import networkx as nx
import numpy as np
from rdkit import Chem

from lib.reactions.mapping import process_reactants, atom_map, bond_map


class Reactions(object):
    def __init__(self, reaction_type, reaction_meta):
        self.reaction_type = reaction_type
        self.reaction_meta = reaction_meta
        self.index = list(range(len(reaction_meta)))
        prob = [_[1] for _ in reaction_meta]
        self.prob = np.array(prob) / sum(prob)
        self.reaction_maps = {}
        self.product_idx = {}
        self.all_reactions = set(self.index)
        for idx in self.index:
            self.reaction_maps[idx] = []
            self.product_idx[idx] = reaction_meta[idx][-1]

    def process(self, reactants_meta):
        for idx in self.index:
            meta = self.reaction_meta[idx]
            reaction = meta[0]
            molecules = []
            for t in self.reaction_type:
                molecules.append(Chem.MolFromSmiles(reactants_meta[t]['smiles']))
            reactants = process_reactants(molecules)  # new copy of molecules every time
            products = reaction.RunReactants(reactants)
            for product in products:
                amap, reacting_atoms = atom_map(product, reaction)
                bmap = bond_map(reactants, product, reaction)
                self.reaction_maps[idx].append((reacting_atoms, amap, bmap))

    def choose(self, idx=None):
        if idx is None:
            idx = np.random.choice(self.index, p=self.prob)
        return idx, self.reaction_maps[idx], self.product_idx[idx]


def allowed_p(reacted_atoms, reactions):
    seen = set()
    while True:
        if seen == reactions.all_reactions:
            break
        idx, reaction_maps, product_idx = reactions.choose()  # choose randomly
        if idx in seen:
            continue
        seen.add(idx)
        for reaction_map in reaction_maps:
            allowed = True
            for ri in reaction_map[0]:
                if set.issubset(set(reaction_map[0][ri]), reacted_atoms[ri]):  # only idle function groups
                    # multi-step reactions are considered as reactions with all reactants in one step.
                    allowed = False
            if allowed:
                return reaction_map, product_idx
    return None, None  # if no available reaction is chosen.


class Reactor(object):
    def __init__(self, reactants_meta, reaction_templates):
        r"""
        :param reactants_meta: dict {type: {smiles: ''}
        :param reaction_templates: dict {reaction_type: [(reaction_object, probability,
        product_idx (None for all)),],...}, if one reaction_type is mapped to some reactions
        """
        self.cg_system = None
        self.aa_system = Chem.RWMol()
        self.meta = nx.Graph()
        self.reaction_templates = {}
        for reaction_type in reaction_templates:
            self.reaction_templates[reaction_type] = Reactions(reaction_type, reaction_templates[reaction_type])
            self.reaction_templates[reaction_type].process(reactants_meta)

    def process(self, cg_system, reactions):
        r"""
        :param cg_system: nx.Graph(global_molecule_id, smiles=smiles)
        :param reactions: list of reactions, [(reaction_type, i, j, k,...),]
        :return: aa_system and meta info
        """
        self.cg_system = cg_system
        self.aa_system = Chem.RWMol()
        self.meta = nx.Graph()
        global_count = 0
        for node in cg_system.nodes:  # add all monomers first
            atom_idx = {}
            reactant = cg_system.nodes[node]
            reactant_molecule = Chem.MolFromSmiles(reactant['smiles'])
            for atom_id in range(reactant_molecule.GetNumAtoms()):
                atom = reactant_molecule.GetAtomWithIdx(atom_id)
                self.aa_system.AddAtom(atom)
                atom_idx[atom_id] = atom_id + global_count
            for bond in reactant_molecule.GetBonds():
                self.aa_system.AddBond(bond.GetBeginAtomIdx() + global_count, bond.GetEndAtomIdx() + global_count,
                                       bond.GetBondType())
            global_count += reactant_molecule.GetNumAtoms()
            self.meta.add_node(node, atom_idx=atom_idx, reacting_map={}, rm_atoms=set())
        for edge in cg_system.edges:
            self.meta.add_edge(*edge)
        for r in reactions:
            reaction_type = r[0]
            reactant_idx = r[1:]
            rxn_tpls = self.reaction_templates.get(reaction_type)
            if rxn_tpls is None:
                raise ValueError(f"Reaction {r} is not defined in reactions!")

            reactants = [self.meta.nodes[_] for _ in reactant_idx]
            key = tuple(sorted(reactant_idx))  # store reacted atoms
            reacted_atoms = {}

            for ri in range(len(reactants)):
                reacted_atoms[ri] = set()

            for ri, rt in enumerate(reactants):  # keep reactant order
                for k in rt['reacting_map']:
                    for at in rt['reacting_map'][k]:
                        reacted_atoms[ri].add(at)

            reaction_map, product_idx = allowed_p(reacted_atoms, rxn_tpls)
            if not reaction_map:
                raise (ValueError(f"{reactant_idx} can not react!"))

            amap, bmap = reaction_map[1], reaction_map[2]
            for ri, rt in enumerate(reactants):
                if rt['reacting_map'].get(key) is None:
                    rt['reacting_map'][key] = set()
                for at in reaction_map[0][ri]:
                    rt['reacting_map'][key].add(at)

            for atom in amap:
                if product_idx is not None:
                    if atom.product_id not in product_idx:
                        reactant = reactants[atom.reactant_id]
                        reactant['rm_atoms'].add(reactant['atom_idx'][atom.reactant_atom_id])
            for b in bmap:
                if b.status == 'deleted':
                    reactant = reactants[b.reactants_id[0]]
                    bi, bj = reactant['atom_idx'][b.reactant_atoms_id[0]], reactant['atom_idx'][b.reactant_atoms_id[1]]
                    self.aa_system.RemoveBond(bi, bj)
                if b.status == 'changed':
                    reactant = reactants[b.reactants_id[0]]  # bond changed, in same reactant
                    bi, bj = reactant['atom_idx'][b.reactant_atoms_id[0]], reactant['atom_idx'][b.reactant_atoms_id[1]]
                    bond = self.aa_system.GetBondBetweenAtoms(bi, bj)
                    bond.SetBondType(b.bond_type)
                if b.status == 'new':
                    reactant0 = reactants[b.reactants_id[0]]
                    reactant1 = reactants[b.reactants_id[1]]
                    bi = reactant0['atom_idx'][b.reactant_atoms_id[0]]
                    bj = reactant1['atom_idx'][b.reactant_atoms_id[1]]
                    self.aa_system.AddBond(bi, bj, b.bond_type)

        rm_all = []
        for m in self.meta.nodes:  # remove atoms in side productions
            molecule = self.meta.nodes[m]
            for idx in molecule['atom_idx'].values():
                atom = self.aa_system.GetAtomWithIdx(idx)
                atom.SetIntProp('molecule_id', m)
            rm_all.extend(list(molecule['rm_atoms']))

        rm_all = sorted(list(set(rm_all)), reverse=True)
        for bi in rm_all:
            self.aa_system.RemoveAtom(bi)
        return Chem.RemoveAllHs(self.aa_system), self.meta
