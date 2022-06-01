from abc import ABCMeta, abstractmethod

import networkx as nx
from rdkit import Chem

from lib.reactions.mapping import process_reactants, atom_map, bond_map

__doc__ = """Deprecated, newer and more generalized method is build. However, this method are more flexible
"dirty" ways for complex systems.
"""

class Reactor(metaclass=ABCMeta):
    def __init__(self, reactants, reaction_templates):
        self.reactants = reactants
        self.reaction_templates = reaction_templates
        self.meta = None
        atom_bond_maps = {}
        for reaction_type in reaction_templates:  # process all reactions
            atom_bond_maps[reaction_type] = []
            reaction = reaction_templates[reaction_type]
            molecules = []
            for rt in reaction_type:
                for t in self.reactants:
                    if self.reactants[t]['reaction_type'] == rt:
                        smiles = self.reactants[t]['smiles']
                molecules.append(Chem.MolFromSmiles(smiles))
            reactants = process_reactants(molecules)
            products = reaction.RunReactants(reactants)
            for product in products:  # record all possible outcomes
                amap, reacting_atoms = atom_map(product, reaction)
                bmap = bond_map(reactants, product, reaction)
                atom_bond_maps[reaction_type].append((reacting_atoms, amap, bmap))
        self.atom_bond_maps = atom_bond_maps
        self.aa_system = None
        self.cg_system = None
        self.reactions = None

    def process(self, cg_system, reactions, *args):
        r"""
        :param cg_system: a cg system graph
        :param reactions: cg reaction data, (reaction_type, cg_i, cg_j)
        :param args: other args pass to self defined process
        :return:
        """
        self.cg_system = cg_system
        self.aa_system = Chem.RWMol()
        self.meta = nx.Graph()
        self.reactions = reactions
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
        return self.react(*args)

    @abstractmethod
    def react(self, *args):
        pass
