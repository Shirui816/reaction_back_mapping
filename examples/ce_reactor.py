from rdkit import Chem
from rdkit.Chem import rdChemReactions

from lib.reactor.reactor import Reactor

reaction_templates = {
    ('A', 'A'): rdChemReactions.ReactionFromSmarts(
        "[C:1]#[N:2].[C:3]#[N:4]>>[C:1]=[N:2][C:3]=[N:4]"
    ),
    ('A', 'A', 'A'): rdChemReactions.ReactionFromSmarts(
        "[C:1]#[N:2].[C:3]#[N:4].[C:5]#[N:6]>>[c:1]1[n:2][c:3][n:4][c:5][n:6]1"
    )
}


class ReactorCE(Reactor):
    def __init__(self, monomers, product_id=0):
        super().__init__(monomers, reaction_templates)
        self.product_id = product_id

    def react(self):
        for r in self.reactions:
            reaction_type = r[0]
            reactant_idx = r[1:]
            atom_bond_map = self.atom_bond_maps[reaction_type]
            reactants = [self.meta.nodes[_] for _ in reactant_idx]
            key = tuple(sorted(reactant_idx))  # store reacted atoms

            reacted_atoms = {}
            for ri in range(len(reactants)):
                reacted_atoms[ri] = set()

            for ri, rt in enumerate(reactants):  # keep reactant order
                for k in rt['reacting_map']:
                    for at in rt['reacting_map'][k]:
                        reacted_atoms[ri].add(at)

            allowed = False
            reaction_map = 'dummy'
            for reaction_map in atom_bond_map:
                allowed = True
                for ri in reaction_map[0]:
                    if set.issubset(set(reaction_map[0][ri]), reacted_atoms[ri]):
                        allowed = False
                if allowed:
                    break
            if not allowed:
                raise (ValueError(f"{reactant_idx} can not react!"))

            amap, bmap = reaction_map[1], reaction_map[2]
            for ri, rt in enumerate(reactants):
                if rt['reacting_map'].get(key) is None:
                    rt['reacting_map'][key] = set()
                for at in reaction_map[0][ri]:
                    rt['reacting_map'][key].add(at)

            for atom in amap:
                if atom.product_id != self.product_id:
                    reactant = reactants[atom.reactant_id]
                    reactant['rm_atoms'].add(reactant['atom_idx'][atom.reactant_atom_id])
            for b in bmap:
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
