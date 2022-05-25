from rdkit import Chem


def divide_into_molecules(aa_system):
    res = []
    for m in Chem.rdmolops.GetMolFrags(aa_system, asMols=True):
        n = Chem.RWMol()
        for atom in m.GetAtoms():
            n.AddAtom(atom)
        for bond in m.GetBonds():
            n.AddBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondType())
        res.append(n)
    return res


def set_molecule_id_for_h(molecule):
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() != 1:
            for nbr_atom in atom.GetNeighbors():
                if nbr_atom.GetAtomicNum() == 1:
                    nbr_atom.SetIntProp("molecule_id", atom.GetIntProp("molecule_id"))
    return molecule
