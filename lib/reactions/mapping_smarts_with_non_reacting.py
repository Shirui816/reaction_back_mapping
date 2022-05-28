from collections import namedtuple

from lib.reactions.mapping import map_reacting_atoms

# TODO: this works where SMARTS contains non-reacting atoms.


atom_info = namedtuple('atom_info',
                       ('reactant_id', 'reactant_atom_id', 'product_id', 'product_atom_id', 'atomic_number'))
bond_info = namedtuple('bond_info',
                       ('product_id', 'product_atoms_id', 'reactants_id', 'reactant_atoms_id', 'bond_type', 'status'))


def map_atoms1(products, reacting_map):  # reacted only
    r"""
    :param products:  list of products
    :param reacting_map: reacting atom map
    :return: atom map and reacting atoms
    """
    amap = []
    ra_dict = {}
    for ip, p in enumerate(products):
        p_idx = ip
        for a in p.GetAtoms():
            p_aidx = a.GetIdx()
            old_mapno = a.GetPropsAsDict().get('old_mapno')
            r_aidx = a.GetPropsAsDict().get("react_atom_idx")
            if old_mapno is not None:  # reacting atoms
                r_idx = reacting_map[old_mapno]
                if ra_dict.get(r_idx) is None:
                    ra_dict[r_idx] = []
                ra_dict[r_idx].append(r_aidx)
                amap.append(atom_info(r_idx, r_aidx, p_idx, p_aidx, a.GetAtomicNum()))  # all Atoms
    return amap, ra_dict


def atom_map1(products, reaction):
    r"""Map atoms between reactants and products
    :param products: product list
    :param reaction: reaction object
    :return: atom map and reacting atoms
    """
    reacting_map = map_reacting_atoms(reaction)
    amap, reacting_atoms = map_atoms1(products, reacting_map)
    return amap, reacting_atoms


def bond_map1(reactants, products, reaction):
    r"""Map bonds between reactants and products
    :param reactants: list of reactants
    :param products:  list of products
    :param reaction:  reaction object
    :return:  list of bond map
    """
    res = []
    amap, reacting_atoms = atom_map1(products, reaction)
    for ir, r in enumerate(reactants):
        for bond in r.GetBonds():  # exist in reactants, but not in production
            a_pid = b_pid = None
            a, b, t = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondType()
            for atom in amap:
                if atom.reactant_id == ir:
                    if atom.reactant_atom_id == a:
                        a_pid = atom.product_id
                        a_paid = atom.product_atom_id
                    if atom.reactant_atom_id == b:
                        b_pid = atom.product_id
                        b_paid = atom.product_atom_id
            if a_pid is None or b_pid is None:
                continue
            if a_pid != b_pid:
                res.append(bond_info(a_pid, (a_paid, b_paid), (ir, ir), (a, b), t, 'deleted'))
            if a_pid == b_pid:
                p_bond = products[a_pid].GetBondBetweenAtoms(a_paid, b_paid)
                if p_bond is None:
                    res.append(bond_info(a_pid, (a_paid, b_paid), (ir, ir), (a, b), t, 'deleted'))

    for ip, p in enumerate(products):
        for bond in p.GetBonds():  # exist in productions, compare with reactants
            a_rid = b_rid = None
            a, b, t = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondType()
            for atom in amap:
                if atom.product_id == ip:
                    if atom.product_atom_id == a:
                        a_rid = atom.reactant_id
                        a_raid = atom.reactant_atom_id
                    if atom.product_atom_id == b:
                        b_rid = atom.reactant_id
                        b_raid = atom.reactant_atom_id
            if a_rid is None or b_rid is None:
                continue
            if a_rid != b_rid:
                res.append(bond_info(ip, (a, b), (a_rid, b_rid), (a_raid, b_raid), t, 'new'))
            if a_rid == b_rid:
                r_bond = reactants[a_rid].GetBondBetweenAtoms(a_raid, b_raid)
                if r_bond is None:
                    res.append(bond_info(ip, (a, b), (a_rid, b_rid), (a_raid, b_raid), t, 'new'))
                else:
                    if r_bond.GetBondType() != t:
                        res.append(bond_info(ip, (a, b), (a_rid, b_rid), (a_raid, b_raid), t, 'changed'))
    return res
