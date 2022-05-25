import networkx as nx


def read_cg_topology(cg_system, monomers):
    cg_sys = nx.Graph()
    for bond in cg_system.data['bond']:
        bond_type, i, j = bond[0], bond[1], bond[2]
        type_i, type_j = cg_system.data['type'][i], cg_system.data['type'][j]
        cg_sys.add_node(i, type=type_i, reaction_type=monomers[type_i]['reaction_type'],
                        smiles=monomers[type_i]['smiles'])
        cg_sys.add_node(j, type=type_j, reaction_type=monomers[type_j]['reaction_type'],
                        smiles=monomers[type_j]['smiles'])
        cg_sys.add_edge(i, j, bond_type=bond_type)
    cg_molecules = [cg_sys.subgraph(c).copy() for c in nx.connected_components(cg_sys)]
    return cg_sys, cg_molecules
