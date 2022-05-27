import os
import re
import warnings

from rdkit import Chem

__doc__ = """
Ported from  asteeves/rdlt
"""

this_dir, this_filename = os.path.split(__file__)
opls_tpl_path = os.path.join(this_dir, "assets", "opls", "STaGE_opls_tomoltemplate_opls.txt")
opls_bonded_itp = os.path.join(this_dir, "assets", "opls", "bonded.itp")  # from gromacs


# lopls_tpl_path = 'lopls_tomoltemplate.txt'


def generate_feature_defn(fpath):
    """Write a feature definition file in RDKit style from the moltemplate
    conversion document. Only need to run this function if the conversion
    document has been changed.
    fpath -- file path of the moltemplate conversion doc
    fdefout -- file path to write the feature definition file
    cdictout -- file path to create a dictionary of atom types to charges
    """
    info = {}
    fdefn = ''
    with open(fpath, 'r') as infile:
        feat_index = 0
        for line in [line.strip() for line in infile if line.strip()]:
            if line[0] != '*':
                el, atomname, typename, patt, lttype, chg, desc = [el.strip() for el in line.split("|")]
                info[lttype] = (atomname, typename, chg, desc)
                # write lttype, SMARTS to feature definintion file
                # NOTE: using feature family to store the atom names is dangerous
                # because rdkit won't assign mutliple features in same family.
                # So I had to assign an index to make unique names [AHS]
                fdefn += \
                    """
                    DefineFeature {0} {1}
                    Family {2}{3}
                    EndFeature""".format(lttype, patt, feat_index, atomname)
                # add new charge dictionary entry
                feat_index += 1
    return info, fdefn


def ff(molecule, **kwargs):
    info, fdefn = generate_feature_defn(opls_tpl_path)

    factory = Chem.ChemicalFeatures.BuildFeatureFactoryFromString(fdefn)
    features = factory.GetFeaturesForMol(molecule)

    [molecule.GetAtomWithIdx(f.GetAtomIds()[0]).SetProp('AtomType', f.GetType()) for f in features]
    defaults = kwargs['defaults']

    for atom in molecule.GetAtoms():
        try:
            # print("Atom {0} has {1}".format(at.GetIdx(), at.GetProp('AtomType')))
            _ = atom.GetProp("AtomType")
        except KeyError:
            # print("Atom {0}:{1} does not have an assigned atom type!".format(atom.GetIdx(), atom.GetSymbol()))
            _m = Chem.RWMol()
            i_a = count = 0
            _m.AddAtom(atom)
            count += 1
            for nbr_atom in atom.GetNeighbors():
                j = count
                _m.AddAtom(nbr_atom)
                count += 1
                _m.AddBond(i_a, j, molecule.GetBondBetweenAtoms(atom.GetIdx(), nbr_atom.GetIdx()).GetBondType())
                for _nbr_atom in nbr_atom.GetNeighbors():
                    if _nbr_atom.GetIdx() != atom.GetIdx():
                        k = count
                        _m.AddAtom(_nbr_atom)
                        count += 1
                        _m.AddBond(j, k,
                                   molecule.GetBondBetweenAtoms(_nbr_atom.GetIdx(), nbr_atom.GetIdx()).GetBondType())
            smiles = Chem.MolToSmiles(_m)
            msg = "---------------------------------------------------------\n" \
                  "Chemical env of atom {0}({1}) is\n\t{2}\nType is not found from database, set manually as\n\t{3}"
            # default = input("Enter default atom type (@atom:xxx): ")
            # at.SetProp("AtomType", default)
            s = "*** NO DEFAULT TYPE ***"
            if not defaults is None:
                if not defaults.get(atom.GetSymbol()) is None:
                    atom.SetProp("AtomType", defaults.get("O"))
                    dt = defaults.get("O")
                    s = info[dt][0] + ", " + info[dt][3]
            print(msg.format(atom.GetIdx(), atom.GetSymbol(), smiles, s))

            # TODO: generate SMARTS of un-typed atoms, set default for all cases.

    sum_of_charge = 0
    for atom in molecule.GetAtoms():
        names = info[atom.GetProp("AtomType")]
        atom.SetProp("ElementType", names[0])
        atom.SetProp("OplsType", names[1])
        atom.SetProp("Charge", names[2])
        sum_of_charge += float(names[2])
    if abs(sum_of_charge) > 0.01:
        warnings.warn("*** Total Charge is %.6f !" % sum_of_charge)

    bonded_itp = open(opls_bonded_itp, 'r').read()
    # TODO: make a reader and hash map, but for small systems re is fine
    bonds = []
    angles = []
    dihedrals = []
    for bond in molecule.GetBonds():  # bonds
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        atom_i = molecule.GetAtomWithIdx(i)
        atom_j = molecule.GetAtomWithIdx(j)

        bonds.append('%s-%s %d %d' % (atom_i.GetProp("ElementType"),
                                      atom_j.GetProp("ElementType"),
                                      i, j))

        for nbr_i in atom_i.GetNeighbors():  # a bond (i, j) and neighbor of i and j define a dihedral
            ni = nbr_i.GetIdx()
            if ni == j:
                continue
            for nbr_j in atom_j.GetNeighbors():
                nj = nbr_j.GetIdx()
                if nj == i or ni == nj:  # avoid 3-ring
                    continue
                dihedral = [nbr_i.GetProp("ElementType"),
                            atom_i.GetProp("ElementType"),
                            atom_j.GetProp("ElementType"),
                            nbr_j.GetProp("ElementType")]
                if re.search(r'\s*'.join(dihedral) + r'\s*[1-9]+', bonded_itp) or re.search(
                        r'\s*'.join(dihedral[::-1]) + r'\s*[1-9]+', bonded_itp):
                    dihedrals.append('%s %d %d %d %d' % ('-'.join(dihedral), ni, i, j, nj))

    for atom in molecule.GetAtoms():  # any atom can be a center of an angle
        for nbr_i in atom.GetNeighbors():
            for nbr_j in atom.GetNeighbors():
                if nbr_i.GetIdx() <= nbr_j.GetIdx():
                    continue
                angle = [nbr_i.GetProp("ElementType"), atom.GetProp("ElementType"), nbr_j.GetProp("ElementType")]
                if re.search(r'\s*'.join(angle) + r'\s*[1-9]+', bonded_itp) or re.search(
                        r'\s*'.join(angle[::-1]) + r'\s*[1-9]+', bonded_itp):
                    angles.append('%s %d %d %d' % ('-'.join(angle), nbr_i.GetIdx(), atom.GetIdx(), nbr_j.GetIdx()))

    return molecule, bonds, angles, dihedrals
