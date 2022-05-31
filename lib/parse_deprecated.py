import argparse

from lib.cg_reader.xml_parser import XmlParser
from lib.topology.cg_top import read_cg_topology

description = """Backmap CG polymer simulations to AA.
Currently only OPLS force field is supported.
Polymer are generated by SMARTS reactions.
"""


def parse():
    arg_parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # monomer_name is dump, it should be found automatically.
    arg_parser.add_argument('-M', '--molecules',
                            nargs='+',
                            help="molecules with <cg_type>,<monomer_name>,<smiles>...",
                            dest="molecules",
                            type=str,
                            metavar='<cg_type>,<molecule_name>,<smiles>;',
                            required=True,
                            )

    arg_parser.add_argument('-X', '--xml',
                            help='Xml file contains CG bond and type.',
                            dest='xml_file',
                            type=str,
                            metavar='<xml_file_name>',
                            required=True,
                            )

    arg_parser.add_argument('-F',
                            '--force_field',
                            nargs=1,
                            help='OPLS only for now.',
                            dest='force_field',
                            type=str,
                            metavar='<force_field_name>',
                            required=True,
                            )

    args = arg_parser.parse_args()
    molecules = {}
    for line in args.molecules:
        molecule = line.split(',')
        while '' in molecule:
            molecule.remove('')
        molecules[molecule[0]] = {}
        molecules[molecule[0]]['smiles'] = molecule[-1]
        molecules[molecule[0]]['reaction_type'] = molecule[1]
    # TODO: autodetect monomer and automatically choose reactions using SMARTS.

    xml = XmlParser(args.xml_file)
    box = (xml.box.lx, xml.box.ly, xml.box.lz, xml.box.xy, xml.box.xz, xml.box.yz)
    box = tuple(map(float, box))
    cg_sys, cg_mols = read_cg_topology(xml, molecules)
    return cg_sys, cg_mols, molecules, box, xml
