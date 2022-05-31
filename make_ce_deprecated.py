from multiprocessing import Pool

from rdkit import Chem

from examples_deprecated.ce_reactor import ReactorCE
from lib.forcefield.opls import ff
from lib.parse_deprecated import parse
from lib.reactor.utils import set_molecule_id_for_h
from lib.writer.xml_writer import write_xml

__doc__ = """CE maker
tests/ce: python ../../make_ce_deprecated.py -X monomerCG.xml -M 'A,A,c1cc(OC#N)ccc1C(C)(C)c1ccc(OC#N)cc1' -F opls
"""

cg_sys, cg_mols, monomers, box, xml = parse()
reactor = ReactorCE(monomers)
dimers = xml.data.get("di")
trimers = xml.data.get('tri')
reactions = []
# for dimer in dimers:
# reactions.append((('A', 'A'), dimer[0], dimer[1]))
for trimer in trimers:
    reactions.append((('A', 'A', 'A'), trimer[0], trimer[1], trimer[2]))
aa_sys, meta = reactor.process(cg_sys, reactions)
aa_mols = list(Chem.rdmolops.GetMolFrags(aa_sys, asMols=True))
[Chem.SanitizeMol(_) for _ in aa_mols]
aa_mols_h = [Chem.AddHs(m) for m in aa_mols]
aa_mols_h = [set_molecule_id_for_h(mh) for mh in aa_mols_h]

defaults = {"O": "@atom:122", "N": "@atom:708"}


def processing(i):
    molecule = aa_mols_h[i]
    plm_h_ff, bonds, angles, dihedrals = ff(molecule, defaults=defaults)
    # TODO: Add SMARTS support, for example, defualts={"O": {SMARTS_1: Type_1}...}
    write_xml(plm_h_ff, box, bonds, angles, dihedrals, '%06d' % i)


if __name__ == "__main__":
    p = Pool(20)
    p.map(processing, range(len(aa_mols_h)))
    p.close()
    p.join()
