from copy import deepcopy
from multiprocessing import Process

from rdkit import Chem

from examples.pi_reactor import ReactorPI
from lib.forcefield.opls import ff
from lib.parse import parse
from lib.reactor.utils import divide_into_molecules
from lib.reactor.utils import set_molecule_id_for_h
from lib.writer.xml_writer import write_xml

description = """PI maker
tests/pi: python ../../make_pi.py -X monomer_pi.xml -M 'A,diamine,c1cc(N)ccc1Oc1ccc(N)cc1' 'B,dianhydride,O=c2oc(=O)c3cc1c(=O)oc(=O)c1cc23' -F opls
Parallel over molecules.
Therefore parallelization is not important if whole system is one huge molecule.
"""

cg_sys, cg_mols, monomers, box, xml = parse()
reactor = ReactorPI(monomers)
reactions = xml.data.get("di")
if not reactions:
    reactions = []
    for edge in cg_sys.edges:
        ti, tj = xml.data['type'][edge[0]], xml.data['type'][edge[1]]
        reaction_type = tuple((monomers[ti]['reaction_type'], monomers[tj]['reaction_type']))
        reactions.append((reaction_type, edge[0], edge[1]))
aa_sys, meta = reactor.process(cg_sys, reactions)
aa_mols = [deepcopy(_) for _ in divide_into_molecules(aa_sys)]
[Chem.SanitizeMol(_) for _ in aa_mols]
aa_mols_h = [Chem.AddHs(m) for m in aa_mols]
aa_mols_h = [set_molecule_id_for_h(mh) for mh in aa_mols_h]
job_lst = aa_mols_h

defaults = {"O": "@atom:178", "N": "@atom239"}


def processing(i):
    molecule = job_lst[i]
    plm_h_ff, bonds, angles, dihedrals = ff(molecule, defaults=defaults)
    write_xml(plm_h_ff, box, bonds, angles, dihedrals, '%06d' % i)


if __name__ == "__main__":
    jobs = [Process(target=processing, args=(i,)) for i in range(len(job_lst))]
    for job in jobs:
        job.start()
    for job in jobs:
        job.join()
