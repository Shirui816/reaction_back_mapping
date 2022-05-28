from copy import deepcopy
from multiprocessing import Process

from rdkit import Chem

from examples.pf_resin_reactor import ReactorPFResin
from lib.forcefield.opls import ff
from lib.parse import parse
from lib.reactor.utils import divide_into_molecules
from lib.reactor.utils import set_molecule_id_for_h
from lib.writer.xml_writer import write_xml

description = """PI maker
tests/pfr: python ../../make_pi.py -X monomer_pi.xml -M 'A,A,Oc1ccccc1' 'B,B,C=O' -F opls
Parallel over molecules.
Therefore parallelization is not important if whole system is one huge molecule.
"""

cg_sys, cg_mols, monomers, box, xml = parse()
reactor = ReactorPFResin(monomers)

reactions = []
for monomer in cg_sys.nodes:
    if cg_sys.nodes[monomer]['type'] == 'B':  # CHO
        if len(cg_sys.adj[monomer]) == 2:
            l, r = list(cg_sys.adj[monomer])
            reactions.append((('A', 'B', 'A'), l, monomer, r))
        elif len(cg_sys.adj[monomer]) == 1:
            l = list(cg_sys.adj[monomer])[0]
            reactions.append((('A', 'B'), l, monomer))
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
    # for i in range(len(job_lst)):
    #    print('------------', i, '----------')
    #    processing(i)

    jobs = [Process(target=processing, args=(i,)) for i in range(len(job_lst))]
    for job in jobs:
        job.start()
    for job in jobs:
        job.join()
