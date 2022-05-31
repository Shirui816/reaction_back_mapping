from copy import deepcopy
from multiprocessing import Process

import numpy as np
from rdkit import Chem

from examples_deprecated.pf_resin_reactor import ReactorPFResin
from lib.forcefield.opls import ff
from lib.parse_deprecated import parse
from lib.reactor.utils import divide_into_molecules
from lib.reactor.utils import set_molecule_id_for_h
from lib.writer.xml_writer import write_xml

__doc__ = """pf resin maker
tests/pfr: python ../../make_pf_resin_deprecated.py -M 'A,A,Oc1ccccc1' 'B,B,C=O' -F opls -X monomerpfr.xml
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
_s = [Chem.SanitizeMol(_) for _ in aa_mols]
assert np.allclose(_s, 0)
aa_mols_h = [Chem.AddHs(m) for m in aa_mols]
aa_mols_h = [set_molecule_id_for_h(mh) for mh in aa_mols_h]
job_lst = aa_mols_h

defaults = {"O": "@atom:178"}


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
