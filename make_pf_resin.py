from copy import deepcopy
from multiprocessing import Process

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdChemReactions

from lib.forcefield.opls import ff
from lib.parse_deprecated import parse
from lib.reactor.reactor import Reactor
from lib.reactor.utils import divide_into_molecules
from lib.reactor.utils import set_molecule_id_for_h
from lib.writer.xml_writer import write_xml

__doc__ = """pf resin maker
tests/pfr: python ../../make_pf_resin_deprecated.py -M 'A,A,Oc1ccccc1' 'B,B,C=O' -F opls -X monomerpfr.xml
Parallel over molecules.
Therefore parallelization is not important if whole system is one huge molecule.
reaction templates:
{reaction_type: {reactants_tuple: tuple, reactions: [(reaction, probability, needed production ids),...]}
for CG molecule A-B-A, different class of reactions may take place, thus reaction_type may be
(A, B, A), (A1, B1, A1),... etc. reactants_tuple is set to identify the molecules/
I'll use reaction_type (tuple) as reactants_tuple if it's blank
E.G., for pf resins, reactions may happen on meta-position or para-position of the Ar ring.
Therefore, for certain monomer combination, probability is used for choosing reaction.
needed production ids is a list of ints or None, if None, atoms from all products are preserved.
e.g., A+B->C+D, if D is unwanted, set needed production ids = [0]
The reaction
"""

reaction_templates = {
    # reaction type may be different from reactants_tuple, represents the class of reactions
    # i.e., same CG monomer combination can have different CLASS of reactions.
    # same CLASS of reactions have probabilities summing to 1, e.g., reaction on meta- or para- positions
    ('A', 'B', 'A'): {  # reaction type, name of reaction.
        'reactants_tuple': ('A', 'B', 'A'),  # tuple of types of CG molecules, if not given, use reaction type instead
        'reactions':
            [
                (rdChemReactions.ReactionFromSmarts(
                    "[O:1][c:2][c:3].[C:8]=[O:9].[O:10][c:11][c:12]>>[O:1][c:2][c:3][C:8][c:12][c:11][O:10].[O:9]"
                ), 1.0, [0])
            ]
    },
    ('A', 'B'): {
        'reactions':
            [
                (rdChemReactions.ReactionFromSmarts(
                    "[O:1][c:2][c:3].[C:8]=[O:9]>>[O:1][c:2][c:3][C:8][O:9]"
                ), 1.0, [0])
            ]
    },
    ('B', 'A'): {
        'reactions':
            [
                (rdChemReactions.ReactionFromSmarts(
                    "[C:8]=[O:9].[O:1][c:2][c:3]>>[O:1][c:2][c:3][C:8][O:9]"
                ), 1.0, [0])
            ]
    }
}

cg_sys, cg_mols, monomers, box, xml = parse()
reactor = Reactor(monomers, reaction_templates)

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

defaults = {"O": "@atom:178", "N": "@atom:239"}


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
