from multiprocessing import Pool

from rdkit import Chem
from rdkit.Chem import rdChemReactions

from lib.forcefield.opls import ff
from lib.parse_deprecated import parse
from lib.reactor.reactor import Reactor
from lib.reactor.utils import set_molecule_id_for_h
from lib.writer.xml_writer import write_xml

__doc__ = """CE maker
tests/ce: python ../../make_ce.py -X monomerCG.xml -M 'A,A,c1cc(OC#N)ccc1C(C)(C)c1ccc(OC#N)cc1' -F opls
(reaction_type: (reaction, probability, needed production ids))
"""

reaction_templates = {
    ('A', 'A'): {
        'reactants_tuple': ('A', 'A'),
        'reactions':
            [
                (rdChemReactions.ReactionFromSmarts(
                    "[C:1]#[N:2].[C:3]#[N:4]>>[C:1]=[N:2][C:3]=[N:4]"
                ), 1.0, None)
            ],
    },
    ('A', 'A', 'A'): {
        'reactions':
            [
                (rdChemReactions.ReactionFromSmarts(
                    "[C:1]#[N:2].[C:3]#[N:4].[C:5]#[N:6]>>[c:1]1[n:2][c:3][n:4][c:5][n:6]1"
                ), 1.0, None)
            ]
    }
}

cg_sys, cg_mols, monomers, box, xml = parse()
reactor = Reactor(monomers, reaction_templates)
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

defaults = {"O": "@atom:122"}


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
