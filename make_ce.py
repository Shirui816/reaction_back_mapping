from rdkit import Chem

from examples import ReactorCE
from lib.forcefield.opls import ff
from lib.parse import parse
from lib.reactor.utils import set_molecule_id_for_h
from lib.writer.xml_writer import write_xml


def processing(job):
    molecule, im, box = job
    plm_h_ff, bonds, angles, dihedrals = ff(molecule)
    write_xml(plm_h_ff, box, bonds, angles, dihedrals, '%06d' % im)


if __name__ == "__main__":
    from multiprocessing.dummy import Pool

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
    [set_molecule_id_for_h(mh) for mh in aa_mols_h]
    boxes = [box] * len(aa_mols_h)
    job_lst = list(zip(aa_mols_h, list(range(len(aa_mols))), boxes))
    p = Pool(20)
    p.map(processing, job_lst)
    p.close()
    p.join()
    # for job in job_lst:
    #    processing(job)
