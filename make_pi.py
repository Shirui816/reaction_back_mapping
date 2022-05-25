from copy import deepcopy
from multiprocessing.dummy import Pool

from rdkit import Chem

from examples.pi_reactor import ReactorPI
from lib.forcefield.opls import ff
from lib.parse import parse
from lib.reactor.utils import divide_into_molecules
from lib.reactor.utils import set_molecule_id_for_h
from lib.writer.xml_writer import write_xml


def processing(job):
    molecule, im, box = job
    plm_h_ff, bonds, angles, dihedrals = ff(molecule)
    write_xml(plm_h_ff, box, bonds, angles, dihedrals, '%06d' % im)


if __name__ == "__main__":
    cg_sys, cg_mols, monomers, box, xml = parse()
    reactor = ReactorPI(monomers)
    reactions = xml.data.get("reaction_2")
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
    boxes = [box] * len(aa_mols_h)
    job_lst = list(zip(aa_mols_h, list(range(len(aa_mols))), boxes))
    # for job in job_lst:
    #    processing(job, box)
    p = Pool(40)
    p.map(processing, job_lst)
    p.close()
    p.join()
