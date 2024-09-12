import os
import sys
err  = sys.stderr
from typing import Optional, Union, List, Tuple
from itertools import combinations

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO

__all__ = ['generate_ssbonded_pdb']

def generate_ssbonded_pdb(
        pdb_input: Union[os.PathLike[str], str],
        pdb_output: Union[os.PathLike[str], str],
        ssbond: str,
        cutoff: Optional[float] = None,
        pair_list: Optional[List[Tuple[int, int]]] = None
) -> List[str]:
    """
    Unit of `cutoff` is Anstrong.
    Assume that there is no CYX residue.
    cutoff is prior than user pair list.
    """
    parser = PDBParser()
    struct = parser.get_structure("", pdb_input)
    io = PDBIO()

    leap_fmt = 'bond {{model}}.{id_x}.SG {{model}}.{id_y}.SG'
    leap_cmds = []

    models = list(struct.get_models())
    model = models[0]
    chains = {c.id: {'seg_id': i, 'num_res': len(list(c.get_residues()))} for i, c in enumerate(model.get_chains())}
    num_res = [(c['seg_id'], c['num_res']) for c in chains.values()]
    num_res.sort()
    num_res = [x[1] for x in num_res]
    cys_list = [r for r in model.get_residues() if r.resname == 'CYS']

    if os.path.exists(ssbond):
        cyx_residues = {}
        f=open(ssbond, 'r')
        for line in f.readlines():
            if line[:2] =='# ': continue
            if line[:2] =='##':
                cid = line.split()[3]
                cyx_residues[cid] = []
                continue

            v = line.split()
            i_res, j_res, i_id, j_id, i_cid, j_cid = int(v[1]), int(v[2]), int(v[4]), int(v[5]), v[6], v[7]
            print (f'ssbond = {i_id} {i_cid} {i_res:4} - {j_id} {j_cid} {j_res:4}', file=err)

            if i_res not in cyx_residues[i_cid]:
                cyx_residues[i_cid].append(i_res)
            if j_res not in cyx_residues[j_cid]:
                cyx_residues[j_cid].append(j_res)

            i_res += sum(num_res[0:i_id], start=0)
            j_res += sum(num_res[0:j_id], start=0)

            s = leap_fmt.format(id_x=i_res, id_y=j_res)
            leap_cmds.append(s)

         # Save
        for cys in cys_list:
            cid = cys.full_id[2]
            resnum = cys.full_id[3][1]
            if cid in cyx_residues and resnum in cyx_residues[cid]:
                cys.resname = 'CYX'

        io.set_structure(struct)
        io.save(pdb_output)
        print ('\n'.join(leap_cmds))
        return leap_cmds       

    #models = list(struct.get_models())
    #model = models[0]
    #chains = {c.id: {'seg_id': i, 'num_res': len(list(c.get_residues()))} for i, c in enumerate(model.get_chains())}
    #num_res = [(c['seg_id'], c['num_res']) for c in chains.values()]
    #num_res.sort()
    #num_res = [x[1] for x in num_res]
    #cys_list = [r for r in model.get_residues() if r.resname == 'CYS']
    sg_list = [a for r in cys_list for a in r.get_atoms() if a.get_id() == 'SG']

    if cutoff is not None:
        pair_list = [(x, y) for x, y in combinations(sg_list, 2) if (x - y) <= cutoff]
    elif pair_list is not None:
        pair_list = list(pair_list)
    else:
        NotImplementedError('Wrong criteria!')

    print(file=sys.stderr)
    for k, v in chains.items():
        print(k, v, file=sys.stderr)

    # Rename CYS to CYX for tleap
    for x, y in pair_list:
        x.get_parent().resname = 'CYX'
        y.get_parent().resname = 'CYX'

        print(file=sys.stderr)
        print(x.get_full_id(), file=sys.stderr)
        print(y.get_full_id(), file=sys.stderr)

        ridx = x.get_parent().get_full_id()[3][1]
        cx = x.get_full_id()[2]
        cx_id = chains[cx]['seg_id']
        ridx += sum(num_res[0:cx_id], start=0)

        ridy = y.get_parent().get_full_id()[3][1]
        cy = y.get_full_id()[2]
        cy_id = chains[cy]['seg_id']
        ridy += sum(num_res[0:cy_id], start=0)

        s = leap_fmt.format(id_x=ridx, id_y=ridy)
        leap_cmds.append(s)

    # Save
    io.set_structure(struct)
    io.save(pdb_output)

    print ('\n'.join(leap_cmds))
    #exit()

    return leap_cmds
