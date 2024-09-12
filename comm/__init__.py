from typing import List, Optional, Union
import os

from . import send_mail
from . import sort_pdb
from .header import *
from .gmail import *

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO

from casp15.utils import use_temp

__all__ = ['submit', 'prepare_pdb', 'Header']


def prepare_pdb(
        head: Header,
        pdb_file: Union[os.PathLike[str], str],
        model_number: int,
        output: str = 'comm.pdb',
        template_pdb_file: Optional[str] = None,
        crf_path: Optional[str] = None
):
    if crf_path is not None:
        with open(crf_path, 'r') as f:
            parents = f.read().strip().split()
        parent_line = 'PARENT ' + ' '.join(parents)
    else:
        parent_line = 'PARENT N/A'

    template_pdb_file = os.path.abspath(template_pdb_file)

    pdb_path = os.path.abspath(pdb_file)
    with use_temp():
        parser = PDBParser()
        struct1 = parser.get_structure('', template_pdb_file)
        struct2 = parser.get_structure('', pdb_path)
        models1 = list(struct1.get_models())
        models2 = list(struct2.get_models())
        sort_pdb.sort_model(models1[0], models2[0])
        io = PDBIO()
        io.set_structure(struct1)
        io.save('tmp.pdb', preserve_atom_numbering=True, write_end=False)
        body = sort_pdb.get_atoms('tmp.pdb')

    out = '\n'.join(head.get_header())
    out += '\n'
    out += 'MODEL  {n:1d}\n'.format(n=model_number)
    out += parent_line
    out += '\n'
    out += body
    out += '\n'
    out += 'END'

    with open(output, 'w') as f:
        f.write(out)


def submit():
    pass


if __name__ == '__main__':
    submit()
