import os
import sys
from typing import Dict, Union, Optional

from .leap import call_tleap
from .disulfide import generate_ssbonded_pdb as gen_ss

__all__ = ['gen_init']


def gen_init(pdb_input: Union[os.PathLike[str], str],
             pdb_output: Union[os.PathLike[str], str],
             prmtop: Union[os.PathLike[str], str],
             inpcrd: Union[os.PathLike[str], str],
             ssbond = 'SSBOND',
             tleap_path: Optional[Union[os.PathLike[str], str]] = None,
             include_dir: Optional[Union[os.PathLike[str], str]] = None,
             ss_kwargs: Optional[Dict] = None):

    cmds = gen_ss(pdb_input, 'ss.pdb', ssbond='SSBOND', **ss_kwargs)

    print('\n  Amber Di-Sulfide Commands:', file=sys.stderr)
    for cmd in cmds:
        print('   ', cmd, file=sys.stderr)

    call_tleap('ss.pdb', pdb_output, prmtop, inpcrd,
               include_dir=include_dir,
               tleap_path=tleap_path,
               ss_cmds=cmds)
    os.remove('ss.pdb')
