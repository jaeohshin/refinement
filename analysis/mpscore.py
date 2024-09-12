import os
from typing import Dict, Union

from casp16.utils import run


def mp_score(target: Union[str, os.PathLike[str]], full=True,
             exec_path: Union[str, os.PathLike[str]] = 'mpscore.py') -> Dict:
    cmd = [exec_path, target]
    if full:
        cmd.append('full')
    result = run(cmd).split()
    return {'mp_score': float(result[5])}
