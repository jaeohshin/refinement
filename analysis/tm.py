import re
import os
from typing import AnyStr, Dict, List, Optional, Tuple, Union

import numpy as np
from casp16.utils import run, use_temp

__all__ = ['tm_score']


def tm_score(mobile: Union[str, os.PathLike[str]],
             target: Union[str, os.PathLike[str]],
             args: Optional[List[AnyStr]] = None,
             exec_path: Union[str, os.PathLike[str]] = 'TMscore') -> Dict:
    """
    To rotate structure 1 from (x, y, z) to (X, Y, Z):
    X[i] = t[0] + u[0][0] * x[i] + u[0][1] * y[i] + u[0][2] * z[i];
    """
    cmd = [exec_path]
    if args is not None:
        for x in args:
            cmd.append(x)
    cmd.append(os.path.abspath(mobile))
    cmd.append(os.path.abspath(target))
    cmd.append('-m')
    cmd.append('matrix.txt')

    with use_temp():
        lines = run(cmd).splitlines()
        lines += open('matrix.txt', 'r').readlines()

    r, rot, vec = tm_parser(lines)

    return {'tm_score': r, 'matrix': rot, 'vector': vec}


def tm_parser(lines: List[str]) -> Tuple[float, np.ndarray, np.ndarray]:
    r = None
    re_score = re.compile(r'TM-score\s*=\s*(\d*\.\d*)')
    rowcount = 0
    matrix = []
    line_it = iter(lines)
    header_check = False
    alignment = []
    for line in line_it:
        if 4 >= rowcount > 0:
            if rowcount >= 2:
                a = list(map(float, line.split()))
                matrix.extend(a[2:5])
                matrix.append(a[1])
            rowcount += 1
        elif not header_check and line.startswith(' * '):
            a = line.split(None, 2)
            if len(a) == 3:
                header_check = a[1]
        elif line.lower().startswith(' -------- rotation matrix'):
            rowcount = 1
        elif line.startswith('(":" denotes'):
            alignment = [next(line_it).rstrip() for i in range(3)]
        else:
            match = re_score.search(line)
            if match is not None:
                r = float(match.group(1))
    matrix = np.array(matrix).reshape((3, 4))
    rot = matrix[:, 0:3]
    vec = matrix[:, 3]
    return r, rot, vec
