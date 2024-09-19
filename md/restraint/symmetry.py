from itertools import combinations, chain, repeat, tee
from typing import List, Tuple

import numpy as np
import openmm as mm
from openmm import app as app
from openmm import unit as u

from .meta import Restraint

__all__ = ['Symmetry', 'corresponding_residues']


def pairwise(iterable):  # Backport
    # pairwise('ABCDEFG') --> AB BC CD DE EF FG
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def corresponding_residues(chain_lengths: List[int], times: int) -> List[Tuple[int, int, int, int]]:
    ncs = chain_lengths
    rep = times
    rp = np.fromiter(chain([0], chain(*repeat(ncs, rep))), dtype=int)
    cs = np.cumsum(rp)
    m0 = [range(a, b) for a, b in pairwise(cs)]
    pairs = list()
    for i in range(len(ncs)):
        for y in zip(*m0[i::len(ncs)]):
            pairs.append(y)
    return pairs


class Symmetry(Restraint):
    def __init__(self, positions: u.Quantity, topology: app.Topology):
        super().__init__()
        self.positions = positions
        self.topology = topology

    def apply(self, system: mm.System, *args, **kwargs):
        """
        weight: float
        cutoff: float
        pairs: List of tuples of equivalent residues.
        """
        weight = float(kwargs['weight'])
        cutoff = float(kwargs['cutoff']) * u.angstrom
        pairs = kwargs['pairs']

        force = mm.CustomCompoundBondForce(4, 'ks*(d1-d2)^2;d1=distance(p1,p2);d2=distance(p3,p4)')
        force.addGlobalParameter('ks', weight * u.kilocalorie_per_mole / u.angstrom ** 2)
        residues = list(self.topology.residues())

        dps = []

        for p1, p2 in combinations(pairs, 2):
            for (r11, r12), (r21, r22) in zip(combinations(p1, 2), combinations(p2, 2)):
                d1 = abs(r11 - r21)
                d2 = abs(r12 - r22)
                if d1 <= 1 and d2 <= 1:
                    continue
                dps.append((r11, r21, r12, r22))

        for r1, r2, r3, r4 in dps:
            res1 = residues[r1]
            res2 = residues[r2]
            res3 = residues[r2]
            res4 = residues[r2]

            a1 = a2 = a3 = a4 = None
            for a1 in res1.atoms():
                if a1.name == 'CA':
                    break
            else:
                Exception('Cannot find CA from residue[{}]'.format(r1))
            for a2 in res2.atoms():
                if a2.name == 'CA':
                    break
            else:
                Exception('Cannot find CA from residue[{}]'.format(r2))
            for a3 in res3.atoms():
                if a3.name == 'CA':
                    break
            else:
                Exception('Cannot find CA from residue[{}]'.format(r3))
            for a4 in res4.atoms():
                if a4.name == 'CA':
                    break
            else:
                Exception('Cannot find CA from residue[{}]'.format(r4))

            i1 = a1.index
            i2 = a2.index
            i3 = a3.index
            i4 = a4.index

            pos1 = self.positions[i1]
            pos2 = self.positions[i2]
            pos3 = self.positions[i3]
            pos4 = self.positions[i4]

            dr1 = (pos1 - pos2).in_units_of(u.angstrom)
            dr1 = np.dot(dr1, dr1)
            dr1 = np.sqrt(dr1)

            dr2 = (pos3 - pos4).in_units_of(u.angstrom)
            dr2 = np.dot(dr2, dr2)
            dr2 = np.sqrt(dr2)

            okay = dr1 <= cutoff or dr2 <= cutoff

            if okay:
                force.addBond((i1, i2, i3, i4))

        system.addForce(force)
