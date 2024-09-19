import numpy as np
import openmm as mm
from openmm import app as app
from openmm import unit as u

from .meta import Restraint

__all__ = ['Distance']


class Distance(Restraint):
    def __init__(self, positions: u.Quantity, topology: app.Topology):
        super().__init__()
        self.positions = positions
        self.topology = topology

        self.atom_indices = dict()
        for atom in self.topology.atoms():
            res_ind = atom.residue.index
            atom_name = atom.name
            self.atom_indices[res_ind, atom_name] = atom.index

    def apply(self, system: mm.System, *args, **kwargs):
        """
        potential: str = 'harmonic', 'flat_bottom'
        parameters: dict
        weight: float (for harmonic)
        r_cut: float (for flat_bottom)
        """
        potential = str(kwargs['potential']).lower()
        params = kwargs['parameters']

        if potential not in ['harmonic', 'flat_bottom']:
            raise NotImplementedError('"{}" is not implemented.'.format(potential))

        force = None
        if potential == 'harmonic':
            force = mm.CustomBondForce('k*(r-r0)^2')

        elif potential == 'flat_bottom':
            rc = float(kwargs['r_cut']) * u.kilocalorie_per_mole / u.angstrom ** 2
            force = mm.CustomBondForce('k*xi^2*step(xi);xi=dr-rc;dr=abs(r-r0)')
            force.addGlobalParameter('rc', rc)

        else:
            raise NotImplementedError(f'Potential "{potential}" is not implemented to distance restraints.')

        force.addPerBondParameter('k')
        force.addPerBondParameter('r0')

        for (res1, res2), (k, dist) in params.items():
            i1 = self.atom_indices[res1, 'CA']
            i2 = self.atom_indices[res2, 'CA']

            pos1 = self.positions[i1]
            pos2 = self.positions[i2]

            # dr = (pos1 - pos2).in_units_of(u.angstrom)
            # dr2 = dr.dot(dr)
            # dist = np.sqrt(dr2) * u.angstrom
            # k = 745.2 / dr2 * u.kilocalorie_per_mole / u.angstrom ** 2

            dist = dist * u.angstrom
            k = k * u.kilocalorie_per_mole / u.angstrom ** 2
            force.addBond(i1, i2, (k, dist))

        system.addForce(force)
