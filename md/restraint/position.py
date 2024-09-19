import openmm as mm
from openmm import app as app
from openmm import unit as u

from .meta import Restraint

__all__ = ['Position']


class Position(Restraint):
    def __init__(self, positions: u.Quantity, topology: app.Topology):
        super().__init__()
        self.positions = positions
        self.topology = topology

    def apply(self, system: mm.System, *args, **kwargs):
        """
        potential: str = 'harmonic', 'flat_bottom'
        mode: str = 'ca', 'backbone', 'heavy'
        weight: float
        r_cut: float
        """
        mode = str(kwargs['mode']).lower()
        potential = str(kwargs['potential']).lower()
        if isinstance(kwargs['positions'], u.Quantity):
            self.positions = kwargs['positions']

        if mode not in ['ca', 'backbone', 'heavy']:
            raise NotImplementedError('"{}" is not implemented.'.format(mode))

        force = None
        k = float(kwargs['weight']) * u.kilocalorie_per_mole / u.angstrom ** 2
        if potential == 'harmonic':
            force = mm.CustomExternalForce('k*periodicdistance(x,y,z,x0,y0,z0)^2')
            force.addGlobalParameter('k', k)
        if potential == 'flat_bottom':
        #    rc = float(kwargs['r_cut']) * u.angstrom ##April 2nd, 16:45.# April 4th, changed rc=2 A, April 8th: rc=3
#            rc = 3.0 * u.angstrom
            ## April 23rd. changed rc=0.3 (=3 Anstrom)
#            force = mm.CustomExternalForce('k*xi^2*step(xi);xi=dr-0.4 ;dr=periodicdistance(x,y,z,x0,y0,z0)')
## Aug. 19th. I changed the flat-bottom potential range to 0.2 from 0.4
            force = mm.CustomExternalForce('k*xi^2*step(xi);xi=dr-0.3 ;dr=periodicdistance(x,y,z,x0,y0,z0)')
            force.addGlobalParameter('k', k)
        #    force.addGlobalParameter('rc', rc)
        force.addPerParticleParameter('x0')
        force.addPerParticleParameter('y0')
        force.addPerParticleParameter('z0')

        if mode == 'ca':
            for atom in self.topology.atoms():
                if not atom.name == 'CA':
                    continue
                index = atom.index
                x0, y0, z0 = self.positions[index]
                force.addParticle(index, (x0, y0, z0))

        elif mode == 'backbone':
            for atom in self.topology.atoms():
                if atom.name not in ('N', 'CA', 'C', 'O'):
                    continue
                if atom.name =='O' and atom.residue.name == 'HOH': 
                    continue
                ## April. 4th. Added above two lines: 'O' atoms in the water molecules are free to move.
                index = atom.index
                x0, y0, z0 = self.positions[index]
                force.addParticle(index, (x0, y0, z0))

        elif mode == 'heavy':
            for atom in self.topology.atoms():
                if atom.residue.name in ('WAT', 'HOH', 'SOL', 'Na', 'Cl', 'Na+', 'Cl-'):
                    continue
                if atom.name[0] == 'H':
                    continue
                index = atom.index
                x0, y0, z0 = self.positions[index]
                force.addParticle(index, (x0, y0, z0))

        else:
            raise NotImplementedError(f'Mode "{mode}" is not implemented for position restraints.')

        system.addForce(force)
