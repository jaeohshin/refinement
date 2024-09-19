import logging
import time
from collections import namedtuple

import openmm as mm
from openmm import app
from openmm import unit as u


class LoggerReporter:

    def __init__(self, logger: logging.Logger, report_interval: int, flush=False):
        self._report_interval = report_interval
        self._logger = logger
        self._flush = flush
        self._has_initialized = False

    def describeNextReport(self, simulation: app.Simulation):
        steps = self._report_interval - simulation.currentStep % self._report_interval
        # steps, positions, velocities, forces, energy,
        return steps, False, False, False, True

    def report(self, simulation: app.Simulation, state: mm.State):
        if not self._has_initialized:
            self._initialize_constants(simulation)
            self._init_clock_time = time.time()
            self._init_simul_time = state.getTime()
            self._init_steps = simulation.currentStep
            self._has_initialized = True

        values = self._construct_values(simulation, state)

        self._log(f"Step: {values.step}   Time: {values.time:.3f} ps")
        self._log(f"  {values.elapsed} elapsed   speed: {values.speed} ns/day")
        self._log(f"  KE: {values.ke:.6g}   kJ/mol PE: {values.pe:.6g} kJ/mol")
        self._log(f"  Temperature: {values.temp:.3f} K")
        self._log(f"  Volume: {values.vol:.3f} nm^3   Density: {values.dens:.3f} g/mL")
        for name, e in values.terms.items():
            self.__log(f"  {name} {e}")

    def _log(self, msg):
        self._logger.log(logging.INFO, msg)
        if self._flush:
            for handler in self._logger.handlers:
                handler.flush()

    def _initialize_constants(self, simulation: app.Simulation):
        system = simulation.system
        # Compute the number of degrees of freedom.
        dof = 0
        for i in range(system.getNumParticles()):
            if system.getParticleMass(i) > 0 * u.dalton:
                dof += 3
        for i in range(system.getNumConstraints()):
            p1, p2, distance = system.getConstraintParameters(i)
            if system.getParticleMass(p1) > 0 * u.dalton or system.getParticleMass(p2) > 0 * u.dalton:
                dof -= 1
        if any(type(system.getForce(i)) == mm.CMMotionRemover for i in range(system.getNumForces())):
            dof -= 3
        self._dof = dof
        # Compute the total system mass.
        self._total_mass = 0 * u.dalton
        for i in range(system.getNumParticles()):
            self._total_mass += system.getParticleMass(i)

    def _construct_values(self, simulation: app.Simulation, state: mm.State):
        values = namedtuple('values',
                            ['step', 'time', 'pe', 'ke', 'temp', 'vol', 'dens', 'speed', 'elapsed', 'terms'])
        vol = state.getPeriodicBoxVolume()
        clock_time = time.time()
        # Step
        values.step = simulation.currentStep
        # Time (ps)
        values.time = state.getTime().value_in_unit(u.picosecond)
        # Potential Energy (kJ/mol)
        values.pe = state.getPotentialEnergy().value_in_unit(u.kilojoules_per_mole)
        # Kinetic Energy (kJ/mol)
        values.ke = state.getKineticEnergy().value_in_unit(u.kilojoules_per_mole)
        # Temperature (K)
        values.temp = (2 * state.getKineticEnergy() / (self._dof * u.MOLAR_GAS_CONSTANT_R)).value_in_unit(u.kelvin)
        # Volume (nm^3)
        values.vol = vol.value_in_unit(u.nanometer ** 3)
        # Density (g/mL)
        values.dens = (self._total_mass / vol).value_in_unit(u.gram / u.item / u.milliliter)
        # Speed (ns/day)
        elapsed_days = (clock_time - self._init_clock_time) / 86400
        elapsed_ns = (state.getTime() - self._init_simul_time).value_in_unit(u.nanosecond)
        if elapsed_days > 0:
            values.speed = '{:3g}'.format(elapsed_ns / elapsed_days)
        else:
            values.speed = '--'
        # Elapsed Time
        elapsed_seconds = time.time() - self._init_clock_time
        elapsed_str = time.strftime('%H:%M:%S', time.gmtime(elapsed_seconds))
        values.elapsed = elapsed_str

        system = simulation.system
        terms = {}
        for i, force in enumerate(system.getForces()):
            force_energy = simulation.context.getState(getEnergy=True, groups={1<<i}).getPotentialEnergy()
            terms[force.__name] = force_energy
        values.terms = terms

        return values
