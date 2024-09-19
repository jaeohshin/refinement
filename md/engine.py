import sys
from cmath import log
import os
from typing import Any, Union, Optional, Dict, List

import openmm as mm
from openmm import app as app
from openmm import unit as u
import mdtraj.reporters
import numpy as np
import re

import parmed

from casp16.pipeline import TaskStopError
from .reporter import *
from .restraint import *

__all__ = ['Constructor', 'Minimize', 'RunMD', 'GenerateState']

def get_connected_atoms(iatm,connect,catms):
    for jatm in connect[iatm-1]:
        if jatm < iatm:
            continue
        catms.add(jatm)
        get_connected_atoms(jatm,connect,catms)

def get_atom_inds(residue):
    inds = {}
    for atom in residue.atoms:
        inds[atom.name] = atom.idx
    return inds

def standard_radian(inp):
    """ shift randian to [-pi,pi] """
    twopi = 2.0*np.pi
    inp -= int(inp/twopi) * twopi
    if inp >= np.pi:
        inp -= twopi
    elif inp < -np.pi:
        inp += twopi
    return inp

def get_torsion_inds(residues, ires):
    residue = residues[ires]
    inds = get_atom_inds(residue)

    # prev_residue
    if 'H2' in inds:
        prev_residue = None
    else:
        try:
            prev_residue = residues[ires-1]
        except:
            prev_residue = None
    if prev_residue is not None:
        prev_inds = get_atom_inds(prev_residue)

    # next_residue
    if 'OXT' in inds:
        next_residue = None
    else:
        try:
            next_residue = residues[ires+1]
        except:
            next_residue = None
    if next_residue is not None:
        next_inds = get_atom_inds(next_residue)

    resname = residue.name
    # phi
    if prev_residue is not None and 'C' in prev_inds:
        phi_inds = ( prev_inds['C'], inds['N'], inds['CA'], inds['C'] )
    else:
        phi_inds = None
    # psi
    if next_residue is not None and 'N' in next_inds:
        psi_inds = ( inds['N'], inds['CA'], inds['C'], next_inds['N'] )
    else:
        psi_inds = None
    # chi1
    


    if resname in ('ARG', 'ASN', 'ASP', 'ASH', 'GLN', 'GLU', 'GLH', 'HIS',
            'HID', 'HIP', 'HIE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'TRP',
            'TYR'):
        ag = 'CG'
    elif resname in ('ILE', 'VAL'):
        ag = 'CG1'
    elif resname == 'THR':
        ag = 'OG1'
    elif resname in ('CYS', 'CYX'):
        ag = 'SG'
    elif resname == 'SER':
        ag = 'OG'
    else:
        ag = None
    if ag is not None:
        chi1_inds = ( inds['N'], inds['CA'], inds['CB'], inds[ag] )
    else:
        chi1_inds = None



    # chi2
    if resname in ('ARG', 'GLN', 'GLU', 'GLH', 'LYS', 'PRO'):
        ad = 'CD'
    elif resname in ('LEU', 'PHE', 'TRP', 'TYR'):
        ad = 'CD1'
    elif resname in ('ASN', 'ASP', 'ASH'):
        ad = 'OD1'
    elif resname in ( 'HIS', 'HIP', 'HIE', 'HID' ):
        ad = 'ND1'
    elif resname == 'MET':
        ad = 'SD'
    elif resname == 'ILE':
        ad = 'CD1'
    else:
        ad = None
    if ad is not None:
        chi2_inds = ( inds['CA'], inds['CB'], inds[ag], inds[ad] )
    else:
        chi2_inds = None

    return phi_inds, psi_inds, chi1_inds, chi2_inds

def read_stap_table(tfile):
    table = []
    with open(tfile,'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#'):
                continue
            table.append(list(map(float,line.split())))
    return (
            np.array(table) *
            u.kilocalorie_per_mole
    ).value_in_unit(u.kilojoule_per_mole)

class Constructor:
    """
    MD system constructor.
    Also prepare restraints.
    """

    def __init__(self,
                 prmtop: Union[os.PathLike[str], str],
                 inpcrd: Union[os.PathLike[str], str],
                 langevin_kwargs: Optional[Dict[str, Any]] = None,
                 system_kwargs: Optional[Dict[str, Any]] = None,
                 constant_pressure=True,
                 restraints: Optional[List[Dict[str, Any]]] = None,
                 stap: bool = False, 
                 stap_weight: float = 0.5,
                 restrain_forces= {'restrain_path': None, 'restrain_weights':None},
                 **kwargs):
        self.prmtop = prmtop
        self.inpcrd = inpcrd
        self.stap = stap
        self.stap_weight = stap_weight

        self.restrain_path = restrain_forces['restrain_path']
        self.restrain_weights = restrain_forces['restrain_weights']

        if langevin_kwargs is not None:
            self.langevin_kwargs = langevin_kwargs
        else:
            self.langevin_kwargs = dict()

        if system_kwargs is not None:
            self.system_kwargs = system_kwargs
        else:
            self.system_kwargs = dict()

        if constant_pressure:
            self.constant_pressure = True

        self.crd = app.AmberInpcrdFile(self.inpcrd)
        self.top = app.AmberPrmtopFile(self.prmtop)

        self.parmed_structure = parmed.load_file (self.prmtop)

        self.system = None
        self.integrator = None

        if restraints is not None:
            self.restraints = restraints
        else:
            self.restraints = dict()

        self.__dict__.update(kwargs)

    def _generate_system(self):
        self.positions = self.crd.getPositions()

        langevin_kwargs = {
            'temperature': 300.00 * u.kelvin,
            'frictionCoeff': 1.0 / u.picosecond,
            'stepSize': 2.0 * u.femtosecond  # dt = 2 fs
        }
        langevin_kwargs.update(self.langevin_kwargs)
        langevin_args = [langevin_kwargs[key] for key in ['temperature', 'frictionCoeff', 'stepSize']]
        self.integrator = mm.LangevinIntegrator(*langevin_args)

        system_kwargs = {
            'nonbondedMethod': app.PME,
            'nonbondedCutoff': 1.0 * u.nanometer,       # 1.2
            'constraints': app.HBonds,
            'rigidWater': True
            #'implicitSolvent': app.OBC2,
            #'rigidWater': False, 
            #'ewaldErrorTolerance': 0.0005
        }
        system_kwargs.update(self.system_kwargs)
        self.system = self.top.createSystem(**system_kwargs)


        # set force group
        self.force_group = {}
        for i in range(self.system.getNumForces()):
            force = self.system.getForce(i)
            m = re.search(r'OpenMM::(\S+)',str(force))
            name = m.group(1)
            if re.search(r'Bond',name):
                name = 'bond'
                grpnum = 0
            elif re.search(r'Angle',name):
                name = 'angle'
                grpnum = 1
            elif re.search(r'Torsion',name):
                name = 'torsion'
                grpnum = 2
            elif re.search(r'Nonbonded',name):
                name = 'nonbond'
                grpnum = 3
            elif re.search(r'GB',name):
                name = 'gb'
                grpnum = 4
            elif re.search(r'CMMotion',name):
                name = 'nonbond'
                grpnum = 3
            else:
                raise
            force.setForceGroup(grpnum)
            self.force_group[name] = grpnum

        if self.constant_pressure:
            self.system.addForce(mm.MonteCarloBarostat(1.0 * u.atmosphere, 300.00 * u.kelvin, 25))

        if self.stap:
            path_stap = f'/gpfs/deepfold/casp/refine/data/stap'
            self.add_stap_force(path_stap, stap_weight=self.stap_weight)

        if self.restrain_path:
            self.add_restrain_force()

    def _apply_restraints(self):
        for res in self.restraints:
            if res['type'] == 'position':
                Position(self.positions, self.top.topology).apply(self.system, **res['kwargs'])
            elif res['type'] == 'distance':
                Distance(self.positions, self.top.topology).apply(self.system, **res['kwargs'])
            elif res['type'] == 'symmetry':
                Symmetry(self.positions, self.top.topology).apply(self.system, **res['kwargs'])
            else:
                NotImplementedError('"{}" is not implemented restraint types.'.format(res['type']))

    def save_system(self, system_output, integrator_output):
        if any(map(lambda x: x is None, [self.system, self.integrator])):
            self._generate_system()
        self._apply_restraints()
        with open(system_output, 'w') as f:
            f.write(mm.XmlSerializer.serialize(self.system))
        with open(integrator_output, 'w') as f:
            f.write(mm.XmlSerializer.serialize(self.integrator))

    def set_force_group(self,fgrp_name):
        # if the name already exists
        fgrp_name = fgrp_name.lower()
        try:
            return self.force_group[fgrp_name]
        except:
            pass

        # determine a new force group number
        fgroup_num = max(self.force_group.values())+1

        self.force_group[fgrp_name] = fgroup_num
        return fgroup_num

    def add_stap_force(self, path_stap, stap_weight=0.5):
        tor_tables = []
        table_ind = {}
        aa_list = (
                'ALA', 'ARG', 'ASN', 'ASP', 'ASH', 'CYS', 'CYX', 'GLN', 'GLU', 'GLH',
                'GLY', 'HIS', 'HID', 'HIP', 'HIE', 'ILE', 'LEU', 'LYS', 'MET', 'PHE',
                'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
        )
        stap_list = ['phipsi', 'phichi1', 'psichi1', 'chi1chi2']
        atom_ind = {}
    
        residues = self.parmed_structure.residues
    
        stap_params = []
        for i, residue in enumerate(residues):
            resname = residue.name
            # convert AMBER-style residue names to standard ones
            if resname == 'ASH':
                resname = 'ASP'
            elif resname in ('HID', 'HIP', 'HIE'):
                resname = 'HIS'
            elif resname == 'CYX':
                resname = 'CYS'
            if resname not in aa_list:
                # not a protein residue
                continue
            (
                phi_inds, psi_inds, chi1_inds, chi2_inds
            ) = get_torsion_inds(residues, i)
            if 'phipsi' in stap_list and phi_inds is not None \
                    and psi_inds is not None:
                key = ( resname, 'phipsi' )
                if key not in table_ind:
                    tfile = f'{path_stap}/phipsi/{resname}.dat'
                    tor_tables.append(read_stap_table(tfile))
                    table_ind[key] = len(tor_tables)-1
                tind = table_ind[key]
                stap_params.append( (phi_inds, psi_inds, tind) )
            if 'phichi1' in stap_list and phi_inds is not None \
                    and chi1_inds is not None:
                key = ( resname, 'phichi1' )
                if key not in table_ind:
                    tfile = f'{path_stap}/phichi1/{resname}.dat'
                    tor_tables.append(read_stap_table(tfile))
                    table_ind[key] = len(tor_tables)-1
                tind = table_ind[key]
                stap_params.append( (phi_inds, chi1_inds, tind) )
            if 'psichi1' in stap_list and psi_inds is not None \
                    and chi1_inds is not None:
                key = ( resname, 'psichi1' )
                if key not in table_ind:
                    tfile = f'{path_stap}/psichi1/{resname}.dat'
                    tor_tables.append(read_stap_table(tfile))
                    table_ind[key] = len(tor_tables)-1
                tind = table_ind[key]
                stap_params.append( (psi_inds, chi1_inds, tind) )
            if 'chi1chi2' in stap_list and chi1_inds is not None \
                    and chi2_inds is not None:
                key = ( resname, 'chi1chi2' )
                if key not in table_ind:
                    tfile = f'{path_stap}/chi1chi2/{resname}.dat'
                    tor_tables.append(read_stap_table(tfile))
                    table_ind[key] = len(tor_tables)-1
                tind = table_ind[key]
                stap_params.append( (chi1_inds, chi2_inds, tind) )
    
        # If stap_params is empty, return immediately-skipping the setup.
        if not stap_params:
            return
    
        tor_tables.append(tor_tables[0]) # dummy because it should be periodic
        tor_tables = np.array(tor_tables) # type, 1st_tor, 2nd_tor (x, y, z)
        tor_tables *= stap_weight
        ntype, n_ang1, n_ang2 = tor_tables.shape
        # xsize, ysize, zsize, values, xmin, xmax, ymin, ymax, zmin, zmax
        pi = np.pi
        stap_table = mm.Continuous3DFunction(ntype, n_ang1, n_ang2,
                tor_tables.T.ravel(), 0.0, float(ntype)-1, -pi, pi, -pi, pi,
                periodic=True)
        force = mm.CustomCompoundBondForce(8,
                f'stap_func(itype,tor1,tor2);'
                f'tor1=dihedral(p1,p2,p3,p4);'
                f'tor2=dihedral(p5,p6,p7,p8);'
        )
        force.addPerBondParameter('itype')
        force.addTabulatedFunction('stap_func', stap_table)
    
        fgroup_num = self.system.getNumForces() + 1
    
        force.setForceGroup(fgroup_num)
        for atoms1, atoms2, itype in stap_params:
            force.addBond(atoms1 + atoms2, [float(itype)])
        self.system.addForce(force)

    def add_restrain_force(self):
        lines = []
        for rest_file in self.restrain_path:
            if not os.path.exists(rest_file):
                print (f'{rest_file} is not exist !!!! check it')
                continue
            f=open(rest_file, 'r')
            lines += f.readlines()
            f.close()

        forces = {}
        tab_dist_funcs = {}
        for line in lines:
            line = line.strip()
            if line == '' or line[0] == '#':
                continue
            c = line.split()
            if len(c) < 3:
                continue
            c[0] = c[0].upper()
            #NOTE: energy: kcal/mol, distance: A
            if c[0] == 'DISTANCE':
                # Flat-bottomed harmonic distance restraint
                # DISTANCE grpname atm1 atm2 force_const lower_d [ upper_d ]
                grpname = c[1]
                atm1 = int(c[2])-1
                atm2 = int(c[3])-1
                fconst = ( float(c[4]) * u.kilocalorie_per_mole / u.angstrom**2 ).value_in_unit(u.kilojoule_per_mole / u.nanometer**2)
                try:
                    weight = self.restrain_weights['distance'][grpname]
                    fconst *= weight
                except:
                    pass
                lower_d = float(c[5]) / 10.0 # nanometer
                try:
                    upper_d = float(c[6]) / 10.0 # nanometer
                except:
                    upper_d = lower_d
                assert( upper_d >= lower_d )
                fkey = ('DISTANCE',grpname)
                if fkey in forces:
                    force = forces[fkey]
                else:
                    force = mm.CustomBondForce('k*( step(-dr0)*dr0^2 + step(dr1)*dr1^2 );dr1=r-r1;dr0=r-r0')
                    forces[fkey] = force
                    force.addPerBondParameter('k')
                    force.addPerBondParameter('r0')
                    force.addPerBondParameter('r1')

                    grp = self.set_force_group(grpname)
                    force.setForceGroup(grp)

                force.addBond(atm1,atm2,(fconst,lower_d,upper_d))

            elif c[0] in ('LORENTZ', 'DISTANCE_LORENTZ'):
                # Flat-bottomed Lorentzian distance restraint
                # LORENTZ grpname atm1 atm2 fconst sigma lower_d [ upper_d ]
                grpname = c[1]
                atm1 = int(c[2])-1
                atm2 = int(c[3])-1
                fconst = ( float(c[4]) * u.kilocalorie_per_mole ).value_in_unit( u.kilojoule_per_mole )
                try:
                    weight = self.restrain_weights['distance'][grpname]
                    fconst *= weight
                except:
                    pass
                sigma = float(c[5]) / 10.0 # nanometer
                lower_d = float(c[6]) / 10.0 # nanometer
                try:
                    upper_d = float(c[7]) / 10.0 # nanometer
                except:
                    upper_d = lower_d
                assert( upper_d >= lower_d )
                fkey = ('LORENTZ',grpname)
                if fkey in forces:
                    force = forces[fkey]
                else:
                    force = mm.CustomBondForce('k*( step(-dr0)*(dr02/(dr02+sigma2)) + step(dr1)*(dr12/(dr12+sigma2)) );dr02=dr0^2;dr12=dr1^2;dr0=r-r0;dr1=r-r1')
                    forces[fkey] = force
                    force.addPerBondParameter('k')
                    force.addPerBondParameter('sigma2')
                    force.addPerBondParameter('r0')
                    force.addPerBondParameter('r1')
                    grp = self.set_force_group(grpname)
                    force.setForceGroup(grp)
                force.addBond(atm1,atm2,(fconst, sigma**2, lower_d,upper_d))

            elif c[0] == 'TABULAR_DISTANCE':
                # TABULAR_DISTANCE grpname xmin xmax yval1 yval2 yval3 ...
                grpname = c[1]
                atm1 = int(c[2])-1
                atm2 = int(c[3])-1
                xmin = float(c[4]) * 0.1 # x min (A -> nm)
                xmax = float(c[5]) * 0.1 # x max (A -> nm)
                yvals = ( np.array(list(map(float,c[6:]))) * u.kilocalorie_per_mole ).value_in_unit(u.kilojoule_per_mole)
                try:
                    weight = self.restrain_weights['distance'][grpname]
                    yvals *= weight
                except:
                    pass
                fkey = ('TABULAR_DISTANCE',grpname)
                self.set_force_group(grpname)
                if fkey in tab_dist_funcs:
                    tab_func = tab_dist_funcs[fkey]
                else:
                    tab_dist_funcs[fkey] = { 'xmin': xmin, 'xmax': xmax,
                            'yvals': [], 'atom_pairs': [] }
                    tab_func = tab_dist_funcs[fkey]
                if tab_func['xmin'] != xmin:
                    self.print_error("'xmin' does not match in the restraint file")
                    self.print_error(f'LINE: {line.rstrip()}')
                if tab_func['xmax'] != xmax:
                    self.print_error("'xmax' does not match in the restraint file")
                    self.print_error(f'LINE: {line.rstrip()}')
                if len(tab_func['yvals']) and len(yvals) != len(tab_func['yvals'][0]):
                    self.print_error("'yvals' size does not match in the restraint file")
                    self.print_error(f'LINE: {line.rstrip()}')
                tab_func['atom_pairs'].append((atm1,atm2))
                tab_func['yvals'].append(yvals)
            elif c[0] == 'ANGLE':
                # Flat-bottomed harmonic angle restraint
                # ANGLE grpname atm1 atm2 atm3 force_const lower_theta [ upper_theta ]
                grpname = c[1]
                atm1 = int(c[2])-1
                atm2 = int(c[3])-1
                atm3 = int(c[4])-1
                fconst = ( float(c[5]) * u.kilocalorie_per_mole ).value_in_unit(u.kilojoule_per_mole)
                try:
                    weight = self.restrain_weights['distance'][grpname]
                    fconst *= weight
                except:
                    pass
                lower_theta = float(c[6]) * np.pi / 180.0 # radian
                try:
                    upper_theta = float(c[7]) * np.pi / 180.0 # radian
                except:
                    upper_theta = lower_theta
                assert( lower_theta >= 0.0 and lower_theta <= np.pi )
                assert( upper_theta >= 0.0 and upper_theta <= np.pi )
                assert( upper_theta >= lower_theta )
                fkey = ('ANGLE',grpname)
                if fkey in forces:
                    force = forces[fkey]
                else:
                    force = mm.CustomAngleForce('k*( step(-dtheta0)*dtheta0^2 + step(dtheta1)*dtheta1^2 );dtheta1=theta-theta1;dtheta0=theta-theta0')
                    forces[fkey] = force
                    force.addPerAngleParameter('k')
                    force.addPerAngleParameter('theta0')
                    force.addPerAngleParameter('theta1')
                    grp = self.set_force_group(grpname)
                    force.setForceGroup(grp)
                force.addAngle(atm1,atm2,atm3,(fconst,lower_theta,upper_theta))
            elif c[0] == 'TORSION':
                # Flat-bottomed harmonic torsion angle restraint
                # TORSION grpname atm1 atm2 atm3 atm4 force_const lower_theta [ upper_theta ]
                grpname = c[1]
                atm1 = int(c[2])-1
                atm2 = int(c[3])-1
                atm3 = int(c[4])-1
                atm4 = int(c[5])-1
                fconst = ( float(c[6]) * u.kilocalorie_per_mole ).value_in_unit(u.kilojoule_per_mole)
                try:
                    weight = self.restrain_weights['torsion'][grpname]
                    fconst *= weight
                except:
                    pass
                lower_theta = float(c[7]) * np.pi / 180.0 # radian
                try:
                    upper_theta = float(c[8]) * np.pi / 180.0 # radian
                except:
                    upper_theta = lower_theta
                lower_theta = standard_radian(lower_theta)
                upper_theta = standard_radian(upper_theta)

                fkey = ('TORSION',grpname)
                if fkey in forces:
                    force = forces[fkey]
                else:
                    twopi = 2.0 * np.pi
                    force = mm.CustomTorsionForce(
                        'k*s*dtheta2;'
                        's = step(sign0 * sign1 * (theta1-theta0));'
                        'sign0=theta-theta0;'
                        'sign1=theta-theta1;'
                        'dtheta2=dtheta*dtheta;'
                        'dtheta = min(dtheta_low, dtheta_high);'
                        'dtheta_low = min('
                        '   abs(theta0-theta),'
                        f'   abs(theta0-(theta-{twopi}))'
                        ');'
                        'dtheta_high = min('
                        '   abs(theta-theta1),'
                        f'   abs((theta+{twopi})-theta1)'
                        ')'
                    )
                    forces[fkey] = force
                    force.addPerTorsionParameter('k')
                    force.addPerTorsionParameter('theta0')
                    force.addPerTorsionParameter('theta1')
                    grp = self.set_force_group(grpname)
                    force.setForceGroup(grp)
                force.addTorsion(atm1, atm2, atm3, atm4,
                        (fconst, lower_theta, upper_theta))
            elif c[0] == 'TORSION_LORENTZ':
                # Flat-bottomed lorentzian torsion angle restraint
                # TORSION_LORENTZ grpname atm1 atm2 atm3 amt4 fconst
                #         sigma_degree lower_degree [ upper_degree ]
                grpname = c[1]
                atm1 = int(c[2])-1
                atm2 = int(c[3])-1
                atm3 = int(c[4])-1
                atm4 = int(c[5])-1
                fconst = ( float(c[6]) * u.kilocalorie_per_mole
                        ).value_in_unit( u.kilojoule_per_mole )
                try:
                    weight = self.restrain_weights['torsion'][grpname]
                    fconst *= weight
                except:
                    pass
                sigma = float(c[7]) * np.pi / 180.0 # radian
                lower_theta = float(c[8]) * np.pi / 180.0 # radian
                try:
                    upper_theta = float(c[9]) * np.pi / 180.0 # radian
                except:
                    upper_theta = lower_theta
                lower_theta = standard_radian(lower_theta)
                upper_theta = standard_radian(upper_theta)
                fkey = ('TORSION_LORENTZ',grpname)
                if fkey in forces:
                    force = forces[fkey]
                else:
                    twopi = 2.0 * np.pi
                    force = mm.CustomTorsionForce(
                        'k*s*(dtheta2/(dtheta2+sigma2));'
                        's = step(sign0 * sign1 * (theta1-theta0));'
                        'sign0=theta-theta0;'
                        'sign1=theta-theta1;'
                        'dtheta2=dtheta*dtheta;'
                        'dtheta = min(dtheta_low, dtheta_high);'
                        'dtheta_low = min('
                        '   abs(theta0-theta),'
                        f'   abs(theta0-(theta-{twopi}))'
                        ');'
                        'dtheta_high = min('
                        '   abs(theta-theta1),'
                        f'   abs((theta+{twopi})-theta1)'
                        ')'
                    )
                    forces[fkey] = force
                    force.addPerTorsionParameter('k')
                    force.addPerTorsionParameter('sigma2')
                    force.addPerTorsionParameter('theta0')
                    force.addPerTorsionParameter('theta1')
                    grp = self.set_force_group(grpname)
                    force.setForceGroup(grp)
                force.addTorsion(atm1, atm2, atm3, atm4,
                        (fconst, sigma**2, lower_theta, upper_theta))
            elif c[0] == 'SYMMETRY_DISTANCE':
                # symmetry distance restraint
                # SYMMETRY_DISTANCE grpname atm1 atm2 atm3 atm4 fconst
                grpname = c[1]
                atm1 = int(c[2])-1
                atm2 = int(c[3])-1
                atm3 = int(c[4])-1
                atm4 = int(c[5])-1
                fconst = ( float(c[6]) * u.kilocalorie_per_mole /
                        u.angstrom**2 ).value_in_unit(
                                u.kilojoule_per_mole / u.nanometer**2 )
                try:
                    weight = self.restrain_weights['symmetry'][grpname]
                    fconst *= weight
                except:
                    pass
                fkey = ('SYMMTRY_DISTANCE',grpname)
                if fkey in forces:
                    force = forces[fkey]
                else:
                    force = mm.CustomCompoundBondForce(4,
                            'k*dd2;dd2=(distance(p1,p2) - distance(p3,p4))^2')
                    forces[fkey] = force
                    force.addPerBondParameter('k')
                    grp = self.set_force_group(grpname)
                    force.setForceGroup(grp)
                force.addBond((atm1, atm2, atm3, atm4), [fconst])
            else:
                raise ValueError('unknown restraint %s'%c[0])

        # add forces into system
        for force in forces.values():
            self.system.addForce(force)

        # tabular distance function
        for i, fkey in enumerate(tab_dist_funcs):
            if fkey[0] == 'TABULAR_DISTANCE':
                grpname = fkey[1]
                tab = tab_dist_funcs[fkey]
                xmin = tab['xmin']
                xmax = tab['xmax']
                ytable = np.array(tab['yvals'])
                atom_pairs = tab['atom_pairs']
                grp = self.set_force_group(grpname)

                ntype, npoint = ytable.shape

                tab_function = mm.Continuous2DFunction(ntype, npoint,
                        ytable.T.ravel(), 0.0, float(ntype-1), xmin, xmax)

                fname = f'dist_2dfunc{i}'
                force = mm.CustomCompoundBondForce(2,
                        f'{fname}(itype,d);'
                        f'd=min({xmax},max({xmin},distance(p1,p2)));'
                )
                force.addPerBondParameter('itype')
                force.addTabulatedFunction(fname, tab_function)
                force.setForceGroup(grp)
                for itype in range(ntype):
                    atom_pair = atom_pairs[itype]
                    force.addBond(atom_pair, [float(itype)])
                self.system.addForce(force)

class GenerateState(Constructor):
    def __init__(self,
                 prmtop: Union[os.PathLike[str], str],
                 inpcrd: Union[os.PathLike[str], str]):
        super().__init__(prmtop, inpcrd)

    def _generate_simulation(self):
        self._platform = mm.Platform_getPlatformByName('CPU')
        self._simulation = app.Simulation(
            self.top.topology,
            self.system,
            self.integrator,
            self._platform
        )
        self._simulation.context.setPositions(self.crd.getPositions())

    def save_state(self, output: Union[os.PathLike[str], str]):
        self._generate_system()
        self._generate_simulation()
        self._simulation.context.setVelocitiesToTemperature(300.00 * u.kelvin)
        self._simulation.saveState(output)


class MDBase:
    """
    Base class for refinement recipes.
    """

    def __init__(self, **kwargs):
        """
        prefix: str
        state_output: str

        prmtop: str
        state_input: str
        system_input: str
        integrator_input: str

        platform_name: str
        properties: Dict[str, str]

        logger: logging.Logger
        """
        # Check kwargs is complete.
        # `logger` is provided from `pipeline`.

        self.__dict__.update(kwargs)

        self.log(f'>> "{self.prefix}"')

        self._load_topology()
        # self._setup_platform()
        # self._setup_simulation()

    def _setup_platform(self):
        self._platform = mm.Platform_getPlatformByName(self.platform_name)

    def _setup_simulation(self):
        self._simulation = app.Simulation(
            topology=self._topology,
            system=self.system_input,
            integrator=self.integrator_input,
            platform=self._platform,
            platformProperties=self.properties,
            # state=self.state_input
        )
        self._simulation.loadState(self.state_input)

    def _setup_reporters(self, reporters):
        self._simulation.reporters.clear()
        self._simulation.reporters.extend(reporters)

    def _save_state(self):
        self.log('   Save a state file...')
        self._simulation.saveState(self.state_output)

    def _load_topology(self):
        top = app.AmberPrmtopFile(self.prmtop)
        self._topology = top.topology

    def _save_pdb(self):
        self.log('   Save a PDB file...')
        app.PDBFile.writeFile(self._topology,
                              self._simulation.context.getState(getPositions=True).getPositions(),
                              open(f'{self.prefix}.pdb', 'w'),
                              keepIds=True
                            )

    def log(self, msg: str, level: int = logging.INFO):
        self.logger.log(level, msg)
        for x in self.logger.handlers:
            x.flush()

    def _setup(self):
        self.log('   Prepare a platform...')
        self._setup_platform()
        self.log('   Prepare a simulation...')
        self._setup_simulation()


class Minimize(MDBase):
    """
    Run energy minimization.
    """

    def __init__(self,
                 prefix,
                 prmtop,
                 state_output,
                 state_input, system_input, integrator_input,
                 platform_name, properties,
                 logger, **kwargs):
        base_kwargs = {
            'prefix': prefix,
            'prmtop': prmtop,
            'state_output': state_output,
            'state_input': state_input,
            'system_input': system_input,
            'integrator_input': integrator_input,
            'platform_name': platform_name,
            'properties': properties,
            'logger': logger
        }
        super().__init__(**base_kwargs)
        self.__dict__.update(kwargs)

        self._setup()
        reporters = [
            app.StateDataReporter(file=sys.stderr, reportInterval=100, step=True, potentialEnergy=True,
                                  elapsedTime=True, separator='\t')
        ]
        self._setup_reporters(reporters)

    def run(self, max_iter: int = 0):
        self.log('   Minimize energy...')
        self._simulation.minimizeEnergy(maxIterations=max_iter)
        self._save_state()
        self._save_pdb()


class RunMD(MDBase):
    """
    Run molecular dynamics simulation.
    """

    def __init__(self,
                 prefix,
                 prmtop,
                 state_output,
                 state_input, system_input, integrator_input,
                 platform_name, properties,
                 traj_out: Union[os.PathLike[str], str], nstxout: int,
                 logger, nstlog: int,
                 gen_vel_to: Optional[float] = None, **kwargs):
        base_kwargs = {
            'prefix': prefix,
            'prmtop': prmtop,
            'state_output': state_output,
            'state_input': state_input,
            'system_input': system_input,
            'integrator_input': integrator_input,
            'platform_name': platform_name,
            'properties': properties,
            'logger': logger
        }
        super().__init__(**base_kwargs)
        self.traj_out = traj_out,
        self.nstxout = int(nstxout)
        self.nstlog = int(nstlog)

        self.gen_vel_to = gen_vel_to
        self.__dict__.update(kwargs)

        self._setup()

        """ Below is change I made  """
#        atom_indices = [atom.index for atom in list(self.top.topology.residues())[10].atoms()] #This

        reporters = [
            mdtraj.reporters.XTCReporter(file=traj_out, reportInterval=self.nstxout, append=False),
#            mdtraj.reporters.XTCReporter(file=traj_out, reportInterval=self.nstxout, atomSubset = 'protein', append=False),
#            mdtraj.reporters.append(StateDataReporter(stdout, 1000, step=True, potentiaEnergy=True, temperature=True)),
#            mdtraj.reporters.PDBReporter('pdb_out.pdb', reportInterval=5000, atom_indices=atom_indices), #and This
            LoggerReporter(self.logger, report_interval=self.nstlog, flush=True)
        ]
        self._setup_reporters(reporters)

        if self.gen_vel_to is not None:
            self._simulation.context.setVelocitiesToTemperature(gen_vel_to)

    def step(self, steps: int):
        self.log(f"   Advance the simulation for {steps} steps...")
        self._simulation.step(steps)
        self._save_state()
        self._save_pdb()

    def run_for(self, maxh: float):
        self.log("")
        self.log(f"   Advance the simulation for {maxh} hours...")
        self._simulation.runForClockTime(maxh * 0.95 * u.hour, None, self.state_output, 10.0 * u.minute)
        self._save_state()
        self._save_pdb()
