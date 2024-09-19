"""
Aug. 19, 2024. Jaeoh Shin.
I will simplify the code for kinase project.
- make the minimization process two steps (instead of four)
- make the equilibration process one step (instead of two)
"""

import logging
import sys
from itertools import combinations
import xml.etree.ElementTree as ET
import argparse

import openmm as mm
from openmm import unit as u
import numpy as np

import parmed
from Bio.PDB.PDBParser import PDBParser

sys.path.append("/gpfs/deepfold/users/jaeohshin2/mdsampling/")

import casp16

logging.getLogger().setLevel(logging.INFO)
logging.getLogger().addHandler(logging.StreamHandler(stream=sys.stderr))


# weights for extra forces

STAP_WEIGHT = 0.0

def build_distance_restraint_v7(input_pdb):
    parser = PDBParser()
    struct = parser.get_structure("", input_pdb)
    model = next(struct.get_models())
    abl = []
    for residue in model.get_residues():
        for atom in filter(lambda a: a.get_id() == 'CA', residue.get_atoms()):
            abl.append({'coord': atom.get_coord(), 'plddt': atom.get_bfactor()})
    good_residue = []
    for i, v in enumerate(abl):
        if v['plddt'] < 90.0:
            continue
        good_residue.append(i)
    dist_dict = {(i, j): np.linalg.norm(abl[i]['coord'] - abl[j]['coord']) for i, j in combinations(good_residue, 2)}

    k0 = 0.0
    cutoff = 1.0
    params = {(i, j): (r, k0 / r ** 2) for (i, j), r in dist_dict.items() if r <= cutoff}
    return params

""" Prepare the system """
def initialize():
    print('\nRun modeller\n', file=sys.stderr)

    # Modeller
    casp16.modeller.gen_init(
        'prep.pdb',
        'init.pdb',
        'init.prmtop',
        'init.rst',
        #include_dir='/gpfs/deepfold/casp/refine/cufix',
        ss_kwargs={
            'cutoff': 3.0
        },
        ssbond = 'SSBOND',
        tleap_path='/user/deepfold/anaconda3/envs/refine/bin/tleap'
    )
    print('\nGenerate an initial state\n', file=sys.stderr)
    gen = casp16.md.GenerateState(
        'init.prmtop',
        'init.rst'
        #'model.prmtop',
        #'model.inpcrd'
    )

    gen.save_state('init.state.xml')
    del gen


def minimize(ngpus: int = 1):
    print('\nminim.1\n', file=sys.stderr)


    gen = casp16.md.Constructor(
        'init.prmtop', 'init.rst',
        stap = True,
        stap_weight = STAP_WEIGHT,
        restraints=[
            {'type': 'position', 'kwargs': {
                'potential': 'harmonic',
                'mode': 'heavy',
                'weight': 100.0,
                'positions': None}}
        ]
    )
    gen.save_system('minim.1.system.xml', 'minim.1.integrator.xml')
    del gen

    handle = casp16.md.Minimize(
        'minim.1',
        'init.prmtop',
        'minim.1.state.xml',
        'init.state.xml', 'minim.1.system.xml', 'minim.1.integrator.xml',
        platform_name='CUDA',
        properties={'CudaPrecision': 'mixed', 'DeviceIndex': ','.join(str(i) for i in range(ngpus))},
        logger=logging.getLogger(), nstlog=10,
    )
    handle.run(max_iter=800)
    del handle

"""
Deleted three minimize() functions.
"""

def equilibrate(ngpus: int = 1):
#    ca_dict = build_distance_restraint_v7('prep.pdb')

    print('\nequil.1\n', file=sys.stderr)
    gen = casp16.md.Constructor(
        'init.prmtop', 'init.rst',
        stap = False,
        stap_weight = STAP_WEIGHT,
        restraints=[
            {'type': 'position', 'kwargs': {
                'potential': 'harmonic',
                'mode': 'backbone',
                'weight': 5.0,
                'positions': None}}
#            ,
#            {'type': 'distance', 'kwargs': {
#                'potential': 'flat_bottom',
#                'r_cut': 2.0,
#                'parameters': ca_dict}}
        ]
#        ,
#        restrain_forces= RESTRAINTS
    )
    gen.save_system('equil.1.system.xml', 'equil.1.integrator.xml')
    del gen

    handle = casp16.md.RunMD(
        'equil.1',
        'init.prmtop',
        'equil.1.state.xml',
        'minim.1.state.xml', 'equil.1.system.xml', 'equil.1.integrator.xml',
        # 'CUDA', {'CudaPrecision': 'mixed', 'DeviceIndex': '0'},
        platform_name='CUDA',
        properties={'CudaPrecision': 'mixed', 'DeviceIndex': ','.join(str(i) for i in range(ngpus))},
        traj_out='equil.1.xtc', nstxout=5000, logger=logging.getLogger(), nstlog=10000,
        gen_vel_to=300.00 * u.kelvin
    )
    handle.step(steps=10000)  
    del handle

"""
Delete the 2nd equilibration
"""


def production(ngpus: int = 1):
#    ca_dict = build_distance_restraint_v7('prep.pdb')

    tree = ET.parse('equil.1.state.xml')
    root = tree.getroot()
    se = root.find('Parameters')
    root.remove(se)
    tree.write('prod.0.state.xml')

    print('\nprod.1\n', file=sys.stderr)
    gen = casp16.md.Constructor(
        'init.prmtop', 'init.rst',
        stap = True,
        stap_weight = STAP_WEIGHT,

        restraints=[
#            {'type': 'position', 'kwargs': {
#                'potential': 'flat_bottom',
#                'mode': 'ca',
#                'weight': 0.0, 
#                'positions': None}}
#            ,
#            {'type': 'distance', 'kwargs': {
#                'potential': 'flat_bottom',
#                'r_cut': 2.0,
#                'parameters': ca_dict}}
        ]
#        ,
#       restrain_forces= RESTRAINTS
    )
    gen.save_system('prod.1.system.xml', 'prod.1.integrator.xml')
    del gen

    handle = casp16.md.RunMD(
        'prod.1',
        'init.prmtop',
        'prod.1.state.xml',
        'prod.0.state.xml', 'prod.1.system.xml', 'prod.1.integrator.xml',
        # 'CUDA', {'CudaPrecision': 'mixed', 'DeviceIndex': '0'},
        platform_name='CUDA',
        properties={'CudaPrecision': 'mixed', 'DeviceIndex': ','.join(str(i) for i in range(ngpus))},
        traj_out='prod.1.xtc', nstxout=5000, logger=logging.getLogger(), nstlog=10000
    )
    handle.step(steps=25000000)  # 50 ns
    handle.step(steps=25000000)  # 50 ns
    handle.step(steps=25000000)  # 50 ns
    handle.step(steps=25000000)  # 50 ns
    del handle


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gpus', type=int, default=1)
    args = parser.parse_args()
    ngpus = args.gpus

    initialize()
    minimize(ngpus)
    equilibrate(ngpus)
    production(ngpus)

if __name__ == '__main__':
    main()

