from typing import Union
from os import PathLike
import Bio.PDB.Model

import warnings
import Bio.PDB.Atom

# Disable PDBConstructionWarning
warnings.filterwarnings('ignore', category=Bio.PDB.Atom.PDBConstructionWarning)


def sort_model(
        target_model: Bio.PDB.Model.Model,
        template_model: Bio.PDB.Model.Model):
    for r1, r2 in zip(target_model.get_residues(), template_model.get_residues()):
        # print(r1)
        atoms = {x.full_id[4][0]: x for x in r2.get_atoms()}
        for a in r1.get_atoms():
            a.set_coord(atoms[a.id].get_coord())
            a.set_bfactor(atoms[a.id].get_bfactor())
            a.set_occupancy(atoms[a.id].get_occupancy())
            # a.set_occupancy(1.0)
            # print(a.full_id, '=>', a.get_coord())


def get_atoms(input_pdb: Union[PathLike[str], str]) -> str:
    with open(input_pdb, 'r') as f:
        lines = f.readlines()
    atoms = []
    for line in lines:
        if line.startswith('ATOM'):
            atoms.append(line)
        elif line.startswith('TER'):
            atoms.append('TER')
    return ''.join(atoms)
