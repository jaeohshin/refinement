import re
import os
from typing import Union, Optional, List
import sys

from casp16.utils import use_temp, cwd, run

__all__ = ['call_tleap']

# Load FF
v6_str = """# Load force fields
source leaprc.protein.ff19SB
source leaprc.water.tip3p
"""

# Load PDB
load_pdb_str = """# Load model
{model} = loadPdb {pdb_in}
# Info on system charge
charge {model}
"""

# Add solvent, counter-ions # {{ ... }}
tip3p_solvate_str = """# Add solvent and counterions
solvateBox {model} TIP3PBOX {{ 10.0 10.0 10.0 }}
addIons {model} Na+ 0
addIons {model} Cl- 0
"""

# Describe the configuration
desc_str = """# Describe about the configuration
desc {model}
"""

# Save the configuration
save_pdb = """# Save the configuration
savePdb {model} {pdb_out}
"""

# Add ions and box
add_ions = """# Add ions
addIons {model} Na+ {n_ions}
addIons {model} Cl- {n_ions}
# Info on system charge
charge {model}
# Set box
# setBox {model} "centers"
"""

# Save the parameters
save_parm_str = """# Save the parameters
saveAmberParm {model} {prmtop} {inpcrd}
savePdb {model} {pdb_out}
"""

# Exit
exit_str = """# Exit
quit
"""


def call_tleap(
        pdb_input: Union[os.PathLike[str], str],
        pdb_output: Union[os.PathLike[str], str],
        prmtop: Union[os.PathLike[str], str],
        inpcrd: Union[os.PathLike[str], str],
        tleap_path: Optional[Union[os.PathLike[str], str]] = None,
        include_dir: Optional[Union[os.PathLike[str], str]] = None,
        ss_cmds: Optional[List[str]] = None
):
    """No hydrogen!!!"""
    pdb_input = os.path.abspath(pdb_input)
    pdb_output = os.path.abspath(pdb_output)
    prmtop = os.path.abspath(prmtop)
    inpcrd = os.path.abspath(inpcrd)

    if include_dir is not None:
        include_dir = os.path.abspath(include_dir)

    model_name = 'protein'

    leap_cmd_1 = []
    leap_cmd_2 = []

    if tleap_path is not None:
        tleap_path = os.path.abspath(tleap_path)
    else:
        tleap_path = 'tleap'

    if include_dir is not None:
        s = f'addPath "{include_dir}"'
        leap_cmd_1.append(s)
        leap_cmd_2.append(s)

    leap_cmd_1.append(v6_str)
    leap_cmd_2.append(v6_str)

    leap_cmd_1.append(load_pdb_str.format(pdb_in=pdb_input, model=model_name))
    leap_cmd_1.append(tip3p_solvate_str.format(model=model_name))
    leap_cmd_1.append(exit_str)

    with use_temp():
        # with cwd('.'):
        tleap_input = '\n'.join(leap_cmd_1)
        with open('leap.in', 'w') as f:
            f.write(tleap_input)
        leap_out = run([tleap_path, '-f', 'leap.in'])

        print(leap_out, file=sys.stderr)

        p = re.compile(r'Volume: (.*) A\^3')
        m = p.findall(leap_out)
        n_wat = int(float(m[0].strip()))
        n_ions = int(6.02e-4 * n_wat)  # 0.15 M

        leap_cmd_2.append(load_pdb_str.format(pdb_in=pdb_input, model=model_name))
        leap_cmd_2.append(desc_str.format(model=model_name))
        for s in ss_cmds:
            leap_cmd_2.append(s.format(model=model_name))
        leap_cmd_2.append(tip3p_solvate_str.format(model=model_name))
        leap_cmd_2.append(add_ions.format(model=model_name, n_ions=n_ions))
        leap_cmd_2.append(save_parm_str.format(model=model_name, pdb_out=pdb_output, prmtop=prmtop, inpcrd=inpcrd))
        leap_cmd_2.append(exit_str)

        tleap_input = '\n'.join(leap_cmd_2)
        with open('leap.in', 'w') as f:
            f.write(tleap_input)
        leap_out = run([tleap_path, '-f', 'leap.in'])

        print(leap_out, file=sys.stderr)
