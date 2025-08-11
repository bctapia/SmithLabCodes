""" """

from pathlib import Path
import shutil
from smithlab import slurm


def write_mcmd(
    gas_path,
    pol_path,
    temp,
    trial,
    press,
    chem_pot,
    iters,
    np,
    mc_props=None,
    md_props=None,
    forcefield="trappe/amber",
):
    """Create all the necessary directories and files for running an MCMD simulation. The structure of the directories looks like:
    ~/MCMD/Polymer/Gas/Temperature/Trial/Pressure
    For example: ~/MCMD/CANAL-Me-Me2F-CN/C3H8/308/Trial1/5
    """
    gas_name = Path(gas_path).name
    pol_name = Path(pol_path).name

    gas_stem = Path(gas_path).stem
    pol_stem = Path(pol_path).stem

    pol_stem_split = pol_stem.split("_")[0]

    folder_path = (
        Path.home()
        / "MCMD"
        / pol_stem_split
        / str(gas_stem).upper()
        / f"{temp}"
        / f"Trial{trial}"
        / f"{press}"
    )

    folder_path.mkdir(parents=True, exist_ok=True)

    # copy all necessary files
    shutil.copy(gas_path, folder_path)
    shutil.copy(pol_path, folder_path)

    with open(folder_path / "run_MCMD.py", "w", encoding="utf-8") as f:

        # header
        f.writelines(
            ["from pysimm.apps import mc_md\n", "from pysimm import system\n", "import numpy\n\n"]
        )

        # general data
        f.write(f"chem_pot = lambda x: {chem_pot[0]}*numpy.log(x)-{chem_pot[1]}\n\n")
        f.write(f"press = {press}  # bar\n")
        f.write(f"temp = {temp}  # K\n")
        f.write(f"iters = {iters}\n")
        f.write(f"polymer_name = '{pol_name}'\n")
        f.write(f"gas_name = '{gas_name}'\n")
        f.write(f"np_val = {np}\n")
        f.write("rst = False\n\n\n")

        f.write("def run(test=False):\n")
        f.write("    frame = system.read_lammps(polymer_name)\n")
        f.write(f"    frame.forcefield = '{forcefield}'\n\n")

        f.write("    gas1 = system.read_lammps(gas_name)\n")
        f.write(f"    frame.forcefield = '{forcefield}'\n\n")

        if mc_props is None:
            mc_props = {}

        mc_props_default = {
            "rigid_type": True,  # only FALSE if modeled as single point (e.g UA CH4)
            "max_ins": 20000,
            "Rcutoff_Low": 0.4,
            "Run_Type": {"steps": 100},
            "CBMC_Info": {"rcut_cbmc": 2.0},
            "Simulation_Length_Info": {"run": 2000000, "coord_freq": 25000, "prop_freq": 100},
            "VDW_Style": {"cut_val": 15.0},
            "Charge_Style": {"cut_val": 15.0},
            "Property_Info": {
                "prop1": "energy_total",
                "prop2": "pressure",
                "prop3": "nmols",
                "prop4": "density",
            },
        }

        def deep_update(default: dict, override: dict) -> dict:
            for key, value in override.items():
                if key in default and isinstance(default[key], dict) and isinstance(value, dict):
                    deep_update(default[key], value)
                else:
                    default[key] = value
            return default

        mc_props_full = deep_update(mc_props_default.copy(), mc_props or {})
        g = lambda k: mc_props_full[k]

        f.write("    mc_props = {\n")
        f.write(
            f"        'rigid_type': {g('rigid_type')},  # only FALSE if modeled as single point (e.g UA CH4)\n"
        )
        f.write(f"        'max_ins': {g('max_ins')},\n")
        f.write(f"        'Chemical_Potential_Info': chem_pot(press),\n")
        f.write(f"        'Temperature_Info': temp,\n")
        f.write(f"        'Rcutoff_Low': {g('Rcutoff_Low')},\n")
        f.write(f"        'Run_Type': {{\n")
        f.write(f"            'steps': {g('Run_Type')['steps']}\n")
        f.write(f"            }},\n")
        f.write(f"        'CBMC_Info': {{\n")
        f.write(f"            'rcut_cbmc': {g('CBMC_Info')['rcut_cbmc']}\n")
        f.write(f"            }},\n")
        f.write(f"        'Simulation_Length_Info': {{\n")
        f.write(f"            'run': {g('Simulation_Length_Info')['run']},\n")
        f.write(f"            'coord_freq': {g('Simulation_Length_Info')['coord_freq']},\n")
        f.write(f"            'prop_freq': {g('Simulation_Length_Info')['prop_freq']}\n")
        f.write(f"            }},\n")
        f.write(f"        'VDW_Style': {{\n")
        f.write(f"            'cut_val': {g('VDW_Style')['cut_val']}\n")
        f.write(f"            }},\n")
        f.write(f"        'Charge_Style': {{\n")
        f.write(f"            'cut_val': {g('Charge_Style')['cut_val']}\n")
        f.write(f"            }},\n")
        f.write(f"        'Property_Info': {{\n")
        f.write(f"            'prop1': '{g('Property_Info')['prop1']}',\n")
        f.write(f"            'prop2': '{g('Property_Info')['prop2']}',\n")
        f.write(f"            'prop3': '{g('Property_Info')['prop3']}',\n")
        f.write(f"            'prop4': '{g('Property_Info')['prop4']}'\n")
        f.write(f"            }}\n")
        f.write("        }\n")

        if md_props is None:
            md_props = {}

        md_props_default = {
            "pressure": {"iso": "iso"},
            "kspace_style": 1e-5,
            "timestep": 0.5,
            "cutoff": 15.0,
            "length": 2000000,
            "thermo": 25000,
            "dump": 25000,
            "print_to_screen": False,
        }

        md_props_full = deep_update(md_props_default.copy(), md_props or {})
        g = lambda k: md_props_full[k]

        f.write("    md_props = {\n")
        f.write(f"        'temp': temp,\n")
        f.write(f"        'pressure': {{\n")
        f.write(f"            'start': press,\n")
        f.write(f"            'iso': '{g('pressure')['iso']}'\n")
        f.write(f"            }},\n")
        f.write(f"        'kspace_style': {g('kspace_style')},\n")
        f.write(f"        'timestep': {g('timestep')},\n")
        f.write(f"        'cutoff': {g('cutoff')},\n")
        f.write(f"        'length': {g('length')},\n")
        f.write(f"        'thermo': {g('thermo')},\n")
        f.write(f"        'dump': {g('dump')},\n")
        f.write(f"        'print_to_screen': {g('print_to_screen')}\n")
        f.write("        }\n")

        f.write(
            "    sim_result = mc_md.mc_md(gas1, frame, mcmd_niter=iters, sim_folder='results', mc_props=mc_props, md_props=md_props, np=np_val, prefix='mpiexec', restart=rst)\n"
        )
        f.write("    sim_result.write_lammps('SYSplusGAS.lmps')\n")
        f.write("    sim_result.write_xyz('SYSplusGAS.xyz')\n\n\n")

        f.write("if __name__ == '__main__':\n")
        f.write("    run(False)\n")

    # write batch script
    commands = """source /etc/profile
module load anaconda/2023a-tensorflow
module load intel/oneapi/compiler/latest
module load intel/oneapi/mpi/latest
module load intel/oneapi/mkl/latest

eval "$(conda shell.bash hook)"
conda activate lammps_test_env

# Set environment variables
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=spread

python run_MCMD.py
"""
    if np == 40:
        partition = "xeon-g6-volta"
    else:
        partition = "xeon-p8"

    slurm.write_batch(
        folder_path / "flex.sh", N=1, n=1, partition=partition, cpus_per_task=np, command=commands
    )
