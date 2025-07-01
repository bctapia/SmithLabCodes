"""
This script has been written to run grand canonical Monte Carlo simulations
The python package pysimm (https://pysimm.org/) is required to run this code
Pysimm has been modified for parallel optimization and is available here: https://github.com/bctapia/pysimm
Extensive parts of this code have been taken from pysimm examples
Written by: Brandon C. Tapia (bctapia@mit.edu) while in the Smith Research Lab at the Massachusetts Institute of Technology
The MIT License is in effect
Please cite this code if used
"""

from pysimm import cassandra
from pysimm import system
from os import path as osp
import numpy
import pandas as pd


def run(test=True):

    # Gas names as they will be referred through simulations, any number of gasses allowed
    gas_names = ["Gas1", "Gas2", "etc..."]

    # The chemical potential is related to the system pressure with this equation
    # See the 'chemical potential to pressure' base code to develop your own calibrations
    chem_pots = [
        lambda x: 2.3112 * numpy.log(x) - 38.8328
    ]  # dependent on temperature and gas inserted

    # Root directory for gas data.
    data_dir = osp.join(
        ".", "gases"
    )  # gases must be in a folder labeled gases, one down from main folder

    # Setup of adsorbate model
    gases = []
    for gn in gas_names:
        gases.append(system.read_lammps(osp.join(data_dir, gn + ".lmps")))
        gases[-1].forcefield = "trappe/amber"  # change depending on gas forcefields you are using

    # Setup of adsorbent model
    frame = system.read_lammps(
        "POLYMER_SYSTEM.lmps"
    )  # update based on your polymer .lmps file name
    frame.forcefield = "trappe/amber"  # change depending on polymer forcefields you are using

    # Constant for loadings calculations
    molec2mmols_g = 1e3 / frame.mass  # converting molecules -> mmol/g

    print("The conversion is:", molec2mmols_g)

    # Setup of the GCMC simulations
    css = cassandra.Cassandra(frame)
    sim_settings = css.read_input("run_props.inp")  # all Cassandra setup is in this run_props file

    # This function sets up the GCMC simulation but does NOT perform it
    def calculate_isotherm_point(gas_name, press):
        run_fldr = osp.join(gas_name, str(press))
        idx = gas_names.index(gas_name)
        css.add_gcmc(
            species=gases[idx],
            is_new=True,
            chem_pot=chem_pots[idx](press),
            out_folder=run_fldr,
            props_file="gcmc.inp",
            **sim_settings,
        )
        css.run()
        full_prp = css.run_queue[0].get_prp()
        step_return = full_prp[0]
        mols_return = full_prp[3]
        loading_return = molec2mmols_g * numpy.average(full_prp[3][int(len(2 * full_prp[3]) / 3) :])
        return [step_return, mols_return, loading_return]

    loadings = []
    mols = []

    # Calculation of adsorption isotherms
    gas_press = numpy.array(
        [0.5, 1, 1.5, 2]
    )  # bar # update prssure based on points you want calculated

    loadings = dict.fromkeys(gas_names)
    for gn in gas_names:
        loadings[gn] = []
        for p in gas_press:
            full_data = calculate_isotherm_point(gn, p)  # now actually running the simulation
            # Extracting important data from the full dataset
            loading_data = full_data[2]
            mols_data = full_data[1]
            step_data = full_data[0]
            loadings[gn].append(loading_data)
            mols.append(mols_data[:, numpy.newaxis])
    # code to print important data into .txt files
    mols = numpy.hstack(mols)
    mols = pd.DataFrame(mols)
    mols.columns = numpy.tile(gas_press, len(gas_names))
    mols = mols.add_prefix("P=")
    steps = pd.DataFrame(step_data)
    mols.insert(0, "Steps", steps)
    mols.to_string("mols_1.txt", index=False)  # mols is the moles inserted for each P

    loadings = pd.DataFrame(loadings)  # loadings is the equilibrium mols transferred into mmol/g
    pressure = pd.DataFrame(gas_press)  # pressure restates the pressures inputted
    loadings.insert(0, "Pressure (bar)", pressure)
    loadings.to_string("loadings_1.txt", index=False)


if __name__ == "__main__":
    run()
