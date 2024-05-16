'''
This script runs multi-gas grand canonical Monte Carlo (GCMC) simulations
required for usage of the MC/MD or MC codes
The python package pysimm (https://pysimm.org/) is required to run this code
Pysimm has been modified for parallel optimization and is available here: https://github.com/bctapia/pysimm
Extensive parts of this code have been taken from pysimm examples
Author: Brandon C. Tapia (bctapia@mit.edu) while in the Smith Research Lab at the Massachusetts Institute of Technology
The MIT License is in effect
Please cite this code if used
'''

from pysimm import system
from pysimm import cassandra
from os import path as osp
import numpy
import pandas as pd

def run(test=True):

    # gas names as they will be referred through simulations 
    gas_names = ['Gas1', 'Gas2'] # update based on your gas names

    # chemical potentials as a function of pressure
    def chempots(x): # see the calibration_curve_SHELL.py script to generate chemical potentials
        species1 = 2.3112*numpy.log(x)-38.8328 # dependent on temperature and gas inserted
        species2 = 2.5456*numpy.log(x)-41.2463 # dependent on temperature and gas inserted
        return [species1, species2]

    # if molecule is more than one sphere rigid MUST = True
    isrigid = [False, True]

     # update pressure based on points you want calculated
    gas_press = [5, 10, 15] # bar

    # max number of molecules inserted (adjust as needed)
    maxins = [1000,1000] # larger values require more memory

    # root directory for adsorbate files
    data_dir = osp.join('.', 'gases') # you need to place the gases .lmps in './gases'

    # setup of adsorbate models
    gases = [] 
    for gn in gas_names:
        gases.append(system.read_lammps(osp.join(data_dir, gn + '.lmps')))
        gases[-1].forcefield = 'trappe/amber' # modify based on your forcefields

    # setup of the adsorbent model
    frame = system.read_lammps('POLYMER_SYSTEM.lmps') # name of your fixed polymer system
    frame.forcefield = 'trappe/amber' # modify based on your forcefields

    
    # constant for loading calculations
    molec2mmols_g = 1e+3/frame.mass # molecules -> mmol/g
    print('CONVERSION (molecules to mmol/g) = ', molec2mmols_g)

    # defining the function that can run the GCMC simulation
    def calculate_isotherm_point_multi(gas_names, press):

        # css must be within the definition to remove interference from the loop structure
        css = cassandra.Cassandra(frame) # setup of the Cassandra object
        sim_settings = css.read_input('run_props.inp') # location of Cassandra parameters

        run_fldr = osp.join('results', str(press)) # all simulation results go into ./results

        # GCMC addition of gases to the fixed polymer frame
        css.add_gcmc(species=gases, # all parameters defined above; do not change this function
                     is_rigid=isrigid,
                     max_ins=maxins, 
                     chem_pot=chempots(press),
                     out_folder = run_fldr, **sim_settings)
        css.run()
        full_prp = css.run_queue[0].get_prp() # getting columns from the prp file created by Cassandra
        # modify this code if you are running more than 2 gases or changed parameters to be reported in .inp 
        step_return = full_prp[0]
        mols1_return = full_prp[6]
        mols2_return = full_prp[7]
        return [step_return, mols1_return, mols2_return]
    
    mols1 = []
    mols2 = []

    # looping through all pressures specified
    for p in gas_press:
        full_data = calculate_isotherm_point_multi(gas_names, p) # running the GCMC code
        step_data = full_data[0]
        mols1_data = full_data[1]
        mols2_data = full_data[2]
        mols1.append(mols1_data[:,numpy.newaxis])
        mols2.append(mols2_data[:,numpy.newaxis])

    # placing all species 1 nmol results in single file 
    mols1 = numpy.hstack(mols1)
    mols1 = pd.DataFrame(mols1)
    mols1.columns = numpy.tile(gas_press, 1)
    mols1 = mols1.add_prefix('P=')
    steps = pd.DataFrame(step_data)
    mols1.insert(0,'Steps', steps )
    mols1.to_string(gas_names[0]+'piece_1.txt', index=False)

    # placing all species 1 nmol results in single file 
    mols2 = numpy.hstack(mols2)
    mols2 = pd.DataFrame(mols2)
    mols2.columns = numpy.tile(gas_press, 1)
    mols2 = mols2.add_prefix('P=')
    steps = pd.DataFrame(step_data)
    mols2.insert(0,'Steps', steps )
    mols2.to_string(gas_names[1]+'piece_1.txt', index=False)

if __name__ == '__main__':
    run(True)