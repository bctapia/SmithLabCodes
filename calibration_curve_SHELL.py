'''
This script has been written to find the chemical potential vs pressure calibration
required for usage of the MC/MD or MC codes 
The python package pysimm (https://pysimm.org/) is required to run this code
Pysimm has been modified for parallel optimization and is available here: https://github.com/bctapia/pysimm
Extensive parts of this code have been taken from pysimm examples
Written by: Brandon C. Tapia (bctapia@mit.edu) while in the Smith Research Lab at the Massachusetts Institute of Technology
The MIT License is in effect
Please cite this code if used
'''

from pysimm import cassandra
from pysimm import system
from os import path as osp
import numpy # this can be removed if 'chem_pot_inp' (line 29) is modified
import pandas as pd


def run(test=False):
    '''
    Code to create a pressure-chemical potential calibration curve required for GCMC simulations
    The code runs GCMC simulations for a reservoir of gas molecules in contact with a fixed system of gas molecules 
    '''

    # Gas name as it will be referred through simulations
    gas_names = ['GAS'] # only one gas allowed

    # Chemical potentials to sample at
    chem_pot_inp = numpy.arange(-26,-15,1) # make sure to input negative numbers; # could also use numpy.linspace() or numpy.arange() if you import numpy

    chem_pots = [lambda x: x] # don't edit this line

    # Subdirectory for insertable gas data
    data_dir = osp.join('.', 'gases')

    # Setup of gas to be inserted
    gases = []
    for gn in gas_names:
        gases.append(system.read_lammps(osp.join(data_dir, gn + '.lmps')))
        gases[-1].forcefield = 'trappe/amber' # change depending on gas forcefields you are using

    # Setup of fixed gas system
    frame = system.read_lammps('GAS_SYSTEM.lmps')
    frame.forcefield = 'trappe/amber' # change depending on gas forcefields you are using

    # Setup of the GCMC simulations
    css = cassandra.Cassandra(frame)
    sim_settings = css.read_input('run_props.inp')

    # GCMC simulations to calculate pressure
    def calculate_isotherm_point(gas_name, cp):

        if cp >= 0:
            print('[WARNING]: chemical potentials should be negative; results may be incorrect')
            cp_flip=cp
            post = 'WARNING_POSITIVE_chem_pot_'
        else:
            cp_flip = -cp # '-' can't be in file names
            post = 'neg_chem_pot'
        
        run_fldr = osp.join(gas_name, str(cp_flip)+post)
        idx = gas_names.index(gas_name)
        # sim_settings.update({'Run_Name':  'gcmc'})
        css.add_gcmc(species=gases[idx], is_new=True, chem_pot=chem_pots[idx](cp),
                     out_folder=run_fldr, props_file='gcmc.inp', max_ins=20000,
                     **sim_settings)
        css.run()
        full_prp = css.run_queue[0].get_prp()
        return [full_prp[0],full_prp[4]]

    pressure = []

    for gn in gas_names:
        for c in chem_pot_inp:
            data = calculate_isotherm_point(gn, c)
            pressure_data = data[1]
            step_data = data[0]
            pressure.append(pressure_data[:, numpy.newaxis])

    pressure = numpy.hstack(pressure)
    pressure = pd.DataFrame(pressure)
    steps = pd.DataFrame(step_data)
    pressure.columns = chem_pot_inp
    pressure = pressure.add_prefix('P (bar) for CP=')
    pressure.insert(0,'Steps', steps )
    pressure.to_string('pressure.txt', index=False)

    chem_pot = pd.DataFrame(chem_pot_inp)
    chem_pot.columns = ["Chemical_Potential"]
    chem_pot.to_string('Chemical_Potential.txt', index=False)
    
    

if __name__ == '__main__':
    run(False)

