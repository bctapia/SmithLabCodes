'''
This script has been written to automate grand canonical Monte Carlo/molecular dynamics iterations
The python package pysimm (https://pysimm.org/) is required to run this code
Pysimm has been modified for parallel optimization and is available here: https://github.com/bctapia/pysimm
Extensive parts of this code have been taken from pysimm examples
Written by: Brandon C. Tapia (bctapia@mit.edu) while in the Smith Research Lab at the Massachusetts Institute of Technology
The MIT License is in effect
Please cite this code if used
'''
# Imports
from pysimm.apps import mc_md 
from pysimm import system
import numpy

# The chemical potential is related to the system pressure with this equation
# See the 'chemical potential to pressure' base code to develop your own calibrations
chem_pot = lambda x: 2.4997*numpy.log(x)-41.9215 # dependent on temperature and gas inserted

press = 1 # bar # change depending on what pressure your system is at
temp = 308 # K # change depending on what temperature your system is at

def run(test=False): 
    frame = system.read_lammps('POLYMER_SYSTEM.lmps') # change name depending on your polymer LAMMPS file
    frame.forcefield = 'trappe/amber' # change depending on polymer forcefields you are using

    gas1 = system.read_lammps('c3h6.lmps') # change name depending on your gas LAMMPS file
    gas1.forcefield = 'trappe/amber' # change depending on gas forcefields you are using
    
    mc_props = {'rigid_type': True, # only FALSE if modeled as single point (e.g UA CH4)
                'max_ins': 20000, # See Cassandra (https://cassandra.nd.edu/) for understanding of these variables
                'Chemical_Potential_Info': chem_pot(press),
                'Temperature_Info': temp,
                'Rcutoff_Low': 0.4,
                'Run_Type': {'steps': 100},
                'CBMC_Info': {'rcut_cbmc': 2.0},
                'Simulation_Length_Info': {'run': 2000000,
                                           'coord_freq': 25000,
                                           'prop_freq': 100},
                'VDW_Style': {'cut_val': 15.0},
                'Charge_Style': {'cut_val': 15.0},
                'Property_Info': {'prop1': 'energy_total',
                                  'prop2': 'pressure',
                                  'prop3': 'nmols',
                                  'prop4': 'density'}}
    
    md_props = {'temp': temp, # See LAMMPS (https://www.lammps.org) for understanding of these variables
                'pressure': {'start': press,
                             'iso': 'iso'},
                'kspace_style': 1e-5,
                'timestep': 0.5,
                'cutoff': 15.0,
                'length': 2000000,
                'thermo': 25000,
                'dump': 25000,
                'print_to_screen': False} # True is fine but slows computation down
    
    # np and prefix are not standard pysimm variables. Added to allow for parallelization
    sim_result = mc_md.mc_md(gas1, frame, mcmd_niter=10, sim_folder='results', mc_props=mc_props, md_props=md_props, np=48, prefix='mpiexec')
    sim_result.write_lammps('SYSplusGAS.lmps') # writing out final systems
    sim_result.write_xyz('SYSplusGAS.xyz') # writing out final systems

if __name__ == '__main__':
    run(False)

