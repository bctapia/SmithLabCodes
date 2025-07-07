import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, curve_fit
from smithlab import lammps, regression


def get_param(prp_in, prop):
    """
    Extracts a parameter from a .prp Cassandra file
    """

    with open(prp_in, "r", encoding="utf-8") as file:
        lines = file.readlines()

    for i, line in enumerate(lines):
        if i>5:
            break
        stripped = line.strip()
        if stripped.startswith("# MC_STEP"):
            columns = stripped.split()

            for j, text in enumerate(columns):
                if text == "MC_STEP":
                    step_idx = j - 1
                if text == prop:
                    p_idx = j - 1
                    break

    data = np.genfromtxt(prp_in, comments="#", unpack=True)

    steps = data[step_idx]
    prop_array = data[p_idx]

    return steps.astype(int), prop_array.astype(float)


def avg_prop(data_in, include=1/3):
    """
    calculates the ensemble average of a specified parameter by calculating the parameter over the "include" fraction
    Example: average_prop(inserted_array, 1/3) calculates the mean and standard deviation from the final 1/3 of inserted_array
    """
    if include > float(1):
        print("include is a fraction")

    first_idx_use = int(np.ceil(len(data_in) * (1 - include)))

    data_use = data_in[first_idx_use:]

    prop_avg = np.average(data_use)
    prop_std = np.std(data_use)

    return prop_avg, prop_std


def loading(prp_in, lammps_in=None, mw=None):
    """
    Computes the gravimetric sorption capacity (loading) [mmol/g] from a Cassandra result.out.prp file
    If lammps_in is specified, the molar mass of the system is automatically determined
    Otherwise, mw must specify the molar mass (g/mol)
    """
    
    if lammps_in:
        mw_use = lammps.get_mw(lammps_in)
    else:
        mw_use = mw

    _, molec_inserted = get_param(prp_in, "Nmols_2")
    
    molec_avg, molec_std = avg_prop(molec_inserted)
    
    molec_2_mmolg = 1E3/mw_use
    
    # the math behind the next line: 
    # moles inserted = molec_inserted/N_A
    # mass of one mol of system: mw/N_A
    # loading [mol/g] = (molec_inserted/N_A)/(mw/N_A) = moles_inserted/mw
    # loading [mmol/g] = (1E3/mw)*moles_inserted = molec_2_mmolg*moles_inserted
    loading_data = molec_2_mmolg * molec_avg
    
    loading_std = molec_2_mmolg * molec_std

    return loading_data, loading_std


def concentration(loading, loading_std, pol_density):
    """
    Computes the volumetric sorption capacity (concentration) (cm^3_STP/cm^3_pol) from loading
    """
    c = 22.414 * loading * pol_density
    c_std = 22.414 * loading_std * pol_density
    return c, c_std


# ===================================WRAPPERS=================================== #
# the following are wrappers to "automate" some tasks using the above functions


def pressure_calib(files_in, include=[0.1, 40], plot=True):
    """
    Computes the chemical potential vs pressure calibration curve from the provided Cassandra files
    Requires that
    """
    chempot_array = np.array([])
    press_array = np.array([])
    press_std_array = np.array([])

    subdirs = [root for root, dirs, files in os.walk(files_in)]

    for subdir in subdirs:
        result_file = os.path.join(subdir, "result.out.prp")
        input_file = os.path.join(subdir, "gcmc.inp")

        if os.path.isfile(result_file) and os.path.isfile(input_file):

            with open(input_file, "r", encoding="utf-8") as file:
                lines = file.readlines()

            for i, line in enumerate(lines):
                stripped = line.strip()
                if stripped.startswith("# Chemical_Potential_Info"):
                    chempot = float(lines[i + 1])
                    chempot_array = np.append(chempot_array, chempot)
                    break
            
            steps, press = get_param(result_file, "Pressure")
            press_avg, press_std = avg_prop(press)
            press_array = np.append(press_array, press_avg)
            press_std_array = np.append(press_std_array, press_std)
    
    chempot_array_idx = np.argsort(chempot_array)
    chempot_array = chempot_array[chempot_array_idx]
    press_array = press_array[chempot_array_idx]
    press_std_array = press_std_array[chempot_array_idx]
    
    minval = np.argmin(abs(press_array-include[0]))
    maxval = np.argmin(abs(press_array-include[1]))
    
    chempot_array_idx = chempot_array_idx[minval:maxval]
    chempot_array = chempot_array[minval:maxval]
    press_array = press_array[minval:maxval]
    press_std_array = press_std_array[minval:maxval]

    x0 = [4,-80]
    #popt = np.polyfit(np.log(press_array), chempot_array, deg=1)
    reg_loss = regression.loss(x0, press_array, press_std_array)
    
    result = minimize(regression.loss, x0, args=(press_array, chempot_array))
    #def log_model(x, a, b):
    #    return a*np.log(x)+b
    result = result.x
    #popt, pcov = curve_fit(log_model, press_std_array, chempot_array)
    if plot:
        fig, ax = plt.subplots(1,2)
        ax[0].plot(press_array, chempot_array, 'o')
        ax[0].plot(press_array, result[0]*np.log(press_array)+result[1], 'o')
        ax[0].set_xlabel(r"$P\;\mathrm{(bar)}$")
        ax[0].set_ylabel(r"$\mu\;\mathrm{(kJ\;mol^{-1}})$")
        ax[1].plot(np.log(press_array), chempot_array, 'o')
        ax[1].plot(np.log(press_array), result[0]*np.log(press_array)+result[1], 'o')
        plt.show()
    return press_array, press_std_array, chempot_array, result.x
    #return press_array, press_std_array, chempot_array, popt


def MC_isotherm(files_in, lammps_in, pol_density=None, mw=None):
    """
    Computes the sorption isotherm from the provided Cassandra files
    """

    pressure_array = []
    sorption_array = []
    sorption_std_array = []
    subdirs = [root for root, dirs, file in os.walk(files_in)]

    for subdir in subdirs:
        result_file = os.path.join(subdir, "result.prp.out")

        if os.path.isfile(result_file):

            loading_data, loading_std = loading(result_file, lammps_in)

            if pol_density:
                loading_data, loading_std = concentration(loading_data, loading_std, pol_density)
        # look for filepath to get the pressure
        # pressure_array.append(float(pressure_val))
        sorption_array.append(loading_data)
        sorption_std_array.append(loading_std)

    return #pressure, sorption_array, sorption_std_array

def GCMC_isotherm():

    return

# TODO: consider adding Gelmin-Rubin test for data to include