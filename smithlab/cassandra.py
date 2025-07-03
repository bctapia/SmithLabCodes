import os
import re
import numpy as np
from smithlab import lammps


def get_param(prp_in, prop, return_step=True):
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
    if return_step:
        return steps, prop_array
    else:
        return prop_array


def avg_prop(data_in, include=1/3, std=True):
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
    print(prop_std)
    if std:
        return prop_avg, prop_std
    else:
        return prop_avg


def loading(prp_in, lammps_in=None, mw=None):
    """
    Computes the gravimetric loading from a Cassandra result.out.prp file
    If lammps_in is specified, the molar mass of the system is automatically determined
    Otherwise, mw must specify the molar mass (g/mol)
    """
    if lammps_in:
        mw_use = lammps.get_mw(lammps_in)
    else:
        mw_use = mw

    molec_inserted = get_param(prp_in, "Nmols_2")
    molec_avg = avg_prop(molec_inserted)





# ===================================WRAPPERS=================================== #
# the following are wrappers to "automate" some tasks using the above functions


def pressure_calib(files_in):
    """
    Computes the chemical potential vs pressure calibration curve from the provided Cassandra files
    Requires that
    """
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
                    break

            steps, press = get_param(result_file, "Pressure")

    return


def MC_isotherm():
    """
    Computes the sorption isotherm from the provided Cassandra files
    """
    return


def GCMC_isotherm():

    return
