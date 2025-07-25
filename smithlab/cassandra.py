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


def avg_prop(data_in, include=float(1 / 3)):
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

    molec_2_mmolg = 1e3 / mw_use

    # the math behind the next line:
    # moles inserted = molec_inserted/N_A
    # mass of one mol of system: mw/N_A
    # loading [mol/g] = (molec_inserted/N_A)/(mw/N_A) = moles_inserted/mw
    # loading [mmol/g] = (1E3/mw)*moles_inserted = molec_2_mmolg*moles_inserted
    loading_data = molec_2_mmolg * molec_avg

    loading_std = molec_2_mmolg * molec_std

    return loading_data, loading_std


# TODO: allow for mw to be passed
def concentration(loading, loading_std, lammps_in=None, pol_density=None):
    """
    Computes the volumetric sorption capacity (concentration) (cm^3_STP/cm^3_pol) from loading
    """
    if pol_density is None:
        density = lammps.get_density(lammps_in)
    else:
        density = pol_density

    c = 22.414 * loading * density
    c_std = 22.414 * loading_std * density
    return c, c_std


def get_max_step(file):

    steps, _ = get_param(file, "MC_STEP")  # prop we get doesn't matter, we just need the steps

    return steps[-1]


# ===================================WRAPPERS===================================
# the following are wrappers to "automate" some tasks using the above functions


def pressure_calib(files_in, include=[0.1, 40], plot=True, steps=2e6):
    """
    Computes the chemical potential vs pressure calibration curve from the provided Cassandra files
    Requires that
    """
    chempot_array = np.array([])
    step_array = np.array([])
    press_array = np.array([])
    press_std_array = np.array([])

    max_step_array = np.array([])

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

            steps_file, press = get_param(result_file, "Pressure")

            max_step_array = np.append(max_step_array, steps_file[-1])
            press_avg, press_std = avg_prop(press)
            press_array = np.append(press_array, press_avg)
            press_std_array = np.append(press_std_array, press_std)

    chempot_array_idx = np.argsort(chempot_array)
    chempot_array = chempot_array[chempot_array_idx]
    max_step_array = max_step_array[chempot_array_idx]
    press_array = press_array[chempot_array_idx]
    press_std_array = press_std_array[chempot_array_idx]

    minval = np.argmin(abs(press_array - include[0]))
    maxval = np.argmin(abs(press_array - include[1]))

    if include[0] < press_array[0]:
        print(f"Warning: {include[0]} is below the minimum pressure in the data ({press_array[0]})")
    elif press_array[minval] > include[0]:
        minval = minval - 1

    if include[1] > press_array[-1]:
        print(
            f"Warning: {include[1]} is above the maximum pressure in the data ({press_array[-1]})"
        )
        maxval = maxval + 1
    elif press_array[maxval] < include[1]:
        maxval = maxval + 1  # TODO: +1 or +2??

    chempot_array_idx = chempot_array_idx[minval:maxval]
    chempot_array = chempot_array[minval:maxval]
    max_step_array = max_step_array[minval:maxval]
    press_array = press_array[minval:maxval]
    press_std_array = press_std_array[minval:maxval]

    for i, step_data in enumerate(max_step_array):
        if step_data < steps:
            print(
                f"Warning: {int(step_data)} steps for mu={chempot_array[i]:.2f} (p={press_array[i]:.2f}) is < {int(steps)} steps"
            )

    # x0 = [4, -80]
    # popt = np.polyfit(np.log(press_array), chempot_array, deg=1)
    # reg_loss = regression.loss(x0, press_array, press_std_array)

    # result = minimize(regression.loss, x0, args=(press_array, chempot_array))
    def log_model(x, a, b):
        return a * np.log(x) + b

    # result = result.x
    # print(press_array)
    result, pcov = curve_fit(log_model, press_array, chempot_array)
    if plot:
        press_array_model = np.linspace(press_array[0], press_array[-1], 500)
        fig, ax = plt.subplots(1, 2)
        ax[0].plot(press_array, chempot_array, "o")
        ax[0].plot(press_array_model, result[0] * np.log(press_array_model) + result[1])
        ax[0].set_xlabel(r"$P\;\mathrm{(bar)}$")
        ax[0].set_ylabel(r"$\mu\;\mathrm{(kJ\;mol^{-1}})$")
        ax[1].plot(np.log(press_array), chempot_array, "o")
        ax[1].plot(np.log(press_array_model), result[0] * np.log(press_array_model) + result[1])
        ax[1].set_xlabel(r"$\mathrm{ln}(P)\;\mathrm{(ln(bar))}$")
        ax[1].set_ylabel(r"$\mu\;\mathrm{(kJ\;mol^{-1}})$")
        plt.tight_layout()
        plt.show()
    # return press_array, press_std_array, chempot_array, result.x
    return result[0], result[1]


# TODO: configure MW once configured in concentration()
def mc_isotherm(files_in, lammps_in, pol_density=None, mw=None, steps=2e6):
    """
    Computes the sorption isotherm from the provided Cassandra files
    """

    pressure_array = np.array([])
    sorption_array = np.array([])
    sorption_std_array = np.array([])

    max_step_array = np.array([])

    subdirs = [root for root, dirs, file in os.walk(files_in)]
    for subdir in subdirs:

        result_file = os.path.join(subdir, "result.out.prp")

        if os.path.isfile(result_file):
            split_path = os.path.normpath(subdir).split(os.path.sep)
            loading_data, loading_std = loading(result_file, lammps_in)

            if pol_density is not None:
                loading_data, loading_std = concentration(
                    loading_data, loading_std, lammps_in=None, pol_density=pol_density
                )
            else:
                loading_data, loading_std = concentration(
                    loading_data, loading_std, lammps_in=lammps_in
                )

            pressure_array = np.append(pressure_array, float(split_path[-1]))
            sorption_array = np.append(sorption_array, loading_data)
            sorption_std_array = np.append(sorption_std_array, loading_std)
            max_step_array = np.append(max_step_array, get_max_step(result_file))

    pressure_array_idx = np.argsort(pressure_array)
    pressure_array = pressure_array[pressure_array_idx]
    sorption_array = sorption_array[pressure_array_idx]
    sorption_std_array = sorption_std_array[pressure_array_idx]
    max_step_array = max_step_array[pressure_array_idx]

    # TODO: is it better to have steps returned (see mcmd_isotherm) to avoid forced printing?
    for i, step_data in enumerate(max_step_array):
        if step_data < steps:
            print(
                f"Warning: {int(step_data)} steps for P={pressure_array[i]:.2f} is < {int(steps)} steps"
            )

    return pressure_array, sorption_array, sorption_std_array


def mcmd_isotherm(files_in, lammps_start, pol_density=None, mw=None):
    """
    files_in will be the Trial folder
    """
    iter_array = np.array([])
    pressure_array = np.array([])
    sorption_array = np.array([])
    sorption_std_array = np.array([])

    max_step_array = np.array([])

    subdirs = [root for root, dirs, file in os.walk(files_in)]

    for subdir in subdirs:

        # print(subdir)
        head, tail = os.path.split(subdir)

        # check if the folder is a number
        pressure = None
        try:
            pressure = float(tail)
        except ValueError:
            continue

        # if it is, check if this folder has a results folder
        res_folder = os.path.join(subdir, "results")
        if os.path.exists(res_folder):

            # find the .lmps file in the current directory:
            lammps_use_path = None
            num_lammps = 0

            for file in os.listdir(subdir):
                if file.endswith(".lmps") and file.startswith(lammps_start):
                    lammps_use_path = os.path.join(subdir, file)
                    num_lammps += 1
            if num_lammps > 1:
                print("multiple .lmps found, unsure what to do, exiting")
                break

            # find the X.gcmc.prp with the highest number
            iter_num = 0
            for file in os.listdir(res_folder):
                if file.endswith(".gcmc.prp"):
                    iter_current = int(file.split(".")[0])
                    iter_num = iter_current if iter_current > iter_num else iter_num
            # TODO: this is going to fail if it cant find the folder. will need to indent probably
            result_file = os.path.join(res_folder, f"{iter_num}.gcmc.prp")

            loading_data, loading_std = loading(result_file, lammps_use_path)

            if pol_density is not None:
                loading_data, loading_std = concentration(
                    loading_data, loading_std, lammps_in=None, pol_density=pol_density
                )
            else:
                loading_data, loading_std = concentration(
                    loading_data, loading_std, lammps_in=lammps_use_path
                )

            pressure_array = np.append(pressure_array, pressure)
            sorption_array = np.append(sorption_array, loading_data)
            sorption_std_array = np.append(sorption_std_array, loading_std)
            max_step_array = np.append(max_step_array, get_max_step(result_file))
            iter_array = np.append(iter_array, iter_num)

    pressure_array_idx = np.argsort(pressure_array)
    pressure_array = pressure_array[pressure_array_idx]
    sorption_array = sorption_array[pressure_array_idx]
    sorption_std_array = sorption_std_array[pressure_array_idx]
    max_step_array = max_step_array[pressure_array_idx]
    iter_array = iter_array[pressure_array_idx]

    print(pressure_array, sorption_array, sorption_std_array, max_step_array, iter_array)


def mcmd_iters(files_in):

    pressure_array = np.array([])
    iter_array = np.array([])
    molec_array = np.array([])
    molec_std_array = np.array([])
    max_step_array = np.array([])

    subdirs = [root for root, dirs, file in os.walk(files_in)]

    for subdir in subdirs:

        head, tail = os.path.split(subdir)

        pressure = None
        try:
            pressure = float(tail)
        except ValueError:
            continue
            
        iter_ind = np.array([])
        molec_ind = np.array([])
        molec_std_ind = np.array([])
        max_step_ind = np.array([])

        res_folder = os.path.join(subdir, "results")
        if os.path.exists(res_folder):
            for file in os.listdir(res_folder):

                if file.endswith(".gcmc.prp"):
                    iter_current = int(file.split(".")[0])
                
                    result_file = os.path.join(res_folder, f"{iter_current}.gcmc.prp")
                
                    _, molec = get_param(result_file, "Nmols_2")
                    molec_avg, molec_std = avg_prop(molec)

                    iter_ind = np.append(iter_ind, iter_current)
                    molec_ind = np.append(molec_ind, molec_avg)
                    molec_std_ind = np.append(molec_std_ind, molec_std)
                    max_step_ind = np.append(max_step_ind, get_max_step(result_file))

            iter_ind_idx = np.argsort(iter_ind)
            iter_ind = iter_ind[iter_ind_idx]
            molec_ind = molec_ind[iter_ind_idx]
            molec_std_ind = molec_std_ind[iter_ind_idx]
            max_step_ind = max_step_ind[iter_ind_idx]

            pressure_array = np.append(pressure_array, pressure)
            iter_array = np.append(iter_array, iter_ind)
            molec_array = np.append(molec_array, molec_ind)
            molec_std_array = np.append(molec_std_array, molec_std_ind)
            max_step_array = np.append(max_step_array, max_step_ind)

    # TODO: will need to perform an argsort to get pressures sorted but for now this is good
    return pressure_array, iter_array, molec_array, molec_std_array, max_step_array

# TODO: consider adding Gelmin-Rubin test for data to include
