import os
import re
import numpy as np


def pressure_calibration(files_in):
    """
    Computes the chemical potential vs pressure calibration curve from the provided Cassandra files
    """
    subdirs = [root for root, dirs, files in os.walk(files_in)]
    # print them if you like

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

            with open(result_file, "r", encoding="utf-8") as file:
                lines = file.readlines()

            for i, line in enumerate(lines):
                stripped = line.strip()

                if stripped.startswith("# MC_STEP"):
                    columns = stripped.split()

                    for j, text in enumerate(columns):
                        if text.startswith("MC_STEP"):
                            step_idx = j - 1
                        if text.startswith("Pressure"):
                            p_idx = j - 1

                    break

            file = np.transpose(np.genfromtxt(result_file))
            steps = file[step_idx]
            press = file[p_idx]
            print(steps)

            # print(chempot)

    return


def MC_isotherm():
    """
    Computes the sorption isotherm from the provided Cassandra files
    """
    return


def GCMC_isotherm():

    return
