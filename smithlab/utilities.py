"""
Helper functioons for the Automatic Charge Calculator III package
https://acc.biodata.ceitec.cz/

Copyright 2025 Brandon C. Tapia
"""
from pathlib import Path
import numpy as np

def change_charge(mol2_in, charge_in, mol2_out):
    """
    Changes the partial charges of a mol2 file to those specified in the charge_in text file.
    """
    # Read the mol2 file
    with open(mol2_in, "r") as f:
        mol2_lines = f.readlines()

    # Read the charge file
    charge_lines = np.loadtxt(charge_in, dtype=str, skiprows=1)