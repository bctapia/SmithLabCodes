"""
Helper functions for the Automatic Charge Calculator III package
https://acc.biodata.ceitec.cz/

Copyright 2025 Brandon C. Tapia
"""
from pathlib import Path
import numpy as np


def change_charge(mol2_in, charge_in, mol2_out):
    with open(mol2_in, "r") as f:
        mol2_lines = f.readlines()

    charges = np.loadtxt(charge_in, dtype=float, skiprows=1)

    in_atom = False

    for i, line in enumerate(mol2_lines):
        if line.startswith("@<TRIPOS>ATOM"):
            in_atom = True
            continue
        if in_atom:
            # atom section ends at next tag or non-atom line
            if line.startswith("@<TRIPOS>") or not line.strip():
                in_atom = False
                continue

            cols = line.split() 
            if not cols or not cols[0].isdigit():
                in_atom = False
                continue

            atom_idx = int(cols[0]) - 1
            q = charges[atom_idx]


            mol2_lines[i] = (
                f"{int(cols[0]):>7d} "
                f"{cols[1]:<8s}"
                f"{float(cols[2]):>10.4f}"
                f"{float(cols[3]):>10.4f}"
                f"{float(cols[4]):>10.4f}"
                f" {cols[5]:<6s}"
                f"{int(cols[6]):>5d} "
                f"{cols[7]:<8s}"
                f"{q:>10.4f}\n"
            )

    # Write once
    with open(mol2_out, "w") as f:
        f.writelines(mol2_lines)
