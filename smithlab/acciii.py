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


def change_charge_lammps(lmp_in, charge_in, lmp_out):
    with open(lmp_in, "r") as f:
        lmp_lines = f.readlines()

    # assumes one charge per atom, same ordering
    charges = np.loadtxt(charge_in, dtype=float, skiprows=1)

    in_atoms = False
    atom_idx = 0

    for i, line in enumerate(lmp_lines):
        stripped = line.strip()

        # Enter Atoms section
        if stripped.startswith("Atoms"):
            in_atoms = True
            continue

        # Skip blank line right after header
        if in_atoms and stripped == "":
            continue

        # Exit Atoms section when next section starts
        if in_atoms and (
            stripped.startswith("Bonds") or
            stripped.startswith("Angles") or
            stripped.startswith("Dihedrals") or
            stripped.startswith("Impropers") or
            stripped.endswith("Coeffs")
        ):
            in_atoms = False
            continue

        # Modify atom lines
        if in_atoms:
            parts = line.split()

            # skip malformed lines
            if len(parts) < 7:
                continue

            # atom_style full:
            # id mol type q x y z ...
            parts[3] = f"{charges[atom_idx]:.6f}"
            atom_idx += 1

            # preserve formatting loosely
            lmp_lines[i] = " ".join(parts) + "\n"

    with open(lmp_out, "w") as f:
        f.writelines(lmp_lines)