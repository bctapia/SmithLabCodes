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


def neutralize_lammps_charge(lmp_in, lmp_out, atom_types=None, charge_col=3, type_col=2):
    """
    Neutralize total charge in a LAMMPS data file by distributing the charge
    correction across selected atom types.

    Assumes atom_style full:
        id mol type q x y z ...

    Parameters
    ----------
    lmp_in : str or Path
        Input LAMMPS data file.
    lmp_out : str or Path
        Output LAMMPS data file.
    atom_types : None, int, or iterable of int
        Atom types over which to distribute the correction.
        If None, all atoms are used.
    charge_col : int
        Zero-based charge column. For atom_style full, q is column 3.
    type_col : int
        Zero-based atom type column. For atom_style full, type is column 2.
    """

    with open(lmp_in, "r") as f:
        lmp_lines = f.readlines()

    if atom_types is not None:
        if isinstance(atom_types, int):
            atom_types = {atom_types}
        else:
            atom_types = set(atom_types)

    atom_line_indices = []
    charges = []
    types = []

    in_atoms = False

    for i, line in enumerate(lmp_lines):
        stripped = line.strip()

        if stripped.startswith("Atoms"):
            in_atoms = True
            continue

        if in_atoms and stripped == "":
            continue

        if in_atoms and (
            stripped.startswith("Bonds") or
            stripped.startswith("Angles") or
            stripped.startswith("Dihedrals") or
            stripped.startswith("Impropers") or
            stripped.endswith("Coeffs")
            or stripped.startswith("Velocities")
        ):
            in_atoms = False
            continue

        if in_atoms:
            parts = line.split()

            if len(parts) <= max(charge_col, type_col):
                continue

            atom_line_indices.append(i)
            types.append(int(parts[type_col]))
            charges.append(float(parts[charge_col]))

    charges = np.array(charges, dtype=float)
    types = np.array(types, dtype=int)

    if len(charges) == 0:
        raise ValueError("No atom lines found in Atoms section.")

    if atom_types is None:
        selected = np.ones(len(charges), dtype=bool)
    else:
        selected = np.isin(types, list(atom_types))

    if not np.any(selected):
        raise ValueError(f"No atoms found with atom_types={atom_types}")

    total_charge = charges.sum()
    correction_per_atom = -total_charge / selected.sum()

    charges[selected] += correction_per_atom

    for local_idx, line_idx in enumerate(atom_line_indices):
        parts = lmp_lines[line_idx].split()
        parts[charge_col] = f"{charges[local_idx]:.6f}"
        lmp_lines[line_idx] = " ".join(parts) + "\n"

    with open(lmp_out, "w") as f:
        f.writelines(lmp_lines)

    return {
        "initial_total_charge": total_charge,
        "final_total_charge": charges.sum(),
        "n_corrected_atoms": int(selected.sum()),
        "correction_per_atom": correction_per_atom,
    }