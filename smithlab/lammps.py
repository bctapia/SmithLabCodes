"""
Copyright 2025 Brandon C. Tapia

MIT License
"""

import numpy as np
from pathlib import Path



def reformat(lammps_in, lammps_out, lammps_ref):
    """
    Reformats the LAMMPS file to replace the style that might have been removed during polymatic
    """

    # TODO: add checking to make sure the reference file really is the correct reference

    with open(lammps_in, "r", encoding="utf-8") as file:
        lines = file.readlines()

    for i, line in enumerate(lines):

        stripped = line.strip()
        columns = stripped.split()

        if stripped.startswith("Masses"):
            preamble = lines[:i]

        if stripped.startswith("Atoms"):
            postamble = lines[i:]
            break

    with open(lammps_ref, "r", encoding="utf-8") as file:
        lines = file.readlines()

    for i, line in enumerate(lines):

        stripped = line.strip()
        columns = stripped.split()

        if stripped.startswith("Masses"):
            data_start = i

        if stripped.startswith("Atoms"):
            data_end = i
            break

    midamble = lines[data_start:data_end]

    with open(lammps_out, "w", encoding="utf-8") as file:
        file.writelines(preamble)
        file.writelines(midamble)
        file.writelines(postamble)

    return


def reorder_sections(lammps_in, lammps_out):
    """
    Sorts the Atoms and Velocities sections of a LAMMPS data file by atom ID.
    """
    atom_index = []
    atom_lines = []
    vel_index = []
    vel_lines = []

    pre_atom = []
    between_sections = []
    post_vel = []

    with open(lammps_in, "r", encoding="utf-8") as file:
        lines = file.readlines()

    in_atoms = False
    in_velocities = False
    atom_section_done = False

    for i, line in enumerate(lines):
        stripped = line.strip()
        columns = stripped.split()

        if stripped.startswith("Atoms"):
            in_atoms = True
            pre_atom = lines[: i + 1]
            continue

        if in_atoms and not atom_section_done:

            if not columns:
                continue

            if columns[0].isalpha():
                in_atoms = False
                atom_section_done = True
                between_start = i

                if stripped.startswith("Velocities"):
                    in_velocities = True
                    between_sections = lines[between_start : i + 1]
                continue

            atom_index.append(int(columns[0]))
            atom_lines.append(line)

        if in_velocities:
            if not columns:
                continue
            if columns[0].isalpha():
                post_vel = lines[i:]
                break
            vel_index.append(int(columns[0]))
            vel_lines.append(line)
            continue

    sorted_atom_idx = np.argsort(atom_index)
    sorted_vel_idx = np.argsort(vel_index)

    sorted_atom_lines = [atom_lines[i] for i in sorted_atom_idx]
    sorted_vel_lines = [vel_lines[i] for i in sorted_vel_idx]

    with open(lammps_out, "w", encoding="utf-8") as file:
        file.writelines(pre_atom)
        file.writelines("\n")
        file.writelines(sorted_atom_lines)
        file.writelines("\n")
        file.writelines(between_sections)
        file.writelines("\n")
        file.writelines(sorted_vel_lines)
        file.writelines("\n")
        file.writelines(post_vel)


def setup_lammps(lammps_in, lammps_out, lammps_ref):
    reformat(lammps_in, lammps_out, lammps_ref)
    reorder_sections(lammps_out, lammps_out)


def get_flucuating_prop(file_in, property, start=0, end=np.inf):

    in_props = False
    prop_array = np.array([])
    step_array = np.array([])

    with open(file_in, "r") as file:
        lines = file.readlines()

    for line in lines:

        if not line.strip():
            continue

        split_line = line.strip().split()

        if split_line[0] == "Step":

            in_props = True
            prop_idx = split_line.index(property)

            continue

        if split_line[0] == "Loop":
            in_props = False

        if in_props:
            step = int(split_line[0])

            if step >= start and step <= end:
                prop = float(split_line[prop_idx])
                step_array = np.append(step_array, step)
                prop_array = np.append(prop_array, prop)

    prop_avg = np.average(prop_array)
    prop_std = np.std(prop_array)

    return step_array, prop_array, prop_avg, prop_std


def polym_stats(lammps_in):

    in_atoms = False
    molecule_num = np.array([])

    with open(lammps_in, "r", encoding="utf-8") as file:
        lines = file.readlines()

    for i, line in enumerate(lines):

        stripped = line.strip()
        columns = stripped.split()

        if not columns:
            continue

        if stripped.startswith("Atoms"):
            in_atoms = True
            continue

        if in_atoms and columns[0].isalpha():
            break

        if in_atoms:
            molecule_num = np.append(molecule_num, int(columns[1]))

    molecule, count = np.unique(molecule_num, return_counts=True)

    return molecule, count


def get_mw(lammps_in):
    """
    Computes the molecular weight of a system
    """

    in_mass = False
    in_atoms = False

    atom_type = np.array([])
    atom_mass = np.array([])
    atoms = np.array([])

    with open(lammps_in, "r", encoding="utf-8") as file:
        lines = file.readlines()

    for i, line in enumerate(lines):
        stripped = line.strip()
        columns = stripped.split()

        if not columns:
            continue

        if stripped.startswith("Masses"):
            in_mass = True
            continue

        if in_mass and columns[0].isalpha():
            in_mass = False

        if in_mass:
            atom_type = np.append(atom_type, int(columns[0]))
            atom_mass = np.append(atom_mass, float(columns[1]))

        if stripped.startswith("Atoms"):
            in_atoms = True
            continue

        if in_atoms and columns[0].isalpha():
            in_atoms = False

        if in_atoms:
            atoms = np.append(atoms, columns[2])

    types, count = np.unique(atoms, return_counts=True)

    mw = 0
    for i, a_type in enumerate(types):
        for j, m_type in enumerate(atom_type):
            if int(a_type) == int(m_type):
                mw += atom_mass[j] * count[i]

    return float(mw)


def get_density(lammps_in):
    """
    Computes the density of a system from a LAMMPS data file.
    """
    na = 6.02214076 * 10**23
    count = 0

    with open(lammps_in, "r", encoding="utf-8") as file:
        lines = file.readlines()

    for line in lines:
        stripped = line.strip()
        columns = stripped.split()

        if stripped.endswith("xhi"):
            x_length = float(columns[1]) - float(columns[0])
            count += 1
        elif stripped.endswith("yhi"):
            y_length = float(columns[1]) - float(columns[0])
            count += 1
        elif stripped.endswith("zhi"):
            z_length = float(columns[1]) - float(columns[0])
            count += 1

        if count == 3:
            break

    mw = get_mw(lammps_in)

    density = (mw / na) / ((x_length * y_length * z_length) / 1e24)

    return density


def append_traj(traj_in, traj_out, time_spacer=100):
    """
    Appends LAMMPS trajectories that each start at timestep = 0.
    Adds a fixed offset between them to avoid overlapping timesteps.
    """
    time_offset = 0

    with open(traj_out, "w", encoding="utf-8") as outfile:
        for traj in traj_in:
            print(f"Processing {traj}")
            with open(traj, "r", encoding="utf-8") as infile:
                lines = infile.readlines()

            timesteps = [int(lines[i+1].strip()) for i, line in enumerate(lines)
                         if line.startswith("ITEM: TIMESTEP")]

            if not timesteps:
                print(f"Warning: no timesteps found in {traj}")
                continue

            last_t = max(timesteps)

            for i, line in enumerate(lines):
                if line.startswith("ITEM: TIMESTEP"):
                    old_t = int(lines[i+1].strip())
                    lines[i+1] = f"{old_t + time_offset}\n"

            outfile.writelines(lines)
            time_offset += last_t + time_spacer


def remove_atoms(lammps_in, lammps_out, types):

    atoms_to_keep = []

    in_atoms = False

    with open(lammps_in, "r", encoding="utf-8") as f:
        lines = f.readlines()

    for line in lines:
        stripped = line.strip()
        if not stripped:
            if in_atoms:
                in_atoms = False
            continue

        if stripped.startswith("Atoms"):
            in_atoms = True
            continue

        if in_atoms:
            if stripped.startswith("#") or not stripped[0].isdigit():
                continue

            cols = stripped.split()
            atom_id = int(cols[0])
            atom_type = int(cols[2])

            if atom_type not in types:
                atoms_to_keep.append(atom_id)

        elif stripped.startswith("Impropers"):
            in_atoms = False
            in_bonds = False
            in_angles = False
            in_dihedrals = False
            in_impropers = True

        if in_atoms:
            if float([columns[2]]) in types:
                # DO SOMETHING HERE
                return

    with open(lammps_out, "w", encoding="utf-8") as out:
        in_atoms = False

        for line in lines:
            stripped = line.strip()

            parts = stripped.split()
            if len(parts) == 2 and parts[1] == "atoms":
                out.write(f"{new_num_atoms} atoms\n")
                continue

            if stripped.startswith("Atoms"):
                in_atoms = True
                out.write(line)
                continue

            if in_atoms:
                if not stripped:
                    in_atoms = False
                    out.write(line)
                    continue

                if stripped.startswith("#") or not stripped[0].isdigit():
                    out.write(line)
                    continue

                cols = stripped.split()
                atom_type = int(cols[2])

                if atom_type in types:
                    continue
                else:
                    out.write(line)
                    continue

            out.write(line)


def change_partial_charge(lammps_in, charge_in, lammps_out, atom_style="full"):
    """
    Replace partial charges in a LAMMPS data file.

    Parameters
    ----------
    lammps_in : str or Path
        Input LAMMPS data file.
    charge_in : str or Path
        File containing charges; charges are read after skipping the first row.
    lammps_out : str or Path
        Output LAMMPS data file.
    atom_style : str
        Either 'full' or 'charge'.
    """

    lammps_in = Path(lammps_in)
    charge_in = Path(charge_in)
    lammps_out = Path(lammps_out)

    if atom_style == "full":
        charge_col = 3
        min_cols = 7
    elif atom_style == "charge":
        charge_col = 2
        min_cols = 6
    else:
        raise ValueError("atom_style must be 'full' or 'charge'")

    charges = np.genfromtxt(charge_in, skip_header=1)
    charges = np.asarray(charges).flatten()

    with open(lammps_in, "r") as f:
        lines = f.readlines()

    out_lines = []
    in_atoms_section = False
    atom_count = 0

    section_headers = {
        "Masses", "Velocities", "Bonds", "Angles", "Dihedrals", "Impropers",
        "Pair", "Bond", "Angle", "Dihedral", "Improper"
    }

    for i, line in enumerate(lines):
        stripped = line.strip()

        if stripped.startswith("Atoms"):
            in_atoms_section = True
            out_lines.append(line)
            continue

        if in_atoms_section and stripped == "":
            out_lines.append(line)
            continue

        if in_atoms_section:
            first_word = stripped.split()[0] if stripped else ""
            if first_word in section_headers:
                in_atoms_section = False
                out_lines.append(line)
                continue

            if not stripped or stripped.startswith("#"):
                out_lines.append(line)
                continue

            parts = line.split()

            if len(parts) < min_cols:
                raise ValueError(
                    f"Unexpected Atoms line format at line {i+1}: {line.rstrip()}"
                )

            if atom_count >= len(charges):
                raise ValueError("More atoms than charges were provided.")

            parts[charge_col] = f"{charges[atom_count]:.5f}"
            atom_count += 1
            out_lines.append(" ".join(parts) + "\n")
        else:
            out_lines.append(line)

    if atom_count != len(charges):
        raise ValueError(
            f"Updated {atom_count} atoms, but found {len(charges)} charges."
        )

    with open(lammps_out, "w") as f:
        f.writelines(out_lines)