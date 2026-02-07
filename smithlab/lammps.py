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

def combine_coeffs(*args, coeffs_out="ff_system.in"):
    """
    Combine multiple LAMMPS coeff include files into one, renumbering type IDs
    so each file occupies a unique block of IDs.

    Supported commands (first token):
      mass
      pair_coeff
      bond_coeff
      angle_coeff
      dihedral_coeff
      improper_coeff

    Notes:
    - For pair_coeff, this handles the common "pair_coeff i j eps sigma"
    """

    in_files = [Path(a) for a in args]

    # Running offsets (how much to add to types from each new file)
    off_mass = off_pair = off_bond = off_angle = off_dihedral = off_improper = 0

    # Track max index seen in *current* file to update offsets after processing it
    out_lines = []

    out_lines.append(f"# read in as: include {coeffs_out}\n")

    for fname in in_files:
        text = fname.read_text(encoding="utf-8").splitlines(True)  # keep line endings

        max_mass = max_pair = max_bond = max_angle = max_dihedral = max_improper = 0

        out_lines.append(f"\n# ================= {fname.name} =================\n")
        out_lines.append(f"# molecule ID filename toff {off_pair} boff {off_bond} aoff {off_angle} doff {off_dihedral} ioff {off_improper}\n")
        for line in text:
            raw = line
            s = line.strip()

            # keep blank lines and comments untouched
            if not s or s.startswith("#"):
                out_lines.append(raw)
                continue

            cols = s.split()
            cmd = cols[0]

            if cmd == "mass":
                t = int(cols[1]) + off_mass
                cols[1] = str(t)
                max_mass = max(max_mass, t)

            elif cmd == "pair_coeff":
                # Common single-type form: pair_coeff <type> eps sigma ...
                # If you use i j form, see note below.
                t = int(cols[1]) + off_pair
                cols[1] = str(t)
                cols[2] = str(t)
                max_pair = max(max_pair, t)

            elif cmd == "bond_coeff":
                t = int(cols[1]) + off_bond
                cols[1] = str(t)
                max_bond = max(max_bond, t)

            elif cmd == "angle_coeff":
                t = int(cols[1]) + off_angle
                cols[1] = str(t)
                max_angle = max(max_angle, t)

            elif cmd == "dihedral_coeff":
                t = int(cols[1]) + off_dihedral
                cols[1] = str(t)
                max_dihedral = max(max_dihedral, t)

            elif cmd == "improper_coeff":
                t = int(cols[1]) + off_improper
                cols[1] = str(t)
                max_improper = max(max_improper, t)

            # If it's some other command, keep it as-is (e.g. pair_modify, kspace_style, etc.)
            else:
                print("Warning: Writing a line that may be a global variable and probably should be defined elsewhere:")
                print(raw)
                out_lines.append(raw)
                continue

            # Reconstruct line with original newline
            newline = "\n" if raw.endswith("\n") else ""
            out_lines.append(" ".join(cols) + newline)
        out_lines.append("\n")

        # After finishing THIS file, update offsets for the next file.
        # Offsets should jump to (current max + 1) so the next file starts after it.
        # If a file had none of that command, max_* stays 0 and we don't move that offset.
        if max_mass:
            off_mass = max_mass
        if max_pair:
            off_pair = max_pair
        if max_bond:
            off_bond = max_bond
        if max_angle:
            off_angle = max_angle
        if max_dihedral:
            off_dihedral = max_dihedral
        if max_improper:
            off_improper = max_improper

        if max_mass != max_pair:
            raise ValueError("Mass and pair types should be equivalent")

    Path(coeffs_out).write_text("".join(out_lines), encoding="utf-8")
