"""
Copyright 2025 Brandon C. Tapia

MIT License
"""

from collections import defaultdict
from pathlib import Path
import re
import numpy as np


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


def add_xyz_atoms_to_lammps_data(
    lammps_in,
    xyz_in,
    lammps_out,
    timestep,
    label_to_type,
    keep_existing_atoms=False,
    atom_style="full",
    charge=0.0,
    expand_box=True,
    padding=0.0,
):
    """
    Add atoms from a selected timestep of an XYZ-like trajectory to a LAMMPS data file.

    Each newly added atom gets its own molecule ID when atom_style is "full" or
    "molecular".

    Parameters
    ----------
    lammps_in : str
        Input LAMMPS data file.

    xyz_in : str
        XYZ-like trajectory file with frames beginning with lines like:

            MC_STEP: 2000000

        followed by atom coordinate lines:

            A0  x y z
            C   x y z

    lammps_out : str
        Output LAMMPS data file.

    timestep : int
        Timestep to extract from xyz_in.

    label_to_type : dict
        Mapping from XYZ labels to LAMMPS atom types.

        Example:
            {"A0": 1}

        This adds all A0 atoms as atom type 1.

    keep_existing_atoms : bool, default True
        If True, keep atoms already in the LAMMPS data file and append new atoms.

        If False, remove existing atoms, velocities, bonds, angles, dihedrals,
        and impropers, then add only selected XYZ atoms.

    atom_style : str, default "full"
        Atom style used for writing the Atoms section.

        Supported:
            "full":
                atom-ID molecule-ID atom-type charge x y z

            "molecular":
                atom-ID molecule-ID atom-type x y z

            "atomic":
                atom-ID atom-type x y z

            "charge":
                atom-ID atom-type charge x y z

    charge : float, default 0.0
        Charge assigned to newly added atoms for atom styles that include charge.

    expand_box : bool, default True
        If True, expand the simulation box to contain all output atoms.

    padding : float, default 0.0
        Extra padding added to the expanded box.

    Notes
    -----
    This function only adds atoms. It does not add bonds, angles, dihedrals,
    or impropers for the newly added atoms.
    """

    section_names = [
        "Masses",
        "Pair Coeffs",
        "Bond Coeffs",
        "Angle Coeffs",
        "Dihedral Coeffs",
        "Improper Coeffs",
        "Atoms",
        "Velocities",
        "Bonds",
        "Angles",
        "Dihedrals",
        "Impropers",
    ]

    count_patterns = {
        "atoms": r"^\s*(\d+)\s+atoms\b",
        "bonds": r"^\s*(\d+)\s+bonds\b",
        "angles": r"^\s*(\d+)\s+angles\b",
        "dihedrals": r"^\s*(\d+)\s+dihedrals\b",
        "impropers": r"^\s*(\d+)\s+impropers\b",
        "atom_types": r"^\s*(\d+)\s+atom\s+types\b",
        "bond_types": r"^\s*(\d+)\s+bond\s+types\b",
        "angle_types": r"^\s*(\d+)\s+angle\s+types\b",
        "dihedral_types": r"^\s*(\d+)\s+dihedral\s+types\b",
        "improper_types": r"^\s*(\d+)\s+improper\s+types\b",
    }

    def strip_comment(line):
        return line.split("#", 1)[0].strip()

    def is_section_header(line):
        clean = strip_comment(line)
        return clean in section_names

    def parse_lammps_data(filename):
        with open(filename, "r") as f:
            raw_lines = [line.rstrip("\n") for line in f]

        header_lines = []
        sections = {name: [] for name in section_names}

        current_section = None
        found_first_section = False

        for line in raw_lines:
            if is_section_header(line):
                current_section = strip_comment(line)
                found_first_section = True
                continue

            if not found_first_section:
                header_lines.append(line)
            else:
                if current_section is not None:
                    if line.strip() != "":
                        sections[current_section].append(line)

        counts = {key: 0 for key in count_patterns}

        for line in header_lines:
            for key, pattern in count_patterns.items():
                m = re.search(pattern, line)
                if m:
                    counts[key] = int(m.group(1))

        box = {
            "xlo": None,
            "xhi": None,
            "ylo": None,
            "yhi": None,
            "zlo": None,
            "zhi": None,
        }

        for line in header_lines:
            parts = line.split()
            if len(parts) >= 4:
                if parts[2] == "xlo" and parts[3] == "xhi":
                    box["xlo"] = float(parts[0])
                    box["xhi"] = float(parts[1])
                elif parts[2] == "ylo" and parts[3] == "yhi":
                    box["ylo"] = float(parts[0])
                    box["yhi"] = float(parts[1])
                elif parts[2] == "zlo" and parts[3] == "zhi":
                    box["zlo"] = float(parts[0])
                    box["zhi"] = float(parts[1])

        return {
            "header_lines": header_lines,
            "sections": sections,
            "counts": counts,
            "box": box,
        }

    def parse_xyz_timestep(filename, target_timestep):
        selected_atoms = []
        in_target_frame = False
        found_target = False

        with open(filename, "r") as f:
            for raw_line in f:
                line = raw_line.strip()

                if line == "":
                    continue

                if line.startswith("MC_STEP"):
                    parts = line.replace(":", " ").split()

                    step_value = None
                    for p in parts:
                        try:
                            step_value = int(p)
                            break
                        except ValueError:
                            pass

                    if step_value is None:
                        raise ValueError(f"Could not parse timestep from line: {line}")

                    if step_value == target_timestep:
                        in_target_frame = True
                        found_target = True
                        selected_atoms = []
                    else:
                        if in_target_frame:
                            break
                        in_target_frame = False

                    continue

                if in_target_frame:
                    parts = line.split()
                    if len(parts) < 4:
                        continue

                    label = parts[0]

                    try:
                        x = float(parts[1])
                        y = float(parts[2])
                        z = float(parts[3])
                    except ValueError:
                        continue

                    selected_atoms.append(
                        {
                            "label": label,
                            "x": x,
                            "y": y,
                            "z": z,
                        }
                    )

        if not found_target:
            raise ValueError(f"Timestep {target_timestep} was not found in {filename}.")

        return selected_atoms

    def count_valid_lines(lines):
        n = 0
        for line in lines:
            clean = strip_comment(line)
            if clean:
                n += 1
        return n

    def max_atom_id(atom_lines):
        max_id = 0

        for line in atom_lines:
            clean = strip_comment(line)
            if not clean:
                continue

            parts = clean.split()
            if len(parts) < 1:
                continue

            try:
                max_id = max(max_id, int(parts[0]))
            except ValueError:
                pass

        return max_id

    def max_molecule_id(atom_lines, atom_style):
        """
        Return max molecule ID from Atoms lines.

        For atom_style='full':
            atom-ID mol-ID atom-type charge x y z

        For atom_style='molecular':
            atom-ID mol-ID atom-type x y z

        For atom_style='atomic' or 'charge', molecule IDs are not present.
        """
        if atom_style not in ["full", "molecular"]:
            return 0

        max_mol = 0

        for line in atom_lines:
            clean = strip_comment(line)
            if not clean:
                continue

            parts = clean.split()
            if len(parts) < 2:
                continue

            try:
                max_mol = max(max_mol, int(parts[1]))
            except ValueError:
                pass

        return max_mol

    def make_atom_line(atom_id, atom_type, x, y, z, mol_id=None):
        if atom_style == "full":
            if mol_id is None:
                raise ValueError("mol_id is required for atom_style='full'.")

            return (
                f"{atom_id:d} {mol_id:d} {atom_type:d} "
                f"{charge:.8f} {x:.8f} {y:.8f} {z:.8f}"
            )

        elif atom_style == "molecular":
            if mol_id is None:
                raise ValueError("mol_id is required for atom_style='molecular'.")

            return (
                f"{atom_id:d} {mol_id:d} {atom_type:d} "
                f"{x:.8f} {y:.8f} {z:.8f}"
            )

        elif atom_style == "atomic":
            return f"{atom_id:d} {atom_type:d} {x:.8f} {y:.8f} {z:.8f}"

        elif atom_style == "charge":
            return (
                f"{atom_id:d} {atom_type:d} {charge:.8f} "
                f"{x:.8f} {y:.8f} {z:.8f}"
            )

        else:
            raise ValueError(
                "Unsupported atom_style. Use one of: "
                "'full', 'molecular', 'atomic', 'charge'."
            )

    def parse_atom_xyz_from_lammps_line(line, atom_style):
        clean = strip_comment(line)
        if not clean:
            return None

        parts = clean.split()

        try:
            if atom_style == "full":
                if len(parts) < 7:
                    return None
                return float(parts[4]), float(parts[5]), float(parts[6])

            elif atom_style == "molecular":
                if len(parts) < 6:
                    return None
                return float(parts[3]), float(parts[4]), float(parts[5])

            elif atom_style == "atomic":
                if len(parts) < 5:
                    return None
                return float(parts[2]), float(parts[3]), float(parts[4])

            elif atom_style == "charge":
                if len(parts) < 6:
                    return None
                return float(parts[3]), float(parts[4]), float(parts[5])

            else:
                return None

        except ValueError:
            return None

    def write_section(f, name, lines):
        clean_lines = []

        for line in lines:
            if strip_comment(line):
                clean_lines.append(line)

        if len(clean_lines) == 0:
            return

        f.write(f"\n{name}\n\n")
        for line in clean_lines:
            f.write(f"{line}\n")

    if atom_style not in ["full", "molecular", "atomic", "charge"]:
        raise ValueError(
            "Unsupported atom_style. Use one of: "
            "'full', 'molecular', 'atomic', 'charge'."
        )

    data = parse_lammps_data(lammps_in)
    sections = data["sections"]
    counts = data["counts"]
    box = data["box"]

    xyz_atoms = parse_xyz_timestep(xyz_in, timestep)

    atoms_to_add = []
    for atom in xyz_atoms:
        label = atom["label"]

        if label in label_to_type:
            atom_type = int(label_to_type[label])

            atoms_to_add.append(
                {
                    "label": label,
                    "type": atom_type,
                    "x": atom["x"],
                    "y": atom["y"],
                    "z": atom["z"],
                }
            )

    if keep_existing_atoms:
        output_atoms = list(sections["Atoms"])
        next_atom_id = max_atom_id(output_atoms) + 1
    else:
        output_atoms = []
        next_atom_id = 1

        sections["Velocities"] = []
        sections["Bonds"] = []
        sections["Angles"] = []
        sections["Dihedrals"] = []
        sections["Impropers"] = []

        counts["bonds"] = 0
        counts["angles"] = 0
        counts["dihedrals"] = 0
        counts["impropers"] = 0

    if atom_style in ["full", "molecular"]:
        next_molecule_id = max_molecule_id(output_atoms, atom_style) + 1
    else:
        next_molecule_id = None

    for atom in atoms_to_add:
        if atom_style in ["full", "molecular"]:
            mol_id = next_molecule_id
            next_molecule_id += 1
        else:
            mol_id = None

        atom_line = make_atom_line(
            next_atom_id,
            atom["type"],
            atom["x"],
            atom["y"],
            atom["z"],
            mol_id=mol_id,
        )

        output_atoms.append(atom_line)
        next_atom_id += 1

    sections["Atoms"] = output_atoms
    counts["atoms"] = count_valid_lines(sections["Atoms"])

    if atoms_to_add:
        max_added_type = max(atom["type"] for atom in atoms_to_add)
        counts["atom_types"] = max(counts["atom_types"], max_added_type)

    if not keep_existing_atoms:
        counts["bonds"] = 0
        counts["angles"] = 0
        counts["dihedrals"] = 0
        counts["impropers"] = 0

    if expand_box and atoms_to_add:
        xs = []
        ys = []
        zs = []

        for line in sections["Atoms"]:
            xyz = parse_atom_xyz_from_lammps_line(line, atom_style)
            if xyz is not None:
                x, y, z = xyz
                xs.append(x)
                ys.append(y)
                zs.append(z)

        if xs:
            if box["xlo"] is None:
                box["xlo"] = min(xs) - padding
            else:
                box["xlo"] = min(box["xlo"], min(xs) - padding)

            if box["xhi"] is None:
                box["xhi"] = max(xs) + padding
            else:
                box["xhi"] = max(box["xhi"], max(xs) + padding)

            if box["ylo"] is None:
                box["ylo"] = min(ys) - padding
            else:
                box["ylo"] = min(box["ylo"], min(ys) - padding)

            if box["yhi"] is None:
                box["yhi"] = max(ys) + padding
            else:
                box["yhi"] = max(box["yhi"], max(ys) + padding)

            if box["zlo"] is None:
                box["zlo"] = min(zs) - padding
            else:
                box["zlo"] = min(box["zlo"], min(zs) - padding)

            if box["zhi"] is None:
                box["zhi"] = max(zs) + padding
            else:
                box["zhi"] = max(box["zhi"], max(zs) + padding)

    if box["xlo"] is None:
        box["xlo"], box["xhi"] = -10.0, 10.0
    if box["ylo"] is None:
        box["ylo"], box["yhi"] = -10.0, 10.0
    if box["zlo"] is None:
        box["zlo"], box["zhi"] = -10.0, 10.0

    with open(lammps_out, "w") as f:
        f.write("LAMMPS data file modified by add_xyz_atoms_to_lammps_data\n\n")

        f.write(f"{counts['atoms']} atoms\n")
        f.write(f"{counts['bonds']} bonds\n")
        f.write(f"{counts['angles']} angles\n")
        f.write(f"{counts['dihedrals']} dihedrals\n")
        f.write(f"{counts['impropers']} impropers\n\n")

        f.write(f"{counts['atom_types']} atom types\n")
        f.write(f"{counts['bond_types']} bond types\n")
        f.write(f"{counts['angle_types']} angle types\n")

        if counts["dihedral_types"] > 0:
            f.write(f"{counts['dihedral_types']} dihedral types\n")

        if counts["improper_types"] > 0:
            f.write(f"{counts['improper_types']} improper types\n")

        f.write("\n")
        f.write(f"{box['xlo']:.8f} {box['xhi']:.8f} xlo xhi\n")
        f.write(f"{box['ylo']:.8f} {box['yhi']:.8f} ylo yhi\n")
        f.write(f"{box['zlo']:.8f} {box['zhi']:.8f} zlo zhi\n")

        write_section(f, "Masses", sections["Masses"])
        write_section(f, "Pair Coeffs", sections["Pair Coeffs"])
        write_section(f, "Bond Coeffs", sections["Bond Coeffs"])
        write_section(f, "Angle Coeffs", sections["Angle Coeffs"])
        write_section(f, "Dihedral Coeffs", sections["Dihedral Coeffs"])
        write_section(f, "Improper Coeffs", sections["Improper Coeffs"])

        write_section(f, "Atoms", sections["Atoms"])
        write_section(f, "Velocities", sections["Velocities"])

        write_section(f, "Bonds", sections["Bonds"])
        write_section(f, "Angles", sections["Angles"])
        write_section(f, "Dihedrals", sections["Dihedrals"])
        write_section(f, "Impropers", sections["Impropers"])

def combine_data_files(
    lammps_array_in,
    lammps_out,
    include_atoms=None,
    include_velocities=None,
):
    """
    Combine multiple LAMMPS data files into one.

    Parameters
    ----------
    lammps_array_in : list[str]
        List of input LAMMPS data file paths.

    lammps_out : str
        Output LAMMPS data file path.

    include_atoms : list[bool] or None
        Whether to include atoms and corresponding bonds, angles, dihedrals,
        impropers from each input file.

        If None, atoms/topology are included for all files.

    include_velocities : list[bool] or None
        Whether to include velocities from each input file.

        If None, velocities are included for all files where atoms are included.

    Notes
    -----
    This function assumes conventional LAMMPS data-file sections such as:

        Masses
        Pair Coeffs
        Bond Coeffs
        Angle Coeffs
        Dihedral Coeffs
        Improper Coeffs
        Atoms
        Velocities
        Bonds
        Angles
        Dihedrals
        Impropers

    It preserves coefficient text after the first type-ID column, while offsetting
    type IDs as needed.

    Atom lines are assumed to have:

        atom-ID molecule-ID atom-type ...

    which is common for atom styles such as full, molecular, charge, etc.

    Bond lines are assumed to have:

        bond-ID bond-type atom1 atom2

    Angle lines:

        angle-ID angle-type atom1 atom2 atom3

    Dihedral lines:

        dihedral-ID dihedral-type atom1 atom2 atom3 atom4

    Improper lines:

        improper-ID improper-type atom1 atom2 atom3 atom4
    """

    lammps_array_in = list(lammps_array_in)
    n_files = len(lammps_array_in)

    if include_atoms is None:
        include_atoms = [True] * n_files

    if include_velocities is None:
        include_velocities = [True] * n_files

    if len(include_atoms) != n_files:
        raise ValueError("include_atoms must have the same length as lammps_array_in.")

    if len(include_velocities) != n_files:
        raise ValueError("include_velocities must have the same length as lammps_array_in.")

    section_names = [
        "Masses",
        "Pair Coeffs",
        "Bond Coeffs",
        "Angle Coeffs",
        "Dihedral Coeffs",
        "Improper Coeffs",
        "Atoms",
        "Velocities",
        "Bonds",
        "Angles",
        "Dihedrals",
        "Impropers",
    ]

    coeff_sections = {
        "Masses": "atom_type",
        "Pair Coeffs": "atom_type",
        "Bond Coeffs": "bond_type",
        "Angle Coeffs": "angle_type",
        "Dihedral Coeffs": "dihedral_type",
        "Improper Coeffs": "improper_type",
    }

    topology_sections = {
        "Bonds": ("bond", "bond_type", 2),
        "Angles": ("angle", "angle_type", 3),
        "Dihedrals": ("dihedral", "dihedral_type", 4),
        "Impropers": ("improper", "improper_type", 4),
    }

    def strip_comment(line):
        return line.split("#", 1)[0].strip()

    def split_comment(line):
        if "#" in line:
            body, comment = line.split("#", 1)
            return body.rstrip(), " # " + comment.strip()
        return line.rstrip(), ""

    def is_section_header(line):
        clean = strip_comment(line)
        return clean in section_names

    def parse_data_file(filename):
        filename = Path(filename)

        with open(filename, "r") as f:
            raw_lines = f.readlines()

        header_lines = []
        sections = {name: [] for name in section_names}

        current_section = None
        found_first_section = False

        for line in raw_lines:
            stripped = line.strip()

            if is_section_header(line):
                current_section = strip_comment(line)
                found_first_section = True
                continue

            if not found_first_section:
                header_lines.append(line.rstrip("\n"))
                continue

            if current_section is None:
                continue

            if stripped == "":
                continue

            if stripped.startswith("#"):
                continue

            sections[current_section].append(line.rstrip("\n"))

        counts = parse_header_counts(header_lines)
        box = parse_box(header_lines)

        return {
            "filename": str(filename),
            "header_lines": header_lines,
            "sections": sections,
            "counts": counts,
            "box": box,
        }

    def parse_header_counts(header_lines):
        patterns = {
            "atoms": r"^\s*(\d+)\s+atoms\b",
            "bonds": r"^\s*(\d+)\s+bonds\b",
            "angles": r"^\s*(\d+)\s+angles\b",
            "dihedrals": r"^\s*(\d+)\s+dihedrals\b",
            "impropers": r"^\s*(\d+)\s+impropers\b",
            "atom_types": r"^\s*(\d+)\s+atom\s+types\b",
            "bond_types": r"^\s*(\d+)\s+bond\s+types\b",
            "angle_types": r"^\s*(\d+)\s+angle\s+types\b",
            "dihedral_types": r"^\s*(\d+)\s+dihedral\s+types\b",
            "improper_types": r"^\s*(\d+)\s+improper\s+types\b",
        }

        counts = {key: 0 for key in patterns}

        for line in header_lines:
            for key, pat in patterns.items():
                m = re.search(pat, line)
                if m:
                    counts[key] = int(m.group(1))

        return counts

    def parse_box(header_lines):
        box = {
            "xlo": None,
            "xhi": None,
            "ylo": None,
            "yhi": None,
            "zlo": None,
            "zhi": None,
        }

        for line in header_lines:
            parts = line.split()
            if len(parts) >= 4:
                if parts[2] == "xlo" and parts[3] == "xhi":
                    box["xlo"] = float(parts[0])
                    box["xhi"] = float(parts[1])
                elif parts[2] == "ylo" and parts[3] == "yhi":
                    box["ylo"] = float(parts[0])
                    box["yhi"] = float(parts[1])
                elif parts[2] == "zlo" and parts[3] == "zhi":
                    box["zlo"] = float(parts[0])
                    box["zhi"] = float(parts[1])

        return box

    def count_nonempty(lines):
        return sum(1 for line in lines if strip_comment(line))

    def max_first_int(lines):
        max_val = 0
        for line in lines:
            clean = strip_comment(line)
            if not clean:
                continue
            parts = clean.split()
            if len(parts) > 0:
                max_val = max(max_val, int(parts[0]))
        return max_val

    def offset_coeff_line(line, type_offset):
        body, comment = split_comment(line)
        parts = body.split()
        if not parts:
            return None

        parts[0] = str(int(parts[0]) + type_offset)
        return " ".join(parts) + comment

    def offset_atom_line(line, atom_offset, mol_offset, atom_type_offset):
        body, comment = split_comment(line)
        parts = body.split()
        if len(parts) < 3:
            raise ValueError(f"Atom line has fewer than 3 columns: {line}")

        parts[0] = str(int(parts[0]) + atom_offset)
        parts[1] = str(int(parts[1]) + mol_offset)
        parts[2] = str(int(parts[2]) + atom_type_offset)

        return " ".join(parts) + comment

    def offset_velocity_line(line, atom_offset):
        body, comment = split_comment(line)
        parts = body.split()
        if len(parts) < 4:
            raise ValueError(f"Velocity line has fewer than 4 columns: {line}")

        parts[0] = str(int(parts[0]) + atom_offset)
        return " ".join(parts) + comment

    def offset_topology_line(line, id_offset, type_offset, atom_offset, n_atoms_in_topology):
        body, comment = split_comment(line)
        parts = body.split()

        expected_min_cols = 2 + n_atoms_in_topology
        if len(parts) < expected_min_cols:
            raise ValueError(
                f"Topology line has fewer than {expected_min_cols} columns: {line}"
            )

        parts[0] = str(int(parts[0]) + id_offset)
        parts[1] = str(int(parts[1]) + type_offset)

        for i in range(n_atoms_in_topology):
            atom_col = 2 + i
            parts[atom_col] = str(int(parts[atom_col]) + atom_offset)

        return " ".join(parts) + comment

    parsed_files = [parse_data_file(fname) for fname in lammps_array_in]

    combined = {name: [] for name in section_names}

    atom_offset = 0
    mol_offset = 0

    bond_offset = 0
    angle_offset = 0
    dihedral_offset = 0
    improper_offset = 0

    atom_type_offset = 0
    bond_type_offset = 0
    angle_type_offset = 0
    dihedral_type_offset = 0
    improper_type_offset = 0

    total_atoms = 0
    total_bonds = 0
    total_angles = 0
    total_dihedrals = 0
    total_impropers = 0

    total_atom_types = 0
    total_bond_types = 0
    total_angle_types = 0
    total_dihedral_types = 0
    total_improper_types = 0

    for file_index, data in enumerate(parsed_files):
        sections = data["sections"]

        use_atoms = include_atoms[file_index]
        use_velocities = include_velocities[file_index] and use_atoms

        current_atom_offset = atom_offset
        current_mol_offset = mol_offset

        current_bond_offset = bond_offset
        current_angle_offset = angle_offset
        current_dihedral_offset = dihedral_offset
        current_improper_offset = improper_offset

        current_atom_type_offset = atom_type_offset
        current_bond_type_offset = bond_type_offset
        current_angle_type_offset = angle_type_offset
        current_dihedral_type_offset = dihedral_type_offset
        current_improper_type_offset = improper_type_offset

        type_offsets = {
            "atom_type": current_atom_type_offset,
            "bond_type": current_bond_type_offset,
            "angle_type": current_angle_type_offset,
            "dihedral_type": current_dihedral_type_offset,
            "improper_type": current_improper_type_offset,
        }

        # Always append type/coefficient sections.
        for section, type_key in coeff_sections.items():
            offset = type_offsets[type_key]
            for line in sections[section]:
                clean = strip_comment(line)
                if clean:
                    combined[section].append(offset_coeff_line(line, offset))

        if use_atoms:
            for line in sections["Atoms"]:
                clean = strip_comment(line)
                if clean:
                    combined["Atoms"].append(
                        offset_atom_line(
                            line,
                            current_atom_offset,
                            current_mol_offset,
                            current_atom_type_offset,
                        )
                    )

            if use_velocities:
                for line in sections["Velocities"]:
                    clean = strip_comment(line)
                    if clean:
                        combined["Velocities"].append(
                            offset_velocity_line(line, current_atom_offset)
                        )

            # Append topology sections only when atoms are included.
            for section, info in topology_sections.items():
                topology_name, type_key, n_atoms_in_topology = info

                if topology_name == "bond":
                    id_offset = current_bond_offset
                elif topology_name == "angle":
                    id_offset = current_angle_offset
                elif topology_name == "dihedral":
                    id_offset = current_dihedral_offset
                elif topology_name == "improper":
                    id_offset = current_improper_offset
                else:
                    raise RuntimeError("Unknown topology type.")

                type_offset = type_offsets[type_key]

                for line in sections[section]:
                    clean = strip_comment(line)
                    if clean:
                        combined[section].append(
                            offset_topology_line(
                                line,
                                id_offset,
                                type_offset,
                                current_atom_offset,
                                n_atoms_in_topology,
                            )
                        )

        # Update type offsets.
        n_atom_types = max_first_int(sections["Masses"])
        n_bond_types = max_first_int(sections["Bond Coeffs"])
        n_angle_types = max_first_int(sections["Angle Coeffs"])
        n_dihedral_types = max_first_int(sections["Dihedral Coeffs"])
        n_improper_types = max_first_int(sections["Improper Coeffs"])

        atom_type_offset += n_atom_types
        bond_type_offset += n_bond_types
        angle_type_offset += n_angle_types
        dihedral_type_offset += n_dihedral_types
        improper_type_offset += n_improper_types

        total_atom_types += n_atom_types
        total_bond_types += n_bond_types
        total_angle_types += n_angle_types
        total_dihedral_types += n_dihedral_types
        total_improper_types += n_improper_types

        if use_atoms:
            n_atoms = count_nonempty(sections["Atoms"])
            n_bonds = count_nonempty(sections["Bonds"])
            n_angles = count_nonempty(sections["Angles"])
            n_dihedrals = count_nonempty(sections["Dihedrals"])
            n_impropers = count_nonempty(sections["Impropers"])

            atom_offset += n_atoms
            bond_offset += n_bonds
            angle_offset += n_angles
            dihedral_offset += n_dihedrals
            improper_offset += n_impropers

            total_atoms += n_atoms
            total_bonds += n_bonds
            total_angles += n_angles
            total_dihedrals += n_dihedrals
            total_impropers += n_impropers

            # Molecule offset is taken from max molecule ID in this file.
            max_mol_id = 0
            for line in sections["Atoms"]:
                clean = strip_comment(line)
                if not clean:
                    continue
                parts = clean.split()
                if len(parts) >= 2:
                    max_mol_id = max(max_mol_id, int(parts[1]))

            mol_offset += max_mol_id

    # Determine output box as the min/max envelope of all input boxes.
    xlo_values = []
    xhi_values = []
    ylo_values = []
    yhi_values = []
    zlo_values = []
    zhi_values = []

    for data in parsed_files:
        box = data["box"]
        if box["xlo"] is not None:
            xlo_values.append(box["xlo"])
            xhi_values.append(box["xhi"])
        if box["ylo"] is not None:
            ylo_values.append(box["ylo"])
            yhi_values.append(box["yhi"])
        if box["zlo"] is not None:
            zlo_values.append(box["zlo"])
            zhi_values.append(box["zhi"])

    if xlo_values:
        xlo, xhi = min(xlo_values), max(xhi_values)
    else:
        xlo, xhi = 0.0, 1.0

    if ylo_values:
        ylo, yhi = min(ylo_values), max(yhi_values)
    else:
        ylo, yhi = 0.0, 1.0

    if zlo_values:
        zlo, zhi = min(zlo_values), max(zhi_values)
    else:
        zlo, zhi = 0.0, 1.0

    def write_section(f, name, lines):
        if len(lines) == 0:
            return

        f.write(f"\n{name}\n\n")
        for line in lines:
            f.write(f"{line}\n")

    with open(lammps_out, "w") as f:
        f.write("LAMMPS data file combined by combine_data_files\n\n")

        f.write(f"{total_atoms} atoms\n")
        f.write(f"{total_bonds} bonds\n")
        f.write(f"{total_angles} angles\n")
        f.write(f"{total_dihedrals} dihedrals\n")
        f.write(f"{total_impropers} impropers\n\n")

        f.write(f"{total_atom_types} atom types\n")
        f.write(f"{total_bond_types} bond types\n")
        f.write(f"{total_angle_types} angle types\n")
        f.write(f"{total_dihedral_types} dihedral types\n")
        f.write(f"{total_improper_types} improper types\n\n")

        f.write(f"{xlo:.8f} {xhi:.8f} xlo xhi\n")
        f.write(f"{ylo:.8f} {yhi:.8f} ylo yhi\n")
        f.write(f"{zlo:.8f} {zhi:.8f} zlo zhi\n")

        write_section(f, "Masses", combined["Masses"])
        write_section(f, "Pair Coeffs", combined["Pair Coeffs"])
        write_section(f, "Bond Coeffs", combined["Bond Coeffs"])
        write_section(f, "Angle Coeffs", combined["Angle Coeffs"])
        write_section(f, "Dihedral Coeffs", combined["Dihedral Coeffs"])
        write_section(f, "Improper Coeffs", combined["Improper Coeffs"])

        write_section(f, "Atoms", combined["Atoms"])

        if len(combined["Velocities"]) > 0:
            write_section(f, "Velocities", combined["Velocities"])

        write_section(f, "Bonds", combined["Bonds"])
        write_section(f, "Angles", combined["Angles"])
        write_section(f, "Dihedrals", combined["Dihedrals"])
        write_section(f, "Impropers", combined["Impropers"])


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


def parse_local_dump(
    filename,
    start_timestep=0,
    selected_types=None,
    outfile_prefix=None
):
    data = defaultdict(list)

    with open(filename) as f:
        lines = iter(f)

        skip_frame = False
        n_entries = 0

        for line in lines:

            if line.startswith("ITEM: TIMESTEP"):

                timestep = int(next(lines).strip())
                skip_frame = timestep < start_timestep

            elif line.startswith("ITEM: NUMBER OF ENTRIES"):

                n_entries = int(next(lines).strip())

            elif line.startswith("ITEM: ENTRIES"):

                for _ in range(n_entries):

                    cols = next(lines).split()

                    type_id = int(cols[0])
                    value = float(cols[1])

                    if skip_frame:
                        continue

                    if selected_types is not None:
                        if type_id not in selected_types:
                            continue

                    data[type_id].append(value)

    if outfile_prefix:

        for type_id, values in data.items():

            outfile = f"{outfile_prefix}_type_{type_id}.dat"

            with open(outfile, "w") as f:

                f.write("value\n")

                for v in values:
                    f.write(f"{v:.8f}\n")

            print(f"Wrote {outfile}")

    return data


def _parse_angles_from_data(data_file, angle_type=None):
    """
    Parse a LAMMPS data file and return angles (optionally filtered by type)
    along with a mapping from atom_id -> atom type name (from Masses comments).
    """
    type_id_to_name = {}
    atom_type_ids = {}
    atom_type_names = {}
    angles = []

    section = None
    section_names = {
        "Masses", "Atoms", "Bonds", "Angles", "Dihedrals", "Impropers",
        "Velocities", "Pair Coeffs", "Bond Coeffs", "Angle Coeffs",
        "Dihedral Coeffs", "Improper Coeffs",
    }

    with open(data_file, "r", encoding="utf-8") as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue

            no_comment = stripped.split("#", 1)[0].strip()
            comment = stripped.split("#", 1)[1].strip() if "#" in stripped else ""

            if no_comment in section_names:
                section = no_comment
                continue

            if section == "Masses":
                cols = no_comment.split()
                if len(cols) >= 2:
                    type_id = int(cols[0])
                    type_id_to_name[type_id] = (
                        comment.split()[0] if comment else str(type_id)
                    )

            elif section == "Atoms":
                cols = no_comment.split()
                if len(cols) >= 7:
                    atom_id = int(cols[0])
                    atom_type_ids[atom_id] = int(cols[2])

            elif section == "Angles":
                cols = no_comment.split()
                if len(cols) >= 5:
                    ang_id = int(cols[0])
                    ang_type = int(cols[1])
                    a1, a2, a3 = map(int, cols[2:5])
                    if angle_type is None or ang_type == angle_type:
                        angles.append((ang_id, a1, a2, a3))

    for atom_id, type_id in atom_type_ids.items():
        atom_type_names[atom_id] = type_id_to_name.get(type_id, str(type_id))

    return angles, atom_type_names


def _bond_angle(p1, p2, p3):
    """Bond angle (degrees) at vertex p2 formed by p1-p2-p3."""
    v1 = p1 - p2
    v2 = p3 - p2

    n1 = np.linalg.norm(v1)
    n2 = np.linalg.norm(v2)

    if n1 == 0 or n2 == 0:
        return np.nan

    cos_theta = np.dot(v1, v2) / (n1 * n2)
    # Clamp for numerical stability
    cos_theta = max(-1.0, min(1.0, cos_theta))
    return np.degrees(np.arccos(cos_theta))


def calc_angles(
    data_file,
    dump_file,
    angle_type,
    output_file=None,
    max_frames=None,
    start_timestep=0,
):
    """
    Compute bond angles for all angles of a given type across a LAMMPS dump.

    Parameters
    ----------
    data_file : str or Path
        LAMMPS data file containing the topology (Masses, Atoms, Angles).
    dump_file : str or Path
        LAMMPS trajectory file (must contain id and either xu/yu/zu or x/y/z).
    angle_type : int
        Angle type ID to extract.
    output_file : str or Path or None, optional
        If provided, results are streamed to this CSV file as they are computed.
        Columns: timestep, angle_id, atom1, atom2, atom3, angle_deg.
    max_frames : int or None, optional
        Maximum number of frames to process. None = all frames.
    start_timestep : int, optional
        Skip frames whose timestep is below this value.

    Returns
    -------
    results : dict
        {
            "angles_topology": list of (ang_id, a1, a2, a3),
            "atom_type_names": {atom_id: type_name},
            "timesteps": np.ndarray of shape (n_frames,),
            "angles": np.ndarray of shape (n_frames, n_angles) in degrees,
        }
    """
    data_file = Path(data_file)
    dump_file = Path(dump_file)

    angles_top, atom_type_names = _parse_angles_from_data(
        data_file, angle_type=angle_type
    )

    if not angles_top:
        raise ValueError(
            f"No angles of type {angle_type} found in {data_file}."
        )

    timesteps = []
    angles_per_frame = []

    out = None
    if output_file is not None:
        out = open(output_file, "w", encoding="utf-8")
        out.write("timestep,angle_id,atom1,atom2,atom3,angle_deg\n")

    try:
        with open(dump_file, "r", encoding="utf-8") as f:
            frame_count = 0

            while True:
                line = f.readline()
                if not line:
                    break
                if not line.startswith("ITEM: TIMESTEP"):
                    continue

                timestep = int(f.readline().strip())

                f.readline()                       # ITEM: NUMBER OF ATOMS
                n_atoms = int(f.readline().strip())

                f.readline(); f.readline(); f.readline(); f.readline()  # box bounds

                atom_header = f.readline().strip().split()
                columns = atom_header[2:]
                id_col = columns.index("id")

                if all(c in columns for c in ["xu", "yu", "zu"]):
                    x_col = columns.index("xu")
                    y_col = columns.index("yu")
                    z_col = columns.index("zu")
                elif all(c in columns for c in ["x", "y", "z"]):
                    x_col = columns.index("x")
                    y_col = columns.index("y")
                    z_col = columns.index("z")
                else:
                    raise ValueError(
                        "Dump must contain either xu yu zu or x y z."
                    )

                skip_frame = timestep < start_timestep
                coords = {}
                for _ in range(n_atoms):
                    cols = f.readline().split()
                    if skip_frame:
                        continue
                    atom_id = int(cols[id_col])
                    coords[atom_id] = np.array(
                        [float(cols[x_col]),
                         float(cols[y_col]),
                         float(cols[z_col])]
                    )

                if skip_frame:
                    continue

                frame_angles = np.empty(len(angles_top))
                for i, (ang_id, a1, a2, a3) in enumerate(angles_top):
                    theta = _bond_angle(
                        coords[a1], coords[a2], coords[a3]
                    )
                    frame_angles[i] = theta

                    if out is not None:
                        out.write(
                            f"{timestep},{ang_id},{a1},{a2},{a3},"
                            f"{theta:.8f}\n"
                        )

                timesteps.append(timestep)
                angles_per_frame.append(frame_angles)

                frame_count += 1
                if max_frames is not None and frame_count >= max_frames:
                    break
    finally:
        if out is not None:
            out.close()
            print(f"Wrote {output_file}")

    return {
        "angles_topology": angles_top,
        "atom_type_names": atom_type_names,
        "timesteps": np.array(timesteps, dtype=int),
        "angles": np.array(angles_per_frame),  # shape: (n_frames, n_angles)
    }


def _parse_dihedrals_from_data(data_file, dihedral_type=None):
    """
    Parse a LAMMPS data file and return dihedrals (optionally filtered by type)
    along with a mapping from atom_id -> atom type name (from Masses comments).
    """
    type_id_to_name = {}
    atom_type_ids = {}
    atom_type_names = {}
    dihedrals = []

    section = None
    section_names = {
        "Masses", "Atoms", "Bonds", "Angles", "Dihedrals", "Impropers",
        "Velocities", "Pair Coeffs", "Bond Coeffs", "Angle Coeffs",
        "Dihedral Coeffs", "Improper Coeffs",
    }

    with open(data_file, "r", encoding="utf-8") as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue

            no_comment = stripped.split("#", 1)[0].strip()
            comment = stripped.split("#", 1)[1].strip() if "#" in stripped else ""

            if no_comment in section_names:
                section = no_comment
                continue

            if section == "Masses":
                cols = no_comment.split()
                if len(cols) >= 2:
                    type_id = int(cols[0])
                    type_id_to_name[type_id] = (
                        comment.split()[0] if comment else str(type_id)
                    )

            elif section == "Atoms":
                cols = no_comment.split()
                if len(cols) >= 7:
                    atom_id = int(cols[0])
                    atom_type_ids[atom_id] = int(cols[2])

            elif section == "Dihedrals":
                cols = no_comment.split()
                if len(cols) >= 6:
                    dih_id = int(cols[0])
                    dih_type = int(cols[1])
                    a1, a2, a3, a4 = map(int, cols[2:6])
                    if dihedral_type is None or dih_type == dihedral_type:
                        dihedrals.append((dih_id, a1, a2, a3, a4))

    for atom_id, type_id in atom_type_ids.items():
        atom_type_names[atom_id] = type_id_to_name.get(type_id, str(type_id))

    return dihedrals, atom_type_names


def _torsion_angle(p1, p2, p3, p4):
    """Signed torsion angle (degrees) for four points."""
    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3

    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)

    n1_norm = np.linalg.norm(n1)
    n2_norm = np.linalg.norm(n2)
    b2_norm = np.linalg.norm(b2)

    if n1_norm == 0 or n2_norm == 0 or b2_norm == 0:
        return np.nan

    n1 /= n1_norm
    n2 /= n2_norm
    m1 = np.cross(n1, b2 / b2_norm)

    x = np.dot(n1, n2)
    y = np.dot(m1, n2)
    return np.degrees(np.arctan2(y, x))


def calc_dihedrals(
    data_file,
    dump_file,
    dihedral_type,
    output_file=None,
    max_frames=None,
    start_timestep=0,
):
    """
    Compute torsion angles for all dihedrals of a given type across a LAMMPS dump.

    Parameters
    ----------
    data_file : str or Path
        LAMMPS data file containing the topology (Masses, Atoms, Dihedrals).
    dump_file : str or Path
        LAMMPS trajectory file (must contain id and either xu/yu/zu or x/y/z).
    dihedral_type : int
        Dihedral type ID to extract.
    output_file : str or Path or None, optional
        If provided, results are streamed to this CSV file as they are computed.
        Columns: timestep, dihedral_id, atom1, atom2, atom3, atom4, angle_deg.
    max_frames : int or None, optional
        Maximum number of frames to process. None = all frames.
    start_timestep : int, optional
        Skip frames whose timestep is below this value.

    Returns
    -------
    results : dict
        {
            "dihedrals": list of (dih_id, a1, a2, a3, a4),
            "atom_type_names": {atom_id: type_name},
            "timesteps": np.ndarray of shape (n_frames,),
            "angles": np.ndarray of shape (n_frames, n_dihedrals) in degrees,
        }
    """
    data_file = Path(data_file)
    dump_file = Path(dump_file)

    dihedrals, atom_type_names = _parse_dihedrals_from_data(
        data_file, dihedral_type=dihedral_type
    )

    if not dihedrals:
        raise ValueError(
            f"No dihedrals of type {dihedral_type} found in {data_file}."
        )

    timesteps = []
    angles_per_frame = []

    out = None
    if output_file is not None:
        out = open(output_file, "w", encoding="utf-8")
        out.write("timestep,dihedral_id,atom1,atom2,atom3,atom4,angle_deg\n")

    try:
        with open(dump_file, "r", encoding="utf-8") as f:
            frame_count = 0

            while True:
                line = f.readline()
                if not line:
                    break
                if not line.startswith("ITEM: TIMESTEP"):
                    continue

                timestep = int(f.readline().strip())

                f.readline()                       # ITEM: NUMBER OF ATOMS
                n_atoms = int(f.readline().strip())

                f.readline(); f.readline(); f.readline(); f.readline()  # box bounds

                atom_header = f.readline().strip().split()
                columns = atom_header[2:]
                id_col = columns.index("id")

                if all(c in columns for c in ["xu", "yu", "zu"]):
                    x_col = columns.index("xu")
                    y_col = columns.index("yu")
                    z_col = columns.index("zu")
                elif all(c in columns for c in ["x", "y", "z"]):
                    x_col = columns.index("x")
                    y_col = columns.index("y")
                    z_col = columns.index("z")
                else:
                    raise ValueError(
                        "Dump must contain either xu yu zu or x y z."
                    )

                # Read atoms (skip storing if we will discard frame)
                skip_frame = timestep < start_timestep
                coords = {}
                for _ in range(n_atoms):
                    cols = f.readline().split()
                    if skip_frame:
                        continue
                    atom_id = int(cols[id_col])
                    coords[atom_id] = np.array(
                        [float(cols[x_col]),
                         float(cols[y_col]),
                         float(cols[z_col])]
                    )

                if skip_frame:
                    continue

                frame_angles = np.empty(len(dihedrals))
                for i, (dih_id, a1, a2, a3, a4) in enumerate(dihedrals):
                    angle = _torsion_angle(
                        coords[a1], coords[a2], coords[a3], coords[a4]
                    )
                    frame_angles[i] = angle

                    if out is not None:
                        out.write(
                            f"{timestep},{dih_id},{a1},{a2},{a3},{a4},"
                            f"{angle:.8f}\n"
                        )

                timesteps.append(timestep)
                angles_per_frame.append(frame_angles)

                frame_count += 1
                if max_frames is not None and frame_count >= max_frames:
                    break
    finally:
        if out is not None:
            out.close()
            print(f"Wrote {output_file}")

    return {
        "dihedrals": dihedrals,
        "atom_type_names": atom_type_names,
        "timesteps": np.array(timesteps, dtype=int),
        "angles": np.array(angles_per_frame),  # shape: (n_frames, n_dihedrals)
    }

### MANIPULATIONS OF LAMMPS FILES ###

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