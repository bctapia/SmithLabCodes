"""
Copyright 2025 Brandon C. Tapia

MIT License
"""

from collections import defaultdict
from pathlib import Path
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