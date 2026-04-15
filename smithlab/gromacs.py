""" """

import subprocess
import os


def lammps_dihedrals(lammps_in, lammps_out, update_title=True):
    """
    removes the dihedral data from the LAMMPS data file
    """
    with open(lammps_in, "r", encoding="utf-8") as f:
        lines = f.readlines()

    in_dihedral = False

    with open(lammps_out, "w", encoding="utf-8") as f:

        for i, line in enumerate(lines):

            stripped = line.strip()
            columns = stripped.split()

            if not columns and not in_dihedral:
                f.writelines("\n")
                continue
            elif not columns or "dihedral" in columns or "dihedrals" in columns:
                continue

            if stripped.startswith("Dihedral"):
                in_dihedral = True
                continue

            if in_dihedral:
                try:
                    float(columns[0])
                except ValueError:
                    in_dihedral = False

            if not in_dihedral:
                if update_title and i == 0:
                    f.writelines("lmp_system\n")
                else:
                    f.writelines(line)


def intermol(intermol, lmp_data, lmp_in, pair_style, dihedral_remove=True, fix=True):
    """
    runs intermol on the data file
    """

    if dihedral_remove:
        os.rename(lmp_data, "lmp_in_temp_renamed.lmps")
        lammps_dihedrals("lmp_in_temp_renamed.lmps", lmp_data)

    cmd = (f'python3 {intermol} --lmp_in {lmp_in} --gromacs -ls "{pair_style}"')

    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    print("===== STDOUT =====")
    print(result.stdout)

    print("===== STDERR =====")
    print(result.stderr)

    if result.returncode != 0:
        print(f"Command failed with return code {result.returncode}")

    if dihedral_remove:
        # switching names back
        os.rename(lmp_data, f"{lmp_data}_no_dihedrals.lmps")
        os.rename("lmp_in_temp_renamed.lmps", f"{lmp_data}")


# intermol does a lot of the work but there are some things we still need to iron out
def fix_harmonic_bonds(top_in, top_out):
    """
    fixes the harmonic potential
    adds the fourier dihedral
    """

    with open(top_in, "r", encoding="utf-8") as f:
        lines = f.readlines()

    in_bonds = False
    skip_next = False
    with open(top_out, "w", encoding="utf-8") as f:

        for i, line in enumerate(lines):

            stripped = line.strip()
            columns = stripped.split()

            if not columns:
                in_bonds = False
                f.writelines("\n")
                continue

            if stripped == "[ bonds ]":
                in_bonds = True
                skip_next = True
                f.writelines(line)
                continue

            if in_bonds and skip_next:
                skip_next = False
                f.writelines(line)
                continue

            if in_bonds:
                columns[4] = float(columns[4]) * 2
                f.writelines(
                    f"     {columns[0]:<8}{columns[1]:<8}{columns[2]:<8}{columns[3]:<17}{columns[4]:<15.8e}\n"
                )
            else:
                f.writelines(line)


def add_fourier_dihedrals(lammps_in, top_in):

    with open(lammps_in, "r", encoding="utf-8") as f:
        lines = f.readlines()

    in_dihedral_coeffs = False
    in_dihedrals = False

    dihedral_coeffs = {}

    specific_dihedral_types = []
    atom_1 = []
    atom_2 = []
    atom_3 = []
    atom_4 = []

    for i, line in enumerate(lines):

        stripped = line.strip()
        columns = stripped.split()

        if not columns:
            continue
        elif stripped.startswith("Dihedral Coeffs"):
            in_dihedral_coeffs = True
            in_dihedrals = False
            continue
        elif stripped.startswith("Atoms"):
            in_dihedral_coeffs = False
            continue
        elif stripped.startswith("Dihedrals"):
            in_dihedral_coeffs = False
            in_dihedrals = True
            continue
        elif stripped.startswith("Impropers"):
            in_dihedrals = False
            continue

        if in_dihedral_coeffs:
            dihedral_type = int(columns[0])
            m = int(columns[1])
            expected_len = 2 + 3 * m
            if len(columns) < expected_len:
                raise ValueError(f"Dihedral type {dihedral_type} says it has {m} Fourier terms, but {len(columns) - 2} coefficient fields were found.")
            
            coeffs = []
            for term_idx in range(m):
                base = 2 + 3 * term_idx
                K = float(columns[base])
                n = int(float(columns[base + 1]))
                d = float(columns[base + 2])
                coeffs.append((K, n, d))

            dihedral_coeffs[dihedral_type] = coeffs

        if in_dihedrals:
            specific_dihedral_types.append(int(columns[1]))
            atom_1.append(columns[2])
            atom_2.append(columns[3])
            atom_3.append(columns[4])
            atom_4.append(columns[5])

    with open(top_in, "r", encoding="utf-8") as f:
        lines = f.readlines()

    insert_index = next((i for i, line in enumerate(lines) if "[ system ]" in line), None)

    if insert_index is None:
        raise ValueError("Could not find '[ system ]' in the topology file.")

    inserted_lines = ["[ dihedrals ]\n;   ai   aj   ak   al   funct   phi   K     mult.\n"]

    for i, specific_d_type in enumerate(specific_dihedral_types):
        if specific_d_type not in dihedral_coeffs:
            raise ValueError(f"Dihedral type {specific_d_type} found in dihedrals section but not in dihedral coefficients")

        # no factor of 2 used here
        for K, n, d in dihedral_coeffs[specific_d_type]:
            line = (
                f"     {atom_1[i]:<6}{atom_2[i]:<6}{atom_3[i]:<6}{atom_4[i]:<6}"
                f"  1  {d:<18.10f}{K * 4.184:<18.10f} {n}\n"
            )
            inserted_lines.append(line)

    inserted_lines.append("\n")
    lines[insert_index:insert_index] = inserted_lines

    with open(top_in, "w", encoding="utf-8") as f:
        f.writelines(lines)
