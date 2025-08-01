"""
"""
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


def intermol(intermol, lmp_in, pair_style, dihedral_remove=True, fix=True):
    """
    runs intermol on the data file
    """

    if dihedral_remove:
        os.rename(lmp_in, "lmp_in_temp_renamed.lmps")
        lammps_dihedrals("lmp_in_temp_renamed.lmps", lmp_in)

    cmd = (
        'python3 /mnt/c/Users/btapi/OneDrive/Documents/GitHub/InterMol/intermol/convert.py '
        '--lmp_in /home/btapia/intermol_test/equil.in --gromacs ' # TODO generalize
        '-ls "pair_style cut/coul/long 15"'  # TODO generalize
    )
    
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    print("===== STDOUT =====")
    print(result.stdout)

    print("===== STDERR =====")
    print(result.stderr)

    if result.returncode != 0:
        print(f"Command failed with return code {result.returncode}")

    if dihedral_remove:
        # switching names back
        os.rename(lmp_in, f"{lmp_in}_no_dihedrals.lmps")
        os.rename("lmp_in_temp_renamed.lmps", f"{lmp_in}") 
        
        

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
                columns[4] = float(columns[4])*2
                f.writelines(f"     {columns[0]:<8}{columns[1]:<8}{columns[2]:<8}{columns[3]:<17}{columns[4]:<15.8e}\n")
            else:
                f.writelines(line)



def add_fourier_dihedrals(lammps_in, top_in):

    with open(lammps_in, "r", encoding="utf-8") as f:
        lines = f.readlines()

    in_dihedral_coeffs = False
    in_dihedrals = False

    dihedral_types = []
    dihedral_K = []
    dihedral_n1 = []
    dihedral_d1 = []

    specific_dihedral_index = []
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
            continue
        elif stripped.startswith("Atoms"):
            in_dihedral_coeffs = False
            continue
        elif stripped.startswith("Dihedrals"):
            in_dihedrals = True
            continue
        elif stripped.startswith("Impropers"):
            in_dihedrals = False
            continue

        if in_dihedral_coeffs:
            if int(columns[1]) != 1:
                print("sorry, only funct=1 is supported right now, cannot continue :(")
                break

            dihedral_types.append(int(columns[0]))
            dihedral_K.append(float(columns[2]))
            dihedral_n1.append(float(columns[3]))
            dihedral_d1.append(float(columns[4]))

        if in_dihedrals:
            specific_dihedral_index.append(columns[0])
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

        index_use = None
        # look through the rtype array:
        for j, d_type in enumerate(dihedral_types):
            if specific_d_type == d_type:
                index_use = j
                break
        
        # no factor of 2 used here
        line = f"     {atom_1[i]:<6}{atom_2[i]:<6}{atom_3[i]:<6}{atom_4[i]:<6}  1  {dihedral_d1[index_use]:<.15} {dihedral_K[index_use]*4.184:<.15} {int(dihedral_n1[index_use])}\n"
        inserted_lines.append(line)

    lines[insert_index:insert_index] = inserted_lines

    with open(top_in, "w", encoding="utf-8") as f:
        f.writelines(lines)