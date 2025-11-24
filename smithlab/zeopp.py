"""
Copyright 2025. Brandon C. Tapia

MIT License
"""

import os
import numpy as np

def write_cuc(lammps_in, cuc_out):
    """
    Writes an CUC file from a LAMMPS data file (atom_style full),
    formatted for use with Zeo++.
    """
    atom_section = False
    atom_lines = []

    with open(lammps_in, "r", encoding="utf-8") as file:
        lines = file.readlines()

    for i, line in enumerate(lines):
        stripped = line.strip()
        columns = stripped.split()

        if stripped.endswith("xhi"):
            x_start = float(columns[0])
            x_length = float(columns[1]) - float(columns[0])

        elif stripped.endswith("yhi"):
            y_start = float(columns[0])
            y_length = float(columns[1]) - float(columns[0])

        elif stripped.endswith("zhi"):
            z_start = float(columns[0])
            z_length = float(columns[1]) - float(columns[0])

        elif stripped.startswith("Atoms"):
            atom_section = True
            continue

        if atom_section:
            if not columns:
                continue
            if columns[0].isalpha():
                break
            atom_lines.append(columns)

    with open(cuc_out, "w", encoding="utf-8") as file:
        file.write("Processing: file_from_smithlab.zeopp.write_cuc\n")
        file.write(f"Unit_cell: {x_length} {y_length} {z_length} 90 90 90\n")
        #print(atom_lines)
        for cols in atom_lines:
            atom_type = cols[2]
            x, y, z = float(cols[4]), float(cols[5]), float(cols[6])
            shift_x = x - x_start
            shift_y = y - y_start
            shift_z = z - z_start
            fractional_x = shift_x / x_length
            fractional_y = shift_y / y_length
            fractional_z = shift_z / z_length
            file.write(f"A{atom_type} {fractional_x} {fractional_y} {fractional_z}\n")


def write_radfile(lammps_in, radfile_out):

    pair_section = False
    identifiers = []
    sigma = []

    with open(lammps_in, "r", encoding="utf-8") as file:
        lines = file.readlines()

    for i, line in enumerate(lines):
        stripped = line.strip()
        columns = stripped.split()

        if stripped.startswith("Pair Coeffs"):
            pair_section = True
            continue

        if pair_section:
            if not columns:
                continue
            if columns[0].isalpha():
                break
            identifiers.append(columns[0])
            sigma.append(columns[2])

    with open(radfile_out, "w", encoding="utf-8") as file:
        for i, identifier in enumerate(identifiers):
            file.write(f"A{identifier} {float(sigma[i])/2}\n") # IS THIS THE SIGMA WE CARE ABOUT??
        file.write("Si 1.35")


def write_massfile(lammps_in, massfile_out):

    mass_section = False
    identifiers = []
    molwt = []

    with open(lammps_in, "r", encoding="utf-8") as file:
        lines = file.readlines()

    for i, line in enumerate(lines):
        stripped = line.strip()
        columns = stripped.split()

        if stripped.startswith("Masses"):
            mass_section = True
            continue

        if mass_section:
            if not columns:
                continue
            if columns[0].isalpha():
                break

            identifiers.append(columns[0])
            molwt.append(columns[1])

    with open(massfile_out, "w", encoding="utf-8") as file:
        for i, identifier in enumerate(identifiers):
            file.write(f"A{identifier} {molwt[i]}\n")
        file.write("Si 28.0855")

def setup_zeopp(lammps_in,
                cuc_out="system.cuc",
                radfile_out="radii.rad",
                massfile_out="molwt.mass"):
    write_cuc(lammps_in, cuc_out)
    write_radfile(lammps_in, radfile_out)
    write_massfile(lammps_in, massfile_out)

def zeopp_command(
        zeopp_loc="./network",
        cuc_file="system.cuc")
        rad_file="radii.rad",
        mass_file="molwt.mass",
        visvcoro="0.2"


# we want to call zeopp as:
# ./network -r radii.rad -mass molwt.mass -resex "resex.out" -visVoro 0.2
# -zvis "zeovis.out" -axs 1.8 "axs.out"

#/home/gridsan/btapia/zeo++-0.3/network -ha -r radii.rad -mass molwt.mass -visVoro 0.2 system.cuc