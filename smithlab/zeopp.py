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
        cuc_file="system.cuc",
        rad_file="radii.rad",
        mass_file="molwt.mass",
        visvcoro="0.2"):
    return


# we want to call zeopp as:
# ./network -r radii.rad -mass molwt.mass -resex "resex.out" -visVoro 0.2
# -zvis "zeovis.out" -axs 1.8 "axs.out"

#/home/gridsan/btapia/zeo++-0.3/network -ha -r radii.rad -mass molwt.mass -visVoro 0.2 system.cuc


def cube2xyz(cube_in, xyz_out, d_spacing=None, d_min=0.0):

    with open(cube_in, "r") as f:

        _ = f.readline()
        _ = f.readline()

        # Line 3: natoms, origin
        parts = f.readline().split()
        natoms = int(parts[0])
        origin = np.array(list(map(float, parts[1:4])), dtype=float)

        # Line 4: nx and x-axis vector
        parts = f.readline().split()
        nx = int(parts[0])
        ax = np.array(list(map(float, parts[1:4])), dtype=float)

        # Line 5: ny and y-axis vector
        parts = f.readline().split()
        ny = int(parts[0])
        by = np.array(list(map(float, parts[1:4])), dtype=float)

        # Line 6: nz and z-axis vector
        parts = f.readline().split()
        nz = int(parts[0])
        cz = np.array(list(map(float, parts[1:4])), dtype=float)

        # Skip atom lines
        for _ in range(natoms):
            f.readline()

        # Read all remaining volumetric values
        vals = []
        for line in f:
            vals.extend(line.split())

        vals = np.array(list(map(float, vals)), dtype=float)
        ngrid = nx * ny * nz

        if vals.size < ngrid:
            raise ValueError(f"Not enough grid values: expected {ngrid}, got {vals.size}")
        if vals.size > ngrid:
            #raise ValueError(f"Too many grid values: expected {ngrid}, got {vals.size}")
            vals = vals[:ngrid]

    # Reshape to (nx, ny, nz); cube order is x outer, y middle, z inner
    vals = vals.reshape((nx, ny, nz))

    # ---- compute strides for downsampling ----
    # Native step lengths (norm of grid vectors)
    step_x = np.linalg.norm(ax)
    step_y = np.linalg.norm(by)
    step_z = np.linalg.norm(cz)

    # Desired spacings: if None, keep native grid (stride = 1)
    def compute_stride(desired, native):
        if desired is None or desired <= 0:
            return 1
        factor = int(round(desired / native))
        if factor < 1:
            factor = 1
        return factor

    dx = d_spacing
    dy = d_spacing
    dz = d_spacing

    stride_x = compute_stride(dx, step_x)
    stride_y = compute_stride(dy, step_y)
    stride_z = compute_stride(dz, step_z)

    # Index arrays we will keep
    ix_idx = np.arange(0, nx, stride_x, dtype=int)
    iy_idx = np.arange(0, ny, stride_y, dtype=int)
    iz_idx = np.arange(0, nz, stride_z, dtype=int)

    # Subsample distances without building the full big coords array first
    vals_sub = vals[np.ix_(ix_idx, iy_idx, iz_idx)]  # shape (nx', ny', nz')

    # Build coords only for kept indices
    ix_grid, iy_grid, iz_grid = np.meshgrid(ix_idx, iy_idx, iz_idx, indexing="ij")

    # r = origin + ix*ax + iy*by + iz*cz
    coords_sub = (
        origin[None, None, None, :]
        + ix_grid[..., None] * ax[None, None, None, :]
        + iy_grid[..., None] * by[None, None, None, :]
        + iz_grid[..., None] * cz[None, None, None, :]
    )  # (nx', ny', nz', 3)

    # Flatten
    coords_flat = coords_sub.reshape(-1, 3)
    vals_flat = vals_sub.reshape(-1)

    # ---- FILTER: keep only distances >= dmin ----
    mask = vals_flat >= d_min
    coords_keep = coords_flat[mask]
    vals_keep = vals_flat[mask][:, None]

    data = np.hstack((coords_keep, vals_keep))

    header = "x y z distance"
    np.savetxt(xyz_out, data, fmt="%.6f", header=header)
    return
