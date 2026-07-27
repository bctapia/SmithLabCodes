"""
Copyright 2026. Brandon C. Tapia

MIT License
"""

import os
import numpy as np
import matplotlib.pyplot as plt


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
        f.readline()
        f.readline()

        parts = f.readline().split()
        natoms = int(parts[0])
        origin = np.array(parts[1:4], dtype=np.float32)

        parts = f.readline().split()
        nx = int(parts[0])
        ax = np.array(parts[1:4], dtype=np.float32)

        parts = f.readline().split()
        ny = int(parts[0])
        by = np.array(parts[1:4], dtype=np.float32)

        parts = f.readline().split()
        nz = int(parts[0])
        cz = np.array(parts[1:4], dtype=np.float32)

        for _ in range(natoms):
            f.readline()

        ngrid = nx * ny * nz

        # ---- Stream values into a preallocated float32 array ----
        vals = np.empty(ngrid, dtype=np.float32)
        i = 0
        for line in f:
            tokens = line.split()
            if not tokens:
                continue
            n = len(tokens)
            if i + n > ngrid:
                n = ngrid - i
                if n <= 0:
                    break
            # np.fromstring is deprecated; use this instead
            vals[i:i + n] = np.asarray(tokens[:n], dtype=np.float32)
            i += n
            if i >= ngrid:
                break

        if i < ngrid:
            raise ValueError(f"Not enough grid values: expected {ngrid}, got {i}")

    # Reshape (view, no copy)
    vals = vals.reshape((nx, ny, nz))

    # ---- Strides for downsampling ----
    step_x = float(np.linalg.norm(ax))
    step_y = float(np.linalg.norm(by))
    step_z = float(np.linalg.norm(cz))

    def compute_stride(desired, native):
        if desired is None or desired <= 0:
            return 1
        return max(1, int(round(desired / native)))

    stride_x = compute_stride(d_spacing, step_x)
    stride_y = compute_stride(d_spacing, step_y)
    stride_z = compute_stride(d_spacing, step_z)

    # Subsample (this is a view if strides divide evenly, else a copy of just the kept slab)
    vals_sub = vals[::stride_x, ::stride_y, ::stride_z]
    nx_s, ny_s, nz_s = vals_sub.shape

    # Free the full-resolution array now that we have the subsample
    del vals

    # ---- Mask FIRST, then compute coords only for kept points ----
    vals_flat = vals_sub.ravel()
    del vals_sub

    mask = vals_flat >= d_min
    n_keep = int(mask.sum())

    if n_keep == 0:
        np.savetxt(xyz_out, np.empty((0, 4)), fmt="%.6f", header="x y z distance")
        return

    kept_vals = vals_flat[mask]
    kept_lin = np.flatnonzero(mask).astype(np.int64)
    del vals_flat, mask

    # Recover (i, j, k) indices in the subsampled grid
    ii_s, jj_s, kk_s = np.unravel_index(kept_lin, (nx_s, ny_s, nz_s))
    del kept_lin

    # Map back to original-grid indices
    ii = (ii_s.astype(np.float32)) * (stride_x)
    jj = (jj_s.astype(np.float32)) * (stride_y)
    kk = (kk_s.astype(np.float32)) * (stride_z)
    del ii_s, jj_s, kk_s

    # Compute coords for kept points only:  r = origin + i*ax + j*by + k*cz
    # Build column-by-column to avoid a (n_keep, 3) intermediate before stacking
    out = np.empty((n_keep, 4), dtype=np.float32)
    out[:, 0] = origin[0] + ii * ax[0] + jj * by[0] + kk * cz[0]
    out[:, 1] = origin[1] + ii * ax[1] + jj * by[1] + kk * cz[1]
    out[:, 2] = origin[2] + ii * ax[2] + jj * by[2] + kk * cz[2]
    out[:, 3] = kept_vals

    np.savetxt(xyz_out, out, fmt="%.6f", header="x y z distance")

def plot_cube(xyz_in, keep_every=1, d_lo=None, d_hi=None, vmin=None, vmax=None, png_out=None):
    if keep_every < 1:
        raise ValueError("keep_every must be >= 1")
    if d_lo is not None and d_hi is not None and d_lo > d_hi:
        raise ValueError("d_lo cannot be greater than d_hi")

    data = np.loadtxt(xyz_in, comments="#", dtype=np.float32)
    x, y, z, d = data.T

    xmin, xmax = np.min(x), np.max(x)
    ymin, ymax = np.min(y), np.max(y)
    zmin, zmax = np.min(z), np.max(z)

    px = 0.02 * (xmax - xmin if xmax > xmin else 1.0)
    py = 0.02 * (ymax - ymin if ymax > ymin else 1.0)
    pz = 0.02 * (zmax - zmin if zmax > zmin else 1.0)

    xmin, xmax = xmin - px, xmax + px
    ymin, ymax = ymin - py, ymax + py
    zmin, zmax = zmin - pz, zmax + pz


    # distance filter
    mask = np.ones(d.shape, dtype=bool)
    if d_lo is not None:
        mask &= (d >= d_lo)
    if d_hi is not None:
        mask &= (d <= d_hi)

    x, y, z, d = x[mask], y[mask], z[mask], d[mask]

    if len(d) == 0:
        print("No points match the distance range.")
        return

    # downsample
    idx = np.arange(0, len(x), keep_every)
    x_t, y_t, z_t, d_t = x[idx], y[idx], z[idx], d[idx]

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection="3d")

    sc = ax.scatter(x_t, y_t, z_t, c=d_t, cmap="viridis", s=2, vmin=vmin, vmax=vmax, alpha=0.8, antialiaseds=True, linewidths=0)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_zlim(zmin, zmax)

    # remove gray panes/background + grid
    ax.grid(False)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.fill = False
        axis.pane.set_facecolor((1, 1, 1, 0))   # transparent
        axis.pane.set_edgecolor((1, 1, 1, 0))   # no pane edge

    ax.tick_params(axis='x', which='both', direction='out', length=0)
    ax.tick_params(axis='y', which='both', direction='out', length=0)
    ax.tick_params(axis='z', which='both', direction='out', length=0)

    ax.set_box_aspect((1, 1, 1))   # cubic axes box
    # optional: removes perspective distortion so cube looks more "true"
    #ax.set_proj_type('ortho')

    corners = np.array([
        [xmin, ymin, zmin], [xmax, ymin, zmin],
        [xmax, ymax, zmin], [xmin, ymax, zmin],
        [xmin, ymin, zmax], [xmax, ymin, zmax],
        [xmax, ymax, zmax], [xmin, ymax, zmax],
    ])

    edges = [
        (0,1), (1,2), (2,3), (3,0),  # bottom
        (4,5), (5,6), (6,7), (7,4),  # top
        (0,4), (1,5), (2,6), (3,7)   # verticals
    ]

    for i, j in edges:
        ax.plot(
            [corners[i,0], corners[j,0]],
            [corners[i,1], corners[j,1]],
            [corners[i,2], corners[j,2]],
            color="black", lw=1.0, antialiased=True
        )


    cbar = plt.colorbar(sc, ax=ax, pad=0.1)
    cbar.set_label("Distance to closest atom (Å)")

    #ax.set_xlabel("x")
    #ax.set_ylabel("y")
    #ax.set_zlabel("z")

    plt.tight_layout()
    if png_out:
        plt.savefig(png_out, dpi=600, bbox_inches="tight")  # export DPI
    #plt.show()
