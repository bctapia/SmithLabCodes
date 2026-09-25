"""
Copyright 2026. Brandon C. Tapia

MIT License
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
from collections import deque


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

        for _ in range(abs(natoms)):
            f.readline()

        # Zeo++ emits a dataset line ("1    1") after the atom block, which the
        # cube spec only prescribes when natoms < 0. Sniff for it: the dataset
        # line is bare integers, data lines are E-formatted floats.
        pos = f.tell()
        tokens = f.readline().split()
        is_dataset_line = (
            0 < len(tokens) <= 10
            and all(t.lstrip("+-").isdigit() for t in tokens)
        )
        if not is_dataset_line:
            f.seek(pos)

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
                raise ValueError(
                    f"Too many grid values: expected {ngrid}, overflowed at "
                    f"index {i} with {n} more on one line. Header/data mismatch."
                )
            vals[i:i + n] = np.asarray(tokens, dtype=np.float32)
            i += n

        if i != ngrid:
            raise ValueError(f"Grid value count mismatch: expected {ngrid}, got {i}")

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

    vals_sub = vals[::stride_x, ::stride_y, ::stride_z]
    nx_s, ny_s, nz_s = vals_sub.shape
    del vals

    # ---- Mask FIRST, then compute coords only for kept points ----
    vals_flat = np.ascontiguousarray(vals_sub).ravel()
    del vals_sub

    mask = vals_flat >= d_min
    n_keep = int(mask.sum())

    if n_keep == 0:
        np.savetxt(xyz_out, np.empty((0, 4)), fmt="%.6f", header="x y z distance")
        return

    kept_vals = vals_flat[mask]
    kept_lin = np.flatnonzero(mask).astype(np.int64)
    del vals_flat, mask

    ii_s, jj_s, kk_s = np.unravel_index(kept_lin, (nx_s, ny_s, nz_s))
    del kept_lin

    ii = ii_s.astype(np.float64) * stride_x
    jj = jj_s.astype(np.float64) * stride_y
    kk = kk_s.astype(np.float64) * stride_z
    del ii_s, jj_s, kk_s

    out = np.empty((n_keep, 4), dtype=np.float64)
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


# ============================================================
# 1. READ ZEO++ DISTANCE GRID
# ============================================================
def cube2npy(cube_in, npy_out, meta_out=None):
    """
    Read a Zeo++ distance-grid cube file and save the full field as a
    3D float32 array.

    Assumes Zeo++ -gridG, i.e. distances/grid vectors are in Angstrom.

    Returns
    -------
    vals : ndarray, shape (nx, ny, nz)
        Distance to nearest atomic surface in Angstrom.

    meta : dict
        Grid/cell information.
    """

    with open(cube_in, "r") as f:
        f.readline()
        f.readline()

        parts = f.readline().split()
        natoms = int(parts[0])
        origin = np.array(parts[1:4], dtype=np.float64)

        raw_n = np.empty(3, dtype=np.int64)
        n = np.empty(3, dtype=np.int64)

        # Rows are voxel vectors
        vox = np.empty((3, 3), dtype=np.float64)

        for ax in range(3):
            parts = f.readline().split()

            raw_n[ax] = int(parts[0])
            n[ax] = abs(raw_n[ax])

            vox[ax] = np.array(
                parts[1:4],
                dtype=np.float64
            )

        # Atom block
        for _ in range(abs(natoms)):
            f.readline()

        # Zeo++ may write a dataset line such as:
        #
        #     1    1
        #
        # Skip it if present.
        pos = f.tell()
        tokens = f.readline().split()

        if not (
            0 < len(tokens) <= 10
            and all(t.lstrip("+-").isdigit() for t in tokens)
        ):
            f.seek(pos)

        ngrid = int(n.prod())

        vals = np.empty(
            ngrid,
            dtype=np.float32
        )

        i = 0

        for line in f:
            tokens = line.split()

            if not tokens:
                continue

            k = len(tokens)

            if i + k > ngrid:
                raise ValueError(
                    f"Too many grid values: expected {ngrid}, "
                    f"overflowed at index {i} with {k} more."
                )

            vals[i:i + k] = np.asarray(
                tokens,
                dtype=np.float32
            )

            i += k

        if i != ngrid:
            raise ValueError(
                f"Expected {ngrid} grid values, got {i}"
            )

    # Cube convention: last axis varies fastest
    vals = vals.reshape(tuple(n))

    # Full lattice vectors
    cell = vox * n[:, None]

    meta = {
        "n": n,
        "raw_n": raw_n,
        "vox": vox,
        "cell": cell,
        "origin": origin,
        "spacing": np.linalg.norm(vox, axis=1),
        "voxel_volume": abs(np.linalg.det(vox)),
        "cell_volume": abs(np.linalg.det(cell)),
        "natoms": abs(natoms),
    }

    np.save(npy_out, vals)

    if meta_out is not None:
        np.savez(meta_out, **meta)

    return vals, meta



# ============================================================
# 2. UNION-FIND
# ============================================================

def _find(parent, x):
    root = x

    while parent[root] != root:
        root = parent[root]

    # Path compression
    while parent[x] != root:
        parent[x], x = root, parent[x]

    return root


def _union(parent, a, b):
    ra = _find(parent, a)
    rb = _find(parent, b)

    if ra != rb:
        parent[max(ra, rb)] = min(ra, rb)

# ============================================================
# 3. PERIODIC WRAPPING
# ============================================================

def _detect_wrapping(n_initial, periodic_edges, roots):
    """
    Determine whether each periodic component has non-zero winding
    around each lattice direction.

    Parameters
    ----------
    n_initial : int
        Number of components before periodic merging.

    periodic_edges : list
        Entries are (u, v, axis).

        u -> v crosses +1 periodic image along `axis`.

    roots : ndarray
        Union-find root for every original component.

    Returns
    -------
    wraps_root : ndarray, shape (n_initial + 1, 3)
        True if the periodic component winds around lattice direction
        a, b, or c.

    Notes
    -----
    This detects actual winding. Merely touching opposite faces is
    not sufficient.
    """

    adj = [
        [] for _ in range(n_initial + 1)
    ]

    for u, v, axis in periodic_edges:

        shift = np.zeros(
            3,
            dtype=np.int32
        )

        shift[axis] = 1

        adj[u].append(
            (v, shift)
        )

        adj[v].append(
            (u, -shift)
        )

    # Periodic-image coordinate assigned to each initial component
    image = np.zeros(
        (n_initial + 1, 3),
        dtype=np.int32
    )

    assigned = np.zeros(
        n_initial + 1,
        dtype=bool
    )

    wraps_root = np.zeros(
        (n_initial + 1, 3),
        dtype=bool
    )

    for start in range(1, n_initial + 1):

        if assigned[start]:
            continue

        if not adj[start]:
            assigned[start] = True
            continue

        root = roots[start]

        assigned[start] = True
        image[start] = 0

        queue = deque([start])

        while queue:

            u = queue.popleft()

            for v, shift in adj[u]:

                expected = image[u] + shift

                if not assigned[v]:

                    image[v] = expected
                    assigned[v] = True
                    queue.append(v)

                else:

                    # Reaching the same component through another
                    # path with a different image coordinate means
                    # the path winds around the periodic cell.
                    mismatch = expected - image[v]

                    wraps_root[root] |= (
                        mismatch != 0
                    )

    return wraps_root

def _merge_periodic_labels(labels, n_initial):
    """
    Merge nonperiodic connected-component labels across PBCs and
    determine true winding in the three lattice directions.

    Returns
    -------
    labels : ndarray
        Periodically merged labels.

    n : int
        Number of periodic pores.

    wraps : ndarray, shape (n + 1, 3)
        wraps[i] = [wrap_a, wrap_b, wrap_c]
    """

    parent = np.arange(
        n_initial + 1,
        dtype=np.int64
    )

    periodic_edges = []

    for axis in range(3):

        lo = np.take(
            labels,
            0,
            axis=axis
        ).ravel()

        hi = np.take(
            labels,
            -1,
            axis=axis
        ).ravel()

        both = (lo > 0) & (hi > 0)

        if not np.any(both):
            continue

        pairs = np.unique(
            np.stack(
                [lo[both], hi[both]],
                axis=1
            ),
            axis=0
        )

        for lo_id, hi_id in pairs:

            lo_id = int(lo_id)
            hi_id = int(hi_id)

            # Moving from high face -> low face is taken as
            # moving +1 lattice image in this direction.
            periodic_edges.append(
                (hi_id, lo_id, axis)
            )

            _union(
                parent,
                lo_id,
                hi_id
            )

    # Collapse union-find
    roots = np.array(
        [
            _find(parent, i)
            for i in range(n_initial + 1)
        ],
        dtype=np.int64
    )

    uniq = np.unique(
        roots[1:]
    )

    n = len(uniq)

    remap = np.zeros(
        n_initial + 1,
        dtype=np.int32
    )

    remap[uniq] = np.arange(
        1,
        n + 1,
        dtype=np.int32
    )

    initial_to_periodic = remap[roots]

    labels = initial_to_periodic[labels]

    # Determine actual winding
    wraps_root = _detect_wrapping(
        n_initial,
        periodic_edges,
        roots
    )

    wraps = np.zeros(
        (n + 1, 3),
        dtype=bool
    )

    for root in uniq:

        pid = remap[root]

        wraps[pid] = wraps_root[root]

    return labels, n, wraps


# ============================================================
# 4. CONNECTED PORES + GLOBAL METRICS
# ============================================================

def find_pores(vals, meta, probe_radius, periodic=True):
    """
    Identify connected probe-accessible pores.

    A voxel is accessible when

        vals >= probe_radius

    Connectivity is 6-neighbor / face connectivity.

    With periodic=True, components are joined through the periodic
    boundaries and true winding is determined.

    Returns
    -------
    labels : int32 ndarray
        0 = blocked
        1 = largest pore
        2 = second-largest pore
        ...

    stats : dict
        Global connectivity metrics and per-pore information.
    """

    open_ = vals >= probe_radius

    # Face connectivity only
    struct = ndimage.generate_binary_structure(
        rank=3,
        connectivity=1
    )

    labels, n_initial = ndimage.label(
        open_,
        structure=struct
    )

    labels = labels.astype(
        np.int32
    )

    # --------------------------------------------------------
    # Periodic merging + winding
    # --------------------------------------------------------

    if periodic and n_initial > 0:

        labels, n, wraps_by_id = _merge_periodic_labels(
            labels,
            n_initial
        )

    else:

        n = n_initial

        wraps_by_id = np.zeros(
            (n + 1, 3),
            dtype=bool
        )

    # --------------------------------------------------------
    # Sort pore IDs by volume
    # --------------------------------------------------------

    counts = np.bincount(
        labels.ravel(),
        minlength=n + 1
    )[1:]

    if n > 0:

        order = np.argsort(
            counts
        )[::-1]

        old_ids = order + 1

        relabel = np.zeros(
            n + 1,
            dtype=np.int32
        )

        relabel[old_ids] = np.arange(
            1,
            n + 1,
            dtype=np.int32
        )

        labels = relabel[labels]

        counts = counts[order]

        # Keep wrapping data in the exact same order
        wraps = wraps_by_id[old_ids]

    else:

        wraps = np.zeros(
            (0, 3),
            dtype=bool
        )

    # --------------------------------------------------------
    # Volumes
    # --------------------------------------------------------

    voxel_volume = meta["voxel_volume"]

    volumes = (
        counts.astype(np.float64)
        * voxel_volume
    )

    V_acc = volumes.sum()

    if V_acc > 0:
        omega = volumes / V_acc
    else:
        omega = np.zeros_like(volumes)

    # --------------------------------------------------------
    # Percolation
    # --------------------------------------------------------
    #
    # Columns correspond to lattice directions a,b,c.
    #
    # A pore is called percolating if it winds around ANY
    # periodic lattice direction.
    #

    percolating = wraps.any(
        axis=1
    )

    has_percolating = bool(
        np.any(percolating)
    )

    # --------------------------------------------------------
    # Global concentration / fragmentation metrics
    # --------------------------------------------------------

    phi = (
        V_acc / meta["cell_volume"]
        if meta["cell_volume"] > 0
        else 0.0
    )

    omega_1 = (
        omega[0]
        if n > 0
        else 0.0
    )

    P_abs = (
        volumes[0] / meta["cell_volume"]
        if n > 0
        else 0.0
    )

    # Simpson / HHI concentration
    lambda_ = (
        np.sum(omega ** 2)
        if n > 0
        else 0.0
    )

    N_eff = (
        1.0 / lambda_
        if lambda_ > 0
        else 0.0
    )

    # --------------------------------------------------------
    # Susceptibility
    # --------------------------------------------------------
    #
    # IMPORTANT:
    # Exclude every PERCOLATING component, not automatically
    # the largest component.
    #

    finite = ~percolating

    finite_volumes = volumes[finite]
    finite_omega = omega[finite]

    V_finite = finite_volumes.sum()

    # Dimensional susceptibility [Angstrom^3]
    chi_volume = (
        np.sum(finite_volumes ** 2)
        / V_finite
        if V_finite > 0
        else 0.0
    )

    # Dimensionless susceptibility:
    #
    # chi* = chi / V_acc
    #
    #      = sum(omega_i^2) / sum(omega_i)
    #
    # over NONPERCOLATING clusters only.
    finite_omega_sum = finite_omega.sum()

    chi_star = (
        np.sum(finite_omega ** 2)
        / finite_omega_sum
        if finite_omega_sum > 0
        else 0.0
    )

    stats = {

        # Input
        "probe_radius": float(probe_radius),

        # Number of connected pores
        "n_pores": int(n),

        # Accessible free-volume fraction
        "phi": float(phi),

        # Largest connected pore / accessible volume
        "omega_1": float(omega_1),

        # Largest connected pore / total cell volume
        "P_abs": float(P_abs),

        # Simpson / HHI concentration
        "lambda": float(lambda_),

        # Effective number of pores
        "N_eff": float(N_eff),

        # Susceptibility
        "chi_volume": float(chi_volume),
        "chi_star": float(chi_star),

        # Percolation
        "has_percolating": has_percolating,
        "n_percolating": int(
            np.count_nonzero(percolating)
        ),

        # Per-pore arrays, sorted largest -> smallest
        "volumes": volumes,
        "omega": omega,

        # columns = a,b,c lattice directions
        "wraps": wraps,

        # True if wraps in any dimension
        "percolating": percolating,
    }

    return labels, stats


# ============================================================
# 5. SCAN MANY PROBE RADII
# ============================================================

def scan_probe_radii(
    vals,
    meta,
    probe_radii,
    periodic=True
):
    """
    Evaluate the connectivity metrics over a series of probe radii.

    Returns
    -------
    results : list of dict
        One stats dictionary for each probe radius.
    """

    results = []

    for r in probe_radii:

        _, stats = find_pores(
            vals,
            meta,
            probe_radius=r,
            periodic=periodic
        )

        results.append(stats)

    return results


# ============================================================
# 6. GRID GEOMETRY HELPERS
# ============================================================

def _flat_to_ijk(index, shape):
    return tuple(
        int(x)
        for x in np.unravel_index(
            index,
            shape
        )
    )


def _ijk_to_cart(ijk, meta):
    """
    Convert grid index to Cartesian position.
    """

    ijk = np.asarray(
        ijk,
        dtype=np.float64
    )

    return (
        meta["origin"]
        + ijk[0] * meta["vox"][0]
        + ijk[1] * meta["vox"][1]
        + ijk[2] * meta["vox"][2]
    )


def _neighbors6_flat(index, shape, periodic=True):
    """
    Return flat indices of 6 face-neighbors.
    """

    nx, ny, nz = shape

    yz = ny * nz

    i = index // yz
    rem = index - i * yz

    j = rem // nz
    k = rem - j * nz

    nbrs = []

    # a / axis 0
    if i > 0:
        nbrs.append(index - yz)
    elif periodic and nx > 1:
        nbrs.append(index + (nx - 1) * yz)

    if i < nx - 1:
        nbrs.append(index + yz)
    elif periodic and nx > 1:
        nbrs.append(index - (nx - 1) * yz)

    # b / axis 1
    if j > 0:
        nbrs.append(index - nz)
    elif periodic and ny > 1:
        nbrs.append(index + (ny - 1) * nz)

    if j < ny - 1:
        nbrs.append(index + nz)
    elif periodic and ny > 1:
        nbrs.append(index - (ny - 1) * nz)

    # c / axis 2
    if k > 0:
        nbrs.append(index - 1)
    elif periodic and nz > 1:
        nbrs.append(index + nz - 1)

    if k < nz - 1:
        nbrs.append(index + 1)
    elif periodic and nz > 1:
        nbrs.append(index - nz + 1)

    # Mainly protects against pathological dimensions of 1 or 2.
    return tuple(
        dict.fromkeys(nbrs)
    )


# ============================================================
# 7. CAVITY / BOTTLENECK MERGE TREE
# ============================================================

def build_bottleneck_merge_tree(
    vals,
    meta,
    min_radius=0.0,
    periodic=True
):
    """
    Construct the 0-dimensional superlevel-set merge tree of the
    distance field.

    Interpretation
    --------------
    As radius r decreases, consider

        A(r) = {x : D(x) >= r}

    Local maxima of D create new connected components ("cavities").

    When two previously separate components first connect, their
    merger occurs at a saddle / bottleneck radius R_t.

    For two merging cavity branches:

        beta = R_t / min(R_i, R_j)

    where R_i and R_j are their cavity maximum radii.

    beta ~ 1  -> weak constriction / broad connection
    beta << 1 -> severe bottleneck

    The returned edges form a MERGE TREE. This captures hierarchical
    bottlenecks but NOT redundant alternative paths or cycles.

    Parameters
    ----------
    vals : ndarray
        Zeo++ distance field.

    meta : dict
        Grid metadata.

    min_radius : float
        Ignore voxels below this distance.

        min_radius=0 analyzes geometrical void space outside atomic
        surfaces.

    periodic : bool
        Apply 6-neighbor periodic adjacency.

    Returns
    -------
    tree : dict
        Contains cavity and merge-throat information.
    """

    shape = vals.shape

    flat = np.asarray(
        vals,
        dtype=np.float32
    ).ravel()

    eligible = np.flatnonzero(
        flat >= min_radius
    )

    if eligible.size == 0:
        return {
            "cavities": [],
            "merge_throats": [],
            "n_cavities": 0,
            "n_merge_throats": 0,
            "min_radius": float(min_radius),
            "periodic": bool(periodic),
        }

    # Descending distance
    #
    # Equal-valued voxels are processed together as a plateau.
    order_local = np.argsort(
        flat[eligible],
        kind="stable"
    )[::-1]

    order = eligible[order_local]

    sorted_values = flat[order]

    N = flat.size

    # -1 means voxel has not yet entered the superlevel set
    parent = np.full(
        N,
        -1,
        dtype=np.int64
    )

    # Only meaningful for roots that already belong to an
    # older superlevel component.
    peak_for_root = np.full(
        N,
        -1,
        dtype=np.int32
    )

    peak_radius = []
    peak_index = []

    death_radius = []
    elder_cavity = []
    death_index = []

    merge_throats = []

    def find_voxel(x):

        root = x

        while parent[root] != root:
            root = parent[root]

        while parent[x] != root:
            parent[x], x = root, parent[x]

        return root

    def union_plateau(a, b):
        """
        Union used only while building a same-valued plateau.
        """

        ra = find_voxel(a)
        rb = find_voxel(b)

        if ra == rb:
            return ra

        if ra < rb:
            parent[rb] = ra
            return ra

        parent[ra] = rb
        return rb

    # --------------------------------------------------------
    # Process one exact distance level at a time
    # --------------------------------------------------------

    start = 0

    while start < len(order):

        t = sorted_values[start]

        end = start + 1

        while (
            end < len(order)
            and sorted_values[end] == t
        ):
            end += 1

        batch = order[start:end]

        # ----------------------------------------------------
        # Activate whole plateau
        # ----------------------------------------------------

        parent[batch] = batch

        # ----------------------------------------------------
        # Join equal-valued neighboring plateau voxels
        # ----------------------------------------------------

        for idx in batch:

            idx = int(idx)

            for nbr in _neighbors6_flat(
                idx,
                shape,
                periodic=periodic
            ):

                if (
                    parent[nbr] >= 0
                    and flat[nbr] == t
                ):
                    union_plateau(
                        idx,
                        nbr
                    )

        # ----------------------------------------------------
        # Find plateau components and which OLDER components
        # each plateau touches.
        # ----------------------------------------------------

        touched = {}
        representative = {}

        for idx in batch:

            idx = int(idx)

            plateau_root = find_voxel(idx)

            if plateau_root not in touched:
                touched[plateau_root] = set()
                representative[plateau_root] = idx

            for nbr in _neighbors6_flat(
                idx,
                shape,
                periodic=periodic
            ):

                # Strictly older component
                if (
                    parent[nbr] >= 0
                    and flat[nbr] > t
                ):
                    touched[plateau_root].add(
                        find_voxel(nbr)
                    )

        # ----------------------------------------------------
        # Interpret each plateau:
        #
        # 0 older neighbors -> new cavity maximum
        # 1 older component -> ordinary growth
        # >1 older components -> bottleneck / merge
        # ----------------------------------------------------

        for plateau_root in representative:

            # Roots may have changed if an earlier plateau at the
            # same value merged some older components.
            older_roots = {
                find_voxel(r)
                for r in touched[plateau_root]
            }

            # -----------------------------------------------
            # New cavity
            # -----------------------------------------------

            if len(older_roots) == 0:

                cid = len(peak_radius)

                peak_radius.append(
                    float(t)
                )

                peak_index.append(
                    representative[plateau_root]
                )

                death_radius.append(None)
                elder_cavity.append(None)
                death_index.append(None)

                peak_for_root[plateau_root] = cid

                continue

            # -----------------------------------------------
            # Plateau simply grows one component
            # -----------------------------------------------

            if len(older_roots) == 1:

                root = next(
                    iter(older_roots)
                )

                parent[plateau_root] = root

                continue

            # -----------------------------------------------
            # Several components merge:
            # choose the oldest/highest cavity as survivor.
            # -----------------------------------------------

            def survivor_key(root):

                cid = int(
                    peak_for_root[root]
                )

                # Higher maximum survives.
                # Cavity ID provides deterministic tie breaking.
                return (
                    peak_radius[cid],
                    -cid
                )

            survivor = max(
                older_roots,
                key=survivor_key
            )

            survivor_cid = int(
                peak_for_root[survivor]
            )

            throat_idx = representative[
                plateau_root
            ]

            throat_ijk = _flat_to_ijk(
                throat_idx,
                shape
            )

            throat_xyz = _ijk_to_cart(
                throat_ijk,
                meta
            )

            # Every other component dies at this saddle.
            for root in older_roots:

                if root == survivor:
                    continue

                child_cid = int(
                    peak_for_root[root]
                )

                Ri = peak_radius[
                    child_cid
                ]

                Rj = peak_radius[
                    survivor_cid
                ]

                Rt = float(t)

                denom = min(
                    Ri,
                    Rj
                )

                beta = (
                    Rt / denom
                    if denom > 0
                    else np.nan
                )

                persistence = (
                    Ri - Rt
                )

                death_radius[
                    child_cid
                ] = Rt

                elder_cavity[
                    child_cid
                ] = survivor_cid

                death_index[
                    child_cid
                ] = throat_idx

                merge_throats.append({
                    "cavity_i": child_cid,
                    "cavity_j": survivor_cid,

                    "R_i": Ri,
                    "R_j": Rj,

                    "R_throat": Rt,
                    "D_throat": 2.0 * Rt,

                    "beta": beta,

                    # Topological prominence of the child cavity
                    "persistence": persistence,

                    # Representative point on the saddle plateau
                    "grid_index": throat_ijk,
                    "position": throat_xyz,
                })

                # Merge losing component into survivor
                parent[root] = survivor

            # Plateau belongs to surviving component
            parent[plateau_root] = survivor

        start = end

    # --------------------------------------------------------
    # Construct cavity output
    # --------------------------------------------------------

    cavities = []

    for cid, radius in enumerate(
        peak_radius
    ):

        idx = peak_index[cid]

        ijk = _flat_to_ijk(
            idx,
            shape
        )

        xyz = _ijk_to_cart(
            ijk,
            meta
        )

        death = death_radius[cid]

        persistence = (
            radius - death
            if death is not None
            else None
        )

        cavities.append({
            "id": cid,

            "R_cavity": radius,
            "D_cavity": 2.0 * radius,

            "grid_index": ijk,
            "position": xyz,

            "death_radius": death,

            "persistence": persistence,

            # Cavity into which this branch eventually merges
            "elder_cavity": elder_cavity[cid],
        })

    return {
        "cavities": cavities,
        "merge_throats": merge_throats,

        "n_cavities": len(cavities),
        "n_merge_throats": len(merge_throats),

        "min_radius": float(min_radius),
        "periodic": bool(periodic),
    }


# ============================================================
# 8. WIDEST-PATH BOTTLENECK BETWEEN TWO CAVITIES
# ============================================================

def merge_tree_widest_path_radius(
    tree,
    cavity_a,
    cavity_b
):
    """
    Return the bottleneck radius connecting two cavities in the
    merge tree.

    For a path gamma,

        R*(A,B) = max_gamma min_{x in gamma} D(x)

    In the maximum-connectivity merge tree this is simply the
    minimum throat radius along the unique tree path.

    Returns None if the cavities are not connected above the
    tree's min_radius.
    """

    n = tree["n_cavities"]

    if not (
        0 <= cavity_a < n
        and 0 <= cavity_b < n
    ):
        raise ValueError(
            "Invalid cavity ID."
        )

    if cavity_a == cavity_b:
        return tree["cavities"][
            cavity_a
        ]["R_cavity"]

    adj = [
        [] for _ in range(n)
    ]

    for edge in tree["merge_throats"]:

        i = edge["cavity_i"]
        j = edge["cavity_j"]
        Rt = edge["R_throat"]

        adj[i].append(
            (j, Rt)
        )

        adj[j].append(
            (i, Rt)
        )

    stack = [
        (
            cavity_a,
            np.inf,
            -1
        )
    ]

    while stack:

        node, bottleneck, parent_node = stack.pop()

        if node == cavity_b:
            return float(
                bottleneck
            )

        for nbr, Rt in adj[node]:

            if nbr == parent_node:
                continue

            stack.append(
                (
                    nbr,
                    min(
                        bottleneck,
                        Rt
                    ),
                    node
                )
            )

    return None