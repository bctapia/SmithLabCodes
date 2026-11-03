import numpy as np
import MDAnalysis as mda
import freud


def c2_bondvec(
    data_file,
    dump_file,
    timestep_fs=1.0,
    max_lag_ps=20000.0,
    stride=1,
    backbone_selection=None,
    origin_stride_ps=1.0,
    output_file="c2_bondvec.csv",
):
    """
    C2(t) = <P2(u(t0).u(t0+t))>
    where u is the unit bond vector. Average over bonds and time origins.

    Args:
        data_file: LAMMPS data file
        dump_file: LAMMPS dump file
            Trajectory MUST be unwrapped. If trajectory is wrapped but image flags are present, ensure MDAnalysis is properly unwrapping the coordinates
        timestep_fs: Timestep MD integration
        max_lag_ps: The maximum time separation of each calculation
        stride: How many timestep_fs to stride between collecting u(t0+t) (=1 means analyze at every timestep)
        origin_stride_ps: How many timestep_fs to stride before starting another u(t0) ensemble average (=1 means average using every timestep)
        backbone_selecton: Analyze the C2 of only specific bond-vectors using MDAnalysis "sel" syntax
        output_file: Path of output file
    """

    def p2(x):
        return 0.5 * (3.0 * x * x - 1.0)

    print("Loading universe...")
    u = mda.Universe(data_file, dump_file, topology_format="DATA", format="LAMMPSDUMP")

    if len(u.bonds) == 0:
        raise RuntimeError("No bonds found in topology. Ensure DATA file has a Bonds section.")

    print("===Topology Info===")
    print("Atoms:", u.atoms.n_atoms)
    print("Frames:", u.trajectory.n_frames)
    print("Bonds:", len(u.bonds))
    print("===================")

    u.trajectory[0]
    t0 = u.trajectory.time
    u.trajectory[1]
    t1 = u.trajectory.time
    dt_step = t1 - t0   # in LAMMPS steps
    dt_ps = dt_step * timestep_fs / 1000.0

    # Bond indices in terms of atom indices into u.atoms.positions
    ij = u.bonds.indices
    i_all = ij[:, 0].astype(np.int64)
    j_all = ij[:, 1].astype(np.int64)

    # Optional: filter bonds so both atoms are in a selection
    if backbone_selection is not None:
        sel = u.select_atoms(backbone_selection)
        sel_idx = sel.indices
        mask = np.isin(i_all, sel_idx) & np.isin(j_all, sel_idx)
        i = i_all[mask]
        j = j_all[mask]
        print(f"Filtered bonds with selection [{backbone_selection}]. {len(i)} bonds remain.")
    else:
        i, j = i_all, j_all
        print(f"Using all bonds: {len(i)} bonds")

    nb = len(i)
    if nb == 0:
        raise RuntimeError("No bonds available after selection.")

    dt_sample_ps = float(dt_ps) * int(stride)
    
    # requested lag in frames
    max_lag_req = int(round(float(max_lag_ps) / dt_sample_ps))
    
    # maximum possible lag given the file length
    n_total = u.trajectory.n_frames
    n_proc = (n_total - 1) // stride + 1   # frames actually visited by [::stride]
    max_lag_possible = u.trajectory.n_frames - 1
    max_lag_possible = n_proc - 1
    max_lag = min(max_lag_req, max_lag_possible)
    

    origin_stride = max(1, int(round(float(origin_stride_ps) / dt_sample_ps)))

    print(f"Sampled dt = {dt_sample_ps:.6g} ps")
    print(f"MAX_LAG_PS = {max_lag_ps} ps -> max_lag = {max_lag} frames")
    print(f"ORIGIN_STRIDE_PS = {origin_stride_ps} ps -> origin_stride = {origin_stride} frames")

    # Accumulators per lag
    c2_sum = np.zeros(max_lag + 1, dtype=np.float64)
    c2_cnt = np.zeros(max_lag + 1, dtype=np.int64)

    # Ring buffer of past unit bond vectors: each entry shape (nb, 3)
    ring = [None] * (max_lag + 1)

    atoms = u.atoms
    ring_step = -1

    print("Streaming trajectory and computing C2...")
    for ts in u.trajectory[::stride]:
        ring_step += 1

        pos = atoms.positions
        vec = pos[j] - pos[i]                       # (nb, 3)
        norm = np.linalg.norm(vec, axis=1)          # (nb,)
        if np.any(norm == 0.0):
            raise RuntimeError(
                f"Zero-length bond vector encountered at frame {ring_step} (t={ts.time})."
            )
        uvec = (vec / norm[:, None]).astype(np.float32)

        # Store current frame in ring
        ring[ring_step % (max_lag + 1)] = uvec

        # For each lag, origin frame = ring_step - lag
        lag_max_here = min(ring_step, max_lag)

        # lag = 0 always valid; origin == current
        # Only count if the origin (ring_step) is on the origin grid
        if (ring_step % origin_stride) == 0:
            dots0 = np.einsum("ij,ij->i", uvec, uvec)
            c2_sum[0] += p2(dots0).mean()
            c2_cnt[0] += 1

        for lag in range(1, lag_max_here + 1):
            origin_idx = ring_step - lag
            if (origin_idx % origin_stride) != 0:
                continue  # this origin is not selected

            u0 = ring[origin_idx % (max_lag + 1)]
            # dot per bond, vectorized over bonds
            dots = np.einsum("ij,ij->i", uvec, u0)   # (nb,)
            c2_sum[lag] += p2(dots).mean()           # average over bonds
            c2_cnt[lag] += 1

        if (ring_step % 200) == 0 and ring_step > 0:
            print(f"Processed sampled frame {ring_step} (t ~ {ring_step * dt_sample_ps:.2f} ps)")

    valid = c2_cnt > 0
    c2 = np.zeros_like(c2_sum)
    c2[valid] = c2_sum[valid] / c2_cnt[valid]

    times_ps = np.arange(max_lag + 1, dtype=np.float64) * dt_sample_ps
    out = np.column_stack((times_ps, c2, c2_cnt))

    np.savetxt(
        output_file,
        out,
        delimiter=",",
        header="time_ps,C2,navg_origins",
        comments="",
    )
    print(f"Done. Wrote {output_file}")

def sq(
    data_file,
    dump_file,
    selection=None,
    q_min_Ainv=0.4,
    q_max_Ainv=3.0,
    bins=200,
    stride=1,
    output_file="sq.csv",
    num_proc=None,
    method="direct"
):
    """
    Compute isotropically-averaged static structure factor S(q) using freud.

    Uses: freud.diffraction.StaticStructureFactorDirect(bins, k_min, k_max) :contentReference[oaicite:2]{index=2}

    Args:
        data_file, dump_file: LAMMPS DATA + dump
        selection: MDAnalysis selection (use same group as Fs)
        q_min_Ainv, q_max_Ainv: q-range (1/Å)
        bins: number of q bins
        stride: analyze every `stride` trajectory frames
        output_file: CSV with columns q_Ainv,Sq,nframes
    """
    if q_max_Ainv <= q_min_Ainv:
        raise ValueError("q_max_Ainv must be > q_min_Ainv")
    if bins < 10:
        raise ValueError("bins should be >= 10")

    print("Loading universe...")
    u = mda.Universe(data_file, dump_file, topology_format="DATA", format="LAMMPSDUMP")
    
    if selection is None:
        ag = u.atoms
        print(f"Using all atoms: {ag.n_atoms}")
    else:
        ag = u.select_atoms(selection)
        print(f"Using selection [{selection}]: {ag.n_atoms} atoms")
    
    N = ag.n_atoms
    if N == 0:
        raise RuntimeError(f"Selection [{selection}] returned 0 atoms")

    print(f"Selected atoms: {N}")
    print(f"Frames: {u.trajectory.n_frames} (stride={stride})")
    freud.parallel.set_num_threads(nthreads=num_proc)
    if method=="direct":
        sf = freud.diffraction.StaticStructureFactorDirect(
            bins=bins, k_min=q_min_Ainv, k_max=q_max_Ainv
        )
    elif method=="debye":
        sf = freud.diffraction.StaticStructureFactorDebye(
            num_k_values=bins, k_min=q_min_Ainv, k_max=q_max_Ainv
        )

    Sq_sum = np.zeros(bins, dtype=np.float64)
    nframes = 0

    for ts in u.trajectory[::stride]:
        # MDAnalysis provides unit cell vectors as a (3,3) matrix here :contentReference[oaicite:3]{index=3}
        box_mda = ts.triclinic_dimensions.astype(np.float64)  # shape (3,3)
        box = freud.box.Box.from_matrix(box_mda)

        # Wrap positions into the primary cell (freud expects positions inside the box)
        r = ag.positions.astype(np.float64)  # Å
        r = box.wrap(r)

        sf.compute((box, r))  # :contentReference[oaicite:4]{index=4}

        # freud stores k values + S(k) after compute
        Sq_sum += sf.S_k
        nframes += 1

        if (nframes % 50) == 0:
            print(f"Processed {nframes} frames")

    if nframes == 0:
        raise RuntimeError("No frames processed (check stride / trajectory).")

    if method == "debye":
        q = np.asarray(sf.k_values, dtype=np.float64).copy()
    else:
        q = np.asarray(sf.bin_centers, dtype=np.float64).copy()

    Sq = Sq_sum / nframes

    # crude "first peak" finder: first local maximum above q>0
    # (skip very small q where finite-size effects dominate)
    q_star = np.nan
    if len(q) >= 3:
        start = np.searchsorted(q, max(q_min_Ainv, 0.2))  # ignore q < 0.2 1/Å by default
        for i in range(max(start, 1), len(q) - 1):
            if Sq[i] > Sq[i - 1] and Sq[i] > Sq[i + 1]:
                q_star = q[i]
                break

    out = np.column_stack([q, Sq, np.full_like(q, nframes, dtype=np.int64)])
    np.savetxt(
        output_file,
        out,
        delimiter=",",
        header="q_Ainv,Sq,nframes_avg",
        comments="",
    )
    print(f"Wrote {output_file}")
    if np.isfinite(q_star):
        print(f"Estimated first-peak q* ≈ {q_star:.4f} 1/Å")
    else:
        print("Did not find a clear first local maximum (try larger q_max_Ainv or more bins).")

    return q, Sq, q_star


def fs_self(
    data_file,
    dump_file,
    q_Ainv=1.0,
    timestep_fs=1.0,
    max_lag_ps=20000.0,
    stride=1,
    selection=None,
    origin_stride_ps=1.0,
    output_file="fs_self.csv",
):
    """
    Self-intermediate scattering function, isotropic/orientation-averaged:
        Fs(q,t) = < sin(q|dr|) / (q|dr|) >
    where dr = r(t0+t) - r(t0). Average over atoms in `selection` and time origins.

    Args:
        data_file: LAMMPS data file (for topology)
        dump_file: LAMMPS dump file
            Trajectory should be unwrapped for meaningful displacements.
        q_Ainv: |q| in 1/Angstrom.
        timestep_fs: Timestep between *dumped* frames in fs (not necessarily MD dt).
        max_lag_ps: Maximum lag time to compute, in ps.
        stride: Analyze every `stride` frames of the dump.
        selection: MDAnalysis selection for atoms to include (e.g. "type 1", "name C*", "all").
        origin_stride_ps: Spacing between time origins, in ps.
        output_file: CSV output path.
    """
    if q_Ainv <= 0:
        raise ValueError("q_Ainv must be > 0.")

    print("Loading universe...")
    u = mda.Universe(data_file, dump_file, topology_format="DATA", format="LAMMPSDUMP")

    print("===Trajectory Info===")
    print("Atoms:", u.atoms.n_atoms)
    print("Frames:", u.trajectory.n_frames)
    print("=====================")

    # Determine dt between stored frames in ps (robust to LAMMPS 'time' being in steps)
    if u.trajectory.n_frames < 2:
        raise RuntimeError("Need at least 2 frames in the trajectory.")

    u.trajectory[0]
    t0 = u.trajectory.time
    u.trajectory[1]
    t1 = u.trajectory.time
    dt_step = t1 - t0                  # whatever units LAMMPS wrote into 'time'
    dt_ps = dt_step * timestep_fs / 1000.0

    dt_sample_ps = float(dt_ps) * int(stride)

    # Requested lag (in sampled frames)
    max_lag_req = int(round(float(max_lag_ps) / dt_sample_ps))

    n_total = u.trajectory.n_frames
    n_proc = (n_total - 1) // stride + 1
    max_lag_possible = n_proc - 1
    max_lag = min(max_lag_req, max_lag_possible)

    origin_stride = max(1, int(round(float(origin_stride_ps) / dt_sample_ps)))
    
    if selection is None:
        atoms = u.atoms
        print(f"Using all atoms: {atoms.n_atoms}")
    else:
        atoms = u.select_atoms(selection)
        print(f"Using selection [{selection}]: {atoms.n_atoms} atoms")

    n_sel = atoms.n_atoms
    if n_sel == 0:
        raise RuntimeError(f"Selection [{selection}] returned 0 atoms.")

    print(f"Selected atoms: {n_sel}")
    print(f"Sampled dt = {dt_sample_ps:.6g} ps")
    print(f"MAX_LAG_PS = {max_lag_ps} ps -> max_lag = {max_lag} frames")
    print(f"ORIGIN_STRIDE_PS = {origin_stride_ps} ps -> origin_stride = {origin_stride} frames")
    print(f"q = {q_Ainv} 1/Å")

    # Accumulators per lag
    fs_sum = np.zeros(max_lag + 1, dtype=np.float64)
    fs_cnt = np.zeros(max_lag + 1, dtype=np.int64)

    # Ring buffer of positions for selected atoms: each entry shape (n_sel, 3)
    ring = [None] * (max_lag + 1)
    ring_step = -1

    def sinc_real(x):
        # returns sin(x)/x with correct x->0 limit
        out = np.ones_like(x, dtype=np.float64)
        mask = x != 0.0
        out[mask] = np.sin(x[mask]) / x[mask]
        return out

    print("Streaming trajectory and computing Fs(q,t)...")
    for ts in u.trajectory[::stride]:
        ring_step += 1

        pos = atoms.positions.astype(np.float32)  # Å
        ring[ring_step % (max_lag + 1)] = pos

        lag_max_here = min(ring_step, max_lag)

        # lag=0: dr=0 => Fs=1, but keep the same origin-grid logic as your C2
        if (ring_step % origin_stride) == 0:
            fs_sum[0] += 1.0
            fs_cnt[0] += 1

        for lag in range(1, lag_max_here + 1):
            origin_idx = ring_step - lag
            if (origin_idx % origin_stride) != 0:
                continue

            r0 = ring[origin_idx % (max_lag + 1)]
            dr = pos - r0                              # (n_sel, 3)
            r = np.linalg.norm(dr, axis=1).astype(np.float64)  # (n_sel,)

            x = q_Ainv * r
            fs_val = sinc_real(x).mean()               # average over atoms
            fs_sum[lag] += fs_val
            fs_cnt[lag] += 1

        if (ring_step % 200) == 0 and ring_step > 0:
            print(f"Processed sampled frame {ring_step} (t ~ {ring_step * dt_sample_ps:.2f} ps)")

    valid = fs_cnt > 0
    fs = np.zeros_like(fs_sum)
    fs[valid] = fs_sum[valid] / fs_cnt[valid]

    times_ps = np.arange(max_lag + 1, dtype=np.float64) * dt_sample_ps
    out = np.column_stack((times_ps, fs, fs_cnt))

    np.savetxt(
        output_file,
        out,
        delimiter=",",
        header="time_ps,Fs_q,navg_origins",
        comments="",
    )
    print(f"Done. Wrote {output_file}")


def msd(
    data_file,
    dump_file,
    timestep_fs=1.0,
    max_lag_ps=20000.0,
    origin_stride_ps=1.0,
    atom_selection=None,
    output_file="msd.csv",
    remove_drift=True,
):
    """
    Calculate whole-chain center-of-mass MSD within an NVT trajectory.

    MSD(tau) = <|R_i(t0 + tau) - R_i(t0)|^2>

    Averages equally over chains and selected time origins.

    Requirements:
        - One LAMMPS molecule ID per complete chain.
        - Dump contains continuously unwrapped xu yu zu coordinates.
        - Fixed simulation box and uniformly spaced trajectory frames.
        - Coordinates in Angstrom.

    Parameters:
        timestep_fs: MD integration timestep, NOT the dump interval.
        atom_selection: Selects chains to include. All atoms belonging
                        to each selected chain are included in its COM.
                        None includes every chain.
        remove_drift: Remove mass-weighted collective translation.

    Returns:
        times_ps, msd_A2, n_origins
    """

    if not np.isfinite(timestep_fs) or timestep_fs <= 0:
        raise ValueError("timestep_fs must be finite and positive.")
    if not np.isfinite(max_lag_ps) or max_lag_ps < 0:
        raise ValueError("max_lag_ps must be finite and nonnegative.")
    if not np.isfinite(origin_stride_ps) or origin_stride_ps <= 0:
        raise ValueError("origin_stride_ps must be finite and positive.")

    u = mda.Universe(
        data_file,
        dump_file,
        topology_format="DATA",
        format="LAMMPSDUMP",
        atom_style="id resid type charge x y z",
        lammps_coordinate_convention="unwrapped",
        dt=timestep_fs / 1000.0,
    )

    try:
        n_frames = len(u.trajectory)

        if n_frames < 2:
            raise ValueError("Trajectory must contain at least two frames.")

        if atom_selection is None:
            selected = u.atoms
        else:
            selected = u.select_atoms(atom_selection)

        if selected.n_atoms == 0:
            raise ValueError("The atom selection returned no atoms.")

        # Expand the selection to complete chains.
        chain_atoms = selected.residues.atoms
        chains = chain_atoms.split("residue")
        n_chains = len(chains)

        if not np.all(
            np.isfinite(chain_atoms.masses) & (chain_atoms.masses > 0)
        ):
            raise ValueError("All chain atoms must have valid positive masses.")

        chain_masses = np.array(
            [chain.masses.sum() for chain in chains],
            dtype=np.float64,
        )

        com = np.empty((n_frames, n_chains, 3), dtype=np.float64)
        steps = np.empty(n_frames, dtype=np.int64)

        u.trajectory[0]
        reference_box = u.trajectory.ts.triclinic_dimensions.copy()

        print(f"Chains: {n_chains}")
        print(f"Frames: {n_frames}")
        print("Computing unwrapped chain COM positions...")

        for iframe, ts in enumerate(u.trajectory):
            steps[iframe] = int(ts.data["step"])

            if not np.allclose(
                ts.triclinic_dimensions,
                reference_box,
                rtol=1e-6,
                atol=1e-5,
            ):
                raise ValueError(
                    "The box changes during this trajectory. "
                    "This function assumes a fixed NVT box."
                )

            for ichain, chain in enumerate(chains):
                com[iframe, ichain] = chain.center_of_mass()

            if iframe > 0 and iframe % 500 == 0:
                print(f"Processed {iframe}/{n_frames} frames")

    finally:
        u.trajectory.close()

    if not np.all(np.isfinite(com)):
        raise ValueError("Nonfinite chain COM coordinates encountered.")

    # Determine frame spacing directly from LAMMPS timestep numbers.
    step_intervals = np.diff(steps)

    if np.any(step_intervals <= 0):
        raise ValueError("LAMMPS timestep numbers must increase.")
    if not np.all(step_intervals == step_intervals[0]):
        raise ValueError("Trajectory frames must be uniformly spaced.")

    dt_ps = float(step_intervals[0]) * timestep_fs / 1000.0

    if remove_drift:
        system_com = np.average(com, axis=1, weights=chain_masses)
        com -= system_com[:, None, :]

    max_lag = min(
        int(np.floor(max_lag_ps / dt_ps)),
        n_frames - 1,
    )
    origin_stride = max(
        1,
        int(round(origin_stride_ps / dt_ps)),
    )

    print(f"Frame spacing: {dt_ps:.6g} ps")
    print(f"Maximum lag: {max_lag * dt_ps:.6g} ps")
    print(f"Actual origin spacing: {origin_stride * dt_ps:.6g} ps")
    print(f"System drift removed: {remove_drift}")

    msd_vals = np.zeros(max_lag + 1, dtype=np.float64)
    n_origins = np.zeros(max_lag + 1, dtype=np.int64)

    for lag in range(max_lag + 1):
        origins = np.arange(0, n_frames - lag, origin_stride)

        dr = com[origins + lag] - com[origins]
        squared_displacements = np.sum(dr**2, axis=2)

        msd_vals[lag] = squared_displacements.mean()
        n_origins[lag] = len(origins)

    times_ps = np.arange(max_lag + 1) * dt_ps

    np.savetxt(
        output_file,
        np.column_stack((times_ps, msd_vals, n_origins)),
        delimiter=",",
        header="time_ps,chain_COM_MSD_A2,n_origins",
        comments="",
        fmt=["%.10g", "%.10g", "%d"],
    )

    print(f"Wrote {output_file}")

    return times_ps, msd_vals, n_origins