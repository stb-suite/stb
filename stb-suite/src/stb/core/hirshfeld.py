"""Hirshfeld-I (iterative) atomic charge partitioning -- the genuinely new
physics behind stb-hirshfeldIons/Analysis (Stages 2/3). Not needed by/not
imported by anything reading SIESTA's own NATIVE Hirshfeld output
(stb-nativecharges, core/siesta_log.py) -- that's a plain text-table
reader, this module is the actual density-partitioning method SIESTA
itself doesn't implement iteratively.

Shared by 2 consumers from the start (hirshfeld_ions.py's informational
pass-0 report, hirshfeld_analysis.py's full iteration), so it lives
directly in core/ rather than waiting for an extract-on-second-use move.

All grid arrays here follow the same convention as core/rho_io.py /
core/grid_export.py: a plain numpy array from sisl's Grid.grid (already
totaled via core.rho_io.read_total, never a raw per-spin-channel read),
paired with the geometry's own 3x3 lattice (Angstrom). Distances use the
same minimum-image (periodic-aware) fractional-coordinate convention as
aimd_analysis.py's compute_rdf: `disp -= np.round(disp)` before converting
to Cartesian, so an atom near a cell edge still gets density contributions
wrapping around from its own nearest periodic image -- both for an
isolated atom's own vacuum box (where the atom is always placed at the
exact fractional center by structure_io.write_isolated_atom_fdf, so this
reduces to a plain, non-wrapped distance) and for the combined production
system (where atoms can sit anywhere, including near a cell face, and the
wrapping genuinely matters).
"""

from __future__ import annotations

import time

import numpy as np

from stb.core.cli import print_progress_line, finish_progress_line


def radial_density_profile(grid_data: np.ndarray, lattice: np.ndarray, center_frac,
                            r_max: float | None = None, n_bins: int = 200):
    """Spherically-averages `grid_data` (a real-space density grid, e.g. an
    isolated atom's own total .RHO) around `center_frac` (fractional
    coordinates of the atom, e.g. (0.5, 0.5, 0.5) for
    write_isolated_atom_fdf's box-centered convention), using minimum
    -image distance. Returns (r_bin_centers, rho_of_r), each a (n_bins,)
    array -- rho_of_r is the mean grid value among voxels whose
    minimum-image distance to `center_frac` falls in that bin (0.0 for a
    bin no voxel reaches, e.g. one near r_max in a cubic box's corner
    region no sphere of that radius fully samples).

    `r_max` defaults to the largest distance actually present on the grid
    (so the profile always covers every voxel) -- for an isolated atom in
    a large vacuum box this comfortably exceeds where the density has
    already decayed to ~0, which is exactly what
    interpolate_reference_density's zero-extrapolation beyond r_max
    relies on.
    """
    nx, ny, nz = grid_data.shape
    fx, fy, fz = np.meshgrid(np.arange(nx) / nx, np.arange(ny) / ny, np.arange(nz) / nz,
                              indexing='ij')
    disp = np.stack([fx, fy, fz], axis=-1) - np.asarray(center_frac, dtype=float)
    disp -= np.round(disp)
    cart = disp @ lattice
    r = np.linalg.norm(cart, axis=-1)

    if r_max is None:
        r_max = float(r.max())
    bin_edges = np.linspace(0.0, r_max, n_bins + 1)
    bin_idx = np.clip(np.digitize(r.ravel(), bin_edges) - 1, 0, n_bins - 1)
    sums = np.bincount(bin_idx, weights=grid_data.ravel(), minlength=n_bins)
    counts = np.bincount(bin_idx, minlength=n_bins)
    rho_of_r = np.divide(sums, counts, out=np.zeros(n_bins), where=counts > 0)
    r_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    # Prepend r=0 and append r=r_max as explicit interpolation anchors,
    # reusing the first/last bin's own average. Without this, EVERY query
    # exactly at r=0 (the atom's own position -- the single most important
    # point of the whole profile, its density peak) or r=r_max falls
    # OUTSIDE [r_centers.min(), r_centers.max()] (bin centers sit strictly
    # inside (0, r_max) by construction), which silently zeroed out the
    # nucleus's own density under callers' left=0.0/right=0.0 np.interp
    # convention -- a real bug caught by this module's own analytic unit
    # test (a single isolated atom's reconstructed population came out
    # measurably short of its true integral; traced to exactly this gap).
    r_out = np.concatenate(([0.0], r_centers, [r_max]))
    rho_out = np.concatenate(([rho_of_r[0]], rho_of_r, [rho_of_r[-1]]))
    return r_out, rho_out


def interpolate_reference_density(target_r, r_neutral, rho_neutral, r_ion, rho_ion,
                                   charge: float, cap: float = 1.0):
    """rho_ref(target_r) = (1-frac)*rho_neutral(target_r) + frac*rho_ion(target_r),
    frac = min(|charge|, cap) -- the Hirshfeld-I reference-density
    interpolation. Both profiles are linearly interpolated (np.interp)
    onto `target_r` independently (they need not share the same r-grid),
    extrapolating to 0.0 beyond either profile's own r_max -- matches the
    physical expectation of the density vanishing far from the nucleus,
    and keeps a slightly different r_max between a species' neutral/ion
    calc from corrupting the blend.
    """
    frac = min(abs(charge), cap)
    neutral_i = np.interp(target_r, r_neutral, rho_neutral, left=0.0, right=0.0)
    ion_i = np.interp(target_r, r_ion, rho_ion, left=0.0, right=0.0)
    return (1.0 - frac) * neutral_i + frac * ion_i


def compute_hirshfeld_charges(rho_real: np.ndarray, lattice: np.ndarray, frac_positions,
                               z_vals, species_profiles: dict, progress_label: str | None = None):
    """One Hirshfeld pass against a FIXED real density `rho_real` (the
    combined/production system's own total .RHO). `species_profiles` maps
    atom index -> (r_array, rho_of_r_array): the CURRENT iteration's
    reference profile for that atom (the neutral-only profile for pass 0
    / plain Hirshfeld, or an interpolate_reference_density blend for a
    later Hirshfeld-I round).

    Promolecule density at each grid voxel = sum of every atom's own
    reference profile evaluated at its minimum-image distance to that
    voxel (an atom near a cell edge still gets contributions "wrapping
    around" from its own nearest periodic image, as it would in the real
    crystal). weight_A(voxel) = rho_A(voxel) / promolecule(voxel) (0 where
    promolecule is 0, i.e. no atom's reference reaches that far -- a true
    vacuum region contributes nothing to any atom's population either
    way). charge_A = Z_A - integral(weight_A * rho_real).

    Two passes over the atom list (accumulate the shared promolecule sum,
    then integrate each atom's own weighted share) keep peak memory at
    O(grid size) instead of O(atoms x grid size) -- each atom's own
    density array is recomputed rather than cached, trading compute for
    memory, deliberately, since a production grid can be large. This is
    the slow, silent part of the whole workflow for a production-sized
    system (many atoms x a fine grid) -- `progress_label`, if given (e.g.
    "Pass 0" or "Iteration 3/20"), prints a live step/elapsed-time line via
    core.cli.print_progress_line over both passes combined (2*n_atoms
    steps) so a long run doesn't look frozen; None (the default) stays
    silent, for callers that don't want stray stderr output (e.g. a unit
    test).

    The per-atom minimum-image GEOMETRY arrays (grid_frac/disp/cart/r --
    the ones whose size is O(grid size), not O(grid size) x 8 bytes worth
    of double precision anyone actually needs for an Angstrom-scale
    distance) are float32, not float64 -- half the memory of the original
    implementation for a real production grid (SIESTA's own .RHO is
    already float32; the density/population arithmetic itself stays
    float64 throughout, only the distance computation is narrowed).
    `del`eted explicitly as soon as each is no longer needed rather than
    relying on CPython's refcounting alone to free it -- this function's
    own peak was observed to OOM-kill a real 216x216x180 (~8.4M point)
    production run otherwise.

    Returns (charges, populations), each a (n_atoms,) array.
    """
    n_atoms = len(frac_positions)
    nx, ny, nz = rho_real.shape
    fx, fy, fz = np.meshgrid(np.arange(nx, dtype=np.float32) / nx,
                              np.arange(ny, dtype=np.float32) / ny,
                              np.arange(nz, dtype=np.float32) / nz,
                              indexing='ij')
    grid_frac = np.stack([fx, fy, fz], axis=-1)
    del fx, fy, fz
    lattice32 = np.asarray(lattice, dtype=np.float32)
    cell_volume = abs(np.linalg.det(lattice))
    voxel_volume = cell_volume / rho_real.size

    def atom_density(i):
        disp = grid_frac - np.asarray(frac_positions[i], dtype=np.float32)
        disp -= np.round(disp)
        cart = disp @ lattice32
        del disp
        r = np.linalg.norm(cart, axis=-1)
        del cart
        r_arr, rho_arr = species_profiles[i]
        return np.interp(r, r_arr, rho_arr, left=0.0, right=0.0)

    start_time = time.monotonic()
    total_steps = 2 * n_atoms

    def report(step, atom_num, phase):
        if progress_label is None:
            return
        elapsed = time.monotonic() - start_time
        print_progress_line(
            f"  {progress_label}: {phase} atom {atom_num}/{n_atoms}, {elapsed:6.1f}s elapsed",
            step, total_steps)

    promolecule = np.zeros(rho_real.shape, dtype=float)
    for i in range(n_atoms):
        promolecule += atom_density(i)
        report(i + 1, i + 1, "building density,")
    safe_promolecule = np.where(promolecule > 0, promolecule, 1.0)

    charges = np.zeros(n_atoms)
    populations = np.zeros(n_atoms)
    for i in range(n_atoms):
        rho_a = atom_density(i)
        weight = np.where(promolecule > 0, rho_a / safe_promolecule, 0.0)
        pop = float(np.sum(weight * rho_real) * voxel_volume)
        populations[i] = pop
        charges[i] = z_vals[i] - pop
        report(n_atoms + i + 1, i + 1, "integrating")

    if progress_label is not None:
        finish_progress_line()

    return charges, populations


def iterate_hirshfeld_i(rho_real: np.ndarray, lattice: np.ndarray, frac_positions, symbols,
                         z_vals, neutral_profiles: dict, cation_profiles: dict,
                         anion_profiles: dict, tol: float = 0.005, max_iter: int = 20,
                         show_progress: bool = False, initial_charges=None,
                         initial_populations=None, on_iteration=None):
    """The Hirshfeld-I convergence loop, following the original iterative-
    Hirshfeld formulation (Bultinck et al., J. Chem. Phys. 126, 144111
    (2007)) literally: EVERY ATOM independently interpolates its own
    reference density between its species' NEUTRAL profile and whichever
    of that species' CATION or ANION profile matches that atom's OWN
    charge sign from the previous round -- so two atoms of the same
    species can legitimately end up leaning toward opposite ion states in
    the same round (e.g. one edge carbon reading slightly positive, the
    bulk of the carbons reading negative). This is why both
    `cation_profiles` and `anion_profiles` (species symbol -> (r_array,
    rho_of_r_array), from radial_density_profile on that species'
    isolated +1/-1 '*.RHO') must always be supplied, never just "the"
    ion for a species -- a coarser variant that picks one sign per
    species from a single vote (this codebase's own earlier
    implementation) cannot represent this and silently mis-references
    every atom on the "wrong" side of that vote.

    Iteration 0 is plain simple Hirshfeld (every atom's reference = its
    own species' neutral profile, frac=0 -- directly comparable to
    stb-nativecharges'/SIESTA's own native Hirshfeld output for the same
    system). Each subsequent round recomputes charges against the fixed
    real density using this per-atom cation-or-anion blend
    (interpolate_reference_density, frac = min(|q|, 1)), until
    max|Delta q| < tol or max_iter rounds are used.

    `show_progress`, if True, forwards a per-round progress_label
    ("Pass 0"/"Iteration N/max_iter") to each compute_hirshfeld_charges
    call -- this loop repeats that function's own already-slow atom scan
    up to max_iter+1 times, so a production run benefits from the live
    step/elapsed-time line even more than a single compute_hirshfeld_charges
    call does.

    `initial_charges`/`initial_populations`, if both given, are used
    directly as iteration 0's result instead of recomputing it -- a caller
    that already ran its own pass-0 compute_hirshfeld_charges call (e.g.
    to print a "simple Hirshfeld" comparison table before iterating) would
    otherwise pay for that identical, expensive full-grid atom scan twice.

    `on_iteration`, if given, is called after EVERY round >= 1 (never for
    pass 0, which has no delta to report) as
    on_iteration(iteration, max_delta, worst_atom_idx, charges, populations)
    -- `worst_atom_idx` is the 0-based index of the single atom that moved
    the most this round (i.e. argmax|new_charges - old_charges|, the same
    atom max_delta itself came from), `charges`/`populations` are this
    round's full, already-updated arrays. Lets a caller print a live,
    per-round convergence line (which atom/species is driving the
    tolerance check, not just the bare number) instead of only a summary
    once the whole loop finishes -- this loop's own inner
    compute_hirshfeld_charges calls are already the slow part (see
    `show_progress`), so surfacing per-round detail as it happens, not
    only at the end, matters for a long production run.

    Returns (charges, populations, history), where history is a list of
    (iteration, max_abs_delta_q) tuples, iteration 0's own entry carrying
    max_abs_delta_q=None (no previous round to compare against).
    """
    n_atoms = len(frac_positions)
    if initial_charges is not None and initial_populations is not None:
        charges, populations = initial_charges, initial_populations
    else:
        species_profiles = {i: neutral_profiles[symbols[i]] for i in range(n_atoms)}
        charges, populations = compute_hirshfeld_charges(
            rho_real, lattice, frac_positions, z_vals, species_profiles,
            progress_label="Pass 0" if show_progress else None)
    history = [(0, None)]

    for iteration in range(1, max_iter + 1):
        species_profiles = {}
        for i in range(n_atoms):
            sym = symbols[i]
            r_n, rho_n = neutral_profiles[sym]
            # Each atom picks its OWN ion state from its OWN charge sign
            # this round -- not a per-species decision. A charge of
            # exactly 0.0 (vanishingly rare in practice) leans cation
            # (frac=0 either way makes the choice moot: rho_ref = rho_n).
            r_i, rho_i = cation_profiles[sym] if charges[i] >= 0 else anion_profiles[sym]
            blended = interpolate_reference_density(r_n, r_n, rho_n, r_i, rho_i, charges[i])
            species_profiles[i] = (r_n, blended)

        new_charges, new_populations = compute_hirshfeld_charges(
            rho_real, lattice, frac_positions, z_vals, species_profiles,
            progress_label=f"Iteration {iteration}/{max_iter}" if show_progress else None)
        delta = np.abs(new_charges - charges)
        worst_atom_idx = int(np.argmax(delta))
        max_delta = float(delta[worst_atom_idx])
        history.append((iteration, max_delta))
        charges, populations = new_charges, new_populations
        if on_iteration is not None:
            on_iteration(iteration, max_delta, worst_atom_idx, charges, populations)
        if max_delta < tol:
            break

    return charges, populations, history
