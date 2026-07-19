#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "2.1.0"

import os
import sys
import argparse
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt
from ase import Atoms
from ase.io import write as ase_write
from pymatgen.io.ase import AseAtomsAdaptor

from stb.core import structure_io, kspace, mace_relax, mace_phonons
from stb.core.cli import color_text, show_intro, print_dual, print_section
from stb.core.deps import require_mace
from stb.phonons_pos import build_band_path, band_path_to_phonopy_format, band_tick_positions, pretty_label

require_mace()

REPORT_FILE = "stb_mlphonons_report.txt"

_SERIES_COLORS = ['tab:blue', 'tab:orange']

# Same 3-format convention already used by stb-ani2traj/stb-mlmd for
# multi-frame trajectories (xsf/pdb natively carry a lattice, extxyz via a
# Lattice= header tag) -- duplicated here rather than imported, same as
# mlmd.py already does, since it's tiny and tied to each tool's own CLI
# --animate-format choices/help text.
OUTPUT_FORMATS = {
    "xsf": ("xsf", ".xsf"),
    "pdb": ("proteindatabank", ".pdb"),
    "xyz": ("extxyz", ".xyz"),
}

# Below this (THz), a "negative" frequency is treated as numerical noise
# (e.g. residual acoustic-sum-rule violation after symmetrize_force_
# constants(), or plain floating-point roundoff at Gamma where the true
# acoustic-branch frequency is exactly 0) rather than a genuine imaginary/
# unstable mode -- verified live: post-ASR-correction, a real run's worst
# Gamma-point residual was ~5e-5 THz, several orders of magnitude below
# this threshold; a genuinely unstable mode is typically at least ~0.1 THz.
_IMAGINARY_MODE_TOL_THZ = -0.01


def print_symmetry_table(phonon, f_out=None):
    """Reports the symmetry reduction Phonopy applied (via spglib) when
    generating the displacement dataset -- how many MACE evaluations were
    actually needed vs. the naive (no symmetry) count. Adapted from the
    now-deleted phonons_ml.py (consolidated into this tool)."""
    sym = phonon.symmetry
    space_group = sym.get_international_table() or "unknown"
    point_group = sym.pointgroup_symbol
    n_ops = len(sym.symmetry_operations['rotations'])
    n_used = len(phonon.dataset['first_atoms'])
    n_naive = phonon.dataset['natom'] * 6

    print_dual(f"Detected symmetry : space group {space_group}, point group {point_group} "
               f"-- {n_ops} operation(s)", f_out)
    print_dual(f"Displacements needed : {n_used} (vs. {n_naive} without symmetry reduction)", f_out)
    if n_naive > 0:
        print_dual(f"Reduction : {100 * (1 - n_used / n_naive):.1f}% fewer MACE evaluations", f_out)


def run_phonon_pipeline(atoms_ref, calc, args, f_out, tag=""):
    """Runs the full single-volume ML phonon pipeline for one calculator:
    optional pre-relax, displacement generation (+ symmetry-reduction
    report), force constants, acoustic-sum-rule correction
    (phonon.symmetrize_force_constants -- eliminates the spurious near-zero
    negative frequency at Gamma that plain finite-difference force
    constants otherwise carry from translational-invariance numerical
    noise, verified live: a real run showed -0.0003 THz at Gamma before
    this, essentially the acoustic sum rule violation this fixes), band
    structure, DOS + species-projected DOS, and thermal properties.

    `tag` is a short label (e.g. "fine-tuned", "foundation (small)") used
    only in printed progress/diagnostic messages, so a caller running this
    twice (fine-tuned vs. foundation comparison) can tell which run is
    which. Returns a dict with everything a caller might want to report,
    plot, or compare -- 'bs'/'band_labels'/'path_connections' are absent if
    the structure is fully vacuum-padded (no meaningful q-path).
    """
    atoms = atoms_ref.copy()
    atoms.calc = calc
    f0 = np.abs(atoms.get_forces()).max()
    print_dual(f"[{tag}] Residual force on input structure: {f0:.4f} eV/Ang", f_out)

    phonon, _, relax_info = mace_phonons.generate_ml_displacements(
        atoms, args.dim, args.distance, calc, relax=args.relax, fmax=args.fmax)
    if relax_info is not None:
        _, f1, converged, steps_used = relax_info
        print_dual(f"[{tag}] After relax: max|F| = {f1:.4f} eV/Ang "
                   f"({'converged' if converged else 'hit step cap, not fully converged'}, "
                   f"{steps_used} steps)", f_out)

    print_symmetry_table(phonon, f_out)

    n_disp = len(phonon.supercells_with_displacements)
    print_dual(f"[{tag}] Displacements needed: {n_disp}", f_out)
    report_every = max(1, n_disp // 10)

    def _progress(i, n):
        if i % report_every == 0 or i == n:
            print(f"  ... [{tag}] {i}/{n}")

    mace_phonons.compute_force_constants(phonon, calc, progress_callback=_progress)
    phonon.symmetrize_force_constants(show_drift=False)

    result = {"phonon": phonon}

    path = build_band_path(phonon.primitive, 1.0, args.vacuum_gap)
    if path is not None:
        kpoints_dict, path_segments, bravais_name = path
        bands, band_labels, path_connections = band_path_to_phonopy_format(
            kpoints_dict, path_segments, args.band_points)
        phonon.run_band_structure(bands, with_group_velocities=False,
                                  path_connections=path_connections,
                                  labels=band_labels, is_legacy_plot=False)
        result["bs"] = phonon.band_structure
        result["band_labels"] = band_labels
        result["path_connections"] = path_connections
        result["bravais_name"] = bravais_name
        result["path_segments"] = path_segments

    # with_eigenvectors + is_mesh_symmetry=False: needed for projected DOS
    # (per-atom character), same mesh call phonons_pos.py itself uses for
    # the same reason -- total DOS comes from the same mesh pass, no
    # redundant computation for the two.
    phonon.run_mesh(args.mesh, with_eigenvectors=True, is_mesh_symmetry=False)
    phonon.run_total_dos()
    dos_dict = phonon.get_total_dos_dict()
    result["dos_freqs"] = dos_dict['frequency_points']
    result["dos"] = dos_dict['total_dos']

    phonon.run_projected_dos()
    pdos_dict = phonon.get_projected_dos_dict()
    proj = pdos_dict['projected_dos']
    species_pdos = {}
    for i, sym in enumerate(phonon.primitive.symbols):
        species_pdos[sym] = species_pdos.get(sym, np.zeros_like(result["dos_freqs"])) + proj[i]
    result["species_pdos"] = species_pdos

    phonon.run_thermal_properties(t_min=args.tmin, t_max=args.tmax, t_step=args.tstep)
    tp = phonon.get_thermal_properties_dict()
    result["temps"] = tp['temperatures']
    result["free_energy"] = tp['free_energy']
    result["entropy"] = tp['entropy']
    result["heat_capacity"] = tp['heat_capacity']

    return result


def _pbc_atomic_shift(positions_new, positions_ref, cell):
    """Per-atom displacement magnitude between two position arrays of the
    same periodic cell, using the minimum-image convention (same convention
    as core/md_traj.py's trajectory unwrapping elsewhere in this suite) --
    NOT a plain Cartesian difference. A plain difference is wrong here: an
    atom sitting near a cell face can have its wrapped fractional
    coordinate flip from ~0.999 to ~0.001 between the two supercells for a
    displacement of a fraction of an Angstrom, which shows up as a jump of
    nearly a full lattice vector (tens of Angstrom) in raw Cartesian
    coordinates -- verified live: an amplitude-calibration probe this
    caused to read ~74 Ang instead of the true ~0.0002 Ang for a real
    unstable-mode test case, silently breaking the requested-amplitude
    scaling below.
    """
    frac_diff = (positions_new - positions_ref) @ np.linalg.inv(cell)
    frac_diff -= np.round(frac_diff)
    cart_diff = frac_diff @ cell
    return np.linalg.norm(cart_diff, axis=1)


def freeze_unstable_mode(phonon, args, species_meta, f_out):
    """If the DOS mesh (already computed with eigenvectors on a full,
    non-symmetry-reduced grid inside run_phonon_pipeline) contains a
    genuine imaginary mode -- below the same numerical-noise tolerance used
    for the band-path warning -- displaces the supercell along that mode's
    eigenvector and writes it out as a new structure, the exact escape-a-
    saddle-point trick stb-phononsPos's own --freeze-unstable-mode already
    does for the SIESTA path (phonon.run_modulations, calibrated to the
    requested Angstrom amplitude via a unit-amplitude probe displacement
    first -- the modulation amplitude parameter isn't itself in Angstrom).
    Unlike phonons_pos.py, no bohr conversion is needed on write: this
    tool's Phonopy object is built directly in Angstrom (see
    core/mace_phonons.py's module docstring), so the modulated supercell's
    positions/cell are already in the right convention for from_pymatgen/
    write_fdf. Returns the written path, or None if the mesh had no mode
    below the noise tolerance.
    """
    freqs = phonon.mesh.frequencies
    min_freq = freqs.min()
    if min_freq >= _IMAGINARY_MODE_TOL_THZ:
        print_dual("No imaginary mode on the sampled mesh -- nothing to freeze.", f_out)
        return None

    q_idx, band_idx = np.unravel_index(np.argmin(freqs), freqs.shape)
    q_point = phonon.mesh.qpoints[q_idx]

    phonon.run_modulations(args.mesh, [[q_point, band_idx, 1.0, 0.0]])
    mod = phonon.modulation
    probe_shift = _pbc_atomic_shift(
        mod.modulated_supercells[0].positions, mod.supercell.positions,
        np.array(mod.supercell.cell)).max()
    scale = args.freeze_amplitude / probe_shift

    phonon.run_modulations(args.mesh, [[q_point, band_idx, scale, 0.0]])
    mod = phonon.modulation
    achieved_shift = _pbc_atomic_shift(
        mod.modulated_supercells[0].positions, mod.supercell.positions,
        np.array(mod.supercell.cell)).max()

    frozen = mod.modulated_supercells[0]
    frozen_atoms = Atoms(numbers=frozen.numbers, positions=frozen.positions,
                         cell=frozen.cell, pbc=True)
    pmg = AseAtomsAdaptor.get_structure(frozen_atoms)
    new_structure = structure_io.from_pymatgen(pmg, species_meta=species_meta)
    out_path = os.path.join(args.output_dir, "frozen_mode.fdf")
    structure_io.write_fdf(new_structure, out_path)

    print_dual(f"Softest mode      : q = {np.array2string(np.array(q_point), precision=4)}, "
               f"band {band_idx}, {min_freq:.4f} THz", f_out)
    print_dual(f"Displacement      : requested {args.freeze_amplitude:.4f} Ang, "
               f"achieved {achieved_shift:.4f} Ang (max atomic shift)", f_out)
    print_dual(f"Structure written : {out_path}", f_out)
    return out_path


def animate_normal_mode(phonon, qpoint, band_index, amplitude, n_frames, dim):
    """Builds a multi-frame ASE trajectory animating one normal mode's
    vibration at fractional q-point `qpoint`, by calling phonon.
    run_modulations once per frame with the phase swept over a full 360
    -degree cycle -- same phonopy modulation machinery freeze_unstable_mode
    uses for a single displaced snapshot, sampled here at many phases
    instead, purely for visualization (OVITO/VMD), not a physical
    perturbation. `dim` is the modulation supercell size and must make
    `qpoint` commensurate (q times dim ~ integer vector) -- Gamma is
    commensurate with any dim, [1, 1, 1] included; a non-Gamma q needs a
    correspondingly larger dim. phonon.get_frequencies(q) diagonalizes the
    dynamical matrix directly at the given q (no mesh/supercell needed for
    this part), so it works for any q regardless of `dim`, letting the
    out-of-range check below run before wasting time on a bad request.

    `amplitude` is calibrated to true Angstrom via the same unit-amplitude
    probe + rescale trick freeze_unstable_mode uses -- phonopy's own
    "amplitude" parameter is an internal eigenvector-normalization scale,
    not itself in Angstrom (verified live: amplitude=1.0 gave a ~0.094 Ang
    max shift for a real optical mode, not 1 Ang), so passing the user's
    requested Angstrom value straight through would silently produce the
    wrong displacement, just like the pre-fix freeze_unstable_mode bug this
    mirrors. The scale is amplitude-only (not phase-dependent, since phase
    is applied as a separate multiplicative factor in phonopy's own
    modulation formula), so one probe covers every frame in the cycle.

    Returns (frames, n_bands_at_q, frequency_thz) -- frames is None (with
    n_bands_at_q still valid) if band_index is out of range, so the caller
    can print a friendlier error.
    """
    frequencies_at_q = phonon.get_frequencies(np.array(qpoint))
    n_bands = len(frequencies_at_q)
    if not (0 <= band_index < n_bands):
        return None, n_bands, None

    phonon.run_modulations(dim, [[qpoint, band_index, 1.0, 0.0]])
    mod = phonon.modulation
    probe_shift = _pbc_atomic_shift(
        mod.modulated_supercells[0].positions, mod.supercell.positions,
        np.array(mod.supercell.cell)).max()
    scale = amplitude / probe_shift if probe_shift > 1e-12 else 0.0

    frames = []
    for phase_deg in np.linspace(0.0, 360.0, n_frames, endpoint=False):
        phonon.run_modulations(dim, [[qpoint, band_index, scale, phase_deg]])
        mod = phonon.modulation
        cell_atoms = mod.modulated_supercells[0]
        frames.append(Atoms(numbers=cell_atoms.numbers, positions=cell_atoms.positions,
                            cell=cell_atoms.cell, pbc=True))

    return frames, n_bands, frequencies_at_q[band_index]


def plot_bands(series, out_path):
    """series: list of (label, result_dict, color). Overlays every series
    on one figure (a single entry for the normal case, two for a fine-tuned
    -vs-foundation comparison -- same q-path by construction, since both
    come from the same structure/--dim/--band-points, so only the first
    series' tick positions are needed for the shared x-axis).
    """
    fig, ax = plt.subplots(figsize=(8, 5))
    for label, result, color in series:
        bs = result["bs"]
        first = True
        for seg_dist, seg_freq in zip(bs.distances, bs.frequencies):
            for band in seg_freq.T:
                ax.plot(seg_dist, band, color=color, linewidth=0.8,
                       label=label if first else None)
                first = False
    ax.axhline(0.0, color='gray', linestyle='--', linewidth=0.8)

    _, result0, _ = series[0]
    bs0 = result0["bs"]
    ticks = band_tick_positions(bs0.distances, result0["band_labels"], result0["path_connections"])
    ax.set_xticks([t[0] for t in ticks])
    ax.set_xticklabels([pretty_label(t[1]) for t in ticks])
    for t in ticks:
        ax.axvline(t[0], color='gray', linewidth=0.5, alpha=0.5)
    ax.set_xlim(bs0.distances[0][0], bs0.distances[-1][-1])
    ax.set_ylabel("Frequency (THz)")
    ax.set_title("Phonon band structure" if len(series) == 1
                 else "Phonon band structure: fine-tuned vs. foundation")
    if len(series) > 1:
        handles, labels_ = ax.get_legend_handles_labels()
        seen = set()
        uniq = [(h, l) for h, l in zip(handles, labels_) if not (l in seen or seen.add(l))]
        ax.legend([h for h, _ in uniq], [l for _, l in uniq], fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_dos(series, out_path):
    """series: list of (label, result_dict, color). Per-species PDOS
    breakdown is only drawn for a single-series (non-comparison) run --
    with two full models overlaid, per-species curves for both would be
    too cluttered to read; the comparison focuses on the total DOS."""
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for label, result, color in series:
        ax.plot(result["dos_freqs"], result["dos"],
                label=f"{label} (total)" if len(series) > 1 else "Total",
                color=color, linewidth=1.5)
        if len(series) == 1:
            for sym, pdos in result["species_pdos"].items():
                ax.plot(result["dos_freqs"], pdos, label=sym, linewidth=0.8, linestyle='--')
    ax.set_xlabel("Frequency (THz)")
    ax.set_ylabel("DOS (states/THz/cell)")
    ax.set_title("Phonon density of states" if len(series) == 1
                 else "Phonon DOS: fine-tuned vs. foundation")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_thermal(series, out_path):
    """series: list of (label, result_dict, color)."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    for label, result, color in series:
        axes[0].plot(result["temps"], result["free_energy"], label=label, color=color)
        axes[1].plot(result["temps"], result["entropy"], label=label, color=color)
        axes[2].plot(result["temps"], result["heat_capacity"], label=label, color=color)
    axes[0].set_xlabel("Temperature (K)")
    axes[0].set_ylabel("Free energy (kJ/mol)")
    axes[0].set_title("Free energy")
    axes[1].set_xlabel("Temperature (K)")
    axes[1].set_ylabel("Entropy (J/K/mol)")
    axes[1].set_title("Entropy")
    axes[2].set_xlabel("Temperature (K)")
    axes[2].set_ylabel("Heat capacity (J/K/mol)")
    axes[2].set_title("Heat capacity (Cv)")
    if len(series) > 1:
        axes[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def run_qha(atoms_ref, calc, args, f_out):
    """Quasi-harmonic approximation: runs the displacement->force-constant
    ->thermal-properties part of the pipeline at --n-volumes isotropically
    -scaled volumes (+/- --strain-range %), then fits phonopy's own
    PhonopyQHA (equation of state + per-temperature thermal expansion/bulk
    modulus/Gibbs free energy) -- adapted directly from the now-deleted
    phonons_qha.py. No foundation-model comparison in this mode (out of
    scope -- would double an already --n-volumes-times-repeated
    calculation); no per-volume PDOS either (QHA's focus is the equation
    of state, not vibrational character). Returns (qha, volumes, energies,
    temperatures_ref).
    """
    from phonopy.api_qha import PhonopyQHA

    v0_input = atoms_ref.get_volume()
    strains_pct = np.linspace(-abs(args.strain_range), abs(args.strain_range), args.n_volumes)
    # PhonopyQHA's t_max must sit a few grid points short of the actual end
    # of the temperature array (needs neighbors for finite-difference
    # derivatives) -- pad the per-volume grid past the user's --tmax.
    t_max_grid = args.tmax + 3 * args.tstep

    print_dual(f"Volume range      : {-abs(args.strain_range):.1f}% to "
               f"{abs(args.strain_range):.1f}% ({args.n_volumes} points, V0={v0_input:.4f} Ang^3)", f_out)
    print_dual(f"{'Strain%':>8} {'Volume(Ang^3)':>14} {'Energy(eV)':>14} "
               f"{'N_disp':>7} {'MinFreq(THz)':>13}", f_out)

    volumes, energies = [], []
    free_energy_list, entropy_list, cv_list = [], [], []
    temperatures_ref = None
    any_unstable = False

    for i, strain_pct in enumerate(strains_pct):
        factor = (1.0 + strain_pct / 100.0) ** (1.0 / 3.0)
        atoms_i = atoms_ref.copy()
        atoms_i.set_cell(np.array(atoms_ref.get_cell()) * factor, scale_atoms=True)
        atoms_i.calc = calc

        def _progress(n_done, n_total, _i=i):
            report_every = max(1, n_total // 5)
            if n_done % report_every == 0 or n_done == n_total:
                print(f"  ... volume {_i + 1}/{len(strains_pct)}: {n_done}/{n_total} forces")

        phonon, relaxed_atoms, _ = mace_phonons.generate_ml_displacements(
            atoms_i, args.dim, args.distance, calc, relax=args.relax, fmax=args.fmax)
        mace_phonons.compute_force_constants(phonon, calc, progress_callback=_progress)
        phonon.symmetrize_force_constants(show_drift=False)

        energy_i = relaxed_atoms.get_potential_energy()
        volume_i = relaxed_atoms.get_volume()

        phonon.run_mesh(args.mesh)
        min_freq_i = phonon.mesh.frequencies.min()
        phonon.run_thermal_properties(t_min=args.tmin, t_max=t_max_grid, t_step=args.tstep)
        tp = phonon.thermal_properties
        if temperatures_ref is None:
            temperatures_ref = tp.temperatures

        vol_dir = os.path.join(args.output_dir, f"vol_{i:02d}")
        os.makedirs(vol_dir, exist_ok=True)
        phonon.save(os.path.join(vol_dir, "phonopy_disp.yaml"), settings={'force_constants': True})

        volumes.append(volume_i)
        energies.append(energy_i)
        free_energy_list.append(tp.free_energy)
        entropy_list.append(tp.entropy)
        cv_list.append(tp.heat_capacity)

        n_disp = len(phonon.dataset['first_atoms'])
        freq_flag = color_text(f"{min_freq_i:13.4f}", 'red' if min_freq_i < _IMAGINARY_MODE_TOL_THZ else 'green')
        if min_freq_i < _IMAGINARY_MODE_TOL_THZ:
            any_unstable = True
        print_dual(f"{strain_pct:8.2f} {volume_i:14.4f} {energy_i:14.6f} "
                   f"{n_disp:7d} {freq_flag}", f_out)

    if any_unstable:
        print_dual(color_text(
            "[WARNING] At least one volume has an imaginary (negative) frequency -- QHA "
            "assumes every sampled point is dynamically stable. Results derived across "
            "the fit may not be physically meaningful.", 'red'), f_out)

    free_energy = np.array(free_energy_list).T
    entropy = np.array(entropy_list).T
    cv = np.array(cv_list).T

    qha = PhonopyQHA(volumes=volumes, electronic_energies=energies,
                     temperatures=temperatures_ref, free_energy=free_energy,
                     cv=cv, entropy=entropy, eos=args.eos, t_max=args.tmax)
    return qha, volumes, energies, temperatures_ref


def _max_frequency(result):
    """Max band-path frequency if a q-path exists (vacuum-padded/0D
    structures have none), else the max mesh frequency -- same fallback
    both run_convergence_check and plot_convergence need."""
    if "bs" in result:
        return np.concatenate(result["bs"].frequencies).max()
    return result["phonon"].mesh.frequencies.max()


def run_convergence_check(atoms_ref, calc, args, conv_dims, f_out):
    """Runs the full single-volume pipeline once per supercell size in
    `conv_dims` (reusing run_phonon_pipeline unchanged, just with --dim
    overridden per point) and reports how the maximum band-path/mesh
    frequency and the heat capacity/entropy at a reference temperature
    change between them -- MACE evaluations are cheap enough that this is
    an actually affordable convergence scan, unlike the real SIESTA-based
    workflow where a scan of this kind rarely is. `conv_dims` is a list of
    (na, nb, nc) tuples; returns a list of (dim_tuple, result_dict) for the
    caller to plot/save.
    """
    import copy

    ref_t = min(args.tmax, 300.0)
    rows = []
    print_dual(f"{'Supercell':>12} {'N_disp':>7} {'MaxFreq(THz)':>13} "
               f"{'Cv@' + str(int(ref_t)) + 'K':>12} {'S@' + str(int(ref_t)) + 'K':>12}", f_out)
    for dim_tuple in conv_dims:
        dim_args = copy.copy(args)
        dim_args.dim = list(dim_tuple)
        tag = "x".join(str(d) for d in dim_tuple)
        result = run_phonon_pipeline(atoms_ref, calc, dim_args, f_out, tag=tag)
        n_disp = len(result["phonon"].dataset['first_atoms'])
        max_freq = _max_frequency(result)
        t_idx = int(np.argmin(np.abs(result["temps"] - ref_t)))
        cv_ref = result["heat_capacity"][t_idx]
        s_ref = result["entropy"][t_idx]
        print_dual(f"{tag:>12} {n_disp:7d} {max_freq:13.4f} {cv_ref:12.4f} {s_ref:12.4f}", f_out)
        rows.append((dim_tuple, result))
    return rows


def plot_convergence(rows, ref_t, out_path):
    """rows: list of (dim_tuple, result_dict) from run_convergence_check.
    x-axis is the supercell's atom-count multiplier (na*nb*nc) so unevenly
    -spaced dim choices (e.g. 1x1x1, 2x2x2, 3x3x3) still land at a
    meaningful, monotonic position instead of just a bare category index.
    """
    labels = ["x".join(str(d) for d in dim_tuple) for dim_tuple, _ in rows]
    n_cells = [int(np.prod(dim_tuple)) for dim_tuple, _ in rows]
    max_freqs = [_max_frequency(result) for _, result in rows]
    cvs = []
    for _, result in rows:
        t_idx = int(np.argmin(np.abs(result["temps"] - ref_t)))
        cvs.append(result["heat_capacity"][t_idx])

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    axes[0].plot(n_cells, max_freqs, marker='o')
    axes[0].set_xticks(n_cells)
    axes[0].set_xticklabels(labels, rotation=30, ha='right')
    axes[0].set_ylabel("Max frequency (THz)")
    axes[0].set_title("Frequency convergence")
    axes[1].plot(n_cells, cvs, marker='o', color='tab:orange')
    axes[1].set_xticks(n_cells)
    axes[1].set_xticklabels(labels, rotation=30, ha='right')
    axes[1].set_ylabel(f"Heat capacity @ {ref_t:.0f} K (J/K/mol)")
    axes[1].set_title("Thermal-property convergence")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Standalone phonon calculation driven entirely by a MACE potential -- the "
                     "MACE-MP-0 foundation model, or a custom model fine-tuned on your own "
                     "SIESTA data via stb-mlffAnalysis (--custom-model). One shot: displacements, "
                     "force constants, band structure, DOS (total + per-species), and thermal "
                     "properties -- or, with --qha, a quasi-harmonic volume scan (thermal "
                     "expansion, bulk modulus vs. T). No SIESTA and no separate analysis tool "
                     "needed. A fast heuristic preview, not a substitute for the real SIESTA-"
                     "based stb-phononsCreate/stb-phononsPos workflow -- ML potentials are "
                     "markedly less reliable for phonons (imaginary modes especially) than DFT.",
        epilog="Example usage:\n"
               "  stb-mlphonons --file structure.fdf --dim 2 2 2\n"
               "  stb-mlphonons --file structure.fdf --custom-model mlff_model.model --dim 3 3 3\n"
               "  stb-mlphonons --file structure.fdf --qha --n-volumes 7 --strain-range 5",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--file", required=True, help="Input structure (.fdf).")
    parser.add_argument("--model", choices=["small", "medium", "large"], default="small",
                        help="MACE-MP-0 model size (default: small). Ignored if --custom-model is given.")
    parser.add_argument("--custom-model", default=None, metavar="PATH",
                        help="Use a custom MACE model file instead of the MACE-MP-0 foundation "
                             "potential -- e.g. one fine-tuned via stb-mlffAnalysis. Overrides --model.")
    parser.add_argument("--skip-foundation-comparison", action="store_true",
                        help="If --custom-model is given, also evaluate the same phonon "
                             "pipeline with the raw (non-fine-tuned) foundation model "
                             "(--model) and overlay both on the bands/DOS/thermal plots, "
                             "to show whether fine-tuning actually changed the phonon "
                             "prediction. On by default; pass this to skip it (--qha never "
                             "does this comparison, regardless of this flag).")
    parser.add_argument("--dim", type=int, nargs=3, default=[2, 2, 2], metavar=("NA", "NB", "NC"),
                        help="Supercell dimensions for the finite-displacement method (default: 2 2 2).")
    parser.add_argument("--distance", type=float, default=0.01,
                        help="Displacement distance, Angstrom (default: 0.01).")
    parser.add_argument("--no-relax", dest="relax", action="store_false",
                        help="Skip the MACE pre-relax (positions only) before generating "
                             "displacements. On by default -- a phonon calculation assumes "
                             "~zero net force at the reference geometry; skipping this on an "
                             "unrelaxed structure is a common cause of spurious imaginary modes.")
    parser.add_argument("--fmax", type=float, default=0.05,
                        help="Force convergence threshold for the pre-relax, eV/Ang (default: 0.05).")
    parser.add_argument("--vacuum-gap", type=float, default=10.0,
                        help="Gap (Ang) used to detect vacuum-padded axes (default: 10.0).")
    parser.add_argument("--band-points", type=int, default=101,
                        help="q-points per band-structure segment (default: 101).")
    parser.add_argument("--mesh", type=int, nargs=3, default=[20, 20, 20], metavar=("NA", "NB", "NC"),
                        help="q-point mesh for DOS/thermal properties (default: 20 20 20).")
    parser.add_argument("--tmin", type=float, default=0.0, help="Min temperature, K (default: 0).")
    parser.add_argument("--tmax", type=float, default=1000.0, help="Max temperature, K (default: 1000).")
    parser.add_argument("--tstep", type=float, default=10.0, help="Temperature step, K (default: 10).")
    parser.add_argument("--qha", action="store_true",
                        help="Quasi-harmonic approximation instead of a single-volume "
                             "calculation: scans --n-volumes isotropically-scaled volumes and "
                             "fits an equation of state (thermal expansion, bulk modulus vs. "
                             "T). No foundation-model comparison in this mode. Off by default.")
    parser.add_argument("--n-volumes", type=int, default=5,
                        help="Number of volume points for --qha, >= 4 (default: 5).")
    parser.add_argument("--strain-range", type=float, default=3.0,
                        help="Volume range around equilibrium for --qha, +/- %% (default: 3.0).")
    parser.add_argument("--eos", choices=["vinet", "murnaghan", "birch_murnaghan"], default="vinet",
                        help="Equation of state for --qha (default: vinet).")
    parser.add_argument("--check-convergence", action="store_true",
                        help="Instead of a single-volume calculation, run the same pipeline "
                             "at multiple supercell sizes (--convergence-dims) and report how "
                             "the max frequency / heat capacity / entropy change between them "
                             "-- MACE is cheap enough to check --dim convergence directly, "
                             "unlike DFT. Not supported together with --qha. Off by default.")
    parser.add_argument("--convergence-dims", type=int, nargs='+', default=None, metavar="N",
                        help="Flat list of supercell sizes for --check-convergence, given as "
                             "consecutive NA NB NC triples (e.g. 1 1 1 2 2 2 3 3 3 tests "
                             "1x1x1, 2x2x2, and 3x3x3) -- must be a multiple of 3 values. "
                             "Default: --dim itself and --dim+1 on every non-vacuum axis "
                             "(two points).")
    parser.add_argument("--freeze-unstable-mode", action="store_true",
                        help="If the phonon mesh has a genuine imaginary (unstable) mode, "
                             "write a new structure displaced along that mode's eigenvector -- "
                             "a common trick to escape a saddle point and find the true "
                             "minimum. Same feature stb-phononsPos already has for the SIESTA "
                             "path. Not supported together with --qha/--check-convergence. "
                             "Off by default.")
    parser.add_argument("--freeze-amplitude", type=float, default=0.05,
                        help="Max atomic displacement (Angstrom) for --freeze-unstable-mode "
                             "(default: 0.05).")
    parser.add_argument("--animate-mode", type=int, default=None, metavar="BAND_INDEX",
                        help="Animate a normal mode's vibration and write it as a multi-frame "
                             "trajectory (xsf/pdb/xyz, same convention as stb-ani2traj/"
                             "stb-mlmd) for viewing in OVITO/VMD -- BAND_INDEX selects which "
                             "mode (0 = lowest frequency) at --animate-qpoint. Not supported "
                             "together with --qha/--check-convergence. Off by default.")
    parser.add_argument("--animate-qpoint", type=float, nargs=3, default=[0.0, 0.0, 0.0],
                        metavar=("QA", "QB", "QC"),
                        help="Fractional q-point to animate a mode at (default: 0 0 0, Gamma).")
    parser.add_argument("--animate-dim", type=int, nargs=3, default=[1, 1, 1],
                        metavar=("NA", "NB", "NC"),
                        help="Supercell size used to build the animated modulation -- must "
                             "make --animate-qpoint commensurate (q x dim ~ integer vector); "
                             "the default 1 1 1 only works exactly at Gamma (default: 1 1 1).")
    parser.add_argument("--animate-amplitude", type=float, default=0.5,
                        help="Max atomic displacement for the animation, Angstrom -- "
                             "deliberately larger than --freeze-amplitude since this is for "
                             "visualization, not a physical perturbation (default: 0.5).")
    parser.add_argument("--animate-frames", type=int, default=30,
                        help="Number of frames in one full vibration cycle (default: 30).")
    parser.add_argument("--animate-format", choices=sorted(OUTPUT_FORMATS), default="xsf",
                        help="Output trajectory format for --animate-mode, same convention as "
                             "stb-ani2traj/stb-mlmd (default: xsf).")
    parser.add_argument("--device", choices=["cpu", "cuda"], default="cpu",
                        help="Device to run the model on (default: cpu).")
    parser.add_argument("-o", "--output-dir", default="mlphonons_out",
                        help="Output directory for all plots/data (default: mlphonons_out).")
    parser.add_argument("--save-data", action="store_true",
                        help="Also write the raw numeric data (.dat) behind each plot. Off by default.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-mlphonons {VERSION}")

    args = parser.parse_args()

    if args.qha and args.n_volumes < 4:
        parser.error("--n-volumes must be at least 4 -- vinet/murnaghan/birch_murnaghan all "
                     "fit 4 parameters (E0, V0, B0, B0'), so need at least 4 volumes.")

    for flag_name, flag_val in (("--freeze-unstable-mode", args.freeze_unstable_mode),
                                 ("--animate-mode", args.animate_mode is not None)):
        if args.qha and flag_val:
            parser.error(f"{flag_name} is not supported together with --qha.")
        if args.check_convergence and flag_val:
            parser.error(f"{flag_name} is not supported together with --check-convergence.")
    if args.qha and args.check_convergence:
        parser.error("--qha is not supported together with --check-convergence.")

    conv_dims_override = None
    if args.convergence_dims is not None:
        if len(args.convergence_dims) % 3 != 0:
            parser.error("--convergence-dims must be a multiple of 3 values (NA NB NC per point).")
        flat = args.convergence_dims
        conv_dims_override = [tuple(flat[i:i + 3]) for i in range(0, len(flat), 3)]
        if len(conv_dims_override) < 2:
            parser.error("--check-convergence needs at least 2 supercell sizes to compare.")

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - ML Phonons",
            "Bands, DOS, and thermal properties, driven entirely by a MACE potential",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("ML PHONONS:", 'bold'))
    print("-" * 60)

    os.makedirs(args.output_dir, exist_ok=True)
    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    def fail(message):
        print_dual(color_text(f"[FAIL] {message}", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_dual(color_text("===== STB-MLPHONONS REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Input file        : {args.file}", f_out)
    mode_label = ("QHA (volume scan)" if args.qha
                  else "supercell convergence check" if args.check_convergence
                  else "single-volume")
    print_dual(f"Mode              : {mode_label}", f_out)
    print_dual(f"Supercell dim     : {args.dim[0]} x {args.dim[1]} x {args.dim[2]}", f_out)
    print_dual(f"Displacement dist.: {args.distance} Ang", f_out)
    print_dual(f"Model             : "
               f"{'custom (' + args.custom_model + ')' if args.custom_model else f'MACE-MP-0 ({args.model})'}", f_out)
    compare = bool(args.custom_model) and not args.skip_foundation_comparison and not args.qha
    if compare:
        print_dual(f"Comparison        : also evaluating with foundation model ({args.model})", f_out)
    if report_path:
        print_dual(f"Report file       : {report_path}", f_out)

    print_section("[1] READING STRUCTURE", f_out)

    if not os.path.isfile(args.file):
        fail(f"Input file '{args.file}' not found.")
    if args.custom_model and not os.path.isfile(args.custom_model):
        fail(f"--custom-model file not found: {args.custom_model}")

    try:
        structure = structure_io.read_fdf(args.file)
    except Exception as e:
        fail(f"Could not read '{args.file}': {e}")

    pmg_structure = structure_io.to_pymatgen(structure)
    atoms = AseAtomsAdaptor.get_atoms(pmg_structure)

    lattice_ang = np.array(atoms.get_cell())
    frac_coords = atoms.get_scaled_positions()
    vacuum_axes = kspace.detect_vacuum_axes(frac_coords, lattice_ang, args.vacuum_gap)
    print_dual(f"[OK] Read {len(atoms)} atom(s), "
               f"dimensionality: {kspace.dimensionality_label(vacuum_axes)}", f_out)

    vacuum_dims_requested = [axis for axis, is_vac in zip('abc', vacuum_axes)
                              if is_vac and args.dim['abc'.index(axis)] > 1]
    if vacuum_dims_requested:
        print_dual(color_text(
            f"[WARNING] --dim requests more than 1 repetition along vacuum-padded "
            f"axis/axes {', '.join(vacuum_dims_requested)} (gap >= {args.vacuum_gap} Ang). "
            "Replicating a supercell across vacuum only multiplies the MACE cost "
            "without adding real periodicity -- consider --dim 1 on that axis.", 'yellow'), f_out)

    model_arg = args.custom_model if args.custom_model else args.model
    calc = mace_relax.get_calculator(model=model_arg, device=args.device)
    if not args.relax:
        print_dual(color_text(
            "[WARNING] --no-relax: generating displacements from the structure as given. "
            "If it isn't already at a MACE-relaxed minimum, expect spurious imaginary "
            "modes downstream.", 'yellow'), f_out)

    if args.check_convergence:
        if conv_dims_override is not None:
            conv_dims = conv_dims_override
        else:
            base = list(args.dim)
            bumped = [d if is_vac else d + 1 for d, is_vac in zip(base, vacuum_axes)]
            conv_dims = [tuple(base), tuple(bumped)]

        print_section("[2] SUPERCELL CONVERGENCE CHECK", f_out)
        ref_t = min(args.tmax, 300.0)
        rows = run_convergence_check(atoms, calc, args, conv_dims, f_out)

        print_section("[3] SUMMARY & FILES", f_out)
        conv_plot = os.path.join(args.output_dir, "convergence.png")
        plot_convergence(rows, ref_t, conv_plot)
        print_dual(f"[OK] Saved {conv_plot}", f_out)
        if args.save_data:
            conv_data = os.path.join(args.output_dir, "convergence.dat")
            with open(conv_data, "w") as f:
                f.write("# dim n_disp max_freq(THz) Cv_ref(J/K/mol) S_ref(J/K/mol)\n")
                for dim_tuple, result in rows:
                    tag = "x".join(str(d) for d in dim_tuple)
                    n_disp = len(result["phonon"].dataset['first_atoms'])
                    max_freq = _max_frequency(result)
                    t_idx = int(np.argmin(np.abs(result["temps"] - ref_t)))
                    f.write(f"{tag} {n_disp} {max_freq:.6f} {result['heat_capacity'][t_idx]:.6f} "
                            f"{result['entropy'][t_idx]:.6f}\n")
            print_dual(f"[OK] Saved {conv_data}", f_out)
        print_dual("Status            : OK", f_out)
        print_dual(f"Convergence plot  : {conv_plot}", f_out)
        if report_path:
            print_dual(f"Report            : {report_path}", f_out)
        if f_out:
            f_out.close()
        print("\n" + "-" * 60)
        print(color_text("ML phonon supercell convergence check complete.\n", 'bold'))
        return

    if args.qha:
        print_section("[2] VOLUME SCAN", f_out)
        qha, volumes, energies, temperatures_ref = run_qha(atoms, calc, args, f_out)

        print_section("[3] EQUILIBRIUM FIT (T=0, static)", f_out)
        e0, b0, b0p, v0_fit = qha.get_bulk_modulus_parameters()
        print_dual(f"E0 (min. energy)  : {e0:.6f} eV", f_out)
        print_dual(f"V0 (eq. volume)   : {v0_fit:.6f} Ang^3", f_out)
        print_dual(f"B0 (bulk modulus) : {b0:.4f} GPa", f_out)
        print_dual(f"B0' (pressure derivative): {b0p:.4f}", f_out)

        print_section("[4] THERMAL EXPANSION SUMMARY", f_out)
        temps_qha = qha.thermal_expansion
        n_t = len(temps_qha)
        t_grid_trunc = temperatures_ref[:n_t]
        ref_idx = int(np.argmin(np.abs(t_grid_trunc - 300.0)))
        last_idx = n_t - 1
        for label, idx in (("T ~ 300 K", ref_idx), (f"T = {t_grid_trunc[last_idx]:.0f} K", last_idx)):
            print_dual(f"At {label}:", f_out)
            print_dual(f"  Volume            : {qha.volume_temperature[idx]:.4f} Ang^3", f_out)
            print_dual(f"  Thermal expansion : {qha.thermal_expansion[idx]:.6e} 1/K", f_out)
            print_dual(f"  Bulk modulus      : {qha.bulk_modulus_temperature[idx]:.4f} GPa", f_out)
            print_dual(f"  Gibbs free energy : {qha.gibbs_temperature[idx]:.6f} eV", f_out)

        qha_plot = os.path.join(args.output_dir, "qha.png")
        fig, axes = plt.subplots(2, 2, figsize=(11, 9))
        axes[0, 0].scatter(volumes, energies)
        axes[0, 0].set_xlabel("Volume (Ang^3)")
        axes[0, 0].set_ylabel("Energy (eV)")
        axes[0, 0].set_title("EOS fit points")
        axes[0, 1].plot(t_grid_trunc, qha.thermal_expansion)
        axes[0, 1].set_xlabel("Temperature (K)")
        axes[0, 1].set_ylabel("alpha (1/K)")
        axes[0, 1].set_title("Thermal expansion")
        axes[1, 0].plot(t_grid_trunc, qha.volume_temperature)
        axes[1, 0].set_xlabel("Temperature (K)")
        axes[1, 0].set_ylabel("Volume (Ang^3)")
        axes[1, 0].set_title("Equilibrium volume")
        axes[1, 1].plot(t_grid_trunc, qha.bulk_modulus_temperature)
        axes[1, 1].set_xlabel("Temperature (K)")
        axes[1, 1].set_ylabel("Bulk modulus (GPa)")
        axes[1, 1].set_title("Bulk modulus")
        fig.tight_layout()
        fig.savefig(qha_plot, dpi=150)
        plt.close(fig)
        print_dual(f"[OK] Saved {qha_plot}", f_out)
        if args.save_data:
            qha_data = os.path.join(args.output_dir, "qha_thermal_expansion.dat")
            np.savetxt(qha_data, np.column_stack([
                t_grid_trunc, qha.thermal_expansion, qha.volume_temperature, qha.bulk_modulus_temperature]),
                      header="T(K) alpha(1/K) Volume(Ang^3) BulkModulus(GPa)")
            print_dual(f"[OK] Saved {qha_data}", f_out)

        print_section("[5] SUMMARY & FILES", f_out)
        print_dual("Status            : OK", f_out)
        print_dual(f"QHA plot          : {qha_plot}", f_out)
        if report_path:
            print_dual(f"Report            : {report_path}", f_out)
        if f_out:
            f_out.close()
        print("\n" + "-" * 60)
        print(color_text("ML QHA complete -- a fast heuristic, not a substitute for DFT phonons.\n", 'bold'))
        return

    print_section("[2] DISPLACEMENTS & FORCE CONSTANTS", f_out)
    result_main = run_phonon_pipeline(atoms, calc, args, f_out, tag="fine-tuned" if args.custom_model else args.model)
    yaml_path = os.path.join(args.output_dir, "phonopy_disp.yaml")
    result_main["phonon"].save(yaml_path, settings={'force_constants': True})
    print_dual(f"[OK] Force constants computed (ASR-corrected). Saved {yaml_path}", f_out)

    result_found = None
    if compare:
        print_dual(f"\n[INFO] Also running the raw foundation model ({args.model}) for comparison...", f_out)
        found_calc = mace_relax.get_calculator(model=args.model, device=args.device)
        result_found = run_phonon_pipeline(atoms, found_calc, args, f_out, tag=f"foundation ({args.model})")

    print_section("[3] BAND STRUCTURE", f_out)
    bands_plot = None
    if "bs" not in result_main:
        print_dual(color_text(
            "Skipped: every axis is vacuum-padded (0D/isolated system) -- a q-path "
            "isn't physically meaningful here.", 'yellow'), f_out)
    else:
        print_dual(f"Bravais lattice   : {result_main['bravais_name']}", f_out)
        path_str = " | ".join("-".join(pretty_label(l) for l in seg) for seg in result_main['path_segments'])
        print_dual(f"Path              : {path_str}", f_out)

        for tag, result in (("Fine-tuned" if args.custom_model else f"MACE-MP-0 ({args.model})", result_main),
                            (f"Foundation ({args.model})", result_found)):
            if result is None:
                continue
            band_min_freq = np.concatenate(result["bs"].frequencies).min()
            print_dual(f"[{tag}] Minimum band-path frequency: {band_min_freq:.4f} THz", f_out)
            if band_min_freq < _IMAGINARY_MODE_TOL_THZ:
                print_dual(color_text(
                    f"[WARNING] [{tag}] Negative (imaginary) frequency found along the band "
                    "path -- sign of an unstable/unconverged structure, or a limitation of "
                    "the ML potential for this system. Try a tighter --fmax, or verify with "
                    "the real SIESTA-based phonon workflow.", 'red'), f_out)
            else:
                print_dual(color_text(f"[{tag}] No imaginary modes found along the band path.", 'green'), f_out)

        series = [(("Fine-tuned" if args.custom_model else f"MACE-MP-0 ({args.model})"), result_main, _SERIES_COLORS[0])]
        if result_found is not None:
            series.append((f"Foundation ({args.model})", result_found, _SERIES_COLORS[1]))

        bands_plot = os.path.join(args.output_dir, "bands.png")
        plot_bands(series, bands_plot)
        print_dual(f"[OK] Saved {bands_plot}", f_out)
        if args.save_data:
            bands_data = os.path.join(args.output_dir, "bands.dat")
            with open(bands_data, "w") as f:
                f.write("# distance(1/Ang) freq_band1 freq_band2 ... (THz, blank line = segment break)\n")
                for seg_dist, seg_freq in zip(result_main["bs"].distances, result_main["bs"].frequencies):
                    for d, row in zip(seg_dist, seg_freq):
                        f.write(f"{d:.6f} " + " ".join(f"{fr:.6f}" for fr in row) + "\n")
                    f.write("\n")
            print_dual(f"[OK] Saved {bands_data}", f_out)

    if args.freeze_unstable_mode:
        print_section("[3b] MODE FREEZE", f_out)
        freeze_unstable_mode(result_main["phonon"], args, structure.species_meta, f_out)

    if args.animate_mode is not None:
        print_section("[3c] MODE ANIMATION", f_out)
        frames, n_bands, freq = animate_normal_mode(
            result_main["phonon"], args.animate_qpoint, args.animate_mode,
            args.animate_amplitude, args.animate_frames, args.animate_dim)
        if frames is None:
            print_dual(color_text(
                f"[WARNING] --animate-mode {args.animate_mode} out of range -- only {n_bands} "
                f"band(s) at q = {args.animate_qpoint}. Skipping animation.", 'red'), f_out)
        else:
            ase_format, ext = OUTPUT_FORMATS[args.animate_format]
            anim_path = os.path.join(args.output_dir, f"mode_{args.animate_mode}{ext}")
            ase_write(anim_path, frames, format=ase_format)
            print_dual(f"q-point           : {args.animate_qpoint}", f_out)
            print_dual(f"Band index        : {args.animate_mode} ({freq:.4f} THz)", f_out)
            print_dual(f"Frames            : {len(frames)} ({args.animate_format} format)", f_out)
            print_dual(f"[OK] Saved {anim_path}", f_out)

    print_section("[4] DENSITY OF STATES", f_out)
    series = [(("Fine-tuned" if args.custom_model else f"MACE-MP-0 ({args.model})"), result_main, _SERIES_COLORS[0])]
    if result_found is not None:
        series.append((f"Foundation ({args.model})", result_found, _SERIES_COLORS[1]))
    dos_plot = os.path.join(args.output_dir, "dos.png")
    plot_dos(series, dos_plot)
    print_dual(f"[OK] Saved {dos_plot}", f_out)
    if result_main["species_pdos"]:
        print_dual(f"Species (PDOS)    : {', '.join(sorted(result_main['species_pdos']))}", f_out)
    if args.save_data:
        dos_data = os.path.join(args.output_dir, "dos.dat")
        cols = [result_main["dos_freqs"], result_main["dos"]]
        header = "freq(THz) total_DOS"
        for sym, pdos in sorted(result_main["species_pdos"].items()):
            cols.append(pdos)
            header += f" {sym}_PDOS"
        np.savetxt(dos_data, np.column_stack(cols), header=header)
        print_dual(f"[OK] Saved {dos_data}", f_out)

    print_section("[5] THERMAL PROPERTIES", f_out)
    thermal_plot = os.path.join(args.output_dir, "thermal.png")
    plot_thermal(series, thermal_plot)
    print_dual(f"[OK] Saved {thermal_plot}", f_out)
    if args.save_data:
        thermal_data = os.path.join(args.output_dir, "thermal.dat")
        np.savetxt(thermal_data, np.column_stack([
            result_main["temps"], result_main["free_energy"], result_main["entropy"], result_main["heat_capacity"]]),
                  header="T(K) FreeEnergy(kJ/mol) Entropy(J/K/mol) HeatCapacity(J/K/mol)")
        print_dual(f"[OK] Saved {thermal_data}", f_out)

    print_section("[6] SUMMARY & FILES", f_out)
    print_dual("Status            : OK", f_out)
    print_dual(f"Force constants   : {yaml_path}", f_out)
    if bands_plot:
        print_dual(f"Band structure    : {bands_plot}", f_out)
    print_dual(f"DOS               : {dos_plot}", f_out)
    print_dual(f"Thermal properties: {thermal_plot}", f_out)
    if report_path:
        print_dual(f"Report            : {report_path}", f_out)

    if f_out:
        f_out.close()

    print("\n" + "-" * 60)
    print(color_text("ML phonons complete -- a fast heuristic, not a substitute for DFT phonons.\n", 'bold'))


if __name__ == "__main__":
    main()
