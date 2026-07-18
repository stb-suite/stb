#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.0.0"

import os
import sys
import argparse
from collections import Counter
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt

from stb.core.deps import require_sisl
sisl = require_sisl()

from stb.core.md_traj import read_frame_lattices, read_md_timestep_fs, unwrap_trajectory
from stb.core.cli import color_text, show_intro, print_dual, print_section

REPORT_FILE = "stb_aimdAnalysis_report.txt"


def cell_min_perpendicular_width(cell):
    """Half the shortest perpendicular distance between opposite faces of
    the cell (V / |cross product of the other two vectors|, minimized over
    the 3 axis pairs) -- the correct upper bound for a minimum-image RDF cutoff
    (unlike just taking half the shortest lattice VECTOR norm, this stays
    correct for a non-orthogonal/skewed cell too)."""
    a, b, c = cell
    vol = abs(np.dot(a, np.cross(b, c)))
    widths = [
        vol / np.linalg.norm(np.cross(b, c)),
        vol / np.linalg.norm(np.cross(c, a)),
        vol / np.linalg.norm(np.cross(a, b)),
    ]
    return min(widths) / 2.0


def compute_rdf(frac_positions, cells, symbols, pair=None, r_max=None, n_bins=200):
    """Average radial distribution function g(r) over all frames.

    `frac_positions` is a list (one per frame) of (natoms, 3) FRACTIONAL
    coordinate arrays -- distances are computed in fractional space (minimum-
    image convention, `disp -= round(disp)`) and converted to Cartesian via
    each frame's own cell, so this stays correct even for a variable-cell
    trajectory. Vectorized over all pairs at once per frame (numpy
    broadcasting) rather than building an ase.Atoms per frame, since a long
    AIMD trajectory can have thousands of frames.

    `pair`, if given, is a (species_a, species_b) tuple restricting the
    histogram to that pair only (species_a == species_b is a valid same-
    species RDF); None means every pair (the total RDF).

    Returns (r_centers, g_r). r_max defaults to the smallest
    cell_min_perpendicular_width found across all frames (the physical limit
    of the minimum-image convention).
    """
    natoms = frac_positions[0].shape[0]
    if r_max is None:
        r_max = min(cell_min_perpendicular_width(c) for c in cells)

    if pair is not None:
        idx_a = [i for i, s in enumerate(symbols) if s == pair[0]]
        idx_b = [i for i, s in enumerate(symbols) if s == pair[1]]
    else:
        idx_a = idx_b = list(range(natoms))

    bin_edges = np.linspace(0.0, r_max, n_bins + 1)
    counts = np.zeros(n_bins)
    n_frames = len(frac_positions)
    density_sum = 0.0

    for frac, cell in zip(frac_positions, cells):
        vol = abs(np.dot(cell[0], np.cross(cell[1], cell[2])))
        pa = frac[idx_a]
        pb = frac[idx_b]
        disp = pa[:, None, :] - pb[None, :, :]
        disp -= np.round(disp)
        cart = disp @ cell
        dist = np.linalg.norm(cart, axis=-1)
        same_species_block = pair is None or pair[0] == pair[1]
        if same_species_block and idx_a is idx_b:
            np.fill_diagonal(dist, np.inf)
        dist = dist[(dist > 1e-6) & (dist < r_max)]
        h, _ = np.histogram(dist, bins=bin_edges)
        counts += h
        density_sum += len(idx_b) / vol

    r_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    dr = bin_edges[1] - bin_edges[0]
    mean_density = density_sum / n_frames
    shell_vol = 4.0 * np.pi * r_centers**2 * dr
    norm = len(idx_a) * n_frames * shell_vol * mean_density
    norm[norm == 0] = np.nan
    g_r = counts / norm
    return r_centers, g_r


def compute_msd(cart_positions, dt_fs, symbols, species=None):
    """Mean-squared displacement vs time, averaged over multiple time
    origins (a sliding window over the whole trajectory, not just t=0 --
    standard MD practice for reducing statistical noise) and over every
    atom of `species` (or all atoms if None).

    `cart_positions` must already be PBC-unwrapped (core.md_traj.
    unwrap_trajectory) -- an MSD computed on wrapped/re-folded coordinates
    would be dominated by box-edge teleports instead of real diffusion.

    Returns (t_fs, msd) with t_fs starting at 0.
    """
    idx = ([i for i, s in enumerate(symbols) if s == species]
           if species is not None else list(range(len(symbols))))
    traj = np.array([p[idx] for p in cart_positions])  # (nframes, nsel, 3)
    n_frames = traj.shape[0]

    msd = np.zeros(n_frames)
    for lag in range(n_frames):
        disp = traj[lag:] - traj[:n_frames - lag]
        msd[lag] = np.mean(np.sum(disp**2, axis=-1))

    t_fs = np.arange(n_frames) * dt_fs
    return t_fs, msd


def fit_diffusion_coefficient(t_fs, msd, fit_start=None, fit_end=None):
    """Linear fit of MSD(t) over [fit_start, fit_end] fs (default: the second
    half of the trajectory, where the regime is diffusive rather than the
    short-time ballistic ramp) via the 3D Einstein relation D = slope / 6.
    Returns (D_cm2_s, slope, (t_lo, t_hi)). slope is in Ang^2/fs; the
    Ang^2/fs -> cm^2/s conversion factor is 1e-16 (Ang^2 -> cm^2 = 1e-16) /
    1e-15 (fs -> s) = 1e-1... done explicitly below rather than folded into
    a magic constant.
    """
    t_hi = t_fs[-1]
    t_lo = fit_start if fit_start is not None else t_hi / 2.0
    t_hi = fit_end if fit_end is not None else t_hi
    mask = (t_fs >= t_lo) & (t_fs <= t_hi)
    if mask.sum() < 2:
        return None, None, (t_lo, t_hi)
    slope, _ = np.polyfit(t_fs[mask], msd[mask], 1)  # Ang^2 / fs
    ang2_to_cm2 = 1.0e-16
    fs_to_s = 1.0e-15
    D_cm2_s = (slope / 6.0) * ang2_to_cm2 / fs_to_s
    return D_cm2_s, slope, (t_lo, t_hi)


def compute_vacf_vdos(cart_positions, dt_fs):
    """Velocity autocorrelation function (VACF) and the vibrational density
    of states (VDOS) derived from it.

    .ANI carries no velocities, only positions -- v(t) is estimated by
    central finite difference, v(t) = (x(t+dt) - x(t-dt)) / (2*dt), the
    standard approximation when only a position trajectory is available.
    This is markedly less accurate than SIESTA's own MD velocities would be
    (it amplifies any positional noise into the derivative) -- treat the
    resulting VDOS as qualitative, not benchmark-grade.

    VACF is normalized (VACF(0) = 1), averaged over atoms and over sliding
    time origins (via FFT-based autocorrelation, O(N log N) instead of the
    O(N^2) direct sum -- matters for long trajectories). VDOS is the real
    part of the FFT of the VACF (Wiener-Khinchin theorem), with the
    frequency axis converted from 1/fs to cm^-1 (1 fs^-1 = 1e15 Hz;
    freq_cm-1 = freq_Hz / (c * 100), c in m/s).

    Returns (freq_cm1, vdos, t_fs, vacf) -- freq_cm1/vdos only cover the
    non-negative-frequency half of the spectrum.
    """
    pos = np.asarray(cart_positions)  # (nframes, natoms, 3)
    vel = np.gradient(pos, dt_fs, axis=0)
    n_frames, natoms, _ = vel.shape

    vacf = np.zeros(n_frames)
    for a in range(natoms):
        for d in range(3):
            v = vel[:, a, d]
            v = v - v.mean()
            n = len(v)
            padded = np.zeros(2 * n)
            padded[:n] = v
            f = np.fft.fft(padded)
            acf = np.fft.ifft(f * np.conj(f)).real[:n]
            acf /= np.arange(n, 0, -1)
            vacf += acf
    vacf /= (natoms * 3)
    if vacf[0] != 0:
        vacf /= vacf[0]

    t_fs = np.arange(n_frames) * dt_fs

    spectrum = np.fft.rfft(vacf)
    freqs_per_fs = np.fft.rfftfreq(n_frames, d=dt_fs)
    c_cm_per_s = 2.99792458e10
    freq_cm1 = (freqs_per_fs * 1.0e15) / c_cm_per_s
    vdos = np.abs(spectrum.real)
    return freq_cm1, vdos, t_fs, vacf


def main():
    parser = argparse.ArgumentParser(
        description="Physical post-processing of a SIESTA AIMD trajectory (<label>.ANI): "
                     "radial distribution function g(r), mean-squared displacement (MSD) "
                     "and diffusion coefficient, and a VACF-derived vibrational density of "
                     "states (VDOS). Complements stb-ani2traj, which only converts the "
                     "trajectory's format for external viewers.",
        epilog="Example usage:\n"
               "  stb-aimdAnalysis --label aimd\n"
               "  stb-aimdAnalysis -l aimd --pair O-H --skip 50\n"
               "  stb-aimdAnalysis -l aimd --fit-start 5000 --fit-end 15000",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("-l", "--label", required=True, help="SystemLabel (required).")
    parser.add_argument("--stride", type=int, default=1,
                        help="Keep only every Nth frame (default: 1, every frame).")
    parser.add_argument("--skip", type=int, default=0,
                        help="Discard the first N frames (post-stride) as thermal "
                             "equilibration, before any statistics are computed (default: 0).")
    parser.add_argument("--pair", default=None, metavar="A-B",
                        help="Restrict the RDF to one species pair, e.g. 'O-H' (default: "
                             "total RDF over all pairs).")
    parser.add_argument("--r-max", type=float, default=None,
                        help="RDF cutoff in Angstrom (default: the minimum-image limit of "
                             "the smallest cell found in the trajectory).")
    parser.add_argument("--n-bins", type=int, default=200, help="Number of RDF bins (default: 200).")
    parser.add_argument("--fit-start", type=float, default=None, metavar="FS",
                        help="MSD linear-fit window start, in fs (default: half the "
                             "trajectory length).")
    parser.add_argument("--fit-end", type=float, default=None, metavar="FS",
                        help="MSD linear-fit window end, in fs (default: end of trajectory).")
    parser.add_argument("--save-data", action="store_true",
                        help="Also write the raw numeric data (.dat) behind each plot. Off by default.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-aimdAnalysis {VERSION}")

    args = parser.parse_args()

    pair = None
    if args.pair:
        parts = args.pair.split("-")
        if len(parts) != 2:
            parser.error(f"--pair must be 'A-B' (e.g. 'O-H'), got '{args.pair}'.")
        pair = tuple(parts)

    ani_file = f"{args.label}.ANI"

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - AIMD Trajectory Analysis",
            "RDF, MSD/diffusion, and VACF-derived VDOS",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("AIMD TRAJECTORY ANALYSIS:", 'bold'))
    print("-" * 60)

    report_path = REPORT_FILE if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    def fail(message):
        print_dual(color_text(f"[FAIL] {message}", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_dual(color_text("===== STB-AIMDANALYSIS REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Input file        : {ani_file}", f_out)
    print_dual(f"Stride / skip     : {args.stride} / {args.skip}", f_out)
    if pair:
        print_dual(f"RDF pair          : {pair[0]}-{pair[1]}", f_out)
    if report_path:
        print_dual(f"Report file       : {report_path}", f_out)

    print_section("[1] READING TRAJECTORY", f_out)

    if not os.path.exists(ani_file):
        fail(f"Input file '{ani_file}' not found.")
    if args.stride < 1:
        fail(f"--stride must be >= 1, got {args.stride}.")

    try:
        sile = sisl.get_sile(ani_file)
        all_frames = sile.read_geometry[:]()
    except Exception as e:
        fail(f"Could not read '{ani_file}': {e}")

    if not all_frames:
        fail(f"'{ani_file}' contains no frames.")

    all_cells, all_steps = read_frame_lattices(args.label, len(all_frames), f_out)
    if all_cells is None:
        fail(f"Could not find a readable '{args.label}.out', '{args.label}.XV' or "
             f"'{args.label}.fdf' -- a lattice is required ('{ani_file}' doesn't carry "
             f"one of its own).")

    frames = all_frames[::args.stride]
    cells = all_cells[::args.stride]

    if args.skip >= len(frames):
        fail(f"--skip {args.skip} discards all {len(frames)} available (post-stride) frame(s).")
    frames = frames[args.skip:]
    cells = cells[args.skip:]

    symbols = [a.symbol for a in frames[0].atoms]
    if pair and (pair[0] not in symbols or pair[1] not in symbols):
        fail(f"--pair {args.pair}: species not found in structure "
             f"(available: {', '.join(sorted(set(symbols)))}).")

    print_dual(f"[OK] Read {len(frames)} frame(s) from {ani_file} "
               f"(stride {args.stride}, skip {args.skip}).", f_out)
    print_dual(f"Atoms per frame   : {len(symbols)}", f_out)
    composition = ", ".join(f"{el}{n}" for el, n in sorted(Counter(symbols).items()))
    print_dual(f"Composition       : {composition}", f_out)

    initial_step, dt_fs = read_md_timestep_fs(args.label)
    if dt_fs is None:
        print_dual(color_text(
            "[WARNING] MD.LengthTimeStep not found in "
            f"'{args.label}.fdf' -- assuming 1 fs per (strided) frame for MSD/VDOS time axes.",
            'yellow'), f_out)
        dt_fs = 1.0
    else:
        dt_fs *= args.stride

    frac_positions = [g.fxyz for g in frames]
    cart_positions_raw = [g.xyz for g in frames]
    cart_positions = unwrap_trajectory(cart_positions_raw, cells)

    print_section("[2] RADIAL DISTRIBUTION FUNCTION (RDF)", f_out)
    r, g_r = compute_rdf(frac_positions, cells, symbols, pair=pair,
                         r_max=args.r_max, n_bins=args.n_bins)
    valid = ~np.isnan(g_r)
    if valid.any():
        peak_idx = np.nanargmax(g_r)
        print_dual(f"[OK] RDF computed over {len(frames)} frame(s), r up to {r[-1]:.3f} Ang.", f_out)
        print_dual(f"First peak        : r = {r[peak_idx]:.3f} Ang (g(r) = {g_r[peak_idx]:.3f})", f_out)
    else:
        print_dual(color_text("[WARNING] RDF could not be computed (no pairs in range).", 'yellow'), f_out)

    pair_tag = f"_{pair[0]}{pair[1]}" if pair else ""
    rdf_plot = f"{args.label}_rdf{pair_tag}.png"
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(r, g_r)
    ax.set_xlabel("r (Ang)")
    ax.set_ylabel("g(r)")
    ax.set_title(f"RDF{' (' + pair[0] + '-' + pair[1] + ')' if pair else ''} -- {args.label}")
    ax.axhline(1.0, color='gray', linestyle='--', linewidth=0.8)
    fig.tight_layout()
    fig.savefig(rdf_plot, dpi=150)
    plt.close(fig)
    print_dual(f"[OK] Saved {rdf_plot}", f_out)
    rdf_data_file = None
    if args.save_data:
        rdf_data_file = f"{args.label}_rdf{pair_tag}.dat"
        np.savetxt(rdf_data_file, np.column_stack([r, g_r]), header="r(Ang) g(r)")
        print_dual(f"[OK] Saved {rdf_data_file}", f_out)

    print_section("[3] MEAN-SQUARED DISPLACEMENT (MSD) & DIFFUSION", f_out)
    species_set = sorted(set(symbols))
    msd_results = {}
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for sp in species_set + [None]:
        t_fs, msd = compute_msd(cart_positions, dt_fs, symbols, species=sp)
        D, slope, (t_lo, t_hi) = fit_diffusion_coefficient(t_fs, msd, args.fit_start, args.fit_end)
        label = sp if sp is not None else "all"
        msd_results[label] = (t_fs, msd, D)
        ax.plot(t_fs, msd, label=label)
        if D is not None:
            print_dual(f"[OK] Diffusion coefficient ({label:<4}): D = {D:.3e} cm^2/s "
                       f"(fit window: {t_lo:.1f}-{t_hi:.1f} fs)", f_out)
        else:
            print_dual(color_text(
                f"[WARNING] Diffusion coefficient ({label}): not enough frames in the "
                f"fit window to fit.", 'yellow'), f_out)
    ax.set_xlabel("t (fs)")
    ax.set_ylabel("MSD (Ang^2)")
    ax.set_title(f"MSD -- {args.label}")
    ax.legend()
    fig.tight_layout()
    msd_plot = f"{args.label}_msd.png"
    fig.savefig(msd_plot, dpi=150)
    plt.close(fig)
    print_dual(f"[OK] Saved {msd_plot}", f_out)
    msd_data_file = None
    if args.save_data:
        msd_data_file = f"{args.label}_msd.dat"
        t_all, msd_all, _ = msd_results["all"]
        np.savetxt(msd_data_file, np.column_stack([t_all, msd_all]), header="t(fs) MSD_all(Ang^2)")
        print_dual(f"[OK] Saved {msd_data_file}", f_out)

    print_section("[4] VELOCITY AUTOCORRELATION / VIBRATIONAL DOS (VDOS)", f_out)
    if len(frames) < 4:
        print_dual(color_text(
            "[WARNING] Trajectory too short for a meaningful VACF/VDOS "
            f"(only {len(frames)} frame(s)) -- skipping.", 'yellow'), f_out)
        vdos_plot = None
        vdos_data_file = None
    else:
        freq_cm1, vdos, t_fs_vacf, vacf = compute_vacf_vdos(cart_positions, dt_fs)
        peak_idx = np.argmax(vdos[1:]) + 1 if len(vdos) > 1 else 0
        print_dual(f"[OK] VDOS computed (finite-difference velocities -- qualitative, "
                   f"not benchmark-grade).", f_out)
        print_dual(f"Dominant peak     : {freq_cm1[peak_idx]:.1f} cm^-1", f_out)

        fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
        axes[0].plot(t_fs_vacf, vacf)
        axes[0].set_xlabel("t (fs)")
        axes[0].set_ylabel("VACF (normalized)")
        axes[0].set_title("Velocity autocorrelation")
        axes[1].plot(freq_cm1, vdos)
        axes[1].set_xlabel("Frequency (cm^-1)")
        axes[1].set_ylabel("VDOS (arb. units)")
        axes[1].set_title("Vibrational density of states")
        fig.suptitle(args.label)
        fig.tight_layout()
        vdos_plot = f"{args.label}_vdos.png"
        fig.savefig(vdos_plot, dpi=150)
        plt.close(fig)
        print_dual(f"[OK] Saved {vdos_plot}", f_out)
        vdos_data_file = None
        if args.save_data:
            vdos_data_file = f"{args.label}_vdos.dat"
            np.savetxt(vdos_data_file, np.column_stack([freq_cm1, vdos]),
                      header="freq(cm^-1) VDOS(arb.units)")
            print_dual(f"[OK] Saved {vdos_data_file}", f_out)

    print_section("[5] SUMMARY & FILES", f_out)
    print_dual("Status            : OK", f_out)
    for f in [rdf_plot, rdf_data_file, msd_plot, msd_data_file, vdos_plot, vdos_data_file]:
        if f:
            print_dual(f"  - {f}", f_out)
    if report_path:
        print_dual(f"Report            : {report_path}", f_out)

    if f_out:
        f_out.close()

    print("\n" + "-" * 60)
    print(color_text("Trajectory analyzed -- physics extracted from the noise.\n", 'bold'))


if __name__ == "__main__":
    main()
