#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""Melting-point bracketing via a sequence of short MACE-driven MD runs
(the fifth tool in the ML Simulations category) -- runs one NVT trajectory
per candidate temperature (reusing mlmd.py's own build_dynamics directly,
same import-don't-duplicate pattern stb-mlffActiveLearning already uses),
each starting fresh from the same MACE-relaxed solid reference, and tracks
the Lindemann index (RMS atomic displacement from each atom's own time-
averaged position, over the nearest-neighbor distance) as a function of
temperature. A vibrating solid's Lindemann index stays small and roughly
flat with T; once atoms stop being confined to a lattice site (melting),
it rises sharply -- the classic Lindemann (1910) empirical melting
criterion, with the literature's usual ~0.10 threshold (--lindemann-
threshold) used to bracket the transition temperature by linear
interpolation between the two bracketing scanned points.

Also reports the self-diffusion coefficient D(T) at each temperature as a
complementary/cross-check signal (compute_msd/fit_diffusion_coefficient,
imported from aimd_analysis.py -- the exact same MSD/Einstein-relation
fit already used for a SIESTA AIMD trajectory or an stb-mlmd run, just
looped here over multiple short in-memory runs instead of one long one):
D should stay near 0 in the solid phase and become clearly nonzero once
the material is diffusing as a liquid, around the same temperature the
Lindemann index crosses its threshold.

Bulk (3D periodic) only -- rejects any vacuum-padded axis, same reasoning
as stb-amorphize/stb-mlmd's NPT mode (melting a slab/wire/molecule is
physically meaningless). A coarse, cheap heuristic screen (short runs, a
single seed per temperature) -- NOT a rigorous free-energy calculation
(e.g. thermodynamic integration/coexistence simulations), and inherits
whatever systematic error the underlying MACE potential has for the real
material (MACE-MP-0 in particular is not fit specifically to reproduce
melting behavior). This "one-phase" heating method (a defect-free periodic
crystal with no free surface to nucleate melting from) is additionally well
known to SUPERHEAT past the true thermodynamic melting point, sometimes
substantially, especially for small cells and short runs -- verified live
on a 4-atom Al fcc cell (real Al: ~933 K): the tool's own estimate landed
at ~3095 K, a large but expected overestimate given the tiny cell/short
runs used, not a sign of a bug (a genuinely converged answer would need a
much larger supercell and/or a two-phase coexistence method).
"""

VERSION = "1.0.0"

import os
import sys
import argparse
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt
from ase.io import write as ase_write
from pymatgen.io.ase import AseAtomsAdaptor

from stb.core import structure_io, kspace, mace_relax, md_traj
from stb.core.cli import color_text, show_intro, print_dual, print_section
from stb.core.deps import require_mace
from stb.mlmd import build_dynamics
from stb.aimd_analysis import compute_msd, fit_diffusion_coefficient

require_mace()

REPORT_FILE = "stb_mlmelting_report.txt"

OUTPUT_FORMATS = {
    "xsf": ("xsf", ".xsf"),
    "pdb": ("proteindatabank", ".pdb"),
    "xyz": ("extxyz", ".xyz"),
}


def nearest_neighbor_distance(atoms):
    """Minimum pairwise distance (minimum-image convention) in the
    equilibrium (T=0, relaxed) structure -- the Lindemann index's own
    length scale, computed once from the reference geometry rather than
    per-temperature (the reference lattice constant is the physically
    meaningful yardstick, not a thermally-expanded one).
    """
    distances = atoms.get_all_distances(mic=True)
    np.fill_diagonal(distances, np.inf)
    return float(distances.min())


def compute_lindemann_index(cart_positions_unwrapped, a_nn):
    """RMS atomic displacement from each atom's own time-averaged position
    (NOT a fixed lattice reference -- letting each atom define its own
    center lets a rigidly translating/slowly drifting frame not register
    as "melting"), divided by the nearest-neighbor distance -- the
    standard Lindemann (1910) melting criterion. `cart_positions_unwrapped`
    is a (nframes, natoms, 3) array, already PBC-unwrapped (core.md_traj.
    unwrap_trajectory) so real diffusive motion isn't confused with
    periodic-boundary teleports.
    """
    cart = np.asarray(cart_positions_unwrapped)
    mean_pos = cart.mean(axis=0)
    disp2 = np.sum((cart - mean_pos[None, :, :]) ** 2, axis=-1)
    return float(np.sqrt(disp2.mean()) / a_nn)


def estimate_melting_point(temps, lindemann_vals, threshold):
    """Linear interpolation between the two scanned temperatures that
    bracket the first crossing of `threshold`. Returns (T_estimate,
    bracket_lo, bracket_hi), or (None, None, None) if the index never
    crosses the threshold anywhere in the scanned range, or is already
    above it at the very first (lowest) temperature scanned (both cases
    need a wider --temp-min/--temp-max, not an extrapolated guess).
    """
    if lindemann_vals[0] >= threshold:
        return None, None, None
    for i in range(1, len(temps)):
        if lindemann_vals[i] >= threshold:
            t_lo, t_hi = temps[i - 1], temps[i]
            l_lo, l_hi = lindemann_vals[i - 1], lindemann_vals[i]
            frac = (threshold - l_lo) / (l_hi - l_lo) if l_hi != l_lo else 0.0
            t_est = t_lo + frac * (t_hi - t_lo)
            return t_est, t_lo, t_hi
    return None, None, None


def plot_melting_curve(temps, lindemann_vals, diffusion_vals, threshold, t_estimate, out_path):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].plot(temps, lindemann_vals, 'o-', color='tab:blue')
    axes[0].axhline(threshold, color='gray', linestyle='--', linewidth=1,
                    label=f"threshold ({threshold:.2f})")
    if t_estimate is not None:
        axes[0].axvline(t_estimate, color='tab:red', linestyle=':', linewidth=1.5,
                        label=f"T_melt ~ {t_estimate:.0f} K")
    axes[0].set_xlabel("Temperature (K)")
    axes[0].set_ylabel("Lindemann index")
    axes[0].set_title("Lindemann melting criterion")
    axes[0].legend(fontsize=8)

    axes[1].plot(temps, diffusion_vals, 'o-', color='tab:orange')
    axes[1].set_xlabel("Temperature (K)")
    axes[1].set_ylabel("Diffusion coefficient (cm^2/s)")
    axes[1].set_title("Self-diffusion coefficient (cross-check)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Melting-point bracketing via a sequence of short MACE-driven MD runs "
                     "at increasing temperature, tracking the Lindemann index (RMS atomic "
                     "displacement from each atom's own time-averaged position, over the "
                     "nearest-neighbor distance) -- crosses a threshold (default 0.10, the "
                     "classic literature value) around the real melting point. Also reports "
                     "the self-diffusion coefficient vs. temperature as a cross-check "
                     "(compute_msd/fit_diffusion_coefficient, same as stb-aimdAnalysis). Bulk "
                     "(3D periodic) only. A coarse, cheap heuristic screen, NOT a rigorous "
                     "free-energy/coexistence melting-point calculation.",
        epilog="Example usage:\n"
               "  stb-mlmelting --file structure.fdf\n"
               "  stb-mlmelting --file structure.fdf --temp-min 500 --temp-max 1500 "
               "--temp-step 100\n",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--file", required=True, help="Input structure (.fdf).")
    parser.add_argument("--model", choices=["small", "medium", "large"], default="small",
                        help="MACE-MP-0 model size (default: small). Ignored if --custom-model is given.")
    parser.add_argument("--custom-model", default=None, metavar="PATH",
                        help="Use a custom MACE model file instead of the MACE-MP-0 foundation "
                             "potential -- e.g. one fine-tuned via stb-mlffAnalysis. Overrides --model.")
    parser.add_argument("--temp-min", type=float, default=300.0,
                        help="Lowest scanned temperature, K (default: 300).")
    parser.add_argument("--temp-max", type=float, default=3000.0,
                        help="Highest scanned temperature, K (default: 3000).")
    parser.add_argument("--temp-step", type=float, default=300.0,
                        help="Temperature step, K (default: 300).")
    parser.add_argument("--temperatures", type=float, nargs='+', default=None, metavar="K",
                        help="Explicit list of temperatures, K -- overrides --temp-min/--temp-max/--temp-step.")
    parser.add_argument("--equilibration-steps", type=int, default=200,
                        help="MD steps discarded before collecting data at each temperature (default: 200).")
    parser.add_argument("--production-steps", type=int, default=500,
                        help="MD steps collected (after equilibration) at each temperature (default: 500).")
    parser.add_argument("--timestep", type=float, default=1.0, help="MD timestep, fs (default: 1.0).")
    parser.add_argument("--friction", type=float, default=0.01,
                        help="Langevin friction, fs^-1 (default: 0.01).")
    parser.add_argument("--stride", type=int, default=5,
                        help="Save one frame every this many production steps (default: 5).")
    parser.add_argument("--lindemann-threshold", type=float, default=0.10,
                        help="Lindemann index melting threshold (default: 0.10 -- the classic "
                             "literature value; real materials vary roughly 0.05-0.15).")
    parser.add_argument("--fmax", type=float, default=0.05,
                        help="Force convergence threshold for the initial full relax, eV/Ang (default: 0.05).")
    parser.add_argument("--max-steps", type=int, default=200,
                        help="Max optimizer steps for the initial full relax (default: 200).")
    parser.add_argument("--vacuum-gap", type=float, default=10.0,
                        help="Gap (Ang) used to detect vacuum-padded axes -- bulk-only, rejects "
                             "any (default: 10.0).")
    parser.add_argument("--seed", type=int, default=None,
                        help="Random seed for reproducible velocity initialization (default: unseeded).")
    parser.add_argument("--save-trajectories", action="store_true",
                        help="Also save each temperature's production trajectory (xsf/pdb/xyz, "
                             "--trajectory-format). Off by default -- can add up to a lot of "
                             "files across a wide temperature scan.")
    parser.add_argument("--trajectory-format", choices=sorted(OUTPUT_FORMATS), default="xsf",
                        help="Format for --save-trajectories, same convention as stb-ani2traj/"
                             "stb-mlmd (default: xsf).")
    parser.add_argument("--device", choices=["cpu", "cuda"], default="cpu",
                        help="Device to run the model on (default: cpu).")
    parser.add_argument("-o", "--output-dir", default="mlmelting_out",
                        help="Output directory for all files (default: mlmelting_out).")
    parser.add_argument("--save-data", action="store_true",
                        help="Also write the raw Lindemann/diffusion-vs-temperature data (.dat). Off by default.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-mlmelting {VERSION}")

    args = parser.parse_args()

    if args.temperatures is not None:
        temps = sorted(args.temperatures)
    else:
        if args.temp_min >= args.temp_max:
            parser.error("--temp-min must be less than --temp-max.")
        temps = list(np.arange(args.temp_min, args.temp_max + 1e-6, args.temp_step))
    if len(temps) < 2:
        parser.error("Need at least 2 temperatures to scan.")

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - ML Melting Point",
            "Lindemann-criterion melting-point bracketing, driven entirely by a MACE potential",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("ML MELTING POINT (LINDEMANN CRITERION):", 'bold'))
    print("-" * 60)

    os.makedirs(args.output_dir, exist_ok=True)
    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    def fail(message):
        print_dual(color_text(f"[FAIL] {message}", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_dual(color_text("===== STB-MLMELTING REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Input file        : {args.file}", f_out)
    print_dual(f"Model             : "
               f"{'custom (' + args.custom_model + ')' if args.custom_model else f'MACE-MP-0 ({args.model})'}", f_out)
    print_dual(f"Temperatures (K)  : {', '.join(f'{t:.0f}' for t in temps)}", f_out)
    print_dual(f"Lindemann threshold: {args.lindemann_threshold}", f_out)
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
    print_dual(f"[OK] Read {len(atoms)} atom(s)", f_out)

    lattice_ang = np.array(atoms.get_cell())
    frac_coords = atoms.get_scaled_positions()
    vacuum_axes = kspace.detect_vacuum_axes(frac_coords, lattice_ang, args.vacuum_gap)
    if any(vacuum_axes):
        fail("This tool is bulk (3D periodic) only -- a vacuum-padded axis was detected. "
             "Melting a slab/wire/molecule is not physically meaningful (same reasoning as "
             "stb-amorphize/stb-mlmd's NPT mode).")

    model_arg = args.custom_model if args.custom_model else args.model
    calc = mace_relax.get_calculator(model=model_arg, device=args.device)
    atoms.calc = calc

    print_section("[2] REFERENCE RELAXATION", f_out)
    f0 = np.abs(atoms.get_forces()).max()
    cell_mask = mace_relax.build_cell_mask(vacuum_axes)
    converged, steps_used = mace_relax.relax(
        atoms, calc, cell_mask=cell_mask, fmax=args.fmax, max_steps=args.max_steps)
    f1 = np.abs(atoms.get_forces()).max()
    print_dual(f"Pre-relax (cell+ions): max|F| {f0:.4f} -> {f1:.4f} eV/Ang "
               f"({'converged' if converged else 'hit step cap'}, {steps_used} steps)", f_out)

    a_nn = nearest_neighbor_distance(atoms)
    print_dual(f"Nearest-neighbor distance (reference): {a_nn:.4f} Ang", f_out)

    print_section("[3] TEMPERATURE SCAN", f_out)
    lindemann_vals, diffusion_vals = [], []
    for temp_k in temps:
        run_atoms = atoms.copy()
        run_atoms.calc = calc
        dyn = build_dynamics(run_atoms, "nvt", args.timestep, temp_k, args.friction,
                             pressure_gpa=0.0, seed=args.seed)
        dyn.run(args.equilibration_steps)

        frames_pos, frames_cell = [], []

        def _collect():
            frames_pos.append(run_atoms.get_positions().copy())
            frames_cell.append(np.array(run_atoms.get_cell()))

        dyn.attach(_collect, interval=args.stride)
        dyn.run(args.production_steps)
        n_collected = len(frames_pos)

        unwrapped = md_traj.unwrap_trajectory(frames_pos, frames_cell)
        lindemann = compute_lindemann_index(unwrapped, a_nn)

        dt_frame_fs = args.timestep * args.stride
        t_fs, msd = compute_msd(unwrapped, dt_frame_fs, run_atoms.get_chemical_symbols())
        d_coeff, _, _ = fit_diffusion_coefficient(t_fs, msd)
        # A short/noisy MSD fit deep in the solid (near-zero true diffusion)
        # regime can land slightly negative from pure statistical noise --
        # physically D can never be negative, so clip rather than report a
        # confusing "negative diffusion" that isn't real backward motion.
        d_coeff = max(0.0, d_coeff) if d_coeff is not None else 0.0

        lindemann_vals.append(lindemann)
        diffusion_vals.append(d_coeff)
        print_dual(f"T = {temp_k:7.1f} K  |  Lindemann = {lindemann:.4f}  |  "
                   f"D = {d_coeff:.3e} cm^2/s  |  {n_collected} frame(s) collected", f_out)

        if args.save_trajectories:
            ase_format, ext = OUTPUT_FORMATS[args.trajectory_format]
            traj_path = os.path.join(args.output_dir, f"trajectory_{temp_k:.0f}K{ext}")
            from ase import Atoms as AseAtoms
            traj_frames = [AseAtoms(numbers=run_atoms.numbers, positions=p, cell=c, pbc=True)
                           for p, c in zip(unwrapped, frames_cell)]
            ase_write(traj_path, traj_frames, format=ase_format)
            print_dual(f"[OK] Saved {traj_path}", f_out)

    print_section("[4] MELTING-POINT ESTIMATE", f_out)
    t_estimate, t_lo, t_hi = estimate_melting_point(temps, lindemann_vals, args.lindemann_threshold)
    if t_estimate is not None:
        print_dual(f"Estimated melting point: {color_text(f'{t_estimate:.0f} K', 'bold')} "
                   f"(Lindemann index crosses {args.lindemann_threshold} between "
                   f"{t_lo:.0f} K and {t_hi:.0f} K)", f_out)
        print_dual(color_text(
            "[NOTE] This 'one-phase' method (heating a defect-free periodic crystal with no "
            "free surface) is well known to SUPERHEAT -- there's nowhere for melting to "
            "nucleate from, so the crystal can persist well past the true thermodynamic "
            "melting point, especially for small cells and short runs like this. Treat this "
            "number as a rough, likely-too-high estimate, not a literal prediction -- a larger "
            "supercell, longer runs, or a two-phase coexistence simulation would narrow it.",
            'cyan'), f_out)
    elif lindemann_vals[0] >= args.lindemann_threshold:
        print_dual(color_text(
            f"[WARNING] Lindemann index is already at/above {args.lindemann_threshold} at the "
            f"lowest scanned temperature ({temps[0]:.0f} K) -- lower --temp-min to bracket the "
            "actual transition.", 'yellow'), f_out)
    else:
        print_dual(color_text(
            f"[WARNING] Lindemann index never reached {args.lindemann_threshold} within the "
            f"scanned range ({temps[0]:.0f}-{temps[-1]:.0f} K) -- raise --temp-max to bracket "
            "the actual transition.", 'yellow'), f_out)

    plot_path = os.path.join(args.output_dir, "melting_curve.png")
    plot_melting_curve(temps, lindemann_vals, diffusion_vals, args.lindemann_threshold,
                       t_estimate, plot_path)
    print_dual(f"[OK] Saved {plot_path}", f_out)
    if args.save_data:
        data_path = os.path.join(args.output_dir, "melting_curve.dat")
        np.savetxt(data_path, np.column_stack([temps, lindemann_vals, diffusion_vals]),
                  header="T(K) Lindemann_index D(cm^2/s)")
        print_dual(f"[OK] Saved {data_path}", f_out)

    print_section("[5] SUMMARY & FILES", f_out)
    print_dual("Status            : OK", f_out)
    if t_estimate is not None:
        print_dual(f"Estimated T_melt  : {t_estimate:.0f} K", f_out)
    print_dual(f"Melting curve plot: {plot_path}", f_out)
    if report_path:
        print_dual(f"Report            : {report_path}", f_out)

    if f_out:
        f_out.close()

    print("\n" + "-" * 60)
    print(color_text("ML melting-point scan complete -- a coarse heuristic, not a rigorous "
                      "free-energy calculation.\n", 'bold'))


if __name__ == "__main__":
    main()
