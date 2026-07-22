#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""Global structure search driven entirely by a MACE potential -- the
fourth tool in the ML Simulations category, and the first without a direct
analog in the SIESTA-based Workflow tools: unlike stb-mlphonons/
stb-mlelastic (fast ML previews of an existing DFT prep+analysis pair),
this explores configuration space itself (which atomic arrangement, at a
fixed cell, has the lowest energy) -- a genuinely new capability, not a
cheaper version of something the suite already does. Two algorithms
(--algorithm):

- **basin-hopping** (default): wraps ase.optimize.basin.BasinHopping
  directly rather than reimplementing the algorithm -- perturb (uniform
  random atomic displacement, --dr), locally relax with MACE (--fmax),
  accept/reject via the Metropolis criterion (probability exp(-dE/kT),
  --temperature in Kelvin, converted to eV via ase.units.kB same as any
  other ASE-facing temperature in this suite), track the best (lowest
  -energy) locally-relaxed structure found across --steps attempts. ASE's
  own implementation only perturbs atomic POSITIONS, never the cell shape/
  volume (it was originally written for clusters, not crystal structure
  prediction).
- **simulated-annealing**: one continuous NVT MD trajectory (Langevin, via
  mlmd.py's own build_dynamics -- same import-don't-duplicate reuse
  stb-mlffActiveLearning/stb-mlmelting already do) whose target temperature
  is cooled from --temp-start down to --temp-end over --steps total MD
  steps (--schedule exponential, the classic SA convention, or linear),
  snapshotting and locally relaxing the running configuration every
  --snapshot-interval steps to find candidate minima along the cooling
  path -- the same "perturb, relax, keep the best" idea as basin hopping,
  just driven by a physically continuous cooling trajectory instead of
  independent random jumps (can find configurations basin hopping's blind
  random perturbation might not reach, at the cost of the schedule itself
  needing some tuning).

Both are appropriate for this tool's primary use case (finding the best
atomic arrangement at a FIXED cell: defect/vacancy site occupation,
disordered/alloy configurations, local amorphous-like arrangements), not
full ab initio crystal structure prediction -- neither touches the cell
shape/volume during the search itself. An optional final full relax
(--no-final-cell-relax to disable) polishes the winning configuration's
cell shape too, via core/mace_relax.py's build_cell_mask (same vacuum
-aware convention as stb-mlrelax/stb-amorphize).

Explicitly NOT rigorous global optimization (no proof of finding the true
global minimum -- both are heuristic searches) and explicitly NOT crystal
structure prediction (cell shape is only touched by the optional final
polish, never during the search itself).
"""

VERSION = "1.1.0"

import os
import sys
import re
import argparse
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt
from ase import units, Atoms
from ase.optimize import FIRE
from ase.optimize.basin import BasinHopping
from ase.io import read as ase_read, write as ase_write
from pymatgen.io.ase import AseAtomsAdaptor

from stb.core import structure_io, kspace, mace_relax
from stb.core.cli import color_text, show_intro, print_dual, print_section
from stb.core.deps import require_mace
from stb.mlmd import build_dynamics

require_mace()

REPORT_FILE = "stb_mlsearch_report.txt"

# Same 3-format convention already used by stb-ani2traj/stb-mlmd/
# stb-mlphonons for multi-frame trajectories -- duplicated here (not
# imported) for the same reason mlmd.py/mlphonons.py already do: tiny, and
# tied to this tool's own CLI --trajectory-format choices/help text.
OUTPUT_FORMATS = {
    "xsf": ("xsf", ".xsf"),
    "pdb": ("proteindatabank", ".pdb"),
    "xyz": ("extxyz", ".xyz"),
}


def parse_basin_log(log_path):
    """Parses ase.optimize.basin.BasinHopping's own logfile ('BasinHopping:
    step %d, energy %f, emin %f' lines, its fixed built-in format) into
    (steps, energies, running_min) arrays for the diagnostics plot --
    avoids subclassing/monkeypatching BasinHopping just to intercept these
    numbers, since they're already written to disk in a stable, parseable
    format. `energies` is the relaxed energy of each ATTEMPTED candidate
    (accepted or not); `running_min` is the best energy found up to and
    including that step.
    """
    steps, energies, emins = [], [], []
    pattern = re.compile(r"BasinHopping:\s*step\s+(-?\d+),\s*energy\s+([\-\d.]+),\s*emin\s+([\-\d.]+)")
    with open(log_path) as f:
        for line in f:
            m = pattern.search(line)
            if m:
                steps.append(int(m.group(1)))
                energies.append(float(m.group(2)))
                emins.append(float(m.group(3)))
    return np.array(steps), np.array(energies), np.array(emins)


def plot_search_history(steps, energies, emins, out_path):
    """Attempted (relaxed) energy per step (scatter) plus the running
    minimum (line) -- lets a user see at a glance whether the search is
    still finding new minima (running-min line still dropping) or has
    plateaued (flat for many consecutive steps), the same kind of read a
    basin-hopping log is normally eyeballed for.
    """
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(steps, energies, 'o', color='tab:blue', alpha=0.5, markersize=4,
            label='Attempted (relaxed) energy')
    ax.plot(steps, emins, '-', color='tab:orange', linewidth=2, label='Running minimum')
    ax.set_xlabel("Basin-hopping step")
    ax.set_ylabel("Energy (eV)")
    ax.set_title("Basin-hopping search history")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_annealing_history(steps, temps, candidate_energies, best_energies, out_path):
    """Two panels: the cooling schedule actually used (temperature vs. MD
    step) and the annealing analog of plot_search_history (each snapshot's
    relaxed candidate energy, plus the running minimum).
    """
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].plot(steps, temps, color='tab:red')
    axes[0].set_xlabel("MD step")
    axes[0].set_ylabel("Target temperature (K)")
    axes[0].set_title("Cooling schedule")

    axes[1].plot(steps, candidate_energies, 'o', color='tab:blue', alpha=0.5, markersize=4,
                label='Candidate (relaxed) energy')
    axes[1].plot(steps, best_energies, '-', color='tab:orange', linewidth=2, label='Running minimum')
    axes[1].set_xlabel("MD step")
    axes[1].set_ylabel("Energy (eV)")
    axes[1].set_title("Simulated-annealing search history")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def run_basin_hopping(atoms, calc, args, f_out):
    """Runs ase.optimize.basin.BasinHopping in place on `atoms` (mutated to
    the best positions found), writing its own progress/diagnostics plus
    search_history.png/.dat and search_trajectory.<ext> into
    args.output_dir. Returns (best_energy, n_new_minima, traj_out) -- same
    shape as run_simulated_annealing's return, so main() can share the
    reporting tail (final polish + summary) between both algorithms.
    """
    print_dual(f"Running {args.steps} step(s) (this can take a while for larger structures)...", f_out)

    log_path = os.path.join(args.output_dir, "basin_hopping.log")
    traj_path = os.path.join(args.output_dir, "lowest.traj")
    lm_traj_path = os.path.join(args.output_dir, "local_minima.traj")
    optimizer_log_path = os.path.join(args.output_dir, "optimizer.log")

    kT = args.temperature * units.kB
    bh = BasinHopping(atoms, temperature=kT, optimizer=FIRE, fmax=args.fmax, dr=args.dr,
                       logfile=log_path, trajectory=traj_path,
                       optimizer_logfile=optimizer_log_path,
                       local_minima_trajectory=lm_traj_path)
    bh.run(args.steps)

    best_positions = np.array(bh.rmin)
    best_energy = float(bh.Emin)
    atoms.set_positions(best_positions)

    steps_arr, energies_arr, emins_arr = parse_basin_log(log_path)
    n_new_minima = int(np.sum(np.diff(emins_arr) < -1e-9)) if len(emins_arr) > 1 else 0

    print_dual(f"[OK] Search complete: {len(steps_arr)} step(s) logged, "
               f"{n_new_minima} new minima found along the way", f_out)

    history_plot = os.path.join(args.output_dir, "search_history.png")
    if len(steps_arr) > 0:
        plot_search_history(steps_arr, energies_arr, emins_arr, history_plot)
        print_dual(f"[OK] Saved {history_plot}", f_out)
    if args.save_data and len(steps_arr) > 0:
        history_data = os.path.join(args.output_dir, "search_history.dat")
        np.savetxt(history_data, np.column_stack([steps_arr, energies_arr, emins_arr]),
                  header="step attempted_energy(eV) running_min(eV)")
        print_dual(f"[OK] Saved {history_data}", f_out)

    traj_out = None
    if os.path.isfile(lm_traj_path):
        frames = ase_read(lm_traj_path, index=':')
        if frames:
            ase_format, ext = OUTPUT_FORMATS[args.trajectory_format]
            traj_out = os.path.join(args.output_dir, f"search_trajectory{ext}")
            ase_write(traj_out, frames, format=ase_format)
            print_dual(f"[OK] Saved {traj_out} ({len(frames)} frame(s), "
                       f"{args.trajectory_format} format)", f_out)

    return best_energy, n_new_minima, traj_out


def run_simulated_annealing(atoms, calc, args, f_out):
    """Simulated annealing: one continuous NVT MD trajectory (Langevin,
    reusing mlmd.py's build_dynamics for the integrator) whose target
    temperature is cooled from --temp-start down to --temp-end over
    --steps total MD steps (--schedule exponential -- the classic SA
    convention -- or linear, updated on the fly via Langevin.
    set_temperature), snapshotting the running configuration every
    --snapshot-interval steps, locally relaxing each snapshot (positions
    only, same --fmax/--max-steps as basin hopping's per-candidate relax),
    and tracking the lowest-energy relaxed candidate seen -- the annealing
    analog of what BasinHopping's own per-step Metropolis+relax loop does
    for random perturbations, just driven by a cooling MD trajectory
    instead of independent random jumps.

    Mutates `atoms` in place to the best positions found. Returns
    (best_energy, n_new_minima, traj_out) -- same shape as
    run_basin_hopping's return.
    """
    print_dual(f"Running {args.steps} MD step(s), cooling {args.temp_start:.0f} K -> "
               f"{args.temp_end:.0f} K ({args.schedule} schedule)...", f_out)

    dyn = build_dynamics(atoms, "nvt", args.timestep, args.temp_start, args.friction,
                         pressure_gpa=0.0, seed=args.seed)

    n_snapshots = max(1, args.steps // args.snapshot_interval)
    best_energy = atoms.get_potential_energy()
    best_positions = atoms.get_positions().copy()
    n_new_minima = 0

    steps_list, temps_list, candidate_e_list, best_e_list = [], [], [], []
    traj_frames = []

    for i in range(n_snapshots):
        frac = i / max(1, n_snapshots - 1)
        if args.schedule == "exponential":
            temp_now = args.temp_start * (args.temp_end / args.temp_start) ** frac
        else:
            temp_now = args.temp_start + (args.temp_end - args.temp_start) * frac

        dyn.set_temperature(temperature_K=temp_now)
        dyn.run(args.snapshot_interval)
        step_now = (i + 1) * args.snapshot_interval

        candidate = atoms.copy()
        candidate.calc = calc
        mace_relax.relax(candidate, calc, cell_mask=None, fmax=args.fmax, max_steps=args.max_steps)
        candidate_energy = candidate.get_potential_energy()

        if candidate_energy < best_energy - 1e-9:
            best_energy = candidate_energy
            best_positions = candidate.get_positions().copy()
            n_new_minima += 1

        steps_list.append(step_now)
        temps_list.append(temp_now)
        candidate_e_list.append(candidate_energy)
        best_e_list.append(best_energy)
        traj_frames.append(Atoms(numbers=candidate.numbers, positions=candidate.positions,
                                 cell=candidate.cell, pbc=True))

    atoms.set_positions(best_positions)

    print_dual(f"[OK] Search complete: {len(steps_list)} snapshot(s) evaluated, "
               f"{n_new_minima} new minima found along the way", f_out)

    history_plot = os.path.join(args.output_dir, "search_history.png")
    plot_annealing_history(steps_list, temps_list, candidate_e_list, best_e_list, history_plot)
    print_dual(f"[OK] Saved {history_plot}", f_out)
    if args.save_data:
        history_data = os.path.join(args.output_dir, "search_history.dat")
        np.savetxt(history_data, np.column_stack([steps_list, temps_list, candidate_e_list, best_e_list]),
                  header="step temperature(K) candidate_energy(eV) running_min(eV)")
        print_dual(f"[OK] Saved {history_data}", f_out)

    traj_out = None
    if traj_frames:
        ase_format, ext = OUTPUT_FORMATS[args.trajectory_format]
        traj_out = os.path.join(args.output_dir, f"search_trajectory{ext}")
        ase_write(traj_out, traj_frames, format=ase_format)
        print_dual(f"[OK] Saved {traj_out} ({len(traj_frames)} frame(s), "
                   f"{args.trajectory_format} format)", f_out)

    return best_energy, n_new_minima, traj_out


def main():
    parser = argparse.ArgumentParser(
        description="Global structure search driven entirely by a MACE potential -- the "
                     "MACE-MP-0 foundation model, or a custom model fine-tuned on your own "
                     "SIESTA data via stb-mlffAnalysis (--custom-model). Two algorithms "
                     "(--algorithm): 'basin-hopping' (default, ase.optimize.basin."
                     "BasinHopping -- random perturbation + local relax + Metropolis accept/"
                     "reject) or 'simulated-annealing' (one continuous cooling NVT MD "
                     "trajectory, periodically snapshotted and locally relaxed). Both search "
                     "for the lowest-energy atomic arrangement at a FIXED cell -- e.g. "
                     "defect/vacancy site occupation, disordered/alloy configurations, local "
                     "amorphous-like rearrangements. NOT crystal structure prediction (the "
                     "cell shape itself is untouched during the search) and NOT a proof of "
                     "finding the true global minimum -- both are heuristic searches.",
        epilog="Example usage:\n"
               "  stb-mlsearch --file structure.fdf --steps 100\n"
               "  stb-mlsearch --file structure.fdf --custom-model mlff_model.model "
               "--dr 0.4 --temperature 1500\n"
               "  stb-mlsearch --file structure.fdf --algorithm simulated-annealing "
               "--steps 2000 --temp-start 1500 --temp-end 50\n",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--file", required=True, help="Input structure (.fdf).")
    parser.add_argument("--model", choices=["small", "medium", "large"], default="small",
                        help="MACE-MP-0 model size (default: small). Ignored if --custom-model is given.")
    parser.add_argument("--custom-model", default=None, metavar="PATH",
                        help="Use a custom MACE model file instead of the MACE-MP-0 foundation "
                             "potential -- e.g. one fine-tuned via stb-mlffAnalysis. Overrides --model.")
    parser.add_argument("--algorithm", choices=["basin-hopping", "simulated-annealing"],
                        default="basin-hopping",
                        help="Search algorithm (default: basin-hopping). 'basin-hopping': "
                             "random perturbation (--dr) + local relax + Metropolis accept/"
                             "reject (--temperature). 'simulated-annealing': one continuous "
                             "cooling NVT MD trajectory (--temp-start/--temp-end/--schedule), "
                             "periodically snapshotted (--snapshot-interval) and locally "
                             "relaxed.")
    parser.add_argument("--steps", type=int, default=50,
                        help="Number of basin-hopping perturb+relax+accept/reject attempts "
                             "(basin-hopping), or total MD steps in the cooling trajectory "
                             "(simulated-annealing) (default: 50).")
    parser.add_argument("--dr", type=float, default=0.3,
                        help="[basin-hopping] Max random atomic displacement per step, "
                             "Angstrom, uniform per Cartesian component (default: 0.3).")
    parser.add_argument("--temperature", type=float, default=1000.0,
                        help="[basin-hopping] Search 'temperature', Kelvin, controlling the "
                             "Metropolis acceptance probability exp(-dE/kT) for an uphill move "
                             "(default: 1000 -- deliberately higher than any physical "
                             "temperature, to encourage escaping shallow local minima; a "
                             "common basin-hopping convention, not a claim about the "
                             "material's real conditions).")
    parser.add_argument("--temp-start", type=float, default=1500.0,
                        help="[simulated-annealing] Starting temperature, Kelvin (default: 1500).")
    parser.add_argument("--temp-end", type=float, default=50.0,
                        help="[simulated-annealing] Final temperature, Kelvin (default: 50 -- "
                             "near-zero, not exactly 0, since the cooling schedule is "
                             "multiplicative for --schedule exponential).")
    parser.add_argument("--schedule", choices=["exponential", "linear"], default="exponential",
                        help="[simulated-annealing] Cooling schedule from --temp-start to "
                             "--temp-end (default: exponential -- the classic simulated-"
                             "annealing convention, cools fast initially then slows down).")
    parser.add_argument("--snapshot-interval", type=int, default=20,
                        help="[simulated-annealing] MD steps between snapshot+local-relax "
                             "candidate evaluations (default: 20).")
    parser.add_argument("--timestep", type=float, default=1.0,
                        help="[simulated-annealing] MD timestep, fs (default: 1.0).")
    parser.add_argument("--friction", type=float, default=0.01,
                        help="[simulated-annealing] Langevin friction, fs^-1 (default: 0.01).")
    parser.add_argument("--fmax", type=float, default=0.05,
                        help="Force convergence threshold for each candidate's local relax, "
                             "eV/Ang (default: 0.05).")
    parser.add_argument("--max-steps", type=int, default=200,
                        help="Max local-optimizer steps per candidate relax (default: 200).")
    parser.add_argument("--seed", type=int, default=None,
                        help="Random seed for reproducible searches (default: unseeded).")
    parser.add_argument("--no-pre-relax", dest="pre_relax", action="store_false",
                        help="Skip the initial full MACE relax (cell + positions) of the "
                             "input structure before starting the search. On by default, so "
                             "the search starts from a sane reference rather than a raw input.")
    parser.add_argument("--no-final-cell-relax", dest="final_cell_relax", action="store_false",
                        help="Skip the final full relax (cell + positions) of the winning "
                             "configuration. On by default -- the search itself only moves "
                             "atomic positions at a fixed cell (see module docstring), so the "
                             "final structure otherwise keeps the pre-search cell shape "
                             "exactly, even if the winning arrangement would relax to a "
                             "slightly different one.")
    parser.add_argument("--vacuum-gap", type=float, default=10.0,
                        help="Gap (Ang) used to detect vacuum-padded axes, for the final cell "
                             "relax's mask (default: 10.0).")
    parser.add_argument("--trajectory-format", choices=sorted(OUTPUT_FORMATS), default="xsf",
                        help="Output format for the full search trajectory (every attempted, "
                             "locally-relaxed candidate), same convention as stb-ani2traj/"
                             "stb-mlmd/stb-mlphonons (default: xsf).")
    parser.add_argument("--device", choices=["cpu", "cuda"], default="cpu",
                        help="Device to run the model on (default: cpu).")
    parser.add_argument("-o", "--output-dir", default="mlsearch_out",
                        help="Output directory for all files (default: mlsearch_out).")
    parser.add_argument("--save-data", action="store_true",
                        help="Also write the raw search-history data (.dat). Off by default.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-mlsearch {VERSION}")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - ML Structure Search",
            "Basin-hopping global search for the lowest-energy arrangement, via MACE",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text(f"ML STRUCTURE SEARCH ({args.algorithm.upper()}):", 'bold'))
    print("-" * 60)

    os.makedirs(args.output_dir, exist_ok=True)
    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    def fail(message):
        print_dual(color_text(f"[FAIL] {message}", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_dual(color_text("===== STB-MLSEARCH REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Input file        : {args.file}", f_out)
    print_dual(f"Model             : "
               f"{'custom (' + args.custom_model + ')' if args.custom_model else f'MACE-MP-0 ({args.model})'}", f_out)
    print_dual(f"Algorithm         : {args.algorithm}", f_out)
    print_dual(f"Steps             : {args.steps}", f_out)
    if args.algorithm == "basin-hopping":
        print_dual(f"Max displacement  : {args.dr} Ang/atom/step", f_out)
        print_dual(f"Search temperature: {args.temperature} K", f_out)
    else:
        print_dual(f"Temperature range : {args.temp_start:.0f} K -> {args.temp_end:.0f} K "
                   f"({args.schedule})", f_out)
        print_dual(f"Snapshot interval : every {args.snapshot_interval} MD step(s)", f_out)
    if args.seed is not None:
        print_dual(f"Random seed       : {args.seed}", f_out)
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
    print_dual(f"Dimensionality    : {kspace.dimensionality_label(vacuum_axes)}", f_out)

    if args.seed is not None:
        np.random.seed(args.seed)

    model_arg = args.custom_model if args.custom_model else args.model
    calc = mace_relax.get_calculator(model=model_arg, device=args.device)
    atoms.calc = calc

    print_section("[2] INITIAL RELAXATION", f_out)
    e0 = atoms.get_potential_energy()
    f0 = np.abs(atoms.get_forces()).max()
    print_dual(f"Input energy      : {e0:.6f} eV ({e0/len(atoms):.6f} eV/atom), "
               f"max|F| = {f0:.4f} eV/Ang", f_out)
    if args.pre_relax:
        cell_mask = mace_relax.build_cell_mask(vacuum_axes)
        converged, steps_used = mace_relax.relax(
            atoms, calc, cell_mask=cell_mask, fmax=args.fmax, max_steps=args.max_steps)
        e1 = atoms.get_potential_energy()
        f1 = np.abs(atoms.get_forces()).max()
        print_dual(f"After pre-relax   : {e1:.6f} eV ({e1/len(atoms):.6f} eV/atom), "
                   f"max|F| = {f1:.4f} eV/Ang ({'converged' if converged else 'hit step cap'}, "
                   f"{steps_used} steps)", f_out)
    else:
        print_dual(color_text(
            "[WARNING] --no-pre-relax: starting the search from the structure as given.", 'yellow'), f_out)

    section_title = "BASIN-HOPPING SEARCH" if args.algorithm == "basin-hopping" else "SIMULATED-ANNEALING SEARCH"
    print_section(f"[3] {section_title}", f_out)
    if args.algorithm == "basin-hopping":
        best_energy, n_new_minima, traj_out = run_basin_hopping(atoms, calc, args, f_out)
    else:
        best_energy, n_new_minima, traj_out = run_simulated_annealing(atoms, calc, args, f_out)

    print_dual(f"Best energy found : {best_energy:.6f} eV ({best_energy/len(atoms):.6f} eV/atom)", f_out)
    improvement = e0 - best_energy
    print_dual(f"Improvement vs. input energy: {improvement:.6f} eV "
               f"({'lower' if improvement > 0 else 'no improvement'})", f_out)
    history_plot = os.path.join(args.output_dir, "search_history.png")

    print_section("[4] FINAL POLISH", f_out)
    if args.final_cell_relax:
        cell_mask = mace_relax.build_cell_mask(vacuum_axes)
        converged, steps_used = mace_relax.relax(
            atoms, calc, cell_mask=cell_mask, fmax=args.fmax, max_steps=args.max_steps)
        final_energy = atoms.get_potential_energy()
        print_dual(f"Final relax (cell+ions): {best_energy:.6f} -> {final_energy:.6f} eV "
                   f"({'converged' if converged else 'hit step cap'}, {steps_used} steps)", f_out)
    else:
        final_energy = best_energy
        print_dual("Skipped (--no-final-cell-relax) -- winning structure keeps the "
                   "pre-search cell shape exactly.", f_out)

    pmg_out = AseAtomsAdaptor.get_structure(atoms)
    best_structure = structure_io.from_pymatgen(pmg_out, species_meta=structure.species_meta)
    best_path = os.path.join(args.output_dir, "best_structure.fdf")
    structure_io.write_fdf(best_structure, best_path)
    print_dual(f"[OK] Saved {best_path}", f_out)

    print_section("[5] SUMMARY & FILES", f_out)
    print_dual("Status            : OK", f_out)
    print_dual(f"Input energy      : {e0:.6f} eV ({e0/len(atoms):.6f} eV/atom)", f_out)
    print_dual(f"Best energy       : {final_energy:.6f} eV ({final_energy/len(atoms):.6f} eV/atom)", f_out)
    print_dual(f"New minima found  : {n_new_minima}", f_out)
    print_dual(f"Best structure    : {best_path}", f_out)
    if os.path.isfile(history_plot):
        print_dual(f"Search history    : {history_plot}", f_out)
    if traj_out:
        print_dual(f"Search trajectory : {traj_out}", f_out)
    if report_path:
        print_dual(f"Report            : {report_path}", f_out)

    if f_out:
        f_out.close()

    print("\n" + "-" * 60)
    print(color_text("ML structure search complete -- a heuristic search, not a proof of "
                      "global optimality.\n", 'bold'))


if __name__ == "__main__":
    main()
