#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""Model-size convergence check for MACE-MP-0 -- the sixth tool in the ML
Simulations category. Runs the exact same reference calculation (a full
relax, cell + positions, via core/mace_relax.py's relax/build_cell_mask --
same vacuum-aware convention as stb-mlrelax) independently with each
requested model size (--models, default small medium large) and any
--custom-models given, then reports how much the answer actually changes
between them: energy/atom, cell volume, lattice parameters, max residual
force, and wall-clock time. Helps decide which model size is "good enough"
for a given structure/property instead of guessing -- e.g. if medium and
large agree to within --energy-tolerance/--lattice-tolerance, small is
probably fine and paying large's much higher cost buys nothing; if they
disagree significantly, size matters here and small shouldn't be trusted.

Not a substitute for validating against real DFT -- this only checks
whether the MACE-MP-0 answer itself has converged with model size, not
whether MACE-MP-0 (at any size) agrees with the real material.
"""

VERSION = "1.0.0"

import os
import sys
import time
import argparse
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt
from pymatgen.io.ase import AseAtomsAdaptor

from stb.core import structure_io, kspace, mace_relax
from stb.core.cli import color_text, show_intro, print_dual, print_section
from stb.core.deps import require_mace

require_mace()

REPORT_FILE = "stb_mlconvergence_report.txt"


def run_one_model(atoms_ref, model_arg, device, args, vacuum_axes, f_out, tag):
    """Runs the reference calculation with one model (a foundation size
    string or a custom model file path) and returns a dict of everything
    the comparison/report needs: {"label", "energy_per_atom", "volume",
    "cellpar", "max_force", "converged", "steps_used", "wall_time_s"}.
    """
    calc = mace_relax.get_calculator(model=model_arg, device=device)
    atoms = atoms_ref.copy()
    atoms.calc = calc

    t0 = time.time()
    converged, steps_used = True, 0
    if args.relax:
        cell_mask = mace_relax.build_cell_mask(vacuum_axes)
        converged, steps_used = mace_relax.relax(
            atoms, calc, cell_mask=cell_mask, fmax=args.fmax, max_steps=args.max_steps)
    energy = atoms.get_potential_energy()
    max_force = np.abs(atoms.get_forces()).max()
    wall_time = time.time() - t0

    result = {
        "label": tag,
        "energy_per_atom": energy / len(atoms),
        "volume": atoms.get_volume(),
        "cellpar": atoms.cell.cellpar(),
        "max_force": max_force,
        "converged": converged,
        "steps_used": steps_used,
        "wall_time_s": wall_time,
    }
    print_dual(f"[{tag}] E/atom = {result['energy_per_atom']:.6f} eV, "
               f"V = {result['volume']:.4f} Ang^3, max|F| = {max_force:.4f} eV/Ang, "
               f"{wall_time:.1f} s ({'converged' if converged else 'hit step cap'}, "
               f"{steps_used} steps)", f_out)
    return result


def plot_convergence(results, energy_tolerance, out_path):
    """3-panel bar/line comparison across models: energy/atom, cell
    volume, and wall-clock time -- the three quantities a user actually
    weighs when picking a model size.
    """
    labels = [r["label"] for r in results]
    energies = [r["energy_per_atom"] for r in results]
    volumes = [r["volume"] for r in results]
    times = [r["wall_time_s"] for r in results]

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))
    axes[0].plot(labels, energies, 'o-', color='tab:blue')
    axes[0].set_ylabel("Energy/atom (eV)")
    axes[0].set_title("Energy convergence")
    axes[0].tick_params(axis='x', rotation=30)

    axes[1].plot(labels, volumes, 'o-', color='tab:green')
    axes[1].set_ylabel("Cell volume (Ang^3)")
    axes[1].set_title("Volume convergence")
    axes[1].tick_params(axis='x', rotation=30)

    axes[2].bar(labels, times, color='tab:orange')
    axes[2].set_ylabel("Wall time (s)")
    axes[2].set_title("Cost")
    axes[2].tick_params(axis='x', rotation=30)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Model-size convergence check for MACE-MP-0 -- runs the same reference "
                     "calculation (full relax, cell + positions) independently with each "
                     "requested model size (--models, default small medium large) and any "
                     "--custom-models, then reports how much energy/atom, cell volume/"
                     "lattice parameters, and wall-clock time actually change between them. "
                     "Helps pick a model size instead of guessing -- NOT a check against real "
                     "DFT, only whether MACE-MP-0's own answer has converged with size.",
        epilog="Example usage:\n"
               "  stb-mlconvergence --file structure.fdf\n"
               "  stb-mlconvergence --file structure.fdf --models small medium "
               "--custom-models mlff_model.model\n",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--file", required=True, help="Input structure (.fdf).")
    parser.add_argument("--models", nargs='+', choices=["small", "medium", "large"],
                        default=["small", "medium", "large"],
                        help="Foundation model sizes to compare (default: small medium large).")
    parser.add_argument("--custom-models", nargs='+', default=None, metavar="PATH",
                        help="Also include one or more custom MACE model files (e.g. fine-"
                             "tuned via stb-mlffAnalysis) in the comparison, labeled by "
                             "filename.")
    parser.add_argument("--no-relax", dest="relax", action="store_false",
                        help="Skip the full relax (cell + positions) -- just a single-point "
                             "energy/force evaluation on the structure as given, for each "
                             "model. On by default (a relax is what most workflows actually "
                             "care about converging).")
    parser.add_argument("--fmax", type=float, default=0.05,
                        help="Force convergence threshold for the relax, eV/Ang (default: 0.05).")
    parser.add_argument("--max-steps", type=int, default=200,
                        help="Max optimizer steps for the relax (default: 200).")
    parser.add_argument("--vacuum-gap", type=float, default=10.0,
                        help="Gap (Ang) used to detect vacuum-padded axes, for the relax's "
                             "cell mask (default: 10.0).")
    parser.add_argument("--energy-tolerance", type=float, default=5.0,
                        help="Flag a pair of consecutive models as NOT converged if their "
                             "energy/atom differs by more than this, meV/atom (default: 5.0).")
    parser.add_argument("--lattice-tolerance", type=float, default=1.0,
                        help="Flag a pair of consecutive models as NOT converged if their "
                             "cell volume differs by more than this, %% (default: 1.0).")
    parser.add_argument("--device", choices=["cpu", "cuda"], default="cpu",
                        help="Device to run the models on (default: cpu).")
    parser.add_argument("-o", "--output-dir", default="mlconvergence_out",
                        help="Output directory for all files (default: mlconvergence_out).")
    parser.add_argument("--save-data", action="store_true",
                        help="Also write the raw per-model comparison data (.dat). Off by default.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-mlconvergence {VERSION}")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - ML Model-Size Convergence",
            "Compares MACE-MP-0 model sizes on the same reference calculation",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("ML MODEL-SIZE CONVERGENCE CHECK:", 'bold'))
    print("-" * 60)

    os.makedirs(args.output_dir, exist_ok=True)
    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    def fail(message):
        print_dual(color_text(f"[FAIL] {message}", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_dual(color_text("===== STB-MLCONVERGENCE REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Input file        : {args.file}", f_out)
    print_dual(f"Models compared   : {', '.join(args.models)}"
               f"{' + custom: ' + ', '.join(args.custom_models) if args.custom_models else ''}", f_out)
    print_dual(f"Relax             : {'yes (cell+positions)' if args.relax else 'no (single-point only)'}", f_out)
    if report_path:
        print_dual(f"Report file       : {report_path}", f_out)

    print_section("[1] READING STRUCTURE", f_out)

    if not os.path.isfile(args.file):
        fail(f"Input file '{args.file}' not found.")
    for cm in (args.custom_models or []):
        if not os.path.isfile(cm):
            fail(f"--custom-models file not found: {cm}")

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

    print_section("[2] PER-MODEL CALCULATION", f_out)
    model_specs = [(m, m) for m in args.models]
    for cm in (args.custom_models or []):
        model_specs.append((cm, os.path.basename(cm)))

    results = []
    for model_arg, tag in model_specs:
        results.append(run_one_model(atoms, model_arg, args.device, args, vacuum_axes, f_out, tag))

    print_section("[3] CONVERGENCE COMPARISON", f_out)
    any_flagged = False
    print_dual(f"{'Pair':<24}{'dE/atom (meV)':>16}{'dVolume (%)':>14}{'Status':>12}", f_out)
    for i in range(1, len(results)):
        r0, r1 = results[i - 1], results[i]
        de_mev = abs(r1["energy_per_atom"] - r0["energy_per_atom"]) * 1000.0
        dvol_pct = 100.0 * abs(r1["volume"] - r0["volume"]) / r0["volume"] if r0["volume"] > 0 else 0.0
        flagged = de_mev > args.energy_tolerance or dvol_pct > args.lattice_tolerance
        any_flagged = any_flagged or flagged
        status = color_text('NOT CONVERGED', 'yellow') if flagged else color_text('OK', 'green')
        pair_label = f"{r0['label']} -> {r1['label']}"
        print_dual(f"{pair_label:<24}{de_mev:16.3f}{dvol_pct:14.3f}{status:>21}", f_out)

    if any_flagged:
        print_dual(color_text(
            "[WARNING] At least one consecutive pair disagrees by more than "
            f"--energy-tolerance ({args.energy_tolerance} meV/atom) and/or --lattice-tolerance "
            f"({args.lattice_tolerance}%) -- results may still be changing with model size; "
            "consider the largest model compared, or a fine-tuned one.", 'yellow'), f_out)
    else:
        print_dual(color_text(
            "[OK] All consecutive pairs agree within tolerance -- the smallest model compared "
            "is likely sufficient for this structure/property.", 'green'), f_out)

    plot_path = os.path.join(args.output_dir, "convergence.png")
    plot_convergence(results, args.energy_tolerance, plot_path)
    print_dual(f"[OK] Saved {plot_path}", f_out)
    if args.save_data:
        data_path = os.path.join(args.output_dir, "convergence.dat")
        with open(data_path, "w") as f:
            f.write("# label energy_per_atom(eV) volume(Ang^3) max_force(eV/Ang) wall_time(s)\n")
            for r in results:
                f.write(f"{r['label']} {r['energy_per_atom']:.8f} {r['volume']:.6f} "
                        f"{r['max_force']:.6f} {r['wall_time_s']:.3f}\n")
        print_dual(f"[OK] Saved {data_path}", f_out)

    print_section("[4] SUMMARY & FILES", f_out)
    print_dual("Status            : OK", f_out)
    fastest_converged = None
    for i in range(1, len(results)):
        r0, r1 = results[i - 1], results[i]
        de_mev = abs(r1["energy_per_atom"] - r0["energy_per_atom"]) * 1000.0
        dvol_pct = 100.0 * abs(r1["volume"] - r0["volume"]) / r0["volume"] if r0["volume"] > 0 else 0.0
        if de_mev <= args.energy_tolerance and dvol_pct <= args.lattice_tolerance:
            fastest_converged = r0["label"]
            break
    if fastest_converged:
        print_dual(f"Recommended       : {fastest_converged} (smallest model already within tolerance)", f_out)
    print_dual(f"Convergence plot  : {plot_path}", f_out)
    if report_path:
        print_dual(f"Report            : {report_path}", f_out)

    if f_out:
        f_out.close()

    print("\n" + "-" * 60)
    print(color_text("ML model-size convergence check complete.\n", 'bold'))


if __name__ == "__main__":
    main()
