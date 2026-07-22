#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""Standalone equation-of-state (E vs V) calculation driven entirely by a
MACE potential -- scans a set of isotropically-scaled cell volumes around
the (optionally pre-relaxed) equilibrium, computes a MACE total energy at
each one (positions relaxed, cell fixed at the scaled volume), and fits an
E(V) equation of state (Birch-Murnaghan/Vinet/Murnaghan/stabilized-jellium,
via ase.eos.EquationOfState) to get V0, E0, B0 (bulk modulus) and B0'
(pressure derivative). No phonons/force-constants at all -- this is the
static-lattice EOS only, the ML Simulations analog of a bare volume-scan
DFT convergence study.

Deliberately distinct from stb-mlelastic: that tool derives the bulk
modulus from the full elastic STIFFNESS TENSOR (a linear stress-strain fit
around a single reference cell), while this one derives it from the
curvature of the total-energy-vs-volume curve itself (B0 = V0 * d^2E/dV^2
at V0, captured by the EOS fit). The two use independent methods and
independent MACE evaluations (isotropic volume scan vs. canonical Voigt
strains) on the same structure -- their B0 values should agree for a
reasonably isotropic material, and comparing them is a useful physical
consistency check on the MACE potential itself, not just a duplicate
feature. Both report bulk modulus in GPa via the same conversion factor
(1 eV/Ang^3 -> GPa, core/eos_fit.py's CONV_EVA3_TO_GPA -- same value as
elastic_analysis.py's own constant of the same name), so the two numbers
are directly comparable without any extra unit bookkeeping.

The EOS-fitting math itself (ase.eos.EquationOfState wrapper, GPa
conversion, R^2) lives in core/eos_fit.py, not here -- extracted once
stb-eosAnalysis (the real-SIESTA analog of this tool) became a second
consumer of the exact same fit.

Bulk (3D periodic) only -- an E-V equation of state assumes a well-defined
cell volume; a vacuum-padded axis (slab/wire/molecule) has no bound volume
to scan (same reasoning as stb-mlmelting/stb-amorphize).

Isotropic volume-scaling convention (factor = (1 + strain_pct/100)**(1/3)
applied to the cell vectors, so --strain-range is a PERCENT VOLUME change,
matching stb-mlphonons's --qha volume scan exactly) is reused so a user
who has already picked a --strain-range for one tool can reuse the same
number in the other with the same meaning.
"""

VERSION = "1.0.2"  # v0-outside-scanned-range warning surfaced (core/eos_fit.py's fit_eos fix)

import os
import sys
import argparse
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt
from pymatgen.io.ase import AseAtomsAdaptor

from stb.core import structure_io, kspace, mace_relax
from stb.core.cli import color_text, show_intro, print_dual, print_section
from stb.core.deps import require_mace
from stb.core.eos_fit import normalize_eos_string, fit_eos

require_mace()

REPORT_FILE = "stb_mleos_report.txt"
_SERIES_COLORS = ['tab:blue', 'tab:orange']


def run_volume_scan(atoms_ref, calc, args, f_out, tag=""):
    """Scans --n-volumes isotropically-scaled volumes (+/- --strain-range %
    around atoms_ref's own volume), relaxing positions only (cell fixed at
    each scaled volume) unless --no-relax. Returns (volumes, energies) in
    ascending-volume order (np.linspace is already monotonic in strain, and
    volume is a monotonic function of the isotropic scale factor, so no
    re-sorting is needed).
    """
    v0_input = atoms_ref.get_volume()
    strains_pct = np.linspace(-abs(args.strain_range), abs(args.strain_range), args.n_volumes)

    print_dual(f"Volume range{tag}      : {-abs(args.strain_range):.1f}% to "
               f"{abs(args.strain_range):.1f}% ({args.n_volumes} points, V0={v0_input:.4f} Ang^3)", f_out)
    print_dual(f"{'Strain%':>8} {'Volume(Ang^3)':>14} {'Energy(eV)':>14} {'Energy/atom(eV)':>16}", f_out)

    volumes, energies = [], []
    n_atoms = len(atoms_ref)
    for strain_pct in strains_pct:
        factor = (1.0 + strain_pct / 100.0) ** (1.0 / 3.0)
        atoms_i = atoms_ref.copy()
        atoms_i.set_cell(np.array(atoms_ref.get_cell()) * factor, scale_atoms=True)
        atoms_i.calc = calc

        if args.relax:
            mace_relax.relax(atoms_i, calc, cell_mask=None, fmax=args.fmax, max_steps=args.max_steps)

        energy_i = atoms_i.get_potential_energy()
        volume_i = atoms_i.get_volume()
        volumes.append(volume_i)
        energies.append(energy_i)

        print_dual(f"{strain_pct:8.2f} {volume_i:14.4f} {energy_i:14.6f} {energy_i / n_atoms:16.6f}", f_out)

    return np.array(volumes), np.array(energies)


def plot_eos(series, out_path):
    """series: list of (label, volumes, energies, fit, color) -- scatter
    the scanned points and overlay each fitted curve; 2 series (fine-tuned
    vs foundation) side by side when comparing, 1 otherwise.
    """
    fig, ax = plt.subplots(figsize=(7, 5))
    for label, volumes, energies, fit, color in series:
        ax.plot(volumes, energies, 'o', color=color, label=f"{label} (data)")
        ax.plot(fit["curve_v"], fit["curve_e"], '-', color=color,
                 label=f"{label} (fit, V0={fit['v0']:.3f} Ang^3)")
    ax.set_xlabel("Volume (Ang^3)")
    ax.set_ylabel("Energy (eV)")
    ax.set_title("Equation of State (MACE)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def report_eos_fit(fit, n_atoms, unit_label, f_out, tag=""):
    v0_str = color_text(f"{fit['v0']:.4f}", 'bold')
    e0_str = color_text(f"{fit['e0']:.6f}", 'bold')
    b0_str = color_text(f"{fit['b0_gpa']:.2f}", 'bold')
    print_dual(f"Equilibrium volume{tag} (V0) : {v0_str} Ang^3 "
               f"({fit['v0'] / n_atoms:.4f} Ang^3/atom)", f_out)
    print_dual(f"Equilibrium energy{tag} (E0) : {e0_str} eV "
               f"({fit['e0'] / n_atoms:.6f} eV/atom)", f_out)
    print_dual(f"Bulk modulus{tag}       (B0) : {b0_str} GPa", f_out)
    if fit["bprime"] is not None:
        bprime_str = color_text(f"{fit['bprime']:.3f}", 'bold')
        print_dual(f"Pressure derivative{tag} (B0'): {bprime_str}", f_out)
    print_dual(f"Fit quality{tag}        (R^2): {fit['r_squared']:.6f}", f_out)
    if fit["v0_outside_range"]:
        print_dual(color_text(
            f"[WARNING] The fitted equilibrium volume{tag} (V0) falls OUTSIDE the range of "
            "scanned volumes -- the scan didn't actually bracket the true minimum; widen "
            "--strain-range and re-run.", 'yellow'), f_out)


def main():
    parser = argparse.ArgumentParser(
        description="Standalone equation-of-state (E vs V) calculation driven entirely by a "
                     "MACE potential -- the MACE-MP-0 foundation model, or a custom model "
                     "fine-tuned on your own SIESTA data via stb-mlffAnalysis (--custom-model). "
                     "Scans isotropically-scaled volumes around equilibrium and fits V0/E0/B0/B0' "
                     "(Birch-Murnaghan by default). Gives an independent, curvature-based bulk "
                     "modulus to cross-check against stb-mlelastic's stress-strain-derived one. "
                     "A fast heuristic, not a substitute for a real SIESTA volume scan.",
        epilog="Example usage:\n"
               "  stb-mleos --file structure.fdf\n"
               "  stb-mleos --file structure.fdf --custom-model mlff_model.model "
               "--n-volumes 11 --strain-range 8.0\n"
               "  stb-mleos --file structure.fdf --eos vinet\n",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--file", required=True, help="Input structure (.fdf).")
    parser.add_argument("--model", choices=["small", "medium", "large"], default="small",
                        help="MACE-MP-0 model size (default: small). Ignored if --custom-model is given.")
    parser.add_argument("--custom-model", default=None, metavar="PATH",
                        help="Use a custom MACE model file instead of the MACE-MP-0 foundation "
                             "potential -- e.g. one fine-tuned via stb-mlffAnalysis. Overrides --model.")
    parser.add_argument("--skip-foundation-comparison", action="store_true",
                        help="If --custom-model is given, also fit the EOS with the raw "
                             "(non-fine-tuned) foundation model (--model) and overlay both curves "
                             "plus a side-by-side V0/E0/B0/B0' table. On by default; pass this to "
                             "skip it.")
    parser.add_argument("--no-pre-relax", dest="pre_relax", action="store_false",
                        help="Skip the full MACE pre-relax (cell + positions) of the reference "
                             "structure before scanning volumes. On by default -- scanning around "
                             "a structure that isn't already at a MACE-relaxed minimum shifts V0 "
                             "away from the true equilibrium and can leave it outside the scanned "
                             "range entirely.")
    parser.add_argument("--no-relax", dest="relax", action="store_false",
                        help="Skip the per-volume atomic-position relax (cell fixed at the scaled "
                             "volume). On by default -- for any structure with more than one atom "
                             "per primitive cell, internal coordinates generally need to relax at "
                             "each fixed strained cell for the energy to be physically meaningful.")
    parser.add_argument("--fmax", type=float, default=0.02,
                        help="Force convergence threshold for both relax stages, eV/Ang (default: 0.02).")
    parser.add_argument("--max-steps", type=int, default=200,
                        help="Max optimizer steps for both relax stages (default: 200).")
    parser.add_argument("--n-volumes", type=int, default=9,
                        help="Number of volumes to scan (default: 9). Needs >= 5 for a robust "
                             "4-parameter EOS fit.")
    parser.add_argument("--strain-range", type=float, default=5.0,
                        help="Max isotropic volume strain, %% (default: 5.0 -- wider than "
                             "stb-mlelastic's linear-regime strain, since an EOS fit wants "
                             "genuine curvature in the E-V curve, not a small-strain linear "
                             "response).")
    parser.add_argument("--eos", choices=["birchmurnaghan", "birch_murnaghan", "vinet", "murnaghan", "sj"],
                        default="birchmurnaghan",
                        help="Equation-of-state form to fit (default: birchmurnaghan). 'sj' = "
                             "stabilized jellium (no B0' reported).")
    parser.add_argument("--vacuum-gap", type=float, default=10.0,
                        help="Gap (Ang) used to detect vacuum-padded axes -- bulk-only, rejects "
                             "any vacuum-padded structure (default: 10.0).")
    parser.add_argument("--device", choices=["cpu", "cuda"], default="cpu",
                        help="Device to run the model on (default: cpu).")
    parser.add_argument("-o", "--output-dir", default="mleos_out",
                        help="Output directory for the plot/data (default: mleos_out).")
    parser.add_argument("--save-data", action="store_true",
                        help="Also write the raw volume-energy data (.dat) behind the plot. Off by default.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-mleos {VERSION}")

    args = parser.parse_args()
    args.eos = normalize_eos_string(args.eos)

    if args.n_volumes < 5:
        parser.error("--n-volumes needs at least 5 points for a robust EOS fit.")

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - ML Equation of State",
            "E vs V curve + equation-of-state fit, driven entirely by a MACE potential",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("ML EQUATION OF STATE:", 'bold'))
    print("-" * 60)

    os.makedirs(args.output_dir, exist_ok=True)
    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    def fail(message):
        print_dual(color_text(f"[FAIL] {message}", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_dual(color_text("===== STB-MLEOS REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Input file        : {args.file}", f_out)
    print_dual(f"Model             : "
               f"{'custom (' + args.custom_model + ')' if args.custom_model else f'MACE-MP-0 ({args.model})'}", f_out)
    print_dual(f"EOS form          : {args.eos}", f_out)
    print_dual(f"Volume scan       : +/-{args.strain_range}% ({args.n_volumes} points)", f_out)
    compare = bool(args.custom_model) and not args.skip_foundation_comparison
    if compare:
        print_dual(f"Comparison        : also fitting with foundation model ({args.model})", f_out)
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
    n_atoms = len(atoms)

    lattice_ang = np.array(atoms.get_cell())
    frac_coords = atoms.get_scaled_positions()
    vacuum_axes = kspace.detect_vacuum_axes(frac_coords, lattice_ang, args.vacuum_gap)
    print_dual(f"[OK] Read {n_atoms} atom(s), "
               f"dimensionality: {kspace.dimensionality_label(vacuum_axes)}", f_out)

    if any(vacuum_axes):
        fail("This tool is bulk (3D periodic) only -- a vacuum-padded axis was detected. "
             "An equation of state assumes a well-defined, bounded cell volume, which a "
             "slab/wire/molecule doesn't have.")

    model_arg = args.custom_model if args.custom_model else args.model
    calc = mace_relax.get_calculator(model=model_arg, device=args.device)

    print_section("[2] REFERENCE PRE-RELAX", f_out)
    atoms.calc = calc
    if args.pre_relax:
        converged, steps_used = mace_relax.relax(atoms, calc, cell_mask=[True] * 6,
                                                  fmax=args.fmax, max_steps=args.max_steps)
        status = color_text("converged", 'green') if converged else color_text("NOT converged", 'red')
        print_dual(f"[OK] Pre-relax (cell + positions): {status} in {steps_used} step(s).", f_out)
    else:
        print_dual("[SKIPPED] --no-pre-relax given; scanning around the input structure as-is.", f_out)
    v0_ref = atoms.get_volume()
    print_dual(f"Reference volume  : {v0_ref:.4f} Ang^3 ({v0_ref / n_atoms:.4f} Ang^3/atom)", f_out)

    print_section("[3] VOLUME SCAN (MACE)", f_out)
    volumes, energies = run_volume_scan(atoms, calc, args, f_out)
    fit = fit_eos(volumes, energies, args.eos)

    plot_series = [("fine-tuned" if args.custom_model else args.model, volumes, energies, fit, _SERIES_COLORS[0])]

    fit_found = None
    if compare:
        print_section("[3b] FOUNDATION MODEL COMPARISON", f_out)
        calc_found = mace_relax.get_calculator(model=args.model, device=args.device)
        atoms_found = atoms.copy()
        atoms_found.calc = calc_found
        if args.pre_relax:
            mace_relax.relax(atoms_found, calc_found, cell_mask=[True] * 6,
                              fmax=args.fmax, max_steps=args.max_steps)
        volumes_found, energies_found = run_volume_scan(atoms_found, calc_found, args, f_out, tag=" (foundation)")
        fit_found = fit_eos(volumes_found, energies_found, args.eos)
        plot_series.append((f"foundation ({args.model})", volumes_found, energies_found, fit_found, _SERIES_COLORS[1]))

    print_section("[4] EQUATION-OF-STATE FIT", f_out)
    report_eos_fit(fit, n_atoms, "", f_out)
    if fit_found is not None:
        print_dual("", f_out)
        report_eos_fit(fit_found, n_atoms, " (foundation)", f_out)
        db0 = fit["b0_gpa"] - fit_found["b0_gpa"]
        print_dual(f"\nB0 difference (fine-tuned - foundation): {db0:+.2f} GPa", f_out)

    print_dual(f"\n[NOTE] Cross-check: compare Bulk modulus (B0) above against stb-mlelastic's own "
               f"'Bulk Modulus (B)' (Hill average, from the stress-strain elastic tensor) on the "
               f"same structure -- an independent-method agreement check on the MACE potential, "
               f"not a duplicate calculation (mlelastic strains the cell along canonical Voigt "
               f"directions; this tool scans the isotropic volume itself).", f_out)

    print_section("[5] SUMMARY & FILES", f_out)
    plot_path = os.path.join(args.output_dir, "eos_curve.png")
    plot_eos(plot_series, plot_path)
    print_dual(f"[OK] Wrote {plot_path}", f_out)

    if args.save_data:
        for label, vols, engs, _, _ in plot_series:
            safe_label = label.replace(" ", "_").replace("(", "").replace(")", "")
            dat_path = os.path.join(args.output_dir, f"eos_data_{safe_label}.dat")
            with open(dat_path, "w") as fh:
                fh.write("# Volume(Ang^3)  Energy(eV)\n")
                for v, e in zip(vols, engs):
                    fh.write(f"{v:.6f}  {e:.6f}\n")
            print_dual(f"[OK] Wrote {dat_path}", f_out)

    print_dual(f"\n{color_text('Done.', 'green')}", f_out)
    if f_out:
        f_out.close()


if __name__ == "__main__":
    main()
