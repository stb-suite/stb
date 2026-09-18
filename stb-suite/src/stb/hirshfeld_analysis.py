#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""Stage 3 (final) of the Hirshfeld-I (iterative) charge-partitioning
workflow. Reads back every ion/<species>/cation//anion/ folder
stb-hirshfeldIons wrote (now assumed already run through SIESTA),
spherically-averages each into a radial density profile, then runs the
actual Hirshfeld-I convergence loop (core.hirshfeld.iterate_hirshfeld_i)
against the combined/production system's fixed real density: following
the original iterative-Hirshfeld formulation (Bultinck et al., J. Chem.
Phys. 126, 144111 (2007)) literally, EACH ATOM independently interpolates
its own reference between its species' neutral profile and whichever of
that species' cation/anion profile matches THAT ATOM's own charge sign
from the previous round -- two atoms of the same species can genuinely
lean toward opposite ion states in the same round. Iterates until
max|Delta q| converges.

This is the only stage that iterates -- no new folders/geometry are
written here, matching the suite's general "an Analysis stage never
generates new geometry" rule. Reports both the iteration-0 (plain simple
Hirshfeld -- directly comparable to stb-nativecharges'/SIESTA's own native
Hirshfeld output for the same system) and the converged Hirshfeld-I result
side by side, plus the convergence history, so the shift between the two
methods is visible rather than only the final number.
"""

VERSION = "2.0.0"  # reads BOTH cation and anion per species and lets each atom
                    # pick its own, per the literal literature formulation
                    # (previously blended toward a single per-species ion sign)

import os
import sys
import time
import json
import argparse
import statistics
from datetime import datetime
from collections import defaultdict

import numpy as np

from stb.core import structure_io, kspace, hirshfeld as hf, siesta_log
from stb.core.deps import require_sisl
from stb.core.cli import color_text, show_intro, print_dual, print_section, print_table
from stb.core.rho_io import find_one_rho, read_total

sisl = require_sisl()

HIRSHFELD_IONS_MANIFEST_FILE = "hirshfeld_ions_manifest.json"
REPORT_FILE = "hirshfeld_analysis_report.txt"
N_RADIAL_BINS = 200
DEFAULT_TOL = 0.005
DEFAULT_MAX_ITER = 20


def load_species_profile(folder, n_bins=N_RADIAL_BINS):
    """Same as hirshfeld_ions.py::load_neutral_profile, generalized to any
    single-isolated-atom folder (neutral or ion) -- reads its own
    structure.fdf for the exact box/center rather than trusting a manifest
    value, then spherically-averages its '*.RHO'.
    """
    rho_path = find_one_rho(folder)
    sile = sisl.get_sile(rho_path)
    grid = read_total(sile)
    structure = structure_io.read_fdf(os.path.join(folder, "structure.fdf"))
    center_frac = structure.atoms[0][1]
    if structure.coord_format == "cartesian":
        center_frac = kspace.to_fractional(center_frac, structure.lattice, True)
    return hf.radial_density_profile(grid.grid, structure.lattice, center_frac, n_bins=n_bins)


def main():
    parser = argparse.ArgumentParser(
        description="Stage 3 of the Hirshfeld-I workflow: reads back every neutral+ion "
                     "isolated-atom folder, then iterates the Hirshfeld-I reference-density "
                     "refinement against the combined/production system's real density until "
                     "convergence.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("-O", "--output-dir", default="hirshfeld_study",
                         help="Same root directory stb-hirshfeldPrep/Ions wrote into "
                              "(default: hirshfeld_study).")
    parser.add_argument("--tol", type=float, default=DEFAULT_TOL, metavar="E-",
                         help=f"Convergence tolerance on max|Delta q| between rounds, in "
                              f"electrons (default: {DEFAULT_TOL:g} e-, matching the method's "
                              f"own prototype).")
    parser.add_argument("--max-iter", type=int, default=DEFAULT_MAX_ITER, metavar="N",
                         help=f"Maximum number of refinement rounds (default: "
                              f"{DEFAULT_MAX_ITER}; the method's own prototype converged in "
                              "8).")
    parser.add_argument("--ref", required=False, default=None,
                         help="Path to a specific SIESTA output file to read each species' "
                              "Z_val (valence charge) from -- same convention as stb-bader's "
                              "own --ref. Overrides auto-detection of combined/'s own '*.out' "
                              "(which assumes that log was left with a .out extension; use "
                              "--ref if yours wasn't).")
    parser.add_argument("--save-report", action="store_true",
                         help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--no-intro", dest="intro", action="store_false",
                         help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version",
                         version=f"stb-hirshfeldAnalysis {VERSION}")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - Hirshfeld-I Analysis (Stage 3/3)",
            "Iterative Hirshfeld-I atomic charges",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    report_path = REPORT_FILE if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    print_dual(color_text("===== STB-HIRSHFELDANALYSIS REPORT =====", 'magenta'), f_out)
    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time      : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Study dir      : {args.output_dir}", f_out)
    print_dual(f"Tolerance      : {args.tol:g} e-", f_out)
    print_dual(f"Max iterations : {args.max_iter}", f_out)

    ions_manifest_path = os.path.join(args.output_dir, HIRSHFELD_IONS_MANIFEST_FILE)
    if not os.path.isfile(ions_manifest_path):
        print_dual(color_text(
            f"[ERROR] '{ions_manifest_path}' not found -- run stb-hirshfeldPrep then "
            "stb-hirshfeldIons first.", 'red'), f_out)
        sys.exit(1)
    with open(ions_manifest_path) as f:
        ions_manifest = json.load(f)
    with open(ions_manifest["hirshfeld_manifest"]) as f:
        manifest = json.load(f)

    try:
        structure = structure_io.read_fdf(manifest["combined_structure"])
    except Exception as e:
        print_dual(color_text(
            f"[ERROR] Could not read the combined structure "
            f"'{manifest['combined_structure']}': {e}", 'red'), f_out)
        sys.exit(1)

    try:
        combined_rho_path = find_one_rho(manifest["combined_dir"])
    except ValueError as e:
        print_dual(color_text(
            f"[ERROR] {e} -- has SIESTA been run in 'combined/' yet?", 'red'), f_out)
        sys.exit(1)
    try:
        combined_sile = sisl.get_sile(combined_rho_path)
        rho_real = read_total(combined_sile).grid
    except Exception as e:
        print_dual(color_text(
            f"[ERROR] Could not read the combined '.RHO' '{combined_rho_path}': {e}", 'red'), f_out)
        sys.exit(1)

    symbols = [sym for sym, _ in structure.atoms]
    positions = np.array([pos for _, pos in structure.atoms])
    frac_positions = kspace.to_fractional(positions, structure.lattice,
                                           structure.coord_format == "cartesian")

    # Z_val (VALENCE charge) -- NOT the atomic number Z used to pick the
    # pseudopotential -- same detected-from-.out-with-tabulated-fallback
    # resolution as stb-bader's own valence_source (see
    # hirshfeld_ions.py's identical block for why this distinction matters:
    # using Z instead of Z_val silently produces nonsense charges). --ref
    # (same flag/behavior as stb-bader's own) takes priority over auto-
    # detecting combined/'s own '*.out'.
    combined_out = args.ref or siesta_log.find_out_file(manifest["combined_dir"], None)
    detected_valence = siesta_log.get_zval_from_output(None, override_path=combined_out) \
        if combined_out else None
    valence_source = {**siesta_log.FALLBACK_VALENCE, **(detected_valence or {})}
    unknown_syms = sorted({sym for sym in symbols if sym not in valence_source})
    if unknown_syms:
        print_dual(color_text(
            f"[ERROR] No valence Z_val (detected or fallback) for species: "
            f"{', '.join(unknown_syms)}.", 'red'), f_out)
        sys.exit(1)
    z_vals = [valence_source[sym] for sym in symbols]
    fallback_syms = sorted(set(symbols) - set((detected_valence or {}).keys()))
    if fallback_syms:
        print_dual(color_text(
            f"[INFO] Z_val for {', '.join(fallback_syms)} taken from the hardcoded fallback "
            "table (not detected in combined/'s own .out) -- verify these match your "
            "pseudopotential's real valence charge if unsure.", 'cyan'), f_out)

    print_section("[1] NEUTRAL + CATION + ANION PROFILES", f_out)
    neutral_profiles, cation_profiles, anion_profiles = {}, {}, {}
    for entry in manifest["species"]:
        sym = entry["symbol"]
        neutral_folder = os.path.join(args.output_dir, entry["neutral_folder"])
        try:
            neutral_profiles[sym] = load_species_profile(neutral_folder)
            print_dual(f"  {sym} (neutral): read from '{neutral_folder}'", f_out)
        except Exception as e:
            print_dual(color_text(
                f"[ERROR] Could not build a neutral profile for species '{sym}' from "
                f"'{neutral_folder}': {e}", 'red'), f_out)
            sys.exit(1)
    for entry in ions_manifest["species"]:
        sym = entry["symbol"]
        cation_folder = os.path.join(args.output_dir, entry["cation_folder"])
        anion_folder = os.path.join(args.output_dir, entry["anion_folder"])
        try:
            cation_profiles[sym] = load_species_profile(cation_folder)
            print_dual(f"  {sym} (cation) : read from '{cation_folder}'", f_out)
        except Exception as e:
            print_dual(color_text(
                f"[ERROR] Could not build a cation profile for species '{sym}' from "
                f"'{cation_folder}': {e} -- has SIESTA been run there yet?", 'red'), f_out)
            sys.exit(1)
        try:
            anion_profiles[sym] = load_species_profile(anion_folder)
            print_dual(f"  {sym} (anion)  : read from '{anion_folder}'", f_out)
        except Exception as e:
            print_dual(color_text(
                f"[ERROR] Could not build an anion profile for species '{sym}' from "
                f"'{anion_folder}': {e} -- has SIESTA been run there yet?", 'red'), f_out)
            sys.exit(1)

    species_profiles_0 = {i: neutral_profiles[symbols[i]] for i in range(len(symbols))}
    charges_0, populations_0 = hf.compute_hirshfeld_charges(
        rho_real, structure.lattice, frac_positions, z_vals, species_profiles_0,
        progress_label="Pass 0")

    print_section("[2] ITERATING HIRSHFELD-I", f_out)
    iter_start = time.monotonic()

    def report_iteration(iteration, max_delta, worst_atom_idx, charges, populations):
        elapsed = time.monotonic() - iter_start
        sym = symbols[worst_atom_idx]
        leaning = "cation-like" if charges[worst_atom_idx] >= 0 else "anion-like"
        status = color_text("CONVERGED", 'green') if max_delta < args.tol else "not yet converged"
        print_dual(
            f"  Iteration {iteration:>2}/{args.max_iter}: max|Delta q| = {max_delta:.6f} e- "
            f"-- atom #{worst_atom_idx + 1} ({sym}), now {charges[worst_atom_idx]:+.4f} e- "
            f"({leaning}, pop. {populations[worst_atom_idx]:.4f} e-), {elapsed:6.1f}s elapsed "
            f"-- {status}",
            f_out)

    charges_final, populations_final, history = hf.iterate_hirshfeld_i(
        rho_real, structure.lattice, frac_positions, symbols, z_vals,
        neutral_profiles, cation_profiles, anion_profiles, tol=args.tol, max_iter=args.max_iter,
        show_progress=True, initial_charges=charges_0, initial_populations=populations_0,
        on_iteration=report_iteration)

    print_section("[3] CONVERGENCE HISTORY", f_out)
    conv_rows = []
    for iteration, max_delta in history:
        delta_str = f"{max_delta:.6f}" if max_delta is not None else "N/A (pass 0)"
        conv_rows.append(([str(iteration), delta_str], None))
    print_table(["Iteration", "max|Delta q| (e-)"], conv_rows, f_out)
    converged = history[-1][1] is not None and history[-1][1] < args.tol
    if converged:
        print_dual(color_text(
            f"Converged after {history[-1][0]} iteration(s) (max|Delta q| < {args.tol:g} e-).",
            'green'), f_out)
    else:
        print_dual(color_text(
            f"[WARNING] Did not confirm convergence within --max-iter {args.max_iter} "
            f"iterations (last max|Delta q| = {history[-1][1]:.6f} e- >= {args.tol:g} e-). "
            "Re-run with a larger --max-iter, or treat the result below as unconverged.",
            'yellow'), f_out)

    print_section("[4] PER-ATOM CHARGES: SIMPLE HIRSHFELD vs. HIRSHFELD-I", f_out)
    headers = ["Idx", "Elem", "Simple Hirshfeld(e-)", "Hirshfeld-I(e-)", "Shift(e-)"]
    rows = []
    for i in range(len(symbols)):
        shift = charges_final[i] - charges_0[i]
        rows.append(([str(i + 1), symbols[i], f"{charges_0[i]:+.4f}",
                       f"{charges_final[i]:+.4f}", f"{shift:+.4f}"], None))
    print_table(headers, rows, f_out)

    print_section("[5] PER-SPECIES SUMMARY (HIRSHFELD-I)", f_out)
    by_species = defaultdict(list)
    for i, sym in enumerate(symbols):
        by_species[sym].append(charges_final[i])
    species_rows = []
    any_genuinely_mixed = False
    for sym in sorted(by_species):
        values = by_species[sym]
        n_cation = sum(1 for c in values if c >= 0)
        n_anion = sum(1 for c in values if c < 0)
        is_mixed = n_cation > 0 and n_anion > 0
        any_genuinely_mixed = any_genuinely_mixed or is_mixed
        mean = statistics.mean(values)
        std = statistics.pstdev(values) if len(values) > 1 else 0.0
        species_rows.append(([sym, str(len(values)), str(n_cation), str(n_anion),
                              f"{mean:+.4f}", f"{std:.4f}"], 'cyan' if is_mixed else None))
    print_table(["Elem", "N", "Cation-like", "Anion-like", "Mean(e-)", "Std(e-)"],
                species_rows, f_out)
    if any_genuinely_mixed:
        print_dual(color_text(
            "Species highlighted above converged with atoms genuinely leaning toward BOTH "
            "ion states -- exactly the case a single-sign-per-species approach cannot "
            "represent correctly.", 'cyan'), f_out)

    print_section("[6] SUMMARY", f_out)
    total_elapsed = time.monotonic() - iter_start
    print_dual(f"Atoms analyzed  : {len(symbols)}", f_out)
    print_dual(f"Iterations run  : {history[-1][0]} (of --max-iter {args.max_iter})", f_out)
    print_dual(f"Converged       : {'yes' if converged else 'NO'}", f_out)
    print_dual(f"Iterating time  : {total_elapsed:.1f}s", f_out)
    print_dual(f"Total charge    : {float(np.sum(charges_final)):+.4f} e- "
               f"(should be close to 0 for a neutral combined system)", f_out)
    if report_path:
        print_dual(f"Report          : {report_path}", f_out)

    if f_out:
        f_out.close()

    print("\n" + "-" * 60)
    print(color_text("Hirshfeld-I Analysis (Stage 3) complete.\n", 'bold'))


if __name__ == "__main__":
    main()
