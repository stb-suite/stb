#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""Aggregates a stb-eosInputs 'vol_*' sweep of real SIESTA calculations and
fits an equation of state (E vs V) -- the real-DFT analog of stb-mleos,
which does the exact same fit but on MACE-computed (volume, energy) pairs
instead of real SIESTA output files. Reuses the fitting math itself
(core/eos_fit.py's fit_eos/normalize_eos_string) rather than duplicating it,
since a Birch-Murnaghan/Vinet/etc. curve fit is the same regardless of
whether the underlying energies came from SIESTA or MACE; only the
data-mining step (reading each folder's calc.out here, vs. an in-memory
MACE evaluation loop there) differs.

Volume comes directly from each folder's own SIESTA output (core/
siesta_log.py's get_outcell, the LAST 'outcell:' block in the file) rather
than re-deriving it from the folder's input .fdf -- one fewer file
dependency, and the authoritative source for whatever cell SIESTA actually
used for that run (matches stb-strainAnalysis's own convention of reading
the cell straight out of the .out file).

Deliberately distinct from stb-elasticAnalysis's own bulk modulus: that one
comes from a linear stress-strain fit around a single reference cell (the
elastic STIFFNESS TENSOR, real SIESTA stress tensors); this one comes from
the curvature of the total-energy-vs-volume curve across several separate
SIESTA calculations. The two are independent methods on independent SIESTA
data -- agreement between them is a genuine physical consistency check on
the DFT setup itself (pseudopotentials, basis, mesh cutoff), not a
duplicate calculation, and the report prints an explicit [NOTE] pointing at
this cross-check (plus stb-mleos, the MACE analog of THIS SAME method, for
a DFT-vs-ML sanity check on top).

Plots use the same gnuplot .dat+.gplot pair convention as this workflow
category's siblings (stb-strainAnalysis/stb-elasticAnalysis/stb-
convergenceAnalysis), not the newer matplotlib convention used by the ML
Simulations tools -- consistent with every other WORKFLOW_TOOLS analysis
stage.

Three features on top of the original single-form fit:
- --eos all: fits every supported EOS form (birchmurnaghan/vinet/murnaghan/sj) on
  the exact same (volume, energy) data and reports V0/E0/B0/B0'/R^2 side by side --
  a robustness check on the fit itself (do different functional forms agree?), not
  a second data-mining pass. Mutually exclusive with --target-pressure (which needs
  one specific form to invert).
- --target-pressure: inverts the fitted EOS (P = -dE/dV, via a numerical derivative
  of the already-fitted curve, core/eos_fit.py's invert_pressure) to report the
  predicted equilibrium volume and linear lattice-parameter scale factor at one or
  more target external pressures (GPa) -- flags [EXTRAPOLATED] whenever the
  requested pressure falls outside the range the scan actually achieved, since that
  prediction then relies on the fitted form's shape well beyond where any real data
  constrained it.
- v0-outside-scanned-range warning: core/eos_fit.py's fit_eos now calls
  eos.fit(warn=False) and checks this explicitly itself (v0_outside_range in the
  returned dict) instead of relying on ase.eos's own warnings.warn() call, which
  bypassed every caller's print_dual-formatted report entirely. Fixes the same
  latent gap in stb-mleos too, since both tools share this fit_eos().

invert_pressure was verified against the analytic Birch-Murnaghan pressure formula
before being wired into the CLI: fit a synthetic exact BM curve (E0=-300 eV, B0=90
GPa, B0'=4.3, V0=160.103 Ang^3), inverted at several target pressures, then
evaluated the analytic P(V) formula at each returned volume -- recovered the
requested pressure to within ~0.0001 GPa in every case (e.g. requesting 5.0 GPa
returned V=152.311 Ang^3, and the analytic formula at that volume gives 4.9997 GPa).
Also confirmed a wildly out-of-range request (50 GPa, when the scan only spans
about -4.6 to +6.4 GPa) is correctly flagged [EXTRAPOLATED] -- catching a real bug
found while designing this check: np.interp CLAMPS to the nearest volume at the
scan boundary for an out-of-range pressure rather than extrapolating, so checking
the RETURNED VOLUME against the scanned range (the first, wrong implementation)
never flags anything; the fix checks the REQUESTED PRESSURE against the range of
pressures the fitted curve actually spans instead.
"""

VERSION = "1.1.0"  # --eos all, --target-pressure, v0-outside-range warning (core/eos_fit.py)

import os
import re
import sys
import glob
import argparse
from datetime import datetime

import numpy as np

from stb.core import structure_io, siesta_log
from stb.core.cli import color_text, show_intro, print_dual, print_section
from stb.core.eos_fit import normalize_eos_string, fit_eos, invert_pressure

REPORT_FILE = "stb_eosAnalysis_report.txt"
_ALL_EOS_FORMS = ["birchmurnaghan", "vinet", "murnaghan", "sj"]


def write_eos_plot(dat_path, gplot_path, volumes, energies, fit):
    """Writes <dat_path> (raw scanned points + the fitted curve, in separate
    column blocks so gnuplot can plot the data as points and the curve as a
    line from the same file) plus a companion .gplot, same convention as
    stb-elasticAnalysis's write_stress_plots/stb-convergenceAnalysis's
    write_curve_plot.
    """
    with open(dat_path, 'w') as f:
        f.write("# Equation of state | Data block, then a blank line, then the fitted curve\n")
        f.write("# 1:Volume(Ang^3) 2:Energy(eV)\n")
        for v, e in zip(volumes, energies):
            f.write(f"{v:.6f}  {e:.6f}\n")
        f.write("\n\n")
        f.write("# Fitted curve\n")
        for v, e in zip(fit["curve_v"], fit["curve_e"]):
            f.write(f"{v:.6f}  {e:.6f}\n")

    dat_name = os.path.basename(dat_path)
    base_name = os.path.splitext(os.path.basename(gplot_path))[0]
    with open(gplot_path, 'w') as f:
        f.writelines([
            '# --- STB Plot Configuration ---\n',
            '# Generated by stb-eosAnalysis\n',
            'set terminal pdfcairo enhanced color font "Arial,14" size 7,5\n',
            f'set output "{base_name}.pdf"\n\n',
            'set title "Equation of State (SIESTA)"\n',
            'set xlabel "Volume (Ang^3)"\n',
            'set ylabel "Energy (eV)"\n',
            'set grid\n',
            'set key top center\n',
            (f'plot "{dat_name}" index 0 using 1:2 with points pt 7 ps 1.5 lc rgb "#2255cc" '
             f'title "SIESTA data", "{dat_name}" index 1 using 1:2 with lines lw 2 '
             f'lc rgb "#cc5522" title "EOS fit"\n'),
        ])


def main():
    parser = argparse.ArgumentParser(
        description="Aggregates a stb-eosInputs 'vol_*' sweep of real SIESTA calculations "
                     "and fits an equation of state (Birch-Murnaghan by default) to get "
                     "V0/E0/B0 (bulk modulus)/B0'.",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Example usage:\n"
               "  stb-eosAnalysis --dir eos_runs --file calc.out\n"
               "  stb-eosAnalysis --dir eos_runs --eos vinet\n",
    )
    parser.add_argument("--dir", default="eos_runs",
                        help="Directory containing the 'vol_*' run folders (default: eos_runs).")
    parser.add_argument("--file", default="calc.out",
                        help="SIESTA output filename inside each run folder (default: calc.out).")
    parser.add_argument("--eos",
                        choices=["birchmurnaghan", "birch_murnaghan", "vinet", "murnaghan", "sj", "all"],
                        default="birchmurnaghan",
                        help="Equation-of-state form to fit (default: birchmurnaghan). 'sj' = "
                             "stabilized jellium (no B0' reported). 'all' fits every form on "
                             "the same data and reports them side by side (the plotted curve "
                             "still uses birchmurnaghan); mutually exclusive with "
                             "--target-pressure.")
    parser.add_argument("--target-pressure", type=float, nargs='+', default=None, metavar="GPA",
                        help="One or more target external pressures (GPa) -- inverts the "
                             "fitted EOS (P = -dE/dV) to report the predicted equilibrium "
                             "volume and linear lattice-parameter scale factor at each one. "
                             "Not available with --eos all (pick one specific form first).")
    parser.add_argument("-o", "--output", default="eos_curve",
                        help="Base filename (no extension) for the .dat/.gplot outputs "
                             "(default: eos_curve).")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-eosAnalysis {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()
    args.eos = normalize_eos_string(args.eos)

    if args.eos == "all" and args.target_pressure:
        parser.error("--target-pressure needs one specific --eos form, not 'all'.")

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("EQUATION-OF-STATE ANALYSIS:", 'bold'))
    print("-" * 60)

    if not os.path.isdir(args.dir):
        sys.exit(color_text(f"[ERROR] Directory '{args.dir}' not found.", 'red'))

    folders = sorted(
        f for f in os.listdir(args.dir)
        if os.path.isdir(os.path.join(args.dir, f)) and f.startswith("vol_")
    )
    if not folders:
        sys.exit(color_text(f"[ERROR] No 'vol_*' folders found in '{args.dir}'.", 'red'))

    report_path = REPORT_FILE if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    print_dual(color_text("===== STB-EOSANALYSIS REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time  : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Directory  : {args.dir}", f_out)
    print_dual(f"Output file: {args.file}", f_out)
    print_dual(f"EOS form   : {args.eos}", f_out)

    print_section("[1] READING FOLDERS", f_out)
    rows = []  # (volume, energy, n_atoms_or_None, scf_converged, max_force)
    n_skipped = 0
    for folder in folders:
        folder_path = os.path.join(args.dir, folder)
        out_path = os.path.join(folder_path, args.file)

        if not os.path.isfile(out_path):
            n_skipped += 1
            print_dual(f"   -> {folder:<20} : {color_text('SKIP', 'yellow')} "
                        f"(missing {args.file})", f_out)
            continue

        energy = siesta_log.get_free_energy(out_path)
        cell = siesta_log.get_outcell(out_path)
        if energy is None or cell is None:
            n_skipped += 1
            print_dual(f"   -> {folder:<20} : {color_text('SKIP', 'yellow')} "
                        f"(could not parse energy/cell)", f_out)
            continue

        volume = abs(np.linalg.det(cell))
        scf_converged, _iterations = siesta_log.get_scf_convergence(out_path)
        max_force = siesta_log.get_max_force(out_path)

        n_atoms = None
        fdf_candidates = glob.glob(os.path.join(folder_path, "*.fdf"))
        if len(fdf_candidates) == 1:
            try:
                n_atoms = len(structure_io.read_fdf(fdf_candidates[0]).atoms)
            except Exception:
                n_atoms = None

        rows.append((volume, energy, n_atoms, scf_converged, max_force))
        status = color_text('OK', 'green') if scf_converged else \
            color_text('OK (SCF not confirmed converged)', 'yellow')
        print_dual(f"   -> {folder:<20} : {status}", f_out)

    if not rows:
        print_dual(color_text("\n[FAIL] No valid data found in any 'vol_*' folder.", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    rows.sort(key=lambda r: r[0])
    atom_counts = {r[2] for r in rows if r[2] is not None}
    if len(atom_counts) > 1:
        print_dual(color_text(
            f"[WARNING] Folders report different atom counts ({sorted(atom_counts)}) -- "
            f"make sure '{args.dir}' only contains one stb-eosInputs sweep.", 'yellow'), f_out)

    print_dual(f"\n  Runs found: {len(rows)} (skipped: {n_skipped})", f_out)

    print_section("[2] VOLUME-ENERGY TABLE", f_out)
    header = f"{'Volume(Ang^3)':>14} {'Energy(eV)':>14} {'E/atom(eV)':>14} {'SCF':<6} {'MaxForce(eV/Ang)':>16}"
    print_dual(header, f_out)
    print_dual("-" * len(header), f_out)
    for volume, energy, n_atoms, scf_converged, max_force in rows:
        e_per_atom_str = f"{energy / n_atoms:14.6f}" if n_atoms else f"{'--':>14}"
        scf_str = "OK" if scf_converged else color_text("WARN", 'yellow')
        force_str = f"{max_force:.6f}" if max_force is not None else "--"
        print_dual(f"{volume:14.4f} {energy:14.6f} {e_per_atom_str} {scf_str:<6} {force_str:>16}", f_out)
    print_dual("-" * len(header), f_out)

    volumes = np.array([r[0] for r in rows])
    energies = np.array([r[1] for r in rows])

    if len(volumes) < 4:
        print_dual(color_text(
            f"\n[FAIL] Only {len(volumes)} valid volume(s) found -- at least 4 are needed "
            "for a 4-parameter EOS fit (5+ recommended for a robust one).", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    n_atoms_ref = rows[0][2]

    if args.eos == "all":
        print_section("[3] EQUATION-OF-STATE FIT (ALL FORMS)", f_out)
        header_all = (f"{'EOS form':<15} {'V0(Ang^3)':>11} {'E0(eV)':>14} "
                      f"{'B0(GPa)':>9} {'B0prime':>9} {'R^2':>10}")
        print_dual(header_all, f_out)
        print_dual("-" * len(header_all), f_out)
        fits_by_form = {}
        for form in _ALL_EOS_FORMS:
            fit_i = fit_eos(volumes, energies, form)
            fits_by_form[form] = fit_i
            bprime_str = f"{fit_i['bprime']:.3f}" if fit_i["bprime"] is not None else "--"
            flag = color_text(" [V0 OUTSIDE SCANNED RANGE]", 'yellow') if fit_i["v0_outside_range"] else ""
            print_dual(f"{form:<15} {fit_i['v0']:11.4f} {fit_i['e0']:14.6f} "
                        f"{fit_i['b0_gpa']:9.2f} {bprime_str:>9} {fit_i['r_squared']:10.6f}{flag}", f_out)
        print_dual("-" * len(header_all), f_out)
        fit = fits_by_form["birchmurnaghan"]
        print_dual(f"\n[INFO] The plotted curve below uses birchmurnaghan's fit -- the table "
                   f"above already compares every form on the exact same (volume, energy) data.", f_out)
    else:
        print_section("[3] EQUATION-OF-STATE FIT", f_out)
        fit = fit_eos(volumes, energies, args.eos)

        v0_str = color_text(f"{fit['v0']:.4f}", 'bold')
        e0_str = color_text(f"{fit['e0']:.6f}", 'bold')
        b0_str = color_text(f"{fit['b0_gpa']:.2f}", 'bold')
        v0_per_atom = f" ({fit['v0'] / n_atoms_ref:.4f} Ang^3/atom)" if n_atoms_ref else ""
        e0_per_atom = f" ({fit['e0'] / n_atoms_ref:.6f} eV/atom)" if n_atoms_ref else ""
        print_dual(f"Equilibrium volume (V0): {v0_str} Ang^3{v0_per_atom}", f_out)
        print_dual(f"Equilibrium energy (E0): {e0_str} eV{e0_per_atom}", f_out)
        print_dual(f"Bulk modulus       (B0): {b0_str} GPa", f_out)
        if fit["bprime"] is not None:
            bprime_str = color_text(f"{fit['bprime']:.3f}", 'bold')
            print_dual(f"Pressure derivative (B0'): {bprime_str}", f_out)
        print_dual(f"Fit quality        (R^2): {fit['r_squared']:.6f}", f_out)
        if fit["r_squared"] < 0.99:
            print_dual(color_text(
                "[WARNING] R^2 < 0.99 -- the scanned strain range may be too wide (leaving "
                "the harmonic/near-harmonic regime the EOS form assumes) or too narrow "
                "(too little curvature to constrain the fit); consider re-running "
                "stb-eosInputs with a different --strain-range.", 'yellow'), f_out)
        if fit["v0_outside_range"]:
            print_dual(color_text(
                "[WARNING] The fitted equilibrium volume (V0) falls OUTSIDE the range of "
                "scanned volumes -- the scan didn't actually bracket the true minimum; "
                "widen --strain-range in stb-eosInputs and re-run.", 'yellow'), f_out)

        if args.target_pressure:
            print_section("[3b] TARGET PRESSURE", f_out)
            header_p = f"{'P(GPa)':>10} {'V(Ang^3)':>12} {'a/a0 (linear)':>14}"
            print_dual(header_p, f_out)
            print_dual("-" * len(header_p), f_out)
            any_extrapolated = False
            for p in args.target_pressure:
                v_at_p, extrapolated = invert_pressure(fit, p)
                any_extrapolated = any_extrapolated or extrapolated
                scale = (v_at_p / fit["v0"]) ** (1.0 / 3.0)
                flag = color_text(" [EXTRAPOLATED]", 'yellow') if extrapolated else ""
                print_dual(f"{p:10.3f} {v_at_p:12.4f} {scale:14.6f}{flag}", f_out)
            print_dual("-" * len(header_p), f_out)
            if any_extrapolated:
                print_dual(color_text(
                    "[WARNING] One or more requested pressures fall outside the range "
                    "actually achieved across the scanned volumes -- those predictions "
                    "rely on the fitted EOS form's shape well beyond where any data "
                    "constrained it; widen --strain-range in stb-eosInputs and re-run "
                    "for a more reliable prediction.", 'yellow'), f_out)

    print_dual(f"\n[NOTE] Cross-check: compare Bulk modulus (B0) above against "
               f"stb-elasticAnalysis's own 'Bulk Modulus (B)' (Hill average, from the "
               f"stress-strain elastic tensor) on the same structure -- an independent "
               f"-method agreement check on the same DFT setup, not a duplicate "
               f"calculation. Also compare against stb-mleos (the MACE analog of THIS "
               f"SAME curvature-based method) for a DFT-vs-ML sanity check.", f_out)

    print_section("[4] SUMMARY & FILES", f_out)
    dat_path = f"{args.output}.dat"
    gplot_path = f"{args.output}.gplot"
    write_eos_plot(dat_path, gplot_path, volumes, energies, fit)
    print_dual(f"[OK] Wrote {dat_path}, {gplot_path} "
                f"(gnuplot {gplot_path} to render {args.output}.pdf)", f_out)
    if report_path:
        print_dual(f"[OK] Report saved to: {report_path}", f_out)

    if f_out:
        f_out.close()


if __name__ == "__main__":
    main()
