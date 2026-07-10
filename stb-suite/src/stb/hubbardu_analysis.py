#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.0.0"

import os
import re
import sys
import json
import argparse
import numpy as np
from stb.core.cli import COLORS, color_text, show_intro

# Reference GGA+U values (eV), for comparison ONLY -- never substituted for
# the linear-response result computed by this tool. Materials Project's
# oxide calibration (docs.materialsproject.org/methodology/materials-methodology/
# calculation-details/gga+u-calculations/hubbard-u-values), itself derived
# from Wang, Maxisch & Ceder, Phys. Rev. B 73, 195107 (2006). Calibrated for
# VASP GGA+U on oxides -- a starting-point sanity check, not a validation:
# the actual U depends on functional, pseudopotential, and basis set.
REFERENCE_U = {
    "V": 3.25, "Cr": 3.7, "Mn": 3.9, "Fe": 5.3,
    "Co": 3.32, "Ni": 6.2, "Mo": 4.38, "W": 6.2,
}

_OCC_RE = re.compile(r'Occupations:\s*([0-9eE+.\-\s]+)')


def parse_occupation(out_path):
    """Last 'Occupations:' line's LAST numeric value (the total shell
    occupation, summed over spin) from a SIESTA .out file produced with
    LDAU.PotentialShift T. Robust to both non-polarized (2 numbers: single
    spin channel + total) and spin-polarized (3 numbers: up, down, total)
    output, since the total is always the last value either way (verified
    against SIESTA's own source, Src/dftu.F, format '(a,/,a,3f12.6)').
    Returns None on any parse error.
    """
    last = None
    try:
        with open(out_path, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                match = _OCC_RE.search(line)
                if match:
                    nums = match.group(1).split()
                    if nums:
                        try:
                            last = float(nums[-1])
                        except ValueError:
                            pass
    except Exception:
        return None
    return last


def linear_fit_r2(alphas, occs):
    """np.polyfit degree-1 fit, returning (slope, intercept, r_squared)."""
    alphas = np.array(alphas)
    occs = np.array(occs)
    slope, intercept = np.polyfit(alphas, occs, 1)
    predicted = slope * alphas + intercept
    ss_res = np.sum((occs - predicted) ** 2)
    ss_tot = np.sum((occs - np.mean(occs)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return float(slope), float(intercept), float(r2)


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Analyzes a stb-hubbardu linear-response sweep and computes the Hubbard U.", 'bold')}
Fits the self-consistent (screened) and frozen-density (bare) occupation
responses to the perturbation strength, then applies the standard
Cococcioni & de Gironcoli formula: U = 1/chi0 - 1/chi (Phys. Rev. B 71,
035105, 2005). Writes a ready-to-use %%block LDAU.proj snippet with the
computed U.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage example:\n"
               "  %(prog)s --dir hubbardu_runs --file calc.out\n"
    )

    parser.add_argument("--dir", type=str, default="hubbardu_runs",
                        help="Directory containing the run folders and run_manifest.json "
                             "written by stb-hubbardu (default: hubbardu_runs).")
    parser.add_argument("--file", type=str, default="calc.out",
                        help="SIESTA output filename inside each run folder (default: calc.out).")
    parser.add_argument("--r2-tolerance", type=float, default=0.98,
                        help="Warn if either fit's R^2 falls below this (default: 0.98) -- a poor "
                             "linear fit means the perturbation range is probably too wide.")
    parser.add_argument("-o", "--output", type=str, default=None,
                        help="Output .fdf snippet filename (default: <species>_LDAU.fdf).")
    parser.add_argument("-v", "--version", action="version", version=f"stb-hubbarduAnalysis {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("Analyze a Hubbard U linear-response sweep:", 'bold'))
    print("-" * 60)

    if not os.path.isdir(args.dir):
        print(color_text(f"Error: Directory '{args.dir}' not found.", 'red'))
        sys.exit(1)

    manifest_path = os.path.join(args.dir, "run_manifest.json")
    if not os.path.isfile(manifest_path):
        print(color_text(
            f"Error: 'run_manifest.json' not found in '{args.dir}' -- this directory "
            "wasn't generated by stb-hubbardu.", 'red'))
        sys.exit(1)

    with open(manifest_path, 'r') as f:
        manifest = json.load(f)
    species, n, l = manifest["species"], manifest["n"], manifest["l"]

    scf_points, frozen_points = [], []
    print(f"\n{color_text('READING RUNS:', 'bold')}")
    for folder, info in manifest["runs"].items():
        out_path = os.path.join(args.dir, folder, args.file)
        occ = parse_occupation(out_path)
        if occ is None:
            print(f"   -> {folder:<24} : {color_text('SKIP', 'yellow')} (missing or unparseable {args.file})")
            continue

        kind, alpha = info["kind"], info["alpha"]
        if kind in ("scf", "reference"):
            scf_points.append((alpha, occ))
        if kind in ("frozen", "reference"):
            frozen_points.append((alpha, occ))
        print(f"   -> {folder:<24} : {color_text('OK', 'green')} (occupation={occ:.6f})")

    if len(scf_points) < 2 or len(frozen_points) < 2:
        print(color_text(
            "\nError: Not enough valid runs to fit a response (need at least 2 scf and "
            "2 frozen points, including the reference).", 'red'))
        sys.exit(1)

    scf_points.sort()
    frozen_points.sort()
    chi, chi_intercept, chi_r2 = linear_fit_r2(*zip(*scf_points))
    chi0, chi0_intercept, chi0_r2 = linear_fit_r2(*zip(*frozen_points))

    print(f"\n{color_text('Self-consistent (screened) response chi:', 'cyan')} "
          f"{chi:.6f} 1/eV  (R^2={chi_r2:.4f})")
    print(f"{color_text('Frozen-density (bare) response chi0:', 'cyan')} "
          f"{chi0:.6f} 1/eV  (R^2={chi0_r2:.4f})")

    for name, r2 in (("self-consistent", chi_r2), ("frozen-density", chi0_r2)):
        if r2 < args.r2_tolerance:
            print(color_text(
                f"[WARNING] The {name} response fit has R^2={r2:.4f}, below the "
                f"{args.r2_tolerance} tolerance -- the perturbation range may be too wide "
                "for the linear regime. Consider narrower --alphas in stb-hubbardu.", 'yellow'))

    if chi == 0 or chi0 == 0:
        print(color_text(
            "\nError: A zero response slope was fitted -- cannot compute U (division by "
            "zero). Check that the occupation actually responds to the perturbation.", 'red'))
        sys.exit(1)

    U = 1.0 / chi0 - 1.0 / chi
    print(f"\n{color_text('Computed U (Cococcioni & de Gironcoli, PRB 71, 035105, 2005):', 'bold')} "
          f"{color_text(f'{U:.4f} eV', 'green')}")

    if species in REFERENCE_U:
        print(color_text(
            f"[INFO] Literature reference (GGA+U, oxides, Wang/Maxisch/Ceder PRB 73, 195107 "
            f"(2006)) for {species}: {REFERENCE_U[species]:.2f} eV -- a sanity-check comparison "
            "only, not a validation (functional/pseudopotential/basis-dependent).", 'cyan'))

    output = args.output or f"{species}_LDAU.fdf"
    with open(output, 'w') as f:
        f.write("%block LDAU.proj\n")
        f.write(f"{species}   1\n")
        f.write(f"n={n}    {l}\n")
        f.write(f"{U:.3f}    0.000\n")
        f.write("0.000    0.000\n")
        f.write("%endblock LDAU.proj\n")

    print(f"\n{color_text('[Saved]', 'cyan')} Ready-to-use DFT+U block -> {output}")


if __name__ == "__main__":
    main()
