#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""Generates the folders for a real-SIESTA equation-of-state (E vs V) scan --
the DFT-workflow analog of stb-mleos, which does the exact same isotropic
-volume-scan idea but with a MACE potential in-memory instead of real SIESTA
calculations to run by hand. Same relationship stb-strain has to
stb-mlelastic's stress-strain method.

Writes one vol_<pct> folder per scanned volume, each containing a copy of
the input .fdf with ONLY %block LatticeVectors rewritten (isotropically
scaled by factor = (1 + strain_pct/100)**(1/3), the same convention
stb-mleos/stb-mlphonons's --qha already use, so a --strain-range value means
the same thing in all three tools) -- everything else (species, basis,
pseudopotentials, SCF block, and critically the atomic positions, which are
expected in FRACTIONAL coordinates so they scale for free with the cell)
preserved verbatim via core/structure_io.py's rewrite_fdf_lattice, same
pattern stb-strain itself uses. The user runs SIESTA in each folder (a
single-point or fixed-cell relaxation -- do NOT let the cell relax, or the
whole point of the scan is lost), then stb-eosAnalysis aggregates the
results.

Bulk (3D periodic) only -- an equation of state assumes a well-defined,
bounded cell volume; a vacuum-padded axis (slab/wire/molecule) has nothing
of the kind (same reasoning as stb-mleos/stb-mlmelting/stb-amorphize).
"""

VERSION = "1.0.0"

import os
import sys
import argparse

import numpy as np

from stb.core import structure_io, kspace
from stb.core.cli import color_text, show_intro


def main():
    parser = argparse.ArgumentParser(
        description="Generates a set of isotropically-scaled-volume SIESTA input folders "
                     "for a real-DFT equation-of-state (E vs V) scan -- run SIESTA in each "
                     "one (fixed cell, positions may relax), then aggregate with "
                     "stb-eosAnalysis. The real-SIESTA analog of stb-mleos's MACE-driven "
                     "in-memory scan.",
        epilog="Example usage:\n"
               "  stb-eosInputs --file structure.fdf\n"
               "  stb-eosInputs --file structure.fdf --n-volumes 11 --strain-range 8.0 "
               "-o eos_runs\n",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--file", required=True, help="Input structure (.fdf). Atomic "
                         "coordinates must be in fractional (ScaledCartesian/Fractional) "
                         "format, same requirement as stb-strain.")
    parser.add_argument("--n-volumes", type=int, default=9,
                        help="Number of volumes to scan (default: 9). Needs >= 5 for a "
                             "robust 4-parameter EOS fit downstream in stb-eosAnalysis.")
    parser.add_argument("--strain-range", type=float, default=5.0,
                        help="Max isotropic volume strain, %% (default: 5.0), spanning "
                             "[-strain-range, +strain-range] -- same meaning as stb-mleos's "
                             "own --strain-range.")
    parser.add_argument("-o", "--output-dir", default="eos_runs",
                        help="Directory to write the 'vol_<pct>' run folders into "
                             "(default: eos_runs).")
    parser.add_argument("--vacuum-gap", type=float, default=10.0,
                        help="Gap (Ang) used to detect vacuum-padded axes -- bulk-only, "
                             "rejects any vacuum-padded structure (default: 10.0).")
    parser.add_argument("-v", "--version", action="version", version=f"stb-eosInputs {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.n_volumes < 5:
        parser.error("--n-volumes needs at least 5 points for a robust downstream EOS fit.")

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("EQUATION-OF-STATE INPUTS:", 'bold'))
    print("-" * 60)

    if not os.path.isfile(args.file):
        sys.exit(color_text(f"[ERROR] Input file '{args.file}' not found.", 'red'))

    try:
        structure = structure_io.read_fdf(args.file)
    except Exception as e:
        sys.exit(color_text(f"[ERROR] Could not read '{args.file}': {e}", 'red'))

    lattice_vectors = structure_io.raw_lattice_vectors(structure)

    positions = np.array([pos for _, pos in structure.atoms])
    is_cartesian = structure.coord_format == 'cartesian'
    frac_coords = kspace.to_fractional(positions, structure.lattice, is_cartesian)
    vacuum_axes = kspace.detect_vacuum_axes(frac_coords, structure.lattice, args.vacuum_gap)
    print(f"[INFO] Detected dimensionality: {kspace.dimensionality_label(vacuum_axes)}")

    if any(vacuum_axes):
        sys.exit(color_text(
            "[ERROR] This tool is bulk (3D periodic) only -- a vacuum-padded axis was "
            "detected. An equation of state assumes a well-defined, bounded cell volume, "
            "which a slab/wire/molecule doesn't have.", 'red'))

    if is_cartesian:
        print(color_text(
            "[WARNING] Atomic coordinates in this file are Cartesian, not fractional -- "
            "scaling only the lattice vectors will NOT scale the atomic positions along "
            "with the cell, giving a physically wrong (fixed-Cartesian-position) volume "
            "scan. Convert to fractional coordinates first (e.g. via stb-translate).",
            'yellow'))

    strains_pct = np.linspace(-abs(args.strain_range), abs(args.strain_range), args.n_volumes)
    v0_input = abs(np.linalg.det(structure.lattice))

    print(f"[INFO] Volume range: {-abs(args.strain_range):.1f}% to "
          f"{abs(args.strain_range):.1f}% ({args.n_volumes} points, V0={v0_input:.4f} Ang^3)")
    print("[INFO] Generating FDF files with scaled volumes...")

    os.makedirs(args.output_dir, exist_ok=True)
    input_basename = os.path.basename(args.file)

    for strain_pct in strains_pct:
        factor = (1.0 + strain_pct / 100.0) ** (1.0 / 3.0)
        new_vectors = lattice_vectors * factor

        strain_prefix = "m" if strain_pct < 0 else ""
        folder = os.path.join(args.output_dir, f"vol_{strain_prefix}{abs(strain_pct):.2f}")
        os.makedirs(folder, exist_ok=True)

        output_fdf = os.path.join(folder, input_basename)
        structure_io.rewrite_fdf_lattice(args.file, new_vectors, output_fdf)
        scaled_volume = v0_input * factor ** 3
        print(f"[OK] Generated: {output_fdf} (strain: {strain_pct:+.2f}%, "
              f"volume: {scaled_volume:.4f} Ang^3)")

    print("\n" + "-" * 60)
    print(f"[INFO] Run SIESTA in each '{args.output_dir}/vol_*' folder with a FIXED cell "
          "(positions may relax, but do not let the cell relax -- that would defeat the "
          "purpose of the scan), then aggregate the results with stb-eosAnalysis "
          f"--dir {args.output_dir}.")
    print(color_text("Complete job!\n", 'bold'))


if __name__ == "__main__":
    main()
