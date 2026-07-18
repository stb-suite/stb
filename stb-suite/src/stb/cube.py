#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.9.1"

import os
import sys
import argparse
from datetime import datetime
import numpy as np

from stb.core.deps import require_sisl
sisl = require_sisl()

from stb.core.cli import COLORS, color_text, show_intro, print_dual, print_section

REPORT_FILE = "stb_cube_report.txt"

# sisl's own read_grid(index=...) string/int convention for selecting a spin
# component -- see rhoSileSiesta.read_grid (sisl.io.siesta.binaries). Only
# meaningful for a simple COLLINEAR polarized calculation (2 components);
# non-collinear/spin-orbit calculations have a genuinely different
# 4-component spin texture (charge, mx, my, mz) that this up/down/diff split
# doesn't apply to -- see detect_spin_configuration below, which is what
# tells the two cases apart.
SPIN_INDEX = {'total': 'total', 'up': 0, 'down': 1, 'diff': 'z'}


def get_cell_matrix(lattice_obj):
    """
    Robustly retrieves the unit cell matrix from sisl object.
    Supports: .cell, .matrix, .vectors, .get_cell()
    """
    attributes = ['cell', 'matrix', 'vectors']

    for attr in attributes:
        if hasattr(lattice_obj, attr):
            return getattr(lattice_obj, attr)

    if hasattr(lattice_obj, 'get_cell') and callable(lattice_obj.get_cell):
        return lattice_obj.get_cell()

    try:
        return np.array(lattice_obj)
    except Exception:
        return None


def detect_spin_configuration(label, f_out=None):
    """Reads <label>.HSX (the SIESTA Hamiltonian file), if present, via sisl
    to determine the actual spin configuration: unpolarized, polarized
    (collinear), non-collinear, or spin-orbit. This is the authoritative
    source for "what kind of spin calculation was this" -- unlike a grid
    file's raw component count alone, which can't tell a 4-component
    non-collinear/spin-orbit run apart from anything else that happens to
    have more than one component.

    Returns a sisl.physics.Spin object, or None if no .HSX is found or it
    can't be read (not fatal -- --spin still works from the grid's own
    component count alone, just without this extra identification/
    cross-check). Verified against a real SIESTA spin-polarized O2 run,
    where this correctly reports Spin{polarized}.
    """
    hsx_file = f"{label}.HSX"
    if not os.path.exists(hsx_file):
        return None
    try:
        H = sisl.get_sile(hsx_file).read_hamiltonian()
        return H.spin
    except Exception as e:
        print_dual(color_text(
            f"[WARNING] Found '{hsx_file}' but could not read its spin "
            f"configuration: {e}", 'yellow'), f_out)
        return None


def main():
    parser = argparse.ArgumentParser(
        description="Convert SIESTA Grid files (VT, VH, RHO, BADER) to Gaussian Cube format.",
        epilog="Example usage:\n"
               "  stb_convert --label graphene --type RHO\n"
               "  stb_convert -l system -t VT -o potential.cube",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("-l", "--label", help="System Label (required).", required=True)
    parser.add_argument("-t", "--type", help="File type to convert.", choices=['VT', 'VH', 'RHO', 'BADER'], required=True)
    parser.add_argument("-o", "--output", help="Custom output filename. Default: label_type.cube", default=None)
    parser.add_argument(
        "--spin", choices=['total', 'up', 'down', 'diff'], default='total',
        help="Which spin channel to convert, for spin-polarized calculations: "
             "'total' (up+down, default), 'up', 'down', or 'diff' (up-down, i.e. "
             "the local magnetization/spin density). Auto-identified and "
             "cross-checked against <label>.HSX when present: ignored (falls back "
             "to 'total') for non-spin-polarized or non-collinear/spin-orbit "
             "calculations, where up/down/diff doesn't apply. Verified against a "
             "real SIESTA spin-polarized O2 run: 'diff' integrates to exactly the "
             "spin moment SIESTA itself reports."
    )
    parser.add_argument(
        "--save-report", action="store_true",
        help=f"Also persist the report to {REPORT_FILE}. Off by default."
    )

    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-cube {VERSION}")

    args = parser.parse_args()

    # Filenames
    input_file = f"{args.label}.{args.type}"
    xv_file = f"{args.label}.XV"
    default_output = (f"{args.label}_{args.type}.cube" if args.spin == 'total'
                      else f"{args.label}_{args.type}_{args.spin}.cube")
    output_file = args.output if args.output else default_output

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - Grid Converter",
            "Professional Grid to Cube Converter",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("GRID CONVERTER:", 'bold'))
    print("-" * 60)

    report_path = REPORT_FILE if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    def fail(message):
        print_dual(color_text(f"[FAIL] {message}", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_dual(color_text("===== STB-CUBE REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Input file        : {input_file}", f_out)
    print_dual(f"Output file       : {output_file}", f_out)
    print_dual(f"Spin channel      : {args.spin} (requested)", f_out)
    if report_path:
        print_dual(f"Report file       : {report_path}", f_out)

    print_section("[1] READING GRID", f_out)

    # --- Check Input ---
    if not os.path.exists(input_file):
        fail(f"Input file '{input_file}' not found.")

    # --- Load Grid (and identify/cross-check spin configuration) ---
    try:
        sile = sisl.get_sile(input_file)
        nspin, _mesh = sile.read_grid_size()
    except Exception as e:
        fail(f"Could not read grid: {e}")

    spin_config = detect_spin_configuration(args.label, f_out)
    spin_kind = None
    if spin_config is not None:
        spin_kind = ("spin-orbit" if spin_config.is_spinorbit else
                     "non-collinear" if spin_config.is_noncolinear else
                     "polarized" if spin_config.is_polarized else
                     "unpolarized")
        print_dual(f"[INFO] Spin configuration (from {args.label}.HSX): {spin_kind}", f_out)
        expected_nspin = {"unpolarized": 1, "polarized": 2,
                          "non-collinear": 4, "spin-orbit": 4}[spin_kind]
        if expected_nspin != nspin:
            print_dual(color_text(
                f"[WARNING] '{args.label}.HSX' reports {spin_kind} ({expected_nspin} "
                f"expected component(s)) but '{input_file}' has {nspin} -- they may not "
                f"belong to the same calculation.", 'yellow'), f_out)

    effective_spin = args.spin
    if nspin == 1:
        if args.spin != 'total':
            print_dual(color_text(
                f"[WARNING] '{input_file}' is not spin-polarized (1 component) -- "
                f"ignoring --spin {args.spin}.", 'yellow'), f_out)
            effective_spin = 'total'
    elif spin_kind in ("non-collinear", "spin-orbit"):
        if args.spin != 'total':
            print_dual(color_text(
                f"[WARNING] --spin {args.spin} assumes a simple collinear polarized "
                f"calculation; '{args.label}.HSX' is {spin_kind} (a 4-component "
                f"charge/mx/my/mz spin texture, not a scalar up/down split) -- falling "
                f"back to 'total'.", 'yellow'), f_out)
            effective_spin = 'total'
    elif nspin > 1 and effective_spin != 'total':
        print_dual(f"[INFO] Spin-polarized grid detected ({nspin} components) -- "
                   f"converting the '{effective_spin}' channel.", f_out)

    try:
        grid = sile.read_grid(index=SPIN_INDEX[effective_spin])
    except Exception as e:
        fail(f"Could not read grid: {e}")

    print_dual(f"[OK] Read the file {input_file}", f_out)

    # --- Load Geometry (Optional) ---
    if os.path.exists(xv_file):
        print_dual(f"[INFO] Found geometry file: {xv_file}. Mapping atoms...", f_out)
        try:
            geom = sisl.get_sile(xv_file).read_geometry()
            grid_cell = get_cell_matrix(grid.lattice)
            # A stale/mismatched .XV (e.g. from a different relaxation step)
            # would silently place atoms at the wrong positions relative to
            # this grid -- same cross-check stb-bader already does before
            # trusting a .RHO + geometry pairing.
            if grid_cell is not None and not np.allclose(grid_cell, geom.cell, atol=1e-3):
                print_dual(color_text(
                    f"[WARNING] Lattice mismatch between the {args.type} grid and "
                    f"'{xv_file}' -- they may not belong to the same calculation. Atom "
                    f"positions in the .cube output may be wrong.", 'yellow'), f_out)
            grid.geometry = geom
        except Exception as e:
            print_dual(color_text(f"[WARNING] Failed to read geometry from .XV: {e}", 'yellow'), f_out)
    else:
        print_dual(color_text(
            f"[WARNING] No .XV file found. Atoms will be missing in the .cube output.", 'yellow'), f_out)

    # --- Metadata Display ---
    print_dual(f"Shape             : {grid.shape}", f_out)

    cell = get_cell_matrix(grid.lattice)
    if cell is not None:
        if np.allclose(cell, np.diag(np.diag(cell)), atol=1e-6):
            print_dual(f"Lattice           : [{cell[0,0]:.2f}, {cell[1,1]:.2f}, {cell[2,2]:.2f}] Ang", f_out)
        else:
            # Non-orthogonal cell (common for this suite's 2D/hexagonal
            # fixtures) -- printing only the diagonal would silently hide
            # real off-diagonal/shear terms.
            print_dual("Lattice (Ang, non-orthogonal):", f_out)
            for row in cell:
                print_dual(f"    [{row[0]:9.4f}, {row[1]:9.4f}, {row[2]:9.4f}]", f_out)

    if args.type.upper() in ("RHO", "BADER"):
        # Approximate integral check -- uses grid.grid/grid.volume as sisl
        # left them after read_grid() (both in Angstrom-based units), so this
        # must run BEFORE the unit-restoration step right below.
        total_val = grid.grid.sum() * grid.volume / (grid.shape[0] * grid.shape[1] * grid.shape[2])
        if effective_spin == 'diff':
            # Verified against a real SIESTA spin-polarized O2 run: this
            # integral matches SIESTA's own reported spin moment |S| exactly.
            print_dual(f"Integral          : {total_val:.4f} e (net spin moment |S| -- "
                       f"compare to SIESTA's own 'spin moment' report)", f_out)
        elif effective_spin != 'total':
            print_dual(f"Integral          : {total_val:.4f} e ({effective_spin}-spin channel only)", f_out)
        else:
            print_dual(f"Integral          : {total_val:.4f} e (Check if close to total electrons)", f_out)

    # --- Restore native SIESTA units for density-type grids before writing ---
    # sisl's read_grid() scales RHO/BADER values by `sile.grid_unit` (the
    # Bohr^-3 -> Ang^-3 conversion, so the Integral check above lines up with
    # grid.volume, itself in Ang^3) -- but its cube writer (grid.write,
    # default unit="Bohr") only converts the LATTICE/atom coordinates to
    # Bohr, never the density values (sisl's own write_grid docstring: "grid
    # data is assumed to be unit-less"). Left uncorrected, every RHO/BADER
    # .cube this tool writes would pair Bohr-unit coordinates with
    # e/Ang^3-valued density -- ~6.75x too large relative to the e/Bohr^3
    # that standard cube-reading software (VESTA, VMD, Multiwfn, ...) assumes
    # goes with Bohr coordinates. Verified directly: sisl's grid_unit for
    # RHO/BADER (6.748334...) is exactly 1/(0.529177 Ang)^3, the Bohr^3 ->
    # Ang^3 conversion; dividing back by it restores e/Bohr^3, matching the
    # true SIESTA-native density paired correctly with Bohr coordinates. VT/VH
    # don't need this -- their grid_unit (13.6057) is a Ry -> eV energy
    # conversion, unrelated to any length-cubed density mismatch.
    if args.type.upper() in ("RHO", "BADER"):
        grid.grid = grid.grid / sile.grid_unit

    print_section("[2] WRITING OUTPUT", f_out)
    try:
        grid.write(output_file)
        print_dual(f"[OK] Saved {output_file}", f_out)
    except Exception as e:
        fail(f"Failed to write output file: {e}")

    print_section("[3] SUMMARY & FILES", f_out)
    print_dual("Status            : OK", f_out)
    print_dual(f"Output            : {output_file}", f_out)
    if report_path:
        print_dual(f"Report            : {report_path}", f_out)

    if f_out:
        f_out.close()

    print("\n" + "-" * 60)
    # Funny closing phrase
    print(color_text("Conversion complete. It's hip to be square (cube).\n", 'bold'))


if __name__ == "__main__":
    main()
