#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "2.0.0"

import sys
import os
import time
import argparse
from datetime import datetime

import numpy as np
from pymatgen.io.ase import AseAtomsAdaptor
from stb.core import structure_io
from stb.core import kspace
from stb.core import mace_relax
from stb.core import citations
from stb.core import structure_checks
from stb.core import symmetry as core_symmetry
from stb.core.ase_view import view_structure_interactive
from stb.core.cli import COLORS, color_text, show_intro, print_dual, print_section, print_table
from stb.core.deps import require_mace

REPORT_FILE = "stb_supercell_report.txt"
BIB_FILE = "references.bib"

# Same default vacuum-gap threshold as stb-fetch/stb-kgrid/stb-mlrelax
# (core/kspace.py's other callers) -- used to detect vacuum-padded axes both
# for --ml-relax-cell's cell mask and the symmetry layer-group check.
VACUUM_GAP_ANG = 10.0


def _fail(message, f_out):
    """Prints a red [ERROR] line, closes the report file if one is open, and
    exits with status 1 -- same single error-exit pattern as
    stacking2D.py/dftu.py/fetch.py's own _fail()."""
    print_dual(color_text(f"[ERROR] {message}", 'red'), f_out)
    if f_out:
        f_out.close()
    sys.exit(1)


def parse_matrix(values: list[float]) -> np.ndarray:
    """Builds a 3x3 supercell matrix from either 3 numbers (diagonal) or 9
    numbers (full matrix, row-major) -- same convention as phonopy's DIM tag.
    """
    if len(values) == 3:
        return np.diag(values).astype(float)
    if len(values) == 9:
        return np.array(values, dtype=float).reshape(3, 3)
    raise ValueError(f"Expected 3 (diagonal) or 9 (full matrix) numbers, got {len(values)}.")


def _describe_structure(pmg_structure, vacuum_axes, f_out):
    """Formula/atoms/dimensionality/cell parameters -- same field set as
    mlrelax.py's own [1] INPUT STRUCTURE section."""
    print_dual(f"Formula        : {pmg_structure.composition.reduced_formula}", f_out)
    print_dual(f"Atoms          : {len(pmg_structure)}", f_out)
    print_dual(f"Dimensionality : {kspace.dimensionality_label(vacuum_axes)}", f_out)
    a, b, c, alpha, beta, gamma = pmg_structure.lattice.parameters
    print_dual(f"Cell a,b,c     : {a:.4f}, {b:.4f}, {c:.4f} Ang", f_out)
    print_dual(f"Cell angles    : {alpha:.2f}, {beta:.2f}, {gamma:.2f} deg", f_out)
    print_dual(f"Cell volume    : {pmg_structure.lattice.volume:.4f} Ang^3", f_out)
    if sum(vacuum_axes) == 0:
        print_dual(f"Density        : {len(pmg_structure) / pmg_structure.lattice.volume:.4f} atoms/Ang^3", f_out)


def _validate_structure(pmg_structure, vacuum_axes, f_out):
    """Shared malformation checklist (core.structure_checks) plus a
    space-group label -- same shape as mlrelax.py's _structure_validation(),
    wrapped in try/except by the caller (a validation failure is reported,
    never fatal)."""
    structure_checks.run_malformation_checks(pmg_structure, vacuum_axes, f_out)
    sg_label = core_symmetry.space_group_label(pmg_structure)
    print_dual(f"Space group    : {sg_label}", f_out)
    return sg_label


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Builds a supercell from a SIESTA FDF structure file.", 'bold')}
Give -d 3 numbers for a diagonal supercell, or 9 for a full row-major matrix.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s -f structure.fdf -d 2 2 2\n"
               "  %(prog)s -f structure.fdf -d 2 0 0 0 2 0 0 0 1 -o big_cell.fdf\n"
               "  %(prog)s -f structure.fdf -d 2 2 2 --ml-relax --save-report --view\n"
    )

    parser.add_argument(
        "-f", "--file", dest="filename", type=str, required=True,
        help="Path to the input structure file (.fdf)."
    )
    parser.add_argument(
        "-d", "--dim", type=float, nargs='+', required=True,
        help="Supercell dimensions: 3 numbers for a diagonal supercell (e.g. 2 2 2),\n"
             "or 9 numbers for a full row-major 3x3 matrix (e.g. 2 0 0 0 2 0 0 0 1)."
    )
    parser.add_argument(
        "-o", "--output", type=str, default="supercell.fdf",
        help="Output .fdf file name (default: supercell.fdf)."
    )
    parser.add_argument("-sp", "--symprec", type=float, default=0.01,
                        help="Symmetry tolerance in Angstroms for the before/after "
                             "symmetry analysis (default: 0.01).")

    parser.add_argument("--ml-relax", action="store_true",
                        help="Pre-relax the built supercell with a MACE potential "
                             "(needs the optional 'ml' extra: pip install stb_suite[ml]) "
                             "before writing it out -- positions only by default. "
                             "Off by default.")
    parser.add_argument("--ml-relax-cell", action="store_true",
                        help="With --ml-relax, also relax the cell -- any vacuum-padded "
                             "axis always stays exactly fixed. Only valid together with "
                             "--ml-relax.")
    parser.add_argument("--model", choices=["small", "medium", "large"], default="small",
                        help="MACE-MP-0 foundation model size for --ml-relax (default: small).")
    parser.add_argument("--custom-model", default=None, metavar="PATH",
                        help="Path to a custom fine-tuned .model file for --ml-relax, "
                             "instead of a MACE-MP-0 foundation size.")
    parser.add_argument("--device", choices=["cpu", "cuda"], default="cpu",
                        help="Device to run the MACE model on (default: cpu).")

    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the full run report (including the symmetry "
                             f"analysis) to {REPORT_FILE}. Off by default.")
    parser.add_argument("--view", action="store_true",
                        help="Open an interactive 3D view (via ASE) comparing the input "
                             "structure and the final supercell (page through frames in "
                             "ase-gui) after writing the output file. Needs a display. "
                             "Off by default.")

    parser.add_argument("-v", "--version", action="version", version=f"stb-supercell {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.ml_relax_cell and not args.ml_relax:
        parser.error("--ml-relax-cell is only valid together with --ml-relax.")
    if (args.custom_model or args.model != "small") and not args.ml_relax:
        parser.error("--model/--custom-model are only valid together with --ml-relax.")
    if args.ml_relax:
        require_mace()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    report_path = REPORT_FILE if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    print_dual(color_text("===== STB-SUPERCELL REPORT =====", 'magenta'), f_out)

    model_desc = f"a custom model ({args.custom_model})" if args.custom_model else f"MACE-MP-0 ({args.model})"

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time      : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Input file     : {args.filename}", f_out)
    print_dual(f"Dimensions     : {args.dim}", f_out)
    print_dual(f"Output file    : {args.output}", f_out)
    print_dual(f"ML pre-relax   : {model_desc if args.ml_relax else 'no'}", f_out)
    if args.ml_relax:
        print_dual(f"Relax cell     : {'yes (vacuum axes fixed)' if args.ml_relax_cell else 'no (positions only)'}", f_out)

    if not os.path.exists(args.filename):
        _fail(f"File '{args.filename}' not found.", f_out)

    try:
        matrix = parse_matrix(args.dim)
    except ValueError as e:
        _fail(str(e), f_out)

    rounded = np.round(matrix)
    if not np.allclose(matrix, rounded, atol=1e-6):
        _fail(f"supercell matrix must have integer entries, got:\n{matrix}", f_out)
    matrix = rounded.astype(int)

    det = int(round(np.linalg.det(matrix)))
    if det == 0:
        _fail("supercell matrix is singular (determinant is zero).", f_out)
    if det < 0:
        print_dual(color_text(
            f"[WARNING] supercell matrix has a negative determinant ({det}) -- "
            "the resulting cell is mirrored (still geometrically valid).", 'yellow'), f_out)

    try:
        structure = structure_io.read_fdf(args.filename)
    except ValueError as e:
        _fail(str(e), f_out)

    pmg_structure = structure_io.to_pymatgen(structure)

    print_section("[1] INPUT STRUCTURE", f_out)
    frac_coords_before = [site.frac_coords for site in pmg_structure]
    vacuum_axes_before = kspace.detect_vacuum_axes(frac_coords_before, pmg_structure.lattice.matrix, VACUUM_GAP_ANG)
    _describe_structure(pmg_structure, vacuum_axes_before, f_out)

    print_section("[2] STRUCTURE VALIDATION (pre-transform)", f_out)
    try:
        _validate_structure(pmg_structure, vacuum_axes_before, f_out)
    except Exception as e:
        print_dual(color_text(f"[WARNING] Structure validation could not complete: {e}", 'yellow'), f_out)

    print_section("[3] SUPERCELL CONSTRUCTION", f_out)
    print_dual("Supercell matrix:", f_out)
    for row in matrix:
        print_dual(f"  {row[0]:3d} {row[1]:3d} {row[2]:3d}", f_out)
    print_dual(f"Determinant (cell multiplication factor): {abs(det)}", f_out)

    supercell = pmg_structure.copy()
    supercell.make_supercell(matrix)

    print_dual(f"Output formula : {supercell.composition.reduced_formula}", f_out)
    print_dual(f"Output atoms   : {len(supercell)}", f_out)
    print_dual(f"Output volume  : {supercell.volume:.4f} Ang^3", f_out)

    ml_relax_info = None
    if args.ml_relax:
        print_section("[4] ML PRE-RELAXATION (MACE)", f_out)
        print_dual(f"Model           : {model_desc} (device={args.device})", f_out)
        print_dual(f"Cell relaxation : "
                   f"{'in cell (vacuum axes fixed)' if args.ml_relax_cell else 'positions only'}", f_out)

        frac_coords_super = [site.frac_coords for site in supercell]
        vacuum_axes_super = kspace.detect_vacuum_axes(frac_coords_super, supercell.lattice.matrix, VACUUM_GAP_ANG)

        model_arg = args.custom_model if args.custom_model else args.model
        try:
            calc = mace_relax.get_calculator(model_arg, device=args.device)
        except ValueError as e:
            _fail(str(e), f_out)
        for line in mace_relax.describe_model(model_arg, calc):
            print_dual(line, f_out)

        atoms = AseAtomsAdaptor.get_atoms(supercell)
        atoms.calc = calc
        e0 = atoms.get_potential_energy()
        f0 = float(np.abs(atoms.get_forces()).max())
        a0, b0, c0, _, _, _ = atoms.cell.cellpar()
        vol0 = atoms.get_volume()

        cell_mask = mace_relax.build_cell_mask(vacuum_axes_super) if args.ml_relax_cell else None
        t0 = time.time()
        converged, steps_used = mace_relax.relax(atoms, calc, cell_mask=cell_mask, fmax=0.05, max_steps=200)
        wall_time = time.time() - t0

        e1 = atoms.get_potential_energy()
        f1 = float(np.abs(atoms.get_forces()).max())
        a1, b1, c1, _, _, _ = atoms.cell.cellpar()

        print_dual(f"Steps used : {steps_used} "
                   f"({'converged' if converged else 'hit step cap, NOT converged'})", f_out)
        print_dual(f"Wall time  : {wall_time:.1f} s", f_out)

        n_atoms = len(atoms)
        rows = [
            (["Energy (eV)", f"{e0:.6f}", f"{e1:.6f}",
              f"{e1 - e0:+.6f} ({(e1 - e0) / n_atoms:+.6f}/atom)"], None),
            (["Max force (eV/Ang)", f"{f0:.4f}", f"{f1:.4f}", f"{f1 - f0:+.4f}"], None),
        ]
        if args.ml_relax_cell:
            rows.append((["Lattice a,b,c (Ang)", f"{a0:.4f}, {b0:.4f}, {c0:.4f}",
                          f"{a1:.4f}, {b1:.4f}, {c1:.4f}",
                          f"max {100 * max(abs(a1 - a0) / a0, abs(b1 - b0) / b0, abs(c1 - c0) / c0):+.2f}%"], None))
            vol1 = atoms.get_volume()
            rows.append((["Volume (Ang^3)", f"{vol0:.4f}", f"{vol1:.4f}",
                          f"{100 * (vol1 - vol0) / vol0:+.2f}%"], None))
        print_table(["Quantity", "Before", "After", "Change"], rows, f_out)

        supercell = AseAtomsAdaptor.get_structure(atoms)
        ml_relax_info = (converged, steps_used, e0, e1)

    final_atoms = AseAtomsAdaptor.get_atoms(supercell)

    print_section("[5] STRUCTURE VALIDATION (post-transform)", f_out)
    frac_coords_final = [site.frac_coords for site in supercell]
    vacuum_axes_final = kspace.detect_vacuum_axes(frac_coords_final, supercell.lattice.matrix, VACUUM_GAP_ANG)
    try:
        _validate_structure(supercell, vacuum_axes_final, f_out)
    except Exception as e:
        print_dual(color_text(f"[WARNING] Structure validation could not complete: {e}", 'yellow'), f_out)

    print_section("[6] SYMMETRY ANALYSIS (BEFORE / AFTER)", f_out)
    print_dual(f"Detailed symmetry analysis (Tolerance: {args.symprec} Ang):", f_out)
    before_info = core_symmetry.symmetry_summary(pmg_structure, args.symprec, VACUUM_GAP_ANG)
    after_info = core_symmetry.symmetry_summary(supercell, args.symprec, VACUUM_GAP_ANG)
    if "Error" in before_info or "Error" in after_info:
        print_dual(color_text("[WARNING] Symmetry analysis failed for at least one structure.", 'yellow'), f_out)
        print_dual(f"  Before: {before_info.get('Error', 'OK')}", f_out)
        print_dual(f"  After : {after_info.get('Error', 'OK')}", f_out)
    else:
        properties = ["Crystal System", "Space Group", "Layer Group", "Point Group", "Hall Symbol"]
        rows = [([prop, str(before_info.get(prop, "N/A")), str(after_info.get(prop, "N/A"))], None)
                for prop in properties]
        print_table(["Property", "Before", "After"], rows, f_out)

    print_section("[7] WRITING OUTPUT FILE", f_out)
    header_comment = [
        "Supercell built by stb-supercell from an input structure.",
        f"Input file: {args.filename}",
        f"Transformation matrix: {matrix.tolist()} (determinant {det}, "
        f"{abs(det)}x multiplication, {len(pmg_structure)} -> {len(supercell)} atoms).",
    ]
    if ml_relax_info is not None:
        converged, steps_used, e0, e1 = ml_relax_info
        header_comment.append(
            f"ML pre-relaxed with {model_desc} "
            f"({'converged' if converged else 'NOT converged'} in {steps_used} step(s), "
            f"E = {e1:.6f} eV, delta E = {e1 - e0:+.6f} eV)."
        )
    new_structure = structure_io.from_pymatgen(supercell, species_meta=structure.species_meta)
    structure_io.write_fdf(new_structure, args.output, header_comment=header_comment)
    print_dual(color_text(f"[OK] Structure written to '{args.output}'.", 'green'), f_out)

    print_section("[8] REFERENCES", f_out)
    bib_entries = [citations.SIESTA, citations.SIESTA_RECENT]
    if args.ml_relax:
        bib_entries.append(citations.MACE)
        if not args.custom_model:
            bib_entries.append(citations.MACE_MP)
    citations.write_bib_file(BIB_FILE, bib_entries)
    print_dual(color_text(
        f"[OK] Citations for the methods used in this run written to '{BIB_FILE}' "
        f"({len(bib_entries)} entries).", 'green'), f_out)

    print_section("[9] SUMMARY & FILES", f_out)
    print_dual("Status         : OK", f_out)
    print_dual(f"Input file     : {args.filename}", f_out)
    print_dual(f"Output file    : {args.output}", f_out)
    print_dual(f"References     : {BIB_FILE}", f_out)
    if report_path:
        print_dual(f"Report         : {report_path}", f_out)

    if f_out:
        f_out.close()

    # --view runs last, after every check/report section above has already
    # printed, so a blocking GUI window never delays or hides them -- shows
    # both frames (input structure vs. final supercell) so the user can page
    # through the actual comparison in ase-gui.
    if args.view:
        input_atoms = AseAtomsAdaptor.get_atoms(pmg_structure)
        view_structure_interactive([input_atoms, final_atoms])


if __name__ == "__main__":
    main()
