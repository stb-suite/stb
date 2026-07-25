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
from pymatgen.core.periodic_table import Element
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
from stb.core.passivation import passivate_dangling_bonds

REPORT_FILE = "stb_passivate_report.txt"
BIB_FILE = "references.bib"

# Same default vacuum-gap threshold as stb-unitcell/stb-supercell/stb-fetch/
# stb-kgrid/stb-mlrelax (core/kspace.py's other callers) -- passivate's usual
# input is a cut slab (vacuum-padded), so unlike stb-crystalbuilder this tool
# must stay vacuum-aware for both --ml-relax-cell's cell mask and the
# before/after Layer Group check.
VACUUM_GAP_ANG = 10.0


def _fail(message, f_out):
    """Prints a red [ERROR] line, closes the report file if one is open, and
    exits with status 1 -- same single error-exit pattern as
    stb-unitcell/stb-supercell/stb-slab/stb-crystalbuilder's own _fail()."""
    print_dual(color_text(f"[ERROR] {message}", 'red'), f_out)
    if f_out:
        f_out.close()
    sys.exit(1)


def _describe_structure(pmg_structure, vacuum_axes, f_out):
    """Formula/atoms/dimensionality/cell parameters -- same field set as
    stb-supercell/stb-unitcell's own [1] INPUT STRUCTURE section."""
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
    space-group label -- same shape as stb-supercell/stb-slab's own
    _validate_structure(), wrapped in try/except by the caller (a validation
    failure is reported, never fatal)."""
    structure_checks.run_malformation_checks(pmg_structure, vacuum_axes, f_out)
    sg_label = core_symmetry.space_group_label(pmg_structure)
    print_dual(f"Space group    : {sg_label}", f_out)
    return sg_label


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Caps dangling bonds on a cut surface with a passivating atom.", 'bold')}
Detects undercoordinated atoms (e.g. from stb-slab) via local coordination
geometry and places a passivant (default H) along each missing-bond
direction. Only single-missing-bond sites are auto-passivated -- sites
missing 2+ bonds are geometrically underdetermined from local coordination
alone and are reported instead of guessed.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s -f slab.fdf\n"
               "  %(prog)s -f slab.fdf --passivant F --bond-length 1.4 -o slab_F.fdf\n"
               "  %(prog)s -f slab.fdf --ml-relax --save-report --view\n"
    )

    parser.add_argument("-f", "--file", dest="filename", type=str, required=True,
                        help="Path to the input structure file (.fdf), typically a cut slab.")
    parser.add_argument("--passivant", type=str, default="H",
                        help="Element to cap dangling bonds with (default: H).")
    parser.add_argument("--cutoff", type=float, default=None,
                        help="Neighbor-search radius in Angstrom, used to determine each atom's "
                             "coordination number. Default: auto-detected as 1.25x the shortest "
                             "pairwise distance in the structure.")
    parser.add_argument("--bond-length", type=float, default=None,
                        help="Passivant bond length in Angstrom. Default: auto, per surface-"
                             "species/passivant pair, as the sum of their pymatgen atomic radii "
                             "(an approximation -- e.g. Si+H = 1.35 Ang vs. the experimental "
                             "Si-H ~1.48 Ang).")
    parser.add_argument("-sp", "--symprec", type=float, default=0.01,
                        help="Symmetry tolerance in Angstroms for the before/after "
                             "symmetry analysis (default: 0.01).")

    parser.add_argument("--ml-relax", action="store_true",
                        help="Pre-relax the passivated structure with a MACE potential "
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

    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the full run report (including the symmetry "
                             f"analysis) to {REPORT_FILE}. Off by default.")
    parser.add_argument("--view", action="store_true",
                        help="Open an interactive 3D view (via ASE) comparing the input "
                             "structure and the passivated structure (page through frames "
                             "in ase-gui) after writing the output file. Needs a display. "
                             "Off by default.")

    parser.add_argument("-o", "--output", type=str, default="passivated.fdf",
                        help="Output .fdf file name (default: passivated.fdf).")
    parser.add_argument("-v", "--version", action="version", version=f"stb-passivate {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    try:
        Element(args.passivant)
    except ValueError as e:
        parser.error(str(e))

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

    print_dual(color_text("===== STB-PASSIVATE REPORT =====", 'magenta'), f_out)

    model_desc = f"a custom model ({args.custom_model})" if args.custom_model else f"MACE-MP-0 ({args.model})"

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time      : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Input file     : {args.filename}", f_out)
    print_dual(f"Passivant      : {args.passivant}", f_out)
    print_dual(f"Cutoff         : {args.cutoff if args.cutoff is not None else 'auto'}", f_out)
    print_dual(f"Bond length    : {args.bond_length if args.bond_length is not None else 'auto'}", f_out)
    print_dual(f"Output file    : {args.output}", f_out)
    print_dual(f"ML pre-relax   : {model_desc if args.ml_relax else 'no'}", f_out)
    if args.ml_relax:
        print_dual(f"Relax cell     : {'yes (vacuum axes fixed)' if args.ml_relax_cell else 'no (positions only)'}", f_out)

    if not os.path.exists(args.filename):
        _fail(f"File '{args.filename}' not found.", f_out)

    try:
        structure = structure_io.read_fdf(args.filename)
    except (FileNotFoundError, ValueError) as e:
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

    print_section("[3] PASSIVATION", f_out)
    new_pmg, report = passivate_dangling_bonds(
        pmg_structure, passivant=args.passivant, cutoff=args.cutoff, bond_length=args.bond_length)

    n_passivated = len(report["passivated"])
    n_unresolved = len(report["unresolved"])
    print_dual(f"Dangling bonds found : {n_passivated + n_unresolved}", f_out)
    print_dual(f"Auto-passivated      : {n_passivated} with {args.passivant}", f_out)
    if n_unresolved:
        print_dual(color_text(
            f"[WARNING] {n_unresolved} atom(s) missing 2+ bonds -- left unpassivated "
            "(geometrically underdetermined from local coordination alone):", 'yellow'), f_out)
        for idx, symbol, pos, deficit in report["unresolved"]:
            print_dual(f"    #{idx + 1:<4} {symbol:<3} deficit={deficit}  at {pos}", f_out)

    print_dual(f"Output formula : {new_pmg.composition.reduced_formula}", f_out)
    print_dual(f"Output atoms   : {len(new_pmg)}", f_out)

    species_meta = structure_io.ensure_species_id(dict(structure.species_meta), args.passivant)

    ml_relax_info = None
    if args.ml_relax:
        print_section("[4] ML PRE-RELAXATION (MACE)", f_out)
        print_dual(f"Model           : {model_desc}", f_out)
        print_dual(f"Cell relaxation : "
                   f"{'in cell (vacuum axes fixed)' if args.ml_relax_cell else 'positions only'}", f_out)

        frac_coords_new = [site.frac_coords for site in new_pmg]
        vacuum_axes_new = kspace.detect_vacuum_axes(frac_coords_new, new_pmg.lattice.matrix, VACUUM_GAP_ANG)

        model_arg = args.custom_model if args.custom_model else args.model
        calc = mace_relax.get_calculator(model_arg)
        for line in mace_relax.describe_model(model_arg, calc):
            print_dual(line, f_out)

        atoms = AseAtomsAdaptor.get_atoms(new_pmg)
        atoms.calc = calc
        e0 = atoms.get_potential_energy()
        f0 = float(np.abs(atoms.get_forces()).max())
        a0, b0, c0, _, _, _ = atoms.cell.cellpar()
        vol0 = atoms.get_volume()

        cell_mask = mace_relax.build_cell_mask(vacuum_axes_new) if args.ml_relax_cell else None
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

        new_pmg = AseAtomsAdaptor.get_structure(atoms)
        ml_relax_info = (converged, steps_used, e0, e1)

    final_atoms = AseAtomsAdaptor.get_atoms(new_pmg)

    print_section("[5] STRUCTURE VALIDATION (post-transform)", f_out)
    frac_coords_final = [site.frac_coords for site in new_pmg]
    vacuum_axes_final = kspace.detect_vacuum_axes(frac_coords_final, new_pmg.lattice.matrix, VACUUM_GAP_ANG)
    try:
        _validate_structure(new_pmg, vacuum_axes_final, f_out)
    except Exception as e:
        print_dual(color_text(f"[WARNING] Structure validation could not complete: {e}", 'yellow'), f_out)

    print_section("[6] SYMMETRY ANALYSIS (BEFORE / AFTER)", f_out)
    print_dual(f"Detailed symmetry analysis (Tolerance: {args.symprec} Ang):", f_out)
    before_info = core_symmetry.symmetry_summary(pmg_structure, args.symprec, VACUUM_GAP_ANG)
    after_info = core_symmetry.symmetry_summary(new_pmg, args.symprec, VACUUM_GAP_ANG)
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
        "Dangling bonds passivated by stb-passivate from an input structure.",
        f"Input file: {args.filename}",
        f"Passivant: {args.passivant}. Dangling bonds found: {n_passivated + n_unresolved}, "
        f"auto-passivated: {n_passivated}, unresolved (2+ missing bonds): {n_unresolved}.",
    ]
    if ml_relax_info is not None:
        converged, steps_used, e0, e1 = ml_relax_info
        header_comment.append(
            f"ML pre-relaxed with {model_desc} "
            f"({'converged' if converged else 'NOT converged'} in {steps_used} step(s), "
            f"E = {e1:.6f} eV, delta E = {e1 - e0:+.6f} eV)."
        )
    new_structure = structure_io.from_pymatgen(new_pmg, species_meta=species_meta)
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
    # both frames (input structure vs. passivated structure) so the user can
    # page through the actual comparison in ase-gui.
    if args.view:
        input_atoms = AseAtomsAdaptor.get_atoms(pmg_structure)
        view_structure_interactive([input_atoms, final_atoms])


if __name__ == "__main__":
    main()
