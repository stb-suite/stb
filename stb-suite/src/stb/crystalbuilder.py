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
from datetime import datetime

import numpy as np
from pymatgen.core import Structure, Lattice
from pymatgen.core.periodic_table import Element
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.groups import SpaceGroup
from pymatgen.io.ase import AseAtomsAdaptor
import argparse
from stb.core import structure_io
from stb.core import kspace
from stb.core import mace_relax
from stb.core import citations
from stb.core import structure_checks
from stb.core import symmetry as core_symmetry
from stb.core.symmetry import reduce_to_unitcell, UNITCELL_MODES
from stb.core.ase_view import view_structure_interactive
from stb.core.cli import COLORS, color_text, show_intro, print_dual, print_section, print_table
from stb.core.deps import require_mace

REPORT_FILE = "stb_crystalbuilder_report.txt"
BIB_FILE = "references.bib"

# Same default vacuum-gap threshold as stb-unitcell/stb-fetch/stb-kgrid/stb-mlrelax
# (core/kspace.py's other callers) -- kept only for the [1] dimensionality line;
# every space group is inherently 3D-periodic, so this is never expected to detect
# vacuum here, unlike stb-unitcell which operates on arbitrary user structures.
VACUUM_GAP_ANG = 10.0


def _fail(message, f_out):
    """Prints a red [ERROR] line, closes the report file if one is open, and
    exits with status 1 -- same single error-exit pattern as stb-unitcell/
    stacking2D.py/supercell.py/slab.py/nanotube.py/defect.py/sqs.py's own
    _fail()."""
    print_dual(color_text(f"[ERROR] {message}", 'red'), f_out)
    if f_out:
        f_out.close()
    sys.exit(1)


def parse_sites(raw_sites):
    """Parses [['Ni','0','0','0'], ['O','0.25','0.25','0.25']] (argparse's
    --site nargs=4, all strings) into [('Ni', [0.0, 0.0, 0.0]), ...].
    """
    sites = []
    for symbol, x, y, z in raw_sites:
        try:
            coords = [float(x), float(y), float(z)]
        except ValueError:
            raise ValueError(f"--site '{symbol} {x} {y} {z}': x/y/z must be numbers.")
        sites.append((symbol, coords))
    return sites


def resolve_spacegroup(spec):
    """Returns spec as an int if it looks like one, else as-is (a symbol string)."""
    try:
        return int(spec)
    except ValueError:
        return spec


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Builds a structure from a space group and Wyckoff positions.", 'bold')}
Give only the symmetrically-distinct sites (e.g. one per Wyckoff letter) --
pymatgen expands the rest via the space group's own symmetry operations.

Bulk (3D) crystals only: the 230 space groups accepted by --spacegroup are
the international crystallographic groups for 3D-periodic crystals -- there
is no 2D/1D analog (layer/rod groups) here. To get a 2D or 1D structure,
build the 3D bulk crystal with this tool first, then cut it down with
stb-slab (a surface slab), stb-nanotube (rolled into a tube/ribbon), or
stb-2Dstacking (an isolated monolayer) -- e.g. build bulk graphite
(space group P6_3/mmc) here, then stb-slab --hkl 0 0 1 it into a graphene
monolayer.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s --spacegroup Fm-3m --a 3.52 --site Ni 0 0 0\n"
               "  %(prog)s --spacegroup 225 --a 3.52 --site Ni 0 0 0 -o fcc_ni.fdf\n"
               "  %(prog)s --spacegroup Fd-3m --a 8.396 \\\n"
               "      --site Fe 0.125 0.125 0.125 --site Fe 0.5 0.5 0.5 \\\n"
               "      --site O 0.2548 0.2548 0.2548\n"
               "  %(prog)s --spacegroup Fm-3m --a 3.52 --site Ni 0 0 0 \\\n"
               "      --reduce primitive --ml-relax --save-report --view\n"
               "\n"
               "Notes:\n"
               "  - Bulk (3D) crystals only -- see the description above for how to get\n"
               "    a 2D/1D structure out of this tool's output.\n"
    )

    parser.add_argument("--spacegroup", type=str, required=True,
                        help="Space group as an international symbol (e.g. 'Fm-3m') or "
                             "number (e.g. '225', must be between 1 and 230).")
    parser.add_argument("--a", type=float, required=True, help="Lattice constant a (Ang).")
    parser.add_argument("--b", type=float, default=None, help="Lattice constant b (Ang). Default: same as --a.")
    parser.add_argument("--c", type=float, default=None, help="Lattice constant c (Ang). Default: same as --a.")
    parser.add_argument("--alpha", type=float, default=90.0, help="Lattice angle alpha, degrees (default: 90).")
    parser.add_argument("--beta", type=float, default=90.0, help="Lattice angle beta, degrees (default: 90).")
    parser.add_argument("--gamma", type=float, default=90.0, help="Lattice angle gamma, degrees (default: 90).")

    parser.add_argument("--site", dest="sites", action="append", nargs=4, required=True,
                        metavar=("SYMBOL", "X", "Y", "Z"),
                        help="One symmetrically-distinct Wyckoff site: element symbol and "
                             "fractional x y z. Repeat --site for each distinct site.")

    parser.add_argument("--symprec", type=float, default=1e-3,
                        help="Symmetry precision, used both for the post-build symmetry "
                             "detection and for --reduce (default: 1e-3, matches "
                             "stb-symmetry/stb-unitcell).")
    parser.add_argument("--angle-tolerance", type=float, default=5.0,
                        help="Symmetry angle tolerance in degrees (default: 5.0, pymatgen's own default).")

    parser.add_argument("--reduce", choices=UNITCELL_MODES, default=None,
                        help="Reduce the space-group-expanded structure to its "
                             "primitive/conventional/refined cell (via stb-unitcell's own "
                             "reduction) before writing. Off by default -- the raw "
                             "space-group-expanded structure is written as-is.")

    parser.add_argument("--ml-relax", action="store_true",
                        help="Pre-relax the final structure with a MACE potential "
                             "(needs the optional 'ml' extra: pip install stb_suite[ml]) "
                             "before writing it out -- positions only by default. "
                             "Off by default.")
    parser.add_argument("--ml-relax-cell", action="store_true",
                        help="With --ml-relax, also relax the cell. Only valid together "
                             "with --ml-relax.")
    parser.add_argument("--model", choices=["small", "medium", "large"], default="small",
                        help="MACE-MP-0 foundation model size for --ml-relax (default: small).")
    parser.add_argument("--custom-model", default=None, metavar="PATH",
                        help="Path to a custom fine-tuned .model file for --ml-relax, "
                             "instead of a MACE-MP-0 foundation size.")

    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the full run report (including the symmetry "
                             f"analysis) to {REPORT_FILE}. Off by default.")
    parser.add_argument("--view", action="store_true",
                        help="Open an interactive 3D view (via ASE) comparing the raw "
                             "space-group-expanded structure and the final structure "
                             "(page through frames in ase-gui) after writing the output "
                             "file. Needs a display. Off by default.")

    parser.add_argument("-o", "--output", type=str, default="crystal.fdf",
                        help="Output .fdf file name (default: crystal.fdf).")
    parser.add_argument("-v", "--version", action="version", version=f"stb-crystalbuilder {VERSION}")
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

    print_dual(color_text("===== STB-CRYSTALBUILDER REPORT =====", 'magenta'), f_out)

    model_desc = f"a custom model ({args.custom_model})" if args.custom_model else f"MACE-MP-0 ({args.model})"

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time      : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Space group    : {args.spacegroup}", f_out)
    b_meta = args.b if args.b is not None else args.a
    c_meta = args.c if args.c is not None else args.a
    print_dual(f"Lattice        : a={args.a} b={b_meta} c={c_meta} "
               f"alpha={args.alpha} beta={args.beta} gamma={args.gamma}", f_out)
    print_dual(f"Symprec        : {args.symprec}", f_out)
    print_dual(f"Angle tolerance: {args.angle_tolerance} deg", f_out)
    print_dual(f"Reduce mode    : {args.reduce if args.reduce else 'no'}", f_out)
    print_dual(f"Output file    : {args.output}", f_out)
    print_dual(f"ML pre-relax   : {model_desc if args.ml_relax else 'no'}", f_out)
    if args.ml_relax:
        print_dual(f"Relax cell     : {'yes' if args.ml_relax_cell else 'no (positions only)'}", f_out)

    print_section("[1] BUILD PARAMETERS", f_out)
    try:
        sites = parse_sites(args.sites)
    except ValueError as e:
        _fail(str(e), f_out)

    for symbol, _ in sites:
        try:
            Element(symbol)
        except ValueError as e:
            _fail(str(e), f_out)

    b = args.b if args.b is not None else args.a
    c = args.c if args.c is not None else args.a
    lattice = Lattice.from_parameters(args.a, b, c, args.alpha, args.beta, args.gamma)

    spacegroup = resolve_spacegroup(args.spacegroup)
    if isinstance(spacegroup, int) and not (1 <= spacegroup <= 230):
        _fail(f"Space group number must be between 1 and 230 (inclusive), got {spacegroup}.", f_out)

    print_dual(f"Requested space group : {args.spacegroup}", f_out)
    print_dual(f"Lattice a,b,c         : {args.a}, {b}, {c} Ang", f_out)
    print_dual(f"Lattice angles        : {args.alpha}, {args.beta}, {args.gamma} deg", f_out)
    rows = [([symbol, f"{coords[0]:.6f}", f"{coords[1]:.6f}", f"{coords[2]:.6f}"], None)
            for symbol, coords in sites]
    print_table(["Site", "x (frac)", "y (frac)", "z (frac)"], rows, f_out)

    species = [symbol for symbol, _ in sites]
    coords = [xyz for _, xyz in sites]

    try:
        structure = Structure.from_spacegroup(spacegroup, lattice, species, coords)
    except ValueError as e:
        _fail(str(e), f_out)

    print_section("[2] STRUCTURE VALIDATION (pre-transform)", f_out)
    frac_coords_before = [site.frac_coords for site in structure]
    vacuum_axes_before = kspace.detect_vacuum_axes(frac_coords_before, structure.lattice.matrix, VACUUM_GAP_ANG)
    print_dual(f"Dimensionality : {kspace.dimensionality_label(vacuum_axes_before)} -- "
               "every space group is inherently 3D-periodic; use stb-slab/stb-nanotube/"
               "stb-2Dstacking on the output to cut a 2D/1D structure from it.", f_out)
    try:
        structure_checks.run_malformation_checks(structure, vacuum_axes_before, f_out)
    except Exception as e:
        print_dual(color_text(f"[WARNING] Structure validation could not complete: {e}", 'yellow'), f_out)

    print_dual(f"Output formula : {structure.composition.reduced_formula}", f_out)
    print_dual(f"Output atoms   : {len(structure)}", f_out)

    print_section("[3] SYMMETRY DETECTION", f_out)
    try:
        sga = SpacegroupAnalyzer(structure, symprec=args.symprec, angle_tolerance=args.angle_tolerance)
        detected_number = sga.get_space_group_number()
        detected_symbol = sga.get_space_group_symbol()
    except Exception as e:
        _fail(f"symmetry detection failed: {e}", f_out)

    print_dual(f"Space group    : {detected_symbol} (No. {detected_number})", f_out)
    print_dual(f"Point group    : {sga.get_point_group_symbol()}", f_out)
    print_dual(f"Crystal system : {sga.get_crystal_system()}", f_out)
    print_dual(f"Hall symbol    : {sga.get_hall()}", f_out)

    try:
        requested_number = SpaceGroup(str(spacegroup)).int_number if isinstance(spacegroup, str) \
            else SpaceGroup.from_int_number(spacegroup).int_number
    except ValueError:
        requested_number = None

    if requested_number is not None and requested_number != detected_number:
        print_dual(color_text(
            f"[NOTE] Requested space group (No. {requested_number}) differs from the "
            f"detected one (No. {detected_number}) -- the given sites likely sit on a "
            "higher-symmetry special position than requested. Still valid, just be aware.",
            'yellow'), f_out)

    final_structure = structure
    reduce_info = None
    if args.reduce:
        print_section("[4] UNIT CELL REDUCTION", f_out)
        try:
            reduced_structure, reduce_sga = reduce_to_unitcell(
                final_structure, args.reduce, symprec=args.symprec, angle_tolerance=args.angle_tolerance)
        except ValueError as e:
            _fail(f"unit cell reduction failed: {e}", f_out)

        output_formula = reduced_structure.composition.reduced_formula
        output_atoms = len(reduced_structure)
        print_dual(f"Mode               : {args.reduce}", f_out)
        print_dual(f"Input formula      : {final_structure.composition.reduced_formula}", f_out)
        print_dual(f"Input atoms        : {len(final_structure)}", f_out)
        print_dual(f"Output formula     : {output_formula}", f_out)
        print_dual(f"Output atoms       : {output_atoms}", f_out)
        print_dual(f"Reduction factor   : {len(final_structure) / output_atoms:.3g}x", f_out)
        if output_atoms == len(final_structure):
            print_dual(color_text(
                f"[NOTE] Structure is already the {args.reduce} cell at this symprec "
                "(no reduction).", 'yellow'), f_out)
        final_structure = reduced_structure
        reduce_info = (output_formula, output_atoms)

    ml_relax_info = None
    if args.ml_relax:
        print_section("[5] ML PRE-RELAXATION (MACE)", f_out)
        print_dual(f"Model           : {model_desc}", f_out)
        print_dual(f"Cell relaxation : {'in cell' if args.ml_relax_cell else 'positions only'}", f_out)

        model_arg = args.custom_model if args.custom_model else args.model
        calc = mace_relax.get_calculator(model_arg)
        for line in mace_relax.describe_model(model_arg, calc):
            print_dual(line, f_out)

        atoms = AseAtomsAdaptor.get_atoms(final_structure)
        atoms.calc = calc
        e0 = atoms.get_potential_energy()
        f0 = float(np.abs(atoms.get_forces()).max())
        a0, b0, c0, _, _, _ = atoms.cell.cellpar()
        vol0 = atoms.get_volume()

        # Every space group is inherently 3D-periodic (no vacuum-padded axis is
        # ever possible from Structure.from_spacegroup), so the cell mask relaxes
        # all 3 axes -- unlike stb-unitcell, which must mask out a vacuum axis for
        # arbitrary user-supplied slabs/wires.
        cell_mask = mace_relax.build_cell_mask((False, False, False)) if args.ml_relax_cell else None
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

        final_structure = AseAtomsAdaptor.get_structure(atoms)
        ml_relax_info = (converged, steps_used, e0, e1)

    print_section("[6] STRUCTURE VALIDATION (post-transform)", f_out)
    frac_coords_final = [site.frac_coords for site in final_structure]
    vacuum_axes_final = kspace.detect_vacuum_axes(frac_coords_final, final_structure.lattice.matrix, VACUUM_GAP_ANG)
    try:
        structure_checks.run_malformation_checks(final_structure, vacuum_axes_final, f_out)
        final_sg_label = core_symmetry.space_group_label(final_structure, symprec=args.symprec)
        print_dual(f"Space group    : {final_sg_label}", f_out)
    except Exception as e:
        print_dual(color_text(f"[WARNING] Structure validation could not complete: {e}", 'yellow'), f_out)

    print_section("[7] SYMMETRY ANALYSIS (BEFORE / AFTER)", f_out)
    print_dual(f"Detailed symmetry analysis (Tolerance: {args.symprec} Ang):", f_out)
    print_dual("Before = raw space-group-expanded structure. After = final structure "
               "(post --reduce/--ml-relax, if given) -- these are EXPECTED to match: "
               "reduction/relaxation preserves the crystal's symmetry, it doesn't "
               "change it.", f_out)
    before_info = core_symmetry.symmetry_summary(structure, args.symprec, VACUUM_GAP_ANG)
    after_info = core_symmetry.symmetry_summary(final_structure, args.symprec, VACUUM_GAP_ANG)
    if "Error" in before_info or "Error" in after_info:
        print_dual(color_text("[WARNING] Symmetry analysis failed for at least one structure.", 'yellow'), f_out)
        print_dual(f"  Before: {before_info.get('Error', 'OK')}", f_out)
        print_dual(f"  After : {after_info.get('Error', 'OK')}", f_out)
    else:
        # No "Layer Group" row here (unlike stb-unitcell/stb-slab/etc.) -- every
        # space group this tool builds from is inherently 3D-periodic, so it would
        # always read "N/A (not 2D-periodic)", never actual information.
        properties = ["Crystal System", "Space Group", "Point Group", "Hall Symbol"]
        rows = [([prop, str(before_info.get(prop, "N/A")), str(after_info.get(prop, "N/A"))], None)
                for prop in properties]
        print_table(["Property", "Before", "After"], rows, f_out)
        if all(before_info.get(prop) == after_info.get(prop) for prop in properties):
            print_dual(color_text("[OK] Before/After match exactly, as expected.", 'green'), f_out)

    print_section("[8] WRITING OUTPUT FILE", f_out)
    header_comment = [
        "Structure built by stb-crystalbuilder from a space group and Wyckoff sites.",
        f"Requested space group: {args.spacegroup} -- detected {detected_symbol} "
        f"(No. {detected_number}).",
        f"Lattice: a={args.a} b={b} c={c} alpha={args.alpha} beta={args.beta} gamma={args.gamma}.",
        "Sites: " + "; ".join(f"{symbol} {coords[0]:.6f} {coords[1]:.6f} {coords[2]:.6f}"
                               for symbol, coords in sites) + ".",
    ]
    if reduce_info is not None:
        output_formula, output_atoms = reduce_info
        header_comment.append(
            f"Reduced to {args.reduce} cell: {len(structure)} -> {output_atoms} atoms "
            f"(reduction factor {len(structure) / output_atoms:.3g}x)."
        )
    if ml_relax_info is not None:
        converged, steps_used, e0, e1 = ml_relax_info
        header_comment.append(
            f"ML pre-relaxed with {model_desc} "
            f"({'converged' if converged else 'NOT converged'} in {steps_used} step(s), "
            f"E = {e1:.6f} eV, delta E = {e1 - e0:+.6f} eV)."
        )

    species_meta = {}
    for symbol in species:
        species_meta = structure_io.ensure_species_id(species_meta, symbol)

    new_structure = structure_io.from_pymatgen(final_structure, species_meta=species_meta)
    structure_io.write_fdf(new_structure, args.output, header_comment=header_comment)
    print_dual(color_text(f"[OK] Structure written to '{args.output}'.", 'green'), f_out)

    print_section("[9] REFERENCES", f_out)
    bib_entries = [citations.SIESTA, citations.SIESTA_RECENT]
    if args.ml_relax:
        bib_entries.append(citations.MACE)
        if not args.custom_model:
            bib_entries.append(citations.MACE_MP)
    citations.write_bib_file(BIB_FILE, bib_entries)
    print_dual(color_text(
        f"[OK] Citations for the methods used in this run written to '{BIB_FILE}' "
        f"({len(bib_entries)} entries).", 'green'), f_out)

    print_section("[10] SUMMARY & FILES", f_out)
    print_dual("Status         : OK", f_out)
    print_dual(f"Output file    : {args.output}", f_out)
    print_dual(f"References     : {BIB_FILE}", f_out)
    if report_path:
        print_dual(f"Report         : {report_path}", f_out)

    if f_out:
        f_out.close()

    # --view runs last, after every check/report section above has already
    # printed, so a blocking GUI window never delays or hides them -- shows
    # both frames (raw built structure vs. final structure) so the user can
    # page through the actual comparison in ase-gui.
    if args.view:
        before_atoms = AseAtomsAdaptor.get_atoms(structure)
        after_atoms = AseAtomsAdaptor.get_atoms(final_structure)
        view_structure_interactive([before_atoms, after_atoms])


if __name__ == "__main__":
    main()
