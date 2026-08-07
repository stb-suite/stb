#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "2.0.0"

import sys
import time
import difflib
import argparse
from datetime import datetime

import numpy as np
from ase.build import molecule
from ase.collections import g2
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

REPORT_FILE = "stb_molecule_report.txt"
BIB_FILE = "references.bib"

# Same default vacuum-gap threshold as stb-unitcell/stb-supercell/stb-fetch/
# stb-kgrid/stb-mlrelax (core/kspace.py's other callers). A molecule built
# here is always vacuum-padded on all 3 axes (0D, isolated) -- unlike
# stb-passivate/stb-crystalbuilder, there is no periodic cell to speak of,
# so symmetry is reported via a POINT group (core_symmetry.point_group_label),
# not a space/layer group.
VACUUM_GAP_ANG = 10.0


def _fail(message, f_out):
    """Prints a red [ERROR] line, closes the report file if one is open, and
    exits with status 1 -- same single error-exit pattern as stb-unitcell/
    stb-supercell/stb-passivate/stb-crystalbuilder's own _fail()."""
    print_dual(color_text(f"[ERROR] {message}", 'red'), f_out)
    if f_out:
        f_out.close()
    sys.exit(1)


def print_available_names():
    names = sorted(g2.names)
    print(color_text(f"Available molecules ({len(names)}):", 'bold'))
    columns = 6
    width = max(len(n) for n in names) + 2
    for i in range(0, len(names), columns):
        row = names[i:i + columns]
        print("  " + "".join(f"{n:<{width}}" for n in row))


def _describe_structure(pmg_structure, vacuum_axes, f_out):
    """Formula/atoms/dimensionality/cell parameters -- same field set as
    stb-supercell/stb-crystalbuilder's own [1] section, minus the atomic
    density line (always skipped here: a molecule's box is vacuum-padded
    on all 3 axes by construction)."""
    print_dual(f"Formula        : {pmg_structure.composition.reduced_formula}", f_out)
    print_dual(f"Atoms          : {len(pmg_structure)}", f_out)
    print_dual(f"Dimensionality : {kspace.dimensionality_label(vacuum_axes)}", f_out)
    a, b, c, alpha, beta, gamma = pmg_structure.lattice.parameters
    print_dual(f"Cell a,b,c     : {a:.4f}, {b:.4f}, {c:.4f} Ang", f_out)


def _validate_structure(pmg_structure, vacuum_axes, f_out):
    """Shared malformation checklist (core.structure_checks) plus a POINT
    -group label -- unlike stb-crystalbuilder/stb-passivate's own
    _validate_structure() (which report a space/layer group), a molecule's
    vacuum box has no physically meaningful periodic cell, so
    core_symmetry.point_group_label (pymatgen's non-periodic
    PointGroupAnalyzer) is used instead, same as stb-fetch's own 0D
    branch."""
    structure_checks.run_malformation_checks(pmg_structure, vacuum_axes, f_out)
    pg_label = core_symmetry.point_group_label(pmg_structure)
    if pg_label:
        print_dual(f"Point group    : {pg_label}", f_out)
    else:
        print_dual(color_text("[WARNING] Point group could not be determined.", 'yellow'), f_out)
    return pg_label


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Builds a reference molecule from ASE's G2 database.", 'bold')}
Places the molecule in a vacuum box, ready as an isolated-molecule
reference (e.g. for stb-cohesive) or a quick starting structure.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s --name H2O\n"
               "  %(prog)s --name CH3OH --vacuum 12 -o methanol.fdf\n"
               "  %(prog)s --name H2O --ml-relax --save-report --view\n"
               "  %(prog)s --list\n"
    )

    parser.add_argument("--name", type=str, default=None,
                        help="Exact (case-sensitive) name from ASE's G2 database, "
                             "e.g. 'H2O', 'CO2', 'C6H6'. See --list for all names. "
                             "Required unless --list is given.")
    parser.add_argument("--list", dest="list_names", action="store_true",
                        help="Print all available molecule names and exit.")
    parser.add_argument("--vacuum", type=float, default=10.0,
                        help="Vacuum padding in Angstrom, added in every direction "
                             "around the molecule (default: 10.0).")

    parser.add_argument("--ml-relax", action="store_true",
                        help="Pre-relax the molecule's geometry with a MACE potential "
                             "(needs the optional 'ml' extra: pip install stb_suite[ml]) "
                             "before writing it out. Positions only -- there is no "
                             "physically meaningful periodic cell to relax for an "
                             "isolated molecule, so no --ml-relax-cell option exists "
                             "here. Off by default.")
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
                        help="Open an interactive 3D view (via ASE) comparing the "
                             "as-built and final molecule (page through frames in "
                             "ase-gui) after writing the output file. Needs a display. "
                             "Off by default.")

    parser.add_argument("-o", "--output", type=str, default="molecule.fdf",
                        help="Output .fdf file name (default: molecule.fdf).")
    parser.add_argument("-v", "--version", action="version", version=f"stb-molecule {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.list_names:
        print_available_names()
        sys.exit(0)

    if args.name is None:
        parser.error("--name is required (unless --list is given).")

    if args.custom_model and not args.ml_relax:
        parser.error("--custom-model is only valid together with --ml-relax.")
    if args.model != "small" and not args.ml_relax:
        parser.error("--model is only valid together with --ml-relax.")
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

    print_dual(color_text("===== STB-MOLECULE REPORT =====", 'magenta'), f_out)

    model_desc = f"a custom model ({args.custom_model})" if args.custom_model else f"MACE-MP-0 ({args.model})"

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time      : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Molecule name  : {args.name}", f_out)
    print_dual(f"Vacuum padding : {args.vacuum} Ang", f_out)
    print_dual(f"Output file    : {args.output}", f_out)
    print_dual(f"ML pre-relax   : {model_desc if args.ml_relax else 'no'}", f_out)

    try:
        atoms = molecule(args.name, vacuum=args.vacuum)
    except KeyError:
        case_insensitive_hit = next(
            (n for n in g2.names if n.lower() == args.name.lower()), None)
        msg = f"'{args.name}' is not a known molecule name (names are case-sensitive)."
        if case_insensitive_hit:
            msg += f" Did you mean '{case_insensitive_hit}'?"
        else:
            suggestions = difflib.get_close_matches(args.name, g2.names, n=3)
            if suggestions:
                msg += f" Did you mean: {', '.join(suggestions)}?"
        msg += " Run with --list to see all available names."
        _fail(msg, f_out)

    pmg_structure = AseAtomsAdaptor.get_structure(atoms)

    print_section("[1] BUILD PARAMETERS", f_out)
    print_dual(f"Requested name : {args.name} (ASE G2 database)", f_out)
    frac_coords_before = pmg_structure.frac_coords
    vacuum_axes_before = kspace.detect_vacuum_axes(frac_coords_before, pmg_structure.lattice.matrix, VACUUM_GAP_ANG)
    _describe_structure(pmg_structure, vacuum_axes_before, f_out)

    print_section("[2] STRUCTURE VALIDATION (pre-transform)", f_out)
    try:
        pg_before = _validate_structure(pmg_structure, vacuum_axes_before, f_out)
    except Exception as e:
        pg_before = None
        print_dual(color_text(f"[WARNING] Structure validation could not complete: {e}", 'yellow'), f_out)

    ml_relax_info = None
    final_pmg = pmg_structure
    if args.ml_relax:
        print_section("[3] ML PRE-RELAXATION (MACE)", f_out)
        print_dual(f"Model           : {model_desc} (device={args.device})", f_out)
        print_dual("Cell relaxation : n/a (isolated molecule, positions only)", f_out)

        model_arg = args.custom_model if args.custom_model else args.model
        try:
            calc = mace_relax.get_calculator(model_arg, device=args.device)
        except ValueError as e:
            _fail(str(e), f_out)
        for line in mace_relax.describe_model(model_arg, calc):
            print_dual(line, f_out)

        atoms.calc = calc
        e0 = atoms.get_potential_energy()
        f0 = float(np.abs(atoms.get_forces()).max())

        t0 = time.time()
        converged, steps_used = mace_relax.relax(atoms, calc, cell_mask=None, fmax=0.05, max_steps=200)
        wall_time = time.time() - t0

        e1 = atoms.get_potential_energy()
        f1 = float(np.abs(atoms.get_forces()).max())

        print_dual(f"Steps used : {steps_used} "
                   f"({'converged' if converged else 'hit step cap, NOT converged'})", f_out)
        print_dual(f"Wall time  : {wall_time:.1f} s", f_out)

        n_atoms = len(atoms)
        rows = [
            (["Energy (eV)", f"{e0:.6f}", f"{e1:.6f}",
              f"{e1 - e0:+.6f} ({(e1 - e0) / n_atoms:+.6f}/atom)"], None),
            (["Max force (eV/Ang)", f"{f0:.4f}", f"{f1:.4f}", f"{f1 - f0:+.4f}"], None),
        ]
        print_table(["Quantity", "Before", "After", "Change"], rows, f_out)

        final_pmg = AseAtomsAdaptor.get_structure(atoms)
        ml_relax_info = (converged, steps_used, e0, e1)

    print_section("[4] STRUCTURE VALIDATION (post-transform)", f_out)
    frac_coords_final = final_pmg.frac_coords
    vacuum_axes_final = kspace.detect_vacuum_axes(frac_coords_final, final_pmg.lattice.matrix, VACUUM_GAP_ANG)
    try:
        pg_after = _validate_structure(final_pmg, vacuum_axes_final, f_out)
    except Exception as e:
        pg_after = None
        print_dual(color_text(f"[WARNING] Structure validation could not complete: {e}", 'yellow'), f_out)

    print_section("[5] SYMMETRY ANALYSIS (BEFORE / AFTER)", f_out)
    print_dual("Point group of the as-built molecule vs. the final (post-ML-relax, if any) "
               "geometry -- an isolated molecule has no periodic space/layer group, so only "
               "the point group applies here.", f_out)
    print_table(["Property", "Before", "After"],
                [(["Point Group", str(pg_before) if pg_before else "N/A", str(pg_after) if pg_after else "N/A"], None)],
                f_out)
    if pg_before and pg_after and pg_before == pg_after:
        print_dual(color_text("[OK] Before/After match exactly.", 'green'), f_out)
    elif pg_before and pg_after:
        print_dual(color_text(
            "[NOTE] Point group changed -- relaxation moved the molecule slightly off its "
            "as-built symmetry (or onto a different one). Not necessarily a problem, just "
            "worth being aware of.", 'yellow'), f_out)

    print_section("[6] WRITING OUTPUT FILE", f_out)
    header_comment = [
        "Reference molecule built by stb-molecule from ASE's G2 database.",
        f"Name: {args.name}. Vacuum padding: {args.vacuum} Ang.",
        f"Point group: {pg_before if pg_before else 'N/A'}.",
    ]
    if ml_relax_info is not None:
        converged, steps_used, e0, e1 = ml_relax_info
        header_comment.append(
            f"ML pre-relaxed with {model_desc} "
            f"({'converged' if converged else 'NOT converged'} in {steps_used} step(s), "
            f"E = {e1:.6f} eV, delta E = {e1 - e0:+.6f} eV)."
        )

    species_meta = {}
    for symbol in dict.fromkeys(site.specie.symbol for site in final_pmg):
        species_meta = structure_io.ensure_species_id(species_meta, symbol)

    new_structure = structure_io.from_pymatgen(final_pmg, species_meta=species_meta, coord_format="cartesian")
    structure_io.write_fdf(new_structure, args.output, header_comment=header_comment)
    print_dual(color_text(f"[OK] Structure written to '{args.output}'.", 'green'), f_out)

    print_section("[7] REFERENCES", f_out)
    bib_entries = [citations.SIESTA, citations.SIESTA_RECENT]
    if args.ml_relax:
        bib_entries.append(citations.MACE)
        if not args.custom_model:
            bib_entries.append(citations.MACE_MP)
    citations.write_bib_file(BIB_FILE, bib_entries)
    print_dual(color_text(
        f"[OK] Citations for the methods used in this run written to '{BIB_FILE}' "
        f"({len(bib_entries)} entries).", 'green'), f_out)

    print_section("[8] SUMMARY & FILES", f_out)
    print_dual("Status         : OK", f_out)
    print_dual(f"Output file    : {args.output}", f_out)
    print_dual(f"References     : {BIB_FILE}", f_out)
    if report_path:
        print_dual(f"Report         : {report_path}", f_out)

    if f_out:
        f_out.close()

    # --view runs last, after every check/report section above has already
    # printed, so a blocking GUI window never delays or hides them -- shows
    # both frames (as-built vs. final molecule) so the user can page
    # through the actual comparison in ase-gui.
    if args.view:
        before_atoms = AseAtomsAdaptor.get_atoms(pmg_structure)
        after_atoms = AseAtomsAdaptor.get_atoms(final_pmg)
        view_structure_interactive([before_atoms, after_atoms])


if __name__ == "__main__":
    main()
