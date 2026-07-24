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
from pymatgen.core.surface import SlabGenerator
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

REPORT_FILE = "stb_slab_report.txt"
BIB_FILE = "references.bib"

# Same default vacuum-gap threshold as stb-fetch/stb-kgrid/stb-mlrelax
# (core/kspace.py's other callers), used to detect vacuum-padded axes on the
# BULK input. A cut slab's own threshold instead tracks whatever vacuum the
# user actually asked for (see min(VACUUM_GAP_ANG, args.min_vacuum_size) at
# its use site below) -- a user-requested vacuum smaller than 10 Ang would
# otherwise not register as "vacuum-padded" at all.
VACUUM_GAP_ANG = 10.0


def slab_metrics(slab):
    """(thickness, vacuum) in Angstrom. reorient_lattice=True (SlabGenerator's
    own default, kept here) puts the vacuum-normal on the c-axis, so this
    matches what kspace.detect_vacuum_axes will later find on the written file.
    """
    cart_z = slab.cart_coords[:, 2]
    thickness = float(cart_z.max() - cart_z.min())
    vacuum = float(slab.lattice.c - thickness)
    return thickness, vacuum


def sort_key(slab, symprec):
    """Non-polar and symmetric terminations first -- pymatgen's own get_slabs()
    order carries no physical meaning, so index 0 needs a meaningful default.
    """
    return (slab.is_polar(), not slab.is_symmetric(symprec=symprec))


def print_termination_table(slabs, symprec):
    """Candidate-termination table for -i/--interactive's manual pick. Prints
    straight to the console, not threaded through the persisted report --
    same convention/rationale as stacking2D.py's own ZSL match table: it's an
    interactive decision aid, not a fact about the run's final outcome (which
    IS captured in the report via the chosen termination's own summary).
    """
    header = (f"{'ID':<4} | {'Formula':<12} | {'Atoms':<6} | {'Thickness(A)':<13} | "
              f"{'Vacuum(A)':<10} | {'Polar':<6} | {'Symmetric':<9}")
    print(color_text(header, 'blue'))
    print("-" * len(header))
    for i, slab in enumerate(slabs):
        thickness, vacuum = slab_metrics(slab)
        polar = slab.is_polar()
        symmetric = slab.is_symmetric(symprec=symprec)
        row = (f"{i:<4} | {slab.composition.reduced_formula:<12} | {len(slab):<6} | "
               f"{thickness:<13.2f} | {vacuum:<10.2f} | "
               f"{'Yes' if polar else 'No':<6} | {'Yes' if symmetric else 'No':<9}")
        if not polar and symmetric:
            print(color_text(row, 'green'))
        else:
            print(row)


def prompt_termination(slabs, hkl, symprec):
    print("\n" + color_text(f"--- Terminations for hkl={hkl} ({len(slabs)} found) ---", 'bold'))
    print_termination_table(slabs, symprec)
    while True:
        user_input = input(color_text(
            f"\nSelect a termination ID (0 to {len(slabs) - 1}) or 'q' to quit: ", 'yellow'))
        if user_input.lower() == 'q':
            sys.exit(0)
        try:
            selected = int(user_input)
            if 0 <= selected < len(slabs):
                return selected
        except ValueError:
            pass


def print_slab_summary(slab, hkl, polar, symmetric, f_out):
    """`polar`/`symmetric` are passed in rather than recomputed here, since
    they're Slab-only properties (Slab.is_polar()/is_symmetric()) -- after
    --ml-relax, `slab` may be a plain pymatgen Structure (AseAtomsAdaptor.
    get_structure() doesn't reconstruct the Slab subclass), so the caller
    caches these from the original, pre-relax Slab object instead."""
    thickness, vacuum = slab_metrics(slab)
    print_dual(f"Formula          : {slab.composition.reduced_formula}", f_out)
    print_dual(f"Total atoms      : {len(slab)}", f_out)
    print_dual(f"Miller index     : {hkl}", f_out)
    print_dual(f"Slab thickness   : {thickness:.2f} Ang", f_out)
    print_dual(f"Vacuum thickness : {vacuum:.2f} Ang", f_out)
    print_dual(f"Polar            : {'Yes' if polar else 'No'}", f_out)
    print_dual(f"Symmetric        : {'Yes' if symmetric else 'No'}", f_out)


def _validate_structure(pmg_structure, vacuum_axes, f_out):
    """Shared malformation checklist (core.structure_checks) plus a
    space-group label -- same shape as mlrelax.py/supercell.py's own
    _validate_structure(), wrapped in try/except by the caller (a validation
    failure is reported, never fatal)."""
    structure_checks.run_malformation_checks(pmg_structure, vacuum_axes, f_out)
    sg_label = core_symmetry.space_group_label(pmg_structure)
    print_dual(f"Space group      : {sg_label}", f_out)
    return sg_label


def write_slab(slab, species_meta, out_path, header_comment):
    new_structure = structure_io.from_pymatgen(slab, species_meta=species_meta)
    structure_io.write_fdf(new_structure, out_path, header_comment=header_comment)


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Cuts a Miller-index slab from a bulk SIESTA FDF structure.", 'bold')}
Adds vacuum along the surface normal; the vacuum axis lands on c.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s -f bulk.fdf --hkl 1 0 0\n"
               "  %(prog)s -f bulk.fdf --hkl 1 1 1 -i\n"
               "  %(prog)s -f bulk.fdf --hkl 1 1 1 --all -o surface.fdf\n"
               "  %(prog)s -f bulk.fdf --hkl 1 1 1 --passivate\n"
               "  %(prog)s -f bulk.fdf --hkl 1 1 1 --ml-relax --save-report --view\n"
    )

    parser.add_argument("-f", "--file", dest="filename", type=str, required=True,
                        help="Path to the input bulk structure file (.fdf).")
    parser.add_argument("--hkl", type=int, nargs=3, required=True, metavar=("H", "K", "L"),
                        help="Miller index of the surface to cut, e.g. --hkl 1 0 0")
    parser.add_argument("--min-slab-size", type=float, default=10.0,
                        help="Minimum slab thickness in Angstrom (default: 10.0).")
    parser.add_argument("--min-vacuum-size", type=float, default=15.0,
                        help="Minimum vacuum thickness in Angstrom (default: 15.0).")
    parser.add_argument("--lll-reduce", action="store_true",
                        help="LLL-reduce the slab lattice (pymatgen SlabGenerator option).")
    parser.add_argument("--center-slab", action="store_true",
                        help="Center the slab in the middle of the vacuum (pymatgen SlabGenerator option).")
    parser.add_argument("--primitive", action="store_true",
                        help="Reduce the bulk to its primitive cell before cutting the slab. "
                             "Off by default so the input cell is used as given.")
    parser.add_argument("--symmetrize", action="store_true",
                        help="Ask pymatgen to symmetrize polar/asymmetric terminations.")
    parser.add_argument("--symprec", type=float, default=0.1,
                        help="Symmetry tolerance, Ang, used for the polar/symmetric diagnostics "
                             "AND the before/after symmetry table (default: 0.1).")

    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("-i", "--interactive", action="store_true",
                           help="Show all terminations found and prompt for one to keep.")
    selection.add_argument("--termination", type=int, default=None,
                                help="Directly select a termination index (0-based). "
                                     "Default: index 0 after sorting non-polar/symmetric terminations first.")
    parser.add_argument("--all", action="store_true",
                        help="Write every termination found, instead of just one.")

    parser.add_argument("--passivate", action="store_true",
                        help="Cap dangling bonds on each written termination with a passivating "
                             "atom (same logic as stb-passivate). Only single-missing-bond sites "
                             "are auto-passivated; sites missing 2+ bonds are reported instead of "
                             "guessed.")
    parser.add_argument("--passivant", type=str, default="H",
                        help="Element to cap dangling bonds with, only with --passivate (default: H).")
    parser.add_argument("--cutoff", type=float, default=None,
                        help="Neighbor-search radius in Angstrom for --passivate. "
                             "Default: auto-detected (see stb-passivate --help).")
    parser.add_argument("--bond-length", type=float, default=None,
                        help="Passivant bond length in Angstrom for --passivate. "
                             "Default: auto per species pair (see stb-passivate --help).")

    parser.add_argument("--ml-relax", action="store_true",
                        help="Pre-relax each written slab with a MACE potential before writing it "
                             "out (needs the optional 'ml' extra: pip install stb_suite[ml]) -- "
                             "positions only by default. Off by default.")
    parser.add_argument("--ml-relax-cell", action="store_true",
                        help="With --ml-relax, also relax the in-plane cell -- the vacuum axis "
                             "always stays exactly fixed. Only valid together with --ml-relax.")
    parser.add_argument("--model", choices=["small", "medium", "large"], default="small",
                        help="MACE-MP-0 foundation model size for --ml-relax (default: small).")
    parser.add_argument("--custom-model", default=None, metavar="PATH",
                        help="Path to a custom fine-tuned .model file for --ml-relax, instead of "
                             "a MACE-MP-0 foundation size.")

    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the full run report (including the symmetry analysis) "
                             f"to {REPORT_FILE}. Off by default.")
    parser.add_argument("--view", action="store_true",
                        help="Open an interactive 3D view (via ASE) of the bulk structure and "
                             "every written slab (page through frames in ase-gui) after writing "
                             "the output file(s). Needs a display. Off by default.")

    parser.add_argument("-o", "--output", type=str, default="slab.fdf",
                        help="Output .fdf file name (default: slab.fdf). When more than one "
                             "termination is written, it is suffixed '_term0', '_term1', ...")
    parser.add_argument("-v", "--version", action="version", version=f"stb-slab {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if not args.passivate:
        if args.passivant != "H":
            parser.error("--passivant is only valid with --passivate.")
        if args.cutoff is not None:
            parser.error("--cutoff is only valid with --passivate.")
        if args.bond_length is not None:
            parser.error("--bond-length is only valid with --passivate.")
    else:
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

    def fail(message):
        print_dual(color_text(f"[ERROR] {message}", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_dual(color_text("===== STB-SLAB REPORT =====", 'magenta'), f_out)

    model_desc = f"a custom model ({args.custom_model})" if args.custom_model else f"MACE-MP-0 ({args.model})"

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time        : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Input file       : {args.filename}", f_out)
    print_dual(f"Miller index     : {tuple(args.hkl)}", f_out)
    print_dual(f"Min. slab size   : {args.min_slab_size} Ang", f_out)
    print_dual(f"Min. vacuum size : {args.min_vacuum_size} Ang", f_out)
    print_dual(f"Passivate        : {'yes (' + args.passivant + ')' if args.passivate else 'no'}", f_out)
    print_dual(f"Output file      : {args.output}", f_out)
    print_dual(f"ML pre-relax     : {model_desc if args.ml_relax else 'no'}", f_out)
    if args.ml_relax:
        print_dual(f"Relax cell       : {'yes (vacuum axis fixed)' if args.ml_relax_cell else 'no (positions only)'}", f_out)

    if not os.path.exists(args.filename):
        fail(f"File '{args.filename}' not found.")

    hkl = tuple(args.hkl)
    if hkl == (0, 0, 0):
        fail("Miller index (0, 0, 0) is not valid.")

    try:
        structure = structure_io.read_fdf(args.filename)
    except (FileNotFoundError, ValueError) as e:
        fail(str(e))

    pmg_structure = structure_io.to_pymatgen(structure)

    print_section("[1] BULK STRUCTURE (INPUT)", f_out)
    frac_coords_bulk = [site.frac_coords for site in pmg_structure]
    vacuum_axes_bulk = kspace.detect_vacuum_axes(frac_coords_bulk, pmg_structure.lattice.matrix, VACUUM_GAP_ANG)
    print_dual(f"Formula          : {pmg_structure.composition.reduced_formula}", f_out)
    print_dual(f"Atoms            : {len(pmg_structure)}", f_out)
    print_dual(f"Dimensionality   : {kspace.dimensionality_label(vacuum_axes_bulk)}", f_out)
    a, b, c, alpha, beta, gamma = pmg_structure.lattice.parameters
    print_dual(f"Cell a,b,c       : {a:.4f}, {b:.4f}, {c:.4f} Ang", f_out)
    print_dual(f"Cell angles      : {alpha:.2f}, {beta:.2f}, {gamma:.2f} deg", f_out)

    print_section("[2] STRUCTURE VALIDATION (bulk, pre-cut)", f_out)
    try:
        _validate_structure(pmg_structure, vacuum_axes_bulk, f_out)
    except Exception as e:
        print_dual(color_text(f"[WARNING] Structure validation could not complete: {e}", 'yellow'), f_out)

    print_section("[3] SLAB GENERATION (TERMINATION SEARCH)", f_out)
    print_dual(f"LLL-reduce       : {'yes' if args.lll_reduce else 'no'}", f_out)
    print_dual(f"Center slab      : {'yes' if args.center_slab else 'no'}", f_out)
    print_dual(f"Reduce to primitive first : {'yes' if args.primitive else 'no'}", f_out)
    print_dual(f"Symmetrize       : {'yes' if args.symmetrize else 'no'}", f_out)
    try:
        gen = SlabGenerator(
            pmg_structure, hkl, args.min_slab_size, args.min_vacuum_size,
            lll_reduce=args.lll_reduce, center_slab=args.center_slab, primitive=args.primitive,
        )
        slabs = gen.get_slabs(symmetrize=args.symmetrize)
    except Exception as e:
        fail(f"Could not generate slabs for hkl={hkl}: {e}")

    if not slabs:
        fail(f"No slabs could be generated for Miller index {hkl}.")

    slabs.sort(key=lambda s: sort_key(s, args.symprec))
    print_dual(f"Terminations found: {len(slabs)}", f_out)

    if args.interactive:
        chosen_indices = [prompt_termination(slabs, hkl, args.symprec)]
    elif args.termination is not None:
        if not (0 <= args.termination < len(slabs)):
            fail(f"--termination {args.termination} is out of range (0 to {len(slabs) - 1}).")
        chosen_indices = [args.termination]
    elif args.all:
        chosen_indices = list(range(len(slabs)))
    else:
        chosen_indices = [0]
    print_dual(f"Selected termination(s): {chosen_indices}", f_out)

    base, ext = os.path.splitext(args.output)
    if not ext:
        ext = ".fdf"
    multi = len(chosen_indices) > 1

    # Vacuum-gap threshold for the CUT slabs: the smaller of the generic
    # 10 Ang default and whatever vacuum the user actually asked for -- a
    # requested --min-vacuum-size below 10 Ang would otherwise not register
    # as vacuum-padded at all.
    slab_vacuum_gap = min(VACUUM_GAP_ANG, args.min_vacuum_size)

    slab_results = {}
    passivation_reports = {}
    # polar/symmetric are Slab-only properties (Slab.is_polar()/
    # is_symmetric()) -- cached here from the original Slab object, since
    # --ml-relax later replaces slab_results[idx] with a plain pymatgen
    # Structure (AseAtomsAdaptor.get_structure() doesn't reconstruct the
    # Slab subclass), which would otherwise raise AttributeError below.
    slab_flags = {}
    for idx in chosen_indices:
        slab = slabs[idx]
        slab_flags[idx] = (slab.is_polar(), slab.is_symmetric(symprec=args.symprec))
        species_meta = structure.species_meta
        if args.passivate:
            slab, report = passivate_dangling_bonds(
                slab, passivant=args.passivant, cutoff=args.cutoff, bond_length=args.bond_length)
            passivation_reports[idx] = report
            species_meta = structure_io.ensure_species_id(dict(species_meta), args.passivant)
        slab_results[idx] = (slab, species_meta)

    ml_relax_info = {}
    if args.ml_relax:
        print_section("[4] ML PRE-RELAXATION (MACE)", f_out)
        print_dual(f"Model           : {model_desc}", f_out)
        print_dual(f"Cell relaxation : "
                   f"{'in-plane only (vacuum axis fixed)' if args.ml_relax_cell else 'positions only'}", f_out)
        model_arg = args.custom_model if args.custom_model else args.model
        calc = mace_relax.get_calculator(model_arg)
        for line in mace_relax.describe_model(model_arg, calc):
            print_dual(line, f_out)

        for idx in chosen_indices:
            slab, species_meta = slab_results[idx]
            print_dual(f"\n{color_text(f'--- [Termination {idx}] ---', 'bold')}", f_out)
            frac_coords_slab = [site.frac_coords for site in slab]
            vacuum_axes_slab = kspace.detect_vacuum_axes(frac_coords_slab, slab.lattice.matrix, slab_vacuum_gap)
            n_atoms = len(slab)
            atoms = AseAtomsAdaptor.get_atoms(slab)
            atoms.calc = calc
            e0 = atoms.get_potential_energy()
            f0 = float(np.abs(atoms.get_forces()).max())
            a0, b0, _ = atoms.cell.cellpar()[:3]

            cell_mask = mace_relax.build_cell_mask(vacuum_axes_slab) if args.ml_relax_cell else None
            t0 = time.time()
            converged, steps_used = mace_relax.relax(atoms, calc, cell_mask=cell_mask, fmax=0.05, max_steps=200)
            wall_time = time.time() - t0

            e1 = atoms.get_potential_energy()
            f1 = float(np.abs(atoms.get_forces()).max())
            a1, b1, _ = atoms.cell.cellpar()[:3]
            relaxed = AseAtomsAdaptor.get_structure(atoms)

            slab_results[idx] = (relaxed, species_meta)
            ml_relax_info[idx] = (converged, steps_used, e0, e1)

            print_dual(f"Steps used : {steps_used} "
                       f"({'converged' if converged else 'hit step cap, NOT converged'})", f_out)
            print_dual(f"Wall time  : {wall_time:.1f} s", f_out)

            rows = [
                (["Energy (eV)", f"{e0:.6f}", f"{e1:.6f}",
                  f"{e1 - e0:+.6f} ({(e1 - e0) / n_atoms:+.6f}/atom)"], None),
                (["Max force (eV/Ang)", f"{f0:.4f}", f"{f1:.4f}", f"{f1 - f0:+.4f}"], None),
            ]
            if args.ml_relax_cell:
                rows.append((["Lattice a, b (Ang)", f"{a0:.4f}, {b0:.4f}", f"{a1:.4f}, {b1:.4f}",
                              f"max {100 * max(abs(a1 - a0) / a0, abs(b1 - b0) / b0):+.2f}%"], None))
            print_table(["Quantity", "Before", "After", "Change"], rows, f_out)

    print_section("[5] STRUCTURE VALIDATION, SYMMETRY & WRITING OUTPUT FILE(S)", f_out)
    output_files = []
    view_atoms = [AseAtomsAdaptor.get_atoms(pmg_structure)]

    for idx in chosen_indices:
        slab, species_meta = slab_results[idx]
        polar, symmetric = slab_flags[idx]
        print_dual(f"\n{color_text(f'--- [Termination {idx}] ---', 'bold')}", f_out)
        print_slab_summary(slab, hkl, polar, symmetric, f_out)

        if idx in passivation_reports:
            report = passivation_reports[idx]
            n_passivated = len(report["passivated"])
            n_unresolved = len(report["unresolved"])
            print_dual(f"Dangling bonds found : {n_passivated + n_unresolved}", f_out)
            print_dual(f"Auto-passivated      : {n_passivated} with {args.passivant}", f_out)
            if n_unresolved:
                print_dual(color_text(
                    f"[WARNING] {n_unresolved} atom(s) missing 2+ bonds -- left unpassivated "
                    "(geometrically underdetermined from local coordination alone):", 'yellow'), f_out)
                for atom_idx, symbol, pos, deficit in report["unresolved"]:
                    print_dual(f"    #{atom_idx + 1:<4} {symbol:<3} deficit={deficit}  at {pos}", f_out)

        frac_coords_slab = [site.frac_coords for site in slab]
        vacuum_axes_slab = kspace.detect_vacuum_axes(frac_coords_slab, slab.lattice.matrix, slab_vacuum_gap)
        try:
            _validate_structure(slab, vacuum_axes_slab, f_out)
        except Exception as e:
            print_dual(color_text(f"[WARNING] Structure validation could not complete: {e}", 'yellow'), f_out)

        print_dual(f"Detailed symmetry analysis (Tolerance: {args.symprec} Ang):", f_out)
        before_info = core_symmetry.symmetry_summary(pmg_structure, args.symprec, VACUUM_GAP_ANG)
        after_info = core_symmetry.symmetry_summary(slab, args.symprec, slab_vacuum_gap)
        if "Error" in before_info or "Error" in after_info:
            print_dual(color_text("[WARNING] Symmetry analysis failed for at least one structure.", 'yellow'), f_out)
            print_dual(f"  Bulk (before)     : {before_info.get('Error', 'OK')}", f_out)
            print_dual(f"  Slab (after)      : {after_info.get('Error', 'OK')}", f_out)
        else:
            properties = ["Crystal System", "Space Group", "Layer Group", "Point Group", "Hall Symbol"]
            rows = [([prop, str(before_info.get(prop, "N/A")), str(after_info.get(prop, "N/A"))], None)
                    for prop in properties]
            print_table(["Property", "Bulk (before)", "Slab (after)"], rows, f_out)

        out_path = f"{base}_term{idx}{ext}" if multi else args.output
        header_comment = [
            f"Slab cut by stb-slab from {args.filename}, Miller index {hkl}.",
            f"Termination {idx} of {len(slabs)} found "
            f"({'polar' if polar else 'non-polar'}, "
            f"{'symmetric' if symmetric else 'asymmetric'}).",
        ]
        if idx in passivation_reports:
            header_comment.append(
                f"Passivated {len(passivation_reports[idx]['passivated'])} dangling bond(s) "
                f"with {args.passivant}.")
        if idx in ml_relax_info:
            converged, steps_used, e0, e1 = ml_relax_info[idx]
            header_comment.append(
                f"ML pre-relaxed with {model_desc} "
                f"({'converged' if converged else 'NOT converged'} in {steps_used} step(s), "
                f"E = {e1:.6f} eV, delta E = {e1 - e0:+.6f} eV)."
            )
        write_slab(slab, species_meta, out_path, header_comment)
        print_dual(color_text(f"[OK] Slab written to '{out_path}'.", 'green'), f_out)
        output_files.append(out_path)
        view_atoms.append(AseAtomsAdaptor.get_atoms(slab))

    print_section("[6] REFERENCES", f_out)
    bib_entries = [citations.SIESTA, citations.SIESTA_RECENT]
    if args.ml_relax:
        bib_entries.append(citations.MACE)
        if not args.custom_model:
            bib_entries.append(citations.MACE_MP)
    citations.write_bib_file(BIB_FILE, bib_entries)
    print_dual(color_text(
        f"[OK] Citations for the methods used in this run written to '{BIB_FILE}' "
        f"({len(bib_entries)} entries).", 'green'), f_out)

    print_section("[7] SUMMARY & FILES", f_out)
    print_dual("Status                    : OK", f_out)
    print_dual(f"Terminations written      : {len(output_files)}", f_out)
    shown_files = output_files[:5]
    for of in shown_files:
        print_dual(f"Output file               : {of}", f_out)
    if len(output_files) > len(shown_files):
        print_dual(f"                            ... and {len(output_files) - len(shown_files)} more", f_out)
    print_dual(f"References                : {BIB_FILE}", f_out)
    if report_path:
        print_dual(f"Report                    : {report_path}", f_out)

    if f_out:
        f_out.close()

    # --view runs last, after every check/report section above has already
    # printed, so a blocking GUI window never delays or hides them -- shows
    # the bulk structure and every written slab so the user can page through
    # the actual comparison in ase-gui.
    if args.view:
        view_structure_interactive(view_atoms)


if __name__ == "__main__":
    main()
