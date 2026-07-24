#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""
STB Monolayer Stacker (ZSL Heteroepitaxy Mode)
Reads two Siesta .fdf files, finds a commensurate supercell using the ZSL algorithm,
and stacks them into a van der Waals heterostructure.
"""

VERSION = "2.2.0"

import os
import sys
import time
import warnings
import argparse
import textwrap
import numpy as np
from time import sleep
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from stb.core import citations, kspace, mace_relax, structure_checks, structure_io
from stb.core.ase_view import view_structure_interactive
from stb.core.deps import require_mace
from stb.core.heterostructure import find_zsl_match, build_stacked_structure
from stb.core.symmetry import layer_group_label

# Suppress Pymatgen warnings for cleaner CLI output
warnings.filterwarnings("ignore")

# ANSI Colors for terminal
from stb.core.cli import COLORS, color_text, print_dual, print_section, print_table, show_intro

REPORT_FILE = "stb_2dstacking_report.txt"
BIB_FILE = "references.bib"

# Same default vacuum-gap threshold as stb-fetch/stb-kgrid (core/kspace.py's
# other callers) -- used to detect the out-of-plane vacuum axis for the
# structure-validation checklist below.
VACUUM_GAP_ANG = 10.0

# Only consumer of this citation in the suite so far -- kept local rather
# than moved into core/citations.py (same "extract on 2nd/3rd use" policy
# documented there).
_BIB_ZSL = ("Zur1984", """@article{Zur1984,
  author  = {Zur, A. and McGill, T. C.},
  title   = {Lattice match: An application to heteroepitaxy},
  journal = {Journal of Applied Physics},
  year    = {1984},
  volume  = {55},
  number  = {2},
  pages   = {378--386},
  doi     = {10.1063/1.333084}
}""")


def _fail(message, f_out):
    """Prints a red [ERROR] line, closes the report file if one is open,
    and exits with status 1 -- same single error-exit pattern as
    dftu.py/fetch.py's own _fail()."""
    print_dual(color_text(f"[ERROR] {message}", 'red'), f_out)
    if f_out:
        f_out.close()
    sys.exit(1)

def parse_fdf_to_pymatgen(filepath, f_out=None):
    """Parses a Siesta .fdf file and returns a Pymatgen Structure object."""
    try:
        return structure_io.to_pymatgen(structure_io.read_fdf(filepath))
    except Exception as e:
        _fail(f"Error parsing {filepath}: {e}", f_out)

def stack_heterostructure(layer1, layer2, gaps, target_vacuum, max_area, max_strain, shift_x=0.0, shift_y=0.0, twist_angle=0.0, strain_mode='top', match_id=0, interactive=False, batch_sym=False):
    """Matches lattices, handles selection, and generates heterostructures.

    Thin wrapper over core.heterostructure's find_zsl_match() (the ZSL
    search, run once) + build_stacked_structure() (called once per
    shift/gap combination requested) -- same public behavior as before
    this was split out, moved once stackingfault.py needed the same two
    pieces separately (run the search once, sweep many more shift/gap
    combinations than batch_sym's fixed 6 points cheaply).
    """
    t_mat1, t_mat2, best_match_data = find_zsl_match(
        layer1, layer2, max_area=max_area, max_strain=max_strain,
        match_id=match_id, interactive=interactive, twist_angle=twist_angle)

    # Define Shifts
    shifts_to_run = {}
    if batch_sym:
        print(color_text("[INFO] Batch Symmetry Mode active. Generating high-symmetry configurations...", 'cyan'))
        if shift_x != 0.0 or shift_y != 0.0:
            print(color_text(f"  -> [WARNING] Custom shifts (-tx {shift_x}, -ty {shift_y}) are IGNORED in --batch_sym mode.", 'yellow'))
        if twist_angle != 0.0:
            print(color_text(f"  -> [WARNING] Using --batch_sym with a twist angle ({twist_angle}°). Shifting may not yield true AA/AB symmetries due to Moiré patterns.", 'yellow'))

        gamma = layer1.lattice.angles[2]
        is_hex = (abs(gamma - 120.0) < 2.0 or abs(gamma - 60.0) < 2.0) and abs(layer1.lattice.a - layer1.lattice.b) < 0.1

        shifts_to_run = {"AA": [0.0, 0.0], "Bridge_X": [0.5, 0.0], "Bridge_Y": [0.0, 0.5], "Center": [0.5, 0.5]}
        if is_hex:
            print(color_text("  -> [INFO] Hexagonal lattice detected. Adding Hex_AB and Hex_BA stackings.", 'green'))
            shifts_to_run.update({"Hex_AB": [1/3, 2/3], "Hex_BA": [2/3, 1/3]})
    else:
        shifts_to_run = {"Custom": [shift_x, shift_y]}

    results = {}
    print(color_text("[INFO] Building heterostructures...", 'cyan'))
    for name, shift in shifts_to_run.items():
        for gap in gaps:
            hetero, _n_layer1_atoms, max_strain_val = build_stacked_structure(
                layer1, layer2, t_mat1, t_mat2, shift[0], shift[1], gap,
                target_vacuum=target_vacuum, strain_mode=strain_mode)

            res_name = f"{name}_gap_{gap:.2f}" if len(gaps) > 1 else name
            results[res_name] = (hetero, max_strain_val)

    return results

def analyze_and_report_symmetry(layer1, layer2, hetero, symprec=0.01, f_out=None):
    """Prints a symmetry comparison table (crystal system, 3D space group of
    the vacuum-padded cell, layer group, point group, Hall symbol) for
    layer1/layer2/the heterostructure -- always folded into the tool's own
    main numbered report (f_out), never a separate file, via the shared
    core.cli.print_table (same convention as core/structure_checks.py's
    validation checklist and stb-mlrelax's before/after table)."""
    def get_details(struct):
        try:
            sga = SpacegroupAnalyzer(struct, symprec=symprec)
            details = {
                "Space Group": f"{sga.get_space_group_symbol()} ({sga.get_space_group_number()})",
                "Point Group": sga.get_point_group_symbol(),
                "Crystal System": str(sga.get_crystal_system()).title(),
                "Hall Symbol": sga.get_hall(),
            }
        except Exception as e:
            return {"Error": f"Analysis failed ({e})"}

        # Space Group above is the 3D space group of the vacuum-padded cell --
        # valid, but not the physically correct classification for a genuinely
        # 2D-periodic structure. Layer Group (spglib's get_layergroup(), same
        # call stb-fetch already uses) accounts for the vacuum direction
        # properly. Additive: never replaces Space Group above.
        try:
            vacuum_axes = kspace.detect_vacuum_axes(struct.frac_coords, struct.lattice.matrix, VACUUM_GAP_ANG)
            if sum(vacuum_axes) == 1:
                label = layer_group_label(struct, vacuum_axes.index(True), symprec)
                details["Layer Group"] = label if label else "N/A (needs spglib >= 2.1.0)"
            else:
                details["Layer Group"] = "N/A (not 2D-periodic)"
        except Exception:
            details["Layer Group"] = "N/A"

        return details

    l1_info = get_details(layer1)
    l2_info = get_details(layer2)
    het_info = get_details(hetero)

    print_dual(f"Detailed symmetry analysis (Tolerance: {symprec} Ang):", f_out)

    if "Error" in l1_info or "Error" in l2_info or "Error" in het_info:
        print_dual(color_text("[WARNING] Symmetry analysis failed for at least one structure.", 'yellow'), f_out)
        print_dual(f"  Layer 1        : {l1_info.get('Error', 'OK')}", f_out)
        print_dual(f"  Layer 2        : {l2_info.get('Error', 'OK')}", f_out)
        print_dual(f"  Heterostructure: {het_info.get('Error', 'OK')}", f_out)
        return

    properties = ["Crystal System", "Space Group", "Layer Group", "Point Group", "Hall Symbol"]
    rows = [([prop, str(l1_info.get(prop, "N/A")), str(l2_info.get(prop, "N/A")), str(het_info.get(prop, "N/A"))], None)
            for prop in properties]
    print_table(["Property", "Layer 1", "Layer 2", "Heterostructure"], rows, f_out)

def main():
    parser = argparse.ArgumentParser(
        description="STB Monolayer Stacking Tool with Geometric Control & Symmetry Analysis.",
        epilog="Example usage:\n"
               "  stb_stacking -l1 bottom.fdf -l2 top.fdf -i\n"
               "  stb_stacking -l1 bottom.fdf -l2 top.fdf --batch_sym -g 3.2",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument("-l1", "--layer1", required=True, help="Path to the bottom monolayer .fdf file")
    parser.add_argument("-l2", "--layer2", required=True, help="Path to the top monolayer .fdf file")
    
    parser.add_argument("-a", "--max_area", type=float, default=150.0, help="Maximum supercell area (default: 150.0)")
    parser.add_argument("-s", "--max_strain", type=float, default=0.05, help="Maximum allowed strain fraction (default: 0.05)")
    
    match_group = parser.add_mutually_exclusive_group()
    match_group.add_argument("-i", "--interactive", action="store_true", help="Interactive mode: shows table to select match ID")
    match_group.add_argument("-id", "--match_id", type=int, default=0, help="Directly select a specific match ID (default: 0)")
    
    parser.add_argument("--batch_sym", action="store_true", help="Generate all high-symmetry stackings automatically")
    
    gap_group = parser.add_mutually_exclusive_group()
    gap_group.add_argument("-g", "--gap", type=float, default=3.2, help="Van der Waals gap in Angstroms (default: 3.2)")
    gap_group.add_argument("--gap_range", type=float, nargs=3, metavar=('START', 'END', 'STEPS'), help="Range for energy curve. Ex: 2.5 4.5 11")
    
    parser.add_argument("--vacuum", type=float, default=None, help="Target vacuum space in Angstroms. Inherits L1 by default.")
    parser.add_argument("-t", "--twist", type=float, default=0.0, help="Twist angle of layer 2 in degrees (default: 0.0)")
    parser.add_argument("-tx", "--shift_x", type=float, default=0.0, help="Fractional shift of layer 2 in X (default: 0.0)")
    parser.add_argument("-ty", "--shift_y", type=float, default=0.0, help="Fractional shift of layer 2 in Y (default: 0.0)")
    parser.add_argument("-sm", "--strain_mode", choices=['top', 'bottom', 'sym'], default='top', help="Strain distribution mode (default: top)")
    
    parser.add_argument("-o", "--output", default="stacked_structure.fdf", help="Output .fdf base file name")
    parser.add_argument("-sp", "--symprec", type=float, default=0.01, help="Symmetry tolerance in Angstroms (default: 0.01)")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the full run report (including the symmetry analysis) "
                             f"to {REPORT_FILE}. Off by default.")
    parser.add_argument("--view", action="store_true",
                        help="Open an interactive 3D view of the final heterostructure(s) (via ASE) "
                             "after writing the output file(s). Needs a display. Off by default.")

    parser.add_argument("--ml-relax", action="store_true",
                        help="Pre-relax each generated heterostructure with a MACE potential "
                             "(needs the optional 'ml' extra: pip install stb_suite[ml]) before "
                             "writing it out -- positions only by default. Off by default.")
    parser.add_argument("--ml-relax-cell", action="store_true",
                        help="With --ml-relax, also relax the in-plane cell -- the vacuum axis "
                             "(c) always stays exactly fixed. Only valid together with --ml-relax.")
    parser.add_argument("--model", choices=["small", "medium", "large"], default="small",
                        help="MACE-MP-0 foundation model size for --ml-relax (default: small).")
    parser.add_argument("--custom-model", default=None, metavar="PATH",
                        help="Path to a custom fine-tuned .model file for --ml-relax, instead of "
                             "a MACE-MP-0 foundation size.")

    parser.add_argument("-v", "--version", action="version", version=f"stb-2Dstacking {VERSION}")
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
            "Siesta ToolBox Suite - Monolayer Stacker",
            "ZSL Heteroepitaxy Algorithm & Geometric Control",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    report_path = REPORT_FILE if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    print_dual(color_text("===== STB-2DSTACKING REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Layer 1 (bottom): {args.layer1}", f_out)
    print_dual(f"Layer 2 (top)   : {args.layer2}", f_out)
    print_dual(f"Max area        : {args.max_area} A^2", f_out)
    print_dual(f"Max strain      : {args.max_strain:.2%}", f_out)
    if args.gap_range:
        print_dual(f"Gap             : range {args.gap_range[0]}-{args.gap_range[1]} A, "
                   f"{int(args.gap_range[2])} point(s)", f_out)
    else:
        print_dual(f"Gap             : {args.gap} A", f_out)
    if args.batch_sym:
        stacking_desc = "batch high-symmetry configurations (--batch_sym)"
    else:
        stacking_desc = f"custom (twist={args.twist} deg, shift=({args.shift_x}, {args.shift_y}))"
    print_dual(f"Stacking mode   : {stacking_desc}", f_out)
    print_dual(f"Strain mode     : {args.strain_mode}", f_out)

    print_section("[1] PARSING INPUT STRUCTURES", f_out)
    layer1 = parse_fdf_to_pymatgen(args.layer1, f_out)
    layer2 = parse_fdf_to_pymatgen(args.layer2, f_out)
    print_dual(f"Layer 1 formula : {layer1.composition.reduced_formula} ({len(layer1)} atoms)", f_out)
    print_dual(f"Layer 2 formula : {layer2.composition.reduced_formula} ({len(layer2)} atoms)", f_out)

    if args.gap_range:
        gaps = np.linspace(args.gap_range[0], args.gap_range[1], int(args.gap_range[2]))
        print_dual(color_text(
            f"[INFO] Energy Curve Mode detected. Generating {len(gaps)} distance points.", 'cyan'), f_out)
    else:
        gaps = [args.gap]

    print_section("[2] LATTICE MATCHING & HETEROSTRUCTURE CONSTRUCTION (ZSL)", f_out)
    print_dual(f"Searching commensurate supercells (max_area={args.max_area} A^2, "
               f"max_strain={args.max_strain:.2%})...", f_out)
    # find_zsl_match()/build_stacked_structure() (via stack_heterostructure()) print their
    # own colored progress/match-table output directly to the console -- shared with
    # stackingfault.py, so it isn't threaded through f_out here (out of scope for this
    # tool alone to change). The parameters and final per-result outcome are still
    # captured in the report via the summary lines below.
    results = stack_heterostructure(
        layer1, layer2, gaps=gaps, target_vacuum=args.vacuum,
        max_area=args.max_area, max_strain=args.max_strain,
        shift_x=args.shift_x, shift_y=args.shift_y, twist_angle=args.twist,
        strain_mode=args.strain_mode, match_id=args.match_id,
        interactive=args.interactive, batch_sym=args.batch_sym
    )
    print_dual(f"Heterostructure(s) built: {len(results)}", f_out)

    model_desc = f"a custom model ({args.custom_model})" if args.custom_model else f"MACE-MP-0 ({args.model})"
    ml_relax_info = {}
    if args.ml_relax:
        print_section("[3] ML PRE-RELAXATION (MACE)", f_out)
        print_dual(f"Model           : {model_desc}", f_out)
        print_dual(f"Cell relaxation : "
                   f"{'in-plane only (vacuum axis fixed)' if args.ml_relax_cell else 'positions only'}", f_out)
        from pymatgen.io.ase import AseAtomsAdaptor
        model_arg = args.custom_model if args.custom_model else args.model
        calc = mace_relax.get_calculator(model_arg)
        for line in mace_relax.describe_model(model_arg, calc):
            print_dual(line, f_out)

        for name, (hetero, applied_strain) in results.items():
            print_dual(f"\n{color_text(f'--- [{name}] ---', 'bold')}", f_out)
            vacuum_axes = kspace.detect_vacuum_axes(hetero.frac_coords, hetero.lattice.matrix, VACUUM_GAP_ANG)
            n_atoms = len(hetero)
            atoms = AseAtomsAdaptor.get_atoms(hetero)
            atoms.calc = calc
            e0 = atoms.get_potential_energy()
            f0 = float(np.abs(atoms.get_forces()).max())
            a0, b0, c0, _, _, gamma0 = atoms.cell.cellpar()

            cell_mask = mace_relax.build_cell_mask(vacuum_axes) if args.ml_relax_cell else None
            t0 = time.time()
            converged, steps_used = mace_relax.relax(atoms, calc, cell_mask=cell_mask, fmax=0.05, max_steps=200)
            wall_time = time.time() - t0

            e1 = atoms.get_potential_energy()
            f1 = float(np.abs(atoms.get_forces()).max())
            a1, b1, c1, _, _, gamma1 = atoms.cell.cellpar()
            relaxed = AseAtomsAdaptor.get_structure(atoms)

            results[name] = (relaxed, applied_strain)
            ml_relax_info[name] = (converged, steps_used, e1, e1 - e0)

            print_dual(f"Steps used : {steps_used} "
                       f"({'converged' if converged else 'hit step cap, NOT converged'})", f_out)
            print_dual(f"Wall time  : {wall_time:.1f} s", f_out)

            # Before/after table, same convention as stb-mlrelax's own
            # [5] BEFORE / AFTER COMPARISON section.
            rows = [
                (["Energy (eV)", f"{e0:.6f}", f"{e1:.6f}",
                  f"{e1 - e0:+.6f} ({(e1 - e0) / n_atoms:+.6f}/atom)"], None),
                (["Max force (eV/Ang)", f"{f0:.4f}", f"{f1:.4f}", f"{f1 - f0:+.4f}"], None),
            ]
            if args.ml_relax_cell:
                rel_change = max(abs(a1 - a0) / a0, abs(b1 - b0) / b0)
                rows.append((["Lattice a, b (Ang)", f"{a0:.4f}, {b0:.4f}", f"{a1:.4f}, {b1:.4f}",
                              f"max {100 * rel_change:+.2f}%"], None))
                rows.append((["Lattice gamma (deg)", f"{gamma0:.2f}", f"{gamma1:.2f}",
                              f"{gamma1 - gamma0:+.2f}"], None))
                rows.append((["Vacuum axis c (Ang)", f"{c0:.4f}", f"{c1:.4f}",
                              "fixed" if c1 == c0 else f"{c1 - c0:+.4f} (SHOULD be 0)"], None))
            print_table(["Quantity", "Before", "After", "Change"], rows, f_out)

    base_out_name = args.output[:-4] if args.output.endswith('.fdf') else args.output
    is_multi_output = args.batch_sym or len(gaps) > 1

    print_section("[4] STRUCTURE VALIDATION & WRITING OUTPUT FILE(S)", f_out)
    output_files = []
    hetero_atoms = []

    for name, (hetero, applied_strain) in results.items():
        print_dual(f"\n{color_text(f'--- [{name}] ---', 'bold')}", f_out)

        print_dual(f"Formula                  : {hetero.composition.reduced_formula}", f_out)
        print_dual(f"Total atoms              : {len(hetero)}", f_out)
        print_dual(f"Lattice C Parameter      : {hetero.lattice.c:.4f} Ang (Vacuum Corrected)", f_out)
        print_dual(f"Max Applied Linear Strain: {applied_strain:.2%}", f_out)

        # Run the full shared validation checklist and the full symmetry table for
        # EVERY generated structure, not just a representative one -- both are short
        # enough (a handful of rows each) that showing them for every --batch_sym/
        # --gap_range result is more useful than a suppressed/summarized view.
        vacuum_axes = kspace.detect_vacuum_axes(hetero.frac_coords, hetero.lattice.matrix, VACUUM_GAP_ANG)
        structure_checks.run_malformation_checks(hetero, vacuum_axes, f_out)
        analyze_and_report_symmetry(layer1, layer2, hetero, symprec=args.symprec, f_out=f_out)

        if is_multi_output:
            out_filename = f"{base_out_name}_{name}.fdf"
        else:
            out_filename = args.output

        header_comment = [
            "Heterostructure built by stb-2Dstacking (ZSL lattice matching).",
            f"Layer 1 (bottom): {args.layer1}",
            f"Layer 2 (top)   : {args.layer2}",
            f"Stacking configuration: {name} (twist={args.twist} deg, strain_mode={args.strain_mode}).",
            f"Max applied linear strain: {applied_strain:.4%}.",
        ]
        if name in ml_relax_info:
            converged, steps_used, energy, delta_e = ml_relax_info[name]
            header_comment.append(
                f"ML pre-relaxed with {model_desc} "
                f"({'converged' if converged else 'NOT converged'} in {steps_used} step(s), "
                f"E = {energy:.6f} eV, delta E = {delta_e:+.6f} eV)."
            )
        new_structure = structure_io.from_pymatgen(hetero)
        structure_io.write_fdf(new_structure, out_filename, header_comment=header_comment)
        print_dual(color_text(f"[OK] Structure written to '{out_filename}'.", 'green'), f_out)
        output_files.append(out_filename)
        hetero_atoms.append(hetero)

    print_section("[5] REFERENCES", f_out)
    bib_entries = [citations.SIESTA, citations.SIESTA_RECENT, _BIB_ZSL]
    if args.ml_relax:
        bib_entries.append(citations.MACE)
        if not args.custom_model:
            bib_entries.append(citations.MACE_MP)
    citations.write_bib_file(BIB_FILE, bib_entries)
    print_dual(color_text(
        f"[OK] Citations for the methods used in this run written to '{BIB_FILE}' "
        f"({len(bib_entries)} entries).", 'green'), f_out)

    print_section("[6] SUMMARY & FILES", f_out)
    print_dual("Status                    : OK", f_out)
    print_dual(f"Heterostructures generated: {len(results)}", f_out)
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

    # --view runs last, after every check/report section above has already printed, so
    # a blocking GUI window never delays or hides them.
    if args.view:
        from pymatgen.io.ase import AseAtomsAdaptor
        atoms_list = [AseAtomsAdaptor.get_atoms(h) for h in hetero_atoms]
        view_structure_interactive(atoms_list)


if __name__ == "__main__":
    main()
