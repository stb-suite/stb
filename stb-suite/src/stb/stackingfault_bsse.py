#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.0.0"

import os
import sys
import argparse
from stb.core import structure_io, adsorption_sites
from stb.core.bsse import make_ghost_variant, strip_config_extra_include
from stb.core.siesta_log import report_quality_diagnostics
from stb.core.cli import color_text, show_intro, print_dual, print_section, print_table, capture_library_noise
from stb.core.pseudopotentials import copy_pseudo
from stb.stackingfault_analysis import read_grid_manifest, RUN_SUBDIR

REPORT_FILE = "stackingfault_bsse_report.txt"


def read_point_theory_flags(point_dir):
    """Returns vdw_corrected: whether this grid point's own config_extra.fdf
    (written by stb-stackingfault's write_grid_folder) forced the DFT-D3
    dispersion correction for THIS point. A counterpoise correction is only
    meaningful if every fragment is evaluated at the exact same level of
    theory as E_site, not just the same geometry -- same reasoning as
    adsorb_bsse.py's read_site_theory_flags, simplified here since
    stb-stackingfault's own config_extra.fdf never carries Spin/dipole
    blocks (only single-point/z-relax + optional D3) -- there's nothing
    else to inherit. Returns False if the point has no config_extra.fdf at
    all (shouldn't happen for a folder stb-stackingfault actually wrote).
    """
    path = os.path.join(point_dir, adsorption_sites.CONFIG_EXTRA_FILE)
    if not os.path.isfile(path):
        return False
    with open(path) as f:
        text = f.read()
    return "DFTD3" in text


def resolve_point_geometry(point_dir, mode):
    """Returns (geometry, status): the FdfStructure to build ghost
    fragments from, and 'ready'/'not_relaxed'/'mismatched'/'missing'.

    Under --mode 1 (SIESTA relaxes z for real), the counterpoise
    correction is only valid at the ACTUAL relaxed geometry -- a finished
    'siesta.XV' must be present (used_relaxed=True), otherwise 'not_relaxed'
    (not fatal, just skipped -- same as adsorb_bsse.py's own site-readiness
    gate). Under --mode 2/3, the point's own structure.fdf IS ALREADY the
    final geometry SIESTA evaluates (MD.Steps is forced to 0 there, nothing
    ever moves), so no relaxation to wait for -- 'ready' as soon as the
    folder exists, whether or not a '.XV' happens to be present too.
    """
    try:
        geometry, used_relaxed = structure_io.read_relaxed_or_input(point_dir)
    except FileNotFoundError:
        return None, 'missing'
    except ValueError:
        return None, 'mismatched'
    if mode == 1 and not used_relaxed:
        return None, 'not_relaxed'
    return geometry, 'ready'


def write_bsse_folders(bsse_dir, point_fdf, n_layer1_atoms, calc_text, pp_path, vdw_corrected):
    """Writes the two ghost-fragment references a Boys-Bernardi counterpoise
    (BSSE) correction of this grid point needs: '<bsse_dir>/bsse_layer1/'
    (real layer 1 + ghost layer 2) and '<bsse_dir>/bsse_layer2/' (ghost
    layer 1 + real layer 2), both at the SAME geometry as the point itself
    (see resolve_point_geometry). `n_layer1_atoms` is the layer1/layer2
    split index -- already known exactly from stb-stackingfault's own
    manifest (core/heterostructure.py::build_stacked_structure always
    writes layer 1's atoms first, layer 2's after, by construction), no
    guessing needed unlike a generic 2-fragment split.

    Unlike stb-stackingfault's own grid folders, these are ALWAYS
    single-point (adsorption_sites.SINGLE_POINT_BLOCK) regardless of which
    --mode produced the real geometry -- a ghost fragment is evaluated
    once, frozen, at the already-resolved geometry; it never needs (or
    should get) --mode 1's own restricted relaxation block. `vdw_corrected`
    (from read_point_theory_flags) matches the ghost fragments' level of
    theory to the point's own D3 status -- without it, a ghost fragment
    missing pairwise dispersion terms E_site actually includes would be a
    different-level-of-theory bug, same class adsorb_bsse.py's own
    vdw_corrected inheritance closes.
    """
    n_total = len(point_fdf.atoms)
    layer1_variant = make_ghost_variant(point_fdf, n_layer1_atoms, n_total)  # ghost layer 2
    layer2_variant = make_ghost_variant(point_fdf, 0, n_layer1_atoms)        # ghost layer 1

    config_extra_content = adsorption_sites.SINGLE_POINT_BLOCK
    if vdw_corrected:
        config_extra_content += adsorption_sites.VDW_CORRECTION_BLOCK

    for sub_dir, variant in [("bsse_layer1", layer1_variant), ("bsse_layer2", layer2_variant)]:
        out_dir = os.path.join(bsse_dir, sub_dir)
        os.makedirs(out_dir, exist_ok=True)
        structure_io.write_fdf(variant, os.path.join(out_dir, "structure.fdf"))
        with open(os.path.join(out_dir, adsorption_sites.CONFIG_EXTRA_FILE), "w") as f:
            f.write(config_extra_content)
        with open(os.path.join(out_dir, "calc.fdf"), "w") as f:
            f.write(structure_io.prepend_include(calc_text, adsorption_sites.CONFIG_EXTRA_FILE))
        present_labels = sorted({symbol for symbol, _ in variant.atoms})
        for label in present_labels:
            real_symbol = label[:-len("_ghost")] if label.endswith("_ghost") else label
            copy_pseudo(pp_path, real_symbol, out_dir, dest_label=label)


def read_original_calc_text(point_calc_fdf_path):
    """Recovers the original --calc template text from a grid point's own
    calc.fdf, stripping the '%include config_extra.fdf' sidecar
    stb-stackingfault's own write_grid_folder added -- needed so
    write_bsse_folders can prepend its OWN (different) config_extra.fdf
    without doubling the include line. Thin wrapper around
    core.bsse.strip_config_extra_include (shared with adsorb_bsse.py).
    """
    with open(point_calc_fdf_path) as f:
        text = f.read()
    return strip_config_extra_include(text, adsorption_sites.CONFIG_EXTRA_FILE)


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Stage 2 of the 2D Stacking Fault workflow: generates the "
        "BSSE (counterpoise) ghost-fragment folders at each grid point's actual geometry, once "
        "stb-stackingfault's own positions/shift_II_JJ/ folders are ready (relaxed, for --mode 1; "
        "immediately, for --mode 2/3).", 'bold')}""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s\n"
               "  %(prog)s --dir sf_run --file calc.out\n"
    )

    parser.add_argument("--dir", type=str, default=RUN_SUBDIR,
                         help=f"Root directory containing 'positions/shift_II_JJ/' and "
                              f"sf_manifest.json (default: '{RUN_SUBDIR}', stb-stackingfault "
                              "Stage 1's own self-contained run folder). BSSE folders are written "
                              "to this same directory's own 'bsse/' subfolder, a sibling of "
                              "'positions/'.")
    parser.add_argument("--file", type=str, default="calc.out",
                         help="SIESTA output filename inside each grid point folder (default: "
                              "calc.out). Only used for the advisory quality diagnostics below "
                              "-- never gates whether BSSE folders get written.")
    parser.add_argument("--force-tolerance", type=float, default=0.05,
                         help="Residual atomic force in eV/Ang (default: 0.05) above which a "
                              "point's relaxation is flagged as possibly incomplete -- advisory "
                              "only, never blocks BSSE folder generation.")
    parser.add_argument("--save-report", action="store_true",
                         help=f"Also persist the report to <dir>/{REPORT_FILE}. Off by default.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-stackingfaultBsse {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    library_warnings = []  # collected via capture_library_noise, reported in [4]

    manifest_result = read_grid_manifest(args.dir)
    if manifest_result is None:
        print(color_text(
            f"[ERROR] No readable sf_manifest.json found in '{args.dir}'. Did you run "
            "stb-stackingfault? (BSSE needs the manifest's n_layer1_atoms -- it can't be "
            "recovered from folder names alone.)", 'red'))
        sys.exit(1)
    grid_rows, _grid_nx, _grid_ny, _scan, _scan_points, mode, n_layer1_atoms = manifest_result
    if n_layer1_atoms is None:
        print(color_text(
            f"[ERROR] sf_manifest.json in '{args.dir}' has no 'n_layer1_atoms' -- it was written "
            "by an older stb-stackingfault version. Re-run Stage 1 to get a current manifest.",
            'red'))
        sys.exit(1)

    bsse_root = os.path.join(args.dir, "bsse")

    report_path = os.path.join(args.dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    print_dual(color_text("===== STACKING-FAULT BSSE (POST-RUN) REPORT =====", 'magenta'), f_out)

    # --- [0] RUN METADATA. ---
    print_section("[0] RUN METADATA", f_out)
    print_dual(f"  Directory       : {args.dir}", f_out)
    print_dual(f"  Output file     : {args.file}", f_out)
    print_dual(f"  Mode (Stage 1)  : {mode}", f_out)
    print_dual(f"  Save report     : {'yes -> ' + report_path if args.save_report else 'no'}", f_out)
    print_dual(color_text(
        "\n  [WHY THIS STAGE EXISTS] SIESTA's localized (PAO) basis lets layer 1 and layer 2 "
        "artificially lower each other's energy by borrowing basis functions across the gap "
        "(BSSE) -- a real effect for a van-der-Waals-bound interlayer interaction, and one that "
        "does NOT cancel out across the grid: eclipsed/high-symmetry registries have atoms "
        "closer together (more basis overlap, more BSSE) than offset ones, so uncorrected BSSE "
        "can distort the gamma-surface's SHAPE, not just shift every point by a constant. "
        "Because stb-stackingfaultAnalysis only ever reports RELATIVE energies (equilibrium, "
        "corrugation = max - min), no isolated-layer reference is needed here (unlike "
        "stb-adsorbBsse) -- each point's own two ghost-fragment energies are enough; the "
        "constant isolated-layer energy they'd otherwise be compared against cancels in any "
        "relative comparison anyway.", 'cyan'), f_out)

    # --- [1] POINT SCAN: which grid points are ready for BSSE. ---
    print_section("[1] POINT SCAN", f_out)
    ready, not_relaxed, mismatched, missing = [], [], [], []
    rows = []
    for label, _i, _j, _shift_x, _shift_y in grid_rows:
        point_dir = os.path.join(args.dir, "positions", label)
        with capture_library_noise(library_warnings, f"sisl (point {label} .XV read)"):
            geometry, status = resolve_point_geometry(point_dir, mode)
        if status == 'missing':
            missing.append(label)
            rows.append(([label, "SKIP (folder not found)"], 'red'))
            continue
        if status == 'mismatched':
            mismatched.append(label)
            rows.append(([label, "SKIP (atom count/species mismatch)"], 'red'))
            continue
        if status == 'not_relaxed':
            not_relaxed.append(label)
            rows.append(([label, "SKIP (no .XV yet -- not relaxed)"], 'yellow'))
            continue
        ready.append((label, point_dir, geometry))
        rows.append(([label, "ready"], 'green'))
        out_path = os.path.join(point_dir, args.file)
        if os.path.exists(out_path):
            report_quality_diagnostics(label, out_path, args.force_tolerance, f_out)

    print_table(["Point", "Status"], rows, f_out)
    print_dual(f"\n  Grid points found  : {len(grid_rows)}", f_out)
    print_dual(f"  Ready              : {len(ready)}", f_out)
    if mode == 1:
        print_dual(f"  Not yet relaxed    : {len(not_relaxed)}"
                    + (f" ({', '.join(not_relaxed)})" if not_relaxed else ""), f_out)
    if mismatched:
        print_dual(color_text(
            f"  Mismatched/stale   : {len(mismatched)} ({', '.join(mismatched)})", 'red'), f_out)
    if missing:
        print_dual(color_text(
            f"  Missing folders    : {len(missing)} ({', '.join(missing)})", 'red'), f_out)

    if not ready:
        print_dual(color_text(
            "\n[ERROR] No grid point is ready yet"
            + (" -- run SIESTA in 'positions/shift_II_JJ/' first." if mode == 1 else "."),
            'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    # --- [2] WRITING BSSE FOLDERS. ---
    print_section("[2] WRITING BSSE FOLDERS", f_out)
    written = []
    for label, point_dir, geometry in ready:
        calc_text = read_original_calc_text(os.path.join(point_dir, "calc.fdf"))
        vdw_corrected = read_point_theory_flags(point_dir)
        bsse_dir = os.path.join(bsse_root, label)
        # The point's own folder already has every real pseudopotential
        # file it needs (copied there by stb-stackingfault) -- reused
        # directly as the pseudopotential source for the ghost variants
        # too, so this stage never needs its own -p/--pseudo-dir flag.
        write_bsse_folders(bsse_dir, geometry, n_layer1_atoms, calc_text, point_dir, vdw_corrected)
        written.append(label)
        print_dual(f"  {color_text('[OK]', 'green')} {bsse_dir}/bsse_layer1/, {bsse_dir}/bsse_layer2/ "
                    f"(vdw: {'yes' if vdw_corrected else 'no'} -- inherited from the point)", f_out)

    # --- [3] SUMMARY & NEXT STEPS. ---
    print_section("[3] SUMMARY & NEXT STEPS", f_out)
    print_dual(f"  BSSE folder pair(s) written : {len(written)} of {len(grid_rows)} grid point(s)",
                f_out)
    if not_relaxed or mismatched or missing:
        print_dual(color_text(
            f"  {len(not_relaxed) + len(mismatched) + len(missing)} grid point(s) skipped -- "
            "re-run stb-stackingfaultBsse once they're ready (already-written folders are left "
            "untouched).", 'yellow'), f_out)
    print_dual(f"  Run SIESTA in every 'bsse/shift_II_JJ/bsse_layer1/' and "
                "'bsse/shift_II_JJ/bsse_layer2/' folder above, then use stb-stackingfaultAnalysis "
                "(Stage 3) for the BSSE-corrected gamma-surface.", f_out)
    if report_path:
        print_dual(f"  Report: {report_path}", f_out)

    # --- [4] LIBRARY WARNINGS: always last. ---
    print_section("[4] LIBRARY WARNINGS", f_out)
    if library_warnings:
        print_dual(color_text(
            "Messages emitted by external libraries (sisl) during this run -- collected here "
            "instead of interleaved with the report above; harmless in almost every case, but "
            "worth a look if a section above looks suspicious.", 'cyan'), f_out)
        for entry in library_warnings:
            print_dual(entry, f_out)
    else:
        print_dual("No library warnings.", f_out)

    if f_out:
        f_out.close()

    print(f"\n{color_text('Success:', 'green')} {len(written)} of {len(grid_rows)} grid point(s) "
          "got BSSE folders written at their geometry.")
    if report_path:
        print(f"Full report: {report_path}")


if __name__ == "__main__":
    main()
