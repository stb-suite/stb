#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "2.1.0"  # multi-cutoff BSSE convergence scan: detects and reports/plots
                    # 'atoms_bsse_check_<cutoff>/' (2+ points, from stb-cohesive's
                    # --bsse-convergence-increment given multiple values), a new
                    # [3b] section, alongside the existing single-point
                    # 'atoms_bsse_check/' comparison (unchanged behavior)

import os
import re
import sys
import glob
import argparse
from stb.core import siesta_log, structure_io
from stb.core.cli import color_text, print_dual, print_section, print_table, show_intro

REPORT_FILE = "cohesive_results.dat"
PLOT_FILE = "cohesive_correction.png"
PLOT_COLORS = ['#2255cc', '#cc5522', '#22aa55']

def get_atom_counts(fdf_path):
    """Number of atoms per chemical species -- thin wrapper around the shared
    core/structure_io.py::atom_counts (moved there once cohesive_energy.py's
    prep side needed the identical logic to decide which isolated-atom
    folders are actually needed)."""
    try:
        structure = structure_io.read_fdf(fdf_path)
    except FileNotFoundError:
        print(color_text(f"[ERROR] Structure file '{fdf_path}' not found.", 'red'))
        sys.exit(1)
    except ValueError as e:
        print(color_text(f"[ERROR] {e}", 'red'))
        sys.exit(1)

    return structure_io.atom_counts(structure)


_SITE_DIR_RE = re.compile(r"^site_.+_x(\d+)$")


def read_bsse_energy(bsse_root, sym, out_file):
    """Reads the BSSE-corrected isolated-atom energy for species `sym` under
    `bsse_root` ('atoms_bsse' or 'atoms_bsse_check'), transparently handling
    both layouts stb-cohesive can produce (see cohesive_energy.py::
    resolve_bsse_sites): a flat 'bsse_root/<sym>/<out_file>' (single site, or
    --no-bsse-multi-site), or a multi-site 'bsse_root/<sym>/site_<wyckoff>_x
    <mult>/<out_file>' per symmetrically distinct site -- averaged, weighted
    by multiplicity, so a species with e.g. 4 octahedral + 2 tetrahedral
    atoms gets a single reference that's actually representative of its real
    population, not just whichever site happened to be generated first.

    Returns (energy, complete): `complete` is False if the flat file, or ANY
    expected multi-site subfolder, is missing/unreadable -- a partial
    weighted average would silently misrepresent the correction, so this
    always drops the whole species rather than average over what's present
    (same "advisory only, never silently wrong" spirit as the rest of the
    suite; the caller decides whether to skip BSSE reporting entirely).
    """
    sym_dir = os.path.join(bsse_root, sym)
    flat_path = os.path.join(sym_dir, out_file)
    if os.path.isfile(flat_path):
        energy = siesta_log.get_free_energy(flat_path)
        return energy, energy is not None

    if not os.path.isdir(sym_dir):
        return None, False

    site_dirs = sorted(
        d for d in os.listdir(sym_dir)
        if _SITE_DIR_RE.match(d) and os.path.isdir(os.path.join(sym_dir, d))
    )
    if not site_dirs:
        return None, False

    total_weight = 0
    weighted_sum = 0.0
    for d in site_dirs:
        mult = int(_SITE_DIR_RE.match(d).group(1))
        e = siesta_log.get_free_energy(os.path.join(sym_dir, d, out_file))
        if e is None:
            return None, False
        weighted_sum += e * mult
        total_weight += mult
    if total_weight == 0:
        return None, False
    return weighted_sum / total_weight, True


def find_reference_cutoff(root_dir):
    """Best-effort recovery of the numeric BSSE cutoff (Ang) used to build
    the ghost-cluster reference(s) under `root_dir` (e.g. 'atoms_bsse/'),
    from the 'bsse_cutoff.txt' sidecar cohesive_energy.py writes alongside
    every reference it generates -- a flat 'atoms_bsse/'-style folder name
    alone carries no cutoff information (unlike the multi-scan
    'atoms_bsse_check_<cutoff>/' naming, where it's already in the folder
    name itself). Returns None if not found/unreadable (e.g. a run from
    before this sidecar file existed) -- advisory only, degrades to "cutoff
    unknown" in the report/plot rather than blocking anything.
    """
    for dirpath, _dirs, files in os.walk(root_dir):
        if "bsse_cutoff.txt" in files:
            try:
                with open(os.path.join(dirpath, "bsse_cutoff.txt")) as f:
                    return float(f.read().strip())
            except (ValueError, OSError):
                return None
    return None


_SCAN_DIR_RE = re.compile(r"atoms_bsse_check_([\d.]+)$")


def discover_scan_points(root_dir):
    """Finds every 'atoms_bsse_check_<cutoff>/' folder directly under
    `root_dir` -- written by stb-cohesive only when --bsse-convergence-
    increment is given 2+ values (a full cutoff scan), as opposed to the
    single flat 'atoms_bsse_check/' this deliberately does NOT match (that
    one-point case is handled separately, unchanged, for backward
    compatibility). Returns a list of (cutoff, dir_path) sorted ascending.
    """
    points = []
    for entry in glob.glob(os.path.join(root_dir, "atoms_bsse_check_*")):
        if not os.path.isdir(entry):
            continue
        m = _SCAN_DIR_RE.search(os.path.basename(entry))
        if m:
            points.append((float(m.group(1)), entry))
    return sorted(points, key=lambda t: t[0])


def plot_bsse_correction(used_species, isolated_energies, bsse_energies, check_energies,
                          e_coh_per_atom, e_coh_per_atom_bsse, e_coh_per_atom_check,
                          out_path, view, scan_curve=None):
    """2- or 3-panel matplotlib figure illustrating the BSSE correction --
    left panel shows WHERE it comes from (per-species isolated-atom
    reference energy SHIFT introduced by each correction, relative to the
    uncorrected value -- not the absolute isolated-atom energies, which are
    typically tens of eV and would dwarf a correction of order 0.1-1 eV on
    the same axis, making it invisible), middle panel shows the NET effect
    (final cohesive energy per atom, absolute, under each method). Built
    entirely from data main() already computed -- no new energy extraction/
    arithmetic here. Only called when a BSSE-corrected reference is
    actually available and complete (see [4] of the report); always saved
    to `out_path`, only shown interactively (plt.show) when `view` is True,
    matching stb-elasticAnalysis's own --view semantics.

    `scan_curve`, if given, is a list of (cutoff, e_coh_per_atom) sorted
    ascending by cutoff (see [3b] BSSE CUTOFF CONVERGENCE SCAN) -- when
    present (2+ points), a third panel plots cohesive energy per atom vs.
    cutoff as a line, the clearest way to see whether the correction has
    actually plateaued or is still trending.
    """
    import numpy as np
    import matplotlib.pyplot as plt

    species = list(used_species.keys())
    methods = []
    values = {}
    if bsse_energies:
        methods.append("BSSE-corrected")
        values["BSSE-corrected"] = [bsse_energies[s] - isolated_energies[s] for s in species]
    if check_energies:
        methods.append("BSSE check")
        values["BSSE check"] = [check_energies[s] - isolated_energies[s] for s in species]

    has_scan = scan_curve and len(scan_curve) >= 2
    if has_scan:
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5))
    else:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    x = np.arange(len(species))
    width = 0.8 / max(len(methods), 1)
    for i, method in enumerate(methods):
        offset = (i - (len(methods) - 1) / 2.0) * width
        # color index +1: keeps PLOT_COLORS[0] (blue) reserved for
        # "Uncorrected" everywhere, matching the right panel's own coloring,
        # even though the uncorrected reference here is just the axhline.
        ax1.bar(x + offset, values[method], width, label=method,
                color=PLOT_COLORS[(i + 1) % len(PLOT_COLORS)])
    ax1.set_xticks(x)
    ax1.set_xticklabels(species)
    ax1.set_ylabel("Isolated-atom energy shift vs. uncorrected (eV)")
    ax1.set_title("BSSE Correction per Species (relative to uncorrected)")
    ax1.axhline(0, color=PLOT_COLORS[0], linestyle='--', linewidth=1.2,
                label="Uncorrected (reference)")
    ax1.grid(True, alpha=0.3, axis='y')
    ax1.legend()

    labels2 = ["Uncorrected"]
    vals2 = [e_coh_per_atom]
    if e_coh_per_atom_bsse is not None:
        labels2.append("BSSE-corrected")
        vals2.append(e_coh_per_atom_bsse)
    if e_coh_per_atom_check is not None:
        labels2.append("BSSE check")
        vals2.append(e_coh_per_atom_check)
    ax2.bar(labels2, vals2, color=PLOT_COLORS[:len(labels2)])
    ax2.set_ylabel("Cohesive energy per atom (eV/atom)")
    ax2.set_title("Cohesive Energy per Atom -- Net Effect")
    ax2.grid(True, alpha=0.3, axis='y')

    if has_scan:
        cutoffs = [c for c, _e in scan_curve]
        energies = [e for _c, e in scan_curve]
        ax3.plot(cutoffs, energies, marker='o', color=PLOT_COLORS[1])
        ax3.set_xlabel("BSSE cutoff (Ang)")
        ax3.set_ylabel("Cohesive energy per atom (eV/atom)")
        ax3.set_title("BSSE Cutoff Convergence Scan")
        ax3.grid(True, alpha=0.3)

    fig.suptitle("BSSE (Counterpoise) Correction")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    if view:
        plt.show()
    else:
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Process cohesive energy results from SIESTA calculations.",
        epilog="Example usage:\n"
               "  stb_cohesive_analysis -o calc.out\n"
               "  stb_cohesive_analysis -o calc.out -d /path/to/results",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("-o", "--out", dest="out_file", type=str, required=True,
                        help="Name of the SIESTA output file (e.g., calc.out)")

    parser.add_argument("-d", "--dir", dest="dir_path", type=str, default="cohesive_runs",
                        help="Path to the folder containing 'structure' and 'atoms' "
                             "directories (default: cohesive_runs, matching stb-cohesive's "
                             "own --output-dir default -- no --dir needed in the common "
                             "case, even run from the parent directory).")

    parser.add_argument("--force-tolerance", dest="force_tolerance", type=float, default=0.05,
                        help="Residual atomic force in eV/Ang (default: 0.05) above which the "
                             "full structure's calc.out is flagged as possibly not relaxed -- "
                             "the cohesive energy would then reflect a strained/off-equilibrium "
                             "geometry rather than the true minimum. Advisory only, never "
                             "blocks the result.")

    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to <dir>/{REPORT_FILE}. Off by default.")

    parser.add_argument("--view", action="store_true",
                        help="Show an interactive matplotlib preview of the BSSE-correction "
                             "figure (see [4] of the report) before finishing. The figure "
                             "itself is always saved when BSSE data is available, regardless "
                             "of this flag. Off by default.")

    parser.add_argument("-v", "--version", action="version", version=f"stb-cohesiveAnalysis {VERSION}")

    parser.add_argument("--no-intro", dest="intro", action="store_false",
                        help="Do not show the introduction")

    args = parser.parse_args()

    if not os.path.isdir(args.dir_path):
        print(color_text(f"[FAIL] Directory '{args.dir_path}' not found.", 'red'))
        sys.exit(1)

    if args.intro == True:
        show_intro([
            "Siesta ToolBox Suite",
            "Cohesive Energy Post-Processing",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    root_dir = args.dir_path
    struct_dir = os.path.join(root_dir, "structure")
    atoms_dir = os.path.join(root_dir, "atoms")
    atoms_bsse_dir = os.path.join(root_dir, "atoms_bsse")
    atoms_bsse_check_dir = os.path.join(root_dir, "atoms_bsse_check")
    bsse_available = os.path.isdir(atoms_bsse_dir)
    check_available = os.path.isdir(atoms_bsse_check_dir)
    scan_points = discover_scan_points(root_dir)
    scan_available = len(scan_points) > 0

    report_path = os.path.join(root_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    print_dual(color_text("===== STB-COHESIVEANALYSIS REPORT (COHESIVE ENERGY) =====", 'magenta'), f_out)

    # --- [0] RUN METADATA ---
    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Results directory      : {root_dir}", f_out)
    print_dual(f"Target log file        : {args.out_file}", f_out)
    print_dual(f"Force tolerance        : {args.force_tolerance} eV/Ang", f_out)
    print_dual(f"BSSE reference          : {'found' if bsse_available else 'not found'}", f_out)
    print_dual(f"BSSE convergence check  : {'found' if check_available else 'not found'}", f_out)
    if scan_available:
        cuts = ", ".join(f"{c:.1f}" for c, _d in scan_points)
        print_dual(f"BSSE cutoff scan        : found, {len(scan_points)} point(s) ({cuts} Ang)", f_out)
    print_dual(f"Save report             : {'yes' if args.save_report else 'no'}", f_out)
    print_dual(f"View (matplotlib)       : {'yes' if args.view else 'no'}", f_out)
    if report_path:
        print_dual(f"Report file             : {report_path}", f_out)

    # --- [1] INPUT STRUCTURE ---
    print_section("[1] INPUT STRUCTURE", f_out)
    if not os.path.exists(struct_dir) or not os.path.exists(atoms_dir):
        print_dual(color_text(
            f"[FAIL] Required directories 'structure' and/or 'atoms' not found in "
            f"'{root_dir}'.", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    struct_fdf = os.path.join(struct_dir, "structure.fdf")
    atom_counts = get_atom_counts(struct_fdf)
    if not atom_counts:
        print_dual(color_text(
            "[FAIL] Could not extract atom counts. Check your structure.fdf", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    total_atoms = sum(atom_counts.values())
    if total_atoms == 0:
        print_dual(color_text(
            "[FAIL] Structure file declares species but has zero atoms. Check your "
            "structure.fdf", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    composition = ", ".join(f"{sym}{n}" for sym, n in atom_counts.items())
    print_dual("Parsed from structure.fdf in 'structure/'. Each species with a nonzero count "
               "needs its own isolated-atom reference in 'atoms/' (see [2]).", f_out)
    print_table(["Quantity", "Value"], [
        (["Composition", composition], None),
        (["Total atoms", str(total_atoms)], None),
    ], f_out)
    used_species = {sym: count for sym, count in atom_counts.items() if count > 0}
    unused_species = [sym for sym, count in atom_counts.items() if count == 0]
    if unused_species:
        print_dual(color_text(
            f"[INFO] Declared but never placed (no isolated-atom calculation needed): "
            f"{', '.join(unused_species)}", 'cyan'), f_out)

    # --- [2] ENERGY EXTRACTION ---
    print_section("[2] ENERGY EXTRACTION", f_out)
    errors = False

    struct_out_path = os.path.join(struct_dir, args.out_file)
    e_bulk = siesta_log.get_free_energy(struct_out_path)
    if e_bulk is None:
        print_dual(color_text(
            f"[WARNING] Could not find '{args.out_file}' or finished calculation for the "
            "full structure.", 'yellow'), f_out)
        errors = True
    else:
        print_dual(f"Full structure energy   : {e_bulk:12.6f} eV", f_out)

        # Numerical-quality diagnostics (best-effort, never blocks): SCF
        # convergence and residual forces, checked ONLY on the full
        # structure -- an isolated atom has zero net force by symmetry, and
        # the BSSE cluster's real atom only feels Pulay-type forces from the
        # ghosts' incomplete basis (not a physical force to "relax" away),
        # so neither is a meaningful place to ask "is this geometry actually
        # at its minimum".
        converged, iterations = siesta_log.get_scf_convergence(struct_out_path)
        if not converged:
            print_dual(color_text(
                f"[WARNING] Could not confirm SCF convergence for the full structure "
                f"('{args.out_file}') -- this result may be unreliable.", 'yellow'), f_out)
        else:
            print_dual(f"SCF convergence          : converged after {iterations} iteration(s)", f_out)

        max_force = siesta_log.get_max_force(struct_out_path)
        if max_force is not None:
            if max_force > args.force_tolerance:
                print_dual(color_text(
                    f"[WARNING] Residual force on the full structure ({max_force:.4f} eV/Ang) "
                    f"exceeds --force-tolerance ({args.force_tolerance} eV/Ang) -- this "
                    "geometry may not be relaxed, so the cohesive energy may reflect a "
                    "strained/off-equilibrium structure rather than the true minimum.",
                    'yellow'), f_out)
            else:
                print_dual(f"Residual force           : {max_force:.4f} eV/Ang (within "
                           f"--force-tolerance {args.force_tolerance} eV/Ang)", f_out)

    # Isolated-atom energies: uncorrected (always attempted), BSSE-corrected
    # and BSSE-check (best-effort, auto-detected from atoms_bsse[_check]/
    # written by stb-cohesive --bsse-correction/--bsse-convergence-check).
    isolated_energies = {}
    sum_isolated_energy = 0.0
    for sym, count in used_species.items():
        sym_dir = os.path.join(atoms_dir, sym)
        e_atom = siesta_log.get_free_energy(os.path.join(sym_dir, args.out_file))
        if e_atom is None:
            print_dual(color_text(
                f"[WARNING] Could not find '{args.out_file}' or results for isolated atom: "
                f"{sym}", 'yellow'), f_out)
            errors = True
        else:
            isolated_energies[sym] = e_atom
            sum_isolated_energy += (e_atom * count)

    bsse_energies = {}
    bsse_complete = bsse_available
    if bsse_available:
        for sym in used_species:
            e_atom_bsse, ok = read_bsse_energy(atoms_bsse_dir, sym, args.out_file)
            if not ok:
                print_dual(color_text(
                    f"[WARNING] Could not find '{args.out_file}' or results for the "
                    f"BSSE-corrected isolated atom: {sym} -- BSSE-corrected cohesive energy "
                    "will be skipped.", 'yellow'), f_out)
                bsse_complete = False
            else:
                bsse_energies[sym] = e_atom_bsse

    check_energies = {}
    check_complete = check_available
    if check_available:
        for sym in used_species:
            e_atom_check, ok = read_bsse_energy(atoms_bsse_check_dir, sym, args.out_file)
            if not ok:
                print_dual(color_text(
                    f"[WARNING] Could not find '{args.out_file}' or results for the BSSE "
                    f"convergence-check isolated atom: {sym} -- the convergence check will "
                    "be skipped.", 'yellow'), f_out)
                check_complete = False
            else:
                check_energies[sym] = e_atom_check

    # Multi-cutoff BSSE convergence scan (best-effort, never blocks): reads
    # every 'atoms_bsse_check_<cutoff>/' point found, dropping (with a
    # warning) any point missing a species -- same "advisory only, whole
    # point dropped rather than a silently wrong partial result" spirit as
    # bsse_energies/check_energies above.
    scan_energies = []  # list of (cutoff, {sym: energy})
    for cutoff, scan_dir in scan_points:
        point_energies = {}
        point_complete = True
        for sym in used_species:
            e_atom_scan, ok = read_bsse_energy(scan_dir, sym, args.out_file)
            if not ok:
                print_dual(color_text(
                    f"[WARNING] Could not find '{args.out_file}' or results for the BSSE "
                    f"convergence-scan reference at cutoff {cutoff:.1f} Ang for: {sym} -- this "
                    "scan point will be skipped.", 'yellow'), f_out)
                point_complete = False
            else:
                point_energies[sym] = e_atom_scan
        if point_complete:
            scan_energies.append((cutoff, point_energies))

    print_dual(
        "\nIsolated-atom (one atom in vacuum) reference energies used in the cohesive-energy "
        "sum below. The 'Delta' columns show the shift each correction introduces relative to "
        "the uncorrected reference -- BSSE removes an artificial over-stabilization of the "
        "isolated atom (it's computed with a smaller effective basis than it has inside the "
        "solid), so a positive Delta here is the expected, physically correct direction.",
        f_out)
    energy_rows = []
    for sym in used_species:
        delta_bsse = f"{bsse_energies[sym] - isolated_energies[sym]:+.6f}" if sym in bsse_energies else "--"
        delta_check = f"{check_energies[sym] - isolated_energies[sym]:+.6f}" if sym in check_energies else "--"
        row = [
            sym,
            f"{isolated_energies[sym]:.6f}" if sym in isolated_energies else "--",
            f"{bsse_energies[sym]:.6f}" if sym in bsse_energies else "--",
            delta_bsse,
            f"{check_energies[sym]:.6f}" if sym in check_energies else "--",
            delta_check,
        ]
        energy_rows.append((row, 'yellow' if sym not in isolated_energies else None))
    print_table(["Species", "Uncorrected (eV)", "BSSE-corrected (eV)", "Delta BSSE (eV)",
                 "BSSE check (eV)", "Delta check (eV)"],
                energy_rows, f_out)

    if errors:
        print_dual(color_text(
            "\n[FAIL] Cannot calculate cohesive energy because some calculations are "
            "missing or incomplete.", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    # --- [3] COHESIVE ENERGY RESULTS ---
    print_section("[3] COHESIVE ENERGY RESULTS", f_out)
    e_coh_total = e_bulk - sum_isolated_energy
    e_coh_per_atom = e_coh_total / total_atoms

    e_coh_per_atom_bsse = None
    if bsse_available and bsse_complete:
        sum_isolated_bsse = sum(bsse_energies[sym] * count for sym, count in used_species.items())
        e_coh_total_bsse = e_bulk - sum_isolated_bsse
        e_coh_per_atom_bsse = e_coh_total_bsse / total_atoms

    e_coh_per_atom_check = None
    if check_available and check_complete and bsse_available and bsse_complete:
        sum_isolated_check = sum(check_energies[sym] * count for sym, count in used_species.items())
        e_coh_total_check = e_bulk - sum_isolated_check
        e_coh_per_atom_check = e_coh_total_check / total_atoms

    print_dual(f"Sum of isolated atoms   : {sum_isolated_energy:12.6f} eV", f_out)
    print_dual(f"Bulk structure energy   : {e_bulk:12.6f} eV", f_out)
    print_dual(
        "\nCohesive energy = (Bulk structure energy - Sum of isolated atoms) / Total atoms. "
        "More negative means a more strongly bound structure. 'Delta vs Uncorrected' isolates "
        "the effect of each correction alone on the final per-atom result.", f_out)

    result_rows = [(["Uncorrected", f"{e_coh_total:.4f}", f"{e_coh_per_atom:.4f}", "(reference)"], None)]
    if e_coh_per_atom_bsse is not None:
        delta_per_atom = e_coh_per_atom_bsse - e_coh_per_atom
        result_rows.append((["BSSE-corrected", f"{e_coh_total_bsse:.4f}", f"{e_coh_per_atom_bsse:.4f}",
                             f"{delta_per_atom:+.4f}"], None))
    if e_coh_per_atom_check is not None:
        delta_check_per_atom = e_coh_per_atom_check - e_coh_per_atom
        result_rows.append((["BSSE check (larger cutoff)", f"{e_coh_total_check:.4f}",
                             f"{e_coh_per_atom_check:.4f}", f"{delta_check_per_atom:+.4f}"], None))
    print_table(["Method", "Total (eV)", "Per atom (eV/atom)", "Delta vs Uncorrected (eV/atom)"],
                result_rows, f_out)

    if e_coh_per_atom_bsse is not None:
        print_dual(f"BSSE correction per atom: {delta_per_atom:+.4f} eV/atom (uncorrected LCAO "
                   "cohesive energies systematically over-bind -- expect this to make the "
                   "energy less negative)", f_out)
    elif not bsse_available:
        print_dual(color_text(
            "[NOTE] This result is NOT corrected for BSSE (Basis Set Superposition Error) -- "
            "a known bias of LCAO cohesive energies: the isolated atom is computed with a "
            "smaller effective basis than it has inside the solid (no neighbors' orbitals "
            "available), which systematically over-binds. Re-run stb-cohesive with "
            "--bsse-correction (the CLI default) for a BSSE-corrected reference.", 'yellow'), f_out)

    if e_coh_per_atom_check is not None:
        check_delta = e_coh_per_atom_check - e_coh_per_atom_bsse
        print_dual(f"BSSE cutoff convergence shift: {check_delta:+.4f} eV/atom (difference "
                   "between --bsse-cutoff and --bsse-cutoff + --bsse-convergence-increment -- "
                   "small means the BSSE correction is converged w.r.t. cutoff; large means "
                   "increase --bsse-cutoff further)", f_out)

    # --- [3b] BSSE CUTOFF CONVERGENCE SCAN ---
    scan_curve = []
    if scan_available:
        print_section("[3b] BSSE CUTOFF CONVERGENCE SCAN", f_out)
        if e_coh_per_atom_bsse is None:
            print_dual(color_text(
                "[INFO] Skipped -- no complete base BSSE-corrected reference to compare "
                "against (see [2]/[3]).", 'cyan'), f_out)
        else:
            print_dual(
                "Cohesive energy per atom (BSSE-corrected) at each scanned cutoff -- from "
                "stb-cohesive's --bsse-convergence-increment given multiple values. The "
                "correction is converged once these stop changing appreciably as the cutoff "
                "grows; if it's still trending by the last point, increase --bsse-cutoff/"
                "--bsse-convergence-increment further.", f_out)

            base_cutoff = find_reference_cutoff(atoms_bsse_dir)
            scan_rows = []
            if base_cutoff is not None:
                scan_curve.append((base_cutoff, e_coh_per_atom_bsse))
                scan_rows.append(([f"{base_cutoff:.1f} Ang (base)", f"{e_coh_per_atom_bsse:.4f}",
                                   "(reference)"], None))
            else:
                print_dual(color_text(
                    "[INFO] Base --bsse-cutoff value not recoverable from 'atoms_bsse/' (an "
                    "older run, generated before this sidecar file existed) -- the scan below "
                    "is shown on its own, without the base point for reference.", 'cyan'), f_out)

            for cutoff, point_energies in scan_energies:
                sum_scan = sum(point_energies[sym] * count for sym, count in used_species.items())
                e_coh_scan = (e_bulk - sum_scan) / total_atoms
                scan_curve.append((cutoff, e_coh_scan))
                delta_base = e_coh_scan - e_coh_per_atom_bsse
                scan_rows.append(([f"{cutoff:.1f} Ang", f"{e_coh_scan:.4f}", f"{delta_base:+.4f}"], None))

            print_table(["Cutoff", "Per atom (eV/atom)", "Delta vs base BSSE (eV/atom)"],
                        scan_rows, f_out)

            if len(scan_curve) >= 2:
                (prev_cutoff, prev_e), (last_cutoff, last_e) = scan_curve[-2], scan_curve[-1]
                print_dual(f"Last step change: {last_e - prev_e:+.4f} eV/atom (cutoff "
                           f"{prev_cutoff:.1f} -> {last_cutoff:.1f} Ang) -- small means "
                           "converged, large means keep increasing the cutoff.", f_out)

    # --- [4] CORRECTION PLOT ---
    print_section("[4] CORRECTION PLOT", f_out)
    plot_path = None
    if bsse_available and bsse_complete:
        plot_path = os.path.join(root_dir, PLOT_FILE)
        plot_bsse_correction(used_species, isolated_energies, bsse_energies, check_energies,
                              e_coh_per_atom, e_coh_per_atom_bsse, e_coh_per_atom_check,
                              plot_path, args.view, scan_curve=scan_curve)
        print_dual(color_text(f"[OK] Saved: {plot_path}", 'green')
                   + (" (shown interactively via --view)" if args.view else ""), f_out)
        if not args.view:
            print_dual("Pass --view to also show this plot interactively before finishing.", f_out)
    else:
        print_dual(color_text(
            "[INFO] Skipped -- no complete BSSE-corrected data available to compare against "
            "the uncorrected reference (see [2]/[3]).", 'cyan'), f_out)

    # --- [5] SUMMARY & FILES ---
    print_section("[5] SUMMARY & FILES", f_out)
    print_dual(f"Cohesive energy per atom (uncorrected)  : {e_coh_per_atom:.4f} eV/atom", f_out)
    if e_coh_per_atom_bsse is not None:
        print_dual(f"Cohesive energy per atom (BSSE-corrected): {e_coh_per_atom_bsse:.4f} eV/atom", f_out)
    if report_path:
        print_dual(f"Report      : {report_path}", f_out)
    else:
        print_dual("Pass --save-report to persist this report to a file.", f_out)
    if plot_path:
        print_dual(f"Plot        : {plot_path}", f_out)

    if f_out:
        f_out.close()

if __name__ == "__main__":
    main()
