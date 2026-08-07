#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""Stage 2 of 2: reads back the finished SIESTA runs from stb-convergence's
own 'convergence_<parameter>_<value>' folders and reports where each swept
parameter (Mesh.CutOff, PAO.EnergyShift, k-grid density) converges.

Folder discovery (discover_parameter_sweeps) auto-detects BOTH layouts
stb-convergence can produce: its default nested layout, <dir>/<parameter>/
convergence_<parameter>_<value>/ (one subfolder per swept parameter -- so
pointing --dir at the top-level output directory picks up every sweep from
one stb-convergence run in a single pass), and a flat layout, <dir>/
convergence_<parameter>_<value>/ directly (pointing --dir at one
parameter's own subfolder, or an older/hand-built sweep). Every parameter
found is analyzed independently and reported in its own numbered
subsection -- this replaces an earlier version of this tool that could
only ever look at exactly one parameter per run (majority-voting across a
mixed folder listing when more than one was present), which is no longer
needed now that the two tools' layouts actually match.

Since stb-convergence (Stage 1) fully relaxes the structure -- both atomic
positions AND the cell -- at every swept value (see its own module
docstring for why: SIESTA's real-space grid "eggbox effect" can bias a
relaxation's converged geometry at a coarse setting even while the SCF
energy at that arbitrary geometry looks fine), checking only the total
energy here would not actually verify Stage 1's own claim. So each sweep
is judged by TWO independent criteria, both read straight from calc.out:
  - Energy: |Delta E/atom| between consecutive values in the
    more-accurate-ward direction, against --tolerance (eV/atom).
  - Structure: |Delta V/V| of the relaxed cell volume (core.siesta_log.
    get_outcell, the actual relaxed cell SIESTA used, not a re-derivation
    from the folder's own structure.fdf), against --volume-tolerance (%).
The recommended value takes whichever of the two demands the tighter
(more accurate-ward) setting, so it satisfies both criteria at once, and a
[NOTE] flags it whenever they disagree -- exactly the scenario stb
-convergence's own docstring warns about (energy looking converged before
the relaxed geometry actually has).

--force-tolerance is NOT a third convergence criterion -- it's a
per-run data-quality flag (residual force still above tolerance, or the
CG relaxation ran its full declared step budget without visibly
stopping early) warning that a given run's own geometry might not be a
genuine relaxed minimum, independent of whether its energy/volume happen
to look converged against neighboring points.

Report/plot style matches the rest of the suite's "stb-standard report"
convention (strain_analysis.py/elastic_analysis.py/etc.): a numbered
[0]...[N] report (always printed), --save-report (persists it to a file),
--save-gnuplot (per-parameter curve data + a companion .gplot script,
energy + volume as a 2-panel multiplot), --view (an on-demand interactive
matplotlib preview of every discovered parameter, never saved to disk),
and -o/--output-dir.
"""

VERSION = "2.0.0"

import os
import re
import sys
import argparse
import numpy as np
from datetime import datetime
from stb.core import structure_io, siesta_log
from stb.core.cli import color_text, show_intro, print_dual, print_section, print_table
from stb.convergence import PARAMETER_DEFAULTS, PARAMETER_ORDER, substitute_numeric_tag, substitute_kgrid_tag

_FOLDER_RE = re.compile(r'^convergence_([a-zA-Z]+)_([\d.]+)$')
_RELAX_STEPS_RE = re.compile(r'MD\.(?:NumCGsteps|Steps)\s+(\d+)', re.IGNORECASE)
REPORT_FILE = "convergence_report.txt"

# Whether accuracy improves as the swept value increases -- derived from
# stb-convergence's own PARAMETER_DEFAULTS['tighter'] rather than a second,
# independently-maintained copy of the same physical fact (true for
# meshcutoff: a finer real-space grid; false for energyshift/kgrid: both
# converge towards a SMALLER value, the opposite direction from
# meshcutoff). Getting this wrong would make the plateau scan below walk
# away from convergence instead of towards it for 2 of the 3 parameters.
INCREASING_IMPROVES_ACCURACY = {p: (d["tighter"] == "larger") for p, d in PARAMETER_DEFAULTS.items()}


def parse_convergence_folder_name(folder_name):
    """Parses 'convergence_meshcutoff_250.0000' into ('meshcutoff', 250.0),
    matching the naming convention stb-convergence writes.
    """
    match = _FOLDER_RE.match(folder_name)
    if not match:
        return None, None
    return match.group(1), float(match.group(2))


def discover_parameter_sweeps(base_dir):
    """Finds every 'convergence_<parameter>_<value>' run folder reachable
    from base_dir, supporting both stb-convergence's default nested layout
    (base_dir/<parameter>/convergence_<parameter>_<value>/, one subfolder
    per swept parameter) and a flat layout (convergence_<parameter>_<value>/
    directly under base_dir -- pointing --dir at a single parameter's own
    subfolder, e.g. 'convergence_runs/meshcutoff', or an older/hand-built
    sweep). Returns {parameter: [folder_path, ...]}, auto-detecting every
    parameter present instead of requiring exactly one. Same flat-or-nested
    duality as siesta_log.find_strain_folders.
    """
    found = {}
    for entry in sorted(os.listdir(base_dir)):
        sub = os.path.join(base_dir, entry)
        if not os.path.isdir(sub):
            continue
        if entry in PARAMETER_ORDER:
            for leaf in sorted(os.listdir(sub)):
                leaf_path = os.path.join(sub, leaf)
                if not os.path.isdir(leaf_path):
                    continue
                param, _value = parse_convergence_folder_name(leaf)
                if param == entry:
                    found.setdefault(param, []).append(leaf_path)
        else:
            param, _value = parse_convergence_folder_name(entry)
            if param is not None:
                found.setdefault(param, []).append(sub)
    return found


def read_declared_relax_steps(folder):
    """Reads MD.NumCGsteps/MD.Steps back out of this folder's own
    config_extra.fdf (stb-convergence's forced-relaxation directive) --
    None if the file is absent (e.g. a hand-built folder), not an error.
    """
    path = os.path.join(folder, "config_extra.fdf")
    if not os.path.isfile(path):
        return None
    try:
        with open(path) as f:
            text = f.read()
    except OSError:
        return None
    m = _RELAX_STEPS_RE.search(text)
    return int(m.group(1)) if m else None


def quality_flag(scf_converged, max_force, force_tolerance, geometry_steps, declared_steps):
    """Short data-quality code for one run's table row -- 'OK', or a comma
    -joined list of concerns (never blocks the convergence scan itself,
    just flags a run whose own geometry/energy might not be trustworthy).
    """
    issues = []
    if not scf_converged:
        issues.append("SCF?")
    if max_force is not None and max_force > force_tolerance:
        issues.append("F>tol")
    if declared_steps is not None and geometry_steps >= declared_steps:
        issues.append("steps")
    return ",".join(issues) if issues else "OK"


def _scan_plateau(triples, tolerance, increasing_is_accurate, relative):
    """Shared plateau scan: `triples` is [(value, quantity, scf_converged), ...]
    sorted ascending by value. Scans from cheap towards accurate (reversed if
    increasing_is_accurate is False) and returns the first value beyond which
    every remaining step's delta stays under `tolerance` AND every remaining
    run has a confirmed SCF cycle -- a single small delta can be a
    coincidence (e.g. k-grid parity, or two SCF-unconverged points landing
    close by chance), not a genuine plateau, so the whole tail must hold.
    `relative=True` expresses each step's delta as a percent of the previous
    quantity (used for cell volume); `relative=False` uses the raw absolute
    delta (used for energy/atom). Returns None if no such point exists.
    """
    scan_source = triples if increasing_is_accurate else list(reversed(triples))
    deltas = []
    prev = None
    for value, quantity, scf_converged in scan_source:
        if prev is None:
            delta = None
        elif relative:
            delta = abs(quantity - prev) / prev * 100.0 if prev else None
        else:
            delta = abs(quantity - prev)
        deltas.append((value, delta, scf_converged))
        prev = quantity

    for i in range(1, len(deltas)):
        tail = deltas[i:]
        if all(d is not None and d < tolerance and scf_ok for _, d, scf_ok in tail):
            return deltas[i][0]
    return None


def scan_energy_convergence(rows, tolerance, increasing_is_accurate):
    """Energy-per-atom convergence: rows are the full per-value tuples from
    analyze_one_sweep (value at index 0, e_per_atom at index 2, scf_converged
    at index 6); tolerance is an absolute eV/atom delta.
    """
    triples = [(r[0], r[2], r[6]) for r in rows]
    return _scan_plateau(triples, tolerance, increasing_is_accurate, relative=False)


def scan_volume_convergence(rows, tolerance_pct, increasing_is_accurate):
    """Relaxed-cell-volume convergence: rows are the full per-value tuples
    from analyze_one_sweep (value at index 0, volume at index 4, scf_converged
    at index 6), pre-filtered by the caller to only rows with a real volume;
    tolerance_pct is a relative (%) delta between consecutive values.
    """
    triples = [(r[0], r[4], r[6]) for r in rows]
    return _scan_plateau(triples, tolerance_pct, increasing_is_accurate, relative=True)


def write_curve_plot(dat_path, parameter_name, rows):
    """Writes <dat_path> plus a companion .gplot next to it -- a 2-panel
    multiplot (E/atom and relaxed cell volume vs. the swept value), same
    <name>.dat + <name>.gplot convention as the rest of this workflow
    category (elastic_analysis.py/strain_analysis.py/eos_analysis.py).
    """
    with open(dat_path, 'w') as f:
        f.write(f"# Convergence curve | Parameter: {parameter_name}\n")
        f.write("# 1:Value 2:Energy(eV) 3:Energy_per_atom(eV) 4:Abs_Delta_E_per_atom(eV) "
                "5:Volume(Ang^3) 6:Abs_Delta_V_over_V(%) 7:SCF_converged\n")
        for value, energy, e_per_atom, de, volume, dv_pct, scf_converged, _mf, _q in rows:
            f.write(f"{value:.6f}  {energy:.6f}  {e_per_atom:.6f}  {(de or 0.0):.6f}  "
                    f"{(volume if volume is not None else 0.0):.6f}  {(dv_pct or 0.0):.6f}  "
                    f"{int(scf_converged)}\n")

    base = os.path.splitext(dat_path)[0]
    base_name = os.path.basename(base)
    gplot_path = f"{base}.gplot"
    dat_name = os.path.basename(dat_path)
    with open(gplot_path, 'w') as f:
        f.writelines([
            '# --- STB Plot Configuration ---\n',
            '# Generated by stb-convergenceAnalysis\n',
            'set terminal pdfcairo enhanced color font "Arial,14" size 7,8\n',
            f'set output "{base_name}.pdf"\n\n',
            'set multiplot layout 2,1\n\n',
            'set grid\n',
            f'set title "Convergence: {parameter_name} (energy)"\n',
            f'set xlabel "{parameter_name}"\n',
            'set ylabel "Energy per atom (eV)"\n',
            (f'plot "{dat_name}" using 1:3 with linespoints lw 2 pt 7 '
             f'lc rgb "#2255cc" title "E/atom"\n\n'),
            f'set title "Convergence: {parameter_name} (relaxed cell volume)"\n',
            f'set xlabel "{parameter_name}"\n',
            'set ylabel "Volume (Ang^3)"\n',
            (f'plot "{dat_name}" using 1:5 with linespoints lw 2 pt 7 '
             f'lc rgb "#cc5522" title "Volume"\n\n'),
            'unset multiplot\n',
        ])
    return gplot_path


def analyze_one_sweep(parameter, folder_paths, out_filename, force_tolerance, f_out):
    """Reads every folder's calc.out + structure.fdf + config_extra.fdf and
    returns (rows, n_skipped) sorted ascending by swept value. Each row is
    (value, energy, n_atoms, e_per_atom, delta_e_per_atom, volume,
    delta_volume_pct, scf_converged, max_force, quality).
    """
    data = []
    n_skipped = 0
    for folder in folder_paths:
        name = os.path.basename(folder)
        _param, value = parse_convergence_folder_name(name)
        out_path = os.path.join(folder, out_filename)
        struct_path = os.path.join(folder, "structure.fdf")

        if value is None or not os.path.isfile(out_path) or not os.path.isfile(struct_path):
            n_skipped += 1
            continue

        energy = siesta_log.get_free_energy(out_path)
        if energy is None:
            n_skipped += 1
            continue

        try:
            structure = structure_io.read_fdf(struct_path)
        except (ValueError, FileNotFoundError):
            n_skipped += 1
            continue
        n_atoms = len(structure.atoms)

        scf_converged, _iterations = siesta_log.get_scf_convergence(out_path)
        max_force = siesta_log.get_max_force(out_path)
        outcell = siesta_log.get_outcell(out_path)
        volume = abs(float(np.linalg.det(outcell))) if outcell is not None else None
        geometry_steps = siesta_log.count_geometry_steps(out_path)
        declared_steps = read_declared_relax_steps(folder)
        quality = quality_flag(scf_converged, max_force, force_tolerance, geometry_steps, declared_steps)

        data.append([value, energy, n_atoms, scf_converged, max_force, volume, quality])

    data.sort(key=lambda d: d[0])

    rows = []
    prev_e_per_atom = None
    prev_volume = None
    for value, energy, n_atoms, scf_converged, max_force, volume, quality in data:
        e_per_atom = energy / n_atoms
        de = None if prev_e_per_atom is None else abs(e_per_atom - prev_e_per_atom)
        dv_pct = None
        if volume is not None and prev_volume not in (None, 0):
            dv_pct = abs(volume - prev_volume) / prev_volume * 100.0
        rows.append((value, energy, e_per_atom, de, volume, dv_pct, scf_converged, max_force, quality))
        prev_e_per_atom = e_per_atom
        if volume is not None:
            prev_volume = volume

    return rows, n_skipped


def print_sweep_section(section_label, parameter, rows, n_skipped, n_found, args, f_out):
    """Prints one '[2.N] <PARAMETER> SWEEP' section: table + dual
    energy/volume convergence verdict + cross-check note. Returns
    (energy_converged, volume_converged, recommended) for the closing
    summary table.
    """
    print_section(section_label, f_out)
    increasing_is_accurate = INCREASING_IMPROVES_ACCURACY.get(parameter, True)
    print_dual(f"Folders found: {n_found} (skipped: {n_skipped} unreadable/incomplete)  |  "
               f"Increasing {parameter} improves accuracy: {increasing_is_accurate}", f_out)

    if not rows:
        print_dual(color_text("[WARNING] No valid data found for this parameter.", 'yellow'), f_out)
        return None, None, None

    header = ["Value", "E/atom(eV)", "|Delta E/atom|", "Volume(Ang^3)", "|Delta V/V|(%)", "Quality"]
    table_rows = []
    n_flagged = 0
    for value, _energy, e_per_atom, de, volume, dv_pct, _scf, _mf, quality in rows:
        if quality != "OK":
            n_flagged += 1
        table_rows.append(([
            f"{value:.4f}",
            f"{e_per_atom:.6f}",
            f"{de:.6f}" if de is not None else "--",
            f"{volume:.4f}" if volume is not None else "--",
            f"{dv_pct:.4f}" if dv_pct is not None else "--",
            quality,
        ], 'yellow' if quality != "OK" else None))
    print_table(header, table_rows, f_out)
    if n_flagged:
        print_dual(color_text(
            f"[WARNING] {n_flagged}/{len(rows)} run(s) flagged (SCF? = SCF convergence not "
            f"confirmed, F>tol = residual force above --force-tolerance "
            f"{args.force_tolerance} eV/Ang, steps = relaxation used its full declared step "
            "budget) -- their geometry may not be a genuine relaxed minimum.", 'yellow'), f_out)

    energy_converged = scan_energy_convergence(rows, args.tolerance, increasing_is_accurate)
    volume_rows = [r for r in rows if r[4] is not None]
    volume_converged = (scan_volume_convergence(volume_rows, args.volume_tolerance, increasing_is_accurate)
                         if len(volume_rows) >= 2 else None)

    direction_hint = "or higher" if increasing_is_accurate else "or lower"
    if energy_converged is not None:
        print_dual(f"{color_text('Energy converged at:', 'green')} {parameter} = "
                   f"{energy_converged:.4f} ({direction_hint}, |Delta E/atom| < {args.tolerance} eV)", f_out)
    else:
        print_dual(color_text(
            f"[WARNING] No point met the {args.tolerance} eV/atom energy tolerance within the "
            "scanned range; consider extending it.", 'yellow'), f_out)

    if volume_converged is not None:
        print_dual(f"{color_text('Structure converged at:', 'green')} {parameter} = "
                   f"{volume_converged:.4f} ({direction_hint}, |Delta V/V| < {args.volume_tolerance}%, "
                   "relaxed cell)", f_out)
    else:
        print_dual(color_text(
            f"[WARNING] No point met the {args.volume_tolerance}% relaxed-volume tolerance within "
            "the scanned range; consider extending it.", 'yellow'), f_out)

    recommended = None
    if energy_converged is not None and volume_converged is not None:
        if energy_converged == volume_converged:
            recommended = energy_converged
            print_dual("Energy and structure agree on the same converged value.", f_out)
        else:
            recommended = (max(energy_converged, volume_converged) if increasing_is_accurate
                           else min(energy_converged, volume_converged))
            print_dual(color_text(
                f"[NOTE] Energy and relaxed structure disagree on where {parameter} converges -- "
                f"using the more accurate-ward of the two ({recommended:.4f}) so BOTH criteria "
                "are satisfied. This is exactly the scenario stb-convergence's own docstring "
                "warns about: the total energy can look converged before the relaxed geometry "
                "actually is.", 'yellow'), f_out)
    elif energy_converged is not None:
        recommended = energy_converged
    elif volume_converged is not None:
        recommended = volume_converged

    return energy_converged, volume_converged, recommended


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Stage 2 of 2: analyzes a stb-convergence sweep and reports where "
                                     "each swept parameter's energy AND relaxed structure converge.", 'bold')}""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s --dir convergence_runs\n"
               "  %(prog)s --dir convergence_runs --tolerance 0.005 --save-report --save-gnuplot --view\n"
               "  %(prog)s --dir convergence_runs --apply calc.fdf\n"
    )

    parser.add_argument("--dir", type=str, default="convergence_runs",
                        help="Directory to scan for 'convergence_*' run folders -- either "
                             "stb-convergence's default nested layout (one <parameter>/ subfolder "
                             "per swept parameter under here) or a flat layout (pointing directly "
                             "at one parameter's own subfolder) (default: convergence_runs).")
    parser.add_argument("--file", type=str, default="calc.out",
                        help="SIESTA output filename inside each run folder (default: calc.out).")
    parser.add_argument("--tolerance", type=float, default=0.01,
                        help="Energy convergence tolerance in eV/atom (default: 0.01).")
    parser.add_argument("--volume-tolerance", type=float, default=0.5,
                        help="Relaxed-cell-volume convergence tolerance, as a percent change "
                             "between consecutive swept values (default: 0.5).")
    parser.add_argument("--force-tolerance", type=float, default=0.05,
                        help="Residual-force data-quality flag threshold in eV/Ang (default: "
                             "0.05) -- NOT a convergence criterion, just flags a run whose CG "
                             "relaxation may not have reached a genuine minimum.")
    parser.add_argument("-o", "--output-dir", type=str, default=".",
                        help="Directory to write <parameter>_convergence.dat/.gplot into with "
                             "--save-gnuplot, and the report into with --save-report (default: "
                             "current directory). Created if it doesn't exist.")
    parser.add_argument("--apply", type=str, default=None, metavar="CALC_FDF",
                        help="Writes every parameter's converged (recommended) value into this "
                             "calc.fdf in place (same tag substitution stb-convergence itself "
                             "uses) -- parameters with no converged value found are skipped.")
    parser.add_argument("--vacuum-gap", type=float, default=10.0,
                        help="--apply with a kgrid sweep only: vacuum-axis detection threshold in "
                             "Angstrom, same meaning as stb-kgrid's --vacuum-gap (default: 10.0).")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to <output-dir>/{REPORT_FILE}. Off by default.")
    parser.add_argument("--save-gnuplot", action="store_true",
                        help="Also write, per discovered parameter, curve data + a companion "
                             ".gplot script (2-panel: energy + relaxed volume). Off by default.")
    parser.add_argument("--view", action="store_true",
                        help="Show an interactive matplotlib preview (energy + volume per "
                             "discovered parameter) before finishing. Off by default.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-convergenceAnalysis {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    os.makedirs(args.output_dir, exist_ok=True)
    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    def fail(message):
        print_dual(color_text(f"[FAIL] {message}", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_dual(color_text(
        "===== STB-CONVERGENCEANALYSIS STAGE 2 REPORT (CONVERGENCE TEST ANALYSIS) =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time         : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", f_out)
    print_dual(f"Directory         : {args.dir}", f_out)
    print_dual(f"Output filename   : {args.file}", f_out)
    print_dual(f"Energy tolerance  : {args.tolerance} eV/atom", f_out)
    print_dual(f"Volume tolerance  : {args.volume_tolerance} %", f_out)
    print_dual(f"Force tolerance   : {args.force_tolerance} eV/Ang (data-quality flag only)", f_out)
    print_dual(f"Output dir        : {args.output_dir}", f_out)
    print_dual(f"Save report       : {'yes' if args.save_report else 'no'}", f_out)
    print_dual(f"Save gnuplot      : {'yes' if args.save_gnuplot else 'no'}", f_out)
    print_dual(f"View (matplotlib) : {'yes' if args.view else 'no'}", f_out)
    if report_path:
        print_dual(f"Report file       : {report_path}", f_out)

    if not os.path.isdir(args.dir):
        fail(f"Directory '{args.dir}' not found.")

    print_section("[1] DISCOVERED PARAMETER SWEEPS", f_out)
    sweeps = discover_parameter_sweeps(args.dir)
    if not sweeps:
        fail(f"No 'convergence_*' folders found under '{args.dir}' (expected either "
             f"'{args.dir}/<parameter>/convergence_<parameter>_<value>/', stb-convergence's "
             f"default layout, or 'convergence_<parameter>_<value>/' directly under '{args.dir}').")

    print_table(["Parameter", "Folders found"], [
        ([p, str(len(sweeps[p]))], None) for p in PARAMETER_ORDER if p in sweeps
    ], f_out)

    results = {}
    for i, parameter in enumerate([p for p in PARAMETER_ORDER if p in sweeps], start=1):
        rows, n_skipped = analyze_one_sweep(parameter, sweeps[parameter], args.file,
                                             args.force_tolerance, f_out)
        energy_c, volume_c, recommended = print_sweep_section(
            f"[2.{i}] {parameter.upper()} SWEEP", parameter, rows, n_skipped,
            len(sweeps[parameter]), args, f_out)
        results[parameter] = {
            "rows": rows, "energy_converged": energy_c, "volume_converged": volume_c,
            "recommended": recommended,
        }

    print_section("[3] OUTPUT FILES", f_out)
    if args.save_gnuplot:
        written = []
        for parameter, r in results.items():
            if not r["rows"]:
                continue
            dat_path = os.path.join(args.output_dir, f"{parameter}_convergence.dat")
            gplot_path = write_curve_plot(dat_path, parameter, r["rows"])
            written += [dat_path, gplot_path]
        if written:
            print_dual(color_text(f"[OK] Curve data + gnuplot script written: {len(written)} "
                                   "file(s).", 'green'), f_out)
            for path in written:
                print_dual(f"  - {path}", f_out)
            print_dual(f"  Render with: cd {args.output_dir} && gnuplot <parameter>_convergence.gplot", f_out)
        else:
            print_dual("No data available to plot.", f_out)
    else:
        print_dual("Not written (off by default -- pass --save-gnuplot to write, per "
                   "parameter, curve data + a companion .gplot script).", f_out)

    print_section("[4] SUMMARY & NEXT STEPS", f_out)
    print_table(["Parameter", "Energy-converged", "Structure-converged", "Recommended"], [
        ([p,
          f"{r['energy_converged']:.4f}" if r['energy_converged'] is not None else "--",
          f"{r['volume_converged']:.4f}" if r['volume_converged'] is not None else "--",
          f"{r['recommended']:.4f}" if r['recommended'] is not None else "NOT CONVERGED"],
         'yellow' if r['recommended'] is None else None)
        for p, r in results.items()
    ], f_out)
    if report_path:
        print_dual(f"Report            : {report_path}", f_out)

    if args.apply:
        print_section("[5] APPLY", f_out)
        try:
            with open(args.apply, 'r') as f:
                calc_text = f.read()
        except OSError as e:
            print_dual(color_text(f"[ERROR] Could not read '{args.apply}': {e}", 'red'), f_out)
            calc_text = None

        if calc_text is not None:
            applied = []
            for parameter, r in results.items():
                value = r["recommended"]
                if value is None:
                    print_dual(color_text(
                        f"  Skipped {parameter}: no converged value to apply.", 'yellow'), f_out)
                    continue
                try:
                    if parameter == "kgrid":
                        first_row_folder = sweeps[parameter][0]
                        struct_path = os.path.join(first_row_folder, "structure.fdf")
                        structure = structure_io.read_fdf(struct_path)
                        calc_text, divisions = substitute_kgrid_tag(calc_text, structure, value, args.vacuum_gap)
                        detail = f"grid {divisions[0]} {divisions[1]} {divisions[2]}"
                    else:
                        calc_text = substitute_numeric_tag(calc_text, parameter, value)
                        detail = f"{value:.4f} Ry"
                except (ValueError, FileNotFoundError) as e:
                    print_dual(color_text(f"  [ERROR] {parameter}: {e}", 'red'), f_out)
                else:
                    applied.append((parameter, value, detail))

            if applied:
                try:
                    with open(args.apply, 'w') as f:
                        f.write(calc_text)
                except OSError as e:
                    print_dual(color_text(f"[ERROR] Could not write '{args.apply}': {e}", 'red'), f_out)
                else:
                    for parameter, value, detail in applied:
                        print_dual(f"{color_text('[Applied]', 'green')} {parameter} = "
                                   f"{value:.4f} ({detail}) -> {args.apply}", f_out)

    if f_out:
        f_out.close()

    if args.view:
        import matplotlib.pyplot as plt
        plottable = {p: r["rows"] for p, r in results.items() if r["rows"]}
        if plottable:
            n = len(plottable)
            fig, axes = plt.subplots(n, 2, figsize=(11, 4 * n), squeeze=False)
            for i, (parameter, rows) in enumerate(plottable.items()):
                values = [r[0] for r in rows]
                e_per_atom = [r[2] for r in rows]
                volumes = [r[4] for r in rows]
                axes[i][0].plot(values, e_per_atom, 'o-', color='#2255cc')
                axes[i][0].set_title(f"{parameter}: Energy per atom")
                axes[i][0].set_xlabel(parameter)
                axes[i][0].set_ylabel("E/atom (eV)")
                axes[i][0].grid(True, alpha=0.3)
                if any(v is not None for v in volumes):
                    axes[i][1].plot(values, volumes, 'o-', color='#cc5522')
                axes[i][1].set_title(f"{parameter}: Relaxed cell volume")
                axes[i][1].set_xlabel(parameter)
                axes[i][1].set_ylabel("Volume (Ang^3)")
                axes[i][1].grid(True, alpha=0.3)
            fig.tight_layout()
            plt.show()


if __name__ == "__main__":
    main()
