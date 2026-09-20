#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

from stb import __version__ as VERSION

import os
import sys
import json
import argparse
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from stb.core.cli import color_text, show_intro, print_dual, print_section, print_table, get_input
from stb.core.siesta_log import get_free_energy, report_quality_diagnostics
from stb.core.phonon_workflow import load_phonon_with_force_constants

REPORT_FILE = "her_stage3.txt"

# Literature gas-phase H2 reference constants (Norskov et al., widely used
# throughout the CHE/HER literature) -- used for the H2 SIDE of every
# thermal-correction mode (even 'local'/'full', which compute the
# ADSORBED H*'s own ZPE/entropy from DFT+phonons): H2's own vibrational
# structure (one stretch mode + well-characterized rotational structure)
# is experimentally very well known, and re-deriving it from a DFT
# harmonic calculation would introduce more xc-functional noise than it
# removes. Fixed at their 298.15 K values regardless of --temp (a
# standard simplification in the field -- the ADSORBATE's own thermal
# terms are what's actually temperature-resolved here).
H2_ZPE_EV = 0.270
H2_TS_298K_EV = 0.400
NORSKOV_STANDARD_OFFSET_EV = 0.240

_FREQ_CONVERSION_THZ = 15.633302  # sqrt(eV / (amu * Ang^2)) -> THz (same constant ASE's own
                                  # vibrational-analysis code uses)
_BOLTZMANN_EV_K = 8.617333262e-5  # eV/K
_H_MASS_AMU = 1.00794
_EV_PER_THZ = 0.00413566733  # E = h*f, h in eV.s, f in THz*1e12 -> eV
_KJ_MOL_TO_EV = 0.01036427
_J_MOL_TO_EV = _KJ_MOL_TO_EV / 1000.0


def find_winning_site_energy(sites_root, out_file):
    """Re-scans 'sites/site_*/' for the lowest-FreeEng site -- same
    logic stb-herRefs already used to pick it (duplicated, not imported:
    each stage in this suite re-reads persisted files/folders rather than
    importing a sibling stage's module, so every stage stays independently
    invocable). Returns (label, energy) or (None, None) if nothing usable
    is found.
    """
    if not os.path.isdir(sites_root):
        return None, None
    site_dirs = sorted(
        d.path for d in os.scandir(sites_root) if d.is_dir() and d.name.startswith("site_")
    )
    best_label, best_energy = None, float('inf')
    for d in site_dirs:
        energy = get_free_energy(os.path.join(d, out_file))
        if energy is not None and energy < best_energy:
            best_energy = energy
            best_label = os.path.basename(d)
    if best_label is None:
        return None, None
    return best_label, best_energy


def read_fa_force(fa_path, atom_index):
    """Reads the force (Fx, Fy, Fz, eV/Ang) on ONE specific atom (0-based
    `atom_index`) from a SIESTA .FA file -- format: first line = atom
    count, then one row per atom '<1-based index> Fx Fy Fz'. Same
    fail-soft "None on any parse issue" contract as core.siesta_log's own
    parsers (this one isn't in core/ since HER's local-ZPE mode is its
    only consumer so far -- extract on second use, same policy as
    everywhere else in this suite).
    """
    try:
        with open(fa_path) as f:
            lines = f.readlines()
    except OSError:
        return None
    if len(lines) < atom_index + 2:
        return None
    try:
        parts = lines[atom_index + 1].split()
        return np.array([float(parts[1]), float(parts[2]), float(parts[3])])
    except (IndexError, ValueError):
        return None


def compute_local_zpe_entropy(zpe_dir, temperature_k, f_out):
    """Reads the 6 disp_NNN/her_zpe.FA force files stb-herRefs
    --zpe-mode local wrote, builds the 3x3 partial Hessian for the H
    atom via central finite difference (substrate held fixed -- a
    decoupled-oscillator approximation, see her_refs.py's own
    write_local_zpe_folders docstring), diagonalizes it, and returns
    (zpe_ev, ts_ev) for the ADSORBED H* -- standard harmonic-oscillator
    formulas (same math already used by Raman's Bose-Einstein weighting
    this suite already has, core/spectrum.py::bose_einstein_weight, just
    applied here to get absolute ZPE/entropy rather than a relative
    Stokes population factor). Returns (None, None) on any read failure.
    """
    meta_path = os.path.join(zpe_dir, "zpe_local_meta.json")
    try:
        with open(meta_path) as f:
            meta = json.load(f)
    except (OSError, ValueError):
        print_dual(color_text(f"    [ERROR] Could not read '{meta_path}'.", 'red'), f_out)
        return None, None

    h_index = meta["h_index"]
    displacement_ang = meta["displacement_ang"]
    order = meta["order"]

    forces = {}
    for i, entry in enumerate(order, start=1):
        fa_path = os.path.join(zpe_dir, f"disp_{i:03d}", "her_zpe.FA")
        force = read_fa_force(fa_path, h_index)
        if force is None:
            print_dual(color_text(f"    [ERROR] Could not read forces from '{fa_path}'.", 'red'), f_out)
            return None, None
        forces[(entry["axis"], entry["sign"])] = force

    Phi = np.zeros((3, 3))
    for axis in range(3):
        f_plus = forces[(axis, 1.0)]
        f_minus = forces[(axis, -1.0)]
        Phi[:, axis] = -(f_plus - f_minus) / (2.0 * displacement_ang)
    Phi = 0.5 * (Phi + Phi.T)

    eigenvalues = np.linalg.eigvalsh(Phi)
    zpe_ev, ts_ev = 0.0, 0.0
    print_dual("    Local (partial-Hessian) H* vibrational modes:", f_out)
    for eig in eigenvalues:
        d_eig = eig / _H_MASS_AMU
        if d_eig > 0:
            freq_thz = _FREQ_CONVERSION_THZ * np.sqrt(d_eig)
            e_mode = freq_thz * _EV_PER_THZ
            print_dual(f"      {freq_thz * 33.35641:>8.2f} cm^-1  ({freq_thz:>6.2f} THz)  "
                        f"E = {e_mode:.4f} eV", f_out)
            zpe_ev += 0.5 * e_mode
            if temperature_k > 0:
                x = e_mode / (_BOLTZMANN_EV_K * temperature_k)
                S = _BOLTZMANN_EV_K * (x / np.expm1(x) - np.log1p(-np.exp(-x)))
                ts_ev += temperature_k * S
        else:
            freq_thz = _FREQ_CONVERSION_THZ * np.sqrt(-d_eig)
            print_dual(color_text(
                f"      {freq_thz:>6.2f} THz (IMAGINARY) -- excluded from ZPE/entropy", 'yellow'), f_out)
    return zpe_ev, ts_ev


def compute_full_zpe_entropy(zpe_dir, system_label, temperature_k, f_out):
    """Full-structure ZPE/entropy via phonopy's own thermal-properties
    pipeline, reusing core.phonon_workflow.load_phonon_with_force_constants
    (the SAME FORCE_SETS-building/loading code Raman/IR/Phonons already
    use) -- zero new phonon code for this mode, only the calling
    convention differs (Gamma-only mesh, since this "supercell" is
    already the full slab, no further periodic replication needed).
    Returns (zpe_ev, ts_ev) or (None, None) on failure.
    """
    try:
        phonon, _internal_to_angstrom, original_dir = load_phonon_with_force_constants(
            zpe_dir, system_label, False, f_out)
    except SystemExit:
        return None, None
    try:
        phonon.run_mesh([1, 1, 1])
        phonon.run_thermal_properties(t_step=10, t_max=temperature_k + 10, t_min=temperature_k)
        tp_dict = phonon.get_thermal_properties_dict()
        temps = np.array(tp_dict['temperatures'])
        idx = int(np.abs(temps - temperature_k).argmin())
        zpe_ev = phonon.get_zero_point_energy() * _KJ_MOL_TO_EV
        entropy_ev_k = tp_dict['entropy'][idx] * _J_MOL_TO_EV
        ts_ev = temperature_k * entropy_ev_k
    finally:
        os.chdir(original_dir)
    return zpe_ev, ts_ev


# Colors for the [4] FINAL RESULT breakdown table AND the gnuplot bar chart
# (write_energy_breakdown_plot below) -- kept in sync so the printed table
# and the saved plot tell exactly the same story: a plain additive
# contribution ('term'), a running checkpoint ('subtotal', Delta-E_H
# corrected -- raw + BSSE), and the final answer ('total', Delta-G_H*
# itself). print_table takes a color NAME (core.cli.COLORS); the plot needs
# an actual hex string per bar -- two small dicts, not one, since the two
# consumers want different value types for the same three kinds.
_BREAKDOWN_TABLE_COLOR = {"term": None, "subtotal": "cyan", "total": "green"}
# Plain hex (no '#', no quotes) -- written into the .dat file as a gnuplot
# packed-RGB integer literal (0xRRGGBB), NOT a quoted color-name string:
# `lc rgb variable` reads an INTEGER from the data column (see `gnuplot
# -e 'help linecolor'`: "rgbcolor variable # integer value is read from
# input file") -- a quoted "#RRGGBB" string in that column silently breaks
# the whole plot ("x range is invalid"), verified live against a real
# gnuplot install before settling on this format.
_BREAKDOWN_PLOT_COLOR = {"term": "2255cc", "subtotal": "888888", "total": "22aa55"}


def write_energy_breakdown_plot(plot_dir, breakdown_rows):
    """Writes a gnuplot .dat + .gplot categorical bar chart of the additive
    energy contributions that build up to Delta-G_H* -- same <name>.dat +
    <name>.gplot convention as the rest of the suite (e.g.
    adsorb_analysis.py's write_curve_plot, xrdrank.py's ranking bar chart).
    Bars are colored by `kind` (term/subtotal/total, see
    _BREAKDOWN_PLOT_COLOR) via gnuplot's `lc rgb variable`, so the running-
    subtotal and final-total bars visually stand out from the plain
    contribution terms -- the same three-way distinction the [4] FINAL
    RESULT table already prints (`breakdown_rows` is that exact list, so
    the plot and the report can never disagree). A y=0 reference line (set
    xzeroaxis) matters here specifically because some contributions (e.g.
    -Delta-TS) can legitimately be negative.

    Explicit x positions (column 1) rather than `set style data histograms`
    -- verified live that `using 2:xtic(3):4 with boxes lc rgb variable`
    (histogram style, implicit x) fails outright ("x range is invalid");
    `using 1:2:4:xtic(3)` (explicit x, then y, then the variable-color
    column, with xtic(label) last) is the form that actually renders.

    `breakdown_rows` is [(label, value_eV, kind), ...]. Returns
    (dat_path, gplot_path).
    """
    os.makedirs(plot_dir, exist_ok=True)
    dat_path = os.path.join(plot_dir, "her_energy_breakdown.dat")
    with open(dat_path, "w") as f:
        f.write("# HER Delta-G_H* energy breakdown\n")
        f.write("# 1:Index 2:Value(eV) 3:Label 4:Color(0xRRGGBB, term/subtotal/total)\n")
        for i, (label, value, kind) in enumerate(breakdown_rows, start=1):
            f.write(f'{i}  {value:.6f}  "{label}"  0x{_BREAKDOWN_PLOT_COLOR[kind]}\n')

    # gnuplot anchors rotated xtic labels right at the bottom of the plot
    # FRAME, not at y=0 -- with an auto-scaled yrange, a bar reaching close
    # to that bottom edge leaves the rotated (diagonal) label nowhere to go
    # but back up INTO the bar itself (verified live: the default auto-range
    # visibly ran the labels through the bar fills). Fix: force an explicit
    # yrange with generous headroom below the lowest bar (and a little above
    # the highest one) so the labels always land in clear space beneath the
    # zero line, scaled to the data's own span rather than a fixed eV
    # constant so it works whether the bars span ~2 eV (typical) or ~80 eV
    # (a pathological fixture/test case).
    values = [value for _label, value, _kind in breakdown_rows]
    y_min, y_max = min(values + [0.0]), max(values + [0.0])
    span = max(y_max - y_min, 1e-6)
    y_bottom = min(y_min, 0.0) - 0.35 * span
    y_top = max(y_max, 0.0) + 0.15 * span

    gplot_path = os.path.join(plot_dir, "her_energy_breakdown.gplot")
    with open(gplot_path, "w") as f:
        f.writelines([
            '# --- STB Plot Configuration ---\n',
            '# Generated by stb-herAnalysis\n',
            'set terminal pdfcairo enhanced color font "Arial,12" size 8,6\n',
            'set output "her_energy_breakdown.pdf"\n\n',
            'set title "HER Delta-G_H* energy breakdown"\n',
            'set ylabel "Energy (eV)"\n',
            'set style fill solid 0.7 border -1\n',
            'set boxwidth 0.6\n',
            'set xtics rotate by -30 right\n',
            f'set yrange [{y_bottom:.6f}:{y_top:.6f}]\n',
            'set grid ytics lt 0 lw 1 lc rgb "#bbbbbb"\n',
            'set xzeroaxis lt -1 lw 2 lc rgb "#000000"\n',
            'unset key\n\n',
            'plot "her_energy_breakdown.dat" using 1:2:4:xtic(3) with boxes lc rgb variable\n',
        ])
    return dat_path, gplot_path


def show_energy_breakdown_matplotlib(breakdown_rows):
    """Opens a blocking matplotlib window with the SAME energy-breakdown
    bar chart write_energy_breakdown_plot saves for gnuplot -- an on-screen
    convenience alongside (not instead of) the gnuplot .dat/.gplot pair,
    which stays the canonical saved-to-disk format for this WORKFLOW_TOOLS
    stage (CLAUDE.md's "Plot conventions differ by category" convention;
    MLSIM_TOOLS write matplotlib PNGs directly, but this isn't one of
    those -- same "--view-plots shows the same plot the tool already
    saves, doesn't replace the saved format" precedent as stb-adsorb's own
    --view-plots flag). Reuses `_BREAKDOWN_PLOT_COLOR` (prefixed with '#'
    for matplotlib, which -- unlike gnuplot's `lc rgb variable` -- takes a
    plain "#RRGGBB" string directly) so the on-screen chart and the saved
    PDF are colored identically. Returns True if the window was shown,
    False if matplotlib couldn't open a display (e.g. a headless SSH
    session with no X11/Wayland forwarding) -- caller reports that as an
    advisory note, not a crash.
    """
    labels = [label for label, _value, _kind in breakdown_rows]
    values = [value for _label, value, _kind in breakdown_rows]
    colors = [f"#{_BREAKDOWN_PLOT_COLOR[kind]}" for _label, _value, kind in breakdown_rows]
    try:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.bar(labels, values, color=colors, edgecolor='black', linewidth=0.5)
        ax.axhline(0, color='black', linewidth=1)
        ax.set_ylabel("Energy (eV)")
        ax.set_title("HER $\\Delta G_{H^*}$ energy breakdown")
        ax.grid(axis='y', linestyle=':', color='#bbbbbb')
        plt.setp(ax.get_xticklabels(), rotation=25, ha='right')
        fig.tight_layout()
        plt.show()
        return True
    except Exception as e:
        print(color_text(f"[NOTE] Could not open a matplotlib window ({e}) -- no display "
                          "available (e.g. a headless/SSH session without X11 forwarding). "
                          "Use --plot for the gnuplot .dat/.gplot files instead.", 'yellow'))
        return False


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Stage 3 of 3: combines every reference energy from stb-herRefs "
        "into the final Delta-G_H* (Gibbs free energy of H adsorption).", 'bold')}
Delta-G_H* = Delta-E_H (BSSE-corrected) + Delta-ZPE - T*Delta-S, the computational hydrogen
electrode (CHE, Norskov et al.) descriptor for HER activity -- Delta-G_H* ~= 0 is the Sabatier-
optimal binding strength.

--zpe-mode full's Delta-ZPE/Delta-S computation is NOT simply ZPE(site) - 0.5*ZPE(H2): that
would include the ENTIRE slab's own vibrational zero-point energy (perturbed by H adsorption),
not just H*'s contribution. The physically correct expression subtracts the clean slab's own
(equally full-structure) ZPE/entropy first:
    Delta-ZPE = [ZPE(site) - ZPE(clean_slab)] - 0.5*ZPE(H2)
    Delta-S   = [S(site) - S(clean_slab)] - 0.5*S(H2)
which is why stb-herRefs --zpe-mode full computes TWO full phonon calculations (site AND clean
slab), not one. --zpe-mode local sidesteps this entirely: it never displaces slab atoms in the
first place (a decoupled H*-only oscillator approximation), so there is no slab ZPE/entropy to
subtract.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage example:\n"
               "  %(prog)s --directory her_study --temp 298.15\n"
    )

    parser.add_argument("-dir", "--directory", type=str, default="her_study",
                         help="Root directory written by stb-her/stb-herRefs (default: her_study).")
    parser.add_argument("--file", type=str, default="calc.out",
                         help="SIESTA output filename inside each folder (default: calc.out).")
    parser.add_argument("--temp", type=float, default=298.15,
                         help="Temperature in Kelvin, for the local/full ZPE modes' entropy "
                              "term (default: 298.15). Has no effect on --zpe-mode standard "
                              "(fixed Norskov offset) or on the H2 gas-phase reference (fixed "
                              "at its own 298.15 K literature value regardless of --temp -- see "
                              "the module docstring's H2_ZPE_EV/H2_TS_298K_EV).")
    parser.add_argument("--force-tolerance", type=float, default=0.05,
                         help="Residual atomic force in eV/Ang (default: 0.05) above which a "
                              "folder is flagged as possibly not relaxed/converged. Advisory only.")
    parser.add_argument("-o", "--output", type=str, default="HER_report",
                         help="Base filename (no extension) for the final report (default: HER_report).")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to <directory>/{REPORT_FILE}. Off by default.")
    parser.add_argument("--plot", dest="plot", action="store_true", default=None,
                         help="Save a gnuplot .dat/.gplot bar chart of the Delta-G_H* energy "
                              "breakdown to <directory>/plot/, skipping the interactive prompt "
                              "this tool otherwise asks at the end. See --no-plot.")
    parser.add_argument("--no-plot", dest="plot", action="store_false",
                         help="Skip the energy-breakdown plot entirely, skipping the interactive "
                              "prompt (see --plot). Use this for non-interactive/scripted runs.")
    parser.add_argument("--show", dest="show", action="store_true", default=None,
                         help="Open a blocking matplotlib window with the same energy-breakdown "
                              "chart on screen, skipping the interactive prompt this tool "
                              "otherwise asks. Independent of --plot/--no-plot (the gnuplot "
                              "files, if any, are still the canonical saved output). Needs a "
                              "display. See --no-show.")
    parser.add_argument("--no-show", dest="show", action="store_false",
                         help="Skip the on-screen matplotlib window entirely, skipping the "
                              "interactive prompt (see --show). Use this for non-interactive/"
                              "scripted/headless runs.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-herAnalysis {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("HER WORKFLOW -- STAGE 3: ANALYSIS", 'bold'))
    print("-" * 60)

    output_root = args.directory
    stage2_report = os.path.join(output_root, "her_stage2.txt")
    if not os.path.isfile(stage2_report):
        print(color_text(f"[ERROR] '{stage2_report}' not found -- run stb-herRefs (Stage 2) first.", 'red'))
        sys.exit(1)

    zpe_mode = None
    with open(stage2_report) as f:
        for line in f:
            if line.startswith("ZPE mode       :"):
                zpe_mode = line.split(":", 1)[1].strip()
                break
    if zpe_mode is None:
        print(color_text(f"[ERROR] Could not recover the ZPE mode from '{stage2_report}'.", 'red'))
        sys.exit(1)

    energy_rows = []  # (label, energy_eV, source_path) -- the [1] recap table at the end

    def read_energy(folder, label, f_out, skip_quality_check=False):
        out_path = os.path.join(output_root, folder, args.file)
        energy = get_free_energy(out_path)
        if energy is None:
            print_dual(color_text(f"[ERROR] Could not read energy from '{out_path}'.", 'red'), f_out)
            sys.exit(1)
        print_dual(f"  {label:<14} = {energy:>12.6f} eV  ({out_path})", f_out)
        if not skip_quality_check:
            report_quality_diagnostics(label, out_path, args.force_tolerance, f_out)
        # 03_slab_deformed/04_slab_ghost (skip_quality_check=True below) are
        # single-point evaluations of the WINNING SITE's own relaxed geometry
        # with H removed/ghosted -- the slab atoms nearest the former H
        # position are never at their own equilibrium there BY CONSTRUCTION
        # (verified live: only the 2-3 atoms closest to the former H carry a
        # large force, ~1 eV/Ang each, while every other atom stays under
        # 0.03 eV/Ang and the NET force is ~0 -- a localized, physically
        # -explainable effect of removing/ghosting H, not a sign of a broken
        # or unconverged calculation). report_quality_diagnostics' force
        # check has no way to know that distinction, so it's skipped
        # entirely for these two folders rather than printing a warning (or
        # even a softened note) that would always fire and never mean
        # anything actionable here.
        energy_rows.append((label, energy, out_path))
        return energy

    report_path = os.path.join(output_root, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None
    print_dual(f"{color_text('===== HER STAGE 3 REPORT (ANALYSIS) =====', 'magenta')}", f_out)

    print_section('[0] RUN METADATA', f_out)
    print_dual(f"Date/time       : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", f_out)
    print_dual(f"Directory       : {output_root}", f_out)
    print_dual(f"SIESTA output   : {args.file}", f_out)
    print_dual(f"ZPE mode        : {zpe_mode}", f_out)
    print_dual(f"Temperature     : {args.temp} K", f_out)
    print_dual(f"Force tolerance : {args.force_tolerance} eV/Ang (advisory)", f_out)
    print_dual(f"Save report     : {'yes -> ' + os.path.join(output_root, REPORT_FILE) if args.save_report else 'no'}",
                f_out)

    print_section('[1] ELECTRONIC ENERGIES (FreeEng)', f_out)
    e_clean = read_energy("00_clean_slab", "E_clean", f_out)
    site_label, e_slab_h = find_winning_site_energy(os.path.join(output_root, "sites"), args.file)
    if site_label is None:
        print_dual(color_text("[ERROR] Could not find a winning site with a readable energy "
                               "under 'sites/'.", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)
    print_dual(f"  {'E_slab+H':<14} = {e_slab_h:>12.6f} eV  (winning site: {site_label})", f_out)
    report_quality_diagnostics(site_label, os.path.join(output_root, "sites", site_label, args.file),
                                args.force_tolerance, f_out)
    energy_rows.append(("E_slab+H", e_slab_h, os.path.join(output_root, "sites", site_label, args.file)))
    e_h2 = read_energy("02_h2_molecule", "E_H2", f_out)
    e_deformed = read_energy("03_slab_deformed", "E_deformed", f_out, skip_quality_check=True)
    e_ghost = read_energy("04_slab_ghost", "E_ghost", f_out, skip_quality_check=True)
    e_h_ghost_slab = read_energy("06_h_ghost_slab", "E_H_ghost", f_out)
    e_h_isolated = read_energy("07_h_isolated", "E_H_iso", f_out)

    print_dual("", f_out)
    print_table(["Term", "Energy (eV)", "Source"],
                [([label, f"{energy:.6f}", path], None) for label, energy, path in energy_rows],
                f_out)

    print_section('[2] BSSE CORRECTION', f_out)
    e_deformation = e_deformed - e_clean
    bsse_slab = e_deformed - e_ghost
    bsse_adsorbate = e_h_isolated - e_h_ghost_slab
    bsse_correction = bsse_slab + bsse_adsorbate
    delta_e_raw = e_slab_h - e_clean - 0.5 * e_h2
    delta_e_corrected = delta_e_raw + bsse_correction
    print_table(["Term", "Value (eV)"], [
        (["Deformation cost (diagnostic, not used below)", f"{e_deformation:+.4f}"], None),
        (["BSSE (slab side)", f"{bsse_slab:+.4f}"], None),
        (["BSSE (H side)", f"{bsse_adsorbate:+.4f}"], None),
        (["BSSE (total)", f"{bsse_correction:+.4f}"], 'cyan'),
        (["Delta-E_H (raw)", f"{delta_e_raw:+.4f}"], None),
        (["Delta-E_H (corrected)", f"{delta_e_corrected:+.4f}"], 'cyan'),
    ], f_out)

    print_section('[3] THERMAL CORRECTION', f_out)
    if zpe_mode == "standard":
        delta_thermo = NORSKOV_STANDARD_OFFSET_EV
        print_dual(f"  Standard Norskov offset: Delta-ZPE - T*Delta-S ~= "
                    f"{NORSKOV_STANDARD_OFFSET_EV:+.4f} eV", f_out)
    elif zpe_mode == "local":
        zpe_dir = os.path.join(output_root, "05_zpe_calc")
        zpe_h_star, ts_h_star = compute_local_zpe_entropy(zpe_dir, args.temp, f_out)
        if zpe_h_star is None:
            if f_out:
                f_out.close()
            sys.exit(1)
        delta_zpe = zpe_h_star - 0.5 * H2_ZPE_EV
        delta_ts = ts_h_star - 0.5 * H2_TS_298K_EV
        delta_thermo = delta_zpe - delta_ts
        print_dual(f"  ZPE(H*, local)  = {zpe_h_star:.4f} eV   1/2 ZPE(H2, lit.) = "
                    f"{0.5 * H2_ZPE_EV:.4f} eV   Delta-ZPE = {delta_zpe:+.4f} eV", f_out)
        print_dual(f"  TS(H*, local)   = {ts_h_star:.4f} eV   1/2 TS(H2, lit.)  = "
                    f"{0.5 * H2_TS_298K_EV:.4f} eV   Delta-TS  = {delta_ts:+.4f} eV", f_out)
        print_dual(f"  Net correction (Delta-ZPE - Delta-TS) : {delta_thermo:+.4f} eV", f_out)
    else:  # full
        print_dual("  Site (slab+H) full phonon calculation:", f_out)
        zpe_site, ts_site = compute_full_zpe_entropy(
            os.path.join(output_root, "05_zpe_calc_site"), "her_zpe_site", args.temp, f_out)
        print_dual("  Clean slab full phonon calculation:", f_out)
        zpe_clean, ts_clean = compute_full_zpe_entropy(
            os.path.join(output_root, "05_zpe_calc_clean"), "her_zpe_clean", args.temp, f_out)
        if zpe_site is None or zpe_clean is None:
            print_dual(color_text("[ERROR] Full-mode phonon calculation(s) incomplete.", 'red'), f_out)
            if f_out:
                f_out.close()
            sys.exit(1)
        delta_zpe = (zpe_site - zpe_clean) - 0.5 * H2_ZPE_EV
        delta_ts = (ts_site - ts_clean) - 0.5 * H2_TS_298K_EV
        delta_thermo = delta_zpe - delta_ts
        print_dual(f"  ZPE(site)={zpe_site:.4f}  ZPE(clean)={zpe_clean:.4f}  "
                    f"1/2 ZPE(H2, lit.)={0.5 * H2_ZPE_EV:.4f}  Delta-ZPE={delta_zpe:+.4f} eV", f_out)
        print_dual(f"  TS(site)={ts_site:.4f}  TS(clean)={ts_clean:.4f}  "
                    f"1/2 TS(H2, lit.)={0.5 * H2_TS_298K_EV:.4f}  Delta-TS={delta_ts:+.4f} eV", f_out)
        print_dual(f"  Net correction (Delta-ZPE - Delta-TS) : {delta_thermo:+.4f} eV", f_out)

    delta_g = delta_e_corrected + delta_thermo

    # Same list feeds the [4] table below AND write_energy_breakdown_plot,
    # so the printed report and the saved plot can never disagree. 'term' =
    # a plain additive contribution, 'subtotal' = a running checkpoint
    # (Delta-E_H corrected -- raw + BSSE, already shown in [2] too),
    # 'total' = the final Delta-G_H* itself.
    breakdown_rows = [
        ("Delta-E_H (raw)", delta_e_raw, "term"),
        ("BSSE correction", bsse_correction, "term"),
        ("Delta-E_H (corrected)", delta_e_corrected, "subtotal"),
    ]
    if zpe_mode == "standard":
        breakdown_rows.append(("Norskov thermal offset", NORSKOV_STANDARD_OFFSET_EV, "term"))
    else:
        breakdown_rows.append(("Delta-ZPE", delta_zpe, "term"))
        breakdown_rows.append(("-Delta-TS", -delta_ts, "term"))
    breakdown_rows.append(("Delta-G_H* (TOTAL)", delta_g, "total"))

    print_dual(f"\n{color_text('[4] FINAL RESULT', 'magenta')}", f_out)
    print_dual("=" * 60, f_out)
    print_dual(f"  Delta-G_H* = {delta_g:+.4f} eV", f_out)
    print_dual("=" * 60, f_out)
    verdict = ("near-optimal (Sabatier principle, |Delta-G_H*| small)" if abs(delta_g) < 0.2
               else "too strong (H* binds too tightly, poor H2 release)" if delta_g < 0
               else "too weak (H* adsorption itself is the bottleneck)")
    print_dual(f"  Qualitative HER assessment: {verdict}", f_out)
    print_dual("", f_out)
    print_dual("  Energy breakdown (each term's contribution to Delta-G_H*, top to bottom):",
                f_out)
    print_table(["Contribution", "Value (eV)"],
                [([label, f"{value:+.4f}"], _BREAKDOWN_TABLE_COLOR[kind])
                 for label, value, kind in breakdown_rows], f_out)

    report_str_path = os.path.join(output_root, f"{args.output}.txt")
    with open(report_str_path, "w") as f_report:
        f_report.write(f"HER Delta-G_H* = {delta_g:+.4f} eV\n")
        f_report.write(f"ZPE mode = {zpe_mode}, T = {args.temp} K\n")
        f_report.write(f"Delta-E_H (corrected) = {delta_e_corrected:+.4f} eV\n")
        f_report.write(f"Thermal correction = {delta_thermo:+.4f} eV\n")

    want_plot = args.plot
    if want_plot is None:
        try:
            answer = get_input(f"\nSave a gnuplot energy-breakdown chart to "
                                f"'{os.path.join(output_root, 'plot')}'? [y/N]: ").strip().lower()
        except EOFError:
            answer = ""
        want_plot = answer in ("y", "yes")
    if want_plot:
        plot_dir = os.path.join(output_root, "plot")
        dat_path, gplot_path = write_energy_breakdown_plot(plot_dir, breakdown_rows)
        print_dual(f"\n{color_text('[Saved]', 'cyan')} {dat_path}, {gplot_path} "
                    f"(cd {plot_dir} && gnuplot {os.path.basename(gplot_path)} -- writes "
                    "her_energy_breakdown.pdf)", f_out)
    else:
        dat_path = gplot_path = None
        print_dual("\nNo plot generated.", f_out)

    want_show = args.show
    if want_show is None:
        try:
            answer = get_input("Show the energy-breakdown chart on screen now (matplotlib)? "
                                "[y/N]: ").strip().lower()
        except EOFError:
            answer = ""
        want_show = answer in ("y", "yes")
    if want_show:
        shown = show_energy_breakdown_matplotlib(breakdown_rows)
        if shown:
            print_dual("Matplotlib window closed.", f_out)

    print_section('[5] SUMMARY & FILES', f_out)
    print_dual(f"Winning site        : {site_label}", f_out)
    print_dual(f"ZPE mode            : {zpe_mode}", f_out)
    print_dual(f"Delta-G_H*          : {delta_g:+.4f} eV", f_out)
    if report_path:
        print_dual(f"Report              : {report_path}", f_out)
    print_dual(f"Files               : {report_str_path}", f_out)
    if dat_path:
        print_dual(f"Plot data           : {dat_path}", f_out)
        print_dual(f"Plot script         : {gplot_path}", f_out)

    if f_out:
        f_out.close()

    print("\n[INFO] Complete job!")
    print("\n" + "-" * 60)
    print(color_text("HER analysis complete.\n", 'bold'))


if __name__ == "__main__":
    main()
