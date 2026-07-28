#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "2.0.0"

import os
import sys
import argparse
import matplotlib.pyplot as plt
from stb.core import structure_io
from stb.core import citations
from stb.core import xrd as xrd_core
from stb.core.symmetry import symmetry_summary
from stb.core.cli import color_text, show_intro, run_with_spinner, print_dual, print_section, print_table
from stb.core.deps import require_pyxtal

REPORT_FILE = "stb_xrd_report.txt"
BIB_FILE = "references.bib"


def plot_matplotlib(xrd, structure, wavelength_label, two_theta_range,
                     sim_profile=None, experimental=None, compare_to=None, similarity=None):
    """Shows an interactive preview: the simulated stick pattern (impulses,
    2-theta vs. normalized intensity), or -- with --compare-to -- an overlay
    of the simulated and experimental profiles instead.

    Deliberately NOT pyxtal's own XRD.plot_pxrd(): checked pyxtal's
    installed source directly and its show=True path only ever calls the
    non-blocking fig.show(), never the blocking plt.show() -- exactly the
    cause of a real, reported bug ("the plot disappears immediately"). Built
    here instead from data this tool already has (xrd.pxrd), ending in a
    blocking plt.show(), the same convention already used correctly by
    stb-workfunction/stb-density's own --view.
    """
    lo, hi = two_theta_range
    fig, ax = plt.subplots(figsize=(10, 5))

    if experimental is not None:
        sim_x, sim_y = sim_profile
        sim_y = sim_y / sim_y.max() if sim_y.max() > 0 else sim_y
        exp_x, exp_y = experimental
        exp_y = exp_y / exp_y.max() if exp_y.max() > 0 else exp_y

        ax.plot(sim_x, sim_y, label="Simulated", color="C0")
        ax.plot(exp_x, exp_y, label=f"Experimental ({compare_to})", color="C1", alpha=0.7)
        ax.set_ylabel("Normalized intensity")
        ax.set_title(f"{structure.composition.reduced_formula} -- Similarity: {similarity:.4f}")
        ax.legend()
    else:
        two_theta = [row[0] for row in xrd.pxrd]
        intensity = [row[5] for row in xrd.pxrd]
        max_i = max(intensity) if intensity else 1.0
        norm_i = [i / max_i * 100.0 for i in intensity]
        ax.vlines(two_theta, 0, norm_i, colors="C0", linewidth=1.5)
        ax.set_ylim(0, 105)
        ax.set_ylabel("Intensity (%)")
        ax.set_title(f"Simulated XRD -- {structure.composition.reduced_formula} ({wavelength_label})")

    ax.set_xlabel("2-theta (deg)")
    ax.set_xlim(lo, hi)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    plt.show()


def main():
    require_pyxtal()

    parser = argparse.ArgumentParser(
        description=f"""{color_text("Simulates a powder XRD pattern from a structure.", 'bold')}
Uses pyxtal's own diffraction engine (not pymatgen's) so the same dependency
already needed by stb-crystalcast covers this too.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s --file structure.fdf --format fdf\n"
               "  %(prog)s --file structure.fdf --format fdf --wavelength MoKa --top 10\n"
               "  %(prog)s --file relaxed.STRUCT_OUT --format struct_out \\\n"
               "      --two-theta-range 10 80 --view\n"
               "  %(prog)s --file structure.fdf --format fdf \\\n"
               "      --compare-to experimental.dat --view --save-report --save-gnuplot\n"
    )

    parser.add_argument("--file", required=True, help="Path to structure file.")
    parser.add_argument("--format", required=True, choices=["fdf", "struct_out"],
                        help="Input file format:\n"
                             "  fdf:        SIESTA structure input (%%block LatticeVectors etc.)\n"
                             "  struct_out: SIESTA post-relaxation output (.STRUCT_OUT)")
    parser.add_argument("--wavelength", type=str, default="CuKa",
                        help="X-ray source: a known name (CuKa, CuKa1, MoKa, CrKa, FeKa, CoKa, "
                             "AgKa, ...) or a wavelength in Ang as a plain number "
                             "(default: CuKa, 1.54184 Ang).")
    parser.add_argument("--two-theta-range", type=float, nargs=2, default=[0.0, 90.0],
                        metavar=("MIN", "MAX"),
                        help="2-theta range in degrees to compute (default: 0 90).")
    parser.add_argument("--top", type=int, default=None,
                        help="Only show the N strongest peaks in the printed table (default: "
                             "show all of them). The saved data file (--save-gnuplot) always has "
                             "all of them.")
    parser.add_argument("--compare-to", type=str, default=None, metavar="PATH",
                        help="Compare the simulated pattern against an experimental one: a plain "
                             "text file with two columns (2theta intensity), whitespace- or "
                             "comma-separated, '#' comments and blank lines allowed. Prints a "
                             "similarity score (pyxtal's cosine-weighted metric, 0-1, higher is "
                             "more similar); combine with --view for a visual overlay.")
    parser.add_argument("-o", "--output-dir", type=str, default=".",
                        help="Directory to write references.bib into (and the report/data+gnuplot "
                             "files, with --save-report/--save-gnuplot) (default: current "
                             "directory). Created if it doesn't exist.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the full run report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--save-gnuplot", action="store_true",
                        help="Also write xrd_pattern.dat (all peaks, pyxtal's own format) and "
                             "xrd_pattern.gplot (a stick-pattern plot script) together. Off by "
                             "default -- this tool used to write the data file unconditionally on "
                             "every run; that's no longer the case.")
    parser.add_argument("--view", action="store_true",
                        help="Show an interactive matplotlib preview of the pattern (or, with "
                             "--compare-to, an overlay of both patterns) before finishing. Off by "
                             "default.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-xrd {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    lo, hi = args.two_theta_range
    if lo >= hi:
        print(color_text("Error: --two-theta-range MIN must be less than MAX.", 'red'))
        sys.exit(1)

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    os.makedirs(args.output_dir, exist_ok=True)

    try:
        structure = structure_io.read_siesta_structure(args.file, args.format)
    except FileNotFoundError:
        print(color_text(f"Error: structure file '{args.file}' not found.", 'red'))
        sys.exit(1)
    except ValueError as e:
        print(color_text(f"Error: {e}", 'red'))
        sys.exit(1)

    try:
        wavelength, wavelength_label = xrd_core.resolve_wavelength(args.wavelength)
    except ValueError as e:
        print(color_text(f"Error: {e}", 'red'))
        sys.exit(1)

    try:
        xrd = xrd_core.compute_pattern(structure, wavelength=wavelength, two_theta_range=(lo, hi))
    except ValueError as e:
        print(color_text(f"Error: {e}", 'red'))
        sys.exit(1)

    rows = sorted(xrd.pxrd, key=lambda row: row[5], reverse=True)
    shown = rows if args.top is None else rows[:args.top]

    sym_info = symmetry_summary(structure)
    lat = structure.lattice

    similarity = None
    sim_profile = None
    experimental = None
    if args.compare_to:
        from pyxtal.XRD import Similarity

        try:
            experimental = xrd_core.read_experimental_pattern(args.compare_to)
        except (FileNotFoundError, ValueError) as e:
            print(color_text(f"Error: {e}", 'red'))
            sys.exit(1)

        sim_profile = xrd.get_profile()
        similarity = run_with_spinner(
            lambda: Similarity(sim_profile, experimental).value, label="Computing similarity")

    dat_path = gplot_path = None
    if args.save_gnuplot:
        dat_path = os.path.join(args.output_dir, "xrd_pattern.dat")
        gplot_path = os.path.join(args.output_dir, "xrd_pattern.gplot")
        xrd.save(dat_path)
        xrd_core.write_xrd_gnuplot(dat_path, gplot_path, structure.composition.reduced_formula,
                                    wavelength_label)

    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    print_dual(color_text("===== STB-XRD REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Input file      : {args.file}", f_out)
    print_dual(f"Format          : {args.format}", f_out)
    print_dual(f"Wavelength      : {wavelength_label} ({wavelength:.5f} Ang)", f_out)
    print_dual(f"2-theta range   : {lo} - {hi} deg", f_out)
    print_dual(f"Top peaks shown : {'all' if args.top is None else args.top}", f_out)
    print_dual(f"Output dir      : {args.output_dir}", f_out)
    print_dual(f"Save gnuplot    : {'yes' if args.save_gnuplot else 'no'}", f_out)
    print_dual(f"View (matplotlib): {'yes' if args.view else 'no'}", f_out)

    print_section("[1] STRUCTURE", f_out)
    print_dual(f"Formula         : {structure.composition.reduced_formula}", f_out)
    print_dual(f"Number of sites : {len(structure)}", f_out)
    if "Error" in sym_info:
        print_dual(color_text(f"[WARNING] Symmetry analysis failed: {sym_info['Error']}", 'yellow'), f_out)
    else:
        print_dual(f"Space Group     : {sym_info['Space Group']}", f_out)
        print_dual(f"Crystal System  : {sym_info['Crystal System']}", f_out)
        print_dual(f"Point Group     : {sym_info['Point Group']}", f_out)
        print_dual(f"Hall Symbol     : {sym_info['Hall Symbol']}", f_out)
        print_dual(f"Layer Group     : {sym_info['Layer Group']}", f_out)
    print_dual(f"Lattice a,b,c   : {lat.a:.4f}, {lat.b:.4f}, {lat.c:.4f} Ang", f_out)
    print_dual(f"Lattice angles  : alpha={lat.alpha:.2f}, beta={lat.beta:.2f}, "
               f"gamma={lat.gamma:.2f} deg", f_out)
    print_dual(f"Cell volume     : {lat.volume:.4f} Ang^3", f_out)
    print_dual(f"Density         : {structure.density:.4f} g/cm^3", f_out)

    print_section("[2] DIFFRACTION PATTERN", f_out)
    strongest = rows[0]
    min_d_row = min(rows, key=lambda r: r[1])
    print_dual(f"Peaks found in range : {len(rows)}", f_out)
    print_dual(f"Strongest peak       : 2theta={strongest[0]:.3f} deg, d={strongest[1]:.4f} Ang, "
               f"hkl=({int(strongest[2])} {int(strongest[3])} {int(strongest[4])}), "
               f"intensity={strongest[5]:.2f}", f_out)
    print_dual(f"Resolution (min d)   : {min_d_row[1]:.4f} Ang "
               f"(2theta={min_d_row[0]:.3f} deg)", f_out)
    if args.top is not None and args.top < len(rows):
        note = f"Showing top {len(shown)} of {len(rows)} peaks by intensity"
        note += " (full list in the data file)." if args.save_gnuplot else "."
        print_dual(note, f_out)
    table_rows = [
        ([f"{two_theta:.3f}", f"{d:.4f}", str(int(h)), str(int(k)), str(int(l)), f"{intensity:.2f}"], None)
        for two_theta, d, h, k, l, intensity in shown
    ]
    print_table(["2-theta (deg)", "d (Ang)", "h", "k", "l", "Intensity"], table_rows, f_out)

    if args.compare_to:
        print_section("[3] EXPERIMENTAL COMPARISON", f_out)
        print_dual(f"Experimental file  : {args.compare_to}", f_out)
        print_dual(f"Experimental points: {experimental.shape[1]}", f_out)
        print_dual(f"Similarity score   : {similarity:.4f} "
                   "(0-1, cosine-weighted, higher is more similar)", f_out)

    print_section("[4] OUTPUT DATA & PLOTS", f_out)
    if dat_path:
        print_dual(color_text(f"[OK] Data written to '{dat_path}' ({len(rows)} peaks).", 'green'), f_out)
        print_dual(color_text(f"[OK] Gnuplot script written to '{gplot_path}'.", 'green'), f_out)
    else:
        print_dual("Not written (off by default -- pass --save-gnuplot to write "
                   "xrd_pattern.dat/xrd_pattern.gplot).", f_out)

    print_section("[5] REFERENCES", f_out)
    bib_entries = [citations.SIESTA, citations.SIESTA_RECENT, citations.PYXTAL]
    citations.write_bib_file(os.path.join(args.output_dir, BIB_FILE), bib_entries)
    print_dual(color_text(
        f"[OK] Citations for the methods used in this run written to "
        f"'{os.path.join(args.output_dir, BIB_FILE)}' ({len(bib_entries)} entries).", 'green'), f_out)

    print_section("[6] SUMMARY & FILES", f_out)
    print_dual("Status         : OK", f_out)
    print_dual(f"References     : {os.path.join(args.output_dir, BIB_FILE)}", f_out)
    if dat_path:
        print_dual(f"Data           : {dat_path}", f_out)
        print_dual(f"Gnuplot script : {gplot_path}", f_out)
    if report_path:
        print_dual(f"Report         : {report_path}", f_out)

    if f_out:
        f_out.close()

    # --view runs last, after the report is fully printed/closed, so a
    # blocking matplotlib window never delays or hides it.
    if args.view:
        plot_matplotlib(xrd, structure, wavelength_label, (lo, hi),
                         sim_profile=sim_profile, experimental=experimental,
                         compare_to=args.compare_to, similarity=similarity)


if __name__ == "__main__":
    main()
