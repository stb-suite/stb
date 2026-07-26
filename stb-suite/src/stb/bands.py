#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "2.0.0"

import argparse
import os
from datetime import datetime
import matplotlib.pyplot as plt
from stb.core import citations
from stb.core.cli import color_text, show_intro, print_dual, print_section, print_table
from stb.core.siesta_bands import _is_gamma, read_data, shift_bands, cbm_vbm, read_eig_mesh

REPORT_FILE = "stb_bands_report.txt"
BIB_FILE = "references.bib"


def plot_gnuplot(high_sym, nspin, output_dir="."):
######################### PDF Plot
    # Build gnuplot-formatted labels (Gamma -> {/Symbol G}) in a local list
    # instead of mutating the caller's high_sym. xticsl below already wraps
    # each label in its own quotes, so no quotes are added here.
    gplot_labels = [
        "{/Symbol G}" if _is_gamma(label) else label for _, label in high_sym
    ]
    # Plot
    fileout=[]
    fileout.append('# Set terminal and output\n')
    fileout.append('set terminal pdfcairo enhanced font "Arial,25" size 8,8\n')
    fileout.append('set output "energy_bands.pdf"\n')
    fileout.append('\n')
    fileout.append('# Title and axis labels\n')
    fileout.append('#set title "Advanced Data Visualization Example" font ",16"\n')
    fileout.append('set ylabel "Energy (eV)" font "Arial,28" #offset -0.1,0\n')
    fileout.append('set xlabel "" font "Arial,28"\n')
    fileout.append('\n')
    fileout.append('# Ranges and scales\n')
    fileout.append('#set xrange [0:20]\n')
    fileout.append('set yrange [-20:20]\n')
    fileout.append('#set logscale y 10  # Set log scale for y-axis\n')
    fileout.append('\n')
    fileout.append('# Grid settings\n')
    fileout.append('#set grid xtics ytics mxtics mytics\n')
    fileout.append('#set mxtics 5  # Minor grid on x-axis\n')
    fileout.append('#set mytics 5  # Minor grid on y-axis\n')
    fileout.append('\n')
    fileout.append('# Tics and formatting\n')
    xticsl="set xtics ("
    for i in range(len(high_sym)):
        xticsl=xticsl+f' "{gplot_labels[i]}" {float(high_sym[i][0])} , '
    xticsl=xticsl[:-2]
    xticsl=xticsl+')  font "Arial,28"\n'
    fileout.append(xticsl)
    fileout.append('set ytics format "%.1f" font "Arial,28" \n')
    fileout.append('set grid xtics front')
    fileout.append('#set no key \n')
    fileout.append('\n')
    fileout.append('# Annotations\n')
    fileout.append('#set label "Critical point" at 10, 100 front font ",12" textcolor rgb "red"\n')
    fileout.append(f'set arrow from {float(high_sym[0][0])},0.0 to {float(high_sym[-1][0])} ,0.0 nohead dt 2 lc rgb "dark-gray" lw 4 back\n')
    for i in range(len(high_sym)-2):
        fileout.append(f'set arrow from {float(high_sym[i+1][0])} ,graph 0 to {float(high_sym[i+1][0])},graph 1 nohead dt 2 lc rgb "dark-gray" lw 4 back\n')
    fileout.append('\n')
    fileout.append('# Color palette and colorbar (useful for pm3d plots)\n')
    fileout.append('#set palette defined (0 "blue", 1 "green", 2 "yellow", 3 "red")\n')
    fileout.append('#set colorbox vertical user origin 0.9, 0.2 size 0.03, 0.6\n')
    fileout.append('#set cblabel "Colorbar (units)" font ",12"\n')
    fileout.append('\n')
    fileout.append('# Multiplot (e.g., additional smaller plots in the same figure)\n')
    fileout.append('#set multiplot layout 1,1 title "Multiplot with Color Palette" font ",14"\n')
    fileout.append('#unset multiplot\n')
    fileout.append('\n')
    fileout.append('# Line styles\n')
    fileout.append('set style line 1 lc rgb "#1f77b4" lt 1 lw 3 pt 7 ps 1.5   # Solid blue line with points\n')
    if nspin == 2:
        fileout.append('set style line 2 lc rgb "#d62728" lt 2 dt 2 lw 3        # Dashed red line (spin down)\n')
    else:
        fileout.append('#set style line 2 lc rgb "#ff7f0e" lt 2 lw 2 pt 5 ps 1.5   # Dashed orange line\n')
    fileout.append('#set style line 3 lc rgb "#2ca02c" lt 3 lw 2               # Solid green line\n')
    fileout.append('\n')
    fileout.append('# Plot commands\n')
    if nspin == 2:
        fileout.append('plot "bands_gnuplot.dat" using 1:2 with lines ls 1 title "Spin Up", '
                        '"" using 1:3 with lines ls 2 title "Spin Down"\n')
    else:
        fileout.append(f'plot "bands_gnuplot.dat" using 1:2 with lines ls 1 title "" \n')
    with open(os.path.join(output_dir, 'bands.gplot'), 'w') as file:
        file.writelines(fileout)
    return

def write_gnuplot_bands(dic_bands, nspin, output_dir="."):
    # define initial key
    key_init = next(iter(dic_bands))
    nbands = dic_bands[key_init].shape[1]
    # write the file
    file_name = os.path.join(output_dir, "bands_gnuplot.dat")
    with open(file_name, 'w') as file:
        for i in range(nbands):
            for key in dic_bands.keys():
                cols = "     ".join(f"{dic_bands[key][s][i]}" for s in range(nspin))
                file.write(f"{key}     {cols}\n")
            file.write("\n")
    return

def mesh_k_formatter(kpoints):
    # kpoints is None (plain mesh index) or an (nk, 3) array of Cartesian
    # (kx, ky, kz) in 1/Ang read from a companion .KP file -- see
    # read_eig_mesh().
    if kpoints is None:
        return lambda k: f"index {k}"
    def fmt(k):
        kx, ky, kz = kpoints[k]
        return f"kx={kx:.6f}, ky={ky:.6f}, kz={kz:.6f} 1/Ang (index {k})"
    return fmt

# --- Report formatting -----------------------------------------------
# Same numbered-section report style as the rest of the suite's newer
# tools ([0] RUN METADATA ... [N] SUMMARY & FILES via print_section/
# print_dual/print_table) -- replaces the old ad hoc plain-text layout.
# The underlying numbers/physics are unchanged; only how they're printed.

def extrema_rows(combined, k_format):
    """Returns print_table rows for one (vbm, cbm, vbm_k, cbm_k,
    indirect_gap, gap_type, direct_gap, direct_k) tuple from cbm_vbm()."""
    vbm, cbm, vbm_k, cbm_k, indirect_gap, gap_type, direct_gap, direct_k = combined
    return [
        (["VBM", f"{vbm:.6f} eV (k = {k_format(vbm_k)})"], None),
        (["CBM", f"{cbm:.6f} eV (k = {k_format(cbm_k)})"], None),
        (["Indirect gap", f"{indirect_gap:.6f} eV (fundamental: CBM - VBM, any k)"], None),
        (["Direct gap", f"{direct_gap:.6f} eV (same-k minimum, k = {k_format(direct_k)})"], None),
        (["Gap type", gap_type], 'yellow' if gap_type == "Metallic" else None),
    ]


def print_spin_channels(result, k_format, f_out):
    """Prints one extrema table per spin channel plus the half-metallic
    flag -- shared by the line ([3]) and mesh ([4]) sections."""
    for s, spin_result in enumerate(result["spins"]):
        spin_name = "Spin up" if s == 0 else "Spin down"
        print_dual(f"{spin_name}:", f_out)
        print_table(["Quantity", "Value"], extrema_rows(spin_result, k_format), f_out)
    half_metallic = result["half_metallic"]
    line = f"Half-metallic : {'Yes' if half_metallic else 'No'}"
    print_dual(color_text(line, 'green') if half_metallic else line, f_out)


def main():
    parser = argparse.ArgumentParser(
        description="Process the band structure data.",
        epilog="Example usage:\n"
               "  stb_bands --label siesta --shift fermi\n"
               "  stb_bands --file siesta.bands --shift fermi\n"
               "  stb_bands --file siesta.bands --shift manual --manual-value 0.5\n"
               "  stb_bands --file siesta.bands --shift fermi --gap-tol 0.05\n"
               "  stb_bands --file siesta.bands --shift fermi --eig-file siesta.EIG\n"
               "  stb_bands --file siesta.bands --shift fermi --save-report",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("--label", type=str, default=None,
                        help="SIESTA output label (e.g. 'siesta'). Shorthand for "
                             "--file <label>.bands, and auto-detects <label>.EIG next "
                             "to it for the mesh (k-grid) gap comparison -- if that "
                             ".EIG isn't there, the mesh comparison is skipped rather "
                             "than treated as an error. Also auto-detects <label>.KP "
                             "(mesh k-point coordinates) alongside the .EIG, if present. "
                             "Mutually exclusive with --file/--eig-file/--kp-file.")

    parser.add_argument("--file",  dest="input_file", type=str, default=None,
                        help="Path to the input file containing band structure data (e.g., siesta.bands). "
                             "Alternative to --label when you need an explicit path.")

    parser.add_argument("--shift", type=str, choices=["vbm", "cbm", "fermi", "manual"], required=True,
                        help="Reference energy shift:\n"
                             "  - 'vbm'    : Valence Band Maximum\n"
                             "  - 'cbm'    : Conduction Band Minimum\n"
                             "  - 'fermi'  : Fermi level\n"
                             "  - 'manual' : Custom shift value (requires --manual-value).")

    parser.add_argument("--manual-value", type=float,
                        help="Custom energy shift value (required if --shift manual is used).")

    parser.add_argument("--gap-tol", type=float, default=0.01,
                        help="Energy tolerance in eV below which a gap is classified as Metallic "
                             "(default: 0.01).")

    parser.add_argument("--eig-file", type=str, default=None,
                        help="Optional SIESTA .EIG file (full k-mesh eigenvalues from the SCF run). "
                             "When given, compares the k-mesh gap against the k-path gap from --file "
                             "-- the high-symmetry path may not pass through the true VBM/CBM. Requires sisl.")

    parser.add_argument("--kp-file", type=str, default=None,
                        help="Optional SIESTA .KP file matching --eig-file's k-mesh (Cartesian "
                             "kx,ky,kz in 1/Ang). When given, mesh VBM/CBM/direct-gap k-points are "
                             "reported as (kx, ky, kz) instead of a bare mesh index. Requires --eig-file.")

    parser.add_argument("-o", "--output-dir", type=str, default=".",
                        help="Directory to write bands_gnuplot.dat and bands.gplot into "
                             "(and stb_bands_report.txt/references.bib, with --save-report) "
                             "(default: current directory). Created if it doesn't exist.")

    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the full run report to {REPORT_FILE}. Off by default.")

    parser.add_argument("-v", "--version", action="version", version=f"stb-bands {VERSION}")

    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()


    # verify the manual value
    if args.shift == "manual" and args.manual_value is None:
        parser.error("--manual-value is required when --shift is set to 'manual'.")

    # Resolve --label into --file/--eig-file/--kp-file, or fall back to explicit paths.
    if args.label:
        if args.input_file or args.eig_file or args.kp_file:
            parser.error("--label cannot be combined with --file/--eig-file/--kp-file.")
        args.input_file = f"{args.label}.bands"
        eig_candidate = f"{args.label}.EIG"
        eig_auto_note = None
        if os.path.isfile(eig_candidate):
            args.eig_file = eig_candidate
            kp_candidate = f"{args.label}.KP"
            if os.path.isfile(kp_candidate):
                args.kp_file = kp_candidate
            else:
                eig_auto_note = (f"No '{kp_candidate}' found next to the label; mesh "
                                  "VBM/CBM/direct-gap k-points will be reported by mesh index only.")
        else:
            eig_auto_note = f"No '{eig_candidate}' found next to the label; mesh (k-grid) gap comparison will be skipped."
    else:
        eig_auto_note = None
        if not args.input_file:
            parser.error("one of --label or --file is required.")
        elif args.kp_file and not args.eig_file:
            parser.error("--kp-file requires --eig-file.")

    if args.intro == True:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    shift_desc = args.shift if args.shift != "manual" else f"manual ({args.manual_value} eV)"

    # Read (and validate) the input data BEFORE creating --output-dir or
    # opening the report file, so a bad input file doesn't leave an empty
    # directory/empty report behind (same fix as dos.py's -o). Everything
    # printed below -- including [0] RUN METADATA -- only happens once
    # this has already succeeded.
    fermi_energy, high_sym, dic_bands, nspin = read_data(args.input_file)
    result = cbm_vbm(fermi_energy, dic_bands, nspin, args.gap_tol)

    os.makedirs(args.output_dir, exist_ok=True)
    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    print_dual(color_text("===== STB-BANDS REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time      : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Input file     : {args.input_file}", f_out)
    print_dual(f"Shift mode     : {shift_desc}", f_out)
    print_dual(f"Gap tolerance  : {args.gap_tol} eV", f_out)
    if args.eig_file:
        print_dual(f"Mesh file      : {args.eig_file}" + (f" (k-points: {args.kp_file})" if args.kp_file else " (no --kp-file: mesh index only)"), f_out)
    elif eig_auto_note:
        print_dual(color_text(f"[INFO] {eig_auto_note}", 'yellow'), f_out)
    print_dual(f"Output dir     : {args.output_dir}", f_out)

    nk = len(dic_bands)
    nbands = next(iter(dic_bands.values())).shape[1]

    print_section("[1] INPUT DATA", f_out)
    print_dual(f"Fermi energy   : {fermi_energy:.6f} eV", f_out)
    print_dual(f"Spin channels  : {nspin} ({'non-polarized' if nspin == 1 else 'polarized'})", f_out)
    print_dual(f"Bands x k-points: {nbands} x {nk}", f_out)
    print_table(["k", "Label"], [([f"{float(k):.6f}", label], None) for k, label in high_sym], f_out)

    print_section("[2] BAND GAP ANALYSIS (k-path)", f_out)
    print_table(["Quantity", "Value"], extrema_rows(result["combined"], result["k_format"]), f_out)

    if nspin == 2:
        print_section("[3] SPIN-POLARIZED ANALYSIS (k-path)", f_out)
        print_spin_channels(result, result["k_format"], f_out)

    mesh_result = None
    mesh_warnings = []
    print_section("[4] MESH VS LINE COMPARISON", f_out)
    if not args.eig_file:
        print_dual("Not available -- no k-mesh eigenvalues were given (--eig-file, or "
                   "<label>.EIG auto-detected next to --label). The indirect gap above "
                   "may be an overestimate if the true VBM/CBM lie off the high-symmetry path.", f_out)
    else:
        fermi_mesh, dic_mesh, nspin_mesh, mesh_kpoints = read_eig_mesh(args.eig_file, args.kp_file)
        if abs(fermi_mesh - fermi_energy) > 1e-3:
            mesh_warnings.append(
                f"Fermi energy mismatch between --file ({fermi_energy:.6f} eV) and "
                f"--eig-file ({fermi_mesh:.6f} eV); the two inputs may be from different calculations."
            )
        if nspin_mesh != nspin:
            mesh_warnings.append(
                f"nspin mismatch between --file (nspin={nspin}) and --eig-file "
                f"(nspin={nspin_mesh}); mesh comparison uses combined values only."
            )
        for w in mesh_warnings:
            print_dual(color_text(f"[WARNING] {w}", 'yellow'), f_out)
        mesh_result = cbm_vbm(fermi_mesh, dic_mesh, nspin_mesh, args.gap_tol,
                               k_format=mesh_k_formatter(mesh_kpoints))
        print_dual(f"Mesh k-points  : {len(dic_mesh)}"
                   + ("" if args.kp_file else " (no --kp-file: mesh index only)"), f_out)
        print_table(["Quantity", "Value"], extrema_rows(mesh_result["combined"], mesh_result["k_format"]), f_out)
        if mesh_result["spins"]:
            print_spin_channels(mesh_result, mesh_result["k_format"], f_out)

        m_indirect = mesh_result["combined"][4]
        line_indirect = result["combined"][4]
        diff = m_indirect - line_indirect
        print_dual(f"Indirect gap diff (mesh - line): {diff:.6f} eV", f_out)
        if abs(diff) < args.gap_tol:
            print_dual(color_text("[OK] Line and mesh indirect gaps agree within tolerance.", 'green'), f_out)
        elif diff < 0:
            print_dual(color_text(
                "[NOTE] Mesh gap is smaller than the line gap: the true VBM/CBM likely lie "
                "off the high-symmetry path. Trust the mesh value for the fundamental gap.", 'yellow'), f_out)
        else:
            print_dual(color_text(
                "[WARNING] Mesh gap is larger than the line gap: the k-mesh may be too coarse "
                "to capture the extrema found on the path. Consider a denser --eig-file mesh.", 'yellow'), f_out)

    vbm, cbm, vbm_k, cbm_k, indirect_gap, gap_type, direct_gap, direct_k = result["combined"]

    if args.shift == "vbm":
        rshift = vbm
    elif args.shift == "cbm":
        rshift = cbm
    elif args.shift == "fermi":
        rshift = fermi_energy
    elif args.shift == "manual":
        rshift = args.manual_value

    shifted_bands = shift_bands(dic_bands, rshift)
    write_gnuplot_bands(shifted_bands, nspin, args.output_dir)
    plot_gnuplot(high_sym, nspin, args.output_dir)

    print_section("[5] WRITING OUTPUT FILES", f_out)
    print_dual(color_text(
        f"[OK] Gnuplot data written to '{os.path.join(args.output_dir, 'bands_gnuplot.dat')}' "
        f"and '{os.path.join(args.output_dir, 'bands.gplot')}' (run gnuplot on the latter for a PDF plot).",
        'green'), f_out)

    print_section("[6] REFERENCES", f_out)
    bib_entries = [citations.SIESTA, citations.SIESTA_RECENT]
    citations.write_bib_file(os.path.join(args.output_dir, BIB_FILE), bib_entries)
    print_dual(color_text(
        f"[OK] Citations for the methods used in this run written to "
        f"'{os.path.join(args.output_dir, BIB_FILE)}' ({len(bib_entries)} entries).", 'green'), f_out)

    if mesh_result is not None:
        best_indirect, best_type = mesh_result["combined"][4], mesh_result["combined"][5]
        best_source = "the k-mesh (denser sampling than the path)"
        half_metallic = mesh_result["half_metallic"] if nspin == 2 else False
    else:
        best_indirect, best_type = indirect_gap, gap_type
        best_source = "the k-path (no k-mesh was provided)"
        half_metallic = result["half_metallic"] if nspin == 2 else False

    print_section("[7] SUMMARY & FILES", f_out)
    print_dual(f"Status         : OK", f_out)
    print_dual(f"Best fundamental (indirect) gap: {best_indirect:.6f} eV ({best_type}), from {best_source}.", f_out)
    if nspin == 2 and half_metallic:
        print_dual(color_text(
            "Half-metallic character detected (one spin channel metallic, the other has a gap).",
            'green'), f_out)
    print_dual(f"Gnuplot data   : {os.path.join(args.output_dir, 'bands_gnuplot.dat')}", f_out)
    print_dual(f"Gnuplot script : {os.path.join(args.output_dir, 'bands.gplot')}", f_out)
    print_dual(f"References     : {os.path.join(args.output_dir, BIB_FILE)}", f_out)
    if report_path:
        print_dual(f"Report         : {report_path}", f_out)

    if f_out:
        f_out.close()

    # Interactive matplotlib preview runs last, after every report/file
    # write above, so a blocking GUI window never delays or hides them.
    plot(shifted_bands, high_sym, nspin)

def plot(dic, custom_ticks, nspin):
    # Organize the data
    x_values = sorted(dic.keys())
    nbands = next(iter(dic.values())).shape[1]
    # Organize the high symmetries points
    tick_positions = [float(t[0]) for t in custom_ticks]
    tick_labels = ["Γ" if _is_gamma(t[1]) else t[1] for t in custom_ticks]
    # plot: spin channels overlaid in a single panel (up: solid blue,
    # down: dashed red), one legend entry per channel rather than per band.
    spin_styles = [("blue", "-", "Spin Up"), ("red", "--", "Spin Down")]
    plt.figure(figsize=(8, 6))
    plt.xticks(tick_positions, tick_labels)
    for s in range(nspin):
        color, linestyle, label = spin_styles[s]
        for i in range(nbands):
            y_vals = [dic[x][s][i] for x in x_values]
            plt.plot(x_values, y_vals, color=color, linestyle=linestyle,
                     label=label if (nspin == 2 and i == 0) else None)
    # vertical lines in High Symmetries
    for pos in tick_positions:
        plt.axvline(x=pos, color='gray', linestyle='--', linewidth=1)
    # plot limits
    plt.ylim(-20, 20)
    plt.ylabel("Energy")
    plt.grid(True)
    if nspin == 2:
        plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
