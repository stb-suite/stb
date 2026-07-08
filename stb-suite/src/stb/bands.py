#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.9.1"

import argparse
import os
from datetime import datetime
import numpy as np
import re
import matplotlib.pyplot as plt
from stb.core.cli import color_text, show_intro
from stb.core.deps import require_sisl

def _is_gamma(label):
    clean_text = re.sub(r'[^a-zA-Z\s]', '', label).lower()
    return "gamma" in clean_text.split()

def plot_gnuplot(high_sym, nspin):
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
    with open('bands.gplot', 'w') as file:
        file.writelines(fileout)
    return

def read_data(file_path = "siesta.bands"):
    # SIESTA .bands format:
    #   line 0: Fermi energy
    #   line 1: kmin kmax
    #   line 2: emin emax
    #   line 3: nbands nspin nk
    #   body:   nk k-blocks (first line has the k value, then as many
    #           continuation lines as needed to reach nbands*nspin values)
    #   footer: a line with the count of high-symmetry points, followed by
    #           that many "position 'LABEL'" lines
    with open(file_path, "r") as f:
        lines = f.readlines()

    fermi_energy = float(lines[0].split()[0])
    nbands, nspin, nk = (int(x) for x in lines[3].split())
    values_per_k = nbands * nspin

    dic_bands = {}
    idx = 4
    for _ in range(nk):
        tokens = lines[idx].split()
        key = float(tokens[0])
        values = tokens[1:]
        idx += 1
        while len(values) < values_per_k:
            values.extend(lines[idx].split())
            idx += 1
        # First nbands values are the spin-up channel, the next nbands (if
        # any) are spin-down -- confirmed against sisl's own bandsSileSiesta
        # ("eb[ik, :, :] = l.reshape(ns, no)"), not interleaved per band.
        dic_bands[key] = np.array(values, dtype=float).reshape(nspin, nbands)

    # Whatever remains is deterministically the high-symmetry footer: we
    # already consumed exactly nk k-blocks by token count, so there is no
    # ambiguity left between a continuation line and the count line.
    while idx < len(lines) and not lines[idx].split():
        idx += 1
    n_high_sym = int(lines[idx].split()[0])
    idx += 1

    high_sym = []
    for _ in range(n_high_sym):
        parts = lines[idx].split()
        label = " ".join(parts[1:]).strip("'\"")
        high_sym.append([parts[0], label])
        idx += 1

    return fermi_energy, high_sym, dic_bands, nspin

def read_eig_mesh(eig_file, kp_file=None):
    # A SIESTA .EIG file holds the eigenvalues at every k-point of the SCF
    # k-mesh (not just the high-symmetry path in .bands), which is what a
    # gap comparison against the path needs. Reuse sisl's own reader
    # (sisl/io/siesta/eig.py) instead of hand-rolling a second parser --
    # it already returns exactly (ns, nk, nb), Ef-shifted.
    sisl = require_sisl()
    sile = sisl.get_sile(eig_file)
    fermi_energy = sile.read_fermi_level()
    eigs = sile.read_data()
    nspin = eigs.shape[0]
    nk = eigs.shape[1]
    # Same (nspin, nbands)-per-k convention as read_data(); add Ef back so
    # values are absolute eV, matching the rest of this module.
    dic_mesh = {ik: eigs[:, ik, :] + fermi_energy for ik in range(nk)}

    kpoints = None
    if kp_file:
        # The .EIG file itself has no k-vectors, only a bare mesh index --
        # the .KP file (same calculation) holds the actual Cartesian
        # (kx, ky, kz), in 1/Ang (sisl converts from the file's raw 1/Bohr).
        kp_sile = sisl.get_sile(kp_file)
        kpoints, _weights = kp_sile.read_data()
        if kpoints.shape[0] != nk:
            raise ValueError(
                f"--kp-file '{kp_file}' has {kpoints.shape[0]} k-points but "
                f"'{eig_file}' has {nk} -- they must be from the same calculation."
            )
    return fermi_energy, dic_mesh, nspin, kpoints

def write_gnuplot_bands(dic_bands, nspin):
    # define initial key
    key_init = next(iter(dic_bands))
    nbands = dic_bands[key_init].shape[1]
    # write the file
    file_name="bands_gnuplot.dat"
    with open(file_name, 'w') as file:
        for i in range(nbands):
            for key in dic_bands.keys():
                cols = "     ".join(f"{dic_bands[key][s][i]}" for s in range(nspin))
                file.write(f"{key}     {cols}\n")
            file.write("\n")
    return

def _band_extrema(values_by_k, fermi_energy, gap_tol):
    # Start value as infinity
    vbm = -np.inf
    cbm = np.inf
    vbm_k = None
    cbm_k = None
    # Direct gap: the smallest CBM(k) - VBM(k) at one and the same k --
    # tracked alongside the (possibly different-k) global VBM/CBM below.
    direct_gap = np.inf
    direct_k = None
    for k, band in values_by_k.items():
        below = band[band <= fermi_energy]
        above = band[band > fermi_energy]
        if below.size > 0:
            local_vbm = np.nanmax(below)
            if local_vbm > vbm:
                vbm = local_vbm
                vbm_k = k
        if above.size > 0:
            local_cbm = np.nanmin(above)
            if local_cbm < cbm:
                cbm = local_cbm
                cbm_k = k
            if below.size > 0:
                local_gap = local_cbm - local_vbm
                if local_gap < direct_gap:
                    direct_gap = local_gap
                    direct_k = k
    if vbm_k is None or cbm_k is None:
        raise ValueError(
            f"No occupied/empty states found on {'both sides' if vbm_k is None and cbm_k is None else ('the occupied side' if vbm_k is None else 'the empty side')} "
            f"of Fermi energy {fermi_energy:.6f} eV -- check that this Fermi energy actually matches the eigenvalue file."
        )
    # Indirect (= fundamental) gap: CBM - VBM regardless of whether they sit
    # at the same k. Always <= direct_gap, since direct_gap is the same
    # quantity restricted to matching k. When they coincide (within
    # gap_tol) the fundamental gap is itself direct.
    indirect_gap = cbm - vbm if cbm > vbm else 0.0  # Avoid negative values
    if indirect_gap < gap_tol:
        gap_type = "Metallic"
    elif direct_k is not None and (direct_gap - indirect_gap) < gap_tol:
        gap_type = "Direct"
    else:
        gap_type = "Indirect"
    return vbm, cbm, vbm_k, cbm_k, indirect_gap, gap_type, direct_gap, direct_k

def _default_k_format(k):
    return f"{k:.6f}"

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

def cbm_vbm(fermi_energy, dic_bands, nspin, gap_tol=0.01, k_format=None):
    # Pure computation, no printing -- callers render the result (console
    # and/or bands_analysis.txt) via write_analysis_report(), so the same
    # numbers are never formatted two different ways in two places.
    k_format = k_format or _default_k_format
    combined = _band_extrema(
        {k: arr.reshape(-1) for k, arr in dic_bands.items()}, fermi_energy, gap_tol
    )
    result = {"combined": combined, "spins": [], "half_metallic": False, "k_format": k_format}
    if nspin == 2:
        for s in range(nspin):
            spin_result = _band_extrema(
                {k: arr[s] for k, arr in dic_bands.items()}, fermi_energy, gap_tol
            )
            result["spins"].append(spin_result)
        gaps = [spin_result[4] for spin_result in result["spins"]]  # indirect gap
        result["half_metallic"] = min(gaps) < gap_tol and max(gaps) >= gap_tol

    return result

# --- Report formatting -----------------------------------------------
# Plain-text report, but laid out like a short write-up: a title block,
# then one category per physically distinct source of data (line vs.
# mesh), a comparison section only when both are available, and a
# one-line summary at the end that states the headline result plainly.

_WIDTH = 74
_LABEL_W = 14

def _rule(char="-"):
    return char * _WIDTH

def _kv(label, value, width=_LABEL_W):
    return f"{label:<{width}}: {value}"

def _extrema_lines(combined, k_format, indent=""):
    vbm, cbm, vbm_k, cbm_k, indirect_gap, gap_type, direct_gap, direct_k = combined
    return [
        indent + _kv("VBM", f"{vbm:.6f} eV  (k = {k_format(vbm_k)})"),
        indent + _kv("CBM", f"{cbm:.6f} eV  (k = {k_format(cbm_k)})"),
        indent + _kv("Indirect gap", f"{indirect_gap:.6f} eV  (fundamental: CBM - VBM, any k)"),
        indent + _kv("Direct gap", f"{direct_gap:.6f} eV  (same-k minimum, k = {k_format(direct_k)})"),
        indent + _kv("Gap type", gap_type),
    ]

def _spin_block(result, k_format):
    lines = []
    for s, spin_result in enumerate(result["spins"]):
        spin_name = "Spin up" if s == 0 else "Spin down"
        lines.append(f"{spin_name}:")
        lines.extend(_extrema_lines(spin_result, k_format, indent="  "))
        lines.append("")
    lines.append(_kv("Half-metallic", "Yes" if result["half_metallic"] else "No"))
    return lines

def write_analysis_report(fermi_energy, result, nspin, mesh_result=None, gap_tol=0.01,
                           mesh_warnings=None, input_file=None, eig_file=None, kp_file=None):
    combined = result["combined"]
    k_format = result["k_format"]
    indirect_gap, gap_type = combined[4], combined[5]

    lines = []
    lines.append(_rule("="))
    lines.append("BAND STRUCTURE ANALYSIS REPORT - STB Suite".center(_WIDTH))
    lines.append(_rule("="))
    lines.append("")
    lines.append(_kv("Generated", datetime.now().strftime("%Y-%m-%d %H:%M:%S"), width=16))
    lines.append(_kv("Fermi energy", f"{fermi_energy:.6f} eV", width=16))
    lines.append(_kv("Spin channels", f"{nspin} ({'non-polarized' if nspin == 1 else 'polarized'})", width=16))
    lines.append(_kv("Gap tolerance", f"{gap_tol:.6f} eV", width=16))

    lines.append("")
    lines.append(_rule())
    lines.append("LINE (K-PATH) RESULTS")
    lines.append(_rule())
    if input_file:
        lines.append(_kv("Source file", input_file))
        lines.append("")
    lines.extend(_extrema_lines(combined, k_format))
    if nspin == 2:
        lines.append("")
        lines.extend(_spin_block(result, k_format))

    lines.append("")
    lines.append(_rule())
    lines.append("MESH (K-GRID) RESULTS")
    lines.append(_rule())
    if mesh_result is None:
        lines.append("Not available -- no k-mesh eigenvalues were given (--eig-file, or")
        lines.append("<label>.EIG auto-detected next to --label). The indirect gap above")
        lines.append("may be an overestimate if the true VBM/CBM lie off the high-symmetry path.")
    else:
        m_combined = mesh_result["combined"]
        m_k_format = mesh_result["k_format"]
        if eig_file:
            lines.append(_kv("Source file", eig_file))
        if kp_file:
            lines.append(_kv("k-points file", kp_file))
        else:
            lines.append("k-points file : not given -- mesh k-points shown by mesh index only.")
        lines.append("")
        lines.extend(_extrema_lines(m_combined, m_k_format))
        if mesh_result["spins"]:
            lines.append("")
            lines.extend(_spin_block(mesh_result, m_k_format))

    if mesh_result is not None:
        m_indirect, m_gap_type = mesh_result["combined"][4], mesh_result["combined"][5]
        diff = m_indirect - indirect_gap
        lines.append("")
        lines.append(_rule())
        lines.append("MESH vs LINE COMPARISON")
        lines.append(_rule())
        for w in (mesh_warnings or []):
            lines.append(f"[WARNING] {w}")
        lines.append(_kv("Indirect gap diff", f"{diff:.6f} eV  (mesh - line)", width=20))
        if abs(diff) < gap_tol:
            lines.append("Line and mesh indirect gaps agree within tolerance.")
        elif diff < 0:
            lines.append("[WARNING] Mesh gap is smaller than the line gap: the true VBM/CBM "
                          "likely lie off the high-symmetry path. Trust the mesh value for the "
                          "fundamental gap.")
        else:
            lines.append("[WARNING] Mesh gap is larger than the line gap: the k-mesh may be too "
                          "coarse to capture the extrema found on the path. Consider a denser "
                          "--eig-file mesh.")

    lines.append("")
    lines.append(_rule())
    lines.append("SUMMARY")
    lines.append(_rule())
    if mesh_result is not None:
        best_indirect, best_type = mesh_result["combined"][4], mesh_result["combined"][5]
        best_source = "the k-mesh (denser sampling than the path)"
        half_metallic = mesh_result["half_metallic"] if nspin == 2 else False
    else:
        best_indirect, best_type = indirect_gap, gap_type
        best_source = "the k-path (no k-mesh was provided)"
        half_metallic = result["half_metallic"] if nspin == 2 else False
    lines.append(f"Best fundamental (indirect) gap estimate: {best_indirect:.6f} eV ({best_type}), from {best_source}.")
    if nspin == 2 and half_metallic:
        lines.append("Half-metallic character detected (one spin channel metallic, the other has a gap).")
    lines.append(_rule("="))

    content = "\n".join(lines) + "\n"
    with open("bands_analysis.txt", "w") as f:
        f.write(content)
    return content

    content = "\n".join(lines) + "\n"
    with open("bands_analysis.txt", "w") as f:
        f.write(content)
    return content

def shift_bands(dic, val):
    return {k: arr - val for k, arr in dic.items()}

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


def main():
    parser = argparse.ArgumentParser(
        description="Process the band structure data.",
        epilog="Example usage:\n"
               "  stb_bands --label siesta --shift fermi\n"
               "  stb_bands --file siesta.bands --shift fermi\n"
               "  stb_bands --file siesta.bands --shift manual --manual-value 0.5\n"
               "  stb_bands --file siesta.bands --shift fermi --gap-tol 0.05\n"
               "  stb_bands --file siesta.bands --shift fermi --eig-file siesta.EIG",
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
        if os.path.isfile(eig_candidate):
            args.eig_file = eig_candidate
            kp_candidate = f"{args.label}.KP"
            if os.path.isfile(kp_candidate):
                args.kp_file = kp_candidate
            else:
                print(f"[INFO] No '{kp_candidate}' found next to the label; "
                      "mesh VBM/CBM/direct-gap k-points will be reported by mesh index only.")
        else:
            print(f"[INFO] No '{eig_candidate}' found next to the label; "
                  "mesh (k-grid) gap comparison will be skipped.")
    elif not args.input_file:
        parser.error("one of --label or --file is required.")
    elif args.kp_file and not args.eig_file:
        parser.error("--kp-file requires --eig-file.")

    if args.intro == True:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("BANDS:", 'bold'))
    print("-"*60)

    # Condition to shift the band structure
    print("[INFO] Reading band structure data ...")
    fermi_energy,high_sym,dic_bands,nspin=read_data(args.input_file)
    result = cbm_vbm(fermi_energy, dic_bands, nspin, args.gap_tol)

    mesh_result = None
    mesh_warnings = []
    if args.eig_file:
        print("[INFO] Reading k-mesh eigenvalues (--eig-file) ...")
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
        mesh_result = cbm_vbm(fermi_mesh, dic_mesh, nspin_mesh, args.gap_tol,
                               k_format=mesh_k_formatter(mesh_kpoints))

    report = write_analysis_report(fermi_energy, result, nspin, mesh_result, args.gap_tol, mesh_warnings,
                                    input_file=args.input_file, eig_file=args.eig_file, kp_file=args.kp_file)
    print("\n" + report)
    print("[INFO] Full report saved to bands_analysis.txt")
    vbm, cbm, vbm_k, cbm_k, indirect_gap, gap_type, direct_gap, direct_k = result["combined"]

    if args.shift == "vbm":
        rshift = vbm
    elif args.shift == "cbm":
        rshift = cbm
    elif args.shift == "fermi":
        rshift = fermi_energy
    elif args.shift == "manual":
        rshift = args.manual_value
    print("[INFO] Write files...")
    print("[WARNING] \n")

    shifted_bands = shift_bands(dic_bands, rshift)
    write_gnuplot_bands(shifted_bands, nspin)
    plot(shifted_bands, high_sym, nspin)
    plot_gnuplot(high_sym, nspin)

    print("\n[INFO] Complete job!") 
    print("\n"+"-"*60)
    print(color_text("Bands found! But still no sign of Metallica.\n\n", 'bold'))

if __name__ == "__main__":
    main()
