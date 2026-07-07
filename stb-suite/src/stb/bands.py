#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.9.1"

import os
import sys
import warnings
import subprocess
from time import sleep
import argparse
import textwrap
from typing import List, Dict
import numpy as np
import re
import argparse
import matplotlib.pyplot as plt
from stb.core.cli import COLORS, color_text, show_intro
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

def read_eig_mesh(eig_file):
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
    return fermi_energy, dic_mesh, nspin

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
    band_gap = cbm - vbm if cbm > vbm else 0.0  # Avoid negative values
    if band_gap < gap_tol:
        gap_type = "Metallic"
    elif vbm_k == cbm_k:
        gap_type = "Direct"
    else:
        gap_type = "Indirect"
    return vbm, cbm, vbm_k, cbm_k, band_gap, gap_type

def cbm_vbm(fermi_energy, high_sym, dic_bands, nspin, gap_tol=0.01):
    combined = _band_extrema(
        {k: arr.reshape(-1) for k, arr in dic_bands.items()}, fermi_energy, gap_tol
    )
    vbm, cbm, vbm_k, cbm_k, band_gap, gap_type = combined
    print(f"[INFO] Fermi: {fermi_energy} \n[INFO] VBM: {vbm:.6f} (k={vbm_k}) \n[INFO] CBM: {cbm:.6f} (k={cbm_k})\n[INFO] Band Gap: {band_gap:.6f} eV ({gap_type})")

    result = {"combined": combined, "spins": [], "half_metallic": False}
    if nspin == 2:
        for s in range(nspin):
            spin_result = _band_extrema(
                {k: arr[s] for k, arr in dic_bands.items()}, fermi_energy, gap_tol
            )
            result["spins"].append(spin_result)
            spin_name = "Up" if s == 0 else "Down"
            svbm, scbm, svbm_k, scbm_k, sgap, sgap_type = spin_result
            print(f"[INFO] Spin {spin_name}: VBM={svbm:.6f} (k={svbm_k}) CBM={scbm:.6f} (k={scbm_k}) Gap={sgap:.6f} eV ({sgap_type})")
        gaps = [spin_result[4] for spin_result in result["spins"]]
        result["half_metallic"] = min(gaps) < gap_tol and max(gaps) >= gap_tol
        if result["half_metallic"]:
            print("[INFO] Half-metallic character detected (one spin channel metallic, the other has a gap)")

    return result

def write_analysis_report(fermi_energy, result, nspin, mesh_result=None, gap_tol=0.01):
    vbm, cbm, vbm_k, cbm_k, band_gap, gap_type = result["combined"]
    lines = [
        "BAND STRUCTURE ANALYSIS REPORT - STB Suite",
        "-"*60,
        f"Fermi energy : {fermi_energy:.6f} eV",
        f"VBM          : {vbm:.6f} eV (k = {vbm_k})",
        f"CBM          : {cbm:.6f} eV (k = {cbm_k})",
        f"Band gap     : {band_gap:.6f} eV",
        f"Gap type     : {gap_type}",
    ]
    if nspin == 2:
        lines.append("")
        lines.append("== Per-spin channel ==")
        for s, (svbm, scbm, svbm_k, scbm_k, sgap, sgap_type) in enumerate(result["spins"]):
            spin_name = "Spin Up" if s == 0 else "Spin Down"
            lines.append(f"{spin_name:<10}: VBM={svbm:.6f} eV (k={svbm_k})  CBM={scbm:.6f} eV (k={scbm_k})  Gap={sgap:.6f} eV ({sgap_type})")
        lines.append("")
        lines.append(f"Half-metallic: {'Yes' if result['half_metallic'] else 'No'}")

    if mesh_result is not None:
        m_vbm, m_cbm, m_vbm_k, m_cbm_k, m_gap, m_gap_type = mesh_result["combined"]
        diff = m_gap - band_gap
        lines.append("")
        lines.append("== Mesh (k-grid) vs Line (k-path) comparison ==")
        lines.append(f"Line gap (k-path)  : {band_gap:.6f} eV ({gap_type})")
        lines.append(f"Mesh gap (k-grid)  : {m_gap:.6f} eV ({m_gap_type})")
        lines.append(f"Difference (mesh - line): {diff:.6f} eV")
        if abs(diff) < gap_tol:
            lines.append("Line and mesh gaps agree within tolerance.")
        elif diff < 0:
            lines.append("[WARNING] Mesh gap is smaller than the line gap: the true VBM/CBM "
                          "likely lie off the high-symmetry path. Trust the mesh value for the "
                          "fundamental gap.")
        else:
            lines.append("[WARNING] Mesh gap is larger than the line gap: the k-mesh may be too "
                          "coarse to capture the extrema found on the path. Consider a denser "
                          "--eig-file mesh.")

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
               "  stb_bands --file siesta.bands --shift fermi\n"
               "  stb_bands --file siesta.bands --shift manual --manual-value 0.5",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("--file",  dest="input_file", type=str, required=True,
                        help="Path to the input file containing band structure data (e.g., siesta.bands).")

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

    parser.add_argument("-v", "--version", action="version", version=f"stb-bands {VERSION}")
    
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()


    # verify the manual value
    if args.shift == "manual" and args.manual_value is None:
        parser.error("--manual-value is required when --shift is set to 'manual'.")
    
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
    print("\n[INFO] Read file ...")
    fermi_energy,high_sym,dic_bands,nspin=read_data(args.input_file)
    print("[INFO] Calculate VBM and CBM ...")
    result = cbm_vbm(fermi_energy, high_sym, dic_bands, nspin, args.gap_tol)

    mesh_result = None
    if args.eig_file:
        print("[INFO] Read --eig-file (k-mesh) ...")
        fermi_mesh, dic_mesh, nspin_mesh = read_eig_mesh(args.eig_file)
        if abs(fermi_mesh - fermi_energy) > 1e-3:
            print(f"[WARNING] Fermi energy mismatch between --file ({fermi_energy:.6f} eV) and "
                  f"--eig-file ({fermi_mesh:.6f} eV); the two inputs may be from different calculations.")
        if nspin_mesh != nspin:
            print(f"[WARNING] nspin mismatch between --file (nspin={nspin}) and --eig-file "
                  f"(nspin={nspin_mesh}); mesh comparison uses combined values only.")
        print("[INFO] Calculate mesh (k-grid) VBM and CBM ...")
        mesh_result = cbm_vbm(fermi_mesh, [], dic_mesh, nspin_mesh, args.gap_tol)

    write_analysis_report(fermi_energy, result, nspin, mesh_result, args.gap_tol)
    print("[INFO] Output saved to bands_analysis.txt")
    vbm, cbm, vbm_k, cbm_k, band_gap, gap_type = result["combined"]

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
