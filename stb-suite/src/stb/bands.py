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

def _is_gamma(label):
    clean_text = re.sub(r'[^a-zA-Z\s]', '', label).lower()
    return "gamma" in clean_text.split()

def plot_gnuplot(high_sym):
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
    fileout.append('#set style line 2 lc rgb "#ff7f0e" lt 2 lw 2 pt 5 ps 1.5   # Dashed orange line\n')
    fileout.append('#set style line 3 lc rgb "#2ca02c" lt 3 lw 2               # Solid green line\n')
    fileout.append('\n')
    fileout.append('# Plot commands\n')
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
        dic_bands[key] = np.array(values, dtype=float)

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
        label = parts[1].strip("'\"") if len(parts) > 1 else ""
        high_sym.append([parts[0], label])
        idx += 1

    return fermi_energy, high_sym, dic_bands

def write_gnuplot_bands(dic_bands):
    # define initial key
    key_init = next(iter(dic_bands))
    # write the file
    file_name="bands_gnuplot.dat"
    with open(file_name, 'w') as file:
        for i in range(len(dic_bands[key_init])):
            for key in dic_bands.keys():
                file.write(f"{key}     {dic_bands[key][i]}\n")
            file.write("\n")
    return

def cbm_vbm(fermi_energy,high_sym,dic_bands):
    # Start value as infinity
    vbm = -np.inf
    cbm = np.inf
    vbm_k = None
    cbm_k = None
    for k, band in dic_bands.items():
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
    if band_gap == 0.0:
        gap_type = "Metallic"
    elif vbm_k == cbm_k:
        gap_type = "Direct"
    else:
        gap_type = "Indirect"
    print(f"[INFO] Fermi: {fermi_energy} \n[INFO] VBM: {vbm:.6f} (k={vbm_k}) \n[INFO] CBM: {cbm:.6f} (k={cbm_k})\n[INFO] Band Gap: {band_gap:.6f} eV ({gap_type})")
    return vbm, cbm, vbm_k, cbm_k, band_gap, gap_type

def write_analysis_report(fermi_energy, vbm, cbm, vbm_k, cbm_k, band_gap, gap_type):
    lines = [
        "BAND STRUCTURE ANALYSIS REPORT - STB Suite",
        "-"*60,
        f"Fermi energy : {fermi_energy:.6f} eV",
        f"VBM          : {vbm:.6f} eV (k = {vbm_k})",
        f"CBM          : {cbm:.6f} eV (k = {cbm_k})",
        f"Band gap     : {band_gap:.6f} eV",
        f"Gap type     : {gap_type}",
    ]
    content = "\n".join(lines) + "\n"
    with open("bands_analysis.txt", "w") as f:
        f.write(content)
    return content

def shift_bands(dic, val):
    return {k: [v - val for v in list] for k, list in dic.items()}

def plot(dic,custom_ticks):
    # Organize the data
    x_values = sorted(dic.keys())
    num_lines = len(next(iter(dic.values())))
    y_series = [[] for _ in range(num_lines)]
    for x in x_values:
        for i in range(num_lines):
            y_series[i].append(dic[x][i])
    # Organize the high symmetries points
    tick_positions = [float(t[0]) for t in custom_ticks]
    tick_labels = ["Γ" if _is_gamma(t[1]) else t[1] for t in custom_ticks]
    # plot
    plt.figure(figsize=(8, 6))
    for i, y_vals in enumerate(y_series):
        plt.xticks(tick_positions, tick_labels)
        plt.plot(x_values, y_vals,color='blue')
    # vertical lines in High Symmetries
    for pos in tick_positions:
        plt.axvline(x=pos, color='gray', linestyle='--', linewidth=1)
    # plot limits
    plt.ylim(-20, 20)
    plt.ylabel("Energy")
    plt.grid(True)
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
    fermi_energy,high_sym,dic_bands=read_data(args.input_file)
    print("[INFO] Calculate VBM and CBM ...")
    vbm, cbm, vbm_k, cbm_k, band_gap, gap_type = cbm_vbm(fermi_energy, high_sym, dic_bands)
    write_analysis_report(fermi_energy, vbm, cbm, vbm_k, cbm_k, band_gap, gap_type)
    print("[INFO] Output saved to bands_analysis.txt")

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
    write_gnuplot_bands(shifted_bands)
    plot(shifted_bands, high_sym)
    plot_gnuplot(high_sym)

    print("\n[INFO] Complete job!") 
    print("\n"+"-"*60)
    print(color_text("Bands found! But still no sign of Metallica.\n\n", 'bold'))

if __name__ == "__main__":
    main()
