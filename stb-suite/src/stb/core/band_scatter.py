"""Shared "scatter a scalar quantity onto a band structure" helpers: one
point per (k, band[, category]), sized/colored by a continuous value.

Extracted from fatbands.py once stb-ipr became a second consumer of the
exact same shape of output (fatbands: weight per orbital/atom/species
category; stb-ipr: a single continuous IPR value per (k, band), no
category at all -- exactly fatbands.py's own "single category" case,
generalized so neither tool re-implements it).
"""

import os

import numpy as np
import matplotlib.pyplot as plt


# Distinct, colorblind-legible colors, cycled if there are more series
# than colors (a plain fallback, not meant to be exhaustive).
PLOT_COLORS = [
    "#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e",
    "#17becf", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22",
]


def write_scalar_data(output_dir, name, rows):
    """rows: list of (k_position, energy, value). One file per `name`,
    long format -- a splot/pm3d grid doesn't fit this data (discrete per
    (k, band) points, not a continuous (k, energy) surface)."""
    path = os.path.join(output_dir, f"{name}.dat")
    with open(path, "w") as f:
        f.write(f"#{'k_position':<14}{'energy(eV)':<16}{'value':<10}\n")
        for k_pos, energy, value in rows:
            f.write(f"{k_pos:<15.6f}{energy:<16.6f}{value:<10.6f}\n")
    return path


def write_scalar_gplot(output_dir, gplot_name, pdf_name, high_sym, names, ps_scale=8.0):
    """One "plot ... with points" term per entry in `names` (each backed
    by its own `write_scalar_data`-written &lt;name&gt;.dat), sharing one set
    of xtics/high-symmetry vertical lines. `_is_gamma` substitution is
    inlined here rather than imported from core.siesta_bands, to avoid a
    dependency cycle since this module is meant to be usable by any
    band-like scatter tool, not just ones that already depend on
    siesta_bands."""
    def is_gamma(label):
        import re
        clean = re.sub(r'[^a-zA-Z\s]', '', label).lower()
        return "gamma" in clean.split()

    gplot_labels = [
        "{/Symbol G}" if is_gamma(label) else label for _, label in high_sym
    ]
    lines = []
    lines.append('set terminal pdfcairo enhanced font "Arial,25" size 8,8\n')
    lines.append(f'set output "{pdf_name}"\n')
    lines.append('set ylabel "Energy (eV)" font "Arial,28"\n')
    lines.append('set xlabel "" font "Arial,28"\n')
    xticsl = "set xtics ("
    for i in range(len(high_sym)):
        xticsl += f' "{gplot_labels[i]}" {float(high_sym[i][0])} , '
    xticsl = xticsl[:-2] + ')  font "Arial,28"\n'
    lines.append(xticsl)
    lines.append('set ytics format "%.1f" font "Arial,28"\n')
    lines.append('set grid xtics front\n')
    lines.append(f'set arrow from {float(high_sym[0][0])},0.0 to {float(high_sym[-1][0])} ,0.0 nohead dt 2 lc rgb "dark-gray" lw 4 back\n')
    for i in range(len(high_sym) - 2):
        lines.append(f'set arrow from {float(high_sym[i+1][0])} ,graph 0 to {float(high_sym[i+1][0])},graph 1 nohead dt 2 lc rgb "dark-gray" lw 4 back\n')

    plot_terms = []
    for i, name in enumerate(names):
        color = PLOT_COLORS[i % len(PLOT_COLORS)]
        plot_terms.append(
            f'"{name}.dat" using 1:2:({ps_scale}*$3) with points '
            f'pt 7 ps variable lc rgb "{color}" title "{name}"'
        )
    lines.append("plot " + ", \\\n     ".join(plot_terms) + "\n")

    path = os.path.join(output_dir, gplot_name)
    with open(path, "w") as f:
        f.writelines(lines)
    return path


def plot_scalar_on_bands(high_sym, k_pos, energy, value, value_label, is_gamma_fn, is_signed=False):
    """Single continuous scalar (e.g. IPR, fatbands' single-category case,
    or a signed quantity like a spin-texture component) scattered onto the
    band structure: marker size and color both encode `value` (same
    redundant-encoding convention as stb-fatbands/stb-stm), with a
    colorbar. `is_gamma_fn` is passed in (rather than imported) so this
    module has no dependency on core.siesta_bands.

    `is_signed` matters for two things: marker SIZE always uses abs(value)
    (a negative size crashes matplotlib -- verified live: IPR/fatbands
    weights are never negative so this never came up before a signed
    quantity, e.g. spin_moment's Sx/Sy/Sz in [-1, 1], was plotted this
    way), and when True the color uses a diverging colormap centered at
    zero instead of a sequential one, so sign is visually legible."""
    tick_positions = [float(t[0]) for t in high_sym]
    tick_labels = ["Γ" if is_gamma_fn(t[1]) else t[1] for t in high_sym]

    plt.figure(figsize=(8, 6))
    plt.xticks(tick_positions, tick_labels)
    for pos in tick_positions:
        plt.axvline(x=pos, color='gray', linestyle='--', linewidth=1)

    value = np.asarray(value)
    size = 10 + 200 * np.abs(value)
    if is_signed:
        vmax = float(np.max(np.abs(value))) if value.size else 1.0
        vmax = vmax if vmax > 0 else 1.0
        sc = plt.scatter(k_pos, energy, s=size, c=value, cmap='coolwarm',
                          vmin=-vmax, vmax=vmax, alpha=0.85)
    else:
        sc = plt.scatter(k_pos, energy, s=size, c=value, cmap='viridis', alpha=0.85)
    plt.colorbar(sc, label=value_label)

    plt.ylabel("Energy (eV)")
    plt.grid(True)
    plt.show()


def plot_multi_series_on_bands(high_sym, series, is_gamma_fn):
    """Multiple named series (fatbands' multi-category case): one
    scatter per series, discrete color, shared legend. `series` is a dict
    name -> list of (k_pos, energy, value)."""
    tick_positions = [float(t[0]) for t in high_sym]
    tick_labels = ["Γ" if is_gamma_fn(t[1]) else t[1] for t in high_sym]

    plt.figure(figsize=(8, 6))
    plt.xticks(tick_positions, tick_labels)
    for pos in tick_positions:
        plt.axvline(x=pos, color='gray', linestyle='--', linewidth=1)

    for i, (name, rows) in enumerate(series.items()):
        if not rows:
            continue
        k_pos, energy, value = zip(*rows)
        value = np.asarray(value)
        color = PLOT_COLORS[i % len(PLOT_COLORS)]
        # abs() for the same reason plot_scalar_on_bands uses it: a
        # negative marker size crashes matplotlib. Dormant for
        # stb-fatbands' current values (orbital weights, always >= 0),
        # but a latent crash for any future signed-quantity consumer.
        plt.scatter(k_pos, energy, s=10 + 200 * np.abs(value), c=color,
                    alpha=0.6, label=name)
    plt.legend()

    plt.ylabel("Energy (eV)")
    plt.grid(True)
    plt.show()
