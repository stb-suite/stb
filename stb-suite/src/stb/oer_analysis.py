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
from pymatgen.core.periodic_table import Element
from stb.core import structure_io
from stb.core.cli import color_text, show_intro, print_dual, print_section, print_table, get_input
from stb.core.siesta_log import get_free_energy, check_scf_and_force, report_quality_diagnostics
from stb.core.phonon_workflow import load_phonon_with_force_constants

REPORT_FILE = "oer_stage4.txt"

# Literature gas-phase H2 reference constant (Norskov et al., widely used throughout the
# CHE/HER/OER literature) -- reused verbatim from HER's her_analysis.py. Fixed at its 298.15 K
# literature value regardless of --temp, same simplification HER already makes: H2's own
# vibrational structure is experimentally very well known, unlike the adsorbed/isolated
# intermediates this module computes from the user's own DFT+phonon data.
H2_ZPE_EV = 0.270
H2_TS_298K_EV = 0.400

# 4 * 1.23 V, the well-established experimental reversible potential for 2H2O -> O2 + 4(H+ + e-)
# (Rossmeisl et al. 2007; Man et al. 2011) -- NOT flagged [UNVERIFIED], unlike a per-species ZPE/TS
# table would be (see the OER plan's Decision 4: no --zpe-mode standard in v1 for exactly this
# reason). G(O2) is DERIVED from this value rather than computed via DFT -- see main()'s
# description for why (Sargeant et al. 2021: >=0.3 eV systematic GGA/O2-triplet error).
FOUR_TIMES_U0_EV = 4.920

_FREQ_CONVERSION_THZ = 15.633302  # sqrt(eV / (amu * Ang^2)) -> THz, same constant HER/ASE use
_BOLTZMANN_EV_K = 8.617333262e-5  # eV/K
_EV_PER_THZ = 0.00413566733  # E = h*f, h in eV.s, f in THz*1e12 -> eV
_KJ_MOL_TO_EV = 0.01036427
_J_MOL_TO_EV = _KJ_MOL_TO_EV / 1000.0


def find_winning_dir(candidates_root, out_file, prefix="site_"):
    """Re-derives which folder (name starting with `prefix`) under
    candidates_root has the lowest FreeEng -- duplicated logic (not
    imported) from oer_refs.py's own find_winning_site, same "each stage
    re-derives, doesn't trust a sibling stage's report" policy already
    established for HER (e.g. her_analysis.py's own
    find_winning_site_energy). Also reused by locate_intermediate_dir below
    (prefix='ooh_star_orient') to pick the winner among OOH* orientation
    candidates sampled AT THE SAME SITE. Returns (basename, dir, energy) or
    (None, None, None).
    """
    if not os.path.isdir(candidates_root):
        return None, None, None
    site_dirs = sorted(
        d.path for d in os.scandir(candidates_root) if d.is_dir() and d.name.startswith(prefix))
    best_dir, best_energy = None, float('inf')
    for d in site_dirs:
        energy = get_free_energy(os.path.join(d, out_file))
        if energy is not None and energy < best_energy:
            best_energy = energy
            best_dir = d
    if best_dir is None:
        return None, None, None
    return os.path.basename(best_dir), best_dir, best_energy


def locate_intermediate_dir(output_root, name, out_file):
    """Returns the final relaxed folder for O*/OOH* ('name' is 'o' or
    'ooh'), matching oer_refs.py's own _locate_intermediate resolution
    (duplicated, not imported): the single 'intermediates/<name>_star/'
    folder if it exists, else (OOH* orientation sampling only) the winner
    among 'intermediates/<name>_star_orient*/' -- MULTIPLE ORIENTATIONS AT
    THE SAME SITE, never a different one. Either way, stb-oerIntermediates
    (Stage 2) derives O*/OOH* from the SAME winning OH* site
    unconditionally.
    """
    single_dir = os.path.join(output_root, "intermediates", f"{name}_star")
    if os.path.isdir(single_dir):
        return single_dir
    intermediates_root = os.path.join(output_root, "intermediates")
    _label, d, _energy = find_winning_dir(intermediates_root, out_file, prefix=f"{name}_star_orient")
    return d


def read_fa_force(fa_path, atom_index):
    """Reads the force (Fx, Fy, Fz, eV/Ang) on ONE atom (0-based
    atom_index) from a SIESTA .FA file -- format: first line = atom
    count, then one row per atom '<1-based index> Fx Fy Fz'. Same
    fail-soft "None on any parse issue" contract as core.siesta_log's own
    parsers. Duplicated from her_analysis.py.
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
    """Generalizes her_analysis.py::compute_local_zpe_entropy from ONE
    atom (3x3 Hessian, uniform mass) to N local atoms (3N x 3N Hessian):
    reads zpe_local_meta.json (local_indices/local_symbols/
    displacement_ang/order/system_label) and, for every displacement
    folder, reads the force on EVERY local atom (not just the one
    displaced) to assemble the full coupled block
    Phi[3i+a, 3j+b] = -(F_j,b(+) - F_j,b(-)) / (2*d) -- captures
    inter-atom coupling (e.g. O-H stretching within OH*/OOH*/H2O), not
    independent per-atom Hessians. Symmetrizes Phi, then builds the
    MASS-WEIGHTED dynamical matrix D[3i+a,3j+b] = Phi[3i+a,3j+b] /
    sqrt(m_i * m_j) (masses via pymatgen.core.periodic_table.Element)
    BEFORE diagonalizing -- the correct generalization of HER's
    per-eigenvalue '/mass' division (only valid there because all 3 DOF
    belonged to one atom of one mass); reduces EXACTLY to HER's formula
    when there's only 1 local atom (sqrt(m*m) == m). Eigenvalues of D are
    already omega^2, no further per-mode mass division needed.

    Returns (zpe_ev, ts_ev) for this ONE species/intermediate -- local
    mode's substrate-frozen (or, for H2O, isolated-molecule) treatment
    needs no further "delta vs. clean" subtraction, exactly as in HER.
    Returns (None, None) on any read failure.
    """
    meta_path = os.path.join(zpe_dir, "zpe_local_meta.json")
    try:
        with open(meta_path) as f:
            meta = json.load(f)
    except (OSError, ValueError):
        print_dual(color_text(f"    [ERROR] Could not read '{meta_path}'.", 'red'), f_out)
        return None, None

    local_indices = meta["local_indices"]
    local_symbols = meta["local_symbols"]
    displacement_ang = meta["displacement_ang"]
    order = meta["order"]
    system_label = meta["system_label"]
    n_local = len(local_indices)
    index_to_local = {atom_index: k for k, atom_index in enumerate(local_indices)}
    # local_symbols is written as the BARE real element by
    # oer_refs.py::write_local_zpe_folders (never a Stage-1/2 fragment
    # label like 'O_ads', which Element() cannot parse) -- resolve any
    # older zpe_local_meta.json written before that fix by falling back
    # to structure_io.real_element against the first displacement
    # folder's own structure.fdf species_meta.
    try:
        masses = np.array([Element(sym).atomic_mass for sym in local_symbols], dtype=float)
    except ValueError:
        disp1_structure = structure_io.read_fdf(os.path.join(zpe_dir, "disp_001", "structure.fdf"))
        masses = np.array([Element(structure_io.real_element(sym, disp1_structure.species_meta)).atomic_mass
                            for sym in local_symbols], dtype=float)

    forces = {}  # (moved_local_k, axis, sign) -> (n_local, 3) array of forces on every local atom
    for i, entry in enumerate(order, start=1):
        fa_path = os.path.join(zpe_dir, f"disp_{i:03d}", f"{system_label}.FA")
        moved_k = index_to_local[entry["atom_index"]]
        row = np.zeros((n_local, 3))
        for j, atom_index in enumerate(local_indices):
            force = read_fa_force(fa_path, atom_index)
            if force is None:
                print_dual(color_text(f"    [ERROR] Could not read forces from '{fa_path}'.", 'red'), f_out)
                return None, None
            row[j] = force
        forces[(moved_k, entry["axis"], entry["sign"])] = row

    dof = 3 * n_local
    Phi = np.zeros((dof, dof))
    for moved_k in range(n_local):
        for axis in range(3):
            f_plus = forces[(moved_k, axis, 1.0)]
            f_minus = forces[(moved_k, axis, -1.0)]
            d_force = -(f_plus - f_minus) / (2.0 * displacement_ang)  # (n_local, 3)
            col = 3 * moved_k + axis
            Phi[:, col] = d_force.reshape(-1)
    Phi = 0.5 * (Phi + Phi.T)

    mass_vec = np.repeat(masses, 3)
    D = Phi / np.sqrt(np.outer(mass_vec, mass_vec))

    eigenvalues = np.linalg.eigvalsh(D)
    zpe_ev, ts_ev = 0.0, 0.0
    print_dual(f"    Local (partial-Hessian, {n_local}-atom) vibrational modes:", f_out)
    for eig in eigenvalues:
        if eig > 0:
            freq_thz = _FREQ_CONVERSION_THZ * np.sqrt(eig)
            e_mode = freq_thz * _EV_PER_THZ
            print_dual(f"      {freq_thz * 33.35641:>8.2f} cm^-1  ({freq_thz:>6.2f} THz)  "
                        f"E = {e_mode:.4f} eV", f_out)
            zpe_ev += 0.5 * e_mode
            if temperature_k > 0:
                x = e_mode / (_BOLTZMANN_EV_K * temperature_k)
                S = _BOLTZMANN_EV_K * (x / np.expm1(x) - np.log1p(-np.exp(-x)))
                ts_ev += temperature_k * S
        else:
            freq_thz = _FREQ_CONVERSION_THZ * np.sqrt(-eig)
            print_dual(color_text(
                f"      {freq_thz:>6.2f} THz (IMAGINARY) -- excluded from ZPE/entropy", 'yellow'), f_out)
    return zpe_ev, ts_ev


def compute_full_zpe_entropy(zpe_dir, system_label, temperature_k, f_out):
    """Full-structure ZPE/entropy via phonopy's own thermal-properties
    pipeline, reusing core.phonon_workflow.load_phonon_with_force_constants
    -- identical to her_analysis.py's own compute_full_zpe_entropy (zero
    changes needed: it already operates on a whole structure's phonon
    calc regardless of atom count). Returns (zpe_ev, ts_ev) or
    (None, None) on failure.
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


# The 5 states of the 4-electron AEM reaction coordinate, in order -- shared
# by the [6] cumulative free-energy table and both plot functions below, so
# the printed numbers and the chart can never disagree.
OER_STATE_LABELS = ["* + H2O", "OH*", "O*", "OOH*", "* + O2"]


def write_oer_free_energy_diagram_plot(plot_dir, g_u0, g_ueq, eta, pds_index):
    """Writes a gnuplot .dat + .gplot free-energy (Norskov-style CHE) diagram:
    cumulative G along the 5-state reaction coordinate (* + H2O -> OH* -> O*
    -> OOH* -> * + O2), one series at U = 0 V and one at the equilibrium
    potential U = 1.23 V (each step there shifts by -eU, so the n-th state's
    G shifts by -n*1.23 eV -- see main()'s own g_ueq construction, which this
    function only plots, never recomputes). At U = 1.23 V the OVERALL
    reaction is flat (G returns to 0 at the last state, by construction of
    G(O2) -- see FOUR_TIMES_U0_EV), and the single steepest uphill segment in
    that curve IS the potential-determining step: `pds_index` (0-based, same
    index main() uses for 'Step N') is shaded to make that visually obvious
    rather than requiring the reader to eyeball slopes.

    Same explicit-x-position / lc rgb "#RRGGBB" (fixed per-series color, not
    the packed-0xRRGGBB-per-row 'lc rgb variable' trick her_analysis.py's
    categorical bar chart needs -- there are only 2 fixed-color line series
    here, not one bar per row) convention as the rest of the suite. Returns
    (dat_path, gplot_path).
    """
    os.makedirs(plot_dir, exist_ok=True)
    dat_path = os.path.join(plot_dir, "oer_free_energy_diagram.dat")
    x_states = [i + 1 for i in range(len(OER_STATE_LABELS))]  # 1-based, matches xtic column below
    hw = _PLATEAU_HALF_WIDTH

    def _plateau_block(values):
        # One line PER STATE (2 points, its own plateau), single blank lines
        # between them -- a lone blank line breaks a gnuplot 'with lines'
        # trace without starting a new index, which is exactly the "one
        # disconnected shelf per state" look this needs.
        lines = []
        for xi, yi in zip(x_states, values):
            lines.append(f"{xi - hw:.6f}  {yi:.6f}\n")
            lines.append(f"{xi + hw:.6f}  {yi:.6f}\n")
            lines.append("\n")
        return lines

    def _connector_block(values):
        lines = []
        for i in range(len(x_states) - 1):
            lines.append(f"{x_states[i] + hw:.6f}  {values[i]:.6f}\n")
            lines.append(f"{x_states[i + 1] - hw:.6f}  {values[i + 1]:.6f}\n")
            lines.append("\n")
        return lines

    with open(dat_path, "w") as f:
        f.write("# OER free-energy (CHE) diagram -- step/plateau form (see module docstring:\n")
        f.write("# each state is a stable G LEVEL, not a point on a slanted line). 4 blocks\n")
        f.write("# (gnuplot 'index' 0-3), each a series of 2-point line segments (breaks on\n")
        f.write("# blank lines): 0=U0 plateaus, 1=U0 dashed connectors, 2=Ueq plateaus,\n")
        f.write("# 3=Ueq dashed connectors. 1:x(Ang, arbitrary reaction-coordinate units) 2:G(eV)\n\n")
        f.writelines(_plateau_block(g_u0))
        f.write("\n")
        f.writelines(_connector_block(g_u0))
        f.write("\n")
        f.writelines(_plateau_block(g_ueq))
        f.write("\n")
        f.writelines(_connector_block(g_ueq))

    values = list(g_u0) + list(g_ueq)
    y_min, y_max = min(values + [0.0]), max(values + [0.0])
    span = max(y_max - y_min, 1e-6)
    y_bottom = min(y_min, 0.0) - 0.30 * span
    y_top = max(y_max, 0.0) + 0.20 * span

    xtics = ", ".join(f'"{label}" {x}' for label, x in zip(OER_STATE_LABELS, x_states))

    gplot_path = os.path.join(plot_dir, "oer_free_energy_diagram.gplot")
    with open(gplot_path, "w") as f:
        f.writelines([
            '# --- STB Plot Configuration ---\n',
            '# Generated by stb-oerAnalysis\n',
            'set terminal pdfcairo enhanced color font "Arial,12" size 9,6\n',
            'set output "oer_free_energy_diagram.pdf"\n\n',
            f'set title "OER free-energy diagram (eta = {eta:+.3f} V, PDS = Step {pds_index + 1})"\n',
            'set xlabel "Reaction coordinate"\n',
            'set ylabel "Free energy G (eV, relative to * + H2O)"\n',
            f'set xtics ({xtics}) rotate by -15 right\n',
            f'set xrange [{x_states[0] - 0.6}:{x_states[-1] + 0.6}]\n',
            f'set yrange [{y_bottom:.6f}:{y_top:.6f}]\n',
            'set grid ytics lt 0 lw 1 lc rgb "#bbbbbb"\n',
            'set xzeroaxis lt -1 lw 1 lc rgb "#000000"\n',
            'set key top left\n',
            # Shade the PDS step (the segment between its two endpoint states)
            # so the steepest-uphill segment of the U=1.23V curve is visually
            # obvious, not just implied by the title text.
            f'set object 1 rect from {pds_index + 1},graph 0 to {pds_index + 2},graph 1 '
            'fc rgb "#ffe8a3" fs solid 0.35 noborder behind\n',
            f'set label 1 "PDS" at {pds_index + 1.5},graph 0.93 center tc rgb "#cc6600" '
            'font ",11" front\n\n',
            # Plateaus (solid, thick) + connectors (dashed, thin), 2 colors --
            # 4 separate 'index' blocks from the SAME datafile, no per-row
            # color column needed (each series is one fixed color, unlike
            # her_analysis.py's categorical bar chart).
            'plot "oer_free_energy_diagram.dat" index 0 using 1:2 with lines lw 4 '
            'lc rgb "#2255cc" title "U = 0 V", \\\n'
            '     "" index 1 using 1:2 with lines lw 1.5 dt 2 lc rgb "#2255cc" notitle, \\\n'
            '     "" index 2 using 1:2 with lines lw 4 lc rgb "#cc2222" '
            'title "U = 1.23 V (eq.)", \\\n'
            '     "" index 3 using 1:2 with lines lw 1.5 dt 2 lc rgb "#cc2222" notitle\n',
        ])
    return dat_path, gplot_path


_PLATEAU_HALF_WIDTH = 0.32  # each state's G value is a level, not a point in time -- see
                            # _build_oer_free_energy_figure's own docstring for why a plateau,
                            # not a line-through-a-point, is the physically correct way to draw it


def _build_oer_free_energy_figure(g_u0, g_ueq, eta, pds_index):
    """Shared matplotlib figure builder for the free-energy diagram -- used
    by BOTH show_oer_free_energy_diagram_matplotlib (plt.show(), needs a
    display) and write_oer_free_energy_diagram_png (fig.savefig(), headless
    -safe) so the on-screen chart and the PNG embedded in the .md report can
    never visually disagree. Returns the Figure (caller decides show/save).

    Drawn as a STEP/PLATEAU diagram (the standard convention in the CHE/OER
    literature, e.g. Nørskov-style volcano/reaction-coordinate plots): each
    state is a stable, well-defined free-energy LEVEL -- not an instant in
    time -- so it's drawn as a short horizontal "shelf" at that G value,
    not a single point on a slanted line. A dashed segment connects
    consecutive plateaus purely to guide the eye from one level to the
    next; the slope of that dashed segment carries no physical meaning
    (there's no intermediate state living partway between OH* and O*, for
    instance) -- solid vs. dashed is exactly how the literature marks that
    distinction visually.
    """
    x_states = np.arange(len(OER_STATE_LABELS))
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.axvspan(pds_index, pds_index + 1, color="#ffe8a3", alpha=0.5, zorder=0)
    # Anchored near the BOTTOM of the axes (blended transform: x in data
    # coordinates, y in axes-fraction) rather than near the top of the data
    # range -- the legend always sits top-left, and the PDS step can be
    # anywhere along the x-axis (step 1 here, but any of 1-4 in general),
    # so a top-anchored label would collide with the legend whenever the
    # PDS happens to fall under it. The bottom is never contested.
    ax.text(pds_index + 0.5, 0.06, "PDS", color="#cc6600", ha='center', va='bottom',
            fontsize=11, fontweight='bold', transform=ax.get_xaxis_transform())

    for values, color, label in [(g_u0, "#2255cc", "U = 0 V"),
                                  (g_ueq, "#cc2222", "U = 1.23 V (eq.)")]:
        for xi, yi in zip(x_states, values):
            ax.plot([xi - _PLATEAU_HALF_WIDTH, xi + _PLATEAU_HALF_WIDTH], [yi, yi],
                     color=color, linewidth=3.5, solid_capstyle='butt', zorder=3)
        for i in range(len(x_states) - 1):
            ax.plot([x_states[i] + _PLATEAU_HALF_WIDTH, x_states[i + 1] - _PLATEAU_HALF_WIDTH],
                     [values[i], values[i + 1]], color=color, linewidth=1.2, linestyle='--',
                     zorder=2)
        ax.plot([], [], color=color, linewidth=3.5, label=label)  # legend-only proxy

    ax.axhline(0, color='black', linewidth=1)
    ax.set_xticks(x_states)
    ax.set_xticklabels(OER_STATE_LABELS, rotation=15, ha='right')
    ax.set_ylabel("Free energy G (eV, relative to * + H2O)")
    ax.set_title(f"OER free-energy diagram ($\\eta$ = {eta:+.3f} V, PDS = Step {pds_index + 1})")
    ax.grid(axis='y', linestyle=':', color='#bbbbbb')
    ax.legend(loc='upper left')
    fig.tight_layout()
    return fig


def show_oer_free_energy_diagram_matplotlib(g_u0, g_ueq, eta, pds_index):
    """Opens a blocking matplotlib window with the SAME free-energy diagram
    write_oer_free_energy_diagram_plot saves for gnuplot -- an on-screen
    convenience alongside (not instead of) the gnuplot .dat/.gplot pair,
    same precedent as her_analysis.py's show_energy_breakdown_matplotlib.
    Returns True if the window was shown, False if matplotlib couldn't open
    a display.
    """
    try:
        fig = _build_oer_free_energy_figure(g_u0, g_ueq, eta, pds_index)
        plt.show()
        return True
    except Exception as e:
        print(color_text(f"[NOTE] Could not open a matplotlib window ({e}) -- no display "
                          "available (e.g. a headless/SSH session without X11 forwarding). "
                          "Use --plot for the gnuplot .dat/.gplot files instead.", 'yellow'))
        return False


def write_oer_free_energy_diagram_png(plot_dir, g_u0, g_ueq, eta, pds_index):
    """Saves the same free-energy diagram as a PNG via fig.savefig -- unlike
    plt.show(), this never needs a live display, so it's safe to call
    unconditionally (used to embed a real chart image in
    write_markdown_report's consolidated .md deliverable, independent of
    whether the user asked for the interactive --show window). Returns the
    PNG path, or None if matplotlib failed to render for any reason (the
    caller still writes the rest of the .md report rather than crashing the
    whole Stage 4 run over a missing image).
    """
    try:
        os.makedirs(plot_dir, exist_ok=True)
        fig = _build_oer_free_energy_figure(g_u0, g_ueq, eta, pds_index)
        png_path = os.path.join(plot_dir, "oer_free_energy_diagram.png")
        fig.savefig(png_path, dpi=150)
        plt.close(fig)
        return png_path
    except Exception as e:
        print(color_text(f"[NOTE] Could not render the free-energy diagram PNG for the "
                          f"Markdown report ({e}).", 'yellow'))
        return None


# Colors for the BSSE-correction chart, same term/total distinction (and
# same rationale) as her_analysis.py's own _BREAKDOWN_TABLE_COLOR/
# _BREAKDOWN_PLOT_COLOR: 'term' = a plain additive contribution (the slab
# or adsorbate component), 'total' = the summed correction actually
# applied to that intermediate's energy in [2]'s table above.
_BSSE_PLOT_COLOR = {"term": "2255cc", "total": "22aa55"}


def write_bsse_correction_plot(plot_dir, bsse_chart_rows):
    """Writes a gnuplot .dat + .gplot categorical bar chart of the BSSE
    correction's slab/adsorbate components and total, per intermediate
    (OH*/O*/OOH*) -- same explicit-x, packed-0xRRGGBB-per-row 'lc rgb
    variable' convention as her_analysis.py's own energy-breakdown chart
    (see write_energy_breakdown_plot there for why explicit x, not
    histogram style). In --bsse-mode shared (the default), all 3
    intermediates' bars are identical by construction (one triad reused
    for all 3) -- the chart still shows all 9 bars so --bsse-mode full's
    genuinely different per-intermediate values render with the exact same
    code path. `bsse_chart_rows` is [(label, value_eV, kind), ...].
    Returns (dat_path, gplot_path).
    """
    os.makedirs(plot_dir, exist_ok=True)
    dat_path = os.path.join(plot_dir, "oer_bsse_correction.dat")
    with open(dat_path, "w") as f:
        f.write("# OER BSSE correction, per intermediate (slab/adsorbate components + total)\n")
        f.write("# 1:Index 2:Value(eV) 3:Label 4:Color(0xRRGGBB, term/total)\n")
        for i, (label, value, kind) in enumerate(bsse_chart_rows, start=1):
            f.write(f'{i}  {value:.6f}  "{label}"  0x{_BSSE_PLOT_COLOR[kind]}\n')

    values = [value for _label, value, _kind in bsse_chart_rows]
    y_min, y_max = min(values + [0.0]), max(values + [0.0])
    span = max(y_max - y_min, 1e-6)
    y_bottom = min(y_min, 0.0) - 0.35 * span
    y_top = max(y_max, 0.0) + 0.15 * span

    gplot_path = os.path.join(plot_dir, "oer_bsse_correction.gplot")
    with open(gplot_path, "w") as f:
        f.writelines([
            '# --- STB Plot Configuration ---\n',
            '# Generated by stb-oerAnalysis\n',
            'set terminal pdfcairo enhanced color font "Arial,12" size 9,6\n',
            'set output "oer_bsse_correction.pdf"\n\n',
            'set title "OER BSSE correction by intermediate (slab / adsorbate / total)"\n',
            'set ylabel "Energy (eV)"\n',
            'set style fill solid 0.7 border -1\n',
            'set boxwidth 0.6\n',
            'set xtics rotate by -30 right\n',
            f'set yrange [{y_bottom:.6f}:{y_top:.6f}]\n',
            'set grid ytics lt 0 lw 1 lc rgb "#bbbbbb"\n',
            'set xzeroaxis lt -1 lw 2 lc rgb "#000000"\n',
            'unset key\n\n',
            'plot "oer_bsse_correction.dat" using 1:2:4:xtic(3) with boxes lc rgb variable\n',
        ])
    return dat_path, gplot_path


def _build_bsse_correction_figure(bsse_chart_rows):
    """Shared matplotlib figure builder for the BSSE-correction chart --
    used by both show_bsse_correction_matplotlib (plt.show()) and
    write_bsse_correction_png (fig.savefig()), same pattern as
    _build_oer_free_energy_figure above.
    """
    labels = [label for label, _value, _kind in bsse_chart_rows]
    values = [value for _label, value, _kind in bsse_chart_rows]
    colors = [f"#{_BSSE_PLOT_COLOR[kind]}" for _label, _value, kind in bsse_chart_rows]
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.bar(labels, values, color=colors, edgecolor='black', linewidth=0.5)
    ax.axhline(0, color='black', linewidth=1)
    ax.set_ylabel("Energy (eV)")
    ax.set_title("OER BSSE correction by intermediate (slab / adsorbate / total)")
    ax.grid(axis='y', linestyle=':', color='#bbbbbb')
    plt.setp(ax.get_xticklabels(), rotation=30, ha='right')
    fig.tight_layout()
    return fig


def show_bsse_correction_matplotlib(bsse_chart_rows):
    """Opens a blocking matplotlib window with the SAME BSSE-correction
    chart write_bsse_correction_plot saves for gnuplot. Returns True if
    shown, False if matplotlib couldn't open a display.
    """
    try:
        _fig = _build_bsse_correction_figure(bsse_chart_rows)
        plt.show()
        return True
    except Exception as e:
        print(color_text(f"[NOTE] Could not open a matplotlib window ({e}) -- no display "
                          "available (e.g. a headless/SSH session without X11 forwarding). "
                          "Use --plot for the gnuplot .dat/.gplot files instead.", 'yellow'))
        return False


def write_bsse_correction_png(plot_dir, bsse_chart_rows):
    """Saves the BSSE-correction chart as a PNG via fig.savefig -- headless
    -safe, used to embed a real chart image in write_markdown_report's
    consolidated .md deliverable (same pattern as
    write_oer_free_energy_diagram_png). Returns the PNG path, or None on
    any rendering failure.
    """
    try:
        os.makedirs(plot_dir, exist_ok=True)
        fig = _build_bsse_correction_figure(bsse_chart_rows)
        png_path = os.path.join(plot_dir, "oer_bsse_correction.png")
        fig.savefig(png_path, dpi=150)
        plt.close(fig)
        return png_path
    except Exception as e:
        print(color_text(f"[NOTE] Could not render the BSSE-correction chart PNG for the "
                          f"Markdown report ({e}).", 'yellow'))
        return None


def _build_bsse_raw_corrected_figure(bsse_summary_rows):
    """Shared matplotlib figure builder for the raw-vs-BSSE-corrected
    comparison: `bsse_summary_rows` is [(label, e_raw, bsse, e_corrected),
    ...] for OH*/O*/OOH*. One small-multiple PANEL per intermediate, each
    with its OWN y-axis range -- NOT a single shared axis. Two reasons a
    shared axis can't work here: (1) BSSE itself (~0.1-0.5 eV) is 3-4
    orders of magnitude smaller than E(raw)/E(corrected) (~-3700 to
    -4100 eV), so a bar chart from 0 would render the correction as a
    literally invisible sliver on top of an enormous bar; (2) the three
    intermediates' OWN raw energies differ from each other by hundreds of
    eV (different atom counts/composition), so even comparing raw values
    ACROSS intermediates on one axis would dwarf the eV-scale BSSE shift
    WITHIN each one. Each panel instead shows exactly the two numbers that
    are meaningfully comparable to each other -- Raw and Corrected for
    THAT intermediate -- as two BARS. Each bar is drawn the normal
    matplotlib way (from 0 up/down to its actual value, e.g. to -3709 eV)
    but the axes' own `ylim` is set to the tight [min-pad, max+pad] window
    around just those two values -- matplotlib clips the bar to what's
    visible rather than removing it, so within that zoomed window each bar
    reads as an ordinary floating bar reaching the panel's floor, with its
    height directly showing the BSSE shift between Raw and Corrected (also
    labeled numerically in the gap between the two bars).
    """
    n = len(bsse_summary_rows)
    fig, axes = plt.subplots(1, n, figsize=(3.4 * n, 5.5), sharex=False)
    if n == 1:
        axes = [axes]
    for ax, (label, e_raw, bsse, e_corr) in zip(axes, bsse_summary_rows):
        ax.bar([0], [e_raw], width=0.6, color="#2255cc", edgecolor='black', linewidth=0.5,
               zorder=2, label="Raw")
        ax.bar([1], [e_corr], width=0.6, color="#22aa55", edgecolor='black', linewidth=0.5,
               zorder=2, label="Corrected")
        mid_y = (e_raw + e_corr) / 2.0
        ax.annotate(f"BSSE\n{bsse:+.4f} eV", xy=(0.5, mid_y), ha='center', va='center',
                    fontsize=9, zorder=3,
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#888888"))
        ax.set_xlim(-0.6, 1.6)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Raw", "Corrected"])
        ax.set_title(label)
        span = max(abs(bsse), 1e-6)
        pad = span * 0.9
        ax.set_ylim(min(e_raw, e_corr) - pad, max(e_raw, e_corr) + pad)
        ax.grid(axis='y', linestyle=':', color='#dddddd', zorder=0)
        ax.set_axisbelow(True)
        ax.ticklabel_format(axis='y', useOffset=False, style='plain')
    axes[0].set_ylabel("Energy (eV)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.suptitle("OER BSSE correction: raw vs. corrected energy, per intermediate", y=1.08)
    fig.tight_layout()
    return fig


def show_bsse_raw_corrected_matplotlib(bsse_summary_rows):
    """Opens a blocking matplotlib window with the SAME raw-vs-corrected
    chart write_bsse_raw_corrected_plot saves for gnuplot. Returns True if
    shown, False if matplotlib couldn't open a display.
    """
    try:
        _fig = _build_bsse_raw_corrected_figure(bsse_summary_rows)
        plt.show()
        return True
    except Exception as e:
        print(color_text(f"[NOTE] Could not open a matplotlib window ({e}) -- no display "
                          "available (e.g. a headless/SSH session without X11 forwarding). "
                          "Use --plot for the gnuplot .dat/.gplot files instead.", 'yellow'))
        return False


def write_bsse_raw_corrected_png(plot_dir, bsse_summary_rows):
    """Saves the raw-vs-corrected chart as a PNG via fig.savefig --
    headless-safe, used to embed a real chart image in
    write_markdown_report's consolidated .md deliverable. Returns the PNG
    path, or None on any rendering failure.
    """
    try:
        os.makedirs(plot_dir, exist_ok=True)
        fig = _build_bsse_raw_corrected_figure(bsse_summary_rows)
        png_path = os.path.join(plot_dir, "oer_bsse_raw_corrected.png")
        fig.savefig(png_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        return png_path
    except Exception as e:
        print(color_text(f"[NOTE] Could not render the raw-vs-corrected chart PNG for the "
                          f"Markdown report ({e}).", 'yellow'))
        return None


def write_bsse_raw_corrected_plot(plot_dir, bsse_summary_rows):
    """Writes a gnuplot .dat + .gplot version of the raw-vs-corrected
    comparison (see _build_bsse_raw_corrected_figure's docstring for why
    each intermediate needs its OWN y-axis range, not one shared axis, and
    for the "ordinary boxes clipped to a tight yrange reads as a floating
    bar" trick used here too): one multiplot PANEL per intermediate, each
    with its own `set yrange` (pre-computed here in Python, same "compute
    the padding once, write it as a literal" convention as every other
    chart in this module) tightly zoomed around just that intermediate's
    Raw/Corrected pair, drawn as two `with boxes` bars, labeled with the
    BSSE value via `set label` in the gap between them.
    Returns (dat_path, gplot_path).
    """
    os.makedirs(plot_dir, exist_ok=True)
    dat_path = os.path.join(plot_dir, "oer_bsse_raw_corrected.dat")
    with open(dat_path, "w") as f:
        f.write("# OER BSSE correction: raw vs. corrected energy, per intermediate\n")
        f.write("# 1:Label 2:E_raw(eV) 3:BSSE(eV) 4:E_corrected(eV)\n")
        for label, e_raw, bsse, e_corr in bsse_summary_rows:
            f.write(f'"{label}"  {e_raw:.6f}  {bsse:.6f}  {e_corr:.6f}\n')

    n = len(bsse_summary_rows)
    gplot_path = os.path.join(plot_dir, "oer_bsse_raw_corrected.gplot")
    lines = [
        '# --- STB Plot Configuration ---\n',
        '# Generated by stb-oerAnalysis\n',
        'set terminal pdfcairo enhanced color font "Arial,11" size ' f'{3.4 * n:.1f},5.5\n',
        'set output "oer_bsse_raw_corrected.pdf"\n\n',
        f'set multiplot layout 1,{n} title '
        '"OER BSSE correction: raw vs. corrected energy, per intermediate"\n',
        'unset key\n',
        'set xrange [-0.6:1.6]\n',
        'set xtics ("Raw" 0, "Corrected" 1)\n',
        'set boxwidth 0.6\n',
        'set style fill solid 0.85 border -1\n',
        'set grid ytics lt 0 lw 1 lc rgb "#dddddd"\n\n',
    ]
    for i, (label, e_raw, bsse, e_corr) in enumerate(bsse_summary_rows):
        span = max(abs(bsse), 1e-6)
        pad = span * 0.9
        y_bottom = min(e_raw, e_corr) - pad
        y_top = max(e_raw, e_corr) + pad
        mid_y = (e_raw + e_corr) / 2.0
        lines.extend([
            'unset label\n',
            f'set label "BSSE\\n{bsse:+.4f} eV" at 0.5,{mid_y:.6f} center front '
            'boxed tc rgb "#444444" font ",10"\n',
            f'set yrange [{y_bottom:.6f}:{y_top:.6f}]\n',
            f'set title "{label}"\n',
            # Ordinary boxes (drawn from y=0 up/down to the real value, e.g.
            # -3709 eV) clipped by the tight yrange above -- reads as a
            # floating bar reaching the panel floor, same trick the
            # matplotlib version uses (see _build_bsse_raw_corrected_figure).
            f'plot "oer_bsse_raw_corrected.dat" every ::{i}::{i} using (0):2 with boxes '
            'lc rgb "#2255cc", \\\n'
            f'     "" every ::{i}::{i} using (1):4 with boxes lc rgb "#22aa55"\n\n',
        ])
    lines.append('unset multiplot\n')
    with open(gplot_path, "w") as f:
        f.writelines(lines)
    return dat_path, gplot_path


def markdown_table(headers, rows):
    """Plain GitHub-flavored-Markdown table (list of str rows already
    formatted, header row, '---' separator) -- the single formatting
    helper write_markdown_report reuses for every section table instead of
    hand-writing '|'-joins repeatedly.
    """
    lines = ["| " + " | ".join(headers) + " |",
              "|" + "|".join(["---"] * len(headers)) + "|"]
    for row in rows:
        lines.append("| " + " | ".join(str(c) for c in row) + " |")
    return "\n".join(lines)


def write_markdown_report(output_root, md_name, ctx):
    """Writes a single self-contained OER_report.md consolidating every
    number this Stage 4 run computed (electronic energies, BSSE, ZPE/TS,
    the 4 reaction steps, eta/PDS, the full free-energy diagram) -- 'todos
    os dados de todas as etapas' the user asked for, built from THIS run's
    own already-verified live data (not by re-parsing/concatenating the
    Stage 2/3 .txt reports, which may not even exist -- e.g. Stage 1 here
    was run without --save-report -- and which Stage 4 already fully
    re-derives itself, same "each stage re-derives, doesn't trust a sibling
    stage's report" policy as find_winning_dir/locate_intermediate_dir
    above). Embeds the free-energy diagram as a PNG (write_oer_free_energy_
    diagram_png, headless-safe) regardless of whether --plot/--show were
    used for the gnuplot/on-screen versions. Returns the .md path.
    """
    md_path = os.path.join(output_root, md_name)
    lines = []
    lines.append(f"# OER Study Report -- `{output_root}`")
    lines.append("")
    lines.append(f"*Generated {ctx['date']} by stb-oerAnalysis (STB-SUITE)*")
    lines.append("")
    lines.append("## 1. Overview")
    lines.append("")
    lines.append(markdown_table(
        ["Property", "Value"],
        [["ZPE mode", ctx['zpe_mode']],
         ["Temperature", f"{ctx['temp']} K"],
         ["Winning OH\\* site", ctx['winning_oh_label']],
         ["O\\* folder (derived)", f"`{ctx['o_dir']}`"],
         ["OOH\\* folder (derived)", f"`{ctx['ooh_dir']}`"]]))
    lines.append("")
    lines.append("## 2. Electronic Energies (FreeEng)")
    lines.append("")
    lines.append(markdown_table(
        ["Term", "Energy (eV)", "Source"],
        [[label, f"{energy:.6f}", f"`{path}`"] for label, energy, path in ctx['energy_rows']]))
    lines.append("")
    lines.append("## 3. BSSE Correction")
    lines.append("")
    if ctx.get('bsse_rc_png_path'):
        rel_png = os.path.relpath(ctx['bsse_rc_png_path'], output_root)
        lines.append(f"![OER BSSE correction: raw vs. corrected]({rel_png})")
        lines.append("")
    if ctx.get('bsse_png_path'):
        rel_png = os.path.relpath(ctx['bsse_png_path'], output_root)
        lines.append(f"![OER BSSE correction: slab/adsorbate components]({rel_png})")
        lines.append("")
    lines.append(markdown_table(
        ["Intermediate", "E (raw, eV)", "BSSE (eV)", "E (corrected, eV)"],
        [["OH*", f"{ctx['e_oh']:.4f}", f"{ctx['bsse_oh']:+.4f}", f"{ctx['e_oh_corr']:.4f}"],
         ["O*", f"{ctx['e_o']:.4f}", f"{ctx['bsse_o']:+.4f}", f"{ctx['e_o_corr']:.4f}"],
         ["OOH*", f"{ctx['e_ooh']:.4f}", f"{ctx['bsse_ooh']:+.4f}", f"{ctx['e_ooh_corr']:.4f}"]]))
    lines.append("")
    lines.append("## 4. Thermal Correction (ZPE / T·S)")
    lines.append("")
    lines.append(markdown_table(
        ["Species", "ZPE (eV)", "TS (eV)", "Delta-ZPE (eV)", "Delta-TS (eV)"],
        [[row[0], f"{row[1]:.4f}", f"{row[2]:.4f}",
          "--" if row[3] is None else f"{row[3]:+.4f}",
          "--" if row[4] is None else f"{row[4]:+.4f}"] for row in ctx['thermal_rows']]))
    lines.append("")
    lines.append("## 5. Gas-Phase O2 Reference (derived, not DFT-computed)")
    lines.append("")
    lines.append(f"G(H2) = {ctx['g_h2']:.4f} eV, G(H2O,l) = {ctx['g_h2o']:.4f} eV")
    lines.append("")
    lines.append(f"**G(O2) = 2·G(H2O) - 2·G(H2) + 4.92 = {ctx['g_o2']:.4f} eV** "
                  "(derived from the experimental 4×1.23 V reversible potential -- "
                  "see `--help` for why no O2 SIESTA reference is used)")
    lines.append("")
    lines.append("## 6. Reaction Steps & Potential-Determining Step (PDS)")
    lines.append("")
    lines.append(markdown_table(
        ["Step", "Reaction", "dG (eV)", "PDS"],
        [[str(i + 1), label.split(": ", 1)[1], f"{dg:+.4f}", "**<-- PDS**" if i == ctx['pds_index'] else ""]
         for i, (label, dg) in enumerate(ctx['steps'])]))
    lines.append("")
    lines.append(f"Sum dG1+dG2+dG3+dG4 = {sum(dg for _l, dg in ctx['steps']):.6f} eV "
                  f"(exactly {FOUR_TIMES_U0_EV} eV by construction)")
    lines.append("")
    lines.append("## 7. Free-Energy Diagram")
    lines.append("")
    if ctx.get('png_path'):
        rel_png = os.path.relpath(ctx['png_path'], output_root)
        lines.append(f"![OER free-energy diagram]({rel_png})")
        lines.append("")
    lines.append(markdown_table(
        ["State", "G at U=0 V (eV)", "G at U=1.23 V eq. (eV)"],
        [[label, f"{g0:+.4f}", f"{geq:+.4f}"]
         for label, g0, geq in zip(OER_STATE_LABELS, ctx['g_u0'], ctx['g_ueq'])]))
    lines.append("")
    lines.append("## 8. Final Result")
    lines.append("")
    lines.append(f"| **Theoretical overpotential eta** | **{ctx['eta']:+.4f} V** |")
    lines.append("|---|---|")
    lines.append(f"| **Potential-determining step** | **Step {ctx['pds_index'] + 1}** |")
    lines.append("")
    lines.append(f"Qualitative assessment: *{ctx['verdict']}*")
    lines.append("")
    lines.append("> **[LIMITATION]** AEM descriptor only (OH\\*/O\\*/OOH\\* at one site) -- does not "
                  "model the lattice oxygen evolution mechanism (LOER). See `stb-oerAnalysis --help`.")
    lines.append("")
    lines.append("## 9. Files")
    lines.append("")
    file_rows = [["Plain-text summary", f"`{ctx['report_str_path']}`"]]
    if ctx.get('report_path'):
        file_rows.append(["Full Stage 4 report (.txt)", f"`{ctx['report_path']}`"])
    if ctx.get('dat_path'):
        file_rows.append(["Gnuplot data/script (free energy)", f"`{ctx['dat_path']}`, `{ctx['gplot_path']}`"])
    if ctx.get('png_path'):
        file_rows.append(["Free-energy diagram (PNG)", f"`{ctx['png_path']}`"])
    if ctx.get('bsse_dat_path'):
        file_rows.append(["Gnuplot data/script (BSSE components)",
                           f"`{ctx['bsse_dat_path']}`, `{ctx['bsse_gplot_path']}`"])
    if ctx.get('bsse_png_path'):
        file_rows.append(["BSSE correction chart, components (PNG)", f"`{ctx['bsse_png_path']}`"])
    if ctx.get('bsse_rc_dat_path'):
        file_rows.append(["Gnuplot data/script (BSSE raw vs. corrected)",
                           f"`{ctx['bsse_rc_dat_path']}`, `{ctx['bsse_rc_gplot_path']}`"])
    if ctx.get('bsse_rc_png_path'):
        file_rows.append(["BSSE correction chart, raw vs. corrected (PNG)", f"`{ctx['bsse_rc_png_path']}`"])
    lines.append(markdown_table(["File", "Path"], file_rows))
    lines.append("")

    with open(md_path, "w") as f:
        f.write("\n".join(lines))
    return md_path


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Stage 4 of 4: combines every reference energy from stb-oerRefs "
        "into Delta-G1..4, the theoretical overpotential eta, and the potential-determining step "
        "(PDS).", 'bold')}
The 4-electron computational hydrogen electrode (CHE) descriptor for OER (Rossmeisl et al. 2007;
Man et al. 2011), via the 3 adsorbed intermediates OH*, O*, OOH*:
    * + H2O(l)  -> OH*  + (H+ + e-)      dG1 = G(OH*)  + 0.5*G(H2) - G(*)    - G(H2O)
    OH*         -> O*   + (H+ + e-)      dG2 = G(O*)   + 0.5*G(H2) - G(OH*)
    O* + H2O(l) -> OOH* + (H+ + e-)      dG3 = G(OOH*) + 0.5*G(H2) - G(O*)   - G(H2O)
    OOH*        -> *    + O2 + (H+ + e-) dG4 = G(*) + G(O2) - G(OOH*) + 0.5*G(H2)
eta = max(dG1..dG4) - 1.23 V; the step with the largest dG is the potential-determining step (PDS).

[IMPORTANT] No O2 SIESTA reference is ever used. DFT/GGA badly describes the O2 triplet ground
state (a systematic error of at least ~0.3 eV across common functionals -- Sargeant et al. 2021,
J. Electroanal. Chem.), so G(O2) is instead DERIVED from the experimental total reaction free
energy of the overall 4-electron reaction (4 * 1.23 V = 4.92 eV):
    G(O2) = 2*G(H2O) - 2*G(H2) + 4.92 eV
By construction dG1+dG2+dG3+dG4 == 4.92 eV exactly, independent of the actual electronic
energies -- this module prints that sum as a live internal sanity check.

No --zpe-mode 'standard' is offered: every intermediate's (and H2O's own) ZPE/entropy is computed
from your own DFT+phonon data (--zpe-mode local/full in stb-oerRefs), never a hardcoded per-species
literature default -- OER would need 4 additional constants beyond HER's single well-established
H2 value, and their exact digits could not be pinned to a verified primary source during this
tool's development.

[LIMITATION] This descriptor assumes the adsorbate evolution mechanism (AEM: OH*/O*/OOH* bound at
a single site) -- it does NOT model the lattice oxygen evolution mechanism (LOER, where lattice
oxygen atoms participate directly, common on some oxide/oxyhydroxide catalysts and linked to
catalyst instability under anodic potential; see Exner, ChemCatChem 2021). Treat eta here as an
AEM-only estimate, not the full mechanistic picture, especially for perovskites or Ru/Ir oxides.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage example:\n"
               "  %(prog)s --directory oer_study --temp 298.15\n"
    )

    parser.add_argument("-dir", "--directory", type=str, default="oer_study",
                         help="Root directory written by stb-oer/stb-oerIntermediates/stb-oerRefs "
                              "(default: oer_study).")
    parser.add_argument("--file", type=str, default="calc.out",
                         help="SIESTA output filename inside each folder (default: calc.out).")
    parser.add_argument("--temp", type=float, default=298.15,
                         help="Temperature in Kelvin, for the local/full ZPE modes' entropy term "
                              "(default: 298.15). Has no effect on H2's own thermal term (fixed at "
                              "its 298.15 K literature value regardless of --temp, see the module "
                              "docstring's H2_ZPE_EV/H2_TS_298K_EV).")
    parser.add_argument("--force-tolerance", type=float, default=0.05,
                         help="Residual atomic force in eV/Ang (default: 0.05) above which a "
                              "folder is flagged as possibly not relaxed/converged. Advisory only.")
    parser.add_argument("-o", "--output", type=str, default="OER_report",
                         help="Base filename (no extension) for the final report (default: OER_report).")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to <directory>/{REPORT_FILE}. Off by default.")
    parser.add_argument("--plot", dest="plot", action="store_true", default=None,
                         help="Save gnuplot .dat/.gplot charts (free-energy (CHE) diagram AND "
                              "BSSE-correction breakdown) to <directory>/plot/, skipping the "
                              "interactive prompt this tool otherwise asks at the end. See "
                              "--no-plot.")
    parser.add_argument("--no-plot", dest="plot", action="store_false",
                         help="Skip both charts entirely, skipping the interactive prompt (see "
                              "--plot). Use this for non-interactive/scripted runs.")
    parser.add_argument("--show", dest="show", action="store_true", default=None,
                         help="Open blocking matplotlib windows with the same two charts on "
                              "screen (free-energy diagram, then BSSE correction), skipping the "
                              "interactive prompt this tool otherwise asks. Independent of "
                              "--plot/--no-plot (the gnuplot files, if any, are still the "
                              "canonical saved output). Needs a display. See --no-show.")
    parser.add_argument("--no-show", dest="show", action="store_false",
                         help="Skip both on-screen matplotlib windows entirely, skipping the "
                              "interactive prompt (see --show). Use this for non-interactive/"
                              "scripted/headless runs.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-oerAnalysis {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("OER WORKFLOW -- STAGE 4: ANALYSIS", 'bold'))
    print("-" * 60)

    output_root = args.directory
    stage2_report = os.path.join(output_root, "oer_stage2.txt")
    stage3_report = os.path.join(output_root, "oer_stage3.txt")
    if not os.path.isfile(stage3_report):
        print(color_text(f"[ERROR] '{stage3_report}' not found -- run stb-oerRefs (Stage 3) first.", 'red'))
        sys.exit(1)
    if not os.path.isfile(stage2_report):
        print(color_text(f"[ERROR] '{stage2_report}' not found -- run stb-oerIntermediates "
                          "(Stage 2) first.", 'red'))
        sys.exit(1)

    zpe_mode = None
    with open(stage3_report) as f:
        for line in f:
            if line.startswith("ZPE mode"):
                zpe_mode = line.split(":", 1)[1].strip()
                break
    if zpe_mode is None:
        print(color_text("[ERROR] Could not recover ZPE mode from the Stage 3 report.", 'red'))
        sys.exit(1)

    sites_root = os.path.join(output_root, "sites")
    winning_oh_label, winning_oh_dir, _e = find_winning_dir(sites_root, args.file)
    o_dir = locate_intermediate_dir(output_root, "o", args.file)
    ooh_dir = locate_intermediate_dir(output_root, "ooh", args.file)
    if winning_oh_dir is None or o_dir is None or ooh_dir is None:
        print(color_text("[ERROR] Could not locate the winning OH*/O*/OOH* folder(s) -- did every "
                          "relaxation finish?", 'red'))
        sys.exit(1)

    energy_rows = []  # (label, energy_eV, source_path) -- the [1] recap table below

    def read_energy(folder, label, f_out):
        out_path = os.path.join(output_root, folder, args.file)
        energy = get_free_energy(out_path)
        if energy is None:
            print_dual(color_text(f"[ERROR] Could not read energy from '{out_path}'.", 'red'), f_out)
            sys.exit(1)
        print_dual(f"  {label:<14} = {energy:>12.6f} eV  ({out_path})", f_out)
        report_quality_diagnostics(label, out_path, args.force_tolerance, f_out)
        energy_rows.append((label, energy, out_path))
        return energy

    def read_energy_at(path, label, f_out):
        out_path = os.path.join(path, args.file)
        energy = get_free_energy(out_path)
        if energy is None:
            print_dual(color_text(f"[ERROR] Could not read energy from '{out_path}'.", 'red'), f_out)
            sys.exit(1)
        print_dual(f"  {label:<14} = {energy:>12.6f} eV  ({out_path})", f_out)
        report_quality_diagnostics(label, out_path, args.force_tolerance, f_out)
        energy_rows.append((label, energy, out_path))
        return energy

    report_path = os.path.join(output_root, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None
    print_dual(f"{color_text('===== OER STAGE 4 REPORT (ANALYSIS) =====', 'magenta')}", f_out)

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
    e_oh = read_energy_at(winning_oh_dir, "E_OH", f_out)
    e_o = read_energy_at(o_dir, "E_O", f_out)
    e_ooh = read_energy_at(ooh_dir, "E_OOH", f_out)
    e_h2 = read_energy("02_h2_molecule", "E_H2", f_out)
    e_h2o = read_energy("03_h2o_molecule", "E_H2O", f_out)
    e_deformed = read_energy("04_slab_deformed", "E_deformed", f_out)

    print_dual("", f_out)
    print_table(["Term", "Energy (eV)", "Source"],
                [([label, f"{energy:.6f}", path], None) for label, energy, path in energy_rows],
                f_out)

    print_section('[2] BSSE CORRECTION', f_out)
    e_deformation = e_deformed - e_clean
    print_dual(f"  Deformation cost (diagnostic, not used below) : {e_deformation:>10.4f} eV", f_out)

    def bsse_shift(prefix, f_out):
        # slab_only is read at THIS triad's own geometry (not the diagnostic
        # '04_slab_deformed', which is always at OH*'s geometry) -- the
        # counterpoise difference must compare the SAME geometry across basis
        # sets, and each intermediate gets its OWN triad at its own geometry
        # (see oer_refs.py's write_bsse_triad docstring for why a single
        # triad shared across all 3 was removed -- it systematically
        # under-corrected OH*/O*).
        e_slab_only = read_energy(f"{prefix}_slab_only", f"{prefix}_slab_only", f_out)
        e_slab_ghost = read_energy(f"{prefix}_slab_ghost", f"{prefix}_slab_gh", f_out)
        e_ads_ghost = read_energy(f"{prefix}_adsorbate_ghost_slab", f"{prefix}_ads_gh", f_out)
        e_isolated = read_energy(f"{prefix}_isolated", f"{prefix}_iso", f_out)
        bsse_slab = e_slab_only - e_slab_ghost
        bsse_ads = e_isolated - e_ads_ghost
        total = bsse_slab + bsse_ads
        print_dual(f"  BSSE ({prefix}): slab={bsse_slab:+.4f} eV  adsorbate={bsse_ads:+.4f} eV  "
                    f"total={total:+.4f} eV", f_out)
        return bsse_slab, bsse_ads, total

    bsse_oh_slab, bsse_oh_ads, bsse_oh = bsse_shift("05_bsse_OH", f_out)
    bsse_o_slab, bsse_o_ads, bsse_o = bsse_shift("05_bsse_O", f_out)
    bsse_ooh_slab, bsse_ooh_ads, bsse_ooh = bsse_shift("05_bsse_OOH", f_out)

    e_oh_corr = e_oh + bsse_oh
    e_o_corr = e_o + bsse_o
    e_ooh_corr = e_ooh + bsse_ooh
    print_dual("", f_out)
    print_table(["Intermediate", "E (raw, eV)", "BSSE (eV)", "E (corrected, eV)"], [
        (["OH*", f"{e_oh:.4f}", f"{bsse_oh:+.4f}", f"{e_oh_corr:.4f}"], None),
        (["O*", f"{e_o:.4f}", f"{bsse_o:+.4f}", f"{e_o_corr:.4f}"], None),
        (["OOH*", f"{e_ooh:.4f}", f"{bsse_ooh:+.4f}", f"{e_ooh_corr:.4f}"], None),
    ], f_out)

    # Flattened (label, value_eV, kind) rows -- 3 bars per intermediate
    # (slab component, adsorbate component, total), same shape
    # her_analysis.py's own energy-breakdown chart uses -- feeds
    # write_bsse_correction_plot/show_bsse_correction_matplotlib below, so
    # the printed numbers and the chart can never disagree.
    bsse_chart_rows = [
        ("OH* (slab)", bsse_oh_slab, "term"), ("OH* (ads.)", bsse_oh_ads, "term"),
        ("OH* (total)", bsse_oh, "total"),
        ("O* (slab)", bsse_o_slab, "term"), ("O* (ads.)", bsse_o_ads, "term"),
        ("O* (total)", bsse_o, "total"),
        ("OOH* (slab)", bsse_ooh_slab, "term"), ("OOH* (ads.)", bsse_ooh_ads, "term"),
        ("OOH* (total)", bsse_ooh, "total"),
    ]
    # (label, e_raw, bsse, e_corrected) -- feeds the raw-vs-corrected
    # small-multiples chart (write_bsse_raw_corrected_plot/
    # show_bsse_raw_corrected_matplotlib), the same numbers as the [2]
    # table above, just paired for that chart's own data format.
    bsse_summary_rows = [
        ("OH*", e_oh, bsse_oh, e_oh_corr),
        ("O*", e_o, bsse_o, e_o_corr),
        ("OOH*", e_ooh, bsse_ooh, e_ooh_corr),
    ]

    print_section('[3] THERMAL CORRECTION', f_out)

    def thermal_term(name, dir_suffix, phonon_label, f_out):
        if zpe_mode == "local":
            zpe_dir = os.path.join(output_root, f"08_zpe_calc_{dir_suffix}")
            zpe, ts = compute_local_zpe_entropy(zpe_dir, args.temp, f_out)
        else:
            zpe_dir = os.path.join(output_root, f"09_zpe_calc_{dir_suffix}")
            zpe, ts = compute_full_zpe_entropy(zpe_dir, phonon_label, args.temp, f_out)
        if zpe is None:
            print_dual(color_text(f"[ERROR] Could not compute {name}'s ZPE/entropy.", 'red'), f_out)
            if f_out:
                f_out.close()
            sys.exit(1)
        print_dual(f"  ZPE({name}) = {zpe:.4f} eV   TS({name}) = {ts:.4f} eV", f_out)
        return zpe, ts

    zpe_h2o, ts_h2o = thermal_term("H2O", "H2O", "oer_zpe_h2o", f_out)
    if zpe_mode == "full":
        zpe_clean, ts_clean = thermal_term("clean slab", "clean", "oer_zpe_clean", f_out)
    else:
        zpe_clean, ts_clean = 0.0, 0.0  # local mode: substrate frozen, delta vs. clean is 0 by construction

    zpe_oh_raw, ts_oh_raw = thermal_term("OH*", "OH", "oer_zpe_oh", f_out)
    zpe_o_raw, ts_o_raw = thermal_term("O*", "O", "oer_zpe_o", f_out)
    zpe_ooh_raw, ts_ooh_raw = thermal_term("OOH*", "OOH", "oer_zpe_ooh", f_out)

    # thermal_term()/compute_local_zpe_entropy()/compute_full_zpe_entropy()
    # all return ts_ev as an ALREADY temperature-integrated T*S term (eV,
    # not an eV/K entropy) -- see compute_local_zpe_entropy's own
    # 'ts_ev += temperature_k * S' line -- exactly like H2_TS_298K_EV
    # (already "T*S at 298.15 K", per its own module-level comment) and
    # HER's her_analysis.py, which never re-multiplies by temperature
    # either (its own delta_ts = ts_h_star - 0.5*H2_TS_298K_EV, no extra
    # '* args.temp'). d_ts_*/ts_h2o below are therefore ALREADY T*S, in
    # eV -- multiplying by args.temp again here was a real, previously
    # undetected bug (inflating every thermal correction by a factor of
    # args.temp, e.g. ~298x at room temperature) fixed together with the
    # zpe_local_meta.json fragment-label fix above.
    d_zpe_oh = (zpe_oh_raw - zpe_clean)
    d_ts_oh = (ts_oh_raw - ts_clean)
    d_zpe_o = (zpe_o_raw - zpe_clean)
    d_ts_o = (ts_o_raw - ts_clean)
    d_zpe_ooh = (zpe_ooh_raw - zpe_clean)
    d_ts_ooh = (ts_ooh_raw - ts_clean)

    # (label, zpe, ts, delta_zpe_or_None, delta_ts_or_None) -- feeds both the
    # print_table recap below AND write_markdown_report's [4] table, so the
    # two can never disagree.
    thermal_rows = [
        ("H2O", zpe_h2o, ts_h2o, None, None),
        ("OH*", zpe_oh_raw, ts_oh_raw, d_zpe_oh, d_ts_oh),
        ("O*", zpe_o_raw, ts_o_raw, d_zpe_o, d_ts_o),
        ("OOH*", zpe_ooh_raw, ts_ooh_raw, d_zpe_ooh, d_ts_ooh),
    ]
    print_dual("", f_out)
    print_table(["Species", "ZPE (eV)", "TS (eV)", "Delta-ZPE (eV)", "Delta-TS (eV)"], [
        ([row[0], f"{row[1]:.4f}", f"{row[2]:.4f}",
          "--" if row[3] is None else f"{row[3]:+.4f}",
          "--" if row[4] is None else f"{row[4]:+.4f}"], None)
        for row in thermal_rows
    ], f_out)

    g_oh = e_oh_corr + d_zpe_oh - d_ts_oh
    g_o = e_o_corr + d_zpe_o - d_ts_o
    g_ooh = e_ooh_corr + d_zpe_ooh - d_ts_ooh
    g_h2 = e_h2 + H2_ZPE_EV - H2_TS_298K_EV
    g_h2o = e_h2o + zpe_h2o - ts_h2o
    g_clean = e_clean

    print_section('[4] GAS-PHASE O2 REFERENCE (DERIVED)', f_out)
    g_o2 = 2.0 * g_h2o - 2.0 * g_h2 + FOUR_TIMES_U0_EV
    print_dual(f"  G(H2)  = {g_h2:.4f} eV   G(H2O,l) = {g_h2o:.4f} eV", f_out)
    print_dual(f"  G(O2) [derived, NOT computed via DFT] = 2*G(H2O) - 2*G(H2) + 4.92 = "
                f"{g_o2:.4f} eV", f_out)

    print_section('[5] REACTION STEPS & POTENTIAL-DETERMINING STEP', f_out)
    dG1 = g_oh + 0.5 * g_h2 - g_clean - g_h2o
    dG2 = g_o + 0.5 * g_h2 - g_oh
    dG3 = g_ooh + 0.5 * g_h2 - g_o - g_h2o
    dG4 = g_clean + g_o2 - g_ooh + 0.5 * g_h2
    steps = [("1: * + H2O -> OH*", dG1), ("2: OH* -> O*", dG2), ("3: O* + H2O -> OOH*", dG3),
             ("4: OOH* -> * + O2", dG4)]
    pds_index = int(np.argmax([dG1, dG2, dG3, dG4]))
    for i, (label, dg) in enumerate(steps):
        marker = color_text(" <-- PDS", 'yellow') if i == pds_index else ""
        print_dual(f"  Step {label:<22} dG = {dg:+.4f} eV{marker}", f_out)
    step_sum = dG1 + dG2 + dG3 + dG4
    sanity_ok = abs(step_sum - FOUR_TIMES_U0_EV) < 1e-6
    print_dual(f"  Sum dG1+dG2+dG3+dG4 = {step_sum:.6f} eV (should equal {FOUR_TIMES_U0_EV} eV "
                f"exactly by construction) -- {'OK' if sanity_ok else color_text('MISMATCH -- BUG', 'red')}",
                f_out)
    print_dual("", f_out)
    print_table(["Step", "Reaction", "dG (eV)"], [
        ([str(i + 1), label.split(": ", 1)[1], f"{dg:+.4f}"], 'yellow' if i == pds_index else None)
        for i, (label, dg) in enumerate(steps)
    ], f_out)

    eta = max(dG1, dG2, dG3, dG4) - 1.23
    print_dual(f"\n{color_text('[6] FINAL RESULT', 'magenta')}", f_out)
    print_dual("=" * 60, f_out)
    print_dual(f"  Theoretical overpotential eta = {eta:+.4f} V", f_out)
    print_dual(f"  Potential-determining step (PDS) = Step {pds_index + 1}", f_out)
    print_dual("=" * 60, f_out)
    verdict = ("excellent OER catalyst (eta close to 0)" if eta < 0.3
               else "moderate OER activity" if eta < 0.6
               else "poor OER catalyst (large overpotential)")
    print_dual(f"  Qualitative assessment: {verdict}", f_out)
    print_dual(color_text(
        "  [LIMITATION] AEM descriptor only (OH*/O*/OOH* at one site) -- does not model the "
        "lattice oxygen evolution mechanism (LOER). See --help.", 'cyan'), f_out)

    # Cumulative free energy along the 5-state reaction coordinate, at U=0
    # and at the equilibrium potential U=1.23V -- each electron-transfer
    # step shifts G by -eU, so the n-th state (n=0..4 electrons transferred
    # so far) shifts by -n*1.23 eV. Feeds both write_oer_free_energy_diagram_plot
    # and show_oer_free_energy_diagram_matplotlib below, so the printed
    # table and the chart can never disagree.
    g_u0 = [0.0, dG1, dG1 + dG2, dG1 + dG2 + dG3, dG1 + dG2 + dG3 + dG4]
    g_ueq = [g - n * 1.23 for n, g in enumerate(g_u0)]
    print_dual("", f_out)
    print_dual("  Cumulative free energy along the reaction coordinate (relative to * + H2O):", f_out)
    print_table(["State", "G at U=0 V (eV)", "G at U=1.23 V eq. (eV)"], [
        ([label, f"{g0:+.4f}", f"{geq:+.4f}"], None)
        for label, g0, geq in zip(OER_STATE_LABELS, g_u0, g_ueq)
    ], f_out)

    report_str_path = os.path.join(output_root, f"{args.output}.txt")
    with open(report_str_path, "w") as f_report:
        f_report.write(f"OER eta = {eta:+.4f} V, PDS = Step {pds_index + 1}\n")
        f_report.write(f"ZPE mode = {zpe_mode}, T = {args.temp} K\n")
        f_report.write(f"dG1={dG1:+.4f}  dG2={dG2:+.4f}  dG3={dG3:+.4f}  dG4={dG4:+.4f} eV\n")

    want_plot = args.plot
    if want_plot is None:
        try:
            answer = get_input(f"\nSave gnuplot charts (free-energy diagram + BSSE correction) to "
                                f"'{os.path.join(output_root, 'plot')}'? [y/N]: ").strip().lower()
        except EOFError:
            answer = ""
        want_plot = answer in ("y", "yes")
    if want_plot:
        plot_dir = os.path.join(output_root, "plot")
        dat_path, gplot_path = write_oer_free_energy_diagram_plot(plot_dir, g_u0, g_ueq, eta, pds_index)
        print_dual(f"\n{color_text('[Saved]', 'cyan')} {dat_path}, {gplot_path} "
                    f"(cd {plot_dir} && gnuplot {os.path.basename(gplot_path)} -- writes "
                    "oer_free_energy_diagram.pdf)", f_out)
        bsse_dat_path, bsse_gplot_path = write_bsse_correction_plot(plot_dir, bsse_chart_rows)
        print_dual(f"{color_text('[Saved]', 'cyan')} {bsse_dat_path}, {bsse_gplot_path} "
                    f"(cd {plot_dir} && gnuplot {os.path.basename(bsse_gplot_path)} -- writes "
                    "oer_bsse_correction.pdf)", f_out)
        bsse_rc_dat_path, bsse_rc_gplot_path = write_bsse_raw_corrected_plot(plot_dir, bsse_summary_rows)
        print_dual(f"{color_text('[Saved]', 'cyan')} {bsse_rc_dat_path}, {bsse_rc_gplot_path} "
                    f"(cd {plot_dir} && gnuplot {os.path.basename(bsse_rc_gplot_path)} -- writes "
                    "oer_bsse_raw_corrected.pdf)", f_out)
    else:
        dat_path = gplot_path = bsse_dat_path = bsse_gplot_path = None
        bsse_rc_dat_path = bsse_rc_gplot_path = None
        print_dual("\nNo plot generated.", f_out)

    want_show = args.show
    if want_show is None:
        try:
            answer = get_input("Show the charts on screen now (matplotlib, free-energy diagram "
                                "then BSSE correction)? [y/N]: ").strip().lower()
        except EOFError:
            answer = ""
        want_show = answer in ("y", "yes")
    if want_show:
        shown = show_oer_free_energy_diagram_matplotlib(g_u0, g_ueq, eta, pds_index)
        if shown:
            print_dual("Matplotlib window closed.", f_out)
        shown = show_bsse_correction_matplotlib(bsse_chart_rows)
        if shown:
            print_dual("Matplotlib window closed.", f_out)
        shown = show_bsse_raw_corrected_matplotlib(bsse_summary_rows)
        if shown:
            print_dual("Matplotlib window closed.", f_out)

    # Consolidated Markdown deliverable ("todos os dados de todas as
    # etapas") -- always written (like OER_report.txt above), not gated
    # behind a --markdown flag: it's cheap, headless-safe (the embedded PNGs
    # use fig.savefig, never plt.show()), and is exactly the kind of
    # write-once summary a user reaches for after --plot/--show finish.
    plot_dir = os.path.join(output_root, "plot")
    png_path = write_oer_free_energy_diagram_png(plot_dir, g_u0, g_ueq, eta, pds_index)
    bsse_png_path = write_bsse_correction_png(plot_dir, bsse_chart_rows)
    bsse_rc_png_path = write_bsse_raw_corrected_png(plot_dir, bsse_summary_rows)
    md_path = write_markdown_report(output_root, f"{args.output}.md", {
        'date': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'zpe_mode': zpe_mode, 'temp': args.temp,
        'winning_oh_label': winning_oh_label, 'o_dir': o_dir, 'ooh_dir': ooh_dir,
        'energy_rows': energy_rows,
        'e_oh': e_oh, 'e_o': e_o, 'e_ooh': e_ooh,
        'bsse_oh': bsse_oh, 'bsse_o': bsse_o, 'bsse_ooh': bsse_ooh,
        'e_oh_corr': e_oh_corr, 'e_o_corr': e_o_corr, 'e_ooh_corr': e_ooh_corr,
        'thermal_rows': thermal_rows,
        'g_h2': g_h2, 'g_h2o': g_h2o, 'g_o2': g_o2,
        'steps': steps, 'pds_index': pds_index, 'eta': eta, 'verdict': verdict,
        'g_u0': g_u0, 'g_ueq': g_ueq,
        'report_str_path': report_str_path, 'report_path': report_path,
        'dat_path': dat_path, 'gplot_path': gplot_path, 'png_path': png_path,
        'bsse_dat_path': bsse_dat_path, 'bsse_gplot_path': bsse_gplot_path,
        'bsse_png_path': bsse_png_path,
        'bsse_rc_dat_path': bsse_rc_dat_path, 'bsse_rc_gplot_path': bsse_rc_gplot_path,
        'bsse_rc_png_path': bsse_rc_png_path,
    })
    print_dual(f"\n{color_text('[Saved]', 'cyan')} {md_path} (consolidated Markdown report, all "
                "stages' data)", f_out)

    print_section('[7] SUMMARY & FILES', f_out)
    print_dual(f"Winning OH* site    : {winning_oh_label}", f_out)
    print_dual(f"O* folder           : {o_dir}", f_out)
    print_dual(f"OOH* folder         : {ooh_dir}", f_out)
    print_dual(f"eta                 : {eta:+.4f} V", f_out)
    print_dual(f"PDS                 : Step {pds_index + 1}", f_out)
    if report_path:
        print_dual(f"Report              : {report_path}", f_out)
    print_dual(f"Files               : {report_str_path}", f_out)
    print_dual(f"Markdown report     : {md_path}", f_out)
    if dat_path:
        print_dual(f"Plot data           : {dat_path}", f_out)
        print_dual(f"Plot script         : {gplot_path}", f_out)
        print_dual(f"BSSE plot data      : {bsse_dat_path}", f_out)
        print_dual(f"BSSE plot script    : {bsse_gplot_path}", f_out)
        print_dual(f"BSSE raw/corr. data : {bsse_rc_dat_path}", f_out)
        print_dual(f"BSSE raw/corr. plot : {bsse_rc_gplot_path}", f_out)

    if f_out:
        f_out.close()

    print("\n[INFO] Complete job!")
    print("\n" + "-" * 60)
    print(color_text("OER analysis complete.\n", 'bold'))


if __name__ == "__main__":
    main()
