#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""Stage 3 (final) of the Hirshfeld-I (iterative) charge-partitioning
workflow. Reads back every ion/<species>/cation//anion/ folder
stb-hirshfeldIons wrote (now assumed already run through SIESTA),
spherically-averages each into a radial density profile, then runs the
actual Hirshfeld-I convergence loop (core.hirshfeld.iterate_hirshfeld_i)
against the combined/production system's fixed real density: following
the original iterative-Hirshfeld formulation (Bultinck et al., J. Chem.
Phys. 126, 144111 (2007)) literally, EACH ATOM independently interpolates
its own reference between its species' neutral profile and whichever of
that species' cation/anion profile matches THAT ATOM's own charge sign
from the previous round -- two atoms of the same species can genuinely
lean toward opposite ion states in the same round. Iterates until
max|Delta q| converges.

This is the only stage that iterates -- no new folders/geometry are
written here, matching the suite's general "an Analysis stage never
generates new geometry" rule.

Report: the combined system (dimensionality, grid, electron count against
sum Z_val), every reference atom checked (electrons in the box and in the
radial profile against Z_val / Z_val -+ 1, radius holding 99.9 % of the
density against the box, SCF convergence), the live iteration, the
convergence history with the contraction ratio and the ESTIMATED distance
from the fixed point (max|Delta q| < tol only says the last step was
small), per-atom charges by every method (SIESTA Voronoi and Hirshfeld,
pass-0 Hirshfeld, Hirshfeld-I, the reference each atom ended on, |q| > 1
flagged), per-species statistics, charge balance, the point-charge dipole
against SIESTA's along the non-periodic axes, and ALL warnings in a last
section. Optional outputs (user's choice): --save-data (per-atom table,
per-iteration charges, extended XYZ with the charges), --png/--show (a
6-panel dashboard).
"""

from stb import __version__ as VERSION

import os
import sys
import time
import json
import argparse
import statistics
from datetime import datetime
from collections import defaultdict

import numpy as np

from stb.core import structure_io, kspace, hirshfeld as hf, siesta_log
from stb.core.deps import require_sisl
from stb.core.cli import color_text, show_intro, print_dual, print_section, print_table
from stb.core.rho_io import find_one_rho, read_total

sisl = require_sisl()

HIRSHFELD_IONS_MANIFEST_FILE = "hirshfeld_ions_manifest.json"
REPORT_FILE = "hirshfeld_analysis_report.txt"
DATA_FILE = "hirshfeld_charges.dat"
HISTORY_FILE = "hirshfeld_iterations.dat"
XYZ_FILE = "hirshfeld_charges.xyz"
PNG_FILE = "hirshfeld_analysis.png"
N_RADIAL_BINS = None  # None: bin width from the grid spacing (core.hirshfeld)
DEFAULT_TOL = 0.005
DEFAULT_MAX_ITER = 20
ELECTRON_COUNT_TOL = 0.02   # e-: reference / combined electron-count checks
REMAINING_ERROR_FLOOR = 1e-3  # e-: no remaining-error warning below this, whatever --tol
BOX_FILL_LIMIT = 0.8        # r(99.9 %) above this fraction of the box inradius -> warning
VACUUM_GAP = 5.0            # Ang: an axis with a gap this large is non-periodic
EANG_TO_DEBYE = 4.803204
AU_TO_DEBYE = 2.541746
STATES = ("neutral", "cation", "anion")
OKABE_ITO = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000"]
MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]


def load_species_profile(folder, n_bins=N_RADIAL_BINS):
    """Radial profile of one isolated-atom folder (neutral or ion) -- reads
    its own structure.fdf for the exact box/center rather than trusting a
    manifest value, then spherically-averages its '*.RHO'."""
    profile, _ = load_reference(folder, n_bins)
    return profile


def load_reference(folder, n_bins=N_RADIAL_BINS):
    """(profile, info) of one isolated-atom folder. info: electrons on the
    grid (what SIESTA put in the box), grid spacing, box inradius (half the
    shortest box height) and SCF convergence of the folder's log."""
    rho_path = find_one_rho(folder)
    grid = read_total(sisl.get_sile(rho_path))
    structure = structure_io.read_fdf(os.path.join(folder, "structure.fdf"))
    lattice = np.asarray(structure.lattice, dtype=float)
    center_frac = structure.atoms[0][1]
    if structure.coord_format == "cartesian":
        center_frac = kspace.to_fractional(center_frac, lattice, True)
    profile = hf.radial_density_profile(grid.grid, lattice, center_frac, n_bins=n_bins)
    volume = abs(np.linalg.det(lattice))
    heights = [volume / np.linalg.norm(np.cross(lattice[(i + 1) % 3], lattice[(i + 2) % 3]))
               for i in range(3)]
    out_file = siesta_log.find_out_file(folder, None)
    scf = siesta_log.get_scf_convergence(out_file) if out_file else None
    info = {
        "grid_electrons": float(np.sum(grid.grid) * volume / grid.grid.size),
        "spacing": grid_spacing(lattice, grid.grid.shape),
        "inradius": 0.5 * min(heights),
        "scf": scf,
    }
    return profile, info


def grid_spacing(lattice, shape):
    """Finest spacing (Ang) along the three cell vectors."""
    return min(float(np.linalg.norm(lattice[i])) / n for i, n in enumerate(shape))


def unwrap_vacuum_axes(frac, vacuum_axes):
    """Shifts the fractional coordinate of every non-periodic axis so its
    largest empty gap sits at the cell boundary -- a molecule or slab cut
    by the boundary is made whole before a dipole is taken."""
    frac = np.array(frac, dtype=float) % 1.0
    for axis, is_vacuum in enumerate(vacuum_axes):
        if not is_vacuum or len(frac) < 2:
            continue
        values = np.sort(frac[:, axis])
        gaps = np.diff(np.concatenate((values, [values[0] + 1.0])))
        start = values[(int(np.argmax(gaps)) + 1) % len(values)]
        frac[:, axis] = (frac[:, axis] - start) % 1.0
    return frac


def fmt(value, spec="+.4f", missing="N/A"):
    return missing if value is None else format(value, spec)


def reference_label(charge):
    """Which reference an atom ended on: neutral blended with the ion of its
    own sign, frac = min(|q|, 1) (core.hirshfeld.interpolate_reference_density)."""
    frac = min(abs(charge), 1.0)
    ion = "cation" if charge >= 0 else "anion"
    label = f"{frac:.2f} {ion}"
    return label + " (capped)" if abs(charge) > 1.0 else label


def make_dashboard(symbols, history, charges_history, charges_0, charges_final,
                   voronoi, native_h, profiles, ref_info, cart, tol, ratio, remaining,
                   png_path, show):
    """2x3 matplotlib dashboard: convergence, charge per iteration, methods
    per species, per-atom Hirshfeld-I vs Hirshfeld, the reference radial
    profiles and a map of the structure coloured by charge. Saved to
    png_path and/or shown. Returns True when something was produced."""
    try:
        import matplotlib
        if not show:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print(color_text("[WARNING] matplotlib is not installed -- skipping the dashboard.", 'yellow'))
        return False
    if show and (matplotlib.get_backend().lower() == "agg"
                 or (sys.platform.startswith("linux") and not (os.environ.get("DISPLAY")
                                                               or os.environ.get("WAYLAND_DISPLAY")))):
        print(color_text("[WARNING] --show needs a display (local desktop or `ssh -X`); "
                         "only saving the figure.", 'yellow'))
        show = False
        plt.switch_backend("Agg")
        if png_path is None:
            png_path = PNG_FILE

    species = sorted(set(symbols), key=symbols.index)
    color = {s: OKABE_ITO[i % len(OKABE_ITO)] for i, s in enumerate(species)}
    marker = {s: MARKERS[i % len(MARKERS)] for i, s in enumerate(species)}
    n_atoms = len(symbols)
    q_hist = np.array(charges_history)            # (n_rounds, n_atoms)
    rounds = np.arange(len(q_hist))

    fig, axes = plt.subplots(2, 3, figsize=(16, 9.5))
    ax_c, ax_q, ax_m, ax_s, ax_p, ax_map = axes.ravel()

    # (a) convergence
    its = [it for it, d in history if d is not None]
    ds = [d for _, d in history if d is not None]
    if ds:
        ax_c.semilogy(its, ds, "o-", color=OKABE_ITO[0], lw=2, ms=6, label=r"max|$\Delta q$|")
    ax_c.axhline(tol, color="0.4", ls="--", lw=1.2, label=f"--tol {tol:g} e-")
    if remaining is not None:
        ax_c.axhline(remaining, color=OKABE_ITO[1], ls=":", lw=1.5,
                     label=f"est. remaining error {remaining:.2g} e-")
    title = "Convergence"
    if ratio is not None:
        title += f"  (contraction ratio {ratio:.2f})"
    ax_c.set_title(title); ax_c.set_xlabel("iteration"); ax_c.set_ylabel("e-")
    ax_c.grid(alpha=0.3, which="both"); ax_c.legend(fontsize=8)

    # (b) charges along the iteration
    if n_atoms <= 24:
        seen = set()
        for i, s in enumerate(symbols):
            ax_q.plot(rounds, q_hist[:, i], "-", marker=marker[s], ms=4, lw=1.5, color=color[s],
                      label=s if s not in seen else None)
            seen.add(s)
        ax_q.set_title("Charge of every atom per iteration")
    else:
        for s in species:
            idx = [i for i, x in enumerate(symbols) if x == s]
            sub = q_hist[:, idx]
            ax_q.plot(rounds, sub.mean(axis=1), "-", marker=marker[s], ms=4, lw=2, color=color[s],
                      label=f"{s} (mean)")
            ax_q.fill_between(rounds, sub.min(axis=1), sub.max(axis=1), color=color[s], alpha=0.2)
        ax_q.set_title("Charge per species per iteration (band = min-max)")
    ax_q.axhline(0, color="0.5", lw=1)
    ax_q.set_xlabel("iteration (0 = simple Hirshfeld)"); ax_q.set_ylabel("charge (e-)")
    ax_q.grid(alpha=0.3); ax_q.legend(fontsize=8)

    # (c) methods per species (mean, error bar = std over the species' atoms)
    methods = []
    if any(v is not None for v in voronoi):
        methods.append(("Voronoi (SIESTA)", voronoi))
    if any(v is not None for v in native_h):
        methods.append(("Hirshfeld (SIESTA)", native_h))
    methods.append(("Hirshfeld (pass 0)", list(charges_0)))
    methods.append(("Hirshfeld-I", list(charges_final)))
    width = 0.8 / len(methods)
    # methods in grays (species own the colours elsewhere), Hirshfeld-I darkest
    method_colors = ["#D0D0D0", "#A8A8A8", "#7A7A7A", "#1A1A1A"][-len(methods):]
    x = np.arange(len(species))
    for k, (name, values) in enumerate(methods):
        means, stds = [], []
        for s in species:
            vals = [values[i] for i, x_ in enumerate(symbols) if x_ == s and values[i] is not None]
            means.append(np.mean(vals) if vals else np.nan)
            stds.append(np.std(vals) if len(vals) > 1 else 0.0)
        ax_m.bar(x + (k - (len(methods) - 1) / 2) * width, means, width * 0.92, yerr=stds,
                 color=method_colors[k], label=name, capsize=3, edgecolor="white", linewidth=1)
    ax_m.axhline(0, color="0.4", lw=1)
    ax_m.set_xticks(x); ax_m.set_xticklabels(species)
    ax_m.set_title("Mean charge per species, by method"); ax_m.set_ylabel("charge (e-)")
    ax_m.grid(alpha=0.3, axis="y"); ax_m.legend(fontsize=8)

    # (d) per-atom Hirshfeld-I against simple Hirshfeld
    for s in species:
        idx = [i for i, x_ in enumerate(symbols) if x_ == s]
        ax_s.scatter(np.asarray(charges_0)[idx], np.asarray(charges_final)[idx], s=50,
                     marker=marker[s], color=color[s], edgecolor="white", linewidth=1, label=s,
                     zorder=3)
    lo = min(np.min(charges_0), np.min(charges_final), 0.0)
    hi = max(np.max(charges_0), np.max(charges_final), 0.0)
    pad = 0.08 * (hi - lo or 1.0)
    ax_s.plot([lo - pad, hi + pad], [lo - pad, hi + pad], "--", color="0.5", lw=1, label="y = x")
    if np.ptp(charges_0) > 1e-6:
        slope = float(np.polyfit(charges_0, charges_final, 1)[0])
        ax_s.text(0.04, 0.95, f"slope {slope:.2f}", transform=ax_s.transAxes, va="top", fontsize=9)
    ax_s.set_xlabel("simple Hirshfeld (e-)"); ax_s.set_ylabel("Hirshfeld-I (e-)")
    ax_s.set_title("Per atom: Hirshfeld-I vs simple Hirshfeld")
    ax_s.grid(alpha=0.3); ax_s.legend(fontsize=8)

    # (e) reference radial profiles 4 pi r^2 rho(r)
    styles = {"neutral": "-", "cation": "--", "anion": ":"}
    r_lim = 0.0
    for s in species:
        for state in STATES:
            r, rho = profiles[state][s]
            if not np.any(rho):
                continue
            ax_p.plot(r, 4 * np.pi * r ** 2 * rho, styles[state], color=color[s], lw=2,
                      label=f"{s} {state}")
            r99 = ref_info[(s, state)].get("r999")
            if r99:
                r_lim = max(r_lim, r99)
    ax_p.set_xlim(0, 1.3 * r_lim if r_lim else None)
    ax_p.set_xlabel(r"r ($\AA$)"); ax_p.set_ylabel(r"4$\pi$r$^2\rho$(r)  (e-/$\AA$)")
    ax_p.set_title("Reference atoms (- neutral, -- cation, : anion)")
    ax_p.grid(alpha=0.3); ax_p.legend(fontsize=7, ncol=2)

    # (f) structure map coloured by the Hirshfeld-I charge, on the plane of
    # the two Cartesian axes with the largest spread of the atoms
    spread = np.ptp(cart, axis=0)
    a, b = sorted(np.argsort(spread)[-2:])
    vmax = max(float(np.max(np.abs(charges_final))), 1e-3)
    for s in species:
        idx = [i for i, x_ in enumerate(symbols) if x_ == s]
        sc = ax_map.scatter(cart[idx, a], cart[idx, b], c=np.asarray(charges_final)[idx],
                            cmap="RdBu_r", vmin=-vmax, vmax=vmax, s=110, marker=marker[s],
                            edgecolor="0.25", linewidth=0.8, label=s, zorder=3)
    if n_atoms <= 30:
        for i in range(n_atoms):
            ax_map.annotate(str(i + 1), (cart[i, a], cart[i, b]), textcoords="offset points",
                            xytext=(6, 6), fontsize=8, color="0.2")
    fig.colorbar(sc, ax=ax_map, label="Hirshfeld-I charge (e-)")
    axis_name = "XYZ"
    ax_map.set_xlabel(f"{axis_name[a]} ($\\AA$)"); ax_map.set_ylabel(f"{axis_name[b]} ($\\AA$)")
    ax_map.set_aspect("equal", adjustable="datalim")
    ax_map.set_title("Structure coloured by Hirshfeld-I charge")
    ax_map.grid(alpha=0.3); ax_map.legend(fontsize=8, loc="best")

    fig.suptitle("Hirshfeld-I charge analysis", fontsize=14)
    fig.tight_layout()
    if png_path:
        fig.savefig(png_path, dpi=150)
    if show:
        plt.show()
    plt.close(fig)
    return png_path or True


def main():
    parser = argparse.ArgumentParser(
        description="Stage 3 of the Hirshfeld-I workflow: reads back every neutral+ion "
                     "isolated-atom folder, then iterates the Hirshfeld-I reference-density "
                     "refinement against the combined/production system's real density until "
                     "convergence.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("-O", "--output-dir", default="hirshfeld_study",
                         help="Same root directory stb-hirshfeldPrep/Ions wrote into "
                              "(default: hirshfeld_study).")
    parser.add_argument("--tol", type=float, default=DEFAULT_TOL, metavar="E-",
                         help=f"Convergence tolerance on max|Delta q| between rounds, in "
                              f"electrons (default: {DEFAULT_TOL:g} e-, matching the method's "
                              f"own prototype). The report also estimates how far the charges "
                              f"still are from the fixed point, which can be 2-3 x this value "
                              f"when the iteration converges slowly.")
    parser.add_argument("--max-iter", type=int, default=DEFAULT_MAX_ITER, metavar="N",
                         help=f"Maximum number of refinement rounds (default: "
                              f"{DEFAULT_MAX_ITER}; the method's own prototype converged in "
                              "8).")
    parser.add_argument("--ref", required=False, default=None,
                         help="Path to a specific SIESTA output file to read each species' "
                              "Z_val (valence charge) from -- same convention as stb-bader's "
                              "own --ref. Overrides auto-detection of combined/'s own '*.out' "
                              "(which assumes that log was left with a .out extension; use "
                              "--ref if yours wasn't).")
    parser.add_argument("--save-report", action="store_true",
                         help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--save-data", action="store_true",
                         help=f"Write {DATA_FILE} (per-atom charges by every method, populations, "
                              f"reference used), {HISTORY_FILE} (every atom's charge per iteration) "
                              f"and {XYZ_FILE} (extended XYZ with the charges as per-atom columns; "
                              f"'ase gui' colours atoms by the Hirshfeld-I charge). Off by default.")
    parser.add_argument("--png", action="store_true",
                         help=f"Save the 6-panel dashboard (convergence, charge per iteration, methods "
                              f"per species, Hirshfeld-I vs Hirshfeld per atom, reference profiles, "
                              f"structure coloured by charge) as {PNG_FILE}. Off by default.")
    parser.add_argument("--show", action="store_true",
                         help="Open the dashboard on screen (matplotlib). Needs a display. Off by "
                              "default.")
    parser.add_argument("--no-intro", dest="intro", action="store_false",
                         help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version",
                         version=f"stb-hirshfeldAnalysis {VERSION}")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - Hirshfeld-I Analysis (Stage 3/3)",
            "Iterative Hirshfeld-I atomic charges",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    report_path = REPORT_FILE if args.save_report else None
    f_out = open(report_path, "w") if report_path else None
    warnings = []

    def warn(message):
        warnings.append(message)

    print_dual(color_text("===== STB-HIRSHFELDANALYSIS REPORT =====", 'magenta'), f_out)
    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time      : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Version        : {VERSION}", f_out)
    print_dual(f"Study dir      : {args.output_dir}", f_out)
    print_dual(f"Tolerance      : {args.tol:g} e-", f_out)
    print_dual(f"Max iterations : {args.max_iter}", f_out)
    print_dual(f"Radial bins    : {hf.BIN_WIDTH_PER_SPACING:g} x the grid spacing of each "
               "reference atom (empty bins dropped)", f_out)

    ions_manifest_path = os.path.join(args.output_dir, HIRSHFELD_IONS_MANIFEST_FILE)
    if not os.path.isfile(ions_manifest_path):
        print_dual(color_text(
            f"[ERROR] '{ions_manifest_path}' not found -- run stb-hirshfeldPrep then "
            "stb-hirshfeldIons first.", 'red'), f_out)
        sys.exit(1)
    with open(ions_manifest_path) as f:
        ions_manifest = json.load(f)
    # Paths are relative to the study folder; an older absolute path that
    # points at ANOTHER study (a copied/moved study) yields to this one's own.
    manifest_path, manifest_note = hf.resolve_study_path(
        args.output_dir, ions_manifest.get("hirshfeld_manifest"), hf.HIRSHFELD_MANIFEST_FILE)
    try:
        manifest, path_notes = hf.load_study_manifest(os.path.dirname(manifest_path) or ".")
    except (OSError, ValueError) as e:
        print_dual(color_text(f"[ERROR] Could not read '{manifest_path}': {e} -- run "
                              "stb-hirshfeldPrep first.", 'red'), f_out)
        sys.exit(1)
    for kind, text in ([manifest_note] if manifest_note else []) + path_notes:
        if kind == "info":
            print_dual(color_text(f"   [INFO] {text}", 'cyan'), f_out)
        else:
            warn(text)

    try:
        structure = structure_io.read_fdf(manifest["combined_structure"])
    except Exception as e:
        print_dual(color_text(
            f"[ERROR] Could not read the combined structure "
            f"'{manifest['combined_structure']}': {e}", 'red'), f_out)
        sys.exit(1)

    try:
        combined_rho_path = find_one_rho(manifest["combined_dir"])
    except ValueError as e:
        print_dual(color_text(
            f"[ERROR] {e} -- has SIESTA been run in 'combined/' yet?", 'red'), f_out)
        sys.exit(1)
    try:
        combined_sile = sisl.get_sile(combined_rho_path)
        rho_real = read_total(combined_sile).grid
    except Exception as e:
        print_dual(color_text(
            f"[ERROR] Could not read the combined '.RHO' '{combined_rho_path}': {e}", 'red'), f_out)
        sys.exit(1)

    symbols = [sym for sym, _ in structure.atoms]
    lattice = np.asarray(structure.lattice, dtype=float)
    positions = np.array([pos for _, pos in structure.atoms])
    frac_positions = kspace.to_fractional(positions, lattice,
                                           structure.coord_format == "cartesian")
    n_atoms = len(symbols)
    species = sorted(set(symbols), key=symbols.index)

    # Z_val (VALENCE charge) -- NOT the atomic number Z used to pick the
    # pseudopotential -- detected from combined/'s .out with the tabulated
    # fallback, same resolution as stb-bader's valence_source (using Z
    # instead of Z_val silently produces nonsense charges). --ref takes
    # priority over auto-detecting combined/'s own '*.out'.
    combined_out = args.ref or siesta_log.find_out_file(manifest["combined_dir"], None)
    if args.ref:
        print_dual(color_text(f"   [INFO] Using external reference file: {args.ref}", 'cyan'), f_out)
    detected_valence = siesta_log.get_zval_from_output(None, override_path=combined_out) \
        if combined_out else None
    valence_source = {**siesta_log.FALLBACK_VALENCE, **(detected_valence or {})}
    unknown_syms = sorted({sym for sym in symbols if sym not in valence_source})
    if unknown_syms:
        print_dual(color_text(
            f"[ERROR] No valence Z_val (detected or fallback) for species: "
            f"{', '.join(unknown_syms)}.", 'red'), f_out)
        sys.exit(1)
    z_vals = [valence_source[sym] for sym in symbols]
    fallback_syms = sorted(set(symbols) - set((detected_valence or {}).keys()))
    if fallback_syms:
        warn(f"Z_val for {', '.join(fallback_syms)} taken from the hardcoded fallback table (not "
             "detected in combined/'s own .out) -- verify these match your pseudopotential's "
             "valence charge.")

    # ---------------------------------------------------------------- [1]
    print_section("[1] COMBINED SYSTEM", f_out)
    vacuum_axes = kspace.detect_vacuum_axes(frac_positions, lattice, VACUUM_GAP)
    dim_label = kspace.dimensionality_label(vacuum_axes)
    cell_volume = abs(np.linalg.det(lattice))
    combined_spacing = grid_spacing(lattice, rho_real.shape)
    n_electrons = float(np.sum(rho_real) * cell_volume / rho_real.size)
    z_total = float(sum(z_vals))
    counts = {s: symbols.count(s) for s in species}
    print_dual(f"Structure      : {manifest['combined_structure']}", f_out)
    print_dual(f"Formula        : {''.join(f'{s}{counts[s]}' for s in species)} ({n_atoms} atoms)", f_out)
    print_dual(f"Cell volume    : {cell_volume:.4f} Ang^3", f_out)
    print_dual(f"Periodicity    : {dim_label} (non-periodic axes: "
               f"{', '.join('abc'[i] for i in range(3) if vacuum_axes[i]) or 'none'}; "
               f"vacuum gap >= {VACUUM_GAP:g} Ang)", f_out)
    print_dual(f"Density grid   : {' x '.join(map(str, rho_real.shape))} "
               f"(spacing {combined_spacing:.4f} Ang), from {combined_rho_path}", f_out)
    print_dual(f"Valence Z_val  : " + ", ".join(f"{s} {valence_source[s]:g}" for s in species)
               + f"  (sum {z_total:g} e-; "
               + ("detected from " + str(combined_out) if detected_valence else "fallback table")
               + ")", f_out)
    print_dual(f"Electrons      : {n_electrons:.4f} e- in the density "
               f"(sum Z_val - electrons = {z_total - n_electrons:+.4f} e-)", f_out)
    if abs(z_total - n_electrons) > ELECTRON_COUNT_TOL:
        warn(f"The combined density holds {n_electrons:.4f} e-, not sum Z_val = {z_total:g} e- "
             f"({z_total - n_electrons:+.4f}). Either the calculation is charged (NetCharge) -- "
             "then the charges add up to that net charge -- or Z_val does not match the "
             "pseudopotentials used.")
    if combined_out:
        conv, n_scf = siesta_log.get_scf_convergence(combined_out)
        print_dual(f"SCF            : {'converged after ' + str(n_scf) + ' iterations' if conv else 'NOT confirmed converged'}"
                   f" ({combined_out})", f_out)
        if not conv:
            warn(f"combined/ SCF convergence not confirmed in '{combined_out}' -- the density "
                 "the charges come from may be unconverged.")
    else:
        print_dual("SCF            : no SIESTA log found in combined/ (N/A)", f_out)

    # ---------------------------------------------------------------- [2]
    print_section("[2] REFERENCE ATOMS (neutral + cation + anion)", f_out)
    profiles = {state: {} for state in STATES}
    label_width = max(len(s) for s in species) + len(" (neutral) ")
    ref_info = {}
    sources = {}

    def read_ref(sym, state, folder):
        try:
            profile, info = load_reference(folder)
        except Exception as e:
            print_dual(color_text(
                f"[ERROR] Could not build the {state} profile of species '{sym}' from "
                f"'{folder}': {e} -- has SIESTA been run there yet?", 'red'), f_out)
            sys.exit(1)
        profiles[state][sym] = profile
        ref_info[(sym, state)] = info
        sources[(sym, state)] = folder
        print_dual(f"  {f'{sym} ({state})':<{label_width}}: read from '{folder}'", f_out)

    for entry in manifest["species"]:
        read_ref(entry["symbol"], "neutral", os.path.join(args.output_dir, entry["neutral_folder"]))
    for entry in ions_manifest["species"]:
        sym = entry["symbol"]
        if entry.get("cation_zero_density", False):
            # Hydrogen's (Z_val=1) cation is a bare proton: exactly zero
            # electron density everywhere -- no SIESTA folder exists, reuse
            # the neutral profile's radial grid with rho=0.
            r_neutral, _ = profiles["neutral"][sym]
            profiles["cation"][sym] = (r_neutral, np.zeros_like(r_neutral))
            ref_info[(sym, "cation")] = {"grid_electrons": 0.0, "spacing": None,
                                          "inradius": None, "scf": None, "analytic": True}
            sources[(sym, "cation")] = None
            print_dual(f"  {f'{sym} (cation)':<{label_width}}: analytical zero density (Z_val <= 1 -- bare "
                       "nucleus, no SIESTA folder)", f_out)
        else:
            read_ref(sym, "cation", os.path.join(args.output_dir, entry["cation_folder"]))
        read_ref(sym, "anion", os.path.join(args.output_dir, entry["anion_folder"]))

    ref_rows = []
    for sym in species:
        zv = valence_source[sym]
        for state, target in (("neutral", zv), ("cation", zv - 1), ("anion", zv + 1)):
            if sym not in profiles[state]:
                continue
            info = ref_info[(sym, state)]
            r, rho = profiles[state][sym]
            n_prof = hf.profile_electron_count(r, rho)
            info["profile_electrons"] = n_prof
            info["r999"] = hf.profile_radius(r, rho)
            flag = None
            if info.get("analytic"):
                ref_rows.append(([sym, state, f"{target:g}", "0 (exact)", "0 (exact)", "-", "-",
                                  "-", "-", "-"], None))
                continue
            n_grid = info["grid_electrons"]
            fill = info["r999"] / info["inradius"] if info["r999"] and info["inradius"] else None
            scf = info["scf"]
            scf_str = "N/A" if scf is None else (f"yes ({scf[1]})" if scf[0] else "NO")
            if abs(n_grid - target) > ELECTRON_COUNT_TOL:
                flag = 'yellow'
                warn(f"{sym} {state} reference holds {n_grid:.4f} e- on its grid, expected "
                     f"{target:g} (Z_val {zv:g}{'' if state == 'neutral' else (' - 1' if state == 'cation' else ' + 1')})"
                     f" -- wrong charge state, or Z_val does not match its pseudopotential "
                     f"('{sources[(sym, state)]}').")
            if abs(n_prof - n_grid) > ELECTRON_COUNT_TOL:
                flag = 'yellow'
                warn(f"{sym} {state}: the radial profile integrates to {n_prof:.4f} e- but the grid "
                     f"holds {n_grid:.4f} e- -- the spherical average is not resolving this atom "
                     "(too coarse a mesh?).")
            if fill is not None and fill > BOX_FILL_LIMIT:
                flag = 'yellow'
                warn(f"{sym} {state}: 99.9 % of its density lies within {info['r999']:.2f} Ang, "
                     f"{fill * 100:.0f} % of the box inradius ({info['inradius']:.2f} Ang) -- the "
                     "density reaches the box edge (an anion's extra electron may be held only by "
                     "the box/basis); rerun Stage 1 with a larger --vacuum.")
            if scf is not None and not scf[0]:
                flag = 'yellow'
                warn(f"{sym} {state}: SCF convergence not confirmed in '{sources[(sym, state)]}'.")
            ref_rows.append(([sym, state, f"{target:g}", f"{n_grid:.4f}", f"{n_prof:.4f}",
                              f"{n_prof - target:+.4f}", f"{info['r999']:.2f}" if info['r999'] else "-",
                              f"{info['inradius']:.2f}", f"{info['spacing']:.4f}", scf_str], flag))
    print_dual("", f_out)
    print_table(["Elem", "State", "Target N", "N on grid", "N radial", "Radial-target",
                 "r99.9%(Ang)", "Box inradius(Ang)", "Spacing(Ang)", "SCF conv."], ref_rows, f_out)
    print_dual("N on grid = electrons SIESTA put in the box; N radial = 4 pi int r^2 rho(r) of the "
               "spherical average actually used as reference; r99.9% = radius holding 99.9 % of it.",
               f_out)

    neutral_profiles, cation_profiles, anion_profiles = (profiles["neutral"], profiles["cation"],
                                                         profiles["anion"])
    species_profiles_0 = {i: neutral_profiles[symbols[i]] for i in range(n_atoms)}
    charges_0, populations_0 = hf.compute_hirshfeld_charges(
        rho_real, lattice, frac_positions, z_vals, species_profiles_0, progress_label="Pass 0")

    # ---------------------------------------------------------------- [3]
    print_section("[3] ITERATING HIRSHFELD-I", f_out)
    iter_start = time.monotonic()
    charges_history = [np.array(charges_0)]
    worst_by_iter = {}

    def report_iteration(iteration, max_delta, worst_atom_idx, charges, populations):
        elapsed = time.monotonic() - iter_start
        charges_history.append(np.array(charges))
        worst_by_iter[iteration] = worst_atom_idx
        sym = symbols[worst_atom_idx]
        leaning = "cation-like" if charges[worst_atom_idx] >= 0 else "anion-like"
        status = color_text("CONVERGED", 'green') if max_delta < args.tol else "not yet converged"
        print_dual(
            f"  Iteration {iteration:>2}/{args.max_iter}: max|Delta q| = {max_delta:.6f} e- "
            f"-- atom #{worst_atom_idx + 1} ({sym}), now {charges[worst_atom_idx]:+.4f} e- "
            f"({leaning}, pop. {populations[worst_atom_idx]:.4f} e-), {elapsed:6.1f}s elapsed "
            f"-- {status}",
            f_out)

    charges_final, populations_final, history = hf.iterate_hirshfeld_i(
        rho_real, lattice, frac_positions, symbols, z_vals,
        neutral_profiles, cation_profiles, anion_profiles, tol=args.tol, max_iter=args.max_iter,
        show_progress=True, initial_charges=charges_0, initial_populations=populations_0,
        on_iteration=report_iteration)
    iter_elapsed = time.monotonic() - iter_start

    # ---------------------------------------------------------------- [4]
    print_section("[4] CONVERGENCE HISTORY", f_out)
    conv_rows = []
    prev = None
    for iteration, max_delta in history:
        if max_delta is None:
            conv_rows.append((["0", "N/A (pass 0 = simple Hirshfeld)", "-", "-"], None))
            continue
        ratio_str = f"{max_delta / prev:.3f}" if prev else "-"
        w = worst_by_iter.get(iteration)
        conv_rows.append(([str(iteration), f"{max_delta:.6f}", ratio_str,
                           f"#{w + 1} ({symbols[w]})" if w is not None else "-"], None))
        prev = max_delta
    print_table(["Iteration", "max|Delta q| (e-)", "Ratio to previous", "Atom that moved most"],
                conv_rows, f_out)
    converged = history[-1][1] is not None and history[-1][1] < args.tol
    ratio, remaining = hf.estimate_remaining_error(history)
    if converged:
        print_dual(color_text(
            f"Converged after {history[-1][0]} iteration(s) (max|Delta q| < {args.tol:g} e-).",
            'green'), f_out)
    else:
        last = history[-1][1]
        warn(f"Did not confirm convergence within --max-iter {args.max_iter} iterations (last "
             f"max|Delta q| = {last if last is not None else float('nan'):.6f} e- >= {args.tol:g} e-). "
             "Re-run with a larger --max-iter, or treat the charges as unconverged.")
        print_dual(color_text("NOT converged within --max-iter (see WARNINGS).", 'yellow'), f_out)
    if ratio is not None and remaining is not None:
        print_dual(f"Contraction ratio      : {ratio:.3f} (each step ~{ratio * 100:.0f} % of the previous)",
                   f_out)
        print_dual(f"Est. remaining error   : {remaining:.2g} e- on the slowest atom, distance to "
                   f"the converged charges (last step x ratio / (1 - ratio))", f_out)
        if remaining > max(args.tol, REMAINING_ERROR_FLOOR):
            suggested = args.tol * (1 - ratio) / ratio
            warn(f"The last step was {history[-1][1]:.2g} e-, but at a contraction ratio of "
                 f"{ratio:.2f} the charges are still ~{remaining:.2g} e- from the converged "
                 f"Hirshfeld-I values ({remaining / args.tol:.1f} x --tol). For the error itself to "
                 f"be below {args.tol:g} e-, rerun with --tol {suggested:.2g} (Stage 3 only, no new SIESTA "
                 "run).")
    elif ratio is not None:
        warn(f"The last rounds did not contract (ratio {ratio:.2f}): the iteration may be "
             "oscillating -- check the per-iteration charges (--save-data / --png).")

    # Native SIESTA charges ('Charge.Voronoi end' / 'Charge.Hirshfeld end',
    # forced by stb-hirshfeldPrep into combined/'s config_extra.fdf), keyed
    # by SIESTA's own 1-based atom index.
    native_voronoi = siesta_log.get_voronoi_charges(combined_out) if combined_out else None
    native_hirsh = siesta_log.get_hirshfeld_charges(combined_out) if combined_out else None
    voronoi_by_id = {d["id"]: d["charge"] for d in native_voronoi} if native_voronoi else {}
    hirsh_by_id = {d["id"]: d["charge"] for d in native_hirsh} if native_hirsh else {}
    voronoi = [voronoi_by_id.get(i + 1) for i in range(n_atoms)]
    native_h = [hirsh_by_id.get(i + 1) for i in range(n_atoms)]

    # ---------------------------------------------------------------- [5]
    print_section("[5] PER-ATOM CHARGES: VORONOI vs. SIMPLE HIRSHFELD vs. HIRSHFELD-I", f_out)
    if not voronoi_by_id:
        print_dual(color_text(
            "[INFO] Native Voronoi charges not found in combined/'s own .out ('Voronoi "
            "Atomic Populations:' block, from 'Charge.Voronoi end') -- Voronoi column "
            "below reads N/A. Pass --ref to point at the right .out if 'combined/' has "
            "more than one, or if it was renamed.", 'yellow'), f_out)
    cart = frac_positions @ lattice
    headers = ["Idx", "Elem", "X(Ang)", "Y(Ang)", "Z(Ang)", "Voronoi(e-)", "Hirsh. SIESTA(e-)",
               "Simple Hirshfeld(e-)", "Hirshfeld-I(e-)", "HI - H(e-)", "Pop. HI(e-)",
               "Reference (neutral + frac ion)"]
    rows = []
    capped = []
    for i in range(n_atoms):
        q = float(charges_final[i])
        if abs(q) > 1.0:
            capped.append(i)
        rows.append(([str(i + 1), symbols[i], f"{cart[i, 0]:.4f}", f"{cart[i, 1]:.4f}",
                      f"{cart[i, 2]:.4f}", fmt(voronoi[i]), fmt(native_h[i]),
                      f"{charges_0[i]:+.4f}", f"{q:+.4f}", f"{q - charges_0[i]:+.4f}",
                      f"{populations_final[i]:.4f}", reference_label(q)],
                     'yellow' if abs(q) > 1.0 else None))
    print_table(headers, rows, f_out)
    print_dual("Charge = Z_val - population (positive = electrons lost). Reference: the atom's "
               "final reference density, (1 - frac) x neutral + frac x ion of its own sign, "
               "frac = min(|q|, 1).", f_out)
    if capped:
        warn(f"{len(capped)} atom(s) ended with |q| > 1 e- ("
             + ", ".join(f"#{i + 1} {symbols[i]} {charges_final[i]:+.3f}" for i in capped[:8])
             + (" ..." if len(capped) > 8 else "")
             + "): the reference is capped at the +-1 ion, so these charges are not fully "
             "self-consistent -- a +-2 ion reference would be needed for them.")

    # ---------------------------------------------------------------- [6]
    print_section("[6] PER-SPECIES SUMMARY (HIRSHFELD-I)", f_out)
    species_rows = []
    any_genuinely_mixed = False
    for sym in sorted(species):
        idx = [i for i, s in enumerate(symbols) if s == sym]
        values = [float(charges_final[i]) for i in idx]
        n_cation = sum(1 for c in values if c >= 0)
        n_anion = sum(1 for c in values if c < 0)
        is_mixed = n_cation > 0 and n_anion > 0
        any_genuinely_mixed = any_genuinely_mixed or is_mixed
        mean = statistics.mean(values)
        std = statistics.pstdev(values) if len(values) > 1 else 0.0
        mean_h = statistics.mean(float(charges_0[i]) for i in idx)
        voro_vals = [voronoi[i] for i in idx if voronoi[i] is not None]
        ratio_hi_h = f"{mean / mean_h:.2f}" if abs(mean_h) > 1e-3 else "-"
        species_rows.append(([sym, str(len(values)), str(n_cation), str(n_anion), f"{mean:+.4f}",
                              f"{std:.4f}", f"{min(values):+.4f}", f"{max(values):+.4f}",
                              fmt(statistics.mean(voro_vals) if voro_vals else None),
                              f"{mean_h:+.4f}", ratio_hi_h], 'cyan' if is_mixed else None))
    print_table(["Elem", "N", "Cation-like", "Anion-like", "Mean(e-)", "Std(e-)", "Min(e-)",
                 "Max(e-)", "Mean Voronoi(e-)", "Mean simple H.(e-)", "HI / H"],
                species_rows, f_out)
    if any_genuinely_mixed:
        print_dual(color_text(
            "Species highlighted above converged with atoms genuinely leaning toward BOTH "
            "ion states -- exactly the case a single-sign-per-species approach cannot "
            "represent correctly.", 'cyan'), f_out)
    print_dual("HI / H: how much Hirshfeld-I enlarges the simple-Hirshfeld charge (simple "
               "Hirshfeld's neutral references are known to understate charge transfer).", f_out)

    # ---------------------------------------------------------------- [7]
    print_section("[7] CHARGE BALANCE AND CROSS-CHECKS", f_out)
    pop_sum = float(np.sum(populations_final))
    print_dual(f"Sum Z_val                  : {z_total:.4f} e-", f_out)
    print_dual(f"Electrons in the density   : {n_electrons:.4f} e-", f_out)
    print_dual(f"Assigned to atoms (HI)     : {pop_sum:.4f} e-  (unassigned: "
               f"{n_electrons - pop_sum:+.4f} e-, density no reference atom reaches)", f_out)
    if abs(n_electrons - pop_sum) > ELECTRON_COUNT_TOL:
        warn(f"{n_electrons - pop_sum:.4f} e- of the combined density lie where no reference "
             "atom reaches and are assigned to no atom -- reference boxes too small?")
    totals = [("Voronoi (SIESTA)", voronoi), ("Hirshfeld (SIESTA)", native_h),
              ("Simple Hirshfeld (pass 0)", list(charges_0)), ("Hirshfeld-I", list(charges_final))]
    for name, values in totals:
        if all(v is not None for v in values):
            print_dual(f"Total charge, {name:<25}: {float(np.sum(values)):+.4f} e-", f_out)
    if all(v is not None for v in native_h):
        diff = np.abs(np.array(native_h) - np.array(charges_0))
        print_dual(f"Simple Hirshfeld here vs SIESTA's own: max |diff| {diff.max():.4f} e-, mean "
                   f"{diff.mean():.4f} e- (SIESTA uses its own free-atom densities from the basis; "
                   "this tool the SCF isolated atoms of neutral/ -- a small difference is expected).",
                   f_out)
        if diff.max() > 0.1:
            warn(f"Simple Hirshfeld (pass 0) differs from SIESTA's native Hirshfeld by up to "
                 f"{diff.max():.3f} e- -- larger than the reference-atom difference alone should "
                 "give; check that neutral/ used the same basis and pseudopotentials as combined/.")
    mean_abs = {name: float(np.mean(np.abs(values))) for name, values in totals
                if all(v is not None for v in values)}
    print_dual("Mean |q| per atom           : " + ", ".join(f"{k} {v:.4f}" for k, v in mean_abs.items()),
               f_out)

    # ---------------------------------------------------------------- [8]
    print_section("[8] DIPOLE CHECK (point charges vs SIESTA, non-periodic axes only)", f_out)
    dipole_rows = []
    if not any(vacuum_axes):
        print_dual("Skipped: the system is periodic along all three axes, where a point-charge "
                   "dipole is not defined.", f_out)
    else:
        frac_u = unwrap_vacuum_axes(frac_positions, vacuum_axes)
        cart_u = frac_u @ lattice
        siesta_dip = siesta_log.get_electric_dipole(combined_out) if combined_out else None
        siesta_debye = siesta_dip * AU_TO_DEBYE if siesta_dip is not None else None
        methods = [(name, values) for name, values in totals if all(v is not None for v in values)]
        for axis in range(3):
            if not vacuum_axes[axis]:
                continue
            u = lattice[axis] / np.linalg.norm(lattice[axis])
            row = [f"along {'abc'[axis]}"]
            for name, values in methods:
                row.append(f"{float(np.dot(np.asarray(values, float) @ cart_u, u)) * EANG_TO_DEBYE:+.3f}")
            row.append(fmt(float(np.dot(siesta_debye, u)) if siesta_debye is not None else None, "+.3f"))
            dipole_rows.append((row, None))
        if all(vacuum_axes):
            row = ["|mu|"]
            for name, values in methods:
                row.append(f"{np.linalg.norm(np.asarray(values, float) @ cart_u) * EANG_TO_DEBYE:.3f}")
            row.append(fmt(float(np.linalg.norm(siesta_debye)) if siesta_debye is not None else None, ".3f"))
            dipole_rows.append((row, None))
        print_table(["Component (Debye)"] + [name for name, _ in methods] + ["SIESTA (density)"],
                    dipole_rows, f_out)
        print_dual("Point charges alone leave out each atom's own dipole, so they need not match "
                   "SIESTA's dipole of the full density; the closer one says which charges carry "
                   "more of the electrostatics. SIESTA's value is read from 'Electric dipole' in "
                   "combined/'s .out." + ("" if siesta_debye is not None else
                                          " (not found there: column N/A)"), f_out)

    # ---------------------------------------------------------------- [9]
    print_section("[9] SUMMARY", f_out)
    big = int(np.argmax(np.abs(charges_final)))
    print_dual(f"Atoms analyzed  : {n_atoms}", f_out)
    print_dual(f"Iterations run  : {history[-1][0]} (of --max-iter {args.max_iter})", f_out)
    print_dual(f"Converged       : {'yes' if converged else 'NO'}"
               + (f" (est. remaining error {remaining:.2g} e-)" if remaining is not None else ""),
               f_out)
    print_dual(f"Iterating time  : {iter_elapsed:.1f}s", f_out)
    print_dual(f"Total charge    : {float(np.sum(charges_final)):+.4f} e- "
               f"(should be close to 0 for a neutral combined system)", f_out)
    print_dual(f"Largest charge  : atom #{big + 1} ({symbols[big]}) {charges_final[big]:+.4f} e- "
               f"(simple Hirshfeld {charges_0[big]:+.4f} e-)", f_out)
    print_dual("Per species     : " + ", ".join(
        f"{s} {statistics.mean(float(charges_final[i]) for i in range(n_atoms) if symbols[i] == s):+.3f}"
        for s in species) + " e- (mean Hirshfeld-I)", f_out)

    written = []
    if args.save_data:
        with open(DATA_FILE, "w") as f:
            f.write("# stb-hirshfeldAnalysis per-atom charges (e-; charge = Z_val - population)\n")
            f.write("# idx elem x(Ang) y(Ang) z(Ang) Z_val voronoi_siesta hirshfeld_siesta "
                    "hirshfeld_pass0 hirshfeld_i pop_hirshfeld_i reference_frac reference_ion\n")
            for i in range(n_atoms):
                q = float(charges_final[i])
                f.write(f"{i + 1:5d} {symbols[i]:>3} {cart[i, 0]:12.6f} {cart[i, 1]:12.6f} "
                        f"{cart[i, 2]:12.6f} {z_vals[i]:6.2f} {fmt(voronoi[i], '+10.5f', '       nan'):>10} "
                        f"{fmt(native_h[i], '+10.5f', '       nan'):>10} {charges_0[i]:+10.5f} "
                        f"{q:+10.5f} {populations_final[i]:10.5f} {min(abs(q), 1.0):6.3f} "
                        f"{'cation' if q >= 0 else 'anion'}\n")
        with open(HISTORY_FILE, "w") as f:
            f.write("# stb-hirshfeldAnalysis: every atom's charge (e-) per iteration "
                    "(0 = simple Hirshfeld); last column max|Delta q|\n")
            f.write("# iter " + " ".join(f"{s}{i + 1}" for i, s in enumerate(symbols)) + " max_dq\n")
            for (it, d), q in zip(history, charges_history):
                f.write(f"{it:4d} " + " ".join(f"{x:+.6f}" for x in q)
                        + f" {d if d is not None else float('nan'):.6f}\n")
        try:
            from ase import Atoms
            from ase.io import write as ase_write
            atoms = Atoms(symbols=symbols, positions=cart, cell=lattice,
                          pbc=[not v for v in vacuum_axes])
            atoms.set_initial_charges(np.asarray(charges_final, float))
            atoms.new_array("hirshfeld_i", np.asarray(charges_final, float))
            atoms.new_array("hirshfeld", np.asarray(charges_0, float))
            if all(v is not None for v in voronoi):
                atoms.new_array("voronoi", np.asarray(voronoi, float))
            ase_write(XYZ_FILE, atoms, format="extxyz")
            written += [DATA_FILE, HISTORY_FILE, XYZ_FILE]
        except ImportError:
            warn(f"ase is not installed -- {XYZ_FILE} not written.")
            written += [DATA_FILE, HISTORY_FILE]
    if args.png or args.show:
        produced = make_dashboard(symbols, history, charges_history, charges_0, charges_final,
                                  voronoi, native_h, profiles, ref_info, cart, args.tol, ratio,
                                  remaining, PNG_FILE if args.png else None, args.show)
        if isinstance(produced, str):
            written.append(produced)
    if report_path:
        written.append(report_path)
    print_dual(f"Files written   : {', '.join(written) if written else 'none (add --save-report, --save-data, --png)'}",
               f_out)

    # ---------------------------------------------------------------- [10]
    print_section("[10] WARNINGS", f_out)
    if warnings:
        for k, message in enumerate(warnings, 1):
            print_dual(color_text(f"  [WARNING {k}] {message}", 'yellow'), f_out)
    else:
        print_dual(color_text("  none", 'green'), f_out)

    if f_out:
        f_out.close()

    print("\n" + "-" * 60)
    print(color_text("Hirshfeld-I Analysis (Stage 3) complete.\n", 'bold'))


if __name__ == "__main__":
    main()
