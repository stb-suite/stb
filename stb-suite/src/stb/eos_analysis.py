#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""Stage 2 of the real-SIESTA equation-of-state workflow (5.15): aggregates the
'vol_*' folders written by stb-eosInputs and fits an equation of state to the
(volume, energy) pairs (the fitting math is core/eos_fit.py).

Volume comes from each folder's own SIESTA output (the LAST 'outcell:' block), not
from the input .fdf: the authoritative cell SIESTA actually used.

What the report adds around the fit itself, all from the same runs:
- the scan record eos_scan.json (when present): planned vs actual volumes, formula,
  reference structure, the input's distance from the fitted equilibrium;
- the level of theory of every folder (settings, k-grid, pseudopotential hashes),
  checked for being identical -- E(V) is one curve only if it is;
- an independent pressure P(V) from SIESTA's own stress tensor, compared with
  -dE/dV of the fitted curve (SIESTA's atom-centred orbitals keep their shape and
  range when the cell is scaled, so the two must agree; an inconsistency points to the
  real-space mesh (egg-box effect), an unconverged SCF or a poor fit), and a second bulk modulus, B = -V dP/dV, from the stress alone;
- uncertainties of V0 and B0 by leave-one-out (jackknife), their sensitivity to
  the fitted window, the agreement between EOS forms, residuals in meV/atom;
- the suggested LatticeConstant for the equilibrium volume;
- a pass/warn table of checks, and every warning again at the end.

Deliberately distinct from stb-elasticAnalysis's bulk modulus (a linear stress-strain
fit around one cell): agreement between the two is a physical consistency check of
the DFT setup.

Plots: --show opens a matplotlib dashboard and --png saves it; --save-gnuplot writes
the gnuplot .dat + .gplot pair (the convention of the workflow analysis stages);
--save-data writes a CSV table and a JSON with every number.
"""

from stb import __version__ as VERSION

import os
import sys
import json
import glob
import hashlib
import argparse
from datetime import datetime

import numpy as np

from stb.core import structure_io, siesta_log, fdf_directives
from stb.core.cli import color_text, show_intro, print_dual, print_section, print_table
from stb.core.eos_fit import normalize_eos_string, fit_eos, invert_pressure

REPORT_FILE = "stb_eosAnalysis_report.txt"
MANIFEST_FILE = "eos_scan.json"
_ALL_EOS_FORMS = ["birchmurnaghan", "vinet", "murnaghan", "sj"]
_EV_A3_TO_GPA = 160.21766208

# thresholds of the checks table
R2_GOOD = 0.9999
RESID_RMS_MEV_ATOM = 0.5          # fit residual rms
PRESSURE_DEV_GPA = 0.5            # |P_SIESTA - P_fit|
PRESSURE_AT_V0_GPA = 0.2          # residual SIESTA pressure at the fitted V0
FORCE_TOL = 0.05                  # eV/Ang, fixed-cell residual force
BPRIME_RANGE = (4.5, 7.0)
FORM_SPREAD_PCT = 1.0             # B0 spread between EOS forms
JACKKNIFE_PCT = 1.0               # B0 jackknife sigma / B0
WINDOW_PCT = 2.0                  # B0 change when refitting the central points only

_THEORY_KEYS = ["xcfunctional", "xcauthors", "paobasistype", "paobasissize", "paoenergyshift",
                "meshcutoff", "spin", "electronictemperature", "dmtolerance", "dftd3"]
_THEORY_LABELS = {"xcfunctional": "XC functional", "xcauthors": "XC authors", "paobasistype": "Basis type",
                  "paobasissize": "Basis size", "paoenergyshift": "PAO energy shift", "meshcutoff": "Mesh cutoff",
                  "spin": "Spin", "electronictemperature": "Electronic temperature",
                  "dmtolerance": "DM tolerance", "dftd3": "DFT-D3"}


# ---- reading --------------------------------------------------------------------

def folder_theory(folder_path):
    """{label: value} of the settings that shape the energy surface, as SIESTA reads them in this
    folder (config_extra.fdf first, then calc.fdf, then the rest), plus the k-grid and the md5 of
    each pseudopotential file."""
    files = sorted(glob.glob(os.path.join(folder_path, "*.fdf")))
    extra = [f for f in files if os.path.basename(f) == "config_extra.fdf"]
    main = [f for f in files if os.path.basename(f) == "calc.fdf"]
    rest = [f for f in files if f not in extra and f not in main]
    texts = []
    for f in extra + main + rest:
        try:
            with open(f, errors="replace") as fh:
                texts.append(fh.read())
        except OSError:
            pass
    d = fdf_directives.parse_directives(*texts)
    theory = {k: " ".join(d[k].split()) for k in _THEORY_KEYS if k in d}
    kg = fdf_directives.parse_kgrid(*texts)
    theory["kgrid"] = "x".join(str(v) for v in kg) if kg else (f"cutoff {d['kgridcutoff']}" if "kgridcutoff" in d else "Gamma")
    for pp in sorted(glob.glob(os.path.join(folder_path, "*.psml")) + glob.glob(os.path.join(folder_path, "*.psf"))):
        with open(pp, "rb") as fh:
            theory[f"pseudo {os.path.basename(pp)}"] = hashlib.md5(fh.read()).hexdigest()[:8]
    return theory


def count_atoms(folder_path):
    for candidate in sorted(glob.glob(os.path.join(folder_path, "*.fdf"))):
        if os.path.basename(candidate) == "config_extra.fdf":
            continue
        try:
            n = len(structure_io.read_fdf(candidate).atoms)
        except Exception:
            continue
        if n:
            return n
    return None


def read_manifest(directory):
    path = os.path.join(directory, MANIFEST_FILE)
    if not os.path.isfile(path):
        return None
    try:
        with open(path) as f:
            m = json.load(f)
        return m if m.get("format") == "stb-eos-scan" else None
    except (OSError, json.JSONDecodeError):
        return None


# ---- numerics -------------------------------------------------------------------

def curve_tools(fit):
    """Cubic-spline interpolators over the fitted curve: E(V), P(V) = -dE/dV (GPa), B(V) = V d2E/dV2 (GPa)."""
    from scipy.interpolate import CubicSpline
    cs = CubicSpline(fit["curve_v"], fit["curve_e"])
    e_of_v = lambda v: cs(v)
    p_of_v = lambda v: -cs(v, 1) * _EV_A3_TO_GPA
    b_of_v = lambda v: np.asarray(v) * cs(v, 2) * _EV_A3_TO_GPA
    return e_of_v, p_of_v, b_of_v


def jackknife(volumes, energies, form):
    """Leave-one-out standard errors of V0 and B0 (None when the refits fail)."""
    n = len(volumes)
    v0s, b0s = [], []
    for i in range(n):
        mask = np.arange(n) != i
        try:
            f = fit_eos(volumes[mask], energies[mask], form)
        except Exception:
            continue
        v0s.append(f["v0"])
        b0s.append(f["b0_gpa"])
    if len(v0s) < n:
        return None
    factor = (n - 1) / n
    return (float(np.sqrt(factor * np.sum((np.array(v0s) - np.mean(v0s)) ** 2))),
            float(np.sqrt(factor * np.sum((np.array(b0s) - np.mean(b0s)) ** 2))))


def stress_bulk_modulus(volumes, pressures, v0):
    """(B from stress, P at V0): polynomial fit of the SIESTA pressure against V -- cubic with 6+ points
    (a quadratic over a +/-8 % window truncates the curvature of P(V) and overestimates B by ~2 %), quadratic
    below -- and B = -V dP/dV at V0."""
    mask = np.isfinite(pressures)
    if mask.sum() < 4:
        return None, None
    c = np.polyfit(volumes[mask], pressures[mask], 3 if mask.sum() >= 6 else 2)
    dp = float(np.polyval(np.polyder(c), v0))
    return float(-v0 * dp), float(np.polyval(c, v0))


# ---- outputs --------------------------------------------------------------------

def write_eos_plot(dat_path, gplot_path, ctx):
    """<dat_path>: block 0 = scanned points (V, E, E/atom, P_SIESTA), block 1 = fitted curve
    (V, E, E/atom, P = -dE/dV); plus the companion .gplot with E(V) and P(V) side by side."""
    n = ctx["n_atoms"] or 1
    with open(dat_path, 'w') as f:
        f.write("# Equation of state | block 0: SIESTA points, block 1: fitted curve (" + ctx["fit_label"] + ")\n")
        f.write("# 1:Volume(Ang^3) 2:Energy(eV) 3:E/atom(eV) 4:P_siesta(GPa)\n")
        for v, e, p in zip(ctx["volumes"], ctx["energies"], ctx["pressures"]):
            f.write(f"{v:.6f}  {e:.6f}  {e / n:.6f}  {p:.5f}\n")
        f.write("\n\n# 1:Volume(Ang^3) 2:Energy(eV) 3:E/atom(eV) 4:P_fit(GPa)\n")
        for v, e, p in zip(ctx["curve_v"], ctx["curve_e"], ctx["curve_p"]):
            f.write(f"{v:.6f}  {e:.6f}  {e / n:.6f}  {p:.5f}\n")
    dat_name = os.path.basename(dat_path)
    base = os.path.splitext(os.path.basename(gplot_path))[0]
    with open(gplot_path, 'w') as f:
        f.writelines([
            '# --- STB Plot Configuration ---\n',
            '# Generated by stb-eosAnalysis\n',
            'set terminal pdfcairo enhanced color font "Arial,12" size 11,4.5\n',
            f'set output "{base}.pdf"\n',
            f'set multiplot layout 1,2 title "Equation of state (SIESTA) -- {ctx["fit_label"]}: V0 = {ctx["v0"]:.4f} Ang^3, B0 = {ctx["b0"]:.2f} GPa"\n\n',
            'set grid\nset key top center\n',
            'set xlabel "Volume (Ang^3)"\nset ylabel "Energy (eV)"\nset title "E(V)"\n',
            (f'plot "{dat_name}" index 0 using 1:2 with points pt 7 ps 1.4 lc rgb "#2255cc" title "SIESTA", '
             f'"{dat_name}" index 1 using 1:2 with lines lw 2 lc rgb "#cc5522" title "fit", '
             f'"+" using (({ctx["v0"]:.6f})):(({ctx["e0"]:.6f})) with points pt 13 ps 2 lc rgb "black" title "V0"\n\n'),
            'set ylabel "Pressure (GPa)"\nset title "P(V)"\nset key top right\n',
            (f'plot "{dat_name}" index 0 using 1:4 with points pt 7 ps 1.4 lc rgb "#2255cc" title "SIESTA stress", '
             f'"{dat_name}" index 1 using 1:4 with lines lw 2 lc rgb "#cc5522" title "-dE/dV of the fit", 0 lc rgb "gray" notitle\n'),
            'unset multiplot\n',
        ])


def make_dashboard(ctx, png_path, show):
    """2x3 matplotlib figure: E(V) with the fit(s), residuals, P(V) from the stress and from the fit,
    B(P), the B0 of every EOS form with the jackknife bar, and a text summary. Saved to `png_path`
    and/or shown on screen. Returns (produced, shown): `shown` is False when no display was available."""
    try:
        import matplotlib
        if not show:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print(color_text("[WARNING] matplotlib is not installed -- skipping the dashboard.", 'yellow'))
        return False, False
    if show and (matplotlib.get_backend().lower() == "agg"
                 or (sys.platform.startswith("linux") and not (os.environ.get("DISPLAY")
                                                               or os.environ.get("WAYLAND_DISPLAY")))):
        print(color_text("[WARNING] --show needs a display (local desktop or `ssh -X`); only saving the figure.", 'yellow'))
        show = False
        plt.switch_backend("Agg")

    n = ctx["n_atoms"] or 1
    unit = "per atom" if ctx["n_atoms"] else "cell"
    v, e = ctx["volumes"] / n, ctx["energies"] / n
    emin = ctx["e0"] / n
    fig, ax = plt.subplots(2, 3, figsize=(15, 8.5))
    blue, orange, green, grey = "#0072B2", "#D55E00", "#009E73", "#666666"

    a = ax[0, 0]
    colors = {"birchmurnaghan": orange, "vinet": green, "murnaghan": "#CC79A7", "sj": "#E69F00"}
    for form, fi in ctx["fits_by_form"].items():
        a.plot(fi["curve_v"] / n, 1000 * (fi["curve_e"] / n - emin), lw=2 if form == ctx["form"] else 1,
               ls="-" if form == ctx["form"] else "--", color=colors.get(form, grey),
               label=("fit: " if form == ctx["form"] else "") + form)
    a.plot(v, 1000 * (e - emin), "o", color=blue, ms=6, label="SIESTA")
    a.axvline(ctx["v0"] / n, color=grey, lw=0.8, ls=":")
    a.plot([ctx["v0"] / n], [0], "k*", ms=11, label=f"V0 = {ctx['v0'] / n:.3f}")
    a.set_xlabel(f"Volume ({unit}, Å³)")
    a.set_ylabel(f"E − E0 (meV {unit})")
    a.set_title("Energy vs volume")
    a.legend(fontsize=8)
    a.grid(alpha=0.3)

    a = ax[0, 1]
    res = ctx["residuals_mev"]
    a.bar(v, res, width=0.8 * np.min(np.diff(np.sort(v))) if len(v) > 1 else 0.5, color=blue)
    a.axhline(0, color="k", lw=0.8)
    lim = max(float(np.max(np.abs(res))) * 1.3, 0.05)
    a.set_ylim(-lim, lim)
    a.set_xlabel(f"Volume ({unit}, Å³)")
    a.set_ylabel(f"E_SIESTA − E_fit (meV {unit})")
    a.set_title(f"Fit residuals (rms {ctx['resid_rms']:.3f} meV {unit})")
    a.grid(alpha=0.3)

    a = ax[0, 2]
    a.plot(ctx["curve_v"] / n, ctx["curve_p"], color=orange, lw=2, label="−dE/dV of the fit")
    a.plot(v, ctx["pressures"], "o", color=blue, ms=6, label="SIESTA stress")
    a.axhline(0, color="k", lw=0.8)
    a.axvline(ctx["v0"] / n, color=grey, lw=0.8, ls=":")
    a.set_xlabel(f"Volume ({unit}, Å³)")
    a.set_ylabel("Pressure (GPa)")
    a.set_title(f"Pressure: max |ΔP| = {ctx['p_dev']:.3f} GPa")
    a.legend(fontsize=8)
    a.grid(alpha=0.3)

    a = ax[1, 0]
    a.plot(ctx["curve_p"], ctx["curve_b"], color=orange, lw=2)
    a.plot([0], [ctx["b0"]], "k*", ms=11, label=f"B0 = {ctx['b0']:.1f} GPa")
    if ctx["b_stress"] is not None:
        a.plot([ctx["p_at_v0"] if ctx["p_at_v0"] is not None else 0], [ctx["b_stress"]], "s", color=blue,
               label=f"B from stress = {ctx['b_stress']:.1f}")
    a.set_xlabel("Pressure (GPa)")
    a.set_ylabel("B = −V dP/dV (GPa)")
    a.set_title("Bulk modulus along the curve")
    a.legend(fontsize=8)
    a.grid(alpha=0.3)

    a = ax[1, 1]
    labels, values, bar_colors = [], [], []
    for f_, fi in ctx["fits_by_form"].items():
        short = f_.replace("birchmurnaghan", "Birch–\nMurnaghan")
        labels.append(short if len(ctx["fits_by_form"]) > 1 else f"{short}\n(all points)")
        values.append(fi["b0_gpa"])
        bar_colors.append(colors.get(f_, grey))
    if ctx["window_b0"] is not None:
        labels.append("central 60 %\nof the points")
        values.append(ctx["window_b0"])
        bar_colors.append("#56B4E9")
    if ctx["b_stress"] is not None:
        labels.append("from the\nstress")
        values.append(ctx["b_stress"])
        bar_colors.append(blue)
    a.bar(range(len(values)), values, color=bar_colors)
    if ctx["jack"]:
        a.axhspan(ctx["b0"] - ctx["jack"][1], ctx["b0"] + ctx["jack"][1], color=grey, alpha=0.25,
                  label=f"leave-one-out ±{ctx['jack'][1]:.2f} GPa")
        a.legend(fontsize=8, loc="upper right")
    span = max(max(values) - min(values), 1.0)
    a.set_ylim(min(values) - 2 * span, max(values) + 2 * span)
    for i, v_ in enumerate(values):
        a.text(i, v_, f"{v_:.2f}", ha="center", va="bottom", fontsize=8)
    a.set_xticks(range(len(values)))
    a.set_xticklabels(labels, fontsize=8)
    a.set_ylabel("B0 (GPa)")
    a.set_title("B0 by EOS form and by method")
    a.grid(alpha=0.3, axis="y")

    a = ax[1, 2]
    a.axis("off")
    lines = [f"{ctx['title']}", "",
             f"Fit ({ctx['fit_label']}):",
             f"  V0 = {ctx['v0']:.4f} Å³" + (f"  ({ctx['v0'] / n:.4f} /atom)" if ctx["n_atoms"] else ""),
             f"  E0 = {ctx['e0']:.6f} eV" + (f"  ({ctx['e0'] / n:.5f} /atom)" if ctx["n_atoms"] else ""),
             f"  B0 = {ctx['b0']:.2f} GPa" + (f" ± {ctx['jack'][1]:.2f}" if ctx["jack"] else ""),
             f"  B0' = {ctx['bprime']:.3f}" if ctx["bprime"] is not None else "  B0' = –",
             f"  R² = {ctx['r2']:.6f}", ""]
    if ctx["b_stress"] is not None:
        lines.append(f"B from stress: {ctx['b_stress']:.2f} GPa")
    if ctx["scale"] is not None:
        lines.append(f"V0/V_input: linear scale {ctx['scale']:.5f}")
    if ctx["new_lc"] is not None:
        lines.append(f"LatticeConstant → {ctx['new_lc']:.5f} Å")
    lines += ["", f"{ctx['n_ok']} of {ctx['n_total']} runs, checks: {ctx['n_pass']} ok / {ctx['n_warn']} warn"]
    a.text(0.02, 0.98, "\n".join(lines), va="top", ha="left", family="monospace", fontsize=10)

    fig.suptitle(f"Equation of state — {ctx['title']}", fontsize=14)
    fig.tight_layout()
    if png_path:
        fig.savefig(png_path, dpi=150)
    if show:
        plt.show()
    plt.close(fig)
    return True, show


# ---- main -----------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Stage 2 of 2 (equation of state): aggregates the 'vol_*' folders written by "
                     "stb-eosInputs and fits an equation of state (Birch-Murnaghan by default) to get "
                     "V0/E0/B0 (bulk modulus)/B0'. Also reads each run's SIESTA pressure (stress tensor) "
                     "and checks it against -dE/dV of the fit, gives a second B0 from the stress alone, "
                     "the uncertainty of V0/B0 (leave-one-out), the residuals, the consistency of the "
                     "level of theory across folders, and the suggested LatticeConstant; every warning "
                     "is repeated in the last section. --show / --png draw a dashboard, --save-gnuplot "
                     "writes the .dat/.gplot pair, --save-data a CSV + JSON.",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Example usage:\n"
               "  stb-eosAnalysis --dir eos_runs --show\n"
               "  stb-eosAnalysis --dir eos_runs --eos all --png --save-gnuplot --save-report\n"
               "  stb-eosAnalysis --dir eos_runs --eos vinet --target-pressure 0 5 10\n",
    )
    parser.add_argument("--dir", default="eos_runs",
                        help="Directory containing the 'vol_*' run folders (default: eos_runs).")
    parser.add_argument("--file", default="calc.out",
                        help="SIESTA output filename inside each run folder (default: calc.out).")
    parser.add_argument("--eos",
                        choices=["birchmurnaghan", "birch_murnaghan", "vinet", "murnaghan", "sj", "all"],
                        default="birchmurnaghan",
                        help="Equation-of-state form to fit (default: birchmurnaghan). 'sj' = stabilized "
                             "jellium (no B0' reported). 'all' fits every form on the same data and "
                             "reports them side by side (the reference fit for everything else is then "
                             "birchmurnaghan); mutually exclusive with --target-pressure.")
    parser.add_argument("--target-pressure", type=float, nargs='+', default=None, metavar="GPA",
                        help="One or more target external pressures (GPa) -- inverts the fitted EOS "
                             "(P = -dE/dV) to report the predicted equilibrium volume and linear "
                             "lattice-parameter scale factor at each one. Not available with --eos all.")
    parser.add_argument("-o", "--output", default="eos_curve",
                        help="Base filename (no extension) of the output files (default: eos_curve).")
    parser.add_argument("--show", action="store_true",
                        help="Open the analysis dashboard (matplotlib, 6 panels) on screen. Needs a display.")
    parser.add_argument("--png", action="store_true",
                        help="Save the dashboard as <output>_summary.png.")
    parser.add_argument("--save-gnuplot", action="store_true",
                        help="Also write <output>.dat and <output>.gplot (E(V) and P(V), data + fit); "
                             "'gnuplot <output>.gplot' renders <output>.pdf. Off by default.")
    parser.add_argument("--save-data", action="store_true",
                        help="Also write <output>_table.csv (one row per volume) and <output>_results.json "
                             "(every fit, check and number of the report). Off by default.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-eosAnalysis {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()
    args.eos = normalize_eos_string(args.eos)

    if args.eos == "all" and args.target_pressure:
        parser.error("--target-pressure needs one specific --eos form, not 'all'.")

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("EQUATION-OF-STATE WORKFLOW -- STAGE 2: ANALYSIS", 'bold'))
    print("-" * 60)

    if not os.path.isdir(args.dir):
        sys.exit(color_text(f"[ERROR] Directory '{args.dir}' not found.", 'red'))
    folders = sorted(f for f in os.listdir(args.dir)
                     if os.path.isdir(os.path.join(args.dir, f)) and f.startswith("vol_"))
    if not folders:
        sys.exit(color_text(f"[ERROR] No 'vol_*' folders found in '{args.dir}'.", 'red'))

    manifest = read_manifest(args.dir)
    plan = {p["folder"]: p for p in (manifest or {}).get("scan", {}).get("points", [])}

    report_path = REPORT_FILE if args.save_report else None
    f_out = open(report_path, "w") if report_path else None
    warnings = []

    def warn(text):
        warnings.append(text)
        return color_text(f"[WARNING] {text}", 'yellow')

    print_dual(color_text("===== STB-EOSANALYSIS REPORT =====", 'magenta'), f_out)

    # ------------------------------------------------------------------ [0]
    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Tool            : stb-eosAnalysis {VERSION}", f_out)
    print_dual(f"Date            : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Directory       : {args.dir}  ({len(folders)} 'vol_*' folder(s))", f_out)
    print_dual(f"SIESTA output   : {args.file}", f_out)
    print_dual(f"EOS form        : {args.eos}" + (" (every form; birchmurnaghan is the reference)" if args.eos == "all" else ""), f_out)
    print_dual(f"Scan record     : {MANIFEST_FILE if manifest else '(not found -- planned-vs-actual checks skipped)'}", f_out)

    # ------------------------------------------------------------------ [1]
    print_section("[1] SCAN & LEVEL OF THEORY", f_out)
    if manifest:
        st, sc, cl = manifest["structure"], manifest["scan"], manifest["calc"]
        print_dual(f"System          : {st['formula']}  ({st['n_atoms']} atoms)"
                   + (f", {st['space_group']} ({st['crystal_system']})" if st.get("space_group") else ""), f_out)
        print_dual(f"Input structure : {st['file']}, V = {st['volume']:.4f} Ang^3 ({st['volume'] / st['n_atoms']:.4f} /atom)", f_out)
        print_dual(f"Planned scan    : {sc['n_volumes']} volumes, +/-{sc['strain_range_pct']:g} % in volume"
                   f"; ions {'relaxed' if cl.get('relax_positions') else 'fixed'}", f_out)
        print_dual(f"Written by      : {manifest.get('tool', '?')} on {manifest.get('created', '?')}", f_out)
    theories = {f: folder_theory(os.path.join(args.dir, f)) for f in folders}
    keys = []
    for th in theories.values():
        for k in th:
            if k not in keys:
                keys.append(k)
    ref_theory = max(theories.values(), key=lambda t: list(theories.values()).count(t))
    rows, differing = [], {}
    for k in keys:
        vals = {f: th.get(k, "(absent)") for f, th in theories.items()}
        top = ref_theory.get(k, "(absent)")
        odd = [f for f, v in vals.items() if v != top]
        if odd:
            differing[k] = [(f, vals[f]) for f in odd]
        if not k.startswith("pseudo"):
            label = _THEORY_LABELS.get(k, "SCF k-grid" if k == "kgrid" else k)
            rows.append(([label, top, "IDENTICAL in all folders" if not odd else f"DIFFERS in {len(odd)} folder(s)"],
                         'yellow' if odd else None))
    pp = [f"{k[7:]} {ref_theory[k]}" for k in keys if k.startswith("pseudo") and k in ref_theory]
    print_table(["Setting", "Value", "Consistency"], rows, f_out)
    if pp:
        print_dual(f"Pseudopotentials: {', '.join(pp)} (md5, first 8)", f_out)
    for k, lst in differing.items():
        print_dual(warn(f"level of theory differs between folders for {k}: "
                        + ", ".join(f"{f} '{v}'" for f, v in lst[:5]) + (" ..." if len(lst) > 5 else "")
                        + " -- E(V) is one curve only if every folder uses the same settings."), f_out)

    # ------------------------------------------------------------------ [2]
    print_section("[2] READING FOLDERS", f_out)
    data = []     # dicts
    skipped = []
    for folder in folders:
        folder_path = os.path.join(args.dir, folder)
        out_path = os.path.join(folder_path, args.file)
        if not os.path.isfile(out_path):
            skipped.append((folder, f"missing {args.file}"))
            print_dual(f"   {folder:<14}: {color_text('SKIP', 'yellow')}  (missing {args.file} -- SIESTA not run here?)", f_out)
            continue
        energy = siesta_log.get_free_energy(out_path)
        cell = siesta_log.get_outcell(out_path)
        if energy is None or cell is None:
            skipped.append((folder, "could not parse energy/cell"))
            print_dual(f"   {folder:<14}: {color_text('SKIP', 'yellow')}  (could not parse the free energy / cell: unfinished run?)", f_out)
            continue
        converged, iterations = siesta_log.get_scf_convergence(out_path)
        stress = siesta_log.get_stress_tensor(out_path)
        pressure = float(-np.trace(stress) / 3 * _EV_A3_TO_GPA) if stress is not None else np.nan
        item = {"folder": folder, "volume": float(abs(np.linalg.det(cell))), "energy": float(energy),
                "scf": bool(converged), "iterations": iterations, "force": siesta_log.get_max_force(out_path),
                "pressure": pressure, "n_atoms": count_atoms(folder_path),
                "version": siesta_log.get_siesta_version(out_path), "plan": plan.get(folder)}
        data.append(item)
        flag = color_text("OK", 'green') if converged else color_text("OK (SCF not confirmed)", 'yellow')
        print_dual(f"   {folder:<14}: {flag}", f_out)
    if not data:
        print_dual(color_text("\n[FAIL] No valid data found in any 'vol_*' folder.", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)
    data.sort(key=lambda d: d["volume"])
    print_dual(f"\nRuns found      : {len(data)} of {len(folders)} folder(s)"
               + (f", {len(skipped)} skipped" if skipped else ""), f_out)
    if skipped:
        print_dual(warn(f"{len(skipped)} folder(s) skipped ({', '.join(f for f, _ in skipped)}); the fit uses the "
                        f"remaining {len(data)}."), f_out)
    n_atoms_set = {d["n_atoms"] for d in data if d["n_atoms"]}
    if len(n_atoms_set) > 1:
        print_dual(warn(f"folders report different atom counts ({sorted(n_atoms_set)}) -- make sure '{args.dir}' "
                        "holds one stb-eosInputs sweep only."), f_out)
    n_atoms = next(iter(n_atoms_set)) if len(n_atoms_set) == 1 else None
    versions = {d["version"] for d in data if d["version"]}
    if versions:
        print_dual(f"SIESTA version  : {', '.join(sorted(versions))}", f_out)

    # planned vs actual
    v_ref = None
    if manifest:
        v_ref = manifest["structure"]["volume"]
        off = [(d["folder"], d["volume"], d["plan"]["volume"]) for d in data
               if d["plan"] and abs(d["volume"] / d["plan"]["volume"] - 1) > 1e-3]
        for folder, v_act, v_plan in off:
            print_dual(warn(f"{folder}: SIESTA ran a cell of {v_act:.4f} Ang^3 but the plan says {v_plan:.4f} -- the "
                            "folder was edited or SIESTA used another cell (a %include problem?)."), f_out)
    else:
        zero = next((d for d in data if d["folder"] == "vol_0.00"), None)
        v_ref = zero["volume"] if zero else None

    # ------------------------------------------------------------------ [3]
    print_section("[3] VOLUME-ENERGY-PRESSURE TABLE", f_out)
    volumes = np.array([d["volume"] for d in data])
    energies = np.array([d["energy"] for d in data])
    pressures = np.array([d["pressure"] for d in data])
    e_unit = n_atoms or 1
    i_min = int(np.argmin(energies))
    rows = []
    for d in data:
        dv = (100 * (d["volume"] / v_ref - 1)) if v_ref else None
        scf_txt = "OK" if d["scf"] else "NOT CONFIRMED"
        rows.append(([d["folder"], f"{d['volume']:.4f}", f"{d['volume'] / e_unit:.4f}",
                      f"{dv:+.2f}" if dv is not None else "-", f"{d['energy']:.6f}",
                      f"{d['energy'] / e_unit:.6f}", f"{1000 * (d['energy'] - energies[i_min]) / e_unit:.2f}",
                      f"{d['pressure']:+.3f}" if np.isfinite(d['pressure']) else "-", scf_txt
                      + (f" ({d['iterations']})" if d["iterations"] else ""),
                      f"{d['force']:.5f}" if d["force"] is not None else "-"],
                     None if d["scf"] else 'yellow'))
    print_table(["Folder", "V (Ang^3)", "V/atom" if n_atoms else "V/cell", "dV (%)", "E (eV)",
                 "E/atom" if n_atoms else "E/cell", "E-Emin (meV/at)" if n_atoms else "E-Emin (meV)",
                 "P_SIESTA (GPa)", "SCF (iter)", "Fmax (eV/Ang)"], rows, f_out)
    print_dual("P_SIESTA = -trace(stress)/3 from each run's own stress tensor (positive = the cell wants to expand).", f_out)
    if len(volumes) < 4:
        print_dual(color_text(f"\n[FAIL] Only {len(volumes)} valid volume(s): a 4-parameter EOS fit needs at least 4 "
                              "(5+ recommended).", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)
    not_conv = [d["folder"] for d in data if not d["scf"]]
    if not_conv:
        print_dual(warn(f"SCF convergence not confirmed in {', '.join(not_conv)}: those energies may be off."), f_out)

    # ------------------------------------------------------------------ [4]
    form = "birchmurnaghan" if args.eos == "all" else args.eos
    fits_by_form = {}
    forms_to_fit = _ALL_EOS_FORMS if args.eos == "all" else [form]
    for fm in forms_to_fit:
        try:
            fits_by_form[fm] = fit_eos(volumes, energies, fm)
        except Exception as e:
            print_dual(warn(f"the {fm} fit failed ({e})."), f_out)
    if form not in fits_by_form:
        print_dual(color_text("[FAIL] The equation-of-state fit failed.", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)
    fit = fits_by_form[form]
    e_of_v, p_of_v, b_of_v = curve_tools(fit)
    v0, e0, b0, bprime, r2 = fit["v0"], fit["e0"], fit["b0_gpa"], fit["bprime"], fit["r_squared"]

    print_section("[4] EQUATION-OF-STATE FIT" + (" (ALL FORMS)" if args.eos == "all" else ""), f_out)
    resid = 1000 * (energies - e_of_v(volumes)) / e_unit
    resid_rms, resid_max = float(np.sqrt(np.mean(resid ** 2))), float(np.max(np.abs(resid)))
    fit_rows = []
    for fm, fi in fits_by_form.items():
        res_f = 1000 * (energies - curve_tools(fi)[0](volumes)) / e_unit
        fit_rows.append(([fm + (" *" if fm == form and args.eos == "all" else ""), f"{fi['v0']:.4f}",
                          f"{fi['v0'] / e_unit:.4f}", f"{fi['e0']:.6f}", f"{fi['b0_gpa']:.2f}",
                          f"{fi['bprime']:.3f}" if fi["bprime"] is not None else "-", f"{fi['r_squared']:.6f}",
                          f"{np.sqrt(np.mean(res_f ** 2)):.4f}",
                          "V0 OUTSIDE SCAN" if fi["v0_outside_range"] else "ok"],
                         'yellow' if fi["v0_outside_range"] else ('bold' if fm == form else None)))
    print_table(["EOS form", "V0 (Ang^3)", "V0/atom" if n_atoms else "V0/cell", "E0 (eV)", "B0 (GPa)", "B0'", "R^2",
                 "rms (meV/at)" if n_atoms else "rms (meV)", "Bracket"], fit_rows, f_out)
    if args.eos == "all":
        print_dual("(*) reference fit for the uncertainty, the pressure check, the plots and the files below.", f_out)
    else:
        print_dual("", f_out)
        print_dual(f"Equilibrium volume   V0 : {color_text(f'{v0:.4f}', 'bold')} Ang^3"
                   + (f"   ({v0 / n_atoms:.4f} Ang^3/atom)" if n_atoms else ""), f_out)
        print_dual(f"Equilibrium energy   E0 : {color_text(f'{e0:.6f}', 'bold')} eV"
                   + (f"   ({e0 / n_atoms:.6f} eV/atom)" if n_atoms else ""), f_out)
        print_dual(f"Bulk modulus         B0 : {color_text(f'{b0:.2f}', 'bold')} GPa", f_out)
        if bprime is not None:
            print_dual(f"Pressure derivative  B0': {color_text(f'{bprime:.3f}', 'bold')}", f_out)

    jack = jackknife(volumes, energies, form)
    if jack:
        print_dual(f"Uncertainty (leave-one-out): V0 +/- {jack[0]:.4f} Ang^3, B0 +/- {jack[1]:.2f} GPa "
                   f"({100 * jack[1] / b0:.2f} % of B0)", f_out)
    # sensitivity to the fitted window: central points only
    n_all = len(volumes)
    n_cent = max(5, int(round(0.6 * n_all)))
    window = None
    if n_all >= 7 and n_cent < n_all:
        lo = (n_all - n_cent) // 2
        try:
            fw = fit_eos(volumes[lo:lo + n_cent], energies[lo:lo + n_cent], form)
            window = (n_cent, fw["v0"], fw["b0_gpa"])
            print_dual(f"Window sensitivity         : refit with the central {n_cent} of {n_all} points -> V0 {fw['v0']:.4f}, "
                       f"B0 {fw['b0_gpa']:.2f} GPa ({fw['b0_gpa'] - b0:+.2f} GPa)", f_out)
        except Exception:
            window = None
    if args.eos == "all":
        b0_all = [fi["b0_gpa"] for fi in fits_by_form.values()]
        print_dual(f"Spread between EOS forms   : B0 {min(b0_all):.2f} .. {max(b0_all):.2f} GPa "
                   f"({100 * (max(b0_all) - min(b0_all)) / b0:.2f} % of B0), V0 {min(fi['v0'] for fi in fits_by_form.values()):.4f} .. "
                   f"{max(fi['v0'] for fi in fits_by_form.values()):.4f} Ang^3", f_out)
    print_dual(f"Fit quality                : R^2 {r2:.6f}; residuals rms {resid_rms:.4f}, max {resid_max:.4f} "
               f"meV/{'atom' if n_atoms else 'cell'}", f_out)
    kappa = 1000.0 / b0
    print_dual(f"Compressibility            : 1/B0 = {kappa:.3f} TPa^-1; B(V0) from the curve = {float(b_of_v(v0)):.2f} GPa", f_out)

    # geometry suggestion
    scale = (v0 / v_ref) ** (1.0 / 3.0) if v_ref else None
    new_lc = None
    if scale is not None:
        lc = (manifest or {}).get("structure", {}).get("latticeconstant") if manifest else None
        print_dual(f"Versus the input structure : V0 / V_input = {v0 / v_ref:.5f} ({100 * (v0 / v_ref - 1):+.2f} % in volume, "
                   f"linear scale {scale:.5f})", f_out)
        if lc:
            new_lc = lc * scale
            print_dual(f"Suggested LatticeConstant  : {lc:.5f} -> {new_lc:.5f} Ang (the equilibrium volume at the "
                       "input's shape)", f_out)

    # ------------------------------------------------------------------ [5]
    print_section("[5] PRESSURE CROSS-CHECK (SIESTA stress vs -dE/dV)", f_out)
    p_fit = p_of_v(volumes)
    have_p = np.isfinite(pressures)
    b_stress, p_at_v0 = stress_bulk_modulus(volumes, pressures, v0)
    p_dev = float(np.nanmax(np.abs(pressures - p_fit))) if have_p.any() else float("nan")
    if have_p.any():
        rows = [([d["folder"], f"{d['volume']:.4f}", f"{pressures[i]:+.4f}", f"{p_fit[i]:+.4f}",
                  f"{pressures[i] - p_fit[i]:+.4f}"], None) for i, d in enumerate(data)]
        print_table(["Folder", "V (Ang^3)", "P_SIESTA (GPa)", "P_fit = -dE/dV (GPa)", "Difference (GPa)"], rows, f_out)
        print_dual(f"Max |P_SIESTA - P_fit|     : {p_dev:.4f} GPa; mean difference {float(np.nanmean(pressures - p_fit)):+.4f} GPa", f_out)
        if p_at_v0 is not None:
            print_dual(f"SIESTA pressure at V0      : {p_at_v0:+.4f} GPa (residual pressure; ~0 if energy and stress agree)", f_out)
        if b_stress is not None:
            print_dual(f"B0 from the stress alone   : {b_stress:.2f} GPa (B = -V dP/dV at V0, polynomial fit of P_SIESTA) "
                       f"vs {b0:.2f} GPa from the energies ({b_stress - b0:+.2f} GPa, {100 * (b_stress - b0) / b0:+.2f} %)", f_out)
        print_dual(color_text("An independent second route to B0 on the same runs: with SIESTA's fixed-shape orbitals the stress "
                              "and -dE/dV must agree, so a disagreement points to the real-space mesh (raise Mesh.Cutoff), "
                              "an unconverged SCF (tighten DM.Tolerance) or a too narrow/noisy scan.", 'cyan'), f_out)
    else:
        print_dual("No stress tensor was found in these outputs -- the pressure cross-check is skipped.", f_out)

    # ------------------------------------------------------------------ [6]
    print_section("[6] TARGET PRESSURE", f_out)
    targets = []
    if args.target_pressure:
        rows, any_extrap = [], False
        for p in args.target_pressure:
            v_at_p, extrapolated = invert_pressure(fit, p)
            any_extrap = any_extrap or extrapolated
            scale_p = (v_at_p / v0) ** (1.0 / 3.0)
            targets.append({"pressure": p, "volume": float(v_at_p), "linear_scale_vs_v0": float(scale_p),
                            "extrapolated": bool(extrapolated)})
            rows.append(([f"{p:.3f}", f"{v_at_p:.4f}", f"{v_at_p / e_unit:.4f}", f"{scale_p:.6f}",
                          "EXTRAPOLATED" if extrapolated else "within the scan"], 'yellow' if extrapolated else None))
        print_table(["P (GPa)", "V (Ang^3)", "V/atom" if n_atoms else "V/cell", "L/L(V0)", "Status"], rows, f_out)
        if any_extrap:
            print_dual(warn("one or more target pressures lie outside the range the scan achieved: those predictions "
                            "rely on the EOS form's shape beyond the data; widen --strain-range in stb-eosInputs."), f_out)
    else:
        print_dual("Not requested (--target-pressure GPA [GPA ...] inverts the fit: the volume at each pressure).", f_out)

    # ------------------------------------------------------------------ [7]
    print_section("[7] CHECKS", f_out)
    n_below, n_above = int(np.sum(volumes < v0)), int(np.sum(volumes > v0))
    checks = []

    def check(name, ok, detail):
        checks.append((name, bool(ok), detail))
    n_planned = (manifest or {}).get("scan", {}).get("n_volumes") or len(folders)
    check("All planned runs available", len(data) == n_planned, f"{len(data)} of {n_planned}")
    check("SCF converged in every run", not not_conv, "yes" if not not_conv else f"not confirmed: {', '.join(not_conv)}")
    check("V0 inside the scanned range", not fit["v0_outside_range"],
          f"V0 {v0:.3f} in [{volumes.min():.3f}, {volumes.max():.3f}]")
    check("Points on both sides of V0 (>= 2 each)", n_below >= 2 and n_above >= 2, f"{n_below} below, {n_above} above")
    check(f"Fit R^2 >= {R2_GOOD}", r2 >= R2_GOOD, f"{r2:.6f}")
    check(f"Residual rms <= {RESID_RMS_MEV_ATOM} meV/{'atom' if n_atoms else 'cell'}", resid_rms <= RESID_RMS_MEV_ATOM, f"{resid_rms:.4f}")
    if bprime is not None:
        check(f"B0' within {BPRIME_RANGE[0]}-{BPRIME_RANGE[1]}", BPRIME_RANGE[0] <= bprime <= BPRIME_RANGE[1], f"{bprime:.3f}")
    if args.eos == "all":
        spread = 100 * (max(b0_all) - min(b0_all)) / b0
        check(f"EOS forms agree on B0 (< {FORM_SPREAD_PCT} %)", spread < FORM_SPREAD_PCT, f"{spread:.2f} %")
    if jack:
        check(f"B0 uncertainty (leave-one-out) < {JACKKNIFE_PCT} %", 100 * jack[1] / b0 < JACKKNIFE_PCT, f"{100 * jack[1] / b0:.2f} %")
    if window:
        check(f"B0 stable against the fitted window (< {WINDOW_PCT} %)", abs(window[2] - b0) / b0 * 100 < WINDOW_PCT,
              f"{100 * (window[2] - b0) / b0:+.2f} %")
    if have_p.any():
        check(f"P_SIESTA agrees with -dE/dV (< {PRESSURE_DEV_GPA} GPa)", p_dev < PRESSURE_DEV_GPA, f"max {p_dev:.3f} GPa")
        if p_at_v0 is not None:
            check(f"Residual pressure at V0 small (< {PRESSURE_AT_V0_GPA} GPa)", abs(p_at_v0) < PRESSURE_AT_V0_GPA, f"{p_at_v0:+.3f} GPa")
        if b_stress is not None:
            check("B0 from stress agrees with B0 from energy (< 5 %)", abs(b_stress - b0) / b0 * 100 < 5.0,
                  f"{100 * (b_stress - b0) / b0:+.2f} %")
    forces = [d["force"] for d in data if d["force"] is not None]
    if forces:
        check(f"Residual forces <= {FORCE_TOL} eV/Ang", max(forces) <= FORCE_TOL, f"max {max(forces):.5f}")
    check("Same level of theory in every folder", not differing, "identical" if not differing else f"differs: {', '.join(differing)}")
    rows = [([name, "OK" if ok else "WARN", detail], None if ok else 'yellow') for name, ok, detail in checks]
    print_table(["Check", "Status", "Detail"], rows, f_out)
    for name, ok, detail in checks:
        if ok:
            continue
        hint = {
            "V0 inside the scanned range": "the scan did not bracket the minimum -- widen --strain-range in stb-eosInputs.",
            "Points on both sides of V0 (>= 2 each)": "the minimum sits at the edge of the scan; re-centre the scan on V0 (use the suggested LatticeConstant).",
        }.get(name)
        if name.startswith("Fit R^2") or name.startswith("Residual rms"):
            hint = "the scan may be too wide (anharmonic) or too narrow (noise): try another --strain-range, and check the SCF settings."
        if name.startswith("P_SIESTA agrees") or name.startswith("Residual pressure") or name.startswith("B0 from stress"):
            hint = "energy and stress disagree: raise Mesh.Cutoff (egg-box effect of the real-space mesh), tighten DM.Tolerance, or widen the scan."
        if name.startswith("Residual forces"):
            hint = "the ions are not at equilibrium at the scaled cells: use stb-eosInputs --relax-positions."
        if name.startswith("Same level"):
            hint = "folders were edited after stb-eosInputs: rerun them with identical settings."
        if name == "SCF converged in every run":
            hint = "rerun the unconverged folders with more MaxSCFIterations / a smaller SCF.Mixer.Weight."
        print_dual(warn(f"{name}: {detail}." + (f" {hint}" if hint else "")), f_out)
    if v_ref and abs(v0 / v_ref - 1) > 0.03:
        print_dual(color_text(f"[NOTE] The input structure was {100 * abs(v0 / v_ref - 1):.1f} % away from the fitted "
                              "equilibrium volume: consider re-running stb-eosInputs from the suggested LatticeConstant.", 'cyan'), f_out)
    structure_info = (manifest or {}).get("structure", {})
    if structure_info.get("crystal_system") and structure_info["crystal_system"] != "cubic":
        print_dual(color_text(f"[NOTE] {structure_info['crystal_system']} crystal scaled isotropically: B0 is the bulk modulus "
                              "at fixed c/a and angles, not the fully relaxed solid's.", 'cyan'), f_out)

    # ------------------------------------------------------------------ [8]
    print_section("[8] SUMMARY & FILES", f_out)
    n_pass = sum(1 for _, ok, _ in checks if ok)
    n_warn = len(checks) - n_pass
    title = (structure_info.get("formula") or os.path.basename(os.path.abspath(args.dir)))
    print_dual(f"Result          : V0 = {v0:.4f} Ang^3, E0 = {e0:.6f} eV, B0 = {b0:.2f} GPa"
               + (f" +/- {jack[1]:.2f}" if jack else "") + (f", B0' = {bprime:.3f}" if bprime is not None else "")
               + f"  ({form}, R^2 {r2:.6f})", f_out)
    print_dual(f"Checks          : {n_pass} ok, {n_warn} warning(s)", f_out)

    curve_p = p_of_v(fit["curve_v"])
    curve_b = b_of_v(fit["curve_v"])
    ctx = {"volumes": volumes, "energies": energies, "pressures": pressures, "n_atoms": n_atoms, "e0": e0, "v0": v0, "b0": b0,
           "bprime": bprime, "r2": r2, "form": form, "fit_label": form, "fits_by_form": fits_by_form,
           "curve_v": fit["curve_v"], "curve_e": fit["curve_e"], "curve_p": curve_p, "curve_b": curve_b,
           "residuals_mev": resid, "resid_rms": resid_rms, "p_dev": p_dev if np.isfinite(p_dev) else 0.0,
           "b_stress": b_stress, "p_at_v0": p_at_v0, "jack": jack, "scale": scale, "new_lc": new_lc,
           "window_b0": (window[2] if window else None), "title": title, "n_ok": len(data), "n_total": len(folders), "n_pass": n_pass, "n_warn": n_warn}

    written = []
    if args.png or args.show:
        png_path = f"{args.output}_summary.png" if args.png else None
        produced, shown = make_dashboard(ctx, png_path, args.show)
        if produced:
            if png_path:
                written.append(png_path)
                print_dual(f"[OK] Dashboard      : {png_path}", f_out)
            if shown:
                print_dual("[OK] Dashboard shown on screen.", f_out)
    if args.save_gnuplot:
        dat_path, gplot_path = f"{args.output}.dat", f"{args.output}.gplot"
        write_eos_plot(dat_path, gplot_path, ctx)
        written += [dat_path, gplot_path]
        print_dual(f"[OK] gnuplot files  : {dat_path}, {gplot_path} (gnuplot {gplot_path} renders {args.output}.pdf)", f_out)
    if args.save_data:
        csv_path, json_path = f"{args.output}_table.csv", f"{args.output}_results.json"
        with open(csv_path, "w") as f:
            f.write("folder,volume_ang3,energy_ev,pressure_siesta_gpa,pressure_fit_gpa,residual_mev,scf_converged,max_force_ev_ang\n")
            for i, d in enumerate(data):
                f.write(f"{d['folder']},{d['volume']:.6f},{d['energy']:.6f},{d['pressure']:.5f},{p_fit[i]:.5f},"
                        f"{resid[i] * e_unit:.5f},{int(d['scf'])},{d['force'] if d['force'] is not None else ''}\n")
        results = {"tool": f"stb-eosAnalysis {VERSION}", "date": datetime.now().isoformat(timespec="seconds"),
                   "directory": args.dir, "form": form, "n_atoms": n_atoms,
                   "fits": {fm: {"v0": fi["v0"], "e0": fi["e0"], "b0_gpa": fi["b0_gpa"], "bprime": fi["bprime"],
                                 "r_squared": fi["r_squared"]} for fm, fi in fits_by_form.items()},
                   "jackknife": ({"v0_sigma": jack[0], "b0_sigma": jack[1]} if jack else None),
                   "window": ({"n_points": window[0], "v0": window[1], "b0_gpa": window[2]} if window else None),
                   "b0_from_stress_gpa": b_stress, "pressure_at_v0_gpa": p_at_v0, "max_pressure_difference_gpa": p_dev,
                   "residual_rms_mev": resid_rms * e_unit, "linear_scale_vs_input": scale, "suggested_latticeconstant": new_lc,
                   "targets": targets,
                   "checks": [{"name": n_, "ok": ok, "detail": dt} for n_, ok, dt in checks],
                   "points": [{"folder": d["folder"], "volume": d["volume"], "energy": d["energy"],
                               "pressure_siesta": (d["pressure"] if np.isfinite(d["pressure"]) else None),
                               "scf": d["scf"], "max_force": d["force"]} for d in data],
                   "warnings": warnings}
        with open(json_path, "w") as f:
            json.dump(results, f, indent=1)
            f.write("\n")
        written += [csv_path, json_path]
        print_dual(f"[OK] Data files     : {csv_path}, {json_path}", f_out)
    if not written and not args.show:
        print_dual("No plot or data file requested: --show / --png (dashboard), --save-gnuplot, --save-data.", f_out)
    if report_path:
        print_dual(f"[OK] Report saved to: {report_path}", f_out)

    # ------------------------------------------------------------------ [9]
    print_section("[9] NEXT STEPS", f_out)
    steps = []
    if new_lc and abs(v0 / v_ref - 1) > 0.01:
        steps.append(f"Re-centre the scan: set LatticeConstant {new_lc:.5f} Ang in the structure and rerun stb-eosInputs "
                     "(a scan centred on V0 constrains the fit best).")
    steps.append("Converge what shapes the curve: PAO.EnergyShift / basis size (the biggest effect on V0 and B0), SCF k-grid, "
                 "Mesh.Cutoff; rerun stb-eosInputs with each and compare this report.")
    steps.append("Cross-check B0: stb-elasticAnalysis (independent stress-strain route) and stb-mlelastic/stb-mlphonons --qha with a model fine-tuned on these settings.")
    if structure_info.get("crystal_system") and structure_info["crystal_system"] != "cubic":
        steps.append("Non-cubic crystal: relax the shape at V0 (cell with fixed volume) before comparing B0 with experiment.")
    for k, text in enumerate(steps, 1):
        print_dual(f"  {k}. {text}", f_out)

    # ------------------------------------------------------------------ [10]
    print_section("[10] WARNINGS", f_out)
    if warnings:
        for k, w in enumerate(warnings, 1):
            print_dual(color_text(f"{k:>2}. {w}", 'yellow'), f_out)
    else:
        print_dual(color_text("None.", 'green'), f_out)

    if f_out:
        f_out.close()

    print("\n[INFO] Complete job!")
    print("\n" + "-" * 60)
    print(color_text("Equation of state fitted.\n", 'bold'))


if __name__ == "__main__":
    main()
