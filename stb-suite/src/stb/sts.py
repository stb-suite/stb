#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""Simulated STS (scanning tunneling spectroscopy) dI/dV(E) curve at one
fixed real-space tip point, from a SIESTA full-BZ .WFSX + .fdf(+ion files).
Conceptually "PDOS projected onto a single real-space point instead of onto
atomic orbitals": for every (k, band, spin) eigenstate, |psi(r_tip)|^2 is
evaluated and Gaussian-broadened around its eigenvalue, then k-weighted and
summed into a continuous curve -- the Tersoff-Hamann-like proxy this
suite's stb-stm already uses for real-space STM images, extended to an
energy axis at one fixed point instead of a spatial map at one fixed
energy window.

Needs a real per-orbital radial basis (from a <label>.fdf with its
.ion/.ion.xml files alongside it -- NOT <label>.XV/.HSX, same requirement
as stb-wfdensity) and a full-BZ-sampled .WFSX -- NOT the band-path .WFSX
stb-fatbands uses, since this is a DOS-like, Brillouin-zone-integrated
quantity. SIESTA writes this kind of WFSX via 'WriteWaveFunctions T' +
an explicit '%block WaveFuncKPoints' list (named <label>.selected.WFSX by
SIESTA's own convention -- confirmed live: WriteWaveFunctions T alone,
without a WaveFuncKPoints block, silently writes nothing).

Known limitations:
 - Not a real STS simulation with a tip-orbital/decay model -- the tip is
   treated as a structureless point probe (same philosophy as stb-stm's
   constant-height/current modes).
 - No metadata distinguishes a full-BZ-mesh WFSX from a band-path one;
   this tool only prints an advisory (based on how uniform/small the
   k-point count and weights look), never blocks.
 - No --shift vbm/cbm (unlike stb-fatbands/stb-wfdensity) -- would need a
   second, expensive full pass over the whole WFSX; use --shift fermi
   with the Fermi-source hierarchy below instead.
 - Cost scales as O(nk * nbands) point-wavefunction evaluations --
   expect this to take a while for a dense mesh/many bands.
 - Only the first spinor component (spinor=0) is evaluated.
 - The --xy/--height tip-position mode assumes the chosen vacuum axis is
   Cartesian-aligned (same caveat stb-stm documents).
 - SIESTA's numerical atomic orbitals (PAOs) have a FINITE, hard cutoff
   radius (confined basis, e.g. ~2.5-3 Ang for a typical DZP basis) --
   |psi(r_tip)|^2 is exactly (not just numerically small) zero beyond
   that radius from every atom, unlike a plane-wave or Gaussian basis.
   Verified live: at a tip height beyond every orbital's cutoff, EVERY
   band gives exactly 0.0, producing a flat, uninformative curve with no
   error. Pick --height/--point within a few Angstrom of the surface
   (check the basis's PAO.EnergyShift-determined cutoff if unsure), or
   this tool will silently return nothing useful.

Output/report style rewritten (v1.0.0 -> v2.0.0) to match the rest of the
Analysis category (stb-wfdensity/stb-density): a numbered [0]...[6] report,
--save-report, --save-gnuplot (the .dat file used to be written
unconditionally with no way to opt out, and no .gplot script existed at
all despite the module being "gnuplot output" in spirit -- now both are
off by default and opted INTO together, with a real .gplot script
written alongside the .dat). --view replaces --no-plot -- the matplotlib
preview is now off by default and opted into, instead of on by default
and opted out of, matching stb-wfdensity's own convention.

Also adds --bands-file/--fermi-file and .out auto-detection for --shift
fermi's Fermi-energy requirement, alongside the pre-existing --fermi:
this tool used to only accept an explicit --fermi value, forcing the user
to look the number up by hand even though every finished SIESTA run
already has it in its own .out log. Reuses
core.siesta_bands.resolve_fermi_energy_hierarchy, the same priority-
ordered resolution (--fermi > --bands-file > --fermi-file > auto-detected
.out, decoupled from --label via core.siesta_log.find_out_file -- many
real SIESTA jobs redirect stdout to a generic name like calc.out instead
of <label>.out) first written for stb-wfdensity's --band vbm/cbm and
extracted to core/ once this tool became a second consumer.

v2.0.0 also fixes two more real bugs, found while comparing this tool
against stb-stm/stb-wfdensity's own already-fixed versions of the exact
same code:
 - --label + --geometry-file together used to be rejected outright, the
   same overly strict validation found and fixed in stb-wfdensity
   (SystemLabel "siesta" with the real input file named calc.fdf is a
   common real-world mismatch -- load_parent() already prefers an
   explicit --geometry-file over <label>.fdf on its own). --label +
   --wfsx is still rejected (--wfsx IS what --label auto-detects, so
   giving both is just ambiguous).
 - check_axis_alignment's returned warning (for the --xy/--height mode)
   was computed but never printed -- a real, silent gap; a sheared cell's
   "height above the surface" warning simply never reached the user.
   Fixed to print it, same as stb-stm already does.
 - the --xy/--height tip height used a naive `xyz[:, axis].max()` for the
   "topmost surface atom" reference -- the exact same bug found and fixed
   in stb-stm via core.kspace.find_surface_reference: a structure whose
   atoms straddle the periodic cell boundary (real vacuum gap in the
   MIDDLE of the cell, not padded after the atoms) silently picks the
   wrong bounding atom. Now reuses find_surface_reference directly,
   consistent with stb-stm's own fix.
"""

VERSION = "2.0.0"

import argparse
import os
import sys
import warnings
from datetime import datetime

import numpy as np

from stb.core import citations, kspace
from stb.core.cli import color_text, show_intro, print_dual, print_section, print_table
from stb.core.deps import require_sisl
from stb.core.siesta_bands import resolve_fermi_energy_hierarchy
from stb.core.siesta_wfsx import resolve_wfsx_path, load_parent, iter_wfsx_states_by_k, read_k_weights
from stb.core.broadening import sigma_ev_from_args
from stb.stm import resolve_axis, check_axis_alignment

REPORT_FILE = "stb_sts_report.txt"
BIB_FILE = "references.bib"


def main():
    parser = argparse.ArgumentParser(
        description="Simulated STS dI/dV(E) spectroscopy curve at a fixed tip point, from a "
                    "full-BZ SIESTA .WFSX.",
        epilog="Example usage:\n"
               "  stb-sts --label siesta --xy 0 0 --height 3.0 --erange -3 3 --sigma 50\n"
               "  stb-sts --label siesta --point 1.2 0.7 12.0 --erange -5 5 --fwhm 100\n"
               "  stb-sts --label siesta --xy 0 0 --height 3.0 --erange -5 5 --sigma 50 "
               "--shift fermi --fermi-file calc.out --save-report --save-gnuplot --view\n",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("--label", type=str, default=None,
                        help="SIESTA output label. Auto-detects <label>.WFSX, falling back to "
                             "<label>.selected.WFSX then <label>.bands.WFSX, and, unless "
                             "--geometry-file overrides it, <label>.fdf (with its .ion/.ion.xml "
                             "files alongside it). Mutually exclusive with --wfsx (not "
                             "--geometry-file -- see its own help: the real fdf is very often NOT "
                             "named <label>.fdf).")
    parser.add_argument("--wfsx", type=str, default=None,
                        help="Explicit path to the .WFSX file (alternative to --label).")
    parser.add_argument("--geometry-file", type=str, default=None,
                        help="Explicit .fdf geometry path (must have the basis's .ion/.ion.xml "
                             "files alongside it). Required if --wfsx is used instead of --label; "
                             "optional (and common) together WITH --label too, since the real fdf "
                             "is often not literally named <label>.fdf (e.g. SystemLabel 'siesta' "
                             "with the actual input file called calc.fdf) -- --label still "
                             "auto-detects the .WFSX in that case, only the geometry source changes.")

    point_group = parser.add_mutually_exclusive_group(required=True)
    point_group.add_argument("--point", type=float, nargs=3, default=None, metavar=("X", "Y", "Z"),
                        help="Absolute Cartesian tip position (Ang).")
    point_group.add_argument("--xy", type=float, nargs=2, default=None, metavar=("X", "Y"),
                        help="In-plane Cartesian position (Ang); needs --height too. "
                             "stb-stm-style: the tip sits --height Ang above the topmost "
                             "surface atom along the vacuum axis.")
    parser.add_argument("--height", type=float, default=None,
                        help="Height (Ang) above the topmost surface atom, for --xy mode.")
    parser.add_argument("--axis", type=int, default=None, choices=[0, 1, 2],
                        help="Surface-normal axis for --xy mode. Auto-detected from the "
                             "geometry's vacuum padding if omitted and exactly one vacuum axis "
                             "is found.")

    parser.add_argument("--erange", type=float, nargs=2, required=True, metavar=("EMIN", "EMAX"),
                        help="Energy window (eV), relative to --shift's reference.")
    parser.add_argument("--npoints", type=int, default=400,
                        help="Number of energy samples (default: 400).")
    parser.add_argument("--shift", type=str, choices=["none", "fermi", "manual"], default="none",
                        help="Energy reference (default: none, i.e. raw eigenvalues as stored "
                             "in the .WFSX). No vbm/cbm option here -- see the module docstring.")
    parser.add_argument("--fermi", type=float, default=None,
                        help="Fermi energy (eV) for --shift fermi. Highest-priority source if given.")
    parser.add_argument("--bands-file", type=str, default=None,
                        help="Companion .bands file to read the Fermi energy from, for --shift "
                             "fermi.")
    parser.add_argument("--fermi-file", type=str, default=None,
                        help="Explicit SIESTA .out log to read the Fermi energy from, for --shift "
                             "fermi -- an alternative to --bands-file for a run with no saved "
                             ".bands file. Not assumed to be named after --label.")
    parser.add_argument("--manual-value", type=float, default=None,
                        help="Custom shift value (eV), required if --shift manual.")

    broadening = parser.add_mutually_exclusive_group(required=True)
    broadening.add_argument("--sigma", type=float, default=None,
                        help="Gaussian broadening standard deviation, in meV.")
    broadening.add_argument("--fwhm", type=float, default=None,
                        help="Gaussian broadening full width at half maximum, in meV.")

    parser.add_argument("-o", "--output-dir", type=str, default=".",
                        help="Directory to write sts.dat (and, with --save-gnuplot/--save-report, "
                             "the .gplot script/report/references.bib) into (default: current "
                             "directory). Created if it doesn't exist.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the full run report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--save-gnuplot", action="store_true",
                        help="Also write a real .gplot script alongside sts.dat. Off by default -- "
                             "this tool used to write the .dat file unconditionally but never wrote "
                             "a .gplot script at all.")
    parser.add_argument("--view", action="store_true",
                        help="Show an interactive matplotlib preview before finishing. Off by "
                             "default (replaces the old --no-plot, which was on by default).")

    parser.add_argument("-v", "--version", action="version", version=f"stb-sts {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.xy is not None and args.height is None:
        parser.error("--height is required when --xy is given.")
    if args.shift == "manual" and args.manual_value is None:
        parser.error("--manual-value is required when --shift is 'manual'.")

    try:
        sigma_eV = sigma_ev_from_args(args.sigma, args.fwhm)
    except ValueError as e:
        parser.error(str(e))

    if args.label:
        if args.wfsx:
            parser.error("--label cannot be combined with --wfsx.")
        args.wfsx = resolve_wfsx_path(args.label, suffixes=(".WFSX", ".selected.WFSX", ".bands.WFSX"))
    elif not args.wfsx:
        parser.error("one of --label or --wfsx is required.")
    elif not args.geometry_file:
        parser.error("--geometry-file is required when using --wfsx instead of --label.")
    # --label + --geometry-file together is valid and common -- see the module docstring
    # (the same fix already made for stb-wfdensity): load_parent() below already prefers
    # an explicit geometry_file over <label>.fdf on its own.

    if args.wfsx is None or not os.path.isfile(args.wfsx):
        tried = (f"'{args.label}.WFSX', '{args.label}.selected.WFSX' and "
                 f"'{args.label}.bands.WFSX'") if args.label else f"'{args.wfsx}'"
        parser.error(f"No .WFSX file found ({tried}).")

    fermi_energy = fermi_source = None
    if args.shift == "fermi":
        fermi_energy, fermi_source = resolve_fermi_energy_hierarchy(
            args.fermi, args.bands_file, args.fermi_file, args.label)
        if fermi_energy is None:
            parser.error(
                "--shift fermi needs a Fermi energy -- none found. Pass --fermi <value>, "
                "--bands-file <path>, --fermi-file <path>, or leave a SIESTA .out log "
                "(<label>.out, or the sole *.out) in the current directory."
            )

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("STS: reading input data", 'bold'))
    print("-" * 60)

    sisl = require_sisl()
    # Silence sisl's own per-call SislInfo notices from wavefunction() --
    # both of the ones it prints here are expected/benign for this tool's
    # use: we deliberately translate the geometry outside its original
    # cell (see the r_tip trick below), and k != Gamma is the normal case
    # for a real k-mesh. Printed once per (k, band) state would otherwise
    # be hundreds of lines of noise for a realistic mesh.
    warnings.filterwarnings("ignore", category=sisl.SislInfo)

    print("[INFO] Resolving geometry (needs real per-orbital basis, from .fdf+ion files) ...")
    geometry, geo_source, _ = load_parent(args.label, None, args.geometry_file, mode="orbitals_required")
    if geometry is None:
        parser.error("No usable geometry found -- stb-sts needs <label>.fdf (with its "
                     ".ion/.ion.xml files alongside it) or --geometry-file.")
    print(f"[INFO] Using '{geo_source}' ({geometry.no} orbitals).")

    axis_warning = None
    if args.point is not None:
        r_tip = np.array(args.point)
        axis = None
        z_top = None
        print(f"[INFO] Tip position: {r_tip} Ang (absolute).")
    else:
        try:
            axis, _vacuum_axes = resolve_axis(geometry, args.axis, geometry.cell)
        except ValueError as e:
            print(color_text(f"[ERROR] {e}", 'red'))
            sys.exit(1)
        axis_warning = check_axis_alignment(geometry.cell, axis)
        if axis_warning:
            print(color_text(f"[WARNING] {axis_warning}", 'yellow'))
        axis_len = float(np.linalg.norm(geometry.cell[axis]))
        # frac_start is the atom immediately below the largest real vacuum gap --
        # NOT necessarily xyz[:, axis].max() (see find_surface_reference's own
        # docstring: a structure whose atoms straddle the periodic cell boundary,
        # with the real vacuum in the middle of the cell, needs the gap-aware
        # reference or "topmost atom" silently picks the wrong one).
        frac_start, _gap_size_frac = kspace.find_surface_reference(geometry.fxyz[:, axis])
        z_top = frac_start * axis_len
        tip_along_axis = z_top + args.height
        r_tip = np.empty(3)
        other_axes = [i for i in range(3) if i != axis]
        r_tip[other_axes[0]] = args.xy[0]
        r_tip[other_axes[1]] = args.xy[1]
        r_tip[axis] = tip_along_axis
        print(f"[INFO] Using axis {axis} as the surface normal; topmost atom at "
              f"{z_top:.3f} Ang; tip at {r_tip} Ang ({args.height:.3f} Ang above surface).")

    print(f"[INFO] Reading k-point weights from '{args.wfsx}' ...")
    k_weights = read_k_weights(args.wfsx, geometry)
    nk = len(k_weights)
    mesh_warning = None
    if nk <= 4 and np.allclose(k_weights, k_weights[0]):
        mesh_warning = (
            f"Only {nk} k-point(s), all equally weighted -- this doesn't look like a "
            "converged full-BZ mesh. stb-sts needs a genuine k-mesh WFSX (WriteWaveFunctions T "
            "+ %block WaveFuncKPoints on a normal SCF run), not a band-path WFSX like "
            "stb-fatbands uses -- the resulting curve may not be physically meaningful.")
        print(color_text(f"[WARNING] {mesh_warning}", 'yellow'))
    k_weights = k_weights / k_weights.sum()

    print(f"[INFO] Evaluating |psi(r_tip)|^2 for every (k, band, spin) state "
          f"({nk} k-points) -- this may take a moment ...")
    eigs = []
    weights2 = []
    kweights_flat = []

    # Fixed for the whole run (depends only on geometry/r_tip, neither of
    # which changes per k/spin) -- computed once here rather than
    # recomputed on every (k, spin) iteration below.
    geo_t = geometry.translate(-r_tip)

    for k_index, block in iter_wfsx_states_by_k(args.wfsx, geometry):
        for spin, state in block.items():
            nbands = state.state.shape[0]
            k_vec = state.info.get("k", (0.0, 0.0, 0.0))
            for b in range(nbands):
                one = state.sub(b, inplace=False)
                grid = sisl.Grid((1, 1, 1), lattice=geometry.lattice, dtype=complex)
                sisl.physics.electron.wavefunction(one.state, grid, geometry=geo_t, k=k_vec, spinor=0)
                psi2 = float(np.abs(grid.grid[0, 0, 0]) ** 2)
                eigs.append(float(np.asarray(one.eig).reshape(-1)[0]))
                weights2.append(psi2)
                kweights_flat.append(k_weights[k_index])

    eigs = np.array(eigs)
    weights2 = np.array(weights2)
    kweights_flat = np.array(kweights_flat)
    print(f"[INFO] Accumulated {len(eigs)} (k, band, spin) contributions.")

    # E stays in the user's requested (possibly Fermi-relative) reference
    # for the output x-axis; E_query is E shifted back to the raw,
    # absolute-eigenvalue reference the accumulated `eigs` are stored in --
    # same "shift the query grid, not the data" convention as stb-coop.
    E = np.linspace(args.erange[0], args.erange[1], args.npoints)
    if args.shift == "fermi":
        E_query = E + fermi_energy
    elif args.shift == "manual":
        E_query = E + args.manual_value
    else:
        E_query = E
    dist = sisl.physics.get_distribution("gaussian", smearing=sigma_eV)

    sts_curve = np.zeros_like(E)
    for eig, psi2, kw in zip(eigs, weights2, kweights_flat):
        sts_curve += kw * psi2 * dist(E_query - eig)

    zero_signal_warning = None
    if not np.any(weights2 > 0):
        zero_signal_warning = (
            "|psi(r_tip)|^2 is exactly zero for every state -- the tip point is likely "
            "beyond every orbital's basis cutoff radius (SIESTA's confined PAOs have a hard "
            "cutoff, unlike a plane-wave/Gaussian basis). Try a smaller --height/closer "
            "--point.")
        print(color_text(f"[WARNING] {zero_signal_warning}", 'yellow'))

    os.makedirs(args.output_dir, exist_ok=True)
    out_path = os.path.join(args.output_dir, "sts.dat")
    with open(out_path, "w") as f:
        f.write(f"# Generated by STB (stb-sts). Tip position: {r_tip} Ang\n")
        f.write(f"# Broadening: sigma = {sigma_eV * 1000:.3f} meV, shift = {args.shift}\n")
        f.write(f"#{'Energy(eV)':<14}{'dIdV(arb.units)':<18}\n")
        for e, v in zip(E, sts_curve):
            f.write(f"{e:<15.6f}{v:<18.6e}\n")

    gplot_path = None
    if args.save_gnuplot:
        gplot_path = out_path.rsplit('.', 1)[0] + ".gplot"
        dat_base = os.path.basename(out_path)
        stem = os.path.splitext(dat_base)[0]
        lines = [
            'set terminal pdfcairo enhanced font "Arial,14" size 8,6\n',
            f'set output "{stem}.pdf"\n',
            'set xlabel "Energy (eV)"\n',
            'set ylabel "dI/dV (arb. units)"\n',
            'set grid\n',
            f'plot "{dat_base}" using 1:2 with lines lw 2 title "Simulated STS at {r_tip} Ang"\n',
        ]
        with open(gplot_path, 'w') as f:
            f.writelines(lines)

    # --- From here on: the numbered, save-able report --------------------
    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    print_dual(color_text("\n===== STB-STS REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time      : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Label          : {args.label}" if args.label else "Label          : (explicit --wfsx/--geometry-file)", f_out)
    print_dual(f"WFSX file      : {args.wfsx}", f_out)
    print_dual(f"Geometry source: {geo_source}", f_out)
    print_dual(f"Tip position   : {r_tip} Ang" + (f" ({args.height:.3f} Ang above axis {axis})" if args.point is None else " (absolute)"), f_out)
    print_dual(f"Energy shift   : {args.shift}" + (f" (Fermi source: {fermi_source})" if args.shift == "fermi" else ""), f_out)
    print_dual(f"Broadening     : sigma = {sigma_eV * 1000:.3f} meV", f_out)
    print_dual(f"Output dir     : {args.output_dir}", f_out)
    print_dual(f"Save gnuplot   : {'yes' if args.save_gnuplot else 'no'}", f_out)
    print_dual(f"View (matplotlib): {'yes' if args.view else 'no'}", f_out)

    print_section("[1] INPUT DATA", f_out)
    print_table(["Quantity", "Value"], [
        (["Orbitals (basis)", f"{geometry.no}"], None),
        (["k-points in .WFSX", f"{nk}"], None),
        (["Energy samples", f"{args.npoints}"], None),
    ], f_out)
    if mesh_warning:
        print_dual(color_text(f"[WARNING] {mesh_warning}", 'yellow'), f_out)

    print_section("[2] TIP POSITION", f_out)
    if args.point is not None:
        print_table(["Quantity", "Value"], [
            (["Mode", "--point (absolute Cartesian)"], None),
            (["Position", f"{r_tip} Ang"], None),
        ], f_out)
    else:
        print_table(["Quantity", "Value"], [
            (["Mode", "--xy + --height (surface-relative)"], None),
            (["Surface-normal axis", f"{axis}"], None),
            (["Topmost surface atom", f"{z_top:.3f} Ang"], None),
            (["Height above surface", f"{args.height:.3f} Ang"], None),
            (["Absolute position", f"{r_tip} Ang"], None),
        ], f_out)
        if axis_warning:
            print_dual(color_text(f"[WARNING] {axis_warning}", 'yellow'), f_out)

    print_section("[3] STS CURVE", f_out)
    print_dual(f"(k, band, spin) contributions: {len(eigs)}", f_out)
    print_dual(f"Peak dI/dV                   : {sts_curve.max():.6e} arb. units "
               f"at E = {E[int(np.argmax(sts_curve))]:.4f} eV", f_out)
    if zero_signal_warning:
        print_dual(color_text(f"[WARNING] {zero_signal_warning}", 'yellow'), f_out)
    else:
        print_dual(color_text("[OK] Non-zero signal recovered across the requested energy window.", 'green'), f_out)

    print_section("[4] OUTPUT DATA & PLOTS", f_out)
    print_dual(color_text(f"[OK] Data written to '{out_path}'.", 'green'), f_out)
    if args.save_gnuplot:
        print_dual(color_text(f"[OK] Gnuplot script written to '{gplot_path}'.", 'green'), f_out)
    else:
        print_dual("Gnuplot script not written (off by default -- pass --save-gnuplot to write it).", f_out)

    print_section("[5] REFERENCES", f_out)
    bib_entries = [citations.SIESTA, citations.SIESTA_RECENT, citations.TERSOFF_HAMANN]
    citations.write_bib_file(os.path.join(args.output_dir, BIB_FILE), bib_entries)
    print_dual(color_text(
        f"[OK] Citations for the methods used in this run written to "
        f"'{os.path.join(args.output_dir, BIB_FILE)}' ({len(bib_entries)} entries).", 'green'), f_out)

    print_section("[6] SUMMARY & FILES", f_out)
    print_dual("Status         : OK", f_out)
    print_dual(f"Data           : {out_path}", f_out)
    if gplot_path:
        print_dual(f"Gnuplot script : {gplot_path}", f_out)
    print_dual(f"References     : {os.path.join(args.output_dir, BIB_FILE)}", f_out)
    if report_path:
        print_dual(f"Report         : {report_path}", f_out)

    if f_out:
        f_out.close()

    # --view runs last, after the report is fully printed/closed, so a
    # blocking matplotlib window never delays or hides it.
    if args.view:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(7, 5))
        plt.plot(E, sts_curve, color="#cc5522", lw=2)
        plt.xlabel("Energy (eV)")
        plt.ylabel("dI/dV (arb. units)")
        plt.title(f"Simulated STS at {r_tip} Ang")
        plt.grid(True)
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    main()
