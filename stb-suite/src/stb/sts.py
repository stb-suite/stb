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
   second, expensive full pass over the whole WFSX; pass --fermi directly.
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
"""

VERSION = "1.0.0"

import argparse
import os
import sys
import warnings
import numpy as np
import matplotlib.pyplot as plt

from stb.core.cli import color_text, show_intro
from stb.core.deps import require_sisl
from stb.core.siesta_wfsx import resolve_wfsx_path, load_parent, iter_wfsx_states_by_k, read_k_weights
from stb.core.broadening import sigma_ev_from_args
from stb.stm import resolve_axis, check_axis_alignment


def main():
    parser = argparse.ArgumentParser(
        description="Simulated STS dI/dV(E) spectroscopy curve at a fixed tip point, from a "
                    "full-BZ SIESTA .WFSX.",
        epilog="Example usage:\n"
               "  stb-sts --label siesta --xy 0 0 --height 3.0 --erange -3 3 --sigma 50\n"
               "  stb-sts --label siesta --point 1.2 0.7 12.0 --erange -5 5 --fwhm 100\n",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("--label", type=str, default=None,
                        help="SIESTA output label. Auto-detects <label>.WFSX, falling back to "
                             "<label>.selected.WFSX then <label>.bands.WFSX, and <label>.fdf "
                             "(with its .ion/.ion.xml files alongside it). Mutually exclusive "
                             "with --wfsx/--geometry-file.")
    parser.add_argument("--wfsx", type=str, default=None,
                        help="Explicit path to the .WFSX file (alternative to --label).")
    parser.add_argument("--geometry-file", type=str, default=None,
                        help="Explicit .fdf geometry path (must have the basis's .ion/.ion.xml "
                             "files alongside it). Required if --wfsx is used instead of --label.")

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
                        help="Fermi energy (eV), required if --shift fermi.")
    parser.add_argument("--manual-value", type=float, default=None,
                        help="Custom shift value (eV), required if --shift manual.")

    broadening = parser.add_mutually_exclusive_group(required=True)
    broadening.add_argument("--sigma", type=float, default=None,
                        help="Gaussian broadening standard deviation, in meV.")
    broadening.add_argument("--fwhm", type=float, default=None,
                        help="Gaussian broadening full width at half maximum, in meV.")

    parser.add_argument("-o", "--output-dir", type=str, default=".",
                        help="Directory to write sts.dat into (default: current directory). "
                             "Created if it doesn't exist.")
    parser.add_argument("--no-plot", dest="plot", action="store_false",
                        help="Skip the matplotlib preview (the data file is still written).")
    parser.add_argument("-v", "--version", action="version", version=f"stb-sts {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.xy is not None and args.height is None:
        parser.error("--height is required when --xy is given.")
    if args.shift == "fermi" and args.fermi is None:
        parser.error("--fermi is required when --shift is 'fermi'.")
    if args.shift == "manual" and args.manual_value is None:
        parser.error("--manual-value is required when --shift is 'manual'.")

    try:
        sigma_eV = sigma_ev_from_args(args.sigma, args.fwhm)
    except ValueError as e:
        parser.error(str(e))

    if args.label:
        if args.wfsx or args.geometry_file:
            parser.error("--label cannot be combined with --wfsx/--geometry-file.")
        args.wfsx = resolve_wfsx_path(args.label, suffixes=(".WFSX", ".selected.WFSX", ".bands.WFSX"))
    elif not args.wfsx:
        parser.error("one of --label or --wfsx is required.")
    elif not args.geometry_file:
        parser.error("--geometry-file is required when using --wfsx instead of --label.")

    if args.wfsx is None or not os.path.isfile(args.wfsx):
        tried = (f"'{args.label}.WFSX', '{args.label}.selected.WFSX' and "
                 f"'{args.label}.bands.WFSX'") if args.label else f"'{args.wfsx}'"
        parser.error(f"No .WFSX file found ({tried}).")

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("STS SIMULATOR:", 'bold'))
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

    if args.point is not None:
        r_tip = np.array(args.point)
        print(f"[INFO] Tip position: {r_tip} Ang (absolute).")
    else:
        try:
            axis, _vacuum_axes = resolve_axis(geometry, args.axis, geometry.cell)
        except ValueError as e:
            print(color_text(f"[ERROR] {e}", 'red'))
            sys.exit(1)
        check_axis_alignment(geometry.cell, axis)
        z_top = float(geometry.xyz[:, axis].max())
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
    if nk <= 4 and np.allclose(k_weights, k_weights[0]):
        print(color_text(
            f"[WARNING] Only {nk} k-point(s), all equally weighted -- this doesn't look like a "
            "converged full-BZ mesh. stb-sts needs a genuine k-mesh WFSX (WriteWaveFunctions T "
            "+ %block WaveFuncKPoints on a normal SCF run), not a band-path WFSX like "
            "stb-fatbands uses -- the resulting curve may not be physically meaningful.",
            'yellow'))
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
        E_query = E + args.fermi
    elif args.shift == "manual":
        E_query = E + args.manual_value
    else:
        E_query = E
    dist = sisl.physics.get_distribution("gaussian", smearing=sigma_eV)

    sts_curve = np.zeros_like(E)
    for eig, psi2, kw in zip(eigs, weights2, kweights_flat):
        sts_curve += kw * psi2 * dist(E_query - eig)

    if not np.any(weights2 > 0):
        print(color_text(
            "[WARNING] |psi(r_tip)|^2 is exactly zero for every state -- the tip point is "
            "likely beyond every orbital's basis cutoff radius (SIESTA's confined PAOs have a "
            "hard cutoff, unlike a plane-wave/Gaussian basis). Try a smaller --height/closer "
            "--point.", 'yellow'))

    os.makedirs(args.output_dir, exist_ok=True)
    out_path = os.path.join(args.output_dir, "sts.dat")
    with open(out_path, "w") as f:
        f.write(f"# Generated by STB. Tip position: {r_tip} Ang\n")
        f.write(f"# Broadening: sigma = {sigma_eV * 1000:.3f} meV, shift = {args.shift}\n")
        f.write(f"#{'Energy(eV)':<14}{'dIdV(arb.units)':<18}\n")
        for e, v in zip(E, sts_curve):
            f.write(f"{e:<15.6f}{v:<18.6e}\n")
    print(f"[INFO] Wrote {out_path}")

    if args.plot:
        plt.figure(figsize=(7, 5))
        plt.plot(E, sts_curve, color="#cc5522", lw=2)
        plt.xlabel("Energy (eV)")
        plt.ylabel("dI/dV (arb. units)")
        plt.title(f"Simulated STS at {r_tip} Ang")
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    print("\n[INFO] Complete job!")
    print("\n" + "-" * 60)


if __name__ == "__main__":
    main()
