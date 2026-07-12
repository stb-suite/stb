#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""Spin texture: <Sx>, <Sy>, <Sz> expectation values scattered onto a band
structure, one point per (k, band), via sisl's EigenstateElectron.
spin_moment(). Only physically meaningful for a non-collinear (nspin=4) or
spin-orbit (nspin=8) SIESTA calculation -- verified live in sisl's own
source: spin_moment() needs Sk() (a real overlap matrix), and the
underlying physics (Pauli-matrix expectation values of a 2-component
spinor) has no meaning for a collinear (nspin=1/2) wavefunction.

Needs the same "hamiltonian_preferred" parent as stb-fatbands/stb-ipr: a
.HSX gives the true overlap-aware spin moment, a .fdf-derived Geometry
falls back to an orthogonal-basis approximation (printed warning).

Unlike stb-fatbands, this tool does NOT require a companion .bands file --
spin_moment() is computed directly from whatever k-points the .WFSX
itself contains (any WFSX: band-path, full-BZ mesh, or even a single
Gamma-only point), and the x-axis is a plain 0-based k-INDEX rather than a
physical k-path arc length (no .bands file means no k-path geometry to
compute that from). For a real SOC band-structure spin-texture plot
(e.g. Rashba splitting, topological surface states), pass a genuine
band-path .WFSX and expect the x-axis to just be evenly spaced by k-index,
not scaled to match stb-bands/stb-fatbands' k-path plots exactly.

Known limitations:
 - Needs nspin=4 (non-collinear) or nspin=8 (spin-orbit); aborts with a
   clear error for any other nspin.
 - Test-fixture caveat (see test/3-analysis/17-spintexture/): the bundled
   fixture is a single isolated O atom run as non-collinear WITHOUT an
   explicit initial spin canting (no %block DM.InitSpin -- its exact
   syntax could not be confirmed reliably during development). It
   converges with a physically sensible, strongly non-zero Sz (~+-1.0,
   as expected for an open-shell atom) and Sx/Sy ~0 (numerical noise,
   ~1e-15 to 1e-30) -- this validates that the numerical pipeline is
   correct, but does NOT demonstrate a genuinely "textured" (canted,
   k-dependent Sx/Sy) spin texture. A real spin-orbit-coupling
   calculation on a suitable heavy-element system would be needed for
   that.
 - Only "diagonal" projection (spin_moment per band), not the orbital-
   resolved "hadamard" projection sisl also supports.
"""

VERSION = "1.0.0"

import argparse
import os
import numpy as np

from stb.core.cli import color_text, show_intro
from stb.core.siesta_bands import _is_gamma, select_band_vbm_cbm
from stb.core.siesta_wfsx import resolve_wfsx_path, load_parent, read_wfsx_states
from stb.core.band_scatter import write_scalar_data, write_scalar_gplot, plot_scalar_on_bands

COMPONENT_NAMES = ("Sx", "Sy", "Sz")


def main():
    parser = argparse.ArgumentParser(
        description="Spin texture (<Sx>,<Sy>,<Sz>) scattered onto a band structure, from a "
                    "non-collinear/SOC SIESTA .WFSX.",
        epilog="Example usage:\n"
               "  stb-spintexture --label siesta --shift fermi --fermi -4.2\n",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("--label", type=str, default=None,
                        help="SIESTA output label. Auto-detects <label>.WFSX, falling back to "
                             "<label>.selected.WFSX then <label>.bands.WFSX, and <label>.HSX "
                             "(optional, upgrades accuracy) then <label>.fdf. Mutually exclusive "
                             "with --wfsx/--hsx-file/--geometry-file.")
    parser.add_argument("--wfsx", type=str, default=None,
                        help="Explicit path to the .WFSX file (alternative to --label).")
    parser.add_argument("--hsx-file", type=str, default=None,
                        help="Optional explicit .HSX path for physically correct overlap-"
                             "weighted spin moments. Falls back to --geometry-file / auto-"
                             "detected .fdf with an accuracy warning if omitted.")
    parser.add_argument("--geometry-file", type=str, default=None,
                        help="Optional explicit .fdf geometry path, used as the fallback "
                             "parent when no .HSX is available.")

    parser.add_argument("--shift", type=str, choices=["vbm", "cbm", "fermi", "manual", "none"],
                        default="none",
                        help="Energy reference (default: none, i.e. raw eigenvalues as stored "
                             "in the .WFSX).")
    parser.add_argument("--manual-value", type=float, default=None,
                        help="Custom shift value (eV), required if --shift manual.")
    parser.add_argument("--fermi", type=float, default=None,
                        help="Fermi energy (eV), required if --shift is fermi/vbm/cbm.")
    parser.add_argument("--gap-tol", type=float, default=0.01,
                        help="Energy tolerance in eV for VBM/CBM classification when --shift "
                             "is vbm/cbm (default: 0.01).")

    parser.add_argument("-o", "--output-dir", type=str, default=".",
                        help="Directory to write spintexture_S{x,y,z}.dat/.gplot into (default: "
                             "current directory). Created if it doesn't exist.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-spintexture {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.shift == "manual" and args.manual_value is None:
        parser.error("--manual-value is required when --shift is 'manual'.")
    if args.shift in ("fermi", "vbm", "cbm") and args.fermi is None:
        parser.error("--fermi is required when --shift is fermi/vbm/cbm.")

    if args.label:
        if args.wfsx or args.hsx_file or args.geometry_file:
            parser.error("--label cannot be combined with --wfsx/--hsx-file/--geometry-file.")
        args.wfsx = resolve_wfsx_path(args.label, suffixes=(".WFSX", ".selected.WFSX", ".bands.WFSX"))
    elif not args.wfsx:
        parser.error("one of --label or --wfsx is required.")

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

    print("\n" + color_text("SPIN TEXTURE:", 'bold'))
    print("-" * 60)

    print("[INFO] Resolving geometry/Hamiltonian ...")
    parent, parent_source, is_approx = load_parent(args.label, args.hsx_file, args.geometry_file)
    if parent is None:
        parser.error("No geometry/Hamiltonian could be resolved -- pass --hsx-file, "
                      "--geometry-file, or use --label with a <label>.HSX/.fdf next to it.")
    if is_approx:
        print(color_text(
            "[WARNING] No .HSX found -- spin moments use an implicit-orthogonal-basis "
            "approximation, not SIESTA's true non-orthogonal overlap-weighted values. Pass "
            "--hsx-file <label>.HSX for physically correct values.", 'yellow'))
    else:
        print(f"[INFO] Using '{parent_source}' for overlap-aware spin moments.")

    print(f"[INFO] Reading wavefunction coefficients from '{args.wfsx}' ...")
    sizes, states_by_k = read_wfsx_states(args.wfsx, parent)

    if sizes.nspin not in (4, 8):
        parser.error(
            f"spin_moment() only makes physical sense for a non-collinear (nspin=4) or "
            f"spin-orbit (nspin=8) calculation -- '{args.wfsx}' has nspin={sizes.nspin}."
        )
    print(f"[INFO] nspin={sizes.nspin} ({'non-collinear' if sizes.nspin == 4 else 'spin-orbit'}), "
          f"{sizes.nk} k-point(s).")

    rshift = 0.0
    if args.shift in ("vbm", "cbm"):
        print(f"[INFO] Searching whole k-mesh for the global {args.shift.upper()} ...")
        k_index, spin, band_index = select_band_vbm_cbm(states_by_k, 1, args.fermi, args.gap_tol, args.shift)
        target_state = states_by_k[k_index][spin]
        rshift = float(np.asarray(target_state.eig)[band_index])
        print(f"[INFO] {args.shift.upper()} = {rshift:.6f} eV.")
    elif args.shift == "fermi":
        rshift = args.fermi
    elif args.shift == "manual":
        rshift = args.manual_value

    print("[INFO] Computing spin moments for every (k, band) state ...")
    rows = {name: [] for name in COMPONENT_NAMES}
    for k_index, block in enumerate(states_by_k):
        state = block[0]
        eigs = np.asarray(state.eig) - rshift
        sm = state.spin_moment().real  # (3, nbands)
        nbands = min(sm.shape[1], len(eigs))
        for b in range(nbands):
            for i, name in enumerate(COMPONENT_NAMES):
                rows[name].append((float(k_index), float(eigs[b]), float(sm[i, b])))

    os.makedirs(args.output_dir, exist_ok=True)
    names = [f"spintexture_{name}" for name in COMPONENT_NAMES]
    for name, key in zip(names, COMPONENT_NAMES):
        path = write_scalar_data(args.output_dir, name, rows[key])
        print(f"[INFO] Wrote {path}")
    # A synthetic "high_sym" with just the k-index range, so the shared
    # gnuplot/matplotlib helpers (built for a real .bands k-path) still
    # work with a plain k-index x-axis -- no .bands file is needed here.
    nk = len(states_by_k)
    high_sym = [["0", "k=0"], [str(nk - 1), f"k={nk - 1}"]] if nk > 1 else [["0", "k=0"], ["0", "k=0"]]
    gplot_path = write_scalar_gplot(args.output_dir, "spintexture.gplot", "spintexture.pdf",
                                    high_sym, names)
    print(f"[INFO] Wrote {gplot_path}")

    for key in COMPONENT_NAMES:
        k_pos, energy, value = zip(*rows[key])
        plot_scalar_on_bands(high_sym, k_pos, energy, value, key, _is_gamma, is_signed=True)

    print("\n[INFO] Complete job!")
    print("\n" + "-" * 60)


if __name__ == "__main__":
    main()
