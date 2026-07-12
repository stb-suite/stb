#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""Inverse Participation Ratio (IPR) scattered onto a band structure: one
value per (k, band), sized/colored on the same k-path plot stb-fatbands
uses for orbital projections -- here the "category" is a single continuous
localization measure instead of an orbital/atom/species weight.

    IPR_q(k, band) = sum_i |psi_i|^(2q)

(sisl's EigenstateElectron.ipr(q), default q=2) is large for a localized
state (weight concentrated on few orbitals) and small for an extended one
(weight spread over many orbitals) -- the standard Anderson-localization
diagnostic. There's no universal absolute scale: the value implicitly
depends on the number of orbitals/cell size, so only compare IPR values
computed on the SAME system/basis, not across different structures.

Needs a real per-orbital radial... actually no -- IPR only needs the
orbital EXPANSION COEFFICIENTS (norm2/hadamard, same as stb-fatbands),
not real-space radial shapes, so it works identically well from a
.HSX-derived Hamiltonian (accurate, overlap-aware) or a .fdf-derived
Geometry approximation (same "hamiltonian_preferred" mode and printed
accuracy warning as stb-fatbands).

Known limitations (mirrors stb-fatbands' own, since this shares its .bands/
.WFSX loading and cross-check machinery):
 - No file metadata ties a given .WFSX to a given .bands file -- same
   k-count/orbital-count guard and --k-tol eigenvalue cross-check as
   stb-fatbands.
 - IPR values are basis-size-dependent, not directly comparable across
   different structures/basis sets.
"""

VERSION = "1.0.0"

import argparse
import os
import numpy as np

from stb.core.cli import color_text, show_intro
from stb.core.siesta_bands import read_data, _is_gamma, shift_bands, cbm_vbm
from stb.core.siesta_wfsx import (
    resolve_wfsx_path, load_parent, _geometry_of, read_wfsx_states,
)
from stb.core.band_scatter import write_scalar_data, write_scalar_gplot, plot_scalar_on_bands


def main():
    parser = argparse.ArgumentParser(
        description="Inverse Participation Ratio (IPR) scattered onto a band structure, from "
                    ".bands + .WFSX.",
        epilog="Example usage:\n"
               "  stb-ipr --label siesta --shift fermi\n"
               "  stb-ipr --label siesta --shift fermi --q 4\n",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("--label", type=str, default=None,
                        help="SIESTA output label. Shorthand for --file <label>.bands, "
                             "auto-detects <label>.bands.WFSX (falling back to <label>.WFSX "
                             "then <label>.selected.WFSX, with a note), and auto-detects "
                             "<label>.HSX (optional, upgrades weight accuracy) then <label>.fdf "
                             "(geometry+basis fallback). Mutually exclusive with "
                             "--file/--wfsx/--hsx-file/--geometry-file.")
    parser.add_argument("--file", dest="input_file", type=str, default=None,
                        help="Explicit path to the .bands file (alternative to --label).")
    parser.add_argument("--wfsx", type=str, default=None,
                        help="Explicit path to the .WFSX file. Required if --file is used "
                             "instead of --label.")
    parser.add_argument("--hsx-file", type=str, default=None,
                        help="Optional explicit .HSX path for physically correct overlap-"
                             "weighted IPR. Falls back to --geometry-file / auto-detected .fdf "
                             "with an accuracy warning if omitted.")
    parser.add_argument("--geometry-file", type=str, default=None,
                        help="Optional explicit .fdf geometry path, used as the fallback "
                             "parent when no .HSX is available. Requires sisl.")

    parser.add_argument("--shift", type=str, choices=["vbm", "cbm", "fermi", "manual"], required=True,
                        help="Reference energy shift, same vocabulary as stb-bands/stb-fatbands.")
    parser.add_argument("--manual-value", type=float,
                        help="Custom energy shift value (required if --shift manual is used).")
    parser.add_argument("--gap-tol", type=float, default=0.01,
                        help="Energy tolerance in eV for VBM/CBM classification when --shift "
                             "is vbm/cbm (default: 0.01).")

    parser.add_argument("--q", type=int, default=2,
                        help="IPR order parameter (default: 2, the standard choice). Passed "
                             "straight to sisl's EigenstateElectron.ipr(q=...).")

    parser.add_argument("--k-tol", type=float, default=0.001,
                        help="Cross-check tolerance (eV) between WFSX-stored eigenvalues and "
                             ".bands eigenvalues at sampled k-points (default: 0.001). Same "
                             "file-order sanity guard as stb-fatbands.")

    parser.add_argument("-o", "--output-dir", type=str, default=".",
                        help="Directory to write ipr.dat and ipr.gplot into (default: current "
                             "directory). Created if it doesn't exist.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-ipr {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.shift == "manual" and args.manual_value is None:
        parser.error("--manual-value is required when --shift is set to 'manual'.")

    if args.label:
        if args.input_file or args.wfsx or args.hsx_file or args.geometry_file:
            parser.error("--label cannot be combined with --file/--wfsx/--hsx-file/--geometry-file.")
        args.input_file = f"{args.label}.bands"
        args.wfsx = resolve_wfsx_path(args.label)
    elif not args.input_file:
        parser.error("one of --label or --file is required.")
    elif not args.wfsx:
        parser.error("--wfsx is required when using --file instead of --label.")

    if not os.path.isfile(args.input_file):
        parser.error(f"'{args.input_file}' not found.")
    if args.wfsx is None or not os.path.isfile(args.wfsx):
        tried = (f"'{args.label}.bands.WFSX', '{args.label}.WFSX' and "
                 f"'{args.label}.selected.WFSX'") if args.label else f"'{args.wfsx}'"
        parser.error(f"No .WFSX file found ({tried}).")

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("IPR (INVERSE PARTICIPATION RATIO):", 'bold'))
    print("-" * 60)

    print("[INFO] Reading band structure data ...")
    fermi_energy, high_sym, dic_bands, nspin = read_data(args.input_file)
    bands_keys = list(dic_bands.keys())

    print("[INFO] Resolving geometry/Hamiltonian ...")
    parent, parent_source, is_approx = load_parent(args.label, args.hsx_file, args.geometry_file)
    if parent is None:
        parser.error("No geometry/Hamiltonian could be resolved -- pass --hsx-file, "
                      "--geometry-file, or use --label with a <label>.HSX/.fdf next to it.")
    if is_approx:
        print(color_text(
            "[WARNING] No .HSX found -- IPR uses an implicit-orthogonal-basis approximation, "
            "not SIESTA's true non-orthogonal overlap-weighted norm. Pass --hsx-file "
            "<label>.HSX for physically correct values.", 'yellow'))
    else:
        print(f"[INFO] Using '{parent_source}' for overlap-aware IPR.")

    geometry = _geometry_of(parent)

    print(f"[INFO] Reading wavefunction coefficients from '{args.wfsx}' ...")
    sizes, states_by_k = read_wfsx_states(args.wfsx, parent)

    if sizes.nk != len(bands_keys):
        parser.error(
            f"k-point count mismatch: '{args.wfsx}' has {sizes.nk} k-points but "
            f"'{args.input_file}' has {len(bands_keys)} -- they must be from the same "
            "calculation/k-path."
        )
    if sizes.no_u != geometry.no:
        parser.error(
            f"Orbital count mismatch: '{args.wfsx}' has {sizes.no_u} orbitals but the "
            f"loaded geometry ('{parent_source}') has {geometry.no} -- they must be from "
            "the same calculation."
        )

    # Sanity cross-check: same file-order-correspondence guard as
    # stb-fatbands (sample first/middle/last k, compare eigenvalues).
    n = len(bands_keys)
    sample_idx = sorted(set([0, n // 2, n - 1]))
    for i in sample_idx:
        k_key = bands_keys[i]
        band_eigs = dic_bands[k_key]
        wfsx_block = states_by_k[i]
        for s in range(nspin):
            state = wfsx_block.get(s)
            if state is None:
                continue
            wfsx_eigs = np.asarray(state.eig)
            m = min(len(wfsx_eigs), band_eigs.shape[1])
            diff = np.abs(wfsx_eigs[:m] - band_eigs[s][:m])
            if diff.size and np.max(diff) > args.k_tol:
                parser.error(
                    f"Eigenvalue mismatch between '{args.wfsx}' and '{args.input_file}' at "
                    f"k-index {i} (spin {s}): max diff {np.max(diff):.6f} eV exceeds --k-tol "
                    f"{args.k_tol} eV -- these files likely don't correspond to the same "
                    "calculation/k-path order."
                )
    print("[INFO] .bands / .WFSX correspondence check passed.")

    result = cbm_vbm(fermi_energy, dic_bands, nspin, args.gap_tol)
    vbm, cbm, _, _, _, _, _, _ = result["combined"]
    if args.shift == "vbm":
        rshift = vbm
    elif args.shift == "cbm":
        rshift = cbm
    elif args.shift == "fermi":
        rshift = fermi_energy
    else:
        rshift = args.manual_value
    shifted_bands = shift_bands(dic_bands, rshift)

    print(f"[INFO] Computing IPR (q={args.q}) for every (k, band, spin) state ...")
    rows = []
    for i, k_key in enumerate(bands_keys):
        wfsx_block = states_by_k[i]
        energies = shifted_bands[k_key]
        for s in range(nspin):
            state = wfsx_block.get(s)
            if state is None:
                continue
            ipr = np.asarray(state.ipr(q=args.q)).real
            nbands_avail = min(len(ipr), energies.shape[1])
            for b in range(nbands_avail):
                rows.append((k_key, float(energies[s][b]), float(ipr[b])))

    os.makedirs(args.output_dir, exist_ok=True)
    data_path = write_scalar_data(args.output_dir, "ipr", rows)
    print(f"[INFO] Wrote {data_path}")
    gplot_path = write_scalar_gplot(args.output_dir, "ipr.gplot", "ipr.pdf", high_sym, ["ipr"])
    print(f"[INFO] Wrote {gplot_path}")

    k_pos, energy, value = zip(*rows)
    plot_scalar_on_bands(high_sym, k_pos, energy, value, f"IPR (q={args.q})", _is_gamma)

    print("\n[INFO] Complete job!")
    print("\n" + "-" * 60)


if __name__ == "__main__":
    main()
