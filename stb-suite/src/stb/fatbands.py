#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""Orbital-projected ("fatbands") band structure: combines a SIESTA .bands
file (k-path, eigenvalues) with a .WFSX wavefunction-coefficient file to
color/size each (k, band) point by its orbital, atom, or species character.

Needs sisl (like stb-bands' --eig-file path), and a .WFSX file written
alongside the same calculation that wrote .bands (`WFS.Write.For.Bands T`
in the fdf -- SIESTA/sisl's own convention names it <label>.bands.WFSX, not
<label>.WFSX). If a companion .HSX (Hamiltonian + overlap) is also
available, weights are computed with SIESTA's true non-orthogonal-basis
overlap matrix (physically correct); otherwise they fall back to an
implicit-orthogonal approximation (raw |c|^2 per orbital) with a printed
warning, since sisl's norm2() silently uses an identity overlap whenever
its `parent` is a bare Geometry rather than a Hamiltonian.

Known limitations (mirrors dos.py's own documented orbital-naming caveats,
since both tools share stb.core.orbitals):
 - l > 3 (g-orbitals and beyond) are excluded from --projection l/ml
   (counted and warned about once), but still contribute fully to
   --projection atom/species, which don't depend on l.
 - non-collinear/SOC (nspin 4/8) inputs are not specially handled beyond
   whatever sisl's own norm2()/Sk() already account for internally; spin
   labeling stays collinear-only (same gap dos.py documents for PDOS.xml).
 - There is no file metadata tying a given .WFSX to a given .bands file.
   This tool trusts file-order correspondence, validated only by a
   k-point/orbital-count guard and a --k-tol eigenvalue cross-check at a
   few sampled k-points -- always let --label auto-detect paired files from
   the same run rather than mixing --file/--wfsx from different
   calculations.
"""

VERSION = "1.0.0"

import argparse
import os
import sys
import numpy as np

from stb.core.cli import color_text, show_intro
from stb.core.siesta_bands import read_data, _is_gamma, shift_bands, cbm_vbm
from stb.core.orbitals import ORBITAL_ORDER, get_orbital_name, get_detailed_orbital_name
from stb.core.siesta_wfsx import (
    resolve_wfsx_path, load_parent, _geometry_of, read_wfsx_states,
)
from stb.core.band_scatter import (
    write_scalar_data, write_scalar_gplot, plot_scalar_on_bands, plot_multi_series_on_bands,
)

CATEGORY_CHOICES = ("l", "ml", "atom", "species")


def build_orbital_table(geometry, projection):
    """Per-global-orbital-index parallel lists: category name (per
    `projection`), atom index, species symbol. Orbitals with l > 3 get
    category=None for l/ml (excluded from those sums), but atom index/
    species are always recorded since --projection atom/species don't
    depend on l at all."""
    categories = []
    atoms = []
    species = []
    n_excluded = 0
    excluded_l = set()
    for ia, io in geometry.iter_orbitals(local=True):
        atom = geometry.atoms[ia]
        orb = atom.orbitals[io]
        atoms.append(ia)
        species.append(atom.symbol)
        if projection == "atom":
            categories.append(str(ia))
        elif projection == "species":
            categories.append(atom.symbol)
        elif projection == "l":
            name = get_orbital_name(orb.l)
            if name is None:
                n_excluded += 1
                excluded_l.add(orb.l)
            categories.append(name)
        else:  # ml
            name = get_detailed_orbital_name(orb.l, orb.m)
            if name is None:
                n_excluded += 1
                excluded_l.add(orb.l)
            categories.append(name)

    if n_excluded and projection in ("l", "ml"):
        print(f"[WARNING] {n_excluded} orbital(s) with l={sorted(excluded_l)} "
              "(l > 3, i.e. g-orbitals or beyond, are not supported) excluded "
              f"from --projection {projection} -- still counted in --projection "
              "atom/species, which don't depend on l.", file=sys.stderr)

    return categories, atoms, species


def _category_sort_key(projection):
    if projection in ("l", "ml"):
        def key(c):
            return (ORBITAL_ORDER.index(c) if c in ORBITAL_ORDER else 999, c)
        return key
    if projection == "atom":
        return lambda c: int(c)
    return lambda c: c  # species: plain alphabetical




def main():
    parser = argparse.ArgumentParser(
        description="Orbital-projected (fatbands) band structure from .bands + .WFSX.",
        epilog="Example usage:\n"
               "  stb-fatbands --label siesta --shift fermi --projection l\n"
               "  stb-fatbands --label siesta --shift fermi --projection species --category O Sn\n"
               "  stb-fatbands --file siesta.bands --wfsx siesta.bands.WFSX --hsx-file siesta.HSX --shift fermi",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("--label", type=str, default=None,
                        help="SIESTA output label. Shorthand for --file <label>.bands, "
                             "auto-detects <label>.bands.WFSX (falling back to <label>.WFSX "
                             "with a note), and auto-detects <label>.HSX (optional, upgrades "
                             "weight accuracy) then <label>.fdf (geometry+basis fallback -- NOT "
                             "<label>.XV, which has no per-orbital basis information). "
                             "Mutually exclusive with --file/--wfsx/--hsx-file/--geometry-file.")

    parser.add_argument("--file", dest="input_file", type=str, default=None,
                        help="Explicit path to the .bands file (alternative to --label).")

    parser.add_argument("--wfsx", type=str, default=None,
                        help="Explicit path to the .WFSX file. Required if --file is used "
                             "instead of --label.")

    parser.add_argument("--hsx-file", type=str, default=None,
                        help="Optional explicit .HSX path (Hamiltonian+overlap) for physically "
                             "correct overlap-weighted norms. Falls back to --geometry-file / "
                             "auto-detected .fdf with an accuracy warning if omitted.")

    parser.add_argument("--geometry-file", type=str, default=None,
                        help="Optional explicit .fdf geometry path (must have the basis's "
                             ".ion/.ion.xml files alongside it), used as the fallback parent "
                             "when no .HSX is available. A bare .XV has no per-orbital basis "
                             "information and will fail the orbital-count check against the "
                             ".WFSX -- use --hsx-file or a .fdf instead. Requires sisl.")

    parser.add_argument("--shift", type=str, choices=["vbm", "cbm", "fermi", "manual"], required=True,
                        help="Reference energy shift, same vocabulary as stb-bands:\n"
                             "  - 'vbm'    : Valence Band Maximum\n"
                             "  - 'cbm'    : Conduction Band Minimum\n"
                             "  - 'fermi'  : Fermi level\n"
                             "  - 'manual' : Custom shift value (requires --manual-value).")

    parser.add_argument("--manual-value", type=float,
                        help="Custom energy shift value (required if --shift manual is used).")

    parser.add_argument("--gap-tol", type=float, default=0.01,
                        help="Energy tolerance in eV for VBM/CBM classification when --shift "
                             "is vbm/cbm (default: 0.01). Same meaning as stb-bands.")

    parser.add_argument("--projection", type=str, choices=CATEGORY_CHOICES, default="l",
                        help="Category to color/size markers by (default: l):\n"
                             "  l:       aggregate by angular momentum (s, p, d, f)\n"
                             "  ml:      detailed orbital (s, px, py, pz, dxy, ...)\n"
                             "  atom:    per atom index\n"
                             "  species: per chemical species\n"
                             "Single choice only -- unlike stb-dos's --type, the plot needs one "
                             "weight per marker; run again for another projection.")

    parser.add_argument("--category", type=str, nargs="+", default=None,
                        help="Restrict output to specific category values (e.g. --projection "
                             "species --category O Sn). Default: all categories found.")

    parser.add_argument("--k-tol", type=float, default=0.001,
                        help="Cross-check tolerance (eV) between WFSX-stored eigenvalues and "
                             ".bands eigenvalues at sampled k-points (default: 0.001). This is a "
                             "file-order sanity guard, not a physics tolerance like --gap-tol -- "
                             "raising it only weakens the check that .bands/.WFSX truly "
                             "correspond to the same calculation.")

    parser.add_argument("-o", "--output-dir", type=str, default=".",
                        help="Directory to write fatbands_<category>.dat and fatbands.gplot "
                             "into (default: current directory). Created if it doesn't exist.")

    parser.add_argument("-v", "--version", action="version", version=f"stb-fatbands {VERSION}")

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
        tried = f"'{args.label}.bands.WFSX' and '{args.label}.WFSX'" if args.label else f"'{args.wfsx}'"
        parser.error(f"No .WFSX file found ({tried}). stb-fatbands needs the wavefunction "
                     "coefficients written for the band-structure k-path (SIESTA fdf option "
                     "'WFS.Write.For.Bands T').")

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("FATBANDS:", 'bold'))
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
            "[WARNING] No .HSX found -- orbital weights use an implicit-orthogonal-basis "
            "approximation (|c|^2), not SIESTA's true non-orthogonal overlap-weighted norm. "
            "Pass --hsx-file <label>.HSX for physically correct weights.", 'yellow'))
    else:
        print(f"[INFO] Using '{parent_source}' for overlap-aware orbital weights.")

    geometry = _geometry_of(parent)

    print(f"[INFO] Reading wavefunction coefficients from '{args.wfsx}' ...")
    sizes, states_by_k = read_wfsx_states(args.wfsx, parent)

    if sizes.nk != len(bands_keys):
        parser.error(
            f"k-point count mismatch: '{args.wfsx}' has {sizes.nk} k-points but "
            f"'{args.input_file}' has {len(bands_keys)} -- they must be from the same "
            "calculation/k-path (see stb-fatbands --help for the k-order assumption this "
            "tool relies on)."
        )
    if sizes.no_u != geometry.no:
        parser.error(
            f"Orbital count mismatch: '{args.wfsx}' has {sizes.no_u} orbitals but the "
            f"loaded geometry ('{parent_source}') has {geometry.no} -- they must be from "
            "the same calculation."
        )

    # Sanity cross-check: sample first/middle/last k, compare WFSX-stored
    # eigenvalues against .bands' own -- a hard error (not a warning) on
    # mismatch, since a silently-mis-zipped plot would look plausible but
    # be wrong.
    n = len(bands_keys)
    sample_idx = sorted(set([0, n // 2, n - 1]))
    for i in sample_idx:
        k_key = bands_keys[i]
        band_eigs = dic_bands[k_key]  # shape (nspin, nbands)
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

    print(f"[INFO] Categorizing {geometry.no} orbitals by --projection {args.projection} ...")
    categories, _atoms, _species = build_orbital_table(geometry, args.projection)
    categories = np.array(categories, dtype=object)

    all_categories = sorted({c for c in categories if c is not None}, key=_category_sort_key(args.projection))
    if args.category:
        missing = set(args.category) - set(all_categories)
        if missing:
            parser.error(f"--category value(s) not found for --projection {args.projection}: "
                         f"{sorted(missing)}. Available: {all_categories}")
        selected_categories = [c for c in all_categories if c in args.category]
    else:
        selected_categories = all_categories

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

    print("[INFO] Computing orbital-projected weights (this may take a moment) ...")
    masks = {cat: (categories == cat) for cat in selected_categories}
    weights_by_category = {cat: [] for cat in selected_categories}

    for i, k_key in enumerate(bands_keys):
        wfsx_block = states_by_k[i]
        energies = shifted_bands[k_key]
        for s in range(nspin):
            state = wfsx_block.get(s)
            if state is None:
                continue
            # norm2(projection="hadamard") is a sesquilinear (basis-resolved
            # inner product) quantity -- real by construction, but sisl
            # returns it as a complex dtype (matching the complex
            # wavefunction coefficients) with a numerically ~0 imaginary
            # part. Take .real explicitly rather than letting float()
            # silently (and noisily -- ComplexWarning) drop it.
            weights = state.norm2(projection="hadamard").real  # (nwf, norb)
            nwf = weights.shape[0]
            nbands_avail = min(nwf, energies.shape[1])
            for cat in selected_categories:
                cat_weight_per_band = weights[:nbands_avail, masks[cat]].sum(axis=1)
                for b in range(nbands_avail):
                    weights_by_category[cat].append(
                        (k_key, float(energies[s][b]), float(cat_weight_per_band[b]))
                    )

    os.makedirs(args.output_dir, exist_ok=True)
    data_names = [f"fatbands_{cat}" for cat in selected_categories]
    for cat, name in zip(selected_categories, data_names):
        path = write_scalar_data(args.output_dir, name, weights_by_category[cat])
        print(f"[INFO] Wrote {path}")
    gplot_path = write_scalar_gplot(args.output_dir, "fatbands.gplot", "fatbands.pdf",
                                    high_sym, data_names)
    print(f"[INFO] Wrote {gplot_path}")

    if len(selected_categories) == 1:
        cat = selected_categories[0]
        k_pos, energy, weight = zip(*weights_by_category[cat])
        plot_scalar_on_bands(high_sym, k_pos, energy, weight, f"{cat} weight", _is_gamma)
    else:
        plot_multi_series_on_bands(high_sym, weights_by_category, _is_gamma)

    print("\n[INFO] Complete job!")
    print("\n" + "-" * 60)


if __name__ == "__main__":
    main()
