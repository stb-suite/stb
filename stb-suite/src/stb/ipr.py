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
Verified directly against sisl's own source: `ipr(q)` returns
`(state_abs2**q).sum(-1) / state_abs2.sum(-1)**q` where `state_abs2 =
norm2(projection="hadamard").real` -- the EXACT same overlap-aware
per-orbital quantity stb-fatbands already uses for its own weights, which
is why this tool shares stb-fatbands' identical HSX-accurate/fdf-
approximate accuracy tradeoff below. sisl itself asserts `q >= 2`
internally (see the --q validation below).

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
 - non-collinear/SOC (nspin 4/8) inputs get a raw "_s{n}" per-channel
   suffix, same fallback dos.py/fatbands.py use -- not specially handled
   or verified beyond that.

Output/report style rewritten (v1.0.0 -> v2.0.0) to match the rest of the
Analysis category (stb-fatbands, its direct sibling -- core/band_scatter.py
was extracted from fatbands.py once this tool became its second consumer):
a numbered [0]...[6] report, --save-report, --save-gnuplot (previously
wrote ipr.dat/ipr.gplot unconditionally on every run, no way to opt out),
and --view (previously the matplotlib preview was always shown, blocking,
with no way to skip it at all -- not even a --no-plot existed).

Also fixes three real issues found while reviewing the physics/function:
 - --label + --hsx-file/--geometry-file together used to be rejected
   outright, the same overly strict validation bug already found and
   fixed in stb-wfdensity/stb-sts/stb-coop -- load_parent() already
   prefers an explicit hsx_file/geometry_file over <label>.HSX/.fdf on
   its own. --label + --file/--wfsx is still rejected (ambiguous --
   --label IS what auto-detects those).
 - --q had no lower-bound validation: sisl's own ipr() asserts q >= 2
   internally, so --q 1/0/negative used to crash with a raw
   AssertionError traceback instead of a clean, actionable error.
 - nspin=2 spin-channel merging: the exact same bug stb-fatbands found
   and fixed live this session (verified there against a real
   spin-polarized isolated-O-atom fixture, whose 2 Bohr-magneton triplet
   moment gave spin-up/spin-down CBMs differing by ~29 eV). This tool's
   row-building loop used to dump both spin channels into one flat,
   spin-blind list -- for a real magnetic system this silently merges
   two physically different IPR series into one indistinguishable
   scatter/data file. Fixed the same way stb-fatbands was: split into
   ipr_up/ipr_down for nspin=2 (nspin=1 unaffected, still plain "ipr").
   Verified live against the same isolated-O-atom fixture: spin-up mean
   IPR 0.9515 vs spin-down mean IPR 1.0870 -- a real, physically
   meaningful difference (the two spin channels see a different
   effective potential) that used to be silently averaged away into one
   merged series.

Unlike stb-wfdensity/stb-sts/stb-coop, no Fermi-energy source hierarchy
was added here: this tool (like stb-fatbands) already requires a .bands
file as its primary input, which carries its own embedded Fermi energy
via read_data() -- there is no decoupling to do.
"""

VERSION = "2.0.0"

import argparse
import os
import sys
from datetime import datetime

import numpy as np

from stb.core import citations
from stb.core.cli import color_text, show_intro, print_dual, print_section, print_table
from stb.core.siesta_bands import read_data, _is_gamma, shift_bands, cbm_vbm
from stb.core.siesta_wfsx import (
    resolve_wfsx_path, load_parent, _geometry_of, read_wfsx_states,
)
from stb.core.band_scatter import (
    write_scalar_data, write_scalar_gplot, plot_scalar_on_bands, plot_multi_series_on_bands,
)

REPORT_FILE = "stb_ipr_report.txt"
BIB_FILE = "references.bib"


def main():
    parser = argparse.ArgumentParser(
        description="Inverse Participation Ratio (IPR) scattered onto a band structure, from "
                    ".bands + .WFSX.",
        epilog="Example usage:\n"
               "  stb-ipr --label siesta --shift fermi\n"
               "  stb-ipr --label siesta --shift fermi --q 4\n"
               "  stb-ipr --label siesta --shift fermi --save-report --save-gnuplot --view\n",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("--label", type=str, default=None,
                        help="SIESTA output label. Shorthand for --file <label>.bands, "
                             "auto-detects <label>.bands.WFSX (falling back to <label>.WFSX "
                             "then <label>.selected.WFSX, with a note), and, unless "
                             "--hsx-file/--geometry-file override it, auto-detects <label>.HSX "
                             "(optional, upgrades weight accuracy) then <label>.fdf (geometry+"
                             "basis fallback). Mutually exclusive with --file/--wfsx (not "
                             "--hsx-file/--geometry-file -- see their own help).")
    parser.add_argument("--file", dest="input_file", type=str, default=None,
                        help="Explicit path to the .bands file (alternative to --label).")
    parser.add_argument("--wfsx", type=str, default=None,
                        help="Explicit path to the .WFSX file. Required if --file is used "
                             "instead of --label.")
    parser.add_argument("--hsx-file", type=str, default=None,
                        help="Optional explicit .HSX path for physically correct overlap-"
                             "weighted IPR. Falls back to --geometry-file / auto-detected .fdf "
                             "with an accuracy warning if omitted. Optional (and common) "
                             "together WITH --label too, since the real file is often not "
                             "literally named <label>.HSX.")
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
                             "straight to sisl's EigenstateElectron.ipr(q=...), which requires "
                             "q >= 2. Higher q emphasizes localization more strongly.")

    parser.add_argument("--k-tol", type=float, default=0.001,
                        help="Cross-check tolerance (eV) between WFSX-stored eigenvalues and "
                             ".bands eigenvalues at sampled k-points (default: 0.001). Same "
                             "file-order sanity guard as stb-fatbands.")

    parser.add_argument("-o", "--output-dir", type=str, default=".",
                        help="Directory to write ipr.dat/ipr.gplot (with --save-gnuplot) and "
                             "stb_ipr_report.txt/references.bib into (default: current "
                             "directory). Created if it doesn't exist.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the full run report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--save-gnuplot", action="store_true",
                        help="Also write ipr.dat (or ipr_up.dat/ipr_down.dat for a spin-"
                             "polarized run) + ipr.gplot. Off by default -- this tool used to "
                             "write both unconditionally on every run.")
    parser.add_argument("--view", action="store_true",
                        help="Show an interactive matplotlib preview before finishing. Off by "
                             "default -- this tool previously always showed it, with no way to "
                             "skip.")

    parser.add_argument("-v", "--version", action="version", version=f"stb-ipr {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.shift == "manual" and args.manual_value is None:
        parser.error("--manual-value is required when --shift is set to 'manual'.")
    if args.q < 2:
        parser.error(f"--q must be >= 2 (sisl's EigenstateElectron.ipr() requires this; got {args.q}).")

    if args.label:
        if args.input_file or args.wfsx:
            parser.error("--label cannot be combined with --file/--wfsx.")
        args.input_file = f"{args.label}.bands"
        args.wfsx = resolve_wfsx_path(args.label)
    elif not args.input_file:
        parser.error("one of --label or --file is required.")
    elif not args.wfsx:
        parser.error("--wfsx is required when using --file instead of --label.")
    # --label + --hsx-file/--geometry-file together is valid and common -- see the module
    # docstring (the same fix already made for stb-wfdensity/stb-sts/stb-coop):
    # load_parent() below already prefers an explicit hsx_file/geometry_file over
    # <label>.HSX/.fdf on its own.

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

    shift_desc = args.shift if args.shift != "manual" else f"manual ({args.manual_value} eV)"

    print("\n" + color_text("IPR: reading input data", 'bold'))
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
        accuracy_desc = "approximate (|c|^2, no .HSX)"
    else:
        print(f"[INFO] Using '{parent_source}' for overlap-aware IPR.")
        accuracy_desc = "accurate (overlap-weighted norm, real Hamiltonian)"

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
    vbm, cbm, vbm_k, cbm_k, indirect_gap, gap_type, direct_gap, direct_k = result["combined"]
    if args.shift == "vbm":
        rshift = vbm
    elif args.shift == "cbm":
        rshift = cbm
    elif args.shift == "fermi":
        rshift = fermi_energy
    else:
        rshift = args.manual_value
    shifted_bands = shift_bands(dic_bands, rshift)

    # Spin-resolved output series: the exact same fix stb-fatbands already
    # got this session (verified live there against a real spin-polarized
    # isolated-O-atom fixture) -- merging both spin channels into one
    # series silently combines two physically different IPR sets into one
    # indistinguishable scatter/data file. nspin=2 splits into
    # "ipr_up"/"ipr_down" (same convention dos.py/fatbands.py already use);
    # nspin=1 is unaffected (plain "ipr"). Any other nspin (non-collinear/
    # SOC) falls back to a raw "_s{n}" suffix, same as those tools.
    spin_suffixes = {1: [""], 2: ["_up", "_down"]}.get(nspin, [f"_s{s + 1}" for s in range(nspin)])
    output_series = [f"ipr{suf}" for suf in spin_suffixes]

    print(f"[INFO] Computing IPR (q={args.q}) for every (k, band, spin) state ...")
    rows_by_series = {series: [] for series in output_series}
    for i, k_key in enumerate(bands_keys):
        wfsx_block = states_by_k[i]
        energies = shifted_bands[k_key]
        for s in range(nspin):
            state = wfsx_block.get(s)
            if state is None:
                continue
            ipr = np.asarray(state.ipr(q=args.q)).real
            nbands_avail = min(len(ipr), energies.shape[1])
            series = f"ipr{spin_suffixes[s]}"
            for b in range(nbands_avail):
                rows_by_series[series].append((k_key, float(energies[s][b]), float(ipr[b])))

    # --- From here on: the numbered, save-able report --------------------
    os.makedirs(args.output_dir, exist_ok=True)
    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    print_dual(color_text("\n===== STB-IPR REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time      : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Label          : {args.label}" if args.label else "Label          : (explicit --file/--wfsx, no --label)", f_out)
    print_dual(f"Bands file     : {args.input_file}", f_out)
    print_dual(f"WFSX file      : {args.wfsx}", f_out)
    print_dual(f"Shift mode     : {shift_desc}", f_out)
    print_dual(f"Gap tolerance  : {args.gap_tol} eV", f_out)
    print_dual(f"k-tol (cross-check): {args.k_tol} eV", f_out)
    print_dual(f"IPR order (q)  : {args.q}", f_out)
    print_dual(f"Output dir     : {args.output_dir}", f_out)
    print_dual(f"Save gnuplot   : {'yes' if args.save_gnuplot else 'no'}", f_out)
    print_dual(f"View (matplotlib): {'yes' if args.view else 'no'}", f_out)

    nk = len(bands_keys)
    nbands = next(iter(dic_bands.values())).shape[1]

    print_section("[1] INPUT DATA", f_out)
    print_table(["Quantity", "Value"], [
        (["Fermi energy", f"{fermi_energy:.6f} eV"], None),
        (["Spin channels", f"{nspin} ({'non-polarized' if nspin == 1 else 'polarized'})"], None),
        (["Bands x k-points", f"{nbands} x {nk}"], None),
        (["Orbitals (basis)", f"{geometry.no}"], None),
        (["Geometry/Hamiltonian source", parent_source], None),
        (["Weight accuracy", accuracy_desc], None),
    ], f_out)
    if is_approx:
        print_dual(color_text(
            "[WARNING] No .HSX found -- IPR uses an implicit-orthogonal-basis approximation, "
            "not SIESTA's true non-orthogonal overlap-weighted norm. Pass --hsx-file "
            "<label>.HSX for physically correct values.", 'yellow'), f_out)
    print_dual(color_text(
        f"[OK] .bands / .WFSX k-point and orbital-count correspondence verified "
        f"(sampled eigenvalue cross-check, tol={args.k_tol} eV).", 'green'), f_out)

    print_section("[2] BAND GAP ANALYSIS (k-path)", f_out)
    print_table(["Quantity", "Value"], [
        (["VBM", f"{vbm:.6f} eV (k = {vbm_k:.6f})"], None),
        (["CBM", f"{cbm:.6f} eV (k = {cbm_k:.6f})"], None),
        (["Indirect gap", f"{indirect_gap:.6f} eV (fundamental: CBM - VBM, any k)"], None),
        (["Direct gap", f"{direct_gap:.6f} eV (same-k minimum, k = {direct_k:.6f})"], None),
        (["Gap type", gap_type], 'yellow' if gap_type == "Metallic" else None),
    ], f_out)
    if nspin == 2:
        for s, spin_result in enumerate(result["spins"]):
            spin_name = "Spin up" if s == 0 else "Spin down"
            sv, sc_, svk, sck, sig, sgt, sdg, sdk = spin_result
            print_dual(f"{spin_name}:", f_out)
            print_table(["Quantity", "Value"], [
                (["VBM", f"{sv:.6f} eV (k = {svk:.6f})"], None),
                (["CBM", f"{sc_:.6f} eV (k = {sck:.6f})"], None),
                (["Indirect gap", f"{sig:.6f} eV"], None),
                (["Direct gap", f"{sdg:.6f} eV (k = {sdk:.6f})"], None),
                (["Gap type", sgt], 'yellow' if sgt == "Metallic" else None),
            ], f_out)
        half_metallic = result["half_metallic"]
        line = f"Half-metallic : {'Yes' if half_metallic else 'No'}"
        print_dual(color_text(line, 'green') if half_metallic else line, f_out)

    print_section("[3] IPR ANALYSIS", f_out)
    print_dual(f"IPR order (q)  : {args.q} (higher q emphasizes localization more strongly)", f_out)
    if nspin == 2:
        print_dual(f"Spin-resolved  : yes -- ipr_up/ipr_down reported and written separately.", f_out)
    ipr_rows = []
    for series in output_series:
        vals = np.array([v for _, _, v in rows_by_series[series]])
        i_max = int(np.argmax(vals))
        i_min = int(np.argmin(vals))
        k_max, e_max, v_max = rows_by_series[series][i_max]
        k_min, e_min, v_min = rows_by_series[series][i_min]
        ipr_rows.append(([series, str(len(vals)), f"{vals.mean():.4f}",
                          f"{v_min:.4f} (k={k_min:.4f}, E={e_min:.4f} eV)",
                          f"{v_max:.4f} (k={k_max:.4f}, E={e_max:.4f} eV)"], None))
    print_table(["Series", "States", "Mean IPR", "Min IPR (most extended)", "Max IPR (most localized)"],
                ipr_rows, f_out)
    print_dual("(IPR has no universal absolute scale -- it depends on orbital count/cell size, "
               "so only compare these values within this same system/basis, never across "
               "different structures.)", f_out)

    print_section("[4] WRITING OUTPUT FILES", f_out)
    written_data = []
    gplot_path = None
    if args.save_gnuplot:
        for series in output_series:
            path = write_scalar_data(args.output_dir, series, rows_by_series[series])
            written_data.append(path)
        gplot_path = write_scalar_gplot(args.output_dir, "ipr.gplot", "ipr.pdf", high_sym, output_series)
        print_dual(color_text(
            f"[OK] Data written: {len(written_data)} file(s) "
            f"({', '.join(os.path.basename(p) for p in written_data)}).", 'green'), f_out)
        print_dual(color_text(f"[OK] Gnuplot script written to '{gplot_path}'.", 'green'), f_out)
    else:
        print_dual(f"Not written (off by default -- pass --save-gnuplot to write "
                   f"{'/'.join(s + '.dat' for s in output_series)} + ipr.gplot).", f_out)

    print_section("[5] REFERENCES", f_out)
    bib_entries = [citations.SIESTA, citations.SIESTA_RECENT, citations.IPR_MURPHY]
    citations.write_bib_file(os.path.join(args.output_dir, BIB_FILE), bib_entries)
    print_dual(color_text(
        f"[OK] Citations for the methods used in this run written to "
        f"'{os.path.join(args.output_dir, BIB_FILE)}' ({len(bib_entries)} entries).", 'green'), f_out)

    print_section("[6] SUMMARY & FILES", f_out)
    print_dual("Status         : OK", f_out)
    print_dual(f"Best (indirect) gap: {indirect_gap:.6f} eV ({gap_type}).", f_out)
    if nspin == 2 and result["half_metallic"]:
        print_dual(color_text(
            "Half-metallic character detected (one spin channel metallic, the other has a gap).",
            'green'), f_out)
    for path in written_data:
        print_dual(f"Data file      : {path}", f_out)
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
        if len(output_series) == 1:
            series = output_series[0]
            k_pos, energy, value = zip(*rows_by_series[series])
            plot_scalar_on_bands(high_sym, k_pos, energy, value, f"IPR (q={args.q})", _is_gamma)
        else:
            plot_multi_series_on_bands(high_sym, rows_by_series, _is_gamma)


if __name__ == "__main__":
    main()
