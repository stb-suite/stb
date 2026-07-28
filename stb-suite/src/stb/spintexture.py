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

Physics verified live: for a normalized 2-component spinor, the Pauli
expectation value magnitude |S| = sqrt(Sx^2+Sy^2+Sz^2) is bounded by 1 --
confirmed directly on this suite's own IsolatedO fixture (max |S| =
1.000001, i.e. 1 to within numerical noise). This tool now reports this
normalization check explicitly (same "verify, don't just assert" pattern
already used by stb-wfdensity's own |psi|^2 normalization check).

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

Output/report style rewritten (v1.0.0 -> v2.0.0) to match the rest of the
Analysis category (stb-wfdensity/stb-sts/stb-coop/stb-ipr/stb-effmass): a
numbered [0]...[6] report, --save-report, --save-gnuplot (previously
wrote spintexture_S{x,y,z}.dat/.gplot unconditionally on every run, no
way to opt out), and --view (previously the matplotlib preview was
always shown, blocking, with no way to skip it at all).

Also fixes two real issues found while reviewing the physics/function:
 - --label + --hsx-file/--geometry-file together used to be rejected
   outright, the same overly strict validation bug already found and
   fixed in stb-wfdensity/stb-sts/stb-coop/stb-effmass -- load_parent()
   already prefers an explicit hsx_file/geometry_file over
   <label>.HSX/.fdf on its own. --label + --wfsx is still rejected
   (ambiguous -- --label IS what auto-detects that).
 - --shift fermi/vbm/cbm's Fermi-energy source decoupled from --label,
   the same priority-ordered hierarchy stb-effmass (--band vbm/cbm) has:
   --fermi (explicit) > --bands-file > --fermi-file > an auto-detected
   .out log, via the shared core.siesta_bands.resolve_fermi_energy_hierarchy
   -- this tool previously only ever accepted an explicit --fermi value.
"""

VERSION = "2.0.0"

import argparse
import os
from datetime import datetime

import numpy as np

from stb.core import citations
from stb.core.cli import color_text, show_intro, print_dual, print_section, print_table
from stb.core.siesta_bands import _is_gamma, select_band_vbm_cbm, resolve_fermi_energy_hierarchy
from stb.core.siesta_wfsx import resolve_wfsx_path, load_parent, read_wfsx_states
from stb.core.band_scatter import write_scalar_data, write_scalar_gplot, plot_scalar_on_bands

REPORT_FILE = "stb_spintexture_report.txt"
BIB_FILE = "references.bib"

COMPONENT_NAMES = ("Sx", "Sy", "Sz")


def main():
    parser = argparse.ArgumentParser(
        description="Spin texture (<Sx>,<Sy>,<Sz>) scattered onto a band structure, from a "
                    "non-collinear/SOC SIESTA .WFSX.",
        epilog="Example usage:\n"
               "  stb-spintexture --label siesta --shift fermi --fermi -4.2\n"
               "  stb-spintexture --label siesta --shift fermi --fermi-file calc.out \\\n"
               "      --save-report --save-gnuplot --view\n",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("--label", type=str, default=None,
                        help="SIESTA output label. Auto-detects <label>.WFSX, falling back to "
                             "<label>.selected.WFSX then <label>.bands.WFSX, and, unless "
                             "--hsx-file/--geometry-file override it, <label>.HSX (optional, "
                             "upgrades accuracy) then <label>.fdf. Mutually exclusive with "
                             "--wfsx (not --hsx-file/--geometry-file -- see their own help).")
    parser.add_argument("--wfsx", type=str, default=None,
                        help="Explicit path to the .WFSX file (alternative to --label).")
    parser.add_argument("--hsx-file", type=str, default=None,
                        help="Optional explicit .HSX path for physically correct overlap-"
                             "weighted spin moments. Falls back to --geometry-file / auto-"
                             "detected .fdf with an accuracy warning if omitted. Optional (and "
                             "common) together WITH --label too, since the real file is often "
                             "not literally named <label>.HSX.")
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
                        help="Fermi energy (eV) for --shift fermi/vbm/cbm. Highest-priority "
                             "source if given.")
    parser.add_argument("--bands-file", type=str, default=None,
                        help="Companion .bands file to read the Fermi energy from, for --shift "
                             "fermi/vbm/cbm.")
    parser.add_argument("--fermi-file", type=str, default=None,
                        help="Explicit SIESTA .out log to read the Fermi energy from, for "
                             "--shift fermi/vbm/cbm -- an alternative to --bands-file for a run "
                             "with no saved .bands file. Not assumed to be named after --label.")
    parser.add_argument("--gap-tol", type=float, default=0.01,
                        help="Energy tolerance in eV for VBM/CBM classification when --shift "
                             "is vbm/cbm (default: 0.01).")

    parser.add_argument("-o", "--output-dir", type=str, default=".",
                        help="Directory to write spintexture_S{x,y,z}.dat/.gplot (with "
                             "--save-gnuplot) and stb_spintexture_report.txt/references.bib "
                             "into (default: current directory). Created if it doesn't exist.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the full run report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--save-gnuplot", action="store_true",
                        help="Also write spintexture_S{x,y,z}.dat + spintexture.gplot. Off by "
                             "default -- this tool used to write both unconditionally on every "
                             "run.")
    parser.add_argument("--view", action="store_true",
                        help="Show an interactive matplotlib preview before finishing. Off by "
                             "default -- this tool previously always showed it, with no way to "
                             "skip.")

    parser.add_argument("-v", "--version", action="version", version=f"stb-spintexture {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.shift == "manual" and args.manual_value is None:
        parser.error("--manual-value is required when --shift is 'manual'.")

    if args.label:
        if args.wfsx:
            parser.error("--label cannot be combined with --wfsx.")
        args.wfsx = resolve_wfsx_path(args.label, suffixes=(".WFSX", ".selected.WFSX", ".bands.WFSX"))
    elif not args.wfsx:
        parser.error("one of --label or --wfsx is required.")
    # --label + --hsx-file/--geometry-file together is valid and common -- see the module
    # docstring (the same fix already made for stb-wfdensity/stb-sts/stb-coop/stb-effmass):
    # load_parent() below already prefers an explicit hsx_file/geometry_file over
    # <label>.HSX/.fdf on its own.

    if args.wfsx is None or not os.path.isfile(args.wfsx):
        tried = (f"'{args.label}.WFSX', '{args.label}.selected.WFSX' and "
                 f"'{args.label}.bands.WFSX'") if args.label else f"'{args.wfsx}'"
        parser.error(f"No .WFSX file found ({tried}).")

    fermi_energy = fermi_source = None
    if args.shift in ("fermi", "vbm", "cbm"):
        fermi_energy, fermi_source = resolve_fermi_energy_hierarchy(
            args.fermi, args.bands_file, args.fermi_file, args.label)
        if fermi_energy is None:
            parser.error(
                f"--shift {args.shift} needs a Fermi energy -- none found. Pass --fermi <value>, "
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

    print("\n" + color_text("SPINTEXTURE: reading input data", 'bold'))
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
        accuracy_desc = "approximate (implicit-orthogonal, no .HSX)"
    else:
        print(f"[INFO] Using '{parent_source}' for overlap-aware spin moments.")
        accuracy_desc = "accurate (overlap-weighted, real Hamiltonian)"

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
    vbm_cbm_desc = None
    shift_desc = args.shift if args.shift != "manual" else f"manual ({args.manual_value} eV)"
    if args.shift in ("vbm", "cbm"):
        print(f"[INFO] Searching whole k-mesh for the global {args.shift.upper()} "
              f"(Fermi = {fermi_energy:.6f} eV, source: {fermi_source}) ...")
        k_index, spin, band_index = select_band_vbm_cbm(states_by_k, 1, fermi_energy, args.gap_tol, args.shift)
        target_state = states_by_k[k_index][spin]
        rshift = float(np.asarray(target_state.eig)[band_index])
        print(f"[INFO] {args.shift.upper()} = {rshift:.6f} eV.")
        vbm_cbm_desc = (f"{args.shift.upper()} = {rshift:.6f} eV, found at k-index {k_index}, "
                        f"band {band_index + 1} (Fermi = {fermi_energy:.6f} eV, source: {fermi_source})")
    elif args.shift == "fermi":
        rshift = fermi_energy
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

    # --- From here on: the numbered, save-able report --------------------
    os.makedirs(args.output_dir, exist_ok=True)
    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    print_dual(color_text("\n===== STB-SPINTEXTURE REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time      : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Label          : {args.label}" if args.label else "Label          : (explicit --wfsx)", f_out)
    print_dual(f"WFSX file      : {args.wfsx}", f_out)
    print_dual(f"Hamiltonian src: {parent_source}", f_out)
    print_dual(f"Shift mode     : {shift_desc}" + (f" (Fermi source: {fermi_source})" if args.shift in ("fermi", "vbm", "cbm") else ""), f_out)
    print_dual(f"Output dir     : {args.output_dir}", f_out)
    print_dual(f"Save gnuplot   : {'yes' if args.save_gnuplot else 'no'}", f_out)
    print_dual(f"View (matplotlib): {'yes' if args.view else 'no'}", f_out)

    print_section("[1] INPUT DATA", f_out)
    print_table(["Quantity", "Value"], [
        (["k-points in .WFSX", f"{sizes.nk}"], None),
        (["Spin channels", f"{sizes.nspin} ({'non-collinear' if sizes.nspin == 4 else 'spin-orbit'})"], None),
        (["Weight accuracy", accuracy_desc], None),
    ], f_out)
    if is_approx:
        print_dual(color_text(
            "[WARNING] No .HSX found -- spin moments use an implicit-orthogonal-basis "
            "approximation, not SIESTA's true non-orthogonal overlap-weighted values. Pass "
            "--hsx-file <label>.HSX for physically correct values.", 'yellow'), f_out)

    print_section("[2] ENERGY REFERENCE", f_out)
    print_dual(f"Shift mode: {shift_desc}", f_out)
    if vbm_cbm_desc:
        print_dual(vbm_cbm_desc, f_out)
    elif args.shift == "fermi":
        print_dual(f"Fermi = {fermi_energy:.6f} eV, source: {fermi_source}.", f_out)

    print_section("[3] SPIN TEXTURE ANALYSIS", f_out)
    sx = np.array([v for _, _, v in rows["Sx"]])
    sy = np.array([v for _, _, v in rows["Sy"]])
    sz = np.array([v for _, _, v in rows["Sz"]])
    S_mag = np.sqrt(sx**2 + sy**2 + sz**2)
    print_table(["Component", "Mean", "Min", "Max"], [
        (["Sx", f"{sx.mean():.6f}", f"{sx.min():.6f}", f"{sx.max():.6f}"], None),
        (["Sy", f"{sy.mean():.6f}", f"{sy.min():.6f}", f"{sy.max():.6f}"], None),
        (["Sz", f"{sz.mean():.6f}", f"{sz.min():.6f}", f"{sz.max():.6f}"], None),
    ], f_out)
    print_dual(f"|S| (Pauli expectation magnitude): mean={S_mag.mean():.6f}, max={S_mag.max():.6f}", f_out)
    if S_mag.max() > 1.02:
        print_dual(color_text(
            "[WARNING] |S| exceeds 1 by more than 2% for at least one state -- this should "
            "never happen for a properly normalized spinor; check the .WFSX/.HSX correspondence.",
            'yellow'), f_out)
    else:
        print_dual(color_text(
            "[OK] |S| stays within the physical bound of 1 (to numerical noise), as expected "
            "for a normalized 2-component spinor.", 'green'), f_out)
    i_max = int(np.argmax(S_mag))
    k_max, e_max, _ = rows["Sz"][i_max]
    print_dual(f"Most spin-polarized state: |S|={S_mag[i_max]:.6f} at k-index={int(k_max)}, "
               f"E={e_max:.6f} eV (Sx={sx[i_max]:.4f}, Sy={sy[i_max]:.4f}, Sz={sz[i_max]:.4f}).", f_out)

    print_section("[4] OUTPUT DATA & PLOTS", f_out)
    names = [f"spintexture_{name}" for name in COMPONENT_NAMES]
    nk = len(states_by_k)
    high_sym = [["0", "k=0"], [str(nk - 1), f"k={nk - 1}"]] if nk > 1 else [["0", "k=0"], ["0", "k=0"]]
    written_data = []
    gplot_path = None
    if args.save_gnuplot:
        for name, key in zip(names, COMPONENT_NAMES):
            path = write_scalar_data(args.output_dir, name, rows[key])
            written_data.append(path)
        gplot_path = write_scalar_gplot(args.output_dir, "spintexture.gplot", "spintexture.pdf",
                                        high_sym, names)
        print_dual(color_text(
            f"[OK] Data written: {', '.join(os.path.basename(p) for p in written_data)}.", 'green'), f_out)
        print_dual(color_text(f"[OK] Gnuplot script written to '{gplot_path}'.", 'green'), f_out)
    else:
        print_dual("Not written (off by default -- pass --save-gnuplot to write "
                   "spintexture_S{x,y,z}.dat + spintexture.gplot).", f_out)

    print_section("[5] REFERENCES", f_out)
    bib_entries = [citations.SIESTA, citations.SIESTA_RECENT]
    citations.write_bib_file(os.path.join(args.output_dir, BIB_FILE), bib_entries)
    print_dual(color_text(
        f"[OK] Citations for the methods used in this run written to "
        f"'{os.path.join(args.output_dir, BIB_FILE)}' ({len(bib_entries)} entries).", 'green'), f_out)

    print_section("[6] SUMMARY & FILES", f_out)
    print_dual("Status         : OK", f_out)
    for path in written_data:
        print_dual(f"Data           : {path}", f_out)
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
        for key in COMPONENT_NAMES:
            k_pos, energy, value = zip(*rows[key])
            plot_scalar_on_bands(high_sym, k_pos, energy, value, key, _is_gamma, is_signed=True)


if __name__ == "__main__":
    main()
