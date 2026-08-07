#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.4.1"  # BSSE text updated: now Stage 3 (stb-adsorbBsse became the new Stage 2, generating
                    # bsse/ folders at each site's RELAXED geometry instead of stb-adsorb's old
                    # pre-relaxation guess); no functional change to this tool's own aggregation.

import os
import re
import sys
import glob
import shutil
import argparse
import numpy as np
import ase
import matplotlib.pyplot as plt
from datetime import datetime
from pymatgen.io.ase import AseAtomsAdaptor
from stb.core import siesta_log, structure_io
from stb.core.siesta_log import check_scf_and_force, report_quality_diagnostics
from stb.core.cli import color_text, show_intro, print_dual, print_section, print_table, capture_library_noise
from stb.core.ase_view import view_structure_interactive

REPORT_FILE = "adsorption_report.txt"
SITES_REPORT_FILE = "adsorption_sites.txt"
RANKING_PLOT_FILE = "adsorption_ranking.png"
_HEIGHT_SUFFIX_RE = re.compile(r'_h[\d.]+$')


def read_site_table(sites_root):
    """Parses the '# SITE_TABLE' section stb-adsorb writes in
    sites/adsorption_sites.txt into {label: (adsorbate_name, height)} --
    the source of truth for which adsorbate/<name>/ reference and which
    height each site belongs to (more robust than re-guessing either from
    the folder name by regex). Returns None if the report (or its table
    section) isn't found, so the caller can fall back to treating every
    site as belonging to a single, unnamed adsorbate -- e.g. a report the
    user deleted, or folders assembled by hand rather than via stb-adsorb.
    """
    report_path = os.path.join(sites_root, SITES_REPORT_FILE)
    if not os.path.isfile(report_path):
        return None
    table = {}
    in_table = False
    with open(report_path) as f:
        for line in f:
            if line.startswith("# SITE_TABLE"):
                in_table = True
                continue
            if not in_table:
                continue
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            label, ads_name = parts[0], parts[1]
            try:
                height = float(parts[2])
            except ValueError:
                continue
            table[label] = (ads_name, height)
    return table or None


def read_bsse_energy(bsse_dir, file_name):
    """Reads the two ghost-fragment reference energies stb-adsorb's
    --bsse-correction (default ON) writes under 'bsse/<site_label>/' -- a
    tree parallel to 'sites/', not nested inside a site's own folder --
    '<bsse_dir>/bsse_slab/<file_name>' and
    '<bsse_dir>/bsse_adsorbate/<file_name>' -- and returns (e_bsse_slab,
    e_bsse_adsorbate), or (None, None) if either is missing/unreadable.
    Best-effort per site (never blocks the uncorrected result for that
    site or any other). Callers are expected to have already confirmed
    '<bsse_dir>/bsse_slab' exists before calling this (main() does, to
    decide whether to even attempt reading BSSE energies for a site) --
    this function itself only worries about whether the energies inside
    are actually readable.
    """
    slab_path = os.path.join(bsse_dir, "bsse_slab", file_name)
    ads_path = os.path.join(bsse_dir, "bsse_adsorbate", file_name)
    e_slab = siesta_log.get_free_energy(slab_path)
    e_ads = siesta_log.get_free_energy(ads_path)
    if e_slab is None or e_ads is None:
        return None, None
    return e_slab, e_ads


def write_curve_plot(dat_path, rows):
    """Writes <dat_path> plus a companion .gplot (site index vs. adsorption
    energy, uncorrected and BSSE-corrected if available) -- same
    <name>.dat + <name>.gplot convention as the rest of the suite (e.g.
    convergence_analysis.py::write_curve_plot, phonons_pos.py's
    write_thermal_plots). Sites are discrete/categorical (not a continuous
    sweep), so this is a scatter, not a line.
    """
    with open(dat_path, 'w') as f:
        f.write("# Adsorption energy per candidate site\n")
        f.write("# 1:SiteIndex 2:E_ads(eV) 3:E_ads_BSSE(eV) 4:Adsorbate 5:Label\n")
        for i, row in enumerate(rows, start=1):
            bsse_str = f"{row.e_ads_bsse:.6f}" if row.e_ads_bsse is not None else "nan"
            f.write(f"{i}  {row.e_ads:.6f}  {bsse_str}  {row.ads_name}  {row.label}\n")

    base = os.path.splitext(dat_path)[0]
    base_name = os.path.basename(base)
    gplot_path = f"{base}.gplot"
    dat_name = os.path.basename(dat_path)
    has_bsse = any(row.e_ads_bsse is not None for row in rows)
    plot_lines = (f'"{dat_name}" using 1:2 with points pt 7 ps 2 lc rgb "#2255cc" title "E_{{ads}}"')
    if has_bsse:
        plot_lines += (f', "{dat_name}" using 1:3 with points pt 9 ps 2 lc rgb "#cc5522" '
                        'title "E_{ads} (BSSE)"')
    with open(gplot_path, 'w') as f:
        f.writelines([
            '# --- STB Plot Configuration ---\n',
            '# Generated by stb-adsorbAnalysis\n',
            'set terminal pdfcairo enhanced color font "Arial,14" size 7,5\n',
            f'set output "{base_name}.pdf"\n\n',
            'set title "Adsorption energy per site"\n',
            'set xlabel "Site index"\n',
            'set ylabel "E_{ads} (eV)"\n',
            'set grid\n',
            'set key top right\n',
            f'plot {plot_lines}\n',
        ])
    return gplot_path


def write_height_curve_plot(family, family_rows, out_dir):
    """Writes an E_ads-vs-height approach curve for one site "family"
    (same site + adsorbate, swept across --height-sweep) -- <family>.dat +
    <family>.gplot, same convention as write_curve_plot above. Only called
    for families with 2+ heights (see main()).
    """
    family_rows = sorted(family_rows, key=lambda r: r.height)
    dat_path = os.path.join(out_dir, f"height_curve_{family}.dat")
    with open(dat_path, 'w') as f:
        f.write(f"# Adsorption energy vs. height for site family '{family}'\n")
        f.write("# 1:Height(Ang) 2:E_ads(eV) 3:E_ads_BSSE(eV)\n")
        for row in family_rows:
            bsse_str = f"{row.e_ads_bsse:.6f}" if row.e_ads_bsse is not None else "nan"
            f.write(f"{row.height:.4f}  {row.e_ads:.6f}  {bsse_str}\n")

    base_name = f"height_curve_{family}"
    gplot_path = os.path.join(out_dir, f"{base_name}.gplot")
    dat_name = os.path.basename(dat_path)
    has_bsse = any(row.e_ads_bsse is not None for row in family_rows)
    plot_lines = (f'"{dat_name}" using 1:2 with linespoints lw 2 pt 7 lc rgb "#2255cc" title "E_{{ads}}"')
    if has_bsse:
        plot_lines += (f', "{dat_name}" using 1:3 with linespoints lw 2 pt 9 lc rgb "#cc5522" '
                        'title "E_{ads} (BSSE)"')
    with open(gplot_path, 'w') as f:
        f.writelines([
            '# --- STB Plot Configuration ---\n',
            '# Generated by stb-adsorbAnalysis\n',
            'set terminal pdfcairo enhanced color font "Arial,14" size 7,5\n',
            f'set output "{base_name}.pdf"\n\n',
            f'set title "Approach curve: {family}"\n',
            'set xlabel "Height (Ang)"\n',
            'set ylabel "E_{ads} (eV)"\n',
            'set grid\n',
            'set key top right\n',
            f'plot {plot_lines}\n',
        ])
    return dat_path, gplot_path


def plot_adsorption_ranking(rows, out_path, show=False):
    """Bar chart of E_ads (and E_ads_BSSE, where available) per site,
    ranked most stable first -- the matplotlib complement to
    write_curve_plot's gnuplot .dat/.gplot pair (kept as-is; this is
    additive, not a replacement), same single-series-bar-chart style as
    adsorb.py::plot_ml_ranking. `rows` is the same list `main()` already
    sorted by `e_ads` (ascending -- most stable first).

    `show=True` (--view-plots) additionally blocks on plt.show() before
    the figure is closed -- same convention as adsorb.py's --view-plots.
    """
    labels = [f"{r.ads_name or 'default'}\n{r.label}" for r in rows]
    e_ads = [r.e_ads for r in rows]
    has_bsse = [r.e_ads_bsse is not None for r in rows]
    x = np.arange(len(rows))
    width = 0.35 if any(has_bsse) else 0.6
    fig, ax = plt.subplots(figsize=(max(6, 0.9 * len(labels)), 5))
    colors = ['tab:green' if i == 0 else 'tab:blue' for i in range(len(rows))]
    ax.bar(x - width / 2 if any(has_bsse) else x, e_ads, width=width, color=colors, label="E_ads")
    if any(has_bsse):
        bsse_vals = [r.e_ads_bsse if r.e_ads_bsse is not None else np.nan for r in rows]
        ax.bar(x + width / 2, bsse_vals, width=width, color='tab:orange', label="E_ads (BSSE)")
        ax.legend(fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("E_ads (eV)")
    ax.set_title("Adsorption energy per site (most stable first)")
    ax.axhline(0.0, color='gray', linestyle='--', linewidth=0.8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    if show:
        plt.show()
    plt.close(fig)


def read_site_geometry_atoms(folder):
    """Returns an ase.Atoms for whatever geometry is available in `folder`:
    prefers a finished SIESTA run's own '<label>.XV' (the actual relaxed
    geometry SIESTA used, converted the same way mlff_analysis.py already
    does for a sisl Geometry), and falls back to 'structure.fdf' (the
    input Stage 1 always writes, via this suite's own core/structure_io.py
    reader) when no .XV exists yet -- e.g. before SIESTA has been run in
    that folder at all. Returns None if neither is readable, so callers
    can skip that frame instead of crashing --view over one missing/
    unrun folder.
    """
    xv_files = sorted(glob.glob(os.path.join(folder, "*.XV")))
    if xv_files:
        from stb.core.deps import require_sisl
        sisl = require_sisl()
        geom = sisl.get_sile(xv_files[0]).read_geometry()
        return ase.Atoms(symbols=[a.symbol for a in geom.atoms],
                          positions=geom.xyz, cell=np.array(geom.cell), pbc=True)
    fdf_path = os.path.join(folder, "structure.fdf")
    if os.path.isfile(fdf_path):
        pmg_structure = structure_io.to_pymatgen(structure_io.read_fdf(fdf_path))
        return AseAtomsAdaptor.get_atoms(pmg_structure)
    return None


class SiteRow:
    """One analyzed site: label, which adsorbate/height it belongs to, and
    its computed energies. A plain attribute-holder (not a dataclass, to
    keep this a light dependency-free module) instead of a bare tuple --
    the field count grew enough (multi-adsorbate + height-sweep) that
    positional tuple unpacking everywhere had become error-prone to edit.
    """
    def __init__(self, label, ads_name, height, e_site, e_ads, e_ads_bsse, scf_ok, max_force):
        self.label = label
        self.ads_name = ads_name
        self.height = height
        self.e_site = e_site
        self.e_ads = e_ads
        self.e_ads_bsse = e_ads_bsse
        self.scf_ok = scf_ok
        self.max_force = max_force


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Computes adsorption energies from an stb-adsorb sweep: "
        "E_ads = E_site - E_clean_slab - E_adsorbate.", 'bold')}""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s\n"
               "  %(prog)s --dir adsorption_run --file calc.out\n"
               "  %(prog)s --apply best_structure.fdf\n"
    )

    parser.add_argument("--dir", type=str, default="adsorption_run",
                         help="Root directory containing 'clean_slab/', 'adsorbate*/' and 'sites/' "
                              "(default: adsorption_run -- stb-adsorb's own default --output-dir, "
                              "so this stage runs against that one's output with no extra flag in "
                              "the common case).")
    parser.add_argument("--file", type=str, default="calc.out",
                         help="SIESTA output filename inside each folder (default: calc.out).")
    parser.add_argument("-o", "--output", type=str, default="adsorption_curve.dat",
                         help="Output data file name (default: adsorption_curve.dat).")
    parser.add_argument("--apply", type=str, default=None, metavar="STRUCTURE_FDF",
                         help="Copy the most stable site's structure.fdf (BSSE-corrected ranking "
                             "when available for every site, else the uncorrected ranking) to "
                             "this path.")
    parser.add_argument("--force-tolerance", type=float, default=0.05,
                         help="Residual atomic force in eV/Ang (default: 0.05, same as "
                              "stb-cohesiveAnalysis) above which a site's calc.out is flagged as "
                              "possibly not relaxed -- E_ads would then reflect a strained/off-"
                              "equilibrium geometry rather than the true adsorption minimum. "
                              "Advisory only, never blocks the result.")
    parser.add_argument("--save-report", action="store_true",
                         help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--view", action="store_true",
                         help="Open clean_slab, every adsorbate reference and every site "
                              "(preferring each folder's own <label>.XV -- the actual relaxed "
                              "geometry, if SIESTA has been run there -- falling back to its "
                              "input structure.fdf otherwise) in ASE's interactive 3D viewer as "
                              "a multi-frame browser, after everything else has finished. "
                              "Needs a display.")
    parser.add_argument("--view-plots", action="store_true",
                         help="Also show the matplotlib adsorption-ranking plot on screen "
                              "(one blocking window), instead of only saving it as a PNG. Off "
                              "by default.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-adsorbAnalysis {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("Analyze an adsorption-energy study:", 'bold'))
    print("-" * 60)

    clean_slab_out = os.path.join(args.dir, "clean_slab", args.file)
    sites_root = os.path.join(args.dir, "sites")
    bsse_root = os.path.join(args.dir, "bsse")

    if not os.path.isdir(sites_root):
        print(color_text(f"[ERROR] '{sites_root}' not found. Did you run stb-adsorb?", 'red'))
        sys.exit(1)

    site_table = read_site_table(sites_root)

    library_warnings = []  # collected via capture_library_noise, reported in [5]
    view_frames = []  # (label, ase.Atoms), only populated when --view is given

    report_path = REPORT_FILE if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    print_dual(f"{color_text('===== ADSORPTION ENERGY REPORT =====', 'magenta')}", f_out)

    print_section('[0] RUN METADATA', f_out)
    print_dual(f"Date/time  : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", f_out)
    print_dual(f"Directory  : {args.dir}", f_out)
    print_dual(f"Output file: {args.file}", f_out)
    if site_table is None:
        print_dual(color_text(
            "[NOTE] No 'sites/adsorption_sites.txt' site table found -- treating every site "
            "as a single unnamed adsorbate (pre-multi-adsorbate/height-sweep layout, or the "
            "report was deleted).", 'yellow'), f_out)

    print_section('[1] REFERENCE ENERGIES', f_out)
    e_clean_slab = siesta_log.get_free_energy(clean_slab_out)
    if e_clean_slab is None:
        print_dual(color_text(f"[ERROR] Could not read energy from '{clean_slab_out}'.", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)
    print_dual(f"E_clean_slab : {e_clean_slab:.6f} eV  ({clean_slab_out})", f_out)
    report_quality_diagnostics("clean_slab", clean_slab_out, args.force_tolerance, f_out)
    if args.view:
        with capture_library_noise(library_warnings, "sisl/ASE (clean_slab geometry)"):
            atoms = read_site_geometry_atoms(os.path.join(args.dir, "clean_slab"))
        if atoms is not None:
            view_frames.append(("clean_slab", atoms))

    # Adsorbate names present, from the site table if available, else a
    # single unnamed one mapping to the legacy 'adsorbate/' folder.
    ads_names = sorted({name for name, _h in site_table.values()}) if site_table else [None]
    adsorbate_energies = {}  # name -> energy
    for name in ads_names:
        candidate_dirs = []
        if name:
            candidate_dirs.append(os.path.join(args.dir, f"adsorbate_{name}"))
        candidate_dirs.append(os.path.join(args.dir, "adsorbate"))
        ads_dir = next((d for d in candidate_dirs if os.path.isdir(d)), None)
        if ads_dir is None:
            print_dual(color_text(
                f"[ERROR] No adsorbate reference folder found for "
                f"'{name or '(default)'}' (tried: {', '.join(candidate_dirs)}).", 'red'), f_out)
            if f_out:
                f_out.close()
            sys.exit(1)
        ads_out_path = os.path.join(ads_dir, args.file)
        e_ads = siesta_log.get_free_energy(ads_out_path)
        if e_ads is None:
            print_dual(color_text(f"[ERROR] Could not read energy from "
                                   f"'{ads_out_path}'.", 'red'), f_out)
            if f_out:
                f_out.close()
            sys.exit(1)
        adsorbate_energies[name] = e_ads
        print_dual(f"E_adsorbate ({name or 'default'}) : {e_ads:.6f} eV  ({ads_dir})", f_out)
        report_quality_diagnostics(f"adsorbate ({name or 'default'})", ads_out_path,
                                    args.force_tolerance, f_out)
        if args.view:
            with capture_library_noise(library_warnings, f"sisl/ASE (adsorbate {name or 'default'} geometry)"):
                atoms = read_site_geometry_atoms(ads_dir)
            if atoms is not None:
                view_frames.append((f"adsorbate({name or 'default'})", atoms))

    print_section('[2] SITE RESULTS: CONFIGURATION COUNT & TABLE', f_out)
    site_dirs = sorted(
        d for d in os.listdir(sites_root)
        if os.path.isdir(os.path.join(sites_root, d)) and d.startswith("site_")
    )
    if not site_dirs:
        print_dual(color_text(f"[ERROR] No 'site_*' folders found in '{sites_root}'.", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    rows = []
    n_skip_missing = 0
    n_skip_unparseable = 0
    n_bsse_found = 0
    incomplete_bsse_labels = []
    scf_warn_labels = []
    force_warn_labels = []
    bsse_scf_warn_labels = []
    table_rows = []
    for label in site_dirs:
        site_dir = os.path.join(sites_root, label)
        ads_name, height = site_table.get(label, (None, None)) if site_table else (None, None)

        if args.view:
            with capture_library_noise(library_warnings, f"sisl/ASE (site {label} geometry)"):
                atoms = read_site_geometry_atoms(site_dir)
            if atoms is not None:
                view_frames.append((label, atoms))

        out_path = os.path.join(site_dir, args.file)
        if not os.path.exists(out_path):
            n_skip_missing += 1
            table_rows.append(([label, ads_name or '--', "--" if height is None else f"{height:.2f}",
                                 "SKIP (missing calc.out)", "--", "--", "--"], 'yellow'))
            continue
        e_site = siesta_log.get_free_energy(out_path)
        if e_site is None:
            n_skip_unparseable += 1
            table_rows.append(([label, ads_name or '--', "--" if height is None else f"{height:.2f}",
                                 "SKIP (unparseable energy)", "--", "--", "--"], 'yellow'))
            continue
        e_ads = e_site - e_clean_slab - adsorbate_energies[ads_name]
        scf_ok, max_force = check_scf_and_force(out_path)
        if not scf_ok:
            scf_warn_labels.append(label)
        if max_force is not None and max_force > args.force_tolerance:
            force_warn_labels.append(label)

        e_ads_bsse = None
        bsse_dir = os.path.join(bsse_root, label)
        if os.path.isdir(os.path.join(bsse_dir, "bsse_slab")):
            e_bsse_slab, e_bsse_ads = read_bsse_energy(bsse_dir, args.file)
            if e_bsse_slab is not None and e_bsse_ads is not None:
                e_ads_bsse = e_site - e_bsse_slab - e_bsse_ads
                n_bsse_found += 1
                bsse_slab_scf, _f1 = check_scf_and_force(
                    os.path.join(bsse_dir, "bsse_slab", args.file))
                bsse_ads_scf, _f2 = check_scf_and_force(
                    os.path.join(bsse_dir, "bsse_adsorbate", args.file))
                if not (bsse_slab_scf and bsse_ads_scf):
                    bsse_scf_warn_labels.append(label)
            else:
                incomplete_bsse_labels.append(label)

        rows.append(SiteRow(label, ads_name, height, e_site, e_ads, e_ads_bsse, scf_ok, max_force))
        row_color = 'yellow' if (not scf_ok or label in force_warn_labels) else None
        table_rows.append(([
            label, ads_name or '--', "--" if height is None else f"{height:.2f}",
            f"{e_ads:.6f}", f"{e_ads_bsse:.6f}" if e_ads_bsse is not None else "--",
            "OK" if scf_ok else "WARN", f"{max_force:.4f}" if max_force is not None else "--",
        ], row_color))

    n_skipped = n_skip_missing + n_skip_unparseable
    n_bsse_absent = len(site_dirs) - n_bsse_found - len(incomplete_bsse_labels)
    print_dual(f"Site folders found  : {len(site_dirs)}", f_out)
    print_dual(f"Read successfully   : {len(rows)}", f_out)
    print_dual(f"Skipped             : {n_skipped}  (missing {args.file}: {n_skip_missing}, "
                f"unparseable energy: {n_skip_unparseable})", f_out)
    if len(ads_names) > 1:
        print_dual("Per-adsorbate breakdown:", f_out)
        for name in ads_names:
            n_found_name = sum(1 for _l in site_dirs
                                if (site_table.get(_l, (None, None))[0] if site_table else None) == name)
            n_read_name = sum(1 for r in rows if r.ads_name == name)
            print_dual(f"  {name or '(default)'}: {n_found_name} found, {n_read_name} read", f_out)
    families = {}
    for label in site_dirs:
        families.setdefault(_HEIGHT_SUFFIX_RE.sub('', label), []).append(label)
    swept_families = {fam: labels for fam, labels in families.items() if len(labels) > 1}
    if swept_families:
        print_dual(f"Height-sweep families: {len(swept_families)} site(s) swept across "
                    f"multiple heights (family -> heights found):", f_out)
        for fam, labels in swept_families.items():
            print_dual(f"  {fam}: {len(labels)} height(s)", f_out)
    print_dual(f"BSSE coverage       : complete {n_bsse_found}, incomplete "
                f"{len(incomplete_bsse_labels)}, absent {n_bsse_absent}  (of {len(site_dirs)} site(s))", f_out)

    print_dual(color_text(
        "\n[BSSE PHYSICS CHECK] E_ads(BSSE) = E_site - E_bsse_slab - E_bsse_adsorbate, where "
        "E_bsse_slab/E_bsse_adsorbate are single-point ghost-fragment energies read from "
        "'bsse/<site>/bsse_slab/' and 'bsse/<site>/bsse_adsorbate/' -- written by stb-adsorbBsse "
        "(Stage 2) at the site's actual RELAXED geometry (from its finished siesta.XV), not the "
        "pre-relaxation guess, since each fragment must be evaluated at the same geometry as the "
        "real, relaxed site for the correction to mean anything (the standard Boys-Bernardi "
        "counterpoise, applied to both the substrate and the adsorbate). Full explanation: "
        "examples/4.8-adsorption/README.md.", 'cyan'), f_out)

    print_dual("", f_out)
    headers = ["Site", "Adsorbate", "Height", "E_ads(eV)", "E_ads_BSSE(eV)", "SCF", "MaxF(eV/A)"]
    print_table(headers, table_rows, f_out)

    if incomplete_bsse_labels:
        print_dual(color_text(
            f"[WARNING] {len(incomplete_bsse_labels)} site(s) have a 'bsse_slab/'/"
            "'bsse_adsorbate/' folder but incomplete/unreadable results -- BSSE-corrected "
            f"energy skipped for: {', '.join(incomplete_bsse_labels)}.", 'yellow'), f_out)
    if scf_warn_labels:
        print_dual(color_text(
            f"[WARNING] {len(scf_warn_labels)} site(s) never confirmed SCF convergence -- "
            f"their E_ads may be unreliable: {', '.join(scf_warn_labels)}.", 'yellow'), f_out)
    if force_warn_labels:
        print_dual(color_text(
            f"[WARNING] {len(force_warn_labels)} site(s) have residual force above "
            f"--force-tolerance ({args.force_tolerance} eV/Ang), possibly not relaxed: "
            f"{', '.join(force_warn_labels)}.", 'yellow'), f_out)
    if bsse_scf_warn_labels:
        print_dual(color_text(
            f"[WARNING] {len(bsse_scf_warn_labels)} site(s)' BSSE ghost-fragment calculation(s) "
            f"never confirmed SCF convergence -- their E_ads_BSSE may be unreliable: "
            f"{', '.join(bsse_scf_warn_labels)}.", 'yellow'), f_out)

    if not rows:
        print_dual(color_text("\n[ERROR] No valid site results found.", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    # Primary ranking always uses the uncorrected E_ads -- it's available
    # for every row, so it's the only metric safe to compare across all
    # sites uniformly (it's already referenced against each row's own
    # adsorbate energy, so comparing across different adsorbates is
    # meaningful too, unlike stb-adsorb --ml-rank's raw MACE energies).
    # BSSE-corrected values are reported alongside (and separately
    # ranked below, when complete for every site) rather than silently
    # mixed into the same sort.
    rows.sort(key=lambda r: r.e_ads)

    print_section('[3] SUMMARY & PLOT', f_out)
    print_dual(f"Sites analyzed : {len(rows)} (skipped: {n_skipped})", f_out)
    best = rows[0]
    verdict = "exothermic (favorable)" if best.e_ads < 0 else "endothermic (unfavorable)"
    print_dual(f"{color_text('Most stable site (uncorrected):', 'green')} {best.label}  "
                f"(E_ads = {best.e_ads:.6f} eV, {verdict})", f_out)

    apply_source_label = best.label
    if n_bsse_found == 0:
        print_dual(color_text(
            "\n[NOTE] No BSSE-corrected results found -- run stb-adsorbBsse (Stage 2) once these "
            "sites have finished relaxing, then run SIESTA in the 'bsse/site_*/' folders it "
            "writes, for a corrected reference.", 'yellow'), f_out)
    else:
        bsse_rows = [r for r in rows if r.e_ads_bsse is not None]
        if len(bsse_rows) == len(rows):
            bsse_ranked = sorted(bsse_rows, key=lambda r: r.e_ads_bsse)
            bsse_best = bsse_ranked[0]
            bsse_verdict = "exothermic (favorable)" if bsse_best.e_ads_bsse < 0 else "endothermic (unfavorable)"
            print_dual(f"{color_text('Most stable site (BSSE-corrected):', 'green')} "
                        f"{bsse_best.label}  (E_ads = {bsse_best.e_ads_bsse:.6f} eV, {bsse_verdict})", f_out)
            shift = bsse_best.e_ads_bsse - bsse_best.e_ads
            print_dual(f"BSSE correction at that site: {shift:+.6f} eV (uncorrected LCAO "
                        "adsorption energies systematically over-bind -- expect this to make "
                        "the energy less negative)", f_out)
            apply_source_label = bsse_best.label
        else:
            print_dual(color_text(
                f"\n[NOTE] BSSE-corrected energy available for only {len(bsse_rows)}/{len(rows)} "
                "site(s) -- not reporting a BSSE-ranked \"most stable\" until all sites have it.",
                'yellow'), f_out)

    if len(ads_names) > 1:
        print_dual(f"\n{color_text('Best site per adsorbate (uncorrected):', 'cyan')}", f_out)
        for name in ads_names:
            per_ads = [r for r in rows if r.ads_name == name]
            if not per_ads:
                continue
            best_per_ads = min(per_ads, key=lambda r: r.e_ads)
            print_dual(f"  {name}: {best_per_ads.label} (E_ads = {best_per_ads.e_ads:.6f} eV)", f_out)
            per_ads_bsse = [r for r in per_ads if r.e_ads_bsse is not None]
            if per_ads_bsse and len(per_ads_bsse) == len(per_ads):
                best_per_ads_bsse = min(per_ads_bsse, key=lambda r: r.e_ads_bsse)
                print_dual(f"     BSSE-corrected: {best_per_ads_bsse.label} "
                            f"(E_ads = {best_per_ads_bsse.e_ads_bsse:.6f} eV)", f_out)

    gplot_path = write_curve_plot(args.output, rows)
    print_dual(f"\n{color_text('[Saved]', 'cyan')} Curve data -> {args.output}, {gplot_path} "
                f"(cd {os.path.dirname(args.output) or '.'} && gnuplot {os.path.basename(gplot_path)})",
                f_out)

    # Height-sweep approach curves: one per "family" (same site + same
    # adsorbate, swept across height) that has 2+ heights among the rows
    # actually analyzed above.
    families = {}
    for row in rows:
        family = _HEIGHT_SUFFIX_RE.sub('', row.label)
        families.setdefault(family, []).append(row)
    out_dir = os.path.dirname(args.output) or "."
    for family, family_rows in families.items():
        if len({r.height for r in family_rows}) < 2:
            continue
        dat_path, height_gplot = write_height_curve_plot(family, family_rows, out_dir)
        print_dual(f"{color_text('[Saved]', 'cyan')} Height curve -> {dat_path}, {height_gplot}", f_out)

    ranking_plot_path = os.path.join(out_dir, RANKING_PLOT_FILE)
    with capture_library_noise(library_warnings, "matplotlib plot_adsorption_ranking"):
        plot_adsorption_ranking(rows, ranking_plot_path, show=args.view_plots)
    print_dual(f"{color_text('[Saved]', 'cyan')} Ranking plot -> {ranking_plot_path}", f_out)

    if report_path:
        print_dual(f"{color_text('[Saved]', 'cyan')} Report     -> {report_path}", f_out)

    # --- [4] APPLY: always present (fixed section numbering), even when
    # --apply wasn't passed -- same "Not requested" pattern as adsorb.py's
    # own [3] ML PRE-SCREENING.
    print_section('[4] APPLY', f_out)
    if not args.apply:
        print_dual("  Not requested (pass --apply <path> to copy the most stable site's "
                    "structure.fdf -- BSSE-corrected ranking when available for every site, "
                    "else the uncorrected ranking -- to that path).", f_out)
    else:
        src = os.path.join(sites_root, apply_source_label, "structure.fdf")
        try:
            shutil.copy(src, args.apply)
        except OSError as e:
            print_dual(color_text(f"[ERROR] Could not copy '{src}' to '{args.apply}': {e}", 'red'), f_out)
        else:
            print_dual(f"{color_text('[Applied]', 'green')} {apply_source_label} -> {args.apply}", f_out)

    # --- [5] LIBRARY WARNINGS: always last. ---
    print_section('[5] LIBRARY WARNINGS', f_out)
    if library_warnings:
        print_dual(color_text(
            "Messages emitted by external libraries (sisl, pymatgen, matplotlib) during this "
            "run -- collected here instead of interleaved with the report above; harmless in "
            "almost every case, but worth a look if a section above looks suspicious.", 'cyan'), f_out)
        for entry in library_warnings:
            print_dual(entry, f_out)
    else:
        print_dual("No library warnings.", f_out)

    if f_out:
        f_out.close()

    if args.view:
        print(f"\n{color_text('--view:', 'cyan')} opening {len(view_frames)} frame(s) in ASE's "
              "interactive viewer (use the frame slider/menu to step through them):")
        for i, (label, _atoms) in enumerate(view_frames):
            print(f"  {i} = {label}")
        if view_frames:
            view_structure_interactive([atoms for _label, atoms in view_frames])
        else:
            print(color_text("  No readable structures found -- nothing to view.", 'yellow'))


if __name__ == "__main__":
    main()
