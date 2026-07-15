#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.3.0"  # SCF convergence + residual-force diagnostics, per-adsorbate BSSE best

import os
import re
import sys
import shutil
import argparse
from datetime import datetime
from stb.core import siesta_log
from stb.core.cli import color_text, show_intro, print_dual

REPORT_FILE = "adsorption_report.txt"
SITES_REPORT_FILE = "adsorption_sites.txt"
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


def read_bsse_energy(site_dir, file_name):
    """Reads the two ghost-fragment reference energies stb-adsorb's
    --bsse-correction (default ON) writes under a site folder --
    'bsse_slab/<file_name>' and 'bsse_adsorbate/<file_name>' -- and returns
    (e_bsse_slab, e_bsse_adsorbate), or (None, None) if either is missing/
    unreadable. Best-effort per site (never blocks the uncorrected result
    for that site or any other).
    """
    slab_path = os.path.join(site_dir, "bsse_slab", file_name)
    ads_path = os.path.join(site_dir, "bsse_adsorbate", file_name)
    if not os.path.isdir(os.path.join(site_dir, "bsse_slab")):
        return None, None
    e_slab = siesta_log.get_free_energy(slab_path)
    e_ads = siesta_log.get_free_energy(ads_path)
    if e_slab is None or e_ads is None:
        return None, None
    return e_slab, e_ads


def check_scf_and_force(out_path):
    """Returns (scf_converged, max_force) for one SIESTA .out file --
    core/siesta_log.py's get_scf_convergence/get_max_force, the same
    numerical-quality diagnostics stb-cohesiveAnalysis already runs on its
    own full-structure reference. Neither call is expensive (single
    sequential file read each), so it's fine to run per folder even for a
    large --all-sites/--height-sweep/--adsorbate sweep.
    """
    scf_ok, _iterations = siesta_log.get_scf_convergence(out_path)
    max_force = siesta_log.get_max_force(out_path)
    return scf_ok, max_force


def report_quality_diagnostics(label, out_path, force_tolerance, f_out):
    """Prints (and persists) a numerical-quality warning for one reference
    folder (clean_slab or an isolated-adsorbate) if its SCF cycle never
    confirmed convergence, or its residual force exceeds force_tolerance --
    silent when both are fine, matching this suite's "advisory only, don't
    clutter a clean run" convention (e.g. convergence_analysis.py's SCF
    gating, cohesive_analysis.py's --force-tolerance check). Site-level
    diagnostics are folded into the [2] SITE RESULTS table instead, since
    there can be many of them.
    """
    scf_ok, max_force = check_scf_and_force(out_path)
    if not scf_ok:
        print_dual(color_text(
            f"  [WARNING] Could not confirm SCF convergence for {label} ('{out_path}') -- "
            "this energy may be unreliable.", 'yellow'), f_out)
    if max_force is not None and max_force > force_tolerance:
        print_dual(color_text(
            f"  [WARNING] Residual force on {label} ({max_force:.4f} eV/Ang) exceeds "
            f"--force-tolerance ({force_tolerance} eV/Ang) -- this geometry may not be "
            "relaxed.", 'yellow'), f_out)


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
               "  %(prog)s --dir . --file calc.out\n"
               "  %(prog)s --dir . --apply best_structure.fdf\n"
    )

    parser.add_argument("--dir", type=str, default=".",
                         help="Root directory containing 'clean_slab/', 'adsorbate*/' and 'sites/' "
                              "(default: current directory).")
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

    if not os.path.isdir(sites_root):
        print(color_text(f"[ERROR] '{sites_root}' not found. Did you run stb-adsorb?", 'red'))
        sys.exit(1)

    site_table = read_site_table(sites_root)

    with open(REPORT_FILE, 'w') as f_out:
        print_dual(f"{color_text('===== ADSORPTION ENERGY REPORT =====', 'magenta')}", f_out)

        print_dual(f"\n{color_text('[0] RUN METADATA', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        print_dual(f"Date/time  : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", f_out)
        print_dual(f"Directory  : {args.dir}", f_out)
        print_dual(f"Output file: {args.file}", f_out)
        if site_table is None:
            print_dual(color_text(
                "[NOTE] No 'sites/adsorption_sites.txt' site table found -- treating every site "
                "as a single unnamed adsorbate (pre-multi-adsorbate/height-sweep layout, or the "
                "report was deleted).", 'yellow'), f_out)

        print_dual(f"\n{color_text('[1] REFERENCE ENERGIES', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        e_clean_slab = siesta_log.get_free_energy(clean_slab_out)
        if e_clean_slab is None:
            print_dual(color_text(f"[ERROR] Could not read energy from '{clean_slab_out}'.", 'red'), f_out)
            sys.exit(1)
        print_dual(f"E_clean_slab : {e_clean_slab:.6f} eV  ({clean_slab_out})", f_out)
        report_quality_diagnostics("clean_slab", clean_slab_out, args.force_tolerance, f_out)

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
                sys.exit(1)
            ads_out_path = os.path.join(ads_dir, args.file)
            e_ads = siesta_log.get_free_energy(ads_out_path)
            if e_ads is None:
                print_dual(color_text(f"[ERROR] Could not read energy from "
                                       f"'{ads_out_path}'.", 'red'), f_out)
                sys.exit(1)
            adsorbate_energies[name] = e_ads
            print_dual(f"E_adsorbate ({name or 'default'}) : {e_ads:.6f} eV  ({ads_dir})", f_out)
            report_quality_diagnostics(f"adsorbate ({name or 'default'})", ads_out_path,
                                        args.force_tolerance, f_out)

        print_dual(f"\n{color_text('[2] SITE RESULTS', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        site_dirs = sorted(
            d for d in os.listdir(sites_root)
            if os.path.isdir(os.path.join(sites_root, d)) and d.startswith("site_")
        )
        if not site_dirs:
            print_dual(color_text(f"[ERROR] No 'site_*' folders found in '{sites_root}'.", 'red'), f_out)
            sys.exit(1)

        rows = []
        n_skipped = 0
        n_bsse_found = 0
        incomplete_bsse_labels = []
        scf_warn_labels = []
        force_warn_labels = []
        bsse_scf_warn_labels = []
        header = (f"{'Site':<30} {'Adsorbate':<12} {'Height':<9} {'E_ads(eV)':<14} "
                  f"{'E_ads_BSSE(eV)':<16}{'SCF':<6}{'MaxF(eV/A)':<12}")
        print_dual(header, f_out)
        print_dual("-" * len(header), f_out)
        for label in site_dirs:
            site_dir = os.path.join(sites_root, label)
            ads_name, height = site_table.get(label, (None, None)) if site_table else (None, None)

            out_path = os.path.join(site_dir, args.file)
            if not os.path.exists(out_path):
                n_skipped += 1
                print_dual(f"{label:<30} {color_text('SKIP', 'yellow')} (missing {args.file})", f_out)
                continue
            e_site = siesta_log.get_free_energy(out_path)
            if e_site is None:
                n_skipped += 1
                print_dual(f"{label:<30} {color_text('SKIP', 'yellow')} (could not parse energy)", f_out)
                continue
            e_ads = e_site - e_clean_slab - adsorbate_energies[ads_name]
            scf_ok, max_force = check_scf_and_force(out_path)
            if not scf_ok:
                scf_warn_labels.append(label)
            if max_force is not None and max_force > args.force_tolerance:
                force_warn_labels.append(label)

            e_ads_bsse = None
            if os.path.isdir(os.path.join(site_dir, "bsse_slab")):
                e_bsse_slab, e_bsse_ads = read_bsse_energy(site_dir, args.file)
                if e_bsse_slab is not None and e_bsse_ads is not None:
                    e_ads_bsse = e_site - e_bsse_slab - e_bsse_ads
                    n_bsse_found += 1
                    bsse_slab_scf, _f1 = check_scf_and_force(
                        os.path.join(site_dir, "bsse_slab", args.file))
                    bsse_ads_scf, _f2 = check_scf_and_force(
                        os.path.join(site_dir, "bsse_adsorbate", args.file))
                    if not (bsse_slab_scf and bsse_ads_scf):
                        bsse_scf_warn_labels.append(label)
                else:
                    incomplete_bsse_labels.append(label)

            rows.append(SiteRow(label, ads_name, height, e_site, e_ads, e_ads_bsse, scf_ok, max_force))
            bsse_str = f"{e_ads_bsse:<16.6f}" if e_ads_bsse is not None else f"{'--':<16}"
            height_str = f"{height:<9.2f}" if height is not None else f"{'--':<9}"
            scf_str = color_text("WARN", 'yellow') if not scf_ok else "OK"
            force_str = f"{max_force:.4f}" if max_force is not None else "--"
            print_dual(f"{label:<30} {(ads_name or '--'):<12} {height_str} {e_ads:<14.6f} "
                        f"{bsse_str}{scf_str:<6}{force_str:<12}", f_out)
        print_dual("-" * len(header), f_out)
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

        print_dual(f"\n{color_text('[3] SUMMARY', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        print_dual(f"Sites analyzed : {len(rows)} (skipped: {n_skipped})", f_out)
        best = rows[0]
        verdict = "exothermic (favorable)" if best.e_ads < 0 else "endothermic (unfavorable)"
        print_dual(f"{color_text('Most stable site (uncorrected):', 'green')} {best.label}  "
                    f"(E_ads = {best.e_ads:.6f} eV, {verdict})", f_out)

        apply_source_label = best.label
        if n_bsse_found == 0:
            print_dual(color_text(
                "\n[NOTE] No BSSE-corrected results found -- re-run stb-adsorb with "
                "--bsse-correction (the CLI default) for a corrected reference.", 'yellow'), f_out)
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

        print_dual(f"{color_text('[Saved]', 'cyan')} Report     -> {REPORT_FILE}", f_out)

        if args.apply:
            print_dual(f"\n{color_text('[4] APPLY', 'magenta')}", f_out)
            print_dual("-" * 60, f_out)
            src = os.path.join(sites_root, apply_source_label, "structure.fdf")
            try:
                shutil.copy(src, args.apply)
            except OSError as e:
                print_dual(color_text(f"[ERROR] Could not copy '{src}' to '{args.apply}': {e}", 'red'), f_out)
            else:
                print_dual(f"{color_text('[Applied]', 'green')} {apply_source_label} -> {args.apply}", f_out)


if __name__ == "__main__":
    main()
