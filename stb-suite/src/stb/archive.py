#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.0.0"

import os
import io
import sys
import glob
import tarfile
import argparse
from datetime import datetime

from stb.core import siesta_log, structure_io
from stb.core.cli import COLORS, color_text, show_intro, print_dual, print_section

REPORT_FILE = "stb_archive_report.txt"

# The complement of stb-clean's default --keep list: inputs (needed to
# resubmit) PLUS the essential outputs stb-clean deliberately deletes (since
# THAT tool preps for a restart, not archiving). .xml catches PDOS.xml
# (and any other stray XML) via its extension, same os.path.splitext-based
# matching as stb-clean's should_delete. Deliberately excludes large
# regenerable intermediates (.HSX, .DM, .RHO, .VT, .VH, MESSAGES, ...) --
# --include adds any of those back in if a specific archive genuinely needs
# them (e.g. .RHO for a Bader/cube analysis someone downstream wants to redo
# without rerunning SIESTA).
DEFAULT_INCLUDE_EXTS = ['.fdf', '.psf', '.psml', '.sh', '.out', '.XV', '.STRUCT_OUT',
                        '.bands', '.DOS', '.xml']


def should_archive(file_path, include_exts):
    return os.path.isfile(file_path) and os.path.splitext(file_path)[1] in include_exts


def build_manifest(path, label, files):
    """Returns the MANIFEST.txt content (list of lines) embedded in every
    archive: what this calculation is (SystemLabel, run type, convergence,
    energy -- the same core.siesta_log fields stb-status reports, read
    quietly here instead of via that tool's own printed report) and exactly
    which files got packaged.
    """
    lines = [
        "STB-ARCHIVE MANIFEST",
        "=" * 60,
        f"Created           : {datetime.now():%Y-%m-%d %H:%M:%S}",
        f"Source directory  : {os.path.abspath(path)}",
    ]
    if label:
        lines.append(f"SystemLabel       : {label}")
    out_file = siesta_log.find_out_file(path, label)
    if out_file:
        dynamics = siesta_log.get_dynamics_type(out_file)
        category = siesta_log.categorize_dynamics(dynamics)
        scf_ok, _ = siesta_log.get_scf_convergence(out_file)
        max_force = siesta_log.get_max_force(out_file)
        free_energy = siesta_log.get_free_energy(out_file)
        lines.append(f"Run type          : {category}" + (f" ({dynamics})" if dynamics else ""))
        if scf_ok is not None:
            lines.append(f"SCF converged     : {'yes' if scf_ok else 'no'}")
        if max_force is not None:
            lines.append(f"Max force (eV/Ang): {max_force:.6f}")
        if free_energy is not None:
            lines.append(f"Free energy (eV)  : {free_energy:.6f}")
    lines.append("")
    lines.append(f"Files included ({len(files)}):")
    for f in files:
        lines.append(f"  - {f}")
    return lines


def archive_one(path, include_exts, label_override, output_override, f_out=None):
    """Packages one directory's matching files (+ a generated MANIFEST.txt)
    into a tar.gz. Returns (ok, archive_path_or_None, n_files)."""
    label = structure_io.find_fdf_system_label(path, label_override)
    files = sorted(f for f in os.listdir(path) if should_archive(os.path.join(path, f), include_exts))

    if not files:
        print_dual(color_text(
            "[WARNING] No matching files found (check --include/--exclude) -- "
            "nothing to archive.", 'yellow'), f_out)
        return False, None, 0

    if output_override:
        archive_path = output_override
    else:
        base = os.path.basename(os.path.abspath(path.rstrip(os.sep))) or "archive"
        parent = os.path.dirname(os.path.abspath(path.rstrip(os.sep))) or "."
        archive_path = os.path.join(parent, f"{base}_archive.tar.gz")

    manifest_bytes = "\n".join(build_manifest(path, label, files)).encode()

    with tarfile.open(archive_path, "w:gz") as tar:
        for f in files:
            tar.add(os.path.join(path, f), arcname=f)
        info = tarfile.TarInfo(name="MANIFEST.txt")
        info.size = len(manifest_bytes)
        tar.addfile(info, io.BytesIO(manifest_bytes))

    print_dual(f"[OK] Archived {len(files)} file(s) + MANIFEST.txt to {archive_path}", f_out)
    return True, archive_path, len(files)


def main():
    parser = argparse.ArgumentParser(
        description="Package a finished SIESTA calculation (inputs + essential outputs) "
                     "into a single .tar.gz for sharing or reproducibility -- the complement "
                     "of stb-clean, which deletes that same class of outputs to prep for a "
                     "restart instead of preserving them. --path can be a glob pattern to "
                     "archive several folders at once (e.g. 'strain_*'), one archive each."
    )
    parser.add_argument(
        "--path", default=".",
        help="Directory to archive, or a glob pattern (e.g. 'strain_*') matching several "
             "directories (default: current directory)."
    )
    parser.add_argument(
        "-l", "--label", default=None,
        help="SystemLabel override, applied to every matched directory. Auto-detected per "
             "directory from its .fdf file(s) if not given (only affects the manifest's "
             "run-type/convergence/energy summary, not which files get archived)."
    )
    parser.add_argument(
        "--include", nargs="+", default=None, metavar="EXT",
        help="Extra extensions to include, on top of the default set "
             f"({' '.join(DEFAULT_INCLUDE_EXTS)})."
    )
    parser.add_argument(
        "--exclude", nargs="+", default=None, metavar="EXT",
        help="Extensions to remove from the default include set."
    )
    parser.add_argument(
        "-o", "--output", default=None,
        help="Archive filename. Default: <foldername>_archive.tar.gz (next to --path). "
             "Ignored in batch mode (multiple directories matched --path)."
    )
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-archive {VERSION}")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - Archive Packager",
            "Bundles a finished calculation for sharing, not deleting",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("ARCHIVE:", 'bold'))
    print("-" * 60)

    include_exts = set(DEFAULT_INCLUDE_EXTS)
    for ext in (args.exclude or []):
        if not ext.startswith('.'):
            parser.error(f"--exclude entries must start with '.' -- got '{ext}' (did you mean '.{ext}'?)")
        include_exts.discard(ext)
    for ext in (args.include or []):
        if not ext.startswith('.'):
            parser.error(f"--include entries must start with '.' -- got '{ext}' (did you mean '.{ext}'?)")
        include_exts.add(ext)

    is_glob = any(c in args.path for c in '*?[')
    matches = sorted(p for p in glob.glob(args.path) if os.path.isdir(p)) if is_glob else [args.path]
    if not matches:
        parser.error(f"No directory matched '{args.path}'.")
    if not is_glob and not os.path.isdir(args.path):
        parser.error(f"'{args.path}' is not a directory.")
    batch_mode = len(matches) > 1

    if batch_mode and args.output:
        parser.error(
            "--output can't be used when --path matches multiple directories "
            f"({len(matches)} matched '{args.path}'); each gets its own auto-named archive.")

    report_path = REPORT_FILE if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    print_dual(color_text("===== STB-ARCHIVE REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    if batch_mode:
        print_dual(f"Target pattern    : '{args.path}' ({len(matches)} director(y/ies) matched)", f_out)
    else:
        print_dual(f"Target directory  : {matches[0]}", f_out)
    print_dual(f"Include extensions: {' '.join(sorted(include_exts))}", f_out)
    if report_path:
        print_dual(f"Report file       : {report_path}", f_out)

    print_section("[1] ARCHIVING", f_out)
    total_ok = total_failed = 0
    archives = []
    for target_dir in matches:
        if batch_mode:
            print_dual(f"\n--- {target_dir} ---", f_out)
        ok, archive_path, _ = archive_one(target_dir, include_exts, args.label, args.output, f_out)
        if ok:
            total_ok += 1
            archives.append(archive_path)
        else:
            total_failed += 1

    print_section("[2] SUMMARY & FILES", f_out)
    print_dual(f"Archived          : {total_ok}", f_out)
    if total_failed:
        print_dual(color_text(f"Skipped (empty)   : {total_failed}", 'yellow'), f_out)
    for a in archives:
        print_dual(f"  - {a}", f_out)
    if report_path:
        print_dual(f"Report            : {report_path}", f_out)

    if f_out:
        f_out.close()

    print("\n" + "-" * 60)
    print(color_text("Archived -- science preserved, not deleted.\n", 'bold'))

    if total_ok == 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
