#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.9.1"

import os
import sys
import io
import glob
import tarfile
import argparse
from datetime import datetime
from stb.core.cli import COLORS, color_text, show_intro, print_dual, print_section

# Purpose: prep a SIESTA calc folder to restart a fresh run after a previous
# one that may or may not have finished OK -- delete every stale output/
# intermediate file (.out, .XV, .STRUCT_OUT, .bands, .DOS, .DM, .RHO,
# .ion(.xml), MESSAGES, ...), keep only what's needed to resubmit
# (pseudopotentials, the .fdf input, the submission script). Deliberately the
# OPPOSITE of an analysis-results archiver: the default --keep list does NOT
# protect outputs other suite tools consume (stb-bands/-dos/-structural,
# core.siesta_log.check_scf_and_force, ...) -- that's by design, since this
# tool's job is to clear them out before the next run. Every removed file is
# archived (tar.gz, next to --path) before deletion unless --no-archive is
# given -- that safety net is what makes "just delete it, no per-file
# confirmation" a reasonable default, so there's no --dry-run: the archive
# already gives you a no-risk way to see (and recover) exactly what was
# removed.
REPORT_FILE = "stb_clean_report.txt"


def should_delete(file_path, allowed_exts):
    return os.path.isfile(file_path) and os.path.splitext(file_path)[1] not in allowed_exts


def collect_targets(path, allowed_exts, recursive):
    """Returns the sorted list of file paths (under `path`) that should be
    removed. Recursive walk skips hidden directories (.git and friends) --
    never something this tool should be reaching into."""
    targets = []
    if recursive:
        for root, dirs, files in os.walk(path):
            dirs[:] = [d for d in dirs if not d.startswith('.')]
            for f in files:
                full = os.path.join(root, f)
                if should_delete(full, allowed_exts):
                    targets.append(full)
    else:
        for f in os.listdir(path):
            full = os.path.join(path, f)
            if should_delete(full, allowed_exts):
                targets.append(full)
    return sorted(targets)


def archive_targets(path, targets, f_out=None):
    """tar.gz's `targets` (paths under `path`) into a timestamped archive in
    `path`'s own parent directory, preserving each file's path relative to
    `path`. Returns the archive path."""
    parent = os.path.dirname(os.path.abspath(path.rstrip(os.sep))) or "."
    base = os.path.basename(os.path.abspath(path.rstrip(os.sep))) or "clean"
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    archive_path = os.path.join(parent, f"{base}_backup_{timestamp}.tar.gz")
    with tarfile.open(archive_path, "w:gz") as tar:
        for target in targets:
            tar.add(target, arcname=os.path.relpath(target, path))
    print_dual(f"[OK] Archived {len(targets)} file(s) to {archive_path}", f_out)
    return archive_path


def remove_targets(targets, no_confirm, f_out=None):
    """Deletes `targets` -- one confirmation for the whole list (not one per
    file), unless --no-confirm. Returns (removed, cancelled, failed) counts;
    a single OSError (permission, etc.) is reported and skipped rather than
    aborting the rest of the list."""
    if not no_confirm:
        print_dual(f"\n{len(targets)} file(s) will be removed:", f_out)
        preview = targets if len(targets) <= 20 else targets[:20] + [f"... and {len(targets) - 20} more"]
        for t in preview:
            print_dual(f"  {t}", f_out)
        answer = input(color_text(f"\nDelete these {len(targets)} file(s)? [y/N] ", 'yellow')).strip().lower()
        if answer != 'y':
            print_dual("[INFO] Cancelled -- no files removed.", f_out)
            return 0, len(targets), 0

    removed = failed = 0
    for target in targets:
        try:
            os.remove(target)
            removed += 1
        except OSError as e:
            print_dual(color_text(f"[WARNING] Could not remove {target}: {e}", 'yellow'), f_out)
            failed += 1
    return removed, 0, failed


def clean_one(path, allowed_exts, recursive, archive, no_confirm, f_out=None):
    """Cleans a single directory: collect -> archive -> remove -> report.
    Returns (removed, failed) counts."""
    targets = collect_targets(path, allowed_exts, recursive)
    if not targets:
        print_dual("[OK] Nothing to remove.", f_out)
        return 0, 0

    if archive:
        archive_targets(path, targets, f_out)

    removed, cancelled, failed = remove_targets(targets, no_confirm, f_out)
    print_dual(f"[OK] Removed {removed}, cancelled {cancelled}, failed {failed} "
               f"(of {len(targets)} candidate(s)).", f_out)
    return removed, failed


def main():
    parser = argparse.ArgumentParser(
        description="Clean a SIESTA calc folder to restart a fresh run: removes every file "
                     "except the extensions in --keep (pseudopotentials, .fdf input, submission "
                     "script by default). This deletes outputs/results too (.out, .XV, .bands, "
                     ".DOS, ...) -- it's meant to prep for a re-run, not to archive results. "
                     "Every removed file is archived (tar.gz) first unless --no-archive."
    )
    parser.add_argument(
        '--path', default='.',
        help="Directory to clean, or a glob pattern (e.g. 'strain_*') matching several "
             "directories to clean one after another (default: current directory)."
    )
    parser.add_argument(
        '--keep', nargs='+', default=['.psml', '.psf', '.fdf', '.sh'],
        help="Extensions to keep, i.e. NOT delete (default: .psml .psf .fdf .sh -- "
             "pseudopotentials, the SIESTA input, and the submission script). Everything "
             "else, including prior run outputs, is removed."
    )
    parser.add_argument(
        '--recursive', action='store_true',
        help="Also clean subdirectories (skips hidden ones, e.g. .git), not just the top "
             "level of --path."
    )
    parser.add_argument(
        '--no-archive', action='store_true',
        help="Skip archiving removed files before deleting them. Archiving (tar.gz, next "
             "to --path) is the default safety net -- turn it off only if you're sure."
    )
    parser.add_argument(
        '--no-confirm', '-n', action='store_true',
        help="Don't ask for confirmation before deleting."
    )
    parser.add_argument(
        '--save-report', action='store_true',
        help=f"Also persist the report to {REPORT_FILE} (next to --path). Off by default."
    )
    parser.add_argument("-v", "--version", action="version",
                        version=f"stb-clean {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("CLEAN:", 'bold'))
    print("-" * 60)

    for ext in args.keep:
        if not ext.startswith('.'):
            parser.error(f"--keep entries must start with '.' -- got '{ext}' (did you mean '.{ext}'?)")
    allowed_exts = set(args.keep)

    is_glob = any(c in args.path for c in '*?[')
    matches = sorted(p for p in glob.glob(args.path) if os.path.isdir(p)) if is_glob else [args.path]
    if not matches:
        parser.error(f"No directory matched '{args.path}'.")
    if not is_glob and not os.path.isdir(args.path):
        parser.error(f"'{args.path}' is not a directory.")
    batch_mode = len(matches) > 1

    report_path = os.path.join(matches[0] if not batch_mode else ".", REPORT_FILE) if args.save_report else None
    # NOTE: buffered in memory, NOT written to disk until every directory has
    # already been cleaned -- report_path sits inside a directory this same
    # run is about to scan for deletion candidates, and a .txt extension
    # isn't in any default --keep list, so opening the real file up front
    # would make the tool delete its own still-being-written report.
    f_out = io.StringIO() if report_path else None

    print_dual(color_text("===== STB-CLEAN REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    if batch_mode:
        print_dual(f"Target pattern    : '{args.path}' ({len(matches)} director(y/ies) matched)", f_out)
    else:
        print_dual(f"Target directory  : {matches[0]}", f_out)
    print_dual(f"Keep extensions   : {' '.join(sorted(allowed_exts))}", f_out)
    print_dual(f"Recursive         : {'yes' if args.recursive else 'no'}", f_out)
    print_dual(f"Archive backup    : {'no (--no-archive)' if args.no_archive else 'yes'}", f_out)
    if report_path:
        print_dual(f"Report file       : {report_path}", f_out)

    print_section("[1] CLEANING", f_out)
    total_removed = total_failed = 0
    for target_dir in matches:
        if batch_mode:
            print_dual(f"\n--- {target_dir} ---", f_out)
        removed, failed = clean_one(
            target_dir, allowed_exts, args.recursive, not args.no_archive,
            args.no_confirm, f_out)
        total_removed += removed
        total_failed += failed

    print_section("[2] SUMMARY & FILES", f_out)
    print_dual(f"Total removed     : {total_removed}", f_out)
    if total_failed:
        print_dual(color_text(f"Total failed      : {total_failed}", 'yellow'), f_out)
    if report_path:
        print_dual(f"Report            : {report_path}", f_out)

    if report_path:
        with open(report_path, "w") as f:
            f.write(f_out.getvalue())

    print("\n" + "-" * 60)
    print(color_text("Cleaned up! May the deleted files rest in pieces.\n", 'bold'))

    if total_failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
