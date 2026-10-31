#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

from stb import __version__ as VERSION

import os
import re
import sys
import io
import glob
import time
import tarfile
import argparse
from datetime import datetime
from stb.core.cli import color_text, show_intro, print_dual, print_table, Sections

# Purpose: prep a SIESTA calc folder to restart a fresh run after a previous one that may or may not have finished
# OK -- delete every stale output/intermediate file (.out, .XV, .STRUCT_OUT, .bands, .DOS, .DM, .RHO, .ion(.xml),
# MESSAGES, ...), keep only what's needed to resubmit (pseudopotentials, the .fdf input, the submission script).
# Deliberately the OPPOSITE of an analysis-results archiver (that is stb-archive, 3.7): the default --keep list does
# NOT protect outputs other suite tools consume -- the job is to clear them out before the next run. Every removed
# file goes into a tar.gz next to the folder first, and the archive is read back and checked before anything is
# deleted (--no-archive skips both). --dry-run reports without touching anything.
#
# Never removed, whatever --keep says: hidden files (and hidden directories with --recursive), earlier stb-clean
# backups, and every file a kept .fdf reads (%include FILE, '%block X < FILE', 'Label < FILE') -- deleting one of
# those would break the very input the folder is being cleaned to rerun.
REPORT_FILE = "stb_clean_report.txt"
DEFAULT_KEEP = ['.psml', '.psf', '.fdf', '.sh']
BACKUP_RE = re.compile(r"_backup_\d{8}_\d{6}(_\d+)?\.tar\.gz$")
RECENT_SECONDS = 120          # an unfinished .out written this recently may belong to a running job

from stb.core.run_files import (  # noqa: F401  (shared with stb-archive)
    CAT_RESTART, CAT_LOG, CAT_JOB, CAT_GEOM, CAT_ELEC, CAT_GRID, CAT_BASIS, CAT_POST, CAT_SCRIPT,
    CAT_OTHER, CATEGORY_ORDER, classify, human_size, file_size, looks_like_job_script, read_fdf_settings, collect_files,
    fdf_key as _fdf_key, fdf_true as _fdf_true,
)


# ---------------------------------------------------------------- one folder ------------------------------------
def plan_folder(path, allowed_exts, recursive):
    """Decides, for every file of the folder, whether it stays (and why) or goes (and what it is)."""
    files = collect_files(path, recursive)
    keep, remove = {}, []
    for f in files:
        name = os.path.basename(f)
        ext = os.path.splitext(f)[1]
        if name.startswith('.'):
            keep[f] = "hidden file (never removed)"
        elif BACKUP_RE.search(name):
            keep[f] = "earlier stb-clean backup (never removed)"
        elif ext in allowed_exts:
            keep[f] = {".fdf": "SIESTA input (.fdf)", ".psf": "pseudopotential", ".psml": "pseudopotential",
                       ".sh": "job script" if looks_like_job_script(f) else "script"}.get(ext, f"kept by --keep ({ext})")
        else:
            remove.append(f)

    # files a kept .fdf reads: followed until nothing new is protected (an included file may include others)
    protected = {}
    while True:
        fdfs = [f for f in keep if f.endswith(".fdf")] + [f for f in protected if f.endswith(".fdf")]
        refs, labels, species = read_fdf_settings(fdfs)
        new = {os.path.normpath(r): src for r, src in refs.items()
               if os.path.normpath(r) in {os.path.normpath(x) for x in remove} and os.path.normpath(r) not in protected}
        if not new:
            break
        protected.update(new)
    for f in list(remove):
        if os.path.normpath(f) in protected:
            remove.remove(f)
            keep[f] = f"read by {os.path.relpath(protected[os.path.normpath(f)], path)}"
    missing_refs = sorted(os.path.relpath(r, path) for r in refs if not os.path.exists(r))
    return {"path": path, "keep": keep, "remove": remove, "labels": labels, "species": species,
            "fdfs": [f for f in keep if f.endswith(".fdf")], "missing_refs": missing_refs,
            "n_protected": len(protected)}


def previous_runs(path, files):
    """The SIESTA logs among `files` (any .out whose head names SIESTA): (relpath, finished, age in s, last line)."""
    runs = []
    for f in files:
        if not f.endswith(".out"):
            continue
        try:
            with open(f, "rb") as fh:
                head = fh.read(4000).decode(errors="replace")
                if "siesta" not in head.lower():
                    continue
                fh.seek(max(0, os.path.getsize(f) - 6000))
                tail = fh.read().decode(errors="replace")
        except OSError:
            continue
        finished = "Job completed" in tail or ">> End of run" in tail
        steps = None
        try:
            moves = re.findall(r"Begin (?:CG|FIRE|Broyden|MD)[^=\n]*=\s*(\d+)", open(f, errors="replace").read())
            steps = (int(moves[-1]) + 1) if moves else None
        except (OSError, ValueError):
            pass
        runs.append((os.path.relpath(f, path), finished, time.time() - os.path.getmtime(f), steps))
    return runs


def restart_rows(plan):
    """What the next run will start from: SIESTA's restart files vs the fdf flags that would read them."""
    label = plan["labels"].get("systemlabel") or "siesta"
    label = label.split()[0]
    rows = []
    removed = {os.path.basename(f) for f in plan["remove"]}
    kept = {os.path.basename(f) for f in plan["keep"]}
    for ext, flag, what, fresh in ((".DM", "DM.UseSaveDM", "density matrix", "a superposition of atomic densities"),
                                   (".XV", "MD.UseSaveXV", "positions + velocities", "the geometry in the .fdf"),
                                   (".CG", "MD.UseSaveCG", "CG optimizer history", "a new CG history")):
        fname = label + ext
        value = plan["labels"].get(_fdf_key(flag))
        on = value is not None and _fdf_true(value)
        flag_txt = f"{flag} {'true' if on else 'false (default)' if value is None else 'false'}"
        if fname in kept:
            state, effect = "kept", (f"next run reads it" if on else f"kept but not read ({flag} off)")
        elif fname in removed:
            state, effect = "removed", (f"starts from {fresh}" if on else f"starts from {fresh} (flag off anyway)")
        else:
            state, effect = "absent", f"starts from {fresh}"
        rows.append(([fname, what, state, flag_txt, effect], None))
    return label, rows


def archive_targets(path, targets):
    """tar.gz of `targets` (relative paths kept) in `path`'s parent directory. Returns the archive path."""
    parent = os.path.dirname(os.path.abspath(path.rstrip(os.sep))) or "."
    base = os.path.basename(os.path.abspath(path.rstrip(os.sep))) or "clean"
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    archive_path = os.path.join(parent, f"{base}_backup_{stamp}.tar.gz")
    n = 1
    while os.path.exists(archive_path):
        archive_path = os.path.join(parent, f"{base}_backup_{stamp}_{n}.tar.gz")
        n += 1
    with tarfile.open(archive_path, "w:gz") as tar:
        for target in targets:
            tar.add(target, arcname=os.path.relpath(target, path), recursive=False)
    return archive_path


def verify_archive(archive_path, path, targets):
    """Reads the archive back: every target must be in it, with its size. Returns None or a reason it failed."""
    try:
        with tarfile.open(archive_path, "r:gz") as tar:
            members = {m.name: m for m in tar.getmembers()}
    except (OSError, tarfile.TarError) as e:
        return f"the archive cannot be read back ({e})"
    for t in targets:
        m = members.get(os.path.relpath(t, path))
        if m is None:
            return f"{os.path.relpath(t, path)} is missing from the archive"
        if m.isfile() and m.size != os.lstat(t).st_size:
            return f"{os.path.relpath(t, path)} has the wrong size in the archive"
    return None


def confirm(n):
    try:
        answer = input(color_text(f"\nDelete these {n} file(s)? [y/N] ", 'yellow')).strip().lower()
    except EOFError:
        answer = ""
    return answer == 'y'


def clean_one(path, args, allowed_exts, sections, warn, f_out, title_suffix=""):
    """Analyses, reports, archives (and checks the archive), removes. Returns a summary record."""
    rec = {"dir": path, "kept": 0, "removed": 0, "freed": 0, "archive": None, "status": "ok"}
    plan = plan_folder(path, allowed_exts, args.recursive)
    keep, remove = plan["keep"], plan["remove"]
    rel = lambda f: os.path.relpath(f, path)  # noqa: E731
    all_files = list(keep) + remove
    total = sum(file_size(f) for f in all_files)
    rm_size = sum(file_size(f) for f in remove)
    rec["kept"] = len(keep)

    # ---- the folder and the run that left these files
    sections(f"FOLDER{title_suffix}")
    print_dual(f"Directory         : {os.path.abspath(path)}", f_out)
    print_dual(f"Files             : {len(all_files)} ({human_size(total)})"
               f"{' -- including subdirectories' if args.recursive else ''}", f_out)
    print_dual(f"Keep / remove     : {len(keep)} kept ({human_size(total - rm_size)}), "
               f"{len(remove)} to remove ({human_size(rm_size)})", f_out)
    label = (plan["labels"].get("systemlabel") or "siesta").split()[0]
    if plan["fdfs"]:
        print_dual(f"SIESTA input      : {', '.join(rel(f) for f in plan['fdfs'])}; SystemLabel {label}"
                   f"{' (default)' if 'systemlabel' not in plan['labels'] else ''}", f_out)
    runs = previous_runs(path, all_files)
    if runs:
        for name, finished, age, steps in runs:
            when = f"last written {age / 60:.0f} min ago" if age < 7200 else f"last written {age / 86400:.1f} days ago" \
                if age > 86400 else f"last written {age / 3600:.1f} h ago"
            extra = f", {steps} geometry steps" if steps else ""
            print_dual(f"Previous run      : {name} -- {'finished normally' if finished else 'did NOT finish'}"
                       f"{extra} ({when})", f_out)
            if not finished and age < RECENT_SECONDS:
                warn(f"{path}: {name} did not finish and was written {age:.0f} s ago -- is the job still running? "
                     f"Cleaning would delete files it is writing.")
            if finished and args.no_archive and not args.dry_run and any(rel(f) == name for f in remove):
                warn(f"{path}: the results of a FINISHED run ({name}) are removed with --no-archive -- no copy kept.")
    elif any(f.endswith(".out") for f in all_files):
        print_dual("Previous run      : no SIESTA log among the .out files", f_out)
    else:
        print_dual("Previous run      : none found (no SIESTA .out log)", f_out)
    if os.path.isfile(os.path.join(path, "0_NORMAL_EXIT")):
        print_dual("Exit marker       : 0_NORMAL_EXIT present (SIESTA ended normally)", f_out)

    # ---- what stays
    sections(f"KEPT{title_suffix}")
    rows = [([rel(f), human_size(file_size(f)), why], 'cyan' if why.startswith("read by") else None)
            for f, why in sorted(keep.items(), key=lambda kv: rel(kv[0]))]
    if rows:
        shown = rows if len(rows) <= 30 else rows[:30]
        print_table(["file", "size", "why it stays"], shown, f_out)
        if len(rows) > 30:
            print_dual(f"... and {len(rows) - 30} more", f_out)
    else:
        print_dual("nothing is kept", f_out)
    if plan["n_protected"]:
        print_dual(f"Protected         : {plan['n_protected']} file(s) read by a kept .fdf, kept although their "
                   f"extension is not in --keep", f_out)
    if not plan["fdfs"]:
        warn(f"{path}: no .fdf input is kept -- there is nothing to rerun from (check --keep).")
    for r in plan["missing_refs"]:
        warn(f"{path}: a kept .fdf reads '{r}', which is not in the folder.")
    psp = {os.path.splitext(os.path.basename(f))[0] for f in keep if f.endswith((".psf", ".psml"))}
    lacking = [s for s in plan["species"] if s not in psp]
    if plan["species"]:
        print_dual(f"Species           : {', '.join(plan['species'])} -- pseudopotential kept for "
                   f"{', '.join(s for s in plan['species'] if s in psp) or 'none'}", f_out)
    if lacking:
        warn(f"{path}: no {'/'.join(lacking)} .psf/.psml is kept in the folder -- fine only if the job script "
             f"or a pseudopotential directory provides it.")

    # ---- what goes
    sections(f"REMOVED{title_suffix}")
    if not remove:
        print_dual("nothing to remove", f_out)
    else:
        by_cat = {}
        for f in remove:
            cat = classify(os.path.basename(f))
            if cat in (CAT_OTHER, CAT_POST) and looks_like_job_script(f):
                cat = CAT_SCRIPT
            by_cat.setdefault(cat, []).append(f)
        rows = []
        for cat in CATEGORY_ORDER:
            fs = by_cat.get(cat)
            if not fs:
                continue
            uniq = list(dict.fromkeys(os.path.basename(f) for f in fs))
            names = ", ".join(uniq[:3]) + (f", +{len(uniq) - 3}" if len(uniq) > 3 else "")
            rows.append(([cat, str(len(fs)), human_size(sum(file_size(f) for f in fs)), names], None))
        rows.append(([f"total", str(len(remove)), human_size(rm_size), ""], 'bold'))
        print_table(["kind", "files", "size", "e.g."], rows, f_out)
        big = sorted(remove, key=file_size, reverse=True)[:5]
        print_dual("", f_out)
        print_dual("Largest           : " + ", ".join(f"{rel(f)} ({human_size(file_size(f))})" for f in big), f_out)
        scripts = [rel(f) for f in by_cat.get(CAT_SCRIPT, [])]
        for s in scripts:
            warn(f"{path}: '{s}' looks like a job script and will be removed -- add its extension to --keep "
                 f"(--keep {' '.join(sorted(allowed_exts))} {os.path.splitext(s)[1] or s}).")

    # ---- what the next run starts from
    sections(f"RESTART{title_suffix}")
    label, rows = restart_rows(plan)
    print_table(["file", "holds", "after cleaning", "fdf flag", "next run"], rows, f_out)
    moved = max((s or 0) for _, _, _, s in runs) if runs else 0
    gone = [n for n in (f"{label}.XV", f"{label}.STRUCT_OUT")
            if n in {os.path.basename(f) for f in remove}]
    if moved > 1 and gone:
        print_dual(color_text(f"Note              : the previous run moved the geometry ({moved} steps); its final "
                              f"geometry is in {' and '.join(gone)}, which go -- the next run starts again from the "
                              f".fdf geometry. Keep .XV, or convert {gone[-1]} to an .fdf first (stb-translate, 3.1).",
                              'cyan'), f_out)
    print_dual("To continue from the previous run instead of starting fresh, keep its restart files, e.g. "
               f"--keep {' '.join(sorted(allowed_exts))} .DM .XV", f_out)

    # ---- do it
    sections(f"ARCHIVE & REMOVAL{title_suffix}")
    rec["freed"] = rm_size
    if not remove:
        print_dual("Nothing to do.", f_out)
        rec["status"] = "nothing to do"
        return rec
    if args.dry_run:
        print_dual(color_text(f"Dry run -- nothing archived, nothing removed ({len(remove)} file(s), "
                              f"{human_size(rm_size)} would go).", 'cyan'), f_out)
        rec["status"], rec["removed"] = "dry run", len(remove)
        return rec
    if not args.no_confirm:
        preview = remove if len(remove) <= 20 else remove[:20]
        print(f"\n{len(remove)} file(s) will be removed:")
        for t in preview:
            print(f"  {rel(t)}")
        if len(remove) > 20:
            print(f"  ... and {len(remove) - 20} more")
        if not confirm(len(remove)):
            print_dual("Cancelled -- nothing archived, nothing removed.", f_out)
            rec["status"], rec["freed"] = "cancelled", 0
            return rec
    if not args.no_archive:
        archive = archive_targets(path, remove)
        problem = verify_archive(archive, path, remove)
        a_size = file_size(archive)
        print_dual(f"Archive           : {archive}", f_out)
        print_dual(f"                    {len(remove)} file(s), {human_size(a_size)} "
                   f"({100 * a_size / max(rm_size, 1):.0f} % of {human_size(rm_size)})", f_out)
        if problem:
            print_dual(color_text(f"Verification      : FAILED -- {problem}; nothing removed", 'red'), f_out)
            warn(f"{path}: archive check failed ({problem}) -- nothing was removed.")
            rec["status"], rec["freed"], rec["archive"] = "FAILED", 0, archive
            return rec
        print_dual(color_text("Verification      : read back, every file present with its size -- ok", 'green'),
                   f_out)
        print_dual(f"Restore with      : tar -xzf {archive} -C {os.path.abspath(path)}", f_out)
        rec["archive"] = archive
    else:
        print_dual(color_text("Archive           : none (--no-archive) -- the removed files cannot be recovered",
                              'yellow'), f_out)
    removed = failed = 0
    for t in remove:
        try:
            os.remove(t)
            removed += 1
        except OSError as e:
            failed += 1
            warn(f"{path}: could not remove {rel(t)}: {e}")
    print_dual(f"Removed           : {removed} file(s), {human_size(rm_size)} freed"
               f"{f'; {failed} failed' if failed else ''}", f_out)
    rec["removed"] = removed
    if failed:
        rec["status"] = "FAILED"
    return rec


def build_parser():
    parser = argparse.ArgumentParser(
        description="Clean a SIESTA calc folder to restart a fresh run: removes every file except the extensions "
                    "in --keep (pseudopotentials, .fdf input, submission script by default) -- outputs/results too "
                    "(.out, .XV, .DM, .bands, ...): it prepares a re-run, it does not archive results (that is "
                    "stb-archive). Never removed: hidden files, earlier stb-clean backups and every file a kept .fdf "
                    "reads (%include, '< file'). Every removed file is archived (tar.gz, next to the folder) and "
                    "the archive is checked before anything is deleted, unless --no-archive. The report says what "
                    "the previous run was, what stays and why, what goes by kind, and what the next run will start "
                    "from (.DM/.XV/.CG vs the fdf UseSave flags)."
    )
    parser.add_argument('--path', default='.',
                        help="Directory to clean, or a glob pattern (e.g. 'strain_*') matching several directories "
                             "to clean one after another (default: current directory).")
    parser.add_argument('--keep', nargs='+', default=list(DEFAULT_KEEP),
                        help="Extensions to keep, i.e. NOT delete (default: .psml .psf .fdf .sh -- pseudopotentials, "
                             "the SIESTA input, and the submission script). Everything else, including prior run "
                             "outputs, is removed. Add .DM .XV to continue from the previous run.")
    parser.add_argument('--recursive', action='store_true',
                        help="Also clean subdirectories (skips hidden ones, e.g. .git), not just the top level of "
                             "--path.")
    parser.add_argument('--dry-run', action='store_true',
                        help="Report what would be kept and removed, and touch nothing (no archive either).")
    parser.add_argument('--no-archive', action='store_true',
                        help="Skip archiving removed files before deleting them. Archiving (tar.gz, next to --path, "
                             "read back before anything is deleted) is the default safety net -- turn it off only "
                             "if you're sure.")
    parser.add_argument('--no-confirm', '-n', action='store_true',
                        help="Don't ask for confirmation before deleting.")
    parser.add_argument('--save-report', action='store_true',
                        help=f"Also persist the report to {REPORT_FILE} (in --path; in the current directory for a "
                             f"glob). Off by default.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-clean {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

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
    batch = len(matches) > 1

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    report_path = os.path.join(matches[0] if not batch else ".", REPORT_FILE) if args.save_report else None
    # Buffered in memory and written only after every folder is cleaned: the report sits inside a folder this run
    # cleans, and a .txt is not kept by default -- the tool would delete its own report.
    f_out = io.StringIO() if report_path else None
    warns = []
    sections = Sections(f_out)

    print("\n" + color_text("CLEAN:", 'bold'))
    print("-" * 60)
    print_dual(color_text("===== STB-CLEAN REPORT =====", 'magenta'), f_out)
    sections("RUN METADATA")
    print_dual(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"stb               : {VERSION}", f_out)
    if is_glob:
        print_dual(f"Target            : '{args.path}' -- {len(matches)} director{'ies' if batch else 'y'} matched",
                   f_out)
    else:
        print_dual(f"Target            : {os.path.abspath(matches[0])}", f_out)
    print_dual(f"Keep extensions   : {' '.join(sorted(allowed_exts))}"
               f"{' (default)' if sorted(allowed_exts) == sorted(DEFAULT_KEEP) else ''}", f_out)
    print_dual("Always kept       : hidden files, earlier stb-clean backups, files a kept .fdf reads", f_out)
    print_dual(f"Recursive         : {'yes (hidden directories skipped)' if args.recursive else 'no'}", f_out)
    mode = "dry run (nothing touched)" if args.dry_run else \
        ("delete without backup (--no-archive)" if args.no_archive else "archive (tar.gz, checked), then delete")
    print_dual(f"Mode              : {mode}{'; no confirmation (-n)' if args.no_confirm and not args.dry_run else ''}",
               f_out)

    records = []
    for d in matches:
        suffix = f" -- {d}" if batch else ""
        records.append(clean_one(d, args, allowed_exts, sections, warns.append, f_out, suffix))

    sections("SUMMARY & FILES")
    if batch:
        rows = [([r["dir"], str(r["kept"]), str(r["removed"]), human_size(r["freed"]),
                  os.path.basename(r["archive"]) if r["archive"] else "-", r["status"]],
                 'red' if r["status"] == "FAILED" else None) for r in records]
        print_table(["folder", "kept", "removed" if not args.dry_run else "to remove", "size", "archive", "status"],
                    rows, f_out)
        print_dual("", f_out)
    n_rm = sum(r["removed"] for r in records)
    freed = sum(r["freed"] for r in records)
    if args.dry_run:
        print_dual(f"Would remove      : {n_rm} file(s), {human_size(freed)} -- rerun without --dry-run", f_out)
    else:
        print_dual(f"Removed           : {n_rm} file(s), {human_size(freed)} freed", f_out)
        for r in records:
            if r["archive"]:
                print_dual(f"Archive           : {r['archive']}", f_out)
    if report_path:
        print_dual(f"Report            : {report_path}", f_out)
    sections("WARNINGS")
    if warns:
        for i, w in enumerate(warns, 1):
            print_dual(color_text(f"  {i}. {w}", 'yellow'), f_out)
    else:
        print_dual(color_text("  none", 'green'), f_out)

    if report_path:
        with open(report_path, "w") as f:
            f.write(f_out.getvalue())

    print("\n" + "-" * 60)
    if n_rm and not args.dry_run:
        print(color_text("Cleaned up! May the deleted files rest in pieces.\n", 'bold'))
    else:
        print(color_text("Nothing removed.\n", 'bold'))

    if any(r["status"] == "FAILED" for r in records):
        sys.exit(1)


if __name__ == "__main__":
    main()
