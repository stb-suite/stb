#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

from stb import __version__ as VERSION

import os
import io
import re
import sys
import glob
import hashlib
import tarfile
import argparse
from datetime import datetime

from stb.core import structure_io
from stb.core.run_files import (
    CAT_RESTART, CAT_LOG, CAT_JOB, CAT_GEOM, CAT_ELEC, CAT_GRID, CAT_BASIS, CAT_POST, CAT_OTHER,
    classify, human_size, file_size, looks_like_job_script, read_fdf_settings, collect_files,
)
from stb.core.cli import color_text, show_intro, print_dual, print_table, Sections

# Purpose: package a finished SIESTA calculation for sharing or for the record -- the complement of stb-clean (3.2),
# which removes outputs to prepare a rerun. Every file of the folder is sorted by kind (core/run_files: SIESTA's naming
# conventions) and either goes into the archive or is left out with a reason; what the input needs to run again
# (the .fdf files, every file they read, the pseudopotentials, the job script) always goes in. The archive holds a
# MANIFEST.txt (what the run was, from stb-status, and every file with its SHA-256) and a SHA256SUMS that
# `sha256sum -c` checks after extraction; it is read back and every file's checksum compared before the report says ok.
# An existing archive is never overwritten.

REPORT_FILE = "stb_archive_report.txt"
OTHER_LIMIT = 1_000_000                  # an unclassified file goes in only below this size (bytes) in 'results' mode
BIG_BINARY = {".HSX", ".HS", ".WFSX", ".WFS", ".TSHS", ".TSDE"}
ARCHIVE_RE = re.compile(r"(_archive(_\d+)?|_backup_\d{8}_\d{6}(_\d+)?)\.tar\.gz$")
CAT_INPUT = "input (needed to run again)"

MODES = {
    "results": "inputs + logs + geometry/trajectory + small electronic outputs (bands, DOS, PDOS, eigenvalues) + "
               "post-processing; no grids, restart files, basis tables or H/S/wavefunction binaries",
    "inputs": "only what the run needs to be repeated: .fdf files, every file they read, pseudopotentials, job scripts",
    "full": "every file (except hidden files and earlier archives)",
}


# ---------------------------------------------------------------------------------------------------- selection --
def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


MULTI_SUFFIXES = (".ion.xml", ".ion.nc", ".grid.nc", ".PDOS.xml", ".DM.nc")


def suffix(name):
    """What --include needs to match a file: a compound SIESTA suffix, else the extension, else the name."""
    for sfx in MULTI_SUFFIXES:
        if name.endswith(sfx):
            return sfx
    return os.path.splitext(name)[1] or os.path.basename(name)


def plan_archive(path, mode, recursive, include, exclude):
    """Every file of the folder: (relpath, kind, size, why) for what goes in and for what is left out."""
    files = [f for f in collect_files(path, recursive) if not os.path.basename(f).startswith(".")
             and not ARCHIVE_RE.search(os.path.basename(f))]
    # inputs: the .fdf files, everything they read (followed through included files), pseudopotentials, scripts
    fdfs = [f for f in files if f.endswith(".fdf")]
    refs, labels, species = read_fdf_settings(fdfs)
    seen = set()
    while True:
        new = [r for r in refs if r.endswith(".fdf") and r not in seen and os.path.isfile(r)]
        if not new:
            break
        seen.update(new)
        refs2, _, _ = read_fdf_settings(fdfs + sorted(seen))
        refs.update(refs2)
    read_by = {os.path.normpath(r): src for r, src in refs.items()}
    missing_refs = sorted(os.path.relpath(r, path) for r in read_by if not os.path.exists(r))

    def matches(name, exts):
        return any(name.endswith(e) for e in exts)

    keep, leave = [], []
    for f in files:
        rel = os.path.relpath(f, path)
        name = os.path.basename(f)
        ext = os.path.splitext(name)[1]
        size = file_size(f)
        kind = classify(name)
        reason_in = None
        if name.endswith(".fdf"):
            reason_in, kind = "SIESTA input", CAT_INPUT
        elif ext in (".psf", ".psml"):
            reason_in, kind = "pseudopotential", CAT_INPUT
        elif os.path.normpath(f) in read_by:
            reason_in, kind = f"read by {os.path.relpath(read_by[os.path.normpath(f)], path)}", CAT_INPUT
        elif kind in (CAT_OTHER, CAT_POST) and looks_like_job_script(f):
            reason_in, kind = "job script", CAT_INPUT
        if matches(name, exclude) and kind != CAT_INPUT:
            leave.append((rel, kind, size, "--exclude"))
            continue
        if reason_in or matches(name, include):
            keep.append((rel, kind, size, reason_in or "--include"))
            continue
        if mode == "full":
            keep.append((rel, kind, size, "full mode"))
            continue
        if mode == "inputs":
            leave.append((rel, kind, size, "inputs mode: not an input"))
            continue
        # results mode
        why_out = None
        if kind == CAT_RESTART and ext != ".XV":
            why_out = "restart file (regenerable by rerunning)"
        elif kind == CAT_GRID:
            why_out = "real-space grid (large)"
        elif kind == CAT_BASIS:
            why_out = "basis table (SIESTA regenerates it)"
        elif ext in BIG_BINARY:
            why_out = "H/S or wavefunction binary (large)"
        elif kind == CAT_OTHER and size >= OTHER_LIMIT:
            why_out = f"unknown kind, over {human_size(OTHER_LIMIT)}"
        if why_out:
            leave.append((rel, kind, size, why_out))
        else:
            keep.append((rel, kind, size, "result"))
    label = (labels.get("systemlabel") or "").split()[0] if labels.get("systemlabel") else None
    present_pseudo = {os.path.splitext(os.path.basename(r[0]))[0] for r in keep
                      if r[0].endswith((".psf", ".psml"))}
    missing_pseudo = [sp for sp in species if sp not in present_pseudo]
    subdirs = [] if recursive else sorted(d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))
                                          and not d.startswith("."))
    return dict(keep=keep, leave=leave, label=label, species=species, missing_refs=missing_refs,
                missing_pseudo=missing_pseudo, n_files=len(files), subdirs=subdirs)


def archive_name(path, output, force):
    """The archive path: -o, else <folder>_archive.tar.gz next to the folder; never an existing file unless --force."""
    if output:
        if os.path.exists(output) and not force:
            return None, f"'{output}' exists (pass --force to overwrite)"
        return output, None
    parent = os.path.dirname(os.path.abspath(path.rstrip(os.sep))) or "."
    base = os.path.basename(os.path.abspath(path.rstrip(os.sep))) or "archive"
    candidate = os.path.join(parent, f"{base}_archive.tar.gz")
    n = 1
    while os.path.exists(candidate) and not force:
        candidate = os.path.join(parent, f"{base}_archive_{n}.tar.gz")
        n += 1
    return candidate, None


def build_manifest(path, plan, status, mode, checksums):
    lines = ["STB-ARCHIVE MANIFEST", "=" * 72,
             f"Created           : {datetime.now():%Y-%m-%d %H:%M:%S}",
             f"stb               : {VERSION}",
             f"Source folder     : {os.path.abspath(path)}",
             f"Mode              : {mode} -- {MODES[mode]}"]
    if status:
        info = status.get("info") or {}
        lines += ["", "THE RUN (stb-status)",
                  f"SystemLabel       : {status.get('label')}",
                  f"State             : {status.get('state')}" + (f" ({status['reason']})" if status.get("reason") else ""),
                  f"Run type          : {status.get('type')}"]
        if info.get("version"):
            lines.append(f"SIESTA            : {info['version']}" + (f", {info['nodes']:g} MPI processes"
                                                                       if info.get("nodes") else ""))
        if status.get("formula"):
            lines.append(f"Atoms             : {status['natoms']} ({status['formula']})")
        for key, name in (("xc", "XC functional"), ("spin", "Spin")):
            if info.get(key):
                lines.append(f"{name:<18}: {info[key]}")
        if info.get("kpoints"):
            lines.append(f"k-points          : {info['kpoints']:g}")
        if info.get("mesh"):
            lines.append(f"Mesh cutoff       : {info['mesh'][0]:g} Ry required, {info['mesh'][1]:.1f} Ry used")
        if status.get("type") == "relaxation":
            lines.append("Relaxation        : " + {True: "converged", False: "NOT converged", None: "unknown"}[
                info.get("relaxed")])
        if status.get("energy") is not None:
            lines.append(f"Total energy      : {status['energy']:.6f} eV")
        if info.get("fermi") is not None:
            lines.append(f"Fermi level       : {info['fermi']:.6f} eV")
        gap = status.get("gap")
        if gap:
            lines.append("Band gap (.EIG)   : " + ("metallic" if gap["metallic"] or gap["gap"] < 0.05 else
                                                 f"{gap['gap']:.4f} eV"))
        if status.get("max_force") is not None:
            lines.append(f"Max force         : {status['max_force']:.6f} eV/Ang")
        if status.get("wall") is not None:
            lines.append(f"Wall time         : {status['wall']:.1f} s")
    lines += ["", f"FILES IN THIS ARCHIVE ({len(plan['keep'])}; SHA-256 also in SHA256SUMS)"]
    for rel, kind, size, why in plan["keep"]:
        lines.append(f"  {checksums[rel]}  {size:>10d}  {rel}   [{kind}; {why}]")
    lines += ["", f"LEFT OUT ({len(plan['leave'])})"]
    for rel, kind, size, why in plan["leave"]:
        lines.append(f"  {size:>10d}  {rel}   [{kind}; {why}]")
    lines += ["", "RESTORE", "  tar -xzf <this archive>",
              f"  cd {os.path.basename(os.path.abspath(path.rstrip(os.sep)))} && sha256sum -c SHA256SUMS"]
    return "\n".join(lines) + "\n"


def write_and_verify(path, archive_path, plan, manifest, sums_text):
    """Writes the archive (one top folder named after the source folder), reads it back and compares every file's
    size and SHA-256. Returns (problems, n_members)."""
    top = os.path.basename(os.path.abspath(path.rstrip(os.sep))) or "archive"
    with tarfile.open(archive_path, "w:gz") as tar:
        for rel, _, _, _ in plan["keep"]:
            tar.add(os.path.join(path, rel), arcname=f"{top}/{rel}", recursive=False)
        for name, text in (("MANIFEST.txt", manifest), ("SHA256SUMS", sums_text)):
            data = text.encode()
            info = tarfile.TarInfo(name=f"{top}/{name}")
            info.size = len(data)
            info.mtime = int(datetime.now().timestamp())
            tar.addfile(info, io.BytesIO(data))
    expected = {f"{top}/{rel}": (size, None) for rel, _, size, _ in plan["keep"]}
    sums = {f"{top}/{line.split('  ', 1)[1]}": line.split("  ", 1)[0] for line in sums_text.splitlines() if line}
    problems, n = [], 0
    try:
        with tarfile.open(archive_path, "r:gz") as tar:
            for m in tar:
                n += 1
                if m.name in expected:
                    fh = tar.extractfile(m)
                    h = hashlib.sha256()
                    for chunk in iter(lambda: fh.read(1 << 20), b""):
                        h.update(chunk)
                    if m.size != expected[m.name][0]:
                        problems.append(f"{m.name}: size {m.size} instead of {expected[m.name][0]}")
                    if h.hexdigest() != sums.get(m.name):
                        problems.append(f"{m.name}: checksum differs")
                    expected.pop(m.name)
    except (OSError, tarfile.TarError) as e:
        return [f"the archive cannot be read back ({e})"], n
    problems += [f"{name} is missing from the archive" for name in expected]
    return problems, n


def by_kind(rows):
    groups = {}
    for rel, kind, size, why in rows:
        g = groups.setdefault(kind, [0, 0, [], set()])
        g[0] += 1
        g[1] += size
        g[2].append(os.path.basename(rel))
        g[3].add(why)
    return groups


# --------------------------------------------------------------------------------------------------- one folder --
def archive_one(path, args, sec, out, warn, buf):
    plan = plan_archive(path, args.mode, args.recursive, args.include or [], args.exclude or [])
    try:
        from stb.status import inspect_folder
        status = inspect_folder(path, args.label)
    except Exception:
        status = None
    result = dict(path=path, state=(status or {}).get("state", "-"), files=len(plan["keep"]),
                  size=sum(r[2] for r in plan["keep"]), archive=None, ok=False)

    sec(f"FOLDER {path}")
    total = sum(r[2] for r in plan["keep"]) + sum(r[2] for r in plan["leave"])
    out(f"Files             : {plan['n_files']} ({human_size(total)})" + (" -- including subfolders" if args.recursive
                                                                          else ""))
    out(f"Going in          : {len(plan['keep'])} files ({human_size(result['size'])}); left out "
        f"{len(plan['leave'])} ({human_size(sum(r[2] for r in plan['leave']))})")
    if status:
        info = status.get("info") or {}
        out(f"The run           : {status.get('state')}, {status.get('type')}"
            + (f", SIESTA {info['version']}" if info.get("version") else "")
            + (f", E = {status['energy']:.4f} eV" if status.get("energy") is not None else ""))
        if status.get("state") not in ("FINISHED", "NO LOG"):
            warn(f"{path}: the run is {status.get('state')}"
                 + (f" ({status['reason']})" if status.get("reason") else "") + " -- archiving an unfinished or failed run.")
        if status.get("type") == "relaxation" and info.get("relaxed") is False:
            warn(f"{path}: the relaxation did NOT converge ('Final (unrelaxed)'); the archive keeps an unrelaxed result.")
        if status.get("state") == "NO LOG":
            out("                    no SIESTA log here: the manifest has no run summary")

    sec("GOING IN")
    groups = by_kind(plan["keep"])
    inputs = [r for r in plan["keep"] if r[1] == CAT_INPUT]
    if inputs:
        print_table(["input", "size", "why"], [([r[0], human_size(r[2]), r[3]], None) for r in inputs], buf)
    rows = []
    for kind, (n, size, names, _) in sorted(groups.items(), key=lambda kv: -kv[1][1]):
        if kind == CAT_INPUT:
            continue
        shown = ", ".join(names[:3]) + (f", +{len(names) - 3}" if len(names) > 3 else "")
        rows.append(([kind, str(n), human_size(size), shown], None))
    if rows:
        print_table(["kind", "files", "size", "e.g."], rows, buf)
    rerun_ok = bool([r for r in inputs if r[0].endswith(".fdf")]) and not plan["missing_refs"] and not plan["missing_pseudo"]
    out("Rerunnable        : " + ("yes -- every .fdf, the files they read and a pseudopotential for every species" if
                                 rerun_ok else "NO"))
    if plan["subdirs"]:
        out(f"Subfolders        : {', '.join(plan['subdirs'])} -- not included (--recursive takes them)")
    if plan["missing_refs"]:
        warn(f"{path}: an .fdf reads {', '.join(plan['missing_refs'])}, which is not in the folder -- the archive "
             "cannot rerun the calculation.")
    if plan["missing_pseudo"]:
        warn(f"{path}: no pseudopotential (.psf/.psml) for {', '.join(plan['missing_pseudo'])} -- the archive cannot "
             "rerun the calculation without it.")
    if not inputs:
        warn(f"{path}: no .fdf input in the folder.")

    sec("LEFT OUT")
    groups = by_kind(plan["leave"])
    if groups:
        rows = [([kind, str(n), human_size(size), "; ".join(sorted(whys))[:70]], None)
                for kind, (n, size, names, whys) in sorted(groups.items(), key=lambda kv: -kv[1][1])]
        print_table(["kind", "files", "size", "why"], rows, buf)
        big = sorted(plan["leave"], key=lambda r: -r[2])[:3]
        out("Largest           : " + ", ".join(f"{os.path.basename(r[0])} ({human_size(r[2])})" for r in big))
        exts = sorted({suffix(r[0]) for r in plan["leave"] if r[3] != "--exclude"})
        if exts:
            out(f"Add back with     : --include {' '.join(exts)}  (or --mode full)")
    else:
        out("Nothing -- every file goes in")

    sec("ARCHIVE")
    if args.dry_run:
        out(f"Dry run           : nothing written ({len(plan['keep'])} files, {human_size(result['size'])} would go in)")
        return result
    if not plan["keep"]:
        out("Nothing to archive")
        warn(f"{path}: no file to archive (check --mode/--include/--exclude).")
        return result
    archive_path, err = archive_name(path, args.output, args.force)
    if err:
        out(color_text(f"Not written       : {err}", 'red'))
        warn(f"{path}: {err}.")
        return result
    if os.path.basename(archive_path) != f"{os.path.basename(os.path.abspath(path.rstrip(os.sep)))}_archive.tar.gz" \
            and not args.output:
        out(f"Name              : an archive with the default name exists -- writing {os.path.basename(archive_path)}"
            " (never overwritten)")
    checksums = {rel: sha256(os.path.join(path, rel)) for rel, _, _, _ in plan["keep"]}
    sums_text = "".join(f"{checksums[rel]}  {rel}\n" for rel, _, _, _ in plan["keep"])
    manifest = build_manifest(path, plan, status, args.mode, checksums)
    problems, n_members = write_and_verify(path, archive_path, plan, manifest, sums_text)
    asize = os.path.getsize(archive_path)
    out(f"Archive           : {archive_path}")
    out(f"                    {n_members} members ({len(plan['keep'])} files + MANIFEST.txt + SHA256SUMS), "
        f"{human_size(asize)} ({100 * asize / max(1, result['size']):.0f} % of {human_size(result['size'])})")
    out(f"SHA-256 (archive) : {sha256(archive_path)}")
    out("Verification      : read back, every file's size and SHA-256 compared -- "
        + ("ok" if not problems else "FAILED: " + "; ".join(problems[:3])))
    top = os.path.basename(os.path.abspath(path.rstrip(os.sep)))
    out(f"Restore with      : tar -xzf {os.path.basename(archive_path)} && cd {top} && sha256sum -c SHA256SUMS")
    if problems:
        warn(f"{path}: the archive check failed -- " + "; ".join(problems[:3]))
    result.update(archive=archive_path, ok=not problems, asize=asize)
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Package a SIESTA calculation folder into a verified .tar.gz for sharing or for the record -- the\n"
                    "complement of stb-clean (3.2). Every file is sorted by kind; what the input needs to run again\n"
                    "(.fdf files, every file they read, pseudopotentials, the job script) always goes in. The archive\n"
                    "holds MANIFEST.txt (the run, from stb-status, and every file with its SHA-256) and SHA256SUMS\n"
                    "(sha256sum -c), and is read back and checked. An existing archive is never overwritten.\n\n"
                    "Modes:\n" + "\n".join(f"  {k:<8} {v}" for k, v in MODES.items()),
        epilog="Example usage:\n"
               "  stb-archive --path si_relax --dry-run       # what would go in, and why\n"
               "  stb-archive --path si_relax                 # results mode -> si_relax_archive.tar.gz\n"
               "  stb-archive --path si_relax --include .RHO  # also the density\n"
               "  stb-archive --path si_relax --mode inputs   # only what reruns it\n"
               "  stb-archive --path 'strain/a_*'             # one archive per folder",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--path", default=".",
                        help="Directory to archive, or a glob pattern (e.g. 'strain_*') matching several directories,\n"
                             "one archive each (default: current directory).")
    parser.add_argument("--mode", choices=list(MODES), default="results", help="What goes in (default: results).")
    parser.add_argument("--include", nargs="+", default=None, metavar="EXT",
                        help="Also include files ending with these (e.g. .RHO .HSX), whatever the mode.")
    parser.add_argument("--exclude", nargs="+", default=None, metavar="EXT",
                        help="Leave out files ending with these (inputs always go in).")
    parser.add_argument("--recursive", action="store_true",
                        help="Also the files of subfolders (hidden ones skipped), in the same archive.")
    parser.add_argument("-l", "--label", default=None,
                        help="SystemLabel override for the run summary (auto-detected from the .fdf).")
    parser.add_argument("-o", "--output", default=None,
                        help="Archive filename. Default: <folder>_archive.tar.gz next to the folder (_1, _2, ... if it\n"
                             "exists). Not with a glob.")
    parser.add_argument("--force", action="store_true", help="Overwrite an existing archive of the same name.")
    parser.add_argument("--dry-run", action="store_true", help="Report what would go in; write nothing.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also write the report to {REPORT_FILE} (current folder). Off by default.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-archive {VERSION}")
    args = parser.parse_args()

    for flag, values in (("--include", args.include), ("--exclude", args.exclude)):
        for ext in values or []:
            if not ext.startswith("."):
                parser.error(f"{flag} entries must start with '.' -- got '{ext}' (did you mean '.{ext}'?)")

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - Archive Packager",
            "Bundles a finished calculation for sharing, not deleting",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    is_glob = any(c in args.path for c in "*?[")
    matches = sorted(p for p in glob.glob(args.path) if os.path.isdir(p)) if is_glob else [args.path]
    if not matches:
        parser.error(f"No directory matched '{args.path}'.")
    if not is_glob and not os.path.isdir(args.path):
        parser.error(f"'{args.path}' is not a directory.")
    batch = len(matches) > 1
    if batch and args.output:
        parser.error(f"--output can't be used when --path matches {len(matches)} directories; each gets its own "
                     "auto-named archive.")

    buf = io.StringIO() if args.save_report else None
    sec = Sections(buf)
    warnings = []

    def out(text=""):
        print_dual(text, buf)

    def warn(message):
        warnings.append(message)

    print("\n" + color_text("ARCHIVE:", 'bold'))
    print("-" * 60)
    out(color_text("===== STB-ARCHIVE REPORT =====", 'magenta'))
    sec("RUN METADATA")
    out(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}")
    out(f"stb               : {VERSION}")
    out(f"Target            : {args.path}" + (f" ({len(matches)} folders, one archive each)" if batch else "")
        + (" (recursive)" if args.recursive else ""))
    out(f"Mode              : {args.mode} -- {MODES[args.mode]}")
    if args.include or args.exclude:
        out(f"Adjusted          : include {' '.join(args.include or []) or '-'}; exclude {' '.join(args.exclude or []) or '-'}")
    out("Always in         : .fdf files, every file they read, pseudopotentials, job scripts")
    out("Never in          : hidden files, earlier archives/backups")
    if args.dry_run:
        out("Dry run           : nothing will be written")

    results = [archive_one(d, args, sec, out, warn, buf) for d in matches]

    sec("SUMMARY & FILES")
    if batch:
        rows = [([r["path"], r["state"], str(r["files"]), human_size(r["size"]),
                  os.path.basename(r["archive"]) if r["archive"] else "-",
                  human_size(r.get("asize", 0)) if r["archive"] else "-",
                  "ok" if r["ok"] else ("dry run" if args.dry_run else "FAILED")],
                 None if r["ok"] or args.dry_run else "red") for r in results]
        print_table(["folder", "run", "files", "size", "archive", "archive size", "check"], rows, buf)
    done = [r for r in results if r["ok"]]
    out(f"Archived          : {len(done)} of {len(results)} folder(s)" + (" (dry run)" if args.dry_run else ""))
    for r in done:
        out(f"  - {r['archive']}")
    if args.save_report:
        out(f"Report            : {REPORT_FILE}")

    sec("WARNINGS")
    if warnings:
        for n_, w in enumerate(warnings, 1):
            out(color_text(f"  {n_}. {w}", 'yellow'))
    else:
        out("  none")
    if buf is not None:
        with open(REPORT_FILE, "w") as fh:
            fh.write(buf.getvalue())

    print("\n" + "-" * 60)
    print(color_text("Archived -- science preserved, not deleted.\n", 'bold'))
    if not args.dry_run and not done:
        sys.exit(1)


if __name__ == "__main__":
    main()
