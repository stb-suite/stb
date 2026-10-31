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
import argparse
from datetime import datetime
import numpy as np

from stb.core.deps import require_sisl
sisl = require_sisl()

from ase import Atoms
from ase.io import write as ase_write, read as ase_read

from stb.core import siesta_log
from stb.core.md_traj import read_static_lattice, unwrap_trajectory
from stb.core.kspace import detect_vacuum_axes, dimensionality_label
from stb.core.cli import color_text, show_intro, print_dual, print_table, Sections

# Purpose: SIESTA's <label>.ANI is a bare multi-frame XYZ -- no cell, no time, no energies. This tool writes the same
# frames as a trajectory that carries the lattice (per frame, from the log, for a variable-cell run), the real time and,
# for extended XYZ, each frame's thermodynamics -- and its report says what the run was (ensemble, time step, how long),
# which MD step each frame is, how the cell, temperature and energy behaved, how the atoms moved (largest step,
# boundary crossings, closest approach), and reads the written file back.
#
# Which MD step is frame k? Checked on real SIESTA 5.4 runs (examples/3.5-stb-ani2traj, WriteCoorStep true): the
# test is done again on every run -- against the per-step coordinates when the log has them (WriteCoorStep), else
# against the initial geometry -- and the thermodynamics of each frame are taken from ITS step, not from frame k's
# index.

REPORT_FILE = "stb_ani2traj_report.txt"

# --out-format -> (ASE format name, default file extension). All three carry a full per-frame lattice.
#   xsf: multi-frame AXSF (ANIMSTEPS + PRIMVEC/PRIMCOORD per frame) -- OVITO, VMD's xsf plugin, XCrySDen.
#   pdb: multi-model PDB (CRYST1 + MODEL/ENDMDL) -- VMD's own format, also OVITO. Coordinates with 3 decimals.
#   xyz: extended XYZ (Lattice="..." per frame) -- OVITO-native; VMD's xyz reader ignores the Lattice tag. The only
#        one that also carries per-frame scalars (energies, temperature, time, ...).
OUTPUT_FORMATS = {
    "xsf": ("xsf", ".xsf"),
    "pdb": ("proteindatabank", ".pdb"),
    "xyz": ("extxyz", ".xyz"),
}
# Position precision of each writer, for the read-back check (Ang)
FORMAT_TOL = {"xsf": 1e-5, "pdb": 2e-3, "xyz": 1e-6}

ENSEMBLES = {"verlet": "NVE (velocity Verlet)", "nose": "NVT (Nose thermostat)",
             "parrinellorahman": "NPH (Parrinello-Rahman, no thermostat)",
             "noseparrinellorahman": "NPT (Nose + Parrinello-Rahman)", "anneal": "annealing (velocity/cell rescaling)"}
AMU_G = 1.66053906660e-24           # g
_OUTCOOR_RE = re.compile(r"outcoor: Atomic coordinates \((\w+)\):\s*\n((?:\s*[-\d.Ee+]+\s+[-\d.Ee+]+\s+[-\d.Ee+]+\s+\d+.*\n)+)")
_MD_STEP_RE = re.compile(r"Begin MD step\s*=\s*(\d+)")


# ------------------------------------------------------------------------------------------------ small helpers --
def human_size(n):
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024.0


def fdf_value(text, key):
    """First occurrence of an fdf label (SIESTA's first-occurrence-wins), case/./_/- insensitive."""
    norm = key.lower().replace(".", "").replace("_", "").replace("-", "")
    for line in text.splitlines():
        line = line.split("#", 1)[0].strip()
        if not line or line.startswith("%"):
            continue
        parts = line.split(None, 1)
        if parts[0].lower().replace(".", "").replace("_", "").replace("-", "") == norm:
            return parts[1].strip() if len(parts) > 1 else ""
    return None


def find_input_fdf(folder, label):
    """The input .fdf of this run: <label>.fdf, else the only .fdf in the folder whose SystemLabel is `label`."""
    direct = os.path.join(folder, f"{label}.fdf")
    if os.path.isfile(direct):
        return direct
    hits = []
    for path in sorted(glob.glob(os.path.join(folder, "*.fdf"))):
        try:
            with open(path, errors="ignore") as fh:
                value = fdf_value(fh.read(), "SystemLabel")
        except OSError:
            continue
        if value and value.split()[0] == label:
            hits.append(path)
    return hits[0] if len(hits) == 1 else None


def read_log(path):
    """MD facts of a SIESTA log: version, finished, the 'redata:' MD settings, per-step data and coordinates."""
    info = dict(version=None, finished=False, dynamics=None, dt_fs=None, first_step=None, last_step=None,
                t_initial=None, t_target=None, steps=[], coords=[], coord_unit=None)
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as fh:
            text = fh.read()
    except OSError:
        return None
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("Version") and ":" in s and info["version"] is None:
            info["version"] = s.split(":", 1)[1].strip()
        elif s.startswith("redata: Dynamics option"):
            info["dynamics"] = s.split("=", 1)[1].strip()
        elif s.startswith("redata: Length of MD time step"):
            try:
                value, unit = s.split("=", 1)[1].split()[:2]
                info["dt_fs"] = float(value) * {"fs": 1.0, "ps": 1e3}.get(unit.lower(), 1.0)
            except ValueError:
                pass
        elif s.startswith("redata: Initial MD time step"):
            info["first_step"] = int(s.split("=", 1)[1])
        elif s.startswith("redata:   Final MD time step") or s.startswith("redata: Final MD time step"):
            info["last_step"] = int(s.split("=", 1)[1])
        elif s.startswith("redata: Initial Temperature of MD run"):
            info["t_initial"] = float(s.split("=", 1)[1].split()[0])
        elif s.startswith("redata: Target Temperature"):
            info["t_target"] = float(s.split("=", 1)[1].split()[0])
        elif "End of run" in s or s.startswith("Job completed"):
            info["finished"] = True
    info["steps"] = siesta_log.get_md_trajectory(path)
    # per-step coordinates (WriteCoorStep true): one 'outcoor' block after each 'Begin MD step'
    markers = list(_MD_STEP_RE.finditer(text))
    for i, m in enumerate(markers):
        end = markers[i + 1].start() if i + 1 < len(markers) else len(text)
        block = _OUTCOOR_RE.search(text, m.end(), end)
        if block:
            info["coord_unit"] = block.group(1).lower()
            rows = [ln.split()[:3] for ln in block.group(2).strip().splitlines()]
            info["coords"].append((int(m.group(1)), np.array(rows, dtype=float)))
    return info


def mic_max(a, b, cell):
    d = (a - b) @ np.linalg.inv(cell)
    d -= np.round(d)
    return float(np.abs(d @ cell).max())


def stats_row(name, values, unit, times_ps=None, digits=4):
    v = np.asarray([x for x in values if x is not None], dtype=float)
    if v.size == 0:
        return None
    drift = ""
    if times_ps is not None and v.size >= 3 and np.ptp(times_ps) > 0:
        drift = f"{np.polyfit(times_ps[:v.size], v, 1)[0]:+.{digits}g}"
    f = f"{{:.{digits}f}}"
    return ([name, unit, f.format(v[0]), f.format(v[-1]), f.format(v.mean()), f.format(v.std()),
             f.format(v.min()), f.format(v.max()), drift], None)


# --------------------------------------------------------------------------------------------------------- main --
def main():
    parser = argparse.ArgumentParser(
        description="Convert a SIESTA AIMD trajectory (<label>.ANI) into a multi-frame trajectory that also carries\n"
                    "the lattice -- for OVITO, VMD, XCrySDen. .ANI is a bare multi-frame XYZ with no cell; the cell\n"
                    "is read per frame from the SIESTA log (correct for variable-cell/NPT runs), else a single one\n"
                    "from <label>.XV or the input .fdf. The report says what the run was (ensemble, time step, how\n"
                    "long), which MD step each frame is, how the cell, temperature and energy behaved (<label>.MDE),\n"
                    "how the atoms moved (largest step, boundary crossings, closest approach), and reads the written\n"
                    "trajectory back.",
        epilog="Example usage:\n"
               "  stb-ani2traj -l water                        # water.ANI -> water_traj.xsf\n"
               "  stb-ani2traj -l si -of xyz --unwrap          # extended XYZ with E_KS, T, E_tot, P, V, time\n"
               "  stb-ani2traj -l run/si -of pdb --skip 100 --stride 5\n"
               "  stb-ani2traj --dry-run                       # the only *.ANI here; report only",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("-l", "--label", default=None,
                        help="SystemLabel of the run; may include a folder (run/si). Default: the only *.ANI here.")
    parser.add_argument(
        "-of", "--out-format", choices=sorted(OUTPUT_FORMATS), default="xsf",
        help="Output trajectory format (default: xsf):\n"
             "  xsf - multi-frame AXSF, read natively by OVITO, VMD's xsf plugin and XCrySDen.\n"
             "  pdb - multi-model PDB (CRYST1), VMD's own default format, also read by OVITO (3 decimals).\n"
             "  xyz - extended XYZ (Lattice tag), OVITO-native; VMD's xyz reader ignores\n"
             "        the Lattice tag (positions show, but without periodic boundaries).\n"
             "        The only format that also carries per-frame E_KS/Temp_ion/E_tot/P/V/Time.")
    parser.add_argument("-o", "--output", default=None,
                        help="Output filename. Default: <label>_traj.<ext>, next to the .ANI.")
    parser.add_argument("--skip", type=int, default=0,
                        help="Drop the first N frames (equilibration) before --stride (default: 0).")
    parser.add_argument("--stride", type=int, default=1,
                        help="Keep only every Nth frame (default: 1, every frame).")
    parser.add_argument(
        "--unwrap", action="store_true",
        help="Remove periodic-boundary 'jumps': each frame's positions are corrected by the minimum-image\n"
             "displacement relative to the previous frame, so atoms that cross the box edge move continuously\n"
             "instead of teleporting. Off by default (positions as SIESTA wrote them).")
    parser.add_argument("--fdf", default=None,
                        help="Input .fdf of the run (MD settings, fallback cell). Default: <label>.fdf, else the\n"
                             "only .fdf here whose SystemLabel is <label>.")
    parser.add_argument("--log", default=None,
                        help="SIESTA output log. Default: <label>.out, else the only *.out next to the .ANI.")
    parser.add_argument("--dry-run", action="store_true", help="Report only; write no trajectory.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also write the report to {REPORT_FILE} (next to the .ANI). Off by default.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-ani2traj {VERSION}")
    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - AIMD Trajectory Converter",
            "Adds back the lattice SIESTA's own .ANI file doesn't carry",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    label = args.label
    if label is None:
        found = sorted(glob.glob("*.ANI"))
        if len(found) != 1:
            print(color_text(f"[FAIL] --label not given and {len(found)} *.ANI file(s) here"
                             + (f" ({', '.join(found)})" if found else "") + " -- name the run with -l.", 'red'))
            sys.exit(1)
        label = found[0][:-4]
    folder = os.path.dirname(label) or "."
    base = os.path.basename(label)
    ase_format, default_ext = OUTPUT_FORMATS[args.out_format]
    ani_file = f"{label}.ANI"
    output_file = args.output if args.output else f"{label}_traj{default_ext}"

    report_path = os.path.join(folder, REPORT_FILE) if args.save_report else None
    buf = io.StringIO() if report_path else None
    sec = Sections(buf)
    warnings = []

    def out(text=""):
        print_dual(text, buf)

    def warn(message):
        warnings.append(message)

    def save_report():
        if buf is not None:
            with open(report_path, "w") as fh:
                fh.write(buf.getvalue())

    def fail(message):
        out(color_text(f"[FAIL] {message}", 'red'))
        save_report()
        sys.exit(1)

    print("\n" + color_text("AIMD TRAJECTORY CONVERTER:", 'bold'))
    print("-" * 60)
    out(color_text("===== STB-ANI2TRAJ REPORT =====", 'magenta'))

    # [0] --------------------------------------------------------------------------------------------------------
    sec("RUN METADATA")
    out(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}")
    out(f"stb / sisl        : {VERSION} / {sisl.__version__}")
    out(f"Input             : {ani_file}" + ("" if args.label else "  (the only *.ANI here)"))
    out(f"Output            : {'none (--dry-run)' if args.dry_run else output_file} ({args.out_format})")
    out(f"Frames kept       : skip {args.skip}, then every {args.stride}")
    out(f"Unwrap PBC        : {'yes' if args.unwrap else 'no'}")

    if args.stride < 1:
        fail(f"--stride must be >= 1, got {args.stride}.")
    if args.skip < 0:
        fail(f"--skip must be >= 0, got {args.skip}.")

    # [1] --------------------------------------------------------------------------------------------------------
    sec("INPUTS")
    if not os.path.exists(ani_file):
        fail(f"Input file '{ani_file}' not found: SIESTA writes it only with WriteMDXmol true (and <label>.MDE "
             "with WriteMDHistory true).")
    try:
        all_frames = sisl.get_sile(ani_file).read_geometry[:]()
    except Exception as e:
        fail(f"Could not read '{ani_file}': {e}")
    if not all_frames:
        fail(f"'{ani_file}' contains no frames.")
    st = os.stat(ani_file)
    nat = len(all_frames[0].atoms)
    out(f"Trajectory        : {ani_file} ({human_size(st.st_size)}, {len(all_frames)} frames x {nat} atoms)")

    fdf_path = args.fdf or find_input_fdf(folder, base)
    fdf_text = ""
    if fdf_path and os.path.isfile(fdf_path):
        with open(fdf_path, errors="ignore") as fh:
            fdf_text = fh.read()
        out(f"Input .fdf        : {fdf_path}")
    else:
        if args.fdf:
            fail(f"--fdf '{args.fdf}' not found.")
        out(f"Input .fdf        : none found ({base}.fdf, or a .fdf with SystemLabel {base}) -- pass --fdf")
    log_path = args.log or siesta_log.find_out_file(folder, base)
    if args.log and not os.path.isfile(args.log):
        fail(f"--log '{args.log}' not found.")
    log = read_log(log_path) if log_path and os.path.isfile(log_path) else None
    if log is None:
        out("SIESTA log        : none found (--log, <label>.out or the only *.out) -- no per-frame cell, no MD settings")
    else:
        out(f"SIESTA log        : {log_path}" + (f" (SIESTA {log['version']})" if log["version"] else "")
            + f", {len(log['steps'])} MD steps, " + ("finished" if log["finished"] else "NOT finished"))
        if not log["finished"]:
            warn(f"{log_path} did not finish: the run may still be going, or it stopped -- the last frames may be "
                 "missing from the log.")
    mde = siesta_log.get_mde_trajectory(label)
    out(f"Thermodynamics    : {base}.MDE " + (f"({len(mde)} rows: T, E_KS, E_tot, V, P per step)" if mde else
                                              "not found (WriteMDHistory true writes it) -- T and E_KS from the log"))
    xv_file = f"{label}.XV"
    xv_geom = None
    if os.path.isfile(xv_file):
        try:
            xv_geom = sisl.get_sile(xv_file).read_geometry()
            out(f"Last geometry     : {xv_file}")
        except Exception:
            pass

    # MD settings: the log's redata block, else the fdf
    dynamics = (log or {}).get("dynamics") or fdf_value(fdf_text, "MD.TypeOfRun")
    dt_fs = (log or {}).get("dt_fs")
    if dt_fs is None and fdf_text:
        raw = fdf_value(fdf_text, "MD.LengthTimeStep")
        if raw:
            try:
                parts = raw.split()
                dt_fs = float(parts[0]) * {"fs": 1.0, "ps": 1e3}.get(parts[1].lower() if len(parts) > 1 else "fs", 1.0)
            except ValueError:
                dt_fs = None
    first_step = (log or {}).get("first_step") or int((fdf_value(fdf_text, "MD.InitialTimeStep") or "1").split()[0])
    run_type = (fdf_value(fdf_text, "MD.TypeOfRun") or "").lower().replace(" ", "")
    if not run_type and dynamics:
        d = dynamics.lower()
        run_type = ("noseparrinellorahman" if "nose" in d and "parrinello" in d else "parrinellorahman"
                    if "parrinello" in d else "nose" if "nose" in d else "verlet" if "verlet" in d else
                    "anneal" if "anneal" in d else "")
    ensemble = ENSEMBLES.get(run_type, dynamics or "unknown")
    out(f"Dynamics          : {dynamics or run_type or 'unknown'} -> {ensemble}")
    out(f"Time step         : {f'{dt_fs:g} fs' if dt_fs else 'unknown (no log line, no MD.LengthTimeStep)'}")
    t_target = (log or {}).get("t_target")
    t_initial = (log or {}).get("t_initial")
    if t_initial is not None or t_target is not None:
        out(f"Temperatures      : initial {t_initial if t_initial is not None else '-'} K"
            + (f", target {t_target:g} K" if t_target is not None and run_type in ("nose", "noseparrinellorahman")
               else ""))

    # [2] --------------------------------------------------------------------------------------------------------
    sec("SYSTEM")
    symbols = [a.symbol for a in all_frames[0].atoms]
    counts = {}
    for s_ in symbols:
        counts[s_] = counts.get(s_, 0) + 1
    formula = "".join(f"{k}{v if v > 1 else ''}" for k, v in counts.items())
    out(f"Atoms             : {nat} ({formula})")

    # ---- cells per frame, and which MD step each frame is ----------------------------------------------------
    steps = (log or {}).get("steps") or []
    step_numbers = [s["step"] for s in steps]
    offset, how = None, None
    coords = (log or {}).get("coords") or []
    cell0 = None
    if steps and steps[0]["cell"] is not None:
        cell0 = steps[0]["cell"]
    else:
        cell0 = read_static_lattice(label, None, fdf_path=fdf_path)
        if cell0 is None and xv_geom is not None:
            cell0 = np.array(xv_geom.cell)
    if cell0 is None:
        fail(f"No lattice: no usable '{log_path or base + '.out'}', '{xv_file}' or input .fdf -- '{ani_file}' carries "
             "none of its own.")
    if coords:
        to_ang = {"ang": 1.0, "bohr": 0.529177210903}.get(log["coord_unit"], None)
        by_step = {s: c for s, c in coords}
        cell_of = {st_["step"]: st_["cell"] for st_ in steps if st_["cell"] is not None}
        scores = {}
        for off in (0, 1, -1):
            devs = []
            for k in range(min(len(all_frames), 20)):
                s = first_step + k + off
                if s in by_step:
                    c = by_step[s]
                    cell_s = cell_of.get(s, cell0)          # a variable cell: each step's own
                    cart = c @ cell_s if to_ang is None else c * to_ang
                    if cart.shape == all_frames[k].xyz.shape:
                        devs.append(mic_max(all_frames[k].xyz, cart, cell_s))
            if devs:
                scores[off] = max(devs)
        if scores:
            best = min(scores, key=scores.get)
            if scores[best] < 1e-3:
                offset, how = best, (f"checked against the log's per-step coordinates (WriteCoorStep), largest "
                                     f"difference {scores[best]:.1e} Ang")
    if offset is None:
        initial = None
        try:
            if fdf_path:
                initial = sisl.get_sile(fdf_path).read_geometry().xyz
        except Exception:
            initial = None
        if initial is not None and initial.shape == all_frames[0].xyz.shape:
            dev0 = mic_max(all_frames[0].xyz, initial, cell0)
            offset = 0 if dev0 < 1e-4 else 1
            how = (f"inferred from the first frame vs the input geometry ({dev0:.1e} Ang apart)")
        else:
            offset = 0
            how = "assumed (no per-step coordinates, no input geometry to compare)"
            warn("which MD step each frame is could not be checked (no WriteCoorStep coordinates, no input .fdf): "
                 "the per-frame energies/temperatures assume frame k = step k.")

    frame_steps = [first_step + k + offset for k in range(len(all_frames))]
    step_index = {s: i for i, s in enumerate(step_numbers)}
    per_frame_cells, cell_source = [], None
    if steps and all(step_index.get(s) is not None and steps[step_index[s]]["cell"] is not None
                     for s in frame_steps[:-1]):
        for s in frame_steps:
            i = step_index.get(s)
            if i is not None and steps[i]["cell"] is not None:
                per_frame_cells.append(steps[i]["cell"])
            elif xv_geom is not None:                       # the last geometry: its cell is the .XV's
                per_frame_cells.append(np.array(xv_geom.cell))
            else:
                per_frame_cells.append(per_frame_cells[-1])
        cell_source = f"per frame, from the log's 'outcell' block of each frame's MD step"
    else:
        per_frame_cells = [np.array(cell0)] * len(all_frames)
        cell_source = "one static cell (log has no usable per-step cell) -- " + (
            f"{xv_file}" if os.path.isfile(xv_file) else str(fdf_path))
        if run_type in ("parrinellorahman", "noseparrinellorahman"):
            warn("a variable-cell run but no per-step cell in the log: every frame gets the same cell, which is wrong "
                 "for this run.")

    lengths = np.linalg.norm(per_frame_cells[0], axis=1)
    volume = abs(np.linalg.det(per_frame_cells[0]))
    mass = sum(a.mass for a in all_frames[0].atoms)
    out(f"Cell (frame 1)    : a = {lengths[0]:.4f}, b = {lengths[1]:.4f}, c = {lengths[2]:.4f} Ang, "
        f"V = {volume:.3f} Ang^3")
    out(f"Density           : {mass * AMU_G / (volume * 1e-24):.4f} g/cm^3 ({nat / volume:.4f} atoms/Ang^3)")
    vac = detect_vacuum_axes(all_frames[0].xyz @ np.linalg.inv(per_frame_cells[0]), per_frame_cells[0], 4.0)
    out(f"Dimensionality    : {dimensionality_label(vac)}")

    # [3] --------------------------------------------------------------------------------------------------------
    sec("FRAMES AND MD STEPS")
    out(f"Frames in .ANI    : {len(all_frames)}")
    if steps:
        out(f"MD steps in log   : {len(steps)} (steps {step_numbers[0]}-{step_numbers[-1]})")
    if mde:
        out(f"Rows in .MDE      : {len(mde)} (steps {mde[0]['step']}-{mde[-1]['step']})")
    rule = ("frame k is the geometry of MD step k" if offset == 0 else
            f"frame k is the geometry AFTER MD step k (= the geometry of step k+{offset})" if offset > 0 else
            f"frame k is the geometry of MD step k{offset}")
    out(f"Frame <-> step    : {rule} -- {how}")
    out("                    each frame gets the cell, energies and temperature of ITS step")
    if steps and len(all_frames) > len(steps) + 1:
        warn(f"{ani_file} has {len(all_frames)} frames but the log only {len(steps)} MD steps: SIESTA APPENDS to an "
             "existing .ANI, so frames of an earlier run are probably at the start of the file.")
    if dt_fs:
        t0 = (frame_steps[0] - first_step) * dt_fs
        t1 = (frame_steps[-1] - first_step) * dt_fs
        out(f"Time span         : {t0:g} to {t1:g} fs ({(t1 - t0) / 1000:.4f} ps; t = 0 at the input geometry)")
    if xv_geom is not None and xv_geom.xyz.shape == all_frames[-1].xyz.shape:
        dev = mic_max(all_frames[-1].xyz, xv_geom.xyz, per_frame_cells[-1])
        out(f"Last frame vs .XV : {dev:.1e} Ang apart" + (" (same geometry)" if dev < 1e-3 else ""))
        if dev > 1e-3:
            warn(f"the last frame of {ani_file} is {dev:.3f} Ang from {xv_file}: the two may come from different "
                 "runs (or the run was restarted).")

    # frames kept
    keep = list(range(len(all_frames)))[args.skip::args.stride]
    if not keep:
        fail(f"--skip {args.skip} leaves no frame (the .ANI has {len(all_frames)}).")
    frames = [all_frames[i] for i in keep]
    cells = [per_frame_cells[i] for i in keep]
    fsteps = [frame_steps[i] for i in keep]
    out(f"Kept              : {len(frames)} frames (indices {keep[0] + 1}..{keep[-1] + 1} of {len(all_frames)}, "
        f"every {args.stride})")

    # per-frame thermodynamics from the frame's own step
    mde_by_step = {r["step"]: r for r in mde}
    log_by_step = {s["step"]: s for s in steps}
    times = [((s - first_step) * dt_fs if dt_fs else None) for s in fsteps]

    def thermo(s, key_mde, key_log=None):
        if s in mde_by_step and mde_by_step[s].get(key_mde) is not None:
            return mde_by_step[s][key_mde]
        if key_log and s in log_by_step:
            return log_by_step[s].get(key_log)
        return None

    # [4] --------------------------------------------------------------------------------------------------------
    sec("LATTICE")
    out(f"Source            : {cell_source}")
    vols = np.array([abs(np.linalg.det(c)) for c in cells])
    lens = np.array([np.linalg.norm(c, axis=1) for c in cells])
    varies = np.ptp(vols) > 1e-6 * vols.mean() or np.ptp(lens) > 1e-6
    if varies:
        rows = [([n, f"{v[0]:.4f}", f"{v[-1]:.4f}", f"{v.mean():.4f}", f"{v.std():.4f}", f"{v.min():.4f}",
                  f"{v.max():.4f}"], None) for n, v in (("a (Ang)", lens[:, 0]), ("b (Ang)", lens[:, 1]),
                                                        ("c (Ang)", lens[:, 2]), ("V (Ang^3)", vols))]
        print_table(["", "first", "last", "mean", "std", "min", "max"], rows, buf)
        out(f"Volume change     : {100 * (vols[-1] / vols[0] - 1):+.2f} % first -> last, {100 * vols.std() / vols.mean():.2f} "
            "% fluctuation")
    else:
        out("Cell              : constant over the frames kept")

    # [5] --------------------------------------------------------------------------------------------------------
    sec("THERMODYNAMICS")
    t_ps = np.array([t / 1000 for t in times]) if dt_fs else None
    series = {
        "T": [thermo(s, "T", "Temp_ion") for s in fsteps],
        "E_KS": [thermo(s, "E_KS", "E_KS") for s in fsteps],
        "E_tot": [thermo(s, "E_tot") for s in fsteps],
        "V": [thermo(s, "volume") for s in fsteps],
        "P": [thermo(s, "pressure") for s in fsteps],
    }
    rows = []
    for name, unit, key, digits in (("T", "K", "T", 2), ("E_KS", "eV", "E_KS", 4), ("E_tot", "eV", "E_tot", 4),
                                    ("V", "Ang^3", "V", 3), ("P", "kBar", "P", 2)):
        r = stats_row(name, series[key], unit, t_ps, digits)
        if r:
            rows.append(r)
    if rows:
        out(f"Per frame kept, from {'the .MDE' if mde else 'the log'} at each frame's own MD step "
            f"(drift = slope of a linear fit, per ps):")
        print_table(["", "unit", "first", "last", "mean", "std", "min", "max", "drift /ps"], rows, buf)
        temps = np.array([x for x in series["T"] if x is not None], float)
        if temps.size > 2:
            ndof = 3 * nat - 3
            expected = np.sqrt(2.0 / ndof)
            out(f"T fluctuation     : {temps.std() / temps.mean():.3f} relative; canonical expectation sqrt(2/(3N-3)) = "
                f"{expected:.3f} for {nat} atoms")
            if t_target is not None and run_type in ("nose", "noseparrinellorahman"):
                out(f"T vs target       : mean {temps.mean():.1f} K vs {t_target:g} K ({temps.mean() - t_target:+.1f} K)")
        etot = np.array([x for x in series["E_tot"] if x is not None], float)
        if etot.size > 2 and dt_fs and run_type == "verlet":
            slope = np.polyfit(t_ps[:etot.size], etot, 1)[0]
            out(f"E_tot drift (NVE) : {1000 * slope / nat:+.3f} meV/atom/ps; fluctuation {1000 * etot.std() / nat:.3f} "
                "meV/atom")
            if abs(1000 * slope / nat) > 5:
                warn(f"the total energy of this NVE run drifts by {1000 * slope / nat:+.1f} meV/atom/ps: the time step "
                     "or the SCF tolerance is too large for good energy conservation.")
        out("                    (stb-aimdAnalysis, menu 2.18, does RDF, MSD and VDOS on the same trajectory)")
    else:
        out("No per-step thermodynamics (no .MDE, no MD steps in the log)")

    # [6] --------------------------------------------------------------------------------------------------------
    sec("MOTION")
    pos_all = [g.xyz for g in all_frames]
    jumps, max_step, max_step_at = 0, 0.0, (0, 0)
    for t in range(1, len(pos_all)):
        inv = np.linalg.inv(per_frame_cells[t])
        raw = (pos_all[t] - pos_all[t - 1]) @ inv
        jumps += int(np.count_nonzero(np.abs(raw) > 0.5))
        mic = raw - np.round(raw)
        dist = np.linalg.norm(mic @ per_frame_cells[t], axis=1)
        i = int(np.argmax(dist))
        if dist[i] > max_step:
            max_step, max_step_at = float(dist[i]), (t, i)
    outside = 0
    for t, p in enumerate(pos_all):
        f = p @ np.linalg.inv(per_frame_cells[t])
        outside += int(np.count_nonzero((f < -1e-6) | (f >= 1 + 1e-6)))
    out(f"Largest step      : {max_step:.4f} Ang between consecutive frames ({symbols[max_step_at[1]]}"
        f"{max_step_at[1] + 1}, frames {max_step_at[0]}->{max_step_at[0] + 1})"
        + (f" -> up to {max_step / dt_fs * 1e5:.0f} m/s" if dt_fs else ""))
    if max_step > 0.3:
        warn(f"an atom moved {max_step:.2f} Ang in one frame: the time step may be too large (or frames are missing).")
    out(f"Boundary jumps    : {jumps} coordinate jumps of more than half a cell between consecutive frames")
    out(f"Outside the cell  : {outside} atom-frames with a fractional coordinate outside [0, 1)")
    if outside and not jumps:
        out("                    positions are not wrapped into the cell, so atoms move continuously already: "
            "--unwrap changes nothing here")
    if jumps and not args.unwrap:
        out("                    --unwrap makes the jumps continuous motion")
    p_first = pos_all[keep[0]]
    p_last = (unwrap_trajectory([pos_all[i] for i in keep], cells)[-1] if len(keep) > 1 else p_first)
    disp = np.linalg.norm(p_last - p_first, axis=1)
    rows = []
    for el in counts:
        idx = [i for i, s_ in enumerate(symbols) if s_ == el]
        rows.append(([el, str(len(idx)), f"{np.sqrt(np.mean(disp[idx] ** 2)):.3f}", f"{disp[idx].max():.3f}"], None))
    out("Displacement, first -> last frame kept (unwrapped):")
    print_table(["species", "atoms", "RMS (Ang)", "max (Ang)"], rows, buf)
    dmin, dmin_at = np.inf, None
    for k_, (g, c) in enumerate(zip(frames, cells)):
        a = Atoms(symbols=symbols, positions=g.xyz, cell=c, pbc=True)
        if len(a) > 1:
            d = a.get_all_distances(mic=True)
            np.fill_diagonal(d, np.inf)
            i, j = np.unravel_index(int(np.argmin(d)), d.shape)
            if d[i, j] < dmin:
                dmin, dmin_at = float(d[i, j]), (k_, i, j)
    if dmin_at:
        k_, i, j = dmin_at
        out(f"Closest approach  : {dmin:.3f} Ang ({symbols[i]}{i + 1}-{symbols[j]}{j + 1}, frame {keep[k_] + 1})")
        if dmin < 0.5:
            warn(f"two atoms came within {dmin:.2f} Ang: the run may have blown up.")

    # [7] --------------------------------------------------------------------------------------------------------
    sec("OUTPUT")
    positions = [g.xyz for g in frames]
    if args.unwrap:
        positions = unwrap_trajectory(positions, cells)
    ase_frames = []
    embedded = []
    for i, g in enumerate(frames):
        atoms = Atoms(symbols=symbols, positions=positions[i], cell=cells[i], pbc=True)
        if args.out_format == "xyz":
            s = fsteps[i]
            for key, value in (("E_KS", thermo(s, "E_KS", "E_KS")), ("Temp_ion", thermo(s, "T", "Temp_ion")),
                               ("E_tot", thermo(s, "E_tot")), ("Volume", thermo(s, "volume")),
                               ("Pressure_kBar", thermo(s, "pressure"))):
                if value is not None:
                    atoms.info[key] = value
            atoms.info["MD_step"] = s
            if times[i] is not None:
                atoms.info["Time"] = times[i]
        ase_frames.append(atoms)
    if args.out_format == "xyz":
        embedded = sorted(set().union(*[set(a.info) for a in ase_frames]))
    if args.dry_run:
        out(f"Dry run           : nothing written ({len(ase_frames)} frames would go to {output_file})")
    else:
        try:
            ase_write(output_file, ase_frames, format=ase_format)
        except Exception as e:
            fail(f"Failed to write output file: {e}")
        out(f"File              : {output_file} ({human_size(os.path.getsize(output_file))}, {len(ase_frames)} frames, "
            f"{args.out_format})")
        out(f"Per frame         : cell" + (", positions unwrapped" if args.unwrap else ", positions as in the .ANI")
            + (f"; {', '.join(embedded)}" if embedded else "; no scalars (only xyz carries them)"))
        try:
            back = ase_read(output_file, index=":", format=ase_format)
            problems = []
            if len(back) != len(ase_frames):
                problems.append(f"{len(back)} frames instead of {len(ase_frames)}")
            pos_err = max(float(np.abs(b.positions - a.positions).max()) for b, a in zip(back, ase_frames))
            cell_err = max(float(np.abs(b.cell.cellpar() - a.cell.cellpar()).max()) for b, a in zip(back, ase_frames))
            if pos_err > FORMAT_TOL[args.out_format]:
                problems.append(f"positions differ by {pos_err:.1e} Ang")
            if cell_err > 1e-2:
                problems.append(f"cell parameters differ by {cell_err:.1e}")
            if any(list(b.get_chemical_symbols()) != symbols for b in back):
                problems.append("species differ")
            if args.out_format == "xyz":
                lost = sorted({k for b, a in zip(back, ase_frames) for k in a.info if k not in b.info})
                if lost:
                    problems.append(f"per-frame values lost: {', '.join(lost)}")
            out("Read-back check   : read again with ASE -- " + ("ok" if not problems else "FAILED: " + "; ".join(problems)))
            out(f"                    {len(back)} frames; positions to {pos_err:.1e} Ang, cell parameters to "
                f"{cell_err:.1e} (the {args.out_format} writer's precision)")
            if problems:
                warn("the trajectory read back differs from what was written: " + "; ".join(problems))
        except Exception as e:
            warn(f"could not read {output_file} back to check it: {e}")

    # [8] --------------------------------------------------------------------------------------------------------
    sec("SUMMARY & FILES")
    out(f"Run               : {formula}, {ensemble}, {len(all_frames)} frames"
        + (f", dt {dt_fs:g} fs" if dt_fs else ""))
    if not args.dry_run:
        out(f"Written           : {output_file} ({len(ase_frames)} frames)")
        viewer = {"xsf": f"ovito {output_file}  |  vmd {output_file}  |  xcrysden --axsf {output_file}",
                  "pdb": f"vmd {output_file}  |  ovito {output_file}",
                  "xyz": f"ovito {output_file}  (VMD reads it without the cell)"}[args.out_format]
        out(f"Open with         : {viewer}")
    if report_path:
        out(f"Report            : {report_path}")

    sec("WARNINGS")
    if warnings:
        for n_, w in enumerate(warnings, 1):
            out(color_text(f"  {n_}. {w}", 'yellow'))
    else:
        out("  none")
    save_report()

    print("\n" + "-" * 60)
    print(color_text("Trajectory converted -- now go make a nice movie of it.\n", 'bold'))


if __name__ == "__main__":
    main()
