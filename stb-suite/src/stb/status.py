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
import csv
import glob
import time
import argparse
from datetime import datetime

import numpy as np

from stb.core import siesta_log, structure_io
from stb.core.cli import color_text, show_intro, print_dual, print_table, Sections

# Purpose: what state is a SIESTA calculation folder in, and is its result usable? One pass over the log (and
# MESSAGES, CLOCK, the .EIG) gives: the state (finished / failed, with SIESTA's own reason / probably running /
# stopped), the run (version, processes, start, wall time), the system (atoms, cell, spin, XC, k-points, mesh), the
# progress (steps done vs requested; for a relaxation, forces and stress against their tolerances and whether SIESTA
# called it relaxed), the SCF (iterations per step, tolerance, failures), the results (energies, Fermi level, gap,
# spin moment, forces, pressure) and which suite tools can use the files present. Several folders (a glob, or
# --recursive) give one table with energies relative to the lowest one per atom.

REPORT_FILE = "stb_status_report.txt"
RECENT_SECONDS = 120            # an unfinished log written this recently probably belongs to a running job

# Output files, the suite tool that reads them, and its menu code
NEXT_TOOLS = [
    ([".bands"], "band structure", "stb-bands (2.1)"),
    ([".PDOS.xml", ".PDOS", ".DOS"], "DOS / PDOS", "stb-dos (2.2)"),
    ([".EIG"], "eigenvalues", "band gap here; stb-dos (2.2)"),
    ([".RHO", ".DRHO", ".BADER", ".VT", ".VH", ".TOCH"], "real-space grids", "stb-cube (3.3), stb-bader (2.6), "
                                                                             "stb-density (2.8)"),
    ([".HSX", ".HS"], "Hamiltonian + overlap", "stb-siesta2wtb (3.4)"),
    ([".ANI"], "MD/relaxation trajectory", "stb-ani2traj (3.5), stb-aimdAnalysis (2.18)"),
    ([".STRUCT_OUT", ".XV"], "final geometry", "stb-translate (3.1)"),
    ([".DM"], "density matrix (restart)", "DM.UseSaveDM; stb-clean (3.2) keeps or removes it"),
    ([".WFSX"], "wavefunctions", "stb-fatbands (2.10), stb-wfdensity (2.12)"),
]

_RE_STEP = re.compile(r"Begin (CG opt\. move|MD step|FIRE opt\. move|Broyden opt\. move)\s*=\s*(\d+)")
_RE_SCF_LINE = re.compile(r"^\s*scf:\s+\d+")


# ------------------------------------------------------------------------------------------------ small helpers --
def human_size(n):
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024.0


def human_time(seconds):
    if seconds is None:
        return "-"
    seconds = float(seconds)
    if seconds < 120:
        return f"{seconds:.1f} s"
    if seconds < 7200:
        return f"{seconds / 60:.1f} min"
    return f"{seconds / 3600:.2f} h"


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


def after_eq(line):
    return line.split("=", 1)[1].strip() if "=" in line else ""


def first_float(text):
    try:
        return float(text.split()[0])
    except (IndexError, ValueError):
        return None


def parse_log(path):
    """Everything this tool reports from a SIESTA log, in one pass. Every field None/[] when absent."""
    info = dict(version=None, nodes=None, start=None, end=None, spin=None, xc=None, kpoints=None, mesh=None,
                dynamics=None, variable_cell=None, max_moves=None, force_tol=None, stress_tol=None, scf_max=None,
                dm_tol=None, scf_abort=None, first_md=None, last_md=None, relaxed=None, steps=[], final={},
                warnings=[], fermi=None, spin_moment=None, dipole=None, n_electrons=None)
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as fh:
            lines = fh.read().splitlines()
    except OSError:
        return None
    current = None
    xc_next = False
    in_final = False
    for i, raw in enumerate(lines):
        s = raw.strip()
        if xc_next:
            info["xc"] = s.replace("xc_check:", "").strip()
            xc_next = False
        if s.startswith("Version") and ":" in s and info["version"] is None:
            info["version"] = s.split(":", 1)[1].strip()
        elif s.startswith("* Running on") and "nodes" in s:
            info["nodes"] = first_float(s.split("on", 1)[1])
        elif s.startswith(">> Start of run:"):
            info["start"] = s.split(":", 1)[1].strip()
        elif s.startswith(">> End of run:"):
            info["end"] = s.split(":", 1)[1].strip()
        elif s.startswith("xc_check: Exchange-correlation functional"):
            xc_next = True
        elif s.startswith("redata: Spin configuration"):
            info["spin"] = after_eq(s)
        elif s.startswith("redata: Max. number of SCF Iter"):
            info["scf_max"] = first_float(after_eq(s))
        elif s.startswith("redata: DM tolerance for SCF"):
            info["dm_tol"] = first_float(after_eq(s))
        elif s.startswith("redata: SCF convergence failure will abort"):
            info["scf_abort"] = True
        elif s.startswith("redata: Dynamics option"):
            info["dynamics"] = after_eq(s)
        elif s.startswith("redata: Variable cell"):
            info["variable_cell"] = after_eq(s).upper().startswith("T")
        elif s.startswith("redata: Maximum number of optimization moves"):
            info["max_moves"] = first_float(after_eq(s))
        elif s.startswith("redata: Force tolerance"):
            info["force_tol"] = first_float(after_eq(s))
        elif s.startswith("redata: Stress tolerance"):
            info["stress_tol"] = first_float(after_eq(s))
        elif s.startswith("redata: Initial MD time step"):
            info["first_md"] = first_float(after_eq(s))
        elif s.startswith("redata:") and "Final MD time step" in s:
            info["last_md"] = first_float(after_eq(s))
        elif s.startswith("siesta: k-grid: Number of k-points"):
            info["kpoints"] = first_float(after_eq(s))
        elif s.startswith("InitMesh: Mesh cutoff (required, used)"):
            parts = after_eq(s).split()
            if len(parts) >= 2:
                info["mesh"] = (float(parts[0]), float(parts[1]))
        elif s.startswith("Total number of electrons:"):
            info["n_electrons"] = first_float(s.split(":", 1)[1])
        elif _RE_STEP.search(s):
            m = _RE_STEP.search(s)
            current = dict(step=int(m.group(2)), scf=0, converged=False, E_KS=None, max_force=None,
                           net_force=None, pressure=None, max_stress=None)
            info["steps"].append(current)
        elif _RE_SCF_LINE.match(raw):
            if current is None:              # a single point: one implicit step
                current = dict(step=1, scf=0, converged=False, E_KS=None, max_force=None, net_force=None,
                               pressure=None, max_stress=None)
                info["steps"].append(current)
            current["scf"] += 1
        elif s.startswith("SCF cycle converged after") and current is not None:
            current["converged"] = True
            n = first_float(s.split("after", 1)[1])
            if n is not None:
                current["scf"] = int(n)
        elif s.startswith("siesta: E_KS(eV) =") and current is not None:
            current["E_KS"] = first_float(after_eq(s))
        elif s.startswith("Tot ") and current is not None:
            vals = s.split()[1:4]
            try:
                current["net_force"] = float(np.linalg.norm([float(v) for v in vals]))
            except ValueError:
                pass
        elif s.startswith("Max ") and current is not None:
            current["max_force"] = first_float(s.split(None, 1)[1])
        elif s.startswith("Stress tensor Voigt") and current is not None:
            try:
                voigt = [float(v) for v in s.split(":", 1)[1].split()[:6]]
                current["max_stress"] = max(abs(v) for v in voigt)          # kBar
            except ValueError:
                pass
        elif s.startswith("siesta: Pressure (static):") and current is not None:
            rest = s.split("Pressure (static):", 1)[1].split()
            if rest:                              # one-line form "... -31.29 kBar"
                current["pressure"] = first_float(rest[0])
            else:                                 # the table form of the final block: the 'kBar' row's Solid column
                for nxt in lines[i + 1:i + 6]:
                    if nxt.strip().endswith("kBar"):
                        current["pressure"] = first_float(nxt.split(":", 1)[1])
                        break
        elif s.startswith("outcoor: Relaxed atomic coordinates"):
            info["relaxed"] = True
        elif s.startswith("outcoor: Final (unrelaxed)"):
            info["relaxed"] = False
        elif s.startswith("siesta: Final energy (eV):"):
            in_final = True
        elif in_final and s.startswith("siesta:") and "=" in s:
            key = s.split(":", 1)[1].split("=")[0].strip()
            info["final"][key] = first_float(after_eq(s))
        elif in_final and not s:
            in_final = False
        elif s.startswith("siesta: Electric dipole (Debye)"):
            try:
                info["dipole"] = np.array([float(v) for v in after_eq(s).split()[:3]])
            except ValueError:
                pass
        if "WARNING" in s and len(info["warnings"]) < 200:
            info["warnings"].append(s)
    info["fermi"] = info["final"].get("Fermi") or siesta_log.get_fermi_energy(path)
    info["spin_moment"] = siesta_log.get_spin_moment(path)
    info["mtime"] = os.path.getmtime(path)
    return info


def read_messages(folder):
    """SIESTA's MESSAGES file: (fatal lines, info lines)."""
    path = os.path.join(folder, "MESSAGES")
    if not os.path.isfile(path):
        return None, None
    fatal, infos = [], []
    with open(path, errors="ignore") as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("FATAL"):
                fatal.append(s.split(":", 1)[1].strip() if ":" in s else s)
            elif s:
                infos.append(s)
    return fatal, infos


def read_clock_seconds(folder):
    """Wall time from SIESTA's CLOCK file: the last 'End of run' (or the last line) in seconds."""
    path = os.path.join(folder, "CLOCK")
    if not os.path.isfile(path):
        return None, False
    last, ended = None, False
    with open(path, errors="ignore") as fh:
        for line in fh:
            parts = line.split()
            if parts and parts[-1].replace(".", "", 1).isdigit():
                last = float(parts[-1])
                if line.startswith("End of run"):
                    ended = True
    return last, ended


def find_run_folders(root):
    """Every non-hidden folder under `root` that holds a SIESTA log or input."""
    hits = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = sorted(d for d in dirnames if not d.startswith("."))
        if any(f.endswith(".out") for f in filenames) or any(f.endswith(".fdf") for f in filenames):
            hits.append(dirpath)
    return hits


def n_moves(steps):
    """Geometry moves done: CG/FIRE/Broyden count the input geometry as move 0."""
    if not steps:
        return 0
    return steps[-1]["step"] if steps[0]["step"] == 0 else len(steps)


def run_type(info):
    d = (info or {}).get("dynamics") or ""
    dl = d.lower()
    if not d or "single" in dl or ("cg" not in dl and "md" not in dl and "verlet" not in dl and "nose" not in dl
                                   and "parrinello" not in dl and "anneal" not in dl and "broyden" not in dl
                                   and "fire" not in dl):
        return "single point"
    if any(k in dl for k in ("cg", "broyden", "fire")) and "optim" in dl:
        return "relaxation"
    if "optim" in dl:
        return "relaxation"
    return "MD"


# ------------------------------------------------------------------------------------------------- one folder --
def inspect_folder(path, explicit_label):
    """Gathers everything about one folder; prints nothing. Returns a dict (see the keys used by the printers)."""
    label = structure_io.find_fdf_system_label(path, explicit_label)
    out_file = siesta_log.find_out_file(path, label)
    r = dict(path=path, label=label, out=out_file, fdf=None, info=None, state="NO LOG", reason=None, type="-",
             wall=None, warnings=[], files=[], size=0, natoms=None, formula=None, energy=None, gap=None)
    for fdf in sorted(glob.glob(os.path.join(path, "*.fdf"))):
        with open(fdf, errors="ignore") as fh:
            value = fdf_value(fh.read(), "SystemLabel")
        if value and label and value.split()[0] == label:
            r["fdf"] = fdf
            break
    if r["fdf"] is None:
        fdfs = sorted(glob.glob(os.path.join(path, "*.fdf")))
        r["fdf"] = fdfs[0] if fdfs else None

    # structure: the .XV (last geometry) or the input
    xv = os.path.join(path, f"{label}.XV") if label else None
    try:
        if xv and os.path.isfile(xv):
            import sisl
            g = sisl.get_sile(xv).read_geometry()
            symbols, cell = [a.symbol for a in g.atoms], np.array(g.cell)
        elif r["fdf"]:
            st = structure_io.read_fdf(r["fdf"])
            symbols, cell = [s_ for s_, _ in st.atoms], np.array(st.lattice)
        else:
            symbols, cell = None, None
    except Exception:
        symbols, cell = None, None
    if symbols:
        counts = {}
        for s_ in symbols:
            counts[s_] = counts.get(s_, 0) + 1
        r["natoms"] = len(symbols)
        r["formula"] = "".join(f"{k}{v if v > 1 else ''}" for k, v in counts.items())
        r["cell"] = cell

    fatal, infos = read_messages(path)
    r["wall"], clock_ended = read_clock_seconds(path)
    if out_file:
        info = parse_log(out_file)
        r["info"] = info
        r["type"] = run_type(info)
        finished = bool(info["end"]) or clock_ended
        age = time.time() - info["mtime"]
        if fatal:
            r["state"], r["reason"] = "FAILED", fatal[0]
        elif finished:
            r["state"] = "FINISHED"
        elif age < RECENT_SECONDS:
            r["state"], r["reason"] = "RUNNING?", f"log written {age:.0f} s ago, no end of run yet"
        else:
            r["state"], r["reason"] = "STOPPED", (f"no end of run, log last written {human_time(age)} ago -- "
                                                 "killed, out of time, or crashed")
        steps = info["steps"]
        if steps:
            last_e = next((s["E_KS"] for s in reversed(steps) if s["E_KS"] is not None), None)
            r["energy"] = info["final"].get("Total") if info["final"].get("Total") is not None else last_e
            r["max_force"] = next((s["max_force"] for s in reversed(steps) if s["max_force"] is not None), None)
            r["pressure"] = next((s["pressure"] for s in reversed(steps) if s["pressure"] is not None), None)
    # what makes two energies comparable: the same composition (element fractions) at the same level of theory
    fdf_text = ""
    if r["fdf"]:
        with open(r["fdf"], errors="ignore") as fh:
            fdf_text = fh.read()
    info_ = r["info"] or {}
    if r["natoms"]:
        fr = {k: round(v / r["natoms"], 4) for k, v in counts.items()}
        r["group"] = (tuple(sorted(fr.items())), info_.get("xc"), info_.get("spin"),
                      (info_.get("mesh") or (None,))[0], (fdf_value(fdf_text, "PAO.BasisSize") or "DZP").upper(),
                      (fdf_value(fdf_text, "PAO.EnergyShift") or "").lower())
    eig = os.path.join(path, f"{label}.EIG") if label else None
    if eig and os.path.isfile(eig):
        r["gap"] = siesta_log.read_eig_gap(eig) or None

    # files
    names = sorted(os.listdir(path))
    r["size"] = sum(os.path.getsize(os.path.join(path, n)) for n in names if os.path.isfile(os.path.join(path, n)))
    for exts, what, tool in NEXT_TOOLS:
        present = [n for n in names if any(n.endswith(e) for e in exts)]
        if present:
            r["files"].append((what, present, tool))

    # warnings
    info = r["info"]
    w = r["warnings"]
    if not out_file:
        w.append("no SIESTA log (<label>.out or the only *.out): not run here, or the log has another name.")
    else:
        if r["state"] == "FAILED":
            w.append(f"SIESTA stopped with an error: {r['reason']}")
        elif r["state"] == "RUNNING?":
            w.append(f"{r['reason']}: probably still running -- results are partial.")
        elif r["state"] == "STOPPED":
            w.append(f"the run did not finish: {r['reason']}.")
        unconverged = [s["step"] for s in info["steps"] if not s["converged"] and s["scf"] > 0]
        if unconverged:
            w.append(f"SCF not converged in {len(unconverged)} step(s) (step {unconverged[0]}"
                     + (f"..{unconverged[-1]}" if len(unconverged) > 1 else "") + ").")
        if r["type"] == "relaxation" and info["relaxed"] is False:
            w.append("relaxation NOT converged: SIESTA wrote 'Final (unrelaxed)' coordinates after "
                     f"{n_moves(info['steps'])} of {info['max_moves']:g} moves -- raise MD.NumCGsteps or restart from "
                     "the .XV." if info["max_moves"] else "relaxation NOT converged ('Final (unrelaxed)').")
        nf = next((s["net_force"] for s in reversed(info["steps"]) if s["net_force"] is not None), None)
        if nf is not None and nf > 0.05:
            w.append(f"net force on the cell {nf:.3f} eV/Ang (should be ~0): mesh too coarse (egg-box effect) or "
                     "an unconverged SCF.")
    return r


# --------------------------------------------------------------------------------------------------- printing --
def print_folder_report(r, out, sec, buf):
    info = r["info"] or {}
    sec("RUN")
    out(f"Folder            : {r['path']}")
    out(f"SystemLabel       : {r['label'] or 'unknown (no .fdf with SystemLabel; pass --label)'}")
    out(f"Input .fdf        : {r['fdf'] or 'none'}")
    out(f"Log               : {r['out'] or 'none'}" + (f" (SIESTA {info['version']})" if info.get("version") else ""))
    state_color = {"FINISHED": "green", "FAILED": "red", "RUNNING?": "yellow", "STOPPED": "yellow"}.get(r["state"])
    out("State             : " + color_text(r["state"], state_color) if state_color else
        f"State             : {r['state']}")
    if r["reason"]:
        out(f"                    {r['reason']}")
    if info:
        if info.get("nodes"):
            out(f"Processes         : {info['nodes']:g} MPI")
        out(f"Started           : {info.get('start') or '-'}")
        out(f"Ended             : {info.get('end') or '-'}")
    out(f"Wall time         : {human_time(r['wall'])}" + (" (CLOCK)" if r["wall"] is not None else ""))
    if not info:
        return

    sec("SYSTEM")
    if r["formula"]:
        out(f"Atoms             : {r['natoms']} ({r['formula']})")
    if r.get("cell") is not None:
        lengths = np.linalg.norm(r["cell"], axis=1)
        out(f"Cell (last)       : a = {lengths[0]:.4f}, b = {lengths[1]:.4f}, c = {lengths[2]:.4f} Ang, "
            f"V = {abs(np.linalg.det(r['cell'])):.4f} Ang^3")
    out(f"Electrons         : {info['n_electrons']:g}" if info.get("n_electrons") is not None else "Electrons         : -")
    out(f"Spin              : {info.get('spin') or '-'}")
    out(f"XC functional     : {info.get('xc') or '-'}")
    out(f"k-points          : {info['kpoints']:g}" if info.get("kpoints") else "k-points          : -")
    if info.get("mesh"):
        out(f"Mesh cutoff       : {info['mesh'][0]:g} Ry required, {info['mesh'][1]:.1f} Ry used")

    sec("PROGRESS")
    steps = info["steps"]
    out(f"Run type          : {r['type']}" + (f" ({info['dynamics']})" if info.get("dynamics") else ""))
    if r["type"] == "relaxation":
        out(f"Moves             : {n_moves(steps)} of at most {info['max_moves']:g} (plus the input geometry)"
            if info.get("max_moves") else f"Moves             : {n_moves(steps)}")
        out(f"Variable cell     : {'yes' if info.get('variable_cell') else 'no'}")
        verdict = {True: "converged ('Relaxed atomic coordinates')", False: "NOT converged ('Final (unrelaxed)')",
                   None: "not finished"}[info.get("relaxed")]
        out(f"Relaxation        : {verdict}")
        last = steps[-1] if steps else {}
        if info.get("force_tol") is not None and last.get("max_force") is not None:
            ok = last["max_force"] <= info["force_tol"]
            out(f"Max force (last)  : {last['max_force']:.4f} eV/Ang vs tolerance {info['force_tol']:g} "
                f"({'ok' if ok else 'above'})")
        if info.get("variable_cell") and info.get("stress_tol") is not None and last.get("max_stress") is not None:
            st_gpa = last["max_stress"] / 10.0
            ok = st_gpa <= info["stress_tol"]
            out(f"Max stress (last) : {st_gpa:.4f} GPa (largest component) vs tolerance {info['stress_tol']:g} GPa "
                f"({'ok' if ok else 'above'})")
    elif r["type"] == "MD":
        requested = (info["last_md"] - info["first_md"] + 1) if info.get("last_md") and info.get("first_md") else None
        out(f"MD steps          : {len(steps)}" + (f" of {requested:g}" if requested else ""))
        mde = siesta_log.get_mde_trajectory(os.path.join(r["path"], r["label"])) if r["label"] else []
        if mde:
            temps = np.array([m["T"] for m in mde])
            out(f"Temperature       : {temps[0]:.1f} -> {temps[-1]:.1f} K (mean {temps.mean():.1f} K, .MDE)")
            etot = [m["E_tot"] for m in mde if m["E_tot"] is not None]
            if len(etot) > 2:
                out(f"E_tot             : {etot[0]:.4f} -> {etot[-1]:.4f} eV ({1000 * (etot[-1] - etot[0]) / r['natoms']:+.3f} "
                    "meV/atom)" if r["natoms"] else "")
    if steps and len(steps) > 1:
        shown = steps if len(steps) <= 8 else steps[:3] + [None] + steps[-4:]
        rows = []
        for s in shown:
            if s is None:
                rows.append((["...", "", "", "", "", ""], None))
                continue
            rows.append(([str(s["step"]), str(s["scf"]) + ("" if s["converged"] else " (not conv.)"),
                          f"{s['E_KS']:.4f}" if s["E_KS"] is not None else "-",
                          f"{s['max_force']:.4f}" if s["max_force"] is not None else "-",
                          f"{s['pressure']:.2f}" if s["pressure"] is not None else "-",
                          f"{s['net_force']:.4f}" if s["net_force"] is not None else "-"],
                         None if s["converged"] else "yellow"))
        print_table(["step", "SCF its", "E_KS (eV)", "max F (eV/Ang)", "P (kBar)", "net F"], rows, buf)

    sec("SCF")
    if steps:
        its = np.array([s["scf"] for s in steps])
        conv = sum(1 for s in steps if s["converged"])
        out(f"Iterations        : {its.sum()} in {len(steps)} step(s); per step min {its.min()}, mean {its.mean():.1f}, "
            f"max {its.max()}")
        out(f"Converged steps   : {conv} of {len(steps)}")
    if info.get("dm_tol") is not None:
        out(f"Tolerance         : dDmax < {info['dm_tol']:g}" + (f", at most {info['scf_max']:g} iterations"
                                                               if info.get("scf_max") else "")
            + (" (failure aborts the run)" if info.get("scf_abort") else ""))

    sec("RESULTS")
    fin = info.get("final") or {}
    if fin:
        for key, name in (("Total", "Total energy"), ("Band Struct.", "Band structure"), ("Fermi", "Fermi level")):
            if fin.get(key) is not None:
                out(f"{name:<18}: {fin[key]:.6f} eV")
    elif r.get("energy") is not None:
        out(f"E_KS (last step)  : {r['energy']:.6f} eV (no final-energy block: the run did not finish)")
    if r["energy"] is not None and r["natoms"]:
        out(f"Energy per atom   : {r['energy'] / r['natoms']:.6f} eV")
    gap = r.get("gap")
    if gap:
        if gap["metallic"] or gap["gap"] < 0.05:
            out(f"Band gap (.EIG)   : metallic on the SCF k-mesh ({gap['n_k']} k-points)")
        else:
            out(f"Band gap (.EIG)   : {gap['gap']:.4f} eV ({'direct' if gap['direct'] else 'indirect'}; VBM "
                f"{gap['vbm']:.4f}, CBM {gap['cbm']:.4f} eV; {gap['n_k']} k-points)")
    if info.get("spin_moment") is not None:
        out(f"Spin moment       : {info['spin_moment']:.5f} mu_B")
    last = steps[-1] if steps else {}
    if last.get("max_force") is not None:
        out(f"Max force         : {last['max_force']:.6f} eV/Ang (net {last['net_force']:.6f})"
            if last.get("net_force") is not None else f"Max force         : {last['max_force']:.6f} eV/Ang")
    if last.get("pressure") is not None:
        out(f"Pressure          : {last['pressure']:.4f} kBar ({last['pressure'] / 10:.4f} GPa)")
    if info.get("dipole") is not None:
        out(f"Electric dipole   : {np.round(info['dipole'], 5).tolist()} Debye")
    if info.get("warnings"):
        uniq = list(dict.fromkeys(info["warnings"]))
        out(f"Log warnings      : {len(info['warnings'])} line(s), {len(uniq)} distinct -- e.g. {uniq[0][:90]}")

    sec("FILES & NEXT STEPS")
    out(f"Folder size       : {human_size(r['size'])} in {len(os.listdir(r['path']))} files")
    if r["files"]:
        rows = [([what, ", ".join(names[:3]) + (f" +{len(names) - 3}" if len(names) > 3 else ""), tool], None)
                for what, names, tool in r["files"]]
        print_table(["output", "files", "read by"], rows, buf)
    else:
        out("No analysis outputs (bands, DOS, grids, trajectory, ...) in this folder")


def main():
    parser = argparse.ArgumentParser(
        description="Status of SIESTA calculation folders: state (finished / failed, with SIESTA's reason / running /\n"
                    "stopped), run (version, processes, wall time), system (atoms, cell, spin, XC, k-points, mesh),\n"
                    "progress (steps vs requested; a relaxation's forces and stress against its tolerances), SCF\n"
                    "(iterations, failures), results (energies, Fermi level, gap, spin moment, forces, pressure)\n"
                    "and which suite tools can use the files present. Several folders (a glob, or --recursive)\n"
                    "give one table, with energies relative to the lowest per atom.",
        epilog="Example usage:\n"
               "  stb-status                               # this folder, the full report\n"
               "  stb-status --path relax                  # another folder\n"
               "  stb-status --path 'strain/a_*'           # several folders: one table, relative energies\n"
               "  stb-status --path project --recursive    # every run folder under project/\n"
               "  stb-status --path project --recursive --csv status.csv",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--path", default=".",
                        help="Directory to inspect, or a glob pattern (e.g. 'strain_*') matching several directories\n"
                             "(default: current directory).")
    parser.add_argument("--recursive", action="store_true",
                        help="Every non-hidden folder under --path that holds a .out or .fdf.")
    parser.add_argument("-l", "--label", default=None,
                        help="SystemLabel override, applied to every directory. Auto-detected per directory from\n"
                             "its .fdf file(s) if not given.")
    parser.add_argument("--details", action="store_true",
                        help="With several folders, also print the full report of each one.")
    parser.add_argument("--csv", default=None, help="Also write the summary table to this CSV file.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also write the report to {REPORT_FILE} (current folder). Off by default.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-status {VERSION}")
    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - Status Checker",
            "What's in this calc folder, and did it finish OK?",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    is_glob = any(c in args.path for c in "*?[")
    if args.recursive:
        if not os.path.isdir(args.path):
            parser.error(f"'{args.path}' is not a directory.")
        matches = find_run_folders(args.path)
    elif is_glob:
        matches = sorted(p for p in glob.glob(args.path) if os.path.isdir(p))
    else:
        if not os.path.isdir(args.path):
            parser.error(f"'{args.path}' is not a directory.")
        matches = [args.path]
    if not matches:
        parser.error(f"No run folder matched '{args.path}'.")
    batch = len(matches) > 1

    buf = io.StringIO() if args.save_report else None
    sec = Sections(buf)

    def out(text=""):
        print_dual(text, buf)

    print("\n" + color_text("STATUS:", 'bold'))
    print("-" * 60)
    out(color_text("===== STB-STATUS REPORT =====", 'magenta'))
    sec("RUN METADATA")
    out(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}")
    out(f"stb               : {VERSION}")
    out(f"Target            : {args.path}" + (" (recursive)" if args.recursive else "")
        + (f" -- {len(matches)} folder(s)" if batch else ""))

    results = [inspect_folder(d, args.label) for d in matches]

    if not batch:
        print_folder_report(results[0], out, sec, buf)
    else:
        sec("FOLDERS")
        # energies relative to the lowest per atom, among folders with the same composition
        best = {}
        comparable = [r for r in results if r["energy"] is not None and r["natoms"] and r.get("group")
                      and r["type"] != "MD"]
        for r in comparable:
            per_atom = r["energy"] / r["natoms"]
            best[r["group"]] = min(best.get(r["group"], per_atom), per_atom)
        rows, table = [], []
        root = os.path.commonpath([os.path.abspath(r["path"]) for r in results])
        for r in results:
            name = os.path.relpath(os.path.abspath(r["path"]), root) if len(results) > 1 else r["path"]
            info = r["info"] or {}
            steps = info.get("steps") or []
            scf_ok = all(s["converged"] for s in steps) if steps else None
            de = None
            if r in comparable:
                de = 1000 * (r["energy"] / r["natoms"] - best[r["group"]])
            gap = r.get("gap")
            gap_s = "-" if not gap else ("metal" if gap["metallic"] or gap["gap"] < 0.05 else f"{gap['gap']:.3f}")
            relax = ""
            if r["type"] == "relaxation":
                relax = {True: " (relaxed)", False: " (not relaxed)", None: ""}[info.get("relaxed")]
            nsteps = n_moves(steps) if r["type"] == "relaxation" else len(steps)
            cells = [name, r["state"], r["type"] + relax, str(nsteps) if steps else "-",
                     "-" if scf_ok is None else ("ok" if scf_ok else "NO"),
                     f"{r['energy']:.4f}" if r["energy"] is not None else "-",
                     f"{de:.1f}" if de is not None else "-",
                     f"{r['max_force']:.3f}" if r.get("max_force") is not None else "-",
                     f"{r['pressure']:.1f}" if r.get("pressure") is not None else "-", gap_s, human_time(r["wall"])]
            color = {"FAILED": "red", "STOPPED": "yellow", "RUNNING?": "yellow"}.get(r["state"])
            if r["state"] == "FINISHED" and (scf_ok is False or (r["type"] == "relaxation" and info.get("relaxed") is False)):
                color = "yellow"
            rows.append((cells, color))
            table.append(cells)
        headers = ["folder", "state", "type", "steps", "SCF", "E (eV)", "dE (meV/atom)", "max F", "P (kBar)",
                   "gap (eV)", "wall"]
        print_table(headers, rows, buf)
        out("dE: energy per atom above the lowest folder with the same composition and level of theory (XC, spin, "
            "mesh cutoff, basis); not for MD")
        out("steps: SCF steps, relaxation moves or MD steps; max F in eV/Ang; P = last pressure")
        counts = {}
        for r in results:
            counts[r["state"]] = counts.get(r["state"], 0) + 1
        out("States            : " + ", ".join(f"{v} {k}" for k, v in sorted(counts.items())))
        if args.csv:
            with open(args.csv, "w", newline="") as fh:
                csv.writer(fh).writerows([headers] + table)
            out(f"CSV               : {args.csv}")
        if args.details:
            for r in results:
                out("")
                out(color_text(f"===== {r['path']} =====", 'magenta'))
                print_folder_report(r, out, sec, buf)

    sec("SUMMARY & FILES")
    if not batch:
        r = results[0]
        verdict = r["state"]
        if r["state"] == "FINISHED":
            problems = [w for w in r["warnings"] if "NOT converged" in w or "SCF not converged" in w]
            verdict += " -- " + ("but see the warnings" if problems else "result usable")
        out(f"Verdict           : {verdict}")
    out(f"Folders inspected : {len(results)}")
    if args.save_report:
        out(f"Report            : {REPORT_FILE}")

    sec("WARNINGS")
    all_w = [(r["path"], w) for r in results for w in r["warnings"]]
    if all_w:
        for n_, (p, w) in enumerate(all_w, 1):
            out(color_text(f"  {n_}. " + (f"{p}: " if batch else "") + w, 'yellow'))
    else:
        out("  none")
    if buf is not None:
        with open(REPORT_FILE, "w") as fh:
            fh.write(buf.getvalue())

    print("\n" + "-" * 60)
    print(color_text("Status check complete.\n", 'bold'))


if __name__ == "__main__":
    main()
