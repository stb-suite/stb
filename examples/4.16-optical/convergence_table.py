#!/usr/bin/env python3
"""Side-by-side summary of several Stage 2 runs of the same crystal (helper of example_4.16.sh).
Usage: convergence_table.py <dir> <run> [<run> ...]

Each <run> is a folder <dir>/<run>/ on which stb-opticalAnalysis already ran, with its report text in
<dir>/<run>_stage2.log. Prints, per run: the setting that differs from the FIRST run (the baseline), the
static dielectric constant, the energy and height of the strongest eps2 peak, the energy of the
energy-loss (plasmon) maximum, and SIESTA's f-sum check -- the quantities whose convergence README
section 5.3 discusses.
"""
import re
import sys

root, runs = sys.argv[1], sys.argv[2:]


def grab(text, pattern, group=1, default=None, flags=re.M):
    m = re.search(pattern, text, flags)
    return m.group(group) if m else default


rows = []
for run in runs:
    t = open(f"{root}/{run}_stage2.log").read()
    eps1 = float(grab(t, r"^xx\s*\|\s*([-+0-9.]+)\s*\|", 1, "nan"))
    peaks = [(float(e), float(v)) for e, v in re.findall(r"^\d+\s*\|\s*([0-9.]+)\s*\|\s*([0-9.]+)\s*\|", t, re.M)]
    main = max(peaks, key=lambda p: p[1]) if peaks else (float("nan"), float("nan"))
    loss = float(grab(t, r"Energy-loss peak \(bulk plasmon\)\s*\|\s*([0-9.]+)", 1, "nan"))
    fsum = float(grab(t, r"^xx\s*\|\s*yes\s*\|\s*([0-9.]+)", 1, "nan"))
    setting = {
        "basis": grab(t, r"^Basis size\s*\|\s*(\S+)"),
        "shift": grab(t, r"^PAO energy shift\s*\|\s*(.+?)\s*\|?$"),
        "cutoff": grab(t, r"^Mesh cutoff\s*\|\s*(.+?)\s*$"),
        "mesh": grab(t, r"^xx\s*\|\s*dir_xx\s*\|\s*(\d+)x"),
        "broaden": grab(t, r"^xx\s*\|\s*dir_xx\s*\|\s*[^|]+\|\s*([0-9.]+)\s*\|"),
        "bands": grab(t, r"^xx\s*\|\s*dir_xx\s*\|\s*[^|]+\|\s*[0-9.]+\s*\|\s*(\d+)\s*\|"),
    }
    rows.append((run, setting, eps1, main, loss, fsum))

base = rows[0][1]
names = {"basis": "basis", "shift": "EnergyShift", "cutoff": "Mesh.Cutoff", "mesh": "Optical.Mesh",
         "broaden": "Broaden (eV)", "bands": "bands"}
print(f"{'run':<12s}{'what differs from production':<30s}{'eps1(0)':>9s}{'main eps2 peak (eV)':>22s}"
      f"{'loss peak (eV)':>16s}{'f-sum':>8s}")
for run, setting, eps1, main, loss, fsum in rows:
    diffs = [f"{names[k]} = {setting[k]}" for k in names if setting[k] != base[k]]
    what = "(baseline)" if run == runs[0] else "; ".join(diffs) or "(same)"
    print(f"{run:<12s}{what:<30s}{eps1:>9.3f}{main[0]:>12.2f} ({main[1]:>5.1f}){loss:>12.2f}{fsum:>8.3f}")
