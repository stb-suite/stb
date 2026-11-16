"""Basis-set convergence by equations of state: layout, defaults and the analysis
shared by stb-convergence --parameter basis (which writes the scans) and
stb-convergenceAnalysis (which fits them).

A basis of numerical atomic orbitals has two independent axes: range
(PAO.EnergyShift) and cardinality (PAO.BasisSize). The total energy has no
plateau in either, so each candidate setting is judged by a physical quantity:
the equilibrium volume V0 (in-plane area for a 2D material) of an E(V) scan of
single points. Three stages, each one a set of such scans:

  mesh  : Mesh.Cutoff, at a tight PAO.EnergyShift -- the real-space mesh the
          other two stages need (without an orbital filter, 350 Ry is not
          enough for every pseudopotential)
  range : PAO.EnergyShift, at the mesh found
  size  : PAO.BasisSize one step up, at the EnergyShift found

Layout: <output-dir>/basis/<stage>_<value>/vol_<pct>/ and one
basis_scan_<stage>.json next to the sub-scans.
"""

import os
import re
import json
import glob

import numpy as np

STAGES = ("mesh", "range", "size")
STAGE_TITLES = {"mesh": "real-space mesh (Mesh.Cutoff)", "range": "orbital range (PAO.EnergyShift)",
                "size": "cardinality (PAO.BasisSize)"}
MESH_VALUES = [350.0, 500.0, 700.0, 1000.0]
SHIFT_VALUES = [0.02, 0.01, 0.005, 0.002, 0.001, 0.0005]
MESH_STAGE_SHIFT = 0.002
SIZE_LADDER = ["SZ", "DZ", "DZP", "TZP", "TZ2P"]
SCAN_POINTS = 9
# Half-width of the scan, % of the scanned measure: +-12 % in volume for a solid, +-8 % in area for a
# 2D material (+-4 % in its lattice vectors), +-4 % in length for a wire.
SCAN_RANGE_PCT = {3: 12.0, 2: 8.0, 1: 4.0}
# The control quantity is V0 (which in 2D is the area times a constant, in 1D the length). Largest
# relative change, in %, between a setting and each of the next two tighter ones for it to count as
# closed. From nine systems (Si, Ge, GaAs, NaCl, Al, Cu, Ag, graphene, hBN): a starting point.
CLOSE_PCT = {"mesh": {3: 0.2, 2: 0.2, 1: 0.1}, "range": {3: 0.5, 2: 0.3, 1: 0.15}, "size": {3: 1.0, 2: 0.6, 1: 0.3}}
# B0 is quoted in the verdict only when the fit residual is below this (meV/atom).
B0_RESIDUAL_MEV = 0.2
# SIESTA's rule of thumb for the largest orbital cutoff radius (Bohr).
RC_RANGE_BOHR = (5.0, 10.0)


def manifest_name(stage):
    return f"basis_scan_{stage}.json"


def sub_scan_name(stage, value):
    if stage == "mesh":
        return f"mesh_{value:g}"
    if stage == "range":
        return f"range_{value:g}"
    return f"size_{value}"


def strain_folder(strain_pct):
    return f"vol_{'m' if strain_pct < 0 else ''}{abs(strain_pct):.2f}"


def next_size(size):
    """The basis size one step above `size` on SIESTA's ladder (None at the top or for an unknown name)."""
    upper = [s.upper() for s in SIZE_LADDER]
    if size.upper() not in upper or size.upper() == upper[-1]:
        return None
    return SIZE_LADDER[upper.index(size.upper()) + 1]


def orbital_radii(out_path):
    """(smallest, largest) orbital cutoff radius in Bohr from the 'rc =' lines of the basis generation."""
    radii = []
    try:
        with open(out_path, errors="ignore") as fh:
            for line in fh:
                match = re.match(r"^\s*rc\s*=\s*([0-9]+\.[0-9]+)\s*$", line)
                if match:
                    radii.append(float(match.group(1)))
                elif radii and line.startswith("siesta: "):
                    break
    except OSError:
        return None
    radii = [r for r in radii if r > 0.0]
    return (min(radii), max(radii)) if radii else None


def read_manifests(base_dir):
    """Every basis_scan_<stage>.json under base_dir, {stage: manifest}."""
    found = {}
    for stage in STAGES:
        path = os.path.join(base_dir, manifest_name(stage))
        if os.path.isfile(path):
            try:
                with open(path) as fh:
                    found[stage] = json.load(fh)
            except (OSError, json.JSONDecodeError):
                pass
    return found


def fit_sub_scan(folder, out_file, n_atoms):
    """Fits E(V) of one sub-scan (its vol_*/ single points). Returns a dict: n (points read), v0, b0,
    rms (meV/atom), radii, filter (True when a log shows the orbital filter on), problems (list of
    strings); v0 is None when fewer than 5 points could be read or the fit failed."""
    from stb.core import siesta_log
    from stb.core.eos_fit import fit_eos
    volumes, energies, radii, problems, filtered = [], [], None, [], False
    folders = sorted(glob.glob(os.path.join(folder, "vol_*")))
    for vol_dir in folders:
        out_path = os.path.join(vol_dir, out_file)
        if not os.path.isfile(out_path):
            continue
        energy = siesta_log.get_total_energy(out_path)
        cell = siesta_log.get_outcell(out_path)
        if energy is None or cell is None:
            problems.append(f"{os.path.basename(vol_dir)}: unfinished")
            continue
        if not siesta_log.get_scf_convergence(out_path)[0]:
            problems.append(f"{os.path.basename(vol_dir)}: SCF not confirmed")
        volumes.append(abs(float(np.linalg.det(cell))))
        energies.append(float(energy))
        if radii is None:
            radii = orbital_radii(out_path)
            try:
                with open(out_path, errors="ignore") as fh:
                    head = fh.read(400000)
                filtered = bool(re.search(r"^\s*FilterCutoff\s+[1-9]|Filter:\s+Number of eigenvalues|"
                                          r"filtering orbs|Filtering orbital", head, re.M | re.I))
            except OSError:
                pass
    result = {"n": len(volumes), "n_folders": len(folders), "v0": None, "b0": None, "rms": None,
              "radii": radii, "filter": filtered, "problems": problems, "outside": False}
    if len(volumes) < 5:
        return result
    order = np.argsort(volumes)
    volumes, energies = np.array(volumes)[order], np.array(energies)[order]
    try:
        fit = fit_eos(volumes, energies, "birchmurnaghan")
    except Exception as error:
        problems.append(f"fit failed: {error}")
        return result
    residual = fit.get("residuals")
    result.update(v0=float(fit["v0"]), b0=float(fit["b0_gpa"]), outside=bool(fit["v0_outside_range"]),
                  rms=(float(1000 * np.sqrt(np.mean(np.asarray(residual) ** 2)) / max(n_atoms, 1))
                       if residual is not None else None))
    return result


def closing_point(v0_list, tolerance_pct):
    """Index of the first setting (list ordered from loosest to tightest) whose V0 differs by less than
    tolerance_pct from the next TWO settings, and how it was verified: (index, 'two') / (index, 'one')
    when only one tighter setting exists and it agrees / (None, None) when nothing closes. Settings
    without a fit (None) cannot be closing points and cannot confirm one."""
    n = len(v0_list)

    def close(i, j):
        return (v0_list[i] is not None and v0_list[j] is not None
                and abs(v0_list[i] / v0_list[j] - 1.0) * 100.0 < tolerance_pct)

    for i in range(n):
        if i + 2 < n and close(i, i + 1) and close(i, i + 2):
            return i, "two"
    if n >= 2 and close(n - 2, n - 1):
        return n - 2, "one"
    return None, None
