#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

from stb import __version__ as VERSION

import os
import re
import sys
import glob
import argparse
import datetime
import numpy as np
from stb.core import structure_io, kspace, fdf_directives
from stb.core.cli import color_text, show_intro, print_dual, print_section, print_table
from stb.core.dielectric import read_epsimg
from stb.core.calc_directives import OPTICAL_DIRECTION_VECTORS, OPTICAL_OFFDIAG_PAIRS
from stb.core.optical_properties import (
    compute_all, correct_2d_perpendicular, correct_2d_parallel, molecular_polarizability,
    derive_from_eps, find_peaks_table, absorption_onset, tauc_direct_gap, zero_crossings,
    free_electron_plasmon_ev, effective_electrons, penn_gap_ev, apply_scissor,
)
from stb.core.optical_setup import reciprocal_spacings, spacing_band, symmetry_analysis
from stb.core.siesta_log import get_scf_convergence, get_optical_info, read_eig_gap
from stb.core.spectrum import read_experimental_spectrum, find_spectrum_peaks, match_peaks

REPORT_FILE = "optical_stage2.txt"
RY_TO_EV = 13.605693122994
HC_EV_NM = 1239.841984  # photon energy (eV) * wavelength (nm)

_LABEL_RE = re.compile(r'^\s*SystemLabel\s+(\S+)', re.IGNORECASE | re.MULTILINE)
_VECTOR_BLOCK_RE = re.compile(
    r'%block\s+Optical\.Vector\s*\n\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)',
    re.IGNORECASE)

# Same direction -> Optical.Vector convention stb-optical (Stage 1) writes
# folders with -- see core/calc_directives.py's own docstring for the
# full diagonal ('xx'/'yy'/'zz')-vs-biaxial ('xy'/'xz'/'yz') derivation.
_AXIS_UNIT_VECTORS = OPTICAL_DIRECTION_VECTORS
_DIAGONAL_AXES = ("xx", "yy", "zz")
_AXIS_ORDER = {"xx": 0, "yy": 1, "zz": 2, "xy": 3, "xz": 4, "yz": 5, "avg": 6}

_QUANTITY_COLUMNS = {
    # name -> [(column_key, dat_column_index_1based, legend_label), ...]
    "eps": [("eps1", 2, "eps1"), ("eps2", 3, "eps2")],
    "n_k": [("n", 4, "n"), ("k", 5, "k")],
    "alpha": [("alpha", 6, "alpha (cm^-1)")],
    "R": [("R", 7, "R")],
    "L": [("L", 8, "L")],
    "sigma1": [("sigma1", 9, "sigma1 (S/m)")],
}
_QUANTITY_TITLES = {
    "eps": "Dielectric function", "n_k": "Refractive index / extinction coefficient",
    "alpha": "Absorption coefficient", "R": "Normal-incidence reflectivity",
    "L": "Energy-loss function", "sigma1": "Real optical conductivity",
}

_PALETTE = ["#2255cc", "#cc2222", "#22aa55", "#aa22cc", "#cc8822", "#22aacc"]

CONFIG_EXTRA_FILE = "config_extra.fdf"


def _read_folder_fdf_text(folder):
    """calc.fdf's own text, plus config_extra.fdf's (if present),
    concatenated -- stb-optical (Stage 1) writes SystemLabel and the rest
    of the user's template directly into calc.fdf, but the forced
    MD.TypeOfRun/MD.Steps and Optical.Mesh/Optical.Vector/Optical.Broaden
    block into a %include'd config_extra.fdf sidecar instead (this
    suite's standard config_extra.fdf convention -- never edited into
    calc.fdf itself, see CLAUDE.md). detect_label/detect_axis below search
    this combined text so they find SystemLabel/Optical.Vector regardless
    of which of the 2 files the directive actually lives in.
    """
    parts = []
    calc_path = os.path.join(folder, "calc.fdf")
    if os.path.isfile(calc_path):
        with open(calc_path) as f:
            parts.append(f.read())
    extra_path = os.path.join(folder, CONFIG_EXTRA_FILE)
    if os.path.isfile(extra_path):
        with open(extra_path) as f:
            parts.append(f.read())
    return "\n".join(parts)


def detect_label(folder):
    """SystemLabel from one direction folder's calc.fdf/config_extra.fdf --
    same regex approach as raman_analysis.py::detect_optical_label/
    core.phonon_workflow.detect_system_label, reimplemented locally per
    -folder (not a single workflow-wide detection) since Stage 2
    aggregates folders that could in principle carry different
    SystemLabels (not the normal case, but cheap to handle correctly
    rather than assume).
    """
    match = _LABEL_RE.search(_read_folder_fdf_text(folder))
    return match.group(1) if match else None


def detect_axis(folder):
    """Direction label ('xx'/'yy'/'zz' diagonal, 'xy'/'xz'/'yz' biaxial)
    from the %block Optical.Vector in one direction folder's calc.fdf/
    config_extra.fdf -- read from the files' own content, NOT inferred
    from the folder name (robust to the user renaming/reorganizing Stage
    1's dir_xx/dir_yy/.../dir_yz folders). Matches both the +vector and
    -vector convention (e.g. -1,0,0 is still "xx") -- physically the SAME
    linear response either way (eps enters via |<i|p.e|j>|^2, quadratic
    in the polarization unit vector e, so a sign flip cannot change the
    result; stb-optical itself only ever writes the +vector convention,
    this just makes detection robust to a hand-edited or externally
    -generated folder using the opposite sign). Returns None if no block
    is found or the vector doesn't match one of the 6 canonical
    directions at all (e.g. a genuine off-axis, non-bisector direction)
    -- callers fall back to the folder's basename as a display label in
    that case.
    """
    match = _VECTOR_BLOCK_RE.search(_read_folder_fdf_text(folder))
    if match is None:
        return None
    vec = tuple(float(v) for v in match.groups())
    # Tolerance-based, not exact-tuple equality: build_optical_block writes
    # the bisector directions ((+/-1/sqrt(2), ...)) rounded to 4 decimal
    # places, an irrational value that would never exactly equal this
    # module's own unrounded constant otherwise. 1e-3 comfortably covers
    # that rounding (worst case ~5e-5) with no risk of matching a genuine
    # off-axis direction by accident.
    for axis, unit_vec in _AXIS_UNIT_VECTORS.items():
        neg_unit_vec = tuple(-v for v in unit_vec)
        if (all(abs(a - b) < 1e-3 for a, b in zip(vec, unit_vec)) or
                all(abs(a - b) < 1e-3 for a, b in zip(vec, neg_unit_vec))):
            return axis
    return None


def find_direction_folders(directory):
    """Every immediate subdirectory of `directory` that contains a
    calc.fdf -- the general "aggregate whatever's present" scan this
    Stage 2 uses instead of hardcoding dir_xx/dir_yy/.../dir_yz (Stage 1's
    own naming is just a convention, not a contract Stage 2 depends on).
    Returns a sorted list of folder paths.
    """
    candidates = []
    for entry in sorted(glob.glob(os.path.join(directory, "*"))):
        if os.path.isdir(entry) and os.path.isfile(os.path.join(entry, "calc.fdf")):
            candidates.append(entry)
    return candidates


def read_direction_result(folder, out_file, label_override=None):
    """Reads one direction folder's .EPSIMG and derives every optical
    quantity via core.optical_properties.compute_all. Returns a dict
    with 'axis' (detected or the folder's own basename as fallback),
    'label', 'scf_ok', plus every array from compute_all -- or None if
    the .EPSIMG couldn't be found/read at all (folder skipped, not a
    fatal error for the whole run -- same advisory-only convention as
    cohesive_analysis.py/oer_analysis.py).
    """
    label = label_override or detect_label(folder) or "siesta"
    axis = detect_axis(folder) or os.path.basename(folder)

    epsimg_path = os.path.join(folder, f"{label}.EPSIMG")
    if not os.path.isfile(epsimg_path):
        return None
    omega, eps2 = read_epsimg(epsimg_path)
    if len(omega) == 0:
        return None

    scf_ok, _iterations = get_scf_convergence(os.path.join(folder, out_file))
    result = compute_all(omega, eps2)
    result["axis"] = axis
    result["label"] = label
    result["scf_ok"] = scf_ok
    result["folder"] = folder
    result["siesta"] = get_optical_info(os.path.join(folder, out_file))
    return result


def trim_energy_range(result, energy_min, energy_max):
    """Filters every array in `result` (in place, returns a new dict) to
    [energy_min, energy_max] -- both optional, None means unbounded on
    that side. Applied identically to every derived quantity so the
    report/CSV/plots all stay consistent with each other.
    """
    omega = result["omega"]
    mask = np.ones_like(omega, dtype=bool)
    if energy_min is not None:
        mask &= omega >= energy_min
    if energy_max is not None:
        mask &= omega <= energy_max
    trimmed = dict(result)
    for key in ("omega", "eps1", "eps2", "n", "k", "alpha", "R", "L", "sigma1"):
        trimmed[key] = result[key][mask]
    return trimmed


def compute_isotropic_average(results):
    """Polycrystalline/isotropic-averaged dielectric response, for direct
    comparison against a typically-polycrystalline or powder experimental
    measurement (which doesn't correspond to any single crystallographic
    direction). Only computed when all 3 of xx/yy/zz are present in
    `results`. Returns a result dict shaped exactly like
    read_direction_result's own output (axis='avg'), or None if xx/yy/zz
    aren't all present.

    Averages eps2 (the fundamental linear-response quantity SIESTA
    actually computes) across the 3 directions, THEN derives every other
    quantity -- including eps1 -- via a fresh core.optical_properties.
    compute_all call on that averaged eps2. This is the physically
    correct order: n/k/alpha/R/L/sigma1 are NONLINEAR functions of
    eps1/eps2, so averaging them directly would not correspond to the
    dielectric response of any real (isotropic) medium, unlike averaging
    eps1/eps2 first. eps1 specifically is not separately averaged and
    overridden -- Kramers-Kronig is a LINEAR transform, so KK(mean(eps2))
    equals mean(KK(eps2)) up to the same small numerical KK error
    documented in core/dielectric.py, making a fresh KK-transform of the
    averaged eps2 equivalent to (and simpler than) averaging 3
    separately-computed eps1 arrays.

    If the 3 directions' energy grids don't match exactly (same length
    and values within 1e-6 -- expected when all 3 came from the same
    calc.fdf template, as stb-optical always writes, but not assumed
    blindly), the non-reference grids are interpolated (np.interp) onto
    the xx direction's own grid before averaging.
    """
    by_axis = {r["axis"]: r for r in results}
    if not all(axis in by_axis for axis in _DIAGONAL_AXES):
        return None

    ref_omega = by_axis["xx"]["omega"]
    eps2_sum = np.zeros_like(ref_omega)
    interpolated = False
    for axis in _DIAGONAL_AXES:
        r = by_axis[axis]
        if len(r["omega"]) == len(ref_omega) and np.allclose(r["omega"], ref_omega, atol=1e-6):
            eps2_sum = eps2_sum + r["eps2"]
        else:
            interpolated = True
            eps2_sum = eps2_sum + np.interp(ref_omega, r["omega"], r["eps2"])
    eps2_avg = eps2_sum / 3.0

    avg_result = compute_all(ref_omega, eps2_avg)
    avg_result["axis"] = "avg"
    avg_result["label"] = by_axis["xx"]["label"]
    avg_result["scf_ok"] = all(by_axis[axis]["scf_ok"] for axis in _DIAGONAL_AXES)
    avg_result["folder"] = "(isotropic average of xx, yy, zz)"
    avg_result["interpolated"] = interpolated
    return avg_result


def reconstruct_offdiagonal(results):
    """Off-diagonal dielectric-tensor components eps_ij(E), for every
    biaxial direction ('xy'/'xz'/'yz') present in `results` whose 2
    matching diagonal directions (see core.calc_directives.
    OPTICAL_OFFDIAG_PAIRS) are ALSO present -- eps_ij = eps_(bisector ij)
    - (eps_ii + eps_jj) / 2, applied separately to eps1(E) and eps2(E)
    (both are components of the SAME symmetric dielectric tensor, so the
    identical linear reconstruction applies to each independently; no
    Kramers-Kronig re-transform needed here, unlike compute_isotropic_
    average's own eps2-first approach, since this is a plain linear
    combination of already-consistent eps1/eps2 pairs, not a change of
    which underlying quantity is being averaged).

    Deliberately does NOT derive n/k/alpha/R/L/sigma1 for these --
    those textbook relations (n = sqrt(eps), etc.) assume a principal
    -axis (diagonal) dielectric response; an off-diagonal tensor
    component has no such standalone refractive-index/absorption
    -coefficient interpretation, so only eps1_ij(E)/eps2_ij(E) themselves
    are physically meaningful outputs here.

    Returns a list of dicts, each {'axis', 'omega', 'eps1', 'eps2',
    'scf_ok', 'interpolated', 'needs'} -- one entry per reconstructable
    biaxial direction found (possibly empty). Energy grids are matched
    onto the biaxial direction's own omega (interpolating the 2 diagonal
    directions to it if needed), same convention as compute_isotropic_
    average.
    """
    by_axis = {r["axis"]: r for r in results}
    reconstructed = []
    for mixed_axis, (need_a, need_b) in OPTICAL_OFFDIAG_PAIRS.items():
        if not (mixed_axis in by_axis and need_a in by_axis and need_b in by_axis):
            continue
        r_mixed, r_a, r_b = by_axis[mixed_axis], by_axis[need_a], by_axis[need_b]
        ref_omega = r_mixed["omega"]
        interpolated = False

        def _on_ref_grid(r):
            nonlocal interpolated
            if len(r["omega"]) == len(ref_omega) and np.allclose(r["omega"], ref_omega, atol=1e-6):
                return r["eps1"], r["eps2"]
            interpolated = True
            return (np.interp(ref_omega, r["omega"], r["eps1"]),
                    np.interp(ref_omega, r["omega"], r["eps2"]))

        eps1_a, eps2_a = _on_ref_grid(r_a)
        eps1_b, eps2_b = _on_ref_grid(r_b)
        eps1_ij = r_mixed["eps1"] - (eps1_a + eps1_b) / 2.0
        eps2_ij = r_mixed["eps2"] - (eps2_a + eps2_b) / 2.0

        reconstructed.append({
            "axis": mixed_axis,
            "omega": ref_omega,
            "eps1": eps1_ij,
            "eps2": eps2_ij,
            "scf_ok": r_mixed["scf_ok"] and r_a["scf_ok"] and r_b["scf_ok"],
            "interpolated": interpolated,
            "needs": (need_a, need_b),
        })
    return reconstructed


def write_direction_dat(dat_path, result):
    """9-column .dat for one direction: E_eV, eps1, eps2, n, k,
    alpha_cm-1, R, L, sigma1_Spm -- same commented-header convention as
    raman_analysis.py::write_spectrum_plot's own .dat writer.
    """
    with open(dat_path, "w") as f:
        f.write(f"# stb-opticalAnalysis -- optical properties along direction {result['axis']}\n")
        f.write("# columns: 1=E(eV) 2=eps1 3=eps2 4=n 5=k 6=alpha(cm^-1) 7=R 8=L 9=sigma1(S/m)\n")
        for i in range(len(result["omega"])):
            f.write(f"{result['omega'][i]:12.6f} {result['eps1'][i]:14.6f} {result['eps2'][i]:14.6f} "
                     f"{result['n'][i]:12.6f} {result['k'][i]:12.6f} {result['alpha'][i]:16.6e} "
                     f"{result['R'][i]:12.6f} {result['L'][i]:14.6e} {result['sigma1'][i]:16.6e}\n")


def write_results_csv(csv_path, results):
    """Long-format CSV, one row per (direction, energy) point -- same
    convention as gqca_analysis.py::write_results_csv.
    """
    with open(csv_path, "w") as f:
        f.write("direction,E_eV,eps1,eps2,n,k,alpha_cm-1,R,L,sigma1_Spm\n")
        for result in results:
            for i in range(len(result["omega"])):
                f.write(f"{result['axis']},{result['omega'][i]:.6f},{result['eps1'][i]:.6f},"
                        f"{result['eps2'][i]:.6f},{result['n'][i]:.6f},{result['k'][i]:.6f},"
                        f"{result['alpha'][i]:.6e},{result['R'][i]:.6f},{result['L'][i]:.6e},"
                        f"{result['sigma1'][i]:.6e}\n")


def write_offdiagonal_dat(dat_path, result):
    """3-column .dat for one reconstructed off-diagonal tensor component:
    E_eV, eps1_ij, eps2_ij -- deliberately not the 9-column format
    write_direction_dat uses (no n/k/alpha/R/L/sigma1, see reconstruct_
    offdiagonal's own docstring for why those aren't meaningful here).
    """
    with open(dat_path, "w") as f:
        f.write(f"# stb-opticalAnalysis -- reconstructed off-diagonal dielectric tensor "
                f"component eps_{result['axis']}\n")
        f.write(f"# eps_{result['axis']} = eps_(bisector {result['axis']}) - "
                f"(eps_{result['needs'][0]} + eps_{result['needs'][1]}) / 2\n")
        f.write(f"# columns: 1=E(eV) 2=eps1_{result['axis']} 3=eps2_{result['axis']}\n")
        for i in range(len(result["omega"])):
            f.write(f"{result['omega'][i]:12.6f} {result['eps1'][i]:14.6f} {result['eps2'][i]:14.6f}\n")


def write_offdiagonal_csv(csv_path, results):
    """Long-format CSV, one row per (off-diagonal direction, energy)
    point: direction,E_eV,eps1_ij,eps2_ij.
    """
    with open(csv_path, "w") as f:
        f.write("direction,E_eV,eps1_ij,eps2_ij\n")
        for result in results:
            for i in range(len(result["omega"])):
                f.write(f"{result['axis']},{result['omega'][i]:.6f},"
                        f"{result['eps1'][i]:.6f},{result['eps2'][i]:.6f}\n")


def write_offdiagonal_gplot(gplot_path, dat_paths_by_axis):
    """One gnuplot script overlaying eps1_ij/eps2_ij for every
    reconstructed off-diagonal direction present -- same 2-curve-per
    -direction convention as write_combined_gplot's own 'eps' quantity,
    kept as a separate function since the off-diagonal .dat files have a
    different (3-column, not 9-column) layout.
    """
    pdf_name = os.path.splitext(os.path.basename(gplot_path))[0] + ".pdf"
    plot_terms = []
    color_idx = 0
    for axis in sorted(dat_paths_by_axis):
        dat_path = dat_paths_by_axis[axis]
        for col, label in ((2, "eps1"), (3, "eps2")):
            color = _PALETTE[color_idx % len(_PALETTE)]
            color_idx += 1
            plot_terms.append(
                f'"{os.path.basename(dat_path)}" using 1:{col} with lines lw 2 lc rgb "{color}" '
                f'title "{label}_{axis}"')
    with open(gplot_path, "w") as f:
        f.writelines([
            '# --- STB Plot Configuration ---\n',
            '# Generated by stb-opticalAnalysis\n',
            'set terminal pdfcairo enhanced color font "Arial,14" size 8,5\n',
            f'set output "{pdf_name}"\n\n',
            'set title "Reconstructed off-diagonal dielectric tensor components"\n',
            'set xlabel "Photon energy (eV)"\n',
            'set ylabel "eps_{ij}"\n',
            'set grid\n',
            'set key top right\n',
            'plot ' + ', \\\n     '.join(plot_terms) + '\n',
        ])


def write_polarizability_dat(dat_path, axis, omega, alpha1_si, alpha2_si, alpha1_ang3, alpha2_ang3):
    """5-column .dat for one direction's molecular polarizability (0D
    correction): E_eV, alpha1_SI, alpha2_SI, alpha1_Ang3, alpha2_Ang3.
    """
    with open(dat_path, "w") as f:
        f.write(f"# stb-opticalAnalysis -- molecular polarizability along direction {axis} "
                "(dilute Clausius-Mossotti correction)\n")
        f.write("# columns: 1=E(eV) 2=alpha1_SI(C*m^2/V) 3=alpha2_SI(C*m^2/V) 4=alpha1(Ang^3) "
                "5=alpha2(Ang^3)\n")
        for i in range(len(omega)):
            f.write(f"{omega[i]:12.6f} {alpha1_si[i]:16.6e} {alpha2_si[i]:16.6e} "
                     f"{alpha1_ang3[i]:14.6f} {alpha2_ang3[i]:14.6f}\n")


def write_polarizability_csv(csv_path, rows):
    """Long-format CSV, one row per (direction, energy) point, for the
    0D molecular-polarizability correction.
    """
    with open(csv_path, "w") as f:
        f.write("direction,E_eV,alpha1_SI,alpha2_SI,alpha1_Ang3,alpha2_Ang3\n")
        for axis, omega, a1si, a2si, a1a3, a2a3 in rows:
            for i in range(len(omega)):
                f.write(f"{axis},{omega[i]:.6f},{a1si[i]:.6e},{a2si[i]:.6e},"
                        f"{a1a3[i]:.6f},{a2a3[i]:.6f}\n")


def write_combined_gplot(gplot_path, dat_paths_by_axis, quantity, experimental=None,
                          experimental_dat_path=None):
    """One gnuplot script overlaying `quantity`'s column(s) across every
    direction present -- e.g. --plot-quantity eps plots eps1 AND eps2
    for every direction (paired quantities that are conventionally shown
    together, same as n_k), while alpha/R/L/sigma1 each plot a single
    curve per direction. Lets the user directly compare anisotropy (e.g.
    eps2 along x vs. z for a 2D material) in one figure.

    `experimental` (optional): (energy_ev, value) arrays of a real
    measured spectrum to overlay -- written to `experimental_dat_path`
    (required together with `experimental`) and plotted against a
    SEPARATE right-hand y-axis (gnuplot y2), same rationale as
    raman_analysis.py::write_spectrum_plot's own experimental overlay:
    simulated and measured values are generally on different scales
    (e.g. absorption coefficient in cm^-1 vs. an arbitrary-units
    spectrophotometer reading) -- a shared single axis would make one
    curve invisible next to the other. Only meaningful for a single-
    column `quantity` (alpha/R/L/sigma1) -- callers must not pass
    `experimental` together with a paired quantity (eps/n_k); this
    function does not itself validate that (see main()'s own guard).
    """
    columns = _QUANTITY_COLUMNS[quantity]
    title = _QUANTITY_TITLES[quantity]
    pdf_name = os.path.splitext(os.path.basename(gplot_path))[0] + ".pdf"

    plot_terms = []
    color_idx = 0
    for axis in sorted(dat_paths_by_axis):
        dat_path = dat_paths_by_axis[axis]
        for _key, col, label in columns:
            color = _PALETTE[color_idx % len(_PALETTE)]
            color_idx += 1
            plot_terms.append(
                f'"{os.path.basename(dat_path)}" using 1:{col} with lines lw 2 lc rgb "{color}" '
                f'axes x1y1 title "{label} ({axis})"' if experimental is not None else
                f'"{os.path.basename(dat_path)}" using 1:{col} with lines lw 2 lc rgb "{color}" '
                f'title "{label} ({axis})"')

    if experimental is not None:
        exp_energy, exp_value = experimental
        with open(experimental_dat_path, "w") as f:
            f.write("# stb-opticalAnalysis -- user-supplied experimental spectrum\n")
            f.write("# columns: 1=E(eV) 2=value (original scale)\n")
            for e, v in zip(exp_energy, exp_value):
                f.write(f"{e:12.6f} {v:14.6f}\n")
        plot_terms.append(
            f'"{os.path.basename(experimental_dat_path)}" using 1:2 with lines lw 2 '
            'lc rgb "#333333" axes x1y2 title "Experimental"')

    with open(gplot_path, "w") as f:
        lines = [
            '# --- STB Plot Configuration ---\n',
            '# Generated by stb-opticalAnalysis\n',
            'set terminal pdfcairo enhanced color font "Arial,14" size 8,5\n',
            f'set output "{pdf_name}"\n\n',
            f'set title "{title}"\n',
            'set xlabel "Photon energy (eV)"\n',
            'set ylabel "' + title + '"\n',
            'set grid\n',
        ]
        if experimental is not None:
            lines += ['set y2label "Experimental (arb. units)"\n', 'set y2tics\n']
        lines += [
            'set key top right\n',
            'plot ' + ', \\\n     '.join(plot_terms) + '\n',
        ]
        f.writelines(lines)


# --- Analyses and dashboard ---------------------------------------------------------


def reference_result(results):
    """The spectrum the single-spectrum tables use: the isotropic average when xx/yy/zz are all
    present, else the first diagonal direction, else the first result."""
    for r in results:
        if r["axis"] == "avg":
            return r
    for r in results:
        if r["axis"] in _DIAGONAL_AXES:
            return r
    return results[0]


def static_properties(result, volume_ang3, n_valence):
    """eps1(0), n(0), R(0), absorption onset, Tauc gap and the Penn gap of one spectrum."""
    omega = result["omega"]
    i0 = int(np.argmin(omega))
    eps1_0, n0 = float(result["eps1"][i0]), float(result["n"][i0])
    onset = absorption_onset(omega, result["eps2"])
    peaks = find_peaks_table(omega, result["eps2"], 1, 0.05, 0.05)
    first_peak = peaks[0]["E"] if peaks else float("nan")
    tauc, fit, window = tauc_direct_gap(omega, result["alpha"], onset, first_peak)
    plasmon = free_electron_plasmon_ev(n_valence, volume_ang3) if n_valence and volume_ang3 else float("nan")
    return {"eps1_0": eps1_0, "n0": n0, "R0": float(((n0 - 1.0) / (n0 + 1.0)) ** 2),
            "onset": onset, "tauc": tauc, "tauc_fit": fit, "tauc_window": window,
            "penn": penn_gap_ev(eps1_0, plasmon) if np.isfinite(plasmon) else float("nan"),
            "plasmon_free": plasmon}


def spectral_features(result, n_peaks, sigma_ev, e_max):
    """Peaks of eps2, the absorption/reflectivity maxima, the loss-function peak (bulk plasmon
    candidate) and the zero crossings of eps1, for one spectrum."""
    omega = result["omega"]
    features = {"eps2_peaks": find_peaks_table(omega, result["eps2"], n_peaks, sigma_ev, 0.05, e_max)}
    mask = omega <= (e_max if e_max else omega.max())
    for key in ("alpha", "R", "eps1"):
        arr = np.where(mask, result[key], -np.inf)
        j = int(np.argmax(arr))
        features[key + "_max"] = (float(omega[j]), float(result[key][j]))
    loss_peaks = find_peaks_table(omega, result["L"], 1, sigma_ev, 0.2)
    features["loss_peak"] = loss_peaks[0] if loss_peaks else None
    features["eps1_zeros"] = [(e, d) for e, d in zero_crossings(omega, result["eps1"], sigma_ev) if e <= omega.max()]
    return features


def interp_at(result, key, energy):
    return float(np.interp(energy, result["omega"], result[key]))


def make_dashboard(results, ref, stats, offdiag_results, scissor_result, neff, n_valence, plot_emax,
                    png_path, show, experimental=None, exp_quantity=None, title=""):
    """3x3 matplotlib figure that completes the analysis visually: eps1/eps2, n/k, absorption
    coefficient with its onset, the direct-gap Tauc plot, reflectivity, energy-loss function with the
    plasmon, optical conductivity, the f-sum-rule accumulation N_eff(E), and either the scissor
    what-if, the reconstructed off-diagonal components, or the anisotropy. Saved to `png_path` and/or
    shown on screen. Returns True when something was produced."""
    try:
        import matplotlib
        if not show:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print(color_text("[WARNING] matplotlib is not installed -- skipping the dashboard.", 'yellow'))
        return False
    if show and (matplotlib.get_backend().lower() == "agg"
                 or (sys.platform.startswith("linux") and not (os.environ.get("DISPLAY")
                                                               or os.environ.get("WAYLAND_DISPLAY")))):
        print(color_text("[WARNING] --show needs a display (local desktop or `ssh -X`); "
                          "only saving the figure.", 'yellow'))
        show = False
        plt.switch_backend("Agg")

    okabe = {"xx": "#0072B2", "yy": "#D55E00", "zz": "#009E73", "avg": "#000000",
             "xy": "#CC79A7", "xz": "#E69F00", "yz": "#56B4E9"}
    plotted = [r for r in results if r["axis"] in _DIAGONAL_AXES or r["axis"] == "avg"] or results
    fig, axes = plt.subplots(3, 3, figsize=(16, 12))
    ax = axes.ravel()

    def style(a, xlabel="Photon energy (eV)", ylabel=""):
        a.set_xlim(0, plot_emax)
        a.set_xlabel(xlabel)
        a.set_ylabel(ylabel)
        a.grid(alpha=0.3)

    def lw_of(r):
        return 2.4 if r["axis"] == "avg" else 1.3

    for r in plotted:
        c, lab = okabe.get(r["axis"], "0.4"), r["axis"]
        ax[0].plot(r["omega"], r["eps1"], color=c, lw=lw_of(r), label=f"eps1 {lab}")
        ax[0].plot(r["omega"], r["eps2"], color=c, lw=lw_of(r), ls="--", label=f"eps2 {lab}")
        ax[1].plot(r["omega"], r["n"], color=c, lw=lw_of(r), label=f"n {lab}")
        ax[1].plot(r["omega"], r["k"], color=c, lw=lw_of(r), ls="--", label=f"k {lab}")
        ax[2].semilogy(r["omega"], np.clip(r["alpha"], 1.0, None), color=c, lw=lw_of(r), label=lab)
        ax[4].plot(r["omega"], r["R"], color=c, lw=lw_of(r), label=lab)
        ax[5].plot(r["omega"], r["L"], color=c, lw=lw_of(r), label=lab)
        ax[6].plot(r["omega"], r["sigma1"] / 1e6, color=c, lw=lw_of(r), label=lab)
    ax[0].axhline(0, color="0.5", lw=0.8)
    style(ax[0], ylabel="dielectric function"); ax[0].set_title("Dielectric function")
    style(ax[1], ylabel="n, k"); ax[1].set_title("Refractive index (solid) / extinction (dashed)")
    ax[1].axhline(0, color="0.5", lw=0.8)
    style(ax[2], ylabel="alpha (cm$^{-1}$)"); ax[2].set_title("Absorption coefficient")
    ax[2].set_ylim(1e2, max(1e6, 2 * float(np.max([r["alpha"].max() for r in plotted]))))
    onset = stats["onset"]
    if np.isfinite(onset):
        ax[2].axvline(onset, color="0.4", ls=":", lw=1.2)
        ax[2].annotate(f"onset {onset:.2f} eV", (onset, ax[2].get_ylim()[1] / 3), rotation=90,
                       va="top", ha="right", fontsize=9, color="0.3")
    style(ax[4], ylabel="R"); ax[4].set_title("Normal-incidence reflectivity"); ax[4].set_ylim(0, 1)
    style(ax[5], ylabel="-Im(1/eps)"); ax[5].set_title("Energy-loss function")
    lp = stats.get("loss_peak")
    if lp:
        ax[5].plot([lp["E"]], [lp["value"]], "o", color="#CC0000")
        ax[5].annotate(f"{lp['E']:.2f} eV", (lp["E"], lp["value"]), textcoords="offset points",
                       xytext=(6, 4), fontsize=9)
    if np.isfinite(stats["plasmon_free"]):
        ax[5].axvline(stats["plasmon_free"], color="0.4", ls=":", lw=1.2,
                      label=f"free-electron {stats['plasmon_free']:.1f} eV")
    style(ax[6], ylabel="sigma1 (10$^6$ S/m)"); ax[6].set_title("Real optical conductivity")
    for a in (ax[0], ax[1], ax[2], ax[4], ax[5], ax[6]):
        a.legend(fontsize=7, ncol=2)

    # Tauc plot of the reference spectrum
    y = (ref["alpha"] * ref["omega"]) ** 2
    ax[3].plot(ref["omega"], y / 1e12, color="#0072B2", lw=1.6)
    if stats["tauc_fit"] is not None and np.isfinite(stats["tauc"]):
        slope, intercept = stats["tauc_fit"]
        xs = np.linspace(stats["tauc"], stats["tauc_window"][1] + 0.3, 20)
        ax[3].plot(xs, (slope * xs + intercept) / 1e12, color="#CC0000", lw=1.3, ls="--",
                   label=f"fit -> E$_g^{{Tauc}}$ = {stats['tauc']:.2f} eV")
        ax[3].axvline(stats["tauc"], color="#CC0000", ls=":", lw=1)
        ax[3].legend(fontsize=8)
    style(ax[3], ylabel=r"($\alpha E$)$^2$ (10$^{12}$)"); ax[3].set_title(f"Direct-gap Tauc plot ({ref['axis']})")
    ax[3].set_xlim(max(stats["onset"] - 0.6, 0) if np.isfinite(stats["onset"]) else 0,
                   min(plot_emax, (stats["onset"] + 3.0) if np.isfinite(stats["onset"]) else plot_emax))
    ymax_t = y[(ref["omega"] >= ax[3].get_xlim()[0]) & (ref["omega"] <= ax[3].get_xlim()[1])]
    if len(ymax_t):
        ax[3].set_ylim(0, 1.1 * float(ymax_t.max()) / 1e12)

    # f-sum accumulation
    if neff is not None:
        ax[7].plot(ref["omega"], neff, color="#0072B2", lw=1.8, label=f"N$_{{eff}}$({ref['axis']})")
        if n_valence:
            ax[7].axhline(n_valence, color="0.4", ls="--", lw=1.2, label=f"valence electrons = {n_valence:g}")
        ax[7].legend(fontsize=8)
    style(ax[7], ylabel="electrons per cell"); ax[7].set_title("f-sum rule: transitions accounted for")
    ax[7].set_xlim(0, max(plot_emax, 30))

    # panel 9: scissor what-if / off-diagonal / anisotropy
    if scissor_result is not None:
        ax[8].plot(ref["omega"], ref["eps2"], color="0.4", lw=1.4, label="eps2 (DFT)")
        ax[8].plot(scissor_result["omega"], scissor_result["eps2"], color="#CC0000", lw=1.6, label="eps2 (scissor)")
        ax[8].plot(ref["omega"], ref["eps1"], color="0.6", lw=1.0, ls="--", label="eps1 (DFT)")
        ax[8].plot(scissor_result["omega"], scissor_result["eps1"], color="#CC79A7", lw=1.2, ls="--", label="eps1 (scissor)")
        ax[8].set_title(f"Scissor what-if ({scissor_result['delta']:+.2f} eV)")
        style(ax[8], ylabel="dielectric function"); ax[8].legend(fontsize=8)
    elif offdiag_results:
        for r in offdiag_results:
            c = okabe.get(r["axis"], "0.4")
            ax[8].plot(r["omega"], r["eps1"], color=c, lw=1.3, label=f"eps1_{r['axis']}")
            ax[8].plot(r["omega"], r["eps2"], color=c, lw=1.3, ls="--", label=f"eps2_{r['axis']}")
        ax[8].set_title("Off-diagonal tensor components (reconstructed)")
        style(ax[8], ylabel="eps$_{ij}$"); ax[8].legend(fontsize=8)
    else:
        diag = [r for r in results if r["axis"] in _DIAGONAL_AXES]
        if len(diag) > 1:
            mean2 = np.mean([r["eps2"] for r in diag], axis=0)
            for r in diag:
                ax[8].plot(r["omega"], r["eps2"] - mean2, color=okabe.get(r["axis"], "0.4"), lw=1.3,
                           label=f"eps2_{r['axis']} - mean")
            ax[8].set_title("Anisotropy: eps2 minus the directional mean")
            style(ax[8], ylabel="Delta eps2"); ax[8].legend(fontsize=8)
        else:
            ax[8].text(0.5, 0.5, "single direction:\nno anisotropy to show", ha="center", va="center")
            ax[8].set_axis_off()

    if experimental is not None and exp_quantity is not None:
        panel = {"alpha": 2, "R": 4, "L": 5, "sigma1": 6}.get(exp_quantity)
        if panel is not None:
            twin = ax[panel].twinx()
            twin.plot(experimental[0], experimental[1], color="#333333", lw=1.4, ls="-.", label="experimental")
            twin.set_ylabel("experimental (arb. units)")
            twin.legend(fontsize=8, loc="lower right")
    fig.suptitle(title or "Optical properties", fontsize=14)
    fig.tight_layout()
    if png_path:
        fig.savefig(png_path, dpi=150)
    if show:
        plt.show()
    plt.close(fig)
    return True


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Stage 2 of 2: aggregates every direction folder written by "
        "stb-optical and derives the full set of linear optical properties.", 'bold')}
Scans --directory for every subfolder containing a calc.fdf, identifies each one's direction from its
own %block Optical.Vector (not the folder name -- robust to renamed/reorganized folders; 'xx'/'yy'/
'zz' diagonal or 'xy'/'xz'/'yz' biaxial, the same convention stb-optical writes), reads its
SystemLabel.EPSIMG, and derives eps1(E)/eps2(E) (Kramers-Kronig transform), n(E)/k(E) (refractive
index/extinction coefficient), the absorption coefficient, normal-incidence reflectivity, the
electron energy-loss function, and the real optical conductivity -- all standard textbook relations
(Wooten/Fox, "Optical Properties of Solids"), not flagged [UNVERIFIED]. Missing or incomplete
direction folders are skipped with a warning, never block the directions that ARE ready (advisory
only, same convention as this suite's other multi-folder analysis stages).

[OFF-DIAGONAL RECONSTRUCTION] Whenever a biaxial direction ('xy'/'xz'/'yz') AND its 2 matching
diagonal directions are all present, also reconstructs the true off-diagonal dielectric-tensor
component eps_ij(E) = eps_(bisector ij)(E) - (eps_ii(E) + eps_jj(E)) / 2 (same bisector-based formula
raman_analysis.py's own raman_tensor_full uses for the Raman tensor's Rxy/Rxz/Ryz) and writes it to
its own dedicated output -- eps1_ij/eps2_ij only, NOT n/k/alpha/R/L/sigma1 (those textbook relations
assume a principal-axis/diagonal response and have no standalone meaning for a genuine off-diagonal
tensor component).

[ASSUMPTION] Reflectivity R(E) is normal-incidence, vacuum(n0=1)/material interface -- no angle-
dependent Fresnel formula or thin-film interference is computed. [KNOWN LIMITATION, 2D/1D/0D inputs]
The raw supercell values above are diluted by the vacuum region and NOT corrected by default -- pass
--dimensionality-correction (see its own --help entry) to additionally restore the intrinsic 2D
dielectric function or extract a 0D molecular polarizability, written ALONGSIDE the raw values, never
replacing them (1D is not yet implemented -- see --dimensionality-correction's own help text; this
correction is only ever applied to the 3 diagonal directions, never to a biaxial/reconstructed one).

When all 3 of xx/yy/zz are present, also computes and reports/plots an isotropic ('avg') average --
eps2 averaged across the 3 directions first, then every other quantity (including eps1) re-derived
from that averaged eps2 -- the standard convention for comparing against a polycrystalline/powder
experimental measurement, which doesn't correspond to any single crystallographic direction.

[REPORT] System/ground state and per-folder Optical settings read back from the folders; data quality
(SCF, SIESTA's f-sum check, the f-sum rule recomputed from eps2 as N_eff(E)); static eps1/n/R, absorption
onset, Tauc direct gap and Penn gap; peaks of eps2 and the plasmon; anisotropy and a symmetry-consistency
check of the tensor; optical constants at chosen energies; a plain-language interpretation. Optional:
--show/--png (9-panel matplotlib dashboard), --scissor (post-processing what-if for the DFT gap error),
--energies, --n-peaks.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s --directory optical_study --plot-quantity alpha\n"
               "  %(prog)s --directory optical_study --plot-quantity alpha "
               "--experimental uvvis.dat\n"
               "  %(prog)s --directory optical_study --plot-quantity all\n"
               "  %(prog)s --directory optical_study --dimensionality-correction --thickness 3.13\n"
    )

    parser.add_argument("--directory", type=str, default="optical_study",
                         help="Root directory written by stb-optical (default: optical_study).")
    parser.add_argument("--file", type=str, default="calc.out",
                         help="SIESTA output filename inside each direction folder, used only for "
                              "the SCF-convergence advisory check (default: calc.out).")
    parser.add_argument("-l", "--label", type=str, default=None,
                         help="SystemLabel override (default: auto-detected per folder).")
    parser.add_argument("--energy-min", type=float, default=None, metavar="EV",
                         help="Trim the reported/plotted energy range from below (default: no trim).")
    parser.add_argument("--energy-max", type=float, default=None, metavar="EV",
                         help="Trim the reported/plotted energy range from above (default: no trim).")
    parser.add_argument("--plot-quantity", type=str, default="eps",
                         choices=["eps", "n_k", "alpha", "R", "L", "sigma1", "all"],
                         help="Which quantity the combined .gplot overlays across every direction "
                              "present (default: eps -- plots eps1 AND eps2 together). 'all' writes "
                              "one .gplot per quantity in a single run (mutually exclusive with "
                              "--experimental -- pick one specific quantity for a direct comparison "
                              "instead).")
    parser.add_argument("--experimental", type=str, default=None, metavar="PATH",
                         help="2-column (energy_eV, value) experimental spectrum (e.g. UV-Vis-NIR "
                              "absorbance/reflectance) to overlay on the plot -- only valid with a "
                              "single-valued --plot-quantity (alpha/R/L/sigma1), not eps/n_k/all. "
                              "Peak positions are matched against the isotropic average (if all 3 "
                              "directions are present) or the first available direction otherwise.")
    parser.add_argument("--peak-prominence", type=float, default=0.01, metavar="FRACTION",
                         help="Minimum peak prominence, as a fraction of that spectrum's own max "
                              "value, for --experimental's peak-finding (default: 0.01).")
    parser.add_argument("--dimensionality-correction", action="store_true",
                         help="Restore the intrinsic (vacuum-independent) response from the raw "
                              "vacuum-padded supercell result. 2D (1 vacuum axis): rescales eps(E) "
                              "via the D-field/E-field continuity argument (Laturia et al., npj 2D "
                              "Mater. Appl. 2, 6, 2018; Yang & Gao, npj 2D Mater. Appl. 5, 78, 2021) "
                              "-- requires --thickness. 0D (3 vacuum axes, an isolated molecule): "
                              "extracts the molecular polarizability alpha(E) via the dilute-limit "
                              "Clausius-Mossotti relation -- no extra input needed (uses the cell "
                              "volume). 1D (2 vacuum axes): NOT implemented -- the literature "
                              "approach (Maxwell-Garnett effective medium) is not a simple "
                              "generalization of the 2D formula and has not been verified for this "
                              "tool; prints a [NOTE] and reports only the raw supercell values. 3D "
                              "(bulk): no-op, nothing to correct. Off by default -- always written "
                              "ALONGSIDE (never replacing) the raw supercell .dat/.csv files.")
    parser.add_argument("--thickness", type=float, default=None, metavar="ANG",
                         help="Physical thickness of the 2D material in Ang, required together "
                              "with --dimensionality-correction for a 2D (1 vacuum axis) input. "
                              "Deliberately not auto-estimated -- 'thickness' of a 2D material is "
                              "not a single well-defined quantity (atomic span? + van der Waals "
                              "radii? bulk interlayer spacing?), left to the user's own judgment.")
    parser.add_argument("-o", "--output", type=str, default="optical_results",
                         help="Base filename (no extension) for the .csv/.dat/.gplot outputs "
                              "(default: optical_results).")
    parser.add_argument("--save-report", action="store_true",
                         help=f"Also persist the report to <directory>/{REPORT_FILE}. Off by default.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-opticalAnalysis {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    parser.add_argument("--show", action="store_true",
                         help="Open the analysis dashboard (matplotlib, 9 panels) on screen. Needs a display.")
    parser.add_argument("--png", action="store_true",
                         help="Save the dashboard as <directory>/<output>_summary.png.")
    parser.add_argument("--plot-emax", type=float, default=None, metavar="EV",
                         help="Upper energy of the dashboard's x axes (default: the data range, at most 20 eV). "
                              "--energy-min/--energy-max instead trim the DATA used by every analysis.")
    parser.add_argument("--energies", type=float, nargs="+", default=[0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0],
                         metavar="EV",
                         help="Photon energies (eV) of the 'optical constants at selected energies' table "
                              "(0 = static limit; default: 0 1 2 3 4 5 10).")
    parser.add_argument("--n-peaks", type=int, default=6, metavar="N",
                         help="Number of eps2 peaks listed in the spectral-features table (default: 6).")
    parser.add_argument("--smoothing", type=float, default=0.05, metavar="EV",
                         help="Gaussian sigma (eV) used only to locate peaks/zero crossings in a spiky "
                              "spectrum (default: 0.05).")
    parser.add_argument("--scissor", type=float, default=None, metavar="EV",
                         help="What-if scissor shift (eV) applied in post-processing to the reference spectrum: "
                              "eps2'(E) = eps2(E - EV), eps1 re-derived by Kramers-Kronig; reports the effect on the "
                              "static constant, edge, peaks and plasmon and overlays it on the dashboard.")
    args = parser.parse_args()

    if args.experimental is not None and args.plot_quantity in ("eps", "n_k", "all"):
        print(color_text(
            "[ERROR] --experimental requires a single-valued --plot-quantity (alpha/R/L/sigma1) "
            "-- eps/n_k are paired curves and 'all' writes multiple plots, both ambiguous to "
            "overlay a single experimental curve onto.", 'red'))
        sys.exit(1)
    if args.experimental is not None and not os.path.isfile(args.experimental):
        print(color_text(f"[ERROR] --experimental file '{args.experimental}' not found.", 'red'))
        sys.exit(1)

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("OPTICAL PROPERTIES WORKFLOW -- STAGE 2: ANALYSIS", 'bold'))
    print("-" * 60)

    if not os.path.isdir(args.directory):
        print(color_text(f"[ERROR] Directory '{args.directory}' not found -- run stb-optical first.", 'red'))
        sys.exit(1)

    folders = find_direction_folders(args.directory)
    if not folders:
        print(color_text(
            f"[ERROR] No direction folder (a subdirectory with a calc.fdf) found under "
            f"'{args.directory}' -- run stb-optical first.", 'red'))
        sys.exit(1)

    report_path = os.path.join(args.directory, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, 'w') if report_path else None

    warnings = []                       # collected, repeated together at the end

    def warn(text):
        warnings.append(text)
        return color_text(f"[WARNING] {text}", 'yellow')

    section_counter = [0]

    def section(title):
        print_section(f"[{section_counter[0]}] {title}", f_out)
        section_counter[0] += 1

    print_dual(f"{color_text('===== OPTICAL STAGE 2 REPORT (ANALYSIS) =====', 'magenta')}", f_out)

    section('RUN METADATA')
    print_dual(f"Tool            : stb-opticalAnalysis {VERSION}", f_out)
    print_dual(f"Date            : {datetime.datetime.now():%Y-%m-%d %H:%M}", f_out)
    print_dual(f"Directory       : {args.directory}", f_out)
    print_dual(f"Folders found   : {len(folders)}", f_out)
    print_dual(f"SIESTA log      : <folder>/{args.file}" + ("  (default)" if args.file == "calc.out" else ""), f_out)
    window = []
    if args.energy_min is not None:
        window.append(f">= {args.energy_min:g} eV")
    if args.energy_max is not None:
        window.append(f"<= {args.energy_max:g} eV")
    print_dual(f"Energy window   : {' and '.join(window) if window else 'whole EPSIMG range (no trim)'}", f_out)

    results = []
    seen_axes = {}
    for folder in folders:
        result = read_direction_result(folder, args.file, args.label)
        if result is None:
            print_dual(color_text(
                f"  [WARNING] Could not read '.EPSIMG' from '{folder}' -- skipped.", 'yellow'), f_out)
            continue
        if args.energy_min is not None or args.energy_max is not None:
            result = trim_energy_range(result, args.energy_min, args.energy_max)
            if len(result["omega"]) == 0:
                print_dual(color_text(
                    f"  [WARNING] '{folder}' has no data points left after --energy-min/"
                    "--energy-max trimming -- skipped.", 'yellow'), f_out)
                continue
        if result["axis"] in seen_axes:
            print_dual(color_text(
                f"  [WARNING] Both '{seen_axes[result['axis']]}' and '{folder}' were detected "
                f"as direction '{result['axis']}' -- keeping the latter, the former's own "
                "output file(s) will be silently overwritten below. Rename/remove one of them "
                "if this is unintended.", 'yellow'), f_out)
        seen_axes[result["axis"]] = folder
        results.append(result)

    if not results:
        print_dual(color_text(
            "\n[ERROR] No usable direction results -- run SIESTA in every "
            f"'{args.directory}/*/' folder, then re-run stb-opticalAnalysis.", 'red'))
        if f_out:
            f_out.close()
        sys.exit(1)

    results.sort(key=lambda r: (_AXIS_ORDER.get(r["axis"], 99), r["axis"]))
    print_dual(f"Usable folders  : {len(results)} ({', '.join(r['axis'] for r in results)})", f_out)

    avg_result = compute_isotropic_average(results)
    if avg_result is not None:
        if avg_result["interpolated"]:
            print_dual(color_text(
                "  [NOTE] Direction energy grids did not match exactly -- interpolated onto "
                "the xx direction's grid before averaging.", 'yellow'), f_out)
        results.append(avg_result)

    # ---- structure, symmetry, ground state
    structure, vacuum_axes, lattice, volume = None, None, None, None
    struct_path = os.path.join(results[0]["folder"], "structure.fdf")
    if os.path.isfile(struct_path):
        structure = structure_io.read_fdf(struct_path)
        lattice = np.array(structure.lattice, dtype=float)
        volume = float(abs(np.linalg.det(lattice)))
        positions = np.array([pos for _, pos in structure.atoms])
        is_cartesian = structure.coord_format == 'cartesian'
        frac_coords = kspace.to_fractional(positions, structure.lattice, is_cartesian)
        vacuum_axes = kspace.detect_vacuum_axes(frac_coords, structure.lattice, 10.0)
    ref = reference_result(results)
    siesta_infos = [r.get("siesta") or {} for r in results if r["axis"] != "avg"]
    n_valence = next((i["nelectrons"] for i in siesta_infos if "nelectrons" in i), None)
    n_valence_known = n_valence is not None
    gap_info = {}
    for r in results:
        if r["axis"] != "avg":
            gap_info = read_eig_gap(os.path.join(r["folder"], f"{r['label']}.EIG"))
            if gap_info:
                break
    metallic = bool(gap_info.get("metallic"))
    is_vacuum_cell = bool(vacuum_axes is not None and any(vacuum_axes))
    volume_density = None if is_vacuum_cell else volume   # free-electron plasmon/Penn need the real electron density
    diag_results = [r for r in results if r["axis"] in _DIAGONAL_AXES]
    plotted_results = [r for r in results if r["axis"] in _DIAGONAL_AXES or r["axis"] == "avg"] or results

    # ================================================================ system & calculation
    section('SYSTEM & CALCULATION (read back from the folders)')
    if structure is not None:
        counts = structure_io.atom_counts(structure)
        formula = " ".join(f"{s}{counts[s]}" for s in structure_io.species_list(structure))
        lengths = [float(np.linalg.norm(v)) for v in lattice]
        print_dual(f"Structure       : {formula}  ({len(structure.atoms)} atoms), cell volume {volume:.4f} Ang^3", f_out)
        print_dual(f"Lattice (Ang)   : a={lengths[0]:.4f}  b={lengths[1]:.4f}  c={lengths[2]:.4f}", f_out)
        print_dual(f"Detected : {kspace.dimensionality_label(vacuum_axes)}", f_out)
        if any(vacuum_axes):
            print_dual(color_text(
                "[KNOWN LIMITATION] Vacuum-padded (2D/slab) input -- see stb-optical --help; "
                "every direction's dielectric response below is diluted by the vacuum region, "
                "no rescaling applied unless --dimensionality-correction is given.", 'yellow'), f_out)
        if n_valence_known:
            print_dual(f"Valence electrons : {n_valence:g} (from SIESTA's optical log)", f_out)
    else:
        print_dual("(structure.fdf not found in the first folder: cell-dependent analyses are skipped)", f_out)

    first_folder = results[0]["folder"]
    folder_texts = []
    for name in ("config_extra.fdf", "calc.fdf"):
        path = os.path.join(first_folder, name)
        if os.path.isfile(path):
            with open(path, 'r', encoding='utf-8', errors='ignore') as fh:
                folder_texts.append(fh.read())
    template_values = fdf_directives.parse_directives(*folder_texts)
    template_kgrid = fdf_directives.parse_kgrid(*folder_texts)
    theory_keys = [("xcfunctional", "XC functional"), ("xcauthors", "XC authors"),
                   ("paobasissize", "Basis size"), ("paoenergyshift", "PAO energy shift"),
                   ("meshcutoff", "Mesh cutoff"), ("spin", "Spin"),
                   ("electronictemperature", "Electronic temperature")]
    theory_rows = [([label, template_values[key]], None) for key, label in theory_keys if key in template_values]
    if template_kgrid:
        theory_rows.append((["SCF k-grid", "x".join(str(v) for v in template_kgrid)], None))
    if theory_rows:
        print_dual("\nGround state (from the first folder's calc.fdf + config_extra.fdf):", f_out)
        print_table(["Setting", "Value"], theory_rows, f_out)

    optical_rows, signature = [], {}
    for r in results:
        if r["axis"] == "avg":
            continue
        info = r.get("siesta") or {}
        mesh = info.get("mesh")
        sp = ""
        if mesh and lattice is not None:
            spacings = reciprocal_spacings(lattice, mesh)
            worst = max(s for i, s in enumerate(spacings) if not (vacuum_axes and vacuum_axes[i])) \
                if not (vacuum_axes and all(vacuum_axes)) else 0.0
            sp = f"{'x'.join(str(m) for m in mesh)} ({spacing_band(worst)})"
        elif mesh:
            sp = "x".join(str(m) for m in mesh)
        broaden_ev = info.get("broaden_ry", float('nan')) * RY_TO_EV if "broaden_ry" in info else float('nan')
        scissor_ev = info.get("scissor_ry", 0.0) * RY_TO_EV if "scissor_ry" in info else float('nan')
        optical_rows.append(([r["axis"], os.path.basename(r["folder"]), sp or "?",
                               f"{broaden_ev:.3f}" if np.isfinite(broaden_ev) else "?",
                               str(info.get("nbands", "?")), f"{scissor_ev:.3f}" if np.isfinite(scissor_ev) else "?",
                               info.get("polarization", "?")], None))
        signature[r["axis"]] = (info.get("mesh"), info.get("nbands"), round(broaden_ev, 4) if np.isfinite(broaden_ev) else None,
                                 round(scissor_ev, 4) if np.isfinite(scissor_ev) else None)
    if optical_rows:
        print_dual("\nOptical settings (as SIESTA ran them):", f_out)
        print_table(["Dir", "Folder", "Optical.Mesh (spacing)", "Broaden (eV)", "Bands", "Scissor (eV)", "Polarization"],
                    optical_rows, f_out)
        if len(set(signature.values())) > 1:
            print_dual(warn("The direction folders were NOT run with identical Optical.Mesh / Broaden / "
                             "NumberOfBands / Scissor: eps(E) of different directions is then not directly "
                             "comparable (anisotropy and the eps_ij reconstruction are affected)."), f_out)
    sym = symmetry_analysis(structure) if structure is not None else None
    if sym is not None:
        print_dual(f"\nSymmetry        : {sym['symbol']} (No. {sym['number']}), {sym['system']}; "
                    f"independent dielectric-tensor components: {sym['independent']} of 6", f_out)

    # ================================================================ data quality
    section('DATA QUALITY & CONVERGENCE')
    quality_rows, neff_by_axis = [], {}
    for r in results:
        omega = r["omega"]
        info = r.get("siesta") or {}
        step = float(np.median(np.diff(omega))) if len(omega) > 1 else float('nan')
        neff = None
        if volume is not None:
            neff = effective_electrons(omega, r["eps2"], volume)
            neff_by_axis[r["axis"]] = neff
        frac = (neff[-1] / n_valence) if (neff is not None and n_valence_known and n_valence) else None
        scf_txt = "yes" if r["scf_ok"] else "NO"
        quality_rows.append(([r["axis"], scf_txt,
                               f"{info['fsum']:.4f}" if "fsum" in info else "-",
                               f"{neff[-1]:.3f}" if neff is not None else "-",
                               f"{frac:.3f}" if frac is not None else "-",
                               f"{omega.min():.2f}-{omega.max():.1f}", str(len(omega)),
                               f"{step:.4f}"], None if r["scf_ok"] else 'yellow'))
        if not r["scf_ok"] and r["axis"] != "avg":
            print_dual(warn(f"SCF convergence not confirmed for direction {r['axis']} ('{r['folder']}'): "
                             "its spectrum may be unreliable."), f_out)
    fsum_low = {}
    for r in results:
        info = r.get("siesta") or {}
        if "fsum" in info and info["fsum"] < 0.9 and r["axis"] != "avg":
            fsum_low.setdefault(round(info["fsum"], 3), []).append(r["axis"])
    for value, axes_low in sorted(fsum_low.items()):
        n_bands_used = max((i.get("nbands", 0) for i in siesta_infos), default=0)
        print_dual(warn(f"SIESTA's f-sum check is {value:.3f} for {', '.join(axes_low)} (1 = complete): the transition "
                         "sum does not account for all valence electrons. More Optical.NumberOfBands helps only if "
                         f"fewer than all bands were used ({n_bands_used} here); beyond that the limit is the basis "
                         "set -- a richer basis (TZP/TZDP, smaller PAO.EnergyShift) is what raises it."), f_out)
    print_table(["Dir", "SCF ok", "f-sum (SIESTA)", "N_eff(E_max)", "N_eff/N_val", "E range (eV)", "Points",
                 "Step (eV)"], quality_rows, f_out)
    print_dual(color_text(
        "N_eff(E) = (2 eps0 m V / pi hbar^2) Int E eps2 dE is the f-sum rule evaluated on the spectrum: the "
        "number of electrons whose transitions are included up to E. It should approach the valence "
        "electron count (SIESTA's own check, f-sum, is that ratio). A value well below 1 means too few "
        "bands/too small an energy range, i.e. a truncated spectrum.", 'cyan'), f_out)
    if args.energy_max is not None:
        print_dual(color_text(
            f"[NOTE] --energy-max {args.energy_max:g} eV trims the data used here: N_eff and the Kramers-Kronig "
            "eps1 are evaluated on the trimmed spectrum, so they are lower/less accurate than on the full range.", 'cyan'), f_out)

    # ================================================================ static response & edge
    section('STATIC RESPONSE & ABSORPTION EDGE')
    stats_by_axis = {}
    static_rows = []
    for r in plotted_results:
        st = static_properties(r, volume_density, n_valence if n_valence_known else None)
        stats_by_axis[r["axis"]] = st
        static_rows.append(([r["axis"], f"{st['eps1_0']:.4f}", f"{st['n0']:.4f}", f"{st['R0']:.4f}",
                              f"{st['onset']:.2f}" if np.isfinite(st['onset']) else "-",
                              f"{st['tauc']:.2f}" if (np.isfinite(st['tauc']) and not metallic) else "-",
                              f"{st['penn']:.2f}" if np.isfinite(st['penn']) else "-"], None))
    print_table(["Dir", "eps1(E->0)", "n(E->0)", "R(E->0)", "eps2 onset (eV)", "Tauc direct gap (eV)",
                 "Penn gap (eV)"], static_rows, f_out)
    ref_stats = stats_by_axis.get(ref["axis"]) or static_properties(ref, volume_density, n_valence if n_valence_known else None)
    if metallic:
        ref_stats = dict(ref_stats, tauc=float('nan'), tauc_fit=None, tauc_window=None)   # a Tauc gap means nothing for a metal
    print_dual(color_text(
        "eps2 onset: first energy where the smoothed eps2 exceeds 5 % of its maximum (the independent-particle "
        "direct absorption edge; no indirect transitions exist in this spectrum). Tauc: extrapolation of "
        "(alpha E)^2 vs E on the rising edge -- indicative, it assumes parabolic bands. Penn gap: "
        "hbar*omega_p / sqrt(eps1(0) - 1), the one-oscillator gap that would give this static constant.", 'cyan'), f_out)
    if gap_info:
        kind = "direct" if gap_info["direct"] else "indirect"
        if metallic:
            print_dual(f"Kohn-Sham band structure (SCF k-mesh, {gap_info['n_k']} points): metallic/semimetallic -- no gap "
                        "(bands reach the Fermi level).", f_out)
        else:
            print_dual(f"Kohn-Sham band gap (SCF k-mesh, {gap_info['n_k']} points): fundamental {gap_info['gap']:.3f} eV "
                        f"({kind}); smallest direct gap {gap_info['direct_min']:.3f} eV.", f_out)
    if metallic:
        print_dual(color_text(
            "[NOTE] Zero Kohn-Sham gap: this is a metal/semimetal. The eps2 'onset' above is NOT a gap: it is set by the "
            "smallest transition energy the optical k-mesh samples and by the broadening (a Dirac point, for instance, "
            "needs a very dense mesh), and SIESTA's optical module has no intraband (Drude) term, so the low-energy "
            "response is incomplete.", 'cyan'), f_out)
    elif np.isfinite(ref_stats["onset"]):
        print_dual(f"Reference spectrum '{ref['axis']}': the absorption edge starts at {ref_stats['onset']:.2f} eV "
                    f"(wavelength {HC_EV_NM / ref_stats['onset']:.0f} nm)."
                    + (f" It matches the smallest direct Kohn-Sham gap ({gap_info['direct_min']:.2f} eV): the edge is direct"
                       + (", and the indirect gap is invisible here." if not gap_info['direct'] else ".")
                       if gap_info and abs(ref_stats['onset'] - gap_info['direct_min']) < 0.3 else ""), f_out)
        if ref_stats["onset"] > 0.05:
            print_dual(f"No absorption below the edge: the spectrum describes a dielectric (non-metallic) response; "
                        f"eps1(0) = {ref_stats['eps1_0']:.2f}.", f_out)

    # ================================================================ spectral features
    feature_emax = min(float(ref["omega"].max()), 30.0)
    features = spectral_features(ref, args.n_peaks, args.smoothing, feature_emax)
    section(f"SPECTRAL FEATURES (reference spectrum: {ref['axis']})")
    if features["eps2_peaks"]:
        peak_rows = []
        for j, pk in enumerate(features["eps2_peaks"], 1):
            e = pk["E"]
            peak_rows.append(([str(j), f"{e:.2f}", f"{pk['value']:.2f}",
                               f"{pk['fwhm']:.2f}" if np.isfinite(pk['fwhm']) else "-",
                               f"{interp_at(ref, 'eps1', e):+.2f}", f"{interp_at(ref, 'n', e):.3f}",
                               f"{interp_at(ref, 'k', e):.3f}", f"{interp_at(ref, 'alpha', e):.3e}",
                               f"{interp_at(ref, 'R', e):.3f}"], None))
        print_dual(f"Peaks of eps2 up to {feature_emax:g} eV (smoothed with sigma = {args.smoothing:g} eV; "
                    "interband critical-point-like structure):", f_out)
        print_table(["#", "E (eV)", "eps2", "FWHM (eV)", "eps1", "n", "k", "alpha (cm^-1)", "R"], peak_rows, f_out)
    else:
        print_dual("No eps2 peak above the prominence threshold in this range.", f_out)
    print_dual("", f_out)
    feat_rows = []
    for key, label, fmt in (("alpha_max", "Absorption coefficient maximum", "{:.4e} cm^-1"),
                            ("R_max", "Reflectivity maximum", "{:.4f}"),
                            ("eps1_max", "Maximum of eps1", "{:.3f}")):
        e, v = features[key]
        feat_rows.append(([label, f"{e:.2f}", fmt.format(v)], None))
    lp = features["loss_peak"]
    if lp:
        feat_rows.append((["Energy-loss peak (bulk plasmon)", f"{lp['E']:.2f}",
                            f"L = {lp['value']:.2f}, FWHM {lp['fwhm']:.2f} eV" if np.isfinite(lp['fwhm'])
                            else f"L = {lp['value']:.2f}"], None))
    print_table(["Feature", "E (eV)", "Value"], feat_rows, f_out)
    if features["eps1_zeros"]:
        print_dual("Zero crossings of eps1: " + ", ".join(
            f"{e:.2f} eV ({'rising' if d > 0 else 'falling'})" for e, d in features["eps1_zeros"][:8]), f_out)
    if is_vacuum_cell:
        print_dual(color_text("[NOTE] Vacuum-padded cell: the free-electron plasmon comparison and the Penn gap are skipped (the cell "
                               "volume includes the vacuum, so the electron density is not the material's).", 'cyan'), f_out)
    if lp and np.isfinite(ref_stats["plasmon_free"]):
        print_dual(f"Free-electron plasma energy for {n_valence:g} valence electrons in this cell: "
                    f"{ref_stats['plasmon_free']:.2f} eV; the loss peak sits at {lp['E']:.2f} eV "
                    f"({lp['E'] / ref_stats['plasmon_free']:.3f} x). A ratio near 1 marks a collective "
                    "valence-electron plasmon.", f_out)

    # ================================================================ anisotropy & symmetry
    section('ANISOTROPY & SYMMETRY CHECKS')
    verdicts = []
    if len(diag_results) >= 2:
        eps_static = [stats_by_axis[r["axis"]]["eps1_0"] for r in diag_results]
        n_static = [stats_by_axis[r["axis"]]["n0"] for r in diag_results]
        mean2 = np.mean([r["eps2"] for r in diag_results], axis=0)
        norm_mean = float(np.linalg.norm(mean2)) or 1.0
        dev_rows = []
        for r in diag_results:
            dev = float(np.linalg.norm(r["eps2"] - mean2) / norm_mean)
            j = int(np.argmax(r["eps2"]))
            dev_rows.append(([r["axis"], f"{stats_by_axis[r['axis']]['eps1_0']:.4f}",
                               f"{stats_by_axis[r['axis']]['n0']:.4f}", f"{r['omega'][j]:.2f}",
                               f"{r['eps2'][j]:.2f}", f"{dev:.2e}"], None))
        print_table(["Dir", "eps1(0)", "n(0)", "E of max eps2 (eV)", "max eps2", "||eps2 - mean|| / ||mean||"],
                    dev_rows, f_out)
        spread = (max(eps_static) - min(eps_static)) / (np.mean(eps_static) or 1.0)
        print_dual(f"Static anisotropy (eps_max - eps_min)/mean : {spread:.3e};  birefringence "
                    f"dn(0) = {max(n_static) - min(n_static):.4f}", f_out)
        verdicts.append(("diag", spread, max(float(np.linalg.norm(r['eps2'] - mean2) / norm_mean) for r in diag_results)))
    else:
        print_dual("Fewer than two diagonal directions: no anisotropy to compare.", f_out)
    if sym is not None and len(diag_results) >= 2:
        locked = [c for cls in sym["classes"] for c in cls if len(cls) > 1 and c in _DIAGONAL_AXES]
        have = [r for r in diag_results if r["axis"] in locked]
        if len(have) >= 2:
            ref_c = have[0]
            worst = max(float(np.linalg.norm(r["eps2"] - ref_c["eps2"]) / (np.linalg.norm(ref_c["eps2"]) or 1.0))
                        for r in have[1:])
            ok = worst < 1e-2
            print_dual(f"Symmetry expects {' = '.join('eps_' + c for c in locked)}: observed relative deviation "
                        f"{worst:.2e} -> " + (color_text("consistent", 'green') if ok else
                                               warn(f"directions that symmetry makes equal differ by {worst:.1e}: check that "
                                                    "all folders used the same settings/structure")), f_out)
    print_dual("", f_out)

    offdiag_results = reconstruct_offdiagonal(results)
    section('OFF-DIAGONAL DIELECTRIC TENSOR (eps_ij RECONSTRUCTION)')
    if offdiag_results:
        print_dual("eps_ij = eps_(bisector ij) - (eps_ii + eps_jj) / 2, applied to eps1(E) and "
                    "eps2(E) independently -- see stb-optical --help for the full derivation.", f_out)
        scale = max((float(np.max(np.abs(r["eps2"]))) for r in diag_results), default=1.0) or 1.0
        off_rows = []
        for r in offdiag_results:
            i0 = int(np.argmin(r["omega"]))
            need_a, need_b = r["needs"]
            jmax = int(np.argmax(np.abs(r["eps2"])))
            rel = float(np.max(np.abs(r["eps2"])) / scale)
            off_rows.append(([f"eps_{r['axis']}", f"dir_{r['axis']}, dir_{need_a}, dir_{need_b}",
                               f"{r['eps1'][i0]:+.4f}", f"{r['eps2'][i0]:+.4f}",
                               f"{r['eps2'][jmax]:+.4f} @ {r['omega'][jmax]:.2f} eV", f"{rel:.2e}"], None))
            if r["interpolated"]:
                print_dual(color_text(
                    f"[NOTE] Direction energy grids did not match exactly -- interpolated "
                    f"onto {r['axis']}'s own grid before reconstructing.", 'yellow'), f_out)
            if not r["scf_ok"]:
                print_dual(warn(f"SCF convergence not confirmed for one of the 3 folders feeding eps_{r['axis']}: "
                                 "this component may be unreliable."), f_out)
        print_table(["Component", "From", "eps1(E->0)", "eps2(E->0)", "largest |eps2_ij|", "relative to max eps2_ii"],
                    off_rows, f_out)
        if sym is not None:
            zero_by_sym = [c for c in sym["zero"] if c in {r["axis"] for r in offdiag_results}]
            if zero_by_sym:
                worst_rel = max(float(np.max(np.abs(r["eps2"])) / scale) for r in offdiag_results if r["axis"] in zero_by_sym)
                ok = worst_rel < 1e-2
                print_dual(f"Symmetry makes {', '.join('eps_' + c for c in zero_by_sym)} vanish: largest observed "
                            f"{worst_rel:.2e} of max eps2 -> " + (color_text("consistent", 'green') if ok else
                            warn("an off-diagonal component that symmetry forbids is not ~0: check the folders")), f_out)
        n_offdiag_possible = len(OPTICAL_OFFDIAG_PAIRS)
        if len(offdiag_results) < n_offdiag_possible:
            missing_pairs = [axis for axis in OPTICAL_OFFDIAG_PAIRS if axis not in
                              {r["axis"] for r in offdiag_results}]
            print_dual(f"({', '.join(missing_pairs)} not reconstructed -- missing the biaxial "
                        "folder and/or one of its 2 matching diagonal folders.)", f_out)
    else:
        print_dual("Not available: it needs a biaxial folder (xy/xz/yz) plus the two matching diagonal folders.", f_out)

    # ================================================================ selected energies
    section(f"OPTICAL CONSTANTS AT SELECTED ENERGIES (reference spectrum: {ref['axis']})")
    sel_rows = []
    for e in args.energies:
        e_eval = max(e, float(ref["omega"].min()))
        sel_rows.append(([f"{e:g}" if e > 0 else "0 (static)", f"{HC_EV_NM / e_eval:.0f}" if e_eval > 0 else "-",
                           f"{interp_at(ref, 'eps1', e_eval):+.3f}", f"{interp_at(ref, 'eps2', e_eval):.3f}",
                           f"{interp_at(ref, 'n', e_eval):.4f}", f"{interp_at(ref, 'k', e_eval):.4f}",
                           f"{interp_at(ref, 'alpha', e_eval):.3e}", f"{interp_at(ref, 'R', e_eval):.4f}",
                           f"{interp_at(ref, 'L', e_eval):.3e}"], None))
    print_table(["E (eV)", "lambda (nm)", "eps1", "eps2", "n", "k", "alpha (cm^-1)", "R", "L"], sel_rows, f_out)

    # ================================================================ dimensionality correction
    corrected_dat_paths = []
    polarizability_rows = []
    if args.dimensionality_correction:
        section('DIMENSIONALITY CORRECTION')
        if structure is None or vacuum_axes is None:
            print_dual(color_text(
                "  [WARNING] Could not read structure.fdf -- --dimensionality-correction "
                "skipped.", 'yellow'), f_out)
        else:
            n_vacuum = sum(vacuum_axes)
            lattice = np.array(structure.lattice, dtype=float)
            cell_volume_ang3 = float(abs(np.dot(lattice[0], np.cross(lattice[1], lattice[2]))))

            if n_vacuum == 0:
                print_dual("  3D (bulk) input -- nothing to correct.", f_out)
            elif n_vacuum == 2:
                print_dual(color_text(
                    "  [NOTE] 1D (wire/tube-like) input -- dimensionality correction is NOT "
                    "implemented in this version. The literature approach (Maxwell-Garnett "
                    "effective-medium theory with cross-sectional volume fraction) is not a "
                    "simple generalization of the 2D formula and has not been verified for "
                    "this tool. Reporting only the raw supercell values above.", 'yellow'), f_out)
            elif n_vacuum == 1:
                if args.thickness is None:
                    print_dual(color_text(
                        "  [ERROR] --dimensionality-correction on a 2D (1 vacuum axis) input "
                        "requires --thickness (Ang).", 'red'), f_out)
                    if f_out:
                        f_out.close()
                    sys.exit(1)
                vacuum_idx = vacuum_axes.index(True)
                cell_length_ang = float(np.linalg.norm(lattice[vacuum_idx]))
                if args.thickness >= cell_length_ang:
                    print_dual(color_text(
                        f"  [ERROR] --thickness ({args.thickness} Ang) must be smaller than "
                        f"the supercell's own length along the vacuum axis ({cell_length_ang:.4f} "
                        "Ang) -- the material can't be thicker than the whole cell.", 'red'), f_out)
                    if f_out:
                        f_out.close()
                    sys.exit(1)
                print_dual(f"  Reference : Laturia, Van de Put, Vandenberghe, npj 2D Mater. "
                            "Appl. 2, 6 (2018); Yang & Gao, npj 2D Mater. Appl. 5, 78 (2021).", f_out)
                print_dual(f"  Cell length (vacuum axis) : {cell_length_ang:.4f} Ang", f_out)
                print_dual(f"  Thickness (user-supplied) : {args.thickness} Ang", f_out)
                corrected_results = []
                for result in results:
                    if result["axis"] not in _DIAGONAL_AXES:
                        continue
                    axis_idx = {"xx": 0, "yy": 1, "zz": 2}[result["axis"]]
                    if axis_idx == vacuum_idx:
                        e1c, e2c = correct_2d_perpendicular(
                            result["eps1"], result["eps2"], cell_length_ang, args.thickness)
                        kind = "perpendicular"
                    else:
                        e1c, e2c = correct_2d_parallel(
                            result["eps1"], result["eps2"], cell_length_ang, args.thickness)
                        kind = "parallel"
                    corrected = derive_from_eps(result["omega"], e1c, e2c)
                    corrected["axis"] = result["axis"]
                    print_dual(f"  Direction {result['axis']} ({kind}) -- eps1_2D(E->0) = "
                                f"{corrected['eps1'][int(np.argmin(corrected['omega']))]:.4f} "
                                f"(raw supercell: {result['eps1'][int(np.argmin(result['omega']))]:.4f})",
                                f_out)
                    corrected_results.append(corrected)
                if corrected_results:
                    corrected_csv_path = os.path.join(args.directory, f"{args.output}_2Dcorrected.csv")
                    write_results_csv(corrected_csv_path, corrected_results)
                    for corrected in corrected_results:
                        dat_path = os.path.join(
                            args.directory, f"{args.output}_{corrected['axis']}_2Dcorrected.dat")
                        write_direction_dat(dat_path, corrected)
                        corrected_dat_paths.append(dat_path)
                    corrected_dat_paths.append(corrected_csv_path)
            elif n_vacuum == 3:
                print_dual(f"  Cell volume : {cell_volume_ang3:.4f} Ang^3", f_out)
                print_dual("  Extracting molecular polarizability (dilute Clausius-Mossotti "
                            "relation).", f_out)
                for result in results:
                    if result["axis"] not in _DIAGONAL_AXES:
                        continue
                    a1si, a2si, a1a3, a2a3 = molecular_polarizability(
                        result["eps1"], result["eps2"], cell_volume_ang3)
                    i0 = int(np.argmin(result["omega"]))
                    print_dual(f"  Direction {result['axis']} -- alpha(E->0) = {a1a3[i0]:.4f} "
                                "Ang^3", f_out)
                    dat_path = os.path.join(
                        args.directory, f"{args.output}_{result['axis']}_polarizability.dat")
                    write_polarizability_dat(dat_path, result["axis"], result["omega"],
                                              a1si, a2si, a1a3, a2a3)
                    corrected_dat_paths.append(dat_path)
                    polarizability_rows.append(
                        (result["axis"], result["omega"], a1si, a2si, a1a3, a2a3))
                if polarizability_rows:
                    pol_csv_path = os.path.join(args.directory, f"{args.output}_polarizability.csv")
                    write_polarizability_csv(pol_csv_path, polarizability_rows)
                    corrected_dat_paths.append(pol_csv_path)


    # ================================================================ scissor what-if
    scissor_result = None
    if args.scissor is not None:
        section(f"SCISSOR WHAT-IF ({args.scissor:+g} eV, post-processing of the {ref['axis']} spectrum)")
        eps2_sc = apply_scissor(ref["omega"], ref["eps2"], args.scissor)
        scissor_result = compute_all(ref["omega"], eps2_sc)
        scissor_result["axis"], scissor_result["delta"] = ref["axis"], args.scissor
        st_sc = static_properties(scissor_result, volume_density, n_valence if n_valence_known else None)
        feat_sc = spectral_features(scissor_result, args.n_peaks, args.smoothing, feature_emax)
        e2_pk = lambda f: (f"{f['eps2_peaks'][0]['E']:.2f}" if f["eps2_peaks"] else "-")
        main_pk = lambda f: (f"{max(f['eps2_peaks'], key=lambda p: p['value'])['E']:.2f}" if f["eps2_peaks"] else "-")
        sc_rows = [
            (["eps1(E->0)", f"{ref_stats['eps1_0']:.3f}", f"{st_sc['eps1_0']:.3f}"], None),
            (["n(E->0)", f"{ref_stats['n0']:.3f}", f"{st_sc['n0']:.3f}"], None),
            (["R(E->0)", f"{ref_stats['R0']:.4f}", f"{st_sc['R0']:.4f}"], None),
            (["eps2 onset (eV)", f"{ref_stats['onset']:.2f}", f"{st_sc['onset']:.2f}"], None),
            (["first eps2 peak (eV)", e2_pk(features), e2_pk(feat_sc)], None),
            (["strongest eps2 peak (eV)", main_pk(features), main_pk(feat_sc)], None),
            (["loss peak (eV)", f"{features['loss_peak']['E']:.2f}" if features['loss_peak'] else "-",
              f"{feat_sc['loss_peak']['E']:.2f}" if feat_sc['loss_peak'] else "-"], None),
        ]
        print_table(["Quantity", "DFT", f"scissor {args.scissor:+g} eV"], sc_rows, f_out)
        print_dual(color_text(
            "eps2'(E) = eps2(E - delta) (zero below delta), the position-matrix-element-conserving scissor "
            "(Levine and Allan, Phys. Rev. Lett. 63, 1719 (1989)); eps1 is re-derived by Kramers-Kronig. It "
            "corrects the gap error of semilocal DFT as a rigid shift, not the excitonic/quasiparticle "
            "physics, so use it to see where the peaks land, not as a prediction.", 'cyan'), f_out)

    # ================================================================ experimental comparison
    experimental, experimental_dat_path = None, None
    if args.experimental is not None:
        exp_energy, exp_value = read_experimental_spectrum(args.experimental)
        experimental = (exp_energy, exp_value)
        experimental_dat_path = os.path.join(args.directory, f"{args.output}_experimental.dat")
        section('EXPERIMENTAL COMPARISON')
        ref_result = next((r for r in results if r["axis"] == "avg"), results[0])
        print_dual(f"Simulated reference : direction '{ref_result['axis']}'", f_out)
        sim_key = _QUANTITY_COLUMNS[args.plot_quantity][0][0]
        sim_peaks = find_spectrum_peaks(ref_result["omega"], ref_result[sim_key], args.peak_prominence)
        exp_peaks = find_spectrum_peaks(exp_energy, exp_value, args.peak_prominence)
        matches = match_peaks(sim_peaks, exp_peaks)
        if matches:
            for exp_e, sim_e, delta in matches:
                if sim_e is None:
                    print_dual(f"  Exp {exp_e:.4f} eV -> no simulated peak to match", f_out)
                else:
                    print_dual(f"  Exp {exp_e:.4f} eV -> Sim {sim_e:.4f} eV (delta = {delta:+.4f} eV)", f_out)
            deltas = [d for _e, s, d in matches if s is not None]
            if deltas:
                print_dual(f"  Mean shift (sim - exp): {np.mean([-d for d in deltas]):+.3f} eV over {len(deltas)} matched peak(s)"
                            " -- a rigid shift of this size would be a scissor candidate (--scissor).", f_out)
        else:
            print_dual("  No experimental peaks found above --peak-prominence threshold.", f_out)

    # ================================================================ interpretation
    section('INTERPRETATION')
    bullets = []
    st = ref_stats
    if metallic:
        bullets.append("Metallic/semimetallic response (zero Kohn-Sham gap): interband absorption starts at the lowest "
                       "transition energy the k-mesh samples, not at a gap; the free-carrier (Drude) response is absent "
                       f"from SIESTA's optical module. eps1(0) = {st['eps1_0']:.2f} is the interband contribution only.")
    elif np.isfinite(st["onset"]) and st["onset"] > 0.05:
        bullets.append(f"Dielectric (gapped) response: no absorption below {st['onset']:.2f} eV "
                       f"({HC_EV_NM / st['onset']:.0f} nm); static dielectric constant eps1(0) = {st['eps1_0']:.2f}, "
                       f"refractive index n(0) = {st['n0']:.2f}, normal-incidence reflectivity R(0) = {st['R0']:.2f}.")
    elif np.isfinite(st["onset"]):
        bullets.append("Absorption starts at ~0 eV: a zero-gap response (no Drude term in SIESTA's optical module).")
    if features["eps2_peaks"]:
        cands = [p for p in features["eps2_peaks"] if not metallic or p["E"] > 1.0] or features["eps2_peaks"]
        main = max(cands, key=lambda p: p["value"])
        bullets.append(f"Strongest interband structure: eps2 = {main['value']:.1f} at {main['E']:.2f} eV "
                       f"(absorption coefficient up to {features['alpha_max'][1]:.2e} cm^-1 at {features['alpha_max'][0]:.2f} eV)."
                       + (" Peaks below 1 eV are set by the k-mesh sampling of the Dirac/Fermi region and are not listed here." if metallic else ""))
    if features["loss_peak"]:
        lp = features["loss_peak"]
        extra = ""
        if np.isfinite(st["plasmon_free"]):
            extra = f", {lp['E'] / st['plasmon_free']:.2f} x the free-electron plasma energy of {st['plasmon_free']:.1f} eV"
        bullets.append(f"Bulk plasmon (energy-loss maximum) at {lp['E']:.2f} eV{extra}; eps1 crosses zero there.")
    if np.isfinite(st["penn"]):
        bullets.append(f"Penn gap {st['penn']:.2f} eV: the single average gap consistent with eps1(0) = {st['eps1_0']:.2f}.")
    fsums = [i["fsum"] for i in siesta_infos if "fsum" in i]
    if fsums:
        worst = min(fsums)
        bullets.append("Sum rule: SIESTA's f-sum check is " + f"{worst:.3f}" + (
            " -- complete." if worst >= 0.95 else " -- the transition sum is incomplete (limited by the number of bands "
            "and, once all are used, by the basis set): more weight at high energy would lower eps1(0)."))
    if verdicts:
        _kind, spread, dev = verdicts[0]
        bullets.append("The three diagonal directions " + ("agree (isotropic response)" if spread < 1e-2 and dev < 1e-2 else
                       f"differ (static anisotropy {spread:.1e}, spectral deviation {dev:.1e})") + ".")
    bullets.append("This is an independent-particle (RPA, no local fields, no excitons) spectrum on top of a semilocal "
                   "ground state: expect the whole spectrum red-shifted by about the DFT gap error and eps1(0) "
                   "overestimated; --scissor shows where the peaks would land.")
    for b in bullets:
        print_dual(f"* {b}", f_out)

    # ================================================================ summary & files
    section('SUMMARY & FILES')

    csv_path = os.path.join(args.directory, f"{args.output}.csv")
    write_results_csv(csv_path, results)
    print_dual(f"CSV (all directions) : {csv_path}", f_out)

    summary_csv = os.path.join(args.directory, f"{args.output}_summary.csv")
    with open(summary_csv, "w") as fh:
        fh.write("direction,eps1_static,n_static,R_static,onset_eV,tauc_direct_gap_eV,penn_gap_eV,N_eff_Emax,fsum_siesta,ks_gap_eV,ks_direct_min_eV\n")
        for r in plotted_results:
            s_r = stats_by_axis[r["axis"]]
            neff_end = neff_by_axis[r["axis"]][-1] if r["axis"] in neff_by_axis else float("nan")
            fsum = (r.get("siesta") or {}).get("fsum", float("nan"))
            fh.write(f"{r['axis']},{s_r['eps1_0']:.6f},{s_r['n0']:.6f},{s_r['R0']:.6f},{s_r['onset']:.4f},"
                     f"{s_r['tauc']:.4f},{s_r['penn']:.4f},{neff_end:.4f},{fsum},"
                     f"{gap_info.get('gap', float('nan')):.4f},{gap_info.get('direct_min', float('nan')):.4f}\n")
    print_dual(f"CSV (summary)        : {summary_csv}", f_out)

    dat_paths_by_axis = {}
    for result in results:
        dat_path = os.path.join(args.directory, f"{args.output}_{result['axis']}.dat")
        write_direction_dat(dat_path, result)
        dat_paths_by_axis[result["axis"]] = dat_path
        print_dual(f"Direction {result['axis']} .dat     : {dat_path}", f_out)

    if scissor_result is not None:
        sc_dat = os.path.join(args.directory, f"{args.output}_{ref['axis']}_scissor.dat")
        write_direction_dat(sc_dat, dict(scissor_result, axis=f"{ref['axis']}_scissor{args.scissor:+g}eV"))
        print_dual(f"Scissor what-if .dat : {sc_dat}", f_out)

    if offdiag_results:
        offdiag_csv_path = os.path.join(args.directory, f"{args.output}_offdiagonal.csv")
        write_offdiagonal_csv(offdiag_csv_path, offdiag_results)
        print_dual(f"CSV (off-diagonal)   : {offdiag_csv_path}", f_out)
        offdiag_dat_paths = {}
        for r in offdiag_results:
            dat_path = os.path.join(args.directory, f"{args.output}_{r['axis']}_offdiag.dat")
            write_offdiagonal_dat(dat_path, r)
            offdiag_dat_paths[r["axis"]] = dat_path
            print_dual(f"eps_{r['axis']} .dat            : {dat_path}", f_out)
        offdiag_gplot_path = os.path.join(args.directory, f"{args.output}_offdiagonal.gplot")
        write_offdiagonal_gplot(offdiag_gplot_path, offdiag_dat_paths)
        print_dual(f"Combined plot (off-diagonal) : {offdiag_gplot_path} "
                    f"(cd {args.directory} && gnuplot {os.path.basename(offdiag_gplot_path)})", f_out)

    for path in corrected_dat_paths:
        print_dual(f"Dimensionality correction : {path}", f_out)

    if args.plot_quantity == "all":
        for quantity in ("eps", "n_k", "alpha", "R", "L", "sigma1"):
            gplot_path = os.path.join(args.directory, f"{args.output}_{quantity}.gplot")
            write_combined_gplot(gplot_path, dat_paths_by_axis, quantity)
            print_dual(f"Combined plot ({quantity}) : {gplot_path}", f_out)
    else:
        gplot_path = os.path.join(args.directory, f"{args.output}.gplot")
        write_combined_gplot(gplot_path, dat_paths_by_axis, args.plot_quantity,
                              experimental, experimental_dat_path)
        print_dual(f"Combined plot ({args.plot_quantity}) : {gplot_path} "
                    f"(cd {args.directory} && gnuplot {os.path.basename(gplot_path)})", f_out)

    if args.png or args.show:
        plot_emax = args.plot_emax if args.plot_emax else min(float(ref["omega"].max()), 20.0)
        png_path = os.path.join(args.directory, f"{args.output}_summary.png") if args.png else None
        if args.show:
            print_dual("Dashboard            : opening on screen (close the window to continue)...", f_out)
        label = os.path.basename(os.path.abspath(args.directory))
        done = make_dashboard(results, ref, ref_stats | {"loss_peak": features["loss_peak"]}, offdiag_results,
                               scissor_result, neff_by_axis.get(ref["axis"]), n_valence if n_valence_known else None,
                               plot_emax, png_path, args.show, experimental,
                               args.plot_quantity if args.experimental else None,
                               title=f"Optical properties -- {label} (reference: {ref['axis']})")
        if done and png_path:
            print_dual(f"Dashboard (PNG)      : {png_path}", f_out)
    else:
        print_dual("Dashboard            : not requested (add --show to open it, --png to save it)", f_out)
    if report_path:
        print_dual(f"Report               : {report_path}", f_out)

    if warnings:
        print_dual(color_text(f"\n{len(warnings)} warning(s) above:", 'yellow'), f_out)
        for w in warnings:
            print_dual(f"  - {w}", f_out)

    if f_out:
        f_out.close()

    print("\n[INFO] Complete job!")
    print("\n" + "-" * 60)
    print(color_text("Optical properties analysis complete!\n", 'bold'))


if __name__ == "__main__":
    main()
