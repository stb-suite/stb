#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.0.0"

import os
import re
import sys
import glob
import argparse
import numpy as np
from stb.core import structure_io, kspace
from stb.core.cli import color_text, show_intro, print_dual, print_section
from stb.core.dielectric import read_epsimg
from stb.core.calc_directives import OPTICAL_DIRECTION_VECTORS, OPTICAL_OFFDIAG_PAIRS
from stb.core.optical_properties import (
    compute_all, correct_2d_perpendicular, correct_2d_parallel, molecular_polarizability,
    derive_from_eps,
)
from stb.core.siesta_log import get_scf_convergence
from stb.core.spectrum import read_experimental_spectrum, find_spectrum_peaks, match_peaks

REPORT_FILE = "optical_stage2.txt"

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
experimental measurement, which doesn't correspond to any single crystallographic direction.""",
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

    print_dual(f"{color_text('===== OPTICAL STAGE 2 REPORT (ANALYSIS) =====', 'magenta')}", f_out)

    print_section('[0] RUN METADATA', f_out)
    print_dual(f"Directory       : {args.directory}", f_out)
    print_dual(f"Folders found   : {len(folders)}", f_out)

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

    avg_result = compute_isotropic_average(results)
    if avg_result is not None:
        if avg_result["interpolated"]:
            print_dual(color_text(
                "  [NOTE] Direction energy grids did not match exactly -- interpolated onto "
                "the xx direction's grid before averaging.", 'yellow'), f_out)
        results.append(avg_result)

    structure = None
    vacuum_axes = None
    struct_path = os.path.join(results[0]["folder"], "structure.fdf")
    if os.path.isfile(struct_path):
        structure = structure_io.read_fdf(struct_path)
        positions = np.array([pos for _, pos in structure.atoms])
        is_cartesian = structure.coord_format == 'cartesian'
        frac_coords = kspace.to_fractional(positions, structure.lattice, is_cartesian)
        vacuum_axes = kspace.detect_vacuum_axes(frac_coords, structure.lattice, 10.0)
        print_section('[1] DIMENSIONALITY', f_out)
        print_dual(f"Detected : {kspace.dimensionality_label(vacuum_axes)}", f_out)
        if any(vacuum_axes):
            print_dual(color_text(
                "[KNOWN LIMITATION] Vacuum-padded (2D/slab) input -- see stb-optical --help; "
                "every direction's dielectric response below is diluted by the vacuum region, "
                "no rescaling applied unless --dimensionality-correction is given.", 'yellow'), f_out)

    print_section('[2] RESULTS BY DIRECTION', f_out)
    for result in results:
        i0 = int(np.argmin(result["omega"]))
        i_alpha_max = int(np.argmax(result["alpha"]))
        i_R_max = int(np.argmax(result["R"]))
        print_dual(f"  Direction {result['axis']} ({result['folder']}):", f_out)
        print_dual(f"    eps1(E->0)      : {result['eps1'][i0]:.4f}", f_out)
        print_dual(f"    n(E->0)         : {result['n'][i0]:.4f}", f_out)
        print_dual(f"    alpha peak      : {result['alpha'][i_alpha_max]:.4e} cm^-1 "
                    f"at {result['omega'][i_alpha_max]:.4f} eV", f_out)
        print_dual(f"    R peak          : {result['R'][i_R_max]:.4f} "
                    f"at {result['omega'][i_R_max]:.4f} eV", f_out)
        if not result["scf_ok"]:
            print_dual(color_text(
                f"    [WARNING] Could not confirm SCF convergence for '{result['folder']}' -- "
                "this direction's spectrum may be unreliable.", 'yellow'), f_out)

    offdiag_results = reconstruct_offdiagonal(results)
    if offdiag_results:
        print_section('[2b] OFF-DIAGONAL DIELECTRIC TENSOR (eps_ij RECONSTRUCTION)', f_out)
        print_dual("eps_ij = eps_(bisector ij) - (eps_ii + eps_jj) / 2, applied to eps1(E) and "
                    "eps2(E) independently -- see stb-optical --help for the full derivation.", f_out)
        for r in offdiag_results:
            i0 = int(np.argmin(r["omega"]))
            need_a, need_b = r["needs"]
            print_dual(f"  eps_{r['axis']} (from dir_{r['axis']}, dir_{need_a}, dir_{need_b}):", f_out)
            print_dual(f"    eps1_{r['axis']}(E->0)  : {r['eps1'][i0]:.4f}", f_out)
            print_dual(f"    eps2_{r['axis']}(E->0)  : {r['eps2'][i0]:.4f}", f_out)
            if r["interpolated"]:
                print_dual(color_text(
                    f"    [NOTE] Direction energy grids did not match exactly -- interpolated "
                    f"onto {r['axis']}'s own grid before reconstructing.", 'yellow'), f_out)
            if not r["scf_ok"]:
                print_dual(color_text(
                    f"    [WARNING] Could not confirm SCF convergence for one of the 3 folders "
                    f"feeding eps_{r['axis']} -- this component may be unreliable.", 'yellow'), f_out)
        n_offdiag_possible = len(OPTICAL_OFFDIAG_PAIRS)
        if len(offdiag_results) < n_offdiag_possible:
            missing_pairs = [axis for axis in OPTICAL_OFFDIAG_PAIRS if axis not in
                              {r["axis"] for r in offdiag_results}]
            print_dual(f"  ({', '.join(missing_pairs)} not reconstructed -- missing the biaxial "
                        "folder and/or one of its 2 matching diagonal folders.)", f_out)

    corrected_dat_paths = []
    polarizability_rows = []
    if args.dimensionality_correction:
        print_section('[2c] DIMENSIONALITY CORRECTION', f_out)
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

    print_section('[3] SUMMARY & FILES', f_out)

    csv_path = os.path.join(args.directory, f"{args.output}.csv")
    write_results_csv(csv_path, results)
    print_dual(f"CSV (all directions) : {csv_path}", f_out)

    dat_paths_by_axis = {}
    for result in results:
        dat_path = os.path.join(args.directory, f"{args.output}_{result['axis']}.dat")
        write_direction_dat(dat_path, result)
        dat_paths_by_axis[result["axis"]] = dat_path
        print_dual(f"Direction {result['axis']} .dat     : {dat_path}", f_out)

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
        experimental = None
        experimental_dat_path = None
        if args.experimental is not None:
            exp_energy, exp_value = read_experimental_spectrum(args.experimental)
            experimental = (exp_energy, exp_value)
            experimental_dat_path = os.path.join(args.directory, f"{args.output}_experimental.dat")

            print_section('[3b] EXPERIMENTAL COMPARISON', f_out)
            ref_result = next((r for r in results if r["axis"] == "avg"), results[0])
            print_dual(f"Simulated reference : direction '{ref_result['axis']}'", f_out)
            sim_key = _QUANTITY_COLUMNS[args.plot_quantity][0][0]
            sim_peaks = find_spectrum_peaks(ref_result["omega"], ref_result[sim_key],
                                             args.peak_prominence)
            exp_peaks = find_spectrum_peaks(exp_energy, exp_value, args.peak_prominence)
            matches = match_peaks(sim_peaks, exp_peaks)
            if matches:
                for exp_e, sim_e, delta in matches:
                    if sim_e is None:
                        print_dual(f"  Exp {exp_e:.4f} eV -> no simulated peak to match", f_out)
                    else:
                        print_dual(f"  Exp {exp_e:.4f} eV -> Sim {sim_e:.4f} eV "
                                    f"(delta = {delta:+.4f} eV)", f_out)
            else:
                print_dual("  No experimental peaks found above --peak-prominence threshold.", f_out)

        gplot_path = os.path.join(args.directory, f"{args.output}.gplot")
        write_combined_gplot(gplot_path, dat_paths_by_axis, args.plot_quantity,
                              experimental, experimental_dat_path)
        print_dual(f"\nCombined plot ({args.plot_quantity}) : {gplot_path} "
                    f"(cd {args.directory} && gnuplot {os.path.basename(gplot_path)})", f_out)
    if report_path:
        print_dual(f"Report                : {report_path}", f_out)

    if f_out:
        f_out.close()

    print("\n[INFO] Complete job!")
    print("\n" + "-" * 60)
    print(color_text("Optical properties analysis complete!\n", 'bold'))


if __name__ == "__main__":
    main()
