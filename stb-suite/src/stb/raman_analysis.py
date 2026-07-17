#!/usr/bin/env python

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
import json
import argparse
from datetime import datetime
import numpy as np
from stb.core.cli import color_text, show_intro, print_dual
from stb.core.dielectric import read_epsimg, static_dielectric
from stb.core.siesta_log import check_scf_and_force
from stb.core.spectrum import (
    bose_einstein_weight, build_lorentzian_spectrum, read_experimental_spectrum,
    find_spectrum_peaks, match_peaks,
)

REPORT_FILE = "raman_stage3.txt"
# 1 THz = 33.35641 cm^-1 -- the conventional Raman-spectroscopy frequency
# unit, used for the report table and spectrum plot even though Phonopy's
# own native unit (and Stage 2's MODE_TABLE) is THz throughout the rest of
# this workflow.
THZ_TO_CM1 = 33.35641

_LABEL_RE = re.compile(r'^\s*SystemLabel\s+(\S+)', re.IGNORECASE | re.MULTILINE)


def detect_optical_label(optical_root):
    """SystemLabel from the first mode_*/calc.fdf under optical_disp/ --
    same regex-based approach as core.phonon_workflow.detect_system_label,
    reimplemented locally (not extracted to core) since it globs a
    different folder-naming pattern (mode_* vs. disp-*) with only one real
    consumer so far.
    """
    mode_dirs = sorted(glob.glob(os.path.join(optical_root, "mode_*")))
    if not mode_dirs:
        return None
    fdf_path = os.path.join(mode_dirs[0], "calc.fdf")
    try:
        with open(fdf_path) as f:
            match = _LABEL_RE.search(f.read())
    except OSError:
        return None
    return match.group(1) if match else None


def read_mode_table(report_path):
    """Parses Stage 2's # MODE_TABLE (label, mode_index, band_index,
    frequency_thz, sign, axis, dir, [derived_from]) into a list of dicts --
    same "read the persisted # XXX_TABLE, don't infer from folder names"
    convention as stb-stackingfaultAnalysis's GRID_TABLE reader. The 8th
    column (derived_from) is optional -- older Stage 2 reports (written
    before stb-ramanModes --skip-degenerate existed) only have 7 columns,
    and every row simply has no derived_from (None).
    """
    rows = []
    in_table = False
    with open(report_path) as f:
        for line in f:
            if line.startswith("# MODE_TABLE"):
                in_table = True
                continue
            if not in_table:
                continue
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 7:
                continue
            rows.append({
                "label": parts[0],
                "mode_index": int(parts[1]),
                "band_index": int(parts[2]),
                "frequency_thz": float(parts[3]),
                "sign": parts[4],
                "axis": parts[5],
                "dir": parts[6],
                "derived_from": (int(parts[7]) if len(parts) >= 8 and parts[7] != "-" else None),
            })
    return rows


def group_by_mode(rows):
    """{mode_index: {"frequency_thz": f, "folders": {(sign, axis): dir},
    "derived_from": int or None}}. A mode written by --skip-degenerate as a
    DERIVED row has no folders at all -- its Raman tensor is reused from
    "derived_from"'s already-computed result instead of read from disk.
    """
    modes = {}
    for row in rows:
        entry = modes.setdefault(row["mode_index"], {
            "frequency_thz": row["frequency_thz"], "folders": {}, "derived_from": None})
        if row["derived_from"] is not None:
            entry["derived_from"] = row["derived_from"]
        if row["sign"] != "DERIVED":
            entry["folders"][(row["sign"], row["axis"])] = row["dir"]
    return modes


def compute_mode_tensor(folders, label, out_file, f_out):
    """Reads <label>.EPSIMG from each of the (sign, axis) folders for one
    mode (6 for a diagonal-only run, 12 with stb-ramanModes --full-tensor),
    gets the static dielectric constant via core.dielectric.
    static_dielectric, and returns (static_eps, n_missing, scf_ok_all) --
    static_eps is a {(sign, axis): value} dict, ready for
    raman_tensor_full's central finite difference. Missing/unparseable
    folders are skipped and counted in n_missing.
    """
    static_eps = {}
    n_missing = 0
    scf_ok_all = True
    for (sign, axis), folder in folders.items():
        epsimg_path = os.path.join(folder, f"{label}.EPSIMG")
        out_path = os.path.join(folder, out_file)
        scf_ok, max_force = check_scf_and_force(out_path)
        if not scf_ok:
            scf_ok_all = False
            print_dual(color_text(
                f"    [WARNING] Could not confirm SCF convergence for {folder} -- "
                "this dielectric value may be unreliable.", 'yellow'), f_out)
        if not os.path.isfile(epsimg_path):
            print_dual(color_text(f"    [SKIP] {epsimg_path} not found.", 'yellow'), f_out)
            n_missing += 1
            continue
        try:
            omega, eps2 = read_epsimg(epsimg_path)
            if len(omega) == 0:
                raise ValueError("no data rows parsed")
            static_eps[(sign, axis)] = static_dielectric(omega, eps2)
        except Exception as e:
            print_dual(color_text(f"    [SKIP] Could not parse {epsimg_path}: {e}", 'yellow'), f_out)
            n_missing += 1
    return static_eps, n_missing, scf_ok_all


# Off-diagonal probe axis -> the diagonal pair it's built from, matching
# stb-ramanModes' _AXES_OFFDIAG face-diagonal directions (i+j)/sqrt(2).
_OFFDIAG_PAIRS = {"xy": ("x", "y"), "xz": ("x", "z"), "yz": ("y", "z")}

# Tensor component -> (row, col) into the 3x3 T0 sidecar matrix, and the
# unit probe vector for each of the 6 possible directions -- same set and
# naming as stb-ramanModes'/core.raman_symmetry's _AXES_DIAGONAL/
# _AXES_OFFDIAG/_PROBE_VECTORS, duplicated here (not imported) to keep
# this Stage 3 module's dependency on Stage 2's/core's internals limited
# to the persisted JSON sidecar contract, not their Python objects.
_INV_SQRT2 = 0.70710678118654752440
_TENSOR_AXIS_INDICES = {"x": (0, 0), "y": (1, 1), "z": (2, 2), "xy": (0, 1), "xz": (0, 2), "yz": (1, 2)}
_PROBE_VECTORS = {
    "x": (1.0, 0.0, 0.0), "y": (0.0, 1.0, 0.0), "z": (0.0, 0.0, 1.0),
    "xy": (_INV_SQRT2, _INV_SQRT2, 0.0), "xz": (_INV_SQRT2, 0.0, _INV_SQRT2),
    "yz": (0.0, _INV_SQRT2, _INV_SQRT2),
}


def reconstruct_reduced_tensor(sidecar_path, folders, label, out_file, delta_ang, f_out):
    """Reads a stb-ramanModes --reduce-components sidecar
    (mode_XX_tensor_form.json: the symmetry-fixed tensor SHAPE T0, and
    which single Optical.Vector direction was actually probed) plus that
    one +/-delta pair of folders, and reconstructs the FULL Raman tensor
    (Rxx, Ryy, Rzz, Rxy, Rxz, Ryz) algebraically: T0 combined with the
    ONE measured scale factor c = R_probe / (n_probe^T T0 n_probe) gives
    every component directly (Rij = c * T0_ij) -- no further finite
    differences needed, since T0's SHAPE (not just its scale) was already
    fixed by symmetry when stb-ramanModes decided this mode qualified for
    reduction (see core.raman_symmetry.tensor_form).

    Returns (r, scf_ok) -- r has the same {'x':...,'xy':...} dict shape
    raman_tensor_full() produces, so callers can treat a reduced mode
    identically to a directly-measured one from here on; scf_ok is
    whether both probe-axis folders confirmed SCF convergence. Returns
    (None, True) if the sidecar or its folders can't be read/parsed
    (caller reports this the same "incomplete" way as any other mode
    with missing/unparseable data).
    """
    try:
        with open(sidecar_path) as f:
            sidecar = json.load(f)
        t0 = np.array(sidecar["T0"])
        probe_axis = sidecar["probe_axis"]
    except (OSError, ValueError, KeyError) as e:
        print_dual(color_text(f"    [SKIP] Could not read {sidecar_path}: {e}", 'yellow'), f_out)
        return None, True

    plus_dir = folders.get(("plus", probe_axis))
    minus_dir = folders.get(("minus", probe_axis))
    if plus_dir is None or minus_dir is None:
        print_dual(color_text(
            f"    [SKIP] Missing +/-delta folder(s) at the reduced probe axis "
            f"'{probe_axis}'.", 'yellow'), f_out)
        return None, True

    scf_ok = True
    static_vals = {}
    for sign, folder in (("plus", plus_dir), ("minus", minus_dir)):
        epsimg_path = os.path.join(folder, f"{label}.EPSIMG")
        out_path = os.path.join(folder, out_file)
        ok, _max_force = check_scf_and_force(out_path)
        scf_ok = scf_ok and ok
        if not ok:
            print_dual(color_text(
                f"    [WARNING] Could not confirm SCF convergence for {folder} -- "
                "this dielectric value may be unreliable.", 'yellow'), f_out)
        if not os.path.isfile(epsimg_path):
            print_dual(color_text(f"    [SKIP] {epsimg_path} not found.", 'yellow'), f_out)
            return None, scf_ok
        try:
            omega, eps2 = read_epsimg(epsimg_path)
            if len(omega) == 0:
                raise ValueError("no data rows parsed")
            static_vals[sign] = static_dielectric(omega, eps2)
        except Exception as e:
            print_dual(color_text(f"    [SKIP] Could not parse {epsimg_path}: {e}", 'yellow'), f_out)
            return None, scf_ok

    r_probe = (static_vals["plus"] - static_vals["minus"]) / (2.0 * delta_ang)
    n = np.array(_PROBE_VECTORS[probe_axis])
    denom = float(n @ t0 @ n)
    if abs(denom) < 1e-12:
        print_dual(color_text(
            f"    [SKIP] Degenerate reduced-tensor calibration in '{sidecar_path}' "
            "(near-zero denominator).", 'yellow'), f_out)
        return None, scf_ok
    c = r_probe / denom

    r = {axis: c * t0[i, j] for axis, (i, j) in _TENSOR_AXIS_INDICES.items()}
    return r, scf_ok


def raman_tensor_full(static_eps, delta_ang):
    """Rxx/Ryy/Rzz (and, when stb-ramanModes was run with --full-tensor,
    Rxy/Rxz/Ryz too) via a central finite difference across +/-delta.

    The diagonal components are the direct finite difference of the
    static dielectric constant measured along x/y/z. The off-diagonal
    ones are recovered from the mixed-direction measurement along
    n=(i+j)/sqrt(2): R_n = n^T R n = (Rii+Rjj)/2 + Rij, so
    Rij = R_n - (Rii+Rjj)/2, using the diagonal components already
    computed above.

    Returns (r, has_offdiag) -- r is a dict with keys 'x','y','z' always,
    plus 'xy','xz','yz' only if that mixed-direction data is present in
    static_eps at all (has_offdiag=True); any individual component is
    None if its own +/-delta pair is incomplete.
    """
    r = {}
    for axis in ("x", "y", "z"):
        plus = static_eps.get(("plus", axis))
        minus = static_eps.get(("minus", axis))
        r[axis] = None if (plus is None or minus is None) else (plus - minus) / (2.0 * delta_ang)

    has_offdiag = any((sign, axis) in static_eps
                       for sign in ("plus", "minus") for axis in _OFFDIAG_PAIRS)
    if has_offdiag:
        for mixed_axis, (i, j) in _OFFDIAG_PAIRS.items():
            plus = static_eps.get(("plus", mixed_axis))
            minus = static_eps.get(("minus", mixed_axis))
            if plus is None or minus is None or r[i] is None or r[j] is None:
                r[mixed_axis] = None
                continue
            r_mixed = (plus - minus) / (2.0 * delta_ang)
            r[mixed_axis] = r_mixed - (r[i] + r[j]) / 2.0
    return r, has_offdiag


def raman_activity(r, has_offdiag):
    """Raman activity (45*a'^2 + 7*gamma'^2, the standard
    isotropic-averaging-over-orientations formula). Uses the FULL
    invariant (gamma'^2 including the off-diagonal Rxy/Rxz/Ryz terms) when
    `has_offdiag` and every off-diagonal component is present; otherwise
    falls back to the diagonal-only approximation (missing the
    6*(Rxy^2+Rxz^2+Ryz^2)/2 anisotropy contribution).

    Returns (activity, is_full, alpha_prime, gamma_prime2) -- activity is
    None (with alpha_prime/gamma_prime2 also None) if any diagonal
    component is missing (nothing usable at all); is_full says which
    gamma'^2 formula was actually used, so the caller can report per-mode
    whether the invariant is exact or approximate. alpha_prime/gamma_prime2
    are returned alongside activity (not just folded into it) so
    depolarization_ratio() below can reuse them without recomputing.
    """
    if any(r[a] is None for a in ("x", "y", "z")):
        return None, False, None, None
    rxx, ryy, rzz = r["x"], r["y"], r["z"]
    alpha_prime = (rxx + ryy + rzz) / 3.0
    diag_aniso = (rxx - ryy) ** 2 + (ryy - rzz) ** 2 + (rzz - rxx) ** 2

    is_full = has_offdiag and all(r.get(a) is not None for a in _OFFDIAG_PAIRS)
    if is_full:
        rxy, rxz, ryz = r["xy"], r["xz"], r["yz"]
        gamma_prime2 = 0.5 * diag_aniso + 3.0 * (rxy ** 2 + rxz ** 2 + ryz ** 2)
    else:
        gamma_prime2 = 0.5 * diag_aniso

    activity = 45.0 * alpha_prime ** 2 + 7.0 * gamma_prime2
    return activity, is_full, alpha_prime, gamma_prime2


def depolarization_ratio(alpha_prime, gamma_prime2):
    """Depolarization ratio rho = 3*gamma'^2 / (45*alpha'^2 + 4*gamma'^2)
    (Long, "The Raman Effect" -- linearly polarized incident light, 90
    degree scattering, isotropic/randomly-oriented sample -- the standard
    convention reported alongside Raman activity in the literature). NOT
    the same denominator as raman_activity()'s 45a'^2+7gamma'^2 -- that's
    the total scattered intensity invariant, this is a different
    combination that isolates the polarized/depolarized intensity ratio.

    Totally symmetric modes (whose alpha_prime != 0) give rho < 3/4
    ("polarized"); every other mode gives rho = 3/4 exactly
    ("depolarized") for an isotropic sample -- a classic diagnostic,
    reported alongside the number.

    Returns None if both alpha_prime and gamma_prime2 are zero (rho is
    undefined, 0/0 -- a mode with no Raman activity at all).
    """
    denom = 45.0 * alpha_prime ** 2 + 4.0 * gamma_prime2
    if denom == 0:
        return None
    return 3.0 * gamma_prime2 / denom


def write_spectrum_plot(dat_path, gplot_path, grid, intensity, all_full, temperature_k=None,
                         experimental=None, experimental_dat_path=None):
    """experimental (optional): (freq_cm1, intensity) arrays of a real
    measured spectrum to overlay -- written to `experimental_dat_path`
    (required together with `experimental`) and plotted against a
    SEPARATE right-hand y-axis (gnuplot y2), since simulated "activity"
    and real experimental intensity are on completely different,
    incomparable scales -- a shared single axis would make one curve
    invisible next to the other.
    """
    scope = "full tensor" if all_full else "diagonal tensor only, approximate"
    if temperature_k is not None:
        scope += f", {temperature_k:.0f} K Bose weighting"
    with open(dat_path, "w") as f:
        f.write(f"# stb-ramanAnalysis -- Raman spectrum ({scope})\n")
        f.write("# columns: 1=frequency(cm^-1) 2=intensity(arb. units)\n")
        for freq, inten in zip(grid, intensity):
            f.write(f"{freq:12.4f} {inten:14.6f}\n")

    if experimental is not None:
        exp_freq, exp_intensity = experimental
        with open(experimental_dat_path, "w") as f:
            f.write("# stb-ramanAnalysis -- user-supplied experimental Raman spectrum\n")
            f.write("# columns: 1=frequency(cm^-1) 2=intensity(arb. units, original scale)\n")
            for freq, inten in zip(exp_freq, exp_intensity):
                f.write(f"{freq:12.4f} {inten:14.6f}\n")

    pdf_name = os.path.splitext(os.path.basename(dat_path))[0] + ".pdf"
    with open(gplot_path, "w") as f:
        lines = [
            '# --- STB Plot Configuration ---\n',
            '# Generated by stb-ramanAnalysis\n',
            'set terminal pdfcairo enhanced color font "Arial,14" size 8,5\n',
            f'set output "{pdf_name}"\n\n',
            f'set title "Raman Spectrum ({scope})"\n',
            'set xlabel "Raman shift (cm^{-1})"\n',
            'set ylabel "Simulated intensity (arb. units)"\n',
            'set grid\n',
        ]
        if experimental is not None:
            lines += [
                'set y2label "Experimental intensity (arb. units)"\n',
                'set y2tics\n',
                'set key top right\n',
                f'plot "{os.path.basename(dat_path)}" using 1:2 with lines lw 2 lc rgb "#2255cc" '
                'axes x1y1 title "Simulated", \\\n'
                f'     "{os.path.basename(experimental_dat_path)}" using 1:2 with lines lw 2 '
                'lc rgb "#cc5522" axes x1y2 title "Experimental"\n',
            ]
        else:
            lines += [
                'unset key\n',
                f'plot "{os.path.basename(dat_path)}" using 1:2 with lines lw 2 lc rgb "#2255cc"\n',
            ]
        f.writelines(lines)


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Stage 3 of 3: analyzes a stb-ramanModes Optical-calculation "
        "sweep and computes the Raman spectrum.", 'bold')}
Reads each mode's folders (+/-delta x, y, z, and, if stb-ramanModes was run
with --full-tensor, also the (x+y)/sqrt(2), (x+z)/sqrt(2), (y+z)/sqrt(2)
mixed directions), extracts the imaginary part of the dielectric function
from SystemLabel.EPSIMG, obtains the static dielectric constant via a
Kramers-Kronig transform (core/dielectric.py), and takes a central finite
difference across +/-delta to get the Raman tensor per mode. Auto-detects
per mode whether the full symmetric tensor (Rxx, Ryy, Rzz, Rxy, Rxz, Ryz)
or just the diagonal (Rxx, Ryy, Rzz) is available from the folders Stage 2
actually wrote -- the full tensor gives the EXACT literature Raman-activity
invariant (45a'^2+7gamma'^2); diagonal-only gives an APPROXIMATION missing
the off-diagonal anisotropy contribution. Writes a Lorentzian-summed
spectrum (.dat/.gplot).""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage example:\n"
               "  %(prog)s --directory raman_study\n"
    )

    parser.add_argument("-dir", "--directory", type=str, default="raman_study",
                        help="Root directory written by stb-raman/stb-ramanModes (default: raman_study).")
    parser.add_argument("-l", "--label", type=str, default=None,
                        help="SystemLabel used in Stage 2's Optical calc.fdf (default: auto-detected).")
    parser.add_argument("--file", type=str, default="calc.out",
                        help="SIESTA output filename inside each mode_*/ folder (default: calc.out).")
    parser.add_argument("--linewidth", type=float, default=10.0,
                        help="Lorentzian linewidth (cm^-1) for the summed spectrum (default: 10.0).")
    parser.add_argument("--temperature", type=float, default=None, metavar="K",
                        help="Weight each mode's spectrum contribution by the Stokes "
                             "Bose-Einstein thermal population factor (n+1) at this temperature "
                             "(Kelvin). Default: no weighting (every mode contributes by its raw "
                             "activity alone, as before this option existed).")
    parser.add_argument("-o", "--output", type=str, default="raman_spectrum",
                        help="Base filename (no extension) for the spectrum .dat/.gplot pair "
                             "(default: raman_spectrum).")
    parser.add_argument("--experimental", type=str, default=None, metavar="FILE",
                        help="Overlay a measured Raman spectrum (2-column text file: frequency "
                             "cm^-1, intensity) on the simulated one, and report a nearest-"
                             "neighbor peak match (simulated vs. experimental peak positions, "
                             "one row per experimental peak). Plotted on a separate y-axis "
                             "(simulated 'activity' and real intensity aren't on the same scale).")
    parser.add_argument("--peak-prominence", type=float, default=0.01, metavar="FRACTION",
                        help="Minimum peak prominence, as a fraction of that spectrum's own max "
                             "intensity, for --experimental's peak-finding (default: 0.01 -- "
                             "filters out sub-1%% noise bumps in both spectra).")
    parser.add_argument("-v", "--version", action="version", version=f"stb-ramanAnalysis {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("RAMAN SPECTRUM WORKFLOW -- STAGE 3: ANALYSIS", 'bold'))
    print("-" * 60)

    stage2_report = os.path.join(args.directory, "raman_stage2.txt")
    if not os.path.isfile(stage2_report):
        print(color_text(f"[ERROR] '{stage2_report}' not found -- run stb-ramanModes (Stage 2) first.", 'red'))
        sys.exit(1)

    if args.experimental is not None and not os.path.isfile(args.experimental):
        print(color_text(f"[ERROR] --experimental file '{args.experimental}' not found.", 'red'))
        sys.exit(1)

    optical_root = os.path.join(args.directory, "optical_disp")
    label = args.label or detect_optical_label(optical_root)
    if not label:
        print(color_text(
            f"[ERROR] Could not detect SystemLabel from '{optical_root}' -- pass -l/--label.", 'red'))
        sys.exit(1)

    # The finite-difference delta used to build the +/-delta pairs -- Stage
    # 2's own MODE_TABLE doesn't persist it, but every folder's calc.fdf
    # doesn't either (the displacement is baked into structure.fdf, not
    # recorded as a number) -- reconstructed here isn't possible, so read
    # it back from Stage 2's own report text instead of re-deriving it.
    stage2_delta = None
    with open(stage2_report) as f:
        for line in f:
            if line.startswith("Displacement      :"):
                try:
                    stage2_delta = float(line.split(":")[1].split()[0])
                except (IndexError, ValueError):
                    pass
                break
    if stage2_delta is None:
        print(color_text(
            f"[ERROR] Could not recover the Stage 2 displacement value from '{stage2_report}'.", 'red'))
        sys.exit(1)

    rows = read_mode_table(stage2_report)
    if not rows:
        print(color_text(f"[ERROR] No # MODE_TABLE rows found in '{stage2_report}'.", 'red'))
        sys.exit(1)
    modes = group_by_mode(rows)

    report_path = os.path.join(args.directory, REPORT_FILE)
    with open(report_path, "w") as f_out:
        print_dual(f"{color_text('===== RAMAN STAGE 3 REPORT (ANALYSIS) =====', 'magenta')}", f_out)

        print_dual(f"\n{color_text('[0] RUN METADATA', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        print_dual(f"Date/time         : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", f_out)
        print_dual(f"Directory         : {args.directory}", f_out)
        print_dual(f"SystemLabel       : {label}", f_out)
        print_dual(f"Modes in table    : {len(modes)}", f_out)
        print_dual(f"Displacement      : {stage2_delta} Ang (from Stage 2)", f_out)

        print_dual(f"\n{color_text('[1] READING RUNS', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)

        mode_results = []  # (mode_index, freq_thz, r, activity, is_full, rho)
        computed = {}  # mode_index -> (r, activity, is_full, rho, scf_ok_all), for --skip-degenerate reuse
        for mode_index in sorted(modes):
            entry = modes[mode_index]
            freq_thz = entry["frequency_thz"]
            n_folders_total = len(entry["folders"])
            print_dual(f"  Mode {mode_index} ({freq_thz:.4f} THz / {freq_thz * THZ_TO_CM1:.2f} cm^-1):", f_out)

            derived_from = entry.get("derived_from")
            if derived_from is not None:
                rep = computed.get(derived_from)
                if rep is None:
                    print_dual(color_text(
                        f"    -> SKIP (representative mode {derived_from} has no usable Raman "
                        "tensor)", 'yellow'), f_out)
                    continue
                r, activity, is_full, rho, scf_ok_all = rep
                print_dual("    [Degenerate partner] Reusing mode "
                            f"{derived_from}'s Raman tensor by symmetry "
                            "(stb-ramanModes --skip-degenerate) -- activity and rho are "
                            "rotational invariants, identical for every partner in a degenerate "
                            "group.", f_out)
                poln_tag = ""
                if rho is not None:
                    poln_tag = "polarized" if rho < 0.75 - 1e-6 else "depolarized"
                offdiag_str = (f"  Rxy={r['xy']:.6f}  Rxz={r['xz']:.6f}  Ryz={r['yz']:.6f}"
                                if is_full else "")
                rho_str = f"  rho~{rho:.4f} ({poln_tag})" if rho is not None else ""
                print_dual(f"    -> Rxx={r['x']:.6f}  Ryy={r['y']:.6f}  Rzz={r['z']:.6f}{offdiag_str}  "
                            f"(1/Ang)  activity~{activity:.6f}{rho_str}  [derived]", f_out)
                mode_results.append((mode_index, freq_thz, r, activity, is_full, rho))
                computed[mode_index] = rep
                continue

            sidecar_path = os.path.join(optical_root, f"mode_{mode_index:02d}_tensor_form.json")
            n_missing = 0
            if os.path.isfile(sidecar_path):
                print_dual("    [Component reduction] Reconstructing the full tensor from a "
                            "single symmetry-fixed probe (stb-ramanModes --reduce-components).", f_out)
                r, scf_ok_all = reconstruct_reduced_tensor(
                    sidecar_path, entry["folders"], label, args.file, stage2_delta, f_out)
                has_offdiag = r is not None
                if r is None:
                    n_missing = n_folders_total
            else:
                static_eps, n_missing, scf_ok_all = compute_mode_tensor(
                    entry["folders"], label, args.file, f_out)
                r, has_offdiag = raman_tensor_full(static_eps, stage2_delta)
            activity, is_full, alpha_prime, gamma_prime2 = (
                (None, False, None, None) if r is None else raman_activity(r, has_offdiag))
            if activity is None:
                print_dual(color_text(
                    f"    -> SKIP (incomplete: {n_missing}/{n_folders_total} folder(s) "
                    "missing/unparseable)", 'yellow'), f_out)
            else:
                rho = depolarization_ratio(alpha_prime, gamma_prime2)
                poln_tag = ""
                if rho is not None:
                    poln_tag = "polarized" if rho < 0.75 - 1e-6 else "depolarized"
                offdiag_str = (f"  Rxy={r['xy']:.6f}  Rxz={r['xz']:.6f}  Ryz={r['yz']:.6f}"
                                if is_full else "")
                rho_str = f"  rho~{rho:.4f} ({poln_tag})" if rho is not None else ""
                print_dual(f"    -> Rxx={r['x']:.6f}  Ryy={r['y']:.6f}  Rzz={r['z']:.6f}{offdiag_str}  "
                            f"(1/Ang)  activity~{activity:.6f}{rho_str}"
                            + ("" if scf_ok_all else color_text("  [some folders unconverged]", 'yellow')), f_out)
                mode_results.append((mode_index, freq_thz, r, activity, is_full, rho))
                computed[mode_index] = (r, activity, is_full, rho, scf_ok_all)

        if not mode_results:
            print_dual(color_text(
                "\n[ERROR] No mode had a complete Raman tensor -- nothing to report/plot. "
                "Make sure SIESTA finished in every optical_disp/mode_*/ folder.", 'red'), f_out)
            sys.exit(1)

        all_full = all(is_full for *_, is_full, _ in mode_results)
        any_full = any(is_full for *_, is_full, _ in mode_results)

        print_dual(f"\n{color_text('[2] RAMAN-ACTIVE MODES SUMMARY', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)

        def _rho_cols(rho):
            if rho is None:
                return f"{'N/A':<10}{'N/A'}"
            return f"{rho:<10.4f}{'polarized' if rho < 0.75 - 1e-6 else 'depolarized'}"

        if any_full:
            print_dual(f"  {'Mode':<6}{'THz':<10}{'cm^-1':<10}{'Rxx':<12}{'Ryy':<12}{'Rzz':<12}"
                        f"{'Rxy':<12}{'Rxz':<12}{'Ryz':<12}{'Activity':<12}{'rho':<10}{'Polarization'}", f_out)
            for mode_index, freq_thz, r, activity, is_full, rho in mode_results:
                offdiag_cols = (f"{r['xy']:<12.6f}{r['xz']:<12.6f}{r['yz']:<12.6f}" if is_full
                                 else f"{'N/A':<12}{'N/A':<12}{'N/A':<12}")
                print_dual(f"  {mode_index:<6}{freq_thz:<10.4f}{freq_thz * THZ_TO_CM1:<10.2f}"
                            f"{r['x']:<12.6f}{r['y']:<12.6f}{r['z']:<12.6f}{offdiag_cols}"
                            f"{activity:<12.6f}{_rho_cols(rho)}", f_out)
        else:
            print_dual(f"  {'Mode':<6}{'THz':<12}{'cm^-1':<12}{'Rxx':<14}{'Ryy':<14}{'Rzz':<14}"
                        f"{'Activity (approx)':<20}{'rho':<10}{'Polarization'}", f_out)
            for mode_index, freq_thz, r, activity, is_full, rho in mode_results:
                print_dual(f"  {mode_index:<6}{freq_thz:<12.4f}{freq_thz * THZ_TO_CM1:<12.2f}"
                            f"{r['x']:<14.6f}{r['y']:<14.6f}{r['z']:<14.6f}"
                            f"{activity:<20.6f}{_rho_cols(rho)}", f_out)

        if all_full:
            print_dual(color_text(
                "\n[NOTE] Activity computed from the FULL Raman tensor (diagonal + off-diagonal "
                "Rxy/Rxz/Ryz), the exact literature invariant (45a^2+7gamma^2).", 'green'), f_out)
        elif any_full:
            print_dual(color_text(
                "\n[NOTE] Mixed scope: some modes have the full tensor (exact invariant), others "
                "only the diagonal one (approximate -- missing the off-diagonal Rxy/Rxz/Ryz "
                "anisotropy contribution). Re-run stb-ramanModes with --full-tensor for the "
                "modes still shown as diagonal-only.", 'yellow'), f_out)
        else:
            print_dual(color_text(
                "\n[NOTE] Activity is APPROXIMATE: computed from the diagonal Raman tensor "
                "(Rxx, Ryy, Rzz) only. The full literature invariant (45a^2+7gamma^2) also needs "
                "the off-diagonal tensor components (Rxy, Rxz, Ryz) -- re-run stb-ramanModes "
                "with --full-tensor to get those.", 'yellow'), f_out)
        print_dual(color_text(
            "[NOTE] rho (depolarization ratio) uses the same diagonal-only/full-tensor scope as "
            "activity above, per mode.", 'cyan'), f_out)

        print_dual(f"\n{color_text('[3] SPECTRUM', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        freqs_cm1 = [freq_thz * THZ_TO_CM1 for _, freq_thz, _, _, _, _ in mode_results]
        activities = [activity for _, _, _, activity, _, _ in mode_results]
        if args.temperature is not None:
            # Imaginary (negative-frequency) modes aren't real oscillators --
            # Bose-Einstein occupation is undefined for them (see
            # bose_einstein_weight's docstring) -- left unweighted (weight 1)
            # rather than raising, since Stage 2 already reports these as
            # [IMAGINARY] and the unweighted spectrum already includes them
            # as-is; thermal weighting just shouldn't touch them.
            n_imaginary_skipped = 0
            weights = []
            for _, freq_thz, _, _, _, _ in mode_results:
                if freq_thz <= 0:
                    weights.append(1.0)
                    n_imaginary_skipped += 1
                else:
                    weights.append(bose_einstein_weight(freq_thz, args.temperature))
            activities = [a * w for a, w in zip(activities, weights)]
            print_dual(f"Thermal weighting : Bose-Einstein Stokes factor (n+1) at "
                        f"{args.temperature:.1f} K applied to every mode's contribution.", f_out)
            if n_imaginary_skipped:
                print_dual(color_text(
                    f"[WARNING] {n_imaginary_skipped} imaginary-frequency mode(s) left unweighted "
                    "(Bose-Einstein occupation is undefined for a non-real vibration).", 'yellow'), f_out)
        grid, intensity = build_lorentzian_spectrum(freqs_cm1, activities, args.linewidth)
        dat_path = os.path.join(args.directory, f"{args.output}.dat")
        gplot_path = os.path.join(args.directory, f"{args.output}.gplot")

        experimental, experimental_dat_path = None, None
        if args.experimental is not None:
            exp_freq, exp_intensity = read_experimental_spectrum(args.experimental)
            experimental = (exp_freq, exp_intensity)
            experimental_dat_path = os.path.join(args.directory, f"{args.output}_experimental.dat")

        write_spectrum_plot(dat_path, gplot_path, grid, intensity, all_full, args.temperature,
                             experimental, experimental_dat_path)
        print_dual(f"{color_text('[Saved]', 'cyan')} {dat_path}, {gplot_path} "
                    f"(cd {args.directory} && gnuplot {os.path.basename(gplot_path)})", f_out)

        extra_files = ""
        if experimental is not None:
            print_dual(f"\n{color_text('[3b] EXPERIMENTAL COMPARISON', 'magenta')}", f_out)
            print_dual("-" * 60, f_out)
            print_dual(f"Experimental file : {args.experimental} ({len(exp_freq)} point(s))", f_out)
            sim_peaks = find_spectrum_peaks(grid, intensity, args.peak_prominence)
            exp_peaks = find_spectrum_peaks(exp_freq, exp_intensity, args.peak_prominence)
            print_dual(f"Peaks found       : {len(sim_peaks)} simulated, {len(exp_peaks)} experimental", f_out)
            if len(exp_peaks) == 0:
                print_dual(color_text(
                    "[WARNING] No experimental peaks found above the prominence threshold -- "
                    "try lowering --peak-prominence.", 'yellow'), f_out)
            elif len(sim_peaks) == 0:
                print_dual(color_text(
                    "[WARNING] No simulated peaks found -- nothing to match against.", 'yellow'), f_out)
            else:
                matches = match_peaks(sim_peaks, exp_peaks)
                print_dual(f"  {'Exp. peak (cm^-1)':<20}{'Sim. peak (cm^-1)':<20}{'Delta (cm^-1)'}", f_out)
                for exp_f, sim_f, delta in matches:
                    print_dual(f"  {exp_f:<20.2f}{sim_f:<20.2f}{delta:+.2f}", f_out)
                mean_abs_delta = float(np.mean([abs(d) for _, _, d in matches]))
                print_dual(f"Mean |delta|      : {mean_abs_delta:.2f} cm^-1 (average across "
                            f"{len(matches)} matched pair(s))", f_out)
            print_dual(f"{color_text('[Saved]', 'cyan')} {experimental_dat_path}", f_out)
            extra_files = f", {experimental_dat_path}"

        print_dual(f"\n{color_text('[4] SUMMARY & FILES', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        print_dual(f"Modes analyzed      : {len(mode_results)}/{len(modes)}", f_out)
        print_dual(f"Report              : {report_path}", f_out)
        print_dual(f"Files               : {dat_path}, {gplot_path}{extra_files}", f_out)

    print("\n[INFO] Complete job!")
    print("\n" + "-" * 60)
    print(color_text("Raman analysis complete.\n", 'bold'))


if __name__ == "__main__":
    main()
