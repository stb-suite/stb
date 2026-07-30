#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
#     Mechanical Properties Analyzer            #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""Stage 2 of 2: reads back the finished SIESTA runs from stb-strain's own
'strain_*' folders (flat or nested -- see find_strain_folders) and fits a
stress-strain curve per direction: initial modulus, peak (UTS), critical
strain, toughness, and optionally the 0.2% offset yield. Multiple
directions found under --dir are compared automatically (no flag needed --
find_strain_folders/collect_strain_groups already report however many are
present). Dimensionality (3D/2D/1D, controlling report units GPa/N/m/nN) is
auto-detected from a real .fdf structure file if one is found alongside the
SIESTA output (same vacuum-axis convention as stb-strain/stb-mlelastic --
see find_structure_fdf), falling back to 3D when none is found (the common,
expected case for hand-built/legacy folders holding only a bare SIESTA
.out) -- pass --dimensionality explicitly to override either way.

Report/plot style matches the rest of the suite's "stb-standard report"
convention (aimd_analysis.py/mlelastic.py/etc.): a numbered [0]...[5]
report (always printed), --save-report (persists it to a file),
--save-gnuplot (writes the curve data + a companion .gplot script --
already-existing write_strain_gplot/write_compare_gplot, just gated behind
this flag instead of unconditional), --view (an on-demand interactive
matplotlib preview, never saved to disk -- gnuplot is the only
persisted-plot mechanism here, matching the rest of this convention), and
-o/--output-dir.
"""

VERSION = "2.0.0"

import os
import sys
import glob
import argparse
import numpy as np
from scipy.stats import linregress
# Use trapezoid for compatibility with new Scipy
from scipy.integrate import trapezoid
from stb.core import siesta_log, structure_io, kspace
from stb.core.cli import COLORS, color_text, show_intro, print_dual, print_section, print_table

REPORT_FILE = "stb_strainAnalysis_report.txt"

# ==========================================
#           CORE LOGIC
# ==========================================

def calculate_yield_stress(strain, stress, modulus, offset=0.002):
    """
    Calculates 0.2% Offset Yield Strength.
    Returns (Yield Stress, Yield Strain).

    Borrowed from macroscopic metallurgy (empirical onset of dislocation-
    glide plasticity) -- only reported when --yield is passed, since a
    perfect, defect-free periodic crystal under affine (Cauchy-Born) strain
    has no dislocation-nucleation mechanism to actually yield through.
    """
    # Offset line: Stress = E * (Strain - 0.002)
    offset_line = modulus * (strain - offset)

    idx = np.where(strain > 0)[0]
    if len(idx) < 2: return 0.0, 0.0

    diff = stress[idx] - offset_line[idx]
    for i in range(len(diff)-1):
        if diff[i] * diff[i+1] < 0:
            return stress[idx][i], strain[idx][i]

    return 0.0, 0.0


def compute_cross_section_area(cell, strained_axis_index):
    """Area (Ang^2) of the cell cross-section perpendicular to the strained
    axis -- the parallelogram formed by the other 2 lattice vectors. Used to
    convert kBar -> axial Force (nN) for 1D mode: SIESTA's stress is already
    normalized by this same area (vacuum included), so stress*area recovers
    the real physical force independent of how much vacuum padding was used
    to build the wire/tube's simulation cell.
    """
    other = [i for i in range(3) if i != strained_axis_index]
    return float(np.linalg.norm(np.cross(cell[other[0]], cell[other[1]])))


def analyze_mechanics(data, direction, is_2d=False, z_height=20.0,
                       is_1d=False, cell_area=None, cross_section=None):
    """
    Calculates the initial slope, peak (UTS), Critical Strain, Toughness.
    data: [Strain%, StrainFrac, Sxx, Syy, Szz, Syz, Sxz, Sxy] (kBar)

    Unit/quantity depends on dimensionality:
    - 3D: Stress (GPa), kBar -> GPa direct conversion.
    - 2D: Stress (N/m), corrected by the vacuum (Z) height.
    - 1D: Force (nN) -- NOT stress. kBar is already normalized by the
      full cell's cross-section (vacuum included); multiplying by that same
      cross-section area recovers the real axial force, independent of the
      arbitrary vacuum padding chosen when the wire/tube cell was built.
      Toughness in this mode is Work/length (nN, dimensionally = J/m).
    """
    quantity = "Force" if is_1d else "Stress"

    if is_1d:
        factor = cell_area * 1e-3  # kBar * Ang^2 -> nN
        unit = "nN"
    elif is_2d:
        # Convert kBar to N/m using Z-height
        # 1 kBar = 1e8 Pa (N/m^2); Z_height in Angstrom (1e-10 m)
        # Stress_2D = Stress_kBar * 1e8 * Z * 1e-10 = Stress_kBar * Z * 0.01
        factor = 0.01 * z_height
        unit = "N/m"
    else:
        # Standard 3D: kBar -> GPa (1 kBar = 0.1 GPa)
        factor = 0.1
        unit = "GPa"

    strain_pct = data[:, 0]
    strain_frac = data[:, 1]
    stress_raw = data[:, 2:8] # xx, yy, zz, yz, xz, xy

    if len(strain_pct) < 2:
        raise ValueError(
            f"Only {len(strain_pct)} valid strain step(s) found for direction '{direction}' -- "
            "at least 2 are needed to fit a slope (check for failed/unconverged SCF runs).")

    # Identify Axial Component
    # xy/xz/yz (biaxial) map to their first axis -- see the [WARNING] printed
    # in main() when this happens; there is no single "biaxial modulus" here,
    # just the stress along the first strained axis.
    d_map = {'xx': 0, 'yy': 1, 'zz': 2, 'x': 0, 'y': 1, 'z': 2,
             'xy': 0, 'xz': 0, 'yz': 1}
    if direction.lower() not in d_map:
        print(f"{color_text('[WARNING]', 'yellow')} Unrecognized direction '{direction}' -- "
              "defaulting to the xx stress component (not a validated mapping).")
    ax_idx = d_map.get(direction.lower(), 0) # Default to xx if unknown

    sigma_axial = stress_raw[:, ax_idx] * factor

    # --- Initial Slope (dQuantity/dStrain) ---
    # Linear fit over a genuinely small-strain window: prefer every point
    # within +-2% (the classic small-strain-elasticity range); but if the
    # sweep's own step size is coarser than that (very common for a large
    # -deformation sweep, e.g. 2%/step up to 40%), fewer than
    # MIN_FIT_POINTS points fall in that window. The OLD behavior fell back
    # to fitting the ENTIRE sweep (including large, manifestly nonlinear/
    # plastic strain) in that case -- not just imprecise but physically
    # wrong: verified live on a real SIESTA 0-40% sweep (2% steps, so only
    # the 0%/2% points fall in +-2%), that fallback gave R^2 as low as
    # ~0.01 (nearly uncorrelated) instead of a real elastic-regime slope.
    # Falling back to the MIN_FIT_POINTS points nearest to zero strain
    # instead (still a genuine small-strain fit, just wider than +-2% when
    # the sampling is coarse) recovers R^2 > 0.99 on the same data.
    MIN_FIT_POINTS = 3
    limit_idx = np.where(np.abs(strain_pct) <= 2.0)[0]
    fit_window_expanded = len(limit_idx) < MIN_FIT_POINTS
    if fit_window_expanded:
        order = np.argsort(np.abs(strain_pct))
        limit_idx = order[:min(MIN_FIT_POINTS, len(strain_pct))]

    slope, intercept, r_val, _, fit_std_err = linregress(
        strain_frac[limit_idx], sigma_axial[limit_idx])
    fit_n_points = len(limit_idx)
    fit_strain_lo = float(strain_pct[limit_idx].min())
    fit_strain_hi = float(strain_pct[limit_idx].max())

    # --- Peak (UTS / max Force) + Critical Strain ---
    uts_idx = int(np.argmax(sigma_axial))
    uts = sigma_axial[uts_idx]
    uts_strain = strain_pct[uts_idx]
    # True peak may lie outside the tested range if it falls on either edge
    # of the (strain-sorted) sweep.
    uts_at_boundary = uts_idx in (0, len(sigma_axial) - 1)

    # --- Toughness (Energy Density for 3D/2D; Work/length for 1D) ---
    toughness = trapezoid(sigma_axial, strain_frac)

    # --- Yield Strength (0.2% Offset) -- only used if --yield is passed ---
    yield_stress, yield_strain = calculate_yield_stress(strain_frac, sigma_axial, slope)

    result = {
        'modulus': slope,
        'quantity': quantity,
        'uts': uts,
        'uts_idx': uts_idx,
        'uts_strain': uts_strain,
        'uts_at_boundary': uts_at_boundary,
        'toughness': toughness,
        'yield': yield_stress,
        'r_squared': r_val**2,
        'unit': unit,
        'slope': slope,
        'intercept': intercept,
        'fit_std_err': fit_std_err,
        'fit_n_points': fit_n_points,
        'fit_strain_lo': fit_strain_lo,
        'fit_strain_hi': fit_strain_hi,
        'fit_window_expanded': fit_window_expanded,
        'strain_pct': strain_pct,
        'strain_frac': strain_frac,
        'axial_stress': sigma_axial,
        'ax_idx': ax_idx,
        'stress_raw': stress_raw,
        'factor': factor,
    }

    if is_1d and cross_section:
        # Conventional Stress (GPa) using a user-supplied physical
        # cross-section instead of the (vacuum-inflated) cell area --
        # Stress_GPa = Force_nN / cross_section_Ang^2 * 100 (derived and
        # cross-checked: reduces to the plain 3D kBar->GPa factor when
        # cross_section == cell_area).
        result['conventional_stress'] = sigma_axial / cross_section * 100.0
        result['cross_section'] = cross_section

    return result


def find_strain_folders(base_dir):
    """Finds every 'strain_*' run folder reachable from base_dir, supporting
    2 layouts: flat (strain_<dir>_<pct> directly under base_dir -- what you
    get by pointing base_dir at a single direction's own subfolder, e.g.
    'strain_runs/x') and nested (stb-strain's own default output layout: one
    <direction>/ subfolder per strain direction under base_dir, each holding
    that direction's strain_<direction>_<pct> folders). This lets --dir point
    at either a single direction's subfolder (this direction only) or the
    top-level output directory itself (every direction found under it,
    compared automatically) without the caller needing to know which layout
    is present.
    """
    direct = [os.path.join(base_dir, d) for d in sorted(os.listdir(base_dir))
              if os.path.isdir(os.path.join(base_dir, d)) and d.startswith('strain_')]
    if direct:
        return direct
    nested = []
    for d in sorted(os.listdir(base_dir)):
        sub = os.path.join(base_dir, d)
        if not os.path.isdir(sub):
            continue
        nested.extend(
            os.path.join(sub, e) for e in sorted(os.listdir(sub))
            if os.path.isdir(os.path.join(sub, e)) and e.startswith('strain_'))
    return nested


def _canonical_direction(direction):
    """Folds the doubled-letter uniaxial convention ('xx'/'yy'/'zz', used by
    some older fixtures) into the canonical single-letter form ('x'/'y'/'z'),
    so both spellings of the same physical direction group together instead
    of being treated as 2 different directions (which could otherwise
    spuriously trigger the "multiple directions found" comparison for what
    is really 1 direction). Genuine biaxial pairs ('xy'/'xz'/'yz') are
    returned unchanged.
    """
    d = direction.lower()
    if len(d) == 2 and d[0] == d[1]:
        return d[0]
    return d


def _is_uniaxial(direction):
    """True for 'x'/'y'/'z' and also the doubled-letter convention 'xx'/'yy'/
    'zz' (older fixtures use this Voigt-style naming for the same uniaxial
    case) -- False only for a genuine biaxial pair like 'xy'/'xz'/'yz'.
    """
    return len(direction) == 1 or direction[0] == direction[1]


def _first_folder_for_direction(folders, direction):
    """First folder in `folders` whose parsed, canonicalized direction
    matches `direction` -- used to read the reference cell for 1D mode's
    cross-section from a folder that actually belongs to that direction,
    not an arbitrary one that may have been strained along a different axis.
    """
    for f in folders:
        d, _ = siesta_log.parse_strain_folder_name(f)
        if d is not None and _canonical_direction(d) == direction:
            return f
    return None


def collect_strain_groups(siesta_out, folders):
    """Scans the given 'strain_*' folders, groups them by detected direction
    (instead of assuming a single direction is present), and reads the Voigt
    stress from each. Returns {direction: sorted_data_array}, data_array
    columns [Strain%, StrainFrac, Sxx, Syy, Szz, Syz, Sxz, Sxy] (kBar).
    `folders` entries may be full/relative paths (e.g. 'strain_runs/strain_x_1.00').
    Direction keys are canonicalized (see _canonical_direction) so 'xx' and
    'x' fixtures group together.
    """
    groups = {}
    for f in folders:
        d, val = siesta_log.parse_strain_folder_name(f)
        if d is None:
            continue
        d = _canonical_direction(d)
        fpath = os.path.join(f, siesta_out)
        if not os.path.exists(fpath):
            continue
        stress = siesta_log.get_stress_voigt_kbar(fpath)
        if not stress:
            continue
        groups.setdefault(d, []).append([val, val / 100.0] + stress)
        print(f"   Reading {f}: Strain {val:>5.2f}% OK")

    return {d: np.array(sorted(rows, key=lambda x: x[0])) for d, rows in groups.items()}


def find_structure_fdf(folder):
    """First *.fdf in `folder` (sorted, excluding "config_extra.fdf" -- the
    fixed name stb-strain writes for its own %block Geometry.Constraints,
    which never has a LatticeVectors/AtomicCoordinates block of its own)
    that structure_io.read_fdf can parse. Returns the parsed FdfStructure,
    or None if no .fdf exists in the folder or none of them parse -- the
    common, expected case for hand-built/legacy strain_* folders that only
    ever hold a bare SIESTA .out. Deliberately does NOT import strain.py's
    own EXTRA_FDF_FILE constant (hardcoded literal instead) -- this tool
    must stay usable on folders stb-strain itself never touched.
    """
    for path in sorted(glob.glob(os.path.join(folder, "*.fdf"))):
        if os.path.basename(path) == "config_extra.fdf":
            continue
        try:
            return structure_io.read_fdf(path)
        except (ValueError, FileNotFoundError):
            continue
    return None


def _dimensionality_label(is_2d, is_1d):
    if is_1d:
        return "1D (Wire/Tube)"
    if is_2d:
        return "2D (Sheet)"
    return "3D (Bulk)"


def print_mechanics_block(direction, results, is_2d, is_1d=False, show_yield=False, f_out=None):
    """Prints one direction's mechanical properties via print_dual (console,
    and also the report file when f_out is given) -- shared by the single
    -direction and multi-direction (auto-compared) report paths so the two
    stay in sync.
    """
    u = results['unit']
    qty = results['quantity']
    modulus_label = f"Initial Slope (d{'F' if is_1d else chr(963)}/d{chr(949)})"
    modulus_str = f"{results['modulus']:.4f} {u}"
    peak_str = f"{results['uts']:.4f} {u}"
    peak_label = f"Peak {qty}"

    print_dual(f"Direction       : {direction.upper()}", f_out)
    print_dual(f"Dimensionality  : {_dimensionality_label(is_2d, is_1d)}", f_out)
    print_dual("-" * 50, f_out)
    print_dual(f"{modulus_label:<16}: {color_text(modulus_str, 'green')} (R²={results['r_squared']:.4f}, "
               f"std.err.={results['fit_std_err']:.4f} {u})", f_out)
    print_dual(f"  Fit window    : {results['fit_n_points']} point(s), strain "
               f"{results['fit_strain_lo']:.4f}% to {results['fit_strain_hi']:.4f}%", f_out)
    if results['fit_window_expanded']:
        print_dual(color_text(
            "  [NOTE] Fewer than 3 strain steps fall within +-2% -- the fit window was "
            "widened to the nearest few points instead of silently fitting the entire "
            "(likely nonlinear/plastic at large strain) sweep.", 'cyan'), f_out)
    print_dual(f"{peak_label:<16}: {color_text(peak_str, 'red')}", f_out)
    print_dual(f"Critical Strain : {results['uts_strain']:.4f} %", f_out)
    if not is_1d:
        print_dual(color_text(
            "[INFO] Single-direction slope under a clamped transverse cell -- for the "
            "rigorous small-strain elastic-tensor Young's Modulus, use stb-elasticAnalysis.",
            'cyan'), f_out)
    if results['uts_at_boundary']:
        print_dual(color_text(
            "[WARNING] Peak occurred at the edge of the tested strain range -- "
            "the true peak may lie beyond --stmin/--stmax; treat it as a lower bound.",
            'yellow'), f_out)
    if is_1d and 'conventional_stress' in results:
        conv_peak = results['conventional_stress'][results['uts_idx']]
        print_dual(f"Peak Stress (conventional, cross-section={results['cross_section']:.4f} "
                   f"Ang^2): {conv_peak:.4f} GPa", f_out)

    toughness_unit = "nN (work/length, ≡ J/m)" if is_1d else ('J/m^2' if is_2d else 'GJ/m^3')
    print_dual(f"Toughness       : {results['toughness']:.4f} {toughness_unit}", f_out)

    if show_yield:
        if results['yield'] > 0:
            print_dual(f"Yield {qty} (0.2%): {results['yield']:.4f} {u}", f_out)
        else:
            print_dual(f"Yield {qty} (0.2%): {color_text('Not detected (Linear)', 'yellow')}", f_out)
        print_dual(color_text(
            "[INFO] '0.2% offset yield' is a macroscopic-plasticity concept borrowed from "
            "metallurgy (dislocation-glide onset) -- a defect-free periodic crystal under "
            "affine strain has no such mechanism; interpret with care.", 'cyan'), f_out)

    if not is_1d and _is_uniaxial(direction):
        transverse = _transverse_stress_lines(results)
        if transverse:
            print_dual("Transverse " + qty + " @ peak (diagnostic, clamped cell): " +
                       ", ".join(transverse), f_out)
            print_dual(color_text(
                "[INFO] Not a Poisson's ratio (that needs transverse STRAIN, not stress) -- "
                "just how much elastic response this direction's clamped transverse cell is "
                "suppressing.", 'cyan'), f_out)


def _transverse_stress_lines(results):
    """The other on-diagonal Voigt normal components (not the analyzed
    axial one), evaluated at the peak (UTS) strain step -- a diagnostic for
    how much the clamped-transverse-cell assumption suppresses real elastic
    response, not a derived Poisson's ratio.
    """
    axis_names = ['xx', 'yy', 'zz']
    ax_idx = results['ax_idx']
    uts_idx = results['uts_idx']
    factor = results['factor']
    stress_raw = results['stress_raw']
    lines = []
    for idx in range(3):
        if idx == ax_idx:
            continue
        val = stress_raw[uts_idx, idx] * factor
        lines.append(f"{axis_names[idx]}={val:.4f} {results['unit']}")
    return lines


def _write_direction_block(f, results, label):
    """Appends one direction's 2 gnuplot 'index' blocks (curve, peak point)
    to an already-open data file handle. Columns match np.savetxt's curve
    file exactly (Strain%, StrainFrac, Axial_Quantity), so both the
    standalone curve file and this helper stay reusable with 'using 2:3'.
    """
    strain_pct = results['strain_pct']
    strain_frac = results['strain_frac']
    sigma_axial = results['axial_stress']
    uts_idx = results['uts_idx']

    f.write(f"# {label}: curve\n")
    f.write(f"# Strain(%) Strain(Frac) {results['quantity']}\n")
    for p, e, s in zip(strain_pct, strain_frac, sigma_axial):
        f.write(f"{p:12.6f}  {e:12.6f}  {s:12.6f}\n")
    f.write("\n\n")
    f.write(f"# {label}: peak point\n")
    f.write(f"# Strain(%) Strain(Frac) {results['quantity']}\n")
    f.write(f"{strain_pct[uts_idx]:12.6f}  {strain_frac[uts_idx]:12.6f}  {sigma_axial[uts_idx]:12.6f}\n")


def write_strain_gplot(curve_filename, gplot_filename, direction, results, unit):
    """Writes the companion .gplot for a single direction's curve, the
    initial-slope fit line, and the peak point -- same data-file + .gplot
    convention used throughout the suite (stb-hubbarduAnalysis, stb-kgrid,
    stb-density, stb-workfunction, stb-xrd). Appends the peak point as a
    second gnuplot 'index' block to the curve file already written by
    np.savetxt (index 0), instead of duplicating the data in a second file.
    """
    qty = results['quantity']
    with open(curve_filename, 'a') as f:
        f.write("\n\n")
        f.write("# index 1: peak point\n")
        f.write(f"# Strain(%) Strain(Frac) {qty}\n")
        uts_idx = results['uts_idx']
        f.write(f"{results['strain_pct'][uts_idx]:12.6f}  "
                f"{results['strain_frac'][uts_idx]:12.6f}  "
                f"{results['axial_stress'][uts_idx]:12.6f}\n")

    # The fit line is only drawn across the strain window it was actually
    # fit over (see analyze_mechanics) -- NOT the full plotted strain range
    # -- so a small-strain elastic fit never visually implies the material
    # follows that straight line all the way out to a large-deformation
    # sweep's peak/plastic regime. Gnuplot has no native per-curve x-domain
    # clip, so this uses the standard "undefined outside the domain"
    # ternary idiom (1/0 -> NaN, silently skipped when plotting).
    fit_lo = results['fit_strain_lo'] / 100.0
    fit_hi = results['fit_strain_hi'] / 100.0
    lines = [
        '# --- STB Plot Configuration ---\n',
        '# Generated by stb-strainAnalysis\n',
        'set terminal pdfcairo enhanced color font "Arial,14" size 7,5\n',
        f'set output "{gplot_filename.rsplit(".", 1)[0]}.pdf"\n\n',
        f'set title "{qty}-Strain Curve ({direction.upper()})"\n',
        'set xlabel "Strain (fraction)"\n',
        f'set ylabel "{qty} ({unit})"\n',
        'set grid\n',
        'set key top left\n',
        f"f(x) = {results['slope']:.8f}*x + {results['intercept']:.8f}\n",
        (f'plot "{curve_filename}" index 0 using 2:3 with linespoints pt 7 ps 1.2 lc rgb "#2255cc" '
         f'title "{qty}-strain data", \\\n'
         f'     (x>={fit_lo:.8f} && x<={fit_hi:.8f}) ? f(x) : 1/0 with lines lc rgb "#2255cc" '
         'dt 2 title "Initial slope fit", \\\n'
         f'     "{curve_filename}" index 1 using 2:3 with points pt 9 ps 2 lc rgb "#cc2222" '
         'title "Peak"\n')
    ]
    with open(gplot_filename, 'w') as f:
        f.writelines(lines)


def write_compare_gplot(data_filename, gplot_filename, per_direction, unit):
    """Same convention as write_strain_gplot, generalized to N directions:
    one combined data file (2 'index' blocks per direction: curve + peak
    point) and one .gplot overlaying every direction's curve + fit, one
    color per direction. Generalizes the fixed-2-branch multi-index pattern
    already used by hubbardu_analysis.py::write_response_plot.
    """
    colors = ['#2255cc', '#cc5522', '#22aa55', '#aa22aa', '#aaaa22', '#22aaaa']
    directions = list(per_direction.keys())
    qty = next(iter(per_direction.values()))['quantity']

    with open(data_filename, 'w') as f:
        for i, direction in enumerate(directions):
            if i > 0:
                f.write("\n\n")
            _write_direction_block(f, per_direction[direction], f"index {2*i}/{2*i+1} ({direction})")

    lines = [
        '# --- STB Plot Configuration ---\n',
        '# Generated by stb-strainAnalysis (auto-compare)\n',
        'set terminal pdfcairo enhanced color font "Arial,14" size 7,5\n',
        f'set output "{gplot_filename.rsplit(".", 1)[0]}.pdf"\n\n',
        f'set title "{qty}-Strain Comparison"\n',
        'set xlabel "Strain (fraction)"\n',
        f'set ylabel "{qty} ({unit})"\n',
        'set grid\n',
        'set key top left\n',
    ]
    plot_parts = []
    for i, direction in enumerate(directions):
        results = per_direction[direction]
        color = colors[i % len(colors)]
        idx0, idx1 = 2 * i, 2 * i + 1
        # Same fit-window clipping as write_strain_gplot (see its own
        # comment) -- each direction's fit line only spans the strain
        # range it was actually fit over, not the full comparison plot.
        fit_lo = results['fit_strain_lo'] / 100.0
        fit_hi = results['fit_strain_hi'] / 100.0
        lines.append(f"f_{i}(x) = {results['slope']:.8f}*x + {results['intercept']:.8f}\n")
        plot_parts.append(
            f'"{data_filename}" index {idx0} using 2:3 with linespoints pt 7 ps 1.2 lc rgb "{color}" '
            f'title "{direction.upper()}"')
        plot_parts.append(
            f'(x>={fit_lo:.8f} && x<={fit_hi:.8f}) ? f_{i}(x) : 1/0 with lines lc rgb "{color}" '
            'dt 2 notitle')
        plot_parts.append(
            f'"{data_filename}" index {idx1} using 2:3 with points pt 9 ps 2 lc rgb "{color}" notitle')
    lines.append('plot ' + ', \\\n     '.join(plot_parts) + '\n')

    with open(gplot_filename, 'w') as f:
        f.writelines(lines)

# ==========================================
#           MAIN
# ==========================================

def main():
    parser = argparse.ArgumentParser(
        description="Mechanical Properties from Stress-Strain Data (Stage 2 of stb-strain).",
        epilog="Example usage:\n"
               "  stb-strainAnalysis --file calc.out\n"
               "  stb-strainAnalysis --file calc.out --dir strain_runs --save-report --save-gnuplot --view\n"
               "  stb-strainAnalysis --file calc.out --dimensionality 1d --cross-section 45.2\n",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("-f", "--file", required=True, help="Siesta output file (e.g., calc.out).")
    parser.add_argument("--dir", default="strain_runs",
                        help="Directory to scan for 'strain_*' run folders (default: strain_runs). "
                             "Accepts either stb-strain's own output layout -- point this at a single "
                             "direction's own subfolder (e.g. strain_runs/x) for that direction alone, "
                             "or at the top-level output directory itself to automatically compare "
                             "every direction found under it -- or a flat directory of 'strain_*' "
                             "folders.")
    parser.add_argument("--dimensionality", choices=["auto", "3d", "2d", "1d"], default="auto",
                        help="Physical dimensionality, controlling report units (GPa/N/m/nN) and how "
                             "the cross-section for 1D Force mode is computed. 'auto' (default) tries "
                             "to read a real, parseable .fdf structure file from the first strain "
                             "folder found and classify it via vacuum-axis detection (same convention "
                             "as stb-strain/stb-mlelastic); if none is found (the common, expected "
                             "case for hand-built/legacy folders holding only a bare SIESTA .out), "
                             "falls back to 3d.")
    parser.add_argument("--vacuum-gap", type=float, default=10.0,
                        help="Vacuum-gap threshold in Ang, used only when --dimensionality auto finds "
                             "a real .fdf structure file to classify (default: 10.0, matching "
                             "stb-strain/stb-mlelastic).")
    parser.add_argument("--thickness", type=float, default=20.0,
                        help="Vacuum height (Z) for 2D conversion (Angstrom).")
    parser.add_argument("--cross-section", type=float, default=None,
                        help="Physical cross-section area (Ang^2) of the wire/tube wall, only used "
                             "when the run is 1D (--dimensionality 1d, or auto resolving to 1d), to "
                             "additionally report a conventional Stress (GPa) -- there is no way to "
                             "auto-derive this (unlike the 2D vacuum height), since a wire/tube wall "
                             "thickness is a physical convention, not a cell property.")
    parser.add_argument("--yield", dest="show_yield", action="store_true",
                        help="Also report the 0.2%% offset Yield -- a macroscopic-plasticity concept "
                             "borrowed from metallurgy; off by default since a defect-free periodic "
                             "crystal under affine DFT strain has no dislocation-nucleation mechanism "
                             "to actually yield through.")
    parser.add_argument("-o", "--output-dir", default=".",
                        help="Directory to write <direction>_curve.dat/.gplot (single direction) or "
                             "comparison_curve.dat/.gplot (multiple directions) into with "
                             f"--save-gnuplot, and {REPORT_FILE} into with --save-report (default: "
                             "current directory). Created if it doesn't exist.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--save-gnuplot", action="store_true",
                        help="Also write the curve data + a companion .gplot script (single "
                             "direction: <direction>_curve.dat/.gplot; multiple directions: "
                             "comparison_curve.dat/.gplot). Off by default.")
    parser.add_argument("--view", action="store_true",
                        help="Show an interactive matplotlib preview (the stress/force-strain curve, "
                             "fit line, and peak) before finishing. Off by default.")
    parser.add_argument("-v", "--version", action="version",
                        version=f"stb-strainAnalysis {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.cross_section is not None and args.dimensionality in ("3d", "2d"):
        parser.error("--cross-section only applies together with --dimensionality 1d (or auto, "
                     "if it resolves to 1d).")

    if args.intro:
        show_intro([
            f"Siesta Tool Box - Mechanical Analysis v{VERSION}",
            "Developed by Dr. Carlos M. O. Bastos",
            "bastoscmo.github.io"
        ])

    os.makedirs(args.output_dir, exist_ok=True)
    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    def fail(message):
        print_dual(color_text(f"[FAIL] {message}", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_dual(color_text(
        "===== STB-STRAINANALYSIS STAGE 2 REPORT (STRESS-STRAIN ANALYSIS) =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Directory         : {args.dir}", f_out)
    print_dual(f"Output filename   : {args.file}", f_out)
    print_dual(f"Dimensionality    : {args.dimensionality}", f_out)
    print_dual(f"Output dir        : {args.output_dir}", f_out)
    print_dual(f"Save report       : {'yes' if args.save_report else 'no'}", f_out)
    print_dual(f"Save gnuplot      : {'yes' if args.save_gnuplot else 'no'}", f_out)
    print_dual(f"View (matplotlib) : {'yes' if args.view else 'no'}", f_out)
    if report_path:
        print_dual(f"Report file       : {report_path}", f_out)

    print_section("[1] INPUT DATA", f_out)
    if not os.path.isdir(args.dir):
        fail(f"Directory '{args.dir}' not found.")
    folders = find_strain_folders(args.dir)
    if not folders:
        fail(f"No 'strain_*' folders found in '{args.dir}'.")
    print_dual(f"Found {len(folders)} strain step folder(s) under '{args.dir}'.", f_out)
    groups = collect_strain_groups(args.file, folders)
    if not groups:
        fail("No stress data found.")
    print_table(["Direction", "Steps", "Strain range (%)"], [
        ([d.upper(), str(len(data)), f"{data[:, 0].min():.2f} to {data[:, 0].max():.2f}"], None)
        for d, data in sorted(groups.items())
    ], f_out)
    for direction in sorted(groups):
        if not _is_uniaxial(direction):
            print_dual(color_text(
                f"[WARNING] Biaxial direction '{direction}' detected -- reporting only the "
                f"stress component along its first axis ('{direction[0]}{direction[0]}'), not "
                "a true biaxial modulus.", 'yellow'), f_out)

    print_section("[2] DIMENSIONALITY DETECTION", f_out)
    if args.dimensionality != "auto":
        dimensionality = args.dimensionality
        print_dual(f"Dimensionality    : {dimensionality.upper()} (manual override)", f_out)
    else:
        structure = find_structure_fdf(folders[0])
        if structure is None:
            print_dual(
                f"[NOTE] No parseable .fdf structure file found in '{folders[0]}' -- assuming "
                "3D bulk. This is the expected case for a hand-built/legacy strain_* folder "
                "holding only a bare SIESTA output; pass --dimensionality explicitly to "
                "override.", f_out)
            dimensionality = "3d"
        else:
            positions = np.array([pos for _, pos in structure.atoms])
            is_cartesian = structure.coord_format == 'cartesian'
            frac_coords = kspace.to_fractional(positions, structure.lattice, is_cartesian)
            vacuum_axes = kspace.detect_vacuum_axes(frac_coords, structure.lattice, args.vacuum_gap)
            print_dual(f"[INFO] Auto-detected from a real structure file found in "
                       f"'{folders[0]}' -- dimensionality: {kspace.dimensionality_label(vacuum_axes)}",
                       f_out)
            if vacuum_axes == [False, False, True]:
                dimensionality = "2d"
            elif vacuum_axes == [True, True, False]:
                dimensionality = "1d"
            elif not any(vacuum_axes):
                dimensionality = "3d"
            else:
                print_dual(color_text(
                    "[WARNING] Vacuum padding detected on an axis this tool doesn't support "
                    "(expects in-plane xy for 2D, or wire-along-c for 1D) -- falling back to "
                    "3D (GPa). Pass --dimensionality explicitly to override.", 'yellow'), f_out)
                dimensionality = "3d"
        print_dual(f"Dimensionality    : {dimensionality.upper()} (auto-detected)", f_out)

    is_2d = dimensionality == "2d"
    is_1d = dimensionality == "1d"

    if args.cross_section is not None and not is_1d:
        fail("--cross-section only applies when the run is 1D -- resolved dimensionality here "
             f"is {dimensionality.upper()}. Pass --dimensionality 1d explicitly if you believe "
             "this is wrong.")

    # Detect Z-height automatically if 2D
    if is_2d and os.path.exists(os.path.join(folders[0], args.file)):
        z_auto = siesta_log.get_cell_height(os.path.join(folders[0], args.file))
        if z_auto > 1.0:
            print_dual(f"[INFO] Detected cell Z-height: {z_auto:.4f} Ang (used for N/m "
                       "conversion)", f_out)
            args.thickness = z_auto

    axis_index = {'x': 0, 'y': 1, 'z': 2}
    results_by_dir = {}
    try:
        for d, data in groups.items():
            cell_area = None
            if is_1d:
                if d[0] not in axis_index:
                    raise ValueError(
                        f"1D mode needs a recognized x/y/z axis to compute a cross-section; "
                        f"direction '{d}' isn't one.")
                # Read the cross-section from THIS direction's own first folder --
                # not a shared/arbitrary one, since a folder strained along a
                # different axis would give the wrong "other 2 vectors" here.
                ref_folder = _first_folder_for_direction(folders, d)
                outcell = siesta_log.get_outcell(os.path.join(ref_folder, args.file))
                if outcell is None:
                    raise ValueError(
                        "1D mode requires an 'outcell: Unit cell vectors' block in the SIESTA "
                        f"output (direction '{d}') to compute the cross-section area; not found.")
                cell_area = compute_cross_section_area(outcell, axis_index[d[0]])
            results_by_dir[d] = analyze_mechanics(data, d, is_2d, args.thickness,
                                                   is_1d, cell_area, args.cross_section)
    except ValueError as e:
        fail(str(e))

    results_by_dir = dict(sorted(results_by_dir.items()))
    u = next(iter(results_by_dir.values()))['unit']
    qty = next(iter(results_by_dir.values()))['quantity']

    print_section("[3] MECHANICAL PROPERTIES", f_out)
    if len(results_by_dir) > 1:
        print_dual(f"{len(results_by_dir)} direction(s) found -- comparing automatically "
                   "(no flag needed). Full detail per direction follows; quick-reference "
                   "table first:", f_out)
        print_table(["Direction", f"Slope ({u})", "R²", f"Peak ({u})", "Crit. Strain (%)",
                     f"Toughness ({'nN' if is_1d else ('J/m^2' if is_2d else 'GJ/m^3')})",
                     "Notes"], [
            ([d.upper(), f"{r['modulus']:.4f}", f"{r['r_squared']:.4f}", f"{r['uts']:.4f}",
              f"{r['uts_strain']:.4f}", f"{r['toughness']:.4f}",
              "edge of range" if r['uts_at_boundary'] else "--"],
             'yellow' if r['uts_at_boundary'] else None)
            for d, r in results_by_dir.items()
        ], f_out)
        print_dual("", f_out)
    for i, (d, r) in enumerate(results_by_dir.items()):
        if i > 0:
            print_dual("", f_out)
        print_mechanics_block(d, r, is_2d, is_1d, args.show_yield, f_out)

    print_section("[4] OUTPUT FILES", f_out)
    written_files = []
    if args.save_gnuplot:
        if len(results_by_dir) == 1:
            detected_dir = next(iter(results_by_dir))
            data = groups[detected_dir]
            results = results_by_dir[detected_dir]
            curve_dat = os.path.join(args.output_dir, f"{detected_dir}_curve.dat")
            curve_gplot = os.path.join(args.output_dir, f"{detected_dir}_curve.gplot")
            header = (f"{qty}-Strain Curve | Dir: {detected_dir} | Unit: {u}\n"
                      f"1:Strain(%) 2:Strain(Frac) 3:{qty}({u})")
            out_data = np.column_stack((data[:, 0], data[:, 1], results['axial_stress']))
            np.savetxt(curve_dat, out_data, fmt="%12.6f", header=header)
            write_strain_gplot(curve_dat, curve_gplot, detected_dir, results, u)
            written_files += [curve_dat, curve_gplot]
        else:
            compare_dat = os.path.join(args.output_dir, "comparison_curve.dat")
            compare_gplot = os.path.join(args.output_dir, "comparison_curve.gplot")
            write_compare_gplot(compare_dat, compare_gplot, results_by_dir, u)
            written_files += [compare_dat, compare_gplot]
        print_dual(color_text(
            f"[OK] Curve data + gnuplot script written: {len(written_files)} file(s).", 'green'), f_out)
        for path in written_files:
            print_dual(f"  - {path}", f_out)
    else:
        print_dual("Not written (off by default -- pass --save-gnuplot to write the curve "
                   "data + a companion .gplot script, renderable via 'gnuplot <script>').", f_out)

    print_section("[5] SUMMARY & NEXT STEPS", f_out)
    print_dual(color_text("[OK] Analysis complete.", 'green'), f_out)
    if report_path:
        print_dual(f"Report saved to   : {report_path}", f_out)
    else:
        print_dual("Pass --save-report to persist this report to a file.", f_out)
    if not args.save_gnuplot:
        print_dual("Pass --save-gnuplot to also write the curve data + a gnuplot script.", f_out)
    if not args.view:
        print_dual("Pass --view to see an interactive matplotlib plot before finishing.", f_out)

    if f_out:
        f_out.close()

    if args.view:
        import matplotlib.pyplot as plt
        # Fit line clipped to the actual fit window (see analyze_mechanics/
        # write_strain_gplot's own comment) -- never drawn across the full
        # strain range, which would visually overstate how far the
        # small-strain elastic fit is supposed to apply.
        if len(results_by_dir) == 1:
            detected_dir = next(iter(results_by_dir))
            results = results_by_dir[detected_dir]
            fig, ax = plt.subplots(figsize=(7, 5))
            ax.plot(results['strain_frac'], results['axial_stress'], 'o-', color='#2255cc',
                    label=f"{qty}-strain data")
            x_fit = np.array([results['fit_strain_lo'], results['fit_strain_hi']]) / 100.0
            ax.plot(x_fit, results['slope'] * x_fit + results['intercept'], '--',
                    color='#2255cc', label="Initial slope fit")
            uts_idx = results['uts_idx']
            ax.plot(results['strain_frac'][uts_idx], results['axial_stress'][uts_idx], 'o',
                    color='#cc2222', markersize=10, label="Peak")
            ax.set_title(f"{qty}-Strain Curve ({detected_dir.upper()})")
        else:
            colors = ['#2255cc', '#cc5522', '#22aa55', '#aa22aa', '#aaaa22', '#22aaaa']
            fig, ax = plt.subplots(figsize=(7, 5))
            for i, (d, results) in enumerate(results_by_dir.items()):
                color = colors[i % len(colors)]
                ax.plot(results['strain_frac'], results['axial_stress'], 'o-', color=color,
                        label=d.upper())
                x_fit = np.array([results['fit_strain_lo'], results['fit_strain_hi']]) / 100.0
                ax.plot(x_fit, results['slope'] * x_fit + results['intercept'], '--', color=color)
                uts_idx = results['uts_idx']
                ax.plot(results['strain_frac'][uts_idx], results['axial_stress'][uts_idx], 'o',
                        color=color, markersize=10)
            ax.set_title(f"{qty}-Strain Comparison")
        ax.set_xlabel("Strain (fraction)")
        ax.set_ylabel(f"{qty} ({u})")
        ax.grid(True, alpha=0.3)
        ax.legend()
        fig.tight_layout()
        plt.show()

if __name__ == "__main__":
    main()
