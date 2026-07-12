#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
#     Mechanical Properties Analyzer            #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.9.1"

import os
import sys
import argparse
import numpy as np
from scipy.stats import linregress
# Use trapezoid for compatibility with new Scipy
from scipy.integrate import trapezoid
from stb.core import siesta_log
from stb.core.cli import COLORS, color_text, show_intro

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
    convert kBar -> axial Force (nN) for --1d: SIESTA's stress is already
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
    - 2D (--2d): Stress (N/m), corrected by the vacuum (Z) height.
    - 1D (--1d): Force (nN) -- NOT stress. kBar is already normalized by the
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
    # Linear fit only in small strain region (-2% to 2%)
    limit_idx = np.where(abs(strain_pct) <= 2.0)
    if len(limit_idx[0]) < 3: limit_idx = range(len(strain_pct))

    slope, intercept, r_val, _, _ = linregress(strain_frac[limit_idx], sigma_axial[limit_idx])

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


def _canonical_direction(direction):
    """Folds the doubled-letter uniaxial convention ('xx'/'yy'/'zz', used by
    some older fixtures) into the canonical single-letter form ('x'/'y'/'z'),
    so both spellings of the same physical direction group together instead
    of being treated as 2 different directions (which could otherwise
    spuriously trigger the "multiple directions found" error, or show up as
    2 separate rows in --compare for what is really 1 direction). Genuine
    biaxial pairs ('xy'/'xz'/'yz') are returned unchanged.
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
    matches `direction` -- used to read the reference cell for --1d's
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


def _dimensionality_label(is_2d, is_1d):
    if is_1d:
        return "1D (Wire/Tube)"
    if is_2d:
        return "2D (Sheet)"
    return "3D (Bulk)"


def build_report(direction, results, is_2d, is_1d=False, show_yield=False):
    """Returns (screen_lines, file_lines) describing one direction's
    mechanical properties -- shared by the single-direction and --compare
    report paths so the two stay in sync.
    """
    u = results['unit']
    qty = results['quantity']
    modulus_label = f"Initial Slope (d{'F' if is_1d else chr(963)}/d{chr(949)})"
    modulus_str = f"{results['modulus']:.2f} {u}"
    peak_str = f"{results['uts']:.2f} {u}"
    peak_label = f"Peak {qty}"

    screen = [
        f"Direction       : {direction.upper()}",
        f"Dimensionality  : {_dimensionality_label(is_2d, is_1d)}",
        "-" * 50,
        f"{modulus_label:<16}: {color_text(modulus_str, 'green')} (R²={results['r_squared']:.4f})",
        f"{peak_label:<16}: {color_text(peak_str, 'red')}",
        f"Critical Strain : {results['uts_strain']:.2f} %",
    ]
    if not is_1d:
        screen.append(color_text(
            "[INFO] Single-direction slope under a clamped transverse cell -- for the "
            "rigorous small-strain elastic-tensor Young's Modulus, use stb-elasticAnalysis.",
            'cyan'))
    if results['uts_at_boundary']:
        screen.append(color_text(
            "[WARNING] Peak occurred at the edge of the tested strain range -- "
            "the true peak may lie beyond --stmin/--stmax; treat it as a lower bound.",
            'yellow'))
    if is_1d and 'conventional_stress' in results:
        conv_peak = results['conventional_stress'][results['uts_idx']]
        screen.append(f"Peak Stress (conventional, cross-section={results['cross_section']:.2f} "
                       f"Ang^2): {conv_peak:.2f} GPa")

    toughness_unit = "nN (work/length, ≡ J/m)" if is_1d else ('J/m^2' if is_2d else 'GJ/m^3')
    screen.append(f"Toughness       : {results['toughness']:.4f} {toughness_unit}")

    if show_yield:
        if results['yield'] > 0:
            screen.append(f"Yield {qty} (0.2%): {results['yield']:.2f} {u}")
        else:
            screen.append(f"Yield {qty} (0.2%): {color_text('Not detected (Linear)', 'yellow')}")
        screen.append(color_text(
            "[INFO] '0.2% offset yield' is a macroscopic-plasticity concept borrowed from "
            "metallurgy (dislocation-glide onset) -- a defect-free periodic crystal under "
            "affine strain has no such mechanism; interpret with care.", 'cyan'))

    if not is_1d and _is_uniaxial(direction):
        transverse = _transverse_stress_lines(results)
        if transverse:
            screen.append("Transverse " + qty + " @ peak (diagnostic, clamped cell): " +
                           ", ".join(transverse))
            screen.append(color_text(
                "[INFO] Not a Poisson's ratio (that needs transverse STRAIN, not stress) -- "
                "just how much elastic response this direction's clamped transverse cell is "
                "suppressing.", 'cyan'))

    file_lines = [
        f"Direction: {direction}",
        f"{modulus_label}: {results['modulus']:.4f} {u}",
        f"{peak_label}: {results['uts']:.4f} {u}",
        f"Critical Strain: {results['uts_strain']:.4f} %",
    ]
    if results['uts_at_boundary']:
        file_lines.append(
            "WARNING: Peak occurred at the edge of the tested strain range -- "
            "the true peak may lie beyond --stmin/--stmax; treat it as a lower bound.")
    if is_1d and 'conventional_stress' in results:
        conv_peak = results['conventional_stress'][results['uts_idx']]
        file_lines.append(f"Peak Stress (conventional, cross-section={results['cross_section']:.4f} "
                           f"Ang^2): {conv_peak:.4f} GPa")
    if show_yield:
        file_lines.append(f"Yield {qty} (0.2%): {results['yield']:.4f} {u}")
    file_lines.append(f"Toughness: {results['toughness']:.4f} {toughness_unit}")
    if not is_1d and _is_uniaxial(direction):
        transverse = _transverse_stress_lines(results)
        if transverse:
            file_lines.append("Transverse " + qty + " @ peak (diagnostic): " + ", ".join(transverse))

    return screen, file_lines


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
        lines.append(f"{axis_names[idx]}={val:.2f} {results['unit']}")
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
        (f'plot "{curve_filename}" index 0 using 2:3 with points pt 7 ps 1.2 lc rgb "#2255cc" '
         f'title "{qty}-strain data", \\\n'
         '     f(x) with lines lc rgb "#2255cc" dt 2 title "Initial slope fit", \\\n'
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
        '# Generated by stb-strainAnalysis --compare\n',
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
        lines.append(f"f_{i}(x) = {results['slope']:.8f}*x + {results['intercept']:.8f}\n")
        plot_parts.append(
            f'"{data_filename}" index {idx0} using 2:3 with points pt 7 ps 1.2 lc rgb "{color}" '
            f'title "{direction.upper()}"')
        plot_parts.append(f'f_{i}(x) with lines lc rgb "{color}" dt 2 notitle')
        plot_parts.append(
            f'"{data_filename}" index {idx1} using 2:3 with points pt 9 ps 2 lc rgb "{color}" notitle')
    lines.append('plot ' + ', \\\n     '.join(plot_parts) + '\n')

    with open(gplot_filename, 'w') as f:
        f.writelines(lines)

# ==========================================
#           MAIN
# ==========================================

def main():
    parser = argparse.ArgumentParser(description="Mechanical Properties from Stress-Strain Data.")
    parser.add_argument("-f", "--file", required=True, help="Siesta output file (e.g., calc.out).")
    parser.add_argument("-o", "--output", default="mechanical_curve.dat", help="Output raw data file.")
    parser.add_argument("--dir", default="strain_runs",
                        help="Directory containing the 'strain_*' run folders (default: strain_runs), "
                             "same convention as stb-strain's --output-dir.")
    parser.add_argument("--2d", dest="is2d", action="store_true", help="Enable 2D units (N/m).")
    parser.add_argument("--thickness", type=float, default=20.0, help="Vacuum height (Z) for 2D conversion (Angstrom).")
    parser.add_argument("--1d", dest="is1d", action="store_true",
                        help="Enable 1D mode: report axial Force (nN) instead of Stress, using the "
                             "cell's own cross-section (auto-detected, vacuum-agnostic) -- for "
                             "wires/nanotubes (stb-nanotube output).")
    parser.add_argument("--cross-section", type=float, default=None,
                        help="Physical cross-section area (Ang^2) of the wire/tube wall, only used "
                             "with --1d, to additionally report a conventional Stress (GPa) -- there "
                             "is no way to auto-derive this (unlike the 2D vacuum height), since a "
                             "wire/tube wall thickness is a physical convention, not a cell property.")
    parser.add_argument("--yield", dest="show_yield", action="store_true",
                        help="Also report the 0.2%% offset Yield -- a macroscopic-plasticity concept "
                             "borrowed from metallurgy; off by default since a defect-free periodic "
                             "crystal under affine DFT strain has no dislocation-nucleation mechanism "
                             "to actually yield through.")
    parser.add_argument("--compare", action="store_true",
                        help="Analyze every strain direction found (e.g. 'strain_x_*' and "
                             "'strain_y_*' side by side) and print/save a comparison table "
                             "instead of requiring a single direction.")
    parser.add_argument("-v", "--version", action="version",
                        version=f"stb-strainAnalysis {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.is1d and args.is2d:
        sys.exit(color_text("[ERROR] --1d and --2d are mutually exclusive.", 'red'))
    if args.cross_section is not None and not args.is1d:
        sys.exit(color_text("[ERROR] --cross-section only applies together with --1d.", 'red'))

    if args.intro:
        show_intro([
            f"Siesta Tool Box - Mechanical Analysis v{VERSION}",
            "Developed by Dr. Carlos M. O. Bastos",
            "bastoscmo.github.io"
        ])

    print(color_text("-> Scanning directories...", 'green'))

    if not os.path.isdir(args.dir):
        sys.exit(color_text(f"[ERROR] Directory '{args.dir}' not found.", 'red'))

    folders = [os.path.join(args.dir, d) for d in os.listdir(args.dir)
               if os.path.isdir(os.path.join(args.dir, d)) and d.startswith('strain_')]
    if not folders:
        sys.exit(color_text(f"[ERROR] No 'strain_*' folders found in '{args.dir}'.", 'red'))

    # Detect Z-height automatically if 2D
    if args.is2d and os.path.exists(os.path.join(folders[0], args.file)):
        z_auto = siesta_log.get_cell_height(os.path.join(folders[0], args.file))
        if z_auto > 1.0:
            print(f"{color_text('[INFO]', 'cyan')} Detected Cell Z-Height: {z_auto:.2f} Ang (used for N/m conversion)")
            args.thickness = z_auto

    print(f"   Found {len(folders)} strain steps.")

    groups = collect_strain_groups(args.file, folders)
    if not groups:
        sys.exit(color_text("[ERROR] No stress data found.", 'red'))

    if len(groups) > 1 and not args.compare:
        dirs_str = ", ".join(sorted(groups))
        sys.exit(color_text(
            f"[ERROR] Multiple strain directions found ({dirs_str}) in '{args.dir}'; "
            "pass --compare to analyze all of them, or run from a directory with only one "
            "direction's 'strain_*' folders.", 'red'))

    for direction, data in groups.items():
        if not _is_uniaxial(direction):
            print(f"{color_text('[WARNING]', 'yellow')} Biaxial direction '{direction}' detected -- "
                  f"reporting only the stress component along its first axis "
                  f"('{direction[0]}{direction[0]}'), not a true biaxial modulus.")

    axis_index = {'x': 0, 'y': 1, 'z': 2}
    results_by_dir = {}
    try:
        for d, data in groups.items():
            cell_area = None
            if args.is1d:
                if d[0] not in axis_index:
                    raise ValueError(
                        f"--1d needs a recognized x/y/z axis to compute a cross-section; "
                        f"direction '{d}' isn't one.")
                # Read the cross-section from THIS direction's own first folder --
                # not a shared/arbitrary one, since a folder strained along a
                # different axis would give the wrong "other 2 vectors" here.
                ref_folder = _first_folder_for_direction(folders, d)
                outcell = siesta_log.get_outcell(os.path.join(ref_folder, args.file))
                if outcell is None:
                    raise ValueError(
                        "--1d requires an 'outcell: Unit cell vectors' block in the SIESTA "
                        f"output (direction '{d}') to compute the cross-section area; not found.")
                cell_area = compute_cross_section_area(outcell, axis_index[d[0]])
            results_by_dir[d] = analyze_mechanics(data, d, args.is2d, args.thickness,
                                                   args.is1d, cell_area, args.cross_section)
    except ValueError as e:
        sys.exit(color_text(f"[ERROR] {e}", 'red'))

    if args.compare:
        # --- COMPARISON TABLE ---
        u = next(iter(results_by_dir.values()))['unit']
        qty = next(iter(results_by_dir.values()))['quantity']
        print("\n" + "="*70)
        print(color_text(f"   MECHANICAL {qty.upper()} COMPARISON", 'bold').center(70))
        print("="*70)
        header = f"{'Direction':<10}{'Slope (' + u + ')':<18}{'Peak (' + u + ')':<16}{'Crit. Strain %':<16}"
        print(header)
        print("-"*70)
        report_lines = [header, "-"*70]
        for direction, results in results_by_dir.items():
            row = (f"{direction.upper():<10}{results['modulus']:<18.2f}"
                   f"{results['uts']:<16.2f}{results['uts_strain']:<16.2f}")
            flag = " (edge of range)" if results['uts_at_boundary'] else ""
            print(row + flag)
            report_lines.append(row + flag)
        print("="*70)

        comparison_file = "mechanical_comparison.txt"
        with open(comparison_file, "w") as f:
            f.write("========================================\n")
            f.write(f"   MECHANICAL {qty.upper()} COMPARISON\n")
            f.write("========================================\n")
            f.write("\n".join(report_lines) + "\n")
        print(f"\n{color_text('[Saved]', 'cyan')} Comparison -> {comparison_file}")

        gplot_base = os.path.splitext(args.output)[0]
        compare_dat = f"{gplot_base}_comparison.dat"
        compare_gplot = f"{gplot_base}_comparison.gplot"
        write_compare_gplot(compare_dat, compare_gplot, results_by_dir, u)
        print(f"{color_text('[Saved]', 'cyan')} Plot data  -> {compare_dat}")
        print(f"{color_text('[Saved]', 'cyan')} Plot script-> {compare_gplot} (gnuplot {compare_gplot})\n")
        return

    # --- SINGLE-DIRECTION REPORT (default path) ---
    detected_dir, data = next(iter(groups.items()))
    results = results_by_dir[detected_dir]
    u = results['unit']
    qty = results['quantity']

    print("\n" + "="*50)
    print(color_text(f"   MECHANICAL {qty.upper()} REPORT", 'bold').center(50))
    print("="*50)
    screen_lines, file_lines = build_report(detected_dir, results, args.is2d, args.is1d, args.show_yield)
    for line in screen_lines:
        print(line)
    print("="*50)

    # --- SAVE FILES ---

    # 1. Curve Data
    header = (f"{qty}-Strain Curve | Dir: {detected_dir} | Unit: {u}\n"
              f"1:Strain(%) 2:Strain(Frac) 3:{qty}({u})")

    out_data = np.column_stack((data[:, 0], data[:, 1], results['axial_stress']))
    np.savetxt(args.output, out_data, fmt="%12.6f", header=header)
    print(f"\n{color_text('[Saved]', 'cyan')} Curve data -> {args.output}")

    # 2. Report Text
    report_file = "mechanical_report.txt"
    with open(report_file, "w") as f:
        f.write("========================================\n")
        f.write(f"      MECHANICAL {qty.upper()} REPORT      \n")
        f.write("========================================\n")
        f.write("\n".join(file_lines) + "\n")
    print(f"{color_text('[Saved]', 'cyan')} Summary    -> {report_file}")

    # 3. Plot (data + gnuplot script)
    gplot_base = os.path.splitext(args.output)[0]
    gplot_file = f"{gplot_base}.gplot"
    write_strain_gplot(args.output, gplot_file, detected_dir, results, u)
    print(f"{color_text('[Saved]', 'cyan')} Plot script-> {gplot_file} (gnuplot {gplot_file})\n")

if __name__ == "__main__":
    main()
