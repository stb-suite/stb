#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Elastic Generator       #
#################################################

VERSION = "1.9.1"

import os
import sys
import shutil
import argparse
import numpy as np
from time import sleep
from stb.core import structure_io, kspace, symmetry
from stb.core.cli import COLORS, color_text, show_intro

# ==========================================
#           HELPER FUNCTIONS
# ==========================================

def get_strain_matrix(direction, delta):
    """Returns the deformation matrix (I + epsilon)."""
    d = direction.lower()

    # Combined energy-method pattern (e.g. 'xx+yy'): sum each submode's own
    # strain tensor -- recursing here (instead of a separate code path)
    # guarantees a combined pattern's deformation can never drift from what
    # its pure submodes mean individually.
    if '+' in d:
        epsilon = sum(get_strain_matrix(part, delta) - np.eye(3) for part in d.split('+'))
        return np.eye(3) + epsilon

    # Robust mapping (accepts x or xx, zx or xz) -- the 6 canonical Voigt
    # modes are built from core/symmetry.py's single source of truth (also
    # used by equivalent_strain_modes/find_mapping_operation), so the
    # deformation this tool actually generates can never drift from what
    # the symmetry analysis assumes it means.
    canonical = {'x': 'xx', 'y': 'yy', 'z': 'zz', 'zy': 'yz', 'zx': 'xz', 'yx': 'xy'}
    voigt_mode = canonical.get(d, d)
    if voigt_mode in symmetry.VOIGT_MODES:
        epsilon = symmetry.strain_tensor(voigt_mode, delta)
    else:
        epsilon = np.zeros((3, 3))
        if d == 'bi':
            epsilon[0, 0] = epsilon[1, 1] = delta
        elif d == 'hydro':
            epsilon[0, 0] = epsilon[1, 1] = epsilon[2, 2] = delta

    return np.eye(3) + epsilon


def direction_axes(direction):
    """Which Cartesian axis indices (0=x, 1=y, 2=z) a --dirs strain mode
    touches -- same direction vocabulary as get_strain_matrix, used to check
    each requested direction against the structure's vacuum-padded axes
    before generating anything.
    """
    d = direction.lower()
    if d in ('x', 'xx'): return {0}
    if d in ('y', 'yy'): return {1}
    if d in ('z', 'zz'): return {2}
    if d in ('yz', 'zy'): return {1, 2}
    if d in ('xz', 'zx'): return {0, 2}
    if d in ('xy', 'yx'): return {0, 1}
    if d == 'bi': return {0, 1}
    if d == 'hydro': return {0, 1, 2}
    return set()


def generate_verify_script():
    """Generates a small helper script to check calculation status."""
    content = r"""#!/bin/bash
for d in strain_*/ ; do
    [ -d "$d" ] && echo "Checking $d..." && tail -n 1 "$d/calc.out" 2>/dev/null
done
"""
    with open("verify_calc.sh", "w") as f:
        f.write(content)
    try:
        os.chmod("verify_calc.sh", 0o755)
    except OSError:
        pass # Ignore if permissions cannot be changed (e.g., Windows)

# ==========================================
#                  MAIN
# ==========================================

def main():
    desc_text = "Siesta Elastic Strain Generator (Batch Mode)"
    # Detailed help description
    help_dirs = (
        "List of directions to generate. \n"
        "Options: [xx, yy, zz, xy, xz, yz, bi, hydro, all]. \n"
        "Example: --dirs xx yy xy\n"
        "'all' auto-reduces to 1 representative direction per symmetry-\n"
        "equivalent group (e.g. just 2 of the 6 for a cubic material) --\n"
        "stb-elasticAnalysis reconstructs the rest by symmetry. Pass every\n"
        "direction explicitly (not 'all') to run all 6 anyway, e.g. to get\n"
        "stb-elasticAnalysis's diagnostic comparison across truly independent\n"
        "DFT calculations instead."
    )

    parser = argparse.ArgumentParser(description=desc_text, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("--file", "-i", required=True, help="Input structural FDF file")
    parser.add_argument("--method", choices=["stress", "energy"], default="stress",
                        help="Physical quantity the generated strains are meant to feed "
                             "(default: stress). 'stress': single-component canonical Voigt "
                             "strains, for stb-elasticAnalysis's stress-strain fit (--method "
                             "stress there too). 'energy': a point-group-selected set of pure "
                             "AND combined (two Voigt components at once, e.g. 'xx+yy') strain "
                             "patterns for an energy-strain (parabolic) fit -- combined patterns "
                             "are physically required for off-diagonal constants there, since a "
                             "pure single-component strain's energy curvature only ever "
                             "determines its own diagonal constant. Always auto-selects (--dirs "
                             "must be left at its 'all' default); needs more DFT calculations "
                             "than --method stress for the same crystal system, since each "
                             "pattern only yields one number (a curvature) instead of a full "
                             "stress column.")
    parser.add_argument("--dirs", nargs='+', default=["all"], help=help_dirs)
    parser.add_argument("--symmetry-method", choices=["basic", "full"], default="basic",
                        help="How 'all' picks which directions to actually run (default: "
                             "basic). 'basic': only reduces when a single point-group operation "
                             "maps one canonical direction exactly onto another (e.g. cubic); "
                             "misses reductions like hexagonal's C11=C22, which come from the "
                             "full elastic tensor's symmetry, not a pairwise direction match. "
                             "'full': picks the minimal set of directions that fully determines "
                             "the point group's independent elastic constants (greedy rank "
                             "selection) -- catches hexagonal/trigonal reductions too. Must "
                             "match the --symmetry-method passed to stb-elasticAnalysis later.")
    parser.add_argument("--max", type=float, default=2.0, help="Max strain %%")
    parser.add_argument("--steps", type=int, default=4, help="Steps per direction")
    parser.add_argument("--output", default="structure.fdf", help="Output filename")
    parser.add_argument("--vacuum-gap", type=float, default=10.0,
                        help="Vacuum gap threshold in Ang used to detect which lattice axes are "
                             "periodic vs. vacuum-padded (default: 10.0), same convention as "
                             "stb-kgrid/stb-strain. Requested directions touching a vacuum-padded "
                             "axis are skipped with a warning instead of generated.")
    parser.add_argument("--symprec", type=float, default=1e-3,
                        help="Symmetry-detection tolerance (default: 1e-3, pymatgen's own default), "
                             "used only for the informational note about triclinic/monoclinic "
                             "structures (this method assumes no normal-shear strain coupling).")
    parser.add_argument("--angle-tolerance", type=float, default=5.0,
                        help="Symmetry angle tolerance in degrees (default: 5.0, pymatgen's own "
                             "default), same use as --symprec.")

    parser.add_argument("-v", "--version", action="version",
                        version=f"stb-elasticInputs {VERSION}")
    # --- NEW: Argument --no-intro ---
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.method == "energy" and args.dirs != ["all"]:
        sys.exit(color_text(
            "[ERROR] --method energy always auto-selects its own strain patterns (pure and "
            "combined) from the structure's point group -- an explicit --dirs list isn't "
            "supported yet. Omit --dirs (or pass --dirs all, the default) to use it.", 'red'))

    # Shows intro if --no-intro flag is NOT used
    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - Elastic Generator",
            "Strain Structure Generator (Batch Mode)",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Integrated Style Refactoring"
        ])

    # Expand 'all' keyword. Only 'all' is eligible for the symmetry-based
    # reduction below -- an explicit list is respected exactly as given
    # (the way to deliberately run redundant directions and get
    # stb-elasticAnalysis's diagnostic/pooling across them).
    used_all = 'all' in args.dirs
    dirs_to_run = []
    if used_all:
        dirs_to_run = ['xx', 'yy', 'zz', 'xy', 'xz', 'yz']
    else:
        dirs_to_run = args.dirs

    try:
        structure = structure_io.read_fdf(args.file)
        lattice = structure_io.raw_lattice_vectors(structure)
        print(f"{color_text('[INFO]', 'green')} Loaded: {args.file}")
    except Exception as e:
        sys.exit(f"{color_text('[ERROR]', 'red')} {e}")

    # Save an UNDEFORMED copy for stb-elasticAnalysis's later symmetry
    # detection. Deliberately NOT read from inside a strain_*/ folder for
    # that: any single-axis normal strain (e.g. 'xx' alone) genuinely lowers
    # the crystal's symmetry by construction (a cubic cell strained along x
    # only has a != b == c, exactly tetragonal, not "cubic with noise") --
    # verified live catching this exact mistake (detected 4/mmm instead of
    # m-3m for a strained folder). This plain, un-strained copy is the only
    # correct reference.
    try:
        shutil.copy(args.file, "reference_structure.fdf")
    except OSError:
        pass  # best-effort; stb-elasticAnalysis just skips symmetry detection without it

    # Detect which axes are periodic vs. vacuum-padded, same convention as
    # stb-kgrid/stb-kpath/stb-strain -- straining a vacuum gap just moves
    # empty space, not the material.
    positions = np.array([pos for _, pos in structure.atoms])
    is_cartesian = structure.coord_format == 'cartesian'
    frac_coords = kspace.to_fractional(positions, structure.lattice, is_cartesian)
    vacuum_axes = kspace.detect_vacuum_axes(frac_coords, structure.lattice, args.vacuum_gap)
    print(f"{color_text('[INFO]', 'green')} Detected dimensionality: {kspace.dimensionality_label(vacuum_axes)}")

    axis_letters = ['x', 'y', 'z']
    periodic_axes = {i for i in range(3) if not vacuum_axes[i]}

    if args.method == "stress":
        # Unlike stb-strain (a single requested direction, block on any
        # vacuum hit), --dirs is a list -- drop only the physically
        # meaningless entries and keep going with the rest, so 'all' on a
        # 2D/1D structure still does something useful instead of failing
        # outright.
        filtered_dirs = []
        for d in dirs_to_run:
            touched = direction_axes(d)
            blocked = sorted(touched - periodic_axes)
            if blocked:
                letters = ', '.join(axis_letters[i] for i in blocked)
                print(f"{color_text('[WARNING]', 'yellow')} Skipping direction '{d}': touches "
                      f"vacuum-padded axis/axes {letters} (detected via --vacuum-gap "
                      f"{args.vacuum_gap} Ang); straining a vacuum gap doesn't correspond to any "
                      "physical deformation.")
                continue
            filtered_dirs.append(d)
        dirs_to_run = filtered_dirs

        if not dirs_to_run:
            periodic_str = ', '.join(axis_letters[i] for i in sorted(periodic_axes)) or \
                'none (structure is fully isolated / 0D)'
            sys.exit(color_text(
                "[ERROR] None of the requested directions are physically valid for this structure "
                f"-- periodic axis/axes available: {periodic_str}.", 'red'))

        # Symmetry-based reduction: only for the 'all' convenience default (an
        # explicit --dirs list is never second-guessed). Keeps 1 representative
        # per group of symmetry-equivalent directions -- stb-elasticAnalysis
        # reconstructs the skipped ones' stiffness columns exactly via the same
        # point-group operation (see core/symmetry.py::find_mapping_operation),
        # so this is a real DFT-cost cut, not an approximation.
        if used_all and args.symmetry_method == "basic":
            try:
                pmg_structure = structure_io.to_pymatgen(structure)
                groups, point_group, _ops = symmetry.equivalent_strain_modes(
                    pmg_structure, args.symprec, args.angle_tolerance)
                kept, seen = [], set()
                for d in dirs_to_run:
                    if d in seen:
                        continue
                    group = next(g for g in groups if d in g)
                    survivors = [m for m in dirs_to_run if m in group]
                    representative = survivors[0]
                    kept.append(representative)
                    skipped = survivors[1:]
                    if skipped:
                        print(f"{color_text('[INFO]', 'cyan')} Skipping {skipped} -- symmetry-equivalent "
                              f"to '{representative}' (point group {point_group}); stb-elasticAnalysis "
                              f"will derive their stiffness columns from the '{representative}' "
                              "calculation instead of running DFT again.")
                    seen.update(group)
                dirs_to_run = kept
            except Exception:
                pass  # optimization only, never blocks -- falls back to running everything
        elif used_all and args.symmetry_method == "full":
            try:
                pmg_structure = structure_io.to_pymatgen(structure)
                basis, point_group = symmetry.symmetry_allowed_basis(
                    pmg_structure, args.symprec, args.angle_tolerance)
                kept = symmetry.select_directions_by_rank(basis, order=dirs_to_run)
                skipped = [d for d in dirs_to_run if d not in kept]
                if skipped:
                    print(f"{color_text('[INFO]', 'cyan')} Point group {point_group} has "
                          f"{len(basis)} independent elastic constant(s); {kept} fully determine "
                          f"them (--symmetry-method full). Skipping {skipped} -- "
                          "stb-elasticAnalysis --symmetry-method full will reconstruct the "
                          "remaining stiffness-matrix entries from the symmetry-allowed fit "
                          "instead of running DFT again.")
                dirs_to_run = kept
            except Exception:
                pass  # optimization only, never blocks -- falls back to running everything

        # Informational note when the structure has genuine normal-shear
        # coupling constants (C14-C16, C24-C26, C34-C36, C45, C46, C56) --
        # stb-elasticAnalysis recovers all of them correctly from the full 6
        # -direction data (see its own docstring), but only if all 6 directions
        # are actually run; this is a reminder for that case, not a limitation
        # warning (see core/symmetry.py::crystal_system docstring for why only
        # these 2 systems are flagged specifically -- they're the only ones
        # where every subclass is guaranteed to have nonzero coupling).
        try:
            pmg_structure = structure_io.to_pymatgen(structure)
            system, point_group = symmetry.crystal_system(pmg_structure, args.symprec, args.angle_tolerance)
            if system in ('triclinic', 'monoclinic'):
                print(color_text(
                    f"[INFO] Detected {system} symmetry (point group {point_group}) -- this crystal "
                    "class always has some nonzero normal-shear coupling constants (C14-C16, C24-C26, "
                    "C34-C36, C45, C46, C56). stb-elasticAnalysis recovers them correctly, but only "
                    "from directions that were actually run: make sure to keep all 6 default "
                    "directions (--dirs all, the default) rather than a partial --dirs list, or those "
                    "coupling terms will read as 0 from missing data instead of their real value.",
                    'cyan'))
        except Exception:
            pass  # symmetry detection is advisory only, never blocks the tool

    else:
        # --method energy: dirs_to_run is currently the 6 canonical modes
        # (from the 'all' expansion above, --dirs is forced to 'all' by the
        # early check) -- discard it and pick from the full 21-pattern pool
        # (6 pure + 15 combined) instead, filtered against the SAME vacuum
        # -padded axes, then reduced by symmetry the same way as --method
        # stress's --symmetry-method full, just with a 1-equation-per
        # -pattern rank criterion instead of 6-per-direction (see
        # core/symmetry.py's energy-method module note for why).
        candidates = symmetry.energy_pattern_candidates()
        filtered_candidates = []
        for name, eps_t in candidates:
            touched = set().union(*(direction_axes(part) for part in name.split('+')))
            blocked = sorted(touched - periodic_axes)
            if blocked:
                continue
            filtered_candidates.append((name, eps_t))

        if not filtered_candidates:
            periodic_str = ', '.join(axis_letters[i] for i in sorted(periodic_axes)) or \
                'none (structure is fully isolated / 0D)'
            sys.exit(color_text(
                "[ERROR] No strain pattern is physically valid for this structure -- periodic "
                f"axis/axes available: {periodic_str}.", 'red'))

        try:
            pmg_structure = structure_io.to_pymatgen(structure)
            basis, point_group = symmetry.symmetry_allowed_basis(
                pmg_structure, args.symprec, args.angle_tolerance)
            kept = symmetry.select_energy_patterns_by_rank(basis, filtered_candidates)
            dirs_to_run = [name for name, _ in kept]
            achieved = len(kept)
            print(f"{color_text('[INFO]', 'cyan')} Point group {point_group} has "
                  f"{len(basis)} independent elastic constant(s); {achieved} strain pattern(s) "
                  f"needed to determine {'all of them' if achieved >= len(basis) else 'the ones reachable without leaving the periodic plane/axis'} "
                  f"(of {len(filtered_candidates)} physically valid candidates, {len(candidates)} in total): {dirs_to_run}.")
        except Exception as e:
            sys.exit(color_text(
                f"[ERROR] --method energy requires symmetry detection to pick its strain "
                f"patterns and it failed: {e}", 'red'))

    print(f"{color_text('[INFO]', 'green')} Modes: {dirs_to_run}")

    count = 0
    strains = np.linspace(-args.max, args.max, args.steps)

    for d in dirs_to_run:
        for s in strains:
            if abs(s) < 1e-9: s = 0.0

            delta = s / 100.0
            def_matrix = get_strain_matrix(d, delta)
            new_lattice = np.dot(lattice, def_matrix)

            # Create folder: e.g. strain_xx_-1.0
            folder = f"strain_{d}_{'m' if s < 0 else ''}{abs(s):.2f}"
            if not os.path.exists(folder): os.makedirs(folder)

            structure_io.rewrite_fdf_lattice(args.file, new_lattice, os.path.join(folder, args.output))
            count += 1

    print(f"\n{color_text('[SUCCESS]', 'green')} Generated {count} structures.")
    generate_verify_script()
    print(f"Run calculations inside each 'strain_...' folder.")

if __name__ == "__main__":
    main()
