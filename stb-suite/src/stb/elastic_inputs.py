#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Elastic Generator       #
#################################################

"""Stage 1 of 2: applies one or more canonical Voigt strains (xx/yy/zz
normal, xy/xz/yz shear, or a --method energy combined pattern like
'xx+yy') to a SIESTA structure and writes one ready-to-run folder per
scanned strain value -- structure + calc.fdf + config_extra.fdf +
pseudopotentials -- mirroring stb-strain's own split-file convention (see
its own module docstring for the general rationale), generalized to
however many directions --dirs actually requests in one run (up to all 6
canonical directions, or 'all''s symmetry-reduced subset -- unlike
stb-strain, which only ever strains one direction per invocation).
Doesn't run SIESTA itself -- run each folder yourself, then use
stb-elasticAnalysis.

Output layout: one <output-dir>/<direction>/strain_<direction>_<pct>/
folder per (direction, strain value), same convention as stb-strain's
<output-dir>/<direction>/... -- but since a single stb-elasticInputs
invocation can generate MULTIPLE directions at once, the report is a
single, whole-run report at <output-dir>/<REPORT_FILE> (not one per
direction), and the undeformed reference_structure.fdf lives at the TOP
of <output-dir> -- exactly what lets 'cd <output-dir> && stb-elasticAnalysis'
keep working with no extra flags needed on the Stage 2 side.

--calc's role, and why its override %include is PREPENDED, not appended:
like stb-strain, --calc is a reusable template copied into every
generated folder (with one %include line added). But where stb-strain's
own include (a %block Geometry.Constraints) is a block name that never
pre-exists in a real calc.fdf, this tool needs to force a pure
SINGLE-POINT SCF (no ionic or cell relaxation at all -- there is no
--relax-mode here, every generated folder is single-point, since the
imposed strain itself must stay exactly at the sampled geometry) by
overriding MD.TypeOfRun/MD.Steps/MD.VariableCell --
directives that ARE commonly already set in a real relaxation calc.fdf
(e.g. this repo's own examples/4.1-strain/calc.fdf has all three). SIESTA's
fdf reader is first-occurrence-wins for duplicate labels (verified in
hubbardu.py::write_run_folder, against real SIESTA 5.4.2 runs) -- so the
override %include is prepended at the very top of the generated calc.fdf,
before the base template's own directives, guaranteeing it wins regardless
of where in the base file those directives happen to appear. Inserting
it AFTER the structure %include instead (as stb-strain's own
insert_include_after_structure does for Geometry.Constraints) would be
silently ignored whenever --calc already sets any of these -- the common
case for a real relaxation template.
"""

VERSION = "2.0.0"

import os
import sys
import shutil
import argparse
import numpy as np
from stb.core import structure_io, kspace, symmetry
from stb.core.cli import color_text, show_intro, print_dual, print_section, print_table
from stb.core.pseudopotentials import BANKS, resolve_pseudo_source, get_required_pseudos, copy_pseudo
from stb.core.structure_io import read_md_state, is_siesta_true, prepend_include

REPORT_FILE = "elastic_stage1.txt"
EXTRA_FDF_FILE = "config_extra.fdf"
REFERENCE_FILE = "reference_structure.fdf"

# ==========================================
#           HELPER FUNCTIONS
# ==========================================

def get_strain_matrix(direction, delta):
    """Returns the deformation matrix (I + epsilon).

    NOTE: mlelastic.py imports this function directly (module-level) --
    keep its name/signature stable.
    """
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

    NOTE: mlelastic.py imports this function directly (module-level) --
    keep its name/signature stable.
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


def print_symmetry_table(dirs_before, kept, point_group, ops, symmetry_method, mapping=None, f_out=None):
    """Pretty table of which of the 6 canonical directions were kept vs.
    suppressed by symmetry reduction, plus a reduced summary of the detected
    point group and its operations (see core/symmetry.py::operations_summary).

    `mapping` ({suppressed_direction: representative_direction}) is given for
    'basic' (an exact pairwise point-group-operation correspondence exists --
    see core/symmetry.py::find_mapping_operation); left None for 'full' (a
    joint least-squares fit over the symmetry-allowed subspace -- there's no
    single representative direction a suppressed one maps onto).
    """
    suppressed = [d for d in dirs_before if d not in kept]
    print_dual("", f_out)
    print_dual("-" * 60, f_out)
    print_dual(color_text(f"DEFORMATION DIRECTIONS (symmetry-method {symmetry_method})", 'cyan').center(60), f_out)
    print_dual("-" * 60, f_out)
    print_dual(f"  Detected symmetry : point group {point_group} -- {symmetry.operations_summary(ops)}", f_out)
    print_dual("-" * 60, f_out)
    print_dual(f"  {'Direction':<12}{'Status':<14}Reconstructed from", f_out)
    print_dual(f"  {'-' * 56}", f_out)
    for d in dirs_before:
        if d in kept:
            status = color_text("RUN".ljust(14), 'green')
            print_dual(f"  {d:<12}{status}--", f_out)
        else:
            source = mapping.get(d, "symmetry-allowed fit") if mapping else "symmetry-allowed fit"
            status = color_text("SUPPRESSED".ljust(14), 'yellow')
            print_dual(f"  {d:<12}{status}{source}", f_out)
    print_dual(f"  {'-' * 56}", f_out)
    if suppressed:
        print_dual(f"  {len(kept)} of {len(dirs_before)} run -- {len(suppressed)} suppressed by "
                   "symmetry (exactly reconstructible by stb-elasticAnalysis)", f_out)
    else:
        print_dual(f"  All {len(dirs_before)} direction(s) run -- no symmetry reduction available "
                   "for this point group.", f_out)
    print_dual("-" * 60, f_out)


# ==========================================
#                  MAIN
# ==========================================

def main():
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

    parser = argparse.ArgumentParser(
        description="Stage 1 of 2: applies one or more canonical Voigt strains to a SIESTA "
                    "structure and writes one ready-to-run folder (structure + calc.fdf + "
                    "config_extra.fdf + pseudopotentials) per scanned strain value, forcing a "
                    "pure single-point SCF at each imposed strain. Doesn't run SIESTA -- run "
                    "each folder yourself, then use stb-elasticAnalysis.",
        epilog="Example usage:\n"
               "  stb-elasticInputs -s structure.fdf -c calc.fdf -p dojo --dirs all\n"
               "  stb-elasticInputs -s structure.fdf -c calc.fdf --dirs xx yy xy --max 1.5 --steps 6 --save-report\n",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument("-s", "--structure", required=True,
                        help="Input structure.fdf -- lattice + atomic coordinates. Only this "
                             "file's %%block LatticeVectors is rewritten per strain step; "
                             "everything else is preserved verbatim. Atomic coordinates MUST be "
                             "fractional (checked in [1] of the report -- Cartesian positions "
                             "would not scale with the strained cell).")
    parser.add_argument("-c", "--calc", required=True,
                        help="calc.fdf template -- SCF/basis/MD settings, normally %%including "
                             "the structure file above by name. Copied into every generated "
                             "folder with one line added: '%%include config_extra.fdf', "
                             "prepended at the very top (see [4] of the report for why) -- this "
                             "forces every generated folder to a pure single-point SCF, "
                             "regardless of what MD settings this template itself contains.")
    parser.add_argument("-p", "--pseudo-dir", default="",
                        help="Pseudopotentials source (optional): a bundled bank "
                             f"({', '.join(BANKS)}) or a folder path. If given, the required "
                             "pseudopotential for every species in the structure is copied into "
                             "every generated folder, so it's immediately ready for SIESTA.")
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
    parser.add_argument("--max", type=float, default=2.0, help="Max strain %% (default: 2.0)")
    parser.add_argument("--steps", type=int, default=4, help="Steps per direction (default: 4)")
    parser.add_argument("-o", "--output-dir", default="elastic_runs",
                        help="Top-level directory (default: elastic_runs). Run folders go into "
                             "its own '<output-dir>/<direction>/strain_<direction>_<val>/' "
                             "subfolder, same convention as stb-strain's own --output-dir -- so "
                             "stb-elasticAnalysis can just 'cd <output-dir>' and find everything "
                             "(including reference_structure.fdf) with no extra flags.")
    parser.add_argument("--vacuum-gap", type=float, default=10.0,
                        help="Vacuum gap threshold in Ang used to detect which lattice axes are "
                             "periodic vs. vacuum-padded (default: 10.0), same convention as "
                             "stb-kgrid/stb-strain. Requested directions touching a vacuum-padded "
                             "axis are skipped with a warning instead of generated.")
    parser.add_argument("--symprec", type=float, default=0.01,
                        help="Symmetry-detection tolerance (default: 0.01, pymatgen's own default), "
                             "used for the symmetry-based direction reduction and the "
                             "informational note about triclinic/monoclinic structures. Loosen "
                             "further (e.g. 0.02-0.05) for a structure relaxed with a looser "
                             "force tolerance -- a tighter value than the relaxation's own "
                             "residual numerical noise can misdetect the true space group.")
    parser.add_argument("--angle-tolerance", type=float, default=5.0,
                        help="Symmetry angle tolerance in degrees (default: 5.0, pymatgen's own "
                             "default), same use as --symprec.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to <output-dir>/{REPORT_FILE}. Off by default.")
    parser.add_argument("-v", "--version", action="version",
                        version=f"stb-elasticInputs {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.method == "energy" and args.dirs != ["all"]:
        sys.exit(color_text(
            "[ERROR] --method energy always auto-selects its own strain patterns (pure and "
            "combined) from the structure's point group -- an explicit --dirs list isn't "
            "supported yet. Omit --dirs (or pass --dirs all, the default) to use it.", 'red'))

    if not os.path.isfile(args.structure):
        sys.exit(color_text(f"[ERROR] Structure file '{args.structure}' not found.", 'red'))
    if not os.path.isfile(args.calc):
        sys.exit(color_text(f"[ERROR] Calc file '{args.calc}' not found.", 'red'))
    if args.pseudo_dir:
        try:
            args.pseudo_dir = resolve_pseudo_source(args.pseudo_dir)
        except ValueError as e:
            sys.exit(color_text(f"[ERROR] {e}", 'red'))

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - Elastic Generator",
            "Strain Structure Generator (Batch Mode)",
            f"Version {VERSION} | University of Brasilia - 2026",
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
        "===== STB-ELASTICINPUTS STAGE 1 REPORT (ELASTIC CONSTANTS PREP) =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Structure file    : {args.structure}", f_out)
    print_dual(f"Calc template     : {args.calc}", f_out)
    print_dual(f"Pseudo source     : {args.pseudo_dir or '(not given)'}", f_out)
    print_dual(f"Method            : {args.method}", f_out)
    print_dual(f"Directions        : {' '.join(args.dirs)}", f_out)
    if args.method == "stress":
        print_dual(f"Symmetry method   : {args.symmetry_method}", f_out)
    print_dual(f"Strain range      : +/-{args.max}%, {args.steps} steps", f_out)
    print_dual(f"Output dir        : {args.output_dir}", f_out)
    print_dual(f"Vacuum-gap        : {args.vacuum_gap} Ang", f_out)
    print_dual(f"Symmetry tol      : symprec={args.symprec}, angle-tolerance={args.angle_tolerance} deg", f_out)
    print_dual(f"Save report       : {'yes' if args.save_report else 'no'}", f_out)
    if report_path:
        print_dual(f"Report file       : {report_path}", f_out)

    print_section("[1] INPUT STRUCTURE", f_out)
    try:
        structure = structure_io.read_fdf(args.structure)
    except (ValueError, FileNotFoundError) as e:
        fail(str(e))

    counts = structure_io.atom_counts(structure)
    species = sorted(sp for sp, n in counts.items() if n > 0)
    composition = ", ".join(f"{sp}{n}" for sp, n in counts.items())
    cell_volume = abs(np.linalg.det(structure.lattice))
    print_table(["Quantity", "Value"], [
        (["Composition", composition], None),
        (["Total atoms", str(sum(counts.values()))], None),
        (["Coordinate format", structure.coord_format], None),
        (["Lattice constant", f"{structure.lattice_constant} Ang"], None),
        (["Cell volume (unstrained)", f"{cell_volume:.4f} Ang^3"], None),
    ], f_out)

    if structure.coord_format == 'cartesian':
        fail("Atomic coordinates in this file are Cartesian, not fractional. stb-elasticInputs "
             "rewrites only %block LatticeVectors and leaves atomic positions untouched -- "
             "fractional coordinates scale automatically with the new strained cell "
             "(physically correct), but Cartesian coordinates would NOT, silently placing "
             "atoms at the wrong positions under the strained lattice. Convert to fractional "
             "coordinates first and re-run.")

    lattice = structure_io.raw_lattice_vectors(structure)

    # Save an UNDEFORMED copy for stb-elasticAnalysis's later symmetry
    # detection. Deliberately NOT read from inside a strain_*/ folder for
    # that: any single-axis normal strain (e.g. 'xx' alone) genuinely lowers
    # the crystal's symmetry by construction (a cubic cell strained along x
    # only has a != b == c, exactly tetragonal, not "cubic with noise") --
    # verified live catching this exact mistake (detected 4/mmm instead of
    # m-3m for a strained folder). This plain, un-strained copy is the only
    # correct reference. Lives at the TOP of --output-dir (not per
    # -direction), matching where stb-elasticAnalysis's own
    # --reference-structure default expects to find it.
    try:
        shutil.copy(args.structure, os.path.join(args.output_dir, REFERENCE_FILE))
    except OSError:
        pass  # best-effort; stb-elasticAnalysis just skips symmetry detection without it

    print_section("[2] DIMENSIONALITY & DIRECTION VALIDATION", f_out)
    # Detect which axes are periodic vs. vacuum-padded, same convention as
    # stb-kgrid/stb-kpath/stb-strain -- straining a vacuum gap just moves
    # empty space, not the material.
    positions = np.array([pos for _, pos in structure.atoms])
    is_cartesian = structure.coord_format == 'cartesian'
    frac_coords = kspace.to_fractional(positions, structure.lattice, is_cartesian)
    vacuum_axes = kspace.detect_vacuum_axes(frac_coords, structure.lattice, args.vacuum_gap)
    print_dual(f"Detected dimensionality: {kspace.dimensionality_label(vacuum_axes)}", f_out)
    print_table(["Axis", "Status"], [
        ([letter, "vacuum" if is_vac else "periodic"], 'yellow' if is_vac else None)
        for letter, is_vac in zip('xyz', vacuum_axes)
    ], f_out)

    axis_letters = ['x', 'y', 'z']
    periodic_axes = {i for i in range(3) if not vacuum_axes[i]}

    # Expand 'all' keyword. Only 'all' is eligible for the symmetry-based
    # reduction below -- an explicit list is respected exactly as given
    # (the way to deliberately run redundant directions and get
    # stb-elasticAnalysis's diagnostic/pooling across them).
    used_all = 'all' in args.dirs
    dirs_to_run = ['xx', 'yy', 'zz', 'xy', 'xz', 'yz'] if used_all else list(args.dirs)

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
                print_dual(color_text(
                    f"[WARNING] Skipping direction '{d}': touches vacuum-padded axis/axes "
                    f"{letters} (detected via --vacuum-gap {args.vacuum_gap} Ang); straining a "
                    "vacuum gap doesn't correspond to any physical deformation.", 'yellow'), f_out)
                continue
            filtered_dirs.append(d)
        dirs_to_run = filtered_dirs

        if not dirs_to_run:
            periodic_str = ', '.join(axis_letters[i] for i in sorted(periodic_axes)) or \
                'none (structure is fully isolated / 0D)'
            fail(f"None of the requested directions are physically valid for this structure -- "
                 f"periodic axis/axes available: {periodic_str}.")
    else:
        # --method energy: dirs_to_run (from the 'all' expansion above,
        # --dirs is forced to 'all' by the early check) is discarded here in
        # favor of the full 21-pattern pool (6 pure + 15 combined), filtered
        # against the SAME vacuum-padded axes.
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
            fail(f"No strain pattern is physically valid for this structure -- periodic "
                 f"axis/axes available: {periodic_str}.")

    print_section("[3] DEFORMATION DIRECTIONS / SYMMETRY REDUCTION", f_out)
    if args.method == "stress":
        # Symmetry-based reduction: only for the 'all' convenience default (an
        # explicit --dirs list is never second-guessed). Keeps 1 representative
        # per group of symmetry-equivalent directions -- stb-elasticAnalysis
        # reconstructs the skipped ones' stiffness columns exactly via the same
        # point-group operation (see core/symmetry.py::find_mapping_operation),
        # so this is a real DFT-cost cut, not an approximation.
        if used_all and args.symmetry_method == "basic":
            try:
                pmg_structure = structure_io.to_pymatgen(structure)
                groups, point_group, ops = symmetry.equivalent_strain_modes(
                    pmg_structure, args.symprec, args.angle_tolerance)
                dirs_before = list(dirs_to_run)
                kept, seen, mapping = [], set(), {}
                for d in dirs_to_run:
                    if d in seen:
                        continue
                    group = next(g for g in groups if d in g)
                    survivors = [m for m in dirs_to_run if m in group]
                    representative = survivors[0]
                    kept.append(representative)
                    for skipped_d in survivors[1:]:
                        mapping[skipped_d] = representative
                    seen.update(group)
                dirs_to_run = kept
                print_symmetry_table(dirs_before, kept, point_group, ops, "basic", mapping=mapping, f_out=f_out)
            except Exception:
                pass  # optimization only, never blocks -- falls back to running everything
        elif used_all and args.symmetry_method == "full":
            try:
                pmg_structure = structure_io.to_pymatgen(structure)
                basis, point_group = symmetry.symmetry_allowed_basis(
                    pmg_structure, args.symprec, args.angle_tolerance)
                _, ops = symmetry.get_point_group_operations(
                    pmg_structure, args.symprec, args.angle_tolerance)
                dirs_before = list(dirs_to_run)
                kept = symmetry.select_directions_by_rank(basis, order=dirs_to_run)
                dirs_to_run = kept
                print_dual(f"{color_text('[INFO]', 'cyan')} Point group {point_group} has "
                           f"{len(basis)} independent elastic constant(s); {len(kept)} direction(s) "
                           "fully determine them (--symmetry-method full).", f_out)
                print_symmetry_table(dirs_before, kept, point_group, ops, "full", f_out=f_out)
            except Exception:
                pass  # optimization only, never blocks -- falls back to running everything
        else:
            print_dual(f"Explicit --dirs given ({', '.join(dirs_to_run)}) -- no symmetry "
                       "reduction applied.", f_out)

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
                print_dual(color_text(
                    f"[INFO] Detected {system} symmetry (point group {point_group}) -- this crystal "
                    "class always has some nonzero normal-shear coupling constants (C14-C16, C24-C26, "
                    "C34-C36, C45, C46, C56). stb-elasticAnalysis recovers them correctly, but only "
                    "from directions that were actually run: make sure to keep all 6 default "
                    "directions (--dirs all, the default) rather than a partial --dirs list, or those "
                    "coupling terms will read as 0 from missing data instead of their real value.",
                    'cyan'), f_out)
        except Exception:
            pass  # symmetry detection is advisory only, never blocks the tool
    else:
        try:
            pmg_structure = structure_io.to_pymatgen(structure)
            basis, point_group = symmetry.symmetry_allowed_basis(
                pmg_structure, args.symprec, args.angle_tolerance)
            _, ops = symmetry.get_point_group_operations(
                pmg_structure, args.symprec, args.angle_tolerance)
            kept = symmetry.select_energy_patterns_by_rank(basis, filtered_candidates)
            dirs_to_run = [name for name, _ in kept]
            achieved = len(kept)
            print_dual(f"{color_text('[INFO]', 'cyan')} Point group {point_group} has "
                       f"{len(basis)} independent elastic constant(s); {achieved} strain pattern(s) "
                       f"needed to determine "
                       f"{'all of them' if achieved >= len(basis) else 'the ones reachable without leaving the periodic plane/axis'} "
                       f"(of {len(filtered_candidates)} physically valid candidate(s), {len(candidates)} in total).", f_out)
            print_symmetry_table([name for name, _ in filtered_candidates], dirs_to_run,
                                 point_group, ops, "energy", f_out=f_out)
        except Exception as e:
            fail(f"--method energy requires symmetry detection to pick its strain patterns "
                 f"and it failed: {e}")

    print_dual(f"\nModes: {dirs_to_run}", f_out)

    print_section("[4] SINGLE-POINT SCF ENFORCEMENT", f_out)
    with open(args.calc) as f:
        original_calc_text = f.read()
    structure_basename = os.path.basename(args.structure)
    calc_basename = os.path.basename(args.calc)
    if "%include" not in original_calc_text or structure_basename not in original_calc_text:
        print_dual(color_text(
            f"[NOTE] Could not confirm '{calc_basename}' references '{structure_basename}' via "
            f"%include -- if your calc.fdf doesn't already include the structure file by this "
            f"exact name, add '%include {structure_basename}' to it (or the equivalent for your "
            "own convention) so SIESTA picks up the strained geometry in each generated folder.",
            'yellow'), f_out)

    before = read_md_state(original_calc_text)
    steps_label = f"{before['steps_key']}={before['steps']}" if before['steps_key'] else "(absent)"
    print_dual("Calc template (current state, before forcing):", f_out)
    print_dual(f"  MD.TypeOfRun={before['typeofrun'] or '(absent)'}  Steps: {steps_label}  "
               f"MD.VariableCell={before['variablecell'] or '(absent)'}", f_out)
    print_dual(
        "Every generated folder is forced to a pure single-point SCF -- NO ionic relaxation, "
        "NO cell relaxation -- regardless of the state above, since the imposed strain must "
        "stay exactly at the sampled geometry for the stress/energy fit to be meaningful. "
        f"Forced via '%include {EXTRA_FDF_FILE}' PREPENDED at the very top of every generated "
        f"{calc_basename} (before your own directives, including the structure %include) -- "
        "SIESTA's fdf reader is first-occurrence-wins for duplicate labels, so this ordering "
        "guarantees the forced values win even if your own template already sets any of them "
        "(the common case for a real relaxation calc.fdf).", f_out)

    extra_fdf_text = (
        "# Auto-generated by stb-elasticInputs -- forces a pure single-point SCF (no ionic or\n"
        "# cell relaxation) at the imposed-strain geometry, regardless of --calc's own settings.\n"
        "MD.TypeOfRun       CG\n"
        "MD.Steps           0\n"
        "MD.VariableCell    false\n"
    )
    print_dual(f"\n{EXTRA_FDF_FILE} (written into every generated folder):", f_out)
    for line in extra_fdf_text.rstrip("\n").split("\n"):
        print_dual(f"  {line}", f_out)
    forced_calc_text = prepend_include(original_calc_text, EXTRA_FDF_FILE)

    print_section("[5] PSEUDOPOTENTIALS", f_out)
    if args.pseudo_dir:
        found, missing = get_required_pseudos(species, args.pseudo_dir)
        print_dual(f"Source            : {args.pseudo_dir}", f_out)
        print_table(["Species", "Status"], [
            ([sp, "MISSING" if sp in missing else "found"], 'yellow' if sp in missing else None)
            for sp in species
        ], f_out)
        if missing:
            print_dual(color_text(
                f"[WARNING] Missing pseudopotential(s) for: {', '.join(missing)} -- these will "
                "need to be added manually to every generated folder.", 'yellow'), f_out)
        else:
            print_dual(color_text(
                "[OK] All required pseudopotentials found -- will be copied into every "
                "generated folder.", 'green'), f_out)
    else:
        print_dual("Not given (pass -p/--pseudo-dir -- a bundled bank or a folder path -- to "
                   "copy the required pseudopotential for every species into every generated "
                   "folder). Pseudopotentials will need to be added manually.", f_out)

    print_section("[6] GENERATED STRAIN FOLDERS", f_out)
    strains = np.linspace(-args.max, args.max, args.steps)

    folder_rows = []
    per_direction_count = {}
    for d in dirs_to_run:
        for s in strains:
            if abs(s) < 1e-9:
                s = 0.0

            delta = s / 100.0
            def_matrix = get_strain_matrix(d, delta)
            new_lattice = np.dot(lattice, def_matrix)

            folder = os.path.join(
                args.output_dir, d, f"strain_{d}_{'m' if s < 0 else ''}{abs(s):.2f}")
            os.makedirs(folder, exist_ok=True)

            output_structure = os.path.join(folder, structure_basename)
            structure_io.rewrite_fdf_lattice(args.structure, new_lattice, output_structure)

            with open(os.path.join(folder, calc_basename), "w") as fh:
                fh.write(forced_calc_text)
            with open(os.path.join(folder, EXTRA_FDF_FILE), "w") as fh:
                fh.write(extra_fdf_text)

            for sym in species:
                copy_pseudo(args.pseudo_dir, sym, folder)

            volume = abs(np.linalg.det(new_lattice)) * structure.lattice_constant ** 3
            folder_rows.append(([folder, d, f"{s:+.2f}",
                                f"{structure_basename}, {calc_basename}, {EXTRA_FDF_FILE}",
                                f"{volume:.4f}"], None))
            per_direction_count[d] = per_direction_count.get(d, 0) + 1

    print_table(["Folder", "Direction", "Strain (%)", "Files", "Volume (Ang^3)"], folder_rows, f_out)

    print_section("[7] SUMMARY & NEXT STEPS", f_out)
    print_dual(f"{len(folder_rows)} folder(s) written under '{args.output_dir}' "
               f"({len(per_direction_count)} direction(s)):", f_out)
    for d, n in per_direction_count.items():
        print_dual(f"  {d}: {n} folder(s)", f_out)
    if report_path:
        print_dual(f"Report            : {report_path}", f_out)
    print_dual(color_text("\nNext steps:", 'yellow'), f_out)
    print_dual(f"  1. Run SIESTA (single-point) in every "
               f"'{args.output_dir}/<direction>/strain_*/' folder.", f_out)
    print_dual(f"  2. Once they're done: cd {args.output_dir} && stb-elasticAnalysis", f_out)

    if f_out:
        f_out.close()

if __name__ == "__main__":
    main()
