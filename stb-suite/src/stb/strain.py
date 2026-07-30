#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""Stage 1 of 2: applies a uniaxial/biaxial Cartesian strain to a SIESTA
structure and writes one '<output-dir>/<direction>/strain_<direction>_<pct>/'
folder per scanned strain value, ready to run -- structure + calc.fdf +
config_extra.fdf + linked pseudopotentials -- for a stress-strain elastic
-constant workflow. Every direction gets its own <output-dir>/<direction>/
subfolder (not flat directly in <output-dir>) so 2 different --stdir runs
into the same --output-dir never mix folders or clobber each other's
--save-report file, and so stb-strainAnalysis can point --dir at one
direction's subfolder alone for a single-direction report, or at
<output-dir> itself with --compare to read every direction found under it
side by side. Stage 2 is stb-strainAnalysis, which aggregates the finished
SIESTA runs back into a stress-strain fit. Doesn't run SIESTA itself -- run
each folder yourself, then use stb-strainAnalysis.

Physics reviewed: the deformation itself (apply_cartesian_strain) uses the
standard small-strain deformation-gradient convention, F = I + strain
-tensor, new_vec = F . lattice_vectors.T -- correct as written, unchanged
by this rewrite.

Takes TWO separate input files -- --structure (lattice + atomic
coordinates) and --calc (SCF/basis/pseudopotential/MD settings, normally
%including the structure file by name) -- same split as her.py/gqca.py's
own Stage 1 tools: --calc is a reusable template copied into every
generated folder (with one %include line added, see --relax-mode below),
while --structure is what actually gets strained and rewritten per
folder. --pseudo-dir (optional, a bundled bank or a folder path, same
convention as stb-inputfile) links the required pseudopotential for every
species into each generated folder too, so every folder is immediately
ready to hand to SIESTA.

--relax-mode selects HOW each generated folder's cell is allowed to
respond to the imposed strain, via SIESTA's own %block
Geometry.Constraints ('stress N' lines, N=1..6 in standard Voigt order
XX/YY/ZZ/YZ/XZ/XY) -- written into a new config_extra.fdf file, included
near the top of calc.fdf right after its %include of the structure file
(the exact position/syntax verified against a real user relaxation
calc.fdf and a real stress-scan calc.fdf, both of which already carry this
same block commented out as a template). Both modes share ONE mechanism
(compute_voigt_status/build_geometry_constraints_block below) rather than
two -- 'cell-fixed' is simply the special case that fixes all 6
components:
 - 'cell-fixed' (relaxed-ion method): all 6 Voigt stress components are
   fixed -- the cell stays exactly at the imposed strain (physically
   equivalent to disabling MD.VariableCell entirely) while atomic
   positions still relax, via whatever MD.TypeOfRun/MD.Steps is already in
   --calc (not touched here -- presumed to be the same calc.fdf used for
   the structure's original relaxation).
 - 'stress-constrained': only the imposed strain direction's own
   diagonal component(s) are fixed, plus anything (diagonal or shear)
   touching a vacuum axis (straining/relaxing a vacuum gap is physically
   meaningless, same reasoning as the vacuum-direction rejection in [2]
   below); every other periodic direction is left free to relax to zero
   stress.
This tool does NOT force MD.VariableCell itself -- Geometry.Constraints
only has an effect when the cell is actually variable, so that precondition
is the user's own responsibility in --calc; [4] below only checks/warns
about it, matching the same reasoning for MD.Steps (a relaxation with zero
steps is a no-op in either mode -- also only warned about, never forced).

Two real gaps fixed here (relative to the original single-file version):
 - --structure must be fractional-coordinate (checked, not just
   documented): rewrite_fdf_lattice only ever replaces %block
   LatticeVectors, leaving atomic positions untouched. That's physically
   correct ONLY for fractional coordinates (they scale for free with the
   new cell) -- Cartesian coordinates would silently end up wrong under
   the strained lattice, since their raw numbers wouldn't follow the
   deformation at all. Previously only a docstring warning; now a hard,
   checked error.
 - --calc's MD.Steps/MD.NumCGsteps count is read (for the [4] warning
   above) recognizing BOTH spellings -- a real relaxation calc.fdf and
   stb-inputfile's own generated template both use 'MD.Steps', not
   'MD.NumCGsteps' (the only one this tool originally recognized).

Output/report style (v2.0.0 -> v3.0.0, the --relax-mode rewrite): a
numbered [0]...[7] report and an opt-in --save-report (off by default),
same shape as her.py/gqca.py's own Stage 1 reports.
"""

VERSION = "3.0.0"

import os
import re
import sys
import argparse
import numpy as np
from stb.core import structure_io, kspace, symmetry
from stb.core.cli import color_text, show_intro, print_dual, print_section, print_table
from stb.core.pseudopotentials import BANKS, resolve_pseudo_source, get_required_pseudos, link_pseudo

REPORT_FILE = "strain_stage1.txt"
EXTRA_FDF_FILE = "config_extra.fdf"

_MD_TYPEOFRUN_VALUE_RE = re.compile(r'MD\.TypeOfRun\s+(\S+)', re.IGNORECASE)
_MD_STEPS_VALUE_RE = re.compile(r'MD\.Steps\s+(\d+)', re.IGNORECASE)
_MD_NUMCGSTEPS_VALUE_RE = re.compile(r'MD\.NumCGsteps\s+(\d+)', re.IGNORECASE)
_MD_VARIABLECELL_VALUE_RE = re.compile(r'MD\.VariableCell\s+(\S+)', re.IGNORECASE)
_SIESTA_TRUE_VALUES = {'t', 'true', '.true.', 'yes'}

_VOIGT_LABELS = {1: 'XX', 2: 'YY', 3: 'ZZ', 4: 'YZ', 5: 'XZ', 6: 'XY'}
_VOIGT_AXES = {1: (0,), 2: (1,), 3: (2,), 4: (1, 2), 5: (0, 2), 6: (0, 1)}  # 0=x, 1=y, 2=z


def determine_strain_type(direction):
    """Determine strain type based on direction input."""
    if len(direction) == 1:  # x, y, z
        return 'uniaxial'
    elif len(direction) == 2:  # xy, xz, yx, yz, zx, zy
        return 'biaxial'
    else:
        raise ValueError("Invalid direction. Use x, y, z for uniaxial or combinations like xy, yz for biaxial.")

def normalize_direction(direction):
    """Normalize direction input (e.g., yx -> xy)."""
    if len(direction) == 2:
        return ''.join(sorted(direction.lower()))
    return direction.lower()


def _strip_fdf_comments(calc_text):
    """Removes everything from the first '#' onward on every line (SIESTA's
    own fdf comment convention -- also used for trailing inline comments
    after a real value, e.g. 'DM.UseSaveDM .true.  #(...)'). Used before
    scanning calc_text for a directive's CURRENT value -- otherwise a
    plain-English comment merely mentioning a directive's name (e.g. '##
    MD.VariableCell should stay true', a real, caught-live bug while
    writing this tool's own example) could be misread as the actual value
    if that comment happens to appear earlier in the file than the real
    directive line."""
    return "\n".join(line.split('#', 1)[0] for line in calc_text.splitlines())


def _read_effective_steps(calc_text):
    """Returns (key_name, value) for the relaxation step count -- prefers
    'MD.Steps' (what a real relaxation calc.fdf and stb-inputfile's own
    generated template both use), falling back to 'MD.NumCGsteps' (an
    older/alternate spelling some templates may still use) if 'MD.Steps'
    is absent. (None, None) if neither is present."""
    calc_text = _strip_fdf_comments(calc_text)
    steps_match = _MD_STEPS_VALUE_RE.search(calc_text)
    if steps_match:
        return 'MD.Steps', steps_match.group(1)
    numcg_match = _MD_NUMCGSTEPS_VALUE_RE.search(calc_text)
    if numcg_match:
        return 'MD.NumCGsteps', numcg_match.group(1)
    return None, None


def _read_md_state(calc_text):
    """Read-only scan of calc_text's CURRENT MD.TypeOfRun/MD.Steps (or
    MD.NumCGsteps)/MD.VariableCell directives, for reporting state in [4]
    below -- this tool never rewrites any of these, only the %include line
    for config_extra.fdf, so this is purely informational/advisory. A
    value is None when the tag is absent (a normal template, not a
    malformed one)."""
    stripped = _strip_fdf_comments(calc_text)
    run_match = _MD_TYPEOFRUN_VALUE_RE.search(stripped)
    varcell_match = _MD_VARIABLECELL_VALUE_RE.search(stripped)
    steps_key, steps_value = _read_effective_steps(calc_text)
    return {
        'typeofrun': run_match.group(1) if run_match else None,
        'steps_key': steps_key,
        'steps': steps_value,
        'variablecell': varcell_match.group(1) if varcell_match else None,
    }


def _is_siesta_true(value):
    """True if value.lower() is one of SIESTA's boolean-true spellings."""
    return value is not None and value.lower() in _SIESTA_TRUE_VALUES


def compute_voigt_status(norm_dir, vacuum_axes, relax_mode):
    """Returns {1..6: reason_str or None} -- a reason means 'fix this Voigt
    stress component in %block Geometry.Constraints'; None means 'leave it
    free to relax'. One shared rule for both --relax-mode choices:
    - 'cell-fixed': all 6 components fixed -- the cell stays exactly at
      the imposed strain (physically equivalent to disabling
      MD.VariableCell entirely, just expressed as constraints instead).
    - 'stress-constrained': fixes only the diagonal component(s) (1/2/3)
      whose axis letter is in norm_dir (the imposed strain itself) OR
      whose axis is vacuum-padded (protects the vacuum gap from spurious
      relaxation); plus any shear component (4/5/6) whose 2 axes include a
      vacuum axis (tilting the cell into the vacuum direction is just as
      physically meaningless as straining along it -- same reasoning as
      the vacuum-direction rejection in [2]). Every other periodic
      direction is left free to relax to zero stress."""
    status = {}
    if relax_mode == 'cell-fixed':
        for voigt_idx in _VOIGT_LABELS:
            status[voigt_idx] = "cell-fixed mode: entire cell locked at the imposed strain"
        return status
    for voigt_idx, axes in _VOIGT_AXES.items():
        if len(axes) == 1:
            axis_letter = 'xyz'[axes[0]]
            if axis_letter in norm_dir:
                status[voigt_idx] = f"imposed strain direction ('{norm_dir}')"
            elif vacuum_axes[axes[0]]:
                status[voigt_idx] = f"vacuum-axis protection ('{axis_letter}' is vacuum-padded)"
            else:
                status[voigt_idx] = None
        else:
            if any(vacuum_axes[a] for a in axes):
                vac_letters = ', '.join('xyz'[a] for a in axes if vacuum_axes[a])
                status[voigt_idx] = f"vacuum-axis protection (involves '{vac_letters}')"
            else:
                status[voigt_idx] = None
    return status


def build_geometry_constraints_block(voigt_status):
    """%block Geometry.Constraints ... %endblock text, one 'stress N  #
    Fixes <LABEL>' line per component with a reason (i.e. status[N] is not
    None) -- same label format as the 2 real reference calc.fdf files this
    was verified against, only WITHOUT the '#' comment markers (active,
    not a menu of options) and only listing components actually fixed (the
    free ones aren't shown at all, so the file never suggests an inactive
    option is available)."""
    lines = ["%block Geometry.Constraints"]
    for voigt_idx in sorted(voigt_status):
        if voigt_status[voigt_idx] is not None:
            lines.append(f"  stress {voigt_idx}  # Fixes {_VOIGT_LABELS[voigt_idx]}")
    lines.append("%endblock Geometry.Constraints")
    return "\n".join(lines) + "\n"


def insert_include_after_structure(calc_text, structure_basename, include_name=EXTRA_FDF_FILE):
    """Inserts '%include <include_name>' right after the existing
    '%include <structure_basename>' line -- the exact position the
    Geometry.Constraints block already occupies (commented out, as a
    template) in the 2 real reference calc.fdf files this was verified
    against. Falls back to prepending at the very top of the file if no
    such line is found (calc.fdf not %including the structure file by
    this exact name already triggers a separate [NOTE] in [4], so this
    never fails silently). --calc itself is otherwise untouched -- no
    other directive (MD.VariableCell/MD.TypeOfRun/MD.Steps) is rewritten
    here."""
    include_line = f"%include {include_name}"
    pattern = re.compile(
        r'^([ \t]*%include[ \t]+' + re.escape(structure_basename) + r'[ \t]*)$',
        re.IGNORECASE | re.MULTILINE)
    new_text, count = pattern.subn(r'\1\n' + include_line, calc_text, count=1)
    if count == 0:
        new_text = include_line + "\n\n" + calc_text
    return new_text


def print_axis_symmetry_table(requested_axis, groups, point_group, ops, vacuum_axes, f_out=None):
    """Pretty table of the 3 Cartesian axes' symmetry-equivalence groups,
    for a uniaxial strain request -- same visual convention as
    stb-elasticInputs' deformation-direction table (see elastic_inputs.py::
    _print_symmetry_table), reusing core/symmetry.py::operations_summary for
    the reduced point-group/operations line. Always shown (even if nothing
    is equivalent), same "give full information regardless" philosophy.

    Public (not underscore-prefixed): stb_suite.py's interactive wrapper
    imports this directly to build its own pre-run "generate the equivalent
    axis/axes too?" prompt, since stb-strain's CLI only ever strains one
    direction per invocation and the wrapper needs the same table before
    deciding whether to loop over extra directions. `f_out`, if given,
    mirrors every line into the numbered [3] report section too (added
    once this became a section of stb-strain's own stage-1 report, not
    just a standalone console preview); the wrapper's own call omits it.

    Returns the list of axis letters equivalent to `requested_axis` (periodic
    ones only -- vacuum-padded axes are excluded even if symmetry-equivalent,
    since they can't physically be strained).
    """
    axis_letters = ['x', 'y', 'z']
    requested_idx = axis_letters.index(requested_axis)
    group = next(g for g in groups if requested_idx in g)

    print_dual("", f_out)
    print_dual("-" * 60, f_out)
    print_dual(color_text(f"AXIS SYMMETRY (uniaxial direction '{requested_axis}')", 'cyan').center(60), f_out)
    print_dual("-" * 60, f_out)
    print_dual(f"  Detected symmetry : point group {point_group} -- {symmetry.operations_summary(ops)}", f_out)
    print_dual("-" * 60, f_out)
    print_dual(f"  {'Axis':<8}{'Status':<14}Equivalent to", f_out)
    print_dual(f"  {'-' * 52}", f_out)
    equivalent = []
    for i, letter in enumerate(axis_letters):
        if i == requested_idx:
            status = color_text("REQUESTED".ljust(14), 'green')
            print_dual(f"  {letter:<8}{status}--", f_out)
        elif vacuum_axes[i]:
            status = color_text("VACUUM".ljust(14), 'yellow')
            print_dual(f"  {letter:<8}{status}(not periodic -- can't be strained)", f_out)
        elif i in group:
            status = color_text("EQUIVALENT".ljust(14), 'cyan')
            print_dual(f"  {letter:<8}{status}{requested_axis}", f_out)
            equivalent.append(letter)
        else:
            print_dual(f"  {letter:<8}{'INDEPENDENT':<14}--", f_out)
    print_dual(f"  {'-' * 52}", f_out)
    if equivalent:
        letters = ', '.join(equivalent)
        verb = "is" if len(equivalent) == 1 else "are"
        print_dual(f"  {letters} {verb} equivalent to '{requested_axis}' by symmetry -- straining "
                   f"{'it' if len(equivalent) == 1 else 'them'} should give the same mechanical "
                   "response; you may not need to compute both.", f_out)
    else:
        print_dual(f"  No other periodic axis is equivalent to '{requested_axis}' for this point group.", f_out)
    print_dual("-" * 60, f_out)
    return equivalent


_CANONICAL_DIRECTIONS = ('x', 'y', 'z', 'xy', 'xz', 'yz')


def print_direction_selection_table(groups, point_group, ops, vacuum_axes, f_out=None):
    """Classifies all 6 canonical stb-strain directions (x, y, z uniaxial; xy,
    xz, yz biaxial) into VACUUM (touches a vacuum-padded axis -- can't be
    strained), REDUNDANT (symmetry-equivalent to an already-listed
    representative), or INDEPENDENT (a representative worth actually running),
    and prints the classification -- same visual convention as
    print_axis_symmetry_table above, but covering the whole direction menu
    upfront (before any direction is chosen) instead of one requested axis at a
    time.

    Biaxial equivalence is derived from the SAME axis `groups`
    equivalent_cartesian_axes already computes for uniaxial directions, no
    separate biaxial symmetry computation needed: a biaxial direction's strain
    tensor is I - e_c(x)e_c, where c is its own EXCLUDED axis (e.g. 'xy' =
    diag(1,1,0) excludes 'z' -- see apply_cartesian_strain, this is a diagonal
    biaxial strain, not a shear). Under a point-group rotation R this
    transforms to I - (Rc)(x)(Rc), so two biaxial directions are
    symmetry-equivalent iff their excluded axes are related by some point
    -group operation -- exactly the criterion `groups` already encodes for
    single axes. Uniaxial and biaxial directions are tracked as separate
    equivalence spaces (never cross-equivalent to each other).

    Returns the ordered list (canonical order) of INDEPENDENT direction
    strings -- the representatives worth offering the user to choose from.
    """
    axis_letters = 'xyz'
    axis_group_of = {}
    for group in groups:
        frozen = frozenset(group)
        for i in group:
            axis_group_of[i] = frozen

    print_dual("", f_out)
    print_dual("-" * 60, f_out)
    print_dual(color_text("STRAIN DIRECTION SELECTION", 'cyan').center(60), f_out)
    print_dual("-" * 60, f_out)
    print_dual(f"  Detected symmetry : point group {point_group} -- {symmetry.operations_summary(ops)}", f_out)
    print_dual("-" * 60, f_out)
    print_dual(f"  {'Direction':<11}{'Type':<10}{'Status':<14}Notes", f_out)
    print_dual(f"  {'-' * 56}", f_out)

    independent = []
    representative_of = {'uniaxial': {}, 'biaxial': {}}
    for direction in _CANONICAL_DIRECTIONS:
        kind = 'uniaxial' if len(direction) == 1 else 'biaxial'
        involved = [axis_letters.index(c) for c in direction]

        if any(vacuum_axes[i] for i in involved):
            vac_letters = ', '.join(axis_letters[i] for i in involved if vacuum_axes[i])
            status = color_text("VACUUM".ljust(14), 'yellow')
            print_dual(f"  {direction:<11}{kind:<10}{status}involves vacuum axis '{vac_letters}'", f_out)
            continue

        if kind == 'uniaxial':
            group_key = axis_group_of[involved[0]]
        else:
            excluded_letter = next(c for c in axis_letters if c not in direction)
            group_key = axis_group_of[axis_letters.index(excluded_letter)]

        seen = representative_of[kind]
        if group_key in seen:
            rep = seen[group_key]
            status = color_text("REDUNDANT".ljust(14), 'cyan')
            print_dual(f"  {direction:<11}{kind:<10}{status}equivalent to '{rep}'", f_out)
        else:
            seen[group_key] = direction
            independent.append(direction)
            status = color_text("INDEPENDENT".ljust(14), 'green')
            print_dual(f"  {direction:<11}{kind:<10}{status}--", f_out)

    print_dual(f"  {'-' * 56}", f_out)
    print_dual(f"  {len(independent)} independent, non-vacuum direction(s): {', '.join(independent)}", f_out)
    print_dual("-" * 60, f_out)
    return independent


def apply_cartesian_strain(lattice_vectors, strain, direction):
    """
    Apply uniaxial or biaxial strain in Cartesian coordinates.

    Args:
        lattice_vectors: 3x3 numpy array of lattice vectors
        strain: strain value (positive for tension, negative for compression)
        direction: strain direction (x, y, z, xy, xz, yz, etc.)

    Returns:
        Strained lattice vectors
    """
    # Normalize and validate direction
    direction = normalize_direction(direction)
    valid_directions = {'x', 'y', 'z', 'xy', 'xz', 'yz'}
    if direction not in valid_directions:
        raise ValueError(f"Invalid direction '{direction}'. Use x, y, z, xy, xz, or yz.")

    # Create strain tensor
    strain_tensor = np.zeros((3, 3))

    if len(direction) == 1:  # Uniaxial
        if direction == 'x':
            strain_tensor[0, 0] = strain
        elif direction == 'y':
            strain_tensor[1, 1] = strain
        elif direction == 'z':
            strain_tensor[2, 2] = strain

    else:  # Biaxial
        if 'x' in direction:
            strain_tensor[0, 0] = strain
        if 'y' in direction:
            strain_tensor[1, 1] = strain
        if 'z' in direction:
            strain_tensor[2, 2] = strain

    # Apply strain transformation: new_vec = (I + ε) · vec
    identity = np.eye(3)
    transformation = identity + strain_tensor

    # Transform each lattice vector
    strained_vectors = np.dot(transformation, lattice_vectors.T).T

    return strained_vectors

def main():
    parser = argparse.ArgumentParser(
        description="Stage 1 of 2: applies strain in Cartesian coordinates to a SIESTA structure "
                   "and writes one ready-to-run folder (structure + calc.fdf + config_extra.fdf + "
                   "pseudopotentials) per scanned strain value. Type (uniaxial/biaxial) is "
                   "inferred from direction; --relax-mode selects how the cell is allowed to "
                   "respond (see its own help). Doesn't run SIESTA -- run each folder yourself, "
                   "then use stb-strainAnalysis.",
        epilog="Example usage:\n"
               "  stb-strain -s structure.fdf -c calc.fdf -p dojo --relax-mode cell-fixed --stdir x --stmin -2 --stmax 2 --step 1\n"
               "  stb-strain -s structure.fdf -c calc.fdf --relax-mode stress-constrained --stdir xy --stmax 5 --save-report\n",
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
                            "folder with one line added: '%%include config_extra.fdf' (the "
                            "Geometry.Constraints block --relax-mode builds, see [4] of the "
                            "report) -- no other directive in this file is modified.")
    parser.add_argument("--relax-mode", required=True, choices=["cell-fixed", "stress-constrained"],
                       help="How the cell responds to the imposed strain, via a %%block "
                            "Geometry.Constraints written to config_extra.fdf: "
                            "'cell-fixed' fixes all 6 Voigt stress components (cell stays "
                            "exactly at the imposed strain; atomic positions still relax, via "
                            "whatever MD.TypeOfRun/MD.Steps is already in --calc). "
                            "'stress-constrained' fixes only the imposed strain direction's own "
                            "component(s) (plus anything touching a vacuum axis); every other "
                            "periodic direction relaxes freely to zero stress. Either way, "
                            "MD.VariableCell must be enabled in your own --calc for "
                            "Geometry.Constraints to have any effect -- checked (not forced) in "
                            "[4] of the report.")
    parser.add_argument("-p", "--pseudo-dir", default="",
                       help="Pseudopotentials source (optional): a bundled bank "
                            f"({', '.join(BANKS)}) or a folder path. If given, the required "
                            "pseudopotential for every species in the structure is linked into "
                            "every generated folder, so it's immediately ready for SIESTA.")
    parser.add_argument("--stdir", required=True,
                       help="Direction of strain: x, y, z for uniaxial; xy, xz, yz, etc. for biaxial.")
    parser.add_argument("--stmin", type=float, default=0,
                       help="Minimum strain percentage (default: 0). Can be negative for compression.")
    parser.add_argument("--stmax", type=float, default=25,
                       help="Maximum strain percentage (default: 25).")
    parser.add_argument("--step", type=float, default=1,
                       help="Strain step percentage (default: 1).")
    parser.add_argument("-o", "--output-dir", default="strain_runs",
                       help="Top-level directory (default: strain_runs). The run folders go into "
                            "its own '<output-dir>/<direction>/strain_<direction>_<val>/' subfolder "
                            "-- so 2 different --stdir runs into the same --output-dir never mix "
                            "their folders, and stb-strainAnalysis can point --dir at one "
                            "direction's subfolder alone, or at <output-dir> itself with "
                            "--compare to analyze every direction found under it.")
    parser.add_argument("--vacuum-gap", type=float, default=10.0,
                       help="Vacuum gap threshold in Ang used to detect which lattice axes are "
                            "periodic vs. vacuum-padded (default: 10.0), same convention as "
                            "stb-kgrid/stb-kpath. Straining a vacuum-padded axis is physically "
                            "meaningless and is refused.")
    parser.add_argument("--symprec", type=float, default=1e-3,
                       help="Symmetry-detection tolerance (default: 1e-3, pymatgen's own default), "
                            "used only for the informational note about symmetry-equivalent "
                            "directions (uniaxial only) -- never blocks the run.")
    parser.add_argument("--angle-tolerance", type=float, default=5.0,
                       help="Symmetry angle tolerance in degrees (default: 5.0, pymatgen's own "
                            "default), same use as --symprec.")
    parser.add_argument("--save-report", action="store_true",
                       help=f"Also persist the report to <output-dir>/{REPORT_FILE}. Off by default.")
    parser.add_argument("-v", "--version", action="version",
                        version=f"stb-strain {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    args = parser.parse_args()

    if not os.path.isfile(args.structure):
        sys.exit(color_text(f"[ERROR] Structure file '{args.structure}' not found.", 'red'))
    if not os.path.isfile(args.calc):
        sys.exit(color_text(f"[ERROR] Calc file '{args.calc}' not found.", 'red'))
    if args.pseudo_dir:
        try:
            args.pseudo_dir = resolve_pseudo_source(args.pseudo_dir)
        except ValueError as e:
            sys.exit(color_text(f"[ERROR] {e}", 'red'))
    if args.stmin > args.stmax:
        sys.exit(color_text("[ERROR] Minimum strain cannot be greater than maximum strain.", 'red'))
    if args.step == 0:
        sys.exit(color_text("[ERROR] Step cannot be zero.", 'red'))
    try:
        strain_type = determine_strain_type(args.stdir)
    except ValueError as e:
        sys.exit(color_text(f"[ERROR] {e}", 'red'))
    norm_dir = normalize_direction(args.stdir)

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("STRAIN WORKFLOW -- STAGE 1: STRESS-STRAIN PREP", 'bold'))
    print("-" * 60)

    # Every strain folder for THIS direction lives under its own
    # <output-dir>/<direction>/ subfolder (not flat in <output-dir> itself) --
    # so stb-strainAnalysis can scan one direction's folders in isolation
    # (--dir <output-dir>/<direction>) or every direction at once
    # (--dir <output-dir> --compare) without them colliding, and so 2 separate
    # directions run into the same --output-dir never overwrite each other's
    # --save-report file.
    direction_dir = os.path.join(args.output_dir, norm_dir)
    os.makedirs(direction_dir, exist_ok=True)
    report_path = os.path.join(direction_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    def fail(message):
        print_dual(color_text(f"[FAIL] {message}", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_dual(color_text("===== STB-STRAIN STAGE 1 REPORT (STRESS-STRAIN PREP) =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Structure file    : {args.structure}", f_out)
    print_dual(f"Calc template     : {args.calc}", f_out)
    print_dual(f"Relax mode        : {args.relax_mode}", f_out)
    print_dual(f"Pseudo source     : {args.pseudo_dir or '(not given)'}", f_out)
    print_dual(f"Direction         : {norm_dir} ({strain_type})", f_out)
    print_dual(f"Strain range      : {args.stmin}% to {args.stmax}%, step {args.step}%", f_out)
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
        fail("Atomic coordinates in this file are Cartesian, not fractional. stb-strain "
             "rewrites only %block LatticeVectors and leaves atomic positions untouched -- "
             "fractional coordinates scale automatically with the new strained cell "
             "(physically correct), but Cartesian coordinates would NOT, silently placing "
             "atoms at the wrong positions under the strained lattice. Convert to fractional "
             "coordinates first and re-run.")

    lattice_vectors = structure_io.raw_lattice_vectors(structure)

    print_section("[2] DIMENSIONALITY & DIRECTION VALIDATION", f_out)
    # Detect which axes are periodic vs. vacuum-padded, same convention
    # as stb-kgrid/stb-kpath/stb-mlrelax -- straining a vacuum gap just
    # moves empty space, not the material.
    positions = np.array([pos for _, pos in structure.atoms])
    is_cartesian = structure.coord_format == 'cartesian'
    frac_coords = kspace.to_fractional(positions, structure.lattice, is_cartesian)
    vacuum_axes = kspace.detect_vacuum_axes(frac_coords, structure.lattice, args.vacuum_gap)
    print_dual(f"Detected dimensionality: {kspace.dimensionality_label(vacuum_axes)}", f_out)
    print_table(["Axis", "Status"], [
        ([letter, "vacuum" if is_vac else "periodic"], 'yellow' if is_vac else None)
        for letter, is_vac in zip('xyz', vacuum_axes)
    ], f_out)

    # Only check letters that are actually valid axes here -- an invalid
    # direction (e.g. 'w') is still caught later by apply_cartesian_strain's
    # own validation, with its original, more specific error message.
    axis_index = {'x': 0, 'y': 1, 'z': 2}
    vacuum_requested = [c for c in dict.fromkeys(norm_dir)
                         if c in axis_index and vacuum_axes[axis_index[c]]]
    if vacuum_requested:
        periodic = [c for c, is_vac in zip('xyz', vacuum_axes) if not is_vac]
        periodic_str = ', '.join(periodic) if periodic else 'none (structure is fully isolated / 0D)'
        fail(
            f"Direction '{norm_dir}' includes vacuum-padded axis/axes "
            f"'{', '.join(vacuum_requested)}' (detected via --vacuum-gap {args.vacuum_gap} Ang); "
            "straining a vacuum gap doesn't correspond to any physical deformation. "
            f"Periodic axis/axes available for this structure: {periodic_str}."
        )

    print_section("[3] AXIS SYMMETRY ADVISORY", f_out)
    # Advisory only (never blocks the run): if the requested uniaxial
    # direction is symmetry-equivalent to another periodic axis, straining
    # that axis instead/too would just repeat the same DFT calculation.
    # Biaxial directions are out of scope for this check -- see
    # core/symmetry.py::equivalent_cartesian_axes docstring.
    if strain_type == 'uniaxial' and norm_dir in ('x', 'y', 'z'):
        try:
            pmg_structure = structure_io.to_pymatgen(structure)
            groups, point_group = symmetry.equivalent_cartesian_axes(
                pmg_structure, args.symprec, args.angle_tolerance)
            _, ops = symmetry.get_point_group_operations(
                pmg_structure, args.symprec, args.angle_tolerance)
            print_axis_symmetry_table(norm_dir, groups, point_group, ops, vacuum_axes, f_out)
        except Exception as e:
            print_dual(color_text(
                f"[NOTE] Could not compute the symmetry-equivalence advisory ({e}) -- "
                "informational only, does not block the run.", 'yellow'), f_out)
    else:
        print_dual("Skipped -- symmetry-equivalence check is uniaxial-only (see "
                   f"core/symmetry.py::equivalent_cartesian_axes docstring). Requested "
                   f"direction is biaxial: {norm_dir}.", f_out)

    print_section("[4] RELAXATION MODE & CELL CONSTRAINTS", f_out)
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

    mode_rationale = {
        "cell-fixed": "the cell is locked exactly at the imposed strain (all 6 Voigt stress "
                      "components fixed); only atomic positions relax, via whatever "
                      "MD.TypeOfRun/MD.Steps is already in --calc.",
        "stress-constrained": "only the imposed strain direction's own stress component(s) are "
                      "fixed; every other periodic direction relaxes freely to zero stress, "
                      "while the vacuum axis (if any) stays protected.",
    }
    print_dual(f"Mode              : {args.relax_mode} -- {mode_rationale[args.relax_mode]}", f_out)

    before = _read_md_state(original_calc_text)
    steps_label = f"{before['steps_key']}={before['steps']}" if before['steps_key'] else "(absent)"
    print_dual(f"Calc template (current state, read-only -- nothing below is forced):", f_out)
    print_dual(f"  MD.TypeOfRun={before['typeofrun'] or '(absent)'}  Steps: {steps_label}  "
               f"MD.VariableCell={before['variablecell'] or '(absent)'}", f_out)
    if not _is_siesta_true(before['variablecell']):
        print_dual(color_text(
            "[WARNING] MD.VariableCell is not enabled in this calc.fdf (or is absent). The "
            "%block Geometry.Constraints written below only has an effect while the cell is "
            "actually variable -- add 'MD.VariableCell true' to your --calc yourself; this tool "
            "does not force it.", 'yellow'), f_out)
    if before['steps'] is None or before['steps'] == '0':
        print_dual(color_text(
            f"[WARNING] Relaxation step count is {steps_label} -- with 0 (or no) steps, "
            "atomic/cell relaxation is a no-op regardless of --relax-mode. Set MD.Steps (or "
            "MD.NumCGsteps) > 0 in your --calc for this mode to actually relax anything.",
            'yellow'), f_out)

    voigt_status = compute_voigt_status(norm_dir, vacuum_axes, args.relax_mode)
    print_table(["Component", "Label", "Status", "Reason"], [
        ([str(i), _VOIGT_LABELS[i], "FIXED" if voigt_status[i] else "free",
          voigt_status[i] or "--"], None if voigt_status[i] is None else 'cyan')
        for i in sorted(voigt_status)
    ], f_out)

    extra_fdf_text = (f"# Auto-generated by stb-strain (relax-mode={args.relax_mode})\n"
                       + build_geometry_constraints_block(voigt_status))
    print_dual(f"\n{EXTRA_FDF_FILE} (written into every generated folder):", f_out)
    for line in extra_fdf_text.rstrip("\n").split("\n"):
        print_dual(f"  {line}", f_out)
    print_dual(f"\n'%include {EXTRA_FDF_FILE}' is inserted right after '%include "
               f"{structure_basename}' in every generated {calc_basename} copy.", f_out)
    forced_calc_text = insert_include_after_structure(original_calc_text, structure_basename)

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
                "[OK] All required pseudopotentials found -- will be linked into every "
                "generated folder.", 'green'), f_out)
    else:
        print_dual("Not given (pass -p/--pseudo-dir -- a bundled bank or a folder path -- to "
                   "link the required pseudopotential for every species into every generated "
                   "folder). Pseudopotentials will need to be added manually.", f_out)

    print_section("[6] GENERATED STRAIN FOLDERS", f_out)

    # Generate strained structures
    n_steps = int(round((args.stmax - args.stmin) / args.step)) + 1
    strain_values = list(np.linspace(args.stmin, args.stmin + (n_steps - 1) * args.step, n_steps) / 100)

    folder_rows = []
    for strain in strain_values:
        # Handle negative strain (compression) in folder name
        strain_prefix = "m" if strain < 0 else ""
        folder = os.path.join(
            direction_dir, f"strain_{norm_dir}_{strain_prefix}{abs((strain * 100)):.2f}")
        os.makedirs(folder, exist_ok=True)

        try:
            new_vectors = apply_cartesian_strain(lattice_vectors, strain, norm_dir)
        except ValueError as e:
            fail(str(e))

        output_structure = os.path.join(folder, structure_basename)
        structure_io.rewrite_fdf_lattice(args.structure, new_vectors, output_structure)

        # calc.fdf (with the %include added above) and config_extra.fdf are
        # identical for every folder -- only the structure changes per strain step.
        with open(os.path.join(folder, calc_basename), "w") as fh:
            fh.write(forced_calc_text)
        with open(os.path.join(folder, EXTRA_FDF_FILE), "w") as fh:
            fh.write(extra_fdf_text)

        for sym in species:
            link_pseudo(args.pseudo_dir, sym, folder)

        volume = abs(np.linalg.det(new_vectors)) * structure.lattice_constant ** 3
        folder_rows.append(([folder, f"{strain * 100:+.2f}",
                            f"{structure_basename}, {calc_basename}, {EXTRA_FDF_FILE}",
                            f"{volume:.4f}"], None))

    print_table(["Folder", "Strain (%)", "Files", "Volume (Ang^3)"], folder_rows, f_out)

    print_section("[7] SUMMARY & NEXT STEPS", f_out)
    print_dual(f"{len(folder_rows)} folder(s) written under '{direction_dir}'.", f_out)
    if report_path:
        print_dual(f"Report            : {report_path}", f_out)
    print_dual(color_text("\nNext steps:", 'yellow'), f_out)
    print_dual(f"  1. Run SIESTA ({args.relax_mode} relaxation) in every "
               f"'{direction_dir}/strain_{norm_dir}_*/' folder.", f_out)
    print_dual(f"  2. Once they're done, run: stb-strainAnalysis --dir {direction_dir} (this "
               f"direction only), or --dir {args.output_dir} --compare (every direction found "
               "under it, side by side).", f_out)

    if f_out:
        f_out.close()

    print("\n[INFO] Complete job!")
    print("\n" + "-" * 60)
    print(color_text("Strain folders ready for Stage 2 (stb-strainAnalysis).\n", 'bold'))

if __name__ == "__main__":
    main()
