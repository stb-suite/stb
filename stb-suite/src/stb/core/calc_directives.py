"""Regex substitute-or-append helpers for forcing a specific directive onto
a calc.fdf template's text (kept as a plain string throughout the suite,
since it's user-supplied and its exact formatting/comments are otherwise
preserved verbatim) -- shared by every WORKFLOW_TOOLS prep tool that needs
to write an INDEPENDENT (non-coupled) folder per candidate/image, where
letting SIESTA relax on its own would silently move the geometry away from
the specific point being sampled.
"""

import re

_MD_RUNTYPE_RE = re.compile(r'MD\.TypeOfRun\s+\S+', re.IGNORECASE)
_MD_NUMCGSTEPS_RE = re.compile(r'MD\.NumCGsteps\s+\d+', re.IGNORECASE)
_MD_FCDISPL_RE = re.compile(r'MD\.FCDispl\s+[-+0-9.eE]+\s*\S*', re.IGNORECASE)
_BORNCHARGE_RE = re.compile(r'BornCharge\s+\S+', re.IGNORECASE)
_POLARIZATION_GRIDS_RE = re.compile(
    r'%block\s+PolarizationGrids.*?%endblock\s+PolarizationGrids',
    re.IGNORECASE | re.DOTALL,
)


def force_single_point(calc_text):
    """Substitutes/appends 'MD.TypeOfRun CG' and 'MD.NumCGsteps 0' onto
    calc_text (regex .subn, append if the tag isn't present at all rather
    than erroring, since an absent tag is a normal template, not a
    malformed one). Moved here from neb.py once stackingfault.py became a
    second consumer -- same reasoning in both callers: each folder (a NEB
    image, a stacking-fault grid point) is an independent, single-point
    evaluation with no inter-folder coupling once written to disk, so
    letting SIESTA relax one on its own would silently move it away from
    the specific point/geometry being sampled. Forced unconditionally, no
    opt-out flag, in both callers.
    """
    new_text, count = _MD_RUNTYPE_RE.subn('MD.TypeOfRun          CG', calc_text)
    if count == 0:
        new_text += "\nMD.TypeOfRun          CG\n"
    new_text, count = _MD_NUMCGSTEPS_RE.subn('MD.NumCGsteps         0', new_text)
    if count == 0:
        new_text += "MD.NumCGsteps         0\n"
    return new_text


def force_born_charge_run(calc_text, fc_displ_bohr, polarization_grid):
    """Substitutes/appends the fdf stanza SIESTA's OWN native Born-
    effective-charge automation needs -- MD.TypeOfRun FC (triggers
    SIESTA's internal finite-displacement engine, distinct from the
    phonopy-external displacement scheme used everywhere else in this
    suite), MD.FCDispl <fc_displ_bohr> bohr, BornCharge T, and a
    %block PolarizationGrids ... %endblock (the Berry-phase k-point grids
    that make macroscopic-polarization/Born-charge evaluation happen at
    all -- SIESTA only computes it when this block is present). Used by
    stb-irModes' bulk (3D periodic) path -- see core/born_charges.py for
    the resulting SystemLabel.BC file this run produces.

    `polarization_grid`: a flat sequence of 9 ints, row-major 3x3 (one
    row per lattice vector -- diagonal entries are 1D line-integral
    k-point counts, off-diagonal entries are the 2D surface-integral
    mesh, per SIESTA's own documentation). Same substitute-or-append
    convention as force_single_point: an absent tag/block is a normal
    template, not a malformed one.

    NOTE: this fdf combination has not been verified empirically against
    a real SIESTA run in this environment (see core/born_charges.py's own
    format caveat) -- flagged during planning as needing confirmation
    that MD.TypeOfRun FC/BornCharge/PolarizationGrids coexist cleanly
    with an otherwise CG/single-point-tuned calc.fdf template.
    """
    new_text, count = _MD_RUNTYPE_RE.subn('MD.TypeOfRun          FC', calc_text)
    if count == 0:
        new_text += "\nMD.TypeOfRun          FC\n"

    fcdispl_line = f"MD.FCDispl            {fc_displ_bohr} bohr"
    new_text, count = _MD_FCDISPL_RE.subn(fcdispl_line, new_text)
    if count == 0:
        new_text += fcdispl_line + "\n"

    new_text, count = _BORNCHARGE_RE.subn('BornCharge            T', new_text)
    if count == 0:
        new_text += "BornCharge            T\n"

    grid = list(polarization_grid)
    if len(grid) != 9:
        raise ValueError(
            f"polarization_grid must have exactly 9 values (3x3, row-major), got {len(grid)}")
    grid_lines = "\n".join(f"  {grid[3 * r]}  {grid[3 * r + 1]}  {grid[3 * r + 2]}" for r in range(3))
    block_text = f"%block PolarizationGrids\n{grid_lines}\n%endblock PolarizationGrids"
    new_text, count = _POLARIZATION_GRIDS_RE.subn(block_text, new_text)
    if count == 0:
        new_text += "\n" + block_text + "\n"

    return new_text
