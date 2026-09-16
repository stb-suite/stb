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
_MD_STEPS_RE = re.compile(r'MD\.Steps\s+\d+', re.IGNORECASE)
_MD_FCDISPL_RE = re.compile(r'MD\.FCDispl\s+[-+0-9.eE]+\s*\S*', re.IGNORECASE)
_BORNCHARGE_RE = re.compile(r'BornCharge\s+\S+', re.IGNORECASE)
_POLARIZATION_GRIDS_RE = re.compile(
    r'%block\s+PolarizationGrids.*?%endblock\s+PolarizationGrids',
    re.IGNORECASE | re.DOTALL,
)


def force_single_point(calc_text):
    """Substitutes/appends 'MD.TypeOfRun CG' and 'MD.Steps 0' onto
    calc_text (regex .subn, append if the tag isn't present at all rather
    than erroring, since an absent tag is a normal template, not a
    malformed one). Uses MD.Steps, not the older MD.NumCGsteps spelling
    (deprecated, no longer honored by current SIESTA) -- see
    structure_io.read_effective_md_steps for the same MD.Steps-preferred
    convention on the reading side. Moved here from neb.py once
    stackingfault.py became a second consumer -- same reasoning in both
    callers: each folder (a NEB image, a stacking-fault grid point) is an
    independent, single-point evaluation with no inter-folder coupling
    once written to disk, so letting SIESTA relax one on its own would
    silently move it away from the specific point/geometry being sampled.
    Forced unconditionally, no opt-out flag, in both callers.
    """
    new_text, count = _MD_RUNTYPE_RE.subn('MD.TypeOfRun          CG', calc_text)
    if count == 0:
        new_text += "\nMD.TypeOfRun          CG\n"
    new_text, count = _MD_STEPS_RE.subn('MD.Steps              0', new_text)
    if count == 0:
        new_text += "MD.Steps              0\n"
    return new_text


_SQRT2_INV = 2.0 ** -0.5

# Shared by stb-optical (writes these folders) and stb-opticalAnalysis
# (detects them back from each folder's own %block Optical.Vector, and
# reconstructs the off-diagonal eps_ij from the biaxial ones) -- a single
# source of truth for the direction-name -> Optical.Vector convention, so
# the two stages can never silently drift apart on it. 'xx'/'yy'/'zz' are
# the pure Cartesian unit axes (Optical.Vector = e_i, giving the diagonal
# dielectric-tensor component eps_ii directly). 'xy'/'xz'/'yz' are the
# BIAXIAL directions: SIESTA's Optical.Vector always computes a single
# -direction response n^T.eps.n, never a genuine off-diagonal component
# directly, so each biaxial direction uses the NORMALIZED BISECTOR of the
# two axes involved (e.g. (e_x+e_y)/sqrt(2) for 'xy') -- the standard
# trick to reach one: for a symmetric dielectric tensor, the response
# along that bisector expands to (eps_ii + eps_jj)/2 + eps_ij, so
# eps_ij = eps_(bisector ij) - (eps_ii + eps_jj) / 2 once eps_ii, eps_jj
# AND the bisector's own response are all known. Same reconstruction
# formula/convention as raman_analysis.py's own raman_tensor_full (Rij
# from a mixed n=(i+j)/sqrt(2) direction), reimplemented here (not
# imported) since it operates on eps(E) spectra rather than a single
# static Raman-tensor number per mode.
OPTICAL_DIRECTION_VECTORS = {
    "xx": (1.0, 0.0, 0.0),
    "yy": (0.0, 1.0, 0.0),
    "zz": (0.0, 0.0, 1.0),
    "xy": (_SQRT2_INV, _SQRT2_INV, 0.0),
    "xz": (_SQRT2_INV, 0.0, _SQRT2_INV),
    "yz": (0.0, _SQRT2_INV, _SQRT2_INV),
}

# Which 2 diagonal directions each biaxial direction needs alongside it
# to be reconstructed into a real eps_ij.
OPTICAL_OFFDIAG_PAIRS = {
    "xy": ("xx", "yy"),
    "xz": ("xx", "zz"),
    "yz": ("yy", "zz"),
}


def build_optical_block(mesh, broaden_ev, axis_vec, nbands=None):
    """The %block Optical.Mesh/Optical.Vector + OpticalCalculation T fdf
    stanza for one Optical.Vector direction -- SIESTA's interband/RPA
    dielectric-function machinery, verified via SIESTA's own official
    tutorial to work for fully 3D periodic bulk crystals (not just
    vacuum-padded systems, unlike %block ExternalElectricField).

    Moved here from raman_modes.py once optical.py (stb-optical) became a
    second consumer -- same "extract on second use" convention as
    force_single_point above. raman_modes.py now imports this instead of
    keeping its own copy.
    """
    lines = [
        "OpticalCalculation T",
        f"Optical.Broaden        {broaden_ev} eV",
        "%block Optical.Mesh",
        f"  {mesh[0]}  {mesh[1]}  {mesh[2]}",
        "%endblock Optical.Mesh",
        "%block Optical.Vector",
        f"  {axis_vec[0]:.4f}  {axis_vec[1]:.4f}  {axis_vec[2]:.4f}",
        "%endblock Optical.Vector",
    ]
    if nbands:
        lines.append(f"Optical.NumberOfBands  {nbands}")
    return "\n".join(lines) + "\n"


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
