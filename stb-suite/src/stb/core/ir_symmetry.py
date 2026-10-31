"""Group-theory classification of which Gamma-point phonon modes are
symmetry-forbidden from being IR active -- lets stb-irModes
(--use-symmetry) skip modes for the non-bulk (dipole-derivative) path
that group theory guarantees are IR-silent, since their dipole
derivative is exactly zero by symmetry, not just numerically small. (For
the bulk Born-effective-charge path, one equilibrium calculation already
covers every mode regardless -- this classification is informational
only there, see stb-irModes' own docstring.)

A mode's degenerate band-group is IR-active iff its irreducible
representation is contained in the VECTOR (translational) representation
itself (x, y, z) -- unlike Raman, which needs the quadratic-function
representation (the symmetric square of the vector representation).
This is genuinely SIMPLER than Raman's rule, not harder: the vector
representation's own character, chi_vec(g), is exactly
core.vibrational_symmetry.acoustic_chi_vector(irreps) -- the same
acoustic-band character Raman's classify_modes() already extracts and
then squares into chi_quad -- reused here directly, no squaring step at
all.

For centrosymmetric point groups this makes Raman- and IR-activity
mutually exclusive (the textbook "rule of mutual exclusion": Raman needs
a gerade/even-parity irrep, IR needs an ungerade/odd-parity one, and a
centrosymmetric group's irreps are strictly one or the other) -- not
hard-coded as a g/u shortcut, it falls out of the general reduction
formula automatically for all 32 point groups, same design stance as
core.raman_symmetry.classify_modes.
"""

from __future__ import annotations

from stb.core.vibrational_symmetry import (
    run_gamma_irreps, acoustic_chi_vector, reduction_multiplicity,
)

IR_ACTIVE_TOL = 1e-3


class IRModeSymmetry:
    """band_indices: the 0-based Phonopy band indices in this degenerate
    group. label: best-effort Mulliken symbol (e.g. "A1"), or None if
    unavailable. is_ir_active: True/False when the group-theory decision
    is certain, None when unknown (never used to justify a skip --
    callers must treat None the same as True, i.e. keep the mode). Kept
    as its own small class, parallel to but not derived from
    core.raman_symmetry.ModeSymmetry -- renaming that class's
    `is_raman_active` field to a generic `is_active` across all of
    Raman's already-working call sites isn't worth it for a 3-line class.
    """

    def __init__(self, band_indices, label, is_ir_active):
        self.band_indices = list(band_indices)
        self.label = label
        self.is_ir_active = is_ir_active


def classify_ir_modes(phonon):
    """Runs phonon.run_irreps([0, 0, 0]) and classifies every Gamma-point
    band-group by IR activity.

    Returns (mode_symmetries, point_group, error) -- identical contract
    to core.raman_symmetry.classify_modes: mode_symmetries is
    list[IRModeSymmetry] (including the acoustic band-groups, same
    reasoning as Raman's classifier -- callers exclude those by band
    index < 3, same as core.phonon_workflow.get_gamma_modes already
    does); point_group is always available; error is None on success, or
    a short human-readable reason when classification wasn't possible at
    all (mode_symmetries is [] in that case -- callers must treat this as
    "IR information unavailable", never as "no modes are IR-active").

    Never raises: every failure mode of phonon.run_irreps, plus anything
    unexpected in this module's own reduction-formula math, is caught and
    reported via `error` instead -- same "opt-in feature must always be
    safe to enable" reasoning as core.raman_symmetry.classify_modes.
    """
    irreps, point_group, error = run_gamma_irreps(phonon)
    if error is not None:
        return [], point_group, error

    try:
        mode_symmetries = _classify_from_irreps(irreps)
    except Exception as e:
        return [], point_group, str(e)

    return mode_symmetries, point_group, None


def _classify_from_irreps(irreps):
    characters = irreps.characters  # (n_groups, |G|)
    band_groups = irreps.band_indices  # list of lists of 0-based band idx
    n_elem = irreps.conventional_rotations.shape[0]
    labels = getattr(irreps, "_ir_labels", None)  # best-effort, display only

    chi_vec = acoustic_chi_vector(irreps)

    mode_symmetries = []
    for i, bset in enumerate(band_groups):
        n_i = reduction_multiplicity(chi_vec, characters[i], n_elem)
        is_active = bool(abs(n_i.real) > IR_ACTIVE_TOL)
        label = labels[i] if labels else None
        mode_symmetries.append(IRModeSymmetry(bset, label, is_active))

    return mode_symmetries
