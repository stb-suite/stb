"""Group-theory classification of which Gamma-point phonon modes are
symmetry-forbidden from being Raman active -- lets stb-ramanModes
(--use-symmetry) skip writing Optical-calculation folders for those modes
entirely, since their Raman tensor is guaranteed exactly zero by symmetry,
not just numerically small.

Works for all 32 crystallographic point groups uniformly (not just a
centrosymmetric g/u shortcut): a mode's degenerate band-group is
Raman-active iff its irreducible representation is contained in the
REDUCIBLE representation spanned by the quadratic basis functions
(x^2, y^2, z^2, xy, xz, yz) -- the symmetric square of the ordinary
(x, y, z) vector representation. Standard result (see e.g. Cotton,
"Chemical Applications of Group Theory"):

    chi_quad(g) = (chi_v(g)^2 + chi_v(g^2)) / 2

and the multiplicity of a representation i (character chi_i) inside the
quadratic representation is the usual reduction-formula inner product

    n_i = (1/|G|) * sum_g chi_quad(g) * conj(chi_i(g))

n_i > 0 (up to numerical tolerance) means Raman-active; n_i == 0 means
symmetry-forbidden.

Implemented entirely from Phonopy's public IrReps API
(phonopy/phonon/irreps.py, verified by reading the installed source):
`characters`/`conventional_rotations`/`band_indices` are all per-POINT-
GROUP-ELEMENT arrays in lockstep order (NOT per conjugacy class -- see
IrReps._get_characters/_get_ground_matrix, which iterate the full little
group of q, one entry per operation). chi_v(g) is obtained directly from
whichever band-group(s) phonopy already assigned to the acoustic
(translational) bands -- no separate table of vector-representation
characters is needed, and no per-point-group Raman-active-label table is
hand-typed anywhere in this module: everything is derived algorithmically
from data Phonopy already computed. Verified against a live MACE-MP-0
silicon calculation during planning: the triply-degenerate optical mode
(11.7 THz, point group m-3m) came out n_i=1.0 (Raman-active), matching
the well-known literature result that Si's Raman mode is T2g; the
acoustic modes came out n_i=0.0 (correctly never Raman-active).

Mulliken labels (T2g, A1g, ...) for the report are read from IrReps'
internal `_ir_labels` attribute (not part of the public API) purely for
display -- never for the active/inactive decision itself, which only uses
the public `characters`/`conventional_rotations`/`band_indices`. If label
matching isn't available at all (e.g. Phonopy's own character-table
matching fails for some numerically-imperfect structure), `run_irreps`
itself raises before this module ever runs -- see classify_modes' single
broad try/except below, which degrades to "symmetry unavailable" rather
than guessing.
"""

from __future__ import annotations

import numpy as np

RAMAN_ACTIVE_TOL = 1e-3


class ModeSymmetry:
    """band_indices: the 0-based Phonopy band indices in this degenerate
    group. label: best-effort Mulliken symbol (e.g. "T2g"), or None if
    unavailable. is_raman_active: True/False when the group-theory
    decision is certain, None when unknown (never used to justify a
    skip -- callers must treat None the same as True, i.e. keep the
    mode).
    """

    def __init__(self, band_indices, label, is_raman_active):
        self.band_indices = list(band_indices)
        self.label = label
        self.is_raman_active = is_raman_active


def classify_modes(phonon):
    """Runs phonon.run_irreps([0, 0, 0]) and classifies every Gamma-point
    band-group by Raman activity.

    Returns (mode_symmetries, point_group, error):
      - mode_symmetries: list[ModeSymmetry], one entry per degenerate
        band-group INCLUDING the acoustic ones (callers exclude those
        the same way core.phonon_workflow.get_gamma_modes already does,
        by band index < 3 -- kept here so chi_v can be derived from
        them).
      - point_group: phonon.symmetry.pointgroup_symbol (e.g. "m-3m"),
        always available even on failure below.
      - error: None on success, or a short human-readable reason string
        (e.g. "non-primitive cell") when classification wasn't possible
        at all -- mode_symmetries is [] in that case; callers must treat
        this as "symmetry information unavailable", never as "no modes
        are Raman-active".

    Never raises: every failure mode of phonon.run_irreps (non-primitive
    cell, unsupported magnetic space group, internal character-table
    matching failure for a numerically-imperfect structure) -- plus
    anything unexpected in this module's own reduction-formula math below
    -- is caught and reported via `error` instead, since --use-symmetry
    must always be safe to enable and fall back to "run every mode"
    rather than crash Stage 2. Deliberately broad (Exception, not just
    RuntimeError): this is an opt-in performance feature, not a
    correctness-critical path, so "never take the whole Stage 2 run down"
    outweighs the usual preference for a narrow except clause.
    """
    point_group = phonon.symmetry.pointgroup_symbol

    try:
        irreps = phonon.run_irreps([0, 0, 0])
        mode_symmetries = _classify_from_irreps(irreps)
    except Exception as e:
        return [], point_group, str(e)

    return mode_symmetries, point_group, None


def _classify_from_irreps(irreps):
    characters = irreps.characters  # (n_groups, |G|)
    rotations = irreps.conventional_rotations  # (|G|, 3, 3)
    band_groups = irreps.band_indices  # list of lists of 0-based band idx
    n_elem = rotations.shape[0]
    labels = getattr(irreps, "_ir_labels", None)  # best-effort, display only

    # The 3 acoustic (translational) bands are always Phonopy's first 3
    # at Gamma -- same convention already used by
    # core.phonon_workflow.get_gamma_modes(exclude_acoustic=True). A
    # band-group is "acoustic" if it's entirely contained in {0, 1, 2}
    # (one 3D group for high-symmetry point groups, up to three separate
    # 1D/2D groups when x, y, z split across different irreps).
    acoustic_group_idx = [i for i, bset in enumerate(band_groups) if max(bset) < 3]

    chi_v = np.zeros(n_elem, dtype=complex)
    for i in acoustic_group_idx:
        chi_v += characters[i]

    # chi_v(g^2): find, for each element g, the index g' whose rotation
    # matrix equals rotations[g] @ rotations[g] (exact integer match).
    # match is never empty in correct operation -- group closure
    # guarantees R^2 is itself an element of the same (finite, closed)
    # point group -- but the explicit message beats a bare IndexError if
    # that assumption is ever wrong (caught by classify_modes' broad
    # except either way).
    sq_index = np.empty(n_elem, dtype=int)
    for g in range(n_elem):
        r_sq = rotations[g] @ rotations[g]
        match = np.flatnonzero((rotations == r_sq).all(axis=(1, 2)))
        assert len(match) > 0, f"rotation element {g} squared isn't in the point group -- unexpected"
        sq_index[g] = match[0]

    chi_quad = (chi_v ** 2 + chi_v[sq_index]) / 2.0

    mode_symmetries = []
    for i, bset in enumerate(band_groups):
        n_i = np.sum(chi_quad * np.conj(characters[i])) / n_elem
        is_active = bool(abs(n_i.real) > RAMAN_ACTIVE_TOL)
        label = labels[i] if labels else None
        mode_symmetries.append(ModeSymmetry(bset, label, is_active))

    return mode_symmetries
