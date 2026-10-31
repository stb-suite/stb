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

from stb.core.vibrational_symmetry import (
    run_gamma_irreps, acoustic_chi_vector, reduction_multiplicity,
)

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
    rotations = irreps.conventional_rotations  # (|G|, 3, 3)
    band_groups = irreps.band_indices  # list of lists of 0-based band idx
    n_elem = rotations.shape[0]
    labels = getattr(irreps, "_ir_labels", None)  # best-effort, display only

    chi_v = acoustic_chi_vector(irreps)

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
        n_i = reduction_multiplicity(chi_quad, characters[i], n_elem)
        is_active = bool(abs(n_i.real) > RAMAN_ACTIVE_TOL)
        label = labels[i] if labels else None
        mode_symmetries.append(ModeSymmetry(bset, label, is_active))

    return mode_symmetries


# --- Component-level tensor-form reduction (stb-ramanModes --reduce-components) ---
#
# Symmetric-tensor basis order used throughout this section: xx, yy, zz,
# xy, xz, yz (matches _OFFDIAG_PAIRS' naming in raman_analysis.py).
_QUAD_INDICES = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
_QUAD_AXIS_NAMES = ["xx", "yy", "zz", "xy", "xz", "yz"]

# The 6 Optical.Vector probe directions stb-ramanModes can write folders
# for (x/y/z always; xy/xz/yz only with --full-tensor) -- same unit
# vectors as raman_modes.py's _AXES_DIAGONAL/_AXES_OFFDIAG, duplicated
# here (not imported) to keep this module independent of the CLI tool.
_INV_SQRT2 = 0.70710678118654752440
_PROBE_VECTORS = {
    "x": np.array([1.0, 0.0, 0.0]), "y": np.array([0.0, 1.0, 0.0]), "z": np.array([0.0, 0.0, 1.0]),
    "xy": np.array([_INV_SQRT2, _INV_SQRT2, 0.0]),
    "xz": np.array([_INV_SQRT2, 0.0, _INV_SQRT2]),
    "yz": np.array([0.0, _INV_SQRT2, _INV_SQRT2]),
}


def _vector_to_symmetric_matrix(v):
    Q = np.zeros((3, 3))
    for k, (i, j) in enumerate(_QUAD_INDICES):
        Q[i, j] = v[k]
        Q[j, i] = v[k]
    return Q


def _symmetric_matrix_to_vector(Q):
    return np.array([Q[i, j] for i, j in _QUAD_INDICES])


def _quad_representation_matrix(rotation):
    """6x6 matrix acting on the (xx,yy,zz,xy,xz,yz) basis induced by the
    3x3 Cartesian rotation `rotation`: a quadratic form (as a 3x3
    symmetric matrix Q) transforms as Q' = R Q R^T under r' = R r. Built
    by direct probing (each of the 6 basis quadratic functions run
    through the transform) rather than a hand-derived index formula --
    simpler to get right, and cross-checked during development: its
    trace exactly reproduces (to 1e-15) the (chi_v(g)^2+chi_v(g^2))/2
    formula already used above for the same chi_quad(g), on a real
    MACE-MP-0 silicon calculation.
    """
    D = np.zeros((6, 6))
    for k in range(6):
        v = np.zeros(6)
        v[k] = 1.0
        transformed = rotation @ _vector_to_symmetric_matrix(v) @ rotation.T
        D[:, k] = _symmetric_matrix_to_vector(transformed)
    return D


def tensor_form(phonon, band_index, full_tensor=False):
    """For a Gamma-point mode that is (a) NOT part of a degenerate
    band-group (its own group has exactly 1 band -- a 1D irrep, no
    "which partner" ambiguity) and (b) has multiplicity exactly 1 in the
    quadratic-function representation (n_i == 1: the Raman tensor's SHAPE
    is fixed up to one overall scale by symmetry alone, no further
    ambiguity), returns a dict describing that fixed shape and the single
    probe direction needed to fully determine it:

        {"T0": 3x3 array (the invariant tensor shape, unit-normalized),
         "probe_axis": str, one of "x"/"y"/"z"/"xy"/"xz"/"yz"}

    Every OTHER component of the true Raman tensor is then exactly
    `measured_value_at_probe_axis / (n_probe^T T0 n_probe) * T0`, i.e. a
    single +/-delta measurement at `probe_axis` is enough to reconstruct
    the entire tensor (stb-ramanModes only needs to write that one pair
    of folders; stb-ramanAnalysis reconstructs the rest algebraically).

    Returns None when the mode doesn't qualify for this reduction:
    degenerate (band-group size != 1), inactive (n_i == 0 -- already
    handled by --use-symmetry, nothing to probe at all), multiplicity
    >= 2 (genuinely needs more than one independent measurement -- not
    reduced by this function, falls back to full probing), or any
    unexpected failure (mirrors classify_modes' own broad safety net --
    this is an opt-in optimization, never worth crashing Stage 2 over).

    `full_tensor`: when False (stb-ramanModes' default diagonal-only
    scope), `probe_axis` is chosen only from {"x","y","z"}; when True,
    also considers the face-diagonal directions {"xy","xz","yz"} (whichever
    has the largest |n^T T0 n|, for the best-conditioned single
    calibration measurement) -- matches whichever probe set the caller
    is actually going to write folders for.
    """
    try:
        irreps = phonon.run_irreps([0, 0, 0])
        band_groups = irreps.band_indices
        group_idx = next((i for i, bset in enumerate(band_groups) if band_index in bset), None)
        if group_idx is None or len(band_groups[group_idx]) != 1:
            return None

        characters = irreps.characters
        rotations = irreps.conventional_rotations
        n_elem = rotations.shape[0]
        chi_i = characters[group_idx]

        projector = np.zeros((6, 6), dtype=complex)
        for g in range(n_elem):
            projector += np.conj(chi_i[g]) * _quad_representation_matrix(rotations[g])
        projector = (projector / n_elem).real

        eigvals, eigvecs = np.linalg.eigh(projector)
        n_i = int(round(float(np.sum(eigvals[eigvals > 0.5]))))
        if n_i != 1:
            return None

        v0 = eigvecs[:, int(np.argmax(eigvals))]
        t0 = _vector_to_symmetric_matrix(v0)
        # Deterministic sign/normalization: make the largest-|value| entry
        # of T0 positive (T0 and -T0 are equally valid basis vectors of a
        # 1-dim real subspace -- fixing a convention keeps reruns stable).
        flat = t0.flatten()
        if flat[np.argmax(np.abs(flat))] < 0:
            t0 = -t0

        axis_names = ("x", "y", "z") if not full_tensor else ("x", "y", "z", "xy", "xz", "yz")
        best_axis, best_overlap = None, 0.0
        for name in axis_names:
            n = _PROBE_VECTORS[name]
            overlap = abs(float(n @ t0 @ n))
            if overlap > best_overlap:
                best_axis, best_overlap = name, overlap
        if best_axis is None or best_overlap < 1e-8:
            return None

        return {"T0": t0, "probe_axis": best_axis}
    except Exception:
        return None
