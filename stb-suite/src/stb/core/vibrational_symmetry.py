"""Shared group-theory boilerplate for classifying Gamma-point phonon
modes by activity selection rule -- factored out once core.ir_symmetry
became a second consumer of the same phonon.run_irreps() plumbing
core.raman_symmetry already had (Raman needs the quadratic-function
representation, IR needs the vector/translational one; both start from
the exact same irreps object and the exact same acoustic-band character
extraction). Kept deliberately minimal: only the genuinely
representation-agnostic pieces live here. Raman's quadratic-specific
squaring (chi_quad = (chi_v^2+chi_v(g^2))/2, needing the g -> g^2 lookup)
and its component-reduction machinery (tensor_form) have no IR analog at
all and stay in core/raman_symmetry.py.
"""

from __future__ import annotations

import numpy as np


def run_gamma_irreps(phonon):
    """Wraps phonon.run_irreps([0, 0, 0]) in a broad try/except -- both
    Raman's and IR's mode classifiers need this exact same fragile call
    (raises on a non-primitive cell, unsupported magnetic space group, or
    an internal character-table matching failure for a numerically-
    imperfect structure). Returns (irreps, point_group, error): irreps is
    None (with error set to a short human-readable reason) on any
    failure; point_group (phonon.symmetry.pointgroup_symbol, e.g. "m-3m")
    is always available, even on failure.
    """
    point_group = phonon.symmetry.pointgroup_symbol
    try:
        irreps = phonon.run_irreps([0, 0, 0])
    except Exception as e:
        return None, point_group, str(e)
    return irreps, point_group, None


def acoustic_chi_vector(irreps):
    """chi_v(g): the character of the vector (translational) representation,
    one complex entry per point-group element -- obtained by summing the
    characters of whichever degenerate band-group(s) Phonopy assigned to
    the 3 acoustic bands (always Phonopy's first 3 bands at Gamma, same
    convention core.phonon_workflow.get_gamma_modes(exclude_acoustic=True)
    already uses). This already IS the vector-representation character
    directly -- translations transform exactly as x, y, z under the point
    group -- so core.ir_symmetry reuses it as-is (IR's selection rule),
    while core.raman_symmetry squares it into chi_quad (Raman's selection
    rule needs the quadratic-function representation instead).

    `characters`/`conventional_rotations`/`band_indices` are all per-
    POINT-GROUP-ELEMENT arrays in lockstep order (not per conjugacy
    class -- see Phonopy's IrReps._get_characters/_get_ground_matrix).
    """
    characters = irreps.characters  # (n_groups, |G|)
    band_groups = irreps.band_indices  # list of lists of 0-based band idx
    n_elem = irreps.conventional_rotations.shape[0]

    acoustic_group_idx = [i for i, bset in enumerate(band_groups) if max(bset) < 3]
    chi_v = np.zeros(n_elem, dtype=complex)
    for i in acoustic_group_idx:
        chi_v += characters[i]
    return chi_v


def reduction_multiplicity(chi_target, chi_i, n_elem):
    """n_i = (1/|G|) * sum_g chi_target(g) * conj(chi_i(g)) -- the
    standard character-theory reduction-formula inner product: how many
    times irrep i (character chi_i) is contained in whichever reducible
    representation chi_target describes. Shared by Raman
    (chi_target = chi_quad, the quadratic-function representation) and
    IR (chi_target = chi_v itself, the vector representation).
    """
    return np.sum(chi_target * np.conj(chi_i)) / n_elem
