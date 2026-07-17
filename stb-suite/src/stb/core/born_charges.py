"""SIESTA SystemLabel.BC (Born effective charge) parsing, shared by
stb-irAnalysis's bulk-path (3D periodic) IR intensity readout.

One small parsing module per external SIESTA output-file TYPE, same
precedent as core/dielectric.py (.EPSIMG) -- kept out of core/siesta_log.py,
which is scoped to .out calculation logs specifically, not SIESTA's other
auxiliary output files.

FORMAT (confirmed against a real SystemLabel.BC, SIESTA MaX-1.2/PSML
build, 2026-07 -- the x/y/z-labeled layout originally assumed here, from
SIESTA's Born-effective-charge tutorial/exercise material, does NOT
match: this SIESTA build's real .BC file has a bare "BC matrix" header
line followed directly by N*3 UNLABELED rows of 3 floats each, one 3x3
tensor per atom, IN FILE ORDER, no per-atom separator:

    BC matrix
           0.0000449      -0.0000023      -0.0000023
          -0.0000023       0.0000449      -0.0000023
          -0.0000023      -0.0000023       0.0000449
           0.0000656       0.0000106       0.0000106
           0.0000106       0.0000656       0.0000106
           0.0000106       0.0000106       0.0000657

read_born_charges() accepts EITHER this unlabeled layout or the
originally-assumed x/y/z-labeled one (SIESTA's own tutorial material
does show labels, so a different build/version may still emit them) --
it locates every run of 3 consecutive data rows that are either all
unlabeled or labeled x/y/z in order (ignoring any other header/label
text around them) and treats each run, IN FILE ORDER, as one atom's Z*
tensor. This degrades safely (a malformed/unexpected file yields None,
not a misparsed tensor); atom_labels stays best-effort-generic ("atom_i")
since there's still no confirmed real example of a per-atom species/index
label in this file.
"""

from __future__ import annotations

import re

import numpy as np

_ROW_RE = re.compile(
    r'^\s*(?:([xyzXYZ])\s+)?([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*$'
)


def read_born_charges(path):
    """Reads a SystemLabel.BC file: (atom_labels, Z_star).

    Z_star.shape == (n_atoms, 3, 3); Z_star[i, tau, beta] is the Born
    effective charge component d F_(i,tau) / d E_beta for atom i (force
    on atom i along Cartesian direction tau, per unit applied field along
    Cartesian direction beta) -- combine with a mode's own per-atom
    eigendisplacement (core.phonon_workflow.mode_eigendisplacement) via
    dmu_beta/dQ = sum_atom sum_tau Z_star[atom,tau,beta] * e[atom,tau] to
    get that mode's IR dipole derivative.

    atom_labels is a list of best-effort per-atom label strings (row
    order in the file -- see the format caveat in this module's
    docstring), same length as Z_star's first axis; only for display, not
    used in any downstream physics.

    Returns None on file-not-found, empty file, or if no valid 3x3 block
    is found at all -- a Born-charge run is always exactly ONE folder per
    stb-irModes run (not a sweep-tolerant loop over many folders like
    core.dielectric's per-displacement reader), so there's no benefit to
    a raise-and-let-caller-catch style here; the caller just reports
    "incomplete" the same way it would for a missing file.
    """
    try:
        with open(path) as f:
            lines = f.readlines()
    except OSError:
        return None

    rows = []  # (axis_letter_or_None, [v0, v1, v2])
    for line in lines:
        m = _ROW_RE.match(line)
        if m is None:
            continue
        axis = m.group(1).lower() if m.group(1) else None
        try:
            vals = [float(m.group(i)) for i in (2, 3, 4)]
        except ValueError:
            continue
        rows.append((axis, vals))

    tensors = []
    i = 0
    axis_order = ('x', 'y', 'z')
    while i + 3 <= len(rows):
        triplet = rows[i:i + 3]
        axes = tuple(axis for axis, _ in triplet)
        if axes == axis_order or axes == (None, None, None):
            tensor = np.array([vals for _, vals in triplet], dtype=float)
            tensors.append(tensor)
            i += 3
        else:
            i += 1

    if not tensors:
        return None

    atom_labels = [f"atom_{i}" for i in range(len(tensors))]
    return atom_labels, np.stack(tensors, axis=0)
