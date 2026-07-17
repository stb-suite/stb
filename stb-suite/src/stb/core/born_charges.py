"""SIESTA SystemLabel.BC (Born effective charge) parsing, shared by
stb-irAnalysis's bulk-path (3D periodic) IR intensity readout.

One small parsing module per external SIESTA output-file TYPE, same
precedent as core/dielectric.py (.EPSIMG) -- kept out of core/siesta_log.py,
which is scoped to .out calculation logs specifically, not SIESTA's other
auxiliary output files.

FORMAT CAVEAT (unlike core/dielectric.py's .EPSIMG parser, which was
verified against SIESTA's own Src/optical.F source): the exact .BC layout
below is reconstructed from SIESTA's official Born-effective-charge
tutorial/exercise material (BornCharge .true. + MD.TypeOfRun FC +
%block PolarizationGrids), not from a real .BC file read in this
environment. The documented example is one 3x3 tensor block per atom,
each block a column header line (x/y/z) followed by 3 rows labeled
x/y/z:

    BC matrix
                x              y              z
    x       2.7661178      0.0000107     -0.0000000
    y      -0.0085767      2.7538482      0.0000000
    z       0.0004767      0.0000000      0.7117852

but the atom-labeling convention (species symbol vs. bare index) was not
confirmed. read_born_charges() is written to be tolerant of that
uncertainty: it locates every run of 3 consecutive x/y/z-labeled data
rows anywhere in the file (ignoring any header/label text around them)
and treats each run, IN FILE ORDER, as one atom's Z* tensor -- this
degrades safely (a malformed/unexpected file yields None, not a
misparsed tensor) but the atom-label strings it returns are best-effort
only. VERIFY AGAINST A REAL SystemLabel.BC FILE before trusting this in
production; adjust the regex/labeling below once one is available.
"""

from __future__ import annotations

import re

import numpy as np

_ROW_RE = re.compile(
    r'^\s*([xyzXYZ])\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*'
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

    rows = []  # (axis_letter, [v0, v1, v2])
    for line in lines:
        m = _ROW_RE.match(line)
        if m is None:
            continue
        axis = m.group(1).lower()
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
        if tuple(axis for axis, _ in triplet) == axis_order:
            tensor = np.array([vals for _, vals in triplet], dtype=float)
            tensors.append(tensor)
            i += 3
        else:
            i += 1

    if not tensors:
        return None

    atom_labels = [f"atom_{i}" for i in range(len(tensors))]
    return atom_labels, np.stack(tensors, axis=0)
