"""Single source of truth for k-space / k-grid math.

Consolidates compute_monkhorts, which used to be duplicated identically
(save for its error handling) in kgrid.py, cohesive_energy.py and
inputfile.py.
"""

from __future__ import annotations

import math

import numpy as np


def compute_monkhorts(cella, cellb, cellc, k_density: float) -> list[int]:
    """Computes the reciprocal lattice vectors and Monkhorst-Pack divisions.

    Raises ValueError if the cell volume is (numerically) zero.
    """
    volume = np.dot(cella, np.cross(cellb, cellc))

    if abs(volume) < 1e-9:
        raise ValueError("Cell volume is zero. Check lattice vectors.")

    b1 = 2 * np.pi * np.cross(cellb, cellc) / volume
    b2 = 2 * np.pi * np.cross(cellc, cella) / volume
    b3 = 2 * np.pi * np.cross(cella, cellb) / volume

    lengths = [np.linalg.norm(b) for b in (b1, b2, b3)]
    divisions = [max(1, math.ceil(length / k_density)) for length in lengths]
    return divisions
