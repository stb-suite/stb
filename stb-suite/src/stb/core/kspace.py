"""Single source of truth for k-space / k-grid math.

Consolidates compute_monkhorts, which used to be duplicated identically
(save for its error handling) in kgrid.py, cohesive_energy.py and
inputfile.py.
"""

from __future__ import annotations

import math

import numpy as np


def compute_monkhorts(cella, cellb, cellc, k_density: float, vacuum_axes=None) -> list[int]:
    """Computes the reciprocal lattice vectors and Monkhorst-Pack divisions.

    `vacuum_axes`, if given, is a length-3 sequence of bools (see
    `detect_vacuum_axes`): axes flagged True are forced to a single
    (Gamma-only) division regardless of k_density, since there is no real
    periodicity to sample along a vacuum-padded direction. Axes flagged False
    (or all axes, if vacuum_axes is None) use the density-based formula.

    Raises ValueError if the cell volume is (numerically) zero, or if
    k_density isn't a positive number (zero overflows math.ceil, negative
    silently returns a meaningless [1, 1, 1] grid).
    """
    if k_density <= 0:
        raise ValueError(f"k_density must be positive, got {k_density}.")

    volume = np.dot(cella, np.cross(cellb, cellc))

    if abs(volume) < 1e-9:
        raise ValueError("Cell volume is zero. Check lattice vectors.")

    b1 = 2 * np.pi * np.cross(cellb, cellc) / volume
    b2 = 2 * np.pi * np.cross(cellc, cella) / volume
    b3 = 2 * np.pi * np.cross(cella, cellb) / volume

    lengths = [np.linalg.norm(b) for b in (b1, b2, b3)]
    if vacuum_axes is None:
        vacuum_axes = (False, False, False)
    divisions = [
        1 if vacuum_axes[i] else max(1, math.ceil(lengths[i] / k_density))
        for i in range(3)
    ]
    return divisions


def _largest_circular_gap(fracs) -> float:
    """Largest empty span (fractional units, in [0, 1]) between points wrapped
    on a periodic ring. A single point (or none) yields 1.0 -- the whole cell
    is "empty" relative to it, which is the correct vacuum reading for an
    isolated atom along every axis.
    """
    fracs = np.sort(np.asarray(fracs, dtype=float) % 1.0)
    if len(fracs) <= 1:
        return 1.0
    gaps = np.diff(fracs)
    wraparound = 1.0 - fracs[-1] + fracs[0]
    return float(max(gaps.max(), wraparound))


def to_fractional(positions, lattice, is_cartesian: bool):
    """Converts atomic positions to fractional coordinates.

    `positions` is already fractional if `is_cartesian` is False (returned
    as-is). Otherwise it is treated as Cartesian Angstrom and converted via
    the inverse of `lattice` (rows are the 3 lattice vectors, Angstrom).
    """
    positions = np.asarray(positions, dtype=float)
    if not is_cartesian:
        return positions
    return positions @ np.linalg.inv(np.asarray(lattice, dtype=float))


def detect_vacuum_axes(frac_coords, lattice, vacuum_gap: float) -> list[bool]:
    """Flags, for each of the 3 lattice directions, whether atoms leave a gap
    of at least `vacuum_gap` Angstrom with nothing in it (wrapped
    periodically). That gap means the direction has no real periodicity to
    sample -- e.g. a slab's out-of-plane axis, a wire's 2 transverse axes, or
    all 3 axes for an isolated molecule in a box -- as opposed to a genuinely
    periodic axis that merely happens to be long.
    """
    frac_coords = np.asarray(frac_coords, dtype=float)
    axes_vacuum = []
    for i in range(3):
        gap_frac = _largest_circular_gap(frac_coords[:, i])
        gap_ang = gap_frac * np.linalg.norm(lattice[i])
        axes_vacuum.append(gap_ang >= vacuum_gap)
    return axes_vacuum
