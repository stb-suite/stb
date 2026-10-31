"""Single source of truth for k-space / k-grid math.

Consolidates compute_monkhorts, which used to be duplicated identically
(save for its error handling) in kgrid.py, cohesive_energy.py and
inputfile.py.
"""

from __future__ import annotations

import math

import numpy as np

from stb.core.cli import print_dual


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


def find_surface_reference(fracs):
    """Returns (frac_start, gap_size_frac): the fractional coordinate of the
    atom immediately BELOW the largest circular (wrapped) gap along one
    axis, and how far into it the search should go -- i.e. "where does the
    real vacuum region start, and how far does it extend" along this axis.

    Added for stb-stm, which (unlike `detect_vacuum_axes` above) needs to
    know not just THAT an axis is vacuum-padded but WHERE the vacuum
    actually is relative to the atoms, to correctly report "height above
    the surface". A naive `xyz[:, axis].max()` (the previous approach)
    silently picks the WRONG bounding atom whenever a structure's atoms are
    stored straddling the periodic cell boundary with the real vacuum
    gap in the middle of the cell instead of padded after the atoms (e.g.
    some externally-fetched database structures, as opposed to ones this
    suite's own stb-slab builds) -- verified live on a real fetched CrS
    monolayer (atoms at fractional z = 0, 0, 0.066, 0.934): the largest
    real gap (~87% of the cell, between the two z=0.066/0.934 atoms) was
    silently missed in favor of the tiny ~7% wraparound sliver beyond the
    naive "topmost" (z=0.934) atom, collapsing the whole search window to
    ~1.5 Ang instead of the genuine ~20 Ang vacuum region.

    `gap_size_frac` is always capped at HALF the identified gap's own
    fractional size. A periodic cell's stacking axis is topologically a
    ring: ANY single compact atomic region surrounded by one vacuum gap
    has TWO faces exposed to that same gap (conventionally "top" and
    "bottom" of the slab), not one -- searching the FULL gap risks the
    outside-in scan eventually crossing into the far face's own LDOS tail
    and reporting a nonsense mix of both surfaces as one. Verified live on
    the same real CrS monolayer: searching the full ~20 Ang gap gave
    wildly unphysical "corrugation" values (18-20 Ang, vs. real STM
    corrugation of sub-Angstrom to a few Ang) at several representative
    --iso values, all fixed by capping at half (~10 Ang) instead. This
    matches this tool's own documented, pre-existing limitation ("only
    images the surface exposed in the +axis direction... a slab with two
    exposed faces only gets its 'top' one imaged this way") -- the cap
    simply enforces that limitation numerically instead of relying on the
    atom arrangement to accidentally already imply it (as a perfectly
    centered single-plane structure like this suite's own graphene fixture
    happens to).
    """
    fracs_sorted = np.sort(np.asarray(fracs, dtype=float) % 1.0)
    if len(fracs_sorted) <= 1:
        return (float(fracs_sorted[0]) if len(fracs_sorted) else 0.0), 0.5

    if len(np.unique(np.round(fracs_sorted, 8))) == 1:
        return float(fracs_sorted[0]), 0.5

    gaps = np.diff(fracs_sorted)
    wraparound = 1.0 - fracs_sorted[-1] + fracs_sorted[0]
    idx = int(np.argmax(gaps))
    if wraparound > gaps[idx]:
        frac_start, gap_size = float(fracs_sorted[-1]), float(wraparound)
    else:
        frac_start, gap_size = float(fracs_sorted[idx]), float(gaps[idx])
    return frac_start, min(gap_size, 0.5)


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


def dimensionality_label(vacuum_axes) -> str:
    """Classifies a system as 0D/1D/2D/3D from which axes were flagged as
    vacuum-padded by detect_vacuum_axes (not from k-grid divisions
    themselves -- a genuinely periodic but large-celled axis can also round
    down to a single division, and shouldn't be mislabeled as vacuum).
    Shared by stb-kgrid and stb-strain (moved here once the second real
    consumer needed it, same extract-on-second-use policy as the rest of
    core/).
    """
    vacuum_count = sum(vacuum_axes)
    if vacuum_count == 3:
        return "0D (e.g., a molecule)"
    elif vacuum_count == 2:
        return "1D (e.g., a nanotube or polymer)"
    elif vacuum_count == 1:
        return "2D (e.g., a slab or surface)"
    return "3D (bulk material)"


def print_density_recommendation(f_out=None) -> None:
    """Prints a k-point density recommendation table -- shared by stb-kgrid's
    own report and stb-suite's interactive 1.2 wrapper (which shows it BEFORE
    prompting for a density, so the choice is informed). Moved here once the
    interactive wrapper became a second consumer, alongside stb-kgrid's own
    CLI report, of what used to be a kgrid.py-local helper.

    Uses print_dual so it participates correctly in --save-report when a
    caller passes f_out; the interactive wrapper calls it with no f_out for
    a plain console-only display.
    """
    print_dual("\nK-Point Density Recommendation Guide", f_out)
    print_dual("-" * 60, f_out)
    print_dual("  Density (1/Ang)        Accuracy Level", f_out)
    print_dual("  ---------------      --------------------------", f_out)
    print_dual("  0.05 - 0.1           High precision", f_out)
    print_dual("  0.10 - 0.30          Medium precision", f_out)
    print_dual("  0.30 - 0.50          Low precision", f_out)
    print_dual("", f_out)
    print_dual("  Tip: for most systems, a density between 0.2 and 0.3 is", f_out)
    print_dual("  generally accurate enough while keeping cost reasonable.", f_out)
