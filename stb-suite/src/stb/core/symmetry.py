"""Shared symmetry-based cell reduction, used by stb-unitcell and stb-fetch.

Wraps pymatgen's SpacegroupAnalyzer to reduce a structure to its primitive
cell, its conventional cell, or a symmetry-refined version of the input cell.
"""

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

UNITCELL_MODES = ("primitive", "conventional", "refined")


def reduce_to_unitcell(structure, mode, symprec=1e-3, angle_tolerance=5.0):
    """Returns (new_structure, sga) for the requested `mode`:
      - primitive: smallest possible cell.
      - conventional: standardized, usually larger, cell.
      - refined: conventional-sized cell with atomic positions snapped to the
        detected symmetry -- cleans up numerical noise from relaxations or
        hand-built/CIF structures without changing which atoms are present.

    `sga` (the SpacegroupAnalyzer built from `structure`) is returned too, so
    callers can report the detected space group/point group/crystal system
    without re-running symmetry detection.

    Note: the output's atom order and coordinate origin are NOT guaranteed to
    match the input, in any mode -- spglib rebuilds the cell from the
    detected symmetry operations from scratch, free to pick any symmetry
    -equivalent origin/ordering.

    Raises ValueError (from SpacegroupAnalyzer) if symmetry detection fails,
    or if `mode` isn't one of UNITCELL_MODES.
    """
    if mode not in UNITCELL_MODES:
        raise ValueError(f"Unknown unitcell mode {mode!r}, must be one of {UNITCELL_MODES}.")

    sga = SpacegroupAnalyzer(structure, symprec=symprec, angle_tolerance=angle_tolerance)
    if mode == "primitive":
        new_structure = sga.get_primitive_standard_structure()
    elif mode == "conventional":
        new_structure = sga.get_conventional_standard_structure()
    else:
        new_structure = sga.get_refined_structure()
    return new_structure, sga
