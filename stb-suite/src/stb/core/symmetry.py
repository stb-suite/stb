"""Shared symmetry-based cell reduction, used by stb-unitcell and stb-fetch.

Wraps pymatgen's SpacegroupAnalyzer to reduce a structure to its primitive
cell, its conventional cell, or a symmetry-refined version of the input cell.
"""

import numpy as np
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

UNITCELL_MODES = ("primitive", "conventional", "refined")


def find_inequivalent_sites(pmg_structure, symprec, filter_species=None):
    """Returns (sites, space_group_label) where `sites` is a list of
    (index, wyckoff_letter, multiplicity) -- one representative atom (0-based
    index) per symmetrically distinct site, via spglib's equivalent-atoms
    mapping. If filter_species is given, only representatives of that
    species are returned (their multiplicity still counts all symmetry
    -equivalent atoms of that site, filtered species or not, since
    spglib never groups atoms of different species together).

    Moved here once stb-hubbardu became a second consumer (alongside
    stb-defect --all-inequivalent-sites) -- both need the same "which sites
    are physically distinct vs. just symmetry copies of each other"
    classification.
    """
    sga = SpacegroupAnalyzer(pmg_structure, symprec=symprec)
    dataset = sga.get_symmetry_dataset()

    groups = {}
    for i, rep in enumerate(dataset.equivalent_atoms):
        groups.setdefault(rep, []).append(i)

    sites = []
    for rep in sorted(groups):
        symbol = pmg_structure[rep].specie.symbol
        if filter_species is not None and symbol != filter_species:
            continue
        sites.append((rep, dataset.wyckoffs[rep], len(groups[rep])))

    space_group = f"{dataset.international} (No. {dataset.number})"
    return sites, space_group


def equivalent_cartesian_axes(pmg_structure, symprec=1e-3, angle_tolerance=5.0):
    """Groups the 3 Cartesian axes (x=0, y=1, z=2) by point-group symmetry:
    axes in the same group are related by a rotation/reflection of the
    structure's point group, so straining along one axis must give an
    identical mechanical response to straining along the other (Neumann's
    principle). Returns (groups, point_group_symbol), where `groups` is a
    list of sets of axis indices partitioning {0, 1, 2}.

    Uses exact axis alignment -- does some point-group operation map unit
    vector i onto (+/-) unit vector j -- rather than the weaker "the linear
    elastic tensor happens to be isotropic in this plane" criterion. This
    matters for hexagonal/trigonal in-plane axes (e.g. x vs y, 60 degrees
    apart): linear elasticity is isotropic there, but this is deliberately
    NOT treated as an equivalence, since the large-strain/ideal-strength
    calculations this is meant to support (stb-strain) are exactly where
    higher-order anharmonic terms make e.g. zigzag vs armchair directions
    genuinely different in real 2D materials (a known effect, e.g. in
    graphene) -- verified live against a synthetic cubic structure (m-3m:
    x/y/z all grouped together) and the suite's hexagonal 2D test fixture
    (-6m2: x and y correctly NOT grouped).
    """
    sga = SpacegroupAnalyzer(pmg_structure, symprec=symprec, angle_tolerance=angle_tolerance)
    point_group = sga.get_point_group_symbol()
    ops = sga.get_point_group_operations(cartesian=True)

    parent = [0, 1, 2]

    def find(i):
        while parent[i] != i:
            i = parent[i]
        return i

    def union(i, j):
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[ri] = rj

    basis = np.eye(3)
    for op in ops:
        for i in range(3):
            mapped = op.rotation_matrix @ basis[i]
            if not np.isclose(np.linalg.norm(mapped), 1.0, atol=1e-3):
                continue
            for j in range(3):
                if i != j and np.allclose(np.abs(mapped), basis[j], atol=1e-3):
                    union(i, j)

    groups = {}
    for i in range(3):
        groups.setdefault(find(i), set()).add(i)
    return list(groups.values()), point_group


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
