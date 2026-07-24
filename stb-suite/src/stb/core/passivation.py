"""Dangling-bond passivation, used by stb-passivate and stb-slab's --passivate.

Caps undercoordinated surface atoms (e.g. from a cut slab) with a passivating
atom (usually H) along the missing-bond direction, determined purely from
local coordination geometry -- no bulk-lattice assumptions hard-coded in.
"""

import statistics
from collections import defaultdict

import numpy as np
from pymatgen.core import Structure
from pymatgen.core.periodic_table import Element


def _auto_cutoff(structure):
    """1.25x the shortest pairwise distance in the structure -- adapts to
    whatever material/bond length is given."""
    dm = structure.distance_matrix
    np.fill_diagonal(dm, np.inf)
    return 1.25 * dm.min()


def _auto_bond_length(species_symbol, passivant):
    """Sum of atomic radii for the two species -- an approximation (e.g. Si +
    H = 1.35 Ang vs. the experimental Si-H ~1.48 Ang), good enough as a
    default; callers can override.
    """
    return Element(species_symbol).atomic_radius + Element(passivant).atomic_radius


def passivate_dangling_bonds(structure, passivant="H", cutoff=None, bond_length=None):
    """Finds undercoordinated atoms and caps single-missing-bond (deficit == 1)
    sites with `passivant`, placed along the direction opposite the sum of
    unit vectors to existing neighbors (the exact missing-bond direction for
    an ideal local geometry, e.g. tetrahedral -- confirmed empirically against
    a real Si(111) slab: recovers the 109.5 degree bond angle exactly).

    Atoms with deficit >= 2 are NOT auto-passivated -- with only 1 or fewer
    existing bond vectors, the remaining directions are geometrically
    underdetermined without extra symmetry assumptions this function doesn't
    make. They're reported instead of guessing.

    Args:
        structure: pymatgen Structure (e.g. a cut slab).
        passivant: element symbol to cap dangling bonds with (default "H").
        cutoff: neighbor-search radius in Angstrom. Default: auto-detected as
            1.25x the shortest pairwise distance in `structure`.
        bond_length: passivant bond length in Angstrom. Default: auto,
            per (surface species, passivant) pair, as the sum of their
            pymatgen atomic radii.

    Returns:
        (new_structure, report) where report is:
        {
            "passivated": [(atom_index, species_symbol, position), ...],
            "unresolved": [(atom_index, species_symbol, position, deficit), ...],
        }
    """
    if cutoff is None:
        cutoff = _auto_cutoff(structure)

    coord_numbers = []
    neighbor_lists = []
    for site in structure:
        neighbors = structure.get_neighbors(site, cutoff)
        coord_numbers.append(len(neighbors))
        neighbor_lists.append(neighbors)

    by_species = defaultdict(list)
    for site, cn in zip(structure, coord_numbers):
        by_species[site.specie.symbol].append(cn)
    reference_coord = {symbol: statistics.mode(cns) for symbol, cns in by_species.items()}

    report = {"passivated": [], "unresolved": []}
    new_sites = []
    for i, (site, cn, neighbors) in enumerate(zip(structure, coord_numbers, neighbor_lists)):
        deficit = reference_coord[site.specie.symbol] - cn
        if deficit <= 0:
            continue
        if deficit >= 2:
            report["unresolved"].append((i, site.specie.symbol, site.coords, deficit))
            continue

        vectors = [(n.coords - site.coords) / np.linalg.norm(n.coords - site.coords) for n in neighbors]
        missing_dir = -np.sum(vectors, axis=0)
        missing_dir_norm = np.linalg.norm(missing_dir)
        if missing_dir_norm < 1e-6:
            # The existing neighbor directions cancel out (near-)exactly, e.g. an
            # atypical/non-physical coordination geometry -- the missing-bond
            # direction is undetermined from local geometry alone, same as the
            # deficit >= 2 case above. Reporting this instead of dividing by
            # ~0 avoids silently placing the passivant at a NaN position, which
            # otherwise crashes downstream symmetry analysis (spglib segfaults
            # on NaN coordinates rather than raising a catchable exception).
            report["unresolved"].append((i, site.specie.symbol, site.coords, deficit))
            continue
        missing_dir /= missing_dir_norm

        length = bond_length if bond_length is not None else _auto_bond_length(site.specie.symbol, passivant)
        new_pos = site.coords + missing_dir * length
        new_sites.append(new_pos)
        report["passivated"].append((i, site.specie.symbol, new_pos))

    new_structure = structure.copy()
    for pos in new_sites:
        new_structure.append(passivant, pos, coords_are_cartesian=True)

    return new_structure, report
