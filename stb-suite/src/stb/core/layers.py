"""Atomic layers of a slab: which atoms form each plane parallel to the surface.
Shared by stb-slab (--layers, the layer report) and stb-adsorb (--fix-layers)."""

import numpy as np


def surface_normal(lattice_matrix):
    """Unit normal of the surface (the a-b plane), on the side c points to."""
    a, b, c = np.asarray(lattice_matrix, dtype=float)
    normal = np.cross(a, b)
    normal /= np.linalg.norm(normal)
    return normal if np.dot(c, normal) >= 0 else -normal


def joined_frac_c(frac_c):
    """Fractional c of every atom with the slab in one piece: every atom above
    the widest empty gap along c -- the vacuum -- is brought down by one c."""
    frac_c = np.array(frac_c, dtype=float)
    frac_c -= np.floor(frac_c)
    order = np.argsort(frac_c)
    heights = frac_c[order]
    gaps = np.diff(np.append(heights, heights[0] + 1.0))
    frac_c[order[int(np.argmax(gaps)) + 1:]] -= 1.0
    return frac_c


# Atoms whose heights differ by less than this (Angstrom) belong to one atomic layer.
LAYER_TOL_ANG = 0.1


def atomic_layers(slab, tol=LAYER_TOL_ANG):
    """Atomic layers of a slab, bottom to top: a list of (height, atom indices),
    height in Angstrom along the surface normal (the layer's mean). Atoms are
    sorted by height and a new layer starts at every step larger than `tol`, so a
    buckled plane (the two halves of a Si(111) bilayer, 0.8 Ang apart) counts as
    two layers."""
    matrix = np.asarray(slab.lattice.matrix, dtype=float)
    c_height = float(np.dot(matrix[2], surface_normal(matrix)))
    heights = joined_frac_c(slab.frac_coords[:, 2]) * c_height
    order = np.argsort(heights, kind="stable")
    groups = [[order[0]]]
    for previous, index in zip(order[:-1], order[1:]):
        if heights[index] - heights[previous] > tol:
            groups.append([])
        groups[-1].append(index)
    return [(float(np.mean(heights[g])), np.array(g)) for g in groups]
