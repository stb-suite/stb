"""Setup-level facts about an optical calculation, shared by stb-optical (Stage 1) and
stb-opticalAnalysis (Stage 2): the k-point spacing a mesh gives, a rule-of-thumb
fine/moderate/coarse band for it, and what crystal symmetry says about the dielectric tensor
(independent components, components locked equal or zero, the smallest direction set).
"""

import numpy as np

from stb.core import structure_io
from stb.core.calc_directives import OPTICAL_OFFDIAG_PAIRS as _OFFDIAG_REQUIRES

_COMPONENTS = ["xx", "yy", "zz", "xy", "xz", "yz"]
_COMPONENT_INDEX = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
# Mesh-spacing bands (1/Ang, 2*pi included) for the advisory in [3]: a rule of thumb from
# converging bulk Si (README of examples/5.14-optical), not a guarantee for any system.
_SPACING_FINE, _SPACING_COARSE = 0.07, 0.12



def reciprocal_spacings(lattice, mesh):
    """|b_i| / n_i in 1/Ang (2*pi included) for each axis: the k-point spacing a mesh gives."""
    recip = 2.0 * np.pi * np.linalg.inv(np.asarray(lattice, dtype=float)).T
    return [float(np.linalg.norm(recip[i]) / mesh[i]) for i in range(3)]


def spacing_band(spacing):
    return "fine" if spacing <= _SPACING_FINE else ("moderate" if spacing <= _SPACING_COARSE else "coarse")


def symmetry_analysis(fdf_structure, symprec=0.01):
    """What crystal symmetry says about the dielectric tensor, in the CARTESIAN frame of
    the cell as given: space group, crystal system, the number of independent components,
    which Cartesian components can be non-zero, which are locked equal, and the smallest
    set of requested directions that determines the whole tensor. The tensor is projected
    onto the subspace invariant under every symmetry operation (Reynolds operator), so it
    holds for any orientation of the cell, not only axis-aligned ones. None if the
    analysis fails (degenerate cell, ...)."""
    try:
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
        sga = SpacegroupAnalyzer(structure_io.to_pymatgen(fdf_structure), symprec=symprec)
        rotations = [op.rotation_matrix for op in sga.get_symmetry_operations(cartesian=True)]
        info = {"symbol": sga.get_space_group_symbol(), "number": sga.get_space_group_number(),
                "system": sga.get_crystal_system(), "n_ops": len(rotations), "symprec": symprec}
    except Exception:
        return None
    projector = np.zeros((6, 6))
    for k, (i, j) in enumerate(_COMPONENT_INDEX):
        basis = np.zeros((3, 3))
        basis[i, j] = basis[j, i] = 1.0
        avg = sum(r @ basis @ r.T for r in rotations) / len(rotations)
        projector[:, k] = [avg[a, b] for a, b in _COMPONENT_INDEX]
    u, sv, _ = np.linalg.svd(projector)
    rank = int(np.sum(sv > 1e-8))
    inv = u[:, :rank]
    tol = 1e-6
    zero = [bool(np.max(np.abs(inv[c])) < tol) for c in range(6)]
    classes = []
    for c in range(6):
        if zero[c]:
            continue
        for cls in classes:
            if np.allclose(inv[c], inv[cls[0]], atol=tol):
                cls.append(c)
                break
        else:
            classes.append([c])
    needed = []
    for cls in classes:
        rep = cls[0]
        needed.append(rep)
        if rep >= 3:                                    # an off-diagonal needs its diagonal pair
            for diag in _OFFDIAG_REQUIRES[_COMPONENTS[rep]]:
                if _COMPONENTS.index(diag) not in needed:
                    needed.append(_COMPONENTS.index(diag))
    info.update({
        "independent": rank,
        "zero": [_COMPONENTS[c] for c in range(6) if zero[c]],
        "classes": [[_COMPONENTS[c] for c in cls] for cls in classes],
        "minimal": [_COMPONENTS[c] for c in sorted(set(needed))],
    })
    return info
