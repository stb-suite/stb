"""Shared MACE-MP-0 relaxation helpers, used by stb-mlrelax, stb-defect's
--ml-rank, and stb-amorphize. Callers must call core.deps.require_mace()
themselves first -- this module assumes mace/torch are already available,
and only imports them lazily inside each function (not at module level),
so merely importing this module never forces the heavy PyTorch/mace
dependency chain to load.
"""

# ASE/Voigt strain order: [xx, yy, zz, yz, xz, xy]. Maps each shear component
# to the pair of axes it mixes, so any vacuum axis also fixes the shears that
# would tilt it into the periodic directions.
_VOIGT_SHEAR_AXES = {3: (1, 2), 4: (0, 2), 5: (0, 1)}


def build_cell_mask(vacuum_axes):
    """Returns a 6-element Voigt strain mask (True = relax, False = fixed) that
    relaxes every direction NOT flagged as vacuum, and fixes the rest -- the
    normal strain along each vacuum axis, plus any shear that mixes a vacuum
    axis with another one. This generalizes cell relaxation to 2D (1 vacuum
    axis, e.g. a slab: in-plane lattice relaxes, vacuum thickness doesn't)
    and 1D (2 vacuum axes, e.g. a wire: only the periodic axis relaxes)
    structures, not just bulk 3D ones.

    Verified live against a graphene slab (vacuum along c): masked relaxation
    from a deliberately-off in-plane lattice constant converged to the real
    graphene value (~2.46-2.50 Ang) while the vacuum axis length changed by
    exactly 0.
    """
    mask = [not v for v in vacuum_axes]
    for ax_a, ax_b in _VOIGT_SHEAR_AXES.values():
        mask.append(not (vacuum_axes[ax_a] or vacuum_axes[ax_b]))
    return mask


def get_calculator(model="small", device="cpu", dtype="float64"):
    """Loads the MACE-MP-0 foundation potential as an ASE calculator.
    Default float64 for geometry optimization (MACE's own guidance is
    unambiguous that float32 is for MD, not geometry optimization) --
    callers doing MD (e.g. stb-amorphize's melt/quench stages) should pass
    dtype="float32" explicitly.
    """
    from mace.calculators import mace_mp
    return mace_mp(model=model, device=device, default_dtype=dtype)


def relax(atoms, calc, cell_mask=None, optimizer="FIRE", fmax=0.05, max_steps=200):
    """Relaxes `atoms` in place (positions, plus cell if `cell_mask` is
    given -- a 6-element Voigt mask from build_cell_mask()). Returns
    (converged, steps_used).
    """
    from ase.optimize import FIRE, BFGS, LBFGS
    from ase.filters import FrechetCellFilter
    optimizers = {"FIRE": FIRE, "BFGS": BFGS, "LBFGS": LBFGS}

    atoms.calc = calc
    target = FrechetCellFilter(atoms, mask=cell_mask) if cell_mask is not None else atoms
    opt = optimizers[optimizer](target, logfile=None)
    converged = opt.run(fmax=fmax, steps=max_steps)
    return converged, opt.nsteps
