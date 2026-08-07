"""Shared MACE-MP-0 relaxation helpers, used by stb-mlrelax, stb-defect's
--ml-rank, stb-amorphize, stb-adsorb's --ml-rank/--ml-prerelax, and
stb-neb's --ml-neb/--ml-prerelax-endpoints (relax_neb() below is the only
function of the four that operates on a whole band of images at once
instead of a single Atoms object -- everything else here, including
--ml-prerelax-endpoints' single-endpoint relax, reuses relax()). Callers
must call core.deps.require_mace() themselves first -- this module assumes
mace/torch are already available, and only imports them lazily inside each
function (not at module level), so merely importing this module never
forces the heavy PyTorch/mace dependency chain to load.
"""

import numpy as np

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


def gpu_available():
    """Returns (available, detail) -- whether a CUDA-capable GPU is actually
    usable by the installed PyTorch build, and a short human-readable reason
    either way (e.g. "NVIDIA GeForce RTX 3080" when available, or "PyTorch
    was not built with CUDA support (a CPU-only build)" when not). Imports
    torch lazily (same policy as the rest of this module -- callers that
    never touch GPU/MACE at all shouldn't pay for the import). Used both to
    validate an explicit --device cuda request (resolve_device below) and,
    standalone, by the interactive stb-suite menu to tell the user up front
    whether choosing "cuda" will actually work on this machine, before they
    commit to it.
    """
    import torch
    if not torch.backends.cuda.is_built():
        return False, "PyTorch was not built with CUDA support (a CPU-only build)"
    if not torch.cuda.is_available():
        return False, ("no CUDA-capable GPU was detected by the installed PyTorch build -- "
                        "check 'nvidia-smi' and your GPU drivers")
    return True, torch.cuda.get_device_name(0)


def resolve_device(device):
    """Validates `device` ('cpu' or 'cuda') against what's actually
    available, raising a clear ValueError -- not a silent fallback to CPU,
    and not letting torch's own often-cryptic runtime error surface deep
    inside a model-loading call instead -- if 'cuda' is requested but
    gpu_available() says no. A silent fallback would quietly give a much
    slower run than the user asked for, with no indication anything was
    wrong. Returns `device` unchanged when valid ('cpu' is always valid,
    trivially).
    """
    if device != "cuda":
        return device
    available, detail = gpu_available()
    if not available:
        raise ValueError(
            f"--device cuda requested, but {detail}. Install a CUDA-enabled PyTorch build "
            "(see https://pytorch.org/get-started/locally/) or pass --device cpu instead."
        )
    return device


def get_calculator(model="small", device="cpu", dtype="float64"):
    """Loads a MACE potential as an ASE calculator: the MACE-MP-0 foundation
    model (`model` = "small"/"medium"/"large", downloaded/cached on first
    use) by default, or a custom model file (`model` = a path to a `.model`
    -- e.g. one fine-tuned by stb-mlffAnalysis on real SIESTA data) when
    `model` resolves to an existing file. This single entry point is shared
    by every ML consumer in the suite (stb-mlrelax, stb-defect --ml-rank,
    stb-amorphize, stb-adsorb/stb-neb --ml-*), so a custom model becomes
    usable everywhere this is called from -- stb-mlrelax, stb-mlmd,
    stb-mlphonons, stb-mlelastic, and stb-mlsearch already expose a
    --custom-model CLI flag for it, the others can gain one the same way
    once there's a real need (same "expose on first genuine use" policy as
    everywhere else in this suite).

    `device` is validated via resolve_device() before anything is loaded --
    every caller in the suite funnels through this one function, so this is
    the single place a bad --device cuda request (no GPU/no CUDA-enabled
    torch) is caught, consistently, everywhere.

    Default float64 for geometry optimization (MACE's own guidance is
    unambiguous that float32 is for MD, not geometry optimization) --
    callers doing MD (e.g. stb-amorphize's melt/quench stages) should pass
    dtype="float32" explicitly.
    """
    device = resolve_device(device)
    import os
    if os.path.isfile(model):
        from mace.calculators import MACECalculator
        return MACECalculator(model_paths=[model], device=device, default_dtype=dtype)
    from mace.calculators import mace_mp
    return mace_mp(model=model, device=device, default_dtype=dtype)


def describe_model(model_arg, calc):
    """Best-effort extra detail about a loaded MACE model (approximate
    parameter count, cutoff radius) -- purely informational, so any failure
    to introspect internal attributes (these aren't part of MACE's public
    API and can differ across mace-torch versions/calculator classes) is
    swallowed silently rather than breaking the run. Moved here from
    mlrelax.py once stacking2D.py's --ml-relax became a second consumer of
    the exact same introspection.
    """
    lines = []
    try:
        models = getattr(calc, "models", None) or [getattr(calc, "model")]
        n_params = sum(p.numel() for p in models[0].parameters())
        lines.append(f"MACE model parameters : ~{n_params:,}")
    except Exception:
        pass
    try:
        r_max = float(models[0].r_max)
        lines.append(f"MACE cutoff radius     : {r_max:.2f} Ang")
    except Exception:
        pass
    return lines


def relax(atoms, calc, cell_mask=None, optimizer="FIRE", fmax=0.05, max_steps=200,
          step_history=None):
    """Relaxes `atoms` in place (positions, plus cell if `cell_mask` is
    given -- a 6-element Voigt mask from build_cell_mask()). Returns
    (converged, steps_used).

    If `step_history` is given (a list), appends one (step_index, energy,
    max|F|) tuple per optimizer iteration via ASE's standard
    Optimizer.attach() observer hook -- lets a caller (stb-mlrelax) plot/
    report the optimizer's own convergence trace without this function
    reimplementing the step loop itself. None (default) is a no-op, so
    every other caller of this shared function (stb-defect --ml-rank,
    stb-amorphize, stb-adsorb, stb-neb, stb-mlsearch, stb-mlconvergence)
    is unaffected.
    """
    from ase.optimize import FIRE, BFGS, LBFGS
    from ase.filters import FrechetCellFilter
    optimizers = {"FIRE": FIRE, "BFGS": BFGS, "LBFGS": LBFGS}

    atoms.calc = calc
    target = FrechetCellFilter(atoms, mask=cell_mask) if cell_mask is not None else atoms
    opt = optimizers[optimizer](target, logfile=None)
    if step_history is not None:
        def _record():
            step_history.append((len(step_history), atoms.get_potential_energy(),
                                  float(np.abs(atoms.get_forces()).max())))
        opt.attach(_record, interval=1)
    converged = opt.run(fmax=fmax, steps=max_steps)
    return converged, opt.nsteps


def relax_neb(images, calc, k=0.1, fmax=0.05, max_steps=200, optimizer="FIRE", climb=True):
    """Runs a climbing-image NEB in place on `images` (a list of ASE Atoms
    -- both endpoints plus already-guessed interior images, e.g. from
    linear or IDPP interpolation). All images must share one exact cell --
    ase.mep.neb.NEB raises NotImplementedError otherwise (no variable-cell
    NEB support in ASE). One `calc` is shared across every image
    (allow_shared_calculator=True -- a MACE ASE calculator holds no
    cross-call mutable state that a second concurrent Atoms object sharing
    it would corrupt; building one calculator per image would just reload
    the same model N times for no benefit).

    Two-stage strategy, standard climbing-image NEB practice: running with
    climb=True from step 0 lets a still badly-shaped band pick the wrong
    interior image as the climbing image (whichever happens to be
    highest-energy before the band has relaxed into its true
    saddle-adjacent shape), which can then converge onto a spurious
    stationary point instead of the real transition state. So stage 1
    (climb=False, a loosened fmax, half the step budget) lets the band
    find its overall shape first; stage 2 sets neb.climb = True (a plain
    mutable attribute on the NEB object) and tightens to the caller's real
    fmax for the remaining step budget, refining only the climbing image
    onto the true saddle while the rest of the band keeps relaxing
    normally.

    Endpoints (images[0]/images[-1]) are never moved by either stage
    (ase.mep.neb.BaseNEB.get_positions/set_positions only touch
    images[1:-1]) -- they're taken as given, the same "already relaxed,
    not this tool's job to move" convention as every other WORKFLOW_TOOLS
    pair in this suite.

    Returns (converged, steps_stage1, steps_stage2, energies) where
    energies is each image's potential energy (eV) after the final stage,
    in image order.
    """
    from ase.optimize import FIRE, BFGS, LBFGS
    from ase.mep.neb import NEB
    optimizers = {"FIRE": FIRE, "BFGS": BFGS, "LBFGS": LBFGS}

    for image in images:
        image.calc = calc

    # method="improvedtangent" pinned explicitly only to silence ASE's
    # UserWarning about the aseneb -> improvedtangent default change in
    # recent versions -- it's already the new default, so this changes no
    # behavior, just avoids spurious warning noise in the persistent report.
    neb = NEB(images, k=k, climb=False, allow_shared_calculator=True,
              method="improvedtangent")
    opt1 = optimizers[optimizer](neb, logfile=None)
    stage1_fmax = max(2 * fmax, 0.1)
    stage1_steps = max(max_steps // 2, 1)
    opt1.run(fmax=stage1_fmax, steps=stage1_steps)

    neb.climb = climb
    opt2 = optimizers[optimizer](neb, logfile=None)
    stage2_steps = max(max_steps - opt1.nsteps, 1)
    converged = opt2.run(fmax=fmax, steps=stage2_steps)

    energies = [image.get_potential_energy() for image in images]
    return converged, opt1.nsteps, opt2.nsteps, energies
