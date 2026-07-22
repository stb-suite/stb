"""Shared equation-of-state (E vs V) fitting helpers, used by stb-mleos (MACE
-driven) and stb-eosAnalysis (real SIESTA data). Deliberately its own module
rather than living in mleos.py: mleos.py calls core.deps.require_mace() at
import time, so anything importing it -- including a plain DFT-workflow tool
that never touches MACE/torch at all -- would be forced to load the heavy
optional 'ml' dependency chain just to fit a curve. This module only needs
ase.eos (a core, always-installed dependency), so both consumers can import
it unconditionally.
"""

import numpy as np

# eV/Ang^3 -> GPa. Same value as elastic_analysis.py's own CONV_EVA3_TO_GPA,
# kept as a local copy rather than imported from there -- core/ modules are
# the shared foundation tools import FROM, not the reverse, and this is a
# single fixed physical conversion factor (1 eV/Ang^3 = 160.21766 GPa),
# not something that could drift between the two copies.
CONV_EVA3_TO_GPA = 160.21766

# CLI-friendly aliases -> ase.eos.EquationOfState's own eos_string spelling.
# 'birch_murnaghan' accepted as an alias since that's the spelling
# stb-mlphonons's --qha --eos flag already uses (phonopy's own convention);
# ASE spells the same EOS without the underscore.
_EOS_ALIASES = {"birch_murnaghan": "birchmurnaghan"}


def normalize_eos_string(name):
    return _EOS_ALIASES.get(name, name)


def fit_eos(volumes, energies, eos_string):
    """Fits ase.eos.EquationOfState and returns a dict of the derived
    quantities plus the fitted curve for plotting. B0/B0' are None for EOS
    forms that don't expose a 4-parameter [E0, B0, B0', V0] fit (only 'sj'
    among the ones offered by either caller -- its own fit_sjeos() doesn't
    populate eos_parameters the same way).

    Calls eos.fit(warn=False) rather than its own default (warn=True) --
    ase's own warning ("The minimum volume of your fit is not in your
    volumes") goes through Python's warnings module, which bypasses every
    caller's own print_dual-formatted report entirely (easy to miss, or to
    print in the wrong place/format). v0_outside_range in the returned dict
    is the same check done explicitly instead, so callers can raise their
    own clean [WARNING] in-report.
    """
    from ase.eos import EquationOfState

    eos = EquationOfState(volumes, energies, eos=eos_string)
    v0, e0, B = eos.fit(warn=False)
    b0_gpa = B * CONV_EVA3_TO_GPA

    bprime = None
    if eos_string != "sj":
        params = np.asarray(eos.eos_parameters)
        if params.size == 4:
            bprime = params[2]

    _, _, _, _, x, y, v, e = eos.getplotdata()

    # Fit-quality R^2 against the actual scanned data points (not the fine
    # plotting grid) -- a real BM/Vinet/etc. curve, not just a passthrough
    # of curve_fit's own internal residual.
    y_at_data = eos.func(v, *eos.eos_parameters) if eos_string != "sj" else eos.fit0(v ** -(1.0 / 3.0))
    ss_res = np.sum((e - y_at_data) ** 2)
    ss_tot = np.sum((e - np.mean(e)) ** 2)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0

    v0_outside_range = not (float(np.min(v)) <= v0 <= float(np.max(v)))

    return {
        "v0": v0, "e0": e0, "b0_gpa": b0_gpa, "bprime": bprime,
        "r_squared": r_squared, "curve_v": x, "curve_e": y,
        "v0_outside_range": v0_outside_range,
    }


def invert_pressure(fit, target_pressure_gpa):
    """Given a fit_eos() result, finds the volume at which the fitted EOS
    predicts external pressure P = target_pressure_gpa (from the standard
    relation P = -dE/dV). Uses a numerical derivative (np.gradient) of the
    already-fitted curve (fit['curve_v']/['curve_e'], a smooth analytic BM/
    Vinet/etc. form evaluated on a fine 100-point grid) rather than an
    analytic per-EOS-form derivative -- one implementation works for every
    EOS string this module supports, at the cost of a small, negligible
    finite-difference error on an already-smooth curve.

    Returns (volume, extrapolated): `extrapolated` is True if
    target_pressure_gpa itself falls outside the range of pressures the
    fitted curve actually spans over its scanned volumes (fit['curve_v']'s
    own range -- exactly the volumes actually scanned). Checking the
    REQUESTED PRESSURE against that range, not the resulting volume against
    fit['curve_v'], matters: np.interp clamps to the nearest endpoint
    outside its xp range rather than extrapolating, so a wildly out-of
    -range pressure would otherwise silently return an in-range (boundary)
    volume and look like a normal interpolation.
    """
    curve_v = fit["curve_v"]
    curve_e = fit["curve_e"]
    dEdV = np.gradient(curve_e, curve_v)  # eV/Ang^3
    pressure_gpa = -dEdV * CONV_EVA3_TO_GPA

    order = np.argsort(pressure_gpa)
    pressure_sorted = pressure_gpa[order]
    volume = float(np.interp(target_pressure_gpa, pressure_sorted, curve_v[order]))
    extrapolated = not (float(np.min(pressure_sorted)) <= target_pressure_gpa <= float(np.max(pressure_sorted)))
    return volume, extrapolated
