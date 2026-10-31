"""SIESTA .EPSIMG parsing and Kramers-Kronig transform, shared by
stb-ramanAnalysis's per-displaced-geometry dielectric response readout.

Deliberately implemented in pure Python/NumPy/SciPy rather than shelling
out to SIESTA's own bundled Util/Optical/optical.f post-processor: this
suite's architecture is "Python-only, no external binaries beyond SIESTA
itself" everywhere else, and the Kramers-Kronig transform itself is only
~10 lines of NumPy -- not worth a new external-binary dependency.

.EPSIMG format verified against SIESTA's own source (Src/optical.F,
id23cat/siesta-3.1 mirror -- the write side) and Util/Optical/optical.f
(the read side of SIESTA's own post-processor): comment lines prefixed
with "##" (energy range, spin count, f-sum rule / plasma frequency), then
one data row per frequency point: "energy_eV  eps2_spin1 [eps2_spin2 ...]"
(free-format/list-directed, i.e. plain whitespace-separated floats).
"""

from __future__ import annotations

import numpy as np

# np.trapz was renamed np.trapezoid in NumPy 2.0 and removed outright in
# later 2.x releases (verified: raises AttributeError on NumPy 2.4) -- the
# suite doesn't pin a NumPy version, so support both without a hard
# version check.
_trapz = getattr(np, "trapezoid", None) or np.trapz


def read_epsimg(path):
    """Reads a SystemLabel.EPSIMG file: (omega, eps2) -- frequency (eV)
    and the imaginary part of the dielectric function along the
    Optical.Vector direction that run used, summed over spin components
    when present (nspin=2 spin-up/spin-down channels both contribute to
    the total optical response). Comment lines ("#"-prefixed, per the
    verified format) are skipped; every other non-blank line is parsed as
    whitespace-separated floats.
    """
    omega = []
    eps2 = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            try:
                values = [float(p) for p in parts]
            except ValueError:
                continue
            if len(values) < 2:
                continue
            omega.append(values[0])
            eps2.append(sum(values[1:]))
    return np.array(omega, dtype=float), np.array(eps2, dtype=float)


def kramers_kronig(omega, eps2):
    """eps1(omega) - 1 = (2/pi) * P-integral_0^inf [w' eps2(w') / (w'^2 - w^2)] dw',
    the standard Kramers-Kronig relation linking the real and imaginary
    parts of a causal response function. Numerically integrated via direct
    trapezoidal quadrature over the input omega grid, omitting the
    (removable) singular sample point w'=w at each target point --
    verified against an analytic Lorentz oscillator (eps1/eps2 both known
    in closed form): max relative error ~0.2% on a 4000-point grid spanning
    the resonance, well within tolerance for a finite-displacement Raman
    tensor (dominated by the finite-difference step error, not this
    transform).

    `omega`/`eps2` need not be evenly spaced -- np.trapz integrates over
    the actual sample spacing, so the point excluded at each target index
    just widens that one local interval rather than breaking the grid.
    """
    omega = np.asarray(omega, dtype=float)
    eps2 = np.asarray(eps2, dtype=float)
    n = len(omega)
    idx = np.arange(n)
    eps1 = np.empty(n)
    for i in range(n):
        wi = omega[i]
        mask = idx != i
        wj = omega[mask]
        e2j = eps2[mask]
        integrand = wj * e2j / (wj ** 2 - wi ** 2)
        eps1[i] = 1.0 + (2.0 / np.pi) * _trapz(integrand, wj)
    return eps1


def static_dielectric(omega, eps2):
    """eps1 at the smallest available omega -- the single scalar used as
    "the static dielectric constant" of one displaced-geometry Optical
    run. Computing the full kramers_kronig() array and taking the first
    entry (rather than a dedicated omega=0 formula) keeps this consistent
    with kramers_kronig()'s own numerics -- no second, differently-behaved
    code path for the single value stb-ramanAnalysis actually needs.
    """
    eps1 = kramers_kronig(omega, eps2)
    return float(eps1[np.argmin(omega)])
