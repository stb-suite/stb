"""Linear optical properties derived from the complex dielectric function
eps(E) = eps1(E) + i*eps2(E), shared by stb-opticalAnalysis.

Physics: standard textbook relations linking a material's complex
dielectric function to its complex refractive index N(E) = n(E) + i*k(E)
(N^2 = eps), and from there to every other linear-response optical
observable computed here -- absorption coefficient, normal-incidence
reflectivity, the electron energy-loss function, and the real part of
the optical conductivity. See e.g. F. Wooten, "Optical Properties of
Solids" (Academic Press, 1972), or M. Fox, "Optical Properties of
Solids" (Oxford, 2010). Unlike this suite's GQCA module, none of these
formulas need an [UNVERIFIED] flag -- they are well-established,
non-contested relations, not a still-developing formalism.

Deliberately independent of any SIESTA-file-parsing concerns -- eps1/eps2
are passed in as plain arrays (the caller gets them from
core.dielectric.read_epsimg + core.dielectric.kramers_kronig), keeping
this module light and unit-testable without any file I/O, same
separation-of-concerns precedent as core/gqca_solver.py.
"""

from __future__ import annotations

import numpy as np
from scipy import constants

from stb.core.dielectric import kramers_kronig

# hc in eV*cm, derived from scipy.constants (not a hardcoded literal) --
# used by absorption_coefficient to convert a photon energy in eV to an
# absorption coefficient in the conventional cm^-1 units.
HC_EV_CM = constants.h * constants.c / constants.e * 100.0


def refractive_index(eps1, eps2):
    """Complex refractive index N(E) = n(E) + i*k(E) from N^2 = eps(E):

        n(E) = sqrt( (sqrt(eps1^2 + eps2^2) + eps1) / 2 )
        k(E) = sqrt( (sqrt(eps1^2 + eps2^2) - eps1) / 2 )

    k >= 0 always (it's a sqrt of a non-negative quantity, since
    sqrt(eps1^2+eps2^2) >= |eps1| unconditionally) -- the sign convention
    that makes N=n+ik describe an absorptive (not amplifying) medium.
    """
    eps1 = np.asarray(eps1, dtype=float)
    eps2 = np.asarray(eps2, dtype=float)
    eps_mag = np.sqrt(eps1 ** 2 + eps2 ** 2)
    n = np.sqrt(np.clip((eps_mag + eps1) / 2.0, 0.0, None))
    k = np.sqrt(np.clip((eps_mag - eps1) / 2.0, 0.0, None))
    return n, k


def absorption_coefficient(k, energy_ev):
    """alpha(E) [cm^-1] = 4*pi*k(E)*E[eV] / hc[eV*cm] -- the standard
    Beer-Lambert absorption coefficient (intensity ~ exp(-alpha*x)),
    from E = hc/lambda and alpha = 4*pi*k/lambda.
    """
    k = np.asarray(k, dtype=float)
    energy_ev = np.asarray(energy_ev, dtype=float)
    return 4.0 * np.pi * k * energy_ev / HC_EV_CM


def reflectivity(n, k):
    """R(E) = ((n-1)^2 + k^2) / ((n+1)^2 + k^2) -- normal-incidence
    Fresnel reflectance at a vacuum(n0=1)/material interface. This
    assumption (vacuum incidence, normal angle) is NOT generalized here
    -- angle-dependent/thin-film-interference reflectivity is out of
    scope for this module.
    """
    n = np.asarray(n, dtype=float)
    k = np.asarray(k, dtype=float)
    return ((n - 1.0) ** 2 + k ** 2) / ((n + 1.0) ** 2 + k ** 2)


def energy_loss_function(eps1, eps2):
    """L(E) = eps2 / (eps1^2 + eps2^2) = Im(-1/eps(E)) -- the electron
    energy-loss function (EELS), whose peaks locate plasmon-like
    excitations (where eps1 crosses zero with small eps2).
    """
    eps1 = np.asarray(eps1, dtype=float)
    eps2 = np.asarray(eps2, dtype=float)
    return eps2 / (eps1 ** 2 + eps2 ** 2)


def optical_conductivity(eps2, energy_ev):
    """sigma1(E) [S/m, SI] = eps0 * omega(E) * eps2(E), the real part of
    the optical conductivity -- derived from sigma(omega) =
    -i*omega*eps0*(eps_r(omega) - 1); its real part is exactly
    omega*eps0*eps2 (the -1 and the imaginary unit only affect the
    imaginary part, sigma2, which is out of scope here).
    omega(E) = E[eV]*e/hbar (angular frequency, rad/s).
    """
    eps2 = np.asarray(eps2, dtype=float)
    energy_ev = np.asarray(energy_ev, dtype=float)
    omega = energy_ev * constants.e / constants.hbar
    return constants.epsilon_0 * omega * eps2


def compute_all(omega_ev, eps2):
    """Single entry point: given one direction's raw (omega_ev, eps2)
    read from a SystemLabel.EPSIMG (core.dielectric.read_epsimg), runs
    the full Kramers-Kronig transform (core.dielectric.kramers_kronig,
    reused as-is -- not reimplemented) to get eps1(E), then every
    derived quantity above, vectorized over the whole energy grid.

    Returns a dict: omega, eps1, eps2, n, k, alpha, R, L, sigma1 (each a
    1D array of the same length as omega_ev).
    """
    omega_ev = np.asarray(omega_ev, dtype=float)
    eps2 = np.asarray(eps2, dtype=float)
    eps1 = kramers_kronig(omega_ev, eps2)
    n, k = refractive_index(eps1, eps2)
    alpha = absorption_coefficient(k, omega_ev)
    R = reflectivity(n, k)
    L = energy_loss_function(eps1, eps2)
    sigma1 = optical_conductivity(eps2, omega_ev)
    return {
        "omega": omega_ev, "eps1": eps1, "eps2": eps2,
        "n": n, "k": k, "alpha": alpha, "R": R, "L": L, "sigma1": sigma1,
    }
