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


def derive_from_eps(omega_ev, eps1, eps2):
    """n/k/alpha/R/L/sigma1 from an ALREADY-KNOWN eps1/eps2 pair -- the
    second half of compute_all, factored out so a caller that already
    has a self-consistent (eps1, eps2) pair from somewhere OTHER than a
    fresh Kramers-Kronig transform (e.g. the dimensionality-correction
    functions below, whose eps1_2D/eps2_2D come out of the correction
    formula already paired) doesn't have to -- and must NOT -- re-run KK
    on it. Returns the same dict shape as compute_all (minus omega/eps1/
    eps2, which the caller already has).
    """
    omega_ev = np.asarray(omega_ev, dtype=float)
    eps1 = np.asarray(eps1, dtype=float)
    eps2 = np.asarray(eps2, dtype=float)
    n, k = refractive_index(eps1, eps2)
    alpha = absorption_coefficient(k, omega_ev)
    R = reflectivity(n, k)
    L = energy_loss_function(eps1, eps2)
    sigma1 = optical_conductivity(eps2, omega_ev)
    return {
        "omega": omega_ev, "eps1": eps1, "eps2": eps2,
        "n": n, "k": k, "alpha": alpha, "R": R, "L": L, "sigma1": sigma1,
    }


def compute_all(omega_ev, eps2):
    """Single entry point: given one direction's raw (omega_ev, eps2)
    read from a SystemLabel.EPSIMG (core.dielectric.read_epsimg), runs
    the full Kramers-Kronig transform (core.dielectric.kramers_kronig,
    reused as-is -- not reimplemented) to get eps1(E), then every
    derived quantity via derive_from_eps, vectorized over the whole
    energy grid.

    Returns a dict: omega, eps1, eps2, n, k, alpha, R, L, sigma1 (each a
    1D array of the same length as omega_ev).
    """
    omega_ev = np.asarray(omega_ev, dtype=float)
    eps2 = np.asarray(eps2, dtype=float)
    eps1 = kramers_kronig(omega_ev, eps2)
    return derive_from_eps(omega_ev, eps1, eps2)


def correct_2d_perpendicular(eps1_sc, eps2_sc, cell_length_ang, thickness_ang):
    """Restores the intrinsic 2D out-of-plane (perpendicular to the
    layer) dielectric function from a vacuum-padded supercell's raw
    computed eps_SC, via the displacement-field (D) continuity argument:

        eps_2D,perp(E) = [1 + (c/t) * (1/eps_SC,perp(E) - 1)]^-1

    c = cell_length_ang (the supercell's own length along the vacuum
    axis), t = thickness_ang (the 2D material's physical thickness,
    user-supplied -- see stb-opticalAnalysis --help for why this isn't
    auto-estimated). eps_SC = eps1_sc + i*eps2_sc is complex, so this is
    computed with numpy complex arithmetic throughout, not separate
    real/imaginary algebra.

    Reference: Yang & Gao, npj 2D Mater. Appl. 5, 78 (2021), "A method
    to restore the intrinsic dielectric functions of 2D materials in
    periodic calculations..."; following Laturia, Van de Put,
    Vandenberghe, npj 2D Mater. Appl. 2, 6 (2018), "Dielectric
    properties of hexagonal boron nitride and transition metal
    dichalcogenides: from monolayer to bulk." Formula confirmed verbatim
    against a secondary source (PMC8376903) quoting Laturia et al.'s own
    Eq. 20 (2018) during this tool's development.

    Returns (eps1_2d, eps2_2d) -- the real/imaginary parts of the
    corrected complex dielectric function, same array shape as input.
    """
    eps1_sc = np.asarray(eps1_sc, dtype=float)
    eps2_sc = np.asarray(eps2_sc, dtype=float)
    eps_sc = eps1_sc + 1j * eps2_sc
    ratio = cell_length_ang / thickness_ang
    eps_2d = 1.0 / (1.0 + ratio * (1.0 / eps_sc - 1.0))
    return eps_2d.real, eps_2d.imag


def correct_2d_parallel(eps1_sc, eps2_sc, cell_length_ang, thickness_ang):
    """Restores the intrinsic 2D in-plane dielectric function via the
    electric-field (E) continuity argument (linear rescaling, unlike
    the perpendicular case's complex inversion):

        eps_2D,para(E) = 1 + (c/t) * (eps_SC,para(E) - 1)

    Same variables/reference as correct_2d_perpendicular. Returns
    (eps1_2d, eps2_2d).
    """
    eps1_sc = np.asarray(eps1_sc, dtype=float)
    eps2_sc = np.asarray(eps2_sc, dtype=float)
    eps_sc = eps1_sc + 1j * eps2_sc
    ratio = cell_length_ang / thickness_ang
    eps_2d = 1.0 + ratio * (eps_sc - 1.0)
    return eps_2d.real, eps_2d.imag


_ANG3_PER_M3 = 1.0e30


def molecular_polarizability(eps1_sc, eps2_sc, cell_volume_ang3):
    """Molecular polarizability alpha(E) of an isolated (0D) system from
    its own vacuum-padded supercell's raw eps_SC, via the dilute-limit
    Clausius-Mossotti relation (eps_SC -> 1 as the box volume -> large,
    the normal regime for a single molecule in a big vacuum box):

        alpha(E) [SI, C*m^2/V] = (eps_SC(E) - 1) * eps0 * V_box

    Also returned in cubic-Angstrom (the conventional unit for comparing
    against tabulated CCSD/experimental molecular polarizabilities):
    alpha_Ang3 = alpha_SI / (4*pi*eps0) * 1e30 (the standard SI<->Gaussian
    polarizability-volume conversion, alpha_SI = 4*pi*eps0*alpha_Gaussian).

    Standard classical electrostatics result (dilute limit of the
    Clausius-Mossotti/Lorentz-Lorenz relation) -- well-established, but
    the exact convention (this dilute-limit form vs. the full Clausius-
    Mossotti relation) has not been separately cross-checked against a
    computational-chemistry-specific primary source during this tool's
    development; flagged here for awareness, not as a correctness doubt.

    Returns (alpha1_SI, alpha2_SI, alpha1_ang3, alpha2_ang3) -- real and
    imaginary parts in both unit conventions.
    """
    eps1_sc = np.asarray(eps1_sc, dtype=float)
    eps2_sc = np.asarray(eps2_sc, dtype=float)
    eps_sc = eps1_sc + 1j * eps2_sc
    alpha_si = (eps_sc - 1.0) * constants.epsilon_0 * cell_volume_ang3 * 1e-30
    alpha_ang3 = alpha_si / (4.0 * np.pi * constants.epsilon_0) * _ANG3_PER_M3
    return alpha_si.real, alpha_si.imag, alpha_ang3.real, alpha_ang3.imag
