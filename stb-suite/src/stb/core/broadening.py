"""Shared Gaussian-broadening convention: sigma/FWHM given in meV, converted
to eV internally. Extracted from convdos.py once stb-sts/stb-coop became
further consumers of the same --sigma/--fwhm mutual-exclusivity and unit
convention -- one source of truth for the conversion constant and the
validation logic, instead of three (or four) independently re-derived
copies.
"""

import numpy as np

# FWHM = 2*sqrt(2*ln2) * sigma, for a Gaussian.
FWHM_TO_SIGMA = 2 * np.sqrt(2 * np.log(2))


def sigma_ev_from_args(sigma_mev, fwhm_mev):
    """Converts a --sigma/--fwhm pair (exactly one of which is expected to
    be non-None, enforced by the caller's argparse mutually-exclusive
    group) into a sigma value in eV. Raises ValueError if the given value
    isn't positive."""
    if sigma_mev is not None:
        if sigma_mev <= 0:
            raise ValueError("--sigma must be positive.")
        sigma_mev_value = sigma_mev
    else:
        if fwhm_mev <= 0:
            raise ValueError("--fwhm must be positive.")
        sigma_mev_value = fwhm_mev / FWHM_TO_SIGMA
    return sigma_mev_value / 1000.0
