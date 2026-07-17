"""Physically-generic spectrum-building helpers shared by the Raman and
IR spectrum workflows -- moved out of raman_analysis.py once
ir_analysis.py became a second consumer of code that never had any
Raman-specific content to begin with (unit-agnostic Lorentzian summation,
generic Bose-Einstein phonon population physics, generic peak-finding/
matching against an experimental spectrum). Left behind in
raman_analysis.py (not moved here): write_spectrum_plot (its .dat/.gplot
text has real Raman-specific wording -- the diagonal/full-tensor scope
caveat has no IR equivalent) and raman_activity/depolarization_ratio
(Raman-tensor-specific physics).
"""

from __future__ import annotations

import numpy as np
from scipy.constants import h as PLANCK_H, k as BOLTZMANN_K
from scipy.signal import find_peaks


def bose_einstein_weight(freq_thz, temperature_k):
    """Stokes-scattering thermal population factor (n+1), where
    n = 1/(exp(hf/kT) - 1) is the Bose-Einstein occupation of the mode at
    this frequency and temperature (f = freq_thz * 1e12 Hz). At T -> 0,
    n -> 0 and the weight -> 1 (spontaneous Stokes scattering is always
    present, even with zero thermally excited phonons); at high T, the
    weight grows (~linearly in the classical limit), matching the
    stronger low-frequency Stokes lines seen experimentally at higher
    temperature. Same scipy.constants h/k pair phonons_pos.py already
    uses for its own Debye-temperature calculation.

    np.expm1 (exp(x)-1, computed directly rather than via exp(x) then
    subtracting 1) keeps this numerically stable for small x (high T or
    low frequency), where a naive exp(x)-1 loses precision to
    cancellation.

    Raises ValueError for freq_thz <= 0: Bose-Einstein occupation is only
    defined for a real (positive-frequency) vibration. An imaginary mode
    (dynamically unstable structure -- Stage 2 already flags these with
    [IMAGINARY]) isn't a real oscillator at all, and naively plugging a
    negative frequency into this formula doesn't degrade gracefully -- it
    blows up near x=0 (small |freq_thz|), verified live against this
    exact fabricated-imaginary-mode scenario during testing (intensities
    came out in the millions, wrong sign and all). Callers must filter
    these out before weighting, not rely on this function to do it softly.
    """
    if freq_thz <= 0:
        raise ValueError(
            f"Bose-Einstein weight is undefined for a non-positive frequency ({freq_thz} THz) "
            "-- this mode is imaginary/unstable, not a real vibration.")
    f_hz = freq_thz * 1e12
    x = (PLANCK_H * f_hz) / (BOLTZMANN_K * temperature_k)
    n = 1.0 / np.expm1(x)
    return n + 1.0


def build_lorentzian_spectrum(freqs_cm1, activities, linewidth_cm1, n_points=2000, pad_cm1=100.0):
    """Sum-of-Lorentzians spectrum over a frequency axis auto-ranged
    around the mode frequencies (padded by `pad_cm1` on each side).
    """
    if not freqs_cm1:
        return np.array([]), np.array([])
    fmin = max(0.0, min(freqs_cm1) - pad_cm1)
    fmax = max(freqs_cm1) + pad_cm1
    grid = np.linspace(fmin, fmax, n_points)
    intensity = np.zeros_like(grid)
    for f0, a in zip(freqs_cm1, activities):
        intensity += a * (linewidth_cm1 ** 2) / ((grid - f0) ** 2 + linewidth_cm1 ** 2)
    return grid, intensity


def read_experimental_spectrum(path):
    """Reads a 2-column (frequency cm^-1, intensity) experimental
    spectrum text file -- '#'-prefixed comment lines and blank lines
    skipped, same lenient whitespace-separated parsing convention as
    core.dielectric.read_epsimg. Returns (freq_cm1, intensity) arrays,
    sorted ascending by frequency (scipy.signal.find_peaks needs an
    ascending x-axis, and a raw export isn't guaranteed to already be
    sorted that way).
    """
    freq, intensity = [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                freq.append(float(parts[0]))
                intensity.append(float(parts[1]))
            except ValueError:
                continue
    freq = np.array(freq, dtype=float)
    intensity = np.array(intensity, dtype=float)
    order = np.argsort(freq)
    return freq[order], intensity[order]


def find_spectrum_peaks(freq, intensity, prominence_fraction=0.01):
    """Local maxima of `intensity` over `freq`, via scipy.signal.find_peaks.
    Works directly on a simulated Lorentzian-sum spectrum (smooth,
    well-separated peaks) and on a real experimental spectrum alike --
    `prominence_fraction` (relative to the spectrum's own peak intensity,
    default 1%) filters out trivial noise-level bumps without needing an
    absolute-intensity threshold, since simulated intensity/activity and
    real experimental intensity are on completely different, incomparable
    scales. Returns peak positions in `freq`'s own units (frequency
    array must already be ascending, same requirement as
    read_experimental_spectrum's output / build_lorentzian_spectrum's
    grid).
    """
    if len(intensity) == 0 or intensity.max() <= 0:
        return np.array([])
    peak_idx, _ = find_peaks(intensity, prominence=prominence_fraction * intensity.max())
    return freq[peak_idx]


def match_peaks(sim_peaks, exp_peaks):
    """Nearest-neighbor peak matching: for each experimental peak, the
    closest simulated peak in frequency (cm^-1). Simple by design --
    appropriate for sparse, well-separated vibrational spectra (Raman/IR),
    unlike a full XRD pattern (stb-xrdrank's normalized-cross-correlation
    machinery is solving a different, denser-pattern problem and isn't
    needed here). Returns a list of (exp_freq, sim_freq, delta) tuples,
    one per experimental peak, in the same order as `exp_peaks`;
    sim_freq/delta are None when there are no simulated peaks to match
    against at all.
    """
    matches = []
    for ef in exp_peaks:
        if len(sim_peaks) == 0:
            matches.append((ef, None, None))
            continue
        idx = int(np.argmin(np.abs(sim_peaks - ef)))
        sf = sim_peaks[idx]
        matches.append((ef, sf, sf - ef))
    return matches
