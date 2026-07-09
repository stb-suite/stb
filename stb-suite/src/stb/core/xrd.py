"""Shared powder-XRD helpers, built on pyxtal's diffraction engine (not
pymatgen's -- pyxtal is already a hard dependency for stb-crystalcast, and
its Similarity class and March-Dollase preferred-orientation support have
no pymatgen equivalent).

Extracted from stb-xrd once stb-xrdrank needed the same pattern-computation
and experimental-file-reading logic, per this repo's own extract-on-second
-use policy (see structure_io.py's module docstring for the same rationale
applied to format readers).
"""

import numpy as np
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.diffraction.xrd import WAVELENGTHS

MIN_EXPERIMENTAL_POINTS = 4  # pyxtal.XRD.Similarity cubic-interpolates each pattern (scipy
                             # interp1d, kind='cubic'), which needs at least 4 points per curve
                             # -- fewer than that raises a raw, unhandled scipy ValueError deep
                             # inside Similarity() instead of a clean error.


def resolve_wavelength(spec):
    """Returns (wavelength_in_ang, label): `spec` matched case-insensitively
    against pymatgen's WAVELENGTHS dict of named X-ray sources (CuKa, MoKa,
    ...), or else parsed directly as a wavelength in Ang. Raises ValueError
    with a clear message (listing the known names) if neither works.
    """
    for name, value in WAVELENGTHS.items():
        if name.lower() == spec.lower():
            return value, name
    try:
        value = float(spec)
    except ValueError:
        raise ValueError(
            f"--wavelength '{spec}' is not a known source name and isn't a number. "
            f"Known names: {', '.join(WAVELENGTHS)}.")
    return value, f"{value:.5f} Ang"


def read_experimental_pattern(path):
    """Reads a plain 2-column (2theta, intensity) text file: whitespace- or
    comma-separated, blank lines and '#' comments skipped. Returns a (2, N)
    array -- the same shape pyxtal.XRD.Similarity expects for both patterns
    it compares (it interpolates internally, so this doesn't need to share
    a grid with a simulated pattern).
    """
    two_theta, intensity = [], []
    with open(path) as f:
        for lineno, raw_line in enumerate(f, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.replace(",", " ").split()
            if len(parts) < 2:
                continue
            try:
                two_theta.append(float(parts[0]))
                intensity.append(float(parts[1]))
            except ValueError:
                raise ValueError(f"'{path}' line {lineno}: expected two numbers, got '{line}'.")
    if len(two_theta) < MIN_EXPERIMENTAL_POINTS:
        raise ValueError(
            f"'{path}' has only {len(two_theta)} data point(s) -- at least "
            f"{MIN_EXPERIMENTAL_POINTS} are needed for the similarity comparison.")
    return np.array([two_theta, intensity])


def compute_pattern(structure, wavelength=1.54184, two_theta_range=(0.0, 90.0)):
    """Builds a pyxtal.XRD.XRD powder pattern for a pymatgen Structure.

    Raises ValueError with a clean message -- instead of letting pyxtal's own
    bare ValueError (it calls max() on an empty list internally) leak out --
    when no reflection falls inside `two_theta_range`.
    """
    from pyxtal.XRD import XRD

    atoms = AseAtomsAdaptor.get_atoms(structure)
    lo, hi = two_theta_range
    try:
        xrd = XRD(atoms, wavelength=wavelength, thetas=[lo, hi])
    except ValueError:
        raise ValueError("no peaks found in the given two-theta range -- try widening it.")
    if len(xrd.pxrd) == 0:
        raise ValueError("no peaks found in the given two-theta range -- try widening it.")
    return xrd
