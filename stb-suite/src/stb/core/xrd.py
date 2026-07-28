"""Shared powder-XRD helpers, built on pyxtal's diffraction engine (not
pymatgen's -- pyxtal is already a hard dependency for stb-crystalcast, and
its Similarity class and March-Dollase preferred-orientation support have
no pymatgen equivalent).

Extracted from stb-xrd once stb-xrdrank needed the same pattern-computation
and experimental-file-reading logic, per this repo's own extract-on-second
-use policy (see structure_io.py's module docstring for the same rationale
applied to format readers).
"""

import os

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
    """Reads a 2-theta/intensity text file: whitespace- or comma-separated,
    blank lines and '#' comments skipped. Intensity is always taken from the
    LAST column, 2-theta from the first -- for a plain 2-column file that's
    the same thing, but it also makes stb-xrd's own `--save-gnuplot` output
    (`xrd_pattern.dat`, 6 columns: 2theta d h k l intensity) a valid
    `--compare-to` input directly, instead of silently reading its `d`
    column (Ang, index 1) as intensity. Returns a (2, N) array -- the same
    shape pyxtal.XRD.Similarity expects for both patterns it compares (it
    interpolates internally, so this doesn't need to share a grid with a
    simulated pattern).
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
                intensity.append(float(parts[-1]))
            except ValueError:
                raise ValueError(f"'{path}' line {lineno}: expected two numbers, got '{line}'.")
    if len(two_theta) < MIN_EXPERIMENTAL_POINTS:
        raise ValueError(
            f"'{path}' has only {len(two_theta)} data point(s) -- at least "
            f"{MIN_EXPERIMENTAL_POINTS} are needed for the similarity comparison.")
    return np.array([two_theta, intensity])


def looks_like_peak_list(two_theta, cv_threshold=0.3):
    """Heuristic: does `two_theta` look like a sparse list of discrete peak
    positions (e.g. an indexed/literature peak table, or stb-xrd's own
    stick-pattern `--save-gnuplot` output) rather than a densely,
    ~evenly-sampled continuous scan (a real diffractometer trace)?

    Real continuous scans have a near-constant step between consecutive
    (sorted) 2-theta points -- coefficient of variation (std/mean) of the
    spacing is close to 0. A peak list has a highly UNEVEN spacing instead:
    near-zero gaps between peaks that happen to sit close together, and
    large empty gaps everywhere else -- verified on stb-xrd's own 386-peak
    Sn3O4 pattern (CV ~1.7) vs. a synthetic uniform 0.02-deg-step scan over
    the same range (CV ~2e-13). `cv_threshold=0.3` sits well below the
    former and well above ordinary instrumental step-size jitter.
    """
    sorted_theta = np.sort(np.asarray(two_theta))
    diffs = np.diff(sorted_theta)
    if len(diffs) < 2 or diffs.mean() == 0:
        return True
    return (diffs.std() / diffs.mean()) > cv_threshold


def broaden_peak_list(two_theta, intensity, two_theta_range, res=0.01, fwhm=0.1):
    """Converts a discrete (2theta, intensity) peak list into a continuously
    sampled Gaussian-broadened profile, using the same defaults
    `pyxtal.XRD.XRD.get_profile()` applies to the simulated pattern
    (method='gaussian', res=0.01 deg, FWHM=0.1 deg) -- so a peak-list
    experimental input becomes directly comparable to the simulated profile
    instead of forcing pyxtal.XRD.Similarity to cubic-interpolate straight
    through the zero-intensity gaps between peaks (which is not the same
    shape as a real broadened line and drags the similarity score down even
    when both patterns come from the exact same structure -- verified live:
    comparing stb-xrd's own stick output for a structure against its own
    simulated pattern for that SAME structure scored only ~0.32-0.46
    unbroadened, vs. ~1.0 once both sides are broadened the same way).
    """
    from pyxtal.XRD import Profile

    lo, hi = two_theta_range
    profile = Profile(method="gaussian", res=res, user_kwargs={"FWHM": fwhm})
    return profile.get_profile(np.asarray(two_theta), np.asarray(intensity), lo, hi)


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


def write_xrd_data(dat_path, rows, total_peaks, structure_label, formula, space_group,
                    wavelength_label, wavelength, two_theta_range):
    """Writes the peak list (2theta, d, h, k, l, intensity) to a plain text
    file with a complete, human-readable header -- unlike pyxtal's own
    XRD.save() (a single terse "wavelength/thetas ..." comment line), this
    is meant to be the ONE place the full peak table lives: stb-xrd itself
    no longer prints the (often huge) peak table to the console/report, on
    the reasoning that it's redundant once it's saved here.

    `rows` is the (already sorted/--top-filtered) subset actually written;
    `total_peaks` is the full count found in range, for the header's own
    "N of TOTAL" note when --top trimmed the list.
    """
    lo, hi = two_theta_range
    header = [
        "# ============================================================\n",
        "# STB-XRD -- Simulated Powder Diffraction Pattern\n",
        "# ============================================================\n",
        f"# Structure       : {structure_label}\n",
        f"# Formula         : {formula}\n",
        f"# Space Group     : {space_group}\n",
        f"# Wavelength      : {wavelength_label} ({wavelength:.5f} Ang)\n",
        f"# 2-theta range   : {lo} - {hi} deg\n",
        f"# Peaks           : {len(rows)} of {total_peaks} (sorted by intensity, descending)\n",
        "# ------------------------------------------------------------\n",
        "# Columns: 2-theta(deg)  d(Ang)  h  k  l  Intensity(normalized)\n",
        "# ============================================================\n",
    ]
    with open(dat_path, "w") as f:
        f.writelines(header)
        for two_theta, d, h, k, l, intensity in rows:
            f.write(f"{two_theta:12.6f} {d:12.6f} {int(h):3d} {int(k):3d} {int(l):3d} {intensity:12.6f}\n")


def write_xrd_gnuplot(dat_path, gplot_path, formula, wavelength_label):
    """Writes a .gplot script plotting the simulated stick pattern (2-theta
    vs. intensity -- columns 1 and 6 of the .dat file pyxtal's own
    XRD.save() writes, per its documented 2theta/d/h/k/l/intensity column
    order) as impulses, the standard powder-XRD "stick pattern" convention.

    Companion to XRD.save() (called directly by stb-xrd, no wrapper needed
    here since pyxtal's own writer already matches this suite's plain
    -text-data-file convention) -- extracted alongside it once stb-xrd
    needed a gnuplot script to go with the data file it already writes.
    """
    dat_basename = os.path.basename(dat_path)
    lines = [
        'set terminal pdfcairo enhanced font "Arial,14" size 8,6\n',
        'set output "xrd_pattern.pdf"\n',
        'set xlabel "2{/Symbol Q} (deg)" font "Arial,16"\n',
        'set ylabel "Intensity (a.u.)" font "Arial,16"\n',
        f'set title "Simulated XRD pattern: {formula} ({wavelength_label})"\n',
        'set grid xtics ytics lt 0 lw 1 lc rgb "#bbbbbb"\n',
        'unset key\n',
        '\n',
        f'plot "{dat_basename}" using 1:6 with impulses lw 2 lc rgb "navy"\n',
    ]
    with open(gplot_path, "w") as f:
        f.writelines(lines)
