#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.9.1"

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from stb.core.cli import color_text, show_intro

# FWHM = 2*sqrt(2*ln2) * sigma, for a Gaussian.
FWHM_TO_SIGMA = 2 * np.sqrt(2 * np.log(2))


def gaussian_kernel_1d(size, sigma):
    """Creates a normalized 1D Gaussian kernel of the given odd size.
    `sigma` is in samples (grid points), not a physical energy unit -- see
    energy_grid_spacing() for the eV -> samples conversion callers are
    expected to do first (main() additionally takes --sigma/--fwhm in meV,
    converted to eV before that)."""
    kernel_1d = np.arange(size) - (size - 1) / 2.0
    kernel_1d = np.exp(-(kernel_1d**2) / (2 * sigma**2))
    kernel_1d /= np.sum(kernel_1d)  # Normalize the kernel
    return kernel_1d


def energy_grid_spacing(energy):
    """Median spacing (eV) between consecutive energy points, so --sigma/
    --fwhm can be given in physically meaningful units instead of grid
    samples -- whose count depends on how densely the input file happens to
    be sampled. Raises ValueError if the energy column isn't increasing
    (spacing would be zero or negative, meaning sigma-in-samples can't be
    computed)."""
    d_energy = float(np.median(np.diff(energy)))
    if d_energy <= 0:
        raise ValueError("Energy column must be sorted in strictly increasing order.")
    return d_energy


def grid_uniformity_warning(energy, d_energy, rel_tol=0.01):
    """The convolution advances one grid sample at a time, so --sigma/--fwhm
    (converted to samples via a single scalar spacing) only means what it
    says if the grid is uniformly spaced -- true by construction for
    stb-dos's own output (a linspace), but not guaranteed for a DOS file
    from elsewhere. Returns a message if consecutive spacings vary by more
    than `rel_tol` of the median spacing, else None."""
    spread = float(np.std(np.diff(energy)))
    if spread > rel_tol * abs(d_energy):
        return (f"energy grid spacing varies by ~{spread:.6f} eV (std) around the "
                f"median {d_energy:.6f} eV -- --sigma/--fwhm's physical broadening "
                "width will drift across the spectrum on a non-uniform grid.")
    return None


def kernel_size_for_sigma(sigma_samples):
    """Smallest odd kernel width covering about 3 standard deviations on
    each side of center (>99% of a Gaussian's mass), used when --size isn't
    given explicitly."""
    half_width = max(1, int(np.ceil(3 * sigma_samples)))
    return 2 * half_width + 1


def read_column_labels(path, n_columns):
    """Best-effort column labels from the input file's leading '#' header
    (e.g. stb-dos writes '#Energy(eV)    s           p           ...').
    Returns None if there's no header line or its label count doesn't match
    the data, so callers can fall back to generic names instead of
    mislabeling columns."""
    try:
        with open(path) as f:
            first_line = f.readline()
    except OSError:
        return None
    if not first_line.lstrip().startswith('#'):
        return None
    labels = first_line.lstrip('#').split()
    if len(labels) != n_columns:
        return None
    return labels


def filter_columns(data, kernel):
    """Gaussian-convolves every column except column 0 (energy).

    Uses numpy's own (vectorized, well-tested) convolution rather than a
    hand-rolled loop. mode='same' zero-pads at the boundaries -- the
    physically sensible choice here, since a DOS is genuinely ~0 outside
    the computed energy window.

    Returns (energy, filtered) where filtered has shape
    (n_energy_points, n_columns - 1).
    """
    energy = data[:, 0]
    filtered = np.empty((data.shape[0], data.shape[1] - 1))
    for i in range(1, data.shape[1]):
        filtered[:, i - 1] = np.convolve(data[:, i], kernel, mode='same')
    return energy, filtered


def print_conservation_check(energy, original, filtered, dos_labels):
    """Prints each column's integrated DOS (trapezoidal) before vs. after
    broadening -- a normalized kernel should barely change it, so this is a
    quick runtime sanity check that nothing went wrong (a large drift here
    usually means the energy window is too narrow for the requested
    --sigma/--fwhm, losing weight off the zero-padded edges)."""
    print("[INFO] DOS conservation check (integral before -> after broadening):")
    for i, label in enumerate(dos_labels):
        before = np.trapz(original[:, i + 1], energy)
        after = np.trapz(filtered[:, i], energy)
        print(f"       {label:<12}: {before:12.4f} -> {after:12.4f}")


def plot_all(energy, original, filtered, dos_labels):
    """One figure, one row of (original, filtered) subplots per DOS column
    -- a single window to close instead of one blocking plt.show() per
    column, and each subplot is titled with its actual DOS column label
    instead of a generic "Original"/"Filtered" that doesn't say which
    column it is."""
    n_cols = filtered.shape[1]
    fig, axes = plt.subplots(n_cols, 2, figsize=(10, 3 * n_cols), squeeze=False)
    for i in range(n_cols):
        label = dos_labels[i]
        axes[i, 0].plot(energy, original[:, i + 1])
        axes[i, 0].set_title(f"Original: {label}")
        axes[i, 0].set_xlabel("Energy")
        axes[i, 0].set_ylabel("DOS")

        axes[i, 1].plot(energy, filtered[:, i])
        axes[i, 1].set_title(f"Filtered: {label}")
        axes[i, 1].set_xlabel("Energy")
        axes[i, 1].set_ylabel("DOS")

    fig.tight_layout()
    plt.show()
    plt.close(fig)


def write_output(path, energy, filtered, energy_label, dos_labels, sigma_mev, sigma_samples, size):
    """Writes energy + all filtered DOS columns. Header line 1 lists every
    actual column (instead of a fixed 'Energy DOS_filtered' regardless of
    how many DOS columns there really are); line 2 records the broadening
    parameters actually used, so the file is self-describing even without
    the run's console log. Column labels stay on line 1 so a re-filtered
    file can still have its columns picked up by read_column_labels()."""
    columns_line = energy_label + " " + " ".join(f"{label}_filtered" for label in dos_labels)
    params_line = f"Gaussian broadening: sigma = {sigma_mev:.3f} meV ({sigma_samples:.3f} samples), kernel size = {size}"
    header = columns_line + "\n" + params_line
    out_array = np.column_stack([energy, filtered])
    np.savetxt(path, out_array, fmt='%.6f', header=header)


def main():
    parser = argparse.ArgumentParser(
        description="Apply Gaussian convolution to broaden a DOS file's columns.",
        epilog="Example usage:\n"
               "  stb-convdos --file dos_total.dat --sigma 50 --out dos_filtered.dat\n"
               "  stb-convdos --file dos_total.dat --fwhm 100 --size 41 --out dos_filtered.dat --no-plot",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--file", dest="input_file", required=True,
                        help="Input DOS file: whitespace-separated columns, energy first.")
    broadening = parser.add_mutually_exclusive_group(required=True)
    broadening.add_argument("--sigma", type=float, default=None,
                        help="Gaussian broadening standard deviation, in meV. Converted internally "
                             "to grid samples using the input file's own energy spacing, so the "
                             "same --sigma gives the same physical broadening regardless of how "
                             "densely the file is sampled. Must be positive. Alternative to --fwhm.")
    broadening.add_argument("--fwhm", type=float, default=None,
                        help="Gaussian broadening full width at half maximum, in meV "
                             "(FWHM = 2.3548 * sigma) -- the width more commonly quoted in "
                             "spectroscopy/DOS literature. Must be positive. Alternative to --sigma.")
    parser.add_argument("--size", type=int, default=None,
                        help="Gaussian kernel width, in samples. Optional -- by default sized "
                             "automatically from --sigma/--fwhm to cover about 3 standard "
                             "deviations on each side. Must be a positive odd number if given "
                             "explicitly (an even size shifts the filtered curve by half a bin, "
                             "since it has no single center sample).")
    parser.add_argument("--out", required=True, dest="outfile",
                        help="Output file for the filtered data.")
    parser.add_argument("--no-plot", dest="plot", action="store_false",
                        help="Skip the before/after plot (all columns are still filtered and written).")
    parser.add_argument("-v", "--version", action="version",
                        version=f"stb-convdos {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    args = parser.parse_args()

    if args.sigma is not None:
        if args.sigma <= 0:
            parser.error("--sigma must be positive.")
        sigma_mev = args.sigma
    else:
        if args.fwhm <= 0:
            parser.error("--fwhm must be positive.")
        sigma_mev = args.fwhm / FWHM_TO_SIGMA
    if args.size is not None and (args.size <= 0 or args.size % 2 == 0):
        parser.error("--size must be a positive odd number.")

    if args.intro == True:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])
    print("\n" + color_text("DOS Convolution Tool:", 'bold'))
    print("-"*60)

    print("\n[INFO] Reading file ...")
    inp_file = np.loadtxt(args.input_file)
    if inp_file.ndim != 2 or inp_file.shape[1] < 2:
        print(color_text("[ERROR] Input file needs at least 2 columns (energy + at least one DOS column).", 'red'))
        sys.exit(1)
    if inp_file.shape[0] < 2:
        print(color_text("[ERROR] Input file needs at least 2 energy points to determine grid spacing.", 'red'))
        sys.exit(1)

    labels = read_column_labels(args.input_file, inp_file.shape[1])
    energy_label = labels[0] if labels else "Energy"
    dos_labels = labels[1:] if labels else [f"DOS_{i}" for i in range(1, inp_file.shape[1])]

    try:
        d_energy = energy_grid_spacing(inp_file[:, 0])
    except ValueError as e:
        print(color_text(f"[ERROR] {e}", 'red'))
        sys.exit(1)

    uniformity_warning = grid_uniformity_warning(inp_file[:, 0], d_energy)
    if uniformity_warning:
        print(color_text(f"[WARNING] {uniformity_warning}", 'yellow'))

    sigma_ev = sigma_mev / 1000.0
    sigma_samples = sigma_ev / d_energy
    size = args.size if args.size is not None else kernel_size_for_sigma(sigma_samples)
    print(f"[INFO] Energy grid spacing: {d_energy:.6f} eV -> sigma = {sigma_mev:.3f} meV "
          f"= {sigma_samples:.3f} samples, kernel size = {size}")
    if size > inp_file.shape[0]:
        print(color_text(f"[ERROR] Kernel size ({size}) is larger than the number of energy "
                          f"points ({inp_file.shape[0]}) in the file -- double-check "
                          "--sigma/--fwhm and the file's energy units (a much bigger value "
                          "than intended, e.g. eV instead of meV, is the usual cause).", 'red'))
        sys.exit(1)

    kernel = gaussian_kernel_1d(size, sigma_samples)
    print("[INFO] Applying Gaussian convolution ...")
    energy, filtered = filter_columns(inp_file, kernel)
    print_conservation_check(energy, inp_file, filtered, dos_labels)

    # Write the actual deliverable before plotting, so a display/backend
    # problem (e.g. no X11 and MPLBACKEND unset) can't cost the already-done
    # convolution work -- the result file is saved either way.
    write_output(args.outfile, energy, filtered, energy_label, dos_labels, sigma_mev, sigma_samples, size)
    print(f"\n[OK] Filtered data written to {args.outfile}")

    if args.plot:
        n_cols = filtered.shape[1]
        if n_cols > 12:
            print(f"[INFO] {n_cols} columns to plot -- the combined figure will be tall; "
                  "pass --no-plot to skip it.")
        print("[INFO] Plotting original vs. filtered ...")
        plot_all(energy, inp_file, filtered, dos_labels)
    print("[INFO] Complete job!")
    print("\n"+"-"*60)
    print(color_text("We’ve convoluted everything… including the soul of the electron.\n\n", 'bold'))

if __name__ == "__main__":
    main()
