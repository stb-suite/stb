#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "2.0.0"

import os
import sys
import glob
import argparse
from datetime import datetime
import numpy as np
# np.trapz was removed in NumPy 2.0; use scipy's trapezoid for compatibility
# with both old and new NumPy (same fix already applied in strain_analysis.py).
from scipy.integrate import trapezoid
import matplotlib.pyplot as plt
from stb.core import citations
from stb.core.cli import color_text, show_intro, print_dual, print_section, print_table

REPORT_FILE = "stb_convdos_report.txt"
BIB_FILE = "references.bib"

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


def plot_all(energy, original, filtered, dos_labels, title=None):
    """One figure, one row of (original, filtered) subplots per DOS column
    -- a single window to close instead of one blocking plt.show() per
    column, and each subplot is titled with its actual DOS column label
    instead of a generic "Original"/"Filtered" that doesn't say which
    column it is."""
    n_cols = filtered.shape[1]
    fig, axes = plt.subplots(n_cols, 2, figsize=(10, 3 * n_cols), squeeze=False)
    if title:
        fig.suptitle(title)
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


# --- Pure per-file computation -------------------------------------------
# Everything above (and this function) is unchanged physics/math from the
# original single-file tool -- only pulled out of main() into its own
# function, with no printing/file writing, so both --file and the new
# --dir (batch) mode below can share the exact same pipeline instead of
# --dir reimplementing it. Raises ValueError on any per-file validation
# failure (caller decides: --file aborts the run, --dir skips the file
# with a warning and continues).

def process_one_file(path, sigma_mev, size_override=None):
    inp_file = np.loadtxt(path)
    if inp_file.ndim != 2 or inp_file.shape[1] < 2:
        raise ValueError("Input file needs at least 2 columns (energy + at least one DOS column).")
    if inp_file.shape[0] < 2:
        raise ValueError("Input file needs at least 2 energy points to determine grid spacing.")

    labels = read_column_labels(path, inp_file.shape[1])
    energy_label = labels[0] if labels else "Energy"
    dos_labels = labels[1:] if labels else [f"DOS_{i}" for i in range(1, inp_file.shape[1])]

    d_energy = energy_grid_spacing(inp_file[:, 0])
    uniformity_warning = grid_uniformity_warning(inp_file[:, 0], d_energy)

    sigma_ev = sigma_mev / 1000.0
    sigma_samples = sigma_ev / d_energy
    size = size_override if size_override is not None else kernel_size_for_sigma(sigma_samples)
    if size > inp_file.shape[0]:
        raise ValueError(
            f"Kernel size ({size}) is larger than the number of energy points "
            f"({inp_file.shape[0]}) in the file -- double-check --sigma/--fwhm and "
            "the file's energy units (a much bigger value than intended, e.g. eV "
            "instead of meV, is the usual cause)."
        )

    kernel = gaussian_kernel_1d(size, sigma_samples)
    energy, filtered = filter_columns(inp_file, kernel)

    conservation = [
        (label, float(trapezoid(inp_file[:, i + 1], energy)), float(trapezoid(filtered[:, i], energy)))
        for i, label in enumerate(dos_labels)
    ]

    return {
        "path": path,
        "energy": energy,
        "original": inp_file,
        "filtered": filtered,
        "energy_label": energy_label,
        "dos_labels": dos_labels,
        "d_energy": d_energy,
        "uniformity_warning": uniformity_warning,
        "sigma_mev": sigma_mev,
        "sigma_samples": sigma_samples,
        "size": size,
        "conservation": conservation,
        "n_points": inp_file.shape[0],
        "n_columns": inp_file.shape[1],
    }


def main():
    parser = argparse.ArgumentParser(
        description="Apply Gaussian convolution to broaden a DOS file's columns.",
        epilog="Example usage:\n"
               "  stb-convdos --file dos_total.dat --sigma 50 --out dos_filtered.dat\n"
               "  stb-convdos --file dos_total.dat --fwhm 100 --size 41 --out dos_filtered.dat --no-plot\n"
               "  stb-convdos --dir dos_output/ --sigma 50 --output-dir dos_output_filtered/",
        formatter_class=argparse.RawTextHelpFormatter
    )
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--file", dest="input_file", default=None,
                        help="Input DOS file: whitespace-separated columns, energy first. "
                             "Requires --out. Alternative to --dir.")
    input_group.add_argument("--dir", dest="input_dir", default=None,
                        help="Recursively broaden every *.dat file found under this directory "
                             "(e.g. a stb-dos output folder: dos_total.dat + dos_per_atom/*.dat "
                             "+ dos_per_species/*.dat) with the same --sigma/--fwhm/--size. Each "
                             "file is written under --output-dir preserving its relative path "
                             "and filename -- a broadened mirror of the whole input tree. A file "
                             "that fails validation is skipped with a warning, not fatal to the "
                             "rest of the batch. Alternative to --file.")
    broadening = parser.add_mutually_exclusive_group(required=True)
    broadening.add_argument("--sigma", type=float, default=None,
                        help="Gaussian broadening standard deviation, in meV. Converted internally "
                             "to grid samples using each input file's own energy spacing, so the "
                             "same --sigma gives the same physical broadening regardless of how "
                             "densely a file is sampled. Must be positive. Alternative to --fwhm.")
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
    parser.add_argument("--out", dest="outfile", default=None,
                        help="Output file for the filtered data. Required with --file, invalid with --dir.")
    parser.add_argument("-o", "--output-dir", dest="output_dir", default=None,
                        help="Root of the broadened output tree for --dir (default: "
                             "'<dir>_filtered', a sibling directory). Invalid with --file "
                             "(use --out there).")
    parser.add_argument("--no-plot", dest="plot", action="store_false",
                        help="Skip the before/after plot (all columns are still filtered and "
                             "written). With --dir, only the first processed file is ever "
                             "plotted regardless of this flag's default; --no-plot turns that "
                             "single plot off too.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the full run report to {REPORT_FILE}. Off by default.")
    parser.add_argument("-v", "--version", action="version",
                        version=f"stb-convdos {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    args = parser.parse_args()

    if args.input_file and not args.outfile:
        parser.error("--out is required with --file.")
    if args.input_file and args.output_dir:
        parser.error("--output-dir is only valid with --dir (use --out with --file).")
    if args.input_dir and args.outfile:
        parser.error("--out is only valid with --file (use --output-dir with --dir).")

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

    # --- Resolve input mode and read/validate BEFORE writing anything,
    # so a bad input can't leave a partial output tree/empty report
    # behind (same convention as bands.py/dos.py/structural.py). ---
    if args.input_file:
        mode = "file"
        try:
            result = process_one_file(args.input_file, sigma_mev, args.size)
        except FileNotFoundError:
            print(color_text(f"[ERROR] File not found: {args.input_file}", 'red'))
            sys.exit(1)
        except ValueError as e:
            print(color_text(f"[ERROR] {e}", 'red'))
            sys.exit(1)
        results = [result]
        skipped = []
        output_dir = None
    else:
        mode = "dir"
        if not os.path.isdir(args.input_dir):
            print(color_text(f"[ERROR] Directory not found: {args.input_dir}", 'red'))
            sys.exit(1)
        files = sorted(glob.glob(os.path.join(args.input_dir, "**", "*.dat"), recursive=True))
        if not files:
            print(color_text(f"[ERROR] No .dat files found under '{args.input_dir}'.", 'red'))
            sys.exit(1)

        results = []
        skipped = []
        for f in files:
            try:
                results.append(process_one_file(f, sigma_mev, args.size))
            except (ValueError, FileNotFoundError, OSError) as e:
                skipped.append((f, str(e)))
        if not results:
            print(color_text(f"[ERROR] None of the {len(files)} .dat file(s) under "
                              f"'{args.input_dir}' could be processed.", 'red'))
            sys.exit(1)

        output_dir = args.output_dir or (args.input_dir.rstrip('/\\') + "_filtered")

    # --- Write output files ---
    written_files = []
    if mode == "file":
        out_dir = os.path.dirname(args.outfile)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        r = results[0]
        write_output(args.outfile, r["energy"], r["filtered"], r["energy_label"], r["dos_labels"],
                     r["sigma_mev"], r["sigma_samples"], r["size"])
        written_files.append(args.outfile)
    else:
        for r in results:
            rel_path = os.path.relpath(r["path"], args.input_dir)
            out_path = os.path.join(output_dir, rel_path)
            out_subdir = os.path.dirname(out_path)
            if out_subdir:
                os.makedirs(out_subdir, exist_ok=True)
            write_output(out_path, r["energy"], r["filtered"], r["energy_label"], r["dos_labels"],
                         r["sigma_mev"], r["sigma_samples"], r["size"])
            written_files.append(out_path)

    report_path = os.path.join(output_dir if mode == "dir" else (os.path.dirname(args.outfile) or "."),
                                REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    print_dual(color_text("===== STB-CONVDOS REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time      : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    if mode == "file":
        print_dual(f"Mode           : single file", f_out)
        print_dual(f"Input file     : {args.input_file}", f_out)
        print_dual(f"Output file    : {args.outfile}", f_out)
    else:
        print_dual(f"Mode           : directory (recursive)", f_out)
        print_dual(f"Input dir      : {args.input_dir}", f_out)
        print_dual(f"Output dir     : {output_dir}", f_out)
    broadening_desc = f"sigma = {args.sigma} meV" if args.sigma is not None else f"fwhm = {args.fwhm} meV"
    print_dual(f"Broadening     : {broadening_desc} (-> sigma = {sigma_mev:.3f} meV)", f_out)
    print_dual(f"Kernel size    : {'explicit ' + str(args.size) if args.size is not None else 'auto (3-sigma rule)'}", f_out)
    print_dual(f"Plot           : {'on' if args.plot else 'off'}"
               + (" (--dir: first processed file only)" if mode == "dir" and args.plot else ""), f_out)

    print_section("[1] INPUT DATA", f_out)
    if mode == "file":
        r = results[0]
        print_table(["Quantity", "Value"], [
            (["Energy points", str(r["n_points"])], None),
            (["Columns", f"{r['n_columns']} ({r['energy_label']} + {len(r['dos_labels'])} DOS)"], None),
            (["Grid spacing", f"{r['d_energy']:.6f} eV"], None),
        ], f_out)
        if r["uniformity_warning"]:
            print_dual(color_text(f"[WARNING] {r['uniformity_warning']}", 'yellow'), f_out)
    else:
        rows = [([os.path.relpath(r["path"], args.input_dir), "OK", str(r["n_columns"]), str(r["n_points"])], None)
                for r in results]
        rows += [([os.path.relpath(p, args.input_dir), f"Skipped ({msg})", "--", "--"], 'yellow')
                 for p, msg in skipped]
        print_table(["File", "Status", "Columns", "Points"], rows, f_out)
        for r in results:
            if r["uniformity_warning"]:
                print_dual(color_text(
                    f"[WARNING] {os.path.relpath(r['path'], args.input_dir)}: {r['uniformity_warning']}",
                    'yellow'), f_out)

    print_section("[2] BROADENING PARAMETERS", f_out)
    sizes = {r["size"] for r in results}
    samples = {round(r["sigma_samples"], 6) for r in results}
    if mode == "file" or (len(sizes) == 1 and len(samples) == 1):
        r = results[0]
        print_table(["Quantity", "Value"], [
            (["Sigma", f"{r['sigma_mev']:.3f} meV ({r['sigma_samples']:.3f} samples)"], None),
            (["Kernel size", f"{r['size']} samples"], None),
        ], f_out)
        if mode == "dir":
            print_dual(f"(same for all {len(results)} processed files -- they share one energy grid)", f_out)
    else:
        rows = [([os.path.relpath(r["path"], args.input_dir), f"{r['sigma_samples']:.3f}", str(r["size"])], None)
                for r in results]
        print_dual(f"Sigma = {sigma_mev:.3f} meV requested for every file; samples/kernel size "
                   "vary because their energy grids differ:", f_out)
        print_table(["File", "Sigma (samples)", "Kernel size"], rows, f_out)

    print_section("[3] CONSERVATION CHECK", f_out)
    print_dual("Integrated DOS (trapezoidal) before -> after broadening -- a normalized", f_out)
    print_dual("kernel should barely change this; a large drift usually means the energy", f_out)
    print_dual("window is too narrow for the requested broadening (weight lost off the", f_out)
    print_dual("zero-padded edges).", f_out)
    if mode == "file":
        rows = [([label, f"{before:.4f}", f"{after:.4f}"], None) for label, before, after in results[0]["conservation"]]
        print_table(["Column", "Before", "After"], rows, f_out)
    else:
        rows = [([os.path.relpath(r["path"], args.input_dir), label, f"{before:.4f}", f"{after:.4f}"], None)
                for r in results for label, before, after in r["conservation"]]
        print_table(["File", "Column", "Before", "After"], rows, f_out)

    print_section("[4] WRITING OUTPUT FILES", f_out)
    for path in written_files:
        print_dual(color_text(f"[OK] Filtered data written to {path}", 'green'), f_out)

    print_section("[5] REFERENCES", f_out)
    bib_dir = output_dir if mode == "dir" else (os.path.dirname(args.outfile) or ".")
    bib_path = os.path.join(bib_dir, BIB_FILE)
    bib_entries = [citations.SIESTA, citations.SIESTA_RECENT]
    citations.write_bib_file(bib_path, bib_entries)
    print_dual(color_text(
        f"[OK] Citations for the methods used in this run written to "
        f"'{bib_path}' ({len(bib_entries)} entries).", 'green'), f_out)

    print_section("[6] SUMMARY & FILES", f_out)
    print_dual("Status         : OK", f_out)
    if skipped:
        print_dual(color_text(f"Skipped        : {len(skipped)} of {len(results) + len(skipped)} "
                               "file(s) (see [1] INPUT DATA)", 'yellow'), f_out)
    for path in written_files:
        print_dual(f"Data file      : {path}", f_out)
    print_dual(f"References     : {bib_path}", f_out)
    if report_path:
        print_dual(f"Report         : {report_path}", f_out)

    if f_out:
        f_out.close()

    # Interactive matplotlib preview runs last, after every report/file
    # write above, so a blocking GUI window never delays or hides them.
    # --dir only ever plots the first processed file (see --no-plot help).
    if args.plot:
        if mode == "file":
            n_cols = results[0]["filtered"].shape[1]
            if n_cols > 12:
                print(f"[INFO] {n_cols} columns to plot -- the combined figure will be tall; "
                      "pass --no-plot to skip it.")
            plot_all(results[0]["energy"], results[0]["original"], results[0]["filtered"],
                     results[0]["dos_labels"])
        else:
            r = results[0]
            print(f"[INFO] Plotting only the first processed file ({os.path.relpath(r['path'], args.input_dir)}) "
                  f"-- pass --no-plot to skip it too.")
            plot_all(r["energy"], r["original"], r["filtered"], r["dos_labels"],
                     title=os.path.relpath(r["path"], args.input_dir))


if __name__ == "__main__":
    main()
