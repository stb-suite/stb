#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.11.0"

import sys
import argparse
import numpy as np
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.diffraction.xrd import WAVELENGTHS
from stb.core import structure_io
from stb.core.cli import color_text, show_intro, run_with_spinner
from stb.core.deps import require_pyxtal


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


MIN_EXPERIMENTAL_POINTS = 4  # pyxtal.XRD.Similarity cubic-interpolates each pattern (scipy
                             # interp1d, kind='cubic'), which needs at least 4 points per curve
                             # -- fewer than that raises a raw, unhandled scipy ValueError deep
                             # inside Similarity() instead of this tool's own clean error.


def read_experimental_pattern(path):
    """Reads a plain 2-column (2theta, intensity) text file for --compare-to:
    whitespace- or comma-separated, blank lines and '#' comments skipped.
    Returns a (2, N) array -- the same shape pyxtal.XRD.Similarity expects
    for both patterns it compares (it interpolates internally, so this
    doesn't need to share a grid with the simulated pattern).
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


def main():
    require_pyxtal()
    from pyxtal.XRD import XRD

    parser = argparse.ArgumentParser(
        description=f"""{color_text("Simulates a powder XRD pattern from a structure.", 'bold')}
Uses pyxtal's own diffraction engine (not pymatgen's) so the same dependency
already needed by stb-crystalcast covers this too.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s --file structure.fdf --format fdf\n"
               "  %(prog)s --file structure.fdf --format fdf --wavelength MoKa --top 10\n"
               "  %(prog)s --file relaxed.STRUCT_OUT --format struct_out \\\n"
               "      --two-theta-range 10 80 --plot\n"
               "  %(prog)s --file structure.fdf --format fdf \\\n"
               "      --compare-to experimental.dat --plot\n"
    )

    parser.add_argument("--file", required=True, help="Path to structure file.")
    parser.add_argument("--format", required=True, choices=["fdf", "struct_out"],
                        help="Input file format:\n"
                             "  fdf:        SIESTA structure input (%%block LatticeVectors etc.)\n"
                             "  struct_out: SIESTA post-relaxation output (.STRUCT_OUT)")
    parser.add_argument("--wavelength", type=str, default="CuKa",
                        help="X-ray source: a known name (CuKa, CuKa1, MoKa, CrKa, FeKa, CoKa, "
                             "AgKa, ...) or a wavelength in Ang as a plain number "
                             "(default: CuKa, 1.54184 Ang).")
    parser.add_argument("--two-theta-range", type=float, nargs=2, default=[0.0, 90.0],
                        metavar=("MIN", "MAX"),
                        help="2-theta range in degrees to compute (default: 0 90).")
    parser.add_argument("--top", type=int, default=None,
                        help="Only show the N strongest peaks in the printed table (default: "
                             "show all of them). The output data file always has all of them.")
    parser.add_argument("--plot", action="store_true",
                        help="Show an interactive plot of the pattern with hkl labels (or, with "
                             "--compare-to, an overlay of both patterns instead).")
    parser.add_argument("--compare-to", type=str, default=None, metavar="PATH",
                        help="Compare the simulated pattern against an experimental one: a plain "
                             "text file with two columns (2theta intensity), whitespace- or "
                             "comma-separated, '#' comments and blank lines allowed. Prints a "
                             "similarity score (pyxtal's cosine-weighted metric, 0-1, higher is "
                             "more similar); combine with --plot for a visual overlay.")
    parser.add_argument("-o", "--output", type=str, default="xrd_pattern.dat",
                        help="Output data file name -- columns are 2theta, d, h, k, l, "
                             "intensity (default: xrd_pattern.dat).")
    parser.add_argument("-v", "--version", action="version", version=f"stb-xrd {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("Simulate a powder XRD pattern:", 'bold'))
    print("-" * 60)

    try:
        structure = structure_io.read_siesta_structure(args.file, args.format)
    except FileNotFoundError:
        print(color_text(f"Error: structure file '{args.file}' not found.", 'red'))
        sys.exit(1)
    except ValueError as e:
        print(color_text(f"Error: {e}", 'red'))
        sys.exit(1)

    try:
        wavelength, wavelength_label = resolve_wavelength(args.wavelength)
    except ValueError as e:
        print(color_text(f"Error: {e}", 'red'))
        sys.exit(1)

    lo, hi = args.two_theta_range
    if lo >= hi:
        print(color_text("Error: --two-theta-range MIN must be less than MAX.", 'red'))
        sys.exit(1)

    print(f"  {color_text('Input file:', 'cyan')} {args.file}")
    print(f"  {color_text('Formula:', 'cyan')} {structure.composition.reduced_formula}")
    print(f"  {color_text('Wavelength:', 'cyan')} {wavelength_label} ({wavelength:.5f} Ang)")
    print(f"  {color_text('2-theta range:', 'cyan')} {lo} - {hi} deg")

    atoms = AseAtomsAdaptor.get_atoms(structure)
    try:
        xrd = XRD(atoms, wavelength=wavelength, thetas=[lo, hi])
    except ValueError:
        # pyxtal's own XRD.__init__ calls max() on an empty list when no
        # reflection falls inside [lo, hi], raising a bare ValueError instead
        # of returning an empty pattern.
        print(color_text(
            "\nError: no peaks found in the given --two-theta-range -- try widening it.", 'red'))
        sys.exit(1)

    if len(xrd.pxrd) == 0:
        print(color_text(
            "\nError: no peaks found in the given --two-theta-range -- try widening it.", 'red'))
        sys.exit(1)

    rows = sorted(xrd.pxrd, key=lambda row: row[5], reverse=True)
    shown = rows if args.top is None else rows[:args.top]

    print(f"\n{color_text(f'Peaks ({len(shown)} of {len(rows)} shown):', 'bold')}")
    print(f"  {'2-theta':>8}  {'d (Ang)':>8}  {'h':>3} {'k':>3} {'l':>3}  {'Intensity':>10}")
    for two_theta, d, h, k, l, intensity in shown:
        print(f"  {two_theta:8.3f}  {d:8.4f}  {int(h):3d} {int(k):3d} {int(l):3d}  {intensity:10.2f}")

    xrd.save(args.output)
    print(f"\n{color_text('Success:', 'green')} pattern data ({len(rows)} peaks) written to "
          f"'{color_text(args.output, 'bold')}'")

    similarity = None
    sim_profile = None
    if args.compare_to:
        from pyxtal.XRD import Similarity

        try:
            experimental = read_experimental_pattern(args.compare_to)
        except (FileNotFoundError, ValueError) as e:
            print(color_text(f"Error: {e}", 'red'))
            sys.exit(1)

        print(f"\n  {color_text('Compared to:', 'cyan')} {args.compare_to}")
        sim_profile = xrd.get_profile()
        similarity = run_with_spinner(
            lambda: Similarity(sim_profile, experimental).value, label="Computing similarity")
        print(f"  {color_text('Similarity score:', 'cyan')} {similarity:.4f} "
              "(0-1, cosine-weighted, higher is more similar)")

    if args.plot:
        if args.compare_to:
            import matplotlib.pyplot as plt

            sim_x, sim_y = sim_profile
            sim_y = sim_y / sim_y.max() if sim_y.max() > 0 else sim_y
            exp_x, exp_y = experimental
            exp_y = exp_y / exp_y.max() if exp_y.max() > 0 else exp_y

            # fig.show() (non-blocking), not the pyplot-level plt.show() --
            # matches pyxtal's own plot_pxrd() below, so both --plot paths in
            # this tool behave the same way instead of one hanging without a
            # running GUI event loop and the other returning immediately.
            fig, ax = plt.subplots(figsize=(10, 5))
            ax.plot(sim_x, sim_y, label="Simulated", color="C0")
            ax.plot(exp_x, exp_y, label=f"Experimental ({args.compare_to})", color="C1", alpha=0.7)
            ax.set_xlabel("2-theta (deg)")
            ax.set_ylabel("Normalized intensity")
            ax.set_title(f"Similarity: {similarity:.4f}")
            ax.set_xlim(lo, hi)
            ax.legend()
            fig.show()
        else:
            xrd.plot_pxrd(show=True)


if __name__ == "__main__":
    main()
