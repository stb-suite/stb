#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.0.0"

import os
import re
import sys
import json
import shutil
import argparse
from stb.core import structure_io
from stb.core.cli import COLORS, color_text, show_intro

# Standard correlated shell per element (periodic-table block classification --
# uncontroversial chemistry, unlike the Hubbard U value itself). Used only as
# a default when --shell isn't given; --shell always overrides.
SHELL_NAMES = {"3d": (3, 2), "4d": (4, 2), "5d": (5, 2), "4f": (4, 3), "5f": (5, 3)}

_3D = ["Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn"]
_4D = ["Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd"]
_5D = ["La", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg"]
_4F = ["Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]
_5F = ["Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"]

DEFAULT_SHELL = {
    **{el: "3d" for el in _3D}, **{el: "4d" for el in _4D}, **{el: "5d" for el in _5D},
    **{el: "4f" for el in _4F}, **{el: "5f" for el in _5F},
}

# Default perturbation strengths (eV) for the linear-response fit -- small
# enough to stay in the linear regime, symmetric around zero so the fit isn't
# skewed by an accidental asymmetry in the response.
DEFAULT_ALPHAS = [-0.15, -0.10, -0.05, 0.05, 0.10, 0.15]

_LABEL_RE = re.compile(r'^\s*SystemLabel\s+(\S+)', re.IGNORECASE | re.MULTILINE)


def get_system_label(calc_text):
    match = _LABEL_RE.search(calc_text)
    return match.group(1) if match else "siesta"


def ldau_proj_block(species, n, l, u, potential_shift):
    """The %block LDAU.proj stanza plus, when potential_shift, the
    LDAU.PotentialShift flag that makes SIESTA treat `u` as a rigid
    perturbation on the shell (Cococcioni & de Gironcoli, PRB 71, 035105,
    2005) instead of the real Hubbard term -- and, as a side effect, makes it
    print the shell's occupation to the .out file (verified against SIESTA's
    own source, Src/dftu.F, where that print statement is gated on this exact
    flag).
    """
    lines = []
    if potential_shift:
        lines.append("LDAU.PotentialShift T\n")
    lines.append("%block LDAU.proj\n")
    lines.append(f"{species}   1\n")
    lines.append(f"n={n}    {l}\n")
    lines.append(f"{u:.3f}    0.000\n")
    lines.append("0.000    0.000\n")
    lines.append("%endblock LDAU.proj\n")
    return "".join(lines)


def write_run_folder(output_dir, folder_name, calc_template, structure_path, extra_fdf_text):
    folder = os.path.join(output_dir, folder_name)
    os.makedirs(folder, exist_ok=True)
    shutil.copy(structure_path, os.path.join(folder, "structure.fdf"))
    with open(os.path.join(folder, "calc.fdf"), 'w') as f:
        f.write(calc_template)
        f.write("\n# --- stb-hubbardu ---\n")
        f.write(extra_fdf_text)
    return folder


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Generates the run folders for a Hubbard U linear-response calculation (Cococcioni & de Gironcoli, Phys. Rev. B 71, 035105, 2005).", 'bold')}
For a chosen species/shell, writes: a 'reference' folder (unperturbed, needed
to seed the frozen-density runs), a 'scf_alpha_<v>' folder per perturbation
strength (full self-consistent response), and a 'frozen_alpha_<v>' folder per
perturbation strength (single-iteration, frozen-density response). Doesn't
run SIESTA itself -- run each folder's calculation yourself, then use
stb-hubbarduAnalysis.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage example:\n"
               "  %(prog)s -s structure.fdf -c calc.fdf --species Mn\n"
               "  %(prog)s -s structure.fdf -c calc.fdf --species Fe --shell 3d "
               "--alphas -0.2 -0.1 0.1 0.2\n"
    )

    parser.add_argument("-s", "--structure", type=str, required=True,
                        help="Input structure.fdf, copied as 'structure.fdf' into every run folder.")
    parser.add_argument("-c", "--calc", type=str, required=True,
                        help="Existing calc.fdf template (with %%include structure.fdf and everything "
                             "else already configured); the DFT+U perturbation lines are appended.")
    parser.add_argument("--species", type=str, required=True,
                        help="Species label to apply the Hubbard U correction to (must be present in "
                             "--structure).")
    parser.add_argument("--shell", type=str, default=None, choices=sorted(SHELL_NAMES),
                        help="Correlated shell (3d/4d/5d/4f/5f). Default: the standard shell for "
                             "--species (transition metals -> (n)d, lanthanides -> 4f, actinides -> 5f).")
    parser.add_argument("--alphas", type=float, nargs='+', default=DEFAULT_ALPHAS,
                        help=f"Perturbation strengths in eV (default: {DEFAULT_ALPHAS}).")
    parser.add_argument("-o", "--output-dir", type=str, default="hubbardu_runs",
                        help="Output directory for the run folders (default: hubbardu_runs).")
    parser.add_argument("-v", "--version", action="version", version=f"stb-hubbardu {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("Generate a Hubbard U linear-response sweep:", 'bold'))
    print("-" * 60)

    if not os.path.exists(args.structure):
        print(color_text(f"Error: File '{args.structure}' not found.", 'red'))
        sys.exit(1)
    if not os.path.exists(args.calc):
        print(color_text(f"Error: File '{args.calc}' not found.", 'red'))
        sys.exit(1)

    try:
        structure = structure_io.read_fdf(args.structure)
    except (FileNotFoundError, ValueError) as e:
        print(color_text(f"Error: {e}", 'red'))
        sys.exit(1)

    if args.species not in structure.species:
        print(color_text(
            f"Error: Species '{args.species}' not found in '{args.structure}' "
            f"(available: {', '.join(structure.species)}).", 'red'))
        sys.exit(1)

    shell_name = args.shell or DEFAULT_SHELL.get(args.species)
    if shell_name is None:
        print(color_text(
            f"Error: No default correlated shell known for '{args.species}' -- "
            f"pass --shell explicitly ({', '.join(sorted(SHELL_NAMES))}).", 'red'))
        sys.exit(1)
    n, l = SHELL_NAMES[shell_name]

    if len(set(args.alphas)) != len(args.alphas):
        print(color_text("Error: --alphas contains duplicate values.", 'red'))
        sys.exit(1)
    if 0.0 in args.alphas:
        print(color_text(
            "Error: 0.0 is not a valid perturbation strength -- alpha=0 is already "
            "covered by the 'reference' folder.", 'red'))
        sys.exit(1)

    with open(args.calc, 'r') as f:
        calc_template = f.read()
    label = get_system_label(calc_template)

    print(f"  {color_text('Species / shell:', 'cyan')} {args.species} ({shell_name}: n={n}, l={l})")
    print(f"  {color_text('Perturbations (eV):', 'cyan')} {args.alphas}")

    os.makedirs(args.output_dir, exist_ok=True)
    manifest = {"species": args.species, "n": n, "l": l, "label": label, "runs": {}}

    ref_extra = ldau_proj_block(args.species, n, l, 0.0, potential_shift=True)
    ref_folder = write_run_folder(args.output_dir, "reference", calc_template, args.structure, ref_extra)
    manifest["runs"]["reference"] = {"kind": "reference", "alpha": 0.0}
    print(f"  {color_text('[OK]', 'green')} {ref_folder} (alpha=0, reference)")

    for alpha in args.alphas:
        scf_name = f"scf_alpha_{alpha:.4f}"
        scf_extra = ldau_proj_block(args.species, n, l, alpha, potential_shift=True)
        scf_folder = write_run_folder(args.output_dir, scf_name, calc_template, args.structure, scf_extra)
        manifest["runs"][scf_name] = {"kind": "scf", "alpha": alpha}
        print(f"  {color_text('[OK]', 'green')} {scf_folder} (self-consistent response)")

        frozen_name = f"frozen_alpha_{alpha:.4f}"
        frozen_extra = (ldau_proj_block(args.species, n, l, alpha, potential_shift=True)
                         + "MaxSCFIterations 1\n" + "DM.UseSaveDM T\n")
        frozen_folder = write_run_folder(args.output_dir, frozen_name, calc_template, args.structure, frozen_extra)
        manifest["runs"][frozen_name] = {"kind": "frozen", "alpha": alpha}
        print(f"  {color_text('[OK]', 'green')} {frozen_folder} (frozen-density response)")

    manifest_path = os.path.join(args.output_dir, "run_manifest.json")
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)

    print(f"\n{color_text('Success:', 'green')} Generated {1 + 2*len(args.alphas)} run(s) in "
          f"'{color_text(args.output_dir, 'bold')}'")
    print(color_text("\nNext steps:", 'yellow'))
    print(f"  1. Run SIESTA in '{args.output_dir}/reference/' first.")
    print(f"  2. Copy the resulting '{label}.DM' from 'reference/' into every "
          f"'{args.output_dir}/frozen_alpha_*' folder.")
    print(f"  3. Run SIESTA in every 'scf_alpha_*' and 'frozen_alpha_*' folder.")
    print(f"  4. Run stb-hubbarduAnalysis --dir {args.output_dir}")


if __name__ == "__main__":
    main()
