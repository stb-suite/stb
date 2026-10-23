#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""Stage 1 of the real-SIESTA equation-of-state workflow (4.18): writes one run
folder per scanned volume -- the DFT analog of stb-mleos, which does the same
isotropic volume scan with a MACE potential in memory.

Every folder `vol_<pct>/` holds the structure with ONLY %block LatticeVectors
rewritten (isotropically scaled by factor = (1 + strain_pct/100)**(1/3), the
convention stb-mleos and stb-mlphonons --qha already use, so --strain-range means
the same thing in all three; atomic positions are fractional and scale for free),
the user's own calc.fdf template UNTOUCHED except for a `%include
config_extra.fdf` at its very top, the config_extra.fdf sidecar with every
directive the scan needs forced (fixed cell, single point -- or fixed-cell ionic
relaxation with --relax-positions --, and the SCF k-grid when the template has
none), and the pseudopotentials (-p: a bundled bank or a folder), the same
conventions as stb-optical and the other workflow prep stages.

`eos_scan.json` next to the folders records the scan (reference structure, planned
volumes, forced directives, level of theory) so that stb-eosAnalysis can compare
what was planned with what SIESTA actually did.

Bulk (3D periodic) only: an equation of state needs a bounded cell volume.
The scan scales the cell isotropically and keeps its shape -- for a non-cubic
crystal it gives E(V) at the input c/a (and b/a), so its B0 is the bulk modulus
of THAT constrained path, not of the fully relaxed solid; the report says so.
"""

from stb import __version__ as VERSION

import os
import re
import sys
import json
import argparse
import datetime

import numpy as np

from stb.core import structure_io, kspace, fdf_directives
from stb.core.cli import color_text, show_intro, print_dual, print_section, print_table
from stb.core.optical_setup import reciprocal_spacings
from stb.core.pseudopotentials import resolve_pseudo_source, copy_pseudo, pseudo_valence

REPORT_FILE = "stb_eosInputs_report.txt"
MANIFEST_FILE = "eos_scan.json"
CONFIG_EXTRA_FILE = "config_extra.fdf"
DEFAULT_KGRID_DENSITY = 0.2     # 1/Ang, same convention as stb-kgrid and the other workflows
RELAX_FORCE_TOL = "0.01 eV/Ang"

_SINGLE_POINT_EXTRA = (
    "# Fixed cell, single point: every folder is one independent sample of E(V).\n"
    "MD.TypeOfRun          CG\n"
    "MD.Steps              0\n"
    "MD.NumCGsteps         0\n"
    "MD.VariableCell       false\n"
)
_RELAX_EXTRA = (
    "# Fixed cell, ionic relaxation: the atoms may move, the cell may NOT (that would undo the scan).\n"
    "MD.TypeOfRun          CG\n"
    "MD.Steps              {steps}\n"
    "MD.NumCGsteps         {steps}\n"
    f"MD.MaxForceTol        {RELAX_FORCE_TOL}\n"
    "MD.VariableCell       false\n"
)

_INCLUDE_RE = re.compile(r"^\s*%include\s+(\S+)", re.I | re.M)
_LATTICE_DEF_RE = re.compile(r"^\s*(%block\s+(LatticeVectors|LatticeParameters)\b|LatticeConstant\b)", re.I | re.M)


def strain_folder_name(strain_pct):
    prefix = "m" if strain_pct < 0 else ""
    return f"vol_{prefix}{abs(strain_pct):.2f}"


def shortest_distance(structure, volume_factor=1.0):
    """Shortest interatomic distance (Ang) of `structure` (periodic), scaled by the linear factor."""
    try:
        from pymatgen.io.ase import AseAtomsAdaptor
        from ase.neighborlist import neighbor_list
        atoms = AseAtomsAdaptor.get_atoms(structure_io.to_pymatgen(structure))
        d = neighbor_list("d", atoms, 7.0)
        return float(d.min()) * volume_factor if len(d) else None
    except Exception:
        return None


def space_group_info(structure, symprec=0.01):
    """(symbol, number, crystal system) from spglib via pymatgen, or None."""
    try:
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
        sga = SpacegroupAnalyzer(structure_io.to_pymatgen(structure), symprec=symprec)
        return sga.get_space_group_symbol(), sga.get_space_group_number(), sga.get_crystal_system()
    except Exception:
        return None


def main():
    parser = argparse.ArgumentParser(
        description="Stage 1 of 2 (equation of state): writes one SIESTA folder per scanned volume. "
                     "Takes an ALREADY-RELAXED structure (--file, fractional coordinates) and your "
                     "calc.fdf template (--calc, default ./calc.fdf), and writes "
                     "'<output-dir>/vol_<pct>/' with the cell scaled isotropically (volume strain "
                     "from -strain-range to +strain-range), your template untouched except for a "
                     "'%include config_extra.fdf' at its top, the config_extra.fdf sidecar with the "
                     "forced directives (fixed cell, single point -- or ionic relaxation with "
                     "--relax-positions --, and the SCF k-grid when the template defines none or "
                     "you ask for one) and the pseudopotentials (-p). Also writes eos_scan.json "
                     "(the plan, read by stb-eosAnalysis). Does not run SIESTA: run each folder, "
                     "then use stb-eosAnalysis.",
        epilog="Example usage:\n"
               "  stb-eosInputs --file structure.fdf -p dojo\n"
               "  stb-eosInputs --file structure.fdf --calc calc.fdf --n-volumes 11 --strain-range 8.0 -p dojo\n"
               "  stb-eosInputs --file structure.fdf --kgrid-density 0.15 --relax-positions -o eos_runs\n",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("-f", "--file", required=True,
                        help="Input structure (.fdf), already relaxed. Atomic coordinates must be in "
                             "fractional format (same requirement as stb-strain). It is written, "
                             "scaled, under the same name into every folder, and your calc.fdf must "
                             "%%include it.")
    parser.add_argument("-c", "--calc", default=None, metavar="FDF",
                        help="SIESTA input template (default: calc.fdf). Copied into every folder with "
                             "'%%include config_extra.fdf' prepended; not edited otherwise. If it does not "
                             "exist and -c was not given, only the scaled structures are written (the "
                             "earlier behavior). If it is the same file as --file (a self-contained input), "
                             "the scaled copy of it is the run input.")
    parser.add_argument("-p", "--pseudo-dir", default="",
                        help="Pseudopotential source: a bundled bank (dojo, virtual_vault) or a folder. "
                             "The .psml/.psf of each species is copied into every folder.")
    parser.add_argument("--n-volumes", type=int, default=9,
                        help="Number of volumes to scan (default: 9). At least 5 for the 4-parameter fit "
                             "of stb-eosAnalysis; an odd number includes the input volume itself.")
    parser.add_argument("--strain-range", type=float, default=5.0,
                        help="Max isotropic volume strain, %% (default: 5.0), spanning "
                             "[-strain-range, +strain-range] -- the same meaning as stb-mleos's "
                             "--strain-range.")
    parser.add_argument("--relax-positions", action="store_true",
                        help="Let the ions relax at every volume (cell fixed): needed for crystals with "
                             "free internal coordinates. Default: single point (positions as given).")
    parser.add_argument("--relax-steps", type=int, default=100,
                        help="Max CG steps per volume with --relax-positions (default: 100).")
    parser.add_argument("--kgrid", type=int, nargs=3, default=None, metavar=("NX", "NY", "NZ"),
                        help="Force this SCF k-grid in every folder (config_extra.fdf), over the template's.")
    parser.add_argument("--kgrid-density", type=float, default=None, metavar="D",
                        help="Force an SCF k-grid of this density, 1/Ang (as stb-kgrid), over the template's. "
                             "Without --kgrid/--kgrid-density the template's own grid is kept; if it has "
                             f"none, one of density {DEFAULT_KGRID_DENSITY:g} is forced (SIESTA would sample "
                             "Gamma only).")
    parser.add_argument("-o", "--output-dir", default="eos_runs",
                        help="Directory for the 'vol_<pct>' folders (default: eos_runs).")
    parser.add_argument("--vacuum-gap", type=float, default=10.0,
                        help="Gap (Ang) used to detect vacuum-padded axes -- bulk only, a vacuum-padded "
                             "structure is rejected (default: 10.0).")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to <output-dir>/{REPORT_FILE}. Off by default.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-eosInputs {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.n_volumes < 5:
        parser.error("--n-volumes needs at least 5 points for a robust downstream EOS fit.")
    if not 0 < args.strain_range < 40:
        parser.error("--strain-range must be between 0 and 40 (%).")
    if args.relax_steps < 1:
        parser.error("--relax-steps must be >= 1.")

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("EQUATION-OF-STATE WORKFLOW -- STAGE 1: VOLUME FOLDERS", 'bold'))
    print("-" * 60)

    if not os.path.isfile(args.file):
        sys.exit(color_text(f"[ERROR] Input file '{args.file}' not found.", 'red'))
    explicit_calc = args.calc is not None
    args.calc = args.calc or "calc.fdf"
    calc_mode = "template"
    if not os.path.isfile(args.calc):
        if explicit_calc:
            sys.exit(color_text(f"[ERROR] Calc file '{args.calc}' not found.", 'red'))
        calc_mode = "structure-only"
    elif os.path.samefile(args.calc, args.file):
        calc_mode = "self-contained"
    pseudo_label = args.pseudo_dir
    if args.pseudo_dir:
        try:
            args.pseudo_dir = resolve_pseudo_source(args.pseudo_dir)
        except ValueError as e:
            sys.exit(color_text(f"[ERROR] {e}", 'red'))
    try:
        structure = structure_io.read_fdf(args.file)
    except Exception as e:
        sys.exit(color_text(f"[ERROR] Could not read '{args.file}': {e}", 'red'))

    species = structure_io.species_list(structure)
    counts = structure_io.atom_counts(structure)
    formula = " ".join(f"{sym}{counts[sym]}" for sym in species)
    n_atoms = len(structure.atoms)
    lattice = structure_io.lattice_only(structure)
    raw_vectors = structure_io.raw_lattice_vectors(structure)
    v0 = float(abs(np.linalg.det(lattice)))
    lengths = [float(np.linalg.norm(v)) for v in lattice]
    angles = [float(np.degrees(np.arccos(np.dot(lattice[(i + 1) % 3], lattice[(i + 2) % 3])
                                          / (lengths[(i + 1) % 3] * lengths[(i + 2) % 3])))) for i in range(3)]
    positions = np.array([pos for _, pos in structure.atoms])
    is_cartesian = structure.coord_format == 'cartesian'
    frac_coords = kspace.to_fractional(positions, structure.lattice, is_cartesian)
    vacuum_axes = kspace.detect_vacuum_axes(frac_coords, structure.lattice, args.vacuum_gap)
    if any(vacuum_axes):
        sys.exit(color_text(
            f"[ERROR] {kspace.dimensionality_label(vacuum_axes)}: this tool is bulk (3D periodic) only -- an "
            "equation of state assumes a well-defined, bounded cell volume, which a slab/wire/molecule "
            "doesn't have.", 'red'))

    os.makedirs(args.output_dir, exist_ok=True)
    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None
    warnings = []

    def warn(text):
        warnings.append(text)
        return color_text(f"[WARNING] {text}", 'yellow')

    calc_template = ""
    if calc_mode == "template":
        with open(args.calc) as f:
            calc_template = f.read()
    elif calc_mode == "self-contained":
        with open(args.file) as f:
            calc_template = f.read()
    template_values = fdf_directives.parse_directives(calc_template)
    template_kgrid = fdf_directives.parse_kgrid(calc_template)

    print_dual(color_text("===== STB-EOSINPUTS REPORT =====", 'magenta'), f_out)

    # ------------------------------------------------------------------ [0]
    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Tool            : stb-eosInputs {VERSION}", f_out)
    print_dual(f"Date            : {datetime.datetime.now():%Y-%m-%d %H:%M}", f_out)
    print_dual(f"Structure       : {args.file}", f_out)
    calc_desc = {"template": f"{args.calc}" + ("  (default)" if not explicit_calc else ""),
                 "self-contained": f"{args.calc} (the structure file is itself the input)",
                 "structure-only": "(none -- ./calc.fdf not found: only the scaled structures are written)"}[calc_mode]
    print_dual(f"Calc template   : {calc_desc}", f_out)
    print_dual(f"Pseudo source   : {(pseudo_label + '  (' + args.pseudo_dir + ')') if pseudo_label and pseudo_label != args.pseudo_dir else (args.pseudo_dir or '(none -- pseudopotentials not copied)')}", f_out)
    print_dual(f"Scan            : {args.n_volumes} volumes, +/-{args.strain_range:g} % in volume", f_out)
    print_dual(f"Ions            : {'relax at every volume (cell fixed), max ' + str(args.relax_steps) + ' CG steps' if args.relax_positions else 'fixed (single point)'}", f_out)
    print_dual(f"Output          : {args.output_dir}/", f_out)

    # ------------------------------------------------------------------ [1]
    print_section("[1] INPUT STRUCTURE", f_out)
    print_dual(f"Formula         : {formula}  ({n_atoms} atom(s), {len(species)} species)", f_out)
    print_dual(f"Lattice (Ang)   : a={lengths[0]:.4f}  b={lengths[1]:.4f}  c={lengths[2]:.4f}", f_out)
    print_dual(f"Angles (deg)    : alpha={angles[0]:.2f}  beta={angles[1]:.2f}  gamma={angles[2]:.2f}", f_out)
    print_dual(f"Volume          : {v0:.4f} Ang^3  ({v0 / n_atoms:.4f} Ang^3/atom)", f_out)
    try:
        density = structure_io.to_pymatgen(structure).density
        print_dual(f"Density         : {density:.4f} g/cm^3", f_out)
    except Exception:
        pass
    d_min = shortest_distance(structure)
    if d_min:
        print_dual(f"Shortest bond   : {d_min:.4f} Ang", f_out)
    print_dual(f"Detected        : {kspace.dimensionality_label(vacuum_axes)}", f_out)
    sg = space_group_info(structure)
    if sg:
        print_dual(f"Space group     : {sg[0]} (No. {sg[1]}), {sg[2]} system", f_out)
        if sg[2] != "cubic":
            print_dual(warn(
                f"The {sg[2]} cell is scaled isotropically, so the scan keeps its shape: the fit gives "
                "the bulk modulus along the constrained path (fixed c/a and angles), not the fully "
                "relaxed solid's. To compare with experiment, relax the shape first at the volume found "
                "(or run stb-elasticInputs for the full tensor)."), f_out)
    rows, total_e = [], 0.0
    for sym in species:
        zval, source = pseudo_valence(sym, args.pseudo_dir)
        if zval is not None:
            total_e += zval * counts[sym]
        rows.append(([sym, str(structure.species_meta.get(sym, {}).get('Z', '?')), str(counts[sym]),
                      f"{zval:g}" if zval is not None else "?",
                      f"{zval * counts[sym]:g}" if zval is not None else "?", source], None))
    print_table(["Species", "Z", "Atoms", "Valence e-", "Total e-", "Valence from"], rows, f_out)
    print_dual(f"Valence electrons : {total_e:g}", f_out)
    if is_cartesian:
        print_dual(warn("Atomic coordinates are Cartesian, not fractional -- scaling only the lattice "
                        "vectors would NOT scale the positions, giving a physically wrong (fixed-Cartesian) "
                        "scan. Convert to fractional coordinates first (e.g. stb-translate)."), f_out)

    # ------------------------------------------------------------------ [2]
    print_section("[2] LEVEL OF THEORY READ FROM THE TEMPLATE", f_out)
    forced_kgrid, kgrid_source = None, None
    if calc_mode == "structure-only":
        print_dual("No calc template: the folders get only the scaled structure, and the level of theory "
                   "is whatever you put around it.", f_out)
    else:
        theory_keys = [("xcfunctional", "XC functional"), ("xcauthors", "XC authors"),
                       ("paobasistype", "Basis type"), ("paobasissize", "Basis size"),
                       ("paoenergyshift", "PAO energy shift"), ("meshcutoff", "Mesh cutoff"),
                       ("spin", "Spin"), ("electronictemperature", "Electronic temperature"),
                       ("dmtolerance", "DM tolerance"), ("maxscfiterations", "Max SCF iterations"),
                       ("dftd3", "DFT-D3")]
        theory_rows = [([label, template_values[key]], None) for key, label in theory_keys if key in template_values]
        if template_kgrid:
            sp = reciprocal_spacings(lattice, template_kgrid)
            theory_rows.append((["SCF k-grid (kgrid.MonkhorstPack)",
                                 "x".join(str(v) for v in template_kgrid)
                                 + f"  (spacing {sp[0]:.3f} {sp[1]:.3f} {sp[2]:.3f} 1/Ang)"], None))
        elif "kgridcutoff" in template_values:
            theory_rows.append((["SCF k-grid (kgrid.Cutoff)", template_values["kgridcutoff"]], None))
        print_table(["Setting", "Value in the template"], theory_rows or [(["(none recognized)", "-"], None)], f_out)

        template_has_kgrid = bool(template_kgrid) or "kgridcutoff" in template_values
        if args.kgrid is not None:
            forced_kgrid, kgrid_source = tuple(args.kgrid), "explicit --kgrid"
        elif args.kgrid_density is not None or not template_has_kgrid:
            dens = args.kgrid_density if args.kgrid_density is not None else DEFAULT_KGRID_DENSITY
            try:
                forced_kgrid = tuple(kspace.compute_monkhorts(lattice[0], lattice[1], lattice[2], dens))
            except ValueError as e:
                sys.exit(color_text(f"[ERROR] {e}", 'red'))
            kgrid_source = (f"--kgrid-density {dens:g} 1/Ang" if args.kgrid_density is not None else
                            f"auto-suggested, density {dens:g} 1/Ang (the template defines no k-grid)")
        if forced_kgrid is not None:
            sp = reciprocal_spacings(lattice, forced_kgrid)
            print_dual(f"SCF k-grid forced : {'x'.join(str(v) for v in forced_kgrid)}  (spacing "
                       f"{sp[0]:.3f} {sp[1]:.3f} {sp[2]:.3f} 1/Ang) -- {kgrid_source}; written into "
                       f"{CONFIG_EXTRA_FILE}" + (f", overriding the template's {'x'.join(str(v) for v in template_kgrid)}"
                                                 if template_kgrid else ""), f_out)
        else:
            print_dual("SCF k-grid forced : no -- the template's own grid is kept (--kgrid / --kgrid-density override it)", f_out)

        includes = [os.path.basename(m) for m in _INCLUDE_RE.findall(calc_template)]
        if calc_mode == "template":
            if os.path.basename(args.file) not in includes:
                print_dual(warn(f"'{args.calc}' does not %include '{os.path.basename(args.file)}' (it includes: "
                                f"{', '.join(includes) if includes else 'nothing'}): SIESTA would NOT use the scaled "
                                "cell written into the folders. Add '%include " + os.path.basename(args.file)
                                + "' to the template."), f_out)
            if _LATTICE_DEF_RE.search(calc_template):
                print_dual(warn(f"'{args.calc}' defines a lattice of its own (LatticeVectors/LatticeParameters/"
                                "LatticeConstant) next to the included structure: whichever SIESTA reads first "
                                "wins, and the scan may silently use the unscaled cell. Keep the cell only in "
                                f"'{args.file}'."), f_out)
        shift = template_values.get("paoenergyshift", "")
        if shift:
            print_dual(color_text(
                f"[NOTE] PAO.EnergyShift {shift}: V0 and B0 depend on the orbital cutoff radii; converge them "
                "(smaller shift = longer orbitals) together with the k-grid and Mesh.Cutoff before quoting a "
                "bulk modulus.", 'cyan'), f_out)
        if "paobasissize" in template_values and template_values["paobasissize"].upper() in ("SZ", "DZ", "SZP"):
            print_dual(warn(f"Basis {template_values['paobasissize']} is small for an equation of state: "
                            "DZP or better is the usual starting point."), f_out)
        if "meshcutoff" not in template_values:
            print_dual(warn("The template defines no Mesh.Cutoff: SIESTA's small default is used, and E(V) "
                            "picks up egg-box noise of the real-space grid."), f_out)

    # ------------------------------------------------------------------ [3]
    print_section("[3] SCAN DEFINITION", f_out)
    strains = np.linspace(-abs(args.strain_range), abs(args.strain_range), args.n_volumes)
    print_dual(f"Volume strain   : {-args.strain_range:+.2f} % .. {args.strain_range:+.2f} %, {args.n_volumes} points "
               f"(step {strains[1] - strains[0]:.3f} %)", f_out)
    print_dual("Cell scaling    : L = L0 * (1 + strain/100)^(1/3) on every lattice vector (shape kept); "
               "fractional positions unchanged", f_out)
    if args.n_volumes % 2 == 0:
        print_dual(color_text("[NOTE] An even number of volumes has no point at the input volume itself; use an "
                              "odd number to keep V0_input in the scan.", 'cyan'), f_out)
    plan = []
    rows = []
    for s in strains:
        factor = (1.0 + s / 100.0) ** (1.0 / 3.0)
        vol = v0 * (1.0 + s / 100.0)
        plan.append({"folder": strain_folder_name(s), "strain_pct": float(s), "linear_factor": float(factor),
                     "volume": float(vol), "volume_per_atom": float(vol / n_atoms),
                     "lengths": [float(x * factor) for x in lengths],
                     "d_min": (d_min * factor if d_min else None)})
        rows.append(([strain_folder_name(s) + "/", f"{s:+.2f}", f"{vol:.4f}", f"{vol / n_atoms:.4f}", f"{factor:.5f}",
                      " ".join(f"{x * factor:.4f}" for x in lengths),
                      f"{d_min * factor:.4f}" if d_min else "-"], None))
    print_table(["Folder", "dV/V0 (%)", "V (Ang^3)", "V/atom", "L/L0", "a b c (Ang)", "d_min (Ang)"], rows, f_out)
    if args.strain_range < 2.0:
        print_dual(warn(f"--strain-range {args.strain_range:g} % is very narrow: the curvature of E(V) is then "
                        "comparable to the numerical noise of the SCF. 4-8 % is the usual window."), f_out)
    elif args.strain_range > 15.0:
        print_dual(warn(f"--strain-range {args.strain_range:g} % reaches far from equilibrium; the EOS forms "
                        "fit the whole range worse and the B0' they return is less reliable. 4-10 % is typical."), f_out)

    # ------------------------------------------------------------------ [4]
    print_section("[4] DIRECTIVES FORCED THROUGH config_extra.fdf (your template is not edited)", f_out)
    extra_text = ""
    if calc_mode != "structure-only":
        header = ("# Auto-generated by stb-eosInputs -- %include'd at the very top of calc.fdf. SIESTA's fdf reader is\n"
                  "# first-occurrence-wins, so these values override any conflicting setting in the template.\n\n")
        kgrid_block = ""
        if forced_kgrid is not None:
            kgrid_block = (f"# SCF k-grid ({kgrid_source}) -- the SAME grid in every folder.\n"
                           f"kgrid.MonkhorstPack   [{forced_kgrid[0]}  {forced_kgrid[1]}  {forced_kgrid[2]}]\n\n")
        md_block = (_RELAX_EXTRA.format(steps=args.relax_steps) if args.relax_positions else _SINGLE_POINT_EXTRA)
        extra_text = header + kgrid_block + md_block

        def template_has(*keys):
            return next((f"{template_values[k]}" for k in keys if k in template_values), "(absent)")
        forced_rows = [
            (["kgrid.MonkhorstPack",
              "x".join(str(v) for v in forced_kgrid) if forced_kgrid else "(not forced)",
              ("x".join(str(v) for v in template_kgrid) if template_kgrid else template_has("kgridcutoff")),
              kgrid_source if forced_kgrid else "the template's own grid"], None),
            (["MD.TypeOfRun", "CG", template_has("mdtypeofrun"), "no molecular dynamics"], None),
            (["MD.Steps", str(args.relax_steps) if args.relax_positions else "0", template_has("mdsteps"),
              "ionic relaxation budget" if args.relax_positions else "single point, no relaxation"], None),
            (["MD.NumCGsteps", str(args.relax_steps) if args.relax_positions else "0", template_has("mdnumcgsteps"),
              "older spelling of the same budget, forced for templates that still use it"], None),
            (["MD.VariableCell", "false", template_has("mdvariablecell"),
              "the cell is the scan variable: it must NEVER relax"], None),
        ]
        if args.relax_positions:
            forced_rows.append((["MD.MaxForceTol", RELAX_FORCE_TOL, template_has("mdmaxforcetol"),
                                 "ionic convergence criterion"], None))
        print_dual(f"Every folder gets its own {CONFIG_EXTRA_FILE}, and calc.fdf starts with '%include "
                   f"{CONFIG_EXTRA_FILE}' followed by your template untouched.", f_out)
        print_table(["Directive", "Forced value", "Template had", "Effect"], forced_rows, f_out)
        print_dual(color_text("Everything else (basis, mesh cutoff, XC, SCF) is whatever your template says; the same "
                              "settings in every folder is what makes E(V) a single curve.", 'cyan'), f_out)
        for key, label in (("mdvariablecell", "MD.VariableCell"),):
            if template_values.get(key, "").lower().strip(".") in ("true", "t", "yes"):
                print_dual(color_text(f"[NOTE] The template had {label} true (a variable-cell relaxation); it is "
                                      "overridden by the forced false above.", 'cyan'), f_out)
    else:
        print_dual("(no calc template: nothing is forced)", f_out)

    # ------------------------------------------------------------------ [5]
    print_section("[5] FOLDERS WRITTEN", f_out)
    structure_name = os.path.basename(args.file)
    run_input = "calc.fdf" if calc_mode == "template" else structure_name
    folder_rows, pseudo_problems = [], []
    for item, s in zip(plan, strains):
        folder = os.path.join(args.output_dir, item["folder"])
        os.makedirs(folder, exist_ok=True)
        factor = item["linear_factor"]
        scaled_path = os.path.join(folder, structure_name)
        structure_io.rewrite_fdf_lattice(args.file, raw_vectors * factor, scaled_path)
        files = [structure_name]
        if calc_mode == "template":
            with open(os.path.join(folder, "calc.fdf"), "w") as f:
                f.write(structure_io.prepend_include(calc_template, CONFIG_EXTRA_FILE))
            with open(os.path.join(folder, CONFIG_EXTRA_FILE), "w") as f:
                f.write(extra_text)
            files += ["calc.fdf", CONFIG_EXTRA_FILE]
        elif calc_mode == "self-contained":
            with open(scaled_path) as f:
                scaled_text = f.read()
            with open(scaled_path, "w") as f:
                f.write(structure_io.prepend_include(scaled_text, CONFIG_EXTRA_FILE))
            with open(os.path.join(folder, CONFIG_EXTRA_FILE), "w") as f:
                f.write(extra_text)
            files += [CONFIG_EXTRA_FILE]
        for sym in species:
            copy_pseudo(args.pseudo_dir, sym, folder)
            for ext in ("psml", "psf"):
                if os.path.exists(os.path.join(folder, f"{sym}.{ext}")):
                    files.append(f"{sym}.{ext}")
                    break
            else:
                if args.pseudo_dir:
                    pseudo_problems.append(f"{item['folder']}: {sym}")
        folder_rows.append(([item["folder"] + "/", f"{s:+.2f}", ", ".join(files)], None))
    print_dual(f"Under {args.output_dir}/:", f_out)
    print_table(["Folder", "dV/V0 (%)", "Files"], folder_rows, f_out)
    for problem in pseudo_problems:
        print_dual(warn(f"no pseudopotential file for {problem} under '{args.pseudo_dir}': add it to the folder "
                        "before running SIESTA."), f_out)
    if not args.pseudo_dir and calc_mode != "structure-only":
        print_dual(warn("no -p/--pseudo-dir: the pseudopotentials were not copied; SIESTA needs <species>.psml/.psf "
                        "in every folder."), f_out)

    manifest = {
        "format": "stb-eos-scan", "version": 1, "tool": f"stb-eosInputs {VERSION}",
        "created": datetime.datetime.now().isoformat(timespec="seconds"),
        "structure": {"file": args.file, "formula": formula, "n_atoms": n_atoms, "species": species,
                      "volume": v0, "lengths": lengths, "angles": angles,
                      "space_group": (sg[0] if sg else None), "crystal_system": (sg[2] if sg else None),
                      "d_min": d_min, "latticeconstant": structure.lattice_constant if hasattr(structure, "lattice_constant") else None},
        "calc": {"mode": calc_mode, "file": (args.calc if calc_mode != "structure-only" else None),
                 "run_input": run_input, "output_expected": "calc.out",
                 "kgrid_forced": list(forced_kgrid) if forced_kgrid else None,
                 "kgrid_template": list(template_kgrid) if template_kgrid else None,
                 "relax_positions": bool(args.relax_positions), "relax_steps": args.relax_steps,
                 "settings": {k: v for k, v in template_values.items()
                              if k in ("xcfunctional", "xcauthors", "paobasistype", "paobasissize", "paoenergyshift",
                                       "meshcutoff", "spin", "electronictemperature", "dmtolerance", "dftd3")}},
        "scan": {"n_volumes": args.n_volumes, "strain_range_pct": args.strain_range, "points": plan},
        "pseudo_source": pseudo_label or None,
    }
    with open(os.path.join(args.output_dir, MANIFEST_FILE), "w") as f:
        json.dump(manifest, f, indent=1)
        f.write("\n")

    # ------------------------------------------------------------------ [6]
    print_section("[6] SUMMARY & NEXT STEPS", f_out)
    print_dual(f"{len(plan)} volume folder(s) written under '{args.output_dir}' ({run_input} is the SIESTA input).", f_out)
    print_dual(f"Scan record          : {os.path.join(args.output_dir, MANIFEST_FILE)} (read by stb-eosAnalysis)", f_out)
    cost = "single-point SCF" if not args.relax_positions else f"SCF + up to {args.relax_steps} CG steps"
    print_dual(f"Cost                 : {len(plan)} independent runs ({cost} each)", f_out)
    if report_path:
        print_dual(f"Report               : {report_path}", f_out)
    print_dual(color_text("\nNext steps:", 'yellow'), f_out)
    print_dual("  1. Run SIESTA in every folder, e.g.:", f_out)
    print_dual(f"       for d in {args.output_dir}/vol_*; do (cd $d && mpirun -np 4 siesta {run_input} > calc.out); done", f_out)
    print_dual("  2. Check that every run finished (SCF converged), then:", f_out)
    print_dual(f"       stb-eosAnalysis --dir {args.output_dir} --eos all --show", f_out)
    print_dual(color_text(
        "\nConvergence checklist (what changes V0 and B0): PAO.EnergyShift and basis size (largest effect "
        "on B0), SCF k-grid, Mesh.Cutoff, electronic temperature, pseudopotential. Re-run this stage with "
        "another value and compare the Stage 2 results.", 'cyan'), f_out)
    if warnings:
        print_dual(color_text(f"\n{len(warnings)} warning(s) above -- read them before running SIESTA:", 'yellow'), f_out)
        for w in warnings:
            print_dual(f"  - {w}", f_out)

    if f_out:
        f_out.close()

    print("\n[INFO] Complete job!")
    print("\n" + "-" * 60)
    print(color_text("Volume folders ready for Stage 2 (stb-eosAnalysis).\n", 'bold'))


if __name__ == "__main__":
    main()
