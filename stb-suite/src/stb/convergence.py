#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""Stage 1 of 2: writes one ready-to-run SIESTA folder (structure + calc.fdf
+ config_extra.fdf + pseudopotentials) per swept value of one or more
convergence-test parameters (Mesh.CutOff, PAO.EnergyShift, k-grid density),
one subfolder per parameter under -o/--output-dir. Doesn't run SIESTA itself
-- run each folder yourself, then use stb-convergenceAnalysis.

Every generated folder fully relaxes the structure (both atomic positions
AND the cell, MD.VariableCell true) at each swept value, rather than running
a fixed-geometry single-point SCF. This is deliberate, and specific to why a
Mesh.CutOff/PAO.EnergyShift/k-grid convergence test matters for SIESTA in
particular: SIESTA evaluates the Hartree/XC potential on a discrete
real-space grid (set by Mesh.CutOff), so forces and stresses carry a small,
grid-position-dependent numerical noise (the "eggbox effect") that plane
-wave codes don't have. At a coarse Mesh.CutOff (or an incomplete PAO basis,
or too sparse a k-mesh) this noise is large enough to bias where a variable
-cell relaxation actually lands -- the SCF total energy at some ARBITRARY
fixed geometry can look fine even while the relaxed lattice parameters/bond
lengths are still drifting with the parameter. Relaxing the full structure
at each swept value is what actually answers "is this setting good enough
to trust the geometry SIESTA gives me", which is normally the point of
running the calculation in the first place -- not just "is the energy at
one arbitrary geometry converged".

'%include config_extra.fdf' is PREPENDED at the very top of the generated
calc.fdf (SIESTA's fdf reader is first-occurrence-wins for duplicate
labels), forcing CG relaxation (MD.Steps = --relax-steps,
MD.VariableCell true) regardless of what --calc's own template says. The
swept parameter's own value is ALSO written into config_extra.fdf (not
regex-substituted into --calc's own body) for the exact same reason -- it's
forced active via the same first-occurrence-wins mechanism, so it doesn't
need to already exist in --calc at all.

Output layout: <output-dir>/<parameter>/convergence_<parameter>_<value>/ --
one subfolder per swept parameter (so 'all' cleanly separates its 3
independent sweeps), same nested-plus-repeated-name convention
stb-elasticInputs already uses for its own <output-dir>/<direction>/
strain_<direction>_<pct>/ folders.

Per-parameter range flags (--meshcutoff-range/--energyshift-range/
--kgrid-range) let each selected parameter get its own custom sweep range
in the SAME invocation -- deliberately not a single shared --min/--max/
--step (an earlier version of this tool had that, and it could only ever
apply to exactly one parameter at a time, forcing multiple separate
invocations -- and multiple separate reports -- whenever more than one
parameter needed a non-default range). Any parameter without its own
-range flag uses its suggested default (see [2] of the report); every
selected parameter, customized or not, is always shown together in the
same single report.

KNOWN, DELIBERATE FOLLOW-UP: stb-convergenceAnalysis (Stage 2) does not yet
understand this per-parameter subfolder layout -- it still expects the old
flat <output-dir>/convergence_<parameter>_<value>/ layout (see its own
_FOLDER_RE + a non-recursive os.listdir(args.dir)). This is Stage 1 only;
Stage 2 gets its own follow-up update separately (same sequencing as this
session's phonons 5.4.1-then-5.4.2 work) -- flagged explicitly in [6] below
rather than left silently broken.
"""

from stb import __version__ as VERSION

import os
import re
import sys
import shutil
import argparse
import numpy as np
from stb.core import structure_io, kspace
from stb.core.cli import color_text, show_intro, print_dual, print_section, print_table
from stb.core.pseudopotentials import BANKS, resolve_pseudo_source, get_required_pseudos, copy_pseudo
from stb.core.structure_io import read_md_state, prepend_include

REPORT_FILE = "convergence_stage1.txt"
EXTRA_FDF_FILE = "config_extra.fdf"

# NOTE: substitute_numeric_tag/substitute_kgrid_tag/build_values/
# TAG_PATTERNS are imported directly by convergence_analysis.py
# (`from stb.convergence import substitute_numeric_tag, substitute_kgrid_tag`)
# -- kept for its own --apply feature (writing a converged value back into a
# real calc.fdf). Stage 1 itself no longer uses substitute_numeric_tag/
# substitute_kgrid_tag (see build_parameter_directive below) -- keep their
# names/signatures/behavior stable regardless.
TAG_PATTERNS = {
    "meshcutoff": re.compile(r'(Mesh\.CutOff\s+)([\d.]+)(\s*Ry)', re.IGNORECASE),
    "energyshift": re.compile(r'(PAO\.EnergyShift\s+)([\d.]+)(\s*Ry)', re.IGNORECASE),
}

# Suggested sweep ranges per parameter -- standard SIESTA convergence-test
# scales, not auto-derived from the structure itself. 'tighter' records
# which direction is MORE converged: meshcutoff is the opposite direction
# from energyshift/kgrid (a real, physical fact -- see
# convergence_analysis.py's own INCREASING_IMPROVES_ACCURACY, which encodes
# the same thing for its own plateau scan). Heavier elements/higher-accuracy
# work may need a wider meshcutoff range than the default; energyshift's
# textbook range is often quoted as a handful of specific values (0.001,
# 0.005, 0.01, 0.02, 0.05 Ry) rather than evenly spaced -- this tool's sweep
# is always linear (min/max/step), so the default below is an evenly spaced
# approximation of that same range, not the exact textbook set.
PARAMETER_DEFAULTS = {
    "meshcutoff":  {"min": 100.0, "max": 400.0, "step": 50.0, "unit": "Ry", "tighter": "larger"},
    "energyshift": {"min": 0.001, "max": 0.05, "step": 0.01, "unit": "Ry", "tighter": "smaller"},
    "kgrid":       {"min": 0.05, "max": 0.30, "step": 0.05, "unit": "1/Ang", "tighter": "smaller"},
}
PARAMETER_ORDER = ["meshcutoff", "energyshift", "kgrid"]
# PAO.EnergyShift acts logarithmically: almost all of the change sits between the
# smallest values, and a linear grid (0.001, 0.011, 0.021 ...) has no point in the
# range that is actually used. Default since 2.0.1: this geometric set.
ENERGYSHIFT_VALUES = [0.02, 0.01, 0.005, 0.002, 0.001, 0.0005]


def substitute_numeric_tag(calc_text, parameter, value):
    """Replaces the numeric value of a 'Tag   NUMBER  Ry' line, keeping the
    tag name and unit untouched. Raises ValueError if the tag isn't found --
    a silently-unmodified calc.fdf would make every run identical.
    """
    pattern = TAG_PATTERNS[parameter]
    new_text, count = pattern.subn(rf'\g<1>{value:.4f}\g<3>', calc_text)
    if count == 0:
        raise ValueError(
            f"Could not find a '{parameter}' tag to substitute in the calc.fdf template."
        )
    return new_text


def substitute_kgrid_tag(calc_text, structure, density, vacuum_gap):
    """Computes a dimension-aware Monkhorst-Pack grid at the given density
    (same logic as stb-kgrid) and substitutes it into 'kgrid.MonkhorstPack'.
    """
    positions = np.array([pos for _, pos in structure.atoms])
    is_cartesian = structure.coord_format == 'cartesian'
    frac_coords = kspace.to_fractional(positions, structure.lattice, is_cartesian)
    vacuum_axes = kspace.detect_vacuum_axes(frac_coords, structure.lattice, vacuum_gap)
    divisions = kspace.compute_monkhorts(
        structure.lattice[0], structure.lattice[1], structure.lattice[2], density, vacuum_axes
    )
    kgrid_str = f"kgrid.MonkhorstPack   [{divisions[0]}  {divisions[1]}  {divisions[2]}]"
    new_text, count = re.subn(r'kgrid\.MonkhorstPack\s+\[.*?\]', kgrid_str, calc_text, flags=re.IGNORECASE)
    if count == 0:
        raise ValueError(
            "Could not find a 'kgrid.MonkhorstPack' tag to substitute in the calc.fdf template."
        )
    return new_text, divisions


def build_values(min_v, max_v, step):
    n_steps = int(round((max_v - min_v) / step)) + 1
    return list(np.linspace(min_v, min_v + (n_steps - 1) * step, n_steps))


def build_parameter_directive(parameter, value, structure, vacuum_gap):
    """The swept parameter's own fdf directive line for THIS value, to be
    forced active via config_extra.fdf (not regex-substituted into --calc's
    body) -- so it doesn't need to already exist in --calc at all. Returns
    (directive_line, detail_string_for_the_report_table).
    """
    if parameter == "kgrid":
        positions = np.array([pos for _, pos in structure.atoms])
        is_cartesian = structure.coord_format == 'cartesian'
        frac_coords = kspace.to_fractional(positions, structure.lattice, is_cartesian)
        vacuum_axes = kspace.detect_vacuum_axes(frac_coords, structure.lattice, vacuum_gap)
        divisions = kspace.compute_monkhorts(
            structure.lattice[0], structure.lattice[1], structure.lattice[2], value, vacuum_axes)
        line = f"kgrid.MonkhorstPack   [{divisions[0]}  {divisions[1]}  {divisions[2]}]"
        detail = f"grid {divisions[0]} {divisions[1]} {divisions[2]}"
    elif parameter == "meshcutoff":
        line = f"Mesh.CutOff           {value:.4f}  Ry"
        detail = f"{value:.4f} Ry"
    elif parameter == "energyshift":
        line = f"PAO.EnergyShift     {value:.4f} Ry"
        detail = f"{value:.4f} Ry"
    else:
        raise ValueError(f"Unknown parameter '{parameter}'.")
    return line, detail


def run_basis_mode(args, parser):
    """--parameter basis: one stage of the basis-set convergence by equations of state (see
    core/basis_scan.py). Writes <output-dir>/basis/<stage>_<value>/vol_<pct>/ single points."""
    import json
    from stb.core import basis_scan as bs
    from stb.core import fdf_directives

    for path, what in ((args.structure, "Structure"), (args.calc, "Calc")):
        if not os.path.isfile(path):
            sys.exit(color_text(f"[ERROR] {what} file '{path}' not found.", 'red'))
    if args.pseudo_dir:
        try:
            args.pseudo_dir = resolve_pseudo_source(args.pseudo_dir)
        except ValueError as e:
            sys.exit(color_text(f"[ERROR] {e}", 'red'))
    with open(args.calc) as f:
        calc_text = f.read()
    template = fdf_directives.parse_directives(calc_text)
    stage = args.basis_stage

    def template_number(key):
        try:
            return float(template[key].split()[0])
        except (KeyError, ValueError, IndexError):
            return None

    template_mesh, template_shift = template_number("meshcutoff"), template_number("paoenergyshift")
    template_size = (template.get("paobasissize") or "DZP").split()[0]
    mesh_fixed = args.mesh_cutoff if args.mesh_cutoff is not None else template_mesh
    if stage == "mesh":
        shift_fixed = args.energyshift if args.energyshift is not None else bs.MESH_STAGE_SHIFT
        values = sorted(set(args.basis_values or bs.MESH_VALUES))
        size_fixed = template_size
    elif stage == "range":
        if mesh_fixed is None:
            parser.error("--basis-stage range needs the Mesh.Cutoff found by the mesh stage: pass --mesh-cutoff "
                         "(the template sets none).")
        values = sorted(set(args.basis_values or bs.SHIFT_VALUES), reverse=True)
        shift_fixed, size_fixed = None, template_size
    else:
        if mesh_fixed is None:
            parser.error("--basis-stage size needs --mesh-cutoff (the template sets no Mesh.Cutoff).")
        shift_fixed = args.energyshift if args.energyshift is not None else template_shift
        if shift_fixed is None:
            parser.error("--basis-stage size needs the PAO.EnergyShift found by the range stage: pass --energyshift.")
        upper = bs.next_size(template_size)
        values = args.basis_sizes or ([template_size, upper] if upper else None)
        if not values:
            parser.error(f"No basis size above '{template_size}' on the ladder {' < '.join(bs.SIZE_LADDER)}: "
                         "give the two sizes to compare with --basis-sizes.")
        size_fixed = None
    if stage != "size" and args.basis_values and any(v <= 0 for v in args.basis_values):
        parser.error("--basis-values must be positive.")

    if args.intro:
        show_intro(["Siesta ToolBox Suite", "A comprehensive toolkit for SIESTA DFT simulations",
                    f"Version {VERSION} | University of Brasilia - 2026", "Developed by Dr. Carlos M. O. Bastos"])
    base = os.path.join(args.output_dir, "basis")
    os.makedirs(base, exist_ok=True)
    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    def fail(message):
        print_dual(color_text(f"[FAIL] {message}", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_dual(color_text("===== STB-CONVERGENCE REPORT (BASIS SET BY EQUATIONS OF STATE) =====", 'magenta'), f_out)
    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Structure file    : {args.structure}", f_out)
    print_dual(f"Calc template     : {args.calc}", f_out)
    print_dual(f"Pseudo source     : {args.pseudo_dir or '(not given)'}", f_out)
    print_dual(f"Stage             : {stage} -- {bs.STAGE_TITLES[stage]}", f_out)
    print_dual(f"Output dir        : {base}", f_out)

    print_section("[1] INPUT STRUCTURE", f_out)
    try:
        structure = structure_io.read_fdf(args.structure)
    except (ValueError, FileNotFoundError) as e:
        fail(str(e))
    if structure.coord_format == 'cartesian':
        fail("Atomic coordinates are Cartesian: scaling the lattice vectors would not move the atoms. Convert to "
             "fractional coordinates first (stb-translate).")
    counts = structure_io.atom_counts(structure)
    species = sorted(sp for sp, n in counts.items() if n > 0)
    n_atoms = sum(counts.values())
    lattice = np.array(structure.lattice, dtype=float) * structure.lattice_constant
    positions = np.array([pos for _, pos in structure.atoms])
    vacuum_axes = kspace.detect_vacuum_axes(positions, lattice, args.vacuum_gap)
    periodic_axes = [i for i, vac in enumerate(vacuum_axes) if not vac]
    n_periodic = len(periodic_axes)
    if n_periodic == 0:
        fail("No periodic axis: this convergence scans the cell (volume, area or length).")
    volume0 = float(abs(np.linalg.det(lattice)))
    measure = {3: "volume", 2: "area", 1: "length"}[n_periodic]
    scan_range = args.scan_range if args.scan_range is not None else bs.SCAN_RANGE_PCT[n_periodic]
    strains = np.linspace(-scan_range, scan_range, args.scan_points)
    print_table(["Quantity", "Value"], [
        (["Composition", ", ".join(f"{sp}{n}" for sp, n in counts.items())], None),
        (["Atoms", str(n_atoms)], None),
        (["Dimensionality", kspace.dimensionality_label(vacuum_axes)], None),
        (["Cell volume", f"{volume0:.4f} Ang^3"], None),
    ], f_out)

    print_section("[2] THE THREE STAGES", f_out)
    print_dual(
        "A numerical-orbital basis has two independent axes, range (PAO.EnergyShift) and cardinality "
        "(PAO.BasisSize), and the total energy has no plateau in either. Each setting is therefore judged by "
        f"the equilibrium {measure} of an E({measure[0].upper()}) scan of single points, one stage at a time:", f_out)
    for name in bs.STAGES:
        print_dual(f"  {'->' if name == stage else '  '} {name:<6}: {bs.STAGE_TITLES[name]}", f_out)
    print_dual("Run this stage, analyse it (stb-convergenceAnalysis), then generate the next one with the value found.",
               f_out)
    caution = fdf_directives.filter_caution(calc_text)
    if caution:
        print_dual(color_text(f"[WARNING] {caution} With the filter on, this convergence tests the filtered basis.",
                              'yellow'), f_out)
    if re.search(r"^\s*%block\s+PAO\.?Basis\b", calc_text, re.M | re.I):
        print_dual(color_text(
            "[WARNING] The template has a %block PAO.Basis: SIESTA takes the orbitals of those species from the "
            "block and ignores PAO.EnergyShift / PAO.BasisSize for them, so the range and size stages would "
            "compare identical bases. Remove the block for this convergence.", 'yellow'), f_out)
    if stage == "mesh" and not args.energyshift:
        print_dual(f"PAO.EnergyShift is held at {shift_fixed:g} Ry in this stage (a tight value, so that the mesh "
                   "is tested with the long orbitals the next stage will use; --energyshift to change).", f_out)

    print_section("[3] SCANS OF THIS STAGE", f_out)
    print_dual(f"{measure.capitalize()} strain : -{scan_range:g} % .. +{scan_range:g} %, {args.scan_points} single points per "
               "scan" + ("" if n_periodic == 3 else f" (only the periodic vector(s) "
                         f"{', '.join('abc'[i] for i in periodic_axes)} are scaled; the vacuum is not touched)"), f_out)
    raw_vectors = np.asarray(structure_io.raw_lattice_vectors(structure), dtype=float)
    structure_basename, calc_basename = os.path.basename(args.structure), os.path.basename(args.calc)
    forced_calc = prepend_include(calc_text, EXTRA_FDF_FILE)
    rows, sub_scans = [], []
    for value in values:
        mesh = value if stage == "mesh" else mesh_fixed
        shift = value if stage == "range" else shift_fixed
        size = value if stage == "size" else size_fixed
        name = bs.sub_scan_name(stage, value)
        extra = ("# Auto-generated by stb-convergence --parameter basis: the settings of this scan and a single\n"
                 "# point at fixed cell, whatever the template says (first occurrence wins).\n"
                 f"PAO.BasisSize       {size}\n"
                 + (f"PAO.EnergyShift     {shift:g} Ry\n" if shift is not None else "")
                 + (f"Mesh.CutOff         {mesh:g} Ry\n" if mesh is not None else "")
                 + "\nMD.TypeOfRun        CG\nMD.Steps            0\nMD.VariableCell     false\n")
        for strain in strains:
            factor = (1.0 + strain / 100.0) ** (1.0 / n_periodic)
            axis_factors = np.array([factor if i in periodic_axes else 1.0 for i in range(3)])
            folder = os.path.join(base, name, bs.strain_folder(strain))
            os.makedirs(folder, exist_ok=True)
            structure_io.rewrite_fdf_lattice(args.structure, raw_vectors * axis_factors[:, None],
                                             os.path.join(folder, structure_basename))
            with open(os.path.join(folder, calc_basename), "w") as fh:
                fh.write(forced_calc)
            with open(os.path.join(folder, EXTRA_FDF_FILE), "w") as fh:
                fh.write(extra)
            for sym in species:
                copy_pseudo(args.pseudo_dir, sym, folder)
        sub_scans.append({"name": name, "value": value, "mesh_cutoff": mesh, "energyshift": shift, "size": size})
        rows.append(([name + "/", str(size), f"{shift:g} Ry" if shift is not None else "(template)",
                      f"{mesh:g} Ry" if mesh is not None else "(template)", str(len(strains))], None))
    print_table(["Scan", "PAO.BasisSize", "PAO.EnergyShift", "Mesh.Cutoff", "Folders"], rows, f_out)
    if stage != "mesh" and mesh_fixed is not None:
        print_dual(f"Mesh.Cutoff {mesh_fixed:g} Ry in every scan"
                   + (" (--mesh-cutoff)" if args.mesh_cutoff is not None else " (the template's; use the value the "
                      "mesh stage found)") + ".", f_out)

    manifest = {"format": "stb-basis-scan", "version": 1, "tool": f"stb-convergence {VERSION}", "stage": stage,
                "structure": args.structure, "calc": args.calc, "n_atoms": n_atoms, "volume": volume0,
                "n_periodic": n_periodic, "periodic_axes": periodic_axes, "measure": measure,
                "scan_range_pct": scan_range, "scan_points": args.scan_points, "sub_scans": sub_scans}
    with open(os.path.join(base, bs.manifest_name(stage)), "w") as fh:
        json.dump(manifest, fh, indent=1)

    print_section("[4] SUMMARY & NEXT STEPS", f_out)
    n_total = len(sub_scans) * len(strains)
    print_dual(f"{n_total} single-point folder(s) written under '{base}' ({len(sub_scans)} scan(s) x {len(strains)}).", f_out)
    print_dual(color_text("\nNext steps:", 'yellow'), f_out)
    print_dual(f"  1. Run SIESTA in every '{base}/{stage}_*/vol_*/' folder (single points).", f_out)
    print_dual(f"  2. stb-convergenceAnalysis --dir {args.output_dir}", f_out)
    following = {"mesh": "  3. Then the range stage at the mesh found:  stb-convergence --parameter basis "
                         "--basis-stage range --mesh-cutoff <Ry> ...",
                 "range": "  3. Then the size stage at the EnergyShift found:  stb-convergence --parameter basis "
                          "--basis-stage size --mesh-cutoff <Ry> --energyshift <Ry> ...",
                 "size": "  3. If the size stage passes, the basis of the template with the values found is the "
                         "converged one for this structure."}[stage]
    print_dual(following, f_out)
    print_dual(color_text(
        "\nThe verdict holds for THIS structure's equilibrium geometry. Energies that involve an isolated atom or "
        "molecule (cohesive, adsorption) need the counterpoise correction and their own test; in 2D the in-plane "
        "lattice parameter does not test how the orbitals decay into the vacuum (work function, adsorption, "
        "interlayer distance).", 'cyan'), f_out)
    if f_out:
        f_out.close()


def main():
    parser = argparse.ArgumentParser(
        description="Stage 1 of 2: writes one ready-to-run SIESTA folder per swept value of one "
                    "or more convergence-test parameters, one subfolder per parameter, fully "
                    "relaxing the structure (positions + cell) at each. Doesn't run SIESTA -- "
                    "run each folder yourself, then use stb-convergenceAnalysis.",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s -s structure.fdf -c calc.fdf -p dojo --parameter all\n"
               "  %(prog)s -s structure.fdf -c calc.fdf -p dojo --parameter all "
               "--meshcutoff-range 150 350 50\n"
               "  %(prog)s -s structure.fdf -c calc.fdf -p dojo --parameter kgrid "
               "--relax-steps 200 --save-report\n"
    )

    parser.add_argument("-s", "--structure", type=str, default="structure.fdf",
                        help="Input structure.fdf, copied as 'structure.fdf' into every "
                             "generated folder (default: structure.fdf).")
    parser.add_argument("-c", "--calc", type=str, default="calc.fdf",
                        help="calc.fdf template (default: calc.fdf) -- SCF/basis/MD settings, normally %%including "
                             "the structure file above by name. Copied into every generated "
                             "folder with one line added: '%%include config_extra.fdf', "
                             "prepended at the very top (see [3] of the report) -- this forces "
                             "full relaxation and the swept parameter's value, regardless of "
                             "what this template itself contains.")
    parser.add_argument("-p", "--pseudo-dir", default="",
                        help="Pseudopotentials source (optional): a bundled bank "
                             f"({', '.join(BANKS)}) or a folder path. If given, the required "
                             "pseudopotential for every species in the structure is copied into "
                             "every generated folder, so it's immediately ready for SIESTA.")
    parser.add_argument("--parameter", nargs='+', required=True,
                        help="Which tag(s) to sweep: meshcutoff, energyshift, kgrid, or 'all' "
                             "for all three (one independent sweep per parameter, each in its "
                             "own subfolder). Or 'basis' (alone): the basis-set convergence by "
                             "equations of state, see --basis-stage.")
    parser.add_argument("--basis-stage", choices=["mesh", "range", "size"], default="mesh",
                        help="--parameter basis: which stage to write (default: mesh). Each stage is a set of\n"
                             "E(V) scans of single points, judged by the equilibrium volume (2D: area):\n"
                             "  mesh : Mesh.Cutoff 350, 500, 700, 1000 Ry at PAO.EnergyShift 0.002 Ry\n"
                             "  range: PAO.EnergyShift 0.02 ... 0.0005 Ry at --mesh-cutoff\n"
                             "  size : the template's PAO.BasisSize and one step up (DZP -> TZP) at\n"
                             "         --mesh-cutoff and --energyshift\n"
                             "Run them in this order, each with the value the previous one found.")
    parser.add_argument("--basis-values", type=float, nargs='+', default=None, metavar="V",
                        help="--parameter basis: the Mesh.Cutoff values (stage mesh) or PAO.EnergyShift values "
                             "(stage range) to scan instead of the defaults.")
    parser.add_argument("--basis-sizes", nargs=2, default=None, metavar=("SIZE", "LARGER"),
                        help="--parameter basis --basis-stage size: the two PAO.BasisSize to compare (default: "
                             "the template's and the next one on SZ < DZ < DZP < TZP < TZ2P).")
    parser.add_argument("--mesh-cutoff", type=float, default=None, metavar="RY",
                        help="--parameter basis, stages range and size: Mesh.Cutoff of every scan (default: the "
                             "template's).")
    parser.add_argument("--energyshift", type=float, default=None, metavar="RY",
                        help="--parameter basis: PAO.EnergyShift held fixed in stage mesh (default 0.002) and "
                             "stage size (default: the template's).")
    parser.add_argument("--scan-points", type=int, default=9,
                        help="--parameter basis: single points per E(V) scan (default: 9, at least 5).")
    parser.add_argument("--scan-range", type=float, default=None, metavar="PCT",
                        help="--parameter basis: half-width of each scan in %% of the volume (default 12), of the "
                             "area for a 2D material (default 8) or of the length for a wire (default 4).")
    parser.add_argument("--meshcutoff-range", type=float, nargs=3, default=None,
                        metavar=("MIN", "MAX", "STEP"),
                        help="Custom Mesh.CutOff sweep range (Ry). Omit to use the suggested "
                             "default (see [2] of the report). Only meaningful if 'meshcutoff' "
                             "is in --parameter.")
    parser.add_argument("--energyshift-range", type=float, nargs=3, default=None,
                        metavar=("MIN", "MAX", "STEP"),
                        help="Custom PAO.EnergyShift sweep range (Ry) (see --meshcutoff-range). "
                             "Only meaningful if 'energyshift' is in --parameter.")
    parser.add_argument("--energyshift-values", type=float, nargs='+', default=None, metavar="RY",
                        help="Explicit PAO.EnergyShift values (Ry) instead of a linear range. Default when "
                             "neither this nor --energyshift-range is given: the geometric set "
                             f"{' '.join(f'{v:g}' for v in ENERGYSHIFT_VALUES)} (the parameter acts "
                             "logarithmically; versions up to 2.0.0 used 0.001 to 0.05 in steps of 0.01, "
                             "which has no point between 0.001 and 0.011).")
    parser.add_argument("--single-point", action="store_true",
                        help="Do not relax: every folder is a single point at the input geometry (MD.Steps 0, "
                             "fixed cell), so only the swept parameter changes. Default: full relaxation "
                             "(positions + cell), which mixes the parameter's effect with the relaxation's "
                             "noise but gives the relaxed volume to judge convergence by.")
    parser.add_argument("--kgrid-range", type=float, nargs=3, default=None,
                        metavar=("MIN", "MAX", "STEP"),
                        help="Custom k-grid density sweep range (1/Ang) (see "
                             "--meshcutoff-range). Only meaningful if 'kgrid' is in --parameter.")
    parser.add_argument("--relax-steps", type=int, default=100,
                        help="Number of relaxation steps (MD.Steps) forced in "
                             "every generated folder (default: 100). Both atomic positions and "
                             "the cell relax (MD.VariableCell true).")
    parser.add_argument("--max-stress-tol", type=float, default=0.02,
                        help="Relaxation mode: MD.MaxStressTol forced in every folder, in GPa (default: 0.02). "
                             "SIESTA's own default, 1 GPa, is ~1 %% of the volume of Si (B0 ~ 88 GPa) -- twice "
                             "the 0.5 %% the analysis judges by -- and stopped most relaxations at the input cell.")
    parser.add_argument("--max-force-tol", type=float, default=0.01,
                        help="Relaxation mode: MD.MaxForceTol forced in every folder, in eV/Ang (default: 0.01).")
    parser.add_argument("--vacuum-gap", type=float, default=10.0,
                        help="--parameter kgrid only: vacuum-axis detection threshold in Angstrom, "
                             "same meaning as stb-kgrid's --vacuum-gap (default: 10.0).")
    parser.add_argument("-o", "--output-dir", type=str, default="convergence_runs",
                        help="Top-level directory (default: convergence_runs). Run folders go "
                             "into its own '<output-dir>/<parameter>/convergence_<parameter>_"
                             "<value>/' subfolder, one tree per swept parameter.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to <output-dir>/{REPORT_FILE}. Off by default.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-convergence {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if "basis" in args.parameter:
        if len(args.parameter) > 1:
            parser.error("--parameter basis is a workflow of its own: give it alone.")
        if args.scan_points < 5:
            parser.error("--scan-points must be at least 5 (the equation of state has 4 parameters).")
        if args.scan_range is not None and not 0 < args.scan_range < 40:
            parser.error("--scan-range must be between 0 and 40 %.")
        run_basis_mode(args, parser)
        return

    # Expand 'all' / validate parameter names, de-duplicated, in a stable
    # (PARAMETER_ORDER) order regardless of how the user listed them.
    requested = []
    for p in args.parameter:
        if p == "all":
            for k in PARAMETER_ORDER:
                if k not in requested:
                    requested.append(k)
        elif p in PARAMETER_ORDER:
            if p not in requested:
                requested.append(p)
        else:
            parser.error(f"--parameter: unrecognized value '{p}' (choose from "
                         f"{', '.join(PARAMETER_ORDER)}, or 'all').")
    args.parameter = requested

    if args.relax_steps <= 0:
        parser.error("--relax-steps must be > 0.")

    range_flags = {
        "meshcutoff": args.meshcutoff_range,
        "energyshift": args.energyshift_range,
        "kgrid": args.kgrid_range,
    }
    for p, custom in range_flags.items():
        if custom is not None and p not in args.parameter:
            parser.error(f"--{p}-range was given but '{p}' is not in --parameter.")

    if args.energyshift_values is not None:
        if "energyshift" not in args.parameter:
            parser.error("--energyshift-values was given but 'energyshift' is not in --parameter.")
        if args.energyshift_range is not None:
            parser.error("--energyshift-values and --energyshift-range are two ways to say the same thing: give one.")
        if any(v <= 0 for v in args.energyshift_values):
            parser.error("--energyshift-values must be positive.")
    ranges = {}
    for p in args.parameter:
        custom = range_flags[p]
        if p == "energyshift" and custom is None:
            given = args.energyshift_values is not None
            ranges[p] = {"values": sorted(set(args.energyshift_values if given else ENERGYSHIFT_VALUES), reverse=True),
                         "source": "custom (--energyshift-values)" if given else "suggested default (geometric)"}
        elif custom is not None:
            mn, mx, st = custom
            if st <= 0:
                parser.error(f"--{p}-range: step must be > 0.")
            if mx <= mn:
                parser.error(f"--{p}-range: max must be greater than min.")
            ranges[p] = {"min": mn, "max": mx, "step": st, "source": f"custom (--{p}-range)"}
        else:
            d = PARAMETER_DEFAULTS[p]
            ranges[p] = {"min": d["min"], "max": d["max"], "step": d["step"],
                         "source": "suggested default"}

    if not os.path.isfile(args.structure):
        sys.exit(color_text(f"[ERROR] Structure file '{args.structure}' not found.", 'red'))
    if not os.path.isfile(args.calc):
        sys.exit(color_text(f"[ERROR] Calc file '{args.calc}' not found.", 'red'))
    if args.pseudo_dir:
        try:
            args.pseudo_dir = resolve_pseudo_source(args.pseudo_dir)
        except ValueError as e:
            sys.exit(color_text(f"[ERROR] {e}", 'red'))

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    os.makedirs(args.output_dir, exist_ok=True)
    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    def fail(message):
        print_dual(color_text(f"[FAIL] {message}", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_dual(color_text(
        "===== STB-CONVERGENCE STAGE 1 REPORT (CONVERGENCE TEST PREP) =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Structure file    : {args.structure}", f_out)
    print_dual(f"Calc template     : {args.calc}", f_out)
    print_dual(f"Pseudo source     : {args.pseudo_dir or '(not given)'}", f_out)
    print_dual(f"Parameter(s)      : {', '.join(args.parameter)}", f_out)
    print_dual(f"Geometry          : " + ("single point at the input geometry (--single-point)" if args.single_point
                                          else f"full relaxation, positions + cell, up to {args.relax_steps} steps"), f_out)
    print_dual(f"Output dir        : {args.output_dir}", f_out)
    print_dual(f"Save report       : {'yes' if args.save_report else 'no'}", f_out)
    if report_path:
        print_dual(f"Report file       : {report_path}", f_out)

    print_section("[1] INPUT STRUCTURE", f_out)
    try:
        structure = structure_io.read_fdf(args.structure)
    except (ValueError, FileNotFoundError) as e:
        fail(str(e))

    counts = structure_io.atom_counts(structure)
    species = sorted(sp for sp, n in counts.items() if n > 0)
    composition = ", ".join(f"{sp}{n}" for sp, n in counts.items())
    print_table(["Quantity", "Value"], [
        (["Composition", composition], None),
        (["Total atoms", str(sum(counts.values()))], None),
        (["Coordinate format", structure.coord_format], None),
    ], f_out)

    print_section("[2] PARAMETER SELECTION", f_out)
    rows = []
    for p in args.parameter:
        r = ranges[p]
        d = PARAMETER_DEFAULTS[p]
        rows.append(([p, r["source"],
                      (", ".join(f"{v:g}" for v in r["values"]) + f" ({d['unit']})") if "values" in r else
                      f"{r['min']:g} to {r['max']:g}, step {r['step']:g} ({d['unit']})",
                      f"{d['tighter']} = more converged"], None))
    print_table(["Parameter", "Range source", "Sweep range", "Convergence direction"], rows, f_out)
    if "energyshift" in args.parameter:
        print_dual(color_text(
            "[NOTE] PAO.EnergyShift: the TOTAL energy has no practical plateau in this parameter. Bases of "
            "different range are not nested, and the total energy still moves by 0.01-1.5 eV/atom between "
            "0.011 and 0.001 Ry (Si 13 meV, graphene 23 meV, NaCl 0.52 eV, Ag 1.28 eV, GaAs 1.46 eV). Judge "
            "this sweep by a physical quantity -- the relaxed volume / lattice parameter (reported by "
            "stb-convergenceAnalysis), or better V0 and B0 from an E(V) scan at each value (stb-eosInputs "
            "with a template per value) -- and an energy DIFFERENCE where one matters.", 'cyan'), f_out)

    print_section("[3] RELAXATION & PARAMETER ENFORCEMENT", f_out)
    with open(args.calc) as f:
        original_calc_text = f.read()
    structure_basename = os.path.basename(args.structure)
    calc_basename = os.path.basename(args.calc)
    if "%include" not in original_calc_text or structure_basename not in original_calc_text:
        print_dual(color_text(
            f"[NOTE] '{calc_basename}' may not %include '{structure_basename}' -- add it if "
            "needed so SIESTA picks up the structure in each generated folder.", 'yellow'), f_out)

    before = read_md_state(original_calc_text)
    steps_label = f"{before['steps_key']}={before['steps']}" if before['steps_key'] else "(absent)"
    print_dual(f"Calc template (before forcing): MD.TypeOfRun={before['typeofrun'] or '(absent)'}  "
               f"Steps: {steps_label}  MD.VariableCell={before['variablecell'] or '(absent)'}", f_out)
    if args.single_point:
        print_dual(f"Forced via '{EXTRA_FDF_FILE}' (prepended -- always wins) in every folder -- a single point "
                   "at the input geometry, so only the swept parameter changes between folders:", f_out)
        forced_rows = [(["MD.TypeOfRun", "CG"], None), (["MD.Steps", "0 (--single-point)"], None),
                       (["MD.VariableCell", "false"], None)]
    else:
        print_dual(f"Forced via '{EXTRA_FDF_FILE}' (prepended -- always wins) in every folder -- the "
                   "full structure (positions + cell) relaxes at each swept value:", f_out)
        forced_rows = [(["MD.TypeOfRun", "CG"], None), (["MD.Steps", str(args.relax_steps)], None),
                       (["MD.VariableCell", "true (full relaxation: positions + cell)"], None),
                       (["MD.MaxStressTol", f"{args.max_stress_tol:g} GPa (SIESTA's default: 1 GPa)"], None),
                       (["MD.MaxForceTol", f"{args.max_force_tol:g} eV/Ang"], None)]
    print_table(["Directive", "Forced value"], forced_rows + [
        (["Swept parameter", f"{', '.join(args.parameter)} -- one value per folder, see [5]"], None),
    ], f_out)
    forced_calc_text = prepend_include(original_calc_text, EXTRA_FDF_FILE)
    from stb.core import fdf_directives
    caution = fdf_directives.filter_caution(original_calc_text)
    if caution:
        print_dual(color_text(f"[WARNING] {caution}", 'yellow'), f_out)
    d3_value = fdf_directives.parse_directives(original_calc_text).get("dftd3", "")
    if d3_value.strip().lower().strip(".") in ("true", "t", "yes", "y"):
        print_dual(color_text(
            "[WARNING] The template turns DFT-D3 on and every folder relaxes the CELL. SIESTA 5.4.2 writes "
            "the D3 contribution to the stress tensor with the wrong sign (energies and forces are right; "
            "measured on Ag and Cu: P(D3) - P(PBE) = +4.75 GPa against -dE_D3/dV = -4.75 GPa), so with that "
            "version the cell relaxes to the wrong volume -- too large -- and the convergence is judged at a "
            "geometry that is not the minimum. Check the SIESTA version; with 5.4.2, converge the parameter "
            "without D3 or get the volume from the energies (stb-eosInputs / stb-eosAnalysis).",
            'yellow'), f_out)

    print_section("[4] PSEUDOPOTENTIALS", f_out)
    if args.pseudo_dir:
        found, missing = get_required_pseudos(species, args.pseudo_dir)
        print_dual(f"Source            : {args.pseudo_dir}", f_out)
        print_table(["Species", "Status"], [
            ([sp, "MISSING" if sp in missing else "found"], 'yellow' if sp in missing else None)
            for sp in species
        ], f_out)
        if missing:
            print_dual(color_text(
                f"[WARNING] Missing pseudopotential(s) for: {', '.join(missing)} -- these will "
                "need to be added manually to every generated folder.", 'yellow'), f_out)
        else:
            print_dual(color_text(
                "[OK] All required pseudopotentials found -- will be copied into every "
                "generated folder.", 'green'), f_out)
    else:
        print_dual("Not given (pass -p/--pseudo-dir -- a bundled bank or a folder path -- to "
                   "copy the required pseudopotential for every species into every generated "
                   "folder). Pseudopotentials will need to be added manually.", f_out)

    print_section("[5] GENERATED CONVERGENCE FOLDERS", f_out)
    print_dual(f"Under '{args.output_dir}/<parameter>/':", f_out)
    folder_rows = []
    per_parameter_count = {}
    for p in args.parameter:
        r = ranges[p]
        values = r["values"] if "values" in r else build_values(r["min"], r["max"], r["step"])
        for value in values:
            folder_name = f"convergence_{p}_{value:.4f}"
            folder = os.path.join(args.output_dir, p, folder_name)
            os.makedirs(folder, exist_ok=True)

            try:
                directive_line, detail = build_parameter_directive(p, value, structure, args.vacuum_gap)
            except ValueError as e:
                fail(str(e))

            if args.single_point:
                extra_fdf_text = (
                    "# Auto-generated by stb-convergence (--single-point) -- forces this folder's own swept\n"
                    "# parameter value and a single point at the input geometry, regardless of --calc's settings.\n"
                    f"{directive_line}\n"
                    "\n"
                    "MD.TypeOfRun       CG\n"
                    "MD.Steps           0\n"
                    "MD.VariableCell    false\n"
                )
            else:
                extra_fdf_text = (
                    "# Auto-generated by stb-convergence -- forces this folder's own swept parameter\n"
                    "# value plus full relaxation (positions + cell), regardless of --calc's own settings.\n"
                    f"{directive_line}\n"
                    "\n"
                    "MD.TypeOfRun       CG\n"
                    f"MD.Steps           {args.relax_steps}\n"
                    "MD.VariableCell    true\n"
                    f"MD.MaxStressTol    {args.max_stress_tol:g} GPa\n"
                    f"MD.MaxForceTol     {args.max_force_tol:g} eV/Ang\n"
                )

            shutil.copy(args.structure, os.path.join(folder, structure_basename))
            with open(os.path.join(folder, calc_basename), "w") as fh:
                fh.write(forced_calc_text)
            with open(os.path.join(folder, EXTRA_FDF_FILE), "w") as fh:
                fh.write(extra_fdf_text)

            for sym in species:
                copy_pseudo(args.pseudo_dir, sym, folder)

            folder_rows.append(([folder_name, p, detail,
                                f"{structure_basename}, {calc_basename}, {EXTRA_FDF_FILE}"], None))
            per_parameter_count[p] = per_parameter_count.get(p, 0) + 1

    print_table(["Folder", "Parameter", "Value", "Files"], folder_rows, f_out)

    print_section("[6] SUMMARY & NEXT STEPS", f_out)
    print_dual(f"{len(folder_rows)} folder(s) written under '{args.output_dir}' "
               f"({len(per_parameter_count)} parameter(s)):", f_out)
    for p, n in per_parameter_count.items():
        print_dual(f"  {p}: {n} folder(s)", f_out)
    if report_path:
        print_dual(f"Report            : {report_path}", f_out)
    print_dual(color_text(
        "\n[NOTE] stb-convergenceAnalysis (Stage 2) doesn't yet support this per-parameter "
        f"layout -- point --dir at one parameter's own subfolder (e.g. '{args.output_dir}/"
        f"{args.parameter[0]}'), not '{args.output_dir}' itself. A Stage 2 update is planned "
        "separately.", 'yellow'), f_out)
    print_dual(color_text("\nNext steps:", 'yellow'), f_out)
    print_dual(f"  1. Run SIESTA in every '{args.output_dir}/<parameter>/convergence_*/' folder.", f_out)
    print_dual("  2. Once they're done: stb-convergenceAnalysis --dir <one parameter's subfolder>", f_out)

    if f_out:
        f_out.close()


if __name__ == "__main__":
    main()
