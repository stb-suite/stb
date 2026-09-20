#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

from stb import __version__ as VERSION

import os
import sys
import argparse
import glob
from collections import Counter
from datetime import datetime
import numpy as np
from phonopy.interface.siesta import read_siesta, get_physical_units
from stb.core.cli import color_text, show_intro, print_dual, print_section, print_table
from stb.core.pseudopotentials import BANKS, resolve_pseudo_source, get_required_pseudos
from stb.core import kspace
from stb.core.phonon_workflow import build_phonon_displacements, write_displacement_folders
from stb.core.structure_io import read_md_state, prepend_include

REPORT_FILE = "raman_stage1.txt"
EXTRA_FDF_FILE = "config_extra.fdf"


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Stage 1 of 3: generates the phonon-displacement folders for a "
        "Raman spectrum workflow.", 'bold')}
Workflow: (1) this tool generates 'raman_study/phonon_disp/disp-*/' -- run
SIESTA in every folder yourself; (2) once they're done, run stb-ramanModes
to build FORCE_SETS from them, identify the Gamma-point vibrational modes,
and generate the +/-delta Optical-calculation displacement folders needed
for the Raman tensor; (3) run SIESTA in those, then stb-ramanAnalysis to
get the Raman-active frequencies/activities and spectrum. This whole
workflow computes its own phonons -- it does NOT need item 4's
stb-phononsCreate/stb-phononsPos run first (though the underlying Phonopy
machinery is shared, via core/phonon_workflow.py).""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s -s structure.fdf -c calc.fdf -p dojo\n"
               "  %(prog)s -s structure.fdf -c calc.fdf -p ~/pseudos -dim 2 2 2 -d 0.02\n"
               "  %(prog)s -s structure.fdf -c calc.fdf -p dojo -dim 5 5 1 --symprec 0.02\n"
               "  %(prog)s -s structure.fdf -c calc.fdf -p dojo -dim 3 3 1 --kgrid-density 0.05\n"
               "\n"
               "Notes:\n"
               "  - --symprec (default 0.01 Ang) controls how many displacement folders Phonopy\n"
               "    actually needs to write here: a HIGHER detected symmetry means FEWER, cheaper\n"
               "    SIESTA runs (see [4] PHONON DISPLACEMENTS' reduction-percent line). But its\n"
               "    real, more important effect is downstream, in stb-ramanModes: that tool\n"
               "    reloads these same force constants with its OWN --symprec and uses whichever\n"
               "    point group IT detects to decide which modes are Raman-active/forbidden and\n"
               "    which ones are truly degenerate (--use-symmetry/--skip-degenerate). Phonopy's\n"
               "    own raw default (1e-5) is far tighter than any real DFT relaxation's residual\n"
               "    noise (e.g. a cell angle of 119.999 deg instead of an exact 120 deg is normal\n"
               "    CG-relaxation noise, not a real symmetry breaking) and can silently DETECT A\n"
               "    LOWER SYMMETRY than the crystal actually has -- this is why 0.01 (pymatgen's\n"
               "    own default) is used here instead. Pass THE SAME --symprec to both stages, or\n"
               "    the mode numbering/classification stb-ramanModes reports may not match what\n"
               "    this stage assumed when it wrote (or skipped) a displacement here. If the\n"
               "    space/point group stb-ramanModes reports doesn't match what you expect for\n"
               "    this crystal, rerun both stages with a looser --symprec (e.g. 0.02-0.1).\n"
               "  - --kgrid-density (default 0.2 1/Ang, same convention as stb-kgrid/stb-strain/\n"
               "    stb-phononsCreate) is unrelated to symmetry -- it sets how finely the ELECTRONIC\n"
               "    structure is sampled for the SUPERCELL force calculation done here, which is\n"
               "    what the phonon frequencies computed by stb-ramanModes are actually built\n"
               "    from. A coarse density under-samples the electronic structure and can leave\n"
               "    the reported frequencies far from experiment -- most dramatically for a\n"
               "    metal/semimetal with a Kohn anomaly (e.g. graphene's G-band/E2g mode, whose\n"
               "    frequency is directly coupled to the Fermi-surface electrons: verified here to\n"
               "    move by >100 cm^-1 just from doubling the k-density, long after the supercell\n"
               "    size itself had stopped mattering). If a computed frequency looks far from the\n"
               "    literature/experimental value, suspect --kgrid-density before --dim (supercell\n"
               "    size) -- in practice the two are NOT interchangeable levers: --kgrid-density\n"
               "    converges the frequency for a semimetal much faster than growing --dim does.\n"
               "  - This flag only tunes the k-grid HERE, for [4]'s big supercell SCF. Stage 2's\n"
               "    own single-point Optical (dielectric-tensor) runs are small, unit-cell-sized\n"
               "    calculations that reuse whatever kgrid.MonkhorstPack is already inside the\n"
               "    --calc file YOU pass to stb-ramanModes -- there is no equivalent\n"
               "    --kgrid-density flag on that stage. To keep the whole workflow at a\n"
               "    consistent, converged k-density end to end, edit that block by hand in the\n"
               "    calc.fdf you give stb-ramanModes (it can be a different file from this one).\n"
    )

    parser.add_argument("-s", "--structure", type=str, default="structure.fdf",
                        help="Input structure file containing the unit cell (default: structure.fdf)")
    parser.add_argument("-c", "--calc", type=str, default="calc.fdf",
                        help="Input calculation parameters file (default: calc.fdf)")
    parser.add_argument("-p", "--pseudo-dir", type=str, default=".",
                        help=f"Pseudopotentials source: a bundled bank ({', '.join(BANKS)}) or a "
                             "folder path (default: current directory).")
    parser.add_argument("-dim", type=int, nargs=3, default=[2, 2, 2],
                        help="Supercell dimensions for the phonon force-constant calculation "
                             "(default: 2 2 2).")
    parser.add_argument("-d", "--distance", type=float, default=0.02,
                        help="Phonon displacement distance in Angstroms (default: 0.02).")
    parser.add_argument("--kgrid", type=int, nargs=3, default=None,
                        help="Explicit Monkhorst-Pack grid for the supercell's forced "
                             "single-point SCF (e.g. --kgrid 2 2 2). Overrides the automatic "
                             "--kgrid-density suggestion. Written into config_extra.fdf (see "
                             "[2]).")
    parser.add_argument("--kgrid-density", type=float, default=0.2,
                        help="Target k-point density (1/Ang), used to auto-suggest the supercell "
                             "k-grid when --kgrid isn't given (default: 0.2, same convention as "
                             "stb-kgrid/stb-strain/stb-phononsCreate).")
    parser.add_argument("--vacuum-gap", type=float, default=10.0,
                        help="Minimum gap (Ang) along an axis to consider it vacuum-padded, "
                             "for the supercell-dimension advisory (default: 10.0)")
    parser.add_argument("--symprec", type=float, default=0.01,
                        help="Symmetry-detection tolerance (Ang) Phonopy uses internally to "
                             "reduce how many displacements are actually needed (default: 0.01, "
                             "pymatgen's own default -- matches the rest of the suite, and "
                             "deliberately NOT Phonopy's own raw default of 1e-5, which is far "
                             "too tight for a real DFT-relaxed structure and can misdetect the "
                             "true point group -- see stb-ramanModes' own --symprec, which is "
                             "what actually drives the Raman-active/inactive classification. "
                             "Loosen further (e.g. 0.02-0.05) for a structure relaxed with a "
                             "looser force tolerance.")
    parser.add_argument("-O", "--output-dir", type=str, default="raman_study",
                        help="Root directory for the whole Raman workflow (default: raman_study) -- "
                             "phonon displacements are written under <output-dir>/phonon_disp/.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to <output-dir>/{REPORT_FILE}. Off by default.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-raman {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    try:
        args.pseudo_dir = resolve_pseudo_source(args.pseudo_dir)
    except ValueError as e:
        print(color_text(f"[ERROR] {e}", 'red'))
        sys.exit(1)

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("RAMAN SPECTRUM WORKFLOW -- STAGE 1: PHONON DISPLACEMENTS", 'bold'))
    print("-" * 60)

    print("\n[INFO] Validating input files ...")
    if not os.path.exists(args.structure):
        print(color_text(f"[ERROR] Structure file '{args.structure}' not found.", 'red'))
        sys.exit(1)
    if not os.path.exists(args.calc):
        print(color_text(f"[ERROR] Calculation file '{args.calc}' not found.", 'red'))
        sys.exit(1)

    print(f"[INFO] Reading unit cell from '{args.structure}' ...")
    try:
        unitcell = read_siesta(args.structure)
    except Exception as e:
        print(color_text(f"[ERROR] Failed to read {args.structure}. Make sure it's properly "
                          f"formatted.\nDetails: {e}", 'red'))
        sys.exit(1)

    output_root = args.output_dir
    phonon_disp_dir = os.path.join(output_root, "phonon_disp")
    existing_disps = sorted(glob.glob(os.path.join(phonon_disp_dir, "disp-*")))
    if existing_disps:
        print(color_text(
            f"\n[CRITICAL ERROR] '{phonon_disp_dir}' already contains {len(existing_disps)} "
            "displacement folder(s) from a previous run.", 'red'))
        print(color_text(
            "Regenerating on top of them can leave stale disp-* folders mixed in with the new "
            f"ones, silently corrupting FORCE_SETS during Stage 2. Remove or move aside "
            f"'{output_root}' and rerun.", 'yellow'))
        sys.exit(1)
    os.makedirs(phonon_disp_dir, exist_ok=True)

    report_path = os.path.join(output_root, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    print_dual(f"{color_text('===== RAMAN STAGE 1 REPORT (PHONON DISPLACEMENTS) =====', 'magenta')}", f_out)

    print_section('[0] RUN METADATA', f_out)
    print_dual(f"Date/time         : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", f_out)
    print_dual(f"Structure file    : {args.structure}", f_out)
    print_dual(f"Calc file         : {args.calc}", f_out)
    print_dual(f"Pseudo source     : {args.pseudo_dir}", f_out)
    print_dual(f"Supercell dim     : {args.dim[0]} x {args.dim[1]} x {args.dim[2]}", f_out)
    print_dual(f"Displacement dist.: {args.distance} Ang", f_out)
    print_dual(f"Output root       : {output_root}", f_out)
    print_dual(f"Vacuum-gap        : {args.vacuum_gap} Ang", f_out)
    print_dual(f"Symmetry tolerance: symprec={args.symprec:g} Ang", f_out)
    print_dual(f"Save report       : {'yes' if args.save_report else 'no'}", f_out)
    if report_path:
        print_dual(f"Report file       : {report_path}", f_out)

    print_section('[1] INPUT STRUCTURE', f_out)
    bohr_to_angstrom = get_physical_units().Bohr
    lattice_ang = np.array(unitcell.cell) * bohr_to_angstrom
    frac_coords = np.array(unitcell.scaled_positions)
    vacuum_axes = kspace.detect_vacuum_axes(frac_coords, lattice_ang, args.vacuum_gap)

    elements = sorted(set(unitcell.symbols))
    counts = Counter(unitcell.symbols)
    composition = ", ".join(f"{sp}{n}" for sp, n in counts.items())
    n_atoms = len(unitcell.symbols)
    cell_volume = abs(np.linalg.det(lattice_ang))
    supercell_atoms = n_atoms * args.dim[0] * args.dim[1] * args.dim[2]
    print_table(["Quantity", "Value"], [
        (["Composition", composition], None),
        (["Total atoms (unit cell)", str(n_atoms)], None),
        (["Cell volume (unit cell)", f"{cell_volume:.4f} Ang^3"], None),
        (["Dimensionality", kspace.dimensionality_label(vacuum_axes)], None),
        (["Supercell", f"{args.dim[0]} x {args.dim[1]} x {args.dim[2]} "
                        f"({supercell_atoms} atoms)"], None),
    ], f_out)

    vacuum_dims_requested = [axis for axis, is_vac in zip('abc', vacuum_axes)
                              if is_vac and args.dim['abc'.index(axis)] > 1]
    if vacuum_dims_requested:
        print_dual(color_text(
            f"[WARNING] -dim requests more than 1 repetition along vacuum-padded "
            f"axis/axes {', '.join(vacuum_dims_requested)} (gap >= {args.vacuum_gap} Ang). "
            "Replicating a supercell across vacuum only multiplies the SIESTA cost without "
            "adding real periodicity -- consider -dim 1 on that axis.", 'yellow'), f_out)

    print_section('[2] SINGLE-POINT SCF ENFORCEMENT', f_out)
    with open(args.calc) as f:
        original_calc_text = f.read()
    structure_basename = os.path.basename(args.structure)
    calc_basename = os.path.basename(args.calc)
    if "%include" not in original_calc_text or structure_basename not in original_calc_text:
        print_dual(color_text(
            f"[NOTE] Could not confirm '{calc_basename}' references '{structure_basename}' via "
            f"%include -- if your calc.fdf doesn't already include the structure file by this "
            f"exact name, add '%include {structure_basename}' to it (or the equivalent for your "
            "own convention) so SIESTA picks up each displaced supercell's exact geometry.",
            'yellow'), f_out)

    before = read_md_state(original_calc_text)
    steps_label = f"{before['steps_key']}={before['steps']}" if before['steps_key'] else "(absent)"
    print_dual("Calc template (current state, before forcing):", f_out)
    print_dual(f"  MD.TypeOfRun={before['typeofrun'] or '(absent)'}  Steps: {steps_label}  "
               f"MD.VariableCell={before['variablecell'] or '(absent)'}", f_out)
    print_dual(
        "Every generated folder is forced to a pure single-point SCF -- NO ionic relaxation, "
        "NO cell relaxation -- regardless of the state above, since each displaced supercell's "
        "force must be measured at EXACTLY the displaced geometry Phonopy generated: if SIESTA "
        "relaxed it back toward equilibrium, the measured force would no longer correspond to "
        "the intended displacement, silently corrupting the finite-difference force constants "
        "(and everything derived from them downstream -- phonon frequencies, eigenvectors, and "
        f"the Raman tensor itself). Forced via '%include {EXTRA_FDF_FILE}' PREPENDED at the very "
        f"top of every generated {calc_basename} (before your own directives, including the "
        "structure %include) -- SIESTA's fdf reader is first-occurrence-wins for duplicate "
        "labels, so this ordering guarantees the forced values win even if your own template "
        "already sets any of them (the common case for a real relaxation calc.fdf reused here).",
        f_out)

    # Supercell k-grid: --calc's own k-grid was tuned for the (smaller) unit
    # cell, so it over-samples the supercell if left as-is -- compute (or
    # take explicitly via --kgrid) a grid sized for the supercell itself, at
    # --kgrid-density (default 0.2, same convention as stb-kgrid/stb-strain/
    # stb-phononsCreate), and put it FIRST in config_extra.fdf (before the MD
    # block below) so it's the first directive SIESTA reads.
    supercell_lattice = np.diag(args.dim) @ lattice_ang
    if args.kgrid is not None:
        kgrid = tuple(args.kgrid)
        kgrid_source = "explicit --kgrid"
    else:
        try:
            kgrid = tuple(kspace.compute_monkhorts(
                supercell_lattice[0], supercell_lattice[1], supercell_lattice[2],
                args.kgrid_density, vacuum_axes))
        except ValueError as e:
            print_dual(color_text(f"[ERROR] {e}", 'red'), f_out)
            if f_out:
                f_out.close()
            sys.exit(1)
        kgrid_source = f"auto-suggested, density={args.kgrid_density:g} 1/Ang"
    print_dual(f"\nSupercell k-grid  : {kgrid[0]} {kgrid[1]} {kgrid[2]} ({kgrid_source})", f_out)

    extra_fdf_text = (
        "# Auto-generated by stb-raman.\n"
        f"# Supercell k-grid ({kgrid_source}) -- takes precedence over --calc's own k-grid,\n"
        "# which was tuned for the smaller unit cell.\n"
        f"kgrid.MonkhorstPack   [{kgrid[0]}  {kgrid[1]}  {kgrid[2]}]\n"
        "\n"
        "# Forces a pure single-point SCF (no ionic or cell relaxation) at each displaced\n"
        "# supercell's exact geometry, regardless of --calc's own settings.\n"
        "MD.TypeOfRun       CG\n"
        "MD.Steps           0\n"
        "MD.VariableCell    false\n"
    )
    print_dual(f"\n{EXTRA_FDF_FILE} (written into every generated folder):", f_out)
    for line in extra_fdf_text.rstrip("\n").split("\n"):
        print_dual(f"  {line}", f_out)
    forced_calc_text = prepend_include(original_calc_text, EXTRA_FDF_FILE)

    print_section('[3] PSEUDOPOTENTIALS', f_out)
    print_dual(f"Source            : {args.pseudo_dir}", f_out)
    pseudos_to_copy, missing = get_required_pseudos(elements, args.pseudo_dir)
    print_table(["Species", "Status"], [
        ([sp, "MISSING" if sp in missing else "found"], 'yellow' if sp in missing else None)
        for sp in elements
    ], f_out)

    if missing:
        print_dual(color_text(
            f"[CRITICAL ERROR] Missing pseudopotentials for: {', '.join(missing)} -- add the "
            f"necessary '{missing[0]}.psf' or '{missing[0]}.psml' file(s) into "
            f"'{args.pseudo_dir}' and rerun.", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)
    print_dual(color_text(
        "[OK] All required pseudopotentials found -- will be copied into every generated folder.",
        'green'), f_out)

    print_section('[4] PHONON DISPLACEMENTS', f_out)
    print(f"[INFO] Generating supercell {args.dim} with {args.distance} Ang displacements ...")
    supercell_matrix = [
        [args.dim[0], 0, 0],
        [0, args.dim[1], 0],
        [0, 0, args.dim[2]]
    ]
    # phonopy's SIESTA interface keeps the structure internally in bohr --
    # distance conversion is handled inside build_phonon_displacements()
    # (see core/phonon_workflow.py for the ~1.89x pitfall this avoids).
    phonon, supercells = build_phonon_displacements(
        unitcell, supercell_matrix, args.distance, symprec=args.symprec)

    n_used = len(phonon.dataset['first_atoms'])
    n_naive = phonon.dataset['natom'] * 6
    reduction_pct = 100.0 * (1.0 - n_used / n_naive) if n_naive else 0.0
    print_dual(f"Displacements needed  : {color_text(str(n_used), 'green')} "
                f"(vs. {n_naive} without symmetry reduction -- {reduction_pct:.0f}% fewer "
                "independent SIESTA runs)", f_out)

    print_dual(f"Building {len(supercells)} displacement folders in '{phonon_disp_dir}' ...", f_out)
    folders, yaml_path = write_displacement_folders(
        phonon_disp_dir, phonon, supercells, args.structure, args.calc, pseudos_to_copy)

    # write_displacement_folders copies --calc verbatim -- overwrite each
    # folder's copy with the forced single-point/supercell-k-grid version
    # from [2], and add the config_extra.fdf sidecar, right after.
    for folder in folders:
        with open(os.path.join(folder, calc_basename), "w") as fh:
            fh.write(forced_calc_text)
        with open(os.path.join(folder, EXTRA_FDF_FILE), "w") as fh:
            fh.write(extra_fdf_text)

    print_dual(f"Saved Phonopy metadata to '{yaml_path}'", f_out)

    print_section('[5] SUMMARY & NEXT STEPS', f_out)
    print_dual(f"Displacement folders : {len(folders)} (disp-001 .. disp-{len(folders):03d})", f_out)
    if report_path:
        print_dual(f"Report               : {report_path}", f_out)
    print_dual(f"Files                : {yaml_path}, {phonon_disp_dir}/disp-*/", f_out)
    print_dual(color_text(
        f"\n[NOTE] '{os.path.basename(args.calc)}' was forced to single-point SCF with its own "
        f"supercell k-grid ({kgrid[0]} {kgrid[1]} {kgrid[2]}, see [2]) in every disp-* folder.",
        'yellow'), f_out)
    print_dual(color_text("\nNext steps:", 'yellow'), f_out)
    print_dual(f"  1. Run SIESTA in every '{phonon_disp_dir}/disp-*/' folder.", f_out)
    print_dual(f"  2. Once they're done, run: stb-ramanModes --directory {output_root}", f_out)

    if f_out:
        f_out.close()

    print("\n[INFO] Complete job!")
    print("\n" + "-" * 60)
    print(color_text("Phonon displacement folders ready for Stage 2 (stb-ramanModes).\n", 'bold'))


if __name__ == "__main__":
    main()
