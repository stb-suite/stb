#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""CLI-only tool (deliberately NOT wired into stb-suite's interactive menu):
one call = one real-DFT NEB step, meant to be called from inside a cluster
submission script's loop (see stb-neb --mode 3/--mode 4's own printed
example), alternating with real SIESTA runs -- there is no live Python
process to keep an optimizer running in memory between cycles, since each
cycle's forces only exist once a real (queued, possibly hours-later) SIESTA
job has actually finished.

Mechanics, verified directly against the installed ASE source before
writing this (not assumed):
- ase.optimize.optimize.Optimizer.__init__(..., restart=<path>) calls
  self.read() if that file exists, else self.initialize() -- native,
  automatic restart.
- ase.optimize.FIRE.step() ends every call with self.dump((self.vel,
  self.dt)) -- so even a single opt.run(fmax=..., steps=1) call correctly
  persists FIRE's velocity/timestep state to disk (JSON, via
  ase.io.jsonio.write_json) for the NEXT separate process invocation to
  pick back up.
- core/mace_relax.py::relax_neb already proves FIRE(neb_object, ...) works
  directly on an ase.mep.neb.NEB (not just a plain Atoms) -- same pattern
  here, except the "calculator" is ase.calculators.singlepoint.
  SinglePointCalculator wrapping each image's already-finished SIESTA
  energy/forces, not a live MACE calculator.

Endpoints (images[0]/images[-1], the already-relaxed site_A/site_B from
stb-nebSites) are never moved -- ase.mep.neb.BaseNEB.get_positions/
set_positions only ever touch the interior images, same guarantee
core/mace_relax.py::relax_neb already documents and relies on.
"""

VERSION = "1.0.0"

import os
import re
import sys
import glob
import argparse
import numpy as np
from pymatgen.io.ase import AseAtomsAdaptor
from stb.core import structure_io, siesta_log
from stb.core.cli import color_text, show_intro, print_dual, print_section
from stb.neb import write_image_folder

CONVERGED_SENTINEL = "NEB_CONVERGED"
_CYCLE_RE = re.compile(r"^cycle_(\d+)$")


def find_latest_cycle(root_dir):
    """Returns (cycle_num, cycle_dir) for the highest-numbered 'cycle_NN'
    directory directly under root_dir, or (None, None) if none exist.
    """
    best = None
    for entry in sorted(os.listdir(root_dir)):
        m = _CYCLE_RE.match(entry)
        full = os.path.join(root_dir, entry)
        if m and os.path.isdir(full):
            n = int(m.group(1))
            if best is None or n > best[0]:
                best = (n, full)
    return best if best else (None, None)


def read_image_results(cycle_dir, out_filename):
    """Reads every 'image_*' folder in cycle_dir: structure (positions),
    energy (out_filename, e.g. calc.out), and forces (<SystemLabel>.FA,
    label auto-detected per folder via structure_io.find_fdf_system_label
    -- the same convention stb-status/stb-archive already use). Returns
    (labels, pmg_structures, energies, forces_list) all sorted by image
    label, or raises ValueError naming the first folder that isn't ready
    yet (missing .out/.FA, or unparseable) -- a NEB step needs EVERY
    image's real forces at once, so a partial cycle can't take a step at
    all, unlike stb-nebAnalysis's per-folder skip-and-continue.
    """
    image_dirs = sorted(
        d for d in os.listdir(cycle_dir)
        if os.path.isdir(os.path.join(cycle_dir, d)) and d.startswith("image_")
    )
    if not image_dirs:
        raise ValueError(f"No 'image_*' folders found in '{cycle_dir}'.")

    labels, structures, energies, forces_list = [], [], [], []
    for label in image_dirs:
        image_dir = os.path.join(cycle_dir, label)
        fdf_path = os.path.join(image_dir, "structure.fdf")
        if not os.path.isfile(fdf_path):
            raise ValueError(f"'{image_dir}' has no structure.fdf.")
        structures.append(structure_io.to_pymatgen(structure_io.read_fdf(fdf_path)))

        out_path = os.path.join(image_dir, out_filename)
        energy = siesta_log.get_free_energy(out_path) if os.path.isfile(out_path) else None
        if energy is None:
            raise ValueError(f"'{image_dir}': could not read energy from '{out_filename}' "
                              "-- has SIESTA finished here yet?")
        energies.append(energy)

        system_label = structure_io.find_fdf_system_label(image_dir)
        fa_path = os.path.join(image_dir, f"{system_label}.FA") if system_label else None
        forces = siesta_log.read_fa_forces(fa_path) if fa_path and os.path.isfile(fa_path) else None
        if forces is None:
            raise ValueError(f"'{image_dir}': could not read forces from "
                              f"'{system_label}.FA' -- has SIESTA finished here yet?")
        forces_list.append(forces)

        labels.append(label)

    return labels, structures, energies, forces_list


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("One real-DFT NEB step: reads the latest cycle_NN/image_*/ "
        "(finished SIESTA single-points), takes exactly one climbing-image NEB step (ASE's own "
        "FIRE optimizer, restart-persisted between calls), and either writes cycle_{N+1}/ or "
        "declares convergence. Meant to be called in a loop from a cluster submission script "
        "-- see stb-neb --mode 3/--mode 4's own printed example. Not available in the "
        "interactive stb-suite menu.", 'bold')}""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s --dir .\n"
               "  %(prog)s --dir . --fmax 0.03 --climb-after 5\n"
    )
    parser.add_argument("--dir", type=str, default=".",
                         help="Root directory containing every 'cycle_NN/' (default: current "
                              "directory).")
    parser.add_argument("--file", type=str, default="calc.out",
                         help="SIESTA output filename inside each image folder (default: "
                              "calc.out).")
    parser.add_argument("--k", type=float, default=0.1,
                         help="NEB spring constant, eV/Ang^2 (default: 0.1, ASE's own default).")
    parser.add_argument("--fmax", type=float, default=0.05,
                         help="Force convergence threshold, eV/Ang (default: 0.05).")
    parser.add_argument("--climb-after", type=int, default=0, metavar="N",
                         help="Only enable the climbing-image correction once the CURRENT cycle "
                              "number is >= N (default: 0, climb from the start -- appropriate "
                              "once a MACE-shaped path already has a sensible shape, e.g. "
                              "stb-neb --mode 3; a larger value is recommended for --mode 4, "
                              "which starts from a raw interpolated path).")
    parser.add_argument("--save-report", action="store_true",
                         help="Also persist this cycle's report to <dir>/neb_cycle_report.txt. "
                              "Off by default.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-nebCycle {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    if not os.path.isdir(args.dir):
        print(color_text(f"[ERROR] '{args.dir}' not found.", 'red'))
        sys.exit(1)

    converged_path = os.path.join(args.dir, CONVERGED_SENTINEL)
    if os.path.isfile(converged_path):
        print(color_text(
            f"[CONVERGED] '{converged_path}' already exists from a previous call -- nothing "
            "to do. Delete it if you want to force another step.", 'yellow'))
        sys.exit(0)

    report_path = os.path.join(args.dir, "neb_cycle_report.txt") if args.save_report else None
    f_out = open(report_path, 'a') if report_path else None

    print_dual(color_text("===== STB-NEBCYCLE =====", 'magenta'), f_out)

    cycle_num, cycle_dir = find_latest_cycle(args.dir)
    if cycle_dir is None:
        print_dual(color_text(
            f"[ERROR] No 'cycle_NN' folder found under '{args.dir}'. Did you run stb-neb "
            "--mode 3 or --mode 4?", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    climb = cycle_num >= args.climb_after

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"  Root dir        : {args.dir}", f_out)
    print_dual(f"  Current cycle   : {os.path.basename(cycle_dir)}", f_out)
    print_dual(f"  Spring k        : {args.k} eV/Ang^2", f_out)
    print_dual(f"  fmax target     : {args.fmax} eV/Ang", f_out)
    print_dual(f"  Climbing image  : {'ON' if climb else 'OFF'} (cycle {cycle_num} "
                f"{'>=' if climb else '<'} --climb-after {args.climb_after})", f_out)

    print_section("[1] READING SIESTA RESULTS", f_out)
    try:
        labels, structures, energies, forces_list = read_image_results(cycle_dir, args.file)
    except ValueError as e:
        print_dual(color_text(f"[ERROR] {e}", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)
    for label, energy in zip(labels, energies):
        print_dual(f"  {color_text('[OK]', 'green')} {os.path.basename(cycle_dir)}/{label} : "
                    f"E = {energy:.6f} eV", f_out)

    from ase.calculators.singlepoint import SinglePointCalculator
    from ase.mep.neb import NEB
    from ase.optimize import FIRE

    images = []
    for pmg, energy, forces in zip(structures, energies, forces_list):
        atoms = AseAtomsAdaptor.get_atoms(pmg)
        atoms.calc = SinglePointCalculator(atoms, energy=energy, forces=forces)
        images.append(atoms)

    neb = NEB(images, k=args.k, climb=climb, method="improvedtangent")
    with np.errstate(invalid="raise"):
        try:
            raw_forces = neb.get_forces()
        except FloatingPointError:
            raw_forces = None
    if raw_forces is None or not np.isfinite(raw_forces).all():
        # ASE's "improvedtangent" tangent estimate normalizes by its own
        # norm (tangent /= np.linalg.norm(tangent)) -- degenerate for two
        # neighboring images sitting at (near-)identical positions/
        # energies, a 0/0 that silently produces NaN forces instead of
        # raising. Verified live: NaN then propagates into the NEXT
        # cycle's written positions (a corrupted structure.fdf) and into
        # FIRE's own persisted velocity state, poisoning every cycle after
        # it. Caught here, immediately, rather than writing a corrupted
        # cycle_{N+1}/ that would only surface as a mysterious failure
        # later.
        print_dual(color_text(
            "[ERROR] NEB force evaluation produced a non-finite (NaN/Inf) result -- likely two "
            "neighboring images have become degenerate (near-identical position or energy), "
            "which the 'improvedtangent' tangent estimate can't normalize (0/0). This usually "
            "means the path has collapsed; inspect the last few cycles' geometries before "
            "continuing (do not delete neb_cycle_state.json and retry blindly -- it may already "
            "hold a poisoned velocity from this same event).", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)
    max_perp_force = float(np.abs(raw_forces).max())

    print_section("[2] NEB STEP", f_out)
    print_dual(f"  Max residual force (this cycle, pre-step) : {max_perp_force:.4f} eV/Ang", f_out)

    state_path = os.path.join(args.dir, "neb_cycle_state.json")
    had_state = os.path.isfile(state_path)
    opt = FIRE(neb, restart=state_path, logfile=None)

    # Deliberately opt.step() directly, NOT opt.run(fmax=..., steps=1):
    # Dynamics.irun() (what run() drives) re-evaluates neb.get_forces() on
    # the JUST-MOVED geometry right after step() to decide/report
    # convergence for its own next iteration -- but every image's
    # "calculator" here is a SinglePointCalculator wrapping one finished
    # SIESTA result, which cannot recompute anything for a new geometry it
    # was never constructed with (verified live: this crashes with
    # ase.calculators.calculator.PropertyNotImplementedError the moment
    # run() tries it). Calling step() directly only ever evaluates forces
    # on the CURRENT, already-known geometry (needed for the gradient),
    # then moves the positions and returns -- no post-move re-evaluation.
    # Convergence is instead checked here, on the PRE-step forces already
    # computed above for the report, exactly equivalent to what run()
    # would have checked before deciding whether a step was even needed.
    if max_perp_force < args.fmax:
        converged = True
    else:
        opt.step()
        converged = False
    print_dual(f"  FIRE step taken (optimizer state: "
                f"{'loaded from previous cycle' if had_state else 'fresh, first cycle'})", f_out)

    print_section("[3] RESULT", f_out)
    if converged:
        with open(converged_path, "w") as f:
            f.write(f"Converged at {os.path.basename(cycle_dir)}, max force below "
                     f"{args.fmax} eV/Ang.\n")
        print_dual(color_text(
            f"[CONVERGED] Max residual force is already below --fmax ({args.fmax} eV/Ang) -- "
            f"wrote {converged_path}. No new cycle folder written; stop the sub.sh loop.",
            'green'), f_out)
    else:
        next_dir = os.path.join(args.dir, f"cycle_{cycle_num + 1:02d}")
        calc_fdf_path = os.path.join(cycle_dir, labels[0], "calc.fdf")
        with open(calc_fdf_path) as f:
            calc_text = f.read()
        for label, atoms in zip(labels, images):
            pmg_new = AseAtomsAdaptor.get_structure(atoms)
            species_meta = structure_io.species_dict(structure_io.read_fdf(
                os.path.join(cycle_dir, label, "structure.fdf")))
            write_image_folder(os.path.join(next_dir, label), pmg_new, calc_text, species_meta,
                                pp_path=os.path.join(cycle_dir, label))
            print_dual(f"  {color_text('[OK]', 'green')} {next_dir}/{label}", f_out)
        print_dual(f"[OK] {len(labels)} image folder(s) written under '{next_dir}' -- run SIESTA "
                    "there next, then call stb-nebCycle again.", f_out)

    if f_out:
        f_out.close()

    print(f"\n{color_text('Success:', 'green')} cycle {cycle_num} processed.")
    if report_path:
        print(f"Report: {report_path}")


if __name__ == "__main__":
    main()
