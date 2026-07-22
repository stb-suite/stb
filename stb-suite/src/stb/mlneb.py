#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""Standalone MACE-driven NEB (nudged elastic band) -- the seventh tool in
the ML Simulations category. stb-neb already has a `--ml-neb` preview mode
(a real climbing-image NEB on MACE-MP-0, used as a fast look before
committing to the real SIESTA-based `image_NN/` folders that tool
generates), but that mode is bolted onto a DFT-oriented workflow: it always
needs `--calc`/pseudopotentials even though the ML path never touches
either, and its `--ml-model` only accepts `small`/`medium`/`large` (no
`--custom-model` support). This tool is the pure-MACE extraction of that
same capability for users who only want the ML-level transition path/
barrier estimate -- e.g. defect/vacancy migration, an adsorbate hopping
between sites -- with no DFT step involved at all.

Deliberately reuses (not duplicates) neb.py's own pure-geometry helpers --
check_composition_match, wrap_into_cell, resolve_lattice_mismatch,
linear_interpolate_images, idpp_refine_images, compute_frozen_indices,
cumulative_reaction_coordinates, check_path_quality, write_path_trajectory,
write_ml_preview_plot -- since interpolating/quality-checking a reaction
path is exactly the same problem whether the images end up going to SIESTA
or straight into a MACE NEB; only write_image_folder (SIESTA input/
pseudopotential specific) is not reused, since this tool doesn't generate
any DFT input. The actual NEB relaxation reuses core/mace_relax.py's
relax_neb (already the same function stb-neb's --ml-neb mode calls) --
climbing-image NEB, two-stage convergence (climb=False first to let the
band find its shape, then climb=True to refine the true saddle point).
"""

VERSION = "1.0.0"

import os
import sys
import argparse
from datetime import datetime

import numpy as np
from pymatgen.io.ase import AseAtomsAdaptor

from stb.core import structure_io, mace_relax
from stb.core.cli import color_text, show_intro, print_dual, print_section
from stb.core.deps import require_mace
from stb.neb import (
    check_composition_match, wrap_into_cell, resolve_lattice_mismatch,
    linear_interpolate_images, idpp_refine_images, compute_frozen_indices,
    cumulative_reaction_coordinates, check_path_quality, write_path_trajectory,
    write_ml_preview_plot,
)

require_mace()

REPORT_FILE = "stb_mlneb_report.txt"


def run_single_neb(initial_pmg, final_pmg_matched, calc, args, f_out, tag=""):
    """Runs the shared endpoint-pre-relax + interpolate (+ optional IDPP)
    + climbing-image NEB pipeline on one (initial, final) pymatgen
    Structure pair (both already wrapped into [0, 1) and lattice-matched --
    see wrap_into_cell/resolve_lattice_mismatch). `args` needs the same
    attribute names main() itself parses (n_images, autosort_tol, idpp,
    idpp_fmax, idpp_steps, pre_relax_endpoints, freeze_substrate,
    freeze_threshold, k, fmax, max_steps) -- any caller's argparse
    Namespace works as long as it carries these, same duck-typing
    convention elastic_analysis.py's emit_elastic_report already uses.
    `tag` labels progress lines only, for a caller running this many times
    (e.g. stb-mldiffusion, one call per candidate defect hop).

    Returns (barrier, dE, ase_images, pmg_images, converged).
    """
    prefix = f"[{tag}] " if tag else ""
    if args.pre_relax_endpoints:
        for label, pmg in (("initial", initial_pmg), ("final", final_pmg_matched)):
            ase_atoms = AseAtomsAdaptor.get_atoms(pmg)
            ase_atoms.calc = calc
            f0 = np.abs(ase_atoms.get_forces()).max()
            converged, steps = mace_relax.relax(ase_atoms, calc, fmax=args.fmax, max_steps=args.max_steps)
            f1 = np.abs(ase_atoms.get_forces()).max()
            relaxed_pmg = AseAtomsAdaptor.get_structure(ase_atoms)
            if label == "initial":
                initial_pmg = relaxed_pmg
            else:
                final_pmg_matched = relaxed_pmg
            print_dual(f"{prefix}[{label}] max|F| {f0:.4f} -> {f1:.4f} eV/Ang "
                       f"({'converged' if converged else 'hit step cap'}, {steps} steps)", f_out)
    else:
        print_dual(color_text(
            f"{prefix}[WARNING] --no-pre-relax-endpoints: using the endpoints as given -- the "
            "reported barrier may be contaminated by residual forces already present.", 'yellow'), f_out)

    pmg_images = linear_interpolate_images(initial_pmg, final_pmg_matched, args.n_images,
                                           autosort_tol=args.autosort_tol)
    print_dual(f"{prefix}[OK] Linear interpolation: {len(pmg_images)} images", f_out)
    if args.idpp:
        pmg_images = idpp_refine_images(pmg_images, fmax=args.idpp_fmax, steps=args.idpp_steps)
        print_dual(f"{prefix}[OK] IDPP refinement applied to interior images", f_out)

    reaction_coords = cumulative_reaction_coordinates(pmg_images)
    check_path_quality(pmg_images, reaction_coords, f_out)

    ase_images = [AseAtomsAdaptor.get_atoms(s) for s in pmg_images]
    if args.freeze_substrate:
        frozen_indices = compute_frozen_indices(pmg_images, threshold=args.freeze_threshold)
        if frozen_indices:
            from ase.constraints import FixAtoms
            for image in ase_images:
                image.set_constraint(FixAtoms(indices=frozen_indices))
        print_dual(f"{prefix}Freezing {len(frozen_indices)}/{len(ase_images[0])} atom(s) with "
                   f"< {args.freeze_threshold} Ang displacement between endpoints", f_out)

    converged, s1, s2, energies = mace_relax.relax_neb(
        ase_images, calc, k=args.k, fmax=args.fmax, max_steps=args.max_steps)
    print_dual(f"{prefix}[{'converged' if converged else 'hit step cap, not fully converged'}] "
               f"after {s1} (stage 1) + {s2} (stage 2, climbing) step(s)", f_out)

    from ase.mep.neb import NEBTools
    barrier, dE = NEBTools(ase_images).get_barrier(fit=True)
    print_dual(f"{prefix}Barrier (forward) : {color_text(f'{barrier:.4f}', 'bold')} eV", f_out)
    print_dual(f"{prefix}Reaction energy   : {dE:.4f} eV (fitted spline over the relaxed band)", f_out)

    return barrier, dE, ase_images, pmg_images, converged


def main():
    parser = argparse.ArgumentParser(
        description="Standalone climbing-image NEB (nudged elastic band) driven entirely by "
                     "a MACE potential -- the MACE-MP-0 foundation model, or a custom model "
                     "fine-tuned on your own SIESTA data via stb-mlffAnalysis "
                     "(--custom-model). Interpolates images between two already-relaxed "
                     "endpoint structures (--initial/--final), optionally refines the path "
                     "with IDPP, runs a real climbing-image NEB with MACE, and reports the "
                     "barrier + reaction energy. No SIESTA and no calc.fdf/pseudopotentials "
                     "needed -- for the real DFT-level path, see stb-neb "
                     "(and its own --ml-neb preview mode, which this tool's machinery was "
                     "extracted from).",
        epilog="Example usage:\n"
               "  stb-mlneb --initial initial.fdf --final final.fdf\n"
               "  stb-mlneb --initial initial.fdf --final final.fdf --custom-model "
               "mlff_model.model --n-images 9 --freeze-substrate\n",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--initial", required=True, help="Initial (relaxed) endpoint structure (.fdf).")
    parser.add_argument("--final", required=True, help="Final (relaxed) endpoint structure (.fdf).")
    parser.add_argument("--model", choices=["small", "medium", "large"], default="small",
                        help="MACE-MP-0 model size (default: small). Ignored if --custom-model is given.")
    parser.add_argument("--custom-model", default=None, metavar="PATH",
                        help="Use a custom MACE model file instead of the MACE-MP-0 foundation "
                             "potential -- e.g. one fine-tuned via stb-mlffAnalysis. Overrides --model.")
    parser.add_argument("--n-images", type=int, default=7,
                        help="Total images including both endpoints (default: 7).")
    parser.add_argument("--autosort-tol", type=float, default=0.5,
                        help="Atom-correspondence tolerance between --initial/--final, Ang "
                             "(default: 0.5, pymatgen's own documented value).")
    parser.add_argument("--no-idpp", dest="idpp", action="store_false",
                        help="Skip IDPP (Image Dependent Pair Potential) refinement of the "
                             "interior images before the MACE NEB. On by default -- a purely "
                             "classical/geometric smoothing that usually gives a much better "
                             "starting guess than plain linear interpolation.")
    parser.add_argument("--idpp-fmax", type=float, default=0.1,
                        help="Force convergence for IDPP refinement (default: 0.1).")
    parser.add_argument("--idpp-steps", type=int, default=100,
                        help="Max steps for IDPP refinement (default: 100).")
    parser.add_argument("--no-pre-relax-endpoints", dest="pre_relax_endpoints", action="store_false",
                        help="Skip the initial MACE relax (positions only) of both endpoints. "
                             "On by default, so the reported barrier isn't contaminated by "
                             "residual forces already present in --initial/--final.")
    parser.add_argument("--freeze-substrate", action="store_true",
                        help="Freeze atoms that barely move between --initial/--final (< "
                             "--freeze-threshold) during the NEB relax -- fewer degrees of "
                             "freedom, and avoids a spectator atom spuriously wandering during "
                             "a fast/loose ML relax. Off by default.")
    parser.add_argument("--freeze-threshold", type=float, default=0.3,
                        help="Displacement threshold for --freeze-substrate, Ang (default: 0.3).")
    parser.add_argument("--k", type=float, default=0.1,
                        help="NEB spring constant, eV/Ang^2 (default: 0.1).")
    parser.add_argument("--fmax", type=float, default=0.05,
                        help="Force convergence threshold for the NEB relax, eV/Ang (default: 0.05).")
    parser.add_argument("--max-steps", type=int, default=200,
                        help="Max steps per NEB stage (default: 200).")
    parser.add_argument("--device", choices=["cpu", "cuda"], default="cpu",
                        help="Device to run the model on (default: cpu).")
    parser.add_argument("-o", "--output-dir", default="mlneb_out",
                        help="Output directory for all files (default: mlneb_out).")
    parser.add_argument("--save-images", action="store_true",
                        help="Also write each relaxed image as its own image_NN.fdf (bare "
                             "structure, no calc/pseudopotentials) -- e.g. to hand off to a "
                             "real SIESTA NEB later. Off by default.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-mlneb {VERSION}")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - ML NEB",
            "Climbing-image nudged elastic band, driven entirely by a MACE potential",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("ML NEB (CLIMBING-IMAGE):", 'bold'))
    print("-" * 60)

    os.makedirs(args.output_dir, exist_ok=True)
    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    def fail(message):
        print_dual(color_text(f"[FAIL] {message}", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_dual(color_text("===== STB-MLNEB REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Initial structure : {args.initial}", f_out)
    print_dual(f"Final structure   : {args.final}", f_out)
    print_dual(f"Model             : "
               f"{'custom (' + args.custom_model + ')' if args.custom_model else f'MACE-MP-0 ({args.model})'}", f_out)
    print_dual(f"N images          : {args.n_images}", f_out)
    print_dual(f"IDPP refinement   : {'yes' if args.idpp else 'no'}", f_out)
    if report_path:
        print_dual(f"Report file       : {report_path}", f_out)

    print_section("[1] READING STRUCTURES", f_out)

    if not os.path.isfile(args.initial):
        fail(f"Initial structure file '{args.initial}' not found.")
    if not os.path.isfile(args.final):
        fail(f"Final structure file '{args.final}' not found.")
    if args.custom_model and not os.path.isfile(args.custom_model):
        fail(f"--custom-model file not found: {args.custom_model}")

    try:
        initial_structure = structure_io.read_fdf(args.initial)
        final_structure = structure_io.read_fdf(args.final)
    except Exception as e:
        fail(f"Could not read structure(s): {e}")

    check_composition_match(initial_structure, final_structure)
    species_meta = structure_io.species_dict(initial_structure)
    print_dual(f"[OK] Composition matches: {len(initial_structure.atoms)} atom(s)", f_out)

    initial_pmg = wrap_into_cell(structure_io.to_pymatgen(initial_structure))
    final_pmg = wrap_into_cell(structure_io.to_pymatgen(final_structure))
    final_pmg_matched = resolve_lattice_mismatch(initial_pmg, final_pmg)

    model_arg = args.custom_model if args.custom_model else args.model
    calc = mace_relax.get_calculator(model=model_arg, device=args.device)

    print_section("[2] NEB PIPELINE (MACE)", f_out)
    barrier, dE, ase_images, pmg_images, converged = run_single_neb(
        initial_pmg, final_pmg_matched, calc, args, f_out)
    print_dual(color_text(
        "[NOTE] An ML-level estimate from a fast surrogate potential, not a DFT-level "
        "barrier -- treat as a screening/preview number, not a final result.", 'yellow'), f_out)

    print_section("[3] SUMMARY & FILES", f_out)
    preview_path = os.path.join(args.output_dir, "neb_preview.png")
    write_ml_preview_plot(ase_images, preview_path)
    print_dual(f"[OK] Saved {preview_path}", f_out)

    traj_path = os.path.join(args.output_dir, "neb_path.xyz")
    write_path_trajectory(pmg_images, traj_path)
    print_dual(f"[OK] Saved {traj_path}", f_out)

    if args.save_images:
        pmg_relaxed = [AseAtomsAdaptor.get_structure(a) for a in ase_images]
        for i, pmg in enumerate(pmg_relaxed):
            fdf_structure = structure_io.from_pymatgen(pmg, species_meta=species_meta,
                                                        coord_format="fractional")
            image_path = os.path.join(args.output_dir, f"image_{i:02d}.fdf")
            structure_io.write_fdf(fdf_structure, image_path)
        print_dual(f"[OK] Saved image_00.fdf .. image_{len(pmg_relaxed)-1:02d}.fdf", f_out)

    print_dual("Status            : OK", f_out)
    print_dual(f"Barrier           : {barrier:.4f} eV", f_out)
    print_dual(f"Reaction energy   : {dE:.4f} eV", f_out)
    if report_path:
        print_dual(f"Report            : {report_path}", f_out)

    if f_out:
        f_out.close()

    print("\n" + "-" * 60)
    print(color_text("ML NEB complete -- a fast heuristic, not a substitute for a DFT-level "
                      "barrier.\n", 'bold'))


if __name__ == "__main__":
    main()
