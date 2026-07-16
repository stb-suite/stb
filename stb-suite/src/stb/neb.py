#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.0.0"

import os
import sys
import argparse
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from stb.core import structure_io
from stb.core.cli import color_text, show_intro, print_dual
from stb.core.pseudopotentials import resolve_pseudo_source, link_pseudo
from stb.core.deps import require_mace
from stb.core.calc_directives import force_single_point

REPORT_FILE = "neb_setup.txt"


def check_composition_match(initial_structure, final_structure):
    """Hard-validates the two endpoint FdfStructures have identical
    composition (same element counts) -- NEB has no notion of a reaction
    that creates or destroys atoms; that's a different problem
    (stb-cohesive-style thermodynamics), not a reaction path. Compares
    FdfStructure.atoms directly (before any pymatgen conversion) since
    that's the cheapest, least ambiguous source of "what atoms are
    actually here". Exits with a clear per-species breakdown if they
    don't match.
    """
    initial_counts = Counter(sym for sym, _ in initial_structure.atoms)
    final_counts = Counter(sym for sym, _ in final_structure.atoms)
    if initial_counts == final_counts:
        return
    all_symbols = sorted(set(initial_counts) | set(final_counts))
    breakdown = ", ".join(
        f"{sym}: {initial_counts.get(sym, 0)} vs {final_counts.get(sym, 0)}"
        for sym in all_symbols
        if initial_counts.get(sym, 0) != final_counts.get(sym, 0)
    )
    print(color_text(
        f"[ERROR] Initial and final structures have different composition ({breakdown}) -- "
        "NEB requires the exact same atoms on both endpoints (a reaction that creates or "
        "destroys atoms is not a reaction path in this sense).", 'red'))
    sys.exit(1)


def wrap_into_cell(pmg_structure):
    """Returns a copy of `pmg_structure` with every fractional coordinate
    wrapped into [0, 1) -- a canonicalization, not a physical change (same
    atoms, same lattice, just relabeling which periodic image each site is
    nominally reported in). Required before linear_interpolate_images'
    deliberate pbc=False interpolation below: two independently-written
    endpoint .fdf files can describe the SAME physical position with
    different wrapping conventions (e.g. frac z=0.95 in one file, z=-0.05
    in the other -- 1.0 apart nominally, but literally the same point once
    wrapped). Interpolating between the UNWRAPPED coordinates would
    silently produce a spuriously huge, wrong path instead of the
    near-zero-length one the two structures actually represent --
    reproduced live: frac z 0.95 -> -0.05 in a 10 Ang cell interpolated
    (unwrapped) as a full 10 Ang traversal through z=7/4.5/2/-0.5, instead
    of a static point once both were wrapped first.
    """
    return Structure(pmg_structure.lattice, pmg_structure.species,
                      pmg_structure.frac_coords % 1.0, coords_are_cartesian=False)


def resolve_lattice_mismatch(initial_pmg, final_pmg, tol=1e-3):
    """Returns a copy of `final_pmg` rebuilt on the INITIAL structure's
    lattice (species/frac_coords kept from final_pmg). Prints a
    [WARNING] if the two lattices differ by more than `tol` Ang in any
    component -- not just stylistic: ase.mep.neb.NEB/idpp_interpolate both
    raise NotImplementedError on any per-image cell mismatch (no
    variable-cell NEB support in ASE), and pymatgen's Structure.interpolate
    raises ValueError on unequal lattices unless interpolate_lattices=True
    -- so a single fixed lattice for the whole band is the only thing
    either downstream library can actually run, confirmed live against
    both APIs. Below `tol`, no message (ordinary floating-point/rounding
    noise between two independently-relaxed endpoint calculations).
    """
    initial_matrix = np.array(initial_pmg.lattice.matrix)
    final_matrix = np.array(final_pmg.lattice.matrix)
    max_diff = float(np.abs(initial_matrix - final_matrix).max())
    if max_diff > tol:
        print(color_text(
            f"[WARNING] Initial and final structures have different lattices (largest "
            f"component difference: {max_diff:.4f} Ang) -- ase.mep.neb.NEB and pymatgen's "
            "interpolation both require every image to share one exact cell (no variable-cell "
            "NEB support in ASE). Adopting the INITIAL structure's lattice for the whole band; "
            "the final structure's atomic (fractional) positions are kept, its own lattice is "
            "discarded.", 'yellow'))
    return Structure(initial_pmg.lattice, final_pmg.species, final_pmg.frac_coords,
                      coords_are_cartesian=False)


def linear_interpolate_images(initial_pmg, final_pmg_matched, n_images, autosort_tol=0.5):
    """Linearly interpolates `n_images` pymatgen Structures (endpoints
    included) between the two endpoints. `autosort_tol` fixes atom-index
    correspondence between the two endpoint .fdf files (a classic NEB
    pitfall if they don't list atoms in the same order) -- pymatgen's own
    documented "usually works well" value is 0.5 Ang, exposed here as
    --autosort-tol for the rare case a very close-packed structure needs a
    tighter value.

    pbc=False is deliberate, not pymatgen's own default (pbc=True):
    Structure.interpolate's default takes the MINIMUM-IMAGE path in
    fractional coordinates, which can silently interpolate "the short way
    around" through a periodic boundary instead of the direct path between
    the two endpoint coordinates as actually given -- reproduced live
    (frac x: 0.1 -> 0.9 interpolated as 0.1 -> -0.1, i.e. through 0/1,
    instead of 0.1 -> 0.9 directly) -- physically wrong for a reaction path
    whose two endpoints are specific, intentional configurations, not just
    "whichever periodic image happens to be closer".

    Precondition: both `initial_pmg` and `final_pmg_matched` must already
    be wrapped into [0, 1) (see wrap_into_cell(), called on both endpoints
    in main() before this) -- with pbc=False, this function takes the
    fractional coordinates completely literally, so two endpoints that
    represent the same physical position via different (unwrapped)
    conventions would otherwise interpolate a spuriously huge path instead
    of a near-zero-length one.
    """
    return initial_pmg.interpolate(final_pmg_matched, nimages=n_images - 1,
                                    interpolate_lattices=False, pbc=False,
                                    autosort_tol=autosort_tol)


def idpp_refine_images(pmg_images, fmax=0.1, steps=100):
    """Refines the interior images of `pmg_images` (endpoints untouched)
    with ASE's Image Dependent Pair Potential method -- a purely
    classical/geometric smoothing (no MACE/chemistry needed) that usually
    gives a much better starting guess than plain linear interpolation,
    especially when the direct path would otherwise pass unphysically
    close to another atom. Builds an explicit ase.mep.neb.NEB object first
    (method="improvedtangent") and hands THAT to idpp_interpolate rather
    than a plain list -- idpp_interpolate accepts either, but a plain list
    makes it build its own NEB(images) internally with ASE's default
    method, which emits a UserWarning about the aseneb -> improvedtangent
    default change; passing our own object avoids the spurious warning
    without changing behavior. traj=None/log=None so this doesn't litter
    the working directory with idpp.traj/idpp.log.
    """
    from ase.mep.neb import NEB, idpp_interpolate
    ase_images = [AseAtomsAdaptor.get_atoms(s) for s in pmg_images]
    neb = NEB(ase_images, method="improvedtangent")
    idpp_interpolate(neb, traj=None, log=None, fmax=fmax, steps=steps)
    return [AseAtomsAdaptor.get_structure(a) for a in ase_images]


def cumulative_reaction_coordinates(pmg_images):
    """[0.0, ...]: cumulative Cartesian L2-norm displacement from image_00,
    one numpy flatten+norm per consecutive image pair. No periodic-image
    (minimum-distance) wrapping -- every image shares one lattice and was
    built by direct (pbc=False) interpolation, so a naive Cartesian
    distance between consecutive images is the physically meaningful one,
    unlike a general min-image pairwise-distance calculation.
    """
    coords = [0.0]
    for prev, curr in zip(pmg_images[:-1], pmg_images[1:]):
        delta = curr.cart_coords - prev.cart_coords
        coords.append(coords[-1] + float(np.linalg.norm(delta)))
    return coords


def compute_frozen_indices(pmg_images, threshold=0.3):
    """Returns the sorted atom indices whose Cartesian displacement between
    the FIRST and LAST image (the two endpoints exactly as returned by the
    interpolation -- never touched by IDPP or --ml-neb, so this reflects
    the user's own initial/final structures, index-matched by
    linear_interpolate_images' autosort_tol) is below `threshold` Ang.
    These atoms are essentially spectators to the reaction and safe to
    freeze during --ml-neb's climbing-image relaxation -- both for speed
    (fewer degrees of freedom for the optimizer) and physical correctness
    (an atom far from the reaction site shouldn't spuriously wander during
    a fast/loose ML relax). Self-regulating: if every atom moves a lot
    between endpoints (e.g. a small-molecule isomerization, no real
    "substrate"), the returned list is simply empty and nothing is frozen.
    """
    initial, final = pmg_images[0], pmg_images[-1]
    frozen = []
    for i in range(len(initial)):
        delta = final.cart_coords[i] - initial.cart_coords[i]
        if np.linalg.norm(delta) < threshold:
            frozen.append(i)
    return frozen


def check_path_quality(pmg_images, reaction_coords, f_out):
    """Prints (and persists) advisory warnings about the interpolated
    path's geometry, based on the step size between consecutive images:
    a near-zero step (two images landing on almost the same geometry --
    wasted SIESTA calculations, or a sign --n-images is larger than the
    path needs) and an unusually large step relative to the mean (often a
    sign of a bad atom correspondence between the two endpoint .fdf files
    -- see --autosort-tol -- rather than a genuinely non-uniform path).
    Silent when the path looks geometrically reasonable, same "advisory
    only, don't clutter a clean run" convention as
    core/siesta_log.py::report_quality_diagnostics.
    """
    steps = [reaction_coords[i + 1] - reaction_coords[i] for i in range(len(reaction_coords) - 1)]
    if not steps:
        return
    mean_step = sum(steps) / len(steps)
    if mean_step <= 0:
        return
    for i, step in enumerate(steps):
        if step < 0.05:
            print_dual(color_text(
                f"  [WARNING] image_{i:02d} and image_{i + 1:02d} are nearly identical "
                f"(step = {step:.4f} Ang) -- consider fewer --n-images, or check that "
                "--initial/--final aren't already the same structure.", 'yellow'), f_out)
        elif step > 3.0 * mean_step:
            print_dual(color_text(
                f"  [WARNING] Unusually large step between image_{i:02d} and image_{i + 1:02d} "
                f"({step:.4f} Ang vs. a mean step of {mean_step:.4f} Ang) -- may indicate a bad "
                "atom correspondence between --initial/--final (see --autosort-tol) rather than "
                "a genuinely non-uniform reaction path.", 'yellow'), f_out)


def write_path_trajectory(pmg_images, out_path):
    """Writes every image (in order) to a single multi-frame extended-XYZ
    file, viewable directly in VESTA/OVITO/ASE-GUI -- a structural
    counterpart to write_ml_preview_plot's energy-only preview, always
    generated (cheap, no extra dependency: ASE writes a Sequence[Atoms] to
    one multi-frame file natively).
    """
    from ase.io import write as ase_write
    ase_images = [AseAtomsAdaptor.get_atoms(s) for s in pmg_images]
    ase_write(out_path, ase_images)


def write_image_folder(out_dir, pmg_structure, calc_text, species_meta, pp_path):
    """Writes structure.fdf + calc.fdf + linked pseudopotentials for one
    image_NN/ folder -- same shape as adsorb.py's write_reference_folder.
    """
    os.makedirs(out_dir, exist_ok=True)
    fdf_structure = structure_io.from_pymatgen(pmg_structure, species_meta=species_meta,
                                                coord_format="fractional")
    structure_io.write_fdf(fdf_structure, os.path.join(out_dir, "structure.fdf"))
    with open(os.path.join(out_dir, "calc.fdf"), "w") as f:
        f.write(calc_text)
    symbols = {site.specie.symbol for site in pmg_structure}
    for sym in sorted(symbols):
        link_pseudo(pp_path, sym, out_dir)


def write_ml_preview_plot(ase_images, out_path):
    """Energy-profile preview of the MACE-MP-0-relaxed band, via ASE's own
    NEBTools.plot_band() -- only called when --ml-neb converged a band.
    """
    from ase.mep.neb import NEBTools
    fig = NEBTools(ase_images).plot_band()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Prepares SIESTA folders for a NEB (Nudged Elastic Band) "
        "reaction-path study: interpolates images between two already-relaxed endpoint "
        "structures, optionally refining the path with MACE-MP-0 before writing single-point "
        "SIESTA folders per image.", 'bold')}""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s -i initial.fdf -f final.fdf -c calc.fdf -n 7\n"
               "  %(prog)s -i initial.fdf -f final.fdf -c calc.fdf --idpp\n"
               "  %(prog)s -i initial.fdf -f final.fdf -c calc.fdf --ml-neb\n"
               "  %(prog)s -i initial.fdf -f final.fdf -c calc.fdf --ml-neb --idpp --ml-k 0.2\n"
    )

    parser.add_argument("-i", "--initial", type=str, required=True,
                         help="Initial (already-relaxed) endpoint structure.fdf.")
    parser.add_argument("-f", "--final", type=str, required=True,
                         help="Final (already-relaxed) endpoint structure.fdf -- same "
                              "composition as --initial.")
    parser.add_argument("-c", "--calc", type=str, required=True,
                         help="calc.fdf template (kgrid, basis, XC, %%include structure.fdf, "
                              "etc.) -- copied into every image_NN/ folder with MD.TypeOfRun/"
                              "MD.NumCGsteps forced to a single-point evaluation.")
    parser.add_argument("-p", "--pseudo-dir", type=str, default="",
                         help="Pseudopotentials source (optional): a bundled bank or a folder path.")

    parser.add_argument("-n", "--n-images", type=int, default=7,
                         help="Total images along the band, endpoints included (default: 7).")
    parser.add_argument("--autosort-tol", type=float, default=0.5,
                         help="Atom-correspondence tolerance (Ang) for linear interpolation "
                              "(default: 0.5, pymatgen's own suggested value).")

    parser.add_argument("--idpp", action="store_true",
                         help="Refine the interior images with ASE's IDPP method after linear "
                              "interpolation -- a better initial guess than plain linear, no "
                              "MACE needed.")
    parser.add_argument("--idpp-fmax", type=float, default=0.1)
    parser.add_argument("--idpp-steps", type=int, default=100)

    parser.add_argument("--ml-neb", action="store_true",
                         help="Run a real climbing-image NEB on MACE-MP-0, live, before writing "
                              "SIESTA folders -- a fast pre-relaxation of the whole band's shape "
                              "(not a replacement for stb-nebAnalysis's DFT-level energies). "
                              "Needs the optional 'ml' extra.")
    parser.add_argument("--ml-model", choices=["small", "medium", "large"], default="small")
    parser.add_argument("--ml-device", choices=["cpu", "cuda"], default="cpu")
    parser.add_argument("--ml-fmax", type=float, default=0.05)
    parser.add_argument("--ml-k", type=float, default=0.1,
                         help="NEB spring constant, eV/Ang^2 (default: 0.1, ASE's own default).")
    parser.add_argument("--ml-max-steps", type=int, default=200)
    parser.add_argument("--ml-freeze-substrate", dest="ml_freeze_substrate", action="store_true",
                         default=True,
                         help="With --ml-neb: freeze atoms whose position barely differs between "
                              "--initial and --final (default: ON) -- fewer degrees of freedom "
                              "for the ML optimizer and avoids spurious drift of atoms that are "
                              "spectators to the reaction. See --ml-freeze-threshold.")
    parser.add_argument("--no-ml-freeze-substrate", dest="ml_freeze_substrate", action="store_false",
                         help="Let every atom relax during --ml-neb, even ones that don't move "
                              "between --initial and --final.")
    parser.add_argument("--ml-freeze-threshold", type=float, default=0.3,
                         help="Displacement threshold in Ang (default: 0.3) below which an atom "
                              "is considered a spectator and frozen by --ml-freeze-substrate.")
    parser.add_argument("--ml-prerelax-endpoints", action="store_true",
                         help="Relax both endpoints' positions (cell fixed) with MACE-MP-0 "
                              "before interpolating -- independent of --ml-neb (a cheap safety "
                              "net when you're not fully certain the endpoints are already "
                              "relaxed, without committing to a full band NEB). Needs the "
                              "optional 'ml' extra.")

    parser.add_argument("-O", "--output-dir", type=str, default=".",
                         help="Root directory (default: current directory) for every image_NN/.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-neb {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    if args.n_images < 3:
        parser.error("-n/--n-images must be >= 3 (2 endpoints + at least 1 interior image).")

    print("\n" + color_text("Prepare a NEB reaction-path study:", 'bold'))
    print("-" * 60)

    if not os.path.exists(args.initial):
        print(color_text(f"[ERROR] Initial structure file '{args.initial}' not found.", 'red'))
        sys.exit(1)
    if not os.path.exists(args.final):
        print(color_text(f"[ERROR] Final structure file '{args.final}' not found.", 'red'))
        sys.exit(1)
    if not os.path.exists(args.calc):
        print(color_text(f"[ERROR] Calc file '{args.calc}' not found.", 'red'))
        sys.exit(1)

    if args.pseudo_dir:
        try:
            args.pseudo_dir = resolve_pseudo_source(args.pseudo_dir)
        except ValueError as e:
            print(color_text(f"[ERROR] {e}", 'red'))
            sys.exit(1)

    initial_structure = structure_io.read_fdf(args.initial)
    final_structure = structure_io.read_fdf(args.final)
    check_composition_match(initial_structure, final_structure)
    species_meta = structure_io.species_dict(initial_structure)

    initial_pmg = structure_io.to_pymatgen(initial_structure)
    final_pmg = structure_io.to_pymatgen(final_structure)

    with open(args.calc) as f:
        calc_text = f.read()
    single_point_calc_text = force_single_point(calc_text)

    output_root = args.output_dir
    os.makedirs(output_root, exist_ok=True)
    report_path = os.path.join(output_root, REPORT_FILE)

    with open(report_path, 'w') as f_out:
        print_dual(f"\n{color_text('[0] RUN METADATA', 'bold')}", f_out)
        print_dual(f"  Initial structure : {args.initial}", f_out)
        print_dual(f"  Final structure   : {args.final}", f_out)
        print_dual(f"  Calc template     : {args.calc}", f_out)
        print_dual(f"  Pseudo dir        : {args.pseudo_dir or '(none)'}", f_out)
        print_dual(f"  Output dir        : {output_root}", f_out)
        print_dual(f"  N images          : {args.n_images}", f_out)
        print_dual(f"  IDPP refinement   : {'yes' if args.idpp else 'no'}", f_out)
        print_dual(f"  ML-NEB (MACE-MP-0): {'yes' if args.ml_neb else 'no'}", f_out)
        if args.ml_neb:
            print_dual(f"  ML-NEB freeze substrate: {'yes' if args.ml_freeze_substrate else 'no'}"
                        + (f" (threshold {args.ml_freeze_threshold} Ang)"
                           if args.ml_freeze_substrate else ""), f_out)
        print_dual(f"  ML pre-relax endpoints: {'yes' if args.ml_prerelax_endpoints else 'no'}", f_out)

        print_dual(f"\n{color_text('[1] INTERPOLATION', 'bold')}", f_out)

        if args.ml_prerelax_endpoints:
            require_mace()
            from stb.core import mace_relax
            print_dual(f"  {color_text('ML pre-relax:', 'cyan')} relaxing both endpoints "
                        "(positions only) with MACE-MP-0 ...", f_out)
            calc_mace_endpoints = mace_relax.get_calculator(model=args.ml_model, device=args.ml_device)
            for label, pmg in (("initial", initial_pmg), ("final", final_pmg)):
                ase_atoms = AseAtomsAdaptor.get_atoms(pmg)
                converged, steps = mace_relax.relax(ase_atoms, calc_mace_endpoints,
                                                     fmax=args.ml_fmax, max_steps=200)
                relaxed_pmg = AseAtomsAdaptor.get_structure(ase_atoms)
                if label == "initial":
                    initial_pmg = relaxed_pmg
                else:
                    final_pmg = relaxed_pmg
                print_dual(f"  {'Converged' if converged else 'Hit step cap, not fully converged'} "
                            f"({label}) after {steps} step(s).", f_out)

        initial_pmg = wrap_into_cell(initial_pmg)
        final_pmg = wrap_into_cell(final_pmg)
        final_pmg_matched = resolve_lattice_mismatch(initial_pmg, final_pmg)

        pmg_images = linear_interpolate_images(initial_pmg, final_pmg_matched, args.n_images,
                                                autosort_tol=args.autosort_tol)
        print_dual(f"  Linear interpolation: {len(pmg_images)} images.", f_out)

        if args.idpp:
            pmg_images = idpp_refine_images(pmg_images, fmax=args.idpp_fmax, steps=args.idpp_steps)
            print_dual("  IDPP refinement: interior images refined with ASE's Image Dependent "
                        "Pair Potential method.", f_out)

        ml_neb_used = False
        if args.ml_neb:
            require_mace()
            from stb.core import mace_relax
            print_dual(f"  {color_text('ML-NEB:', 'cyan')} running a real climbing-image NEB "
                        "on MACE-MP-0 ...", f_out)
            ase_images = [AseAtomsAdaptor.get_atoms(s) for s in pmg_images]
            if args.ml_freeze_substrate:
                frozen_indices = compute_frozen_indices(pmg_images, threshold=args.ml_freeze_threshold)
                if frozen_indices:
                    from ase.constraints import FixAtoms
                    for atoms in ase_images:
                        atoms.set_constraint(FixAtoms(indices=frozen_indices))
                print_dual(f"  Freezing {len(frozen_indices)}/{len(ase_images[0])} atom(s) with "
                            f"< {args.ml_freeze_threshold} Ang displacement between endpoints.", f_out)
            calc_mace_neb = mace_relax.get_calculator(model=args.ml_model, device=args.ml_device)
            converged, s1, s2, energies = mace_relax.relax_neb(
                ase_images, calc_mace_neb, k=args.ml_k, fmax=args.ml_fmax,
                max_steps=args.ml_max_steps)
            pmg_images = [AseAtomsAdaptor.get_structure(a) for a in ase_images]
            ml_neb_used = True
            print_dual(f"  {'Converged' if converged else 'Hit step cap, not fully converged'} "
                        f"after {s1} (stage 1) + {s2} (stage 2, climbing) step(s).", f_out)

            from ase.mep.neb import NEBTools
            barrier, dE = NEBTools(ase_images).get_barrier(fit=True)
            print_dual(f"  ML barrier estimate: {barrier:.4f} eV (forward), reaction energy: "
                        f"{dE:.4f} eV (fitted spline over the relaxed band).", f_out)
            print_dual(color_text(
                "  Note: an ML-level estimate from a fast surrogate potential, not a DFT "
                "adsorption/reaction barrier -- stb-nebAnalysis will report the DFT-level "
                "value from each image_NN/'s SIESTA single-point once you run them.", 'yellow'), f_out)

            preview_path = os.path.join(output_root, "neb_ml_preview.png")
            write_ml_preview_plot(ase_images, preview_path)
            print_dual(f"  {color_text('[Saved]', 'cyan')} {preview_path}", f_out)

        reaction_coords = cumulative_reaction_coordinates(pmg_images)
        check_path_quality(pmg_images, reaction_coords, f_out)

        trajectory_path = os.path.join(output_root, "neb_path.xyz")
        write_path_trajectory(pmg_images, trajectory_path)
        print_dual(f"  {color_text('[Saved]', 'cyan')} {trajectory_path} (all images, viewable in "
                    "VESTA/OVITO/ASE-GUI)", f_out)

        print_dual(f"\n{color_text('[2] IMAGE FOLDERS', 'bold')}", f_out)
        report_rows = []  # (label, index, reaction_coord, dir)
        for i, pmg_image in enumerate(pmg_images):
            label = f"image_{i:02d}"
            image_dir = os.path.join(output_root, label)
            write_image_folder(image_dir, pmg_image, single_point_calc_text, species_meta,
                                args.pseudo_dir)
            print_dual(f"  {color_text('[OK]', 'green')} {image_dir}", f_out)
            report_rows.append((label, i, reaction_coords[i], image_dir))

        print_dual(f"\n{color_text('[3] SUMMARY', 'bold')}", f_out)
        print_dual(f"  {len(pmg_images)} image folder(s) written under '{output_root}'.", f_out)
        print_dual("  Run SIESTA (single-point: MD.NumCGsteps forced to 0) in every image_NN/ "
                    "folder, then use stb-nebAnalysis.", f_out)

        f_out.write(f"\n# ML_NEB_USED: {'yes' if ml_neb_used else 'no'}\n")
        f_out.write("# IMAGE_TABLE -- parsed by stb-nebAnalysis, do not reorder the columns\n")
        f_out.write(f"# {'label':<14}{'index':<8}{'reaction_coord':<18}{'dir'}\n")
        for label, index, reaction_coord, image_dir in report_rows:
            f_out.write(f"{label:<16}{index:<8}{reaction_coord:<18.6f}{image_dir}\n")

    print(f"\n{color_text('Success:', 'green')} {len(pmg_images)} image folder(s) written under "
          f"'{output_root}'.")
    print(f"Full report: {report_path}")


if __name__ == "__main__":
    main()
