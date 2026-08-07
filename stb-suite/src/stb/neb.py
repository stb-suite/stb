#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.0.0"

import os
import sys
import json
import argparse
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from stb.core import structure_io
from stb.core import neb_manifest
from stb.core.cli import color_text, show_intro, print_dual
from stb.core.pseudopotentials import resolve_pseudo_source, copy_pseudo
from stb.core.deps import require_mace
from stb.core.calc_directives import force_single_point

REPORT_FILE = "neb_setup.txt"
MACE_RESULT_FILE = "neb_mace_result.json"

_MODE_DESCRIPTIONS = {
    1: "100% MACE-MP-0, JSON result, no SIESTA",
    2: "100% MACE-MP-0, then 1 single-point SIESTA per image",
    3: "MACE-MP-0 + a few real-DFT NEB refinement cycles (default)",
    4: "100% real-DFT NEB from a plain interpolated path",
}


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


def resolve_manifest_pair(initial_path, final_path, initial_structure, final_structure):
    """Best-effort: if BOTH --initial/--final are directories containing a
    neb_manifest.json (written by stb-nebSites), loads them, cross-checks
    them against EACH OTHER (neb_manifest.validate_manifest_pair) AND
    against the structures actually read back by read_relaxed_or_input --
    the second check catches a stale manifest (structure.fdf/siesta.XV
    regenerated or hand-edited since stb-nebSites wrote it), which a
    manifest-vs-manifest-only comparison can't see.

    Returns (manifest_initial, manifest_final, pair_ids_match) when both
    sides have a manifest and everything checks out. Returns None if
    EITHER side has no manifest at all -- not an error, just "no proof
    available" (an older stb-nebSites run predating this feature, or a
    hand-built --initial/--final) -- falls back to the existing distance
    -based --autosort-tol behavior. Exits (clean [ERROR], same convention
    as check_composition_match) for every case where a manifest genuinely
    IS present but can't be trusted.
    """
    paths = {"--initial": (initial_path, initial_structure), "--final": (final_path, final_structure)}
    if not all(os.path.isdir(p) for p, _ in paths.values()):
        return None
    manifest_paths = {label: os.path.join(p, neb_manifest.MANIFEST_FILENAME)
                       for label, (p, _) in paths.items()}
    if not all(os.path.isfile(mp) for mp in manifest_paths.values()):
        return None

    manifests = {}
    for label, mp in manifest_paths.items():
        try:
            manifests[label] = neb_manifest.load_manifest(mp)
        except ValueError as e:
            print(color_text(f"[ERROR] {e}", 'red'))
            sys.exit(1)

    for label, (_p, structure) in paths.items():
        actual_sequence = [sym for sym, _ in structure.atoms]
        if actual_sequence != manifests[label]["species_sequence"]:
            print(color_text(
                f"[ERROR] {label}'s read-back structure does not match its own "
                f"'{neb_manifest.MANIFEST_FILENAME}' -- the folder's structure.fdf/siesta.XV was "
                "likely regenerated or hand-edited after stb-nebSites wrote the manifest. Delete "
                "the stale manifest (or re-run stb-nebSites) before retrying.", 'red'))
            sys.exit(1)

    try:
        pair_ids_match = neb_manifest.validate_manifest_pair(manifests["--initial"], manifests["--final"])
    except ValueError as e:
        print(color_text(f"[ERROR] {e}", 'red'))
        sys.exit(1)
    return manifests["--initial"], manifests["--final"], pair_ids_match


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
    interpolation -- never touched by IDPP or the MACE path-shaping stage,
    so this reflects the user's own initial/final structures, index-matched
    by linear_interpolate_images' autosort_tol) is below `threshold` Ang.
    These atoms are essentially spectators to the reaction and safe to
    freeze during the MACE climbing-image relaxation (modes 1/2/3) -- both for speed
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


def check_endpoint_displacement(initial_pmg, final_pmg, f_out, threshold, idpp_used):
    """Advisory diagnostic: reports the max/mean per-atom Cartesian
    displacement between --initial and --final, and -- only when a
    meaningful FRACTION of the whole structure moved, not just the
    reacting atom(s) itself -- recommends --idpp. Deliberately gated on
    BOTH an absolute minimum (>= 3 atoms) AND a relative fraction (> 10%
    of the structure) so an ordinary single- or few-atom reaction hop
    (e.g. this suite's own initial.fdf/final.fdf H-diffusion fixture:
    1/9 atoms move ~2.5 Ang, the reacting H itself -- a completely normal,
    expected displacement for any real reaction path) never triggers this;
    only a broad multi-atom rearrangement does. This exact shape (many
    atoms, not just the reactive site, displaced beyond `threshold`) is
    what was live-verified to correlate with an unphysical ~47 eV barrier
    from plain linear interpolation with no --idpp on real user data (16/33
    atoms > 0.3 Ang) -- independent of, and in addition to, atom-order
    correctness (this can fire even when the correspondence is manifest
    -PROVEN correct). `threshold` reuses --ml-freeze-threshold's own
    already-documented default (0.3 Ang) rather than inventing a new flag.
    """
    deltas = np.linalg.norm(final_pmg.cart_coords - initial_pmg.cart_coords, axis=1)
    n_moved = int(np.sum(deltas > threshold))
    print_dual(f"  Endpoint displacement: max {deltas.max():.3f} Ang, mean {deltas.mean():.3f} Ang, "
                f"{n_moved}/{len(deltas)} atom(s) moved > {threshold} Ang.", f_out)
    if n_moved >= 3 and (n_moved / len(deltas)) > 0.1 and not idpp_used:
        print_dual(color_text(
            f"  [WARNING] {n_moved}/{len(deltas)} atoms moved more than {threshold} Ang between "
            "--initial and --final -- a broad rearrangement, not just the reacting atom(s). Plain "
            "linear interpolation over this large a change can pass unphysically close to another "
            "atom in an intermediate image. Consider --idpp (ASE's Image Dependent Pair Potential "
            "refinement) for a substantially better starting path.", 'yellow'), f_out)


def check_path_quality(pmg_images, reaction_coords, f_out, manifest_proven=False):
    """Prints (and persists) advisory warnings about the interpolated
    path's geometry, based on the step size between consecutive images:
    a near-zero step (two images landing on almost the same geometry --
    wasted SIESTA calculations, or a sign --n-images is larger than the
    path needs) and an unusually large step relative to the mean (often a
    sign of a bad atom correspondence between the two endpoint .fdf files
    -- see --autosort-tol -- rather than a genuinely non-uniform path,
    unless `manifest_proven` is True -- see below). Silent when the path
    looks geometrically reasonable, same "advisory only, don't clutter a
    clean run" convention as core/siesta_log.py::report_quality_diagnostics.
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
            hint = ("may indicate a genuinely non-uniform reaction path (atom correspondence is "
                    "already PROVEN via neb_manifest.json, not a distance-matching guess)"
                    if manifest_proven else
                    "may indicate a bad atom correspondence between --initial/--final (see "
                    "--autosort-tol) rather than a genuinely non-uniform reaction path")
            print_dual(color_text(
                f"  [WARNING] Unusually large step between image_{i:02d} and image_{i + 1:02d} "
                f"({step:.4f} Ang vs. a mean step of {mean_step:.4f} Ang) -- {hint}.", 'yellow'), f_out)


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
    """Writes structure.fdf + calc.fdf + copied pseudopotentials for one
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
        copy_pseudo(pp_path, sym, out_dir)


def write_ml_preview_plot(ase_images, out_path):
    """Energy-profile preview of the MACE-MP-0-relaxed band, via ASE's own
    NEBTools.plot_band() -- only called after the MACE path-shaping stage
    (modes 1/2/3) has run.
    """
    from ase.mep.neb import NEBTools
    fig = NEBTools(ase_images).plot_band()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def _write_mode1_json(output_root, ase_images, reaction_coords, converged, barrier, dE, args, f_out):
    """Mode 1's entire output: a single JSON (MACE_RESULT_FILE) with every
    image's energy/positions/symbols plus the fitted barrier/reaction
    energy -- stb-nebAnalysis reads this directly, no calc.out anywhere,
    since mode 1 never touches SIESTA at all. `ase_images` already carry
    their MACE-computed energies (get_potential_energy() below reads the
    same cached result NEBTools.get_barrier() itself just used).
    """
    backward = barrier - dE
    payload = {
        "mode": 1,
        "k": args.ml_k,
        "model": args.ml_model,
        "converged": bool(converged),
        "barrier_forward_eV": float(barrier),
        "barrier_backward_eV": float(backward),
        "reaction_energy_eV": float(dE),
        "images": [
            {
                "index": i,
                "label": f"image_{i:02d}",
                "reaction_coord": reaction_coords[i],
                "energy_eV": float(atoms.get_potential_energy()),
                "symbols": atoms.get_chemical_symbols(),
                "positions": atoms.positions.tolist(),
                "cell": np.asarray(atoms.cell).tolist(),
            }
            for i, atoms in enumerate(ase_images)
        ],
    }
    json_path = os.path.join(output_root, MACE_RESULT_FILE)
    with open(json_path, "w") as f:
        json.dump(payload, f, indent=2)

    print_dual(f"\n{color_text('[2] MACE RESULT (JSON)', 'bold')}", f_out)
    print_dual(f"  {color_text('[Saved]', 'cyan')} {json_path} ({len(ase_images)} image(s), "
                f"barrier {barrier:.4f} eV forward / {backward:.4f} eV backward, reaction "
                f"energy {dE:.4f} eV)", f_out)

    print_dual(f"\n{color_text('[3] SUMMARY', 'bold')}", f_out)
    print_dual("  Mode 1: no SIESTA folders written -- the path was fully converged in "
                "MACE-MP-0. Run stb-nebAnalysis on this directory to read the JSON directly.",
                f_out)

    f_out.write("\n# MODE: 1\n")
    f_out.write(f"# MACE_RESULT_FILE: {MACE_RESULT_FILE}\n")


_CLIMB_AFTER_SUGGESTION = {3: 0, 4: 5}


def _print_cluster_submission_snippet(f_out, mode):
    """Prints the sub.sh-ready loop (modes 3/4 only): run SIESTA in every
    image_* of the current cycle, call stb-nebCycle (a separate, CLI-only
    tool -- deliberately not part of the interactive stb-suite menu) to
    take one real-DFT NEB step, and stop once it writes the NEB_CONVERGED
    sentinel. --climb-after differs by mode: mode 3 already handed off a
    MACE-shaped (already climbing-image-converged) path, so climbing from
    cycle 0 is safe; mode 4 starts from a raw interpolated path, where
    climbing too early risks locking onto the wrong image as the saddle
    (same reasoning core/mace_relax.py::relax_neb's own two-stage
    climb=False-then-True approach exists for, here spread across
    separate cluster-queued cycles instead of one live process).
    """
    climb_after = _CLIMB_AFTER_SUGGESTION[mode]
    print_dual(f"\n{color_text('[4] CLUSTER SUBMISSION', 'bold')}", f_out)
    print_dual("  stb-nebCycle is a separate CLI-only tool (not in the interactive stb-suite "
                "menu) that reads the latest finished cycle_NN/ and takes one real-DFT NEB "
                "step -- run it in a loop from your own submission script, alternating with "
                "real SIESTA runs:", f_out)
    snippet = f"""
for cycle in $(seq -w 0 29); do
    cdir="cycle_${{cycle}}"
    [ -d "$cdir" ] || break
    for img in "$cdir"/image_*; do
        (cd "$img" && siesta < calc.fdf > calc.out)   # adjust to your real queue wrapper
    done
    stb-nebCycle --dir . --fmax 0.05 --climb-after {climb_after}
    [ -f NEB_CONVERGED ] && {{ echo "NEB converged at cycle $cycle"; break; }}
done
"""
    print_dual(snippet, f_out)


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
               "  %(prog)s -i initial.fdf -f final.fdf -c calc.fdf --mode 1\n"
               "  %(prog)s -i initial.fdf -f final.fdf -c calc.fdf --mode 2 --idpp --ml-k 0.2\n"
               "  %(prog)s -i initial.fdf -f final.fdf -c calc.fdf --mode 4\n"
    )

    parser.add_argument("-i", "--initial", type=str, required=True,
                         help="Initial (already-relaxed) endpoint: a structure.fdf file, or a "
                              "directory containing one (e.g. a stb-nebSites site_A/ folder) -- "
                              "a directory's own finished siesta.XV, if present, is preferred "
                              "over its structure.fdf guess.")
    parser.add_argument("-f", "--final", type=str, required=True,
                         help="Final (already-relaxed) endpoint -- same file-or-directory "
                              "convention as --initial, same composition.")
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
                              "(default: 0.5, pymatgen's own suggested value for two endpoints "
                              "that may have a different atom order). Pass 0 to skip distance "
                              "-based matching entirely and use --initial/--final's atom order "
                              "directly -- the right choice (and more robust for a small/"
                              "densely-packed cell) when both endpoints already share a "
                              "guaranteed matching order, e.g. a stb-nebSites site_A/site_B pair. "
                              "If matching fails at this tolerance, stb-neb automatically retries "
                              "once with --autosort-tol 0 (a [WARNING], not an [ERROR]) before "
                              "giving up -- this is usually already the right correspondence, so "
                              "most runs never need this flag set explicitly at all.")

    parser.add_argument("--idpp", action="store_true",
                         help="Refine the interior images with ASE's IDPP method after linear "
                              "interpolation -- a better initial guess than plain linear, no "
                              "MACE needed.")
    parser.add_argument("--idpp-fmax", type=float, default=0.1)
    parser.add_argument("--idpp-steps", type=int, default=100)

    parser.add_argument("--mode", type=int, choices=[1, 2, 3, 4], default=3,
                         help="How to build the path and what to write (default: 3):\n"
                              "  1 = 100%% MACE-MP-0 (full climbing-image NEB) -> a single "
                              "neb_mace_result.json for stb-nebAnalysis, no SIESTA folders at all.\n"
                              "  2 = 100%% MACE-MP-0, then ONE single-point SIESTA per image "
                              "-> image_NN/ (this is the old --ml-neb behavior).\n"
                              "  3 = MACE-MP-0 shapes the path, then a FEW real-DFT NEB "
                              "refinement cycles (via the separate stb-nebCycle CLI tool, run "
                              "in your own submission-script loop) -> cycle_00/image_NN/ + a "
                              "printed loop snippet.\n"
                              "  4 = 100%% real-DFT NEB from a plain interpolated path, no MACE "
                              "at all -> cycle_00/image_NN/ + a printed loop snippet (same "
                              "stb-nebCycle tool as mode 3, more cycles expected).")
    parser.add_argument("--ml-model", choices=["small", "medium", "large"], default="small")
    parser.add_argument("--ml-device", choices=["cpu", "cuda"], default="cpu")
    parser.add_argument("--ml-fmax", type=float, default=0.05)
    parser.add_argument("--ml-k", type=float, default=0.1,
                         help="NEB spring constant, eV/Ang^2 (default: 0.1, ASE's own default).")
    parser.add_argument("--ml-max-steps", type=int, default=200)
    parser.add_argument("--ml-freeze-substrate", dest="ml_freeze_substrate", action="store_true",
                         default=True,
                         help="With --mode 1/2/3: freeze atoms whose position barely differs "
                              "between --initial and --final (default: ON) -- fewer degrees of "
                              "freedom for the MACE optimizer and avoids spurious drift of atoms "
                              "that are spectators to the reaction. See --ml-freeze-threshold.")
    parser.add_argument("--no-ml-freeze-substrate", dest="ml_freeze_substrate", action="store_false",
                         help="Let every atom relax during the MACE path-shaping stage (modes "
                              "1/2/3), even ones that don't move between --initial and --final.")
    parser.add_argument("--ml-freeze-threshold", type=float, default=0.3,
                         help="Displacement threshold in Ang (default: 0.3) below which an atom "
                              "is considered a spectator and frozen by --ml-freeze-substrate.")
    parser.add_argument("--ml-prerelax-endpoints", action="store_true",
                         help="Relax both endpoints' positions (cell fixed) with MACE-MP-0 "
                              "before interpolating -- independent of --mode (a cheap safety "
                              "net when you're not fully certain the endpoints are already "
                              "relaxed). Needs the optional 'ml' extra.")

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
        print(color_text(f"[ERROR] Initial structure/folder '{args.initial}' not found.", 'red'))
        sys.exit(1)
    if not os.path.exists(args.final):
        print(color_text(f"[ERROR] Final structure/folder '{args.final}' not found.", 'red'))
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

    # --initial/--final each accept a bare .fdf file (used exactly as given)
    # or a directory -- e.g. one of stb-nebSites' own site_A/site_B folders
    # -- in which case the folder's finished 'siesta.XV' is preferred over
    # its pre-relaxation 'structure.fdf' guess (structure_io.py::
    # read_relaxed_or_input, shared with stb-adsorbBsse's own "prefer the
    # relaxed geometry" need). used_relaxed is False only when a directory
    # was given with no .XV in it yet -- flagged in [0] below, not fatal
    # (useful to preview the path before the real SIESTA relaxation
    # finishes), unlike stb-adsorbBsse, which skips a not-yet-relaxed site
    # outright: a NEB image folder from an unrelaxed guess is still a
    # legitimate (if weaker) starting point, worth generating regardless.
    try:
        initial_structure, initial_relaxed = structure_io.read_relaxed_or_input(args.initial)
        final_structure, final_relaxed = structure_io.read_relaxed_or_input(args.final)
    except (FileNotFoundError, ValueError) as e:
        print(color_text(f"[ERROR] {e}", 'red'))
        sys.exit(1)
    check_composition_match(initial_structure, final_structure)
    manifest_pair = resolve_manifest_pair(args.initial, args.final, initial_structure, final_structure)
    if manifest_pair is not None:
        args.autosort_tol = 0.0
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
        print_dual(f"  Initial structure : {args.initial}"
                    + ("" if initial_relaxed else " (NOT YET RELAXED)"), f_out)
        print_dual(f"  Final structure   : {args.final}"
                    + ("" if final_relaxed else " (NOT YET RELAXED)"), f_out)
        if not initial_relaxed or not final_relaxed:
            unrelaxed = ", ".join(
                p for p, relaxed in ((args.initial, initial_relaxed), (args.final, final_relaxed))
                if not relaxed)
            print_dual(color_text(
                f"  [WARNING] {unrelaxed}: no finished 'siesta.XV' found in this directory yet -- "
                "using the pre-relaxation structure.fdf guess instead. stb-neb officially expects "
                "already-relaxed endpoints (Section 1 of examples/4.8-adsorption/README.md's "
                "reasoning applies here too); the band below is still generated, useful to preview "
                "the path, but re-run once SIESTA has actually relaxed this endpoint.", 'yellow'), f_out)
        print_dual(f"  Atom correspondence: "
                    + (f"PROVEN via {neb_manifest.MANIFEST_FILENAME} (stb-nebSites site_A/site_B pair)"
                       if manifest_pair is not None
                       else f"distance-based matching (--autosort-tol {args.autosort_tol})"), f_out)
        print_dual(f"  Calc template     : {args.calc}", f_out)
        print_dual(f"  Pseudo dir        : {args.pseudo_dir or '(none)'}", f_out)
        print_dual(f"  Output dir        : {output_root}", f_out)
        print_dual(f"  N images          : {args.n_images}", f_out)
        print_dual(f"  IDPP refinement   : {'yes' if args.idpp else 'no'}", f_out)
        print_dual(f"  Mode              : {args.mode} ({_MODE_DESCRIPTIONS[args.mode]})", f_out)
        mace_used_this_mode = args.mode in (1, 2, 3)
        if mace_used_this_mode:
            print_dual(f"  MACE-MP-0 path shaping: yes (model={args.ml_model}, k={args.ml_k}, "
                        f"fmax={args.ml_fmax})", f_out)
            print_dual(f"  ML-NEB freeze substrate: {'yes' if args.ml_freeze_substrate else 'no'}"
                        + (f" (threshold {args.ml_freeze_threshold} Ang)"
                           if args.ml_freeze_substrate else ""), f_out)
        else:
            print_dual("  MACE-MP-0 path shaping: no (mode 4 -- real DFT from the start)", f_out)
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

        check_endpoint_displacement(initial_pmg, final_pmg_matched, f_out,
                                     threshold=args.ml_freeze_threshold, idpp_used=args.idpp)

        if manifest_pair is not None:
            _manifest_i, _manifest_f, pair_ids_match = manifest_pair
            print_dual(color_text(
                "  [OK] Atom order PROVEN via neb_manifest.json -- skipping distance-based "
                "matching entirely (forcing --autosort-tol 0).", 'green'), f_out)
            if not pair_ids_match:
                print_dual(color_text(
                    "  [WARNING] --initial/--final's manifests have different pair_id -- they "
                    "don't look like they came from the same stb-nebSites invocation (their "
                    "species_sequence still matches exactly, so interpolation proceeds normally; "
                    "double-check these really are the intended site_A/site_B pair).", 'yellow'), f_out)

        try:
            pmg_images = linear_interpolate_images(initial_pmg, final_pmg_matched, args.n_images,
                                                    autosort_tol=args.autosort_tol)
        except ValueError as e:
            if args.autosort_tol != 0:
                # Distance-based matching failed at the requested tolerance. Rather than making
                # the user re-run with --autosort-tol 0 by hand (the fix is almost always this,
                # per the guidance below -- and this exact failure was reported live, twice, by
                # a user who kept hitting it through the interactive menu without setting the
                # advanced option), fall back automatically to index-based matching (--autosort-
                # tol 0 semantics: use --initial/--final's own atom order directly, no distance
                # matching at all). This is safe to do silently-but-reported because a genuinely
                # wrong fallback correspondence would show up as an implausibly large jump in
                # check_path_quality() below, which already warns on exactly that.
                print_dual(color_text(
                    f"[WARNING] Could not match atoms at --autosort-tol {args.autosort_tol} Ang "
                    f"({e}). Falling back to using --initial/--final's own atom order directly "
                    "(equivalent to --autosort-tol 0) -- the right correspondence whenever both "
                    "endpoints share a guaranteed matching atom order (e.g. a stb-nebSites "
                    "site_A/site_B pair, or any pair built from the same base structure). If "
                    "--initial/--final were instead built independently and may have a genuinely "
                    "different atom order, check the path-quality section below for an implausibly "
                    "large per-image displacement -- that would be the sign this fallback guessed "
                    "wrong, not this warning alone.", 'yellow'), f_out)
                try:
                    pmg_images = linear_interpolate_images(initial_pmg, final_pmg_matched,
                                                            args.n_images, autosort_tol=0)
                except ValueError:
                    pmg_images = None
            else:
                pmg_images = None

            if pmg_images is None:
                print_dual(color_text(
                    f"[ERROR] Could not match atoms between --initial and --final for interpolation: "
                    f"{e}. This usually means --autosort-tol ({args.autosort_tol} Ang) is too tight "
                    "for how much some atoms actually moved between the two (already-relaxed) "
                    "endpoints -- e.g. real substrate reconstruction near an adsorption/defect site "
                    "can easily exceed pymatgen's own suggested 0.5 Ang default.", 'red'), f_out)
                print_dual(color_text(
                    "  If --initial/--final came from stb-nebSites (site_A/site_B): they already "
                    "share the EXACT same atom order (same base slab, same construction) -- pass "
                    "--autosort-tol 0 to use that order directly instead of re-matching by distance. "
                    "This is usually the right fix (and more robust than just raising the tolerance: "
                    "verified live that even a much larger tolerance can still fail on a small/"
                    "densely-packed periodic cell, where multiple atoms of the same species sit "
                    "within any reasonably large distance of each other, making the distance-based "
                    "match genuinely ambiguous no matter how loose the tolerance is).", 'yellow'), f_out)
                print_dual(color_text(
                    "  If the two endpoints were built independently (not from the same "
                    "stb-nebSites pair) and may have a different atom order: try a moderately "
                    "larger --autosort-tol first (e.g. 1.0 Ang); if that still fails, --initial/"
                    "--final may not actually be the same physical system (same atom count/"
                    "species), or the correspondence is genuinely ambiguous and needs to be fixed "
                    "by hand.", 'yellow'), f_out)
                if f_out:
                    f_out.close()
                sys.exit(1)
        print_dual(f"  Linear interpolation: {len(pmg_images)} images.", f_out)

        if args.idpp:
            pmg_images = idpp_refine_images(pmg_images, fmax=args.idpp_fmax, steps=args.idpp_steps)
            print_dual("  IDPP refinement: interior images refined with ASE's Image Dependent "
                        "Pair Potential method.", f_out)

        ml_neb_used = False
        ase_images, mace_converged, mace_barrier, mace_dE = None, None, None, None
        if mace_used_this_mode:
            require_mace()
            from stb.core import mace_relax
            print_dual(f"  {color_text('MACE-MP-0:', 'cyan')} running a real climbing-image NEB "
                        "to shape the path ...", f_out)
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
            mace_converged, s1, s2, energies = mace_relax.relax_neb(
                ase_images, calc_mace_neb, k=args.ml_k, fmax=args.ml_fmax,
                max_steps=args.ml_max_steps)
            pmg_images = [AseAtomsAdaptor.get_structure(a) for a in ase_images]
            ml_neb_used = True
            print_dual(f"  {'Converged' if mace_converged else 'Hit step cap, not fully converged'} "
                        f"after {s1} (stage 1) + {s2} (stage 2, climbing) step(s).", f_out)

            from ase.mep.neb import NEBTools
            mace_barrier, mace_dE = NEBTools(ase_images).get_barrier(fit=True)
            print_dual(f"  MACE barrier estimate: {mace_barrier:.4f} eV (forward), reaction "
                        f"energy: {mace_dE:.4f} eV (fitted spline over the relaxed band).", f_out)
            print_dual(color_text(
                "  Note: an ML-level estimate from a fast surrogate potential, not a DFT "
                "adsorption/reaction barrier.", 'yellow'), f_out)

            preview_path = os.path.join(output_root, "neb_ml_preview.png")
            write_ml_preview_plot(ase_images, preview_path)
            print_dual(f"  {color_text('[Saved]', 'cyan')} {preview_path}", f_out)

        reaction_coords = cumulative_reaction_coordinates(pmg_images)
        check_path_quality(pmg_images, reaction_coords, f_out, manifest_proven=(manifest_pair is not None))

        trajectory_path = os.path.join(output_root, "neb_path.xyz")
        write_path_trajectory(pmg_images, trajectory_path)
        print_dual(f"  {color_text('[Saved]', 'cyan')} {trajectory_path} (all images, viewable in "
                    "VESTA/OVITO/ASE-GUI)", f_out)

        if args.mode == 1:
            _write_mode1_json(output_root, ase_images, reaction_coords, mace_converged,
                               mace_barrier, mace_dE, args, f_out)
        else:
            images_root = output_root if args.mode == 2 else os.path.join(output_root, "cycle_00")
            print_dual(f"\n{color_text('[2] IMAGE FOLDERS', 'bold')}", f_out)
            report_rows = []  # (label, index, reaction_coord, dir)
            for i, pmg_image in enumerate(pmg_images):
                label = f"image_{i:02d}"
                image_dir = os.path.join(images_root, label)
                write_image_folder(image_dir, pmg_image, single_point_calc_text, species_meta,
                                    args.pseudo_dir)
                print_dual(f"  {color_text('[OK]', 'green')} {image_dir}", f_out)
                report_rows.append((label, i, reaction_coords[i], image_dir))

            print_dual(f"\n{color_text('[3] SUMMARY', 'bold')}", f_out)
            print_dual(f"  {len(pmg_images)} image folder(s) written under '{images_root}'.", f_out)
            if args.mode == 2:
                print_dual("  Run SIESTA (single-point: MD.NumCGsteps forced to 0) in every "
                            "image_NN/ folder, then use stb-nebAnalysis.", f_out)
            else:
                print_dual("  Run SIESTA (single-point) in every cycle_00/image_NN/ folder, then "
                            "start the refinement loop below -- see [4].", f_out)
                _print_cluster_submission_snippet(f_out, args.mode)

            f_out.write(f"\n# MODE: {args.mode}\n")
            f_out.write(f"# ML_NEB_USED: {'yes' if ml_neb_used else 'no'}\n")
            f_out.write("# IMAGE_TABLE -- parsed by stb-nebAnalysis, do not reorder the columns\n")
            f_out.write(f"# {'label':<14}{'index':<8}{'reaction_coord':<18}{'dir'}\n")
            for label, index, reaction_coord, image_dir in report_rows:
                f_out.write(f"{label:<16}{index:<8}{reaction_coord:<18.6f}{image_dir}\n")

    if args.mode == 1:
        print(f"\n{color_text('Success:', 'green')} MACE-MP-0 NEB complete -- "
              f"{os.path.join(output_root, MACE_RESULT_FILE)} written, no SIESTA folders "
              "(mode 1).")
    else:
        print(f"\n{color_text('Success:', 'green')} {len(pmg_images)} image folder(s) written "
              f"under '{images_root}'.")
    print(f"Full report: {report_path}")


if __name__ == "__main__":
    main()
