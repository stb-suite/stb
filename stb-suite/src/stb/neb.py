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
import time
import argparse
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from stb.core import structure_io
from stb.core import neb_manifest
from stb.core import adsorption_sites
from stb.core.cli import (color_text, show_intro, print_dual, print_section,
                           print_progress_line, finish_progress_line, capture_library_noise)
from stb.core.pseudopotentials import resolve_pseudo_source, copy_pseudo
from stb.core.deps import require_mace

REPORT_FILE = "neb_setup.txt"
MACE_RESULT_FILE = "neb_mace_result.json"
RUN_SUBDIR = "neb_run"

MODE_DESCRIPTIONS = {
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
    neb_manifest.json (written by any tool that proves atom-index
    correspondence between two structures -- see core/neb_manifest.py),
    loads them, cross-checks them against EACH OTHER
    (neb_manifest.validate_manifest_pair) AND against the structures
    actually read back by read_relaxed_or_input -- the second check catches
    a stale manifest (structure.fdf/siesta.XV regenerated or hand-edited
    since the manifest was written), which a manifest-vs-manifest-only
    comparison can't see.

    Returns (manifest_initial, manifest_final, pair_ids_match) when both
    sides have a manifest and everything checks out. Returns None if
    EITHER side has no manifest at all -- not an error, just "no proof
    available" (a hand-built --initial/--final, the common case) -- falls
    back to the existing distance-based --autosort-tol behavior. Exits
    (clean [ERROR], same convention as check_composition_match) for every
    case where a manifest genuinely IS present but can't be trusted.
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
                "likely regenerated or hand-edited after the manifest was written. Delete the "
                "stale manifest before retrying.", 'red'))
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


def require_lattice_match(initial_pmg, final_pmg, f_out, tol=1e-3):
    """Hard-validates that the two endpoints' lattices agree within `tol`
    Ang in every matrix component -- not just stylistic:
    ase.mep.neb.NEB/idpp_interpolate both raise NotImplementedError on any
    per-image cell mismatch (no variable-cell NEB support in ASE), and
    pymatgen's Structure.interpolate raises ValueError on unequal lattices
    unless interpolate_lattices=True -- so a single fixed lattice for the
    whole band is the only thing either downstream library can actually
    run, confirmed live against both APIs.

    Above `tol`, this is now a hard [ERROR] + exit(1) (previously only a
    [WARNING] that silently rebuilt --final onto --initial's lattice): a
    mismatch this large means the user's own two "already-relaxed"
    endpoints don't actually share a cell, which stb-neb should surface
    rather than paper over -- re-relaxing --final on --initial's own cell
    is the user's job, not this tool's.

    Below `tol` (ordinary floating-point/rounding noise between two
    independently-relaxed endpoint calculations, not a real difference),
    silently returns a copy of `final_pmg` rebuilt bit-for-bit onto
    initial_pmg's lattice (species/frac_coords kept from final_pmg) --
    still needed even for a tiny numerical difference, since pymatgen/
    ASE's own lattice-equality checks are tighter than this tool's `tol`.
    """
    initial_matrix = np.array(initial_pmg.lattice.matrix)
    final_matrix = np.array(final_pmg.lattice.matrix)
    max_diff = float(np.abs(initial_matrix - final_matrix).max())
    if max_diff > tol:
        print_dual(color_text(
            f"[ERROR] Initial and final structures have different lattices (largest component "
            f"difference: {max_diff:.4f} Ang, tolerance: {tol} Ang) -- ase.mep.neb.NEB and "
            "pymatgen's interpolation both require every image to share one exact cell (no "
            "variable-cell NEB support in ASE). stb-neb requires --initial/--final to already "
            "share the same lattice (it no longer silently overrides one with the other) -- "
            "re-relax --final with --initial's own cell (e.g. a fixed-cell relaxation, or copy "
            "--initial's %block LatticeVectors into --final's structure) before retrying.", 'red'),
            f_out)
        f_out.close()
        sys.exit(1)
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


def write_image_folder(out_dir, pmg_structure, calc_text, species_meta, pp_path,
                        force_spin=False, force_vdw=False, force_dipole=False):
    """Writes structure.fdf + config_extra.fdf + calc.fdf + copied
    pseudopotentials for one image_NN/ folder. config_extra.fdf is always
    written and always %include'd at the very top of calc.fdf
    (structure_io.prepend_include) -- same "config_extra.fdf as the
    standard mechanism for forcing calculation directives" convention as
    core/adsorption_sites.py's write_reference_folder, but WITHOUT that
    function's FIXED_CELL_BLOCK (meaningless here -- an image's cell is
    fixed by construction, shared with every other image on the band, not
    something that needs a directive to enforce) and without calling
    write_reference_folder itself (that function is slab/pymatgen
    -adsorption-specific). Unlike force_spin/force_vdw/force_dipole,
    SINGLE_POINT_BLOCK is unconditional (no opt-out): a NEB image is one
    independent, uncoupled sample point on the band, so letting SIESTA
    relax it on its own would silently move it off the interpolated path
    -- caught the calc.fdf template's OWN MD.TypeOfRun/MD.Steps
    always losing to config_extra.fdf's (first-occurrence-wins), so
    `calc_text` here is the caller's ORIGINAL --calc text, unmodified.
    """
    os.makedirs(out_dir, exist_ok=True)
    fdf_structure = structure_io.from_pymatgen(pmg_structure, species_meta=species_meta,
                                                coord_format="fractional")
    structure_io.write_fdf(fdf_structure, os.path.join(out_dir, "structure.fdf"))
    with open(os.path.join(out_dir, adsorption_sites.CONFIG_EXTRA_FILE), "w") as f:
        f.write(adsorption_sites.SINGLE_POINT_BLOCK)
        if force_spin:
            f.write(adsorption_sites.SPIN_POLARIZED_BLOCK)
        if force_dipole:
            f.write(adsorption_sites.DIPOLE_CORRECTION_BLOCK)
        if force_vdw:
            f.write(adsorption_sites.VDW_CORRECTION_BLOCK)
    with open(os.path.join(out_dir, "calc.fdf"), "w") as f:
        f.write(structure_io.prepend_include(calc_text, adsorption_sites.CONFIG_EXTRA_FILE))
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


def _print_library_warnings_section(f_out, library_warnings):
    """Prints the shared '[9] LIBRARY WARNINGS' section. Must be called
    BEFORE any '# MODE:'/'# IMAGE_TABLE' marker lines are written to f_out
    -- stb-nebAnalysis's read_image_table() treats every non-blank,
    non-'#'-prefixed line after '# IMAGE_TABLE' as a table row, with no
    explicit end-of-table marker, so the marker block must stay the very
    last thing written to the report file.
    """
    print_section('[9] LIBRARY WARNINGS', f_out)
    if library_warnings:
        print_dual(color_text(
            "  Messages emitted by external libraries (MACE/torch/ASE/pymatgen) during this "
            "run -- collected here instead of interleaved with the report above; harmless in "
            "almost every case (import-time notices, deprecation-style warnings), but worth a "
            "glance.", 'cyan'), f_out)
        for entry in library_warnings:
            print_dual(f"  {entry}", f_out)
    else:
        print_dual("  No library warnings.", f_out)


def _write_mode1_json(run_root, ase_images, reaction_coords, converged, barrier, dE, args, f_out,
                       run_start, library_warnings):
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
    json_path = os.path.join(run_root, MACE_RESULT_FILE)
    with open(json_path, "w") as f:
        json.dump(payload, f, indent=2)

    print_section('[6] MACE RESULT (JSON)', f_out)
    print_dual(f"  {color_text('[Saved]', 'cyan')} {json_path} ({len(ase_images)} image(s), "
                f"barrier {barrier:.4f} eV forward / {backward:.4f} eV backward, reaction "
                f"energy {dE:.4f} eV)", f_out)

    elapsed = time.monotonic() - run_start
    print_section('[8] SUMMARY', f_out)
    print_dual("  Mode 1: no SIESTA folders written -- the path was fully converged in "
                f"MACE-MP-0. Run stb-nebAnalysis --dir {run_root} to read the JSON directly.",
                f_out)
    print_dual(f"  Total elapsed     : {elapsed:.1f}s", f_out)

    _print_library_warnings_section(f_out, library_warnings)

    f_out.write("\n# MODE: 1\n")
    f_out.write(f"# MACE_RESULT_FILE: {MACE_RESULT_FILE}\n")


_CLIMB_AFTER_SUGGESTION = {3: 0, 4: 5}


def _print_cluster_submission_snippet(f_out, mode, run_root, max_cycles, siesta_exe, mpirun_np,
                                       conda_env, cycle_fmax):
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

    max_cycles/siesta_exe/mpirun_np/conda_env/cycle_fmax only shape this
    PRINTED snippet -- stb-neb itself never executes SIESTA in modes 3/4,
    so none of these affect this run's own behavior. cycle_fmax is
    stb-nebCycle's OWN --fmax (real-DFT NEB force-convergence threshold,
    eV/Ang) -- a completely different value from this tool's own
    --ml-fmax (the MACE-MP-0 path-shaping/pre-relax target), just
    printed here since stb-nebCycle is CLI-only and never runs as part
    of this call. Cycle numbers are built via
    `printf "%02d"`, not `seq -w`: `seq -w`'s zero-padding widens to match
    the LARGEST number in the range (3 digits once max_cycles >= 100,
    e.g. "005"), which would then fail to match the 2-digit-MINIMUM
    `cycle_NN` folders stb-nebCycle actually writes via Python's own
    f"cycle_{{n:02d}}" (minimum-width, so "cycle_05" stays 2 digits
    regardless of how many cycles exist overall) -- printf "%02d" mirrors
    that exact minimum-width semantics at any --max-cycles value.
    """
    climb_after = _CLIMB_AFTER_SUGGESTION[mode]
    print_section('[7] CLUSTER SUBMISSION', f_out)
    print_dual(f"  stb-nebCycle is a separate CLI-only tool (not in the interactive stb-suite "
                f"menu) that reads the latest finished cycle_NN/ and takes one real-DFT NEB "
                f"step -- run it in a loop from your own submission script, alternating with "
                f"real SIESTA runs (every cycle_NN/ lives under '{run_root}'):", f_out)

    if conda_env:
        conda_block = f'source ~/miniconda3/etc/profile.d/conda.sh\nconda activate {conda_env}\n'
    else:
        conda_block = ('# TODO: activate your own environment (e.g. `conda activate <env>`) '
                        'if needed\n')
    mpirun_prefix = f"mpirun -np {mpirun_np} " if mpirun_np else ""

    snippet = f"""
{conda_block}
cd "{run_root}"
for cyc in $(seq 0 {max_cycles - 1}); do
    cycle=$(printf "%02d" "$cyc")
    cdir="cycle_${{cycle}}"
    [ -d "$cdir" ] || break
    for img in "$cdir"/image_*; do
        (cd "$img" && {mpirun_prefix}{siesta_exe} calc.fdf --out calc.out)
    done
    stb-nebCycle --dir . --fmax {cycle_fmax} --climb-after {climb_after} --no-intro
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
                              "directory containing one -- a directory's own finished siesta.XV, "
                              "if present, is preferred over its structure.fdf guess.")
    parser.add_argument("-f", "--final", type=str, required=True,
                         help="Final (already-relaxed) endpoint -- same file-or-directory "
                              "convention as --initial, same composition.")
    parser.add_argument("-c", "--calc", type=str, required=True,
                         help="calc.fdf template (kgrid, basis, XC, %%include structure.fdf, "
                              "etc.) -- copied into every image_NN/ folder with MD.TypeOfRun/"
                              "MD.Steps forced to a single-point evaluation.")
    parser.add_argument("-p", "--pseudo-dir", type=str, default="",
                         help="Pseudopotentials source (optional): a bundled bank or a folder path.")

    parser.add_argument("--force-spin", dest="force_spin", action="store_true", default=True,
                         help="Force Spin polarized (via config_extra.fdf, overriding --calc's "
                              "own Spin tag if any) in every image_NN/ (default: ON) -- a "
                              "reacting atom commonly leaves the system with a net magnetic "
                              "moment, which a spin-restricted calculation cannot represent at "
                              "all; costs nothing for a genuinely closed-shell path (converges "
                              "to zero moment). No effect with --mode 1 (no SIESTA folders "
                              "written).")
    parser.add_argument("--no-force-spin", dest="force_spin", action="store_false",
                         help="Leave --calc's own Spin tag untouched in every image_NN/.")
    parser.add_argument("--force-vdw", dest="force_vdw", action="store_true", default=True,
                         help="Force the Grimme DFT-D3 dispersion (van der Waals) correction "
                              "(default: ON) in every image_NN/'s config_extra.fdf. No effect "
                              "with --mode 1.")
    parser.add_argument("--no-force-vdw", dest="force_vdw", action="store_false",
                         help="Leave --calc's own dispersion settings untouched.")
    parser.add_argument("--force-dipole", dest="force_dipole", action="store_true", default=True,
                         help="Force the slab dipole correction (default: ON) in every "
                              "image_NN/'s config_extra.fdf -- only meaningful for a genuinely "
                              "one-sided slab/surface reaction path; harmless (evaluates near "
                              "zero) otherwise. No effect with --mode 1.")
    parser.add_argument("--no-force-dipole", dest="force_dipole", action="store_false",
                         help="Leave --calc's own dipole settings untouched.")

    parser.add_argument("-n", "--n-images", type=int, default=7,
                         help="Total images along the band, endpoints included (default: 7).")
    parser.add_argument("--autosort-tol", type=float, default=0.5,
                         help="Atom-correspondence tolerance (Ang) for linear interpolation "
                              "(default: 0.5, pymatgen's own suggested value for two endpoints "
                              "that may have a different atom order). Pass 0 to skip distance "
                              "-based matching entirely and use --initial/--final's atom order "
                              "directly -- the right choice (and more robust for a small/"
                              "densely-packed cell) when both endpoints already share a "
                              "guaranteed matching order (e.g. a neb_manifest.json-proven pair). "
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

    parser.add_argument("--max-cycles", type=int, default=30,
                         help="Cluster-submission snippet only (--mode 3/4): max number of "
                              "stb-nebCycle refinement cycles the printed loop will attempt "
                              "before giving up (default: 30). Has no effect on this run itself.")
    parser.add_argument("--siesta-exe", type=str, default="siesta",
                         help="Cluster-submission snippet only (--mode 3/4): SIESTA executable "
                              "name or path for the printed loop to call (default: 'siesta', "
                              "assumed already on PATH).")
    parser.add_argument("--mpirun-np", type=int, default=None, metavar="N",
                         help="Cluster-submission snippet only (--mode 3/4): prefix the SIESTA "
                              "call in the printed loop with 'mpirun -np N'. Omit for no mpirun "
                              "(default).")
    parser.add_argument("--conda-env", type=str, default="",
                         help="Cluster-submission snippet only (--mode 3/4): conda environment "
                              "name to activate at the top of the printed loop. Omit to print a "
                              "generic placeholder reminder instead (default).")
    parser.add_argument("--cycle-fmax", type=float, default=0.05,
                         help="Cluster-submission snippet only (--mode 3/4): the looped "
                              "stb-nebCycle call's own --fmax (real-DFT NEB force-convergence "
                              "threshold, eV/Ang, default: 0.05) -- distinct from this tool's own "
                              "--ml-fmax (the MACE-MP-0 path-shaping target).")

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
    # or a directory -- in which case the folder's finished 'siesta.XV' is preferred over
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

    output_root = args.output_dir
    run_root = os.path.join(output_root, RUN_SUBDIR)
    os.makedirs(run_root, exist_ok=True)
    report_path = os.path.join(run_root, REPORT_FILE)

    run_start = time.monotonic()
    library_warnings = []

    composition_counts = Counter(sym for sym, _ in initial_structure.atoms)
    composition_str = " ".join(f"{sym}{n}" for sym, n in sorted(composition_counts.items()))
    n_atoms_total = sum(composition_counts.values())

    with open(report_path, 'w') as f_out:
        print_section('[0] RUN METADATA', f_out)
        print_dual(f"  Started           : {time.strftime('%Y-%m-%d %H:%M:%S')}", f_out)
        print_dual(f"  Composition       : {composition_str} ({n_atoms_total} atoms)", f_out)
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
                    + (f"PROVEN via {neb_manifest.MANIFEST_FILENAME}"
                       if manifest_pair is not None
                       else f"distance-based matching (--autosort-tol {args.autosort_tol})"), f_out)
        print_dual(f"  Calc template     : {args.calc}", f_out)
        print_dual(f"  Pseudo dir        : {args.pseudo_dir or '(none)'}", f_out)
        print_dual(f"  Output dir        : {output_root}", f_out)
        print_dual(f"  Run folder        : {run_root} (every generated file/folder lives here)", f_out)
        print_dual(f"  N images          : {args.n_images}", f_out)
        print_dual(f"  IDPP refinement   : {'yes' if args.idpp else 'no'}", f_out)
        print_dual(f"  Mode              : {args.mode} ({MODE_DESCRIPTIONS[args.mode]})", f_out)
        print_dual(f"  Force spin/vdW/dipole: {'ON' if args.force_spin else 'off'}/"
                    f"{'ON' if args.force_vdw else 'off'}/{'ON' if args.force_dipole else 'off'} "
                    "(config_extra.fdf per image_NN/)", f_out)
        if args.mode != 1:
            print_dual("  Single-point (MD.TypeOfRun CG / MD.Steps 0): ON, unconditional "
                        "(config_extra.fdf per image_NN/ -- each image is one independent, "
                        "uncoupled point on the band; SIESTA must not relax it off the path)",
                        f_out)
        if args.mode == 1 and not (args.force_spin and args.force_vdw and args.force_dipole):
            print_dual(color_text(
                "  [NOTE] --no-force-spin/--no-force-vdw/--no-force-dipole have no effect with "
                "--mode 1 -- no SIESTA folders (and therefore no config_extra.fdf) are written "
                "in this mode.", 'yellow'), f_out)
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

        if args.ml_prerelax_endpoints:
            print_section('[1] ML PRE-RELAX ENDPOINTS', f_out)
            prerelax_max_steps = 200
            with capture_library_noise(library_warnings, "MACE import"):
                require_mace()
            from stb.core import mace_relax
            with capture_library_noise(library_warnings, "MACE calculator setup"):
                calc_mace_endpoints = mace_relax.get_calculator(model=args.ml_model, device=args.ml_device)
            for label, pmg in (("initial", initial_pmg), ("final", final_pmg)):
                print_dual(f"  Relaxing the {label} endpoint (positions only, MACE-MP-0) ...", f_out)
                ase_atoms = AseAtomsAdaptor.get_atoms(pmg)

                def _on_step(step, energy, max_force, _label=label):
                    print_progress_line(
                        f"    {_label}: step {step}/{prerelax_max_steps}, E = {energy:.4f} eV, "
                        f"max|F| = {max_force:.4f} eV/Ang (target {args.ml_fmax})",
                        step, prerelax_max_steps)

                with capture_library_noise(library_warnings, f"MACE pre-relax ({label})"):
                    converged, steps = mace_relax.relax(ase_atoms, calc_mace_endpoints,
                                                         fmax=args.ml_fmax, max_steps=prerelax_max_steps,
                                                         on_step=_on_step)
                finish_progress_line()
                relaxed_pmg = AseAtomsAdaptor.get_structure(ase_atoms)
                if label == "initial":
                    initial_pmg = relaxed_pmg
                else:
                    final_pmg = relaxed_pmg
                print_dual(f"  {'Converged' if converged else 'Hit step cap, not fully converged'} "
                            f"({label}) after {steps} step(s).", f_out)

        print_section('[2] ENDPOINT CHECKS', f_out)
        initial_pmg = wrap_into_cell(initial_pmg)
        final_pmg = wrap_into_cell(final_pmg)
        final_pmg_matched = require_lattice_match(initial_pmg, final_pmg, f_out)
        print_dual(color_text(
            "  [OK] Endpoints share the same composition and lattice (validated above).",
            'green'), f_out)

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
                    "don't look like they came from the same tool invocation (their "
                    "species_sequence still matches exactly, so interpolation proceeds normally; "
                    "double-check these really are the intended endpoint pair).", 'yellow'), f_out)

        print_section('[3] PATH INTERPOLATION', f_out)
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
                    "endpoints share a guaranteed matching atom order (e.g. any pair built from "
                    "the same base structure). If "
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
                    "  If --initial/--final were built from the same base structure (e.g. two "
                    "derived copies of one reference calculation): they already share the EXACT "
                    "same atom order -- pass --autosort-tol 0 to use that order directly instead "
                    "of re-matching by distance. This is usually the right fix (and more robust "
                    "than just raising the tolerance: verified live that even a much larger "
                    "tolerance can still fail on a small/densely-packed periodic cell, where "
                    "multiple atoms of the same species sit within any reasonably large distance "
                    "of each other, making the distance-based match genuinely ambiguous no matter "
                    "how loose the tolerance is).", 'yellow'), f_out)
                print_dual(color_text(
                    "  If the two endpoints were built independently and may have a different "
                    "atom order: try a moderately larger --autosort-tol first (e.g. 1.0 Ang); if "
                    "that still fails, --initial/--final may not actually be the same physical "
                    "system (same atom count/species), or the correspondence is genuinely "
                    "ambiguous and needs to be fixed by hand.", 'yellow'), f_out)
                if f_out:
                    f_out.close()
                sys.exit(1)
        print_dual(f"  [OK] Linear interpolation: {len(pmg_images)} images "
                    f"(image_00 .. image_{len(pmg_images) - 1:02d}).", f_out)

        if args.idpp:
            print_dual(f"  Refining interior images with ASE's IDPP method "
                        f"(fmax={args.idpp_fmax}, max {args.idpp_steps} steps) ...", f_out)
            with capture_library_noise(library_warnings, "IDPP refinement"):
                pmg_images = idpp_refine_images(pmg_images, fmax=args.idpp_fmax, steps=args.idpp_steps)
            print_dual("  [OK] IDPP refinement applied to interior images.", f_out)

        ml_neb_used = False
        ase_images, mace_converged, mace_barrier, mace_dE = None, None, None, None
        if mace_used_this_mode:
            print_section('[4] MACE-MP-0 PATH SHAPING', f_out)
            mace_start = time.monotonic()
            with capture_library_noise(library_warnings, "MACE import"):
                require_mace()
            from stb.core import mace_relax
            print_dual(f"  Model             : MACE-MP-0 ({args.ml_model}, device={args.ml_device})", f_out)
            print_dual(f"  Spring constant k : {args.ml_k} eV/Ang^2", f_out)
            print_dual(f"  Target fmax       : {args.ml_fmax} eV/Ang (max {args.ml_max_steps} steps)", f_out)
            ase_images = [AseAtomsAdaptor.get_atoms(s) for s in pmg_images]
            if args.ml_freeze_substrate:
                frozen_indices = compute_frozen_indices(pmg_images, threshold=args.ml_freeze_threshold)
                if frozen_indices:
                    from ase.constraints import FixAtoms
                    for atoms in ase_images:
                        atoms.set_constraint(FixAtoms(indices=frozen_indices))
                print_dual(f"  Freezing {len(frozen_indices)}/{len(ase_images[0])} atom(s) with "
                            f"< {args.ml_freeze_threshold} Ang displacement between endpoints.", f_out)

            with capture_library_noise(library_warnings, "MACE calculator setup"):
                calc_mace_neb = mace_relax.get_calculator(model=args.ml_model, device=args.ml_device)

            _stage_names = {1: "shaping", 2: "climbing"}

            def _on_step(stage, step, residual):
                stage_steps = max(args.ml_max_steps // 2, 1) if stage == 1 else args.ml_max_steps
                print_progress_line(
                    f"    stage {stage}/2 ({_stage_names[stage]}): step {step}, "
                    f"max residual force {residual:.4f} eV/Ang (target {args.ml_fmax})",
                    step, stage_steps)

            print_dual("  Running the climbing-image NEB (stage 1: shape the band, "
                        "stage 2: climb to the saddle) ...", f_out)
            with capture_library_noise(library_warnings, "MACE NEB relax"):
                mace_converged, s1, s2, energies = mace_relax.relax_neb(
                    ase_images, calc_mace_neb, k=args.ml_k, fmax=args.ml_fmax,
                    max_steps=args.ml_max_steps, on_step=_on_step)
            finish_progress_line()
            pmg_images = [AseAtomsAdaptor.get_structure(a) for a in ase_images]
            ml_neb_used = True
            mace_elapsed = time.monotonic() - mace_start
            print_dual(f"  [{'OK' if mace_converged else 'WARNING'}] "
                        f"{'Converged' if mace_converged else 'Hit step cap, not fully converged'} "
                        f"after {s1} (stage 1) + {s2} (stage 2, climbing) step(s), "
                        f"{mace_elapsed:.1f}s elapsed.", f_out)

            from ase.mep.neb import NEBTools
            mace_barrier, mace_dE = NEBTools(ase_images).get_barrier(fit=True)
            print_dual(f"  MACE barrier estimate: {mace_barrier:.4f} eV (forward), reaction "
                        f"energy: {mace_dE:.4f} eV (fitted spline over the relaxed band).", f_out)
            print_dual(color_text(
                "  Note: an ML-level estimate from a fast surrogate potential, not a DFT "
                "adsorption/reaction barrier.", 'yellow'), f_out)

            preview_path = os.path.join(run_root, "neb_ml_preview.png")
            write_ml_preview_plot(ase_images, preview_path)
            print_dual(f"  {color_text('[Saved]', 'cyan')} {preview_path}", f_out)

        print_section('[5] PATH QUALITY', f_out)
        reaction_coords = cumulative_reaction_coordinates(pmg_images)
        print_dual(f"  Total path length : {reaction_coords[-1]:.4f} Ang "
                    f"(cumulative Cartesian displacement, image_00 to image_{len(pmg_images) - 1:02d})",
                    f_out)
        check_path_quality(pmg_images, reaction_coords, f_out, manifest_proven=(manifest_pair is not None))

        trajectory_path = os.path.join(run_root, "neb_path.xyz")
        write_path_trajectory(pmg_images, trajectory_path)
        print_dual(f"  {color_text('[Saved]', 'cyan')} {trajectory_path} (all images, viewable in "
                    "VESTA/OVITO/ASE-GUI)", f_out)

        if args.mode == 1:
            _write_mode1_json(run_root, ase_images, reaction_coords, mace_converged,
                               mace_barrier, mace_dE, args, f_out, run_start, library_warnings)
        else:
            images_root = run_root if args.mode == 2 else os.path.join(run_root, "cycle_00")
            print_section('[6] IMAGE FOLDERS', f_out)
            report_rows = []  # (label, index, reaction_coord, dir)
            for i, pmg_image in enumerate(pmg_images):
                label = f"image_{i:02d}"
                image_dir = os.path.join(images_root, label)
                write_image_folder(image_dir, pmg_image, calc_text, species_meta,
                                    args.pseudo_dir, force_spin=args.force_spin,
                                    force_vdw=args.force_vdw, force_dipole=args.force_dipole)
                print_dual(f"  {color_text('[OK]', 'green')} {image_dir}", f_out)
                report_rows.append((label, i, reaction_coords[i], image_dir))
            print_dual(f"  {len(pmg_images)} image folder(s) written under '{images_root}'.", f_out)

            if args.mode in (3, 4):
                _print_cluster_submission_snippet(f_out, args.mode, run_root, args.max_cycles,
                                                   args.siesta_exe, args.mpirun_np, args.conda_env,
                                                   args.cycle_fmax)

            elapsed = time.monotonic() - run_start
            print_section('[8] SUMMARY', f_out)
            print_dual(f"  Mode              : {args.mode} ({MODE_DESCRIPTIONS[args.mode]})", f_out)
            print_dual(f"  Images written    : {len(pmg_images)} under '{images_root}'", f_out)
            print_dual(f"  Total elapsed     : {elapsed:.1f}s", f_out)
            if args.mode == 2:
                print_dual("  Next step: run SIESTA (single-point: MD.Steps forced to 0) in "
                            f"every image_NN/ folder, then run stb-nebAnalysis --dir {run_root}.", f_out)
            else:
                print_dual("  Next step: run SIESTA (single-point) in every cycle_00/image_NN/ "
                            "folder, then start the refinement loop -- see [7].", f_out)

            _print_library_warnings_section(f_out, library_warnings)

            f_out.write(f"\n# MODE: {args.mode}\n")
            f_out.write(f"# ML_NEB_USED: {'yes' if ml_neb_used else 'no'}\n")
            f_out.write("# IMAGE_TABLE -- parsed by stb-nebAnalysis, do not reorder the columns\n")
            f_out.write(f"# {'label':<14}{'index':<8}{'reaction_coord':<18}{'dir'}\n")
            for label, index, reaction_coord, image_dir in report_rows:
                f_out.write(f"{label:<16}{index:<8}{reaction_coord:<18.6f}{image_dir}\n")

    if args.mode == 1:
        print(f"\n{color_text('Success:', 'green')} MACE-MP-0 NEB complete -- "
              f"{os.path.join(run_root, MACE_RESULT_FILE)} written, no SIESTA folders "
              "(mode 1).")
    else:
        print(f"\n{color_text('Success:', 'green')} {len(pmg_images)} image folder(s) written "
              f"under '{images_root}'.")
    print(f"Full report: {report_path}")


if __name__ == "__main__":
    main()
