#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""Stage 1 of the NEB / Reaction Path workflow (item 4.9): given a slab and
an adsorbate/molecule, enumerates the symmetrically distinct adsorption
sites (same AdsorbateSiteFinder + --symprec machinery as stb-adsorb) and
writes TWO endpoint folders -- site_A/, site_B/ -- for the two candidates
the user picks. Those are the "already-relaxed endpoints" stb-neb (now
Stage 2, code 4.9.2) expects, once you've run SIESTA in each.

Distinct from stb-adsorb: no isolated-adsorbate reference, no E_ads
ranking, no --all-sites/--ml-rank/--height-sweep/--both-sides -- this tool
only cares about two concrete points on the slab, not ranking many of
them. Distinct from stb-mldiffusion: that one builds a vacancy-hop pair in
BULK via a nearest-neighbor-ATOM shell search and runs the NEB itself with
MACE; this tool builds an ADSORBATE-hop pair on a SLAB via symmetry-
reduced site enumeration and only writes real SIESTA input folders --
running the actual NEB is stb-neb/stb-nebAnalysis's job (Stages 2/3).
"""

VERSION = "1.0.0"

import os
import sys
import uuid
import argparse
import numpy as np
import matplotlib.pyplot as plt
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.adsorption import AdsorbateSiteFinder, plot_slab, get_rot
from stb.core import structure_io
from stb.core import neb_manifest
from stb.core.cli import color_text, show_intro, print_dual, print_section, print_table, capture_library_noise
from stb.core.pseudopotentials import resolve_pseudo_source
from stb.core.ase_view import view_structure_interactive
from stb.core.adsorption_sites import (
    resolve_adsorbate, resolve_slab_orientation, write_reference_folder,
    min_adsorbate_image_distance, cluster_candidate_coords, _MIN_LATERAL_IMAGE_SEPARATION_ANG,
)

PREP_REPORT_FILE = "neb_sites_prep_report.txt"

# Dot color per site type -- shared by both plotting functions below, and
# consistent between them (a given type always looks the same on either
# plot).
_SITE_TYPE_COLORS = {"ontop": "tab:red", "bridge": "tab:blue", "hollow": "tab:green"}

# plot_slab draws its own repeated (default repeat=5) substrate-atom
# circles at zorder = 2*idx/2*idx+1, where idx runs over ALL drawn copies
# -- for an N-atom cell that's up to 25*N atoms, so zorder can reach
# ~50*N. A low, fixed zorder (e.g. 100) is fine for a 2-atom toy cell but
# gets silently drawn UNDER the substrate for anything bigger (reported
# live: dots/numbers invisible on a real, bigger structure) -- both
# plotting functions below use this same safely-large constant instead so
# our own markers/labels always end up on top, regardless of cell size.
_FOREGROUND_ZORDER = 10_000


def write_candidate_sites_plot(pmg_structure, candidates, out_path, show=False):
    """Saves a top-view PNG of the slab with EVERY enumerated candidate
    site marked with a small dot and its 0-based index (as a text label
    next to the dot) -- the exact same numbering [2]'s printed table (and
    --site-a/--site-b) use, so picking an endpoint pair can be done
    visually instead of only by reading coordinates out of the table. Dot
    color encodes the site type (ontop/bridge/hollow, see
    _SITE_TYPE_COLORS), with one legend entry per type actually present
    among `candidates`.

    Deliberately NOT core/adsorption_sites.py::write_site_plot (which
    stb-adsorb still uses for its own generic 'sanity check' overview): that
    one draws pymatgen's OWN find_adsorption_sites() candidates at its own
    default settings, not necessarily matching this run's --height/
    --symprec/--site-type -- for a tool whose whole point is picking two
    specific numbered candidates, showing a plot with a different site set
    (and no numbers at all) would be actively misleading rather than just
    unhelpful. Same plot_slab(..., adsorption_sites=False) + get_rot()
    rotation write_chosen_sites_plot uses, see that function's docstring
    for why the rotation is needed.
    """
    fig, ax = plt.subplots(figsize=(6, 6))
    plot_slab(pmg_structure, ax, adsorption_sites=False)
    symm_op = get_rot(pmg_structure)
    seen_types = []
    for i, st, coord in candidates:
        x, y = symm_op.operate(coord)[:2]
        color = _SITE_TYPE_COLORS.get(st, "tab:gray")
        label = st if st not in seen_types else None
        seen_types.append(st)
        ax.scatter([x], [y], marker='o', s=40, c=color, edgecolors='k', linewidths=0.6,
                   zorder=_FOREGROUND_ZORDER, label=label)
        ax.annotate(str(i), (x, y), textcoords="offset points", xytext=(5, 5), fontsize=8,
                    fontweight='bold', color=color, zorder=_FOREGROUND_ZORDER + 1)
    n_types = len(set(seen_types))
    if n_types:
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.02), ncol=n_types, fontsize=9)
    ax.set_title("Candidate sites (numbered -- see [2]/[3] table)")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)


def write_chosen_sites_plot(pmg_structure, coord_a, coord_b, out_path, show=False):
    """Saves a top-view PNG of the slab with ONLY the two CHOSEN endpoints
    marked and labeled 'A'/'B' -- distinct from write_candidate_sites_plot's
    candidate_sites.png (every enumerated candidate, numbered), this one is
    specific to the actual --site-a/--site-b pair a run picked, the same
    two points site_A/site_B get written from.

    plot_slab (adsorption_sites=False here, since we're drawing our own
    markers, not pymatgen's candidate-site ones) internally reorients the
    slab so the surface normal aligns with Cartesian z before drawing --
    `get_rot(pmg_structure)` is the exact same rotation, applied by hand to
    `coord_a`/`coord_b` (already-Cartesian points from
    AdsorbateSiteFinder.find_adsorption_sites, in the ORIGINAL unrotated
    frame) so they land in the same 2D axes pymatgen's own site markers
    would -- verified against pymatgen's plot_slab source, which applies
    this identical get_rot() rotation to its own adsorption-site markers
    before plotting them.
    """
    fig, ax = plt.subplots(figsize=(6, 6))
    plot_slab(pmg_structure, ax, adsorption_sites=False)
    symm_op = get_rot(pmg_structure)
    (xa, ya), (xb, yb) = symm_op.operate(coord_a)[:2], symm_op.operate(coord_b)[:2]
    ax.scatter([xa], [ya], marker='*', s=400, c='tab:red', edgecolors='k', linewidths=1.0,
               zorder=_FOREGROUND_ZORDER, label='site A')
    ax.scatter([xb], [yb], marker='*', s=400, c='tab:blue', edgecolors='k', linewidths=1.0,
               zorder=_FOREGROUND_ZORDER, label='site B')
    ax.annotate('A', (xa, ya), textcoords="offset points", xytext=(8, 8), fontsize=12,
                fontweight='bold', color='tab:red', zorder=_FOREGROUND_ZORDER + 1)
    ax.annotate('B', (xb, yb), textcoords="offset points", xytext=(8, 8), fontsize=12,
                fontweight='bold', color='tab:blue', zorder=_FOREGROUND_ZORDER + 1)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.02), ncol=2, fontsize=9)
    ax.set_title("Chosen NEB endpoint sites")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)


def find_candidate_sites(finder, site_types, height, symprec, lattice):
    """Returns (candidates, raw_counts, reduced_counts) where `candidates`
    is a single flat list of (index, site_type, coord) -- one combined,
    0-based enumeration across every requested site type, in the order
    site_types itself lists them (ontop, then bridge, then hollow for
    --site-type all) -- this is the numbering --site-a/--site-b index
    into, since a diffusion-path endpoint pair may legitimately span two
    different site types (e.g. ontop -> hollow). `raw_counts`/
    `reduced_counts` (per type) are for the [2] report table, same
    "how many candidates, how many are symmetry-equivalent" probe
    stb-adsorb's own [2] section already does.

    Every candidate (across all site types) is re-wrapped into the same
    periodic image via cluster_candidate_coords, anchored on the first one
    found -- pymatgen's own symm_reduce/put_coord_inside wrap each
    candidate independently with no shared reference point, so two
    symmetrically-equivalent candidates can land in different periodic
    images and appear far apart despite being the same physical site (see
    cluster_candidate_coords's own docstring).
    """
    raw_counts, reduced_counts = {}, {}
    reduced_by_type = {}
    for st in site_types:
        raw = finder.find_adsorption_sites(distance=height, symm_reduce=0, positions=[st])
        reduced = finder.find_adsorption_sites(distance=height, symm_reduce=symprec, positions=[st])
        raw_counts[st], reduced_counts[st] = len(raw[st]), len(reduced[st])
        reduced_by_type[st] = reduced[st]

    anchor = next((c for st in site_types for c in reduced_by_type[st]), None)
    candidates = []
    for st in site_types:
        coords = reduced_by_type[st]
        if anchor is not None:
            coords = cluster_candidate_coords(coords, lattice, reference=anchor)
        for coord in coords:
            candidates.append((len(candidates), st, coord))
    return candidates, raw_counts, reduced_counts


def nearest_candidate_pair(candidates, lattice):
    """Returns (i, j, distance) for the closest pair of candidate sites by
    minimum-image Cartesian distance between their (already-at-height)
    coordinates -- the site-level analog of stb-mldiffusion's atom-level
    find_neighbor_shell, printed as a suggestion only (never auto-selected
    -- the user picks --site-a/--site-b explicitly). Minimum-image (not
    plain norm) as defense in depth on top of find_candidate_sites' own
    clustering: even if some future caller feeds in un-clustered
    candidates, a genuinely close pair through the periodic boundary won't
    be missed just because their representatives happen to sit in
    different periodic images (same frac_diff -= round(frac_diff) idiom as
    cluster_candidate_coords). None if fewer than 2 candidates exist.
    """
    if len(candidates) < 2:
        return None
    coords = np.array([c for _i, _st, c in candidates])
    inv_lattice = np.linalg.inv(lattice)
    frac = coords @ inv_lattice
    diffs = frac[:, None, :] - frac[None, :, :]
    diffs -= np.round(diffs)
    dists = np.linalg.norm(diffs @ lattice, axis=-1)
    np.fill_diagonal(dists, np.inf)
    i, j = np.unravel_index(np.argmin(dists), dists.shape)
    return int(min(i, j)), int(max(i, j)), float(dists[i, j])


def candidate_table_rows(candidates, suggested_pair):
    rows = []
    suggested = {suggested_pair[0], suggested_pair[1]} if suggested_pair else set()
    for i, st, coord in candidates:
        mark = " (suggested pair)" if i in suggested else ""
        rows.append(([str(i), st, f"{coord[0]:.4f}, {coord[1]:.4f}, {coord[2]:.4f}{mark}"], None))
    return rows


def center_slab_in_vacuum(pmg_structure, tol_ang=1.0):
    """Recenters the slab along c (assumed already relabeled so vacuum is on
    c, see resolve_slab_orientation) so it sits in the middle of the vacuum
    gap instead of touching a periodic boundary.

    A freshly-cut slab (pymatgen's own SlabGenerator, stb-slab) conventionally
    starts near frac z=0 with all vacuum stacked above it. During a real
    SIESTA relaxation, an atom that starts at frac z close to 0 can easily
    drift slightly negative and get wrapped to frac z close to 1 (the TOP of
    the cell) instead -- which then reads as an enormous, spurious
    displacement to anything comparing this relaxed geometry against another
    snapshot of the same atom (e.g. stb-neb matching site_A against site_B).
    Shifting the whole slab -- a rigid translation along c, then wrapping
    into the cell -- to be centered removes this risk entirely without
    changing the physical structure (bond lengths/angles are untouched, only
    the coordinate origin along c moves).

    Returns (structure, shifted, shift_ang). `shifted` is False (structure
    returned unchanged, not even copied) when the slab is already within
    `tol_ang` of centered, to avoid gratuitously rewriting coordinates for a
    structure that doesn't need it.
    """
    frac_z = pmg_structure.frac_coords[:, 2]
    slab_center = (frac_z.min() + frac_z.max()) / 2.0
    shift = 0.5 - slab_center
    c_length = np.linalg.norm(pmg_structure.lattice.matrix[2])
    shift_ang = shift * c_length
    if abs(shift_ang) < tol_ang:
        return pmg_structure, False, shift_ang

    pmg_structure = pmg_structure.copy()
    pmg_structure.translate_sites(
        list(range(len(pmg_structure))), [0.0, 0.0, shift], frac_coords=True, to_unit_cell=True)
    return pmg_structure, True, shift_ang


def min_adsorbate_slab_distance(pmg_structure, n_substrate):
    """Minimum periodic distance between any slab atom (index < n_substrate)
    and any adsorbate atom (index >= n_substrate) -- catches a --height
    that places the adsorbate unphysically close to (or inside) the slab.
    Same small generic helper as adsorb.py's own copy (not worth a
    cross-import for a 6-line distance-matrix slice, same "plain generic
    helper" bar stb-adsorb's own build_sweep_values duplicate already sets).
    """
    n_total = len(pmg_structure)
    if n_substrate == 0 or n_substrate >= n_total:
        return None
    dm = pmg_structure.distance_matrix
    return float(dm[:n_substrate, n_substrate:].min())


def write_neb_manifest(site_dir, pair_id, site_label, adsorbate_name, n_substrate, n_adsorbate):
    """Re-reads the just-written structure.fdf (the authoritative on-disk
    atom order AFTER write_fdf's own species-grouping) and records it as
    site_dir/neb_manifest.json, so a later stb-neb run on this exact folder
    can PROVE its atom order instead of guessing via --autosort-tol.
    """
    written = structure_io.read_fdf(os.path.join(site_dir, "structure.fdf"))
    species_sequence = [symbol for symbol, _ in written.atoms]
    return neb_manifest.write_manifest(
        site_dir, pair_id=pair_id, site_label=site_label,
        species_sequence=species_sequence, n_substrate=n_substrate,
        n_adsorbate=n_adsorbate, adsorbate=adsorbate_name)


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Stage 1 of the NEB / Reaction Path workflow: enumerates "
        "symmetrically distinct adsorption sites on a slab and writes two endpoint folders "
        "(site_A/, site_B/) for the pair you pick -- the inputs stb-neb (Stage 2) expects "
        "once you've relaxed them with SIESTA.", 'bold')}""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s -s slab.fdf -c calc.fdf --adsorbate H\n"
               "  %(prog)s -s slab.fdf -c calc.fdf --adsorbate H --site-a 0 --site-b 2\n"
               "  %(prog)s -s slab.fdf -c calc.fdf --adsorbate H2O --site-type ontop\n"
    )

    parser.add_argument("-s", "--structure", type=str, default="structure.fdf",
                         help="Input slab/2D structure.fdf (vacuum along c). Default: structure.fdf.")
    parser.add_argument("-c", "--calc", type=str, default="calc.fdf",
                         help="calc.fdf template already configured for the slab -- copied as-is "
                              "(cell fixed via a config_extra.fdf sidecar) into 'site_A/'/'site_B/'. "
                              "Default: calc.fdf.")
    parser.add_argument("-p", "--pseudo-dir", type=str, default="",
                         help="Pseudopotentials source (optional): a bundled bank or a folder path.")

    parser.add_argument("--adsorbate", type=str, help="Element symbol (single atom) or ASE G2 "
                         "molecule name (e.g. H2O) -- the ONE species hopping between the two "
                         "chosen sites. Comma-separated lists (stb-adsorb's multi-adsorbate "
                         "comparison) aren't meaningful here.")
    parser.add_argument("--list", action="store_true", help="List the available G2 molecule names and exit.")

    parser.add_argument("--site-type", choices=["ontop", "bridge", "hollow", "all"], default="all",
                         help="Which family/families of candidate sites to enumerate (default: all).")
    parser.add_argument("--height", type=float, default=2.0,
                         help="Adsorption distance in Ang above the surface (default: 2.0), shared "
                              "by every candidate site and both endpoints.")
    parser.add_argument("--symprec", type=float, default=0.01,
                         help="Symmetry-reduction tolerance for site-finding (default: 0.01, "
                              "pymatgen's own default).")
    parser.add_argument("--force-spin", dest="force_spin", action="store_true", default=True,
                         help="Force Spin polarized (via config_extra.fdf, overriding --calc's own "
                              "Spin tag if any) in site_A/site_B (default: ON) -- a single adsorbate "
                              "atom (or a molecule containing one) bonded to a slab commonly leaves "
                              "the combined system with a net magnetic moment (most simply, an odd "
                              "total valence-electron count, which a spin-restricted calculation "
                              "cannot represent at all), same reasoning stb-adsorb already applies "
                              "to its isolated-adsorbate reference. Costs nothing for a genuinely "
                              "closed-shell site (converges to zero moment).")
    parser.add_argument("--no-force-spin", dest="force_spin", action="store_false",
                         help="Leave --calc's own Spin tag untouched in site_A/site_B.")
    parser.add_argument("--vacuum-gap", type=float, default=10.0,
                         help="Vacuum-axis detection threshold in Ang (default: 10.0), same "
                              "convention as stb-adsorb/stb-kgrid/stb-slab.")

    parser.add_argument("--site-a", type=int, default=None, metavar="INDEX",
                         help="0-based index (from the [2] candidate table) of the FIRST endpoint. "
                              "Both --site-a and --site-b must be given together; if either is "
                              "missing, the candidate table is printed and the tool exits so you "
                              "can re-run with both.")
    parser.add_argument("--site-b", type=int, default=None, metavar="INDEX",
                         help="0-based index of the SECOND endpoint. See --site-a.")

    parser.add_argument("-O", "--output-dir", type=str, default="nebsites_run",
                         help="Root directory (default: nebsites_run) for 'site_A' and 'site_B'.")
    parser.add_argument("--save-report", action="store_true",
                         help=f"Also persist the full report narrative to {PREP_REPORT_FILE} (under "
                              "--output-dir). Off by default.")
    parser.add_argument("--view-plots", action="store_true",
                         help="Also show the candidate-site layout plot and (once --site-a/"
                              "--site-b are given) the chosen-sites plot on screen (one "
                              "blocking window each), instead of only saving them as PNGs. "
                              "Off by default.")
    parser.add_argument("--view", action="store_true",
                         help="Once --site-a/--site-b are given: open the clean slab, site_A "
                              "and site_B in ASE's interactive 3D viewer as a multi-frame "
                              "browser, after everything else has finished. Needs a display. "
                              "A separate flag from --view-plots, which is matplotlib.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-nebSites {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    if args.list:
        from ase.collections import g2
        print(color_text("Available ASE G2 molecule names:", 'bold'))
        for name in sorted(g2.names):
            print(f"  {name}")
        sys.exit(0)

    if not args.adsorbate:
        parser.error("--adsorbate is required (unless --list).")
    if (args.site_a is None) != (args.site_b is None):
        parser.error("--site-a and --site-b must be given together (or neither, to see the "
                      "candidate table first).")
    if args.site_a is not None and args.site_a == args.site_b:
        parser.error("--site-a and --site-b must be different indices.")

    if not os.path.exists(args.structure):
        print(color_text(f"[ERROR] Structure file '{args.structure}' not found.", 'red'))
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

    structure = structure_io.read_fdf(args.structure)
    try:
        structure, relabel_note = resolve_slab_orientation(structure, args.vacuum_gap)
    except ValueError as e:
        print(color_text(f"[ERROR] stb-nebSites {e}", 'red'))
        sys.exit(1)
    pmg_structure = structure_io.to_pymatgen(structure)
    pmg_structure, slab_recentered, recenter_shift_ang = center_slab_in_vacuum(pmg_structure)
    slab_species_meta = structure_io.species_dict(structure)
    n_substrate = len(pmg_structure)

    try:
        mol = resolve_adsorbate(args.adsorbate)
    except ValueError as e:
        print(color_text(f"[ERROR] {e}", 'red'))
        sys.exit(1)

    with open(args.calc) as f:
        calc_text = f.read()

    finder = AdsorbateSiteFinder(pmg_structure)
    site_types = ["ontop", "bridge", "hollow"] if args.site_type == "all" else [args.site_type]

    output_root = args.output_dir
    os.makedirs(output_root, exist_ok=True)
    report_path = os.path.join(output_root, PREP_REPORT_FILE) if args.save_report else None
    f_out = open(report_path, 'w') if report_path else None
    library_warnings = []
    view_frames = [("clean_slab", AseAtomsAdaptor.get_atoms(pmg_structure))]  # (label, ase.Atoms), for --view

    print_dual(color_text("===== STB-NEBSITES REPORT =====", 'magenta'), f_out)

    # --- [0] RUN METADATA ---
    print_section("[0] RUN METADATA", f_out)
    print_dual(f"  Structure       : {args.structure}", f_out)
    print_dual(f"  Calc template   : {args.calc}", f_out)
    print_dual(f"  Pseudo dir      : {args.pseudo_dir or '(none)'}", f_out)
    print_dual(f"  Output dir      : {output_root}", f_out)
    print_dual(f"  Adsorbate       : {args.adsorbate} ({len(mol)} atom(s))", f_out)
    print_dual(f"  Site type       : {args.site_type}", f_out)
    print_dual(f"  Height          : {args.height:.2f} Ang", f_out)
    print_dual(f"  symprec         : {args.symprec}", f_out)
    print_dual(f"  Force spin      : {'yes -> Spin polarized' if args.force_spin else 'no (--calc as-is)'}",
                f_out)
    print_dual(f"  Save report     : {'yes -> ' + report_path if args.save_report else 'no'}", f_out)

    # --- [1] REFERENCE SLAB ---
    print_section("[1] REFERENCE SLAB", f_out)
    if relabel_note:
        print_dual(color_text(f"  [INFO] {relabel_note}", 'yellow'), f_out)
    else:
        print_dual("  Vacuum axis already on c -- no relabeling needed.", f_out)
    if slab_recentered:
        print_dual(color_text(
            f"  [INFO] Slab was near a cell boundary along c (uneven vacuum split) -- shifted by "
            f"{recenter_shift_ang:+.3f} Ang to center it in the vacuum gap. This avoids a real "
            "risk: an atom starting near frac z=0 can drift slightly negative during SIESTA "
            "relaxation and wrap around to frac z~1 (the top of the cell) instead, which reads "
            "as an enormous, spurious displacement when stb-neb later matches this site against "
            "its pair.", 'yellow'), f_out)
    else:
        print_dual("  Slab already centered in the vacuum gap (no boundary-wrap risk).", f_out)

    lateral_dist = min_adsorbate_image_distance(mol.cart_coords, pmg_structure.lattice.matrix)
    if lateral_dist < _MIN_LATERAL_IMAGE_SEPARATION_ANG:
        print_dual(color_text(
            f"  [WARNING] '{args.adsorbate}' on this slab: only {lateral_dist:.2f} Ang from its "
            f"own periodic image in the ab-plane (recommended >= "
            f"{_MIN_LATERAL_IMAGE_SEPARATION_ANG:.0f} Ang) -- the slab's in-plane supercell may "
            "be too small, letting the adsorbate spuriously interact with its own lateral copies. "
            "Use a bigger in-plane supercell (stb-supercell) if so.", 'yellow'), f_out)

    # --- [2] CANDIDATE SITES: FINDING & COUNT ---
    print_section("[2] CANDIDATE SITES: FINDING & COUNT", f_out)
    with capture_library_noise(library_warnings, "pymatgen AdsorbateSiteFinder"):
        candidates, raw_counts, reduced_counts = find_candidate_sites(
            finder, site_types, args.height, args.symprec, pmg_structure.lattice.matrix)

    count_rows = [([st, str(raw_counts[st]), str(reduced_counts[st]),
                     str(raw_counts[st] - reduced_counts[st])], None) for st in site_types]
    if len(site_types) > 1:
        total_raw, total_reduced = sum(raw_counts.values()), sum(reduced_counts.values())
        count_rows.append((["TOTAL", str(total_raw), str(total_reduced),
                             str(total_raw - total_reduced)], 'cyan'))
    print_table(["Site type", "Raw candidates", "After symm. reduction", "Reduced away (equivalent)"],
                count_rows, f_out)

    if candidates:
        with capture_library_noise(library_warnings, "pymatgen plot_slab (numbered candidates)"):
            plot_path = os.path.join(output_root, "candidate_sites.png")
            write_candidate_sites_plot(pmg_structure, candidates, plot_path, show=args.view_plots)
        print_dual(f"  {color_text('[Saved]', 'cyan')} {plot_path} (every candidate above, numbered "
                    "with the SAME index the table below and --site-a/--site-b use)", f_out)

    if len(candidates) < 2:
        print_dual(color_text(
            f"[ERROR] Only {len(candidates)} candidate site(s) found for '{args.site_type}' -- a "
            "diffusion-path pair needs at least 2. Try --site-type all, or a bigger --symprec.",
            'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    suggested_pair = nearest_candidate_pair(candidates, pmg_structure.lattice.matrix)

    # --- [3] SITE SELECTION ---
    print_section("[3] SITE SELECTION", f_out)
    print_table(["Index", "Type", "Coordinate (Ang)"], candidate_table_rows(candidates, suggested_pair), f_out)
    if suggested_pair:
        si, sj, sdist = suggested_pair
        print_dual(f"  {color_text('[Suggested]', 'cyan')} closest candidate pair: site {si} <-> "
                    f"site {sj} ({sdist:.2f} Ang apart) -- a plausible short hop; not chosen "
                    "automatically, --site-a/--site-b is always your call.", f_out)

    if args.site_a is None:
        print_dual(color_text(
            "\n[STOP] Pass --site-a/--site-b (0-based indices from the table above) to pick the "
            "two endpoints, then re-run.", 'yellow'), f_out)
        if f_out:
            f_out.close()
        print(f"\n{len(candidates)} candidate site(s) listed above -- re-run with --site-a/--site-b.")
        sys.exit(1)

    n_candidates = len(candidates)
    for label, idx in (("--site-a", args.site_a), ("--site-b", args.site_b)):
        if not (0 <= idx < n_candidates):
            print_dual(color_text(
                f"[ERROR] {label} {idx} out of range: only {n_candidates} candidate site(s) found.",
                'red'), f_out)
            if f_out:
                f_out.close()
            sys.exit(1)

    _idx_a, type_a, coord_a = candidates[args.site_a]
    _idx_b, type_b, coord_b = candidates[args.site_b]
    chosen_dist = float(np.linalg.norm(np.array(coord_a) - np.array(coord_b)))
    print_dual(f"  Chosen: site {args.site_a} ({type_a}) <-> site {args.site_b} ({type_b}), "
                f"{chosen_dist:.2f} Ang apart.", f_out)

    # --- [4] WRITING ENDPOINT FOLDERS ---
    print_section("[4] WRITING ENDPOINT FOLDERS", f_out)

    with capture_library_noise(library_warnings, "pymatgen plot_slab (chosen sites)"):
        chosen_plot_path = os.path.join(output_root, "chosen_sites.png")
        write_chosen_sites_plot(pmg_structure, coord_a, coord_b, chosen_plot_path, show=args.view_plots)
    print_dual(f"  {color_text('[Saved]', 'cyan')} {chosen_plot_path} (top-view of exactly the two "
                "chosen endpoints, 'A'/'B' labeled)", f_out)

    endpoint_dirs = {}
    pair_id = uuid.uuid4().hex[:12]
    n_adsorbate = len(mol)
    for label, tag, coord in (("site_A", "A", coord_a), ("site_B", "B", coord_b)):
        ads_struct = finder.add_adsorbate(mol, coord)
        min_dist = min_adsorbate_slab_distance(ads_struct, n_substrate)
        if min_dist is not None and min_dist < 0.7:
            print_dual(color_text(
                f"  [WARNING] {label}: closest slab-adsorbate distance is only {min_dist:.3f} Ang "
                "-- likely overlapping atoms (height too small for this site). Check before "
                "running SIESTA.", 'yellow'), f_out)
        site_dir = os.path.join(output_root, label)
        write_reference_folder(site_dir, ads_struct, calc_text, slab_species_meta, args.pseudo_dir,
                                force_spin=args.force_spin)
        manifest_path = write_neb_manifest(site_dir, pair_id, label, args.adsorbate,
                                            n_substrate, n_adsorbate)
        print_dual(f"  {color_text('[OK]', 'green')} {site_dir}", f_out)
        print_dual(f"  {color_text('[OK]', 'green')} {manifest_path} (atom-order proof for stb-neb)",
                    f_out)
        endpoint_dirs[label] = site_dir
        view_frames.append((f"{label} ({tag})", AseAtomsAdaptor.get_atoms(ads_struct)))

    # --- [5] SUMMARY & NEXT STEPS ---
    print_section("[5] SUMMARY & NEXT STEPS", f_out)
    print_dual(f"  2 endpoint folder(s) written under '{output_root}'.", f_out)
    if report_path:
        print_dual(f"  Report          : {report_path}", f_out)
    print_dual(f"  Run SIESTA (a real relaxation, not single-point) in '{endpoint_dirs['site_A']}' "
                f"and '{endpoint_dirs['site_B']}', then:", f_out)
    print_dual(f"    stb-neb --initial {endpoint_dirs['site_A']} --final {endpoint_dirs['site_B']} "
                "-c <calc.fdf>", f_out)
    print_dual("  (stb-neb accepts a folder directly and prefers its finished siesta.XV over the "
                "pre-relaxation structure.fdf -- see stb-neb --help.)", f_out)
    print_dual("  Each folder also carries a neb_manifest.json recording the exact atom order -- "
                "stb-neb will detect it and use it as a PROVEN atom correspondence, skipping "
                "distance-based --autosort-tol matching entirely.", f_out)

    # --- [6] LIBRARY WARNINGS ---
    print_section("[6] LIBRARY WARNINGS", f_out)
    if library_warnings:
        print_dual(color_text(
            "Messages emitted by external libraries (pymatgen) during this run -- collected here "
            "instead of interleaved with the report above; harmless in almost every case, but "
            "worth a look if a section above looks suspicious.", 'cyan'), f_out)
        for entry in library_warnings:
            print_dual(entry, f_out)
    else:
        print_dual("No library warnings.", f_out)

    if f_out:
        f_out.close()

    print(f"\n{color_text('Success:', 'green')} 2 endpoint folder(s) written under '{output_root}' "
          f"(site_A={type_a}, site_B={type_b}).")
    if report_path:
        print(f"Full report: {report_path}")

    if args.view:
        print(f"\n{color_text('--view:', 'cyan')} opening {len(view_frames)} frame(s) in ASE's "
              "interactive viewer (use the frame slider/menu to step through them):")
        for i, (label, _atoms) in enumerate(view_frames):
            print(f"  {i} = {label}")
        view_structure_interactive([atoms for _label, atoms in view_frames])


if __name__ == "__main__":
    main()
