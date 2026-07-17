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
import argparse
from fractions import Fraction
import numpy as np
from stb.core import structure_io, kspace
from stb.core.cli import color_text, show_intro, print_dual
from stb.core.deps import require_icet
from stb.core.pseudopotentials import resolve_pseudo_source, link_pseudo

REPORT_FILE = "gqca_stage1.txt"

_RUNTYPE_RE = re.compile(r'MD\.TypeOfRun\s+\S+', re.IGNORECASE)
_VARIABLECELL_RE = re.compile(r'MD\.VariableCell\s+\S+', re.IGNORECASE)


def force_full_relaxation(calc_text, is_2d):
    """Substitutes/appends 'MD.TypeOfRun CG' always, plus 'MD.VariableCell T'
    ONLY for a genuinely 3D (bulk, no vacuum axis) input. Real (not
    single-point) relaxation is structurally required here -- see
    stb-gqca --help: the 3 representative cluster structures (AA/AB/BB)
    have, in general, different equilibrium lattice constants (atomic size
    mismatch between A and B), and evaluating them at a fixed/parent
    lattice constant would contaminate the mixing energy with spurious
    lattice-mismatch strain unrelated to real chemical/configurational
    energetics -- same class of bug already caught and fixed in this
    suite's HER (H2 molecule) and OER (O*/OOH* intermediates) workflows.

    [KNOWN LIMITATION, 2D inputs] Cell relaxation for a 2D (vacuum-padded)
    input is NOT restricted to the in-plane lattice vectors here --
    MD.VariableCell is left OFF entirely for 2D inputs (positions-only
    relaxation, at the parent/input in-plane lattice constant) rather than
    risk an unconfirmed SIESTA %block Geometry.Constraints stress-
    component-fixing syntax (researched during this tool's development;
    the manual documents 'cellside'/'cellangle' keywords but historically
    marks them "not implemented yet", and no verbatim confirmed example of
    the alternative stress-tensor-component-fixing block syntax could be
    found). This means a 2D alloy's in-plane lattice mismatch between the
    two pure end-members is NOT captured -- a real accuracy limitation,
    documented here and in the Stage 1 report rather than silently
    assumed away. Revisit once SIESTA's exact 2D-cell-constraint syntax is
    confirmed against a real installation.
    """
    new_text, count = _RUNTYPE_RE.subn('MD.TypeOfRun          CG', calc_text)
    if count == 0:
        new_text += "\nMD.TypeOfRun          CG\n"
    if is_2d:
        return new_text
    new_text2, count = _VARIABLECELL_RE.subn('MD.VariableCell       T', new_text)
    if count == 0:
        new_text2 += "MD.VariableCell       T\n"
    return new_text2


def min_sublattice_distance(pmg_structure, site_indices, start_radius=6.0, max_radius=50.0):
    """Shortest periodic distance between any two sublattice sites --
    defines the ClusterSpace cutoff (a tight shell around the nearest-
    neighbor distance) without a hardcoded/lattice-specific constant.

    Uses a periodic neighbor search (pymatgen's own get_neighbors, radius
    doubled until at least one same-sublattice neighbor is found), NOT a
    plain site-to-site distance_matrix -- verified live during this
    tool's development that distance_matrix silently breaks whenever the
    sublattice has only ONE atom per primitive cell (a common case: e.g.
    rocksalt's cation sublattice, zincblende's cation/anion sublattices),
    since that atom's nearest neighbors are periodic IMAGES of itself,
    which a same-cell-only site-to-site matrix never sees at all (its
    only diagonal entry is a trivial self-distance of exactly 0, not a
    meaningful periodic nearest-neighbor distance).
    """
    radius = start_radius
    while radius <= max_radius:
        min_dist = None
        for i in site_indices:
            for neighbor in pmg_structure.get_neighbors(pmg_structure[i], r=radius):
                if neighbor.index in site_indices:
                    d = neighbor.nn_distance
                    if min_dist is None or d < min_dist:
                        min_dist = d
        if min_dist is not None:
            return min_dist
        radius *= 2
    raise ValueError(f"Could not find any two sublattice sites within {max_radius} Ang of each "
                      "other -- unexpectedly sparse --sublattice selection.")


def select_pair_orbit(cluster_space):
    """Returns (cv_index, radius, multiplicity) for the order-2 (pair)
    orbit of smallest radius. `cv_index` is the row position in
    cluster_space.to_dataframe() -- CRITICAL to get this from the
    dataframe (or get_multiplicities(), indexed identically), not from
    cluster_space.orbit_list.orbits' own list position: the latter is
    0-indexed EXCLUDING the zerolet term, while get_cluster_vector()'s
    returned array and the dataframe both INCLUDE the zerolet at position
    0 -- verified live during this tool's development that orbits[0]
    (order=1, the singlet) corresponds to dataframe/cluster-vector row 1,
    not row 0, i.e. every orbit_list.orbits position is off-by-one versus
    the cluster-vector index actually needed to read cv[cv_index] in
    build_ab_structure. Using the dataframe directly sidesteps this
    mismatch entirely.

    ClusterSpace only generates orbits over sites that actually have >1
    allowed species (i.e. the disordered sublattice) -- verified live on
    GaAs, MoS2, and a simple-cubic test fixture during this tool's
    development: every order-2 orbit found was already a sublattice-
    sublattice pair, no extra filtering needed. Smallest radius selects
    the nearest-neighbor shell if the cutoff let more than one order-2
    orbit through.
    """
    df = cluster_space.to_dataframe()
    pair_rows = df[df["order"] == 2]
    if pair_rows.empty:
        raise ValueError("No pair (order-2) orbit found on the target sublattice within the "
                          "computed cutoff -- this shouldn't happen for a valid --sublattice choice.")
    pair_rows = pair_rows.sort_values("radius")
    row = pair_rows.iloc[0]
    cv_index = int(pair_rows.index[0])
    return cv_index, float(row["radius"]), int(row["multiplicity"])


def sublattice_fraction(cluster_space, sublattice_symbol, species_b):
    """Fraction of the icet-detected PRIMITIVE cell's atoms that belong to
    the target sublattice -- needed to convert a LOCAL cluster composition
    (e.g. 1 B atom in a 2-site pair, i.e. 50% of the sublattice) into the
    all-atom fraction icet's own enumerate_structures/concentration_restrictions
    API expects (a fraction of ALL atoms, not just the disordered
    sublattice). Same Fraction-based bookkeeping pattern as
    sqs.py::minimal_scaling, adapted (not imported -- self-contained
    workflow) since this is a different consumer with a different exact
    need (a fixed 50% sublattice target, not an arbitrary user composition).
    """
    n_primitive = len(cluster_space.primitive_structure)
    sublattices = cluster_space.get_sublattices(cluster_space.primitive_structure)
    target_symbols = {sublattice_symbol, species_b}
    matching = [sl for sl in sublattices if set(sl.chemical_symbols) == target_symbols]
    if not matching:
        raise ValueError(f"Could not find the {sublattice_symbol}/{species_b} sublattice in "
                          "icet's own sublattice detection -- unexpected structure/cutoffs combination.")
    return Fraction(len(matching[0].indices), n_primitive)


_MAX_AB_CANDIDATES = 2000  # hard cap, independent of cell size -- see build_ab_structure's
                           # docstring for why open-ended cell-size expansion is unsafe here.


def build_ab_structure(cluster_space, cv_index, sub_fraction, species_b, vacuum_gap, f_out):
    """Returns (best_structure, best_deviation) -- the n_j=1 (AB)
    representative structure and its achieved deviation from the ideal
    -1.0 pair correlation (0.0 for a perfectly AB-ordered result; see
    the [IMPORTANT PHYSICS] note below for why this is not always 0 --
    `best_deviation` is persisted by main() into the Stage 1 report so
    stb-gqcaAnalysis can carry the same caveat forward to Delta_1).

    Among every symmetrically distinct ordered structure at exactly 50% B
    symmetrically distinct ordered structure at exactly 50% B (of the
    sublattice, at the SMALLEST cell size reaching that composition
    exactly) that icet's Hart-Forcade enumeration finds, picks the one
    whose cluster-vector component at `cv_index` (the pair orbit's own
    position in get_cluster_vector()'s returned array -- see
    select_pair_orbit's docstring for why this must come from the
    dataframe/multiplicities indexing, not orbit_list.orbits' own list
    position) is closest to -1.0 (the "every pair is unlike, AB-ordered"
    extreme, e.g. a rock-salt-type alternation) -- the structure that
    most purely realizes "one A, one B" in the LOCAL pair environment the
    GQCA formalism needs, as opposed to merely having the right GLOBAL
    50:50 composition (most candidates at that composition do NOT purely
    realize this -- a random 50:50 mixture would average to pair-vector
    0, not -1).

    [IMPORTANT PHYSICS] -1.0 is NOT always achievable. Verified live
    during this tool's development: a simple-cubic sublattice (bipartite
    -- 2-colorable so every nearest neighbor is unlike) reaches exactly
    -1.0 at the smallest possible cell. An FCC sublattice (e.g. rocksalt
    NaCl's cation sublattice) does NOT -- FCC's own close-packed planes
    contain triangles (an odd cycle), so no 2-coloring can make every
    nearest-neighbor pair unlike (the same geometric frustration behind
    the classic FCC antiferromagnetic Ising model never reaching perfect
    Neel order). For such lattices, this function returns the BEST
    achievable ordering at the smallest exact-composition cell, not a
    perfect one -- documented as a genuine, inherent limitation of the
    pair-cluster approximation on frustrated lattices (not a bug to fix),
    flagged via a report warning whenever the achieved value isn't -1.0.

    Deliberately searches ONLY the single smallest cell size -- expressed
    as icet's own `sizes` unit (a MULTIPLE OF THE PRIMITIVE CELL, i.e.
    `sizes=[2]` means a 2x-primitive supercell, NOT "2 atoms total";
    verified live during this tool's development after an earlier version
    of this function conflated the two and silently over-searched by a
    factor of n_primitive -- e.g. asked for an 18-atom cell on MoS2's
    3-atom primitive when a 6-atom cell was the true smallest exact-
    composition size). The correct smallest multiplier is
    `(target_b_fraction * n_primitive)`'s Fraction denominator: the
    smallest `m` for which `target_b_fraction * m * n_primitive` (the
    exact B-atom count at that size) is an integer. Rather than expanding
    to larger cells hunting for a better -1.0 match -- verified live that
    candidate counts explode combinatorially with cell size (5 -> 94 ->
    1552 -> 8000+ for NaCl's sublattice at sizes 4/8/12/16 in the old,
    over-sized units) WITHOUT reliably improving the achieved extremum
    (the best value actually got WORSE at larger sizes in that live
    test) for a frustrated lattice -- an open-ended size search is both
    unsafe (unbounded runtime) and not even reliably beneficial.
    `_MAX_AB_CANDIDATES` is a hard per-call safety cap independent of
    cell size, for lattices where even the smallest exact-composition
    cell is itself large (an unusual sublattice-fraction edge case, not
    exercised in this tool's own development but capped defensively
    regardless).

    [2D inputs] icet's enumeration is pbc-unaware in this suite's actual
    pipeline BY DEFAULT (a pymatgen-derived ase.Atoms always carries
    pbc=[True,True,True], even for a vacuum-padded 2D structure -- see
    detect_structure_vacuum's own docstring), so an unrestricted search
    can tile the supercell PARTLY ALONG the vacuum direction instead of
    purely in-plane -- verified live that this produces a physically
    nonsensical result (a MoS2 vacuum thickness ballooning from ~33 Ang
    to ~99 Ang, with the vacuum lattice vector no longer axis-aligned).
    Fixed here by setting `pbc=[True,True,False]` (vacuum axis only) on
    a COPY of the primitive structure before calling enumerate_structures
    -- confirmed live this restricts Hart-Forcade to the correct, smaller
    2D search space (1 candidate at cell size 2x primitive for MoS2,
    down from 3 in an earlier unrestricted-plus-post-filter version of
    this function) and returns an exactly axis-aligned, un-reshaped
    vacuum vector (unlike Niggli-reduced 3D-search results, whose vacuum
    vector can come out tilted even when its LENGTH is preserved). An
    earlier attempt at this exact fix appeared to fail (zero candidates
    at every size up to 24x primitive) -- root cause turned out to be
    the sizes-unit bug described above (the size passed was ~3x too
    large for MoS2's primitive, an infeasible request under the pbc
    restriction), not the pbc restriction itself; re-verified live after
    fixing the sizes calculation that pbc restriction alone (no
    niggli_reduce override needed -- icet's own default already turns
    Niggli reduction off whenever pbc isn't True on all 3 axes) works
    correctly.

    Enumerates from cluster_space.primitive_structure, NOT the raw input
    atoms -- verified live during this tool's development that
    enumerate_structures requires its `structure` argument to already be
    exactly as large as `chemical_symbols` (icet's own primitive-reduced
    list, via cluster_space.chemical_symbols); passing a non-primitive
    input structure (e.g. the user's own supercell) alongside the
    primitive-sized chemical_symbols list raises an opaque icet-internal
    error instead of a clear one.
    """
    from icet.tools import enumerate_structures

    prim = cluster_space.primitive_structure
    n_primitive = len(prim)
    target_b_fraction = Fraction(1, 2) * sub_fraction
    size_fraction = target_b_fraction * n_primitive
    size_multiplier = size_fraction.denominator

    prim_lattice = np.array(prim.cell[:])
    prim_frac = prim.get_scaled_positions()
    vacuum_axes = kspace.detect_vacuum_axes(prim_frac, prim_lattice, vacuum_gap)
    is_2d = any(vacuum_axes)

    search_structure = prim
    if is_2d:
        search_structure = prim.copy()
        search_structure.pbc = [not v for v in vacuum_axes]

    best_structure, best_deviation = None, None
    n_evaluated = 0
    for candidate in enumerate_structures(
            structure=search_structure, sizes=[size_multiplier],
            chemical_symbols=cluster_space.chemical_symbols,
            concentration_restrictions={species_b: (float(target_b_fraction), float(target_b_fraction))}):
        cv = cluster_space.get_cluster_vector(candidate)
        deviation = abs(cv[cv_index] - (-1.0))
        if best_deviation is None or deviation < best_deviation:
            best_structure, best_deviation = candidate, deviation
        n_evaluated += 1
        if n_evaluated >= _MAX_AB_CANDIDATES:
            print_dual(color_text(
                f"  [WARNING] Hit the {_MAX_AB_CANDIDATES}-candidate safety cap while searching "
                f"for the AB structure at cell size {size_multiplier}x primitive -- using the best "
                "match found so far, which may not be the true optimum at this size.", 'yellow'), f_out)
            break

    if best_structure is None:
        raise ValueError("icet's enumeration found no candidate structure at 50% B composition "
                          "on the target sublattice -- unexpected for a valid --sublattice/--species-b.")
    if best_deviation >= 1e-6:
        print_dual(color_text(
            f"  [NOTE] The best achievable AB-ordered pair correlation at cell size "
            f"{size_multiplier}x primitive deviates by {best_deviation:.4f} from the ideal -1.0 "
            "(every-pair-unlike) extreme -- expected on a geometrically frustrated sublattice "
            "(e.g. FCC/HCP/triangular; see build_ab_structure's own docstring), not necessarily a "
            f"bug. Evaluated {n_evaluated} candidate(s).", 'yellow'), f_out)
        _report_larger_cell_sensitivity(search_structure, cluster_space, cv_index, species_b,
                                         target_b_fraction, size_multiplier, best_deviation, f_out)
    return best_structure, best_deviation


_SENSITIVITY_PROBE_MAX_CANDIDATES = 500  # small/cheap on purpose -- an informational sensitivity
                                          # check, not a second exhaustive search (see
                                          # _report_larger_cell_sensitivity's own docstring for
                                          # why an open-ended search here is unsafe, same
                                          # combinatorial-explosion risk as build_ab_structure's
                                          # own size search).


def _report_larger_cell_sensitivity(search_structure, cluster_space, cv_index, species_b,
                                     target_b_fraction, size_multiplier, best_deviation, f_out):
    """[Delta_1 sensitivity check] Only called when the smallest cell's
    AB-ordering deviation is non-negligible (a geometrically frustrated
    sublattice). Runs ONE bounded, cheap probe at the next cell size
    (2x the smallest exact-composition size -- always still exactly
    reaches the target composition, since size_multiplier was chosen as
    the SMALLEST integer solution, so any integer multiple of it is also
    exact) to check whether a larger cell could realize a meaningfully
    better AB ordering, capped at `_SENSITIVITY_PROBE_MAX_CANDIDATES`
    (much smaller than build_ab_structure's own `_MAX_AB_CANDIDATES` --
    this is a quick sanity probe, not a real search, and deliberately
    stays bounded rather than expanding further: candidate counts explode
    combinatorially with cell size for frustrated lattices, the same
    runaway-search risk documented at length in build_ab_structure's own
    docstring). Purely informational -- never changes which structure
    Stage 1 actually writes to cluster_n1/ (that stays the smallest-cell
    result, for consistency/reproducibility and to keep Stage 1's runtime
    bounded and predictable regardless of what this probe finds).
    """
    from icet.tools import enumerate_structures

    probe_size = size_multiplier * 2
    best_probe_deviation = None
    n_probed = 0
    for candidate in enumerate_structures(
            structure=search_structure, sizes=[probe_size],
            chemical_symbols=cluster_space.chemical_symbols,
            concentration_restrictions={species_b: (float(target_b_fraction), float(target_b_fraction))}):
        cv = cluster_space.get_cluster_vector(candidate)
        deviation = abs(cv[cv_index] - (-1.0))
        if best_probe_deviation is None or deviation < best_probe_deviation:
            best_probe_deviation = deviation
        n_probed += 1
        if n_probed >= _SENSITIVITY_PROBE_MAX_CANDIDATES:
            break

    if best_probe_deviation is None:
        print_dual(color_text(
            f"  [NOTE] Delta_1 sensitivity check: no candidate found at the next cell size "
            f"({probe_size}x primitive) within {_SENSITIVITY_PROBE_MAX_CANDIDATES} probes -- "
            "inconclusive, not necessarily a problem.", 'cyan'), f_out)
        return

    improvement = best_deviation - best_probe_deviation
    if improvement > 0.02:
        print_dual(color_text(
            f"  [NOTE] Delta_1 sensitivity check: a {probe_size}x-primitive cell achieves a "
            f"meaningfully better AB ordering (deviation {best_probe_deviation:.4f} vs. "
            f"{best_deviation:.4f} at the size actually used, from {n_probed} probed candidate(s)) "
            "-- if Delta_1's accuracy matters for your use case, consider manually building and "
            "relaxing a larger cluster_n1 cell yourself (Stage 1 always uses the smallest exact-"
            "composition cell, for predictable runtime).", 'cyan'), f_out)
    else:
        print_dual(color_text(
            f"  [NOTE] Delta_1 sensitivity check: a {probe_size}x-primitive cell does not "
            f"meaningfully improve the AB ordering (deviation {best_probe_deviation:.4f} vs. "
            f"{best_deviation:.4f} at the size actually used, from {n_probed} probed candidate(s)) "
            "-- the deviation at the smallest cell already appears representative, not an "
            "artifact of under-sizing.", 'cyan'), f_out)


def build_pure_structure(base_fdf_structure, site_indices, symbol):
    """Trivial n_j=0/n_j=2 representative: every sublattice site set to
    one symbol, everything else (lattice, other species) unchanged.
    """
    new_atoms = [(symbol if i in site_indices else sym, pos)
                 for i, (sym, pos) in enumerate(base_fdf_structure.atoms)]
    present = {sym for sym, _ in new_atoms}
    species_meta = dict(base_fdf_structure.species_meta)
    species_meta = structure_io.ensure_species_id(species_meta, symbol)
    species_meta = {k: v for k, v in species_meta.items() if k in present}
    species = [s for s in list(base_fdf_structure.species) + [symbol] if s in present]
    species = list(dict.fromkeys(species))
    return structure_io.FdfStructure(
        lattice=base_fdf_structure.lattice, lattice_constant=base_fdf_structure.lattice_constant,
        species=species, species_meta=species_meta, atoms=new_atoms,
        coord_format=base_fdf_structure.coord_format, raw_lines=[],
    )


def detect_structure_vacuum(fdf_structure, vacuum_gap):
    """Detects (vacuum_axes, is_2d) FRESH from `fdf_structure`'s own
    lattice/coordinates -- NOT reused from the original input structure's
    own detection. Important for the AB (n_j=1) structure specifically:
    it's built from icet's own primitive-cell reduction (optionally
    Niggli-reduced), which can reorder or reshape the lattice vectors
    relative to the original input -- verified during this tool's
    development that icet's enumerate_structures defaults to Niggli
    reduction whenever ASE reports the structure as periodic in all 3
    directions (true even for a vacuum-padded 2D structure, since ASE
    Atoms objects converted from a pymatgen Structure always carry
    pbc=[True,True,True] regardless of physical dimensionality). Reusing
    a vacuum-axis mask computed from the ORIGINAL structure's axis order
    against a possibly-reordered AB cell would silently apply the k-grid's
    vacuum exclusion to the wrong axis.
    """
    positions = np.array([pos for _, pos in fdf_structure.atoms])
    is_cartesian = fdf_structure.coord_format == 'cartesian'
    frac_coords = kspace.to_fractional(positions, fdf_structure.lattice, is_cartesian)
    vacuum_axes = kspace.detect_vacuum_axes(frac_coords, fdf_structure.lattice, vacuum_gap)
    return vacuum_axes, any(vacuum_axes)


def write_cluster_folder(out_dir, fdf_structure, calc_text_template, pp_path, kdensity, vacuum_gap):
    """Writes structure.fdf + calc.fdf (relaxation directives forced,
    based on THIS structure's own freshly-detected dimensionality -- see
    detect_structure_vacuum) + linked pseudos + a density-based k-grid
    (recomputed per folder -- every n_j structure has a different cell,
    and potentially a different axis order) for one representative-
    cluster folder. Returns (divisions, is_2d).
    """
    vacuum_axes, is_2d = detect_structure_vacuum(fdf_structure, vacuum_gap)
    calc_text = force_full_relaxation(calc_text_template, is_2d)

    os.makedirs(out_dir, exist_ok=True)
    structure_io.write_fdf(fdf_structure, os.path.join(out_dir, "structure.fdf"))
    with open(os.path.join(out_dir, "calc.fdf"), "w") as f:
        f.write(calc_text)
    symbols = sorted({sym for sym, _ in fdf_structure.atoms})
    for sym in symbols:
        link_pseudo(pp_path, sym, out_dir)
    divisions = kspace.compute_monkhorts(
        fdf_structure.lattice[0], fdf_structure.lattice[1], fdf_structure.lattice[2],
        kdensity, vacuum_axes=vacuum_axes)
    with open(os.path.join(out_dir, "calc.fdf"), "a") as f:
        f.write(f"\nkgrid.MonkhorstPack   [{divisions[0]}  {divisions[1]}  {divisions[2]}]\n")
    return divisions, is_2d


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Stage 1 of 2: generates the 3 ordered pair-cluster-representative "
        "structures (AA/AB/BB) for a substitutional alloy A(1-x)Bx.", 'bold')}
Implements the pair-cluster (N=2) special case of the Generalized Quasi-Chemical Approximation
(GQCA, Sher, van Schilfgaarde, Chen, Chen, Phys. Rev. B 36, 4279, 1987) -- the classic Guggenheim
quasi-chemical approximation, applicable to ANY lattice/dimensionality (2D or 3D alike; verified on
both a 3D zincblende cell and a 2D transition-metal-dichalcogenide cell during this tool's
development). Disorders one existing sublattice (--sublattice) with a second species (--species-b,
typically absent from the input structure) and writes 3 folders: the two pure end-members
(cluster_n0/, cluster_n2/) and the ordered 50:50 compound that most purely realizes an "every
neighbor pair is unlike" (AB) local environment along the shortest sublattice-sublattice bond
(cluster_n1/) -- NOT a random 50:50 mixture, whose bonds would average to a mix of AA/AB/BB.

All 3 folders are written for a REAL relaxation (MD.TypeOfRun CG always; MD.VariableCell T too for
a 3D/bulk input) -- the 3 representative structures generally have different equilibrium lattice
constants (atomic size mismatch between A and B), and single-point evaluation at a fixed lattice
constant would contaminate the mixing energy with spurious strain. [KNOWN LIMITATION] For a 2D
(vacuum-padded) input, cell relaxation is left OFF entirely (positions-only, at the input's own
in-plane lattice constant) rather than risk an unconfirmed SIESTA cell-constraint syntax for
excluding the vacuum axis -- see force_full_relaxation's docstring. Doesn't run SIESTA -- run each
folder yourself, then use stb-gqcaAnalysis.

[WHY ONLY 3 STRUCTURES / ONLY EVER 50:50] The pair cluster (N=2) has exactly 3 possible local
compositions (0, 1, or 2 B atoms per pair) -- there is no "user-chosen x" knob here, and
cluster_n1/ is always built at the cluster-internal 50:50 composition (the ONLY interior point a
2-site cluster has), regardless of the alloy composition you actually care about. Your real target
composition x (e.g. Mo0.7W0.3S2) is swept entirely in Stage 2 (stb-gqcaAnalysis --x-min/--x-max),
which combines these 3 fixed reference energies via the mass-action law to get x_j(x,T) at ANY x in
[0,1] -- these 3 DFT calculations are the complete, sufficient input for the whole x/T thermodynamic
surface, not a coarse sample of it.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage example:\n"
               "  %(prog)s -f mos2.fdf --sublattice Mo --species-b W -c calc.fdf -p dojo\n"
    )

    parser.add_argument("-f", "--file", type=str, required=True,
                         help="Input structure.fdf (2D or 3D).")
    parser.add_argument("--sublattice", type=str, required=True,
                         help="Existing species whose sites become the disordered alloy sublattice.")
    parser.add_argument("--species-b", type=str, required=True,
                         help="Second alloy species (typically absent from the input structure).")
    parser.add_argument("-c", "--calc", type=str, required=True,
                         help="calc.fdf template for the 3 cluster relaxations.")
    parser.add_argument("-p", "--pseudo-dir", type=str, default="",
                         help="Pseudopotentials source: a bundled bank or a folder path.")
    parser.add_argument("--kdensity", type=float, default=0.04,
                         help="k-point density for the per-folder Monkhorst-Pack grid, same "
                              "convention as stb-kgrid (default: 0.04).")
    parser.add_argument("--vacuum-gap", type=float, default=10.0,
                         help="Vacuum-axis detection threshold in Ang (default: 10.0).")
    parser.add_argument("-O", "--output-dir", type=str, default="gqca_study",
                         help="Root directory for the whole GQCA workflow (default: gqca_study).")
    parser.add_argument("-v", "--version", action="version", version=f"stb-gqca {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("GQCA WORKFLOW -- STAGE 1: CLUSTER STRUCTURE GENERATION", 'bold'))
    print("-" * 60)

    if not os.path.exists(args.file):
        print(color_text(f"[ERROR] Structure file '{args.file}' not found.", 'red'))
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

    if args.sublattice == args.species_b:
        print(color_text("[ERROR] --sublattice and --species-b must be different species.", 'red'))
        sys.exit(1)

    require_icet()
    from pymatgen.io.ase import AseAtomsAdaptor

    structure = structure_io.read_fdf(args.file)
    pmg_structure = structure_io.to_pymatgen(structure)
    site_indices = [i for i, site in enumerate(pmg_structure) if site.specie.symbol == args.sublattice]
    if not site_indices:
        print(color_text(f"[ERROR] Species '{args.sublattice}' not found in the structure.", 'red'))
        sys.exit(1)

    other_symbols = {site.specie.symbol for i, site in enumerate(pmg_structure) if i not in site_indices}
    if args.species_b in other_symbols:
        print(color_text(
            f"[WARNING] --species-b '{args.species_b}' is already present elsewhere in the "
            "structure (not on --sublattice) -- the alloy substitution is still well-defined, but "
            "double-check this is intentional.", 'yellow'))

    input_vacuum_axes, input_is_2d = detect_structure_vacuum(structure, args.vacuum_gap)

    with open(args.calc) as f:
        calc_text = f.read()

    output_root = args.output_dir
    os.makedirs(output_root, exist_ok=True)
    report_path = os.path.join(output_root, REPORT_FILE)

    with open(report_path, 'w') as f_out:
        print_dual(f"{color_text('===== GQCA STAGE 1 REPORT (CLUSTER STRUCTURE GENERATION) =====', 'magenta')}", f_out)

        print_dual(f"\n{color_text('[0] RUN METADATA', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        print_dual(f"Structure       : {args.file}", f_out)
        print_dual(f"Sublattice      : {args.sublattice} ({len(site_indices)} site(s))", f_out)
        print_dual(f"Species B       : {args.species_b}", f_out)
        print_dual(f"Calc template   : {args.calc}", f_out)
        print_dual(f"Output root     : {output_root}", f_out)

        print_dual(f"\n{color_text('[1] DIMENSIONALITY', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        print_dual(f"Detected (input) : {kspace.dimensionality_label(input_vacuum_axes)}", f_out)
        if input_is_2d:
            print_dual(color_text(
                "[KNOWN LIMITATION] 2D input -- cell relaxation left OFF (positions-only, at each "
                "structure's own in-plane lattice constant). See stb-gqca --help / "
                "force_full_relaxation's docstring for why. Each of the 3 cluster folders below has "
                "its own dimensionality re-detected independently (the AB structure's cell can be "
                "reshaped by icet's own primitive-cell reduction).", 'yellow'), f_out)
        else:
            print_dual("Cell treatment   : MD.VariableCell T (full cell relaxes)", f_out)

        ase_atoms = AseAtomsAdaptor.get_atoms(pmg_structure)
        from icet import ClusterSpace
        try:
            min_dist = min_sublattice_distance(pmg_structure, site_indices)
            cutoff = min_dist * 1.05
            chemical_symbols = [
                [args.sublattice, args.species_b] if i in site_indices else [site.specie.symbol]
                for i, site in enumerate(pmg_structure)
            ]
            cluster_space = ClusterSpace(structure=ase_atoms, cutoffs=[cutoff], chemical_symbols=chemical_symbols)
            cv_index, orbit_radius, orbit_multiplicity = select_pair_orbit(cluster_space)
        except ValueError as e:
            print_dual(color_text(f"[ERROR] {e}", 'red'), f_out)
            sys.exit(1)

        print_dual(f"\n{color_text('[2] PAIR ORBIT', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        print_dual(f"Shortest sublattice-sublattice distance : {min_dist:.4f} Ang", f_out)
        print_dual(f"Orbit radius    : {orbit_radius:.4f} Ang", f_out)
        print_dual(f"Multiplicity    : {orbit_multiplicity} "
                    "(icet-derived, informational -- degeneracy g_j=[1,2,1] used by the solver "
                    "is the trivial pair-cluster combinatorics, not this multiplicity)", f_out)

        try:
            sub_fraction = sublattice_fraction(cluster_space, args.sublattice, args.species_b)
        except ValueError as e:
            print_dual(color_text(f"[ERROR] {e}", 'red'), f_out)
            sys.exit(1)

        print_dual(f"\n{color_text('[3] CLUSTER STRUCTURES', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        print_dual("cluster_n1 (AB) is always built at the fixed cluster-internal 50:50 "
                    "composition (the pair's only interior point) -- your real target alloy "
                    "composition x is swept entirely in Stage 2, not chosen here. See stb-gqca "
                    "--help for why these 3 fixed references are sufficient.", f_out)

        report_rows = []
        for n_j, symbol in ((0, args.sublattice), (2, args.species_b)):
            fdf_struct = build_pure_structure(structure, site_indices, symbol)
            out_dir = os.path.join(output_root, f"cluster_n{n_j}")
            divisions, folder_is_2d = write_cluster_folder(
                out_dir, fdf_struct, calc_text, args.pseudo_dir, args.kdensity, args.vacuum_gap)
            print_dual(f"  {color_text('[OK]', 'green')} {out_dir} ({len(fdf_struct.atoms)} atoms, "
                        f"k-grid {divisions}, {'2D positions-only' if folder_is_2d else '3D full-cell'} "
                        "relaxation)", f_out)
            report_rows.append((n_j, out_dir, len(fdf_struct.atoms), divisions))

        try:
            ab_atoms, ab_deviation = build_ab_structure(cluster_space, cv_index, sub_fraction,
                                                          args.species_b, args.vacuum_gap, f_out)
        except ValueError as e:
            print_dual(color_text(f"[ERROR] {e}", 'red'), f_out)
            sys.exit(1)
        ab_pmg = AseAtomsAdaptor.get_structure(ab_atoms)
        species_meta = structure_io.species_dict(structure)
        for sym in (args.sublattice, args.species_b):
            species_meta = structure_io.ensure_species_id(species_meta, sym)
        ab_fdf = structure_io.from_pymatgen(ab_pmg, species_meta=species_meta, coord_format="fractional")
        out_dir = os.path.join(output_root, "cluster_n1")
        divisions, folder_is_2d = write_cluster_folder(
            out_dir, ab_fdf, calc_text, args.pseudo_dir, args.kdensity, args.vacuum_gap)
        print_dual(f"  {color_text('[OK]', 'green')} {out_dir} ({len(ab_fdf.atoms)} atoms, "
                    f"k-grid {divisions}, {'2D positions-only' if folder_is_2d else '3D full-cell'} "
                    "relaxation, AB-ordered)", f_out)
        report_rows.append((1, out_dir, len(ab_fdf.atoms), divisions))

        print_dual(f"\n{color_text('[4] SUMMARY & NEXT STEPS', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        print_dual(f"3 cluster folder(s) written under '{output_root}'.", f_out)
        print_dual(f"Report               : {report_path}", f_out)
        print_dual(color_text("\nNext steps:", 'yellow'), f_out)
        print_dual(f"  1. Run SIESTA in every '{output_root}/cluster_n*/' folder (a full relaxation).", f_out)
        print_dual(f"  2. Once they're done, run: stb-gqcaAnalysis --directory {output_root}", f_out)

        f_out.write("\nSublattice     : " + args.sublattice + "\n")
        f_out.write("Species B      : " + args.species_b + "\n")
        f_out.write(f"AB deviation   : {ab_deviation:.6f}\n")

    print("\n[INFO] Complete job!")
    print("\n" + "-" * 60)
    print(color_text("Cluster structure folders ready for Stage 2 (stb-gqcaAnalysis).\n", 'bold'))


if __name__ == "__main__":
    main()
