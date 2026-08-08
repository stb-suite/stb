"""Shared slab/adsorbate-site geometry helpers, extracted from adsorb.py
once neb_sites.py became a second consumer (same extract-on-second-use
policy as the rest of core/). These are the pieces of adsorb.py that are
about SLAB + CANDIDATE-SITE geometry in general -- vacuum-axis validation,
resolving an adsorbate name to a pymatgen Molecule, and writing one
fixed-cell reference folder -- not the E_ads-specific machinery (isolated
-adsorbate reference, BSSE, height sweeps, ML ranking), which stays in
adsorb.py since neb_sites.py doesn't need any of it.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from pymatgen.core import Molecule
from pymatgen.core.periodic_table import Element
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.adsorption import plot_slab
from stb.core import structure_io, kspace
from stb.core.cli import color_text
from stb.core.pseudopotentials import copy_pseudo

# Every folder written via write_reference_folder() gets this via a
# %include sidecar (same prepend_include convention as hubbardu.py/
# elastic_inputs.py/phonons_create.py) rather than editing --calc in
# place: the slab's in-plane lattice constant is set by the underlying
# bulk periodicity, so a single local adsorbate/molecule has no physical
# business changing it, and an uncontrolled cell relaxation would also
# change the vacuum spacing along c. See
# examples/4.8-adsorption/README.md for the full rationale (written for
# stb-adsorb, applies identically to neb_sites.py's site_A/site_B).
FIXED_CELL_BLOCK = (
    "# Auto-generated -- keeps the cell fixed (see examples/4.8-adsorption/\n"
    "# README.md for why).\n"
    "MD.VariableCell false\n"
)
# Opt-in, appended to config_extra.fdf only when a caller passes
# write_reference_folder(..., force_spin=True) -- a single adsorbate atom (or
# a molecule containing one, e.g. NO) bonded to a slab very commonly leaves
# the combined slab+adsorbate system with a net magnetic moment (most simply,
# an odd total valence-electron count, which a spin-restricted calculation
# cannot represent correctly at all), same reasoning adsorb.py's own
# force_spin_polarized() already applies to its isolated-adsorbate reference.
# Unlike that function (which edits calc_text's Spin tag directly), this goes
# through config_extra.fdf so it overrides any Spin tag already present in
# the caller's own --calc template (SIESTA's fdf reader is first-occurrence
# -wins for a duplicate label, and config_extra.fdf is %include'd at the very
# top of the written calc.fdf -- see structure_io.prepend_include). Costs
# nothing for a genuinely closed-shell site (converges to zero moment).
SPIN_POLARIZED_BLOCK = (
    "# Auto-generated -- forces spin-polarized SCF (see write_reference_folder's\n"
    "# force_spin docstring for why).\n"
    "Spin                polarized\n"
)
# Opt-in, appended to config_extra.fdf only when a caller passes
# write_reference_folder(..., force_dipole=True) -- adsorbing on only ONE face of
# a slab breaks whatever inversion/mirror symmetry the clean slab had along the
# surface normal, giving the cell a net dipole moment along a PERIODIC direction;
# without the correction, that spurious periodic-image field contaminates the
# total energy -- same reasoning stb-her/stb-oer already apply unconditionally to
# their own one-sided-adsorbate site folders (their own force_dipole_correction,
# which edits calc_text directly instead of going through config_extra.fdf, since
# neither tool has this file at all). Free when the true dipole is already zero
# (a symmetric/clean slab): the correction itself evaluates to zero, so there's
# no downside to applying it whenever the caller actually wants it.
DIPOLE_CORRECTION_BLOCK = (
    "# Auto-generated -- forces the slab dipole correction (see write_reference_folder's\n"
    "# force_dipole docstring for why).\n"
    "Slab.DipoleCorrection      .true.\n"
)
# Opt-in, appended to config_extra.fdf only when a caller passes
# write_reference_folder(..., force_vdw=True) -- standard GGA functionals miss
# dispersion (van der Waals) forces entirely, which matter for a physisorption
# -like or weakly-bound adsorption site (her.py/oer.py already print an
# advisory [NOTE] recommending a vdW correction for exactly this reason, but
# neither forces one). `DFTD3 .true.` is SIESTA's own Grimme DFT-D3 pairwise
# correction added on top of whatever XC functional the user's --calc template
# already declares -- confirmed against a real, successfully-completed SIESTA
# run in this exact workflow (MESSAGES: "DFT-D3: loading default parameters
# for functional PBE"), the only vdW mechanism with any real evidence of use
# anywhere in this suite (SIESTA's OTHER vdW mechanism, a native non-local
# XC.Functional VDW/XC.Authors DRSLL-style functional, has no evidence of use
# anywhere and is NOT what this forces). Uses SIESTA's own DFTD3.* defaults
# (BJ damping, PBE-parametrized dispersion coefficients) -- no sub-parameters
# forced here, same "just the boolean toggle" scope as inputfile.py/
# cohesive_energy.py's own --dispersion/--d3 flags elsewhere in this suite.
VDW_CORRECTION_BLOCK = (
    "# Auto-generated -- forces the Grimme DFT-D3 dispersion correction (see\n"
    "# write_reference_folder's force_vdw docstring for why).\n"
    "DFTD3                   .true.\n"
)
CONFIG_EXTRA_FILE = "config_extra.fdf"

_MIN_LATERAL_IMAGE_SEPARATION_ANG = 8.0

# Neighboring in-plane cell offsets to probe -- the 8 nearest (n1, n2) pairs
# around the origin. Sufficient here because a real "too small" supercell
# always shows up in one of the immediate neighbors (the adsorbate would
# have to be larger than the ENTIRE cell for a further-out image to be the
# true minimum).
_NEIGHBOR_OFFSETS = [(n1, n2) for n1 in (-1, 0, 1) for n2 in (-1, 0, 1) if (n1, n2) != (0, 0)]


def resolve_adsorbate(name):
    """Returns a pymatgen Molecule for `name`: a molecule from ASE's bundled
    G2 database (same technique as stb-molecule.py, which doesn't expose an
    importable function -- ase.build.molecule is called directly here), or a
    single chemical element (one atom at the origin). Tries G2 first since
    G2 names are compound formulas ('H2O', 'CH4', ...) that never collide
    with a bare element symbol.
    """
    from ase.collections import g2
    from ase.build import molecule as ase_build_molecule

    if name in g2.names:
        atoms = ase_build_molecule(name)
        return AseAtomsAdaptor.get_molecule(atoms)

    try:
        Element(name)
    except ValueError:
        raise ValueError(
            f"'{name}' is not a recognized element symbol or ASE G2 molecule name. "
            "Pass --list to see the available G2 molecules."
        )
    return Molecule([name], [[0.0, 0.0, 0.0]])


def reorient_vacuum_to_c(structure, vacuum_axes):
    """Permutes the lattice vectors (and the matching fractional-coordinate
    components) so the single vacuum-padded axis becomes c -- a pure
    relabeling, not a literal 3D rotation: AdsorbateSiteFinder only cares
    about which lattice VECTOR is "out of plane" (mi_vec = cross(a, b)), not
    any particular Cartesian alignment, same convention stb-slab already
    writes (SlabGenerator's reorient_lattice=True, slab.py:21-29). Restores
    a right-handed lattice (positive determinant) if the permutation itself
    would flip handedness, by swapping the two in-plane axes -- the physical
    crystal is identical either way (same "relabeling, not a different
    crystal" spirit as core/symmetry.py's reduce_to_unitcell docstring).

    `vacuum_axes` must have exactly one True entry (caller's responsibility
    -- see resolve_slab_orientation). Returns (new_structure, new_vacuum_axes).
    """
    vacuum_axis = vacuum_axes.index(True)
    order = [i for i in range(3) if i != vacuum_axis] + [vacuum_axis]
    new_lattice = structure.lattice[order]
    if np.linalg.det(new_lattice) < 0:
        order[0], order[1] = order[1], order[0]
        new_lattice = structure.lattice[order]

    positions = np.array([pos for _, pos in structure.atoms])
    is_cartesian = structure.coord_format == 'cartesian'
    frac_coords = kspace.to_fractional(positions, structure.lattice, is_cartesian)
    new_frac_coords = frac_coords[:, order]
    new_atoms = [(sym, new_frac_coords[i]) for i, (sym, _pos) in enumerate(structure.atoms)]

    new_structure = structure_io.FdfStructure(
        lattice=new_lattice,
        lattice_constant=structure.lattice_constant,
        species=structure.species,
        species_meta=structure.species_meta,
        atoms=new_atoms,
        coord_format="fractional",
        raw_lines=[],
    )
    return new_structure, [False, False, True]


def resolve_slab_orientation(structure, vacuum_gap):
    """Validates the input has a single well-defined vacuum axis, and
    reorients it to c (see reorient_vacuum_to_c) if it isn't already there.
    Exits with a clear error only for the cases that can't be auto-fixed:
    no vacuum axis at all (a bulk 3D structure, not a slab) or vacuum on 2-3
    axes (a wire or an isolated molecule, no single well-defined surface to
    adsorb onto). Raises ValueError instead of calling sys.exit() itself --
    both callers (adsorb.py, neb_sites.py) print their own [ERROR]-prefixed
    message and exit with their own tool's exact wording/report structure.
    """
    positions = np.array([pos for _, pos in structure.atoms])
    is_cartesian = structure.coord_format == 'cartesian'
    frac_coords = kspace.to_fractional(positions, structure.lattice, is_cartesian)
    vacuum_axes = kspace.detect_vacuum_axes(frac_coords, structure.lattice, vacuum_gap)
    if sum(vacuum_axes) != 1:
        detected = ', '.join(axis for axis, is_vac in zip('abc', vacuum_axes) if is_vac) or 'none'
        raise ValueError(
            "needs a slab/2D material with vacuum along exactly one axis; detected vacuum "
            f"axis/axes: {detected}. A bulk 3D structure (no vacuum) or a wire/molecule "
            "(vacuum on 2-3 axes) doesn't have a single well-defined surface to adsorb onto."
        )

    relabel_note = None
    if not vacuum_axes[2]:
        old_axis = 'abc'[vacuum_axes.index(True)]
        structure, vacuum_axes = reorient_vacuum_to_c(structure, vacuum_axes)
        relabel_note = (
            f"Input structure's vacuum axis was '{old_axis}', not c -- relabeled the lattice "
            "vectors so vacuum is on c (the convention this tool/AdsorbateSiteFinder assume). "
            "This is a pure relabeling, not a rotation: the physical structure is unchanged, "
            "but every written structure.fdf below uses this relabeled a/b/c order, not the "
            "input file's original one."
        )

    return structure, relabel_note


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

    Moved here from neb_sites.py once adsorb.py became a second consumer --
    same extract-on-second-use policy as the rest of core/.
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


def cluster_candidate_coords(coords, lattice, reference=None):
    """Re-wraps candidate Cartesian coordinates (adsorption sites, all near
    the same slab surface) into a single, mutually consistent periodic image
    -- minimum-image relative to `reference` (or coords[0] if not given) --
    instead of pymatgen's own independent per-point floor-wrap
    (AdsorbateSiteFinder's internal put_coord_inside), which can leave two
    symmetrically-equivalent candidates on opposite sides of the cell purely
    from which side of a fractional-coordinate boundary they happened to
    floor to (verified live: AdsorbateSiteFinder.symm_reduce keeps whichever
    representative of an equivalence orbit it encounters first, with no
    preference for proximity to anything else -- so a borderline atom at
    e.g. frac x=-0.0003 floors to x~0.9997, the opposite corner from an
    equivalent neighbor at x~0.0003). Same `frac_diff -= round(frac_diff)`
    minimum-image idiom as core/md_traj.py::unwrap_trajectory /
    mlphonons.py::_pbc_atomic_shift.

    NOTE: this only guarantees the points end up MUTUALLY close (all within
    half a cell of `reference`) -- `reference` itself can sit anywhere in
    [0,1), so the whole clustered group can end up straddling or sitting
    outside the conventional [0,1) cell. A follow-up attempt to also force
    the group's centroid into [0,1) via an extra whole-lattice-vector shift
    was tried and reverted (2026-08-06): it fixed that for some cases but
    regressed others (a cluster whose true position straddles the x=0/x=1
    seam cannot be BOTH mutually close AND fully inside a fixed [0,1) box --
    doing so provably requires shifting which cell is highlighted in the
    plot, not just the candidate coordinates, a bigger change deliberately
    not made). Kept simple by design: mutual closeness only.

    A pure lattice-vector translation per point -- the physical site
    (chemistry, local coordination) is completely unchanged, only which
    periodic copy of it is reported/plotted/used. Returns `coords`
    unchanged if empty.
    """
    if not coords:
        return coords
    inv_lattice = np.linalg.inv(lattice)
    ref_frac = np.asarray(reference if reference is not None else coords[0]) @ inv_lattice
    wrapped = []
    for c in coords:
        frac = np.asarray(c) @ inv_lattice
        diff = frac - ref_frac
        diff -= np.round(diff)
        wrapped.append((ref_frac + diff) @ lattice)
    return wrapped


def write_reference_folder(out_dir, pmg_structure, calc_text, species_meta, pp_path,
                            force_spin=False, force_dipole=False, force_vdw=False):
    """Writes structure.fdf + calc.fdf + config_extra.fdf + copied
    pseudopotentials for one reference/candidate folder (adsorb.py's
    clean_slab/, adsorbate/, sites/site_*/; neb_sites.py's site_A/,
    site_B/). calc.fdf is '%include config_extra.fdf' followed by the
    untouched --calc template (structure_io.prepend_include, same
    convention as hubbardu.py) -- config_extra.fdf itself always forces
    MD.VariableCell false (FIXED_CELL_BLOCK), and additionally forces
    Spin polarized (SPIN_POLARIZED_BLOCK) when `force_spin=True` -- opt-in
    per caller, since it's only clearly universally correct for a folder
    that actually places an adsorbate atom onto a slab (see
    SPIN_POLARIZED_BLOCK's own comment), not e.g. a bare clean-slab
    reference -- the slab dipole correction (DIPOLE_CORRECTION_BLOCK)
    when `force_dipole=True` -- opt-in since it's only meaningful for a
    genuine 2D-periodic-plus-vacuum slab, not e.g. an isolated molecule in
    an all-around vacuum box -- and the DFT-D3 dispersion correction
    (VDW_CORRECTION_BLOCK) when `force_vdw=True` -- unlike the other two,
    meaningful for every folder this function writes (slab, isolated
    molecule, or combined site alike), so callers that want it on
    unconditionally (stb-adsorb, see its own force_vdw=True call sites)
    just always pass True, same "hardcoded True at the call site, no
    opt-out flag" style already used there for force_dipole. Returns the
    FdfStructure written, in case a caller needs the exact same geometry
    object again without re-reading the file.
    """
    os.makedirs(out_dir, exist_ok=True)
    fdf_structure = structure_io.from_pymatgen(pmg_structure, species_meta=species_meta, coord_format="fractional")
    structure_io.write_fdf(fdf_structure, os.path.join(out_dir, "structure.fdf"))
    with open(os.path.join(out_dir, CONFIG_EXTRA_FILE), "w") as f:
        f.write(FIXED_CELL_BLOCK)
        if force_spin:
            f.write(SPIN_POLARIZED_BLOCK)
        if force_dipole:
            f.write(DIPOLE_CORRECTION_BLOCK)
        if force_vdw:
            f.write(VDW_CORRECTION_BLOCK)
    with open(os.path.join(out_dir, "calc.fdf"), "w") as f:
        f.write(structure_io.prepend_include(calc_text, CONFIG_EXTRA_FILE))
    symbols = {site.specie.symbol for site in pmg_structure}
    for sym in sorted(symbols):
        copy_pseudo(pp_path, sym, out_dir)
    return fdf_structure


def generate_systematic_orientations(mol, n_polar=1, n_azimuthal=1):
    """Returns a list of `n_polar * n_azimuthal` rotated copies of `mol` (a
    pymatgen `Molecule`, rotated about its own center of mass) systematically
    covering orientation space -- deliberately NOT random (unlike
    `mladsorb.py::pick_best_orientation`'s `Rotation.random`), so the same
    `n_polar`/`n_azimuthal` always produces the exact same orientation set,
    reproducibly, with no seed needed.

    `AdsorbateSiteFinder.add_adsorbate` always translates whichever atom(s)
    sit at the MOST NEGATIVE z of the molecule IT'S GIVEN to touch the
    surface, then aligns that same z-axis to the site's surface normal (see
    its own docstring, quoted in write_site_plot above) -- so pre-rotating
    the molecule so a DIFFERENT reference direction becomes its "most
    negative z" is exactly how to control which part of the molecule ends
    up facing the surface. `n_polar` reference directions (in the
    molecule's OWN, unrotated frame) are distributed evenly over a sphere
    via a Fibonacci-sphere lattice (deterministic, no clustering at the
    poles the way a naive lat/long grid would have) -- for a 3-atom
    molecule like H2O, different directions land the rotated molecule's O
    atom, or one of the H atoms, or a point in between, at the new
    "most-negative-z" position, which is what ends up touching the
    surface. Each of those is then combined with `n_azimuthal` evenly
    spaced spins (0 to 360 deg, exclusive) about that same new z-axis --
    the remaining physically meaningful degree of freedom once "which part
    touches down" is fixed (in-plane rotation / tilt of the rest of the
    molecule relative to the site's local symmetry).

    `n_polar=1, n_azimuthal=1` (the default) returns `[mol]` unchanged --
    the existing, single-orientation behavior every caller had before this
    function existed.
    """
    n_polar = max(1, n_polar)
    n_azimuthal = max(1, n_azimuthal)
    if n_polar == 1 and n_azimuthal == 1:
        return [mol]

    from scipy.spatial.transform import Rotation

    center = mol.center_of_mass
    centered = mol.cart_coords - center

    if n_polar == 1:
        directions = [np.array([0.0, 0.0, -1.0])]
    else:
        golden_angle = np.pi * (3.0 - np.sqrt(5.0))
        directions = []
        for i in range(n_polar):
            z = 1.0 - 2.0 * i / (n_polar - 1)
            radius = np.sqrt(max(0.0, 1.0 - z * z))
            theta = golden_angle * i
            directions.append(np.array([radius * np.cos(theta), radius * np.sin(theta), z]))

    azimuth_angles = np.linspace(0.0, 2.0 * np.pi, n_azimuthal, endpoint=False)
    down = np.array([0.0, 0.0, -1.0])

    orientations = []
    for direction in directions:
        base_rot, _rssd = Rotation.align_vectors([down], [direction])
        for phi in azimuth_angles:
            spin = Rotation.from_rotvec([0.0, 0.0, phi])
            combined = spin * base_rot
            rotated_coords = combined.apply(centered) + center
            orientations.append(Molecule(mol.species, rotated_coords))
    return orientations


def deduplicate_orientations(entries, n_substrate, rmsd_tol=0.3, energy_tol=0.01):
    """Greedy duplicate filter for a list of relaxed orientation candidates
    at the SAME site: `entries` is a list of `(energy, ase.Atoms)` already
    sorted ascending by energy (lowest/best first); returns the list of
    INDICES into `entries` to keep (always keeps index 0, the best one).

    RMSD is computed on the adsorbate atoms only (`atoms.positions[
    n_substrate:]`) -- the substrate is identical across every candidate at
    one site (same relaxation, same fixed atoms), so no Kabsch alignment is
    needed, just a direct per-atom Cartesian distance (same atom count/
    order across all candidates, guaranteed by construction: every
    candidate here comes from the SAME `AdsorbateSiteFinder.add_adsorbate`
    call shape, just a different pre-rotated molecule). A candidate is
    dropped only if it's within BOTH `rmsd_tol` (Ang) of AND `energy_tol`
    (eV) of an already-kept one -- requiring both avoids two failure modes
    a single-metric check would have: merging two genuinely different local
    minima that happen to share a similar energy by coincidence, or merging
    two near-identical geometries whose energy differs only by MACE's own
    numerical noise.

    Does not correct for periodic wraparound (no minimum-image convention)
    -- two candidates starting at the same site coordinate and relaxing
    independently aren't expected to drift across a cell boundary
    differently enough for this to matter in practice; documented as a
    known simplification, not handled.
    """
    kept_indices = []
    kept_ads_positions = []
    kept_energies = []
    for i, (energy, atoms) in enumerate(entries):
        ads_positions = atoms.get_positions()[n_substrate:]
        is_duplicate = False
        for other_energy, other_positions in zip(kept_energies, kept_ads_positions):
            if abs(energy - other_energy) >= energy_tol:
                continue
            rmsd = np.sqrt(np.mean(np.sum((ads_positions - other_positions) ** 2, axis=1)))
            if rmsd < rmsd_tol:
                is_duplicate = True
                break
        if not is_duplicate:
            kept_indices.append(i)
            kept_ads_positions.append(ads_positions)
            kept_energies.append(energy)
    return kept_indices


def wrap_markers_into_cell(cart_sites, lattice_matrix, symm_op):
    """For each candidate site (a Cartesian, unrotated 3-vector), returns the
    ROTATED xy projection (via `symm_op`, same rotation plot_slab/
    write_site_plot already apply to every marker) of its unique periodic
    image whose in-plane fractional coordinate lies in [0, 1) x [0, 1) --
    i.e. strictly inside the single unit cell plot_slab's own outline
    highlights. Solves `[alpha, beta] = inv([a2d, b2d]) @ xy` (the point's
    coordinate in the in-plane lattice-vector basis) then wraps both via
    `% 1.0`, mirroring `neb.py::wrap_into_cell`'s `frac_coords % 1.0` idiom
    but done in this plot's own 2D ROTATED frame directly (candidates here
    are bare Cartesian points, not full pymatgen Structures, so that 3D
    helper doesn't apply as-is).

    Deliberately one marker per candidate, always inside the cell -- NOT
    guaranteed mutually close for two symmetrically-equivalent candidates
    whose true position straddles a cell seam (e.g. frac x=0.001 and
    x=0.999 are both, correctly, already inside [0,1) and each stays
    exactly there, even though they're physically adjacent across the
    boundary and so can end up looking far apart, near opposite edges of
    the same drawn cell) -- cluster_candidate_coords' own docstring already
    documents this exact tradeoff as unavoidable for a single point per
    candidate; this function picks "inside the cell" as the priority,
    cluster_candidate_coords picks "mutually close" for its own (different)
    callers. Only used for this plot -- cluster_candidate_coords itself is
    unchanged and still used for the real distance/symmetry logic
    elsewhere (main()'s own site selection, the "closest candidate pair"
    suggestion), where a single canonical coordinate genuinely is needed.
    """
    a2d = np.asarray(symm_op.operate(lattice_matrix[0])[:2])
    b2d = np.asarray(symm_op.operate(lattice_matrix[1])[:2])
    inv_ab = np.linalg.inv(np.column_stack([a2d, b2d]))
    wrapped = []
    for site in cart_sites:
        xy0 = np.asarray(symm_op.operate(site)[:2])
        alpha, beta = inv_ab @ xy0
        alpha, beta = alpha % 1.0, beta % 1.0
        wrapped.append((alpha * a2d + beta * b2d).tolist())
    return wrapped


def write_site_plot(pmg_structure, out_path, show=False):
    """Saves a top-view PNG of the slab with every ontop/bridge/hollow
    adsorption site marked (pymatgen.analysis.adsorption.plot_slab for the
    slab itself, already confirmed to work on a plain Structure -- no Slab
    wrapper needed here). The adsorption-site markers come from
    AdsorbateSiteFinder.find_adsorption_sites() with pymatgen's own
    defaults, not the caller's own --height/--symprec/--site-type, so treat
    it as a quick sanity check of the site layout, not a literal preview of
    what will be written below.

    Markers are drawn HERE, manually, rather than via plot_slab's own
    `adsorption_sites=True` convenience (which builds its own internal
    AdsorbateSiteFinder and plots its raw result) -- pymatgen's own
    AdsorbateSiteFinder.symm_reduce()/put_coord_inside() wrap each candidate
    independently with no shared reference point, so two symmetrically
    -equivalent candidates can land in different periodic images and appear
    far apart on this plot despite being the same physical site (verified
    live on a real user structure). wrap_markers_into_cell() wraps every
    candidate into the SAME, single highlighted unit cell (per explicit
    request: only markers inside the cell, none of the periodic clutter a
    prior version of this function drew across the whole plotted window)
    -- see its own docstring for the one tradeoff this implies (two
    genuinely-adjacent-across-a-cell-seam candidates can still land near
    opposite edges of the cell, since each is wrapped independently).
    Marker style (color/marker/markersize/mew/zorder) copied exactly from
    plot_slab's own adsorption_sites=True code path, so the plot looks
    identical apart from this fix.

    `show=True` (--view-plots) additionally blocks on plt.show() before the
    figure is closed -- same blocking-preview convention as every other
    --view in this suite.
    """
    from pymatgen.analysis.adsorption import AdsorbateSiteFinder, get_rot

    fig, ax = plt.subplots(figsize=(6, 6))
    plot_slab(pmg_structure, ax, adsorption_sites=False)
    ads_sites = AdsorbateSiteFinder(pmg_structure).find_adsorption_sites()["all"]
    if ads_sites:
        symm_op = get_rot(pmg_structure)
        xy = wrap_markers_into_cell(ads_sites, pmg_structure.lattice.matrix, symm_op)
        if xy:
            ax.plot(*zip(*xy, strict=True), color="k", marker="x", markersize=10, mew=1,
                    linestyle="", zorder=10000)
    ax.set_title("Candidate adsorption sites (ontop/bridge/hollow)")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)


def min_adsorbate_image_distance(ads_cart_coords, slab_lattice_matrix):
    """Minimum distance between any adsorbate atom and the nearest periodic
    image of any adsorbate atom, across the SLAB's in-plane lattice vectors
    (a, b) only -- c is the vacuum/surface-normal axis, irrelevant here.
    `ads_cart_coords` is just the adsorbate molecule's own (relative)
    Cartesian coordinates -- translation-invariant (depends only on the
    lattice + the adsorbate's own internal geometry), so this can run once
    per adsorbate rather than once per site. Checks whether the slab's
    lateral supercell is large enough that the adsorbate doesn't spuriously
    interact with its own periodic copy -- the in-plane analog of
    --vacuum-box's check for an isolated reference.
    """
    a_vec, b_vec = slab_lattice_matrix[0], slab_lattice_matrix[1]
    min_dist = None
    for n1, n2 in _NEIGHBOR_OFFSETS:
        shift = n1 * a_vec + n2 * b_vec
        diffs = ads_cart_coords[:, None, :] - (ads_cart_coords[None, :, :] + shift)
        local_min = float(np.linalg.norm(diffs, axis=-1).min())
        if min_dist is None or local_min < min_dist:
            min_dist = local_min
    return min_dist
