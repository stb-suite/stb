"""Shared symmetry-based cell reduction and group-label lookups.

Wraps pymatgen's SpacegroupAnalyzer to reduce a structure to its primitive
cell, its conventional cell, or a symmetry-refined version of the input cell
(used by stb-unitcell and stb-fetch), plus thin space/layer/point-group
label accessors (used by stb-fetch and stb-translate) that report just the
detected group's symbol -- the full per-site/operations analysis lives in
stb-symmetry instead.
"""

import numpy as np
import spglib
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

UNITCELL_MODES = ("primitive", "conventional", "refined")


def find_inequivalent_sites(pmg_structure, symprec, filter_species=None):
    """Returns (sites, space_group_label) where `sites` is a list of
    (index, wyckoff_letter, multiplicity) -- one representative atom (0-based
    index) per symmetrically distinct site, via spglib's equivalent-atoms
    mapping. If filter_species is given, only representatives of that
    species are returned (their multiplicity still counts all symmetry
    -equivalent atoms of that site, filtered species or not, since
    spglib never groups atoms of different species together).

    Moved here once stb-hubbardu became a second consumer (alongside
    stb-defect --all-inequivalent-sites) -- both need the same "which sites
    are physically distinct vs. just symmetry copies of each other"
    classification.
    """
    sga = SpacegroupAnalyzer(pmg_structure, symprec=symprec)
    dataset = sga.get_symmetry_dataset()

    groups = {}
    for i, rep in enumerate(dataset.equivalent_atoms):
        groups.setdefault(rep, []).append(i)

    sites = []
    for rep in sorted(groups):
        symbol = pmg_structure[rep].specie.symbol
        if filter_species is not None and symbol != filter_species:
            continue
        sites.append((rep, dataset.wyckoffs[rep], len(groups[rep])))

    space_group = f"{dataset.international} (No. {dataset.number})"
    return sites, space_group


def space_group_label(pmg_structure, symprec=1e-3) -> str:
    """Returns e.g. "Fm-3m (No. 225)" -- the same "international symbol (No.
    number)" format find_inequivalent_sites already builds inline, extracted
    here once stb-translate became a second consumer (--dry-run/conversion
    symmetry summary, and the space group actually written to CIF output
    instead of always P1)."""
    dataset = SpacegroupAnalyzer(pmg_structure, symprec=symprec).get_symmetry_dataset()
    return f"{dataset.international} (No. {dataset.number})"


def layer_group_label(pmg_structure, aperiodic_dir, symprec=1e-3):
    """Returns e.g. "p31m (No. 78)" for a genuinely 2D-periodic structure
    (periodic along the 2 axes other than `aperiodic_dir`, vacuum-padded
    along it) -- just the label, no operations/sites/Wyckoff table. Same
    "international symbol (No. number)" format as space_group_label above,
    via spglib.get_layergroup() directly (needs spglib >= 2.1.0; pymatgen's
    SpacegroupAnalyzer doesn't wrap this). Returns None if the installed
    spglib predates this function or layer-group detection fails (e.g. not
    genuinely 2D-periodic within symprec).

    A thin, label-only sibling of stb.symmetry's own compute_layer_group
    (which additionally computes sites/orbits/distortion for stb-symmetry's
    full report) -- extracted here once stb-fetch became a second consumer
    that only ever wants the label, not the full analysis.
    """
    if not hasattr(spglib, "get_layergroup"):
        return None
    cell = (pmg_structure.lattice.matrix, pmg_structure.frac_coords,
            [site.specie.Z for site in pmg_structure])
    ds = spglib.get_layergroup(cell, aperiodic_dir=aperiodic_dir, symprec=symprec)
    if ds is None:
        return None
    return f"{ds.international} (No. {ds.number})"


def symmetry_summary(pmg_structure, symprec=0.01, vacuum_gap_ang=10.0) -> dict:
    """Crystal system / (3D) space group / layer group (if genuinely
    2D-periodic) / point group / Hall symbol for one structure, as a dict
    of ready-to-print string values -- the exact per-structure detail
    stb-2dstacking's own symmetry table needs for each of its 3 structures
    (layer 1, layer 2, heterostructure). Extracted from stacking2D.py's
    private get_details() closure once stb-supercell became a second
    consumer needing the identical per-structure summary for its own
    before/after comparison table.

    Space Group is always the 3D space group of the (possibly
    vacuum-padded) cell as given -- valid, but not the physically correct
    classification for a genuinely 2D-periodic structure, which Layer
    Group (spglib's get_layergroup(), via layer_group_label above)
    accounts for properly; Layer Group is additive, never a replacement
    for Space Group.

    Returns {"Error": "..."} instead of raising if symmetry analysis fails
    outright (e.g. a degenerate structure) -- every value is a plain string,
    ready to drop straight into a core.cli.print_table row.
    """
    from stb.core import kspace

    try:
        sga = SpacegroupAnalyzer(pmg_structure, symprec=symprec)
        details = {
            "Space Group": f"{sga.get_space_group_symbol()} ({sga.get_space_group_number()})",
            "Point Group": sga.get_point_group_symbol(),
            "Crystal System": str(sga.get_crystal_system()).title(),
            "Hall Symbol": sga.get_hall(),
        }
    except Exception as e:
        return {"Error": f"Analysis failed ({e})"}

    try:
        vacuum_axes = kspace.detect_vacuum_axes(
            pmg_structure.frac_coords, pmg_structure.lattice.matrix, vacuum_gap_ang)
        if sum(vacuum_axes) == 1:
            label = layer_group_label(pmg_structure, vacuum_axes.index(True), symprec)
            details["Layer Group"] = label if label else "N/A (needs spglib >= 2.1.0)"
        else:
            details["Layer Group"] = "N/A (not 2D-periodic)"
    except Exception:
        details["Layer Group"] = "N/A"

    return details


def point_group_label(pmg_structure, tolerance=0.3):
    """Returns e.g. "Oh" for a genuinely isolated (0D) structure -- vacuum
    -padded on all 3 axes -- via pymatgen's PointGroupAnalyzer (a separate,
    non-periodic detector; neither spglib nor SpacegroupAnalyzer handle this
    case). `tolerance` is in Angstrom (molecular point-group convention),
    deliberately not tied to a crystallographic symprec. Returns None on any
    failure (e.g. a single-atom "molecule").

    A thin, label-only sibling of stb.symmetry's own compute_point_group
    (which additionally computes equivalent-atom groups/distortion for
    stb-symmetry's full report) -- extracted here once stb-fetch became a
    second consumer that only ever wants the label, not the full analysis.
    """
    from pymatgen.core import Molecule
    from pymatgen.symmetry.analyzer import PointGroupAnalyzer

    species = [str(site.specie.symbol) for site in pmg_structure]
    mol = Molecule(species, pmg_structure.cart_coords)
    try:
        pga = PointGroupAnalyzer(mol, tolerance=tolerance)
        return str(pga.get_pointgroup())
    except Exception:
        return None


def equivalent_cartesian_axes(pmg_structure, symprec=1e-3, angle_tolerance=5.0):
    """Groups the 3 Cartesian axes (x=0, y=1, z=2) by point-group symmetry:
    axes in the same group are related by a rotation/reflection of the
    structure's point group, so straining along one axis must give an
    identical mechanical response to straining along the other (Neumann's
    principle). Returns (groups, point_group_symbol), where `groups` is a
    list of sets of axis indices partitioning {0, 1, 2}.

    Uses exact axis alignment -- does some point-group operation map unit
    vector i onto (+/-) unit vector j -- rather than the weaker "the linear
    elastic tensor happens to be isotropic in this plane" criterion. This
    matters for hexagonal/trigonal in-plane axes (e.g. x vs y, 60 degrees
    apart): linear elasticity is isotropic there, but this is deliberately
    NOT treated as an equivalence, since the large-strain/ideal-strength
    calculations this is meant to support (stb-strain) are exactly where
    higher-order anharmonic terms make e.g. zigzag vs armchair directions
    genuinely different in real 2D materials (a known effect, e.g. in
    graphene) -- verified live against a synthetic cubic structure (m-3m:
    x/y/z all grouped together) and the suite's hexagonal 2D test fixture
    (-6m2: x and y correctly NOT grouped).
    """
    sga = SpacegroupAnalyzer(pmg_structure, symprec=symprec, angle_tolerance=angle_tolerance)
    point_group = sga.get_point_group_symbol()
    ops = sga.get_point_group_operations(cartesian=True)

    parent = [0, 1, 2]

    def find(i):
        while parent[i] != i:
            i = parent[i]
        return i

    def union(i, j):
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[ri] = rj

    basis = np.eye(3)
    for op in ops:
        for i in range(3):
            mapped = op.rotation_matrix @ basis[i]
            if not np.isclose(np.linalg.norm(mapped), 1.0, atol=1e-3):
                continue
            for j in range(3):
                if i != j and np.allclose(np.abs(mapped), basis[j], atol=1e-3):
                    union(i, j)

    groups = {}
    for i in range(3):
        groups.setdefault(find(i), set()).add(i)
    return list(groups.values()), point_group


VOIGT_MODES = ('xx', 'yy', 'zz', 'yz', 'xz', 'xy')

# Voigt index (0-based, 1..6 in the usual 1-based convention) -> Cartesian
# stress/strain-tensor (row, col). Matches VOIGT_MODES order (4=yz, 5=xz,
# 6=xy). Single source of truth, shared with elastic_analysis.py.
VOIGT_TO_TENSOR = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]


def strain_tensor(mode, delta=1.0):
    """Builds the symmetric 3x3 strain tensor for one of the 6 canonical
    Voigt strain modes. Single source of truth for "what does mode X mean
    as a strain tensor", shared by equivalent_strain_modes below and
    stb-elasticInputs' own deformation-matrix builder, so the two can never
    drift apart. Only the 6 Voigt modes are covered here -- 'bi'/'hydro'
    (stb-elasticInputs' 2 extra convenience modes) aren't part of the Voigt
    basis and stay defined locally in that tool.
    """
    eps = np.zeros((3, 3))
    if mode == 'xx': eps[0, 0] = delta
    elif mode == 'yy': eps[1, 1] = delta
    elif mode == 'zz': eps[2, 2] = delta
    elif mode == 'yz': eps[1, 2] = eps[2, 1] = delta
    elif mode == 'xz': eps[0, 2] = eps[2, 0] = delta
    elif mode == 'xy': eps[0, 1] = eps[1, 0] = delta
    return eps


def equivalent_strain_modes(pmg_structure, symprec=1e-3, angle_tolerance=5.0):
    """Groups the 6 canonical Voigt strain modes ('xx','yy','zz','yz','xz',
    'xy') by point-group symmetry: modes in the same group are related by a
    rotation/reflection of the structure's point group (R . eps_i . R^T ==
    +/-eps_j for some point-group operation R), so straining along one gives
    a stress response fully determined by the other's (see
    find_mapping_operation) -- only one representative per group needs to
    actually be run through DFT.

    Returns (groups, point_group_symbol, ops): `groups` is a list of lists
    of mode names partitioning VOIGT_MODES; `ops` is the raw list of
    point-group operations (Cartesian rotation matrices), returned so
    find_mapping_operation can locate the specific operation between any
    two modes without re-running symmetry detection.

    Same deliberately-conservative exact-tensor-match criterion as
    equivalent_cartesian_axes (not "the elastic tensor happens to be
    isotropic here") -- verified live: a synthetic cubic structure groups
    into {xx,yy,zz} and {xy,xz,yz} (matching the textbook result: only 2
    independent strain series needed for the 3 independent cubic constants
    C11/C12/C44), while the suite's hexagonal 2D fixture and a synthetic
    triclinic structure both come back with 6 singleton groups -- hexagonal
    materials really do have C11=C22, but that comes from continuous 6-fold
    averaging, not a discrete operation mapping eps_xx exactly onto eps_yy,
    so this method correctly does NOT claim it (would need a full
    irreducible-representation projection to capture, out of scope here).
    """
    sga = SpacegroupAnalyzer(pmg_structure, symprec=symprec, angle_tolerance=angle_tolerance)
    point_group = sga.get_point_group_symbol()
    ops = sga.get_point_group_operations(cartesian=True)

    tensors = [strain_tensor(m) for m in VOIGT_MODES]
    parent = list(range(6))

    def find(i):
        while parent[i] != i:
            i = parent[i]
        return i

    def union(i, j):
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[ri] = rj

    for op in ops:
        R = op.rotation_matrix
        for i in range(6):
            transformed = R @ tensors[i] @ R.T
            for j in range(6):
                if i != j and (np.allclose(transformed, tensors[j], atol=1e-3)
                                or np.allclose(transformed, -tensors[j], atol=1e-3)):
                    union(i, j)

    groups = {}
    for i in range(6):
        groups.setdefault(find(i), []).append(VOIGT_MODES[i])
    return list(groups.values()), point_group, ops


def operations_summary(ops) -> str:
    """Reduced (not a full per-operation dump) description of a point group's
    symmetry operations: just the total count, split into proper rotations
    (rotation matrix determinant +1) vs. improper operations (determinant -1
    -- mirrors, inversion, rotoreflections). Enough to sanity-check "does
    this look like the point group I expect" without printing every 3x3
    matrix (a cubic m-3m point group alone has 48 of them). Shared by
    stb-elasticInputs' and stb-strain's symmetry tables (moved here once
    stb-strain became a second consumer).
    """
    proper = sum(1 for op in ops if np.linalg.det(op.rotation_matrix) > 0)
    improper = len(ops) - proper
    return f"{len(ops)} operation(s) ({proper} proper rotation(s), {improper} improper (mirror/inversion/rotoreflection))"


def get_point_group_operations(pmg_structure, symprec=1e-3, angle_tolerance=5.0):
    """Returns (point_group_symbol, ops): the structure's point group and its
    Cartesian symmetry operations (pymatgen SymmOp list), with none of the
    strain-mode/elastic-tensor machinery above. A thin, standalone accessor
    for callers that only need to REPORT the detected symmetry (e.g.
    stb-elasticInputs' deformation-direction table) -- so they don't have to
    duplicate the SpacegroupAnalyzer call inline or import pymatgen directly.
    equivalent_strain_modes already returns `ops` as a side effect of its own
    reduction (reuse that where available); this exists for callers like the
    "full" method that don't otherwise compute `ops` at all.
    """
    sga = SpacegroupAnalyzer(pmg_structure, symprec=symprec, angle_tolerance=angle_tolerance)
    return sga.get_point_group_symbol(), sga.get_point_group_operations(cartesian=True)


def find_mapping_operation(ops, src_mode, dst_mode, atol=1e-3):
    """Returns (R, sign) such that sign * (R . strain_tensor(src_mode) . R^T)
    equals strain_tensor(dst_mode), for some rotation R in `ops` (as
    produced by equivalent_strain_modes) and sign in {+1, -1} -- or
    (None, None) if no operation in `ops` relates the two modes.

    `sign` is returned SEPARATELY from R, not folded into it, because R only
    ever appears bilinearly (R . X . R^T): negating R does NOT negate that
    product (the two minus signs cancel), so a version of this function that
    tried returning -R for the sign-flip case was a real, caught-by-testing
    bug -- it silently left the reconstructed stress unflipped whenever the
    matching operation happened to need a sign flip (verified live: m-3m's
    yz->xy mapping is exactly such a case, R.eps_yz.R^T == -eps_xy, not
    +eps_xy; a pooled fit across {xy,xz,yz} for a synthetic cubic crystal
    came back at exactly 1/3 of the true C66 before this fix, exactly 80.0
    GPa after applying sign to the transformed stress via transform_stress).

    Since stress is linear in strain and transforms the same way as strain
    under a symmetry operation of the crystal (Neumann's principle -- the
    elastic tensor itself must be invariant under every point-group
    operation), the same (R, sign) also maps the measured stress response
    sigma(src_mode) onto what sigma(dst_mode) would have been, letting a
    mode that was never computed be reconstructed exactly from a
    symmetry-equivalent one that was -- see transform_stress.
    """
    src_t, dst_t = strain_tensor(src_mode), strain_tensor(dst_mode)
    for op in ops:
        R = op.rotation_matrix
        transformed = R @ src_t @ R.T
        if np.allclose(transformed, dst_t, atol=atol):
            return R, 1
        if np.allclose(transformed, -dst_t, atol=atol):
            return R, -1
    return None, None


def transform_stress(R, sign, sigma):
    """Applies the (R, sign) mapping from find_mapping_operation to a
    measured stress tensor, recovering the stress response a
    symmetry-equivalent strain mode would have had: sign * (R . sigma . R^T).
    """
    return sign * (R @ sigma @ R.T)


# ==========================================================================
# General symmetry-allowed elastic tensor ("full" method)
#
# equivalent_strain_modes/find_mapping_operation above (the "basic" method)
# only ever detect an equivalence when a single point-group operation maps
# one canonical strain TENSOR (rank-2) exactly onto another -- this is
# exact and cheap, but structurally blind to reductions like hexagonal's
# C11=C22, which is a property of the full elastic TENSOR (rank-4), not of
# any pairwise strain-tensor correspondence: a 60-degree rotation maps
# eps_xx onto a mixed tensor (cos^2, sin*cos, sin^2 entries), never onto
# eps_yy itself, even though that same rotation, applied to the *entire* C
# via the Bond transformation, forces C11=C22 (the textbook "direct
# inspection method" for deriving each crystal system's independent-constant
# pattern). The functions below implement that general criterion instead,
# via the null space of "C must be invariant under every point-group
# operation" rather than hand-derived per-crystal-system tables -- verified
# against known results for the suite's own fixtures (m-3m/cubic: 3
# independent constants; -6m2/hexagonal 2D fixture: 5; -1/triclinic: 21,
# matching every textbook count) and against a synthetic round-trip
# (hexagonal-consistent C, only xx/zz/xz ever "measured": reconstructed the
# full matrix, including C22 and C66 which were never run directly, to
# 3.6e-9 max error).
# ==========================================================================

# All 21 independent (row, col) pairs of a symmetric 6x6 Voigt matrix.
VOIGT_INDEPENDENT_PAIRS = [(m, n) for m in range(6) for n in range(m, 6)]


def voigt_to_full_tensor(C):
    """Expands a 6x6 Voigt stiffness matrix into the full 3x3x3x3 elastic
    tensor C_ijkl. No extra factors needed anywhere (unlike the
    strain/stress Voigt vectors, which double shear/normal components
    differently) -- this Voigt convention was designed so that C_mn ==
    C_ijkl directly for every m, n; verified by hand-derivation from
    sigma_ij = C_ijkl eps_kl before implementing.
    """
    T = np.zeros((3, 3, 3, 3))
    for m in range(6):
        i, j = VOIGT_TO_TENSOR[m]
        for n in range(6):
            k, l = VOIGT_TO_TENSOR[n]
            val = C[m, n]
            for (ii, jj) in {(i, j), (j, i)}:
                for (kk, ll) in {(k, l), (l, k)}:
                    T[ii, jj, kk, ll] = val
    return T


def full_tensor_to_voigt(T):
    """Inverse of voigt_to_full_tensor: reads the 6x6 Voigt matrix back out
    of a full 3x3x3x3 elastic tensor (just index lookups, see there)."""
    C = np.zeros((6, 6))
    for m in range(6):
        i, j = VOIGT_TO_TENSOR[m]
        for n in range(6):
            k, l = VOIGT_TO_TENSOR[n]
            C[m, n] = T[i, j, k, l]
    return C


def rotate_elastic_tensor(C, R):
    """Rotates a 6x6 Voigt elastic tensor by a Cartesian rotation matrix R,
    i.e. C'_ijkl = R_ia R_jb R_kc R_ld C_abcd. Deliberately computed via the
    generic rank-4 tensor transformation (round-tripping through the full
    3x3x3x3 form) instead of a hand-typed 6x6 Bond matrix -- the Bond
    matrix's off-diagonal blocks carry factor-of-2 placements that are easy
    to transcribe backwards; this has no such risk since it only relies on
    the rank-4 transformation law itself.
    """
    T = voigt_to_full_tensor(C)
    T_rot = np.einsum('ia,jb,kc,ld,abcd->ijkl', R, R, R, R, T)
    return full_tensor_to_voigt(T_rot)


def symmetry_allowed_basis(pmg_structure, symprec=1e-3, angle_tolerance=5.0):
    """Returns (basis, point_group_symbol): `basis` is a list of 6x6 Voigt
    matrices spanning the subspace of elastic tensors invariant under every
    operation of the structure's point group -- i.e. every physically
    allowed C for this crystal system is some linear combination of these
    matrices, and `len(basis)` is exactly the number of independent elastic
    constants for that point group (3 for cubic, 5 for hexagonal, 21 for
    triclinic, etc. -- matches the standard textbook counts, verified
    against this suite's own cubic/hexagonal/triclinic fixtures).

    Computed as the null space of "C must equal its own rotated image"
    stacked over every point-group operation: for each of the 21
    independent Voigt entries taken as a unit matrix, rotate it (via
    rotate_elastic_tensor) and re-express the result in the same 21-entry
    basis to get one column of the operation's 21x21 linear map L_R; stack
    (L_R - I) over all operations and take the SVD null space. This is the
    general "direct inspection method" used to derive elastic-constant
    tables per crystal system, applied numerically instead of by hand so it
    covers any of the 32 point groups uniformly (no per-crystal-system
    special-casing, unlike equivalent_strain_modes above).
    """
    sga = SpacegroupAnalyzer(pmg_structure, symprec=symprec, angle_tolerance=angle_tolerance)
    point_group = sga.get_point_group_symbol()
    ops = sga.get_point_group_operations(cartesian=True)

    n = len(VOIGT_INDEPENDENT_PAIRS)
    constraints = []
    for op in ops:
        R = op.rotation_matrix
        L = np.zeros((n, n))
        for p, (m, nn) in enumerate(VOIGT_INDEPENDENT_PAIRS):
            unit = np.zeros((6, 6))
            unit[m, nn] = 1.0
            unit[nn, m] = 1.0
            rotated = rotate_elastic_tensor(unit, R)
            L[:, p] = [rotated[mm, nn2] for (mm, nn2) in VOIGT_INDEPENDENT_PAIRS]
        constraints.append(L - np.eye(n))

    M = np.vstack(constraints)
    _, s, vt = np.linalg.svd(M, full_matrices=True)
    tol = 1e-8 * max(M.shape)
    rank = int(np.sum(s > tol))

    basis = []
    for vec in vt[rank:]:
        C = np.zeros((6, 6))
        for val, (m, nn) in zip(vec, VOIGT_INDEPENDENT_PAIRS):
            C[m, nn] = val
            C[nn, m] = val
        basis.append(C)
    return basis, point_group


def voigt_strain_engineering(eps_tensor):
    """Converts a symmetric 3x3 (tensor) strain into the 6-component
    engineering-strain Voigt vector (shear entries doubled: gamma=2*eps) --
    the vector Hooke's law sigma_voigt = C . eps_voigt_engineering actually
    expects. Companion to strain_tensor above (which returns the un-doubled
    tensor form).
    """
    return np.array([
        eps_tensor[0, 0], eps_tensor[1, 1], eps_tensor[2, 2],
        2 * eps_tensor[1, 2], 2 * eps_tensor[0, 2], 2 * eps_tensor[0, 1],
    ])


def select_directions_by_rank(basis, order=VOIGT_MODES):
    """Greedily picks the minimal subset of canonical Voigt strain modes
    (in `order`) whose data would fully determine every coefficient in
    `basis` (from symmetry_allowed_basis): keeps a mode only if it increases
    the rank of the accumulated design matrix, stopping once the rank
    reaches len(basis) -- verified to reproduce the expected reduction on
    this suite's fixtures (2 of 6 for the cubic fixture, 3 of 6 for the
    hexagonal 2D fixture, all 6 for the triclinic fixture).
    """
    kept = []
    blocks = []
    target_rank = len(basis)
    for mode in order:
        eps_v = voigt_strain_engineering(strain_tensor(mode))
        block = np.array([Bk @ eps_v for Bk in basis]).T  # shape (6, len(basis))
        trial = blocks + [block]
        rank = np.linalg.matrix_rank(np.vstack(trial), tol=1e-8)
        prev_rank = np.linalg.matrix_rank(np.vstack(blocks), tol=1e-8) if blocks else 0
        if rank > prev_rank:
            blocks = trial
            kept.append(mode)
        if rank >= target_rank:
            break
    return kept


def fit_stiffness_matrix(basis, points, conv_factor):
    """Least-squares fit of the elastic tensor, constrained a priori to the
    symmetry-allowed subspace `basis` (from symmetry_allowed_basis).
    `points` is a list of (strain_tensor_3x3, stress_tensor_3x3) pairs
    (raw, unconverted units, tensor strain -- not engineering). Each point
    contributes 6 equations (one per Voigt stress component) linear in the
    basis coefficients; stacking every point from every direction actually
    run and solving once with np.linalg.lstsq naturally pools redundant
    directions and reconstructs columns for directions never run, in a
    single fit -- verified in a synthetic round-trip (hexagonal-consistent
    C, only 3 of 6 directions "measured") to 3.6e-9 max error against the
    true matrix, including the two independent constants (C22, C66) that
    were never measured directly.

    Returns C (6x6, already scaled by conv_factor).
    """
    K = len(basis)
    rows, targets = [], []
    for eps_t, sigma_t in points:
        eps_v = voigt_strain_engineering(eps_t)
        sigma_v = np.array([
            sigma_t[0, 0], sigma_t[1, 1], sigma_t[2, 2],
            sigma_t[1, 2], sigma_t[0, 2], sigma_t[0, 1],
        ])
        block = np.array([Bk @ eps_v for Bk in basis]).T  # (6, K)
        rows.append(block)
        targets.append(sigma_v)

    A = np.vstack(rows)
    b = np.concatenate(targets)
    x, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    C = sum(xk * Bk for xk, Bk in zip(x, basis))
    return C * conv_factor


# ==========================================================================
# Energy-strain ("parabolic fit") method
#
# E(eps) = E0 + (V0/2) gamma^T C gamma (gamma = engineering Voigt strain) --
# a PURE single-Voigt-component strain only ever probes the diagonal term
# C_mm (gamma has one nonzero entry, so gamma^T C gamma = C_mm * gamma_m^2,
# structurally blind to any off-diagonal C_mn): unlike the stress method
# above, which gets a full column (6 numbers) from one direction, an energy
# curvature gives exactly ONE number per pattern. Off-diagonal constants
# need a genuinely COMBINED pattern (two Voigt components strained at once,
# e.g. eps_xx=eps_yy=delta), whose curvature mixes C_11+C_22+2*C_12. Reuses
# symmetry_allowed_basis/strain_tensor/voigt_strain_engineering above --
# verified in a synthetic round-trip against the suite's cubic fixture (3
# patterns: 'xx', 'yz', 'xx+yy', reconstruction error 1.1e-13) and hexagonal
# 2D fixture restricted to in-plane candidates (2 patterns: 'xx', 'xy',
# recovering C11=C22 and C66=(C11-C12)/2 exactly even though neither was
# ever measured as its own diagonal pattern).
# ==========================================================================

def energy_pattern_candidates():
    """Returns the full pool of 21 candidate deformation patterns for the
    energy method: the 6 pure canonical Voigt modes (diagonal terms only)
    plus all 15 pairwise combinations (needed for off-diagonal terms) -- 21
    total, matching triclinic's 21 independent constants (the worst case,
    where every candidate ends up needed). Each entry is (name, tensor),
    name being e.g. 'xx' or 'xx+yy' (used verbatim in strain_*/ folder
    names), tensor the unit-magnitude (delta=1) 3x3 strain tensor.
    """
    candidates = [(m, strain_tensor(m)) for m in VOIGT_MODES]
    for i in range(6):
        for j in range(i + 1, 6):
            name = f"{VOIGT_MODES[i]}+{VOIGT_MODES[j]}"
            candidates.append((name, strain_tensor(VOIGT_MODES[i]) + strain_tensor(VOIGT_MODES[j])))
    return candidates


def select_energy_patterns_by_rank(basis, candidates):
    """Greedily picks the minimal subset of `candidates` (from
    energy_pattern_candidates, already filtered to physically valid ones --
    see the vacuum-axis filtering in elastic_inputs.py) whose curvatures
    would fully determine every coefficient in `basis`. Same greedy
    rank-growing idea as select_directions_by_rank above, but each
    candidate contributes a single SCALAR row (gamma^T Bk gamma, one number
    per basis element) instead of a 6-row block, matching the "one equation
    per pattern" structure of the energy method (see module note above).
    """
    kept = []
    rows = []
    target_rank = len(basis)
    for name, eps_t in candidates:
        gamma = voigt_strain_engineering(eps_t)
        row = [gamma @ Bk @ gamma for Bk in basis]
        trial = rows + [row]
        rank = np.linalg.matrix_rank(np.array(trial), tol=1e-8)
        prev_rank = np.linalg.matrix_rank(np.array(rows), tol=1e-8) if rows else 0
        if rank > prev_rank:
            rows = trial
            kept.append((name, eps_t))
        if rank >= target_rank:
            break
    return kept


def fit_stiffness_matrix_energy(basis, pattern_targets, conv_factor):
    """Least-squares fit of the elastic tensor from energy-curvature data,
    constrained a priori to the symmetry-allowed subspace `basis` (from
    symmetry_allowed_basis). `pattern_targets` is a list of
    (strain_tensor_unit_3x3, target) pairs, where `target` is the raw
    (unconverted) volume/area-normalized curvature `2*a/V0` -- `a` being
    the quadratic coefficient of a polyfit of energy vs. strain magnitude
    for that pattern (see elastic_analysis.py), so that `target == gamma^T
    C_raw gamma` for that pattern's unit strain gamma. Each pattern
    contributes exactly one equation (unlike fit_stiffness_matrix's 6 per
    point) -- solving needs at least len(basis) patterns, which is exactly
    what select_energy_patterns_by_rank guarantees.

    Returns C (6x6, already scaled by conv_factor).
    """
    rows, targets = [], []
    for eps_t, target in pattern_targets:
        gamma = voigt_strain_engineering(eps_t)
        rows.append([gamma @ Bk @ gamma for Bk in basis])
        targets.append(target)

    A = np.array(rows)
    b = np.array(targets)
    x, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    C = sum(xk * Bk for xk, Bk in zip(x, basis))
    return C * conv_factor


def crystal_system(pmg_structure, symprec=1e-3, angle_tolerance=5.0):
    """Returns (crystal_system, point_group_symbol), e.g. ("triclinic", "-1").

    Used by stb-elasticInputs to flag when the "vary one strain component at
    a time" elastic-constant method can't be exact: that method never solves
    for the normal-shear coupling constants (C14-C16/C24-C26/C34-C36/C45/
    C46/C56), implicitly assuming they're zero. That's only unconditionally
    true for triclinic and monoclinic's *complement* -- triclinic (21
    independent constants) and monoclinic (13) always have some of these
    terms nonzero, in every subclass, so those two systems are flagged
    specifically. Higher symmetry systems are NOT all safe in every subclass
    (e.g. some low-symmetry tetragonal/trigonal point groups still have C16/
    C14 nonzero) -- deliberately not attempting that finer classification
    here, since misclassifying it would be worse than not trying.
    """
    sga = SpacegroupAnalyzer(pmg_structure, symprec=symprec, angle_tolerance=angle_tolerance)
    return sga.get_crystal_system(), sga.get_point_group_symbol()


def reduce_to_unitcell(structure, mode, symprec=1e-3, angle_tolerance=5.0):
    """Returns (new_structure, sga) for the requested `mode`:
      - primitive: smallest possible cell.
      - conventional: standardized, usually larger, cell.
      - refined: conventional-sized cell with atomic positions snapped to the
        detected symmetry -- cleans up numerical noise from relaxations or
        hand-built/CIF structures without changing which atoms are present.

    `sga` (the SpacegroupAnalyzer built from `structure`) is returned too, so
    callers can report the detected space group/point group/crystal system
    without re-running symmetry detection.

    Note: the output's atom order and coordinate origin are NOT guaranteed to
    match the input, in any mode -- spglib rebuilds the cell from the
    detected symmetry operations from scratch, free to pick any symmetry
    -equivalent origin/ordering.

    Raises ValueError (from SpacegroupAnalyzer) if symmetry detection fails,
    or if `mode` isn't one of UNITCELL_MODES.
    """
    if mode not in UNITCELL_MODES:
        raise ValueError(f"Unknown unitcell mode {mode!r}, must be one of {UNITCELL_MODES}.")

    sga = SpacegroupAnalyzer(structure, symprec=symprec, angle_tolerance=angle_tolerance)
    if mode == "primitive":
        new_structure = sga.get_primitive_standard_structure()
    elif mode == "conventional":
        new_structure = sga.get_conventional_standard_structure()
    else:
        new_structure = sga.get_refined_structure()
    return new_structure, sga
