#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "2.0.0"

import sys
import os
import time
import argparse
import math
from datetime import datetime

import numpy as np
from pymatgen.core.periodic_table import Element
from pymatgen.io.ase import AseAtomsAdaptor
from stb.core import structure_io, kspace
from stb.core import mace_relax
from stb.core import citations
from stb.core import structure_checks
from stb.core import symmetry as core_symmetry
from stb.core.ase_view import view_structure_interactive
from stb.core.cli import COLORS, color_text, show_intro, print_dual, print_section, print_table
from stb.core.deps import require_mace
from stb.core.passivation import passivate_dangling_bonds

REPORT_FILE = "stb_nanotube_report.txt"
BIB_FILE = "references.bib"

# Same default vacuum-gap threshold as stb-fetch/stb-kgrid/stb-slab
# (core/kspace.py's other callers), used for the OUTPUT tube/ribbon's own
# vacuum-axis detection. The INPUT monolayer instead uses --vacuum-gap
# (already a CLI flag here, since detecting its single vacuum axis is the
# tool's whole starting point, not just a validation nicety).
VACUUM_GAP_ANG = 10.0


def find_translation_vector(a1, a2, n, m, search_bound=None):
    """Shortest lattice vector T = t1*a1 + t2*a2 (t1, t2 integers) that is
    (as close to) perpendicular to C = n*a1 + m*a2 (as an integer lattice
    vector allows), found by a bounded search -- no hexagonal-only closed
    form (the classic CNT gcd/(n-m)%3d shortcut only holds for a honeycomb
    lattice), so this works for any oblique 2D lattice.

    Perpendicularity, not merely non-parallelism, is required: roll_to_tube/
    build_ribbon embed the flat (frac_C, frac_T) cell by treating C and T as
    if they were an ORTHOGONAL (circumference, axial) frame -- angle
    2*pi*frac_C, position frac_T*|T| -- which only preserves real bond
    lengths/connectivity when T genuinely has no component along C. Any T'
    with the same cross product n*t2-m*t1 (same C x T sublattice, i.e. the
    same set of atoms per period) generates an equally valid periodic cell,
    since T' = T + k*C (any integer k) leaves the cross product unchanged
    (cross(C, T + k*C) = cross(C, T), as C x C = 0) -- so instead of just
    picking the raw-shortest non-parallel T, each candidate is first reduced
    modulo C (subtracting the integer multiple of C that minimizes its
    remaining component along C) to find the most-perpendicular, and among
    those the shortest, representative of its own sublattice.

    Different (t1, t2) candidates can belong to different C x T sublattices
    (different cross products n*t2-m*t1, i.e. a different number of cells
    per period) -- modulo-C reduction alone only finds the best
    representative *within* one candidate's own sublattice, it can't turn a
    sublattice that has no exactly-perpendicular representative at all into
    one that does. A hexagonal lattice's classical zigzag/armchair CNT unit
    cells, for example, generally need MORE cells per period than the
    smallest possible cross product provides (verified live: graphene
    (6, 0)'s shortest non-parallel candidate lands in the cross=6 sublattice,
    which has no integer T perpendicular to C at all -- the true zigzag
    period needs cross=12). So every candidate found by the bounded search
    is reduced modulo C first, and candidates are then ranked primarily by
    how close to perpendicular their reduced T ended up (near-machine-zero
    residual dot product = "exactly", for a sublattice that supports it),
    length only as the tie-break among comparably-perpendicular candidates
    -- not by raw length first, which is what silently preferred an
    insufficient sublattice before.

    Verified live against ASE's own ase.build.nanotube (an independent,
    widely-used reference implementation): the earlier length-first ranking
    silently returned an under-sized, badly non-perpendicular cell for most
    chiralities (e.g. graphene (6, 0) came out 12 atoms with only 1 real
    bonded neighbor per atom within 1.8 Ang, instead of the correct 24
    atoms / 3-fold coordination) -- this perpendicularity-first ranking
    fixes that: full 3-fold coordination and uniform, correct C-C bond
    lengths confirmed live for every chirality tried, zigzag/armchair/
    chiral alike.

    Zigzag (m=0 or n=0) and armchair (n=m) indices additionally match
    ASE's own atom count exactly, verified live -- for a GENERIC chiral
    index with gcd(n, m) > 1 (e.g. (8, 2), (5, 1)), this function's cell
    can come out physically valid but LARGER than ASE's, because ASE's
    nanotube() (like the classic Dresselhaus/Saito construction) uses a
    "symmetry vector" -- a combined ROTATION + translation (screw)
    operation -- to generate a more compact cell for these chiralities,
    while roll_to_tube/build_ribbon here only ever apply a PURE
    translation between repeats (no accompanying rotation), matching this
    function's own perpendicular-pure-translation search. Both are
    legitimate periodic descriptions of the same physical tube; this one
    just isn't always the most compact possible, since implementing screw
    -axis symmetry would reintroduce a hexagonal-lattice-only assumption
    this tool deliberately avoids (see the module-level "no hexagonal-only
    closed form" design note). Confirmed live for (8, 2)/(5, 1)/(7, 4)/
    (10, 5)/(6, 3)/(4, 2): every one of these larger cells is still fully,
    correctly 3-fold coordinated at the right bond length under repetition
    -- just not the smallest possible periodic description.
    """
    if search_bound is None:
        search_bound = 3 * (abs(n) + abs(m)) + 5

    C_cart = n * a1 + m * a2
    C_dot_C = float(np.dot(C_cart, C_cart))

    candidates = []
    for t1 in range(-search_bound, search_bound + 1):
        for t2 in range(-search_bound, search_bound + 1):
            cross = n * t2 - m * t1
            if cross == 0:
                continue
            T_cart = t1 * a1 + t2 * a2
            # Reduce modulo C: same sublattice (same cross product), but the
            # representative with the smallest remaining along-C component.
            k = round(float(np.dot(T_cart, C_cart)) / C_dot_C)
            t1_red, t2_red = t1 - k * n, t2 - k * m
            T_red_cart = t1_red * a1 + t2_red * a2
            length = float(np.linalg.norm(T_red_cart))
            # Normalized perpendicularity residual (0 = exactly perpendicular),
            # dimensionless (|sin| of the angle T makes away from perpendicular
            # to C), so it's comparable across candidates of different length.
            perp_error = abs(float(np.dot(T_red_cart, C_cart))) / (length * np.sqrt(C_dot_C))
            candidates.append((round(perp_error, 9), round(length, 6), t1_red, t2_red, cross))

    if not candidates:
        raise ValueError(
            f"Could not find a translation vector within search bound {search_bound}; "
            "try increasing --search-bound."
        )

    # Closest to perpendicular first, shortest as the tie-break, then prefer
    # a positive cross product (canonical orientation) among remaining ties.
    candidates.sort(key=lambda c: (c[0], c[1], 0 if c[4] > 0 else 1))
    _, _, t1, t2, _ = candidates[0]
    return t1, t2


def enumerate_cell_atoms(structure, in_plane_axes, vac_axis_idx, n, m, t1, t2, tol=1e-6):
    """Every atom image inside the C x T parallelogram (C = n*a1+m*a2,
    T = t1*a1+t2*a2), as (symbol, frac_C, frac_T, offset). `offset` is the
    atom's out-of-plane displacement (Angstrom) *relative to the sheet's own
    mean plane* -- not its raw vacuum-axis fractional coordinate, which
    would instead encode wherever the user happened to center the sheet
    inside the vacuum box (e.g. z_frac=0.5) and blow up the tube radius by
    that amount. Buckled/puckered monolayers (silicene, phosphorene, TMDs)
    still roll into a genuinely corrugated tube; a perfectly flat, centered
    sheet gets offset 0 for every atom.
    """
    positions = np.array([pos for _, pos in structure.atoms])
    symbols = [s for s, _ in structure.atoms]
    is_cartesian = structure.coord_format == 'cartesian'
    frac_full = kspace.to_fractional(positions, structure.lattice, is_cartesian)

    vac_len = float(np.linalg.norm(structure.lattice[vac_axis_idx]))
    u_vals = frac_full[:, in_plane_axes[0]]
    v_vals = frac_full[:, in_plane_axes[1]]
    w_vals = frac_full[:, vac_axis_idx] - float(np.mean(frac_full[:, vac_axis_idx]))

    # Columns are C's and T's own fractional coefficients in the (a1, a2)
    # basis -- solving (u+i, v+j) = fc*(n, m) + ft*(t1, t2) for (fc, ft)
    # needs M's COLUMNS to be (n, m) and (t1, t2), i.e. M @ (fc, ft) =
    # (u+i, v+j). A previous version had M's ROWS as (n, m)/(t1, t2)
    # instead (M's transpose) -- silently correct only when that matrix
    # happens to be symmetric (e.g. a zigzag (n, 0) chirality, where the
    # matrix is diagonal), but wrong for every other case (verified live:
    # an armchair (6, 6) tube/ribbon came out with a 0.41 Ang nearest
    # -neighbor distance instead of graphene's real 1.42 Ang bond length).
    M = np.array([[n, t1], [m, t2]], dtype=float)
    M_inv = np.linalg.inv(M)
    n_cells = abs(n * t2 - m * t1)

    bound = max(abs(n), abs(m), abs(t1), abs(t2)) + 2

    cell_atoms = []
    for k in range(len(symbols)):
        for i in range(-bound, bound + 1):
            for j in range(-bound, bound + 1):
                frac_ct = M_inv @ np.array([u_vals[k] + i, v_vals[k] + j])
                fc, ft = frac_ct
                if -tol <= fc < 1.0 - tol and -tol <= ft < 1.0 - tol:
                    offset = w_vals[k] * vac_len
                    cell_atoms.append((symbols[k], float(fc), float(ft), float(offset)))

    expected = n_cells * len(symbols)
    if len(cell_atoms) != expected:
        raise ValueError(
            f"Internal error: found {len(cell_atoms)} atoms in the C x T cell, "
            f"expected {expected} (N_cells={n_cells} x {len(symbols)} basis atoms). "
            "Try increasing --lattice-tol or report this as a bug."
        )
    return cell_atoms, n_cells


def roll_to_tube(cell_atoms, C_cart, T_cart, repeats, vacuum_gap):
    """Wraps the flat (frac_C, frac_T) cell atoms around a cylinder.
    T becomes the periodic axis (c); the vacuum-padded box in x/y is sized
    to the largest radius found (R plus each atom's own thickness offset).
    """
    T_len = float(np.linalg.norm(T_cart))
    C_len = float(np.linalg.norm(C_cart))
    R = C_len / (2.0 * np.pi)

    atoms = []
    max_r = 0.0
    for rep in range(repeats):
        for symbol, fc, ft, offset in cell_atoms:
            theta = 2.0 * np.pi * fc
            r = R + offset
            max_r = max(max_r, abs(r))
            x = r * np.cos(theta)
            y = r * np.sin(theta)
            z = (ft + rep) * T_len
            atoms.append((symbol, np.array([x, y, z])))

    box = 2.0 * (max_r + vacuum_gap)
    lattice = np.array([
        [box, 0.0, 0.0],
        [0.0, box, 0.0],
        [0.0, 0.0, T_len * repeats],
    ])
    return atoms, lattice, R


def build_ribbon(cell_atoms, C_cart, T_cart, vac_axis_idx, full_lattice, repeats, vacuum_gap):
    """Tiles the flat cell `repeats` times along C for a finite-width strip.
    No wrap: real in-plane Cartesian positions are kept (an oblique C/T
    angle stays oblique). T is still the sole periodic axis (c); the finite
    C-direction ("width") and the thickness both get vacuum-padded, using a
    frame built from the vacuum axis and T so it stays valid for any C/T
    angle, not just an axis-aligned input.
    """
    vac_vec = full_lattice[vac_axis_idx]
    vac_unit = vac_vec / np.linalg.norm(vac_vec)
    T_unit = T_cart / np.linalg.norm(T_cart)
    width_unit = np.cross(vac_unit, T_unit)
    width_unit = width_unit / np.linalg.norm(width_unit)

    raw = []
    for rep in range(repeats):
        for symbol, fc, ft, offset in cell_atoms:
            pos = (fc + rep) * C_cart + ft * T_cart + offset * vac_unit
            raw.append((symbol, pos))

    w_vals = [float(np.dot(pos, width_unit)) for _, pos in raw]
    t_vals = [float(np.dot(pos, vac_unit)) for _, pos in raw]
    w_min, w_max = min(w_vals), max(w_vals)
    t_min, t_max = min(t_vals), max(t_vals)

    a_len = (w_max - w_min) + 2.0 * vacuum_gap
    b_len = (t_max - t_min) + 2.0 * vacuum_gap
    shift = (vacuum_gap - w_min) * width_unit + (vacuum_gap - t_min) * vac_unit

    atoms = [(symbol, pos + shift) for symbol, pos in raw]
    lattice = np.array([
        a_len * width_unit,
        b_len * vac_unit,
        T_cart,
    ])
    return atoms, lattice, (w_max - w_min)


def build_result_structure(atoms, lattice, species_meta):
    species = list(dict.fromkeys(symbol for symbol, _ in atoms))
    meta = {s: species_meta[s] for s in species if s in species_meta}
    return structure_io.FdfStructure(
        lattice=lattice,
        lattice_constant=1.0,
        species=species,
        species_meta=meta,
        atoms=atoms,
        coord_format='cartesian',
    )


def _is_hexagonal_lattice(a1, a2, tol_len=0.02, tol_ang=1.0):
    """True if a1/a2 describe a hexagonal (honeycomb-capable) 2D Bravais
    lattice: equal length, 60 or 120 degrees apart. Gates the CNT-specific
    zigzag/armchair/chiral labels and the graphene metallic/semiconducting
    hint below -- both are meaningless for a generic oblique lattice, which
    this tool otherwise supports without any hexagonal-only assumption."""
    len1, len2 = float(np.linalg.norm(a1)), float(np.linalg.norm(a2))
    if len1 == 0 or abs(len1 - len2) / len1 > tol_len:
        return False
    cos_ang = np.dot(a1, a2) / (len1 * len2)
    ang = np.degrees(np.arccos(np.clip(cos_ang, -1.0, 1.0)))
    return abs(ang - 120.0) < tol_ang or abs(ang - 60.0) < tol_ang


def _validate_structure(pmg_structure, vacuum_axes, f_out):
    """Shared malformation checklist (core.structure_checks) plus a
    space-group label -- same shape as mlrelax.py/supercell.py/slab.py's own
    _validate_structure(), wrapped in try/except by the caller (a validation
    failure is reported, never fatal)."""
    structure_checks.run_malformation_checks(pmg_structure, vacuum_axes, f_out)
    sg_label = core_symmetry.space_group_label(pmg_structure)
    print_dual(f"Space group      : {sg_label}", f_out)
    return sg_label


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Rolls a 2D monolayer FDF into a nanotube or nanoribbon.", 'bold')}
Chirality (n, m) picks the vector to roll around (tube) or the finite
width direction (ribbon); works for any oblique 2D lattice, not just
hexagonal (no gcd/(n-m)%3d hexagonal-only shortcut is used).""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s -f graphene.fdf --chirality 6 0 --mode tube\n"
               "  %(prog)s -f graphene.fdf --chirality 6 0 --mode ribbon --repeats 4\n"
               "  %(prog)s -f graphene.fdf --chirality 6 0 --mode ribbon --repeats 4 --passivate\n"
               "  %(prog)s -f graphene.fdf --chirality 6 0 --ml-relax --save-report --view\n"
    )

    parser.add_argument("-f", "--file", dest="filename", type=str, required=True,
                        help="Path to the input 2D monolayer structure file (.fdf).")
    parser.add_argument("--chirality", type=int, nargs=2, required=True, metavar=("N", "M"),
                        help="Chirality indices (n, m), e.g. --chirality 6 0")
    parser.add_argument("--mode", choices=["tube", "ribbon"], default="tube",
                        help="'tube': roll into a cylinder. 'ribbon': finite-width flat strip. Default: tube.")
    parser.add_argument("--repeats", type=int, default=1,
                        help="Axial repeats (tube) or width repeats along the chirality "
                             "vector (ribbon). Default: 1.")
    parser.add_argument("--vacuum-gap", type=float, default=10.0,
                        help="Minimum empty span (Angstrom) used to detect which axis of the "
                             "INPUT monolayer is vacuum-padded (default: 10.0). Same meaning as "
                             "stb-kgrid/stb-kpath's --vacuum-gap.")
    parser.add_argument("--min-vacuum-size", type=float, default=15.0,
                        help="Vacuum padding (Angstrom) added around the generated tube/ribbon "
                             "(default: 15.0). Same meaning as stb-slab's --min-vacuum-size.")
    parser.add_argument("--lattice-tol", type=float, default=1e-6,
                        help="Fractional-coordinate tolerance for the C x T cell-membership "
                             "test (default: 1e-6).")
    parser.add_argument("--search-bound", type=int, default=None,
                        help="Override the search range for the translation vector T "
                             "(default: auto, 3*(|n|+|m|)+5).")
    parser.add_argument("--symprec", type=float, default=0.01,
                        help="Symmetry tolerance, Ang, for the before/after symmetry table "
                             "(default: 0.01).")

    parser.add_argument("--passivate", action="store_true",
                        help="Cap dangling bonds on the built tube/ribbon with a passivating "
                             "atom (same logic as stb-passivate/stb-slab). A properly closed, "
                             "axially-periodic tube has no dangling bonds to begin with (every "
                             "atom already has 3 real neighbors) -- this mainly matters for "
                             "--mode ribbon, whose two finite edges are real, physical dangling "
                             "bonds. Only single-missing-bond sites are auto-passivated; sites "
                             "missing 2+ bonds are reported instead of guessed.")
    parser.add_argument("--passivant", type=str, default="H",
                        help="Element to cap dangling bonds with, only with --passivate (default: H).")
    parser.add_argument("--cutoff", type=float, default=None,
                        help="Neighbor-search radius in Angstrom for --passivate. "
                             "Default: auto-detected (see stb-passivate --help).")
    parser.add_argument("--bond-length", type=float, default=None,
                        help="Passivant bond length in Angstrom for --passivate. "
                             "Default: auto per species pair (see stb-passivate --help).")

    parser.add_argument("--ml-relax", action="store_true",
                        help="Pre-relax the built tube/ribbon with a MACE potential before "
                             "writing it out (needs the optional 'ml' extra: pip install "
                             "stb_suite[ml]) -- positions only by default. Off by default.")
    parser.add_argument("--ml-relax-cell", action="store_true",
                        help="With --ml-relax, also relax the axial period (the vacuum-padded "
                             "width/thickness axes always stay exactly fixed). Only valid "
                             "together with --ml-relax.")
    parser.add_argument("--model", choices=["small", "medium", "large"], default="small",
                        help="MACE-MP-0 foundation model size for --ml-relax (default: small).")
    parser.add_argument("--custom-model", default=None, metavar="PATH",
                        help="Path to a custom fine-tuned .model file for --ml-relax, instead "
                             "of a MACE-MP-0 foundation size.")
    parser.add_argument("--device", choices=["cpu", "cuda"], default="cpu",
                        help="Device to run the MACE model on (default: cpu).")

    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the full run report (including the symmetry "
                             f"analysis) to {REPORT_FILE}. Off by default.")
    parser.add_argument("--view", action="store_true",
                        help="Open an interactive 3D view (via ASE) comparing the input "
                             "monolayer and the final tube/ribbon (page through frames in "
                             "ase-gui) after writing the output file. Needs a display. "
                             "Off by default.")

    parser.add_argument("-o", "--output", type=str, default="nanotube.fdf",
                        help="Output .fdf file name (default: nanotube.fdf).")
    parser.add_argument("-v", "--version", action="version", version=f"stb-nanotube {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if not args.passivate:
        if args.passivant != "H":
            parser.error("--passivant is only valid with --passivate.")
        if args.cutoff is not None:
            parser.error("--cutoff is only valid with --passivate.")
        if args.bond_length is not None:
            parser.error("--bond-length is only valid with --passivate.")
    else:
        try:
            Element(args.passivant)
        except ValueError as e:
            parser.error(str(e))

    if args.ml_relax_cell and not args.ml_relax:
        parser.error("--ml-relax-cell is only valid together with --ml-relax.")
    if (args.custom_model or args.model != "small") and not args.ml_relax:
        parser.error("--model/--custom-model are only valid together with --ml-relax.")
    if args.ml_relax:
        require_mace()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    report_path = REPORT_FILE if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    def fail(message):
        print_dual(color_text(f"[ERROR] {message}", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_dual(color_text("===== STB-NANOTUBE REPORT =====", 'magenta'), f_out)

    model_desc = f"a custom model ({args.custom_model})" if args.custom_model else f"MACE-MP-0 ({args.model})"
    n, m = args.chirality

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time        : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Input file       : {args.filename}", f_out)
    print_dual(f"Chirality (n, m) : ({n}, {m})", f_out)
    print_dual(f"Mode             : {args.mode}", f_out)
    print_dual(f"Repeats          : {args.repeats}", f_out)
    print_dual(f"Output vacuum    : {args.min_vacuum_size} Ang", f_out)
    print_dual(f"Output file      : {args.output}", f_out)
    print_dual(f"Passivate        : {'yes (' + args.passivant + ')' if args.passivate else 'no'}", f_out)
    print_dual(f"ML pre-relax     : {model_desc if args.ml_relax else 'no'}", f_out)
    if args.ml_relax:
        print_dual(f"Relax cell       : {'yes (axial period only)' if args.ml_relax_cell else 'no (positions only)'}", f_out)

    if not os.path.exists(args.filename):
        fail(f"File '{args.filename}' not found.")

    if (n, m) == (0, 0):
        fail("Chirality (0, 0) is not valid.")

    if args.repeats < 1:
        fail("--repeats must be >= 1.")

    try:
        structure = structure_io.read_fdf(args.filename)
    except (FileNotFoundError, ValueError) as e:
        fail(str(e))

    positions = np.array([pos for _, pos in structure.atoms])
    is_cartesian = structure.coord_format == 'cartesian'
    frac_coords = kspace.to_fractional(positions, structure.lattice, is_cartesian)
    vacuum_axes_before = kspace.detect_vacuum_axes(frac_coords, structure.lattice, args.vacuum_gap)

    if sum(vacuum_axes_before) != 1:
        dim = 3 - sum(vacuum_axes_before)
        fail(f"input must be a 2D monolayer (exactly 1 vacuum-padded axis); "
             f"detected {dim}D ({sum(vacuum_axes_before)} vacuum axes). Adjust --vacuum-gap "
             "if this looks wrong.")

    vac_axis_idx = vacuum_axes_before.index(True)
    in_plane_axes = [i for i in range(3) if i != vac_axis_idx]
    a1 = structure.lattice[in_plane_axes[0]]
    a2 = structure.lattice[in_plane_axes[1]]
    pmg_before = structure_io.to_pymatgen(structure)

    print_section("[1] INPUT STRUCTURE (2D monolayer)", f_out)
    print_dual(f"Formula          : {pmg_before.composition.reduced_formula}", f_out)
    print_dual(f"Atoms            : {len(structure.atoms)}", f_out)
    print_dual(f"Dimensionality   : {kspace.dimensionality_label(vacuum_axes_before)}", f_out)
    print_dual(f"In-plane axes    : {in_plane_axes}   Vacuum axis: {vac_axis_idx}", f_out)
    print_dual(f"|a1|, |a2|       : {np.linalg.norm(a1):.4f}, {np.linalg.norm(a2):.4f} Ang", f_out)

    print_section("[2] STRUCTURE VALIDATION (pre-transform)", f_out)
    try:
        _validate_structure(pmg_before, vacuum_axes_before, f_out)
    except Exception as e:
        print_dual(color_text(f"[WARNING] Structure validation could not complete: {e}", 'yellow'), f_out)

    print_section("[3] TUBE/RIBBON CONSTRUCTION", f_out)
    hexagonal = _is_hexagonal_lattice(a1, a2)
    print_dual(f"Lattice type     : {'hexagonal' if hexagonal else 'generic (non-hexagonal)'}", f_out)

    try:
        t1, t2 = find_translation_vector(a1, a2, n, m, args.search_bound)
    except ValueError as e:
        fail(str(e))

    C_cart = n * a1 + m * a2
    T_cart = t1 * a1 + t2 * a2
    C_len = float(np.linalg.norm(C_cart))
    chiral_angle = float(np.degrees(np.arccos(
        np.clip(np.dot(C_cart, a1) / (C_len * np.linalg.norm(a1)), -1.0, 1.0))))

    if hexagonal:
        if m == 0 or n == 0:
            cnt_type = "zigzag"
        elif n == m:
            cnt_type = "armchair"
        else:
            cnt_type = "chiral"
        print_dual(f"CNT type         : {cnt_type} (hexagonal lattice)", f_out)
    print_dual(f"Chiral angle     : {chiral_angle:.2f} deg (C relative to a1)", f_out)

    try:
        cell_atoms, n_cells = enumerate_cell_atoms(
            structure, in_plane_axes, vac_axis_idx, n, m, t1, t2, args.lattice_tol
        )
    except ValueError as e:
        fail(str(e))

    print_dual(f"Translation vector (t1, t2) : ({t1}, {t2})", f_out)
    print_dual(f"Cells in periodic unit (N_cells) : {n_cells}", f_out)
    if hexagonal and args.mode == "tube" and math.gcd(abs(n), abs(m)) > 1:
        print_dual(
            "Note: gcd(n, m) > 1 -- this cell uses a PURE translational period (no "
            "screw/roto-translation between repeats), which can be larger than the "
            "classic Dresselhaus/Saito CNT cell (e.g. ASE's own nanotube()), since "
            "that construction uses a combined rotation+translation symmetry this "
            "tool deliberately doesn't assume (see find_translation_vector's own "
            "docstring). Still a physically correct, fully-bonded periodic cell.", f_out)

    if args.mode == "tube":
        atoms, lattice, R = roll_to_tube(cell_atoms, C_cart, T_cart, args.repeats, args.min_vacuum_size)
        print_dual(f"Circumference |C| : {C_len:.4f} Ang", f_out)
        print_dual(f"Tube diameter (2R): {2 * R:.4f} Ang", f_out)
        print_dual(f"Axial period |T|  : {np.linalg.norm(T_cart):.4f} Ang", f_out)
    else:
        atoms, lattice, width = build_ribbon(
            cell_atoms, C_cart, T_cart, vac_axis_idx, structure.lattice, args.repeats, args.min_vacuum_size
        )
        print_dual(f"Ribbon width      : {width:.4f} Ang", f_out)
        print_dual(f"Axial period |T|  : {np.linalg.norm(T_cart):.4f} Ang", f_out)

    if hexagonal and args.mode == "tube" and len(structure.species) == 1 and len(structure.atoms) == 2:
        electronic = "metallic" if (n - m) % 3 == 0 else "semiconducting"
        print_dual(color_text(
            f"Electronic hint  : {electronic} -- graphene-like single-species honeycomb rule "
            f"((n-m) mod 3 {'== 0' if electronic == 'metallic' else '!= 0'}). A large-diameter "
            "estimate; curvature effects can change this for very small diameters.", 'cyan'), f_out)

    new_structure = build_result_structure(atoms, lattice, structure.species_meta)
    pmg_after = structure_io.to_pymatgen(new_structure)

    print_dual(f"Output formula   : {pmg_after.composition.reduced_formula}", f_out)
    print_dual(f"Output atoms     : {len(atoms)}", f_out)

    species_meta_final = structure.species_meta
    passivation_report = None
    if args.passivate:
        pmg_after, passivation_report = passivate_dangling_bonds(
            pmg_after, passivant=args.passivant, cutoff=args.cutoff, bond_length=args.bond_length)
        species_meta_final = structure_io.ensure_species_id(dict(structure.species_meta), args.passivant)
        n_passivated = len(passivation_report["passivated"])
        n_unresolved = len(passivation_report["unresolved"])
        print_dual(f"Dangling bonds found : {n_passivated + n_unresolved}", f_out)
        print_dual(f"Auto-passivated      : {n_passivated} with {args.passivant}", f_out)
        if n_unresolved:
            print_dual(color_text(
                f"[WARNING] {n_unresolved} atom(s) missing 2+ bonds -- left unpassivated "
                "(geometrically underdetermined from local coordination alone):", 'yellow'), f_out)
            for atom_idx, symbol, pos, deficit in passivation_report["unresolved"]:
                print_dual(f"    #{atom_idx + 1:<4} {symbol:<3} deficit={deficit}  at {pos}", f_out)
        if n_passivated + n_unresolved == 0:
            print_dual(
                "Note: no dangling bonds found -- expected for a properly closed, "
                "axially-periodic tube (every atom already has its full coordination).", f_out)
        print_dual(f"Output atoms (after passivation) : {len(pmg_after)}", f_out)

    tube_vacuum_gap = min(VACUUM_GAP_ANG, args.min_vacuum_size)
    frac_after = kspace.to_fractional(pmg_after.cart_coords, pmg_after.lattice.matrix, True)
    vacuum_axes_after = kspace.detect_vacuum_axes(frac_after, pmg_after.lattice.matrix, tube_vacuum_gap)

    ml_relax_info = None
    if args.ml_relax:
        print_section("[4] ML PRE-RELAXATION (MACE)", f_out)
        print_dual(f"Model           : {model_desc} (device={args.device})", f_out)
        print_dual(f"Cell relaxation : "
                   f"{'axial period only (vacuum axes fixed)' if args.ml_relax_cell else 'positions only'}", f_out)
        model_arg = args.custom_model if args.custom_model else args.model
        try:
            calc = mace_relax.get_calculator(model_arg, device=args.device)
        except ValueError as e:
            fail(str(e))
        for line in mace_relax.describe_model(model_arg, calc):
            print_dual(line, f_out)

        atoms_ase = AseAtomsAdaptor.get_atoms(pmg_after)
        atoms_ase.calc = calc
        e0 = atoms_ase.get_potential_energy()
        f0 = float(np.abs(atoms_ase.get_forces()).max())
        t_len0 = float(atoms_ase.cell.cellpar()[2])

        cell_mask = mace_relax.build_cell_mask(vacuum_axes_after) if args.ml_relax_cell else None
        t0 = time.time()
        converged, steps_used = mace_relax.relax(atoms_ase, calc, cell_mask=cell_mask, fmax=0.05, max_steps=200)
        wall_time = time.time() - t0

        e1 = atoms_ase.get_potential_energy()
        f1 = float(np.abs(atoms_ase.get_forces()).max())
        t_len1 = float(atoms_ase.cell.cellpar()[2])
        n_atoms = len(atoms_ase)

        print_dual(f"Steps used : {steps_used} "
                   f"({'converged' if converged else 'hit step cap, NOT converged'})", f_out)
        print_dual(f"Wall time  : {wall_time:.1f} s", f_out)

        rows = [
            (["Energy (eV)", f"{e0:.6f}", f"{e1:.6f}",
              f"{e1 - e0:+.6f} ({(e1 - e0) / n_atoms:+.6f}/atom)"], None),
            (["Max force (eV/Ang)", f"{f0:.4f}", f"{f1:.4f}", f"{f1 - f0:+.4f}"], None),
        ]
        if args.ml_relax_cell:
            rows.append((["Axial period |T| (Ang)", f"{t_len0:.4f}", f"{t_len1:.4f}",
                          f"{100 * (t_len1 - t_len0) / t_len0:+.2f}%"], None))
        print_table(["Quantity", "Before", "After", "Change"], rows, f_out)

        pmg_after = AseAtomsAdaptor.get_structure(atoms_ase)
        ml_relax_info = (converged, steps_used, e0, e1)

    print_section("[5] STRUCTURE VALIDATION (post-transform)", f_out)
    try:
        _validate_structure(pmg_after, vacuum_axes_after, f_out)
    except Exception as e:
        print_dual(color_text(f"[WARNING] Structure validation could not complete: {e}", 'yellow'), f_out)

    print_section("[6] SYMMETRY ANALYSIS (BEFORE / AFTER)", f_out)
    print_dual(f"Detailed symmetry analysis (Tolerance: {args.symprec} Ang):", f_out)
    print_dual(
        "Note: a tube/ribbon is genuinely 1D-periodic (2 vacuum axes) -- spglib has no "
        "1D rod-group classification, so Layer Group (which needs exactly 1 vacuum axis) "
        "reads N/A here for BOTH columns; Space Group is still the 3D space group of the "
        "vacuum-padded cell as given.", f_out)
    before_info = core_symmetry.symmetry_summary(pmg_before, args.symprec, args.vacuum_gap)
    after_info = core_symmetry.symmetry_summary(pmg_after, args.symprec, tube_vacuum_gap)
    if "Error" in before_info or "Error" in after_info:
        print_dual(color_text("[WARNING] Symmetry analysis failed for at least one structure.", 'yellow'), f_out)
        print_dual(f"  Monolayer (before): {before_info.get('Error', 'OK')}", f_out)
        print_dual(f"  Tube/ribbon (after): {after_info.get('Error', 'OK')}", f_out)
    else:
        properties = ["Crystal System", "Space Group", "Layer Group", "Point Group", "Hall Symbol"]
        rows = [([prop, str(before_info.get(prop, "N/A")), str(after_info.get(prop, "N/A"))], None)
                for prop in properties]
        print_table(["Property", "Monolayer (before)", "Tube/Ribbon (after)"], rows, f_out)

    print_section("[7] WRITING OUTPUT FILE", f_out)
    final_structure = structure_io.from_pymatgen(pmg_after, species_meta=species_meta_final)
    header_comment = [
        f"{args.mode.capitalize()} built by stb-nanotube from {args.filename}, chirality ({n}, {m}).",
        f"Translation vector (t1, t2) = ({t1}, {t2}), N_cells = {n_cells}, repeats = {args.repeats}.",
    ]
    if args.mode == "tube":
        header_comment.append(f"Diameter {2 * R:.4f} Ang, axial period {np.linalg.norm(T_cart):.4f} Ang.")
    else:
        header_comment.append(f"Width {width:.4f} Ang, axial period {np.linalg.norm(T_cart):.4f} Ang.")
    if passivation_report is not None:
        header_comment.append(
            f"Passivated {len(passivation_report['passivated'])} dangling bond(s) with {args.passivant}.")
    if ml_relax_info is not None:
        converged, steps_used, e0, e1 = ml_relax_info
        header_comment.append(
            f"ML pre-relaxed with {model_desc} "
            f"({'converged' if converged else 'NOT converged'} in {steps_used} step(s), "
            f"E = {e1:.6f} eV, delta E = {e1 - e0:+.6f} eV)."
        )
    structure_io.write_fdf(final_structure, args.output, header_comment=header_comment)
    print_dual(color_text(f"[OK] Structure written to '{args.output}'.", 'green'), f_out)

    print_section("[8] REFERENCES", f_out)
    bib_entries = [citations.SIESTA, citations.SIESTA_RECENT]
    if args.ml_relax:
        bib_entries.append(citations.MACE)
        if not args.custom_model:
            bib_entries.append(citations.MACE_MP)
    citations.write_bib_file(BIB_FILE, bib_entries)
    print_dual(color_text(
        f"[OK] Citations for the methods used in this run written to '{BIB_FILE}' "
        f"({len(bib_entries)} entries).", 'green'), f_out)

    print_section("[9] SUMMARY & FILES", f_out)
    print_dual("Status         : OK", f_out)
    print_dual(f"Input file     : {args.filename}", f_out)
    print_dual(f"Output file    : {args.output}", f_out)
    print_dual(f"References     : {BIB_FILE}", f_out)
    if report_path:
        print_dual(f"Report         : {report_path}", f_out)

    if f_out:
        f_out.close()

    # --view runs last, after every check/report section above has already
    # printed, so a blocking GUI window never delays or hides them -- shows
    # the input monolayer and the final tube/ribbon so the user can page
    # through the actual comparison in ase-gui.
    if args.view:
        input_atoms = AseAtomsAdaptor.get_atoms(pmg_before)
        final_atoms = AseAtomsAdaptor.get_atoms(pmg_after)
        view_structure_interactive([input_atoms, final_atoms])


if __name__ == "__main__":
    main()
