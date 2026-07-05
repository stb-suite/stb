#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.9.1"

import sys
import os
import argparse
import numpy as np
from stb.core import structure_io, kspace
from stb.core.cli import COLORS, color_text, show_intro


def find_translation_vector(a1, a2, n, m, search_bound=None):
    """Shortest lattice vector T = t1*a1 + t2*a2 (t1, t2 integers) that is
    not parallel to C = n*a1 + m*a2, found by a bounded search -- no
    hexagonal-only closed form (the classic CNT gcd/(n-m)%3d shortcut only
    holds for a honeycomb lattice), so this works for any oblique 2D
    lattice. `n*t2 - m*t1` is the exact integer cross product of C and T in
    the (a1, a2) basis: zero iff they are parallel.
    """
    if search_bound is None:
        search_bound = 3 * (abs(n) + abs(m)) + 5

    candidates = []
    for t1 in range(-search_bound, search_bound + 1):
        for t2 in range(-search_bound, search_bound + 1):
            cross = n * t2 - m * t1
            if cross == 0:
                continue
            length = float(np.linalg.norm(t1 * a1 + t2 * a2))
            candidates.append((length, t1, t2, cross))

    if not candidates:
        raise ValueError(
            f"Could not find a translation vector within search bound {search_bound}; "
            "try increasing --search-bound."
        )

    # Shortest first; among near-ties prefer a positive cross product (canonical orientation).
    candidates.sort(key=lambda c: (round(c[0], 6), 0 if c[3] > 0 else 1))
    _, t1, t2, _ = candidates[0]
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

    M = np.array([[n, m], [t1, t2]], dtype=float)
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
    parser.add_argument("-o", "--output", type=str, default="nanotube.fdf",
                        help="Output .fdf file name (default: nanotube.fdf).")
    parser.add_argument("-v", "--version", action="version", version=f"stb-nanotube {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("Roll a 2D monolayer into a nanotube/nanoribbon:", 'bold'))
    print("-" * 60)

    if not os.path.exists(args.filename):
        print(color_text(f"Error: File '{args.filename}' not found.", 'red'))
        sys.exit(1)

    n, m = args.chirality
    if (n, m) == (0, 0):
        print(color_text("Error: Chirality (0, 0) is not valid.", 'red'))
        sys.exit(1)

    if args.repeats < 1:
        print(color_text("Error: --repeats must be >= 1.", 'red'))
        sys.exit(1)

    try:
        structure = structure_io.read_fdf(args.filename)
    except (FileNotFoundError, ValueError) as e:
        print(color_text(f"Error: {e}", 'red'))
        sys.exit(1)

    positions = np.array([pos for _, pos in structure.atoms])
    is_cartesian = structure.coord_format == 'cartesian'
    frac_coords = kspace.to_fractional(positions, structure.lattice, is_cartesian)
    vacuum_axes = kspace.detect_vacuum_axes(frac_coords, structure.lattice, args.vacuum_gap)

    if sum(vacuum_axes) != 1:
        dim = 3 - sum(vacuum_axes)
        print(color_text(
            f"Error: input must be a 2D monolayer (exactly 1 vacuum-padded axis); "
            f"detected {dim}D ({sum(vacuum_axes)} vacuum axes). Adjust --vacuum-gap "
            "if this looks wrong.", 'red'))
        sys.exit(1)

    vac_axis_idx = vacuum_axes.index(True)
    in_plane_axes = [i for i in range(3) if i != vac_axis_idx]
    a1 = structure.lattice[in_plane_axes[0]]
    a2 = structure.lattice[in_plane_axes[1]]

    print(f"  {color_text('Input formula:', 'cyan')} "
          f"{structure_io.to_pymatgen(structure).composition.reduced_formula}")
    print(f"  {color_text('Input atoms:', 'cyan')} {len(structure.atoms)}")
    print(f"  {color_text('In-plane axes:', 'cyan')} {in_plane_axes}  "
          f"{color_text('Vacuum axis:', 'cyan')} {vac_axis_idx}")
    print(f"  {color_text('Chirality (n, m):', 'cyan')} ({n}, {m})")

    try:
        t1, t2 = find_translation_vector(a1, a2, n, m, args.search_bound)
    except ValueError as e:
        print(color_text(f"Error: {e}", 'red'))
        sys.exit(1)

    C_cart = n * a1 + m * a2
    T_cart = t1 * a1 + t2 * a2

    try:
        cell_atoms, n_cells = enumerate_cell_atoms(
            structure, in_plane_axes, vac_axis_idx, n, m, t1, t2, args.lattice_tol
        )
    except ValueError as e:
        print(color_text(f"Error: {e}", 'red'))
        sys.exit(1)

    print(f"  {color_text('Translation vector (t1, t2):', 'cyan')} ({t1}, {t2})")
    print(f"  {color_text('Cells in periodic unit (N_cells):', 'cyan')} {n_cells}")

    if args.mode == "tube":
        atoms, lattice, R = roll_to_tube(cell_atoms, C_cart, T_cart, args.repeats, args.min_vacuum_size)
        print(f"\n  {color_text('Circumference |C|:', 'cyan')} {np.linalg.norm(C_cart):.4f} Ang")
        print(f"  {color_text('Tube diameter (2R):', 'cyan')} {2 * R:.4f} Ang")
        print(f"  {color_text('Axial period |T|:', 'cyan')} {np.linalg.norm(T_cart):.4f} Ang")
    else:
        atoms, lattice, width = build_ribbon(
            cell_atoms, C_cart, T_cart, vac_axis_idx, structure.lattice, args.repeats, args.min_vacuum_size
        )
        print(f"\n  {color_text('Ribbon width:', 'cyan')} {width:.4f} Ang")
        print(f"  {color_text('Axial period |T|:', 'cyan')} {np.linalg.norm(T_cart):.4f} Ang")

    new_structure = build_result_structure(atoms, lattice, structure.species_meta)
    print(f"  {color_text('Output formula:', 'cyan')} "
          f"{structure_io.to_pymatgen(new_structure).composition.reduced_formula}")
    print(f"  {color_text('Output atoms:', 'cyan')} {len(atoms)}")

    structure_io.write_fdf(new_structure, args.output)
    print(f"\n{color_text('Success:', 'green')} Structure written to '{color_text(args.output, 'bold')}'")


if __name__ == "__main__":
    main()
