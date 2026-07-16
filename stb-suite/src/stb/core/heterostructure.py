"""ZSL commensurate-supercell lattice matching + rigid bilayer construction.

Extracted from stacking2D.py's stack_heterostructure() once stackingfault.py
became a second real consumer: a stacking-fault (generalized stacking fault
energy / gamma-surface) grid sweep needs to build many bilayer structures
(one per lateral-shift grid point) from the SAME two monolayers, and the
expensive part -- the ZSL commensurate-supercell search -- only needs to run
ONCE for the whole sweep, not once per grid point. Splitting stack_
heterostructure() into find_zsl_match() (expensive, once) and
build_stacked_structure() (cheap, call many times with the already-found
match) lets both stacking2D.py (which still calls both once per shift/gap
requested) and stackingfault.py (which calls find_zsl_match() once and
build_stacked_structure() N**2 times for an N x N grid) share this without
either re-running the ZSL search unnecessarily.
"""

import sys

import numpy as np
from pymatgen.core import Structure, Lattice
from pymatgen.core.operations import SymmOp
from pymatgen.analysis.interfaces.zsl import ZSLGenerator

from stb.core.cli import color_text


def find_zsl_match(layer1, layer2, max_area=150.0, max_strain=0.05, match_id=0,
                    interactive=False, twist_angle=0.0):
    """Runs the ZSL commensurate-supercell search ONCE and returns
    (t_mat1, t_mat2, best_match_data) -- the two 3x3 supercell
    transformation matrices for the selected match, plus its area/strain
    diagnostics. Expensive (the ZSL search itself); a caller building many
    structures from the same layer pair (e.g. a shift/gap sweep) should
    call this ONCE and reuse the result via build_stacked_structure(), not
    search again per structure.

    Mutates `layer2` in place (applies the twist rotation, if any) --
    callers must pass this SAME (possibly now-rotated) `layer2` object into
    every subsequent build_stacked_structure() call, not the original
    pre-twist structure, exactly as stack_heterostructure() below does.

    Exits (sys.exit) with a clear error if no lattice match is found at
    all, or if `match_id` is out of range for the non-interactive path --
    same behavior as the original, non-split stack_heterostructure().
    """
    if twist_angle != 0.0:
        print(color_text(f"[INFO] Applying initial twist angle of {twist_angle}° to Layer 2...", 'cyan'))
        op = SymmOp.from_axis_angle_and_translation([0, 0, 1], twist_angle)
        layer2.apply_operation(op, fractional=False)

    print(color_text("[INFO] Searching for commensurate supercells (ZSL Algorithm)...", 'cyan'))
    zsl = ZSLGenerator(max_area=max_area, max_length_tol=max_strain, max_angle_tol=0.1)

    raw_matches = list(zsl(layer1.lattice.matrix[:2], layer2.lattice.matrix[:2]))
    if not raw_matches:
        print(color_text("\n[ERROR] No lattice match found! Try increasing --max_area or --max_strain.", 'red'))
        sys.exit(1)

    evaluated_matches = []
    for m in raw_matches:
        t_mat1 = np.eye(3); t_mat1[:2, :2] = m.film_transformation
        t_mat2 = np.eye(3); t_mat2[:2, :2] = m.substrate_transformation

        t_l1 = layer1.copy(); t_l1.make_supercell(t_mat1)
        t_l2 = layer2.copy(); t_l2.make_supercell(t_mat2)

        # Calculate Linear Strain
        strain_a = abs(t_l1.lattice.a / t_l2.lattice.a - 1.0) * 100
        strain_b = abs(t_l1.lattice.b / t_l2.lattice.b - 1.0) * 100
        max_len_strain = max(strain_a, strain_b)

        # Calculate Angular Strain (Un-rotation penalty)
        vec1 = t_l1.lattice.matrix[0]
        vec2 = t_l2.lattice.matrix[0]
        cos_angle = np.clip(np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2)), -1.0, 1.0)
        angle_diff = np.degrees(np.arccos(cos_angle))

        evaluated_matches.append({
            'match': m,
            'strain': max_len_strain,
            'angle_strain': angle_diff,
            'area': m.match_area
        })

    evaluated_matches.sort(key=lambda x: (x['strain'], x['angle_strain']))

    if interactive:
        display_limit = 15
        print("\n" + color_text(f"--- ZSL Commensurate Supercells Found ({len(evaluated_matches)} matches) ---", 'bold'))
        header = f"{'ID':<4} | {'Area (Å²)':<10} | {'Strain (%)':<10} | {'Ang. Strain (°)':<15} | {'Matrix L1'}"
        print(color_text(header, 'blue'))
        print("-" * 75)

        for i, data in enumerate(evaluated_matches[:display_limit]):
            f_mat = np.round(data['match'].film_transformation).astype(int)
            mat_str = f"[{f_mat[0][0]:2d}, {f_mat[0][1]:2d}], [{f_mat[1][0]:2d}, {f_mat[1][1]:2d}]"
            row_text = f"{i:<4} | {data['area']:<10.2f} | {data['strain']:<10.2f} | {data['angle_strain']:<15.2f} | {mat_str}"
            if data['strain'] < 1.0 and data['angle_strain'] < 1.0:
                print(color_text(row_text, 'green'))
            else:
                print(row_text)

        while True:
            try:
                user_input = input(color_text(f"\nSelect a Match ID (0 to {len(evaluated_matches)-1}) or 'q' to quit: ", 'yellow'))
                if user_input.lower() == 'q':
                    sys.exit(0)
                selected_id = int(user_input)
                if 0 <= selected_id < len(evaluated_matches):
                    match_id = selected_id
                    break
            except ValueError:
                pass
    else:
        if match_id >= len(evaluated_matches) or match_id < 0:
            print(color_text(f"[ERROR] Selected match_id {match_id} is out of range.", 'red'))
            sys.exit(1)

    best_match_data = evaluated_matches[match_id]
    best_match = best_match_data['match']

    t_mat1 = np.eye(3); t_mat1[:2, :2] = best_match.film_transformation
    t_mat2 = np.eye(3); t_mat2[:2, :2] = best_match.substrate_transformation

    print(color_text(f"\n[INFO] Selected Match ID {match_id} | Area: {best_match_data['area']:.2f} Å² | Angular Strain: {best_match_data['angle_strain']:.2f}°", 'green'))
    if twist_angle != 0.0 and best_match_data['angle_strain'] > 1.0:
        print(color_text(f"[WARNING] High Angular Strain! This match will 'un-rotate' your structure by ~{best_match_data['angle_strain']:.1f}° to force PBC fit.", 'yellow'))
        print(color_text("          Increase --max_area to find a true Moiré supercell with 0° Angular Strain.", 'yellow'))

    return t_mat1, t_mat2, best_match_data


def build_stacked_structure(layer1, layer2, t_mat1, t_mat2, shift_x, shift_y, gap,
                             target_vacuum=None, strain_mode='top'):
    """Builds ONE stacked bilayer structure for a given lateral shift of
    layer2 (fractional, in layer2's own primitive cell -- applied before
    supercelling, so the natural periodicity of "how much shift is unique
    before repeating" is set by layer2's own lattice, matching
    stacking2D.py's existing -tx/-ty single-shift convention) and a given
    interlayer gap, using an ALREADY-FOUND ZSL match (t_mat1, t_mat2 from
    find_zsl_match()) -- cheap, safe to call many times in a shift/gap
    sweep. Returns (hetero_structure, n_layer1_atoms, max_strain_val):
    the first n_layer1_atoms atoms of hetero_structure are layer1, the
    rest are layer2 (no z-gap heuristic needed -- this is true by
    construction) -- max_strain_val doesn't depend on shift_x/shift_y/gap
    (only on t_mat1/t_mat2/strain_mode), recomputed here for convenience
    since it's cheap.
    """
    layer1_supercell = layer1.copy(); layer1_supercell.make_supercell(t_mat1)
    layer2_supercell = layer2.copy()
    layer2_supercell.translate_sites(range(len(layer2_supercell)), [shift_x, shift_y, 0.0])
    layer2_supercell.make_supercell(t_mat2)

    # Strain logic
    if strain_mode == 'top':
        base_lattice = layer1_supercell.lattice
    elif strain_mode == 'bottom':
        base_lattice = layer2_supercell.lattice
    else:
        base_lattice = Lattice((layer1_supercell.lattice.matrix + layer2_supercell.lattice.matrix) / 2.0)

    strain_a = max(abs(base_lattice.a / layer1_supercell.lattice.a - 1.0), abs(base_lattice.a / layer2_supercell.lattice.a - 1.0))
    strain_b = max(abs(base_lattice.b / layer1_supercell.lattice.b - 1.0), abs(base_lattice.b / layer2_supercell.lattice.b - 1.0))
    max_strain_val = max(strain_a, strain_b)

    # Flawless mapping of rotated atoms to the new lattice box
    l1_frac = [s.frac_coords for s in layer1_supercell]
    l1_cart_strained = [base_lattice.get_cartesian_coords(f) for f in l1_frac]

    l2_frac = [s.frac_coords for s in layer2_supercell]
    l2_cart_strained = [base_lattice.get_cartesian_coords(f) for f in l2_frac]

    z_max_l1 = max([c[2] for c in l1_cart_strained])
    z_min_l1 = min([c[2] for c in l1_cart_strained])
    z_min_l2 = min([c[2] for c in l2_cart_strained])

    z_shift_cart = z_max_l1 - z_min_l2 + gap
    l2_cart_shifted = [c + np.array([0, 0, z_shift_cart]) for c in l2_cart_strained]

    calc_vacuum = target_vacuum if target_vacuum is not None else layer1_supercell.lattice.c - (z_max_l1 - z_min_l1)
    new_z_max = max([c[2] for c in l2_cart_shifted])
    new_c_length = (new_z_max - z_min_l1) + calc_vacuum

    new_matrix = base_lattice.matrix.copy()
    new_matrix[2] = [0, 0, new_c_length]
    new_lattice = Lattice(new_matrix)

    all_species = [s.specie for s in layer1_supercell] + [s.specie for s in layer2_supercell]
    all_coords_cart = l1_cart_strained + l2_cart_shifted

    hetero = Structure(new_lattice, all_species, all_coords_cart, coords_are_cartesian=True)
    n_layer1_atoms = len(layer1_supercell)

    return hetero, n_layer1_atoms, max_strain_val
