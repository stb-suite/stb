#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.17.0"

import os
import sys
import argparse
import warnings
import statistics
import multiprocessing
import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

# --- Warnings Configuration ---
warnings.filterwarnings("ignore", category=UserWarning, module="pybader")
warnings.filterwarnings("ignore", category=DeprecationWarning)

# --- Library Imports ---
from stb.core.deps import require_sisl, require_pybader
sisl = require_sisl()
PyBaderCalc = require_pybader().Bader
from pybader.io import cube as cube_io

# ================= ANSI COLORS =================
from stb.core.cli import COLORS, color_text, show_intro

# ================= SCIENTIFIC DATA =================

# FALLBACK DICTIONARY
# Contains standard chemical valence (s+p for main group, s+d for transition).
# NOTE: DFT pseudopotentials (especially with semi-core states) may differ.
#       This is used ONLY as a per-species fallback when the .out file can't
#       be parsed, or doesn't mention a species the structure actually has.
FALLBACK_VALENCE = {
    # Period 1
    "H": 1.0,  "He": 2.0,
    # Period 2
    "Li": 1.0, "Be": 2.0, "B": 3.0,  "C": 4.0,  "N": 5.0,  "O": 6.0,  "F": 7.0,  "Ne": 8.0,
    # Period 3
    "Na": 1.0, "Mg": 2.0, "Al": 3.0, "Si": 4.0, "P": 5.0,  "S": 6.0,  "Cl": 7.0, "Ar": 8.0,
    # Period 4
    "K": 1.0,  "Ca": 2.0, "Sc": 3.0, "Ti": 4.0, "V": 5.0,  "Cr": 6.0, "Mn": 7.0, "Fe": 8.0,
    "Co": 9.0, "Ni": 10.0,"Cu": 11.0,"Zn": 12.0,"Ga": 3.0, "Ge": 4.0, "As": 5.0, "Se": 6.0,
    "Br": 7.0, "Kr": 8.0,
    # Period 5
    "Rb": 1.0, "Sr": 2.0, "Y": 3.0,  "Zr": 4.0, "Nb": 5.0, "Mo": 6.0, "Tc": 7.0, "Ru": 8.0,
    "Rh": 9.0, "Pd": 10.0,"Ag": 11.0,"Cd": 12.0,"In": 3.0, "Sn": 4.0, "Sb": 5.0, "Te": 6.0,
    "I": 7.0,  "Xe": 8.0,
    # Period 6 (Lanthanides usually 3, but can vary in DFT)
    "Cs": 1.0, "Ba": 2.0,
    "La": 3.0, "Ce": 4.0, "Pr": 3.0, "Nd": 3.0, "Pm": 3.0, "Sm": 3.0, "Eu": 2.0, "Gd": 3.0,
    "Tb": 3.0, "Dy": 3.0, "Ho": 3.0, "Er": 3.0, "Tm": 3.0, "Yb": 2.0, "Lu": 3.0,
    "Hf": 4.0, "Ta": 5.0, "W": 6.0,  "Re": 7.0, "Os": 8.0, "Ir": 9.0, "Pt": 10.0,"Au": 11.0,
    "Hg": 12.0,"Tl": 3.0, "Pb": 4.0, "Bi": 5.0, "Po": 6.0, "At": 7.0, "Rn": 8.0,
    # Period 7 (Actinides)
    "Fr": 1.0, "Ra": 2.0,
    "Ac": 3.0, "Th": 4.0, "Pa": 5.0, "U": 6.0,  "Np": 5.0, "Pu": 4.0, "Am": 3.0, "Cm": 3.0,
    "Bk": 3.0, "Cf": 3.0, "Es": 3.0, "Fm": 3.0, "Md": 3.0, "No": 2.0, "Lr": 3.0,
    "Rf": 4.0, "Db": 5.0, "Sg": 6.0, "Bh": 7.0, "Hs": 8.0, "Mt": 9.0, "Ds": 10.0,"Rg": 11.0,
    "Cn": 12.0,"Nh": 3.0, "Fl": 4.0, "Mc": 5.0, "Lv": 6.0, "Ts": 7.0, "Og": 8.0
}

def get_zval_from_output(label, override_path=None):
    """
    Reads Z_val from the .out file. Supports Siesta 4.x and 5.x formats.
    """

    if override_path:
        target_file = override_path
        if not os.path.exists(target_file):
            print(color_text(f"   [ERROR] Reference file '{target_file}' (via --ref) not found.", 'red'))
            return None
        print(f"   [INFO] Using external reference file: {color_text(target_file, 'bold')}")
    else:
        target_file = f"{label}.out"
        if not os.path.exists(target_file):
            return None

    dynamic_valence = {}
    current_label = None

    try:
        with open(target_file, 'r') as f:
            for line in f:
                # ---------------- LABEL DETECTION ----------------

                # PRIORITY 1: Siesta 5 Standard ("atom: Called for C")
                # This is the safest method as it appears right before the processing block
                if "atom: Called for" in line:
                    parts = line.split()
                    try:
                        # Find "for" and get the next element
                        idx = parts.index("for")
                        current_label = parts[idx+1]
                        # Clean cleanup (e.g. "C(Z=6)" -> "C")
                        current_label = current_label.split('(')[0]
                    except ValueError:
                        pass

                # PRIORITY 2: Old Standard ("Species number: ... Label: C")
                elif "Species number:" in line and "Label:" in line:
                    parts = line.split()
                    try:
                        lbl_candidate = parts[-1]
                        current_label = lbl_candidate
                    except IndexError:
                        continue

                # ---------------- ZVAL DETECTION ----------------

                # Check if we have an active Label
                if current_label:
                    zval = None

                    # CASE 1: Standard Vna ("Vna: chval, zval: ...")
                    if "Vna: chval, zval:" in line:
                        parts = line.split()
                        try:
                            zval = float(parts[-1])
                        except ValueError: pass

                    # CASE 2: PSML/Pseudopotential Generation
                    # "Valence charge in pseudo generation:    4.00000"
                    elif "Valence charge in pseudo generation:" in line:
                        parts = line.split()
                        try:
                            zval = float(parts[-1])
                        except ValueError: pass

                    if zval is not None:
                        dynamic_valence[current_label] = zval

    except Exception as e:
        print(color_text(f"[WARN] Error reading {target_file}: {e}", 'yellow'))
        return None

    return dynamic_valence if dynamic_valence else None

# Below this raw Bader population (electrons), an atom is flagged as suspicious rather
# than reported at face value. A real atom essentially never integrates to exactly/near
# zero in a normal calculation -- it usually means PyBader found no density feature of
# its own to anchor a basin on (e.g. a pseudopotential with deep semicore states frozen
# into the core, leaving too little valence density near that nucleus to resolve), and
# the atom's whole region got folded into a neighbor's basin instead.
ZERO_POPULATION_TOL = 0.01

# Symmetrically-equivalent atoms (same Wyckoff site) are physically required to carry
# the same charge; a spread larger than this across one such group -- while everything
# else about the calculation looks normal -- is a red flag for numerical noise, a
# resolution problem localized to part of the cell, or a mis-set-up structure, rather
# than a real electronic-structure effect.
SYMMETRY_CHARGE_TOL = 0.1

# ================= HELPERS =================

def get_speed_kwargs(speed_mode):
    """Maps --speed onto PyBader's own bundled 'speed' profile (on-grid
    charge assignment instead of the default, slower near-grid method --
    see its default config.ini's [speed] section) instead of a made-up
    method name that PyBader doesn't actually implement.
    """
    if speed_mode == 'fast':
        return {'method': 'ongrid', 'refine_mode': ('changed', 3), 'speed_flag': True}
    return {}


def resolve_threads(requested):
    """--threads if given; otherwise $SLURM_CPUS_PER_TASK if set (a SLURM
    job's cgroup can cap it well below multiprocessing.cpu_count(), which
    reports the whole host's core count); otherwise every core available.
    """
    if requested is not None:
        return requested
    slurm_cpus = os.environ.get("SLURM_CPUS_PER_TASK")
    if slurm_cpus:
        try:
            return int(slurm_cpus)
        except ValueError:
            pass
    return multiprocessing.cpu_count()


def read_spin_density(file_rho, geometry, cube_path):
    """Returns a PyBader-ready spin-density array for a spin-polarized .RHO
    file (its second grid component), or None for a non-spin-polarized run
    OR if anything about processing the spin component fails -- this is an
    optional add-on to the main charge analysis, so any failure here just
    degrades to "no spin reported" rather than aborting the whole run.
    sisl raises when the second grid component doesn't exist, which is how
    the two cases are told apart; there is no separate flag/header field
    for it. Re-uses PyBader's own cube reader (via a real round-trip
    through sisl's cube writer) instead of hand-rolling the bohr/angstrom
    and charge-density unit conversions PyBader's Bader object expects.
    """
    try:
        spin_grid = sisl.get_sile(file_rho).read_grid(index=1)
    except Exception:
        return None
    try:
        spin_grid.set_geometry(geometry)
        spin_grid.write(cube_path)
        density, _, _, _ = cube_io.read(cube_path)
        return density['charge']
    except Exception as e:
        print(color_text(
            f"   [WARN] Spin-density grid detected but failed to process ({e}) -- "
            "continuing with charge-only analysis.", 'yellow'))
        if os.path.exists(cube_path):
            os.remove(cube_path)
        return None


def find_symmetry_groups(geometry, physical_idx, symprec=1e-3, angle_tolerance=5.0):
    """Returns a list of groups of positions into `physical_idx` (i.e. into
    atoms_data, which is built in the same order) that the detected space
    group treats as symmetrically equivalent -- singleton "groups" (no
    equivalent partner) are dropped, since there's nothing to cross-check.
    Returns an empty list if symmetry detection fails for any reason (e.g. a
    relaxed structure with no exact symmetry left at this tolerance) -- this
    is an optional cross-check, not something worth aborting the run over.
    """
    try:
        lattice = Lattice(geometry.cell)
        species = [geometry.atoms[i].symbol for i in physical_idx]
        coords = [geometry.xyz[i] for i in physical_idx]
        pmg = Structure(lattice, species, coords, coords_are_cartesian=True)
        sga = SpacegroupAnalyzer(pmg, symprec=symprec, angle_tolerance=angle_tolerance)
        sym_struct = sga.get_symmetrized_structure()
        return [list(g) for g in sym_struct.equivalent_indices if len(g) > 1]
    except Exception:
        return []

# ================= MAIN LOGIC =================

def solve_bader(label, output_file=None, speed_mode='normal', ref_file=None,
                 threads=None, vacuum_tol=None, keep_cube=True, export_volumes=False):

    file_rho = f"{label}.RHO"
    file_xv = f"{label}.XV"
    file_fdf = f"{label}.fdf"
    file_cube = f"{label}.cube"
    file_spin_cube = f"{label}_spin.cube"

    if not output_file:
        output_file = f"{label}_BADER.txt"

    if not os.path.exists(file_rho):
        print(color_text(f"[ERROR] Grid file '{file_rho}' not found.", 'red'))
        sys.exit(1)

    print(f"[INFO] System: {color_text(label, 'bold')} | Mode: {color_text(speed_mode.upper(), 'cyan')}")
    if speed_mode == 'fast':
        print(color_text(
            "   [CAUTION] --speed fast trades basin-boundary accuracy for speed (PyBader's "
            "on-grid method instead of near-grid) -- prefer 'normal' for a final result.",
            'yellow'))

    # --- Z_val Setup ---
    print(f"0. [Setup] Configuring valence charges...")

    detected_valence = get_zval_from_output(label, override_path=ref_file)

    if detected_valence:
        print(f"   {color_text('[SUCCESS]', 'green')} Z_vals detected: {detected_valence}")
        source_name = f"Siesta Output ({ref_file if ref_file else label+'.out'})"
    else:
        detected_valence = {}
        # --- HIGH VISIBILITY WARNING BLOCK ---
        warn_msg = [
            "\n" + "!" * 65,
            "[WARNING] .out FILE NOT FOUND OR UNREADABLE! USING DEFAULT VALUES.",
            "Please ensure your pseudopotential Z_val matches standard values.",
            "Use --ref to point to a valid .out file.",
            "!" * 65 + "\n"
        ]
        # Print with RED Background and WHITE text
        for line in warn_msg:
            print(f"{COLORS['bg_red']}{COLORS['white']}{COLORS['bold']}{line.center(65)}{COLORS['reset']}")

        source_name = "Hardcoded Defaults (FALLBACK)"

    # Detected values (from the actual pseudopotential used) always win; FALLBACK_VALENCE
    # only fills in species the .out parse didn't mention -- so a PARTIAL parse (some
    # elements detected, others missed, e.g. an unrecognized log-format for one species)
    # doesn't get silently treated as fully successful, and the missed species aren't just
    # dropped from the analysis if a reasonable (if imprecise) fallback exists for them.
    valence_source = {**FALLBACK_VALENCE, **detected_valence}

    cube_files = []
    try:
        # --- STEP 1: SISL ---
        try:
            print(f"1. [SISL] Reading geometry and charge density...")

            if os.path.exists(file_xv):
                geometry = sisl.get_sile(file_xv).read_geometry()
            elif os.path.exists(file_fdf):
                geometry = sisl.get_sile(file_fdf).read_geometry()
            else:
                print(color_text("[ERROR] No geometry file (.XV or .fdf) found.", 'red'))
                sys.exit(1)

            # Atoms with no real element assigned (Z<=0 -- e.g. a floating dummy site with
            # no basis) are silently dropped by sisl's own cube writer, so they're excluded
            # here too to stay consistent with what PyBader will actually see. This does NOT
            # affect SIESTA's own ghost/BSSE atoms: sisl represents those as `AtomGhost`,
            # whose `.Z` is the real (positive) element number -- verified empirically against
            # this installed sisl, sisl's cube writer keeps them, and so does this check.
            physical_idx = [i for i, atom in enumerate(geometry.atoms) if atom.Z > 0]
            dummy_count = len(geometry.atoms) - len(physical_idx)
            if dummy_count:
                print(f"   {color_text('[INFO]', 'cyan')} {dummy_count} dummy/no-element site(s) "
                      "(Z<=0) excluded from the cube file and this analysis.")

            if detected_valence:
                structure_syms = {geometry.atoms[i].symbol for i in physical_idx}
                fallback_syms = sorted(s for s in structure_syms
                                       if s not in detected_valence and s in FALLBACK_VALENCE)
                if fallback_syms:
                    print(color_text(
                        f"   [WARN] Element(s) {', '.join(fallback_syms)} not detected in the "
                        ".out parse -- using the hardcoded fallback Z_val for them instead "
                        "(may not match the actual pseudopotential used).", 'yellow'))

            rho_grid = sisl.get_sile(file_rho).read_grid()
            rho_cell = rho_grid.cell.copy()
            rho_grid.set_geometry(geometry)
            rho_grid.write(file_cube)
            cube_files.append(file_cube)
            density, lattice, cube_atoms, file_info = cube_io.read(file_cube)

            spin_density = read_spin_density(file_rho, geometry, file_spin_cube)
            if spin_density is not None and spin_density.shape == density['charge'].shape:
                density['spin'] = spin_density
                cube_files.append(file_spin_cube)
                print(f"   {color_text('[INFO]', 'cyan')} Spin-polarized grid detected -- "
                      "reporting per-atom net spin (magnetic moment) too.")
            elif spin_density is not None:
                print(color_text(
                    "   [WARN] Spin-density grid shape doesn't match the charge grid -- "
                    "skipping spin-resolved analysis.", 'yellow'))

        except Exception as e:
            print(color_text(f"[ERROR] SISL processing failed: {e}", 'red'))
            sys.exit(1)

        # --- Cross-checks: two independent signals that the .RHO really belongs with
        # this geometry. The cell check compares the .RHO file's OWN lattice (captured
        # above, before set_geometry() overwrote it) against the geometry's -- this catches
        # a genuinely wrong pairing (e.g. a different relaxation step with a different
        # cell). The atom check compares the geometry against what PyBader will actually
        # see in the cube file, so a sisl/PyBader disagreement doesn't silently misattribute
        # charges. Neither check (nor both together) can catch a same-cell, different-
        # atomic-position mismatch -- a .RHO grid carries no independent atomic-position
        # record to compare against; that class of mistake is not detectable from these
        # files alone.
        if not np.allclose(rho_cell, geometry.cell, atol=1e-3):
            print(color_text(
                "[ERROR] Lattice mismatch between the .RHO grid and the geometry file -- "
                "they likely don't belong to the same calculation.", 'red'))
            sys.exit(1)

        if len(cube_atoms) != len(physical_idx):
            print(color_text(
                f"[ERROR] Atom count mismatch: the geometry has {len(physical_idx)} real "
                f"atom(s) (excluding dummy/no-element sites) but the written cube file has "
                f"{len(cube_atoms)} -- refusing to continue with possibly misaligned charges.",
                'red'))
            sys.exit(1)
        cube_z = [int(z) for z in file_info['elements']]
        geo_z = [geometry.atoms[i].Z for i in physical_idx]
        if cube_z != geo_z:
            print(color_text(
                "[ERROR] Atom order/species mismatch between the geometry and the written "
                "cube file -- refusing to continue with possibly misaligned charges.", 'red'))
            sys.exit(1)

        # --- STEP 2: PyBader ---
        n_threads = resolve_threads(threads)
        print(f"2. [PyBader] Starting calculation on {n_threads} threads...")

        # 'output': None stops PyBader's own __call__ from silently pickling the whole
        # calculation to a 'bader.p' file in the working directory -- this tool builds
        # its own report, so that file would just be unexplained clutter.
        extra_kwargs = {'threads': n_threads, 'vacuum_tol': vacuum_tol, 'output': None}
        extra_kwargs.update(get_speed_kwargs(speed_mode))
        if 'spin' in density:
            extra_kwargs['spin_flag'] = True
        if export_volumes:
            # 'atoms' (not 'volumes') so this stays valid even under --speed fast, which
            # deletes bader_volumes but keeps atoms_volumes (see PyBader's own __call__).
            # [-2] is PyBader's own sentinel for "every atom, plus the vacuum bucket if
            # vacuum_tol is set".
            extra_kwargs['export_mode'] = ('atoms', [-2])

        try:
            bader_job = PyBaderCalc(density, lattice, cube_atoms, file_info, **extra_kwargs)
            bader_job()
            raw_populations = bader_job.atoms_charge
            raw_spins = bader_job.atoms_spin if bader_job.spin_bool else None
            raw_volumes = bader_job.atoms_volume
            raw_surf_dist = bader_job.atoms_surface_distance
        except Exception as e:
            print(color_text(f"[ERROR] PyBader calculation failed: {e}", 'red'))
            sys.exit(1)

        if export_volumes:
            print(f"   {color_text('[INFO]', 'cyan')} Exported {len(cube_atoms)} per-atom "
                  f"Bader volume(s) as 'Bader-atoms-<N>.cube' in the current directory.")

        # --- STEP 3 & 4: Analysis and Output ---
        try:
            print("3. [Analysis] compiling results...")

            total_theory = 0.0
            total_raw_known = 0.0
            atoms_data = []
            unknown_syms = set()
            suspicious_zero_ids = []

            n_atoms = min(len(physical_idx), len(raw_populations))
            if len(physical_idx) != len(raw_populations):
                print(color_text(
                    f"   [WARN] {len(physical_idx)} real atom(s) in the geometry but PyBader "
                    f"only returned {len(raw_populations)} population(s) -- only the first "
                    f"{n_atoms} will be reported. Check that .RHO/.XV/.fdf really belong to "
                    "the same run.", 'red'))
            if raw_spins is not None and len(raw_spins) != len(raw_populations):
                print(color_text(
                    "   [WARN] PyBader returned a different number of spin values than charge "
                    "values -- disabling the spin column for this run.", 'yellow'))
                raw_spins = None
            if len(raw_volumes) != len(raw_populations) or len(raw_surf_dist) != len(raw_populations):
                print(color_text(
                    "   [WARN] PyBader returned a different number of volume/surface-distance "
                    "values than charge values -- disabling those columns for this run.", 'yellow'))
                raw_volumes = raw_surf_dist = None

            for pos in range(n_atoms):
                orig_i = physical_idx[pos]
                atom = geometry.atoms[orig_i]
                sym = atom.symbol
                z_val = valence_source.get(sym)
                spin_val = raw_spins[pos] if raw_spins is not None else None
                volume_val = raw_volumes[pos] if raw_volumes is not None else None
                surf_dist_val = raw_surf_dist[pos] if raw_surf_dist is not None else None

                if z_val is None:
                    unknown_syms.add(sym)
                else:
                    total_theory += z_val
                    total_raw_known += raw_populations[pos]

                atoms_data.append({'id': orig_i + 1, 'sym': sym, 'z_val': z_val,
                                    'pop_raw': raw_populations[pos], 'spin': spin_val,
                                    'volume': volume_val, 'surf_dist': surf_dist_val})
                if raw_populations[pos] < ZERO_POPULATION_TOL:
                    suspicious_zero_ids.append(orig_i + 1)

            if suspicious_zero_ids:
                ids_str = ', '.join(str(i) for i in suspicious_zero_ids)
                print(color_text(
                    f"   [WARN] Atom(s) {ids_str} got essentially zero Bader population "
                    f"(< {ZERO_POPULATION_TOL} e-) -- PyBader found almost no density of its "
                    "own to assign there. This is rarely a trustworthy result; common causes "
                    "are a pseudopotential with no resolvable density near the nucleus (e.g. "
                    "deep semicore states frozen into the core) or the atom's whole region "
                    "being folded into a neighboring basin. Treat this atom's charge with real "
                    "suspicion, not as a precise result.", 'red'))

            if unknown_syms:
                print(color_text(
                    f"   [WARN] Element(s) {', '.join(sorted(unknown_syms))} not found in the "
                    "Z_val dictionary! Their net charge cannot be computed and they're excluded "
                    "from the correction-factor/total below (use --ref to point at a matching "
                    ".out file).", 'red'))

            # Unit Correction -- computed from atoms with a known Z_val only, so an
            # unrecognized element can't silently corrupt the correction factor for every
            # other atom. This is a single global scalar assuming a uniform cause; if the
            # real deviation is spatially localized (a coarse mesh region, one bad atom,
            # a wrong file pairing), it will instead smear that localized error evenly
            # across every atom's reported charge -- a large factor is a cue to
            # investigate further, not a guarantee the correction is physically right.
            ratio = total_raw_known / total_theory if total_theory > 0 else 1.0
            correction_factor = 1.0 / ratio if abs(ratio - 1.0) > 0.1 else 1.0

            if correction_factor != 1.0:
                print(color_text(
                    f"   [INFO] Unit mismatch corrected. Factor: {correction_factor:.4f} "
                    "(applied uniformly to every atom -- assumes a global cause; see the "
                    "note at the end of the report).", 'cyan'))

            # Symmetry cross-check -- atoms the detected space group treats as equivalent
            # should carry the same charge; a run-time deviation here is a red flag that's
            # invisible from looking at any atom's own row in isolation.
            inconsistent_groups = []
            for group in find_symmetry_groups(geometry, physical_idx):
                group = [p for p in group if p < n_atoms and atoms_data[p]['z_val'] is not None]
                if len(group) < 2:
                    continue
                charges = [atoms_data[p]['z_val'] - atoms_data[p]['pop_raw'] * correction_factor
                           for p in group]
                if max(charges) - min(charges) > SYMMETRY_CHARGE_TOL:
                    ids = [atoms_data[p]['id'] for p in group]
                    inconsistent_groups.append((ids, charges))

            if inconsistent_groups:
                for ids, charges in inconsistent_groups:
                    pairs = ', '.join(f"#{i}={c:+.3f}" for i, c in zip(ids, charges))
                    print(color_text(
                        f"   [WARN] Symmetry-equivalent atoms disagree on net charge: {pairs} "
                        f"(spread > {SYMMETRY_CHARGE_TOL} e-) -- these sites should be identical "
                        "by symmetry; investigate before trusting either value.", 'red'))

            # --- STEP 4: Output ---
            has_spin = raw_spins is not None

            has_volume = raw_volumes is not None

            out_lines = []
            out_lines.append(f"BADER CHARGE ANALYSIS REPORT - STB Suite v{VERSION}")
            out_lines.append(f"System: {label}")
            out_lines.append(f"Z_val Source: {source_name}")
            out_lines.append(f"PyBader: method={bader_job.method}, refine_method={bader_job.refine_method}, "
                              f"threads={bader_job.threads}, vacuum_tol={bader_job.vacuum_tol}")

            vol_header = f" {'Vol(A^3)':<10} {'SurfD(A)':<10}" if has_volume else ""
            spin_header = f" {'Spin (µB)':<10}" if has_spin else ""
            header = (f"{'Idx':<4} {'Elem':<5} {'Pop(e-)':<12} {'Z_val':<8} {'Net Charge':<12} "
                      f"{'State':<15}{vol_header}{spin_header}")
            sep = "-" * len(header)
            eq = "=" * len(header)

            out_lines.append(eq)
            out_lines.append(header)
            out_lines.append(sep)

            print("\n" + header)
            print(sep)

            total_final = 0.0

            for data in atoms_data:
                pop = data['pop_raw'] * correction_factor
                total_final += pop

                if data['z_val'] is None:
                    z_str = f"{'N/A':<8}"
                    net_str = f"{'N/A':<12}"
                    state, color = "Unknown Z_val", 'yellow'
                else:
                    net = data['z_val'] - pop
                    z_str = f"{data['z_val']:<8.2f}"
                    net_str = f"{net:<+12.4f}"
                    if net > 0.05: state, color = "Donor (+)", 'red'
                    elif net < -0.05: state, color = "Acceptor (-)", 'blue'
                    else: state, color = "Neutral", 'reset'

                vol_str = f" {data['volume']:<10.4f} {data['surf_dist']:<10.4f}" if has_volume else ""
                spin_str = f" {data['spin']:<+10.4f}" if has_spin else ""

                line = (f"{data['id']:<4} {data['sym']:<5} {pop:<12.4f} {z_str} {net_str} "
                        f"{state:<15}{vol_str}{spin_str}")
                out_lines.append(line)
                print(f"{data['id']:<4} {data['sym']:<5} {pop:<12.4f} {z_str} "
                      f"{color_text(net_str, 'bold')} {color_text(state, color)}{vol_str}{spin_str}")

            theory_note = "" if not unknown_syms else \
                f" (excludes {len(unknown_syms)} unknown element(s): {', '.join(sorted(unknown_syms))})"
            footer = [
                sep,
                f"Total Integrated: {total_final:.4f} (Target: {total_theory:.2f}{theory_note})",
                eq,
            ]
            out_lines.extend(footer)
            for l in footer: print(l)

            # Per-species summary -- a quick way to spot a species whose charge is wildly
            # inconsistent across the cell without reading every row (population standard
            # deviation only, since a single outlier atom is exactly what this is meant to
            # surface -- a full statistical treatment isn't the point here).
            by_species = {}
            for data in atoms_data:
                if data['z_val'] is None:
                    continue
                net = data['z_val'] - data['pop_raw'] * correction_factor
                by_species.setdefault(data['sym'], []).append(net)

            if by_species:
                species_header = f"{'Elem':<5} {'N':<4} {'Mean(e-)':<12} {'Std(e-)':<10}"
                species_sep = "-" * len(species_header)
                species_lines = ["", "Per-species net charge summary:", species_header, species_sep]
                for sym in sorted(by_species):
                    values = by_species[sym]
                    mean = statistics.mean(values)
                    std = statistics.pstdev(values) if len(values) > 1 else 0.0
                    species_lines.append(f"{sym:<5} {len(values):<4} {mean:<+12.4f} {std:<10.4f}")
                out_lines.extend(species_lines)
                print("\n" + "\n".join(species_lines))

            if suspicious_zero_ids:
                ids_str = ', '.join(str(i) for i in suspicious_zero_ids)
                out_lines.append(
                    f"\nWARN: Atom(s) {ids_str} got essentially zero Bader population (< "
                    f"{ZERO_POPULATION_TOL} e-) -- treat with suspicion, see console output "
                    "for likely causes."
                )

            if inconsistent_groups:
                for ids, charges in inconsistent_groups:
                    pairs = ', '.join(f"#{i}={c:+.3f}" for i, c in zip(ids, charges))
                    out_lines.append(
                        f"WARN: Symmetry-equivalent atoms disagree on net charge: {pairs} "
                        f"(spread > {SYMMETRY_CHARGE_TOL} e-)."
                    )

            limitations_note = (
                "Note: Bader analysis has known limitations this tool cannot detect or correct "
                "for -- non-nuclear attractors (basins not centered on any atom, common in "
                "ionic/metallic systems) are always folded into the nearest atom by PyBader; "
                "Z_val is only as accurate as the .out parse or the FALLBACK_VALENCE table "
                "(semicore-inclusive pseudopotentials can differ from the 'standard valence' "
                "assumed there); and --speed fast trades basin-boundary accuracy for speed."
            )
            out_lines.append(limitations_note)
            print(color_text("\n  " + limitations_note, 'yellow'))

        except Exception as e:
            print(color_text(f"[ERROR] Analysis/report generation failed: {e}", 'red'))
            sys.exit(1)

        try:
            with open(output_file, "w") as f:
                f.write("\n".join(out_lines))
            print(f"\n{color_text('[OK]', 'green')} Results saved to: {color_text(output_file, 'bold')}")
        except IOError:
            print(color_text(f"[ERROR] Could not save file {output_file}", 'red'))
            sys.exit(1)

    finally:
        if not keep_cube:
            for f in cube_files:
                if os.path.exists(f):
                    os.remove(f)

# ================= EXECUTION =================

def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Computes Bader (AIM) atomic charges from a SIESTA charge-density grid.", 'bold')}
Reads <label>.RHO (charge density) plus <label>.XV or <label>.fdf (geometry),
converts to a .cube file via sisl, and partitions it into atomic basins with
PyBader. Z_val (valence electron count per species) is parsed automatically
from <label>.out when available (or --ref points at a specific SIESTA
output); species it doesn't mention fall back to a hardcoded periodic-table
guess, per species, with a warning. If the .RHO file is spin-polarized, a
per-atom net spin (magnetic moment) is reported alongside the charge
automatically. Dummy/no-element sites (Z<=0) are excluded to match
PyBader's own handling; SIESTA's own ghost/BSSE atoms are unaffected.
Also reports each atom's Bader volume and minimum surface distance, a
per-species mean/std summary, and cross-checks symmetry-equivalent atoms
against each other -- on top of the per-atom warnings for a suspiciously
low (near-zero) population.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s --label siesta\n"
               "  %(prog)s --label siesta --ref relax/siesta.out -o siesta_bader.txt\n"
               "  %(prog)s --label siesta --speed fast --threads 8\n"
               "  %(prog)s --label slab --vacuum-tol 1e-3\n"
               "  %(prog)s --label siesta --export-volumes\n"
    )
    parser.add_argument("-l", "--label", required=True, help="SystemLabel used in Siesta")
    parser.add_argument("-o", "--output", required=False, help="Output filename (default: <label>_BADER.txt)")
    parser.add_argument("--ref", required=False, default=None,
                         help="Path to a specific .out file to read Z_val from (overrides <label>.out)")
    parser.add_argument("--speed", choices=['normal', 'fast'], default='normal',
                         help="'fast' uses PyBader's own bundled speed profile (on-grid charge "
                              "assignment instead of near-grid) -- quicker, less precise basin "
                              "boundaries (default: normal)")
    parser.add_argument("--threads", type=int, default=None, metavar="N",
                         help="Worker threads for PyBader, >= 1 (default: $SLURM_CPUS_PER_TASK if "
                              "set, else every CPU core -- multiprocessing.cpu_count() can "
                              "overcount under a job scheduler's cgroup limits)")
    parser.add_argument("--vacuum-tol", type=float, default=None, metavar="TOL",
                         help="Charge-density threshold (e/Ang^3) below which a voxel is assigned "
                              "to a separate 'vacuum' bucket instead of the nearest atom -- useful "
                              "for slabs/wires/isolated molecules with empty space in the cell "
                              "(default: disabled, matches PyBader's own default)")
    parser.add_argument("--no-cube", dest="keep_cube", action="store_false", default=True,
                         help="Delete the intermediate .cube file(s) after the calculation "
                              "(kept by default -- they can be large for dense grids)")
    parser.add_argument("--export-volumes", action="store_true",
                         help="Also write each atom's individual Bader volume as its own "
                              "'Bader-atoms-<N>.cube' file (plus a vacuum one if --vacuum-tol "
                              "is set), for visual inspection in VESTA/VMD -- one file per atom, "
                              "each as large as the main grid, so this is opt-in (default: off)")
    parser.add_argument("-v", "--version", action="version",
                        version=f"stb-bader {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.threads is not None and args.threads < 1:
        parser.error("--threads must be at least 1.")

    if args.intro:
        show_intro(["Siesta ToolBox Suite - Bader Analysis", f"Version {VERSION} | University of Brasilia"])

    solve_bader(args.label, args.output, args.speed, args.ref,
                threads=args.threads, vacuum_tol=args.vacuum_tol, keep_cube=args.keep_cube,
                export_volumes=args.export_volumes)

    print("\n" + "-" * 60)
    print(color_text("Electron counting is like accounting, but the currency is negative.\n", 'bold'))

if __name__ == "__main__":
    main()
