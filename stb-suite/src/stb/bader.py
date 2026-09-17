#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "2.0.1"

import os
import sys
import argparse
import warnings
import statistics
import multiprocessing
from datetime import datetime
import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

# --- Warnings Configuration ---
warnings.filterwarnings("ignore", category=UserWarning, module="pybader")
warnings.filterwarnings("ignore", category=DeprecationWarning)

# --- Library Imports ---
from stb.core.deps import require_sisl, require_pybader, read_sisl_geometry_xv_or_fdf
sisl = require_sisl()
PyBaderCalc = require_pybader().Bader
from pybader.io import cube as cube_io

# ================= ANSI COLORS =================
from stb.core import citations
from stb.core.cli import COLORS, color_text, show_intro, print_dual, print_section, print_table
# FALLBACK_VALENCE/get_zval_from_output moved to core/siesta_log.py once
# hirshfeld_ions.py/hirshfeld_analysis.py became further consumers of the
# exact same "real, pseudopotential-consistent per-species Z_val" need
# (extract-on-second-use, same policy as the rest of core/).
from stb.core.siesta_log import FALLBACK_VALENCE, get_zval_from_output

REPORT_FILE = "stb_bader_report.txt"
BIB_FILE = "references.bib"

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

# A whole species reporting a LARGE, nearly-uniform net charge across every
# atom of that species (not just one outlier) is a distinct failure
# signature from either check above: a pseudopotential with no resolvable
# valence density maximum at that nucleus (no NLCC, deep frozen semicore
# states) makes Bader donate that species' basin volume to its neighbors
# wholesale, so every atom of the species fails identically rather than
# randomly. Real, physically large and uniform charges do happen (e.g. an
# alkali halide) -- these thresholds are a coarse heuristic, not proof of a
# bug, and may need tuning against more cases than the one that motivated
# them (a Te pseudopotential without NLCC reading a uniform +2.22 e across
# every Te atom in a slab; see
# new_functions/PROPOSAL_hirshfeld-I_implementation.md).
SPECIES_CHARGE_MAGNITUDE_TOL = 1.0    # e-: mean |net charge| above this is "large"
SPECIES_CHARGE_UNIFORMITY_TOL = 0.15  # e-: pstdev below this is "suspiciously uniform"

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


def build_bader_inputs(grid, geometry, physical_idx, cube_filename, output_dir):
    """Builds PyBader's 4 required inputs -- (density, lattice, atoms,
    file_info), exactly the tuple `pybader.io.cube.read()` would return --
    directly from an already-read, in-memory sisl `Grid` + `Geometry`,
    with NO round-trip through a written `.cube` file's text.

    **Two real, verified bugs this sidesteps entirely** (both in the
    write-via-sisl / read-via-pybader round-trip this function replaces,
    not in pybader's own cube.read() in isolation):

    1. **Crash for any grid whose nz isn't a multiple of 6.** sisl's cube
       writer (`Grid.write`) streams the ENTIRE flattened grid as one
       continuous run of values, 6 per line, wrapping across (x,y) row
       boundaries -- a valid Gaussian-cube convention, but NOT the one
       `pybader.io.cube.read()` assumes: it expects every (x,y) row to
       start on its own fresh line and end with its own (possibly short)
       trailing line before the next row begins. Confirmed via
       `bugs_report/BUG_REPORT_stb-bader_cube-roundtrip.md`: on a
       256-point z-grid (256 %% 6 = 4), this raises "could not broadcast
       input array from shape (6,) into shape (4,)".
    2. **Silent ~6.75x charge overestimate for every grid where nz IS a
       multiple of 6** (i.e. whenever bug #1 doesn't crash) -- NOT
       documented in the bug report above, found and verified numerically
       while fixing it (test/3-analysis/6-bader/Sn3O4.RHO, a real
       committed fixture: nz=90, 90%%6=0, so bug #1 never triggers here).
       sisl's `read_grid()` scales RHO/BADER values by `sile.grid_unit`
       (~6.748334, the Bohr^-3 -> Ang^-3 conversion -- see stb-cube.py's
       own `[WARNING]`-free, already-correct precedent for this exact
       factor) so `grid.grid` in memory is already e/Ang^3. sisl's cube
       WRITER never undoes that scaling (its own docstring: "grid data is
       assumed to be unit-less"), so the written file pairs Bohr-unit
       coordinates with e/Ang^3-valued density -- but `pybader.io.cube.
       read()` then multiplies by `ang_to_bohr**3` (numerically identical
       to `sile.grid_unit`) on the assumption the file holds native
       e/Bohr^3 values, compounding the error: verified directly on
       Sn3O4.RHO, the OLD round-trip's total integrated charge came out
       485.88 e against `grid.grid.sum() * voxel_volume` = 72.00 e from
       the same in-memory grid -- a 485.88 / 72.00 = 6.748 ratio, matching
       `sile.grid_unit` exactly. PyBader's own `charge_sum` (pybader/
       utils.py: `charge[atom_num] += density[i]`, then `*= voxel_volume`)
       confirms it wants plain density (e/Ang^3), not a pre-scaled
       quantity -- despite the `Bader` class docstring's own misleading
       "(rho * lattice volume) units" wording.

    `grid.grid` (whatever sisl's default `read_grid()` scaling left it in
    -- e/Ang^3 for a RHO/BADER file) is used AS-IS, matching what the
    (bug-free) round-trip would have produced after pybader's own
    `*= ang_to_bohr**3` undoes a correctly-applied write-side
    `/= sile.grid_unit` -- i.e. this is not a new/different convention,
    just the same one reached without ever touching disk or pybader's
    own (buggy) reader. `lattice`/`atoms` need no Bohr round-trip either:
    sisl's `grid.cell` is already the full cell matrix in Ang (cube.read()
    only reconstructs this from a per-voxel Bohr vector because that's
    literally how the file format stores it), and wrapping `geometry.fxyz`
    into [0, 1) before converting to cartesian replicates cube.read()'s
    own wrap-into-cell step exactly.
    """
    density = {'charge': np.array(grid.grid, dtype=np.float64)}
    lattice = np.array(grid.cell, dtype=np.float64)
    atoms = (geometry.fxyz[physical_idx] % 1.0) @ lattice
    file_info = {
        'filename': cube_filename,
        'prefix': os.path.join(output_dir, ''),
        'file_type': 'cube',
        'write_function': cube_io.write,
        'elements': np.array([geometry.atoms[i].Z for i in physical_idx], dtype=np.int64),
        'voxel_offset': np.array([.5, .5, .5]),
    }
    return density, lattice, atoms, file_info


def write_native_cube(grid, sile, cube_path):
    """Writes `grid` (already read via `sile.read_grid(...)`, still
    carrying sisl's default e/Ang^3-scaled RHO/BADER values in `.grid`) to
    `cube_path` in genuine, standard-compliant Gaussian-cube units
    (e/Bohr^3, correctly paired with the writer's own Bohr-unit
    coordinates) -- the same `grid.grid / sile.grid_unit` restoration
    stb-cube.py already applies (see its own long comment for the full
    derivation), needed here too so a user who opens this file in VESTA/
    VMD/Multiwfn (via --keep-cube) sees a correctly-scaled density, not
    one ~6.75x too large. Mutates `grid.grid` in place (division creates
    and rebinds a new array, so this must run on a copy if the caller
    still needs the original e/Ang^3 values afterward -- see
    build_bader_inputs, called BEFORE this in every caller here, which
    is why that ordering matters).
    """
    grid.grid = grid.grid / sile.grid_unit
    grid.write(cube_path)


def read_spin_density(file_rho, geometry, physical_idx, cube_path):
    """Returns a PyBader-ready NET SPIN (magnetization) density array for a
    spin-polarized .RHO file, or None for a non-spin-polarized run OR if
    anything about processing the spin component fails -- this is an
    optional add-on to the main charge analysis, so any failure here just
    degrades to "no spin reported" rather than aborting the whole run.

    Uses sisl's own `index='z'` combination (== up - down) -- NOT
    `index=1` directly. **Real, verified bug fixed here**: for a
    genuinely spin-polarized (nspin=2) SIESTA .RHO, sisl's raw
    `read_grid(index=1)` returns the DOWN-spin channel alone, not a
    "spin density" -- confirmed against a real spin-polarized O2
    calculation (test/6-utils/3-cube/o2.RHO, a textbook triplet ground
    state where SIESTA's own log reports `|S| = 2.0`): the old
    `index=1` reading integrates to 5.0 e (the down channel alone)
    instead of the correct 2.0 e net moment. `index='z'` (sisl's own
    up-minus-down combination) is the same convention already used
    correctly by stb-cube's own `SPIN_INDEX` and by stb-density's
    analogous fix. `index='z'` raises (not `index=1`'s silent wrong
    answer) on a non-spin-polarized file (only 1 component, sisl needs 2
    to form the [1, -1] combination), which is how the two cases are
    told apart -- there is no separate flag/header field for it.

    Returns the in-memory `grid.grid` array directly (see
    build_bader_inputs' own docstring for why this -- not a round-trip
    through pybader's cube reader -- is the correct, bug-free way to get
    a PyBader-ready density array); still writes `cube_path` (correctly
    unit-restored, see write_native_cube) so --keep-cube keeps an
    inspectable spin-density .cube alongside the charge one.
    """
    try:
        spin_sile = sisl.get_sile(file_rho)
        spin_grid = spin_sile.read_grid(index='z')
    except Exception:
        return None
    try:
        spin_density = np.array(spin_grid.grid, dtype=np.float64)
        spin_grid.set_geometry(geometry)
        write_native_cube(spin_grid, spin_sile, cube_path)
        return spin_density
    except Exception as e:
        print(color_text(
            f"   [WARN] Spin-density grid detected but failed to process ({e}) -- "
            "continuing with charge-only analysis.", 'yellow'))
        if os.path.exists(cube_path):
            os.remove(cube_path)
        return None


def find_symmetry_groups(geometry, physical_idx, symprec=0.01, angle_tolerance=5.0):
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

def compute_bader_charges(label, output_dir, speed_mode='normal', ref_file=None,
                          threads=None, vacuum_tol=None, keep_cube=True, export_volumes=False):
    """Reads <label>.RHO + <label>.XV/.fdf, converts to a .cube file via sisl,
    and partitions it into atomic (Bader) basins with PyBader. Pure
    computation plus the step-by-step progress narration (SISL/PyBader are
    slow, I/O-bound steps -- printing live as they run is useful feedback,
    same convention as other tools' "[INFO] Reading structure file..."
    lines); the actual formatted report is built separately by
    print_bader_report() from the results dict returned here. Aborts (exit 1)
    on any unrecoverable error, same fail-fast behavior as before.
    """
    file_rho = f"{label}.RHO"
    file_cube = os.path.join(output_dir, f"{label}.cube")
    file_spin_cube = os.path.join(output_dir, f"{label}_spin.cube")

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
    fallback_syms = []
    dummy_count = 0
    try:
        # --- STEP 1: SISL ---
        try:
            print(f"1. [SISL] Reading geometry and charge density...")

            geometry, _ = read_sisl_geometry_xv_or_fdf(label)
            if geometry is None:
                print(color_text("[ERROR] No geometry file (.XV or .fdf) found.", 'red'))
                sys.exit(1)

            # Atoms with no real element assigned (Z<=0 -- e.g. a floating dummy site with
            # no basis) are excluded from the Bader-ready atom list -- a dummy site has no
            # physical charge for PyBader to partition. This does NOT affect SIESTA's own
            # ghost/BSSE atoms: sisl represents those as `AtomGhost`, whose `.Z` is the real
            # (positive) element number -- verified empirically against this installed sisl.
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

            # index='total' (sisl's own up+down combination for a spin
            # -polarized .RHO, or just the single component for a
            # non-polarized one) -- NOT the default read_grid() (index=0),
            # which for a genuinely spin-polarized file is only the raw
            # UP-spin channel, not the total charge. Same fix/verification
            # as read_spin_density's index='z' above and stb-density's
            # analogous one: on a real spin-polarized O2 run, the old
            # index=0 reading integrated to 7.0 e (up channel alone)
            # instead of the correct 12.0 e (2 O atoms x 6 valence e- each).
            sile = sisl.get_sile(file_rho)
            rho_grid = sile.read_grid(index='total')
            rho_cell = rho_grid.cell.copy()

            # Built directly from the in-memory grid/geometry -- no round-trip through a
            # written .cube file's text (see build_bader_inputs' own docstring for the 2
            # real, verified bugs -- a crash and a silent ~6.75x charge overestimate --
            # this sidesteps entirely). The on-disk .cube (below) is still written, purely
            # for --keep-cube's user-facing "open this in VESTA/VMD" convenience.
            density, lattice, cube_atoms, file_info = build_bader_inputs(
                rho_grid, geometry, physical_idx, f"{label}.cube", output_dir)

            rho_grid.set_geometry(geometry)
            write_native_cube(rho_grid, sile, file_cube)
            cube_files.append(file_cube)

            spin_density = read_spin_density(file_rho, geometry, physical_idx, file_spin_cube)
            has_spin_grid = False
            if spin_density is not None and spin_density.shape == density['charge'].shape:
                density['spin'] = spin_density
                cube_files.append(file_spin_cube)
                has_spin_grid = True
                print(f"   {color_text('[INFO]', 'cyan')} Spin-polarized grid detected -- "
                      "reporting per-atom net spin (magnetic moment) too.")
            elif spin_density is not None:
                print(color_text(
                    "   [WARN] Spin-density grid shape doesn't match the charge grid -- "
                    "skipping spin-resolved analysis.", 'yellow'))

        except Exception as e:
            print(color_text(f"[ERROR] SISL processing failed: {e}", 'red'))
            sys.exit(1)

        # --- Cross-check: the .RHO file's OWN lattice (captured above, before
        # set_geometry() overwrote it) against the geometry's -- catches a genuinely wrong
        # pairing (e.g. a different relaxation step with a different cell). The former atom
        # -count/species checks here (comparing the geometry against a written-then-reread
        # cube file) are gone along with the round-trip they guarded against -- cube_atoms/
        # file_info['elements'] are now built directly FROM physical_idx (build_bader_inputs),
        # so they cannot diverge from it. This cell check alone (nor any check) still can't
        # catch a same-cell, different-atomic-position mismatch -- a .RHO grid carries no
        # independent atomic-position record to compare against; not detectable from these
        # files alone.
        if not np.allclose(rho_cell, geometry.cell, atol=1e-3):
            print(color_text(
                "[ERROR] Lattice mismatch between the .RHO grid and the geometry file -- "
                "they likely don't belong to the same calculation.", 'red'))
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
        n_exported = 0
        if export_volumes:
            # 'atoms' (not 'volumes') so this stays valid even under --speed fast, which
            # deletes bader_volumes but keeps atoms_volumes (see PyBader's own __call__).
            # [-2] is PyBader's own sentinel for "every atom, plus the vacuum bucket if
            # vacuum_tol is set". 'prefix' routes the exported Bader-atoms-<N>.cube files
            # into --output-dir instead of always landing in the current directory --
            # PyBader's own writer just string-concatenates prefix + filename, so a
            # trailing separator is required.
            extra_kwargs['export_mode'] = ('atoms', [-2])
            extra_kwargs['prefix'] = os.path.join(output_dir, '')
            n_exported = len(cube_atoms) + (1 if vacuum_tol else 0)

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
                  f"Bader volume(s) as 'Bader-atoms-<N>.cube' in {output_dir}.")

        # --- STEP 3: Analysis ---
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

            total_final = sum(d['pop_raw'] * correction_factor for d in atoms_data)

            by_species = {}
            for data in atoms_data:
                if data['z_val'] is None:
                    continue
                net = data['z_val'] - data['pop_raw'] * correction_factor
                by_species.setdefault(data['sym'], []).append(net)

            # Species-wide large-and-uniform-charge check -- see
            # SPECIES_CHARGE_MAGNITUDE_TOL/SPECIES_CHARGE_UNIFORMITY_TOL's own
            # comment above for why this is a distinct signature from both
            # suspicious_zero_ids (one atom, ~zero) and inconsistent_groups
            # (symmetry-equivalent atoms that disagree). Needs >=2 atoms of
            # the species to even demonstrate "uniform" -- a single-atom
            # species can't, that's the other two checks' job.
            suspicious_species = []
            for sym, charges in by_species.items():
                if len(charges) < 2:
                    continue
                mean_abs = statistics.mean(abs(c) for c in charges)
                std = statistics.pstdev(charges)
                if mean_abs > SPECIES_CHARGE_MAGNITUDE_TOL and std < SPECIES_CHARGE_UNIFORMITY_TOL:
                    suspicious_species.append((sym, mean_abs, std, len(charges)))

            if suspicious_species:
                for sym, mean_abs, std, n in suspicious_species:
                    print(color_text(
                        f"   [WARN] Species '{sym}' shows a large ({mean_abs:.2f} e-) and nearly "
                        f"uniform (std {std:.3f} e-) net charge across all {n} atoms -- a common "
                        "signature of a pseudopotential with no resolvable valence density "
                        "maximum at that nucleus (e.g. no NLCC, deep frozen semicore states), "
                        "which makes Bader donate that whole species' basin volume to its "
                        "neighbors. Cross-check with a non-topological charge-partitioning "
                        "method before trusting this species' charge.", 'red'))

        except Exception as e:
            print(color_text(f"[ERROR] Analysis failed: {e}", 'red'))
            sys.exit(1)
    finally:
        if not keep_cube:
            for f in cube_files:
                if os.path.exists(f):
                    os.remove(f)

    return {
        'label': label, 'output_dir': output_dir, 'speed_mode': speed_mode, 'ref_file': ref_file,
        'keep_cube': keep_cube, 'export_volumes': export_volumes, 'n_exported': n_exported,
        'source_name': source_name, 'detected_valence': detected_valence,
        'fallback_syms': fallback_syms, 'dummy_count': dummy_count,
        'method': bader_job.method, 'refine_method': bader_job.refine_method,
        'threads': bader_job.threads, 'vacuum_tol': bader_job.vacuum_tol,
        'has_spin': raw_spins is not None, 'has_volume': raw_volumes is not None,
        'atoms_data': atoms_data, 'unknown_syms': unknown_syms,
        'suspicious_zero_ids': suspicious_zero_ids, 'inconsistent_groups': inconsistent_groups,
        'suspicious_species': suspicious_species,
        'total_theory': total_theory, 'total_final': total_final,
        'correction_factor': correction_factor, 'by_species': by_species,
        'cube_files': cube_files if keep_cube else [],
    }


# --- Report formatting --------------------------------------------------
# Same numbered-section report style as the rest of the suite's newer tools
# ([0] RUN METADATA ... [N] SUMMARY & FILES via print_section/print_dual/
# print_table) -- replaces the old ad hoc out_lines/print() duplication.
# The underlying numbers/physics (compute_bader_charges, above) are
# unchanged; only how they're printed -- and several diagnostics that used
# to be console-only (Z_val source detail, dummy-atom exclusion, unit
# -correction factor) are now also persisted to the report file itself.

def _fmt(value, prec=4):
    return f"{value:.{prec}f}" if value is not None else "N/A"


def print_bader_report(results, args, report_path, f_out):
    print_dual(color_text("===== STB-BADER REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time      : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Label          : {results['label']}", f_out)
    print_dual(f"Output dir     : {results['output_dir']}", f_out)
    print_dual(f"Speed mode     : {results['speed_mode']}"
               + (" (CAUTION: trades basin-boundary accuracy for speed)"
                  if results['speed_mode'] == 'fast' else ""), f_out)
    print_dual(f"Reference file : {results['ref_file'] if results['ref_file'] else results['label']+'.out (default)'}", f_out)
    print_dual(f"Vacuum tol.    : {results['vacuum_tol'] if results['vacuum_tol'] else 'disabled'}", f_out)
    print_dual(f"Threads        : {results['threads']}", f_out)
    print_dual(f"Keep .cube     : {'yes' if results['keep_cube'] else 'no'}", f_out)
    print_dual(f"Export volumes : {'yes (' + str(results['n_exported']) + ' file(s))' if results['export_volumes'] else 'no'}", f_out)

    print_section("[1] VALENCE (Z_val) SETUP", f_out)
    print_dual(f"Source: {results['source_name']}", f_out)
    if results['dummy_count']:
        print_dual(f"{results['dummy_count']} dummy/no-element site(s) (Z<=0) excluded from "
                   "the cube file and this analysis.", f_out)
    species_zval = sorted({d['sym']: d['z_val'] for d in results['atoms_data']}.items())
    rows = []
    for sym, z_val in species_zval:
        if z_val is None:
            rows.append(([sym, "N/A", "unknown"], 'red'))
        elif sym in results['detected_valence']:
            rows.append(([sym, f"{z_val:.2f}", "detected (.out)"], None))
        elif sym in results['fallback_syms']:
            rows.append(([sym, f"{z_val:.2f}", "hardcoded fallback"], 'yellow'))
        else:
            rows.append(([sym, f"{z_val:.2f}", "hardcoded default"], None))
    print_table(["Species", "Z_val", "Source"], rows, f_out)
    if results['unknown_syms']:
        print_dual(color_text(
            f"[WARNING] Element(s) {', '.join(sorted(results['unknown_syms']))} not found in "
            "the Z_val dictionary -- their net charge cannot be computed; excluded from the "
            "correction-factor/total below. Use --ref to point at a matching .out file.", 'red'), f_out)

    print_section("[2] PYBADER CONFIGURATION", f_out)
    print_table(["Setting", "Value"], [
        (["Method", results['method']], None),
        (["Refine method", results['refine_method']], None),
        (["Threads", str(results['threads'])], None),
        (["Vacuum tolerance", str(results['vacuum_tol']) if results['vacuum_tol'] else "disabled"], None),
        (["Spin-polarized", "yes" if results['has_spin'] else "no"], None),
    ], f_out)

    print_section("[3] PER-ATOM BADER POPULATIONS", f_out)
    headers = ["Idx", "Elem", "Pop(e-)", "Z_val", "Net Charge", "State"]
    if results['has_volume']:
        headers += ["Vol(Å³)", "SurfD(Å)"]
    if results['has_spin']:
        headers += ["Spin(µB)"]
    inconsistent_ids = {i for ids, _ in results['inconsistent_groups'] for i in ids}
    suspicious_syms = {sym for sym, _, _, _ in results['suspicious_species']}
    rows = []
    for data in results['atoms_data']:
        pop = data['pop_raw'] * results['correction_factor']
        if data['z_val'] is None:
            z_str, net_str, state = "N/A", "N/A", "Unknown Z_val"
        else:
            net = data['z_val'] - pop
            z_str, net_str = f"{data['z_val']:.2f}", f"{net:+.4f}"
            state = "Donor (+)" if net > 0.05 else "Acceptor (-)" if net < -0.05 else "Neutral"
        cells = [str(data['id']), data['sym'], f"{pop:.4f}", z_str, net_str, state]
        if results['has_volume']:
            cells += [_fmt(data['volume']), _fmt(data['surf_dist'])]
        if results['has_spin']:
            cells.append(f"{data['spin']:+.4f}" if data['spin'] is not None else "N/A")
        # red = this atom itself is individually anomalous (near-zero population);
        # yellow = this atom's symmetry group disagrees, OR its whole species reads
        # large-and-uniform -- the atom itself isn't the anomaly, its species is.
        color = 'red' if data['id'] in results['suspicious_zero_ids'] else \
                'yellow' if data['id'] in inconsistent_ids or data['sym'] in suspicious_syms else None
        rows.append((cells, color))
    print_table(headers, rows, f_out)
    theory_note = "" if not results['unknown_syms'] else \
        f" (excludes {len(results['unknown_syms'])} unknown element(s): " \
        f"{', '.join(sorted(results['unknown_syms']))})"
    print_dual(f"Total Integrated: {results['total_final']:.4f} "
               f"(Target: {results['total_theory']:.2f}{theory_note})", f_out)

    print_section("[4] PER-SPECIES SUMMARY", f_out)
    if results['by_species']:
        rows = []
        for sym in sorted(results['by_species']):
            values = results['by_species'][sym]
            mean = statistics.mean(values)
            std = statistics.pstdev(values) if len(values) > 1 else 0.0
            rows.append(([sym, str(len(values)), f"{mean:+.4f}", f"{std:.4f}"],
                         'yellow' if sym in suspicious_syms else None))
        print_table(["Elem", "N", "Mean(e-)", "Std(e-)"], rows, f_out)
    else:
        print_dual("No species with a known Z_val -- nothing to summarize.", f_out)

    print_section("[5] DIAGNOSTICS & WARNINGS", f_out)
    if results['correction_factor'] != 1.0:
        print_dual(f"[INFO] Unit mismatch corrected. Factor: {results['correction_factor']:.4f} "
                   "(applied uniformly to every atom -- a single global scalar; if the real "
                   "deviation is spatially localized, this smears it evenly across every atom's "
                   "reported charge instead of fixing it exactly).", f_out)
    if results['suspicious_zero_ids']:
        ids_str = ', '.join(str(i) for i in results['suspicious_zero_ids'])
        print_dual(color_text(
            f"[WARNING] Atom(s) {ids_str} got essentially zero Bader population "
            f"(< {ZERO_POPULATION_TOL} e-) -- treat with suspicion (common causes: a "
            "pseudopotential with no resolvable density near the nucleus, or the atom's "
            "region folded into a neighbor's basin).", 'red'), f_out)
    if results['inconsistent_groups']:
        for ids, charges in results['inconsistent_groups']:
            pairs = ', '.join(f"#{i}={c:+.3f}" for i, c in zip(ids, charges))
            print_dual(color_text(
                f"[WARNING] Symmetry-equivalent atoms disagree on net charge: {pairs} "
                f"(spread > {SYMMETRY_CHARGE_TOL} e-) -- these sites should be identical by "
                "symmetry; investigate before trusting either value.", 'red'), f_out)
    if results['suspicious_species']:
        for sym, mean_abs, std, n in results['suspicious_species']:
            print_dual(color_text(
                f"[WARNING] Species '{sym}' shows a large ({mean_abs:.2f} e-) and nearly "
                f"uniform (std {std:.3f} e-) net charge across all {n} atoms -- a common "
                "signature of a pseudopotential with no resolvable valence density maximum "
                "at that nucleus (e.g. no NLCC, deep frozen semicore states), which makes "
                "Bader donate that whole species' basin volume to its neighbors. Cross-check "
                "with a non-topological charge-partitioning method before trusting this "
                "species' charge.", 'red'), f_out)
    print_dual(
        "Bader analysis has known limitations this tool cannot detect or correct for: "
        "non-nuclear attractors (basins not centered on any atom, common in ionic/metallic "
        "systems) are always folded into the nearest atom by PyBader; Z_val is only as "
        "accurate as the .out parse or the FALLBACK_VALENCE table (semicore-inclusive "
        "pseudopotentials can differ from the 'standard valence' assumed there); and "
        "--speed fast trades basin-boundary accuracy for speed. See the example README "
        "for the underlying theory and a fuller discussion of each limitation.", f_out)

    print_section("[6] REFERENCES", f_out)
    bib_entries = [citations.SIESTA, citations.SIESTA_RECENT]
    citations.write_bib_file(os.path.join(args.output_dir, BIB_FILE), bib_entries)
    print_dual(color_text(
        f"[OK] Citations for the methods used in this run written to "
        f"'{os.path.join(args.output_dir, BIB_FILE)}' ({len(bib_entries)} entries).", 'green'), f_out)
    print_dual("Bader charge partitioning via PyBader, implementing the grid-based algorithm "
               "of Tang, Sanville & Henkelman (2009); underlying theory: Bader's Atoms-in-"
               "Molecules (AIM) quantum theory (1990) -- see the example README for both.", f_out)

    print_section("[7] SUMMARY & FILES", f_out)
    print_dual("Status         : OK", f_out)
    if results['cube_files']:
        print_dual(f"Cube file(s)   : {', '.join(results['cube_files'])}", f_out)
    else:
        print_dual("Cube file(s)   : deleted (--no-cube)", f_out)
    if results['export_volumes']:
        print_dual(f"Exported vols. : {results['n_exported']} Bader-atoms-<N>.cube file(s) in "
                   f"{results['output_dir']}", f_out)
    print_dual(f"References     : {os.path.join(args.output_dir, BIB_FILE)}", f_out)
    if report_path:
        print_dual(f"Report         : {report_path}", f_out)


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
               "  %(prog)s --label siesta --ref relax/siesta.out --save-report\n"
               "  %(prog)s --label siesta --speed fast --threads 8\n"
               "  %(prog)s --label slab --vacuum-tol 1e-3\n"
               "  %(prog)s --label siesta --export-volumes\n"
    )
    parser.add_argument("-l", "--label", required=True, help="SystemLabel used in Siesta")
    parser.add_argument("-o", "--output-dir", type=str, default=".",
                        help="Directory to write the .cube file(s)/references.bib into (and "
                             f"{REPORT_FILE}, with --save-report; also any --export-volumes "
                             "files) (default: current directory). Created if it doesn't exist.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the full run report to {REPORT_FILE}. Off by default.")
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
                              "is set) inside --output-dir, for visual inspection in VESTA/VMD "
                              "-- one file per atom, each as large as the main grid, so this is "
                              "opt-in (default: off)")
    parser.add_argument("-v", "--version", action="version",
                        version=f"stb-bader {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.threads is not None and args.threads < 1:
        parser.error("--threads must be at least 1.")

    if args.intro:
        show_intro(["Siesta ToolBox Suite - Bader Analysis", f"Version {VERSION} | University of Brasilia"])

    os.makedirs(args.output_dir, exist_ok=True)

    # Remove a stale <label>_BADER.txt left by an older version of this tool
    # -- it always wrote this unconditionally to the current directory (there
    # was no --output-dir concept before this migration), so that's where a
    # leftover one would be, regardless of --output-dir now.
    stale_path = f"{args.label}_BADER.txt"
    if os.path.exists(stale_path):
        os.remove(stale_path)

    results = compute_bader_charges(
        args.label, args.output_dir, args.speed, args.ref,
        threads=args.threads, vacuum_tol=args.vacuum_tol, keep_cube=args.keep_cube,
        export_volumes=args.export_volumes)

    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    print_bader_report(results, args, report_path, f_out)

    if f_out:
        f_out.close()

    print("\n" + "-" * 60)
    print(color_text("Electron counting is like accounting, but the currency is negative.\n", 'bold'))

if __name__ == "__main__":
    main()
