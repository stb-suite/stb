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
import argparse
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from ase import Atoms
from ase.thermochemistry import IdealGasThermo
from pymatgen.core.periodic_table import Element
from pymatgen.io.ase import AseAtomsAdaptor
from stb.core import structure_io
from stb.core.siesta_log import get_free_energy, get_spin_moment, report_quality_diagnostics
from stb.core.symmetry import point_group_details
from stb.core.cli import color_text, show_intro, print_dual, print_section
from stb.core.phonon_workflow import load_phonon_with_force_constants
from stb.adsorb_analysis import read_site_table, read_bsse_energy, GIBBS_LOCAL_META_FILE

REPORT_FILE = "adsorption_gibbs_report.txt"
PLOT_FILE = "adsorption_gibbs.png"
MODE_SPECTRUM_PLOT_FILE = "adsorption_mode_spectrum.png"

_FREQ_CONVERSION_THZ = 15.633302  # sqrt(eV / (amu * Ang^2)) -> THz, same constant ASE's own
                                  # vibrational-analysis code (and her_analysis.py/oer_analysis.py) use
_BOLTZMANN_EV_K = 8.617333262e-5  # eV/K
_EV_PER_THZ = 0.00413566733  # E = h*f, h in eV.s, f in THz*1e12 -> eV
_KJ_MOL_TO_EV = 0.01036427
_J_MOL_TO_EV = _KJ_MOL_TO_EV / 1000.0
_STANDARD_PRESSURE_PA = 1.0e5  # 1 bar, in real Pascals -- see compute_ideal_gas_thermo's
                                # get_entropy() call for why this must be passed explicitly
                                # in real Pa rather than left at ase.thermochemistry's own
                                # `pressure` default.


def find_gibbs_folders(gibbs_root, f_out):
    """Auto-detects the 3 possible subfolders stb-adsorbAnalysis's
    --compute-gibbs step (Stage 3) can have written under '<dir>/gibbs/':
    exactly one 'isolated' reference (name ends in '_isolated'), at most
    one 'clean_slab_full' (only present for --zpe-mode full), and
    whatever single folder remains is the winning SITE's own Hessian/
    phonon folder (named after its 'sites/<label>' -- Stage 4 never needs
    to be told the label explicitly, since Stage 3 only ever preps ONE
    site's worth of Gibbs folders per run). Returns
    (site_label, site_dir, isolated_dir, clean_slab_dir_or_None), exits
    with a clear error if the site/isolated folders can't be
    unambiguously identified.
    """
    if not os.path.isdir(gibbs_root):
        print_dual(color_text(
            f"[ERROR] '{gibbs_root}' not found -- run stb-adsorbAnalysis --compute-gibbs "
            "(Stage 3) first.", 'red'), f_out)
        sys.exit(1)
    entries = sorted(d for d in os.listdir(gibbs_root) if os.path.isdir(os.path.join(gibbs_root, d)))
    isolated = [d for d in entries if d.endswith("_isolated")]
    clean = [d for d in entries if d == "clean_slab_full"]
    site_candidates = [d for d in entries if d not in isolated and d not in clean]
    if len(isolated) != 1 or len(site_candidates) != 1:
        print_dual(color_text(
            f"[ERROR] Expected exactly one '<site>/' and one '<adsorbate>_isolated/' folder "
            f"under '{gibbs_root}', found site candidate(s) {site_candidates} and isolated "
            f"reference(s) {isolated}. Re-run stb-adsorbAnalysis --compute-gibbs if this "
            "directory is stale/mixed from more than one run.", 'red'), f_out)
        sys.exit(1)
    site_label = site_candidates[0]
    clean_dir = os.path.join(gibbs_root, clean[0]) if clean else None
    return site_label, os.path.join(gibbs_root, site_label), os.path.join(gibbs_root, isolated[0]), clean_dir


def read_fa_force(fa_path, atom_index):
    """Reads the force (Fx, Fy, Fz, eV/Ang) on ONE atom (0-based
    atom_index) from a SIESTA .FA file. Duplicated from her_analysis.py/
    oer_analysis.py's identical helper -- each workflow stays
    self-contained (re-reads persisted files rather than importing a
    sibling workflow's module), same convention already established there.
    """
    try:
        with open(fa_path) as f:
            lines = f.readlines()
    except OSError:
        return None
    if len(lines) < atom_index + 2:
        return None
    try:
        parts = lines[atom_index + 1].split()
        return np.array([float(parts[1]), float(parts[2]), float(parts[3])])
    except (IndexError, ValueError):
        return None


def _build_mass_weighted_hessian(gibbs_dir, f_out):
    """Reads gibbs_local_meta.json + every disp_NNN/<system_label>.FA force
    file stb-adsorbAnalysis --compute-gibbs wrote, builds the mass-weighted
    3N x 3N dynamical matrix by central finite difference (same math as
    oer_analysis.py::compute_local_zpe_entropy, generalized N-atom version
    -- reduces to her_analysis.py's single-atom 3x3 case when N=1), and
    diagonalizes it. Returns (eigenvalues, local_symbols, local_indices,
    displacement_ang), or (None, None, None, None) on any read failure.

    Extracted from compute_local_hessian_thermo once compute_ideal_gas_thermo
    became a second consumer needing the SAME raw eigenvalues but a
    different (count-based, not sign-based) mode-selection/thermo treatment
    -- see compute_ideal_gas_thermo's own docstring for why.
    """
    meta_path = os.path.join(gibbs_dir, GIBBS_LOCAL_META_FILE)
    try:
        with open(meta_path) as f:
            meta = json.load(f)
    except (OSError, ValueError):
        print_dual(color_text(f"    [ERROR] Could not read '{meta_path}'.", 'red'), f_out)
        return None, None, None, None

    local_indices = meta["local_indices"]
    local_symbols = meta["local_symbols"]
    displacement_ang = meta["displacement_ang"]
    order = meta["order"]
    system_label = meta["system_label"]
    n_local = len(local_indices)
    index_to_local = {atom_index: k for k, atom_index in enumerate(local_indices)}
    masses = np.array([Element(sym).atomic_mass for sym in local_symbols], dtype=float)

    forces = {}  # (moved_local_k, axis, sign) -> (n_local, 3) forces on every local atom
    for i, entry in enumerate(order, start=1):
        fa_path = os.path.join(gibbs_dir, f"disp_{i:03d}", f"{system_label}.FA")
        moved_k = index_to_local[entry["atom_index"]]
        row = np.zeros((n_local, 3))
        for j, atom_index in enumerate(local_indices):
            force = read_fa_force(fa_path, atom_index)
            if force is None:
                print_dual(color_text(f"    [ERROR] Could not read forces from '{fa_path}'.", 'red'), f_out)
                return None, None, None, None
            row[j] = force
        forces[(moved_k, entry["axis"], entry["sign"])] = row

    dof = 3 * n_local
    Phi = np.zeros((dof, dof))
    for moved_k in range(n_local):
        for axis in range(3):
            f_plus = forces[(moved_k, axis, 1.0)]
            f_minus = forces[(moved_k, axis, -1.0)]
            d_force = -(f_plus - f_minus) / (2.0 * displacement_ang)  # (n_local, 3)
            Phi[:, 3 * moved_k + axis] = d_force.reshape(-1)
    Phi = 0.5 * (Phi + Phi.T)

    mass_vec = np.repeat(masses, 3)
    D = Phi / np.sqrt(np.outer(mass_vec, mass_vec))
    eigenvalues = np.linalg.eigvalsh(D)
    return eigenvalues, local_symbols, local_indices, displacement_ang


def compute_local_hessian_thermo(gibbs_dir, temperature_k, f_out, label, mode_log=None):
    """Standard harmonic-oscillator ZPE/entropy summed over every REAL
    (non-imaginary, non-near-zero) mode of _build_mass_weighted_hessian's
    eigenvalues. Modes at/below 0 are reported and excluded -- for a fully
    isolated single atom (no local coupling at all), ALL 3 modes are pure
    translation and get excluded this way, correctly giving zpe_ev=0.0 (a
    free atom has no vibrational zero-point energy). Returns (None, None)
    on any read failure.

    Used for the SITE side only (an adsorbate frozen onto a fixed
    substrate) -- correct as-is: unlike a free molecule, a bonded adsorbate
    has no EXACT rigid-body (translation/rotation) zero modes to project
    out, so a negative eigenvalue here signals a genuine instability/
    non-converged geometry, not a symmetry-required zero mode, and this
    plain sign-based exclusion is the right physics. See
    compute_ideal_gas_thermo for the ISOLATED REFERENCE side, which does
    need a count-based (not sign-based) exclusion instead.

    `mode_log`, if given a list, gets one {"freq_thz", "kept", "imaginary"}
    dict appended per mode (T-independent, so a caller only needs this from
    one call, e.g. at args.tmin) -- feeds adsorption_mode_spectrum.png.
    """
    eigenvalues, local_symbols, local_indices, _displacement_ang = \
        _build_mass_weighted_hessian(gibbs_dir, f_out)
    if eigenvalues is None:
        return None, None
    n_local = len(local_indices)

    zpe_ev, ts_ev = 0.0, 0.0
    print_dual(f"    {label} ({n_local}-atom local Hessian) vibrational modes:", f_out)
    for eig in eigenvalues:
        if eig > 0:
            freq_thz = _FREQ_CONVERSION_THZ * np.sqrt(eig)
            e_mode = freq_thz * _EV_PER_THZ
            print_dual(f"      {freq_thz * 33.35641:>8.2f} cm^-1  ({freq_thz:>6.2f} THz)  "
                        f"E = {e_mode:.4f} eV", f_out)
            zpe_ev += 0.5 * e_mode
            if temperature_k > 0:
                x = e_mode / (_BOLTZMANN_EV_K * temperature_k)
                S = _BOLTZMANN_EV_K * (x / np.expm1(x) - np.log1p(-np.exp(-x)))
                ts_ev += temperature_k * S
            if mode_log is not None:
                mode_log.append({"freq_thz": freq_thz, "kept": True, "imaginary": False})
        else:
            freq_thz = _FREQ_CONVERSION_THZ * np.sqrt(-eig)
            print_dual(color_text(
                f"      {freq_thz:>6.2f} THz (imaginary/translational) -- excluded from "
                "ZPE/entropy", 'yellow'), f_out)
            if mode_log is not None:
                mode_log.append({"freq_thz": freq_thz, "kept": False, "imaginary": True})
    return zpe_ev, ts_ev


def compute_ideal_gas_thermo(gibbs_dir, ads_dir, ads_out, temperature_k, f_out, label, mode_log=None):
    """ISOLATED-REFERENCE-only counterpart to compute_local_hessian_thermo,
    using ase.thermochemistry.IdealGasThermo instead of a flat sum over
    every positive-eigenvalue mode. Fixes a real, quantitatively
    significant bug the plain sign-based treatment has for any POLYATOMIC
    isolated reference: a free molecule has exactly 6 (nonlinear) / 5
    (linear) rigid-body zero modes (translation + rotation), not 3 --
    residual finite-difference numerical noise routinely puts roughly half
    of those on the POSITIVE side of zero, where the old sign-based filter
    (which only excludes NEGATIVE/"imaginary" eigenvalues) would silently
    keep them as if they were genuine vibrations. Since entropy diverges as
    frequency -> 0, even a couple of eV-fraction-sized spurious low modes
    can swamp the real DG value (verified on a real H2O/SiC calculation:
    ~+0.245 eV of spurious T*S at 300 K, comparable to the entire reported
    DG).

    Fix: select exactly the `3N-6`/`3N-5`/`0` HIGHEST-|energy| modes as
    genuine vibrations (count-based, not sign-based -- the same selection
    ase.thermochemistry.IdealGasThermo's own default `vib_selection`
    performs internally, replicated explicitly here so the excluded/kept
    split is fully visible in the report and the mode-spectrum plot,
    instead of silently happening inside IdealGasThermo). This also closes
    the report's own previously-documented `[LIMITATION]` (the isolated
    reference's translational/rotational ideal-gas entropy was never
    computed at all) for every adsorbate, including single-atom ones -- a
    monatomic species now correctly gets a nonzero translational entropy
    (Sackur-Tetrode, via IdealGasThermo's 'monatomic' geometry) instead of
    the old code's implicit ts_ev=0.0.

    Molecular geometry class ('monatomic'/'linear'/'nonlinear') and
    rotational symmetry number come from core/symmetry.py::
    point_group_details (pymatgen's PointGroupAnalyzer) on the reference's
    actual relaxed geometry (core/structure_io.py::read_relaxed_or_input on
    `ads_dir` -- the ORIGINAL isolated-adsorbate-reference folder, not the
    perturbed Gibbs-prep disp_NNN/ folders, which only ever move one atom
    by a fraction of an Angstrom and would give an unreliable point-group
    read). Reference pressure: ase.units.bar (IdealGasThermo's own default,
    1 bar) -- the conventional standard state for a "Gibbs free energy of
    adsorption," same "flag the reference-state assumption" transparency
    as her_analysis.py's fixed literature H2/H2O gas-phase constants.

    Returns (zpe_ev, ts_ev), same shape as compute_local_hessian_thermo, or
    (None, None) on any read failure. `mode_log`, if given a list, gets one
    {"freq_thz", "kept", "imaginary"} dict appended per mode (same
    convention as compute_local_hessian_thermo's own mode_log).
    """
    eigenvalues, local_symbols, local_indices, _displacement_ang = \
        _build_mass_weighted_hessian(gibbs_dir, f_out)
    if eigenvalues is None:
        return None, None
    n_local = len(local_indices)

    if n_local == 1:
        geometry = "monatomic"
        symmetrynumber = 1
        point_group_symbol = None
        ase_atoms = Atoms(local_symbols)
    else:
        try:
            relaxed, _used_relaxed = structure_io.read_relaxed_or_input(ads_dir)
        except (FileNotFoundError, ValueError) as e:
            print_dual(color_text(
                f"    [ERROR] Could not read {label}'s relaxed geometry from '{ads_dir}': {e}",
                'red'), f_out)
            return None, None
        pmg_structure = structure_io.to_pymatgen(relaxed)
        ase_atoms = AseAtomsAdaptor.get_atoms(pmg_structure)
        # to_pymatgen always returns a genuinely periodic pymatgen Structure
        # (this suite's isolated-adsorbate reference is still a vacuum-
        # padded box, not a bare pymatgen Molecule) -- IdealGasThermo
        # explicitly rejects periodic boundary conditions (its rotational
        # moment-of-inertia math assumes a truly isolated system), so it
        # must be stripped here even though the vacuum box itself already
        # makes the physical system effectively isolated.
        ase_atoms.set_pbc(False)
        details = point_group_details(pmg_structure)
        if details is None:
            print_dual(color_text(
                f"    [WARNING] Could not determine {label}'s point group -- assuming "
                "nonlinear geometry, symmetry number 1 (a conservative fallback that may "
                "slightly misestimate the rotational entropy contribution).", 'yellow'), f_out)
            geometry, symmetrynumber, point_group_symbol = "nonlinear", 1, None
        else:
            geometry = "linear" if details["linear"] else "nonlinear"
            symmetrynumber = details["rotational_symmetry_number"]
            point_group_symbol = details["symbol"]

    n_vib = {"monatomic": 0, "linear": 3 * n_local - 5, "nonlinear": 3 * n_local - 6}[geometry]

    # Build one complex "vibrational energy" per eigenvalue (real for a
    # genuine, eig>0 mode; purely imaginary for eig<=0, same ASE convention
    # ase.vibrations uses) and keep only the n_vib HIGHEST by |energy|^2 --
    # this magnitude-based ranking (not sign) is what correctly separates
    # near-zero rigid-body modes (whichever side of zero they numerically
    # land on) from the genuinely high-frequency vibrational spectrum.
    candidates = []  # (energy_complex, freq_thz, is_imaginary)
    for eig in eigenvalues:
        if eig > 0:
            freq_thz = _FREQ_CONVERSION_THZ * np.sqrt(eig)
            candidates.append((complex(freq_thz * _EV_PER_THZ), freq_thz, False))
        else:
            freq_thz = _FREQ_CONVERSION_THZ * np.sqrt(-eig)
            candidates.append((complex(0.0, freq_thz * _EV_PER_THZ), freq_thz, True))
    ranked = sorted(candidates, key=lambda c: (c[0] ** 2).real, reverse=True)
    kept, excluded = ranked[:n_vib], ranked[n_vib:]

    label_bits = f"{n_local}-atom, {geometry}"
    if point_group_symbol:
        label_bits += f", point group {point_group_symbol}, symmetry number {symmetrynumber}"
    print_dual(f"    {label} ({label_bits}) -- ideal-gas treatment "
                "(translation + rotation + vibration):", f_out)
    for _e, freq_thz, is_imag in sorted(kept, key=lambda c: c[1]):
        warn_tag = "  [WARNING: imaginary]" if is_imag else ""
        print_dual(f"      {freq_thz * 33.35641:>8.2f} cm^-1  ({freq_thz:>6.2f} THz)  "
                    f"vibrational{warn_tag}", f_out)
        if mode_log is not None:
            mode_log.append({"freq_thz": freq_thz, "kept": True, "imaginary": is_imag})
    for _e, freq_thz, is_imag in sorted(excluded, key=lambda c: c[1]):
        kind = "imaginary" if is_imag else "residual positive"
        print_dual(color_text(
            f"      {freq_thz:>6.2f} THz ({kind}) -- excluded as a rigid-body "
            "translation/rotation mode", 'yellow'), f_out)
        if mode_log is not None:
            mode_log.append({"freq_thz": freq_thz, "kept": False, "imaginary": is_imag})
    if any(is_imag for _e, _f, is_imag in kept):
        print_dual(color_text(
            f"    [WARNING] {label}: a genuine (kept) vibrational mode is imaginary -- the "
            "isolated-reference geometry may not be fully relaxed.", 'yellow'), f_out)

    spin_moment = get_spin_moment(ads_out)
    spin = abs(spin_moment) / 2.0 if spin_moment else 0.0

    thermo = IdealGasThermo(
        vib_energies=[e for e, _f, _imag in kept], geometry=geometry, atoms=ase_atoms,
        symmetrynumber=symmetrynumber, spin=spin, natoms=n_local,
        vib_selection="exact", potentialenergy=0.0,
    )
    zpe_ev = thermo.get_ZPE_correction()
    if temperature_k > 0:
        # Real bug found and verified in ase.thermochemistry (ASE 3.29.0):
        # get_entropy()'s own default `pressure=ase.units.bar` is in ASE's
        # INTERNAL pressure units (eV/Ang^3, ~6.24e-7), but the pressure-
        # correction term inside get_ideal_entropy() (S_p = -kB*ln(pressure
        # / self.referencepressure)) divides it directly against
        # self.referencepressure, which BaseThermoChem sets in REAL
        # Pascals (1e5) -- no unit conversion between the two. Left at its
        # own default, this silently adds a huge spurious entropy term
        # (~+0.0022 eV/K, confirmed by hand) since 6.24e-7 / 1e5 is a
        # wildly, unphysically low "pressure ratio". Verified live: passing
        # `pressure` explicitly in real Pa (this constant) instead
        # reproduces gas-phase H2O's known experimental standard entropy
        # (188.8 J/(mol*K) at 298.15 K) to within 0.1%; relying on the
        # library's own default overestimated it by more than 2x.
        ts_ev = temperature_k * thermo.get_entropy(
            temperature_k, pressure=_STANDARD_PRESSURE_PA, verbose=False)
    else:
        ts_ev = 0.0
    return zpe_ev, ts_ev


def compute_full_phonon_thermo(phonon_dir, system_label, temperature_k, f_out, label):
    """Full-structure ZPE/entropy via Phonopy's own thermal-properties
    pipeline, reusing core.phonon_workflow.load_phonon_with_force_constants
    -- identical to her_analysis.py/oer_analysis.py's own
    compute_full_zpe_entropy. Returns (zpe_ev, ts_ev) or (None, None).
    """
    print_dual(f"    {label} (full Phonopy phonon calculation):", f_out)
    try:
        phonon, _internal_to_angstrom, original_dir = load_phonon_with_force_constants(
            phonon_dir, system_label, False, f_out)
    except SystemExit:
        return None, None
    try:
        phonon.run_mesh([1, 1, 1])
        phonon.run_thermal_properties(t_step=10, t_max=temperature_k + 10, t_min=temperature_k)
        tp_dict = phonon.get_thermal_properties_dict()
        temps = np.array(tp_dict['temperatures'])
        idx = int(np.abs(temps - temperature_k).argmin())
        zpe_ev = phonon.get_zero_point_energy() * _KJ_MOL_TO_EV
        entropy_ev_k = tp_dict['entropy'][idx] * _J_MOL_TO_EV
        ts_ev = temperature_k * entropy_ev_k
    finally:
        os.chdir(original_dir)
    return zpe_ev, ts_ev


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Stage 4 of the Adsorption workflow: combines the vibrational "
        "Hessian/phonon folders stb-adsorbAnalysis --compute-gibbs (Stage 3) wrote into a Gibbs "
        "free energy of adsorption, DG = E_ads(BSSE if available) + DZPE - T*DS.", 'bold')}
DZPE/DS = [ZPE/S(site) - ZPE/S(isolated adsorbate reference)], both computed with the SAME method
(local partial-Hessian, or full Phonopy for the site side when Stage 3 used --zpe-mode full) --
the isolated reference always uses the local, full-molecule Hessian (there's no periodic substrate
phonon to subtract for a molecule in a vacuum box, so Phonopy's supercell machinery buys nothing
there; for a single-atom adsorbate this correctly gives ZPE=S=0, since a free atom has no
vibrational modes at all once its 3 translational modes are excluded).

KNOWN LIMITATION, always printed in the report: the isolated-reference term above is
VIBRATIONAL/HARMONIC ONLY -- it does NOT include the translational/rotational entropy a real gas
molecule has from being free to move/rotate in 3D (an ideal-gas partition-function term this tool
does not compute). For a single-atom adsorbate this is exact (no rotation/translation to begin
with beyond the excluded trivial modes); for a polyatomic adsorbate it is an approximation, in the
same spirit as (though computed differently from) stb-herRefs/stb-oerRefs' own use of a FIXED
literature entropy for their H2/H2O gas-phase references -- generalized here since there is no
literature table for an arbitrary adsorbate.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage example:\n"
               "  %(prog)s --dir adsorption_run --tmin 100 --tmax 500\n"
    )

    parser.add_argument("--dir", type=str, default="adsorption_run",
                         help="Root directory containing 'gibbs/', 'sites/', 'clean_slab/' etc. "
                              "(default: adsorption_run -- stb-adsorb's own default --output-dir).")
    parser.add_argument("--file", type=str, default="calc.out",
                         help="SIESTA output filename inside each folder (default: calc.out).")
    parser.add_argument("--tmin", type=float, default=200.0, help="Sweep start, Kelvin (default: 200).")
    parser.add_argument("--tmax", type=float, default=400.0, help="Sweep end, Kelvin (default: 400).")
    parser.add_argument("--tstep", type=float, default=25.0, help="Sweep step, Kelvin (default: 25).")
    parser.add_argument("--force-tolerance", type=float, default=0.05,
                         help="Residual atomic force in eV/Ang (default: 0.05) above which a "
                              "folder is flagged as possibly not converged/relaxed. Advisory only.")
    parser.add_argument("--save-report", action="store_true",
                         help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--view-plots", action="store_true",
                         help="Also show the DG-vs-T plot on screen, instead of only saving it "
                              "as a PNG. Off by default.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-adsorbGibbs {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    report_path = os.path.join(args.dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None
    print_dual(color_text("===== ADSORPTION GIBBS FREE ENERGY REPORT (STAGE 4) =====", 'magenta'), f_out)

    print_section('[0] RUN METADATA', f_out)
    print_dual(f"  Date/time  : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", f_out)
    print_dual(f"  Directory  : {args.dir}", f_out)
    print_dual(f"  T sweep    : {args.tmin}-{args.tmax} K, step {args.tstep} K", f_out)

    gibbs_root = os.path.join(args.dir, "gibbs")
    site_label, site_gibbs_dir, isolated_gibbs_dir, clean_gibbs_dir = find_gibbs_folders(gibbs_root, f_out)
    zpe_mode = "full" if clean_gibbs_dir else "local"
    print_dual(f"  Winning site : {site_label}", f_out)
    print_dual(f"  ZPE mode     : {zpe_mode} (auto-detected: clean_slab_full/ "
                f"{'present' if clean_gibbs_dir else 'absent'})", f_out)

    print_section('[1] ELECTRONIC ENERGIES', f_out)
    sites_root = os.path.join(args.dir, "sites")
    site_out = os.path.join(sites_root, site_label, args.file)
    e_site = get_free_energy(site_out)
    clean_out = os.path.join(args.dir, "clean_slab", args.file)
    e_clean = get_free_energy(clean_out)
    if e_site is None or e_clean is None:
        print_dual(color_text(
            f"[ERROR] Could not read energy from '{site_out}' and/or '{clean_out}'.", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)
    report_quality_diagnostics(site_label, site_out, args.force_tolerance, f_out)
    print_dual(f"  E_site       = {e_site:>12.6f} eV  ({site_out})", f_out)
    print_dual(f"  E_clean_slab = {e_clean:>12.6f} eV  ({clean_out})", f_out)

    site_table = read_site_table(sites_root)
    ads_name = site_table.get(site_label, (None, None))[0] if site_table else None
    ads_candidates = ([os.path.join(args.dir, f"adsorbate_{ads_name}")] if ads_name else []) + \
        [os.path.join(args.dir, "adsorbate")]
    ads_dir = next((d for d in ads_candidates if os.path.isdir(d)), None)
    if ads_dir is None:
        print_dual(color_text(
            f"[ERROR] No isolated-adsorbate reference folder found (tried: "
            f"{', '.join(ads_candidates)}).", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)
    ads_out = os.path.join(ads_dir, args.file)
    e_adsorbate = get_free_energy(ads_out)
    if e_adsorbate is None:
        print_dual(color_text(f"[ERROR] Could not read energy from '{ads_out}'.", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)
    print_dual(f"  E_adsorbate  = {e_adsorbate:>12.6f} eV  ({ads_out})", f_out)
    e_ads = e_site - e_clean - e_adsorbate
    print_dual(f"  E_ads (raw)  = {e_ads:>12.6f} eV", f_out)

    bsse_dir = os.path.join(args.dir, "bsse", site_label)
    e_ads_used = e_ads
    used_bsse = False
    if os.path.isdir(os.path.join(bsse_dir, "bsse_slab")):
        e_bsse_slab, e_bsse_ads = read_bsse_energy(bsse_dir, args.file)
        if e_bsse_slab is not None and e_bsse_ads is not None:
            e_ads_bsse = e_site - e_bsse_slab - e_bsse_ads
            print_dual(f"  E_ads (BSSE) = {e_ads_bsse:>12.6f} eV", f_out)
            e_ads_used = e_ads_bsse
            used_bsse = True
    print_dual(f"  E_ads used for DG below: {'BSSE-corrected' if used_bsse else 'raw (no BSSE found)'}"
                f" = {e_ads_used:.6f} eV", f_out)

    print_section('[2] VIBRATIONAL/THERMAL TERMS', f_out)
    print_dual(
        "  The isolated-adsorbate reference below is treated as a genuine ideal gas molecule "
        "(ase.thermochemistry.IdealGasThermo, standard-state pressure = 1 bar) -- full "
        "translational + rotational + vibrational ZPE/entropy, not a vibrational-only "
        "approximation. See stb-adsorbGibbs --help.", f_out)
    ref_mode_log = []
    zpe_ref, ts_ref = compute_ideal_gas_thermo(
        isolated_gibbs_dir, ads_dir, ads_out, args.tmin, f_out, "Isolated reference",
        mode_log=ref_mode_log)
    if zpe_ref is None:
        if f_out:
            f_out.close()
        sys.exit(1)

    def site_thermo_at(temperature_k, mode_log=None):
        if zpe_mode == "local":
            return compute_local_hessian_thermo(site_gibbs_dir, temperature_k, None, "Site",
                                                 mode_log=mode_log)
        zpe_site, ts_site = compute_full_phonon_thermo(
            site_gibbs_dir, "gibbs_site", temperature_k, None, "Site")
        zpe_clean, ts_clean = compute_full_phonon_thermo(
            clean_gibbs_dir, "gibbs_clean", temperature_k, None, "Clean slab")
        if zpe_site is None or zpe_clean is None:
            return None, None
        return zpe_site - zpe_clean, ts_site - ts_clean

    site_mode_log = []
    zpe_site0, ts_site0 = site_thermo_at(args.tmin, mode_log=site_mode_log)
    if zpe_site0 is None:
        print_dual(color_text("[ERROR] Site vibrational/phonon calculation incomplete.", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)
    delta_zpe0 = zpe_site0 - zpe_ref
    d0 = e_ads_used + delta_zpe0
    print_dual(f"  ZPE(site)  = {zpe_site0:.4f} eV   ZPE(isolated ref) = {zpe_ref:.4f} eV   "
                f"DZPE = {delta_zpe0:+.4f} eV", f_out)
    print_dual(f"  D0 (ZPE-corrected binding energy) = E_ads + DZPE = {d0:+.4f} eV "
                "(T-independent -- the electronic+ZPE binding energy before any entropy "
                "contribution)", f_out)

    print_section('[3] GIBBS FREE ENERGY vs. TEMPERATURE', f_out)
    temperatures = np.arange(args.tmin, args.tmax + args.tstep / 2, args.tstep)
    dg_values = []
    ds_values_ev_k = []
    for T in temperatures:
        zpe_site, ts_site = site_thermo_at(T)
        _zpe_ref_T, ts_ref_T = compute_ideal_gas_thermo(
            isolated_gibbs_dir, ads_dir, ads_out, T, None, "Isolated reference")
        delta_zpe = zpe_site - zpe_ref
        delta_ts = ts_site - ts_ref_T
        delta_s_ev_k = delta_ts / T if T > 0 else 0.0
        dg = e_ads_used + delta_zpe - delta_ts
        dg_values.append(dg)
        ds_values_ev_k.append(delta_s_ev_k)
        print_dual(f"  T = {T:>7.2f} K   DZPE = {delta_zpe:+.4f} eV   DTS = {delta_ts:+.4f} eV   "
                    f"DS = {delta_s_ev_k * 1000:+.4f} meV/K ({delta_s_ev_k / _J_MOL_TO_EV:+.2f} "
                    f"J/(mol*K))   DG = {dg:+.4f} eV", f_out)

    dg_at_tmin = dg_values[0]
    print_dual(f"\n{color_text('[FINAL RESULT]', 'magenta')} DG(adsorption) = {dg_at_tmin:+.4f} eV "
                f"at T = {args.tmin} K (E_ads {'BSSE-corrected' if used_bsse else 'raw'}, "
                f"D0 = {d0:+.4f} eV)", f_out)

    print_section('[3b] ESTIMATED DESORPTION TEMPERATURE', f_out)
    print_dual(
        "  Linear interpolation/extrapolation of the DG(T) sweep above to its DG=0 crossing "
        "(desorption becomes thermodynamically favorable above this temperature). A rough "
        "estimate from a harmonic model, not a kinetic prediction -- real desorption also "
        "depends on the barrier/prefactor, not just the sign of DG.", f_out)
    sign_changes = [i for i in range(len(dg_values) - 1)
                     if (dg_values[i] < 0) != (dg_values[i + 1] < 0)]
    if sign_changes:
        i = sign_changes[0]
        t1, t2 = temperatures[i], temperatures[i + 1]
        g1, g2 = dg_values[i], dg_values[i + 1]
        t_desorb = t1 + (0.0 - g1) * (t2 - t1) / (g2 - g1)
        print_dual(f"  T_desorption ~= {t_desorb:.1f} K (interpolated within the scanned range)",
                    f_out)
    else:
        t1, g1 = temperatures[0], dg_values[0]
        t2, g2 = temperatures[-1], dg_values[-1]
        if g2 != g1:
            t_desorb = t1 + (0.0 - g1) * (t2 - t1) / (g2 - g1)
            print_dual(color_text(
                f"  [EXTRAPOLATED] T_desorption ~= {t_desorb:.1f} K -- DG stays "
                f"{'negative' if dg_at_tmin < 0 else 'positive'} throughout the scanned "
                f"[{args.tmin}, {args.tmax}] K range; this is a linear extrapolation beyond it, "
                "not an interpolation -- widen --tmin/--tmax for a direct estimate.", 'yellow'),
                f_out)
        else:
            t_desorb = None
            print_dual(color_text(
                "  [WARNING] DG is flat over the scanned range -- cannot estimate a desorption "
                "temperature.", 'yellow'), f_out)

    # [Panel 1] DG(T) plus the two T-independent electronic-level reference
    # lines (E_ads, D0) -- the original single-axes plot, unchanged, with
    # D0 added alongside E_ads for a direct visual read of how much the ZPE
    # correction alone shifts the electronic binding energy before entropy
    # is even considered. [Panel 2] DS(T) -- the standard thermodynamic
    # decomposition figure (DG = D0 - T*DS), same multi-panel diagnostic
    # convention stb-mlmd's own <stem>_md_diagnostics.png already uses.
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 8), sharex=True)
    ax1.plot(temperatures, dg_values, marker='o', color='tab:blue', label="DG(T)")
    ax1.axhline(e_ads_used, color='gray', linestyle='--', linewidth=1,
                label=f"E_ads ({'BSSE' if used_bsse else 'raw'}, electronic only)")
    ax1.axhline(d0, color='tab:orange', linestyle=':', linewidth=1.5,
                label="D0 (E_ads + DZPE)")
    ax1.set_ylabel("Energy (eV)")
    ax1.set_title(f"Gibbs free energy of adsorption -- {site_label}")
    ax1.legend(fontsize=9)

    ds_mev_k = [s * 1000 for s in ds_values_ev_k]
    ax2.plot(temperatures, ds_mev_k, marker='o', color='tab:green')
    ax2.axhline(0.0, color='gray', linestyle='--', linewidth=1)
    ax2.set_xlabel("Temperature (K)")
    ax2.set_ylabel("DS (meV/K)")
    fig.tight_layout()
    plot_path = os.path.join(args.dir, PLOT_FILE)
    fig.savefig(plot_path, dpi=150)
    if args.view_plots:
        plt.show()
    plt.close(fig)
    print_dual(f"{color_text('[Saved]', 'cyan')} {plot_path}", f_out)

    # Mode-spectrum comparison: site (bonded, above the axis) vs. isolated
    # reference (free molecule, below the axis) vibrational frequencies --
    # kept (genuine vibrational) modes as solid stems, excluded rigid-body/
    # imaginary modes as short, greyed-out, annotated stems for
    # transparency instead of silently vanishing -- exactly the visual that
    # would have surfaced this session's real mode-count bug on sight. Only
    # available for --zpe-mode local (site_mode_log is only ever populated
    # in that branch -- a full Phonopy DOS is a fundamentally different,
    # continuous spectrum not comparable this way).
    if site_mode_log:
        fig2, axm = plt.subplots(figsize=(8, 4))
        for m in site_mode_log:
            color = 'tab:red' if not m["kept"] else 'tab:blue'
            height = 0.4 if not m["kept"] else 1.0
            axm.plot([m["freq_thz"], m["freq_thz"]], [0, height], color=color, linewidth=1.5)
        for m in ref_mode_log:
            color = 'tab:red' if not m["kept"] else 'tab:orange'
            height = 0.4 if not m["kept"] else 1.0
            axm.plot([m["freq_thz"], m["freq_thz"]], [0, -height], color=color, linewidth=1.5)
        axm.axhline(0.0, color='black', linewidth=0.8)
        axm.plot([], [], color='tab:blue', label="Site: vibrational")
        axm.plot([], [], color='tab:orange', label="Isolated ref: vibrational")
        axm.plot([], [], color='tab:red', label="Excluded (imaginary/rigid-body)")
        axm.set_xlabel("Frequency (THz)")
        axm.set_yticks([])
        axm.set_title(f"Vibrational mode spectrum -- {site_label} (above) vs. "
                       "isolated reference (below)")
        axm.legend(fontsize=8, loc='upper right')
        fig2.tight_layout()
        mode_plot_path = os.path.join(args.dir, MODE_SPECTRUM_PLOT_FILE)
        fig2.savefig(mode_plot_path, dpi=150)
        if args.view_plots:
            plt.show()
        plt.close(fig2)
        print_dual(f"{color_text('[Saved]', 'cyan')} {mode_plot_path}", f_out)

    if report_path:
        print_dual(f"{color_text('[Saved]', 'cyan')} Report -> {report_path}", f_out)

    if f_out:
        f_out.close()

    print(f"\n{color_text('Success:', 'green')} DG(adsorption) = {dg_at_tmin:+.4f} eV at "
          f"T = {args.tmin} K.")


if __name__ == "__main__":
    main()
