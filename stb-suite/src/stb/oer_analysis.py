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
from datetime import datetime
import numpy as np
from pymatgen.core.periodic_table import Element
from stb.core.cli import color_text, show_intro, print_dual
from stb.core.siesta_log import get_free_energy, check_scf_and_force, report_quality_diagnostics
from stb.core.phonon_workflow import load_phonon_with_force_constants

REPORT_FILE = "oer_stage4.txt"

# Literature gas-phase H2 reference constant (Norskov et al., widely used throughout the
# CHE/HER/OER literature) -- reused verbatim from HER's her_analysis.py. Fixed at its 298.15 K
# literature value regardless of --temp, same simplification HER already makes: H2's own
# vibrational structure is experimentally very well known, unlike the adsorbed/isolated
# intermediates this module computes from the user's own DFT+phonon data.
H2_ZPE_EV = 0.270
H2_TS_298K_EV = 0.400

# 4 * 1.23 V, the well-established experimental reversible potential for 2H2O -> O2 + 4(H+ + e-)
# (Rossmeisl et al. 2007; Man et al. 2011) -- NOT flagged [UNVERIFIED], unlike a per-species ZPE/TS
# table would be (see the OER plan's Decision 4: no --zpe-mode standard in v1 for exactly this
# reason). G(O2) is DERIVED from this value rather than computed via DFT -- see main()'s
# description for why (Sargeant et al. 2021: >=0.3 eV systematic GGA/O2-triplet error).
FOUR_TIMES_U0_EV = 4.920

_FREQ_CONVERSION_THZ = 15.633302  # sqrt(eV / (amu * Ang^2)) -> THz, same constant HER/ASE use
_BOLTZMANN_EV_K = 8.617333262e-5  # eV/K
_EV_PER_THZ = 0.00413566733  # E = h*f, h in eV.s, f in THz*1e12 -> eV
_KJ_MOL_TO_EV = 0.01036427
_J_MOL_TO_EV = _KJ_MOL_TO_EV / 1000.0


def find_winning_dir(candidates_root, out_file):
    """Re-derives which 'site_*' folder under candidates_root has the
    lowest FreeEng -- duplicated logic (not imported) from oer_refs.py's
    own find_winning_site, same "each stage re-derives, doesn't trust a
    sibling stage's report" policy already established for HER (e.g.
    her_analysis.py's own find_winning_site_energy). Returns
    (basename, dir, energy) or (None, None, None).
    """
    if not os.path.isdir(candidates_root):
        return None, None, None
    site_dirs = sorted(
        d.path for d in os.scandir(candidates_root) if d.is_dir() and d.name.startswith("site_"))
    best_dir, best_energy = None, float('inf')
    for d in site_dirs:
        energy = get_free_energy(os.path.join(d, out_file))
        if energy is not None and energy < best_energy:
            best_energy = energy
            best_dir = d
    if best_dir is None:
        return None, None, None
    return os.path.basename(best_dir), best_dir, best_energy


def locate_intermediate_dir(output_root, name, strategy, out_file):
    """Returns the final relaxed folder for O*/OOH* ('name' is 'o' or
    'ooh'), matching oer_refs.py's own _locate_intermediate resolution
    (duplicated, not imported): 'derived' -> the single
    'intermediates/<name>_star/' folder; 'search' -> the winner among
    'intermediates/<name>_candidates/site_*/'.
    """
    if strategy == "derived":
        return os.path.join(output_root, "intermediates", f"{name}_star")
    candidates_root = os.path.join(output_root, "intermediates", f"{name}_candidates")
    _label, d, _energy = find_winning_dir(candidates_root, out_file)
    return d


def read_fa_force(fa_path, atom_index):
    """Reads the force (Fx, Fy, Fz, eV/Ang) on ONE atom (0-based
    atom_index) from a SIESTA .FA file -- format: first line = atom
    count, then one row per atom '<1-based index> Fx Fy Fz'. Same
    fail-soft "None on any parse issue" contract as core.siesta_log's own
    parsers. Duplicated from her_analysis.py.
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


def compute_local_zpe_entropy(zpe_dir, temperature_k, f_out):
    """Generalizes her_analysis.py::compute_local_zpe_entropy from ONE
    atom (3x3 Hessian, uniform mass) to N local atoms (3N x 3N Hessian):
    reads zpe_local_meta.json (local_indices/local_symbols/
    displacement_ang/order/system_label) and, for every displacement
    folder, reads the force on EVERY local atom (not just the one
    displaced) to assemble the full coupled block
    Phi[3i+a, 3j+b] = -(F_j,b(+) - F_j,b(-)) / (2*d) -- captures
    inter-atom coupling (e.g. O-H stretching within OH*/OOH*/H2O), not
    independent per-atom Hessians. Symmetrizes Phi, then builds the
    MASS-WEIGHTED dynamical matrix D[3i+a,3j+b] = Phi[3i+a,3j+b] /
    sqrt(m_i * m_j) (masses via pymatgen.core.periodic_table.Element)
    BEFORE diagonalizing -- the correct generalization of HER's
    per-eigenvalue '/mass' division (only valid there because all 3 DOF
    belonged to one atom of one mass); reduces EXACTLY to HER's formula
    when there's only 1 local atom (sqrt(m*m) == m). Eigenvalues of D are
    already omega^2, no further per-mode mass division needed.

    Returns (zpe_ev, ts_ev) for this ONE species/intermediate -- local
    mode's substrate-frozen (or, for H2O, isolated-molecule) treatment
    needs no further "delta vs. clean" subtraction, exactly as in HER.
    Returns (None, None) on any read failure.
    """
    meta_path = os.path.join(zpe_dir, "zpe_local_meta.json")
    try:
        with open(meta_path) as f:
            meta = json.load(f)
    except (OSError, ValueError):
        print_dual(color_text(f"    [ERROR] Could not read '{meta_path}'.", 'red'), f_out)
        return None, None

    local_indices = meta["local_indices"]
    local_symbols = meta["local_symbols"]
    displacement_ang = meta["displacement_ang"]
    order = meta["order"]
    system_label = meta["system_label"]
    n_local = len(local_indices)
    index_to_local = {atom_index: k for k, atom_index in enumerate(local_indices)}
    masses = np.array([Element(sym).atomic_mass for sym in local_symbols], dtype=float)

    forces = {}  # (moved_local_k, axis, sign) -> (n_local, 3) array of forces on every local atom
    for i, entry in enumerate(order, start=1):
        fa_path = os.path.join(zpe_dir, f"disp_{i:03d}", f"{system_label}.FA")
        moved_k = index_to_local[entry["atom_index"]]
        row = np.zeros((n_local, 3))
        for j, atom_index in enumerate(local_indices):
            force = read_fa_force(fa_path, atom_index)
            if force is None:
                print_dual(color_text(f"    [ERROR] Could not read forces from '{fa_path}'.", 'red'), f_out)
                return None, None
            row[j] = force
        forces[(moved_k, entry["axis"], entry["sign"])] = row

    dof = 3 * n_local
    Phi = np.zeros((dof, dof))
    for moved_k in range(n_local):
        for axis in range(3):
            f_plus = forces[(moved_k, axis, 1.0)]
            f_minus = forces[(moved_k, axis, -1.0)]
            d_force = -(f_plus - f_minus) / (2.0 * displacement_ang)  # (n_local, 3)
            col = 3 * moved_k + axis
            Phi[:, col] = d_force.reshape(-1)
    Phi = 0.5 * (Phi + Phi.T)

    mass_vec = np.repeat(masses, 3)
    D = Phi / np.sqrt(np.outer(mass_vec, mass_vec))

    eigenvalues = np.linalg.eigvalsh(D)
    zpe_ev, ts_ev = 0.0, 0.0
    print_dual(f"    Local (partial-Hessian, {n_local}-atom) vibrational modes:", f_out)
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
        else:
            freq_thz = _FREQ_CONVERSION_THZ * np.sqrt(-eig)
            print_dual(color_text(
                f"      {freq_thz:>6.2f} THz (IMAGINARY) -- excluded from ZPE/entropy", 'yellow'), f_out)
    return zpe_ev, ts_ev


def compute_full_zpe_entropy(zpe_dir, system_label, temperature_k, f_out):
    """Full-structure ZPE/entropy via phonopy's own thermal-properties
    pipeline, reusing core.phonon_workflow.load_phonon_with_force_constants
    -- identical to her_analysis.py's own compute_full_zpe_entropy (zero
    changes needed: it already operates on a whole structure's phonon
    calc regardless of atom count). Returns (zpe_ev, ts_ev) or
    (None, None) on failure.
    """
    try:
        phonon, _internal_to_angstrom, original_dir = load_phonon_with_force_constants(
            zpe_dir, system_label, False, f_out)
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
        description=f"""{color_text("Stage 4 of 4: combines every reference energy from stb-oerRefs "
        "into Delta-G1..4, the theoretical overpotential eta, and the potential-determining step "
        "(PDS).", 'bold')}
The 4-electron computational hydrogen electrode (CHE) descriptor for OER (Rossmeisl et al. 2007;
Man et al. 2011), via the 3 adsorbed intermediates OH*, O*, OOH*:
    * + H2O(l)  -> OH*  + (H+ + e-)      dG1 = G(OH*)  + 0.5*G(H2) - G(*)    - G(H2O)
    OH*         -> O*   + (H+ + e-)      dG2 = G(O*)   + 0.5*G(H2) - G(OH*)
    O* + H2O(l) -> OOH* + (H+ + e-)      dG3 = G(OOH*) + 0.5*G(H2) - G(O*)   - G(H2O)
    OOH*        -> *    + O2 + (H+ + e-) dG4 = G(*) + G(O2) - G(OOH*) + 0.5*G(H2)
eta = max(dG1..dG4) - 1.23 V; the step with the largest dG is the potential-determining step (PDS).

[IMPORTANT] No O2 SIESTA reference is ever used. DFT/GGA badly describes the O2 triplet ground
state (a systematic error of at least ~0.3 eV across common functionals -- Sargeant et al. 2021,
J. Electroanal. Chem.), so G(O2) is instead DERIVED from the experimental total reaction free
energy of the overall 4-electron reaction (4 * 1.23 V = 4.92 eV):
    G(O2) = 2*G(H2O) - 2*G(H2) + 4.92 eV
By construction dG1+dG2+dG3+dG4 == 4.92 eV exactly, independent of the actual electronic
energies -- this module prints that sum as a live internal sanity check.

No --zpe-mode 'standard' is offered: every intermediate's (and H2O's own) ZPE/entropy is computed
from your own DFT+phonon data (--zpe-mode local/full in stb-oerRefs), never a hardcoded per-species
literature default -- OER would need 4 additional constants beyond HER's single well-established
H2 value, and their exact digits could not be pinned to a verified primary source during this
tool's development.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage example:\n"
               "  %(prog)s --directory oer_study --temp 298.15\n"
    )

    parser.add_argument("-dir", "--directory", type=str, default="oer_study",
                         help="Root directory written by stb-oer/stb-oerIntermediates/stb-oerRefs "
                              "(default: oer_study).")
    parser.add_argument("--file", type=str, default="calc.out",
                         help="SIESTA output filename inside each folder (default: calc.out).")
    parser.add_argument("--temp", type=float, default=298.15,
                         help="Temperature in Kelvin, for the local/full ZPE modes' entropy term "
                              "(default: 298.15). Has no effect on H2's own thermal term (fixed at "
                              "its 298.15 K literature value regardless of --temp, see the module "
                              "docstring's H2_ZPE_EV/H2_TS_298K_EV).")
    parser.add_argument("--force-tolerance", type=float, default=0.05,
                         help="Residual atomic force in eV/Ang (default: 0.05) above which a "
                              "folder is flagged as possibly not relaxed/converged. Advisory only.")
    parser.add_argument("-o", "--output", type=str, default="OER_report",
                         help="Base filename (no extension) for the final report (default: OER_report).")
    parser.add_argument("-v", "--version", action="version", version=f"stb-oerAnalysis {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("OER WORKFLOW -- STAGE 4: ANALYSIS", 'bold'))
    print("-" * 60)

    output_root = args.directory
    stage2_report = os.path.join(output_root, "oer_stage2.txt")
    stage3_report = os.path.join(output_root, "oer_stage3.txt")
    if not os.path.isfile(stage3_report):
        print(color_text(f"[ERROR] '{stage3_report}' not found -- run stb-oerRefs (Stage 3) first.", 'red'))
        sys.exit(1)
    if not os.path.isfile(stage2_report):
        print(color_text(f"[ERROR] '{stage2_report}' not found -- run stb-oerIntermediates "
                          "(Stage 2) first.", 'red'))
        sys.exit(1)

    o_strategy, ooh_strategy = None, None
    with open(stage2_report) as f:
        for line in f:
            if line.startswith("O* strategy"):
                o_strategy = line.split(":", 1)[1].strip()
            elif line.startswith("OOH* strategy"):
                ooh_strategy = line.split(":", 1)[1].strip()

    bsse_mode, zpe_mode = None, None
    with open(stage3_report) as f:
        for line in f:
            if line.startswith("BSSE mode"):
                bsse_mode = line.split(":", 1)[1].strip()
            elif line.startswith("ZPE mode"):
                zpe_mode = line.split(":", 1)[1].strip()
    if o_strategy is None or ooh_strategy is None or bsse_mode is None or zpe_mode is None:
        print(color_text("[ERROR] Could not recover strategy/BSSE mode/ZPE mode from the Stage "
                          "2/3 reports.", 'red'))
        sys.exit(1)

    sites_root = os.path.join(output_root, "sites")
    winning_oh_label, winning_oh_dir, _e = find_winning_dir(sites_root, args.file)
    o_dir = locate_intermediate_dir(output_root, "o", o_strategy, args.file)
    ooh_dir = locate_intermediate_dir(output_root, "ooh", ooh_strategy, args.file)
    if winning_oh_dir is None or o_dir is None or ooh_dir is None:
        print(color_text("[ERROR] Could not locate the winning OH*/O*/OOH* folder(s) -- did every "
                          "relaxation finish?", 'red'))
        sys.exit(1)

    def read_energy(folder, label, f_out):
        out_path = os.path.join(output_root, folder, args.file)
        energy = get_free_energy(out_path)
        if energy is None:
            print_dual(color_text(f"[ERROR] Could not read energy from '{out_path}'.", 'red'), f_out)
            sys.exit(1)
        print_dual(f"  {label:<14} = {energy:>12.6f} eV  ({out_path})", f_out)
        report_quality_diagnostics(label, out_path, args.force_tolerance, f_out)
        return energy

    def read_energy_at(path, label, f_out):
        out_path = os.path.join(path, args.file)
        energy = get_free_energy(out_path)
        if energy is None:
            print_dual(color_text(f"[ERROR] Could not read energy from '{out_path}'.", 'red'), f_out)
            sys.exit(1)
        print_dual(f"  {label:<14} = {energy:>12.6f} eV  ({out_path})", f_out)
        report_quality_diagnostics(label, out_path, args.force_tolerance, f_out)
        return energy

    report_path = os.path.join(output_root, REPORT_FILE)
    with open(report_path, "w") as f_out:
        print_dual(f"{color_text('===== OER STAGE 4 REPORT (ANALYSIS) =====', 'magenta')}", f_out)

        print_dual(f"\n{color_text('[0] RUN METADATA', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        print_dual(f"Date/time       : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", f_out)
        print_dual(f"Directory       : {output_root}", f_out)
        print_dual(f"O* strategy     : {o_strategy}   OOH* strategy: {ooh_strategy}", f_out)
        print_dual(f"BSSE mode       : {bsse_mode}", f_out)
        print_dual(f"ZPE mode        : {zpe_mode}", f_out)
        print_dual(f"Temperature     : {args.temp} K", f_out)

        print_dual(f"\n{color_text('[1] ELECTRONIC ENERGIES (FreeEng)', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        e_clean = read_energy("00_clean_slab", "E_clean", f_out)
        e_oh = read_energy_at(winning_oh_dir, "E_OH", f_out)
        e_o = read_energy_at(o_dir, "E_O", f_out)
        e_ooh = read_energy_at(ooh_dir, "E_OOH", f_out)
        e_h2 = read_energy("02_h2_molecule", "E_H2", f_out)
        e_h2o = read_energy("03_h2o_molecule", "E_H2O", f_out)
        e_deformed = read_energy("04_slab_deformed", "E_deformed", f_out)

        print_dual(f"\n{color_text('[2] BSSE CORRECTION', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        e_deformation = e_deformed - e_clean
        print_dual(f"  Deformation cost (diagnostic, not used below) : {e_deformation:>10.4f} eV", f_out)

        def bsse_shift(prefix, f_out):
            # slab_only is read at THIS triad's own geometry (not the diagnostic
            # '04_slab_deformed', which is always at OH*'s geometry) -- the counterpoise
            # difference must compare the SAME geometry across basis sets, and in
            # --bsse-mode shared the triad geometry is OOH*'s, not OH*'s.
            e_slab_only = read_energy(f"{prefix}_slab_only", f"{prefix}_slab_only", f_out)
            e_slab_ghost = read_energy(f"{prefix}_slab_ghost", f"{prefix}_slab_gh", f_out)
            e_ads_ghost = read_energy(f"{prefix}_adsorbate_ghost_slab", f"{prefix}_ads_gh", f_out)
            e_isolated = read_energy(f"{prefix}_isolated", f"{prefix}_iso", f_out)
            bsse_slab = e_slab_only - e_slab_ghost
            bsse_ads = e_isolated - e_ads_ghost
            total = bsse_slab + bsse_ads
            print_dual(f"  BSSE ({prefix}): slab={bsse_slab:+.4f} eV  adsorbate={bsse_ads:+.4f} eV  "
                        f"total={total:+.4f} eV", f_out)
            return total

        if bsse_mode == "shared":
            shared_correction = bsse_shift("05_bsse", f_out)
            bsse_oh = bsse_o = bsse_ooh = shared_correction
        else:
            bsse_oh = bsse_shift("05_bsse_OH", f_out)
            bsse_o = bsse_shift("05_bsse_O", f_out)
            bsse_ooh = bsse_shift("05_bsse_OOH", f_out)

        e_oh_corr = e_oh + bsse_oh
        e_o_corr = e_o + bsse_o
        e_ooh_corr = e_ooh + bsse_ooh
        print_dual(f"  E_OH (corrected)  = {e_oh_corr:.4f} eV   E_O (corrected)  = {e_o_corr:.4f} "
                    f"eV   E_OOH (corrected) = {e_ooh_corr:.4f} eV", f_out)

        print_dual(f"\n{color_text('[3] THERMAL CORRECTION', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)

        def thermal_term(name, dir_suffix, phonon_label, f_out):
            if zpe_mode == "local":
                zpe_dir = os.path.join(output_root, f"08_zpe_calc_{dir_suffix}")
                zpe, ts = compute_local_zpe_entropy(zpe_dir, args.temp, f_out)
            else:
                zpe_dir = os.path.join(output_root, f"09_zpe_calc_{dir_suffix}")
                zpe, ts = compute_full_zpe_entropy(zpe_dir, phonon_label, args.temp, f_out)
            if zpe is None:
                print_dual(color_text(f"[ERROR] Could not compute {name}'s ZPE/entropy.", 'red'), f_out)
                sys.exit(1)
            print_dual(f"  ZPE({name}) = {zpe:.4f} eV   TS({name}) = {ts:.4f} eV", f_out)
            return zpe, ts

        zpe_h2o, ts_h2o = thermal_term("H2O", "H2O", "oer_zpe_h2o", f_out)
        if zpe_mode == "full":
            zpe_clean, ts_clean = thermal_term("clean slab", "clean", "oer_zpe_clean", f_out)
        else:
            zpe_clean, ts_clean = 0.0, 0.0  # local mode: substrate frozen, delta vs. clean is 0 by construction

        zpe_oh_raw, ts_oh_raw = thermal_term("OH*", "OH", "oer_zpe_oh", f_out)
        zpe_o_raw, ts_o_raw = thermal_term("O*", "O", "oer_zpe_o", f_out)
        zpe_ooh_raw, ts_ooh_raw = thermal_term("OOH*", "OOH", "oer_zpe_ooh", f_out)

        d_zpe_oh = (zpe_oh_raw - zpe_clean)
        d_ts_oh = (ts_oh_raw - ts_clean)
        d_zpe_o = (zpe_o_raw - zpe_clean)
        d_ts_o = (ts_o_raw - ts_clean)
        d_zpe_ooh = (zpe_ooh_raw - zpe_clean)
        d_ts_ooh = (ts_ooh_raw - ts_clean)

        g_oh = e_oh_corr + d_zpe_oh - args.temp * d_ts_oh
        g_o = e_o_corr + d_zpe_o - args.temp * d_ts_o
        g_ooh = e_ooh_corr + d_zpe_ooh - args.temp * d_ts_ooh
        g_h2 = e_h2 + H2_ZPE_EV - args.temp * H2_TS_298K_EV
        g_h2o = e_h2o + zpe_h2o - args.temp * ts_h2o
        g_clean = e_clean

        print_dual(f"\n{color_text('[4] GAS-PHASE O2 REFERENCE (DERIVED)', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        g_o2 = 2.0 * g_h2o - 2.0 * g_h2 + FOUR_TIMES_U0_EV
        print_dual(f"  G(H2)  = {g_h2:.4f} eV   G(H2O,l) = {g_h2o:.4f} eV", f_out)
        print_dual(f"  G(O2) [derived, NOT computed via DFT] = 2*G(H2O) - 2*G(H2) + 4.92 = "
                    f"{g_o2:.4f} eV", f_out)

        print_dual(f"\n{color_text('[5] REACTION STEPS & POTENTIAL-DETERMINING STEP', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        dG1 = g_oh + 0.5 * g_h2 - g_clean - g_h2o
        dG2 = g_o + 0.5 * g_h2 - g_oh
        dG3 = g_ooh + 0.5 * g_h2 - g_o - g_h2o
        dG4 = g_clean + g_o2 - g_ooh + 0.5 * g_h2
        steps = [("1: * + H2O -> OH*", dG1), ("2: OH* -> O*", dG2), ("3: O* + H2O -> OOH*", dG3),
                 ("4: OOH* -> * + O2", dG4)]
        pds_index = int(np.argmax([dG1, dG2, dG3, dG4]))
        for i, (label, dg) in enumerate(steps):
            marker = color_text(" <-- PDS", 'yellow') if i == pds_index else ""
            print_dual(f"  Step {label:<22} dG = {dg:+.4f} eV{marker}", f_out)
        step_sum = dG1 + dG2 + dG3 + dG4
        sanity_ok = abs(step_sum - FOUR_TIMES_U0_EV) < 1e-6
        print_dual(f"  Sum dG1+dG2+dG3+dG4 = {step_sum:.6f} eV (should equal {FOUR_TIMES_U0_EV} eV "
                    f"exactly by construction) -- {'OK' if sanity_ok else color_text('MISMATCH -- BUG', 'red')}",
                    f_out)

        eta = max(dG1, dG2, dG3, dG4) - 1.23
        print_dual(f"\n{color_text('[6] FINAL RESULT', 'magenta')}", f_out)
        print_dual("=" * 60, f_out)
        print_dual(f"  Theoretical overpotential eta = {eta:+.4f} V", f_out)
        print_dual(f"  Potential-determining step (PDS) = Step {pds_index + 1}", f_out)
        print_dual("=" * 60, f_out)
        verdict = ("excellent OER catalyst (eta close to 0)" if eta < 0.3
                   else "moderate OER activity" if eta < 0.6
                   else "poor OER catalyst (large overpotential)")
        print_dual(f"  Qualitative assessment: {verdict}", f_out)

        report_str_path = os.path.join(output_root, f"{args.output}.txt")
        with open(report_str_path, "w") as f_report:
            f_report.write(f"OER eta = {eta:+.4f} V, PDS = Step {pds_index + 1}\n")
            f_report.write(f"BSSE mode = {bsse_mode}, ZPE mode = {zpe_mode}, T = {args.temp} K\n")
            f_report.write(f"dG1={dG1:+.4f}  dG2={dG2:+.4f}  dG3={dG3:+.4f}  dG4={dG4:+.4f} eV\n")

        print_dual(f"\n{color_text('[7] SUMMARY & FILES', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        print_dual(f"Winning OH* site    : {winning_oh_label}", f_out)
        print_dual(f"eta                 : {eta:+.4f} V", f_out)
        print_dual(f"Report              : {report_path}", f_out)
        print_dual(f"Files               : {report_str_path}", f_out)

    print("\n[INFO] Complete job!")
    print("\n" + "-" * 60)
    print(color_text("OER analysis complete.\n", 'bold'))


if __name__ == "__main__":
    main()
