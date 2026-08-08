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
from pymatgen.core.periodic_table import Element
from stb.core.siesta_log import get_free_energy, report_quality_diagnostics
from stb.core.cli import color_text, show_intro, print_dual, print_section
from stb.core.phonon_workflow import load_phonon_with_force_constants
from stb.adsorb_analysis import read_site_table, read_bsse_energy, GIBBS_LOCAL_META_FILE

REPORT_FILE = "adsorption_gibbs_report.txt"
PLOT_FILE = "adsorption_gibbs.png"

_FREQ_CONVERSION_THZ = 15.633302  # sqrt(eV / (amu * Ang^2)) -> THz, same constant ASE's own
                                  # vibrational-analysis code (and her_analysis.py/oer_analysis.py) use
_BOLTZMANN_EV_K = 8.617333262e-5  # eV/K
_EV_PER_THZ = 0.00413566733  # E = h*f, h in eV.s, f in THz*1e12 -> eV
_KJ_MOL_TO_EV = 0.01036427
_J_MOL_TO_EV = _KJ_MOL_TO_EV / 1000.0


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


def compute_local_hessian_thermo(gibbs_dir, temperature_k, f_out, label):
    """Reads gibbs_local_meta.json + every disp_NNN/<system_label>.FA
    force file stb-adsorbAnalysis --compute-gibbs wrote, builds the
    mass-weighted 3N x 3N dynamical matrix by central finite difference
    (same math as oer_analysis.py::compute_local_zpe_entropy, generalized
    N-atom version -- reduces to her_analysis.py's single-atom 3x3 case
    when N=1), diagonalizes it, and returns (zpe_ev, ts_ev) -- standard
    harmonic-oscillator ZPE/entropy summed over every REAL (non-imaginary,
    non-near-zero) mode. Modes at/below 0 are reported and excluded --
    for a fully isolated single atom (no local coupling at all), ALL 3
    modes are pure translation and get excluded this way, correctly
    giving zpe_ev=0.0 (a free atom has no vibrational zero-point energy).
    Returns (None, None) on any read failure.
    """
    meta_path = os.path.join(gibbs_dir, GIBBS_LOCAL_META_FILE)
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

    forces = {}  # (moved_local_k, axis, sign) -> (n_local, 3) forces on every local atom
    for i, entry in enumerate(order, start=1):
        fa_path = os.path.join(gibbs_dir, f"disp_{i:03d}", f"{system_label}.FA")
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
            Phi[:, 3 * moved_k + axis] = d_force.reshape(-1)
    Phi = 0.5 * (Phi + Phi.T)

    mass_vec = np.repeat(masses, 3)
    D = Phi / np.sqrt(np.outer(mass_vec, mass_vec))
    eigenvalues = np.linalg.eigvalsh(D)

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
        else:
            freq_thz = _FREQ_CONVERSION_THZ * np.sqrt(-eig)
            print_dual(color_text(
                f"      {freq_thz:>6.2f} THz (imaginary/translational) -- excluded from "
                "ZPE/entropy", 'yellow'), f_out)
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
    print_dual(color_text(
        "  [LIMITATION] The isolated-adsorbate reference term below is VIBRATIONAL ONLY -- no "
        "translational/rotational ideal-gas entropy. Exact for a single-atom adsorbate; an "
        "approximation for a polyatomic one. See stb-adsorbGibbs --help.", 'yellow'), f_out)
    zpe_ref, ts_ref = compute_local_hessian_thermo(isolated_gibbs_dir, args.tmin, f_out, "Isolated reference")
    if zpe_ref is None:
        if f_out:
            f_out.close()
        sys.exit(1)

    def site_thermo_at(temperature_k):
        if zpe_mode == "local":
            return compute_local_hessian_thermo(site_gibbs_dir, temperature_k, None, "Site")
        zpe_site, ts_site = compute_full_phonon_thermo(
            site_gibbs_dir, "gibbs_site", temperature_k, None, "Site")
        zpe_clean, ts_clean = compute_full_phonon_thermo(
            clean_gibbs_dir, "gibbs_clean", temperature_k, None, "Clean slab")
        if zpe_site is None or zpe_clean is None:
            return None, None
        return zpe_site - zpe_clean, ts_site - ts_clean

    zpe_site0, ts_site0 = site_thermo_at(args.tmin)
    if zpe_site0 is None:
        print_dual(color_text("[ERROR] Site vibrational/phonon calculation incomplete.", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)
    print_dual(f"  ZPE(site)  = {zpe_site0:.4f} eV   ZPE(isolated ref) = {zpe_ref:.4f} eV   "
                f"DZPE = {zpe_site0 - zpe_ref:+.4f} eV", f_out)

    print_section('[3] GIBBS FREE ENERGY vs. TEMPERATURE', f_out)
    temperatures = np.arange(args.tmin, args.tmax + args.tstep / 2, args.tstep)
    dg_values = []
    for T in temperatures:
        zpe_site, ts_site = site_thermo_at(T)
        _zpe_ref_T, ts_ref_T = compute_local_hessian_thermo(isolated_gibbs_dir, T, None, "Isolated reference")
        delta_zpe = zpe_site - zpe_ref
        delta_ts = ts_site - ts_ref_T
        dg = e_ads_used + delta_zpe - delta_ts
        dg_values.append(dg)
        print_dual(f"  T = {T:>7.2f} K   DZPE = {delta_zpe:+.4f} eV   DTS = {delta_ts:+.4f} eV   "
                    f"DG = {dg:+.4f} eV", f_out)

    dg_at_tmin = dg_values[0]
    print_dual(f"\n{color_text('[FINAL RESULT]', 'magenta')} DG(adsorption) = {dg_at_tmin:+.4f} eV "
                f"at T = {args.tmin} K (E_ads {'BSSE-corrected' if used_bsse else 'raw'})", f_out)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(temperatures, dg_values, marker='o', color='tab:blue')
    ax.axhline(e_ads_used, color='gray', linestyle='--', linewidth=1,
               label=f"E_ads ({'BSSE' if used_bsse else 'raw'}, electronic only)")
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel("DG (eV)")
    ax.set_title(f"Gibbs free energy of adsorption -- {site_label}")
    ax.legend(fontsize=9)
    fig.tight_layout()
    plot_path = os.path.join(args.dir, PLOT_FILE)
    fig.savefig(plot_path, dpi=150)
    if args.view_plots:
        plt.show()
    plt.close(fig)
    print_dual(f"{color_text('[Saved]', 'cyan')} {plot_path}", f_out)

    if report_path:
        print_dual(f"{color_text('[Saved]', 'cyan')} Report -> {report_path}", f_out)

    if f_out:
        f_out.close()

    print(f"\n{color_text('Success:', 'green')} DG(adsorption) = {dg_at_tmin:+.4f} eV at "
          f"T = {args.tmin} K.")


if __name__ == "__main__":
    main()
