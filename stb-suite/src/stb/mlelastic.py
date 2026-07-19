#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""Standalone elastic-constants calculation driven entirely by a MACE
potential -- the ML Simulations analog of the SIESTA-based stb-elasticInputs
/stb-elasticAnalysis prep+analysis pair, same relationship stb-mlphonons has
to stb-phononsCreate/stb-phononsPos. Unlike that pair, this is a single
self-contained tool (no folders to hand off to SIESTA and reload later):
MACE stress is cheap enough to compute in-memory for every canonical Voigt
strain direction in one run, so there's no DFT-cost reason to do the
symmetry-based direction-reduction stb-elasticInputs' 'all' default performs,
and no separate --method energy path either (that method exists in the DFT
tool specifically to sidestep SIESTA's real-space-grid "eggbox effect" on
stress, which MACE has no analog of -- the direct stress-strain method is
already reliable here).

Deliberately reuses (not duplicates) elastic_analysis.py's post-fit numerics
-- compute_stiffness_matrix, direction_fit_diagnostics, tensor_symmetry_check,
check_stability_and_report, emit_elastic_report, and its physical-unit
constants/tables -- since a stress-strain elastic-constant fit is the exact
same linear algebra regardless of whether the (strain, stress) pairs came
from SIESTA or MACE; only the DATA MINING step (reading calc.out per
strain_*/ folder there) differs, replaced here by compute_ml_elastic_data's
in-memory MACE evaluation loop. Also reuses elastic_inputs.py's
get_strain_matrix/direction_axes for the exact same deformation-matrix/
vacuum-blocking convention stb-elasticInputs itself uses, so a strain
labeled 'xy' means bit-for-bit the same deformation in both tools.

Plots are new matplotlib code (plot_stress_strain), not elastic_analysis.py's
write_stress_plots (gnuplot .dat+.gplot pairs, the older convention) --
matching the newer matplotlib convention this session's ML tools already use
(mlphonons.py, mlmd.py, aimd_analysis.py, mlff_analysis.py).

Four features on top of the original single-model, single-strain-range
calculation:
  - Foundation-model comparison (--custom-model + overlay/table against the
    raw foundation model), same pattern as stb-mlphonons/stb-mlffAnalysis.
  - Sound velocities + Debye temperature (3D only), derived from the
    Voigt-Reuss-Hill bulk/shear moduli check_stability_and_report already
    computes plus the structure's own mass density -- zero extra MACE cost.
  - Directional Young's modulus anisotropy surface (3D only), from the full
    compliance tensor built out of the fitted C_sym.
  - --check-convergence: re-runs the fit at multiple --max strain
    magnitudes to check the linear-regime assumption actually holds.
"""

VERSION = "1.1.0"

import os
import sys
import copy
import argparse
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 -- registers the '3d' projection
from pymatgen.io.ase import AseAtomsAdaptor
from scipy.constants import hbar, k as BOLTZMANN_K, atomic_mass

from stb.core import structure_io, kspace, mace_relax
from stb.core import symmetry as core_symmetry
from stb.core.cli import color_text, show_intro, print_dual, print_section
from stb.core.deps import require_mace
from stb.elastic_inputs import get_strain_matrix, direction_axes
from stb.elastic_analysis import (
    compute_stiffness_matrix, direction_fit_diagnostics, tensor_symmetry_check,
    emit_elastic_report,
    DIRECTION_TO_COLUMN, SHEAR_DIRECTIONS, UNIT_LABELS,
    CONV_EVA3_TO_GPA, CONV_EVA2_TO_NM, CONV_EVA1_TO_NN,
)

require_mace()

REPORT_FILE = "stb_mlelastic_report.txt"
_SERIES_COLORS = ['tab:blue', 'tab:orange']


def compute_ml_elastic_data(atoms_ref, calc, args, periodic_axes, f_out, tag=""):
    """For every canonical Voigt direction not blocked by a vacuum-padded
    axis, strains the reference cell at --steps points in [-max, max]%
    (get_strain_matrix, imported from elastic_inputs.py -- the identical
    deformation convention stb-elasticInputs uses), relaxes atomic
    positions only at fixed (strained) cell unless --no-relax (the
    internal-strain relaxation real elastic-constant methodology needs --
    e.g. diamond-structure C44 depends on it), and records the resulting
    MACE stress tensor. This is the ML-driven, in-memory equivalent of
    running one SIESTA calculation per strain_*/ folder stb-elasticInputs
    would generate. `tag` (e.g. "fine-tuned") only labels progress lines,
    for a caller running this twice (model comparison).

    Returns data ({mode: {'eps': [...], 'stress': [3x3 array, ...]}}),
    exactly the format compute_stiffness_matrix/direction_fit_diagnostics/
    tensor_symmetry_check (imported from elastic_analysis.py) already
    expect -- both DFT- and ML-sourced elastic-constant fits share the same
    downstream numerics, only the data-mining step differs.
    """
    prefix = f"[{tag}] " if tag else ""
    axis_letters = ['x', 'y', 'z']
    dirs_to_run = []
    for d in core_symmetry.VOIGT_MODES:
        touched = direction_axes(d)
        blocked = sorted(touched - periodic_axes)
        if blocked:
            letters = ', '.join(axis_letters[i] for i in blocked)
            print_dual(color_text(
                f"{prefix}[WARNING] Skipping direction '{d}': touches vacuum-padded axis/axes "
                f"{letters} -- straining a vacuum gap doesn't correspond to any physical "
                "deformation.", 'yellow'), f_out)
            continue
        dirs_to_run.append(d)

    strains_pct = np.linspace(-args.max, args.max, args.steps)
    ref_lattice = np.array(atoms_ref.get_cell())
    scaled_positions = atoms_ref.get_scaled_positions()

    data = {}
    for d in dirs_to_run:
        eps_list, stress_list = [], []
        print_dual(f"{prefix}[INFO] Direction '{d}': {len(strains_pct)} strain point(s)", f_out)
        for s in strains_pct:
            s_clean = 0.0 if abs(s) < 1e-9 else s
            delta = s_clean / 100.0
            def_matrix = get_strain_matrix(d, delta)
            new_lattice = ref_lattice @ def_matrix

            deformed = atoms_ref.copy()
            deformed.set_cell(new_lattice, scale_atoms=False)
            deformed.set_scaled_positions(scaled_positions)
            deformed.calc = calc

            if args.relax and len(deformed) > 1:
                mace_relax.relax(deformed, calc, cell_mask=None, fmax=args.fmax,
                                  max_steps=args.max_steps)

            stress = deformed.get_stress(voigt=False)
            eps_list.append(delta)
            stress_list.append(stress)

        data[d] = {'eps': eps_list, 'stress': stress_list}

    return data, dirs_to_run


def run_elastic_pipeline(atoms_ref, calc, args, vacuum_axes, periodic_axes, dimensionality,
                          f_out, tag=""):
    """Full single-model elastic-constants pipeline: optional pre-relax
    (cell+ions), per-direction strain-relax-stress data generation
    (compute_ml_elastic_data), stiffness-matrix fit (compute_stiffness_matrix,
    imported from elastic_analysis.py), and numerical-quality diagnostics
    (direction_fit_diagnostics/tensor_symmetry_check, same module). `tag` is
    a short label (e.g. "fine-tuned", "foundation (small)") used only in
    printed progress messages, letting a caller running this twice (model
    comparison, or a strain-convergence scan) tell which run produced which
    numbers.

    Returns a dict: {"C_sym", "data", "dirs_run", "fit_diagnostics",
    "symmetry_check", "conv_factor", "missing_directions", "atoms_relaxed"}.
    """
    prefix = f"[{tag}] " if tag else ""
    atoms = atoms_ref.copy()
    atoms.calc = calc

    if args.pre_relax:
        cell_mask = mace_relax.build_cell_mask(vacuum_axes)
        f0 = np.abs(atoms.get_forces()).max()
        converged, steps_used = mace_relax.relax(
            atoms, calc, cell_mask=cell_mask, fmax=args.fmax, max_steps=args.max_steps)
        f1 = np.abs(atoms.get_forces()).max()
        print_dual(f"{prefix}Pre-relax (cell+ions): max|F| {f0:.4f} -> {f1:.4f} eV/Ang "
                   f"({'converged' if converged else 'hit step cap, not fully converged'}, "
                   f"{steps_used} steps)", f_out)
    else:
        f0 = np.abs(atoms.get_forces()).max()
        print_dual(color_text(
            f"{prefix}[WARNING] --no-pre-relax: straining the structure as given (residual "
            f"max|F| = {f0:.4f} eV/Ang).", 'yellow'), f_out)

    # Geometric dilution correction: MACE's stress (like SIESTA's) is a
    # volumetric quantity computed over the FULL cell, vacuum included --
    # 2D needs the cell height (Lz) to recover per-area stiffness, 1D the
    # transverse cross-section (Lx*Ly), same physical correction
    # stb-elasticAnalysis applies from a parsed .out file, computed here
    # directly from the (already pre-relaxed) reference cell instead.
    unit_label = UNIT_LABELS[dimensionality]
    cell_now = np.array(atoms.get_cell())
    if dimensionality == "2d":
        conv_factor = cell_now[2, 2] * CONV_EVA2_TO_NM
        print_dual(f"{prefix}Normalization: cell height (Lz) = {cell_now[2,2]:.4f} Ang", f_out)
    elif dimensionality == "1d":
        cross_section = np.linalg.norm(cell_now[0]) * np.linalg.norm(cell_now[1])
        conv_factor = cross_section * CONV_EVA1_TO_NN
        print_dual(f"{prefix}Normalization: vacuum cross-section (Lx*Ly) = "
                   f"{cross_section:.4f} Ang^2", f_out)
    else:
        conv_factor = CONV_EVA3_TO_GPA

    data, dirs_run = compute_ml_elastic_data(atoms, calc, args, periodic_axes, f_out, tag=tag)

    C, filled_by_symmetry, missing_directions = compute_stiffness_matrix(
        data, conv_factor, None, unit_label, f_out)
    C_sym = 0.5 * (C + C.T)
    fit_diagnostics = direction_fit_diagnostics(data, conv_factor, args.fit_quality_tolerance)
    symmetry_check = tensor_symmetry_check(C)

    return {
        "C_sym": C_sym, "data": data, "dirs_run": dirs_run,
        "fit_diagnostics": fit_diagnostics, "symmetry_check": symmetry_check,
        "conv_factor": conv_factor, "missing_directions": missing_directions,
        "atoms_relaxed": atoms,
    }


def plot_stress_strain(series, unit_label, out_path):
    """series: list of (label, data, fit_diagnostics, conv_factor, color).
    One grid panel per canonical Voigt direction measured by the FIRST
    series (assumed to be the superset when comparing 2 models -- both
    always attempt every direction the structure's dimensionality allows),
    stress-strain data + linear fit overlaid, one series' worth of points
    +fit per panel in its own color. Matplotlib equivalent of
    elastic_analysis.py's write_stress_plots (gnuplot .dat+.gplot pairs).
    Returns False (nothing written) if the first series has no direction.
    """
    label0, data0, fit0, _, _ = series[0]
    fit_by_mode0 = {r['mode']: r for r in fit0}
    modes = [m for m in core_symmetry.VOIGT_MODES if m in data0 and m in fit_by_mode0]
    if not modes:
        return False

    ncols = min(3, len(modes))
    nrows = int(np.ceil(len(modes) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 3.8 * nrows), squeeze=False)
    for idx, mode in enumerate(modes):
        ax = axes[idx // ncols][idx % ncols]
        eng_factor = 2.0 if mode in SHEAR_DIRECTIONS else 1.0
        i, j = core_symmetry.VOIGT_TO_TENSOR[DIRECTION_TO_COLUMN[mode]]

        for label, data, fit_diagnostics, conv_factor, color in series:
            fit_by_mode = {r['mode']: r for r in fit_diagnostics}
            if mode not in data or mode not in fit_by_mode:
                continue
            fit = fit_by_mode[mode]
            S = np.array(data[mode]['stress'])
            x = np.array(data[mode]['eps']) * eng_factor
            y = conv_factor * S[:, i, j]

            ax.scatter(x, y, color=color, label=f"{label} data")
            if x.max() > x.min():
                xs = np.linspace(x.min(), x.max(), 20)
                ax.plot(xs, fit['slope'] * xs + fit['residual_stress'], color=color,
                        linestyle='--', label=f"{label} fit (C={fit['slope']:.1f})")

        ax.set_title(mode.upper(), fontsize=9)
        ax.set_xlabel("Engineering strain (fraction)")
        ax.set_ylabel(f"Stress ({unit_label})")
        ax.legend(fontsize=6)

    for idx in range(len(modes), nrows * ncols):
        axes[idx // ncols][idx % ncols].axis('off')
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return True


def compute_sound_velocities(C_sym, atoms):
    """Elastic (bulk) wave velocities (longitudinal Vp, shear Vs, average
    Vm) and Debye temperature from the Voigt-Reuss-Hill bulk/shear moduli
    (same averages check_stability_and_report already computes for the 3D
    case) plus the structure's own mass density -- standard derived
    properties any elastic-constants tool reports, at zero extra MACE cost
    (pure post-processing of the already-fitted C_sym). 3D only (the
    isotropic-average formulas below assume a 3D solid); callers should
    only call this when dimensionality == "3d".

    Debye temperature uses the standard formula theta_D = (hbar/k_B) *
    (6*pi^2*n)^(1/3) * Vm, where n is the ATOM number density (not unit
    cells or formula units) -- consistent with using the average sound
    velocity Vm derived from the same isotropic elastic continuum
    approximation.

    Returns None if the moduli/density aren't physically sane (e.g. an
    unstable structure with non-positive bulk/shear modulus) -- the same
    situation check_stability_and_report already flags as UNSTABLE.
    """
    B_v = (C_sym[0, 0] + C_sym[1, 1] + C_sym[2, 2]
           + 2 * (C_sym[0, 1] + C_sym[0, 2] + C_sym[1, 2])) / 9.0
    G_v = ((C_sym[0, 0] + C_sym[1, 1] + C_sym[2, 2])
           - (C_sym[0, 1] + C_sym[0, 2] + C_sym[1, 2])
           + 3 * (C_sym[3, 3] + C_sym[4, 4] + C_sym[5, 5])) / 15.0
    try:
        S_mat = np.linalg.inv(C_sym)
        B_r = 1.0 / ((S_mat[0, 0] + S_mat[1, 1] + S_mat[2, 2])
                     + 2 * (S_mat[0, 1] + S_mat[0, 2] + S_mat[1, 2]))
        G_r = 15.0 / (4 * (S_mat[0, 0] + S_mat[1, 1] + S_mat[2, 2])
                      - 4 * (S_mat[0, 1] + S_mat[0, 2] + S_mat[1, 2])
                      + 3 * (S_mat[3, 3] + S_mat[4, 4] + S_mat[5, 5]))
    except np.linalg.LinAlgError:
        B_r, G_r = B_v, G_v
    B_h = (B_v + B_r) / 2.0
    G_h = (G_v + G_r) / 2.0
    if B_h <= 0 or G_h <= 0:
        return None

    mass_kg = np.sum(atoms.get_masses()) * atomic_mass
    volume_m3 = atoms.get_volume() * 1e-30
    if volume_m3 <= 0:
        return None
    rho = mass_kg / volume_m3

    B_pa, G_pa = B_h * 1e9, G_h * 1e9
    v_p = np.sqrt((B_pa + 4.0 / 3.0 * G_pa) / rho)
    v_s = np.sqrt(G_pa / rho)
    v_m = (1.0 / 3.0 * (2.0 / v_s ** 3 + 1.0 / v_p ** 3)) ** (-1.0 / 3.0)

    n_density = len(atoms) / volume_m3
    theta_d = (hbar / BOLTZMANN_K) * (6 * np.pi ** 2 * n_density) ** (1.0 / 3.0) * v_m

    return {"rho": rho, "v_p": v_p, "v_s": v_s, "v_m": v_m, "theta_D": theta_d}


def compliance_voigt_to_full(S):
    """Converts a 6x6 Voigt COMPLIANCE matrix (engineering-strain
    convention) into the full 3x3x3x3 compliance tensor s_ijkl, applying
    the standard engineering-to-tensor-strain correction factor: 1 for a
    normal-normal entry, 2 for one shear index, 4 for two shear indices.
    This is NOT the same factor rule a STIFFNESS Voigt-to-tensor conversion
    uses (stiffness has no such factor) -- kept as its own small helper
    rather than risking reusing a stiffness-oriented tensor conversion for
    compliance data.
    """
    s_full = np.zeros((3, 3, 3, 3))
    for m, (i, j) in enumerate(core_symmetry.VOIGT_TO_TENSOR):
        for n, (k, l) in enumerate(core_symmetry.VOIGT_TO_TENSOR):
            factor = (2.0 if m >= 3 else 1.0) * (2.0 if n >= 3 else 1.0)
            val = S[m, n] / factor
            for (a, b) in ((i, j), (j, i)):
                for (c, d) in ((k, l), (l, k)):
                    s_full[a, b, c, d] = val
    return s_full


def directional_youngs_modulus(s_full, direction):
    """1/E(n) = s_ijkl * n_i * n_j * n_k * n_l (Einstein summation, unit
    Cartesian direction n) -- the standard formula for the direction
    -dependent Young's modulus of an anisotropic elastic tensor (e.g. Nye,
    'Physical Properties of Crystals'). Returns NaN if the direction
    somehow gives a non-positive compliance sum (shouldn't happen for a
    physically stable C_sym).
    """
    n = np.asarray(direction, dtype=float)
    n = n / np.linalg.norm(n)
    inv_e = np.einsum('ijkl,i,j,k,l->', s_full, n, n, n, n)
    return 1.0 / inv_e if inv_e > 0 else np.nan


def plot_anisotropy_surface(C_sym, out_path, n_theta=60, n_phi=120):
    """Directional Young's modulus E(n) sampled over a full spherical grid
    (n_theta x n_phi points) and plotted as a 3D "material property
    surface" -- radius in every direction equals E(n) itself (GPa), same
    visual convention as tools like ELATE use to show elastic anisotropy at
    a glance. Colored by E magnitude (viridis) for readability alongside
    the pure shape. Returns (E_min, E_max), or (None, None) if C_sym isn't
    invertible (e.g. a rank-deficient tensor from missing directions --
    caller should already be guarding against this via missing_directions,
    but a structurally singular C_sym must never crash the whole run).
    """
    try:
        S = np.linalg.inv(C_sym)
    except np.linalg.LinAlgError:
        return None, None
    s_full = compliance_voigt_to_full(S)

    thetas = np.linspace(0, np.pi, n_theta)
    phis = np.linspace(0, 2 * np.pi, n_phi)
    E = np.zeros((n_theta, n_phi))
    for it, theta in enumerate(thetas):
        for ip, phi in enumerate(phis):
            n = np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])
            E[it, ip] = directional_youngs_modulus(s_full, n)

    X = E * np.outer(np.sin(thetas), np.cos(phis))
    Y = E * np.outer(np.sin(thetas), np.sin(phis))
    Z = E * np.outer(np.cos(thetas), np.ones_like(phis))

    e_min, e_max = float(np.nanmin(E)), float(np.nanmax(E))
    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection='3d')
    norm = plt.Normalize(e_min, e_max)
    ax.plot_surface(X, Y, Z, facecolors=plt.cm.viridis(norm(E)), rstride=1, cstride=1,
                     antialiased=True, shade=False)
    ax.set_xlabel("E_x (GPa)")
    ax.set_ylabel("E_y (GPa)")
    ax.set_zlabel("E_z (GPa)")
    ax.set_title(f"Directional Young's modulus (GPa)\nmin={e_min:.1f}, max={e_max:.1f}, "
                 f"ratio={e_max/e_min:.3f}" if e_min > 0 else "Directional Young's modulus (GPa)")
    mappable = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
    mappable.set_array(E)
    fig.colorbar(mappable, ax=ax, shrink=0.6, label="E (GPa)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return e_min, e_max


def plot_strain_convergence(rows, modes, unit_label, out_path):
    """rows: list of (max_strain_pct, result_dict) from a --check-convergence
    scan. One line per diagonal Voigt constant actually present (`modes`),
    stiffness constant vs. the --max strain magnitude used to fit it --
    a flat line means the linear-regime assumption holds up to that strain;
    drift signals anharmonic contamination at larger strains.
    """
    fig, ax = plt.subplots(figsize=(7, 5))
    strains = [r[0] for r in rows]
    for mode in modes:
        idx = DIRECTION_TO_COLUMN[mode]
        vals = [r[1]["C_sym"][idx, idx] for r in rows]
        ax.plot(strains, vals, marker='o', label=f"C_{mode}")
    ax.set_xlabel("Max strain magnitude (%)")
    ax.set_ylabel(f"Diagonal stiffness constant ({unit_label})")
    ax.set_title("Strain-amplitude convergence")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Standalone elastic-constants calculation driven entirely by a MACE "
                     "potential -- the MACE-MP-0 foundation model, or a custom model "
                     "fine-tuned on your own SIESTA data via stb-mlffAnalysis "
                     "(--custom-model). Strains the reference cell along every physically "
                     "valid canonical Voigt direction, computes MACE stress at each strain "
                     "point, and fits the full stiffness matrix -- no SIESTA and no separate "
                     "strain_*/ folders needed. Also reports sound velocities/Debye "
                     "temperature and a directional Young's-modulus anisotropy surface (3D "
                     "only), or, with --check-convergence, how the fit changes across "
                     "multiple strain magnitudes. A fast heuristic preview, not a substitute "
                     "for the real SIESTA-based stb-elasticInputs/stb-elasticAnalysis "
                     "workflow.",
        epilog="Example usage:\n"
               "  stb-mlelastic --file structure.fdf\n"
               "  stb-mlelastic --file structure.fdf --custom-model mlff_model.model "
               "--max 1.0 --steps 6\n"
               "  stb-mlelastic --file structure.fdf --check-convergence\n",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--file", required=True, help="Input structure (.fdf).")
    parser.add_argument("--model", choices=["small", "medium", "large"], default="small",
                        help="MACE-MP-0 model size (default: small). Ignored if --custom-model is given.")
    parser.add_argument("--custom-model", default=None, metavar="PATH",
                        help="Use a custom MACE model file instead of the MACE-MP-0 foundation "
                             "potential -- e.g. one fine-tuned via stb-mlffAnalysis. Overrides --model.")
    parser.add_argument("--skip-foundation-comparison", action="store_true",
                        help="If --custom-model is given, also evaluate the same "
                             "elastic-constants pipeline with the raw (non-fine-tuned) "
                             "foundation model (--model) and overlay both on the stress-"
                             "strain plot plus a side-by-side constants table. On by "
                             "default; pass this to skip it (--check-convergence never does "
                             "this comparison, regardless of this flag).")
    parser.add_argument("--no-pre-relax", dest="pre_relax", action="store_false",
                        help="Skip the full MACE pre-relax (cell + positions) of the reference "
                             "structure before straining. On by default -- straining around a "
                             "structure that isn't already at a MACE-relaxed minimum shows up "
                             "as a non-negligible residual stress at zero strain in the fit "
                             "diagnostics, contaminating every constant.")
    parser.add_argument("--no-relax", dest="relax", action="store_false",
                        help="Skip the per-strain-point atomic-position relax (cell fixed at "
                             "the strained shape). On by default -- real elastic-constant "
                             "methodology relaxes internal coordinates at fixed strained cell "
                             "(matters for e.g. diamond-structure C44); skip only for a single "
                             "-atom-per-primitive-cell material, or a quick sanity check.")
    parser.add_argument("--fmax", type=float, default=0.02,
                        help="Force convergence threshold for both relax stages, eV/Ang "
                             "(default: 0.02 -- tighter than stb-mlrelax's 0.05 default, since "
                             "elastic constants are a second derivative and more sensitive to "
                             "residual forces).")
    parser.add_argument("--max-steps", type=int, default=200,
                        help="Max optimizer steps for both relax stages (default: 200).")
    parser.add_argument("--max", type=float, default=1.0,
                        help="Max strain magnitude, %% (default: 1.0 -- smaller than "
                             "stb-elasticInputs' 2.0 default, to stay safely in the linear "
                             "regime given MACE's own accuracy limits).")
    parser.add_argument("--steps", type=int, default=4,
                        help="Strain points per direction, spanning [-max, +max] (default: 4).")
    parser.add_argument("--check-convergence", action="store_true",
                        help="Instead of a single fit, re-run the pipeline at multiple --max "
                             "strain magnitudes (--convergence-strains) and report how each "
                             "diagonal constant changes -- checks whether the linear-regime "
                             "assumption actually holds at the strain used. No foundation-"
                             "model comparison in this mode. Off by default.")
    parser.add_argument("--convergence-strains", type=float, nargs='+', default=None,
                        metavar="PCT",
                        help="Max strain magnitudes (%%) to compare for --check-convergence "
                             "(default: half, the given --max, and double it).")
    parser.add_argument("--dimensionality", choices=["auto", "3d", "2d", "1d"], default="auto",
                        help="Physical dimensionality, controlling report units and which "
                             "constants are physically meaningful (default: auto-detected from "
                             "vacuum padding). Manual override assumes the same axis convention "
                             "as stb-elasticAnalysis: in-plane xy for 2D, wire-along-c for 1D.")
    parser.add_argument("--vacuum-gap", type=float, default=10.0,
                        help="Gap (Ang) used to detect vacuum-padded axes (default: 10.0).")
    parser.add_argument("--fit-quality-tolerance", type=float, default=0.995,
                        help="Minimum R^2 for a direction's stress-strain fit before it's "
                             "flagged as low-quality (default: 0.995).")
    parser.add_argument("--skip-anisotropy-surface", action="store_true",
                        help="Skip the directional Young's-modulus anisotropy surface plot "
                             "(3D only) -- it's an O(n_theta*n_phi) tensor contraction, cheap, "
                             "but this opts out if not wanted. Off by default.")
    parser.add_argument("--device", choices=["cpu", "cuda"], default="cpu",
                        help="Device to run the model on (default: cpu).")
    parser.add_argument("-o", "--output-dir", default="mlelastic_out",
                        help="Output directory for the plot/data (default: mlelastic_out).")
    parser.add_argument("--save-data", action="store_true",
                        help="Also write the raw stress-strain data (.dat) behind the plot. Off by default.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-mlelastic {VERSION}")

    args = parser.parse_args()

    if args.check_convergence and args.convergence_strains is not None \
            and len(args.convergence_strains) < 2:
        parser.error("--convergence-strains needs at least 2 values to compare.")

    # emit_elastic_report (imported from elastic_analysis.py) is generic
    # over any caller's argparse Namespace as long as it carries these exact
    # attribute names -- 'method'/'symmetry_method' only matter there for a
    # DFT-specific [INFO] branch that never triggers here (filled_by_symmetry
    # is always empty below, since we never use symmetry-based direction
    # reduction/pooling), and 'plot_dir' just labels our own --output-dir.
    args.method = "stress"
    args.symmetry_method = "basic"
    args.eggbox_tolerance = 0.0
    args.plot_dir = args.output_dir

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - ML Elastic Constants",
            "Stress-strain stiffness matrix, driven entirely by a MACE potential",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("ML ELASTIC CONSTANTS:", 'bold'))
    print("-" * 60)

    os.makedirs(args.output_dir, exist_ok=True)
    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    def fail(message):
        print_dual(color_text(f"[FAIL] {message}", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_dual(color_text("===== STB-MLELASTIC REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Input file        : {args.file}", f_out)
    mode_label = "strain-amplitude convergence check" if args.check_convergence else "single fit"
    print_dual(f"Mode              : {mode_label}", f_out)
    print_dual(f"Model             : "
               f"{'custom (' + args.custom_model + ')' if args.custom_model else f'MACE-MP-0 ({args.model})'}", f_out)
    if not args.check_convergence:
        print_dual(f"Strain range      : +/-{args.max}% ({args.steps} points/direction)", f_out)
    compare = bool(args.custom_model) and not args.skip_foundation_comparison \
        and not args.check_convergence
    if compare:
        print_dual(f"Comparison        : also evaluating with foundation model ({args.model})", f_out)
    if report_path:
        print_dual(f"Report file       : {report_path}", f_out)

    print_section("[1] READING STRUCTURE", f_out)

    if not os.path.isfile(args.file):
        fail(f"Input file '{args.file}' not found.")
    if args.custom_model and not os.path.isfile(args.custom_model):
        fail(f"--custom-model file not found: {args.custom_model}")

    try:
        structure = structure_io.read_fdf(args.file)
    except Exception as e:
        fail(f"Could not read '{args.file}': {e}")

    pmg_structure = structure_io.to_pymatgen(structure)
    atoms = AseAtomsAdaptor.get_atoms(pmg_structure)

    lattice_ang = np.array(atoms.get_cell())
    frac_coords = atoms.get_scaled_positions()
    vacuum_axes = kspace.detect_vacuum_axes(frac_coords, lattice_ang, args.vacuum_gap)
    vacuum_count = sum(vacuum_axes)
    print_dual(f"[OK] Read {len(atoms)} atom(s), "
               f"dimensionality: {kspace.dimensionality_label(vacuum_axes)}", f_out)

    if args.dimensionality != "auto":
        dimensionality = args.dimensionality
        print_dual(f"Dimensionality    : {dimensionality.upper()} (manual override)", f_out)
    elif vacuum_count == 3:
        fail("Detected a 0D (isolated molecule) structure -- elastic constants aren't "
             "defined for an isolated molecule.")
    elif vacuum_count == 0:
        dimensionality = "3d"
    elif vacuum_count == 1 and vacuum_axes == [False, False, True]:
        dimensionality = "2d"
    elif vacuum_count == 2 and vacuum_axes == [True, True, False]:
        dimensionality = "1d"
    else:
        print_dual(color_text(
            "[WARNING] Vacuum padding detected on an axis this tool doesn't support (expects "
            "in-plane xy for 2D, or wire-along-c for 1D) -- falling back to 3D (GPa). Pass "
            "--dimensionality explicitly to override.", 'yellow'), f_out)
        dimensionality = "3d"
    args.dimensionality = dimensionality
    unit_label = UNIT_LABELS[dimensionality]
    print_dual(f"Mode              : {dimensionality.upper()} (units: {unit_label})", f_out)

    periodic_axes = {i for i in range(3) if not vacuum_axes[i]}

    model_arg = args.custom_model if args.custom_model else args.model
    calc = mace_relax.get_calculator(model=model_arg, device=args.device)

    if args.check_convergence:
        strains = args.convergence_strains if args.convergence_strains else \
            [args.max / 2.0, args.max, args.max * 2.0]

        print_section("[2] STRAIN-AMPLITUDE CONVERGENCE CHECK", f_out)
        rows = []
        for max_strain in strains:
            sub_args = copy.copy(args)
            sub_args.max = max_strain
            tag = f"max={max_strain:.2f}%"
            print_dual(f"\n--- {tag} ---", f_out)
            result = run_elastic_pipeline(atoms, calc, sub_args, vacuum_axes, periodic_axes,
                                           dimensionality, f_out, tag=tag)
            rows.append((max_strain, result))

        modes_present = rows[0][1]["dirs_run"] if rows else []
        print_section("[3] SUMMARY & FILES", f_out)
        print_dual(f"{'MaxStrain%':>12} " + " ".join(f"{('C_'+m):>12}" for m in modes_present), f_out)
        for max_strain, result in rows:
            vals = " ".join(f"{result['C_sym'][DIRECTION_TO_COLUMN[m], DIRECTION_TO_COLUMN[m]]:12.4f}"
                             for m in modes_present)
            print_dual(f"{max_strain:12.2f} {vals}", f_out)

        conv_plot = os.path.join(args.output_dir, "strain_convergence.png")
        plot_strain_convergence(rows, modes_present, unit_label, conv_plot)
        print_dual(f"[OK] Saved {conv_plot}", f_out)
        if args.save_data:
            conv_data = os.path.join(args.output_dir, "strain_convergence.dat")
            with open(conv_data, "w") as f:
                f.write("# max_strain(%) " + " ".join(f"C_{m}({unit_label})" for m in modes_present) + "\n")
                for max_strain, result in rows:
                    vals = " ".join(
                        f"{result['C_sym'][DIRECTION_TO_COLUMN[m], DIRECTION_TO_COLUMN[m]]:.6f}"
                        for m in modes_present)
                    f.write(f"{max_strain:.4f} {vals}\n")
            print_dual(f"[OK] Saved {conv_data}", f_out)
        print_dual("Status            : OK", f_out)
        print_dual(f"Convergence plot  : {conv_plot}", f_out)
        if report_path:
            print_dual(f"Full report       : {report_path}", f_out)
        if f_out:
            f_out.close()
        print("\n" + "-" * 60)
        print(color_text("ML elastic strain-convergence check complete.\n", 'bold'))
        return

    print_section("[2] REFERENCE RELAXATION & STRESS-STRAIN DATA (MACE)", f_out)
    result_main = run_elastic_pipeline(
        atoms, calc, args, vacuum_axes, periodic_axes, dimensionality, f_out,
        tag="fine-tuned" if args.custom_model else args.model)
    if not result_main["dirs_run"]:
        fail("No direction is physically valid for this structure (all touch a vacuum-padded axis).")

    result_found = None
    if compare:
        print_dual(f"\n[INFO] Also running the raw foundation model ({args.model}) for comparison...", f_out)
        found_calc = mace_relax.get_calculator(model=args.model, device=args.device)
        result_found = run_elastic_pipeline(
            atoms, found_calc, args, vacuum_axes, periodic_axes, dimensionality, f_out,
            tag=f"foundation ({args.model})")

    plot_path = os.path.join(args.output_dir, "stress_strain.png")
    series = [(("Fine-tuned" if args.custom_model else f"MACE-MP-0 ({args.model})"),
               result_main["data"], result_main["fit_diagnostics"], result_main["conv_factor"],
               _SERIES_COLORS[0])]
    if result_found is not None:
        series.append((f"Foundation ({args.model})", result_found["data"],
                        result_found["fit_diagnostics"], result_found["conv_factor"],
                        _SERIES_COLORS[1]))
    plotted = plot_stress_strain(series, unit_label, plot_path)
    if plotted:
        print_dual(f"[OK] Saved {plot_path}", f_out)
    if args.save_data:
        data_path = os.path.join(args.output_dir, "stress_strain.dat")
        with open(data_path, "w") as f:
            f.write("# mode strain(frac) eng_strain(frac) sigma_xx sigma_yy sigma_zz "
                    f"sigma_yz sigma_xz sigma_xy [{unit_label}]\n")
            for mode in result_main["dirs_run"]:
                d = result_main["data"][mode]
                eng_factor = 2.0 if mode in SHEAR_DIRECTIONS else 1.0
                for eps, s in zip(d['eps'], d['stress']):
                    row = [result_main["conv_factor"] * s[i, j] for (i, j) in core_symmetry.VOIGT_TO_TENSOR]
                    f.write(f"{mode:>4} {eps:12.6f} {eps*eng_factor:12.6f} "
                            + " ".join(f"{v:12.6f}" for v in row) + "\n")
        print_dual(f"[OK] Saved {data_path}", f_out)

    emit_elastic_report(args, result_main["C_sym"], unit_label, [], result_main["missing_directions"],
                         f_out, fit_diagnostics=result_main["fit_diagnostics"],
                         symmetry_check=result_main["symmetry_check"], plot_files=None,
                         report_path=report_path)
    if plotted:
        print_dual(f"Plot              : {plot_path}", f_out)

    if result_found is not None:
        print_section("[6] FOUNDATION MODEL COMPARISON", f_out)
        modes_present = [m for m in core_symmetry.VOIGT_MODES if m in result_main["dirs_run"]]
        print_dual(f"{'':>10} " + " ".join(f"{('C_'+m):>10}" for m in modes_present), f_out)
        for label, result in (("Fine-tuned", result_main), (f"Foundation ({args.model})", result_found)):
            vals = " ".join(f"{result['C_sym'][DIRECTION_TO_COLUMN[m], DIRECTION_TO_COLUMN[m]]:10.2f}"
                             for m in modes_present)
            print_dual(f"{label:>10} {vals}", f_out)

    if dimensionality == "3d" and result_main["missing_directions"]:
        print_section("[7] SOUND VELOCITIES & DEBYE TEMPERATURE", f_out)
        print_dual(color_text(
            f"[WARNING] Skipped (and [8] DIRECTIONAL ANISOTROPY below) -- direction(s) "
            f"{', '.join(result_main['missing_directions'])} have no data (vacuum-padded axis, "
            "or --dimensionality forced 3D on a lower-dimensional structure). The isotropic "
            "3D-average formulas need a complete, non-degenerate 6-direction stiffness "
            "tensor -- computing them from a rank-deficient one would silently give wrong "
            "numbers, not just imprecise ones.", 'yellow'), f_out)
    elif dimensionality == "3d":
        print_section("[7] SOUND VELOCITIES & DEBYE TEMPERATURE", f_out)
        sound = compute_sound_velocities(result_main["C_sym"], result_main["atoms_relaxed"])
        if sound is None:
            print_dual(color_text(
                "[WARNING] Could not compute sound velocities/Debye temperature (non-positive "
                "bulk/shear modulus -- same instability check_stability_and_report already "
                "flagged).", 'yellow'), f_out)
        else:
            print_dual(f"Density (rho)     : {sound['rho']:.2f} kg/m^3", f_out)
            print_dual(f"Longitudinal (Vp) : {sound['v_p']:.2f} m/s", f_out)
            print_dual(f"Shear (Vs)        : {sound['v_s']:.2f} m/s", f_out)
            print_dual(f"Average (Vm)      : {sound['v_m']:.2f} m/s", f_out)
            theta_d_str = f"{sound['theta_D']:.1f}"
            print_dual(f"Debye temperature : {color_text(theta_d_str, 'bold')} K", f_out)

        anisotropy_plot = None
        if not args.skip_anisotropy_surface:
            print_section("[8] DIRECTIONAL ANISOTROPY", f_out)
            candidate_plot = os.path.join(args.output_dir, "anisotropy_surface.png")
            e_min, e_max = plot_anisotropy_surface(result_main["C_sym"], candidate_plot)
            if e_min is None:
                print_dual(color_text(
                    "[WARNING] Could not compute the anisotropy surface (singular stiffness "
                    "matrix).", 'yellow'), f_out)
            else:
                anisotropy_plot = candidate_plot
                print_dual(f"Young's modulus range: {e_min:.2f} - {e_max:.2f} {unit_label}", f_out)
                if e_min > 0:
                    print_dual(f"Anisotropy ratio (E_max/E_min): {e_max/e_min:.3f} "
                               "(1.0 = perfectly isotropic)", f_out)
                print_dual(f"[OK] Saved {anisotropy_plot}", f_out)

        print_section("[9] SUMMARY & FILES (ML-SPECIFIC)", f_out)
        if sound is not None:
            print_dual(f"Vp={sound['v_p']:.0f} m/s  Vs={sound['v_s']:.0f} m/s  "
                       f"theta_D={sound['theta_D']:.1f} K", f_out)
        if anisotropy_plot:
            print_dual(f"Anisotropy surface: {anisotropy_plot}", f_out)

    if f_out:
        f_out.close()

    print("\n" + "-" * 60)
    print(color_text("ML elastic constants complete -- a fast heuristic, not a substitute "
                      "for DFT elastic constants.\n", 'bold'))


if __name__ == "__main__":
    main()
