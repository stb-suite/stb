#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""stb-mlelastic (7.4): elastic constants from a MACE potential, by the stress-strain method.

The ML Simulations analog of the SIESTA-based stb-elasticInputs/stb-elasticAnalysis pair (as stb-mlphonons is of
5.4): every canonical Voigt strain not blocked by a vacuum axis, at --steps points in +-(--max) %, the atoms relaxed
at each fixed strained cell, the model's stress at each point, a linear fit per direction. One self-contained run,
no folders: a MACE stress costs milliseconds. The linear algebra of the fit is elastic_analysis.py's
(compute_stiffness_matrix, direction_fit_diagnostics, tensor_symmetry_check) and the strain convention
elastic_inputs.py's (get_strain_matrix, direction_axes), so 'xy' means the same deformation in both tools.

On top of the fit the report gives: each model's reference state (relaxed cell, residual stress); a per-direction
table (fit quality, residual stress, internal relaxation); the strained cells against the model card's trained
range; the tensor projected on the crystal's point group (which constants are independent, how far the raw fit is
from that symmetry); Born stability; Voigt/Reuss/Hill moduli and the usual ratios; clamped-ion vs relaxed-ion
constants; Christoffel sound velocities along [100]/[110]/[111] and the Debye temperature; the directional Young's
modulus and linear compressibility; fine-tuned vs foundation. Two modes besides the single fit: --check-convergence
(the fit at several strain amplitudes). Sections numbered in order, every warning in the last one, library messages
only counted (same report style as stb-mlmd/stb-mlphonons).
"""

from stb import __version__ as VERSION

import os
import sys
import copy
import time
import argparse
from datetime import datetime

import numpy as np

from stb.core import structure_io, kspace, model_card
from stb.core.cli import (color_text, show_intro, print_dual, print_table,
                          capture_library_noise, display_available, Sections, run_gnuplot)

REPORT_FILE = "stb_mlelastic_report.txt"
PLOT_DIR = "plot"
EV_A3_TO_GPA = 160.21766208
AMU_A3_TO_G_CM3 = 1.66053906660
HBAR, KB, AMU = 1.054571817e-34, 1.380649e-23, 1.66053906660e-27
VOIGT = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
VNAME = ["11", "22", "33", "44", "55", "66"]
SYM_TOL = 0.01                 # relative deviation of the raw fit from the point-group projection that warns
_COLORS = ["#2a6fdb", "#d9480f", "#2f9e44", "#862e9c"]
SOUND_DIRS = (("[100]", (1, 0, 0)), ("[110]", (1, 1, 0)), ("[111]", (1, 1, 1)))




def fmt(v, spec="{:.2f}", none="-"):
    return none if v is None or (isinstance(v, float) and not np.isfinite(v)) else spec.format(v)


# ---- tensors ---------------------------------------------------------------------

def voigt_to_tensor(C):
    """6x6 Voigt STIFFNESS (no factors) -> c_ijkl."""
    c = np.zeros((3, 3, 3, 3))
    for m, (i, j) in enumerate(VOIGT):
        for n, (k, l) in enumerate(VOIGT):
            for a, b in ((i, j), (j, i)):
                for p, q in ((k, l), (l, k)):
                    c[a, b, p, q] = C[m, n]
    return c


def tensor_to_voigt(c):
    return np.array([[c[i, j, k, l] for (k, l) in VOIGT] for (i, j) in VOIGT])


def compliance_voigt_to_full(S):
    """6x6 Voigt COMPLIANCE (engineering shear) -> s_ijkl: factor 1/2 per shear index."""
    s = np.zeros((3, 3, 3, 3))
    for m, (i, j) in enumerate(VOIGT):
        for n, (k, l) in enumerate(VOIGT):
            val = S[m, n] / ((2.0 if m >= 3 else 1.0) * (2.0 if n >= 3 else 1.0))
            for a, b in ((i, j), (j, i)):
                for p, q in ((k, l), (l, k)):
                    s[a, b, p, q] = val
    return s


def point_group(atoms, symprec=1e-3):
    """(space group, point group, crystal system, [Cartesian rotations]) of a cell, from spglib; None if it fails."""
    try:
        import spglib
        cell = np.array(atoms.get_cell())
        ds = spglib.get_symmetry_dataset((cell, atoms.get_scaled_positions(), atoms.numbers), symprec=symprec)
        if ds is None:
            return None
        sg = ds.international if hasattr(ds, "international") else ds["international"]
        num = ds.number if hasattr(ds, "number") else ds["number"]
        pg = ds.pointgroup if hasattr(ds, "pointgroup") else ds["pointgroup"]
        rots = ds.rotations if hasattr(ds, "rotations") else ds["rotations"]
        uniq, cart = set(), []
        for R in rots:
            key = R.tobytes()
            if key in uniq:
                continue
            uniq.add(key)
            Rc = cell.T @ R @ np.linalg.inv(cell.T)
            # a cell that is symmetric only within the tolerance gives a slightly non-orthogonal Rc: use the nearest
            # rotation (polar decomposition), so the projection makes symmetry-equal entries equal
            u, _sv, vt = np.linalg.svd(Rc)
            cart.append(u @ vt)
        system = ("triclinic" if num <= 2 else "monoclinic" if num <= 15 else "orthorhombic" if num <= 74 else
                  "tetragonal" if num <= 142 else "trigonal" if num <= 167 else "hexagonal" if num <= 194 else "cubic")
        return f"{sg} ({num})", pg, system, cart
    except Exception:
        return None


def project_on_point_group(C, rotations):
    """Average of C over the point group, in the cell's own Cartesian frame: the tensor with exactly the symmetry
    of the crystal. Works for any orientation of the cell."""
    c = voigt_to_tensor(C)
    acc = np.zeros_like(c)
    for R in rotations:
        acc += np.einsum("ia,jb,kc,ld,abcd->ijkl", R, R, R, R, c, optimize=True)
    return tensor_to_voigt(acc / len(rotations))


def independent_constants(P, tol, same=None):
    """Groups the non-zero upper-triangle entries of a symmetry-projected tensor by value: [(names, value)].
    Entries below `tol` are zero; entries equal within `same` (default 1e-4 of the largest) are one constant -- the
    projection makes symmetry-equal entries equal to ~1e-6 (a cell symmetric only within spglib's tolerance), while a
    looser `same` would merge entries that are only numerically close (a P1 cell of a nearly cubic crystal, whose
    entries differ by ~1e-3, once came out as '3 independent constants')."""
    if same is None:
        same = 1e-4 * max(1e-12, float(np.abs(P).max()))
    groups = []
    for m in range(6):
        for n in range(m, 6):
            v = P[m, n]
            if abs(v) <= tol:
                continue
            name = f"C{m + 1}{n + 1}"
            for g in groups:
                if abs(g[1] - v) <= same:
                    g[0].append(name)
                    break
            else:
                groups.append([[name], v])
    return groups


def vrh(C):
    """Voigt, Reuss, Hill bulk and shear moduli and the derived E, nu (3D). Reuss = None if C is singular."""
    Bv = (C[0, 0] + C[1, 1] + C[2, 2] + 2 * (C[0, 1] + C[0, 2] + C[1, 2])) / 9.0
    Gv = ((C[0, 0] + C[1, 1] + C[2, 2]) - (C[0, 1] + C[0, 2] + C[1, 2]) + 3 * (C[3, 3] + C[4, 4] + C[5, 5])) / 15.0
    out = {"V": (Bv, Gv)}
    try:
        S = np.linalg.inv(C)
        Br = 1.0 / ((S[0, 0] + S[1, 1] + S[2, 2]) + 2 * (S[0, 1] + S[0, 2] + S[1, 2]))
        Gr = 15.0 / (4 * (S[0, 0] + S[1, 1] + S[2, 2]) - 4 * (S[0, 1] + S[0, 2] + S[1, 2])
                     + 3 * (S[3, 3] + S[4, 4] + S[5, 5]))
        out["R"] = (Br, Gr)
        out["H"] = ((Bv + Br) / 2, (Gv + Gr) / 2)
    except np.linalg.LinAlgError:
        out["R"] = out["H"] = None
    res = {}
    for k, v in out.items():
        if v is None:
            res[k] = None
            continue
        B, G = v
        E = 9 * B * G / (3 * B + G) if (3 * B + G) else float("nan")
        nu = (3 * B - 2 * G) / (2 * (3 * B + G)) if (3 * B + G) else float("nan")
        res[k] = {"B": B, "G": G, "E": E, "nu": nu}
    return res


def christoffel(C, rho, n):
    """Sound velocities (m/s, sorted) along unit direction n: eigenvalues of Gamma_ik = c_ijkl n_j n_l / rho."""
    c = voigt_to_tensor(C) * 1e9
    n = np.asarray(n, float) / np.linalg.norm(n)
    G = np.einsum("ijkl,j,l->ik", c, n, n) / rho
    ev = np.linalg.eigvalsh(G)
    return np.sqrt(np.clip(ev, 0, None)), ev.min() < 0


def debye_from_vrh(B, G, rho, n_density):
    """Isotropic (Hill) longitudinal, transverse and mean velocities and the Debye temperature."""
    vl = np.sqrt((B + 4.0 / 3.0 * G) * 1e9 / rho)
    vt = np.sqrt(G * 1e9 / rho)
    vm = (1.0 / 3.0 * (2.0 / vt ** 3 + 1.0 / vl ** 3)) ** (-1.0 / 3.0)
    theta = HBAR / KB * (6 * np.pi ** 2 * n_density) ** (1.0 / 3.0) * vm
    return vl, vt, vm, theta


def directional(C, n):
    """Young's modulus E(n) and linear compressibility beta(n) (1/TPa) along unit direction n (3D)."""
    s = compliance_voigt_to_full(np.linalg.inv(C))
    n = np.asarray(n, float) / np.linalg.norm(n)
    inv_e = np.einsum("ijkl,i,j,k,l->", s, n, n, n, n)
    beta = np.einsum("ijkk,i,j->", s, n, n) * 1000.0
    return (1.0 / inv_e if inv_e > 0 else float("nan")), beta


def directional_grid(C, n_theta=61, n_phi=121):
    """E and beta on a sphere; returns thetas, phis, E, beta (or None if C is singular)."""
    try:
        s = compliance_voigt_to_full(np.linalg.inv(C))
    except np.linalg.LinAlgError:
        return None
    th = np.linspace(0, np.pi, n_theta)
    ph = np.linspace(0, 2 * np.pi, n_phi)
    T, P = np.meshgrid(th, ph, indexing="ij")
    N = np.stack([np.sin(T) * np.cos(P), np.sin(T) * np.sin(P), np.cos(T)], axis=-1)
    inv_e = np.einsum("ijkl,...i,...j,...k,...l->...", s, N, N, N, N, optimize=True)
    E = np.where(inv_e > 0, 1.0 / np.where(inv_e > 0, inv_e, 1), np.nan)
    beta = np.einsum("ijkk,...i,...j->...", s, N, N, optimize=True) * 1000.0
    return th, ph, E, beta, N


def youngs_2d(C, theta):
    """In-plane Young's modulus of a 2D sheet along angle theta (from C11, C22, C12, C66)."""
    K = np.array([[C[0, 0], C[0, 1], C[0, 5]], [C[0, 1], C[1, 1], C[1, 5]], [C[0, 5], C[1, 5], C[5, 5]]])
    S = np.linalg.inv(K)
    c, s = np.cos(theta), np.sin(theta)
    inv = S[0, 0] * c ** 4 + S[1, 1] * s ** 4 + (2 * S[0, 1] + S[2, 2]) * c ** 2 * s ** 2 \
        + 2 * S[0, 2] * c ** 3 * s + 2 * S[1, 2] * c * s ** 3
    return 1.0 / inv


# ---- the calculation ---------------------------------------------------------------

def strain_scan(atoms_ref, calc, args, periodic_axes, relax, label, noise):
    """Strains the reference cell along every Voigt direction not blocked by a vacuum axis, at --steps points in
    +-(--max) %, relaxes the atoms at fixed strained cell (if `relax`), records the stress. Returns (data,
    dirs_run, skipped, points) -- data in elastic_analysis' format, points one dict per strained cell (volume per
    atom, pressure, max |stress|, shortest distance, internal relaxation, relax steps)."""
    from stb.core import symmetry as core_symmetry, mace_relax
    from stb.elastic_inputs import get_strain_matrix, direction_axes
    dirs_run, skipped = [], []
    for d in core_symmetry.VOIGT_MODES:
        (skipped if direction_axes(d) - periodic_axes else dirs_run).append(d)
    strains = np.linspace(-args.max, args.max, args.steps)
    ref_cell = np.array(atoms_ref.get_cell())
    frac = atoms_ref.get_scaled_positions()
    data, points = {}, []
    for d in dirs_run:
        eps_l, st_l = [], []
        for s_ in strains:
            delta = 0.0 if abs(s_) < 1e-9 else s_ / 100.0
            cell = ref_cell @ get_strain_matrix(d, delta)
            a = atoms_ref.copy()
            a.set_cell(cell, scale_atoms=False)
            a.set_scaled_positions(frac)
            a.calc = calc
            pos0 = a.get_positions().copy()
            steps = 0
            with capture_library_noise(noise, label):
                if relax and len(a) > 1:
                    _conv, steps = mace_relax.relax(a, calc, cell_mask=None, fmax=args.fmax, max_steps=args.max_steps)
                stress = a.get_stress(voigt=False)
            df = (a.get_positions() - pos0) @ np.linalg.inv(cell)
            df -= np.round(df)
            shift = float(np.linalg.norm(df @ cell, axis=1).max())
            eps_l.append(delta)
            st_l.append(stress)
            points.append({"mode": d, "strain": delta, "volume": a.get_volume() / len(a),
                           "pressure": -np.trace(stress) / 3 * EV_A3_TO_GPA,
                           "max_stress": float(np.abs(stress).max() * EV_A3_TO_GPA),
                           "d_min": model_card.shortest_distance(a), "shift": shift, "steps": steps})
        data[d] = {"eps": eps_l, "stress": st_l}
    return data, dirs_run, skipped, points


def run_model(atoms_in, calc, args, vacuum_axes, periodic_axes, dimensionality, label, noise, clamped=False):
    """One model: pre-relax (cell + ions), the strain scan (relaxed ions, and clamped ions too when `clamped`),
    the fit and its diagnostics. Returns a dict with everything the report needs."""
    from stb.core import mace_relax
    from stb.elastic_analysis import (compute_stiffness_matrix, direction_fit_diagnostics, tensor_symmetry_check,
                                      UNIT_LABELS, CONV_EVA3_TO_GPA, CONV_EVA2_TO_NM, CONV_EVA1_TO_NN)
    t0 = time.monotonic()
    atoms = atoms_in.copy()
    atoms.calc = calc
    res = {"label": label}
    with capture_library_noise(noise, label):
        cell0 = np.array(atoms.get_cell())
        v0 = atoms.get_volume()
        stress0 = atoms.get_stress(voigt=False) * EV_A3_TO_GPA
        f0 = float(np.linalg.norm(atoms.get_forces(), axis=1).max())
        if args.pre_relax:
            converged, steps = mace_relax.relax(atoms, calc, cell_mask=mace_relax.build_cell_mask(vacuum_axes),
                                                fmax=args.fmax, max_steps=args.max_steps)
        else:
            converged, steps = None, 0
        f1 = float(np.linalg.norm(atoms.get_forces(), axis=1).max())
        stress1 = atoms.get_stress(voigt=False) * EV_A3_TO_GPA
    res.update(f0=f0, f1=f1, converged=converged, relax_steps=steps, cell0=cell0, cell1=np.array(atoms.get_cell()),
               v0=v0, v1=atoms.get_volume(), stress0=stress0, stress1=stress1, atoms=atoms)
    cell = res["cell1"]
    if dimensionality == "2d":
        conv = cell[2, 2] * CONV_EVA2_TO_NM
        res["norm_note"] = f"cell height Lz = {cell[2, 2]:.4f} Ang (stress x Lz: per-area stiffness)"
    elif dimensionality == "1d":
        xs = np.linalg.norm(cell[0]) * np.linalg.norm(cell[1])
        conv = xs * CONV_EVA1_TO_NN
        res["norm_note"] = f"cross-section Lx*Ly = {xs:.4f} Ang^2 (stress x area: axial force)"
    else:
        conv = CONV_EVA3_TO_GPA
        res["norm_note"] = None
    res["conv"], res["unit"] = conv, UNIT_LABELS[dimensionality]

    data, dirs, skipped, points = strain_scan(atoms, calc, args, periodic_axes, args.relax, label, noise)
    res.update(data=data, dirs=dirs, skipped=skipped, points=points)
    sink = open(os.devnull, "w")
    with capture_library_noise(noise, label):
        C, _filled, missing = compute_stiffness_matrix(data, conv, None, res["unit"], sink)
    res["C_raw"], res["C"], res["missing"] = C, 0.5 * (C + C.T), missing
    res["fit"] = direction_fit_diagnostics(data, conv, args.fit_quality_tolerance)
    res["tsym"] = tensor_symmetry_check(C)
    res["n_eval"] = len(points)
    if clamped:
        data_c, _d, _s, points_c = strain_scan(atoms, calc, args, periodic_axes, False, label, noise)
        with capture_library_noise(noise, label):
            Cc, _f, _m = compute_stiffness_matrix(data_c, conv, None, res["unit"], sink)
        res["C_clamped"] = 0.5 * (Cc + Cc.T)
        res["n_eval"] += len(points_c)
    sink.close()
    pg = point_group(atoms)
    res["pg"] = pg
    if pg is not None:
        res["C_proj"] = project_on_point_group(res["C"], pg[3])
    res["seconds"] = time.monotonic() - t0
    return res


# ---- derived quantities -----------------------------------------------------------

def stability(C, dimensionality, system):
    """[(criterion, value text, passes)], overall verdict. Necessary and sufficient: all eigenvalues of C > 0
    (Mouhat & Coudert 2014); the crystal-system conditions are its explicit form."""
    rows = []
    if dimensionality == "1d":
        rows.append(("C33 > 0", f"{C[2, 2]:.2f}", C[2, 2] > 0))
        return rows, C[2, 2] > 0
    if dimensionality == "2d":
        c11, c22, c12, c66 = C[0, 0], C[1, 1], C[0, 1], C[5, 5]
        rows += [("C11 > 0", f"{c11:.2f}", c11 > 0), ("C66 > 0", f"{c66:.2f}", c66 > 0),
                 ("C11 C22 - C12^2 > 0", f"{c11 * c22 - c12 ** 2:.2f}", c11 * c22 - c12 ** 2 > 0)]
        return rows, all(r[2] for r in rows)
    ev = np.linalg.eigvalsh(C)
    rows.append(("smallest eigenvalue of C > 0", f"{ev.min():.2f}", ev.min() > 0))
    if system == "cubic":
        rows += [("C11 - C12 > 0", f"{C[0, 0] - C[0, 1]:.2f}", C[0, 0] - C[0, 1] > 0),
                 ("C11 + 2 C12 > 0", f"{C[0, 0] + 2 * C[0, 1]:.2f}", C[0, 0] + 2 * C[0, 1] > 0),
                 ("C44 > 0", f"{C[3, 3]:.2f}", C[3, 3] > 0)]
    elif system in ("hexagonal", "tetragonal"):
        rows += [("C11 > |C12|", f"{C[0, 0] - abs(C[0, 1]):.2f}", C[0, 0] > abs(C[0, 1])),
                 ("2 C13^2 < C33 (C11 + C12)", f"{C[2, 2] * (C[0, 0] + C[0, 1]) - 2 * C[0, 2] ** 2:.2f}",
                  2 * C[0, 2] ** 2 < C[2, 2] * (C[0, 0] + C[0, 1])),
                 ("C44 > 0", f"{C[3, 3]:.2f}", C[3, 3] > 0)]
        if system == "tetragonal":
            rows.append(("C66 > 0", f"{C[5, 5]:.2f}", C[5, 5] > 0))
    return rows, ev.min() > 0


# ---- plots and data ------------------------------------------------------------------

def plot_stress_strain(series, unit, path):
    """One panel per direction; data and fit of every model."""
    import matplotlib.pyplot as plt
    from stb.elastic_analysis import DIRECTION_TO_COLUMN, SHEAR_DIRECTIONS
    from stb.core import symmetry as core_symmetry
    modes = [m for m in core_symmetry.VOIGT_MODES if m in series[0]["data"]]
    if not modes:
        return None
    nc = min(3, len(modes))
    nr = int(np.ceil(len(modes) / nc))
    fig, axes = plt.subplots(nr, nc, figsize=(4.4 * nc, 3.5 * nr), squeeze=False)
    for k, mode in enumerate(modes):
        ax = axes[k // nc][k % nc]
        i, j = VOIGT[DIRECTION_TO_COLUMN[mode]]
        ef = 2.0 if mode in SHEAR_DIRECTIONS else 1.0
        for s_, res in enumerate(series):
            if mode not in res["data"]:
                continue
            fit = {r["mode"]: r for r in res["fit"]}.get(mode)
            x = np.array(res["data"][mode]["eps"]) * ef * 100
            y = res["conv"] * np.array(res["data"][mode]["stress"])[:, i, j]
            c = _COLORS[s_ % 4]
            ax.plot(x, y, "o", color=c, ms=5, label=f"{res['label']}")
            if fit and x.max() > x.min():
                xs = np.linspace(x.min(), x.max(), 20)
                ax.plot(xs, fit["slope"] * xs / 100 + fit["residual_stress"], "--", color=c, lw=1.2,
                        label=f"C = {fit['slope']:.1f}")
        ax.set_title(f"{mode}  (C{VNAME[DIRECTION_TO_COLUMN[mode]]})", fontsize=10)
        ax.set_xlabel("engineering strain (%)")
        ax.set_ylabel(f"stress σ_{mode} ({unit})")
        ax.axhline(0, color="grey", lw=0.5)
        ax.legend(fontsize=7)
    for k in range(len(modes), nr * nc):
        axes[k // nc][k % nc].axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    return fig


def plot_directional(series, path):
    """E(n) and beta(n) in two Cartesian planes (xy; z and x+y), polar, every model."""
    import matplotlib.pyplot as plt
    t = np.linspace(0, 2 * np.pi, 361)
    # Cartesian planes of the cell's own frame (for a cubic cell they are (001) and (1-10); for a hexagonal cell with c
    # along z, the basal plane and a plane containing c)
    planes = (("xy plane: x at 0°, y at 90°", np.stack([np.cos(t), np.sin(t), 0 * t], 1)),
              ("plane of z and x+y: z at 0°, x+y at 90°", np.stack([np.sin(t) / np.sqrt(2), np.sin(t) / np.sqrt(2),
                                                                    np.cos(t)], 1)))
    fig, axes = plt.subplots(1, 4, figsize=(16, 5.0), subplot_kw={"projection": "polar"})
    for p, (title, N) in enumerate(planes):
        for s_, res in enumerate(series):
            vals = np.array([directional(res["C"], n) for n in N])
            axes[2 * p].plot(t, vals[:, 0], color=_COLORS[s_ % 4], lw=1.8, label=res["label"])
            axes[2 * p + 1].plot(t, vals[:, 1], color=_COLORS[s_ % 4], lw=1.8, label=res["label"])
        axes[2 * p].set_title(f"Young's modulus E (GPa)\n{title}", fontsize=9, pad=14)
        axes[2 * p + 1].set_title(f"linear compressibility β (1/TPa)\n{title}", fontsize=9, pad=14)
    for ax in axes:                     # room around the curve: a constant beta would sit on the rim
        top = max(np.nanmax(line.get_ydata()) for line in ax.get_lines())
        ax.set_ylim(0, 1.12 * top)
    h, lab = axes[0].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=len(series), fontsize=9, frameon=False)
    fig.tight_layout(rect=(0, 0.06, 1, 0.97))
    fig.savefig(path, dpi=150)
    return fig


def plot_surface(C, path, label):
    import matplotlib.pyplot as plt
    g = directional_grid(C)
    if g is None:
        return None
    th, ph, E, _b, N = g
    X, Y, Z = (E[..., None] * N).transpose(2, 0, 1)
    fig = plt.figure(figsize=(6.5, 5.6))
    ax = fig.add_subplot(111, projection="3d")
    norm = plt.Normalize(np.nanmin(E), np.nanmax(E))
    ax.plot_surface(X, Y, Z, facecolors=plt.cm.viridis(norm(E)), rstride=1, cstride=1, shade=False, antialiased=True)
    for setter, lab in ((ax.set_xlabel, "E_x"), (ax.set_ylabel, "E_y"), (ax.set_zlabel, "E_z")):
        setter(f"{lab} (GPa)")
    ax.set_title(f"Directional Young's modulus, {label}\nmin {np.nanmin(E):.1f}, max {np.nanmax(E):.1f} GPa", fontsize=10)
    m = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
    m.set_array(E)
    fig.colorbar(m, ax=ax, shrink=0.6, label="E (GPa)")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    return fig


def plot_inplane(series, unit, path):
    import matplotlib.pyplot as plt
    t = np.linspace(0, 2 * np.pi, 361)
    fig, ax = plt.subplots(figsize=(5.5, 5), subplot_kw={"projection": "polar"})
    for s_, res in enumerate(series):
        ax.plot(t, youngs_2d(res["C"], t), color=_COLORS[s_ % 4], lw=1.8, label=res["label"])
    ax.set_ylim(0, 1.12 * max(np.nanmax(line.get_ydata()) for line in ax.get_lines()))
    ax.set_title(f"In-plane Young's modulus ({unit}); a hexagonal sheet is isotropic in-plane", fontsize=10)
    ax.legend(fontsize=7, loc="lower left", bbox_to_anchor=(-0.15, -0.15))
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    return fig


def plot_convergence(rows, groups_names, unit, path):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    x = [r["max"] for r in rows]
    for names in groups_names:
        m, n = int(names[0][1]) - 1, int(names[0][2]) - 1
        axes[0].plot(x, [r["res"]["C"][m, n] for r in rows], "o-", label=" = ".join(names[:3]))
    axes[0].set_xlabel("maximum strain (%)")
    axes[0].set_ylabel(f"constant ({unit})")
    axes[0].set_title("independent constants against the strain amplitude")
    axes[0].legend(fontsize=7)
    ref = rows[0]["res"]["C"]
    for names in groups_names:
        m, n = int(names[0][1]) - 1, int(names[0][2]) - 1
        axes[1].plot(x, [100 * (r["res"]["C"][m, n] / ref[m, n] - 1) for r in rows], "o-",
                     label=" = ".join(names[:3]))
    axes[1].axhline(0, color="grey", lw=0.6)
    axes[1].set_xlabel("maximum strain (%)")
    axes[1].set_ylabel(f"change vs {x[0]:g} % (%)")
    axes[1].set_title("relative change")
    axes[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    return fig


def write_stress_strain_dat(path, res):
    from stb.elastic_analysis import SHEAR_DIRECTIONS
    with open(path, "w") as f:
        f.write(f"# {res['label']}: mode strain(frac) eng_strain(frac) sigma_xx sigma_yy sigma_zz sigma_yz sigma_xz "
                f"sigma_xy [{res['unit']}]\n")
        for mode in res["dirs"]:
            d = res["data"][mode]
            ef = 2.0 if mode in SHEAR_DIRECTIONS else 1.0
            for eps, s_ in zip(d["eps"], d["stress"]):
                f.write(f"{mode:>4} {eps:12.6f} {eps * ef:12.6f} "
                        + " ".join(f"{res['conv'] * s_[i, j]:12.6f}" for (i, j) in VOIGT) + "\n")
            f.write("\n\n")


def write_tensor_dat(path, res):
    with open(path, "w") as f:
        f.write(f"# {res['label']}: stiffness C ({res['unit']}), Voigt order 11 22 33 23 13 12\n")
        for row in res["C"]:
            f.write(" ".join(f"{v:12.4f}" for v in row) + "\n")
        if "C_proj" in res:
            f.write("\n# projected on the point group\n")
            for row in res["C_proj"]:
                f.write(" ".join(f"{v:12.4f}" for v in row) + "\n")
        try:
            S = np.linalg.inv(res["C"])
            f.write(f"\n# compliance S (1/{res['unit']})\n")
            for row in S:
                f.write(" ".join(f"{v:14.6e}" for v in row) + "\n")
        except np.linalg.LinAlgError:
            pass


def write_directional_dat(path, series):
    t = np.linspace(0, 2 * np.pi, 361)
    with open(path, "w") as f:
        f.write("# angle(deg) then per model: E_xy(GPa) beta_xy(1/TPa) E_z,x+y(GPa) beta_z,x+y(1/TPa) (Cartesian planes) -- models: "
                + ", ".join(r["label"] for r in series) + "\n")
        for a in t:
            n1 = (np.cos(a), np.sin(a), 0)
            n2 = (np.sin(a) / np.sqrt(2), np.sin(a) / np.sqrt(2), np.cos(a))
            vals = []
            for r in series:
                e1, b1 = directional(r["C"], n1)
                e2, b2 = directional(r["C"], n2)
                vals += [e1, b1, e2, b2]
            f.write(f"{np.degrees(a):8.2f} " + " ".join(f"{v:12.5f}" for v in vals) + "\n")


def gnuplot_files(plot_dir, series, three_d):
    """stress_strain.gplot (one panel per direction) and, in 3D, directional.gplot. Returns the script names."""
    from stb.elastic_analysis import SHEAR_DIRECTIONS, DIRECTION_TO_COLUMN
    os.makedirs(plot_dir, exist_ok=True)
    gp_colors = ["#2a6fdb", "#d9480f", "#2f9e44", "#862e9c"]
    for k, res in enumerate(series):
        write_stress_strain_dat(os.path.join(plot_dir, f"stress_strain_{k}.dat"), res)
    modes = series[0]["dirs"]
    nc = min(3, len(modes))
    nr = int(np.ceil(len(modes) / nc))
    L = ["# --- STB Plot Configuration ---\n# Generated by stb-mlelastic: gnuplot stress_strain.gplot\n",
         f'set terminal pdfcairo enhanced color font "Arial,9" size {4.2 * nc:.1f},{3.4 * nr:.1f}\n',
         'set output "stress_strain.pdf"\nset grid\nset key top left font ",7"\n',
         f'set multiplot layout {nr},{nc}\nset xlabel "engineering strain (%)"\n']
    for b, mode in enumerate(modes):
        i, j = VOIGT[DIRECTION_TO_COLUMN[mode]]
        col = 4 + DIRECTION_TO_COLUMN[mode]
        L.append(f'set title "{mode} (C{VNAME[DIRECTION_TO_COLUMN[mode]]})"\n'
                 f'set ylabel "stress ({series[0]["unit"]})"\n')
        parts = []
        for k, res in enumerate(series):
            fit = {r["mode"]: r for r in res["fit"]}.get(mode)
            parts.append(f'"stress_strain_{k}.dat" index {b} using ($3*100):{col} with points pt 7 ps 0.6 '
                         f'lc rgb "{gp_colors[k % 4]}" title "{res["label"]}"')
            if fit:
                parts.append(f'{fit["slope"]:.6f}*x/100 + {fit["residual_stress"]:.6f} with lines dt 2 '
                             f'lc rgb "{gp_colors[k % 4]}" title "C = {fit["slope"]:.1f}"')
        L.append("plot " + ", \\\n     ".join(parts) + "\n")
    L.append("unset multiplot\n")
    with open(os.path.join(plot_dir, "stress_strain.gplot"), "w") as f:
        f.writelines(L)
    scripts = ["stress_strain.gplot"]
    if three_d:
        write_directional_dat(os.path.join(plot_dir, "directional.dat"), series)
        L = ["# --- STB Plot Configuration ---\n# Generated by stb-mlelastic: gnuplot directional.gplot\n",
             'set terminal pdfcairo enhanced color font "Arial,9" size 11,5.5\nset output "directional.pdf"\n',
             "set polar\nset angles degrees\nset size square\nset grid polar 30\nunset border\nunset xtics\n"
             "unset ytics\nset key bottom right font \",7\"\nset multiplot layout 1,2\n"]
        for p, title in ((0, "Young's modulus E (GPa), xy plane"), (2, "Young's modulus E (GPa), plane of z and x+y")):
            L.append(f'set title "{title}"\n')
            parts = [f'"directional.dat" using 1:{2 + 4 * k + p} with lines lw 2 lc rgb "{gp_colors[k % 4]}" '
                     f'title "{res["label"]}"' for k, res in enumerate(series)]
            L.append("plot " + ", \\\n     ".join(parts) + "\n")
        L.append("unset multiplot\n")
        with open(os.path.join(plot_dir, "directional.gplot"), "w") as f:
            f.writelines(L)
        scripts.append("directional.gplot")
    return scripts




# ---- report pieces -------------------------------------------------------------------

def report_reference(res, f_out, warn, dimensionality):
    a0, a1 = np.linalg.norm(res["cell0"], axis=1), np.linalg.norm(res["cell1"], axis=1)
    print_dual(f"--- {res['label']}", f_out)
    if res["converged"] is None:
        print_dual(f"Pre-relax         : skipped (--no-pre-relax); max |F| {res['f0']:.4f} eV/Ang", f_out)
    else:
        print_dual(f"Pre-relax         : cell + atoms, max |F| {res['f0']:.4f} -> {res['f1']:.4f} eV/Ang in "
                   f"{res['relax_steps']} step(s) ({'converged' if res['converged'] else 'NOT converged: step cap'}, "
                   f"--fmax)", f_out)
        if not res["converged"]:
            warn(f"{res['label']}: the pre-relax hit --max-steps without converging; the reference is not at a "
                 "minimum and every constant inherits its residual stress.")
    print_dual(f"Lattice a, b, c   : {a0[0]:.5f} {a0[1]:.5f} {a0[2]:.5f} -> {a1[0]:.5f} {a1[1]:.5f} {a1[2]:.5f} Ang "
               f"(volume {100 * (res['v1'] / res['v0'] - 1):+.3f} %)", f_out)
    p0, p1 = -np.trace(res["stress0"]) / 3, -np.trace(res["stress1"]) / 3
    print_dual(f"Stress            : pressure {p0:+.3f} -> {p1:+.3f} GPa; largest component after "
               f"{np.abs(res['stress1']).max():.3f} GPa", f_out)
    if res["norm_note"]:
        print_dual(f"Normalization     : {res['norm_note']}", f_out)
    if abs(res["v1"] / res["v0"] - 1) > 0.01 and dimensionality == "3d":
        warn(f"{res['label']}: its own equilibrium volume is {100 * (res['v1'] / res['v0'] - 1):+.2f} % from the input "
             "structure; the constants are at the model's minimum, not at the input lattice.")
    if np.abs(res["stress1"]).max() > 0.2:
        why = ("--no-pre-relax: the input cell is not stress-free" if res["converged"] is None
               else "the pre-relax did not reach a stress-free cell")
        warn(f"{res['label']}: {np.abs(res['stress1']).max():.2f} GPa of stress remains at zero strain ({why}); it "
             "appears as the fits' intercepts.")


def report_directions(res, f_out, warn, args):
    from stb.elastic_analysis import DIRECTION_TO_COLUMN
    fits = {r["mode"]: r for r in res["fit"]}
    rows = []
    for d in res["dirs"]:
        f = fits.get(d)
        pts = [p for p in res["points"] if p["mode"] == d]
        smax = max(abs(p["strain"]) for p in pts) * 100
        shift = max(p["shift"] for p in pts)
        rows.append(([d, f"C{VNAME[DIRECTION_TO_COLUMN[d]]}", f"{f['slope']:.2f}" if f else "-",
                      f"{f['r_squared']:.6f}" if f else "-", f"{f['residual_stress']:+.3f}" if f else "-",
                      f"{max(p['max_stress'] for p in pts):.2f}", f"{shift:.4f}",
                      f"{shift / smax:.4f}" if smax else "-", str(sum(p['steps'] for p in pts)),
                      "check" if f and f["flagged"] else "ok"], "yellow" if f and f["flagged"] else None))

    flagged = [r["mode"] for r in res["fit"] if r["flagged"]]
    if flagged:
        worst = min(r["r_squared"] for r in res["fit"] if r["flagged"])
        warn(f"{res['label']}: fit(s) {', '.join(flagged)} flagged (lowest R^2 {worst:.4f}, --fit-quality-tolerance "
             f"{args.fit_quality_tolerance}; or an intercept above 10 % of the stress at the largest strain): check the "
             "stress-strain plot, a tighter --fmax, the pre-relax, and --check-convergence.")
    print_dual(f"--- {res['label']}: {len(res['dirs'])} directions x {args.steps} strains, {res['n_eval']} relaxed "
               f"cells, {res['seconds']:.1f} s", f_out)
    print_table(["dir.", "column", f"diagonal ({res['unit']})", "R^2", "intercept", "max|σ| GPa", "max shift Ang",
                 "Ang per %", "relax steps", "fit"], rows, f_out)
    if res["skipped"]:
        print_dual(f"Skipped           : {', '.join(res['skipped'])} -- touch a vacuum axis (no physical deformation)",
                   f_out)


def report_domain(series, card, f_out, warn, dimensionality):
    dom = (card or {}).get("domain") or {}
    rows = []
    for res in series:
        pts = res["points"]
        vols = [p["volume"] for p in pts]
        press = [p["pressure"] for p in pts]
        dmin = min(p["d_min"] for p in pts if p["d_min"] is not None) if any(p["d_min"] for p in pts) else None
        rows.append((res["label"], min(vols), max(vols), min(press), max(press), dmin))
    print_dual("Strained cells (all directions, relaxed atoms):", f_out)
    print_table(["model", "V/atom min", "V/atom max", "pressure min GPa", "pressure max GPa", "shortest dist. Ang"],
                [([r[0], f"{r[1]:.4f}", f"{r[2]:.4f}", f"{r[3]:+.3f}", f"{r[4]:+.3f}", fmt(r[5], "{:.4f}")], None)
                 for r in rows], f_out)
    if not card:
        print_dual("No model card (foundation model): no trained range to compare with.", f_out)
        return
    r = rows[0]
    m, pr, dd = dom.get("measure"), dom.get("pressure_GPa"), dom.get("d_min_ang")
    tab = []
    if m and dom.get("measure_name", "").startswith("volume") and dimensionality == "3d":
        ok = m["min"] <= r[1] and r[2] <= m["max"]
        tab.append(([dom["measure_name"], f"{m['min']:.3f} .. {m['max']:.3f}", f"{r[1]:.3f} .. {r[2]:.3f}",
                     "inside" if ok else "OUTSIDE"], None if ok else "yellow"))
        if not ok:
            warn(f"Model domain: the strained cells' volume per atom ({r[1]:.3f}-{r[2]:.3f}) leaves the trained range "
                 f"{m['min']:.3f}-{m['max']:.3f}.")
    if pr:
        ok = pr["min"] <= r[3] and r[4] <= pr["max"]
        tab.append((["pressure (GPa)", f"{pr['min']:.2f} .. {pr['max']:.2f}", f"{r[3]:+.2f} .. {r[4]:+.2f}",
                     "inside" if ok else "OUTSIDE"], None if ok else "yellow"))
        if not ok:
            warn(f"Model domain: the pressure along the strain scan ({r[3]:+.2f} .. {r[4]:+.2f} GPa) leaves the trained "
                 f"range {pr['min']:.2f} .. {pr['max']:.2f} GPa.")
    if dd and r[5] is not None:
        ok = r[5] >= dd["min"]
        tab.append((["shortest distance (Ang)", f">= {dd['min']:.3f}", f"min {r[5]:.3f}", "inside" if ok else "OUTSIDE"],
                    None if ok else "yellow"))
    if tab:
        print_dual(f"\n{res_label(series[0])} against its model card's trained range:", f_out)
        print_table(["quantity", "trained", "this run", "verdict"], tab, f_out)


def res_label(res):
    return res["label"]


def print_matrix(M, unit, f_out, title, fmtv="{:9.2f}"):
    print_dual(f"{title} ({unit}), Voigt order 1=xx 2=yy 3=zz 4=yz 5=xz 6=xy:", f_out)
    print_dual("        " + "".join(f"{'j=' + str(j + 1):>10}" for j in range(6)), f_out)
    for i in range(6):
        print_dual(f"  i={i + 1} |" + "".join(f" {fmtv.format(v)}" for v in M[i]), f_out)


def report_tensor(res, f_out, warn, dimensionality):
    C = res["C"]
    print_matrix(C, res["unit"], f_out, f"Stiffness C, {res['label']}")
    mx, i, j, flagged = res["tsym"]
    scale = np.abs(C).max()
    print_dual(f"\nC vs C^T          : max |C_ij - C_ji| = {mx:.3f} {res['unit']} (C{i + 1}{j + 1}) "
               f"[{'CHECK' if flagged else 'OK'}] -- reported C is (C + C^T)/2", f_out)
    if flagged:
        warn(f"{res['label']}: C is not symmetric (|C{i + 1}{j + 1} - C{j + 1}{i + 1}| = {mx:.2f} {res['unit']}, above 1 % "
             "of the largest constant): the two directions disagree; check the pre-relax and the fits.")
    if res.get("pg") is not None and "C_proj" in res:
        sg, pgs, system, rots = res["pg"]
        P = res["C_proj"]
        dev = np.abs(C - P).max()
        print_dual(f"Crystal symmetry  : {sg}, point group {pgs} ({system}, {len(rots)} rotations, spglib on the "
                   "pre-relaxed cell)", f_out)
        print_dual(f"Projection        : C averaged over the point group; max |C - C_sym| = {dev:.3f} {res['unit']} "
                   f"({100 * dev / scale:.2f} % of the largest constant)", f_out)
        if dev > SYM_TOL * scale:
            warn(f"{res['label']}: the fitted tensor departs from the crystal's {pgs} symmetry by {100 * dev / scale:.1f} % "
                 "of its largest constant -- noisy fits, an unrelaxed reference, or a structure whose symmetry spglib "
                 "found too high or too low.")
        groups = independent_constants(P, 1e-3 * scale)
        if len(rots) == 1:
            print_dual(f"Independent       : no symmetry found (point group {pgs}): all 21 constants are free "
                       "(the matrix above)", f_out)
            warn(f"{res['label']}: spglib finds no symmetry in the pre-relaxed cell (point group {pgs}): no projection, "
                 "no symmetry check. A structure that should be symmetric was not relaxed to its symmetric minimum "
                 "(pre-relax off, or stopped early).")
        else:
            # a linear relation that is not an equality: hexagonal (and trigonal, transversely isotropic) symmetry ties
            # C66 to (C11 - C12)/2 -- the projection obeys it exactly, so it is not an independent constant
            dependent = set()
            first = {g[0][0]: g[1] for g in groups}
            if "C11" in first and "C12" in first:
                half = (first["C11"] - first["C12"]) / 2
                for k, (names, v) in enumerate(groups):
                    if names[0] == "C66" and abs(v - half) <= 1e-4 * scale:
                        dependent.add(k)
            print_dual(f"Independent       : {len(groups) - len(dependent)} non-zero independent constant(s) under "
                       f"{pgs}:", f_out)
            for k, (names, v) in enumerate(groups):
                tail = "  = (C11 - C12)/2, not independent" if k in dependent else ""
                print_dual(f"  {' = '.join(names):<36} {v:9.2f} {res['unit']}{tail}", f_out)
    if dimensionality == "3d":
        try:
            print_dual("", f_out)
            print_matrix(np.linalg.inv(C) * 1000.0, "1/TPa", f_out, "Compliance S = C^-1", "{:9.3f}")
        except np.linalg.LinAlgError:
            warn(f"{res['label']}: the stiffness matrix is singular; no compliance, moduli or directional properties.")


def report_moduli(res, f_out, warn, system):
    C = res["C"]
    v = vrh(C)
    rows = []
    for k, name in (("V", "Voigt (upper bound)"), ("R", "Reuss (lower bound)"), ("H", "Hill (average)")):
        if v.get(k) is None:
            continue
        x = v[k]
        rows.append(([name, f"{x['B']:.2f}", f"{x['G']:.2f}", f"{x['E']:.2f}", f"{x['nu']:.3f}"], None))
    print_table(["average", "B (GPa)", "G (GPa)", "E (GPa)", "Poisson ν"], rows, f_out)
    h, vv, r = v.get("H"), v["V"], v.get("R")
    out = {"vrh": v}
    if h:
        k = h["G"] / h["B"]
        print_dual(f"\nPugh ratio B/G    : {h['B'] / h['G']:.3f} (above ~1.75: ductile, below: brittle; Pugh 1954)", f_out)
        if system == "cubic":
            cp = C[0, 1] - C[3, 3]
            za = 2 * C[3, 3] / (C[0, 0] - C[0, 1]) if C[0, 0] != C[0, 1] else float("nan")
            print_dual(f"Cauchy pressure   : C12 - C44 = {cp:+.2f} GPa (negative: directional, covalent bonding; "
                       "positive: metallic)", f_out)
            print_dual(f"Zener ratio       : A = 2 C44 / (C11 - C12) = {za:.3f} (1 = isotropic)", f_out)
            out.update(cauchy=cp, zener=za)
        if r:
            au = 5 * vv["G"] / r["G"] + vv["B"] / r["B"] - 6
            print_dual(f"Universal anisot. : A^U = 5 G_V/G_R + B_V/B_R - 6 = {au:.3f} (0 = isotropic; "
                       "Ranganathan & Ostoja-Starzewski 2008)", f_out)
            out["AU"] = au
        hv = 2 * (k * k * h["G"]) ** 0.585 - 3
        print_dual(f"Vickers hardness  : {hv:.1f} GPa (empirical, Chen et al. 2011: 2 (k^2 G)^0.585 - 3, k = G/B)", f_out)
        out.update(pugh=h["B"] / h["G"], hv=hv)
        try:
            S = np.linalg.inv(C)
            beta = [1000 * S[i, :3].sum() for i in range(3)]
            print_dual(f"Linear compress.  : β_x {beta[0]:.3f}, β_y {beta[1]:.3f}, β_z {beta[2]:.3f} 1/TPa "
                       f"(volume: {sum(beta):.3f} = 1000/B_R)", f_out)
            out["beta"] = beta
        except np.linalg.LinAlgError:
            pass
        if h["G"] <= 0 or h["B"] <= 0:
            warn(f"{res['label']}: non-positive Hill modulus -- the structure is mechanically unstable.")
    return out


def report_2d(res, f_out):
    C = res["C"]
    c11, c22, c12, c66 = C[0, 0], C[1, 1], C[0, 1], C[5, 5]
    den = c11 * c22 - c12 ** 2
    u = res["unit"]
    rows = [(["Young's modulus E_x = (C11 C22 - C12^2)/C22", f"{den / c22:.2f} {u}" if c22 else "-"], None),
            (["Young's modulus E_y = (C11 C22 - C12^2)/C11", f"{den / c11:.2f} {u}" if c11 else "-"], None),
            (["Poisson ratio ν_xy = C12/C22", f"{c12 / c22:.3f}" if c22 else "-"], None),
            (["Poisson ratio ν_yx = C12/C11", f"{c12 / c11:.3f}" if c11 else "-"], None),
            (["shear modulus G_xy = C66", f"{c66:.2f} {u}"], None),
            (["layer modulus γ = (C11 + C22 + 2 C12)/4", f"{(c11 + c22 + 2 * c12) / 4:.2f} {u}"], None)]
    print_table(["quantity", "value"], rows, f_out)
    try:
        t = np.linspace(0, 2 * np.pi, 721)
        e = youngs_2d(C, t)
        print_dual(f"In-plane E(θ)     : {e.min():.2f} .. {e.max():.2f} {u} (ratio {e.max() / e.min():.3f}; "
                   "1 = in-plane isotropic, as for a hexagonal sheet)", f_out)
    except np.linalg.LinAlgError:
        pass


def report_clamped(res, f_out, warn):
    Cr, Cc = res["C"], res["C_clamped"]
    rows = []
    groups = independent_constants(res.get("C_proj", Cr), 1e-3 * np.abs(Cr).max())
    for names, _v in groups:
        m, n = int(names[0][1]) - 1, int(names[0][2]) - 1
        d = 100 * (Cr[m, n] / Cc[m, n] - 1) if abs(Cc[m, n]) > 1e-6 else float("nan")
        rows.append(([" = ".join(names[:3]) + (" ..." if len(names) > 3 else ""), f"{Cc[m, n]:.2f}", f"{Cr[m, n]:.2f}",
                      f"{d:+.1f} %"], None))
    print_dual("Clamped ions: the atoms follow the strain homogeneously (no relaxation inside the cell).\n"
               "Relaxed ions: the atoms relax at each strained cell -- the physical, measured constants.", f_out)
    print_table(["constant", f"clamped ({res['unit']})", f"relaxed ({res['unit']})", "relaxation"], rows, f_out)
    big = [r[0][0] for r in rows if r[0][3] != "nan %" and abs(float(r[0][3].rstrip(" %"))) > 5]
    if big:
        print_dual(f"Internal relaxation matters for {', '.join(big)}: a cell with 2+ atoms per primitive cell whose "
                   "atoms are not all at inversion centres (diamond, zincblende, wurtzite...) relaxes inside the cell "
                   "under shear.", f_out)


def report_sound(res, f_out, warn):
    C = res["C"]
    atoms = res["atoms"]
    rho = atoms.get_masses().sum() * AMU / (atoms.get_volume() * 1e-30)
    rows = []
    for name, n in SOUND_DIRS:
        v, neg = christoffel(C, rho, n)
        rows.append(([name, f"{v[2]:.0f}", f"{v[1]:.0f}", f"{v[0]:.0f}"], "yellow" if neg else None))
        if neg:
            warn(f"{res['label']}: an imaginary sound velocity along {name}: the tensor is unstable for that wave.")
    print_dual(f"Density           : {rho:.2f} kg/m^3 ({rho / 1000:.4f} g/cm^3) at the pre-relaxed cell", f_out)
    print_dual("Christoffel velocities along Cartesian directions (m/s): the speeds of the three sound waves, "
               "eigenvalues of c_ijkl n_j n_l / rho:", f_out)
    print_table(["direction", "quasi-longitudinal", "fast transverse", "slow transverse"], rows, f_out)
    h = vrh(C).get("H")
    if h and h["B"] > 0 and h["G"] > 0:
        vl, vt, vm, th = debye_from_vrh(h["B"], h["G"], rho, len(atoms) / (atoms.get_volume() * 1e-30))
        print_dual(f"Isotropic (Hill)  : v_l {vl:.0f}, v_t {vt:.0f}, mean v_m {vm:.0f} m/s "
                   "(v_m^-3 = (v_l^-3 + 2 v_t^-3)/3)", f_out)
        print_dual(f"Debye temperature : {color_text(f'{th:.1f}', 'bold')} K "
                   "(θ_D = ħ/k_B (6π² n)^(1/3) v_m, n = atoms per volume; Anderson 1963)", f_out)
        return {"rho": rho, "vl": vl, "vt": vt, "vm": vm, "theta": th}
    return {"rho": rho}


def report_directional(res, f_out):
    g = directional_grid(res["C"])
    if g is None:
        return None
    th, ph, E, beta, N = g
    i_min = np.unravel_index(np.nanargmin(E), E.shape)
    i_max = np.unravel_index(np.nanargmax(E), E.shape)

    def nice(n):
        n = np.abs(n) / np.abs(n).max()
        for cand, lab in (((1, 0, 0), "<100>"), ((1, 1, 0), "<110>"), ((1, 1, 1), "<111>")):
            cs = np.array(sorted(n, reverse=True))
            c = np.array(cand, float)
            if np.allclose(cs / np.linalg.norm(cs), c / np.linalg.norm(c), atol=0.03):
                return lab
        return "(" + " ".join(f"{x:.2f}" for x in n) + ")"
    rows = [(["Young's modulus E (GPa)", f"{np.nanmin(E):.2f} along {nice(N[i_min])}",
              f"{np.nanmax(E):.2f} along {nice(N[i_max])}", f"{np.nanmax(E) / np.nanmin(E):.3f}"], None),
            (["linear compressibility β (1/TPa)", f"{beta.min():.3f}", f"{beta.max():.3f}",
              f"{beta.max() / beta.min():.3f}" if beta.min() > 0 else "-"], None)]
    vals = [directional(res["C"], n) for _, n in SOUND_DIRS]
    print_table(["property", "minimum", "maximum", "max/min"], rows, f_out)
    print_dual("Along the Cartesian axes: " + ", ".join(f"E{name} {e:.2f} GPa" for (name, _), (e, _b) in
                                                       zip(SOUND_DIRS, vals)), f_out)
    return {"E_min": float(np.nanmin(E)), "E_max": float(np.nanmax(E))}


def report_comparison(series, f_out, system):
    a, b = series[0], series[1]
    groups = independent_constants(a.get("C_proj", a["C"]), 1e-3 * np.abs(a["C"]).max())
    rows = []
    for names, _v in groups:
        m, n = int(names[0][1]) - 1, int(names[0][2]) - 1
        x, y = a["C"][m, n], b["C"][m, n]
        rows.append(([" = ".join(names[:3]) + (" ..." if len(names) > 3 else ""), f"{x:.2f}", f"{y:.2f}",
                      f"{100 * (y / x - 1):+.1f} %" if abs(x) > 1e-9 else "-"], None))
    va, vb = vrh(a["C"]).get("H"), vrh(b["C"]).get("H")
    if va and vb:
        for key, name in (("B", "B Hill (GPa)"), ("G", "G Hill (GPa)"), ("E", "E Hill (GPa)"), ("nu", "Poisson ν")):
            rows.append(([name, f"{va[key]:.3f}" if key == "nu" else f"{va[key]:.2f}",
                          f"{vb[key]:.3f}" if key == "nu" else f"{vb[key]:.2f}",
                          f"{100 * (vb[key] / va[key] - 1):+.1f} %"], None))
    la = np.linalg.norm(a["cell1"], axis=1)
    lb = np.linalg.norm(b["cell1"], axis=1)
    rows.append((["lattice a (Ang, own minimum)", f"{la[0]:.5f}", f"{lb[0]:.5f}", f"{100 * (lb[0] / la[0] - 1):+.2f} %"],
                 None))
    print_table(["quantity", a["label"], b["label"], "foundation vs fine-tuned"], rows, f_out)
    print_dual("Each model is strained about its own relaxed cell. A large difference is what fine-tuning taught the "
               "model: the foundation model was never told about this level of theory.", f_out)


# ---- parser -------------------------------------------------------------------------

def build_parser():
    parser = argparse.ArgumentParser(
        description="Elastic constants by the stress-strain method, driven entirely by a MACE potential -- the "
                    "MACE-MP-0 foundation model, or a custom model fine-tuned on your own SIESTA data via "
                    "stb-mlffAnalysis (--custom-model, always with its model card). Strains the relaxed cell along "
                    "every canonical Voigt direction not blocked by a vacuum axis, relaxes the atoms at each strained "
                    "cell, fits the stiffness tensor and reports: each model's reference state, fit quality per "
                    "direction, the strained cells against the card's trained range, the tensor projected on the "
                    "crystal's point group, Born stability, Voigt/Reuss/Hill moduli, clamped- vs relaxed-ion "
                    "constants, Christoffel sound velocities and the Debye temperature, directional Young's modulus "
                    "and compressibility, and fine-tuned vs foundation; every warning in one final section. Or, with "
                    "--check-convergence, the fit at several strain amplitudes. A fast preview, not a substitute for "
                    "the SIESTA-based stb-elasticInputs/stb-elasticAnalysis workflow.",
        epilog="Example usage:\n"
               "  stb-mlelastic --file structure.fdf\n"
               "  stb-mlelastic --file structure.fdf --custom-model mlff_model.model --custom-model-card mlff_model_card.json "
               "--device cuda --steps 7 --save-report\n"
               "  stb-mlelastic --file structure.fdf --check-convergence --convergence-strains 0.25 0.5 1 2 3\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--file", required=True, help="Input structure (.fdf).")
    parser.add_argument("--model", choices=["small", "medium", "large"], default="medium",
                        help="MACE-MP-0 model size (default: medium). Ignored if --custom-model is given.")
    parser.add_argument("--custom-model", default=None, metavar="PATH",
                        help="Use a custom MACE model file instead of the MACE-MP-0 foundation "
                             "potential -- e.g. one fine-tuned via stb-mlffAnalysis. Overrides --model.")
    model_card.add_card_argument(parser, "--custom-model")
    parser.add_argument("--skip-foundation-comparison", action="store_true",
                        help="With --custom-model the same calculation is repeated with the raw foundation model "
                             "(--model) and compared constant by constant; this flag skips it. --check-convergence "
                             "never compares.")
    parser.add_argument("--no-pre-relax", dest="pre_relax", action="store_false",
                        help="Skip the pre-relax of the reference structure (cell + atoms). On by default: the "
                             "constants are second derivatives about a stress-free minimum; a residual stress shows "
                             "up as the fits' intercepts.")
    parser.add_argument("--no-relax", dest="relax", action="store_false",
                        help="Do not relax the atoms at each strained cell (clamped-ion constants only). On by "
                             "default: the measured constants are relaxed-ion ones (diamond's C44 depends on it).")
    parser.add_argument("--no-clamped", dest="clamped", action="store_false",
                        help="Skip the clamped-ion scan that is compared with the relaxed-ion constants (it doubles "
                             "the number of stress evaluations; run only when the cell has 2+ atoms).")
    parser.add_argument("--fmax", type=float, default=0.001,
                        help="Force threshold of both relaxations, eV/Ang (default: 0.001). Tight on purpose: under "
                             "shear the atoms of a diamond-type cell relax by only ~0.008 Ang per %% strain, and a "
                             "residual force f leaves them ~f/k off (k ~ 13 eV/Ang^2): 0.02 made silicon's C44 3 %% "
                             "low and noisy, 0.001 is converged to 0.3 %%. With MACE the extra steps cost "
                             "milliseconds.")
    parser.add_argument("--max-steps", type=int, default=200,
                        help="Max optimizer steps of both relaxations (default: 200).")
    parser.add_argument("--max", type=float, default=1.0,
                        help="Max strain magnitude, %% (default: 1.0).")
    parser.add_argument("--steps", type=int, default=5,
                        help="Strain points per direction, spanning [-max, +max] (default: 5, which includes the "
                             "unstrained cell).")
    parser.add_argument("--check-convergence", action="store_true",
                        help="Instead of one fit, repeat it at several --max strain magnitudes "
                             "(--convergence-strains) and compare every independent constant: does the linear "
                             "regime hold? No foundation comparison in this mode.")
    parser.add_argument("--convergence-strains", type=float, nargs='+', default=None, metavar="PCT",
                        help="Max strain magnitudes (%%) for --check-convergence (default: half, the given --max, "
                             "and double it).")
    parser.add_argument("--dimensionality", choices=["auto", "3d", "2d", "1d"], default="auto",
                        help="Physical dimensionality, controlling report units and which constants are physically "
                             "meaningful (default: auto-detected from vacuum padding). Manual override assumes the "
                             "same axis convention as stb-elasticAnalysis: in-plane xy for 2D, wire-along-c for 1D.")
    parser.add_argument("--vacuum-gap", type=float, default=10.0,
                        help="Gap (Ang) used to detect vacuum-padded axes (default: 10.0).")
    parser.add_argument("--fit-quality-tolerance", type=float, default=0.995,
                        help="Minimum R^2 of a direction's stress-strain fit before it is flagged (default: 0.995).")
    parser.add_argument("--skip-anisotropy-surface", action="store_true",
                        help="Do not draw the 3D Young's-modulus surface (anisotropy_surface.png). The directional "
                             "numbers and the polar plots are kept.")
    parser.add_argument("--device", choices=["cpu", "cuda"], default="cpu",
                        help="Device to run the model on (default: cpu).")
    parser.add_argument("-o", "--output-dir", default="mlelastic_out",
                        help="Output directory (default: mlelastic_out).")
    parser.add_argument("--show", action="store_true",
                        help="Open the figures on screen at the end (needs a display; without one it is a warning). "
                             "The PNGs are written either way.")
    parser.add_argument("--save-gnuplot", action="store_true",
                        help=f"Write .dat + .gplot into <output-dir>/{PLOT_DIR}/ and run gnuplot there (one PDF per "
                             "figure) when it is installed. Off by default.")
    parser.add_argument("--save-data", action="store_true",
                        help="Write the data behind the report: stress_strain.dat (every strained cell's stress), "
                             "elastic_tensor.dat (C, its projection, S), directional.dat. Off by default.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-mlelastic {VERSION}")
    return parser


# ---- main -----------------------------------------------------------------------------

def main():
    parser = build_parser()
    args = parser.parse_args()
    card = model_card.setup_from_args(args.custom_model, args.custom_model_card, "--custom-model", error=parser.error)
    if args.check_convergence and args.convergence_strains is not None and len(args.convergence_strains) < 2:
        parser.error("--convergence-strains needs at least 2 values to compare.")
    if args.steps < 2:
        parser.error("--steps must be at least 2 (a fit needs two strains).")
    if args.max <= 0:
        parser.error("--max must be positive.")

    import matplotlib
    if not args.show:
        matplotlib.use("Agg")

    if args.intro:
        show_intro(["Siesta ToolBox Suite - ML Elastic Constants",
                    "Stress-strain stiffness tensor, driven entirely by a MACE potential",
                    f"Version {VERSION} | University of Brasilia - 2026",
                    "Developed by Dr. Carlos M. O. Bastos"])

    os.makedirs(args.output_dir, exist_ok=True)
    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None
    warns, noise, figures, files = [], [], [], []
    t_start = time.monotonic()
    sections = Sections(f_out)

    def warn(message):
        warns.append(message)

    def fail(message):
        print_dual(color_text(f"[FAIL] {message}", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    import importlib.util
    from stb.core.deps import require_mace
    if importlib.util.find_spec("mace") is None or importlib.util.find_spec("torch") is None:
        require_mace()
    with capture_library_noise(noise, "import"):
        require_mace()
        from pymatgen.io.ase import AseAtomsAdaptor
        from stb.core import mace_relax
        from stb.elastic_analysis import UNIT_LABELS
        import torch
        import ase
        import mace
    model_card.collect_warnings(warns)

    print("\n" + color_text("ML ELASTIC CONSTANTS:", 'bold'))
    print("-" * 60)
    print_dual(color_text("===== STB-MLELASTIC REPORT =====", 'magenta'), f_out)

    compare = bool(args.custom_model) and not args.skip_foundation_comparison and not args.check_convergence
    sections("RUN METADATA")
    print_dual(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Software          : stb {VERSION}, mace-torch {getattr(mace, '__version__', '?')}, torch "
               f"{torch.__version__}, ase {ase.__version__}", f_out)
    if args.device == "cuda":
        try:
            mace_relax.resolve_device("cuda")
            dev_label = f"cuda ({torch.cuda.get_device_name(0)})"
        except Exception as e:
            fail(f"--device cuda: {e}")
    else:
        dev_label = "cpu"
    print_dual(f"Device            : {dev_label}", f_out)
    print_dual(f"Input file        : {args.file}", f_out)
    print_dual(f"Mode              : {'strain-amplitude convergence check' if args.check_convergence else 'single fit'}",
               f_out)
    if not args.check_convergence:
        print_dual(f"Strains           : {args.steps} points in +-{args.max:g} % per direction "
                   f"({', '.join(f'{x:+.3g}' for x in np.linspace(-args.max, args.max, args.steps))} %)", f_out)
    print_dual(f"Relaxation        : reference {'cell + atoms' if args.pre_relax else 'none (--no-pre-relax)'}; strained "
               f"cells {'atoms relaxed (relaxed-ion constants)' if args.relax else 'atoms clamped (--no-relax)'}; "
               f"--fmax {args.fmax} eV/Ang, --max-steps {args.max_steps}", f_out)
    print_dual(f"Model             : "
               f"{'custom (' + args.custom_model + ')' if args.custom_model else f'MACE-MP-0 ({args.model})'}", f_out)
    if compare:
        print_dual(f"Comparison        : also evaluating with the foundation model ({args.model})", f_out)
    if report_path:
        print_dual(f"Report file       : {report_path}", f_out)

    sections("STRUCTURE")
    if not os.path.isfile(args.file):
        fail(f"Input file '{args.file}' not found.")
    if args.custom_model and not os.path.isfile(args.custom_model):
        fail(f"--custom-model file not found: {args.custom_model}")
    try:
        structure = structure_io.read_fdf(args.file)
    except Exception as e:
        fail(f"Could not read '{args.file}': {e}")
    atoms = AseAtomsAdaptor.get_atoms(structure_io.to_pymatgen(structure))
    vacuum_axes = kspace.detect_vacuum_axes(atoms.get_scaled_positions(), np.array(atoms.get_cell()), args.vacuum_gap)
    nvac = sum(vacuum_axes)
    print_dual(f"[OK] Read {len(atoms)} atom(s), dimensionality: {kspace.dimensionality_label(vacuum_axes)}", f_out)
    if args.dimensionality != "auto":
        dimensionality = args.dimensionality
        print_dual(f"Dimensionality    : {dimensionality.upper()} (manual override)", f_out)
    elif nvac == 3:
        fail("Detected a 0D (isolated molecule) structure -- elastic constants aren't defined for an isolated molecule.")
    elif nvac == 0:
        dimensionality = "3d"
    elif nvac == 1 and list(vacuum_axes) == [False, False, True]:
        dimensionality = "2d"
    elif nvac == 2 and list(vacuum_axes) == [True, True, False]:
        dimensionality = "1d"
    else:
        warn("Vacuum padding on an axis this tool does not support (2D expects the vacuum along z, 1D a wire along c): "
             "treated as 3D (GPa). Pass --dimensionality to override.")
        dimensionality = "3d"
    unit = UNIT_LABELS[dimensionality]
    periodic_axes = {i for i in range(3) if not vacuum_axes[i]}
    cp = atoms.cell.cellpar()
    print_dual(f"Formula           : {atoms.get_chemical_formula(mode='hill')}", f_out)
    print_dual(f"Cell              : a {cp[0]:.5f}, b {cp[1]:.5f}, c {cp[2]:.5f} Ang; alpha {cp[3]:.2f}, beta "
               f"{cp[4]:.2f}, gamma {cp[5]:.2f} deg", f_out)
    vol = atoms.get_volume()
    print_dual(f"Volume            : {vol:.3f} Ang^3 ({vol / len(atoms):.4f} Ang^3/atom)"
               + (f", density {atoms.get_masses().sum() / vol * AMU_A3_TO_G_CM3:.4f} g/cm^3" if dimensionality == "3d"
                  else " (includes vacuum)"), f_out)
    print_dual(f"Shortest distance : {fmt(model_card.shortest_distance(atoms), '{:.4f}')} Ang", f_out)
    pg0 = point_group(atoms)
    if pg0:
        print_dual(f"Symmetry          : {pg0[0]}, point group {pg0[1]} ({pg0[2]})", f_out)
    print_dual(f"Constants         : {dimensionality.upper()}, in {unit}"
               + {"2d": " (stress x cell height: per-area stiffness)", "1d": " (axial force of the wire)"}.get(
                   dimensionality, ""), f_out)

    sections("MODEL")
    model_arg = args.custom_model if args.custom_model else args.model
    t0 = time.monotonic()
    with capture_library_noise(noise, "model"):
        calc = mace_relax.get_calculator(model=model_arg, device=args.device)
    if card is not None:
        model_card.print_model_card(card, args.custom_model_card, args.custom_model, f_out, warn)
        for msg in model_card.check_structure(card, atoms):
            warn(f"Model card check: {msg}")
    else:
        print_dual(f"Model             : MACE-MP-0 {args.model} (foundation, no card): trained on the Materials Project "
                   "PBE/PBE+U relaxations -- not on your SIESTA level of theory", f_out)
    print_dual(f"Loaded in         : {time.monotonic() - t0:.1f} s; float64 (stresses are compared at the 0.01 GPa level)",
               f_out)
    label_main = "fine-tuned" if args.custom_model else f"MACE-MP-0 {args.model}"
    clamped = args.clamped and args.relax and len(atoms) > 1

    def finish(status_lines):
        sections("SUMMARY & FILES")
        for line in status_lines:
            print_dual(line, f_out)
        for kind, path in files:
            print_dual(f"{kind:<18}: {path}", f_out)
        if report_path:
            print_dual(f"Report            : {report_path}", f_out)
        not_asked = [f for f, on in (("--save-gnuplot", args.save_gnuplot), ("--save-data", args.save_data),
                                     ("--show", args.show), ("--save-report", args.save_report)) if not on]
        if not_asked:
            print_dual(f"Also available    : {' '.join(not_asked)}", f_out)
        print_dual(f"Wall time         : {time.monotonic() - t_start:.1f} s", f_out)
        if args.show and not display_available():
            warn("--show needs a display (none found, or a non-interactive matplotlib backend): the figures were "
                 "written as PNG files only.")
        sections("WARNINGS")
        if warns:
            for i, w in enumerate(warns, 1):
                print_dual(color_text(f"  {i}. {w}", 'yellow'), f_out)
        else:
            print_dual(color_text("  none", 'green'), f_out)
        if noise:
            print_dual(f"Library messages  : {len(noise)} distinct message(s) from mace/torch/ase captured, not shown",
                       f_out)
        if f_out:
            f_out.close()
        print("\n" + "-" * 60)
        print(color_text("ML elastic constants complete -- a fast preview, not a substitute for DFT elastic "
                         "constants.\n", 'bold'))
        import matplotlib.pyplot as plt
        if args.show and display_available() and figures:
            plt.show()
        for fig in figures:
            plt.close(fig)

    def gnuplot(series, three_d):
        plot_dir = os.path.join(args.output_dir, PLOT_DIR)
        scripts = gnuplot_files(plot_dir, series, three_d)
        pdfs, msg = run_gnuplot(plot_dir, scripts)
        files.append(("gnuplot folder", f"{plot_dir}/ ({', '.join(scripts)})"))
        for p in pdfs:
            files.append(("gnuplot figure", p))
        if msg:
            warn(msg)

    # ---- convergence mode
    if args.check_convergence:
        strains = sorted(args.convergence_strains or [args.max / 2.0, args.max, args.max * 2.0])
        sections("STRAIN-AMPLITUDE CONVERGENCE CHECK")
        rows = []
        for mx in strains:
            sub = copy.copy(args)
            sub.max = mx
            print(f"  ... max strain {mx:g} %", flush=True)
            res = run_model(atoms, calc, sub, vacuum_axes, periodic_axes, dimensionality, f"{mx:g} %", noise)
            rows.append({"max": mx, "res": res})
        ref = rows[0]["res"]
        groups = independent_constants(ref.get("C_proj", ref["C"]), 1e-3 * np.abs(ref["C"]).max())
        names = [g[0] for g in groups]
        hdr = ["max strain %"] + [" = ".join(n[:2]) + (" ..." if len(n) > 2 else "") for n in names] + ["min R^2"]
        if dimensionality == "3d":
            hdr += ["B Hill", "G Hill"]
        tab = []
        for r in rows:
            C = r["res"]["C"]
            cells = [f"{r['max']:g}"] + [f"{C[int(n[0][1]) - 1, int(n[0][2]) - 1]:.2f}" for n in names]
            cells.append(f"{min(f['r_squared'] for f in r['res']['fit']):.5f}")
            if dimensionality == "3d":
                h = vrh(C).get("H")
                cells += [fmt(h and h["B"]), fmt(h and h["G"])]
            tab.append((cells, None))
        print_dual(f"Strain amplitudes : {', '.join(f'{x:g}' for x in strains)} %, {args.steps} points per direction, "
                   f"{unit}", f_out)
        print_table(hdr, tab, f_out)
        spread = []
        for n in names:
            m_, n_ = int(n[0][1]) - 1, int(n[0][2]) - 1
            v = np.array([r["res"]["C"][m_, n_] for r in rows])
            spread.append((n[0], 100 * (v.max() - v.min()) / abs(v).mean() if abs(v).mean() > 1e-9 else 0.0))
        worst = max(spread, key=lambda x: x[1])
        print_dual(f"\nLargest change    : {worst[0]} varies by {worst[1]:.2f} % over the scan", f_out)
        print_dual("Reading           : a constant that drifts with the amplitude is non-linear there; one that scatters "
                   "at the SMALLEST amplitude is limited by the precision of the stress.", f_out)
        if worst[1] > 2.0:
            warn(f"Convergence: {worst[0]} changes by {worst[1]:.1f} % between {strains[0]:g} and {strains[-1]:g} % strain; "
                 "use an amplitude where the constants are flat.")
        p = os.path.join(args.output_dir, "strain_convergence.png")
        figures.append(plot_convergence(rows, names, unit, p))
        files.append(("Convergence plot", p))
        if args.save_data:
            p = os.path.join(args.output_dir, "strain_convergence.dat")
            with open(p, "w") as f:
                f.write("# max_strain(%) " + " ".join(f"{n[0]}({unit})" for n in names) + " min_R2\n")
                for r in rows:
                    C = r["res"]["C"]
                    f.write(f"{r['max']:.4f} " + " ".join(f"{C[int(n[0][1]) - 1, int(n[0][2]) - 1]:.5f}" for n in names)
                            + f" {min(f_['r_squared'] for f_ in r['res']['fit']):.6f}\n")
            files.append(("Data", p))
        finish(["Status            : OK" + (f" ({len(warns)} warning(s), see the last section)" if warns else ""),
                f"Results           : largest change {worst[0]} {worst[1]:.2f} % over {strains[0]:g}-{strains[-1]:g} %"])
        return

    # ---- single fit
    print(f"  ... {label_main}: pre-relax and strain scan", flush=True)
    main_res = run_model(atoms, calc, args, vacuum_axes, periodic_axes, dimensionality, label_main, noise,
                         clamped=clamped)
    if not main_res["dirs"]:
        fail("No direction is physically valid for this structure (all touch a vacuum-padded axis).")
    series = [main_res]
    if compare:
        print(f"  ... foundation model ({args.model}): pre-relax and strain scan", flush=True)
        with capture_library_noise(noise, "foundation model"):
            fcalc = mace_relax.get_calculator(model=args.model, device=args.device)
        series.append(run_model(atoms, fcalc, args, vacuum_axes, periodic_axes, dimensionality,
                                f"MACE-MP-0 {args.model}", noise))

    sections("REFERENCE STATE")
    print_dual("Each model's cell is relaxed (cell + atoms) to its own stress-free minimum before straining.", f_out)
    for res in series:
        report_reference(res, f_out, warn, dimensionality)

    sections("STRESS-STRAIN FITS")
    print_dual("Per direction: the fitted diagonal constant, R^2 and intercept (stress at zero strain, which should be "
               "~0) of the linear fit, the largest stress reached, and how far the atoms relaxed inside the cell "
               "(largest shift from the homogeneous strain; ~0 for atoms at inversion centres).", f_out)
    for res in series:
        report_directions(res, f_out, warn, args)
    p = os.path.join(args.output_dir, "stress_strain.png")
    fig = plot_stress_strain(series, unit, p)
    if fig is not None:
        figures.append(fig)
        files.append(("Stress-strain plot", p))

    sections("TRAINED RANGE")
    report_domain(series, card, f_out, warn, dimensionality)

    sections("STIFFNESS TENSOR & CRYSTAL SYMMETRY")
    report_tensor(main_res, f_out, warn, dimensionality)
    if main_res["missing"]:
        print_dual(f"\nNo data for       : {', '.join(main_res['missing'])} (vacuum axes) -- those columns are zero",
                   f_out)

    sections("MECHANICAL STABILITY")
    system = main_res["pg"][2] if main_res.get("pg") else None
    for res in series:
        rows, ok = stability(res["C"], dimensionality, system)
        print_dual(f"--- {res['label']} (Born criteria; Mouhat & Coudert 2014)", f_out)
        print_table(["criterion", f"value ({unit})", "verdict"],
                    [([c, v, "PASS" if p_ else "FAIL"], None if p_ else "red") for c, v, p_ in rows], f_out)
        print_dual(f"VERDICT           : {color_text('STABLE', 'green') if ok else color_text('UNSTABLE', 'red')}", f_out)
        res["stable"] = ok
        if not ok:
            warn(f"{res['label']}: mechanically UNSTABLE by the Born criteria -- the structure is not a minimum of this "
                 "model under some strain, or the fit is unreliable.")

    summary = {}
    if dimensionality == "3d" and not main_res["missing"]:
        sections("ELASTIC MODULI")
        print_dual("Polycrystalline averages: Voigt assumes uniform strain, Reuss uniform stress; Hill is their mean "
                   "(Hill 1952).", f_out)
        summary = report_moduli(main_res, f_out, warn, system)
        if "C_clamped" in main_res:
            sections("INTERNAL RELAXATION (clamped vs relaxed ions)")
            report_clamped(main_res, f_out, warn)
        sections("SOUND VELOCITIES & DEBYE TEMPERATURE")
        summary["sound"] = report_sound(main_res, f_out, warn)
        sections("DIRECTIONAL ANISOTROPY")
        summary["dir"] = report_directional(main_res, f_out)
        p = os.path.join(args.output_dir, "directional_moduli.png")
        figures.append(plot_directional(series, p))
        files.append(("Directional plot", p))
        if not args.skip_anisotropy_surface:
            p = os.path.join(args.output_dir, "anisotropy_surface.png")
            fig = plot_surface(main_res["C"], p, main_res["label"])
            if fig is not None:
                figures.append(fig)
                files.append(("Anisotropy surface", p))
    elif dimensionality == "2d":
        sections("IN-PLANE MODULI")
        report_2d(main_res, f_out)
        p = os.path.join(args.output_dir, "inplane_youngs.png")
        figures.append(plot_inplane(series, unit, p))
        files.append(("In-plane E(θ) plot", p))
    elif dimensionality == "3d":
        warn(f"No data for {', '.join(main_res['missing'])}: the 3D moduli, sound velocities and directional properties "
             "need all six directions and are skipped.")

    if len(series) > 1:
        sections("FINE-TUNED vs FOUNDATION")
        report_comparison(series, f_out, system)

    if args.save_data:
        p = os.path.join(args.output_dir, "stress_strain.dat")
        write_stress_strain_dat(p, main_res)
        files.append(("Data", p))
        p = os.path.join(args.output_dir, "elastic_tensor.dat")
        write_tensor_dat(p, main_res)
        files.append(("Data", p))
        if dimensionality == "3d" and not main_res["missing"]:
            p = os.path.join(args.output_dir, "directional.dat")
            write_directional_dat(p, series)
            files.append(("Data", p))
    if args.save_gnuplot:
        gnuplot(series, dimensionality == "3d" and not main_res["missing"])

    lines = ["Status            : OK" + (f" ({len(warns)} warning(s), see the last section)" if warns else "")]
    h = (summary.get("vrh") or {}).get("H")
    if h:
        th = (summary.get("sound") or {}).get("theta")
        lines.append(f"Results           : B {h['B']:.1f}, G {h['G']:.1f}, E {h['E']:.1f} GPa, ν {h['nu']:.3f}"
                     + (f", θ_D {th:.0f} K" if th else "") + f", {'stable' if main_res['stable'] else 'UNSTABLE'}")
    finish(lines)


if __name__ == "__main__":
    main()
