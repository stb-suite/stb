#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

from stb import __version__ as VERSION

import os
import glob
import re
import sys
import shutil
import argparse
from datetime import datetime

import numpy as np
from pymatgen.io.ase import AseAtomsAdaptor

from stb.core import structure_io
from stb.core.adsorption_sites import SINGLE_POINT_BLOCK, CONFIG_EXTRA_FILE
from stb.core.pseudopotentials import BANKS, get_required_pseudos, resolve_pseudo_source
from stb.core.cli import color_text, show_intro, print_dual, print_section, print_table, capture_library_noise
from stb.core.fdf_directives import parse_directives, parse_kgrid
from stb.core.kspace import detect_vacuum_axes, dimensionality_label
from stb.core import mlff_workflow
from stb.core.atomic_env import EnvironmentSet

REPORT_FILE = "stb_mlff_report.txt"

# Every training folder is ONE fixed sample point: SIESTA must evaluate it, not
# move it. A relaxation template (MD.TypeOfRun CG + MD.Steps > 0) would drift the
# geometry towards equilibrium, and MD.VariableCell T would relax a strained cell
# back to the reference volume -- either way the labels end up describing a
# different structure than the one sampled. Forced through config_extra.fdf
# (first-occurrence-wins, see structure_io.prepend_include), same mechanism as
# the other workflows' single-point folders.
_MLFF_EXTRA_BLOCK = (
    SINGLE_POINT_BLOCK
    + "MD.NumCGsteps         0\n"
    + "MD.VariableCell       false\n"
)


def build_rattled_configs(atoms, coord_format, n_configs, stdevs, seed):
    """Returns a list of `n_configs` position arrays (same coordinate system
    as `coord_format`, matching the reference calc.fdf's own
    AtomicCoordinatesFormat), each a Gaussian-displaced ("rattled") copy of
    `atoms`'s equilibrium positions. Cycles through `stdevs` (one or more
    displacement amplitudes) round-robin across the `n_configs` -- using
    several amplitudes instead of just one is standard MLIP training-set
    practice, so the fitted potential sees more than one noise scale.
    """
    configs = []
    for i in range(n_configs):
        stdev = stdevs[i % len(stdevs)]
        a = atoms.copy()
        a.rattle(stdev=stdev, seed=seed + i)
        if coord_format == "fractional":
            configs.append(a.get_scaled_positions(wrap=False))
        else:
            configs.append(a.get_positions())
    return configs


def build_strained_configs(atoms, coord_format, raw_lattice, strains_pct, stdevs, seed):
    """Returns one (positions, raw_lattice) pair per entry of `strains_pct`
    (volumetric strain, %): the cell scaled isotropically by
    (1 + strain/100)**(1/3) -- same convention as stb-eosInputs/stb-mleos --
    and the atoms carried with it (same fractional coordinates), then rattled
    with the `stdevs` amplitudes cycled round-robin.

    Rattled rather than left ideal: at ideal positions of a high-symmetry
    crystal every force vanishes by symmetry, so such a folder would teach
    the model energy and stress only. Rattled, it also adds forces at a
    volume the position-only configurations never reach. Without at least a
    few of these, the training set has a single volume and the fitted model
    has nothing to learn the volume response (E(V), bulk modulus, pressure)
    from. `raw_lattice` is the reference cell in the file's own
    LatticeConstant basis (structure_io.raw_lattice_vectors).

    Only the periodic axes are scaled -- a slab's vacuum or a molecule's box
    keeps its size, as in stb-mlffActiveLearning --sampling-method strain --,
    each by (1 + strain/100)**(1/n_periodic), so the strain is still the
    requested change of the cell volume (an area strain for a 2D cell). With
    no periodic axis at all the cell is scaled isotropically.
    """
    from stb.core.kspace import detect_vacuum_axes
    vac = detect_vacuum_axes(atoms.get_scaled_positions(), np.asarray(atoms.cell), 5.0)
    periodic = [not v for v in vac]
    if not any(periodic):
        periodic = [True, True, True]
    n_per = sum(periodic)
    configs = []
    scaled0 = atoms.get_scaled_positions(wrap=False)
    for i, strain in enumerate(strains_pct):
        f_axis = (1.0 + strain / 100.0) ** (1.0 / n_per)
        factor = np.array([f_axis if p else 1.0 for p in periodic])[:, None]   # per lattice vector (row)
        a = atoms.copy()
        a.set_cell(np.asarray(atoms.get_cell()) * factor, scale_atoms=False)
        a.set_scaled_positions(scaled0)
        a.rattle(stdev=stdevs[i % len(stdevs)], seed=seed + 10000 + i)
        if coord_format == "fractional":
            positions = a.get_scaled_positions(wrap=False)
        else:
            positions = a.get_positions()
        configs.append((positions, np.asarray(raw_lattice) * factor))
    return configs


def sample_aimd_configs(label, coord_format, lattice_constant, n_samples, f_out=None):
    """Returns a list of (positions, raw_lattice_or_None) pairs, one per
    sampled AIMD frame -- `n_samples` frames spaced evenly across
    `<label>.ANI`. raw_lattice is None when the frame's cell matches the
    reference structure closely enough not to need rewriting (same
    LatticeConstant-normalized convention rewrite_fdf_lattice expects,
    i.e. already divided by `lattice_constant`), otherwise the frame's own
    cell (for a variable-cell/NPT run, so each sampled snapshot keeps its
    real cell instead of silently reusing the reference one).
    """
    from stb.core.deps import require_sisl
    sisl = require_sisl()
    from stb.core.md_traj import read_frame_lattices

    ani_file = f"{label}.ANI"
    if not os.path.exists(ani_file):
        raise FileNotFoundError(f"'{ani_file}' not found.")

    sile = sisl.get_sile(ani_file)
    all_frames = sile.read_geometry[:]()
    if not all_frames:
        raise ValueError(f"'{ani_file}' contains no frames.")

    all_cells, _ = read_frame_lattices(label, len(all_frames), f_out)
    if all_cells is None:
        raise ValueError(
            f"Could not find a readable '{label}.out', '{label}.XV' or '{label}.fdf' "
            f"-- a lattice is required ('{ani_file}' doesn't carry one of its own).")

    n_samples = min(n_samples, len(all_frames))
    indices = np.linspace(0, len(all_frames) - 1, n_samples).round().astype(int)
    indices = sorted(set(indices.tolist()))

    configs = []
    for idx in indices:
        cell = np.array(all_cells[idx])
        cart = all_frames[idx].xyz
        if coord_format == "fractional":
            pos = cart @ np.linalg.inv(cell)
        else:
            pos = cart
        raw_lattice = cell / lattice_constant
        configs.append((pos, raw_lattice))
    return configs


TEST_SEED_OFFSET = 777000
TEST_DIR = "mlff_test"


def build_test_configs(atoms, coord_format, raw_lattice, n_test, stdevs, strains_pct, seed):
    """`n_test` held-out configurations drawn like the training ones but with an
    independent seed: the kinds -- one rattle per --stdev, one strained cell per
    --volume-strains entry -- are cycled round-robin. Returns a list of
    (positions, raw_lattice or None, kind, origin). The point is a test set
    from the same distribution the model is meant for, never trained on, so
    stb-mlffValidate can tell how a round really does and whether the next one
    helped."""
    kinds = [("rattle", s) for s in stdevs] + [("strain", p) for p in strains_pct]
    out = []
    for i in range(n_test):
        kind, value = kinds[i % len(kinds)]
        stdev = stdevs[i % len(stdevs)]
        if kind == "rattle":
            pos = build_rattled_configs(atoms, coord_format, 1, [value], seed + TEST_SEED_OFFSET + i)[0]
            out.append((pos, None, "rattle", f"rattle {value:g} Ang"))
        else:
            pos, lat = build_strained_configs(atoms, coord_format, raw_lattice, [value], [stdev],
                                              seed + TEST_SEED_OFFSET + i)[0]
            out.append((pos, lat, "strain", f"volume {value:+g} % (rattle {stdev:g} Ang)"))
    return out


def write_config_folder(output_root, index, ref_fdf, ref_basename, positions,
                         raw_lattice, pseudos_to_copy, extra_block="", prefix="mlff_config",
                         calc_template=None):
    """Writes mlff_config_NNN/<ref_basename> (positions rewritten via
    rewrite_fdf_positions, lattice ALSO rewritten via rewrite_fdf_lattice
    when raw_lattice is given -- e.g. a variable-cell AIMD snapshot) plus
    the copied pseudopotentials, ready for the user to run SIESTA in, plus a
    config_extra.fdf forcing a single-point evaluation (%include'd at the top
    of the .fdf; see _MLFF_EXTRA_BLOCK). `extra_block` (e.g. a supercell's own
    k-grid, stb-mlffActiveLearning --supercell) is written first in that file,
    so it wins over the template. `prefix` names the folder (<prefix>_NNN;
    the held-out test set uses "test"). With `calc_template` (the suite's usual
    structure.fdf + calc.fdf pair), `ref_fdf` is the structure file -- rewritten
    as above -- and the calc.fdf is copied next to it with the config_extra.fdf
    %include on top, its own '%include <structure>' left as it is. Returns the
    folder path.
    """
    folder = os.path.join(output_root, f"{prefix}_{index:03d}")
    os.makedirs(folder, exist_ok=True)
    out_fdf = os.path.join(folder, ref_basename)

    structure_io.rewrite_fdf_positions(ref_fdf, positions, out_fdf)
    if raw_lattice is not None:
        structure_io.rewrite_fdf_lattice(out_fdf, raw_lattice, out_fdf)

    with open(os.path.join(folder, CONFIG_EXTRA_FILE), "w") as f:
        f.write(extra_block + _MLFF_EXTRA_BLOCK)
    main_fdf = out_fdf
    if calc_template is not None:
        main_fdf = os.path.join(folder, os.path.basename(calc_template))
        shutil.copy(calc_template, main_fdf)
    with open(main_fdf) as f:
        fdf_text = f.read()
    if f"%include {CONFIG_EXTRA_FILE}" not in fdf_text:
        with open(main_fdf, "w") as f:
            f.write(structure_io.prepend_include(fdf_text, CONFIG_EXTRA_FILE))

    for pseudo_path in pseudos_to_copy:
        shutil.copy(pseudo_path, os.path.join(folder, os.path.basename(pseudo_path)))

    return folder


# --------------------------------------------------------------------------
# Report helpers: what the reference is, how the template will run, and how
# far each written configuration is from the reference.

VACUUM_GAP = 5.0          # Ang of empty space that marks a non-periodic axis
NEIGHBOUR_CUTOFF = 6.0    # Ang, for nearest-neighbour distances
CLOSE_WARN, CLOSE_FAIL = 0.80, 0.65   # d_min / reference d_min that flags a config

# (shown label, normalized fdf label(s), what SIESTA does when it is absent)
_SETTINGS = [
    ("SystemLabel", ("systemlabel",), "siesta (SIESTA default)"),
    ("XC.functional", ("xcfunctional",), "SIESTA default (LDA)"),
    ("XC.authors", ("xcauthors",), "SIESTA default (PZ)"),
    ("PAO.BasisSize", ("paobasissize",), "SIESTA default (DZP)"),
    ("PAO.EnergyShift", ("paoenergyshift",), "SIESTA default"),
    ("Mesh.Cutoff", ("meshcutoff",), "SIESTA default"),
    ("ElectronicTemperature", ("electronictemperature",), "SIESTA default"),
    ("Spin", ("spin", "spinpolarized"), "not set (non-polarized)"),
    ("DM.Tolerance", ("scfdmtolerance", "dmtolerance"), "SIESTA default"),
    ("MaxSCFIterations", ("scfmaxiterations", "maxscfiterations"), "SIESTA default"),
]
# Template directives that config_extra.fdf overrides in every folder.
_OVERRIDDEN = [("MD.TypeOfRun", "mdtypeofrun"), ("MD.Steps", "mdsteps"),
               ("MD.NumCGsteps", "mdnumcgsteps"), ("MD.VariableCell", "mdvariablecell")]


def fdf_atoms(path):
    """ASE Atoms of the structure in an .fdf (as SIESTA will read it)."""
    return AseAtomsAdaptor.get_atoms(structure_io.to_pymatgen(structure_io.read_fdf(path)))


def pair_min_distances(atoms, cutoff=NEIGHBOUR_CUTOFF):
    """{(A, B): shortest A-B distance} over periodic images, within `cutoff`."""
    from ase.neighborlist import neighbor_list
    i, j, d = neighbor_list("ijd", atoms, cutoff)
    sym = atoms.get_chemical_symbols()
    out = {}
    for a, b, x in zip(i, j, d):
        key = tuple(sorted((sym[a], sym[b])))
        out[key] = min(out.get(key, np.inf), float(x))
    return out


def min_distance(atoms, cutoff=NEIGHBOUR_CUTOFF):
    pairs = pair_min_distances(atoms, cutoff)
    return min(pairs.values()) if pairs else None


def displacements(ref, atoms):
    """Per-atom displacement (Ang) from the reference: the fractional
    difference taken to the nearest image and expressed in the configuration's
    own cell, so a strained cell's uniform scaling is not counted as motion."""
    df = atoms.get_scaled_positions(wrap=False) - ref.get_scaled_positions(wrap=False)
    df -= np.round(df)
    return np.linalg.norm(df @ np.asarray(atoms.cell), axis=1)


def config_stats(ref, atoms, ref_dmin):
    u = displacements(ref, atoms)
    dmin = min_distance(atoms)
    return {"n_atoms": len(atoms), "volume": float(atoms.get_volume()),
            "dV_pct": 100.0 * (atoms.get_volume() / ref.get_volume() - 1.0),
            "u_rms": float(np.sqrt(np.mean(u ** 2))), "u_max": float(u.max()),
            "d_min": dmin, "d_ratio": (dmin / ref_dmin) if (dmin and ref_dmin) else None}


def calc_settings(texts):
    """[(label, value, is_default)] for the settings worth showing, read
    first-occurrence-wins over `texts` in SIESTA's reading order."""
    d = parse_directives(*texts)
    rows = []
    for label, keys, default in _SETTINGS:
        value = next((d[k] for k in keys if k in d), None)
        rows.append((label, value if value is not None else default, value is None))
    return rows, d


def report_reference(atoms, ref_file, calc_template, structure, f_out):
    """Prints the reference: composition, cell, density, dimensionality,
    nearest-neighbour distances. Returns (vacuum_axes, pair distances)."""
    print_dual(f"Structure source  : {ref_file}" + ("" if calc_template else " (self-contained .fdf)"), f_out)
    if calc_template:
        print_dual(f"Calc template     : {calc_template} (%include {os.path.basename(ref_file)})", f_out)
    counts = {}
    for sym in atoms.get_chemical_symbols():
        counts[sym] = counts.get(sym, 0) + 1
    print_dual(f"Formula           : {atoms.get_chemical_formula(mode='hill')} "
               f"({len(atoms)} atoms: {', '.join(f'{k} {v}' for k, v in counts.items())})", f_out)
    a, b, c, al, be, ga = atoms.cell.cellpar()
    print_dual(f"Cell a, b, c      : {a:.4f}  {b:.4f}  {c:.4f} Ang", f_out)
    print_dual(f"Cell angles       : {al:.2f}  {be:.2f}  {ga:.2f} deg", f_out)
    vol = atoms.get_volume()
    mass = float(np.sum(atoms.get_masses()))
    print_dual(f"Volume            : {vol:.3f} Ang^3 ({vol / len(atoms):.3f} Ang^3/atom), "
               f"density {mass / vol * 1.66053907:.4f} g/cm^3", f_out)
    print_dual(f"Coordinates       : {structure.coord_format}, LatticeConstant "
               f"{structure.lattice_constant:g}", f_out)
    vac = detect_vacuum_axes(atoms.get_scaled_positions(), np.asarray(atoms.cell), VACUUM_GAP)
    dim = dimensionality_label(vac)
    print_dual(f"Dimensionality    : {dim}" +
               (f" -- vacuum along {', '.join('abc'[i] for i, v in enumerate(vac) if v)}" if any(vac) else ""), f_out)
    pairs = pair_min_distances(atoms)
    if pairs:
        txt = ", ".join(f"{p[0]}-{p[1]} {d:.3f}" for p, d in sorted(pairs.items(), key=lambda x: x[1]))
        print_dual(f"Shortest distances: {txt} Ang", f_out)
    else:
        print_dual(f"Shortest distances: none within {NEIGHBOUR_CUTOFF:g} Ang", f_out)
    return vac, pairs


def report_settings(texts, atoms, vac, pseudo_label, pseudos, f_out):
    """Prints the calculation settings every folder inherits, the k-grid and
    what config_extra.fdf overrides. Returns a list of warning strings."""
    warnings = []
    rows, d = calc_settings(texts)
    width = max(18, max(len(label) for label, _, _ in rows))
    for label, value, is_default in rows:
        print_dual(f"{label:<{width}}: {value}" + ("  (not in template)" if is_default else ""), f_out)
    kgrid = parse_kgrid(*texts)
    lengths = atoms.cell.lengths()
    if kgrid:
        ka = "  ".join(f"{k * l:.1f}" for k, l in zip(kgrid, lengths))
        print_dual(f"{'k-grid':<{width}}: {kgrid[0]} x {kgrid[1]} x {kgrid[2]}  (k x a = {ka} Ang)", f_out)
        thin = [("abc"[i]) for i in range(3) if not vac[i] and kgrid[i] * lengths[i] < 15.0]
        if thin:
            warnings.append(f"k-grid is coarse along {', '.join(thin)} (k x a < 15 Ang): forces/stress "
                            "may not be converged -- the model learns whatever SIESTA gives it.")
    elif "kgridcutoff" in d:
        print_dual(f"{'k-grid':<{width}}: kgrid.Cutoff {d['kgridcutoff']}", f_out)
    else:
        print_dual(f"{'k-grid':<{width}}: none (Gamma only)", f_out)
        if not all(vac):
            warnings.append("No k-grid in the template: SIESTA samples Gamma only, usually far "
                            "from converged for a periodic cell this small.")
    forced = [f"{lab} {d[k]}" for lab, k in _OVERRIDDEN if k in d]
    print_dual(f"{'Run type':<{width}}: single point in every folder (config_extra.fdf)" +
               (f" -- template's {', '.join(forced)} overridden" if forced else ""), f_out)
    print_dual(f"{'Pseudopotentials':<{width}}: {pseudo_label} -> "
               + ", ".join(f"{os.path.basename(p)} ({os.path.getsize(p) / 1024:.0f} KiB)" for p in pseudos), f_out)
    return warnings


DUP_THRESHOLD = 0.25      # novelty (typical-spacing units) below which a config is a near-duplicate
ML_FORCE_WARN = mlff_workflow.UNPHYSICAL_FMAX   # eV/Ang: foundation-model max force of an unphysical geometry
ML_ENERGY_WARN = mlff_workflow.UNPHYSICAL_DE    # meV/atom above the reference, same purpose
TRAJ_FORMATS = {"xsf": "xsf", "pdb": "proteindatabank", "xyz": "extxyz"}


def config_novelty(stats, frames):
    """Fills stats[i]['novelty']: how far config i's atomic environments are
    from those of the training configs written before it (a test config:
    from every training config and the earlier test ones), in units of the
    whole set's typical spacing (core.atomic_env, the same fingerprint as
    stb-mlffActiveLearning's novelty filter). A value near 0 = a near-copy
    of another config: one more SIESTA run that teaches the model nothing,
    or a test config that does not test anything new. Returns the unit."""
    everything = EnvironmentSet()
    for a in frames:
        everything.add(a)
    scale = everything.typical_spacing()
    seen = EnvironmentSet()
    for st, a in zip(stats, frames):
        st["novelty"] = seen.novelty(a, scale) if (scale and seen.n_configs) else None
        if st["novelty"] is not None and not np.isfinite(st["novelty"]):
            st["novelty"] = None
        seen.add(a)
    return scale


def ml_screen(frames, model, device):
    """Foundation-model single points on every frame (frames[0] = reference):
    list of {'e_atom' eV/atom, 'f_max' eV/Ang, 'f_rms' eV/Ang, 'p_gpa' or None}."""
    from stb.core.deps import require_mace
    require_mace()
    from stb.core.mace_relax import get_calculator
    calc = get_calculator(model, device, allow_uncarded=True)
    out = []
    for a in frames:
        b = a.copy()
        b.calc = calc
        e = b.get_potential_energy()
        f = np.linalg.norm(b.get_forces(), axis=1)
        try:
            p = -float(np.mean(b.get_stress(voigt=True)[:3])) / 0.006241509074   # eV/Ang^3 -> GPa
        except Exception:
            p = None
        out.append({"e_atom": float(e / len(b)), "f_max": float(f.max()),
                    "f_rms": float(np.sqrt(np.mean(f ** 2))), "p_gpa": p})
    return out


def plot_coverage_png(out_path, ref, frames, stats, ref_dmin, with_ml, show):
    """The coverage plots as one matplotlib figure (the gnuplot files stay the
    workflow's plot convention; this is what --show displays, like Stages 2-4):
    pair-distance distribution, displacement vs volume, shortest distance vs
    displacement, and the foundation model's max force (--ml-screen) or the
    novelty of each config."""
    import matplotlib.pyplot as plt
    from ase.neighborlist import neighbor_list
    bins = np.arange(0.0, NEIGHBOUR_CUTOFF + 1e-9, 0.05)
    centers = 0.5 * (bins[1:] + bins[:-1])

    def density(atoms_list):
        hist, n = np.zeros(len(centers)), 0
        for a in atoms_list:
            hist += np.histogram(neighbor_list("d", a, NEIGHBOUR_CUTOFF), bins=bins)[0]
            n += len(a)
        return hist / max(n * 0.05, 1e-12)

    sets = {w: [i for i, st in enumerate(stats) if st["set"] == w] for w in ("train", "test")}
    style = {"train": dict(color="tab:blue", marker="o"), "test": dict(color="tab:orange", marker="s")}
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    ax = axes[0, 0]
    ax.plot(centers, density([ref]), "k-", lw=1, label="reference")
    for w, ii in sets.items():
        if ii:
            ax.plot(centers, density([frames[i] for i in ii]), "-", color=style[w]["color"], label=w)
    ax.set_xlabel("Pair distance (Ang)"); ax.set_ylabel("Pairs per atom per Ang")
    ax.set_title("Pair-distance distribution"); ax.legend(fontsize=8)
    for ax, (xk, yk, xl, yl, title) in zip(
            [axes[0, 1], axes[1, 0]],
            [("dV_pct", "u_rms", "Volume change (%)", "RMS displacement (Ang)", "Displacement vs volume"),
             ("u_rms", "d_min", "RMS displacement (Ang)", "Shortest distance (Ang)", "Shortest distance vs displacement")]):
        for w, ii in sets.items():
            if ii:
                ax.scatter([stats[i][xk] for i in ii], [stats[i][yk] if stats[i][yk] is not None else np.nan for i in ii],
                           s=30, label=w, **style[w])
        ax.set_xlabel(xl); ax.set_ylabel(yl); ax.set_title(title); ax.legend(fontsize=8)
    if ref_dmin:
        for frac, ls in ((CLOSE_WARN, "--"), (CLOSE_FAIL, ":")):
            axes[1, 0].axhline(frac * ref_dmin, color="tab:red", ls=ls, lw=0.8)
    ax = axes[1, 1]
    if with_ml:
        for w, ii in sets.items():
            if ii:
                ax.scatter([stats[i]["d_min"] or np.nan for i in ii], [stats[i]["ml_f_max"] for i in ii],
                           s=30, label=w, **style[w])
        ax.axhline(ML_FORCE_WARN, color="tab:red", ls="--", lw=0.8, label="unphysical")
        ax.set_xlabel("Shortest distance (Ang)"); ax.set_ylabel("Foundation-model max |F| (eV/Ang)")
        ax.set_title("Foundation-model screen")
    else:
        for w, ii in sets.items():
            if ii:
                ax.scatter(ii, [stats[i].get("novelty") if stats[i].get("novelty") is not None else np.nan for i in ii],
                           s=30, label=w, **style[w])
        ax.set_xlabel("Config #"); ax.set_ylabel("Novelty (typical spacings)"); ax.set_title("Novelty per config")
    ax.legend(fontsize=8)
    fig.suptitle("stb-mlff -- training/test set coverage", fontsize=13)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    if show:
        plt.show()
    plt.close(fig)


def write_coverage_plots(out_dir, ref, frames, stats, ref_dmin, with_ml):
    """gnuplot .dat + .gplot pairs (the workflow tools' plot convention) in
    `out_dir`: the pair-distance distribution of the reference, the training
    and the test set, and per-config scatter plots (displacement vs volume,
    shortest distance vs displacement, and with --ml-screen the foundation
    model's max force vs shortest distance). Returns the .gplot paths."""
    os.makedirs(out_dir, exist_ok=True)
    from ase.neighborlist import neighbor_list
    bins = np.arange(0.0, NEIGHBOUR_CUTOFF + 1e-9, 0.05)
    centers = 0.5 * (bins[1:] + bins[:-1])

    def density(atoms_list):
        if not atoms_list:
            return np.full(len(centers), np.nan)
        hist = np.zeros(len(centers))
        n = 0
        for a in atoms_list:
            hist += np.histogram(neighbor_list("d", a, NEIGHBOUR_CUTOFF), bins=bins)[0]
            n += len(a)
        return hist / (n * 0.05)

    sets = {w: [a for a, st in zip(frames, stats) if st["set"] == w] for w in ("train", "test")}
    dist = np.column_stack([centers, density([ref]), density(sets["train"]), density(sets["test"])])
    d_dat = os.path.join(out_dir, "coverage_distances.dat")
    with open(d_dat, "w") as f:
        f.write("# Pair-distance distribution (neighbours per atom per Ang), stb-mlff\n")
        f.write("# 1:r(Ang) 2:reference 3:training 4:test\n")
        for row in dist:
            f.write("  ".join("NaN" if not np.isfinite(v) else f"{v:.6f}" for v in row) + "\n")
    d_gp = os.path.join(out_dir, "coverage_distances.gplot")
    has_test = bool(sets["test"])
    # The reference's neighbour shells as thin vertical lines (its histogram is
    # a set of spikes that would dwarf the sampled distributions).
    from ase.neighborlist import neighbor_list as _nl
    shells = sorted({round(float(x), 2) for x in _nl("d", ref, NEIGHBOUR_CUTOFF)})
    with open(d_gp, "w") as f:
        f.writelines([
            "# --- STB Plot Configuration ---\n", "# Generated by stb-mlff\n", "reset\n",
            'set terminal pdfcairo enhanced color font "Arial,14" size 7,5\n',
            'set output "coverage_distances.pdf"\n\n',
            'set title "Interatomic distances sampled (stb-mlff)"\n',
            'set xlabel "Distance (Ang)"\n', 'set ylabel "Neighbours per atom per Ang"\n',
            'set grid lc rgb "#dddddd"\n', "set key top left opaque\n",
            f"set xrange [0:{NEIGHBOUR_CUTOFF:g}]\n",
            "# reference neighbour shells\n",
        ] + [f'set arrow from {d:.2f}, graph 0 to {d:.2f}, graph 1 nohead dt 3 lw 1 lc rgb "#888888"\n'
             for d in shells] + [
            f'plot "coverage_distances.dat" using 1:3 with lines lw 2 lc rgb "#2255cc" title "training"'
            + (', "" using 1:4 with lines lw 2 dt 2 lc rgb "#cc5522" title "test"' if has_test else "")
            + ', NaN with lines dt 3 lw 1 lc rgb "#888888" title "reference shells"\n',
        ])

    c_dat = os.path.join(out_dir, "coverage_configs.dat")
    cols = ["dV_pct", "u_rms", "u_max", "d_min", "novelty"]
    with open(c_dat, "w") as f:
        f.write("# One row per configuration; index 0 = training, index 1 = test, stb-mlff\n")
        f.write("# 1:dV(%) 2:u_rms(Ang) 3:u_max(Ang) 4:d_min(Ang) 5:novelty 6:ML_Fmax(eV/Ang) "
                "7:ML_dE(meV/atom) 8:folder\n")
        for k, which in enumerate(("train", "test")):
            if k:
                f.write("\n\n")
            for st in stats:
                if st["set"] != which:
                    continue
                vals = [st.get(c) for c in cols] + [st.get("ml_f_max"), st.get("ml_dE")]
                f.write("  ".join("NaN" if v is None else f"{v:.6f}" for v in vals)
                        + f"  {os.path.basename(st['folder'])}\n")
    c_gp = os.path.join(out_dir, "coverage_configs.gplot")
    panels = 3 if with_ml else 2

    def series(using):
        t = f'"coverage_configs.dat" index 0 using {using} with points pt 7 ps 1.2 lc rgb "#2255cc" title "training"'
        if has_test:
            t += f', "" index 1 using {using} with points pt 5 ps 1.2 lc rgb "#cc5522" title "test"'
        return t
    lines = ["# --- STB Plot Configuration ---\n", "# Generated by stb-mlff\n", "reset\n",
             f'set terminal pdfcairo enhanced color font "Arial,12" size {5 * panels},4.8\n',
             'set output "coverage_configs.pdf"\n',
             f"set multiplot layout 1,{panels} title \"Sampled configurations (stb-mlff)\"\n",
             'set grid lc rgb "#dddddd"\nset key outside top center horizontal\nset offsets graph 0.05, graph 0.05, graph 0.05, graph 0.05\n\n',
             'set xlabel "Volume change (%)"\nset ylabel "RMS displacement from reference (Ang)"\n',
             f"plot {series('1:2')}\n\n",
             'set xlabel "Max displacement (Ang)"\nset ylabel "Shortest distance (Ang)"\n']
    if ref_dmin:
        lines.append(f'plot {series("3:4")}, {ref_dmin * CLOSE_WARN:.4f} with lines dt 2 lc rgb "#bb9900" '
                     f'title "{int(CLOSE_WARN * 100)}% ref. d_{{min}}", {ref_dmin * CLOSE_FAIL:.4f} '
                     f'with lines dt 2 lc rgb "#cc0000" title "{int(CLOSE_FAIL * 100)}%"\n\n')
    else:
        lines.append(f"plot {series('3:4')}\n\n")
    if with_ml:
        lines += ['set xlabel "Shortest distance (Ang)"\nset ylabel "Foundation-model max force (eV/Ang)"\n',
                  f"plot {series('4:6')}\n\n"]
    lines.append("unset multiplot\n")
    with open(c_gp, "w") as f:
        f.writelines(lines)
    return [d_gp, c_gp]


def includes_file(calc_path, name):
    """True when `calc_path` has an '%include <name>' line (case-insensitive
    keyword, comments ignored) -- the structure.fdf + calc.fdf pair only works
    if the calc file actually pulls the structure in."""
    pattern = re.compile(r"^\s*%include\s+(\S+)", re.I)
    with open(calc_path) as f:
        for line in f:
            m = pattern.match(line.split("#", 1)[0])
            if m and os.path.basename(m.group(1).strip("'\"")) == name:
                return True
    return False


def resolve_reference(structure, calc, single):
    """(structure_path, calc_path_or_None, error_or_None) for the two input
    styles: the suite's usual pair -s structure.fdf -c calc.fdf (calc.fdf
    %includes the structure file), or one self-contained .fdf (--file)."""
    if single:
        if not os.path.isfile(single):
            return single, None, f"Reference file '{single}' not found."
        return single, None, None
    for label, path in (("Structure", structure), ("Calc", calc)):
        if not os.path.isfile(path):
            return structure, calc, (f"{label} file '{path}' not found (pass -s/--structure and "
                                     "-c/--calc, or one self-contained .fdf with --file).")
    if os.path.basename(structure) == os.path.basename(calc):
        return structure, calc, ("--structure and --calc must have different file names (both are "
                                 "copied into every folder).")
    if not includes_file(calc, os.path.basename(structure)):
        return structure, calc, (f"'{calc}' has no '%include {os.path.basename(structure)}' -- add "
                                 "it, or pass one self-contained .fdf with --file.")
    return structure, calc, None


def main():
    parser = argparse.ArgumentParser(
        description="Stage 1 (Prep) of the Custom ML Force Field workflow: generates "
                     "training configurations (rattled displacements and/or AIMD-sampled "
                     "snapshots) from a reference SIESTA calculation, ready for the user to "
                     "run through SIESTA, plus an optional held-out test set (--n-test). "
                     "Creates mlff_workflow.json, through which the next stages -- Train "
                     "(stb-mlffAnalysis), Validate (stb-mlffValidate), Active learning "
                     "(stb-mlffActiveLearning) -- find what the earlier ones did.",
        epilog="Example usage:\n"
               "  stb-mlff -s structure.fdf -c calc.fdf --n-configs 10 --stdev 0.05 0.1\n"
               "  stb-mlff --n-configs 30 --stdev 0.05 0.1 --volume-strains -5 -3 -1 1 3 5\n"
               "  stb-mlff --n-configs 6 --from-aimd aimd --n-samples 6\n"
               "  stb-mlff --n-configs 30 --stdev 0.05 0.1 --volume-strains -3 3 --n-test 8\n"
               "  stb-mlff --file si.fdf --n-configs 10     # one self-contained .fdf instead",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("-s", "--structure", default="structure.fdf",
                        help="Reference structure (default: structure.fdf): lattice, species and "
                             "coordinates, rewritten in every generated configuration.")
    parser.add_argument("-c", "--calc", default="calc.fdf",
                        help="Calculation template (default: calc.fdf), with '%%include "
                             "<structure>' and the basis/mesh cutoff/SCF settings, copied verbatim "
                             "into every folder (a single-point config_extra.fdf is %%included on top).")
    parser.add_argument("--file", default=None,
                        help="Instead of -s/-c: one self-contained reference .fdf (structure and "
                             "settings in the same file), reused verbatim with new coordinates.")
    parser.add_argument("--n-configs", type=int, default=10,
                        help="Number of rattled configurations to generate (default: 10). "
                             "Ignored if --from-aimd is given without also wanting rattled "
                             "configs (see --n-configs 0 to disable rattling entirely).")
    parser.add_argument("--stdev", type=float, nargs="+", default=[0.05],
                        help="One or more rattling displacement standard deviations, in "
                             "Angstrom (default: 0.05). Cycled round-robin across --n-configs.")
    parser.add_argument("--volume-strains", type=float, nargs="+", default=[], metavar="PCT",
                        help="Also write one configuration per volumetric strain, in %% (e.g. "
                             "-5 -3 -1 1 3 5): the cell scaled isotropically by "
                             "(1 + PCT/100)^(1/3), atoms rattled with --stdev. Without them "
                             "every configuration has the reference volume and the fitted "
                             "model cannot learn the volume response (E(V), bulk modulus, "
                             "pressure). Default: none. Bulk only -- in a slab this also "
                             "scales the vacuum.")
    parser.add_argument("--pseudo-dir", default=".",
                        help="Pseudopotentials source: a bundled bank "
                             f"({', '.join(BANKS)}) or a folder with <element>.psf/.psml "
                             "(default: current directory).")
    parser.add_argument("--from-aimd", default=None, metavar="LABEL",
                        help="Also sample starting geometries from an existing AIMD "
                             "trajectory ('<LABEL>.ANI' + '.out'), in addition to (or instead "
                             "of, with --n-configs 0) rattled configurations -- more "
                             "physically realistic diversity than rattling alone.")
    parser.add_argument("--n-samples", type=int, default=5,
                        help="Number of evenly-spaced AIMD frames to sample when --from-aimd "
                             "is given (default: 5).")
    parser.add_argument("--n-test", type=int, default=0,
                        help="Also write N held-out test configurations to "
                             f"{TEST_DIR}/test_NNN/ (default: 0): drawn like the training ones "
                             "(each --stdev, each --volume-strains entry, cycled) with an "
                             "independent seed, never trained on. Stage 3 (stb-mlffValidate) "
                             "scores every round against them -- needed to know when to stop.")
    parser.add_argument("-o", "--output-dir", default=mlff_workflow.RUN_DIR,
                        help="Run folder (default: " + mlff_workflow.RUN_DIR + "): the training configs go to "
                             "<run>/configs_round_01/mlff_config_NNN/ (they train round 1; Stage 4's picks go "
                             "to configs_round_02/, ...), the test set to <run>/mlff_test/, plus the workflow "
                             "file and the next stages' round_NN/ folders.")
    parser.add_argument("--workflow", default=None, metavar="JSON",
                        help=f"Workflow state file shared by the 4 stages (default: "
                             f"<output-dir>/{mlff_workflow.WORKFLOW_FILE}). Created here; an "
                             "existing one is copied to .bak first.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for rattling (default: 42).")
    parser.add_argument("--d3", action="store_true",
                        help="Add Grimme's DFT-D3 dispersion correction (SIESTA's DFTD3 .true., its own "
                             "defaults) to every folder's config_extra.fdf, so the labels -- and the model -- are "
                             "PBE+D3. Recorded in the workflow file: Stage 4's picks and Stage 3's --md-check "
                             "frames get it too, and Stage 2 checks every folder has it. For vdW-bonded "
                             "systems (layers, molecular crystals, physisorption). Default: whatever the "
                             "template says (off if nothing).")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to <run folder>/{REPORT_FILE}. Off by default.")
    parser.add_argument("--ml-screen", action="store_true",
                        help="Before any SIESTA run, evaluate every config with a MACE-MP foundation "
                             "model (needs the 'ml' extra): its max force and energy above the "
                             "reference flag unphysical geometries, and show how far off the base "
                             "model starts. Informational only -- nothing is filtered.")
    parser.add_argument("--ml-model", default="medium",
                        help="Model for --ml-screen: small/medium/large or a .model file (default: medium).")
    parser.add_argument("--ml-device", default="cpu", choices=["cpu", "cuda"],
                        help="Device for --ml-screen (default: cpu).")
    parser.add_argument("--dup-threshold", type=float, default=DUP_THRESHOLD,
                        help="Flag a config whose atomic environments are all closer than this to "
                             "an earlier one (units of the set's typical spacing; default: "
                             f"{DUP_THRESHOLD:g}; 0 disables).")
    parser.add_argument("--export-traj", nargs="?", const="mlff_configs.xsf", default=None, metavar="FILE",
                        help="Also write reference + every config as one multi-frame file in the run "
                             "folder (default name: mlff_configs.xsf; .xsf, .pdb or .xyz -- extended "
                             "XYZ carrying folder/set/kind per frame), e.g. for VESTA/OVITO.")
    parser.add_argument("--no-plots", dest="plots", action="store_false",
                        help="Do not write the coverage plots (gnuplot .dat + .gplot in <run folder>/coverage/).")
    parser.add_argument("--show", action="store_true",
                        help="Also display the coverage figure on screen (it is always saved as "
                             "<run folder>/coverage/mlff_coverage.png, unless --no-plots).")
    parser.add_argument("--view", action="store_true",
                        help="Open the written configurations in ASE's interactive viewer (frame 0 "
                             "= reference, then training, then test configs, in the report's order). "
                             "Needs a display.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-mlff {VERSION}")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - Custom ML Force Field (Prep)",
            "Generate SIESTA training configurations for MACE fine-tuning",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("CUSTOM ML FORCE FIELD -- PREP:", 'bold'))
    print("-" * 60)

    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    if report_path:
        os.makedirs(args.output_dir, exist_ok=True)   # the report lives in the run folder, like Stages 2-4'
    f_out = open(report_path, "w") if report_path else None
    warns = []           # every warning, printed together in the last section
    library_noise = []   # third-party prints/warnings, counted there too

    def warn(msg):
        warns.append(msg)

    def fail(message):
        print_dual(color_text(f"[FAIL] {message}", 'red'), f_out)
        if warns:
            print_section("WARNINGS SO FAR", f_out)
            for k, w in enumerate(warns, 1):
                print_dual(color_text(f"  {k}. {w}", 'yellow'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_dual(color_text("===== STB-MLFF REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    if args.file:
        print_dual(f"Reference file    : {args.file}", f_out)
    else:
        print_dual(f"Structure file    : {args.structure}", f_out)
        print_dual(f"Calc file         : {args.calc}", f_out)
    print_dual(f"Rattled configs   : {args.n_configs} (stdev {args.stdev})", f_out)
    if args.volume_strains:
        print_dual(f"Volume strains    : {', '.join(f'{x:+g}' for x in args.volume_strains)} %", f_out)
    if args.from_aimd:
        print_dual(f"AIMD samples      : {args.n_samples} (from '{args.from_aimd}.ANI')", f_out)
    print_dual(f"Test configs      : {args.n_test}" + (f" (in {TEST_DIR}/)" if args.n_test else ""), f_out)
    print_dual(f"Run folder        : {args.output_dir}", f_out)
    print_dual(f"Seed              : {args.seed} (test set: {args.seed} + {TEST_SEED_OFFSET})", f_out)
    print_dual(f"Pseudo source     : {args.pseudo_dir}", f_out)
    if report_path:
        print_dual(f"Report file       : {report_path}", f_out)

    print_section("[1] REFERENCE STRUCTURE", f_out)

    ref_file, calc_template, error = resolve_reference(args.structure, args.calc, args.file)
    if error:
        fail(error)

    try:
        ref_structure = structure_io.read_fdf(ref_file)
    except Exception as e:
        fail(f"Could not read '{ref_file}': {e}")
    if not ref_structure.atoms:
        fail(f"No atoms read from '{ref_file}' -- it must hold %block AtomicCoordinatesAndAtomicSpecies.")

    unique_elements = sorted(ref_structure.species)
    pmg_structure = structure_io.to_pymatgen(ref_structure)
    ase_atoms = AseAtomsAdaptor.get_atoms(pmg_structure)
    ref_basename = os.path.basename(ref_file)
    vacuum_axes, ref_pairs = report_reference(ase_atoms, ref_file, calc_template, ref_structure, f_out)
    ref_dmin = min(ref_pairs.values()) if ref_pairs else None
    is_bulk = not any(vacuum_axes)

    try:
        import contextlib
        import io
        with contextlib.redirect_stdout(io.StringIO()):   # its bank banner; [2] names the source
            pseudo_path = resolve_pseudo_source(args.pseudo_dir)
    except ValueError as e:
        fail(str(e))
    pseudos_to_copy, missing = get_required_pseudos(unique_elements, pseudo_path)
    if missing:
        fail(f"Missing pseudopotential(s) for: {', '.join(missing)} (searched '{args.pseudo_dir}').")

    print_section("[2] CALCULATION SETTINGS (inherited by every folder)", f_out)
    texts = []
    for path in ([calc_template] if calc_template else []) + [ref_file]:
        with open(path) as f:
            texts.append(f.read())
    pseudo_label = (f"bundled bank '{args.pseudo_dir}'" if args.pseudo_dir in BANKS else f"'{args.pseudo_dir}'")
    plan_warnings = report_settings(texts, ase_atoms, vacuum_axes, pseudo_label, pseudos_to_copy, f_out)
    system_label = parse_directives(*texts).get("systemlabel", "siesta")
    tmpl_d3 = parse_directives(*texts).get("dftd3")
    if args.d3:
        d3 = {"on": True, "forced": True, "source": "--d3 (config_extra.fdf)"}
        print_dual(f"{'Dispersion':<21}: DFT-D3 forced in every folder (config_extra.fdf, --d3; SIESTA's own "
                   "DFTD3 defaults) -- the model will be PBE+D3", f_out)
        if tmpl_d3 is not None and not mlff_workflow.fdf_true(tmpl_d3):
            print_dual(f"{'':<21}  (the template's 'DFTD3 {tmpl_d3}' is overridden: config_extra.fdf is read "
                       "first)", f_out)
    elif tmpl_d3 is not None and mlff_workflow.fdf_true(tmpl_d3):
        d3 = {"on": True, "forced": False, "source": "template"}
        print_dual(f"{'Dispersion':<21}: DFT-D3 from the template (DFTD3 {tmpl_d3}) -- the model will be PBE+D3",
                   f_out)
    else:
        d3 = {"on": False, "forced": False, "source": "template" if tmpl_d3 is not None else "none"}
        print_dual(f"{'Dispersion':<21}: none (no DFT-D3; --d3 adds it -- for layered or vdW-bonded systems)",
                   f_out)
    d3_block = mlff_workflow.DFTD3_BLOCK if args.d3 else ""

    print_section("[3] SAMPLING PLAN", f_out)
    n_aimd = args.n_samples if args.from_aimd else 0
    plan = []
    if args.n_configs > 0:
        per = [sum(1 for i in range(args.n_configs) if args.stdev[i % len(args.stdev)] == sd) for sd in args.stdev]
        plan.append((["rattle", str(args.n_configs),
                      ", ".join(f"{sd:g} Ang x{n}" for sd, n in zip(args.stdev, per))
                      + f" (RMS shift ~ {', '.join(f'{sd * 3 ** 0.5:.3f}' for sd in args.stdev)} Ang)"], None))
    if args.volume_strains:
        plan.append((["strain", str(len(args.volume_strains)),
                      "volume " + ", ".join(f"{x:+g}" for x in args.volume_strains) + " % (+ rattle)"], None))
    if args.from_aimd:
        plan.append((["aimd", f"<= {n_aimd}", f"evenly spaced frames of {args.from_aimd}.ANI"], None))
    if args.n_test:
        plan.append((["test", str(args.n_test), "same kinds, cycled, independent seed"], None))
    if plan:
        print_table(["Kind", "Count", "Parameters"], plan, f_out)
    if ref_dmin:
        for sd in args.stdev:
            ratio = sd * 3 ** 0.5 / ref_dmin
            if ratio > 0.15:
                plan_warnings.append(f"--stdev {sd:g}: RMS shift {sd * 3 ** 0.5:.3f} Ang is {100 * ratio:.0f} % of "
                                     f"the shortest bond ({ref_dmin:.3f} Ang) -- expect very close atoms "
                                     "(see the d_min column) and hard SCF convergence.")
            elif ratio < 0.01:
                plan_warnings.append(f"--stdev {sd:g}: RMS shift is < 1 % of the shortest bond -- forces "
                                     "will be tiny and nearly harmonic; add a larger amplitude too.")
    if is_bulk and not args.volume_strains and not args.from_aimd:
        plan_warnings.append("No --volume-strains: every config has the reference volume, so the model "
                             "cannot learn E(V), the bulk modulus or pressure.")
    if not is_bulk and args.volume_strains:
        print_dual("Volume strains    : periodic axes only (the vacuum keeps its size) -- an area strain "
                   "for a 2D cell", f_out)
    if args.n_test == 0:
        plan_warnings.append("No --n-test: Stage 3 (stb-mlffValidate) needs a held-out test set "
                             "to say whether a round helped.")
    elif args.n_test < 5:
        plan_warnings.append(f"Only {args.n_test} test config(s): Stage 3's verdict will be noisy (>= 5 advised).")
    if 0 < args.n_configs + len(args.volume_strains) + n_aimd < 20:
        plan_warnings.append(f"{args.n_configs + len(args.volume_strains) + n_aimd} training configs = as many "
                             "energy labels: enough for a quick test, usually not for a production model "
                             "(20-50+ is typical for one phase).")
    if args.n_test:
        n_kinds = len(args.stdev) + len(args.volume_strains)
        eos_ok = len(args.volume_strains) >= 5 and args.n_test >= n_kinds
        print_dual("Test-set E(V)     : " + ("yes -- Stage 3 fits it on the test strains" if eos_ok else
                   f"no -- Stage 3 needs 5+ --volume-strains and --n-test >= {max(n_kinds, len(args.stdev) + 5)} "
                   "(the test kinds are cycled, one config each)"), f_out)
    for w in plan_warnings:
        warn(f"Plan: {w}")
    print_dual(color_text(f"Plan              : {len(plan_warnings)} issue(s) -- see [9]", 'yellow') if plan_warnings
               else color_text("[OK] No issues found in the plan.", 'green'), f_out)

    os.makedirs(args.output_dir, exist_ok=True)

    if args.n_test < 0:
        fail(f"--n-test must be >= 0 (got {args.n_test}).")
    if args.export_traj and os.path.splitext(args.export_traj)[1].lstrip(".").lower() not in TRAJ_FORMATS:
        fail(f"--export-traj must end in {', '.join('.' + k for k in TRAJ_FORMATS)} (got '{args.export_traj}').")
    if args.ml_screen:
        from stb.core.mace_relax import resolve_device
        try:
            resolve_device(args.ml_device)
        except ValueError as e:
            fail(str(e))
    workflow_path = args.workflow or os.path.join(args.output_dir, mlff_workflow.WORKFLOW_FILE)
    wf_base = mlff_workflow.base_dir(workflow_path)

    print_section("[4] GENERATING CONFIGURATIONS", f_out)
    folders = []
    entries = []
    index = 0
    n_ref = len(ref_structure.atoms)
    # Training configs in <run>/configs_round_01/ (they train round 1; Stage 4's picks go to
    # configs_round_02/, ...), the test set in <run>/mlff_test/.
    train_root = mlff_workflow.configs_dir(args.output_dir, 1)
    os.makedirs(train_root, exist_ok=True)
    before = set(glob.glob(os.path.join(train_root, "mlff_config_*")) +
                 glob.glob(os.path.join(args.output_dir, TEST_DIR, "test_*")))

    def _entry(folder, kind, origin, which="train"):
        entries.append({"folder": folder, "set": which, "kind": kind, "origin": origin,
                        "stage": "prep", "round": 0, "n_atoms": n_ref})

    if args.n_configs > 0:
        rattled = build_rattled_configs(
            ase_atoms, ref_structure.coord_format, args.n_configs, args.stdev, args.seed)
        for i, positions in enumerate(rattled):
            index += 1
            folder = write_config_folder(
                train_root, index, ref_file, ref_basename, positions, None, pseudos_to_copy,
                extra_block=d3_block, calc_template=calc_template)
            folders.append(folder)
            _entry(folder, "rattle", f"rattle {args.stdev[i % len(args.stdev)]:g} Ang")
        print_dual(f"[OK] Generated {len(rattled)} rattled configuration(s).", f_out)

    if args.volume_strains:
        bad = [x for x in args.volume_strains if x <= -100.0]
        if bad:
            fail(f"--volume-strains must be > -100 % (got {bad}).")
        strained = build_strained_configs(
            ase_atoms, ref_structure.coord_format,
            structure_io.raw_lattice_vectors(ref_structure),
            args.volume_strains, args.stdev, args.seed)
        for i, (positions, raw_lattice) in enumerate(strained):
            index += 1
            folder = write_config_folder(
                train_root, index, ref_file, ref_basename, positions, raw_lattice, pseudos_to_copy,
                extra_block=d3_block, calc_template=calc_template)
            folders.append(folder)
            _entry(folder, "strain", f"volume {args.volume_strains[i]:+g} % "
                                     f"(rattle {args.stdev[i % len(args.stdev)]:g} Ang)")
        print_dual(f"[OK] Generated {len(strained)} strained configuration(s) "
                   f"(volume {', '.join(f'{x:+g}' for x in args.volume_strains)} %).", f_out)

    if args.from_aimd:
        try:
            aimd_configs = sample_aimd_configs(
                args.from_aimd, ref_structure.coord_format, ref_structure.lattice_constant,
                args.n_samples, f_out)
        except Exception as e:
            fail(f"Could not sample from AIMD trajectory '{args.from_aimd}': {e}")
        for i, (positions, raw_lattice) in enumerate(aimd_configs):
            index += 1
            folder = write_config_folder(
                train_root, index, ref_file, ref_basename, positions, raw_lattice, pseudos_to_copy,
                extra_block=d3_block, calc_template=calc_template)
            folders.append(folder)
            _entry(folder, "aimd", f"AIMD frame {i + 1}/{len(aimd_configs)} of '{args.from_aimd}'")
        print_dual(f"[OK] Generated {len(aimd_configs)} AIMD-sampled configuration(s).", f_out)

    if not folders:
        fail("No configurations were generated (--n-configs 0, no --volume-strains and no --from-aimd).")

    test_folders = []
    if args.n_test:
        tests = build_test_configs(
            ase_atoms, ref_structure.coord_format, structure_io.raw_lattice_vectors(ref_structure),
            args.n_test, args.stdev, args.volume_strains, args.seed)
        test_root = os.path.join(args.output_dir, TEST_DIR)
        os.makedirs(test_root, exist_ok=True)
        for i, (positions, raw_lattice, kind, origin) in enumerate(tests, start=1):
            folder = write_config_folder(test_root, i, ref_file, ref_basename, positions,
                                         raw_lattice, pseudos_to_copy, prefix="test",
                                         extra_block=d3_block, calc_template=calc_template)
            test_folders.append(folder)
            _entry(folder, kind, origin, which="test")
        print_dual(f"[OK] Generated {len(test_folders)} held-out test configuration(s) in "
                   f"{test_root}/ (never trained on).", f_out)

    print_dual(f"[OK] Each folder forces a single-point run via {CONFIG_EXTRA_FILE} "
               "(no relaxation, no cell change)" + (" and DFT-D3 (--d3)" if args.d3 else "") + ".", f_out)
    stale = sorted(before - {os.path.normpath(f) for f in folders + test_folders}
                   - set(folders + test_folders))
    if stale:
        warn(f"Folders: {len(stale)} folder(s) from an earlier Stage 1 run were not rewritten and are not in the new "
             f"workflow file ({', '.join(os.path.relpath(d) for d in stale[:4])}{' ...' if len(stale) > 4 else ''}) "
             "-- delete them so nobody runs SIESTA there for nothing.")

    bak = mlff_workflow.backup(workflow_path)
    if bak:
        warn(f"Workflow: {workflow_path} already existed -- copied to {bak} and started fresh.")
    state = mlff_workflow.new_state(ref_file, n_ref, unique_elements, wf_base)
    if calc_template:
        state["reference"]["calc"] = mlff_workflow.rel(calc_template, wf_base)
    # A folder is stored relative to the workflow file; a bundled bank name
    # (dojo, virtual_vault) as given.
    state["reference"]["pseudo_dir"] = (args.pseudo_dir if args.pseudo_dir in BANKS
                                        else mlff_workflow.rel(pseudo_path, wf_base))
    state["reference"]["dftd3"] = d3
    mlff_workflow.add_configs(state, entries, wf_base)
    mlff_workflow.record(state, "prep", "stb-mlff", sys.argv[1:])
    mlff_workflow.save(state, workflow_path)

    # ---- per-configuration geometry, read back from what was written
    print_section("[5] CONFIGURATIONS", f_out)
    rows, stats, frames = [], [], [ase_atoms]
    for e in entries:
        path = os.path.join(e["folder"], ref_basename)
        try:
            atoms = fdf_atoms(path)
        except Exception as err:
            fail(f"Could not read back {path}: {err}")
        st = config_stats(ase_atoms, atoms, ref_dmin)
        st.update(e)
        stats.append(st)
        atoms.info["folder"], atoms.info["kind"] = e["folder"], e["kind"]
        frames.append(atoms)
        color = None
        if st["d_ratio"] is not None and st["d_ratio"] < CLOSE_FAIL:
            color = 'red'
        elif st["d_ratio"] is not None and st["d_ratio"] < CLOSE_WARN:
            color = 'yellow'
        rows.append(([str(len(frames) - 1), mlff_workflow.rel(e["folder"], args.output_dir), e["set"], e["kind"],
                      e["origin"], f"{st['volume']:.2f}", f"{st['dV_pct']:+.2f}",
                      f"{st['u_rms']:.3f}", f"{st['u_max']:.3f}",
                      f"{st['d_min']:.3f}" if st["d_min"] else "-"], color))
    nov_unit = config_novelty(stats, frames[1:]) if args.dup_threshold > 0 else None
    ml = None
    if args.ml_screen:
        print_dual(f"[INFO] Foundation-model screen ({args.ml_model}, {args.ml_device}) of "
                   f"{len(frames)} structures ...", f_out)
        try:
            with capture_library_noise(library_noise, "foundation-model screen"):
                ml = ml_screen(frames, args.ml_model, args.ml_device)
        except SystemExit:
            raise
        except Exception as err:
            warn(f"ML screen: --ml-screen failed ({err}) -- skipped.")
        if ml:
            for st, r in zip(stats, ml[1:]):
                st["ml_f_max"] = r["f_max"]
                st["ml_dE"] = 1000.0 * (r["e_atom"] - ml[0]["e_atom"])
                st["ml_p"] = r["p_gpa"]
    for i, ((cells, color), st) in enumerate(zip(rows, stats)):
        nov = st.get("novelty")
        cells.append("-" if nov is None else f"{nov:.2f}")
        if nov is not None and nov < args.dup_threshold and color is None:
            color = 'yellow'
        if ml:
            cells += [f"{st['ml_dE']:+.0f}", f"{st['ml_f_max']:.2f}"]
            if st["ml_f_max"] > ML_FORCE_WARN or st["ml_dE"] > ML_ENERGY_WARN:
                color = 'red'
        rows[i] = (cells, color)
    headers = ["#", "Folder", "Set", "Kind", "Detail", "V (Ang^3)", "dV %", "u_rms", "u_max", "d_min", "nov"]
    if ml:
        headers += ["dE_ML", "Fmax_ML"]
    print_table(headers, rows, f_out)
    print_dual(f"u = displacement from the reference (Ang, nearest image, in the config's own cell); "
               f"d_min = shortest interatomic distance (Ang; reference {ref_dmin:.3f})." if ref_dmin else
               "u = displacement from the reference (Ang).", f_out)
    if nov_unit:
        print_dual(f"nov = novelty vs the training configs before it (a test config: vs all training + earlier "
                   f"test), in units of the set's typical environment spacing ({nov_unit:.3f} Ang); "
                   f"< {args.dup_threshold:g} = near-duplicate.", f_out)
    if ml:
        print_dual(f"dE_ML = foundation-model ({args.ml_model}) energy above the reference, meV/atom; "
                   "Fmax_ML = its largest atomic force, eV/Ang (before any SIESTA run).", f_out)
    print_dual(f"Yellow: d_min < {CLOSE_WARN:.0%} of the reference's or a near-duplicate; red: d_min < "
               f"{CLOSE_FAIL:.0%}" + (f", Fmax_ML > {ML_FORCE_WARN:g} eV/Ang or dE_ML > {ML_ENERGY_WARN:g} "
                                      "meV/atom" if ml else "") + ".", f_out)

    print_section("[6] COVERAGE & CHECKS", f_out)
    check_warnings = []

    def _rng(vals, fmt):
        return f"{fmt.format(min(vals))} .. {fmt.format(max(vals))}" if vals else "-"

    cov = []
    for which in ("train", "test"):
        sub = [x for x in stats if x["set"] == which]
        if not sub:
            continue
        cov.append(([which, str(len(sub)), _rng([x["u_rms"] for x in sub], "{:.3f}"),
                     _rng([x["u_max"] for x in sub], "{:.3f}"), _rng([x["dV_pct"] for x in sub], "{:+.2f}"),
                     _rng([x["d_min"] for x in sub if x["d_min"]], "{:.3f}")], None))
    print_table(["Set", "N", "u_rms (Ang)", "u_max (Ang)", "dV (%)", "d_min (Ang)"], cov, f_out)
    train = [x for x in stats if x["set"] == "train"]
    test = [x for x in stats if x["set"] == "test"]
    n_lab_atoms = sum(x["n_atoms"] for x in train)
    print_dual(f"Training labels   : {len(train)} energies, {3 * n_lab_atoms} force components"
               + (f", {6 * len(train)} stress components" if is_bulk else " (stress: only meaningful along periodic axes)"), f_out)
    if train and test:
        lo_v, hi_v = min(x["dV_pct"] for x in train), max(x["dV_pct"] for x in train)
        hi_u = max(x["u_max"] for x in train)
        out = [x for x in test if x["dV_pct"] < lo_v - 0.05 or x["dV_pct"] > hi_v + 0.05 or x["u_max"] > 1.2 * hi_u]
        if out:
            check_warnings.append(f"{len(out)} test config(s) lie outside the training range (volume or "
                                  "displacement): Stage 3 will measure extrapolation there.")
        else:
            print_dual(color_text("[OK] Every test config lies inside the training range (volume and displacement).", 'green'), f_out)
    close = [x for x in stats if x["d_ratio"] is not None and x["d_ratio"] < CLOSE_FAIL]
    near = [x for x in stats if x["d_ratio"] is not None and CLOSE_FAIL <= x["d_ratio"] < CLOSE_WARN]
    if close:
        check_warnings.append(f"{len(close)} config(s) have atoms closer than {CLOSE_FAIL:.0%} of the reference's "
                              f"shortest distance ({', '.join(os.path.basename(x['folder']) for x in close)}): "
                              "SIESTA may not converge, and such geometries are far outside what MD visits -- "
                              "consider a smaller --stdev.")
    if near:
        check_warnings.append(f"{len(near)} config(s) have a short contact ({CLOSE_FAIL:.0%}-{CLOSE_WARN:.0%} of the "
                              "reference's shortest distance) -- fine for robustness, check SCF convergence.")
    if not close and not near:
        print_dual(color_text(f"[OK] No close contacts (every d_min >= {CLOSE_WARN:.0%} of the reference's).", 'green'), f_out)
    if nov_unit:
        dups = [x for x in stats if x.get("novelty") is not None and x["novelty"] < args.dup_threshold]
        dup_test = [x for x in dups if x["set"] == "test"]
        if dups:
            check_warnings.append(
                f"{len(dups)} near-duplicate config(s) (nov < {args.dup_threshold:g}: "
                f"{', '.join(os.path.basename(x['folder']) for x in dups)}) -- extra SIESTA runs that add "
                "little; " + ("a test config this close to a training one does not test generalization; "
                              if dup_test else "") + "use more --stdev values or a different --seed.")
        else:
            novs = [x["novelty"] for x in stats if x.get("novelty") is not None]
            print_dual(color_text(f"[OK] No near-duplicates (novelty {min(novs):.2f} .. {max(novs):.2f})."
                                  if novs else "[OK] No near-duplicates.", 'green'), f_out)
    elif args.dup_threshold > 0:
        print_dual("Near-duplicate check: skipped (fewer than 2 configs).", f_out)
    if ml:
        r0 = ml[0]
        print_dual(f"Foundation model  : {args.ml_model} -- reference max force {r0['f_max']:.3f} eV/Ang"
                   + (f", pressure {r0['p_gpa']:+.2f} GPa" if (r0['p_gpa'] is not None and is_bulk) else ""), f_out)
        if r0["f_max"] > 0.5:
            check_warnings.append(f"The foundation model sees {r0['f_max']:.2f} eV/Ang on the REFERENCE itself: "
                                  "either the reference is not relaxed, or the base model disagrees with it -- "
                                  "the fine-tuning will have to fix that.")
        bad = [x for x in stats if x["ml_f_max"] > ML_FORCE_WARN or x["ml_dE"] > ML_ENERGY_WARN]
        if bad:
            check_warnings.append(f"{len(bad)} config(s) look unphysical to the foundation model "
                                  f"(Fmax > {ML_FORCE_WARN:g} eV/Ang or dE > {ML_ENERGY_WARN:g} meV/atom: "
                                  f"{', '.join(os.path.basename(x['folder']) for x in bad)}) -- check them "
                                  "(--view) before spending SIESTA time; a smaller --stdev usually fixes it.")
        else:
            print_dual(color_text(f"[OK] Foundation-model max force {min(x['ml_f_max'] for x in stats):.2f} .. "
                                  f"{max(x['ml_f_max'] for x in stats):.2f} eV/Ang, dE "
                                  f"{min(x['ml_dE'] for x in stats):+.0f} .. {max(x['ml_dE'] for x in stats):+.0f} "
                                  "meV/atom: no unphysical geometry.", 'green'), f_out)
    for w in check_warnings:
        warn(f"Checks: {w}")
    if check_warnings:
        print_dual(color_text(f"Checks            : {len(check_warnings)} issue(s) -- see [9]", 'yellow'), f_out)

    # ---- optional outputs: coverage plots, multi-frame file; per-config numbers into the workflow file
    plot_files = []
    if args.plots:
        try:
            plot_files = write_coverage_plots(os.path.join(args.output_dir, "coverage"), ase_atoms,
                                              frames[1:], stats, ref_dmin, bool(ml))
        except Exception as err:
            warn(f"Plots: coverage plots not written: {err}")
        show = args.show
        if show:
            import matplotlib
            if matplotlib.get_backend().lower() in ("agg", "pdf", "ps", "svg", "pgf", "cairo", "template"):
                warn("Plots: --show needs a display (no interactive matplotlib backend here) -- the figure is only saved.")
                show = False
        try:
            png_path = os.path.join(args.output_dir, "coverage", "mlff_coverage.png")
            plot_coverage_png(png_path, ase_atoms, frames[1:], stats, ref_dmin, bool(ml), show)
            plot_files.append(png_path)
        except Exception as err:
            warn(f"Plots: coverage figure not written: {err}")
    traj_path = None
    if args.export_traj:
        from ase.io import write as ase_write
        traj_path = (args.export_traj if os.path.dirname(args.export_traj)
                     else os.path.join(args.output_dir, args.export_traj))
        ext = os.path.splitext(traj_path)[1].lstrip(".").lower()
        out_frames = []
        for k, a in enumerate(frames):
            b = a.copy()
            b.info = {"frame": k, "folder": "reference" if k == 0 else mlff_workflow.rel(stats[k - 1]["folder"], args.output_dir),
                      "set": "reference" if k == 0 else stats[k - 1]["set"],
                      "kind": "reference" if k == 0 else stats[k - 1]["kind"]}
            out_frames.append(b)
        try:
            ase_write(traj_path, out_frames, format=TRAJ_FORMATS[ext])
        except Exception as err:
            warn(f"Export: could not write {traj_path}: {err}")
            traj_path = None
    by_folder = {mlff_workflow.rel(x["folder"], wf_base): x for x in stats}
    for c in state["configs"]:
        x = by_folder.get(c["folder"])
        if x is None:
            continue
        c["geometry"] = {k: (round(x[k], 5) if isinstance(x[k], float) else x[k])
                         for k in ("volume", "dV_pct", "u_rms", "u_max", "d_min", "novelty")}
        if ml:
            c["ml_screen"] = {"model": args.ml_model, "f_max": round(x["ml_f_max"], 4),
                              "dE_meV_atom": round(x["ml_dE"], 2)}
    mlff_workflow.save(state, workflow_path)

    print_section("[7] SUMMARY & FILES", f_out)
    print_dual("Status            : OK" + (f" ({len(warns)} warning(s), see [9])" if warns else ""), f_out)
    print_dual(f"Configurations    : {len(folders)} training" +
               (f" + {len(test_folders)} test" if test_folders else "") +
               f" = {len(folders) + len(test_folders)} SIESTA single points, "
               f"{sum(x['n_atoms'] for x in stats)} atoms in total", f_out)
    print_dual(f"Run folder        : {args.output_dir}/", f_out)
    print_dual(f"  {os.path.relpath(train_root, args.output_dir)}/mlff_config_001 .. {os.path.basename(folders[-1])}"
               "  (training, round 1)", f_out)
    if test_folders:
        print_dual(f"  {TEST_DIR}/test_001 .. {os.path.basename(test_folders[-1])}  (held-out test)", f_out)
    main_fdf = os.path.basename(calc_template) if calc_template else ref_basename
    files = [main_fdf] + ([ref_basename] if calc_template else []) + [CONFIG_EXTRA_FILE] + \
            [os.path.basename(p) for p in pseudos_to_copy]
    print_dual(f"Each folder       : {', '.join(files)}", f_out)
    print_dual(f"Workflow file     : {workflow_path} (+ each config's geometry numbers"
               + (" and ML screen" if ml else "") + ")", f_out)
    gp = [p for p in plot_files if p.endswith(".gplot")]
    png = [p for p in plot_files if p.endswith(".png")]
    if png:
        print_dual(f"Coverage figure   : {png[0]}", f_out)
    if gp:
        print_dual(f"Coverage plots    : {', '.join(gp)}", f_out)
        print_dual(f"                    (cd {os.path.dirname(gp[0])} && gnuplot "
                   f"{' '.join(os.path.basename(p) for p in gp)}  -> .pdf)", f_out)
    if traj_path:
        print_dual(f"All structures    : {traj_path} ({len(frames)} frames: 0 = reference, then the # column)", f_out)
    if report_path:
        print_dual(f"Report            : {report_path}", f_out)

    print_section("[8] NEXT STEPS", f_out)
    globs = f"{train_root}/mlff_config_*" + (f" {args.output_dir}/{TEST_DIR}/test_*" if test_folders else "")
    run_here = [f"{train_root}/ ({len(folders)} folders)"] + (
        [f"{os.path.join(args.output_dir, TEST_DIR)}/ ({len(test_folders)} folders)"] if test_folders else [])
    warn(f"SIESTA: run it in {' and '.join(run_here)} before Stage 2 (stb-mlffAnalysis).")
    print_dual(color_text(f"Run SIESTA in     : {' and '.join(run_here)}", 'cyan'), f_out)
    print_dual("1. Run SIESTA in every folder (serial example; on a cluster, one job per folder):", f_out)
    print_dual(f'     for d in {globs}; do (cd "$d" && siesta < {main_fdf} > {system_label}.out); done', f_out)
    print_dual(f"   Each run writes {system_label}.out/.FA/.XV -- Stage 2 needs all three.", f_out)
    print_dual("2. Train a committee: stb-mlffAnalysis (4.17.2) -- it finds the configs through the workflow file.", f_out)
    print_dual("3. Score it on the test set: stb-mlffValidate (4.17.3); then stb-mlffActiveLearning (4.17.4).", f_out)
    tips = [t for t, used in (("--view (ASE viewer)", args.view), ("--export-traj (one multi-frame file)", args.export_traj),
                              ("--ml-screen (foundation-model check before SIESTA)", args.ml_screen),
                              ("--show (display the coverage figure)", args.show)) if not used]
    if tips:
        print_dual(f"Also available    : {'; '.join(tips)}.", f_out)

    print_section("[9] WARNINGS", f_out)
    if warns:
        for k, w in enumerate(warns, 1):
            print_dual(color_text(f"{k:>2}. {w}", 'yellow'), f_out)
    else:
        print_dual(color_text("None.", 'green'), f_out)
    if library_noise:
        print_dual(f"Library messages  : {len(library_noise)} informational message(s) from mace/torch captured "
                   "(model loading, dtype).", f_out)

    if f_out:
        f_out.close()

    print("\n" + "-" * 60)
    print(color_text(f"Run SIESTA in each {train_root}/mlff_config_NNN/ folder" +
                     (f" and in {args.output_dir}/{TEST_DIR}/test_NNN/" if test_folders else "") +
                     ", then train with stb-mlffAnalysis (4.17.2).\n", 'bold'))

    if args.view:
        from stb.core.ase_view import view_structure_interactive
        print(color_text(f"[INFO] Opening ASE viewer: {len(frames)} frames (0 = reference, then the "
                         "# column of [5]).", 'cyan'))
        view_structure_interactive(frames)


if __name__ == "__main__":
    main()
