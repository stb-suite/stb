#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""Climbing-image NEB (nudged elastic band) driven entirely by a MACE potential: the barrier of a
transition between two relaxed endpoint structures, with no SIESTA. The pure-MACE counterpart of
the SIESTA NEB workflow (stb-neb, 5.9): for a quick ML-level estimate -- a vacancy or interstitial
hop, an adsorbate moving between sites, a small rearrangement -- or to prepare a good band before
the DFT NEB.

Pipeline: both endpoints are relaxed with the model (positions only); the band is interpolated
linearly between them (atoms matched by distance, --autosort-tol) and smoothed with IDPP; a
two-stage climbing-image NEB relaxes it (core/mace_relax.relax_neb: the band first finds its
shape, then the highest image climbs to the saddle); the energy profile is fitted with ASE's
force-consistent spline (ase.utils.forcecurve.fit_images).

The report, in the style of the other ML tools (numbered sections, every warning in a last
section, library noise only counted): run settings and device, the structures and the atoms
that move between the endpoints, the model card in full, the endpoint relaxations, the starting
path (step between images, shortest distance per image), the NEB convergence, the band image by
image, the barrier (spline and highest image, forward and backward), the saddle geometry, the
band against the card's trained range, and -- with a fine-tuned model -- the same calculation
with the foundation model. --show displays the figures, --save-gnuplot writes .dat + .gplot into
<output-dir>/plot/ (and the PDFs when gnuplot is installed), --save-data the band as a table,
--save-images every relaxed image as an .fdf.

The geometry helpers are neb.py's own (composition check, wrapping, interpolation, IDPP, frozen
atoms), so a path here and a path written for SIESTA by stb-neb are built the same way. Unlike
stb-neb, a lattice mismatch between the endpoints is not fatal: the band adopts the initial
lattice and the report warns.
"""

from stb import __version__ as VERSION

import os
import re
import sys
import time
import shutil
import argparse
import contextlib
from datetime import datetime

import numpy as np

from stb.core import structure_io, model_card
from stb.core.cli import (color_text, show_intro, print_dual, print_table, capture_library_noise,
                          display_available, Sections, run_gnuplot)
from stb.core.deps import require_mace

REPORT_FILE = "stb_mlneb_report.txt"
PLOT_DIR = "plot"
MOVER_THRESHOLD = 0.1          # Ang: an atom that moves more than this between the endpoints is listed
SPLINE_VS_IMAGE_TOL = 0.05     # eV: spline and highest-image barriers further apart -> the band is too coarse
E_BELOW_TRAINED_TOL = 0.02     # eV/atom below the lowest training energy (custom model) -> flagged
_COLORS = ["#2a6fdb", "#d9480f", "#2f9e44", "#862e9c"]


def resolve_lattice_mismatch(initial_pmg, final_pmg, f_out=None, tol=1e-3, warn=None):
    """`final_pmg` rebuilt on the INITIAL lattice (species and fractional coordinates kept).
    ase.mep.neb.NEB and pymatgen's interpolation both need one exact cell for every image (no
    variable-cell NEB in ASE); below `tol` (noise between two separately relaxed endpoints) the
    rebuild is silent, above it a warning says the final lattice was discarded. The warning goes
    to `warn(message)` when given, else it is printed."""
    from pymatgen.core import Structure
    d = float(np.abs(np.array(initial_pmg.lattice.matrix) - np.array(final_pmg.lattice.matrix)).max())
    if d > tol:
        msg = (f"Initial and final structures have different lattices (largest component difference {d:.4f} "
               "Ang). ase.mep.neb.NEB and pymatgen's interpolation both require every image to share one exact "
               "cell, so the INITIAL structure's lattice is adopted for the whole band; the final structure keeps "
               "its fractional positions and loses its own lattice. Relax both endpoints in the same cell.")
        if warn:
            warn(msg)
        else:
            print_dual(color_text(f"[WARNING] {msg}", 'yellow'), f_out)
    return Structure(initial_pmg.lattice, final_pmg.species, final_pmg.frac_coords, coords_are_cartesian=False)


def _quiet(noise, label):
    return capture_library_noise(noise, label) if noise is not None else contextlib.nullcontext()


def max_force(atoms):
    return float(np.linalg.norm(atoms.get_forces(), axis=1).max())


def relax_endpoint(label, pmg, calc, args, noise=None):
    """One endpoint relaxed with the model (positions only). Returns a dict with its energy and
    forces before/after, the optimizer's verdict and how far the atoms moved."""
    from pymatgen.io.ase import AseAtomsAdaptor
    from stb.core import mace_relax
    at = AseAtomsAdaptor.get_atoms(pmg)
    at.calc = calc
    with _quiet(noise, label):
        e0, f0 = at.get_potential_energy(), max_force(at)
    pos0 = at.get_positions().copy()
    converged, steps = None, 0
    if args.pre_relax_endpoints:
        with _quiet(noise, label):
            converged, steps = mace_relax.relax(at, calc, fmax=args.fmax, max_steps=args.max_steps)
    with _quiet(noise, label):
        e1, f1 = at.get_potential_energy(), max_force(at)
    return {"label": label, "e0": e0, "f0": f0, "e1": e1, "f1": f1, "converged": converged, "steps": steps,
            "moved": float(np.linalg.norm(at.get_positions() - pos0, axis=1).max()),
            "pmg": AseAtomsAdaptor.get_structure(at)}


def neb_band(initial_pmg, final_pmg, calc, args, noise=None, label="model"):
    """The whole pipeline on one (initial, final) pair, both wrapped into [0, 1) and on one
    lattice: endpoint relaxations, interpolation (+ IDPP), climbing-image NEB, spline fit. Prints
    nothing; returns a dict with everything the report needs.

    `args` needs n_images, autosort_tol, idpp, idpp_fmax, idpp_steps, pre_relax_endpoints,
    freeze_substrate, freeze_threshold, k, fmax, max_steps (any argparse Namespace that carries
    them)."""
    from pymatgen.io.ase import AseAtomsAdaptor
    from ase.constraints import FixAtoms
    from ase.utils.forcecurve import fit_images
    from stb.core import mace_relax
    from stb.neb import linear_interpolate_images, idpp_refine_images, compute_frozen_indices

    t0 = time.monotonic()
    ends = [relax_endpoint("initial", initial_pmg, calc, args, noise),
            relax_endpoint("final", final_pmg, calc, args, noise)]
    # The interpolation takes fractional coordinates literally (stb.neb.linear_interpolate_images, pbc=False). Two
    # endpoints wrapped into the cell separately can hold one atom on opposite faces (0.98 and 0.01): it would cross
    # the whole cell. Every atom of the final endpoint goes to the periodic image nearest to its initial position.
    shift = np.round(ends[1]["pmg"].frac_coords - ends[0]["pmg"].frac_coords)
    crossed = int(np.any(shift != 0, axis=1).sum())
    if crossed:
        from pymatgen.core import Structure
        fin = ends[1]["pmg"]
        ends[1]["pmg"] = Structure(fin.lattice, fin.species, fin.frac_coords - shift)
    with _quiet(noise, label):
        guess = linear_interpolate_images(ends[0]["pmg"], ends[1]["pmg"], args.n_images,
                                          autosort_tol=args.autosort_tol)
        if args.idpp:
            guess = idpp_refine_images(guess, fmax=args.idpp_fmax, steps=args.idpp_steps)
    guess_atoms = [AseAtomsAdaptor.get_atoms(s) for s in guess]
    guess_dmin = [model_card.shortest_distance(a) for a in guess_atoms]
    guess_steps = [float(np.linalg.norm(guess_atoms[i + 1].positions - guess_atoms[i].positions))
                   for i in range(len(guess_atoms) - 1)]

    images = [a.copy() for a in guess_atoms]
    frozen = []
    if args.freeze_substrate:
        frozen = compute_frozen_indices(guess, threshold=args.freeze_threshold)
        if frozen:
            for image in images:
                image.set_constraint(FixAtoms(indices=frozen))
    history = []
    with _quiet(noise, label):
        converged, s1, s2, energies = mace_relax.relax_neb(
            images, calc, k=args.k, fmax=args.fmax, max_steps=args.max_steps,
            on_step=lambda stage, step, res: history.append((stage, step, res)))
        fit = fit_images(images)
        forces = [max_force(a) for a in images]
        e_abs = np.array([a.get_potential_energy() for a in images], float)
    energies = np.array(fit.energies, float)
    path = np.array(fit.path, float)
    e_fit = np.array(fit.fit_energies, float)
    s_fit = np.array(fit.fit_path, float)
    k_top = int(np.argmax(e_fit))
    interior = list(range(1, len(images) - 1))
    i_top = max(interior, key=lambda i: energies[i]) if interior else int(np.argmax(energies))
    ref = images[0].get_positions()
    disp = [float(np.linalg.norm(a.get_positions() - ref, axis=1).max()) for a in images]
    return {"label": label, "ends": ends, "crossed": crossed, "guess": guess, "guess_dmin": guess_dmin, "guess_steps": guess_steps,
            "frozen": frozen, "images": images, "converged": bool(converged), "s1": s1, "s2": s2,
            "history": history, "residual": history[-1][2] if history else None,
            "energies": energies, "e_abs": e_abs, "path": path, "s_fit": s_fit, "e_fit": e_fit, "forces": forces,
            "dmin": [model_card.shortest_distance(a) for a in images], "disp": disp,
            "barrier": float(e_fit[k_top] - energies[0]), "barrier_img": float(energies[i_top] - energies[0]),
            "dE": float(energies[-1] - energies[0]), "s_saddle": float(s_fit[k_top]), "i_top": i_top,
            "seconds": time.monotonic() - t0}


# ---- analysis ---------------------------------------------------------------------------

def movers(initial_atoms, final_atoms, threshold=MOVER_THRESHOLD):
    """(index, symbol, displacement) of every atom that moves more than `threshold` between the
    two endpoints (minimum image), largest first."""
    from ase.geometry import find_mic
    d, _ = find_mic(final_atoms.get_positions() - initial_atoms.get_positions(), initial_atoms.cell,
                    initial_atoms.pbc)
    n = np.linalg.norm(d, axis=1)
    sym = initial_atoms.get_chemical_symbols()
    return [(int(i), sym[i], float(n[i])) for i in np.argsort(-n) if n[i] > threshold], float(n.sum())


def neighbor_distances(atoms, index, n=4):
    """The `n` shortest distances (Ang, minimum image) from atom `index`."""
    d = atoms.get_distances(index, list(range(len(atoms))), mic=True)
    d[index] = np.inf
    return np.sort(d)[:n]


def domain_rows(res, card, warn):
    """The band against the model card's trained range: rows for print_table, warnings to warn()."""
    dom = (card or {}).get("domain") or {}
    rows = []
    vol = res["images"][0].get_volume() / len(res["images"][0])
    m = dom.get("measure")
    if m and (dom.get("measure_name") or "").startswith("volume"):
        ok = m["min"] * 0.98 <= vol <= m["max"] * 1.02
        rows.append((["volume per atom (Ang^3)", f"{m['min']:.3f} .. {m['max']:.3f}", f"{vol:.3f}",
                      "inside" if ok else "OUTSIDE"], None if ok else "yellow"))
    dd = dom.get("d_min_ang")
    dmin = min(d for d in res["dmin"] if d is not None)
    if dd:
        ok = dmin >= 0.95 * dd["min"]
        rows.append((["shortest distance (Ang)", f">= {dd['min']:.3f}", f"min {dmin:.3f} (image "
                      f"{int(np.argmin([d if d is not None else np.inf for d in res['dmin']]))})",
                      "inside" if ok else "OUTSIDE"], None if ok else "yellow"))
        if not ok:
            warn(f"Model domain: the band brings two atoms to {dmin:.3f} Ang, below the shortest trained distance "
                 f"({dd['min']:.3f} Ang); the model is extrapolating there.")
    fm = dom.get("force_max_eV_per_ang")
    if fm:
        f0 = max(res["ends"][0]["f0"], res["ends"][1]["f0"])
        rows.append((["largest force (eV/Ang)", f"<= {fm['max']:.2f}", f"{f0:.3f} (endpoints as read)",
                      "inside" if f0 <= fm["max"] else "OUTSIDE"], None if f0 <= fm["max"] else "yellow"))
    er = dom.get("energy_per_atom_eV")
    if er:
        n = len(res["images"][0])
        epa = res["e_abs"] / n
        low = epa.min() < er["min"] - E_BELOW_TRAINED_TOL
        rows.append((["energy per atom (eV)", f"{er['min']:.4f} .. {er['max']:.4f}",
                      f"{epa.min():.4f} .. {epa.max():.4f}", "BELOW the lowest trained" if low else "inside"],
                     "yellow" if low else None))
        if low:
            warn(f"Model domain: the band's energy per atom goes {1000 * (er['min'] - epa.min()):.0f} meV below the "
                 "lowest training configuration -- a sign of a hole in the learned surface; check with SIESTA.")
    return rows


# ---- figures and files --------------------------------------------------------------------

def plot_band(series, card_dmin, fmax, path, siesta=None):
    """neb_profile.png: (a) energy profile, (b) NEB convergence, (c) shortest distance per image."""
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.4))
    for k, res in enumerate(series):
        c = _COLORS[k % 4]
        e0 = res["energies"][0]
        axs[0].plot(res["s_fit"], res["e_fit"] - e0, "-", color=c, lw=1.6)
        axs[0].plot(res["path"], res["energies"] - e0, "o", color=c, ms=6,
                    label=f"{res['label']}: {res['barrier']:.3f} eV")
        axs[0].plot([res["s_saddle"]], [res["barrier"]], "*", color=c, ms=13, mec="k", mew=0.5)
        steps = np.arange(1, len(res["history"]) + 1)
        r = np.array([h[2] for h in res["history"]])
        axs[1].semilogy(steps, r, "-", color=c, lw=1.4, label=res["label"])
        n1 = sum(1 for h in res["history"] if h[0] == 1)
        if n1 and n1 < len(r):
            axs[1].axvline(n1 + 0.5, color=c, ls=":", lw=1)
        idx = np.arange(len(res["dmin"]))
        axs[2].plot(idx, res["dmin"], "o-", color=c, label=res["label"])
    if siesta:
        axs[0].plot([p[0] for p in siesta], [p[1] for p in siesta], "kx", ms=10, mew=2,
                    label="SIESTA, same geometry")
    axs[0].axhline(0, color="0.6", lw=0.8)
    axs[0].set_xlabel("reaction coordinate (Ang)")
    axs[0].set_ylabel("E - E(initial) (eV)")
    axs[0].set_title("(a) energy profile (points: images, line: spline)")
    axs[0].legend(fontsize=8)
    axs[1].axhline(fmax, color="0.4", ls="--", lw=1, label=f"--fmax {fmax:g}")
    axs[1].set_xlabel("optimizer step (dotted: climbing starts)")
    axs[1].set_ylabel("NEB residual force (eV/Ang)")
    axs[1].set_title("(b) convergence")
    axs[1].legend(fontsize=8)
    if card_dmin:
        axs[2].axhline(card_dmin, color="0.4", ls="--", lw=1, label="shortest trained distance")
    axs[2].set_xlabel("image")
    axs[2].set_ylabel("shortest distance (Ang)")
    axs[2].set_title("(c) closest contact along the band")
    axs[2].legend(fontsize=8)
    for ax in axs:
        ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    return fig


def write_band_dat(path, res):
    with open(path, "w") as f:
        f.write(f"# stb-mlneb band, {res['label']}\n# image  s(Ang)  E-E0(eV)  max|F|(eV/Ang)  d_min(Ang)  "
                "max_disp_from_initial(Ang)\n")
        e0 = res["energies"][0]
        for i in range(len(res["images"])):
            f.write(f"{i:3d} {res['path'][i]:10.5f} {res['energies'][i] - e0:12.6f} {res['forces'][i]:10.5f} "
                    f"{res['dmin'][i]:9.4f} {res['disp'][i]:9.4f}\n")


def gnuplot_files(plot_dir, series, fmax):
    """neb_profile.gplot: the energy profile (images + spline) and the convergence of each model."""
    os.makedirs(plot_dir, exist_ok=True)
    for k, res in enumerate(series):
        write_band_dat(os.path.join(plot_dir, f"band_{k}.dat"), res)
        with open(os.path.join(plot_dir, f"spline_{k}.dat"), "w") as f:
            f.write("# s(Ang) E-E0(eV)\n")
            for s_, e_ in zip(res["s_fit"], res["e_fit"] - res["energies"][0]):
                f.write(f"{s_:.5f} {e_:.6f}\n")
        with open(os.path.join(plot_dir, f"convergence_{k}.dat"), "w") as f:
            f.write("# step stage residual(eV/Ang)\n")
            for j, (stage, _, r) in enumerate(res["history"], 1):
                f.write(f"{j} {stage} {r:.6g}\n")
    L = ["# --- STB Plot Configuration ---\n# Generated by stb-mlneb: gnuplot neb_profile.gplot\n",
         'set terminal pdfcairo enhanced color font "Arial,9" size 10,3.8\nset output "neb_profile.pdf"\n',
         "set grid\nset multiplot layout 1,2\n",
         'set title "Energy profile"\nset xlabel "reaction coordinate (Ang)"\nset ylabel "E - E(initial) (eV)"\n'
         "set key top left font \",8\"\n"]
    parts = []
    for k, res in enumerate(series):
        c = _COLORS[k % 4]
        parts.append(f'"spline_{k}.dat" using 1:2 with lines lw 2 lc rgb "{c}" notitle')
        parts.append(f'"band_{k}.dat" using 2:3 with points pt 7 ps 0.8 lc rgb "{c}" '
                     f'title "{res["label"]}: {res["barrier"]:.3f} eV"')
    L.append("plot " + ", \\\n     ".join(parts) + "\n")
    L.append('set title "NEB convergence"\nset xlabel "optimizer step"\nset ylabel "residual force (eV/Ang)"\n'
             "set logscale y\nset key top right font \",8\"\n")
    parts = [f'"convergence_{k}.dat" using 1:3 with lines lw 2 lc rgb "{_COLORS[k % 4]}" title "{res["label"]}"'
             for k, res in enumerate(series)]
    parts.append(f'{fmax} with lines dt 2 lc rgb "#666666" title "--fmax"')
    L.append("plot " + ", \\\n     ".join(parts) + "\nunset multiplot\n")
    with open(os.path.join(plot_dir, "neb_profile.gplot"), "w") as f:
        f.writelines(L)
    return ["neb_profile.gplot"]


# ---- report pieces --------------------------------------------------------------------------

def report_endpoints(res, f_out, warn):
    rows = []
    for e in res["ends"]:
        verdict = ("as given" if e["converged"] is None else "converged" if e["converged"] else "STEP CAP")
        rows.append(([e["label"], f"{e['f0']:.4f}", f"{e['f1']:.4f}", str(e["steps"]), verdict, f"{e['moved']:.4f}",
                      f"{e['e0'] - e['e1']:.4f}", f"{e['e1']:.5f}"], "yellow" if e["converged"] is False else None))
        if e["converged"] is False:
            warn(f"{res['label']}: the {e['label']} endpoint hit --max-steps without converging (max |F| "
                 f"{e['f1']:.3f} eV/Ang); the barrier carries its residual forces.")
        if e["converged"] is not None and e["moved"] > 0.3:
            warn(f"{res['label']}: the {e['label']} endpoint moved by up to {e['moved']:.2f} Ang during its relaxation "
                 "-- it was not a minimum of this model; check that it is still the configuration you meant.")
    print_dual(f"--- {res['label']}", f_out)
    print_table(["endpoint", "max|F| before", "max|F| after", "steps", "relaxation", "largest move Ang",
                 "energy gained eV", "E (eV)"], rows, f_out)
    if res["ends"][0]["converged"] is None:
        warn(f"{res['label']}: --no-pre-relax-endpoints -- the endpoints are not relaxed with this model, so the "
             "barrier is measured from a point that is not its minimum.")
    print_dual(f"Endpoint energies : E(final) - E(initial) = {res['ends'][1]['e1'] - res['ends'][0]['e1']:+.4f} eV",
               f_out)


def report_guess(res, args, f_out, warn):
    steps = res["guess_steps"]
    mean = float(np.mean(steps)) if steps else 0.0
    print_dual(f"--- {res['label']}: {len(res['guess'])} images, linear interpolation"
               + (f" + IDPP (--idpp-fmax {args.idpp_fmax}, --idpp-steps {args.idpp_steps})" if args.idpp else ""),
               f_out)
    if res.get("crossed"):
        print_dual(f"Periodic images   : {res['crossed']} atom(s) sit on opposite faces of the cell in the two endpoints; "
                   "each goes to the image nearest to its initial position (the shortest path)", f_out)
    rows = []
    for i in range(len(res["guess"])):
        st = f"{steps[i - 1]:.4f}" if i else "-"
        flag = ""
        if i and steps[i - 1] < 0.05:
            flag = "nearly identical to the previous"
        elif i and mean and steps[i - 1] > 3 * mean:
            flag = "large jump"
        rows.append(([str(i), st, f"{res['guess_dmin'][i]:.4f}", flag or "ok"], "yellow" if flag else None))
    print_table(["image", "step from previous (Ang)", "shortest distance (Ang)", "check"], rows, f_out)
    if any(s < 0.05 for s in steps):
        warn(f"{res['label']}: consecutive images of the starting path are nearly identical -- fewer --n-images, or "
             "--initial and --final are almost the same structure.")
    if mean and any(s > 3 * mean for s in steps):
        warn(f"{res['label']}: one step of the starting path is more than 3x the mean step -- often a wrong atom "
             "correspondence between --initial and --final (see --autosort-tol) rather than a real feature.")
    if res["frozen"]:
        print_dual(f"Frozen atoms      : {len(res['frozen'])}/{len(res['images'][0])} (move < {args.freeze_threshold} "
                   "Ang between the endpoints; --freeze-substrate)", f_out)


def report_band(res, args, f_out, warn):
    n1 = sum(1 for h in res["history"] if h[0] == 1)
    n2 = len(res["history"]) - n1
    print_dual(f"--- {res['label']}: {'converged' if res['converged'] else 'NOT converged (step cap)'}; "
               f"{n1} shaping + {n2} climbing step(s) of --max-steps {args.max_steps}; residual force "
               f"{res['residual']:.4f} eV/Ang (target --fmax {args.fmax}); {res['seconds']:.1f} s", f_out)
    if not res["converged"]:
        warn(f"{res['label']}: the NEB stopped at --max-steps with a residual force of {res['residual']:.3f} eV/Ang "
             f"(> --fmax {args.fmax}); the barrier is not converged -- raise --max-steps.")
    e0 = res["energies"][0]
    rows = []
    for i in range(len(res["images"])):
        role = ("initial" if i == 0 else "final" if i == len(res["images"]) - 1 else
                "climbing (saddle)" if i == res["i_top"] else "")
        rows.append(([str(i), f"{res['path'][i]:.4f}", f"{res['energies'][i] - e0:+.4f}", f"{res['forces'][i]:.4f}",
                      f"{res['dmin'][i]:.4f}", f"{res['disp'][i]:.4f}", role], "green" if role.startswith("climb")
                     else None))
    print_table(["image", "s (Ang)", "E - E0 (eV)", "max|F| (eV/Ang)", "d_min (Ang)", "max shift from initial",
                 "role"], rows, f_out)
    print_dual("max|F| is the raw force of the model; the NEB converges the force perpendicular to the band (and, on "
               "the climbing image, the full force), so a raw force above --fmax on an interior image is normal.",
               f_out)
    E = res["energies"]
    minima = [i for i in range(1, len(E) - 1) if E[i] < E[i - 1] and E[i] < E[i + 1]]
    if minima:
        warn(f"{res['label']}: the band has an intermediate minimum at image(s) {', '.join(map(str, minima))}: the "
             "transition goes through a metastable state. Relax that image and run two NEBs, one per step.")
    if res["i_top"] in (0, len(E) - 1) or res["barrier"] <= 1e-4:
        warn(f"{res['label']}: no barrier along this band (the highest point is an endpoint): the transition is "
             "downhill with this model, or the endpoints are not distinct minima.")


def report_barrier(res, f_out, warn):
    b, bi, dE = res["barrier"], res["barrier_img"], res["dE"]
    print_dual(f"--- {res['label']}", f_out)
    rows = [(["forward  (initial -> final)", f"{b:.4f}", f"{bi:.4f}"], None),
            (["backward (final -> initial)", f"{b - dE:.4f}", f"{bi - dE:.4f}"], None),
            (["reaction energy E(final) - E(initial)", f"{dE:+.4f}", f"{dE:+.4f}"], None)]
    print_table(["eV", "spline (reported)", "highest image"], rows, f_out)
    print_dual(f"Saddle            : s = {res['s_saddle']:.4f} of {res['path'][-1]:.4f} Ang "
               f"({100 * res['s_saddle'] / res['path'][-1]:.1f} % of the path); climbing image {res['i_top']}"
               if res["path"][-1] > 0 else "Saddle            : path of zero length", f_out)
    if abs(b - bi) > SPLINE_VS_IMAGE_TOL:
        warn(f"{res['label']}: the spline barrier ({b:.3f} eV) and the highest image ({bi:.3f} eV) differ by "
             f"{abs(b - bi):.3f} eV -- the band is too coarse near the saddle; more --n-images.")


def report_saddle(res, mv, f_out):
    if not mv:
        print_dual("No atom moves more than 0.1 Ang between the endpoints.", f_out)
        return
    ims = res["images"]
    i, s_ = mv[0][0], mv[0][1]
    rows = []
    for k, lab in ((0, "initial"), (res["i_top"], "saddle"), (len(ims) - 1, "final")):
        d = neighbor_distances(ims[k], i)
        rows.append(([lab, f"{k}", " ".join(f"{x:.3f}" for x in d)], "green" if lab == "saddle" else None))
    print_dual(f"Atom {i + 1} ({s_}, the largest mover, {mv[0][2]:.3f} Ang between the endpoints): four nearest "
               "neighbours along the band", f_out)
    print_table(["state", "image", "nearest distances (Ang)"], rows, f_out)
    print_dual(f"Saddle image      : max|F| {res['forces'][res['i_top']]:.4f} eV/Ang, shortest distance "
               f"{res['dmin'][res['i_top']]:.4f} Ang, largest shift from the initial image "
               f"{res['disp'][res['i_top']]:.3f} Ang", f_out)


def report_comparison(series, f_out):
    hdr = ["", *[r["label"] for r in series]]
    keys = (("barrier, spline (eV)", lambda r: f"{r['barrier']:.4f}"),
            ("barrier, highest image (eV)", lambda r: f"{r['barrier_img']:.4f}"),
            ("backward barrier (eV)", lambda r: f"{r['barrier'] - r['dE']:.4f}"),
            ("reaction energy (eV)", lambda r: f"{r['dE']:+.4f}"),
            ("saddle at s (Ang)", lambda r: f"{r['s_saddle']:.4f}"),
            ("path length (Ang)", lambda r: f"{r['path'][-1]:.4f}"),
            ("shortest distance on band (Ang)", lambda r: f"{min(r['dmin']):.4f}"),
            ("NEB steps", lambda r: f"{len(r['history'])} ({'converged' if r['converged'] else 'step cap'})"),
            ("time (s)", lambda r: f"{r['seconds']:.1f}"))
    print_table(hdr, [([k, *[f(r) for r in series]], None) for k, f in keys], f_out)
    b0, b1 = series[0]["barrier"], series[1]["barrier"]
    print_dual(f"Difference        : the foundation model's barrier is {b1 - b0:+.4f} eV "
               f"({100 * (b1 - b0) / b0:+.1f} %) from the fine-tuned one" if abs(b0) > 1e-6 else "", f_out)


# ---- committee uncertainty -------------------------------------------------------------------

def committee_members(card, card_path, workflow_path, explicit):
    """(members, note): the committee as [(model, card)] -- --committee-models (each with the
    mlff_model_card.json next to it), else the card's own committee section, whose paths are relative
    to the workflow run folder: --mlff-workflow's folder, or the folder above the card's (best_model/
    lives in the run folder). `note` says where they came from, or why there are none."""
    if explicit:
        out = []
        for m in explicit:
            c = os.path.join(os.path.dirname(os.path.abspath(m)), model_card.CARD_FILE)
            if not os.path.isfile(m) or not os.path.isfile(c):
                return [], f"--committee-models: '{m}' or its {model_card.CARD_FILE} not found"
            out.append((m, c))
        return out, f"{len(out)} model(s) from --committee-models"
    members = ((card or {}).get("committee") or {}).get("members") or []
    if len(members) < 2:
        return [], "the model card records no committee (trained with --committee 1)"
    bases = []
    if workflow_path:
        bases.append(os.path.dirname(os.path.abspath(workflow_path)))
    if card_path:
        bases += [os.path.dirname(os.path.dirname(os.path.abspath(card_path))),
                  os.path.dirname(os.path.abspath(card_path))]
    for b in bases:
        found = [(os.path.join(b, m["model"]), os.path.join(b, m["card"])) for m in members
                 if m.get("model") and m.get("card")]
        if found and all(os.path.isfile(x) and os.path.isfile(y) for x, y in found):
            return found, f"{len(found)} members of round {(card.get('committee') or {}).get('round')} from {b}"
    return [], ("the card's committee members (" + ", ".join(m.get("model", "?") for m in members) + ") were not "
                "found next to it -- pass --mlff-workflow (the run folder) or --committee-models")


def committee_band(res, members, device, noise):
    """Every member on every relaxed image of `res`: energies relative to image 0, the force spread
    per atom (sqrt of the mean squared deviation from the committee mean, stb-mlffValidate's
    definition) and each member's barrier at the climbing image."""
    from stb.core import mace_relax
    e_rel, forces = [], []
    for model, card in members:
        model_card.register(model, card, echo=False)
        model_card.mark_announced(model)
        with _quiet(noise, "committee"):
            calc = mace_relax.get_calculator(model=model, device=device)
            e, f = [], []
            for a in res["images"]:
                b = a.copy()
                b.calc = calc
                e.append(b.get_potential_energy())
                f.append(b.get_forces())
        e = np.array(e)
        e_rel.append(e - e[0])
        forces.append(np.array(f))
    e_rel, forces = np.array(e_rel), np.array(forces)          # (m, n_img), (m, n_img, N, 3)
    mean_f = forces.mean(axis=0)
    spread_atom = np.sqrt(np.mean(np.sum((forces - mean_f) ** 2, axis=3), axis=0))   # (n_img, N)
    i_top = res["i_top"]
    return {"n": len(members), "e_rel": e_rel, "e_std": e_rel.std(axis=0),
            "f_spread": spread_atom.max(axis=1), "f_spread_atom": spread_atom.argmax(axis=1),
            "barriers": e_rel[:, i_top] - e_rel[:, 0], "i_top": i_top}


def report_committee(com, res, card, f_out, warn):
    cal = ((card or {}).get("validation") or {}).get("calibration") or {}
    k = cal.get("error_over_spread")
    ref = cal.get("spread_mean_meV")
    rows = []
    for i in range(len(res["images"])):
        sp = 1000 * com["f_spread"][i]
        cells = [str(i), f"{res['energies'][i] - res['energies'][0]:+.4f}", f"{com['e_std'][i]:.4f}", f"{sp:.1f}",
                 str(com["f_spread_atom"][i] + 1), f"{k * sp:.0f}" if k else "-", f"{sp / ref:.1f}x" if ref else "-"]
        rows.append((cells, "green" if i == com["i_top"] else None))
    print_dual(f"{com['n']} committee members on the {res['label']} band (each member's energies relative to its own "
               "image 0):", f_out)
    print_table(["image", "E - E0 (eV)", "E spread (eV)", "F spread (meV/Ang)", "atom", "est. F error (meV/Ang)",
                 "vs test set"], rows, f_out)
    b = com["barriers"]
    print_dual(f"Barrier per member: {', '.join(f'{x:.3f}' for x in b)} eV (spread {b.std():.3f} eV, range "
               f"{b.max() - b.min():.3f} eV) at the climbing image", f_out)
    if k and ref:
        print_dual(f"Calibration       : on the card's test set the true force error is ~{k:.1f} x the spread, and the "
                   f"mean spread there is {ref:.1f} meV/Ang ('vs test set' = this image's spread over that mean)", f_out)
    else:
        print_dual("Calibration       : the card has no committee calibration (stb-mlffValidate fills it): the spread is "
                   "only a lower bound of the error", f_out)
    top = 1000 * com["f_spread"].max()
    if ref and top > 3 * ref:
        warn(f"Committee: the force spread reaches {top:.0f} meV/Ang on the band ({top / ref:.1f}x the test set's mean): "
             "the band leaves what the training set constrains -- check with SIESTA (--siesta-check) and add these "
             "images to training (--to-mlff).")
    if b.std() > 0.05:
        warn(f"Committee: the members' barriers differ by {b.std():.3f} eV (std): the barrier is not determined by the "
             "training data.")
    print_dual("Seeds fine-tuned from the same foundation model are nearly the same model: the spread is a LOWER "
               "bound of the error (7.1 benchmark: it tracked the error per configuration, 3-15x smaller).", f_out)


# ---- SIESTA check and export to the MLFF workflow ----------------------------------------------

SIESTA_DIR = "siesta_check"
IMAGE_STAMP = "stb_mlneb_image.xyz"
SINGLE_POINT_EXTRA = ("# stb-mlneb: single-point evaluation of a NEB image (no relaxation)\n"
                      "MD.TypeOfRun          CG\nMD.Steps              0\nMD.NumCGsteps         0\n"
                      "MD.VariableCell       false\n")


def siesta_template(args, state, wf_base):
    """The SIESTA inputs every check/export folder starts from. With --mlff-workflow: the workflow's
    reference calc.fdf, pseudopotentials and dispersion, the k-grid scaled to keep the reference's
    k x a; else --calc and --pseudo-dir, the calc's own k-grid (or --kgrid). Raises ValueError."""
    from stb.core.fdf_directives import parse_kgrid
    from stb.core.pseudopotentials import BANKS, resolve_pseudo_source
    from stb.core import mlff_workflow
    ref_cell, extra = None, ""
    if state is not None:
        rec = state.get("reference") or {}
        if not rec.get("calc"):
            raise ValueError("the workflow's reference is one self-contained .fdf (no separate calc.fdf): pass --calc")
        calc = mlff_workflow.resolve(rec["calc"], wf_base)
        ref = structure_io.read_fdf(mlff_workflow.resolve(rec["file"], wf_base))
        ref_cell = np.linalg.norm(np.array(structure_io.lattice_only(ref)), axis=1)
        pseudo_dir = rec.get("pseudo_dir") or "."
        if pseudo_dir not in BANKS:
            pseudo_dir = mlff_workflow.resolve(pseudo_dir, wf_base)
        extra = mlff_workflow.dftd3_block(state)
        source = f"the workflow's reference ({rec['calc']}, pseudopotentials {rec.get('pseudo_dir')})"
    else:
        if not args.calc:
            raise ValueError("--siesta-check needs --calc (or --mlff-workflow)")
        calc, pseudo_dir = args.calc, args.pseudo_dir
        source = f"--calc {args.calc}, pseudopotentials {args.pseudo_dir}"
    if args.calc and state is not None:
        calc = args.calc
        source = f"--calc {args.calc} (k-grid scaled from the workflow's reference cell)"
    if not os.path.isfile(calc):
        raise ValueError(f"calc template '{calc}' not found")
    with open(calc) as f:
        text = f.read()
    inc = [m for m in re.findall(r"^\s*%include\s+(\S+)", text, re.M | re.I) if m != "config_extra.fdf"]
    if not inc:
        raise ValueError(f"'{calc}' does not %include a structure file: the images need a calc.fdf without geometry "
                         "that %includes one (the suite's structure.fdf + calc.fdf pair)")
    return {"calc": calc, "struct": inc[0], "pseudo": resolve_pseudo_source(pseudo_dir), "pseudo_label": pseudo_dir,
            "kgrid": parse_kgrid(text), "ref_cell": ref_cell, "extra": extra, "source": source}


def image_kgrid(tpl, atoms, override):
    if override:
        return list(override), "--kgrid"
    if tpl["kgrid"] is None:
        return None, "none in the template (SIESTA samples Gamma)"
    if tpl["ref_cell"] is None:
        return list(tpl["kgrid"]), "the template's own"
    lengths = atoms.cell.lengths()
    kg = [max(1, int(round(k * lr / l))) for k, lr, l in zip(tpl["kgrid"], tpl["ref_cell"], lengths)]
    return kg, f"the reference's {'x'.join(map(str, tpl['kgrid']))} scaled to this cell (same k x a)"


def write_siesta_folder(folder, atoms, tpl, species_meta, kgrid_override=None):
    """One single-point SIESTA folder for `atoms`: the structure file the template %includes, the
    calc.fdf with config_extra.fdf on top (k-grid, NumberOfAtoms, dispersion, single point), the
    pseudopotentials, and the image's coordinates as a stamp to recognise it later."""
    from pymatgen.io.ase import AseAtomsAdaptor
    from ase.io import write as ase_write
    from stb.core.pseudopotentials import get_required_pseudos
    os.makedirs(folder, exist_ok=True)
    st = structure_io.from_pymatgen(AseAtomsAdaptor.get_structure(atoms), species_meta=species_meta,
                                    coord_format="fractional")
    structure_io.write_fdf(st, os.path.join(folder, tpl["struct"]))
    kg, _ = image_kgrid(tpl, atoms, kgrid_override)
    extra = "# Written by stb-mlneb -- %included first, so it wins over calc.fdf\n"
    if kg:
        extra += f"kgrid.MonkhorstPack [{kg[0]} {kg[1]} {kg[2]}]\n"
    extra += f"NumberOfAtoms {len(atoms)}\n\n" + tpl["extra"] + SINGLE_POINT_EXTRA
    with open(os.path.join(folder, "config_extra.fdf"), "w") as f:
        f.write(extra)
    with open(tpl["calc"]) as f:
        text = f.read()
    with open(os.path.join(folder, os.path.basename(tpl["calc"])), "w") as f:
        f.write(structure_io.prepend_include(text, "config_extra.fdf"))
    pseudos, missing = get_required_pseudos(sorted(set(atoms.get_chemical_symbols())), tpl["pseudo"])
    if missing:
        raise ValueError(f"no pseudopotential for {', '.join(missing)} in '{tpl['pseudo_label']}'")
    for p in pseudos:
        shutil.copy(p, folder)
    ase_write(os.path.join(folder, IMAGE_STAMP), atoms.copy(), format="extxyz")
    return folder


def read_siesta_folder(folder, atoms=None, tol=0.02):
    """'missing' | 'stale' (written for another geometry) | 'pending' (not finished) | dict with the
    SIESTA energy and largest force of a finished folder."""
    from ase.io import read as ase_read
    from stb.core import siesta_log
    if not os.path.isdir(folder):
        return "missing"
    stamp = os.path.join(folder, IMAGE_STAMP)
    if atoms is not None and os.path.isfile(stamp):
        old = ase_read(stamp)
        if len(old) != len(atoms) or np.abs(old.get_positions() - atoms.get_positions()).max() > tol:
            return "stale"
    label = structure_io.find_fdf_system_label(folder)
    out = siesta_log.find_out_file(folder, label) if label else None
    if not out or not os.path.isfile(out):
        return "pending"
    with open(out, errors="replace") as f:
        if "End of run" not in f.read():
            return "pending"
    e = siesta_log.get_free_energy(out)
    fa = siesta_log.read_fa_forces(os.path.join(folder, f"{label}.FA"))
    if e is None:
        return "pending"
    return {"E": e, "fmax": float(np.linalg.norm(fa, axis=1).max()) if fa is not None else None, "out": out}


def check_images(res, which):
    n = len(res["images"])
    return list(range(n)) if which == "all" else sorted({0, res["i_top"], n - 1})


def report_siesta(res, results, f_out, warn, args):
    """The model's band against SIESTA single points on the same geometries."""
    done = {i: r for i, r in results.items() if isinstance(r, dict)}
    rows = []
    for i, r in sorted(results.items()):
        model = res["energies"][i] - res["energies"][0]
        if isinstance(r, dict) and 0 in done:
            s_ = r["E"] - done[0]["E"]
            rows.append(([str(i), f"{model:+.4f}", f"{s_:+.4f}", f"{model - s_:+.4f}",
                          f"{r['fmax']:.3f}" if r["fmax"] is not None else "-", "done"],
                         "green" if i == res["i_top"] else None))
        else:
            rows.append(([str(i), f"{model:+.4f}", "-", "-", "-", r if isinstance(r, str) else "done"], "yellow"))
    print_dual(f"SIESTA single points on the {res['label']} band (E - E of image 0; SIESTA max|F| tells how far the "
               "geometry is from a SIESTA stationary point):", f_out)
    print_table(["image", "model (eV)", "SIESTA (eV)", "model - SIESTA", "SIESTA max|F|", "status"], rows, f_out)
    n = len(res["images"]) - 1
    if all(i in done for i in (0, res["i_top"], n)):
        sb = done[res["i_top"]]["E"] - done[0]["E"]
        sbb = done[res["i_top"]]["E"] - done[n]["E"]
        mb, mbb = res["barrier_img"], res["barrier_img"] - res["dE"]
        print_dual(f"Barrier on this path: model {mb:.4f} eV forward / {mbb:.4f} backward; SIESTA {sb:.4f} / "
                   f"{sbb:.4f} eV ({100 * (mb - sb) / sb:+.0f} % model vs SIESTA)" if abs(sb) > 1e-6 else
                   f"Barrier on this path: model {mb:.4f} eV, SIESTA {sb:.4f} eV", f_out)
        print_dual("SIESTA at the model's geometries is not SIESTA's own barrier: its endpoints and saddle are not "
                   "relaxed with SIESTA (see max|F|); the difference to the model is still the model's error there.",
                   f_out)
        dev = abs(mb - sb)
        if dev > max(0.05, 0.5 * abs(sb)):
            warn(f"SIESTA check: on its own path the model's barrier is {mb:.3f} eV, SIESTA's {sb:.3f} eV -- the model is "
                 "not reliable for this transition. Add the images to the training set (--to-mlff), retrain "
                 "(stb-mlffAnalysis) and run the NEB again.")
        elif dev > max(0.05, 0.2 * abs(sb)):
            warn(f"SIESTA check: the model's barrier ({mb:.3f} eV) is {dev:.3f} eV ({100 * dev / abs(sb):.0f} %) from "
                 f"SIESTA's ({sb:.3f} eV) on the same path: right in kind, not yet in value -- one more round with these "
                 "images (--to-mlff) would close it, or use SIESTA's number.")
    pending = [i for i, r in results.items() if r == "pending"]
    if pending:
        warn(f"SIESTA check: {len(pending)} folder(s) still need SIESTA ({args.output_dir}/{SIESTA_DIR}/): run it there, "
             "then run stb-mlneb again with the same options to compare.")


def export_to_mlff(args, state, wf_path, res, tpl, species_meta, results, warn, f_out):
    """Writes the band's images (and --to-mlff-rattle copies) as training configs of the next round of
    the MLFF workflow: <run>/configs_round_NN/neb_<tag>/image_II[_rattle_K]/, recorded in the JSON
    (set train, kind neb / neb-rattle). An image already labelled by --siesta-check (finished, same
    geometry) is copied with its SIESTA output, so it is not computed twice. Returns the folders that
    still need SIESTA."""
    from stb.core import mlff_workflow
    base = mlff_workflow.base_dir(wf_path)
    latest = max([r.get("round", 0) for r in state.get("rounds", [])] or [0])
    tag = re.sub(r"[^A-Za-z0-9_.-]+", "_", os.path.basename(os.path.abspath(args.output_dir))) or "neb"
    root = os.path.join(mlff_workflow.configs_dir(base, latest + 1), tag if tag.startswith("neb") else f"neb_{tag}")
    known = {c["folder"] for c in state["configs"]}
    entries, pending, copied, skipped = [], [], 0, 0
    rng = np.random.default_rng(args.to_mlff_seed)
    for i, image in enumerate(res["images"]):
        jobs = [(f"image_{i:02d}", image, "neb", f"NEB image {i} of {len(res['images']) - 1}, {res['label']} band")]
        for c in range(args.to_mlff_rattle):
            a = image.copy()
            a.set_positions(a.get_positions() + rng.normal(0.0, args.to_mlff_stdev, (len(a), 3)))
            a.wrap()
            jobs.append((f"image_{i:02d}_rattle_{c + 1}", a, "neb-rattle",
                         f"NEB image {i} rattled {args.to_mlff_stdev:g} Ang"))
        for name, atoms, kind, origin in jobs:
            dest = os.path.join(root, name)
            if mlff_workflow.rel(dest, base) in known:
                skipped += 1
                continue
            r = results.get(i) if kind == "neb" else None
            if isinstance(r, dict):
                shutil.copytree(os.path.dirname(r["out"]), dest, dirs_exist_ok=True)
                copied += 1
            else:
                write_siesta_folder(dest, atoms, tpl, species_meta, args.kgrid)
                pending.append(dest)
            entries.append({"folder": dest, "set": "train", "kind": kind, "origin": origin, "stage": "neb",
                            "round": latest, "n_atoms": len(atoms),
                            "neb": {"image": i, "initial": args.initial, "final": args.final,
                                    "model": res["label"], "E_rel_model_eV": round(float(
                                        res["energies"][i] - res["energies"][0]), 5)}})
    if entries:
        mlff_workflow.backup(wf_path)
        mlff_workflow.add_configs(state, entries, base)
        mlff_workflow.record(state, "neb", "stb-mlneb", sys.argv[1:])
        mlff_workflow.save(state, wf_path)
    print_dual(f"Workflow          : {wf_path} (latest round {latest}; these configs train round {latest + 1})", f_out)
    print_dual(f"Written           : {len(entries)} training config(s) in {os.path.relpath(root)}/ -- {copied} already "
               f"labelled by --siesta-check (copied with their output), {len(pending)} need SIESTA"
               + (f"; {skipped} already in the workflow, left as they are" if skipped else ""), f_out)
    print_dual(f"k-grid            : {image_kgrid(tpl, res['images'][0], args.kgrid)[1]}; template: {tpl['source']}", f_out)
    if pending:
        print_dual(f"Next              : run SIESTA in each folder that needs it, e.g.\n  for d in {os.path.relpath(root)}/*/; "
                   "do (cd $d && mpirun -np 4 siesta < calc.fdf > siesta.out); done", f_out)
    print_dual(f"Then              : stb-mlffAnalysis --workflow {wf_path} --init-from previous  (round {latest + 1}), "
               "stb-mlffValidate, and this NEB again with the new model", f_out)
    if pending:
        warn(f"--to-mlff: {len(pending)} new training folder(s) need SIESTA before stb-mlffAnalysis can use them "
             f"({os.path.relpath(root)}/).")
    return pending


# ---- command line ---------------------------------------------------------------------------

def build_parser():
    parser = argparse.ArgumentParser(
        description="Climbing-image NEB (nudged elastic band) driven entirely by a MACE potential -- the MACE-MP-0 "
                    "foundation model, or a custom model fine-tuned on your own SIESTA data via stb-mlffAnalysis "
                    "(--custom-model, always with its model card). Relaxes the two endpoint structures (--initial/"
                    "--final, same composition and cell), interpolates the band (linear + IDPP), relaxes it with a "
                    "two-stage climbing-image NEB and reports: the atoms that move, the endpoint relaxations, the "
                    "starting path, the NEB convergence, the band image by image, the barrier forward and backward "
                    "(spline and highest image), the saddle geometry, the band against the model card's trained "
                    "range and, with a custom model, the same NEB with the foundation model and the uncertainty of its "
                    "committee; every warning in one final section. Closing the loop with SIESTA: --siesta-check writes "
                    "single-point folders for the endpoints and the saddle (then reads them on the next run and compares "
                    "the model with SIESTA on the same geometries), and --to-mlff adds the images to the training set of "
                    "the 7.1 workflow (--mlff-workflow) for the next round. For the full DFT barrier use stb-neb (5.9).",
        epilog="Example usage:\n"
               "  stb-mlneb --initial initial.fdf --final final.fdf\n"
               "  stb-mlneb --initial initial.fdf --final final.fdf --custom-model mlff_model.model "
               "--custom-model-card mlff_model_card.json --device cuda --n-images 9 --show --save-gnuplot\n"
               "  stb-mlneb --initial initial.fdf --final final.fdf --freeze-substrate --save-images --save-report\n"
               "  stb-mlneb --initial i.fdf --final f.fdf --custom-model best_model/mlff_model_compiled.model "
               "--custom-model-card best_model/mlff_model_card.json --mlff-workflow mlff_runs/mlff_workflow.json "
               "--siesta-check            # then run SIESTA, run it again, and add --to-mlff\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--initial", required=True, help="Initial (relaxed) endpoint structure (.fdf).")
    parser.add_argument("--final", required=True, help="Final (relaxed) endpoint structure (.fdf).")
    parser.add_argument("--model", choices=["small", "medium", "large"], default="medium",
                        help="MACE-MP-0 model size (default: medium). With --custom-model, the foundation model of "
                             "the comparison.")
    parser.add_argument("--custom-model", default=None, metavar="PATH",
                        help="Use a custom MACE model file instead of the MACE-MP-0 foundation potential -- e.g. one "
                             "fine-tuned via stb-mlffAnalysis. Overrides --model.")
    model_card.add_card_argument(parser, "--custom-model")
    parser.add_argument("--skip-foundation-comparison", action="store_true",
                        help="With --custom-model, do not repeat the NEB with the foundation model (--model). The "
                             "comparison is on by default: it shows what the fine-tuning changed in the barrier.")
    parser.add_argument("--mlff-workflow", default=None, metavar="JSON",
                        help="The mlff_workflow.json of the 7.1 run that trained --custom-model: locates its committee "
                             "(for the uncertainty), and gives --siesta-check/--to-mlff the training level of theory "
                             "(calc.fdf, pseudopotentials, dispersion, k x a).")
    parser.add_argument("--no-committee", dest="committee", action="store_false",
                        help="Do not evaluate the committee of the custom model on the band. On by default when the card "
                             "records a committee whose model files are found (next to the card's run folder, or in "
                             "--mlff-workflow's): the force spread per image and each member's barrier.")
    parser.add_argument("--committee-models", nargs="+", default=None, metavar="PATH",
                        help="The committee's model files explicitly (each with its mlff_model_card.json next to it), "
                             "instead of the card's own list.")
    parser.add_argument("--siesta-check", action="store_true",
                        help=f"Write single-point SIESTA folders for images of the relaxed band into <output-dir>/"
                             f"{SIESTA_DIR}/image_NN/ (--siesta-check-images). Run SIESTA there and run stb-mlneb again "
                             "with the same options: finished folders whose geometry matches the band are read and "
                             "the model's profile is compared with SIESTA's on the same geometries. Needs --calc or "
                             "--mlff-workflow.")
    parser.add_argument("--siesta-check-images", choices=["ends-saddle", "all"], default="ends-saddle",
                        help="Which images --siesta-check writes: the two endpoints and the climbing image (default), "
                             "or every image.")
    parser.add_argument("--calc", default=None,
                        help="calc.fdf template for --siesta-check/--to-mlff (no geometry; it must %%include the "
                             "structure file, which each folder gets). Default with --mlff-workflow: its reference "
                             "calc.fdf.")
    parser.add_argument("--pseudo-dir", default=".",
                        help="Pseudopotentials for --siesta-check with --calc: a folder or a bank name (default: .).")
    parser.add_argument("--kgrid", type=int, nargs=3, default=None, metavar=("KA", "KB", "KC"),
                        help="k-grid of the check/export folders. Default: with --mlff-workflow the reference's grid "
                             "scaled to this cell (same k x a as the training set); with --calc its own grid.")
    parser.add_argument("--to-mlff", action="store_true",
                        help="Add the band's images to the training set of the --mlff-workflow run, as configs of its "
                             "next round (<run>/configs_round_NN/neb_<output-dir>/), with the workflow's level of "
                             "theory. Images already labelled by --siesta-check are copied with their output. Then "
                             "stb-mlffAnalysis trains the next round on them.")
    parser.add_argument("--to-mlff-rattle", type=int, default=0, metavar="N",
                        help="With --to-mlff, also N rattled copies of every image (default: 0): the neighbourhood of "
                             "the path, not only the path itself.")
    parser.add_argument("--to-mlff-stdev", type=float, default=0.05,
                        help="Gaussian displacement of --to-mlff-rattle copies, Ang (default: 0.05).")
    parser.add_argument("--to-mlff-seed", type=int, default=7,
                        help="Random seed of --to-mlff-rattle (default: 7).")
    parser.add_argument("--n-images", type=int, default=7,
                        help="Total images including both endpoints (default: 7, at least 3).")
    parser.add_argument("--autosort-tol", type=float, default=0.5,
                        help="Atom-correspondence tolerance between --initial/--final, Ang (default: 0.5, pymatgen's "
                             "own documented value).")
    parser.add_argument("--no-idpp", dest="idpp", action="store_false",
                        help="Skip the IDPP (Image Dependent Pair Potential) refinement of the interior images before "
                             "the NEB. On by default: a geometric smoothing that usually gives a much better start "
                             "than linear interpolation.")
    parser.add_argument("--idpp-fmax", type=float, default=0.1, help="Force convergence for IDPP (default: 0.1).")
    parser.add_argument("--idpp-steps", type=int, default=100, help="Max steps for IDPP (default: 100).")
    parser.add_argument("--no-pre-relax-endpoints", dest="pre_relax_endpoints", action="store_false",
                        help="Skip the relaxation (positions only) of both endpoints with the model. On by default, "
                             "so the barrier is measured from the model's own minima.")
    parser.add_argument("--freeze-substrate", action="store_true",
                        help="Freeze the atoms that move less than --freeze-threshold between --initial and --final "
                             "during the NEB: fewer degrees of freedom, no spectator atom wandering. Off by default.")
    parser.add_argument("--freeze-threshold", type=float, default=0.3,
                        help="Displacement threshold for --freeze-substrate, Ang (default: 0.3).")
    parser.add_argument("--k", type=float, default=0.1, help="NEB spring constant, eV/Ang^2 (default: 0.1).")
    parser.add_argument("--fmax", type=float, default=0.05,
                        help="Force convergence of the endpoint relaxations and of the NEB, eV/Ang (default: 0.05). "
                             "The shaping stage stops at max(2 x fmax, 0.1); the climbing stage at fmax.")
    parser.add_argument("--max-steps", type=int, default=200,
                        help="Optimizer steps of each endpoint relaxation, and of the NEB in total: the shaping stage "
                             "gets half, the climbing stage the rest (default: 200).")
    parser.add_argument("--device", choices=["cpu", "cuda"], default="cpu",
                        help="Device to run the model on (default: cpu).")
    parser.add_argument("-o", "--output-dir", default="mlneb_out", help="Output directory (default: mlneb_out).")
    parser.add_argument("--show", action="store_true",
                        help="Open the figure on screen at the end (needs a display; without one it is a warning). "
                             "neb_profile.png is written either way.")
    parser.add_argument("--save-gnuplot", action="store_true",
                        help=f"Write .dat + .gplot into <output-dir>/{PLOT_DIR}/ and run gnuplot there (neb_profile.pdf) "
                             "when it is installed. Off by default.")
    parser.add_argument("--save-data", action="store_true",
                        help="Write the band as a table: neb_band.dat (and neb_band_foundation.dat with the "
                             "comparison): image, reaction coordinate, energy, max force, shortest distance, shift. "
                             "Off by default.")
    parser.add_argument("--save-images", action="store_true",
                        help="Also write each relaxed image as image_NN.fdf (bare structure, no calc or "
                             "pseudopotentials) -- e.g. to start a SIESTA NEB (stb-neb) from this band. Off by default.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-mlneb {VERSION}")
    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()
    card = model_card.setup_from_args(args.custom_model, args.custom_model_card, "--custom-model", error=parser.error)
    if args.n_images < 3:
        parser.error("--n-images must be at least 3 (two endpoints and one image between them).")
    if args.to_mlff and not args.mlff_workflow:
        parser.error("--to-mlff needs --mlff-workflow (the run whose training set gets the images).")
    if args.mlff_workflow and not os.path.isfile(args.mlff_workflow):
        parser.error(f"--mlff-workflow file '{args.mlff_workflow}' not found.")

    import matplotlib
    if not args.show:
        matplotlib.use("Agg")

    if args.intro:
        show_intro(["Siesta ToolBox Suite - ML NEB",
                    "Climbing-image nudged elastic band, driven entirely by a MACE potential",
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
    if importlib.util.find_spec("mace") is None or importlib.util.find_spec("torch") is None:
        require_mace()
    with capture_library_noise(noise, "import"):
        require_mace()
        from pymatgen.io.ase import AseAtomsAdaptor
        from stb.core import mace_relax
        from stb.neb import check_composition_match, wrap_into_cell
        import torch
        import ase
        import mace
    model_card.collect_warnings(warns)

    print("\n" + color_text("ML NEB (CLIMBING-IMAGE):", 'bold'))
    print("-" * 60)
    print_dual(color_text("===== STB-MLNEB REPORT =====", 'magenta'), f_out)

    compare = bool(args.custom_model) and not args.skip_foundation_comparison
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
    print_dual(f"Endpoints         : {args.initial} -> {args.final}", f_out)
    print_dual(f"Band              : {args.n_images} images (endpoints included), linear interpolation"
               f"{' + IDPP' if args.idpp else ''}, atoms matched within --autosort-tol {args.autosort_tol} Ang", f_out)
    print_dual(f"NEB               : improved tangent, spring k {args.k} eV/Ang^2, FIRE; shaping to "
               f"{max(2 * args.fmax, 0.1):g} then climbing image to --fmax {args.fmax} eV/Ang, --max-steps "
               f"{args.max_steps} in total" + ("; substrate frozen below "
                                                f"{args.freeze_threshold} Ang" if args.freeze_substrate else ""), f_out)
    print_dual(f"Endpoints relaxed : {'yes, positions only, to --fmax' if args.pre_relax_endpoints else 'no (--no-pre-relax-endpoints)'}",
               f_out)
    print_dual(f"Model             : "
               f"{'custom (' + args.custom_model + ')' if args.custom_model else f'MACE-MP-0 ({args.model})'}", f_out)
    if compare:
        print_dual(f"Comparison        : the same NEB with the foundation model (MACE-MP-0 {args.model})", f_out)
    if report_path:
        print_dual(f"Report file       : {report_path}", f_out)

    sections("STRUCTURES")
    for path, flag in ((args.initial, "--initial"), (args.final, "--final")):
        if not os.path.isfile(path):
            fail(f"{flag} file '{path}' not found.")
    if args.custom_model and not os.path.isfile(args.custom_model):
        fail(f"--custom-model file not found: {args.custom_model}")
    try:
        initial_structure = structure_io.read_fdf(args.initial)
        final_structure = structure_io.read_fdf(args.final)
    except Exception as e:
        fail(f"Could not read structure(s): {e}")
    try:
        check_composition_match(initial_structure, final_structure)
    except SystemExit:
        fail("--initial and --final have a different composition: a NEB needs the same atoms at both ends.")
    species_meta = structure_io.species_dict(initial_structure)

    # The workflow and the SIESTA template are checked here, before any model is run: a missing calc.fdf, a template
    # without a structure %include or a missing pseudopotential must not cost a whole band first.
    state, tpl = None, None
    if args.mlff_workflow:
        from stb.core import mlff_workflow
        try:
            state = mlff_workflow.load(args.mlff_workflow)
        except ValueError as e:
            fail(f"--mlff-workflow: {e}")
    if args.siesta_check or args.to_mlff:
        try:
            tpl = siesta_template(args, state, os.path.dirname(os.path.abspath(args.mlff_workflow))
                                  if args.mlff_workflow else None)
        except ValueError as e:
            fail(f"--siesta-check/--to-mlff: {e}")
        from stb.core.pseudopotentials import get_required_pseudos
        elements = sorted({el.symbol for el in structure_io.to_pymatgen(initial_structure).composition.elements})
        _, missing = get_required_pseudos(elements, tpl["pseudo"])
        if missing:
            fail(f"--siesta-check/--to-mlff: no pseudopotential for {', '.join(missing)} in '{tpl['pseudo_label']}'")
    initial_pmg = wrap_into_cell(structure_io.to_pymatgen(initial_structure))
    final_pmg = resolve_lattice_mismatch(initial_pmg, wrap_into_cell(structure_io.to_pymatgen(final_structure)),
                                         warn=warn)
    a_i, a_f = AseAtomsAdaptor.get_atoms(initial_pmg), AseAtomsAdaptor.get_atoms(final_pmg)
    cp = a_i.cell.cellpar()
    vol = a_i.get_volume()
    print_dual(f"[OK] Same composition: {len(a_i)} atom(s), {a_i.get_chemical_formula(mode='hill')}", f_out)
    print_dual(f"Cell              : a {cp[0]:.5f}, b {cp[1]:.5f}, c {cp[2]:.5f} Ang; alpha {cp[3]:.2f}, beta "
               f"{cp[4]:.2f}, gamma {cp[5]:.2f} deg; {vol:.3f} Ang^3 ({vol / len(a_i):.4f} Ang^3/atom)", f_out)
    print_dual(f"Shortest distance : initial {model_card.shortest_distance(a_i):.4f}, final "
               f"{model_card.shortest_distance(a_f):.4f} Ang", f_out)
    mv, total = movers(a_i, a_f)
    print_dual(f"Atoms that move   : {len(mv)} by more than {MOVER_THRESHOLD} Ang between the endpoints (as read; "
               f"sum of all displacements {total:.3f} Ang)", f_out)
    if mv:
        print_table(["atom (fdf order)", "species", "displacement (Ang)"],
                    [([str(i + 1), s_, f"{d:.4f}"], None) for i, s_, d in mv[:8]], f_out)
        if len(mv) > 8:
            print_dual(f"  ... and {len(mv) - 8} more", f_out)
    else:
        warn("No atom moves more than 0.1 Ang between --initial and --final: there is no transition to compute.")
    if mv and mv[0][2] > 4.0:
        warn(f"Atom {mv[0][0] + 1} moves {mv[0][2]:.2f} Ang between the endpoints: check that the two files list the "
             "atoms in the same order (a long hop is real only if you meant it); the band is built atom by atom.")

    sections("MODEL")
    model_arg = args.custom_model if args.custom_model else args.model
    t0 = time.monotonic()
    if card is not None:
        model_card.mark_announced(args.custom_model)
    with capture_library_noise(noise, "model"):
        calc = mace_relax.get_calculator(model=model_arg, device=args.device)
    if card is not None:
        model_card.print_model_card(card, args.custom_model_card, args.custom_model, f_out, warn)
        for msg in model_card.check_structure(card, a_i):
            warn(f"Model card check: {msg}")
        kinds = ((card.get("training_data") or {}).get("kinds") or {})
        if kinds:
            print_dual(f"Training kinds    : {', '.join(f'{k} {v}' for k, v in kinds.items())} -- the card cannot "
                       "say whether defects or transition states were in the training set; a barrier is only as "
                       "good as the training data near the saddle", f_out)
    else:
        print_dual(f"Model             : MACE-MP-0 {args.model} (foundation, no card): trained on the Materials Project "
                   "PBE/PBE+U relaxations -- not on your SIESTA level of theory", f_out)
    print_dual(f"Loaded in         : {time.monotonic() - t0:.1f} s", f_out)
    label_main = "fine-tuned" if args.custom_model else f"MACE-MP-0 {args.model}"

    print(f"  ... {label_main}: endpoints, band and climbing-image NEB", flush=True)
    series = [neb_band(initial_pmg, final_pmg, calc, args, noise, label_main)]
    if compare:
        print(f"  ... MACE-MP-0 {args.model}: the same NEB", flush=True)
        with capture_library_noise(noise, "foundation model"):
            calc_f = mace_relax.get_calculator(model=args.model, device=args.device)
        series.append(neb_band(initial_pmg, final_pmg, calc_f, args, noise, f"MACE-MP-0 {args.model}"))

    sections("ENDPOINTS")
    for res in series:
        report_endpoints(res, f_out, warn)
    sections("STARTING PATH")
    for res in series:
        report_guess(res, args, f_out, warn)
    sections("NEB BAND")
    for res in series:
        report_band(res, args, f_out, warn)
    sections("BARRIER")
    for res in series:
        report_barrier(res, f_out, warn)
    print_dual(color_text("The spline is fitted to the image energies AND the forces along the band (ASE "
                          "fit_images), so it can peak between two images; the highest image is a lower bound.",
                          'cyan'), f_out)
    sections("SADDLE POINT")
    main_res = series[0]
    mv_relaxed, _ = movers(main_res["images"][0], main_res["images"][-1])
    report_saddle(main_res, mv_relaxed, f_out)
    sections("TRAINED RANGE")
    if card is not None:
        rows = domain_rows(main_res, card, warn)
        print_dual("The band against the model card's trained range:", f_out)
        print_table(["quantity", "trained", "this band", "verdict"], rows, f_out)
    else:
        print_dual("No model card (foundation model): no trained range to compare with.", f_out)
    if compare:
        sections("FINE-TUNED vs FOUNDATION")
        report_comparison(series, f_out)
        print_dual("Each model relaxes its own endpoints and band in the same cell; the foundation model is not "
                   "trained on your level of theory, so a large difference is expected -- it shows what the "
                   "fine-tuning changed, not which one is right. A SIESTA NEB (stb-neb) decides that.", f_out)

    com = None
    if args.custom_model or args.committee_models:
        sections("COMMITTEE UNCERTAINTY")
        members, note = ([], "--no-committee") if not args.committee else committee_members(
            card, args.custom_model_card, args.mlff_workflow, args.committee_models)
        if len(members) >= 2:
            print(f"  ... committee: {len(members)} models on the band", flush=True)
            print_dual(f"Committee         : {note}", f_out)
            com = committee_band(main_res, members, args.device, noise)
            report_committee(com, main_res, card, f_out, warn)
        else:
            print_dual(f"Not evaluated     : {note}.", f_out)

    results = {}
    if args.siesta_check:
        sections("SIESTA CHECK")
        written = 0
        for i in check_images(main_res, args.siesta_check_images):
            folder = os.path.join(args.output_dir, SIESTA_DIR, f"image_{i:02d}")
            r = read_siesta_folder(folder, main_res["images"][i])
            if r == "missing":
                write_siesta_folder(folder, main_res["images"][i], tpl, species_meta, args.kgrid)
                r, written = "pending", written + 1
            elif r == "stale":
                warn(f"SIESTA check: {folder} was written for another geometry (an earlier band) -- left untouched and "
                     "not compared; delete it, or use another --output-dir, to check this band.")
            results[i] = r
        print_dual(f"Template          : {tpl['source']}; k-grid {image_kgrid(tpl, main_res['images'][0], args.kgrid)[1]}",
                   f_out)
        if written:
            print_dual(f"Written           : {written} single-point folder(s) in {args.output_dir}/{SIESTA_DIR}/ -- run "
                       "SIESTA in them, then this command again:", f_out)
            print_dual(f"  for d in {args.output_dir}/{SIESTA_DIR}/image_*/; do (cd $d && mpirun -np 4 siesta < "
                       f"{os.path.basename(tpl['calc'])} > siesta.out); done", f_out)
        report_siesta(main_res, results, f_out, warn, args)
    if args.to_mlff:
        sections("EXPORT TO THE MLFF WORKFLOW")
        export_to_mlff(args, state, args.mlff_workflow, main_res, tpl, species_meta, results, warn, f_out)

    # ---- files
    p = os.path.join(args.output_dir, "neb_profile.png")
    dd = ((card or {}).get("domain") or {}).get("d_min_ang") or {}
    done = {i: r for i, r in results.items() if isinstance(r, dict)}
    siesta_pts = ([(main_res["path"][i], done[i]["E"] - done[0]["E"]) for i in sorted(done)] if 0 in done else None)
    figures.append(plot_band(series, dd.get("min"), args.fmax, p, siesta_pts))
    files.append(("Figure", p))
    from ase.io import write as ase_write
    for k, res in enumerate(series):
        suffix = "" if k == 0 else "_foundation"
        p = os.path.join(args.output_dir, f"neb_path{suffix}.xyz")
        ase_write(p, [a.copy() for a in res["images"]])
        files.append(("Relaxed band", p))
        if args.save_data:
            p = os.path.join(args.output_dir, f"neb_band{suffix}.dat")
            write_band_dat(p, res)
            files.append(("Data", p))
    if args.save_images:
        for i, a in enumerate(main_res["images"]):
            st = structure_io.from_pymatgen(AseAtomsAdaptor.get_structure(a), species_meta=species_meta,
                                            coord_format="fractional")
            structure_io.write_fdf(st, os.path.join(args.output_dir, f"image_{i:02d}.fdf"))
        files.append(("Images", f"{args.output_dir}/image_00.fdf .. image_{len(main_res['images']) - 1:02d}.fdf"))
    if args.save_gnuplot:
        plot_dir = os.path.join(args.output_dir, PLOT_DIR)
        scripts = gnuplot_files(plot_dir, series, args.fmax)
        pdfs, msg = run_gnuplot(plot_dir, scripts)
        files.append(("gnuplot folder", f"{plot_dir}/ ({', '.join(scripts)})"))
        for p in pdfs:
            files.append(("gnuplot figure", p))
        if msg:
            warn(msg)

    sections("SUMMARY & FILES")
    status = "OK" if all(r["converged"] for r in series) else "NEB not converged"
    print_dual(f"Status            : {status}" + (f" ({len(warns)} warning(s), see the last section)" if warns else ""),
               f_out)
    for res in series:
        print_dual(f"{res['label']:<18}: barrier {res['barrier']:.4f} eV forward, {res['barrier'] - res['dE']:.4f} "
                   f"eV backward, reaction energy {res['dE']:+.4f} eV", f_out)
    if com is not None:
        b = com["barriers"]
        print_dual(f"Committee         : {com['n']} members, barrier {b.min():.4f}-{b.max():.4f} eV, largest force spread "
                   f"{1000 * com['f_spread'].max():.0f} meV/Ang (image {int(np.argmax(com['f_spread']))})", f_out)
    done = {i: r for i, r in results.items() if isinstance(r, dict)}
    n_last = len(main_res["images"]) - 1
    if all(i in done for i in (0, main_res["i_top"], n_last)):
        print_dual(f"SIESTA check      : barrier on the model's path {done[main_res['i_top']]['E'] - done[0]['E']:.4f} eV "
                   f"(model {main_res['barrier_img']:.4f} eV)", f_out)
    elif args.siesta_check:
        print_dual(f"SIESTA check      : {sum(1 for r in results.values() if r == 'pending')} folder(s) waiting for SIESTA",
                   f_out)
    for kind, path in files:
        print_dual(f"{kind:<18}: {path}", f_out)
    if report_path:
        print_dual(f"Report            : {report_path}", f_out)
    not_asked = [f for f, on in (("--show", args.show), ("--save-gnuplot", args.save_gnuplot),
                                 ("--save-data", args.save_data), ("--save-images", args.save_images),
                                 ("--save-report", args.save_report), ("--siesta-check", args.siesta_check),
                                 ("--to-mlff", args.to_mlff)) if not on]
    if not_asked:
        print_dual(f"Also available    : {' '.join(not_asked)}", f_out)
    print_dual(f"Wall time         : {time.monotonic() - t_start:.1f} s", f_out)
    print_dual("Reading           : an ML-level barrier, as good as the training data near the saddle; confirm a "
               "number you will use with a SIESTA NEB (stb-neb, 5.9) or SIESTA single points on these images.", f_out)
    if args.show and not display_available():
        warn("--show needs a display (none found, or a non-interactive matplotlib backend): the figure was written as "
             "a PNG file only.")
    sections("WARNINGS")
    if warns:
        for i, w in enumerate(warns, 1):
            print_dual(color_text(f"  {i}. {w}", 'yellow'), f_out)
    else:
        print_dual(color_text("  none", 'green'), f_out)
    if noise:
        print_dual(f"Library messages  : {len(noise)} distinct message(s) from mace/torch/ase captured, not shown", f_out)
    if f_out:
        f_out.close()
    print("\n" + "-" * 60)
    print(color_text("ML NEB complete -- a fast estimate, not a substitute for a DFT-level barrier.\n", 'bold'))
    import matplotlib.pyplot as plt
    if args.show and display_available() and figures:
        plt.show()
    for fig in figures:
        plt.close(fig)


if __name__ == "__main__":
    main()
