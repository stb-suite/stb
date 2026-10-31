#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""stb-mladsorb (7.6): adsorption-site screening driven entirely by a MACE potential -- the
MACE-MP-0 foundation model, or a model fine-tuned on SIESTA data with the 7.1 workflow
(--custom-model, always with its model card).

What one run does:
  1. the clean slab: layers, the model's forces on it as given and, with --relax-depth D, its
     relaxation (atoms within D Ang of the top surface free, deeper ones fixed);
  2. the references: each adsorbate isolated in a box (relaxed) and, with --bulk-reference, the
     chemical potential of each element from its bulk (positions + cell relaxed);
  3. the candidate sites (pymatgen's AdsorbateSiteFinder): ontop/bridge/hollow, reduced by the
     slab's symmetry; each site is described by the atom right below it (a hollow above a
     second-layer atom is not the one above a fourth-layer atom: hcp vs fcc on a metal, T4 vs H3
     on Si(111)) and counted in the cell (its multiplicity);
  4. one relaxation per (adsorbate, site, height): the adsorbate and the free substrate move;
     each relaxed structure is then CLASSIFIED by where the adsorbate ended up -- an adsorbate
     placed on a bridge that slides into a hollow is reported as a hollow -- and the runs are
     grouped into distinct minima (equivalent end site and the same energy within 0.02 eV), with
     their multiplicity and Boltzmann population at --temperature;
  5. the neighbours of the most stable minimum: how far each other minimum is (nearest copy of its
     site, and in the .fdf written for it) -- the input of a diffusion barrier, which is a stb-mlneb
     (7.5) run between two of those .fdf files (with its committee, SIESTA check and export to
     training); this tool computes no barrier itself;
  6. optional: a binding curve E_ads(h) at FIXED anchor height (--height-sweep), the committee of
     the fine-tuned model on every minimum (spread of E_ads and of the site ranking), and the same
     screening with the foundation model;
  7. every structure against the model card's trained range; all warnings in the last section.

E_ads = E(slab + adsorbate) - E(slab) - E(isolated adsorbate), all with the same model. With a
fine-tuned model the isolated adsorbate is a structure the model never saw (an atom or molecule
in vacuum: its energy is little more than the model's E0s), so the report also gives the energy
relative to the best site (independent of the reference) and, with --bulk-reference, relative to
the elements' bulk phases (in the training domain of a model trained on that bulk).

Implementation notes. The isolated adsorbate is evaluated with a calculator of its own: a
TorchScript-compiled MACE model on CUDA, once called on a structure with no neighbour within the
cutoff (an isolated atom), fails on every later call ("Global alloc not supported yet").
"""

from stb import __version__ as VERSION

import os
import sys
import time
import argparse
import contextlib
from datetime import datetime

import numpy as np

from stb.core import structure_io, model_card
from stb.core.cli import (color_text, show_intro, print_dual, print_table, capture_library_noise,
                          display_available, Sections, run_gnuplot, print_progress_line, finish_progress_line)
from stb.core.deps import require_mace

REPORT_FILE = "stb_mladsorb_report.txt"
FIGURE_FILE = "mladsorb.png"
PLOT_DIR = "plot"
LAYER_TOL = 0.35          # Ang: atoms whose heights differ by less belong to one layer
SITE_MATCH = 0.6          # Ang: in-plane distance under which a relaxed adsorbate sits AT a site
BELOW_MATCH = 0.5         # Ang: in-plane distance under which a substrate atom is 'below' a site
SYM_MATCH = 0.3           # Ang: in-plane tolerance of the symmetry mapping of sites
ENERGY_DUP = 0.02         # eV: two runs ending at one site within this are one minimum
TOO_CLOSE = 0.8           # x the sum of covalent radii: a held adsorbate closer than this to an atom is unphysical
KB_EV = 8.617333262e-5
_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#17becf"]
TYPE_MARKERS = {"ontop": "o", "bridge": "s", "hollow": "^"}


def _quiet(noise, label):
    return capture_library_noise(noise, label) if noise is not None else contextlib.nullcontext()


def max_force(atoms, idx=None):
    f = np.linalg.norm(atoms.get_forces(), axis=1)
    if idx is not None:
        f = f[list(idx)] if len(idx) else np.zeros(1)
    return float(f.max())


# ---- surface geometry ---------------------------------------------------------------------

class SurfaceFrame:
    """In-plane / normal geometry of a slab with its vacuum along c: heights along the surface
    normal n = a x b, in-plane vectors periodic in a and b (minimum image), 2D plot coordinates."""

    def __init__(self, cell):
        cell = np.asarray(cell, float)
        self.a, self.b, self.c = cell
        n = np.cross(self.a, self.b)
        self.area = float(np.linalg.norm(n))
        self.n = n / self.area
        self.e1 = self.a / np.linalg.norm(self.a)
        self.e2 = np.cross(self.n, self.e1)

    def height(self, xyz):
        return np.asarray(xyz, float) @ self.n

    def _uv(self, xyz):
        p = np.asarray(xyz, float)
        par = p - (p @ self.n) * self.n
        sol, *_ = np.linalg.lstsq(np.array([self.a, self.b]).T, par, rcond=None)
        return sol

    def inplane(self, p_from, p_to):
        """Shortest in-plane vector (Cartesian) from one point to another, periodic in a and b."""
        d = self._uv(p_to) - self._uv(p_from)
        d -= np.round(d)
        best = None
        for i in (-1, 0, 1):
            for j in (-1, 0, 1):
                v = (d[0] + i) * self.a + (d[1] + j) * self.b
                if best is None or v @ v < best @ best:
                    best = v
        return best

    def dist(self, p, q):
        return float(np.linalg.norm(self.inplane(p, q)))

    def uv_frac(self, xyz):
        """In-plane fractional coordinates of a point, in [0, 1)."""
        return self._uv(xyz) % 1.0

    def xy(self, xyz):
        p = np.atleast_2d(np.asarray(xyz, float))
        return np.c_[p @ self.e1, p @ self.e2]


def find_layers(heights, tol=LAYER_TOL):
    """Layers from the top: [{'idx': [...], 'h': mean height}], layer 1 the outermost."""
    order = np.argsort(-np.asarray(heights))
    layers = []
    for i in order:
        if layers and layers[-1]["last"] - heights[i] <= tol:
            layers[-1]["idx"].append(int(i))
            layers[-1]["last"] = heights[i]
        else:
            layers.append({"idx": [int(i)], "last": heights[i]})
    return [{"idx": L["idx"], "h": float(np.mean([heights[i] for i in L["idx"]]))} for L in layers]


def layer_index(layers, n):
    out = np.zeros(n, int)
    for k, L in enumerate(layers, 1):
        out[L["idx"]] = k
    return out


def describe_below(frame, pos, sub_pos, layer_of, symbols, max_layer=4):
    """(layer, text) of the topmost substrate atom right below an in-plane position (layers 1-4)."""
    best = None
    for i, p in enumerate(sub_pos):
        if layer_of[i] > max_layer:
            continue
        if frame.dist(pos, p) < BELOW_MATCH and (best is None or layer_of[i] < layer_of[best]):
            best = i
    if best is None:
        return None, f"nothing below (L1-L{max_layer})"
    return int(layer_of[best]), f"above L{layer_of[best]} {symbols[best]}{best + 1}"


def find_sites(slab_pmg, frame, site_types, distance, symprec, layer_of, symbols):
    """The candidate sites: `slots` (symmetry-reduced, the starting points), `all_sites` (every
    site in the cell, each mapped to its slot by the slab's symmetry operations that keep the
    surface normal), and the number of such operations."""
    from pymatgen.analysis.adsorption import AdsorbateSiteFinder
    from pymatgen.core.operations import SymmOp
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    finder = AdsorbateSiteFinder(slab_pmg)
    reduced = finder.find_adsorption_sites(distance=distance, symm_reduce=symprec, positions=site_types)
    full = finder.find_adsorption_sites(distance=distance, symm_reduce=0, positions=site_types)
    sub_pos = slab_pmg.cart_coords
    ops = [SymmOp.from_rotation_and_translation(np.eye(3), np.zeros(3))]
    try:
        sga = SpacegroupAnalyzer(slab_pmg, symprec=symprec)
        ops += [op for op in sga.get_symmetry_operations(cartesian=True)
                if (op.rotation_matrix @ frame.n) @ frame.n > 0.99]
    except Exception:
        pass

    def equivalent(p, q):
        return any(frame.dist(op.operate(p), q) < SYM_MATCH for op in ops)

    # pymatgen's labels and its list of every site follow a Delaunay triangulation, which is degenerate on a honeycomb
    # or a square net: one position can come out as a 'bridge' in one copy of the cell and a 'hollow' in another, and
    # whole sets of equivalent sites can be missing (on Si(111) 2x2 only 4 of the 12 B2 bridges). So the classes are
    # built from positions: the reduced sites are merged when equivalent (a class keeps the label of its first site),
    # and every site of the cell is generated as the orbit of its class under the slab's symmetry operations -- the
    # multiplicity is the orbit's size.
    slots = []
    for st in site_types:
        for c in reduced.get(st, []):
            c = np.array(c, float)
            if not any(equivalent(c, sl["coord"]) for sl in slots):
                slots.append({"type": st, "coord": c})
    # Numbered by geometry, not by pymatgen's order (which follows the triangulation and changes when the same slab is
    # written with its atoms wrapped differently): type, then the layer of the atom right below, then the in-plane
    # distances to the nearest atoms of the top two layers -- the same surface always gets the same numbers.
    top2 = [i for i in range(len(sub_pos)) if layer_of[i] <= 2]

    def key(sl):
        below, _ = describe_below(frame, sl["coord"], sub_pos, layer_of, symbols)
        d = sorted(round(frame.dist(sl["coord"], sub_pos[i]), 1) for i in top2)[:4]
        return (site_types.index(sl["type"]), below if below is not None else 99, d)

    slots.sort(key=key)
    for k, sl in enumerate(slots, 1):
        sl["slot"] = k

    def in_cell(p):
        u, v = frame.uv_frac(p)
        return u * frame.a + v * frame.b + frame.height(p) * frame.n

    all_sites = []
    for sl in slots:
        orbit = []
        for op in ops:
            q = in_cell(op.operate(sl["coord"]))
            if not any(frame.dist(q, x) < SYM_MATCH for x in orbit):
                orbit.append(q)
        all_sites += [{"type": sl["type"], "coord": q, "slot": sl["slot"]} for q in orbit]
    missed = 0
    for st in site_types:
        for c in full.get(st, []):
            if not any(frame.dist(np.array(c, float), x["coord"]) < SYM_MATCH for x in all_sites):
                missed += 1
                all_sites.append({"type": st, "coord": np.array(c, float), "slot": None})
    for sl in slots:
        sl["mult"] = sum(1 for s in all_sites if s["slot"] == sl["slot"])
        sl["below_layer"], sl["below"] = describe_below(frame, sl["coord"], sub_pos, layer_of, symbols)
        sl["uv"] = frame.uv_frac(sl["coord"])
    return finder, slots, all_sites, len(ops) - 1


def classify(frame, pos, all_sites):
    """Where a relaxed adsorbate (its anchor atom) sits: the nearest site of the cell."""
    d = [frame.dist(pos, s["coord"]) for s in all_sites]
    k = int(np.argmin(d))
    s = all_sites[k]
    return {"slot": s["slot"] if d[k] < SITE_MATCH else None, "type": s["type"], "dist": float(d[k]),
            "near_slot": s["slot"]}


def site_name(slot, slots_by_id, cls=None):
    if slot is None:
        if cls is None:
            return "off-site"
        return f"off-site (near #{cls['near_slot']} {cls['type']}, {cls['dist']:.2f} A)"
    return f"#{slot} {slots_by_id[slot]['type']}"


# ---- one model: slab, references, relaxations ------------------------------------------------

def prepare_slab(slab_given, calc, args, noise, label):
    """The clean slab with the model: energy and forces as given, then (--relax-depth > 0) relaxed
    with every atom deeper than the depth fixed."""
    from ase.constraints import FixAtoms
    from stb.core import mace_relax
    at = slab_given.copy()
    frame = SurfaceFrame(at.cell)
    h = frame.height(at.positions)
    top = h.max()
    free = [i for i in range(len(at)) if args.relax_depth > 0 and h[i] >= top - args.relax_depth - 1e-6]
    fixed = [i for i in range(len(at)) if i not in set(free)]
    at.calc = calc
    with _quiet(noise, label):
        e0 = at.get_potential_energy()
        f_all = np.linalg.norm(at.get_forces(), axis=1)
    out = {"e_given": e0, "f_given": float(f_all.max()), "f_given_atom": int(f_all.argmax()),
           "f_given_free": float(f_all[free].max()) if free else None, "free": free, "fixed": fixed,
           "converged": None, "steps": 0, "moved": 0.0}
    if free:
        if fixed:
            at.set_constraint(FixAtoms(indices=fixed))
        p0 = at.get_positions().copy()
        t0 = time.monotonic()
        with _quiet(noise, label):
            out["converged"], out["steps"] = mace_relax.relax(at, calc, fmax=args.fmax, max_steps=args.max_steps)
        out["seconds"] = time.monotonic() - t0
        out["moved"] = float(np.linalg.norm(at.get_positions() - p0, axis=1).max())
    with _quiet(noise, label):
        out["e"] = at.get_potential_energy()
        out["f_after"] = max_force(at, free) if free else out["f_given"]
    out["atoms"] = at
    out["given"] = slab_given
    return out


def isolated_reference(name, mol, calc_factory, args, noise, label):
    """The adsorbate alone in a cubic box, relaxed, with a calculator of its own (see the module
    docstring: an isolated atom breaks a compiled model on CUDA for every later call)."""
    from pymatgen.io.ase import AseAtomsAdaptor
    from stb.adsorb import isolated_adsorbate_structure, molecule_extent
    from stb.core import mace_relax
    iso = AseAtomsAdaptor.get_atoms(structure_io.to_pymatgen(isolated_adsorbate_structure(mol, args.vacuum_box)))
    with _quiet(noise, label):
        calc = calc_factory()
        iso.calc = calc
        if len(iso) > 1:
            conv, steps = mace_relax.relax(iso, calc, fmax=args.fmax, max_steps=args.max_steps)
        else:
            conv, steps = True, 0
        e = iso.get_potential_energy()
    return {"name": name, "e": e, "converged": conv, "steps": steps, "n": len(iso),
            "extent": molecule_extent(mol), "atoms": iso}


def bulk_reference(path, calc, args, noise, label):
    """Chemical potential of an element: its bulk relaxed (positions + cell) with the model, E/N."""
    from pymatgen.io.ase import AseAtomsAdaptor
    from stb.core import mace_relax
    at = AseAtomsAdaptor.get_atoms(structure_io.to_pymatgen(structure_io.read_fdf(path)))
    v0 = at.get_volume() / len(at)
    with _quiet(noise, label):
        conv, steps = mace_relax.relax(at, calc, cell_mask=[True] * 6, fmax=min(args.fmax, 0.01),
                                       max_steps=max(args.max_steps, 300))
        e = at.get_potential_energy()
    return {"file": path, "element": at.get_chemical_symbols()[0], "n": len(at), "v0": v0,
            "v": at.get_volume() / len(at), "mu": e / len(at), "converged": conv, "steps": steps, "atoms": at}


def start_structure(finder, mol, coord, calc, args, rng, noise, label):
    """The adsorbate placed at `coord` (pymatgen), optionally the best of --n-orientations random
    rotations by a single-point energy (no relaxation)."""
    from pymatgen.core import Molecule
    from pymatgen.io.ase import AseAtomsAdaptor
    if args.n_orientations <= 1 or len(mol) == 1:
        return AseAtomsAdaptor.get_atoms(finder.add_adsorbate(mol, coord)), None
    from scipy.spatial.transform import Rotation
    best = None
    energies = []
    for _ in range(args.n_orientations):
        rot = Rotation.random(rng=rng)
        cm = mol.center_of_mass
        trial = Molecule(mol.species, rot.apply(mol.cart_coords - cm) + cm)
        at = AseAtomsAdaptor.get_atoms(finder.add_adsorbate(trial, coord))
        at.calc = calc
        with _quiet(noise, label):
            e = at.get_potential_energy()
        energies.append(e)
        if best is None or e < best[0]:
            best = (e, at)
    return best[1], (min(energies), max(energies))


def analyse_config(at, start, n_sub, frame, slab_ref, h_top, all_sites):
    """Where the adsorbate ended up and how it is bound."""
    from ase.data import covalent_radii, atomic_numbers
    from ase.geometry import find_mic
    pos = at.get_positions()
    ads = list(range(n_sub, len(at)))
    hs = frame.height(pos[ads])
    anchor = ads[int(np.argmin(hs))]
    cls = classify(frame, pos[anchor], all_sites)
    sym = at.get_chemical_symbols()
    d = at.get_distances(anchor, list(range(n_sub)), mic=True)
    r_a = covalent_radii[atomic_numbers[sym[anchor]]]
    cut = np.array([1.25 * (r_a + covalent_radii[atomic_numbers[sym[j]]]) for j in range(n_sub)])
    nn = np.sort(d[d < cut])
    dsub, _ = find_mic(pos[:n_sub] - slab_ref.get_positions(), at.cell, at.pbc)
    return {"anchor": anchor, "anchor_sym": sym[anchor], "cls": cls,
            "shift": frame.dist(start.get_positions()[anchor], pos[anchor]),
            "height": float(frame.height(pos[anchor]) - h_top), "cn": int(len(nn)),
            "nn": [float(x) for x in nn[:4]], "d_nearest": float(d.min()),
            "sub_disp": float(np.linalg.norm(dsub, axis=1).max()),
            "dmin": model_card.shortest_distance(at)}


def relax_config(start, fixed, calc, args, noise, label, extra_fixed=()):
    from ase.constraints import FixAtoms
    from stb.core import mace_relax
    at = start.copy()
    fx = sorted(set(fixed) | set(extra_fixed))
    if fx:
        at.set_constraint(FixAtoms(indices=fx))
    t0 = time.monotonic()
    with _quiet(noise, label):
        conv, steps = mace_relax.relax(at, calc, fmax=args.fmax, max_steps=args.max_steps)
        e = at.get_potential_energy()
        fres = float(np.linalg.norm(at.get_forces(), axis=1).max())
    return at, e, bool(conv), int(steps), fres, time.monotonic() - t0


def group_minima(configs):
    """Distinct minima: runs of one adsorbate that end at the same site (same slot, or within
    SITE_MATCH in-plane when off-site) with energies within ENERGY_DUP. Sorted by E_ads per
    adsorbate order of appearance, most stable first."""
    minima = []
    for c in sorted(configs, key=lambda r: r["e_ads"]):
        hit = None
        for m in minima:
            if m["name"] != c["name"] or abs(c["e_ads"] - m["e_ads"]) > ENERGY_DUP:
                continue
            if c["end_slot"] is not None and c["end_slot"] == m["end_slot"]:
                hit = m
            elif c["end_slot"] is None and m["end_slot"] is None and c["frame"].dist(
                    c["atoms"].positions[c["anchor"]], m["atoms"].positions[m["anchor"]]) < SITE_MATCH:
                hit = m
            if hit:
                break
        if hit:
            hit["starts"].append(c["id"])
        else:
            m = dict(c)
            m["starts"] = [c["id"]]
            minima.append(m)
    return minima


def boltzmann(minima, slots_by_id, temperature):
    """Population of each minimum of one adsorbate: multiplicity x exp(-dE/kT), normalized."""
    if not minima:
        return
    e0 = minima[0]["e_ads"]
    w = []
    for m in minima:
        mult = slots_by_id[m["end_slot"]]["mult"] if m["end_slot"] is not None else 1
        m["mult"] = mult
        w.append(mult * np.exp(-(m["e_ads"] - e0) / (KB_EV * temperature)))
    w = np.array(w)
    for m, p in zip(minima, w / w.sum()):
        m["pop"] = float(p)


def screen(label, calc, calc_factory, slab_given, ctx, args, adsorbates, noise, rng):
    """The whole screening with one model. The candidate sites are found once, on the first
    model's relaxed slab (ctx['sites'] is None then), and reused by the others so their runs start
    from the same sites."""
    from pymatgen.io.ase import AseAtomsAdaptor
    from pymatgen.analysis.adsorption import AdsorbateSiteFinder
    res = {"label": label}
    t0 = time.monotonic()
    print(f"  ... {label}: clean slab", flush=True)
    slab = prepare_slab(slab_given, calc, args, noise, label)
    res["slab"] = slab
    frame = SurfaceFrame(slab["atoms"].cell)
    n_sub = len(slab["atoms"])
    heights = frame.height(slab["atoms"].positions)
    layers = find_layers(heights)
    layer_of = layer_index(layers, n_sub)
    symbols = slab["atoms"].get_chemical_symbols()
    res.update(frame=frame, layers=layers, layer_of=layer_of, h_top=layers[0]["h"])
    slab_pmg = AseAtomsAdaptor.get_structure(slab["atoms"])
    if ctx.get("sites") is None:
        finder, slots, all_sites, nops = find_sites(slab_pmg, frame, ctx["site_types"], args.height,
                                                    args.symprec, layer_of, symbols)
        ctx["sites"] = (slots, all_sites, nops)
    else:
        finder = AdsorbateSiteFinder(slab_pmg)
    slots, all_sites, _ = ctx["sites"]
    res["iso"] = {}
    for name, mol in adsorbates:
        print(f"  ... {label}: isolated {name}", flush=True)
        res["iso"][name] = isolated_reference(name, mol, calc_factory, args, noise, label)
    res["bulk"] = {}
    for path in args.bulk_reference or []:
        b = bulk_reference(path, calc, args, noise, label)
        res["bulk"][b["element"]] = b

    def mu_sum(mol):
        try:
            return sum(res["bulk"][str(sp.symbol)]["mu"] for sp in mol.species)
        except KeyError:
            return None

    if not slots:
        res["configs"], res["minima"], res["sweep"] = [], [], []
        return res
    jobs = [(name, mol, sl, h) for name, mol in adsorbates for sl in slots for h in ctx["start_heights"]]
    configs = []
    for k, (name, mol, sl, h) in enumerate(jobs, 1):
        print_progress_line(f"  ... {label}: relaxation {k}/{len(jobs)} ({name} on #{sl['slot']} {sl['type']}, "
                            f"start {h:.2f} A)", k, len(jobs))
        coord = sl["coord"] + (h - args.height) * frame.n
        start, orient = start_structure(finder, mol, coord, calc, args, rng, noise, label)
        start.set_constraint()
        at, e, conv, steps, fres, secs = relax_config(start, slab["fixed"], calc, args, noise, label)
        info = analyse_config(at, start, n_sub, frame, slab["atoms"], res["h_top"], all_sites)
        mu = mu_sum(mol)
        configs.append({"id": k, "name": name, "start_slot": sl["slot"], "start_type": sl["type"], "h0": h,
                        "atoms": at, "e": e, "e_ads": e - slab["e"] - res["iso"][name]["e"],
                        "e_form": (e - slab["e"] - mu) if mu is not None else None,
                        "converged": conv, "steps": steps, "fres": fres, "seconds": secs, "orient": orient,
                        "end_slot": info["cls"]["slot"], "end_type": info["cls"]["type"], "frame": frame, **info})
    finish_progress_line()
    res["configs"] = configs
    minima = group_minima(configs)
    slots_by_id = {s["slot"]: s for s in slots}
    for name, _ in adsorbates:
        boltzmann([m for m in minima if m["name"] == name], slots_by_id, args.temperature)
    for m in minima:
        best = min(x["e_ads"] for x in minima if x["name"] == m["name"])
        m["de"] = m["e_ads"] - best
    res["minima"] = minima
    res["sweep"] = []
    if args.height_sweep is not None:
        from stb.adsorb import build_sweep_values
        from ase.data import covalent_radii, atomic_numbers
        hs = build_sweep_values(*args.height_sweep)
        jobs = [(name, mol, sl, h) for name, mol in adsorbates for sl in slots for h in hs]
        for k, (name, mol, sl, h) in enumerate(jobs, 1):
            print_progress_line(f"  ... {label}: binding curve {k}/{len(jobs)} ({name} on #{sl['slot']}, "
                                f"{h:.2f} A)", k, len(jobs))
            coord = sl["coord"] + (h - args.height) * frame.n
            start, _ = start_structure(finder, mol, coord, calc, args, rng, noise, label)
            start.set_constraint()
            ads = list(range(n_sub, len(start)))
            anchor = ads[int(np.argmin(frame.height(start.positions[ads])))]
            at, e, conv, steps, fres, secs = relax_config(start, slab["fixed"], calc, args, noise, label,
                                                          extra_fixed=[anchor])
            d = at.get_distances(anchor, list(range(n_sub)), mic=True)
            j = int(np.argmin(d))
            limit = TOO_CLOSE * (covalent_radii[atomic_numbers[at[anchor].symbol]]
                                 + covalent_radii[atomic_numbers[at[j].symbol]])
            res["sweep"].append({"name": name, "slot": sl["slot"], "type": sl["type"], "h0": h,
                                 "height": float(frame.height(at.positions[anchor]) - res["h_top"]),
                                 "e_ads": e - slab["e"] - res["iso"][name]["e"], "converged": conv, "steps": steps,
                                 "d_nearest": float(d[j]), "too_close": bool(d[j] < limit), "atoms": at})
        finish_progress_line()
    res["seconds"] = time.monotonic() - t0
    return res


# ---- committee --------------------------------------------------------------------------------

def committee_minima(res, members, device, noise, adsorbates):
    """Every committee member on the main model's geometries: clean slab, isolated adsorbates
    (each with its own calculator) and every distinct minimum. Returns E_ads per member and
    minimum, the energy spread, and the force spread per minimum (max over the moving atoms)."""
    from stb.core import mace_relax
    mins = res["minima"]
    moving = [i for i in range(len(mins[0]["atoms"])) if i not in set(res["slab"]["fixed"])] if mins else []
    e_ads, e_form, forces = [], [], []
    for model, card in members:
        model_card.register(model, card, echo=False)
        model_card.mark_announced(model)
        with _quiet(noise, "committee"):
            calc = mace_relax.get_calculator(model=model, device=device)
            s = res["slab"]["atoms"].copy()
            s.calc = calc
            e_slab = s.get_potential_energy()
            mu = {}
            for el, b in res["bulk"].items():
                bb = b["atoms"].copy()
                bb.calc = calc
                mu[el] = bb.get_potential_energy() / len(bb)
            e_iso = {}
            for name, _ in adsorbates:
                c_iso = mace_relax.get_calculator(model=model, device=device)
                a = res["iso"][name]["atoms"].copy()
                a.calc = c_iso
                e_iso[name] = a.get_potential_energy()
            row, rowf, fr = [], [], []
            for m in mins:
                a = m["atoms"].copy()
                a.set_constraint()
                a.calc = calc
                e = a.get_potential_energy()
                row.append(e - e_slab - e_iso[m["name"]])
                syms = [str(sp) for sp in a.get_chemical_symbols()[len(s):]]
                rowf.append(e - e_slab - sum(mu[x] for x in syms) if all(x in mu for x in syms) else np.nan)
                fr.append(a.get_forces()[moving])
        e_ads.append(row)
        e_form.append(rowf)
        forces.append(fr)
    e_ads = np.array(e_ads)                                   # (n_members, n_minima)
    e_form = np.array(e_form)
    spread = []
    for j in range(len(mins)):
        f = np.array([forces[k][j] for k in range(len(members))])
        s_atom = np.sqrt(np.mean(np.sum((f - f.mean(axis=0)) ** 2, axis=2), axis=0))
        spread.append(float(s_atom.max()))
    return {"n": len(members), "e_ads": e_ads, "e_std": e_ads.std(axis=0), "e_form": e_form,
            "form_std": e_form.std(axis=0), "f_spread": np.array(spread)}


# ---- neighbours of the most stable minimum ---------------------------------------------------

def neighbours(res, name):
    """For every other minimum of `name`: the in-plane distance from minimum 1 to the nearest copy of
    its site anywhere in the cell, and to the adsorbate as it sits in that minimum's own structure
    (the .fdf written for it, minimum image) -- a stb-mlneb hop between the two files is the latter."""
    frame = res["frame"]
    mins = [m for m in res["minima"] if m["name"] == name]
    if len(mins) < 2:
        return []
    a = mins[0]
    pa = a["atoms"].positions[a["anchor"]]
    out = []
    for r, m in enumerate(mins[1:], 2):
        sites = [x for x in res["all_sites"] if x["slot"] == m["end_slot"]] if m["end_slot"] is not None else []
        d_site = min(frame.dist(pa, x["coord"]) for x in sites) if sites else None
        out.append({"rank": r, "min": m, "d_site": d_site, "d_file": frame.dist(pa, m["atoms"].positions[m["anchor"]])})
    return out


# ---- figures and files --------------------------------------------------------------------------

def _cell_outline(frame):
    corners = np.array([[0, 0, 0], frame.a, frame.a + frame.b, frame.b, [0, 0, 0]], float)
    return frame.xy(corners)


def _slot_color(slot):
    return _COLORS[(slot - 1) % len(_COLORS)] if slot else "0.5"


def _in_cell_xy(frame, xyz):
    """2D plot coordinates of a point brought into the cell (in-plane fractional in [0, 1))."""
    u, v = frame.uv_frac(xyz)
    return frame.xy(u * frame.a + v * frame.b)[0]


def panel_top(ax, res, slots, all_sites):
    """Top view of the cell and a margin: substrate (top four layers), every site (colour = its
    symmetry-distinct site), the distinct minima (stars, numbered by rank), the moves of the runs that
    left their starting site (arrows)."""
    frame = res["frame"]
    at = res["slab"]["atoms"]
    lay = res["layer_of"]
    sizes = {1: 130, 2: 75, 3: 40, 4: 25}
    shades = {1: "0.15", 2: "0.45", 3: "0.65", 4: "0.82"}
    out = _cell_outline(frame)
    lo, hi = out.min(axis=0) - 0.8, out.max(axis=0) + 0.8
    shifts = [frame.xy(i * frame.a + j * frame.b)[0] for i in (-1, 0, 1) for j in (-1, 0, 1)]

    def replicas(p):
        pts = np.array([p + s_ for s_ in shifts])
        return pts[np.all((pts >= lo) & (pts <= hi), axis=1)]

    for L in (4, 3, 2, 1):
        idx = [i for i in range(len(at)) if lay[i] == L]
        if not idx:
            continue
        pts = np.vstack([replicas(_in_cell_xy(frame, at.positions[i])) for i in idx])
        ax.scatter(pts[:, 0], pts[:, 1], s=sizes[L], c=shades[L], edgecolors="k", linewidths=0.3,
                   label=f"substrate L{L}", zorder=1 + (5 - L))
    ax.plot(out[:, 0], out[:, 1], "k-", lw=0.8)
    for sl in slots:
        pts = [replicas(_in_cell_xy(frame, x["coord"])) for x in all_sites if x["slot"] == sl["slot"]]
        pts = np.vstack([p for p in pts if len(p)]) if any(len(p) for p in pts) else np.empty((0, 2))
        ax.scatter(pts[:, 0], pts[:, 1], marker=TYPE_MARKERS.get(sl["type"], "x"), s=45, facecolors="none",
                   edgecolors=_slot_color(sl["slot"]), linewidths=1.5, zorder=8,
                   label=f"#{sl['slot']} {sl['type']} ({sl['below']})")
    for c in res["configs"]:
        if c["shift"] > SITE_MATCH:
            p1 = _in_cell_xy(frame, c["atoms"].positions[c["anchor"]])
            v = frame.inplane(c["atoms"].positions[c["anchor"]], next(x["coord"] for x in slots
                                                                     if x["slot"] == c["start_slot"]))
            p0 = p1 + np.array([v @ frame.e1, v @ frame.e2])
            ax.annotate("", xy=p1, xytext=p0, arrowprops=dict(arrowstyle="->", color="0.25", lw=1.1), zorder=9)
    for r, m in enumerate(res["minima"], 1):
        p = _in_cell_xy(frame, m["atoms"].positions[m["anchor"]])
        ax.scatter([p[0]], [p[1]], marker="*", s=150, c="gold", edgecolors="k", alpha=0.75, zorder=10,
                   label="minima (rank)" if r == 1 else None)
        ax.annotate(str(r), p, xytext=(5, 5), textcoords="offset points", fontsize=9, weight="bold", zorder=11)
    ax.set_aspect("equal")
    ax.set_xlim(lo[0], hi[0])
    ax.set_ylim(lo[1], hi[1])
    ax.set_xlabel("x (Ang)")
    ax.set_ylabel("y (Ang)")
    ax.set_title("(a) top view: sites, minima (stars), moves (arrows)")
    ax.legend(fontsize=6, loc="upper left", bbox_to_anchor=(1.01, 1.0))


def panel_minima(ax, res, found, com, slots_by_id):
    mins = res["minima"]
    labels = [f"{m['name']}\n{site_name(m['end_slot'], slots_by_id)}" for m in mins]
    x = np.arange(len(mins))
    e = [m["e_ads"] for m in mins]
    w = 0.38 if found else 0.6
    err = com["e_std"] if com is not None else None
    ax.bar(x - (w / 2 if found else 0), e, w, color=[_slot_color(m["end_slot"]) for m in mins], yerr=err,
           capsize=4, edgecolor="k", label=res["label"] + (" (error bars: committee std)" if com is not None else ""))
    if found:
        fe = []
        for m in mins:
            same = [f for f in found["minima"] if f["name"] == m["name"] and f["end_slot"] == m["end_slot"]
                    and m["end_slot"] is not None]
            fe.append(same[0]["e_ads"] if same else np.nan)
        ax.bar(x + w / 2, fe, w, color="white", edgecolor="k", hatch="//", label=found["label"])
    for xi, ei in zip(x, e):
        ax.annotate(f"{ei:.3f}", (xi - (w / 2 if found else 0), ei / 2), ha="center", va="center", fontsize=7,
                    color="white", weight="bold", rotation=90)
    ax.axhline(0, color="0.5", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel("E_ads (eV)")
    ax.set_title("(b) distinct minima (colour = site, as in (a))")
    ax.legend(fontsize=7)


def panel_configs(ax, res, slots_by_id):
    cs = sorted(res["configs"], key=lambda c: c["id"])
    y = np.arange(len(cs))
    ax.barh(y, [c["e_ads"] for c in cs], color=[_slot_color(c["end_slot"]) for c in cs], edgecolor="k")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{c['name']} #{c['start_slot']} {c['start_type']} @{c['h0']:.1f} -> "
                        f"{site_name(c['end_slot'], slots_by_id)}" for c in cs], fontsize=7)
    ax.invert_yaxis()
    ax.axvline(0, color="0.5", lw=0.8)
    ax.set_xlabel("E_ads (eV)")
    ax.set_title("(c) every relaxation: start -> end site (colour = end site)")


def panel_sweep(ax, res):
    groups = {}
    for p in res["sweep"]:
        groups.setdefault((p["name"], p["slot"], p["type"]), []).append(p)
    for (name, slot, st), pts in sorted(groups.items()):
        pts = sorted((p for p in pts if not p["too_close"]), key=lambda p: p["height"])
        ax.plot([p["height"] for p in pts], [p["e_ads"] for p in pts], marker=TYPE_MARKERS.get(st, "o"),
                color=_slot_color(slot), label=f"{name} #{slot} {st}")
    ax.axhline(0, color="0.5", lw=0.8)
    ax.set_xlabel("anchor height above layer 1 (Ang, held fixed)")
    ax.set_ylabel("E_ads (eV)")
    ax.set_title("(d) binding curve at fixed height (points too close to an atom left out)")
    ax.legend(fontsize=7)


def make_figure(res, found, com, slots, all_sites, path):
    import matplotlib.pyplot as plt
    slots_by_id = {s["slot"]: s for s in slots}
    panels = [lambda ax: panel_top(ax, res, slots, all_sites),
              lambda ax: panel_minima(ax, res, found, com, slots_by_id),
              lambda ax: panel_configs(ax, res, slots_by_id)]
    if res["sweep"]:
        panels.append(lambda ax: panel_sweep(ax, res))
    nrow = (len(panels) + 1) // 2
    fig = plt.figure(figsize=(14, 5.2 * nrow))
    grid = fig.add_gridspec(nrow, 2)
    for k, p in enumerate(panels):
        last_alone = k == len(panels) - 1 and len(panels) % 2 == 1
        p(fig.add_subplot(grid[k // 2, :] if last_alone else grid[k // 2, k % 2]))
    fig.suptitle(f"stb-mladsorb: {res['label']}", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    return fig


def write_configs_dat(path, res):
    with open(path, "w") as f:
        f.write(f"# stb-mladsorb, {res['label']}: every relaxation\n# id adsorbate start_slot start_type start_dist(A) "
                "end_slot end_type shift(A) height(A) CN d_nearest(A) E_ads(eV) E_form(eV) steps converged "
                "max|F|(eV/A)\n")
        for c in sorted(res["configs"], key=lambda c: c["id"]):
            ef = f"{c['e_form']:.6f}" if c["e_form"] is not None else "nan"
            f.write(f"{c['id']} {c['name']} {c['start_slot']} {c['start_type']} {c['h0']:.3f} "
                    f"{c['end_slot'] if c['end_slot'] is not None else 0} {c['end_type']} {c['shift']:.4f} "
                    f"{c['height']:.4f} {c['cn']} {c['d_nearest']:.4f} {c['e_ads']:.6f} {ef} {c['steps']} "
                    f"{int(c['converged'])} {c['fres']:.5f}\n")


def write_minima_dat(path, res, com=None):
    below = {s_["slot"]: s_["below_layer"] for s_ in res.get("slots", [])}
    with open(path, "w") as f:
        f.write(f"# stb-mladsorb, {res['label']}: distinct minima, most stable first\n# rank adsorbate end_slot "
                "end_type E_ads(eV) dE(eV) E_form(eV) multiplicity population height(A) CN n_starts "
                "committee_std(eV) below(layer of the atom right below the site, 0 = none in L1-L4)\n")
        for r, m in enumerate(res["minima"], 1):
            ef = f"{m['e_form']:.6f}" if m["e_form"] is not None else "nan"
            cs = f"{com['e_std'][r - 1]:.5f}" if com is not None else "nan"
            f.write(f"{r} {m['name']} {m['end_slot'] if m['end_slot'] is not None else 0} {m['end_type']} "
                    f"{m['e_ads']:.6f} {m['de']:.6f} {ef} {m['mult']} {m['pop']:.6f} {m['height']:.4f} {m['cn']} "
                    f"{len(m['starts'])} {cs} {below.get(m['end_slot']) or 0}\n")


def write_sweep_dat(path, res):
    with open(path, "w") as f:
        f.write("# stb-mladsorb binding curves at fixed anchor height; one block per (adsorbate, site)\n"
                "# height(A) E_ads(eV) start_distance(A) d_nearest(A) too_close(1: left out of the minimum)\n")
        groups = {}
        for p in res["sweep"]:
            groups.setdefault((p["name"], p["slot"], p["type"]), []).append(p)
        for (name, slot, st), pts in sorted(groups.items()):
            f.write(f"# {name} #{slot} {st}\n")
            for p in sorted(pts, key=lambda p: p["height"]):
                f.write(f"{p['height']:.4f} {p['e_ads']:.6f} {p['h0']:.3f} {p['d_nearest']:.4f} {int(p['too_close'])}\n")
            f.write("\n\n")


def gnuplot_files(plot_dir, res, found, slots, all_sites):
    os.makedirs(plot_dir, exist_ok=True)
    frame = res["frame"]
    scripts = []
    head = "# --- STB Plot Configuration ---\n# Generated by stb-mladsorb: gnuplot {name}\n" \
           'set terminal pdfcairo enhanced color font "Arial,9" size {w},{h}\nset output "{pdf}"\nset grid\n'
    # top view
    at = res["slab"]["atoms"]
    with open(os.path.join(plot_dir, "substrate.dat"), "w") as f:
        f.write("# x(A) y(A) layer   (top four layers, cell plus its neighbours)\n")
        for i, p in enumerate(at.positions):
            if res["layer_of"][i] <= 4:
                for si in (-1, 0, 1):
                    for sj in (-1, 0, 1):
                        q = frame.xy(p + si * frame.a + sj * frame.b)[0]
                        f.write(f"{q[0]:.4f} {q[1]:.4f} {res['layer_of'][i]}\n")
    with open(os.path.join(plot_dir, "sites.dat"), "w") as f:
        f.write("# x(A) y(A) slot type\n")
        for s in all_sites:
            q = frame.xy(s["coord"])[0]
            f.write(f"{q[0]:.4f} {q[1]:.4f} {s['slot'] or 0} {s['type']}\n")
    with open(os.path.join(plot_dir, "minima_xy.dat"), "w") as f:
        f.write("# x(A) y(A) rank\n")
        for r, m in enumerate(res["minima"], 1):
            q = frame.xy(m["atoms"].positions[m["anchor"]])[0]
            f.write(f"{q[0]:.4f} {q[1]:.4f} {r}\n")
    out = _cell_outline(frame)
    with open(os.path.join(plot_dir, "cell.dat"), "w") as f:
        f.writelines(f"{p[0]:.4f} {p[1]:.4f}\n" for p in out)
    lo, hi = out.min(axis=0) - 1.0, out.max(axis=0) + 1.0
    L = [head.format(name="top_view.gplot", w=6, h=5, pdf="top_view.pdf"),
         f"set size ratio -1\nset xrange [{lo[0]:.2f}:{hi[0]:.2f}]\nset yrange [{lo[1]:.2f}:{hi[1]:.2f}]\n"
         'set xlabel "x (Ang)"\nset ylabel "y (Ang)"\nset title "Top view: substrate, sites, minima"\n'
         "set key outside right font \",7\"\n",
         'plot "substrate.dat" using 1:($3==2?$2:1/0) with points pt 7 ps 1.2 lc rgb "#777777" title "L2", \\\n'
         '     "substrate.dat" using 1:($3==1?$2:1/0) with points pt 7 ps 2.0 lc rgb "#222222" title "L1", \\\n'
         '     "cell.dat" using 1:2 with lines lc rgb "black" notitle, \\\n'
         '     "sites.dat" using 1:2:3 with points pt 6 ps 1.2 lc variable title "sites (colour: slot)", \\\n'
         '     "minima_xy.dat" using 1:2 with points pt 3 ps 2.5 lc rgb "#d4a017" title "minima", \\\n'
         '     "minima_xy.dat" using 1:2:3 with labels offset 1,1 notitle\n']
    with open(os.path.join(plot_dir, "top_view.gplot"), "w") as f:
        f.writelines(L)
    scripts.append("top_view.gplot")
    # minima
    write_minima_dat(os.path.join(plot_dir, "minima.dat"), res)
    parts = ['"minima.dat" using 1:5:xtic(sprintf("%s #%d %s", stringcolumn(2), $3, stringcolumn(4))) '
             f'with boxes fs solid 0.7 lc rgb "{_COLORS[0]}" title "{res["label"]}"']
    if found:
        write_minima_dat(os.path.join(plot_dir, "minima_foundation.dat"), found)
        parts.append(f'"minima_foundation.dat" using ($1+0.3):5 with boxes fs solid 0.5 lc rgb "{_COLORS[1]}" '
                     f'title "{found["label"]} (its own ranking)"')
    L = [head.format(name="site_energies.gplot", w=6, h=4, pdf="site_energies.pdf"),
         'set style fill solid\nset boxwidth 0.3\nset ylabel "E_{ads} (eV)"\nset title "Distinct minima"\n'
         "set xtics rotate by -30 font \",7\"\nset key bottom right font \",7\"\n",
         "plot " + ", \\\n     ".join(parts) + "\n"]
    with open(os.path.join(plot_dir, "site_energies.gplot"), "w") as f:
        f.writelines(L)
    scripts.append("site_energies.gplot")
    if res["sweep"]:
        write_sweep_dat(os.path.join(plot_dir, "binding_curve.dat"), res)
        blocks = sorted({(p["name"], p["slot"], p["type"]) for p in res["sweep"]})
        L = [head.format(name="binding_curve.gplot", w=6, h=4, pdf="binding_curve.pdf"),
             'set xlabel "anchor height above layer 1 (Ang)"\nset ylabel "E_{ads} (eV)"\n'
             'set title "Binding curve at fixed height"\nset key top right font ",7"\n',
             "plot " + ", \\\n     ".join(
                 f'"binding_curve.dat" index {i} using 1:($5==0?$2:1/0) with linespoints pt 7 '
                 f'lc rgb "{_slot_color(slot)}" title "{name} #{slot} {st}"' for i, (name, slot, st) in enumerate(blocks))
             + "\n"]
        with open(os.path.join(plot_dir, "binding_curve.gplot"), "w") as f:
            f.writelines(L)
        scripts.append("binding_curve.gplot")
    return scripts


def write_xyz(path, frames):
    """Extended XYZ: each frame with its energy, forces and fixed atoms; per-frame metadata in the
    comment line; per-atom 'adsorbate' (1) and 'layer' (substrate layer, 0 for the adsorbate)."""
    from ase.io import write as ase_write
    ase_write(path, frames, format="extxyz")


def _frame(at, n_sub, layer_of, **info):
    from ase.calculators.singlepoint import SinglePointCalculator
    b = at.copy()
    e, f = at.get_potential_energy(), at.get_forces(apply_constraint=False)
    # pymatgen's site properties come along as arrays (initial_magmoms = nan, surface_properties): drop them -- a nan
    # magmom also makes ASE treat the frame read back as changed, and refuse its stored forces
    for k in list(b.arrays):
        if k not in ("numbers", "positions"):
            del b.arrays[k]
    b.wrap()
    b.calc = SinglePointCalculator(b, energy=e, forces=f)
    b.arrays["adsorbate"] = np.array([0] * n_sub + [1] * (len(b) - n_sub), int)
    lay = np.zeros(len(b), int)
    lay[:n_sub] = layer_of[:n_sub]
    b.arrays["layer"] = lay
    b.info.update({k: v for k, v in info.items() if v is not None})
    return b


# ---- report pieces ----------------------------------------------------------------------------

def report_slab(res, args, f_out, warn, given_path):
    at = res["slab"]["atoms"]
    fr = res["frame"]
    sl = res["slab"]
    comp = {}
    for s in at.get_chemical_symbols():
        comp[s] = comp.get(s, 0) + 1
    thick = res["layers"][0]["h"] - res["layers"][-1]["h"]
    cn = abs(fr.c @ fr.n)
    print_dual(f"File              : {given_path}", f_out)
    print_dual(f"Atoms             : {len(at)} ({', '.join(f'{k}{v}' for k, v in comp.items())})", f_out)
    la, lb = np.linalg.norm(fr.a), np.linalg.norm(fr.b)
    gamma = np.degrees(np.arccos(fr.a @ fr.b / (la * lb)))
    print_dual(f"Surface cell      : a = {la:.4f}, b = {lb:.4f} Ang, gamma = {gamma:.2f} deg, area {fr.area:.3f} Ang^2 "
               f"({fr.area / len(res['layers'][0]['idx']):.3f} Ang^2 per top-layer atom)", f_out)
    print_dual(f"Thickness         : {thick:.3f} Ang in {len(res['layers'])} layers; vacuum {cn - thick:.2f} Ang "
               f"(c along the normal {cn:.3f} Ang)", f_out)
    rows = []
    free = set(sl["free"])
    hg = fr.height(sl["given"].positions)
    hg_top = hg[res["layers"][0]["idx"]].mean()
    for k, L in enumerate(res["layers"], 1):
        syms = {}
        for i in L["idx"]:
            syms[at[i].symbol] = syms.get(at[i].symbol, 0) + 1
        state = "free" if all(i in free for i in L["idx"]) else ("fixed" if not any(i in free for i in L["idx"])
                                                                   else "partly free")
        gap = f"{res['layers'][k - 2]['h'] - L['h']:.3f}" if k > 1 else "-"
        rows.append(([f"L{k}", f"{hg[L['idx']].mean() - hg_top:+.3f}", f"{L['h'] - res['h_top']:+.3f}", gap,
                      str(len(L["idx"])),
                      " ".join(f"{s}{n}" for s, n in syms.items()), state], None))
    print_table(["layer", "height as given (Ang)", "height relaxed (Ang)", "spacing (Ang)", "atoms", "elements",
                 "in the relaxations"], rows, f_out)
    print_dual(f"As given          : E = {sl['e_given']:.6f} eV, largest force {sl['f_given']:.4f} eV/Ang "
               f"(atom {sl['f_given_atom'] + 1}, L{res['layer_of'][sl['f_given_atom']]})", f_out)
    if sl["free"]:
        print_dual(f"Relaxed           : --relax-depth {args.relax_depth:g} Ang -> {len(sl['free'])} atom(s) free, "
                   f"{len(sl['fixed'])} fixed; E = {sl['e']:.6f} eV (relaxation {sl['e'] - sl['e_given']:+.4f} eV), "
                   f"max|F| on the free atoms {sl['f_given_free']:.4f} -> {sl['f_after']:.4f} eV/Ang, "
                   f"{sl['steps']} steps ({'converged' if sl['converged'] else 'HIT THE STEP CAP'}), "
                   f"largest move {sl['moved']:.3f} Ang", f_out)
        if not sl["converged"]:
            warn(f"Clean slab: its relaxation hit --max-steps {args.max_steps} (max|F| {sl['f_after']:.3f} eV/Ang); "
                 "every E_ads carries the residual -- raise --max-steps.")
        if not sl["fixed"]:
            warn("Clean slab: --relax-depth frees every atom of the slab (nothing is fixed); the slab can drift and "
                 "both surfaces relax. Use a depth smaller than the slab thickness.")
    else:
        print_dual("Relaxed           : no (--relax-depth 0): the substrate stays as given in every run; only the "
                   "adsorbate moves", f_out)
        top = res["layers"][0]["idx"]
        f_top = max(np.linalg.norm(sl["atoms"].get_forces(apply_constraint=False)[top], axis=1))
        if f_top > 0.1:
            warn(f"Clean slab: the model's force on the top layer as given is {f_top:.2f} eV/Ang -- the surface is not "
                 "at the model's minimum, and with --relax-depth 0 it cannot respond to the adsorbate: the E_ads mix "
                 "binding with a strained surface. Use --relax-depth (e.g. 3-5 Ang, the top layers).")


def report_references(res, args, f_out, warn, custom):
    rows = []
    for name, r in res["iso"].items():
        rows.append(([name, str(r["n"]), f"{r['extent']:.2f}", f"{args.vacuum_box:.1f}", f"{r['e']:.6f}",
                      f"{r['steps']}" + ("" if r["converged"] else " (cap)")], None))
        if r["extent"] > 0.5 * args.vacuum_box:
            warn(f"Reference: '{name}' spans {r['extent']:.2f} Ang, more than half of --vacuum-box "
                 f"({args.vacuum_box:.1f} Ang): it may see its periodic images. Raise --vacuum-box.")
    print_dual("Isolated adsorbate (alone in a cubic box, relaxed with the model):", f_out)
    print_table(["adsorbate", "atoms", "extent (Ang)", "box (Ang)", "E (eV)", "steps"], rows, f_out)
    if custom:
        print_dual("Note              : an atom or molecule in vacuum is far outside a model fine-tuned on condensed "
                   "SIESTA configurations -- its energy is little more than the model's E0s, so E_ads relative to it "
                   "carries that error. The energy differences between sites (dE) do not depend on it; "
                   "--bulk-reference gives a reference inside the training domain.", f_out)
    if res["bulk"]:
        rows = []
        for el, b in res["bulk"].items():
            rows.append(([el, b["file"], str(b["n"]), f"{b['v0']:.3f} -> {b['v']:.3f}", f"{b['mu']:.6f}",
                          f"{b['steps']}" + ("" if b["converged"] else " (cap)")], None))
            if abs(b["v"] / b["v0"] - 1) > 0.03:
                warn(f"Bulk reference {el}: the model's equilibrium volume differs by {100 * (b['v'] / b['v0'] - 1):+.1f} "
                     f"% from the file's ({b['file']}): check that the bulk is the phase you mean.")
        print_dual("Bulk reference (positions and cell relaxed with the model; mu = E/N):", f_out)
        print_table(["element", "file", "atoms", "volume/atom (Ang^3)", "mu (eV/atom)", "steps"], rows, f_out)
        print_dual("E_form = E(slab+adsorbate) - E(slab) - sum mu: the energy to take the adsorbate's atoms from their "
                   "bulk phases (an adatom formation energy).", f_out)


def report_sites(slots, all_sites, nops, args, f_out, warn):
    print_dual(f"Found             : {len(all_sites)} site(s) in the cell, {len(slots)} symmetry-distinct "
               f"(--symprec {args.symprec:g}; {nops} symmetry operation(s) keep the surface normal)", f_out)
    rows = []
    for s in slots:
        rows.append(([f"#{s['slot']}", s["type"], s["below"], str(s["mult"]), f"{s['uv'][0]:.4f} {s['uv'][1]:.4f}"], None))
    print_table(["site", "type", "what is below", "in the cell", "in-plane fractional (a, b)"], rows, f_out)
    unmapped = sum(1 for s in all_sites if s["slot"] is None)
    if unmapped:
        warn(f"Sites: {unmapped} site(s) of the cell could not be mapped to a symmetry-distinct one (the slab's symmetry "
             "was not found at --symprec); multiplicities and populations ignore them.")
    print_dual("Below             : the topmost atom within 0.5 Ang (in-plane) of the site, layers counted from the top "
               "(on fcc(111) a hollow above L2 is hcp and one above L3 is fcc; on Si(111) a hollow above L2 is T4 and "
               "one above L4 is H3).", f_out)


def report_configs(res, slots_by_id, args, f_out, warn):
    rows = []
    for c in sorted(res["configs"], key=lambda c: c["id"]):
        moved = c["end_slot"] != c["start_slot"]
        cells = [str(c["id"]), c["name"], f"#{c['start_slot']} {c['start_type']}", f"{c['h0']:.2f}",
                 site_name(c["end_slot"], slots_by_id, c["cls"]), f"{c['shift']:.2f}", f"{c['height']:.3f}",
                 f"{c['cn']}", f"{c['d_nearest']:.3f}", f"{c['e_ads']:.4f}"]
        if res["slab"]["free"]:
            cells.append(f"{c['sub_disp']:.3f}")
        cells += [f"{c['steps']}" + ("" if c["converged"] else " cap"), f"{c['fres']:.3f}"]
        rows.append((cells, "yellow" if not c["converged"] else ("cyan" if moved else None)))
    head = ["run", "ads", "start", "dist (A)", "ended at", "moved (A)", "height (A)", "CN", "d_nn (A)", "E_ads (eV)"]
    if res["slab"]["free"]:
        head.append("substrate (A)")
    head += ["steps", "max|F|"]
    print_table(head, rows, f_out)
    print_dual("Columns           : start = starting site and pymatgen distance; ended at = nearest site of the cell "
               f"to the adsorbate's lowest atom (within {SITE_MATCH} Ang, else off-site); moved = in-plane shift of that "
               "atom; height = above layer 1; CN = substrate neighbours within 1.25 x the covalent radii; substrate = "
               "largest substrate displacement from the clean slab. Cyan rows: the adsorbate left its starting site.",
               f_out)
    n_moved = sum(1 for c in res["configs"] if c["end_slot"] != c["start_slot"])
    if n_moved:
        print_dual(f"Moved             : {n_moved} of {len(res['configs'])} run(s) ended at another site -- the "
                   "starting site is not a minimum of the model there.", f_out)
    caps = [c["id"] for c in res["configs"] if not c["converged"]]
    if caps:
        warn(f"Relaxations: run(s) {', '.join(map(str, caps))} hit --max-steps {args.max_steps}; their energies are "
             "upper bounds. Raise --max-steps.")
    for c in res["configs"]:
        if c["d_nearest"] < 1.0:
            warn(f"Run {c['id']}: the adsorbate ended {c['d_nearest']:.2f} Ang from a substrate atom -- unphysical; the "
                 "model may be extrapolating.")


def report_minima(res, slots_by_id, args, f_out, warn):
    rows = []
    has_form = any(m["e_form"] is not None for m in res["minima"])
    for r, m in enumerate(res["minima"], 1):
        below = slots_by_id[m["end_slot"]]["below"] if m["end_slot"] is not None else "-"
        cells = [str(r), m["name"], site_name(m["end_slot"], slots_by_id, m["cls"]), below, f"{m['e_ads']:.4f}",
                 f"{m['de']:+.4f}"]
        if has_form:
            cells.append(f"{m['e_form']:.4f}" if m["e_form"] is not None else "-")
        cells += [str(m["mult"]), f"{100 * m['pop']:.1f}", f"{m['height']:.3f}", str(m["cn"]),
                  " ".join(f"{x:.3f}" for x in m["nn"][:3]) or "-", ",".join(map(str, m["starts"]))]
        rows.append((cells, "green" if m["de"] == 0 else None))
    head = ["rank", "ads", "site", "what is below", "E_ads (eV)", "dE (eV)"]
    if has_form:
        head.append("E_form (eV)")
    head += ["in cell", f"pop. {args.temperature:g} K (%)", "height (A)", "CN", "bonds (A)", "from runs"]
    print_table(head, rows, f_out)
    print_dual(f"Grouping          : runs of one adsorbate that end at symmetry-equivalent sites within {ENERGY_DUP} eV "
               "are one minimum; dE is relative to the most stable minimum of the same adsorbate; the population is "
               "multiplicity x exp(-dE/kT), normalized (configurational only: no vibrational entropy).", f_out)
    for name in dict.fromkeys(m["name"] for m in res["minima"]):
        ms = [m for m in res["minima"] if m["name"] == name]
        if len(ms) > 1:
            gap = ms[1]["e_ads"] - ms[0]["e_ads"]
            print_dual(f"{name:<18}: most stable {site_name(ms[0]['end_slot'], slots_by_id)} "
                       f"({slots_by_id[ms[0]['end_slot']]['below'] if ms[0]['end_slot'] else 'off-site'}), "
                       f"{gap:.4f} eV below the next minimum ({site_name(ms[1]['end_slot'], slots_by_id)})", f_out)
            if gap < 0.05:
                warn(f"Minima: for {name} the two most stable minima are only {1000 * gap:.0f} meV apart -- below what "
                     "an ML potential resolves reliably; check both with SIESTA.")
        else:
            print_dual(f"{name:<18}: a single minimum, {site_name(ms[0]['end_slot'], slots_by_id)}", f_out)
        missing = [s for s in slots_by_id.values() if s["slot"] not in {m["end_slot"] for m in ms}]
        if missing:
            miss = ", ".join(f"#{s['slot']} {s['type']}" for s in missing)
            print_dual(f"{'':<18}  no run ended at {miss}: not a minimum of the model (or not reached from these "
                       "starts)", f_out)


def report_sweep(res, f_out, warn):
    groups = {}
    for p in res["sweep"]:
        groups.setdefault((p["name"], p["slot"], p["type"]), []).append(p)
    rows = []
    for (name, slot, st), pts in sorted(groups.items()):
        pts.sort(key=lambda p: p["height"])
        bad = [p for p in pts if p["too_close"]]
        if bad:
            worst = min(bad, key=lambda p: p["d_nearest"])
            warn(f"Binding curve {name} #{slot}: at {len(bad)} height(s) the held adsorbate is closer than {TOO_CLOSE} x "
                 f"the covalent bond length to a substrate atom (down to {worst['d_nearest']:.2f} Ang, E_ads "
                 f"{worst['e_ads']:.1f} eV) -- those points are left out of the minimum: no model is trained there.")
        good = [p for p in pts if not p["too_close"]]
        if len(good) < 1:
            rows.append(([name, f"#{slot} {st}", "-", "-", "-", "-", "every point too close"], "yellow"))
            continue
        h = np.array([p["height"] for p in good])
        e = np.array([p["e_ads"] for p in good])
        k = int(np.argmin(e))
        h_min, e_min, note = h[k], e[k], "lowest point"
        if 0 < k < len(e) - 1:
            c = np.polyfit(h[k - 1:k + 2], e[k - 1:k + 2], 2)
            if c[0] > 0:
                h_min = -c[1] / (2 * c[0])
                e_min = np.polyval(c, h_min)
                note = "parabola through the 3 lowest points"
        elif k == 0:
            note = "at the LOWEST height scanned: extend the sweep down"
            warn(f"Binding curve {name} #{slot}: the minimum is at the lowest height scanned -- extend --height-sweep.")
        else:
            note = "at the HIGHEST height scanned: extend the sweep up"
        e_txt = " ".join(f"{p['e_ads']:.3f}" if not p["too_close"] else "(close)" for p in pts)
        rows.append(([name, f"#{slot} {st}", " ".join(f"{p['height']:.2f}" for p in pts), e_txt,
                      f"{h_min:.3f}", f"{e_min:.4f}", note], None))
    print_table(["ads", "site", "heights (Ang)", "E_ads (eV)", "h at min (Ang)", "E_ads at min (eV)", "how"], rows, f_out)
    print_dual("The anchor atom (the adsorbate's lowest) is held at each height, laterally too; the rest of the "
               "adsorbate and the free substrate relax. The height is measured from layer 1. '(close)': the anchor is "
               f"closer than {TOO_CLOSE} x the covalent bond length to a substrate atom -- left out.", f_out)


def report_committee(com, res, card, slots_by_id, f_out, warn):
    cal = ((card or {}).get("validation") or {}).get("calibration") or {}
    k = cal.get("error_over_spread")
    rows = []
    has_form = bool(np.isfinite(com["e_form"]).any())
    for j, m in enumerate(res["minima"]):
        e = com["e_ads"][:, j]
        first = [i for i, x in enumerate(res["minima"]) if x["name"] == m["name"]][0]
        de = com["e_ads"][:, j] - com["e_ads"][:, first]
        sp = 1000 * com["f_spread"][j]
        cells = [str(j + 1), m["name"], site_name(m["end_slot"], slots_by_id), f"{m['e_ads']:.4f}",
                 " ".join(f"{x:.3f}" for x in e), f"{e.std():.3f}", " ".join(f"{x:+.3f}" for x in de),
                 f"{de.std():.3f}"]
        if has_form:
            cells += [" ".join(f"{x:.3f}" for x in com["e_form"][:, j]), f"{com['form_std'][j]:.3f}"]
        cells += [f"{sp:.0f}", f"{k * sp:.0f}" if k else "-"]
        rows.append((cells, None))
    print_dual(f"{com['n']} committee members on the main model's geometries -- E_ads against each member's own slab and "
               "isolated adsorbate (it carries each member's E0s), dE against each member's own minimum 1"
               + (", E_form against each member's own bulk" if has_form else "") + ":", f_out)
    head = ["rank", "ads", "site", "main E_ads", "members E_ads (eV)", "std", "members dE (eV)", "std"]
    if has_form:
        head += ["members E_form (eV)", "std"]
    head += ["F spread (meV/Ang)", "est. F error"]
    print_table(head, rows, f_out)
    for name in dict.fromkeys(m["name"] for m in res["minima"]):
        idx = [j for j, m in enumerate(res["minima"]) if m["name"] == name]
        if len(idx) < 2:
            continue
        best = [idx[int(np.argmin(com["e_ads"][i, idx]))] for i in range(com["n"])]
        agree = sum(1 for b in best if b == idx[0])
        de = com["e_ads"][:, idx[1]] - com["e_ads"][:, idx[0]]
        print_dual(f"{name:<18}: {agree} of {com['n']} members put the same minimum first; gap to the second "
                   f"{', '.join(f'{x:.3f}' for x in de)} eV (std {de.std():.3f} eV)", f_out)
        if agree < com["n"]:
            warn(f"Committee: for {name} not every member ranks the same minimum first -- the site preference is not "
                 "determined by the training data; check with SIESTA.")
        elif de.std() > 0.5 * abs(de.mean()) and abs(de.mean()) > 0:
            warn(f"Committee: for {name} the gap between the two best minima varies by {de.std():.3f} eV across members, "
                 "comparable with the gap itself.")
    ref = cal.get("spread_mean_meV")
    if k:
        print_dual(f"Calibration       : on the card's test set the true force error is ~{k:.1f} x the spread"
                   + (f"; the mean spread there is {ref:.1f} meV/Ang" if ref else ""), f_out)
    top = 1000 * com["f_spread"].max()
    if ref and top > 3 * ref:
        j = int(np.argmax(com["f_spread"]))
        warn(f"Committee: the force spread reaches {top:.0f} meV/Ang on minimum {j + 1} ({top / ref:.0f}x the test "
             "set's mean): the model is far outside what its training set constrains there -- the site energies and "
             "the ranking need a SIESTA check (and the structures belong in the training set).")
    if com["e_std"].max() > 0.1:
        warn(f"Committee: the members' E_ads differ by up to {com['e_std'].max():.2f} eV (std) on one minimum.")
    print_dual("Seeds fine-tuned from the same foundation model are nearly the same model: the spread is a LOWER bound "
               "of the error (7.5: a committee agreed within 0.09 eV on a barrier 0.47 eV too high).", f_out)


def report_foundation(res, found, slots_by_id, f_out, warn):
    by_id = {c["id"]: c for c in found["configs"]}
    rows = []
    has_form = any(c["e_form"] is not None for c in res["configs"])
    best_a = {n: min(c["e_ads"] for c in res["configs"] if c["name"] == n) for n in {c["name"] for c in res["configs"]}}
    best_b = {n: min(c["e_ads"] for c in found["configs"] if c["name"] == n) for n in {c["name"] for c in found["configs"]}}
    for c in sorted(res["configs"], key=lambda c: c["id"]):
        f = by_id.get(c["id"])
        if f is None:
            continue
        cells = [str(c["id"]), c["name"], f"#{c['start_slot']} {c['start_type']}",
                 site_name(c["end_slot"], slots_by_id), f"{c['e_ads'] - best_a[c['name']]:+.4f}"]
        if has_form:
            cells.append(f"{c['e_form']:.4f}")
        cells += [site_name(f["end_slot"], slots_by_id), f"{f['e_ads'] - best_b[f['name']]:+.4f}"]
        if has_form:
            cells += [f"{f['e_form']:.4f}", f"{f['e_form'] - c['e_form']:+.4f}"]
        rows.append((cells, None if c["end_slot"] == f["end_slot"] else "cyan"))
    print_dual(f"The same runs with {found['label']} (its own clean slab, references and relaxations); dE against each "
               "model's own most stable run" + (", E_form against each model's own bulk" if has_form else "") + ":",
               f_out)
    head = ["run", "ads", "start", f"{res['label']}: end", "dE (eV)"] + (["E_form (eV)"] if has_form else []) + \
           [f"{found['label']}: end", "dE (eV)"] + (["E_form (eV)", "difference"] if has_form else [])
    print_table(head, rows, f_out)
    for name in dict.fromkeys(m["name"] for m in res["minima"]):
        a = [m for m in res["minima"] if m["name"] == name]
        b = [m for m in found["minima"] if m["name"] == name]
        if not (a and b):
            continue
        same = a[0]["end_slot"] == b[0]["end_slot"]
        print_dual(f"{name:<18}: most stable {site_name(a[0]['end_slot'], slots_by_id)} ({res['label']}) vs "
                   f"{site_name(b[0]['end_slot'], slots_by_id)} ({found['label']}) -- "
                   f"{'same site' if same else 'DIFFERENT site'}", f_out)
        if len(a) > 1 and len(b) > 1:
            print_dual(f"{'':<18}  gap to the second minimum {a[1]['e_ads'] - a[0]['e_ads']:.4f} vs "
                       f"{b[1]['e_ads'] - b[0]['e_ads']:.4f} eV", f_out)
        if not same:
            warn(f"Foundation comparison: for {name} the fine-tuning changed which site is most stable -- confirm with "
                 "SIESTA before trusting either.")
    print_dual("E_ads of two models are not compared: each is measured from the model's own isolated adsorbate, i.e. "
               "its own E0s; compare dE and, with --bulk-reference, E_form.", f_out)


def report_neighbours(res, slots_by_id, args, f_out):
    for name in dict.fromkeys(m["name"] for m in res["minima"]):
        nb = neighbours(res, name)
        if not nb:
            continue
        a = [m for m in res["minima"] if m["name"] == name][0]
        a_file = f"rank01_{a['end_type']}_site{a['end_slot'] or 0}_{name}.fdf"
        rows = []
        for n in nb:
            m = n["min"]
            rows.append(([str(n["rank"]), site_name(m["end_slot"], slots_by_id),
                          slots_by_id[m["end_slot"]]["below"] if m["end_slot"] is not None else "-",
                          f"{m['de']:+.4f}", f"{n['d_site']:.3f}" if n["d_site"] is not None else "-",
                          f"{n['d_file']:.3f}"], None))
        print_dual(f"From minimum 1 of {name}, {site_name(a['end_slot'], slots_by_id)}, to every other minimum "
                   "(in-plane distances):", f_out)
        print_table(["rank", "site", "what is below", "dE (eV)", "nearest copy of the site (Ang)",
                     "in its .fdf (Ang)"], rows, f_out)
        near = min(nb, key=lambda n: n["d_file"])
        b_file = (f"rank{near['rank']:02d}_{near['min']['end_type']}_site{near['min']['end_slot'] or 0}_{name}.fdf")
        print_dual("Diffusion barrier : a climbing-image NEB with stb-mlneb (7.5) between two of these minima -- with "
                   "its committee, SIESTA check (--siesta-check) and export to training (--to-mlff):", f_out)
        print_dual(f"  stb-mlneb --initial {args.output_dir}/{a_file} --final {args.output_dir}/{b_file} --autosort-tol 0 "
                   + (f"--custom-model {args.custom_model} --custom-model-card {args.custom_model_card}"
                      if args.custom_model else f"--model {args.model}"), f_out)
        print_dual(f"  (--autosort-tol 0: the .fdf files list the same atoms in the same order, and an adsorbate that hops "
                   "farther than the default 0.5 Ang cannot be matched by distance. The minimum nearest to minimum 1 as "
                   f"written in its .fdf, {near['d_file']:.2f} Ang. A minimum whose "
                   "nearest site copy is closer than its .fdf sits at another copy; such a short hop is real only if "
                   "the adsorbate stays there when relaxed next to minimum 1 -- not checked here)", f_out)


def report_domain(card, res, f_out, warn):
    dom = (card or {}).get("domain") or {}
    rows = []
    trained = set(dom.get("elements") or card.get("elements") or [])
    present = sorted({s for c in res["configs"] for s in c["atoms"].get_chemical_symbols()})
    miss = [s for s in present if trained and s not in trained]
    rows.append((["elements", ", ".join(sorted(trained)) or "-", ", ".join(present),
                  "OUTSIDE" if miss else "inside"], "yellow" if miss else None))
    rows.append((["dimensionality", str(dom.get("dimensionality") or "-"), "2D slab; 0D isolated adsorbate",
                  "OUTSIDE" if dom.get("dimensionality") not in (None, "2D") else "inside"],
                 "yellow" if dom.get("dimensionality") not in (None, "2D") else None))
    dd = dom.get("d_min_ang")
    if dd:
        dmin = min(c["dmin"] for c in res["configs"] if c["dmin"])
        ok = dmin >= 0.95 * dd["min"]
        rows.append((["shortest distance (Ang)", f">= {dd['min']:.3f}", f"{dmin:.3f}", "inside" if ok else "OUTSIDE"],
                     None if ok else "yellow"))
        if not ok:
            warn(f"Model domain: a relaxed structure has two atoms {dmin:.3f} Ang apart, below the shortest trained "
                 f"distance ({dd['min']:.3f} Ang).")
    fm = dom.get("force_max_eV_per_ang")
    if fm:
        f0 = res["slab"]["f_given"]
        rows.append((["largest force (eV/Ang)", f"<= {fm['max']:.2f}", f"{f0:.3f} (slab as given)",
                      "inside" if f0 <= fm["max"] else "OUTSIDE"], None if f0 <= fm["max"] else "yellow"))
    er = dom.get("energy_per_atom_eV")
    if er:
        epa = [c["e"] / len(c["atoms"]) for c in res["configs"]]
        lo = min(epa) < er["min"] - 0.02
        rows.append((["energy per atom (eV)", f"{er['min']:.4f} .. {er['max']:.4f}", f"{min(epa):.4f} .. {max(epa):.4f}",
                      "BELOW the lowest trained" if lo else ("above (a surface costs energy)" if min(epa) > er["max"]
                                                              else "inside")], "yellow" if lo else None))
        if lo:
            warn("Model domain: a relaxed structure is lower in energy per atom than anything in the training set -- a "
                 "hole in the learned surface; check with SIESTA.")
    print_table(["quantity", "trained", "here", "verdict"], rows, f_out)
    kinds = ((card.get("training_data") or {}).get("kinds") or {})
    if kinds:
        print_dual(f"Training kinds    : {', '.join(f'{k} {v}' for k, v in kinds.items())} -- a surface, an adsorbate "
                   "and its bonds to the surface are only as good as the training data near them; the committee "
                   "spread and a SIESTA check of the best sites are the way to find out", f_out)


# ---- main ---------------------------------------------------------------------------------------

def build_parser():
    p = argparse.ArgumentParser(
        description="Adsorption-site screening driven entirely by a MACE potential -- the MACE-MP-0 foundation model, or "
                    "a model fine-tuned on your own SIESTA data with the 7.1 workflow (--custom-model, always with its "
                    "model card). Finds the candidate ontop/bridge/hollow sites of the slab (symmetry-reduced, each "
                    "described by the atom below it), relaxes the adsorbate on each (and the top layers of the slab with "
                    "--relax-depth), reports where each run ENDED (a bridge start that slides into a hollow is a hollow), "
                    "groups the runs into distinct minima with their multiplicity and Boltzmann population, and the "
                    "adsorption energy relative to the isolated adsorbate and, with --bulk-reference, to the elements' "
                    "bulk, and how far the other minima are from the most stable one (the input of a diffusion "
                    "barrier with stb-mlneb, 7.5). Optional: binding curves at fixed height, the committee's spread on "
                    "every minimum and the same screening with the foundation model. Every structure is checked against the model card's trained range; all "
                    "warnings in the last section. No SIESTA input is written; for the DFT pipeline use stb-adsorb (5.8).",
        epilog="Example usage:\n"
               "  stb-mladsorb --slab slab.fdf --adsorbate H\n"
               "  stb-mladsorb --slab slab.fdf --adsorbate O,H --site-type hollow --relax-depth 3\n"
               "  stb-mladsorb --slab si111.fdf --adsorbate Si --relax-depth 4 --bulk-reference si_bulk.fdf "
               "--custom-model best_model/mlff_model_compiled.model --custom-model-card "
               "best_model/mlff_model_card.json --device cuda --show --save-data --save-gnuplot --save-xyz\n"
               "  stb-mladsorb --slab slab.fdf --adsorbate H2O --n-orientations 8 --seed 1\n",
        formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("--slab", required=True, help="Slab/2D-material structure (.fdf), vacuum along one axis.")
    p.add_argument("--adsorbate", required=True,
                   help="Adsorbate(s): an ASE G2 molecule name (e.g. 'H2O') or an element symbol, comma-separated for "
                        "several (e.g. 'H,O').")
    p.add_argument("--model", choices=["small", "medium", "large"], default="medium",
                   help="MACE-MP-0 model size (default: medium). With --custom-model, the foundation model of the "
                        "comparison.")
    p.add_argument("--custom-model", default=None, metavar="PATH",
                   help="A custom MACE model file instead of the MACE-MP-0 foundation potential -- e.g. one fine-tuned "
                        "via stb-mlffAnalysis. Overrides --model.")
    model_card.add_card_argument(p, "--custom-model")
    p.add_argument("--skip-foundation-comparison", action="store_true",
                   help="With --custom-model, do not repeat the screening with the foundation model (--model). The "
                        "comparison is on by default: it shows whether the fine-tuning changed the site ranking.")
    p.add_argument("--mlff-workflow", default=None, metavar="JSON",
                   help="The mlff_workflow.json of the 7.1 run that trained --custom-model: locates its committee.")
    p.add_argument("--no-committee", dest="committee", action="store_false",
                   help="Do not evaluate the committee of the custom model on the minima. On by default when the card "
                        "records a committee whose model files are found (next to the card's run folder, or in "
                        "--mlff-workflow's).")
    p.add_argument("--committee-models", nargs="+", default=None, metavar="PATH",
                   help="The committee's model files explicitly (each with its mlff_model_card.json next to it).")
    p.add_argument("--site-type", choices=["all", "ontop", "bridge", "hollow"], default="all",
                   help="Which candidate site type(s) to screen (default: all).")
    p.add_argument("--height", type=float, default=2.0,
                   help="Starting distance of the adsorbate from the surface, Ang -- pymatgen's AdsorbateSiteFinder "
                        "definition (default: 2.0).")
    p.add_argument("--start-heights", type=float, nargs="+", default=None, metavar="H",
                   help="Several starting distances (each one a run per site), e.g. 1.5 2.0 2.5: a relaxation can end "
                        "in a different minimum depending on where it starts. Default: --height only.")
    p.add_argument("--height-sweep", type=float, nargs=3, default=None, metavar=("MIN", "MAX", "STEP"),
                   help="Binding curve per site: the adsorbate's lowest atom held at each height (the rest relaxes), "
                        "E_ads(h) and its minimum. On top of the free relaxations.")
    p.add_argument("--relax-depth", type=float, default=0.0, metavar="ANG",
                   help="Substrate atoms within this depth (Ang) of the top surface relax -- in the clean slab first, "
                        "then with the adsorbate; deeper atoms stay fixed (default: 0 = the substrate is fixed as "
                        "given and only the adsorbate moves). A slab cut from the bulk is not at the model's minimum: "
                        "use the top layers (3-5 Ang).")
    p.add_argument("--bulk-reference", nargs="+", default=None, metavar="FDF",
                   help="Bulk structure(s) of the adsorbate's element(s), one element per file: mu = E/N after relaxing "
                        "positions and cell with the model, and E_form = E(slab+ads) - E(slab) - sum mu is reported "
                        "next to E_ads.")
    p.add_argument("--temperature", type=float, default=300.0,
                   help="Temperature of the site populations, K (default: 300).")
    p.add_argument("--n-orientations", type=int, default=1,
                   help="For a multi-atom adsorbate, screen this many random orientations per site by a single-point "
                        "energy and relax the best (default: 1 = the G2 orientation). No effect on a single atom.")
    p.add_argument("--seed", type=int, default=None, help="Random seed of --n-orientations (default: unseeded).")
    p.add_argument("--symprec", type=float, default=0.01,
                   help="Symmetry tolerance of the site reduction (default: 0.01).")
    p.add_argument("--vacuum-gap", type=float, default=10.0,
                   help="Gap (Ang) that identifies the slab's vacuum axis (default: 10.0).")
    p.add_argument("--vacuum-box", type=float, default=20.0,
                   help="Cubic box (Ang) of the isolated-adsorbate reference (default: 20.0).")
    p.add_argument("--fmax", type=float, default=0.05,
                   help="Force convergence of every relaxation and of the NEB, eV/Ang (default: 0.05).")
    p.add_argument("--max-steps", type=int, default=300,
                   help="Optimizer steps per relaxation (default: 300).")
    p.add_argument("--top-k", type=int, default=None,
                   help="Write only the N most stable minima as .fdf (default: all minima).")
    p.add_argument("--device", choices=["cpu", "cuda"], default="cpu", help="Device of the model (default: cpu).")
    p.add_argument("-o", "--output-dir", default="mladsorb_out", help="Output directory (default: mladsorb_out).")
    p.add_argument("--show", action="store_true",
                   help=f"Open the figure on screen at the end (needs a display). {FIGURE_FILE} is written either way.")
    p.add_argument("--save-data", action="store_true",
                   help="Write the tables as .dat: every relaxation (mladsorb_configs.dat), the minima "
                        "(mladsorb_minima.dat, with the committee spread) and the binding curves. Off by default.")
    p.add_argument("--save-gnuplot", action="store_true",
                   help=f"Write .dat + .gplot into <output-dir>/{PLOT_DIR}/ (top view, site energies, binding curve) "
                        "and run gnuplot there when it is installed. Off by default.")
    p.add_argument("--save-xyz", action="store_true",
                   help="Write extended XYZ files: every relaxed run (mladsorb_runs.xyz), the minima "
                        "(mladsorb_minima.xyz), the clean slab and the binding curves -- energy, forces, fixed atoms, "
                        "per-atom 'adsorbate' and 'layer', and the run's data in each frame's comment line. Off by "
                        "default.")
    p.add_argument("--save-report", action="store_true", help=f"Also write the report to {REPORT_FILE}.")
    p.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    p.add_argument("-v", "--version", action="version", version=f"stb-mladsorb {VERSION}")
    return p


def main():
    parser = build_parser()
    args = parser.parse_args()
    card = model_card.setup_from_args(args.custom_model, args.custom_model_card, "--custom-model", error=parser.error)
    if args.relax_depth < 0:
        parser.error("--relax-depth must be >= 0.")
    if args.mlff_workflow and not os.path.isfile(args.mlff_workflow):
        parser.error(f"--mlff-workflow file '{args.mlff_workflow}' not found.")

    import matplotlib
    if not args.show:
        matplotlib.use("Agg")

    if args.intro:
        show_intro(["Siesta ToolBox Suite - ML Adsorption Screening",
                    "Adsorption sites and minima with a MACE potential",
                    f"Version {VERSION} | University of Brasilia - 2026",
                    "Developed by Dr. Carlos M. O. Bastos"])

    os.makedirs(args.output_dir, exist_ok=True)
    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None
    warns, noise, files = [], [], []
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
        from stb.adsorb import resolve_slab_orientation, parse_adsorbates
        import torch
        import ase
        import mace
    model_card.collect_warnings(warns)

    print("\n" + color_text("ML ADSORPTION-SITE SCREENING:", 'bold'))
    print("-" * 60)
    print_dual(color_text("===== STB-MLADSORB REPORT =====", 'magenta'), f_out)

    if not os.path.isfile(args.slab):
        fail(f"Slab file '{args.slab}' not found.")
    for b in args.bulk_reference or []:
        if not os.path.isfile(b):
            fail(f"--bulk-reference file '{b}' not found.")
    adsorbates = parse_adsorbates(args.adsorbate)
    names = [n for n, _ in adsorbates]
    compare = bool(args.custom_model) and not args.skip_foundation_comparison
    start_heights = args.start_heights or [args.height]
    site_types = ["ontop", "bridge", "hollow"] if args.site_type == "all" else [args.site_type]

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
    print_dual(f"Slab              : {args.slab}", f_out)
    print_dual(f"Adsorbate(s)      : {', '.join(names)}", f_out)
    print_dual(f"Sites             : {', '.join(site_types)}; start at {', '.join(f'{h:g}' for h in start_heights)} Ang"
               + (f"; {args.n_orientations} orientations each" if args.n_orientations > 1 else ""), f_out)
    print_dual("Substrate         : " + (f"top {args.relax_depth:g} Ang relaxed, deeper atoms fixed" if args.relax_depth
                                         else "fixed as given (only the adsorbate relaxes)"), f_out)
    print_dual(f"Relaxations       : FIRE to max|F| < {args.fmax:g} eV/Ang, at most {args.max_steps} steps", f_out)
    extras = []
    if args.height_sweep:
        extras.append("binding curves {:g}-{:g} Ang step {:g}".format(*args.height_sweep))
    if args.bulk_reference:
        extras.append("bulk reference")
    if compare:
        extras.append(f"foundation comparison ({args.model})")
    print_dual(f"Also              : {', '.join(extras) if extras else '-'}", f_out)
    print_dual(f"Temperature       : {args.temperature:g} K (site populations)", f_out)

    sections("MODEL")
    model_arg = args.custom_model if args.custom_model else args.model
    t0 = time.monotonic()
    if card is not None:
        model_card.mark_announced(args.custom_model)
    with capture_library_noise(noise, "model"):
        calc = mace_relax.get_calculator(model=model_arg, device=args.device)
    if card is not None:
        model_card.print_model_card(card, args.custom_model_card, args.custom_model, f_out, warn)
    else:
        print_dual(f"Model             : MACE-MP-0 {args.model} (foundation, no card): trained on the Materials Project "
                   "PBE/PBE+U relaxations -- not on your SIESTA level of theory", f_out)
    print_dual(f"Loaded in         : {time.monotonic() - t0:.1f} s", f_out)
    label_main = "fine-tuned" if args.custom_model else f"MACE-MP-0 {args.model}"

    try:
        structure = structure_io.read_fdf(args.slab)
        structure, relabel_note = resolve_slab_orientation(structure, args.vacuum_gap)
    except Exception as e:
        fail(f"Could not use '{args.slab}' as a slab: {e}")
    if relabel_note:
        warn(f"Slab: {relabel_note}")
    slab_pmg = structure_io.to_pymatgen(structure)
    species_meta = structure_io.species_dict(structure)
    slab_given = AseAtomsAdaptor.get_atoms(slab_pmg)
    slab_given.set_constraint()
    ctx = {"sites": None, "site_types": site_types, "start_heights": start_heights}
    rng = np.random.default_rng(args.seed)

    def factory(model):
        def make():
            if os.path.isfile(str(model)):
                model_card.mark_announced(model)
            return mace_relax.get_calculator(model=model, device=args.device)
        return make

    print(f"  ... {label_main}: screening", flush=True)
    timing = []
    res = screen(label_main, calc, factory(model_arg), slab_given, ctx, args, adsorbates, noise, rng)
    timing.append((f"{label_main} screening", res["seconds"]))
    slots, all_sites, nops = ctx["sites"]
    res["all_sites"] = all_sites
    res["slots"] = slots
    slots_by_id = {s["slot"]: s for s in slots}

    sections("CLEAN SLAB")
    report_slab(res, args, f_out, warn, args.slab)

    sections("REFERENCES")
    report_references(res, args, f_out, warn, bool(args.custom_model))

    sections("CANDIDATE SITES")
    if not slots:
        fail(f"No {'/'.join(site_types)} sites found on this slab.")
    report_sites(slots, all_sites, nops, args, f_out, warn)

    sections("RELAXATIONS")
    print_dual(f"{len(res['configs'])} relaxation(s) with {label_main} in {res['seconds']:.1f} s:", f_out)
    report_configs(res, slots_by_id, args, f_out, warn)

    sections("DISTINCT MINIMA")
    report_minima(res, slots_by_id, args, f_out, warn)
    if any(len([m for m in res["minima"] if m["name"] == n]) > 1 for n in names):
        print_dual("", f_out)
        report_neighbours(res, slots_by_id, args, f_out)
    print_dual(color_text("[NOTE] ML-level energies from a fast surrogate potential: use them to choose which site(s) "
                          "to compute with SIESTA (stb-adsorb, 5.8), not as the final answer.", 'yellow'), f_out)

    if res["sweep"]:
        sections("BINDING CURVES (FIXED HEIGHT)")
        report_sweep(res, f_out, warn)

    com = None
    if args.custom_model and args.committee and res["minima"]:
        from stb.mlneb import committee_members
        members, note = committee_members(card, args.custom_model_card, args.mlff_workflow, args.committee_models)
        sections("COMMITTEE")
        if members:
            print_dual(f"Members           : {note}", f_out)
            print(f"  ... committee: {len(members)} members on {len(res['minima'])} minima", flush=True)
            t_c = time.monotonic()
            com = committee_minima(res, members, args.device, noise, adsorbates)
            timing.append(("committee", time.monotonic() - t_c))
            report_committee(com, res, card, slots_by_id, f_out, warn)
        else:
            print_dual(f"Not evaluated     : {note}", f_out)

    found = None
    if compare:
        sections("FOUNDATION MODEL COMPARISON")
        with capture_library_noise(noise, "foundation model"):
            found_calc = mace_relax.get_calculator(model=args.model, device=args.device)
        found = screen(f"MACE-MP-0 {args.model}", found_calc, factory(args.model), slab_given, ctx, args, adsorbates,
                       noise, np.random.default_rng(args.seed))
        found["all_sites"] = all_sites
        found["slots"] = slots
        timing.append((f"{found['label']} screening", found["seconds"]))
        report_foundation(res, found, slots_by_id, f_out, warn)

    if card is not None:
        sections("TRAINED RANGE")
        report_domain(card, res, f_out, warn)

    # ---- files
    fig_path = os.path.join(args.output_dir, FIGURE_FILE)
    fig = make_figure(res, found, com, slots, all_sites, fig_path)
    files.append(("Figure", fig_path))
    write_list = res["minima"][:args.top_k] if args.top_k is not None else res["minima"]
    for r, m in enumerate(write_list, 1):
        wrapped = m["atoms"].copy()
        wrapped.wrap()
        st = structure_io.from_pymatgen(AseAtomsAdaptor.get_structure(wrapped), species_meta=species_meta,
                                        coord_format="fractional")
        p = os.path.join(args.output_dir, f"rank{r:02d}_{m['end_type']}_site{m['end_slot'] or 0}_{m['name']}.fdf")
        structure_io.write_fdf(st, p)
    if write_list:
        files.append(("Minima (.fdf)", f"{args.output_dir}/rank01_*.fdf .. rank{len(write_list):02d}_*.fdf"))
    if args.save_data:
        p = os.path.join(args.output_dir, "mladsorb_configs.dat")
        write_configs_dat(p, res)
        files.append(("Data", p))
        p = os.path.join(args.output_dir, "mladsorb_minima.dat")
        write_minima_dat(p, res, com)
        files.append(("Data", p))
        if res["sweep"]:
            p = os.path.join(args.output_dir, "mladsorb_binding_curve.dat")
            write_sweep_dat(p, res)
            files.append(("Data", p))
        if found is not None:
            p = os.path.join(args.output_dir, "mladsorb_configs_foundation.dat")
            write_configs_dat(p, found)
            files.append(("Data", p))
    if args.save_xyz:
        n_sub = len(res["slab"]["atoms"])
        lay = res["layer_of"]
        frames = [_frame(c["atoms"], n_sub, lay, run=c["id"], adsorbate_name=c["name"], start_site=c["start_slot"],
                         start_type=c["start_type"], start_distance=c["h0"], end_site=c["end_slot"] or 0,
                         end_type=c["end_type"], E_ads=c["e_ads"], E_form=c["e_form"], height=c["height"],
                         CN=c["cn"], converged=c["converged"], steps=c["steps"])
                  for c in sorted(res["configs"], key=lambda c: c["id"])]
        p = os.path.join(args.output_dir, "mladsorb_runs.xyz")
        write_xyz(p, frames)
        files.append(("Extended XYZ", f"{p} ({len(frames)} relaxed runs)"))
        frames = [_frame(m["atoms"], n_sub, lay, rank=r, adsorbate_name=m["name"], end_site=m["end_slot"] or 0,
                         end_type=m["end_type"], E_ads=m["e_ads"], dE=m["de"], E_form=m["e_form"],
                         multiplicity=m["mult"], population=m["pop"],
                         committee_std=float(com["e_std"][r - 1]) if com is not None else None)
                  for r, m in enumerate(res["minima"], 1)]
        p = os.path.join(args.output_dir, "mladsorb_minima.xyz")
        write_xyz(p, frames)
        files.append(("Extended XYZ", f"{p} ({len(frames)} minima)"))
        p = os.path.join(args.output_dir, "mladsorb_slab.xyz")
        write_xyz(p, [_frame(res["slab"]["atoms"], n_sub, lay, relax_depth=args.relax_depth)])
        files.append(("Extended XYZ", f"{p} (clean slab)"))
        if res["sweep"]:
            frames = [_frame(q["atoms"], n_sub, lay, adsorbate_name=q["name"], site=q["slot"], height=q["height"],
                             E_ads=q["e_ads"]) for q in res["sweep"]]
            p = os.path.join(args.output_dir, "mladsorb_binding_curve.xyz")
            write_xyz(p, frames)
            files.append(("Extended XYZ", f"{p} ({len(frames)} constrained runs)"))
    if args.save_gnuplot:
        plot_dir = os.path.join(args.output_dir, PLOT_DIR)
        scripts = gnuplot_files(plot_dir, res, found, slots, all_sites)
        pdfs, msg = run_gnuplot(plot_dir, scripts)
        files.append(("gnuplot folder", f"{plot_dir}/ ({', '.join(scripts)})"))
        for p in pdfs:
            files.append(("gnuplot figure", p))
        if msg:
            warn(msg)

    sections("SUMMARY & FILES")
    print_dual("Status            : OK" + (f" ({len(warns)} warning(s), see the last section)" if warns else ""), f_out)
    for name in names:
        ms = [m for m in res["minima"] if m["name"] == name]
        if not ms:
            continue
        m = ms[0]
        below = slots_by_id[m["end_slot"]]["below"] if m["end_slot"] else "off-site"
        txt = (f"{name:<18}: most stable {site_name(m['end_slot'], slots_by_id)} ({below}), E_ads {m['e_ads']:.4f} eV"
               + (f", E_form {m['e_form']:.4f} eV" if m["e_form"] is not None else ""))
        if len(ms) > 1:
            txt += f"; next {site_name(ms[1]['end_slot'], slots_by_id)} +{ms[1]['de']:.4f} eV"
        print_dual(txt, f_out)
    if com is not None:
        print_dual(f"Committee         : E_ads std up to {com['e_std'].max():.4f} eV over the minima", f_out)
    for kind, path in files:
        print_dual(f"{kind:<18}: {path}", f_out)
    if report_path:
        print_dual(f"Report            : {report_path}", f_out)
    not_asked = [f for f, on in (("--show", args.show), ("--save-data", args.save_data),
                                 ("--save-gnuplot", args.save_gnuplot), ("--save-xyz", args.save_xyz),
                                 ("--save-report", args.save_report),
                                 ("--relax-depth", args.relax_depth > 0),
                                 ("--bulk-reference", bool(args.bulk_reference))) if not on]
    if not_asked:
        print_dual(f"Also available    : {' '.join(not_asked)}", f_out)
    print_dual(f"Wall time         : {time.monotonic() - t_start:.1f} s ("
               + ", ".join(f"{k} {v:.0f} s" for k, v in timing) + ")", f_out)
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
    print(color_text("ML adsorption screening complete -- a fast estimate, not a substitute for DFT.\n", 'bold'))
    import matplotlib.pyplot as plt
    if args.show and display_available():
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    main()
