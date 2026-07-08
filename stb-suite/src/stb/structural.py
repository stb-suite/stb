#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.9.1"

import os
import sys
import argparse
import logging
import warnings
from datetime import datetime
import numpy as np
from scipy.spatial import ConvexHull, QhullError
from pymatgen.analysis.local_env import (
    MinimumDistanceNN, CrystalNN, BrunnerNNRelative, EconNN)
from ase.io import read as ase_read
from pymatgen.io.ase import AseAtomsAdaptor
from stb.core import structure_io
from stb.core.cli import color_text, show_intro

def warn_handler(message, category, filename, lineno, file=None, line=None):
    log_message = f"{category.__name__}: {message} (File: {filename}, Line: {lineno})"
    logging.warning(log_message)
    print("[WARNING] Warning detected! Check warnings.log for details.")

def read_structure(path, fmt):
    """Returns a pymatgen Structure for one of this suite's own SIESTA
    formats: .fdf (structure input, via the shared core/structure_io.py
    parser) or .STRUCT_OUT (post-relaxation output, via ASE -- there's no
    shared reader for that format elsewhere in the suite)."""
    if fmt == "fdf":
        fdf_structure = structure_io.read_fdf(path)
        return structure_io.to_pymatgen(fdf_structure)
    atoms = ase_read(path, format="struct_out")
    return AseAtomsAdaptor.get_structure(atoms)

def _safe_mean(values):
    values = [v for v in values if v is not None]
    return float(np.mean(values)) if values else None

def _bond_length_distortion(distances):
    """Baur's bond-length distortion index (%): mean absolute deviation of
    each bond length from this site's own mean bond length. Needs at least
    2 bonds; returns None otherwise (distortion is undefined for a single
    bond)."""
    if len(distances) < 2:
        return None
    d_avg = np.mean(distances)
    if d_avg == 0:
        return None
    return float(100.0 * np.mean([abs(d - d_avg) / d_avg for d in distances]))

def _bond_angle_variance(angles):
    """Bond-angle variance (deg²): variance of this site's L-center-L
    angles around their own mean. Two deliberate departures from Robinson
    et al.'s original (octahedron/tetrahedron) formula, both because the
    coordination number here can be any -- often fractional, "effective"
    -- value, not just 6 or 4:
      1. Uses the site's own observed mean angle as the reference instead
         of a theoretical ideal polyhedral angle (not well-defined for an
         arbitrary/fractional CN).
      2. Includes ALL pairwise ligand-center-ligand angles, not just the
         12 "cis" (~90°) angles Robinson's octahedral formula uses --
         there's no general way to separate "cis" from "trans" pairs
         outside an assumed octahedral topology. For a real 6-coordinate
         site this pulls in the ~180° "trans" angles too, inflating the
         variance well above textbook octahedral BAV values.
    Not directly comparable to literature BAV numbers as a result; still
    a valid, self-consistent distortion metric for comparing sites/species
    within one report. Needs at least 2 angles (3 neighbors); returns None
    otherwise."""
    if len(angles) < 2:
        return None
    theta_avg = np.mean(angles)
    n = len(angles)
    return float(sum((theta - theta_avg) ** 2 for theta in angles) / (n - 1))

def _polyhedron_volume(neighbor_coords):
    """Volume (Å³) of the convex hull of a site's neighbor positions --
    needs at least 4 non-coplanar neighbors to bound a 3D volume; returns
    None otherwise (too few neighbors, or they're coplanar/collinear)."""
    if len(neighbor_coords) < 4:
        return None
    try:
        return float(ConvexHull(np.array(neighbor_coords)).volume)
    except QhullError:
        return None

def _classify_connectivity(shared_ligand_count):
    """How two same-species coordination polyhedra are linked, by how many
    ligand atoms they have in common: 0 = not connected (not reported),
    1 = corner-sharing, 2 = edge-sharing, 3+ = face-sharing."""
    if shared_ligand_count <= 0:
        return None
    if shared_ligand_count == 1:
        return "corner-sharing"
    if shared_ligand_count == 2:
        return "edge-sharing"
    return "face-sharing"

def compute_rdf(structure, r_max=10.0, dr=0.05):
    """Partial (per species pair, canonical order spA <= spB) and total
    radial distribution functions g(r) out to r_max, via a periodic-aware
    neighbor search (Structure.get_all_neighbors -- correctly includes
    contributions from periodic images beyond the first unit cell).

    Standard number-density normalization:
        g_AB(r) = <n_AB(r, r+dr)> / (4*pi*r^2*dr*rho_B)
    averaged over all A-type reference atoms; for a heteronuclear pair
    (A != B) only A-as-reference is accumulated (B-as-reference for the
    same physical pair is skipped) to avoid double counting the same
    proximity event with mismatched normalization -- by construction
    g_AB(r) computed from A's side and g_BA(r) computed from B's side
    describe the same physical distribution, so only one is needed. For a
    homonuclear pair (A == A), both (i,j) and (j,i) orderings legitimately
    contribute (each A atom's own neighbor count), which this naturally
    includes. Total g(r) uses every atom as both reference and target.

    Returns (r_centers, g_total, g_by_pair) where g_by_pair is
    {(spA, spB): array} for spA <= spB.
    """
    n_bins = max(1, int(np.ceil(r_max / dr)))
    bin_edges = np.linspace(0, n_bins * dr, n_bins + 1)
    r_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    volume = structure.lattice.volume

    species_of = [str(site.specie.symbol) for site in structure]
    species_list = sorted(set(species_of))
    species_count = {sp: species_of.count(sp) for sp in species_list}

    all_neighbors = structure.get_all_neighbors(r_max)

    total_hist = np.zeros(n_bins)
    directional_hist = {}

    for i in range(len(structure)):
        sp_i = species_of[i]
        for nbr in all_neighbors[i]:
            d = nbr.nn_distance
            if d <= 0 or d >= r_max:
                continue
            b = min(int(d / dr), n_bins - 1)
            total_hist[b] += 1
            sp_j = str(nbr.specie.symbol)
            ref_species = sp_i if sp_i <= sp_j else sp_j
            if sp_i != ref_species:
                continue
            key = (sp_i, sp_j) if sp_i <= sp_j else (sp_j, sp_i)
            directional_hist.setdefault(key, np.zeros(n_bins))
            directional_hist[key][b] += 1

    shell_volumes = 4.0 * np.pi * r_centers ** 2 * dr
    n_total = len(structure)
    rho_total = n_total / volume
    g_total = total_hist / (shell_volumes * rho_total * n_total)

    g_by_pair = {}
    for (sp_a, sp_b), hist in directional_hist.items():
        n_a = species_count[sp_a]
        rho_b = species_count[sp_b] / volume
        g_by_pair[(sp_a, sp_b)] = hist / (shell_volumes * rho_b * n_a)

    return r_centers, g_total, g_by_pair

def write_rdf_file(path, r_centers, g_total, g_by_pair):
    """Writes the full g(r) curve (total + every species pair) as columns,
    for external plotting -- the main report only quotes each curve's
    first peak."""
    pairs = sorted(g_by_pair)
    header_cols = ["r(Ang)", "g_total(r)"] + [f"g_{a}-{b}(r)" for a, b in pairs]
    header = "  ".join(f"{c:>14}" for c in header_cols)
    with open(path, "w") as f:
        f.write("# " + header + "\n")
        for k, r in enumerate(r_centers):
            row = [r, g_total[k]] + [g_by_pair[pair][k] for pair in pairs]
            f.write("  ".join(f"{v:14.6f}" for v in row) + "\n")

def _first_peak(r_centers, g):
    """(r, g(r)) at the highest point of the curve, or (None, None) if it's
    all zero (no pair found within r_max)."""
    if len(g) == 0 or not np.any(g):
        return None, None
    idx = int(np.argmax(g))
    return float(r_centers[idx]), float(g[idx])

def compute_ecn(structure, mode, atoms_position=None):
    """Pure computation -- no printing, no file I/O. Returns a results dict
    that format_report() turns into the console/file report."""
    with warnings.catch_warnings():
        # CrystalNN (used throughout this function) unconditionally warns
        # twice per neighbor pair that "no oxidation states specified",
        # then again that it's falling back to covalent/atomic radii --
        # pure noise here: SIESTA structures never carry oxidation states
        # (guessing them reliably isn't possible for the covalent/metallic/
        # mixed systems this suite supports -- the same reasoning behind
        # not implementing bond-valence-sum-based metrics), and the
        # covalent/atomic-radius fallback CrystalNN uses instead is exactly
        # the radius data every other coordination-number method here
        # already uses. So this is an expected, permanent, non-actionable
        # warning, not a defect -- filtered by exact message so it's scoped
        # to only these two, and restored on return via catch_warnings().
        # CrystalNN's separate distance_cutoffs refinement (which also
        # consults this same radius lookup) is deliberately left fully
        # enabled: disabling it outright to dodge the warning changes
        # computed coordination numbers by up to ~0.9, a real accuracy
        # loss, not a cosmetic one -- verified by hand before choosing this
        # approach over passing distance_cutoffs=None.
        warnings.filterwarnings("ignore", message="No oxidation states specified.*")
        warnings.filterwarnings("ignore", message="CrystalNN: cannot locate.*")
        return _compute_ecn_impl(structure, mode, atoms_position)

def _compute_ecn_impl(structure, mode, atoms_position=None):
    lattice = structure.lattice
    results = {
        "mode": mode,
        "lattice": {
            "a": lattice.a, "b": lattice.b, "c": lattice.c,
            "alpha": lattice.alpha, "beta": lattice.beta, "gamma": lattice.gamma,
            "volume": lattice.volume,
            "density": float(structure.density),
            "vectors": lattice.matrix,
        },
    }

    # Coordination-number methods. use_weights=True makes every one of
    # them -- EconNN included, its get_cn() respects the flag exactly
    # like the rest even though its formula (Hoppe's ECoN) is defined
    # as a continuous quantity: verified this atom-by-atom, e.g.
    # use_weights=False gives a plain integer count (6) where
    # use_weights=True gives the real ECoN value (5.98) -- return a
    # genuinely "effective" (continuous, neighbor-weighted) coordination
    # number. pymatgen's NearNeighbors.get_cn() defaults to
    # use_weights=False (a plain integer neighbor count) for all four.
    #
    # JmolNN was dropped: its reference distance is a fixed Jmol
    # bonding-radius lookup table, not this atom's own closest-neighbor
    # distance, so its weight isn't on the same scale as the other four
    # and it was consistently the outlier of the five -- worst agreement
    # with the actual local geometry, not just "different."
    methods = {
        "MinDistNN": MinimumDistanceNN(),
        # weighted_cn=True must match the use_weights=True passed to
        # get_cn() below -- CrystalNN raises ValueError otherwise,
        # unlike the other methods (which don't require a matching
        # constructor flag).
        "CrystalNN": CrystalNN(weighted_cn=True),
        "BrunnerNN": BrunnerNNRelative(),
        "EconNN": EconNN()
    }
    ecn_results = {method: [] for method in methods}

    pos_atomics = [(i+1, str(site.specie.symbol), site.coords) for i, site in enumerate(structure)]
    results["atomic_positions"] = pos_atomics
    species_of_all = [p[1] for p in pos_atomics]
    all_species_list = sorted(set(species_of_all))

    # Whole-structure, --mode-independent analysis: minimum same-species
    # distance and coordination-polyhedron connectivity both describe the
    # extended network, not just whichever atoms --list happens to name,
    # so they always use every atom regardless of --mode/--list. This
    # means a second, separate CrystalNN pass over all atoms even in
    # "list" mode's small subset -- deliberately not fused with the
    # --mode-scoped neighbor loop further down, to avoid changing what
    # that loop's pooled distance/angle statistics mean per mode.
    same_species_min_distance = {}
    for sp in all_species_list:
        idx = [i for i, s in enumerate(species_of_all) if s == sp]
        if len(idx) < 2:
            continue
        min_d = min(structure.get_distance(idx[a], idx[b])
                     for a in range(len(idx)) for b in range(a + 1, len(idx)))
        same_species_min_distance[sp] = float(min_d)
    results["same_species_min_distance"] = same_species_min_distance

    connectivity_cnn = CrystalNN()
    connectivity_by_species = {}
    try:
        full_neighbor_ids = {}
        for i in range(len(structure)):
            nbrs = connectivity_cnn.get_nn_info(structure, i)
            full_neighbor_ids[i] = {
                (n["site_index"], tuple(round(x) for x in n.get("image", (0, 0, 0))))
                for n in nbrs
            }
        for sp in all_species_list:
            idx = [i for i, s in enumerate(species_of_all) if s == sp]
            counts = {"corner-sharing": 0, "edge-sharing": 0, "face-sharing": 0}
            for a in range(len(idx)):
                for b in range(a + 1, len(idx)):
                    shared = len(full_neighbor_ids[idx[a]] & full_neighbor_ids[idx[b]])
                    label = _classify_connectivity(shared)
                    if label is not None:
                        counts[label] += 1
            connectivity_by_species[sp] = counts
    except Exception as e:
        print(f"[WARNING] Failed to compute polyhedron connectivity: {e}")
    results["connectivity_by_species"] = connectivity_by_species

    if mode == "mean":
        for i in range(len(structure)):
            for method_name, method in methods.items():
                try:
                    ecn_results[method_name].append(method.get_cn(structure, i, use_weights=True))
                except Exception as e:
                    ecn_results[method_name].append(None)
                    print(f"[WARNING] {method_name} failed for atom {i+1}: {e}")

        # Per-species averages -- a single structure-wide average mixes
        # chemically distinct sites (e.g. a cation's coordination with
        # an anion's) into one physically meaningless number. Reported
        # per species, plus an overall figure for reference.
        species_by_index = [pos_atomics[i][1] for i in range(len(structure))]
        species_list = sorted(set(species_by_index))

        cn_per_species = {}
        for sp in species_list:
            sp_idx = [i for i, s in enumerate(species_by_index) if s == sp]
            cn_per_species[sp] = {
                "n_atoms": len(sp_idx),
                "values": {
                    method_name: _safe_mean([values[i] for i in sp_idx])
                    for method_name, values in ecn_results.items()
                },
            }
        results["cn_per_species"] = cn_per_species
        results["cn_overall"] = {
            method_name: _safe_mean(values) for method_name, values in ecn_results.items()
        }

    elif mode == "list" and atoms_position:
        for i in atoms_position:
            for method_name, method in methods.items():
                try:
                    ecn_results[method_name].append(method.get_cn(structure, i-1, use_weights=True))
                except Exception as e:
                    ecn_results[method_name].append(None)
                    print(f"[WARNING] {method_name} failed for atom {i}: {e}")

        cn_per_atom = []
        for k, atom_index in enumerate(atoms_position):
            atom_id, species, position = pos_atomics[atom_index - 1]
            cn_per_atom.append({
                "atom_id": atom_id,
                "species": species,
                "position": position,
                "values": {method: values[k] for method, values in ecn_results.items()},
            })
        results["cn_per_atom"] = cn_per_atom

    # Bond distances/angles/polyhedron shape, all from one CrystalNN
    # neighbor list per atom -- reused for everything below instead of
    # adding yet another per-method comparison (the 5-method table above
    # is specifically about comparing coordination-number *algorithms*;
    # this geometric analysis just needs one consistent neighbor set, the
    # common choice in the literature).
    #
    # Distances/angles are pooled by species pair/triplet across the
    # structure (a single average lumping e.g. Sn-O with any Sn-Sn/O-O
    # found isn't meaningful on its own); a bond or angle can be counted
    # from more than one of its atoms' neighbor lists in "mean" mode (n=
    # makes that visible) -- this only skews the pooled average if the
    # neighbor relationship isn't symmetric between the atoms involved.
    #
    # Bond-length distortion (BLD), bond-angle variance (BAV), and
    # polyhedron volume are inherently per-atom (they characterize one
    # coordination polyhedron), then averaged per species like ECN.
    cnn = CrystalNN()
    distances_by_pair = {}
    all_distances = []
    angles_by_triplet = {}
    all_angles = []
    bld_by_atom = {}
    bav_by_atom = {}
    volume_by_atom = {}

    indices = range(len(structure)) if mode == "mean" else [i - 1 for i in atoms_position]

    for i in indices:
        try:
            neighbors = cnn.get_nn_info(structure, i)
            center_coords = structure[i].coords

            site_distances = []
            for neighbor in neighbors:
                dist = neighbor["site"].distance(structure[i])
                site_distances.append(dist)
                all_distances.append(dist)
                pair = tuple(sorted((pos_atomics[i][1], str(neighbor["site"].specie.symbol))))
                distances_by_pair.setdefault(pair, []).append(dist)
            bld = _bond_length_distortion(site_distances)
            if bld is not None:
                bld_by_atom[i] = bld
            volume = _polyhedron_volume([n["site"].coords for n in neighbors])
            if volume is not None:
                volume_by_atom[i] = volume

            site_angles = []
            for a in range(len(neighbors)):
                for b in range(a + 1, len(neighbors)):
                    v1 = neighbors[a]["site"].coords - center_coords
                    v2 = neighbors[b]["site"].coords - center_coords
                    cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
                    angle_deg = float(np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0))))
                    site_angles.append(angle_deg)
                    all_angles.append(angle_deg)
                    ligand_pair = tuple(sorted((str(neighbors[a]["site"].specie.symbol),
                                                 str(neighbors[b]["site"].specie.symbol))))
                    triplet = (ligand_pair[0], pos_atomics[i][1], ligand_pair[1])
                    angles_by_triplet.setdefault(triplet, []).append(angle_deg)
            bav = _bond_angle_variance(site_angles)
            if bav is not None:
                bav_by_atom[i] = bav
        except Exception as e:
            print(f"[WARNING] Failed to compute distances/angles for atom {i+1}: {e}")

    results["bond_distances_by_pair"] = {
        pair: (float(np.mean(d)), len(d)) for pair, d in distances_by_pair.items()
    }
    results["bond_distances_overall"] = (
        (float(np.mean(all_distances)), len(all_distances)) if all_distances else None
    )
    results["bond_angles_by_triplet"] = {
        triplet: (float(np.mean(a)), len(a)) for triplet, a in angles_by_triplet.items()
    }
    results["bond_angles_overall"] = (
        (float(np.mean(all_angles)), len(all_angles)) if all_angles else None
    )

    if mode == "mean":
        distortion_per_species = {}
        for sp in species_list:
            sp_idx = [i for i, s in enumerate(species_by_index) if s == sp]
            distortion_per_species[sp] = {
                "n_atoms": len(sp_idx),
                "bld": _safe_mean([bld_by_atom.get(i) for i in sp_idx]),
                "bav": _safe_mean([bav_by_atom.get(i) for i in sp_idx]),
                "volume": _safe_mean([volume_by_atom.get(i) for i in sp_idx]),
            }
        results["distortion_per_species"] = distortion_per_species
    else:
        distortion_per_atom = []
        for atom_index in atoms_position:
            i = atom_index - 1
            distortion_per_atom.append({
                "atom_id": pos_atomics[i][0],
                "species": pos_atomics[i][1],
                "bld": bld_by_atom.get(i),
                "bav": bav_by_atom.get(i),
                "volume": volume_by_atom.get(i),
            })
        results["distortion_per_atom"] = distortion_per_atom

    return results

# --- Report formatting --------------------------------------------------
_WIDTH = 74

def _rule(char="-"):
    return char * _WIDTH

def _fmt(value, width=7, prec=3):
    return f"{value:{width}.{prec}f}" if value is not None else "N/A".rjust(width)

def format_report(results, source_file, fmt, rdf_summary=None):
    lat = results["lattice"]
    lines = []
    lines.append(_rule("="))
    lines.append("STRUCTURAL PROPERTIES REPORT - STB Suite".center(_WIDTH))
    lines.append(_rule("="))
    lines.append("")
    lines.append(f"Generated        : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"Source file      : {source_file}  (format: {fmt})")
    lines.append(f"Analysis mode    : {results['mode']}")

    lines.append("")
    lines.append(_rule())
    lines.append("LATTICE")
    lines.append(_rule())
    lines.append(f"a = {lat['a']:.3f} Å   b = {lat['b']:.3f} Å   c = {lat['c']:.3f} Å")
    lines.append(f"alpha = {lat['alpha']:.2f}°   beta = {lat['beta']:.2f}°   gamma = {lat['gamma']:.2f}°")
    lines.append(f"Volume = {lat['volume']:.3f} Å³   Density = {lat['density']:.3f} g/cm³")
    lines.append("")
    lines.append("Lattice vectors (Å):")
    for i, vec in enumerate(lat["vectors"]):
        lines.append(f"  a_{i+1}: {vec[0]:12.6f}  {vec[1]:12.6f}  {vec[2]:12.6f}")

    if results["mode"] == "mean":
        lines.append("")
        lines.append(_rule())
        lines.append("EFFECTIVE COORDINATION NUMBER (weighted), PER SPECIES")
        lines.append(_rule())
        for sp, data in results["cn_per_species"].items():
            lines.append(f"{sp} ({data['n_atoms']} atoms):")
            for method, avg in data["values"].items():
                lines.append(f"  {method:<15}: {_fmt(avg)}")
            lines.append("")
        lines.append("Overall average:")
        for method, avg in results["cn_overall"].items():
            lines.append(f"  {method:<15}: {_fmt(avg)}")

    elif results["mode"] == "list":
        lines.append("")
        lines.append(_rule())
        lines.append("EFFECTIVE COORDINATION NUMBER (weighted), SPECIFIED ATOMS")
        lines.append(_rule())
        for atom in results["cn_per_atom"]:
            pos = atom["position"]
            lines.append(f"Atom {atom['atom_id']} ({atom['species']}), position: "
                         f"{pos[0]:.6f}  {pos[1]:.6f}  {pos[2]:.6f}")
            for method, val in atom["values"].items():
                lines.append(f"  {method:<15}: {_fmt(val)}")
            lines.append("")

    lines.append(_rule())
    lines.append("AVERAGE BOND DISTANCE, PER SPECIES PAIR")
    lines.append(_rule())
    if results["bond_distances_overall"] is not None:
        for pair in sorted(results["bond_distances_by_pair"]):
            avg, n = results["bond_distances_by_pair"][pair]
            label = f"{pair[0]}-{pair[1]}"
            lines.append(f"  {label:<10}: {avg:.4f} Å  (n={n})")
        overall_avg, overall_n = results["bond_distances_overall"]
        lines.append(f"  {'Overall':<10}: {overall_avg:.4f} Å  (n={overall_n})")
    else:
        lines.append("  No distances could be computed.")

    lines.append("")
    lines.append(_rule())
    lines.append("AVERAGE BOND ANGLE, PER LIGAND-CENTER-LIGAND TRIPLET")
    lines.append(_rule())
    if results["bond_angles_overall"] is not None:
        for triplet in sorted(results["bond_angles_by_triplet"]):
            avg, n = results["bond_angles_by_triplet"][triplet]
            label = f"{triplet[0]}-{triplet[1]}-{triplet[2]}"
            lines.append(f"  {label:<12}: {avg:7.3f}°  (n={n})")
        overall_avg, overall_n = results["bond_angles_overall"]
        lines.append(f"  {'Overall':<12}: {overall_avg:7.3f}°  (n={overall_n})")
    else:
        lines.append("  No angles could be computed (fewer than 2 neighbors per site).")

    lines.append("")
    lines.append(_rule())
    lines.append("COORDINATION POLYHEDRON DISTORTION")
    lines.append(_rule())
    lines.append("BLD: bond-length distortion (%, Baur). BAV: bond-angle variance (deg²),")
    lines.append("from ALL ligand-center-ligand angles (not just Robinson's 12 'cis' angles")
    lines.append("for an ideal octahedron -- there's no general cis/trans split for an")
    lines.append("arbitrary/fractional CN), around the site's own mean angle (not a")
    lines.append("theoretical ideal) -- not directly comparable to textbook BAV values, but")
    lines.append("self-consistent for comparing sites/species within this report. Volume:")
    lines.append("convex hull of the neighbor positions (Å³, needs >= 4 non-coplanar neighbors).")
    lines.append("")
    if results["mode"] == "mean":
        for sp, data in results["distortion_per_species"].items():
            lines.append(f"{sp} ({data['n_atoms']} atoms):")
            lines.append(f"  BLD   : {_fmt(data['bld'])} %")
            lines.append(f"  BAV   : {_fmt(data['bav'])} deg²")
            lines.append(f"  Volume: {_fmt(data['volume'])} Å³")
            lines.append("")
    else:
        for atom in results["distortion_per_atom"]:
            lines.append(f"Atom {atom['atom_id']} ({atom['species']}):")
            lines.append(f"  BLD   : {_fmt(atom['bld'])} %")
            lines.append(f"  BAV   : {_fmt(atom['bav'])} deg²")
            lines.append(f"  Volume: {_fmt(atom['volume'])} Å³")
            lines.append("")

    lines.append(_rule())
    lines.append("SAME-SPECIES MINIMUM DISTANCE (whole structure, independent of --mode)")
    lines.append(_rule())
    if results["same_species_min_distance"]:
        for sp in sorted(results["same_species_min_distance"]):
            label = f"{sp}-{sp}"
            lines.append(f"  {label:<10}: {results['same_species_min_distance'][sp]:.4f} Å")
    else:
        lines.append("  Fewer than 2 atoms of any single species -- nothing to compare.")

    lines.append("")
    lines.append(_rule())
    lines.append("COORDINATION POLYHEDRON CONNECTIVITY (shared ligands between same-")
    lines.append("species centers; whole structure, independent of --mode)")
    lines.append(_rule())
    if results["connectivity_by_species"]:
        for sp, counts in results["connectivity_by_species"].items():
            lines.append(f"{sp}-{sp}:")
            lines.append(f"  corner-sharing: {counts['corner-sharing']} pair(s)")
            lines.append(f"  edge-sharing  : {counts['edge-sharing']} pair(s)")
            lines.append(f"  face-sharing  : {counts['face-sharing']} pair(s)")
            lines.append("")
    else:
        lines.append("  Not available.")
        lines.append("")

    if rdf_summary is not None:
        lines.append(_rule())
        lines.append("RADIAL DISTRIBUTION FUNCTION g(r)")
        lines.append(_rule())
        lines.append(f"Full curve written to {rdf_summary['rdf_file']} "
                     f"(r_max = {rdf_summary['r_max']:.1f} Å, bin width = {rdf_summary['dr']:.2f} Å).")
        lines.append("First peak (highest g(r) in this range):")
        for label, (r, g) in rdf_summary["peaks"].items():
            if r is None:
                lines.append(f"  {label:<10}: no pair found within r_max")
            else:
                lines.append(f"  {label:<10}: r = {r:6.3f} Å   g(r) = {g:8.3f}")
        lines.append("")

    lines.append(_rule())
    lines.append("ATOMIC POSITIONS (Cartesian, Å)")
    lines.append(_rule())
    for atom_id, species, coords in results["atomic_positions"]:
        lines.append(f"{atom_id:>4}  {species:<3}  {coords[0]:12.6f}  {coords[1]:12.6f}  {coords[2]:12.6f}")

    lines.append(_rule("="))
    return "\n".join(lines) + "\n"

def main():
    parser = argparse.ArgumentParser(
        description="Compute ECN and structural properties from a SIESTA structure file.",
        epilog="Example usage:\n"
               "  stb-structural --file structure.fdf --format fdf --mode mean\n"
               "  stb-structural --file siesta.STRUCT_OUT --format struct_out --mode list --list 1,4,5",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--file", required=True, help="Path to structure file.")
    parser.add_argument("--format", required=True, choices=["fdf", "struct_out"],
                        help="Input file format:\n"
                             "  fdf:        SIESTA structure input (%%block LatticeVectors etc.)\n"
                             "  struct_out: SIESTA post-relaxation output (.STRUCT_OUT)")
    parser.add_argument("--mode", choices=["list", "mean"], required=True, help="Calculation mode: list or mean")
    parser.add_argument("--list", type=str, help="List of atom indices (comma-separated, 1-based). Example: 1,4,5,7 - Required for 'list' mode")
    parser.add_argument("-o", "--output-dir", type=str, default=".",
                        help="Directory to write structural_information.dat, warnings.log, and "
                             "rdf.dat into (default: current directory). Created if it doesn't exist.")
    parser.add_argument("--no-rdf", dest="rdf", action="store_false",
                        help="Skip the radial distribution function g(r) (rdf.dat is not written).")
    parser.add_argument("--rdf-rmax", type=float, default=10.0,
                        help="Cutoff radius for g(r), in Å (default: 10.0). Larger values capture "
                             "more coordination shells but take longer for large structures.")
    parser.add_argument("-v", "--version", action="version",
                        version=f"stb-structural {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    args = parser.parse_args()

    if args.mode == "list" and not args.list:
        parser.error("--list is required when --mode is 'list'")
    if args.rdf_rmax <= 0:
        parser.error("--rdf-rmax must be positive.")

    os.makedirs(args.output_dir, exist_ok=True)

    # Configure logger for warnings (done here, not at module level, so importing
    # stb.structural has no side effect of creating warnings.log on disk).
    # filemode='w' so stale warnings from a previous run in the same directory
    # don't linger forever (logging.basicConfig defaults to append mode).
    logging.basicConfig(filename=os.path.join(args.output_dir, "warnings.log"),
                         level=logging.WARNING, format="%(message)s", filemode='w')
    warnings.showwarning = warn_handler

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("STRUCTURAL PROPERTIES:", 'bold'))
    print("-"*60)

    print("\n[INFO] Reading structure file...")
    atoms_position = list(map(int, args.list.strip('[]').split(','))) if args.list else None

    try:
        structure = read_structure(args.file, args.format)
    except FileNotFoundError:
        print(color_text(f"[ERROR] Structure file '{args.file}' not found.", 'red'))
        sys.exit(1)
    except ValueError as e:
        print(color_text(f"[ERROR] {e}", 'red'))
        sys.exit(1)

    if atoms_position:
        invalid = [i for i in atoms_position if i < 1 or i > len(structure)]
        if invalid:
            print(color_text(f"[ERROR] Atom index/indices {invalid} out of range "
                              f"(structure has {len(structure)} atoms, 1-based).", 'red'))
            sys.exit(1)

    print("[INFO] Computing coordination numbers and bond distances...")
    results = compute_ecn(structure, args.mode, atoms_position)

    rdf_summary = None
    if args.rdf:
        print(f"[INFO] Computing radial distribution function g(r) (r_max = {args.rdf_rmax:.1f} Å)...")
        r_centers, g_total, g_by_pair = compute_rdf(structure, r_max=args.rdf_rmax)
        rdf_path = os.path.join(args.output_dir, "rdf.dat")
        write_rdf_file(rdf_path, r_centers, g_total, g_by_pair)
        peaks = {"Total": _first_peak(r_centers, g_total)}
        for pair in sorted(g_by_pair):
            peaks[f"{pair[0]}-{pair[1]}"] = _first_peak(r_centers, g_by_pair[pair])
        rdf_summary = {"rdf_file": rdf_path, "r_max": args.rdf_rmax, "dr": 0.05, "peaks": peaks}
    else:
        # Remove a stale rdf.dat from a previous (non-"--no-rdf") run in
        # this same directory, so it can't be mistaken for output of this run.
        stale_rdf = os.path.join(args.output_dir, "rdf.dat")
        if os.path.exists(stale_rdf):
            os.remove(stale_rdf)

    report = format_report(results, args.file, args.format, rdf_summary)
    print("\n" + report)

    out_path = os.path.join(args.output_dir, "structural_information.dat")
    with open(out_path, "w") as f:
        f.write(report)

    print(f"[INFO] Job complete! Results saved to {out_path}")
    print("-"*60)

if __name__ == "__main__":
    main()
