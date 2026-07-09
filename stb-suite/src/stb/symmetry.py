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
from datetime import datetime
import numpy as np
import spglib
from pymatgen.core.operations import SymmOp
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from stb.core import kspace, structure_io
from stb.core.cli import color_text, show_intro

# Same default as stb-kgrid/stb-mlrelax/etc. (core/kspace.py's other
# callers): minimum empty span along an axis, wrapped periodically, to
# treat it as vacuum-padded rather than genuinely periodic.
VACUUM_GAP_ANG = 10.0

def compute_layer_group(structure, aperiodic_dir, symprec=1e-3):
    """Layer-group analysis for a genuinely 2D-periodic structure (periodic
    along the 2 axes other than `aperiodic_dir`, vacuum-padded along it) --
    the physically correct symmetry classification for a slab, unlike the
    ordinary 3D space group SPACE GROUP reports (which treats the vacuum as
    just an unusually tall periodic cell and can report a group that
    depends on the arbitrary vacuum thickness, not the real 2D symmetry).

    Uses spglib.get_layergroup() directly (needs spglib >= 2.1.0, its
    Python interface for the 80 layer groups; >= 2.2.0 for a Hall-number
    bugfix -- see pyproject.toml's spglib pin) since pymatgen's
    SpacegroupAnalyzer doesn't wrap this yet. Returns None if spglib
    can't determine a layer group (e.g. not genuinely 2D-periodic within
    symprec) or if the installed spglib predates this function.
    """
    if not hasattr(spglib, "get_layergroup"):
        return None

    cell = (structure.lattice.matrix, structure.frac_coords,
            [site.specie.Z for site in structure])
    ds = spglib.get_layergroup(cell, aperiodic_dir=aperiodic_dir, symprec=symprec)
    if ds is None:
        return None

    # ds.equivalent_atoms[i] is the representative atom index of atom i's
    # orbit (0-based, same convention as compute_symmetry's 3D path, but
    # spglib returns it flat here rather than pre-grouped into lists like
    # pymatgen's ss.equivalent_indices).
    orbit_members = {}
    for i, rep in enumerate(ds.equivalent_atoms):
        orbit_members.setdefault(rep, []).append(i)

    sites = []
    for i, site in enumerate(structure):
        rep = ds.equivalent_atoms[i]
        wyckoff = f"{len(orbit_members[rep])}{ds.wyckoffs[i]}"
        sites.append({
            "atom_id": i + 1,
            "species": str(site.specie.symbol),
            "wyckoff": wyckoff,
            "site_symmetry": ds.site_symmetry_symbols[i],
            "orbit": rep + 1,
        })

    orbits = [
        {"wyckoff": f"{len(indices)}{ds.wyckoffs[rep]}",
         "species": str(structure[rep].specie.symbol), "n_atoms": len(indices),
         "example_atom_id": indices[0] + 1, "site_symmetry": ds.site_symmetry_symbols[rep]}
        for rep, indices in orbit_members.items()
    ]

    # Translation components come back as noisy near-integers/near-zeros
    # (e.g. 9.99999980e-01) rather than exact 0/1 -- round before %1 so
    # as_xyz_str() doesn't render spurious "+1"/"+0" terms.
    op_objects = [
        SymmOp.from_rotation_and_translation(rot, np.round(trans, 6) % 1.0)
        for rot, trans in zip(ds.rotations, ds.translations)
    ]

    distortion = _wyckoff_distortion(structure, op_objects)
    for site in sites:
        site["distortion"] = distortion[site["atom_id"] - 1]

    return {
        "symbol": ds.international,
        "number": ds.number,
        "hall_symbol": ds.hall,
        "hall_number": ds.hall_number,
        "point_group": ds.pointgroup,
        "aperiodic_axis": "abc"[aperiodic_dir],
        "n_distinct_sites": len(orbit_members),
        "sites": sites,
        "orbits": orbits,
        "symmetry_operations": [op.as_xyz_str() for op in op_objects],
        # For --write-refined: the standardized (layer-group-consistent)
        # cell spglib derived internally -- same fields/convention as the
        # ordinary 3D SpglibDataset's std_lattice/std_positions/std_types,
        # but standardized against the layer group instead of the (for a
        # vacuum-padded structure, not physically meaningful) 3D space
        # group core/symmetry.py::reduce_to_unitcell() would otherwise use.
        "std_lattice": ds.std_lattice,
        "std_positions": ds.std_positions,
        "std_types": ds.std_types,
    }

def compute_symmetry(structure, symprec=1e-3, angle_tolerance=5.0):
    """Pure computation -- no printing, no file I/O. Returns a results dict
    that format_report() turns into the console/file report.

    Uses pymatgen's SpacegroupAnalyzer (wrapping spglib) -- the same
    approach already shared by stb-unitcell/stb-fetch via
    core/symmetry.py::reduce_to_unitcell(), instead of calling spglib
    directly with a hand-rolled crystal-system classifier: SpacegroupAnalyzer
    provides the crystal system (and Pearson symbol, lattice type, per-orbit
    Wyckoff multiplicity) natively, so there's nothing left to duplicate.
    """
    sga = SpacegroupAnalyzer(structure, symprec=symprec, angle_tolerance=angle_tolerance)
    ss = sga.get_symmetrized_structure()
    # site_symmetry_symbols is the per-atom point-group symmetry AT that
    # site (e.g. "m-3m" for a rock-salt ion sitting on a high-symmetry
    # position) -- same original atom order as ss.wyckoff_letters, verified
    # against ds.wyckoffs (which matches ss.wyckoff_letters exactly).
    dataset = sga.get_symmetry_dataset()
    site_symmetry = dataset.site_symmetry_symbols

    # ss.wyckoff_symbols is per ORBIT (e.g. "4a" for a 4-atom equivalence
    # class), not per atom -- ss.equivalent_indices maps each orbit to its
    # member atom indices (0-based, original atom order preserved). Rebuild
    # a per-atom "multiplicity+letter" label by looking up which orbit each
    # atom belongs to.
    orbit_of_atom = {}
    for orbit_idx, indices in enumerate(ss.equivalent_indices):
        for i in indices:
            orbit_of_atom[i] = orbit_idx

    sites = []
    for i, site in enumerate(structure):
        orbit_idx = orbit_of_atom[i]
        wyckoff = f"{len(ss.equivalent_indices[orbit_idx])}{ss.wyckoff_letters[i]}"
        sites.append({
            "atom_id": i + 1,
            "species": str(site.specie.symbol),
            "wyckoff": wyckoff,
            "site_symmetry": site_symmetry[i],
            "orbit": orbit_idx + 1,
            "frac_coords": site.frac_coords,
        })

    lattice = structure.lattice
    ops = sga.get_symmetry_operations()
    distortion = _wyckoff_distortion(structure, ops)
    for site in sites:
        site["distortion"] = distortion[site["atom_id"] - 1]

    reduced_formula, z = structure.composition.get_reduced_formula_and_factor()
    vacuum_axes = kspace.detect_vacuum_axes(structure.frac_coords, lattice.matrix, VACUUM_GAP_ANG)

    # Layer-group detection only makes sense for a genuinely 2D-periodic
    # structure -- exactly one aperiodic (vacuum) axis. Two or more vacuum
    # axes (a wire or an isolated molecule) has no layer-group analogue;
    # zero means it's not vacuum-padded at all.
    layer_group = None
    if sum(vacuum_axes) == 1:
        layer_group = compute_layer_group(structure, vacuum_axes.index(True), symprec)

    return {
        "vacuum_axes": vacuum_axes,
        "layer_group": layer_group,
        "space_group_symbol": sga.get_space_group_symbol(),
        "space_group_number": sga.get_space_group_number(),
        "hall_symbol": sga.get_hall(),
        "hall_number": dataset.hall_number,
        "point_group": sga.get_point_group_symbol(),
        "crystal_system": sga.get_crystal_system(),
        "lattice_type": sga.get_lattice_type(),
        "pearson_symbol": sga.get_pearson_symbol(),
        "setting_choice": dataset.choice,
        "symprec": symprec,
        "angle_tolerance": angle_tolerance,
        "lattice": {
            "a": lattice.a, "b": lattice.b, "c": lattice.c,
            "alpha": lattice.alpha, "beta": lattice.beta, "gamma": lattice.gamma,
            "volume": lattice.volume,
            "vectors": lattice.matrix,
            "reduced_formula": reduced_formula,
            "z": z,
        },
        "sites": sites,
        "n_distinct_sites": len(ss.equivalent_indices),
        "orbits": [
            {"wyckoff": ss.wyckoff_symbols[k], "species": str(structure[indices[0]].specie.symbol),
             "n_atoms": len(indices), "example_atom_id": indices[0] + 1,
             "site_symmetry": site_symmetry[indices[0]]}
            for k, indices in enumerate(ss.equivalent_indices)
        ],
        "symmetry_operations": [op.as_xyz_str() for op in ops],
    }

def _wyckoff_distortion(structure, ops):
    """Per-atom distortion (Å) from exact symmetry: applies every DETECTED
    symmetry operation to each atom's position and measures the distance to
    the nearest same-species atom that operation should map it onto exactly
    if the structure had perfect symmetry; an atom's distortion is the
    largest such distance over all non-identity operations.

    This is bounded by construction to roughly <= the symprec that produced
    `ops` (that's the definition of symprec: the tolerance within which
    spglib accepted this operation set as valid), so it gives a finer,
    per-atom breakdown *within* that tolerance budget rather than just the
    pass/fail group determination.

    (An earlier, simpler attempt used SpacegroupAnalyzer.get_refined_structure()
    -- its docstring claims "atoms moved to the expected symmetry positions"
    -- but empirically it returns the input positions essentially unchanged
    (verified: every atom's distance to its refined counterpart came out
    exactly 0.0 even for a deliberately noisy test structure), so it
    couldn't be used for this. Recomputing the residuals directly from the
    operations themselves, as done here, does produce the expected
    non-zero, noise-correlated values -- verified against a structure with
    known injected per-atom noise.)
    """
    lattice = structure.lattice
    species_of = [str(s.specie.symbol) for s in structure]
    indices_by_species = {}
    for i, sp in enumerate(species_of):
        indices_by_species.setdefault(sp, []).append(i)

    distortion = []
    for i, site in enumerate(structure):
        candidates = indices_by_species[species_of[i]]
        max_dist = 0.0
        for op in ops:
            p2 = op.operate(site.frac_coords) % 1.0
            best = min(
                lattice.get_distance_and_image(p2, structure[j].frac_coords)[0]
                for j in candidates
            )
            max_dist = max(max_dist, best)
        distortion.append(max_dist)
    return distortion

# Default tolerance sweep for --scan-symprec: log-spaced from tight (typical
# DFT-relaxation numerical noise) to loose enough to reveal near-symmetry a
# structure might have lost to that same noise.
DEFAULT_SYMPREC_SCAN = (1e-5, 1e-4, 1e-3, 1e-2, 5e-2, 1e-1)

def scan_symprec(structure, symprec_values=DEFAULT_SYMPREC_SCAN, angle_tolerance=5.0):
    """Runs symmetry detection at each symprec in symprec_values (ascending)
    and returns a list of (symprec, space_group_symbol, space_group_number)
    -- lets a caller see whether a structure's *apparent* symmetry is being
    hidden by small numerical noise (e.g. from a DFT relaxation) that a
    tighter tolerance is picking up as a genuine distortion."""
    results = []
    for sp in symprec_values:
        sga = SpacegroupAnalyzer(structure, symprec=sp, angle_tolerance=angle_tolerance)
        results.append((sp, sga.get_space_group_symbol(), sga.get_space_group_number()))
    return results

def scan_symprec_layer(structure, aperiodic_dir, symprec_values=DEFAULT_SYMPREC_SCAN):
    """Layer-group analogue of scan_symprec(): runs spglib.get_layergroup()
    at each symprec, for a vacuum-padded structure where the layer group
    (not the 3D space group) is the physically meaningful one to watch for
    tolerance sensitivity. An entry's symbol/number is None if
    get_layergroup() couldn't determine one at that tolerance."""
    results = []
    for sp in symprec_values:
        lg = compute_layer_group(structure, aperiodic_dir, sp)
        results.append((sp, lg["symbol"] if lg else None, lg["number"] if lg else None))
    return results

def _describe_symprec_scan(scan_results, label="space group"):
    """One-line summary of the scan: the first tolerance (if any) at which
    the group changes from the tightest-tolerance result. Only reports the
    first transition, not every subsequent one, if it changes more than
    once across the scanned range. A None symbol/number (detection failed
    at that tolerance, only possible for the layer-group scan) counts as
    its own distinct state, same as any other change."""
    first_symbol, first_number = scan_results[0][1], scan_results[0][2]
    for symprec, symbol, number in scan_results[1:]:
        if number != first_number:
            first_str = f"{first_symbol} (No. {first_number})" if first_number is not None else "not determined"
            cur_str = f"{symbol} (No. {number})" if number is not None else "not determined"
            return f"Symmetry changes at symprec >= {symprec:g}: {first_str} -> {cur_str}."
    return f"No change in detected {label} across the scanned tolerance range."

def compare_symmetry(results_a, results_b, label_a, label_b):
    """Compares two compute_symmetry() results (e.g. a structure before and
    after a SIESTA relaxation) -- same space group or not, and how their
    symmetry-operation counts differ. Both must have been computed with the
    same symprec/angle_tolerance for the comparison to mean anything.

    If both structures have a detected layer group (e.g. comparing a 2D
    slab before/after relaxation), that's compared too -- the 3D space
    group comparison alone would only describe the padded supercells,
    the same limitation SPACE GROUP has for a single structure.
    """
    layer_comparison = None
    lg_a, lg_b = results_a["layer_group"], results_b["layer_group"]
    if lg_a is not None and lg_b is not None:
        layer_comparison = {
            "symbol_a": lg_a["symbol"], "number_a": lg_a["number"],
            "symbol_b": lg_b["symbol"], "number_b": lg_b["number"],
            "n_ops_a": len(lg_a["symmetry_operations"]), "n_ops_b": len(lg_b["symmetry_operations"]),
            "same": lg_a["number"] == lg_b["number"],
        }
    return {
        "label_a": label_a, "label_b": label_b,
        "symbol_a": results_a["space_group_symbol"], "number_a": results_a["space_group_number"],
        "symbol_b": results_b["space_group_symbol"], "number_b": results_b["space_group_number"],
        "n_ops_a": len(results_a["symmetry_operations"]), "n_ops_b": len(results_b["symmetry_operations"]),
        "same_space_group": results_a["space_group_number"] == results_b["space_group_number"],
        "layer_comparison": layer_comparison,
    }

def _describe_comparison(cmp):
    if cmp["same_space_group"]:
        lines = [f"Symmetry PRESERVED between the two structures: both "
                f"{cmp['symbol_a']} (No. {cmp['number_a']}, {cmp['n_ops_a']} ops)."]
    else:
        direction = "gained" if cmp["n_ops_b"] > cmp["n_ops_a"] else "lost"
        lines = [f"Symmetry CHANGED between the two structures: "
                f"{cmp['symbol_a']} (No. {cmp['number_a']}, {cmp['n_ops_a']} ops) -> "
                f"{cmp['symbol_b']} (No. {cmp['number_b']}, {cmp['n_ops_b']} ops) -- "
                f"{direction} symmetry operations."]
    lc = cmp["layer_comparison"]
    if lc is not None:
        if lc["same"]:
            lines.append(f"Layer group PRESERVED: both {lc['symbol_a']} "
                         f"(No. {lc['number_a']}, {lc['n_ops_a']} ops).")
        else:
            direction = "gained" if lc["n_ops_b"] > lc["n_ops_a"] else "lost"
            lines.append(f"Layer group CHANGED: {lc['symbol_a']} (No. {lc['number_a']}, "
                         f"{lc['n_ops_a']} ops) -> {lc['symbol_b']} (No. {lc['number_b']}, "
                         f"{lc['n_ops_b']} ops) -- {direction} symmetry operations.")
    return "\n".join(lines)

# --- Report formatting --------------------------------------------------
_WIDTH = 74

def _rule(char="-"):
    return char * _WIDTH

def format_report(results, source_file, fmt, show_operations=True, symprec_scan=None,
                  layer_symprec_scan=None, comparison=None):
    lat = results["lattice"]
    lines = []
    lines.append(_rule("="))
    lines.append("CRYSTAL SYMMETRY REPORT - STB Suite".center(_WIDTH))
    lines.append(_rule("="))
    lines.append("")
    lines.append(f"Generated        : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"Source file      : {source_file}  (format: {fmt})")
    lines.append(f"Symmetry precision: symprec={results['symprec']:g}, "
                 f"angle_tolerance={results['angle_tolerance']:g}°")

    n_vacuum = sum(results["vacuum_axes"])
    if n_vacuum:
        axis_names = [name for name, is_vacuum in zip("abc", results["vacuum_axes"]) if is_vacuum]
        lines.append("")
        lines.append(f"NOTE: vacuum gap (>= {VACUUM_GAP_ANG:g} Å) detected along axis/axes "
                     f"{', '.join(axis_names)} (e.g. a slab or wire). 3D space-group detection")
        lines.append("isn't physically meaningful for a vacuum-padded structure -- spglib treats")
        lines.append("the vacuum as ordinary empty space, not a broken periodicity, so SPACE GROUP")
        if results["layer_group"] is not None:
            lines.append("below describes the padded supercell, not the real 2D symmetry -- see")
            lines.append("LAYER GROUP for that (spglib's dedicated 2D detection, not a heuristic).")
        elif n_vacuum == 1:
            lines.append("below describes the padded supercell, not the real layer symmetry.")
            lines.append("Dedicated layer-group detection was attempted (exactly one vacuum axis)")
            lines.append("but failed -- the structure may not be genuinely 2D-periodic within")
            lines.append("--symprec, or the installed spglib predates get_layergroup() (>= 2.1.0).")
        else:
            lines.append("below describes the padded supercell, not the real rod/point symmetry.")
            lines.append("Dedicated layer-group detection needs exactly one vacuum axis (this")
            lines.append(f"structure has {n_vacuum}, e.g. a wire or an isolated molecule) so it")
            lines.append("wasn't attempted -- no equivalent rod/point-group detection is available.")

    lines.append("")
    lines.append(_rule())
    lines.append("SPACE GROUP")
    lines.append(_rule())
    lines.append(f"Space group      : {results['space_group_symbol']} (No. {results['space_group_number']})")
    lines.append(f"Hall symbol      : {results['hall_symbol']} (Hall No. {results['hall_number']})")
    lines.append(f"Point group      : {results['point_group']}")
    lines.append(f"Crystal system   : {results['crystal_system']}")
    lines.append(f"Lattice type     : {results['lattice_type']}")
    lines.append(f"Pearson symbol   : {results['pearson_symbol']}")
    if results["setting_choice"]:
        lines.append(f"Setting choice   : {results['setting_choice']}")
    lines.append(f"Symmetry operations: {len(results['symmetry_operations'])}")

    if results["layer_group"] is not None:
        lg = results["layer_group"]
        lines.append("")
        lines.append(_rule())
        lines.append(f"LAYER GROUP (aperiodic axis: {lg['aperiodic_axis']}) -- the physically correct")
        lines.append("symmetry for this vacuum-padded structure (spglib's dedicated layer-group")
        lines.append("detection, not a heuristic filter of the SPACE GROUP above)")
        lines.append(_rule())
        lines.append(f"Layer group      : {lg['symbol']} (No. {lg['number']})")
        lines.append(f"Hall symbol      : {lg['hall_symbol']} (Hall No. {lg['hall_number']})")
        lines.append(f"Point group      : {lg['point_group']}")
        lines.append(f"Symmetry operations: {len(lg['symmetry_operations'])}")
        lines.append("")
        lines.append(f"Symmetrically distinct sites: {lg['n_distinct_sites']}")
        lines.append(f"{'Wyckoff':<10}{'Species':<10}{'n atoms':<10}{'Site symmetry':<16}{'Example atom'}")
        for orbit in lg["orbits"]:
            lines.append(f"{orbit['wyckoff']:<10}{orbit['species']:<10}{orbit['n_atoms']:<10}"
                         f"{orbit['site_symmetry']:<16}{orbit['example_atom_id']}")
        lines.append("")
        lines.append(f"{'Atom':>4}  {'Sp.':<3}  {'Wyckoff':<8}  {'Site symmetry':<16}{'Orbit':<6}Distortion(Å)")
        for site in lg["sites"]:
            lines.append(f"{site['atom_id']:>4}  {site['species']:<3}  {site['wyckoff']:<8}  "
                         f"{site['site_symmetry']:<16}{site['orbit']:<6}{site['distortion']:.4f}")
        if show_operations:
            lines.append("")
            lines.append(f"Symmetry operations ({len(lg['symmetry_operations'])}), in x,y,z notation:")
            for i, op_str in enumerate(lg["symmetry_operations"]):
                lines.append(f"{i+1:>4}: {op_str}")

    lines.append("")
    lines.append(_rule())
    lines.append("LATTICE")
    lines.append(_rule())
    lines.append(f"a = {lat['a']:.3f} Å   b = {lat['b']:.3f} Å   c = {lat['c']:.3f} Å")
    lines.append(f"alpha = {lat['alpha']:.2f}°   beta = {lat['beta']:.2f}°   gamma = {lat['gamma']:.2f}°")
    lines.append(f"Volume = {lat['volume']:.3f} Å³")
    lines.append(f"Formula          : {lat['reduced_formula']} (Z = {lat['z']:g} formula units in this cell)")
    lines.append("")
    lines.append("Lattice vectors (Å):")
    for i, vec in enumerate(lat["vectors"]):
        lines.append(f"  a_{i+1}: {vec[0]:12.6f}  {vec[1]:12.6f}  {vec[2]:12.6f}")

    if symprec_scan is not None:
        lines.append("")
        lines.append(_rule())
        lines.append("SYMPREC SENSITIVITY SCAN")
        lines.append(_rule())
        lines.append(f"{'symprec':<12}{'space group'}")
        for symprec, symbol, number in symprec_scan:
            group_str = f"{symbol} (No. {number})" if number is not None else "not determined"
            lines.append(f"{symprec:<12g}{group_str}")
        lines.append("")
        lines.append(_describe_symprec_scan(symprec_scan))

    if layer_symprec_scan is not None:
        lines.append("")
        lines.append(_rule())
        lines.append("LAYER GROUP SYMPREC SENSITIVITY SCAN")
        lines.append(_rule())
        lines.append(f"{'symprec':<12}{'layer group'}")
        for symprec, symbol, number in layer_symprec_scan:
            group_str = f"{symbol} (No. {number})" if number is not None else "not determined"
            lines.append(f"{symprec:<12g}{group_str}")
        lines.append("")
        lines.append(_describe_symprec_scan(layer_symprec_scan, label="layer group"))

    lines.append("")
    lines.append(_rule())
    lines.append(f"SYMMETRICALLY DISTINCT SITES: {results['n_distinct_sites']}")
    lines.append(_rule())
    lines.append(f"{'Wyckoff':<10}{'Species':<10}{'n atoms':<10}{'Site symmetry':<16}{'Example atom'}")
    for orbit in results["orbits"]:
        lines.append(f"{orbit['wyckoff']:<10}{orbit['species']:<10}{orbit['n_atoms']:<10}"
                     f"{orbit['site_symmetry']:<16}{orbit['example_atom_id']}")

    lines.append("")
    lines.append(_rule())
    lines.append("ATOMIC SITES (fractional coordinates; Distortion = how far this atom sits")
    lines.append("from where the detected symmetry operations would exactly place it --")
    lines.append("bounded by roughly symprec, 0 means the operations map it exactly)")
    lines.append(_rule())
    lines.append(f"{'Atom':>4}  {'Sp.':<3}  {'Wyckoff':<8}  {'Frac. x':>10}  {'Frac. y':>10}  "
                 f"{'Frac. z':>10}  {'Orbit':<6}Distortion(Å)")
    for site in results["sites"]:
        fc = site["frac_coords"]
        lines.append(f"{site['atom_id']:>4}  {site['species']:<3}  {site['wyckoff']:<8}  "
                     f"{fc[0]:10.6f}  {fc[1]:10.6f}  {fc[2]:10.6f}  {site['orbit']:<6}"
                     f"{site['distortion']:.4f}")

    if show_operations:
        lines.append("")
        lines.append(_rule())
        lines.append(f"SYMMETRY OPERATIONS ({len(results['symmetry_operations'])}), in x,y,z notation")
        lines.append(_rule())
        for i, op_str in enumerate(results["symmetry_operations"]):
            lines.append(f"{i+1:>4}: {op_str}")

    if comparison is not None:
        lines.append("")
        lines.append(_rule())
        lines.append("SYMMETRY COMPARISON")
        lines.append(_rule())
        sg_a = f"{comparison['symbol_a']} (No. {comparison['number_a']})"
        sg_b = f"{comparison['symbol_b']} (No. {comparison['number_b']})"
        lines.append(f"{'':<20}{comparison['label_a']:<28}{comparison['label_b']}")
        lines.append(f"{'Space group':<20}{sg_a:<28}{sg_b}")
        lines.append(f"{'Symmetry operations':<20}{comparison['n_ops_a']:<28}{comparison['n_ops_b']}")
        lc = comparison["layer_comparison"]
        if lc is not None:
            lg_a = f"{lc['symbol_a']} (No. {lc['number_a']})"
            lg_b = f"{lc['symbol_b']} (No. {lc['number_b']})"
            lines.append("")
            lines.append(f"{'Layer group':<20}{lg_a:<28}{lg_b}")
            lines.append(f"{'Layer group ops':<20}{lc['n_ops_a']:<28}{lc['n_ops_b']}")
        lines.append("")
        lines.append(_describe_comparison(comparison))

    lines.append(_rule("="))
    return "\n".join(lines) + "\n"

def main():
    parser = argparse.ArgumentParser(
        description="Analyze crystal symmetry (space group, Wyckoff positions, symmetry "
                     "operations) from a SIESTA structure file.",
        epilog="Example usage:\n"
               "  stb-symmetry --file structure.fdf --format fdf\n"
               "  stb-symmetry --file siesta.STRUCT_OUT --format struct_out --symprec 1e-4",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--file", required=True, help="Path to structure file.")
    parser.add_argument("--format", required=True, choices=["fdf", "struct_out"],
                        help="Input file format:\n"
                             "  fdf:        SIESTA structure input (%%block LatticeVectors etc.)\n"
                             "  struct_out: SIESTA post-relaxation output (.STRUCT_OUT)")
    parser.add_argument("--symprec", type=float, default=1e-3,
                        help="Symmetry-detection tolerance in Å (default: 1e-3).")
    parser.add_argument("--angle-tolerance", type=float, default=5.0,
                        help="Symmetry-detection angle tolerance in degrees (default: 5.0).")
    parser.add_argument("--no-operations", dest="show_operations", action="store_false",
                        help="Skip the full symmetry-operations list in the report -- keeps "
                             "just the count in the SPACE GROUP section. Useful for "
                             "high-symmetry structures where the full list is long.")
    parser.add_argument("--scan-symprec", action="store_true",
                        help="Also detect symmetry at a range of tolerances "
                             f"({', '.join(f'{v:g}' for v in DEFAULT_SYMPREC_SCAN)}) and report "
                             "where the space group changes -- reveals symmetry that's present "
                             "but hidden by small numerical noise (e.g. from a DFT relaxation) "
                             "at the --symprec you asked for.")
    parser.add_argument("--compare-to", type=str, metavar="PATH",
                        help="Also analyze a second structure file and report whether it "
                             "shares the same space group -- e.g. compare a pre-relaxation "
                             ".fdf against its post-relaxation .STRUCT_OUT to check whether "
                             "symmetry survived. Analyzed with the same --symprec/"
                             "--angle-tolerance as --file.")
    parser.add_argument("--compare-format", choices=["fdf", "struct_out"],
                        help="File format of --compare-to (required if --compare-to is given).")
    parser.add_argument("--write-refined", type=str, metavar="PATH",
                        help="Also write the symmetry-refined structure (conventional cell, "
                             "positions snapped to the detected symmetry) to this .fdf path "
                             "-- the same reduction stb-unitcell --mode refined uses.")
    parser.add_argument("-o", "--output-dir", type=str, default=".",
                        help="Directory to write symmetry.dat into (default: current "
                             "directory). Created if it doesn't exist.")
    parser.add_argument("-v", "--version", action="version",
                        version=f"stb-symmetry {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    args = parser.parse_args()

    if args.symprec <= 0:
        parser.error("--symprec must be positive.")
    if args.angle_tolerance <= 0:
        parser.error("--angle-tolerance must be positive.")
    if args.compare_to and not args.compare_format:
        parser.error("--compare-format is required when --compare-to is given.")

    os.makedirs(args.output_dir, exist_ok=True)

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("CRYSTAL SYMMETRY ANALYSIS:", 'bold'))
    print("-"*60)

    print("\n[INFO] Reading structure file...")
    try:
        structure = structure_io.read_siesta_structure(args.file, args.format)
    except FileNotFoundError:
        print(color_text(f"[ERROR] Structure file '{args.file}' not found.", 'red'))
        sys.exit(1)
    except ValueError as e:
        print(color_text(f"[ERROR] {e}", 'red'))
        sys.exit(1)

    print("[INFO] Detecting symmetry...")
    try:
        results = compute_symmetry(structure, symprec=args.symprec, angle_tolerance=args.angle_tolerance)
    except Exception as e:
        print(color_text(f"[ERROR] Symmetry detection failed: {e}", 'red'))
        sys.exit(1)

    symprec_scan = None
    layer_symprec_scan = None
    if args.scan_symprec:
        print("[INFO] Scanning symprec sensitivity...")
        try:
            symprec_scan = scan_symprec(structure, angle_tolerance=args.angle_tolerance)
        except Exception as e:
            print(color_text(f"[WARNING] --scan-symprec failed: {e}", 'yellow'))
        if sum(results["vacuum_axes"]) == 1:
            print("[INFO] Scanning layer-group symprec sensitivity...")
            try:
                layer_symprec_scan = scan_symprec_layer(structure, results["vacuum_axes"].index(True))
            except Exception as e:
                print(color_text(f"[WARNING] --scan-symprec (layer group) failed: {e}", 'yellow'))

    comparison = None
    if args.compare_to:
        print("\n[INFO] Reading comparison structure file...")
        try:
            structure2 = structure_io.read_siesta_structure(args.compare_to, args.compare_format)
        except FileNotFoundError:
            print(color_text(f"[ERROR] Comparison structure file '{args.compare_to}' not found.", 'red'))
            sys.exit(1)
        except ValueError as e:
            print(color_text(f"[ERROR] {e}", 'red'))
            sys.exit(1)

        print("[INFO] Detecting symmetry of comparison structure...")
        try:
            results2 = compute_symmetry(structure2, symprec=args.symprec, angle_tolerance=args.angle_tolerance)
        except Exception as e:
            print(color_text(f"[ERROR] Symmetry detection failed for comparison structure: {e}", 'red'))
            sys.exit(1)
        comparison = compare_symmetry(results, results2, args.file, args.compare_to)

    if args.write_refined:
        n_vacuum = sum(results["vacuum_axes"])
        if n_vacuum == 1 and results["layer_group"] is not None:
            # Refine against the layer group instead of the (for a
            # vacuum-padded structure, not physically meaningful) 3D space
            # group -- uses the same std_lattice/std_positions/std_types
            # convention as an ordinary SpglibDataset, just standardized
            # against the 2D symmetry instead.
            print(f"\n[INFO] Writing layer-group-refined structure to {args.write_refined}...")
            try:
                from pymatgen.core import Element, Lattice, Structure
                lg = results["layer_group"]
                species = [Element.from_Z(int(z)).symbol for z in lg["std_types"]]
                refined_structure = Structure(Lattice(lg["std_lattice"]), species, lg["std_positions"])
                structure_io.write_fdf(structure_io.from_pymatgen(refined_structure), args.write_refined)
            except Exception as e:
                print(color_text(f"[WARNING] --write-refined failed: {e}", 'yellow'))
        elif n_vacuum > 0:
            print(color_text(
                "\n[ERROR] --write-refined skipped: refining needs either a genuinely 3D bulk "
                "structure (0 vacuum axes) or a successfully-detected layer group (exactly 1 "
                "vacuum axis) to standardize against -- this structure has "
                f"{n_vacuum} vacuum axis/axes and " +
                ("no layer group was determined" if n_vacuum == 1 else
                 "no equivalent rod/point-group refinement is available") +
                ". Refining against the 3D space group instead could reshape the cell/vacuum "
                "axis into something that no longer represents the intended structure.", 'red'))
        else:
            print(f"\n[INFO] Writing symmetry-refined structure to {args.write_refined}...")
            try:
                from stb.core.symmetry import reduce_to_unitcell
                refined_structure, _ = reduce_to_unitcell(
                    structure, "refined", symprec=args.symprec, angle_tolerance=args.angle_tolerance)
                structure_io.write_fdf(structure_io.from_pymatgen(refined_structure), args.write_refined)
            except Exception as e:
                print(color_text(f"[WARNING] --write-refined failed: {e}", 'yellow'))

    report = format_report(results, args.file, args.format, show_operations=args.show_operations,
                           symprec_scan=symprec_scan, layer_symprec_scan=layer_symprec_scan,
                           comparison=comparison)
    print("\n" + report)

    out_path = os.path.join(args.output_dir, "symmetry.dat")
    with open(out_path, "w") as f:
        f.write(report)

    print(f"[INFO] Job complete! Results saved to {out_path}")
    print("-"*60)

if __name__ == "__main__":
    main()
