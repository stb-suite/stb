#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "2.0.0"

import xml.etree.ElementTree as ET
import numpy as np
import os
import pandas as pd
import argparse
import sys
from datetime import datetime
import matplotlib.pyplot as plt

from stb.core import citations
from stb.core.cli import color_text, show_intro, print_dual, print_section, print_table
from stb.core.orbitals import ORBITAL_MAP, ORBITAL_ORDER, get_orbital_name, get_detailed_orbital_name
from stb.core.siesta_bands import read_data as read_bands_data, read_eig_mesh, cbm_vbm

REPORT_FILE = "stb_dos_report.txt"
BIB_FILE = "references.bib"


def parse_data_string(data_str):
    """
    Parses a space/newline-separated string of numbers into a numpy array.
    """
    if data_str is None:
        return np.array([])
    try:
        data = np.array([float(val) for val in data_str.strip().split()])
        return data
    except Exception as e:
        return np.array([])


def resolve_label_and_dir(filename, label_arg):
    """Returns (label, dirname) used to locate a companion <label>.bands/
    <label>.EIG next to `filename`, for --shift vbm/cbm. If `label_arg` is
    given, it's used directly; otherwise the label is derived from
    `filename`'s own basename (stripping a trailing '.PDOS.xml' or '.xml'),
    so --label is never strictly required when the PDOS.xml already
    follows SIESTA's own <label>.PDOS.xml naming.
    """
    dirname = os.path.dirname(filename) or "."
    if label_arg:
        return label_arg, dirname
    base = os.path.basename(filename)
    if base.endswith(".PDOS.xml"):
        label = base[:-len(".PDOS.xml")]
    elif base.endswith(".xml"):
        label = base[:-len(".xml")]
    else:
        label = base
    return label, dirname


def find_vbm_cbm_source(label, dirname):
    """Returns ('bands', path) or ('eig', path) or (None, None) -- the
    hierarchy requested for --shift vbm/cbm: <label>.bands first (the same
    k-path file stb-bands itself defaults to), <label>.EIG only if no
    .bands is present (the fuller SCF k-mesh, same fallback stb-bands
    --eig-file already relies on for the cases a k-path alone misses).
    """
    bands_path = os.path.join(dirname, f"{label}.bands")
    if os.path.isfile(bands_path):
        return "bands", bands_path
    eig_path = os.path.join(dirname, f"{label}.EIG")
    if os.path.isfile(eig_path):
        return "eig", eig_path
    return None, None


def vbm_cbm_from_source(source_kind, source_path):
    """Returns (fermi_source, nspin_source, vbm, cbm) for a companion
    <label>.bands or <label>.EIG file, via the exact same
    core.siesta_bands parsing/physics (read_data/read_eig_mesh + cbm_vbm)
    stb-bands itself uses -- so the VBM/CBM value used here to shift the
    DOS plot always matches what `stb-bands --shift vbm/cbm` would report
    for the same file, not a second, potentially-diverging computation.
    """
    if source_kind == "bands":
        fermi_source, _high_sym, dic_bands, nspin_source = read_bands_data(source_path)
    else:
        fermi_source, dic_bands, nspin_source, _kpoints = read_eig_mesh(source_path)
    result = cbm_vbm(fermi_source, dic_bands, nspin_source)
    vbm, cbm = result["combined"][0], result["combined"][1]
    return fermi_source, nspin_source, vbm, cbm


def estimate_vbm_cbm_from_dos(energies, dos_total, fermi_energy, threshold_frac):
    """Heuristic VBM/CBM estimate directly from the (broadened) total DOS
    curve, for when no companion .bands/.EIG file is available: walks
    outward from the Fermi energy on each side until the DOS first rises
    above `threshold_frac * max(dos_total)` -- the point where the curve
    leaves the near-zero gap region. Deliberately less precise than a
    .bands/.EIG-derived VBM/CBM (the true band edge is blurred by the
    PDOS's own energy broadening/smearing), which is why this is only ever
    used as an explicit, opt-in (--estimate-from-dos) fallback -- never
    silently substituted for the real thing.

    `energies` must be ascending (SIESTA's own PDOS.xml energy grid
    convention); `dos_total` is the DOS summed over every atom, orbital,
    and spin channel, same length as `energies`.
    """
    energies = np.asarray(energies)
    dos_total = np.asarray(dos_total)
    threshold = threshold_frac * dos_total.max()
    idx_ef = int(np.searchsorted(energies, fermi_energy))

    vbm = None
    for i in range(idx_ef - 1, -1, -1):
        if dos_total[i] > threshold:
            vbm = energies[i]
            break
    cbm = None
    for i in range(idx_ef, len(energies)):
        if dos_total[i] > threshold:
            cbm = energies[i]
            break

    if vbm is None or cbm is None:
        side = ("either side" if vbm is None and cbm is None
                else "the occupied side" if vbm is None else "the empty side")
        raise ValueError(
            f"Could not estimate a VBM/CBM from the DOS: no energy point with DOS "
            f"above {threshold_frac} * max(DOS) was found on {side} of the Fermi "
            "energy. Try a smaller --dos-threshold-frac."
        )
    return vbm, cbm


# --- Parsing (pure -- no printing, raises on fatal errors) -------------

def parse_pdos_xml(input_file, projection_mode):
    """Parses a PDOS.xml file into a plain dict, with no side effects
    (printing/file I/O) so the caller can decide what to report and when
    -- fatal problems raise FileNotFoundError/ET.ParseError/ValueError
    instead of calling sys.exit() directly, mirroring bands.py's own
    "read + validate everything before printing/writing anything" order.
    Non-fatal issues (missing <fermi_energy>/<nspin>, a data-length
    mismatch on one orbital, ...) are collected as (message, color) pairs
    in the returned 'messages' list instead of being printed immediately,
    so main() can render them under a single [1] INPUT DATA section.
    """
    messages = []

    tree = ET.parse(input_file)  # raises ET.ParseError / FileNotFoundError
    root = tree.getroot()

    # --- Fermi energy (for the automatic shift) ---
    fermi_energy_element = root.find('fermi_energy')
    e_fermi = 0.0
    if fermi_energy_element is not None:
        e_fermi = float(fermi_energy_element.text.strip())
    else:
        messages.append(("<fermi_energy> tag not found. Using 0.0 eV as default Fermi level.", 'yellow'))

    # nspin: each <orbital>'s <data> holds nspin values per energy point,
    # interleaved (e0-spin1, e0-spin2, e1-spin1, ...), not one block of
    # num_energy_points values -- matches sisl's own reader
    # (pdosSileSiesta.read_data(): "DOS = ...reshape(-1, nspin)"). A fixed
    # nspin=1 assumption here silently drops all PDOS data on any
    # spin-polarized/non-collinear file (every orbital's data length
    # would mismatch num_energy_points and get skipped).
    nspin_element = root.find('nspin')
    nspin = 1
    if nspin_element is not None:
        try:
            nspin = int(nspin_element.text.strip())
        except (TypeError, ValueError):
            messages.append(("<nspin> tag could not be parsed. Assuming nspin=1.", 'yellow'))
    else:
        messages.append(("<nspin> tag not found. Assuming nspin=1.", 'yellow'))

    if nspin == 1:
        spin_suffixes = ['']
    elif nspin == 2:
        spin_suffixes = ['_up', '_down']
        messages.append((f"Detected nspin=2 (spin-polarized); exporting {spin_suffixes} columns per orbital.", None))
    else:
        spin_suffixes = [f'_s{s + 1}' for s in range(nspin)]
        messages.append((f"nspin={nspin} (non-collinear/SOC) -- spin components exported "
                          f"as raw {spin_suffixes} without physical relabeling (SIESTA orders these "
                          "as total/Mx/My/Mz); consult your SIESTA documentation for the exact "
                          "convention.", 'yellow'))

    # --- Energy grid ---
    energy_values_element = root.find('energy_values')
    if energy_values_element is None:
        raise ValueError("<energy_values> tag not found. Cannot proceed.")

    energy_values = parse_data_string(energy_values_element.text.strip())
    num_energy_points = len(energy_values)
    if num_energy_points == 0:
        raise ValueError("No energy points found in <energy_values>. Cannot proceed.")

    # --- Orbital data ---
    all_orbital_tags = root.findall('orbital')
    if not all_orbital_tags:
        raise ValueError("No <orbital> tags found in the XML file.")

    atom_data = {}
    all_species = set()
    processed_atoms_count = 0
    skipped_lm_orbitals = 0
    skipped_l_values = set()

    for orbital in all_orbital_tags:
        try:
            atom_index = int(orbital.attrib.get('atom_index', -1))
            atom_species = orbital.attrib.get('species', 'Unknown')
            l_attr = orbital.attrib.get('l')
            m_val = int(orbital.attrib.get('m', 999))  # 999: invalid-m flag

            if atom_index == -1:
                messages.append(("Orbital found with no 'atom_index' attribute. Skipping.", 'yellow'))
                continue

            if l_attr is None:
                messages.append(("Orbital found with no 'l' attribute. Skipping.", 'yellow'))
                continue
            l_val = int(l_attr)

            orbital_name = None
            if projection_mode == 'l':
                orbital_name = get_orbital_name(l_val)
            elif projection_mode == 'ml':
                orbital_name = get_detailed_orbital_name(l_val, m_val)

            # Skip if orbital is not one we want to process (l > 3, i.e.
            # g-orbitals and beyond, or an m outside the (l,m) map).
            # Counted (not silently dropped) so the excluded fraction of
            # the basis is visible -- a summed 'total' DOS built only
            # from s/p/d/f will otherwise look complete while quietly
            # missing any g-orbital contribution.
            if orbital_name is None:
                skipped_lm_orbitals += 1
                skipped_l_values.add(l_val)
                continue

            if atom_index not in atom_data:
                atom_data[atom_index] = {'species': atom_species}
                all_species.add(atom_species)
                processed_atoms_count += 1

            data_element = orbital.find('data')
            data_text = data_element.text if data_element is not None else None

            raw_pdos_data = parse_data_string(data_text)
            expected_len = num_energy_points * nspin

            if len(raw_pdos_data) == expected_len:
                # nspin values per energy point, interleaved -- see the
                # nspin-detection comment above. Shape (ne, nspin) even
                # for nspin=1, so downstream aggregation/output code is
                # uniform regardless of spin polarization.
                orbital_pdos_data = raw_pdos_data.reshape(num_energy_points, nspin)
                if orbital_name not in atom_data[atom_index]:
                    atom_data[atom_index][orbital_name] = np.zeros((num_energy_points, nspin))
                atom_data[atom_index][orbital_name] += orbital_pdos_data
            else:
                messages.append((f"Data mismatch for atom {atom_index}, l={l_val}: expected "
                                  f"{expected_len} points (nspin={nspin}), found "
                                  f"{len(raw_pdos_data)}. Skipping orbital.", 'yellow'))

        except Exception as e:
            messages.append((f"Error processing orbital {orbital.attrib.get('index', 'N/A')}: {e}", 'yellow'))

    if not atom_data:
        raise ValueError("No valid atom data was processed.")

    if skipped_lm_orbitals:
        messages.append((f"skipped {skipped_lm_orbitals} orbital(s) with l={sorted(skipped_l_values)} "
                          "(l > 3, i.e. g-orbitals or beyond, are not supported) or an unrecognized m -- "
                          "excluded from all output DOS.", 'yellow'))

    return {
        'e_fermi': e_fermi,
        'nspin': nspin,
        'spin_suffixes': spin_suffixes,
        'energy_values': energy_values,
        'num_energy_points': num_energy_points,
        'num_orbital_tags': len(all_orbital_tags),
        'atom_data': atom_data,
        'all_species': all_species,
        'processed_atoms_count': processed_atoms_count,
        'messages': messages,
    }


def _orbital_array(atom_entry, orb_name, num_energy_points, nspin):
    # (num_energy_points, nspin) array, zeros if this atom has no data for
    # orb_name (e.g. an atom with no d-orbitals in the basis when other
    # atoms do).
    arr = atom_entry.get(orb_name)
    if arr is None:
        return np.zeros((num_energy_points, nspin))
    return arr


def build_output_columns(parsed):
    """Column layout (sorted_columns/spin_columns/header) plus the
    grand-total DOS (summed over every atom, orbital, and spin channel
    individually per spin_column) -- computed unconditionally so it's
    available both for --shift vbm/cbm's DOS-estimate fallback and for
    the always-possible --view matplotlib preview, not only when 'total'
    is one of the requested --type outputs.
    """
    atom_data = parsed['atom_data']
    num_energy_points = parsed['num_energy_points']
    nspin = parsed['nspin']
    spin_suffixes = parsed['spin_suffixes']

    all_orbital_names = set()
    for idx in atom_data:
        all_orbital_names.update(atom_data[idx].keys())
    all_orbital_names.discard('species')

    sorted_columns = [orb for orb in ORBITAL_ORDER if orb in all_orbital_names]
    for orb in sorted(all_orbital_names):
        if orb not in sorted_columns:
            sorted_columns.append(orb)

    # Expand each base orbital column (e.g. 'p') into one column per spin
    # channel (e.g. 'p_up', 'p_down') -- a no-op for nspin=1, where
    # spin_suffixes == [''].
    spin_columns = [f"{col}{suf}" for col in sorted_columns for suf in spin_suffixes]

    header_parts = [f"#{'Energy(eV)':<14}"]
    header_parts.extend([f'{col:<14}' for col in spin_columns])
    header_str = "\t".join(header_parts) + "\n"

    grand_total = {col: np.zeros(num_energy_points) for col in spin_columns}
    for atom_index in atom_data:
        for orb_name in sorted_columns:
            arr = _orbital_array(atom_data[atom_index], orb_name, num_energy_points, nspin)
            for s, suf in enumerate(spin_suffixes):
                grand_total[f"{orb_name}{suf}"] += arr[:, s]

    return {
        'sorted_columns': sorted_columns,
        'spin_columns': spin_columns,
        'all_df_columns': ['Energy(eV)'] + spin_columns,
        'header_str': header_str,
        'grand_total': grand_total,
    }


# --- Energy shift resolution -------------------------------------------

def resolve_shift(parsed, columns, shift_str, label, input_file,
                   allow_dos_estimate, dos_threshold_frac):
    """Resolves --shift into a numeric value, without printing anything --
    returns (shift_value, shift_desc, messages), where messages is a list
    of (text, color) pairs for the caller to render under [2] ENERGY
    SHIFT. Raises ValueError on a fatal, unresolvable request (e.g.
    --shift vbm/cbm with neither a companion .bands/.EIG nor
    --estimate-from-dos)."""
    e_fermi = parsed['e_fermi']
    nspin = parsed['nspin']
    messages = []

    shift_key = shift_str.strip().lower()
    if shift_key == 'fermi':
        shift_value = e_fermi
        messages.append((f"Using automatic Fermi energy shift: {shift_value} eV", None))
        return shift_value, 'fermi', messages

    if shift_key in ('vbm', 'cbm'):
        label_resolved, dirname = resolve_label_and_dir(input_file, label)
        source_kind, source_path = find_vbm_cbm_source(label_resolved, dirname)
        if source_kind is not None:
            fermi_source, nspin_source, vbm, cbm = vbm_cbm_from_source(source_kind, source_path)
            shift_value = vbm if shift_key == 'vbm' else cbm
            messages.append((f"Using {shift_key.upper()} shift from '{source_path}' "
                              f"({'k-path .bands' if source_kind == 'bands' else 'k-mesh .EIG'}): "
                              f"{shift_value:.6f} eV", None))
            if abs(fermi_source - e_fermi) > 1e-3:
                messages.append((f"Fermi energy mismatch between '{input_file}' ({e_fermi:.6f} eV) "
                                  f"and '{source_path}' ({fermi_source:.6f} eV); the two files may be "
                                  "from different calculations.", 'yellow'))
            if nspin_source != nspin:
                messages.append((f"nspin mismatch between '{input_file}' (nspin={nspin}) and "
                                  f"'{source_path}' (nspin={nspin_source}).", 'yellow'))
            return shift_value, shift_key, messages

        if allow_dos_estimate:
            grand_total_flat = np.sum(list(columns['grand_total'].values()), axis=0)
            try:
                vbm_est, cbm_est = estimate_vbm_cbm_from_dos(
                    parsed['energy_values'], grand_total_flat, e_fermi, dos_threshold_frac)
            except ValueError:
                raise
            shift_value = vbm_est if shift_key == 'vbm' else cbm_est
            messages.append((f"no '{label_resolved}.bands' or '{label_resolved}.EIG' found next to "
                              f"'{input_file}' -- estimating {shift_key.upper()} directly from the DOS "
                              f"(--estimate-from-dos): {shift_value:.6f} eV. This is an approximation, "
                              "blurred by the PDOS's own energy broadening -- prefer a companion "
                              ".bands/.EIG file when available.", 'yellow'))
            return shift_value, shift_key, messages

        raise ValueError(
            f"no '{label_resolved}.bands' or '{label_resolved}.EIG' found next to '{input_file}' -- "
            f"cannot compute an exact {shift_key.upper()}. Pass --estimate-from-dos to fall back to "
            "an approximation from the DOS itself, or provide a matching .bands/.EIG file (see --label)."
        )

    try:
        shift_value = float(shift_str)
    except ValueError:
        raise ValueError(f"Invalid shift value '{shift_str}'. Must be 'fermi', 'vbm', 'cbm', or a number.")
    messages.append((f"Using manual energy shift: {shift_value} eV", None))
    return shift_value, 'manual', messages


# --- Output writing ------------------------------------------------------

def _save_dat(path, energies_shifted, data_dict, all_df_columns, header_str):
    data_dict = dict(data_dict)
    data_dict['Energy(eV)'] = energies_shifted
    df = pd.DataFrame(data_dict)
    with open(path, 'w') as f:
        f.write(header_str)
    df.to_csv(path, sep='\t', index=False, header=False, mode='a',
              columns=all_df_columns, float_format='%14.6E')


def write_total_dos(output_dir, energies_shifted, columns):
    path = os.path.join(output_dir, "dos_total.dat")
    _save_dat(path, energies_shifted, columns['grand_total'], columns['all_df_columns'], columns['header_str'])
    return path


def write_per_atom_dos(output_dir, energies_shifted, parsed, columns):
    """Writes one .dat per atom, and also returns `curves` -- an ordered
    {label: {col_name: 1D array}} of each atom's OWN per-orbital/spin
    breakdown (the exact same columns its own .dat has, e.g. 's', 'p',
    'd' -- NOT summed together) -- used by the 'atom' category view/
    gnuplot below to overlay every atom's every orbital as its own curve,
    so e.g. a 3-orbital atom contributes 3 separate lines, not 1."""
    output_dir_atoms = os.path.join(output_dir, "dos_per_atom")
    os.makedirs(output_dir_atoms, exist_ok=True)
    atom_data = parsed['atom_data']
    num_energy_points = parsed['num_energy_points']
    nspin = parsed['nspin']
    spin_suffixes = parsed['spin_suffixes']

    paths = []
    curves = {}
    for atom_index in sorted(atom_data.keys()):
        species = atom_data[atom_index]['species']
        atom_dos_data = {}
        for col in columns['sorted_columns']:
            arr = _orbital_array(atom_data[atom_index], col, num_energy_points, nspin)
            for s, suf in enumerate(spin_suffixes):
                atom_dos_data[f"{col}{suf}"] = arr[:, s]
        label = f"{species}_{atom_index}"
        path = os.path.join(output_dir_atoms, f"{label}.dat")
        _save_dat(path, energies_shifted, atom_dos_data, columns['all_df_columns'], columns['header_str'])
        paths.append(path)
        curves[label] = atom_dos_data
    return paths, curves


def write_per_species_dos(output_dir, energies_shifted, parsed, columns):
    """Writes one .dat per species, and also returns `curves` -- an
    ordered {label: {col_name: 1D array}} of each species' OWN
    per-orbital/spin breakdown (summed over its atoms, but NOT over
    orbitals) -- same purpose as write_per_atom_dos's own `curves`, one
    level up (per species instead of per atom)."""
    output_dir_species = os.path.join(output_dir, "dos_per_species")
    os.makedirs(output_dir_species, exist_ok=True)
    atom_data = parsed['atom_data']
    all_species = parsed['all_species']
    num_energy_points = parsed['num_energy_points']
    nspin = parsed['nspin']
    spin_suffixes = parsed['spin_suffixes']

    species_dos = {name: {col: np.zeros(num_energy_points) for col in columns['spin_columns']}
                    for name in sorted(all_species)}
    for atom_index in atom_data:
        species = atom_data[atom_index]['species']
        if species in species_dos:
            for col in columns['sorted_columns']:
                arr = _orbital_array(atom_data[atom_index], col, num_energy_points, nspin)
                for s, suf in enumerate(spin_suffixes):
                    species_dos[species][f"{col}{suf}"] += arr[:, s]

    paths = []
    curves = {}
    for species_name in species_dos:
        path = os.path.join(output_dir_species, f"dos_{species_name}.dat")
        _save_dat(path, energies_shifted, species_dos[species_name], columns['all_df_columns'], columns['header_str'])
        paths.append(path)
        curves[species_name] = species_dos[species_name]
    return paths, curves


def write_dos_gnuplot(dat_path, all_df_columns):
    """Writes a generic gnuplot script next to `dat_path`: plots every
    non-energy column via gnuplot's own `columnheader(i)`, reading the
    column names straight from the .dat file's own '#Energy(eV) s  p ...'
    header line. Used for the 'total' category, where there is only ever
    one .dat file, so the per-orbital breakdown (s, p, d, ...) is itself
    the most useful single-file plot -- see write_category_gnuplot for
    the 'atom'/'species' categories, which instead overlay one collapsed
    curve per file."""
    ncols = len(all_df_columns)
    dat_dir = os.path.dirname(dat_path)
    dat_base = os.path.basename(dat_path)
    stem = os.path.splitext(dat_base)[0]
    gplot_path = os.path.join(dat_dir, f"{stem}.gplot")
    lines = [
        'set terminal pdfcairo enhanced font "Arial,14" size 8,6\n',
        f'set output "{stem}.pdf"\n',
        'set xlabel "Energy (eV)"\n',
        'set ylabel "DOS (states/eV)"\n',
        'set key outside\n',
        'set grid\n',
        f'plot for [i=2:{ncols}] "{dat_base}" using 1:i with lines lw 2 title columnheader(i)\n',
    ]
    with open(gplot_path, 'w') as f:
        f.writelines(lines)
    return gplot_path


def write_category_gnuplot(dat_paths, labels, spin_columns, gplot_path):
    """Writes ONE gnuplot script for a whole category (dos_per_atom/ or
    dos_per_species/) -- one curve per (file, orbital/spin column) pair,
    e.g. a 3-orbital atom contributes 3 separate curves ('Sn_1: s',
    'Sn_1: p', 'Sn_1: d'), NOT one curve summed across its orbitals.
    Column positions are looked up by index (every .dat in this suite
    shares the same `spin_columns` order, column 1 is Energy) rather than
    gnuplot's own columnheader(), so each curve's title can combine both
    the file's label and that column's name. Placed inside the category's
    own folder, referencing each .dat by its bare basename (assumes
    gnuplot is run from that same folder, the same convention
    write_dos_gnuplot/bands.py's plot_gnuplot already use)."""
    stem = os.path.splitext(os.path.basename(gplot_path))[0]
    plot_entries = ", \\\n     ".join(
        f'"{os.path.basename(p)}" using 1:{col_idx} with lines lw 2 title "{label}: {col_name}"'
        for p, label in zip(dat_paths, labels)
        for col_idx, col_name in enumerate(spin_columns, start=2)
    )
    lines = [
        'set terminal pdfcairo enhanced font "Arial,14" size 8,6\n',
        f'set output "{stem}.pdf"\n',
        'set xlabel "Energy (eV)"\n',
        'set ylabel "DOS (states/eV)"\n',
        'set key outside\n',
        'set grid\n',
        f'plot {plot_entries}\n',
    ]
    with open(gplot_path, 'w') as f:
        f.writelines(lines)
    return gplot_path


def plot_category_matplotlib(energies_shifted, curves, title):
    """One matplotlib figure, one curve per (label, 1D array) in `curves`
    -- shared by the 'total' category (one curve per orbital/spin column
    of the grand total) and the 'atom'/'species' categories (one curve
    per (entity, orbital/spin column) pair, e.g. 'Sn_1: s'/'Sn_1: p' --
    never summed across orbitals). Does NOT call plt.show() itself --
    main() opens one figure per active category, then shows them all at
    once at the very end, after every report line and file has already
    been written, so a blocking GUI window never delays or hides them
    (same ordering as bands.py's own interactive preview)."""
    plt.figure(figsize=(8, 6))
    for label, curve in curves.items():
        plt.plot(energies_shifted, curve, label=label)
    plt.axvline(x=0.0, color='gray', linestyle='--', linewidth=1)
    plt.xlabel("Energy (eV)")
    plt.ylabel("DOS (states/eV)")
    plt.title(title)
    plt.legend()
    plt.grid(True)


def main():

    parser = argparse.ArgumentParser(
        description="Parse a PDOS.xml file and generate Gnuplot-ready .dat files.",
        epilog="Example usage:\n"
               "  stb-dos siesta.PDOS.xml --type total species --projection ml\n"
               "Put the filename first: --type takes one or more values, so anything\n"
               "written after it (including the filename) is consumed as another type.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "filename",
        type=str,
        nargs='?',
        default=None,
        help="The input .PDOS.xml file to process. Give this before --type (see epilog). "
             "Optional if --label is given (resolves to '<label>.PDOS.xml')."
    )

    parser.add_argument(
        "--label",
        type=str,
        default=None,
        help="SIESTA label. Shorthand for filename='<label>.PDOS.xml' if filename isn't "
             "given directly, and (with --shift vbm/cbm) used to locate a companion "
             "'<label>.bands'/'<label>.EIG' next to it. One of filename or --label is "
             "required. If filename is given explicitly without --label, the label is "
             "still derived automatically from its own name for --shift vbm/cbm."
    )

    parser.add_argument(
        "--type",
        nargs='+',
        choices=['total', 'atom', 'species'],
        default=['total', 'atom', 'species'],
        help="Type(s) of DOS to output.\n"
             "  total:   Sum of all atoms.\n"
             "  atom:    One file for each atom.\n"
             "  species: One file for each chemical species (e.g., C, N, B).\n"
             "You can select multiple, e.g., --type total species (default: all three)"
    )

    parser.add_argument(
        "--shift",
        type=str,
        default='fermi',
        help="Energy shift to apply. \n"
             "  'fermi': Automatically shift by the Fermi energy (default).\n"
             "  'vbm'/'cbm': Shift by the Valence Band Maximum / Conduction Band\n"
             "    Minimum, computed from a companion '<label>.bands' (preferred) or\n"
             "    '<label>.EIG' file -- see --label. Falls back to an approximate\n"
             "    estimate from the DOS itself only with --estimate-from-dos.\n"
             "  '0.0':   Use an absolute energy scale (no shift).\n"
             "  '-1.23': Apply a manual shift of -1.23 eV."
    )

    parser.add_argument(
        "--estimate-from-dos", action="store_true",
        help="With --shift vbm/cbm: if no companion '<label>.bands'/'<label>.EIG' is "
             "found, fall back to an approximate VBM/CBM estimated directly from the "
             "DOS itself (blurred by its own energy broadening -- less precise than a "
             "real .bands/.EIG). Off by default: without it, a missing .bands/.EIG is "
             "a clear error, not a silent approximation."
    )

    parser.add_argument(
        "--dos-threshold-frac", type=float, default=0.01,
        help="Only used by --shift vbm/cbm's --estimate-from-dos fallback: fraction of "
             "the peak total DOS above which the curve is considered 'occupied'/'empty' "
             "again past the gap region (default: 0.01)."
    )

    parser.add_argument(
        "--projection",
        type=str,
        choices=['l', 'ml'],
        default='l',
        help="Orbital projection detail level.\n"
             "  l:  Aggregate by angular momentum (s, p, d, f). (default)\n"
             "  ml: Project by magnetic quantum number (s, px, py, pz, dxy, etc.)."
    )

    parser.add_argument("-o", "--output-dir", type=str, default=".",
                        help="Directory to write dos_total.dat, dos_per_atom/, and "
                             "dos_per_species/ into (default: current directory). "
                             "Created if it doesn't exist.")

    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the full run report to {REPORT_FILE}. Off by default.")

    parser.add_argument("--save-gnuplot", action="store_true",
                        help="Also write a .gplot gnuplot script next to every .dat file "
                             "this run generates (dos_total.dat, each dos_per_atom/*.dat, "
                             "each dos_per_species/*.dat). Off by default.")

    parser.add_argument("--view", action="store_true",
                        help="Show an interactive matplotlib preview of the total DOS "
                             "before exiting. Off by default.")

    parser.add_argument("-v", "--version", action="version",
                        version=f"stb-dos {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if not args.filename and not args.label:
        parser.error("one of filename or --label is required.")
    if not args.filename:
        args.filename = f"{args.label}.PDOS.xml"

    if args.intro == True:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    try:
        parsed = parse_pdos_xml(args.filename, args.projection)
        columns = build_output_columns(parsed)
        shift_value, shift_key, shift_messages = resolve_shift(
            parsed, columns, args.shift, args.label, args.filename,
            args.estimate_from_dos, args.dos_threshold_frac)
    except ET.ParseError as e:
        print(f"Error parsing XML file '{args.filename}': {e}", file=sys.stderr)
        print("The file might be corrupted or not well-formed XML.", file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError:
        print(f"Error: File not found at '{args.filename}'", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)

    energies_shifted = parsed['energy_values'] - shift_value

    try:
        os.makedirs(args.output_dir, exist_ok=True)
    except OSError as e:
        print(f"Error: Could not create output directory '{args.output_dir}': {e}", file=sys.stderr)
        sys.exit(1)

    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    print_dual(color_text("===== STB-DOS REPORT =====", 'magenta'), f_out)

    label_resolved, _dirname = resolve_label_and_dir(args.filename, args.label)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time      : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Input file     : {args.filename}", f_out)
    print_dual(f"Label          : {label_resolved} (used for .bands/.EIG lookup with --shift vbm/cbm)", f_out)
    print_dual(f"DOS type(s)    : {', '.join(args.type)}", f_out)
    print_dual(f"Projection     : {args.projection}", f_out)
    print_dual(f"Shift mode     : {args.shift}", f_out)
    print_dual(f"Output dir     : {args.output_dir}", f_out)

    nspin_desc = {1: "non-polarized", 2: "spin-polarized"}.get(parsed['nspin'], "non-collinear/SOC")
    print_section("[1] INPUT DATA", f_out)
    print_table(["Quantity", "Value"], [
        (["Fermi energy", f"{parsed['e_fermi']:.6f} eV"], None),
        (["Spin channels", f"{parsed['nspin']} ({nspin_desc})"], None),
        (["Energy points", f"{parsed['num_energy_points']} "
          f"({parsed['energy_values'].min():.3f} to {parsed['energy_values'].max():.3f} eV)"], None),
        (["Orbital tags", f"{parsed['num_orbital_tags']}"], None),
        (["Atoms processed", f"{parsed['processed_atoms_count']}"], None),
        (["Species found", f"{sorted(parsed['all_species'])}"], None),
    ], f_out)
    for text, color in parsed['messages']:
        print_dual(color_text(f"[WARNING] {text}", color) if color else text, f_out)

    print_section("[2] ENERGY SHIFT", f_out)
    for text, color in shift_messages:
        print_dual(color_text(f"[WARNING] {text}", color) if color else text, f_out)

    print_section("[3] WRITING OUTPUT FILES", f_out)
    written_files = []
    gplot_files = []
    category_curves = {}  # category -> {label: 1D array}, for --view

    if 'total' in args.type:
        path = write_total_dos(args.output_dir, energies_shifted, columns)
        written_files.append(path)
        print_dual(color_text(f"[OK] Saved Total DOS to {path}", 'green'), f_out)
        category_curves['total'] = columns['grand_total']
        if args.save_gnuplot:
            gplot_path = write_dos_gnuplot(path, columns['all_df_columns'])
            gplot_files.append(gplot_path)
            print_dual(color_text(f"[OK] Gnuplot script written to '{gplot_path}'", 'green'), f_out)

    if 'atom' in args.type:
        atom_paths, atom_curves = write_per_atom_dos(args.output_dir, energies_shifted, parsed, columns)
        written_files.extend(atom_paths)
        output_dir_atoms = os.path.join(args.output_dir, "dos_per_atom")
        print_dual(color_text(f"[OK] Saved DOS per atom to '{output_dir_atoms}' directory "
                               f"({len(atom_paths)} files).", 'green'), f_out)
        category_curves['atom'] = {
            f"{label}: {col}": arr for label, cols in atom_curves.items() for col, arr in cols.items()
        }
        if args.save_gnuplot:
            gplot_path = os.path.join(output_dir_atoms, "dos_per_atom.gplot")
            write_category_gnuplot(atom_paths, list(atom_curves.keys()), columns['spin_columns'], gplot_path)
            gplot_files.append(gplot_path)
            print_dual(color_text(f"[OK] Gnuplot script (one curve per atom/orbital) written to "
                                   f"'{gplot_path}'.", 'green'), f_out)

    if 'species' in args.type:
        species_paths, species_curves = write_per_species_dos(args.output_dir, energies_shifted, parsed, columns)
        written_files.extend(species_paths)
        output_dir_species = os.path.join(args.output_dir, "dos_per_species")
        print_dual(color_text(f"[OK] Saved DOS per species to '{output_dir_species}' directory "
                               f"({len(species_paths)} files).", 'green'), f_out)
        category_curves['species'] = {
            f"{label}: {col}": arr for label, cols in species_curves.items() for col, arr in cols.items()
        }
        if args.save_gnuplot:
            gplot_path = os.path.join(output_dir_species, "dos_per_species.gplot")
            write_category_gnuplot(species_paths, list(species_curves.keys()), columns['spin_columns'], gplot_path)
            gplot_files.append(gplot_path)
            print_dual(color_text(f"[OK] Gnuplot script (one curve per species/orbital) written to "
                                   f"'{gplot_path}'.", 'green'), f_out)

    print_section("[4] REFERENCES", f_out)
    bib_entries = [citations.SIESTA, citations.SIESTA_RECENT]
    citations.write_bib_file(os.path.join(args.output_dir, BIB_FILE), bib_entries)
    print_dual(color_text(
        f"[OK] Citations for the methods used in this run written to "
        f"'{os.path.join(args.output_dir, BIB_FILE)}' ({len(bib_entries)} entries).", 'green'), f_out)

    print_section("[5] SUMMARY & FILES", f_out)
    print_dual("Status         : OK", f_out)
    for path in written_files:
        print_dual(f"Data file      : {path}", f_out)
    for path in gplot_files:
        print_dual(f"Gnuplot script : {path}", f_out)
    print_dual(f"References     : {os.path.join(args.output_dir, BIB_FILE)}", f_out)
    if report_path:
        print_dual(f"Report         : {report_path}", f_out)

    if f_out:
        f_out.close()

    # Interactive matplotlib preview runs last, after every report/file
    # write above, so a blocking GUI window never delays or hides them.
    # One figure per category actually written (total/atom/species),
    # shown together with a single plt.show() call.
    if args.view:
        category_titles = {
            'total': "Total density of states",
            'atom': "Density of states per atom",
            'species': "Density of states per species",
        }
        for category, curves in category_curves.items():
            plot_category_matplotlib(energies_shifted, curves, category_titles[category])
        plt.show()


if __name__ == "__main__":
    main()
