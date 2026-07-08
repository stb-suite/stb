#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################
#  Modified by "Parceiro de Programacao" (Gemini) #
#  to include (l,m) orbital projection support    #
#################################################

VERSION = "1.9.1" # Updated version

import xml.etree.ElementTree as ET
import numpy as np
import os
import pandas as pd
import argparse
import sys
from time import sleep


from stb.core.cli import COLORS, color_text, show_intro

# --- NEW: Map for (l,m) detailed projections ---
# This defines the standard SIESTA order for real spherical harmonics
ORBITAL_MAP = {
    0: {0: 's'},
    1: {-1: 'py', 0: 'pz', 1: 'px'},
    2: {-2: 'dxy', -1: 'dyz', 0: 'dz2', 1: 'dxz', 2: 'dx2-y2'},
    3: {-3: 'f-3', -2: 'f-2', -1: 'f-1', 0: 'f0', 1: 'f1', 2: 'f2', 3: 'f3'} # Using simple f names
}

# --- NEW: Sort order for output columns ---
ORBITAL_ORDER = [
    's', 
    'py', 'pz', 'px',
    'p', # For 'l' mode
    'dxy', 'dyz', 'dz2', 'dxz', 'dx2-y2',
    'd', # For 'l' mode
    'f-3', 'f-2', 'f-1', 'f0', 'f1', 'f2', 'f3',
    'f' # For 'l' mode
]



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
        print(f"Warning: Could not parse data string. Error: {e}", file=sys.stderr)
        return np.array([])

def get_orbital_name(l_val):
    """Maps angular momentum number 'l' to its name (s, p, d, f)."""
    l_map = {0: 's', 1: 'p', 2: 'd', 3: 'f'}
    return l_map.get(l_val, None) # Return None if not s,p,d,f

# --- NEW: Function to get detailed (l,m) orbital names ---
def get_detailed_orbital_name(l_val, m_val):
    """Maps (l, m) to orbital name (s, px, py, pz, dxy, ...)."""
    if l_val in ORBITAL_MAP:
        return ORBITAL_MAP[l_val].get(m_val, None) # Return None if m is invalid
    return None # Return None if l is invalid

# --- MODIFIED: Added 'projection_mode' argument ---
def process_pdos_xml(input_file, dos_types, shift_str, projection_mode):
    """
    Main function to parse the PDOS.xml file and generate output files.
    """
    try:
        tree = ET.parse(input_file)
        root = tree.getroot()

        # --- 1. Get Fermi Energy (for automatic shift) ---
        fermi_energy_element = root.find('fermi_energy')
        e_fermi = 0.0
        if fermi_energy_element is not None:
            e_fermi = float(fermi_energy_element.text.strip())
        else:
            print("Warning: <fermi_energy> tag not found. Using 0.0 eV as default Fermi level.", file=sys.stderr)

        # --- 1b. Get nspin: each <orbital>'s <data> holds nspin values per
        # energy point, interleaved (e0-spin1, e0-spin2, e1-spin1, ...), not
        # one block of num_energy_points values -- matches sisl's own reader
        # (pdosSileSiesta.read_data(): "DOS = ...reshape(-1, nspin)"). A
        # fixed nspin=1 assumption here silently drops all PDOS data on any
        # spin-polarized/non-collinear file (every orbital's data length
        # would mismatch num_energy_points and get skipped).
        nspin_element = root.find('nspin')
        nspin = 1
        if nspin_element is not None:
            try:
                nspin = int(nspin_element.text.strip())
            except (TypeError, ValueError):
                print("Warning: <nspin> tag could not be parsed. Assuming nspin=1.", file=sys.stderr)
        else:
            print("Warning: <nspin> tag not found. Assuming nspin=1.", file=sys.stderr)

        if nspin == 1:
            spin_suffixes = ['']
        elif nspin == 2:
            spin_suffixes = ['_up', '_down']
            print(f"Detected nspin=2 (spin-polarized); exporting {spin_suffixes} columns per orbital.")
        else:
            spin_suffixes = [f'_s{s + 1}' for s in range(nspin)]
            print(f"Warning: nspin={nspin} (non-collinear/SOC) -- spin components exported "
                  f"as raw {spin_suffixes} without physical relabeling (SIESTA orders these "
                  "as total/Mx/My/Mz); consult your SIESTA documentation for the exact "
                  "convention.", file=sys.stderr)

        # --- 2. Determine Energy Shift ---
        shift_value = 0.0
        if shift_str.lower() == 'fermi':
            shift_value = e_fermi
            print(f"Using automatic Fermi energy shift: {shift_value} eV")
        else:
            try:
                shift_value = float(shift_str)
                print(f"Using manual energy shift: {shift_value} eV")
            except ValueError:
                print(f"Error: Invalid shift value '{shift_str}'. Must be 'fermi' or a number.", file=sys.stderr)
                sys.exit(1)

        # --- 3. Get Energy Values ---
        energy_values_element = root.find('energy_values')
        if energy_values_element is None:
            print("Error: <energy_values> tag not found. Cannot proceed.", file=sys.stderr)
            sys.exit(1)

        energies_str = energy_values_element.text.strip()
        energy_values = parse_data_string(energies_str)
        energies_shifted = energy_values - shift_value
        num_energy_points = len(energies_shifted)

        if num_energy_points == 0:
            print("Error: No energy points found. Cannot proceed.", file=sys.stderr)
            sys.exit(1)

        print(f"Found {num_energy_points} energy points.")

        # --- 4. Find and Process Orbital Data ---
        all_orbital_tags = root.findall('orbital')

        if not all_orbital_tags:
            print("Error: No <orbital> tags found in the XML file.", file=sys.stderr)
            sys.exit(1)

        print(f"Found {len(all_orbital_tags)} <orbital> tags to process...")
        print(f"Using orbital projection mode: '{projection_mode}'")


        atom_data = {}
        all_species = set()
        processed_atoms_count = 0
        skipped_lm_orbitals = 0
        skipped_l_values = set()

        # --- MODIFIED: This loop is updated for dynamic orbital handling ---
        for orbital in all_orbital_tags:
            try:
                atom_index = int(orbital.attrib.get('atom_index', -1))
                atom_species = orbital.attrib.get('species', 'Unknown')
                l_val = int(orbital.attrib.get('l', -1))
                m_val = int(orbital.attrib.get('m', 999)) # Get m value, 999 as invalid flag

                if atom_index == -1:
                    print(f"Warning: Orbital found with no 'atom_index' attribute. Skipping.", file=sys.stderr)
                    continue

                # --- MODIFIED: Choose orbital name based on projection mode ---
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

                # --- MODIFIED: Initialize atom data dynamically ---
                if atom_index not in atom_data:
                    atom_data[atom_index] = {'species': atom_species}
                    all_species.add(atom_species)
                    processed_atoms_count += 1
                
                data_element = orbital.find('data')
                data_text = None
                if data_element is not None:
                    data_text = data_element.text

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
                    print(f"Warning: Data mismatch for atom {atom_index}, l={l_val}. Skipping orbital.", file=sys.stderr)
                    print(f"Expected {expected_len} points (nspin={nspin}), found {len(raw_pdos_data)}", file=sys.stderr)

            except Exception as e:
                print(f"Error processing orbital {orbital.attrib.get('index', 'N/A')}: {e}", file=sys.stderr)
        
        if not atom_data:
            print("Error: No valid atom data was processed.", file=sys.stderr)
            sys.exit(1)
            
        print(f"Successfully processed data for {processed_atoms_count} atoms.")
        print(f"Found species: {sorted(list(all_species))}")
        if skipped_lm_orbitals:
            print(f"Warning: skipped {skipped_lm_orbitals} orbital(s) with l={sorted(skipped_l_values)} "
                  "(l > 3, i.e. g-orbitals or beyond, are not supported) or an unrecognized m -- "
                  "excluded from all output DOS.", file=sys.stderr)

        # --- 5. Prepare and Write Output Data ---
        
        # --- MODIFIED: Output headers and columns are now DYNAMIC ---
        
        # 5a. Find all unique orbital columns that were processed
        all_orbital_names = set()
        for idx in atom_data:
            all_orbital_names.update(atom_data[idx].keys())
        all_orbital_names.remove('species') # Not a data column

        # 5b. Sort the columns for a clean output file
        sorted_columns = [orb for orb in ORBITAL_ORDER if orb in all_orbital_names]
        # Add any remaining orbitals not in the predefined order (just in case)
        for orb in sorted(list(all_orbital_names)):
            if orb not in sorted_columns:
                sorted_columns.append(orb)

        # 5b2. Expand each base orbital column (e.g. 'p') into one column per
        # spin channel (e.g. 'p_up', 'p_down') -- a no-op for nspin=1, where
        # spin_suffixes == [''].
        spin_columns = [f"{col}{suf}" for col in sorted_columns for suf in spin_suffixes]

        def _orbital_array(atom_entry, orb_name):
            # (num_energy_points, nspin) array, zeros if this atom has no
            # data for orb_name (e.g. an atom with no d-orbitals in the
            # basis when other atoms do).
            arr = atom_entry.get(orb_name)
            if arr is None:
                return np.zeros((num_energy_points, nspin))
            return arr

        print(f"Will generate files with columns: {['Energy(eV)'] + spin_columns}")

        # 5c. Create dynamic header and column list for pandas
        header_parts = [f"#{'Energy(eV)':<14}"]
        header_parts.extend([f'{col:<12}' for col in spin_columns])
        header_str = "\t".join(header_parts) + "\n"

        all_df_columns = ['Energy(eV)'] + spin_columns
        float_format_str = '%14.6E'

        # --- Mode 1: Total DOS ---
        if 'total' in dos_types:
            # Initialize a dictionary for total DOS
            total_dos = {col: np.zeros(num_energy_points) for col in spin_columns}

            for atom_index in atom_data:
                for orb_name in sorted_columns:
                    arr = _orbital_array(atom_data[atom_index], orb_name)
                    for s, suf in enumerate(spin_suffixes):
                        total_dos[f"{orb_name}{suf}"] += arr[:, s]

            total_dos['Energy(eV)'] = energies_shifted
            df_total = pd.DataFrame(total_dos)

            output_file_total = "dos_total.dat"
            with open(output_file_total, 'w') as f:
                f.write(header_str)
            df_total.to_csv(output_file_total, sep='\t', index=False, header=False, mode='a',
                            columns=all_df_columns, # Use dynamic columns
                            float_format=float_format_str)
            print(f"Saved Total DOS to {output_file_total}")

        # --- Mode 2: DOS per Atom ---
        if 'atom' in dos_types:
            output_dir_atoms = "dos_per_atom"
            if not os.path.exists(output_dir_atoms):
                os.makedirs(output_dir_atoms)

            for atom_index in sorted(atom_data.keys()):
                species = atom_data[atom_index]['species']

                # Build data for this atom's DataFrame
                atom_dos_data = {'Energy(eV)': energies_shifted}
                for col in sorted_columns:
                    arr = _orbital_array(atom_data[atom_index], col)
                    for s, suf in enumerate(spin_suffixes):
                        atom_dos_data[f"{col}{suf}"] = arr[:, s]

                df_atom = pd.DataFrame(atom_dos_data)

                output_file_atom = os.path.join(output_dir_atoms, f"{species}_{atom_index}.dat")
                with open(output_file_atom, 'w') as f:
                    f.write(header_str)
                df_atom.to_csv(output_file_atom, sep='\t', index=False, header=False, mode='a',
                               columns=all_df_columns, # Use dynamic columns
                               float_format=float_format_str)

            print(f"Saved DOS per atom to '{output_dir_atoms}' directory.")

        # --- Mode 3: DOS per Species ---
        if 'species' in dos_types:
            output_dir_species = "dos_per_species"
            if not os.path.exists(output_dir_species):
                os.makedirs(output_dir_species)

            species_dos = {}
            for species_name in sorted(list(all_species)):
                # Initialize dict for each species with all possible columns
                species_dos[species_name] = {col: np.zeros(num_energy_points) for col in spin_columns}

            for atom_index in atom_data:
                species = atom_data[atom_index]['species']
                if species in species_dos:
                    for col in sorted_columns:
                        arr = _orbital_array(atom_data[atom_index], col)
                        for s, suf in enumerate(spin_suffixes):
                            species_dos[species][f"{col}{suf}"] += arr[:, s]

            for species_name in species_dos:
                # Build DataFrame for this species
                species_dos_data = species_dos[species_name]
                species_dos_data['Energy(eV)'] = energies_shifted
                
                df_species = pd.DataFrame(species_dos_data)
                
                output_file_species = os.path.join(output_dir_species, f"dos_{species_name}.dat")
                with open(output_file_species, 'w') as f:
                    f.write(header_str)
                df_species.to_csv(output_file_species, sep='\t', index=False, header=False, mode='a',
                                  columns=all_df_columns, # Use dynamic columns
                                  float_format=float_format_str)
                
            print(f"Saved DOS per species to '{output_dir_species}' directory.")

    except ET.ParseError as e:
        print(f"Error parsing XML file '{input_file}': {e}", file=sys.stderr)
        print("The file might be corrupted or not well-formed XML.", file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError:
        print(f"Error: File not found at '{input_file}'", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)

def main():

    parser = argparse.ArgumentParser(
        description="Parse a PDOS.xml file and generate Gnuplot-ready .dat files.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument(
        "filename",
        type=str,
        help="The input .PDOS.xml file to process."
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
             "  '0.0':   Use an absolute energy scale (no shift).\n"
             "  '-1.23': Apply a manual shift of -1.23 eV."
    )

    # --- NEW: Added --projection argument ---
    parser.add_argument(
        "--projection",
        type=str,
        choices=['l', 'ml'],
        default='l',
        help="Orbital projection detail level.\n"
             "  l:  Aggregate by angular momentum (s, p, d, f). (default)\n"
             "  ml: Project by magnetic quantum number (s, px, py, pz, dxy, etc.)."
    )

    parser.add_argument("-v", "--version", action="version",
                        version=f"stb-dos {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()


    if args.intro == True:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("Density of States:", 'bold'))
    print("-"*60)
    
    # --- MODIFIED: Pass args.projection to the processing function ---
    process_pdos_xml(args.filename, args.type, args.shift, args.projection)
    
if __name__ == "__main__":
    main()
