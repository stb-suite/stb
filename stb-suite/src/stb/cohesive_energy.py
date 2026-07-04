#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################
    
VERSION = "1.9.1"    

import os
import sys
import re
import shutil
import argparse
from time import sleep
import numpy as np
from stb.core import structure_io, kspace
from stb.core.cli import COLORS, color_text, show_intro

# Base Template
CALC_TEMPLATE = """
## ===================================================================
## SYSTEM DEFINITION
## ===================================================================
SystemLabel      siesta
SystemName       siesta

%include structure.fdf

## ===================================================================
## BASIS SET DEFINITION
## ===================================================================
PAO.BasisType       split
PAO.BasisSize       DZP
PAO.EnergyShift     0.02 Ry

## ===================================================================
## K-POINT SAMPLING (BRILLOUIN ZONE)
## ===================================================================
kgrid.MonkhorstPack   [1  1  1]
Mesh.CutOff           320  Ry   
FilterCutoff          150  Ry

## ===================================================================
## EXCHANGE-CORRELATION (XC) FUNCTIONAL
## ===================================================================
XC.Functional       GGA
XC.Authors          PBE

## ===================================================================
## SPIN POLARIZATION
## ===================================================================
Spin                non-polarized

## ===================================================================
## SELF-CONSISTENT-FIELD (SCF)
## ===================================================================
MaxSCFIterations        300
SCF.Mixer.Weight        0.1
SCF.DM.Tolerance        1.0d-5  eV
SCF.Mixer.History       6
ElectronicTemperature   300 K
Diag.ParallelOverK      .true.

## ===================================================================
## STARTING THE CALCULATION
## ===================================================================
DM.UseSaveDM            .true.
"""

def parse_structure_fdf(filename):
    """Parses .fdf for lattice vectors and species info"""
    try:
        structure = structure_io.read_fdf(filename)
        return structure.lattice, structure_io.species_dict(structure)
    except FileNotFoundError:
        print(color_text(f"[ERROR] Structure file '{filename}' not found.", 'red'))
        sys.exit(1)
    except Exception as e:
        print(color_text(f"[ERROR] {e}", 'red'))
        sys.exit(1)

def generate_isolated_atom_fdf(symbol, z_num, out_path):
    """Creates a structure.fdf for a single isolated atom in a large box"""
    content = f"""# Isolated {symbol} atom for cohesive energy
NumberOfSpecies    1
NumberofAtoms      1

%block ChemicalSpeciesLabel
 1   {z_num}   {symbol}
%endblock ChemicalSpeciesLabel 

LatticeConstant 1.00 Ang 

AtomicCoordinatesFormat  Fractional 

%block LatticeVectors
 20.000000   0.000000   0.000000 
  0.000000  20.000000   0.000000 
  0.000000   0.000000  20.000000 
%endblock LatticeVectors

%block AtomicCoordinatesAndAtomicSpecies
  0.500000000   0.500000000   0.500000000   1  
%endblock AtomicCoordinatesAndAtomicSpecies
"""
    with open(out_path, 'w') as f:
        f.write(content)
    return

def link_pseudo(pp_path, symbol, target_dir):
    """Symlinks the pseudopotential if the directory is provided"""
    if not pp_path:
        return
    psml_file = f"{symbol}.psml"
    src = os.path.join(os.path.abspath(pp_path), psml_file)
    dst = os.path.join(target_dir, psml_file)
    
    if os.path.exists(src):
        try:
            os.symlink(src, dst)
        except FileExistsError:
            pass
    else:
        print(f"[WARNING] Pseudopotential '{psml_file}' not found in {pp_path}")
    return

def main():
    parser = argparse.ArgumentParser(
        description="Prepare folder structure for cohesive energy calculations.",
        epilog="Example usage:\n"
               "  stb_cohesive -s structure.fdf -k 0.15\n"
               "  stb_cohesive -s structure.fdf --spin -p /path/to/pseudos",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("-s", "--structure", dest="structure", type=str, required=True, 
                        help="Path to the initial structure.fdf file")
    
    parser.add_argument("-p", "--pp-path", dest="pp_path", type=str, default="", 
                        help="Path to the pseudopotentials folder (optional)")
    
    parser.add_argument("-k", "--k-density", dest="k_density", type=float, default=0.2, 
                        help="K-point density in 1/Angstrom (default: 0.2)")
    
    parser.add_argument("--spin", dest="spin", action="store_true", 
                        help="Set the full structure calculation to spin polarized")
    
    parser.add_argument("-v", "--version", action="version", version=f"stb-cohesive {VERSION}")

    parser.add_argument("--no-intro", dest="intro", action="store_false", 
                        help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro == True:
        show_intro([
            "Siesta ToolBox Suite",
            "Cohesive Energy Workflow setup",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("COHESIVE ENERGY:", 'bold'))
    print("-"*60)
    
    # Extract structure data
    print("\n[INFO] Read structure file ...")
    lattice, species = parse_structure_fdf(args.structure)
    print(f"[INFO] Detected species: {', '.join(species.keys())}")
    
    # Calculate K-grid for the full structure
    print("[INFO] Calculate Monkhorst-Pack grid ...")
    try:
        kgrid_divs = kspace.compute_monkhorts(lattice[0], lattice[1], lattice[2], args.k_density)
    except ValueError as e:
        print(color_text(f"[ERROR] {e}", 'red'))
        sys.exit(1)
    print(f"[INFO] Calculated K-grid for full structure: {kgrid_divs[0]} {kgrid_divs[1]} {kgrid_divs[2]} (density = {args.k_density})")

    # 1. Setup full structure directory directly in current folder
    print("\n[INFO] Setting up full structure directory ...")
    struct_dir = "structure"
    os.makedirs(struct_dir, exist_ok=True)
    
    shutil.copy(args.structure, os.path.join(struct_dir, "structure.fdf"))
    
    struct_calc = CALC_TEMPLATE
    if args.spin:
        struct_calc = struct_calc.replace("Spin                non-polarized", "Spin                polarized")
        print("[INFO] Spin polarization ENABLED for the full structure.")
    else:
        print("[INFO] Spin polarization DISABLED for the full structure.")
    
    # Apply calculated K-grid
    kgrid_str = f"kgrid.MonkhorstPack   [{kgrid_divs[0]}  {kgrid_divs[1]}  {kgrid_divs[2]}]"
    struct_calc = re.sub(r'kgrid\.MonkhorstPack\s+\[.*?\]', kgrid_str, struct_calc)
        
    with open(os.path.join(struct_dir, "calc.fdf"), 'w') as f:
        f.write(struct_calc)
        
    for sym in species.keys():
        link_pseudo(args.pp_path, sym, struct_dir)

    # 2. Setup isolated atoms directories directly in current folder
    print("\n[INFO] Setting up isolated atoms directories ...")
    atoms_root = "atoms"
    os.makedirs(atoms_root, exist_ok=True)
    
    for sym, data in species.items():
        print(f"[INFO] Setting up isolated {sym} ...")
        atom_dir = os.path.join(atoms_root, sym)
        os.makedirs(atom_dir, exist_ok=True)
        
        # Isolated atom structure
        generate_isolated_atom_fdf(sym, data['Z'], os.path.join(atom_dir, "structure.fdf"))
        
        # Calc file for isolated atom (Always polarized + Gamma point only)
        atom_calc = CALC_TEMPLATE.replace("Spin                non-polarized", "Spin                polarized")
        atom_calc = re.sub(r'kgrid\.MonkhorstPack\s+\[.*?\]', 'kgrid.MonkhorstPack   [1  1  1]', atom_calc)
        
        with open(os.path.join(atom_dir, "calc.fdf"), 'w') as f:
            f.write(atom_calc)
            
        link_pseudo(args.pp_path, sym, atom_dir)

    print("\n[INFO] Complete job!") 
    print("\n"+"-"*60)
    print(color_text("Setup complete! Folders 'structure' and 'atoms' are ready.\n\n", 'bold'))

if __name__ == "__main__":
    main()
