#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Elastic Generator       #
#################################################

VERSION = "2.1.2"

import os
import sys
import argparse
import numpy as np
from time import sleep
from stb.core import structure_io
from stb.core.cli import COLORS, color_text, show_intro

# ==========================================
#           HELPER FUNCTIONS
# ==========================================

def get_strain_matrix(direction, delta):
    """Returns the deformation matrix (I + epsilon)."""
    epsilon = np.zeros((3, 3))
    d = direction.lower()
    
    # Robust mapping (accepts x or xx, zx or xz)
    if d in ['x', 'xx']:     
        epsilon[0, 0] = delta
    elif d in ['y', 'yy']:   
        epsilon[1, 1] = delta
    elif d in ['z', 'zz']:   
        epsilon[2, 2] = delta
    
    elif d in ['yz', 'zy']:  
        epsilon[1, 2] = epsilon[2, 1] = delta
    elif d in ['xz', 'zx']:  
        epsilon[0, 2] = epsilon[2, 0] = delta
    elif d in ['xy', 'yx']:  
        epsilon[0, 1] = epsilon[1, 0] = delta
    
    elif d == 'bi':    
        epsilon[0, 0] = epsilon[1, 1] = delta
    elif d == 'hydro': 
        epsilon[0, 0] = epsilon[1, 1] = epsilon[2, 2] = delta
            
    return np.eye(3) + epsilon

def generate_verify_script():
    """Generates a small helper script to check calculation status."""
    content = r"""#!/bin/bash
for d in strain_*/ ; do
    [ -d "$d" ] && echo "Checking $d..." && tail -n 1 "$d/calc.out" 2>/dev/null
done
"""
    with open("verify_calc.sh", "w") as f:
        f.write(content)
    try:
        os.chmod("verify_calc.sh", 0o755)
    except OSError:
        pass # Ignore if permissions cannot be changed (e.g., Windows)

# ==========================================
#                  MAIN
# ==========================================

def main():
    desc_text = "Siesta Elastic Strain Generator (Batch Mode)"
    # Detailed help description
    help_dirs = (
        "List of directions to generate. \n"
        "Options: [xx, yy, zz, xy, xz, yz, bi, hydro, all]. \n"
        "Example: --dirs xx yy xy"
    )

    parser = argparse.ArgumentParser(description=desc_text, formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("--file", "-i", required=True, help="Input structural FDF file")
    parser.add_argument("--dirs", nargs='+', default=["all"], help=help_dirs)
    parser.add_argument("--max", type=float, default=2.0, help="Max strain %%")
    parser.add_argument("--steps", type=int, default=4, help="Steps per direction")
    parser.add_argument("--output", default="structure.fdf", help="Output filename")
    
    parser.add_argument("-v", "--version", action="version",
                        version=f"stb-elasticInputs {VERSION}")
    # --- NEW: Argument --no-intro ---
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    
    args = parser.parse_args()

    # Shows intro if --no-intro flag is NOT used
    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - Elastic Generator",
            "Strain Structure Generator (Batch Mode)",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Integrated Style Refactoring"
        ])

    # Expand 'all' keyword
    dirs_to_run = []
    if 'all' in args.dirs:
        dirs_to_run = ['xx', 'yy', 'zz', 'xy', 'xz', 'yz']
    else:
        dirs_to_run = args.dirs

    try:
        lattice = structure_io.raw_lattice_vectors(args.file)
        print(f"{color_text('[INFO]', 'green')} Loaded: {args.file}")
    except Exception as e:
        sys.exit(f"{color_text('[ERROR]', 'red')} {e}")

    print(f"{color_text('[INFO]', 'green')} Modes: {dirs_to_run}")
    
    count = 0
    strains = np.linspace(-args.max, args.max, args.steps)
    
    for d in dirs_to_run:
        for s in strains:
            if abs(s) < 1e-9: s = 0.0
            
            delta = s / 100.0
            def_matrix = get_strain_matrix(d, delta)
            new_lattice = np.dot(lattice, def_matrix)
            
            # Create folder: e.g. strain_xx_-1.0
            folder = f"strain_{d}_{'m' if s < 0 else ''}{abs(s):.2f}"
            if not os.path.exists(folder): os.makedirs(folder)
            
            structure_io.rewrite_fdf_lattice(args.file, new_lattice, os.path.join(folder, args.output))
            count += 1
            
    print(f"\n{color_text('[SUCCESS]', 'green')} Generated {count} structures.")
    generate_verify_script()
    print(f"Run calculations inside each 'strain_...' folder.")

if __name__ == "__main__":
    main()
