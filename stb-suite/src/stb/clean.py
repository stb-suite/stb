#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################
    
VERSION = "1.9.1"  

import os
import sys
import warnings
import subprocess
from time import sleep
import argparse
import textwrap
from typing import List, Dict
import numpy as np
import re
import argparse
from stb.core.cli import COLORS, color_text, show_intro

def should_delete(file, allowed_exts):
    return os.path.isfile(file) and os.path.splitext(file)[1] not in allowed_exts

def main():
    parser = argparse.ArgumentParser(description="Remove all files except those with specified extensions.")
    
    parser.add_argument(
        '--keep', nargs='+', default=['.psml', '.psf', '.fdf', '.sh'],
        help="List of extensions to keep (e.g. .fdf .sh)"
    )
    parser.add_argument(
        '--dry-run', action='store_true',
        help="Only print the files that would be removed"
    )
    parser.add_argument(
        '--no-confirm', '-n', action='store_true',
        help="Do not ask for confirmation before deleting files"
    )
    parser.add_argument(
        '--path', default='.',
        help="Directory to clean (default: current directory)"
    )
    parser.add_argument("-v", "--version", action="version",
                        version=f"stb-clean {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    allowed_exts = set(args.keep)
    
    if args.intro == True:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])
    print("\n" + color_text("Clean:", 'bold'))
    print("-"*60)
    print("Remove all files except those with specified extensions.\n")

    for file in os.listdir(args.path):
        full_path = os.path.join(args.path, file)
        if should_delete(full_path, allowed_exts):
            if args.dry_run:
                print(f"[Dry-run] Would remove: {file}")
            elif args.no_confirm:
                os.remove(full_path)
                print(f"[INFO] Removed: {file}")
            else:
                answer = input(f"[CONFIRM] Delete {file}? [y/N] ").strip().lower()
                if answer == 'y':
                    os.remove(full_path)
                    print(f"[INFO] Removed: {file}")
                else:
                    print(f"[INFO] Skipped: {file}")
    
    print("\n[INFO] Complete job!") 
    print("\n"+"-"*60)
    print(color_text("Cleaned up! May the deleted files rest in pieces.\n\n", 'bold'))

if __name__ == "__main__":
    main()

