#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Elastic Properties      #
#     Analyzer Module (2D/3D Support)           #
#     Integrated Style Refactoring              #
#################################################

VERSION = "1.9.1" # Atualizado versão para refletir a mudança

import os
import sys
import argparse
import re
import numpy as np
from time import sleep
from scipy.stats import linregress
from stb.core import siesta_log
from stb.core.cli import COLORS, color_text, show_intro

# ==========================================
#           HELPERS & CONFIG
# ==========================================

CONV_EVA3_TO_GPA = 160.21766
CONV_EVA2_TO_NM = 16.021766
REPORT_FILE = "mechanical_properties.txt"

def print_dual(text, file_handle=None):
    """Prints to stdout with color, writes to file without color."""
    print(text)
    if file_handle:
        # Regex to strip ANSI escape codes
        clean_text = re.sub(r'\x1b\[[0-9;]*m', '', text)
        file_handle.write(clean_text + "\n")

# ==========================================
#           DATA MINING LOGIC
# ==========================================

def calculate_slope(x, y, factor):
    if len(x) < 2: return 0.0
    slope, _, _, _, _ = linregress(x, y)
    return slope * factor

# ==========================================
#           REPORTING LOGIC
# ==========================================

def check_stability_and_report(C, is_2d, unit_label, f_out):
    print_dual(f"\n{color_text('[4] STABILITY AND PROPERTIES', 'magenta')}", f_out)
    print_dual("-" * 60, f_out)
    
    passes = False
    
    if is_2d:
        # 2D Constants Extraction
        c11, c22, c12, c66 = C[0,0], C[1,1], C[0,1], C[5,5]
        
        # In-plane Properties (2D Orthorhombic Formulas)
        denom = (c11 * c22 - c12**2)
        Ex_2d = denom / c22 if c22 > 0 else 0
        Ey_2d = denom / c11 if c11 > 0 else 0
        
        # Poisson's Ratio
        v12 = c12 / c22 if c22 > 0 else 0 # Strain in x due to stress in y
        v21 = c12 / c11 if c11 > 0 else 0 # Strain in y due to stress in x
        
        print_dual(f"In-Plane Stiffness (Layer Modulus - Ex): {color_text(f'{Ex_2d:.2f}', 'bold')} {unit_label}", f_out)
        print_dual(f"In-Plane Stiffness (Layer Modulus - Ey): {color_text(f'{Ey_2d:.2f}', 'bold')} {unit_label}", f_out)
        print_dual(f"Poisson's Ratio (v_21)                 : {color_text(f'{v21:.3f}', 'bold')}", f_out)
        print_dual(f"Shear Modulus (G_xy = C66)             : {color_text(f'{c66:.2f}', 'bold')} {unit_label}", f_out)
        print_dual("", f_out)

        # 2D Born Criteria
        cond1 = c11 > 0
        cond2 = c66 > 0
        cond3 = denom > 0
        passes = cond1 and cond2 and cond3
        
        print_dual(f"2D Stability Criteria (Born):", f_out)
        print_dual(f"  C11 > 0         : {color_text('PASS', 'green') if cond1 else color_text('FAIL', 'red')}", f_out)
        print_dual(f"  C66 > 0         : {color_text('PASS', 'green') if cond2 else color_text('FAIL', 'red')}", f_out)
        print_dual(f"  C11*C22 > C12^2 : {color_text('PASS', 'green') if cond3 else color_text('FAIL', 'red')}", f_out)

    else:
        # 3D Case (Hill Average)
        # Voigt
        B_v = (C[0,0]+C[1,1]+C[2,2] + 2*(C[0,1]+C[0,2]+C[1,2])) / 9.0
        G_v = ((C[0,0]+C[1,1]+C[2,2]) - (C[0,1]+C[0,2]+C[1,2]) + 3*(C[3,3]+C[4,4]+C[5,5])) / 15.0
        # Reuss
        try:
            S_mat = np.linalg.inv(C)
            B_r = 1.0 / ((S_mat[0,0]+S_mat[1,1]+S_mat[2,2]) + 2*(S_mat[0,1]+S_mat[0,2]+S_mat[1,2]))
            G_r = 15.0 / (4*(S_mat[0,0]+S_mat[1,1]+S_mat[2,2]) - 4*(S_mat[0,1]+S_mat[0,2]+S_mat[1,2]) + 3*(S_mat[3,3]+S_mat[4,4]+S_mat[5,5]))
        except:
            B_r, G_r = B_v, G_v
        
        # Hill
        B_h = (B_v + B_r) / 2.0
        G_h = (G_v + G_r) / 2.0
        E_h = (9 * B_h * G_h) / (3 * B_h + G_h)
        nu_h = (3 * B_h - 2 * G_h) / (2 * (3 * B_h + G_h))

        print_dual(f"Young's Modulus (E)  : {color_text(f'{E_h:.2f}', 'bold')} {unit_label}", f_out)
        print_dual(f"Bulk Modulus    (B)  : {color_text(f'{B_h:.2f}', 'bold')} {unit_label}", f_out)
        print_dual(f"Poisson's Ratio (v)  : {color_text(f'{nu_h:.3f}', 'bold')}", f_out)
        print_dual("", f_out)

        eigenvals = np.linalg.eigvals(C)
        passes = np.all(eigenvals > 0.01)

    if passes:
        print_dual(f"\n{color_text('VERDICT:', 'bold')} {color_text('STABLE', 'green')}", f_out)
    else:
        print_dual(f"\n{color_text('VERDICT:', 'bold')} {color_text('UNSTABLE', 'red')}", f_out)

# ==========================================
#                  MAIN
# ==========================================

def main():
    parser = argparse.ArgumentParser(description="Siesta Elastic Properties Analyzer (STB Suite)")
    
    # --- CHANGED: Explicit flag for filename ---
    parser.add_argument("-f", "--file", type=str, default="calc.out", 
                        help="Siesta log/output filename inside strain folders (default: calc.out)")
    
    parser.add_argument("--2d", dest="is2d", action="store_true", help="Enable 2D mode (N/m)")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - Elastic Analyzer",
            "Elastic Constants & Stability Calculator",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Integrated Style Refactoring"
        ])

    print("\n" + color_text("ELASTIC ANALYZER:", 'bold'))
    print("-" * 60)

    # --- Setup ---
    print(f"[INFO] Target Log File: {color_text(args.file, 'yellow')}")
    
    if args.is2d:
        unit_label = "N/m"
        print(f"[INFO] Mode: {color_text('2D Material', 'cyan')} (Units: N/m)")
    else:
        unit_label = "GPa"
        print(f"[INFO] Mode: {color_text('3D Bulk', 'blue')} (Units: GPa)")

    # --- Data Mining ---
    print(f"[INFO] Scanning for 'strain_*' folders...")
    folders = sorted([f for f in os.listdir('.') if os.path.isdir(f) and f.startswith("strain_")])
    data = {}

    Lz = 1.0
    found_lz = False

    # Pre-scan for Lz (Unit Cell Height) if 2D
    for folder in folders:
        path = os.path.join(folder, args.file)
        if os.path.exists(path):
            val = siesta_log.get_cell_height(path)
            if val > 1.0:
                Lz = val
                found_lz = True
                break
    
    if args.is2d and found_lz:
        print(f"[INFO] Normalization by Cell Height (Lz): {Lz:.2f} Ang")
    
    CONV_FACTOR = Lz * CONV_EVA2_TO_NM if args.is2d else CONV_EVA3_TO_GPA

    with open(REPORT_FILE, "w") as f_out:
        print_dual(f"{color_text('===== ELASTIC PROPERTIES REPORT =====', 'magenta')}", f_out)
    
        loaded_count = 0
        print(f"\n{color_text('READING FOLDERS:', 'bold')}")
    
        for folder in folders:
            direction, val = siesta_log.parse_strain_folder_name(folder)
            if direction is None: continue
            strain = val / 100.0

            # Read the specified file
            S = siesta_log.get_stress_tensor(os.path.join(folder, args.file))
        
            # --- MODIFICAÇÃO: Imprimir status da pasta ---
            msg = f"   -> {folder:<25} : "
            if S is not None:
                if direction not in data: data[direction] = {'eps': [], 'stress': []}
                data[direction]['eps'].append(strain)
                data[direction]['stress'].append(S)
                loaded_count += 1
                msg += color_text("OK", 'green')
            else:
                msg += color_text("FAIL", 'red') + " (No stress data found)"
        
            print(msg)
            # ---------------------------------------------
            
        print(f"\n[OK] Loaded {loaded_count} calculations from {len(folders)} folders.")

        if not data:
            print(f"{color_text('[FAIL]', 'red')} No valid data found in strain folders.")
            print(f"       Check if '{args.file}' exists and has 'Stress tensor'.")
            sys.exit(1)


        # --- C_ij Calculation ---
        C = np.zeros((6, 6))

        # Eixo X (C11, C21, C31) - Tenta 'xx' ou 'x'
        if 'xx' in data:
            e, S = data['xx']['eps'], np.array(data['xx']['stress'])
            for i in range(3): C[i, 0] = calculate_slope(e, S[:, i, i], CONV_FACTOR)
        elif 'x' in data:
            e, S = data['x']['eps'], np.array(data['x']['stress'])
            for i in range(3): C[i, 0] = calculate_slope(e, S[:, i, i], CONV_FACTOR)

        # Eixo Y (C12, C22, C32) - Tenta 'yy' ou 'y'
        if 'yy' in data:
            e, S = data['yy']['eps'], np.array(data['yy']['stress'])
            for i in range(3): C[i, 1] = calculate_slope(e, S[:, i, i], CONV_FACTOR)
        elif 'y' in data:
            e, S = data['y']['eps'], np.array(data['y']['stress'])
            for i in range(3): C[i, 1] = calculate_slope(e, S[:, i, i], CONV_FACTOR)

        # Eixo Z (C13, C23, C33) - Tenta 'zz' ou 'z'
        if 'zz' in data:
            e, S = data['zz']['eps'], np.array(data['zz']['stress'])
            for i in range(3): C[i, 2] = calculate_slope(e, S[:, i, i], CONV_FACTOR)
        elif 'z' in data:
            e, S = data['z']['eps'], np.array(data['z']['stress'])
            for i in range(3): C[i, 2] = calculate_slope(e, S[:, i, i], CONV_FACTOR)

        # C44 (Shear YZ)
        if 'yz' in data: 
            C[3, 3] = calculate_slope(data['yz']['eps'], np.array(data['yz']['stress'])[:, 1, 2], CONV_FACTOR) / 2.0

        # C55 (Shear XZ ou ZX) - Suas pastas usam 'zx'
        if 'zx' in data: 
            C[4, 4] = calculate_slope(data['zx']['eps'], np.array(data['zx']['stress'])[:, 0, 2], CONV_FACTOR) / 2.0
        elif 'xz' in data:
            C[4, 4] = calculate_slope(data['xz']['eps'], np.array(data['xz']['stress'])[:, 0, 2], CONV_FACTOR) / 2.0

        # C66 (Shear XY)
        if 'xy' in data: 
            C[5, 5] = calculate_slope(data['xy']['eps'], np.array(data['xy']['stress'])[:, 0, 1], CONV_FACTOR) / 2.0

        C_sym = 0.5 * (C + C.T)

        # --- CONDITIONAL OUTPUT ---
        if args.is2d:
            # 2D REPORT
            print_dual(f"\n{color_text('[1] 2D STIFFNESS MATRIX (' + unit_label + ')', 'magenta')}", f_out)
            print_dual("      xx       yy       xy", f_out)
        
            indices_2d = [0, 1, 5]
            labels_2d = ["xx", "yy", "xy"]
        
            for r, label_r in zip(indices_2d, labels_2d):
                row_str = f"{label_r} | "
                for c in indices_2d:
                    val = C_sym[r, c]
                    if abs(val) < 0.01: val = 0.0
                    row_str += f"{val:8.2f} "
                print_dual(row_str, f_out)
            
            print_dual(f"\n{color_text('[2] RELEVANT CONSTANTS (2D)', 'magenta')}", f_out)
            labels_show = ["C11", "C22", "C12", "C66"]
            idx_show = [(0,0), (1,1), (0,1), (5,5)]
        
            buf = ""
            for l, (r, c) in zip(labels_show, idx_show):
                buf += f"{color_text(l, 'cyan')}: {color_text(f'{C_sym[r,c]:6.2f}', 'bold')} {unit_label}   "
            print_dual(buf, f_out)

        else:
            # 3D REPORT
            if abs(C_sym[2,2]) < 0.1 * abs(C_sym[0,0]) and abs(C_sym[0,0]) > 10.0:
                print_dual(f"\n{color_text('[WARNING] Material seems to be 2D (C33 << C11). Consider using --2d', 'yellow')}", f_out)
            
            print_dual(f"\n{color_text('[1] 3D STIFFNESS MATRIX (' + unit_label + ')', 'magenta')}", f_out)
            header = "      " + "".join([f"{f'j={j+1}':<9}" for j in range(6)])
            print_dual(f"{header}", f_out)
            for i in range(6):
                print_dual(f"i={i+1} | " + "".join([f"{C_sym[i,j]:9.2f}" for j in range(6)]), f_out)
            
            print_dual(f"\n{color_text('[2] MAIN CONSTANTS', 'magenta')}", f_out)
            labels = ["C11", "C22", "C33", "C44", "C55", "C66", "C12", "C13", "C23"]
            idx = [(0,0), (1,1), (2,2), (3,3), (4,4), (5,5), (0,1), (0,2), (1,2)]
            buf = ""
            for k, (l, (r, c)) in enumerate(zip(labels, idx)):
                buf += f"{l}: {C_sym[r,c]:7.1f} {unit_label}  "
                if (k+1)%3==0: 
                    print_dual(buf, f_out); buf=""
            if buf: print_dual(buf, f_out)

        # Stability and Properties Check
        check_stability_and_report(C_sym, args.is2d, unit_label, f_out)
    
    print("-" * 60)
    print(f"[DONE] Report saved to: {color_text(REPORT_FILE, 'green')}")
    print(color_text("\nScience is organized knowledge. Wisdom is organized life.", 'bold'))

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print(f"\n{color_text('[INFO] Operation cancelled by user.', 'red')}")
        sys.exit(0)
