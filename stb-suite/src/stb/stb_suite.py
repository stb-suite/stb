#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################
    
VERSION = "1.9.1"  

import os
import subprocess
from time import sleep
import argparse
import textwrap
from typing import List, Dict, Callable

try:
    import readline
    readline.parse_and_bind("tab: complete")
except ImportError:
    pass

from stb.core.cli import COLORS, color_text, show_intro, get_input, get_float_input, get_int_input
from stb.core.pseudopotentials import BANKS

def prompt_pseudo_source(optional: bool = True) -> str:
    """Shared pseudopotential-source prompt for every wrapper below that
    needs one (phonons, cohesive energy, input file, Hubbard U prep):
    a bundled bank (see core/pseudopotentials.py) or a custom path. Returns
    the raw string to pass straight through as the tool's -p/--pp-path/
    --pseudo-dir value (each tool resolves it itself); empty string only if
    `optional` and the user skips."""
    bank_list = list(BANKS.items())
    print(f"\n{color_text('Pseudopotential source:', 'yellow')}")
    for i, (name, description) in enumerate(bank_list, 1):
        print(f"  {color_text(str(i), 'cyan')} = Bundled: {description} ({name})")
    print(f"  {color_text(str(len(bank_list) + 1), 'cyan')} = Custom path")
    prompt = f"Select option (1-{len(bank_list) + 1}"
    prompt += ", or Enter to skip): " if optional else "): "
    while True:
        choice = get_input(prompt).strip()
        if not choice and optional:
            return ""
        if choice.isdigit() and 1 <= int(choice) <= len(bank_list):
            return bank_list[int(choice) - 1][0]
        if choice == str(len(bank_list) + 1):
            path = os.path.expanduser(get_input("Custom pseudopotentials folder path: ").strip())
            if os.path.isdir(path):
                return path
            print(color_text(f"Path not found: '{path}'", 'red'))
            continue
        print(color_text("Invalid choice.", 'red'))


def show_main_menu() -> None:
    """Displays the main category menu"""
    print("\n" + color_text("STB-SUITE Main Menu:", 'bold'))
    print("-"*60)
    print(f"{color_text('1.', 'yellow')} {color_text('Inputs', 'blue')}\n    Tools to set up a SIESTA run (input file, k-grid, k-path)\n")
    print(f"{color_text('2.', 'yellow')} {color_text('Structures', 'blue')}\n    Tools to build, generate, or transform structure files (stacking, supercells, slabs, defects, SQS, etc.)\n")
    print(f"{color_text('3.', 'yellow')} {color_text('Analysis', 'blue')}\n    Tools to analyze simulation results (bands, DOS, structures, charge density, etc.)\n")
    print(f"{color_text('4.', 'yellow')} {color_text('Workflow', 'blue')}\n    Complete prep + analysis pipelines for a specific property (strain, elastic constants, cohesive energy, phonons)\n")
    print(f"{color_text('5.', 'yellow')} {color_text('Utils', 'blue')}\n    Helper tools for file management and conversion\n")
    print(f"{color_text('0.', 'yellow')} {color_text('Exit', 'red')}")
    print("-"*60)
    print(color_text("Tip: you can also type a tool code directly (e.g. 4.1.2) to jump straight to it.", 'yellow'))

def show_sub_menu(title: str, tools_dict: Dict) -> None:
    """Displays a sub-menu for a specific tool category"""
    print("\n" + "="*60)
    print(color_text(f"--- {title} ---", 'cyan').center(68))
    print("="*60 + "\n")
    
    for key, info in tools_dict.items():
        menu_title = color_text(info['title'], 'blue')
        desc = textwrap.fill(info['description'], width=55, subsequent_indent='    ')
        print(f"{color_text(str(key)+'.', 'yellow')} {menu_title}\n    {desc}\n")
    
    print(f"{color_text('0.', 'yellow')} {color_text('Back to Main Menu', 'red')}")
    print("-"*60)

def run_tool(tool_name: str, args: List[str]) -> None:
    """Executes a suite tool as a subprocess"""
    try:
        cmd = [f"{tool_name}"] + args
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(color_text(f"\nError running {tool_name}: {e}", 'red'))
    except FileNotFoundError:
        print(color_text(f"\nTool {tool_name} not found!", 'red'))
        print(color_text(f"Make sure {tool_name} is in your system's PATH.", 'yellow'))
    input("\nPress Enter to continue...")

# ==========================================================
# TOOL FUNCTIONS
# ==========================================================


def run_phonon_postprocessing() -> None:
    """Interface for the Phonon Post-Processing (phonons_post.py)"""
    print("\n" + "="*60)
    print(color_text("PHONON POST-PROCESSING", 'bold').center(60))
    print("="*60 + "\n")
    
    # 1. Diretório
    phonon_dir = get_input("Phonon runs directory [default: phonon_runs]: ").strip()
    if not phonon_dir:
        phonon_dir = "phonon_runs"
        
    # 2. System Label
    sys_label = get_input("SystemLabel used in calculations [default: siesta]: ").strip()
    if not sys_label:
        sys_label = "siesta"
        
    # 3. Malha Q (Mesh)
    mesh_input = get_input("\nQ-point mesh (e.g. '20 20 20') [default: 20 20 20]: ").strip()
    if not mesh_input:
        m_x, m_y, m_z = 20, 20, 20
    else:
        try:
            dims = [int(x) for x in mesh_input.split()]
            if len(dims) == 3:
                m_x, m_y, m_z = dims
            else:
                print(color_text("Please provide exactly 3 integers. Using default 20 20 20.", 'yellow'))
                m_x, m_y, m_z = 20, 20, 20
        except ValueError:
            print(color_text("Invalid input format. Using default 20 20 20.", 'yellow'))
            m_x, m_y, m_z = 20, 20, 20

    # 4. Temperaturas
    print(f"\n{color_text('Thermal Properties Settings:', 'yellow')}")
    tmin = get_float_input("Minimum temperature (K) [default: 0]: ", 0.0)
    tmax = get_float_input("Maximum temperature (K) [default: 1000]: ", 1000.0)
    tstep = get_float_input("Temperature step (K) [default: 10]: ", 10.0)
        

    args = [
        "-dir", phonon_dir,
        "-l", sys_label,
        "-m", str(m_x), str(m_y), str(m_z),
        "--tmin", str(tmin),
        "--tmax", str(tmax),
        "--tstep", str(tstep),
        "--no-intro"
    ]

    print(color_text("\nStarting Phonon post-processing...", 'green'))
    run_tool("stb-phononsPos", args)
    
    


def run_phonon_generator() -> None:
    """Interface for the Phonon Displacement Generator (phonons_create.py)"""
    print("\n" + "="*60)
    print(color_text("PHONON DISPLACEMENT GENERATOR", 'bold').center(60))
    print("="*60 + "\n")
    
    # 1. Obter arquivo de estrutura
    structure_file = get_input("Input structure file [default: structure.fdf]: ").strip()
    if not structure_file:
        structure_file = "structure.fdf"
        
    # 2. Obter arquivo de cálculo
    calc_file = get_input("Calculation parameters file [default: calc.fdf]: ").strip()
    if not calc_file:
        calc_file = "calc.fdf"
        
    # 3. Definir dimensões da supercélula (entrada única separada por espaços)
    dim_input = get_input("\nSupercell dimensions (e.g. '2 2 2') [default: 2 2 2]: ").strip()
    
    if not dim_input:
        dim_x, dim_y, dim_z = 2, 2, 2
    else:
        try:
            dims = [int(x) for x in dim_input.split()]
            if len(dims) == 3:
                dim_x, dim_y, dim_z = dims
            else:
                print(color_text("Please provide exactly 3 integers. Using default 2 2 2.", 'yellow'))
                dim_x, dim_y, dim_z = 2, 2, 2
        except ValueError:
            print(color_text("Invalid input format. Using default 2 2 2.", 'yellow'))
            dim_x, dim_y, dim_z = 2, 2, 2
    
    # 4. Definir distância de deslocamento
    distance = get_float_input("\nDisplacement distance in Å [default: 0.01]: ", 0.01)
    
    # 5. Diretório dos pseudopotenciais
    pseudo_dir = prompt_pseudo_source(optional=True)
    if not pseudo_dir:
        pseudo_dir = "."
        
    # 6. Preparar e executar o script
    args = [
        "-s", structure_file,
        "-c", calc_file,
        "-dim", str(dim_x), str(dim_y), str(dim_z),
        "-d", str(distance),
        "-p", pseudo_dir,
        "--no-intro"
    ]

    print(color_text("\nGenerating phonon displacement folders...", 'green'))
    run_tool("stb-phononsCreate", args)


def run_cohesive_setup() -> None:
    """Interface for the Cohesive Energy Setup (cohesive_energy.py)"""
    print("\n" + "="*60)
    print(color_text("COHESIVE ENERGY SETUP", 'bold').center(60))
    print("="*60 + "\n")
    
    # 1. Obter arquivo de estrutura
    struct_file = get_input("Input structure FDF file (-s): ").strip()
    while not os.path.isfile(struct_file):
        print(color_text("File not found!", 'red'))
        struct_file = get_input("Input structure FDF file (-s): ").strip()
        
    # 2. Obter densidade K
    k_density = get_float_input("K-point density (default: 0.2): ", 0.2)
    while k_density <= 0:
        print(color_text("Density must be a positive number!", 'red'))
        k_density = get_float_input("K-point density (default: 0.2): ", 0.2)

    # 3. Obter caminho do PP
    pp_path = prompt_pseudo_source(optional=True)
    
    # 4. Spin polarization
    spin_choice = get_input("Enable spin polarization for full structure? (y/N): ").strip().lower()

    # 5. Isolated-atom vacuum box size
    vacuum = get_float_input("Isolated-atom vacuum box side in Ang (default: 20.0): ", 20.0)
    while vacuum <= 0:
        print(color_text("Vacuum box side must be a positive number!", 'red'))
        vacuum = get_float_input("Isolated-atom vacuum box side in Ang (default: 20.0): ", 20.0)

    args = [
        "-s", struct_file,
        "-k", str(k_density),
        "--vacuum", str(vacuum),
        "--no-intro"
    ]

    if pp_path:
        args.extend(["-p", pp_path])
    if spin_choice in ['y', 'yes']:
        args.append("--spin")

    # Executa o script. Se não tiver os atalhos globais configurados,
    # pode alterar "stb_cohesive" para "python cohesive_energy.py" 
    run_tool("stb-cohesive", args)


def run_cohesive_analysis() -> None:
    """Interface for the Cohesive Energy Analysis (cohesive_analysis.py)"""
    print("\n" + "="*60)
    print(color_text("COHESIVE ENERGY ANALYSIS", 'bold').center(60))
    print("="*60 + "\n")
    
    # 1. Obter nome do ficheiro de output
    out_file = get_input("SIESTA output file name (e.g., calc.out) [-o]: ").strip()
    while not out_file:
        print(color_text("File name cannot be empty!", 'red'))
        out_file = get_input("SIESTA output file name [-o]: ").strip()
        
    # 2. Obter o diretório alvo
    dir_path = get_input("Path to results folder containing 'structure' and 'atoms' (default: current dir) [-d]: ").strip()
    
    args = [
        "-o", out_file,
        "--no-intro"
    ]
    
    if dir_path:
        args.extend(["-d", dir_path])
        
    # Executa o script. Se não tiver atalho configurado, altere para "python cohesive_analysis.py"
    run_tool("stb-cohesiveAnalysis", args)


def run_2d_stacker() -> None:
    """Interface for the Monolayer Stacker (stb.stacking2D:main)"""
    print("\n" + "="*60)
    print(color_text("2D MONOLAYER STACKER", 'bold').center(60))
    print("="*60 + "\n")
    
    # 1. Get Input Files
    layer1 = get_input("Bottom Monolayer FDF file (-l1): ").strip()
    while not os.path.isfile(layer1):
        print(color_text("File not found!", 'red'))
        layer1 = get_input("Bottom Monolayer FDF file (-l1): ").strip()

    layer2 = get_input("Top Monolayer FDF file (-l2): ").strip()
    while not os.path.isfile(layer2):
        print(color_text("File not found!", 'red'))
        layer2 = get_input("Top Monolayer FDF file (-l2): ").strip()

    # 2. Basic numerical parameters
    max_area = get_float_input("Max supercell area in Å² (default: 150.0): ", 150.0)
    max_strain = get_float_input("Max allowed strain fraction (default: 0.05): ", 0.05)

    # 3. Van der Waals Gap Option (Numbered Menu)
    print(f"\n{color_text('Select Van der Waals Gap Option:', 'yellow')}")
    print(f"  {color_text('1', 'cyan')} = Default (3.2 Å)")
    print(f"  {color_text('2', 'cyan')} = Manual (Single value)")
    print(f"  {color_text('3', 'cyan')} = Range (Multiple values for Energy Curve)")
    gap_choice = get_input("Select option (1-3) [default: 1]: ").strip()
    
    gap_args = []
    if gap_choice == '2':
        val = get_float_input("Enter Gap in Å: ", 3.2)
        gap_args = ["-g", str(val)]
    elif gap_choice == '3':
        g_start = get_float_input("Start Gap in Å: ", 3.0)
        g_end = get_float_input("End Gap in Å: ", 4.0)
        # O script nativo usa np.linspace, logo precisa do número total de pontos (steps)
        g_pts = int(get_float_input("Number of points/steps (e.g., 11): ", 11))
        gap_args = ["--gap_range", str(g_start), str(g_end), str(g_pts)]
    else:
        gap_args = ["-g", "3.2"]

    # 4. Stacking & Symmetry Mode (Numbered Menu)
    print(f"\n{color_text('Select Stacking Mode:', 'yellow')}")
    print(f"  {color_text('1', 'cyan')} = Default (twist=0.0, tx=0.0, ty=0.0)")
    print(f"  {color_text('2', 'cyan')} = Manual (Define twist, tx, ty)")
    print(f"  {color_text('3', 'cyan')} = High-Symmetry Points (Batch Mode)")
    stack_mode = get_input("Select option (1-3) [default: 1]: ").strip()
    
    batch_sym = False
    twist = 0.0
    shift_x = 0.0
    shift_y = 0.0
    
    if stack_mode == '3':
        batch_sym = True
    elif stack_mode == '2':
        twist = get_float_input("Initial twist angle in degrees (default: 0.0): ", 0.0)
        shift_x = get_float_input("Fractional shift for layer 2 in X axis [-tx] (default: 0.0): ", 0.0)
        shift_y = get_float_input("Fractional shift for layer 2 in Y axis [-ty] (default: 0.0): ", 0.0)

    # 5. Strain Distribution Mode (Numbered Menu)
    print(f"\n{color_text('Select Strain Distribution Mode:', 'yellow')}")
    print(f"  {color_text('1', 'cyan')} = Top (Strain layer 2 to match layer 1) [Default]")
    print(f"  {color_text('2', 'cyan')} = Bottom (Strain layer 1 to match layer 2)")
    print(f"  {color_text('3', 'cyan')} = Sym (Symmetric strain on both layers)")
    sm_choice = get_input("Select mode (1-3) [default: 1]: ").strip()
    
    sm_map = {'1': 'top', '2': 'bottom', '3': 'sym'}
    strain_mode = sm_map.get(sm_choice, 'top')

    # 6. Build Command and Execute ONCE
    args = [
        "-l1", layer1,
        "-l2", layer2,
        "-i", # Always included as requested
        "-a", str(max_area),
        "-s", str(max_strain),
        "-sm", strain_mode,
        "--no-intro"
    ]

    # Adiciona os argumentos de Gap (seja -g ou --gap_range)
    args.extend(gap_args)

    # Adiciona os argumentos de Simetria e Empilhamento
    if batch_sym:
        args.append("--batch_sym")
    else:
        args.extend(["-t", str(twist), "-tx", str(shift_x), "-ty", str(shift_y)])

    print(color_text(f"\n--- Running 2D Stacker ---", 'green'))
    run_tool("stb-2Dstacking", args)


def run_grid_to_cube() -> None:
    """Interface for the Grid to Cube Converter (cube.py)"""
    print("\n" + "="*60)
    print(color_text("SIESTA GRID TO CUBE CONVERTER", 'bold').center(60))
    print("="*60 + "\n")
    
    # 1. Obter o SystemLabel
    label = get_input("Enter the Siesta SystemLabel (e.g., siesta): ").strip()
    while not label:
        print(color_text("Label cannot be empty!", 'red'))
        label = get_input("Enter the Siesta SystemLabel: ")

    # 2. Escolher o Tipo de Ficheiro
    print(f"\n{color_text('Select Grid Type to Convert:', 'yellow')}")
    print(f"  {color_text('1', 'cyan')} = RHO (Charge Density)")
    print(f"  {color_text('2', 'cyan')} = VT  (Total Potential)")
    print(f"  {color_text('3', 'cyan')} = VH  (Hartree Potential)")
    print(f"  {color_text('4', 'cyan')} = BADER (Charge Analysis Grid)")
    
    type_map = {'1': 'RHO', '2': 'VT', '3': 'VH', '4': 'BADER'}
    choice = get_input("Select type (1-4) [default: 1]: ").strip()
    
    selected_type = type_map.get(choice, 'RHO') # Default is RHO
    
    print(f"\nTarget File: {color_text(f'{label}.{selected_type}', 'cyan')}")

    # 3. Execução
    args = ["--label", label, "--type", selected_type, "--no-intro"]

    print(color_text("\nConverting to Cube format...", 'green'))
    run_tool("stb-cube", args)



def run_density_plotter() -> None:
    """Interface for the Charge Density Plotter (density.py)"""
    print("\n" + "="*60)
    print(color_text("CHARGE DENSITY PLOTTER (RHO)", 'bold').center(60))
    print("="*60 + "\n")
    
    # 1. Obter o SystemLabel
    print(f"This tool reads the {color_text('.RHO', 'yellow')} file generated by Siesta.")
    label = get_input("Enter the Siesta SystemLabel (e.g., siesta): ").strip()
    while not label:
        print(color_text("Label cannot be empty!", 'red'))
        label = get_input("Enter the Siesta SystemLabel: ")

    # 2. Escolher Modo (2D, 3D ou Perfil)
    print(f"\n{color_text('Plot Mode:', 'yellow')}")
    print(f"  {color_text('1', 'cyan')} = 2D Slice (Planar Cut)")
    print(f"  {color_text('2', 'cyan')} = 3D Volume (Point Cloud)")
    print(f"  {color_text('3', 'cyan')} = Planar-Averaged 1D Profile")
    mode_choice = get_input("Select mode (1-3) [default: 1]: ").strip()

    # Argumentos base
    args = ["--label", label, "--no-intro"]

    # Lógica para 3D vs 2D vs Perfil
    if mode_choice == '2':
        # Modo 3D
        args.append("--3d")
        print(color_text("\nSelected: Full 3D Volume export.", 'cyan'))

        iso_min_str = get_input(
            "Minimum |density| (e/Ang^3) to keep, for a manageable point cloud (blank = no filter): "
        ).strip()
        if iso_min_str:
            args.extend(["--iso-min", iso_min_str])

        cube_choice = get_input("Also write a Gaussian .cube file for VESTA/VMD/Avogadro? (y/N): ", 'green').strip().lower()
        if cube_choice == 'y':
            args.append("--cube")
    elif mode_choice == '3':
        # Modo Perfil
        args.append("--profile")
        print(color_text("\nSelected: Planar-Averaged 1D Profile.", 'cyan'))

        print(f"Choose the axis the profile varies {color_text('ALONG', 'bold')}:")
        print(f"  {color_text('0', 'cyan')} = X")
        print(f"  {color_text('1', 'cyan')} = Y")
        print(f"  {color_text('2', 'cyan')} = Z")
        axis = get_int_input("Select axis (0-2) [default: 2]: ", 2)
        if axis not in [0, 1, 2]: axis = 2
        args.extend(["--axis", str(axis)])
    else:
        # Modo 2D (Default)
        print(color_text("\nSelected: 2D Slice Configuration", 'cyan'))

        # Escolher Eixo Normal
        print(f"Choose the axis {color_text('NORMAL', 'bold')} to the cut plane:")
        print(f"  {color_text('0', 'cyan')} = X (Cut YZ plane)")
        print(f"  {color_text('1', 'cyan')} = Y (Cut XZ plane)")
        print(f"  {color_text('2', 'cyan')} = Z (Cut XY plane - Standard)")
        axis = get_int_input("Select axis (0-2) [default: 2]: ", 2)
        if axis not in [0, 1, 2]: axis = 2
        args.extend(["--axis", str(axis)])

        # Escolher Posição (Opcional)
        pos_str = get_input("Position in Angstrom (Press Enter for center): ").strip()
        if pos_str:
            args.extend(["--pos", pos_str])

        contour_choice = get_input("Overlay contour lines on the map? (y/N): ", 'green').strip().lower()
        if contour_choice == 'y':
            args.append("--contour")

    # Densidade de spin ou diferença de densidade (mutuamente com o modo acima)
    spin_choice = get_input("Plot the spin density instead of total charge? (y/N): ", 'green').strip().lower()
    if spin_choice == 'y':
        args.append("--spin")

    rho2 = get_input(
        "Second .RHO file to subtract (Delta rho = rho1 - rho2), blank to skip: "
    ).strip()
    if rho2:
        args.extend(["--rho2", rho2])

    # Escala de cor manual (apenas relevante para os modos 2D/3D, não para o perfil)
    if mode_choice != '3':
        vrange_choice = get_input(
            "Fix the colorbar range manually instead of auto? (y/N): ", 'green'
        ).strip().lower()
        if vrange_choice == 'y':
            vmin_str = get_input("Colorbar minimum (e/Ang^3), blank to leave auto: ").strip()
            if vmin_str:
                args.extend(["--vmin", vmin_str])
            vmax_str = get_input("Colorbar maximum (e/Ang^3), blank to leave auto: ").strip()
            if vmax_str:
                args.extend(["--vmax", vmax_str])

    print(color_text("\nProcessing Density...", 'green'))
    run_tool("stb-density", args)


def run_workfunction_calculator() -> None:
    """Interface for the Work Function Calculator (workfunction.py)"""
    print("\n" + "="*60)
    print(color_text("WORK FUNCTION CALCULATOR", 'bold').center(60))
    print("="*60 + "\n")
    
    # 1. Obter o SystemLabel
    label = get_input("Enter the Siesta SystemLabel (e.g., siesta): ").strip()
    while not label:
        print(color_text("Label cannot be empty!", 'red'))
        label = get_input("Enter the Siesta SystemLabel: ")

    # 2. Arquivo de Potencial/Grid
    print(f"\n{color_text('Grid/Potential File:', 'yellow')}")
    print("Default logic: looks for {label}.VT (Total Potential)")
    grid_file = get_input(f"Grid filename (default: {label}.VT): ").strip()
    
    # 3. Arquivo para Energia de Fermi
    print(f"\n{color_text('Fermi Energy Source:', 'yellow')}")
    print("Default logic: looks for Fermi energy in {label}.out")
    fermi_file = get_input(f"Output filename (default: {label}.out): ").strip()
    
    # 4. Eixo de integração
    print(f"\n{color_text('Integration Axis (Planar Average):', 'yellow')}")
    print(f"  {color_text('0', 'cyan')} = x")
    print(f"  {color_text('1', 'cyan')} = y")
    print(f"  {color_text('2', 'cyan')} = z (standard for slabs)")
    axis_str = get_input(
        "Select axis (0-2), blank to auto-detect from label.XV/.fdf if available: "
    ).strip()
    axis_choice = None
    if axis_str:
        try:
            axis_choice = int(axis_str)
            if axis_choice not in (0, 1, 2):
                print(color_text("Invalid axis. Auto-detecting instead.", 'red'))
                axis_choice = None
        except ValueError:
            print(color_text("Invalid axis. Auto-detecting instead.", 'red'))
            axis_choice = None

    # Montar argumentos
    args = ["--label", label, "--no-intro"]
    if axis_choice is not None:
        args.extend(["--axis", str(axis_choice)])

    if grid_file:
        args.extend(["--grid", grid_file])

    if fermi_file:
        # Nota: O workfunction.py usa --fdf para ler o arquivo de saída (Fermi)
        args.extend(["--file", fermi_file])

    print(color_text("\nRunning Work Function analysis...", 'green'))
    run_tool("stb-workfunction", args)


def run_bader_calculator() -> None:
    """Interface for the Bader Charge Analysis (bader.py)"""
    print("\n" + "="*60)
    print(color_text("BADER CHARGE ANALYSIS", 'bold').center(60))
    print("="*60 + "\n")
    
    # 1. Obter o SystemLabel
    label = get_input("Enter the Siesta SystemLabel (e.g., siesta): ").strip()
    while not label:
        print(color_text("Label cannot be empty!", 'red'))
        label = get_input("Enter the Siesta SystemLabel: ")

    # 2. Configurações de Arquivos
    output_file = get_input(f"Output filename (default: {label}_BADER.txt): ").strip()
    
    # --- NOVO: Pergunta sobre o arquivo de referência (.out) ---
    print(f"\n{color_text('Reference Output File (for Z_val detection):', 'cyan')}")
    print("If your .out file has a different name or path, specify it below.")
    print("Otherwise, leave blank to look for the default file.")
    ref_file = get_input(f"Path to .out file (default: {label}.out): ").strip()

    # 3. Configuração de Velocidade
    print(f"\n{color_text('Select speed mode:', 'yellow')}")
    print(f"  {color_text('1.', 'yellow')} Normal (Precise)")
    print(f"  {color_text('2.', 'yellow')} Fast (Less refined edges)")
    speed_choice = get_input("Choice (1-2, default: 1): ", 'green')
    speed_mode = 'fast' if speed_choice == '2' else 'normal'

    # 4. Vácuo (slabs/fios/moléculas)
    vacuum_tol = get_input(
        "Vacuum tolerance in e/Ang^3, for slabs/wires/isolated molecules (blank to disable): "
    ).strip()

    # 5. Threads
    threads_str = get_input("Worker threads (blank for auto -- $SLURM_CPUS_PER_TASK or all cores): ").strip()

    # 6. Arquivo .cube intermediário
    keep_cube_choice = get_input("Keep the intermediate .cube file(s) after the run? (Y/n): ", 'green').strip().lower()

    # 7. Exportar volumes de Bader
    export_volumes_choice = get_input(
        "Export each atom's Bader volume as its own .cube file, for VESTA/VMD? (y/N): ", 'green'
    ).strip().lower()

    # Argumentos básicos
    args = ["--label", label, "--speed", speed_mode, "--no-intro"]

    # Argumentos opcionais
    if output_file:
        args.extend(["--output", output_file])

    # --- NOVO: Adiciona a flag --ref se o usuário digitou algo ---
    if ref_file:
        args.extend(["--ref", ref_file])

    if vacuum_tol:
        args.extend(["--vacuum-tol", vacuum_tol])
    if threads_str:
        args.extend(["--threads", threads_str])
    if keep_cube_choice == 'n':
        args.append("--no-cube")
    if export_volumes_choice == 'y':
        args.append("--export-volumes")

    run_tool("stb-bader", args)
    

def run_elastic_generator() -> None:
    """Interface for the Elastic Constants Generator (elastic_inputs.py)"""
    print("\n" + "="*60)
    print(color_text("ELASTIC CONSTANTS GENERATOR", 'bold').center(60))
    print("="*60 + "\n")
    
    # 1. Obter ficheiro de estrutura
    input_file = get_input("Input structure file (fdf/poscar): ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input structure file: ")
    
    # 2. Obter deformação máxima e passos
    max_strain = get_float_input("\nMax strain % (default: 2.0): ", 2.0)
    steps = get_int_input("Number of steps per direction (default: 4): ", 4)
    
    # 3. MENU DE DIREÇÕES
    print("\n" + "-"*60)
    print(color_text("SELECT DEFORMATION MODE", 'cyan').center(60))
    print("-"*60)
    print(f"[{color_text('1', 'yellow')}] Full 3D Tensor  (xx, yy, zz, xy, xz, yz) -> Standard 3D")
    print(f"[{color_text('2', 'yellow')}] Normal Strains  (xx, yy, zz)             -> Bulk Modulus")
    print(f"[{color_text('3', 'yellow')}] Shear Strains   (xy, xz, yz)             -> Shear Modulus")
    print(f"[{color_text('4', 'yellow')}] 2D In-Plane     (xx, yy, xy)             -> Graphene/Monolayers")
    print(f"[{color_text('5', 'yellow')}] Uniaxial Z-Only (xx)                     -> Nanowires/Tubes")
    print("-" * 60)
    
    mode = get_input("Select mode (1-5): ")
    
    # Mapeamento da escolha para as strings que o elastic_inputs.py entende
    dirs_map = {
        '1': ["xx", "yy", "zz", "xy", "zx", "yz"],
        '2': ["xx", "yy", "zz"],
        '3': ["xy", "zx", "yz"],
        '4': ["xx", "yy", "xy"],
        '5': ["xx"]
    }
    
    # Se a escolha for inválida, assume o padrão (1)
    selected_dirs = dirs_map.get(mode, dirs_map['1'])
    print(f"Selected directions: {color_text(str(selected_dirs), 'green')}\n")
    
    # Monta a lista de argumentos base
    args = [
        "--file", input_file,
        "--max", str(max_strain),
        "--steps", str(steps),
        "--no-intro",
        "--dirs" # Adiciona a flag --dirs
    ]

    # Adiciona as direções escolhidas à lista de argumentos
    args.extend(selected_dirs)

    run_tool("stb-elasticInputs", args)

def run_elastic_analyzer() -> None:
    """Interface for the Elastic Properties Analyzer """
    print("\n" + "="*60)
    print(color_text("ELASTIC PROPERTIES ANALYZER", 'bold').center(60))
    print("="*60 + "\n")
    
    args = []

    # --- NOVO: Solicita o nome do arquivo de output ---
    print(color_text("Enter the Siesta output filename located inside strain folders.", 'yellow'))
    output_filename = get_input("Filename (default: calc.out): ").strip()
    
    if not output_filename:
        output_filename = "calc.out"
    
    # Passa o argumento -f/--file para o script elastic_analysis.py
    args.extend(["--file", output_filename,"--no-intro"])
    print(f"Targeting file: {color_text(output_filename, 'cyan')}\n")
    # --------------------------------------------------
    
    print(f"Is this a {color_text('2D material', 'cyan')}? (affects stiffness units N/m vs GPa)")
    is_2d = get_input("Enable 2D analysis? (y/N): ").lower()
    if is_2d == 'y' or is_2d == 'yes':
        args.append("--2d")
        print(color_text("-> 2D Mode Enabled", 'green'))
    
    print(color_text("\nRunning analysis in current directory...", 'yellow'))
    run_tool("stb-elasticAnalysis", args)

def run_input_generator() -> None:
    """Interface for the Input File Generator (stb-inputfile)"""
    print("\n" + "="*60)
    print(color_text("INPUT FILE GENERATOR (stb-inputfile)", 'bold').center(60))
    print("="*60 + "\n")
    
    # --- Validação do ficheiro de entrada ---
    # Esta linha agora terá Tab-completion!
    input_file = get_input("Input structure file (e.g., struct.fdf): ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input structure file: ")
    
    # --- PONTO 1: Menu Numérico para o tipo de cálculo ---
    mode_list = [
        'total_energy', 'total_energy+d3',
        'relax', 'relax+d3',
        'aimd', 'aimd+d3',
        'bands', 'bands+d3'
    ]
    
    print(f"\n{color_text('Available calculation modes:', 'yellow')}")
    for i, mode in enumerate(mode_list, 1):
        print(f"  {color_text(str(i)+'.', 'yellow')} {mode}")

    choice = 0
    max_choice = len(mode_list)
    
    while not (1 <= choice <= max_choice):
        choice = get_int_input(f"\nSelect calculation mode (1-{max_choice}): ")
        if not (1 <= choice <= max_choice):
            print(color_text(f"Invalid choice! Please select between 1 and {max_choice}.", 'red'))
            
    calc_type = mode_list[choice - 1]
    print(f"Selected mode: {color_text(calc_type, 'cyan')}") 

    # --- PONTO 2: Validação do caminho do Pseudopotencial ---
    
    args = [
        input_file, 
        "--type", calc_type,
        "--no-intro"
    ]
    
    pp_path = prompt_pseudo_source(optional=True)
    if pp_path:
        args.extend(["--pp-path", pp_path])
    else:
        print(color_text("Skipping pseudopotential path.", 'yellow'))

    run_tool("stb-inputfile", args)

def run_kgrid_generator() -> None:
    """Interface for the K-Grid Generator (stb-kgrid)"""
    print("\n" + "="*60)
    print(color_text("K-GRID GENERATOR (stb-kgrid)", 'bold').center(60))
    print("="*60 + "\n")
    
    # Esta linha agora terá Tab-completion!
    input_file = get_input("Input structure file (fdf/poscar/cif/fhi): ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input structure file: ")

    type_list = ['fdf', 'poscar', 'cif', 'fhi']
    print(f"\n{color_text('Available file types:', 'yellow')}")
    for i, t in enumerate(type_list, 1):
        print(f"  {color_text(str(i)+'.', 'yellow')} {t}")

    max_choice = len(type_list)
    type_choice = get_int_input(f"\nSelect file type (1-{max_choice}) [default: 1 = fdf]: ", 1)
    while not (1 <= type_choice <= max_choice):
        print(color_text(f"Invalid choice! Please select between 1 and {max_choice}.", 'red'))
        type_choice = get_int_input(f"Select file type (1-{max_choice}) [default: 1 = fdf]: ", 1)
    file_type = type_list[type_choice - 1]
    print(f"Selected type: {color_text(file_type, 'cyan')}")

    density = get_float_input("\nK-point density (e.g., 0.2) [default: 0.2]: ", 0.2)
    while density <= 0:
        print(color_text("Density must be a positive number!", 'red'))
        density = get_float_input("K-point density (e.g., 0.2) [default: 0.2]: ", 0.2)

    args = [
        "--file", input_file,
        "--type", file_type,
        "--density", str(density),
        "--no-intro"
    ]
    
    run_tool("stb-kgrid", args)

def run_kpath_generator() -> None:
    """Interface for the K-Path Generator (stb-kpath)"""
    print("\n" + "="*60)
    print(color_text("K-PATH GENERATOR (stb-kpath)", 'bold').center(60))
    print("="*60 + "\n")
    
    # Esta linha agora terá Tab-completion!
    input_file = get_input("Input structure file (fdf): ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input structure file: ")

    precision = get_float_input("\nBravais-lattice detection tolerance / eps (default: 0.0002): ", 0.0002)

    args = [
        "--file", input_file,
        "--prec", str(precision),
        "--no-intro"
    ]

    run_tool("stb-kpath", args)

def run_supercell_generator() -> None:
    """Interface for the Supercell Builder (stb-supercell)"""
    print("\n" + "="*60)
    print(color_text("SUPERCELL BUILDER (stb-supercell)", 'bold').center(60))
    print("="*60 + "\n")

    # Esta linha agora terá Tab-completion!
    input_file = get_input("Input structure file (fdf): ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input structure file: ")

    print(f"\n{color_text('Supercell dimensions:', 'yellow')}")
    print("  Enter 3 numbers for a diagonal supercell (e.g. '2 2 2')")
    print("  or 9 numbers for a full row-major 3x3 matrix (e.g. '2 0 0 0 2 0 0 0 1')")
    dim_values = get_input("Dimensions: ").split()
    while len(dim_values) not in (3, 9):
        print(color_text("Please enter exactly 3 or 9 numbers.", 'red'))
        dim_values = get_input("Dimensions: ").split()

    output_file = get_input("\nOutput file name [default: supercell.fdf]: ").strip()
    if not output_file:
        output_file = "supercell.fdf"

    args = [
        "--file", input_file,
        "--dim", *dim_values,
        "--output", output_file,
        "--no-intro"
    ]

    run_tool("stb-supercell", args)

def run_slab_generator() -> None:
    """Interface for the Slab Builder (stb-slab)"""
    print("\n" + "="*60)
    print(color_text("SLAB BUILDER (stb-slab)", 'bold').center(60))
    print("="*60 + "\n")

    input_file = get_input("Input bulk structure file (fdf): ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input bulk structure file (fdf): ")

    hkl_values = []
    while len(hkl_values) != 3 or hkl_values == [0, 0, 0]:
        hkl_input = get_input("\nMiller index, 3 integers (e.g. '1 0 0'): ").split()
        try:
            hkl_values = [int(v) for v in hkl_input]
        except ValueError:
            hkl_values = []
        if len(hkl_values) != 3:
            print(color_text("Please enter exactly 3 integers.", 'red'))
        elif hkl_values == [0, 0, 0]:
            print(color_text("Miller index (0, 0, 0) is not valid.", 'red'))

    min_slab_size = get_float_input("\nMinimum slab thickness in Ang [default: 10.0]: ", 10.0)
    min_vacuum_size = get_float_input("Minimum vacuum thickness in Ang [default: 15.0]: ", 15.0)

    primitive_choice = get_input("\nReduce to primitive cell before cutting? (y/N): ").strip().lower()
    center_choice = get_input("Center the slab in the vacuum? (y/N): ").strip().lower()
    symmetrize_choice = get_input("Try to symmetrize polar/asymmetric terminations? (y/N): ").strip().lower()

    print(f"\n{color_text('Select Termination Mode:', 'yellow')}")
    print(f"  {color_text('1', 'cyan')} = First best termination only (non-polar/symmetric preferred) [Default]")
    print(f"  {color_text('2', 'cyan')} = Interactive (show all terminations, choose one)")
    print(f"  {color_text('3', 'cyan')} = All terminations (writes one file per termination)")
    mode_choice = get_input("Select option (1-3) [default: 1]: ").strip()

    passivate_choice = get_input("\nPassivate dangling bonds on the cut surface? (y/N): ").strip().lower()
    passivant = None
    cutoff_str = ""
    bond_length_str = ""
    if passivate_choice in ('y', 'yes'):
        passivant = get_input("Passivant element [default: H]: ").strip() or "H"
        cutoff_str = get_input("Neighbor cutoff, Ang (blank = auto-detect): ").strip()
        bond_length_str = get_input("Bond length, Ang (blank = auto per species pair): ").strip()

    output_file = get_input("\nOutput file name [default: slab.fdf]: ").strip()
    if not output_file:
        output_file = "slab.fdf"

    args = [
        "--file", input_file,
        "--hkl", *[str(v) for v in hkl_values],
        "--min-slab-size", str(min_slab_size),
        "--min-vacuum-size", str(min_vacuum_size),
        "--output", output_file,
        "--no-intro"
    ]

    if primitive_choice in ('y', 'yes'):
        args.append("--primitive")
    if center_choice in ('y', 'yes'):
        args.append("--center-slab")
    if symmetrize_choice in ('y', 'yes'):
        args.append("--symmetrize")

    if mode_choice == '2':
        args.append("--interactive")
    elif mode_choice == '3':
        args.append("--all")

    if passivant is not None:
        args.extend(["--passivate", "--passivant", passivant])
        if cutoff_str:
            args.extend(["--cutoff", cutoff_str])
        if bond_length_str:
            args.extend(["--bond-length", bond_length_str])

    run_tool("stb-slab", args)


def run_nanotube_generator() -> None:
    """Interface for the Nanotube/Nanoribbon Builder (stb-nanotube)"""
    print("\n" + "="*60)
    print(color_text("NANOTUBE/NANORIBBON BUILDER (stb-nanotube)", 'bold').center(60))
    print("="*60 + "\n")

    input_file = get_input("Input 2D monolayer structure file (fdf): ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input 2D monolayer structure file (fdf): ")

    chirality = []
    while len(chirality) != 2 or chirality == [0, 0]:
        chir_input = get_input("\nChirality indices, 2 integers (e.g. '6 0'): ").split()
        try:
            chirality = [int(v) for v in chir_input]
        except ValueError:
            chirality = []
        if len(chirality) != 2:
            print(color_text("Please enter exactly 2 integers.", 'red'))
        elif chirality == [0, 0]:
            print(color_text("Chirality (0, 0) is not valid.", 'red'))

    print(f"\n{color_text('Select Mode:', 'yellow')}")
    print(f"  {color_text('1', 'cyan')} = Nanotube (roll into a cylinder) [Default]")
    print(f"  {color_text('2', 'cyan')} = Nanoribbon (finite-width flat strip)")
    mode_choice = get_input("Select option (1-2) [default: 1]: ").strip()
    mode = 'ribbon' if mode_choice == '2' else 'tube'

    repeats = get_int_input(
        "\nRepeats (axial for tube, width for ribbon) [default: 1]: ", 1)
    while repeats < 1:
        print(color_text("Repeats must be >= 1!", 'red'))
        repeats = get_int_input("Repeats [default: 1]: ", 1)

    min_vacuum_size = get_float_input("\nVacuum padding around the structure in Ang [default: 15.0]: ", 15.0)

    output_file = get_input("\nOutput file name [default: nanotube.fdf]: ").strip()
    if not output_file:
        output_file = "nanotube.fdf"

    args = [
        "--file", input_file,
        "--chirality", *[str(v) for v in chirality],
        "--mode", mode,
        "--repeats", str(repeats),
        "--min-vacuum-size", str(min_vacuum_size),
        "--output", output_file,
        "--no-intro"
    ]

    run_tool("stb-nanotube", args)


def run_defect_generator() -> None:
    """Interface for the Point Defect Generator (stb-defect)"""
    print("\n" + "="*60)
    print(color_text("POINT DEFECT GENERATOR (stb-defect)", 'bold').center(60))
    print("="*60 + "\n")

    input_file = get_input("Input structure file (fdf): ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input structure file (fdf): ")

    print(f"\n{color_text('Select Defect Type:', 'yellow')}")
    print(f"  {color_text('1', 'cyan')} = Vacancy (remove an atom)")
    print(f"  {color_text('2', 'cyan')} = Substitution (replace an atom's species)")
    print(f"  {color_text('3', 'cyan')} = Interstitial (add a new atom)")
    type_choice = get_input("Select option (1-3): ").strip()
    type_map = {'1': 'vacancy', '2': 'substitution', '3': 'interstitial'}
    defect_type = type_map.get(type_choice, 'vacancy')

    args = ["--file", input_file, "--type", defect_type, "--no-intro"]

    if defect_type in ('vacancy', 'substitution'):
        print(f"\n{color_text('Select Site:', 'yellow')}")
        print(f"  {color_text('1', 'cyan')} = By atom index (e.g. '3,7')")
        print(f"  {color_text('2', 'cyan')} = By nearest position to a target point")
        print(f"  {color_text('3', 'cyan')} = Every symmetrically distinct site (auto, via spglib)")
        site_choice = get_input("Select option (1-3) [default: 1]: ").strip()

        if site_choice == '2':
            coords = get_input("Target position, 3 numbers (e.g. '0.5 0.5 0.5'): ").split()
            fmt_choice = get_input("Format: (f)ractional or (c)artesian [default: f]: ").strip().lower()
            fmt = 'cartesian' if fmt_choice == 'c' else 'fractional'
            filter_sp = get_input("Restrict search to one element (optional, Enter to skip): ").strip()
            args.extend(["--nearest", *coords, "--nearest-format", fmt])
            if filter_sp:
                args.extend(["--filter-species", filter_sp])
        elif site_choice == '3':
            filter_sp = get_input("Restrict to one element (optional, Enter for all species): ").strip()
            args.append("--all-inequivalent-sites")
            if filter_sp:
                args.extend(["--filter-species", filter_sp])
            ml_rank = get_input(
                "Rank sites by MACE-relaxed energy? Needs the optional 'ml' extra (y/N): "
            ).strip().lower()
            if ml_rank in ('y', 'yes'):
                args.append("--ml-rank")
        else:
            index_str = get_input("Atom index/indices, comma-separated (e.g. '3,7'): ").strip()
            args.extend(["--index", index_str])

        if defect_type == 'substitution':
            new_species = get_input("New element symbol: ").strip()
            args.extend(["--new-species", new_species])
    else:
        coords = get_input("Position for the new atom, 3 numbers (e.g. '0.5 0.5 0.5'): ").split()
        fmt_choice = get_input("Format: (f)ractional or (c)artesian [default: f]: ").strip().lower()
        fmt = 'cartesian' if fmt_choice == 'c' else 'fractional'
        species = get_input("Element symbol of the new atom: ").strip()
        args.extend(["--position", *coords, "--position-format", fmt, "--species", species])

    output_file = get_input("\nOutput file name [default: defect.fdf]: ").strip()
    if not output_file:
        output_file = "defect.fdf"
    args.extend(["--output", output_file])

    run_tool("stb-defect", args)


def run_sqs_generator() -> None:
    """Interface for the SQS Generator (stb-sqs)"""
    print("\n" + "="*60)
    print(color_text("SPECIAL QUASIRANDOM STRUCTURE GENERATOR (stb-sqs)", 'bold').center(60))
    print("="*60 + "\n")

    input_file = get_input("Input structure file (fdf): ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input structure file (fdf): ")

    sublattice = get_input("\nSpecies to disorder (e.g. 'Ni'): ").strip()
    composition = get_input("Target composition, e.g. 'Ni:0.5,Fe:0.5': ").strip()
    scaling_str = get_input("Scaling factor (positive integer, blank = auto-detect): ").strip()
    while scaling_str and (not scaling_str.isdigit() or int(scaling_str) < 1):
        print(color_text("Scaling must be blank (auto-detect) or a positive integer!", 'red'))
        scaling_str = get_input("Scaling factor (positive integer, blank = auto-detect): ").strip()
    scaling = int(scaling_str) if scaling_str else None

    print(f"\n{color_text('Select Search Method:', 'yellow')}")
    print(f"  {color_text('1', 'cyan')} = Monte Carlo [Default]")
    print(f"  {color_text('2', 'cyan')} = Enumeration (exact, only feasible for small cells)")
    method_choice = get_input("Select option (1-2) [default: 1]: ").strip()
    method = 'enumeration' if method_choice == '2' else 'monte_carlo'

    temperature = get_float_input("\nMonte Carlo starting temperature [default: 1.0]: ", 1.0)

    output_file = get_input("\nOutput file name [default: sqs.fdf]: ").strip()
    if not output_file:
        output_file = "sqs.fdf"

    args = [
        "--file", input_file,
        "--sublattice", sublattice,
        "--composition", composition,
        "--method", method,
        "--temperature", str(temperature),
        "--output", output_file,
        "--no-intro"
    ]
    if scaling is not None:
        args.extend(["--scaling", str(scaling)])

    run_tool("stb-sqs", args)


def run_unitcell_generator() -> None:
    """Interface for the Unit Cell Finder (stb-unitcell)"""
    print("\n" + "="*60)
    print(color_text("UNIT CELL FINDER (stb-unitcell)", 'bold').center(60))
    print("="*60 + "\n")

    input_file = get_input("Input structure file (fdf): ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input structure file (fdf): ")

    print(f"\n{color_text('Select Mode:', 'yellow')}")
    print(f"  {color_text('1', 'cyan')} = Primitive (smallest cell) [Default]")
    print(f"  {color_text('2', 'cyan')} = Conventional (standardized, usually larger)")
    print(f"  {color_text('3', 'cyan')} = Refined (conventional cell, positions snapped to symmetry)")
    mode_choice = get_input("Select option (1-3) [default: 1]: ").strip()
    mode = {'2': 'conventional', '3': 'refined'}.get(mode_choice, 'primitive')

    symprec = get_float_input("\nSymmetry precision [default: 1e-3]: ", 1e-3)

    output_file = get_input("\nOutput file name [default: unitcell.fdf]: ").strip()
    if not output_file:
        output_file = "unitcell.fdf"

    args = [
        "--file", input_file,
        "--mode", mode,
        "--symprec", str(symprec),
        "--output", output_file,
        "--no-intro"
    ]

    run_tool("stb-unitcell", args)


def run_crystalbuilder_generator() -> None:
    """Interface for the Crystal Builder (stb-crystalbuilder)"""
    print("\n" + "="*60)
    print(color_text("CRYSTAL BUILDER (stb-crystalbuilder)", 'bold').center(60))
    print("="*60 + "\n")

    spacegroup = get_input("Space group (symbol e.g. 'Fm-3m' or number e.g. '225'): ").strip()

    a = get_float_input("\nLattice constant a (Ang): ")
    b_str = get_input("Lattice constant b (Ang) [default: same as a]: ").strip()
    c_str = get_input("Lattice constant c (Ang) [default: same as a]: ").strip()
    alpha = get_float_input("Lattice angle alpha (deg) [default: 90]: ", 90.0)
    beta = get_float_input("Lattice angle beta (deg) [default: 90]: ", 90.0)
    gamma = get_float_input("Lattice angle gamma (deg) [default: 90]: ", 90.0)

    args = [
        "--spacegroup", spacegroup,
        "--a", str(a),
        "--alpha", str(alpha), "--beta", str(beta), "--gamma", str(gamma),
    ]
    if b_str:
        args.extend(["--b", b_str])
    if c_str:
        args.extend(["--c", c_str])

    print(f"\n{color_text('Enter each symmetrically-distinct Wyckoff site:', 'yellow')}")
    print("  (element symbol + fractional x y z, e.g. 'Ni 0 0 0'; blank line to finish)")
    while True:
        site = get_input("Site (blank to finish): ").strip()
        if not site:
            break
        parts = site.split()
        if len(parts) != 4:
            print(color_text("Expected 4 values: SYMBOL X Y Z.", 'red'))
            continue
        args.extend(["--site", *parts])

    if "--site" not in args:
        print(color_text("No sites given -- aborting.", 'red'))
        return

    output_file = get_input("\nOutput file name [default: crystal.fdf]: ").strip()
    if not output_file:
        output_file = "crystal.fdf"
    args.extend(["--output", output_file, "--no-intro"])

    run_tool("stb-crystalbuilder", args)


def run_crystalcast_generator() -> None:
    """Interface for the Random Crystal Generator (stb-crystalcast)"""
    print("\n" + "="*60)
    print(color_text("RANDOM CRYSTAL GENERATOR (stb-crystalcast)", 'bold').center(60))
    print("="*60 + "\n")

    mode = get_input(
        "Mode: generate / analyze / substitute / subgroup / supergroup "
        "[default: generate]: ").strip().lower()

    if mode == "analyze":
        input_file = get_input("Structure file to analyze (.fdf): ").strip()
        if not input_file:
            print(color_text("No file given -- aborting.", 'red'))
            return
        run_tool("stb-crystalcast", ["--analyze", "-f", input_file, "--no-intro"])
        return

    if mode == "substitute":
        input_file = get_input("Structure file to modify (.fdf): ").strip()
        if not input_file:
            print(color_text("No file given -- aborting.", 'red'))
            return
        print(f"\n{color_text('Enter substitutions:', 'yellow')}")
        print("  (OLD:NEW element pairs, e.g. 'Cl:F'; blank line to finish)")
        subs = []
        while True:
            entry = get_input("Substitution (blank to finish): ").strip()
            if not entry:
                break
            subs.append(entry)
        if not subs:
            print(color_text("No substitutions given -- aborting.", 'red'))
            return
        output_file = get_input("Output file name [default: crystalcast.fdf]: ").strip() or "crystalcast.fdf"
        run_tool("stb-crystalcast", ["--substitute", *subs, "-f", input_file,
                                      "-o", output_file, "--no-intro"])
        return

    if mode in ("subgroup", "supergroup"):
        input_file = get_input("Structure file to transform (.fdf): ").strip()
        if not input_file:
            print(color_text("No file given -- aborting.", 'red'))
            return
        args = [f"--{mode}", "-f", input_file]
        if mode == "supergroup":
            while True:
                target = get_input("Target space group number (required): ").strip()
                if not target:
                    print(color_text("--supergroup requires a target space group -- aborting.", 'red'))
                    return
                try:
                    int(target)
                    break
                except ValueError:
                    print(color_text("Please enter a valid space group number.", 'red'))
            args.extend(["--target-group", target])
            d_tol = get_float_input("Max displacement tolerance --d-tol [default: 1.0]: ", 1.0)
            args.extend(["--d-tol", str(d_tol)])
        else:
            while True:
                target = get_input("Target space group number [default: auto-search]: ").strip()
                if not target:
                    break
                try:
                    int(target)
                    args.extend(["--target-group", target])
                    break
                except ValueError:
                    print(color_text("Please enter a valid space group number, or leave blank.", 'red'))
            eps = get_float_input("Perturbation --eps [default: 0.05]: ", 0.05)
            args.extend(["--eps", str(eps)])
        count = get_int_input("Number of candidates to keep [default: 1]: ", 1)
        args.extend(["--count", str(count)])
        output_file = get_input("Output file name [default: crystalcast.fdf]: ").strip() or "crystalcast.fdf"
        args.extend(["--output", output_file, "--no-intro"])
        run_tool("stb-crystalcast", args)
        return

    dim_str = get_input("Dimension: 3=bulk, 2=layer, 1=rod/wire, 0=cluster [default: 3]: ").strip()
    dim = dim_str if dim_str else "3"

    molecular = False
    if dim != "0":
        molecular_str = get_input(
            "Molecular crystal? Packs whole rigid molecules instead of bare atoms (y/n) "
            "[default: n]: ").strip().lower()
        molecular = molecular_str == "y"

    dim_labels = {"3": "space group", "2": "layer group", "1": "rod group", "0": "point group"}
    group_label = dim_labels.get(dim, "symmetry group")
    group = get_input(f"{group_label.capitalize()} (number, or symbol for dim 3/0): ").strip()

    if molecular:
        print(f"\n{color_text('Enter the composition:', 'yellow')}")
        print("  (molecule name from pyxtal's bundled collection, e.g. 'H2O', or a "
              ".xyz/.gjf/.g03/.json path, + how many of it, e.g. 'H2O 4'; blank line to finish)")
        print("  Run stb-crystalcast --list-molecules to see all bundled names.")
    else:
        print(f"\n{color_text('Enter the composition:', 'yellow')}")
        print("  (element symbol + how many atoms of it, e.g. 'Ni 4'; blank line to finish)")
    species = []
    num_ions = []
    while True:
        entry = get_input("Species (blank to finish): ").strip()
        if not entry:
            break
        parts = entry.split()
        if len(parts) != 2:
            print(color_text("Expected 2 values: SYMBOL COUNT.", 'red'))
            continue
        species.append(parts[0])
        num_ions.append(parts[1])

    if not species:
        print(color_text("No species given -- aborting.", 'red'))
        return

    args = [
        "--dim", dim,
        "--group", group,
        "--species", *species,
        "--num-ions", *num_ions,
    ]
    if molecular:
        args.append("--molecular")

    if dim == "2":
        thickness_str = get_input("Layer thickness in Ang [default: automatic]: ").strip()
        if thickness_str:
            args.extend(["--thickness", thickness_str])
    elif dim == "1":
        area_str = get_input("Rod cross-sectional area in Ang^2 [default: automatic]: ").strip()
        if area_str:
            args.extend(["--area", area_str])
    elif dim == "0":
        vacuum = get_float_input("Vacuum padding in Ang [default: 10.0]: ", 10.0)
        args.extend(["--vacuum", str(vacuum)])

    if dim != "0":
        while True:
            lattice_str = get_input(
                "Fix the cell? 'A B C ALPHA BETA GAMMA' [default: estimate from volume factor]: "
            ).strip()
            if not lattice_str:
                break
            tokens = lattice_str.split()
            try:
                if len(tokens) != 6:
                    raise ValueError
                [float(t) for t in tokens]
            except ValueError:
                print(color_text(
                    "Expected exactly 6 numbers: A B C ALPHA BETA GAMMA (or leave blank).", 'red'))
                continue
            args.extend(["--lattice", *tokens])
            break

    sites_str = get_input(
        "Pre-assign Wyckoff sites? One per species, e.g. '4a 4b,4c' [default: fully random]: "
    ).strip()
    if sites_str:
        args.extend(["--sites", *sites_str.split()])

    volume_factor = get_float_input("\nVolume factor [default: 1.1]: ", 1.1)
    args.extend(["--volume-factor", str(volume_factor)])

    count = get_int_input("Number of structures to generate [default: 1]: ", 1)
    args.extend(["--count", str(count)])

    seed_str = get_input("Random seed [default: none]: ").strip()
    if seed_str:
        args.extend(["--seed", seed_str])

    ml_rank_str = get_input(
        "Rank candidates by MACE-MP-0 relaxed energy? Needs the optional 'ml' extra (y/n) "
        "[default: n]: ").strip().lower()
    if ml_rank_str == "y":
        args.append("--ml-rank")

    output_file = get_input("\nOutput file name [default: crystalcast.fdf]: ").strip()
    if not output_file:
        output_file = "crystalcast.fdf"
    args.extend(["--output", output_file, "--no-intro"])

    run_tool("stb-crystalcast", args)


def run_fetch_generator() -> None:
    """Interface for the Structure Fetcher (stb-fetch)"""
    print("\n" + "="*60)
    print(color_text("STRUCTURE FETCHER (stb-fetch)", 'bold').center(60))
    print("="*60 + "\n")

    # Keep in sync with fetch.py's CURATED_OPTIMADE_PROVIDERS.
    curated_optimade_providers = [
        ("twodmatpedia", "2D Materials Encyclopedia -- confirmed working"),
        ("jarvis", "NIST JARVIS-DFT (2D + 3D materials)"),
        ("aflow", "AFLOW"),
        ("oqmd", "Open Quantum Materials Database"),
        ("alexandria.pbe", "Alexandria (PBE)"),
        ("odbx", "Open Database of Xtals"),
        ("mcloud.mc2d", "Materials Cloud 2D Structures Database"),
        ("mcloud.2dtopo", "Materials Cloud 2D Topological Database"),
        ("cmr", "DTU Computational Materials Repository, incl. C2DB -- currently unreliable"),
    ]

    print(f"{color_text('Select Source:', 'yellow')}")
    print(f"  {color_text('1', 'cyan')} = Crystallography Open Database (COD) -- no account needed")
    print(f"  {color_text('2', 'cyan')} = Materials Project -- needs a free API key")
    print(f"  {color_text('3', 'cyan')} = OPTIMADE -- many databases (2D materials, C2DB-adjacent, ...)")
    source_choice = get_input("Select option (1-3) [default: 1]: ").strip()
    source = {'2': 'materials-project', '3': 'optimade'}.get(source_choice, 'cod')

    args = ["--source", source]

    if source == "optimade":
        print(f"\n{color_text('Select Provider:', 'yellow')}")
        for i, (alias, desc) in enumerate(curated_optimade_providers):
            print(f"  {color_text(str(i + 1), 'cyan')} = {alias} -- {desc}")
        print(f"  {color_text('0', 'cyan')} = Other (type an alias or URL)")
        provider_choice = get_input(f"Select option (0-{len(curated_optimade_providers)}) [default: 1]: ").strip()
        if provider_choice == "0":
            provider = get_input("Provider alias or OPTIMADE base URL: ").strip()
        elif provider_choice.isdigit() and 1 <= int(provider_choice) <= len(curated_optimade_providers):
            provider = curated_optimade_providers[int(provider_choice) - 1][0]
        else:
            provider = curated_optimade_providers[0][0]
        args.extend(["--provider", provider])

    print(f"\n{color_text('Select Query:', 'yellow')}")
    print(f"  {color_text('1', 'cyan')} = By exact id")
    print(f"  {color_text('2', 'cyan')} = By chemical formula")
    query_choice = get_input("Select option (1-2) [default: 1]: ").strip()

    if query_choice == "2":
        formula = get_input("Formula (e.g. 'Fe3O4'): ").strip()
        args.extend(["--formula", formula])
        if source == "materials-project":
            most_stable = get_input("Auto-pick the most stable match? (y/n) [default: n]: ").strip().lower()
            if most_stable == "y":
                args.append("--most-stable")
    else:
        if source == "materials-project":
            material_id = get_input("Materials Project id (e.g. 'mp-19306'): ").strip()
            args.extend(["--material-id", material_id])
        elif source == "optimade":
            optimade_id = get_input("Id within the chosen provider (e.g. '2dm-2127'): ").strip()
            args.extend(["--optimade-id", optimade_id])
        else:
            cod_id = get_input("COD numeric id (e.g. '1010369'): ").strip()
            args.extend(["--cod-id", cod_id])

    if source == "materials-project":
        api_key = get_input("API key (blank = use PMG_MAPI_KEY env var): ").strip()
        if api_key:
            args.extend(["--api-key", api_key])

    print(f"\n{color_text('Reduce to unit cell?', 'yellow')}")
    print(f"  {color_text('0', 'cyan')} = No, keep the structure as fetched [Default]")
    print(f"  {color_text('1', 'cyan')} = Primitive (smallest cell)")
    print(f"  {color_text('2', 'cyan')} = Conventional (standardized, usually larger)")
    print(f"  {color_text('3', 'cyan')} = Refined (conventional cell, positions snapped to symmetry)")
    unitcell_choice = get_input("Select option (0-3) [default: 0]: ").strip()
    unitcell_map = {'1': 'primitive', '2': 'conventional', '3': 'refined'}
    if unitcell_choice in unitcell_map:
        args.extend(["--unitcell", unitcell_map[unitcell_choice]])

    output_file = get_input("\nOutput file name [default: fetched.fdf]: ").strip()
    if not output_file:
        output_file = "fetched.fdf"
    args.extend(["--output", output_file, "--no-intro"])

    run_tool("stb-fetch", args)


def run_passivate_generator() -> None:
    """Interface for the Surface Passivator (stb-passivate)"""
    print("\n" + "="*60)
    print(color_text("SURFACE PASSIVATOR (stb-passivate)", 'bold').center(60))
    print("="*60 + "\n")

    input_file = get_input("Input structure file (fdf, typically a cut slab): ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input structure file (fdf): ")

    passivant = get_input("\nPassivant element [default: H]: ").strip()
    if not passivant:
        passivant = "H"

    cutoff_str = get_input("Neighbor cutoff, Ang (blank = auto-detect): ").strip()
    bond_length_str = get_input("Bond length, Ang (blank = auto per species pair): ").strip()

    output_file = get_input("\nOutput file name [default: passivated.fdf]: ").strip()
    if not output_file:
        output_file = "passivated.fdf"

    args = ["--file", input_file, "--passivant", passivant, "--output", output_file, "--no-intro"]
    if cutoff_str:
        args.extend(["--cutoff", cutoff_str])
    if bond_length_str:
        args.extend(["--bond-length", bond_length_str])

    run_tool("stb-passivate", args)


def run_molecule_generator() -> None:
    """Interface for the Reference Molecule Builder (stb-molecule)"""
    print("\n" + "="*60)
    print(color_text("REFERENCE MOLECULE BUILDER (stb-molecule)", 'bold').center(60))
    print("="*60 + "\n")

    print("Names are case-sensitive (e.g. 'H2O', not 'h2o'). Type 'list' to see all names.")
    name = get_input("Molecule name: ").strip()
    while name.lower() == "list" or not name:
        run_tool("stb-molecule", ["--list", "--no-intro"])
        name = get_input("Molecule name: ").strip()

    vacuum = get_float_input("\nVacuum padding in Ang [default: 10.0]: ", 10.0)

    output_file = get_input("\nOutput file name [default: molecule.fdf]: ").strip()
    if not output_file:
        output_file = "molecule.fdf"

    args = ["--name", name, "--vacuum", str(vacuum), "--output", output_file, "--no-intro"]

    run_tool("stb-molecule", args)


def run_mlrelax_generator() -> None:
    """Interface for the ML Pre-Relaxation (stb-mlrelax)"""
    print("\n" + "="*60)
    print(color_text("ML PRE-RELAXATION (stb-mlrelax)", 'bold').center(60))
    print("="*60 + "\n")
    print(color_text(
        "Needs the optional 'ml' extra (pip install stb_suite[ml] -- PyTorch + "
        "mace-torch). Fast heuristic pre-relax, not a substitute for a real DFT "
        "relaxation.", 'yellow'))

    input_file = get_input("\nInput structure file (fdf): ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input structure file (fdf): ")

    relax_cell = get_input(
        "\nAlso relax the cell? Only for bulk (no vacuum) structures (y/N): "
    ).strip().lower()

    model = get_input("Model size, small/medium/large [default: small]: ").strip()
    if not model:
        model = "small"

    fmax = get_float_input("Force convergence, eV/Ang [default: 0.05]: ", 0.05)

    output_file = get_input("\nOutput file name [default: relaxed.fdf]: ").strip()
    if not output_file:
        output_file = "relaxed.fdf"

    args = [
        "--file", input_file,
        "--model", model,
        "--fmax", str(fmax),
        "--output", output_file,
        "--no-intro"
    ]
    if relax_cell in ('y', 'yes'):
        args.append("--relax-cell")

    run_tool("stb-mlrelax", args)


def run_amorphize_generator() -> None:
    """Interface for the Amorphous Structure Generator (stb-amorphize)"""
    print("\n" + "="*60)
    print(color_text("AMORPHOUS STRUCTURE GENERATOR (stb-amorphize)", 'bold').center(60))
    print("="*60 + "\n")
    print(color_text(
        "Needs the optional 'ml' extra (pip install stb_suite[ml]). Melt-quench MD "
        "with MACE-MP-0 -- a fast heuristic starting guess, bulk (3D periodic) "
        "structures only. Can take a few minutes.", 'yellow'))

    input_file = get_input("\nInput structure file (fdf, typically a supercell): ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input structure file (fdf): ")

    melt_temp = get_float_input("\nMelt temperature, K [default: 3000]: ", 3000.0)
    melt_steps = get_int_input("Melt steps [default: 500]: ", 500)
    quench_temp = get_float_input("Quench target temperature, K [default: 300]: ", 300.0)
    quench_steps = get_int_input("Quench steps [default: 1000]: ", 1000)

    final_relax = get_input(
        "\nDo a final static relax after the quench? (Y/n): "
    ).strip().lower()

    output_file = get_input("\nOutput file name [default: amorphous.fdf]: ").strip()
    if not output_file:
        output_file = "amorphous.fdf"

    args = [
        "--file", input_file,
        "--melt-temp", str(melt_temp),
        "--melt-steps", str(melt_steps),
        "--quench-temp", str(quench_temp),
        "--quench-steps", str(quench_steps),
        "--output", output_file,
        "--no-intro"
    ]
    if final_relax in ('n', 'no'):
        args.append("--no-final-relax")

    run_tool("stb-amorphize", args)


def run_convergence_generator() -> None:
    """Interface for the Convergence Test Prep (stb-convergence)"""
    print("\n" + "="*60)
    print(color_text("CONVERGENCE TEST GENERATOR", 'bold').center(60))
    print("="*60 + "\n")

    struct_file = get_input("Input structure file (-s): ")
    while not os.path.isfile(struct_file):
        print(color_text("File not found!", 'red'))
        struct_file = get_input("Input structure file (-s): ")

    calc_file = get_input("Calc.fdf template file (-c): ")
    while not os.path.isfile(calc_file):
        print(color_text("File not found!", 'red'))
        calc_file = get_input("Calc.fdf template file (-c): ")

    print(f"\n{color_text('Select Parameter to Sweep:', 'yellow')}")
    print(f"  {color_text('1', 'cyan')} = Mesh.CutOff (Ry)")
    print(f"  {color_text('2', 'cyan')} = PAO.EnergyShift (Ry)")
    print(f"  {color_text('3', 'cyan')} = K-grid density (1/Ang)")
    param_choice = get_input("Select option (1-3): ").strip()
    param_map = {'1': 'meshcutoff', '2': 'energyshift', '3': 'kgrid'}
    parameter = param_map.get(param_choice, 'meshcutoff')

    min_value = get_float_input("\nMinimum value: ")
    max_value = get_float_input("Maximum value: ")
    step_value = get_float_input("Step: ")

    args = [
        "--structure", struct_file,
        "--calc", calc_file,
        "--parameter", parameter,
        "--min", str(min_value),
        "--max", str(max_value),
        "--step", str(step_value),
        "--no-intro"
    ]

    if parameter == 'kgrid':
        vacuum_gap = get_float_input("\nVacuum-axis detection threshold in Ang [default: 10.0]: ", 10.0)
        args.extend(["--vacuum-gap", str(vacuum_gap)])

    output_dir = get_input("\nOutput directory [default: convergence_runs]: ").strip()
    if not output_dir:
        output_dir = "convergence_runs"
    args.extend(["--output-dir", output_dir])

    run_tool("stb-convergence", args)


def run_convergence_analyzer() -> None:
    """Interface for the Convergence Test Analysis (stb-convergenceAnalysis)"""
    print("\n" + "="*60)
    print(color_text("CONVERGENCE TEST ANALYZER", 'bold').center(60))
    print("="*60 + "\n")

    run_dir = get_input("Directory with 'convergence_*' folders [default: convergence_runs]: ").strip()
    if not run_dir:
        run_dir = "convergence_runs"

    output_filename = get_input("SIESTA output filename inside each folder [default: calc.out]: ").strip()
    if not output_filename:
        output_filename = "calc.out"

    tolerance = get_float_input("Convergence tolerance in eV/atom [default: 0.001]: ", 0.001)

    args = [
        "--dir", run_dir,
        "--file", output_filename,
        "--tolerance", str(tolerance),
        "--no-intro"
    ]

    run_tool("stb-convergenceAnalysis", args)


def run_hubbardu_prep() -> None:
    """Interface for the Hubbard U Linear-Response Reference Prep (stb-hubbardu)"""
    print("\n" + "="*60)
    print(color_text("HUBBARD U (LINEAR RESPONSE) - STAGE 1: REFERENCE PREP", 'bold').center(60))
    print("="*60 + "\n")

    struct_file = get_input("Input structure file (-s): ")
    while not os.path.isfile(struct_file):
        print(color_text("File not found!", 'red'))
        struct_file = get_input("Input structure file (-s): ")

    calc_file = get_input("Calc.fdf template file (-c): ")
    while not os.path.isfile(calc_file):
        print(color_text("File not found!", 'red'))
        calc_file = get_input("Calc.fdf template file (-c): ")

    species = get_input("Species to correct (e.g. Mn): ").strip()
    while not species:
        print(color_text("Species cannot be empty!", 'red'))
        species = get_input("Species to correct (e.g. Mn): ").strip()

    args = ["--structure", struct_file, "--calc", calc_file, "--species", species, "--no-intro"]

    atom_index = get_input(
        "Atom index to perturb (1-based, counting all atoms in the structure) "
        "[optional -- only needed if the species appears more than once]: "
    ).strip()
    if atom_index:
        args.extend(["--atom-index", atom_index])

    print(f"\n{color_text('Correlated shell:', 'yellow')} (blank = auto-detect from species)")
    shell = get_input("  Shell (3d/4d/5d/4f/5f): ").strip()
    if shell:
        args.extend(["--shell", shell])

    j_value = get_float_input("Exchange J (eV) [default: 0.0]: ", 0.0)
    args.extend(["--j", str(j_value)])

    pseudo_dir = prompt_pseudo_source(optional=True)
    if pseudo_dir:
        args.extend(["--pseudo-dir", pseudo_dir])

    output_dir = get_input("\nOutput directory [default: hubbardu_runs]: ").strip()
    if not output_dir:
        output_dir = "hubbardu_runs"
    args.extend(["--output-dir", output_dir])

    run_tool("stb-hubbardu", args)


def run_hubbardu_alphas() -> None:
    """Interface for the Hubbard U Linear-Response Perturbation Prep (stb-hubbarduAlphas)"""
    print("\n" + "="*60)
    print(color_text("HUBBARD U (LINEAR RESPONSE) - STAGE 2: PERTURBATION PREP", 'bold').center(60))
    print("="*60 + "\n")

    run_dir = get_input("Directory with stage 1's 'reference/' folder [default: hubbardu_runs]: ").strip()
    if not run_dir:
        run_dir = "hubbardu_runs"

    args = ["--dir", run_dir, "--no-intro"]

    alphas_str = get_input(
        "Perturbation strengths in eV, space-separated "
        "[default: -0.10 -0.05 0.05 0.10]: "
    ).strip()
    if alphas_str:
        args.extend(["--alphas"] + alphas_str.split())

    frozen_iter = get_int_input(
        "MaxSCFIterations for the frozen-density runs [default: 2]: ", 2
    )
    args.extend(["--frozen-iterations", str(frozen_iter)])

    run_tool("stb-hubbarduAlphas", args)


def run_hubbardu_analysis() -> None:
    """Interface for the Hubbard U Linear-Response Analysis (stb-hubbarduAnalysis)"""
    print("\n" + "="*60)
    print(color_text("HUBBARD U (LINEAR RESPONSE) - ANALYSIS", 'bold').center(60))
    print("="*60 + "\n")

    run_dir = get_input("Directory with the stb-hubbardu run folders [default: hubbardu_runs]: ").strip()
    if not run_dir:
        run_dir = "hubbardu_runs"

    output_filename = get_input("SIESTA output filename inside each folder [default: calc.out]: ").strip()
    if not output_filename:
        output_filename = "calc.out"

    args = ["--dir", run_dir, "--file", output_filename, "--no-intro"]

    r2_tolerance = get_float_input("R^2 tolerance for the fit warnings [default: 0.98]: ", 0.98)
    args.extend(["--r2-tolerance", str(r2_tolerance)])

    output_file = get_input(
        "Output .fdf snippet filename [optional, press Enter for <species>_LDAU.fdf]: "
    ).strip()
    if output_file:
        args.extend(["--output", output_file])

    run_tool("stb-hubbarduAnalysis", args)


def run_dftu_generator() -> None:
    """Interface for the DFT+U / Hubbard Block Generator (stb-dftu)"""
    print("\n" + "="*60)
    print(color_text("DFT+U / HUBBARD BLOCK GENERATOR (stb-dftu)", 'bold').center(60))
    print("="*60 + "\n")

    print(f"{color_text('What do you want to do?', 'yellow')}")
    print(f"  {color_text('1', 'cyan')} = Generate a %block from the Materials Project table (verify a .fdf) (default)")
    print(f"  {color_text('2', 'cyan')} = Generate a %block LDAU.proj by hand (type each species' U yourself)")
    mode_choice = get_input("Select option (1-2) [default: 1]: ").strip() or '1'

    if mode_choice == '1':
        fdf_path = get_input(".fdf file to read species from: ").strip()
        while not os.path.isfile(fdf_path):
            print(color_text("File not found!", 'red'))
            fdf_path = get_input(".fdf file to read species from: ").strip()

        args = ["--fdf", fdf_path, "--use-reference", "--no-intro"]
        output = get_input("\nAlso save to a file? [optional, press Enter to skip]: ").strip()
        if output:
            args.extend(["--output", output])
        run_tool("stb-dftu", args)
        return

    print(f"\n{color_text('Enter one species at a time (blank species to finish):', 'yellow')}")
    species, u_values, j_values, shells = [], [], [], []
    while True:
        sp = get_input("Species (blank to finish): ").strip()
        if not sp:
            break
        u = get_float_input(f"  U (eV) for {sp}: ")
        j = get_float_input(f"  J (eV) for {sp} [default: 0.0]: ", 0.0)
        shell = get_input(f"  Shell for {sp} (3d/4d/5d/4f/5f) [blank = auto-detect]: ").strip()

        species.append(sp)
        u_values.append(str(u))
        j_values.append(str(j))
        if shell:
            shells.append(shell)

    if not species:
        print(color_text("No species given -- aborting.", 'red'))
        return

    args = ["--species", *species, "--u", *u_values, "--j", *j_values, "--no-intro"]
    if shells:
        if len(shells) != len(species):
            print(color_text(
                "Shell given for some but not all species -- skipping --shell "
                "(auto-detecting for every species instead).", 'yellow'))
        else:
            args.extend(["--shell", *shells])

    output = get_input("\nAlso save to a file? [optional, press Enter to skip]: ").strip()
    if output:
        args.extend(["--output", output])

    run_tool("stb-dftu", args)


def run_dos_parser() -> None:
    """Interface for the PDOS XML Parser (stb-dos)"""
    print("\n" + "="*60)
    print(color_text("PDOS XML PARSER (stb-dos)", 'bold').center(60))
    print("="*60 + "\n")
    
    # Esta linha agora terá Tab-completion!
    input_file = get_input("Input PDOS.xml file: ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input PDOS.xml file: ")

    type_list = ['total', 'atom', 'species']
    print(f"\n{color_text('Available types:', 'yellow')} {', '.join(type_list)}")
    dos_types_str = get_input(f"DOS types (space-separated, default: total atom species): ")
    if not dos_types_str.strip():
        dos_types = ['total', 'atom', 'species']
    else:
        dos_types = dos_types_str.split()
        if not all(t in type_list for t in dos_types):
            print(color_text("Input contains invalid types. Using default.", 'yellow'))
            dos_types = ['total', 'atom', 'species']

    shift = get_input("\nEnergy shift ('fermi', '0.0', or a number, default: fermi): ")
    if not shift.strip():
        shift = 'fermi'
    
    print(f"\n{color_text('Select projection mode:', 'yellow')}")
    print(f"  {color_text('1.', 'yellow')} l (s, p, d, f) [Default]")
    print(f"  {color_text('2.', 'yellow')} ml (s, px, py, pz, dxy...)")

    choice = 0
    while not (1 <= choice <= 2):
        # Usamos get_int_input com default = 1
        choice = get_int_input(f"\nSelect mode (1-2) [default: 1]: ", 1) 
        if not (1 <= choice <= 2):
            print(color_text(f"Invalid choice! Please select 1 or 2.", 'red'))
            
    projection_mode = 'l' if choice == 1 else 'ml'
    print(f"Selected mode: {color_text(projection_mode, 'cyan')}")

    output_dir = get_input("\nOutput directory [default: current directory]: ").strip()
    if not output_dir:
        output_dir = "."

    args = [
        input_file, # Positional argument
        "--shift", shift,
        "--type"
    ]

    args.extend(dos_types)
    args.extend(["--projection", projection_mode])
    args.extend(["--output-dir", output_dir])
    args.append("--no-intro")

    run_tool("stb-dos", args)

def run_strain_generator() -> None:
    """Interface for the Strain Generator (stb-strain)"""
    print("\n" + "="*60)
    print(color_text("STRAIN GENERATOR (stb-strain)", 'bold').center(60))
    print("="*60 + "\n")
    
    # Esta linha agora terá Tab-completion!
    input_file = get_input("Input FDF file: ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input FDF file: ")
    
    direction = get_input("Strain direction (x,y,z,xy,xz,yz): ").lower()
    while not all(c in 'xyz' for c in direction) or len(direction) not in (1, 2):
        print(color_text("Invalid direction! Use x,y,z for uniaxial or xy,xz,yz for biaxial", 'red'))
        direction = get_input("Strain direction (x,y,z,xy,xz,yz): ").lower()
    
    stmin = get_float_input("Minimum strain % (default 0): ", 0.0)
    stmax = get_float_input("Maximum strain % (default 25): ", 25.0)
    while stmax <= stmin:
        print(color_text("Maximum strain must be greater than minimum strain!", 'red'))
        stmax = get_float_input("Maximum strain % (default 25): ", 25.0)
    
    step = get_float_input("Step % (default 1): ", 1.0)
    while step <= 0:
        print(color_text("Step must be positive!", 'red'))
        step = get_float_input("Step % (default 1): ", 1.0)
    
    args = [
        "--file", input_file,
        "--stdir", direction,
        "--stmin", str(stmin),
        "--stmax", str(stmax),
        "--step", str(step),
        "--no-intro"
    ]
    
    run_tool("stb-strain", args)

def run_strain_post_processor() -> None:
    """Interface for the Strain Post-Processing (strain_analysis.py)"""
    print("\n" + "="*60)
    print(color_text("STRAIN POST-PROCESSING ANALYZER", 'bold').center(60))
    print("="*60 + "\n")
    print(color_text("This tool analyzes 'strain_*' folders in the current directory.", 'yellow'))
    
    args = []

    # 1. Pergunta qual o nome do ficheiro de output
    print(color_text("Enter the Siesta output filename located inside strain folders.", 'yellow'))
    siesta_out = get_input("Filename (e.g., calc.out): ").strip()
    while not siesta_out:
        print(color_text("Filename is required!", 'red'))
        siesta_out = get_input("Filename (e.g., calc.out): ").strip()
        
    
    # Adiciona argumentos obrigatórios
    args.extend(["--file", siesta_out, "--no-intro"])
    
     
    # 3. --- NOVO: Opção 2D ---
    print(f"\nIs this a {color_text('2D material', 'cyan')}? (Calculates units in N/m)")
    is_2d = get_input("Enable 2D analysis? (y/N): ").lower()
    if is_2d == 'y' or is_2d == 'yes':
        args.append("--2d")
        print(color_text("-> 2D Mode Enabled (N/m)", 'green'))
    # -------------------------

    print(color_text("\nRunning analysis...", 'yellow'))
    run_tool("stb-strainAnalysis", args)


def run_bands_analyzer() -> None:
    """Interface for the Bands Analyzer (stb-bands)"""
    print("\n" + "="*60)
    print(color_text("BANDS ANALYZER (stb-bands)", 'bold').center(60))
    print("="*60 + "\n")
    
    # Esta linha agora terá Tab-completion!
    label = get_input("SIESTA label (e.g. 'siesta' for siesta.bands [+ siesta.EIG if present]): ")
    while not os.path.isfile(f"{label}.bands"):
        print(color_text(f"File '{label}.bands' not found!", 'red'))
        label = get_input("SIESTA label: ")
    if os.path.isfile(f"{label}.EIG"):
        print(color_text(f"-> Found '{label}.EIG', mesh (k-grid) gap comparison will be included.", 'green'))
    else:
        print(color_text(f"-> No '{label}.EIG' found, mesh (k-grid) gap comparison will be skipped.", 'yellow'))

    shift_options = {
        '1': ('vbm', "Valence Band Maximum"),
        '2': ('cbm', "Conduction Band Minimum"),
        '3': ('fermi', "Fermi level"),
        '4': ('manual', "Custom value")
    }
    
    print("\nEnergy reference options:")
    for key, (_, desc) in shift_options.items():
        print(f" {color_text(key, 'yellow')}. {desc}")
    
    choice = get_input("\nSelect reference (1-4): ")
    while choice not in shift_options:
        print(color_text("Invalid choice!", 'red'))
        choice = get_input("Select reference (1-4): ")
    
    shift_type, _ = shift_options[choice]
    args = ["--label", label, "--shift", shift_type, "--no-intro"]

    if shift_type == "manual":
        manual_value = get_float_input("Enter custom shift value: ")
        args.extend(["--manual-value", str(manual_value)])
    
    run_tool("stb-bands", args)

def run_dos_convolution() -> None:
    """Interface for the DOS Processor (Convolution) (stb-convdos)"""
    print("\n" + "="*60)
    print(color_text("DOS PROCESSOR (CONVOLUTION) (stb-convdos)", 'bold').center(60))
    print("="*60 + "\n")
    
    # Esta linha agora terá Tab-completion!
    input_file = get_input("Input DOS file (e.g., dos_total.dat): ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input DOS file: ")
    
    # Esta linha agora terá Tab-completion!
    out_file = get_input("Output file (default: dos_filtered.dat): ", 'green') or "dos_filtered.dat"

    broadening_choice = get_input("Specify broadening as (1) sigma or (2) FWHM, in meV [default: 1]: ").strip()
    if broadening_choice == '2':
        fwhm = get_float_input("FWHM, in meV (default: 118): ", 118.0)
        while fwhm <= 0:
            print(color_text("FWHM must be positive!", 'red'))
            fwhm = get_float_input("FWHM in meV (default: 118): ", 118.0)
        broadening_flag, broadening_value = "--fwhm", fwhm
    else:
        sigma = get_float_input("Sigma, in meV (default: 50): ", 50.0)
        while sigma <= 0:
            print(color_text("Sigma must be positive!", 'red'))
            sigma = get_float_input("Sigma in meV (default: 50): ", 50.0)
        broadening_flag, broadening_value = "--sigma", sigma

    size_str = get_input("Kernel size in samples (optional, blank = auto-sized from broadening): ").strip()

    plot_choice = get_input("Show before/after plot? (Y/n): ").strip().lower()

    args = [
        "--file", input_file,
        broadening_flag, str(broadening_value),
        "--out", out_file,
        "--no-intro"
    ]

    if size_str:
        args.extend(["--size", size_str])
    if plot_choice in ['n', 'no']:
        args.append("--no-plot")

    run_tool("stb-convdos", args)

def run_structure_analyzer() -> None:
    """Interface for the Structure Analyzer (stb-structural)"""
    print("\n" + "="*60)
    print(color_text("STRUCTURE ANALYZER (stb-structural)", 'bold').center(60))
    print("="*60 + "\n")
    
    # Esta linha agora terá Tab-completion!
    input_file = get_input("Input structure file (SIESTA .fdf or .STRUCT_OUT): ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input structure file: ")

    formats = ['fdf', 'struct_out']
    print(f"\n{color_text('Select input file format:', 'yellow')}")
    for i, fmt in enumerate(formats, 1):
        print(f"  {color_text(str(i)+'.', 'yellow')} {fmt}")
    choice = 0
    while not (1 <= choice <= len(formats)):
        choice = get_int_input(f"\nSelect format (1-{len(formats)}): ")
        if not (1 <= choice <= len(formats)):
            print(color_text(f"Invalid choice! Please select between 1 and {len(formats)}.", 'red'))
    format_type = formats[choice - 1]
    print(f"Selected format: {color_text(format_type, 'cyan')}")

    modes = ['mean', 'list']
    print(f"\n{color_text('Select analysis mode:', 'yellow')}")
    for i, m in enumerate(modes, 1):
        print(f"  {color_text(str(i)+'.', 'yellow')} {m}")
    choice = 0
    while not (1 <= choice <= len(modes)):
        choice = get_int_input(f"\nSelect mode (1-{len(modes)}): ")
        if not (1 <= choice <= len(modes)):
            print(color_text(f"Invalid choice! Please select between 1 and {len(modes)}.", 'red'))
    mode = modes[choice - 1]

    output_dir = get_input("Output directory [default: current directory]: ").strip() or "."

    args = ["--file", input_file, "--mode", mode, "--format", format_type, "--output-dir", output_dir, "--no-intro"]

    if mode == "list":
        atom_list = get_input("Enter atom indices (comma-separated, e.g. 1,4,5): ")
        args.extend(["--list", atom_list])

    rdf_choice = get_input("Compute radial distribution function g(r)? (Y/n): ").strip().lower()
    if rdf_choice in ['n', 'no']:
        args.append("--no-rdf")
    else:
        rdf_rmax = get_float_input("RDF cutoff radius in Ang (default: 10.0): ", 10.0)
        while rdf_rmax <= 0:
            print(color_text("Cutoff radius must be a positive number!", 'red'))
            rdf_rmax = get_float_input("RDF cutoff radius in Ang (default: 10.0): ", 10.0)
        args.extend(["--rdf-rmax", str(rdf_rmax)])

    run_tool("stb-structural", args)

def run_xrd_analyzer() -> None:
    """Interface for the Powder XRD Simulator (stb-xrd)"""
    print("\n" + "="*60)
    print(color_text("POWDER XRD SIMULATOR (stb-xrd)", 'bold').center(60))
    print("="*60 + "\n")

    input_file = get_input("Input structure file (SIESTA .fdf or .STRUCT_OUT): ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input structure file: ")

    formats = ['fdf', 'struct_out']
    print(f"\n{color_text('Select input file format:', 'yellow')}")
    for i, fmt in enumerate(formats, 1):
        print(f"  {color_text(str(i)+'.', 'yellow')} {fmt}")
    choice = 0
    while not (1 <= choice <= len(formats)):
        choice = get_int_input(f"\nSelect format (1-{len(formats)}): ")
        if not (1 <= choice <= len(formats)):
            print(color_text(f"Invalid choice! Please select between 1 and {len(formats)}.", 'red'))
    format_type = formats[choice - 1]

    wavelength = get_input("X-ray source name (CuKa, MoKa, ...) or wavelength in Ang [default: CuKa]: ").strip()
    if not wavelength:
        wavelength = "CuKa"

    two_theta_min = get_float_input("2-theta range minimum in deg [default: 0]: ", 0.0)
    two_theta_max = get_float_input("2-theta range maximum in deg [default: 90]: ", 90.0)

    top_str = get_input("Show only the N strongest peaks [default: show all]: ").strip()

    compare_to = get_input(
        "Compare to an experimental pattern file (2theta, intensity columns) [default: skip]: "
    ).strip()
    while compare_to and not os.path.isfile(compare_to):
        print(color_text("File not found!", 'red'))
        compare_to = get_input(
            "Compare to an experimental pattern file [default: skip]: ").strip()

    plot_choice = get_input("Show an interactive plot? (y/N): ").strip().lower()

    output_file = get_input("Output data file name [default: xrd_pattern.dat]: ").strip() or "xrd_pattern.dat"

    args = ["--file", input_file, "--format", format_type, "--wavelength", wavelength,
            "--two-theta-range", str(two_theta_min), str(two_theta_max),
            "--output", output_file, "--no-intro"]
    if top_str:
        args.extend(["--top", top_str])
    if compare_to:
        args.extend(["--compare-to", compare_to])
    if plot_choice == "y":
        args.append("--plot")

    run_tool("stb-xrd", args)

def run_xrdsearch_generator() -> None:
    """Interface for the XRD-Guided Structure Search prep stage (stb-xrdsearch)"""
    print("\n" + "="*60)
    print(color_text("XRD STRUCTURE SEARCH -- PREP (stb-xrdsearch)", 'bold').center(60))
    print("="*60 + "\n")

    print(f"{color_text('Enter the composition:', 'yellow')}")
    print("  (element symbol + how many atoms of it, e.g. 'Ni 4'; blank line to finish)")
    species = []
    num_ions = []
    while True:
        entry = get_input("Species (blank to finish): ").strip()
        if not entry:
            break
        parts = entry.split()
        if len(parts) != 2:
            print(color_text("Expected 2 values: SYMBOL COUNT.", 'red'))
            continue
        species.append(parts[0])
        num_ions.append(parts[1])

    if not species:
        print(color_text("No species given -- aborting.", 'red'))
        return

    groups = get_input("Space groups to try (comma- or space-separated, e.g. '225,227,229'): ").strip()
    while not groups:
        print(color_text("At least one space group is required.", 'red'))
        groups = get_input("Space groups to try: ").strip()

    count_per_group = get_int_input("Candidates per space group [default: 1]: ", 1)
    output_dir = get_input("Output folder [default: xrd_search]: ").strip() or "xrd_search"

    args = [
        "--species", *species,
        "--num-ions", *num_ions,
        "--groups", groups,
        "--count-per-group", str(count_per_group),
        "--output-dir", output_dir,
        "--no-intro",
    ]
    run_tool("stb-xrdsearch", args)

def run_xrdrank_analyzer() -> None:
    """Interface for the XRD-Guided Structure Search analysis stage (stb-xrdrank)"""
    print("\n" + "="*60)
    print(color_text("XRD STRUCTURE SEARCH -- RANK (stb-xrdrank)", 'bold').center(60))
    print("="*60 + "\n")

    input_dir = get_input("Folder of candidate structures (e.g. from stb-xrdsearch): ").strip()
    while not os.path.isdir(input_dir):
        print(color_text("Folder not found!", 'red'))
        input_dir = get_input("Folder of candidate structures: ").strip()

    experimental = get_input("Experimental XRD pattern file (2theta, intensity columns): ").strip()
    while not os.path.isfile(experimental):
        print(color_text("File not found!", 'red'))
        experimental = get_input("Experimental XRD pattern file: ").strip()

    wavelength = get_input("X-ray source name (CuKa, MoKa, ...) or wavelength in Ang [default: CuKa]: ").strip()
    if not wavelength:
        wavelength = "CuKa"

    top_str = get_input("Show only the N best matches [default: show all]: ").strip()
    output_file = get_input("Output ranking file name [default: xrd_rank.txt]: ").strip() or "xrd_rank.txt"

    args = ["--input-dir", input_dir, "--experimental", experimental,
            "--wavelength", wavelength, "--output", output_file, "--no-intro"]
    if top_str:
        args.extend(["--top", top_str])

    run_tool("stb-xrdrank", args)

def run_file_translator() -> None:
    """Interface for the File Translator (stb-translate)"""
    print("\n" + "="*60)
    print(color_text("FILE TRANSLATOR (stb-translate)", 'bold').center(60))
    print("="*60 + "\n")
    
    # Formatos suportados (o 'cif' foi adicionado à lista de saída)
    input_formats = ['fdf','poscar', 'cif', 'siesta', 'xyz', 'fhi', 'dftb', 'xsf']
    output_formats = ['cif', 'xyz', 'poscar', 'fdf', 'dftb', 'xsf', 'fhi'] # Adicionei 'cif' aqui
    
    input_file = get_input("Input file path: ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input file path: ")

    print(f"\n{color_text('Supported input formats:', 'yellow')}")
    for i, fmt in enumerate(input_formats, 1):
        print(f"  {color_text(str(i)+'.', 'yellow')} {fmt}")

    choice_in = 0
    max_in = len(input_formats)
    while not (1 <= choice_in <= max_in):
        choice_in = get_int_input(f"\nSelect input format (1-{max_in}): ")
        if not (1 <= choice_in <= max_in):
            print(color_text(f"Invalid choice! Please select between 1 and {max_in}.", 'red'))
    
    in_format = input_formats[choice_in - 1]
    print(f"Selected input format: {color_text(in_format, 'cyan')}")

    out_file = get_input("\nOutput file path: ")
    
    print(f"\n{color_text('Supported output formats:', 'yellow')}")
    for i, fmt in enumerate(output_formats, 1):
        print(f"  {color_text(str(i)+'.', 'yellow')} {fmt}")

    choice_out = 0
    max_out = len(output_formats)
    while not (1 <= choice_out <= max_out):
        choice_out = get_int_input(f"\nSelect output format (1-{max_out}): ")
        if not (1 <= choice_out <= max_out):
            print(color_text(f"Invalid choice! Please select between 1 and {max_out}.", 'red'))
            
    out_format = output_formats[choice_out - 1]
    print(f"Selected output format: {color_text(out_format, 'cyan')}")

    # ##### NOVO BLOCO: Seleção do Formato de Coordenadas #####
    print(f"\n{color_text('Select output coordinate format:', 'yellow')}")
    print(f"  {color_text('1.', 'yellow')} Cartesian (Angstroms)")
    print(f"  {color_text('2.', 'yellow')} Direct (Fractional)")
    print(f"  {color_text('3.', 'yellow')} Default (Use input format or output's default)")

    coord_choice = 0
    # Usamos default=3 para que pressionar Enter selecione a opção "Default"
    while not (1 <= coord_choice <= 3):
        coord_choice = get_int_input(f"\nSelect format (1-3) [default: 3]: ", 3) 
        if not (1 <= coord_choice <= 3):
            print(color_text(f"Invalid choice! Please select between 1 and 3.", 'red'))

    coord_format_value = None # Valor a ser passado para o argumento
    
    if coord_choice == 1:
        coord_format_value = "cartesian"
        print(f"Selected coordinate format: {color_text('Cartesian', 'cyan')}")
    elif coord_choice == 2:
        coord_format_value = "direct"
        print(f"Selected coordinate format: {color_text('Direct', 'cyan')}")
    else:
        # coord_format_value permanece None
        print(f"Selected coordinate format: {color_text('Default', 'cyan')}")
    # ##### FIM DO NOVO BLOCO #####

    # Construção dos argumentos base
    args = [
        "--in-format", in_format,
        "--in-file", input_file,
        "--out-format", out_format,
        "--out-file", out_file,
        "--no-intro"
    ]
    
    if coord_format_value:
        args.extend(["--coord-format", coord_format_value])


    if in_format == "xyz":
        print(color_text("\nXYZ format requires a separate lattice file.", 'yellow'))
        # Esta linha agora terá Tab-completion!
        lattice_file = get_input("Lattice vectors file (required for XYZ): ")
        while not os.path.isfile(lattice_file):
            print(color_text("File not found!", 'red'))
            lattice_file = get_input("Lattice vectors file: ")
        args.extend(["--lattice", lattice_file])
    
    run_tool("stb-translate", args)




def run_clean_tool() -> None:
    """Interactive interface for the Clean Files tool (stb-clean)"""
    print("\n" + "="*60)
    print(color_text("CLEAN FILES TOOL (stb-clean)", 'bold').center(60))
    print("="*60 + "\n")

    # Esta linha agora terá Tab-completion!
    path = get_input("Directory to clean (default: current): ").strip()
    if path == "":
        path = "."
    while not os.path.isdir(path):
        print(color_text("Directory not found!", 'red'))
        path = get_input("Enter a valid directory: ").strip()

    default_exts = ['.psml', '.psf', '.fdf', '.sh']
    print(f"\nExtensions to keep (space-separated, default: {' '.join(default_exts)}):")
    ext_input = get_input("Extensions: ").strip()
    if ext_input:
        extensions = ext_input.split()
    else:
        extensions = default_exts

    confirm_choice = get_input("Skip confirmation and delete directly? [y/N]: ").strip().lower()
    no_confirm = confirm_choice == 'y'

    dry_run_choice = get_input("Perform a dry run (show what would be deleted)? [y/N]: ").strip().lower()
    dry_run = dry_run_choice == 'y'

    args = ["--path", path, "--keep"] + extensions
    if no_confirm:
        args.append("--no-confirm")
    if dry_run:
        args.append("--dry-run")
    
    args.append("--no-intro")

    print()
    run_tool("stb-clean", args)

    if not dry_run:
        print("\n" + color_text("Cleanup complete. Your folder is now cleaner than my browser history.", "green"))

def run_symmetry_analyzer() -> None:
    """Interface for the Symmetry Analyzer (stb-symmetry)"""
    print("\n" + "="*60)
    print(color_text("SYMMETRY ANALYZER (stb-symmetry)", 'bold').center(60))
    print("="*60 + "\n")
    
    # Esta linha agora terá Tab-completion!
    input_file = get_input("Input structure file (SIESTA .fdf or .STRUCT_OUT): ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input structure file: ")

    formats = ['fdf', 'struct_out']
    print(f"\n{color_text('Select input file format:', 'yellow')}")
    for i, fmt in enumerate(formats, 1):
        print(f"  {color_text(str(i)+'.', 'yellow')} {fmt}")
    choice = 0
    while not (1 <= choice <= len(formats)):
        choice = get_int_input(f"\nSelect format (1-{len(formats)}): ")
        if not (1 <= choice <= len(formats)):
            print(color_text(f"Invalid choice! Please select between 1 and {len(formats)}.", 'red'))
    format_type = formats[choice - 1]
    print(f"Selected format: {color_text(format_type, 'cyan')}")

    output_dir = get_input("Output directory [default: current directory]: ").strip() or "."

    args = ["--file", input_file, "--format", format_type, "--output-dir", output_dir, "--no-intro"]

    scan_choice = get_input("Scan a range of symprec tolerances to check for hidden symmetry? (y/N): ").strip().lower()
    if scan_choice in ['y', 'yes']:
        args.append("--scan-symprec")

    ops_choice = get_input("Include the full symmetry-operations list in the report? (Y/n): ").strip().lower()
    if ops_choice in ['n', 'no']:
        args.append("--no-operations")

    compare_file = get_input("Compare to a second structure file (e.g. pre/post-relaxation)? "
                              "[optional, press Enter to skip]: ").strip()
    if compare_file:
        if not os.path.isfile(compare_file):
            print(color_text("File not found -- skipping comparison.", 'red'))
        else:
            print(f"\n{color_text('Select format of the comparison file:', 'yellow')}")
            for i, fmt in enumerate(formats, 1):
                print(f"  {color_text(str(i)+'.', 'yellow')} {fmt}")
            choice = 0
            while not (1 <= choice <= len(formats)):
                choice = get_int_input(f"\nSelect format (1-{len(formats)}): ")
                if not (1 <= choice <= len(formats)):
                    print(color_text(f"Invalid choice! Please select between 1 and {len(formats)}.", 'red'))
            args.extend(["--compare-to", compare_file, "--compare-format", formats[choice - 1]])

    write_refined_path = get_input("Write the symmetry-refined structure to a file? "
                                    "[optional, press Enter to skip]: ").strip()
    if write_refined_path:
        args.extend(["--write-refined", write_refined_path])

    run_tool("stb-symmetry", args)

def run_wantibexos_interface() -> None:
    """Interface for the Wantibexos (stb-siesta2wtb)"""
    print("\n" + "="*60)
    print(color_text("WANTIBEXOS INTERFACE (stb-siesta2wtb)", 'bold').center(60))
    print("="*60 + "\n")
    
    # Esta linha agora terá Tab-completion!
    input_file = get_input("Input FDF file: ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input FDF file: ")
    
    # Esta linha agora terá Tab-completion!
    output_file = get_input("SIESTA output file (optional): ")
    fermi = get_input("Manual Fermi level (optional): ")
    
    args = ["--input", input_file]
    if output_file:
        # Validação extra para o ficheiro opcional
        if os.path.isfile(output_file):
            args.extend(["--output", output_file])
        else:
            print(color_text(f"Warning: Output file '{output_file}' not found, skipping.", 'yellow'))
    if fermi:
        args.extend(["--fermi-level", fermi])
    
    args.append("--no-intro")
    
    run_tool("stb-siesta2wtb", args)

# ==========================================================
# SUB-MENU LOGIC
# ==========================================================

# Define the tool dictionaries
INPUT_TOOLS = {
    1: {'title': "Input File Generator (stb-inputfile)",
        'description': "Create a 'calc.fdf' input file from a structure file.",
        'func': run_input_generator},
    2: {'title': "K-Grid Generator (stb-kgrid)",
        'description': "Suggest a Monkhorst-Pack grid (k-points) based on desired density.",
        'func': run_kgrid_generator},
    3: {'title': "K-Path Generator (stb-kpath)",
        'description': "Generate a high-symmetry k-path for band structure calculations.",
        'func': run_kpath_generator},
    4: {'title': "DFT+U / Hubbard Block Generator (stb-dftu)",
        'description': "Generate a ready-to-use %block LDAU.proj snippet from a user-supplied U (and J).",
        'func': run_dftu_generator},
    5: {'title': "Structure Fetcher (stb-fetch)",
        'description': "Fetch a structure from Materials Project or COD and write it as .fdf.",
        'func': run_fetch_generator},
    6: {'title': "ML Pre-Relaxation (stb-mlrelax)",
        'description': "Fast pre-relaxation with the MACE-MP-0 potential (needs the optional 'ml' extra).",
        'func': run_mlrelax_generator},
       }


# Tools that build, generate, or transform a structure file, independent of any
# SIESTA-specific input setup -- split out from INPUT_TOOLS once it grew to 17
# entries mixing the two concerns.
STRUCTURE_TOOLS = {
    1: {'title': "2D Monolayer Stacker (stb-2Dstacking)",
        'description': "Stacks two monolayers into a heterostructure using the ZSL algorithm.",
        'func': run_2d_stacker},
    2: {'title': "Supercell Builder (stb-supercell)",
        'description': "Build a supercell from a structure file using a user-defined transformation matrix.",
        'func': run_supercell_generator},
    3: {'title': "Slab Builder (stb-slab)",
        'description': "Cut a Miller-index slab from a bulk structure, with vacuum.",
        'func': run_slab_generator},
    4: {'title': "Nanotube/Nanoribbon Builder (stb-nanotube)",
        'description': "Roll a 2D monolayer into a nanotube or nanoribbon, given (n, m) chirality.",
        'func': run_nanotube_generator},
    5: {'title': "Point Defect Generator (stb-defect)",
        'description': "Introduce a vacancy, substitution, or interstitial defect.",
        'func': run_defect_generator},
    6: {'title': "SQS Generator (stb-sqs)",
        'description': "Generate a Special Quasirandom Structure for a substitutional alloy.",
        'func': run_sqs_generator},
    7: {'title': "Unit Cell Finder (stb-unitcell)",
        'description': "Find the primitive or conventional unit cell of a structure.",
        'func': run_unitcell_generator},
    8: {'title': "Crystal Builder (stb-crystalbuilder)",
        'description': "Build a structure from a space group and Wyckoff positions.",
        'func': run_crystalbuilder_generator},
    9: {'title': "Surface Passivator (stb-passivate)",
        'description': "Cap dangling bonds on a cut surface with a passivating atom (e.g. H).",
        'func': run_passivate_generator},
    10: {'title': "Reference Molecule Builder (stb-molecule)",
         'description': "Build a reference molecule (e.g. H2O, CO2) from ASE's G2 database.",
         'func': run_molecule_generator},
    11: {'title': "Amorphous Structure Generator (stb-amorphize)",
         'description': "Melt-quench MD with MACE-MP-0 to build an amorphous starting guess (needs the optional 'ml' extra).",
         'func': run_amorphize_generator},
    12: {'title': "Random Crystal Generator (stb-crystalcast)",
         'description': "Cast random bulk/layer/rod/cluster structures (atomic or molecular, "
                         "optionally ML-ranked) from a symmetry group and composition; or "
                         "analyze/substitute/subgroup/supergroup an existing structure.",
         'func': run_crystalcast_generator},
        }


ANALYSIS_TOOLS = {
    1: {'title': "Bands Analyzer (stb-bands)",
        'description': "Analyze .bands files and calculate band gaps.",
        'func': run_bands_analyzer},
    2: {'title': "PDOS XML Parser (stb-dos)",
        'description': "Extract data from PDOS.xml by total, atom, and species.",
        'func': run_dos_parser},
    3: {'title': "DOS Processor (Convolution) (stb-convdos)",
        'description': "Apply Gaussian convolution to Density of States (DOS) files.",
        'func': run_dos_convolution},
    4: {'title': "Structure Analyzer (stb-structural)",
        'description': "Calculate ECN and analyze structural properties.",
        'func': run_structure_analyzer},
    5: {'title': "Symmetry Analyzer (stb-symmetry)",
        'description': "Analyze the symmetry of crystal structures.",
        'func': run_symmetry_analyzer},
    6: {'title': 'Bader Charge Analysis',
        'description': 'Calculate atomic charges using the Bader AIM method from .RHO and .XV files.',
        'func': run_bader_calculator},
    7: {'title': "Work Function Calculator", 'description': "Calculate Work Function from electrostatic potential (.VT).",
        'func': run_workfunction_calculator},
    8: {'title': "Density Plotter (RHO)", 'description': "Export 2D Charge Density Maps or 3D Clouds.",
        'func': run_density_plotter},
    9: {'title': "Powder XRD Simulator (stb-xrd)",
        'description': "Simulate a powder XRD pattern (peak table + optional plot) from a "
                        "structure, optionally compared against an experimental pattern.",
        'func': run_xrd_analyzer},
       }


WORKFLOW_TOOLS = {
    1: {'title': "Stress-Strain",
        'description': "Generate strained structures, then extract stress-strain curves and mechanical properties.",
        'stages': {
            1: {'title': "Stage 1 - Prep (stb-strain)",
                'description': "Generate strained structures for calculations.",
                'func': run_strain_generator},
            2: {'title': "Stage 2 - Analysis (stb-strainAnalysis)",
                'description': "Extract stress-strain curves from strain_* folders.",
                'func': run_strain_post_processor},
        }},
    2: {'title': "Elastic Constants",
        'description': "Generate deformed structures, then compute the stiffness matrix and elastic moduli.",
        'stages': {
            1: {'title': "Stage 1 - Prep (stb-elasticInputs)",
                'description': 'Generates deformed structures to calculate elastic constants.',
                'func': run_elastic_generator},
            2: {'title': "Stage 2 - Analysis (stb-elasticAnalysis)",
                'description': 'Calculates Stiffness Matrix, Young Modulus and Stability from outputs.',
                'func': run_elastic_analyzer},
        }},
    3: {'title': "Cohesive Energy",
        'description': "Set up bulk + isolated-atom calculations, then compute the cohesive energy per atom.",
        'stages': {
            1: {'title': "Stage 1 - Prep (stb-cohesive)",
                'description': "Prepare folder structure and inputs for cohesive energy calculations.",
                'func': run_cohesive_setup},
            2: {'title': "Stage 2 - Analysis (stb_cohesive_analysis)",
                'description': "Process and calculate the final cohesive energy per atom.",
                'func': run_cohesive_analysis},
        }},
    4: {'title': "Phonons",
        'description': "Generate displaced supercells via Phonopy, then post-process into thermal properties.",
        'stages': {
            1: {'title': "Stage 1 - Prep (stb-phononsCreate)",
                'description': "Automate SIESTA phonon displacement folders using Phonopy.",
                'func': run_phonon_generator},
            2: {'title': "Stage 2 - Analysis (stb-phononsPos)",
                'description': "Extract forces, generate FORCE_SETS, and calculate thermal properties.",
                'func': run_phonon_postprocessing},
        }},
    5: {'title': "Convergence Tests",
        'description': "Sweep Mesh.CutOff, k-grid density, or PAO.EnergyShift and check total-energy convergence.",
        'stages': {
            1: {'title': "Stage 1 - Prep (stb-convergence)",
                'description': "Generate a sweep of calc.fdf variants for one parameter.",
                'func': run_convergence_generator},
            2: {'title': "Stage 2 - Analysis (stb-convergenceAnalysis)",
                'description': "Extract energies from convergence_* folders and report the converged value.",
                'func': run_convergence_analyzer},
        }},
    6: {'title': "Structure Solution (XRD)",
        'description': "Cast candidate structures across a set of space groups, then rank them by "
                        "similarity to an experimental powder XRD pattern.",
        'stages': {
            1: {'title': "Stage 1 - Prep (stb-xrdsearch)",
                'description': "Cast candidate structures across a set of space groups for a given composition.",
                'func': run_xrdsearch_generator},
            2: {'title': "Stage 2 - Analysis (stb-xrdrank)",
                'description': "Rank candidate structures by similarity to an experimental XRD pattern.",
                'func': run_xrdrank_analyzer},
        }},
    7: {'title': "Hubbard U (Linear Response)",
        'description': "Compute a first-principles Hubbard U via Cococcioni & de Gironcoli's "
                        "linear-response method, ending in a ready-to-use DFT+U fdf block.",
        'stages': {
            1: {'title': "Stage 1 - Reference Prep (stb-hubbardu)",
                'description': "Generate the 'reference' folder -- run SIESTA in it before stage 2.",
                'func': run_hubbardu_prep},
            2: {'title': "Stage 2 - Perturbation Prep (stb-hubbarduAlphas)",
                'description': "Generate the scf/frozen run folders, auto-copying the reference DM.",
                'func': run_hubbardu_alphas},
            3: {'title': "Stage 3 - Analysis (stb-hubbarduAnalysis)",
                'description': "Fit the occupation responses, compute U, and write the DFT+U block.",
                'func': run_hubbardu_analysis},
        }},
       }


UTILITY_TOOLS = {
    1: {'title': "File Translator (stb-translate)",
        'description': "Convert between file formats (CIF, POSCAR, fdf, xyz...).",
        'func': run_file_translator},
    2: {'title': "Clean File Tools (stb-clean)",
        'description': "Clean the directory of calculation files (except essential ones).",
        'func': run_clean_tool},
    3: {'title': "Grid to Cube Converter", 'description': "Convert Siesta Grid files (VT, RHO, VH) to Gaussian .cube.", 'func': run_grid_to_cube},
    4: {'title': "Wantibexos Interface (stb-siesta2wtb)",
        'description': "Convert SIESTA Hamiltonian to Wantibexos format.",
        'func': run_wantibexos_interface},
}


def _flatten_tool_codes() -> Dict[str, Callable]:
    """Builds a flat {"1.1": func, ..., "2.1": func, ..., "4.1.2": func, ...,
    "5.4": func} lookup across all 5 categories, from the dicts above -- lets
    the main menu jump straight to a tool via a dotted code instead of
    navigating level by level.
    """
    codes: Dict[str, Callable] = {}
    for key, info in INPUT_TOOLS.items():
        codes[f"1.{key}"] = info['func']
    for key, info in STRUCTURE_TOOLS.items():
        codes[f"2.{key}"] = info['func']
    for key, info in ANALYSIS_TOOLS.items():
        codes[f"3.{key}"] = info['func']
    for prop_key, prop_info in WORKFLOW_TOOLS.items():
        for stage_key, stage_info in prop_info['stages'].items():
            codes[f"4.{prop_key}.{stage_key}"] = stage_info['func']
    for key, info in UTILITY_TOOLS.items():
        codes[f"5.{key}"] = info['func']
    return codes


TOOL_CODES = _flatten_tool_codes()


def run_sub_menu(title: str, tools_dict: Dict) -> None:
    """Handles the logic for showing and running a sub-menu.

    Recurses one level deeper when an entry has a 'stages' dict instead of a
    'func' (this is how Workflow's category -> property -> stage nesting
    works, with no separate menu function needed).
    """
    while True:
        show_sub_menu(title, tools_dict)
        try:
            choice_str = get_input(f"\nSelect an option (0-{len(tools_dict)}): ")

            if choice_str == '0':
                break # Go back to the main menu

            try:
                choice = int(choice_str)
            except ValueError:
                choice = float(choice_str)

            if choice in tools_dict:
                entry = tools_dict[choice]
                if 'stages' in entry:
                    run_sub_menu(entry['title'], entry['stages'])
                else:
                    entry['func']() # Run the selected tool
            else:
                print(color_text(f"\nInvalid choice! Please select between 0 and {len(tools_dict)}.", 'red'))
                sleep(1)
                
        except ValueError:
            print(color_text("\nPlease enter a valid number!", 'red'))
            sleep(1)
        except KeyboardInterrupt:
            break # Go back to the main menu

# ==========================================================
# MAIN FUNCTION
# ==========================================================

def main():
    """Main function to run the STB-SUITE interface"""
    show_intro([
        "Siesta ToolBox Suite",
        "A comprehensive toolkit for SIESTA DFT simulations",
        f"Version {VERSION} | University of Brasilia - 2025",
        "Developed by Dr. Carlos M. O. Bastos"
    ])

    while True:
        show_main_menu()
        
        try:
            choice = get_input("\nSelect an option (0-5, or a tool code like 4.1.2): ")

            if choice in TOOL_CODES:
                TOOL_CODES[choice]()
            elif choice == '1':
                run_sub_menu("Inputs", INPUT_TOOLS)
            elif choice == '2':
                run_sub_menu("Structures", STRUCTURE_TOOLS)
            elif choice == '3':
                run_sub_menu("Analysis", ANALYSIS_TOOLS)
            elif choice == '4':
                run_sub_menu("Workflow", WORKFLOW_TOOLS)
            elif choice == '5':
                run_sub_menu("Utils", UTILITY_TOOLS)
            elif choice == '0':
                print(color_text("\nThank you for using STB-SUITE!", 'cyan'))
                break
            else:
                print(color_text("\nInvalid choice! Please select between 0 and 5, or a valid tool code.", 'red'))
                sleep(1)
                
        except ValueError:
            print(color_text("\nPlease enter a valid number!", 'red'))
            sleep(1)
        except KeyboardInterrupt:
            print(color_text("\n\nOperation cancelled by user.", 'yellow'))
            break

if __name__ == "__main__":
    main()
