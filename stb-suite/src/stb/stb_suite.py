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

import numpy as np

try:
    import readline
    readline.parse_and_bind("tab: complete")
except ImportError:
    pass

from stb.core.cli import COLORS, color_text, show_intro, get_input, get_float_input, get_int_input
from stb.core.pseudopotentials import BANKS
from stb.core import structure_io, kspace, symmetry
from stb import strain

def prompt_pseudo_source(optional: bool = True) -> str:
    """Shared pseudopotential-source prompt for every wrapper below that
    needs one (phonons, cohesive energy, input file, Hubbard U prep):
    a bundled bank (see core/pseudopotentials.py) or a custom path. Returns
    the raw string to pass straight through as the tool's -p/--pp-path/
    --pseudo-dir value (each tool resolves it itself); empty string only if
    `optional` and the user skips."""
    bank_list = list(BANKS.items())
    print(f"\n{color_text('Pseudopotential source:', 'yellow')}")
    for i, (name, info) in enumerate(bank_list, 1):
        print(f"  {color_text(str(i), 'cyan')} = Bundled: {info['description']} ({name})")
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
    print(f"{color_text('5.', 'yellow')} {color_text('ML Simulations', 'blue')}\n    Run MD (and more, over time) with a trained MACE potential -- the foundation model or your own fine-tuned one\n")
    print(f"{color_text('6.', 'yellow')} {color_text('Utils', 'blue')}\n    Helper tools for file management and conversion\n")
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

def run_tool(tool_name: str, args: List[str], pause: bool = True) -> None:
    """Executes a suite tool as a subprocess.

    `pause=False` skips the "Press Enter to continue..." block -- for callers
    that invoke run_tool() several times in a row (e.g. run_strain_generator
    looping over multiple symmetry-equivalent directions in one go), so the
    user isn't interrupted after every individual subprocess call, only once
    at the very end.
    """
    try:
        cmd = [f"{tool_name}"] + args
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(color_text(f"\nError running {tool_name}: {e}", 'red'))
    except FileNotFoundError:
        print(color_text(f"\nTool {tool_name} not found!", 'red'))
        print(color_text(f"Make sure {tool_name} is in your system's PATH.", 'yellow'))
    if pause:
        input("\nPress Enter to continue...")

# ==========================================================
# TOOL FUNCTIONS
# ==========================================================


def run_phonon_postprocessing() -> None:
    """Interface for the Phonon Post-Processing (phonons_pos.py)"""
    print("\n" + "="*60)
    print(color_text("PHONON POST-PROCESSING", 'bold').center(60))
    print("="*60 + "\n")

    # 1. Diretório
    phonon_dir = get_input("Phonon runs directory [default: phonon_runs]: ").strip()
    if not phonon_dir:
        phonon_dir = "phonon_runs"

    # 2. System Label -- left blank means "auto-detect from calc.fdf" (the
    # CLI's own default); NOT the same as forcing "siesta" -- an ML-sourced
    # directory (e.g. one produced by stb-mlphonons, see ML Simulations) has
    # no SystemLabel at all and must never get -l passed.
    sys_label = get_input(
        "SystemLabel used in calculations [blank = auto-detect from calc.fdf, "
        "or N/A for an ML-sourced directory]: ").strip()

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

    # 5. Análises extra (Tiers 2-3) -- lista numerada em vez de exigir digitar
    # as palavras-chave de cabeça (um "y" de hábito antigo, ou um typo, batia
    # em nenhuma keyword e era descartado em silêncio -- nenhum erro, nenhum
    # aviso, --bands nunca chegava a ser passado. Reproduzido manualmente:
    # digitar "y" aqui rodava o resto normalmente e só faltava band structure
    # no resumo final, sem nada que apontasse a causa).
    EXTRA_ANALYSIS_OPTIONS = [
        ("bands", "Band structure (--bands)"),
        ("dos", "Total DOS (--dos)"),
        ("pdos", "Projected DOS (--pdos)"),
        ("thermal", "Thermal displacements (--thermal-displacements)"),
        ("freeze", "Freeze unstable mode, if one was found (--freeze-unstable-mode)"),
    ]
    print(f"\n{color_text('Additional analyses (Tiers 2-3) -- select any that apply:', 'yellow')}")
    for i, (_key, desc) in enumerate(EXTRA_ANALYSIS_OPTIONS, start=1):
        print(f"  {color_text(str(i), 'cyan')} = {desc}")
    ALL_ANALYSES_KEYS = {"bands", "dos", "pdos", "thermal"}  # freeze (5) excluded --
    # a different kind of action (writes a structure file), not just a plot/report.
    raw_choice = get_input(
        f"Select by number (1-{len(EXTRA_ANALYSIS_OPTIONS)}), space/comma-separated, "
        "'all' for 1-4, or blank for none: ").strip().lower()

    selected = set()
    unrecognized = []
    option_keys = {key for key, _ in EXTRA_ANALYSIS_OPTIONS}
    for tok in raw_choice.replace(',', ' ').split():
        if tok == 'all':
            selected |= ALL_ANALYSES_KEYS
        elif tok.isdigit() and 1 <= int(tok) <= len(EXTRA_ANALYSIS_OPTIONS):
            selected.add(EXTRA_ANALYSIS_OPTIONS[int(tok) - 1][0])
        elif tok in option_keys:  # old-style keyword still accepted
            selected.add(tok)
        else:
            unrecognized.append(tok)

    if unrecognized:
        print(color_text(
            f"[WARNING] Ignored unrecognized option(s): {', '.join(unrecognized)} "
            f"(expected a number 1-{len(EXTRA_ANALYSIS_OPTIONS)}).", 'yellow'))
    if selected:
        print(color_text(f"Selected: {', '.join(sorted(selected))}", 'green'))

    want_bands = 'bands' in selected
    want_dos = 'dos' in selected
    want_pdos = 'pdos' in selected
    want_thermal = 'thermal' in selected
    want_freeze = 'freeze' in selected

    # Advanced settings (rarely-touched -- gated so the essential flow above
    # stays short; CLI defaults apply untouched when skipped).
    plot_dir, band_points, vacuum_gap, freeze_amplitude, symprec = (
        "phonon_plots", 101, 10.0, 0.05, 0.01)
    advanced_items = "plot directory, symmetry tolerance"
    if want_bands:
        advanced_items += ", band-path resolution, vacuum-gap"
    if want_freeze:
        advanced_items += ", freeze amplitude"
    show_advanced = get_input(f"\nConfigure advanced settings ({advanced_items})? [y/N]: ").strip().lower()
    if show_advanced == 'y':
        plot_dir = get_input("Gnuplot .dat/.gplot output directory (default: phonon_plots): ").strip()
        if not plot_dir:
            plot_dir = "phonon_plots"
        symprec = get_float_input(
            "Symmetry-detection tolerance Ang, used for force-constant symmetrization "
            "and the band-path (default: 0.01 -- NOT Phonopy's/ASE's own much tighter "
            "raw defaults; loosen, don't tighten, if unsure -- a much tighter value "
            "than Stage 1 used can fail to reconstruct force constants): ", 0.01)
        if want_bands:
            band_points = get_int_input(
                "Q-points per band-structure path segment (default: 101): ", 101)
            vacuum_gap = get_float_input(
                "Vacuum gap threshold in Ang, for --bands' periodicity detection "
                "(default: 10.0): ", 10.0)
        if want_freeze:
            freeze_amplitude = get_float_input(
                "Target max atomic displacement for the frozen mode, Ang "
                "(default: 0.05): ", 0.05)

    save_report = get_input("\nAlso save a text report to file? (y/N): ").strip().lower() == 'y'
    save_gnuplot = get_input(
        "Also save gnuplot .dat + .gplot scripts for every computed quantity? (y/N): "
    ).strip().lower() == 'y'
    view_plots = get_input(
        "View plots interactively via matplotlib now? (y/N): ").strip().lower() == 'y'

    args = [
        "-dir", phonon_dir,
        "-m", str(m_x), str(m_y), str(m_z),
        "--tmin", str(tmin),
        "--tmax", str(tmax),
        "--tstep", str(tstep),
        "--plot-dir", plot_dir,
        "--symprec", str(symprec),
        "--no-intro"
    ]
    if sys_label and sys_label.upper() != "N/A":
        args += ["-l", sys_label]
    if want_bands:
        args += ["--bands", "--band-points", str(band_points), "--vacuum-gap", str(vacuum_gap)]
    if want_dos:
        args.append("--dos")
    if want_pdos:
        args.append("--pdos")
    if want_thermal:
        args.append("--thermal-displacements")
    if want_freeze:
        args += ["--freeze-unstable-mode", "--freeze-amplitude", str(freeze_amplitude)]
    if save_report:
        args.append("--save-report")
    if save_gnuplot:
        args.append("--save-gnuplot")
    if view_plots:
        args.append("--view")

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

    # 3b. Malha de k-points da supercélula -- preview the 0.2 1/Ang-density
    # suggestion as the bracketed default, so the user sees real numbers
    # instead of a blind "blank to auto-suggest".
    suggested_kgrid = None
    try:
        _struct = structure_io.read_fdf(structure_file)
        _pos = np.array([p for _, p in _struct.atoms])
        _frac = kspace.to_fractional(_pos, _struct.lattice, _struct.coord_format == 'cartesian')
        _vac = kspace.detect_vacuum_axes(_frac, _struct.lattice, 10.0)
        _sc_lattice = np.diag([dim_x, dim_y, dim_z]) @ np.array(_struct.lattice)
        suggested_kgrid = kspace.compute_monkhorts(
            _sc_lattice[0], _sc_lattice[1], _sc_lattice[2], 0.2, _vac)
    except Exception:
        suggested_kgrid = None

    if suggested_kgrid:
        default_str = f"{suggested_kgrid[0]} {suggested_kgrid[1]} {suggested_kgrid[2]}"
        kgrid_choice = get_input(f"\nSupercell k-grid [{default_str}]: ").strip()
        if not kgrid_choice:
            kgrid_choice = default_str
    else:
        kgrid_choice = get_input("\nSupercell k-grid 'X Y Z' (blank to auto-suggest): ").strip()
    kgrid = None
    if kgrid_choice:
        try:
            kgrid_vals = [int(x) for x in kgrid_choice.split()]
            if len(kgrid_vals) == 3:
                kgrid = kgrid_vals
            else:
                print(color_text("Please provide exactly 3 integers. Auto-suggesting instead.", 'yellow'))
        except ValueError:
            print(color_text("Invalid input format. Auto-suggesting instead.", 'yellow'))

    # 4. Definir distância de deslocamento
    distance = get_float_input("\nDisplacement distance in Å [default: 0.01]: ", 0.01)

    # 5. Diretório dos pseudopotenciais
    pseudo_dir = prompt_pseudo_source(optional=True)
    if not pseudo_dir:
        pseudo_dir = "."

    # 6. Checagem MACE opcional (Tier "ML") -- pergunta direto, não atrás do
    # gate avançado, já que é um recurso com valor próprio, não uma
    # tolerância rara de mexer.
    ml_choice = get_input(
        "\nRun a MACE-MP-0 pre-flight check before generating displacements? "
        "Catches an unrelaxed reference structure (a common cause of spurious "
        "imaginary modes) before spending SIESTA time on it (y/N): ").strip().lower()
    ml_prerelax = ml_choice in ('y', 'yes')
    ml_model, ml_fmax = "small", 0.05
    if ml_prerelax:
        ml_model = get_input(
            "  MACE-MP-0 model size [small/medium/large, default: small]: ").strip().lower()
        if ml_model not in ("small", "medium", "large"):
            ml_model = "small"
        ml_fmax = get_float_input(
            "  Force threshold above which a relax is offered, eV/Ang [default: 0.05]: ", 0.05)

    # Advanced settings (rarely-touched -- gated so the essential flow above
    # stays short; CLI defaults apply untouched when skipped).
    vacuum_gap, symprec, ml_device = 10.0, 0.01, "cpu"
    advanced_items = "vacuum-gap, symmetry tolerance"
    if ml_prerelax:
        advanced_items += ", ML device"
    show_advanced = get_input(f"\nConfigure advanced settings ({advanced_items})? [y/N]: ").strip().lower()
    if show_advanced == 'y':
        vacuum_gap = get_float_input(
            "Vacuum gap threshold in Ang, for the supercell-dimension advisory "
            "(default: 10.0): ", 10.0)
        symprec = get_float_input(
            "Symmetry-detection tolerance Phonopy uses to reduce displacements, in Ang "
            "(default: 0.01 -- NOT Phonopy's own much tighter raw default, which can "
            "misdetect the true space group on a real relaxed structure): ", 0.01)
        if ml_prerelax:
            device_choice = get_input("ML device [cpu/cuda, default: cpu]: ").strip().lower()
            ml_device = device_choice if device_choice in ("cpu", "cuda") else "cpu"

    save_report = get_input("\nAlso save a text report to file? (y/N): ").strip().lower() == 'y'

    # 7. Preparar e executar o script
    args = [
        "-s", structure_file,
        "-c", calc_file,
        "-dim", str(dim_x), str(dim_y), str(dim_z),
        "-d", str(distance),
        "-p", pseudo_dir,
        "--vacuum-gap", str(vacuum_gap),
        "--symprec", str(symprec),
        "--no-intro"
    ]
    if kgrid is not None:
        args += ["--kgrid", str(kgrid[0]), str(kgrid[1]), str(kgrid[2])]
    if ml_prerelax:
        args += ["--ml-prerelax", "--ml-model", ml_model, "--ml-device", ml_device,
                  "--ml-fmax", str(ml_fmax)]
    if save_report:
        args.append("--save-report")

    print(color_text("\nGenerating phonon displacement folders...", 'green'))
    run_tool("stb-phononsCreate", args)


def run_cohesive_setup() -> None:
    """Interface for the Cohesive Energy Setup (cohesive_energy.py)"""
    print("\n" + "="*60)
    print(color_text("COHESIVE ENERGY SETUP", 'bold').center(60))
    print("="*60)
    print(color_text(
        "Sets up a full-structure calculation plus one isolated-atom reference per species. "
        "By default it also builds a BSSE (counterpoise) -corrected reference per species -- "
        "LCAO cohesive energies otherwise systematically over-bind, since the isolated atom "
        "is computed with a smaller effective basis than it has inside the solid.", 'cyan'))
    print()

    struct_file = get_input("Input structure FDF file (-s): ").strip()
    while not os.path.isfile(struct_file):
        print(color_text("File not found!", 'red'))
        struct_file = get_input("Input structure FDF file (-s): ").strip()

    k_density = get_float_input("K-point density (default: 0.2): ", 0.2)
    while k_density <= 0:
        print(color_text("Density must be a positive number!", 'red'))
        k_density = get_float_input("K-point density (default: 0.2): ", 0.2)

    pp_path = prompt_pseudo_source(optional=True)

    spin_choice = get_input("Enable spin polarization for full structure? (y/N): ").strip().lower()

    dispersion_choice = get_input(
        "\nEnable DFT-D3 dispersion correction for every calculation? (y/N -- relevant for "
        "layered/molecular materials, e.g. 2D materials, where van der Waals interactions "
        "contribute non-negligibly): ").strip().lower()
    dispersion = dispersion_choice in ('y', 'yes')

    vacuum = get_float_input("\nIsolated-atom vacuum box side in Ang (default: 20.0): ", 20.0)
    while vacuum <= 0:
        print(color_text("Vacuum box side must be a positive number!", 'red'))
        vacuum = get_float_input("Isolated-atom vacuum box side in Ang (default: 20.0): ", 20.0)

    # BSSE correction (default ON, matching the CLI default)
    bsse_choice = get_input(
        "\nGenerate BSSE (counterpoise) -corrected references too? (Y/n): ").strip().lower()
    bsse_correction = bsse_choice not in ('n', 'no')

    bsse_multi_site = True
    bsse_convergence_check = False
    if bsse_correction:
        multi_site_choice = get_input(
            "  If a species occupies multiple distinct crystallographic sites (e.g. "
            "octahedral vs. tetrahedral), build one BSSE reference per site, weighted by "
            "multiplicity? (Y/n): ").strip().lower()
        bsse_multi_site = multi_site_choice not in ('n', 'no')

        conv_check_choice = get_input(
            "  Also run a BSSE cutoff convergence check (2nd, larger cutoff -- roughly "
            "triples the isolated-atom-type calculation count)? (y/N): ").strip().lower()
        bsse_convergence_check = conv_check_choice in ('y', 'yes')

    # Advanced settings (rarely-touched tolerances -- gated so the essential
    # flow above stays short; CLI defaults apply untouched when skipped).
    vacuum_gap, output_dir, bsse_cutoff, symprec, bsse_convergence_increment = 10.0, "cohesive_runs", 4.0, 0.01, [2.0]
    advanced_items = "vacuum-gap, output directory"
    if bsse_correction:
        advanced_items += ", BSSE cutoff radius/symmetry tolerance"
    if bsse_convergence_check:
        advanced_items += ", convergence-check increment"
    show_advanced = get_input(f"\nConfigure advanced settings ({advanced_items})? [y/N]: ").strip().lower()
    if show_advanced == 'y':
        vacuum_gap = get_float_input(
            "Vacuum gap threshold in Ang, for detecting periodic vs. vacuum-padded axes "
            "(default: 10.0): ", 10.0)
        output_dir = get_input("Output root directory (default: cohesive_runs): ").strip()
        if not output_dir:
            output_dir = "cohesive_runs"
        if bsse_correction:
            bsse_cutoff = get_float_input(
                "BSSE ghost-neighbor cutoff radius in Ang (default: 4.0 -- reaches 1 shell "
                "for longer-bonded solids, several for short-bonded covalent ones like "
                "graphene/diamond; check the printed ghost-neighbor count): ", 4.0)
            symprec = get_float_input(
                "Symmetry-detection tolerance, for multi-site BSSE (default: 0.01): ", 0.01)
            if bsse_convergence_check:
                inc_input = get_input(
                    "BSSE convergence-check cutoff increment(s) in Ang -- one value for a "
                    "single before/after check, or several space-separated values (e.g. "
                    "'2 4 6') for a full convergence scan (default: 2.0): ").strip()
                if inc_input:
                    try:
                        bsse_convergence_increment = [float(x) for x in inc_input.split()]
                    except ValueError:
                        print(color_text("Invalid input, using default (2.0).", 'red'))
                        bsse_convergence_increment = [2.0]

    save_report = get_input("\nAlso save a report to file? [y/N]: ").strip().lower() == 'y'

    args = [
        "-s", struct_file,
        "-k", str(k_density),
        "--vacuum", str(vacuum),
        "--vacuum-gap", str(vacuum_gap),
        "--output-dir", output_dir,
        "--no-intro",
    ]
    if dispersion:
        args.append("--dispersion")
    args.append("--bsse-correction" if bsse_correction else "--no-bsse-correction")
    if bsse_correction:
        args.extend(["--bsse-cutoff", str(bsse_cutoff), "--symprec", str(symprec)])
        args.append("--bsse-multi-site" if bsse_multi_site else "--no-bsse-multi-site")
        if bsse_convergence_check:
            args.extend(["--bsse-convergence-check", "--bsse-convergence-increment"])
            args.extend(str(v) for v in bsse_convergence_increment)

    if pp_path:
        args.extend(["-p", pp_path])
    if spin_choice in ['y', 'yes']:
        args.append("--spin")
    if save_report:
        args.append("--save-report")

    # Recap -- always shown, no confirmation gate
    bsse_desc = "OFF"
    if bsse_correction:
        bsse_desc = (f"ON (cutoff={bsse_cutoff} Ang, "
                     f"multi-site={'on' if bsse_multi_site else 'off'})")
    summary_rows = [
        ("Structure file", struct_file),
        ("K-point density", f"{k_density} 1/Ang"),
        ("Pseudopotentials", pp_path or "(none)"),
        ("Spin (full structure)", "polarized" if spin_choice in ['y', 'yes'] else "non-polarized"),
        ("Dispersion (D3)", "ON" if dispersion else "OFF"),
        ("Isolated-atom vacuum box", f"{vacuum} Ang"),
        ("BSSE correction", bsse_desc),
    ]
    if bsse_correction and bsse_convergence_check:
        incs = ", ".join(f"+{v:g}" for v in bsse_convergence_increment)
        summary_rows.append(("BSSE convergence check", f"ON ({incs} Ang)"))
    summary_rows.append(("Vacuum-gap threshold", f"{vacuum_gap} Ang"))
    summary_rows.append(("Output directory", output_dir))
    summary_rows.append(("Save report", "yes" if save_report else "no"))
    _print_config_summary("CONFIGURATION SUMMARY", summary_rows)

    run_tool("stb-cohesive", args)


def run_cohesive_analysis() -> None:
    """Interface for the Cohesive Energy Analysis (cohesive_analysis.py)"""
    print("\n" + "="*60)
    print(color_text("COHESIVE ENERGY ANALYSIS", 'bold').center(60))
    print("="*60)
    print(color_text(
        "Extracts energies from 'structure'/'atoms' (and, if present, 'atoms_bsse') and "
        "reports the cohesive energy -- auto-detects a BSSE-corrected reference, no flag "
        "needed.", 'cyan'))
    print()

    out_file = get_input("SIESTA output file name (e.g., calc.out) [-o]: ").strip()
    while not out_file:
        print(color_text("File name cannot be empty!", 'red'))
        out_file = get_input("SIESTA output file name [-o]: ").strip()

    dir_path = get_input(
        "Path to results folder containing 'structure' and 'atoms' (default: cohesive_runs) [-d]: "
    ).strip()

    force_tolerance = get_float_input(
        "Force tolerance for the 'is the full structure relaxed' check, in eV/Ang "
        "(default: 0.05): ", 0.05)

    # Save report / view plot?
    save_report = get_input("\nAlso save a text report to file? (y/N): ").strip().lower() == 'y'
    view = get_input(
        "Show an interactive matplotlib plot before finishing? (y/N): ").strip().lower() == 'y'

    args = [
        "-o", out_file,
        "--force-tolerance", str(force_tolerance),
        "--no-intro",
    ]
    if dir_path:
        args.extend(["-d", dir_path])
    if save_report:
        args.append("--save-report")
    if view:
        args.append("--view")

    # Recap -- always shown, no confirmation gate
    summary_rows = [
        ("SIESTA output filename", out_file),
        ("Results directory", dir_path or "cohesive_runs (default)"),
        ("Force tolerance", f"{force_tolerance} eV/Ang"),
        ("Save report", "yes" if save_report else "no"),
        ("View (matplotlib)", "yes" if view else "no"),
    ]
    _print_config_summary("CONFIGURATION SUMMARY", summary_rows)

    run_tool("stb-cohesiveAnalysis", args)


def run_adsorb_setup() -> None:
    """Interface for the Adsorption Prep (adsorb.py)"""
    print("\n" + "="*60)
    print(color_text("ADSORPTION SETUP", 'bold').center(60))
    print("="*60)
    print(color_text(
        "Sets up a clean-slab reference, an isolated-adsorbate reference, and one folder "
        "per candidate adsorption site -- run SIESTA in each, then use Stage 2 (Analysis) "
        "for the real DFT adsorption energy.", 'cyan'))
    print()

    struct_file = get_input("Input slab/2D structure FDF file, vacuum along c (-s): ").strip()
    while not os.path.isfile(struct_file):
        print(color_text("File not found!", 'red'))
        struct_file = get_input("Input slab/2D structure FDF file (-s): ").strip()

    calc_file = get_input("Calc.fdf template already configured for this slab (-c): ").strip()
    while not os.path.isfile(calc_file):
        print(color_text("File not found!", 'red'))
        calc_file = get_input("Calc.fdf template file (-c): ").strip()

    pp_path = prompt_pseudo_source(optional=True)

    print("\nAdsorbate: an element symbol (e.g. 'O') or an ASE G2 molecule name (e.g. 'H2O'). "
          "Comma-separated for more than one (e.g. 'O,N,C') to compare several species at the "
          "same site(s). Type 'list' to see all G2 names.")
    adsorbate = get_input("Adsorbate(s): ").strip()
    while adsorbate.lower() == "list" or not adsorbate:
        run_tool("stb-adsorb", ["--list", "--no-intro"])
        adsorbate = get_input("Adsorbate(s): ").strip()

    ml_prerelax_choice = get_input(
        "\nPre-relax the isolated adsorbate(s) with MACE-MP-0 before generating their reference? "
        "Needs the optional 'ml' extra (y/N): ").strip().lower()
    ml_prerelax = ml_prerelax_choice in ('y', 'yes')

    print(f"\n{color_text('Site type:', 'yellow')}")
    print(f"  {color_text('1', 'cyan')} = ontop")
    print(f"  {color_text('2', 'cyan')} = bridge")
    print(f"  {color_text('3', 'cyan')} = hollow")
    print(f"  {color_text('4', 'cyan')} = all")
    site_choice = get_input("Select option (1-4) [default: 1]: ").strip()
    site_map = {'1': 'ontop', '2': 'bridge', '3': 'hollow', '4': 'all'}
    site_type = site_map.get(site_choice, 'ontop')

    height_sweep = None
    height_sweep_choice = get_input(
        "\nSweep across multiple heights (approach curve), instead of one fixed height? "
        "(y/N): ").strip().lower()
    if height_sweep_choice in ('y', 'yes'):
        h_min = get_float_input("  Minimum height, Ang [default: 1.5]: ", 1.5)
        h_max = get_float_input("  Maximum height, Ang [default: 3.0]: ", 3.0)
        h_step = get_float_input("  Step, Ang [default: 0.5]: ", 0.5)
        height_sweep = (h_min, h_max, h_step)
        height = h_min  # unused when height_sweep is set, kept for the summary/args fallback
    else:
        height = get_float_input("\nAdsorption height above the surface, Ang [default: 2.0]: ", 2.0)

    all_sites_choice = get_input(
        "\nWrite a folder for every symmetrically distinct site (Y), or just one (n)? [Y/n]: "
    ).strip().lower()
    all_sites = all_sites_choice not in ('n', 'no')

    ml_rank = False
    top_k = None
    if all_sites:
        ml_rank_choice = get_input(
            "\nPre-screen sites with a MACE-MP-0 relax before writing SIESTA folders? Needs "
            "the optional 'ml' extra (y/N): ").strip().lower()
        ml_rank = ml_rank_choice in ('y', 'yes')
        if ml_rank:
            top_k_str = get_input(
                "  Only keep the N best-ranked sites (blank = keep all): ").strip()
            top_k = int(top_k_str) if top_k_str.isdigit() else None
    else:
        site_index = get_int_input("Which site (0-based index) [default: 0]: ", 0)

    both_sides = False
    if site_type != 'all' and not ml_rank and height_sweep is None:
        both_sides_choice = get_input(
            "\nAdsorb on both faces (free-standing 2D material)? (y/N): ").strip().lower()
        both_sides = both_sides_choice in ('y', 'yes')

    bsse_choice = get_input(
        "\nGenerate BSSE (counterpoise) -corrected references per site too? (Y/n -- doubles "
        "each site's SIESTA folder count; LCAO adsorption energies are otherwise "
        "systematically over-bound): ").strip().lower()
    bsse_correction = bsse_choice not in ('n', 'no')

    # Advanced settings (rarely-touched -- gated so the essential flow above
    # stays short; CLI defaults apply untouched when skipped).
    symprec, vacuum_gap, vacuum_box, output_dir = 0.01, 10.0, 20.0, "."
    show_advanced = get_input(
        "\nConfigure advanced settings (symmetry tolerance, vacuum-gap, vacuum box, "
        "output directory)? [y/N]: ").strip().lower()
    if show_advanced == 'y':
        symprec = get_float_input("Site symmetry-reduction tolerance [default: 0.01]: ", 0.01)
        vacuum_gap = get_float_input(
            "Vacuum-axis detection threshold in Ang [default: 10.0]: ", 10.0)
        vacuum_box = get_float_input(
            "Isolated-adsorbate vacuum box side in Ang [default: 20.0]: ", 20.0)
        output_dir = get_input("Output root directory [default: current directory]: ").strip()
        if not output_dir:
            output_dir = "."

    args = [
        "-s", struct_file,
        "-c", calc_file,
        "--adsorbate", adsorbate,
        "--site-type", site_type,
        "--symprec", str(symprec),
        "--vacuum-gap", str(vacuum_gap),
        "--vacuum-box", str(vacuum_box),
        "--output-dir", output_dir,
        "--no-intro",
    ]
    if height_sweep is not None:
        args.extend(["--height-sweep", str(height_sweep[0]), str(height_sweep[1]), str(height_sweep[2])])
    else:
        args.extend(["--height", str(height)])
    if pp_path:
        args.extend(["-p", pp_path])
    if ml_prerelax:
        args.append("--ml-prerelax")
    if all_sites:
        args.append("--all-sites")
        if ml_rank:
            args.append("--ml-rank")
            if top_k is not None:
                args.extend(["--top-k", str(top_k)])
    else:
        args.extend(["--site-index", str(site_index)])
    if both_sides:
        args.append("--both-sides")
    args.append("--bsse-correction" if bsse_correction else "--no-bsse-correction")

    summary_rows = [
        ("Structure file", struct_file),
        ("Calc template", calc_file),
        ("Pseudopotentials", pp_path or "(none)"),
        ("Adsorbate(s)", adsorbate),
        ("ML pre-relax adsorbate", "ON" if ml_prerelax else "OFF"),
        ("Site type", site_type),
        ("Sites", "all symmetrically distinct" if all_sites else f"index {site_index}"),
        ("ML pre-screen", f"ON (top {top_k or 'all'})" if ml_rank else "OFF"),
        ("Both faces", "yes" if both_sides else "no"),
        ("BSSE correction", "ON" if bsse_correction else "OFF"),
        ("Height", f"sweep {height_sweep[0]}-{height_sweep[1]} step {height_sweep[2]} Ang"
                   if height_sweep is not None else f"{height} Ang"),
        ("Output directory", output_dir),
    ]
    _print_config_summary("CONFIGURATION SUMMARY", summary_rows)

    run_tool("stb-adsorb", args)


def run_adsorb_analysis() -> None:
    """Interface for the Adsorption Analysis (adsorb_analysis.py)"""
    print("\n" + "="*60)
    print(color_text("ADSORPTION ANALYSIS", 'bold').center(60))
    print("="*60)
    print(color_text(
        "Reads 'clean_slab/', 'adsorbate/' and every 'sites/site_*/' folder and reports "
        "E_ads = E_site - E_clean_slab - E_adsorbate, ranked most stable first.", 'cyan'))
    print()

    dir_path = get_input("Root directory with 'clean_slab'/'adsorbate'/'sites' [default: .]: ").strip()
    if not dir_path:
        dir_path = "."

    out_file = get_input("SIESTA output filename inside each folder [default: calc.out]: ").strip()
    if not out_file:
        out_file = "calc.out"

    force_tolerance = get_float_input(
        "Force tolerance for the 'is this site relaxed' check, in eV/Ang "
        "(default: 0.05): ", 0.05)

    args = ["--dir", dir_path, "--file", out_file,
            "--force-tolerance", str(force_tolerance), "--no-intro"]

    apply_target = get_input(
        "\nCopy the most stable site's structure.fdf to a production file? "
        "Path to write, or leave blank to skip: ").strip()
    if apply_target:
        args.extend(["--apply", apply_target])

    run_tool("stb-adsorbAnalysis", args)


def run_neb_setup() -> None:
    """Interface for the NEB Prep (neb.py)"""
    print("\n" + "="*60)
    print(color_text("NEB SETUP", 'bold').center(60))
    print("="*60)
    print(color_text(
        "Interpolates a reaction-path band between two already-relaxed endpoint structures "
        "and writes one single-point SIESTA folder per image -- run SIESTA in each, then use "
        "Stage 2 (Analysis) for the DFT-level energy profile.", 'cyan'))
    print()

    initial_file = get_input("Initial (already-relaxed) endpoint structure FDF file (-i): ").strip()
    while not os.path.isfile(initial_file):
        print(color_text("File not found!", 'red'))
        initial_file = get_input("Initial endpoint structure FDF file (-i): ").strip()

    final_file = get_input("Final (already-relaxed) endpoint structure FDF file (-f): ").strip()
    while not os.path.isfile(final_file):
        print(color_text("File not found!", 'red'))
        final_file = get_input("Final endpoint structure FDF file (-f): ").strip()

    calc_file = get_input("Calc.fdf template (-c): ").strip()
    while not os.path.isfile(calc_file):
        print(color_text("File not found!", 'red'))
        calc_file = get_input("Calc.fdf template file (-c): ").strip()

    pp_path = prompt_pseudo_source(optional=True)

    n_images = get_int_input("\nTotal images along the band, endpoints included [default: 7]: ", 7)

    idpp_choice = get_input(
        "\nRefine the interpolated path with ASE's IDPP method (better initial guess, no "
        "MACE needed)? (y/N): ").strip().lower()
    idpp = idpp_choice in ('y', 'yes')

    ml_neb_choice = get_input(
        "\nRun a real climbing-image NEB on MACE-MP-0 before writing SIESTA folders? Needs "
        "the optional 'ml' extra (y/N): ").strip().lower()
    ml_neb = ml_neb_choice in ('y', 'yes')
    ml_k = 0.1
    ml_max_steps = 200
    ml_freeze_substrate = True
    ml_freeze_threshold = 0.3
    if ml_neb:
        ml_k = get_float_input("  NEB spring constant, eV/Ang^2 [default: 0.1]: ", 0.1)
        ml_max_steps = get_int_input("  Max optimizer steps [default: 200]: ", 200)
        freeze_choice = get_input(
            "  Freeze atoms that barely move between endpoints (faster, avoids spurious "
            "drift)? [Y/n]: ").strip().lower()
        ml_freeze_substrate = freeze_choice not in ('n', 'no')
        if ml_freeze_substrate:
            ml_freeze_threshold = get_float_input(
                "    Displacement threshold below which an atom is frozen, Ang [default: 0.3]: ", 0.3)

    ml_prerelax_choice = get_input(
        "\nPre-relax both endpoints with MACE-MP-0 before interpolating? Independent of the "
        "ML-NEB choice above -- needs the optional 'ml' extra (y/N): ").strip().lower()
    ml_prerelax_endpoints = ml_prerelax_choice in ('y', 'yes')

    # Advanced settings (rarely-touched -- gated so the essential flow above
    # stays short; CLI defaults apply untouched when skipped).
    autosort_tol, output_dir = 0.5, "."
    show_advanced = get_input(
        "\nConfigure advanced settings (atom-correspondence tolerance, output directory)? "
        "[y/N]: ").strip().lower()
    if show_advanced == 'y':
        autosort_tol = get_float_input(
            "Atom-correspondence tolerance for interpolation, Ang [default: 0.5]: ", 0.5)
        output_dir = get_input("Output root directory [default: current directory]: ").strip()
        if not output_dir:
            output_dir = "."

    args = [
        "-i", initial_file,
        "-f", final_file,
        "-c", calc_file,
        "-n", str(n_images),
        "--autosort-tol", str(autosort_tol),
        "--output-dir", output_dir,
        "--no-intro",
    ]
    if pp_path:
        args.extend(["-p", pp_path])
    if idpp:
        args.append("--idpp")
    if ml_neb:
        args.extend(["--ml-neb", "--ml-k", str(ml_k), "--ml-max-steps", str(ml_max_steps)])
        if ml_freeze_substrate:
            args.extend(["--ml-freeze-substrate", "--ml-freeze-threshold", str(ml_freeze_threshold)])
        else:
            args.append("--no-ml-freeze-substrate")
    if ml_prerelax_endpoints:
        args.append("--ml-prerelax-endpoints")

    summary_rows = [
        ("Initial structure", initial_file),
        ("Final structure", final_file),
        ("Calc template", calc_file),
        ("Pseudopotentials", pp_path or "(none)"),
        ("N images", n_images),
        ("IDPP refinement", "ON" if idpp else "OFF"),
        ("ML-NEB (MACE-MP-0)", f"ON (k={ml_k}, max steps={ml_max_steps})" if ml_neb else "OFF"),
        ("ML-NEB freeze substrate", (f"ON (threshold {ml_freeze_threshold} Ang)" if ml_freeze_substrate else "OFF")
                                     if ml_neb else "n/a"),
        ("ML pre-relax endpoints", "ON" if ml_prerelax_endpoints else "OFF"),
        ("Output directory", output_dir),
    ]
    _print_config_summary("CONFIGURATION SUMMARY", summary_rows)

    run_tool("stb-neb", args)


def run_neb_analysis() -> None:
    """Interface for the NEB Analysis (neb_analysis.py)"""
    print("\n" + "="*60)
    print(color_text("NEB ANALYSIS", 'bold').center(60))
    print("="*60)
    print(color_text(
        "Reads every 'image_NN/' folder and reports the DFT energy profile along the band, "
        "with an approximate reaction-barrier estimate from the highest-energy image.", 'cyan'))
    print()

    dir_path = get_input("Root directory with every 'image_NN/' [default: .]: ").strip()
    if not dir_path:
        dir_path = "."

    out_file = get_input("SIESTA output filename inside each folder [default: calc.out]: ").strip()
    if not out_file:
        out_file = "calc.out"

    force_tolerance = get_float_input(
        "Force tolerance for the 'is this image single-point-converged' check, in eV/Ang "
        "(default: 0.05): ", 0.05)

    args = ["--dir", dir_path, "--file", out_file,
            "--force-tolerance", str(force_tolerance), "--no-intro"]

    apply_target = get_input(
        "\nCopy the highest-energy image's structure.fdf (transition-state guess) to a "
        "production file? Path to write, or leave blank to skip: ").strip()
    if apply_target:
        args.extend(["--apply", apply_target])

    run_tool("stb-nebAnalysis", args)


def run_stackingfault_setup() -> None:
    """Interface for the Stacking Fault Prep (stackingfault.py)"""
    print("\n" + "="*60)
    print(color_text("STACKING FAULT (GAMMA-SURFACE) SETUP", 'bold').center(60))
    print("="*60)
    print(color_text(
        "Rigidly slides one layer of a bilayer across a 2D grid of lateral offsets and writes "
        "one single-point SIESTA folder per grid point -- run SIESTA in each, then use Stage 2 "
        "(Analysis) for the equilibrium stacking and the stacking-fault energy landscape.", 'cyan'))
    print()

    layer1_file = get_input("Bottom monolayer FDF file (-l1): ").strip()
    while not os.path.isfile(layer1_file):
        print(color_text("File not found!", 'red'))
        layer1_file = get_input("Bottom monolayer FDF file (-l1): ").strip()

    print(color_text(
        "\nTop monolayer: pass the SAME file as above for the canonical use case (a material "
        "sliding against itself, e.g. graphite ABA vs. ABC stacking); a different file studies "
        "an interlayer/heterostructure sliding energy landscape instead.", 'cyan'))
    layer2_file = get_input("Top monolayer FDF file (-l2): ").strip()
    while not os.path.isfile(layer2_file):
        print(color_text("File not found!", 'red'))
        layer2_file = get_input("Top monolayer FDF file (-l2): ").strip()

    calc_file = get_input("Calc.fdf template (-c): ").strip()
    while not os.path.isfile(calc_file):
        print(color_text("File not found!", 'red'))
        calc_file = get_input("Calc.fdf template file (-c): ").strip()

    pp_path = prompt_pseudo_source(optional=True)

    grid_n = get_int_input("\nGrid resolution, N x N shifts [default: 7]: ", 7)
    gap = get_float_input("Interlayer gap, Ang (fixed across the grid, unless --ml-relax-gap "
                           "below) [default: 3.2]: ", 3.2)

    ml_prerelax_choice = get_input(
        "\nPre-relax each monolayer with MACE-MP-0 before stacking? Needs the optional 'ml' "
        "extra (y/N): ").strip().lower()
    ml_prerelax_layers = ml_prerelax_choice in ('y', 'yes')

    ml_relax_gap_choice = get_input(
        "\nOptimize the interlayer gap per grid point with MACE-MP-0, instead of the fixed "
        "gap above? Needs the optional 'ml' extra (y/N): ").strip().lower()
    ml_relax_gap = ml_relax_gap_choice in ('y', 'yes')
    ml_relax_gap_window = 1.0
    if ml_relax_gap:
        ml_relax_gap_window = get_float_input(
            "  Scan window, +/- Ang around the gap above [default: 1.0]: ", 1.0)

    ml_preview_choice = get_input(
        "\nGenerate a fast MACE-MP-0 preview of the gamma-surface (no SIESTA) before writing "
        "the grid? Needs the optional 'ml' extra (y/N): ").strip().lower()
    ml_preview = ml_preview_choice in ('y', 'yes')

    # Advanced settings (rarely-touched -- gated so the essential flow above
    # stays short; CLI defaults apply untouched when skipped).
    max_area, max_strain, match_id, vacuum, strain_mode, twist, output_dir = \
        150.0, 0.05, 0, None, "top", 0.0, "."
    show_advanced = get_input(
        "\nConfigure advanced settings (ZSL match tolerance, vacuum, strain mode, twist, "
        "output directory)? [y/N]: ").strip().lower()
    if show_advanced == 'y':
        max_area = get_float_input("Max ZSL supercell area, Ang^2 [default: 150.0]: ", 150.0)
        max_strain = get_float_input("Max ZSL match strain fraction [default: 0.05]: ", 0.05)
        match_id = get_int_input("Which ZSL match to use, 0-based [default: 0]: ", 0)
        vacuum_str = get_input("Target vacuum, Ang [default: inherit layer 1]: ").strip()
        vacuum = float(vacuum_str) if vacuum_str else None
        strain_mode_choice = get_input("Strain mode top/bottom/sym [default: top]: ").strip().lower()
        strain_mode = strain_mode_choice if strain_mode_choice in ('top', 'bottom', 'sym') else 'top'
        twist = get_float_input("Twist angle, degrees (fixed) [default: 0.0]: ", 0.0)
        output_dir = get_input("Output root directory [default: current directory]: ").strip()
        if not output_dir:
            output_dir = "."

    args = [
        "-l1", layer1_file,
        "-l2", layer2_file,
        "-c", calc_file,
        "-n", str(grid_n),
        "-g", str(gap),
        "-a", str(max_area),
        "-s", str(max_strain),
        "-id", str(match_id),
        "-sm", strain_mode,
        "-t", str(twist),
        "--output-dir", output_dir,
        "--no-intro",
    ]
    if pp_path:
        args.extend(["-p", pp_path])
    if vacuum is not None:
        args.extend(["--vacuum", str(vacuum)])
    if ml_prerelax_layers:
        args.append("--ml-prerelax-layers")
    if ml_relax_gap:
        args.extend(["--ml-relax-gap", "--ml-relax-gap-window", str(ml_relax_gap_window)])
    if ml_preview:
        args.append("--ml-preview")

    summary_rows = [
        ("Layer 1", layer1_file),
        ("Layer 2", layer2_file),
        ("Calc template", calc_file),
        ("Pseudopotentials", pp_path or "(none)"),
        ("Grid", f"{grid_n} x {grid_n}"),
        ("Gap (fixed)", f"{gap} Ang"),
        ("Twist (fixed)", f"{twist} deg"),
        ("ML pre-relax layers", "ON" if ml_prerelax_layers else "OFF"),
        ("ML relax gap", f"ON (window +/-{ml_relax_gap_window} Ang)" if ml_relax_gap else "OFF"),
        ("ML preview", "ON" if ml_preview else "OFF"),
        ("Output directory", output_dir),
    ]
    _print_config_summary("CONFIGURATION SUMMARY", summary_rows)

    run_tool("stb-stackingfault", args)


def run_stackingfault_analysis() -> None:
    """Interface for the Stacking Fault Analysis (stackingfault_analysis.py)"""
    print("\n" + "="*60)
    print(color_text("STACKING FAULT (GAMMA-SURFACE) ANALYSIS", 'bold').center(60))
    print("="*60)
    print(color_text(
        "Reads every 'shift_II_JJ/' folder and reports the equilibrium stacking, the "
        "highest-energy registry, the corrugation energy, and the full 2D gamma-surface map.",
        'cyan'))
    print()

    dir_path = get_input("Root directory with every 'shift_II_JJ/' [default: .]: ").strip()
    if not dir_path:
        dir_path = "."

    out_file = get_input("SIESTA output filename inside each folder [default: calc.out]: ").strip()
    if not out_file:
        out_file = "calc.out"

    force_tolerance = get_float_input(
        "Force tolerance for the 'is this point single-point-converged' check, in eV/Ang "
        "(default: 0.05): ", 0.05)

    args = ["--dir", dir_path, "--file", out_file,
            "--force-tolerance", str(force_tolerance), "--no-intro"]

    apply_target = get_input(
        "\nCopy the equilibrium (lowest-energy) point's structure.fdf to a production file? "
        "Path to write, or leave blank to skip: ").strip()
    if apply_target:
        args.extend(["--apply", apply_target])

    run_tool("stb-stackingfaultAnalysis", args)


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

    # 6. Output file name
    output_file = get_input("\nOutput file name [default: stacked_structure.fdf]: ").strip()
    if not output_file:
        output_file = "stacked_structure.fdf"

    # 7. Build Command and Execute ONCE
    args = [
        "-l1", layer1,
        "-l2", layer2,
        "-i", # Always included as requested
        "-a", str(max_area),
        "-s", str(max_strain),
        "-sm", strain_mode,
        "-o", output_file,
        "--no-intro"
    ]

    # Adiciona os argumentos de Gap (seja -g ou --gap_range)
    args.extend(gap_args)

    # Adiciona os argumentos de Simetria e Empilhamento
    if batch_sym:
        args.append("--batch_sym")
    else:
        args.extend(["-t", str(twist), "-tx", str(shift_x), "-ty", str(shift_y)])

    ml_relax = get_input(
        "\nPre-relax the final structure(s) with a MACE ML potential before writing? "
        "(needs the optional 'ml' extra) (y/N): "
    ).strip().lower()
    if ml_relax == 'y':
        args.append("--ml-relax")

        ml_relax_cell = get_input(
            "Also relax the in-plane cell? The vacuum axis (c) always stays fixed (y/N): "
        ).strip().lower()
        if ml_relax_cell == 'y':
            args.append("--ml-relax-cell")

        custom_model = get_input(
            "Use a custom MACE model instead (e.g. one fine-tuned via stb-mlffAnalysis)? "
            "Path, or Enter to use a MACE-MP-0 foundation model: ").strip()
        if custom_model:
            while not os.path.isfile(custom_model):
                print(color_text("File not found!", 'red'))
                custom_model = get_input("Custom model path: ").strip()
            args.extend(["--custom-model", custom_model])
        else:
            model = get_input("Model size, small/medium/large [default: small]: ").strip()
            if not model:
                model = "small"
            args.extend(["--model", model])

    save_report = get_input("Also save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")

    view_choice = get_input("View the structure interactively via ASE before finishing? (y/N): ").strip().lower()
    if view_choice == 'y':
        args.append("--view")

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

    # 3. Canal de spin (só importa para cálculos spin-polarizados; ignorado
    # com aviso caso o arquivo não seja spin-polarizado)
    print(f"\n{color_text('Spin channel (only matters for spin-polarized calculations):', 'yellow')}")
    print(f"  {color_text('1', 'cyan')} = Total (up + down)")
    print(f"  {color_text('2', 'cyan')} = Up")
    print(f"  {color_text('3', 'cyan')} = Down")
    print(f"  {color_text('4', 'cyan')} = Diff (up - down, i.e. spin/magnetization density)")
    spin_map = {'1': 'total', '2': 'up', '3': 'down', '4': 'diff'}
    spin_choice = get_input("Select spin channel (1-4) [default: 1]: ").strip()
    selected_spin = spin_map.get(spin_choice, 'total')

    # 4. Nome de saída (opcional)
    default_output = (f"{label}_{selected_type}.cube" if selected_spin == 'total'
                      else f"{label}_{selected_type}_{selected_spin}.cube")
    output = get_input(f"Output filename (default: {default_output}): ").strip()

    # 5. Execução
    args = ["--label", label, "--type", selected_type, "--spin", selected_spin, "--no-intro"]
    if output:
        args.extend(["--output", output])

    print(color_text("\nConverting to Cube format...", 'green'))
    run_tool("stb-cube", args)


def run_ani2traj_converter() -> None:
    """Interface for the AIMD Trajectory Converter (ani2traj.py)"""
    print("\n" + "="*60)
    print(color_text("SIESTA AIMD TRAJECTORY CONVERTER", 'bold').center(60))
    print("="*60 + "\n")
    print("Converts <label>.ANI (SIESTA's AIMD trajectory, no lattice of its own) into a")
    print("multi-frame format that also carries the cell -- for OVITO, VMD, etc. The")
    print("lattice is read from <label>.XV (or <label>.fdf as a fallback).\n")

    # 1. SystemLabel
    label = get_input("Enter the Siesta SystemLabel (e.g., siesta): ").strip()
    while not label:
        print(color_text("Label cannot be empty!", 'red'))
        label = get_input("Enter the Siesta SystemLabel: ")

    # 2. Formato de saída
    print(f"\n{color_text('Select output trajectory format:', 'yellow')}")
    print(f"  {color_text('1', 'cyan')} = XSF (multi-frame AXSF -- OVITO and VMD's xsf plugin)")
    print(f"  {color_text('2', 'cyan')} = PDB (multi-model, CRYST1 -- VMD's own default format)")
    print(f"  {color_text('3', 'cyan')} = XYZ (extended, Lattice tag -- OVITO-native; no PBC in VMD)")
    format_map = {'1': 'xsf', '2': 'pdb', '3': 'xyz'}
    ext_map = {'xsf': '.xsf', 'pdb': '.pdb', 'xyz': '.xyz'}
    choice = get_input("Select format (1-3) [default: 1]: ").strip()
    selected_format = format_map.get(choice, 'xsf')

    # 3. Stride (opcional)
    stride = get_int_input("Keep every Nth frame (stride) [default: 1]: ", 1)

    # 4. Nome de saída (opcional)
    default_output = f"{label}_traj{ext_map[selected_format]}"
    output = get_input(f"Output filename (default: {default_output}): ").strip()

    args = ["--label", label, "--out-format", selected_format, "--stride", str(stride), "--no-intro"]
    if output:
        args.extend(["--output", output])

    print(color_text("\nConverting AIMD trajectory...", 'green'))
    run_tool("stb-ani2traj", args)


def run_aimd_analyzer() -> None:
    """Interface for the AIMD Trajectory Analysis tool (aimd_analysis.py)"""
    print("\n" + "="*60)
    print(color_text("AIMD TRAJECTORY ANALYSIS", 'bold').center(60))
    print("="*60 + "\n")
    print("Radial distribution function (RDF), mean-squared displacement (MSD) /")
    print("diffusion coefficient, and a VACF-derived vibrational density of states")
    print("(VDOS) from a SIESTA AIMD trajectory (<label>.ANI + .out) or a generic")
    print("ASE-readable trajectory (e.g. one written by stb-mlmd).\n")

    print(f"{color_text('Input source:', 'yellow')}")
    print(f"  {color_text('1', 'cyan')} = SIESTA SystemLabel (<label>.ANI, default)")
    print(f"  {color_text('2', 'cyan')} = Generic trajectory file (e.g. stb-mlmd output)")
    source_choice = get_input("Select option (1-2) [default: 1]: ").strip()

    args = ["--no-intro"]
    if source_choice == '2':
        trajectory = get_input("Trajectory file path: ").strip()
        while not os.path.isfile(trajectory):
            print(color_text("File not found!", 'red'))
            trajectory = get_input("Trajectory file path: ").strip()
        args.extend(["--trajectory", trajectory])
        dt = get_input(
            "Real time between saved frames, fs (blank = auto-detect from embedded "
            "'Time' info, only works for an extended-xyz file): ").strip()
        if dt:
            args.extend(["--dt", dt])
    else:
        label = get_input("Enter the Siesta SystemLabel (e.g., siesta): ").strip()
        while not label:
            print(color_text("Label cannot be empty!", 'red'))
            label = get_input("Enter the Siesta SystemLabel: ")
        args.extend(["--label", label])

        print(f"{color_text('Note:', 'yellow')} the real .fdf input file is very often NOT "
              f"named '{label}.fdf' (SystemLabel is chosen independently of the input "
              "filename) -- needed here for the true MD timestep (falls back to 1 fs/frame "
              "with a warning if not given).")
        default_fdf = f"{label}.fdf"
        if os.path.isfile(default_fdf):
            geometry_file = get_input(
                f"Path to the real .fdf input [default: {default_fdf}]: ").strip() or default_fdf
        else:
            geometry_file = get_input(
                "Path to the real .fdf input (blank to skip, assumes 1 fs/frame): ").strip()
        if geometry_file:
            args.extend(["--geometry-file", geometry_file])

    stride = get_int_input("Keep every Nth frame (stride) [default: 1]: ", 1)
    skip = get_int_input("Discard the first N frames as equilibration (--skip) [default: 0]: ", 0)

    pair = get_input("Restrict the RDF to one species pair, e.g. 'O-H' (default: all pairs): ").strip()

    list_atoms = get_input(
        "List every atom's index/species/coordinates, to help pick --track-atom/"
        "--track-pair? (y/N): ").strip().lower()
    if list_atoms == 'y':
        run_tool("stb-aimdAnalysis", args + ["--list-atoms"], pause=False)

    track_atom = get_input(
        "0-based atom index to track its own displacement over time (blank to skip): ").strip()
    track_pair = get_input(
        "Two 0-based atom indices to track their relative distance, e.g. '0-1' "
        "(blank to skip): ").strip()

    output_dir = get_input("Output directory (default: '.'): ") or "."
    args.extend(["--output-dir", output_dir])

    save_report = get_input("Also save a text report to file? (y/N): ").strip().lower()
    save_gnuplot = get_input("Also save the data + gnuplot scripts? (y/N): ").strip().lower()
    view_choice = get_input("Show the matplotlib plot before finishing? (y/N): ").strip().lower()

    args.extend(["--stride", str(stride), "--skip", str(skip)])
    if pair:
        args.extend(["--pair", pair])
    if track_atom:
        args.extend(["--track-atom", track_atom])
    if track_pair:
        args.extend(["--track-pair", track_pair])
    if save_report == 'y':
        args.append("--save-report")
    if save_gnuplot == 'y':
        args.append("--save-gnuplot")
    if view_choice == 'y':
        args.append("--view")

    print(color_text("\nAnalyzing AIMD trajectory...", 'green'))
    run_tool("stb-aimdAnalysis", args)


def run_status_checker() -> None:
    """Interface for the Status Checker (status.py)"""
    print("\n" + "="*60)
    print(color_text("SIESTA STATUS CHECKER", 'bold').center(60))
    print("="*60 + "\n")
    print("Reports run type, SCF convergence, max force, final energy, and which key")
    print("output files are present in a calc folder.\n")

    path = get_input("Directory (or glob pattern, e.g. strain_*) to inspect (default: current): ").strip()
    if path == "":
        path = "."
    is_glob = any(c in path for c in '*?[')
    while not is_glob and not os.path.isdir(path):
        print(color_text("Directory not found!", 'red'))
        path = get_input("Enter a valid directory (or glob pattern): ").strip()
        is_glob = any(c in path for c in '*?[')

    label = get_input("SystemLabel override (leave blank to auto-detect): ").strip()

    args = ["--path", path, "--no-intro"]
    if label:
        args.extend(["--label", label])

    run_tool("stb-status", args)


def run_archive_packager() -> None:
    """Interface for the Archive Packager (archive.py)"""
    print("\n" + "="*60)
    print(color_text("SIESTA ARCHIVE PACKAGER", 'bold').center(60))
    print("="*60 + "\n")
    print("Packages a finished calculation (inputs + essential outputs) into a single")
    print(".tar.gz for sharing/reproducibility -- the complement of stb-clean, which")
    print("deletes that same class of outputs to prep for a restart instead.\n")

    path = get_input("Directory (or glob pattern, e.g. strain_*) to archive (default: current): ").strip()
    if path == "":
        path = "."
    is_glob = any(c in path for c in '*?[')
    while not is_glob and not os.path.isdir(path):
        print(color_text("Directory not found!", 'red'))
        path = get_input("Enter a valid directory (or glob pattern): ").strip()
        is_glob = any(c in path for c in '*?[')

    extra_include = get_input(
        "Extra extensions to include, space-separated (e.g. .RHO), or Enter to skip: ").strip()
    extra_exclude = get_input(
        "Extensions to exclude from the defaults, space-separated, or Enter to skip: ").strip()

    args = ["--path", path, "--no-intro"]
    if extra_include:
        args.extend(["--include"] + extra_include.split())
    if extra_exclude:
        args.extend(["--exclude"] + extra_exclude.split())

    run_tool("stb-archive", args)


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

    # Densidade de spin: agora detectada e reportada automaticamente sempre que o
    # .RHO for spin-polarizado -- esta pergunta só decide se o usuário quer ver
    # APENAS o spin (pulando a carga), não se o spin é ou não processado.
    spin_choice = get_input(
        "Process ONLY the spin density, skipping charge (only if the .RHO is "
        "spin-polarized -- otherwise spin is detected and shown automatically "
        "alongside charge)? (y/N): ", 'green').strip().lower()
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

    output_dir = get_input("\nOutput directory [default: current directory]: ").strip()
    if output_dir:
        args.extend(["--output-dir", output_dir])

    save_report = get_input("Also save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")

    save_gnuplot = get_input("Also save the data + gnuplot script(s)? (y/N): ").strip().lower()
    if save_gnuplot == 'y':
        args.append("--save-gnuplot")

    view_choice = get_input("Show the matplotlib preview(s) before finishing? (y/N): ").strip().lower()
    if view_choice == 'y':
        args.append("--view")

    print(color_text("\nProcessing Density...", 'green'))
    run_tool("stb-density", args)

def run_stm_analyzer() -> None:
    """Interface for the STM Simulator (stb-stm)"""
    print("\n" + "="*60)
    print(color_text("STM SIMULATOR (stb-stm)", 'bold').center(60))
    print("="*60 + "\n")

    print(f"This tool reads the {color_text('.LDOS', 'yellow')} file generated by Siesta's "
          f"{color_text('%block LocalDensityOfStates', 'yellow')}.")
    label = get_input("Enter the Siesta SystemLabel (e.g., siesta): ").strip()
    while not label:
        print(color_text("Label cannot be empty!", 'red'))
        label = get_input("Enter the Siesta SystemLabel: ")

    if os.path.isfile(f"{label}.LDOS"):
        print(color_text(f"-> Found '{label}.LDOS'.", 'green'))
    elif os.path.isfile(f"{label}.STM.LDOS"):
        print(color_text(f"-> Found '{label}.STM.LDOS'.", 'green'))
    else:
        print(color_text(f"-> No '{label}.LDOS'/'{label}.STM.LDOS' found -- stb-stm will fail "
                          "without one. Add '%block LocalDensityOfStates' to the fdf and rerun SIESTA.",
                          'yellow'))

    print(f"\n{color_text('STM Mode:', 'yellow')}")
    print(f"  {color_text('1', 'cyan')} = Constant-current (classic STM image)")
    print(f"  {color_text('2', 'cyan')} = Constant-height (LDOS map at a fixed height)")
    mode_choice = get_input("Select mode (1-2) [default: 1]: ").strip()

    args = ["--label", label, "--no-intro"]

    # Defaults below must match stm.py's own DEFAULT_Z/DEFAULT_ISO.
    if mode_choice == '2':
        args.append("--mode")
        args.append("height")
        z = get_float_input(
            "Height above the topmost surface atom, in Angstrom [default: 3.0]: ", 3.0)
        args.extend(["--z", str(z)])
    else:
        args.append("--mode")
        args.append("current")
        iso = get_float_input(
            "LDOS threshold (e/Bohr^3) defining the constant-current contour "
            "[default: 0.001]: ", 0.001)
        args.extend(["--iso", str(iso)])

    axis_str = get_input("Surface-normal axis, 0/1/2 (blank = auto-detect from vacuum padding): ").strip()
    if axis_str:
        args.extend(["--axis", axis_str])

    output_dir = get_input("Output directory (default: '.'): ") or "."
    args.extend(["--output-dir", output_dir])

    save_report = get_input("\nAlso save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")
    save_gnuplot = get_input("Also save the data + gnuplot script? (y/N): ").strip().lower()
    if save_gnuplot == 'y':
        args.append("--save-gnuplot")
    view_choice = get_input("Show the matplotlib plot before finishing? (y/N): ").strip().lower()
    if view_choice == 'y':
        args.append("--view")

    run_tool("stb-stm", args)


def run_wfdensity_analyzer() -> None:
    """Interface for the Wavefunction Density tool (stb-wfdensity)"""
    print("\n" + "="*60)
    print(color_text("WAVEFUNCTION DENSITY (stb-wfdensity)", 'bold').center(60))
    print("="*60 + "\n")

    print(f"This tool needs a {color_text('.WFSX', 'yellow')} file (found via the SystemLabel) and a "
          f"{color_text('.fdf', 'yellow')} geometry (with its .ion/.ion.xml files alongside it) -- "
          "asked separately below, since the real fdf is very often NOT named <label>.fdf (e.g. "
          "SystemLabel 'siesta' with the actual input file called calc.fdf).")
    label = get_input("Enter the Siesta SystemLabel (used to find <label>.WFSX etc.): ").strip()
    while not label:
        print(color_text("Label cannot be empty!", 'red'))
        label = get_input("Enter the Siesta SystemLabel: ")

    if os.path.isfile(f"{label}.WFSX"):
        print(color_text(f"-> Found '{label}.WFSX'.", 'green'))
    elif os.path.isfile(f"{label}.selected.WFSX"):
        print(color_text(f"-> Found '{label}.selected.WFSX'.", 'green'))
    elif os.path.isfile(f"{label}.bands.WFSX"):
        print(color_text(f"-> Found '{label}.bands.WFSX'.", 'green'))
    else:
        print(color_text(f"-> No '{label}.WFSX'/'{label}.selected.WFSX'/'{label}.bands.WFSX' "
                          "found -- stb-wfdensity will fail without one.", 'yellow'))

    default_fdf = f"{label}.fdf"
    fdf_prompt_default = default_fdf if os.path.isfile(default_fdf) else None
    if fdf_prompt_default:
        geometry_file = get_input(
            f"Path to the .fdf structure (with its .ion/.ion.xml files alongside it) "
            f"[default: {default_fdf}]: ").strip() or default_fdf
    else:
        print(color_text(f"-> No '{default_fdf}' found -- enter the real input file's path below "
                          "(e.g. calc.fdf).", 'yellow'))
        geometry_file = get_input(
            "Path to the .fdf structure (with its .ion/.ion.xml files alongside it): ").strip()
        while not geometry_file:
            print(color_text("A geometry file is required!", 'red'))
            geometry_file = get_input("Path to the .fdf structure: ").strip()

    print(f"\n{color_text('Band selection:', 'yellow')}")
    print(f"  {color_text('1', 'cyan')} = Band index N at a chosen k-index")
    print(f"  {color_text('2', 'cyan')} = VBM (needs Fermi energy)")
    print(f"  {color_text('3', 'cyan')} = CBM (needs Fermi energy)")
    choice = get_input("Select (1-3) [default: 1]: ").strip()

    args = ["--label", label, "--geometry-file", geometry_file, "--no-intro"]

    if choice in ('2', '3'):
        args.extend(["--band", "vbm" if choice == '2' else "cbm"])
        print(f"\n{color_text('Fermi energy source:', 'yellow')}")
        print(f"  {color_text('1', 'cyan')} = Enter a value directly")
        print(f"  {color_text('2', 'cyan')} = Read from a .bands file")
        print(f"  {color_text('3', 'cyan')} = Read from a SIESTA .out log (any filename, not just <label>.out)")
        print(f"  {color_text('4', 'cyan')} = Auto-detect a .out in this directory [default]")
        fermi_choice = get_input("Select (1-4) [default: 4]: ").strip() or '4'
        if fermi_choice == '1':
            fermi = get_float_input("Fermi energy (eV): ")
            args.extend(["--fermi", str(fermi)])
        elif fermi_choice == '2':
            bands_file = get_input("Path to the .bands file: ").strip()
            args.extend(["--bands-file", bands_file])
        elif fermi_choice == '3':
            out_file = get_input("Path to the .out log (any filename): ").strip()
            args.extend(["--fermi-file", out_file])
        # fermi_choice == '4': pass nothing, stb-wfdensity auto-detects a .out itself.
    else:
        k_index = get_int_input("k-index (0-based) [default: 0]: ", 0)
        args.extend(["--k-index", str(k_index)])
        band_n = get_int_input("Band number (1-based): ", 1)
        args.extend(["--band", str(band_n)])

    spacing = get_float_input("Grid spacing in Angstrom (default: 0.3): ", 0.3)
    args.extend(["--spacing", str(spacing)])

    print(f"\n{color_text('Cut axis (for the 2D slice / 1D profile):', 'yellow')}")
    print(f"  {color_text('0', 'cyan')} = a   {color_text('1', 'cyan')} = b   "
          f"{color_text('2', 'cyan')} = c [default]")
    axis = get_input("Select (0-2) [default: 2]: ").strip() or "2"
    args.extend(["--axis", axis])

    print(f"\n{color_text('Cut mode:', 'yellow')}")
    print(f"  {color_text('1', 'cyan')} = 2D slice (perpendicular to the cut axis) [default]")
    print(f"  {color_text('2', 'cyan')} = 1D profile (planar-averaged along the cut axis)")
    mode_choice = get_input("Select (1-2) [default: 1]: ").strip()
    if mode_choice == '2':
        args.append("--profile")
    else:
        print(f"\n{color_text('Slice position along the cut axis:', 'yellow')}")
        print(f"  {color_text('1', 'cyan')} = Auto-detect the |psi|^2 density peak [default]")
        print(f"  {color_text('2', 'cyan')} = Enter a position manually (Angstrom)")
        pos_choice = get_input("Select (1-2) [default: 1]: ").strip()
        if pos_choice == '2':
            pos = get_float_input("Position along the cut axis (Ang): ")
            args.extend(["--pos", str(pos)])

    output_dir = get_input("Output directory (default: '.'): ") or "."
    args.extend(["--output-dir", output_dir])

    save_report = get_input("\nAlso save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")
    save_gnuplot = get_input("Also save the slice/profile data + gnuplot script? (y/N): ").strip().lower()
    if save_gnuplot == 'y':
        args.append("--save-gnuplot")
    view_choice = get_input("Show the matplotlib plot before finishing? (y/N): ").strip().lower()
    if view_choice == 'y':
        args.append("--view")

    run_tool("stb-wfdensity", args)


def run_sts_analyzer() -> None:
    """Interface for the STS Spectroscopy Simulator (stb-sts)"""
    print("\n" + "="*60)
    print(color_text("STS SPECTROSCOPY SIMULATOR (stb-sts)", 'bold').center(60))
    print("="*60 + "\n")

    print(f"This tool needs a full-BZ {color_text('.WFSX', 'yellow')} file (not a band-path one, "
          f"found via the SystemLabel) and a {color_text('.fdf', 'yellow')} geometry (with its "
          ".ion/.ion.xml files alongside it) -- asked separately below, since the real fdf is "
          "very often NOT named <label>.fdf (e.g. SystemLabel 'siesta' with the actual input "
          "file called calc.fdf).")
    label = get_input("Enter the Siesta SystemLabel (used to find <label>.WFSX etc.): ").strip()
    while not label:
        print(color_text("Label cannot be empty!", 'red'))
        label = get_input("Enter the Siesta SystemLabel: ")

    if os.path.isfile(f"{label}.WFSX"):
        print(color_text(f"-> Found '{label}.WFSX'.", 'green'))
    elif os.path.isfile(f"{label}.selected.WFSX"):
        print(color_text(f"-> Found '{label}.selected.WFSX'.", 'green'))
    elif os.path.isfile(f"{label}.bands.WFSX"):
        print(color_text(f"-> Found '{label}.bands.WFSX'.", 'green'))
    else:
        print(color_text(f"-> No '{label}.WFSX'/'{label}.selected.WFSX'/'{label}.bands.WFSX' "
                          "found -- stb-sts will fail without one.", 'yellow'))

    default_fdf = f"{label}.fdf"
    fdf_prompt_default = default_fdf if os.path.isfile(default_fdf) else None
    if fdf_prompt_default:
        geometry_file = get_input(
            f"Path to the .fdf structure (with its .ion/.ion.xml files alongside it) "
            f"[default: {default_fdf}]: ").strip() or default_fdf
    else:
        print(color_text(f"-> No '{default_fdf}' found -- enter the real input file's path below "
                          "(e.g. calc.fdf).", 'yellow'))
        geometry_file = get_input(
            "Path to the .fdf structure (with its .ion/.ion.xml files alongside it): ").strip()
        while not geometry_file:
            print(color_text("A geometry file is required!", 'red'))
            geometry_file = get_input("Path to the .fdf structure: ").strip()

    args = ["--label", label, "--geometry-file", geometry_file, "--no-intro"]

    print(f"\n{color_text('Tip position:', 'yellow')}")
    print(f"  {color_text('1', 'cyan')} = X,Y + height above the topmost surface atom (slab-style)")
    print(f"  {color_text('2', 'cyan')} = Absolute Cartesian X,Y,Z")
    choice = get_input("Select (1-2) [default: 1]: ").strip()

    if choice == '2':
        x = get_float_input("X (Ang): ")
        y = get_float_input("Y (Ang): ")
        z = get_float_input("Z (Ang): ")
        args.extend(["--point", str(x), str(y), str(z)])
    else:
        x = get_float_input("X (Ang): ")
        y = get_float_input("Y (Ang): ")
        args.extend(["--xy", str(x), str(y)])
        height = get_float_input("Height above the topmost surface atom (Ang): ")
        args.extend(["--height", str(height)])

    emin = get_float_input("Energy window minimum (eV): ")
    emax = get_float_input("Energy window maximum (eV): ")
    args.extend(["--erange", str(emin), str(emax)])

    sigma = get_float_input("Gaussian broadening sigma, in meV (default: 50): ", 50.0)
    args.extend(["--sigma", str(sigma)])

    print(f"\n{color_text('Energy reference:', 'yellow')}")
    print(f"  {color_text('1', 'cyan')} = Raw eigenvalues, no shift [default]")
    print(f"  {color_text('2', 'cyan')} = Shift by the Fermi energy")
    print(f"  {color_text('3', 'cyan')} = Shift by a custom value")
    shift_choice = get_input("Select (1-3) [default: 1]: ").strip()

    if shift_choice == '2':
        args.extend(["--shift", "fermi"])
        print(f"\n{color_text('Fermi energy source:', 'yellow')}")
        print(f"  {color_text('1', 'cyan')} = Enter a value directly")
        print(f"  {color_text('2', 'cyan')} = Read from a .bands file")
        print(f"  {color_text('3', 'cyan')} = Read from a SIESTA .out log (any filename, not just <label>.out)")
        print(f"  {color_text('4', 'cyan')} = Auto-detect a .out in this directory [default]")
        fermi_choice = get_input("Select (1-4) [default: 4]: ").strip() or '4'
        if fermi_choice == '1':
            fermi = get_float_input("Fermi energy (eV): ")
            args.extend(["--fermi", str(fermi)])
        elif fermi_choice == '2':
            bands_file = get_input("Path to the .bands file: ").strip()
            args.extend(["--bands-file", bands_file])
        elif fermi_choice == '3':
            out_file = get_input("Path to the .out log (any filename): ").strip()
            args.extend(["--fermi-file", out_file])
        # fermi_choice == '4': pass nothing, stb-sts auto-detects a .out itself.
    elif shift_choice == '3':
        args.extend(["--shift", "manual"])
        manual_value = get_float_input("Custom shift value (eV): ")
        args.extend(["--manual-value", str(manual_value)])

    output_dir = get_input("Output directory (default: '.'): ") or "."
    args.extend(["--output-dir", output_dir])

    save_report = get_input("\nAlso save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")
    save_gnuplot = get_input("Also save the sts.dat data + gnuplot script? (y/N): ").strip().lower()
    if save_gnuplot == 'y':
        args.append("--save-gnuplot")
    view_choice = get_input("Show the matplotlib plot before finishing? (y/N): ").strip().lower()
    if view_choice == 'y':
        args.append("--view")

    run_tool("stb-sts", args)


def run_coop_analyzer() -> None:
    """Interface for the COOP/COHP Bonding Analysis tool (stb-coop)"""
    print("\n" + "="*60)
    print(color_text("COOP/COHP BONDING ANALYSIS (stb-coop)", 'bold').center(60))
    print("="*60 + "\n")

    print(f"This tool needs a full-BZ {color_text('.WFSX', 'yellow')} file (found via the "
          f"SystemLabel) and a {color_text('.HSX', 'yellow')} file (mandatory -- no approximate "
          "fallback) -- asked separately below, since the real file is very often NOT named "
          "<label>.HSX.")
    label = get_input("Enter the Siesta SystemLabel (used to find <label>.WFSX etc.): ").strip()
    while not label:
        print(color_text("Label cannot be empty!", 'red'))
        label = get_input("Enter the Siesta SystemLabel: ")

    default_hsx = f"{label}.HSX"
    if os.path.isfile(default_hsx):
        hsx_file = get_input(f"Path to the .HSX file [default: {default_hsx}]: ").strip() or default_hsx
    else:
        print(color_text(f"-> No '{default_hsx}' found -- enter the real .HSX path below.", 'yellow'))
        hsx_file = get_input("Path to the .HSX file: ").strip()
        while not hsx_file:
            print(color_text("An .HSX file is required!", 'red'))
            hsx_file = get_input("Path to the .HSX file: ").strip()

    print(f"\n{color_text('Quantity:', 'yellow')}")
    print(f"  {color_text('1', 'cyan')} = COOP (overlap population)")
    print(f"  {color_text('2', 'cyan')} = COHP (Hamilton population)")
    q_choice = get_input("Select (1-2) [default: 1]: ").strip()

    args = ["--label", label, "--hsx-file", hsx_file,
            "--quantity", "cohp" if q_choice == '2' else "coop", "--no-intro"]

    print(f"\n{color_text('Atom pair selection:', 'yellow')}")
    print(f"  {color_text('1', 'cyan')} = By atom index (e.g. 0 1)")
    print(f"  {color_text('2', 'cyan')} = By species (e.g. Sn O)")
    pair_choice = get_input("Select (1-2) [default: 1]: ").strip()

    if pair_choice == '2':
        sp_a = get_input("Species A: ").strip()
        sp_b = get_input("Species B: ").strip()
        args.extend(["--pair-species", sp_a, sp_b])
    else:
        i = get_int_input("Atom index I (0-based): ", 0)
        j = get_int_input("Atom index J (0-based): ", 1)
        args.extend(["--pair", str(i), str(j)])

    emin = get_float_input("Energy window minimum (eV): ")
    emax = get_float_input("Energy window maximum (eV): ")
    args.extend(["--erange", str(emin), str(emax)])

    sigma = get_float_input("Gaussian broadening sigma, in meV (default: 200): ", 200.0)
    args.extend(["--sigma", str(sigma)])

    print(f"\n{color_text('Energy reference:', 'yellow')}")
    print(f"  {color_text('1', 'cyan')} = Raw eigenvalues, no shift [default]")
    print(f"  {color_text('2', 'cyan')} = Shift by the Fermi energy")
    print(f"  {color_text('3', 'cyan')} = Shift by a custom value")
    shift_choice = get_input("Select (1-3) [default: 1]: ").strip()

    if shift_choice == '2':
        args.extend(["--shift", "fermi"])
        print(f"\n{color_text('Fermi energy source:', 'yellow')}")
        print(f"  {color_text('1', 'cyan')} = Enter a value directly")
        print(f"  {color_text('2', 'cyan')} = Read from a .bands file")
        print(f"  {color_text('3', 'cyan')} = Read from a SIESTA .out log (any filename, not just <label>.out)")
        print(f"  {color_text('4', 'cyan')} = Auto-detect a .out in this directory [default]")
        fermi_choice = get_input("Select (1-4) [default: 4]: ").strip() or '4'
        if fermi_choice == '1':
            fermi = get_float_input("Fermi energy (eV): ")
            args.extend(["--fermi", str(fermi)])
        elif fermi_choice == '2':
            bands_file = get_input("Path to the .bands file: ").strip()
            args.extend(["--bands-file", bands_file])
        elif fermi_choice == '3':
            out_file = get_input("Path to the .out log (any filename): ").strip()
            args.extend(["--fermi-file", out_file])
        # fermi_choice == '4': pass nothing, stb-coop auto-detects a .out itself.
    elif shift_choice == '3':
        args.extend(["--shift", "manual"])
        manual_value = get_float_input("Custom shift value (eV): ")
        args.extend(["--manual-value", str(manual_value)])

    bond_order = get_input("\nAlso compute the Mulliken bond order cross-check? (y/N): ").strip().lower()
    if bond_order == 'y':
        args.append("--bond-order")
        default_dm = f"{label}.DM"
        if os.path.isfile(default_dm):
            dm_file = get_input(f"Path to the .DM density matrix [default: {default_dm}]: ").strip() or default_dm
        else:
            print(color_text(f"-> No '{default_dm}' found -- enter the real .DM path below.", 'yellow'))
            dm_file = get_input("Path to the .DM density matrix: ").strip()
        if dm_file and dm_file != default_dm:
            args.extend(["--dm-file", dm_file])
        # dm_file blank or matching the default: let stb-coop's own <label>.DM
        # auto-detect (or its clean error if none is found) handle it.

    output_dir = get_input("Output directory (default: '.'): ") or "."
    args.extend(["--output-dir", output_dir])

    save_report = get_input("\nAlso save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")
    save_gnuplot = get_input("Also save the coop/cohp.dat data + gnuplot script? (y/N): ").strip().lower()
    if save_gnuplot == 'y':
        args.append("--save-gnuplot")
    view_choice = get_input("Show the matplotlib plot before finishing? (y/N): ").strip().lower()
    if view_choice == 'y':
        args.append("--view")

    run_tool("stb-coop", args)


def run_ipr_analyzer() -> None:
    """Interface for the IPR Analyzer (stb-ipr)"""
    print("\n" + "="*60)
    print(color_text("IPR ANALYZER (stb-ipr)", 'bold').center(60))
    print("="*60 + "\n")

    label = get_input("SIESTA label (e.g. 'siesta' for siesta.bands + siesta.bands.WFSX): ")
    while not os.path.isfile(f"{label}.bands"):
        print(color_text(f"File '{label}.bands' not found!", 'red'))
        label = get_input("SIESTA label: ")

    if os.path.isfile(f"{label}.bands.WFSX"):
        print(color_text(f"-> Found '{label}.bands.WFSX'.", 'green'))
    elif os.path.isfile(f"{label}.WFSX"):
        print(color_text(f"-> Found '{label}.WFSX'.", 'green'))
    elif os.path.isfile(f"{label}.selected.WFSX"):
        print(color_text(f"-> Found '{label}.selected.WFSX'.", 'green'))
    else:
        print(color_text(f"-> No '{label}.bands.WFSX'/'{label}.WFSX'/'{label}.selected.WFSX' "
                          "found -- stb-ipr will fail without one.", 'yellow'))

    args = ["--label", label, "--no-intro"]

    default_hsx = f"{label}.HSX"
    if os.path.isfile(default_hsx):
        hsx_file = get_input(
            f"Path to the .HSX file for overlap-aware (accurate) IPR "
            f"[default: {default_hsx}]: ").strip() or default_hsx
        args.extend(["--hsx-file", hsx_file])
    else:
        print(color_text(f"-> No '{default_hsx}' found -- IPR will use an approximate "
                          "orthogonal-basis fallback unless you point at one below.", 'yellow'))
        hsx_file = get_input(
            "Path to the .HSX file (blank to use the approximate fallback): ").strip()
        if hsx_file:
            args.extend(["--hsx-file", hsx_file])
        else:
            default_fdf = f"{label}.fdf"
            geometry_file = get_input(
                f"Path to the .fdf geometry for the fallback [default: {default_fdf}, "
                "blank for auto-detect]: ").strip()
            if geometry_file:
                args.extend(["--geometry-file", geometry_file])

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
    args.extend(["--shift", shift_type])

    if shift_type == "manual":
        manual_value = get_float_input("Enter custom shift value: ")
        args.extend(["--manual-value", str(manual_value)])

    q = get_int_input("IPR order parameter q (default: 2, must be >= 2): ", 2)
    args.extend(["--q", str(q)])

    output_dir = get_input("Output directory (default: '.'): ") or "."
    args.extend(["--output-dir", output_dir])

    save_report = get_input("\nAlso save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")
    save_gnuplot = get_input("Also save the data + gnuplot script? (y/N): ").strip().lower()
    if save_gnuplot == 'y':
        args.append("--save-gnuplot")
    view_choice = get_input("Show the matplotlib plot before finishing? (y/N): ").strip().lower()
    if view_choice == 'y':
        args.append("--view")

    run_tool("stb-ipr", args)


def run_effmass_analyzer() -> None:
    """Interface for the Effective Mass / Velocity Analyzer (stb-effmass)"""
    print("\n" + "="*60)
    print(color_text("EFFECTIVE MASS / VELOCITY (stb-effmass)", 'bold').center(60))
    print("="*60 + "\n")

    print(f"This tool needs a full-BZ {color_text('.WFSX', 'yellow')} file (found via the "
          f"SystemLabel) and a {color_text('.HSX', 'yellow')} file (mandatory -- no approximate "
          "fallback) -- asked separately below, since the real file is very often NOT named "
          "<label>.HSX.")
    label = get_input("Enter the Siesta SystemLabel (used to find <label>.WFSX etc.): ").strip()
    while not label:
        print(color_text("Label cannot be empty!", 'red'))
        label = get_input("Enter the Siesta SystemLabel: ")

    default_hsx = f"{label}.HSX"
    if os.path.isfile(default_hsx):
        hsx_file = get_input(f"Path to the .HSX file [default: {default_hsx}]: ").strip() or default_hsx
    else:
        print(color_text(f"-> No '{default_hsx}' found -- enter the real .HSX path below "
                          "(stb-effmass needs a real Hamiltonian, no approximate fallback).", 'yellow'))
        hsx_file = get_input("Path to the .HSX file: ").strip()
        while not hsx_file:
            print(color_text("An .HSX file is required!", 'red'))
            hsx_file = get_input("Path to the .HSX file: ").strip()

    args = ["--label", label, "--hsx-file", hsx_file, "--no-intro"]

    print(f"\n{color_text('Band selection:', 'yellow')}")
    print(f"  {color_text('1', 'cyan')} = Band index N at a chosen k-index")
    print(f"  {color_text('2', 'cyan')} = VBM (needs Fermi energy)")
    print(f"  {color_text('3', 'cyan')} = CBM (needs Fermi energy)")
    choice = get_input("Select (1-3) [default: 1]: ").strip()

    if choice in ('2', '3'):
        args.extend(["--band", "vbm" if choice == '2' else "cbm"])
        print(f"\n{color_text('Fermi energy source:', 'yellow')}")
        print(f"  {color_text('1', 'cyan')} = Enter a value directly")
        print(f"  {color_text('2', 'cyan')} = Read from a .bands file")
        print(f"  {color_text('3', 'cyan')} = Read from a SIESTA .out log (any filename, not just <label>.out)")
        print(f"  {color_text('4', 'cyan')} = Auto-detect a .out in this directory [default]")
        fermi_choice = get_input("Select (1-4) [default: 4]: ").strip() or '4'
        if fermi_choice == '1':
            fermi = get_float_input("Fermi energy (eV): ")
            args.extend(["--fermi", str(fermi)])
        elif fermi_choice == '2':
            bands_file = get_input("Path to the .bands file: ").strip()
            args.extend(["--bands-file", bands_file])
        elif fermi_choice == '3':
            out_file = get_input("Path to the .out log (any filename): ").strip()
            args.extend(["--fermi-file", out_file])
        # fermi_choice == '4': pass nothing, stb-effmass auto-detects a .out itself.
    else:
        k_index = get_int_input("k-index (0-based) [default: 0]: ", 0)
        args.extend(["--k-index", str(k_index)])
        band_n = get_int_input("Band number (1-based): ", 1)
        args.extend(["--band", str(band_n)])

    output_dir = get_input("Output directory (default: '.'): ") or "."
    args.extend(["--output-dir", output_dir])

    save_report = get_input("\nAlso save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")
    save_gnuplot = get_input("Also save the data + gnuplot script? (y/N): ").strip().lower()
    if save_gnuplot == 'y':
        args.append("--save-gnuplot")
    view_choice = get_input("Show the matplotlib plot before finishing? (y/N): ").strip().lower()
    if view_choice == 'y':
        args.append("--view")

    run_tool("stb-effmass", args)


def run_spintexture_analyzer() -> None:
    """Interface for the Spin Texture Analyzer (stb-spintexture)"""
    print("\n" + "="*60)
    print(color_text("SPIN TEXTURE ANALYZER (stb-spintexture)", 'bold').center(60))
    print("="*60 + "\n")

    print(f"This tool only makes sense for a {color_text('non-collinear/spin-orbit', 'yellow')} "
          "calculation (nspin=4 or 8).")
    label = get_input("Enter the Siesta SystemLabel (used to find <label>.WFSX etc.): ").strip()
    while not label:
        print(color_text("Label cannot be empty!", 'red'))
        label = get_input("Enter the Siesta SystemLabel: ")

    args = ["--label", label, "--no-intro"]

    default_hsx = f"{label}.HSX"
    if os.path.isfile(default_hsx):
        hsx_file = get_input(
            f"Path to the .HSX file for overlap-aware (accurate) spin moments "
            f"[default: {default_hsx}]: ").strip() or default_hsx
        args.extend(["--hsx-file", hsx_file])
    else:
        print(color_text(f"-> No '{default_hsx}' found -- spin moments will use an approximate "
                          "orthogonal-basis fallback unless you point at one below.", 'yellow'))
        hsx_file = get_input(
            "Path to the .HSX file (blank to use the approximate fallback): ").strip()
        if hsx_file:
            args.extend(["--hsx-file", hsx_file])

    shift_options = {
        '1': ('none', "Raw eigenvalues (no shift)"),
        '2': ('vbm', "Valence Band Maximum"),
        '3': ('cbm', "Conduction Band Minimum"),
        '4': ('fermi', "Fermi level"),
        '5': ('manual', "Custom value"),
    }
    print("\nEnergy reference options:")
    for key, (_, desc) in shift_options.items():
        print(f" {color_text(key, 'yellow')}. {desc}")
    choice = get_input("\nSelect reference (1-5) [default: 1]: ").strip() or '1'
    while choice not in shift_options:
        print(color_text("Invalid choice!", 'red'))
        choice = get_input("Select reference (1-5) [default: 1]: ").strip() or '1'

    shift_type, _ = shift_options[choice]
    args.extend(["--shift", shift_type])

    if shift_type == "manual":
        manual_value = get_float_input("Enter custom shift value: ")
        args.extend(["--manual-value", str(manual_value)])
    elif shift_type in ("fermi", "vbm", "cbm"):
        print(f"\n{color_text('Fermi energy source:', 'yellow')}")
        print(f"  {color_text('1', 'cyan')} = Enter a value directly")
        print(f"  {color_text('2', 'cyan')} = Read from a .bands file")
        print(f"  {color_text('3', 'cyan')} = Read from a SIESTA .out log (any filename, not just <label>.out)")
        print(f"  {color_text('4', 'cyan')} = Auto-detect a .out in this directory [default]")
        fermi_choice = get_input("Select (1-4) [default: 4]: ").strip() or '4'
        if fermi_choice == '1':
            fermi = get_float_input("Fermi energy (eV): ")
            args.extend(["--fermi", str(fermi)])
        elif fermi_choice == '2':
            bands_file = get_input("Path to the .bands file: ").strip()
            args.extend(["--bands-file", bands_file])
        elif fermi_choice == '3':
            out_file = get_input("Path to the .out log (any filename): ").strip()
            args.extend(["--fermi-file", out_file])
        # fermi_choice == '4': pass nothing, stb-spintexture auto-detects a .out itself.

    output_dir = get_input("Output directory (default: '.'): ") or "."
    args.extend(["--output-dir", output_dir])

    save_report = get_input("\nAlso save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")
    save_gnuplot = get_input("Also save the data + gnuplot script? (y/N): ").strip().lower()
    if save_gnuplot == 'y':
        args.append("--save-gnuplot")
    view_choice = get_input("Show the matplotlib plot before finishing? (y/N): ").strip().lower()
    if view_choice == 'y':
        args.append("--view")

    run_tool("stb-spintexture", args)


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

    save_report = get_input("\nAlso save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")

    save_gnuplot = get_input("Also save the data + gnuplot script? (y/N): ").strip().lower()
    if save_gnuplot == 'y':
        args.append("--save-gnuplot")

    view_choice = get_input("Show the matplotlib plot before finishing? (y/N): ").strip().lower()
    if view_choice == 'y':
        args.append("--view")

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
    output_dir = get_input("Output directory [default: current directory]: ").strip() or "."

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

    save_report = get_input("\nAlso save a text report to file? (y/N): ").strip().lower()

    # Argumentos básicos
    args = ["--label", label, "--speed", speed_mode, "--output-dir", output_dir, "--no-intro"]

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
    if save_report == 'y':
        args.append("--save-report")

    run_tool("stb-bader", args)
    

def _print_config_summary(title: str, rows: list) -> None:
    """Boxed recap of every EFFECTIVE setting (custom or default) right before dispatch --
    shared by run_elastic_generator/run_elastic_analyzer so the user has full visibility into
    what's about to run, even for values they never touched (no confirmation gate: purely
    informational, matches the "show but always proceed" choice for this workflow)."""
    print("\n" + "-"*60)
    print(color_text(title, 'cyan').center(60))
    print("-"*60)
    for label, value in rows:
        print(f"  {label:<26}: {color_text(str(value), 'green')}")
    print("-"*60)


def run_elastic_generator() -> None:
    """Interface for the Elastic Constants Generator (elastic_inputs.py)"""
    print("\n" + "="*60)
    print(color_text("ELASTIC CONSTANTS GENERATOR", 'bold').center(60))
    print("="*60)
    print(color_text(
        "Generates deformed structures for SIESTA calculations. By default it strains only "
        "the minimal symmetry-distinct set of directions -- stb-elasticAnalysis reconstructs "
        "the rest EXACTLY from the point group, so this is a real DFT-cost cut, not an "
        "approximation.", 'cyan'))
    print()

    # 1. Structure file
    input_file = get_input("Input structure.fdf (lattice + atomic coordinates): ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input structure.fdf: ")

    # 1b. calc.fdf template
    default_calc = "calc.fdf"
    if os.path.isfile(default_calc):
        calc_file = get_input(
            "Input calc.fdf -- SCF/basis/pseudopotential/MD settings (normally %include-ing "
            f"the structure file above) [default: {default_calc}]: ").strip() or default_calc
    else:
        print(color_text(f"-> No '{default_calc}' found -- enter the real calc.fdf path below.", 'yellow'))
        calc_file = get_input(
            "Input calc.fdf (SCF/basis/pseudopotential/MD settings, normally %include-ing the "
            "structure file): ").strip()
    while not os.path.isfile(calc_file):
        print(color_text("File not found!", 'red'))
        calc_file = get_input("Input calc.fdf: ").strip()

    # 1c. Pseudopotentials (optional)
    pseudo_dir = prompt_pseudo_source(optional=True)

    # 2. Sampling
    max_strain = get_float_input("\nMax strain % (default: 2.0): ", 2.0)
    steps = get_int_input("Number of steps per direction (default: 4): ", 4)

    # 3. Fitting method (stress-strain linear vs. energy-strain parabolic)
    print("\n" + "-"*60)
    print(color_text("FITTING METHOD", 'cyan').center(60))
    print("-"*60)
    print(f"[{color_text('1', 'yellow')}] Stress-strain (default) -> linear fit, sigma = C . epsilon")
    print(f"[{color_text('2', 'yellow')}] Energy-strain           -> parabolic fit (independent cross-check)")
    method_choice = get_input("Select method (1-2, default 1): ").strip()
    method = "energy" if method_choice == '2' else "stress"

    # 4. Direction/symmetry strategy
    symmetry_method = None  # only meaningful (and only ever read by elastic_inputs.py) for
                             # --dirs all + --method stress -- --method energy always uses its
                             # own point-group rank selection regardless of this flag's value,
                             # so it's never asked for that method.
    if method == "energy":
        print(color_text(
            "\n--method energy always auto-selects its own pure+combined strain patterns from "
            "the structure's point group -- no direction menu needed.", 'cyan'))
        direction_summary = "energy patterns, auto-selected by point group"
        dirs_args = ["--dirs", "all"]
    else:
        # All 6 canonical directions are always requested ("--dirs all") --
        # stb-elasticInputs auto-detects which ones touch a vacuum-padded
        # axis (skipped, unrelated to symmetry -- already based on whatever
        # dimensionality the structure turns out to have) and reconstructs
        # every symmetry-equivalent one exactly. This choice only controls
        # HOW that symmetry reduction is done -- there's nothing left for a
        # manual direction-subset menu to add.
        print("\n" + "-"*60)
        print(color_text("SELECT DEFORMATION MODE", 'cyan').center(60))
        print("-"*60)
        print(color_text(
            "All 6 canonical directions are requested; vacuum-padded axes are skipped "
            "automatically according to the structure's own dimensionality. This choice "
            "loses no information either way -- stb-elasticAnalysis reconstructs every "
            "skipped direction's stiffness column exactly.", 'cyan'))
        print(f"[{color_text('1', 'yellow')}] Basic (default) -> pairwise point-group match "
              "(misses hexagonal/trigonal-type reductions)")
        print(f"[{color_text('2', 'yellow')}] Full             -> symmetry-allowed-subspace rank "
              "selection (catches those too)")
        print("-" * 60)

        mode = get_input("Select mode (1-2, default 1): ").strip()
        dirs_args = ["--dirs", "all"]
        symmetry_method = "full" if mode == '2' else "basic"
        direction_summary = f"all 6, symmetry-reduced ({symmetry_method})"

    # 5. Advanced settings (rarely-touched tolerances -- gated so the essential
    # flow above stays short; CLI defaults apply untouched when skipped).
    vacuum_gap, symprec, angle_tolerance, output_dir = 10.0, 0.01, 5.0, "elastic_runs"
    show_advanced = get_input(
        "\nConfigure advanced settings (vacuum-gap, symmetry tolerances, output directory)? "
        "[y/N]: ").strip().lower()
    if show_advanced == 'y':
        vacuum_gap = get_float_input(
            "Vacuum gap threshold in Ang, for detecting periodic vs. vacuum-padded axes "
            "(default: 10.0): ", 10.0)
        symprec = get_float_input("Symmetry-detection tolerance (default: 0.01): ", 0.01)
        angle_tolerance = get_float_input("Symmetry angle tolerance in degrees (default: 5.0): ", 5.0)
        output_dir = get_input("Output directory (default: elastic_runs): ").strip()
        if not output_dir:
            output_dir = "elastic_runs"

    save_report = get_input("\nAlso save a report to file? [y/N]: ").strip().lower() == 'y'

    args = [
        "--structure", input_file,
        "--calc", calc_file,
        "--max", str(max_strain),
        "--steps", str(steps),
        "--method", method,
        "--vacuum-gap", str(vacuum_gap),
        "--symprec", str(symprec),
        "--angle-tolerance", str(angle_tolerance),
        "--output-dir", output_dir,
        "--no-intro",
    ] + dirs_args
    if pseudo_dir:
        args.extend(["-p", pseudo_dir])
    if symmetry_method is not None:
        args.extend(["--symmetry-method", symmetry_method])
    if save_report:
        args.append("--save-report")

    # 6. Recap -- always shown, no confirmation gate
    summary_rows = [
        ("Structure file", input_file),
        ("Calc template", calc_file),
        ("Pseudo source", pseudo_dir or "(not given)"),
        ("Strain range", f"+/-{max_strain}%, {steps} steps"),
        ("Fitting method", method),
        ("Direction strategy", direction_summary),
        ("Vacuum-gap threshold", f"{vacuum_gap} Ang"),
        ("Symmetry tolerance", f"symprec={symprec}, angle={angle_tolerance} deg"),
        ("Output directory", output_dir),
        ("Save report", "yes" if save_report else "no"),
    ]
    _print_config_summary("CONFIGURATION SUMMARY", summary_rows)

    run_tool("stb-elasticInputs", args)

def run_elastic_analyzer() -> None:
    """Interface for the Elastic Properties Analyzer """
    print("\n" + "="*60)
    print(color_text("ELASTIC PROPERTIES ANALYZER", 'bold').center(60))
    print("="*60)
    print(color_text(
        "Fits the stiffness matrix from your strain_*/ results. Numerical-quality "
        "diagnostics (eggbox cross-check, fit R^2, tensor symmetry) run automatically and "
        "never block the report -- they only flag suspicious data.", 'cyan'))
    print()

    # 1. Directory with the strain_* run folders (same default as
    # stb-elasticInputs' own --output-dir, so the two stages chain with no
    # extra typing in the common case).
    run_dir = get_input("Directory with elastic run folders [default: elastic_runs]: ").strip()
    if not run_dir:
        run_dir = "elastic_runs"

    # 2. Output/log filename
    output_filename = get_input("Siesta output filename inside strain folders (default: calc.out): ").strip()
    if not output_filename:
        output_filename = "calc.out"

    # 3. Fitting method
    print("\n" + "-"*60)
    print(color_text("FITTING METHOD", 'cyan').center(60))
    print("-"*60)
    print(f"[{color_text('1', 'yellow')}] Stress-strain (default) -> linear fit, sigma = C . epsilon")
    print(f"[{color_text('2', 'yellow')}] Energy-strain           -> parabolic fit (independent cross-check)")
    method_choice = get_input("Select method (1-2, default 1): ").strip()
    method = "energy" if method_choice == '2' else "stress"

    # 4. Reference (undeformed) structure -- feeds dimensionality auto-detection,
    # symmetry detection, and --method energy's reference volume/area/length.
    # Left blank by default -- the CLI itself now defaults to
    # <run_dir>/reference_structure.fdf.
    ref_structure = get_input(
        f"\nUndeformed reference structure file [default: {run_dir}/reference_structure.fdf]: ").strip()

    # 5. Dimensionality -- auto-detected by default, same convention as
    # stb-strainAnalysis' own interactive wrapper.
    dimensionality = get_input(
        "\nDimensionality, auto/3d/2d/1d [default: auto]: ").strip()

    # 6. Symmetry reduction (only meaningful for --method stress; energy always
    # uses the general symmetry-allowed-subspace fit regardless of this flag)
    symmetry_method = None
    if method == "stress":
        sym_choice = get_input(
            "\nSymmetry reduction -- [1] basic (default) or [2] full "
            "(also catches hexagonal/trigonal reductions)? : ").strip()
        symmetry_method = "full" if sym_choice == '2' else "basic"

    # 7. Advanced settings (rarely-touched tolerances -- gated so the essential
    # flow above stays short; CLI defaults apply untouched when skipped).
    vacuum_gap, symprec, angle_tolerance = 10.0, 0.01, 5.0
    eggbox_tolerance, fit_quality_tolerance, plot_dir = 5.0, 0.995, "elastic_plots"
    show_advanced = get_input(
        "\nConfigure advanced settings (vacuum-gap, symmetry/quality tolerances, plot "
        "directory)? [y/N]: ").strip().lower()
    if show_advanced == 'y':
        if not dimensionality or dimensionality == "auto":
            vacuum_gap = get_float_input(
                "Vacuum gap threshold in Ang, for dimensionality auto-detection "
                "(default: 10.0): ", 10.0)
        symprec = get_float_input("Symmetry-detection tolerance (default: 0.01): ", 0.01)
        angle_tolerance = get_float_input("Symmetry angle tolerance in degrees (default: 5.0): ", 5.0)
        if method == "stress":
            eggbox_tolerance = get_float_input(
                "Eggbox cross-check tolerance % (default: 5.0): ", 5.0)
            fit_quality_tolerance = get_float_input(
                "Fit-quality R^2 tolerance (default: 0.995): ", 0.995)
        plot_dir = get_input("Plot data directory (default: elastic_plots): ").strip()
        if not plot_dir:
            plot_dir = "elastic_plots"

    if dimensionality and dimensionality != "auto":
        dimensionality_summary = dimensionality.upper() + " (manual override)"
    else:
        dimensionality_summary = f"auto-detect (vacuum-gap={vacuum_gap} Ang)"

    # 8. Save report?
    save_report = get_input("\nAlso save a text report to file? (y/N): ").strip().lower() == 'y'

    # 9. Save the raw fitting data + gnuplot scripts, and/or preview on screen?
    save_gnuplot = get_input(
        "Also save the fitting data + a gnuplot script? (y/N): ").strip().lower() == 'y'
    view = get_input(
        "Show an interactive matplotlib plot before finishing? (y/N): ").strip().lower() == 'y'

    args = [
        "--dir", run_dir,
        "--file", output_filename,
        "--method", method,
        "--vacuum-gap", str(vacuum_gap),
        "--symprec", str(symprec),
        "--angle-tolerance", str(angle_tolerance),
        "--plot-dir", plot_dir,
        "--no-intro",
    ]
    if ref_structure:
        args.extend(["--reference-structure", ref_structure])
    if dimensionality:
        args.extend(["--dimensionality", dimensionality])
    if symmetry_method is not None:
        args.extend(["--symmetry-method", symmetry_method])
    if method == "stress":
        args.extend(["--eggbox-tolerance", str(eggbox_tolerance),
                     "--fit-quality-tolerance", str(fit_quality_tolerance)])
    if save_report:
        args.append("--save-report")
    if save_gnuplot:
        args.append("--save-gnuplot")
    if view:
        args.append("--view")

    # 10. Recap -- always shown, no confirmation gate
    summary_rows = [
        ("Run directory", run_dir),
        ("Log filename", output_filename),
        ("Fitting method", method),
        ("Reference structure", ref_structure or f"{run_dir}/reference_structure.fdf (default)"),
        ("Dimensionality", dimensionality_summary),
    ]
    if symmetry_method is not None:
        summary_rows.append(("Symmetry reduction", symmetry_method))
    summary_rows.append(("Symmetry tolerance", f"symprec={symprec}, angle={angle_tolerance} deg"))
    if method == "stress":
        summary_rows.append(("Quality diagnostics", f"eggbox={eggbox_tolerance}%, "
                              f"fit-quality R^2={fit_quality_tolerance}"))
    summary_rows.append(("Plot data directory", plot_dir))
    summary_rows.append(("Save report", "yes" if save_report else "no"))
    summary_rows.append(("Save gnuplot", "yes" if save_gnuplot else "no"))
    summary_rows.append(("View (matplotlib)", "yes" if view else "no"))
    _print_config_summary("CONFIGURATION SUMMARY", summary_rows)

    print(color_text(f"\nRunning analysis in '{run_dir}'...", 'yellow'))
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
    mode_list = ['total_energy', 'relax', 'aimd', 'bands']

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

    # --- PONTO 1b: DFT-D3, independente do modo de cálculo ---
    d3_choice = get_input("Enable the DFT-D3 dispersion correction? (y/N): ").strip().lower()

    # --- PONTO 1c: Spin polarization, independente do modo de cálculo ---
    spin_choice = get_input("Spin-polarized calculation? (y/N, default: non-polarized): ").strip().lower()

    # --- PONTO 2: Validação do caminho do Pseudopotencial ---

    args = [
        input_file,
        "--type", calc_type,
        "--no-intro"
    ]

    if d3_choice == 'y':
        args.append("--d3")

    if spin_choice == 'y':
        args.append("--spin-polarized")

    pp_path = prompt_pseudo_source(optional=True)
    if pp_path:
        args.extend(["--pp-path", pp_path])
    else:
        print(color_text("Skipping pseudopotential path.", 'yellow'))

    save_report = get_input("Also save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")

    view_choice = get_input("View the structure interactively via ASE before finishing? (y/N): ").strip().lower()
    if view_choice == 'y':
        args.append("--view")

    run_tool("stb-inputfile", args)

def run_kgrid_generator() -> None:
    """Interface for the K-Grid Generator (stb-kgrid)"""
    print("\n" + "="*60)
    print(color_text("K-GRID GENERATOR (stb-kgrid)", 'bold').center(60))
    print("="*60 + "\n")

    # Esta linha agora terá Tab-completion!
    input_file = get_input("Input structure file (.fdf): ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input structure file: ")

    # Show the density recommendation guide BEFORE asking, so the choice is informed.
    kspace.print_density_recommendation()

    density = get_float_input("\nK-point density (e.g., 0.2) [default: 0.2]: ", 0.2)
    while density <= 0:
        print(color_text("Density must be a positive number!", 'red'))
        density = get_float_input("K-point density (e.g., 0.2) [default: 0.2]: ", 0.2)

    args = [
        "--file", input_file,
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

    # Advanced settings (rarely-touched -- gated so the essential flow above
    # stays short; CLI defaults apply untouched when skipped).
    vacuum_gap, eps, symprec, angle_tolerance = 10.0, 0.0002, 0.01, 5.0
    show_advanced = get_input(
        "\nConfigure advanced settings (vacuum-gap, Bravais-lattice tolerance, "
        "symmetry tolerances)? [y/N]: ").strip().lower()
    if show_advanced == 'y':
        vacuum_gap = get_float_input(
            "Vacuum gap threshold in Ang, for detecting periodic vs. vacuum-padded axes "
            "(default: 10.0): ", 10.0)
        eps = get_float_input("Bravais-lattice detection tolerance / eps (default: 0.0002): ", 0.0002)
        symprec = get_float_input("Symmetry-detection tolerance (default: 0.01): ", 0.01)
        angle_tolerance = get_float_input("Symmetry angle tolerance in degrees (default: 5.0): ", 5.0)

    args = [
        "--file", input_file,
        "--vacuum-gap", str(vacuum_gap),
        "--prec", str(eps),
        "--symprec", str(symprec),
        "--angle-tolerance", str(angle_tolerance),
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

    ml_relax = get_input(
        "\nPre-relax the built supercell with a MACE ML potential before writing? "
        "(needs the optional 'ml' extra) (y/N): "
    ).strip().lower()
    if ml_relax == 'y':
        args.append("--ml-relax")

        ml_relax_cell = get_input(
            "Also relax the cell? Any vacuum-padded axis always stays fixed (y/N): "
        ).strip().lower()
        if ml_relax_cell == 'y':
            args.append("--ml-relax-cell")

        custom_model = get_input(
            "Use a custom MACE model instead (e.g. one fine-tuned via stb-mlffAnalysis)? "
            "Path, or Enter to use a MACE-MP-0 foundation model: ").strip()
        if custom_model:
            while not os.path.isfile(custom_model):
                print(color_text("File not found!", 'red'))
                custom_model = get_input("Custom model path: ").strip()
            args.extend(["--custom-model", custom_model])
        else:
            model = get_input("Model size, small/medium/large [default: small]: ").strip()
            if not model:
                model = "small"
            args.extend(["--model", model])

    save_report = get_input("Also save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")

    view_choice = get_input("View the structure interactively via ASE before finishing? (y/N): ").strip().lower()
    if view_choice == 'y':
        args.append("--view")

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

    ml_relax = get_input(
        "\nPre-relax each written slab with a MACE ML potential before writing? "
        "(needs the optional 'ml' extra) (y/N): "
    ).strip().lower()
    if ml_relax == 'y':
        args.append("--ml-relax")

        ml_relax_cell = get_input(
            "Also relax the in-plane cell? The vacuum axis always stays fixed (y/N): "
        ).strip().lower()
        if ml_relax_cell == 'y':
            args.append("--ml-relax-cell")

        custom_model = get_input(
            "Use a custom MACE model instead (e.g. one fine-tuned via stb-mlffAnalysis)? "
            "Path, or Enter to use a MACE-MP-0 foundation model: ").strip()
        if custom_model:
            while not os.path.isfile(custom_model):
                print(color_text("File not found!", 'red'))
                custom_model = get_input("Custom model path: ").strip()
            args.extend(["--custom-model", custom_model])
        else:
            model = get_input("Model size, small/medium/large [default: small]: ").strip()
            if not model:
                model = "small"
            args.extend(["--model", model])

    save_report = get_input("Also save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")

    view_choice = get_input("View the structure interactively via ASE before finishing? (y/N): ").strip().lower()
    if view_choice == 'y':
        args.append("--view")

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

    passivate_choice = get_input(
        "\nPassivate dangling bonds on the built tube/ribbon? A closed tube has none; "
        "a ribbon's two edges do (y/N): "
    ).strip().lower()
    passivant = None
    cutoff_str = ""
    bond_length_str = ""
    if passivate_choice in ('y', 'yes'):
        passivant = get_input("Passivant element [default: H]: ").strip() or "H"
        cutoff_str = get_input("Neighbor cutoff, Ang (blank = auto-detect): ").strip()
        bond_length_str = get_input("Bond length, Ang (blank = auto per species pair): ").strip()

    default_output = "nanoribbon.fdf" if mode == "ribbon" else "nanotube.fdf"
    output_file = get_input(f"\nOutput file name [default: {default_output}]: ").strip()
    if not output_file:
        output_file = default_output

    args = [
        "--file", input_file,
        "--chirality", *[str(v) for v in chirality],
        "--mode", mode,
        "--repeats", str(repeats),
        "--min-vacuum-size", str(min_vacuum_size),
        "--output", output_file,
        "--no-intro"
    ]

    if passivant is not None:
        args.extend(["--passivate", "--passivant", passivant])
        if cutoff_str:
            args.extend(["--cutoff", cutoff_str])
        if bond_length_str:
            args.extend(["--bond-length", bond_length_str])

    ml_relax = get_input(
        "\nPre-relax the built tube/ribbon with a MACE ML potential before writing? "
        "(needs the optional 'ml' extra) (y/N): "
    ).strip().lower()
    if ml_relax == 'y':
        args.append("--ml-relax")

        ml_relax_cell = get_input(
            "Also relax the axial period? The vacuum axes always stay fixed (y/N): "
        ).strip().lower()
        if ml_relax_cell == 'y':
            args.append("--ml-relax-cell")

        custom_model = get_input(
            "Use a custom MACE model instead (e.g. one fine-tuned via stb-mlffAnalysis)? "
            "Path, or Enter to use a MACE-MP-0 foundation model: ").strip()
        if custom_model:
            while not os.path.isfile(custom_model):
                print(color_text("File not found!", 'red'))
                custom_model = get_input("Custom model path: ").strip()
            args.extend(["--custom-model", custom_model])
        else:
            model = get_input("Model size, small/medium/large [default: small]: ").strip()
            if not model:
                model = "small"
            args.extend(["--model", model])

    save_report = get_input("Also save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")

    view_choice = get_input("View the structure interactively via ASE before finishing? (y/N): ").strip().lower()
    if view_choice == 'y':
        args.append("--view")

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

    if "--ml-rank" not in args:
        ml_relax = get_input(
            "\nPre-relax the output structure(s) with a MACE ML potential before writing? "
            "(needs the optional 'ml' extra) (y/N): "
        ).strip().lower()
        if ml_relax == 'y':
            args.append("--ml-relax")

            ml_relax_cell = get_input(
                "Also relax the cell? Any vacuum-padded axis always stays fixed (y/N): "
            ).strip().lower()
            if ml_relax_cell == 'y':
                args.append("--ml-relax-cell")

            custom_model = get_input(
                "Use a custom MACE model instead (e.g. one fine-tuned via stb-mlffAnalysis)? "
                "Path, or Enter to use a MACE-MP-0 foundation model: ").strip()
            if custom_model:
                while not os.path.isfile(custom_model):
                    print(color_text("File not found!", 'red'))
                    custom_model = get_input("Custom model path: ").strip()
                args.extend(["--custom-model", custom_model])
            else:
                model = get_input("Model size, small/medium/large [default: small]: ").strip()
                if not model:
                    model = "small"
                args.extend(["--model", model])

    save_report = get_input("Also save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")

    view_choice = get_input("View the structure interactively via ASE before finishing? (y/N): ").strip().lower()
    if view_choice == 'y':
        args.append("--view")

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

    ml_relax = get_input(
        "\nPre-relax the SQS structure found with a MACE ML potential before writing? "
        "(needs the optional 'ml' extra) (y/N): "
    ).strip().lower()
    if ml_relax == 'y':
        args.append("--ml-relax")

        ml_relax_cell = get_input(
            "Also relax the cell? Any vacuum-padded axis always stays fixed (y/N): "
        ).strip().lower()
        if ml_relax_cell == 'y':
            args.append("--ml-relax-cell")

        custom_model = get_input(
            "Use a custom MACE model instead (e.g. one fine-tuned via stb-mlffAnalysis)? "
            "Path, or Enter to use a MACE-MP-0 foundation model: ").strip()
        if custom_model:
            while not os.path.isfile(custom_model):
                print(color_text("File not found!", 'red'))
                custom_model = get_input("Custom model path: ").strip()
            args.extend(["--custom-model", custom_model])
        else:
            model = get_input("Model size, small/medium/large [default: small]: ").strip()
            if not model:
                model = "small"
            args.extend(["--model", model])

    save_report = get_input("Also save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")

    view_choice = get_input("View the structure interactively via ASE before finishing? (y/N): ").strip().lower()
    if view_choice == 'y':
        args.append("--view")

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

    symprec = get_float_input("\nSymmetry precision [default: 0.01]: ", 0.01)

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

    ml_relax = get_input(
        "\nPre-relax the reduced unit cell with a MACE ML potential before writing? "
        "(needs the optional 'ml' extra) (y/N): "
    ).strip().lower()
    if ml_relax == 'y':
        args.append("--ml-relax")

        ml_relax_cell = get_input(
            "Also relax the cell? Any vacuum-padded axis always stays fixed (y/N): "
        ).strip().lower()
        if ml_relax_cell == 'y':
            args.append("--ml-relax-cell")

        custom_model = get_input(
            "Use a custom MACE model instead (e.g. one fine-tuned via stb-mlffAnalysis)? "
            "Path, or Enter to use a MACE-MP-0 foundation model: ").strip()
        if custom_model:
            while not os.path.isfile(custom_model):
                print(color_text("File not found!", 'red'))
                custom_model = get_input("Custom model path: ").strip()
            args.extend(["--custom-model", custom_model])
        else:
            model = get_input("Model size, small/medium/large [default: small]: ").strip()
            if not model:
                model = "small"
            args.extend(["--model", model])

    save_report = get_input("Also save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")

    view_choice = get_input("View the structure interactively via ASE before finishing? (y/N): ").strip().lower()
    if view_choice == 'y':
        args.append("--view")

    run_tool("stb-unitcell", args)


def run_crystalbuilder_generator() -> None:
    """Interface for the Crystal Builder (stb-crystalbuilder)"""
    print("\n" + "="*60)
    print(color_text("CRYSTAL BUILDER (stb-crystalbuilder)", 'bold').center(60))
    print("="*60 + "\n")
    print(color_text(
        "Bulk (3D periodic) structures only -- the 230 space groups accepted "
        "below are the international classification for 3D crystals, with no "
        "2D/1D analog. For a slab, monolayer, nanotube, or wire: build the 3D "
        "bulk crystal here first, then cut it down with stb-slab/stb-nanotube/"
        "stb-2Dstacking.", 'yellow'))

    spacegroup = get_input("\nSpace group (symbol e.g. 'Fm-3m' or number e.g. '225'): ").strip()

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

    print(f"\n{color_text('Reduce to a smaller cell?', 'yellow')}")
    print(f"  {color_text('1', 'cyan')} = No reduction [Default]")
    print(f"  {color_text('2', 'cyan')} = Primitive (smallest cell)")
    print(f"  {color_text('3', 'cyan')} = Conventional (standardized, usually larger)")
    print(f"  {color_text('4', 'cyan')} = Refined (conventional cell, positions snapped to symmetry)")
    reduce_choice = get_input("Select option (1-4) [default: 1]: ").strip()
    reduce_mode = {'2': 'primitive', '3': 'conventional', '4': 'refined'}.get(reduce_choice)
    if reduce_mode:
        args.extend(["--reduce", reduce_mode])

    ml_relax = get_input(
        "\nPre-relax the final structure with a MACE ML potential before writing? "
        "(needs the optional 'ml' extra) (y/N): "
    ).strip().lower()
    if ml_relax == 'y':
        args.append("--ml-relax")

        ml_relax_cell = get_input(
            "Also relax the cell? (y/N): "
        ).strip().lower()
        if ml_relax_cell == 'y':
            args.append("--ml-relax-cell")

        custom_model = get_input(
            "Use a custom MACE model instead (e.g. one fine-tuned via stb-mlffAnalysis)? "
            "Path, or Enter to use a MACE-MP-0 foundation model: ").strip()
        if custom_model:
            while not os.path.isfile(custom_model):
                print(color_text("File not found!", 'red'))
                custom_model = get_input("Custom model path: ").strip()
            args.extend(["--custom-model", custom_model])
        else:
            model = get_input("Model size, small/medium/large [default: small]: ").strip()
            if not model:
                model = "small"
            args.extend(["--model", model])

    save_report = get_input("\nAlso save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")

    view_choice = get_input("View the structure interactively via ASE before finishing? (y/N): ").strip().lower()
    if view_choice == 'y':
        args.append("--view")

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
        "Mode: generate / substitute / subgroup / supergroup "
        "[default: generate]: ").strip().lower()

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
        args = ["--substitute", *subs, "-f", input_file]
        save_report = get_input("\nAlso save a text report to file? (y/N): ").strip().lower()
        if save_report == 'y':
            args.append("--save-report")
        view_choice = get_input("View the input and final structure interactively via ASE before finishing? (y/N): ").strip().lower()
        if view_choice == 'y':
            args.append("--view")
        output_file = get_input("Output file name [default: crystalcast.fdf]: ").strip() or "crystalcast.fdf"
        args.extend(["-o", output_file, "--no-intro"])
        run_tool("stb-crystalcast", args)
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
        save_report = get_input("\nAlso save a text report to file? (y/N): ").strip().lower()
        if save_report == 'y':
            args.append("--save-report")
        view_choice = get_input("View the input and final structure(s) interactively via ASE before finishing? (y/N): ").strip().lower()
        if view_choice == 'y':
            args.append("--view")
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
        "Rank candidates by MACE relaxed energy? Needs the optional 'ml' extra (y/n) "
        "[default: n]: ").strip().lower()
    if ml_rank_str == "y":
        args.append("--ml-rank")
        custom_model = get_input(
            "Use a custom MACE model instead (e.g. one fine-tuned via stb-mlffAnalysis)? "
            "Path, or Enter to use a MACE-MP-0 foundation model: ").strip()
        if custom_model:
            while not os.path.isfile(custom_model):
                print(color_text("File not found!", 'red'))
                custom_model = get_input("Custom model path: ").strip()
            args.extend(["--custom-model", custom_model])
        else:
            model = get_input("Model size, small/medium/large [default: small]: ").strip()
            if not model:
                model = "small"
            args.extend(["--model", model])

    save_report = get_input("\nAlso save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")

    view_choice = get_input("View the final structure(s) interactively via ASE before finishing? (y/N): ").strip().lower()
    if view_choice == 'y':
        args.append("--view")

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
    args.extend(["--output", output_file])

    save_report = get_input("Also save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")

    view_choice = get_input("View the structure interactively via ASE before finishing? (y/N): ").strip().lower()
    if view_choice == 'y':
        args.append("--view")

    args.append("--no-intro")

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

    args = ["--file", input_file, "--passivant", passivant, "--no-intro"]
    if cutoff_str:
        args.extend(["--cutoff", cutoff_str])
    if bond_length_str:
        args.extend(["--bond-length", bond_length_str])

    ml_relax = get_input(
        "\nPre-relax the passivated structure with a MACE ML potential before writing? "
        "(needs the optional 'ml' extra) (y/N): "
    ).strip().lower()
    if ml_relax == 'y':
        args.append("--ml-relax")

        ml_relax_cell = get_input(
            "Also relax the cell? Any vacuum-padded axis always stays fixed (y/N): "
        ).strip().lower()
        if ml_relax_cell == 'y':
            args.append("--ml-relax-cell")

        custom_model = get_input(
            "Use a custom MACE model instead (e.g. one fine-tuned via stb-mlffAnalysis)? "
            "Path, or Enter to use a MACE-MP-0 foundation model: ").strip()
        if custom_model:
            while not os.path.isfile(custom_model):
                print(color_text("File not found!", 'red'))
                custom_model = get_input("Custom model path: ").strip()
            args.extend(["--custom-model", custom_model])
        else:
            model = get_input("Model size, small/medium/large [default: small]: ").strip()
            if not model:
                model = "small"
            args.extend(["--model", model])

    save_report = get_input("\nAlso save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")

    view_choice = get_input("View the structure interactively via ASE before finishing? (y/N): ").strip().lower()
    if view_choice == 'y':
        args.append("--view")

    output_file = get_input("\nOutput file name [default: passivated.fdf]: ").strip()
    if not output_file:
        output_file = "passivated.fdf"
    args.extend(["--output", output_file])

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

    args = ["--name", name, "--vacuum", str(vacuum), "--no-intro"]

    ml_relax = get_input(
        "\nPre-relax the molecule's geometry with a MACE ML potential before writing? "
        "(needs the optional 'ml' extra) (y/N): "
    ).strip().lower()
    if ml_relax == 'y':
        args.append("--ml-relax")

        custom_model = get_input(
            "Use a custom MACE model instead (e.g. one fine-tuned via stb-mlffAnalysis)? "
            "Path, or Enter to use a MACE-MP-0 foundation model: ").strip()
        if custom_model:
            while not os.path.isfile(custom_model):
                print(color_text("File not found!", 'red'))
                custom_model = get_input("Custom model path: ").strip()
            args.extend(["--custom-model", custom_model])
        else:
            model = get_input("Model size, small/medium/large [default: small]: ").strip()
            if not model:
                model = "small"
            args.extend(["--model", model])

    save_report = get_input("\nAlso save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")

    view_choice = get_input("View the molecule interactively via ASE before finishing? (y/N): ").strip().lower()
    if view_choice == 'y':
        args.append("--view")

    output_file = get_input("\nOutput file name [default: molecule.fdf]: ").strip()
    if not output_file:
        output_file = "molecule.fdf"
    args.extend(["--output", output_file])

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

    custom_model = get_input(
        "\nUse a custom MACE model instead (e.g. one fine-tuned via stb-mlffAnalysis)? "
        "Path, or Enter to use a MACE-MP-0 foundation model: ").strip()

    args = ["--file", input_file, "--no-intro"]
    if custom_model:
        while not os.path.isfile(custom_model):
            print(color_text("File not found!", 'red'))
            custom_model = get_input("Custom model path: ").strip()
        args.extend(["--custom-model", custom_model])
    else:
        model = get_input("Model size, small/medium/large [default: small]: ").strip()
        if not model:
            model = "small"
        args.extend(["--model", model])

    fmax = get_float_input("Force convergence, eV/Ang [default: 0.05]: ", 0.05)

    optimizer = get_input("Optimizer, FIRE/BFGS/LBFGS [default: FIRE]: ").strip()
    if not optimizer:
        optimizer = "FIRE"

    output_file = get_input("\nOutput file name [default: relaxed.fdf]: ").strip()
    if not output_file:
        output_file = "relaxed.fdf"

    args.extend(["--fmax", str(fmax), "--optimizer", optimizer, "--output", output_file])
    if relax_cell in ('y', 'yes'):
        args.append("--relax-cell")

    save_report = get_input("Also save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")

    view_choice = get_input("View the structure interactively via ASE before finishing? (y/N): ").strip().lower()
    if view_choice == 'y':
        args.append("--view")

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

    save_data = get_input(
        "\nSave simulation data (.dat + gnuplot) for energy/temperature vs. time? (y/N): "
    ).strip().lower()
    save_traj = get_input(
        "Save a trajectory to view the atoms in OVITO/VMD? (y/N): "
    ).strip().lower()

    args = [
        "--file", input_file,
        "--melt-temp", str(melt_temp),
        "--melt-steps", str(melt_steps),
        "--quench-temp", str(quench_temp),
        "--quench-steps", str(quench_steps),
        "--no-intro"
    ]

    if save_data == 'y':
        args.append("--save-data")
    if save_traj == 'y':
        args.append("--save-traj")
        traj_format = get_input("Trajectory format, xsf/pdb/xyz [default: xsf]: ").strip()
        if not traj_format:
            traj_format = "xsf"
        args.extend(["--traj-format", traj_format])
    if save_data == 'y' or save_traj == 'y':
        stride = get_int_input("Sample every N MD steps [default: 10]: ", 10)
        args.extend(["--stride", str(stride)])

    final_relax = get_input(
        "\nDo a final static relax after the quench? (Y/n): "
    ).strip().lower()
    if final_relax in ('n', 'no'):
        args.append("--no-final-relax")

    custom_model = get_input(
        "\nUse a custom MACE model instead (e.g. one fine-tuned via stb-mlffAnalysis)? "
        "Path, or Enter to use a MACE-MP-0 foundation model: ").strip()
    if custom_model:
        while not os.path.isfile(custom_model):
            print(color_text("File not found!", 'red'))
            custom_model = get_input("Custom model path: ").strip()
        args.extend(["--custom-model", custom_model])
    else:
        model = get_input("Model size, small/medium/large [default: small]: ").strip()
        if not model:
            model = "small"
        args.extend(["--model", model])

    save_report = get_input("\nAlso save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")

    view_choice = get_input("View the input and final structure interactively via ASE before finishing? (y/N): ").strip().lower()
    if view_choice == 'y':
        args.append("--view")

    output_file = get_input("\nOutput file name [default: amorphous.fdf]: ").strip()
    if not output_file:
        output_file = "amorphous.fdf"
    args.extend(["--output", output_file])

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

    curve_output = get_input("Output curve data file name [default: convergence_curve.dat]: ").strip()
    if not curve_output:
        curve_output = "convergence_curve.dat"

    args = [
        "--dir", run_dir,
        "--file", output_filename,
        "--tolerance", str(tolerance),
        "--output", curve_output,
        "--no-intro"
    ]

    apply_target = get_input("\nApply the converged value to a calc.fdf now? "
                              "Path to the file, or leave blank to skip: ").strip()
    if apply_target:
        args.extend(["--apply", apply_target])
        vacuum_gap = get_float_input(
            "Vacuum-axis detection threshold in Ang, only used if the sweep "
            "was kgrid [default: 10.0]: ", 10.0)
        args.extend(["--vacuum-gap", str(vacuum_gap)])

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


def run_raman_prep() -> None:
    """Interface for the Raman Spectrum Stage 1 (stb-raman)"""
    print("\n" + "="*60)
    print(color_text("RAMAN SPECTRUM - STAGE 1: PHONON DISPLACEMENTS", 'bold').center(60))
    print("="*60 + "\n")

    structure_file = get_input("Input structure file [default: structure.fdf]: ").strip()
    if not structure_file:
        structure_file = "structure.fdf"

    calc_file = get_input("Calculation parameters file [default: calc.fdf]: ").strip()
    if not calc_file:
        calc_file = "calc.fdf"

    dim_input = get_input("\nSupercell dimensions for the phonon calc (e.g. '2 2 2') "
                           "[default: 2 2 2]: ").strip()
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

    distance = get_float_input("\nPhonon displacement distance in Ang [default: 0.02]: ", 0.02)

    pseudo_dir = prompt_pseudo_source(optional=True)
    if not pseudo_dir:
        pseudo_dir = "."

    output_dir = get_input("\nOutput root directory [default: raman_study]: ").strip()
    if not output_dir:
        output_dir = "raman_study"

    vacuum_gap = 10.0
    show_advanced = get_input("\nConfigure advanced settings (vacuum-gap)? [y/N]: ").strip().lower()
    if show_advanced == 'y':
        vacuum_gap = get_float_input(
            "Vacuum gap threshold in Ang, for the supercell-dimension advisory "
            "(default: 10.0): ", 10.0)

    args = [
        "-s", structure_file, "-c", calc_file,
        "-dim", str(dim_x), str(dim_y), str(dim_z),
        "-d", str(distance), "-p", pseudo_dir,
        "--vacuum-gap", str(vacuum_gap),
        "-O", output_dir, "--no-intro"
    ]

    print(color_text("\nGenerating Raman workflow phonon displacement folders...", 'green'))
    run_tool("stb-raman", args)


def run_raman_modes() -> None:
    """Interface for the Raman Spectrum Stage 2 (stb-ramanModes)"""
    print("\n" + "="*60)
    print(color_text("RAMAN SPECTRUM - STAGE 2: MODES & OPTICAL DISPLACEMENTS", 'bold').center(60))
    print("="*60 + "\n")

    run_dir = get_input("Directory written by Stage 1 [default: raman_study]: ").strip()
    if not run_dir:
        run_dir = "raman_study"

    calc_file = get_input("Calc.fdf template for the Optical calculations: ")
    while not os.path.isfile(calc_file):
        print(color_text("File not found!", 'red'))
        calc_file = get_input("Calc.fdf template for the Optical calculations: ")

    args = ["--directory", run_dir, "--calc", calc_file, "--no-intro"]

    modes_str = get_input(
        "\nMode indices to process, space-separated [optional, default: every "
        "non-acoustic mode]: ").strip()
    if modes_str:
        args.extend(["--modes"] + modes_str.split())

    full_tensor_choice = get_input(
        "\nCompute the FULL Raman tensor (Rxx,Ryy,Rzz,Rxy,Rxz,Ryz) instead of just the "
        "diagonal? Doubles the folder count per mode (12 instead of 6) (y/N): ").strip().lower()
    if full_tensor_choice in ('y', 'yes'):
        args.append("--full-tensor")

    if not modes_str:
        use_symmetry_choice = get_input(
            "\nUse group theory to skip modes that are symmetry-forbidden from being "
            "Raman-active (saves SIESTA time, needs a primitive cell) (y/N): ").strip().lower()
        if use_symmetry_choice in ('y', 'yes'):
            args.append("--use-symmetry")

    reduce_components_choice = get_input(
        "\nProbe only ONE direction (instead of all of them) for modes whose Raman tensor "
        "shape is already fixed by symmetry -- saves SIESTA time, needs a primitive cell "
        "(y/N): ").strip().lower()
    if reduce_components_choice in ('y', 'yes'):
        args.append("--reduce-components")

    skip_degenerate_choice = get_input(
        "\nFor degenerate modes (2+ bands sharing one irrep, e.g. E/T representations), "
        "compute only ONE representative and reuse its Raman activity/depolarization ratio "
        "for the rest (both are rotational invariants, identical for every degenerate "
        "partner) -- saves SIESTA time, needs a primitive cell (y/N): ").strip().lower()
    if skip_degenerate_choice in ('y', 'yes'):
        args.append("--skip-degenerate")

    displacement = get_float_input("\nFinite-difference displacement in Ang [default: 0.02]: ", 0.02)
    args.extend(["--displacement", str(displacement)])

    animations_choice = get_input(
        "\nExport a looping animation (.axsf, viewable in XCrySDen/VESTA) of each selected "
        "mode's eigendisplacement (y/N): ").strip().lower()
    if animations_choice in ('y', 'yes'):
        args.append("--export-animations")
        animation_frames = get_int_input("Frames per mode animation [default: 20]: ", 20)
        args.extend(["--animation-frames", str(animation_frames)])

    optical_mesh, optical_broaden = [10, 10, 10], 0.2
    show_advanced = get_input(
        "\nConfigure advanced settings (Optical mesh/broadening/frequency range/pseudopotential "
        "override/vacuum-gap/rotational-mode tolerance)? [y/N]: ").strip().lower()
    if show_advanced == 'y':
        mesh_input = get_input("Optical.Mesh k-grid (e.g. '10 10 10') [default: 10 10 10]: ").strip()
        if mesh_input:
            try:
                mesh_vals = [int(x) for x in mesh_input.split()]
                if len(mesh_vals) == 3:
                    optical_mesh = mesh_vals
            except ValueError:
                pass
        optical_broaden = get_float_input("Optical.Broaden in eV [default: 0.2]: ", 0.2)
        freq_min = get_input("Skip modes below this frequency in THz [optional]: ").strip()
        if freq_min:
            args.extend(["--freq-min", freq_min])
        freq_max = get_input("Skip modes above this frequency in THz [optional]: ").strip()
        if freq_max:
            args.extend(["--freq-max", freq_max])
        print(color_text(
            "\nPseudopotential source (blank = reuse whatever Stage 1 already copied into its "
            "phonon_disp/disp-*/ folders -- the normal case):", 'yellow'))
        pseudo_dir = prompt_pseudo_source(optional=True)
        if pseudo_dir:
            args.extend(["--pseudo-dir", pseudo_dir])
        vacuum_gap = get_float_input(
            "\nVacuum gap threshold in Ang, to detect an isolated (0D) molecule -- it has 3 "
            "extra trivial (free-rotation) modes at Gamma that don't count as real vibrations "
            "(default: 10.0): ", 10.0)
        args.extend(["--vacuum-gap", str(vacuum_gap)])
        rotational_mode_tol = get_float_input(
            "0D only: frequency threshold in THz below which a mode is treated as free rotation "
            "rather than a real vibration -- lower it if your molecule has genuine low-frequency "
            "skeletal/torsional modes (default: 2.0): ", 2.0)
        args.extend(["--rotational-mode-tol", str(rotational_mode_tol)])

    args.extend(["--optical-mesh", str(optical_mesh[0]), str(optical_mesh[1]), str(optical_mesh[2])])
    args.extend(["--optical-broaden", str(optical_broaden)])

    print(color_text("\nBuilding phonon force constants and Optical displacement folders...", 'green'))
    run_tool("stb-ramanModes", args)


def run_raman_analysis() -> None:
    """Interface for the Raman Spectrum Stage 3 (stb-ramanAnalysis)"""
    print("\n" + "="*60)
    print(color_text("RAMAN SPECTRUM - STAGE 3: ANALYSIS", 'bold').center(60))
    print("="*60 + "\n")

    run_dir = get_input("Directory with the stb-raman/stb-ramanModes run [default: raman_study]: ").strip()
    if not run_dir:
        run_dir = "raman_study"

    output_filename = get_input("SIESTA output filename inside each mode_*/ folder "
                                 "[default: calc.out]: ").strip()
    if not output_filename:
        output_filename = "calc.out"

    args = ["--directory", run_dir, "--file", output_filename, "--no-intro"]

    linewidth = get_float_input("\nLorentzian linewidth for the spectrum, cm^-1 [default: 10.0]: ", 10.0)
    args.extend(["--linewidth", str(linewidth)])

    temperature = get_input(
        "\nWeight the spectrum by the Bose-Einstein thermal population factor at a given "
        "temperature in K [optional, blank = no weighting]: ").strip()
    if temperature:
        args.extend(["--temperature", temperature])

    experimental_file = get_input(
        "\nOverlay a measured Raman spectrum (2-column text file: frequency cm^-1, intensity) "
        "[optional, blank = skip]: ").strip()
    if experimental_file:
        args.extend(["--experimental", experimental_file])
        peak_prominence = get_float_input(
            "Minimum peak prominence, as a fraction of each spectrum's own max intensity "
            "[default: 0.01]: ", 0.01)
        args.extend(["--peak-prominence", str(peak_prominence)])

    run_tool("stb-ramanAnalysis", args)


def run_ir_prep() -> None:
    """Interface for the IR Spectrum Stage 1 (stb-ir)"""
    print("\n" + "="*60)
    print(color_text("IR SPECTRUM - STAGE 1: PHONON DISPLACEMENTS", 'bold').center(60))
    print("="*60 + "\n")

    structure_file = get_input("Input structure file [default: structure.fdf]: ").strip()
    if not structure_file:
        structure_file = "structure.fdf"

    calc_file = get_input("Calculation parameters file [default: calc.fdf]: ").strip()
    if not calc_file:
        calc_file = "calc.fdf"

    dim_input = get_input("\nSupercell dimensions for the phonon calc (e.g. '2 2 2') "
                           "[default: 2 2 2]: ").strip()
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

    distance = get_float_input("\nPhonon displacement distance in Ang [default: 0.02]: ", 0.02)

    pseudo_dir = prompt_pseudo_source(optional=True)
    if not pseudo_dir:
        pseudo_dir = "."

    output_dir = get_input("\nOutput root directory [default: ir_study]: ").strip()
    if not output_dir:
        output_dir = "ir_study"

    vacuum_gap = 10.0
    show_advanced = get_input("\nConfigure advanced settings (vacuum-gap)? [y/N]: ").strip().lower()
    if show_advanced == 'y':
        vacuum_gap = get_float_input(
            "Vacuum gap threshold in Ang, for the supercell-dimension advisory and the Stage 2 "
            "dipole-path preview (default: 10.0): ", 10.0)

    args = [
        "-s", structure_file, "-c", calc_file,
        "-dim", str(dim_x), str(dim_y), str(dim_z),
        "-d", str(distance), "-p", pseudo_dir,
        "--vacuum-gap", str(vacuum_gap),
        "-O", output_dir, "--no-intro"
    ]

    print(color_text("\nGenerating IR workflow phonon displacement folders...", 'green'))
    run_tool("stb-ir", args)


def run_ir_modes() -> None:
    """Interface for the IR Spectrum Stage 2 (stb-irModes)"""
    print("\n" + "="*60)
    print(color_text("IR SPECTRUM - STAGE 2: MODES & DISPLACEMENTS", 'bold').center(60))
    print("="*60 + "\n")

    run_dir = get_input("Directory written by Stage 1 [default: ir_study]: ").strip()
    if not run_dir:
        run_dir = "ir_study"

    calc_file = get_input("Calc.fdf template for the follow-up calculation(s): ")
    while not os.path.isfile(calc_file):
        print(color_text("File not found!", 'red'))
        calc_file = get_input("Calc.fdf template for the follow-up calculation(s): ")

    args = ["--directory", run_dir, "--calc", calc_file, "--no-intro"]

    modes_str = get_input(
        "\nMode indices to process, space-separated [optional, default: every "
        "non-acoustic mode]: ").strip()
    if modes_str:
        args.extend(["--modes"] + modes_str.split())

    if not modes_str:
        use_symmetry_choice = get_input(
            "\nUse group theory to skip modes that are symmetry-forbidden from being "
            "IR-active (saves SIESTA time on the non-bulk path, needs a primitive cell) "
            "(y/N): ").strip().lower()
        if use_symmetry_choice in ('y', 'yes'):
            args.append("--use-symmetry")

    skip_degenerate_choice = get_input(
        "\nFor degenerate modes (2+ bands sharing one irrep), compute only ONE representative "
        "and reuse its IR intensity for the rest (rotational invariant, identical for every "
        "degenerate partner) -- saves SIESTA time on the non-bulk path only, needs a primitive "
        "cell (y/N): ").strip().lower()
    if skip_degenerate_choice in ('y', 'yes'):
        args.append("--skip-degenerate")

    displacement = get_float_input("\nFinite-difference displacement in Ang [default: 0.02]: ", 0.02)
    args.extend(["--displacement", str(displacement)])

    animations_choice = get_input(
        "\nExport a looping animation (.axsf, viewable in XCrySDen/VESTA) of each selected "
        "mode's eigendisplacement (y/N): ").strip().lower()
    if animations_choice in ('y', 'yes'):
        args.append("--export-animations")
        animation_frames = get_int_input("Frames per mode animation [default: 20]: ", 20)
        args.extend(["--animation-frames", str(animation_frames)])

    show_advanced = get_input(
        "\nConfigure advanced settings (frequency range/pseudopotential override/vacuum-gap/"
        "rotational-mode tolerance/bulk Born-charge settings)? [y/N]: ").strip().lower()
    if show_advanced == 'y':
        freq_min = get_input("Skip modes below this frequency in THz [optional]: ").strip()
        if freq_min:
            args.extend(["--freq-min", freq_min])
        freq_max = get_input("Skip modes above this frequency in THz [optional]: ").strip()
        if freq_max:
            args.extend(["--freq-max", freq_max])
        print(color_text(
            "\nPseudopotential source (blank = reuse whatever Stage 1 already copied into its "
            "phonon_disp/disp-*/ folders -- the normal case):", 'yellow'))
        pseudo_dir = prompt_pseudo_source(optional=True)
        if pseudo_dir:
            args.extend(["--pseudo-dir", pseudo_dir])
        vacuum_gap = get_float_input(
            "\nVacuum gap threshold in Ang -- also used to auto-select the bulk vs. non-bulk "
            "IR path (default: 10.0): ", 10.0)
        args.extend(["--vacuum-gap", str(vacuum_gap)])
        rotational_mode_tol = get_float_input(
            "0D only: frequency threshold in THz below which a mode is treated as free rotation "
            "rather than a real vibration (default: 2.0): ", 2.0)
        args.extend(["--rotational-mode-tol", str(rotational_mode_tol)])
        print(color_text(
            "\nBulk (3D) path only -- ignored for a non-bulk structure:", 'yellow'))
        fc_displ = get_float_input("MD.FCDispl in Bohr [default: 0.02]: ", 0.02)
        args.extend(["--fc-displ", str(fc_displ)])
        pol_grid_input = get_input(
            "PolarizationGrids, 9 integers row-major 3x3 [default: 10 4 4 4 10 4 4 4 10]: "
        ).strip()
        if pol_grid_input:
            try:
                pol_grid_vals = [int(x) for x in pol_grid_input.split()]
                if len(pol_grid_vals) == 9:
                    args.extend(["--polarization-grid"] + [str(v) for v in pol_grid_vals])
                else:
                    print(color_text("Expected 9 integers -- using default.", 'yellow'))
            except ValueError:
                print(color_text("Invalid input -- using default.", 'yellow'))

    print(color_text("\nBuilding phonon force constants and IR displacement folder(s)...", 'green'))
    run_tool("stb-irModes", args)


def run_ir_analysis() -> None:
    """Interface for the IR Spectrum Stage 3 (stb-irAnalysis)"""
    print("\n" + "="*60)
    print(color_text("IR SPECTRUM - STAGE 3: ANALYSIS", 'bold').center(60))
    print("="*60 + "\n")

    run_dir = get_input("Directory with the stb-ir/stb-irModes run [default: ir_study]: ").strip()
    if not run_dir:
        run_dir = "ir_study"

    output_filename = get_input("SIESTA output filename inside each displacement folder "
                                 "[default: calc.out]: ").strip()
    if not output_filename:
        output_filename = "calc.out"

    args = ["--directory", run_dir, "--file", output_filename, "--no-intro"]

    linewidth = get_float_input("\nLorentzian linewidth for the spectrum, cm^-1 [default: 10.0]: ", 10.0)
    args.extend(["--linewidth", str(linewidth)])

    temperature = get_input(
        "\nWeight the spectrum by the Bose-Einstein thermal population factor at a given "
        "temperature in K [optional, blank = no weighting]: ").strip()
    if temperature:
        args.extend(["--temperature", temperature])

    experimental_file = get_input(
        "\nOverlay a measured IR spectrum (2-column text file: frequency cm^-1, intensity) "
        "[optional, blank = skip]: ").strip()
    if experimental_file:
        args.extend(["--experimental", experimental_file])
        peak_prominence = get_float_input(
            "Minimum peak prominence, as a fraction of each spectrum's own max intensity "
            "[default: 0.01]: ", 0.01)
        args.extend(["--peak-prominence", str(peak_prominence)])

    run_tool("stb-irAnalysis", args)


def run_her_prep() -> None:
    """Interface for the HER Stage 1 (stb-her)"""
    print("\n" + "="*60)
    print(color_text("HER WORKFLOW - STAGE 1: ADSORPTION SITES", 'bold').center(60))
    print("="*60 + "\n")

    structure_file = get_input("Relaxed slab/2D structure file [default: structure.fdf]: ").strip()
    if not structure_file:
        structure_file = "structure.fdf"

    calc_file = get_input("Calc.fdf template for the site relaxations: ")
    while not os.path.isfile(calc_file):
        print(color_text("File not found!", 'red'))
        calc_file = get_input("Calc.fdf template for the site relaxations: ")

    pseudo_dir = prompt_pseudo_source(optional=True)

    site_type = get_input(
        "\nSite type to search: ontop/bridge/hollow/all [default: all]: ").strip().lower()
    if site_type not in ("ontop", "bridge", "hollow", "all"):
        site_type = "all"

    height = get_float_input("\nH adsorption height in Ang [default: 1.5]: ", 1.5)

    both_sides_choice = get_input(
        "\nAdsorb on both faces (free-standing 2D material with vacuum on both sides) "
        "(y/N): ").strip().lower()

    output_dir = get_input("\nOutput root directory [default: her_study]: ").strip()
    if not output_dir:
        output_dir = "her_study"

    args = [
        "-s", structure_file, "-c", calc_file,
        "--site-type", site_type, "--height", str(height),
        "-O", output_dir, "--no-intro"
    ]
    if pseudo_dir:
        args.extend(["-p", pseudo_dir])
    if both_sides_choice in ('y', 'yes'):
        args.append("--both-sides")

    show_advanced = get_input(
        "\nConfigure advanced settings (symprec/vacuum-gap)? [y/N]: ").strip().lower()
    if show_advanced == 'y':
        symprec = get_float_input(
            "Symmetry-reduction tolerance for site-finding [default: 0.01]: ", 0.01)
        args.extend(["--symprec", str(symprec)])
        vacuum_gap = get_float_input(
            "Vacuum-axis detection threshold in Ang [default: 10.0]: ", 10.0)
        args.extend(["--vacuum-gap", str(vacuum_gap)])

    print(color_text("\nSearching H-adsorption sites...", 'green'))
    run_tool("stb-her", args)


def run_her_refs() -> None:
    """Interface for the HER Stage 2 (stb-herRefs)"""
    print("\n" + "="*60)
    print(color_text("HER WORKFLOW - STAGE 2: REFERENCES & ZPE PREP", 'bold').center(60))
    print("="*60 + "\n")

    run_dir = get_input("Directory written by Stage 1 [default: her_study]: ").strip()
    if not run_dir:
        run_dir = "her_study"

    output_filename = get_input("SIESTA output filename inside each site folder "
                                 "[default: calc.out]: ").strip()
    if not output_filename:
        output_filename = "calc.out"

    pseudo_dir = prompt_pseudo_source(optional=True)

    zpe_mode = get_input(
        "\nZPE/entropy mode: standard/local/full [default: local]: ").strip().lower()
    if zpe_mode not in ("standard", "local", "full"):
        zpe_mode = "local"

    args = ["--directory", run_dir, "--file", output_filename, "--zpe-mode", zpe_mode, "--no-intro"]
    if pseudo_dir:
        args.extend(["-p", pseudo_dir])

    show_advanced = get_input(
        "\nConfigure advanced settings (displacement/supercell/vacuum-box)? [y/N]: ").strip().lower()
    if show_advanced == 'y':
        displacement = get_float_input("Finite-difference displacement in Ang [default: 0.015]: ", 0.015)
        args.extend(["--displacement", str(displacement)])
        vacuum_box = get_float_input("H2 reference vacuum box side in Ang [default: 15.0]: ", 15.0)
        args.extend(["--vacuum-box", str(vacuum_box)])
        if zpe_mode == "full":
            supercell_input = get_input(
                "Supercell dimensions for the full phonon calc (e.g. '1 1 1') "
                "[default: 1 1 1]: ").strip()
            if supercell_input:
                try:
                    dims = [int(x) for x in supercell_input.split()]
                    if len(dims) == 3:
                        args.extend(["--supercell", str(dims[0]), str(dims[1]), str(dims[2])])
                except ValueError:
                    pass

    print(color_text("\nBuilding winning-site references and ZPE folders...", 'green'))
    run_tool("stb-herRefs", args)


def run_her_analysis() -> None:
    """Interface for the HER Stage 3 (stb-herAnalysis)"""
    print("\n" + "="*60)
    print(color_text("HER WORKFLOW - STAGE 3: ANALYSIS", 'bold').center(60))
    print("="*60 + "\n")

    run_dir = get_input("Directory with the stb-her/stb-herRefs run [default: her_study]: ").strip()
    if not run_dir:
        run_dir = "her_study"

    output_filename = get_input("SIESTA output filename inside each folder "
                                 "[default: calc.out]: ").strip()
    if not output_filename:
        output_filename = "calc.out"

    temperature = get_float_input("\nTemperature in K [default: 298.15]: ", 298.15)

    force_tolerance = get_float_input(
        "Force tolerance for the 'is this folder relaxed/converged' check, in eV/Ang "
        "(default: 0.05): ", 0.05)

    args = ["--directory", run_dir, "--file", output_filename, "--temp", str(temperature),
            "--force-tolerance", str(force_tolerance), "--no-intro"]

    run_tool("stb-herAnalysis", args)


def run_oer_prep() -> None:
    """Interface for the OER Stage 1 (stb-oer)"""
    print("\n" + "="*60)
    print(color_text("OER WORKFLOW - STAGE 1: ADSORPTION SITES (OH*)", 'bold').center(60))
    print("="*60 + "\n")

    structure_file = get_input("Relaxed slab/2D structure file [default: structure.fdf]: ").strip()
    if not structure_file:
        structure_file = "structure.fdf"

    calc_file = get_input("Calc.fdf template for the site relaxations: ")
    while not os.path.isfile(calc_file):
        print(color_text("File not found!", 'red'))
        calc_file = get_input("Calc.fdf template for the site relaxations: ")

    pseudo_dir = prompt_pseudo_source(optional=True)

    site_type = get_input(
        "\nSite type to search: ontop/bridge/hollow/all [default: all]: ").strip().lower()
    if site_type not in ("ontop", "bridge", "hollow", "all"):
        site_type = "all"

    height = get_float_input("\nOH adsorption height in Ang [default: 1.8]: ", 1.8)

    both_sides_choice = get_input(
        "\nAdsorb on both faces (free-standing 2D material with vacuum on both sides) "
        "(y/N): ").strip().lower()

    output_dir = get_input("\nOutput root directory [default: oer_study]: ").strip()
    if not output_dir:
        output_dir = "oer_study"

    args = [
        "-s", structure_file, "-c", calc_file,
        "--site-type", site_type, "--height", str(height),
        "-O", output_dir, "--no-intro"
    ]
    if pseudo_dir:
        args.extend(["-p", pseudo_dir])
    if both_sides_choice in ('y', 'yes'):
        args.append("--both-sides")

    show_advanced = get_input(
        "\nConfigure advanced settings (O-H bond length/symprec/vacuum-gap)? [y/N]: ").strip().lower()
    if show_advanced == 'y':
        oh_bond_length = get_float_input(
            "O-H bond length of the adsorbing OH group, in Ang [default: 0.970]: ", 0.970)
        args.extend(["--oh-bond-length", str(oh_bond_length)])
        symprec = get_float_input(
            "Symmetry-reduction tolerance for site-finding [default: 0.01]: ", 0.01)
        args.extend(["--symprec", str(symprec)])
        vacuum_gap = get_float_input(
            "Vacuum-axis detection threshold in Ang [default: 10.0]: ", 10.0)
        args.extend(["--vacuum-gap", str(vacuum_gap)])

    print(color_text("\nSearching OH-adsorption sites...", 'green'))
    run_tool("stb-oer", args)


def run_oer_intermediates() -> None:
    """Interface for the OER Stage 2 (stb-oerIntermediates)"""
    print("\n" + "="*60)
    print(color_text("OER WORKFLOW - STAGE 2: O*/OOH* INTERMEDIATES", 'bold').center(60))
    print("="*60 + "\n")

    run_dir = get_input("Directory written by Stage 1 [default: oer_study]: ").strip()
    if not run_dir:
        run_dir = "oer_study"

    output_filename = get_input("SIESTA output filename inside each site folder "
                                 "[default: calc.out]: ").strip()
    if not output_filename:
        output_filename = "calc.out"

    pseudo_dir = prompt_pseudo_source(optional=True)

    o_strategy = get_input(
        "\nO* starting geometry: derived (from the winning OH* site, cheap)/search (independent "
        "site search, more rigorous) [default: derived]: ").strip().lower()
    if o_strategy not in ("derived", "search"):
        o_strategy = "derived"

    ooh_strategy = get_input(
        "OOH* starting geometry: derived/search [default: derived]: ").strip().lower()
    if ooh_strategy not in ("derived", "search"):
        ooh_strategy = "derived"

    ml_prerelax_choice = get_input(
        "\nPre-relax O*/OOH*'s adsorbate atoms with MACE-MP-0 before writing the CG-relaxation "
        "folder(s) (substrate fixed)? Needs the optional 'ml' extra (y/N): ").strip().lower()
    ml_prerelax = ml_prerelax_choice in ('y', 'yes')

    args = ["--directory", run_dir, "--file", output_filename,
            "--o-strategy", o_strategy, "--ooh-strategy", ooh_strategy, "--no-intro"]
    if pseudo_dir:
        args.extend(["-p", pseudo_dir])
    if ml_prerelax:
        args.append("--ml-prerelax")

    if o_strategy == "search":
        print(color_text("\nO* search settings:", 'yellow'))
        o_site_type = get_input("  Site type: ontop/bridge/hollow/all [default: all]: ").strip().lower()
        if o_site_type not in ("ontop", "bridge", "hollow", "all"):
            o_site_type = "all"
        args.extend(["--o-site-type", o_site_type])
        o_height = get_float_input("  O adsorption height in Ang [default: 1.5]: ", 1.5)
        args.extend(["--o-height", str(o_height)])

    if ooh_strategy == "search":
        print(color_text("\nOOH* search settings:", 'yellow'))
        ooh_site_type = get_input("  Site type: ontop/bridge/hollow/all [default: all]: ").strip().lower()
        if ooh_site_type not in ("ontop", "bridge", "hollow", "all"):
            ooh_site_type = "all"
        args.extend(["--ooh-site-type", ooh_site_type])
        ooh_height = get_float_input("  OOH adsorption height in Ang [default: 2.0]: ", 2.0)
        args.extend(["--ooh-height", str(ooh_height)])

    show_advanced = get_input(
        "\nConfigure advanced settings (OOH* starting bond lengths/bend angle)? [y/N]: ").strip().lower()
    if show_advanced == 'y':
        oo_bond_length = get_float_input(
            "Illustrative starting O-O bond length for OOH*, in Ang [default: 1.45]: ", 1.45)
        args.extend(["--oo-bond-length", str(oo_bond_length)])
        ooh_oh_bond_length = get_float_input(
            "Illustrative starting O-H bond length for OOH*'s new H, in Ang [default: 0.970]: ", 0.970)
        args.extend(["--ooh-oh-bond-length", str(ooh_oh_bond_length)])
        ooh_bend_deg = get_float_input(
            "Illustrative starting O-O-H bend angle for OOH*, in degrees [default: 100.0]: ", 100.0)
        args.extend(["--ooh-bend-deg", str(ooh_bend_deg)])

    print(color_text("\nBuilding O*/OOH* intermediate geometries...", 'green'))
    run_tool("stb-oerIntermediates", args)


def run_oer_refs() -> None:
    """Interface for the OER Stage 3 (stb-oerRefs)"""
    print("\n" + "="*60)
    print(color_text("OER WORKFLOW - STAGE 3: REFERENCES, BSSE & ZPE PREP", 'bold').center(60))
    print("="*60 + "\n")

    run_dir = get_input("Directory written by Stage 1/2 [default: oer_study]: ").strip()
    if not run_dir:
        run_dir = "oer_study"

    output_filename = get_input("SIESTA output filename inside each folder "
                                 "[default: calc.out]: ").strip()
    if not output_filename:
        output_filename = "calc.out"

    pseudo_dir = prompt_pseudo_source(optional=True)

    bsse_mode = get_input(
        "\nBSSE mode: shared (one triad at OOH*'s geometry, cheaper)/full (one triad per "
        "intermediate) [default: shared]: ").strip().lower()
    if bsse_mode not in ("shared", "full"):
        bsse_mode = "shared"

    zpe_mode = get_input(
        "\nZPE/entropy mode: local/full [default: local]: ").strip().lower()
    if zpe_mode not in ("local", "full"):
        zpe_mode = "local"

    args = ["--directory", run_dir, "--file", output_filename, "--bsse-mode", bsse_mode,
            "--zpe-mode", zpe_mode, "--no-intro"]
    if pseudo_dir:
        args.extend(["-p", pseudo_dir])

    show_advanced = get_input(
        "\nConfigure advanced settings (displacement/supercell/vacuum-box)? [y/N]: ").strip().lower()
    if show_advanced == 'y':
        displacement = get_float_input("Finite-difference displacement in Ang [default: 0.015]: ", 0.015)
        args.extend(["--displacement", str(displacement)])
        vacuum_box = get_float_input("H2/H2O reference vacuum box side in Ang [default: 15.0]: ", 15.0)
        args.extend(["--vacuum-box", str(vacuum_box)])
        if zpe_mode == "full":
            supercell_input = get_input(
                "Supercell dimensions for the full phonon calc (e.g. '1 1 1') "
                "[default: 1 1 1]: ").strip()
            if supercell_input:
                try:
                    dims = [int(x) for x in supercell_input.split()]
                    if len(dims) == 3:
                        args.extend(["--supercell", str(dims[0]), str(dims[1]), str(dims[2])])
                except ValueError:
                    pass

    print(color_text("\nBuilding references, BSSE and ZPE folders...", 'green'))
    run_tool("stb-oerRefs", args)


def run_oer_analysis() -> None:
    """Interface for the OER Stage 4 (stb-oerAnalysis)"""
    print("\n" + "="*60)
    print(color_text("OER WORKFLOW - STAGE 4: ANALYSIS", 'bold').center(60))
    print("="*60 + "\n")

    run_dir = get_input("Directory with the stb-oer/stb-oerIntermediates/stb-oerRefs run "
                         "[default: oer_study]: ").strip()
    if not run_dir:
        run_dir = "oer_study"

    output_filename = get_input("SIESTA output filename inside each folder "
                                 "[default: calc.out]: ").strip()
    if not output_filename:
        output_filename = "calc.out"

    temperature = get_float_input("\nTemperature in K [default: 298.15]: ", 298.15)

    force_tolerance = get_float_input(
        "Force tolerance for the 'is this folder relaxed/converged' check, in eV/Ang "
        "(default: 0.05): ", 0.05)

    args = ["--directory", run_dir, "--file", output_filename, "--temp", str(temperature),
            "--force-tolerance", str(force_tolerance), "--no-intro"]

    run_tool("stb-oerAnalysis", args)


def run_gqca_prep() -> None:
    """Interface for the GQCA Stage 1 (stb-gqca)"""
    print("\n" + "="*60)
    print(color_text("GQCA WORKFLOW - STAGE 1: CLUSTER STRUCTURE GENERATION", 'bold').center(60))
    print("="*60 + "\n")
    print(color_text(
        "Builds the 3 ordered pair-cluster-representative structures (AA/AB/BB) for a "
        "substitutional alloy A(1-x)Bx -- works on any lattice, 2D or 3D alike.", 'cyan'))
    print()

    structure_file = get_input("Input structure file (2D or 3D) [default: structure.fdf]: ").strip()
    if not structure_file:
        structure_file = "structure.fdf"

    sublattice = get_input("\nExisting species whose sites become the alloy sublattice "
                            "(e.g. Mo): ").strip()
    species_b = get_input("Second alloy species (e.g. W): ").strip()

    calc_file = get_input("\nCalc.fdf template for the 3 cluster relaxations: ")
    while not os.path.isfile(calc_file):
        print(color_text("File not found!", 'red'))
        calc_file = get_input("Calc.fdf template for the 3 cluster relaxations: ")

    pseudo_dir = prompt_pseudo_source(optional=True)

    output_dir = get_input("\nOutput root directory [default: gqca_study]: ").strip()
    if not output_dir:
        output_dir = "gqca_study"

    args = [
        "-f", structure_file, "--sublattice", sublattice, "--species-b", species_b,
        "-c", calc_file, "-O", output_dir, "--no-intro"
    ]
    if pseudo_dir:
        args.extend(["-p", pseudo_dir])

    show_advanced = get_input(
        "\nConfigure advanced settings (k-point density/vacuum-gap)? [y/N]: ").strip().lower()
    if show_advanced == 'y':
        kdensity = get_float_input("k-point density, same convention as stb-kgrid "
                                    "[default: 0.04]: ", 0.04)
        args.extend(["--kdensity", str(kdensity)])
        vacuum_gap = get_float_input(
            "Vacuum-axis detection threshold in Ang [default: 10.0]: ", 10.0)
        args.extend(["--vacuum-gap", str(vacuum_gap)])

    print(color_text("\nGenerating AA/AB/BB cluster structures...", 'green'))
    run_tool("stb-gqca", args)


def run_gqca_analysis() -> None:
    """Interface for the GQCA Stage 2 (stb-gqcaAnalysis)"""
    print("\n" + "="*60)
    print(color_text("GQCA WORKFLOW - STAGE 2: MASS-ACTION ANALYSIS", 'bold').center(60))
    print("="*60 + "\n")

    run_dir = get_input("Directory written by Stage 1 [default: gqca_study]: ").strip()
    if not run_dir:
        run_dir = "gqca_study"

    output_filename = get_input("SIESTA output filename inside each cluster_n* folder "
                                 "[default: calc.out]: ").strip()
    if not output_filename:
        output_filename = "calc.out"

    print(f"\n{color_text('Temperature:', 'yellow')}")
    print(f"  {color_text('1', 'cyan')} = Single fixed temperature (default)")
    print(f"  {color_text('2', 'cyan')} = Temperature sweep")
    temp_choice = get_input("Select option (1-2) [default: 1]: ").strip()

    args = ["--directory", run_dir, "--file", output_filename, "--no-intro"]
    if temp_choice == '2':
        temp_min = get_float_input("  T min (K): ")
        temp_max = get_float_input("  T max (K): ")
        temp_points = get_int_input("  Number of T points [default: 5]: ", 5)
        args.extend(["--temp-min", str(temp_min), "--temp-max", str(temp_max),
                     "--temp-points", str(temp_points)])
    else:
        temp = get_float_input("  Temperature in K [default: 300.0]: ", 300.0)
        args.extend(["--temp", str(temp)])

    show_advanced = get_input(
        "\nConfigure advanced settings (composition grid/force tolerance/output name)? "
        "[y/N]: ").strip().lower()
    if show_advanced == 'y':
        x_min = get_float_input("  x min [default: 0.0]: ", 0.0)
        x_max = get_float_input("  x max [default: 1.0]: ", 1.0)
        x_points = get_int_input("  Number of x points [default: 21]: ", 21)
        args.extend(["--x-min", str(x_min), "--x-max", str(x_max), "--x-points", str(x_points)])
        force_tolerance = get_float_input(
            "  Force tolerance for the 'is this folder relaxed' check, in eV/Ang "
            "[default: 0.05]: ", 0.05)
        args.extend(["--force-tolerance", str(force_tolerance)])
        output = get_input("  Base filename for .csv/.dat/.gplot [default: gqca_results]: ").strip()
        if output:
            args.extend(["--output", output])

    run_tool("stb-gqcaAnalysis", args)


def run_optical_prep() -> None:
    """Interface for the Optical Properties Stage 1 (stb-optical)"""
    print("\n" + "="*60)
    print(color_text("OPTICAL PROPERTIES - STAGE 1: DIRECTION FOLDERS", 'bold').center(60))
    print("="*60 + "\n")
    print(color_text(
        "Writes one SIESTA folder per requested Cartesian polarization direction, ready for a "
        "linear-response Optical.Calculate run -- no new geometry, same already-relaxed structure "
        "in every folder.", 'cyan'))
    print()

    structure_file = get_input("Already-relaxed input structure file [default: structure.fdf]: ").strip()
    if not structure_file:
        structure_file = "structure.fdf"

    calc_file = get_input("\nCalc.fdf template for the Optical calculations: ")
    while not os.path.isfile(calc_file):
        print(color_text("File not found!", 'red'))
        calc_file = get_input("Calc.fdf template for the Optical calculations: ")

    pseudo_dir = prompt_pseudo_source(optional=True)

    directions_str = get_input(
        "\nDirections to compute, space-separated (x/y/z) [default: x y z]: ").strip()
    directions = directions_str.split() if directions_str else ["x", "y", "z"]

    output_dir = get_input("\nOutput root directory [default: optical_study]: ").strip()
    if not output_dir:
        output_dir = "optical_study"

    args = [
        "-f", structure_file, "-c", calc_file, "--directions", *directions,
        "-O", output_dir, "--no-intro"
    ]
    if pseudo_dir:
        args.extend(["-p", pseudo_dir])

    show_advanced = get_input(
        "\nConfigure advanced settings (Optical.Mesh/Optical.Broaden/Optical.NumberOfBands)? "
        "[y/N]: ").strip().lower()
    if show_advanced == 'y':
        mesh_str = get_input("Optical.Mesh, 3 integers [default: 10 10 10]: ").strip()
        if mesh_str:
            mesh_vals = mesh_str.split()
            if len(mesh_vals) == 3:
                args.extend(["--optical-mesh", *mesh_vals])
            else:
                print(color_text("Expected 3 integers -- using default.", 'yellow'))
        broaden = get_float_input("Optical.Broaden in eV [default: 0.2]: ", 0.2)
        args.extend(["--optical-broaden", str(broaden)])
        nbands_str = get_input("Optical.NumberOfBands [optional, blank = SIESTA default]: ").strip()
        if nbands_str:
            args.extend(["--optical-nbands", nbands_str])

    print(color_text("\nWriting direction folder(s)...", 'green'))
    run_tool("stb-optical", args)


def run_optical_analysis() -> None:
    """Interface for the Optical Properties Stage 2 (stb-opticalAnalysis)"""
    print("\n" + "="*60)
    print(color_text("OPTICAL PROPERTIES - STAGE 2: ANALYSIS", 'bold').center(60))
    print("="*60 + "\n")

    run_dir = get_input("Directory written by Stage 1 [default: optical_study]: ").strip()
    if not run_dir:
        run_dir = "optical_study"

    output_filename = get_input("SIESTA output filename inside each direction folder "
                                 "[default: calc.out]: ").strip()
    if not output_filename:
        output_filename = "calc.out"

    print(f"\n{color_text('Quantity to plot (overlaid across every direction present):', 'yellow')}")
    print(f"  {color_text('1', 'cyan')} = eps1 + eps2, dielectric function (default)")
    print(f"  {color_text('2', 'cyan')} = n + k, refractive index / extinction coefficient")
    print(f"  {color_text('3', 'cyan')} = alpha, absorption coefficient")
    print(f"  {color_text('4', 'cyan')} = R, normal-incidence reflectivity")
    print(f"  {color_text('5', 'cyan')} = L, energy-loss function")
    print(f"  {color_text('6', 'cyan')} = sigma1, real optical conductivity")
    print(f"  {color_text('7', 'cyan')} = all -- write every quantity's .gplot in one run")
    quantity_choice = get_input("Select option (1-7) [default: 1]: ").strip()
    quantity_map = {'1': 'eps', '2': 'n_k', '3': 'alpha', '4': 'R', '5': 'L', '6': 'sigma1', '7': 'all'}
    plot_quantity = quantity_map.get(quantity_choice, 'eps')

    args = ["--directory", run_dir, "--file", output_filename,
            "--plot-quantity", plot_quantity, "--no-intro"]

    if plot_quantity not in ('eps', 'n_k', 'all'):
        experimental_file = get_input(
            "\nOverlay a measured spectrum (2-column text file: energy eV, value) "
            "[optional, blank = skip]: ").strip()
        if experimental_file:
            args.extend(["--experimental", experimental_file])
            peak_prominence = get_float_input(
                "Minimum peak prominence, as a fraction of the spectrum's own max "
                "[default: 0.01]: ", 0.01)
            args.extend(["--peak-prominence", str(peak_prominence)])

    dimcorr_choice = get_input(
        "\nRestore the intrinsic (vacuum-independent) response for a 2D/0D structure "
        "(needs --thickness for 2D) [y/N]: ").strip().lower()
    if dimcorr_choice in ('y', 'yes'):
        args.append("--dimensionality-correction")
        thickness = get_input(
            "Physical thickness of the 2D material in Ang [blank if the structure is 0D/bulk]: "
        ).strip()
        if thickness:
            args.extend(["--thickness", thickness])

    show_advanced = get_input(
        "\nConfigure advanced settings (energy range/output name)? [y/N]: ").strip().lower()
    if show_advanced == 'y':
        energy_min = get_input("Energy min in eV [optional, blank = no trim]: ").strip()
        if energy_min:
            args.extend(["--energy-min", energy_min])
        energy_max = get_input("Energy max in eV [optional, blank = no trim]: ").strip()
        if energy_max:
            args.extend(["--energy-max", energy_max])
        output = get_input("Base filename for .csv/.dat/.gplot [default: optical_results]: ").strip()
        if output:
            args.extend(["--output", output])

    run_tool("stb-opticalAnalysis", args)


def run_mlff_prep() -> None:
    """Interface for the Custom ML Force Field Stage 1 (stb-mlff)"""
    print("\n" + "="*60)
    print(color_text("CUSTOM ML FORCE FIELD - STAGE 1: PREP", 'bold').center(60))
    print("="*60 + "\n")
    print(color_text(
        "Generates SIESTA training configurations (rattled displacements and/or AIMD-sampled "
        "snapshots) from a reference calculation -- run SIESTA in each folder, then use Stage "
        "2 to fine-tune a MACE-MP foundation model on the results.", 'cyan'))
    print()

    calc_file = get_input("Reference calc.fdf (basis/mesh cutoff/pseudos/SCF settings reused as-is): ").strip()
    while not os.path.isfile(calc_file):
        print(color_text("File not found!", 'red'))
        calc_file = get_input("Reference calc.fdf: ").strip()

    pseudo_dir = prompt_pseudo_source(optional=True) or "."

    n_configs = get_int_input("Number of rattled configurations [default: 10]: ", 10)
    stdev_str = get_input("Rattling standard deviation(s) in Ang, space-separated [default: 0.05]: ").strip()
    stdevs = stdev_str.split() if stdev_str else ["0.05"]

    args = ["--file", calc_file, "--pseudo-dir", pseudo_dir,
            "--n-configs", str(n_configs), "--stdev", *stdevs, "--no-intro"]

    use_aimd = get_input(
        "\nAlso sample starting geometries from an existing AIMD trajectory (<label>.ANI)? "
        "[y/N]: ").strip().lower()
    if use_aimd == 'y':
        aimd_label = get_input("AIMD SystemLabel: ").strip()
        n_samples = get_int_input("Number of AIMD frames to sample [default: 5]: ", 5)
        args.extend(["--from-aimd", aimd_label, "--n-samples", str(n_samples)])

    output_dir = get_input("\nOutput directory [default: .]: ").strip() or "."
    args.extend(["--output-dir", output_dir])

    save_report = get_input("Also save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")

    print(color_text("\nGenerating training configurations...", 'green'))
    run_tool("stb-mlff", args)


def run_mlff_analysis() -> None:
    """Interface for the Custom ML Force Field Stage 2 (stb-mlffAnalysis)"""
    print("\n" + "="*60)
    print(color_text("CUSTOM ML FORCE FIELD - STAGE 2: FINE-TUNING", 'bold').center(60))
    print("="*60 + "\n")
    print(color_text(
        "Aggregates the finished SIESTA calculations from Stage 1 (energy + forces, "
        "best-effort stress) and fine-tunes a MACE-MP foundation model on them via "
        "'mace_run_train'.", 'cyan'))
    print()

    path = get_input("Glob pattern for the config folders [default: mlff_config_*]: ").strip()
    if not path:
        path = "mlff_config_*"

    print(f"\n{color_text('Foundation model to fine-tune from:', 'yellow')}")
    print(f"  {color_text('1', 'cyan')} = small (default -- fastest)")
    print(f"  {color_text('2', 'cyan')} = medium")
    print(f"  {color_text('3', 'cyan')} = large")
    model_choice = get_input("Select option (1-3) [default: 1]: ").strip()
    model_map = {'1': 'small', '2': 'medium', '3': 'large'}
    foundation_model = model_map.get(model_choice, 'small')

    epochs = get_int_input("Max training epochs [default: 100]: ", 100)
    batch_size = get_int_input("Batch size [default: 5]: ", 5)

    device = get_input("Device (cpu/cuda) [default: auto-detect]: ").strip().lower()

    name = get_input("Experiment/model name [default: mlff_model]: ").strip() or "mlff_model"
    work_dir = get_input("Working directory [default: .]: ").strip() or "."

    args = ["--path", path, "--foundation-model", foundation_model,
            "--epochs", str(epochs), "--batch-size", str(batch_size),
            "--name", name, "--work-dir", work_dir, "--no-intro"]
    if device in ("cpu", "cuda"):
        args.extend(["--device", device])

    export_lammps = get_input(
        "\nAlso export a LAMMPS-loadable model (needs a LAMMPS build with the MACE pair "
        "style)? (y/N): ").strip().lower()
    if export_lammps == 'y':
        args.append("--export-lammps")

    save_report = get_input("Also save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")

    print(color_text("\nFine-tuning the model (this can take a while)...", 'green'))
    run_tool("stb-mlffAnalysis", args)


def run_mlff_active_learning() -> None:
    """Interface for the Custom ML Force Field Stage 3 (stb-mlffActiveLearning)"""
    print("\n" + "="*60)
    print(color_text("CUSTOM ML FORCE FIELD - STAGE 3: ACTIVE LEARNING", 'bold').center(60))
    print("="*60 + "\n")
    print(color_text(
        "Uses the fine-tuned model itself to screen a fresh batch of candidate "
        "configurations and suggests the ones it's least confident about (largest "
        "predicted force) as new SIESTA calculations to run -- NOT rigorous uncertainty "
        "quantification, a cheap single-model heuristic.", 'cyan'))
    print()

    model_path = get_input("Fine-tuned model file (from stb-mlffAnalysis): ").strip()
    while not os.path.isfile(model_path):
        print(color_text("File not found!", 'red'))
        model_path = get_input("Fine-tuned model file: ").strip()

    calc_file = get_input("Reference calc.fdf (same one stb-mlff used): ").strip()
    while not os.path.isfile(calc_file):
        print(color_text("File not found!", 'red'))
        calc_file = get_input("Reference calc.fdf: ").strip()

    pseudo_dir = prompt_pseudo_source(optional=True) or "."

    print(f"\n{color_text('Candidate sampling method:', 'yellow')}")
    print(f"  {color_text('1', 'cyan')} = Rattle -- independent, static displaced snapshots (default, cheap)")
    print(f"  {color_text('2', 'cyan')} = MD -- short NVT trajectory with the model itself "
          "(stb-mlmd's own dynamics, more physically representative)")
    method_choice = get_input("Select option (1-2) [default: 1]: ").strip()
    sampling_method = {'1': 'rattle', '2': 'md'}.get(method_choice, 'rattle')

    args = ["--model", model_path, "--file", calc_file, "--pseudo-dir", pseudo_dir,
            "--sampling-method", sampling_method, "--no-intro"]

    if sampling_method == "rattle":
        n_candidates = get_int_input("Number of candidates to screen [default: 50]: ", 50)
        stdev_str = get_input("Rattling standard deviation(s) in Ang, space-separated "
                              "[default: 0.15]: ").strip()
        stdevs = stdev_str.split() if stdev_str else ["0.15"]
        args.extend(["--n-candidates", str(n_candidates), "--stdev", *stdevs])
    else:
        md_steps = get_int_input("Total MD steps to run [default: 500]: ", 500)
        md_temperature = get_float_input("NVT temperature, K [default: 800]: ", 800.0)
        stride = get_int_input("Score one frame every Nth MD step [default: 10]: ", 10)
        args.extend(["--md-steps", str(md_steps), "--md-temperature", str(md_temperature),
                     "--stride", str(stride)])

    top_k = get_int_input("Number of highest-scoring candidates to write out [default: 5]: ", 5)
    output_dir = get_input("Output directory [default: .]: ").strip() or "."
    args.extend(["--top-k", str(top_k), "--output-dir", output_dir])

    save_report = get_input("Also save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")

    print(color_text("\nScreening candidate configurations...", 'green'))
    run_tool("stb-mlffActiveLearning", args)


def run_eos_inputs_generator() -> None:
    """Interface for the Equation of State Prep (stb-eosInputs)"""
    print("\n" + "="*60)
    print(color_text("EQUATION OF STATE - STAGE 1: PREP", 'bold').center(60))
    print("="*60 + "\n")

    struct_file = get_input("Input structure file (bulk, 3D periodic, fractional coords): ")
    while not os.path.isfile(struct_file):
        print(color_text("File not found!", 'red'))
        struct_file = get_input("Input structure file: ")

    n_volumes = get_int_input("Number of volumes to scan [default: 9, min 5]: ", 9)
    strain_range = get_float_input("Max isotropic volume strain, %% [default: 5.0]: ", 5.0)
    output_dir = get_input("Output directory [default: eos_runs]: ").strip() or "eos_runs"

    args = ["--file", struct_file, "--n-volumes", str(n_volumes),
            "--strain-range", str(strain_range), "--output-dir", output_dir, "--no-intro"]

    run_tool("stb-eosInputs", args)


def run_eos_analysis_generator() -> None:
    """Interface for the Equation of State Analysis (stb-eosAnalysis)"""
    print("\n" + "="*60)
    print(color_text("EQUATION OF STATE - STAGE 2: ANALYSIS", 'bold').center(60))
    print("="*60 + "\n")

    run_dir = get_input("Directory with 'vol_*' folders [default: eos_runs]: ").strip() or "eos_runs"
    output_filename = get_input("SIESTA output filename inside each folder [default: calc.out]: ").strip() \
        or "calc.out"
    eos = get_input("EOS form, birchmurnaghan/vinet/murnaghan/sj [default: birchmurnaghan]: ").strip()
    curve_output = get_input("Base filename for the .dat/.gplot outputs [default: eos_curve]: ").strip() \
        or "eos_curve"

    args = ["--dir", run_dir, "--file", output_filename, "--output", curve_output, "--no-intro"]
    if eos:
        args.extend(["--eos", eos])

    run_tool("stb-eosAnalysis", args)


def run_mlmd_generator() -> None:
    """Interface for ML Molecular Dynamics (stb-mlmd)"""
    print("\n" + "="*60)
    print(color_text("ML MOLECULAR DYNAMICS (stb-mlmd)", 'bold').center(60))
    print("="*60 + "\n")
    print(color_text(
        "Needs the optional 'ml' extra (pip install stb_suite[ml]). NVE/NVT/NPT MD "
        "driven by MACE-MP-0 or a custom fine-tuned model -- a fast heuristic, not a "
        "substitute for ab initio MD.", 'yellow'))

    input_file = get_input("\nInput structure file (fdf): ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input structure file (fdf): ")

    custom_model = get_input(
        "\nUse a custom MACE model instead (e.g. one fine-tuned via stb-mlffAnalysis)? "
        "Path, or Enter to use a MACE-MP-0 foundation model: ").strip()

    args = ["--file", input_file, "--no-intro"]
    if custom_model:
        while not os.path.isfile(custom_model):
            print(color_text("File not found!", 'red'))
            custom_model = get_input("Custom model path: ").strip()
        args.extend(["--custom-model", custom_model])
    else:
        model = get_input("Model size, small/medium/large [default: small]: ").strip()
        args.extend(["--model", model or "small"])

    print(f"\n{color_text('Ensemble:', 'yellow')}")
    print(f"  {color_text('1', 'cyan')} = NVT (Langevin thermostat, default)")
    print(f"  {color_text('2', 'cyan')} = NVE (microcanonical, no thermostat)")
    print(f"  {color_text('3', 'cyan')} = NPT (Berendsen thermostat+barostat, bulk only)")
    ensemble_choice = get_input("Select option (1-3) [default: 1]: ").strip()
    ensemble = {'1': 'nvt', '2': 'nve', '3': 'npt'}.get(ensemble_choice, 'nvt')
    args.extend(["--ensemble", ensemble])

    temperature = get_float_input("Temperature, K [default: 300]: ", 300.0)
    args.extend(["--temperature", str(temperature)])
    if ensemble == "npt":
        pressure = get_float_input("Pressure, GPa [default: 0.0]: ", 0.0)
        args.extend(["--pressure", str(pressure)])

    equilibration_steps = get_int_input(
        "Equilibration steps to run and discard before production [default: 0]: ", 0)
    n_steps = get_int_input("Number of PRODUCTION MD steps [default: 1000]: ", 1000)
    timestep = get_float_input("Timestep, fs [default: 1.0]: ", 1.0)
    args.extend(["--equilibration-steps", str(equilibration_steps),
                 "--n-steps", str(n_steps), "--timestep", str(timestep)])

    show_advanced = get_input(
        "\nConfigure advanced settings (thermostat/barostat coupling, compressibility)? "
        "[y/N]: ").strip().lower()
    if show_advanced == 'y':
        if ensemble in ("nvt", "npt"):
            taut = get_float_input("Thermostat coupling time, fs [default: 100]: ", 100.0)
            args.extend(["--taut", str(taut)])
        if ensemble == "npt":
            taup = get_float_input("Barostat coupling time, fs [default: 500]: ", 500.0)
            compressibility = get_float_input(
                "Compressibility, eV/Ang^3 [default: 4.57e-5 -- water's value]: ", 4.57e-5)
            args.extend(["--taup", str(taup), "--compressibility", str(compressibility)])

    output = get_input("\nTrajectory filename [default: <input>_md_traj.xsf]: ").strip()
    if output:
        args.extend(["--output", output])

    save_data = get_input("Also save the raw energy/temperature diagnostics data (.dat)? (y/N): ").strip().lower()
    if save_data == 'y':
        args.append("--save-data")

    run_tool("stb-mlmd", args)


def run_mlphonons_generator() -> None:
    """Interface for ML Phonons (stb-mlphonons)"""
    print("\n" + "="*60)
    print(color_text("ML PHONONS (stb-mlphonons)", 'bold').center(60))
    print("="*60 + "\n")
    print(color_text(
        "Needs the optional 'ml' extra (pip install stb_suite[ml]). Standalone phonon "
        "calculation (bands, DOS, thermal properties) driven entirely by a MACE potential -- "
        "no SIESTA, no separate analysis tool. A fast heuristic preview, not a substitute "
        "for the real SIESTA-based Phonons workflow.", 'yellow'))

    input_file = get_input("\nInput structure file (fdf): ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input structure file (fdf): ")

    custom_model = get_input(
        "\nUse a custom MACE model instead (e.g. one fine-tuned via stb-mlffAnalysis)? "
        "Path, or Enter to use a MACE-MP-0 foundation model: ").strip()

    args = ["--file", input_file, "--no-intro"]
    if custom_model:
        while not os.path.isfile(custom_model):
            print(color_text("File not found!", 'red'))
            custom_model = get_input("Custom model path: ").strip()
        args.extend(["--custom-model", custom_model])
    else:
        model = get_input("Model size, small/medium/large [default: small]: ").strip()
        args.extend(["--model", model or "small"])

    dim_str = get_input("Supercell dimensions, 3 integers [default: 2 2 2]: ").strip()
    dim = dim_str.split() if dim_str else ["2", "2", "2"]
    args.extend(["--dim", *dim])

    output_dir = get_input("Output directory [default: mlphonons_out]: ").strip()
    if output_dir:
        args.extend(["--output-dir", output_dir])

    print(f"\n{color_text('Run mode:', 'yellow')}")
    print(f"  {color_text('1', 'cyan')} = Single-volume (bands, DOS, thermal properties) (default)")
    print(f"  {color_text('2', 'cyan')} = Quasi-harmonic approximation (QHA), volume scan")
    print(f"  {color_text('3', 'cyan')} = Supercell (--dim) convergence check")
    mode_choice = get_input("Select option (1-3) [default: 1]: ").strip() or '1'

    if mode_choice == '2':
        args.append("--qha")
        n_volumes = get_input("Number of volumes to scan, >= 4 [default: 5]: ").strip()
        args.extend(["--n-volumes", n_volumes or "5"])
        strain_range = get_input(
            "Max isotropic strain around the reference volume, %% [default: 3.0]: ").strip()
        args.extend(["--strain-range", strain_range or "3.0"])
        eos = get_input(
            "Equation of state, vinet/murnaghan/birch_murnaghan [default: vinet]: ").strip()
        args.extend(["--eos", eos or "vinet"])
    elif mode_choice == '3':
        args.append("--check-convergence")
        conv_dims_str = get_input(
            "Supercell sizes to compare, as space-separated NA NB NC triples "
            "(e.g. 1 1 1 2 2 2 3 3 3) [default: --dim and --dim+1 on every "
            "non-vacuum axis]: ").strip()
        if conv_dims_str:
            args.extend(["--convergence-dims", *conv_dims_str.split()])
    else:
        if custom_model:
            skip_compare = get_input(
                "\nAlso compare against the MACE-MP-0 foundation model on the same plots? "
                "(Y/n): ").strip().lower()
            if skip_compare == 'n':
                args.append("--skip-foundation-comparison")

        freeze = get_input(
            "\nIf a genuine imaginary (unstable) mode is found, write a structure displaced "
            "along it (--freeze-unstable-mode)? (y/N): ").strip().lower()
        if freeze == 'y':
            args.append("--freeze-unstable-mode")
            freeze_amp = get_input("Displacement amplitude, Ang [default: 0.05]: ").strip()
            if freeze_amp:
                args.extend(["--freeze-amplitude", freeze_amp])

        animate = get_input(
            "\nAnimate a normal mode's vibration as a trajectory for OVITO/VMD "
            "(--animate-mode)? (y/N): ").strip().lower()
        if animate == 'y':
            band_idx = get_input(
                "Band index to animate (0 = lowest frequency at the q-point): ").strip()
            args.extend(["--animate-mode", band_idx or "0"])
            qpoint_str = get_input(
                "Fractional q-point, 3 numbers [default: 0 0 0, Gamma]: ").strip()
            if qpoint_str:
                args.extend(["--animate-qpoint", *qpoint_str.split()])
            animate_dim_str = get_input(
                "Modulation supercell size, 3 integers [default: 1 1 1, only exact at "
                "Gamma]: ").strip()
            if animate_dim_str:
                args.extend(["--animate-dim", *animate_dim_str.split()])
            animate_amp = get_input("Animation amplitude, Ang [default: 0.5]: ").strip()
            if animate_amp:
                args.extend(["--animate-amplitude", animate_amp])
            animate_frames = get_input("Number of frames [default: 30]: ").strip()
            if animate_frames:
                args.extend(["--animate-frames", animate_frames])
            animate_fmt = get_input("Output format, xsf/pdb/xyz [default: xsf]: ").strip()
            args.extend(["--animate-format", animate_fmt or "xsf"])

    save_data = get_input("\nAlso save the raw plot data (.dat files)? (y/N): ").strip().lower()
    if save_data == 'y':
        args.append("--save-data")
    save_report = get_input("Also save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")

    print(color_text("\nComputing ML phonons...", 'green'))
    run_tool("stb-mlphonons", args)


def run_mlelastic_generator() -> None:
    """Interface for ML Elastic Constants (stb-mlelastic)"""
    print("\n" + "="*60)
    print(color_text("ML ELASTIC CONSTANTS (stb-mlelastic)", 'bold').center(60))
    print("="*60 + "\n")
    print(color_text(
        "Needs the optional 'ml' extra (pip install stb_suite[ml]). Standalone stiffness "
        "matrix (stress-strain method) driven entirely by a MACE potential -- no SIESTA, no "
        "separate strain_*/ folders. A fast heuristic preview, not a substitute for the real "
        "SIESTA-based Elastic workflow.", 'yellow'))

    input_file = get_input("\nInput structure file (fdf): ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input structure file (fdf): ")

    custom_model = get_input(
        "\nUse a custom MACE model instead (e.g. one fine-tuned via stb-mlffAnalysis)? "
        "Path, or Enter to use a MACE-MP-0 foundation model: ").strip()

    args = ["--file", input_file, "--no-intro"]
    if custom_model:
        while not os.path.isfile(custom_model):
            print(color_text("File not found!", 'red'))
            custom_model = get_input("Custom model path: ").strip()
        args.extend(["--custom-model", custom_model])
    else:
        model = get_input("Model size, small/medium/large [default: small]: ").strip()
        args.extend(["--model", model or "small"])

    max_strain = get_input("Max strain magnitude, %% [default: 1.0]: ").strip()
    if max_strain:
        args.extend(["--max", max_strain])
    steps = get_input("Strain points per direction [default: 4]: ").strip()
    if steps:
        args.extend(["--steps", steps])

    dimensionality = get_input(
        "Dimensionality, auto/3d/2d/1d [default: auto]: ").strip()
    if dimensionality:
        args.extend(["--dimensionality", dimensionality])

    output_dir = get_input("Output directory [default: mlelastic_out]: ").strip()
    if output_dir:
        args.extend(["--output-dir", output_dir])

    check_convergence = get_input(
        "\nInstead run a strain-amplitude convergence check (multiple --max values)? "
        "(y/N): ").strip().lower()
    if check_convergence == 'y':
        args.append("--check-convergence")
        conv_strains = get_input(
            "Strain magnitudes to compare, space-separated %% [default: half/same/double "
            "the max above]: ").strip()
        if conv_strains:
            args.extend(["--convergence-strains", *conv_strains.split()])
    elif custom_model:
        skip_compare = get_input(
            "\nAlso compare against the MACE-MP-0 foundation model on the same plot/table? "
            "(Y/n): ").strip().lower()
        if skip_compare == 'n':
            args.append("--skip-foundation-comparison")

    if check_convergence != 'y':
        skip_aniso = get_input(
            "\nSkip the directional Young's-modulus anisotropy surface (3D only)? "
            "(y/N): ").strip().lower()
        if skip_aniso == 'y':
            args.append("--skip-anisotropy-surface")

    save_data = get_input("\nAlso save the raw stress-strain data (.dat)? (y/N): ").strip().lower()
    if save_data == 'y':
        args.append("--save-data")
    save_report = get_input("Also save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")

    print(color_text("\nComputing ML elastic constants...", 'green'))
    run_tool("stb-mlelastic", args)


def run_mlsearch_generator() -> None:
    """Interface for ML Structure Search (stb-mlsearch)"""
    print("\n" + "="*60)
    print(color_text("ML STRUCTURE SEARCH (stb-mlsearch)", 'bold').center(60))
    print("="*60 + "\n")
    print(color_text(
        "Needs the optional 'ml' extra (pip install stb_suite[ml]). Global structure "
        "search for the lowest-energy atomic arrangement at a fixed cell -- e.g. defect/"
        "vacancy site occupation, disordered/alloy configurations. NOT crystal structure "
        "prediction (cell shape untouched during the search) and NOT a proof of global "
        "optimality.", 'yellow'))

    input_file = get_input("\nInput structure file (fdf): ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input structure file (fdf): ")

    custom_model = get_input(
        "\nUse a custom MACE model instead (e.g. one fine-tuned via stb-mlffAnalysis)? "
        "Path, or Enter to use a MACE-MP-0 foundation model: ").strip()

    args = ["--file", input_file, "--no-intro"]
    if custom_model:
        while not os.path.isfile(custom_model):
            print(color_text("File not found!", 'red'))
            custom_model = get_input("Custom model path: ").strip()
        args.extend(["--custom-model", custom_model])
    else:
        model = get_input("Model size, small/medium/large [default: small]: ").strip()
        args.extend(["--model", model or "small"])

    print(f"\n{color_text('Algorithm:', 'yellow')}")
    print(f"  {color_text('1', 'cyan')} = Basin hopping (default)")
    print(f"  {color_text('2', 'cyan')} = Simulated annealing (MD-based, reuses stb-mlmd's engine)")
    algo_choice = get_input("Select option (1-2) [default: 1]: ").strip() or '1'

    if algo_choice == '2':
        args.extend(["--algorithm", "simulated-annealing"])
        steps = get_input("Total MD steps in the cooling trajectory [default: 50]: ").strip()
        if steps:
            args.extend(["--steps", steps])
        temp_start = get_input("Starting temperature, K [default: 1500]: ").strip()
        if temp_start:
            args.extend(["--temp-start", temp_start])
        temp_end = get_input("Final temperature, K [default: 50]: ").strip()
        if temp_end:
            args.extend(["--temp-end", temp_end])
        schedule = get_input("Cooling schedule, exponential/linear [default: exponential]: ").strip()
        if schedule:
            args.extend(["--schedule", schedule])
        snapshot_interval = get_input("MD steps between snapshot evaluations [default: 20]: ").strip()
        if snapshot_interval:
            args.extend(["--snapshot-interval", snapshot_interval])
    else:
        steps = get_input("Number of basin-hopping steps [default: 50]: ").strip()
        if steps:
            args.extend(["--steps", steps])
        dr = get_input("Max random displacement per step, Ang [default: 0.3]: ").strip()
        if dr:
            args.extend(["--dr", dr])
        temperature = get_input("Search temperature, K [default: 1000]: ").strip()
        if temperature:
            args.extend(["--temperature", temperature])

    seed = get_input("\nRandom seed for reproducibility [default: unseeded]: ").strip()
    if seed:
        args.extend(["--seed", seed])

    output_dir = get_input("Output directory [default: mlsearch_out]: ").strip()
    if output_dir:
        args.extend(["--output-dir", output_dir])

    save_data = get_input("\nAlso save the raw search-history data (.dat)? (y/N): ").strip().lower()
    if save_data == 'y':
        args.append("--save-data")
    save_report = get_input("Also save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")

    print(color_text("\nRunning ML structure search...", 'green'))
    run_tool("stb-mlsearch", args)


def run_mlmelting_generator() -> None:
    """Interface for ML Melting Point (stb-mlmelting)"""
    print("\n" + "="*60)
    print(color_text("ML MELTING POINT (stb-mlmelting)", 'bold').center(60))
    print("="*60 + "\n")
    print(color_text(
        "Needs the optional 'ml' extra (pip install stb_suite[ml]). Brackets a bulk "
        "material's melting point via a sequence of short MACE MD runs, tracking the "
        "Lindemann index vs. temperature. Bulk (3D periodic) only. A coarse heuristic, "
        "known to superheat past the true melting point (small, defect-free periodic "
        "cell -- no free surface to nucleate melting from).", 'yellow'))

    input_file = get_input("\nInput structure file (fdf): ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input structure file (fdf): ")

    custom_model = get_input(
        "\nUse a custom MACE model instead (e.g. one fine-tuned via stb-mlffAnalysis)? "
        "Path, or Enter to use a MACE-MP-0 foundation model: ").strip()

    args = ["--file", input_file, "--no-intro"]
    if custom_model:
        while not os.path.isfile(custom_model):
            print(color_text("File not found!", 'red'))
            custom_model = get_input("Custom model path: ").strip()
        args.extend(["--custom-model", custom_model])
    else:
        model = get_input("Model size, small/medium/large [default: small]: ").strip()
        args.extend(["--model", model or "small"])

    temp_choice = get_input(
        "\nTemperatures to scan: [1] range (min/max/step), [2] explicit list "
        "[default: 1]: ").strip() or '1'
    if temp_choice == '2':
        temps_str = get_input("Temperatures, space-separated K: ").strip()
        if temps_str:
            args.extend(["--temperatures", *temps_str.split()])
    else:
        temp_min = get_input("Min temperature, K [default: 300]: ").strip()
        if temp_min:
            args.extend(["--temp-min", temp_min])
        temp_max = get_input("Max temperature, K [default: 3000]: ").strip()
        if temp_max:
            args.extend(["--temp-max", temp_max])
        temp_step = get_input("Temperature step, K [default: 300]: ").strip()
        if temp_step:
            args.extend(["--temp-step", temp_step])

    eq_steps = get_input("Equilibration steps per temperature [default: 200]: ").strip()
    if eq_steps:
        args.extend(["--equilibration-steps", eq_steps])
    prod_steps = get_input("Production steps per temperature [default: 500]: ").strip()
    if prod_steps:
        args.extend(["--production-steps", prod_steps])

    output_dir = get_input("Output directory [default: mlmelting_out]: ").strip()
    if output_dir:
        args.extend(["--output-dir", output_dir])

    save_data = get_input("\nAlso save the raw Lindemann/diffusion data (.dat)? (y/N): ").strip().lower()
    if save_data == 'y':
        args.append("--save-data")
    save_report = get_input("Also save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")

    print(color_text("\nRunning ML melting-point scan...", 'green'))
    run_tool("stb-mlmelting", args)


def run_mlconvergence_generator() -> None:
    """Interface for ML Model-Size Convergence (stb-mlconvergence)"""
    print("\n" + "="*60)
    print(color_text("ML MODEL-SIZE CONVERGENCE (stb-mlconvergence)", 'bold').center(60))
    print("="*60 + "\n")
    print(color_text(
        "Needs the optional 'ml' extra (pip install stb_suite[ml]). Compares MACE-MP-0 "
        "model sizes (small/medium/large) on the same reference calculation (a full "
        "relax by default) -- helps decide whether a bigger, slower model actually "
        "changes the answer for this structure.", 'yellow'))

    input_file = get_input("\nInput structure file (fdf): ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input structure file (fdf): ")

    args = ["--file", input_file, "--no-intro"]

    models_str = get_input(
        "\nModel sizes to compare, space-separated (small/medium/large) "
        "[default: small medium large]: ").strip()
    if models_str:
        args.extend(["--models", *models_str.split()])

    custom_str = get_input(
        "Also include custom model file(s)? Space-separated paths, or Enter to skip: ").strip()
    if custom_str:
        args.extend(["--custom-models", *custom_str.split()])

    energy_tol = get_input("Energy tolerance, meV/atom [default: 5.0]: ").strip()
    if energy_tol:
        args.extend(["--energy-tolerance", energy_tol])
    lattice_tol = get_input("Lattice/volume tolerance, %% [default: 1.0]: ").strip()
    if lattice_tol:
        args.extend(["--lattice-tolerance", lattice_tol])

    output_dir = get_input("Output directory [default: mlconvergence_out]: ").strip()
    if output_dir:
        args.extend(["--output-dir", output_dir])

    save_data = get_input("\nAlso save the raw comparison data (.dat)? (y/N): ").strip().lower()
    if save_data == 'y':
        args.append("--save-data")
    save_report = get_input("Also save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")

    print(color_text("\nRunning ML model-size convergence check...", 'green'))
    run_tool("stb-mlconvergence", args)


def run_mlneb_generator() -> None:
    """Interface for ML NEB (stb-mlneb)"""
    print("\n" + "="*60)
    print(color_text("ML NEB (stb-mlneb)", 'bold').center(60))
    print("="*60 + "\n")
    print(color_text(
        "Needs the optional 'ml' extra (pip install stb_suite[ml]). Standalone climbing-"
        "image NEB (reaction/migration barrier estimate) driven entirely by a MACE "
        "potential -- no SIESTA, no calc.fdf needed. Requires two already-relaxed "
        "endpoint structures with the exact same composition (e.g. a vacancy/adsorbate "
        "at two different sites).", 'yellow'))

    initial_file = get_input("\nInitial (relaxed) endpoint structure (fdf): ")
    while not os.path.isfile(initial_file):
        print(color_text("File not found!", 'red'))
        initial_file = get_input("Initial structure file (fdf): ")
    final_file = get_input("Final (relaxed) endpoint structure (fdf): ")
    while not os.path.isfile(final_file):
        print(color_text("File not found!", 'red'))
        final_file = get_input("Final structure file (fdf): ")

    custom_model = get_input(
        "\nUse a custom MACE model instead (e.g. one fine-tuned via stb-mlffAnalysis)? "
        "Path, or Enter to use a MACE-MP-0 foundation model: ").strip()

    args = ["--initial", initial_file, "--final", final_file, "--no-intro"]
    if custom_model:
        while not os.path.isfile(custom_model):
            print(color_text("File not found!", 'red'))
            custom_model = get_input("Custom model path: ").strip()
        args.extend(["--custom-model", custom_model])
    else:
        model = get_input("Model size, small/medium/large [default: small]: ").strip()
        args.extend(["--model", model or "small"])

    n_images = get_input("Total images including both endpoints [default: 7]: ").strip()
    if n_images:
        args.extend(["--n-images", n_images])

    freeze = get_input(
        "\nFreeze substrate atoms that barely move between endpoints? (y/N): ").strip().lower()
    if freeze == 'y':
        args.append("--freeze-substrate")

    output_dir = get_input("Output directory [default: mlneb_out]: ").strip()
    if output_dir:
        args.extend(["--output-dir", output_dir])

    save_images = get_input("\nAlso save each relaxed image as its own .fdf? (y/N): ").strip().lower()
    if save_images == 'y':
        args.append("--save-images")
    save_report = get_input("Also save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")

    print(color_text("\nRunning ML NEB...", 'green'))
    run_tool("stb-mlneb", args)


def run_mldiffusion_generator() -> None:
    """Interface for ML Vacancy Diffusion (stb-mldiffusion)"""
    print("\n" + "="*60)
    print(color_text("ML VACANCY DIFFUSION (stb-mldiffusion)", 'bold').center(60))
    print("="*60 + "\n")
    print(color_text(
        "Needs the optional 'ml' extra (pip install stb_suite[ml]). Given a PERFECT bulk "
        "structure and the index of one atom to remove, auto-detects the neighboring hop "
        "targets and screens the vacancy migration barrier for each via MACE-driven NEB.",
        'yellow'))

    input_file = get_input("\nPerfect bulk reference structure (fdf): ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input structure file (fdf): ")

    vacancy_index = get_input("0-based index of the atom to remove (vacancy site): ").strip()
    while not vacancy_index:
        print(color_text("This is required.", 'red'))
        vacancy_index = get_input("0-based index of the atom to remove (vacancy site): ").strip()

    custom_model = get_input(
        "\nUse a custom MACE model instead (e.g. one fine-tuned via stb-mlffAnalysis)? "
        "Path, or Enter to use a MACE-MP-0 foundation model: ").strip()

    args = ["--file", input_file, "--vacancy-index", vacancy_index, "--no-intro"]
    if custom_model:
        while not os.path.isfile(custom_model):
            print(color_text("File not found!", 'red'))
            custom_model = get_input("Custom model path: ").strip()
        args.extend(["--custom-model", custom_model])
    else:
        model = get_input("Model size, small/medium/large [default: small]: ").strip()
        args.extend(["--model", model or "small"])

    n_images = get_input("Total images including both endpoints [default: 5]: ").strip()
    if n_images:
        args.extend(["--n-images", n_images])

    output_dir = get_input("Output directory [default: mldiffusion_out]: ").strip()
    if output_dir:
        args.extend(["--output-dir", output_dir])

    save_report = get_input("\nAlso save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")

    print(color_text("\nRunning ML vacancy migration-barrier screening...", 'green'))
    run_tool("stb-mldiffusion", args)


def run_mlgcmc_generator() -> None:
    """Interface for ML Monte Carlo / GCMC (stb-mlgcmc)"""
    print("\n" + "="*60)
    print(color_text("ML MONTE CARLO / GCMC (stb-mlgcmc)", 'bold').center(60))
    print("="*60 + "\n")
    print(color_text(
        "Needs the optional 'ml' extra (pip install stb_suite[ml]). Canonical/grand-"
        "canonical Monte Carlo adsorption (insert/move/delete a single-atom adsorbate "
        "in a rigid host framework) driven by a MACE potential.", 'yellow'))

    host_file = get_input("\nHost framework structure (fdf, treated as rigid): ")
    while not os.path.isfile(host_file):
        print(color_text("File not found!", 'red'))
        host_file = get_input("Host framework structure (fdf): ")

    adsorbate = get_input("Adsorbate element symbol (single atom, e.g. Ar): ").strip()
    while not adsorbate:
        print(color_text("This is required.", 'red'))
        adsorbate = get_input("Adsorbate element symbol: ").strip()

    args = ["--host", host_file, "--adsorbate", adsorbate, "--no-intro"]

    print(f"\n{color_text('Ensemble:', 'yellow')}")
    print(f"  {color_text('1', 'cyan')} = Grand-canonical (insert/move/delete at a chemical potential) (default)")
    print(f"  {color_text('2', 'cyan')} = Canonical (fixed adsorbate count, move only)")
    ensemble_choice = get_input("Select option (1-2) [default: 1]: ").strip() or '1'

    if ensemble_choice == '2':
        args.extend(["--ensemble", "canonical"])
        n_initial = get_input("Number of adsorbate atoms (fixed): ").strip()
        while not n_initial:
            print(color_text("This is required for the canonical ensemble.", 'red'))
            n_initial = get_input("Number of adsorbate atoms (fixed): ").strip()
        args.extend(["--n-initial", n_initial])
    else:
        mu = get_input(
            "Chemical potential, eV (or Enter to scan multiple values): ").strip()
        if mu:
            args.extend(["--mu", mu])
        else:
            mu_scan = get_input("Chemical potentials to scan, space-separated eV: ").strip()
            if mu_scan:
                args.extend(["--mu-scan", *mu_scan.split()])

    temperature = get_input("\nTemperature, K [default: 300]: ").strip()
    if temperature:
        args.extend(["--temperature", temperature])
    steps = get_input("Total MC steps [default: 2000]: ").strip()
    if steps:
        args.extend(["--steps", steps])
    equilibration_steps = get_input(
        "Equilibration steps to discard before collecting statistics [default: 500 -- "
        "must be less than the total steps above]: ").strip()
    if equilibration_steps:
        args.extend(["--equilibration-steps", equilibration_steps])

    output_dir = get_input("Output directory [default: mlgcmc_out]: ").strip()
    if output_dir:
        args.extend(["--output-dir", output_dir])

    save_data = get_input("\nAlso save the raw history/isotherm data (.dat)? (y/N): ").strip().lower()
    if save_data == 'y':
        args.append("--save-data")
    save_report = get_input("Also save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")

    print(color_text("\nRunning ML Monte Carlo...", 'green'))
    run_tool("stb-mlgcmc", args)


def run_mladsorb_generator() -> None:
    """Interface for ML Adsorption-Site Screening (stb-mladsorb)"""
    print("\n" + "="*60)
    print(color_text("ML ADSORPTION-SITE SCREENING (stb-mladsorb)", 'bold').center(60))
    print("="*60 + "\n")
    print(color_text(
        "Needs the optional 'ml' extra (pip install stb_suite[ml]). Enumerates candidate "
        "ontop/bridge/hollow adsorption sites, relaxes each (substrate fixed) with MACE, "
        "and ranks them by a genuine ML adsorption energy -- no SIESTA and no calc.fdf "
        "needed.", 'yellow'))

    slab_file = get_input("\nSlab/2D-material structure (fdf, vacuum along one axis): ")
    while not os.path.isfile(slab_file):
        print(color_text("File not found!", 'red'))
        slab_file = get_input("Slab structure file (fdf): ")

    adsorbate = get_input(
        "Adsorbate(s): ASE G2 molecule name or element symbol, comma-separated for "
        "multiple (e.g. 'H2O' or 'H,O'): ").strip()
    while not adsorbate:
        print(color_text("This is required.", 'red'))
        adsorbate = get_input("Adsorbate(s): ").strip()

    custom_model = get_input(
        "\nUse a custom MACE model instead (e.g. one fine-tuned via stb-mlffAnalysis)? "
        "Path, or Enter to use a MACE-MP-0 foundation model: ").strip()

    args = ["--slab", slab_file, "--adsorbate", adsorbate, "--no-intro"]
    if custom_model:
        while not os.path.isfile(custom_model):
            print(color_text("File not found!", 'red'))
            custom_model = get_input("Custom model path: ").strip()
        args.extend(["--custom-model", custom_model])
        skip_compare = get_input(
            "Also compare against the MACE-MP-0 foundation model on the same sites? "
            "(Y/n): ").strip().lower()
        if skip_compare == 'n':
            args.append("--skip-foundation-comparison")
    else:
        model = get_input("Model size, small/medium/large [default: small]: ").strip()
        args.extend(["--model", model or "small"])

    site_type = get_input(
        "Site type, all/ontop/bridge/hollow [default: all]: ").strip()
    if site_type:
        args.extend(["--site-type", site_type])
    height = get_input("Adsorbate height above the surface, Ang [default: 2.0]: ").strip()
    if height:
        args.extend(["--height", height])
    top_k = get_input("Only keep the N best sites (Enter for all): ").strip()
    if top_k:
        args.extend(["--top-k", top_k])

    n_orientations = get_input(
        "\nFor a multi-atom adsorbate, sample this many random orientations/site "
        "[default: 1, no sampling]: ").strip()
    if n_orientations:
        args.extend(["--n-orientations", n_orientations])

    diffusion = get_input(
        "Also estimate the adsorbate's surface-diffusion barrier between its 2 best "
        "sites? (y/N): ").strip().lower()
    if diffusion == 'y':
        args.append("--diffusion-barrier")

    output_dir = get_input("\nOutput directory [default: mladsorb_out]: ").strip()
    if output_dir:
        args.extend(["--output-dir", output_dir])

    save_data = get_input("Also save the raw ranking data (.dat)? (y/N): ").strip().lower()
    if save_data == 'y':
        args.append("--save-data")
    save_report = get_input("Also save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")

    print(color_text("\nRunning ML adsorption-site screening...", 'green'))
    run_tool("stb-mladsorb", args)


def run_mleos_generator() -> None:
    """Interface for ML Equation of State (stb-mleos)"""
    print("\n" + "="*60)
    print(color_text("ML EQUATION OF STATE (stb-mleos)", 'bold').center(60))
    print("="*60 + "\n")
    print(color_text(
        "Needs the optional 'ml' extra (pip install stb_suite[ml]). Scans isotropically-"
        "scaled volumes around equilibrium and fits an E(V) equation of state (Birch-"
        "Murnaghan by default) via MACE -- gives an independent, curvature-based bulk "
        "modulus to cross-check against stb-mlelastic's stress-strain-derived one.", 'yellow'))

    struct_file = get_input("\nBulk (3D periodic) structure file (fdf): ")
    while not os.path.isfile(struct_file):
        print(color_text("File not found!", 'red'))
        struct_file = get_input("Structure file (fdf): ")

    custom_model = get_input(
        "\nUse a custom MACE model instead (e.g. one fine-tuned via stb-mlffAnalysis)? "
        "Path, or Enter to use a MACE-MP-0 foundation model: ").strip()

    args = ["--file", struct_file, "--no-intro"]
    if custom_model:
        while not os.path.isfile(custom_model):
            print(color_text("File not found!", 'red'))
            custom_model = get_input("Custom model path: ").strip()
        args.extend(["--custom-model", custom_model])
        skip_compare = get_input(
            "Also fit with the MACE-MP-0 foundation model for comparison? (Y/n): ").strip().lower()
        if skip_compare == 'n':
            args.append("--skip-foundation-comparison")
    else:
        model = get_input("Model size, small/medium/large [default: small]: ").strip()
        args.extend(["--model", model or "small"])

    eos = get_input(
        "EOS form, birchmurnaghan/vinet/murnaghan/sj [default: birchmurnaghan]: ").strip()
    if eos:
        args.extend(["--eos", eos])
    n_volumes = get_input("Number of volumes to scan [default: 9, min 5]: ").strip()
    if n_volumes:
        args.extend(["--n-volumes", n_volumes])
    strain_range = get_input("Max volume strain, %% [default: 5.0]: ").strip()
    if strain_range:
        args.extend(["--strain-range", strain_range])

    output_dir = get_input("\nOutput directory [default: mleos_out]: ").strip()
    if output_dir:
        args.extend(["--output-dir", output_dir])

    save_data = get_input("Also save the raw volume-energy data (.dat)? (y/N): ").strip().lower()
    if save_data == 'y':
        args.append("--save-data")
    save_report = get_input("Also save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")

    print(color_text("\nRunning ML equation-of-state calculation...", 'green'))
    run_tool("stb-mleos", args)


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
        output = get_input("\nAlso save the block to a file? [optional, press Enter to skip]: ").strip()
        if output:
            args.extend(["--output", output])
        save_report = get_input("Also save a text report to file? (y/N): ").strip().lower()
        if save_report == 'y':
            args.append("--save-report")
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

    output = get_input("\nAlso save the block to a file? [optional, press Enter to skip]: ").strip()
    if output:
        args.extend(["--output", output])

    save_report = get_input("Also save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")

    run_tool("stb-dftu", args)


def run_dos_parser() -> None:
    """Interface for the PDOS XML Parser (stb-dos)"""
    print("\n" + "="*60)
    print(color_text("PDOS XML PARSER (stb-dos)", 'bold').center(60))
    print("="*60 + "\n")

    # Esta linha agora terá Tab-completion!
    label_or_path = get_input(
        "SIESTA label (e.g. 'siesta' for siesta.PDOS.xml [+ auto-detected "
        "siesta.bands/siesta.EIG]), or a direct .xml path: "
    )
    while True:
        if label_or_path.endswith(".xml") and os.path.isfile(label_or_path):
            input_file = label_or_path
            label = os.path.basename(input_file)
            if label.endswith(".PDOS.xml"):
                label = label[:-len(".PDOS.xml")]
            elif label.endswith(".xml"):
                label = label[:-len(".xml")]
            break
        elif os.path.isfile(f"{label_or_path}.PDOS.xml"):
            input_file = f"{label_or_path}.PDOS.xml"
            label = label_or_path
            break
        else:
            print(color_text(f"Neither '{label_or_path}' nor "
                              f"'{label_or_path}.PDOS.xml' was found!", 'red'))
            label_or_path = get_input("SIESTA label (or .xml path): ")

    has_bands = os.path.isfile(f"{label}.bands")
    has_eig = os.path.isfile(f"{label}.EIG")
    if has_bands:
        print(color_text(f"-> Found '{label}.bands', available as a VBM/CBM reference.", 'green'))
    elif has_eig:
        print(color_text(f"-> Found '{label}.EIG' (no '{label}.bands'), available as a "
                          "VBM/CBM reference.", 'green'))
    else:
        print(color_text(f"-> No '{label}.bands'/'{label}.EIG' found; --shift vbm/cbm would "
                          "need --estimate-from-dos.", 'yellow'))

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

    shift_options = {
        '1': ('fermi', "Fermi level"),
        '2': ('vbm', "Valence Band Maximum"),
        '3': ('cbm', "Conduction Band Minimum"),
        '4': ('manual', "Custom value"),
    }

    print(f"\n{color_text('Energy shift reference:', 'yellow')}")
    for key, (_, desc) in shift_options.items():
        print(f"  {color_text(key, 'yellow')}. {desc}" + (" [Default]" if key == '1' else ""))

    choice = get_input("\nSelect reference (1-4) [default: 1]: ").strip()
    if not choice:
        choice = '1'
    while choice not in shift_options:
        print(color_text("Invalid choice! Please select 1-4.", 'red'))
        choice = get_input("Select reference (1-4) [default: 1]: ").strip() or '1'

    shift_key, _ = shift_options[choice]
    if shift_key == 'manual':
        shift = str(get_float_input("Enter custom shift value (eV): "))
    else:
        shift = shift_key

    estimate_from_dos = False
    if shift.strip().lower() in ('vbm', 'cbm'):
        if has_bands:
            print(color_text(f"-> Will use '{label}.bands' to compute {shift.strip().upper()}.", 'green'))
        elif has_eig:
            print(color_text(f"-> Will use '{label}.EIG' to compute {shift.strip().upper()}.", 'green'))
        else:
            answer = get_input(
                "No '.bands'/'.EIG' found. Estimate VBM/CBM from the DOS itself instead? (y/N): "
            ).strip().lower()
            if answer == 'y':
                estimate_from_dos = True
            else:
                print(color_text("Falling back to --shift fermi.", 'yellow'))
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

    save_report = get_input("\nAlso save a text report to file? (y/N): ").strip().lower()
    save_gnuplot = get_input(
        "Also save a gnuplot (.gplot) script for each output category "
        "(total/atom/species)? (y/N): "
    ).strip().lower()
    view_plot = get_input(
        "Show an interactive matplotlib preview of the DOS before exiting? (y/N): "
    ).strip().lower()

    args = [
        input_file, # Positional argument
        "--shift", shift,
        "--type"
    ]

    args.extend(dos_types)
    args.extend(["--projection", projection_mode])
    args.extend(["--output-dir", output_dir])
    if estimate_from_dos:
        args.append("--estimate-from-dos")
    if save_report == 'y':
        args.append("--save-report")
    if save_gnuplot == 'y':
        args.append("--save-gnuplot")
    if view_plot == 'y':
        args.append("--view")
    args.append("--no-intro")

    run_tool("stb-dos", args)

def run_strain_generator() -> None:
    """Interface for the Strain Generator (stb-strain)"""
    print("\n" + "="*60)
    print(color_text("STRAIN GENERATOR (stb-strain)", 'bold').center(60))
    print("="*60)
    print(color_text(
        "Generates strained structures along one or more Cartesian directions. Vacuum-padded "
        "axes and symmetry-equivalent directions (both uniaxial and biaxial) are detected "
        "up front, so you only ever choose among the physically valid, non-redundant "
        "directions -- or run all of them at once.", 'cyan'))
    print()

    default_structure = "structure.fdf"
    if os.path.isfile(default_structure):
        structure_file = get_input(
            f"Input structure.fdf -- lattice + atomic coordinates [default: {default_structure}]: "
        ).strip() or default_structure
    else:
        structure_file = get_input("Input structure.fdf (lattice + atomic coordinates): ").strip()
    while not os.path.isfile(structure_file):
        print(color_text("File not found!", 'red'))
        structure_file = get_input("Input structure.fdf: ").strip()

    default_calc = "calc.fdf"
    if os.path.isfile(default_calc):
        calc_file = get_input(
            "Input calc.fdf -- SCF/basis/MD settings (normally %include-ing the structure file "
            f"above) [default: {default_calc}]: ").strip() or default_calc
    else:
        print(color_text(f"-> No '{default_calc}' found -- enter the real calc.fdf path below.", 'yellow'))
        calc_file = get_input(
            "Input calc.fdf (SCF/basis/MD settings, normally %include-ing the structure file): "
        ).strip()
    while not os.path.isfile(calc_file):
        print(color_text("File not found!", 'red'))
        calc_file = get_input("Input calc.fdf: ").strip()

    print(color_text(
        "\nRelax mode -- how the cell responds to the imposed strain (via a %block "
        "Geometry.Constraints written to config_extra.fdf):", 'cyan'))
    print("  1) cell-fixed          -- cell locked at the imposed strain; only atomic positions relax")
    print("  2) stress-constrained  -- only the imposed direction's stress fixed; other directions relax freely")
    relax_mode = ""
    while relax_mode not in ("cell-fixed", "stress-constrained"):
        choice = get_input("Choose 1 or 2: ").strip()
        relax_mode = {"1": "cell-fixed", "2": "stress-constrained"}.get(choice, choice)
        if relax_mode not in ("cell-fixed", "stress-constrained"):
            print(color_text("Invalid choice! Enter 1 or 2.", 'red'))

    pseudo_dir = prompt_pseudo_source(optional=True)

    # Advanced settings up front (rarely-touched tolerances -- gated so the
    # essential flow below stays short) -- collected here, not after the
    # direction question, because the axis-symmetry preview table below needs
    # the actual vacuum-gap/symmetry tolerances to be accurate.
    vacuum_gap, symprec, angle_tolerance, output_dir = 10.0, 0.01, 5.0, "strain_runs"
    show_advanced = get_input(
        "\nConfigure advanced settings (vacuum-gap, symmetry tolerances, output directory)? "
        "[y/N]: ").strip().lower()
    if show_advanced == 'y':
        vacuum_gap = get_float_input(
            "Vacuum gap threshold in Ang, for detecting periodic vs. vacuum-padded axes "
            "(default: 10.0): ", 10.0)
        symprec = get_float_input("Symmetry-detection tolerance (default: 0.01): ", 0.01)
        angle_tolerance = get_float_input("Symmetry angle tolerance in degrees (default: 5.0): ", 5.0)
        output_dir = get_input("Output directory (default: strain_runs): ").strip()
        if not output_dir:
            output_dir = "strain_runs"

    # Dimensionality: detected ONCE, upfront -- before the direction menu
    # below, which needs vacuum_axes to filter out vacuum-padded directions.
    # Also feeds the CONFIGURATION SUMMARY recap, so the structure file is
    # only ever read once.
    structure = None
    vacuum_axes = None
    dimensionality_label = "unknown (could not detect)"
    try:
        structure = structure_io.read_fdf(structure_file)
        positions = np.array([pos for _, pos in structure.atoms])
        is_cartesian = structure.coord_format == 'cartesian'
        frac_coords = kspace.to_fractional(positions, structure.lattice, is_cartesian)
        vacuum_axes = kspace.detect_vacuum_axes(frac_coords, structure.lattice, vacuum_gap)
        dimensionality_label = kspace.dimensionality_label(vacuum_axes)
        print(f"\n{color_text('[INFO]', 'green')} Detected dimensionality: {dimensionality_label}")
    except Exception:
        print(color_text(
            "\n[WARNING] Could not detect dimensionality (advisory only) -- "
            "continuing without it.", 'yellow'))

    # Detect vacuum-padded axes and symmetry-equivalent directions (both
    # uniaxial and biaxial -- see strain.py::print_direction_selection_table)
    # up front, so the direction menu below only ever offers directions that
    # are physically valid and not redundant with one another. Falls back to
    # the old free-text prompt (validated, but with no symmetry/vacuum
    # filtering) if the structure couldn't be read or the symmetry pre-check
    # fails for any reason.
    directions_to_run = None
    if structure is not None:
        try:
            pmg_structure = structure_io.to_pymatgen(structure)
            groups, point_group = symmetry.equivalent_cartesian_axes(
                pmg_structure, symprec, angle_tolerance)
            _, ops = symmetry.get_point_group_operations(pmg_structure, symprec, angle_tolerance)
            independent = strain.print_direction_selection_table(groups, point_group, ops, vacuum_axes)
            if independent:
                uniaxial_independent = [d for d in independent if len(d) == 1]
                biaxial_independent = [d for d in independent if len(d) == 2]

                # Individual directions first, then grouped "run several at
                # once" choices -- ALL UNIAXIAL / ALL BIAXIAL (only when that
                # kind has at least one independent direction) and, always
                # last, ALL of the above (uniaxial + biaxial together).
                menu_entries = [(d, [d]) for d in independent]
                if uniaxial_independent:
                    menu_entries.append((
                        f"ALL independent UNIAXIAL directions ({', '.join(uniaxial_independent)})",
                        uniaxial_independent))
                if biaxial_independent:
                    menu_entries.append((
                        f"ALL independent BIAXIAL directions ({', '.join(biaxial_independent)})",
                        biaxial_independent))
                menu_entries.append((
                    f"ALL independent directions -- uniaxial + biaxial ({len(independent)} total)",
                    independent))

                print()
                for i, (label, _) in enumerate(menu_entries, 1):
                    print(f"  {i}) {label}")
                n_options = len(menu_entries)
                choice = get_input(f"\nSelect a direction [1-{n_options}]: ").strip()
                while not choice.isdigit() or not (1 <= int(choice) <= n_options):
                    print(color_text(f"Invalid choice! Enter a number from 1 to {n_options}.", 'red'))
                    choice = get_input(f"Select a direction [1-{n_options}]: ").strip()
                choice = int(choice)
                directions_to_run = list(menu_entries[choice - 1][1])
            else:
                print(color_text(
                    "[WARNING] No independent, non-vacuum direction is available -- "
                    "falling back to manual direction entry.", 'yellow'))
        except Exception:
            print(color_text(
                "[WARNING] Could not run the symmetry/vacuum pre-check (advisory only) -- "
                "falling back to manual direction entry.", 'yellow'))

    if directions_to_run is None:
        direction = get_input("\nStrain direction (x,y,z,xy,xz,yz): ").lower()
        while not all(c in 'xyz' for c in direction) or len(direction) not in (1, 2):
            print(color_text("Invalid direction! Use x,y,z for uniaxial or xy,xz,yz for biaxial", 'red'))
            direction = get_input("Strain direction (x,y,z,xy,xz,yz): ").lower()
        directions_to_run = [direction]

    stmin = get_float_input("\nMinimum strain % (default 0): ", 0.0)
    stmax = get_float_input("Maximum strain % (default 25): ", 25.0)
    while stmax <= stmin:
        print(color_text("Maximum strain must be greater than minimum strain!", 'red'))
        stmax = get_float_input("Maximum strain % (default 25): ", 25.0)

    step = get_float_input("Step % (default 1): ", 1.0)
    while step <= 0:
        print(color_text("Step must be positive!", 'red'))
        step = get_float_input("Step % (default 1): ", 1.0)

    save_report = get_input("\nAlso save a report to file? [y/N]: ").strip().lower() == 'y'

    # Recap -- always shown, no confirmation gate
    summary_rows = [
        ("Structure file", structure_file),
        ("Calc template", calc_file),
        ("Relax mode", relax_mode),
        ("Pseudo source", pseudo_dir or "(not given)"),
        ("Dimensionality", dimensionality_label),
        ("Direction(s) to generate", ", ".join(directions_to_run)),
        ("Strain range", f"{stmin}% to {stmax}%, step {step}%"),
        ("Output directory", output_dir),
        ("Vacuum-gap threshold", f"{vacuum_gap} Ang"),
        ("Symmetry tolerance", f"symprec={symprec}, angle={angle_tolerance} deg"),
        ("Save report", "yes" if save_report else "no"),
    ]
    _print_config_summary("CONFIGURATION SUMMARY", summary_rows)

    for d in directions_to_run:
        args = [
            "-s", structure_file,
            "-c", calc_file,
            "--relax-mode", relax_mode,
            "--stdir", d,
            "--stmin", str(stmin),
            "--stmax", str(stmax),
            "--step", str(step),
            "--output-dir", output_dir,
            "--vacuum-gap", str(vacuum_gap),
            "--symprec", str(symprec),
            "--angle-tolerance", str(angle_tolerance),
            "--no-intro",
        ]
        if pseudo_dir:
            args.extend(["-p", pseudo_dir])
        if save_report:
            args.append("--save-report")
        run_tool("stb-strain", args, pause=(d == directions_to_run[-1]))

def run_strain_post_processor() -> None:
    """Interface for the Strain Post-Processing (strain_analysis.py)"""
    print("\n" + "="*60)
    print(color_text("STRAIN POST-PROCESSING ANALYZER", 'bold').center(60))
    print("="*60 + "\n")
    print(color_text(
        "Analyzes 'strain_*' folders inside a run directory -- reads either a single "
        "direction's own subfolder or the top-level output directory (every direction "
        "found under it is compared automatically, no flag needed). Dimensionality is "
        "auto-detected from a real structure file when one is found alongside the SIESTA "
        "output; falls back to 3D otherwise.", 'yellow'))

    args = []

    # 1. Diretório com as pastas strain_*
    run_dir = get_input("Directory with 'strain_*' folders [default: strain_runs]: ").strip()
    if not run_dir:
        run_dir = "strain_runs"
    args.extend(["--dir", run_dir])

    # 2. Nome do ficheiro de output dentro de cada pasta
    print(color_text("Enter the Siesta output filename located inside strain folders.", 'yellow'))
    siesta_out = get_input("Filename (e.g., calc.out): ").strip()
    while not siesta_out:
        print(color_text("Filename is required!", 'red'))
        siesta_out = get_input("Filename (e.g., calc.out): ").strip()

    args.extend(["--file", siesta_out, "--no-intro"])

    # 3. Dimensionalidade -- auto-detected by default, same convention as
    # stb-mlelastic's own interactive wrapper.
    dimensionality = get_input(
        "\nDimensionality, auto/3d/2d/1d [default: auto]: ").strip()
    if dimensionality:
        args.extend(["--dimensionality", dimensionality])

    # 4. Advanced settings (rarely-touched tolerances), gated -- same pattern
    # as stb-strain's own interactive wrapper.
    show_advanced = get_input(
        "\nConfigure advanced settings (vacuum-gap for auto-detection, physical "
        "cross-section for a 1D run)? [y/N]: ").strip().lower()
    if show_advanced == 'y':
        vacuum_gap = get_float_input(
            "Vacuum gap threshold in Ang, used only for auto-detection (default: 10.0): ", 10.0)
        args.extend(["--vacuum-gap", str(vacuum_gap)])
        cross_section = get_input(
            "Physical cross-section area in Ang^2, for a conventional GPa reading if this "
            "turns out to be 1D (blank = skip): ").strip()
        if cross_section:
            args.extend(["--cross-section", cross_section])

    # 5. Yield strength (opt-in -- borrowed metallurgy concept, off by default)
    show_yield = get_input(
        "\nAlso report the 0.2% offset Yield strength? (y/N): ").lower()
    if show_yield in ('y', 'yes'):
        args.append("--yield")

    # 6. Output directory + stb-standard report/plot options
    output_dir = get_input("\nOutput directory [default: current directory]: ").strip()
    if output_dir:
        args.extend(["--output-dir", output_dir])

    save_report = get_input("Also save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")
    save_gnuplot = get_input(
        "Also save the curve data + a gnuplot script? (y/N): ").strip().lower()
    if save_gnuplot == 'y':
        args.append("--save-gnuplot")
    view_choice = get_input(
        "Show an interactive matplotlib plot before finishing? (y/N): ").strip().lower()
    if view_choice == 'y':
        args.append("--view")

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

    output_dir = get_input("Output directory (default: '.'): ") or "."
    args.extend(["--output-dir", output_dir])

    save_report = get_input("\nAlso save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")

    run_tool("stb-bands", args)

def run_fatbands_analyzer() -> None:
    """Interface for the Fatbands Analyzer (stb-fatbands)"""
    print("\n" + "="*60)
    print(color_text("FATBANDS ANALYZER (stb-fatbands)", 'bold').center(60))
    print("="*60 + "\n")

    label = get_input("SIESTA label (e.g. 'siesta' for siesta.bands + siesta.bands.WFSX): ")
    while not os.path.isfile(f"{label}.bands"):
        print(color_text(f"File '{label}.bands' not found!", 'red'))
        label = get_input("SIESTA label: ")

    if os.path.isfile(f"{label}.bands.WFSX"):
        print(color_text(f"-> Found '{label}.bands.WFSX'.", 'green'))
    elif os.path.isfile(f"{label}.WFSX"):
        print(color_text(f"-> Found '{label}.WFSX'.", 'green'))
    elif os.path.isfile(f"{label}.selected.WFSX"):
        print(color_text(f"-> Found '{label}.selected.WFSX'.", 'green'))
    else:
        print(color_text(f"-> No '{label}.bands.WFSX'/'{label}.WFSX'/'{label}.selected.WFSX' "
                          "found -- stb-fatbands will fail without one. Re-run SIESTA with "
                          "'WFS.Write.For.Bands T' in the fdf.", 'yellow'))

    if os.path.isfile(f"{label}.HSX"):
        print(color_text(f"-> Found '{label}.HSX', overlap-aware (accurate) weights will be used.", 'green'))
    else:
        print(color_text(f"-> No '{label}.HSX' found -- weights will use an approximate "
                          "orthogonal-basis fallback.", 'yellow'))

    shift_options = {
        '1': ('fermi', "Fermi level"),
        '2': ('vbm', "Valence Band Maximum"),
        '3': ('cbm', "Conduction Band Minimum"),
        '4': ('manual', "Custom value")
    }

    print("\nEnergy reference options:")
    for key, (_, desc) in shift_options.items():
        print(f" {color_text(key, 'yellow')}. {desc}" + (" [Default]" if key == '1' else ""))

    choice = get_input("\nSelect reference (1-4) [default: 1]: ").strip() or '1'
    while choice not in shift_options:
        print(color_text("Invalid choice!", 'red'))
        choice = get_input("Select reference (1-4) [default: 1]: ").strip() or '1'

    shift_type, _ = shift_options[choice]
    args = ["--label", label, "--shift", shift_type, "--no-intro"]

    if shift_type == "manual":
        manual_value = get_float_input("Enter custom shift value: ")
        args.extend(["--manual-value", str(manual_value)])

    projection_options = {
        '1': "l",
        '2': "ml",
        '3': "atom",
        '4': "species",
        '5': "species_l",
    }
    print("\nOrbital projection options:")
    print(f" {color_text('1', 'yellow')}. l (s, p, d, f)")
    print(f" {color_text('2', 'yellow')}. ml (s, px, py, pz, dxy, ...)")
    print(f" {color_text('3', 'yellow')}. atom")
    print(f" {color_text('4', 'yellow')}. species")
    print(f" {color_text('5', 'yellow')}. species_l (species AND s/p/d/f combined, e.g. 'Sn-s', 'O-p') [Default]")
    proj_choice = get_input("Select projection (1-5) [default: 5]: ").strip() or '5'
    while proj_choice not in projection_options:
        print(color_text("Invalid choice!", 'red'))
        proj_choice = get_input("Select projection (1-5) [default: 5]: ").strip() or '5'
    args.extend(["--projection", projection_options[proj_choice]])

    output_dir = get_input("Output directory (default: '.'): ") or "."
    args.extend(["--output-dir", output_dir])

    save_report = get_input("\nAlso save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")
    save_gnuplot = get_input("Also save the data + gnuplot script? (y/N): ").strip().lower()
    if save_gnuplot == 'y':
        args.append("--save-gnuplot")
    view_choice = get_input("Show the matplotlib plot before finishing? (y/N): ").strip().lower()
    if view_choice == 'y':
        args.append("--view")

    run_tool("stb-fatbands", args)

def run_dos_convolution() -> None:
    """Interface for the DOS Processor (Convolution) (stb-convdos)"""
    print("\n" + "="*60)
    print(color_text("DOS PROCESSOR (CONVOLUTION) (stb-convdos)", 'bold').center(60))
    print("="*60 + "\n")

    print(f"{color_text('Process:', 'yellow')}")
    print(f"  {color_text('1.', 'yellow')} A single file [Default]")
    print(f"  {color_text('2.', 'yellow')} A whole folder, recursively (e.g. a stb-dos output folder)")
    mode_choice = get_input("\nSelect (1-2) [default: 1]: ").strip()

    args = []
    if mode_choice == '2':
        # Esta linha agora terá Tab-completion!
        input_dir = get_input("Input folder (e.g. the output of stb-dos): ")
        while not os.path.isdir(input_dir):
            print(color_text("Folder not found!", 'red'))
            input_dir = get_input("Input folder: ")
        output_dir = get_input(
            "Output folder [default: '<folder>_filtered']: "
        ).strip()
        args.extend(["--dir", input_dir])
        if output_dir:
            args.extend(["--output-dir", output_dir])
    else:
        # Esta linha agora terá Tab-completion!
        input_file = get_input("Input DOS file (e.g., dos_total.dat): ")
        while not os.path.isfile(input_file):
            print(color_text("File not found!", 'red'))
            input_file = get_input("Input DOS file: ")

        # Esta linha agora terá Tab-completion!
        out_file = get_input("Output file (default: dos_filtered.dat): ", 'green') or "dos_filtered.dat"
        args.extend(["--file", input_file, "--out", out_file])

    broadening_choice = get_input("\nSpecify broadening as (1) sigma or (2) FWHM, in meV [default: 1]: ").strip()
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
    args.extend([broadening_flag, str(broadening_value)])

    size_str = get_input("Kernel size in samples (optional, blank = auto-sized from broadening): ").strip()
    if size_str:
        args.extend(["--size", size_str])

    plot_prompt = ("Show before/after plot of the first processed file? (Y/n): " if mode_choice == '2'
                    else "Show before/after plot? (Y/n): ")
    plot_choice = get_input(plot_prompt).strip().lower()
    if plot_choice in ['n', 'no']:
        args.append("--no-plot")

    save_report = get_input("\nAlso save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")

    args.append("--no-intro")

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

    save_report = get_input("\nAlso save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")

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
        "Compare to an experimental pattern file (2theta first column, intensity LAST "
        "column -- e.g. stb-xrd's own xrd_pattern.dat works directly) [default: skip]: "
    ).strip()
    while compare_to and not os.path.isfile(compare_to):
        print(color_text("File not found!", 'red'))
        compare_to = get_input(
            "Compare to an experimental pattern file [default: skip]: ").strip()

    raw_experimental = False
    if compare_to:
        print(color_text(
            "\nNote: computing the similarity score can take ~20-30s "
            "(pyxtal's own comparison routine has no faster path).", 'yellow'))
        print("A sparse discrete peak list (not a continuous scan) is normally "
              "auto-detected and Gaussian-broadened before comparing, so it's on "
              "the same footing as the simulated pattern.")
        raw_choice = get_input(
            "Force comparing it exactly as-is instead, with no auto-broadening? "
            "(y/N): ").strip().lower()
        raw_experimental = raw_choice == 'y'

    args = ["--file", input_file, "--format", format_type, "--wavelength", wavelength,
            "--two-theta-range", str(two_theta_min), str(two_theta_max), "--no-intro"]
    if top_str:
        args.extend(["--top", top_str])
    if compare_to:
        args.extend(["--compare-to", compare_to])
        if raw_experimental:
            args.append("--raw-experimental")

    save_report = get_input("\nAlso save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")
    save_gnuplot = get_input("Also save the data + gnuplot script? (y/N): ").strip().lower()
    if save_gnuplot == 'y':
        args.append("--save-gnuplot")
    view_choice = get_input("Show the matplotlib preview before finishing? (y/N): ").strip().lower()
    if view_choice == 'y':
        args.append("--view")

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
    output_formats = ['cif', 'xyz', 'poscar', 'fdf', 'dftb', 'xsf', 'fhi', 'lammps']
    
    input_file = get_input("Input file path: ")
    while not os.path.isfile(input_file):
        print(color_text("File not found!", 'red'))
        input_file = get_input("Input file path: ")

    print(f"\n{color_text('Supported input formats:', 'yellow')}")
    print(f"  {color_text('0.', 'yellow')} Auto-detect from filename")
    for i, fmt in enumerate(input_formats, 1):
        print(f"  {color_text(str(i)+'.', 'yellow')} {fmt}")

    max_in = len(input_formats)
    choice_in = get_int_input(f"\nSelect input format (0-{max_in}) [default: 0]: ", 0)
    while not (0 <= choice_in <= max_in):
        print(color_text(f"Invalid choice! Please select between 0 and {max_in}.", 'red'))
        choice_in = get_int_input(f"Select input format (0-{max_in}) [default: 0]: ", 0)

    if choice_in == 0:
        from stb.translate import guess_format_from_filename
        in_format = guess_format_from_filename(input_file)
        if in_format is None:
            print(color_text(f"Could not auto-detect the format from '{input_file}'.", 'red'))
            while not (1 <= choice_in <= max_in):
                choice_in = get_int_input(f"Select input format (1-{max_in}): ")
                if not (1 <= choice_in <= max_in):
                    print(color_text(f"Invalid choice! Please select between 1 and {max_in}.", 'red'))
            in_format = input_formats[choice_in - 1]
            print(f"Selected input format: {color_text(in_format, 'cyan')}")
        else:
            print(f"Auto-detected input format: {color_text(in_format, 'cyan')}")
    else:
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
        print(color_text("\nXYZ format requires the 3 lattice vectors (Cartesian, in Angstrom).", 'yellow'))
        vec_values = []
        for vec_name in ("A", "B", "C"):
            for axis in ("x", "y", "z"):
                vec_values.append(get_float_input(f"Vector {vec_name}.{axis}: "))
        args.extend(["--lattice-vectors"] + [str(v) for v in vec_values])

    # ##### NOVO BLOCO: Visualização 3D opcional #####
    view_choice = get_input("\nView the structure in 3D? (y/N): ").strip().lower()
    if view_choice == "y":
        window_choice = get_input("  Open interactive window? (Y/n): ").strip().lower()
        if window_choice != "n":
            args.append("--view")
        image_path = get_input("  Save a snapshot image to (path, or Enter to skip): ").strip()
        if image_path:
            args.extend(["--view-image", image_path])
    # ##### FIM DO NOVO BLOCO #####

    # ##### NOVO BLOCO: Salvar relatório em arquivo (opcional) #####
    report_choice = get_input("\nSave the report to a file? (y/N): ").strip().lower()
    if report_choice == "y":
        args.append("--save-report")
    # ##### FIM DO NOVO BLOCO #####

    run_tool("stb-translate", args)




def run_clean_tool() -> None:
    """Interactive interface for the Clean Files tool (stb-clean)"""
    print("\n" + "="*60)
    print(color_text("CLEAN FILES TOOL (stb-clean)", 'bold').center(60))
    print("="*60 + "\n")
    print("Preps a SIESTA calc folder to restart a fresh run: removes prior run outputs")
    print("(logs, .out, .XV, .bands, .DOS, ...), keeps only the extensions below")
    print("(inputs/submission script by default). This deletes results too -- it's for")
    print("restarting, not archiving. Everything removed is backed up (tar.gz) first,")
    print("unless you skip that below.\n")

    # Esta linha agora terá Tab-completion!
    path = get_input("Directory (or glob pattern, e.g. strain_*) to clean (default: current): ").strip()
    if path == "":
        path = "."
    is_glob = any(c in path for c in '*?[')
    while not is_glob and not os.path.isdir(path):
        print(color_text("Directory not found!", 'red'))
        path = get_input("Enter a valid directory (or glob pattern): ").strip()
        is_glob = any(c in path for c in '*?[')

    default_exts = ['.psml', '.psf', '.fdf', '.sh']
    print(f"\nExtensions to keep (space-separated, default: {' '.join(default_exts)}):")
    ext_input = get_input("Extensions: ").strip()
    if ext_input:
        extensions = ext_input.split()
    else:
        extensions = default_exts

    recursive_choice = get_input("Clean subdirectories too (recursive)? [y/N]: ").strip().lower()
    recursive = recursive_choice == 'y'

    archive_choice = get_input("Create a backup archive before deleting? [Y/n]: ").strip().lower()
    no_archive = archive_choice == 'n'

    confirm_choice = get_input("Skip confirmation and delete directly? [y/N]: ").strip().lower()
    no_confirm = confirm_choice == 'y'

    report_choice = get_input("Save the report to a file? [y/N]: ").strip().lower()
    save_report = report_choice == 'y'

    args = ["--path", path, "--keep"] + extensions
    if recursive:
        args.append("--recursive")
    if no_archive:
        args.append("--no-archive")
    if no_confirm:
        args.append("--no-confirm")
    if save_report:
        args.append("--save-report")

    args.append("--no-intro")

    print()
    run_tool("stb-clean", args)

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

    save_report = get_input("\nAlso save a text report to file? (y/N): ").strip().lower()
    if save_report == 'y':
        args.append("--save-report")

    view_choice = get_input("View the structure interactively via ASE before finishing? (y/N): ").strip().lower()
    if view_choice == 'y':
        args.append("--view")

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
                         "substitute/subgroup/supergroup an existing structure.",
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
    10: {'title': "Fatbands Analyzer (stb-fatbands)",
        'description': "Orbital-projected / character-colored band structure from .bands + .WFSX.",
        'func': run_fatbands_analyzer},
    11: {'title': "STM Simulator (stb-stm)",
        'description': "Simulate an STM image (Tersoff-Hamann approximation) from a SIESTA .LDOS grid.",
        'func': run_stm_analyzer},
    12: {'title': "Wavefunction Density (stb-wfdensity)",
        'description': "Real-space |psi|^2 of a single (k, band) eigenstate from .WFSX -> Gaussian cube.",
        'func': run_wfdensity_analyzer},
    13: {'title': "STS Spectroscopy Simulator (stb-sts)",
        'description': "Simulated dI/dV(E) spectroscopy curve at a fixed tip point from a full-BZ .WFSX.",
        'func': run_sts_analyzer},
    14: {'title': "COOP/COHP Bonding Analysis (stb-coop)",
        'description': "Energy-resolved Crystal Orbital Overlap/Hamilton Population curves for selected atom pairs.",
        'func': run_coop_analyzer},
    15: {'title': "IPR Analyzer (stb-ipr)",
        'description': "Inverse Participation Ratio (localization measure) scattered onto a band structure.",
        'func': run_ipr_analyzer},
    16: {'title': "Effective Mass / Velocity (stb-effmass)",
        'description': "Effective mass tensor and band velocity for one chosen (k, band) eigenstate.",
        'func': run_effmass_analyzer},
    17: {'title': "Spin Texture Analyzer (stb-spintexture)",
        'description': "Spin moments (<Sx>,<Sy>,<Sz>) scattered onto a band structure from a non-collinear/SOC .WFSX.",
        'func': run_spintexture_analyzer},
    18: {'title': "AIMD Trajectory Analysis (stb-aimdAnalysis)",
        'description': "RDF, MSD/diffusion coefficient, and VACF-derived VDOS from an AIMD trajectory.",
        'func': run_aimd_analyzer},
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
        'description': "Generate displaced supercells via Phonopy, then post-process into thermal "
                        "properties. For a fast ML-only preview (no SIESTA) instead, including "
                        "QHA, see ML Simulations > ML Phonons (stb-mlphonons).",
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
    8: {'title': "Adsorption",
        'description': "Prepare a clean slab, an isolated adsorbate, and candidate adsorption "
                        "sites, then compute the real DFT adsorption energy per site.",
        'stages': {
            1: {'title': "Stage 1 - Prep (stb-adsorb)",
                'description': "Generate clean-slab/adsorbate/site folders ready for SIESTA.",
                'func': run_adsorb_setup},
            2: {'title': "Stage 2 - Analysis (stb-adsorbAnalysis)",
                'description': "Compute E_ads = E_site - E_clean_slab - E_adsorbate per site.",
                'func': run_adsorb_analysis},
        }},
    9: {'title': "NEB / Reaction Path",
        'description': "Interpolate a reaction-path band between two relaxed endpoints "
                        "(optionally refined with MACE-MP-0), then compute the DFT-level "
                        "energy profile and an approximate barrier.",
        'stages': {
            1: {'title': "Stage 1 - Prep (stb-neb)",
                'description': "Generate one single-point image_NN/ folder per band image.",
                'func': run_neb_setup},
            2: {'title': "Stage 2 - Analysis (stb-nebAnalysis)",
                'description': "Compute the energy profile and forward/backward barrier.",
                'func': run_neb_analysis},
        }},
    10: {'title': "2D Stacking Fault (Gamma-Surface)",
        'description': "Slide one 2D layer across a grid of lateral offsets, then compute the "
                        "equilibrium stacking, corrugation energy, and full gamma-surface map.",
        'stages': {
            1: {'title': "Stage 1 - Prep (stb-stackingfault)",
                'description': "Generate one single-point shift_II_JJ/ folder per grid point.",
                'func': run_stackingfault_setup},
            2: {'title': "Stage 2 - Analysis (stb-stackingfaultAnalysis)",
                'description': "Compute the equilibrium stacking and corrugation energy.",
                'func': run_stackingfault_analysis},
        }},
    11: {'title': "Raman Spectrum",
        'description': "Self-contained: computes its own phonons, then the +/-delta Optical "
                        "(dielectric-response) displacements needed for the diagonal Raman "
                        "tensor, then the approximate Raman-active frequencies and spectrum.",
        'stages': {
            1: {'title': "Stage 1 - Phonon Displacements (stb-raman)",
                'description': "Generate one phonon_disp/disp-NNN/ folder per finite-difference "
                                "displacement.",
                'func': run_raman_prep},
            2: {'title': "Stage 2 - Modes & Optical Displacements (stb-ramanModes)",
                'description': "Build FORCE_SETS, get the Gamma-point modes, and generate one "
                                "optical_disp/mode_*/ folder per +/-delta x/y/z displacement.",
                'func': run_raman_modes},
            3: {'title': "Stage 3 - Analysis (stb-ramanAnalysis)",
                'description': "Compute the diagonal Raman tensor per mode and the approximate "
                                "Raman spectrum.",
                'func': run_raman_analysis},
        }},
    12: {'title': "IR Spectrum",
        'description': "Self-contained: computes its own phonons, then either the +/-delta "
                        "dipole displacements (0D/1D/2D) or a single equilibrium "
                        "Born-effective-charge run (3D bulk), then the IR-active frequencies "
                        "and spectrum. Path auto-selected by structure dimensionality.",
        'stages': {
            1: {'title': "Stage 1 - Phonon Displacements (stb-ir)",
                'description': "Generate one phonon_disp/disp-NNN/ folder per finite-difference "
                                "displacement.",
                'func': run_ir_prep},
            2: {'title': "Stage 2 - Modes & Dipole/Born-Charge Displacements (stb-irModes)",
                'description': "Build FORCE_SETS, get the Gamma-point modes, and generate "
                                "either one dipole_disp/mode_*/ folder pair per +/-delta mode "
                                "(non-bulk) or a single born_charge_disp/equilibrium/ folder "
                                "(bulk).",
                'func': run_ir_modes},
            3: {'title': "Stage 3 - Analysis (stb-irAnalysis)",
                'description': "Compute dmu/dQ per mode (dipole finite-difference or "
                                "Born-charge x eigendisplacement projection) and the IR "
                                "spectrum.",
                'func': run_ir_analysis},
        }},
    13: {'title': "Hydrogen Evolution Reaction (HER)",
        'description': "Self-contained: finds the most stable H-adsorption site on a slab, "
                        "then the H2/BSSE/ZPE references needed for the computational hydrogen "
                        "electrode descriptor Delta-G_H* (Norskov et al.).",
        'stages': {
            1: {'title': "Stage 1 - Adsorption Sites (stb-her)",
                'description': "Find every symmetrically distinct H-adsorption site and write "
                                "one relaxation folder per site.",
                'func': run_her_prep},
            2: {'title': "Stage 2 - References & ZPE Prep (stb-herRefs)",
                'description': "Pick the winning site, build the H2/deformed-slab/BSSE-ghost "
                                "references, and the ZPE calculation folder(s).",
                'func': run_her_refs},
            3: {'title': "Stage 3 - Analysis (stb-herAnalysis)",
                'description': "Combine every reference energy into the final Delta-G_H*.",
                'func': run_her_analysis},
        }},
    14: {'title': "Oxygen Evolution Reaction (OER)",
        'description': "Self-contained: finds the most stable OH-adsorption site on a slab, derives "
                        "(or independently searches for) the O*/OOH* intermediates from it, then the "
                        "H2/H2O/BSSE/ZPE references needed for the 4-electron CHE descriptor "
                        "(Delta-G1..4, overpotential eta, potential-determining step; Rossmeisl et "
                        "al. 2007, Man et al. 2011).",
        'stages': {
            1: {'title': "Stage 1 - Adsorption Sites (stb-oer)",
                'description': "Find every symmetrically distinct OH-adsorption site and write "
                                "one relaxation folder per site.",
                'func': run_oer_prep},
            2: {'title': "Stage 2 - O*/OOH* Intermediates (stb-oerIntermediates)",
                'description': "Derive (or independently search for) O*/OOH* from the winning OH* "
                                "site and write their own CG-relaxation folder(s).",
                'func': run_oer_intermediates},
            3: {'title': "Stage 3 - References, BSSE & ZPE Prep (stb-oerRefs)",
                'description': "Build the H2/H2O/BSSE/ZPE reference folders from the final relaxed "
                                "OH*/O*/OOH* geometries.",
                'func': run_oer_refs},
            4: {'title': "Stage 4 - Analysis (stb-oerAnalysis)",
                'description': "Combine every reference energy into Delta-G1..4, eta, and the PDS.",
                'func': run_oer_analysis},
        }},
    15: {'title': "Generalized Quasi-Chemical Approximation (GQCA)",
        'description': "Self-contained: builds the 3 ordered pair-cluster-representative "
                        "structures (AA/AB/BB) for a substitutional alloy A(1-x)Bx, then solves "
                        "the mass-action equilibrium to get Delta-H_mix/Delta-S_mix/Delta-G_mix"
                        "(x,T) and short-range-order parameters (Sher et al., PRB 36, 4279, "
                        "1987). Works on 2D and 3D structures alike.",
        'stages': {
            1: {'title': "Stage 1 - Structure Generation (stb-gqca)",
                'description': "Generate the 3 ordered pair-representative structures (AA/AB/BB), "
                                "fully relaxed, ready for SIESTA.",
                'func': run_gqca_prep},
            2: {'title': "Stage 2 - Analysis (stb-gqcaAnalysis)",
                'description': "Solve the mass-action equilibrium over an x/T grid and report "
                                "Delta-H_mix, Delta-S_mix, Delta-G_mix, and SRO parameters.",
                'func': run_gqca_analysis},
        }},
    16: {'title': "Optical Properties",
        'description': "Self-contained: writes one SIESTA folder per Cartesian polarization "
                        "direction for a linear-response Optical.Calculate run on an already-"
                        "relaxed structure, then derives eps1(E)/eps2(E), n(E)/k(E), the "
                        "absorption coefficient, normal-incidence reflectivity, the energy-loss "
                        "function, and the real optical conductivity from the resulting "
                        "SystemLabel.EPSIMG file(s).",
        'stages': {
            1: {'title': "Stage 1 - Direction Folders (stb-optical)",
                'description': "Write one direction folder per requested Optical.Vector axis, "
                                "ready for SIESTA.",
                'func': run_optical_prep},
            2: {'title': "Stage 2 - Analysis (stb-opticalAnalysis)",
                'description': "Aggregate every direction folder present and derive the full "
                                "set of linear optical properties.",
                'func': run_optical_analysis},
        }},
    17: {'title': "Custom ML Force Field (MACE fine-tuning)",
        'description': "Generate training configurations (rattled and/or AIMD-sampled), then "
                        "fine-tune a MACE-MP foundation model on the resulting SIESTA "
                        "energy/forces/stress -- a fast potential specialized to your own "
                        "material instead of the generic foundation model. Optionally exports "
                        "a LAMMPS-loadable model.",
        'stages': {
            1: {'title': "Stage 1 - Prep (stb-mlff)",
                'description': "Generate rattled/AIMD-sampled training configurations, ready "
                                "for SIESTA.",
                'func': run_mlff_prep},
            2: {'title': "Stage 2 - Analysis (stb-mlffAnalysis)",
                'description': "Aggregate the finished SIESTA calculations and fine-tune a "
                                "MACE-MP foundation model on them.",
                'func': run_mlff_analysis},
            3: {'title': "Stage 3 - Active Learning (stb-mlffActiveLearning)",
                'description': "Use the fine-tuned model to suggest new configurations to "
                                "label with SIESTA (largest predicted force), for another "
                                "round of refinement.",
                'func': run_mlff_active_learning},
        }},
    18: {'title': "Equation of State",
        'description': "Generate isotropically-scaled-volume structures, then fit an E vs V "
                        "equation of state (V0/E0/B0/B0') from the resulting SIESTA "
                        "energies. For a fast ML-only preview (no SIESTA) instead, see ML "
                        "Simulations > ML Equation of State (stb-mleos).",
        'stages': {
            1: {'title': "Stage 1 - Prep (stb-eosInputs)",
                'description': "Generate isotropically-scaled-volume structures for a real-SIESTA EOS scan.",
                'func': run_eos_inputs_generator},
            2: {'title': "Stage 2 - Analysis (stb-eosAnalysis)",
                'description': "Aggregate the vol_* folders and fit an equation of state (bulk modulus, etc).",
                'func': run_eos_analysis_generator},
        }},
       }


# Tools that RUN a simulation using a trained MACE potential (the MACE-MP-0
# foundation model, or a custom one fine-tuned via stb-mlffAnalysis) instead
# of driving a real SIESTA calculation -- distinct from WORKFLOW_TOOLS (which
# is about generating/analyzing DFT calculations), and from stb-mlrelax/
# stb-amorphize/etc. in STRUCTURE_TOOLS (one-shot structure preprocessing).
# Category 5 in the main menu, between Workflow (4) and Utils (6).
MLSIM_TOOLS = {
    1: {'title': "ML Molecular Dynamics (stb-mlmd)",
        'description': "NVE/NVT/NPT molecular dynamics driven by a MACE potential -- fast "
                        "exploration with the foundation model or your own fine-tuned one.",
        'func': run_mlmd_generator},
    2: {'title': "ML Phonons (stb-mlphonons)",
        'description': "Standalone phonon bands/DOS/thermal properties driven entirely by a "
                        "MACE potential -- no SIESTA, no separate analysis tool needed.",
        'func': run_mlphonons_generator},
    3: {'title': "ML Elastic Constants (stb-mlelastic)",
        'description': "Standalone stiffness matrix (stress-strain method) driven entirely by "
                        "a MACE potential -- no SIESTA, no separate strain_*/ folders needed.",
        'func': run_mlelastic_generator},
    4: {'title': "ML Structure Search (stb-mlsearch)",
        'description': "Basin-hopping global search (perturb + relax + accept/reject) for the "
                        "lowest-energy atomic arrangement at a fixed cell, driven by MACE.",
        'func': run_mlsearch_generator},
    5: {'title': "ML Melting Point (stb-mlmelting)",
        'description': "Brackets a bulk material's melting point via a sequence of short MACE "
                        "MD runs, tracking the Lindemann index vs. temperature.",
        'func': run_mlmelting_generator},
    6: {'title': "ML Model-Size Convergence (stb-mlconvergence)",
        'description': "Compares MACE-MP-0 model sizes (small/medium/large) on the same "
                        "reference calculation -- helps pick a size instead of guessing.",
        'func': run_mlconvergence_generator},
    7: {'title': "ML NEB (stb-mlneb)",
        'description': "Standalone climbing-image NEB (reaction barrier estimate) driven "
                        "entirely by a MACE potential -- no SIESTA, no calc.fdf needed.",
        'func': run_mlneb_generator},
    8: {'title': "ML Vacancy Diffusion (stb-mldiffusion)",
        'description': "Auto-detects neighbor hop targets for a vacancy and screens the "
                        "migration barrier for each via MACE-driven NEB.",
        'func': run_mldiffusion_generator},
    9: {'title': "ML Monte Carlo / GCMC (stb-mlgcmc)",
        'description': "Canonical/grand-canonical Monte Carlo adsorption (insert/move/"
                        "delete a single-atom adsorbate) driven by a MACE potential.",
        'func': run_mlgcmc_generator},
    10: {'title': "ML Adsorption-Site Screening (stb-mladsorb)",
        'description': "Enumerates candidate ontop/bridge/hollow sites, relaxes each with "
                        "MACE, and ranks them by a genuine ML adsorption energy.",
        'func': run_mladsorb_generator},
    11: {'title': "ML Equation of State (stb-mleos)",
        'description': "E vs V curve + equation-of-state fit (V0/E0/B0/B0') driven by MACE -- "
                        "an independent bulk-modulus cross-check against stb-mlelastic.",
        'func': run_mleos_generator},
}


UTILITY_TOOLS = {
    1: {'title': "File Translator (stb-translate)",
        'description': "Convert between file formats (CIF, POSCAR, fdf, xyz...).",
        'func': run_file_translator},
    2: {'title': "Clean File Tools (stb-clean)",
        'description': "Prep a calc folder to restart a fresh run (removes prior outputs, "
                        "keeps only inputs/submission script).",
        'func': run_clean_tool},
    3: {'title': "Grid to Cube Converter", 'description': "Convert Siesta Grid files (VT, RHO, VH) to Gaussian .cube.", 'func': run_grid_to_cube},
    4: {'title': "Wantibexos Interface (stb-siesta2wtb)",
        'description': "Convert SIESTA Hamiltonian to Wantibexos format.",
        'func': run_wantibexos_interface},
    5: {'title': "AIMD Trajectory Converter (stb-ani2traj)",
        'description': "Convert a SIESTA AIMD .ANI trajectory (no lattice of its own) into "
                        "a multi-frame format with the cell, for OVITO/VMD.",
        'func': run_ani2traj_converter},
    6: {'title': "Status Checker (stb-status)",
        'description': "Quick summary of a calc folder: run type, SCF convergence, max "
                        "force, final energy, and which key output files are present.",
        'func': run_status_checker},
    7: {'title': "Archive Packager (stb-archive)",
        'description': "Package a finished calc (inputs + essential outputs) into a "
                        ".tar.gz for sharing/reproducibility -- the complement of stb-clean.",
        'func': run_archive_packager},
}


def _flatten_tool_codes() -> Dict[str, Callable]:
    """Builds a flat {"1.1": func, ..., "2.1": func, ..., "4.1.2": func, ...,
    "5.1": func, ..., "6.4": func} lookup across all 6 categories, from the
    dicts above -- lets the main menu jump straight to a tool via a dotted
    code instead of navigating level by level.
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
    for key, info in MLSIM_TOOLS.items():
        codes[f"5.{key}"] = info['func']
    for key, info in UTILITY_TOOLS.items():
        codes[f"6.{key}"] = info['func']
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
            choice = get_input("\nSelect an option (0-6, or a tool code like 4.1.2): ")

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
                run_sub_menu("ML Simulations", MLSIM_TOOLS)
            elif choice == '6':
                run_sub_menu("Utils", UTILITY_TOOLS)
            elif choice == '0':
                print(color_text("\nThank you for using STB-SUITE!", 'cyan'))
                break
            else:
                print(color_text("\nInvalid choice! Please select between 0 and 6, or a valid tool code.", 'red'))
                sleep(1)
                
        except ValueError:
            print(color_text("\nPlease enter a valid number!", 'red'))
            sleep(1)
        except KeyboardInterrupt:
            print(color_text("\n\nOperation cancelled by user.", 'yellow'))
            break

if __name__ == "__main__":
    main()
