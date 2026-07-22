#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################
    
VERSION = "1.11.0" # SCF/force quality diagnostics, multi-site BSSE averaging, atoms_bsse_check/ convergence report

import os
import re
import sys
import argparse
from stb.core import siesta_log, structure_io
from stb.core.cli import color_text, show_intro

def get_atom_counts(fdf_path):
    """Number of atoms per chemical species -- thin wrapper around the shared
    core/structure_io.py::atom_counts (moved there once cohesive_energy.py's
    prep side needed the identical logic to decide which isolated-atom
    folders are actually needed)."""
    try:
        structure = structure_io.read_fdf(fdf_path)
    except FileNotFoundError:
        print(color_text(f"[ERROR] Structure file '{fdf_path}' not found.", 'red'))
        sys.exit(1)
    except ValueError as e:
        print(color_text(f"[ERROR] {e}", 'red'))
        sys.exit(1)

    return structure_io.atom_counts(structure)


_SITE_DIR_RE = re.compile(r"^site_.+_x(\d+)$")


def read_bsse_energy(bsse_root, sym, out_file):
    """Reads the BSSE-corrected isolated-atom energy for species `sym` under
    `bsse_root` ('atoms_bsse' or 'atoms_bsse_check'), transparently handling
    both layouts stb-cohesive can produce (see cohesive_energy.py::
    resolve_bsse_sites): a flat 'bsse_root/<sym>/<out_file>' (single site, or
    --no-bsse-multi-site), or a multi-site 'bsse_root/<sym>/site_<wyckoff>_x
    <mult>/<out_file>' per symmetrically distinct site -- averaged, weighted
    by multiplicity, so a species with e.g. 4 octahedral + 2 tetrahedral
    atoms gets a single reference that's actually representative of its real
    population, not just whichever site happened to be generated first.

    Returns (energy, complete): `complete` is False if the flat file, or ANY
    expected multi-site subfolder, is missing/unreadable -- a partial
    weighted average would silently misrepresent the correction, so this
    always drops the whole species rather than average over what's present
    (same "advisory only, never silently wrong" spirit as the rest of the
    suite; the caller decides whether to skip BSSE reporting entirely).
    """
    sym_dir = os.path.join(bsse_root, sym)
    flat_path = os.path.join(sym_dir, out_file)
    if os.path.isfile(flat_path):
        energy = siesta_log.get_free_energy(flat_path)
        return energy, energy is not None

    if not os.path.isdir(sym_dir):
        return None, False

    site_dirs = sorted(
        d for d in os.listdir(sym_dir)
        if _SITE_DIR_RE.match(d) and os.path.isdir(os.path.join(sym_dir, d))
    )
    if not site_dirs:
        return None, False

    total_weight = 0
    weighted_sum = 0.0
    for d in site_dirs:
        mult = int(_SITE_DIR_RE.match(d).group(1))
        e = siesta_log.get_free_energy(os.path.join(sym_dir, d, out_file))
        if e is None:
            return None, False
        weighted_sum += e * mult
        total_weight += mult
    if total_weight == 0:
        return None, False
    return weighted_sum / total_weight, True

def main():
    parser = argparse.ArgumentParser(
        description="Process cohesive energy results from SIESTA calculations.",
        epilog="Example usage:\n"
               "  stb_cohesive_analysis -o calc.out\n"
               "  stb_cohesive_analysis -o calc.out -d /path/to/results",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("-o", "--out", dest="out_file", type=str, required=True, 
                        help="Name of the SIESTA output file (e.g., calc.out)")
    
    parser.add_argument("-d", "--dir", dest="dir_path", type=str, default=".",
                        help="Path to the folder containing 'structure' and 'atoms' directories (default: current directory)")

    parser.add_argument("--force-tolerance", dest="force_tolerance", type=float, default=0.05,
                        help="Residual atomic force in eV/Ang (default: 0.05) above which the "
                             "full structure's calc.out is flagged as possibly not relaxed -- "
                             "the cohesive energy would then reflect a strained/off-equilibrium "
                             "geometry rather than the true minimum. Advisory only, never "
                             "blocks the result.")

    parser.add_argument("--save-report", action="store_true",
                        help="Also persist the results to cohesive_results.dat. Off by default.")

    parser.add_argument("-v", "--version", action="version", version=f"stb-cohesiveAnalysis {VERSION}")

    parser.add_argument("--no-intro", dest="intro", action="store_false",
                        help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro == True:
        show_intro([
            "Siesta ToolBox Suite",
            "Cohesive Energy Post-Processing",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("COHESIVE ANALYSIS:", 'bold'))
    print("-" * 60)
    
    root_dir = args.dir_path
    struct_dir = os.path.join(root_dir, "structure")
    atoms_dir = os.path.join(root_dir, "atoms")
    atoms_bsse_dir = os.path.join(root_dir, "atoms_bsse")
    atoms_bsse_check_dir = os.path.join(root_dir, "atoms_bsse_check")

    print("\n[INFO] Checking directories ...")
    if not os.path.exists(struct_dir) or not os.path.exists(atoms_dir):
        print(color_text(f"[ERROR] Required directories 'structure' and/or 'atoms' not found in '{root_dir}'.", 'red'))
        sys.exit(1)

    # atoms_bsse/ and atoms_bsse_check/ (written by stb-cohesive
    # --bsse-correction/--bsse-convergence-check, the former the default)
    # are auto-detected -- no flag needed here, and their absence never
    # blocks the uncorrected result.
    bsse_available = os.path.isdir(atoms_bsse_dir)
    if bsse_available:
        print(f"{color_text('[INFO]', 'cyan')} Detected '{atoms_bsse_dir}' -- will also report a "
              "BSSE (counterpoise) -corrected cohesive energy.")
    check_available = os.path.isdir(atoms_bsse_check_dir)
    if check_available:
        print(f"{color_text('[INFO]', 'cyan')} Detected '{atoms_bsse_check_dir}' -- will also "
              "report a BSSE cutoff convergence check.")

    # 1. Parse Structure for Atom Counts
    struct_fdf = os.path.join(struct_dir, "structure.fdf")
    print("[INFO] Read structure file ...")
    atom_counts = get_atom_counts(struct_fdf)

    if not atom_counts:
        print(color_text("[ERROR] Could not extract atom counts. Check your structure.fdf", 'red'))
        sys.exit(1)

    total_atoms = sum(atom_counts.values())
    if total_atoms == 0:
        print(color_text("[ERROR] Structure file declares species but has zero atoms. Check your structure.fdf", 'red'))
        sys.exit(1)
    print(f"[INFO] Total atoms in cell: {total_atoms}")
    for sym, count in atom_counts.items():
        note = "" if count > 0 else "  (declared but never placed -- no isolated-atom calculation needed)"
        print(f"       -> {sym}: {count} atoms{note}")

    # Species that actually need an isolated-atom reference energy -- a
    # species declared in ChemicalSpeciesLabel but with 0 real atoms always
    # contributes exactly 0 to the cohesive-energy sum regardless of its own
    # energy, so it must never be required/block the analysis (see
    # core/structure_io.py::atom_counts docstring).
    used_species = {sym: count for sym, count in atom_counts.items() if count > 0}

    # 2. Extract Energies
    print("\n[INFO] Extracting energies ...")
    errors = False

    # Bulk Energy
    struct_out_path = os.path.join(struct_dir, args.out_file)
    e_bulk = siesta_log.get_free_energy(struct_out_path)
    if e_bulk is None:
        print(color_text(f"[WARNING] Could not find '{args.out_file}' or finished calculation for the full structure.", 'yellow'))
        errors = True
    else:
        print(f"[INFO] Energy (Full Structure): {e_bulk:12.6f} eV")

        # Numerical-quality diagnostics (best-effort, never blocks): SCF
        # convergence and residual forces, checked ONLY on the full
        # structure -- an isolated atom has zero net force by symmetry, and
        # the BSSE cluster's real atom only feels Pulay-type forces from the
        # ghosts' incomplete basis (not a physical force to "relax" away),
        # so neither is a meaningful place to ask "is this geometry actually
        # at its minimum".
        converged, iterations = siesta_log.get_scf_convergence(struct_out_path)
        if not converged:
            print(color_text(
                f"[WARNING] Could not confirm SCF convergence for the full structure "
                f"('{args.out_file}') -- this result may be unreliable.", 'yellow'))
        else:
            print(f"[INFO] SCF converged after {iterations} iteration(s).")

        max_force = siesta_log.get_max_force(struct_out_path)
        if max_force is not None:
            if max_force > args.force_tolerance:
                print(color_text(
                    f"[WARNING] Residual force on the full structure ({max_force:.4f} eV/Ang) "
                    f"exceeds --force-tolerance ({args.force_tolerance} eV/Ang) -- this geometry "
                    "may not be relaxed, so the cohesive energy may reflect a strained/off-"
                    "equilibrium structure rather than the true minimum.", 'yellow'))
            else:
                print(f"[INFO] Residual force on the full structure: {max_force:.4f} eV/Ang "
                      f"(within --force-tolerance {args.force_tolerance} eV/Ang).")

    # Isolated Atoms Energy (uncorrected)
    isolated_energies = {}
    sum_isolated_energy = 0.0

    for sym, count in used_species.items():
        sym_dir = os.path.join(atoms_dir, sym)
        e_atom = siesta_log.get_free_energy(os.path.join(sym_dir, args.out_file))

        if e_atom is None:
            print(color_text(f"[WARNING] Could not find '{args.out_file}' or results for isolated atom: {sym}", 'yellow'))
            errors = True
        else:
            isolated_energies[sym] = e_atom
            sum_isolated_energy += (e_atom * count)
            print(f"[INFO] Energy (Isolated {sym}):   {e_atom:12.6f} eV")

    # 2b. BSSE-corrected (counterpoise) isolated-atom energies, if the prep
    # stage generated them -- best-effort, never blocks the uncorrected
    # result above: any missing/incomplete species just drops BSSE-corrected
    # reporting entirely (same "advisory only" spirit as the rest of the
    # suite's numerical-quality diagnostics). read_bsse_energy transparently
    # handles both the flat (single site) and multi-site (weighted-average)
    # layouts stb-cohesive --bsse-multi-site can produce.
    bsse_energies = {}
    bsse_complete = bsse_available
    if bsse_available:
        for sym in used_species:
            e_atom_bsse, ok = read_bsse_energy(atoms_bsse_dir, sym, args.out_file)
            if not ok:
                print(color_text(
                    f"[WARNING] Could not find '{args.out_file}' or results for the BSSE-corrected "
                    f"isolated atom: {sym} -- BSSE-corrected cohesive energy will be skipped.", 'yellow'))
                bsse_complete = False
            else:
                bsse_energies[sym] = e_atom_bsse
                print(f"[INFO] Energy (Isolated {sym}, BSSE-corrected): {e_atom_bsse:12.6f} eV")

    # 2c. BSSE cutoff convergence check, if the prep stage generated it
    # (stb-cohesive --bsse-convergence-check) -- same best-effort spirit,
    # purely informational (no pass/fail threshold: how much residual is
    # acceptable is application-dependent).
    check_energies = {}
    check_complete = check_available
    if check_available:
        for sym in used_species:
            e_atom_check, ok = read_bsse_energy(atoms_bsse_check_dir, sym, args.out_file)
            if not ok:
                print(color_text(
                    f"[WARNING] Could not find '{args.out_file}' or results for the BSSE "
                    f"convergence-check isolated atom: {sym} -- the convergence check will be "
                    "skipped.", 'yellow'))
                check_complete = False
            else:
                check_energies[sym] = e_atom_check
                print(f"[INFO] Energy (Isolated {sym}, BSSE check): {e_atom_check:12.6f} eV")

    # 3. Calculate Cohesive Energy
    if errors:
        print(color_text("\n[ERROR] Cannot calculate cohesive energy because some calculations are missing or incomplete.", 'red'))
        sys.exit(1)

    # Cohesive Energy Formula: E_coh = (E_bulk - Sum(E_iso)) / N_atoms
    # Negative value means a bound, stable structure.
    print("\n[INFO] Calculating Cohesive Energy ...")
    e_coh_total =  e_bulk - sum_isolated_energy
    e_coh_per_atom = e_coh_total / total_atoms

    e_coh_total_bsse = None
    e_coh_per_atom_bsse = None
    if bsse_available and bsse_complete:
        sum_isolated_bsse = sum(bsse_energies[sym] * count for sym, count in used_species.items())
        e_coh_total_bsse = e_bulk - sum_isolated_bsse
        e_coh_per_atom_bsse = e_coh_total_bsse / total_atoms

    e_coh_per_atom_check = None
    if check_available and check_complete and bsse_available and bsse_complete:
        sum_isolated_check = sum(check_energies[sym] * count for sym, count in used_species.items())
        e_coh_total_check = e_bulk - sum_isolated_check
        e_coh_per_atom_check = e_coh_total_check / total_atoms

    print("\n[INFO] FINAL RESULTS:")
    print(f"       Sum of Isolated Atoms: {sum_isolated_energy:12.6f} eV")
    print(f"       Bulk Structure Energy: {e_bulk:12.6f} eV")
    print(f"       Total Cohesive Energy (uncorrected):     {e_coh_total:10.4f} eV")
    print(f"       Cohesive Energy per Atom (uncorrected):  {e_coh_per_atom:10.4f} eV/atom")
    if e_coh_per_atom_bsse is not None:
        delta_per_atom = e_coh_per_atom_bsse - e_coh_per_atom
        print(f"       Total Cohesive Energy (BSSE-corrected):     {e_coh_total_bsse:10.4f} eV")
        print(f"       Cohesive Energy per Atom (BSSE-corrected):  {e_coh_per_atom_bsse:10.4f} eV/atom")
        print(f"       BSSE correction per atom: {delta_per_atom:+.4f} eV/atom (uncorrected LCAO "
              "cohesive energies systematically over-bind -- expect this to make the energy "
              "less negative)")
    elif not bsse_available:
        print(color_text(
            "\n[NOTE] This result is NOT corrected for BSSE (Basis Set Superposition Error) -- "
            "a known bias of LCAO cohesive energies: the isolated atom is computed with a "
            "smaller effective basis than it has inside the solid (no neighbors' orbitals "
            "available), which systematically over-binds. Re-run stb-cohesive with "
            "--bsse-correction (the CLI default) for a BSSE-corrected reference.", 'yellow'))

    if e_coh_per_atom_check is not None:
        check_delta = e_coh_per_atom_check - e_coh_per_atom_bsse
        print(f"       Cohesive Energy per Atom (BSSE check, larger cutoff): {e_coh_per_atom_check:10.4f} eV/atom")
        print(f"       BSSE cutoff convergence shift: {check_delta:+.4f} eV/atom (difference "
              "between --bsse-cutoff and --bsse-cutoff + --bsse-convergence-increment -- small "
              "means the BSSE correction is converged w.r.t. cutoff; large means increase "
              "--bsse-cutoff further)")

    # 4. Save results to a .dat file (opt-in)
    if args.save_report:
        print("\n[INFO] Write files ...")
        out_file_path = os.path.join(root_dir, "cohesive_results.dat")
        try:
            with open(out_file_path, 'w') as f_out:
                f_out.write("==================================================\n")
                f_out.write("             COHESIVE ENERGY RESULTS              \n")
                f_out.write("==================================================\n\n")

                f_out.write(f"Total atoms in cell: {total_atoms}\n")
                for sym, count in atom_counts.items():
                    note = "" if count > 0 else "  (declared but never placed)"
                    f_out.write(f"  -> {sym}: {count} atoms{note}\n")
                f_out.write("-" * 50 + "\n")

                f_out.write(f"Energy (Full Structure): {e_bulk:16.6f} eV\n")
                for sym in used_species:
                    f_out.write(f"Energy (Isolated {sym}):   {isolated_energies[sym]:16.6f} eV\n")
                f_out.write("-" * 50 + "\n")

                f_out.write(f"Sum of Isolated Atoms:   {sum_isolated_energy:16.6f} eV\n\n")
                f_out.write(f"Total Cohesive Energy (uncorrected):       {e_coh_total:12.4f} eV\n")
                f_out.write(f"Cohesive Energy per Atom (uncorrected):    {e_coh_per_atom:12.4f} eV/atom\n")
                if e_coh_per_atom_bsse is not None:
                    f_out.write("\n" + "-" * 50 + "\n")
                    f_out.write("BSSE (counterpoise) -corrected reference:\n")
                    for sym in used_species:
                        f_out.write(f"Energy (Isolated {sym}, BSSE-corrected):   {bsse_energies[sym]:16.6f} eV\n")
                    f_out.write(f"Total Cohesive Energy (BSSE-corrected):    {e_coh_total_bsse:12.4f} eV\n")
                    f_out.write(f"Cohesive Energy per Atom (BSSE-corrected): {e_coh_per_atom_bsse:12.4f} eV/atom\n")
                    if e_coh_per_atom_check is not None:
                        f_out.write(f"Cohesive Energy per Atom (BSSE check, larger cutoff): "
                                     f"{e_coh_per_atom_check:12.4f} eV/atom\n")
                        f_out.write(f"BSSE cutoff convergence shift: "
                                     f"{e_coh_per_atom_check - e_coh_per_atom_bsse:+.4f} eV/atom\n")
                elif not bsse_available:
                    f_out.write("\nNOTE: not corrected for BSSE (Basis Set Superposition Error) -- "
                                 "re-run stb-cohesive with --bsse-correction for a corrected reference.\n")
                f_out.write("\n==================================================\n")

            print(f"[INFO] Results saved to: {out_file_path}")
        except Exception as e:
            print(color_text(f"[ERROR] Failed to save results to file: {e}", 'red'))
            sys.exit(1)

    print("\n[INFO] Complete job!") 
    print("\n"+"-"*60)
    print(color_text("Cohesive energy analysis complete!\n\n", 'bold'))

if __name__ == "__main__":
    main()
