#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.0.0"

import os
import re
import sys
import argparse
import numpy as np
from stb.core import siesta_log, structure_io, gqca_solver
from stb.core.cli import color_text, show_intro, print_dual

REPORT_FILE = "gqca_stage2.txt"

_SUBLATTICE_RE = re.compile(r"^Sublattice\s*:\s*(\S+)\s*$", re.MULTILINE)
_SPECIES_B_RE = re.compile(r"^Species B\s*:\s*(\S+)\s*$", re.MULTILINE)
_AB_DEVIATION_RE = re.compile(r"^AB deviation\s*:\s*(\S+)\s*$", re.MULTILINE)


def read_stage1_metadata(directory):
    """Reads back --sublattice/--species-b/the achieved AB pair-correlation
    deviation from Stage 1's own persisted report (gqca.py::main appends
    these lines at the very end of gqca_stage1.txt) -- same "each stage
    re-reads/re-derives, never imports a sibling stage's module"
    convention as every other WORKFLOW_TOOLS pair/trio in this suite.
    Returns (sublattice, species_b, ab_deviation) or (None, None, None)
    if the report is missing/unparsable. `ab_deviation` is 0.0 (not None)
    for reports predating this field (backward-compatible silent
    default), since a missing value here is far more likely an older
    report than a real deviation.
    """
    report_path = os.path.join(directory, "gqca_stage1.txt")
    if not os.path.isfile(report_path):
        return None, None, None
    with open(report_path, 'r', encoding='utf-8', errors='ignore') as f:
        text = f.read()
    sub_match = _SUBLATTICE_RE.search(text)
    b_match = _SPECIES_B_RE.search(text)
    if not sub_match or not b_match:
        return None, None, None
    dev_match = _AB_DEVIATION_RE.search(text)
    ab_deviation = float(dev_match.group(1)) if dev_match else 0.0
    return sub_match.group(1), b_match.group(1), ab_deviation


def read_cluster_energy(cluster_dir, out_file, sublattice, species_b):
    """Reads one cluster_n{j} folder's normalized per-site energy.

    [PER-SITE NORMALIZATION] cluster_n0/cluster_n2 (the pure end-members)
    are built directly from the user's own input cell (whatever size),
    while cluster_n1 (AB) generally comes from icet's own primitive-cell-
    based enumeration -- a DIFFERENT cell size/atom count (see gqca.py::
    build_ab_structure). The raw FreeEng values are therefore NOT on a
    consistent per-site basis across the 3 folders: dividing by each
    folder's OWN count of sublattice+species-B atoms (re-derived from
    that folder's own structure.fdf, not assumed/copied from another
    folder) is what makes e_0/e_1/e_2 comparable at all before they're
    fed to gqca_solver.mixing_energies.

    Returns (e_per_site, n_sites, raw_free_energy) or (None, None, None)
    if the energy or structure file is missing/unreadable.
    """
    out_path = os.path.join(cluster_dir, out_file)
    free_eng = siesta_log.get_free_energy(out_path)
    if free_eng is None:
        return None, None, None
    struct_path = os.path.join(cluster_dir, "structure.fdf")
    try:
        structure = structure_io.read_fdf(struct_path)
    except (FileNotFoundError, ValueError):
        return None, None, None
    counts = structure_io.atom_counts(structure)
    n_sites = counts.get(sublattice, 0) + counts.get(species_b, 0)
    if n_sites == 0:
        return None, None, None
    return free_eng / n_sites, n_sites, free_eng


def write_results_csv(csv_path, rows):
    """Long-format CSV, one row per (x,T) grid point -- x, T, eta, x_0,
    x_1, x_2, H_mix, S_mix, G_mix, SRO.
    """
    with open(csv_path, "w") as f:
        f.write("x,T_K,eta,x_0,x_1,x_2,H_mix_eV,S_mix_eV_per_K,G_mix_eV,SRO\n")
        for row in rows:
            x, T, eta, x_j, H, S, G, sro = row
            eta_str = "inf" if eta == float("inf") else f"{eta:.8g}"
            f.write(f"{x:.6f},{T:.4f},{eta_str},{x_j[0]:.8f},{x_j[1]:.8f},{x_j[2]:.8f},"
                    f"{H:.8f},{S:.8e},{G:.8f},{sro:.8f}\n")


def write_spectrum_plot(dat_path, gplot_path, rows_at_t, t_ref):
    """Delta-H_mix/Delta-S_mix/Delta-G_mix vs x at a single reference
    temperature `t_ref` -- same .dat/.gplot gnuplot-driven convention as
    the rest of this suite's spectrum-style outputs (e.g. raman_analysis.
    py::write_spectrum_plot), not an embedded matplotlib figure.
    """
    with open(dat_path, "w") as f:
        f.write(f"# stb-gqcaAnalysis -- mixing thermodynamics vs composition at T={t_ref:.2f} K\n")
        f.write("# columns: 1=x 2=H_mix(eV) 3=S_mix(eV/K) 4=G_mix(eV) 5=SRO\n")
        for x, T, eta, x_j, H, S, G, sro in rows_at_t:
            f.write(f"{x:12.6f} {H:14.6f} {S:16.8e} {G:14.6f} {sro:12.6f}\n")

    pdf_name = os.path.splitext(os.path.basename(dat_path))[0] + ".pdf"
    with open(gplot_path, "w") as f:
        f.writelines([
            '# --- STB Plot Configuration ---\n',
            '# Generated by stb-gqcaAnalysis\n',
            'set terminal pdfcairo enhanced color font "Arial,14" size 8,5\n',
            f'set output "{pdf_name}"\n\n',
            f'set title "GQCA mixing thermodynamics at T = {t_ref:.0f} K"\n',
            'set xlabel "x (alloy composition, B fraction)"\n',
            'set ylabel "Energy (eV per sublattice site)"\n',
            'set grid\n',
            'set key top right\n',
            f'plot "{os.path.basename(dat_path)}" using 1:2 with lines lw 2 lc rgb "#2255cc" '
            'title "{/Symbol D}H_{mix}", \\\n'
            f'     "{os.path.basename(dat_path)}" using 1:4 with lines lw 2 lc rgb "#cc2222" '
            'title "{/Symbol D}G_{mix}"\n',
        ])


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Stage 2 of 2: solves the GQCA pair-cluster mass-action "
        "equilibrium and reports alloy mixing thermodynamics.", 'bold')}
Reads the 3 relaxed cluster_n{{0,1,2}}/ folders written by stb-gqca (AA/AB/BB pair-cluster
representatives), normalizes each to a per-sublattice-site energy (the 3 folders generally have
different cell sizes -- see read_cluster_energy's own docstring), and solves the mass-action
equilibrium (Sher, van Schilfgaarde, Chen, Chen, Phys. Rev. B 36, 4279, 1987) over a composition x /
temperature T grid to report Delta-H_mix, Delta-S_mix, Delta-G_mix, and a short-range-order (SRO)
parameter.

[UNVERIFIED] The mixing-entropy formula (-kB * Sum_j x_j*ln(x_j/g_j), the standard "gas of
independent clusters" form) and the SRO-parameter normalization could not be pinned to a fully
confirmed closed form from the primary literature during this tool's development -- see
core/gqca_solver.py's own docstrings. The mixing-enthalpy and mass-action-law equations ARE
confirmed against a secondary source that restates the original 1987 formalism. Every report/CSV
column derived from an [UNVERIFIED] quantity is documented as such below.

Also runs a discrete-grid common-tangent (miscibility gap) search on Delta-G_mix(x) at every
sampled T -- kept for completeness, but [STRUCTURAL LIMITATION] Delta-G_mix(x) from this
independent-cluster formalism is PROVABLY CONVEX for any Delta_j (a convex-analysis theorem, not
a numerical coincidence), so a gap is essentially never found here regardless of the input
energies -- this model captures short-range order, not genuine phase separation. See
core/gqca_solver.py::common_tangent_gaps' own docstring for the full argument.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s --directory gqca_study --temp 800\n"
               "  %(prog)s --directory gqca_study --temp-min 300 --temp-max 1500 --temp-points 5\n"
    )

    parser.add_argument("--directory", type=str, default="gqca_study",
                         help="Root directory written by stb-gqca (default: gqca_study).")
    parser.add_argument("--file", type=str, default="calc.out",
                         help="SIESTA output filename inside each cluster_n* folder (default: calc.out).")
    parser.add_argument("--x-min", type=float, default=0.0, help="Composition grid minimum (default: 0.0).")
    parser.add_argument("--x-max", type=float, default=1.0, help="Composition grid maximum (default: 1.0).")
    parser.add_argument("--x-points", type=int, default=21,
                         help="Number of composition grid points (default: 21).")
    parser.add_argument("--temp", type=float, default=None,
                         help="Single fixed temperature in K. Mutually exclusive with --temp-min/"
                              "--temp-max/--temp-points.")
    parser.add_argument("--temp-min", type=float, default=None, help="Temperature sweep minimum (K).")
    parser.add_argument("--temp-max", type=float, default=None, help="Temperature sweep maximum (K).")
    parser.add_argument("--temp-points", type=int, default=5,
                         help="Number of temperature sweep points (default: 5).")
    parser.add_argument("--force-tolerance", type=float, default=0.05,
                         help="Residual atomic force in eV/Ang (default: 0.05) above which a cluster "
                              "folder's calc.out is flagged as possibly not relaxed. Advisory only, "
                              "never blocks the result.")
    parser.add_argument("-o", "--output", type=str, default="gqca_results",
                         help="Base filename (no extension) for the .csv/.dat/.gplot outputs "
                              "(default: gqca_results).")
    parser.add_argument("-v", "--version", action="version", version=f"stb-gqcaAnalysis {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("GQCA WORKFLOW -- STAGE 2: MASS-ACTION ANALYSIS", 'bold'))
    print("-" * 60)

    if args.temp is not None and (args.temp_min is not None or args.temp_max is not None):
        print(color_text("[ERROR] --temp is mutually exclusive with --temp-min/--temp-max.", 'red'))
        sys.exit(1)
    if args.temp is None and (args.temp_min is None or args.temp_max is None):
        print(color_text("[ERROR] Provide either --temp, or both --temp-min and --temp-max.", 'red'))
        sys.exit(1)
    if args.x_points < 2:
        print(color_text("[ERROR] --x-points must be >= 2.", 'red'))
        sys.exit(1)

    if not os.path.isdir(args.directory):
        print(color_text(f"[ERROR] Directory '{args.directory}' not found -- run stb-gqca first.", 'red'))
        sys.exit(1)

    sublattice, species_b, ab_deviation = read_stage1_metadata(args.directory)
    if sublattice is None:
        print(color_text(
            f"[ERROR] Could not read Sublattice/Species B from "
            f"'{os.path.join(args.directory, 'gqca_stage1.txt')}' -- run stb-gqca first.", 'red'))
        sys.exit(1)

    N = 2
    g = gqca_solver.degeneracies(N)

    report_path = os.path.join(args.directory, REPORT_FILE)
    with open(report_path, 'w') as f_out:
        print_dual(f"{color_text('===== GQCA STAGE 2 REPORT (MASS-ACTION ANALYSIS) =====', 'magenta')}", f_out)

        print_dual(f"\n{color_text('[0] RUN METADATA', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        print_dual(f"Directory       : {args.directory}", f_out)
        print_dual(f"Sublattice      : {sublattice}", f_out)
        print_dual(f"Species B       : {species_b}", f_out)
        print_dual(f"Cluster type    : pair (N=2), degeneracies g_j = {g}", f_out)

        print_dual(f"\n{color_text('[1] CLUSTER ENERGIES', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)

        e_per_site = [None, None, None]
        missing = []
        for n_j in (0, 1, 2):
            cluster_dir = os.path.join(args.directory, f"cluster_n{n_j}")
            out_path = os.path.join(cluster_dir, args.file)
            e_norm, n_sites, raw = read_cluster_energy(cluster_dir, args.file, sublattice, species_b)
            if e_norm is None:
                missing.append(cluster_dir)
                print_dual(color_text(
                    f"  [ERROR] Could not read '{args.file}' or 'structure.fdf' from '{cluster_dir}'.",
                    'red'), f_out)
                continue
            e_per_site[n_j] = e_norm
            print_dual(f"  cluster_n{n_j}: FreeEng = {raw:12.6f} eV over {n_sites} sublattice site(s) "
                        f"-> {e_norm:12.6f} eV/site", f_out)
            siesta_log.report_quality_diagnostics(f"cluster_n{n_j}", out_path, args.force_tolerance, f_out)

        if missing:
            print_dual(color_text(
                "\n[ERROR] Cannot proceed -- one or more cluster_n* results are missing or "
                "incomplete. Run SIESTA in every folder listed above, then re-run stb-gqcaAnalysis.",
                'red'), f_out)
            sys.exit(1)

        delta_e = gqca_solver.mixing_energies(e_per_site, N)
        print_dual(f"\nMixing energies Delta_j (excess relative to linear e_0/e_2 interpolation, "
                    f"local cluster composition):", f_out)
        for n_j in (0, 1, 2):
            print_dual(f"  Delta_{n_j} = {delta_e[n_j]:12.6f} eV/site", f_out)

        if ab_deviation >= 1e-6:
            print_dual(color_text(
                f"\n[NOTE] cluster_n1's pair correlation deviated by {ab_deviation:.4f} from the "
                "ideal -1.0 (every-pair-unlike) AB ordering at Stage 1 (see gqca_stage1.txt) -- "
                "the sublattice is geometrically frustrated (e.g. FCC/HCP/triangular), so this "
                "reference structure's local environment is not purely n_j=1. Delta_1 (and every "
                "downstream x_1/H_mix/S_mix/G_mix/SRO value derived from it) is therefore a less "
                "clean estimate than for a bipartite sublattice -- qualitatively meaningful, but "
                "not a rigorous quantitative bound.", 'yellow'), f_out)

        print_dual(f"\n{color_text('[2] COMPOSITION / TEMPERATURE GRID', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)

        x_points = args.x_points
        x_grid = [args.x_min + i * (args.x_max - args.x_min) / (x_points - 1) for i in range(x_points)]
        if args.temp is not None:
            t_grid = [args.temp]
            print_dual(f"x grid          : {args.x_min} .. {args.x_max}, {x_points} point(s)", f_out)
            print_dual(f"T               : {args.temp} K (fixed)", f_out)
        else:
            t_points = args.temp_points
            if t_points < 1:
                print_dual(color_text("[ERROR] --temp-points must be >= 1.", 'red'), f_out)
                sys.exit(1)
            if t_points == 1:
                t_grid = [args.temp_min]
            else:
                t_grid = [args.temp_min + i * (args.temp_max - args.temp_min) / (t_points - 1)
                          for i in range(t_points)]
            print_dual(f"x grid          : {args.x_min} .. {args.x_max}, {x_points} point(s)", f_out)
            print_dual(f"T grid          : {args.temp_min} .. {args.temp_max} K, {t_points} point(s)", f_out)

        rows = []
        for T in t_grid:
            for x in x_grid:
                eta = gqca_solver.solve_eta(x, T, g, delta_e, N)
                x_j = gqca_solver.cluster_fractions(eta, T, g, delta_e, N)
                H = gqca_solver.mixing_enthalpy(x_j, delta_e)
                S = gqca_solver.mixing_entropy(x_j, g)
                G = gqca_solver.free_energy(H, S, T)
                sro = gqca_solver.sro_parameter(x_j, x, N)
                rows.append((x, T, eta, x_j, H, S, G, sro))

        print_dual(f"Computed {len(rows)} grid point(s).", f_out)

        print_dual(f"\n{color_text('[3] SANITY DIAGNOSTICS', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        # Pure end-member limits reconfirmed against the REAL data just read (not a
        # standalone unit test) -- always runs, data-independent gate.
        t_check = t_grid[0]
        x_j_x0 = gqca_solver.cluster_fractions(gqca_solver.solve_eta(0.0, t_check, g, delta_e, N),
                                                 t_check, g, delta_e, N)
        x_j_x1 = gqca_solver.cluster_fractions(gqca_solver.solve_eta(1.0, t_check, g, delta_e, N),
                                                 t_check, g, delta_e, N)
        ok_x0 = abs(x_j_x0[0] - 1.0) < 1e-9 and abs(x_j_x0[1]) < 1e-9 and abs(x_j_x0[2]) < 1e-9
        ok_x1 = abs(x_j_x1[2] - 1.0) < 1e-9 and abs(x_j_x1[0]) < 1e-9 and abs(x_j_x1[1]) < 1e-9
        status_x0 = color_text("OK", 'green') if ok_x0 else color_text("FAIL", 'red')
        status_x1 = color_text("OK", 'green') if ok_x1 else color_text("FAIL", 'red')
        print_dual(f"  x=0 pure-A limit (x_j = [1,0,0]) : {status_x0}", f_out)
        print_dual(f"  x=1 pure-B limit (x_j = [0,0,1]) : {status_x1}", f_out)

        print_dual(f"\n{color_text('[4] MISCIBILITY GAP (COMMON-TANGENT)', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)
        print_dual(color_text(
            "[STRUCTURAL LIMITATION] Within this independent-cluster (Sher-GQCA) formalism, "
            "Delta_G_mix(x) is PROVABLY CONVEX at every composition for ANY Delta_j (a convex-"
            "analysis theorem, not a numerical coincidence -- see core/gqca_solver.py::"
            "common_tangent_gaps' own docstring), so a gap is essentially never found here "
            "regardless of how unfavorable the cluster energies are. This model captures short-"
            "range order but NOT genuine phase separation (that needs cluster-cluster "
            "correlations beyond this independent-cluster treatment). A 'no gap' result below is "
            "therefore an intrinsic property of the model, NOT a physical claim that this alloy "
            "is fully miscible -- do not report it as such.", 'yellow'), f_out)
        print_dual("Discrete-grid common-tangent construction on Delta_G_mix(x) at each sampled T "
                    "(kept for completeness/future formalism extensions) -- resolution limited by "
                    "--x-points (x_alpha/x_beta can only land on a sampled grid point).", f_out)

        phase_rows = []
        for T in t_grid:
            rows_at_t = [r for r in rows if r[1] == T]
            xs_t = np.array([r[0] for r in rows_at_t])
            gs_t = np.array([r[6] for r in rows_at_t])
            gaps = gqca_solver.common_tangent_gaps(xs_t, gs_t)
            if gaps:
                phase_rows.extend((T, xa, xb) for xa, xb in gaps)
                gaps_str = ", ".join(f"[{xa:.4f}, {xb:.4f}]" for xa, xb in gaps)
                print_dual(f"  T = {T:8.2f} K : miscibility gap(s) {gaps_str}", f_out)
            else:
                print_dual(f"  T = {T:8.2f} K : fully miscible (no gap) over the sampled x range", f_out)

        phase_csv_path = None
        if phase_rows:
            phase_csv_path = os.path.join(args.directory, f"{args.output}_phase_diagram.csv")
            with open(phase_csv_path, "w") as f:
                f.write("T_K,x_alpha,x_beta\n")
                for T, xa, xb in phase_rows:
                    f.write(f"{T:.4f},{xa:.6f},{xb:.6f}\n")

        print_dual(f"\n{color_text('[5] SUMMARY & FILES', 'magenta')}", f_out)
        print_dual("-" * 60, f_out)

        csv_path = os.path.join(args.directory, f"{args.output}.csv")
        write_results_csv(csv_path, rows)
        print_dual(f"CSV (full grid) : {csv_path}", f_out)
        if phase_csv_path:
            print_dual(f"Phase diagram (T vs. x_alpha/x_beta) : {phase_csv_path}", f_out)

        t_ref = t_grid[0]
        rows_at_t = [r for r in rows if r[1] == t_ref]
        dat_path = os.path.join(args.directory, f"{args.output}.dat")
        gplot_path = os.path.join(args.directory, f"{args.output}.gplot")
        write_spectrum_plot(dat_path, gplot_path, rows_at_t, t_ref)
        print_dual(f"Plot (T={t_ref:.0f} K) : {dat_path}, {gplot_path} "
                    f"(cd {args.directory} && gnuplot {os.path.basename(gplot_path)})", f_out)
        print_dual(f"Report          : {report_path}", f_out)

        print_dual(color_text(
            "\n[UNVERIFIED] S_mix and SRO use provisional formulas -- see core/gqca_solver.py's "
            "mixing_entropy/sro_parameter docstrings.", 'yellow'), f_out)

    print("\n[INFO] Complete job!")
    print("\n" + "-" * 60)
    print(color_text("GQCA mass-action analysis complete!\n", 'bold'))


if __name__ == "__main__":
    main()
