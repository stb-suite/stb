#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

from stb import __version__ as VERSION

import os
import re
import sys
import argparse
import datetime
import numpy as np
from stb.core import siesta_log, structure_io, gqca_solver
from stb.core.cli import color_text, show_intro, print_dual, print_section, print_table

REPORT_FILE = "gqca_stage2.txt"

KB_EV_PER_K = 8.617333262e-5
EV_TO_KJ_MOL = 96.485332
EV_TO_KCAL_MOL = 23.060549
DEFAULT_T_MIN, DEFAULT_T_MAX, DEFAULT_T_POINTS = 300.0, 1500.0, 5

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


# ---------------------------------------------------------------------------
# Level of theory / convergence / geometry readers (all degrade to None)
# ---------------------------------------------------------------------------

_THEORY_KEYS = [
    ("xcfunctional", "XC functional"), ("xcauthors", "XC authors"),
    ("paobasistype", "Basis type"), ("paobasissize", "Basis size"),
    ("paoenergyshift", "PAO energy shift"), ("meshcutoff", "Mesh cutoff"),
    ("spin", "Spin"), ("electronictemperature", "Electronic temperature"),
    ("dmtolerance", "DM tolerance"), ("mdtypeofrun", "MD run type"),
    ("mdvariablecell", "Variable cell"), ("mdmaxforcetol", "Max force tol"),
    ("mdmaxstresstol", "Max stress tol"),
]
_KGRID_LIST_RE = re.compile(r'kgrid[._\-]?monkhorstpack\s*\[\s*(\d+)\s+(\d+)\s+(\d+)\s*\]', re.IGNORECASE)
_KGRID_BLOCK_RE = re.compile(r'%block\s+kgrid[._\-]?monkhorstpack\s*\n(.*?)%endblock', re.IGNORECASE | re.DOTALL)


def _fdf_key(label):
    """fdf labels ignore case and '.', '_', '-'."""
    return re.sub(r'[._\-]', '', label.lower())


def read_level_of_theory(cluster_dir):
    """{normalized key: value} + 'kgrid' for one cluster folder, merging
    config_extra.fdf then calc.fdf with SIESTA's first-occurrence-wins rule
    (block contents and comments skipped). Empty dict if neither file exists."""
    texts = []
    for name in ("config_extra.fdf", "calc.fdf"):
        path = os.path.join(cluster_dir, name)
        if os.path.isfile(path):
            with open(path, 'r', encoding='utf-8', errors='ignore') as f:
                texts.append(f.read())
    values, in_block = {}, False
    for text in texts:
        for raw in text.splitlines():
            line = raw.split('#', 1)[0].strip()
            if not line:
                continue
            low = line.lower()
            if low.startswith('%block'):
                in_block = True
                continue
            if low.startswith('%endblock'):
                in_block = False
                continue
            if in_block or low.startswith('%'):
                continue
            parts = line.split(None, 1)
            if len(parts) == 2:
                values.setdefault(_fdf_key(parts[0]), parts[1].strip())
    combined = "\n".join(texts)
    m = _KGRID_LIST_RE.search(combined)
    if m:
        values['kgrid'] = "x".join(m.groups())
    else:
        mb = _KGRID_BLOCK_RE.search(combined)
        if mb:
            try:
                rows = [r.split() for r in mb.group(1).strip().splitlines()][:3]
                values['kgrid'] = "x".join(str(int(float(r[i]))) for i, r in enumerate(rows))
            except (ValueError, IndexError):
                pass
    nsteps = values.get('mdnumcgsteps') or values.get('mdsteps')
    if nsteps:
        values['mdnumcgsteps'] = nsteps
    return values


def _lattice_metrics(lattice):
    """(a, b, c, volume) of a 3x3 lattice matrix in Ang, or None."""
    if lattice is None:
        return None
    lat = np.asarray(lattice, dtype=float)
    if lat.shape != (3, 3):
        return None
    return (*[float(np.linalg.norm(v)) for v in lat], abs(float(np.linalg.det(lat))))


def read_cluster_details(cluster_dir, out_file, n_sites):
    """Convergence + geometry facts of one finished cluster folder. Every
    field is None when the log doesn't carry it (a fabricated/minimal calc.out
    must not crash the analysis)."""
    out_path = os.path.join(cluster_dir, out_file)
    scf_ok, scf_iter = siesta_log.get_scf_convergence(out_path)
    try:
        with open(out_path, 'r', encoding='utf-8', errors='ignore') as f:
            tail = f.read()
        finished = "Job completed" in tail
    except OSError:
        finished = False
    stress = siesta_log.get_stress_voigt_kbar(out_path)
    try:
        input_lat = structure_io.read_fdf(os.path.join(cluster_dir, "structure.fdf")).lattice
    except (FileNotFoundError, ValueError):
        input_lat = None
    relaxed_lat = siesta_log.get_outcell(out_path)
    metrics_in = _lattice_metrics(input_lat)
    metrics_out = _lattice_metrics(relaxed_lat) if relaxed_lat is not None else None
    return {
        'scf_ok': scf_ok, 'scf_iter': scf_iter, 'finished': finished,
        'moves': siesta_log.count_geometry_steps(out_path),
        'max_force': siesta_log.get_max_force(out_path),
        'stress': stress[:3] if stress else None,
        'input': metrics_in, 'relaxed': metrics_out,
        'v_site_input': metrics_in[3] / n_sites if metrics_in else None,
        'v_site': metrics_out[3] / n_sites if metrics_out else None,
    }


def load_relaxed_atoms(cluster_dir, out_file):
    """(ase.Atoms, source) with the relaxed geometry (positions AND cell) of a
    cluster folder: from the SIESTA *.XV (sisl) when present, else the input
    structure.fdf with the last outcell of the log. (None, None) on failure."""
    import glob
    try:
        from ase import Atoms
    except ImportError:
        return None, None
    xv_files = sorted(glob.glob(os.path.join(cluster_dir, "*.XV")))
    if xv_files:
        try:
            from stb.core.deps import require_sisl
            geom = require_sisl().get_sile(xv_files[0]).read_geometry()
            atoms = Atoms(symbols=[a.symbol for a in geom.atoms], positions=np.asarray(geom.xyz),
                          cell=np.asarray(geom.cell), pbc=True)
            return atoms, f"relaxed ({os.path.basename(xv_files[0])})"
        except Exception:
            pass
    try:
        from pymatgen.io.ase import AseAtomsAdaptor
        fdf = structure_io.read_fdf(os.path.join(cluster_dir, "structure.fdf"))
        atoms = AseAtomsAdaptor.get_atoms(structure_io.to_pymatgen(fdf))
        outcell = siesta_log.get_outcell(os.path.join(cluster_dir, out_file))
        if outcell is not None:
            atoms.set_cell(np.asarray(outcell), scale_atoms=True)
            return atoms, "input positions + relaxed cell (no .XV found)"
        return atoms, "input structure (no .XV / outcell found)"
    except Exception:
        return None, None


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


def write_spectrum_plot(dat_path, gplot_path, rows, t_grid):
    """Delta-H_mix / Delta-G_mix vs x for EVERY sampled temperature (one
    gnuplot 'index' block per T) -- same .dat/.gplot gnuplot-driven convention
    as the rest of this suite's workflow outputs."""
    with open(dat_path, "w") as f:
        f.write("# stb-gqcaAnalysis -- mixing thermodynamics vs composition, one block per T\n")
        f.write("# columns: 1=x 2=H_mix(eV) 3=S_mix(eV/K) 4=G_mix(eV) 5=SRO\n")
        for T in t_grid:
            f.write(f"# T = {T:.2f} K\n")
            for x, T_r, eta, x_j, H, S, G, sro in rows:
                if T_r == T:
                    f.write(f"{x:12.6f} {H:14.6f} {S:16.8e} {G:14.6f} {sro:12.6f}\n")
            f.write("\n\n")

    pdf_name = os.path.splitext(os.path.basename(dat_path))[0] + ".pdf"
    dat = os.path.basename(dat_path)
    clauses = []
    for i, T in enumerate(t_grid):
        clauses.append(f'"{dat}" index {i} using 1:2 with lines lw 2 dt 1 lc {i + 1} '
                       f'title "{{/Symbol D}}H {T:.0f} K"')
        clauses.append(f'"{dat}" index {i} using 1:4 with lines lw 2 dt 2 lc {i + 1} '
                       f'title "{{/Symbol D}}G {T:.0f} K"')
    with open(gplot_path, "w") as f:
        f.writelines([
            '# --- STB Plot Configuration ---\n',
            '# Generated by stb-gqcaAnalysis (solid = DeltaH_mix, dashed = DeltaG_mix)\n',
            'set terminal pdfcairo enhanced color font "Arial,14" size 8,5\n',
            f'set output "{pdf_name}"\n\n',
            'set title "GQCA mixing thermodynamics"\n',
            'set xlabel "x (alloy composition, B fraction)"\n',
            'set ylabel "Energy (eV per sublattice site)"\n',
            'set grid\n',
            'set key outside right\n',
            'plot ' + ', \\\n     '.join(clauses) + '\n',
        ])


def make_dashboard(rows, t_grid, t_ref, delta_e, e_per_site, details, sublattice, species_b,
                    png_path, show):
    """One 2x3 matplotlib figure that completes the analysis visually:
    DeltaH(x), DeltaG(x) (every T), cluster fractions vs the random
    reference (at t_ref), SRO(x), cluster energy diagram, and the relaxed
    volume per site vs composition (Vegard check). Saved to png_path when
    given and/or shown on screen. Returns True when something was produced."""
    try:
        import matplotlib
        if not show:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print(color_text("[WARNING] matplotlib is not installed -- skipping the dashboard.", 'yellow'))
        return False
    if show and (matplotlib.get_backend().lower() == "agg"
                 or (sys.platform.startswith("linux") and not (os.environ.get("DISPLAY")
                                                               or os.environ.get("WAYLAND_DISPLAY")))):
        print(color_text("[WARNING] --show needs a display (local desktop or `ssh -X`); "
                          "only saving the figure.", 'yellow'))
        show = False
        plt.switch_backend("Agg")

    okabe = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000"]
    colors = {T: okabe[i % len(okabe)] for i, T in enumerate(t_grid)}
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.5))
    ax_h, ax_g, ax_f, ax_s, ax_e, ax_v = axes.ravel()

    for T in t_grid:
        pts = [r for r in rows if r[1] == T]
        xs = [r[0] for r in pts]
        ax_h.plot(xs, [r[4] * 1000 for r in pts], color=colors[T], lw=2, label=f"{T:.0f} K")
        ax_g.plot(xs, [r[6] * 1000 for r in pts], color=colors[T], lw=2, label=f"{T:.0f} K")
        ax_s.plot(xs, [r[7] for r in pts], color=colors[T], lw=2, label=f"{T:.0f} K")
    ax_h.set_title(r"$\Delta H_{mix}$"); ax_h.set_ylabel("meV / site")
    ax_g.set_title(r"$\Delta G_{mix} = \Delta H - T\Delta S$"); ax_g.set_ylabel("meV / site")
    ax_s.set_title("Short-range order (SRO)"); ax_s.set_ylabel("SRO  (>0 clustering, <0 ordering)")
    ax_s.axhline(0, color="0.5", lw=1)
    for ax in (ax_h, ax_g, ax_s):
        ax.set_xlabel(f"x  ({sublattice}$_{{1-x}}${species_b}$_x$)")
        ax.grid(alpha=0.3); ax.legend(title="T", fontsize=8)

    pts = [r for r in rows if r[1] == t_ref]
    xs = np.array([r[0] for r in pts])
    labels = [f"{sublattice}-{sublattice} (AA)", f"{sublattice}-{species_b} (AB)", f"{species_b}-{species_b} (BB)"]
    for j, col in zip(range(3), okabe[:3]):
        ax_f.plot(xs, [r[3][j] for r in pts], color=col, lw=2, label=labels[j])
        ax_f.plot(xs, [gqca_solver.random_reference_fractions(x, 2)[j] for x in xs], color=col,
                  lw=1.2, ls="--")
    ax_f.set_title(f"Pair fractions x$_j$ at {t_ref:.0f} K  (dashed = random alloy)")
    ax_f.set_xlabel("x"); ax_f.set_ylabel("pair fraction"); ax_f.grid(alpha=0.3); ax_f.legend(fontsize=8)

    e_rel = [(e - (e_per_site[0] * (1 - n / 2) + e_per_site[2] * (n / 2))) * 1000 for n, e in enumerate(e_per_site)]
    bars = ax_e.bar(["AA\n(n$_j$=0)", "AB\n(n$_j$=1)", "BB\n(n$_j$=2)"], [d * 1000 for d in delta_e],
                    color=okabe[:3])
    for b, d in zip(bars, delta_e):
        ax_e.annotate(f"{d * 1000:+.1f}", (b.get_x() + b.get_width() / 2, d * 1000),
                      ha="center", va="bottom" if d >= 0 else "top", fontsize=10)
    ax_e.axhline(0, color="0.4", lw=1)
    ax_e.set_title(r"Pair excess energies $\Delta_j$"); ax_e.set_ylabel("meV / site"); ax_e.grid(alpha=0.3, axis="y")

    v = [d['v_site'] if d and d.get('v_site') else None for d in details]
    if all(val is not None for val in v):
        ax_v.plot([0, 0.5, 1], v, "o", color=okabe[1], ms=9, label="relaxed (DFT)")
        ax_v.plot([0, 1], [v[0], v[2]], "--", color="0.4", label="Vegard (linear)")
        exc = (v[1] - 0.5 * (v[0] + v[2])) / (0.5 * (v[0] + v[2])) * 100
        ax_v.set_title(f"Relaxed volume per site  (AB excess {exc:+.2f} %)")
        ax_v.set_xlabel("x (cluster composition)"); ax_v.set_ylabel(r"V / site ($\AA^3$)")
        ax_v.grid(alpha=0.3); ax_v.legend(fontsize=8)
    else:
        ax_v.text(0.5, 0.5, "relaxed cell not found\nin the calc.out files", ha="center", va="center")
        ax_v.set_axis_off()

    fig.suptitle(f"GQCA pair-cluster analysis -- {sublattice}$_{{1-x}}${species_b}$_x$", fontsize=14)
    fig.tight_layout()
    if png_path:
        fig.savefig(png_path, dpi=150)
    if show:
        plt.show()
    plt.close(fig)
    return True


def _fmt_mev(ev):
    return f"{ev * 1000:+.2f}"


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Stage 2 of 2: solves the GQCA pair-cluster mass-action "
        "equilibrium and reports alloy mixing thermodynamics.", 'bold')}
Reads the 3 relaxed cluster_n{{0,1,2}}/ folders written by stb-gqca (AA/AB/BB pair-cluster
representatives; SIESTA output read from <folder>/calc.out by default, level of theory read back
from each folder's own calc.fdf + config_extra.fdf), normalizes each to a per-sublattice-site energy
(the 3 folders generally have different cell sizes -- see read_cluster_energy's own docstring), and
solves the mass-action equilibrium (Sher, van Schilfgaarde, Chen, Chen, Phys. Rev. B 36, 4279, 1987)
over a composition x / temperature T grid to report Delta-H_mix, Delta-S_mix, Delta-G_mix, and a
short-range-order (SRO) parameter.

The report covers: level-of-theory consistency across the 3 folders, per-folder convergence and
relaxed geometry (cell, volume per site, residual force/stress), pair excess energies in several
units, an equivalent regular-solution estimate (Omega, T_c), tables at chosen target compositions,
sanity checks, and a plain-language interpretation. Optional on-screen/PNG dashboard (--show/--png),
interactive 3D viewer of the 3 relaxed structures (--view) and CIF/XSF/... export (--export-relaxed).

[UNVERIFIED] The closed form of the mixing entropy (-(kB/N) * Sum_j x_j*ln(x_j/g_j) per site, the
standard "gas of independent clusters" form -- exact ideal-solution entropy in the random limit)
and the SRO-parameter normalization could not be pinned to a fully confirmed closed form from the
primary literature during this tool's development -- see core/gqca_solver.py's own docstrings
(which also document the per-site vs per-cluster normalization). The mixing-enthalpy and mass-action-law equations ARE
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
               "  %(prog)s                                   # defaults: ./gqca_study, T = 300..1500 K (5 points)\n"
               "  %(prog)s --directory gqca_study --temp 800\n"
               "  %(prog)s --temp-min 300 --temp-max 1500 --temp-points 5 --x-target 0.25 0.5 0.75\n"
               "  %(prog)s --show --view                     # dashboard + 3D viewer on screen\n"
               "  %(prog)s --png --export-relaxed cif        # figure + relaxed structures to files\n"
    )

    parser.add_argument("--directory", type=str, default="gqca_study",
                         help="Root directory written by stb-gqca (default: gqca_study).")
    parser.add_argument("--file", type=str, default="calc.out",
                         help="SIESTA output filename inside each cluster_n* folder (default: calc.out).")
    parser.add_argument("--x-min", type=float, default=0.0, help="Composition grid minimum (default: 0.0).")
    parser.add_argument("--x-max", type=float, default=1.0, help="Composition grid maximum (default: 1.0).")
    parser.add_argument("--x-points", type=int, default=21,
                         help="Number of composition grid points (default: 21).")
    parser.add_argument("--x-target", type=float, nargs="+", default=[0.5],
                         help="Composition(s) x for the detailed per-temperature tables (default: 0.5).")
    parser.add_argument("--temp", type=float, default=None,
                         help="Single fixed temperature in K. Mutually exclusive with --temp-min/"
                              "--temp-max/--temp-points.")
    parser.add_argument("--temp-min", type=float, default=None,
                         help=f"Temperature sweep minimum (K). Without any temperature option the sweep "
                              f"{DEFAULT_T_MIN:.0f}..{DEFAULT_T_MAX:.0f} K is used.")
    parser.add_argument("--temp-max", type=float, default=None, help="Temperature sweep maximum (K).")
    parser.add_argument("--temp-points", type=int, default=DEFAULT_T_POINTS,
                         help=f"Number of temperature sweep points (default: {DEFAULT_T_POINTS}).")
    parser.add_argument("--force-tolerance", type=float, default=0.05,
                         help="Residual atomic force in eV/Ang (default: 0.05) above which a cluster "
                              "folder's calc.out is flagged as possibly not relaxed. Advisory only, "
                              "never blocks the result.")
    parser.add_argument("-o", "--output", type=str, default="gqca_results",
                         help="Base filename (no extension) for the .csv/.dat/.gplot/.png outputs "
                              "(default: gqca_results).")
    parser.add_argument("--show", action="store_true",
                         help="Open the analysis dashboard (matplotlib) on screen. Needs a display.")
    parser.add_argument("--png", action="store_true",
                         help="Save the dashboard as <directory>/<output>_summary.png.")
    parser.add_argument("--plot-temp", type=float, default=None,
                         help="Temperature (K) used by the single-T dashboard panels (default: the lowest "
                              "sampled T; the nearest sampled T is used).")
    parser.add_argument("--view", action="store_true",
                         help="Open each relaxed cluster structure in ASE's interactive 3D viewer, one "
                              "after the other (needs a display).")
    parser.add_argument("--export-relaxed", nargs="?", const="cif", default=None,
                         choices=["cif", "xsf", "xyz", "vasp"], metavar="FORMAT",
                         help="Write the 3 relaxed structures (positions + relaxed cell) to "
                              "<directory>/relaxed/ as cif (default), xsf, xyz or vasp.")
    parser.add_argument("--save-report", action="store_true",
                         help=f"Also persist the report to <directory>/{REPORT_FILE}. Off by default.")
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

    default_sweep = False
    if args.temp is not None and (args.temp_min is not None or args.temp_max is not None):
        print(color_text("[ERROR] --temp is mutually exclusive with --temp-min/--temp-max.", 'red'))
        sys.exit(1)
    if args.temp is None and (args.temp_min is None) != (args.temp_max is None):
        print(color_text("[ERROR] Provide both --temp-min and --temp-max (or --temp, or neither "
                          "for the default sweep).", 'red'))
        sys.exit(1)
    if args.temp is None and args.temp_min is None:
        args.temp_min, args.temp_max, default_sweep = DEFAULT_T_MIN, DEFAULT_T_MAX, True
    if args.x_points < 2:
        print(color_text("[ERROR] --x-points must be >= 2.", 'red'))
        sys.exit(1)
    if any(not (0.0 <= xt <= 1.0) for xt in args.x_target):
        print(color_text("[ERROR] --x-target values must lie in [0, 1].", 'red'))
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

    report_path = os.path.join(args.directory, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, 'w') if report_path else None

    def abort(msg):
        print_dual(color_text(msg, 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_dual(f"{color_text('===== GQCA STAGE 2 REPORT (MASS-ACTION ANALYSIS) =====', 'magenta')}", f_out)

    # ------------------------------------------------------------------ [0]
    print_section('[0] RUN METADATA', f_out)
    print_dual(f"Tool            : stb-gqcaAnalysis {VERSION}", f_out)
    print_dual(f"Date            : {datetime.datetime.now():%Y-%m-%d %H:%M}", f_out)
    print_dual(f"Directory       : {args.directory}", f_out)
    print_dual(f"SIESTA output   : <cluster>/{args.file}" + ("  (default)" if args.file == "calc.out" else ""), f_out)
    print_dual(f"Alloy           : {sublattice}(1-x){species_b}(x) on the {sublattice} sublattice", f_out)
    print_dual(f"Sublattice (A)  : {sublattice}", f_out)
    print_dual(f"Species B       : {species_b}", f_out)
    print_dual(f"Cluster type    : pair (N=2), degeneracies g_j = {g}", f_out)

    # ------------------------------------------------------------------ [1]
    names = {0: f"AA ({sublattice}-{sublattice})", 1: f"AB ({sublattice}-{species_b})",
             2: f"BB ({species_b}-{species_b})"}
    cluster_dirs = {n: os.path.join(args.directory, f"cluster_n{n}") for n in (0, 1, 2)}

    print_section('[1] LEVEL OF THEORY (read back from each folder)', f_out)
    theory = {n: read_level_of_theory(cluster_dirs[n]) for n in (0, 1, 2)}
    if not any(theory.values()):
        print_dual("  (no calc.fdf / config_extra.fdf found in the cluster folders)", f_out)
    else:
        theory_rows, mismatches = [], []
        for key, label in _THEORY_KEYS + [("mdnumcgsteps", "CG steps (max)"), ("kgrid", "k-grid")]:
            vals = [theory[n].get(key) for n in (0, 1, 2)]
            if all(v is None for v in vals):
                continue
            differ = len({v for v in vals if v is not None}) > 1 or any(v is None for v in vals)
            expected_diff = key == "kgrid"
            if differ and not expected_diff:
                mismatches.append(label)
            theory_rows.append(([label] + [v if v is not None else "-" for v in vals]
                                 + ["differs (k-grid follows each cell)" if (differ and expected_diff)
                                    else ("DIFFERS" if differ else "")],
                                 'yellow' if (differ and not expected_diff) else None))
        print_table(["Setting", "cluster_n0", "cluster_n1", "cluster_n2", "Consistency"], theory_rows, f_out)
        if mismatches:
            print_dual(color_text(
                f"[WARNING] The 3 folders do NOT share the same level of theory ({', '.join(mismatches)}). "
                "Mixing energies are only meaningful when every cluster is computed identically -- "
                "re-run the folders that differ.", 'yellow'), f_out)
        else:
            print_dual(color_text("All 3 folders share the same level of theory (k-grid aside, which "
                                   "is density-based per cell).", 'green'), f_out)

    # ------------------------------------------------------------------ [2]
    print_section('[2] CLUSTER ENERGIES, CONVERGENCE & RELAXED GEOMETRY', f_out)
    e_per_site = [None, None, None]
    n_sites_of, raw_of, details = {}, {}, [None, None, None]
    missing = []
    for n_j in (0, 1, 2):
        e_norm, n_sites, raw = read_cluster_energy(cluster_dirs[n_j], args.file, sublattice, species_b)
        if e_norm is None:
            missing.append(cluster_dirs[n_j])
            print_dual(color_text(
                f"  [ERROR] Could not read '{args.file}' or 'structure.fdf' from '{cluster_dirs[n_j]}'.",
                'red'), f_out)
            continue
        e_per_site[n_j] = e_norm
        n_sites_of[n_j], raw_of[n_j] = n_sites, raw
        details[n_j] = read_cluster_details(cluster_dirs[n_j], args.file, n_sites)

    if missing:
        abort("\n[ERROR] Cannot proceed -- one or more cluster_n* results are missing or "
              "incomplete. Run SIESTA in every folder listed above, then re-run stb-gqcaAnalysis.")

    print_table(["Cluster", "Pair", "Sites", "FreeEng (eV)", "E/site (eV)"],
                [([f"cluster_n{n}", names[n], str(n_sites_of[n]), f"{raw_of[n]:.6f}",
                   f"{e_per_site[n]:.6f}"], None) for n in (0, 1, 2)], f_out)

    def _yn(flag, ok="yes", bad="NO"):
        return ok if flag else bad

    print_dual("\nConvergence:", f_out)
    conv_rows = []
    for n in (0, 1, 2):
        d = details[n]
        force_txt = f"{d['max_force']:.4f}" if d['max_force'] is not None else "-"
        stress_txt = " ".join(f"{v:+.2f}" for v in d['stress']) if d['stress'] else "-"
        flag = (not d['scf_ok']) or (d['max_force'] is not None and d['max_force'] > args.force_tolerance)
        conv_rows.append(([f"cluster_n{n}", _yn(d['finished'], "yes", "no"),
                            _yn(d['scf_ok']) + (f" ({d['scf_iter']} it.)" if d['scf_iter'] else ""),
                            str(d['moves']) if d['moves'] else "-", force_txt, stress_txt],
                           'yellow' if flag else None))
    print_table(["Cluster", "Job completed", "SCF converged", "Geom. moves", "Max force (eV/Ang)",
                 "Stress xx yy zz (kbar)"], conv_rows, f_out)
    for n in (0, 1, 2):
        siesta_log.report_quality_diagnostics(f"cluster_n{n}", os.path.join(cluster_dirs[n], args.file),
                                               args.force_tolerance, f_out)

    print_dual("\nRelaxed geometry (cell lengths in Ang; volume per sublattice site in Ang^3):", f_out)
    geo_rows = []
    for n in (0, 1, 2):
        d = details[n]
        rel, inp = d['relaxed'], d['input']
        cell_txt = f"{rel[0]:.3f} {rel[1]:.3f} {rel[2]:.3f}" if rel else "-"
        v_txt = f"{d['v_site']:.3f}" if d['v_site'] else "-"
        dv = (f"{(d['v_site'] / d['v_site_input'] - 1) * 100:+.2f} %"
              if d['v_site'] and d['v_site_input'] else "-")
        geo_rows.append(([f"cluster_n{n}", f"{inp[0]:.3f} {inp[1]:.3f} {inp[2]:.3f}" if inp else "-",
                           cell_txt, v_txt, dv], None))
    print_table(["Cluster", "Input cell a b c", "Relaxed cell a b c", "V/site relaxed", "dV vs input"],
                geo_rows, f_out)

    v = [d['v_site'] for d in details]
    if all(val is not None for val in v):
        v_mid = 0.5 * (v[0] + v[2])
        excess = (v[1] / v_mid - 1) * 100
        mismatch = ((v[2] / v[0]) ** (1.0 / 3.0) - 1) * 100
        print_dual(f"\nEffective lattice mismatch  (V_BB/V_AA)^(1/3) - 1 : {mismatch:+.2f} %", f_out)
        print_dual(f"AB volume vs Vegard (linear) average          : {excess:+.2f} %  "
                    f"(V_AB = {v[1]:.3f} vs {v_mid:.3f} Ang^3/site)", f_out)
    else:
        print_dual(color_text("\n[NOTE] The relaxed cell ('outcell') was not found in every calc.out -- "
                               "volume/Vegard diagnostics skipped.", 'cyan'), f_out)

    # ------------------------------------------------------------------ [3]
    print_section('[3] MIXING ENERGIES', f_out)
    delta_e = gqca_solver.mixing_energies(e_per_site, N)
    print_dual("Delta_j = excess of each pair cluster over the linear e_0/e_2 interpolation, per "
                "sublattice site:", f_out)
    print_table(["j", "Pair", "Delta_j (eV)", "Delta_j (meV)", "kJ/mol", "kcal/mol"],
                [([str(n), names[n], f"{delta_e[n]:+.6f}", f"{delta_e[n] * 1000:+.2f}",
                   f"{delta_e[n] * EV_TO_KJ_MOL:+.2f}", f"{delta_e[n] * EV_TO_KCAL_MOL:+.2f}"], None)
                 for n in (0, 1, 2)], f_out)

    d1 = float(delta_e[1])
    omega = 2.0 * d1
    print_dual("\nRegular-solution equivalent (random-mixing limit of this model, H = Omega x (1-x)):", f_out)
    print_dual(f"  Omega = 2*Delta_1            : {omega * 1000:+.2f} meV/site  "
                f"({omega * EV_TO_KCAL_MOL:+.2f} kcal/mol, {omega * EV_TO_KJ_MOL:+.2f} kJ/mol)", f_out)
    print_dual(f"  H_mix(x=0.5, random limit)   : {omega / 4 * 1000:+.2f} meV/site", f_out)
    if omega > 0:
        print_dual(f"  T_c = Omega / (2 kB)         : {omega / (2 * KB_EV_PER_K):.0f} K  "
                    "(critical temperature of the equivalent REGULAR solution -- an estimate; "
                    "this model itself cannot show a gap, see [7])", f_out)
    else:
        print_dual("  T_c                          : none (Omega <= 0: no phase-separation tendency "
                    "in the regular-solution picture)", f_out)

    if ab_deviation >= 1e-6:
        print_dual(color_text(
            f"\n[NOTE] cluster_n1's pair correlation deviated by {ab_deviation:.4f} from the "
            "ideal -1.0 (every-pair-unlike) AB ordering at Stage 1 (see gqca_stage1.txt) -- "
            "the sublattice is geometrically frustrated (e.g. FCC/HCP/triangular), so this "
            "reference structure's local environment is not purely n_j=1. Delta_1 (and every "
            "downstream x_1/H_mix/S_mix/G_mix/SRO value derived from it) is therefore a less "
            "clean estimate than for a bipartite sublattice -- qualitatively meaningful, but "
            "not a rigorous quantitative bound.", 'yellow'), f_out)

    # ------------------------------------------------------------------ [4]
    print_section('[4] COMPOSITION / TEMPERATURE GRID', f_out)
    x_points = args.x_points
    x_grid = [args.x_min + i * (args.x_max - args.x_min) / (x_points - 1) for i in range(x_points)]
    if args.temp is not None:
        t_grid = [args.temp]
        print_dual(f"x grid          : {args.x_min} .. {args.x_max}, {x_points} point(s)", f_out)
        print_dual(f"T               : {args.temp} K (fixed)", f_out)
    else:
        t_points = args.temp_points
        if t_points < 1:
            abort("[ERROR] --temp-points must be >= 1.")
        t_grid = ([args.temp_min] if t_points == 1 else
                  [args.temp_min + i * (args.temp_max - args.temp_min) / (t_points - 1)
                   for i in range(t_points)])
        print_dual(f"x grid          : {args.x_min} .. {args.x_max}, {x_points} point(s)", f_out)
        print_dual(f"T grid          : {args.temp_min} .. {args.temp_max} K, {t_points} point(s)"
                    + ("  (default sweep; use --temp or --temp-min/--temp-max to change)" if default_sweep else ""),
                    f_out)

    def point(x, T):
        eta = gqca_solver.solve_eta(x, T, g, delta_e, N)
        x_j = gqca_solver.cluster_fractions(eta, T, g, delta_e, N)
        H = gqca_solver.mixing_enthalpy(x_j, delta_e)
        S = gqca_solver.mixing_entropy(x_j, g)
        return (x, T, eta, x_j, H, S, gqca_solver.free_energy(H, S, T), gqca_solver.sro_parameter(x_j, x, N))

    rows = [point(x, T) for T in t_grid for x in x_grid]
    print_dual(f"Computed {len(rows)} grid point(s).", f_out)

    print_dual("\nResults at the target composition(s) (energies in meV per sublattice site, "
                "S in meV/K):", f_out)
    for xt in args.x_target:
        print_dual(f"\n  x = {xt:.3f}   [random-alloy reference: x_0={(1 - xt) ** 2:.3f}, "
                    f"x_1={2 * xt * (1 - xt):.3f}, x_2={xt ** 2:.3f}]", f_out)
        trows = []
        for T in t_grid:
            _, _, eta, x_j, H, S, G, sro = point(xt, T)
            trows.append(([f"{T:.0f}", f"{x_j[0]:.4f}", f"{x_j[1]:.4f}", f"{x_j[2]:.4f}",
                            f"{H * 1000:+.3f}", f"{S * 1000:.4f}", f"{-T * S * 1000:+.3f}",
                            f"{G * 1000:+.3f}", f"{sro:+.4f}", f"{Delta_over_kT(N * d1, T):.2f}"], None))
        print_table(["T (K)", "x_0", "x_1", "x_2", "H_mix", "S_mix", "-T*S", "G_mix", "SRO",
                     "N*Delta_1/kT"], trows, f_out)

    # ------------------------------------------------------------------ [5]
    print_section('[5] SANITY DIAGNOSTICS', f_out)
    # Pure end-member limits reconfirmed against the REAL data just read (not a
    # standalone unit test) -- always runs, data-independent gate.
    t_check = t_grid[0]
    x_j_x0 = gqca_solver.cluster_fractions(gqca_solver.solve_eta(0.0, t_check, g, delta_e, N),
                                             t_check, g, delta_e, N)
    x_j_x1 = gqca_solver.cluster_fractions(gqca_solver.solve_eta(1.0, t_check, g, delta_e, N),
                                             t_check, g, delta_e, N)
    ok_x0 = abs(x_j_x0[0] - 1.0) < 1e-9 and abs(x_j_x0[1]) < 1e-9 and abs(x_j_x0[2]) < 1e-9
    ok_x1 = abs(x_j_x1[2] - 1.0) < 1e-9 and abs(x_j_x1[0]) < 1e-9 and abs(x_j_x1[1]) < 1e-9
    ok_sum = all(abs(sum(r[3]) - 1.0) < 1e-9 and min(r[3]) > -1e-12 for r in rows)
    ok_comp = all(abs(sum((j / N) * r[3][j] for j in range(3)) - r[0]) < 1e-6 for r in rows)
    ok_zero = abs(delta_e[0]) < 1e-12 and abs(delta_e[2]) < 1e-12
    ok_conv = all(d['scf_ok'] and (d['max_force'] is None or d['max_force'] <= args.force_tolerance)
                  for d in details)

    def status(flag, warn=False):
        return color_text("OK", 'green') if flag else color_text("WARN" if warn else "FAIL", 'yellow' if warn else 'red')

    print_dual(f"  x=0 pure-A limit (x_j = [1,0,0]) : {status(ok_x0)}", f_out)
    print_dual(f"  x=1 pure-B limit (x_j = [0,0,1]) : {status(ok_x1)}", f_out)
    print_dual(f"  sum(x_j) = 1 and x_j >= 0        : {status(ok_sum)}", f_out)
    print_dual(f"  sum(n_j/N * x_j) = x (mass bal.) : {status(ok_comp)}", f_out)
    print_dual(f"  Delta_0 = Delta_2 = 0 (reference): {status(ok_zero)}", f_out)
    print_dual(f"  all clusters SCF-converged/relaxed: {status(ok_conv, warn=True)}", f_out)

    # ------------------------------------------------------------------ [6]
    print_section('[6] INTERPRETATION', f_out)
    if abs(d1) < 1e-4:
        print_dual(f"* Delta_1 = {d1 * 1000:+.2f} meV/site ~ 0: unlike pairs cost nothing -> ideal-solution "
                    "behaviour, no preference for ordering or clustering.", f_out)
    elif d1 > 0:
        print_dual(f"* Delta_1 = {d1 * 1000:+.2f} meV/site > 0: AB (unlike) pairs are UNFAVOURABLE relative to "
                    "AA/BB -> tendency to CLUSTERING / phase separation (SRO > 0 at low T). Typical of "
                    "size-mismatched alloys.", f_out)
    else:
        print_dual(f"* Delta_1 = {d1 * 1000:+.2f} meV/site < 0: AB pairs are FAVOURABLE -> tendency to ORDERING / "
                    "compound formation (SRO < 0 at low T).", f_out)
    kt_lo, kt_hi = KB_EV_PER_K * t_grid[0], KB_EV_PER_K * t_grid[-1]
    print_dual(f"* N*|Delta_1|/kT (cluster energy over kT) goes from {N * abs(d1) / kt_lo:.1f} at "
                f"{t_grid[0]:.0f} K to {N * abs(d1) / kt_hi:.1f} at {t_grid[-1]:.0f} K: the larger it is, "
                "the stronger the deviation from random mixing.", f_out)
    xt0 = args.x_target[0]
    sro_lo = point(xt0, t_grid[0])[7]
    sro_hi = point(xt0, t_grid[-1])[7]
    print_dual(f"* SRO at x={xt0:g}: {sro_lo:+.3f} ({t_grid[0]:.0f} K) -> {sro_hi:+.3f} ({t_grid[-1]:.0f} K); "
                "it tends to 0 (random) as T grows.", f_out)
    g_vals = [point(xt0, T)[6] for T in t_grid]
    print_dual(f"* Delta_G_mix(x={xt0:g}) is {'negative' if max(g_vals) < 0 else 'not always negative'} over the "
                "sampled temperatures: " + ", ".join(f"{T:.0f} K: {G * 1000:+.1f} meV" for T, G in zip(t_grid, g_vals)) + ".",
                f_out)
    if omega > 0:
        print_dual(f"* Equivalent regular solution: Omega = {omega * 1000:.1f} meV, T_c ~ "
                    f"{omega / (2 * KB_EV_PER_K):.0f} K. Below T_c a real alloy of this system would tend to "
                    "phase-separate (thermodynamically); epitaxial/kinetically-frozen samples can stay mixed.", f_out)
    print_dual(color_text("* Compare Omega and H_mix(x=0.5) with experiment or SQS/cluster-expansion "
                           "literature for the same system (and level of theory) before drawing "
                           "quantitative conclusions.", 'cyan'), f_out)

    # ------------------------------------------------------------------ [7]
    print_section('[7] MISCIBILITY GAP (COMMON-TANGENT)', f_out)
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

    # ------------------------------------------------------------------ [8]
    print_section('[8] SUMMARY & FILES', f_out)

    csv_path = os.path.join(args.directory, f"{args.output}.csv")
    write_results_csv(csv_path, rows)
    print_dual(f"CSV (full grid)     : {csv_path}", f_out)

    clusters_csv = os.path.join(args.directory, f"{args.output}_clusters.csv")
    with open(clusters_csv, "w") as f:
        f.write("cluster,pair,sites,free_eng_eV,e_per_site_eV,delta_eV,v_per_site_A3,max_force_eV_A,scf_converged\n")
        for n in (0, 1, 2):
            d = details[n]
            v_txt = "" if d['v_site'] is None else f"{d['v_site']:.6f}"
            force_txt = "" if d['max_force'] is None else f"{d['max_force']:.6f}"
            f.write(f"cluster_n{n},{names[n].replace(' ', '')},{n_sites_of[n]},{raw_of[n]:.8f},"
                    f"{e_per_site[n]:.8f},{delta_e[n]:.8f},{v_txt},{force_txt},{int(d['scf_ok'])}\n")
    print_dual(f"CSV (per cluster)   : {clusters_csv}", f_out)
    if phase_csv_path:
        print_dual(f"Phase diagram       : {phase_csv_path}", f_out)

    dat_path = os.path.join(args.directory, f"{args.output}.dat")
    gplot_path = os.path.join(args.directory, f"{args.output}.gplot")
    write_spectrum_plot(dat_path, gplot_path, rows, t_grid)
    print_dual(f"Gnuplot (all T)     : {dat_path}, {gplot_path} "
                f"(cd {args.directory} && gnuplot {os.path.basename(gplot_path)})", f_out)

    t_ref = min(t_grid, key=lambda T: abs(T - args.plot_temp)) if args.plot_temp is not None else t_grid[0]
    if args.png or args.show:
        png_path = os.path.join(args.directory, f"{args.output}_summary.png") if args.png else None
        if args.show:
            print_dual("Dashboard          : opening on screen (close the window to continue)...", f_out)
        done = make_dashboard(rows, t_grid, t_ref, delta_e, e_per_site, details, sublattice, species_b,
                               png_path, args.show)
        if done and png_path:
            print_dual(f"Dashboard (PNG)     : {png_path}", f_out)
    else:
        print_dual("Dashboard          : not requested (add --show to open it, --png to save it)", f_out)

    if args.export_relaxed:
        out_dir = os.path.join(args.directory, "relaxed")
        os.makedirs(out_dir, exist_ok=True)
        ext = {"cif": "cif", "xsf": "xsf", "xyz": "xyz", "vasp": "vasp"}[args.export_relaxed]
        from ase.io import write as ase_write
        for n in (0, 1, 2):
            atoms, source = load_relaxed_atoms(cluster_dirs[n], args.file)
            if atoms is None:
                print_dual(color_text(f"[WARNING] Could not load the relaxed structure of cluster_n{n}.", 'yellow'), f_out)
                continue
            target = os.path.join(out_dir, f"cluster_n{n}.{ext}")
            ase_write(target, atoms, format=args.export_relaxed if ext != "vasp" else "vasp")
            print_dual(f"Relaxed structure   : {target}  [{source}]", f_out)

    if report_path:
        print_dual(f"Report              : {report_path}", f_out)

    print_dual(color_text(
        "\n[UNVERIFIED] S_mix and SRO use provisional formulas -- see core/gqca_solver.py's "
        "mixing_entropy/sro_parameter docstrings.", 'yellow'), f_out)

    if f_out:
        f_out.close()

    if args.view and sys.platform.startswith("linux") and not (os.environ.get("DISPLAY")
                                                              or os.environ.get("WAYLAND_DISPLAY")):
        print(color_text("[WARNING] --view needs a display (local desktop or `ssh -X`); skipped. "
                          "Use --export-relaxed to write the structures to files instead.", 'yellow'))
    elif args.view:
        from stb.core.ase_view import view_structure_interactive
        for n in (0, 1, 2):
            atoms, source = load_relaxed_atoms(cluster_dirs[n], args.file)
            if atoms is None:
                print(color_text(f"[WARNING] No structure to show for cluster_n{n}.", 'yellow'))
                continue
            print(f"\n[VIEW] cluster_n{n} -- {names[n]} -- {source} (close the window for the next one)")
            view_structure_interactive(atoms)

    print("\n[INFO] Complete job!")
    print("\n" + "-" * 60)
    print(color_text("GQCA mass-action analysis complete!\n", 'bold'))


def Delta_over_kT(delta, T):
    return delta / (KB_EV_PER_K * T)


if __name__ == "__main__":
    main()
