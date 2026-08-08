#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.5.0"  # New: [2b] PHYSICAL DIAGNOSTICS (dipole/spin moment/bond-length change),
                    # [5] SUGGESTED NEXT ANALYSES, and an opt-in [6] GIBBS FREE ENERGY (DG) PREP
                    # step (--compute-gibbs) that writes the displacement folders stb-adsorbGibbs
                    # (the new Stage 4) needs -- see their own docstrings/CLAUDE.md for the design.

import os
import re
import sys
import json
import glob
import shutil
import argparse
import numpy as np
import ase
import matplotlib.pyplot as plt
from datetime import datetime
from pymatgen.io.ase import AseAtomsAdaptor
from stb.core import siesta_log, structure_io
from stb.core.siesta_log import check_scf_and_force, report_quality_diagnostics
from stb.core.cli import color_text, show_intro, print_dual, print_section, print_table, capture_library_noise
from stb.core.ase_view import view_structure_interactive
from stb.core.pseudopotentials import copy_pseudo
from stb.core.adsorption_sites import (
    CONFIG_EXTRA_FILE, SPIN_POLARIZED_BLOCK, DIPOLE_CORRECTION_BLOCK, VDW_CORRECTION_BLOCK,
)
from stb.core.phonon_workflow import build_phonon_displacements, write_displacement_folders
from stb.adsorb import min_adsorbate_slab_distance

REPORT_FILE = "adsorption_report.txt"
SITES_REPORT_FILE = "adsorption_sites.txt"
RANKING_PLOT_FILE = "adsorption_ranking.png"
GIBBS_LOCAL_META_FILE = "gibbs_local_meta.json"
_HEIGHT_SUFFIX_RE = re.compile(r'_h[\d.]+$')
_LABEL_RE = re.compile(r'SystemLabel\s+\S+', re.IGNORECASE)

# Gibbs-prep displacement folders (disp_NNN/ or disp-NNN/) are independent,
# single-point finite-difference evaluations at a fixed, displaced geometry --
# letting SIESTA relax one on its own would silently move the atom back and
# invalidate the Hessian. Same config_extra.fdf %include mechanism (see
# structure_io.prepend_include) as adsorb_bsse.py's own
# BSSE_SINGLE_POINT_BLOCK for its ghost-fragment folders (duplicated, not
# imported -- same self-contained-workflow reasoning as force_system_label
# below), used in place of editing calc.fdf's flat text directly (the
# original design here) so this stage's folders are laid out the same way
# as the rest of the Adsorption workflow's (stb-adsorb, stb-adsorbBsse).
GIBBS_SINGLE_POINT_BLOCK = (
    "# Auto-generated -- keeps the cell fixed and forces a single-point SCF\n"
    "# (0 CG steps) so this displaced geometry isn't relaxed away.\n"
    "MD.VariableCell false\n"
    "MD.TypeOfRun       CG\n"
    "MD.Steps           0\n"
    "MD.NumCGsteps      0\n"
)

# Local-mode (partial-Hessian) displacement order (axis, sign) -- same
# convention/order as her_refs.py's/oer_refs.py's own _LOCAL_DISPLACEMENTS,
# applied here per local atom (see write_local_hessian_folders).
_LOCAL_DISPLACEMENTS = [(0, 1.0), (0, -1.0), (1, 1.0), (1, -1.0), (2, 1.0), (2, -1.0)]


def read_site_table(sites_root):
    """Parses the '# SITE_TABLE' section stb-adsorb writes in
    sites/adsorption_sites.txt into {label: (adsorbate_name, height)} --
    the source of truth for which adsorbate/<name>/ reference and which
    height each site belongs to (more robust than re-guessing either from
    the folder name by regex). Returns None if the report (or its table
    section) isn't found, so the caller can fall back to treating every
    site as belonging to a single, unnamed adsorbate -- e.g. a report the
    user deleted, or folders assembled by hand rather than via stb-adsorb.
    """
    report_path = os.path.join(sites_root, SITES_REPORT_FILE)
    if not os.path.isfile(report_path):
        return None
    table = {}
    in_table = False
    with open(report_path) as f:
        for line in f:
            if line.startswith("# SITE_TABLE"):
                in_table = True
                continue
            if not in_table:
                continue
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            label, ads_name = parts[0], parts[1]
            try:
                height = float(parts[2])
            except ValueError:
                continue
            table[label] = (ads_name, height)
    return table or None


def read_bsse_energy(bsse_dir, file_name):
    """Reads the two ghost-fragment reference energies stb-adsorb's
    --bsse-correction (default ON) writes under 'bsse/<site_label>/' -- a
    tree parallel to 'sites/', not nested inside a site's own folder --
    '<bsse_dir>/bsse_slab/<file_name>' and
    '<bsse_dir>/bsse_adsorbate/<file_name>' -- and returns (e_bsse_slab,
    e_bsse_adsorbate), or (None, None) if either is missing/unreadable.
    Best-effort per site (never blocks the uncorrected result for that
    site or any other). Callers are expected to have already confirmed
    '<bsse_dir>/bsse_slab' exists before calling this (main() does, to
    decide whether to even attempt reading BSSE energies for a site) --
    this function itself only worries about whether the energies inside
    are actually readable.
    """
    slab_path = os.path.join(bsse_dir, "bsse_slab", file_name)
    ads_path = os.path.join(bsse_dir, "bsse_adsorbate", file_name)
    e_slab = siesta_log.get_free_energy(slab_path)
    e_ads = siesta_log.get_free_energy(ads_path)
    if e_slab is None or e_ads is None:
        return None, None
    return e_slab, e_ads


def write_curve_plot(dat_path, rows):
    """Writes <dat_path> plus a companion .gplot (site index vs. adsorption
    energy, uncorrected and BSSE-corrected if available) -- same
    <name>.dat + <name>.gplot convention as the rest of the suite (e.g.
    convergence_analysis.py::write_curve_plot, phonons_pos.py's
    write_thermal_plots). Sites are discrete/categorical (not a continuous
    sweep), so this is a scatter, not a line.
    """
    with open(dat_path, 'w') as f:
        f.write("# Adsorption energy per candidate site\n")
        f.write("# 1:SiteIndex 2:E_ads(eV) 3:E_ads_BSSE(eV) 4:Adsorbate 5:Label\n")
        for i, row in enumerate(rows, start=1):
            bsse_str = f"{row.e_ads_bsse:.6f}" if row.e_ads_bsse is not None else "nan"
            f.write(f"{i}  {row.e_ads:.6f}  {bsse_str}  {row.ads_name}  {row.label}\n")

    base = os.path.splitext(dat_path)[0]
    base_name = os.path.basename(base)
    gplot_path = f"{base}.gplot"
    dat_name = os.path.basename(dat_path)
    has_bsse = any(row.e_ads_bsse is not None for row in rows)
    plot_lines = (f'"{dat_name}" using 1:2 with points pt 7 ps 2 lc rgb "#2255cc" title "E_{{ads}}"')
    if has_bsse:
        plot_lines += (f', "{dat_name}" using 1:3 with points pt 9 ps 2 lc rgb "#cc5522" '
                        'title "E_{ads} (BSSE)"')
    with open(gplot_path, 'w') as f:
        f.writelines([
            '# --- STB Plot Configuration ---\n',
            '# Generated by stb-adsorbAnalysis\n',
            'set terminal pdfcairo enhanced color font "Arial,14" size 7,5\n',
            f'set output "{base_name}.pdf"\n\n',
            'set title "Adsorption energy per site"\n',
            'set xlabel "Site index"\n',
            'set ylabel "E_{ads} (eV)"\n',
            'set grid\n',
            'set key top right\n',
            f'plot {plot_lines}\n',
        ])
    return gplot_path


def write_height_curve_plot(family, family_rows, out_dir):
    """Writes an E_ads-vs-height approach curve for one site "family"
    (same site + adsorbate, swept across --height-sweep) -- <family>.dat +
    <family>.gplot, same convention as write_curve_plot above. Only called
    for families with 2+ heights (see main()).
    """
    family_rows = sorted(family_rows, key=lambda r: r.height)
    dat_path = os.path.join(out_dir, f"height_curve_{family}.dat")
    with open(dat_path, 'w') as f:
        f.write(f"# Adsorption energy vs. height for site family '{family}'\n")
        f.write("# 1:Height(Ang) 2:E_ads(eV) 3:E_ads_BSSE(eV)\n")
        for row in family_rows:
            bsse_str = f"{row.e_ads_bsse:.6f}" if row.e_ads_bsse is not None else "nan"
            f.write(f"{row.height:.4f}  {row.e_ads:.6f}  {bsse_str}\n")

    base_name = f"height_curve_{family}"
    gplot_path = os.path.join(out_dir, f"{base_name}.gplot")
    dat_name = os.path.basename(dat_path)
    has_bsse = any(row.e_ads_bsse is not None for row in family_rows)
    plot_lines = (f'"{dat_name}" using 1:2 with linespoints lw 2 pt 7 lc rgb "#2255cc" title "E_{{ads}}"')
    if has_bsse:
        plot_lines += (f', "{dat_name}" using 1:3 with linespoints lw 2 pt 9 lc rgb "#cc5522" '
                        'title "E_{ads} (BSSE)"')
    with open(gplot_path, 'w') as f:
        f.writelines([
            '# --- STB Plot Configuration ---\n',
            '# Generated by stb-adsorbAnalysis\n',
            'set terminal pdfcairo enhanced color font "Arial,14" size 7,5\n',
            f'set output "{base_name}.pdf"\n\n',
            f'set title "Approach curve: {family}"\n',
            'set xlabel "Height (Ang)"\n',
            'set ylabel "E_{ads} (eV)"\n',
            'set grid\n',
            'set key top right\n',
            f'plot {plot_lines}\n',
        ])
    return dat_path, gplot_path


def plot_adsorption_ranking(rows, out_path, show=False):
    """Bar chart of E_ads (and E_ads_BSSE, where available) per site,
    ranked most stable first -- the matplotlib complement to
    write_curve_plot's gnuplot .dat/.gplot pair (kept as-is; this is
    additive, not a replacement), same single-series-bar-chart style as
    adsorb.py::plot_ml_ranking. `rows` is the same list `main()` already
    sorted by `e_ads` (ascending -- most stable first).

    `show=True` (--view-plots) additionally blocks on plt.show() before
    the figure is closed -- same convention as adsorb.py's --view-plots.
    """
    labels = [f"{r.ads_name or 'default'}\n{r.label}" for r in rows]
    e_ads = [r.e_ads for r in rows]
    has_bsse = [r.e_ads_bsse is not None for r in rows]
    x = np.arange(len(rows))
    width = 0.35 if any(has_bsse) else 0.6
    fig, ax = plt.subplots(figsize=(max(6, 0.9 * len(labels)), 5))
    colors = ['tab:green' if i == 0 else 'tab:blue' for i in range(len(rows))]
    ax.bar(x - width / 2 if any(has_bsse) else x, e_ads, width=width, color=colors, label="E_ads")
    if any(has_bsse):
        bsse_vals = [r.e_ads_bsse if r.e_ads_bsse is not None else np.nan for r in rows]
        ax.bar(x + width / 2, bsse_vals, width=width, color='tab:orange', label="E_ads (BSSE)")
        ax.legend(fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("E_ads (eV)")
    ax.set_title("Adsorption energy per site (most stable first)")
    ax.axhline(0.0, color='gray', linestyle='--', linewidth=0.8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    if show:
        plt.show()
    plt.close(fig)


def read_site_geometry_atoms(folder):
    """Returns an ase.Atoms for whatever geometry is available in `folder`:
    prefers a finished SIESTA run's own '<label>.XV' (the actual relaxed
    geometry SIESTA used, converted the same way mlff_analysis.py already
    does for a sisl Geometry), and falls back to 'structure.fdf' (the
    input Stage 1 always writes, via this suite's own core/structure_io.py
    reader) when no .XV exists yet -- e.g. before SIESTA has been run in
    that folder at all. Returns None if neither is readable, so callers
    can skip that frame instead of crashing --view over one missing/
    unrun folder.
    """
    xv_files = sorted(glob.glob(os.path.join(folder, "*.XV")))
    if xv_files:
        from stb.core.deps import require_sisl
        sisl = require_sisl()
        geom = sisl.get_sile(xv_files[0]).read_geometry()
        return ase.Atoms(symbols=[a.symbol for a in geom.atoms],
                          positions=geom.xyz, cell=np.array(geom.cell), pbc=True)
    fdf_path = os.path.join(folder, "structure.fdf")
    if os.path.isfile(fdf_path):
        pmg_structure = structure_io.to_pymatgen(structure_io.read_fdf(fdf_path))
        return AseAtomsAdaptor.get_atoms(pmg_structure)
    return None


def force_system_label(calc_text, label):
    """Substitutes/appends SystemLabel -- every derived Gibbs-prep folder
    needs its own distinct label so SIESTA's per-run files (.FA, .DM, ...)
    never collide between folders, and so stb-adsorbGibbs knows exactly
    which '<label>.FA' filename to read back. Duplicated from her_refs.py/
    oer_refs.py's own force_system_label (each workflow re-reads/derives
    its own folders rather than importing a sibling workflow's module --
    same "self-contained stage" convention HER/OER already follow; the
    within-Adsorption-workflow cross-imports elsewhere, e.g.
    adsorb_bsse.py importing read_site_table from this module, don't apply
    here since that would create a real import cycle with adsorb_bsse.py).
    """
    new_text, count = _LABEL_RE.subn(f'SystemLabel {label}', calc_text)
    if count == 0:
        new_text += f"\nSystemLabel {label}\n"
    return new_text


def read_site_theory_flags(folder):
    """Returns (spin_polarized, dipole_corrected, vdw_corrected): whether
    `folder`'s own config_extra.fdf (written by stb-adsorb's
    write_reference_folder) forced Spin polarized / Slab.DipoleCorrection /
    the DFT-D3 dispersion correction. The Gibbs-prep displacement folders
    below must inherit the SAME level of theory as the structure they're
    perturbing -- forgetting this was a real, user-reported bug already
    fixed once for stb-adsorbBsse's own ghost fragments (see CLAUDE.md);
    the same care applies here, including for DFT-D3 (forced
    unconditionally on every stb-adsorb folder, see adsorb.py's own
    force_vdw=True call sites). Returns (False, False, False) if `folder`
    has no config_extra.fdf at all.
    """
    path = os.path.join(folder, CONFIG_EXTRA_FILE)
    if not os.path.isfile(path):
        return False, False, False
    with open(path) as f:
        text = f.read()
    return "Spin" in text, "Slab.DipoleCorrection" in text, "DFTD3" in text


def read_original_calc_text(calc_fdf_path):
    """Strips the '%include config_extra.fdf' sidecar line
    structure_io.prepend_include added when stb-adsorb first wrote this
    folder's own calc.fdf, recovering the original --calc template text --
    needed so the Gibbs-prep folders below can prepend their OWN config_
    extra.fdf (single-point enforcement + whatever Spin/dipole flags
    read_site_theory_flags found) without doubling the include line.
    Duplicated from adsorb_bsse.py's identical helper (see
    force_system_label's docstring for why: importing it directly would
    create an adsorb_bsse.py <-> adsorb_analysis.py import cycle, since
    adsorb_bsse.py already imports read_site_table from this module). A
    no-op if the file doesn't start with that exact prefix (defensive;
    e.g. a hand-edited calc.fdf, or the isolated-adsorbate reference,
    which never gets a config_extra.fdf include of its own spin override
    -- see adsorb.py's force_spin_polarized).
    """
    with open(calc_fdf_path) as f:
        text = f.read()
    prefix = f"%include {CONFIG_EXTRA_FILE}\n\n"
    if text.startswith(prefix):
        return text[len(prefix):]
    return text


def write_gibbs_folder(out_dir, fdf_structure, calc_text, pp_path, config_extra_content=None):
    """Writes structure.fdf + calc.fdf + copied pseudos for one Gibbs-prep
    displacement folder. Same minimal writer as her_refs.py/oer_refs.py's
    own write_folder (duplicated, not imported -- same reasoning as
    force_system_label above). When `config_extra_content` is given, it's
    written to config_extra.fdf and calc.fdf becomes '%include
    config_extra.fdf' + calc_text (structure_io.prepend_include) instead of
    calc_text written bare -- same split adsorb_bsse.py's write_bsse_folders
    uses for its own ghost-fragment folders.
    """
    os.makedirs(out_dir, exist_ok=True)
    structure_io.write_fdf(fdf_structure, os.path.join(out_dir, "structure.fdf"))
    if config_extra_content is not None:
        with open(os.path.join(out_dir, CONFIG_EXTRA_FILE), "w") as f:
            f.write(config_extra_content)
        calc_text = structure_io.prepend_include(calc_text, CONFIG_EXTRA_FILE)
    with open(os.path.join(out_dir, "calc.fdf"), "w") as f:
        f.write(calc_text)
    symbols = sorted({symbol for symbol, _ in fdf_structure.atoms})
    for symbol in symbols:
        copy_pseudo(pp_path, symbol, out_dir)


def write_local_hessian_folders(gibbs_dir, base_structure, local_indices, displacement_ang,
                                 calc_text, pp_path, system_label, config_extra_content=None):
    """Writes 6*len(local_indices) single-point folders (each local atom
    displaced +/-x/+/-y/+/-z by `displacement_ang` from its position in
    `base_structure`, every other atom -- including OTHER local atoms --
    held fixed) for the 'local' partial-Hessian Gibbs-prep mode: a
    decoupled-oscillator approximation for whichever atoms are "local"
    (the adsorbate's own atoms, substrate frozen, for the site side; every
    atom, for the isolated-adsorbate reference side, since there's no
    substrate to freeze there at all). Same N-atom generalization
    oer_refs.py::write_local_zpe_folders already established over
    her_refs.py's original single-atom (H*) version -- duplicated here
    (not imported) for the same self-contained-workflow reason as
    force_system_label above. Writes a JSON sidecar (gibbs_local_meta.json)
    with the local atom indices/symbols (per-atom mass lookup in Stage 4),
    displacement_ang, the exact (atom_index, axis, sign) order, and the
    SystemLabel used, so stb-adsorbGibbs never has to guess a folder-name
    -to-DOF mapping or a '<label>.FA' filename.
    """
    os.makedirs(gibbs_dir, exist_ok=True)
    inv_lattice = np.linalg.inv(base_structure.lattice)
    local_symbols = [base_structure.atoms[i][0] for i in local_indices]

    order = []
    for atom_index in local_indices:
        for axis, sign in _LOCAL_DISPLACEMENTS:
            order.append({"atom_index": atom_index, "axis": axis, "sign": sign})

    for i, entry in enumerate(order, start=1):
        atom_index, axis, sign = entry["atom_index"], entry["axis"], entry["sign"]
        symbol, frac = base_structure.atoms[atom_index]
        cart = frac @ base_structure.lattice
        delta_cart = np.zeros(3)
        delta_cart[axis] = sign * displacement_ang
        new_frac = (cart + delta_cart) @ inv_lattice
        new_atoms = list(base_structure.atoms)
        new_atoms[atom_index] = (symbol, new_frac)
        disp_structure = structure_io.FdfStructure(
            lattice=base_structure.lattice, lattice_constant=base_structure.lattice_constant,
            species=base_structure.species, species_meta=base_structure.species_meta,
            atoms=new_atoms, coord_format=base_structure.coord_format, raw_lines=[],
        )
        disp_dir = os.path.join(gibbs_dir, f"disp_{i:03d}")
        write_gibbs_folder(disp_dir, disp_structure, calc_text, pp_path,
                            config_extra_content=config_extra_content)

    with open(os.path.join(gibbs_dir, GIBBS_LOCAL_META_FILE), "w") as f:
        json.dump({
            "local_indices": list(local_indices),
            "local_symbols": local_symbols,
            "displacement_ang": displacement_ang,
            "order": order,
            "system_label": system_label,
        }, f)


class SiteRow:
    """One analyzed site: label, which adsorbate/height it belongs to, and
    its computed energies. A plain attribute-holder (not a dataclass, to
    keep this a light dependency-free module) instead of a bare tuple --
    the field count grew enough (multi-adsorbate + height-sweep) that
    positional tuple unpacking everywhere had become error-prone to edit.
    """
    def __init__(self, label, ads_name, height, e_site, e_ads, e_ads_bsse, scf_ok, max_force,
                 dipole=None, spin_moment=None, bond_change=None,
                 e_bsse_slab=None, e_bsse_ads=None):
        self.label = label
        self.ads_name = ads_name
        self.height = height
        self.e_site = e_site
        self.e_ads = e_ads
        self.e_ads_bsse = e_ads_bsse
        self.scf_ok = scf_ok
        self.max_force = max_force
        self.dipole = dipole              # |Electric dipole|, a.u., or None
        self.spin_moment = spin_moment    # net magnetic moment, Bohr magnetons, or None
        self.bond_change = bond_change    # relaxed - initial closest adsorbate-slab distance, Ang, or None
        self.e_bsse_slab = e_bsse_slab    # raw ghost-fragment energy, real slab + ghost adsorbate, eV, or None
        self.e_bsse_ads = e_bsse_ads      # raw ghost-fragment energy, ghost slab + real adsorbate, eV, or None


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Computes adsorption energies from an stb-adsorb sweep: "
        "E_ads = E_site - E_clean_slab - E_adsorbate.", 'bold')}""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s\n"
               "  %(prog)s --dir adsorption_run --file calc.out\n"
               "  %(prog)s --apply best_structure.fdf\n"
    )

    parser.add_argument("--dir", type=str, default="adsorption_run",
                         help="Root directory containing 'clean_slab/', 'adsorbate*/' and 'sites/' "
                              "(default: adsorption_run -- stb-adsorb's own default --output-dir, "
                              "so this stage runs against that one's output with no extra flag in "
                              "the common case).")
    parser.add_argument("--file", type=str, default="calc.out",
                         help="SIESTA output filename inside each folder (default: calc.out).")
    parser.add_argument("-o", "--output", type=str, default="adsorption_curve.dat",
                         help="Output data file name (default: adsorption_curve.dat).")
    parser.add_argument("--apply", type=str, default=None, metavar="STRUCTURE_FDF",
                         help="Copy the most stable site's structure.fdf (BSSE-corrected ranking "
                             "when available for every site, else the uncorrected ranking) to "
                             "this path.")
    parser.add_argument("--force-tolerance", type=float, default=0.05,
                         help="Residual atomic force in eV/Ang (default: 0.05, same as "
                              "stb-cohesiveAnalysis) above which a site's calc.out is flagged as "
                              "possibly not relaxed -- E_ads would then reflect a strained/off-"
                              "equilibrium geometry rather than the true adsorption minimum. "
                              "Advisory only, never blocks the result.")
    parser.add_argument("--save-report", action="store_true",
                         help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--view", action="store_true",
                         help="Open clean_slab, every adsorbate reference and every site "
                              "(preferring each folder's own <label>.XV -- the actual relaxed "
                              "geometry, if SIESTA has been run there -- falling back to its "
                              "input structure.fdf otherwise) in ASE's interactive 3D viewer as "
                              "a multi-frame browser, after everything else has finished. "
                              "Needs a display.")
    parser.add_argument("--view-plots", action="store_true",
                         help="Also show the matplotlib adsorption-ranking plot on screen "
                              "(one blocking window), instead of only saving it as a PNG. Off "
                              "by default.")
    parser.add_argument("--compute-gibbs", action="store_true",
                         help="Also write the vibrational-Hessian displacement folders (under "
                              "'<dir>/gibbs/') stb-adsorbGibbs (Stage 4) needs to compute a "
                              "Gibbs free energy of adsorption, for the winning site (BSSE-"
                              "ranked if complete, else uncorrected) and its isolated-adsorbate "
                              "reference. Off by default -- this is the one exception in this "
                              "tool to 'analysis never generates new folders', by explicit "
                              "design (see CLAUDE.md): the alternative would be a whole extra "
                              "prep-only stage duplicating everything [0]-[4] already do.")
    parser.add_argument("--zpe-mode", choices=["local", "full"], default="local",
                         help="Vibrational Hessian mode for the SITE side, when --compute-gibbs "
                              "(default: local). 'local': partial Hessian, only the adsorbate's "
                              "own atoms displaced (substrate frozen) -- fast, a decoupled-"
                              "oscillator approximation, no Phonopy involved. 'full': a real "
                              "Phonopy phonon calculation of the ENTIRE site structure, minus "
                              "the clean slab's own (also full-structure) phonon calculation -- "
                              "rigorous, but far more SIESTA jobs. The isolated-adsorbate "
                              "reference ALWAYS uses the 'local' (full-Hessian-of-the-whole-"
                              "molecule) method regardless of this flag -- there's no periodic "
                              "substrate phonon to subtract for an isolated molecule in a vacuum "
                              "box, so Phonopy's supercell machinery buys nothing there.")
    parser.add_argument("--displacement", type=float, default=0.015,
                         help="Finite-difference displacement in Ang for the Hessian/Phonopy "
                              "displacements (default: 0.015, same as stb-herRefs/stb-oerRefs).")
    parser.add_argument("--gibbs-supercell", type=int, nargs=3, default=[1, 1, 1],
                         help="Supercell dimensions for --zpe-mode full's site+clean-slab phonon "
                              "calculations (default: 1 1 1 -- appropriate when the slab's own "
                              "in-plane supercell is already large enough).")
    parser.add_argument("-v", "--version", action="version", version=f"stb-adsorbAnalysis {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("Analyze an adsorption-energy study:", 'bold'))
    print("-" * 60)

    clean_slab_out = os.path.join(args.dir, "clean_slab", args.file)
    sites_root = os.path.join(args.dir, "sites")
    bsse_root = os.path.join(args.dir, "bsse")

    if not os.path.isdir(sites_root):
        print(color_text(f"[ERROR] '{sites_root}' not found. Did you run stb-adsorb?", 'red'))
        sys.exit(1)

    site_table = read_site_table(sites_root)

    library_warnings = []  # collected via capture_library_noise, reported in [5]
    view_frames = []  # (label, ase.Atoms), only populated when --view is given

    report_path = REPORT_FILE if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    print_dual(f"{color_text('===== ADSORPTION ENERGY REPORT =====', 'magenta')}", f_out)

    print_section('[0] RUN METADATA', f_out)
    print_dual(f"Date/time  : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", f_out)
    print_dual(f"Directory  : {args.dir}", f_out)
    print_dual(f"Output file: {args.file}", f_out)
    if site_table is None:
        print_dual(color_text(
            "[NOTE] No 'sites/adsorption_sites.txt' site table found -- treating every site "
            "as a single unnamed adsorbate (pre-multi-adsorbate/height-sweep layout, or the "
            "report was deleted).", 'yellow'), f_out)

    print_section('[1] REFERENCE ENERGIES', f_out)
    e_clean_slab = siesta_log.get_free_energy(clean_slab_out)
    if e_clean_slab is None:
        print_dual(color_text(f"[ERROR] Could not read energy from '{clean_slab_out}'.", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)
    print_dual(f"E_clean_slab : {e_clean_slab:.6f} eV  ({clean_slab_out})", f_out)
    report_quality_diagnostics("clean_slab", clean_slab_out, args.force_tolerance, f_out)
    clean_slab_fdf_path = os.path.join(args.dir, "clean_slab", "structure.fdf")
    n_substrate = len(structure_io.read_fdf(clean_slab_fdf_path).atoms) \
        if os.path.isfile(clean_slab_fdf_path) else 0
    if args.view:
        with capture_library_noise(library_warnings, "sisl/ASE (clean_slab geometry)"):
            atoms = read_site_geometry_atoms(os.path.join(args.dir, "clean_slab"))
        if atoms is not None:
            view_frames.append(("clean_slab", atoms))

    # Adsorbate names present, from the site table if available, else a
    # single unnamed one mapping to the legacy 'adsorbate/' folder.
    ads_names = sorted({name for name, _h in site_table.values()}) if site_table else [None]
    adsorbate_energies = {}  # name -> energy
    adsorbate_dirs = {}  # name -> folder, reused by [6] GIBBS FREE ENERGY (DG) PREP below
    for name in ads_names:
        candidate_dirs = []
        if name:
            candidate_dirs.append(os.path.join(args.dir, f"adsorbate_{name}"))
        candidate_dirs.append(os.path.join(args.dir, "adsorbate"))
        ads_dir = next((d for d in candidate_dirs if os.path.isdir(d)), None)
        if ads_dir is None:
            print_dual(color_text(
                f"[ERROR] No adsorbate reference folder found for "
                f"'{name or '(default)'}' (tried: {', '.join(candidate_dirs)}).", 'red'), f_out)
            if f_out:
                f_out.close()
            sys.exit(1)
        ads_out_path = os.path.join(ads_dir, args.file)
        e_ads = siesta_log.get_free_energy(ads_out_path)
        if e_ads is None:
            print_dual(color_text(f"[ERROR] Could not read energy from "
                                   f"'{ads_out_path}'.", 'red'), f_out)
            if f_out:
                f_out.close()
            sys.exit(1)
        adsorbate_energies[name] = e_ads
        adsorbate_dirs[name] = ads_dir
        print_dual(f"E_adsorbate ({name or 'default'}) : {e_ads:.6f} eV  ({ads_dir})", f_out)
        report_quality_diagnostics(f"adsorbate ({name or 'default'})", ads_out_path,
                                    args.force_tolerance, f_out)
        if args.view:
            with capture_library_noise(library_warnings, f"sisl/ASE (adsorbate {name or 'default'} geometry)"):
                atoms = read_site_geometry_atoms(ads_dir)
            if atoms is not None:
                view_frames.append((f"adsorbate({name or 'default'})", atoms))

    print_section('[2] SITE RESULTS: CONFIGURATION COUNT & TABLE', f_out)
    site_dirs = sorted(
        d for d in os.listdir(sites_root)
        if os.path.isdir(os.path.join(sites_root, d)) and d.startswith("site_")
    )
    if not site_dirs:
        print_dual(color_text(f"[ERROR] No 'site_*' folders found in '{sites_root}'.", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    rows = []
    n_skip_missing = 0
    n_skip_unparseable = 0
    n_bsse_found = 0
    incomplete_bsse_labels = []
    scf_warn_labels = []
    force_warn_labels = []
    bsse_scf_warn_labels = []
    table_rows = []
    for label in site_dirs:
        site_dir = os.path.join(sites_root, label)
        ads_name, height = site_table.get(label, (None, None)) if site_table else (None, None)

        if args.view:
            with capture_library_noise(library_warnings, f"sisl/ASE (site {label} geometry)"):
                atoms = read_site_geometry_atoms(site_dir)
            if atoms is not None:
                view_frames.append((label, atoms))

        out_path = os.path.join(site_dir, args.file)
        if not os.path.exists(out_path):
            n_skip_missing += 1
            table_rows.append(([label, ads_name or '--', "--" if height is None else f"{height:.2f}",
                                 "SKIP (missing calc.out)", "--", "--", "--"], 'yellow'))
            continue
        e_site = siesta_log.get_free_energy(out_path)
        if e_site is None:
            n_skip_unparseable += 1
            table_rows.append(([label, ads_name or '--', "--" if height is None else f"{height:.2f}",
                                 "SKIP (unparseable energy)", "--", "--", "--"], 'yellow'))
            continue
        e_ads = e_site - e_clean_slab - adsorbate_energies[ads_name]
        scf_ok, max_force = check_scf_and_force(out_path)
        if not scf_ok:
            scf_warn_labels.append(label)
        if max_force is not None and max_force > args.force_tolerance:
            force_warn_labels.append(label)

        e_ads_bsse = None
        e_bsse_slab = None
        e_bsse_ads = None
        bsse_dir = os.path.join(bsse_root, label)
        if os.path.isdir(os.path.join(bsse_dir, "bsse_slab")):
            e_bsse_slab, e_bsse_ads = read_bsse_energy(bsse_dir, args.file)
            if e_bsse_slab is not None and e_bsse_ads is not None:
                e_ads_bsse = e_site - e_bsse_slab - e_bsse_ads
                n_bsse_found += 1
                bsse_slab_scf, _f1 = check_scf_and_force(
                    os.path.join(bsse_dir, "bsse_slab", args.file))
                bsse_ads_scf, _f2 = check_scf_and_force(
                    os.path.join(bsse_dir, "bsse_adsorbate", args.file))
                if not (bsse_slab_scf and bsse_ads_scf):
                    bsse_scf_warn_labels.append(label)
            else:
                incomplete_bsse_labels.append(label)

        # [2b] PHYSICAL DIAGNOSTICS -- best-effort, never fatal to the main
        # E_ads result: dipole/spin come straight from calc.out (None if
        # SIESTA never printed the line, e.g. a non-polarized run for spin);
        # bond_change compares the closest adsorbate-slab distance at the
        # pre-relaxation guess (structure.fdf) vs. the relaxed geometry
        # (preferring .XV, same as read_site_geometry_atoms's own default).
        dipole_vec = siesta_log.get_electric_dipole(out_path)
        dipole_mag = float(np.linalg.norm(dipole_vec)) if dipole_vec is not None else None
        spin_moment = siesta_log.get_spin_moment(out_path)
        bond_change = None
        if n_substrate > 0:
            try:
                initial_pmg = structure_io.to_pymatgen(
                    structure_io.read_fdf(os.path.join(site_dir, "structure.fdf")))
                dist_initial = min_adsorbate_slab_distance(initial_pmg, n_substrate)
                relaxed_atoms = read_site_geometry_atoms(site_dir)
                if relaxed_atoms is not None and dist_initial is not None:
                    relaxed_pmg = AseAtomsAdaptor.get_structure(relaxed_atoms)
                    dist_relaxed = min_adsorbate_slab_distance(relaxed_pmg, n_substrate)
                    if dist_relaxed is not None:
                        bond_change = dist_relaxed - dist_initial
            except Exception:
                bond_change = None

        rows.append(SiteRow(label, ads_name, height, e_site, e_ads, e_ads_bsse, scf_ok, max_force,
                             dipole=dipole_mag, spin_moment=spin_moment, bond_change=bond_change,
                             e_bsse_slab=e_bsse_slab, e_bsse_ads=e_bsse_ads))
        row_color = 'yellow' if (not scf_ok or label in force_warn_labels) else None
        table_rows.append(([
            label, ads_name or '--', "--" if height is None else f"{height:.2f}",
            f"{e_ads:.6f}", f"{e_ads_bsse:.6f}" if e_ads_bsse is not None else "--",
            "OK" if scf_ok else "WARN", f"{max_force:.4f}" if max_force is not None else "--",
        ], row_color))

    n_skipped = n_skip_missing + n_skip_unparseable
    n_bsse_absent = len(site_dirs) - n_bsse_found - len(incomplete_bsse_labels)
    print_dual(f"Site folders found  : {len(site_dirs)}", f_out)
    print_dual(f"Read successfully   : {len(rows)}", f_out)
    print_dual(f"Skipped             : {n_skipped}  (missing {args.file}: {n_skip_missing}, "
                f"unparseable energy: {n_skip_unparseable})", f_out)
    if len(ads_names) > 1:
        print_dual("Per-adsorbate breakdown:", f_out)
        for name in ads_names:
            n_found_name = sum(1 for _l in site_dirs
                                if (site_table.get(_l, (None, None))[0] if site_table else None) == name)
            n_read_name = sum(1 for r in rows if r.ads_name == name)
            print_dual(f"  {name or '(default)'}: {n_found_name} found, {n_read_name} read", f_out)
    families = {}
    for label in site_dirs:
        families.setdefault(_HEIGHT_SUFFIX_RE.sub('', label), []).append(label)
    swept_families = {fam: labels for fam, labels in families.items() if len(labels) > 1}
    if swept_families:
        print_dual(f"Height-sweep families: {len(swept_families)} site(s) swept across "
                    f"multiple heights (family -> heights found):", f_out)
        for fam, labels in swept_families.items():
            print_dual(f"  {fam}: {len(labels)} height(s)", f_out)
    print_dual(f"BSSE coverage       : complete {n_bsse_found}, incomplete "
                f"{len(incomplete_bsse_labels)}, absent {n_bsse_absent}  (of {len(site_dirs)} site(s))", f_out)

    print_dual(color_text(
        "\n[BSSE PHYSICS CHECK] E_ads(BSSE) = E_site - E_bsse_slab - E_bsse_adsorbate, where "
        "E_bsse_slab/E_bsse_adsorbate are single-point ghost-fragment energies read from "
        "'bsse/<site>/bsse_slab/' and 'bsse/<site>/bsse_adsorbate/' -- written by stb-adsorbBsse "
        "(Stage 2) at the site's actual RELAXED geometry (from its finished siesta.XV), not the "
        "pre-relaxation guess, since each fragment must be evaluated at the same geometry as the "
        "real, relaxed site for the correction to mean anything (the standard Boys-Bernardi "
        "counterpoise, applied to both the substrate and the adsorbate). Full explanation: "
        "examples/4.8-adsorption/README.md.", 'cyan'), f_out)

    print_dual("", f_out)
    headers = ["Site", "Adsorbate", "Height", "E_ads(eV)", "E_ads_BSSE(eV)", "SCF", "MaxF(eV/A)"]
    print_table(headers, table_rows, f_out)

    if incomplete_bsse_labels:
        print_dual(color_text(
            f"[WARNING] {len(incomplete_bsse_labels)} site(s) have a 'bsse_slab/'/"
            "'bsse_adsorbate/' folder but incomplete/unreadable results -- BSSE-corrected "
            f"energy skipped for: {', '.join(incomplete_bsse_labels)}.", 'yellow'), f_out)
    if scf_warn_labels:
        print_dual(color_text(
            f"[WARNING] {len(scf_warn_labels)} site(s) never confirmed SCF convergence -- "
            f"their E_ads may be unreliable: {', '.join(scf_warn_labels)}.", 'yellow'), f_out)
    if force_warn_labels:
        print_dual(color_text(
            f"[WARNING] {len(force_warn_labels)} site(s) have residual force above "
            f"--force-tolerance ({args.force_tolerance} eV/Ang), possibly not relaxed: "
            f"{', '.join(force_warn_labels)}.", 'yellow'), f_out)
    if bsse_scf_warn_labels:
        print_dual(color_text(
            f"[WARNING] {len(bsse_scf_warn_labels)} site(s)' BSSE ghost-fragment calculation(s) "
            f"never confirmed SCF convergence -- their E_ads_BSSE may be unreliable: "
            f"{', '.join(bsse_scf_warn_labels)}.", 'yellow'), f_out)

    if not rows:
        print_dual(color_text("\n[ERROR] No valid site results found.", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_section('[2b] PHYSICAL DIAGNOSTICS', f_out)
    print_dual(color_text(
        "Dipole/spin come straight from each site's calc.out (Slab.DipoleCorrection and "
        "Spin polarized are already forced by stb-adsorb by default, so these should be "
        "non-trivial for a genuine one-sided adsorption site -- '--' means SIESTA never "
        "printed the line, e.g. a non-polarized run for the spin column). Bond change compares "
        "the closest adsorbate-slab distance at the pre-relaxation guess vs. the relaxed "
        "geometry -- positive means the adsorbate moved AWAY from the slab during relaxation, "
        "negative means it moved closer (the common case for a real chemisorption bond).",
        'cyan'), f_out)
    diag_rows = [([
        r.label, f"{r.dipole:.4f}" if r.dipole is not None else "--",
        f"{r.spin_moment:+.4f}" if r.spin_moment is not None else "--",
        f"{r.bond_change:+.4f}" if r.bond_change is not None else "--",
    ], None) for r in rows]
    print_table(["Site", "Dipole(a.u.)", "MagMoment(muB)", "Bond Change(Ang)"], diag_rows, f_out)

    bsse_breakdown_rows = [r for r in rows if r.e_ads_bsse is not None]
    if bsse_breakdown_rows:
        print_section('[2c] BSSE BREAKDOWN', f_out)
        print_dual(color_text(
            "Raw ghost-fragment energies behind each site's E_ads_BSSE above (see the "
            "[BSSE PHYSICS CHECK] note below the [2] table for the full formula) -- Shift is "
            "E_ads_BSSE - E_ads (the counterpoise correction itself; positive means the "
            "correction made binding LESS favorable, the expected direction), and Shift(%) is "
            "that same shift as a fraction of |E_ads| (large values flag a site where BSSE is a "
            "significant fraction of the whole binding energy, not just a small correction).",
            'cyan'), f_out)
        bsse_breakdown_table = []
        for r in bsse_breakdown_rows:
            shift = r.e_ads_bsse - r.e_ads
            shift_pct = (shift / abs(r.e_ads) * 100.0) if r.e_ads != 0 else float('nan')
            bsse_breakdown_table.append(([
                r.label, f"{r.e_bsse_slab:.6f}" if r.e_bsse_slab is not None else "--",
                f"{r.e_bsse_ads:.6f}" if r.e_bsse_ads is not None else "--",
                f"{shift:+.6f}", f"{shift_pct:+.2f}",
            ], None))
        print_table(["Site", "E_bsse_slab(eV)", "E_bsse_adsorbate(eV)", "Shift(eV)", "Shift(%)"],
                     bsse_breakdown_table, f_out)

    # Primary ranking always uses the uncorrected E_ads -- it's available
    # for every row, so it's the only metric safe to compare across all
    # sites uniformly (it's already referenced against each row's own
    # adsorbate energy, so comparing across different adsorbates is
    # meaningful too, unlike stb-adsorb --ml-rank's raw MACE energies).
    # BSSE-corrected values are reported alongside (and separately
    # ranked below, when complete for every site) rather than silently
    # mixed into the same sort.
    rows.sort(key=lambda r: r.e_ads)

    print_section('[3] SUMMARY & PLOT', f_out)
    print_dual(f"Sites analyzed : {len(rows)} (skipped: {n_skipped})", f_out)
    best = rows[0]
    verdict = "exothermic (favorable)" if best.e_ads < 0 else "endothermic (unfavorable)"
    print_dual(f"{color_text('Most stable site (uncorrected):', 'green')} {best.label}  "
                f"(E_ads = {best.e_ads:.6f} eV, {verdict})", f_out)

    apply_source_label = best.label
    if n_bsse_found == 0:
        print_dual(color_text(
            "\n[NOTE] No BSSE-corrected results found -- run stb-adsorbBsse (Stage 2) once these "
            "sites have finished relaxing, then run SIESTA in the 'bsse/site_*/' folders it "
            "writes, for a corrected reference.", 'yellow'), f_out)
    else:
        bsse_rows = [r for r in rows if r.e_ads_bsse is not None]
        if len(bsse_rows) == len(rows):
            bsse_ranked = sorted(bsse_rows, key=lambda r: r.e_ads_bsse)
            bsse_best = bsse_ranked[0]
            bsse_verdict = "exothermic (favorable)" if bsse_best.e_ads_bsse < 0 else "endothermic (unfavorable)"
            print_dual(f"{color_text('Most stable site (BSSE-corrected):', 'green')} "
                        f"{bsse_best.label}  (E_ads = {bsse_best.e_ads_bsse:.6f} eV, {bsse_verdict})", f_out)
            shift = bsse_best.e_ads_bsse - bsse_best.e_ads
            print_dual(f"BSSE correction at that site: {shift:+.6f} eV (uncorrected LCAO "
                        "adsorption energies systematically over-bind -- expect this to make "
                        "the energy less negative)", f_out)
            apply_source_label = bsse_best.label
        else:
            print_dual(color_text(
                f"\n[NOTE] BSSE-corrected energy available for only {len(bsse_rows)}/{len(rows)} "
                "site(s) -- not reporting a BSSE-ranked \"most stable\" until all sites have it.",
                'yellow'), f_out)

    if len(ads_names) > 1:
        print_dual(f"\n{color_text('Best site per adsorbate (uncorrected):', 'cyan')}", f_out)
        for name in ads_names:
            per_ads = [r for r in rows if r.ads_name == name]
            if not per_ads:
                continue
            best_per_ads = min(per_ads, key=lambda r: r.e_ads)
            print_dual(f"  {name}: {best_per_ads.label} (E_ads = {best_per_ads.e_ads:.6f} eV)", f_out)
            per_ads_bsse = [r for r in per_ads if r.e_ads_bsse is not None]
            if per_ads_bsse and len(per_ads_bsse) == len(per_ads):
                best_per_ads_bsse = min(per_ads_bsse, key=lambda r: r.e_ads_bsse)
                print_dual(f"     BSSE-corrected: {best_per_ads_bsse.label} "
                            f"(E_ads = {best_per_ads_bsse.e_ads_bsse:.6f} eV)", f_out)

    gplot_path = write_curve_plot(args.output, rows)
    print_dual(f"\n{color_text('[Saved]', 'cyan')} Curve data -> {args.output}, {gplot_path} "
                f"(cd {os.path.dirname(args.output) or '.'} && gnuplot {os.path.basename(gplot_path)})",
                f_out)

    # Height-sweep approach curves: one per "family" (same site + same
    # adsorbate, swept across height) that has 2+ heights among the rows
    # actually analyzed above.
    families = {}
    for row in rows:
        family = _HEIGHT_SUFFIX_RE.sub('', row.label)
        families.setdefault(family, []).append(row)
    out_dir = os.path.dirname(args.output) or "."
    for family, family_rows in families.items():
        if len({r.height for r in family_rows}) < 2:
            continue
        dat_path, height_gplot = write_height_curve_plot(family, family_rows, out_dir)
        print_dual(f"{color_text('[Saved]', 'cyan')} Height curve -> {dat_path}, {height_gplot}", f_out)

    ranking_plot_path = os.path.join(out_dir, RANKING_PLOT_FILE)
    with capture_library_noise(library_warnings, "matplotlib plot_adsorption_ranking"):
        plot_adsorption_ranking(rows, ranking_plot_path, show=args.view_plots)
    print_dual(f"{color_text('[Saved]', 'cyan')} Ranking plot -> {ranking_plot_path}", f_out)

    if report_path:
        print_dual(f"{color_text('[Saved]', 'cyan')} Report     -> {report_path}", f_out)

    # --- [4] APPLY: always present (fixed section numbering), even when
    # --apply wasn't passed -- same "Not requested" pattern as adsorb.py's
    # own [3] ML PRE-SCREENING.
    print_section('[4] APPLY', f_out)
    if not args.apply:
        print_dual("  Not requested (pass --apply <path> to copy the most stable site's "
                    "structure.fdf -- BSSE-corrected ranking when available for every site, "
                    "else the uncorrected ranking -- to that path).", f_out)
    else:
        src = os.path.join(sites_root, apply_source_label, "structure.fdf")
        try:
            shutil.copy(src, args.apply)
        except OSError as e:
            print_dual(color_text(f"[ERROR] Could not copy '{src}' to '{args.apply}': {e}", 'red'), f_out)
        else:
            print_dual(f"{color_text('[Applied]', 'green')} {apply_source_label} -> {args.apply}", f_out)

    # --- [5] SUGGESTED NEXT ANALYSES: always present, references the same
    # apply_source_label winning site [3]/[4] already settled on. Every
    # suggestion is conditional ("if you also generated <file>...") since
    # none of these SIESTA outputs are guaranteed to exist -- each needs an
    # extra flag/block in the user's own --calc template that stb-adsorb
    # never forces (unlike Slab.DipoleCorrection/Spin polarized). ---
    print_section('[5] SUGGESTED NEXT ANALYSES', f_out)
    winning_site_dir = os.path.join(sites_root, apply_source_label)
    print_dual(f"  All paths below are relative to the winning site's own folder: "
                f"{winning_site_dir}", f_out)
    print_dual(color_text(
        "\n  - Bader charges/spin (if you saved '<label>.RHO'):\n"
        "      stb-bader --label <label>\n"
        "  - Density of states (if you enabled PDOS output, '<label>.PDOS.xml'):\n"
        "      stb-dos <label>.PDOS.xml --type total species\n"
        "  - Work function (if you saved '<label>.VT' -- SaveTotalPotential/"
        "SaveElectrostaticPotential):\n"
        "      stb-workfunction -l <label>\n"
        "      (also a direct dipole-correction consistency check, since [2b] above already "
        "reports the dipole this site's Slab.DipoleCorrection is meant to fix)\n"
        "  - Bonding/antibonding character (if you saved a full-BZ '<label>.WFSX' -- "
        "WriteWaveFunctions T + an explicit %block WaveFuncKPoints -- plus '<label>.HSX'):\n"
        "      stb-coop --label <label> --quantity coop --pair <adsorbate-atom> "
        "<surface-atom> --erange <E1> <E2> --sigma <meV>", 'cyan'), f_out)

    # --- [6] GIBBS FREE ENERGY (DG) PREP: opt-in, --compute-gibbs. ---
    print_section('[6] GIBBS FREE ENERGY (DG) PREP', f_out)
    if not args.compute_gibbs:
        print_dual("  Not requested (pass --compute-gibbs to write the vibrational-Hessian "
                    "displacement folders stb-adsorbGibbs, Stage 4, needs for a ZPE/entropy-"
                    "corrected Gibbs free energy of adsorption -- computed for the winning site "
                    "above).", f_out)
    else:
        winning_row = next((r for r in rows if r.label == apply_source_label), None)
        ads_name = winning_row.ads_name if winning_row else None
        ads_dir = adsorbate_dirs.get(ads_name)
        if ads_dir is None:
            print_dual(color_text(
                f"[ERROR] No isolated-adsorbate reference folder found for '{ads_name}' -- "
                "cannot build the reference Hessian.", 'red'), f_out)
        else:
            try:
                site_relaxed, site_used_relaxed = structure_io.read_relaxed_or_input(winning_site_dir)
            except (FileNotFoundError, ValueError) as e:
                print_dual(color_text(f"[ERROR] {e}", 'red'), f_out)
                site_relaxed, site_used_relaxed = None, False
            try:
                ads_relaxed, ads_used_relaxed = structure_io.read_relaxed_or_input(ads_dir)
            except (FileNotFoundError, ValueError) as e:
                print_dual(color_text(f"[ERROR] {e}", 'red'), f_out)
                ads_relaxed, ads_used_relaxed = None, False

            if not site_used_relaxed:
                print_dual(color_text(
                    f"[ERROR] '{winning_site_dir}' has no finished siesta.XV yet -- a Hessian "
                    "needs the actual relaxed geometry, not the pre-relaxation guess. Run SIESTA "
                    "there first, then re-run with --compute-gibbs.", 'red'), f_out)
            elif not ads_used_relaxed:
                print_dual(color_text(
                    f"[ERROR] '{ads_dir}' has no finished siesta.XV yet -- the isolated-adsorbate "
                    "reference must also be relaxed before its own Hessian means anything.",
                    'red'), f_out)
            else:
                gibbs_root = os.path.join(args.dir, "gibbs")
                gibbs_site_dir = os.path.join(gibbs_root, apply_source_label)
                gibbs_ads_dir = os.path.join(gibbs_root, f"{ads_name or 'default'}_isolated")

                site_calc_text = force_system_label(
                    read_original_calc_text(os.path.join(winning_site_dir, "calc.fdf")), "gibbs_site")
                spin_site, dipole_site, vdw_site = read_site_theory_flags(winning_site_dir)
                site_config_extra = GIBBS_SINGLE_POINT_BLOCK
                if spin_site:
                    site_config_extra += SPIN_POLARIZED_BLOCK
                if dipole_site:
                    site_config_extra += DIPOLE_CORRECTION_BLOCK
                if vdw_site:
                    site_config_extra += VDW_CORRECTION_BLOCK

                # read_original_calc_text (not a plain open().read()) strips the
                # site's own '%include config_extra.fdf' prefix -- this stage writes
                # its OWN config_extra.fdf (GIBBS_SINGLE_POINT_BLOCK + whatever
                # Spin/dipole/vdW flags read_site_theory_flags found) into every
                # disp_NNN/ folder below instead, same mechanism stb-adsorb/
                # stb-adsorbBsse already use for their own folders -- no dangling
                # include, since the file is actually written this time. DFTD3 is
                # read back via read_site_theory_flags(ads_dir) (the isolated
                # reference DOES get its own config_extra.fdf for this, from
                # adsorb.py's force_vdw=True) -- Slab.DipoleCorrection is
                # deliberately never re-applied here even if somehow present,
                # matching write_reference_folder's own adsorbate_dir call (no
                # force_dipole): a boxed molecule isn't a slab. Spin is forced
                # UNCONDITIONALLY instead of gated on read_site_theory_flags:
                # adsorb.py's force_spin_polarized() bakes 'Spin polarized'
                # directly into the isolated reference's calc.fdf body
                # (unconditionally, regardless of --force-spin), a mechanism
                # read_site_theory_flags can't see at all (it only scans a
                # folder's own config_extra.fdf, which the isolated reference
                # never gets one of for Spin) -- so relying on that flag here
                # would silently never fire. Adding it to config_extra.fdf too is
                # redundant with the body text (SIESTA's first-occurrence-wins
                # duplicate-tag rule makes this harmless) but removes the
                # single-point-of-failure on that separate mechanism staying in
                # sync, matching the explicit config_extra.fdf-based precaution
                # already used for every other Spin/dipole/vdW flag in this
                # workflow.
                _spin_ads_flag, _dipole_ads, vdw_ads = read_site_theory_flags(ads_dir)
                ads_calc_text = force_system_label(
                    read_original_calc_text(os.path.join(ads_dir, "calc.fdf")), "gibbs_isolated")
                ads_config_extra = GIBBS_SINGLE_POINT_BLOCK + SPIN_POLARIZED_BLOCK
                if vdw_ads:
                    ads_config_extra += VDW_CORRECTION_BLOCK

                n_total_site = len(site_relaxed.atoms)
                site_local_indices = list(range(n_substrate, n_total_site))
                ads_local_indices = list(range(len(ads_relaxed.atoms)))

                print_dual(f"  Winning site  : {apply_source_label}  ({n_total_site - n_substrate} "
                            "adsorbate atom(s))", f_out)
                print_dual(f"  Isolated ref  : {ads_dir}  ({len(ads_local_indices)} atom(s))", f_out)
                print_dual(f"  ZPE mode (site): {args.zpe_mode}  (isolated reference always uses "
                            "'local' -- see --zpe-mode --help)", f_out)

                write_local_hessian_folders(gibbs_ads_dir, ads_relaxed, ads_local_indices,
                                             args.displacement, ads_calc_text, ads_dir, "gibbs_isolated",
                                             config_extra_content=ads_config_extra)
                print_dual(f"  {color_text('[OK]', 'green')} {gibbs_ads_dir}/disp_001.."
                            f"disp_{6 * len(ads_local_indices):03d} (isolated-reference Hessian, "
                            "local mode)", f_out)

                if args.zpe_mode == "local":
                    write_local_hessian_folders(gibbs_site_dir, site_relaxed, site_local_indices,
                                                 args.displacement, site_calc_text, winning_site_dir,
                                                 "gibbs_site", config_extra_content=site_config_extra)
                    print_dual(f"  {color_text('[OK]', 'green')} {gibbs_site_dir}/disp_001.."
                                f"disp_{6 * len(site_local_indices):03d} (site Hessian, adsorbate "
                                "atoms only, substrate frozen)", f_out)
                else:  # full
                    print_dual(color_text(
                        "  [NOTE] 'full' mode needs a full phonon calculation of BOTH the winning "
                        "site AND the clean slab (stb-adsorbGibbs subtracts the clean slab's own "
                        "phonon ZPE/entropy, same reasoning as stb-herRefs/stb-oerRefs) -- roughly "
                        "double the cost of a single full phonon run.", 'yellow'), f_out)
                    from phonopy.interface.siesta import read_siesta
                    supercell_matrix = [[args.gibbs_supercell[0], 0, 0],
                                         [0, args.gibbs_supercell[1], 0],
                                         [0, 0, args.gibbs_supercell[2]]]

                    site_fdf_path = os.path.join(gibbs_site_dir, "_reference.fdf")
                    os.makedirs(gibbs_site_dir, exist_ok=True)
                    structure_io.write_fdf(site_relaxed, site_fdf_path)
                    site_unitcell = read_siesta(site_fdf_path)
                    site_phonon, site_supercells = build_phonon_displacements(
                        site_unitcell, supercell_matrix, args.displacement)
                    site_folders, _site_yaml = write_displacement_folders(
                        gibbs_site_dir, site_phonon, site_supercells, "structure.fdf",
                        os.path.join(winning_site_dir, "calc.fdf"), [])
                    site_symbols = sorted({sym for sym, _ in site_relaxed.atoms})
                    for d in site_folders:
                        for sym in site_symbols:
                            copy_pseudo(winning_site_dir, sym, d)
                        with open(os.path.join(d, CONFIG_EXTRA_FILE), "w") as f:
                            f.write(site_config_extra)
                        with open(os.path.join(d, "calc.fdf"), "w") as f:
                            f.write(structure_io.prepend_include(site_calc_text, CONFIG_EXTRA_FILE))
                    print_dual(f"  {color_text('[OK]', 'green')} {len(site_folders)} displacement "
                                f"folder(s) under {gibbs_site_dir}/", f_out)

                    clean_slab_dir = os.path.join(args.dir, "clean_slab")
                    gibbs_clean_dir = os.path.join(gibbs_root, "clean_slab_full")
                    clean_fdf_path = os.path.join(gibbs_clean_dir, "_reference.fdf")
                    os.makedirs(gibbs_clean_dir, exist_ok=True)
                    clean_template = structure_io.read_fdf(clean_slab_fdf_path)
                    structure_io.write_fdf(clean_template, clean_fdf_path)
                    clean_unitcell = read_siesta(clean_fdf_path)
                    clean_phonon, clean_supercells = build_phonon_displacements(
                        clean_unitcell, supercell_matrix, args.displacement)
                    # clean_slab/ normally gets Slab.DipoleCorrection + DFTD3 forced
                    # unconditionally by stb-adsorb (never Spin -- that's the winning
                    # SITE's own choice, unrelated to the bare clean slab) -- read its
                    # own config_extra.fdf rather than assuming, same as the site/
                    # isolated-reference sides above. Real gap found and fixed while
                    # doing this refactor: this branch used to never call
                    # read_site_theory_flags(clean_slab_dir) at all, so the clean-slab
                    # phonon reference silently missed DipoleCorrection/DFTD3 that the
                    # actual clean_slab/ folder carries.
                    spin_clean, dipole_clean, vdw_clean = read_site_theory_flags(clean_slab_dir)
                    clean_calc_text = force_system_label(
                        read_original_calc_text(os.path.join(clean_slab_dir, "calc.fdf")), "gibbs_clean")
                    clean_config_extra = GIBBS_SINGLE_POINT_BLOCK
                    if spin_clean:
                        clean_config_extra += SPIN_POLARIZED_BLOCK
                    if dipole_clean:
                        clean_config_extra += DIPOLE_CORRECTION_BLOCK
                    if vdw_clean:
                        clean_config_extra += VDW_CORRECTION_BLOCK
                    clean_folders, _clean_yaml = write_displacement_folders(
                        gibbs_clean_dir, clean_phonon, clean_supercells, "structure.fdf",
                        os.path.join(clean_slab_dir, "calc.fdf"), [])
                    clean_symbols = sorted({sym for sym, _ in clean_template.atoms})
                    for d in clean_folders:
                        for sym in clean_symbols:
                            copy_pseudo(clean_slab_dir, sym, d)
                        with open(os.path.join(d, CONFIG_EXTRA_FILE), "w") as f:
                            f.write(clean_config_extra)
                        with open(os.path.join(d, "calc.fdf"), "w") as f:
                            f.write(structure_io.prepend_include(clean_calc_text, CONFIG_EXTRA_FILE))
                    print_dual(f"  {color_text('[OK]', 'green')} {len(clean_folders)} displacement "
                                f"folder(s) under {gibbs_clean_dir}/", f_out)

                print_dual(f"\n  Run SIESTA in every folder written above, then use: "
                            f"stb-adsorbGibbs --dir {args.dir}", f_out)

    # --- [7] LIBRARY WARNINGS: always last. ---
    print_section('[7] LIBRARY WARNINGS', f_out)
    if library_warnings:
        print_dual(color_text(
            "Messages emitted by external libraries (sisl, pymatgen, matplotlib) during this "
            "run -- collected here instead of interleaved with the report above; harmless in "
            "almost every case, but worth a look if a section above looks suspicious.", 'cyan'), f_out)
        for entry in library_warnings:
            print_dual(entry, f_out)
    else:
        print_dual("No library warnings.", f_out)

    if f_out:
        f_out.close()

    if args.view:
        print(f"\n{color_text('--view:', 'cyan')} opening {len(view_frames)} frame(s) in ASE's "
              "interactive viewer (use the frame slider/menu to step through them):")
        for i, (label, _atoms) in enumerate(view_frames):
            print(f"  {i} = {label}")
        if view_frames:
            view_structure_interactive([atoms for _label, atoms in view_frames])
        else:
            print(color_text("  No readable structures found -- nothing to view.", 'yellow'))


if __name__ == "__main__":
    main()
