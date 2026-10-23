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
from stb.core import structure_io, kspace, fdf_directives, siesta_log
from stb.core.optical_setup import (reciprocal_spacings, spacing_band, symmetry_analysis,
                                    _SPACING_COARSE, _SPACING_FINE)
from stb.core.cli import color_text, show_intro, print_dual, print_section, print_table
from stb.core.calc_directives import build_optical_block, OPTICAL_DIRECTION_VECTORS, OPTICAL_OFFDIAG_PAIRS
from stb.core.adsorption_sites import SINGLE_POINT_BLOCK
from stb.core.pseudopotentials import resolve_pseudo_source, copy_pseudo, pseudo_valence

REPORT_FILE = "optical_stage1.txt"
DEFAULT_KGRID_DENSITY = 0.2     # 1/Ang, same convention as stb-kgrid and the other workflows
CONFIG_EXTRA_FILE = "config_extra.fdf"

# See core/calc_directives.py's own docstring for the full diagonal-vs-
# -biaxial derivation -- shared with stb-opticalAnalysis (Stage 2), which
# performs the actual eps_ij reconstruction from the biaxial folders this
# Stage 1 writes.
_AXIS_VECTORS = OPTICAL_DIRECTION_VECTORS
_OFFDIAG_REQUIRES = OPTICAL_OFFDIAG_PAIRS


def detect_structure_vacuum(fdf_structure, vacuum_gap):
    """Detects (vacuum_axes, is_2d) from `fdf_structure`'s own lattice/
    coordinates -- same pattern as gqca.py's own detect_structure_vacuum,
    duplicated here (not extracted to core/) since this is currently the
    second, not third, consumer of this exact 4-line combination and the
    two callers' surrounding context differs enough (gqca.py re-detects
    per representative-cluster structure; here it's a single detection
    on the one input structure, reused informationally in the report and
    NOT used to change which/how many direction folders get written --
    see stb-optical --help for why no direction is skipped for 2D
    inputs).
    """
    positions = np.array([pos for _, pos in fdf_structure.atoms])
    is_cartesian = fdf_structure.coord_format == 'cartesian'
    frac_coords = kspace.to_fractional(positions, fdf_structure.lattice, is_cartesian)
    vacuum_axes = kspace.detect_vacuum_axes(frac_coords, fdf_structure.lattice, vacuum_gap)
    return vacuum_axes, any(vacuum_axes)


# --- Report helpers ---------------------------------------------------------------

_MD_EXTRA_BLOCK = (
    "MD.NumCGsteps         0\n"
    "MD.VariableCell       false\n"
)


def write_direction_folder(out_dir, fdf_structure, calc_text, extra_fdf_text, pp_path):
    """Writes structure.fdf (unchanged -- no new geometry, this workflow
    never touches atomic positions/lattice) + calc.fdf (the user's own
    --calc template, untouched except for the %include at its very top --
    see structure_io.prepend_include) + config_extra.fdf (the forced
    single-point directives plus this one direction's own Optical.* block,
    `extra_fdf_text`, varying per folder) + copied pseudos for one
    direction folder. Same config_extra.fdf convention as every other
    WORKFLOW_TOOLS prep stage in this suite (see CLAUDE.md's
    "config_extra.fdf" convention note) -- never edits the user's own
    template text directly, and SIESTA's fdf reader is first-occurrence
    -wins for a duplicate label, so %include'ing this ahead of the
    template correctly overrides any conflicting MD.*/Optical.* setting
    the template itself might already carry. The ground-state SCF k-grid
    is forced into config_extra.fdf only when it is REQUIRED -- the
    template defines none (SIESTA would then sample Gamma only) -- or when
    the user asks for one (--kgrid / --kgrid-density); a template's own
    grid is otherwise left alone (see main()).
    """
    os.makedirs(out_dir, exist_ok=True)
    structure_io.write_fdf(fdf_structure, os.path.join(out_dir, "structure.fdf"))
    with open(os.path.join(out_dir, "calc.fdf"), "w") as f:
        f.write(calc_text)
    with open(os.path.join(out_dir, CONFIG_EXTRA_FILE), "w") as f:
        f.write(extra_fdf_text)
    symbols = sorted({sym for sym, _ in fdf_structure.atoms})
    for sym in symbols:
        copy_pseudo(pp_path, sym, out_dir)


def main():
    parser = argparse.ArgumentParser(
        description=f"""{color_text("Stage 1 of 2: writes one SIESTA folder per requested "
        "Cartesian/biaxial polarization direction, ready for a linear-response Optical.Calculate run.", 'bold')}
Takes an ALREADY-RELAXED structure (--file) and a calc.fdf template (--calc, default: calc.fdf in the
current directory), and for each requested direction (--directions, default xx yy zz) writes
'<output-dir>/dir_<axis>/' with the SAME structure (no new geometry -- unlike this suite's other
WORKFLOW_TOOLS prep stages, there's nothing to deform/search/enumerate here). Every forced directive
-- MD.TypeOfRun CG / MD.Steps 0 / MD.NumCGsteps 0 / MD.VariableCell false (single-point, no
relaxation), that folder's own Optical.Mesh/Optical.Vector/Optical.Broaden/Optical.NumberOfBands block
(SIESTA's own linear-response/RPA dielectric-function machinery -- see the SIESTA manual's Optical.*
directives) and, when it is REQUIRED, the ground-state SCF k-grid -- is written into a config_extra.fdf
sidecar, %include'd at the very top of calc.fdf, never edited into your own --calc template text
directly (this suite's standard config_extra.fdf convention -- see CLAUDE.md). The SCF k-grid is
forced when your template defines none (SIESTA would sample Gamma only; a grid is suggested at
density 0.2 1/Ang, vacuum axes getting 1 division) or when you ask for one (--kgrid X Y Z, or
--kgrid-density D); a template that already has a grid keeps it. Doesn't run SIESTA -- run each folder yourself, then use stb-opticalAnalysis (which
aggregates every direction folder present into one combined report/spectrum, so run it once against
the whole --output-dir root, not once per folder).

[DIAGONAL vs. BIAXIAL DIRECTIONS] 'xx'/'yy'/'zz' set Optical.Vector to the pure Cartesian unit axis
(e_x/e_y/e_z) -- the diagonal dielectric-tensor component eps_ii comes directly from that one folder.
'xy'/'xz'/'yz' set Optical.Vector to the NORMALIZED BISECTOR of the two axes involved (e.g. xy ->
(e_x+e_y)/sqrt(2)) -- SIESTA's Optical.Vector always computes a single-direction response n^T.eps.n,
never a genuine off-diagonal component directly, so the bisector is the standard trick to reach one:
for a symmetric dielectric tensor, response along that bisector expands to
(eps_ii + eps_jj)/2 + eps_ij, so once eps_ii, eps_jj AND the bisector's own response are all known,
  eps_ij = eps_(bisector ij) - (eps_ii + eps_jj) / 2
recovers the true off-diagonal component. Requesting 'xy' without also requesting 'xx'/'yy' in the
same (or an earlier) run still writes a valid folder -- it's simply not reconstructible into eps_xy
until the 2 diagonal folders it needs also exist; the Stage 1 report flags this per direction.

Stage 2 (stb-opticalAnalysis) auto-detects each folder's direction from its own Optical.Vector value
(not the folder name) and, whenever a biaxial direction's 2 matching diagonal folders are ALSO
present in the same run, applies the reconstruction formula above automatically -- eps_ij(E) is
written to its own dedicated output alongside the raw per-direction eps_ii(E)/eps_jj(E) results, no
extra flag needed. A biaxial folder without its matching diagonal pair is still read normally by
Stage 2 (its raw eps_(bisector)(E) response is real data), just not reconstructible into eps_ij yet.

[KNOWN LIMITATION, 2D/slab inputs] No direction is skipped or specially treated for a vacuum-padded
input -- every requested direction (including out-of-plane and any biaxial one involving it) is still
written and computed. For a slab, the periodic-supercell dielectric response along ANY direction is
diluted by the vacuum region (a well-known methodological caveat for slab dielectric/polarizability
DFT calculations) -- no vacuum-thickness rescaling is applied here. Interpret out-of-plane 2D results
with this in mind.""",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s -f relaxed.fdf -p dojo\n"
               "  %(prog)s -f relaxed.fdf -p dojo --directions xx yy zz xy xz yz\n"
    )

    parser.add_argument("-f", "--file", type=str, required=True,
                         help="Already-relaxed input structure.fdf.")
    parser.add_argument("-c", "--calc", type=str, default="calc.fdf",
                         help="calc.fdf template for the Optical calculations (default: calc.fdf).")
    parser.add_argument("-p", "--pseudo-dir", type=str, default="",
                         help="Pseudopotentials source: a bundled bank or a folder path.")
    parser.add_argument("--directions", type=str, nargs='+', default=["xx", "yy", "zz"],
                         choices=["xx", "yy", "zz", "xy", "xz", "yz"],
                         help="Direction(s) to compute (default: xx yy zz). 'xx'/'yy'/'zz' are the "
                              "diagonal Cartesian components (Optical.Vector = a unit axis, "
                              "eps_ii read directly). 'xy'/'xz'/'yz' are BIAXIAL directions "
                              "(Optical.Vector = the normalized bisector of the 2 axes, e.g. "
                              "(1,1,0)/sqrt(2) for xy) -- the off-diagonal eps_ij component is only "
                              "reconstructible from these once the matching xx/yy/zz pair also "
                              "exists; see the description above for the exact formula.")
    parser.add_argument("--optical-mesh", type=int, nargs=3, default=[10, 10, 10], metavar=("NX", "NY", "NZ"),
                         help="Optical.Mesh k-point sampling for the optical-transition integral "
                              "(default: 10 10 10).")
    parser.add_argument("--optical-broaden", type=float, default=0.2, metavar="EV",
                         help="Optical.Broaden Gaussian broadening in eV (default: 0.2).")
    parser.add_argument("--optical-nbands", type=int, default=None, metavar="N",
                         help="Optical.NumberOfBands (default: SIESTA's own default -- all bands).")
    parser.add_argument("--kgrid", type=int, nargs=3, default=None, metavar=("NX", "NY", "NZ"),
                         help="Ground-state SCF k-grid forced into every folder's config_extra.fdf "
                              "(e.g. --kgrid 8 8 8). Overrides the template's own grid. Without --kgrid "
                              "or --kgrid-density the template's grid is kept; if the template has none, "
                              "a grid is suggested automatically (density 0.2 1/Ang).")
    parser.add_argument("--kgrid-density", type=float, default=None, metavar="D",
                         help="Suggest the SCF k-grid from a target density in 1/Ang (same convention as "
                              "stb-kgrid; vacuum axes get 1 division) and force it into config_extra.fdf, "
                              "overriding the template's grid. Ignored when --kgrid is given.")
    parser.add_argument("--vacuum-gap", type=float, default=10.0,
                         help="Vacuum-axis detection threshold in Ang (vacuum axes get a single SCF "
                              "k-point division and drive the dimensionality notes; default: 10.0).")
    parser.add_argument("-O", "--output-dir", type=str, default="optical_study",
                         help="Root directory for the whole optical-properties workflow "
                              "(default: optical_study).")
    parser.add_argument("--save-report", action="store_true",
                         help=f"Also persist the report to <output-dir>/{REPORT_FILE}. Off by default.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-optical {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("OPTICAL PROPERTIES WORKFLOW -- STAGE 1: DIRECTION FOLDERS", 'bold'))
    print("-" * 60)

    if not os.path.exists(args.file):
        print(color_text(f"[ERROR] Structure file '{args.file}' not found.", 'red'))
        sys.exit(1)
    if not os.path.exists(args.calc):
        print(color_text(f"[ERROR] Calc file '{args.calc}' not found. The calc.fdf template defaults to "
                          "./calc.fdf -- pass -c/--calc to point at another file.", 'red'))
        sys.exit(1)
    if args.pseudo_dir:
        try:
            args.pseudo_dir = resolve_pseudo_source(args.pseudo_dir)
        except ValueError as e:
            print(color_text(f"[ERROR] {e}", 'red'))
            sys.exit(1)

    structure = structure_io.read_fdf(args.file)
    vacuum_axes, is_2d = detect_structure_vacuum(structure, args.vacuum_gap)

    with open(args.calc) as f:
        calc_template = f.read()
    template_values = fdf_directives.parse_directives(calc_template)
    template_kgrid = fdf_directives.parse_kgrid(calc_template)

    output_root = args.output_dir
    os.makedirs(output_root, exist_ok=True)
    report_path = os.path.join(output_root, REPORT_FILE) if args.save_report else None

    f_out = open(report_path, "w") if report_path else None

    print_dual(f"{color_text('===== OPTICAL STAGE 1 REPORT (DIRECTION FOLDERS) =====', 'magenta')}", f_out)

    species = structure_io.species_list(structure)
    counts = structure_io.atom_counts(structure)
    formula = " ".join(f"{sym}{counts[sym]}" for sym in species)
    lattice = structure_io.lattice_only(structure)
    cell_volume = float(abs(np.linalg.det(lattice)))
    n_atoms = len(structure.atoms)
    lengths = [float(np.linalg.norm(v)) for v in lattice]
    angles = [float(np.degrees(np.arccos(np.dot(lattice[(i + 1) % 3], lattice[(i + 2) % 3])
                                          / (lengths[(i + 1) % 3] * lengths[(i + 2) % 3])))) for i in range(3)]
    warnings = []                       # collected, then repeated together in [8]

    def warn(text):
        warnings.append(text)
        return color_text(f"[WARNING] {text}", 'yellow')

    # ------------------------------------------------------------------ [0]
    print_section('[0] RUN METADATA', f_out)
    print_dual(f"Tool            : stb-optical {VERSION}", f_out)
    print_dual(f"Date            : {datetime.datetime.now():%Y-%m-%d %H:%M}", f_out)
    print_dual(f"Structure       : {args.file}", f_out)
    print_dual(f"Formula         : {formula}  ({n_atoms} atom(s), {len(species)} species: {', '.join(species)})", f_out)
    print_dual(f"Cell volume     : {cell_volume:.4f} Ang^3", f_out)
    print_dual(f"Calc template   : {args.calc}" + ("  (default)" if args.calc == "calc.fdf" else ""), f_out)
    print_dual(f"Directions      : {' '.join(args.directions)}", f_out)
    print_dual(f"Pseudo source   : {args.pseudo_dir if args.pseudo_dir else '(none -- pseudos not copied)'}", f_out)
    print_dual(f"Output root     : {output_root}", f_out)

    # ------------------------------------------------------------------ [1]
    print_section('[1] INPUT STRUCTURE', f_out)
    print_dual(f"Lattice (Ang)   : a={lengths[0]:.4f}  b={lengths[1]:.4f}  c={lengths[2]:.4f}", f_out)
    print_dual(f"Angles (deg)    : alpha={angles[0]:.2f}  beta={angles[1]:.2f}  gamma={angles[2]:.2f}", f_out)
    print_dual(f"Atoms per Ang^3 : {n_atoms / cell_volume:.4f}", f_out)
    rows, total_electrons = [], 0.0
    for sym in species:
        zval, source = pseudo_valence(sym, args.pseudo_dir)
        if zval is not None:
            total_electrons += zval * counts[sym]
        rows.append(([sym, str(structure.species_meta.get(sym, {}).get('Z', '?')), str(counts[sym]),
                       f"{zval:g}" if zval is not None else "?",
                       f"{zval * counts[sym]:g}" if zval is not None else "?", source], None))
    print_table(["Species", "Z", "Atoms", "Valence e-", "Total e-", "Valence from"], rows, f_out)
    n_occupied = int(round(total_electrons / 2.0))
    print_dual(f"Valence electrons : {total_electrons:g}  ->  {n_occupied} occupied band(s) in a "
                "spin-unpolarized calculation", f_out)
    print_dual(f"Detected : {kspace.dimensionality_label(vacuum_axes)}", f_out)
    n_vacuum = sum(vacuum_axes)
    if n_vacuum == 1:
        print_dual(color_text(
            "[KNOWN LIMITATION] Vacuum-padded (2D/slab) input -- the periodic-supercell "
            "dielectric response along ANY direction (including out-of-plane) is diluted by "
            "the vacuum region, a known methodological caveat for slab dielectric/"
            "polarizability DFT calculations. No vacuum-thickness rescaling is applied. Every "
            "requested direction is still computed (no direction is skipped for 2D inputs). "
            "Stage 2 can restore the intrinsic response: stb-opticalAnalysis "
            "--dimensionality-correction (2D also needs --thickness).", 'yellow'), f_out)
    elif n_vacuum == 2:
        print_dual(color_text(
            "[KNOWN LIMITATION] Vacuum on 2 axes (1D wire/chain): the periodic-supercell response is "
            "diluted by the vacuum; Stage 2's --dimensionality-correction does not implement 1D, so only "
            "the raw supercell values are reported.", 'yellow'), f_out)
    elif n_vacuum == 3:
        print_dual(color_text(
            "[KNOWN LIMITATION] Vacuum on all 3 axes (isolated molecule/cluster): the periodic-supercell "
            "response is the dilute-gas limit and scales with the cell volume. Stage 2's "
            "--dimensionality-correction extracts the molecular polarizability alpha(E) from it.", 'yellow'), f_out)

    # SCF k-grid: forced into config_extra.fdf when it is REQUIRED (the template defines none --
    # SIESTA would sample Gamma only) or when asked for (--kgrid / --kgrid-density); a template's
    # own grid is otherwise kept as is.
    template_has_kgrid = bool(template_kgrid) or "kgridcutoff" in template_values
    forced_kgrid, kgrid_source = None, None
    if args.kgrid is not None:
        forced_kgrid, kgrid_source = tuple(args.kgrid), "explicit --kgrid"
    elif args.kgrid_density is not None or not template_has_kgrid:
        density = args.kgrid_density if args.kgrid_density is not None else DEFAULT_KGRID_DENSITY
        try:
            forced_kgrid = tuple(kspace.compute_monkhorts(lattice[0], lattice[1], lattice[2], density,
                                                          vacuum_axes=vacuum_axes))
        except ValueError as e:
            print(color_text(f"[ERROR] {e}", 'red'))
            sys.exit(1)
        kgrid_source = (f"--kgrid-density {density:g} 1/Ang" if args.kgrid_density is not None else
                        f"auto-suggested, density {density:g} 1/Ang (the template defines no k-grid)")

    # ------------------------------------------------------------------ [2]
    print_section('[2] GROUND STATE READ FROM YOUR calc.fdf (the optical response sits on top of it)', f_out)
    theory_keys = [("xcfunctional", "XC functional"), ("xcauthors", "XC authors"),
                   ("paobasistype", "Basis type"), ("paobasissize", "Basis size"),
                   ("paoenergyshift", "PAO energy shift"), ("meshcutoff", "Mesh cutoff"),
                   ("spin", "Spin"), ("spinpolarized", "SpinPolarized"),
                   ("electronictemperature", "Electronic temperature"), ("dmtolerance", "DM tolerance"),
                   ("maxscfiterations", "Max SCF iterations")]
    theory_rows = [([label, template_values[key]], None) for key, label in theory_keys if key in template_values]
    if template_kgrid:
        kg_spacing = reciprocal_spacings(lattice, template_kgrid)
        theory_rows.append((["SCF k-grid (kgrid.MonkhorstPack)",
                             "x".join(str(v) for v in template_kgrid)
                             + f"  (spacing {kg_spacing[0]:.3f} {kg_spacing[1]:.3f} {kg_spacing[2]:.3f} 1/Ang)"], None))
    elif "kgridcutoff" in template_values:
        theory_rows.append((["SCF k-grid (kgrid.Cutoff)", template_values["kgridcutoff"]], None))
    print_table(["Setting", "Value in the template"], theory_rows or [(["(none recognized)", "-"], None)], f_out)
    if forced_kgrid is not None:
        kg_sp = reciprocal_spacings(lattice, forced_kgrid)
        print_dual(f"SCF k-grid forced     : {'x'.join(str(v) for v in forced_kgrid)}  "
                    f"(spacing {kg_sp[0]:.3f} {kg_sp[1]:.3f} {kg_sp[2]:.3f} 1/Ang) -- {kgrid_source}; "
                    f"written into config_extra.fdf" + (f", overriding the template's "
                    f"{'x'.join(str(v) for v in template_kgrid)}" if template_kgrid else ""), f_out)
    else:
        print_dual("SCF k-grid forced     : no -- the template's own grid is kept "
                    "(use --kgrid / --kgrid-density to override it)", f_out)

    spin_value = (template_values.get("spin", "") + " " + template_values.get("spinpolarized", "")).lower()
    if any(word in spin_value for word in ("polarized", "true", ".true.")) and "non" not in spin_value:
        print_dual(color_text(
            "[NOTE] Spin-polarized template: stb-optical does not special-case spin; check the "
            "Stage 2 results against the SIESTA manual's treatment of Optical with spin.", 'cyan'), f_out)
    if "paobasissize" in template_values and template_values["paobasissize"].upper() in ("SZ", "DZ", "SZP"):
        print_dual(color_text(
            f"[NOTE] Basis {template_values['paobasissize']}: optical transitions involve the "
            "unoccupied (conduction) bands, which a minimal/small basis describes poorly; DZP or TZP "
            "with a small PAO.EnergyShift is the usual starting point.", 'cyan'), f_out)

    # ------------------------------------------------------------------ [3]
    print_section('[3] OPTICAL PARAMETERS (applied identically to every direction)', f_out)
    mesh = args.optical_mesh
    n_opt_k = int(np.prod(mesh))
    spacing = reciprocal_spacings(lattice, mesh)
    mesh_rows = []
    for axis_i, (name, n_i, sp_i) in enumerate(zip("abc", mesh, spacing)):
        vac = vacuum_axes[axis_i]
        mesh_rows.append(([f"axis {name}", str(n_i), f"{sp_i:.4f}", "vacuum axis" if vac else spacing_band(sp_i)], None))
    print_table(["Optical.Mesh", "Divisions", "Spacing (1/Ang)", "Assessment"], mesh_rows, f_out)
    print_dual(f"Optical k-points  : {n_opt_k} total (independent of the SCF k-grid; the optical "
                "integral is evaluated on this mesh)", f_out)
    coarse = [name for axis_i, (name, sp_i) in enumerate(zip("abc", spacing))
              if not vacuum_axes[axis_i] and spacing_band(sp_i) == "coarse"]
    if coarse:
        print_dual(warn(f"Optical.Mesh is coarse along {', '.join(coarse)} (spacing > {_SPACING_COARSE} 1/Ang): "
                         "eps2(E) will be spiky and peak positions/heights unreliable. Converge by "
                         "re-running with a denser --optical-mesh until the spectrum stops changing."), f_out)
    vac_dense = [name for axis_i, (name, n_i) in enumerate(zip("abc", mesh)) if vacuum_axes[axis_i] and n_i > 1]
    if vac_dense:
        print_dual(color_text(
            f"[NOTE] Optical.Mesh has more than 1 division along the vacuum axis ({', '.join(vac_dense)}): "
            "there is no dispersion along it, so those k-points only cost time; 1 is enough.", 'cyan'), f_out)

    print_dual("", f_out)
    optical_rows = [
        (["Optical.Broaden", f"{args.optical_broaden} eV",
          "Gaussian width applied to every interband transition before summing into eps2(E)"], None),
        (["Optical.NumberOfBands",
          str(args.optical_nbands) if args.optical_nbands else "(SIESTA default)",
          "bands entering the transition sum; must include conduction bands"], None),
        (["Optical.PolarizationType", "polarized",
          "makes SIESTA honour Optical.Vector (the default would silently average all directions)"], None),
    ]
    print_table(["Parameter", "Value", "Meaning"], optical_rows, f_out)
    if args.optical_nbands:
        n_cond = args.optical_nbands - n_occupied
        print_dual(f"Bands included    : {args.optical_nbands} = {n_occupied} occupied + {max(n_cond, 0)} unoccupied", f_out)
        if n_cond <= 0:
            print_dual(warn(f"Optical.NumberOfBands ({args.optical_nbands}) does not exceed the "
                             f"{n_occupied} occupied bands: no conduction states, no interband absorption. "
                             "Raise --optical-nbands (or leave it to SIESTA's default)."), f_out)
        elif n_cond < 4:
            print_dual(warn(f"Only {n_cond} unoccupied band(s) are included: the spectrum will be "
                             "truncated above the first few transitions. Raise --optical-nbands for a "
                             "wide energy range."), f_out)
    if args.optical_broaden < 0.05:
        print_dual(warn(f"Optical.Broaden = {args.optical_broaden} eV is very small: unless the optical mesh is "
                         "very dense the spectrum is a forest of single-transition spikes."), f_out)
    elif args.optical_broaden > 0.3:
        print_dual(color_text(
            f"[NOTE] Optical.Broaden = {args.optical_broaden} eV washes out fine structure (E1/E2-like "
            "peaks merge); 0.05-0.2 eV is typical for a converged mesh.", 'cyan'), f_out)

    # ------------------------------------------------------------------ [4]
    print_section('[4] DIRECTIVES FORCED THROUGH config_extra.fdf (your template is not edited)', f_out)
    print_dual(f"Every direction folder gets its own {CONFIG_EXTRA_FILE}, and calc.fdf starts with "
                f"'%include {CONFIG_EXTRA_FILE}' followed by your template untouched. SIESTA's fdf reader "
                "is first-occurrence-wins for a duplicate label, so the forced values win over anything "
                "the template says.", f_out)

    def template_has(*keys):
        return next((f"{template_values[k]}" for k in keys if k in template_values), "(absent)")

    forced_rows = [
        (["kgrid.MonkhorstPack",
          "x".join(str(v) for v in forced_kgrid) if forced_kgrid else "(not forced)",
          ("x".join(str(v) for v in template_kgrid) if template_kgrid else
           template_has("kgridcutoff")),
          kgrid_source if forced_kgrid else "the template's own SCF grid is used"], None),
        (["MD.TypeOfRun", "CG", template_has("mdtypeofrun"), "single-point run"], None),
        (["MD.Steps", "0", template_has("mdsteps"), "no relaxation step"], None),
        (["MD.NumCGsteps", "0", template_has("mdnumcgsteps"), "belt and braces: no CG step whatever the template says"], None),
        (["MD.VariableCell", "false", template_has("mdvariablecell"), "the input cell is evaluated as given"], None),
        (["OpticalCalculation", "T", template_has("opticalcalculation"), "switches the optical module on"], None),
        (["Optical.PolarizationType", "polarized", template_has("opticalpolarizationtype"),
          "required for Optical.Vector to be used"], None),
        (["Optical.Broaden", f"{args.optical_broaden} eV", template_has("opticalbroaden"), "from --optical-broaden"], None),
        (["Optical.Mesh", " ".join(str(n) for n in mesh), "block set" if "opticalmesh" in template_values else "(absent)",
          "from --optical-mesh"], None),
        (["Optical.Vector", "per direction, see [5]", "(absent)", "one polarization direction per folder"], None),
        (["Optical.NumberOfBands", str(args.optical_nbands) if args.optical_nbands else "(not forced)",
          template_has("opticalnumberofbands"), "from --optical-nbands"], None),
    ]
    print_table(["Directive", "Forced value", "Template had", "Effect"], forced_rows, f_out)
    print_dual(color_text(
        "The rest of the ground-state SCF (structure, basis, mesh cutoff, ...) is whatever your template "
        "says: every folder repeats the same SCF and only the polarization differs.", 'cyan'), f_out)

    # ------------------------------------------------------------------ [5]
    print_section('[5] DIRECTIONS & DIELECTRIC-TENSOR MAPPING', f_out)
    print_dual(f"{'Direction':<10} {'Kind':<12} {'Optical.Vector':<24} {'Gives access to':<28} {'Needs also':<12}", f_out)
    print_dual("-" * 60, f_out)
    for axis in args.directions:
        axis_vec = _AXIS_VECTORS[axis]
        vec_str = f"{axis_vec[0]:.4f} {axis_vec[1]:.4f} {axis_vec[2]:.4f}"
        if axis in _OFFDIAG_REQUIRES:
            need_a, need_b = _OFFDIAG_REQUIRES[axis]
            have = (need_a in args.directions, need_b in args.directions)
            needs_str = f"{need_a}+{need_b}"
            if not all(have):
                missing = [n for n, ok in zip((need_a, need_b), have) if not ok]
                needs_str += color_text(f" (missing {'/'.join(missing)} in this run)", 'yellow')
            print_dual(f"{axis:<10} {'biaxial':<12} {vec_str:<24} "
                        f"{'eps_' + axis + ' (needs reconstruction)':<28} {needs_str}", f_out)
        else:
            print_dual(f"{axis:<10} {'diagonal':<12} {vec_str:<24} {'eps_' + axis + ' (direct)':<28} {'--':<12}", f_out)
    if any(axis in _OFFDIAG_REQUIRES for axis in args.directions):
        print_dual(color_text(
            "\nReconstruction formula for any biaxial direction above: "
            "eps_ij = eps_(bisector ij) - (eps_ii + eps_jj) / 2 -- stb-opticalAnalysis applies this "
            "automatically once all 3 folders involved (ii, jj, ij) are present and done; no extra "
            "flag needed.", 'cyan'), f_out)

    sym = symmetry_analysis(structure)
    print_dual("\nWhat crystal symmetry says about the tensor (in the Cartesian frame of your cell):", f_out)
    if sym is None:
        print_dual("  (symmetry analysis unavailable for this structure)", f_out)
    else:
        print_dual(f"  Space group           : {sym['symbol']} (No. {sym['number']}), {sym['system']} "
                    f"-- {sym['n_ops']} operations at symprec {sym['symprec']} Ang", f_out)
        print_dual(f"  Independent components: {sym['independent']} of 6", f_out)
        eq = "; ".join("=".join(f"eps_{c}" for c in cls) for cls in sym['classes'])
        print_dual(f"  Non-zero components   : {eq}", f_out)
        if sym['zero']:
            print_dual(f"  Zero by symmetry      : {', '.join('eps_' + c for c in sym['zero'])}", f_out)
        print_dual(f"  Minimal direction set : {' '.join(sym['minimal'])}", f_out)
        redundant = [a for a in args.directions if a not in sym['minimal']]
        absent = [c for c in sym['minimal'] if c not in args.directions]
        if absent:
            print_dual(warn(f"The requested directions do not include {' '.join(absent)}, which symmetry "
                             "does not give from the others."), f_out)
        elif redundant:
            print_dual(color_text(
                f"[NOTE] {' '.join(redundant)} {'follows' if len(redundant) == 1 else 'follow'} from symmetry, "
                "so " + ("it is" if len(redundant) == 1 else "they are") + " redundant; computing "
                + ("it" if len(redundant) == 1 else "them") + " is a consistency check (identical results "
                "expected), not new information. "
                f"--directions {' '.join(sym['minimal'])} would save {len(redundant)} of "
                f"{len(args.directions)} runs.", 'cyan'), f_out)

    # ------------------------------------------------------------------ [6]
    print_section('[6] FOLDERS WRITTEN', f_out)

    forced_calc_text = structure_io.prepend_include(calc_template, CONFIG_EXTRA_FILE)
    kgrid_block = ""
    if forced_kgrid is not None:
        kgrid_block = (f"# Ground-state SCF k-grid ({kgrid_source}) -- takes precedence over --calc's own grid.\n"
                       f"kgrid.MonkhorstPack   [{forced_kgrid[0]}  {forced_kgrid[1]}  {forced_kgrid[2]}]\n\n")
    folder_rows, pseudo_problems = [], []
    for axis in args.directions:
        axis_vec = _AXIS_VECTORS[axis]
        optical_block = build_optical_block(args.optical_mesh, args.optical_broaden, axis_vec,
                                             args.optical_nbands)
        extra_fdf_text = (kgrid_block + SINGLE_POINT_BLOCK + _MD_EXTRA_BLOCK + "\n" + optical_block)
        out_dir = os.path.join(output_root, f"dir_{axis}")
        write_direction_folder(out_dir, structure, forced_calc_text, extra_fdf_text, args.pseudo_dir)
        kind = "biaxial" if axis in _OFFDIAG_REQUIRES else "diagonal"
        files = ["structure.fdf", "calc.fdf", CONFIG_EXTRA_FILE]
        for sym_label in species:
            for ext in ("psml", "psf"):
                if os.path.exists(os.path.join(out_dir, f"{sym_label}.{ext}")):
                    files.append(f"{sym_label}.{ext}")
                    break
            else:
                if args.pseudo_dir:
                    pseudo_problems.append(f"{out_dir}: {sym_label}")
        folder_rows.append(([f"dir_{axis}/", axis, kind,
                              f"{axis_vec[0]:.4f} {axis_vec[1]:.4f} {axis_vec[2]:.4f}",
                              ", ".join(files)], None))
    print_dual(f"Under {output_root}/:", f_out)
    print_table(["Folder", "Dir", "Kind", "Optical.Vector", "Files"], folder_rows, f_out)
    for problem in pseudo_problems:
        print_dual(warn(f"no pseudopotential file for {problem}: add it to the folder (or rerun with -p) "
                         "before running SIESTA."), f_out)

    # ------------------------------------------------------------------ [7]
    print_section('[7] SUMMARY & NEXT STEPS', f_out)
    n_diag = sum(1 for a in args.directions if a not in _OFFDIAG_REQUIRES)
    n_off = len(args.directions) - n_diag
    print_dual(f"{len(args.directions)} direction folder(s) written under '{output_root}' "
                f"({n_diag} diagonal, {n_off} biaxial).", f_out)
    print_dual(f"Cost                 : {len(args.directions)} independent SCF + optical runs "
                f"(optical mesh {n_opt_k} k-points each).", f_out)
    if report_path:
        print_dual(f"Report               : {report_path}", f_out)
    print_dual(color_text("\nNext steps:", 'yellow'), f_out)
    print_dual(f"  1. Run SIESTA in every folder, e.g.:", f_out)
    print_dual(f"       for d in {output_root}/dir_*; do (cd $d && mpirun -np 4 siesta calc.fdf > calc.out); done", f_out)
    print_dual("  2. Check that each folder finished and wrote <SystemLabel>.EPSIMG, then run:", f_out)
    extra_flags = ""
    if is_2d:
        extra_flags = " --dimensionality-correction --thickness <Ang>"
    elif sum(vacuum_axes) == 3:
        extra_flags = " --dimensionality-correction"
    print_dual(f"       stb-opticalAnalysis --directory {output_root} --energy-max 20{extra_flags}", f_out)
    print_dual("     (--plot-quantity all writes eps, n/k, alpha, R, loss and sigma plots; --energy-max "
                "trims the energy range of the reports and plots)", f_out)
    if n_off:
        print_dual("  3. stb-opticalAnalysis reconstructs eps_ij(E) automatically for any biaxial "
                    "direction whose matching diagonal pair is also present -- see Section [5] above.",
                    f_out)
    print_dual(color_text(
        "\nConvergence checklist (what changes the spectrum): Optical.Mesh (peak positions/heights), "
        "PAO.EnergyShift and basis size (conduction bands), SCF k-grid, Optical.NumberOfBands (energy "
        "range). Re-run Stage 1 with a different value and compare the Stage 2 results.", 'cyan'), f_out)
    if warnings:
        print_dual(color_text(f"\n{len(warnings)} warning(s) above -- read them before running SIESTA:", 'yellow'), f_out)
        for w in warnings:
            print_dual(f"  - {w}", f_out)

    if f_out:
        f_out.close()

    print("\n[INFO] Complete job!")
    print("\n" + "-" * 60)
    print(color_text("Direction folders ready for Stage 2 (stb-opticalAnalysis).\n", 'bold'))


if __name__ == "__main__":
    main()
