#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""Real-space |psi|^2 of one chosen electronic eigenstate (band N at
k-point K), from a SIESTA .WFSX file, written as a Gaussian .cube file
(primary output, for VESTA/VMD/Avogadro) plus an optional 2D-slice/1D-profile
gnuplot output. Useful for visualizing a single Bloch state -- e.g. "what
does the CBM look like at Gamma", "is this defect state localized" -- as a
diagnostic, NOT a physical ground-state observable (that already exists via
stb-density's .RHO-based total charge density, which sums over all occupied
states/k with proper weights).

Needs a real per-orbital radial basis -- verified live that this MUST come
from a <label>.fdf file with its .ion/.ion.xml files sitting alongside it,
never from a <label>.HSX/.XV: a Geometry built from .HSX has every orbital's
radial cutoff set to -1 (no shape at all, just enough bookkeeping for the
Hamiltonian's linear algebra), and a .XV-derived Geometry has no orbitals at
all. .HSX is not needed at all for this tool (unlike stb-fatbands, which
needs it for norm2/Sk); reading the .WFSX with a bare .fdf-derived Geometry
as `parent` gives numerically identical wavefunction densities.

Known limitations:
 - The given .fdf must reflect the EXACT basis set (PAO.Basis/
   PAO.BasisSize/etc.) used in the calculation that wrote the .WFSX --
   only the orbital COUNT is guarded, not basis identity; a coincidentally
   equal-count-but-different basis would silently misassign orbitals.
 - --band vbm/cbm materializes the whole .WFSX in memory to search for the
   global extremum -- fine for a band-path or modest mesh, potentially
   slow/memory-heavy for a large full-BZ WFSX; use --k-index/--band N
   directly for large meshes.
 - Only the first spinor component (spinor=0) is evaluated --
   non-collinear/SOC wavefunctions are not fully represented.
 - No file metadata ties a given .WFSX to a given .fdf beyond the
   orbital-count guard (same trust model as stb-fatbands).
 - SIESTA's numerical atomic orbitals (PAOs) have a finite, hard cutoff
   radius (confined basis) -- the density is exactly zero beyond a few
   Angstrom from every atom (unlike a plane-wave/Gaussian basis), so an
   overly large grid mostly shows empty space by construction, not a
   numerical artifact.
"""

VERSION = "1.0.0"

import argparse
import os
import numpy as np

from stb.core.cli import color_text, show_intro
from stb.core.deps import require_sisl
from stb.core.siesta_bands import read_data as read_bands_data, select_band_vbm_cbm
from stb.core.siesta_wfsx import resolve_wfsx_path, load_parent, read_wfsx_states
from stb.core.grid_export import write_data_file, write_profile_data_file, check_planar_orthogonality


def main():
    parser = argparse.ArgumentParser(
        description="Real-space |psi|^2 of one chosen (k, band) eigenstate from a "
                    "SIESTA .WFSX file -> Gaussian .cube.",
        epilog="Example usage:\n"
               "  stb-wfdensity --label siesta --k-index 0 --band 1\n"
               "  stb-wfdensity --label siesta --band vbm --fermi -4.2\n",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("--label", type=str, default=None,
                        help="SIESTA output label. Auto-detects <label>.WFSX, falling back to "
                             "<label>.selected.WFSX (written by WriteWaveFunctions T + %%block "
                             "WaveFuncKPoints) then <label>.bands.WFSX, with a note -- and "
                             "<label>.fdf (with its .ion/.ion.xml files alongside it -- NOT "
                             "<label>.XV/.HSX, neither of which carries real orbital shapes). "
                             "Mutually exclusive with --wfsx/--geometry-file.")
    parser.add_argument("--wfsx", type=str, default=None,
                        help="Explicit path to the .WFSX file (alternative to --label).")
    parser.add_argument("--geometry-file", type=str, default=None,
                        help="Explicit .fdf geometry path (must have the basis's .ion/.ion.xml "
                             "files alongside it). Required if --wfsx is used instead of --label.")

    k_group = parser.add_mutually_exclusive_group()
    k_group.add_argument("--k-index", type=int, default=None,
                        help="0-based index into the .WFSX's own k-point order (default: 0 if "
                             "neither --k-index nor --k-point is given).")
    k_group.add_argument("--k-point", type=float, nargs=3, default=None, metavar=("KX", "KY", "KZ"),
                        help="Match by k-point vector instead of index (same units as the "
                             ".WFSX's own k, matched via sisl's read_eigenstate).")
    parser.add_argument("--spin", type=int, default=0, choices=[0, 1],
                        help="Spin channel (default: 0; must be 0 for a non-polarized "
                             "calculation).")
    parser.add_argument("--band", type=str, default=None, required=True,
                        help="Which band at the chosen k: an integer N (1-based), or 'vbm'/'cbm' "
                             "(requires --fermi or --bands-file).")

    parser.add_argument("--fermi", type=float, default=None,
                        help="Fermi energy (eV), required if --band is vbm/cbm and "
                             "--bands-file isn't given.")
    parser.add_argument("--bands-file", type=str, default=None,
                        help="Optional companion .bands file to read the Fermi energy from "
                             "automatically, for --band vbm/cbm. Convenience only -- no "
                             "correspondence guard against --wfsx like stb-fatbands has.")
    parser.add_argument("--gap-tol", type=float, default=0.01,
                        help="Energy tolerance in eV for VBM/CBM classification (default: 0.01). "
                             "Same meaning as stb-bands/stb-fatbands.")

    grid_group = parser.add_mutually_exclusive_group()
    grid_group.add_argument("--spacing", type=float, default=0.3,
                        help="Grid spacing in Angstrom (default: 0.3). Mutually exclusive with "
                             "--shape.")
    grid_group.add_argument("--shape", type=int, nargs=3, default=None, metavar=("NX", "NY", "NZ"),
                        help="Explicit grid shape (voxel counts), alternative to --spacing.")

    parser.add_argument("--cube-only", action="store_true",
                        help="Skip the optional 2D-slice/1D-profile gnuplot output, write only "
                             "the .cube file.")
    parser.add_argument("--axis", type=int, default=2, choices=[0, 1, 2],
                        help="Axis for the optional gnuplot slice/profile (default: 2).")
    parser.add_argument("--profile", action="store_true",
                        help="Write a planar-averaged 1D profile instead of a 2D slice, for the "
                             "optional gnuplot output.")

    parser.add_argument("-o", "--output-dir", type=str, default=".",
                        help="Directory to write the .cube (and optional slice/profile) into "
                             "(default: current directory). Created if it doesn't exist.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-wfdensity {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.label:
        if args.wfsx or args.geometry_file:
            parser.error("--label cannot be combined with --wfsx/--geometry-file.")
        args.wfsx = resolve_wfsx_path(args.label, suffixes=(".WFSX", ".selected.WFSX", ".bands.WFSX"))
    elif not args.wfsx:
        parser.error("one of --label or --wfsx is required.")
    elif not args.geometry_file:
        parser.error("--geometry-file is required when using --wfsx instead of --label.")

    if args.wfsx is None or not os.path.isfile(args.wfsx):
        tried = (f"'{args.label}.WFSX', '{args.label}.selected.WFSX' and "
                 f"'{args.label}.bands.WFSX'") if args.label else f"'{args.wfsx}'"
        parser.error(f"No .WFSX file found ({tried}).")

    band_vbm_cbm = args.band.lower() in ("vbm", "cbm")
    if not band_vbm_cbm:
        try:
            band_n = int(args.band)
            if band_n < 1:
                raise ValueError
        except ValueError:
            parser.error("--band must be a positive integer (1-based), or 'vbm'/'cbm'.")
    if band_vbm_cbm and args.fermi is None and not args.bands_file:
        parser.error("--band vbm/cbm requires --fermi or --bands-file.")
    if args.k_point is not None and band_vbm_cbm:
        parser.error("--band vbm/cbm searches the whole k-mesh; --k-index/--k-point don't apply.")

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("WAVEFUNCTION DENSITY:", 'bold'))
    print("-" * 60)

    sisl = require_sisl()

    print("[INFO] Resolving geometry (needs real per-orbital basis, from .fdf+ion files) ...")
    geometry, geo_source, _ = load_parent(args.label, None, args.geometry_file, mode="orbitals_required")
    if geometry is None:
        parser.error("No usable geometry found -- stb-wfdensity needs <label>.fdf (with its "
                     ".ion/.ion.xml files alongside it) or --geometry-file. A bare .XV/.HSX "
                     "cannot provide the real orbital shapes wavefunction evaluation needs.")
    print(f"[INFO] Using '{geo_source}' ({geometry.no} orbitals).")

    print(f"[INFO] Reading wavefunction coefficients from '{args.wfsx}' ...")
    sizes, states_by_k = read_wfsx_states(args.wfsx, geometry)

    if sizes.no_u != geometry.no:
        parser.error(
            f"Orbital count mismatch: '{args.wfsx}' has {sizes.no_u} orbitals but the "
            f"loaded geometry ('{geo_source}') has {geometry.no} -- they must be from "
            "the same calculation/basis set."
        )

    nspin = sizes.nspin if sizes.nspin == 2 else 1

    if band_vbm_cbm:
        fermi_energy = args.fermi
        if fermi_energy is None:
            fermi_energy, _, _, _ = read_bands_data(args.bands_file)
        print(f"[INFO] Searching whole k-mesh for the global {args.band.upper()} "
              f"(Fermi = {fermi_energy:.6f} eV) ...")
        k_index, spin, band_index = select_band_vbm_cbm(
            states_by_k, nspin, fermi_energy, args.gap_tol, args.band.lower()
        )
        print(f"[INFO] {args.band.upper()} found at k-index {k_index}, spin {spin}, "
              f"band {band_index + 1} (1-based).")
        full_state = states_by_k[k_index][spin]
    else:
        band_index = band_n - 1
        spin = args.spin
        if args.k_point is not None:
            k_index = None
            for i, block in enumerate(states_by_k):
                state = block.get(spin)
                if state is not None and np.allclose(state.info.get("k"), args.k_point, atol=1e-4):
                    k_index = i
                    break
            if k_index is None:
                parser.error(f"No k-point matching {args.k_point} (spin {spin}) found in "
                             f"'{args.wfsx}'.")
            full_state = states_by_k[k_index][spin]
        else:
            k_index = args.k_index if args.k_index is not None else 0
            if k_index < 0 or k_index >= len(states_by_k):
                parser.error(f"--k-index {k_index} out of range (0-{len(states_by_k) - 1}).")
            block = states_by_k[k_index]
            if spin not in block:
                parser.error(f"--spin {spin} not present at k-index {k_index}.")
            full_state = block[spin]
        nbands_here = full_state.state.shape[0]
        if band_index < 0 or band_index >= nbands_here:
            parser.error(f"--band {band_n} out of range (1-{nbands_here}) at this k-point.")

    es_one = full_state.sub(band_index, inplace=False)
    k_vec = es_one.info.get("k", (0.0, 0.0, 0.0))
    eig = float(np.asarray(es_one.eig).reshape(-1)[0])
    print(f"[INFO] Selected state: k={np.asarray(k_vec)}, spin={spin}, "
          f"band={band_index + 1}, eigenvalue={eig:.6f} eV.")

    shape_or_spacing = tuple(args.shape) if args.shape else args.spacing
    grid = sisl.Grid(shape_or_spacing, geometry=geometry, dtype=complex)
    print(f"[INFO] Evaluating wavefunction on a {grid.grid.shape[0]}x{grid.grid.shape[1]}x"
          f"{grid.grid.shape[2]} grid (this may take a moment) ...")
    sisl.physics.electron.wavefunction(es_one.state, grid, geometry=geometry, k=k_vec,
                                       spinor=0, eta=True)

    density = np.abs(grid.grid) ** 2
    norm_check = float(density.sum() * grid.dvolume)
    print(f"[INFO] Normalization check (sum|psi|^2 * dvol, should be ~1): {norm_check:.6f}")
    if abs(norm_check - 1.0) > 0.05:
        print(color_text(
            "[WARNING] Normalization deviates from 1 by more than 5% -- the grid may be too "
            "coarse (orbitals clipped) or the geometry/basis may not match the .WFSX. Consider "
            "a finer --spacing.", 'yellow'))

    os.makedirs(args.output_dir, exist_ok=True)
    label_tag = f"k{k_index}_b{band_index + 1}"
    cube_path = os.path.join(args.output_dir, f"wfdensity_{label_tag}.cube")
    grid.grid = density
    grid.set_geometry(geometry)
    grid.write(cube_path)
    print(f"[INFO] Cube file saved to '{color_text(cube_path, 'green')}' "
          "(open in VESTA/VMD/Avogadro for real 3D isosurfaces)")

    if not args.cube_only:
        lattice = geometry.cell
        origin = grid.origin
        if args.profile:
            out_path = os.path.join(args.output_dir, f"wfdensity_{label_tag}_profile.dat")
            write_profile_data_file(density, lattice, args.axis, out_path)
        else:
            dim_size = density.shape[args.axis]
            slice_idx = dim_size // 2
            angle = check_planar_orthogonality(lattice, args.axis)
            if abs(angle - 90.0) > 1.0:
                print(color_text(
                    f"[WARNING] The cut plane is skewed ({angle:.1f} deg) -- see stb-density's "
                    "own note on hexagonal/skewed cells.", 'yellow'))
            out_path = os.path.join(args.output_dir, f"wfdensity_{label_tag}_slice.dat")
            write_data_file(density, lattice, origin, out_path, mode='slice',
                            slice_idx=slice_idx, axis=args.axis)
        print(f"[INFO] Wrote {out_path}")

    print("\n[INFO] Complete job!")
    print("\n" + "-" * 60)


if __name__ == "__main__":
    main()
