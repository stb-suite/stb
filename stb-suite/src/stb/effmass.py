#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""Effective mass tensor and band velocity for one chosen electronic
eigenstate (band N at k-point K), from a SIESTA .WFSX + .HSX, via sisl's
EigenstateElectron.effective_mass()/.velocity() (analytic k-derivatives of
the Hamiltonian, sisl.physics.electron.derivative()).

Needs a real sisl.physics.Hamiltonian (.HSX) -- unlike stb-wfdensity/
stb-sts, there is no approximate fallback: the analytic dH/dk, d2H/dk2
derivatives that effective_mass()/velocity() need only exist on a real
Hamiltonian (same requirement as stb-coop's COOP/COHP).

IMPORTANT (verified live, matches sisl's own docstring warning): the 2nd-
order curvature correction inside effective_mass() uses energy differences
between ALL bands at that k-point, so it must be computed on the FULL
multi-band eigenstate, never after `.sub(band_index)` -- doing so silently
gives wildly wrong values (verified: sub-then-calculate differed from the
correct full-state calculation by orders of magnitude on a real test
system, e.g. 3.4e7 vs 3.6e4 for the same Voigt component). This tool always
calls effective_mass()/velocity() on the full state and only indexes the
band of interest afterward.

Unit conversion for effective_mass(): sisl returns it in its own internal
"(ps/Ang)^2" convention. Derived and cross-checked independently via two
methods (hand algebra using sisl's own `_velocity_const = 1/hbar[eV*ps]`,
and scipy.constants) -- both agree to 6 significant figures:
    m*[m0] = 3.301296e-6 * effective_mass_sisl_output
Sanity-checked against a real value from this suite's own Sn3O4 test
system: 134638.3 (sisl units) -> 0.444 m0, a physically ordinary value
(most real materials: 0.05-5 m0).

velocity() is documented by sisl itself as "Ang/ps"; converted here to
km/s (1 Ang/ps = 1e-10 m / 1e-12 s = 100 m/s = 0.1 km/s).

Known limitations:
 - Only a single (k, band) point, not swept along a path -- effective
   mass away from a true band extremum (VBM/CBM/saddle point) is not a
   very meaningful quantity physically; --band vbm/cbm is the intended
   common use, --k-index/--band N for anything else.
 - Directions along a non-periodic (vacuum-padded) axis will read exactly
   0 (sisl's own documented behavior: "Since some directions may not be
   periodic there will be zeros").
 - Berry curvature correction to the velocity is NOT included (sisl's own
   documented limitation of .velocity() -- see its docstring), and this
   tool does not implement Berry curvature at all (out of scope).
 - Effective mass (but NOT velocity) is unavailable for non-collinear/
   SOC Hamiltonians (nspin=4/8): verified live that sisl's effective_mass()
   crashes with a bare "TypeError: 'NoneType' object is not subscriptable"
   for those cases (its 2nd-order k-derivative needs the overlap matrix's
   own 2nd derivative, ddSk, which isn't implemented there in this sisl
   version) -- this tool detects nspin up front and reports velocity only
   in that case, with a warning, rather than letting that traceback surface.
"""

VERSION = "1.0.0"

import argparse
import os
import numpy as np

from stb.core.cli import color_text, show_intro
from stb.core.deps import require_sisl
from stb.core.siesta_bands import read_data as read_bands_data, select_band_vbm_cbm
from stb.core.siesta_wfsx import resolve_wfsx_path, load_parent, read_wfsx_states

# Derived + cross-checked (hand algebra using sisl's own _velocity_const,
# and independently via scipy.constants): see module docstring.
EFFMASS_TO_M0 = 3.301296009086774e-06
ANG_PS_TO_KM_S = 0.1

VOIGT_LABELS = ("xx", "yy", "zz", "yz", "xz", "xy")


def main():
    parser = argparse.ArgumentParser(
        description="Effective mass tensor (Voigt) and band velocity for one chosen (k, band) "
                    "eigenstate, from a SIESTA .WFSX + .HSX.",
        epilog="Example usage:\n"
               "  stb-effmass --label siesta --band vbm --fermi -4.2\n"
               "  stb-effmass --label siesta --k-index 0 --band 1\n",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("--label", type=str, default=None,
                        help="SIESTA output label. Auto-detects <label>.WFSX, falling back to "
                             "<label>.selected.WFSX then <label>.bands.WFSX, and <label>.HSX "
                             "(mandatory -- no approximate fallback). Mutually exclusive with "
                             "--wfsx/--hsx-file.")
    parser.add_argument("--wfsx", type=str, default=None,
                        help="Explicit path to the .WFSX file (alternative to --label).")
    parser.add_argument("--hsx-file", type=str, default=None,
                        help="Explicit .HSX path. Required if --wfsx is used instead of --label.")

    k_group = parser.add_mutually_exclusive_group()
    k_group.add_argument("--k-index", type=int, default=None,
                        help="0-based index into the .WFSX's own k-point order (default: 0 if "
                             "neither --k-index nor --k-point is given).")
    k_group.add_argument("--k-point", type=float, nargs=3, default=None, metavar=("KX", "KY", "KZ"),
                        help="Match by k-point vector instead of index (matched via sisl's "
                             "read_eigenstate).")
    parser.add_argument("--spin", type=int, default=0, choices=[0, 1],
                        help="Spin channel (default: 0; must be 0 for a non-polarized "
                             "calculation).")
    parser.add_argument("--band", type=str, default=None, required=True,
                        help="Which band at the chosen k: an integer N (1-based), or 'vbm'/'cbm' "
                             "(requires --fermi).")
    parser.add_argument("--fermi", type=float, default=None,
                        help="Fermi energy (eV), required if --band is vbm/cbm.")
    parser.add_argument("--gap-tol", type=float, default=0.01,
                        help="Energy tolerance in eV for VBM/CBM classification (default: 0.01). "
                             "Same meaning as stb-bands/stb-fatbands/stb-wfdensity.")

    parser.add_argument("-o", "--output-dir", type=str, default=".",
                        help="Directory to write effmass.txt into (default: current directory). "
                             "Created if it doesn't exist.")
    parser.add_argument("-v", "--version", action="version", version=f"stb-effmass {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.label:
        if args.wfsx or args.hsx_file:
            parser.error("--label cannot be combined with --wfsx/--hsx-file.")
        args.wfsx = resolve_wfsx_path(args.label, suffixes=(".WFSX", ".selected.WFSX", ".bands.WFSX"))
    elif not args.wfsx:
        parser.error("one of --label or --wfsx is required.")
    elif not args.hsx_file:
        parser.error("--hsx-file is required when using --wfsx instead of --label.")

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
    if band_vbm_cbm and args.fermi is None:
        parser.error("--band vbm/cbm requires --fermi.")
    if args.k_point is not None and band_vbm_cbm:
        parser.error("--band vbm/cbm searches the whole k-mesh; --k-index/--k-point don't apply.")

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("EFFECTIVE MASS / BAND VELOCITY:", 'bold'))
    print("-" * 60)

    sisl = require_sisl()

    print("[INFO] Resolving Hamiltonian (mandatory -- effective mass needs a real dH/dk) ...")
    H, source, _ = load_parent(args.label, args.hsx_file, None, mode="hamiltonian_required")
    if H is None:
        parser.error("No .HSX found -- stb-effmass needs --hsx-file, or --label with a "
                     "<label>.HSX next to it. There is no approximate fallback.")
    print(f"[INFO] Using '{source}' ({H.geometry.no} orbitals).")

    print(f"[INFO] Reading wavefunction coefficients from '{args.wfsx}' ...")
    sizes, states_by_k = read_wfsx_states(args.wfsx, H)
    nspin = sizes.nspin if sizes.nspin == 2 else 1

    if band_vbm_cbm:
        print(f"[INFO] Searching whole k-mesh for the global {args.band.upper()} "
              f"(Fermi = {args.fermi:.6f} eV) ...")
        k_index, spin, band_index = select_band_vbm_cbm(
            states_by_k, nspin, args.fermi, args.gap_tol, args.band.lower()
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

    k_vec = full_state.info.get("k", (0.0, 0.0, 0.0))
    eig = float(np.asarray(full_state.eig)[band_index])
    print(f"[INFO] Selected state: k={np.asarray(k_vec)}, spin={spin}, "
          f"band={band_index + 1}, eigenvalue={eig:.6f} eV.")

    print("[INFO] Computing effective mass / velocity on the FULL multi-band state "
          "(never after isolating one band -- see module docstring) ...")
    v_full = full_state.velocity().real  # (3, nbands), Ang/ps
    v_kms = v_full[:, band_index] * ANG_PS_TO_KM_S

    # Verified live: sisl's effective_mass() (2nd-order k-derivative, needs
    # the overlap matrix's own 2nd derivative, ddSk) crashes with a bare
    # "TypeError: 'NoneType' object is not subscriptable" for a
    # non-collinear/SOC Hamiltonian (nspin=4/8) in this sisl version --
    # ddSk simply isn't implemented for that case. velocity() (1st-order)
    # works fine there. Rather than let that raw traceback surface, detect
    # it up front and degrade to velocity-only for nspin=4/8.
    em_voigt = None
    if sizes.nspin in (4, 8):
        print(color_text(
            f"[WARNING] Effective mass (2nd-order k-derivative) is not supported by sisl for "
            f"non-collinear/spin-orbit Hamiltonians (nspin={sizes.nspin}) in this environment "
            "-- reporting velocity only.", 'yellow'))
    else:
        em_voigt_full = full_state.effective_mass().real  # (6, nbands), (ps/Ang)^2
        em_voigt = em_voigt_full[:, band_index] * EFFMASS_TO_M0

    lines = []
    lines.append("=" * 60)
    lines.append("EFFECTIVE MASS / BAND VELOCITY - STB Suite")
    lines.append("=" * 60)
    lines.append(f"Source (.WFSX)   : {args.wfsx}")
    lines.append(f"Source (.HSX)    : {source}")
    lines.append(f"k-point          : {np.asarray(k_vec)}")
    lines.append(f"Spin             : {spin}")
    lines.append(f"Band (1-based)   : {band_index + 1}")
    lines.append(f"Eigenvalue       : {eig:.6f} eV")
    lines.append("-" * 60)
    if em_voigt is not None:
        lines.append("Effective mass tensor (Voigt, in units of the free-electron mass m0):")
        for label, value in zip(VOIGT_LABELS, em_voigt):
            lines.append(f"  m*_{label:<3} = {value:12.6f} m0")
    else:
        lines.append(f"Effective mass tensor: N/A (unsupported for nspin={sizes.nspin}, see warning above)")
    lines.append("-" * 60)
    lines.append("Band velocity (km/s):")
    for axis, value in zip("xyz", v_kms):
        lines.append(f"  v_{axis} = {value:12.6f} km/s")
    lines.append("=" * 60)
    lines.append("Note: a component reading exactly 0 along a vacuum-padded (non-periodic) "
                 "axis is expected, not an error.")
    report = "\n".join(lines) + "\n"

    print("\n" + report)

    os.makedirs(args.output_dir, exist_ok=True)
    out_path = os.path.join(args.output_dir, "effmass.txt")
    with open(out_path, "w") as f:
        f.write(report)
    print(f"[INFO] Wrote {out_path}")

    print("\n[INFO] Complete job!")
    print("\n" + "-" * 60)


if __name__ == "__main__":
    main()
