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

**A real, significant physics finding (v2.0.0), verified live**: sisl's
own `effective_mass()` returns each Voigt component as an INDEPENDENT,
ELEMENT-WISE reciprocal of the corresponding curvature component
(`velocity_const**2 / d2E/dk_i dk_j`, computed separately per Voigt
index) -- NOT a proper matrix inversion of the full 3x3 curvature tensor.
This means m*_xx/yy/zz are only rigorously "the effective mass along x/y/z"
when the off-diagonal terms (yz, xz, xy) are negligible. Verified live on
this suite's own Sn3O4 VBM: the off-diagonal curvature terms were
COMPARABLE TO OR LARGER THAN the diagonal ones, and properly assembling
the full curvature tensor from all 6 Voigt values, inverting THAT matrix,
and diagonalizing gave principal effective masses of -1.320/-0.844/-0.343
m0 -- substantially different (30-100%+) from sisl's own naive per-axis
values (-0.629/-0.757/-0.515 m0) for the exact same state. This is common,
not a rare edge case: away from high-symmetry k-points/directions
(this VBM is not at Gamma), the curvature tensor is generically NOT
diagonal in the Cartesian basis. This tool now reports BOTH: sisl's own
per-axis values (for continuity/comparison) AND the properly diagonalized
principal effective masses + their principal-axis directions, with an
explicit flag when the two disagree meaningfully. The full-tensor
constant (velocity_const**2) is derived empirically from sisl's own
public-API output (matching a returned Voigt component against its own
un-inverted curvature) rather than reaching into sisl's private
`_velocity_const` attribute, so this has no dependency on sisl internals.

Known limitations:
 - Only a single (k, band) point, not swept along a path -- effective
   mass away from a true band extremum (VBM/CBM/saddle point) is not a
   very meaningful quantity physically; --band vbm/cbm is the intended
   common use, --k-index/--band N for anything else.
 - Directions along a non-periodic (vacuum-padded) axis will read exactly
   0 (sisl's own documented behavior: "Since some directions may not be
   periodic there will be zeros") -- the full curvature tensor is then
   singular along that axis and the principal-mass diagonalization is
   skipped (reported, not silently wrong).
 - Berry curvature correction to the velocity is NOT included (sisl's own
   documented limitation of .velocity() -- see its docstring), and this
   tool does not implement Berry curvature at all (out of scope).
 - Effective mass (but NOT velocity) is unavailable for ANY spin-resolved
   Hamiltonian -- collinear spin-polarized (nspin=2) AND non-collinear/SOC
   (nspin=4/8) alike. Verified live: sisl's own `Hamiltonian.ddPk()` (the
   2nd-order k-derivative machinery `effective_mass()`/`derivative(2)`
   need) simply does not accept a `spin` argument at all in this sisl
   version (0.16.4) -- confirmed directly from its signature -- while
   `dPk()` (1st-order, used by `velocity()`) does. A real user hit this
   live on their own spin-polarized calculation: `effective_mass()`
   crashed with `TypeError: SparseOrbitalBZ._ddPk() got an unexpected
   keyword argument 'spin'`. This tool previously only guarded nspin=4/8
   (assuming nspin=2 worked); it now detects ANY nspin != 1 up front and
   degrades to velocity-only with a clear warning, plus a defensive
   try/except around the actual call as a second line of defense in case
   a differently-behaving sisl version raises this for some other reason.
"""

VERSION = "2.0.0"

import argparse
import os
from datetime import datetime

import numpy as np

from stb.core import citations
from stb.core.cli import color_text, show_intro, print_dual, print_section, print_table
from stb.core.deps import require_sisl
from stb.core.siesta_bands import select_band_vbm_cbm, resolve_fermi_energy_hierarchy
from stb.core.siesta_wfsx import resolve_wfsx_path, load_parent, read_wfsx_states

REPORT_FILE = "stb_effmass_report.txt"
BIB_FILE = "references.bib"

# Derived + cross-checked (hand algebra using sisl's own _velocity_const,
# and independently via scipy.constants): see module docstring.
EFFMASS_TO_M0 = 3.301296009086774e-06
ANG_PS_TO_KM_S = 0.1

VOIGT_LABELS = ("xx", "yy", "zz", "yz", "xz", "xy")

# If the largest off-diagonal curvature term exceeds this fraction of the
# mean |diagonal| term, sisl's own naive per-axis Voigt values are flagged
# as potentially misleading (see the module docstring's live verification:
# a real VBM state showed off-diagonal terms LARGER than the diagonal
# ones, a 30-100%+ difference once properly diagonalized).
OFFDIAG_WARN_RATIO = 0.2


def compute_principal_effective_mass(full_state, band_index, em_voigt_sisl, ddv_voigt):
    """Returns (principal_masses_m0, principal_axes, offdiag_ratio) or
    (None, None, None) if the full 3x3 curvature tensor is singular (a
    non-periodic/vacuum-padded axis -- see the module docstring).

    Assembles the FULL symmetric curvature tensor from all 6 Voigt
    components (not just the diagonal ones sisl's own per-axis values
    use independently), inverts THAT matrix, and diagonalizes it --
    the properly matrix-inverted, basis-independent principal effective
    masses, as opposed to sisl's own element-wise-reciprocal Voigt
    values. The velocity_const**2 scale factor is derived empirically
    from sisl's own already-computed em_voigt_sisl/ddv_voigt pair
    (`em = velocity_const**2 / ddv` per component) rather than reaching
    into sisl's private `_velocity_const` attribute.
    """
    C = np.array([
        [ddv_voigt[0], ddv_voigt[5], ddv_voigt[4]],
        [ddv_voigt[5], ddv_voigt[1], ddv_voigt[3]],
        [ddv_voigt[4], ddv_voigt[3], ddv_voigt[2]],
    ])
    diag = np.abs(ddv_voigt[:3])
    offdiag = np.abs(ddv_voigt[3:])
    mean_diag = diag.mean()
    offdiag_ratio = float(offdiag.max() / mean_diag) if mean_diag > 0 else None

    try:
        vc2_candidates = [em_voigt_sisl[i] * ddv_voigt[i] for i in range(6) if abs(ddv_voigt[i]) > 1e-12]
        if not vc2_candidates:
            return None, None, offdiag_ratio
        vc2 = float(np.median(vc2_candidates))
        Minv_full = np.linalg.inv(C) * vc2
    except np.linalg.LinAlgError:
        return None, None, offdiag_ratio

    evals, evecs = np.linalg.eigh(Minv_full)
    return evals * EFFMASS_TO_M0, evecs, offdiag_ratio


def main():
    parser = argparse.ArgumentParser(
        description="Effective mass tensor (Voigt) and band velocity for one chosen (k, band) "
                    "eigenstate, from a SIESTA .WFSX + .HSX.",
        epilog="Example usage:\n"
               "  stb-effmass --label siesta --band vbm --fermi -4.2\n"
               "  stb-effmass --label siesta --k-index 0 --band 1\n"
               "  stb-effmass --label siesta --band vbm --fermi-file calc.out \\\n"
               "      --save-report --save-gnuplot --view\n",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("--label", type=str, default=None,
                        help="SIESTA output label. Auto-detects <label>.WFSX, falling back to "
                             "<label>.selected.WFSX then <label>.bands.WFSX, and, unless "
                             "--hsx-file overrides it, <label>.HSX (mandatory -- no approximate "
                             "fallback). Mutually exclusive with --wfsx (not --hsx-file -- see "
                             "its own help).")
    parser.add_argument("--wfsx", type=str, default=None,
                        help="Explicit path to the .WFSX file (alternative to --label).")
    parser.add_argument("--hsx-file", type=str, default=None,
                        help="Explicit .HSX path. Required if --wfsx is used instead of --label; "
                             "optional (and common) together WITH --label too, since the real "
                             "file is often not literally named <label>.HSX.")

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
                             "(requires --fermi, --bands-file, --fermi-file, or an auto-detected "
                             ".out log -- see below).")

    parser.add_argument("--fermi", type=float, default=None,
                        help="Fermi energy (eV) for --band vbm/cbm. Highest-priority source if given.")
    parser.add_argument("--bands-file", type=str, default=None,
                        help="Companion .bands file to read the Fermi energy from, for --band "
                             "vbm/cbm.")
    parser.add_argument("--fermi-file", type=str, default=None,
                        help="Explicit SIESTA .out log to read the Fermi energy from, for --band "
                             "vbm/cbm -- an alternative to --bands-file for a run with no saved "
                             ".bands file. Not assumed to be named after --label.")
    parser.add_argument("--gap-tol", type=float, default=0.01,
                        help="Energy tolerance in eV for VBM/CBM classification (default: 0.01). "
                             "Same meaning as stb-bands/stb-fatbands/stb-wfdensity.")

    parser.add_argument("-o", "--output-dir", type=str, default=".",
                        help="Directory to write effmass.dat/velocity.dat (with --save-gnuplot) "
                             "and stb_effmass_report.txt/references.bib into (default: current "
                             "directory). Created if it doesn't exist.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the full run report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--save-gnuplot", action="store_true",
                        help="Also write effmass.dat/velocity.dat + a real .gplot bar-chart "
                             "script. Off by default.")
    parser.add_argument("--view", action="store_true",
                        help="Show a matplotlib bar-chart preview (effective mass + velocity "
                             "components) before finishing. Off by default.")

    parser.add_argument("-v", "--version", action="version", version=f"stb-effmass {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    if args.label:
        if args.wfsx:
            parser.error("--label cannot be combined with --wfsx.")
        args.wfsx = resolve_wfsx_path(args.label, suffixes=(".WFSX", ".selected.WFSX", ".bands.WFSX"))
    elif not args.wfsx:
        parser.error("one of --label or --wfsx is required.")
    elif not args.hsx_file:
        parser.error("--hsx-file is required when using --wfsx instead of --label.")
    # --label + --hsx-file together is valid and common -- see the module docstring (the same
    # fix already made for stb-wfdensity/stb-sts/stb-coop): load_parent() below already prefers
    # an explicit hsx_file over <label>.HSX on its own.

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
    if args.k_point is not None and band_vbm_cbm:
        parser.error("--band vbm/cbm searches the whole k-mesh; --k-index/--k-point don't apply.")

    fermi_energy = fermi_source = None
    if band_vbm_cbm:
        fermi_energy, fermi_source = resolve_fermi_energy_hierarchy(
            args.fermi, args.bands_file, args.fermi_file, args.label)
        if fermi_energy is None:
            parser.error(
                "--band vbm/cbm needs a Fermi energy -- none found. Pass --fermi <value>, "
                "--bands-file <path>, --fermi-file <path>, or leave a SIESTA .out log "
                "(<label>.out, or the sole *.out) in the current directory."
            )

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite",
            "A comprehensive toolkit for SIESTA DFT simulations",
            f"Version {VERSION} | University of Brasilia - 2025",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("EFFMASS: reading input data", 'bold'))
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
              f"(Fermi = {fermi_energy:.6f} eV, source: {fermi_source}) ...")
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

    k_vec = np.asarray(full_state.info.get("k", (0.0, 0.0, 0.0)))
    eig = float(np.asarray(full_state.eig)[band_index])
    print(f"[INFO] Selected state: k={k_vec}, spin={spin}, "
          f"band={band_index + 1}, eigenvalue={eig:.6f} eV.")

    print("[INFO] Computing effective mass / velocity on the FULL multi-band state "
          "(never after isolating one band -- see module docstring) ...")
    v_full = full_state.velocity().real  # (3, nbands), Ang/ps
    v_kms = v_full[:, band_index] * ANG_PS_TO_KM_S

    # Verified live: sisl's effective_mass()/derivative(2) (2nd-order
    # k-derivative machinery) crashes for ANY spin-resolved Hamiltonian in
    # this sisl version -- ddPk() simply has no 'spin' parameter at all
    # (confirmed directly from its signature), while dPk() (1st-order,
    # used by velocity() above) does. This affects collinear
    # spin-polarized (nspin=2) as much as non-collinear/SOC (nspin=4/8) --
    # a real user hit exactly this on their own spin-polarized
    # calculation. Rather than let the raw TypeError traceback surface,
    # detect it up front (nspin != 1) and degrade to velocity-only, with a
    # try/except as a second line of defense in case some other sisl
    # version/configuration raises it for a different reason.
    em_voigt = principal_masses = principal_axes = offdiag_ratio = None
    ddv_voigt = None
    if sizes.nspin != 1:
        print(color_text(
            f"[WARNING] Effective mass (2nd-order k-derivative) is not supported by sisl for "
            f"spin-resolved Hamiltonians (nspin={sizes.nspin}) in this environment -- reporting "
            "velocity only.", 'yellow'))
    else:
        try:
            em_voigt_full = full_state.effective_mass().real  # (6, nbands), (ps/Ang)^2
            em_voigt_sisl_full = em_voigt_full[:, band_index]
            em_voigt = em_voigt_sisl_full * EFFMASS_TO_M0
            _dv, ddv_full = full_state.derivative(2)
            ddv_voigt = ddv_full.real[:, band_index]
            principal_masses, principal_axes, offdiag_ratio = compute_principal_effective_mass(
                full_state, band_index, em_voigt_sisl_full, ddv_voigt)
        except TypeError as e:
            print(color_text(
                f"[WARNING] Effective mass calculation failed unexpectedly ({e}) -- reporting "
                "velocity only.", 'yellow'))

    # --- From here on: the numbered, save-able report --------------------
    os.makedirs(args.output_dir, exist_ok=True)
    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    print_dual(color_text("\n===== STB-EFFMASS REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time      : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Label          : {args.label}" if args.label else "Label          : (explicit --wfsx/--hsx-file)", f_out)
    print_dual(f"WFSX file      : {args.wfsx}", f_out)
    print_dual(f"Hamiltonian src: {source}", f_out)
    print_dual(f"Band selection : {args.band}" + (f" (Fermi source: {fermi_source})" if band_vbm_cbm else ""), f_out)
    print_dual(f"Output dir     : {args.output_dir}", f_out)
    print_dual(f"Save gnuplot   : {'yes' if args.save_gnuplot else 'no'}", f_out)
    print_dual(f"View (matplotlib): {'yes' if args.view else 'no'}", f_out)

    print_section("[1] INPUT DATA", f_out)
    print_table(["Quantity", "Value"], [
        (["Orbitals (basis)", f"{H.geometry.no}"], None),
        (["k-points in .WFSX", f"{sizes.nk}"], None),
        (["Spin channels", f"{sizes.nspin} ({'non-polarized' if sizes.nspin == 1 else 'polarized' if sizes.nspin == 2 else 'non-collinear/SOC'})"], None),
    ], f_out)

    print_section("[2] STATE SELECTION", f_out)
    print_table(["Quantity", "Value"], [
        (["k-index", f"{k_index}"], None),
        (["k-vector", f"{k_vec}"], None),
        (["Spin", f"{spin}"], None),
        (["Band (1-based)", f"{band_index + 1}"], None),
        (["Eigenvalue", f"{eig:.6f} eV"], None),
    ], f_out)
    if band_vbm_cbm:
        print_dual(f"Selected via global {args.band.upper()} search "
                   f"(Fermi = {fermi_energy:.6f} eV, source: {fermi_source}, "
                   f"--gap-tol {args.gap_tol} eV).", f_out)

    print_section("[3] EFFECTIVE MASS (per-axis Voigt, sisl's own convention)", f_out)
    if em_voigt is not None:
        print_table(["Component", "Value (m0)"], [
            ([f"m*_{label}", f"{value:12.6f}"], None) for label, value in zip(VOIGT_LABELS, em_voigt)
        ], f_out)
        print_dual("(Each component is sisl's own INDEPENDENT, element-wise reciprocal of the "
                   "corresponding curvature term -- rigorously the effective mass along x/y/z "
                   "only when the off-diagonal terms below are small. See [4] for the properly "
                   "diagonalized principal effective masses.)", f_out)
    else:
        print_dual(f"N/A (unsupported for nspin={sizes.nspin}, see warning above).", f_out)

    print_section("[4] EFFECTIVE MASS (principal, full tensor)", f_out)
    if principal_masses is not None:
        print_dual(f"Off-diagonal/diagonal curvature ratio: {offdiag_ratio:.3f}", f_out)
        if offdiag_ratio > OFFDIAG_WARN_RATIO:
            print_dual(color_text(
                "[WARNING] Off-diagonal curvature is a significant fraction of (or exceeds) the "
                "diagonal curvature -- sisl's own per-axis Voigt values in [3] may be misleading "
                "for this state. Use the principal effective masses below instead.", 'yellow'), f_out)
        else:
            print_dual(color_text(
                "[OK] Off-diagonal curvature is small -- the per-axis Voigt values in [3] and "
                "the principal effective masses below should agree closely.", 'green'), f_out)
        print_table(["Principal mass (m0)", "Direction (unit vector)"], [
            ([f"{principal_masses[i]:12.6f}", f"{np.array2string(principal_axes[:, i], precision=4)}"], None)
            for i in range(3)
        ], f_out)
    elif em_voigt is not None:
        print_dual("N/A -- the full curvature tensor is singular (a non-periodic/vacuum-padded "
                   "axis at this k-point; see the module docstring).", f_out)
    else:
        print_dual(f"N/A (unsupported for nspin={sizes.nspin}, see warning above).", f_out)

    print_section("[5] BAND VELOCITY", f_out)
    print_table(["Component", "Value (km/s)"], [
        ([f"v_{axis}", f"{value:12.6f}"], None) for axis, value in zip("xyz", v_kms)
    ], f_out)
    print_dual(f"|v| = {float(np.linalg.norm(v_kms)):.6f} km/s", f_out)
    print_dual("(A component reading exactly 0 along a vacuum-padded/non-periodic axis is "
               "expected, not an error.)", f_out)

    print_section("[6] OUTPUT DATA & PLOTS", f_out)
    dat_paths = []
    gplot_path = None
    if args.save_gnuplot:
        vel_path = os.path.join(args.output_dir, "velocity.dat")
        with open(vel_path, "w") as f:
            f.write("#component\tvalue(km/s)\n")
            for axis, value in zip("xyz", v_kms):
                f.write(f"v_{axis}\t{value:.6f}\n")
        dat_paths.append(vel_path)

        gplot_lines = [
            'set terminal pdfcairo enhanced font "Arial,14" size 10,5\n',
        ]
        if em_voigt is not None:
            em_path = os.path.join(args.output_dir, "effmass.dat")
            with open(em_path, "w") as f:
                f.write("#component\tvalue(m0)\n")
                for label, value in zip(VOIGT_LABELS, em_voigt):
                    f.write(f"m*_{label}\t{value:.6f}\n")
            dat_paths.append(em_path)
            gplot_lines += [
                'set multiplot layout 1,2\n',
                'set style data histogram\n',
                'set style fill solid 0.7\n',
                'set boxwidth 0.6\n',
                'set grid ytics\n',
                'unset key\n',
                'set title "Effective mass (m0)"\n',
                'set ylabel "m* (m0)"\n',
                'plot "effmass.dat" using 2:xtic(1) lc rgb "#1f77b4"\n',
                'set title "Band velocity (km/s)"\n',
                'set ylabel "v (km/s)"\n',
                'plot "velocity.dat" using 2:xtic(1) lc rgb "#d62728"\n',
                'unset multiplot\n',
            ]
        else:
            gplot_lines += [
                'set style data histogram\n',
                'set style fill solid 0.7\n',
                'set boxwidth 0.6\n',
                'set grid ytics\n',
                'unset key\n',
                'set title "Band velocity (km/s)"\n',
                'set ylabel "v (km/s)"\n',
                'plot "velocity.dat" using 2:xtic(1) lc rgb "#d62728"\n',
            ]
        gplot_path = os.path.join(args.output_dir, "effmass.gplot")
        gplot_lines.insert(1, f'set output "effmass.pdf"\n')
        with open(gplot_path, 'w') as f:
            f.writelines(gplot_lines)
        print_dual(color_text(
            f"[OK] Data written: {', '.join(os.path.basename(p) for p in dat_paths)}.", 'green'), f_out)
        print_dual(color_text(f"[OK] Gnuplot script written to '{gplot_path}'.", 'green'), f_out)
    else:
        print_dual("Not written (off by default -- pass --save-gnuplot to write "
                   "effmass.dat/velocity.dat + effmass.gplot).", f_out)

    print_section("[7] REFERENCES", f_out)
    bib_entries = [citations.SIESTA, citations.SIESTA_RECENT]
    citations.write_bib_file(os.path.join(args.output_dir, BIB_FILE), bib_entries)
    print_dual(color_text(
        f"[OK] Citations for the methods used in this run written to "
        f"'{os.path.join(args.output_dir, BIB_FILE)}' ({len(bib_entries)} entries).", 'green'), f_out)

    print_section("[8] SUMMARY & FILES", f_out)
    print_dual("Status         : OK", f_out)
    for path in dat_paths:
        print_dual(f"Data           : {path}", f_out)
    if gplot_path:
        print_dual(f"Gnuplot script : {gplot_path}", f_out)
    print_dual(f"References     : {os.path.join(args.output_dir, BIB_FILE)}", f_out)
    if report_path:
        print_dual(f"Report         : {report_path}", f_out)

    if f_out:
        f_out.close()

    # --view runs last, after the report is fully printed/closed, so a
    # blocking matplotlib window never delays or hides it.
    if args.view:
        import matplotlib.pyplot as plt
        if em_voigt is not None:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
            ax1.bar([f"m*_{l}" for l in VOIGT_LABELS], em_voigt, color="#1f77b4")
            ax1.set_ylabel("m* (m0)")
            ax1.set_title("Effective mass (per-axis Voigt)")
            ax1.grid(True, axis='y')
        else:
            fig, ax2 = plt.subplots(figsize=(6, 5))
        ax2.bar([f"v_{a}" for a in "xyz"], v_kms, color="#d62728")
        ax2.set_ylabel("v (km/s)")
        ax2.set_title("Band velocity")
        ax2.grid(True, axis='y')
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    main()
