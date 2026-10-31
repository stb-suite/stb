#!/usr/bin/env python
"""
SIESTA ToolBox - HAM_WANTIBEXOS Interface
Convert SIESTA Hamiltonian to Wantibexos code
"""

import os
import sys
import re
import io
import time
import argparse
import textwrap
from datetime import datetime
from time import sleep
import numpy as np
from typing import Tuple
from stb.core import siesta_log
from stb.core.cli import COLORS, color_text, show_intro, print_dual, print_table, Sections
from stb.core.kspace import detect_vacuum_axes, dimensionality_label
from stb.core.deps import require_sisl
sisl = require_sisl()

from stb import __version__ as VERSION

REPORT_FILE = "stb_wtb_report.txt"
VACUUM_GAP = 4.0              # Ang of empty space that make an axis non-periodic (2D/1D/0D)
BYTES_PER_LINE = 84           # a tb-*.ham matrix-element line, for the --dry-run size estimate
NEAR_EF = 5.0                 # eV: the energy window the read-back check is judged on
NEAR_EF_TOL = 1e-3            # eV: largest eigenvalue change accepted inside that window

def parse_arguments() -> argparse.Namespace:
    """Handle command-line arguments with argparse"""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''
        Convert a SIESTA Hamiltonian (<label>.HSX, written with SaveHS true) to the tight-binding
        files of WanTiBEXOS: tb-<NP|SP|SOC>.ham (Fermi level, lattice, then every H and S element
        between the orbitals of the cell and of each neighbouring cell) and basis_set-<..>.basis
        (one line per orbital: species, position, l, m, spin). The report describes the system,
        the basis, the Hamiltonian's range and the Fermi level, and checks the written .ham by
        rebuilding H(k) and S(k) from it and comparing the eigenvalues with sisl's.
        '''),
        epilog=textwrap.dedent('''
        Example usage:
          stb-siesta2wtb -i calc.fdf -o calc.out           # Fermi level from the SIESTA log
          stb-siesta2wtb -i calc.fdf -f -4.02              # Fermi level given by hand (eV)
          stb-siesta2wtb -i calc.fdf -o calc.out --dry-run # report only, write nothing
        '''))

    parser.add_argument('-i', '--input', required=True,
                        help='Input FDF file from SIESTA calculation')
    parser.add_argument('-o', '--output', default=None,
                        help='SIESTA output file for Fermi energy extraction')
    parser.add_argument('-f', '--fermi-level', type=float, default=None,
                        help='Manual Fermi energy setting (optional)')
    parser.add_argument('--no-check', dest='check', action='store_false',
                        help='Skip the read-back check of the written tb-*.ham (rebuilding H(k), S(k)\n'
                             'from the file and comparing eigenvalues with sisl).')
    parser.add_argument('--dry-run', action='store_true',
                        help='Report only: read the Hamiltonian and say what would be written; write nothing.')
    parser.add_argument('--save-report', action='store_true',
                        help=f'Also write the report to {REPORT_FILE} (current folder). Off by default.')
    parser.add_argument('-q', '--quiet', action='store_false',
                        help='Suppress terminal animations and colors')
    parser.add_argument("-v", "--version", action="version",
                        version=f"stb-siesta2wtb {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    return parser.parse_args()

def get_system_label(fdf_path: str) -> str:
    """Extract SystemLabel from SIESTA .fdf file with fallback to 'siesta'"""
    # NOTE: This function isn't strictly necessary if sisl can
    # find the .HSX file from the .fdf, but we'll keep it.
    default_label = "siesta"
    try:
        with open(fdf_path, 'r') as f:
            for line in f:
                # Strip comments and whitespace
                line = re.sub(r'#.*', '', line).strip()
                if not line:
                    continue
                
                # Split into key-value pairs
                key_value = re.split(r'\s*=\s*|\s+', line, maxsplit=1)
                if key_value[0].lower() == 'systemlabel':
                    return key_value[1].split()[0] if len(key_value) > 1 else default_label
    except Exception as e:
        pass
    return default_label

from typing import Tuple

def spin_mapping(spin_obj) -> Tuple[str, int]:
    """Map SIESTA spin types to Wantibexos conventions (code, multiplicity)."""
    # This function was correct. It looks for the keyword
    # in the full string (e.g., 'unpolarized' in '<Spin.unpolarized>')
    spin_str = str(spin_obj).lower()
    mapping = {
        'unpolarized': ('NP', 1),
        'polarized': ('SP', 2),
        'non-colinear': ('SOC', 2),
        'spin-orbit': ('SOC', 2)
    }
    for key in mapping:
        if key in spin_str:
            return mapping[key]
    raise ValueError(f"Unsupported spin type in: '{spin_str}'")



def format_float(value: float) -> str:
    """Format floating point numbers for TB files"""
    return f"{value:12.6f}"


def write_basis(hamiltonian, spin_suffix: str, spin_factor: int) -> None:
    """Write basis set information to file"""
    filename = f"basis_set-{spin_suffix}.basis"
    with open(filename, 'w') as f:
        f.write("# bindex aspecie ax ay az l m spin\n")
        
        for spin in range(spin_factor):
            io = 0
            for ia, atom in enumerate(hamiltonian.geometry.atoms):
                xyz = hamiltonian.geometry.xyz[ia, :]
                for orbital in atom:
                    spin_val = 1 if spin == 0 else -1
                    if spin_suffix == 'SOC':
                        spin_val = 1 - 2*spin  # 1 for up, -1 for down
                    f.write(f"{io+1 + spin*hamiltonian.no:6d} {atom.tag:6} "
                            f"{format_float(xyz[0])} {format_float(xyz[1])} "
                            f"{format_float(xyz[2])} {orbital.l:2} "
                            f"{orbital.m:2} {spin_val:3}\n")
                    io += 1

def write_hamiltonian(hamiltonian, spin_suffix: str, fermi: float) -> None:
    """Write Hamiltonian data in TB format"""
    filename = f"tb-{spin_suffix}.ham"
    nbasis = hamiltonian.no
    ncell = np.prod(hamiltonian.nsc)

    with open(filename, 'w') as f:
        # Header section
        f.write(f"{spin_suffix}\n")
        f.write(f"{format_float(0.0)}\n")  # Scissors operator
        f.write(f"{format_float(fermi)}\n")
        
        # Cell vectors
        for vec in hamiltonian.cell:
            f.write(" ".join([format_float(v) for v in vec]) + "\n")
        
        # Dimensions
        f.write(f"{hamiltonian.no * (2 if spin_suffix != 'NP' else 1)}\n")
        f.write(f"{ncell}\n")
        f.write("# rcell_x rcell_y rcell_z i j ReH ImH S\n")
        
        # Hamiltonian data handling based on spin type
        if spin_suffix == 'NP':
            # Unpolarized case
            for icell in range(ncell):
                for i in range(nbasis):
                    for j in range(nbasis):
                        sc_index = hamiltonian.geometry.o2sc(j + icell*nbasis)
                        H = hamiltonian[i, j + icell*nbasis][0]
                        S = hamiltonian.S[i, j + icell*nbasis]
                        
                        line = (f"{' '.join(format_float(v) for v in sc_index)} "
                                f"{i+1} {j+1} {format_float(H)} "
                                f"{format_float(0.0)} {format_float(S)}\n")
                        f.write(line)
        
        elif spin_suffix == 'SP':
            # Spin polarized case
            S = np.zeros((ncell+1, 2*nbasis+2, 2*nbasis+2))
            reH = np.zeros((ncell+1, 2*nbasis+2, 2*nbasis+2))
            imH = np.zeros((ncell+1, 2*nbasis+2, 2*nbasis+2)) # Although not used here, it's good for consistency
            a = np.zeros((ncell+1, 2*nbasis+2))
            b = np.zeros((ncell+1, 2*nbasis+2))
            c = np.zeros((ncell+1, 2*nbasis+2))

            # Populate matrices
            for icell in range(ncell):
                for i in range(nbasis):
                    for j in range(nbasis):
                        # Up-up block
                        reH[icell+1, i+1, j+1] = hamiltonian[i, j + icell*nbasis][0]
                        S[icell+1, i+1, j+1] = hamiltonian.S[i, j + icell*nbasis]
                        
                        # Down-down block
                        reH[icell+1, i+1+nbasis, j+1+nbasis] = hamiltonian[i, j + icell*nbasis][1]
                        S[icell+1, i+1+nbasis, j+1+nbasis] = hamiltonian.S[i, j + icell*nbasis]
                        
                        # Cell coordinates
                        sc = hamiltonian.geometry.o2sc(j + icell*nbasis)
                        a[icell+1, j+1] = sc[0]
                        b[icell+1, j+1] = sc[1]
                        c[icell+1, j+1] = sc[2]
                        a[icell+1, j+1+nbasis] = sc[0]
                        b[icell+1, j+1+nbasis] = sc[1]
                        c[icell+1, j+1+nbasis] = sc[2]

            # Write to file
            for icell in range(1, ncell+1):
                for i in range(1, 2*nbasis+1):
                    for j in range(1, 2*nbasis+1):
                        line = (f"{format_float(a[icell,j])} {format_float(b[icell,j])} "
                                f"{format_float(c[icell,j])} {i} {j} "
                                f"{format_float(reH[icell,i,j])} {format_float(imH[icell,i,j])} " # ImH will be 0.0
                                f"{format_float(S[icell,i,j])}\n")
                        f.write(line)
        
        elif spin_suffix == 'SOC':
            # Non-colinear and spin-orbit cases
            S = np.zeros((ncell+1, 2*nbasis+2, 2*nbasis+2))
            reH = np.zeros((ncell+1, 2*nbasis+2, 2*nbasis+2))
            imH = np.zeros((ncell+1, 2*nbasis+2, 2*nbasis+2))
            a = np.zeros((ncell+1, 2*nbasis+2))
            b = np.zeros((ncell+1, 2*nbasis+2))
            c = np.zeros((ncell+1, 2*nbasis+2))

            for icell in range(ncell):
                for i in range(nbasis):
                    for j in range(nbasis):
                        sc = hamiltonian.geometry.o2sc(j + icell*nbasis)
                        idx = icell+1
                        
                        # Common cell coordinates
                        a[idx, j+1] = sc[0]
                        b[idx, j+1] = sc[1]
                        c[idx, j+1] = sc[2]
                        a[idx, j+1+nbasis] = sc[0]
                        b[idx, j+1+nbasis] = sc[1]
                        c[idx, j+1+nbasis] = sc[2]

                        # Matrix elements
                        H = hamiltonian[i, j + icell*nbasis]
                        S_val = hamiltonian.S[i, j + icell*nbasis]
                        
                        # Up-up block
                        reH[idx, i+1, j+1] = H[0]
                        S[idx, i+1, j+1] = S_val
                        
                        # Down-down block
                        reH[idx, i+1+nbasis, j+1+nbasis] = H[1]
                        S[idx, i+1+nbasis, j+1+nbasis] = S_val
                        
                        # CORRECTION APPLIED HERE
                        # siesta2wtb.py shows that H[4] and H[5] are the imaginary parts
                        # of the diagonal blocks (up-up, dn-dn) for SOC.
                        # Non-colinear (NC) doesn't have H[4] or H[5].
                        
                        # sisl returns 4 elements for NC (H[0..3])
                        # and 8 elements for SOC (H[0..7])
                        if len(H) > 4:  # Spin-orbit (SOC) case
                            
                            # Up-up block (imaginary)
                            imH[idx, i+1, j+1] = H[4] # <- THIS LINE WAS ADDED/CORRECTED

                            # Up-down block (real and imaginary)
                            reH[idx, i+1, j+1+nbasis] = H[2]
                            imH[idx, i+1, j+1+nbasis] = H[3]
                            
                            # Down-up block (real and imaginary)
                            reH[idx, i+1+nbasis, j+1] = H[6]
                            imH[idx, i+1+nbasis, j+1] = H[7]
                            
                            # Down-down block (imaginary)
                            imH[idx, i+1+nbasis, j+1+nbasis] = H[5]

                        elif len(H) > 2:  # Non-colinear (NC) case
                            # (Has 4 elements, H[0], H[1], H[2], H[3])
                            
                            # Up-down block (real)
                            reH[idx, i+1, j+1+nbasis] = H[2]
                            
                            # Down-up block (real)
                            reH[idx, i+1+nbasis, j+1] = H[3]

            # Write to file
            for icell in range(1, ncell+1):
                for i in range(1, 2*nbasis+1):
                    for j in range(1, 2*nbasis+1):
                        line = (f"{format_float(a[icell,j])} {format_float(b[icell,j])} "
                                f"{format_float(c[icell,j])} {i} {j} "
                                f"{format_float(reH[icell,i,j])} "
                                f"{format_float(imH[icell,i,j])} "
                                f"{format_float(S[icell,i,j])}\n")
                        f.write(line)



# --------------------------------------------------------------------------------------------- report helpers --
def human_size(n):
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024.0


def cell_parameters(cell):
    lengths = np.linalg.norm(cell, axis=1)

    def angle(i, j):
        return np.degrees(np.arccos(np.clip(np.dot(cell[i], cell[j]) / (lengths[i] * lengths[j]), -1, 1)))
    return lengths, (angle(1, 2), angle(0, 2), angle(0, 1)), abs(np.linalg.det(cell))


def read_log_facts(path):
    """What the SIESTA log says that matters here; every field None when absent."""
    facts = dict(version=None, spin=None, electrons=None, finished=False)
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                s = line.strip()
                if s.startswith("Version") and ":" in s and facts["version"] is None:
                    facts["version"] = s.split(":", 1)[1].strip()
                elif "redata: Spin configuration" in s and "=" in s:
                    facts["spin"] = s.split("=", 1)[1].strip()
                elif s.startswith("Total number of electrons:"):
                    try:
                        facts["electrons"] = float(s.split(":", 1)[1].split()[0])
                    except (IndexError, ValueError):
                        pass
                elif "End of run" in s or s.startswith("Job completed"):
                    facts["finished"] = True
    except OSError:
        return None
    facts["fermi"] = siesta_log.get_fermi_energy(path)
    facts["spin_moment"] = siesta_log.get_spin_moment(path)
    return facts


def shell_summary(atom):
    """'2s(2z) 2p(2z) 3d(pol)': the radial shells of one species' basis, from sisl's AtomicOrbitals."""
    shells = {}
    for orb in atom.orbitals:
        n, l = getattr(orb, "n", None), getattr(orb, "l", None)
        if l is None:
            return f"{atom.no} orbitals"
        key = (n, int(l), bool(getattr(orb, "P", False)))
        shells.setdefault(key, set()).add(int(getattr(orb, "zeta", 1)))
    parts = []
    for (n, l, pol), zetas in sorted(shells.items(), key=lambda kv: (kv[0][2], kv[0][0] or 0, kv[0][1])):
        name = f"{n if n is not None else ''}{'spdfg'[l]}"
        parts.append(f"{name}(pol)" if pol else f"{name}({len(zetas)}z)")
    return " ".join(parts)


def read_ham_file(path):
    """Header and matrix-element columns of a tb-*.ham as written by write_hamiltonian."""
    with open(path) as f:
        head = [next(f) for _ in range(9)]
    header = dict(code=head[0].strip(), scissor=float(head[1]), fermi=float(head[2]),
                  cell=np.array([[float(x) for x in head[i].split()] for i in (3, 4, 5)]),
                  dim=int(head[6]), ncell=int(head[7]))
    try:
        import pandas as pd
        data = pd.read_csv(path, skiprows=9, sep=r"\s+", header=None, engine="c").to_numpy()
    except ImportError:
        data = np.loadtxt(path, skiprows=9)
    return header, data


def eigen_from_file(data, dim, cell, kfrac):
    """Eigenvalues of H(k) c = E S(k) c, with H(k) = sum_R H(R) exp(i k.R) built from the file's lines."""
    import scipy.linalg as sl
    k = np.asarray(kfrac, dtype=float) @ (2 * np.pi * np.linalg.inv(cell).T)
    phase = np.exp(1j * (data[:, :3] @ k))
    i = data[:, 3].astype(int) - 1
    j = data[:, 4].astype(int) - 1
    hk = np.zeros((dim, dim), complex)
    sk = np.zeros((dim, dim), complex)
    np.add.at(hk, (i, j), (data[:, 5] + 1j * data[:, 6]) * phase)
    np.add.at(sk, (i, j), data[:, 7] * phase)
    herm = float(max(np.abs(hk - hk.conj().T).max(), np.abs(sk - sk.conj().T).max()))
    e = sl.eigh((hk + hk.conj().T) / 2, (sk + sk.conj().T) / 2, eigvals_only=True)
    return np.sort(e), herm


def eigen_from_sisl(ham, kfrac):
    if ham.spin.is_polarized:
        return np.sort(np.concatenate([ham.eigh(k=kfrac, spin=0), ham.eigh(k=kfrac, spin=1)]))
    return np.sort(ham.eigh(k=kfrac))


def main():
    """Main workflow controller"""

    args = parse_arguments()

    # Hides the intro if --quiet or --no-intro is used
    # Note the quiet logic is inverted (action='store_false')
    # So args.quiet will be 'False' if -q is passed (same for args.intro/--no-intro)
    show_banner = args.quiet and args.intro

    if show_banner:
        show_intro([
            "Siesta ToolBox Suite",
            "SIESTA-to-Wantibexos Hamiltonian Interface",
            f"Version {VERSION} | University of Brasilia - 2026",
            ""
        ])
        print("\n" + color_text("SIESTA-WANTIBEXOS Interface:", 'bold'))
        print("-"*60)

    buf = io.StringIO() if args.save_report else None
    sec = Sections(buf)
    warnings = []

    def out(text=""):
        print_dual(text, buf)

    def warn(message):
        warnings.append(message)

    def finish(code=0):
        if buf is not None:
            with open(REPORT_FILE, "w") as fh:
                fh.write(buf.getvalue())
        if code:
            sys.exit(code)

    def fail(message):
        out(color_text(f"[FAIL] {message}", 'red'))
        finish(1)

    out(color_text("===== STB-SIESTA2WTB REPORT =====", 'magenta'))
    folder = os.path.dirname(os.path.abspath(args.input))
    system_label = get_system_label(args.input)

    # [0] --------------------------------------------------------------------------------------------------------
    sec("RUN METADATA")
    out(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}")
    out(f"stb / sisl        : {VERSION} / {sisl.__version__}")
    out(f"Input             : {args.input} (SystemLabel {system_label})")
    if args.fermi_level is not None:
        out(f"Fermi level       : {args.fermi_level} eV, given with -f")
    elif args.output:
        out(f"Fermi level       : read from {args.output} (-o)")
    else:
        out("Fermi level       : none given (-o or -f) -- 0.0 eV is written")
    out("Mode              : " + ("dry run (nothing written)" if args.dry_run else
                                  "write" + (" + read-back check" if args.check else " (no read-back check)")))

    # [1] --------------------------------------------------------------------------------------------------------
    sec("SIESTA RUN")
    if not os.path.isfile(args.input):
        fail(f"Input file '{args.input}' not found.")
    hsx = os.path.join(folder, f"{system_label}.HSX")
    if os.path.isfile(hsx):
        st = os.stat(hsx)
        out(f"Hamiltonian       : {os.path.relpath(hsx)} ({human_size(st.st_size)}, written "
            f"{datetime.fromtimestamp(st.st_mtime):%Y-%m-%d %H:%M})")
    else:
        out(f"Hamiltonian       : no {system_label}.HSX next to the input -- sisl will look for another H source")
    try:
        with open(args.input, errors="ignore") as fh:
            savehs = re.search(r"^\s*SaveHS\s+(\.?true\.?|t|yes)\b", fh.read(), re.I | re.M)
        out(f"SaveHS in input   : {'yes' if savehs else 'not in this file (maybe in an %include, or a default)'}")
    except OSError:
        pass
    log_path = args.output or siesta_log.find_out_file(folder, system_label)
    log = read_log_facts(log_path) if log_path and os.path.isfile(log_path) else None
    if args.output and log is None:
        fail(f"SIESTA output '{args.output}' (-o) not found.")
    if log is None:
        out("Log               : none found (-o, <label>.out or the only *.out) -- no cross-checks with SIESTA")
    else:
        out(f"Log               : {log_path}" + (f" (SIESTA {log['version']})" if log["version"] else "")
            + ("" if args.output else "  (found, not given with -o: used for the report only)"))
        out(f"Finished          : {'yes' if log['finished'] else 'NO'}")
        if not log["finished"]:
            warn(f"{log_path} did not finish: the Hamiltonian may be from an unconverged SCF step.")
        if log["electrons"] is not None:
            out(f"Electrons         : {log['electrons']:.4f} (valence)")
        if log["spin"]:
            out(f"Spin (log)        : {log['spin']}")
        if log["fermi"] is not None:
            out(f"Fermi level (log) : {log['fermi']:.6f} eV")
        if log["spin_moment"] is not None:
            out(f"Spin moment (log) : {log['spin_moment']:.5f} mu_B")

    t0 = time.time()
    try:
        geom = sisl.get_sile(args.input).read_geometry()
        ham = sisl.get_sile(args.input).read_hamiltonian(geometry=geom)
    except Exception as e:
        fail(f"Could not read the Hamiltonian through '{args.input}': {e} (the run needs SaveHS true).")
    if ham is None:
        fail(f"No Hamiltonian found for '{args.input}' (no {system_label}.HSX): rerun SIESTA with SaveHS true.")
    out(f"Read              : geometry + H + S in {time.time() - t0:.1f} s")
    # sisl reads E_F from an HSX of version >= 1 (SIESTA >= 4.1) and shifts H by -E_F S: the matrix elements that
    # go into the file are then referred to the run's Fermi level, whatever Fermi level the header gets.
    h_zero = None
    try:
        if not os.path.isfile(hsx):
            raise FileNotFoundError(hsx)
        hsx_sile = sisl.get_sile(hsx)
        if (getattr(hsx_sile, "version", 0) or 0) >= 1:
            h_zero = float(hsx_sile.read_fermi_level())
            out(f"Energy zero of H  : the run's Fermi level -- {system_label}.HSX (version {hsx_sile.version}) stores "
                f"E_F = {h_zero:.6f} eV and sisl shifts H by -E_F S on reading")
        else:
            out(f"Energy zero of H  : SIESTA's absolute scale ({system_label}.HSX version 0 has no E_F)")
    except Exception:
        out("Energy zero of H  : unknown (the .HSX could not be inspected)")

    # [2] --------------------------------------------------------------------------------------------------------
    sec("SYSTEM")
    g = ham.geometry
    lengths, angles, volume = cell_parameters(np.asarray(ham.cell))
    symbols = [a.symbol for a in g.atoms]
    counts = {}
    for s_ in symbols:
        counts[s_] = counts.get(s_, 0) + 1
    out(f"Atoms             : {len(symbols)} ({''.join(f'{k}{v if v > 1 else str()}' for k, v in counts.items())})")
    out(f"Cell lengths      : a = {lengths[0]:.4f}, b = {lengths[1]:.4f}, c = {lengths[2]:.4f} Ang")
    out(f"Cell angles       : alpha = {angles[0]:.2f}, beta = {angles[1]:.2f}, gamma = {angles[2]:.2f} deg")
    out(f"Cell volume       : {volume:.4f} Ang^3")
    vac = detect_vacuum_axes(g.xyz @ np.linalg.inv(ham.cell), np.asarray(ham.cell), VACUUM_GAP)
    dim_label = dimensionality_label(vac)
    out(f"Dimensionality    : {dim_label}"
        + (f" -- vacuum along {', '.join(n for n, v in zip('abc', vac) if v)}" if any(vac) else ""))
    rows = []
    for sp in g.atoms.atom:
        n_at = sum(1 for a in g.atoms if a.tag == sp.tag)
        rmax = max((o.R for o in sp.orbitals), default=-1)
        rows.append(([sp.tag, str(sp.Z), str(n_at), str(sp.no), shell_summary(sp),
                      f"{rmax:.2f}" if rmax > 0 else "-"], None))
    print_table(["species", "Z", "atoms", "orbitals", "basis shells", "r_c (Ang)"], rows, buf)
    out(f"Orbitals          : {ham.no} in the cell")

    # [3] --------------------------------------------------------------------------------------------------------
    sec("HAMILTONIAN")
    try:
        spin_suffix, spin_factor = spin_mapping(str(ham.spin))
    except ValueError as e:
        fail(str(e))
    ncell = int(np.prod(ham.nsc))
    dim = ham.no * (2 if spin_suffix != 'NP' else 1)
    out(f"Spin (sisl)       : {ham.spin} -> WanTiBEXOS type {spin_suffix}, {dim} x {dim} matrices")
    out(f"Orthogonal basis  : {'yes' if ham.orthogonal else 'no -- the overlap S goes into the file with H'}")
    out(f"Neighbour cells   : nsc = {ham.nsc[0]} x {ham.nsc[1]} x {ham.nsc[2]} = {ncell} cells (all R with an "
        "orbital overlapping one of the cell)")
    dense = ham.no * ham.no * ncell
    out(f"Stored elements   : {ham.nnz} of {dense} ({100 * ham.nnz / dense:.1f} %) per spin component")
    try:
        hc = ham.tocsr(0).tocoo()
        dist = np.linalg.norm(g.axyz(np.arange(g.na_s))[g.o2a(hc.col)] - g.xyz[g.o2a(hc.row)], axis=1)
        out(f"Longest hopping   : {dist.max():.3f} Ang (orbital centres; the range of the basis)")
        diag = ham.tocsr(0).diagonal()[:ham.no]
        out(f"On-site energies  : {diag.min():.3f} to {diag.max():.3f} eV"
            + (" (from E_F, the zero of H)" if h_zero is not None else " (SIESTA's absolute scale)"))
    except Exception as e:
        warn(f"could not analyse the Hamiltonian's range: {e}")
    try:
        s_gamma = np.linalg.eigvalsh(ham.Sk(k=[0, 0, 0]).toarray())
        out(f"S(Gamma) spectrum : {s_gamma.min():.2e} to {s_gamma.max():.3f} (smallest eigenvalue: how close the "
            "basis is to linear dependence)")
        if s_gamma.min() < 1e-3:
            warn(f"the overlap S(Gamma) has an eigenvalue of {s_gamma.min():.1e}: the basis is nearly linearly "
                 "dependent, and the 6 decimals of the file perturb the highest bands (the read-back check measures by how much).")
    except Exception:
        pass
    if ham.spin.is_polarized and spin_suffix == 'SP':
        out("Layout            : up block (orbitals 1..N), then down block (N+1..2N); no up-down coupling")
    elif spin_suffix == 'SOC':
        out("Layout            : 2N x 2N with up-down blocks (non-collinear or spin-orbit)")

    # [4] --------------------------------------------------------------------------------------------------------
    sec("FERMI LEVEL AND BAND EDGES")
    fermi = 0.0
    if args.fermi_level is not None:
        fermi = args.fermi_level
        out(f"Fermi level       : {fermi:.6f} eV (given with -f)")
        if log and log["fermi"] is not None and abs(log["fermi"] - fermi) > 1e-4:
            out(f"                    the log says {log['fermi']:.6f} eV ({fermi - log['fermi']:+.4f} eV difference)")
    elif args.output is not None:
        fermi = siesta_log.get_fermi_energy(args.output)
        if fermi is None:
            fail(f"Could not find the Fermi energy in '{args.output}'.")
        out(f"Fermi level       : {fermi:.6f} eV (from {args.output})")
    else:
        out("Fermi level       : 0.0 eV (no -o, no -f)")
        if h_zero is None:
            hint = (f" -- {log_path} has E_F = {log['fermi']:.4f} eV: pass -o {log_path}" if log and log["fermi"]
                    is not None else "")
            warn("no Fermi level given (-o or -f) and H is on SIESTA's absolute scale: the header's 0.0 eV is not "
                 "the Fermi level" + hint + ".")
    h_ref = 0.0 if h_zero is not None else fermi          # the energy that is E_F in the H written to the file
    if h_zero is not None:
        if abs(fermi) > 1e-6:
            out(f"Reference         : H in the file is already referred to E_F (0 eV = {h_zero:.4f} eV absolute); the "
                f"header adds E_F = {fermi:.4f} eV")
            warn(f"the matrix elements are already referred to the Fermi level (sisl shifted them by -E_F S, E_F = "
                 f"{h_zero:.4f} eV from {system_label}.HSX), and the header also carries E_F = {fermi:.4f} eV: a reader "
                 f"that subtracts the header's E_F from H moves every band by {-fermi:+.4f} eV (with -f 0 the header "
                 "would match H).")
        else:
            out("Reference         : header E_F = 0.0 matches H, which is already referred to E_F")
    eig = os.path.join(folder, f"{system_label}.EIG")
    gap = siesta_log.read_eig_gap(eig) if os.path.isfile(eig) else {}
    if gap:
        if gap["metallic"] or gap["gap"] < 0.05:
            out(f"SIESTA .EIG       : metallic on the SCF k-mesh ({gap['n_k']} k-points; gap {gap['gap']:.3f} eV, "
                "below the k-mesh/smearing resolution)")
        else:
            out(f"SIESTA .EIG       : VBM {gap['vbm']:.4f}, CBM {gap['cbm']:.4f} eV -> gap {gap['gap']:.4f} eV "
                f"({'direct' if gap['direct'] else 'indirect'}; smallest direct {gap['direct_min']:.4f} eV; "
                f"{gap['n_k']} k-points)")
            if h_zero is None and not (gap["vbm"] - 1e-3 <= fermi <= gap["cbm"] + 1e-3) and (
                    args.fermi_level is not None or args.output):
                warn(f"the Fermi level written ({fermi:.4f} eV) is not inside the gap [{gap['vbm']:.4f}, "
                     f"{gap['cbm']:.4f}] eV of the SIESTA run.")
    else:
        out("SIESTA .EIG       : not found -- no band edges from the SCF k-mesh")
    occupied = None
    if log and log["electrons"] is not None:
        occupied = log["electrons"] / (2.0 if spin_suffix == 'NP' else 1.0)
        out(f"Occupied states   : {occupied:g} of {dim} per k ({log['electrons']:g} electrons"
            + (", 2 per state)" if spin_suffix == 'NP' else ", 1 per spin-orbital)"))

    # [5] --------------------------------------------------------------------------------------------------------
    sec("OUTPUT FILES")
    basis_file = f"basis_set-{spin_suffix}.basis"
    ham_file = f"tb-{spin_suffix}.ham"
    n_lines = ncell * dim * dim
    if args.dry_run:
        out(f"Would write       : {basis_file} ({dim} orbital lines)")
        out(f"                    {ham_file} ({n_lines} matrix-element lines, ~{human_size(n_lines * BYTES_PER_LINE)})")
    else:
        t0 = time.time()
        write_basis(ham, spin_suffix, spin_factor)
        t_basis = time.time() - t0
        t0 = time.time()
        write_hamiltonian(ham, spin_suffix, fermi)
        t_ham = time.time() - t0
        out(f"Basis             : {basis_file} ({dim} orbitals: index, species, x y z (Ang), l, m, spin; "
            f"{human_size(os.path.getsize(basis_file))}, {t_basis:.1f} s)")
        out(f"Hamiltonian       : {ham_file} ({human_size(os.path.getsize(ham_file))}, {t_ham:.1f} s)")
        out(f"Header            : type {spin_suffix}, scissor 0.0, E_F {fermi:.6f} eV, 3 lattice vectors (Ang), "
            f"{dim} orbitals, {ncell} cells")
        out(f"Body              : {n_lines} lines 'R_x R_y R_z i j ReH ImH S' (R in Ang, H in eV), every i, j of "
            "every cell, zeros included")

    # [6] --------------------------------------------------------------------------------------------------------
    sec("READ-BACK CHECK")
    if args.dry_run or not args.check:
        out("Skipped           : " + ("--dry-run" if args.dry_run else "--no-check"))
    else:
        try:
            t0 = time.time()
            header, data = read_ham_file(ham_file)
            problems = []
            if header["code"] != spin_suffix:
                problems.append(f"type {header['code']}")
            if header["dim"] != dim or header["ncell"] != ncell:
                problems.append(f"dimensions {header['dim']} x {header['ncell']}")
            if data.shape[0] != n_lines:
                problems.append(f"{data.shape[0]} lines instead of {n_lines}")
            if abs(header["fermi"] - fermi) > 1e-5:
                problems.append(f"E_F {header['fermi']}")
            if np.abs(header["cell"] - np.asarray(ham.cell)).max() > 1e-5:
                problems.append("lattice vectors")
            nonzero = int(np.count_nonzero((data[:, 5] != 0) | (data[:, 6] != 0) | (data[:, 7] != 0)))
            out(f"File read         : {data.shape[0]} lines in {time.time() - t0:.1f} s; header "
                + ("as written" if not problems else "DIFFERS: " + ", ".join(problems)))
            out(f"Non-zero lines    : {nonzero} ({100 * nonzero / max(1, data.shape[0]):.1f} %)")
            kpts = [[0, 0, 0], [0.5, 0, 0], [1 / 3, 1 / 3, 0], [0.13, 0.27, 0.31]]
            kpts = [[0.0 if vac[ax] else kf[ax] for ax in range(3)] for kf in kpts]
            kpts = [k for i_, k in enumerate(kpts) if k not in kpts[:i_]]
            rows, worst_near, worst_all, herm_max = [], 0.0, 0.0, 0.0
            for kf in kpts:
                e_file, herm = eigen_from_file(data, dim, header["cell"], kf)
                e_ref = eigen_from_sisl(ham, kf)
                de = np.abs(e_file - e_ref)
                near = np.abs(e_ref - h_ref) < NEAR_EF
                below = int(np.count_nonzero(e_file < h_ref))
                d_near = float(de[near].max()) if near.any() else 0.0
                worst_near, worst_all = max(worst_near, d_near), max(worst_all, float(de.max()))
                herm_max = max(herm_max, herm)
                rows.append(([f"({kf[0]:.3f}, {kf[1]:.3f}, {kf[2]:.3f})", str(below), f"{d_near * 1000:.3f}",
                              f"{de.max() * 1000:.1f}", f"{e_ref[np.argmax(de)] - h_ref:+.1f}"], None))
            out(f"Eigenvalues       : H(k) c = E S(k) c from the file vs sisl on the SIESTA H, S (k in reduced units"
                + (", 0 along vacuum axes" if any(vac) else "") + ")")
            print_table(["k", "states < E_F", f"max |dE| within {NEAR_EF:g} eV of E_F (meV)", "max |dE| all (meV)",
                         "worst band at E - E_F (eV)"], rows, buf)
            out(f"Hermiticity       : max |H - H^+|, |S - S^+| from the file = {herm_max:.1e}")
            ok = worst_near < NEAR_EF_TOL and not problems
            out("Verdict           : " + (f"ok -- within {NEAR_EF:g} eV of E_F the file reproduces the SIESTA bands to "
                                         f"{worst_near * 1000:.3f} meV" if ok else "FAILED"))
            if worst_all > 10 * max(worst_near, 1e-6) and worst_all > 1e-3:
                out(f"                    the largest changes ({worst_all * 1000:.0f} meV) are in high, unphysical bands: "
                    "the file's 6 decimals perturb S where the basis is nearly linearly dependent")
            if occupied is not None and rows:
                gamma_below = int(rows[0][0][1])
                insulator = bool(gap) and not gap["metallic"] and gap["gap"] >= 0.05
                if insulator and abs(gamma_below - occupied) > 0.5:
                    warn(f"at Gamma {gamma_below} states of the file's H lie below its Fermi level, but the system has "
                         f"{occupied:g} occupied states -- check the Fermi level.")
            if not ok:
                warn("the tb file does not reproduce the SIESTA eigenvalues near E_F "
                     f"(max {worst_near * 1000:.2f} meV)" + (": " + ", ".join(problems) if problems else "") + ".")
        except Exception as e:
            warn(f"could not read {ham_file} back to check it: {e}")

    # [7] --------------------------------------------------------------------------------------------------------
    sec("FOR WANTIBEXOS")
    sysdim = {0: "3D", 1: "2D", 2: "1D", 3: "0D"}[sum(vac)]
    out(f"PARAMS_FILE       : {ham_file}")
    out(f"SYSDIM            : {sysdim} (from the vacuum axes)")
    if occupied is not None:
        out(f"Valence states    : {occupied:g} occupied of {dim} -- NBANDSV can be at most {int(occupied)}, "
            f"NBANDSC at most {dim - int(occupied)}")
    out("Energies          : eV; " + (f"H referred to E_F (0 eV = {h_zero:.4f} eV absolute), header E_F {fermi:.4f} eV"
                                      if h_zero is not None else f"H on SIESTA's absolute scale, header E_F {fermi:.4f} eV"))

    # [8] --------------------------------------------------------------------------------------------------------
    sec("SUMMARY & FILES")
    if args.dry_run:
        out(f"Dry run           : nothing written ({ham_file} would have {n_lines} lines)")
    else:
        out(f"Written           : {ham_file}, {basis_file}")
    out(f"System            : {''.join(f'{k}{v if v > 1 else str()}' for k, v in counts.items())}, {spin_suffix}, "
        f"{dim} orbitals, {ncell} cells, E_F {fermi:.4f} eV")
    if buf is not None:
        out(f"Report            : {REPORT_FILE}")

    sec("WARNINGS")
    if warnings:
        for n_, w in enumerate(warnings, 1):
            out(color_text(f"  {n_}. {w}", 'yellow'))
    else:
        out("  none")
    finish(0)

    print("\n"+"-"*60)
    if show_banner:
        print(color_text("Saving Hamiltonians before they collapse their own wavefunctions.\n\n", 'bold'))


if __name__ == "__main__":
    main()
