"""Single source of truth for parsing SIESTA calculation logs (.out files).

Consolidates parsers that used to live independently in elastic_analysis.py,
strain_analysis.py, workfunction.py, wantibexos.py and cohesive_analysis.py:
Fermi energy, cell height (Lz), stress tensor (two flavors -- see note below)
and free energy. get_electric_dipole() was added directly here (ir_analysis.py
is its only consumer so far) rather than via an extract-on-second-use move,
since it's exactly the kind of "one labeled quantity from a .out log" parser
this module already exists to hold.

Error contract: every function here returns None (or 1.0 for
get_cell_height, a multiplicative normalization factor rather than a
validity check) instead of raising, on file-not-found, pattern-not-found, or
any parse error. This matches how every existing caller already uses these
values -- scanning many strain_* folders while tolerating incomplete/missing
calculations, or doing an explicit `if value is None` check before deciding
whether to abort. Callers are responsible for checking the result.

Note on stress tensor: get_stress_tensor() (matrix block with a Voigt-line
fallback, eV/Angstrom^3) and get_stress_voigt_kbar() (Voigt line only, raw
kBar) are deliberately NOT unified into one function. Their two callers
(elastic_analysis.py and strain_analysis.py) each already do downstream
strain/stress math expecting a specific shape and unit; forcing one
representation would mean rewriting that math too, which is out of scope for
a parser consolidation.

Note on check_scf_and_force()/report_quality_diagnostics(): the "returns
None instead of raising" contract above describes the pure PARSING
functions (get_fermi_energy, get_free_energy, get_scf_convergence,
get_max_force, ...). These two are reporting helpers built on top of those
parsers -- report_quality_diagnostics in particular has real side effects
(printing to stdout, writing to a report file handle) and is not itself a
parser. Moved here from adsorb_analysis.py once stb-nebAnalysis needed the
identical SCF-convergence/residual-force diagnostic, the same
extract-on-second-use policy as structure_io.py/core/symmetry.py.
"""

from __future__ import annotations

import re

import numpy as np

from stb.core.cli import color_text, print_dual

_STRAIN_FOLDER_RE = re.compile(r"strain_([a-zA-Z0-9]+)_(m?)(\d+\.\d+)")
_EV_A3_TO_GPA = 160.21766  # 1 eV/Angstrom^3 in GPa


def parse_strain_folder_name(folder_name: str) -> tuple[str | None, float | None]:
    """Parses 'strain_xx_1.00' / 'strain_xy_m2.00' into (direction, signed value in %)."""
    match = _STRAIN_FOLDER_RE.search(folder_name)
    if not match:
        return None, None
    direction = match.group(1)
    value = float(match.group(3))
    if match.group(2) == 'm':
        value = -value
    return direction, value


def get_fermi_energy(path: str) -> float | None:
    """Extracts Fermi energy (eV) from a SIESTA .out file.

    Tries, in a single pass over the file: the SIESTA 5.x summary line
    ("... Fermi = ..."), the classic format ("Fermi energy: ... eV"), and the
    per-iteration SCF table (column 'Ef'/'Ef(eV)'). Whichever match appears
    last in the file wins, which naturally favors the final summary over
    intermediate SCF steps.
    """
    E_f = None
    scf_ef_idx = None
    try:
        with open(path, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if not parts:
                    continue

                if "Fermi" in parts and "=" in parts:
                    try:
                        eq_idx = parts.index("=")
                        if eq_idx + 1 < len(parts):
                            E_f = float(parts[eq_idx + 1])
                    except ValueError:
                        pass

                elif "Fermi" in line and "energy:" in line:
                    try:
                        colon_idx = -1
                        for i, p in enumerate(parts):
                            if p.endswith(':'):
                                colon_idx = i
                        if colon_idx != -1 and colon_idx + 1 < len(parts):
                            E_f = float(parts[colon_idx + 1].replace('eV', ''))
                    except ValueError:
                        pass

                if "iscf" in line and ("Ef" in line or "Ef(eV)" in line):
                    for i, p in enumerate(parts):
                        if "Ef" in p:
                            scf_ef_idx = i
                            break
                elif scf_ef_idx is not None and parts[0] == 'scf:':
                    try:
                        data_idx = scf_ef_idx + 1
                        if data_idx < len(parts):
                            E_f = float(parts[data_idx])
                    except ValueError:
                        pass
        return E_f
    except Exception:
        return None


def get_cell_height(path: str) -> float:
    """Norm of the c lattice vector (Angstrom), from the LAST 'outcell:' block in the file.

    Used to normalize 2D properties by cell height. Defaults to 1.0 if not
    found -- a deliberate no-op fallback (this is a multiplicative factor,
    not a validity check). Uses the last occurrence, not the first: SIESTA
    reprints this block on every relaxation step, and the first one is the
    unrelaxed starting geometry.
    """
    z_len = 1.0
    try:
        with open(path, 'r', errors='ignore') as f:
            lines = f.readlines()
        for i in range(len(lines) - 1, -1, -1):
            if "outcell: Unit cell vectors" in lines[i]:
                nums = _parse_float_line(lines[i + 3])
                if nums:
                    z_len = float(np.linalg.norm(nums))
                break
    except Exception:
        pass
    return z_len


def get_outcell(path: str) -> np.ndarray | None:
    """Full 3x3 cell matrix (Angstrom, rows a/b/c), from the LAST
    'outcell:' block in the file -- same block and same last-occurrence
    convention as get_cell_height, generalized to all 3 rows instead of
    just the norm of c. Returns None if the block isn't found or fails to
    parse.
    """
    try:
        with open(path, 'r', errors='ignore') as f:
            lines = f.readlines()
        for i in range(len(lines) - 1, -1, -1):
            if "outcell: Unit cell vectors" in lines[i]:
                rows = [_parse_float_line(lines[i + offset]) for offset in (1, 2, 3)]
                if all(rows):
                    return np.array(rows, dtype=float)
                return None
    except Exception:
        return None
    return None


def get_stress_tensor(path: str) -> np.ndarray | None:
    """Stress tensor (3x3, eV/Angstrom^3) from a SIESTA .out file.

    Tries the full "siesta: Stress tensor (static)" matrix block first (uses
    the LAST such block in the file); falls back to reconstructing the
    matrix from the "Stress tensor Voigt" line (converted from kBar) if the
    matrix block isn't found or fails to parse.
    """
    try:
        with open(path, 'r', errors='ignore') as f:
            lines = f.readlines()

        idx_matrix, idx_voigt = -1, -1
        for i, line in enumerate(lines):
            if "siesta: Stress tensor (static)" in line:
                idx_matrix = i
            if "Stress tensor Voigt" in line:
                idx_voigt = i

        if idx_matrix != -1:
            try:
                rows = []
                offset = 1
                while len(rows) < 3 and offset < 10:
                    nums = _parse_float_line(lines[idx_matrix + offset].strip())
                    if nums:
                        rows.append(nums)
                    offset += 1
                if len(rows) == 3:
                    return np.array(rows)
            except Exception:
                pass

        if idx_voigt != -1:
            try:
                parts = lines[idx_voigt].split(':')[-1].split()
                v_kbar = [float(x) for x in parts]
                factor = 0.1 / _EV_A3_TO_GPA
                v = [x * factor for x in v_kbar]
                return np.array([
                    [v[0], v[5], v[4]],
                    [v[5], v[1], v[3]],
                    [v[4], v[3], v[2]],
                ])
            except Exception:
                pass
        return None
    except Exception:
        return None


def get_stress_voigt_kbar(path: str) -> list[float] | None:
    """Last 'Stress tensor Voigt ... (kbar)' line as [xx, yy, zz, yz, xz, xy], raw kBar."""
    last_stress = None
    try:
        with open(path, 'r') as f:
            for line in f:
                if "Stress tensor Voigt" in line and "(kbar)" in line:
                    parts = line.split(":")[-1].split()
                    if len(parts) >= 6:
                        last_stress = [float(x) for x in parts[:6]]
    except Exception:
        return None
    return last_stress


def get_free_energy(path: str) -> float | None:
    """Last 'siesta: FreeEng' value (eV) from a SIESTA .out file."""
    energy = None
    try:
        with open(path, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if "siesta: FreeEng" in line and "=" in line:
                    try:
                        energy = float(line.split('=')[1].split()[0])
                    except (IndexError, ValueError):
                        pass
    except Exception:
        return None
    return energy


def get_electric_dipole(path: str) -> np.ndarray | None:
    """Last 'siesta: Electric dipole (a.u.) = X  Y  Z' line from a SIESTA
    .out file, as a (3,) array in atomic units. SIESTA only prints this
    line for non-bulk systems (molecule/wire/slab -- any structure with
    at least one non-periodic direction); a genuinely 3D-periodic bulk
    calculation never prints it at all (a naive dipole moment is
    gauge-ambiguous there without a Berry-phase treatment -- see
    core.born_charges for the bulk-appropriate quantity instead). Returns
    None if the line is never found or fails to parse -- same contract as
    every other function in this module.

    Parses tolerantly (splits on the first '=', then takes the first 3
    numeric tokens via _parse_float_line) rather than assuming a fixed
    column width or that exactly 3 numbers always follow -- the precise
    spacing of this line has not been verified against a real SIESTA run
    in this environment; a format surprise should degrade to None (caller
    reports "not found"), not a silently wrong parse.
    """
    dipole = None
    try:
        with open(path, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if "Electric dipole" in line and "a.u." in line and "=" in line:
                    tail = line.split("=", 1)[1]
                    nums = _parse_float_line(tail)
                    if nums:
                        dipole = np.array(nums, dtype=float)
    except Exception:
        return None
    return dipole


def get_scf_convergence(path: str) -> tuple[bool, int | None]:
    """Returns (converged, iterations) from the LAST "SCF cycle converged
    after N iterations" line in a SIESTA .out file. `converged` is False
    (with `iterations` None) if that line never appears -- a conservative
    reading: rather than pattern-match SIESTA's various non-convergence
    messages (which vary by version and this module has never had a
    verified example of), treat "never confirmed converged" as the signal
    worth flagging, same spirit as this module's other "None on anything
    uncertain" functions.
    """
    converged = False
    iterations = None
    try:
        with open(path, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if "SCF cycle converged after" in line:
                    try:
                        iterations = int(line.split("after")[1].split()[0])
                        converged = True
                    except (IndexError, ValueError):
                        pass
    except Exception:
        return False, None
    return converged, iterations


def get_max_force(path: str) -> float | None:
    """Last 'Max' residual atomic force (eV/Ang) from a SIESTA .out file's
    "siesta: Atomic forces" block -- e.g. "   Max    0.689325" (optionally
    followed by "constrained" on the same line, when atoms are fixed; the
    unconstrained and constrained values are numerically identical unless
    constraints are actually declared). Same "last occurrence = final state"
    convention as get_free_energy/get_stress_tensor.
    """
    max_force = None
    try:
        with open(path, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 2 and parts[0] == "Max":
                    try:
                        max_force = float(parts[1])
                    except ValueError:
                        pass
    except Exception:
        return None
    return max_force


def check_scf_and_force(out_path: str) -> tuple[bool, float | None]:
    """Returns (scf_converged, max_force) for one SIESTA .out file -- a
    thin wrapper over get_scf_convergence/get_max_force above. Neither
    call is expensive (single sequential file read each), so it's fine to
    run per folder even for a large sweep.
    """
    scf_ok, _iterations = get_scf_convergence(out_path)
    max_force = get_max_force(out_path)
    return scf_ok, max_force


def report_quality_diagnostics(label: str, out_path: str, force_tolerance: float, f_out) -> None:
    """Prints (and persists, via core.cli.print_dual) an advisory
    numerical-quality warning for one folder if its SCF cycle never
    confirmed convergence, or its residual force exceeds
    force_tolerance -- silent when both are fine, matching this suite's
    "advisory only, don't clutter a clean run" convention (e.g.
    convergence_analysis.py's SCF gating, cohesive_analysis.py's
    --force-tolerance check).
    """
    scf_ok, max_force = check_scf_and_force(out_path)
    if not scf_ok:
        print_dual(color_text(
            f"  [WARNING] Could not confirm SCF convergence for {label} ('{out_path}') -- "
            "this energy may be unreliable.", 'yellow'), f_out)
    if max_force is not None and max_force > force_tolerance:
        print_dual(color_text(
            f"  [WARNING] Residual force on {label} ({max_force:.4f} eV/Ang) exceeds "
            f"--force-tolerance ({force_tolerance} eV/Ang) -- this geometry may not be "
            "relaxed.", 'yellow'), f_out)


def _parse_float_line(line: str) -> list[float] | None:
    """First 3 numeric tokens on the line, skipping any non-numeric ones.

    SIESTA prefixes some repeated blocks (e.g. the last "Stress tensor
    (static)" block in a relaxation run) with a leading "siesta:" label on
    each data line -- skip it and any other stray token instead of assuming
    the first 3 tokens are always the numbers.
    """
    nums: list[float] = []
    for tok in line.replace(',', ' ').split():
        try:
            nums.append(float(tok))
        except ValueError:
            continue
        if len(nums) == 3:
            break
    return nums if len(nums) == 3 else None
