"""Single source of truth for parsing SIESTA calculation logs (.out files).

Consolidates parsers that used to live independently in elastic_analysis.py,
strain_analysis.py, workfunction.py, wantibexos.py and cohesive_analysis.py:
Fermi energy, cell height (Lz), stress tensor (two flavors -- see note below)
and free energy.

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
"""

from __future__ import annotations

import re

import numpy as np

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
        for i in range(len(lines) - 1, 0, -1):
            if "outcell: Unit cell vectors" in lines[i]:
                nums = _parse_float_line(lines[i + 3])
                if nums:
                    z_len = float(np.linalg.norm(nums))
                break
    except Exception:
        pass
    return z_len


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


def _parse_float_line(line: str) -> list[float] | None:
    try:
        parts = line.replace(',', ' ').split()
        nums = [float(x) for x in parts[:3]]
        return nums if len(nums) >= 3 else None
    except Exception:
        return None
