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

import os
import re
import glob

import numpy as np

from stb.core.cli import color_text, print_dual

_STRAIN_FOLDER_RE = re.compile(r"strain_([a-zA-Z0-9]+)_(m?)(\d+\.\d+)")
_EV_A3_TO_GPA = 160.21766  # 1 eV/Angstrom^3 in GPa


def find_out_file(directory: str, label: str | None) -> str | None:
    """<directory>/<label>.out if `label` is known and that file exists,
    else the sole *.out in `directory` if there's exactly one (many SIESTA
    jobs redirect stdout to a generic name like calc.out, not
    <label>.out), else None if there's no match or the directory has
    several *.out files and no way to disambiguate. Shared by stb-status
    and stb-archive, which both need to locate "the" .out log for a
    directory without assuming its filename.
    """
    if label:
        candidate = os.path.join(directory, f"{label}.out")
        if os.path.isfile(candidate):
            return candidate
    matches = glob.glob(os.path.join(directory, "*.out"))
    return matches[0] if len(matches) == 1 else None


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


def find_strain_folders(base_dir: str) -> list[str]:
    """Finds every 'strain_*' run folder reachable from base_dir, supporting
    2 layouts: flat (strain_<dir>_<pct> directly under base_dir -- what you
    get by pointing base_dir at a single direction's own subfolder, e.g.
    'strain_runs/x') and nested (stb-strain's/stb-elasticInputs' own default
    output layout: one <direction>/ subfolder per strain direction under
    base_dir, each holding that direction's strain_<direction>_<pct>
    folders). This lets --dir (or the analog CWD-relative scan) point at
    either a single direction's subfolder (this direction only) or the
    top-level output directory itself (every direction found under it)
    without the caller needing to know which layout is present. Moved here
    from strain_analysis.py once elastic_analysis.py became a second
    consumer needing the identical flat-or-nested folder discovery.
    """
    direct = [os.path.join(base_dir, d) for d in sorted(os.listdir(base_dir))
              if os.path.isdir(os.path.join(base_dir, d)) and d.startswith('strain_')]
    if direct:
        return direct
    nested = []
    for d in sorted(os.listdir(base_dir)):
        sub = os.path.join(base_dir, d)
        if not os.path.isdir(sub):
            continue
        nested.extend(
            os.path.join(sub, e) for e in sorted(os.listdir(sub))
            if os.path.isdir(os.path.join(sub, e)) and e.startswith('strain_'))
    return nested


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


_DYNAMICS_OPTION_RE = re.compile(r"redata: Dynamics option\s*=\s*(.+)")
_GEOMETRY_STEP_RE = re.compile(r"Begin (?:\S+ opt\. move|MD step)\s*=\s*(\d+)")


def get_dynamics_type(path: str) -> str | None:
    """Raw 'redata: Dynamics option' string from a SIESTA .out log, e.g.
    'Single-point calculation', 'CG coord. optimization', 'Verlet MD run'.
    None if the file can't be read or doesn't have that line (older/atypical
    logs). Use categorize_dynamics() to turn this into a coarse
    single-point/relaxation/aimd bucket.
    """
    try:
        with open(path, 'r', errors='ignore') as f:
            for line in f:
                m = _DYNAMICS_OPTION_RE.search(line)
                if m:
                    return m.group(1).strip()
    except Exception:
        pass
    return None


def categorize_dynamics(dynamics_type: str | None) -> str:
    """Buckets get_dynamics_type()'s raw string into 'single-point',
    'relaxation' (CG/Broyden/FIRE coord. optimization), 'aimd' (Verlet/Nose/
    Parrinello-Rahman/Anneal MD), or 'unknown' (None, or a string that
    doesn't match any known SIESTA dynamics-option wording).
    """
    if dynamics_type is None:
        return "unknown"
    lower = dynamics_type.lower()
    if "single-point" in lower:
        return "single-point"
    if "opt." in lower or "optimization" in lower:
        return "relaxation"
    if any(kw in lower for kw in ("md run", "nose", "parrinello", "anneal")):
        return "aimd"
    return "unknown"


def count_geometry_steps(path: str) -> int:
    """Number of 'Begin MD step = N' / 'Begin <optimizer> opt. move = N'
    markers in a SIESTA .out log -- the number of MD steps for an AIMD run,
    or the number of geometry-optimization moves for a relaxation. 0 for a
    single-point run (neither marker appears) or if the file can't be read.
    """
    try:
        with open(path, 'r', errors='ignore') as f:
            text = f.read()
    except Exception:
        return 0
    return len(_GEOMETRY_STEP_RE.findall(text))


_MD_STEP_RE = re.compile(r"Begin MD step\s*=\s*(\d+)")


def get_md_trajectory(path: str) -> list[dict]:
    """Per-MD-step data from an AIMD .out log: splits the file into chunks
    at each 'Begin MD step = N' marker and extracts, from each chunk, the
    outcell 3x3 cell matrix (Angstrom; present every step regardless of
    whether the cell actually varies -- SIESTA reprints it either way), the
    KS energy ('siesta: E_KS(eV) =', eV) and the ionic temperature
    ('siesta: Temp_ion =', K). Chunk-scoped (not a single pass over the
    whole file) so a value from one MD step can never leak into another's
    entry.

    Returns a list of {'step': int, 'cell': np.ndarray|None,
    'E_KS': float|None, 'Temp_ion': float|None}, one per step found, in step
    order -- [] if no 'Begin MD step' marker exists at all (e.g. a
    single-point or plain geometry-relaxation run, not an AIMD one) or the
    file can't be read. Any individual field stays None if that step's
    chunk doesn't contain it, rather than aborting the whole parse -- built
    for stb-ani2traj's per-frame lattice/energy/temperature enrichment,
    which degrades gracefully field by field.
    """
    try:
        with open(path, 'r', errors='ignore') as f:
            text = f.read()
    except Exception:
        return []

    markers = list(_MD_STEP_RE.finditer(text))
    if not markers:
        return []

    steps = []
    for i, m in enumerate(markers):
        start = m.end()
        end = markers[i + 1].start() if i + 1 < len(markers) else len(text)
        chunk_lines = text[start:end].splitlines()

        cell = None
        for j, line in enumerate(chunk_lines):
            if "outcell: Unit cell vectors" in line:
                rows = [_parse_float_line(chunk_lines[j + off]) for off in (1, 2, 3)]
                if all(rows):
                    cell = np.array(rows, dtype=float)
                break

        e_ks = None
        temp_ion = None
        for line in chunk_lines:
            stripped = line.strip()
            if e_ks is None and stripped.startswith("siesta: E_KS(eV)"):
                try:
                    e_ks = float(stripped.split("=", 1)[1].split()[0])
                except (IndexError, ValueError):
                    pass
            elif temp_ion is None and stripped.startswith("siesta: Temp_ion"):
                try:
                    temp_ion = float(stripped.split("=", 1)[1].split()[0])
                except (IndexError, ValueError):
                    pass

        steps.append({'step': int(m.group(1)), 'cell': cell, 'E_KS': e_ks, 'Temp_ion': temp_ion})

    return steps


def get_mde_trajectory(label: str) -> list[dict]:
    """Per-MD-step thermodynamic data straight from SIESTA's own dedicated
    '<SystemLabel>.MDE' file -- a small, clean table (one row per MD step:
    Step, T (K), E_KS (eV), E_tot (eV), Vol (Ang^3), P (kBar)) SIESTA
    already writes for exactly this purpose, instead of regex-scraping it
    back out of scattered .out log lines the way get_md_trajectory() does
    for the cell/E_KS/Temp_ion it needs. Always named after SystemLabel,
    like .XV/.ANI/.HSX/.WFSX -- unlike the INPUT .fdf (see
    core/md_traj.py::read_md_timestep_fs's own docstring), so no
    --geometry-file-style override is needed here.

    Returns a list of {'step': int, 'T': float, 'E_KS': float,
    'E_tot': float|None, 'volume': float|None, 'pressure': float|None},
    one per data row, in file order -- [] if '<label>.MDE' doesn't exist
    or has no parseable rows. Vol/P are read positionally (columns 5/6)
    rather than by matching the header text (whose column labels contain
    embedded spaces, e.g. 'T (K)', making a naive whitespace-split
    header-name match unreliable) but degrade to None row-by-row if a row
    has fewer than 6 columns -- e.g. a fixed-shape-cell run without a
    barostat can omit Vol/P in older SIESTA versions.
    """
    path = f"{label}.MDE"
    if not os.path.exists(path):
        return []
    rows = []
    try:
        with open(path, 'r', errors='ignore') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split()
                try:
                    step = int(parts[0])
                    temp = float(parts[1])
                    e_ks = float(parts[2])
                except (IndexError, ValueError):
                    continue
                e_tot = None
                volume = None
                pressure = None
                try:
                    e_tot = float(parts[3])
                    volume = float(parts[4])
                    pressure = float(parts[5])
                except (IndexError, ValueError):
                    pass
                rows.append({'step': step, 'T': temp, 'E_KS': e_ks, 'E_tot': e_tot,
                             'volume': volume, 'pressure': pressure})
    except Exception:
        return []
    return rows


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


def get_spin_moment(path: str) -> float | None:
    """Last net magnetic (spin) moment, in Bohr magnetons, from a SIESTA
    .out file of a spin-polarized run. SIESTA only ever prints this for a
    'Spin polarized' calculation; a non-polarized run has no such line at
    all, so `None` here can mean either "not found" or "not a
    spin-polarized run" -- callers can't tell the two apart from this
    return value alone (same ambiguity get_electric_dipole already has
    for a bulk-vs-non-bulk run), which is fine for this suite's usage
    (adsorb_analysis.py already knows independently, from its own
    force_spin/config_extra.fdf bookkeeping, whether a site was run
    spin-polarized at all).

    Tries two known SIESTA wordings, in order (last matching line wins,
    same convention as every other parser here): the classic
    "siesta:    Total spin polarization (Qup-Qdown) =    X.XXXXXX" line,
    and a newer "Total spin moment:   X.XXXXXXX" wording seen in some
    SIESTA versions. **Neither has been verified against a real
    spin-polarized SIESTA .out in this environment** (no such fixture
    exists in this repository) -- flagged here deliberately; if a live
    run doesn't match, this returns None (never raises) rather than a
    wrong number, but the exact substring may need adjusting once
    checked against a real spin-polarized output.
    """
    moment = None
    try:
        with open(path, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if "Total spin polarization" in line and "=" in line:
                    try:
                        moment = float(line.split("=", 1)[1].split()[0])
                    except (IndexError, ValueError):
                        pass
                elif "Total spin moment" in line and ":" in line:
                    try:
                        moment = float(line.split(":", 1)[1].split()[0])
                    except (IndexError, ValueError):
                        pass
    except Exception:
        return None
    return moment


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


def read_fa_forces(fa_path: str) -> np.ndarray | None:
    """Reads ALL atomic forces (eV/Ang) from a SIESTA .FA file -- format:
    first line = atom count, then one row per atom '<1-based index> Fx Fy
    Fz'. Returns an (natoms, 3) array, or None on any read/parse failure
    (same fail-soft contract as every other parser in this module).

    Promoted here from mlff_analysis.py once neb_cycle.py became a second
    consumer -- extract-on-second-use, same policy as the rest of core/.
    mlff_analysis.py's own read_fa_forces() (which itself generalized
    her_analysis.py::read_fa_force's single-atom version) now just imports
    this instead of keeping a local copy.
    """
    try:
        with open(fa_path) as f:
            lines = f.readlines()
    except OSError:
        return None
    try:
        n = int(lines[0].split()[0])
    except (IndexError, ValueError):
        return None
    if len(lines) < n + 1:
        return None
    try:
        forces = np.array([[float(x) for x in lines[i + 1].split()[1:4]] for i in range(n)])
    except (IndexError, ValueError):
        return None
    return forces


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
