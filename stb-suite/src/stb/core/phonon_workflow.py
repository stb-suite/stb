"""Shared Phonopy displacement-generation and force-constant-loading logic
for SIESTA phonon calculations.

Extracted from stb-phononsCreate (build_phonon_displacements/
write_displacement_folders) and stb-phononsPos (load_phonon_with_force_
constants, detect_system_label) once stb-raman became a second, real
consumer -- same extract-on-second-use policy as core/heterostructure.py.
Everything specific to one tool's own extra features (ML pre-relax,
symmetry table, thermal properties, band structure, DOS/PDOS, mode
freezing) stays in phonons_create.py/phonons_pos.py -- only the genuinely
shared Phonopy-object plumbing moved here.

get_gamma_modes() is new code: nothing in the suite previously exposed
Gamma-point eigenvectors as reusable data (phonons_pos.py only fetches
eigenvectors transiently, via run_mesh(with_eigenvectors=True), for
PDOS/thermal-displacement purposes).
"""

import os
import re
import glob
import shutil
import subprocess

import numpy as np
from phonopy import Phonopy
from phonopy.interface.siesta import write_siesta, get_physical_units

from stb.core.cli import color_text, print_dual

_LABEL_RE = re.compile(r'^\s*SystemLabel\s+(\S+)', re.IGNORECASE | re.MULTILINE)


def detect_system_label(phonon_dir):
    """Auto-detect SystemLabel from the calc.fdf/structure.fdf copied into
    the first disp-* folder, same regex-based approach as
    hubbardu.get_system_label/wantibexos.get_system_label. Returns None
    (not "siesta") when nothing is found, so callers can tell
    "not detected" apart from a genuine label of "siesta". Moved here
    verbatim from phonons_pos.py once stb-ramanModes needed the identical
    auto-detection.
    """
    disp_dirs = sorted(glob.glob(os.path.join(phonon_dir, "disp-*")))
    if not disp_dirs:
        return None
    for fdf_path in glob.glob(os.path.join(disp_dirs[0], "*.fdf")):
        try:
            with open(fdf_path) as f:
                match = _LABEL_RE.search(f.read())
        except OSError:
            continue
        if match:
            return match.group(1)
    return None


def build_phonon_displacements(unitcell, supercell_matrix, distance_ang):
    """Constructs the Phonopy object and generates the displaced
    supercells via generate_displacements(). `distance_ang` is in
    Angstrom, converted internally to Bohr: phonopy's SIESTA interface
    keeps cell/positions internally in bohr (not Angstrom, unlike every
    other calculator interface), so distance has to be converted to that
    same unit or the real Cartesian displacement ends up ~1.89x smaller
    than requested -- verified numerically against this tool's own output
    (same pitfall documented in phonons_create.py's original inline
    version of this code).

    Returns (phonon, supercells) -- phonon.dataset has the
    symmetry-reduced displacement pattern (spglib, via Phonopy itself),
    supercells is the list of displaced structures (some entries may be
    None, same convention as phonopy.supercells_with_displacements).
    """
    bohr_to_angstrom = get_physical_units().Bohr
    phonon = Phonopy(unitcell, supercell_matrix=supercell_matrix, calculator="siesta")
    phonon.generate_displacements(distance=distance_ang / bohr_to_angstrom)
    supercells = phonon.supercells_with_displacements
    return phonon, supercells


def write_displacement_folders(output_root, phonon, supercells, structure_filename,
                                calc_path, pseudos_to_copy):
    """Writes disp-NNN/ (displaced structure + calc.fdf copied + pseudos)
    and saves phonopy_disp.yaml. Returns (folders, yaml_path) -- folders
    is the list of disp-NNN/ paths actually written (skipping any None
    entry in `supercells`).
    """
    folders = []
    for i, scell in enumerate(supercells):
        if scell is None:
            continue

        folder_name = os.path.join(output_root, f"disp-{i + 1:03d}")
        os.makedirs(folder_name, exist_ok=True)

        disp_struct_path = os.path.join(folder_name, os.path.basename(structure_filename))
        write_siesta(disp_struct_path, scell)

        shutil.copy(calc_path, os.path.join(folder_name, os.path.basename(calc_path)))

        for pseudo_path in pseudos_to_copy:
            pseudo_filename = os.path.basename(pseudo_path)
            shutil.copy(pseudo_path, os.path.join(folder_name, pseudo_filename))

        folders.append(folder_name)

    yaml_path = os.path.join(output_root, "phonopy_disp.yaml")
    phonon.save(yaml_path)
    return folders, yaml_path


def load_phonon_with_force_constants(phonon_dir, system_label, has_embedded_fc, f_out=None):
    """Builds FORCE_SETS from disp-*/<system_label>.FA (subprocess
    phonopy-init --siesta -f, falling back to plain `phonopy` on older
    installs -- same fallback phonons_pos.py already used) and loads the
    Phonopy object with force constants ready. Skips straight to loading
    when `has_embedded_fc` (an ML-sourced phonopy_disp.yaml, e.g. from
    stb-phononsML, already has force constants embedded -- no SIESTA
    .FA/FORCE_SETS involved at all).

    Leaves the process chdir'd into `phonon_dir` on return (phonopy.load
    needs phonopy_disp.yaml/FORCE_SETS to be resolved relative to the
    current directory) -- callers that write further output there
    (thermal properties, bands, ...) rely on this, same as the original
    inline code in phonons_pos.py. Returns
    (phonon, internal_to_angstrom, original_dir); callers must
    os.chdir(original_dir) themselves once done.

    Raises SystemExit(1) (after printing a report-visible error via
    print_dual) on any failure -- same fail-loud behavior as the code
    this was extracted from, since a partial/corrupted phonon load isn't
    recoverable mid-run.
    """
    import phonopy  # local import: phonons_pos.py already guards this at
    # module level with a friendly install hint if missing; core modules
    # shouldn't duplicate that guard for every caller.

    if has_embedded_fc:
        print_dual("Force constants already embedded in phonopy_disp.yaml "
                    "(ML-computed, e.g. via stb-phononsML) -- skipping SIESTA "
                    ".FA/FORCE_SETS extraction.", f_out)
    else:
        fa_files_abs = sorted(glob.glob(os.path.join(phonon_dir, "disp-*", f"{system_label}.FA")))

        if not fa_files_abs:
            print_dual(color_text(
                f"[ERROR] No {system_label}.FA files found in {phonon_dir}/disp-*", 'red'), f_out)
            print_dual(color_text("Make sure SIESTA calculations finished successfully.", 'yellow'), f_out)
            raise SystemExit(1)

        fa_files_rel = [os.path.relpath(f, phonon_dir) for f in fa_files_abs]
        print_dual(f"Force files found : {len(fa_files_rel)} (label '{system_label}')", f_out)

        # phonopy >=4 moved the FORCE_SETS-creation flags ("-f"/"--siesta")
        # out of the main `phonopy` CLI into a separate `phonopy-init`
        # command; older installs (<4) only have `phonopy`. Prefer the new
        # command when present, fall back otherwise.
        force_sets_bin = "phonopy-init" if shutil.which("phonopy-init") else "phonopy"
        cmd = [force_sets_bin, "--siesta", "-f"] + fa_files_rel

        try:
            subprocess.check_call(cmd, cwd=phonon_dir, stdout=subprocess.DEVNULL)
            print_dual(color_text("FORCE_SETS       : generated successfully", 'green'), f_out)
        except subprocess.CalledProcessError as e:
            print_dual(color_text(f"[ERROR] Failed to generate FORCE_SETS. Phonopy error: {e}", 'red'), f_out)
            raise SystemExit(1)

    print("[INFO] Initializing Phonopy API and loading FORCE_SETS ...")
    original_dir = os.getcwd()
    try:
        os.chdir(phonon_dir)
        phonon = phonopy.load("phonopy_disp.yaml")
    except Exception as e:
        print_dual(color_text(f"[ERROR] Could not load phonopy data: {e}", 'red'), f_out)
        os.chdir(original_dir)
        raise SystemExit(1)

    # phonon's internal cell/positions are in bohr (SIESTA-interface
    # convention, calculator="siesta") for a stb-phononsCreate-sourced run,
    # but already in native Angstrom (calculator=None) for an ML-sourced one
    # -- internal_to_angstrom is what actually converts *this* phonon
    # object's internal units to real Angstrom, 1.0 (no-op) for the ML case.
    # Using the real Bohr-to-Angstrom constant unconditionally here was a
    # real, previously-fixed bug (see phonons_pos.py's original note) --
    # preserved exactly here.
    true_bohr_to_angstrom = get_physical_units().Bohr
    internal_to_angstrom = 1.0 if has_embedded_fc else true_bohr_to_angstrom
    return phonon, internal_to_angstrom, original_dir


def get_gamma_modes(phonon, exclude_acoustic=True, extra_trivial_tol_thz=None):
    """Runs phonon.run_qpoints([[0, 0, 0]]) and returns
    (frequencies, mode_indices, n_extra_trivial) at Gamma: frequencies in
    THz, mode_indices the corresponding 0-based Phonopy band indices --
    usable directly as the `band_index` in phonon.run_modulations([[q_point,
    band_index, amplitude, phase]]) (see displace_along_mode below), rather
    than exposing raw eigenvectors and reimplementing Phonopy's own
    mass-weighting/normalization math here.

    `exclude_acoustic=True` (default) drops the first 3 modes (pure
    translation, omega ~ 0 at Gamma) -- these aren't real Raman-active
    vibrations. This is the right, and ONLY, trivial-mode count for a
    genuinely 3D/2D/1D (at least one real periodic direction) structure --
    free rotation isn't a zero-energy motion there, so no other mode is
    ever trivial no matter how soft it is.

    `extra_trivial_tol_thz`: for a genuinely isolated (0D) structure only,
    a real molecule has 3 MORE trivial modes at Gamma beyond translation --
    free rotation (3 for a non-linear molecule, 2 for a linear one, since
    spinning about the molecular axis isn't a real degree of freedom for
    point masses on a line). Verified live with a relaxed H2O molecule via
    MACE-MP-0: bands 0-2 (translation, ~-1.9 to -1.5 THz -- numerical box
    noise) and bands 3-5 (rotation, ~-0.01 to 0.22 THz -- distinctly
    closer to zero) are both far below the genuine vibrational modes at
    44.5/111.6/115.2 THz. When given (in THz), after the first 3 bands, up
    to 3 MORE contiguous bands are also dropped as long as their
    |frequency| stays below this tolerance -- naturally stops at 2 extra
    for a linear molecule (the 3rd real mode is a genuine, much higher-
    frequency vibration) and at 3 for a non-linear one. The hard cap of 3
    extra is what makes this safe to call unconditionally once a caller
    has confirmed 0D (e.g. via core.kspace.detect_vacuum_axes): it can
    never eat into a real vibrational mode of a 3D/2D/1D structure, because
    those callers must never pass this parameter at all -- there's no
    "maybe rotational" ambiguity to resolve there in the first place.
    """
    phonon.run_qpoints([[0, 0, 0]])
    frequencies = phonon.qpoints.frequencies[0]
    mode_indices = np.arange(len(frequencies))
    if not exclude_acoustic:
        return frequencies, mode_indices, 0

    n_trivial = 3
    if extra_trivial_tol_thz is not None:
        for i in range(3, min(6, len(frequencies))):
            if abs(frequencies[i]) < extra_trivial_tol_thz:
                n_trivial += 1
            else:
                break

    return frequencies[n_trivial:], mode_indices[n_trivial:], n_trivial - 3


def displace_along_mode(phonon, band_index, amplitude_ang, internal_to_angstrom, sign=1.0):
    """Displaces the Gamma-point primitive cell along phonon mode
    `band_index` (0-based Phonopy band index at q=[0,0,0], e.g. from
    get_gamma_modes) by `sign * amplitude_ang` (Angstrom).

    Uses phonon.run_modulations -- the SAME machinery phonons_pos.py's
    --freeze-unstable-mode already uses to turn a phonon mode into a real
    Cartesian displacement pattern (mass-weighting and normalization
    handled by Phonopy itself, not reimplemented here) -- with dimension
    [1, 1, 1] (no supercell expansion: Gamma has no cell-to-cell phase
    variation to visualize, so the "modulated supercell" IS the primitive
    cell). run_modulations' amplitude parameter isn't itself in Angstrom
    (same note as the freeze-unstable-mode code this mirrors), so a
    unit-amplitude probe is run first and the real displacement rescaled
    linearly from it.

    Returns a PhonopyAtoms (phonon.modulation.modulated_supercells[0]) --
    same object type/unit-convention as build_phonon_displacements's
    output, ready for phonopy.interface.siesta.write_siesta.
    """
    q_point = [0, 0, 0]
    phonon.run_modulations([1, 1, 1], [[q_point, band_index, 1.0, 0.0]])
    mod = phonon.modulation
    probe_shift = np.linalg.norm(
        mod.modulated_supercells[0].positions - mod.supercell.positions,
        axis=1).max() * internal_to_angstrom
    if probe_shift == 0:
        raise ValueError(f"Mode {band_index} has zero displacement amplitude at a unit "
                          "probe -- cannot calibrate a real Angstrom displacement.")
    scale = sign * amplitude_ang / probe_shift

    phonon.run_modulations([1, 1, 1], [[q_point, band_index, scale, 0.0]])
    mod = phonon.modulation
    return mod.modulated_supercells[0]
