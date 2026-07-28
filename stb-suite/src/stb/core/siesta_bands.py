"""Shared SIESTA .bands-file parsing and k-path/VBM-CBM helpers.

Extracted from bands.py once stb-fatbands became a second consumer of this
logic (both tools need the same k-path positions, high-symmetry ticks, and
VBM/CBM/shift machinery to plot on/reference the same x-axis). Pure
functions operating on an already-parsed dic_bands -- no plotting or
report-formatting here, those stay tool-specific in bands.py/fatbands.py.
"""

import re
import numpy as np
from stb.core.deps import require_sisl
from stb.core.siesta_log import find_out_file, get_fermi_energy


def resolve_fermi_energy_hierarchy(fermi, bands_file, fermi_file, label, dirname="."):
    """Resolves a Fermi energy in priority order: an explicit value >
    a companion .bands file > an explicit .out log > an auto-detected .out
    log in `dirname` (via find_out_file, NOT assumed to be <label>.out --
    many real SIESTA jobs redirect stdout to a generic name instead).
    Returns (fermi_energy, source_description) or (None, None) if nothing
    was found.

    First written for stb-wfdensity's --band vbm/cbm; extracted here once
    stb-sts's --shift fermi became a second consumer of the identical
    hierarchy, so both tools resolve a Fermi energy the same way instead
    of stb-sts only ever accepting an explicit --fermi value.
    """
    if fermi is not None:
        return fermi, "--fermi (explicit value)"
    if bands_file:
        fermi_energy, _, _, _ = read_data(bands_file)
        return fermi_energy, f"'{bands_file}' (--bands-file)"
    if fermi_file:
        fermi_energy = get_fermi_energy(fermi_file)
        return fermi_energy, f"'{fermi_file}' (--fermi-file)"
    out_file = find_out_file(dirname, label)
    if out_file:
        fermi_energy = get_fermi_energy(out_file)
        if fermi_energy is not None:
            return fermi_energy, f"'{out_file}' (auto-detected .out)"
    return None, None


def _is_gamma(label):
    clean_text = re.sub(r'[^a-zA-Z\s]', '', label).lower()
    return "gamma" in clean_text.split()


def read_data(file_path="siesta.bands"):
    # SIESTA .bands format:
    #   line 0: Fermi energy
    #   line 1: kmin kmax
    #   line 2: emin emax
    #   line 3: nbands nspin nk
    #   body:   nk k-blocks (first line has the k value, then as many
    #           continuation lines as needed to reach nbands*nspin values)
    #   footer: a line with the count of high-symmetry points, followed by
    #           that many "position 'LABEL'" lines
    with open(file_path, "r") as f:
        lines = f.readlines()

    fermi_energy = float(lines[0].split()[0])
    nbands, nspin, nk = (int(x) for x in lines[3].split())
    values_per_k = nbands * nspin

    dic_bands = {}
    idx = 4
    for _ in range(nk):
        tokens = lines[idx].split()
        key = float(tokens[0])
        values = tokens[1:]
        idx += 1
        while len(values) < values_per_k:
            values.extend(lines[idx].split())
            idx += 1
        # First nbands values are the spin-up channel, the next nbands (if
        # any) are spin-down -- confirmed against sisl's own bandsSileSiesta
        # ("eb[ik, :, :] = l.reshape(ns, no)"), not interleaved per band.
        dic_bands[key] = np.array(values, dtype=float).reshape(nspin, nbands)

    # Whatever remains is deterministically the high-symmetry footer: we
    # already consumed exactly nk k-blocks by token count, so there is no
    # ambiguity left between a continuation line and the count line.
    while idx < len(lines) and not lines[idx].split():
        idx += 1
    n_high_sym = int(lines[idx].split()[0])
    idx += 1

    high_sym = []
    for _ in range(n_high_sym):
        parts = lines[idx].split()
        label = " ".join(parts[1:]).strip("'\"")
        high_sym.append([parts[0], label])
        idx += 1

    return fermi_energy, high_sym, dic_bands, nspin


def shift_bands(dic, val):
    return {k: arr - val for k, arr in dic.items()}


def read_eig_mesh(eig_file, kp_file=None):
    """Reads a SIESTA .EIG file (eigenvalues at every k-point of the SCF
    k-mesh, not just a high-symmetry path) via sisl. Returns (fermi_energy,
    dic_mesh, nspin, kpoints) in the same (nspin, nbands)-per-k shape as
    read_data() above, so cbm_vbm() works on either source unchanged.

    Extracted from bands.py (--eig-file mesh-vs-line gap comparison) once
    stb-dos became a second consumer, needing the same reader to resolve a
    --shift vbm/cbm reference from a companion .EIG when no .bands file is
    available.
    """
    sisl = require_sisl()
    sile = sisl.get_sile(eig_file)
    fermi_energy = sile.read_fermi_level()
    eigs = sile.read_data()
    nspin = eigs.shape[0]
    nk = eigs.shape[1]
    # Same (nspin, nbands)-per-k convention as read_data(); add Ef back so
    # values are absolute eV, matching the rest of this module.
    dic_mesh = {ik: eigs[:, ik, :] + fermi_energy for ik in range(nk)}

    kpoints = None
    if kp_file:
        # The .EIG file itself has no k-vectors, only a bare mesh index --
        # the .KP file (same calculation) holds the actual Cartesian
        # (kx, ky, kz), in 1/Ang (sisl converts from the file's raw 1/Bohr).
        kp_sile = sisl.get_sile(kp_file)
        kpoints, _weights = kp_sile.read_data()
        if kpoints.shape[0] != nk:
            raise ValueError(
                f"--kp-file '{kp_file}' has {kpoints.shape[0]} k-points but "
                f"'{eig_file}' has {nk} -- they must be from the same calculation."
            )
    return fermi_energy, dic_mesh, nspin, kpoints


def _band_extrema(values_by_k, fermi_energy, gap_tol):
    # Start value as infinity
    vbm = -np.inf
    cbm = np.inf
    vbm_k = None
    cbm_k = None
    # Direct gap: the smallest CBM(k) - VBM(k) at one and the same k --
    # tracked alongside the (possibly different-k) global VBM/CBM below.
    direct_gap = np.inf
    direct_k = None
    for k, band in values_by_k.items():
        below = band[band <= fermi_energy]
        above = band[band > fermi_energy]
        if below.size > 0:
            local_vbm = np.nanmax(below)
            if local_vbm > vbm:
                vbm = local_vbm
                vbm_k = k
        if above.size > 0:
            local_cbm = np.nanmin(above)
            if local_cbm < cbm:
                cbm = local_cbm
                cbm_k = k
            if below.size > 0:
                local_gap = local_cbm - local_vbm
                if local_gap < direct_gap:
                    direct_gap = local_gap
                    direct_k = k
    if vbm_k is None or cbm_k is None:
        raise ValueError(
            f"No occupied/empty states found on {'both sides' if vbm_k is None and cbm_k is None else ('the occupied side' if vbm_k is None else 'the empty side')} "
            f"of Fermi energy {fermi_energy:.6f} eV -- check that this Fermi energy actually matches the eigenvalue file."
        )
    # Indirect (= fundamental) gap: CBM - VBM regardless of whether they sit
    # at the same k. Always <= direct_gap, since direct_gap is the same
    # quantity restricted to matching k. When they coincide (within
    # gap_tol) the fundamental gap is itself direct.
    indirect_gap = cbm - vbm if cbm > vbm else 0.0  # Avoid negative values
    if indirect_gap < gap_tol:
        gap_type = "Metallic"
    elif direct_k is not None and (direct_gap - indirect_gap) < gap_tol:
        gap_type = "Direct"
    else:
        gap_type = "Indirect"
    return vbm, cbm, vbm_k, cbm_k, indirect_gap, gap_type, direct_gap, direct_k


def cbm_vbm(fermi_energy, dic_bands, nspin, gap_tol=0.01, k_format=None):
    # Pure computation, no printing -- callers render the result (console
    # and/or bands_analysis.txt) via write_analysis_report(), so the same
    # numbers are never formatted two different ways in two places.
    k_format = k_format or (lambda k: f"{k:.6f}")
    combined = _band_extrema(
        {k: arr.reshape(-1) for k, arr in dic_bands.items()}, fermi_energy, gap_tol
    )
    result = {"combined": combined, "spins": [], "half_metallic": False, "k_format": k_format}
    if nspin == 2:
        for s in range(nspin):
            spin_result = _band_extrema(
                {k: arr[s] for k, arr in dic_bands.items()}, fermi_energy, gap_tol
            )
            result["spins"].append(spin_result)
        gaps = [spin_result[4] for spin_result in result["spins"]]  # indirect gap
        result["half_metallic"] = min(gaps) < gap_tol and max(gaps) >= gap_tol

    return result


def select_band_vbm_cbm(states_by_k, nspin, fermi_energy, gap_tol, which):
    """Builds a dic_bands-shaped dict (key -> (nspin, nbands) eigenvalue
    array, key meaning irrelevant to cbm_vbm) directly from a .WFSX's own
    states_by_k (as returned by core.siesta_wfsx.read_wfsx_states), so
    cbm_vbm() can be reused to find the global VBM/CBM without needing a
    companion .bands file at all -- then maps the result back to
    (k_index, spin, band_index). Extracted from stb-wfdensity once
    stb-effmass became a second consumer of the exact same "find VBM/CBM
    directly from a WFSX" need."""
    dic_bands = {}
    for k_index, block in enumerate(states_by_k):
        nbands = min(len(block[s].eig) for s in block)
        arr = np.empty((nspin, nbands))
        for s in range(nspin):
            arr[s] = np.asarray(block[s].eig)[:nbands]
        dic_bands[k_index] = arr

    result = cbm_vbm(fermi_energy, dic_bands, nspin, gap_tol)
    vbm, cbm, vbm_k, cbm_k, _, _, _, _ = result["combined"]
    target_k = vbm_k if which == "vbm" else cbm_k
    target_val = vbm if which == "vbm" else cbm

    arr = dic_bands[target_k]
    spin, band = np.unravel_index(np.argmin(np.abs(arr - target_val)), arr.shape)
    return int(target_k), int(spin), int(band)
