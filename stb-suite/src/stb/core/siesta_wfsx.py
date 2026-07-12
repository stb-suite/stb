"""Shared SIESTA .WFSX (wavefunction coefficients) loading helpers.

Extracted from fatbands.py once stb-wfdensity/stb-sts/stb-coop became
further consumers of the same "resolve label/file/wfsx/hsx/geometry into a
usable sisl parent object" logic.

Three different `mode`s exist for `load_parent()` because the *kind* of
parent object needed genuinely differs by tool, not just "accurate vs
approximate":

- "hamiltonian_preferred" (stb-fatbands): tries .HSX (a real
  sisl.physics.Hamiltonian, needed for the true non-orthogonal-basis
  overlap-weighted norm2/Sk projection), falling back to a .fdf(+ion
  files)-derived Geometry (orthogonal-basis approximation, printed warning)
  if no .HSX is found.
- "hamiltonian_required" (stb-coop): only .HSX -- COOP needs Sk() and COHP
  needs Hk(), both of which require a real Hamiltonian; there is no
  meaningful Geometry-only fallback for either.
- "orbitals_required" (stb-wfdensity, stb-sts): only .fdf(+ion files) --
  NEVER tries .HSX at all. Verified live: a Geometry read from .HSX has
  every atom's orbitals with maxR() == -1.0 (no radial shape, just the
  linear-algebra bookkeeping norm2/Sk need) -- useless for
  sisl.physics.electron.wavefunction()'s real-space evaluation, which
  requires atom.maxR() > 0. Only a Geometry read from .fdf with its
  .ion/.ion.xml files sitting alongside carries real radial functions.
  Also verified live that .HSX isn't even needed as a stepping stone here:
  reading the .WFSX with a bare .fdf-derived Geometry as `parent` (no
  Hamiltonian at all) gives numerically identical wavefunction densities.
"""

import os

from stb.core.deps import require_sisl


def resolve_wfsx_path(label, suffixes=(".bands.WFSX", ".WFSX", ".selected.WFSX")):
    """Tries each of `suffixes`, in order, returning the first
    <label><suffix> that exists (with a printed note for anything after the
    first). Returns None if none exist.

    Three real SIESTA/sisl WFSX naming conventions exist, confirmed live in
    this environment, one per fdf mechanism that triggers a wavefunction
    write:
    - "<label>.bands.WFSX": WFS.Write.For.Bands T + %block BandLines.
    - "<label>.selected.WFSX": WriteWaveFunctions T + %block WaveFuncKPoints
      (an explicit k-point list -- this is what a full-BZ-mesh WFSX for
      stb-sts/stb-coop actually looks like; WriteWaveFunctions T alone,
      with no WaveFuncKPoints block, silently writes nothing).
    - "<label>.WFSX": not produced by either mechanism above in testing,
      but kept as a fallback name in case of a differently-configured run
      or a user-renamed file.

    stb-fatbands wants (".bands.WFSX", ".WFSX", ".selected.WFSX") (its
    typical input is a band-path WFSX). stb-wfdensity/stb-sts/stb-coop want
    (".WFSX", ".selected.WFSX", ".bands.WFSX") -- their typical input is a
    full-BZ-mesh WFSX.
    """
    for i, suffix in enumerate(suffixes):
        candidate = f"{label}{suffix}"
        if os.path.isfile(candidate):
            if i > 0:
                print(f"[INFO] No '{label}{suffixes[0]}' found; using '{candidate}' instead.")
            return candidate
    return None


def load_parent(label, hsx_file, geometry_file, mode="hamiltonian_preferred"):
    """Returns (parent, source_path, is_orthogonal_approx).

    `mode` selects which kind of parent this tool actually needs -- see
    the module docstring above. `is_orthogonal_approx` is only meaningful
    for "hamiltonian_preferred" (True when the .fdf-Geometry fallback was
    used instead of a real Hamiltonian); it is always False for the other
    two modes, since neither has an "approximate but still usable"
    fallback path -- either the required file was found, or `parent` is
    None and the caller must error out.
    """
    sisl = require_sisl()

    if mode == "orbitals_required":
        if geometry_file:
            if not os.path.isfile(geometry_file):
                print(f"[WARNING] --geometry-file '{geometry_file}' not found.")
                return None, None, False
            geometry = sisl.get_sile(geometry_file).read_geometry()
            return geometry, geometry_file, False
        if label:
            candidate = f"{label}.fdf"
            if os.path.isfile(candidate):
                geometry = sisl.get_sile(candidate).read_geometry()
                return geometry, candidate, False
        return None, None, False

    if hsx_file is None and label:
        candidate = f"{label}.HSX"
        if os.path.isfile(candidate):
            hsx_file = candidate

    if hsx_file:
        if not os.path.isfile(hsx_file):
            print(f"[WARNING] --hsx-file '{hsx_file}' not found.")
            return None, None, False
        H = sisl.get_sile(hsx_file).read_hamiltonian()
        return H, hsx_file, False

    if mode == "hamiltonian_required":
        return None, None, False

    # mode == "hamiltonian_preferred": fall back to an approximate Geometry.
    if geometry_file:
        if not os.path.isfile(geometry_file):
            print(f"[WARNING] --geometry-file '{geometry_file}' not found.")
            return None, None, True
        geometry = sisl.get_sile(geometry_file).read_geometry()
        return geometry, geometry_file, True

    if label:
        candidate = f"{label}.fdf"
        if os.path.isfile(candidate):
            geometry = sisl.get_sile(candidate).read_geometry()
            return geometry, candidate, True

    return None, None, True


def _geometry_of(parent):
    return getattr(parent, "geometry", parent)


def read_wfsx_states(wfsx_path, parent):
    """Groups wfsxSileSiesta.yield_eigenstate()'s flat (k, spin) iteration
    (k-major, spin-minor -- confirmed from sisl's own
    `product(range(nk), itt_spin)` iteration order) into one dict per
    k-point, keyed by spin index (0 for nspin=1). Returns (sizes,
    states_by_k) where states_by_k[i] == {spin: EigenstateElectron}.

    Materializes the whole file in memory -- fine for a moderate band-path
    WFSX (stb-fatbands, stb-wfdensity's --band vbm/cbm search), but use
    `iter_wfsx_states_by_k` instead for a potentially large full-BZ mesh
    (stb-sts, stb-coop).
    """
    sisl = require_sisl()
    wfsx = sisl.get_sile(wfsx_path, parent=parent)
    sizes = wfsx.read_sizes()

    states_by_k = []
    current = None
    for state in wfsx.yield_eigenstate():
        spin = state.info.get("spin", 0)
        if spin == 0:
            if current is not None:
                states_by_k.append(current)
            current = {}
        current[spin] = state
    if current is not None:
        states_by_k.append(current)

    return sizes, states_by_k


def iter_wfsx_states_by_k(wfsx_path, parent):
    """Generator version of read_wfsx_states -- yields one (k_index,
    {spin: EigenstateElectron}) at a time without materializing the whole
    file, for tools that iterate a potentially large full-BZ mesh."""
    sisl = require_sisl()
    wfsx = sisl.get_sile(wfsx_path, parent=parent)
    k_index = -1
    current = None
    for state in wfsx.yield_eigenstate():
        spin = state.info.get("spin", 0)
        if spin == 0:
            if current is not None:
                yield k_index, current
            k_index += 1
            current = {}
        current[spin] = state
    if current is not None:
        yield k_index, current


def read_k_weights(wfsx_path, parent):
    """Per-k-point BZ weights from the .WFSX file's own header info,
    needed to properly average any DOS-like quantity (STS, COOP/COHP)
    across a full k-mesh."""
    sisl = require_sisl()
    wfsx = sisl.get_sile(wfsx_path, parent=parent)
    _k, weight, _nwf = wfsx.read_info()
    return weight
