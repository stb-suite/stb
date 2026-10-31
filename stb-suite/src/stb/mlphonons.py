#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""stb-mlphonons (7.3): phonons from a MACE potential, by finite displacements (phonopy).

Three modes: one volume (bands, DOS, thermal properties, optionally against the raw foundation
model), --qha (quasi-harmonic volume scan: thermal expansion, B(T), Cp, Grueneisen) and
--check-convergence (the same calculation in several supercells, compared q-point by q-point).

The report is built in numbered sections with every warning collected in a final WARNINGS section
and the messages of MACE/torch/phonopy only counted. Figures: PNGs in the output folder, --show
opens them, --save-gnuplot writes .dat + .gplot into <output>/plot/ and runs gnuplot there.

Units: phonopy works in the PRIMITIVE cell (primitive_matrix 'auto', core/mace_phonons.py), so
its thermal properties are per primitive cell; the report says so and also gives them per atom.
The QHA multiplies them by mace_phonons.cell_factor() to put them on the input cell, the cell
whose volumes and energies it fits -- it once mixed the two (8-atom diamond cell: free energy 4x
too small, thermal expansion 4x too small) and printed the static B0 in eV/Ang^3 as GPa.
"""

from stb import __version__ as VERSION

import os
import sys
import copy
import time
import argparse
from datetime import datetime

import numpy as np
from ase import Atoms
from ase.io import write as ase_write
from pymatgen.io.ase import AseAtomsAdaptor

from stb.core import structure_io, kspace, model_card
from stb.core.cli import (color_text, show_intro, print_dual, print_table,
                          capture_library_noise, display_available, Sections, run_gnuplot)
from stb.core.phonon_workflow import IMAGINARY_MODE_TOL_THZ as _IMAGINARY_MODE_TOL_THZ

REPORT_FILE = "stb_mlphonons_report.txt"
PLOT_DIR = "plot"

THZ_TO_CM = 33.35641
THZ_TO_MEV = 4.135667696
H_OVER_KB = 47.99243           # K per THz
EV_A3_TO_GPA = 160.21766208
KJMOL_TO_MEV = 10.3642688      # meV per kJ/mol
R_GAS = 8.314462618            # J/K/mol = k_B per particle
AMU_A3_TO_G_CM3 = 1.66053906660
SOUND_Q = 0.01                 # 1/Ang: |q| at which the acoustic slopes are taken
REPORT_TEMPS = (100, 300, 500, 800, 1000, 1500, 2000)
SMALL_WIDTH = 10.0             # Ang: a supercell thinner than this is flagged

_SERIES_COLORS = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown']
_GP_COLORS = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#E69F00', '#56B4E9']

# Same 3-format convention already used by stb-ani2traj/stb-mlmd for
# multi-frame trajectories (xsf/pdb natively carry a lattice, extxyz via a
# Lattice= header tag) -- duplicated here rather than imported, same as
# mlmd.py already does, since it's tiny and tied to each tool's own CLI
# --animate-format choices/help text.
OUTPUT_FORMATS = {
    "xsf": ("xsf", ".xsf"),
    "pdb": ("proteindatabank", ".pdb"),
    "xyz": ("extxyz", ".xyz"),
}


def fmt(v, spec="{:.2f}", none="-"):
    return none if v is None or (isinstance(v, float) and not np.isfinite(v)) else spec.format(v)




def perpendicular_widths(cell):
    """Distance between opposite faces of a cell, per axis (Ang)."""
    cell = np.asarray(cell, float)
    v = abs(np.linalg.det(cell))
    return [v / np.linalg.norm(np.cross(cell[(i + 1) % 3], cell[(i + 2) % 3])) for i in range(3)]


def model_range(calc, card):
    """(r_max, interactions) of the MACE model behind `calc`: from the card when there is one, else
    from the loaded torch model. (None, None) when neither says."""
    arch = (card or {}).get("model", {}).get("architecture", {}) if card else {}
    r, n = arch.get("r_max"), arch.get("num_interactions")
    if r is None or n is None:
        try:
            mdl = calc.models[0]
            r = r if r is not None else float(mdl.r_max)
            n = n if n is not None else int(mdl.num_interactions)
        except Exception:
            pass
    return r, n


# ---- the calculation -------------------------------------------------------------

def run_phonon_pipeline(atoms_ref, calc, args, f_out, tag, noise, quiet=False):
    """One model, one volume: optional pre-relax (positions only), symmetry-reduced displacements,
    forces, force constants (residual forces of the undisplaced supercell subtracted, acoustic sum
    rule imposed), band structure, mesh (DOS, PDOS, thermal properties). Prints its own diagnostic
    lines unless `quiet`. Returns a dict with everything the report, the comparison and the plots
    need ('bs' & co. absent for a 0D system)."""
    from stb.core import mace_phonons
    from stb.phonons_pos import build_band_path, band_path_to_phonopy_format

    def say(text):
        if not quiet:
            print_dual(text, f_out)

    t0 = time.monotonic()
    atoms = atoms_ref.copy()
    atoms.calc = calc
    with capture_library_noise(noise, tag):
        f0 = float(np.linalg.norm(atoms.get_forces(), axis=1).max())
        phonon, relaxed, relax_info = mace_phonons.generate_ml_displacements(
            atoms, args.dim, args.distance, calc, relax=args.relax, fmax=args.fmax)
    res = {"tag": tag, "phonon": phonon, "f0": f0, "relax_info": relax_info,
           "n_disp": len(phonon.supercells_with_displacements)}
    shift = float(np.abs(relaxed.positions - atoms_ref.positions).max())
    res["relax_shift"] = shift

    sym = phonon.symmetry
    res["spacegroup"] = sym.get_international_table() or "unknown"
    res["pointgroup"] = sym.pointgroup_symbol
    res["n_ops"] = len(sym.symmetry_operations['rotations'])
    n_naive = phonon.dataset['natom'] * 6

    say(f"Input forces      : max |F| {f0:.4f} eV/Ang")
    if relax_info is not None:
        _, f1, converged, steps = relax_info
        say(f"Pre-relax         : max |F| {f1:.4f} eV/Ang after {steps} step(s) "
            f"({'converged' if converged else 'NOT converged: hit the step cap'}, --fmax {args.fmax}); "
            f"largest atomic shift {shift:.4f} Ang")
    else:
        say("Pre-relax         : skipped (--no-relax)")
    say(f"Symmetry          : {res['spacegroup']} (point group {res['pointgroup']}, {res['n_ops']} operations "
        "in the supercell)")
    say(f"Displacements     : {res['n_disp']} of {args.distance} Ang (vs. {n_naive} without symmetry: "
        f"{100 * (1 - res['n_disp'] / n_naive):.1f} % fewer evaluations)")

    with capture_library_noise(noise, tag):
        res["res_max"] = mace_phonons.compute_force_constants(phonon, calc)
        fc = phonon.force_constants
        res["asr_drift"] = float(np.abs(fc.sum(axis=1)).max())
        phonon.symmetrize_force_constants(show_drift=False)
    res["fc_seconds"] = time.monotonic() - t0
    say(f"Residual forces   : max {res['res_max']:.2e} eV/Ang per atom in the undisplaced supercell -- "
        "subtracted from every displacement's")
    say(f"Acoustic sum rule : largest violation {res['asr_drift']:.2e} eV/Ang^2 before symmetrization "
        "(imposed afterwards)")
    say(f"Force constants   : {res['fc_seconds']:.1f} s ({tag})")

    path = build_band_path(phonon.primitive, 1.0, args.vacuum_gap)
    if path is not None:
        kpoints_dict, path_segments, bravais_name = path
        bands, band_labels, path_connections = band_path_to_phonopy_format(
            kpoints_dict, path_segments, args.band_points)
        with capture_library_noise(noise, tag):
            phonon.run_band_structure(bands, with_group_velocities=False, path_connections=path_connections,
                                      labels=band_labels, is_legacy_plot=False)
        res.update(bs=phonon.band_structure, band_labels=band_labels, path_connections=path_connections,
                   bravais_name=bravais_name, path_segments=path_segments, kpoints=kpoints_dict,
                   band_freqs=np.concatenate(phonon.band_structure.frequencies))

    # with_eigenvectors + no mesh symmetry: the projected DOS needs per-atom character
    with capture_library_noise(noise, tag):
        phonon.run_mesh(args.mesh, with_eigenvectors=True, is_mesh_symmetry=False)
        phonon.run_total_dos()
        dos = phonon.get_total_dos_dict()
        phonon.run_projected_dos()
        proj = phonon.get_projected_dos_dict()['projected_dos']
        phonon.run_thermal_properties(t_min=args.tmin, t_max=args.tmax, t_step=args.tstep)
        tp = phonon.get_thermal_properties_dict()
    res["mesh_freqs"] = np.array(phonon.mesh.frequencies)
    res["mesh_q"] = np.array(phonon.mesh.qpoints)
    res["dos_freqs"], res["dos"] = dos['frequency_points'], dos['total_dos']
    species_pdos = {}
    for i, s_ in enumerate(phonon.primitive.symbols):
        species_pdos[s_] = species_pdos.get(s_, np.zeros_like(res["dos_freqs"])) + proj[i]
    res["species_pdos"] = species_pdos
    res["temps"], res["free_energy"] = tp['temperatures'], tp['free_energy']
    res["entropy"], res["heat_capacity"] = tp['entropy'], tp['heat_capacity']
    zpe = getattr(phonon.thermal_properties, "zero_point_energy", None)
    res["zpe_kjmol"] = float(zpe) if zpe is not None else None
    res["n_prim"] = len(phonon.primitive)
    res["seconds"] = time.monotonic() - t0
    return res


def freqs_at(phonon, q):
    return np.sort(np.asarray(phonon.get_frequencies(np.asarray(q, float))))


MAX_LISTED_BRANCHES = 18


def freq_list(f, spec="{:6.2f}"):
    """Frequencies as text; a large cell lists the lowest 9 and the highest 6 branches."""
    f = list(f)
    if len(f) <= MAX_LISTED_BRANCHES:
        return " ".join(spec.format(x) for x in f)
    return (" ".join(spec.format(x) for x in f[:9]) + f" ... ({len(f) - 15} more) ... "
            + " ".join(spec.format(x) for x in f[-6:]))


def special_points(res):
    """Special points in the order the path visits them: [(label, frac)]."""
    seen = []
    for seg in res.get("path_segments") or []:
        for lab in seg:
            if lab not in [s for s, _ in seen]:
                seen.append((lab, res["kpoints"][lab]))
    return seen


def gamma_irreps(phonon, noise):
    """[(label or None, [band indices], frequency)] at Gamma, from phonopy's irreducible
    representations; None when phonopy cannot label this point group."""
    try:
        with capture_library_noise(noise, "irreps"):
            phonon.set_irreps([0, 0, 0], degeneracy_tolerance=1e-3)
        ir = phonon.irreps
        labels = getattr(ir, "_ir_labels", None)
        sets = getattr(ir, "_degenerate_sets", None)
        freqs = np.asarray(ir.frequencies)
        if sets is None:
            return None
        return [((labels[i] if labels and i < len(labels) else None), list(s), float(freqs[s[0]]))
                for i, s in enumerate(sets)]
    except Exception:
        return None


def sound_velocities(res):
    """Acoustic slopes from Gamma toward every special point: [(label, [v1, v2, v3] m/s or None,
    [linear?])]. f = v |q| for a linear branch (|q| without 2 pi), so v = 100 f[THz]/|q|[1/Ang] m/s.
    A branch whose frequency grows 4x when |q| doubles is quadratic (a flexural ZA mode) and has
    no sound velocity."""
    phonon = res["phonon"]
    rec = np.linalg.inv(np.asarray(phonon.primitive.cell, float)).T     # rows: b_i without 2 pi
    out = []
    for lab, frac in special_points(res):
        frac = np.asarray(frac, float)
        qc = frac @ rec
        if np.linalg.norm(qc) < 1e-9:
            continue
        unit = frac / np.linalg.norm(qc)
        f1 = freqs_at(phonon, unit * SOUND_Q)[:3]
        f2 = freqs_at(phonon, unit * 2 * SOUND_Q)[:3]
        vel, lin = [], []
        for a, b in zip(f1, f2):
            linear = a > 1e-6 and b / a < 3.0
            lin.append(linear)
            vel.append(100.0 * a / SOUND_Q if linear else None)
        out.append((lab, vel, lin))
    return out


def mesh_stats(res):
    """Frequency statistics from the uniform mesh: imaginary modes, softest mode, Debye temperature
    from the second moment of the DOS (theta_D = h/k_B sqrt(5/3 <f^2>)), mean frequency and the
    gaps between bands."""
    f = res["mesh_freqs"]
    q = res["mesh_q"]
    gamma = np.all(np.abs(q) < 1e-9, axis=1)
    f_off = f.copy()
    f_off[gamma, :3] = np.nan            # the three acoustic modes at Gamma are zero by construction
    imag_q = int(np.nansum(np.nanmin(f_off, axis=1) < _IMAGINARY_MODE_TOL_THZ))
    i_min = np.unravel_index(np.nanargmin(f_off), f.shape)
    pos = f[f > 0]
    st = {"imag_q": imag_q, "n_q": len(q), "f_min": float(np.nanmin(f_off)),
          "q_min": q[i_min[0]], "band_min": int(i_min[1]), "f_max": float(f.max()),
          "f_mean": float(pos.mean()) if len(pos) else None,
          "debye_K": float(H_OVER_KB * np.sqrt(5.0 / 3.0 * np.mean(pos ** 2))) if len(pos) else None}
    gaps = []
    for b in range(f.shape[1] - 1):
        lo, hi = float(f[:, b].max()), float(f[:, b + 1].min())
        if hi - lo > 0.05:
            gaps.append((b, lo, hi))
    st["gaps"] = gaps
    return st


def thermal_rows(res):
    """[(T, F, S, Cv, U)] at REPORT_TEMPS inside the computed range (per primitive cell)."""
    t = np.asarray(res["temps"])
    rows = []
    for T in REPORT_TEMPS:
        if T < t[0] - 1e-9 or T > t[-1] + 1e-9:
            continue
        i = int(np.argmin(np.abs(t - T)))
        if abs(t[i] - T) > 0.5 * (t[1] - t[0] if len(t) > 1 else 1):
            continue
        F, S, C = res["free_energy"][i], res["entropy"][i], res["heat_capacity"][i]
        rows.append((float(t[i]), F, S, C, F + t[i] * S / 1000.0))
    return rows


def value_at(res, key, T):
    t = np.asarray(res["temps"])
    i = int(np.argmin(np.abs(t - T)))
    return float(res[key][i]) if abs(t[i] - T) <= max(1.0, 0.5 * (t[1] - t[0] if len(t) > 1 else 1)) else None


def _pbc_atomic_shift(positions_new, positions_ref, cell):
    """Per-atom displacement magnitude between two position arrays of the
    same periodic cell, using the minimum-image convention (same convention
    as core/md_traj.py's trajectory unwrapping elsewhere in this suite) --
    NOT a plain Cartesian difference. A plain difference is wrong here: an
    atom sitting near a cell face can have its wrapped fractional
    coordinate flip from ~0.999 to ~0.001 between the two supercells for a
    displacement of a fraction of an Angstrom, which shows up as a jump of
    nearly a full lattice vector (tens of Angstrom) in raw Cartesian
    coordinates -- verified live: an amplitude-calibration probe this
    caused to read ~74 Ang instead of the true ~0.0002 Ang for a real
    unstable-mode test case, silently breaking the requested-amplitude
    scaling below.
    """
    frac_diff = (positions_new - positions_ref) @ np.linalg.inv(cell)
    frac_diff -= np.round(frac_diff)
    cart_diff = frac_diff @ cell
    return np.linalg.norm(cart_diff, axis=1)


def freeze_unstable_mode(phonon, args, species_meta, f_out):
    """If the DOS mesh (already computed with eigenvectors on a full,
    non-symmetry-reduced grid inside run_phonon_pipeline) contains a
    genuine imaginary mode -- below the same numerical-noise tolerance used
    for the band-path warning -- displaces the supercell along that mode's
    eigenvector and writes it out as a new structure, the exact escape-a-
    saddle-point trick stb-phononsPos's own --freeze-unstable-mode already
    does for the SIESTA path (phonon.run_modulations, calibrated to the
    requested Angstrom amplitude via a unit-amplitude probe displacement
    first -- the modulation amplitude parameter isn't itself in Angstrom).
    Unlike phonons_pos.py, no bohr conversion is needed on write: this
    tool's Phonopy object is built directly in Angstrom (see
    core/mace_phonons.py's module docstring), so the modulated supercell's
    positions/cell are already in the right convention for from_pymatgen/
    write_fdf. Returns the written path, or None if the mesh had no mode
    below the noise tolerance.
    """
    freqs = phonon.mesh.frequencies
    min_freq = freqs.min()
    if min_freq >= _IMAGINARY_MODE_TOL_THZ:
        print_dual("No imaginary mode on the sampled mesh -- nothing to freeze.", f_out)
        return None

    q_idx, band_idx = np.unravel_index(np.argmin(freqs), freqs.shape)
    q_point = phonon.mesh.qpoints[q_idx]

    phonon.run_modulations(args.mesh, [[q_point, band_idx, 1.0, 0.0]])
    mod = phonon.modulation
    probe_shift = _pbc_atomic_shift(
        mod.modulated_supercells[0].positions, mod.supercell.positions,
        np.array(mod.supercell.cell)).max()
    scale = args.freeze_amplitude / probe_shift

    phonon.run_modulations(args.mesh, [[q_point, band_idx, scale, 0.0]])
    mod = phonon.modulation
    achieved_shift = _pbc_atomic_shift(
        mod.modulated_supercells[0].positions, mod.supercell.positions,
        np.array(mod.supercell.cell)).max()

    frozen = mod.modulated_supercells[0]
    frozen_atoms = Atoms(numbers=frozen.numbers, positions=frozen.positions,
                         cell=frozen.cell, pbc=True)
    pmg = AseAtomsAdaptor.get_structure(frozen_atoms)
    new_structure = structure_io.from_pymatgen(pmg, species_meta=species_meta)
    out_path = os.path.join(args.output_dir, "frozen_mode.fdf")
    structure_io.write_fdf(new_structure, out_path)

    print_dual(f"Softest mode      : q = {np.array2string(np.array(q_point), precision=4)}, "
               f"band {band_idx}, {min_freq:.4f} THz", f_out)
    print_dual(f"Displacement      : requested {args.freeze_amplitude:.4f} Ang, "
               f"achieved {achieved_shift:.4f} Ang (max atomic shift)", f_out)
    print_dual(f"Structure written : {out_path}", f_out)
    return out_path


ANIMATE_OPTICAL = "optical"


def _animate_mode_arg(value):
    """--animate-mode value: a band index (int >= 0) or 'optical' (the first optical mode)."""
    if value.strip().lower() == ANIMATE_OPTICAL:
        return ANIMATE_OPTICAL
    try:
        band = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected a band index or '{ANIMATE_OPTICAL}', got {value!r}")
    if band < 0:
        raise argparse.ArgumentTypeError(f"band index must be >= 0, got {band}")
    return band


def animate_normal_mode(phonon, qpoint, band_index, amplitude, n_frames, dim):
    """Builds a multi-frame ASE trajectory animating one normal mode's
    vibration at fractional q-point `qpoint`, by calling phonon.
    run_modulations once per frame with the phase swept over a full 360
    -degree cycle -- same phonopy modulation machinery freeze_unstable_mode
    uses for a single displaced snapshot, sampled here at many phases
    instead, purely for visualization (OVITO/VMD), not a physical
    perturbation. `dim` is the modulation supercell size and must make
    `qpoint` commensurate (q times dim ~ integer vector) -- Gamma is
    commensurate with any dim, [1, 1, 1] included; a non-Gamma q needs a
    correspondingly larger dim. phonon.get_frequencies(q) diagonalizes the
    dynamical matrix directly at the given q (no mesh/supercell needed for
    this part), so it works for any q regardless of `dim`, letting the
    out-of-range check below run before wasting time on a bad request.

    `amplitude` is calibrated to true Angstrom via the same unit-amplitude
    probe + rescale trick freeze_unstable_mode uses -- phonopy's own
    "amplitude" parameter is an internal eigenvector-normalization scale,
    not itself in Angstrom (verified live: amplitude=1.0 gave a ~0.094 Ang
    max shift for a real optical mode, not 1 Ang), so passing the user's
    requested Angstrom value straight through would silently produce the
    wrong displacement, just like the pre-fix freeze_unstable_mode bug this
    mirrors. The scale is amplitude-only (not phase-dependent, since phase
    is applied as a separate multiplicative factor in phonopy's own
    modulation formula), so one probe covers every frame in the cycle.

    Returns (frames, n_bands_at_q, frequency_thz) -- frames is None (with
    n_bands_at_q still valid) if band_index is out of range, so the caller
    can print a friendlier error.
    """
    frequencies_at_q = phonon.get_frequencies(np.array(qpoint))
    n_bands = len(frequencies_at_q)
    if not (0 <= band_index < n_bands):
        return None, n_bands, None

    phonon.run_modulations(dim, [[qpoint, band_index, 1.0, 0.0]])
    mod = phonon.modulation
    probe_shift = _pbc_atomic_shift(
        mod.modulated_supercells[0].positions, mod.supercell.positions,
        np.array(mod.supercell.cell)).max()
    scale = amplitude / probe_shift if probe_shift > 1e-12 else 0.0

    frames = []
    for phase_deg in np.linspace(0.0, 360.0, n_frames, endpoint=False):
        phonon.run_modulations(dim, [[qpoint, band_index, scale, phase_deg]])
        mod = phonon.modulation
        cell_atoms = mod.modulated_supercells[0]
        frames.append(Atoms(numbers=cell_atoms.numbers, positions=cell_atoms.positions,
                            cell=cell_atoms.cell, pbc=True))

    return frames, n_bands, frequencies_at_q[band_index]


# ---- figures (matplotlib) --------------------------------------------------------

def _ticks(res):
    from stb.phonons_pos import band_tick_positions, pretty_label
    ticks = band_tick_positions(res["bs"].distances, res["band_labels"], res["path_connections"])
    return [(t[0], pretty_label(t[1])) for t in ticks]


def _draw_bands(ax, series, title):
    for label, res, color in series:
        first = True
        for seg_d, seg_f in zip(res["bs"].distances, res["bs"].frequencies):
            for band in seg_f.T:
                ax.plot(seg_d, band, color=color, linewidth=0.9, label=label if first else None)
                first = False
    ax.axhline(0.0, color='gray', linestyle='--', linewidth=0.8)
    res0 = series[0][1]
    for x, _ in _ticks(res0):
        ax.axvline(x, color='gray', linewidth=0.5, alpha=0.5)
    ax.set_xticks([x for x, _ in _ticks(res0)])
    ax.set_xticklabels([lab for _, lab in _ticks(res0)])
    ax.set_xlim(res0["bs"].distances[0][0], res0["bs"].distances[-1][-1])
    ax.set_ylabel("Frequency (THz)")
    ax.set_title(title)
    sec = ax.secondary_yaxis('right', functions=(lambda f: f * THZ_TO_CM, lambda c: c / THZ_TO_CM))
    sec.set_ylabel("cm$^{-1}$")
    if len(series) > 1:
        ax.legend(fontsize=8, loc='lower right')


def plot_bands(series, out_path):
    """series: [(label, result, color)] on one shared q-path (same structure, same --band-points)."""
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(9, 5.2))
    _draw_bands(ax, series, "Phonon band structure" if len(series) == 1
                else "Phonon band structure: " + " vs. ".join(s[0] for s in series))
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    return fig


def plot_dos(series, out_path):
    """Total DOS of every series; the per-species PDOS for a single series."""
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    for label, res, color in series:
        ax.plot(res["dos_freqs"], res["dos"], label=f"{label} (total)" if len(series) > 1 else "Total",
                color=color, linewidth=1.5)
        if len(series) == 1 and len(res["species_pdos"]) > 1:
            for k, (s_, pdos) in enumerate(sorted(res["species_pdos"].items())):
                ax.plot(res["dos_freqs"], pdos, label=s_, linewidth=0.9, linestyle='--',
                        color=_SERIES_COLORS[(k + 2) % len(_SERIES_COLORS)])
    ax.set_xlabel("Frequency (THz)")
    ax.set_ylabel(f"DOS (states/THz per {series[0][1]['n_prim']}-atom primitive cell)")
    ax.set_title("Phonon density of states")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    return fig


def plot_thermal(series, out_path):
    """F, S, Cv per atom (so models and cells compare directly), with the Dulong-Petit limit."""
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4))
    for label, res, color in series:
        n = res["n_prim"]
        t = np.asarray(res["temps"])
        axes[0].plot(t, np.asarray(res["free_energy"]) / n * KJMOL_TO_MEV, label=label, color=color)
        axes[1].plot(t, np.asarray(res["entropy"]) / n / R_GAS, label=label, color=color)
        axes[2].plot(t, np.asarray(res["heat_capacity"]) / n / R_GAS, label=label, color=color)
    axes[2].axhline(3.0, color='gray', linestyle='--', linewidth=0.8, label="Dulong-Petit (3 k$_B$)")
    for ax, (yl, ti) in zip(axes, (("F$_{vib}$ (meV/atom)", "Vibrational free energy"),
                                   ("S (k$_B$/atom)", "Entropy"), ("C$_V$ (k$_B$/atom)", "Heat capacity"))):
        ax.set_xlabel("Temperature (K)")
        ax.set_ylabel(yl)
        ax.set_title(ti)
        ax.grid(alpha=0.3)
    axes[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    return fig


def plot_qha(q, out_path):
    """2x3: E(V) + F(V,T) curves, V(T), linear expansion, B(T), Cp and Cv, Grueneisen."""
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.5))
    ax = axes[0, 0]
    v, e = np.asarray(q["volumes"]), np.asarray(q["energies"])
    ax.plot(v, (e - e.min()) * 1000 / q["n_atoms"], 'o', color='k', label="E static (points)")
    for k, (T, vv, ff) in enumerate(q["fv_curves"]):
        c = _SERIES_COLORS[k % len(_SERIES_COLORS)]
        ax.plot(vv, ff - ff.min(), color=c, linewidth=1.2, label=f"F, {T:.0f} K")
        ax.plot(vv[np.argmin(ff)], 0.0, 'v', color=c, markersize=7)
    ax.set_xlabel("Volume (Ang$^3$/cell)")
    ax.set_ylabel("each curve minus its minimum (meV/atom)")
    ax.set_title("E(V) and F(V,T): the minimum moves with T")
    ax.legend(fontsize=7)
    panels = ((axes[0, 1], q["V_T"] / q["n_atoms"], "V (Ang$^3$/atom)", "Equilibrium volume"),
              (axes[0, 2], q["alpha_L"] * 1e6, r"$\alpha_L$ (10$^{-6}$ K$^{-1}$)", "Linear thermal expansion"),
              (axes[1, 0], q["B_T"], "B$_T$ (GPa)", "Isothermal bulk modulus"),
              (axes[1, 2], q["gamma_T"], r"$\gamma$", "Average Grueneisen parameter"))
    for ax, y, yl, ti in panels:
        ax.plot(q["T"], y, color='tab:blue')
        ax.set_xlabel("Temperature (K)")
        ax.set_ylabel(yl)
        ax.set_title(ti)
        ax.grid(alpha=0.3)
    axes[0, 2].axhline(0.0, color='gray', linewidth=0.6)
    ax = axes[1, 1]
    ax.plot(q["T"], q["Cp"] / q["n_atoms"] / R_GAS, color='tab:red', label="C$_P$")
    ax.plot(q["T"], q["Cv_V0"] / q["n_atoms"] / R_GAS, color='tab:blue', linestyle='--', label="C$_V$ at V$_0$")
    ax.axhline(3.0, color='gray', linestyle=':', linewidth=0.8)
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel("k$_B$/atom")
    ax.set_title("Heat capacity")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    return fig


def plot_convergence(rows, ref, out_path):
    """Bands of every supercell on the reference's path, band-path RMSE vs the largest, Cv/S."""
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6), gridspec_kw={"width_ratios": [2.2, 1, 1]})
    have_bs = "bs" in rows[-1]["res"]
    if have_bs:
        _draw_bands(axes[0], [(r["label"], r["res"], _SERIES_COLORS[k % len(_SERIES_COLORS)])
                              for k, r in enumerate(rows)], "Bands per supercell")
        axes[0].legend(fontsize=8, loc='lower right')
    natoms = [r["n_sc"] for r in rows]
    labels = [r["label"] for r in rows]
    axes[1].plot(natoms, [r["rmse"] for r in rows], 'o-')
    axes[1].set_ylabel(f"RMSE vs {rows[-1]['label']} (THz)")
    axes[1].set_title("Frequencies")
    axes[1].axhline(ref, color='gray', linestyle='--', linewidth=0.8)
    axes[2].plot(natoms, [r["cv300"] for r in rows], 'o-', color='tab:blue', label="C$_V$")
    axes[2].set_ylabel("C$_V$ at 300 K (k$_B$/atom)", color='tab:blue')
    ax2 = axes[2].twinx()
    ax2.plot(natoms, [r["s300"] for r in rows], 's-', color='tab:orange', label="S")
    ax2.set_ylabel("S at 300 K (k$_B$/atom)", color='tab:orange')
    axes[2].set_title("Thermal properties")
    for ax in axes[1:]:
        ax.set_xscale('log')
        ax.set_xticks(natoms)
        ax.set_xticklabels([f"{lab}\n{n}" for lab, n in zip(labels, natoms)], fontsize=8)
        ax.minorticks_off()
        ax.set_xlabel("supercell / atoms")
        ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    return fig


# ---- data files and gnuplot -------------------------------------------------------

def write_bands_dat(path, res):
    with open(path, "w") as f:
        f.write("# distance(1/Ang) freq_band1 freq_band2 ... (THz; blank line = segment break)\n")
        for seg_d, seg_f in zip(res["bs"].distances, res["bs"].frequencies):
            for d, row in zip(seg_d, seg_f):
                f.write(f"{d:.6f} " + " ".join(f"{x:.6f}" for x in row) + "\n")
            f.write("\n")


def write_dos_dat(path, res):
    cols, header = [res["dos_freqs"], res["dos"]], f"freq(THz) total_DOS(states/THz per {res['n_prim']}-atom cell)"
    for s_, pdos in sorted(res["species_pdos"].items()):
        cols.append(pdos)
        header += f" {s_}_PDOS"
    np.savetxt(path, np.column_stack(cols), header=header)


def write_thermal_dat(path, res):
    n = res["n_prim"]
    t, F, S, C = (np.asarray(res[k]) for k in ("temps", "free_energy", "entropy", "heat_capacity"))
    np.savetxt(path, np.column_stack([t, F, S, C, F + t * S / 1000.0, F / n * KJMOL_TO_MEV, S / n / R_GAS, C / n / R_GAS]),
               header=f"per {n}-atom primitive cell (mol = Avogadro's number of cells): T(K) F(kJ/mol) S(J/K/mol) "
                      "Cv(J/K/mol) U(kJ/mol) | per atom: F(meV) S(kB) Cv(kB)")


def _gp_label(lab):
    return "{/Symbol G}" if lab == "Γ" else lab.replace("_", "\\\\_")


def _gp_header(name, size="8,5"):
    return ['# --- STB Plot Configuration ---\n', f'# Generated by stb-mlphonons: gnuplot {name}.gplot\n',
            f'set terminal pdfcairo enhanced color font "Arial,10" size {size}\n', f'set output "{name}.pdf"\n',
            'set encoding utf8\nset grid\nset key font ",8"\n']


def write_gnuplot_bands(plot_dir, series):
    """bands.gplot over bands_<k>.dat (one per series)."""
    res0 = series[0][1]
    ticks = _ticks(res0)
    nb = res0["bs"].frequencies[0].shape[1]
    L = _gp_header("bands", "9,5.2")
    L.append('set xtics (' + ", ".join(f'"{_gp_label(lab)}" {x:.6f}' for x, lab in ticks) + ')\n')
    L.append('set grid xtics\nset ylabel "Frequency (THz)"\nset y2label "cm^{-1}"\nset ytics nomirror\nset y2tics 100 format "%.0f"\n'
             f'set link y2 via y*{THZ_TO_CM} inverse y/{THZ_TO_CM}\n'
             f'set xrange [{res0["bs"].distances[0][0]:.6f}:{res0["bs"].distances[-1][-1]:.6f}]\n'
             'set key bottom right\nset title "Phonon band structure"\n')
    parts = []
    for k, (label, _res, _c) in enumerate(series):
        fn = f"bands_{k}.dat"
        parts.append(f'for [i=2:{nb + 1}] "{fn}" using 1:i with lines lw 1.2 lc rgb "{_GP_COLORS[k % 6]}" '
                     f'title (i==2 ? "{label}" : "")')
    L.append("plot " + ", \\\n     ".join(parts) + ', 0 with lines dt 2 lc rgb "gray" notitle\n')
    for k, (_l, res, _c) in enumerate(series):
        write_bands_dat(os.path.join(plot_dir, f"bands_{k}.dat"), res)
    with open(os.path.join(plot_dir, "bands.gplot"), "w") as f:
        f.writelines(L)
    return "bands.gplot"


def write_gnuplot_dos(plot_dir, series):
    L = _gp_header("dos", "7.5,4.5")
    L.append('set xlabel "Frequency (THz)"\nset ylabel "DOS (states/THz/cell)"\nset title "Phonon density of states"\n')
    parts = []
    for k, (label, res, _c) in enumerate(series):
        fn = f"dos_{k}.dat"
        write_dos_dat(os.path.join(plot_dir, fn), res)
        parts.append(f'"{fn}" using 1:2 with lines lw 2 lc rgb "{_GP_COLORS[k % 6]}" title "{label} (total)"')
        if len(series) == 1 and len(res["species_pdos"]) > 1:
            for j, s_ in enumerate(sorted(res["species_pdos"])):
                parts.append(f'"{fn}" using 1:{j + 3} with lines lw 1.2 dt 2 lc rgb "{_GP_COLORS[(j + 2) % 6]}" '
                             f'title "{s_}"')
    L.append("plot " + ", \\\n     ".join(parts) + "\n")
    with open(os.path.join(plot_dir, "dos.gplot"), "w") as f:
        f.writelines(L)
    return "dos.gplot"


def write_gnuplot_thermal(plot_dir, series):
    L = _gp_header("thermal", "14,4.4")
    L.append('set multiplot layout 1,3\nset xlabel "Temperature (K)"\n')
    for col, yl, ti in ((6, "F_{vib} (meV/atom)", "Vibrational free energy"), (7, "S (k_B/atom)", "Entropy"),
                        (8, "C_V (k_B/atom)", "Heat capacity")):
        parts = [f'"thermal_{k}.dat" using 1:{col} with lines lw 2 lc rgb "{_GP_COLORS[k % 6]}" title "{lab}"'
                 for k, (lab, _r, _c) in enumerate(series)]
        if col == 8:
            parts.append('3 with lines dt 2 lc rgb "gray" title "Dulong-Petit"')
        key = {6: "top right", 7: "top left", 8: "bottom right"}[col]
        L.append(f'set key {key}\nset ylabel "{yl}"\nset title "{ti}"\nplot ' + ", ".join(parts) + "\n")
    L.append("unset multiplot\n")
    for k, (_l, res, _c) in enumerate(series):
        write_thermal_dat(os.path.join(plot_dir, f"thermal_{k}.dat"), res)
    with open(os.path.join(plot_dir, "thermal.gplot"), "w") as f:
        f.writelines(L)
    return "thermal.gplot"


def write_qha_dat(path, q):
    np.savetxt(path, np.column_stack([q["T"], q["V_T"], q["V_T"] / q["n_atoms"], q["alpha_V"], q["alpha_L"], q["B_T"],
                                      q["Cp"], q["Cv_V0"], q["gamma_T"], q["G_T"]]),
               header=f"input cell of {q['n_atoms']} atoms: T(K) V(Ang^3) V/atom(Ang^3) alpha_V(1/K) alpha_L(1/K) "
                      "B_T(GPa) Cp(J/K/mol) Cv_at_V0(J/K/mol) gamma G(eV)")


def write_gnuplot_qha(plot_dir, q):
    write_qha_dat(os.path.join(plot_dir, "qha.dat"), q)
    v, e = np.asarray(q["volumes"]), np.asarray(q["energies"])
    np.savetxt(os.path.join(plot_dir, "qha_ev.dat"), np.column_stack([v, (e - e.min()) * 1000 / q["n_atoms"]]),
               header="V(Ang^3/cell) E-Emin(meV/atom), static energies of the scan")
    with open(os.path.join(plot_dir, "qha_fv.dat"), "w") as f:
        f.write("# V(Ang^3/cell) F-min_V F (meV/atom); one block per temperature (index i = block)\n")
        for T, vv, ff in q["fv_curves"]:
            f.write(f"# T = {T:.0f} K\n")
            for a, b in zip(vv, ff - ff.min()):
                f.write(f"{a:.6f} {b:.6f}\n")
            f.write("\n\n")
    L = _gp_header("qha", "15,8.5")
    L.append('set multiplot layout 2,3\n')
    L.append('set xlabel "Volume (Ang^3/cell)"\nset ylabel "each curve minus its minimum (meV/atom)"\n'
             'set title "E(V) and F(V,T): the minimum moves with T"\nset key top center\n'
             'plot for [i=0:' + str(len(q["fv_curves"]) - 1) + '] "qha_fv.dat" index i using 1:2 with lines lw 1.5 '
             'title columnheader(1), "qha_ev.dat" using 1:2 with points pt 7 lc rgb "black" title "E static"\n')
    L.append('set xlabel "Temperature (K)"\n')
    for col, yl, ti in (("3", "V (Ang^3/atom)", "Equilibrium volume"),
                        ("($5*1e6)", "{/Symbol a}_L (10^{-6} K^{-1})", "Linear thermal expansion"),
                        ("6", "B_T (GPa)", "Isothermal bulk modulus")):
        L.append(f'set ylabel "{yl}"\nset title "{ti}"\nplot "qha.dat" using 1:{col} with lines lw 2 '
                 f'lc rgb "{_GP_COLORS[0]}" notitle\n')
    L.append(f'set key bottom right\nset ylabel "k_B/atom"\nset title "Heat capacity"\n'
             f'plot "qha.dat" using 1:($7/{q["n_atoms"]}/{R_GAS}) with lines lw 2 lc rgb "{_GP_COLORS[1]}" title "C_P", '
             f'"" using 1:($8/{q["n_atoms"]}/{R_GAS}) with lines lw 2 dt 2 lc rgb "{_GP_COLORS[0]}" title "C_V at V_0", '
             '3 with lines dt 3 lc rgb "gray" notitle\n')
    L.append(f'set ylabel "{{/Symbol g}}"\nset title "Average Grueneisen parameter"\n'
             f'plot "qha.dat" using 1:9 with lines lw 2 lc rgb "{_GP_COLORS[0]}" notitle\n')
    L.append("unset multiplot\n")
    temps = " ".join(f"{T:.0f}" for T, _v, _f in q["fv_curves"])
    L = [x.replace('title columnheader(1)', f'title sprintf("F, %s K", word("{temps}", i+1))') for x in L]
    with open(os.path.join(plot_dir, "qha.gplot"), "w") as f:
        f.writelines(L)
    return "qha.gplot"


def write_convergence_dat(path, rows):
    with open(path, "w") as f:
        f.write("# dim n_atoms_supercell n_disp max_freq(THz) Cv_300K(J/K/mol) S_300K(J/K/mol) "
                "rmse_vs_largest(THz) maxdiff_vs_largest(THz) Cv_300K(kB/atom) S_300K(kB/atom) ZPE(meV/atom)\n")
        for r in rows:
            f.write(f"{r['label']} {r['n_sc']} {r['res']['n_disp']} {r['fmax']:.6f} {r['cv300_mol']:.6f} "
                    f"{r['s300_mol']:.6f} {r['rmse']:.6f} {r['maxdiff']:.6f} {r['cv300']:.6f} {r['s300']:.6f} "
                    f"{fmt(r['zpe'], '{:.4f}', 'nan')}\n")


def write_gnuplot_convergence(plot_dir, rows):
    write_convergence_dat(os.path.join(plot_dir, "convergence.dat"), rows)
    L = _gp_header("convergence", "16,4.6")
    have_bs = "bs" in rows[-1]["res"]
    L.append('set multiplot layout 1,3\n' if have_bs else 'set multiplot layout 1,2\n')
    if have_bs:
        res0 = rows[-1]["res"]
        nb = res0["bs"].frequencies[0].shape[1]
        L.append('set xtics (' + ", ".join(f'"{_gp_label(lab)}" {x:.6f}' for x, lab in _ticks(res0)) + ')\n')
        L.append('set ylabel "Frequency (THz)"\nset title "Bands per supercell"\nset key bottom right\n')
        parts = []
        for k, r in enumerate(rows):
            fn = f"conv_bands_{r['label']}.dat"
            write_bands_dat(os.path.join(plot_dir, fn), r["res"])
            parts.append(f'for [i=2:{nb + 1}] "{fn}" using 1:i with lines lw 1.2 lc rgb "{_GP_COLORS[k % 6]}" '
                         f'title (i==2 ? "{r["label"]}" : "")')
        L.append("plot " + ", \\\n     ".join(parts) + "\n")
        L.append('set xtics autofreq\n')
    L.append('set logscale x\nset xlabel "atoms in the supercell"\nset key top right\n')
    L.append(f'set ylabel "RMSE vs {rows[-1]["label"]} (THz)"\nset title "Frequencies"\n'
             f'plot "convergence.dat" using 2:7 with linespoints pt 7 lw 2 lc rgb "{_GP_COLORS[0]}" notitle\n')
    L.append('set ylabel "C_V at 300 K (k_B/atom)"\nset y2label "S at 300 K (k_B/atom)"\nset ytics nomirror\n'
             'set y2tics\nset title "Thermal properties"\n'
             f'plot "convergence.dat" using 2:9 with linespoints pt 7 lw 2 lc rgb "{_GP_COLORS[0]}" title "C_V", '
             f'"" using 2:10 axes x1y2 with linespoints pt 5 lw 2 lc rgb "{_GP_COLORS[1]}" title "S"\n')
    L.append("unset multiplot\n")
    with open(os.path.join(plot_dir, "convergence.gplot"), "w") as f:
        f.writelines(L)
    return "convergence.gplot"




# ---- report pieces ---------------------------------------------------------------

def report_highsym(series, f_out, noise):
    """Frequencies at every special point of the path, for every series, and the irreps at Gamma."""
    res0 = series[0][1]
    pts = special_points(res0)
    nb = 3 * res0["n_prim"]
    print_dual(f"Frequencies in THz, sorted, {nb} branches per q-point "
               f"({res0['n_prim']}-atom primitive cell; reduced q in its reciprocal basis).", f_out)
    rows = []
    for lab, frac in pts:
        from stb.phonons_pos import pretty_label
        for k, (label, res, _c) in enumerate(series):
            f = freqs_at(res["phonon"], frac)
            rows.append(([pretty_label(lab) if k == 0 else "", "[" + " ".join(f"{x:.3f}" for x in frac) + "]" if k == 0 else "",
                          label, freq_list(f)], None))
    print_table(["q", "reduced q", "model", "frequencies (THz)"], rows, f_out)
    for label, res, _c in series:
        irr = gamma_irreps(res["phonon"], noise)
        if not irr:
            continue
        acoustic = {lab_ for lab_, idx, fr in irr if abs(fr) < 0.05 and lab_}
        print_dual(f"\nGamma-point modes ({label}): irreducible representations (phonopy), "
                   "IR-active = the representation of a vector (that of the acoustic modes)", f_out)
        rows = []
        for lab_, idx, fr in irr:
            kind = "acoustic" if abs(fr) < 0.05 else ("optical, IR-active" if lab_ in acoustic else "optical")
            rows.append(([lab_ or "?", str(len(idx)), f"{fr:8.3f}", f"{fr * THZ_TO_CM:8.1f}", f"{fr * THZ_TO_MEV:7.2f}",
                          kind], None))
        print_table(["irrep", "deg.", "THz", "cm^-1", "meV", "character"], rows, f_out)


def report_stability_and_bands(label, res, st, f_out, warn):
    from stb.phonons_pos import pretty_label
    if "bs" in res:
        bmin, bmax = float(res["band_freqs"].min()), float(res["band_freqs"].max())
        print_dual(f"[{label}] Band path : {bmin:.3f} .. {bmax:.3f} THz ({bmax * THZ_TO_CM:.1f} cm^-1)", f_out)
        if bmin < _IMAGINARY_MODE_TOL_THZ:
            print_dual(color_text(f"[{label}] Negative (imaginary) frequency found along the band path: "
                                  f"{bmin:.3f} THz", 'red'), f_out)
        else:
            print_dual(color_text(f"[{label}] No imaginary modes found along the band path "
                                  f"(minimum {bmin:.4f} THz, tolerance {_IMAGINARY_MODE_TOL_THZ} THz).", 'green'), f_out)
    print_dual(f"[{label}] Mesh {st['n_q']} q : {st['imag_q']} q-point(s) with an imaginary mode; softest "
               f"{st['f_min']:.3f} THz (band {st['band_min']}, q = [{' '.join(f'{x:.3f}' for x in st['q_min'])}])", f_out)
    if st["imag_q"]:
        warn(f"{label}: {st['imag_q']} of {st['n_q']} mesh q-points have an imaginary mode (softest {st['f_min']:.3f} THz "
             f"at q = [{' '.join(f'{x:.3f}' for x in st['q_min'])}]). Either the structure is not at a minimum of this "
             "model (tighter --fmax, relax the cell first), the supercell is too small (--check-convergence), "
             "or a genuine instability (--freeze-unstable-mode).")
    if st["gaps"]:
        for b, lo, hi in st["gaps"]:
            print_dual(f"[{label}] Gap       : {lo:.3f} .. {hi:.3f} THz ({hi - lo:.3f} THz) between branches "
                       f"{b + 1} and {b + 2}" + (" (acoustic-optical)" if b == 2 else ""), f_out)
    sv = sound_velocities(res) if "bs" in res else []
    if sv:
        rows = []
        for lab, vel, lin in sv:
            cells = [f"{v:7.0f}" if v is not None else "quadratic" for v in vel]
            rows.append(([f"Γ→{pretty_label(lab)}"] + cells, None))
        print_dual(f"[{label}] Sound velocities (m/s), slopes of the three acoustic branches at |q| = {SOUND_Q} 1/Ang:", f_out)
        print_table(["direction", "TA1", "TA2", "LA"], rows, f_out)


def report_dos(label, res, st, f_out):
    tot = np.trapezoid(res["dos"], res["dos_freqs"]) if hasattr(np, "trapezoid") else np.trapz(res["dos"], res["dos_freqs"])
    print_dual(f"[{label}] Max frequency {st['f_max']:.3f} THz ({st['f_max'] * THZ_TO_CM:.1f} cm^-1, "
               f"{st['f_max'] * THZ_TO_MEV:.1f} meV); mean {fmt(st['f_mean'], '{:.3f}')} THz; Debye temperature "
               f"{fmt(st['debye_K'], '{:.0f}')} K (second moment of the DOS)", f_out)
    print_dual(f"[{label}] DOS integrates to {tot:.2f} states (3N = {3 * res['n_prim']})", f_out)
    if res["species_pdos"]:
        rows = []
        f = res["dos_freqs"]
        integ = (lambda y: np.trapezoid(y, f)) if hasattr(np, "trapezoid") else (lambda y: np.trapz(y, f))
        for s_, p in sorted(res["species_pdos"].items()):
            w = integ(p)
            mean = integ(p * f) / w if w > 0 else float('nan')
            hi = f[np.nonzero(p > 0.05 * p.max())[0][-1]] if p.max() > 0 else float('nan')
            rows.append(([s_, f"{w:.2f}", f"{100 * w / tot:.1f} %", f"{mean:.2f}", f"{hi:.2f}"], None))
        print_dual(f"Species (PDOS)    : {', '.join(sorted(res['species_pdos']))}", f_out)
        print_table(["species", "states", "share", "mean f (THz)", "upper edge (THz)"], rows, f_out)


def report_thermal(series, f_out):
    for label, res, _c in series:
        n = res["n_prim"]
        print_dual(f"\n[{label}] per {n}-atom primitive cell (mol = mol of cells) and per atom:", f_out)
        if res.get("zpe_kjmol") is not None:
            print_dual(f"Zero-point energy : {res['zpe_kjmol']:.4f} kJ/mol = {res['zpe_kjmol'] / n * KJMOL_TO_MEV:.2f} "
                       "meV/atom", f_out)
        rows = []
        for T, F, S, C, U in thermal_rows(res):
            rows.append(([f"{T:.0f}", f"{F:.3f}", f"{S:.3f}", f"{C:.3f}", f"{U:.3f}",
                          f"{F / n * KJMOL_TO_MEV:.2f}", f"{S / n / R_GAS:.3f}", f"{C / n / R_GAS:.3f}",
                          f"{100 * C / n / R_GAS / 3:.0f} %"], None))
        print_table(["T (K)", "F kJ/mol", "S J/K/mol", "Cv J/K/mol", "U kJ/mol", "F meV/at", "S kB/at", "Cv kB/at",
                     "Cv/3kB"], rows, f_out)


def report_comparison(main, found, f_out):
    """Fine-tuned vs foundation on the same q-path, the same mesh and the same temperatures."""
    rows = []
    if "bs" in main and "bs" in found:
        d = main["band_freqs"] - found["band_freqs"]
        rows.append((["band path", f"{np.sqrt(np.mean(d ** 2)):.3f}", f"{np.abs(d).max():.3f}", f"{d.mean():+.3f}"], None))
    d = main["mesh_freqs"] - found["mesh_freqs"]
    rows.append(([f"mesh ({len(main['mesh_q'])} q)", f"{np.sqrt(np.mean(d ** 2)):.3f}", f"{np.abs(d).max():.3f}",
                  f"{d.mean():+.3f}"], None))
    nb = d.shape[1]
    if nb > 3:
        rows.append((["  acoustic branches", f"{np.sqrt(np.mean(d[:, :3] ** 2)):.3f}", f"{np.abs(d[:, :3]).max():.3f}",
                      f"{d[:, :3].mean():+.3f}"], None))
        rows.append((["  optical branches", f"{np.sqrt(np.mean(d[:, 3:] ** 2)):.3f}", f"{np.abs(d[:, 3:]).max():.3f}",
                      f"{d[:, 3:].mean():+.3f}"], None))
    print_table(["frequencies, fine-tuned - foundation", "RMSE THz", "max |d| THz", "mean THz"], rows, f_out)
    fm, ff = main["mesh_freqs"].max(), found["mesh_freqs"].max()
    print_dual(f"Highest frequency : fine-tuned {fm:.3f} THz, foundation {ff:.3f} THz ({100 * (ff / fm - 1):+.1f} %)", f_out)
    rows = []
    n = main["n_prim"]
    for T in (300.0, float(main["temps"][-1])):
        vals = [(value_at(r, "heat_capacity", T), value_at(r, "entropy", T), value_at(r, "free_energy", T)) for r in (main, found)]
        if None in vals[0] + vals[1]:
            continue
        rows.append(([f"{T:.0f}", f"{vals[0][0] / n / R_GAS:.3f} / {vals[1][0] / n / R_GAS:.3f}",
                      f"{vals[0][1] / n / R_GAS:.3f} / {vals[1][1] / n / R_GAS:.3f}",
                      f"{vals[0][2] / n * KJMOL_TO_MEV:.2f} / {vals[1][2] / n * KJMOL_TO_MEV:.2f}"], None))
    if rows:
        print_table(["T (K)", "Cv kB/atom (ft / found.)", "S kB/atom", "F meV/atom"], rows, f_out)
    print_dual("A large difference means the fine-tuning changed the curvature of the energy surface: the "
               "foundation model was never told about this level of theory.", f_out)


# ---- QHA ----------------------------------------------------------------------------

def run_qha(atoms_ref, calc, args, f_out, noise, warn, sections):
    """Quasi-harmonic approximation: the single-volume pipeline (no bands, no PDOS) at --n-volumes
    isotropically scaled volumes, then phonopy's PhonopyQHA on the INPUT cell: per-primitive-cell
    F/S/Cv multiplied by mace_phonons.cell_factor(). Also the mode Grueneisen parameters at the
    special points, from the frequencies of the scan."""
    from phonopy.api_qha import PhonopyQHA
    from stb.core import mace_phonons
    from stb.phonons_pos import build_band_path, pretty_label

    sections("VOLUME SCAN")
    n_atoms = len(atoms_ref)
    v0_input = atoms_ref.get_volume()
    strains = np.linspace(-abs(args.strain_range), abs(args.strain_range), args.n_volumes)
    # PhonopyQHA needs a few grid points beyond t_max for its finite differences
    t_max_grid = args.tmax + 3 * args.tstep
    print_dual(f"Volume range      : {-abs(args.strain_range):.1f} % to {abs(args.strain_range):.1f} % of the input "
               f"volume ({args.n_volumes} points, V_input {v0_input:.4f} Ang^3 = {v0_input / n_atoms:.4f} Ang^3/atom), "
               "isotropic, positions relaxed at each volume" + ("" if args.relax else " (--no-relax: not relaxed)"), f_out)

    volumes, energies, F, S, C, fmins, fmaxs, zpes, hs = [], [], [], [], [], [], [], [], []
    temps, factor, pts = None, None, None
    rows = []
    for i, strain in enumerate(strains):
        t0 = time.monotonic()
        a_i = atoms_ref.copy()
        a_i.set_cell(np.array(atoms_ref.get_cell()) * (1.0 + strain / 100.0) ** (1.0 / 3.0), scale_atoms=True)
        a_i.calc = calc
        with capture_library_noise(noise, f"QHA volume {i + 1}"):
            phonon, relaxed, _ = mace_phonons.generate_ml_displacements(
                a_i, args.dim, args.distance, calc, relax=args.relax, fmax=args.fmax)
            mace_phonons.compute_force_constants(phonon, calc)
            phonon.symmetrize_force_constants(show_drift=False)
            e_i = relaxed.get_potential_energy()
            v_i = relaxed.get_volume()
            phonon.run_mesh(args.mesh)
            phonon.run_thermal_properties(t_min=args.tmin, t_max=t_max_grid, t_step=args.tstep)
        tp = phonon.thermal_properties
        if factor is None:
            factor = mace_phonons.cell_factor(phonon)
            temps = np.asarray(tp.temperatures)
            path = build_band_path(phonon.primitive, 1.0, args.vacuum_gap)
            pts = [(lab, path[0][lab]) for lab in dict.fromkeys(l_ for seg in path[1] for l_ in seg)] if path else []
            print_dual(f"Cell factor       : {factor:g} (input cell {n_atoms} atoms / primitive {len(phonon.primitive)}): "
                       "phonopy's per-primitive-cell F, S, Cv are multiplied by it to match the volumes and energies",
                       f_out)
            print_dual(f"{'Strain%':>8} {'V (Ang^3)':>11} {'V/atom':>8} {'a/a_in':>8} {'E-Emin':>9} {'f_min':>8} "
                       f"{'f_max':>8} {'ZPE':>8} {'F_vib(300K)':>11} {'time':>6}", f_out)
            print_dual(f"{'':8} {'':11} {'':8} {'':8} {'meV/atom':>9} {'THz':>8} {'THz':>8} {'meV/at':>8} "
                       f"{'meV/atom':>11} {'s':>6}", f_out)
        f_i = np.asarray(phonon.mesh.frequencies)
        volumes.append(v_i)
        energies.append(e_i)
        F.append(np.asarray(tp.free_energy) * factor)
        S.append(np.asarray(tp.entropy) * factor)
        C.append(np.asarray(tp.heat_capacity) * factor)
        q_nonzero = np.any(np.abs(np.asarray(phonon.mesh.qpoints)) > 1e-9, axis=1)
        fmins.append(float(f_i[q_nonzero].min()) if q_nonzero.any() else float(f_i.min()))
        fmaxs.append(float(f_i.max()))
        zpes.append(float(getattr(tp, "zero_point_energy", np.nan)) * factor / n_atoms * KJMOL_TO_MEV)
        hs.append([freqs_at(phonon, frac) for _, frac in pts])
        vol_dir = os.path.join(args.output_dir, f"vol_{i:02d}")
        os.makedirs(vol_dir, exist_ok=True)
        with capture_library_noise(noise, f"QHA volume {i + 1}"):
            phonon.save(os.path.join(vol_dir, "phonopy_disp.yaml"), settings={'force_constants': True})
        rows.append((strain, v_i, e_i, time.monotonic() - t0))

    e_min = min(energies)
    i300 = int(np.argmin(np.abs(temps - 300.0)))
    for k, (strain, v_i, e_i, dt) in enumerate(rows):
        unstable = fmins[k] < _IMAGINARY_MODE_TOL_THZ
        line = (f"{strain:8.2f} {v_i:11.4f} {v_i / n_atoms:8.4f} {(v_i / v0_input) ** (1 / 3):8.5f} "
                f"{(e_i - e_min) * 1000 / n_atoms:9.3f} {fmins[k]:8.3f} {fmaxs[k]:8.3f} {zpes[k]:8.2f} "
                f"{F[k][i300] / n_atoms * KJMOL_TO_MEV:11.2f} {dt:6.1f}")
        print_dual(color_text(line, 'red') if unstable else line, f_out)
    if any(fm < _IMAGINARY_MODE_TOL_THZ for fm in fmins):
        warn("QHA: at least one volume has an imaginary mode on the mesh -- the quasi-harmonic free energy assumes "
             "every sampled volume is dynamically stable; narrow --strain-range.")
    i_emin = int(np.argmin(energies))
    if i_emin in (0, len(energies) - 1):
        warn(f"QHA: the lowest static energy is at the edge of the scan ({strains[i_emin]:+.1f} %): the model's "
             "equilibrium volume is outside or at the border of the sampled range. Relax the cell with this model "
             "(stb-mlrelax --relax-cell) or widen --strain-range.")

    with capture_library_noise(noise, "QHA fit"):
        qha = PhonopyQHA(volumes=volumes, electronic_energies=energies, temperatures=temps,
                         free_energy=np.array(F).T, cv=np.array(C).T, entropy=np.array(S).T,
                         eos=args.eos, t_max=args.tmax)
        e0, b0, b0p, v0 = qha.get_bulk_modulus_parameters()
        V_T = np.asarray(qha.volume_temperature)
        n_t = min(len(V_T), len(qha.thermal_expansion), len(qha.bulk_modulus_temperature),
                  len(qha.heat_capacity_P_polyfit), len(qha.gruneisen_temperature), len(qha.gibbs_temperature))
        T = temps[:n_t]
        q = {"T": T, "V_T": V_T[:n_t], "alpha_V": np.asarray(qha.thermal_expansion)[:n_t],
             "B_T": np.asarray(qha.bulk_modulus_temperature)[:n_t], "Cp": np.asarray(qha.heat_capacity_P_polyfit)[:n_t],
             "gamma_T": np.asarray(qha.gruneisen_temperature)[:n_t], "G_T": np.asarray(qha.gibbs_temperature)[:n_t]}
    q["alpha_L"] = q["alpha_V"] / 3.0
    q["Cv_V0"] = np.array([np.interp(v0, volumes, np.array(C)[:, k]) for k in range(n_t)])
    q.update(volumes=volumes, energies=energies, n_atoms=n_atoms, e0=e0, v0=v0, b0_gpa=b0 * EV_A3_TO_GPA, b0p=b0p)
    vv = np.linspace(min(volumes), max(volumes), 60)
    q["fv_curves"] = []
    for Tc in (0, 300, 600, 1000, 1500, 2000):
        if Tc > T[-1]:
            continue
        k = int(np.argmin(np.abs(temps - Tc)))
        ftot = np.array(energies) + np.array(F)[:, k] / 96.485332        # kJ/mol -> eV
        coef = np.polyfit(volumes, ftot, 3)
        q["fv_curves"].append((float(temps[k]), vv, (np.polyval(coef, vv) - e_min) * 1000 / n_atoms))

    sections("EQUILIBRIUM FIT (T=0, static)")
    a_ratio = (v0 / v0_input) ** (1 / 3)
    print_dual(f"Equation of state : {args.eos}, {len(volumes)} volumes, static energies (no vibrations)", f_out)
    print_dual(f"E0 (min. energy)  : {e0:.6f} eV per cell ({e0 / n_atoms:.6f} eV/atom)", f_out)
    print_dual(f"V0 (eq. volume)   : {v0:.4f} Ang^3 per cell = {v0 / n_atoms:.4f} Ang^3/atom "
               f"({100 * (v0 / v0_input - 1):+.3f} % vs the input; lattice x {a_ratio:.5f})", f_out)
    print_dual(f"B0 (bulk modulus) : {q['b0_gpa']:.2f} GPa", f_out)
    print_dual(f"B0' (pressure derivative): {b0p:.3f}", f_out)
    if abs(v0 / v0_input - 1) > 0.01:
        warn(f"QHA: the model's static equilibrium volume is {100 * (v0 / v0_input - 1):+.2f} % from the input structure. "
             "The scan is centred on the input; relax the cell with the same model first for a scan centred on V0.")

    sections("THERMAL EXPANSION SUMMARY")
    print_dual(f"Zero-point motion : V(0 K) {q['V_T'][0] / n_atoms:.4f} Ang^3/atom vs static {v0 / n_atoms:.4f} "
               f"({100 * (q['V_T'][0] / v0 - 1):+.3f} %)", f_out)
    rows_t = []
    for Tr in (0, 50, 100, 200, 300, 500, 800, 1000, 1500, 2000):
        if Tr > T[-1]:
            continue
        k = int(np.argmin(np.abs(T - Tr)))
        rows_t.append(([f"{T[k]:.0f}", f"{q['V_T'][k] / n_atoms:.4f}", f"{(q['V_T'][k] / v0) ** (1 / 3):.5f}",
                        f"{q['alpha_V'][k] * 1e6:.3f}", f"{q['alpha_L'][k] * 1e6:.3f}", f"{q['B_T'][k]:.2f}",
                        f"{q['Cp'][k] / n_atoms / R_GAS:.3f}", f"{q['gamma_T'][k]:.3f}",
                        f"{(q['G_T'][k] - e_min) * 1000 / n_atoms:.2f}"], None))
    print_table(["T (K)", "V Ang^3/at", "a/a0", "alpha_V 1e-6/K", "alpha_L 1e-6/K", "B_T GPa", "Cp kB/at", "gamma",
                 "G-Emin meV/at"], rows_t, f_out)
    print_dual("Thermal expansion : alpha_V = (1/V) dV/dT; linear alpha_L = alpha_V / 3 (isotropic scan).", f_out)
    for Tr in (300.0, float(T[-1])):
        k = int(np.argmin(np.abs(T - Tr)))
        print_dual(f"At T = {T[k]:.0f} K:", f_out)
        print_dual(f"  Volume            : {q['V_T'][k]:.4f} Ang^3 per cell", f_out)
        print_dual(f"  Thermal expansion : {q['alpha_V'][k]:.4e} 1/K volumetric, {q['alpha_L'][k]:.4e} 1/K linear", f_out)
        print_dual(f"  Bulk modulus      : {q['B_T'][k]:.2f} GPa", f_out)
        print_dual(f"  Gibbs free energy : {q['G_T'][k]:.6f} eV per cell", f_out)
    neg = T[q["alpha_V"] < 0]
    if len(neg):
        print_dual(f"Negative expansion: below {neg.max():.0f} K (minimum alpha_L {q['alpha_L'].min() * 1e6:.3f} 1e-6/K at "
                   f"{T[np.argmin(q['alpha_L'])]:.0f} K) -- modes that soften under compression (gamma < 0, e.g. transverse "
                   "acoustic) win at low temperature", f_out)
    v_lo, v_hi = min(volumes), max(volumes)
    out = T[(q["V_T"] < v_lo) | (q["V_T"] > v_hi)]
    if len(out):
        warn(f"QHA: V(T) leaves the scanned volume range above {out.min():.0f} K -- the free energy is extrapolated "
             "there; widen --strain-range.")
    elif q["V_T"][-1] > v_lo + 0.85 * (v_hi - v_lo):
        warn(f"QHA: at {T[-1]:.0f} K V(T) is in the top 15 % of the scanned range; widen --strain-range to keep the "
             "fit away from the edge.")

    if pts:
        sections("MODE GRUENEISEN PARAMETERS")
        print_dual("gamma = -d ln f / d ln V, from a straight-line fit of ln f against ln V over the scan; "
                   "acoustic modes at Gamma are left out.", f_out)
        lnv = np.log(volumes)
        rows_g = []
        for j, (lab, frac) in enumerate(pts):
            fr = np.array([h[j] for h in hs])                # (n_vol, n_band)
            vals = []
            for b in range(fr.shape[1]):
                if np.any(fr[:, b] <= 0.05):
                    vals.append("-")
                    continue
                vals.append(f"{-np.polyfit(lnv, np.log(fr[:, b]), 1)[0]:+.2f}")
            mid = fr[len(fr) // 2]
            rows_g.append(([pretty_label(lab), freq_list(mid), freq_list(vals, "{:>5}")], None))
        print_table(["q", "f at the middle volume (THz)", "gamma per branch"], rows_g, f_out)
        q["mode_gruneisen"] = rows_g
    return q


# ---- convergence ------------------------------------------------------------------

def run_convergence(atoms_ref, calc, args, conv_dims, f_out, noise, warn, sections):
    """The single-volume pipeline per supercell, compared on the largest supercell's q-path and at
    the special points (Gamma is commensurate with every supercell, so it says nothing: the band
    path and the zone-boundary points do)."""
    sections("SUPERCELL CONVERGENCE CHECK")
    print_dual(f"Supercells        : {', '.join('x'.join(map(str, d)) for d in conv_dims)}; reference = the largest; "
               f"tolerance {args.convergence_tol} THz (band-path RMSE)", f_out)
    rows = []
    for dim in conv_dims:
        dargs = copy.copy(args)
        dargs.dim = list(dim)
        label = "x".join(str(d) for d in dim)
        print(f"  ... supercell {label}", flush=True)
        res = run_phonon_pipeline(atoms_ref, calc, dargs, f_out, label, noise, quiet=True)
        sc = res["phonon"].supercell
        n = res["n_prim"]
        r = {"label": label, "dim": dim, "res": res, "n_sc": len(sc), "width": min(perpendicular_widths(sc.cell)),
             "fmax": float(res["mesh_freqs"].max()), "zpe": res["zpe_kjmol"] / n * KJMOL_TO_MEV if res["zpe_kjmol"] else None}
        r["cv300_mol"] = value_at(res, "heat_capacity", 300.0) or float('nan')
        r["s300_mol"] = value_at(res, "entropy", 300.0) or float('nan')
        r["cv300"], r["s300"] = r["cv300_mol"] / n / R_GAS, r["s300_mol"] / n / R_GAS
        rows.append(r)
    ref = rows[-1]["res"]
    for r in rows:
        if "bs" in ref:
            d = r["res"]["band_freqs"] - ref["band_freqs"]
        else:
            d = r["res"]["mesh_freqs"] - ref["mesh_freqs"]
        r["rmse"], r["maxdiff"] = float(np.sqrt(np.mean(d ** 2))), float(np.abs(d).max())
        r["imag"] = mesh_stats(r["res"])["imag_q"]

    table = [([r["label"], str(r["n_sc"]), str(r["res"]["n_disp"]), f"{r['width']:.2f}", f"{r['res']['seconds']:.1f}",
               f"{r['rmse']:.3f}" if r is not rows[-1] else "ref", f"{r['maxdiff']:.3f}" if r is not rows[-1] else "ref",
               fmt(r["zpe"], "{:.2f}"), f"{r['cv300']:.4f}", f"{r['s300']:.4f}", str(r["imag"])],
              'yellow' if r is not rows[-1] and r["rmse"] > args.convergence_tol else None) for r in rows]
    print_table(["supercell", "atoms", "disp.", "width Ang", "time s", "RMSE THz", "max|d| THz", "ZPE meV/at",
                 "Cv kB/at 300K", "S kB/at 300K", "imag. q"], table, f_out)
    print_dual("width = smallest distance between opposite faces of the supercell; force constants beyond half of it "
               "are folded back onto nearer images.", f_out)

    pts = special_points(ref) if "bs" in ref else []
    if pts:
        from stb.phonons_pos import pretty_label
        print_dual("\nFrequencies at the special points (THz), per supercell:", f_out)
        trows = []
        for lab, frac in pts:
            for k, r in enumerate(rows):
                fr = freqs_at(r["res"]["phonon"], frac)
                d = fr - freqs_at(ref["phonon"], frac)
                trows.append(([pretty_label(lab) if k == 0 else "", r["label"], freq_list(fr),
                               f"{np.abs(d).max():.3f}" if r is not rows[-1] else "ref"], None))
        print_table(["q", "supercell", "frequencies", "max|d|"], trows, f_out)

    ok = [r for r in rows[:-1] if r["rmse"] <= args.convergence_tol]
    if len(rows) > 1:
        if ok:
            first = ok[0]
            print_dual(color_text(f"Converged         : {first['label']} ({first['n_sc']} atoms) is within "
                                  f"{args.convergence_tol} THz of {rows[-1]['label']} (RMSE {first['rmse']:.3f} THz)", 'green'), f_out)
        else:
            warn(f"Convergence: no supercell smaller than {rows[-1]['label']} is within {args.convergence_tol} THz of it "
                 f"(best RMSE {min(r['rmse'] for r in rows[:-1]):.3f} THz) -- the largest one is not proven converged "
                 "either: add a larger supercell.")
        if rows[-2]["rmse"] > args.convergence_tol:
            print_dual(f"Last step         : {rows[-2]['label']} -> {rows[-1]['label']} still changes the bands by "
                       f"{rows[-2]['rmse']:.3f} THz RMSE", f_out)
    return rows


# ---- main -------------------------------------------------------------------------

def build_parser():
    parser = argparse.ArgumentParser(
        description="Standalone phonon calculation driven entirely by a MACE potential -- the "
                     "MACE-MP-0 foundation model, or a custom model fine-tuned on your own "
                     "SIESTA data via stb-mlffAnalysis (--custom-model). One shot: displacements, "
                     "force constants, band structure, DOS (total + per-species), and thermal "
                     "properties -- or, with --qha, a quasi-harmonic volume scan (thermal "
                     "expansion, bulk modulus vs. T). No SIESTA and no separate analysis tool "
                     "needed. A fast heuristic preview, not a substitute for the real SIESTA-"
                     "based stb-phononsCreate/stb-phononsPos workflow -- ML potentials are "
                     "markedly less reliable for phonons (imaginary modes especially) than DFT.",
        epilog="Example usage:\n"
               "  stb-mlphonons --file structure.fdf --dim 2 2 2\n"
               "  stb-mlphonons --file structure.fdf --custom-model mlff_model.model --custom-model-card mlff_model_card.json --dim 3 3 3\n"
               "  stb-mlphonons --file structure.fdf --qha --n-volumes 7 --strain-range 5 --save-gnuplot\n"
               "  stb-mlphonons --file structure.fdf --check-convergence --convergence-dims 1 1 1 2 2 2 3 3 3 --show",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--file", required=True, help="Input structure (.fdf).")
    parser.add_argument("--model", choices=["small", "medium", "large"], default="medium",
                        help="MACE-MP-0 model size (default: medium). Ignored if --custom-model is given.")
    parser.add_argument("--custom-model", default=None, metavar="PATH",
                        help="Use a custom MACE model file instead of the MACE-MP-0 foundation "
                             "potential -- e.g. one fine-tuned via stb-mlffAnalysis. Overrides --model.")
    model_card.add_card_argument(parser, "--custom-model")
    parser.add_argument("--skip-foundation-comparison", action="store_true",
                        help="If --custom-model is given, also evaluate the same phonon "
                             "pipeline with the raw (non-fine-tuned) foundation model "
                             "(--model) and overlay both on the bands/DOS/thermal plots, "
                             "to show whether fine-tuning actually changed the phonon "
                             "prediction. On by default; pass this to skip it (--qha never "
                             "does this comparison, regardless of this flag).")
    parser.add_argument("--dim", type=int, nargs=3, default=[2, 2, 2], metavar=("NA", "NB", "NC"),
                        help="Supercell dimensions for the finite-displacement method (default: 2 2 2).")
    parser.add_argument("--distance", type=float, default=0.01,
                        help="Displacement distance, Angstrom (default: 0.01).")
    parser.add_argument("--no-relax", dest="relax", action="store_false",
                        help="Skip the MACE pre-relax (positions only) before generating "
                             "displacements. On by default -- a phonon calculation assumes "
                             "~zero net force at the reference geometry; skipping this on an "
                             "unrelaxed structure is a common cause of spurious imaginary modes.")
    parser.add_argument("--fmax", type=float, default=0.05,
                        help="Force convergence threshold for the pre-relax, eV/Ang (default: 0.05).")
    parser.add_argument("--vacuum-gap", type=float, default=10.0,
                        help="Gap (Ang) used to detect vacuum-padded axes (default: 10.0).")
    parser.add_argument("--band-points", type=int, default=101,
                        help="q-points per band-structure segment (default: 101).")
    parser.add_argument("--mesh", type=int, nargs=3, default=[20, 20, 20], metavar=("NA", "NB", "NC"),
                        help="q-point mesh for DOS/thermal properties (default: 20 20 20).")
    parser.add_argument("--tmin", type=float, default=0.0, help="Min temperature, K (default: 0).")
    parser.add_argument("--tmax", type=float, default=1000.0, help="Max temperature, K (default: 1000).")
    parser.add_argument("--tstep", type=float, default=10.0, help="Temperature step, K (default: 10).")
    parser.add_argument("--qha", action="store_true",
                        help="Quasi-harmonic approximation instead of a single-volume "
                             "calculation: scans --n-volumes isotropically-scaled volumes and "
                             "fits an equation of state (thermal expansion, bulk modulus, Cp and "
                             "Grueneisen parameters vs. T, mode Grueneisen parameters). No "
                             "foundation-model comparison in this mode. Off by default.")
    parser.add_argument("--n-volumes", type=int, default=5,
                        help="Number of volume points for --qha, >= 4 (default: 5).")
    parser.add_argument("--strain-range", type=float, default=3.0,
                        help="Volume range around equilibrium for --qha, +/- %% (default: 3.0).")
    parser.add_argument("--eos", choices=["vinet", "murnaghan", "birch_murnaghan"], default="vinet",
                        help="Equation of state for --qha (default: vinet).")
    parser.add_argument("--check-convergence", action="store_true",
                        help="Instead of a single-volume calculation, run the same pipeline "
                             "at multiple supercell sizes (--convergence-dims) and compare them "
                             "with the largest one along the band path and at the special points "
                             "(RMSE, largest difference), plus heat capacity, entropy and ZPE -- "
                             "MACE is cheap enough to check --dim convergence directly, unlike "
                             "DFT. Not supported together with --qha. Off by default.")
    parser.add_argument("--convergence-dims", type=int, nargs='+', default=None, metavar="N",
                        help="Flat list of supercell sizes for --check-convergence, given as "
                             "consecutive NA NB NC triples (e.g. 1 1 1 2 2 2 3 3 3 tests "
                             "1x1x1, 2x2x2, and 3x3x3) -- must be a multiple of 3 values. "
                             "Default: --dim itself and --dim+1 on every non-vacuum axis "
                             "(two points).")
    parser.add_argument("--convergence-tol", type=float, default=0.1,
                        help="Band-path RMSE (THz) against the largest supercell under which a "
                             "smaller one counts as converged (default: 0.1).")
    parser.add_argument("--freeze-unstable-mode", action="store_true",
                        help="If the phonon mesh has a genuine imaginary (unstable) mode, "
                             "write a new structure displaced along that mode's eigenvector -- "
                             "a common trick to escape a saddle point and find the true "
                             "minimum. Same feature stb-phononsPos already has for the SIESTA "
                             "path. Not supported together with --qha/--check-convergence. "
                             "Off by default.")
    parser.add_argument("--freeze-amplitude", type=float, default=0.05,
                        help="Max atomic displacement (Angstrom) for --freeze-unstable-mode "
                             "(default: 0.05).")
    parser.add_argument("--animate-mode", type=_animate_mode_arg, nargs="?", const=ANIMATE_OPTICAL,
                        default=None, metavar="BAND_INDEX",
                        help="Animate a normal mode's vibration and write it as a multi-frame "
                             "trajectory (xsf/pdb/xyz, same convention as stb-ani2traj/"
                             "stb-mlmd) for viewing in OVITO/VMD -- BAND_INDEX selects which "
                             "mode (0 = lowest frequency) at --animate-qpoint; with no value, or "
                             "'optical', the first optical mode (band 3: bands 0-2 are acoustic, "
                             "at Gamma a rigid translation). Not supported "
                             "together with --qha/--check-convergence. Off by default.")
    parser.add_argument("--animate-qpoint", type=float, nargs=3, default=[0.0, 0.0, 0.0],
                        metavar=("QA", "QB", "QC"),
                        help="Fractional q-point to animate a mode at (default: 0 0 0, Gamma).")
    parser.add_argument("--animate-dim", type=int, nargs=3, default=[1, 1, 1],
                        metavar=("NA", "NB", "NC"),
                        help="Supercell size used to build the animated modulation -- must "
                             "make --animate-qpoint commensurate (q x dim ~ integer vector); "
                             "the default 1 1 1 only works exactly at Gamma (default: 1 1 1).")
    parser.add_argument("--animate-amplitude", type=float, default=0.5,
                        help="Max atomic displacement for the animation, Angstrom -- "
                             "deliberately larger than --freeze-amplitude since this is for "
                             "visualization, not a physical perturbation (default: 0.5).")
    parser.add_argument("--animate-frames", type=int, default=30,
                        help="Number of frames in one full vibration cycle (default: 30).")
    parser.add_argument("--animate-format", choices=sorted(OUTPUT_FORMATS), default="xsf",
                        help="Output trajectory format for --animate-mode, same convention as "
                             "stb-ani2traj/stb-mlmd (default: xsf).")
    parser.add_argument("--device", choices=["cpu", "cuda"], default="cpu",
                        help="Device to run the model on (default: cpu).")
    parser.add_argument("-o", "--output-dir", default="mlphonons_out",
                        help="Output directory for all plots/data (default: mlphonons_out).")
    parser.add_argument("--show", action="store_true",
                        help="Open the figures on screen at the end (needs a display; without one "
                             "it is a warning, not an error). The PNGs are written either way.")
    parser.add_argument("--save-gnuplot", action="store_true",
                        help=f"Write every figure as .dat + .gplot into <output-dir>/{PLOT_DIR}/ and run "
                             "gnuplot there (one PDF per figure) when it is installed. Off by default.")
    parser.add_argument("--save-data", action="store_true",
                        help="Also write the raw numeric data (.dat) behind each plot. Off by default.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-mlphonons {VERSION}")
    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()
    card = model_card.setup_from_args(args.custom_model, args.custom_model_card, "--custom-model", error=parser.error)

    if args.qha and args.n_volumes < 4:
        parser.error("--n-volumes must be at least 4 -- vinet/murnaghan/birch_murnaghan all "
                     "fit 4 parameters (E0, V0, B0, B0'), so need at least 4 volumes.")
    for flag_name, flag_val in (("--freeze-unstable-mode", args.freeze_unstable_mode),
                                ("--animate-mode", args.animate_mode is not None)):
        if args.qha and flag_val:
            parser.error(f"{flag_name} is not supported together with --qha.")
        if args.check_convergence and flag_val:
            parser.error(f"{flag_name} is not supported together with --check-convergence.")
    if args.qha and args.check_convergence:
        parser.error("--qha is not supported together with --check-convergence.")
    conv_dims_override = None
    if args.convergence_dims is not None:
        if len(args.convergence_dims) % 3 != 0:
            parser.error("--convergence-dims must be a multiple of 3 values (NA NB NC per point).")
        flat = args.convergence_dims
        conv_dims_override = [tuple(flat[i:i + 3]) for i in range(0, len(flat), 3)]
        if len(conv_dims_override) < 2:
            parser.error("--check-convergence needs at least 2 supercell sizes to compare.")

    import matplotlib
    if not args.show:
        matplotlib.use("Agg")

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - ML Phonons",
            "Bands, DOS, and thermal properties, driven entirely by a MACE potential",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    os.makedirs(args.output_dir, exist_ok=True)
    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None
    warns, noise, figures, files = [], [], [], []
    t_start = time.monotonic()
    sections = Sections(f_out)

    def warn(message):
        warns.append(message)

    def fail(message):
        print_dual(color_text(f"[FAIL] {message}", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    import importlib.util
    from stb.core.deps import require_mace
    if importlib.util.find_spec("mace") is None or importlib.util.find_spec("torch") is None:
        require_mace()          # prints the install hint and exits
    with capture_library_noise(noise, "import"):
        require_mace()
        from pymatgen.io.ase import AseAtomsAdaptor
        from stb.core import mace_relax
        import torch
        import ase
        import mace
        import phonopy
    model_card.collect_warnings(warns)

    print("\n" + color_text("ML PHONONS:", 'bold'))
    print("-" * 60)
    print_dual(color_text("===== STB-MLPHONONS REPORT =====", 'magenta'), f_out)

    mode_label = ("QHA (volume scan)" if args.qha else "supercell convergence check" if args.check_convergence
                  else "single-volume")
    compare = bool(args.custom_model) and not args.skip_foundation_comparison and not args.qha and not args.check_convergence
    sections("RUN METADATA")
    print_dual(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Software          : stb {VERSION}, phonopy {phonopy.__version__}, mace-torch "
               f"{getattr(mace, '__version__', '?')}, torch {torch.__version__}, ase {ase.__version__}", f_out)
    if args.device == "cuda":
        try:
            mace_relax.resolve_device("cuda")
            dev_label = f"cuda ({torch.cuda.get_device_name(0)})"
        except Exception as e:
            fail(f"--device cuda: {e}")
    else:
        dev_label = "cpu"
    print_dual(f"Device            : {dev_label}", f_out)
    print_dual(f"Input file        : {args.file}", f_out)
    print_dual(f"Mode              : {mode_label}", f_out)
    print_dual(f"Supercell dim     : {args.dim[0]} x {args.dim[1]} x {args.dim[2]}", f_out)
    print_dual(f"Displacement dist.: {args.distance} Ang", f_out)
    print_dual(f"q-mesh / path     : {' x '.join(map(str, args.mesh))} (DOS, thermal), {args.band_points} points per "
               "band segment", f_out)
    print_dual(f"Temperatures      : {args.tmin:g} .. {args.tmax:g} K, step {args.tstep:g}", f_out)
    print_dual(f"Model             : "
               f"{'custom (' + args.custom_model + ')' if args.custom_model else f'MACE-MP-0 ({args.model})'}", f_out)
    if compare:
        print_dual(f"Comparison        : also evaluating with foundation model ({args.model})", f_out)
    if report_path:
        print_dual(f"Report file       : {report_path}", f_out)

    # ---- structure
    sections("STRUCTURE")
    if not os.path.isfile(args.file):
        fail(f"Input file '{args.file}' not found.")
    if args.custom_model and not os.path.isfile(args.custom_model):
        fail(f"--custom-model file not found: {args.custom_model}")
    try:
        structure = structure_io.read_fdf(args.file)
    except Exception as e:
        fail(f"Could not read '{args.file}': {e}")
    atoms = AseAtomsAdaptor.get_atoms(structure_io.to_pymatgen(structure))
    n_atoms = len(atoms)
    vacuum_axes = kspace.detect_vacuum_axes(atoms.get_scaled_positions(), np.array(atoms.get_cell()), args.vacuum_gap)
    bulk = not any(vacuum_axes)
    cp = atoms.cell.cellpar()
    vol = atoms.get_volume()
    print_dual(f"[OK] Read {n_atoms} atom(s), dimensionality: {kspace.dimensionality_label(vacuum_axes)}", f_out)
    print_dual(f"Formula           : {atoms.get_chemical_formula(mode='hill')}", f_out)
    print_dual(f"Cell              : a {cp[0]:.4f}, b {cp[1]:.4f}, c {cp[2]:.4f} Ang; alpha {cp[3]:.2f}, beta {cp[4]:.2f}, "
               f"gamma {cp[5]:.2f} deg", f_out)
    print_dual(f"Volume            : {vol:.3f} Ang^3 ({vol / n_atoms:.4f} Ang^3/atom)"
               + (f", density {atoms.get_masses().sum() / vol * AMU_A3_TO_G_CM3:.4f} g/cm^3" if bulk else " (includes vacuum)"),
               f_out)
    d_min = model_card.shortest_distance(atoms)
    print_dual(f"Shortest distance : {fmt(d_min, '{:.4f}')} Ang", f_out)
    sc_cell = np.array(atoms.get_cell()) * np.array(args.dim)[:, None]
    widths = perpendicular_widths(sc_cell)
    print_dual(f"Supercell         : {' x '.join(map(str, args.dim))} -> {n_atoms * int(np.prod(args.dim))} atoms; "
               f"face-to-face widths {', '.join(f'{w:.2f}' for w in widths)} Ang", f_out)
    vac_dims = [ax for ax, v in zip('abc', vacuum_axes) if v and args.dim['abc'.index(ax)] > 1]
    if vac_dims and not args.check_convergence:
        warn(f"--dim repeats the vacuum-padded axis/axes {', '.join(vac_dims)}: that only multiplies the cost "
             "without adding periodicity -- use --dim 1 there.")
    periodic_widths = [w for w, v in zip(widths, vacuum_axes) if not v]
    if periodic_widths and min(periodic_widths) < SMALL_WIDTH and not args.check_convergence:
        warn(f"the supercell is only {min(periodic_widths):.2f} Ang thick along a periodic axis: force constants "
             "farther than half of that are folded onto nearer images. Check with --check-convergence.")

    # ---- model
    sections("MODEL")
    model_arg = args.custom_model if args.custom_model else args.model
    t0 = time.monotonic()
    with capture_library_noise(noise, "model"):
        calc = mace_relax.get_calculator(model=model_arg, device=args.device)
    if card is not None:
        model_card.print_model_card(card, args.custom_model_card, args.custom_model, f_out, warn)
        for msg in model_card.check_structure(card, atoms):
            warn(f"Model card check: {msg}")
    else:
        print_dual(f"Model             : MACE-MP-0 {args.model} (foundation, no card): trained on the Materials Project "
                   "PBE/PBE+U relaxations -- not on your SIESTA level of theory", f_out)
    r_max, n_int = model_range(calc, card)
    if r_max:
        print_dual(f"Model range       : r_max {r_max:g} Ang" + (f" x {n_int} interactions = receptive field "
                   f"{r_max * n_int:g} Ang (force constants can reach that far)" if n_int else ""), f_out)
    print_dual(f"Loaded in         : {time.monotonic() - t0:.1f} s", f_out)
    if not args.relax:
        warn("--no-relax: displacements are generated from the structure as given; if it is not at a minimum of "
             "this model, expect spurious imaginary modes.")

    label_main = "Fine-tuned" if args.custom_model else f"MACE-MP-0 ({args.model})"

    def finish(status_lines):
        """Plots on screen, gnuplot, summary, warnings, close."""
        sections("SUMMARY & FILES")
        for line in status_lines:
            print_dual(line, f_out)
        for kind, path in files:
            print_dual(f"{kind:<18}: {path}", f_out)
        if report_path:
            print_dual(f"Report            : {report_path}", f_out)
        not_asked = [f for f, on in (("--save-gnuplot", args.save_gnuplot), ("--save-data", args.save_data),
                                     ("--show", args.show), ("--save-report", args.save_report)) if not on]
        if not_asked:
            print_dual(f"Also available    : {' '.join(not_asked)}", f_out)
        print_dual(f"Wall time         : {time.monotonic() - t_start:.1f} s", f_out)
        if args.show and not display_available():
            warn("--show needs a display (none found, or a non-interactive matplotlib backend): the figures were "
                 "written as PNG files only.")
        sections("WARNINGS")
        if warns:
            for i, w in enumerate(warns, 1):
                print_dual(color_text(f"  {i}. {w}", 'yellow'), f_out)
        else:
            print_dual(color_text("  none", 'green'), f_out)
        if noise:
            print_dual(f"Library messages  : {len(noise)} distinct message(s) from mace/torch/phonopy captured, not shown", f_out)
        if f_out:
            f_out.close()
        print("\n" + "-" * 60)
        print(color_text("ML phonons complete -- a fast heuristic, not a substitute for DFT phonons.\n", 'bold'))
        import matplotlib.pyplot as plt
        if args.show and display_available() and figures:
            plt.show()
        for fig in figures:
            plt.close(fig)

    def gnuplot(writers):
        plot_dir = os.path.join(args.output_dir, PLOT_DIR)
        os.makedirs(plot_dir, exist_ok=True)
        scripts = [w() for w in writers]
        pdfs, msg = run_gnuplot(plot_dir, scripts)
        files.append(("gnuplot folder", f"{plot_dir}/ ({', '.join(scripts)})"))
        for p in pdfs:
            files.append(("gnuplot figure", p))
        if msg:
            warn(msg)

    # ---- convergence mode
    if args.check_convergence:
        if conv_dims_override is not None:
            conv_dims = conv_dims_override
        else:
            base = list(args.dim)
            conv_dims = [tuple(base), tuple(d if v else d + 1 for d, v in zip(base, vacuum_axes))]
        rows = run_convergence(atoms, calc, args, conv_dims, f_out, noise, warn, sections)
        conv_plot = os.path.join(args.output_dir, "convergence.png")
        figures.append(plot_convergence(rows, args.convergence_tol, conv_plot))
        files.append(("Convergence plot", conv_plot))
        if args.save_data:
            p = os.path.join(args.output_dir, "convergence.dat")
            write_convergence_dat(p, rows)
            files.append(("Data", p))
        if args.save_gnuplot:
            gnuplot([lambda: write_gnuplot_convergence(os.path.join(args.output_dir, PLOT_DIR), rows)])
        finish(["Status            : OK" if not warns else f"Status            : OK ({len(warns)} warning(s), see the "
                                                              "last section)"])
        return

    # ---- QHA mode
    if args.qha:
        q = run_qha(atoms, calc, args, f_out, noise, warn, sections)
        qha_plot = os.path.join(args.output_dir, "qha.png")
        figures.append(plot_qha(q, qha_plot))
        files.append(("QHA plot", qha_plot))
        if args.save_data:
            p = os.path.join(args.output_dir, "qha_thermal_expansion.dat")
            write_qha_dat(p, q)
            files.append(("Data", p))
        if args.save_gnuplot:
            gnuplot([lambda: write_gnuplot_qha(os.path.join(args.output_dir, PLOT_DIR), q)])
        k300 = int(np.argmin(np.abs(q["T"] - 300)))
        finish(["Status            : OK" + (f" ({len(warns)} warning(s), see the last section)" if warns else ""),
                f"Results           : V0 {q['v0'] / n_atoms:.4f} Ang^3/atom, B0 {q['b0_gpa']:.1f} GPa, "
                f"alpha_L({q['T'][k300]:.0f} K) {q['alpha_L'][k300] * 1e6:.2f} 1e-6/K, B_T {q['B_T'][k300]:.1f} GPa"])
        return

    # ---- single volume
    sections("DISPLACEMENTS & FORCE CONSTANTS")
    print_dual(f"--- {label_main}", f_out)
    main_res = run_phonon_pipeline(atoms, calc, args, f_out, label_main, noise)
    ph = main_res["phonon"]
    print_dual(f"Primitive cell    : {len(ph.primitive)} atoms (phonopy, spglib); supercell {len(ph.supercell)} atoms", f_out)
    yaml_path = os.path.join(args.output_dir, "phonopy_disp.yaml")
    with capture_library_noise(noise, "save"):
        ph.save(yaml_path, settings={'force_constants': True})
    files.append(("Force constants", yaml_path))
    print_dual(f"[OK] Force constants (ASR-corrected) saved to {yaml_path}", f_out)
    found_res = None
    if compare:
        print_dual(f"\n--- Foundation ({args.model})", f_out)
        with capture_library_noise(noise, "foundation model"):
            found_calc = mace_relax.get_calculator(model=args.model, device=args.device)
        found_res = run_phonon_pipeline(atoms, found_calc, args, f_out, f"Foundation ({args.model})", noise)

    series = [(label_main, main_res, _SERIES_COLORS[0])]
    if found_res is not None:
        series.append((f"Foundation ({args.model})", found_res, _SERIES_COLORS[1]))
    stats = {lab: mesh_stats(res) for lab, res, _c in series}

    if "bs" in main_res:
        sections("FREQUENCIES AT HIGH-SYMMETRY POINTS")
        report_highsym(series, f_out, noise)

    sections("BAND STRUCTURE & STABILITY")
    if "bs" not in main_res:
        print_dual(color_text("Skipped: every axis is vacuum-padded (0D/isolated system) -- a q-path "
                              "isn't physically meaningful here.", 'yellow'), f_out)
    else:
        from stb.phonons_pos import pretty_label
        print_dual(f"Bravais lattice   : {main_res['bravais_name']}", f_out)
        print_dual("Path              : " + " | ".join("-".join(pretty_label(lab) for lab in seg)
                                                         for seg in main_res['path_segments']), f_out)
    for lab, res, _c in series:
        report_stability_and_bands(lab, res, stats[lab], f_out, warn)
    if "bs" in main_res:
        p = os.path.join(args.output_dir, "bands.png")
        figures.append(plot_bands(series, p))
        files.append(("Band structure", p))
        if args.save_data:
            p = os.path.join(args.output_dir, "bands.dat")
            write_bands_dat(p, main_res)
            files.append(("Data", p))
            if found_res is not None:
                p = os.path.join(args.output_dir, "bands_foundation.dat")
                write_bands_dat(p, found_res)
                files.append(("Data", p))

    if args.freeze_unstable_mode:
        sections("MODE FREEZE")
        out = freeze_unstable_mode(ph, args, structure.species_meta, f_out)
        if out:
            files.append(("Frozen-mode cell", out))
    if args.animate_mode is not None:
        sections("MODE ANIMATION")
        if args.animate_mode == ANIMATE_OPTICAL:
            # bands are sorted by frequency: the 3 lowest are the acoustic branches
            args.animate_mode = 3
            print_dual("Mode              : first optical mode (band 3; bands 0-2 are acoustic)", f_out)
        frames, n_bands, freq = animate_normal_mode(ph, args.animate_qpoint, args.animate_mode,
                                                    args.animate_amplitude, args.animate_frames, args.animate_dim)
        if frames is None and n_bands == 3:
            warn("--animate-mode: the primitive cell has one atom, so there is no optical mode -- "
                 "only the 3 acoustic bands (0-2). Skipping animation.")
            print_dual(color_text("No optical mode (1-atom primitive cell): animation skipped", 'red'), f_out)
        elif frames is None:
            warn(f"--animate-mode {args.animate_mode} out of range -- only {n_bands} band(s) at q = "
                 f"{args.animate_qpoint}. Skipping animation.")
            print_dual(color_text(f"--animate-mode {args.animate_mode} out of range ({n_bands} bands): skipped", 'red'), f_out)
        else:
            ase_format, ext = OUTPUT_FORMATS[args.animate_format]
            anim_path = os.path.join(args.output_dir, f"mode_{args.animate_mode}{ext}")
            ase_write(anim_path, frames, format=ase_format)
            print_dual(f"q-point           : {args.animate_qpoint}", f_out)
            print_dual(f"Band index        : {args.animate_mode} ({freq:.4f} THz)", f_out)
            print_dual(f"Frames            : {len(frames)} ({args.animate_format} format)", f_out)
            print_dual(f"[OK] Saved {anim_path}", f_out)
            files.append(("Mode animation", anim_path))

    sections("DENSITY OF STATES")
    for lab, res, _c in series:
        report_dos(lab, res, stats[lab], f_out)
    p = os.path.join(args.output_dir, "dos.png")
    figures.append(plot_dos(series, p))
    files.append(("DOS", p))
    if args.save_data:
        p = os.path.join(args.output_dir, "dos.dat")
        write_dos_dat(p, main_res)
        files.append(("Data", p))

    sections("THERMAL PROPERTIES")
    print_dual("Harmonic, at the input volume. Classical limit: Cv -> 3 k_B per atom (Dulong-Petit).", f_out)
    report_thermal(series, f_out)
    p = os.path.join(args.output_dir, "thermal.png")
    figures.append(plot_thermal(series, p))
    files.append(("Thermal properties", p))
    if args.save_data:
        p = os.path.join(args.output_dir, "thermal.dat")
        write_thermal_dat(p, main_res)
        files.append(("Data", p))

    if found_res is not None:
        sections("FINE-TUNED vs FOUNDATION")
        report_comparison(main_res, found_res, f_out)

    if args.save_gnuplot:
        pd = os.path.join(args.output_dir, PLOT_DIR)
        writers = [lambda: write_gnuplot_dos(pd, series), lambda: write_gnuplot_thermal(pd, series)]
        if "bs" in main_res:
            writers.insert(0, lambda: write_gnuplot_bands(pd, series))
        gnuplot(writers)

    st = stats[label_main]
    summary = ["Status            : OK" + (f" ({len(warns)} warning(s), see the last section)" if warns else "")]
    if "bs" in main_res:
        g = freqs_at(ph, [0, 0, 0])
        summary.append(f"Results           : highest frequency {st['f_max']:.3f} THz, Gamma optical "
                       f"{', '.join(sorted({f'{x:.2f}' for x in g[3:]}, key=float)) or '-'} THz, "
                       f"Debye T {fmt(st['debye_K'], '{:.0f}')} K, imaginary q {st['imag_q']}/{st['n_q']}")
    finish(summary)


if __name__ == "__main__":
    main()
