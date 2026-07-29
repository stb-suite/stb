#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""Physical post-processing of an AIMD/MD trajectory: radial distribution
function g(r), mean-squared displacement (MSD) and diffusion coefficient, a
VACF-derived vibrational density of states (VDOS), per-atom displacement /
atom-pair relative-distance tracking, and a thermodynamic time series
(energy/temperature/volume/pressure). Input is either a SIESTA AIMD run
(<label>.ANI + .out) or a generic ASE-readable trajectory (e.g. one written
by stb-mlmd). Complements stb-ani2traj, which only converts a SIESTA
trajectory's format for external viewers.

Physics notes:
 - RDF (compute_rdf): minimum-image fractional-coordinate distances,
   re-derived per frame from that frame's own cell (so a variable-cell/
   NPT run stays correct), normalized against the standard ideal-gas
   shell-volume density.
 - MSD/diffusion (compute_msd/fit_diffusion_coefficient): sliding-time-
   origin averaging (not just t=0), 3D Einstein relation D = slope/6, and
   the Ang^2/fs -> cm^2/s conversion (1e-16/1e-15, explicit rather than a
   folded magic constant).
 - VACF/VDOS (compute_vacf_vdos): central-difference velocities (the only
   option -- .ANI carries no velocities), FFT-based autocorrelation
   (Wiener-Khinchin), frequency axis fs^-1 -> cm^-1 conversion via the
   physical constant c = 2.99792458e10 cm/s.

Output/report style matches the rest of the Analysis category (stb-
wfdensity/stb-sts/stb-coop/stb-ipr/stb-effmass/stb-spintexture): a
numbered [0]...[10] report, --save-report, --save-gnuplot (writes a .dat +
.gplot pair for every computed quantity instead of unconditional
matplotlib PNGs), --view (an on-demand interactive matplotlib preview),
and -o/--output-dir.

--track-atom N (0-based atom index) reports the Cartesian displacement of
one specific atom from its own initial position, frame by frame -- built
on the already-PBC-unwrapped trajectory (core.md_traj.unwrap_trajectory)
the MSD/VDOS analyses already use, so a real diffusing/hopping atom's
displacement stays physically continuous instead of resetting every time
it crosses a periodic boundary. Useful e.g. for watching one diffusing
adatom/dopant/defect during a run instead of only the species-averaged
MSD.

--track-pair I-J (two 0-based atom indices, e.g. '0-5') reports the
relative separation between two SPECIFIC atoms every frame, using the
SAME minimum-image convention as compute_rdf (fractional-coordinate
difference, wrapped to [-0.5, 0.5), converted through that frame's own
cell) rather than the unwrapped trajectory -- deliberately different from
--track-atom's convention, because what matters physically for a bond
length/hydrogen-bond distance/reaction coordinate is the atoms' TRUE
instantaneous separation. Unwrapping each atom independently would not
give this: two atoms can each accumulate their own independent drift/
unwrap path that does not cancel back down to their true short
separation, especially over a long trajectory. Verified on this tool's
own 2-atom O2 fixture: --track-pair 0-1 reproduces the same O-O
bond-length range as the RDF's own first peak (~1.0-1.3 Ang).

--geometry-file <path>: --label mode's MD timestep (MD.InitialTimeStep/
MD.LengthTimeStep) and lattice fallback both need the real SIESTA input
.fdf, whose filename is chosen by the user and can differ from the
SystemLabel inside it (e.g. SystemLabel 'siesta' with the real input
called calc.fdf) -- unlike .XV/.ANI/.HSX/.WFSX, which SIESTA itself
always names after SystemLabel. Without this override, a renamed/missing
.fdf silently degrades to "assume 1 fs per frame" with only a warning --
the same gap closed elsewhere in this suite via an explicit, --label-
decoupled file path (stb-sts/stb-coop/stb-ipr/stb-effmass/stb-
spintexture's own --geometry-file). core/md_traj.py's read_static_lattice/
read_frame_lattices/read_md_timestep_fs take an optional fdf_path
override instead of always reconstructing "<label>.fdf" internally;
ani2traj.py's own calls are unaffected (fdf_path defaults to None,
preserving the plain <label>.fdf guess there).

--list-atoms is a standalone early-exit mode: prints every atom's 0-based
index, species, and Cartesian coordinates (from the first frame only, so
it's fast regardless of trajectory length) then exits immediately,
without running RDF/MSD/VDOS/anything else. Lets --track-atom/
--track-pair's required indices (and which specific atoms they refer to,
spatially) be picked directly instead of guessed or hand-counted from the
.fdf. Deliberately a separate opt-in mode rather than a table always
printed in [1] INPUT DATA -- a real structure can easily have hundreds of
atoms, which would make every run's report unusably long by default. The
interactive stb-suite menu asks "List every atom's index/species/
coordinates?" (y/N) right before prompting for --track-atom/--track-pair
-- answering 'y' runs this and prints its output there; 'n' (or Enter)
skips it, same opt-in default as the CLI flag itself.

[7] THERMODYNAMIC TIME SERIES reports and plots energy/temperature/
volume/pressure vs time in one 4-panel figure (2x2, matplotlib via
--view; gnuplot via --save-gnuplot -- a `set multiplot layout 2,2`
script, write_multipanel_gplot, plus an individual .dat per quantity so
each can also be plotted standalone). For --label mode, Energy/
Temperature/Pressure are read from SIESTA's own dedicated '<label>.MDE'
file (core.siesta_log.get_mde_trajectory) -- a small, clean per-step
table (Step, T, E_KS, E_tot, Vol, P) SIESTA already writes for exactly
this purpose, rather than re-scraping scattered .out log lines the way
get_md_trajectory()'s cell/E_KS/Temp_ion does. Unlike the INPUT .fdf,
'.MDE' is always named after SystemLabel (like .XV/.ANI), so no
--geometry-file-style override is needed for it. Both E_tot (total
energy, kinetic+potential) and E_KS (the electronic/potential-like energy
alone) are plotted together on the Energy panel when both are available
-- E_tot trades between kinetic and potential energy by design, so only
the total is meaningful to check for conservation/drift (the same
distinction stb-mlmd's own energy tracking relies on). The Energy panel's
y-axis (matplotlib and gnuplot alike) is eV/atom, not the raw absolute
total -- a system-size-independent scale, comparable across different
structures/runs; the absolute eV values are still written as extra
columns in <stem>_energy.dat (and still printed in the [7] report table)
for completeness, just not what's plotted. Volume is always available
regardless of input source (computed directly from each frame's own
cell); for a generic --trajectory input, Energy/Temperature fall back to
each frame's own embedded 'Epot'/'Temp' info if present (e.g. from
stb-mlmd --out-format xyz) and Pressure is never available (no known
source). Verified against a real 500-step SIESTA Nose-thermostat (NVT)
AIMD run (8-atom SiC supercell, target 500 K): Volume exactly constant
(std = 0.0000 Ang^3, a fixed-cell NVT run), Temperature fluctuating
around the 500 K target (mean 503.1 K, as expected for a small 8-atom
system's canonical-ensemble fluctuations), and E_tot's std (0.0029 eV)
vs. E_KS's std (0.4158 eV) -- ~140x smaller -- confirming E_tot, not
E_KS, is the physically appropriate "conserved-ish" quantity to watch.
"""

VERSION = "2.0.0"

import os
import sys
import argparse
from collections import Counter
from datetime import datetime

import numpy as np

from stb.core.deps import require_sisl
sisl = require_sisl()

from stb.core import citations, siesta_log
from stb.core.md_traj import read_frame_lattices, read_md_timestep_fs, unwrap_trajectory
from stb.core.cli import color_text, show_intro, print_dual, print_section, print_table

REPORT_FILE = "stb_aimdAnalysis_report.txt"
BIB_FILE = "references.bib"


def cell_min_perpendicular_width(cell):
    """Half the shortest perpendicular distance between opposite faces of
    the cell (V / |cross product of the other two vectors|, minimized over
    the 3 axis pairs) -- the correct upper bound for a minimum-image RDF cutoff
    (unlike just taking half the shortest lattice VECTOR norm, this stays
    correct for a non-orthogonal/skewed cell too)."""
    a, b, c = cell
    vol = abs(np.dot(a, np.cross(b, c)))
    widths = [
        vol / np.linalg.norm(np.cross(b, c)),
        vol / np.linalg.norm(np.cross(c, a)),
        vol / np.linalg.norm(np.cross(a, b)),
    ]
    return min(widths) / 2.0


def compute_rdf(frac_positions, cells, symbols, pair=None, r_max=None, n_bins=200):
    """Average radial distribution function g(r) over all frames.

    `frac_positions` is a list (one per frame) of (natoms, 3) FRACTIONAL
    coordinate arrays -- distances are computed in fractional space (minimum-
    image convention, `disp -= round(disp)`) and converted to Cartesian via
    each frame's own cell, so this stays correct even for a variable-cell
    trajectory. Vectorized over all pairs at once per frame (numpy
    broadcasting) rather than building an ase.Atoms per frame, since a long
    AIMD trajectory can have thousands of frames.

    `pair`, if given, is a (species_a, species_b) tuple restricting the
    histogram to that pair only (species_a == species_b is a valid same-
    species RDF); None means every pair (the total RDF).

    Returns (r_centers, g_r). r_max defaults to the smallest
    cell_min_perpendicular_width found across all frames (the physical limit
    of the minimum-image convention).
    """
    natoms = frac_positions[0].shape[0]
    if r_max is None:
        r_max = min(cell_min_perpendicular_width(c) for c in cells)

    if pair is not None:
        idx_a = [i for i, s in enumerate(symbols) if s == pair[0]]
        idx_b = [i for i, s in enumerate(symbols) if s == pair[1]]
    else:
        idx_a = idx_b = list(range(natoms))

    bin_edges = np.linspace(0.0, r_max, n_bins + 1)
    counts = np.zeros(n_bins)
    n_frames = len(frac_positions)
    density_sum = 0.0

    for frac, cell in zip(frac_positions, cells):
        vol = abs(np.dot(cell[0], np.cross(cell[1], cell[2])))
        pa = frac[idx_a]
        pb = frac[idx_b]
        disp = pa[:, None, :] - pb[None, :, :]
        disp -= np.round(disp)
        cart = disp @ cell
        dist = np.linalg.norm(cart, axis=-1)
        same_species_block = pair is None or pair[0] == pair[1]
        if same_species_block and idx_a is idx_b:
            np.fill_diagonal(dist, np.inf)
        dist = dist[(dist > 1e-6) & (dist < r_max)]
        h, _ = np.histogram(dist, bins=bin_edges)
        counts += h
        density_sum += len(idx_b) / vol

    r_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    dr = bin_edges[1] - bin_edges[0]
    mean_density = density_sum / n_frames
    shell_vol = 4.0 * np.pi * r_centers**2 * dr
    norm = len(idx_a) * n_frames * shell_vol * mean_density
    norm[norm == 0] = np.nan
    g_r = counts / norm
    return r_centers, g_r


def compute_msd(cart_positions, dt_fs, symbols, species=None):
    """Mean-squared displacement vs time, averaged over multiple time
    origins (a sliding window over the whole trajectory, not just t=0 --
    standard MD practice for reducing statistical noise) and over every
    atom of `species` (or all atoms if None).

    `cart_positions` must already be PBC-unwrapped (core.md_traj.
    unwrap_trajectory) -- an MSD computed on wrapped/re-folded coordinates
    would be dominated by box-edge teleports instead of real diffusion.

    Returns (t_fs, msd) with t_fs starting at 0.
    """
    idx = ([i for i, s in enumerate(symbols) if s == species]
           if species is not None else list(range(len(symbols))))
    traj = np.array([p[idx] for p in cart_positions])  # (nframes, nsel, 3)
    n_frames = traj.shape[0]

    msd = np.zeros(n_frames)
    for lag in range(n_frames):
        disp = traj[lag:] - traj[:n_frames - lag]
        msd[lag] = np.mean(np.sum(disp**2, axis=-1))

    t_fs = np.arange(n_frames) * dt_fs
    return t_fs, msd


def fit_diffusion_coefficient(t_fs, msd, fit_start=None, fit_end=None):
    """Linear fit of MSD(t) over [fit_start, fit_end] fs (default: the second
    half of the trajectory, where the regime is diffusive rather than the
    short-time ballistic ramp) via the 3D Einstein relation D = slope / 6.
    Returns (D_cm2_s, slope, (t_lo, t_hi)). slope is in Ang^2/fs; the
    Ang^2/fs -> cm^2/s conversion factor is 1e-16 (Ang^2 -> cm^2 = 1e-16) /
    1e-15 (fs -> s) = 1e-1... done explicitly below rather than folded into
    a magic constant.
    """
    t_hi = t_fs[-1]
    t_lo = fit_start if fit_start is not None else t_hi / 2.0
    t_hi = fit_end if fit_end is not None else t_hi
    mask = (t_fs >= t_lo) & (t_fs <= t_hi)
    if mask.sum() < 2:
        return None, None, (t_lo, t_hi)
    slope, _ = np.polyfit(t_fs[mask], msd[mask], 1)  # Ang^2 / fs
    ang2_to_cm2 = 1.0e-16
    fs_to_s = 1.0e-15
    D_cm2_s = (slope / 6.0) * ang2_to_cm2 / fs_to_s
    return D_cm2_s, slope, (t_lo, t_hi)


def compute_vacf_vdos(cart_positions, dt_fs):
    """Velocity autocorrelation function (VACF) and the vibrational density
    of states (VDOS) derived from it.

    .ANI carries no velocities, only positions -- v(t) is estimated by
    central finite difference, v(t) = (x(t+dt) - x(t-dt)) / (2*dt), the
    standard approximation when only a position trajectory is available.
    This is markedly less accurate than SIESTA's own MD velocities would be
    (it amplifies any positional noise into the derivative) -- treat the
    resulting VDOS as qualitative, not benchmark-grade.

    VACF is normalized (VACF(0) = 1), averaged over atoms and over sliding
    time origins (via FFT-based autocorrelation, O(N log N) instead of the
    O(N^2) direct sum -- matters for long trajectories). VDOS is the real
    part of the FFT of the VACF (Wiener-Khinchin theorem), with the
    frequency axis converted from 1/fs to cm^-1 (1 fs^-1 = 1e15 Hz;
    freq_cm-1 = freq_Hz / (c * 100), c in m/s).

    Returns (freq_cm1, vdos, t_fs, vacf) -- freq_cm1/vdos only cover the
    non-negative-frequency half of the spectrum.
    """
    pos = np.asarray(cart_positions)  # (nframes, natoms, 3)
    vel = np.gradient(pos, dt_fs, axis=0)
    n_frames, natoms, _ = vel.shape

    vacf = np.zeros(n_frames)
    for a in range(natoms):
        for d in range(3):
            v = vel[:, a, d]
            v = v - v.mean()
            n = len(v)
            padded = np.zeros(2 * n)
            padded[:n] = v
            f = np.fft.fft(padded)
            acf = np.fft.ifft(f * np.conj(f)).real[:n]
            acf /= np.arange(n, 0, -1)
            vacf += acf
    vacf /= (natoms * 3)
    if vacf[0] != 0:
        vacf /= vacf[0]

    t_fs = np.arange(n_frames) * dt_fs

    spectrum = np.fft.rfft(vacf)
    freqs_per_fs = np.fft.rfftfreq(n_frames, d=dt_fs)
    c_cm_per_s = 2.99792458e10
    freq_cm1 = (freqs_per_fs * 1.0e15) / c_cm_per_s
    vdos = np.abs(spectrum.real)
    return freq_cm1, vdos, t_fs, vacf


def compute_atom_displacement(cart_positions, dt_fs, atom_index):
    """Cartesian displacement of ONE atom from its own initial position,
    frame by frame. `cart_positions` must already be PBC-unwrapped (same
    input compute_msd uses), so a real diffusing/hopping atom's
    displacement stays physically continuous instead of resetting every
    time it crosses a periodic boundary.

    Returns (t_fs, traj, disp, disp_mag): `traj` is the atom's own
    unwrapped (natoms-independent) Cartesian trajectory (nframes, 3), for
    reporting its initial/final absolute position; `disp` is traj - traj[0]
    (nframes, 3); `disp_mag` is its norm (nframes,).
    """
    traj = np.array([p[atom_index] for p in cart_positions])
    disp = traj - traj[0]
    disp_mag = np.linalg.norm(disp, axis=1)
    t_fs = np.arange(len(cart_positions)) * dt_fs
    return t_fs, traj, disp, disp_mag


def compute_relative_distance(frac_positions, cells, dt_fs, atom_i, atom_j):
    """Relative separation vector/distance between two SPECIFIC atoms,
    frame by frame, using the SAME minimum-image convention as compute_rdf
    (fractional-coordinate difference, wrapped to [-0.5, 0.5), converted
    through that frame's own cell -- correct even for a variable-cell run).

    Deliberately NOT based on the unwrapped trajectory compute_atom_
    displacement uses: what matters physically for a bond length/hydrogen
    -bond distance/reaction coordinate is the atoms' TRUE instantaneous
    separation, which minimum-image gives directly. Unwrapping the two
    atoms independently would not give this in general -- each atom can
    accumulate its own separate drift/unwrap path that does not cancel back
    down to the true short separation.

    Returns (t_fs, vec, dist): `vec` is the j-i displacement vector every
    frame (nframes, 3, Angstrom); `dist` is its norm (nframes,).
    """
    n_frames = len(frac_positions)
    vec = np.zeros((n_frames, 3))
    dist = np.zeros(n_frames)
    for k in range(n_frames):
        disp = frac_positions[k][atom_j] - frac_positions[k][atom_i]
        disp -= np.round(disp)
        cart = disp @ cells[k]
        vec[k] = cart
        dist[k] = np.linalg.norm(cart)
    t_fs = np.arange(n_frames) * dt_fs
    return t_fs, vec, dist


def write_xy_data(path, header, columns):
    """Plain np.savetxt wrapper -- writes `columns` (a list of 1D arrays,
    all the same length, first one the shared x-axis) side by side as a
    whitespace-separated .dat file with a one-line header."""
    np.savetxt(path, np.column_stack(columns), header=header)


def write_xy_gplot(gplot_path, dat_path, pdf_path, xlabel, ylabel, title, series):
    """Writes a .gplot script plotting one or more y-series of a plain
    multi-column x-y .dat file. `series` is a list of (column_index,
    curve_label) tuples (1-based; column 1 is always the shared x-axis).

    The plain-XY analog of core/xrd.py's/core/band_scatter.py's own per
    -tool .gplot writers -- kept local here rather than promoted to
    core/ since aimd_analysis.py is still the only consumer of a generic
    multi-column time-series plot in this suite (same extract-on-second
    -use policy as the rest of core/).
    """
    dat_basename = os.path.basename(dat_path)
    lines = [
        'set terminal pdfcairo enhanced font "Arial,14" size 8,6\n',
        f'set output "{os.path.basename(pdf_path)}"\n',
        f'set xlabel "{xlabel}" font "Arial,16"\n',
        f'set ylabel "{ylabel}" font "Arial,16"\n',
        f'set title "{title}"\n',
        'set grid xtics ytics lt 0 lw 1 lc rgb "#bbbbbb"\n',
        ('set key outside\n' if len(series) > 1 else 'unset key\n'),
        '\n',
    ]
    plot_terms = ", ".join(
        f'"{dat_basename}" using 1:{col} with lines lw 2 title "{label}"'
        for col, label in series
    )
    lines.append(f'plot {plot_terms}\n')
    with open(gplot_path, "w") as f:
        f.writelines(lines)


def write_multipanel_gplot(gplot_path, pdf_path, panels):
    """Writes a single .gplot script laying out up to 4 panels (`set
    multiplot layout 2,2`), one per (dat_path, xlabel, ylabel, title,
    column_index) tuple in `panels` -- used for the combined energy/
    temperature/volume/pressure figure, since each quantity has its own
    units/scale and its own separate .dat file, so it's more readable on
    its own small axis than overlaid on one shared plot.
    """
    lines = [
        'set terminal pdfcairo enhanced font "Arial,11" size 10,8\n',
        f'set output "{os.path.basename(pdf_path)}"\n',
        'set multiplot layout 2,2\n',
        'set grid xtics ytics lt 0 lw 1 lc rgb "#bbbbbb"\n',
        'unset key\n',
        '\n',
    ]
    for dat_path, xlabel, ylabel, title, col in panels:
        dat_basename = os.path.basename(dat_path)
        lines.append(f'set xlabel "{xlabel}" font "Arial,13"\n')
        lines.append(f'set ylabel "{ylabel}" font "Arial,13"\n')
        lines.append(f'set title "{title}"\n')
        lines.append(f'plot "{dat_basename}" using 1:{col} with lines lw 2 lc rgb "#1f77b4"\n')
        lines.append('\n')
    lines.append('unset multiplot\n')
    with open(gplot_path, "w") as f:
        f.writelines(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Physical post-processing of an AIMD/MD trajectory -- a SIESTA AIMD run "
                     "(<label>.ANI) or a generic ASE-readable trajectory (e.g. one written by "
                     "stb-mlmd): radial distribution function g(r), mean-squared displacement "
                     "(MSD) and diffusion coefficient, a VACF-derived vibrational density of "
                     "state (VDOS), single-atom displacement tracking, and atom-pair relative "
                     "distance tracking. Complements stb-ani2traj, which only converts a SIESTA "
                     "trajectory's format for external viewers.",
        epilog="Example usage:\n"
               "  stb-aimdAnalysis --label aimd\n"
               "  stb-aimdAnalysis -l aimd --pair O-H --skip 50\n"
               "  stb-aimdAnalysis -l aimd --fit-start 5000 --fit-end 15000\n"
               "  stb-aimdAnalysis -l aimd --track-atom 0 --track-pair 0-1\n"
               "  stb-aimdAnalysis -l aimd --save-report --save-gnuplot --view\n"
               "  stb-aimdAnalysis --trajectory si_bulk_md_traj.xyz --dt 10",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("-l", "--label", default=None,
                        help="SIESTA SystemLabel -- reads '<label>.ANI'. Mutually exclusive with --trajectory.")
    parser.add_argument("--geometry-file", type=str, default=None, metavar="PATH",
                        help="Explicit path to the real SIESTA .fdf input file (--label mode "
                             "only). MD.InitialTimeStep/MD.LengthTimeStep (the real MD timestep) "
                             "and, if '<label>.XV'/'<label>.out' aren't usable, the lattice "
                             "itself are read from here. Auto-detecting '<label>.fdf' almost "
                             "never works in practice -- the input file's own name is chosen by "
                             "the user and is very often NOT <label>.fdf (e.g. SystemLabel "
                             "'siesta' with the real input called calc.fdf) -- pass this "
                             "explicitly instead of relying on that guess.")
    parser.add_argument("--trajectory", default=None, metavar="PATH",
                        help="Read a generic ASE-readable multi-frame trajectory instead of a "
                             "SIESTA .ANI -- e.g. one written by stb-mlmd (xsf/pdb/xyz). "
                             "Mutually exclusive with --label.")
    parser.add_argument("--dt", type=float, default=None, metavar="FS",
                        help="Real time between saved frames in the --trajectory file, in fs "
                             "(--trajectory mode only). Auto-detected if the file is an "
                             "extended-xyz with per-frame 'Time' info (e.g. stb-mlmd --out-format "
                             "xyz); required otherwise (xsf/pdb carry no such data).")
    parser.add_argument("--stride", type=int, default=1,
                        help="Keep only every Nth frame (default: 1, every frame).")
    parser.add_argument("--skip", type=int, default=0,
                        help="Discard the first N frames (post-stride) as thermal "
                             "equilibration, before any statistics are computed (default: 0).")
    parser.add_argument("--pair", default=None, metavar="A-B",
                        help="Restrict the RDF to one species pair, e.g. 'O-H' (default: "
                             "total RDF over all pairs).")
    parser.add_argument("--r-max", type=float, default=None,
                        help="RDF cutoff in Angstrom (default: the minimum-image limit of "
                             "the smallest cell found in the trajectory).")
    parser.add_argument("--n-bins", type=int, default=200, help="Number of RDF bins (default: 200).")
    parser.add_argument("--fit-start", type=float, default=None, metavar="FS",
                        help="MSD linear-fit window start, in fs (default: half the "
                             "trajectory length).")
    parser.add_argument("--fit-end", type=float, default=None, metavar="FS",
                        help="MSD linear-fit window end, in fs (default: end of trajectory).")
    parser.add_argument("--track-atom", type=int, default=None, metavar="N",
                        help="0-based atom index whose own Cartesian displacement (from its "
                             "initial position, PBC-unwrapped -- same convention as the MSD) is "
                             "tracked and reported over time (optional).")
    parser.add_argument("--track-pair", default=None, metavar="I-J",
                        help="Two 0-based atom indices, e.g. '0-5', whose relative distance "
                             "(minimum-image convention, same as the RDF) is tracked and "
                             "reported over time (optional).")
    parser.add_argument("--list-atoms", action="store_true",
                        help="Print a table of every atom's 0-based index, species, and "
                             "Cartesian coordinates (Ang, from the first frame only) to help "
                             "pick --track-atom/--track-pair indices, then exit immediately -- "
                             "no RDF/MSD/VDOS/report is computed. Off by default (a large "
                             "structure can have hundreds of atoms).")

    parser.add_argument("-o", "--output-dir", type=str, default=".",
                        help="Directory to write <stem>_*.dat/.gplot (with --save-gnuplot) and "
                             f"{REPORT_FILE}/{BIB_FILE} into (default: current directory). "
                             "Created if it doesn't exist.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--save-gnuplot", action="store_true",
                        help="Also write a .dat + .gplot pair for every computed quantity "
                             "(RDF, MSD, VACF, VDOS, and the atom-displacement/pair-distance "
                             "series if --track-atom/--track-pair are used). Off by default -- "
                             "this tool previously wrote matplotlib PNGs unconditionally on "
                             "every run with no gnuplot output at all.")
    parser.add_argument("--view", action="store_true",
                        help="Show an interactive matplotlib preview before finishing. Off by "
                             "default -- this tool previously always generated (and saved) the "
                             "matplotlib plots, with no way to skip it.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-aimdAnalysis {VERSION}")

    args = parser.parse_args()

    if not args.label and not args.trajectory:
        parser.error("one of --label or --trajectory is required.")
    if args.label and args.trajectory:
        parser.error("--label and --trajectory are mutually exclusive.")
    if args.geometry_file and not args.label:
        parser.error("--geometry-file only applies to --label mode (a --trajectory input "
                     "carries its own per-frame cell, no .fdf involved).")
    if args.geometry_file and not os.path.isfile(args.geometry_file):
        parser.error(f"--geometry-file '{args.geometry_file}' not found.")

    pair = None
    if args.pair:
        parts = args.pair.split("-")
        if len(parts) != 2:
            parser.error(f"--pair must be 'A-B' (e.g. 'O-H'), got '{args.pair}'.")
        pair = tuple(parts)

    track_pair = None
    if args.track_pair:
        parts = args.track_pair.split("-")
        if len(parts) != 2:
            parser.error(f"--track-pair must be 'I-J' (e.g. '0-5'), got '{args.track_pair}'.")
        try:
            track_pair = (int(parts[0]), int(parts[1]))
        except ValueError:
            parser.error(f"--track-pair indices must be integers, got '{args.track_pair}'.")
        if track_pair[0] == track_pair[1]:
            parser.error("--track-pair requires two DIFFERENT atom indices.")

    # Output filenames are keyed off `stem` regardless of input mode --
    # args.label only resolves the SIESTA <label>.ANI/.out/.XV/.fdf input set.
    if args.label:
        ani_file = f"{args.label}.ANI"
        stem = args.label
    else:
        ani_file = args.trajectory
        stem = os.path.splitext(os.path.basename(args.trajectory))[0]

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - AIMD Trajectory Analysis",
            "RDF, MSD/diffusion, VACF-derived VDOS, atom/pair tracking",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    if args.list_atoms:
        if not os.path.exists(ani_file):
            parser.error(f"Input file '{ani_file}' not found.")
        if args.label:
            try:
                first_frame = sisl.get_sile(ani_file).read_geometry[0]()
            except Exception as e:
                parser.error(f"Could not read '{ani_file}': {e}")
            list_symbols = [a.symbol for a in first_frame.atoms]
            list_positions = first_frame.xyz
        else:
            from ase.io import read as ase_read
            try:
                first_frame = ase_read(args.trajectory, index=0)
            except Exception as e:
                parser.error(f"Could not read '{args.trajectory}': {e}")
            list_symbols = first_frame.get_chemical_symbols()
            list_positions = first_frame.get_positions()

        print(color_text(f"\n{stem}: {len(list_symbols)} atom(s), coordinates from the first "
                         "frame (Ang):", 'bold'))
        print_table(["Index", "Species", "X", "Y", "Z"], [
            ([str(i), sp, f"{pos[0]:.4f}", f"{pos[1]:.4f}", f"{pos[2]:.4f}"], None)
            for i, (sp, pos) in enumerate(zip(list_symbols, list_positions))
        ])
        sys.exit(0)

    print("\n" + color_text("AIMD TRAJECTORY ANALYSIS:", 'bold'))
    print("-" * 60)

    os.makedirs(args.output_dir, exist_ok=True)
    report_path = os.path.join(args.output_dir, REPORT_FILE) if args.save_report else None
    f_out = open(report_path, "w") if report_path else None

    def fail(message):
        print_dual(color_text(f"[FAIL] {message}", 'red'), f_out)
        if f_out:
            f_out.close()
        sys.exit(1)

    print_dual(color_text("===== STB-AIMDANALYSIS REPORT =====", 'magenta'), f_out)

    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Input file        : {ani_file}", f_out)
    if args.geometry_file:
        print_dual(f"Geometry file     : {args.geometry_file}", f_out)
    print_dual(f"Stride / skip     : {args.stride} / {args.skip}", f_out)
    if pair:
        print_dual(f"RDF pair          : {pair[0]}-{pair[1]}", f_out)
    if args.track_atom is not None:
        print_dual(f"Track atom        : index {args.track_atom}", f_out)
    if track_pair:
        print_dual(f"Track pair        : atoms {track_pair[0]}-{track_pair[1]}", f_out)
    print_dual(f"Output dir        : {args.output_dir}", f_out)
    print_dual(f"Save gnuplot      : {'yes' if args.save_gnuplot else 'no'}", f_out)
    print_dual(f"View (matplotlib) : {'yes' if args.view else 'no'}", f_out)
    if report_path:
        print_dual(f"Report file       : {report_path}", f_out)

    print_section("[1] INPUT DATA", f_out)

    if not os.path.exists(ani_file):
        fail(f"Input file '{ani_file}' not found.")
    if args.stride < 1:
        fail(f"--stride must be >= 1, got {args.stride}.")

    if args.label:
        try:
            sile = sisl.get_sile(ani_file)
            all_frames = sile.read_geometry[:]()
        except Exception as e:
            fail(f"Could not read '{ani_file}': {e}")

        if not all_frames:
            fail(f"'{ani_file}' contains no frames.")

        all_cells, all_steps = read_frame_lattices(
            args.label, len(all_frames), f_out, fdf_path=args.geometry_file)
        if all_cells is None:
            fdf_desc = args.geometry_file if args.geometry_file else f"{args.label}.fdf (auto-detected)"
            fail(f"Could not find a readable '{args.label}.out', '{args.label}.XV' or "
                 f"'{fdf_desc}' -- a lattice is required ('{ani_file}' doesn't carry "
                 f"one of its own). If the real input file isn't named '{args.label}.fdf', "
                 f"pass --geometry-file <path>.")

        all_mde = siesta_log.get_mde_trajectory(args.label)
        mde_usable = bool(all_mde) and len(all_mde) >= len(all_frames)

        frames = all_frames[::args.stride]
        cells = all_cells[::args.stride]
        thermo_steps = all_mde[::args.stride] if mde_usable else None

        if args.skip >= len(frames):
            fail(f"--skip {args.skip} discards all {len(frames)} available (post-stride) frame(s).")
        frames = frames[args.skip:]
        cells = cells[args.skip:]
        if thermo_steps is not None:
            thermo_steps = thermo_steps[args.skip:]

        symbols = [a.symbol for a in frames[0].atoms]
        if pair and (pair[0] not in symbols or pair[1] not in symbols):
            fail(f"--pair {args.pair}: species not found in structure "
                 f"(available: {', '.join(sorted(set(symbols)))}).")

        print_dual(f"[OK] Read {len(frames)} frame(s) from {ani_file} "
                   f"(stride {args.stride}, skip {args.skip}).", f_out)

        if all_mde and not mde_usable:
            print_dual(color_text(
                f"[WARNING] '{args.label}.MDE' has {len(all_mde)} row(s) but {len(all_frames)} "
                f"frame(s) were read from '{ani_file}' -- skipping the thermodynamic time "
                "series (row/frame counts must match).", 'yellow'), f_out)

        initial_step, dt_fs = read_md_timestep_fs(args.label, fdf_path=args.geometry_file)
        if dt_fs is None:
            fdf_desc = args.geometry_file if args.geometry_file else f"{args.label}.fdf (auto-detected)"
            print_dual(color_text(
                f"[WARNING] MD.LengthTimeStep not found in '{fdf_desc}' -- assuming 1 fs per "
                "(strided) frame for MSD/VDOS time axes. If the real input file isn't named "
                f"'{args.label}.fdf', pass --geometry-file <path>.",
                'yellow'), f_out)
            dt_fs = 1.0
        else:
            dt_fs *= args.stride

        frac_positions = [g.fxyz for g in frames]
        cart_positions_raw = [g.xyz for g in frames]

    else:
        # Generic ASE-readable trajectory (e.g. stb-mlmd's xsf/pdb/xyz output)
        # -- no sisl/SIESTA involved at all. Every ASE multi-frame reader
        # already carries a per-frame cell (that's the whole point of these
        # 3 formats, same as stb-ani2traj's own choice), so there's no
        # separate lattice-source lookup needed here.
        from ase.io import read as ase_read
        try:
            all_frames = ase_read(args.trajectory, index=':')
        except Exception as e:
            fail(f"Could not read '{args.trajectory}': {e}")

        if not all_frames:
            fail(f"'{args.trajectory}' contains no frames.")

        all_cells = [np.array(f.cell) for f in all_frames]
        frames = all_frames[::args.stride]
        cells = all_cells[::args.stride]

        if args.skip >= len(frames):
            fail(f"--skip {args.skip} discards all {len(frames)} available (post-stride) frame(s).")
        frames = frames[args.skip:]
        cells = cells[args.skip:]

        symbols = frames[0].get_chemical_symbols()
        if pair and (pair[0] not in symbols or pair[1] not in symbols):
            fail(f"--pair {args.pair}: species not found in structure "
                 f"(available: {', '.join(sorted(set(symbols)))}).")

        print_dual(f"[OK] Read {len(frames)} frame(s) from {args.trajectory} "
                   f"(stride {args.stride}, skip {args.skip}).", f_out)

        if args.dt is not None:
            dt_fs = args.dt * args.stride
        elif len(frames) >= 2 and 'Time' in frames[0].info and 'Time' in frames[1].info:
            dt_fs = frames[1].info['Time'] - frames[0].info['Time']
            print_dual(f"[INFO] Auto-detected dt = {dt_fs:.4f} fs from embedded per-frame "
                       "'Time' info.", f_out)
        else:
            fail("--dt is required for --trajectory input without embedded per-frame 'Time' "
                 "info (only an extended-xyz written with Time, e.g. stb-mlmd --out-format "
                 "xyz, carries this automatically).")

        # No .MDE (that's SIESTA-specific) -- Energy/Temperature come from each
        # frame's own embedded info instead, if the writer included them (e.g.
        # stb-mlmd --out-format xyz embeds 'Epot'/'Temp', same convention as
        # 'Time' above). No known source for Pressure in this generic path.
        thermo_steps = [
            {'E_tot': None, 'E_KS': f.info.get('Epot'), 'T': f.info.get('Temp'), 'pressure': None}
            for f in frames
        ] if any('Epot' in f.info or 'Temp' in f.info for f in frames) else None

        frac_positions = [f.get_scaled_positions(wrap=False) for f in frames]
        cart_positions_raw = [f.get_positions() for f in frames]

    natoms = len(symbols)
    if args.track_atom is not None and not (0 <= args.track_atom < natoms):
        fail(f"--track-atom {args.track_atom} out of range (structure has {natoms} atoms, "
             f"valid indices are 0-{natoms - 1}).")
    if track_pair and (not (0 <= track_pair[0] < natoms) or not (0 <= track_pair[1] < natoms)):
        fail(f"--track-pair {args.track_pair} out of range (structure has {natoms} atoms, "
             f"valid indices are 0-{natoms - 1}).")

    composition = ", ".join(f"{el}{n}" for el, n in sorted(Counter(symbols).items()))
    print_table(["Quantity", "Value"], [
        (["Frames analyzed", f"{len(frames)}"], None),
        (["Atoms per frame", f"{natoms}"], None),
        (["Composition", composition], None),
        (["Timestep (post-stride)", f"{dt_fs:.4f} fs"], None),
    ], f_out)

    cart_positions = unwrap_trajectory(cart_positions_raw, cells)

    print_section("[2] RADIAL DISTRIBUTION FUNCTION (RDF)", f_out)
    r, g_r = compute_rdf(frac_positions, cells, symbols, pair=pair,
                         r_max=args.r_max, n_bins=args.n_bins)
    valid = ~np.isnan(g_r)
    if valid.any():
        rdf_peak_idx = np.nanargmax(g_r)
        print_dual(f"[OK] RDF computed over {len(frames)} frame(s), r up to {r[-1]:.3f} Ang.", f_out)
        print_dual(f"First peak        : r = {r[rdf_peak_idx]:.3f} Ang (g(r) = {g_r[rdf_peak_idx]:.3f})", f_out)
    else:
        print_dual(color_text("[WARNING] RDF could not be computed (no pairs in range).", 'yellow'), f_out)

    print_section("[3] MEAN-SQUARED DISPLACEMENT (MSD) & DIFFUSION", f_out)
    species_set = sorted(set(symbols))
    msd_results = {}
    msd_rows = []
    for sp in species_set + [None]:
        t_fs, msd = compute_msd(cart_positions, dt_fs, symbols, species=sp)
        D, slope, (t_lo, t_hi) = fit_diffusion_coefficient(t_fs, msd, args.fit_start, args.fit_end)
        label = sp if sp is not None else "all"
        msd_results[label] = (t_fs, msd, D)
        if D is not None:
            msd_rows.append(([label, f"{D:.3e}", f"{t_lo:.1f}-{t_hi:.1f} fs"], None))
        else:
            msd_rows.append(([label, "n/a", f"{t_lo:.1f}-{t_hi:.1f} fs"], 'yellow'))
    print_table(["Species", "D (cm^2/s)", "Fit window"], msd_rows, f_out)
    if any(D is None for _, _, D in msd_results.values()):
        print_dual(color_text(
            "[WARNING] Not enough frames in the fit window for at least one series -- pass "
            "--fit-start/--fit-end or a longer trajectory.", 'yellow'), f_out)

    print_section("[4] VELOCITY AUTOCORRELATION / VIBRATIONAL DOS (VDOS)", f_out)
    vacf_result = None
    if len(frames) < 4:
        print_dual(color_text(
            "[WARNING] Trajectory too short for a meaningful VACF/VDOS "
            f"(only {len(frames)} frame(s)) -- skipping.", 'yellow'), f_out)
    else:
        freq_cm1, vdos, t_fs_vacf, vacf = compute_vacf_vdos(cart_positions, dt_fs)
        vacf_result = (freq_cm1, vdos, t_fs_vacf, vacf)
        vdos_peak_idx = np.argmax(vdos[1:]) + 1 if len(vdos) > 1 else 0
        print_dual(f"[OK] VDOS computed (finite-difference velocities -- qualitative, "
                   f"not benchmark-grade).", f_out)
        print_dual(f"Dominant peak     : {freq_cm1[vdos_peak_idx]:.1f} cm^-1", f_out)

    print_section("[5] SINGLE-ATOM DISPLACEMENT TRACKING", f_out)
    atom_disp_result = None
    if args.track_atom is not None:
        t_disp, traj_atom, disp, disp_mag = compute_atom_displacement(
            cart_positions, dt_fs, args.track_atom)
        species_tracked = symbols[args.track_atom]
        i_max = int(np.argmax(disp_mag))
        print_dual(f"Tracked atom      : index {args.track_atom} (species {species_tracked})", f_out)
        print_table(["Quantity", "Value"], [
            (["Initial position (Ang)",
              f"({traj_atom[0][0]:.4f}, {traj_atom[0][1]:.4f}, {traj_atom[0][2]:.4f})"], None),
            (["Final position (Ang)",
              f"({traj_atom[-1][0]:.4f}, {traj_atom[-1][1]:.4f}, {traj_atom[-1][2]:.4f})"], None),
            (["Net displacement vector (Ang)",
              f"({disp[-1][0]:+.4f}, {disp[-1][1]:+.4f}, {disp[-1][2]:+.4f})"], None),
            (["Net displacement magnitude (Ang)", f"{disp_mag[-1]:.4f}"], None),
            (["Max displacement magnitude (Ang)",
              f"{disp_mag[i_max]:.4f} at t = {t_disp[i_max]:.1f} fs"], None),
            (["Mean displacement magnitude (Ang)", f"{disp_mag.mean():.4f}"], None),
        ], f_out)
        atom_disp_result = (t_disp, traj_atom, disp, disp_mag, species_tracked)
    else:
        print_dual("Not requested (pass --track-atom N, 0-based, to track one atom's own "
                   "Cartesian displacement over time -- use --list-atoms to see valid indices "
                   "and their coordinates).", f_out)

    print_section("[6] ATOM-PAIR RELATIVE DISTANCE", f_out)
    pair_dist_result = None
    if track_pair:
        i_idx, j_idx = track_pair
        t_pair, vec_pair, dist_pair = compute_relative_distance(
            frac_positions, cells, dt_fs, i_idx, j_idx)
        print_dual(f"Tracked pair      : atom {i_idx} ({symbols[i_idx]}) -- "
                   f"atom {j_idx} ({symbols[j_idx]})", f_out)
        delta = dist_pair[-1] - dist_pair[0]
        print_table(["Quantity", "Value"], [
            (["Initial distance (Ang)", f"{dist_pair[0]:.4f}"], None),
            (["Final distance (Ang)", f"{dist_pair[-1]:.4f}"], None),
            (["Mean distance (Ang)", f"{dist_pair.mean():.4f}"], None),
            (["Std dev (Ang)", f"{dist_pair.std():.4f}"], None),
            (["Min distance (Ang)", f"{dist_pair.min():.4f}"], None),
            (["Max distance (Ang)", f"{dist_pair.max():.4f}"], None),
            (["Net change, final-initial (Ang)", f"{delta:+.4f}"], None),
        ], f_out)
        rel_change = abs(delta) / dist_pair[0] if dist_pair[0] > 1e-9 else 0.0
        if rel_change > 0.3:
            print_dual(color_text(
                f"[WARNING] Distance changed by {rel_change * 100:.1f}% over the trajectory -- "
                "check for a possible bond-breaking/forming or reactive event.", 'yellow'), f_out)
        pair_dist_result = (t_pair, vec_pair, dist_pair, i_idx, j_idx)
    else:
        print_dual("Not requested (pass --track-pair I-J, 0-based indices e.g. '0-5', to track "
                   "the minimum-image distance between two specific atoms over time -- use "
                   "--list-atoms to see valid indices and their coordinates).", f_out)

    print_section("[7] THERMODYNAMIC TIME SERIES (ENERGY / TEMPERATURE / VOLUME / PRESSURE)", f_out)
    t_thermo = np.arange(len(frames)) * dt_fs
    volumes = np.array([abs(np.dot(c[0], np.cross(c[1], c[2]))) for c in cells])

    def _extract_thermo(key):
        if thermo_steps is None:
            return None
        vals = [s.get(key) for s in thermo_steps]
        return np.array(vals, dtype=float) if all(v is not None for v in vals) else None

    e_tot = _extract_thermo('E_tot')
    e_pot = _extract_thermo('E_KS')
    temperature = _extract_thermo('T')
    pressure = _extract_thermo('pressure')

    thermo_rows = [(["Volume (Ang^3)", f"{volumes.mean():.4f}", f"{volumes.std():.4f}",
                     f"{volumes.min():.4f}", f"{volumes.max():.4f}"], None)]
    for row_label, arr in [("Energy, total (eV)", e_tot), ("Energy, potential (eV)", e_pot),
                           ("Temperature (K)", temperature), ("Pressure (kBar)", pressure)]:
        if arr is not None:
            thermo_rows.append(([row_label, f"{arr.mean():.4f}", f"{arr.std():.4f}",
                                 f"{arr.min():.4f}", f"{arr.max():.4f}"], None))
        else:
            thermo_rows.append(([row_label, "n/a", "n/a", "n/a", "n/a"], 'yellow'))
    for row_label, arr in [("Energy, total (eV/atom)", e_tot), ("Energy, potential (eV/atom)", e_pot)]:
        if arr is not None:
            per_atom = arr / natoms
            thermo_rows.append(([row_label, f"{per_atom.mean():.6f}", f"{per_atom.std():.6f}",
                                 f"{per_atom.min():.6f}", f"{per_atom.max():.6f}"], None))
        else:
            thermo_rows.append(([row_label, "n/a", "n/a", "n/a", "n/a"], 'yellow'))
    print_table(["Quantity", "Mean", "Std dev", "Min", "Max"], thermo_rows, f_out)

    if thermo_steps is None:
        if args.label:
            print_dual(f"Energy/Temperature/Pressure not available ('{args.label}.MDE' not "
                       "found, or its row count doesn't match the number of frames read). "
                       "Volume is always available (computed directly from each frame's own "
                       "cell).", f_out)
        else:
            print_dual("Energy/Temperature not available (the --trajectory input has no "
                       "embedded 'Epot'/'Temp' per-frame info, e.g. from stb-mlmd --out-format "
                       "xyz). Pressure is never available for a generic --trajectory input. "
                       "Volume is always available (computed directly from each frame's own "
                       "cell).", f_out)
    elif pressure is None and args.label:
        print_dual("Pressure not available in this '.MDE' file (older SIESTA versions can omit "
                   "it for a fixed-shape-cell run with no barostat).", f_out)

    print_section("[8] WRITING OUTPUT FILES (GNUPLOT)", f_out)
    written_files = []
    if args.save_gnuplot:
        pair_tag = f"_{pair[0]}{pair[1]}" if pair else ""
        rdf_dat = os.path.join(args.output_dir, f"{stem}_rdf{pair_tag}.dat")
        rdf_gplot = os.path.join(args.output_dir, f"{stem}_rdf{pair_tag}.gplot")
        write_xy_data(rdf_dat, "r(Ang) g(r)", [r, g_r])
        write_xy_gplot(rdf_gplot, rdf_dat, f"{stem}_rdf{pair_tag}.pdf", "r (Ang)", "g(r)",
                       f"RDF{' (' + pair[0] + '-' + pair[1] + ')' if pair else ''} -- {stem}",
                       [(2, "g(r)")])
        written_files += [rdf_dat, rdf_gplot]

        msd_dat = os.path.join(args.output_dir, f"{stem}_msd.dat")
        msd_gplot = os.path.join(args.output_dir, f"{stem}_msd.gplot")
        series_labels = species_set + ["all"]
        t_common = msd_results[series_labels[0]][0]
        msd_cols = [t_common] + [msd_results[lbl][1] for lbl in series_labels]
        msd_header = "t(fs) " + " ".join(f"MSD_{lbl}(Ang^2)" for lbl in series_labels)
        msd_series_spec = [(i + 2, lbl) for i, lbl in enumerate(series_labels)]
        write_xy_data(msd_dat, msd_header, msd_cols)
        write_xy_gplot(msd_gplot, msd_dat, f"{stem}_msd.pdf", "t (fs)", "MSD (Ang^2)",
                       f"MSD -- {stem}", msd_series_spec)
        written_files += [msd_dat, msd_gplot]

        if vacf_result is not None:
            freq_cm1, vdos, t_fs_vacf, vacf = vacf_result
            vacf_dat = os.path.join(args.output_dir, f"{stem}_vacf.dat")
            vacf_gplot = os.path.join(args.output_dir, f"{stem}_vacf.gplot")
            write_xy_data(vacf_dat, "t(fs) VACF(normalized)", [t_fs_vacf, vacf])
            write_xy_gplot(vacf_gplot, vacf_dat, f"{stem}_vacf.pdf", "t (fs)", "VACF (normalized)",
                           f"Velocity autocorrelation -- {stem}", [(2, "VACF")])
            written_files += [vacf_dat, vacf_gplot]

            vdos_dat = os.path.join(args.output_dir, f"{stem}_vdos.dat")
            vdos_gplot = os.path.join(args.output_dir, f"{stem}_vdos.gplot")
            write_xy_data(vdos_dat, "freq(cm^-1) VDOS(arb.units)", [freq_cm1, vdos])
            write_xy_gplot(vdos_gplot, vdos_dat, f"{stem}_vdos.pdf", "Frequency (cm^-1)",
                           "VDOS (arb. units)", f"Vibrational DOS -- {stem}", [(2, "VDOS")])
            written_files += [vdos_dat, vdos_gplot]

        if atom_disp_result is not None:
            t_disp, traj_atom, disp, disp_mag, species_tracked = atom_disp_result
            disp_dat = os.path.join(args.output_dir, f"{stem}_disp_atom{args.track_atom}.dat")
            disp_gplot = os.path.join(args.output_dir, f"{stem}_disp_atom{args.track_atom}.gplot")
            write_xy_data(disp_dat, "t(fs) dx(Ang) dy(Ang) dz(Ang) |d|(Ang)",
                         [t_disp, disp[:, 0], disp[:, 1], disp[:, 2], disp_mag])
            write_xy_gplot(disp_gplot, disp_dat, f"{stem}_disp_atom{args.track_atom}.pdf",
                           "t (fs)", "Displacement (Ang)",
                           f"Atom {args.track_atom} ({species_tracked}) displacement -- {stem}",
                           [(2, "dx"), (3, "dy"), (4, "dz"), (5, "|d|")])
            written_files += [disp_dat, disp_gplot]

        if pair_dist_result is not None:
            t_pair, vec_pair, dist_pair, i_idx, j_idx = pair_dist_result
            dist_dat = os.path.join(args.output_dir, f"{stem}_dist_{i_idx}_{j_idx}.dat")
            dist_gplot = os.path.join(args.output_dir, f"{stem}_dist_{i_idx}_{j_idx}.gplot")
            write_xy_data(dist_dat, "t(fs) dx(Ang) dy(Ang) dz(Ang) |r|(Ang)",
                         [t_pair, vec_pair[:, 0], vec_pair[:, 1], vec_pair[:, 2], dist_pair])
            write_xy_gplot(dist_gplot, dist_dat, f"{stem}_dist_{i_idx}_{j_idx}.pdf",
                           "t (fs)", "Distance (Ang)",
                           f"Atom {i_idx}-{j_idx} relative distance -- {stem}",
                           [(2, "dx"), (3, "dy"), (4, "dz"), (5, "|r|")])
            written_files += [dist_dat, dist_gplot]

        thermo_panels = []
        vol_dat = os.path.join(args.output_dir, f"{stem}_volume.dat")
        write_xy_data(vol_dat, "t(fs) Volume(Ang^3)", [t_thermo, volumes])
        written_files.append(vol_dat)
        thermo_panels.append((vol_dat, "t (fs)", "Volume (Ang^3)", f"Volume -- {stem}", 2))

        if e_tot is not None or e_pot is not None:
            # Per-atom columns come FIRST (the plotted/primary series, both
            # here and in the matplotlib panel below) -- the y-axis scale a
            # reader actually wants is eV/atom, so a system-size-independent
            # value is comparable across different structures/runs; the raw
            # absolute totals are still written as extra columns afterward
            # for completeness, not for plotting.
            energy_cols = [t_thermo]
            energy_header_parts = ["t(fs)"]
            energy_series = []
            col_i = 2
            for label_, arr, tag in [("E_tot", e_tot, "Energy_total(eV_per_atom)"),
                                      ("E_pot", e_pot, "Energy_potential(eV_per_atom)")]:
                if arr is not None:
                    energy_cols.append(arr / natoms)
                    energy_header_parts.append(tag)
                    energy_series.append((col_i, label_))
                    col_i += 1
            for label_, arr, tag in [("E_tot", e_tot, "Energy_total(eV)"),
                                      ("E_pot", e_pot, "Energy_potential(eV)")]:
                if arr is not None:
                    energy_cols.append(arr)
                    energy_header_parts.append(tag)
            energy_dat = os.path.join(args.output_dir, f"{stem}_energy.dat")
            write_xy_data(energy_dat, " ".join(energy_header_parts), energy_cols)
            written_files.append(energy_dat)
            thermo_panels.append((energy_dat, "t (fs)", "Energy (eV/atom)", f"Energy -- {stem}",
                                  energy_series[0][0]))
            if len(energy_series) > 1:
                energy_gplot = os.path.join(args.output_dir, f"{stem}_energy.gplot")
                write_xy_gplot(energy_gplot, energy_dat, f"{stem}_energy.pdf", "t (fs)",
                               "Energy (eV/atom)", f"Energy -- {stem}", energy_series)
                written_files.append(energy_gplot)

        if temperature is not None:
            temp_dat = os.path.join(args.output_dir, f"{stem}_temperature.dat")
            write_xy_data(temp_dat, "t(fs) T(K)", [t_thermo, temperature])
            written_files.append(temp_dat)
            thermo_panels.append((temp_dat, "t (fs)", "Temperature (K)",
                                  f"Temperature -- {stem}", 2))

        if pressure is not None:
            press_dat = os.path.join(args.output_dir, f"{stem}_pressure.dat")
            write_xy_data(press_dat, "t(fs) P(kBar)", [t_thermo, pressure])
            written_files.append(press_dat)
            thermo_panels.append((press_dat, "t (fs)", "Pressure (kBar)",
                                  f"Pressure -- {stem}", 2))

        if len(thermo_panels) > 1:
            thermo_gplot = os.path.join(args.output_dir, f"{stem}_thermo.gplot")
            write_multipanel_gplot(thermo_gplot, f"{stem}_thermo.pdf", thermo_panels)
            written_files.append(thermo_gplot)

        print_dual(color_text(
            f"[OK] Data + gnuplot scripts written: {len(written_files)} file(s).", 'green'), f_out)
        for path in written_files:
            print_dual(f"  - {path}", f_out)
    else:
        print_dual("Not written (off by default -- pass --save-gnuplot to write "
                   f"{stem}_{{rdf,msd,vacf,vdos,volume,energy,temperature,pressure}}.dat/.gplot "
                   "(whichever are available) plus a combined 4-panel {stem}_thermo.gplot, plus "
                   "the atom-displacement/pair-distance pair if --track-atom/--track-pair are "
                   "used).", f_out)

    print_section("[9] REFERENCES", f_out)
    bib_path = None
    if args.label:
        bib_path = os.path.join(args.output_dir, BIB_FILE)
        bib_entries = [citations.SIESTA, citations.SIESTA_RECENT]
        citations.write_bib_file(bib_path, bib_entries)
        print_dual(color_text(
            f"[OK] Citations for the methods used in this run written to '{bib_path}' "
            f"({len(bib_entries)} entries).", 'green'), f_out)
    else:
        print_dual("No SIESTA-specific references for a generic --trajectory input -- cite "
                   "whichever tool produced it instead (e.g. stb-mlmd's own MACE/foundation-"
                   "model references).", f_out)

    print_section("[10] SUMMARY & FILES", f_out)
    print_dual("Status            : OK", f_out)
    if valid.any():
        print_dual(f"RDF first peak    : r = {r[rdf_peak_idx]:.3f} Ang", f_out)
    if vacf_result is not None:
        freq_cm1, vdos, _, _ = vacf_result
        peak_idx_v = np.argmax(vdos[1:]) + 1 if len(vdos) > 1 else 0
        print_dual(f"VDOS dominant peak: {freq_cm1[peak_idx_v]:.1f} cm^-1", f_out)
    if atom_disp_result is not None:
        print_dual(f"Atom {args.track_atom} net displacement: {atom_disp_result[3][-1]:.4f} Ang", f_out)
    if pair_dist_result is not None:
        print_dual(f"Pair {track_pair[0]}-{track_pair[1]} distance: "
                   f"{pair_dist_result[2][0]:.4f} -> {pair_dist_result[2][-1]:.4f} Ang", f_out)
    print_dual(f"Volume            : mean = {volumes.mean():.4f} Ang^3 "
               f"(std = {volumes.std():.4f})", f_out)
    if e_tot is not None:
        print_dual(f"Energy (total)    : mean = {e_tot.mean():.4f} eV "
                   f"(std = {e_tot.std():.4f}), {e_tot.mean() / natoms:.6f} eV/atom", f_out)
    if temperature is not None:
        print_dual(f"Temperature       : mean = {temperature.mean():.2f} K "
                   f"(std = {temperature.std():.2f})", f_out)
    if pressure is not None:
        print_dual(f"Pressure          : mean = {pressure.mean():.3f} kBar "
                   f"(std = {pressure.std():.3f})", f_out)
    for path in written_files:
        print_dual(f"  - {path}", f_out)
    if bib_path:
        print_dual(f"References        : {bib_path}", f_out)
    if report_path:
        print_dual(f"Report            : {report_path}", f_out)

    if f_out:
        f_out.close()

    # --view runs last, after the report is fully printed/closed, so a
    # blocking matplotlib window never delays or hides it.
    if args.view:
        import matplotlib.pyplot as plt

        pair_tag_view = f" ({pair[0]}-{pair[1]})" if pair else ""
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.plot(r, g_r)
        ax.set_xlabel("r (Ang)")
        ax.set_ylabel("g(r)")
        ax.set_title(f"RDF{pair_tag_view} -- {stem}")
        ax.axhline(1.0, color='gray', linestyle='--', linewidth=0.8)
        fig.tight_layout()

        fig, ax = plt.subplots(figsize=(7, 4.5))
        for lbl in species_set + ["all"]:
            t_fs, msd, _ = msd_results[lbl]
            ax.plot(t_fs, msd, label=lbl)
        ax.set_xlabel("t (fs)")
        ax.set_ylabel("MSD (Ang^2)")
        ax.set_title(f"MSD -- {stem}")
        ax.legend()
        fig.tight_layout()

        if vacf_result is not None:
            freq_cm1, vdos, t_fs_vacf, vacf = vacf_result
            fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
            axes[0].plot(t_fs_vacf, vacf)
            axes[0].set_xlabel("t (fs)")
            axes[0].set_ylabel("VACF (normalized)")
            axes[0].set_title("Velocity autocorrelation")
            axes[1].plot(freq_cm1, vdos)
            axes[1].set_xlabel("Frequency (cm^-1)")
            axes[1].set_ylabel("VDOS (arb. units)")
            axes[1].set_title("Vibrational density of states")
            fig.suptitle(stem)
            fig.tight_layout()

        if atom_disp_result is not None:
            t_disp, _, disp, disp_mag, species_tracked = atom_disp_result
            fig, ax = plt.subplots(figsize=(7, 4.5))
            ax.plot(t_disp, disp[:, 0], label="dx")
            ax.plot(t_disp, disp[:, 1], label="dy")
            ax.plot(t_disp, disp[:, 2], label="dz")
            ax.plot(t_disp, disp_mag, label="|d|", color='k', linewidth=2)
            ax.set_xlabel("t (fs)")
            ax.set_ylabel("Displacement (Ang)")
            ax.set_title(f"Atom {args.track_atom} ({species_tracked}) displacement -- {stem}")
            ax.legend()
            fig.tight_layout()

        if pair_dist_result is not None:
            t_pair, vec_pair, dist_pair, i_idx, j_idx = pair_dist_result
            fig, ax = plt.subplots(figsize=(7, 4.5))
            ax.plot(t_pair, vec_pair[:, 0], label="dx", alpha=0.6)
            ax.plot(t_pair, vec_pair[:, 1], label="dy", alpha=0.6)
            ax.plot(t_pair, vec_pair[:, 2], label="dz", alpha=0.6)
            ax.plot(t_pair, dist_pair, label="|r|", color='k', linewidth=2)
            ax.set_xlabel("t (fs)")
            ax.set_ylabel("Distance (Ang)")
            ax.set_title(f"Atom {i_idx}-{j_idx} relative distance -- {stem}")
            ax.legend()
            fig.tight_layout()

        fig, axes = plt.subplots(2, 2, figsize=(11, 8))
        axes[0, 0].plot(t_thermo, volumes, color='#1f77b4')
        axes[0, 0].set_xlabel("t (fs)")
        axes[0, 0].set_ylabel("Volume (Ang^3)")
        axes[0, 0].set_title("Volume")

        # Plotted per-atom (not the raw absolute total) -- a system-size
        # -independent y-axis scale, same reasoning as the .dat/.gplot output
        # above; the report table [7] still prints the absolute eV values too.
        if e_tot is not None or e_pot is not None:
            if e_tot is not None:
                axes[0, 1].plot(t_thermo, e_tot / natoms, label="E_tot")
            if e_pot is not None:
                axes[0, 1].plot(t_thermo, e_pot / natoms, label="E_pot")
            axes[0, 1].legend()
        else:
            axes[0, 1].text(0.5, 0.5, "not available", ha='center', va='center',
                            transform=axes[0, 1].transAxes, color='gray')
        axes[0, 1].set_xlabel("t (fs)")
        axes[0, 1].set_ylabel("Energy (eV/atom)")
        axes[0, 1].set_title("Energy")

        if temperature is not None:
            axes[1, 0].plot(t_thermo, temperature, color='#d62728')
        else:
            axes[1, 0].text(0.5, 0.5, "not available", ha='center', va='center',
                            transform=axes[1, 0].transAxes, color='gray')
        axes[1, 0].set_xlabel("t (fs)")
        axes[1, 0].set_ylabel("Temperature (K)")
        axes[1, 0].set_title("Temperature (thermostat)")

        if pressure is not None:
            axes[1, 1].plot(t_thermo, pressure, color='#2ca02c')
        else:
            axes[1, 1].text(0.5, 0.5, "not available", ha='center', va='center',
                            transform=axes[1, 1].transAxes, color='gray')
        axes[1, 1].set_xlabel("t (fs)")
        axes[1, 1].set_ylabel("Pressure (kBar)")
        axes[1, 1].set_title("Pressure")

        fig.suptitle(f"Thermodynamic time series -- {stem}")
        fig.tight_layout()

        plt.show()


if __name__ == "__main__":
    main()
