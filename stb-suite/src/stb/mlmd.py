#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

"""stb-mlmd (5.1): molecular dynamics driven by a MACE potential.

The report is built in numbered sections (run, structure, model and its card, MD setup,
equilibration, production, thermodynamics, ensemble checks, structure along the run and the
model's trained domain, summary) with every warning collected in one final section and the
messages of MACE/torch only counted. Plots are opt-in: --png/--show (a matplotlib dashboard)
and --save-gnuplot (.dat files + one .gplot); --save-data writes the .dat files alone.

NVE energy conservation is measured as a drift in meV/atom/ps (a linear fit of E_tot(t)): the
relative drift used before divided by the absolute total energy, which on SIESTA's scale
(-24826 eV for 216 Si atoms) is so large that the 1 % warning would need a 248 eV drift.
"""

from stb import __version__ as VERSION

import os
import sys
import time
import argparse
from datetime import datetime

import numpy as np

from stb.core import structure_io, kspace
from stb.core import model_card
from stb.core.cli import (color_text, show_intro, print_dual, print_section, print_table,
                          capture_library_noise)

REPORT_FILE = "stb_mlmd_report.txt"

# Same 3 formats/extensions stb-ani2traj already writes a SIESTA AIMD
# trajectory into -- reused here so an MLMD trajectory opens in the same
# viewers (OVITO/VMD) the same way, no separate format story for the suite.
OUTPUT_FORMATS = {
    "xsf": ("xsf", ".xsf"),
    "pdb": ("proteindatabank", ".pdb"),
    "xyz": ("extxyz", ".xyz"),
}

# A good NVE run at 1 fs drifts well below this (same line as stb-mlffAnalysis's stability check).
NVE_DRIFT_WARN = 1.0            # meV/atom/ps
WATER_COMPRESSIBILITY = 4.57e-5  # 1/bar, the --compressibility placeholder default
N_BLOCKS = 5                    # block averaging of the production series
RDF_MAX_FRAMES = 100
RDF_DR = 0.02                   # Ang
MSD_MAX_LAGS = 200
EV_TO_MEV = 1000.0
AMU_A3_TO_G_CM3 = 1.66053906660
KB_EV = 8.617333262e-5


def build_dynamics(atoms, ensemble, timestep_fs, temperature_k, friction_fs, pressure_gpa, seed,
                   taut_fs=100.0, taup_fs=500.0, compressibility_bar=4.57e-5):
    """Returns a configured ASE MD driver, already given Maxwell-Boltzmann
    initial velocities (net momentum removed) at `temperature_k`. Timestep/
    friction/pressure are given in physically intuitive units (fs, fs^-1,
    GPa) and converted to ASE's internal unit system here -- NOT done by
    simply passing the raw numbers through (a real, verified bug in
    stb-amorphize's own MD call: passing timestep in fs directly as ASE
    time units silently runs at timestep/units.fs =~ 10.18x the intended
    step size; fixed there alongside adding this tool, see amorphize.py's
    own comment).

    `compressibility_bar` (1/bar, NPT only) defaults to water's value (the
    same generic placeholder stb-amorphize uses) -- NOT material-specific;
    --compressibility lets a caller override it with something closer to the
    real material's 1/bulk-modulus (1/(B[GPa] * 1e4) in 1/bar) for a
    better-calibrated barostat response. It is converted to ASE's internal
    1/(eV/Ang^3) here: passed raw (as this function once did, labelled
    eV/Ang^3), 4.57e-5 is ~1.6e6 x too small and the barostat never moves the
    cell -- an 8-atom Si cell at 10 GPa changed its volume by 1e-6.
    """
    from ase import units
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary

    rng = np.random.default_rng(seed) if seed is not None else None
    MaxwellBoltzmannDistribution(atoms, temperature_K=temperature_k, rng=rng)
    Stationary(atoms)

    dt = timestep_fs * units.fs
    if ensemble == "nve":
        from ase.md.verlet import VelocityVerlet
        return VelocityVerlet(atoms, timestep=dt)
    if ensemble == "nvt":
        from ase.md.langevin import Langevin
        return Langevin(atoms, timestep=dt, temperature_K=temperature_k,
                         friction=friction_fs / units.fs, rng=rng)
    from ase.md.nptberendsen import NPTBerendsen
    return NPTBerendsen(atoms, timestep=dt, temperature_K=temperature_k,
                        pressure_au=pressure_gpa * units.GPa,
                        taut=taut_fs * units.fs, taup=taup_fs * units.fs,
                        compressibility_au=compressibility_bar / units.bar)


# ---- statistics -------------------------------------------------------------------

def block_se(x, n_blocks=N_BLOCKS):
    """Standard error of the mean of a correlated series by block averaging (None if too short)."""
    x = np.asarray(x, float)
    if len(x) < 2 * n_blocks:
        return None
    means = [b.mean() for b in np.array_split(x, n_blocks)]
    return float(np.std(means, ddof=1) / np.sqrt(n_blocks))


def slope_per_ps(t_fs, y):
    """Linear-fit slope of y(t) per ps."""
    if len(t_fs) < 3:
        return None
    return float(np.polyfit(np.asarray(t_fs) / 1000.0, np.asarray(y, float), 1)[0])


def thirds_shift(x):
    """(mean of the last third - mean of the first third, combined block SE) -- a drift test."""
    x = np.asarray(x, float)
    n3 = len(x) // 3
    if n3 < N_BLOCKS:
        return None, None
    a, b = x[:n3], x[-n3:]
    sa, sb = block_se(a), block_se(b)
    se = np.hypot(sa, sb) if sa is not None and sb is not None else None
    return float(b.mean() - a.mean()), se


def fmt(v, spec="{:.2f}", none="-"):
    return none if v is None or (isinstance(v, float) and not np.isfinite(v)) else spec.format(v)


def pm(mean, err, spec="{:.2f}"):
    return f"{spec.format(mean)} +/- {spec.format(err)}" if err is not None else spec.format(mean)


# ---- structure analysis -----------------------------------------------------------

def perpendicular_widths(cell):
    cell = np.asarray(cell)
    vol = abs(np.linalg.det(cell))
    return [vol / np.linalg.norm(np.cross(cell[(i + 1) % 3], cell[(i + 2) % 3])) for i in range(3)]


def compute_rdf(frames, max_frames=RDF_MAX_FRAMES, dr=RDF_DR):
    """Total radial distribution function averaged over up to `max_frames` frames of a bulk cell:
    dict(r, g, n_r (running coordination), r_peak, g_peak, r_min, coordination, n_frames, rmax)."""
    from ase.neighborlist import neighbor_list
    sel = frames[::max(1, len(frames) // max_frames)]
    rmax = min(8.0, 0.5 * min(min(perpendicular_widths(a.cell)) for a in sel) - 1e-6)
    if rmax < 2.0:
        return None
    edges = np.arange(0.0, rmax + dr, dr)
    hist = np.zeros(len(edges) - 1)
    rho = 0.0
    n = len(sel[0])
    for a in sel:
        d = neighbor_list("d", a, rmax)
        hist += np.histogram(d, bins=edges)[0]
        rho += n / a.get_volume()
    rho /= len(sel)
    r = 0.5 * (edges[1:] + edges[:-1])
    shell = 4.0 / 3.0 * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
    g = hist / (len(sel) * n * rho * shell)
    n_r = np.cumsum(hist) / (len(sel) * n)
    out = {"r": r, "g": g, "n_r": n_r, "n_frames": len(sel), "rmax": rmax,
           "r_peak": None, "g_peak": None, "r_min": None, "coordination": None}
    if g.max() <= 0:
        return out
    ip = int(np.argmax(g))
    out["r_peak"], out["g_peak"] = float(r[ip]), float(g[ip])
    # first minimum after the first peak, on a lightly smoothed curve
    gs = np.convolve(g, np.ones(5) / 5, mode="same")
    hi = min(len(r) - 1, int(np.searchsorted(r, 1.6 * r[ip])))
    if hi > ip + 2:
        im = ip + 1 + int(np.argmin(gs[ip + 1:hi]))
        out["r_min"], out["coordination"] = float(r[im]), float(n_r[im])
    return out


def compute_msd(positions, cells, symbols, dt_fs, max_lags=MSD_MAX_LAGS):
    """Mean-square displacement with every time origin, from the unwrapped saved frames (centre-of-
    mass drift removed), overall and per element; diffusion coefficients from a linear fit over the
    second half of the lag range; Lindemann-like ratio sqrt(<|u - <u>|^2>) / d_nn."""
    from stb.core.md_traj import unwrap_trajectory
    nf = len(positions)
    if nf < 4:
        return None
    u = np.array(unwrap_trajectory([np.asarray(p) for p in positions], [np.asarray(c) for c in cells]))
    u -= u.mean(axis=1, keepdims=True) - u[:1].mean(axis=1, keepdims=True)
    lags = np.unique(np.linspace(1, nf // 2, min(max_lags, nf // 2)).astype(int))
    symbols = np.asarray(symbols)
    elements = list(dict.fromkeys(symbols))
    msd = {"all": [], **{el: [] for el in elements}}
    for lag in lags:
        d2 = ((u[lag:] - u[:-lag]) ** 2).sum(axis=2)        # (origins, atoms)
        msd["all"].append(float(d2.mean()))
        for el in elements:
            msd[el].append(float(d2[:, symbols == el].mean()))
    t = lags * dt_fs
    half = t >= t[len(t) // 2] if len(t) > 3 else np.ones(len(t), bool)
    diff = {}
    for key, vals in msd.items():
        vals = np.asarray(vals)
        if half.sum() >= 3 and len(t) >= 6:
            slope = np.polyfit(t[half], vals[half], 1)[0]            # Ang^2/fs
            diff[key] = float(slope / 6.0 * 0.1)                     # cm^2/s (1 Ang^2/fs = 0.1 cm^2/s)
        else:
            diff[key] = None
    mean_pos = u.mean(axis=0)
    rms_vib = float(np.sqrt(((u - mean_pos) ** 2).sum(axis=2).mean()))
    total = u[-1] - u[0]
    return {"t_fs": t, "msd": {k: np.asarray(v) for k, v in msd.items()}, "D_cm2_s": diff,
            "rms_about_mean": rms_vib, "rms_total": float(np.sqrt((total ** 2).sum(axis=1).mean())),
            "max_total": float(np.sqrt((total ** 2).sum(axis=1)).max()), "elements": elements}


def nearest_neighbour_distance(atoms):
    from ase.neighborlist import neighbor_list
    try:
        i, d = neighbor_list("id", atoms, 4.5)
    except Exception:
        return None
    if not len(d):
        return None
    per_atom = np.full(len(atoms), np.inf)
    np.minimum.at(per_atom, i, d)
    per_atom = per_atom[np.isfinite(per_atom)]
    return float(np.median(per_atom)) if len(per_atom) else None


# ---- model card -------------------------------------------------------------------

def _card_value(d, *keys):
    for k in keys:
        if not isinstance(d, dict):
            return None
        d = d.get(k)
    return d


def print_model_card(card, card_path, model_path, f_out, warn):
    """The card that came with --custom-model, section by section."""
    m = card.get("model") or {}
    com = card.get("committee") or {}
    arch = m.get("architecture") or {}
    lt = card.get("level_of_theory") or {}
    td = card.get("training_data") or {}
    dom = card.get("domain") or {}
    rec = card.get("recipe") or {}
    met = card.get("metrics") or {}
    val = card.get("validation") or {}
    chk = card.get("checks") or {}
    use = card.get("usage") or {}
    prov = card.get("provenance") or {}
    sel = card.get("selection") or {}

    def sub(title):
        print_dual("\n" + color_text(title, 'bold'), f_out)

    sub("Model card -- identity")
    digest = model_card.sha256_file(model_path)
    matched = next((f for f in m.get("files", []) if f.get("sha256") == digest), {})
    print_dual(f"Card file         : {card_path} (format {card.get('format')} v{card.get('version')}, "
               f"written {card.get('created', '?')}, updated {card.get('updated', '-')}, stb {card.get('stb_version', '?')})", f_out)
    print_dual(f"Model             : {m.get('name', '?')} -- {model_path}", f_out)
    print_dual(f"SHA-256           : {digest[:16]}... = card's '{matched.get('file', '?')}' "
               f"({matched.get('kind', '?')}) [verified]", f_out)
    print_dual(f"Committee         : member {m.get('member', '?')} of {com.get('size', '?')} (seed {m.get('seed', '?')}), "
               f"round {com.get('round', '?')}"
               + (f"; committee force spread mean {com['spread']['F_spread_mean_meV']:.1f} meV/Ang"
                  if (com.get("spread") or {}).get("F_spread_mean_meV") is not None else ""), f_out)
    if sel:
        print_dual(f"Selected as best  : by {sel.get('selected_by', '?')} on {sel.get('time', '?')} -- "
                   f"{sel.get('basis') or sel.get('criterion', '?')}", f_out)
    found = m.get("foundation") or {}
    print_dual(f"Architecture      : MACE, r_max {fmt(arch.get('r_max'), '{:g}')} Ang, "
               f"{fmt(arch.get('num_interactions'), '{}')} interactions, {fmt(arch.get('n_params'), '{:,}')} parameters; "
               f"fine-tuned from {found.get('name', '?')} ({found.get('size', '?')}), start: {m.get('init_from', '?')}", f_out)

    sub("Model card -- level of theory (what the model reproduces)")
    s = lt.get("settings") or {}
    print_dual(f"Code              : {lt.get('code', 'SIESTA')} {lt.get('siesta_version') or ''}".rstrip(), f_out)
    keys = list(s)
    for i in range(0, len(keys), 3):
        print_dual(("Settings          : " if i == 0 else "                    ")
                   + "; ".join(f"{k} {s[k]}" for k in keys[i:i + 3]), f_out)
    d3 = lt.get("dftd3") or {}
    print_dual(f"DFT-D3            : {'included in the labels' if d3.get('on') else 'not in the labels'}"
               f" (source: {d3.get('source', '?')})", f_out)
    pp = lt.get("pseudopotentials_md5") or {}
    if pp:
        print_dual("Pseudopotentials  : " + ", ".join(f"{k} (md5 {v})" for k, v in pp.items()), f_out)
    print_dual(f"Energy scale      : {_card_value(card, 'energy', 'scale') or 'SIESTA absolute total energy'}; E0 "
               + ", ".join(f"{el} {e:.4f} eV" for el, e in (_card_value(card, 'energy', 'e0_eV') or {}).items())
               + f" ({_card_value(card, 'energy', 'e0_source') or '?'})", f_out)

    sub("Model card -- training data and trained domain")
    kinds = td.get("kinds") or {}
    apc = td.get("atoms_per_config") or {}
    print_dual(f"Training set      : {td.get('n_configs', '?')} SIESTA configs ({td.get('n_train', '?')} train / "
               f"{td.get('n_valid', '?')} valid), kinds " + ", ".join(f"{k} {v}" for k, v in kinds.items())
               + f"; {apc.get('min', '?')}-{apc.get('max', '?')} atoms per config", f_out)
    lab = td.get("labels") or {}
    print_dual("Labels            : " + ", ".join(k for k in ("energy", "forces", "stress") if lab.get(k))
               + ("; stress fitted" if lab.get("stress_used_in_training") else ""), f_out)
    rows = []
    for name, key, spec in [(dom.get("measure_name") or "cell size", "measure", "{:.3f}"),
                            ("shortest distance (Ang)", "d_min_ang", "{:.3f}"),
                            ("max force per config (eV/Ang)", "force_max_eV_per_ang", "{:.2f}"),
                            ("pressure (GPa)", "pressure_GPa", "{:.2f}"),
                            ("energy per atom (eV)", "energy_per_atom_eV", "{:.4f}"),
                            ("atoms per config", "n_atoms", "{}")]:
        r = dom.get(key)
        if r:
            rows.append(([name, f"{spec.format(r['min'])} .. {spec.format(r['max'])}"], None))
    rows.insert(0, (["elements / dimensionality", f"{', '.join(dom.get('elements') or card.get('elements') or [])} / "
                     f"{dom.get('dimensionality', '?')}"], None))
    print_table(["Trained range", "Value"], rows, f_out)

    sub("Model card -- recipe and accuracy")
    print_dual(f"Recipe            : {rec.get('epochs_max', '?')} epochs (kept {rec.get('kept_epoch', '?')}), batch "
               f"{rec.get('batch_size', '?')}, lr {rec.get('lr', '?')}, loss {rec.get('loss', '?')}, E0s {rec.get('e0s', '?')}, "
               f"{rec.get('device', '?')} {fmt(rec.get('training_seconds'), '{:.0f}')} s", f_out)
    rows = []
    for label, d in [("train", met.get("train")), ("valid (own split)", met.get("valid")),
                     ("out-of-sample (committee)", met.get("out_of_sample_committee")),
                     ("test set (stb-mlffValidate)", val.get("overall")),
                     ("foundation model, raw", met.get("foundation_model"))]:
        if d:
            rows.append(([label, fmt(d.get("n"), "{}"), fmt(d.get("F_rmse_meV"), "{:.1f}"),
                          fmt(d.get("E_rmse_rel_meV"), "{:.2f}"), fmt(d.get("S_rmse_GPa"), "{:.3f}")], None))
    print_table(["Errors vs SIESTA", "N", "F RMSE (meV/Ang)", "E rel (meV/atom)", "S RMSE (GPa)"], rows, f_out)
    if met.get("fine_tuning_gain_pct") is not None:
        print_dual(f"Fine-tuning gain  : {met['fine_tuning_gain_pct']:.0f} % lower force error than the raw foundation model", f_out)
    for g, d in (val.get("groups") or {}).items():
        print_dual(f"  test, {g:<10}: F {fmt(d.get('F_rmse_meV'), '{:.1f}')} meV/Ang over {d.get('n', '?')} configs", f_out)
    ph = val.get("phonons") or {}
    if ph.get("rmse_THz") is not None:
        print_dual(f"Phonons vs SIESTA : RMSE {ph['rmse_THz']:.2f} THz (acoustic {fmt(ph.get('acoustic_rmse_THz'))}, optical "
                   f"{fmt(ph.get('optical_rmse_THz'))}), imaginary q-points: model {ph.get('model_imaginary_q', '?')}", f_out)
    for run, pts in (val.get("md_check") or {}).items():
        if pts:
            f_vals = [p[1] for p in pts if p and p[1] is not None]
            print_dual(f"MD check {run:<9}: F {min(f_vals):.1f}-{max(f_vals):.1f} meV/Ang along {len(f_vals)} SIESTA-labelled frames", f_out)
    cal = val.get("calibration") or {}
    if cal.get("error_over_spread") is not None:
        print_dual(f"Committee spread  : true error ~ {cal['error_over_spread']:.1f} x spread "
                   f"(per-config correlation {fmt(cal.get('config_corr'), '{:+.2f}')})", f_out)
    if val.get("verdict"):
        print_dual(f"Verdict (4.17.3)  : {str(val['verdict']).upper()} (validated {val.get('time', '?')}, "
                   f"{val.get('n_test', '?')} test configs)", f_out)

    sub("Model card -- stability, speed and usage")
    st = chk.get("stability") or {}
    sp = chk.get("speed") or {}
    if st:
        print_dual(f"Stability (4.17.2): {st.get('status', '?')} -- {st.get('n_atoms', '?')} atoms, T mean "
                   f"{fmt(st.get('T_mean'), '{:.0f}')} K, NVE drift {fmt(st.get('drift_meV_atom_ps'), '{:.3f}')} meV/atom/ps", f_out)
    if sp:
        print_dual(f"Speed (4.17.2)    : {fmt(sp.get('ms_per_step'), '{:.1f}')} ms/step for {sp.get('n_atoms', '?')} atoms "
                   f"on {sp.get('device', '?')}", f_out)
    rd = use.get("recommended_dtype") or {}
    if rd:
        print_dual("Recommended dtype : " + ", ".join(f"{k.replace('_', ' ')} {v}" for k, v in rd.items()), f_out)
    for note in use.get("notes") or []:
        print_dual(f"  - {note}", f_out)
    sw = prov.get("software") or {}
    if sw:
        print_dual("Trained with      : " + ", ".join(f"{k} {v}" for k, v in sw.items()), f_out)

    # what the card itself says should worry a user
    if not val.get("overall"):
        warn("Model card: the model has no test-set validation (stb-mlffValidate, 4.17.3) -- its accuracy is known "
             "only on its own training data.")
    elif val.get("verdict") and str(val["verdict"]).lower() not in ("stop",):
        warn(f"Model card: the last validation verdict was {str(val['verdict']).upper()}, not STOP -- the workflow "
             "suggested more training data.")
    if st and st.get("status") not in (None, "OK"):
        warn(f"Model card: this member FAILED Stage 2's MD stability test ({st.get('status')}: {st.get('reason') or '-'}).")
    import re
    my_seed = m.get("seed")
    for w in prov.get("training_warnings") or []:
        seeds = {int(x) for x in re.findall(r"\bseed (\d+)", w)}
        if seeds and my_seed is not None and my_seed not in seeds:
            continue                     # a round-level note about another committee member
        warn(f"Model card (recorded at training): {w}")


# ---- plots ------------------------------------------------------------------------

def _show_possible():
    import matplotlib
    if matplotlib.get_backend().lower() in ("agg", "pdf", "ps", "svg", "pgf", "cairo", "template"):
        return False
    return not (sys.platform.startswith("linux") and not (os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY")))


def make_dashboard(ctx, png_path, show):
    """2x3 matplotlib figure of the run: temperature, energies, pressure, cell volume or shortest
    distance, g(r) and MSD. Saved to `png_path` and/or shown."""
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    blue, orange, green, purple, grey = "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#777777"
    t = ctx["t_ps"]
    fig, ax = plt.subplots(2, 3, figsize=(15, 8.2))
    fig.suptitle(ctx["title"], fontsize=11)

    a = ax[0, 0]
    a.plot(t, ctx["T"], color=blue, lw=1.0, label="T(t)")
    if ctx["ensemble"] != "nve":
        a.axhline(ctx["T_target"], color=grey, ls="--", lw=0.8, label=f"target {ctx['T_target']:g} K")
    a.axhline(np.mean(ctx["T"]), color=orange, lw=1.2, label=f"mean {np.mean(ctx['T']):.1f} K")
    a.set_xlabel("time (ps)"); a.set_ylabel("temperature (K)"); a.set_title("Temperature"); a.legend(fontsize=8)

    a = ax[0, 1]
    ep = EV_TO_MEV * (ctx["Epot_atom"] - np.mean(ctx["Epot_atom"]))
    et = EV_TO_MEV * (ctx["Etot_atom"] - np.mean(ctx["Etot_atom"]))
    a.plot(t, ep, color=blue, lw=1.0, label="E_pot")
    a.plot(t, et, color=orange, lw=1.4, label="E_tot = E_pot + E_kin")
    a.set_xlabel("time (ps)"); a.set_ylabel("energy - mean (meV/atom)")
    a.set_title("Energies" + (" (E_tot is conserved in NVE)" if ctx["ensemble"] == "nve" else ""))
    a.legend(fontsize=8)

    a = ax[0, 2]
    if ctx["P"] is not None:
        a.plot(t, ctx["P"], color=green, lw=1.0, label="P(t)")
        a.axhline(np.mean(ctx["P"]), color=orange, lw=1.2, label=f"mean {np.mean(ctx['P']):.2f} GPa")
        if ctx["ensemble"] == "npt":
            a.axhline(ctx["P_target"], color=grey, ls="--", lw=0.8, label=f"target {ctx['P_target']:g} GPa")
        a.set_ylabel("pressure (GPa)"); a.legend(fontsize=8)
    else:
        a.text(0.5, 0.5, "no pressure\n(vacuum-padded cell)", ha="center", va="center", transform=a.transAxes)
    a.set_xlabel("time (ps)"); a.set_title("Pressure (virial + kinetic)")

    a = ax[1, 0]
    if ctx["ensemble"] == "npt":
        a.plot(t, ctx["V_atom"], color=purple, lw=1.0, label="V/atom")
        a.axhline(ctx["V_atom"][0], color=grey, ls="--", lw=0.8, label="start")
        a.set_ylabel("volume per atom (Ang^3)"); a.set_title("Cell volume")
    else:
        a.plot(t, ctx["d_min"], color=purple, lw=1.0, label="shortest distance")
        if ctx.get("d_min_trained") is not None:
            a.axhline(ctx["d_min_trained"], color=orange, ls="--", lw=1.0, label="shortest in training")
        a.set_ylabel("shortest distance (Ang)"); a.set_title("Closest pair")
    a.set_xlabel("time (ps)"); a.legend(fontsize=8)

    a = ax[1, 1]
    rdf = ctx.get("rdf")
    if rdf:
        a.plot(rdf["r"], rdf["g"], color=blue, lw=1.2)
        if rdf["r_peak"]:
            a.axvline(rdf["r_peak"], color=grey, ls=":", lw=0.8)
        if rdf["r_min"]:
            a.axvline(rdf["r_min"], color=orange, ls="--", lw=0.8,
                      label=f"1st min {rdf['r_min']:.2f} Ang, N = {rdf['coordination']:.2f}")
            a.legend(fontsize=8)
        a.set_xlabel("r (Ang)"); a.set_ylabel("g(r)")
    else:
        a.text(0.5, 0.5, "g(r) not computed", ha="center", va="center", transform=a.transAxes)
    a.set_title("Radial distribution function")

    a = ax[1, 2]
    msd = ctx.get("msd")
    if msd:
        colors = [blue, orange, green, purple, grey]
        keys = msd["elements"] if len(msd["elements"]) > 1 else ["all"]
        for i, k in enumerate(keys):
            a.plot(msd["t_fs"] / 1000, msd["msd"][k], color=colors[i % len(colors)], lw=1.2,
                   label=f"{k}: D = {msd['D_cm2_s'][k]:.2e} cm^2/s" if msd["D_cm2_s"][k] is not None else k)
        a.set_xlabel("lag time (ps)"); a.set_ylabel("MSD (Ang^2)"); a.legend(fontsize=8)
    else:
        a.text(0.5, 0.5, "MSD not computed", ha="center", va="center", transform=a.transAxes)
    a.set_title("Mean-square displacement")

    fig.tight_layout()
    if png_path:
        fig.savefig(png_path, dpi=150)
    if show:
        plt.show()
    plt.close(fig)


def write_data_files(stem, ctx):
    """<stem>_md_thermo.dat, <stem>_md_rdf.dat, <stem>_md_msd.dat. Returns the paths written."""
    paths = []
    thermo = f"{stem}_md_thermo.dat"
    n = ctx["n_atoms"]
    with open(thermo, "w") as f:
        f.write(f"# stb-mlmd {ctx['ensemble'].upper()} production run, {n} atoms, energies per atom\n")
        f.write("# 1:time(fs) 2:T(K) 3:Epot(eV/atom) 4:Ekin(eV/atom) 5:Etot(eV/atom) 6:P(GPa) 7:V(Ang^3/atom) "
                "8:a(Ang) 9:b(Ang) 10:c(Ang) 11:d_min(Ang) 12:Fmax(eV/Ang)\n")
        P = ctx["P"] if ctx["P"] is not None else np.full(len(ctx["t_fs"]), np.nan)
        for i in range(len(ctx["t_fs"])):
            f.write(f"{ctx['t_fs'][i]:.3f} {ctx['T'][i]:.3f} {ctx['Epot_atom'][i]:.8f} {ctx['Ekin_atom'][i]:.8f} "
                    f"{ctx['Etot_atom'][i]:.8f} {P[i]:.5f} {ctx['V_atom'][i]:.6f} {ctx['abc'][i][0]:.5f} "
                    f"{ctx['abc'][i][1]:.5f} {ctx['abc'][i][2]:.5f} {ctx['d_min'][i]:.4f} {ctx['Fmax'][i]:.4f}\n")
    paths.append(thermo)
    rdf = ctx.get("rdf")
    if rdf:
        p = f"{stem}_md_rdf.dat"
        with open(p, "w") as f:
            f.write(f"# total g(r) averaged over {rdf['n_frames']} frames; N(r) = running coordination number\n")
            f.write("# 1:r(Ang) 2:g(r) 3:N(r)\n")
            for r, g, nr in zip(rdf["r"], rdf["g"], rdf["n_r"]):
                f.write(f"{r:.4f} {g:.6f} {nr:.6f}\n")
        paths.append(p)
    msd = ctx.get("msd")
    if msd:
        p = f"{stem}_md_msd.dat"
        keys = ["all"] + (msd["elements"] if len(msd["elements"]) > 1 else [])
        with open(p, "w") as f:
            f.write("# mean-square displacement, every time origin, centre-of-mass drift removed\n")
            f.write("# 1:lag(fs) " + " ".join(f"{i + 2}:MSD_{k}(Ang^2)" for i, k in enumerate(keys)) + "\n")
            for j, t in enumerate(msd["t_fs"]):
                f.write(f"{t:.3f} " + " ".join(f"{msd['msd'][k][j]:.6f}" for k in keys) + "\n")
        paths.append(p)
    return paths


def write_gnuplot(stem, ctx):
    """<stem>_md.gplot: the dashboard as a 2x3 gnuplot multiplot -> <stem>_md.pdf."""
    path = f"{stem}_md.gplot"
    thermo, rdf, msd = f"{stem}_md_thermo.dat", f"{stem}_md_rdf.dat", f"{stem}_md_msd.dat"
    L = ['# --- STB Plot Configuration ---\n', '# Generated by stb-mlmd: gnuplot ' + os.path.basename(path) + '\n',
         'set terminal pdfcairo enhanced color font "Arial,10" size 15,8.2\n',
         f'set output "{os.path.basename(stem)}_md.pdf"\n',
         'set multiplot layout 2,3 title "' + ctx["title"].replace("_", "\\\\_") + '"\n', 'set grid\nset key top right font ",8"\n\n',
         'set xlabel "time (ps)"\nset ylabel "temperature (K)"\nset title "Temperature"\n',
         f'plot "{os.path.basename(thermo)}" using ($1/1000):2 with lines lw 1 lc rgb "#0072B2" title "T(t)"'
         + (f', {ctx["T_target"]:g} with lines dt 2 lc rgb "#777777" title "target"\n' if ctx["ensemble"] != "nve" else "\n"),
         f'set ylabel "energy - mean (meV/atom)"\nset title "Energies"\n'
         f'plot "{os.path.basename(thermo)}" using ($1/1000):(($3-({np.mean(ctx["Epot_atom"]):.8f}))*1000) with lines lw 1 lc rgb "#0072B2" title "E_{{pot}}", '
         f'"" using ($1/1000):(($5-({np.mean(ctx["Etot_atom"]):.8f}))*1000) with lines lw 1.5 lc rgb "#D55E00" title "E_{{tot}}"\n',
         'set ylabel "pressure (GPa)"\nset title "Pressure"\n'
         + (f'plot "{os.path.basename(thermo)}" using ($1/1000):6 with lines lw 1 lc rgb "#009E73" title "P(t)"\n'
            if ctx["P"] is not None else 'unset key\nset yrange [-1:1]\nplot 0 with lines lc rgb "white" notitle\nset autoscale y\nset key\n')]
    if ctx["ensemble"] == "npt":
        L.append(f'set ylabel "volume per atom (Ang^3)"\nset title "Cell volume"\n'
                 f'plot "{os.path.basename(thermo)}" using ($1/1000):7 with lines lw 1 lc rgb "#CC79A7" title "V/atom"\n')
    else:
        L.append(f'set ylabel "shortest distance (Ang)"\nset title "Closest pair"\n'
                 f'plot "{os.path.basename(thermo)}" using ($1/1000):11 with lines lw 1 lc rgb "#CC79A7" title "d_{{min}}"'
                 + (f', {ctx["d_min_trained"]:.4f} with lines dt 2 lc rgb "#D55E00" title "shortest in training"\n'
                    if ctx.get("d_min_trained") is not None else "\n"))
    if ctx.get("rdf"):
        L.append(f'set xlabel "r (Ang)"\nset ylabel "g(r)"\nset title "Radial distribution function"\n'
                 f'plot "{os.path.basename(rdf)}" using 1:2 with lines lw 1.5 lc rgb "#0072B2" title "g(r)"\n')
    else:
        L.append('set title "g(r) not computed"\nset yrange [-1:1]\nplot 0 with lines lc rgb "white" notitle\nset autoscale y\n')
    if ctx.get("msd"):
        L.append(f'set xlabel "lag time (ps)"\nset ylabel "MSD (Ang^2)"\nset title "Mean-square displacement"\n'
                 f'plot "{os.path.basename(msd)}" using ($1/1000):2 with lines lw 1.5 lc rgb "#0072B2" title "all atoms"\n')
    else:
        L.append('set title "MSD not computed"\nset yrange [-1:1]\nplot 0 with lines lc rgb "white" notitle\nset autoscale y\n')
    L.append('unset multiplot\n')
    with open(path, "w") as f:
        f.writelines(L)
    return path


# ---- main -------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Molecular dynamics driven by a MACE potential -- the MACE-MP-0 "
                     "foundation model, or a custom model fine-tuned on your own SIESTA "
                     "data via stb-mlffAnalysis (--custom-model, always with its model card). "
                     "NVE (microcanonical), NVT (Langevin thermostat), or NPT (Berendsen "
                     "thermostat+barostat, bulk-only) ensembles. The report covers the model and "
                     "its card, thermodynamic averages with error bars, ensemble checks (NVE drift "
                     "in meV/atom/ps, temperature fluctuations, pressure), the structure along the "
                     "run (shortest distance, g(r), MSD and diffusion, Lindemann ratio) and whether "
                     "the run stayed inside the model's trained range; every warning goes to one "
                     "final section. A fast heuristic exploration tool, not a substitute for a "
                     "real ab initio MD run.",
        epilog="Example usage:\n"
               "  stb-mlmd --file structure.fdf --ensemble nvt --temperature 300 --n-steps 5000\n"
               "  stb-mlmd --file structure.fdf --ensemble nvt --equilibration-steps 200 --n-steps 5000 --png --save-gnuplot\n"
               "  stb-mlmd --file structure.fdf --custom-model mlff_model.model --custom-model-card mlff_model_card.json --ensemble npt \\\n"
               "                           --temperature 500 --pressure 0.0 --compressibility 1e-6 --n-steps 10000 --device cuda",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--file", required=True, help="Input structure (.fdf).")
    parser.add_argument("--model", choices=["small", "medium", "large"], default="medium",
                        help="MACE-MP-0 model size (default: medium). Ignored if --custom-model is given.")
    parser.add_argument("--custom-model", default=None, metavar="PATH",
                        help="Use a custom MACE model file instead of the MACE-MP-0 foundation "
                             "potential -- e.g. one fine-tuned via stb-mlffAnalysis. Overrides --model.")
    model_card.add_card_argument(parser, "--custom-model")
    parser.add_argument("--ensemble", choices=["nve", "nvt", "npt"], default="nvt",
                        help="MD ensemble (default: nvt). npt requires a bulk (no vacuum-padded "
                             "axis) structure -- cell dynamics on a slab/wire/molecule would be "
                             "physically meaningless.")
    parser.add_argument("--temperature", type=float, default=300.0,
                        help="Target temperature, K (default: 300). In NVE only the initial velocities.")
    parser.add_argument("--pressure", type=float, default=0.0,
                        help="Target pressure, GPa (npt only, default: 0.0 -- ambient).")
    parser.add_argument("--timestep", type=float, default=1.0, help="MD timestep, fs (default: 1.0).")
    parser.add_argument("--n-steps", type=int, default=1000,
                        help="Number of PRODUCTION MD steps, after --equilibration-steps (default: 1000).")
    parser.add_argument("--equilibration-steps", type=int, default=0,
                        help="Steps to run and discard BEFORE production (default: 0). Verified live "
                             "that a freshly Maxwell-Boltzmann-seeded structure relaxes for ~100-200 "
                             "steps before settling -- letting that transient into the saved "
                             "trajectory/statistics biases them towards the initial condition.")
    parser.add_argument("--friction", type=float, default=0.01,
                        help="Langevin friction coefficient, fs^-1 (nvt only, default: 0.01).")
    parser.add_argument("--taut", type=float, default=100.0,
                        help="Thermostat coupling time, fs (npt only -- Berendsen, default: 100).")
    parser.add_argument("--taup", type=float, default=500.0,
                        help="Barostat coupling time, fs (npt only, default: 500).")
    parser.add_argument("--compressibility", type=float, default=WATER_COMPRESSIBILITY,
                        help="Compressibility, 1/bar (npt only, default: 4.57e-5 -- water's "
                             "value, a generic placeholder; override with something closer to "
                             "1/bulk-modulus for the actual material, 1/(B[GPa] x 1e4) -- e.g. Si, "
                             "B ~ 98 GPa: 1e-6 -- for a better-calibrated barostat response).")
    parser.add_argument("--stride", type=int, default=10,
                        help="Save one trajectory frame (and one sample of every quantity) every Nth MD "
                             "step (default: 10).")
    parser.add_argument("--vacuum-gap", type=float, default=10.0,
                        help="Gap (Ang) used to detect vacuum-padded axes (default: 10.0, "
                             "matches stb-kgrid/stb-mlrelax).")
    parser.add_argument("-of", "--out-format", choices=sorted(OUTPUT_FORMATS), default="xsf",
                        help="Trajectory format (default: xsf). Same 3 choices as stb-ani2traj.")
    parser.add_argument("-o", "--output", default=None,
                        help="Trajectory filename. Default: <input stem>_md_traj.<ext>")
    parser.add_argument("--output-structure", default="md_final.fdf",
                        help="Final-frame structure output (default: md_final.fdf).")
    parser.add_argument("--device", choices=["cpu", "cuda"], default="cpu",
                        help="Device to run the model on (default: cpu).")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for initial velocities.")
    parser.add_argument("--png", nargs="?", const="", default=None, metavar="FILE",
                        help="Save the run's dashboard (temperature, energies, pressure, volume or shortest "
                             "distance, g(r), MSD) as a PNG (default name: <input stem>_md_dashboard.png). "
                             "Off by default.")
    parser.add_argument("--show", action="store_true",
                        help="Display the dashboard on screen (needs a display). Off by default.")
    parser.add_argument("--save-gnuplot", action="store_true",
                        help="Write the data files (as --save-data) plus <input stem>_md.gplot, the same "
                             "six panels for gnuplot ('gnuplot <stem>_md.gplot' -> <stem>_md.pdf). Off by default.")
    parser.add_argument("--save-data", action="store_true",
                        help="Write the data behind the report: <stem>_md_thermo.dat (time series: T, "
                             "energies per atom, P, volume, cell, shortest distance, max force), "
                             "<stem>_md_rdf.dat and <stem>_md_msd.dat. Off by default.")
    parser.add_argument("--save-report", action="store_true",
                        help=f"Also persist the report to {REPORT_FILE}. Off by default.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    parser.add_argument("-v", "--version", action="version", version=f"stb-mlmd {VERSION}")

    args = parser.parse_args()
    card = model_card.setup_from_args(args.custom_model, args.custom_model_card, "--custom-model", error=parser.error)

    ase_format, default_ext = OUTPUT_FORMATS[args.out_format]
    stem = os.path.splitext(os.path.basename(args.file))[0]
    traj_path = args.output if args.output else f"{stem}_md_traj{default_ext}"
    png_path = (args.png or f"{stem}_md_dashboard.png") if args.png is not None else None

    if args.intro:
        show_intro([
            "Siesta ToolBox Suite - ML Molecular Dynamics",
            "MD driven by MACE-MP-0 or your own fine-tuned potential",
            f"Version {VERSION} | University of Brasilia - 2026",
            "Developed by Dr. Carlos M. O. Bastos"
        ])

    print("\n" + color_text("ML MOLECULAR DYNAMICS:", 'bold'))
    print("-" * 60)

    report_path = REPORT_FILE if args.save_report else None
    f_out = open(report_path, "w") if report_path else None
    warns, library_noise, card_msgs = [], [], []
    t_start = time.monotonic()

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
    with capture_library_noise(library_noise, "import"):
        require_mace()
        from pymatgen.io.ase import AseAtomsAdaptor
        from stb.core import mace_relax
        import torch
        import ase
        import mace

    print_dual(color_text("===== STB-MLMD REPORT =====", 'magenta'), f_out)

    # ---- [0]
    print_section("[0] RUN METADATA", f_out)
    print_dual(f"Date/time         : {datetime.now():%Y-%m-%d %H:%M:%S}", f_out)
    print_dual(f"Software          : stb {VERSION}, mace-torch {getattr(mace, '__version__', '?')}, "
               f"torch {torch.__version__}, ase {ase.__version__}", f_out)
    dev = args.device
    if dev == "cuda":
        try:
            mace_relax.resolve_device("cuda")
            dev_label = f"cuda ({torch.cuda.get_device_name(0)})"
        except Exception as e:
            fail(f"--device cuda: {e}")
    else:
        dev_label = "cpu"
    print_dual(f"Device            : {dev_label}", f_out)
    print_dual(f"Input file        : {args.file}", f_out)
    print_dual(f"Ensemble          : {args.ensemble.upper()}", f_out)
    print_dual(f"Temperature       : {args.temperature:.1f} K" + (" (initial velocities only)" if args.ensemble == "nve" else ""), f_out)
    if args.ensemble == "npt":
        print_dual(f"Pressure          : {args.pressure:.3f} GPa", f_out)
    print_dual(f"Timestep / steps  : {args.timestep} fs / {args.n_steps} "
               f"(+ {args.equilibration_steps} equilibration)", f_out)
    print_dual(f"Model             : "
               f"{'custom (' + args.custom_model + ')' if args.custom_model else f'MACE-MP-0 ({args.model})'}", f_out)
    if report_path:
        print_dual(f"Report file       : {report_path}", f_out)

    # ---- [1]
    print_section("[1] STRUCTURE", f_out)
    if not os.path.isfile(args.file):
        fail(f"Input file '{args.file}' not found.")
    if args.custom_model and not os.path.isfile(args.custom_model):
        fail(f"--custom-model file not found: {args.custom_model}")
    try:
        structure = structure_io.read_fdf(args.file)
    except Exception as e:
        fail(f"Could not read '{args.file}': {e}")

    positions = np.array([pos for _, pos in structure.atoms])
    is_cartesian = structure.coord_format == 'cartesian'
    frac_coords = kspace.to_fractional(positions, structure.lattice, is_cartesian)
    vacuum_axes = kspace.detect_vacuum_axes(frac_coords, structure.lattice, args.vacuum_gap)
    bulk = not any(vacuum_axes)
    print_dual(f"[OK] Read {len(structure.atoms)} atom(s), "
               f"dimensionality: {kspace.dimensionality_label(vacuum_axes)}", f_out)
    if args.ensemble == "npt" and not bulk:
        fail("--ensemble npt requires a bulk (fully periodic) structure -- a vacuum-padded "
             "axis was detected. Cell dynamics on a slab/wire/molecule would be physically "
             "meaningless. Use --ensemble nve/nvt instead.")

    pmg_structure = structure_io.to_pymatgen(structure)
    atoms = AseAtomsAdaptor.get_atoms(pmg_structure)
    n_atoms = len(atoms)
    # An empty slab of >= model_card.VACUUM_GAP (5 Ang) along an axis -- the model card's own criterion -- means
    # the cell is not a condensed bulk even when --vacuum-gap (10 Ang, for the NPT guard) calls it 3D: an O2
    # molecule in a 10 Ang box once got a density, a g(r) peak of 398 and a "solid-like" Lindemann ratio.
    gaps5 = model_card.periodic_measure(atoms)[2]
    condensed = bulk and not any(gaps5)
    cp = atoms.cell.cellpar()
    v0 = atoms.get_volume()
    mass = float(atoms.get_masses().sum())
    d_min0 = model_card.shortest_distance(atoms)
    d_nn0 = nearest_neighbour_distance(atoms)
    comp = {}
    for s_ in atoms.get_chemical_symbols():
        comp[s_] = comp.get(s_, 0) + 1
    print_dual(f"Formula           : {atoms.get_chemical_formula(mode='hill')} ("
               + ", ".join(f"{k} {v}" for k, v in comp.items()) + ")", f_out)
    print_dual(f"Cell              : a {cp[0]:.4f}, b {cp[1]:.4f}, c {cp[2]:.4f} Ang; "
               f"alpha {cp[3]:.2f}, beta {cp[4]:.2f}, gamma {cp[5]:.2f} deg", f_out)
    print_dual(f"Volume            : {v0:.3f} Ang^3 ({v0 / n_atoms:.4f} Ang^3/atom)"
               + (f", density {mass / v0 * AMU_A3_TO_G_CM3:.4f} g/cm^3" if condensed else " (includes vacuum)"), f_out)
    if bulk and not condensed:
        print_dual(f"Empty space       : gaps >= {model_card.VACUUM_GAP:g} Ang along {sum(gaps5)} axis/axes (below "
                   f"--vacuum-gap {args.vacuum_gap:g}) -- treated as a molecule/slab in a box: no density, pressure, "
                   "g(r), heat capacity or Lindemann ratio", f_out)
        if args.ensemble == "npt":
            warn(f"NPT: the cell has empty gaps >= {model_card.VACUUM_GAP:g} Ang along {sum(gaps5)} axis/axes -- a "
                 "barostat on a molecule or slab in a box squeezes the vacuum, not the material.")
    print_dual(f"Shortest distance : {fmt(d_min0, '{:.4f}')} Ang (median nearest neighbour {fmt(d_nn0, '{:.4f}')} Ang)", f_out)
    n_dof = 3 * n_atoms - 3
    print_dual(f"Degrees of freedom: {n_dof} (3N - 3: net momentum removed)", f_out)

    # ---- [2]
    print_section("[2] MODEL", f_out)
    model_arg = args.custom_model if args.custom_model else args.model
    if card is not None:
        model_card.mark_announced(args.custom_model)
        print_model_card(card, args.custom_model_card, args.custom_model, f_out, warn)
        print_dual("", f_out)
        known = set((card.get("domain") or {}).get("elements") or card.get("elements") or [])
        absent = sorted(set(atoms.get_chemical_symbols()) - known) if known else []
        if absent:
            fail(f"the model was trained only on {', '.join(sorted(known))} and cannot evaluate "
                 f"{', '.join(absent)} (see its model card): train a model with those elements in the SIESTA "
                 "data (4.17) or use the MACE-MP-0 foundation model (--model).")
    else:
        print_dual(f"Model             : MACE-MP-0 foundation model, size {args.model} (trained on the Materials "
                   "Project's PBE/PBE+U relaxations)", f_out)
        print_dual("Energy scale      : MACE-MP-0's own (fitted to Materials Project energies) -- not SIESTA's; use "
                   "energy differences only", f_out)
        print_dual("Model card        : none (a foundation model; a fine-tuned model always comes with one)", f_out)
    model_card.collect_warnings(card_msgs)
    t0 = time.monotonic()
    with capture_library_noise(library_noise, "model loading"):
        try:
            atoms.calc = mace_relax.get_calculator(model=model_arg, device=args.device, dtype="float32")
        except Exception as e:
            fail(f"Could not load the model: {e}")
    print_dual(f"Loaded            : in {time.monotonic() - t0:.1f} s, float32 (MD; the card recommends "
               f"{_card_value(card, 'usage', 'recommended_dtype', 'molecular_dynamics') or 'float32'})"
               if card else f"Loaded            : in {time.monotonic() - t0:.1f} s, float32", f_out)

    # ---- [3]
    print_section("[3] MD SETUP", f_out)
    with capture_library_noise(library_noise, "MD setup"):
        dyn = build_dynamics(atoms, args.ensemble, args.timestep, args.temperature,
                             args.friction, args.pressure, args.seed,
                             taut_fs=args.taut, taup_fs=args.taup, compressibility_bar=args.compressibility)
    integ = {"nve": "velocity Verlet (microcanonical)",
             "nvt": f"Langevin, friction {args.friction:g} fs^-1 (~ {1 / args.friction:.0f} fs relaxation), canonical",
             "npt": f"Berendsen thermostat (tau_T {args.taut:g} fs) + barostat (tau_P {args.taup:g} fs)"}[args.ensemble]
    print_dual(f"Integrator        : {integ}", f_out)
    if args.ensemble == "npt":
        b_implied = 1.0 / (args.compressibility * 1e4)
        print_dual(f"Compressibility   : {args.compressibility:.3e} 1/bar (barostat response of B ~ {b_implied:.1f} GPa)", f_out)
        if args.compressibility == WATER_COMPRESSIBILITY:
            warn("NPT: --compressibility is the 4.57e-5 1/bar placeholder (water) -- for a solid give 1/(B[GPa] x 1e4) "
                 "(Si: ~1e-6): the barostat then responds on its own time scale instead of over- or under-shooting.")
        print_dual("Note              : Berendsen coupling gives the right averages but suppresses fluctuations -- "
                   "do not read heat capacities or compressibilities from them.", f_out)
    print_dual(f"Initial velocities: Maxwell-Boltzmann at {args.temperature:g} K, net momentum removed, seed "
               f"{args.seed if args.seed is not None else 'random'}", f_out)
    total_ps = (args.n_steps + args.equilibration_steps) * args.timestep / 1000.0
    print_dual(f"Length            : {args.equilibration_steps} + {args.n_steps} steps x {args.timestep} fs = "
               f"{args.equilibration_steps * args.timestep / 1000:.3f} + {args.n_steps * args.timestep / 1000:.3f} ps "
               f"({total_ps:.3f} ps in total)", f_out)
    print_dual(f"Sampling          : every {args.stride} steps ({args.stride * args.timestep:g} fs) -> "
               f"{args.n_steps // args.stride} frames", f_out)
    if args.timestep > 2.0:
        warn(f"MD setup: a {args.timestep:g} fs timestep is long for atomistic MD (light atoms vibrate in ~10 fs) -- "
             "check the NVE drift before trusting it.")
    if args.ensemble != "nve" and args.equilibration_steps == 0:
        warn("MD setup: no --equilibration-steps -- the first ~100-200 steps after the Maxwell-Boltzmann start are a "
             "transient and are inside the production averages.")

    # ---- [4]
    if args.equilibration_steps > 0:
        print_section("[4] EQUILIBRATION (discarded)", f_out)
        print_dual(f"Running {args.equilibration_steps} step(s) before production starts -- "
                   "a freshly Maxwell-Boltzmann-seeded structure needs to relax before its "
                   "trajectory/statistics are representative (verified live: ~100-200 steps "
                   "of E_pot relaxation from the initial condition in a 64-atom test case).", f_out)
        t0 = time.monotonic()
        with capture_library_noise(library_noise, "equilibration"):
            dyn.run(args.equilibration_steps)
        print_dual(f"[OK] Equilibration done. E = {atoms.get_potential_energy():.4f} eV, "
                   f"T = {atoms.get_temperature():.1f} K ({time.monotonic() - t0:.1f} s)", f_out)
    else:
        print_section("[4] EQUILIBRATION", f_out)
        print_dual("None requested (--equilibration-steps 0): production starts from the initial velocities.", f_out)

    # ---- [5]
    print_section("[5] PRODUCTION RUN", f_out)
    e0 = atoms.get_potential_energy()
    print_dual(f"Initial state     : E = {e0:.4f} eV, T = {atoms.get_temperature():.1f} K", f_out)

    frames, t_fs, epots, ekins, temps, press, vols, abcs, fmaxs = [], [], [], [], [], [], [], [], []
    n_chunks = args.n_steps // args.stride
    progress_every = max(n_chunks // 10, 1)
    rows = []
    t_prod = time.monotonic()
    pressure_ok = condensed
    from ase import units
    for chunk in range(1, n_chunks + 1):
        with capture_library_noise(library_noise, "production"):
            dyn.run(args.stride)
            epot = atoms.get_potential_energy()
            f = atoms.get_forces()
            p_now = None
            if pressure_ok:
                try:
                    s_ = atoms.get_stress(include_ideal_gas=True, voigt=True)
                    p_now = -float(np.mean(s_[:3])) / units.GPa
                except Exception:
                    pressure_ok = False
        step = chunk * args.stride
        ekin = atoms.get_kinetic_energy()
        temp = atoms.get_temperature()
        t_fs.append(step * args.timestep)
        epots.append(epot)
        ekins.append(ekin)
        temps.append(temp)
        press.append(p_now)
        vols.append(atoms.get_volume())
        abcs.append(atoms.cell.cellpar())
        fmaxs.append(float(np.linalg.norm(f, axis=1).max()))
        frame = atoms.copy()
        if args.out_format == "xyz":
            frame.info['Epot'] = epot
            frame.info['Temp'] = temp
            frame.info['Time'] = step * args.timestep
        frames.append(frame)
        if chunk % progress_every == 0 or chunk == n_chunks:
            rows.append(([f"{step}/{args.n_steps}", f"{step * args.timestep / 1000:.3f}", f"{temp:.1f}",
                          f"{epot / n_atoms:.5f}", fmt(p_now, "{:+.3f}"), f"{vols[-1] / n_atoms:.4f}",
                          f"{fmaxs[-1]:.2f}"], None))
    rest = args.n_steps - n_chunks * args.stride
    if rest:
        with capture_library_noise(library_noise, "production"):
            dyn.run(rest)
    wall = time.monotonic() - t_prod
    if rows:
        print_table(["Step", "t (ps)", "T (K)", "E_pot/atom (eV)", "P (GPa)", "V/atom (Ang^3)", "Fmax (eV/Ang)"],
                    rows, f_out)
    ef = atoms.get_potential_energy()
    print_dual(f"Final state       : E = {ef:.4f} eV, T = {atoms.get_temperature():.1f} K", f_out)
    print_dual(f"[OK] Collected {len(frames)} frame(s) (stride {args.stride}).", f_out)
    ms_step = 1000 * wall / max(args.n_steps, 1)
    ns_day = (86400 / (wall / max(args.n_steps, 1))) * args.timestep * 1e-6 if wall > 0 else None
    print_dual(f"Speed             : {ms_step:.1f} ms/step = {1000 * ms_step / n_atoms:.0f} us/atom/step on {dev_label} "
               f"-> {fmt(ns_day, '{:.3f}')} ns/day at {args.timestep:g} fs ({wall:.0f} s for the production run)", f_out)

    if not frames:
        fail("No frame was collected (--n-steps smaller than --stride).")

    t_fs = np.array(t_fs)
    epots, ekins, temps = np.array(epots), np.array(ekins), np.array(temps)
    etots = epots + ekins
    vols = np.array(vols)
    P = np.array(press, float) if pressure_ok and all(p is not None for p in press) else None
    epa, eka, eta = epots / n_atoms, ekins / n_atoms, etots / n_atoms
    t_ps = t_fs / 1000.0

    # ---- [6]
    print_section("[6] THERMODYNAMICS (production averages)", f_out)
    print_dual(f"Samples           : {len(frames)} frames over {t_ps[-1] - t_ps[0] + args.stride * args.timestep / 1000:.3f} ps; "
               f"+/- = standard error from {N_BLOCKS} block averages", f_out)
    rows = [(["Temperature (K)", pm(temps.mean(), block_se(temps)), f"{temps.std():.2f}",
              f"{temps.min():.1f} .. {temps.max():.1f}"], None),
            (["E_pot (eV/atom)", pm(epa.mean(), block_se(epa), "{:.6f}"), f"{EV_TO_MEV * epa.std():.3f} meV",
              f"{epa.min():.5f} .. {epa.max():.5f}"], None),
            (["E_kin (eV/atom)", pm(eka.mean(), block_se(eka), "{:.6f}"), f"{EV_TO_MEV * eka.std():.3f} meV",
              f"{eka.min():.5f} .. {eka.max():.5f}"], None),
            (["E_tot (eV/atom)", pm(eta.mean(), block_se(eta), "{:.6f}"), f"{EV_TO_MEV * eta.std():.3f} meV",
              f"{eta.min():.5f} .. {eta.max():.5f}"], None)]
    if P is not None:
        rows.append((["Pressure (GPa)", pm(P.mean(), block_se(P), "{:.3f}"), f"{P.std():.3f}",
                      f"{P.min():.2f} .. {P.max():.2f}"], None))
    va = vols / n_atoms
    if args.ensemble == "npt":
        rows.append((["Volume (Ang^3/atom)", pm(va.mean(), block_se(va), "{:.4f}"), f"{va.std():.4f}",
                      f"{va.min():.4f} .. {va.max():.4f}"], None))
        dens = mass / vols * AMU_A3_TO_G_CM3
        rows.append((["Density (g/cm^3)", pm(dens.mean(), block_se(dens), "{:.4f}"), f"{dens.std():.4f}",
                      f"{dens.min():.4f} .. {dens.max():.4f}"], None))
    print_table(["Quantity", "Mean", "Std (fluct.)", "Range"], rows, f_out)
    print_dual(f"Temperature       : mean = {temps.mean():.1f} K, std = {temps.std():.1f} K "
               f"(target {args.temperature:.1f} K)" if args.ensemble != "nve" else
               f"Temperature       : mean = {temps.mean():.1f} K, std = {temps.std():.1f} K "
               f"(started at {args.temperature:.1f} K)", f_out)
    if args.ensemble == "npt":
        abc = np.array(abcs)
        print_dual(f"Cell volume       : mean = {vols.mean():.3f} Ang^3, std = {vols.std():.3f} Ang^3 "
                   f"({100 * (vols.mean() / v0 - 1):+.3f} % vs the input cell)", f_out)
        print_dual(f"Lattice (mean)    : a {abc[:, 0].mean():.4f} +/- {fmt(block_se(abc[:, 0]), '{:.4f}')}, "
                   f"b {abc[:, 1].mean():.4f} +/- {fmt(block_se(abc[:, 1]), '{:.4f}')}, "
                   f"c {abc[:, 2].mean():.4f} +/- {fmt(block_se(abc[:, 2]), '{:.4f}')} Ang", f_out)
        print_dual(f"                    alpha {abc[:, 3].mean():.3f}, beta {abc[:, 4].mean():.3f}, "
                   f"gamma {abc[:, 5].mean():.3f} deg; vs input a {100 * (abc[:, 0].mean() / cp[0] - 1):+.3f} %, "
                   f"b {100 * (abc[:, 1].mean() / cp[1] - 1):+.3f} %, c {100 * (abc[:, 2].mean() / cp[2] - 1):+.3f} %", f_out)
    if args.ensemble == "nvt" and condensed and n_atoms >= 8:
        cv = float(np.var(etots) / (KB_EV * temps.mean()) ** 2 / n_atoms)   # C_V / (N k_B)
        print_dual(f"Heat capacity     : {cv:.2f} kB/atom from the total-energy fluctuations (classical harmonic "
                   "solid: 3, Dulong-Petit; noisy unless the run is long)", f_out)

    # ---- [7]
    print_section("[7] ENSEMBLE CHECKS", f_out)
    sig_t_exp = args.temperature * np.sqrt(2.0 / n_dof) if n_dof > 0 else None
    if args.ensemble == "nve":
        drift = slope_per_ps(t_fs, eta * EV_TO_MEV)
        e_fluct = EV_TO_MEV * eta.std()
        p_fluct = EV_TO_MEV * epa.std()
        print_dual(f"Energy drift      : {drift:+.4f} meV/atom/ps (linear fit of E_tot(t); good: |drift| < "
                   f"{NVE_DRIFT_WARN:g})", f_out)
        print_dual(f"E_tot fluctuation : {e_fluct:.4f} meV/atom (std), {e_fluct / max(p_fluct, 1e-12):.3f} x the E_pot "
                   f"fluctuation ({p_fluct:.3f} meV/atom) -- small means the timestep integrates the motion well", f_out)
        print_dual(f"E_tot start -> end: {etots[0]:.6f} -> {etots[-1]:.6f} eV ({EV_TO_MEV * (eta[-1] - eta[0]):+.4f} meV/atom)", f_out)
        print_dual(f"Temperature       : {temps.mean():.1f} K mean (started at {args.temperature:g} K; from a perfect "
                   "lattice NVE equipartition halves it)", f_out)
        if drift is not None and abs(drift) > NVE_DRIFT_WARN:
            warn(f"NVE: total energy drifts by {drift:+.3f} meV/atom/ps (> {NVE_DRIFT_WARN:g}) -- use a smaller --timestep, "
                 "or the model's forces are not smooth here.")
        if e_fluct > 0.1 * p_fluct and p_fluct > 0:
            warn(f"NVE: E_tot fluctuates by {e_fluct / p_fluct:.2f} x the potential-energy fluctuation -- the timestep is "
                 "too long for these vibrations.")
    else:
        ratio = temps.std() / sig_t_exp if sig_t_exp else None
        print_dual(f"Thermostat        : <T> {temps.mean():.2f} K vs target {args.temperature:g} K "
                   f"({temps.mean() - args.temperature:+.2f} K, block SE {fmt(block_se(temps), '{:.2f} K')})", f_out)
        print_dual(f"T fluctuation     : {temps.std():.2f} K vs canonical T sqrt(2/N_dof) = {fmt(sig_t_exp, '{:.2f}')} K "
                   f"(ratio {fmt(ratio)}" + ("; Berendsen does not sample canonical fluctuations, so this is no check here)" if args.ensemble == "npt" else ")"), f_out)
        se_t = block_se(temps)
        if se_t is not None and abs(temps.mean() - args.temperature) > max(3 * se_t, 0.02 * args.temperature):
            warn(f"Thermostat: <T> = {temps.mean():.1f} K is off the {args.temperature:g} K target by more than 3 standard "
                 "errors -- run longer or equilibrate first.")
        if args.ensemble == "nvt" and ratio is not None and not 0.6 < ratio < 1.6:
            warn(f"Thermostat: the temperature fluctuation is {ratio:.2f} x the canonical value -- the run is too short "
                 "for the Langevin friction, or not equilibrated.")
        if args.ensemble == "npt" and P is not None:
            se_p = block_se(P)
            print_dual(f"Barostat          : <P> {P.mean():+.3f} GPa vs target {args.pressure:g} GPa "
                       f"({P.mean() - args.pressure:+.3f}, block SE {fmt(se_p, '{:.3f} GPa')})", f_out)
            if se_p is not None and abs(P.mean() - args.pressure) > max(3 * se_p, 0.5):
                warn(f"Barostat: <P> = {P.mean():.2f} GPa is off the {args.pressure:g} GPa target -- the cell is still "
                     "relaxing (longer --equilibration-steps) or --compressibility is far from 1/B.")
    for name, series, unit in [("Temperature", temps, "K"), ("E_pot", epa * EV_TO_MEV, "meV/atom")] + \
            ([("Volume", va, "Ang^3/atom")] if args.ensemble == "npt" else []):
        shift, se = thirds_shift(series)
        if shift is None:
            continue
        flag = se is not None and abs(shift) > 3 * se and (name != "E_pot" or args.ensemble != "nve")
        print_dual(f"Stationarity      : {name} last third - first third = {shift:+.4f} {unit} "
                   f"(combined SE {fmt(se, '{:.4f}')})" + (" -- DRIFTING" if flag else ""), f_out)
        if flag:
            warn(f"Stationarity: {name} still drifts between the first and last third of the run ({shift:+.3g} {unit}, "
                 f"> 3 standard errors) -- longer --equilibration-steps.")

    # ---- [8]
    print_section("[8] STRUCTURE ALONG THE RUN & MODEL DOMAIN", f_out)
    t0 = time.monotonic()
    d_mins = np.array([model_card.shortest_distance(a, cutoff=min(6.0, max(3.0, 1.5 * (d_min0 or 2.0)))) or np.nan
                       for a in frames])
    i_dmin = int(np.nanargmin(d_mins)) if np.any(np.isfinite(d_mins)) else None
    if i_dmin is not None:
        print_dual(f"Shortest distance : min {d_mins[i_dmin]:.4f} Ang at t = {t_ps[i_dmin]:.3f} ps, mean "
                   f"{np.nanmean(d_mins):.4f} Ang (start {fmt(d_min0, '{:.4f}')})", f_out)
        if d_min0 and d_mins[i_dmin] < 0.6 * d_min0:
            warn(f"Structure: two atoms came to {d_mins[i_dmin]:.2f} Ang (start {d_min0:.2f}) -- the run may have blown up.")
    rdf = compute_rdf(frames) if condensed else None
    if rdf and rdf["r_peak"]:
        print_dual(f"g(r)              : first peak {rdf['r_peak']:.3f} Ang (g = {rdf['g_peak']:.2f}), first minimum "
                   f"{fmt(rdf['r_min'], '{:.3f}')} Ang, coordination {fmt(rdf['coordination'], '{:.2f}')} "
                   f"({rdf['n_frames']} frames, r <= {rdf['rmax']:.2f} Ang)", f_out)
    elif not condensed:
        print_dual("g(r)              : not computed for a vacuum-padded cell (the density normalization needs a bulk cell)", f_out)
    msd = compute_msd([a.positions for a in frames], [np.asarray(a.cell) for a in frames],
                      atoms.get_chemical_symbols(), args.stride * args.timestep)
    lind = None
    if msd:
        d_nn = rdf["r_peak"] if rdf and rdf["r_peak"] else d_nn0
        lind = msd["rms_about_mean"] / d_nn if d_nn and condensed and n_atoms >= 8 else None
        for key in (["all"] + (msd["elements"] if len(msd["elements"]) > 1 else [])):
            print_dual(f"MSD {key:<14}: {msd['msd'][key][-1]:.4f} Ang^2 at lag {msd['t_fs'][-1] / 1000:.3f} ps, "
                       + (f"D = {msd['D_cm2_s'][key]:.3e} cm^2/s (fit over the second half of the lags)"
                          if msd['D_cm2_s'][key] is not None else "D not fitted (too few frames: lower --stride)"), f_out)
        print_dual(f"Displacement      : RMS about the mean positions {msd['rms_about_mean']:.4f} Ang; start->end RMS "
                   f"{msd['rms_total']:.4f}, max {msd['max_total']:.4f} Ang", f_out)
        if lind is not None:
            state_ = ("solid-like: atoms vibrate about fixed sites" if lind < 0.10 else
                      "large vibrations or diffusion: near melting, a liquid, or a defect hopping" if lind < 0.15 else
                      "liquid-like / diffusive")
            print_dual(f"Lindemann ratio   : {lind:.3f} (RMS about the mean / nearest-neighbour distance {d_nn:.3f} Ang) -- "
                       f"{state_}", f_out)
        if lind is not None and lind < 0.10:
            print_dual("                    (in a solid the MSD levels off at the vibrational amplitude: the D above "
                       "is a fit to that plateau, not diffusion)", f_out)
        d_all = msd["D_cm2_s"]["all"]
        if d_all is not None and d_all > 1e-5 and (lind or 0) > 0.15:
            print_dual("                    D > 1e-5 cm^2/s with a large Lindemann ratio: the system flows (a liquid's D "
                       "is ~1e-5 to 1e-4 cm^2/s)", f_out)
    print_dual(f"(structure analysis: {time.monotonic() - t0:.1f} s)", f_out)

    # trained domain of the model along the run
    if card is not None:
        dom = card.get("domain") or {}
        rows = []

        def status(ok):
            return "inside" if ok else "OUTSIDE"
        if dom.get("measure") and dom.get("measure_name") and condensed:
            lo, hi = dom["measure"]["min"], dom["measure"]["max"]
            run_lo, run_hi = va.min(), va.max()
            ok = lo * (1 - model_card.MEASURE_TOL) <= run_lo and run_hi <= hi * (1 + model_card.MEASURE_TOL)
            rows.append(([dom["measure_name"], f"{lo:.3f} .. {hi:.3f}", f"{run_lo:.3f} .. {run_hi:.3f}", status(ok)],
                         None if ok else 'yellow'))
            if not ok:
                warn(f"Model domain: the {dom['measure_name']} reached {run_lo:.3f}-{run_hi:.3f}, outside the trained "
                     f"{lo:.3f}-{hi:.3f} -- the model is extrapolating in cell size.")
        dd = dom.get("d_min_ang")
        if dd and i_dmin is not None:
            ok = d_mins[i_dmin] >= model_card.D_MIN_FACTOR * dd["min"]
            n_bad = int(np.sum(d_mins < model_card.D_MIN_FACTOR * dd["min"]))
            rows.append((["shortest distance (Ang)", f">= {dd['min']:.3f}", f"min {d_mins[i_dmin]:.3f}",
                          status(ok) + (f" ({n_bad} frames)" if n_bad else "")], None if ok else 'yellow'))
            if not ok:
                warn(f"Model domain: in {n_bad} of {len(frames)} frames two atoms are closer than any pair in the training "
                     f"data ({d_mins[i_dmin]:.3f} < {dd['min']:.3f} Ang).")
        fm = dom.get("force_max_eV_per_ang")
        if fm:
            fmx = max(fmaxs)
            ok = fmx <= model_card.FORCE_FACTOR * fm["max"]
            rows.append((["max force (eV/Ang)", f"<= {fm['max']:.2f}", f"max {fmx:.2f}", status(ok)], None if ok else 'yellow'))
            if not ok:
                warn(f"Model domain: predicted forces reach {fmx:.1f} eV/Ang, above {model_card.FORCE_FACTOR:g} x the largest "
                     f"training force ({fm['max']:.1f}).")
        er = dom.get("energy_per_atom_eV")
        if er:
            ok = epa.max() <= er["max"] + 0.5 * (er["max"] - er["min"])
            rows.append((["E_pot per atom (eV)", f"{er['min']:.4f} .. {er['max']:.4f}",
                          f"{epa.min():.4f} .. {epa.max():.4f}", status(ok)], None if ok else 'yellow'))
            if not ok:
                warn("Model domain: the potential energy per atom rose well above anything in the training data -- the run "
                     "visits configurations the model never saw (e.g. a melt from a crystal-only training set).")
        pr = dom.get("pressure_GPa")
        if pr and P is not None:
            ok = pr["min"] - 2 <= P.mean() <= pr["max"] + 2
            rows.append((["pressure (GPa)", f"{pr['min']:.2f} .. {pr['max']:.2f}", f"mean {P.mean():+.2f}", status(ok)],
                         None if ok else 'yellow'))
        na = dom.get("n_atoms")
        if na:
            rows.append((["atoms per cell", f"{na['min']} .. {na['max']}", f"{n_atoms}",
                          "larger (fine: the model is local, r_max "
                          f"{fmt(_card_value(card, 'model', 'architecture', 'r_max'), '{:g}')} Ang)" if n_atoms > na["max"] else "inside"],
                         None))
        print_dual("\n" + color_text("Model domain along the run (model card's trained range):", 'bold'), f_out)
        print_table(["Quantity", "Trained", "This run", "Status"], rows, f_out)
    model_card.collect_warnings(None)
    for m in card_msgs:
        warn(f"Model card check: {m.split(': ', 1)[-1]}")

    # ---- outputs
    title = (f"stb-mlmd {args.ensemble.upper()} {args.temperature:g} K"
             + (f", {args.pressure:g} GPa" if args.ensemble == "npt" else "")
             + f" -- {atoms.get_chemical_formula(mode='hill')} ({n_atoms} atoms), "
             + (os.path.basename(args.custom_model) if args.custom_model else f"MACE-MP-0 {args.model}"))
    ctx = {"title": title, "ensemble": args.ensemble, "T_target": args.temperature, "P_target": args.pressure,
           "n_atoms": n_atoms, "t_fs": t_fs, "t_ps": t_ps, "T": temps, "Epot_atom": epa, "Ekin_atom": eka,
           "Etot_atom": eta, "P": P, "V_atom": va, "abc": [a[:3] for a in abcs], "d_min": d_mins,
           "Fmax": np.array(fmaxs), "rdf": rdf, "msd": msd,
           "d_min_trained": ((card or {}).get("domain") or {}).get("d_min_ang", {}).get("min") if card else None}
    data_paths, gplot_path = [], None
    if args.save_data or args.save_gnuplot:
        data_paths = write_data_files(stem, ctx)
    if args.save_gnuplot:
        gplot_path = write_gnuplot(stem, ctx)
    show = args.show
    if show and not _show_possible():
        warn("Plots: --show needs a display (local desktop or `ssh -X`) -- the figure was not shown"
             + ("; it is saved with --png." if not png_path else "."))
        show = False
    plotted = False
    if png_path or show:
        try:
            with capture_library_noise(library_noise, "plots"):
                make_dashboard(ctx, png_path, show)
            plotted = bool(png_path)
        except Exception as e:
            warn(f"Plots: dashboard not made ({type(e).__name__}: {e}).")

    from ase.io import write as ase_write
    try:
        ase_write(traj_path, frames, format=ase_format)
    except Exception as e:
        fail(f"Failed to write trajectory: {e}")
    new_pmg = AseAtomsAdaptor.get_structure(atoms)
    new_structure = structure_io.from_pymatgen(new_pmg, species_meta=structure.species_meta)
    structure_io.write_fdf(new_structure, args.output_structure)

    # ---- [9]
    print_section("[9] SUMMARY & FILES", f_out)
    print_dual(f"Status            : OK ({len(warns)} warning(s), see [10])" if warns else "Status            : OK", f_out)
    key = [f"<T> {temps.mean():.1f} K"]
    if args.ensemble == "nve":
        key.append(f"drift {slope_per_ps(t_fs, eta * EV_TO_MEV):+.4f} meV/atom/ps")
    if P is not None:
        key.append(f"<P> {P.mean():+.3f} GPa")
    if args.ensemble == "npt":
        key.append(f"<V> {va.mean():.4f} Ang^3/atom ({100 * (vols.mean() / v0 - 1):+.3f} %)")
    if lind is not None:
        key.append(f"Lindemann {lind:.3f}")
    print_dual("Results           : " + ", ".join(key), f_out)
    print_dual(f"Time              : {time.monotonic() - t_start:.0f} s in total, {ms_step:.1f} ms/step", f_out)
    print_dual(f"Trajectory        : {traj_path} ({len(frames)} frame(s), {args.out_format} format)", f_out)
    print_dual(f"Final structure   : {args.output_structure}", f_out)
    if plotted:
        print_dual(f"Dashboard         : {png_path}", f_out)
    for p in data_paths:
        print_dual(f"Data              : {p}", f_out)
    if gplot_path:
        print_dual(f"Gnuplot           : {gplot_path}  (gnuplot {gplot_path} -> {stem}_md.pdf)", f_out)
    if report_path:
        print_dual(f"Report            : {report_path}", f_out)
    extras = [x for x, on in [("--png (dashboard)", not png_path), ("--show", not args.show),
                              ("--save-gnuplot (gnuplot panels)", not args.save_gnuplot),
                              ("--save-data (.dat files)", not (args.save_data or args.save_gnuplot)),
                              ("--save-report", not report_path)] if on]
    if extras:
        print_dual("Also available    : " + ", ".join(extras), f_out)

    # ---- [10]
    print_section("[10] WARNINGS", f_out)
    if warns:
        for i, w in enumerate(warns, 1):
            print_dual(color_text(f"{i:2d}. {w}", 'yellow'), f_out)
    else:
        print_dual(color_text("None.", 'green'), f_out)
    unique_noise = list(dict.fromkeys(n.split("] ", 1)[-1] for n in library_noise))
    if unique_noise:
        print_dual(f"Library messages  : {len(unique_noise)} distinct informational message(s) from mace/torch/ase captured "
                   "(model loading, dtype, deprecations) -- not part of the results.", f_out)

    if f_out:
        f_out.close()

    print("\n" + "-" * 60)
    print(color_text("MD run complete -- a fast heuristic, not a substitute for ab initio MD.\n", 'bold'))


if __name__ == "__main__":
    main()
