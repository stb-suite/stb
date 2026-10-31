#!/usr/bin/env python3
"""Figures of the 3.2 README (figures/), from this folder's own data:

    run_contents.png   what the finished SIESTA run in siesta_run/ is made of: size and number of files per kind,
                       as stb-clean classifies them (+ the inputs it keeps)
    rerun.png          the experiment of siesta_rerun.sh: SCF iterations per geometry step and wall time when the
                       cleaned folder is run again, starting fresh vs keeping .DM and .XV

    python3 theory_figures.py      (stb_suite, matplotlib; seconds)
"""
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from stb import clean
from stb.core.pseudopotentials import bank_path

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
KEEP_C, REMOVE_C = "#2a78d6", "#eb6834"          # categorical slots 1 and 2 (validated default pair)
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e6e5e1"


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelcolor=INK)
    ax.grid(axis="x", color=GRID, lw=0.8)
    ax.set_axisbelow(True)


def fig_contents():
    plan = clean.plan_folder(os.path.join(HERE, "siesta_run"), set(clean.DEFAULT_KEEP), False)
    rows = []
    kept = list(plan["keep"])
    psml = clean.file_size(os.path.join(bank_path("dojo"), "Si.psml"))   # copied from the bundled bank, not shipped
    rows.append(("kept: inputs (fdf, Si.psml, scripts, bandlines.txt)", sum(map(clean.file_size, kept)) + psml,
                 len(kept) + 1, KEEP_C))
    by = {}
    for f in plan["remove"]:
        cat = clean.classify(os.path.basename(f))
        if cat in (clean.CAT_OTHER, clean.CAT_POST) and clean.looks_like_job_script(f):
            cat = clean.CAT_SCRIPT                      # as the report sorts it
        by.setdefault(cat, []).append(f)
    for cat in reversed(clean.CATEGORY_ORDER):
        if cat in by:
            rows.insert(0, (cat, sum(map(clean.file_size, by[cat])), len(by[cat]), REMOVE_C))
    rows = rows[::-1]
    fig, ax = plt.subplots(figsize=(9, 4.6))
    y = range(len(rows))
    ax.barh(list(y), [r[1] / 1024 for r in rows], color=[r[3] for r in rows], height=0.62, edgecolor="white", lw=2)
    for i, r in zip(y, rows):
        ax.text(r[1] / 1024 + 25, i, f"{clean.human_size(r[1])}  ({r[2]} file{'s' if r[2] > 1 else ''})",
                va="center", fontsize=9, color=INK)
    ax.set_yticks(list(y))
    ax.set_yticklabels([r[0] for r in rows], fontsize=9.5)
    ax.invert_yaxis()
    ax.set_xlabel("size (KB)", color=INK)
    ax.set_xlim(0, max(r[1] for r in rows) / 1024 * 1.35)
    style(ax)
    n_all = sum(r[2] for r in rows)
    size_all = clean.human_size(sum(r[1] for r in rows))
    ax.set_title(f"A finished SIESTA run (Si: CG relaxation, bands, PDOS): {n_all} files, {size_all}\n"
                 "blue: kept (default --keep)   orange: removed", fontsize=10.5, color=INK,
                 loc="left")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "run_contents.png"), dpi=150)
    plt.close(fig)


def scf_per_step(out):
    t = open(out).read()
    blocks = re.split(r"Begin CG opt\. move", t)[1:]
    return [len(re.findall(r"^\s+scf:\s+\d+", b, re.M)) for b in blocks]


def wall(clock):
    return float(open(clock).read().split("End of run")[-1].split()[0])


def fig_rerun():
    ref = os.path.join(HERE, "reference_results")
    fresh = scf_per_step(os.path.join(ref, "rerun_fresh", "calc.out"))
    cont = scf_per_step(os.path.join(ref, "rerun_continue", "calc.out"))
    w_f, w_c = wall(os.path.join(ref, "rerun_fresh", "CLOCK")), wall(os.path.join(ref, "rerun_continue", "CLOCK"))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 3.9), gridspec_kw={"width_ratios": [3, 1.1]})
    x = range(1, len(fresh) + 1)
    a1.bar([i - 0.2 for i in x], fresh, width=0.38, color=REMOVE_C, edgecolor="white", lw=1,
           label=f"default clean -- starts fresh: {len(fresh)} steps, {sum(fresh)} SCF iterations")
    a1.bar([i + 0.2 for i in range(1, len(cont) + 1)], cont, width=0.38, color=KEEP_C, edgecolor="white", lw=1,
           label=f"--keep ... .DM .XV -- continues: {len(cont)} step, {sum(cont)} SCF iterations")
    a1.set_xlabel("geometry step (CG move)", color=INK)
    a1.set_ylabel("SCF iterations", color=INK)
    a1.set_xticks(list(x))
    a1.legend(fontsize=8.5, frameon=False, loc="upper right")
    a1.set_ylim(0, max(fresh) * 1.45)
    style(a1)
    a1.grid(axis="x", visible=False)
    a1.grid(axis="y", color=GRID, lw=0.8)
    a1.set_title("(a) the rerun, step by step", fontsize=10.5, color=INK, loc="left")
    a2.bar([0, 1], [w_f, w_c], color=[REMOVE_C, KEEP_C], width=0.6, edgecolor="white", lw=2)
    for i, w in enumerate((w_f, w_c)):
        a2.text(i, w + 0.4, f"{w:.1f} s", ha="center", fontsize=9.5, color=INK)
    a2.set_xticks([0, 1])
    a2.set_xticklabels(["fresh", "continue"], color=INK)
    a2.set_ylabel("wall time (s, 4 MPI)", color=INK)
    a2.set_ylim(0, w_f * 1.2)
    style(a2)
    a2.grid(axis="x", visible=False)
    a2.grid(axis="y", color=GRID, lw=0.8)
    a2.set_title("(b) cost", fontsize=10.5, color=INK, loc="left")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "rerun.png"), dpi=150)
    plt.close(fig)
    print("fresh", fresh, sum(fresh), w_f, "| continue", cont, sum(cont), w_c)


if __name__ == "__main__":
    os.makedirs(FIG, exist_ok=True)
    fig_contents()
    fig_rerun()
    print("figures written to", FIG)
