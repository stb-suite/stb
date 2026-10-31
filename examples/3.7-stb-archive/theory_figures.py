"""Figure of the 3.7 example (stb-archive), from stb-archive's own selection (plan_archive) on runs/si_relax:

  figures/contents.png   the folder by kind of file: what goes into the archive and what is left out, for each mode

    python3 theory_figures.py        # seconds
"""
import os
import shutil
import tempfile
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from stb.archive import plan_archive
from stb.core.pseudopotentials import bank_path

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)
work = tempfile.mkdtemp()
try:
    run = os.path.join(work, "si_relax")
    shutil.copytree(os.path.join(HERE, "runs", "si_relax"), run)
    shutil.copy(os.path.join(bank_path("dojo"), "Si.psml"), run)
    plans = {m: plan_archive(run, m, True, [], []) for m in ("inputs", "results", "full")}
finally:
    shutil.rmtree(work)

kinds = {}
for rel, kind, size, why in plans["full"]["keep"]:
    kinds[kind] = kinds.get(kind, 0) + size
order = sorted(kinds, key=lambda k: kinds[k])
fig, ax = plt.subplots(figsize=(9.5, 4.4))
y = np.arange(len(order))
h = 0.26
for k, (mode, color) in enumerate((("inputs", "#2e7d32"), ("results", "#1f6f8b"), ("full", "#999999"))):
    inside = {}
    for rel, kind, size, why in plans[mode]["keep"]:
        inside[kind] = inside.get(kind, 0) + size
    vals = [inside.get(kd, 0) / 1024 for kd in order]
    total = sum(r[2] for r in plans[mode]["keep"]) / 1024
    ax.barh(y + (1 - k) * h, vals, height=h, color=color, label=f"--mode {mode}: {len(plans[mode]['keep'])} files, "
                                                                f"{total:.0f} KB")
ax.set_yticks(y)
ax.set_yticklabels(order, fontsize=8)
ax.set_xlabel("size in the archive (KB, before compression)")
ax.set_title("runs/si_relax (+ analysis/): what each mode puts into the archive, by kind of file")
ax.legend(frameon=False, fontsize=8, loc="lower right")
ax.grid(axis="x", alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "contents.png"), dpi=150)
print({m: (len(p["keep"]), sum(r[2] for r in p["keep"])) for m, p in plans.items()})
print("figures/contents.png")
