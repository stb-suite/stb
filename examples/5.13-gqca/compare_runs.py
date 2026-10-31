#!/usr/bin/env python3
"""Side-by-side summary of several finished GQCA Stage 2 runs (helper of
example_5.13.sh). Usage: compare_runs.py <output_dir> <run> [<run> ...]

Each <run> is a folder <output_dir>/<run>/ on which stb-gqcaAnalysis already
ran (so gqca_results_clusters.csv exists) plus its report text
<output_dir>/<run>_stage2.log. Prints, per run: the basis read back from the
folders, the relaxed end-member lattice constants (zincblende conventional
cell: a = (4 * V/site)^(1/3)), their mismatch, Delta_1, and the
regular-solution equivalents Omega = 2*Delta_1 and T_c = Delta_1/kB.
"""
import csv
import re
import sys

KB = 8.617333262e-5  # eV/K

out, runs = sys.argv[1], sys.argv[2:]
cols = []
for run in runs:
    rows = {r["cluster"]: r for r in csv.DictReader(open(f"{out}/{run}/gqca_results_clusters.csv"))}
    txt = open(f"{out}/{run}_stage2.log").read()
    basis = re.search(r"^Basis size\s*\|\s*(\S+)", txt, re.M).group(1)
    mesh = re.search(r"^Mesh cutoff\s*\|\s*(\d+)", txt, re.M).group(1)
    pair_b = re.search(r"2 \| BB \((\S+?)-", txt).group(1)
    a_aa = (4 * float(rows["cluster_n0"]["v_per_site_A3"])) ** (1 / 3)
    a_bb = (4 * float(rows["cluster_n2"]["v_per_site_A3"])) ** (1 / 3)
    d1 = float(rows["cluster_n1"]["delta_eV"])
    cols.append((run, pair_b, a_aa, a_bb, 100 * (a_bb / a_aa - 1), 1000 * d1,
                 2000 * d1, d1 / KB, basis, mesh))

labels = ["run", "basis / mesh cutoff (Ry)", "a(GaAs) relaxed (Ang)", "a(B-As) relaxed (Ang)", "mismatch (%)",
          "Delta_1 (meV/site)", "Omega = 2*Delta_1 (meV)", "T_c = Delta_1/kB (K)"]
fmt = [lambda c: c[0], lambda c: f"{c[8]} / {c[9]}", lambda c: f"{c[2]:.3f}", lambda c: f"{c[3]:.3f} ({c[1]})",
       lambda c: f"{c[4]:+.2f}", lambda c: f"{c[5]:.1f}", lambda c: f"{c[6]:.1f}", lambda c: f"{c[7]:.0f}"]
w = max(len(lab) for lab in labels) + 2
print(" " * w + "".join(f"{fmt[0](c):>16s}" for c in cols))
for lab, f in zip(labels[1:], fmt[1:]):
    print(f"{lab:<{w}s}" + "".join(f"{f(c):>16s}" for c in cols))
print(f"{'(experiment, 300 K)':<{w}s}" + "  GaAs 5.653, InAs 6.058, AlAs 5.661 Ang"
      "  (Vurgaftman et al., J. Appl. Phys. 89, 5815, 2001)")
