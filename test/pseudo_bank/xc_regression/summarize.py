#!/usr/bin/env python
"""Summary of a run of the bank regression: python summarize.py CASES

A case passes when SIESTA wrote 0_NORMAL_EXIT. Failures are sorted by what calc.out/calc.err say:
polarization (perturbative polarization orbital), ranks (more ranks than orbitals/atoms), ekb0
('IMAX = 0'), scalar (no scalar-relativistic part), basis-charge (sr3plus without synthetic-atoms),
nan, timeout (exit code 124), other. Writes CASES/results.tsv; exit code 1 when any case failed.
"""
import os
import re
import sys
from collections import Counter, defaultdict

CATEGORIES = [
    ("polarization", r"Generation of polarization orbital failed"),
    ("ranks", r"too many processors|more processors than|Number of nodes.*larger"),
    ("ekb0", r"IMAX = 0"),
    ("scalar", r"scalar-relativistic part is not available"),
    ("basis-charge", r"Total charge in occupied basis states|nsemx too small"),
    ("nan", r"\bNaN\b"),
]


def classify(run):
    if os.path.exists(os.path.join(run, "0_NORMAL_EXIT")):
        return "OK", ""
    text = ""
    for name in ("calc.out", "calc.err"):
        try:
            text += open(os.path.join(run, name), errors="replace").read()
        except OSError:
            pass
    try:
        if open(os.path.join(run, "exit_code")).read().strip() == "124":
            return "timeout", ""
    except OSError:
        pass
    for cat, pattern in CATEGORIES:
        if re.search(pattern, text, re.I):
            return cat, ""
    line = next((ln.strip() for ln in text.splitlines() if re.search(r"error|stop|abort", ln, re.I)), "")
    return "other", line[:120]


def main():
    root = sys.argv[1] if len(sys.argv) > 1 else "cases"
    rows = []
    for dirpath, dirs, files in os.walk(root):
        base = os.path.basename(dirpath)
        if re.fullmatch(r"np\d+", base):
            case = os.path.dirname(dirpath)
            bank = os.path.basename(os.path.dirname(case))
            status, detail = classify(dirpath)
            rows.append((bank, os.path.basename(case), base[2:], status, detail))
    rows.sort(key=lambda r: (r[0], int(r[2]), r[1]))
    with open(os.path.join(root, "results.tsv"), "w") as f:
        f.write("bank\tcase\tranks\tstatus\tdetail\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    by = defaultdict(Counter)
    for bank, _case, ranks, status, _d in rows:
        by[(bank, ranks)][status] += 1
    print(f"{'bank':<28} {'ranks':>5} {'cases':>6} {'OK':>5}  failures")
    for (bank, ranks), cnt in sorted(by.items(), key=lambda kv: (kv[0][0], int(kv[0][1]))):
        fails = ", ".join(f"{k} {v}" for k, v in sorted(cnt.items()) if k != "OK")
        print(f"{bank:<28} {ranks:>5} {sum(cnt.values()):>6} {cnt['OK']:>5}  {fails or '-'}")
    failed = [r for r in rows if r[3] != "OK"]
    if failed:
        print("\nFailed cases:")
        for bank, case, ranks, status, detail in failed:
            print(f"  {bank}/{case} np{ranks}: {status} {detail}")
    print(f"\n{len(rows) - len(failed)}/{len(rows)} OK")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
