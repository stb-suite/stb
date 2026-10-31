#!/bin/bash
# =====================================================================================================
#  make_reference.sh -- regenerates reference_results/<case>/ of the 3.5 example (stb-ani2traj)
#
#  Every case runs on a copy of siesta_runs/ (three real, finished SIESTA 5.4 AIMD runs: 8 water molecules NVT,
#  8-atom Si NVE, 8-atom Si NPT) in a scratch folder; only the reports and the first lines of the written
#  trajectories are kept (scratch paths written as <work>, dates as <date>). Seconds, no SIESTA.
#
#      source ~/.anaconda3/etc/profile.d/conda.sh && conda activate stb    # any environment with stb_suite
#      bash make_reference.sh
# =====================================================================================================
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
REF="$HERE/reference_results"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
clean() { sed 's/\x1b\[[0-9;]*m//g' | sed "s#$WORK#<work>#g" | sed -E 's/^(Date\/time *: ).*/\1<date>/'; }
report() {       # report <case> <run folder> <stb-ani2traj args...>: a fresh copy of the run, the plain-text report
    local case="$1" run="$2"; shift 2
    rm -rf "$WORK/$case"; mkdir -p "$REF/$case"; cp -a "$HERE/siesta_runs/$run" "$WORK/$case"
    echo "  [run]  $case: stb-ani2traj $*"
    (cd "$WORK/$case" && stb-ani2traj --no-intro "$@" 2>&1) | clean > "$REF/$case/report.txt" || true
}

echo "== 1. Water, NVT: the whole report (xsf)"
report water water_nvt -l water
head -n 14 "$WORK/water/water_traj.xsf" > "$REF/water/xsf_head.txt"

echo "== 2. Si, NPT: a variable cell (extended xyz with per-frame thermodynamics)"
report si_npt si_npt -l si -of xyz
head -n 2 "$WORK/si_npt/si_traj.xyz" > "$REF/si_npt/xyz_head.txt"

echo "== 3. Si, NVE: energy conservation"
report si_nve si_nve -l si -of xyz

echo "== 4. --unwrap on positions SIESTA never wrapped"
report si_nve_unwrap si_nve -l si -of xyz --unwrap

echo "== 5. Equilibration and stride, pdb"
report water_pdb water_nvt -l water -of pdb --skip 100 --stride 5
head -n 6 "$WORK/water_pdb/water_traj.pdb" > "$REF/water_pdb/pdb_head.txt"

echo "== 6. Guard rails"
report guard_nolog si_npt -l si --dry-run
rm -rf "$WORK/guard_nolog"; cp -a "$HERE/siesta_runs/si_npt" "$WORK/guard_nolog"; rm "$WORK/guard_nolog/calc.out"
(cd "$WORK/guard_nolog" && stb-ani2traj --no-intro -l si --dry-run 2>&1) | clean > "$REF/guard_nolog/report.txt" || true
rm -rf "$WORK/guard_appended"; mkdir -p "$REF/guard_appended"; cp -a "$HERE/siesta_runs/si_nve" "$WORK/guard_appended"
cat "$WORK/guard_appended/si.ANI" "$WORK/guard_appended/si.ANI" > "$WORK/guard_appended/x" && mv "$WORK/guard_appended/x" "$WORK/guard_appended/si.ANI"
(cd "$WORK/guard_appended" && stb-ani2traj --no-intro -l si --dry-run 2>&1) | clean > "$REF/guard_appended/report.txt" || true
report guard_missing water_nvt -l nothing

echo "== 7. The interactive menu (stb-suite, code 3.5)"
rm -rf "$WORK/menu"; cp -a "$HERE/siesta_runs/si_npt" "$WORK/menu"; mkdir -p "$REF/menu"
# label (Enter: the only *.ANI), format 3 = xyz, skip 0, stride 1, no unwrap, default output, no dry run, no report
(cd "$WORK/menu" && printf '3.5\n\n3\n0\n1\nn\n\nn\nn\n\n0\n' | stb-suite 2>&1) | clean \
    | sed -n '/STB-ANI2TRAJ REPORT/,/^\[9\] WARNINGS/p' > "$REF/menu/report.txt"

echo "done -- reference_results/<case>/"
