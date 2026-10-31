#!/bin/bash
# =====================================================================================================
#  make_reference.sh -- regenerates reference_results/<case>/ of the 3.3 example (stb-cube)
#
#  Every case runs on a copy of siesta_runs/ (three real, finished SIESTA 5.4 runs: bulk Si, a water molecule,
#  spin-polarized O2) in a scratch folder; only the reports and the first lines of each cube are kept (scratch
#  paths written as <work>). Seconds, no SIESTA. The SIESTA runs themselves come from siesta_runs/run_all.sh.
#
#      source ~/.anaconda3/etc/profile.d/conda.sh && conda activate stb    # any environment with stb_suite
#      bash make_reference.sh
# =====================================================================================================
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
REF="$HERE/reference_results"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
cp -a "$HERE/siesta_runs/si" "$HERE/siesta_runs/h2o" "$HERE/siesta_runs/o2" "$WORK/"

report() {       # report <case> <folder> <stb-cube args...>: runs in $WORK/<folder>, stores the plain-text report
    local case="$1" dir="$2"; shift 2
    mkdir -p "$REF/$case"
    echo "  [run]  $case: stb-cube $*"
    (cd "$WORK/$dir" && stb-cube --no-intro "$@" 2>&1) | sed 's/\x1b\[[0-9;]*m//g' | sed "s#$WORK#<work>#g" \
        | sed -E 's/^(Date\/time *: ).*/\1<date>/; s/(written )[0-9-]+ [0-9:]+\)/\1<date>)/' > "$REF/$case/report.txt" || true
}
header() {       # header <case> <cube>: the cube's header (comment lines, grid, atoms) and its first data line
    head -n "$(( $(awk 'NR==3 {print ($1 < 0 ? -$1 : $1) + 7}' "$WORK/$2") ))" "$WORK/$2" > "$REF/$1/cube_header.txt"
}

echo "== 1. Bulk Si: the valence density, the whole report"
report si_rho si -l si -t RHO
header si_rho si/si_RHO.cube

echo "== 2. Units: the same grid written the naive way (sisl's own writer, e/Ang^3 next to Bohr)"
mkdir -p "$REF/units"
(cd "$WORK/si" && python3 - <<'PY'
import sisl, numpy as np
from ase.io.cube import read_cube_data
g = sisl.get_sile("si.RHO").read_grid()
g.geometry = sisl.get_sile("si.XV").read_geometry()
g.write("naive.cube")                                     # what a plain sisl conversion writes
for name in ("naive.cube", "si_RHO.cube"):
    data, atoms = read_cube_data(name)
    dv_bohr = atoms.get_volume() / data.size / 0.529177210903**3
    print(f"{name:12s}  max {data.max():8.5f}   sum x voxel (Bohr^3) = {data.sum() * dv_bohr:8.4f} e")
PY
) > "$REF/units/integrals.txt"

echo "== 3. Every grid of the Si run"
for t in RHO DRHO BADER TOCH IOCH VT VH VNA; do report si_all/$t si -l si -t "$t" --dry-run; done

echo "== 4. A molecule in a box: water"
report h2o_rho  h2o -l h2o -t RHO
report h2o_drho h2o -l h2o -t DRHO
report h2o_vh   h2o -l h2o -t VH --potential-unit Ha
header h2o_vh h2o/h2o_VH.cube

echo "== 5. Spin: O2 triplet"
for s in total up down diff; do report o2_$s o2 -l o2 -t RHO --spin "$s"; done

echo "== 6. Guard rails"
report guard_missing si -l si -t LDOS
report guard_spin si -l si -t RHO --spin diff --dry-run
cp -a "$WORK/si" "$WORK/si_noxv"; rm "$WORK/si_noxv/si.XV"
report guard_noxv si_noxv -l si -t RHO --dry-run
cp -a "$WORK/si" "$WORK/si_stale"; cp "$WORK/o2/o2.XV" "$WORK/si_stale/si.XV"
report guard_stale si_stale -l si -t RHO --dry-run
cp -a "$WORK/si" "$WORK/si_two"; cp "$WORK/si_two/si.RHO" "$WORK/si_two/si_old.RHO"
report guard_ambiguous si_two -t RHO

echo "== 7. From the parent folder, label found automatically, report saved"
report label_path . -l si/si -t VT --potential-unit Ry --save-report
report autolabel h2o -t DRHO --dry-run

echo "== 8. The interactive menu (stb-suite, code 3.3)"
mkdir -p "$REF/menu"
# type 10 = VH, label h2o, spin total, unit Ha, default output, no dry run, no report, Enter, 0 = exit
(cd "$WORK/h2o" && printf '3.3\n10\nh2o\n1\nHa\n\nn\nn\n\n0\n' | stb-suite 2>&1) | sed 's/\x1b\[[0-9;]*m//g' \
    | sed "s#$WORK#<work>#g" | sed -E 's/^(Date\/time *: ).*/\1<date>/; s/(written )[0-9-]+ [0-9:]+\)/\1<date>)/' \
    | sed -n '/STB-CUBE REPORT/,/hip to be square/p' > "$REF/menu/report.txt"

echo "done -- reference_results/<case>/"
