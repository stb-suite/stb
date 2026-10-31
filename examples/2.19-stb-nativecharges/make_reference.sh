#!/bin/bash
# =====================================================================================================
#  make_reference.sh -- regenerates reference_results/<case>/report.txt of the 2.19 example (stb-nativecharges)
#
#  Every case runs on a copy of runs/ (real SIESTA 5.4 runs with Mulliken, Hirshfeld and Voronoi charges in the log;
#  run_example.sh) in a scratch folder. Dates become <date>. Seconds, no SIESTA.
#
#      source ~/.anaconda3/etc/profile.d/conda.sh && conda activate stb    # any environment with stb_suite
#      bash make_reference.sh
# =====================================================================================================
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
REF="$HERE/reference_results"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
cp -a "$HERE/runs/." "$WORK/"
clean() { sed 's/\x1b\[[0-9;]*m//g' | sed "s#$WORK#<work>#g" | sed -E 's/^(Date\/time *: ).*/\1<date>/'; }
report() {       # report <case> <stb-nativecharges args...>: runs in $WORK, stores the plain-text report
    local case="$1"; shift
    mkdir -p "$REF/$case"
    echo "  [run]  $case: stb-nativecharges $*"
    (cd "$WORK" && stb-nativecharges --no-intro "$@" 2>&1) | clean > "$REF/$case/report.txt" || true
}

echo "== 1. Water, DZP: the whole report"
report h2o_dzp --path h2o_basis/DZP
echo "== 2. The basis series: Mulliken vs Hirshfeld vs Voronoi"
report h2o_basis --path 'h2o_basis/*' --csv h2o_basis_charges.csv
cp "$WORK/h2o_basis_charges.csv" "$REF/h2o_basis/"
report h2o_sz --path h2o_basis/SZ
echo "== 3. An ionic solid: NaCl"
report nacl --path nacl
echo "== 4. Spin: the O2 triplet"
report o2 --path o2
echo "== 5. A log without charge blocks"
mkdir -p "$WORK/nocharges"; cp "$WORK/nacl/nacl.XV" "$WORK/nocharges/"
sed '/^Charge\./d; /^WriteMullikenPop/d' "$WORK/nacl/calc.fdf" > "$WORK/nocharges/calc.fdf"
python3 - "$WORK/nacl/calc.out" "$WORK/nocharges/calc.out" <<'PY'
# the same log with the three population blocks cut out, as SIESTA writes it without the flags
import sys, re
text = open(sys.argv[1]).read()
for h in ("Mulliken", "Hirshfeld", "Voronoi"):
    text = re.sub(rf"\n{h} Atomic Populations:\n.*?\n Total[^\n]*\n", "\n", text, flags=re.S)
open(sys.argv[2], "w").write(text)
PY
report guard_nocharges --path nocharges
echo "== 6. The interactive menu (stb-suite, code 2.19)"
mkdir -p "$REF/menu"
# path glob, no label, details n, no csv, no report, Enter, 0 = exit
(cd "$WORK" && printf '2.19\nh2o_basis/*\n\nn\n\nn\n\n0\n' | stb-suite 2>&1) | clean \
    | sed -n '/STB-NATIVECHARGES REPORT/,/report complete/p' > "$REF/menu/report.txt"
echo "done -- reference_results/<case>/"
