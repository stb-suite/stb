#!/bin/bash
# =====================================================================================================
#  make_reference.sh -- regenerates reference_results/<case>/report.txt of the 3.6 example (stb-status)
#
#  Every case runs on a copy of project/ (real, finished SIESTA 5.4 runs; see run_project.sh) in a scratch folder.
#  Two states are made there from a copy of relax/: a run that is still going (log cut, written just now) and one
#  that stopped (the same, written hours ago). Dates and ages are replaced by <date>/<age>. Seconds, no SIESTA.
#
#      source ~/.anaconda3/etc/profile.d/conda.sh && conda activate stb    # any environment with stb_suite
#      bash make_reference.sh
# =====================================================================================================
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
REF="$HERE/reference_results"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
cp -a "$HERE/project" "$WORK/project"
clean() { sed 's/\x1b\[[0-9;]*m//g' | sed "s#$WORK#<work>#g" | sed -E 's/^(Date\/time *: ).*/\1<date>/; s/written [0-9.]+ (s|min|h) ago/written <age> ago/g'; }
report() {       # report <case> <stb-status args...>: runs in $WORK, stores the plain-text report
    local case="$1"; shift
    mkdir -p "$REF/$case"
    echo "  [run]  $case: stb-status $*"
    (cd "$WORK" && stb-status --no-intro "$@" 2>&1) | clean > "$REF/$case/report.txt" || true
}
cut_run() {      # cut_run <dest>: relax/ stopped half-way -- log cut after move 4, no end of run anywhere
    rm -rf "$WORK/$1"; cp -a "$WORK/project/relax" "$WORK/$1"
    local n; n=$(grep -n "Begin CG opt. move =      4" "$WORK/$1/calc.out" | cut -d: -f1)
    head -n $((n + 30)) "$WORK/$1/calc.out" > "$WORK/$1/x" && mv "$WORK/$1/x" "$WORK/$1/calc.out"
    rm -f "$WORK/$1/MESSAGES" "$WORK/$1/0_NORMAL_EXIT"
    grep -v "End of run" "$WORK/$1/CLOCK" > "$WORK/$1/x" && mv "$WORK/$1/x" "$WORK/$1/CLOCK"
}

echo "== 1. One folder: a converged relaxation"
report relax --path project/relax
echo "== 2. A single point with bands and the density: what can be done next"
report bulk_sp --path project/bulk_sp
echo "== 3. A relaxation that ran out of moves"
report relax_short --path project/relax_short
echo "== 4. A failed SCF"
report scf_fail --path project/scf_fail
echo "== 5. Running, or stopped?"
cut_run running
report running --path running
cut_run stopped; touch -d "3 hours ago" "$WORK/stopped/calc.out"
report stopped --path stopped
echo "== 6. A spin-polarized metal, and an MD run"
report fe_spin --path project/fe_spin
report md --path project/md
echo "== 7. A series: relative energies"
report strain --path 'project/strain/a_*'
echo "== 8. The whole project"
report project --path project --recursive --csv project_status.csv
cp "$WORK/project_status.csv" "$REF/project/"

echo "== 9. The interactive menu (stb-suite, code 3.6)"
mkdir -p "$REF/menu"
# path, recursive y, no label, no csv, no report, Enter, 0 = exit
(cd "$WORK" && printf '3.6\nproject\ny\n\n\nn\n\n0\n' | stb-suite 2>&1) | clean \
    | sed -n '/STB-STATUS REPORT/,/Status check complete/p' > "$REF/menu/report.txt"
echo "done -- reference_results/<case>/"
