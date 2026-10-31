#!/bin/bash
# =====================================================================================================
#  make_reference.sh -- regenerates reference_results/<case>/ of the 3.4 example (stb-siesta2wtb)
#
#  Every case runs on a copy of siesta_runs/ (two real, finished SIESTA 5.4 runs with SaveHS: an hBN monolayer and
#  ferromagnetic bcc Fe) in a scratch folder; only the reports and the first lines of the written files are kept
#  (scratch paths written as <work>, dates as <date>). About 30 s, no SIESTA (siesta_runs/run_all.sh made the runs).
#
#      source ~/.anaconda3/etc/profile.d/conda.sh && conda activate stb    # any environment with stb_suite
#      bash make_reference.sh
# =====================================================================================================
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
REF="$HERE/reference_results"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
clean() { sed 's/\x1b\[[0-9;]*m//g' | sed "s#$WORK#<work>#g" \
          | sed -E 's/^(Date\/time *: ).*/\1<date>/; s/(written )[0-9-]+ [0-9:]+\)/\1<date>)/; s/ in [0-9.]+ s/ in <t> s/; s/, [0-9.]+ s\)/, <t> s)/'; }
report() {       # report <case> <run folder> <stb-siesta2wtb args...>: a fresh copy of the run, the plain-text report
    local case="$1" run="$2"; shift 2
    rm -rf "$WORK/$case"; mkdir -p "$REF/$case"; cp -a "$HERE/siesta_runs/$run" "$WORK/$case"
    echo "  [run]  $case: stb-siesta2wtb $*"
    (cd "$WORK/$case" && stb-siesta2wtb --no-intro "$@" 2>&1) | clean > "$REF/$case/report.txt" || true
}
heads() {        # heads <case> <spin code>: the first lines of both written files
    head -n 14 "$WORK/$1/tb-$2.ham" > "$REF/$1/tb_head.txt"
    head -n 8 "$WORK/$1/basis_set-$2.basis" > "$REF/$1/basis_head.txt"
}

echo "== 1. hBN monolayer (NP), Fermi level from the log"
report hbn hbn -i calc.fdf -o calc.out
heads hbn NP

echo "== 2. bcc Fe (SP)"
report fe fe -i calc.fdf -o calc.out
heads fe SP
awk 'NR > 9 && $4 == 20 && $5 == 20 && $1 == 0 && $2 == 0 && $3 == 0' "$WORK/fe/tb-SP.ham" > "$REF/fe/down_block_line.txt"

echo "== 3. The Fermi level: none, by hand, 0"
report hbn_nofermi hbn -i calc.fdf --dry-run
report hbn_f0 hbn -i calc.fdf -f 0 --no-check
head -n 3 "$WORK/hbn_f0/tb-NP.ham" > "$REF/hbn_f0/tb_head.txt"

echo "== 4. Guard rails"
report guard_nohsx hbn -i calc.fdf -o calc.out
rm -rf "$WORK/guard_nohsx"; cp -a "$HERE/siesta_runs/hbn" "$WORK/guard_nohsx"; rm "$WORK/guard_nohsx/hbn.HSX"
(cd "$WORK/guard_nohsx" && stb-siesta2wtb --no-intro -i calc.fdf -o calc.out 2>&1) | clean > "$REF/guard_nohsx/report.txt" || true
report guard_badlog hbn -i calc.fdf -o siesta.out

echo "== 5. The interactive menu (stb-suite, code 3.4)"
rm -rf "$WORK/menu"; cp -a "$HERE/siesta_runs/hbn" "$WORK/menu"; mkdir -p "$REF/menu"
# fdf, log, no manual E_F, no dry run, check yes (Enter), no report, Enter, 0 = exit
(cd "$WORK/menu" && printf '3.4\ncalc.fdf\ncalc.out\n\nn\n\nn\n\n0\n' | stb-suite 2>&1) | clean \
    | sed -n '/STB-SIESTA2WTB REPORT/,/^\[9\] WARNINGS/p' > "$REF/menu/report.txt"

echo "done -- reference_results/<case>/"
