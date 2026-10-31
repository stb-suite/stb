#!/bin/bash
# Silicon's elastic constants with SIESTA at the level of theory the fine-tuned model was trained on (README section
# 5.11), through the SIESTA elastic workflow (5.2: stb-elasticInputs -> SIESTA -> stb-elasticAnalysis). Needs a
# 'siesta' executable and mpirun; the pseudopotential is the suite's PseudoDojo Si.psml (-p dojo, md5 878e1f00 as in
# the model card); about 6 minutes on 6 cores (37 runs of 8 atoms).
#
#   bash siesta_elastic.sh [work_dir]            run everything, then collect into reference_results/siesta_elastic/
#   bash siesta_elastic.sh [work_dir] collect    only collect the finished runs
#
# work_dir defaults to ./siesta_run; NP (default 6) sets the number of MPI processes. Finished runs are skipped.
set -e
cd "$(dirname "${BASH_SOURCE[0]}")"
HERE="$PWD"
WORK="$(mkdir -p "${1:-siesta_run}" && cd "${1:-siesta_run}" && pwd)"
STEP="${2:-all}"
NP="${NP:-6}"
OUT="$HERE/reference_results/siesta_elastic"

prepare() {
    cd "$WORK"
    cp "$HERE/si8.fdf" structure.fdf
    cp "$OUT/calc.fdf" calc.fdf
    # 1. xx and xy strains, -1 .. +1 %, 7 points: single points (clamped ions)
    [ -d clamped ] || stb-elasticInputs -s structure.fdf -c calc.fdf -p dojo --dirs all --max 1 --steps 7 \
        -o clamped --no-intro < /dev/null > /dev/null
    # 2. the same cells with the atoms relaxed (relaxed ions, what the model computes): config_extra.fdf replaced
    if [ ! -d relaxed ]; then
        cp -r clamped relaxed
        for d in relaxed/*/strain_*; do cp "$OUT/config_extra_relaxed.fdf" "$d/config_extra.fdf"; done
    fi
    # 3. isotropic strains, -2 .. +2 % in lattice, 9 points: B from E(V) and from the pressure
    [ -d hydro ] || stb-elasticInputs -s structure.fdf -c calc.fdf -p dojo --dirs hydro --max 2 --steps 9 \
        -o hydro --no-intro < /dev/null > /dev/null
}

run() {
    cd "$WORK"
    for d in clamped/*/strain_* relaxed/*/strain_* hydro/*/strain_*; do
        grep -q "End of run" "$d/siesta.out" 2>/dev/null && continue
        echo "$d: $(date +%H:%M:%S)"
        (cd "$d" && OMP_NUM_THREADS=1 mpirun -np "$NP" siesta < calc.fdf > siesta.out 2> siesta.err) \
            || echo "$d: FAILED (see $d/siesta.err)"
    done
}

collect() {
    cd "$WORK"
    for s in clamped relaxed; do
        stb-elasticAnalysis -d "$s" -f siesta.out --no-intro < /dev/null | sed 's/\x1b\[[0-9;]*m//g' \
            > "$OUT/elastic_analysis_$s.txt"
    done
    python3 "$HERE/siesta_compare.py" --collect "$WORK"
}

case "$STEP" in
    collect) collect ;;
    *) prepare; run; collect ;;
esac
