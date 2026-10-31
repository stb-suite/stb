#!/bin/bash
# Runs the three SIESTA AIMD calculations of the 6.5 example (SIESTA 5.x + mpirun; ~20 min, the three in parallel on 4 processes each).
# Pseudopotentials: the PseudoDojo .psml files of the suite's bundled bank (copied in; not shipped).
#     bash run_all.sh [NPROC per run]
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
NP=${1:-4}
BANK="$(python3 -c 'from stb.core.pseudopotentials import bank_path; print(bank_path("dojo"))')"
for d in water_nvt si_nve si_npt; do
    (
    cd "$HERE/$d"
    for el in $(awk '/%block ChemicalSpeciesLabel/{f=1;next} /%endblock/{f=0} f{print $3}' calc.fdf); do cp "$BANK/$el.psml" .; done
    echo "[siesta] $d"
    OMP_NUM_THREADS=1 mpirun -np "$NP" siesta calc.fdf < /dev/null > calc.out 2> siesta.err
    ) &
done
wait

# Keep only what the example uses: the input, the log, the trajectory and its companions.
prune() {   # prune <folder> <files to keep...>
    local d="$1"; shift
    (cd "$HERE/$d" && for f in * ; do case " $* " in *" $f "*) ;; *) rm -rf "$f";; esac; done)
}
prune water_nvt calc.fdf calc.out water.ANI water.MDE water.XV
prune si_nve    calc.fdf calc.out si.ANI si.MDE si.XV
prune si_npt    calc.fdf calc.out si.ANI si.MDE si.XV
