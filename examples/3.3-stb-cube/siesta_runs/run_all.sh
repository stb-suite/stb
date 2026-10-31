#!/bin/bash
# Runs the three SIESTA calculations of the 6.3 example (SIESTA 5.x + mpirun; a few minutes on 4 processes).
# Pseudopotentials: the PseudoDojo .psml files of the suite's bundled bank (copied in by this script; not shipped).
#     bash run_all.sh [NPROC]
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
NP=${1:-4}
BANK="$(python3 -c 'from stb.core.pseudopotentials import bank_path; print(bank_path("dojo"))')"
for d in si h2o o2; do
    cd "$HERE/$d"
    for el in $(awk '/%block ChemicalSpeciesLabel/{f=1;next} /%endblock/{f=0} f{print $3}' calc.fdf); do cp "$BANK/$el.psml" .; done
    echo "[siesta] $d"
    OMP_NUM_THREADS=1 mpirun -np "$NP" siesta calc.fdf < /dev/null > calc.out 2> siesta.err
done

# Keep only what the example uses: the input, the log, the geometry and the grids (o2.HSX for the spin check).
prune() {   # prune <folder> <files to keep...>
    local d="$1"; shift
    (cd "$HERE/$d" && for f in * ; do case " $* " in *" $f "*) ;; *) rm -rf "$f";; esac; done)
}
prune si  calc.fdf calc.out si.XV si.RHO si.DRHO si.BADER si.TOCH si.IOCH si.VT si.VH si.VNA
prune h2o calc.fdf calc.out h2o.XV h2o.RHO h2o.DRHO h2o.VH
prune o2  calc.fdf calc.out o2.XV o2.RHO o2.HSX
