#!/bin/bash
# Runs the two SIESTA calculations of the 6.4 example (SIESTA 5.x + mpirun; about a minute on 4 processes).
# Pseudopotentials: the PseudoDojo .psml files of the suite's bundled bank (copied in; not shipped).
#     bash run_all.sh [NPROC]
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
NP=${1:-4}
BANK="$(python3 -c 'from stb.core.pseudopotentials import bank_path; print(bank_path("dojo"))')"
for d in hbn fe; do
    cd "$HERE/$d"
    for el in $(awk '/%block ChemicalSpeciesLabel/{f=1;next} /%endblock/{f=0} f{print $3}' calc.fdf); do cp "$BANK/$el.psml" .; done
    echo "[siesta] $d"
    OMP_NUM_THREADS=1 mpirun -np "$NP" siesta calc.fdf < /dev/null > calc.out 2> siesta.err
done

# Keep only what the example uses: the input, the log, the geometry, the Hamiltonian, eigenvalues and bands.
prune() {   # prune <folder> <files to keep...>
    local d="$1"; shift
    (cd "$HERE/$d" && for f in * ; do case " $* " in *" $f "*) ;; *) rm -rf "$f";; esac; done)
}
prune hbn  calc.fdf calc.out hbn.XV hbn.HSX hbn.EIG hbn.KP hbn.bands
prune fe   calc.fdf calc.out fe.XV fe.HSX fe.EIG fe.KP fe.bands
