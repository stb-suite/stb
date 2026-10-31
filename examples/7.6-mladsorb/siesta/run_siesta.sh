#!/bin/bash
# Runs SIESTA in every folder given (skips finished ones): NP MPI processes (default 6), one folder after another.
#   bash siesta/run_siesta.sh <folder> [<folder> ...]
# siesta must be on PATH (e.g. conda install -c conda-forge siesta). On the machine of the example each relaxation took
# 30-90 minutes on 6 cores.
NP="${NP:-6}"
for d in "$@"; do
    grep -q "End of run" "$d/siesta.out" 2>/dev/null && { echo "$d: finished, skipped"; continue; }
    echo "$d: $(date +%H:%M:%S)"
    (cd "$d" && OMP_NUM_THREADS=1 mpirun -np "$NP" siesta < calc.fdf > siesta.out 2> siesta.err)
done
echo "all done $(date +%H:%M:%S)"
