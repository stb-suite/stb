#!/bin/bash
# Runs every case under CASES with NRANKS MPI ranks, each in <case>/np<NRANKS>/.
#   bash run_cases.sh CASES NRANKS [TIMEOUT_S]
# Environment: SIESTA (default 'siesta'), MPIRUN (default 'mpirun'). With NRANKS = 1 SIESTA runs
# without mpirun unless MPIRUN_SERIAL=1. Cases already finished (0_NORMAL_EXIT) are skipped.
CASES="${1:?cases folder}"
NP="${2:?number of ranks}"
TIMEOUT="${3:-600}"
SIESTA="${SIESTA:-siesta}"
MPIRUN="${MPIRUN:-mpirun}"
n=0
for case in $(find "$CASES" -name calc.fdf -not -path '*/np*/*' -printf '%h\n' | sort); do
    run="$case/np$NP"
    [ -f "$run/0_NORMAL_EXIT" ] && continue
    mkdir -p "$run"
    cp "$case"/calc.fdf "$case"/*.psml "$run"/
    ( cd "$run" || exit
      if [ "$NP" -eq 1 ] && [ "${MPIRUN_SERIAL:-0}" != 1 ]; then
          timeout "$TIMEOUT" "$SIESTA" < calc.fdf > calc.out 2> calc.err
      else
          timeout "$TIMEOUT" "$MPIRUN" -np "$NP" "$SIESTA" < calc.fdf > calc.out 2> calc.err
      fi
      echo $? > exit_code )
    n=$((n+1))
done
echo "$n case(s) run with $NP rank(s)"
