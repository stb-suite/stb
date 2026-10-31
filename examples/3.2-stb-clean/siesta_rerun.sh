#!/bin/bash
# =====================================================================================================
#  siesta_rerun.sh -- the experiment of section 6.4: run SIESTA again in two cleaned copies of siesta_run/
#
#    rerun_fresh/     stb-clean with the default --keep        -> the next run starts from scratch
#    rerun_continue/  stb-clean --keep ... .DM .XV             -> it continues from the relaxed geometry/DM
#
#  Needs SIESTA (5.x) and mpirun; ~25 s on 4 cores. The calc.out and CLOCK of each go to reference_results/.
#      bash siesta_rerun.sh [NPROC]
# =====================================================================================================
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
NP="${1:-4}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
PSML="$(python3 -c 'from stb.core.pseudopotentials import bank_path; print(bank_path("dojo"))')/Si.psml"
command -v siesta > /dev/null || { echo "siesta not found on PATH"; exit 1; }
for case in fresh continue; do
    d="$WORK/$case"; cp -a "$HERE/siesta_run" "$d"; cp "$PSML" "$d/"
    if [ $case = fresh ]; then keep=(.fdf .psf .psml .sh); else keep=(.fdf .psf .psml .sh .DM .XV); fi
    stb-clean --no-intro --path "$d" -n --no-archive --keep "${keep[@]}" > /dev/null
    echo "  [siesta] rerun_$case (kept: ${keep[*]})"
    (cd "$d" && OMP_NUM_THREADS=1 mpirun -np "$NP" siesta calc.fdf < /dev/null > calc.out 2> siesta.err)
    mkdir -p "$HERE/reference_results/rerun_$case"
    cp "$d/calc.out" "$d/CLOCK" "$HERE/reference_results/rerun_$case/"
done
echo "done -- reference_results/rerun_fresh/, reference_results/rerun_continue/"
