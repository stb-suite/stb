#!/bin/bash
# Does an archive made with --mode inputs really reproduce the run? Archives runs/si_relax with --mode inputs, unpacks
# it in a scratch folder, runs SIESTA there and compares the total energy and the relaxed cell with the original
# calc.out. Needs SIESTA 5.x + mpirun (~15 s). Result in reference_results/rerun/.
#     bash rerun_from_archive.sh [NPROC]
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
NP=${1:-4}
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
PSML="$(python3 -c 'from stb.core.pseudopotentials import bank_path; print(bank_path("dojo"))')/Si.psml"
cp -a "$HERE/runs/si_relax" "$WORK/si_relax"; cp "$PSML" "$WORK/si_relax/"
(cd "$WORK" && stb-archive --no-intro --path si_relax --mode inputs > archive_report.txt 2>&1)
mkdir "$WORK/unpacked" && tar -xzf "$WORK/si_relax_archive.tar.gz" -C "$WORK/unpacked"
ls "$WORK/unpacked/si_relax"
(cd "$WORK/unpacked/si_relax" && OMP_NUM_THREADS=1 mpirun -np "$NP" siesta calc.fdf < /dev/null > calc.out 2> siesta.err)
mkdir -p "$HERE/reference_results/rerun"
python3 - "$HERE/runs/si_relax/calc.out" "$WORK/unpacked/si_relax/calc.out" "$WORK/unpacked/si_relax" <<'PY' \
    | tee "$HERE/reference_results/rerun/comparison.txt"
import sys, os
from stb.status import parse_log
a, b, folder = parse_log(sys.argv[1]), parse_log(sys.argv[2]), sys.argv[3]
print(f"files unpacked from the inputs-only archive: {', '.join(sorted(f for f in os.listdir(folder) if f in ('calc.fdf','structure.fdf','bandlines.txt','Si.psml','submit.slurm','MANIFEST.txt','SHA256SUMS')))}")
for name, info in (("original", a), ("rerun   ", b)):
    print(f"{name}: E_total {info['final']['Total']:.6f} eV, {len(info['steps']) - 1} CG moves, "
          f"relaxed: {info['relaxed']}, Fermi {info['fermi']:.6f} eV")
print(f"difference: {1000 * (b['final']['Total'] - a['final']['Total']):+.3f} meV")
PY
