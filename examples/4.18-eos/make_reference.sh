#!/bin/bash
# Regenerates reference_results/si/: every variant of the Si equation-of-state study of README
# section 5 -- Stage 1 (stb-eosInputs) + SIESTA in each of its 11 volume folders. ~10 minutes on
# 2 cores (mpirun -np 2). The pseudopotential (-p dojo) is copied into each folder and removed
# afterwards, with the SIESTA by-products, to keep the repository small.
#
#   bash make_reference.sh [work_dir]      (default: ./reference_run, then trimmed into reference_results/si)
set -e
cd "$(dirname "${BASH_SOURCE[0]}")"
WORK="${1:-reference_run}"
mkdir -p "$WORK"

# name | k-grid | Mesh.Cutoff (Ry) | basis size | PAO.EnergyShift (Ry)
VARIANTS="
base 10 320 DZP 0.02
k04  4  320 DZP 0.02
k06  6  320 DZP 0.02
k08  8  320 DZP 0.02
k14  14 320 DZP 0.02
m150 10 150 DZP 0.02
m200 10 200 DZP 0.02
m450 10 450 DZP 0.02
m600 10 600 DZP 0.02
dzp_s010 10 320 DZP 0.01
dzp_s005 10 320 DZP 0.005
tzp_s020 10 320 TZP 0.02
tzp_s010 10 320 TZP 0.01
tzp_s005 10 320 TZP 0.005
"
echo "$VARIANTS" | while read name k mesh basis shift; do
    [ -z "$name" ] && continue
    [ -d "$WORK/$name/vol_0.00" ] && [ -s "$WORK/$name/vol_8.00/calc.out" ] && { echo "$name: done"; continue; }
    mkdir -p "$WORK/$name"
    sed "s/^kgrid.MonkhorstPack.*/kgrid.MonkhorstPack     [$k  $k  $k]/; s/^Mesh.Cutoff.*/Mesh.Cutoff             $mesh Ry/; \
         s/^PAO.BasisSize.*/PAO.BasisSize           $basis/; s/^PAO.EnergyShift.*/PAO.EnergyShift         $shift Ry/" \
        calc.fdf > "$WORK/$name/calc_template.fdf"
    cp structure.fdf "$WORK/$name/"
    (cd "$WORK/$name" && stb-eosInputs -f structure.fdf -c calc_template.fdf -p dojo --n-volumes 11 \
         --strain-range 8.0 -o . --save-report --no-intro < /dev/null > stage1.log 2>&1)
    for d in "$WORK/$name"/vol_*; do
        (cd "$d" && mpirun -np 2 siesta calc.fdf < /dev/null > calc.out 2> calc.err)
    done
    echo "$name: done"
done

# trim into reference_results/si
for name in $(ls "$WORK"); do
    out="reference_results/si/$name"; mkdir -p "$out"
    cp "$WORK/$name/eos_scan.json" "$WORK/$name/stb_eosInputs_report.txt" "$WORK/$name/calc_template.fdf" "$out/"
    for d in "$WORK/$name"/vol_*; do
        v=$(basename "$d"); mkdir -p "$out/$v"
        cp "$d/structure.fdf" "$d/calc.fdf" "$d/config_extra.fdf" "$d/calc.out" "$out/$v/"
    done
done
echo "reference_results/si/ written"
