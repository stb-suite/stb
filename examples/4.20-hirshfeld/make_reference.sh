#!/bin/bash
# Regenerates reference_results/ of the 4.20 example (README section 6): every SIESTA run behind the worked
# example, with the stb-hirshfeldPrep / stb-hirshfeldIons / stb-hirshfeldAnalysis commands of the tutorial.
#
#   1. relaxations: H2O and N2O (CG, fixed cell), NaCl and MgO (CG, variable cell);
#   2. one Hirshfeld-I study per system: stb-hirshfeldPrep, SIESTA in combined/ and every neutral/<El>/,
#      stb-hirshfeldIons, SIESTA in every ions/<El>/<cation|anion>/, stb-hirshfeldAnalysis --save-data
#      --save-report --png (into <study>/analysis/), and stb-bader on combined/ (into <study>/bader/);
#   3. the H2O studies: MeshCutoff 200-600 Ry, the isolated-atom box 10-28 Ang, the basis (DZP/TZP x
#      PAO.EnergyShift) -- each a complete study of its own.
#
# About 1.5 h on 12 cores (3 SIESTA runs at a time, mpirun -np 4 each).
#
#   bash make_reference.sh [work_dir]     (default: ./reference_run, then trimmed into reference_results/)
#   NP=4 JOBS=3 bash make_reference.sh     MPI ranks per SIESTA run, SIESTA runs at a time
set -e
cd "$(dirname "${BASH_SOURCE[0]}")"
HERE="$(pwd)"
WORK="$(realpath -m "${1:-reference_run}")"
export NP="${NP:-4}"
JOBS="${JOBS:-3}"
DOJO="$(python3 -c 'from stb.core.pseudopotentials import bank_path; print(bank_path("dojo"))')"
mkdir -p "$WORK"

siesta_in() {   # $1 = folder; runs calc.fdf there unless a finished calc.out exists
    if grep -q "End of run" "$1/calc.out" 2>/dev/null; then return 0; fi
    (cd "$1" && OMP_NUM_THREADS=1 mpirun -np "$NP" siesta calc.fdf < /dev/null > calc.out 2> calc.err) || true
    # SIESTA 5.4.2 once aborted ('double free') after a converged SCF of the H atom at 600 Ry with 4 ranks;
    # the same input finished with 2
    if ! grep -q "End of run" "$1/calc.out"; then
        (cd "$1" && OMP_NUM_THREADS=1 mpirun -np 2 siesta calc.fdf < /dev/null > calc.out 2> calc.err) || true
    fi
    grep -q "End of run" "$1/calc.out" || { echo "SIESTA failed in $1" >&2; return 1; }
}
export -f siesta_in

siesta_all() {  # every folder given, JOBS at a time
    printf '%s\n' "$@" | xargs -P "$JOBS" -I{} bash -c 'siesta_in "$1"' _ {}
}

relax() {   # $1 = name, $2 = start structure, $3 = template, $4 = variable cell (yes/no), $5.. = species
    local d="$WORK/relax/$1"
    mkdir -p "$d"
    cp "$2" "$d/structure.fdf"
    cp "$3" "$d/template.fdf"
    for s in "${@:5}"; do cp "$DOJO/$s.psml" "$d/"; done
    { echo "MD.TypeOfRun        CG"; echo "MD.Steps            300"; echo "MD.MaxForceTol      0.01 eV/Ang"
      if [ "$4" = yes ]; then echo "MD.VariableCell     true"; echo "MD.MaxStressTol     0.05 GPa"
      else echo "MD.VariableCell     false"; fi
      echo "%include template.fdf"; } > "$d/calc.fdf"
    siesta_in "$d"
    stb-translate --no-intro -if siesta -i "$d/siesta.STRUCT_OUT" -of fdf -o "$d/relaxed.fdf" < /dev/null > /dev/null
    echo "relax $1: $(grep -c 'SCF cycle converged' "$d/calc.out") CG steps"
}

study() {   # $1 = study name, $2 = structure, $3 = template, $4 = MeshCutoff (Ry), $5 = vacuum box (Ang)
    local d="$WORK/$1"
    mkdir -p "$d"
    if [ ! -f "$d/hirshfeld_manifest.json" ]; then
        (cd "$d" && stb-hirshfeldPrep -s "$2" --calc "$3" -p dojo -O . --mesh-cutoff "$4" --vacuum "$5" \
            --save-report --no-intro < /dev/null > prep.log 2>&1)
    fi
    siesta_all "$d/combined" "$d"/neutral/*
    if [ ! -f "$d/hirshfeld_ions_manifest.json" ]; then
        (cd "$d" && stb-hirshfeldIons -O . -p dojo --save-report --no-intro < /dev/null > ions.log 2>&1)
    fi
    siesta_all "$d"/ions/*/*
    mkdir -p "$d/analysis"
    (cd "$d/analysis" && stb-hirshfeldAnalysis -O .. --save-data --save-report --png --no-intro \
        < /dev/null > analysis.log 2>&1)
    echo "study $1: $(grep -m1 -o 'Converged *: [a-z]*' "$d/analysis/hirshfeld_analysis_report.txt")"
}

bader() {   # $1 = study name: Bader charges of the same combined/ density
    local d="$WORK/$1"
    [ -f "$d/bader/stb_bader_report.txt" ] && return 0
    (cd "$d/combined" && stb-bader -l siesta --ref calc.out -o ../bader --save-report --no-cube --no-intro \
        < /dev/null > ../bader/bader.log 2>&1) || mkdir -p "$d/bader"
}

# ---- 1. relaxations ------------------------------------------------------------------------------------------
relax h2o h2o.fdf calc.fdf no O H
relax n2o n2o.fdf calc.fdf no N O
relax nacl nacl.fdf calc_bulk.fdf yes Na Cl
relax mgo mgo.fdf calc_bulk.fdf yes Mg O
R="$WORK/relax"

# ---- 2. the four systems -------------------------------------------------------------------------------------
study h2o "$R/h2o/relaxed.fdf" "$HERE/calc.fdf" 400 20
study n2o "$R/n2o/relaxed.fdf" "$HERE/calc.fdf" 400 20
study nacl "$R/nacl/relaxed.fdf" "$HERE/calc_bulk.fdf" 400 20
study mgo "$R/mgo/relaxed.fdf" "$HERE/calc_bulk.fdf" 400 20
for s in h2o n2o nacl mgo; do
    mkdir -p "$WORK/$s/bader"; bader $s
    # the same analysis converged to 1e-6 e: the reference the remaining-error estimate is checked against
    mkdir -p "$WORK/$s/analysis_tight"
    (cd "$WORK/$s/analysis_tight" && stb-hirshfeldAnalysis -O .. --tol 1e-6 --max-iter 300 --save-data \
        --save-report --no-intro < /dev/null > analysis.log 2>&1)
done

# ---- 3. H2O studies ------------------------------------------------------------------------------------------
for m in 200 300 600; do study h2o_mesh$m "$R/h2o/relaxed.fdf" "$HERE/calc.fdf" $m 20; done
for v in 10 14 28; do study h2o_box$v "$R/h2o/relaxed.fdf" "$HERE/calc.fdf" 400 $v; done
for b in "DZP 0.02" "DZP 0.005" "TZP 0.01" "TZP 0.005" "SZP 0.01"; do
    set -- $b
    tag="h2o_$(echo "$1" | tr 'A-Z' 'a-z')_s$(echo "$2" | tr -d '.')"
    sed "s/^PAO.BasisSize.*/PAO.BasisSize           $1/; s/^PAO.EnergyShift.*/PAO.EnergyShift         $2 Ry/" \
        calc.fdf > "$WORK/calc_$tag.fdf"
    study "$tag" "$R/h2o/relaxed.fdf" "$WORK/calc_$tag.fdf" 400 20
done

# ---- trim into reference_results/ ----------------------------------------------------------------------------
OUT="$HERE/reference_results"
rm -rf "$OUT"
mkdir -p "$OUT/relax"
for d in "$WORK"/relax/*; do
    n=$(basename "$d"); mkdir -p "$OUT/relax/$n"
    cp "$d/calc.fdf" "$d/template.fdf" "$d/structure.fdf" "$d/relaxed.fdf" "$d/calc.out" "$OUT/relax/$n/"
done
for d in "$WORK"/h2o* "$WORK"/n2o "$WORK"/nacl "$WORK"/mgo; do
    [ -d "$d/combined" ] || continue
    n=$(basename "$d"); mkdir -p "$OUT/$n/analysis"
    cp "$d"/hirshfeld_manifest.json "$d"/hirshfeld_ions_manifest.json "$d"/hirshfeld_prep_report.txt \
       "$d"/hirshfeld_ions_report.txt "$OUT/$n/"
    cp "$d"/analysis/hirshfeld_* "$OUT/$n/analysis/"
    [ -d "$d/analysis_tight" ] && { mkdir -p "$OUT/$n/analysis_tight"; cp "$d"/analysis_tight/hirshfeld_*.dat \
        "$d"/analysis_tight/hirshfeld_analysis_report.txt "$OUT/$n/analysis_tight/"; }
    [ -f "$d/bader/stb_bader_report.txt" ] && { mkdir -p "$OUT/$n/bader"; cp "$d/bader/stb_bader_report.txt" "$OUT/$n/bader/"; }
    for f in "$d"/combined "$d"/neutral/* "$d"/ions/*/*; do
        r=${f#"$d"/}; mkdir -p "$OUT/$n/$r"
        cp "$f"/{structure.fdf,calc.fdf,config_extra.fdf,calc.out} "$OUT/$n/$r/"
        # the density grids of the four systems only (SIESTA's density is exactly zero beyond the orbital radii:
        # ~100 kB each, compressed); the studies keep their inputs, logs and analysis results
        case "$n" in h2o|n2o|nacl|mgo) xz -T0 -c "$f/siesta.RHO" > "$OUT/$n/$r/siesta.RHO.xz";; esac
    done
done
echo "reference_results/ written"
