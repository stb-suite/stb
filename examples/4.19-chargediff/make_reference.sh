#!/bin/bash
# Regenerates reference_results/ of the 4.19 example (README section 6): every SIESTA run behind the worked
# example, with the stb-chargediffPrep / stb-chargediffAnalysis commands of the tutorial.
#
#   1. relaxations (CG, fixed cell): the water dimer, ammonia borane, Li on graphene, and the isolated
#      monomers H2O, NH3, BH3 (for the deformation energy);
#   2. one charge-density-difference set per case: stb-chargediffPrep, SIESTA in combined/ and every
#      fragment folder, the fragments once more WITHOUT ghost atoms (build_inputs.py alone: the uncorrected,
#      BSSE-affected energy), stb-chargediffAnalysis --save-data --save-report;
#   3. the water-dimer studies: MeshCutoff 150-600 Ry, basis DZP/TZP x PAO.EnergyShift, O...O distance 2.5-5 Ang.
#
# About 1.5 h on 2 cores (mpirun -np 2), most of it the Li/graphene relaxation.
#
#   bash make_reference.sh [work_dir]     (default: ./reference_run, then trimmed into reference_results/)
set -e
cd "$(dirname "${BASH_SOURCE[0]}")"
HERE="$(pwd)"
WORK="$(realpath -m "${1:-reference_run}")"
NP="${NP:-2}"
DOJO="$(python3 -c 'from stb.core.pseudopotentials import bank_path; print(bank_path("dojo"))')"
mkdir -p "$WORK"

siesta_in() {   # $1 = folder; runs calc.fdf there unless a finished calc.out exists
    if grep -q "End of run" "$1/calc.out" 2>/dev/null; then return; fi
    (cd "$1" && mpirun -np "$NP" siesta calc.fdf < /dev/null > calc.out 2> calc.err)
}

relax() {   # $1 = name, $2 = start structure, $3 = template, $4.. = species; extra directives in $EXTRA
    local d="$WORK/relax/$1"
    mkdir -p "$d"
    cp "$2" "$d/structure.fdf"
    cp "$3" "$d/template.fdf"
    for s in "${@:4}"; do cp "$DOJO/$s.psml" "$d/"; done
    { echo "MD.TypeOfRun        CG"; echo "MD.Steps            300"; echo "MD.MaxForceTol      0.01 eV/Ang"
      echo "MD.VariableCell     false"; printf "%b" "${EXTRA:-}"; echo "%include template.fdf"; } > "$d/calc.fdf"
    siesta_in "$d"
    # the relaxed geometry, atom order kept (fragment splits by count depend on it)
    stb-translate --no-intro -if siesta -i "$d/siesta.STRUCT_OUT" -of fdf -o "$d/relaxed.fdf" < /dev/null > /dev/null
    echo "relax $1: $(grep -c 'SCF cycle converged' "$d/calc.out") CG steps"
}

cdd() {   # $1 = set name, $2 = structure, $3 = fragment sizes, $4 = template, $5 = MeshCutoff (Ry), $6.. = prep flags
    local d="$WORK/$1"
    if [ ! -f "$d/chargediff_manifest.json" ]; then
        stb-chargediffPrep -s "$2" --fragment-sizes "$3" --calc "$4" --mesh-cutoff "$5" -p dojo -o "$d" \
            "${@:6}" --save-report --no-intro < /dev/null > /dev/null
    fi
    for f in "$d"/combined "$d"/fragment[0-9]; do siesta_in "$f"; done
    for f in "$d"/fragment[0-9]; do
        [ -d "${f}_alone" ] || python3 "$HERE/build_inputs.py" alone "$f" "${f}_alone"
        siesta_in "${f}_alone"
    done
    (cd "$d" && stb-chargediffAnalysis --dir . -o analysis --save-data --save-report --no-intro < /dev/null > /dev/null)
    echo "cdd $1: done"
}

# ---- 1. relaxations ------------------------------------------------------------------------------------------
relax water_dimer water_dimer.fdf calc.fdf O H
relax nh3bh3 nh3bh3.fdf calc.fdf N B H
EXTRA="Slab.DipoleCorrection .true.\n" relax li_graphene li_graphene.fdf calc_graphene.fdf C Li
R="$WORK/relax"
python3 build_inputs.py monomer "$R/water_dimer/relaxed.fdf" 1 3 "$WORK/h2o_start.fdf"
python3 build_inputs.py monomer "$R/nh3bh3/relaxed.fdf" 1 4 "$WORK/nh3_start.fdf"
python3 build_inputs.py monomer "$R/nh3bh3/relaxed.fdf" 5 4 "$WORK/bh3_start.fdf"
relax h2o "$WORK/h2o_start.fdf" calc.fdf O H
relax nh3 "$WORK/nh3_start.fdf" calc.fdf N H
relax bh3 "$WORK/bh3_start.fdf" calc.fdf B H

# ---- 2. the three cases --------------------------------------------------------------------------------------
cdd water "$R/water_dimer/relaxed.fdf" 3,3 calc.fdf 320
cdd nh3bh3 "$R/nh3bh3/relaxed.fdf" 4,4 calc.fdf 320
cdd li_graphene "$R/li_graphene/relaxed.fdf" 18,1 calc_graphene.fdf 300 --force-spin --force-dipole
# flat graphene before Li (the start structure, same settings) for the deformation energy of the sheet
G="$WORK/graphene_flat"
mkdir -p "$G"
python3 build_inputs.py monomer li_graphene.fdf 1 18 "$G/structure.fdf"
cp "$DOJO/C.psml" "$G/"
{ echo "MD.TypeOfRun CG"; echo "MD.Steps 0"; echo "MD.VariableCell false"; echo "Slab.DipoleCorrection .true."
  echo "MeshCutoff 300 Ry"; echo "%include template.fdf"; } > "$G/calc.fdf"
cp calc_graphene.fdf "$G/template.fdf"
siesta_in "$G"

# ---- 3. water-dimer studies ----------------------------------------------------------------------------------
for m in 150 200 450 600; do cdd water_mesh$m "$R/water_dimer/relaxed.fdf" 3,3 calc.fdf $m; done
# basis: the same geometry, another basis (PAO.BasisSize x PAO.EnergyShift)
for b in "DZP 0.02" "DZP 0.005" "TZP 0.01" "TZP 0.005"; do
    set -- $b
    tag="water_$(echo "$1" | tr 'A-Z' 'a-z')_s$(echo "$2" | tr -d '.')"
    sed "s/^PAO.BasisSize.*/PAO.BasisSize           $1/; s/^PAO.EnergyShift.*/PAO.EnergyShift         $2 Ry/" \
        calc.fdf > "$WORK/calc_$tag.fdf"
    cdd "$tag" "$R/water_dimer/relaxed.fdf" 3,3 "$WORK/calc_$tag.fdf" 320
done
# distance: the donor moved rigidly along O1->O2
for r in 2.50 2.60 2.70 2.90 3.00 3.20 3.50 4.00 5.00; do
    python3 build_inputs.py scan "$R/water_dimer/relaxed.fdf" $r "$WORK/dimer_r$r.fdf"
    cdd water_r$r "$WORK/dimer_r$r.fdf" 3,3 calc.fdf 320
done

# ---- trim into reference_results/ ----------------------------------------------------------------------------
OUT="$HERE/reference_results"
rm -rf "$OUT"
mkdir -p "$OUT/relax"
for d in "$WORK"/relax/*; do
    n=$(basename "$d"); mkdir -p "$OUT/relax/$n"
    cp "$d/calc.fdf" "$d/template.fdf" "$d/structure.fdf" "$d/relaxed.fdf" "$d/calc.out" "$OUT/relax/$n/"
done
mkdir -p "$OUT/graphene_flat"; cp "$G"/{calc.fdf,template.fdf,structure.fdf,calc.out} "$OUT/graphene_flat/"
for d in "$WORK"/water* "$WORK"/nh3bh3 "$WORK"/li_graphene; do
    n=$(basename "$d"); mkdir -p "$OUT/$n"
    cp "$d"/chargediff_manifest.json "$d"/chargediff_prep_report.txt "$OUT/$n/"
    mkdir -p "$OUT/$n/analysis"; cp "$d"/analysis/*results.json "$d"/analysis/*atoms.csv "$d"/analysis/*report.txt "$OUT/$n/analysis/"
    for f in "$d"/combined "$d"/fragment*; do
        b=$(basename "$f"); mkdir -p "$OUT/$n/$b"
        cp "$f"/{structure.fdf,calc.fdf,config_extra.fdf,calc.out} "$OUT/$n/$b/"
        # the density grids of the three main cases only (compressed): the studies keep their analysis results
        case "$n/$b" in */*_alone) ;; water/*|nh3bh3/*|li_graphene/*) [ -f "$f/siesta.RHO" ] && xz -T0 -c "$f/siesta.RHO" > "$OUT/$n/$b/siesta.RHO.xz";; esac
    done
done
echo "reference_results/ written"
