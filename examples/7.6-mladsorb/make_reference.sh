#!/bin/bash
# =====================================================================================================
#  make_reference.sh -- regenerates reference_results/ of the 7.6 example (stb-mladsorb)
#
#  Every case runs into ./reference_run/<case>/ (skipped when its report exists; FORCE=1 reruns), then the
#  report and the small data files are copied to ./reference_results/<case>/. GPU by default (about 25 minutes on
#  one RTX 3060); DEVICE=cpu works, several times slower.
#
#  The model is the round-3 silicon model of the 7.5 example (../7.5-mlneb/model_round3/, not copied). Its committee
#  (the three round-3 members) is not shipped -- 11 MB each: set MLFF_RUN to the run folder of the 7.1/5.7 silicon
#  workflow (the folder with mlff_workflow.json) to evaluate it in the main case; without it the main case runs with
#  --no-committee.
#
#  The SIESTA relaxations of section 7 are not run here (hours on a CPU): siesta/ holds their inputs and
#  reference_results/siesta/ their outputs; siesta_compare.py compares them with the model.
#
#      source ~/.anaconda3/etc/profile.d/conda.sh && conda activate stb    # any environment with stb_suite[ml]
#      MLFF_RUN=~/run/stb/5.7/mlff_si/mlff_runs bash make_reference.sh
# =====================================================================================================
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
DEVICE="${DEVICE:-cuda}"
RUN="$HERE/reference_run"
OUT="$HERE/reference_results"
R2="$HERE/../7.2-mlmd/model"           # round 2 of the 7.1 silicon run (the model of 7.2-7.4)
R3="$HERE/../7.5-mlneb/model_round3"   # round 3: round 2 + the vacancy-hop band of 5.7
mkdir -p "$RUN" "$OUT"
export MPLBACKEND=Agg

r2=(--custom-model "$R2/mlff_model_compiled.model" --custom-model-card "$R2/mlff_model_card.json")
r3=(--custom-model "$R3/mlff_model_compiled.model" --custom-model-card "$R3/mlff_model_card.json")
common=(--no-intro --device "$DEVICE" --adsorbate Si --bulk-reference "$HERE/si_bulk.fdf" --save-data --save-report)
quick=(--skip-foundation-comparison --no-committee)
committee=(--no-committee)
if [ -n "$MLFF_RUN" ] && [ -f "$MLFF_RUN/mlff_workflow.json" ]; then
    committee=(--mlff-workflow "$MLFF_RUN/mlff_workflow.json")
fi

run() {   # run <case> <slab> <args...>
    local case="$1" slab="$2"; shift 2
    if [ -f "$RUN/$case/stb_mladsorb_report.txt" ] && [ -z "$FORCE" ]; then
        echo "  [skip] $case (report exists; FORCE=1 to rerun)"; return
    fi
    echo "  [run]  $case: stb-mladsorb --slab $slab $*"
    local t0=$SECONDS
    (cd "$RUN" && stb-mladsorb --slab "$HERE/$slab" "${common[@]}" -o "$case" "$@" > "$case.log" 2>&1) || {
        echo "  [FAIL] $case -- see $RUN/$case.log"; exit 1; }
    echo "         $((SECONDS - t0)) s"
}

echo "== 1. The main case: round-3 model, top 4 Ang relaxed, committee + foundation comparison, every output"
run main si111_2x2.fdf "${r3[@]}" --relax-depth 4 "${committee[@]}" --save-gnuplot --save-xyz

echo "== 2. How much of the substrate relaxes"
run depth_0 si111_2x2.fdf "${r3[@]}" "${quick[@]}" --relax-depth 0
run depth_1 si111_2x2.fdf "${r3[@]}" "${quick[@]}" --relax-depth 1
run depth_7 si111_2x2.fdf "${r3[@]}" "${quick[@]}" --relax-depth 7

echo "== 3. Starting distances and binding curves at fixed height"
run start_heights si111_2x2.fdf "${r3[@]}" "${quick[@]}" --relax-depth 4 --start-heights 1.5 2.0 2.5 3.0
run binding_curve si111_2x2.fdf "${r3[@]}" "${quick[@]}" --relax-depth 4 --height-sweep 1.0 3.0 0.25 \
    --save-gnuplot --save-xyz

echo "== 4. Slab thickness and coverage"
run thick_4bl si111_2x2_4bl.fdf "${r3[@]}" "${quick[@]}" --relax-depth 4
run cell_3x3 si111_3x3.fdf "${r3[@]}" "${quick[@]}" --relax-depth 4

echo "== 5. Other models: round 2, and the three MACE-MP-0 foundation models"
run round2 si111_2x2.fdf "${r2[@]}" "${quick[@]}" --relax-depth 4
for m in small medium large; do run "foundation_$m" si111_2x2.fdf --model "$m" --relax-depth 4; done

echo "== collect into reference_results/"
for d in "$RUN"/*/; do
    c=$(basename "$d")
    mkdir -p "$OUT/$c"
    for f in stb_mladsorb_report.txt mladsorb.png mladsorb_*.dat rank*.fdf mladsorb_minima.xyz mladsorb_slab.xyz; do
        [ -e "$d/$f" ] || compgen -G "$d/$f" > /dev/null && cp $d/$f "$OUT/$c/" 2>/dev/null || true
    done
    [ -d "$d/plot" ] && { rm -rf "$OUT/$c/plot"; cp -r "$d/plot" "$OUT/$c/plot"; }
done
cp "$OUT/main/mladsorb.png" figures/main_dashboard.png
echo "done -- python3 results_figures.py draws the summary figures; python3 siesta_compare.py the SIESTA comparison"
