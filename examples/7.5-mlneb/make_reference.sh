#!/bin/bash
# =====================================================================================================
#  make_reference.sh -- regenerates reference_results/ of the 7.5 example (stb-mlneb)
#
#  Every case runs into ./reference_run/<case>/ (skipped when its report exists), then the reports and
#  data files are copied to ./reference_results/<case>/. GPU by default (about 20 minutes on one RTX 3060);
#  DEVICE=cpu works, but each 63-atom NEB then takes 15-30 minutes.
#
#  The two SIESTA-check cases (round2, round3) need no SIESTA: their single points were run once, in the
#  real loop of section 6 (reference_results/loop_round*/siesta_check/), and are copied into the folders
#  stb-mlneb writes -- after checking that the geometry is the same (the stamp stb_mlneb_image.xyz, 0.02
#  Ang, the tool's own tolerance). The second run of the same command then reads them.
#
#      source ~/.anaconda3/etc/profile.d/conda.sh && conda activate stb    # any environment with stb_suite[ml]
#      bash make_reference.sh
# =====================================================================================================
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
DEVICE="${DEVICE:-cuda}"
RUN="$HERE/reference_run"
OUT="$HERE/reference_results"
R2="$HERE/../7.2-mlmd/model"          # round 2 of the 7.1 silicon run (the model of 7.2-7.4)
R3="$HERE/model_round3"               # round 3: round 2 + 14 configs of the NEB band (section 6.5)
mkdir -p "$RUN" "$OUT"
export MPLBACKEND=Agg

r2=(--custom-model "$R2/mlff_model_compiled.model" --custom-model-card "$R2/mlff_model_card.json")
r3=(--custom-model "$R3/mlff_model_compiled.model" --custom-model-card "$R3/mlff_model_card.json")
ends=(--initial "$HERE/vac_initial.fdf" --final "$HERE/vac_final.fdf")
common=(--no-intro --device "$DEVICE" --save-data --save-report)

run() {   # run <case> <args...>
    local case="$1"; shift
    if [ -f "$RUN/$case/stb_mlneb_report.txt" ] && [ -z "$FORCE" ]; then
        echo "  [skip] $case (report exists; FORCE=1 to rerun)"; return
    fi
    echo "  [run]  $case: stb-mlneb $*"
    local t0=$SECONDS
    (cd "$RUN" && stb-mlneb "${ends[@]}" "${common[@]}" -o "$case" "$@" > "$case.log" 2>&1) || {
        echo "  [FAIL] $case -- see $RUN/$case.log"; exit 1; }
    echo "         $((SECONDS - t0)) s"
}

label_from() {   # label_from <case> <loop dir>: copy finished SIESTA outputs into matching folders
    python3 - "$RUN/$1/siesta_check" "$2/siesta_check" <<'EOF'
import os, shutil, sys
import numpy as np
from ase.io import read
new, old = sys.argv[1], sys.argv[2]
for name in sorted(os.listdir(new)):
    src = os.path.join(old, name)
    if not os.path.isdir(src):
        print(f"         {name}: no finished run in {old}"); continue
    a = read(os.path.join(new, name, "stb_mlneb_image.xyz"))
    b = read(os.path.join(src, "stb_mlneb_image.xyz"))
    d = np.abs(a.get_positions() - b.get_positions()).max()
    if d > 0.02:
        print(f"         {name}: geometry differs by {d:.3f} Ang from the finished run -- left unlabelled"); continue
    for f in ("siesta.out", "siesta.FA"):
        shutil.copy(os.path.join(src, f), os.path.join(new, name, f))
    print(f"         {name}: labelled (geometry within {d:.1e} Ang of the finished run)")
EOF
}

echo "== 1. Foundation models (MACE-MP-0 small / medium / large)"
for m in small medium large; do run "foundation_$m" --model "$m"; done

echo "== 2. Round-2 model: NEB + SIESTA check of every image (pass 1 writes, pass 2 reads)"
if [ ! -f "$RUN/round2/stb_mlneb_report.txt" ] || [ -n "$FORCE" ]; then
    rm -rf "$RUN/round2"
    run round2_pass1 "${r2[@]}" --skip-foundation-comparison --no-committee \
        --siesta-check --siesta-check-images all --calc "$HERE/calc.fdf" --pseudo-dir dojo --kgrid 3 3 3
    mv "$RUN/round2_pass1" "$RUN/round2"
    label_from round2 "$OUT/loop_round2"
    FORCE=1 run round2 "${r2[@]}" --skip-foundation-comparison --no-committee \
        --siesta-check --siesta-check-images all --calc "$HERE/calc.fdf" --pseudo-dir dojo --kgrid 3 3 3 --save-gnuplot
fi

echo "== 3. Round-3 model: NEB + foundation comparison + SIESTA check of the ends and the saddle"
if [ ! -f "$RUN/round3/stb_mlneb_report.txt" ] || [ -n "$FORCE" ]; then
    rm -rf "$RUN/round3"
    run round3_pass1 "${r3[@]}" --no-committee --siesta-check --calc "$HERE/calc.fdf" --pseudo-dir dojo --kgrid 3 3 3
    mv "$RUN/round3_pass1" "$RUN/round3"
    label_from round3 "$OUT/loop_round3"
    FORCE=1 run round3 "${r3[@]}" --no-committee --siesta-check --calc "$HERE/calc.fdf" --pseudo-dir dojo \
        --kgrid 3 3 3 --save-gnuplot
fi

echo "== 4. Number of images and IDPP (round-3 model)"
for n in 5 9 13; do run "images_$n" "${r3[@]}" --skip-foundation-comparison --no-committee --n-images "$n"; done
run no_idpp "${r3[@]}" --skip-foundation-comparison --no-committee --no-idpp

echo "== 5. A tighter force threshold (round-3 model, 7 and 9 images)"
run fmax_0.01 "${r3[@]}" --skip-foundation-comparison --no-committee --fmax 0.01
run fmax_0.01_images_9 "${r3[@]}" --skip-foundation-comparison --no-committee --fmax 0.01 --n-images 9

echo "== 6. Endpoints taken as they are (no pre-relaxation, round-3 model)"
run no_prerelax "${r3[@]}" --skip-foundation-comparison --no-committee --no-pre-relax-endpoints

echo "== collect into reference_results/"
for d in "$RUN"/*/; do
    c=$(basename "$d")
    [ -f "$d/stb_mlneb_report.txt" ] || continue
    mkdir -p "$OUT/$c"
    cp "$d"/stb_mlneb_report.txt "$d"/neb_band*.dat "$d"/neb_profile.png "$d"/neb_path*.xyz "$OUT/$c/"
    [ -d "$d/plot" ] && { rm -rf "$OUT/$c/plot"; cp -r "$d/plot" "$OUT/$c/plot"; }
done
cp "$OUT/round2/neb_profile.png" figures/round2_profile.png
echo "done -- python3 results_figures.py draws the summary figures"
