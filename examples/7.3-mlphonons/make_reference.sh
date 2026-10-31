#!/bin/bash
# Regenerates reference_results/: the phonon calculations of README section 5 -- bulk silicon with the fine-tuned
# model of the 7.2 example (../7.2-mlmd/model/, trained by the 7.1 workflow on SIESTA PBE/DZP data) and the raw
# MACE-MP-0 foundation models for comparison. Everything runs on the CPU: about 10 minutes in total on a 12-core
# desktop (one single-volume 2x2x2 calculation takes ~15 s, the 4x4x4 one ~30 s, the 7-volume QHA ~40 s).
#
#   bash make_reference.sh [work_dir]            run every case, then copy the results into reference_results/
#   bash make_reference.sh [work_dir] collect    only copy the finished cases
#
# work_dir defaults to ./reference_run. Each case is skipped when its report already exists, so the script can be
# re-run after an interruption. Only reports, data files, gnuplot scripts and the small force-constant files are
# copied into reference_results/; the PNG figures stay in the work folder.
set -e
cd "$(dirname "${BASH_SOURCE[0]}")"
HERE="$PWD"
WORK="$(mkdir -p "${1:-reference_run}" && cd "${1:-reference_run}" && pwd)"
STEP="${2:-all}"
export MPLBACKEND=Agg
DEVICE="${DEVICE:-cpu}"

MODEL="$HERE/../7.2-mlmd/model/mlff_model_compiled.model"
CARD="$HERE/../7.2-mlmd/model/mlff_model_card.json"
FT="--custom-model $MODEL --custom-model-card $CARD"

# case <name> <structure> <stb-mlphonons flags...>
case_() {
    local name=$1 file=$2; shift 2
    [ -s "$WORK/$name/stb_mlphonons_report.txt" ] && { echo "$name: done"; return; }
    case "$file" in /*) ;; *) file="$HERE/$file" ;; esac
    mkdir -p "$WORK/$name"
    echo "$name: running ($(date +%H:%M:%S))"
    (cd "$WORK/$name" && stb-mlphonons --file "$file" --device "$DEVICE" -o . --save-report --save-data \
        --no-intro "$@" < /dev/null > log.txt 2>&1) || echo "$name: FAILED (see $WORK/$name/log.txt)"
}

run_all() {
    # 1. the reference calculation: 2x2x2 supercell, fine-tuned model against the medium foundation model,
    #    the first optical mode at Gamma animated, gnuplot files
    case_ base si8.fdf $FT --animate-mode --save-gnuplot
    # 2. the same in a 4x4x4 supercell (512 atoms): converged sound velocities and bands
    case_ base_444 si8.fdf $FT --dim 4 4 4 --skip-foundation-comparison
    # 3. the three foundation models alone
    for m in small medium large; do
        case_ foundation_$m si8.fdf --model $m
    done
    # 4. supercell convergence, 1x1x1 .. 4x4x4
    case_ conv si8.fdf $FT --check-convergence --convergence-dims 1 1 1 2 2 2 3 3 3 4 4 4 --save-gnuplot
    # 5. quasi-harmonic approximation, 7 volumes over +-4 %, up to 1500 K
    case_ qha si8.fdf $FT --qha --n-volumes 7 --strain-range 4 --tmax 1500 --save-gnuplot
    #    the foundation model: first on the same volumes (centred on the fine-tuned model's lattice), then centred on
    #    its own static minimum -- its cell relaxed with stb-mlrelax (a = 5.4554 Ang instead of 5.4802)
    case_ qha_foundation si8.fdf --model medium --qha --n-volumes 7 --strain-range 4 --tmax 1500
    if [ ! -s "$WORK/si8_mp0.fdf" ]; then
        (cd "$WORK" && stb-mlrelax -f "$HERE/si8.fdf" --relax-cell --model medium --fmax 0.001 --device "$DEVICE" \
            -o si8_mp0.fdf --no-intro < /dev/null > mlrelax_mp0.log 2>&1) || echo "stb-mlrelax: FAILED"
    fi
    case_ qha_foundation_v0 "$WORK/si8_mp0.fdf" --model medium --qha --n-volumes 7 --strain-range 4 --tmax 1500
    # 6. displacement amplitude
    for d in 0.001 0.003 0.01 0.03 0.1 0.2; do
        case_ dist_$d si8.fdf $FT --distance $d --skip-foundation-comparison
    done
    # 7. a structure that is not at a minimum: no pre-relax, the default pre-relax, a tight pre-relax
    case_ rattled_norelax si8_rattled.fdf $FT --no-relax --skip-foundation-comparison --freeze-unstable-mode
    case_ rattled_relax   si8_rattled.fdf $FT --skip-foundation-comparison
    case_ rattled_tight   si8_rattled.fdf $FT --skip-foundation-comparison --fmax 0.001
    # 8. elastic constants from stress-strain with the same model: a cross-check of the sound velocities
    if [ ! -s "$WORK/elastic/stb_mlelastic_report.txt" ]; then
        mkdir -p "$WORK/elastic"
        echo "elastic: running ($(date +%H:%M:%S))"
        (cd "$WORK/elastic" && stb-mlelastic --file "$HERE/si8.fdf" $FT --device "$DEVICE" -o . --save-report \
            --no-intro < /dev/null > log.txt 2>&1) || echo "elastic: FAILED (see $WORK/elastic/log.txt)"
    else
        echo "elastic: done"
    fi
    # 9. the SIESTA phonons of the 7.1 example against the model, and the force-noise experiment
    if [ ! -s "$HERE/reference_results/siesta_compare.json" ]; then
        echo "siesta_compare: running ($(date +%H:%M:%S))"
        python3 "$HERE/siesta_compare.py" --device "$DEVICE" > "$WORK/siesta_compare.log" 2>&1 \
            || echo "siesta_compare: FAILED (see $WORK/siesta_compare.log)"
    else
        echo "siesta_compare: done"
    fi
}

collect() {
    for d in "$WORK"/*/; do
        name=$(basename "$d")
        out="$HERE/reference_results/$name"
        if [ -s "$d/stb_mlphonons_report.txt" ]; then
            mkdir -p "$out"
            cp "$d"/stb_mlphonons_report.txt "$d"/*.dat "$out/" 2>/dev/null || true
            [ -f "$d/phonopy_disp.yaml" ] && [ "$name" = base ] && cp "$d/phonopy_disp.yaml" "$out/"
            cp "$d"/mode_*.xsf "$out/" 2>/dev/null || true
            if [ -d "$d/plot" ]; then
                mkdir -p "$out/plot"
                cp "$d"/plot/*.gplot "$d"/plot/*.dat "$out/plot/" 2>/dev/null || true
            fi
        elif [ -s "$d/stb_mlelastic_report.txt" ]; then
            mkdir -p "$out"
            cp "$d"/stb_mlelastic_report.txt "$out/"
        fi
    done
    echo "reference_results/ written"
}

case "$STEP" in
    collect) collect ;;
    *) run_all; collect ;;
esac
