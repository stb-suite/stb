#!/bin/bash
# Regenerates reference_results/: the elastic-constant calculations of README section 5 -- bulk silicon with the
# fine-tuned model of the 7.2 example (../7.2-mlmd/model/, trained by the 7.1 workflow on SIESTA PBE/DZP data), the
# raw MACE-MP-0 foundation models, and two non-cubic cases with the foundation model (hcp magnesium, graphene).
# About 15 minutes on one RTX 3060 (DEVICE=cuda, the default here); DEVICE=cpu works too (similar time for the
# 8-atom cells, longer for the 216-atom one).
#
#   bash make_reference.sh [work_dir]            run every case, then copy the results into reference_results/
#   bash make_reference.sh [work_dir] collect    only copy the finished cases
#
# work_dir defaults to ./reference_run. Each case is skipped when its report already exists. Reports, data files
# and gnuplot scripts are copied into reference_results/; the PNG figures stay in the work folder.
set -e
cd "$(dirname "${BASH_SOURCE[0]}")"
HERE="$PWD"
WORK="$(mkdir -p "${1:-reference_run}" && cd "${1:-reference_run}" && pwd)"
STEP="${2:-all}"
export MPLBACKEND=Agg
DEVICE="${DEVICE:-cuda}"

MODEL="$HERE/../7.2-mlmd/model/mlff_model_compiled.model"
CARD="$HERE/../7.2-mlmd/model/mlff_model_card.json"
FT="--custom-model $MODEL --custom-model-card $CARD"

# case <name> <structure> <stb-mlelastic flags...>
case_() {
    local name=$1 file=$2; shift 2
    [ -s "$WORK/$name/stb_mlelastic_report.txt" ] && { echo "$name: done"; return; }
    case "$file" in /*) ;; *) file="$HERE/$file" ;; esac
    mkdir -p "$WORK/$name"
    echo "$name: running ($(date +%H:%M:%S))"
    (cd "$WORK/$name" && stb-mlelastic --file "$file" --device "$DEVICE" -o . --save-report --save-data \
        --no-intro "$@" < /dev/null > log.txt 2>&1) || echo "$name: FAILED (see $WORK/$name/log.txt)"
}

run_all() {
    # 1. the reference calculation: 7 strains in +-1 %, fine-tuned model against the medium foundation model
    case_ base si8.fdf $FT --steps 7 --save-gnuplot
    # 2. strain amplitude, 0.25 .. 3 %
    case_ conv si8.fdf $FT --steps 7 --check-convergence --convergence-strains 0.25 0.5 1 2 3
    # 3. the force threshold of the internal relaxation
    for f in 0.02 0.005 0.001 0.0002; do
        case_ fmax_$f si8.fdf $FT --steps 5 --fmax $f --skip-foundation-comparison --no-clamped
    done
    # 4. the three foundation models alone
    for m in small medium large; do
        case_ foundation_$m si8.fdf --model $m --steps 7
    done
    # 5. a 216-atom supercell (stb-supercell -d 3 3 3): the constants do not depend on the cell
    if [ ! -s "$WORK/si216.fdf" ]; then
        (cd "$WORK" && stb-supercell -f "$HERE/si8.fdf" -d 3 3 3 -o si216.fdf --no-intro < /dev/null > /dev/null 2>&1)
    fi
    case_ si216 "$WORK/si216.fdf" $FT --steps 5 --skip-foundation-comparison --no-clamped
    # 6. a structure that is not at a minimum: without and with the pre-relax
    case_ rattled_noprerelax si8_rattled.fdf $FT --steps 5 --skip-foundation-comparison --no-clamped --no-pre-relax
    case_ rattled si8_rattled.fdf $FT --steps 5 --skip-foundation-comparison --no-clamped
    # 7. hexagonal: hcp magnesium, foundation model
    case_ mg_hcp mg_hcp.fdf --model medium --steps 7 --save-gnuplot
    # 8. 2D: graphene, foundation model, strains up to 0.5 %
    case_ graphene graphene.fdf --model medium --steps 7 --max 0.5
}

collect() {
    for d in "$WORK"/*/; do
        name=$(basename "$d")
        [ -s "$d/stb_mlelastic_report.txt" ] || continue
        out="$HERE/reference_results/$name"
        mkdir -p "$out"
        cp "$d"/stb_mlelastic_report.txt "$d"/*.dat "$out/" 2>/dev/null || true
        if [ -d "$d/plot" ]; then
            mkdir -p "$out/plot"
            cp "$d"/plot/*.gplot "$d"/plot/*.dat "$out/plot/" 2>/dev/null || true
        fi
    done
    # the tool's own figures shown in the README
    mkdir -p "$HERE/figures"
    for pair in base/stress_strain.png:si_stress_strain.png base/directional_moduli.png:si_directional.png \
                base/anisotropy_surface.png:si_anisotropy_surface.png mg_hcp/directional_moduli.png:mg_directional.png \
                graphene/inplane_youngs.png:graphene_inplane.png conv/strain_convergence.png:si_strain_convergence.png; do
        src="$WORK/${pair%%:*}"; dst="$HERE/figures/${pair##*:}"
        [ -f "$src" ] && cp "$src" "$dst"
    done
    echo "reference_results/ written"
}

case "$STEP" in
    collect) collect ;;
    *) run_all; collect ;;
esac
