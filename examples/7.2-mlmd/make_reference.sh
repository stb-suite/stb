#!/bin/bash
# Regenerates reference_results/: the molecular-dynamics runs of README section 5 -- bulk silicon with the
# fine-tuned model in model/ (and the raw MACE-MP medium foundation model for comparison), on a GPU.
# About 4 h on one RTX 3060 for both lanes run one after the other (~2.5 h with the two lanes side by
# side; 216 atoms cost ~75 ms/step in NVT/NVE and ~140 ms/step in NPT). Trajectories stay in the work
# folder; only reports, data files, gnuplot scripts and dashboards are copied into reference_results/.
#
#   bash make_reference.sh [work_dir] [lane]     lane = A, B or all (default all); work_dir default ./reference_run
#   bash make_reference.sh [work_dir] collect    copy the finished runs into reference_results/
#   bash make_reference.sh [work_dir] siesta     SIESTA single points on the final frames (section 5.3; ~2 h each
#                                                on 6 cores, mpirun -np ${NP:-6}), then siesta_check.py
#
# Each case is skipped when its report already exists, so the script can be re-run after an interruption.
set -e
cd "$(dirname "${BASH_SOURCE[0]}")"
HERE="$PWD"
WORK="$(mkdir -p "${1:-reference_run}" && cd "${1:-reference_run}" && pwd)"
LANE="${2:-all}"
export OMP_NUM_THREADS=2 MPLBACKEND=Agg

FT="--custom-model $HERE/model/mlff_model_compiled.model --custom-model-card $HERE/model/mlff_model_card.json"
FOUND="--model medium"
SI_BETA=1.0e-6          # 1/bar: 1/(B0 * 1e4) with B0 ~ 98 GPa, silicon's own compressibility

structures() {
    cd "$WORK"
    [ -f si216.fdf ] && return
    cp "$HERE/si8.fdf" .
    stb-supercell -f si8.fdf -d 2 2 2 -o si64.fdf  --no-intro < /dev/null > /dev/null
    stb-supercell -f si8.fdf -d 3 3 3 -o si216.fdf --no-intro < /dev/null > /dev/null
    stb-supercell -f si8.fdf -d 4 4 4 -o si512.fdf --no-intro < /dev/null > /dev/null
}

# case <name> <structure> <stb-mlmd flags...>
case_() {
    local name=$1 file=$2; shift 2
    [ -s "$WORK/$name/stb_mlmd_report.txt" ] && { echo "$name: done"; return; }
    mkdir -p "$WORK/$name"
    echo "$name: running ($(date +%H:%M))"
    (cd "$WORK/$name" && stb-mlmd --file "$WORK/$file" --device cuda --save-report --save-gnuplot --png \
        -o traj.xsf --no-intro "$@" < /dev/null > mlmd.log 2>&1) || echo "$name: FAILED (see $WORK/$name/mlmd.log)"
}

lane_A() {
    # 1. the reference run: NVT at 300 K, 10 ps
    case_ nvt300 si216.fdf $FT --ensemble nvt --temperature 300 --equilibration-steps 1000 --n-steps 10000 --seed 1
    # 2. energy conservation vs timestep: NVE from the equilibrated structure, 5 ps each
    for dt_n in 0.5:10000 1:5000 2:2500 3:1667 4:1250 5:1000; do
        dt=${dt_n%%:*}; n=${dt_n##*:}
        case_ nve_dt$dt nvt300/md_final.fdf $FT --ensemble nve --temperature 300 --timestep $dt --n-steps $n \
            --stride $(python -c "print(max(1, round(10/$dt)))") --seed 2
    done
    # 3. equipartition: the perfect crystal given 600 K of kinetic energy, NVE
    case_ nve_equipartition si216.fdf $FT --ensemble nve --temperature 600 --n-steps 2000 --stride 2 --seed 3
    # 4. finite size: NVT at 300 K for 8, 64 and 512 atoms
    case_ size_N8   si8.fdf   $FT --ensemble nvt --temperature 300 --equilibration-steps 1000 --n-steps 20000 --seed 4
    case_ size_N64  si64.fdf  $FT --ensemble nvt --temperature 300 --equilibration-steps 1000 --n-steps 10000 --seed 4
    case_ size_N512 si512.fdf $FT --ensemble nvt --temperature 300 --equilibration-steps 1000 --n-steps 5000 --seed 4
    # 5. thermal expansion, the hottest point
    case_ npt_T1500 si216.fdf $FT --ensemble npt --temperature 1500 --compressibility $SI_BETA \
        --equilibration-steps 2000 --n-steps 6000 --seed 5
    # 6. far outside the training data: 2500 K at 0 GPa (silicon melts at 1687 K)
    case_ melt_ft si216.fdf $FT --ensemble npt --temperature 2500 --compressibility $SI_BETA \
        --equilibration-steps 3000 --n-steps 6000 --seed 6
}

lane_B() {
    # 1b. the same NVT run with the raw foundation model
    case_ nvt300_foundation si216.fdf $FOUND --ensemble nvt --temperature 300 --equilibration-steps 1000 --n-steps 10000 --seed 1
    # 5. thermal expansion: NPT at 0 GPa, silicon's own compressibility
    for T in 300 600 900 1200; do
        case_ npt_T$T si216.fdf $FT --ensemble npt --temperature $T --compressibility $SI_BETA \
            --equilibration-steps 2000 --n-steps 6000 --seed 5
    done
    # 5b. the barostat with the placeholder compressibility (water's, 45x silicon's)
    case_ npt_T300_waterbeta si216.fdf $FT --ensemble npt --temperature 300 --n-steps 4000 --stride 5 --seed 5
    case_ npt_T300_sibeta    si216.fdf $FT --ensemble npt --temperature 300 --compressibility $SI_BETA \
        --n-steps 4000 --stride 5 --seed 5
    # 6b. the melt with the foundation model
    case_ melt_foundation si216.fdf $FOUND --ensemble npt --temperature 2500 --compressibility $SI_BETA \
        --equilibration-steps 3000 --n-steps 6000 --seed 6
}

collect() {
    for d in "$WORK"/*/; do
        name=$(basename "$d"); [ -s "$d/stb_mlmd_report.txt" ] || continue
        out="$HERE/reference_results/$name"; mkdir -p "$out"
        cp "$d"/stb_mlmd_report.txt "$d"/*_md_thermo.dat "$d"/*_md_rdf.dat "$d"/*_md_msd.dat \
           "$d"/*_md.gplot "$d"/md_final.fdf "$out/" 2>/dev/null || true
    done
    echo "reference_results/ written"
}

# SIESTA single point on the final frame of a finished run: <frame name> <case>
siesta_frame() {
    local frame=$1 case=$2 d="$WORK/siesta_check/$1"
    [ -s "$d/siesta.FA" ] && { echo "siesta $frame: done"; return; }
    [ -s "$WORK/$case/md_final.fdf" ] || { echo "siesta $frame: $case not finished"; return; }
    mkdir -p "$d"
    cp "$WORK/$case/md_final.fdf" "$d/structure.fdf"
    cp "$HERE/siesta_calc.fdf" "$d/calc.fdf"
    cp "$HERE/siesta_config_extra.fdf" "$d/config_extra.fdf"
    python3 -c "from stb.core.pseudopotentials import resolve_pseudo_source as r, copy_pseudo as c; c(r('dojo'), 'Si', '$d')" > /dev/null
    echo "siesta $frame: running ($(date +%H:%M))"
    (cd "$d" && mpirun -np "${NP:-6}" siesta calc.fdf < /dev/null > siesta.out 2> siesta.err)
}

siesta_checks() {
    siesta_frame nvt300_final nvt300
    siesta_frame npt1500_final npt_T1500
    siesta_frame melt_ft_final melt_ft
    for d in "$WORK"/siesta_check/*/; do
        [ -s "$d/siesta.FA" ] || continue
        out="$HERE/reference_results/siesta_check/$(basename "$d")"; mkdir -p "$out"
        cp "$d"/structure.fdf "$d"/calc.fdf "$d"/config_extra.fdf "$d"/siesta.out "$d"/siesta.FA "$out/"
    done
    python3 "$HERE/siesta_check.py" --device cuda
}

structures
case "$LANE" in
    siesta) siesta_checks ;;
    A) lane_A ;;
    B) lane_B ;;
    collect) collect ;;
    *) lane_A; lane_B; collect ;;
esac
