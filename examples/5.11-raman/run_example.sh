#!/bin/bash
# Runs every SIESTA calculation of the 5.11 example (runs/), stage by stage, with the suite's own commands:
#   runs/si/               bulk Si (inputs/si): Stage 1 (one displacement + the undisplaced cell), Stage 2 with the full
#                          tensor for one T2g partner, displacement 0.02 Ang, Optical.Mesh 32^3   -- the main case
#   runs/si_d001/, si_d004/  the same with displacement 0.01 and 0.04 Ang (linearity of the finite difference)
#   runs/si_mesh16/        the same with Optical.Mesh 16^3 (convergence of eps1(0) and of the tensor)
#   runs/graphene/         the graphene monolayer of inputs/graphene: Stage 1 in a 2x2x1 supercell, Stage 2 with the
#                          full tensor for ALL three optical modes (B2g and both E2g partners), Optical.Mesh 60x60x1
# SIESTA 5.x + mpirun. On 2 cores: about 50 min for silicon, 25 min for graphene. Pseudopotentials from the suite's
# banks (PseudoDojo for Si, virtual_vault for C), copied in by the tools and removed afterwards. Each SIESTA folder
# keeps only its inputs, calc.out, the forces (.FA) or the dielectric function (siesta.EPSIMG, gzipped): the
# 1.4 GB siesta.contrib.EPSIMG of every optical run is deleted as soon as it finishes.
#     bash run_example.sh [NPROC]
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
R="$HERE/runs"
NP=${1:-2}

siesta_in() {     # siesta_in <folder>...: runs SIESTA in every folder, then keeps only what the stages read
    for d in "$@"; do
        [ -f "$d/calc.fdf" ] || continue
        echo "[siesta] ${d#$R/}"
        (cd "$d" && OMP_NUM_THREADS=1 mpirun -np "$NP" siesta < calc.fdf > calc.out 2> siesta.err)
        rm -f "$d/siesta.contrib.EPSIMG"
    done
}

slim() {          # slim <study>: drops everything SIESTA wrote that the stages do not read
    find "$1" -mindepth 2 -type f \( -name "*.psml" -o -name "*.psf" -o -name "*.DM" -o -name "*.HSX" -o -name "*.ion*" \
        -o -name "*.EIG" -o -name "*.KP" -o -name "*.XV" -o -name "*.STRUCT_OUT" -o -name "*.BONDS*" -o -name "*.ORB_INDX" \
        -o -name "*.alloc" -o -name "*.times" -o -name "*.bib" -o -name "*.xyz" -o -name "fdf*.log" -o -name "*ENTHALPY" \
        -o -name "TIMES" -o -name "CLOCK" -o -name "MESSAGES" -o -name "OUTVARS.yml" -o -name "PARALLEL_DIST" \
        -o -name "NON_TRIMMED_KP_LIST" -o -name "FORCE_STRESS" -o -name "0_NORMAL_EXIT" -o -name "siesta.err" \
        -o -name "*.nc" \) -delete
    find "$1" -name "siesta.EPSIMG" -exec gzip -9 {} \;
}

study() {         # study <name> <inputs> <stage-1 args> -- <stage-2 args>
    local name=$1 inputs=$2; shift 2
    local s1=() s2=()
    while [ "$1" != "--" ]; do s1+=("$1"); shift; done; shift
    s2=("$@")
    rm -rf "$R/$name"; mkdir -p "$R/$name"
    cp "$HERE/inputs/$inputs/structure.fdf" "$HERE/inputs/$inputs/calc.fdf" "$R/$name/"
    (cd "$R/$name" && stb-raman -s structure.fdf -c calc.fdf "${s1[@]}" --no-intro > /dev/null)
    siesta_in "$R/$name"/raman_study/phonon_disp/disp-* "$R/$name/raman_study/phonon_disp/perfect_supercell"
    (cd "$R/$name" && stb-ramanModes --directory raman_study -c calc.fdf "${s2[@]}" --no-intro > /dev/null)
    siesta_in "$R/$name"/raman_study/optical_disp/mode_*_*
    slim "$R/$name"
}

SI_S1=(-p dojo -dim 1 1 1 --kgrid 11 11 11)
SI_S2=(--kgrid 11 11 11 --full-tensor --use-symmetry --skip-degenerate --optical-broaden 0.1 --optical-nbands 16)
study si        si "${SI_S1[@]}" -- "${SI_S2[@]}" --optical-mesh 32 32 32 --export-animations
study si_mesh16 si "${SI_S1[@]}" -- "${SI_S2[@]}" --optical-mesh 16 16 16
study si_d001   si "${SI_S1[@]}" -- "${SI_S2[@]}" --optical-mesh 32 32 32 --displacement 0.01
study si_d004   si "${SI_S1[@]}" -- "${SI_S2[@]}" --optical-mesh 32 32 32 --displacement 0.04
study graphene  graphene -p virtual_vault -dim 2 2 1 -- --full-tensor --optical-mesh 60 60 1 --optical-broaden 0.1 \
    --export-animations
echo "done -- runs/"
