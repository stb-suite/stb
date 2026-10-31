#!/bin/bash
# =====================================================================================================
#  make_reference.sh -- regenerates reference_results/<case>/ of the 5.11 example (Raman workflow)
#
#  Every case runs on a copy of runs/ (finished SIESTA 5.4 calculations, made by run_example.sh) in a scratch
#  folder: the three stages' reports, and the spectra. Dates become <date>, the scratch path <work>. Seconds, no
#  SIESTA.
#
#      source ~/.anaconda3/etc/profile.d/conda.sh && conda activate stb    # any environment with stb_suite
#      bash make_reference.sh
# =====================================================================================================
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
REF="$HERE/reference_results"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
clean() { sed 's/\x1b\[[0-9;]*m//g' | sed "s#$WORK#<work>#g" | sed -E 's/^(Date\/time *: ).*/\1<date>/'; }
copy_run() { rm -rf "$2"; mkdir -p "$2"; cp -a "$HERE/runs/$1/." "$2/"; find "$2" -name "siesta.EPSIMG.gz" -exec gunzip -f {} \; ; }
report() {       # report <case> <file> <dir> <command...>: runs in <dir>, stores the plain-text report
    local case="$1" file="$2" dir="$3"; shift 3
    mkdir -p "$REF/$case"
    echo "  [run]  $case/$file: $*"
    (cd "$dir" && "$@" --no-intro < /dev/null 2>&1) | clean > "$REF/$case/$file" || true
}
SI2=(--kgrid 11 11 11 --full-tensor --use-symmetry --skip-degenerate --optical-mesh 32 32 32 --optical-broaden 0.1
     --optical-nbands 16)

echo "== 1. Stage 1 on bulk Si"
mkdir -p "$WORK/stage1"; cp "$HERE/inputs/si/structure.fdf" "$HERE/inputs/si/calc.fdf" "$WORK/stage1/"
report si stage1.txt "$WORK/stage1" stb-raman -s structure.fdf -c calc.fdf -p dojo -dim 1 1 1 --kgrid 11 11 11

echo "== 2. Stage 2 on the finished phonons (and without --kgrid)"
mkdir -p "$WORK/stage2/raman_study"; cp "$HERE/inputs/si/calc.fdf" "$WORK/stage2/"
cp -a "$HERE/runs/si/raman_study/phonon_disp" "$WORK/stage2/raman_study/"
report si stage2.txt "$WORK/stage2" stb-ramanModes --directory raman_study -c calc.fdf "${SI2[@]}"
rm -rf "$WORK/stage2/raman_study/optical_disp"
report si stage2_nokgrid.txt "$WORK/stage2" stb-ramanModes --directory raman_study -c calc.fdf -p dojo --full-tensor \
    --use-symmetry --skip-degenerate

echo "== 3. Stage 3 on the finished Si runs"
copy_run si "$WORK/si"
report si stage3.txt "$WORK/si" stb-ramanAnalysis --directory raman_study --save-gnuplot
cp "$WORK/si/raman_study/plot/raman_spectrum.dat" "$REF/si/"

echo "== 4. The diagonal-only trap (the same 6 x/y/z runs, analysed alone)"
mkdir -p "$WORK/diag/raman_study"; cp "$HERE/inputs/si/calc.fdf" "$WORK/diag/"
cp -a "$HERE/runs/si/raman_study/phonon_disp" "$WORK/diag/raman_study/"
(cd "$WORK/diag" && stb-ramanModes --directory raman_study -c calc.fdf -p dojo --kgrid 11 11 11 --use-symmetry \
    --skip-degenerate --optical-mesh 32 32 32 --optical-broaden 0.1 --optical-nbands 16 --no-intro > /dev/null 2>&1)
for f in "$WORK"/diag/raman_study/optical_disp/mode_01_*; do
    cp "$WORK/si/raman_study/optical_disp/$(basename "$f")/calc.out" \
       "$WORK/si/raman_study/optical_disp/$(basename "$f")/siesta.EPSIMG" "$f/"
done
report si_diagonal stage3.txt "$WORK/diag" stb-ramanAnalysis --directory raman_study

echo "== 5. Convergence: displacement and Optical.Mesh"
for s in si_d001 si_d004 si_mesh16; do
    copy_run "$s" "$WORK/$s"
    report "$s" stage3.txt "$WORK/$s" stb-ramanAnalysis --directory raman_study
done

echo "== 6. Graphene"
copy_run graphene "$WORK/graphene"
mkdir -p "$REF/graphene"; cp "$WORK/graphene/raman_study/raman_stage2.txt" "$REF/graphene/stage2.txt"
report graphene stage3.txt "$WORK/graphene" stb-ramanAnalysis --directory raman_study --save-gnuplot
cp "$WORK/graphene/raman_study/plot/raman_spectrum.dat" "$REF/graphene/"

echo "== 7. The interactive menu (stb-suite, code 5.11.3)"
mkdir -p "$REF/menu"
(cd "$WORK/si" && printf '5.11.3\nraman_study\n\n\n\n\nn\nn\n\n0\n' | stb-suite 2>&1) | clean \
    | sed -n '/RAMAN STAGE 3 REPORT/,/Raman analysis complete/p' > "$REF/menu/report.txt"
echo "done -- reference_results/<case>/"
