#!/bin/bash
# Guided example: 7.2 -- ML Molecular Dynamics (stb-mlmd; code 7.2 in the stb-suite menu).
#
# Not an automated test (see test/7-mlsimulations/2-mlmd/test.sh for that) -- a commented walk-through: it runs
# real commands, one case at a time, into its own output/<case>/ folder, and shows the piece of output that proves
# what just happened. Pauses between sections so you can read before moving on. Safe to re-run any time -- it
# always starts by wiping its own output/.
#
# Every MD here is real and runs on the CPU in about 10 minutes in total (8- and 64-atom silicon cells, short
# runs), with the fine-tuned model shipped in model/ (trained by the 7.1 workflow on SIESTA PBE/DZP data of bulk
# silicon) and, in Case 7, the MACE-MP-0 medium foundation model (downloaded once on first use). Set DEVICE=cuda
# to run on a GPU. The long runs of README section 5 (216-512 atoms, 5-10 ps, on a GPU) are NOT rerun here: Case 10
# reads their reports from reference_results/ (make_reference.sh regenerates them, ~4 h on one GPU).
#
# Needs the 'ml' extra: pip install stb_suite[ml]

set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

export MPLBACKEND=Agg           # headless-safe; the dashboards are still written with --png
DEVICE="${DEVICE:-cpu}"
MODEL="$DIR/model/mlff_model_compiled.model"
CARD="$DIR/model/mlff_model_card.json"
FT="--custom-model $MODEL --custom-model-card $CARD"

OUT="$DIR/output"
rm -rf "$OUT"
mkdir -p "$OUT"

pause() {
    echo
    read -p "  [Press Enter to continue] " -r
    echo
}

banner() {
    echo "=================================================================="
    echo " $1"
    echo "=================================================================="
}

strip_colors() { sed 's/\x1b\[[0-9;]*m//g'; }

# section N of a report: from "[N] TITLE" up to the next "[M] " header
section() {   # $1 = file, $2 = section number
    strip_colors < "$1" | awk -v n="[$2]" 'index($0, n) == 1 {on = 1; print; next}
        on && /^\[[0-9]+\] / {exit} on {print}'
}

# mlmd <case folder> <structure> <flags...>: runs stb-mlmd there, report saved, log kept
mlmd() {
    local case=$1 file=$2; shift 2
    mkdir -p "$OUT/$case"
    (cd "$OUT/$case" && stb-mlmd --file "$file" --device "$DEVICE" --save-report --no-intro "$@" < /dev/null > log.txt 2>&1)
}


banner "Welcome: molecular dynamics with a machine-learned potential"
cat <<'EOF'
stb-mlmd integrates Newton's equations for every atom with the forces of a MACE potential -- the MACE-MP-0
foundation model, or a model you fine-tuned on YOUR SIESTA calculations with the 7.1 workflow, which always comes
with its model card. NVE (constant energy), NVT (Langevin thermostat) or NPT (Berendsen thermostat + barostat).

The report says what the trajectory means -- averages with error bars, the ensemble checks, the structure (g(r),
coordination, MSD, Lindemann ratio) -- and whether to trust it: the model card in full, and a check of the run
against the range of structures the model was trained on. Every warning goes to the last section.

This tour (bulk silicon, CPU): a first NVT run and its report (1), the optional outputs (2), equipartition (3),
energy conservation vs timestep (4), NPT and the barostat's compressibility (5), a structure outside the trained
range (6), the foundation model (7), the guard rails (8), the interactive menu (9), and the real GPU runs of the
README (10). Theory: README section 1. Worked example with the numbers: README section 5.
EOF
pause


banner "Case 1: a first NVT run -- 8 Si atoms at 300 K with the fine-tuned model and its card"
cat <<'EOF'
si8.fdf is the conventional diamond cell at SIESTA's relaxed lattice constant (5.48024 Ang). The model comes with
its card: --custom-model AND --custom-model-card, always both. 200 equilibration steps are discarded, then 1500
production steps of 1 fs. Report section [2] prints the card -- here its first blocks: identity (SHA-256 checked
against the file), level of theory, trained range.
EOF
echo "\$ stb-mlmd --file si8.fdf $FT --ensemble nvt --temperature 300 --equilibration-steps 200 --n-steps 1500 --seed 1 --png --save-gnuplot --save-report"
mlmd case1_nvt300 "$DIR/si8.fdf" $FT --ensemble nvt --temperature 300 --equilibration-steps 200 --n-steps 1500 \
    --seed 1 --png --save-gnuplot
R1="$OUT/case1_nvt300/stb_mlmd_report.txt"
section "$R1" 2 | sed -n '1,40p'
echo "   ... (recipe and accuracy, stability and usage follow in the report)"
pause
cat <<'EOF'
Sections [6]-[7]: averages with block-average errors (+/-), the heat capacity from the energy fluctuations, and
the checks. With 8 atoms the temperature fluctuates by ~90 K at any instant: sqrt(2/N_dof) = sqrt(2/21) = 31 %.
EOF
section "$R1" 6
section "$R1" 7
pause
cat <<'EOF'
Section [8]: the structure along the run (g(r): 4 neighbours at 2.37 Ang, a solid-like Lindemann ratio) and the
domain table -- every range of the card's training data against this run. [10] collects every warning.
EOF
section "$R1" 8
section "$R1" 10
pause


banner "Case 2: the optional outputs -- dashboard, data files, gnuplot"
cat <<'EOF'
Nothing but the report, the trajectory and md_final.fdf is written unless asked. Case 1 asked for --png (six-panel
dashboard) and --save-gnuplot (three .dat files + a .gplot script that draws the same panels into a PDF).
EOF
ls -1 "$OUT/case1_nvt300"
echo
echo "--- head of si8_md_thermo.dat (one line per saved frame, energies per atom) ---"
head -4 "$OUT/case1_nvt300/si8_md_thermo.dat"
if command -v gnuplot > /dev/null; then
    (cd "$OUT/case1_nvt300" && gnuplot si8_md.gplot 2> /dev/null) && echo && echo "[OK] gnuplot si8_md.gplot -> $(ls "$OUT/case1_nvt300"/si8_md.pdf)"
else
    echo "(gnuplot not installed: 'gnuplot si8_md.gplot' would draw si8_md.pdf)"
fi
pause


banner "Case 3: equipartition -- a perfect crystal given 600 K of kinetic energy, no thermostat"
cat <<'EOF'
--supercell 2 2 2 repeats si8.fdf into the 64-atom cell before the run (the same cell, and the same atom order, as
stb-supercell -d 2 2 2 would write; report [1] says 8 -> 64 atoms). NVE from the perfect lattice with velocities
drawn at 600 K: all the energy starts as kinetic, and within a few vibrational periods (~65 fs) equipartition
shares it with the potential energy -- the temperature drops to about half, 300 K. This is the transient
--equilibration-steps throws away.
EOF
echo "\$ stb-mlmd --file si8.fdf --supercell 2 2 2 $FT --ensemble nve --temperature 600 --n-steps 300 --stride 3 --seed 3"
mlmd case3_equipartition "$DIR/si8.fdf" $FT --supercell 2 2 2 --ensemble nve --temperature 600 --n-steps 300 \
    --stride 3 --seed 3 --save-data
strip_colors < "$OUT/case3_equipartition/stb_mlmd_report.txt" | grep -E "^Supercell|^\[OK\] Supercell|^Degrees of freedom" | sed 's/^/  /'
echo
section "$OUT/case3_equipartition/stb_mlmd_report.txt" 5 | sed -n '1,16p'
echo
echo "--- T(t) over the first 150 fs (from si8_md_thermo.dat) ---"
awk '!/^#/ && $1 <= 150 {printf "  t = %5.0f fs   T = %6.1f K\n", $1, $2}' "$OUT/case3_equipartition/si8_md_thermo.dat" | awk 'NR % 5 == 1'
pause


banner "Case 4: energy conservation vs timestep (NVE, 1 ps each)"
cat <<'EOF'
Velocity Verlet conserves a "shadow" energy: the total energy oscillates by ~(omega dt)^2/4 but does not drift as
long as dt is well below the stability limit (omega dt = 2, ~20 fs for silicon's 15.5 THz optical phonon). Three
NVE runs of 1 ps from the equilibrated structure of Case 1, at 1, 3 and 6 fs. Watch the drift (meV/atom/ps) and the
total-energy fluctuation relative to the potential-energy one.
EOF
printf "  %-6s %-12s %-62s\n" "dt" "steps" "report [7]"
for dt_n in 1:1000 3:333 6:167; do
    dt=${dt_n%%:*}; n=${dt_n##*:}
    mlmd case4_nve_dt$dt "$OUT/case1_nvt300/md_final.fdf" $FT --ensemble nve --temperature 300 --timestep $dt \
        --n-steps $n --stride 1 --seed 2
    r="$OUT/case4_nve_dt$dt/stb_mlmd_report.txt"
    drift=$(strip_colors < "$r" | grep "^Energy drift" | sed 's/Energy drift *: //')
    fluc=$(strip_colors < "$r" | grep "^E_tot fluctuation" | sed 's/E_tot fluctuation *: //; s/ -- .*//')
    printf "  %-6s %-12s drift %s\n  %-6s %-12s fluct %s\n" "$dt fs" "$n" "$drift" "" "" "$fluc"
done
echo
echo "[10] of the 6 fs run:"
section "$OUT/case4_nve_dt6/stb_mlmd_report.txt" 10 | sed -n '3,6p'
pause


banner "Case 5: NPT -- the barostat and its compressibility"
cat <<'EOF'
Berendsen's barostat moves the cell at a rate set by --compressibility / --taup. The default compressibility is
water's (4.57e-5 1/bar), a placeholder 45x larger than silicon's (1/(B x 1e4) = 1.0e-6 with B ~ 98 GPa): with it the
cell responds as if the coupling time were ~11 fs instead of 500 fs. Two runs at 300 K, 0 GPa: the placeholder
(warned) and silicon's own value. Both give a lattice constant above SIESTA's static 5.48024 Ang -- thermal
expansion -- but the placeholder barostat believes silicon is as soft as B ~ 2 GPa and overreacts: the volume
fluctuates ~2-3x more. (8 atoms and 1.5 ps are too few for a converged volume; section 5.7 of the README has the
216-atom runs.)
EOF
for tag in placeholder:4.57e-5 silicon:1.0e-6; do
    name=${tag%%:*}; beta=${tag##*:}
    mlmd case5_npt_$name "$DIR/si8.fdf" $FT --ensemble npt --temperature 300 --compressibility $beta \
        --equilibration-steps 1000 --n-steps 1500 --seed 5
    r="$OUT/case5_npt_$name/stb_mlmd_report.txt"
    echo "--- compressibility $beta 1/bar ($name) ---"
    strip_colors < "$r" | grep -E "^Compressibility|^Lattice \(mean\)|^Barostat|^Volume \(Ang|Stationarity      : Volume" | sed 's/^/  /'
done
echo
echo "[10] of the placeholder run:"
section "$OUT/case5_npt_placeholder/stb_mlmd_report.txt" 10 | sed -n '3,5p'
pause


banner "Case 6: outside the trained range -- silicon squeezed to 5.10 Ang"
cat <<'EOF'
The card says the model was trained on volumes of 18.7-22.4 Ang^3/atom and shortest distances down to 1.85 Ang. A
cell compressed to a = 5.10 Ang (16.6 Ang^3/atom, ~ -20 % in volume, tens of GPa) is outside both: nothing was
measured there. The run is not refused -- the model gives numbers -- but the domain table says OUTSIDE and the
model card's warnings go to [10].
EOF
mkdir -p "$OUT/case6_outside"
sed 's/5\.48023800/5.10000000/g' "$DIR/si8.fdf" > "$OUT/case6_outside/si8_squeezed.fdf"
mlmd case6_outside "$OUT/case6_outside/si8_squeezed.fdf" $FT --ensemble nvt --temperature 300 --equilibration-steps 200 \
    --n-steps 300 --seed 1
r="$OUT/case6_outside/stb_mlmd_report.txt"
section "$r" 8 | sed -n '/Model domain/,$p'
section "$r" 10
pause


banner "Case 7: the same run with the raw MACE-MP-0 foundation model"
cat <<'EOF'
--model medium instead of the custom model: no card, no domain table, an energy on the foundation model's own scale.
The foundation model was never told about this SIESTA setup -- on the card's held-out test set its forces are ~18x
worse (517 vs 29 meV/Ang). At the same cell its pressure differs because its own equilibrium lattice differs.
(First use downloads the model into ~/.cache/mace.)
EOF
mlmd case7_foundation "$DIR/si8.fdf" --model medium --ensemble nvt --temperature 300 --equilibration-steps 200 \
    --n-steps 1500 --seed 1
for r in "$OUT/case1_nvt300/stb_mlmd_report.txt" "$OUT/case7_foundation/stb_mlmd_report.txt"; do
    echo "--- $(basename "$(dirname "$r")") ---"
    strip_colors < "$r" | grep -E "^Model             :|^Model card        :|^SHA-256|^Pressure \(GPa\)|^Temperature \(K\)|^g\(r\)|^Lindemann" | sed 's/^/  /'
done
pause


banner "Case 8: guard rails"
cat <<'EOF'
(a) a custom model without its card; (b) an element the model never saw (Ge, from a Si-only model);
(c) NPT on a slab (25 Ang of cell along z, ~20 Ang of vacuum). Each stops with a clean message, no traceback.
EOF
mkdir -p "$OUT/case8_guards"
cd "$OUT/case8_guards"
echo "--- (a) --custom-model without --custom-model-card ---"
stb-mlmd --file "$DIR/si8.fdf" --custom-model "$MODEL" --n-steps 10 --no-intro < /dev/null 2>&1 | strip_colors | grep -E "ERROR|error" | head -3 || true
echo "--- (b) germanium with a silicon-only model ---"
sed 's/ 1 14 Si/ 1 32 Ge/' "$DIR/si8.fdf" > ge8.fdf
stb-mlmd --file ge8.fdf $FT --n-steps 10 --no-intro < /dev/null 2>&1 | strip_colors | grep -E "ERROR|cannot evaluate" | head -3 || true
echo "--- (c) NPT on a slab ---"
awk '/^ 0.00000000 0.00000000 5.48023800/ {print " 0.00000000 0.00000000 25.00000000"; next} {print}' "$DIR/si8.fdf" | \
    python3 -c "
import sys
lines = sys.stdin.read().splitlines()
out, inside = [], False
for l in lines:
    if 'AtomicCoordinatesAndAtomicSpecies' in l and l.startswith('%block'): inside = True; out.append(l); continue
    if 'AtomicCoordinatesAndAtomicSpecies' in l and l.startswith('%endblock'): inside = False
    if inside:
        x, y, z, s = l.split(); z = float(z) * 5.48023800 / 25.0; l = f' {x} {y} {z:.8f} {s}'
    out.append(l)
print('\n'.join(out))" > slab.fdf
stb-mlmd --file slab.fdf $FT --ensemble npt --compressibility 1e-6 --n-steps 10 --no-intro < /dev/null 2>&1 | strip_colors | grep -E "ERROR|requires a bulk" | head -3 || true
cd "$DIR"
pause


banner "Case 9: the same tool from the interactive menu (stb-suite, code 7.2)"
cat <<'EOF'
The menu asks the same questions and builds the same command line. Answers piped in: structure, the custom model
(Enter accepts the card next to it), a 2x2x2 supercell (64 atoms), NVT, 300 K, 100 + 300 steps of 1 fs, no advanced
settings, a trajectory name, no plots, save the report.
EOF
mkdir -p "$OUT/case9_menu"
cd "$OUT/case9_menu"
printf '7.2\n%s\n%s\n\n2 2 2\n1\n300\n100\n300\n1.0\nn\nmenu_traj.xsf\nn\nn\nn\nn\ny\n\n0\n' "$DIR/si8.fdf" "$MODEL" | \
    stb-suite > menu_log.txt 2>&1 || true
ls -1
strip_colors < stb_mlmd_report.txt | grep -E "^Supercell|^Ensemble|^SHA-256|^Status|^Results" | sed 's/^/  /'
cd "$DIR"
pause


banner "Case 10: the real runs of README section 5 -- 216-512 atoms, 5-10 ps, on a GPU"
cat <<'EOF'
reference_results/ holds the reports and data of the GPU runs (made by make_reference.sh; trajectories left out).
results_figures.py reads them back and prints the tables of README section 5 -- equipartition, timestep, finite
size, thermal expansion vs experiment, the melt -- and redraws the figures in figures/.
EOF
python3 "$DIR/results_figures.py"
pause

banner "Done"
cat <<EOF
Everything this tour wrote is in $OUT/. README.md has the theory (section 1), what each report section means
(section 3) and the worked example on the GPU runs (section 5).
EOF
