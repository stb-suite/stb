#!/bin/bash
# Guided example: Workflow 7.1 -- Custom ML Force Field (all four stages: stb-mlff 7.1.1,
# stb-mlffAnalysis 7.1.2, stb-mlffValidate 7.1.3, stb-mlffActiveLearning 7.1.4).
#
# Not an automated test (see test/7-mlsimulations/1-mlff/{prep,analysis,validate,active-learning}/
# test.sh for that) -- a commented walk-through that runs in about a minute: it runs real
# commands into its own output/ folder and shows the piece of output that proves what just
# happened. Pauses between sections so you can read before moving on. Safe to re-run: it always
# starts by wiping its own output/.
#
# Two kinds of cases:
#   - LIVE: run here, in seconds -- Stage 1 on the graphene 3x3 cell, labelling its folders with
#     the finished SIESTA runs, and one quick turn of the whole loop (train -> validate ->
#     active learning, dry run) on the 2-atom primitive-cell data with a TINY recipe (2 models,
#     5 epochs: its numbers are meaningless, the mechanics are the point).
#   - REFERENCE: the full run of README section 8 (graphene 3x3, 18 atoms: 30 + 6 SIESTA single
#     points, 3-model committees, 100 epochs, ~1 h of CPU training) is NOT repeated -- its real
#     reports are in reference_results/graphene_3x3/ and are shown here.
# Nothing is fabricated: every SIESTA label (SIESTA 5.4.2, PBE, DZP) and every report was produced
# with the commands quoted. This script never runs SIESTA.
#
# Needs the 'ml' extra (pip install stb_suite[ml]) for Case 8 only.

set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

# --show needs a display; Agg keeps the script headless-safe (every figure is still saved).
export MPLBACKEND=Agg

OUT="$DIR/output"
rm -rf "$OUT"
mkdir -p "$OUT"
REF="$DIR/reference_results"
G3="$REF/graphene_3x3/mlff_runs"
P1="$REF/primitive_too_hard"

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

# "run SIESTA": copy the finished outputs of the same folder from reference_results/, after
# checking that the structure written here is the one SIESTA was run on.
label_from_reference() {   # $1 = run folder here, $2 = the same run folder in reference_results
    local n=0 bad=0
    for d in $(cd "$2" && ls -d configs_round_*/mlff_config_* mlff_test/test_* 2>/dev/null); do
        [ -d "$1/$d" ] || continue
        if python3 - "$1/$d/structure.fdf" "$2/$d/structure.fdf" <<'EOF'
import sys, numpy as np
from stb.core import structure_io
a, b = (structure_io.to_pymatgen(structure_io.read_fdf(f)) for f in sys.argv[1:3])
ok = len(a) == len(b) and np.allclose(a.lattice.matrix, b.lattice.matrix, atol=1e-6) and \
     np.allclose(a.cart_coords, b.cart_coords, atol=1e-5)
sys.exit(0 if ok else 1)
EOF
        then
            cp "$2/$d"/siesta.out "$2/$d"/siesta.FA "$2/$d"/siesta.XV "$1/$d/"
            n=$((n + 1))
        else
            echo "  [!] $d differs from the reference run -- not labelled"
            bad=$((bad + 1))
        fi
    done
    echo "  labelled $n folder(s) with the finished SIESTA runs of reference_results/ ($bad differ)"
}


banner "Welcome: a SIESTA-quality force field for graphene, in four stages"
cat <<'EOF'
A machine-learned interatomic potential (MLIP) reproduces the SIESTA potential-energy surface
-- energies, forces, stresses -- far more cheaply, but only where it was trained. This workflow
fine-tunes the MACE-MP foundation model on SIESTA single points, in a loop:

  Stage 1 (stb-mlff, 7.1.1)                 training + held-out test configurations
  SIESTA                                      one single point per folder (you run it)
  Stage 2 (stb-mlffAnalysis, 7.1.2)         fine-tune a committee of models -> round N
  Stage 3 (stb-mlffValidate, 7.1.3)         score every round on the test set; phonons, E(V)
  Stage 4 (stb-mlffActiveLearning, 7.1.4)   let the committee pick what SIESTA labels next

All four talk through one file, mlff_runs/mlff_workflow.json. Cases 1-3 and 8-10 run live (a
minute in total); Cases 4-7 show the real reports of the full run of README section 8.
EOF
pause


banner "Case 1 (live): Stage 1 -- the plan and the configurations"
cat <<'EOF'
structure.fdf is the relaxed graphene 3x3 supercell (stb-supercell -d 3 3 1 of graphene.fdf),
calc.fdf the SIESTA template it is %included in (k 10x10x1 = 30x30 per primitive cell, README
8.1). 16 rattles cycling 4 amplitudes, 5 area strains (graphene is 2D: only the in-plane axes
are scaled), 9 held-out test configs drawn the same way with another seed.
EOF
mkdir -p "$OUT/run"
cp structure.fdf calc.fdf "$OUT/run/"
cd "$OUT/run"
echo "\$ stb-mlff -s structure.fdf -c calc.fdf --pseudo-dir virtual_vault --n-configs 16 \\"
echo "      --stdev 0.015 0.03 0.045 0.06 --volume-strains -3 -1.5 1.5 3 4.5 --n-test 9 --save-report"
stb-mlff -s structure.fdf -c calc.fdf --pseudo-dir virtual_vault --n-configs 16 \
    --stdev 0.015 0.03 0.045 0.06 --volume-strains -3 -1.5 1.5 3 4.5 --n-test 9 \
    --save-report --no-intro > stage1.log 2>&1
section mlff_runs/stb_mlff_report.txt 2 | grep -E "k-grid|Run type|Dispersion"
section mlff_runs/stb_mlff_report.txt 3
cat <<'EOF'

The reference run added --ml-screen (a MACE-MP evaluation of every config before SIESTA; it
needs a minute to load the model). Its columns dE_ML / Fmax_ML turn each amplitude into an
energy and, as 3/2 kT per atom, a temperature (README 8.2):
EOF
section "$G3/stb_mlff_report.txt" 5 | sed -n '4p;6,13p' | awk -F'|' '{printf "%-34s|%-22s|%-7s|%s\n", $2, $5, $12, $13}'
pause


banner "Case 2 (reference): what too large a rattle looks like, before paying SIESTA"
cat <<'EOF'
The first attempt for this example used --stdev up to 0.08 Ang. The foundation-model screen
flagged four configs above 15 eV/Ang -- the workflow's one definition of 'unphysical', shared by
the four stages -- BEFORE a single SIESTA run, and the amplitudes were reduced:
EOF
section "$REF/graphene_3x3/stage1_first_try_stdev008.txt" 9 | head -4
pause


banner "Case 3 (live): SIESTA -- where to run it, and what Stage 2 needs from it"
cat <<'EOF'
Stage 1's last section names the folders: mlff_runs/configs_round_01/ (training, round 1) and
mlff_runs/mlff_test/ (held-out test). Each holds calc.fdf (your template, verbatim, with
'%include config_extra.fdf' on top), structure.fdf (the sampled geometry) and config_extra.fdf
(forces a single point: SIESTA reads the first occurrence of a label). On a cluster: one job per
folder, ~30-60 s each on two cores. Here, the finished runs are copied in from reference_results/
after checking the geometries are identical (same seeds -> same configurations).
EOF
cat mlff_runs/configs_round_01/mlff_config_001/config_extra.fdf
echo
echo "\$ for d in mlff_runs/configs_round_01/mlff_config_* mlff_runs/mlff_test/test_*; do"
echo "      (cd \$d && mpirun -np 2 siesta < calc.fdf > siesta.out); done"
label_from_reference mlff_runs "$G3"
pause


banner "Case 4 (reference): Stage 2, round 1 -- a committee of 3 fine-tuned models"
cat <<'EOF'
$ stb-mlffAnalysis --committee 3 --device cpu --save-data --save-report     (~25 min on 2 cores)
Every folder on one level of theory, the per-element E0 fitted on SIESTA's energy scale, three
members fine-tuned from MACE-MP medium (seeds 123-125, own train/valid split each), every model
evaluated on every config, and a short MD with each model (stability + speed):
EOF
section "$G3/round_01/stb_mlffAnalysis_report.txt" 3 | tail -5
section "$G3/round_01/stb_mlffAnalysis_report.txt" 4 | sed -n '/Out-of-sample pool/,/^ALL/p'
section "$G3/round_01/stb_mlffAnalysis_report.txt" 7 | sed -n '3,6p'
pause


banner "Case 5 (reference): Stage 3, round 1 -- test set, phonons against SIESTA, verdict"
cat <<'EOF'
$ stb-mlffValidate --phonon-ref phonons/phonon_runs --target-force 30 --save-data --save-report
The 9 test configs were never trained on; +/- is a bootstrap over them. --phonon-ref runs the
model on the SAME displaced 5x5 supercells as a finished SIESTA phonon calculation (5.4).
EOF
section "$G3/stb_mlffValidate_report_round1.txt" 3 | sed -n '3,6p'
section "$G3/stb_mlffValidate_report_round1.txt" 5 | grep -E "Averaging|True error|Correlation"
section "$G3/stb_mlffValidate_report_round1.txt" 6 | sed -n '3,6p'
section "$G3/stb_mlffValidate_report_round1.txt" 10 | grep -E "^\[CONTINUE|^\[STOP"
pause


banner "Case 6 (reference): Stage 4 -- molecular dynamics with the model, the committee picks"
cat <<'EOF'
$ stb-mlffActiveLearning --sampling-method md --auto-sampling --top-k 6 --save-data --save-report
NVT trajectories run WITH THE MODEL (minutes, no SIESTA); every frame scored by the committee's
force spread; the unphysical, the ones above the energy cap and the already-known (novelty)
dropped; a diverse batch of 6 written to configs_round_02/. Each energy is also a temperature.
EOF
section "$G3/stb_mlffActiveLearning_report.txt" 1 | grep -E "Calibration|Training energies|Energy cap"
section "$G3/stb_mlffActiveLearning_report.txt" 2 | sed -n '/Auto-sampling scan/,/^Chosen/p'
section "$G3/stb_mlffActiveLearning_report.txt" 4 | grep -E "^Energy above|^Diversity"
pause


banner "Case 7 (reference): round 2 and the MD check -- STOP"
cat <<'EOF'
SIESTA on the 6 picks; stb-mlffAnalysis --init-from previous (round 2); then
$ stb-mlffValidate --phonon-ref phonons/phonon_runs --target-force 30 --compare 1 2 --save-report
and an MD check: 2 ps at 600 K run with round 2's model, 6 frames labelled by SIESTA.
EOF
section "$G3/stb_mlffValidate_report.txt" 3 | sed -n '3,7p'
section "$G3/stb_mlffValidate_report.txt" 6 | sed -n '3,7p'
section "$G3/stb_mlffValidate_report.txt" 8 | grep -E "^Run|^Trend"
section "$G3/stb_mlffValidate_report.txt" 9 | sed -n '/^Overall/,/^Configs/p'
section "$G3/stb_mlffValidate_report.txt" 10 | grep -E "^\[STOP|^\[CONTINUE"
pause


banner "Case 8 (live): the whole loop in a minute, and the trap of sampling too hard"
cat <<'EOF'
reference_results/primitive_too_hard/ is a real run on the 2-atom primitive cell whose two
active-learning rounds chose 0.15 Ang rattles (the old rule: the highest committee spread still
'physical'). Their picks sit 1.2-2.4 eV/atom above equilibrium (~10^4 K), and the test force
error GREW from 21 to 37 meV/Ang. Its validation of the three rounds:
EOF
section "$P1/stb_mlffValidate_report_rounds123.txt" 3 | sed -n '3,7p'
strip_colors < "$P1/stb_mlffValidate_report_rounds123.txt" | grep -E "^\[REGRESSED|^Learning curve|GROWS"
cat <<'EOF'

Now one turn of the loop on its round-1 data (14 training + 8 test configs of 2 atoms), with a
TINY recipe -- 2 models, 5 epochs, no MD stability test, no foundation comparison -- so it
runs in seconds. Its errors are far above a real recipe's; watch the mechanics, and above all
the energy cap of Stage 4 (--dry-run: nothing written).
EOF
mkdir -p "$OUT/primitive/mlff_runs"
cp "$P1/structure.fdf" "$P1/calc.fdf" "$OUT/primitive/"
cp -r "$P1/mlff_runs/configs_round_01" "$P1/mlff_runs/mlff_test" "$OUT/primitive/mlff_runs/"
cp "$P1/mlff_workflow.after_stage1.json" "$OUT/primitive/mlff_runs/mlff_workflow.json"
cd "$OUT/primitive"
if python3 -c "import mace" 2>/dev/null; then
    echo "\$ stb-mlffAnalysis --committee 2 --epochs 5 --device cpu --no-stability --skip-foundation-comparison"
    stb-mlffAnalysis --committee 2 --epochs 5 --device cpu --no-stability --skip-foundation-comparison \
        --no-intro > stage2.log 2>&1
    section stage2.log 3 | grep -E "^Seed|^12[34] "
    echo "\$ stb-mlffValidate --skip-foundation-comparison"
    stb-mlffValidate --skip-foundation-comparison --no-intro > stage3.log 2>&1 || true
    section stage3.log 3 | sed -n '3,5p'
    echo "\$ stb-mlffActiveLearning --auto-sampling --top-k 5 --dry-run      (rattle, the default)"
    stb-mlffActiveLearning --auto-sampling --top-k 5 --dry-run --no-intro > stage4_dry.log 2>&1 || true
    section stage4_dry.log 1 | grep -E "Training energies|Energy cap"
    section stage4_dry.log 2 | sed -n '/Auto-sampling scan/,/^Chosen/p'
    section stage4_dry.log 4 | grep -E "^Energy above"
    cat <<'EOF'

The old rule took the most uncertain setting still under 15 eV/Ang (0.15 Ang, ~8000 K here);
the cap -- the training set's 90th-percentile energy x --frontier 1.5 -- keeps the next round one
step beyond the labelled range instead (README section 9).
EOF
else
    echo "  [skipped: needs the 'ml' extra -- pip install stb_suite[ml]]"
fi
cd "$OUT/run"
pause


banner "Case 9 (live): the phonon reference itself -- a k-grid lesson (README section 8.1)"
cat <<'EOF'
A model can only be as good as its labels. The Gamma E2g phonon of graphene by frozen phonon in
the primitive cell (atoms displaced by +u and -u), SIESTA at four k-grids:
EOF
python3 "$DIR/e2g_frozen_phonon.py" "$REF/kgrid_e2g"
cat <<'EOF'

A grid divisible by 3 contains the Dirac point K exactly; a COARSE one (15x15) then exaggerates
graphene's Kohn anomaly and makes the frequency depend on the amplitude. The primitive-cell run
of Case 8 used 15x15; this example's labels use 30x30 (within ~0.5 THz of 45x45).
EOF
pause


banner "Case 10 (live): the interactive menu path (stb-suite, shortcut 7.1.1) gives the same folders"
cat <<'EOF'
The same Stage 1 as Case 1, answered through the menu: structure and calc files (Enter =
defaults), pseudopotentials 2 = virtual_vault, 16 configs, the 4 amplitudes, the 5 strains, 9 test
configs, no D3, no AIMD, run folder menu_runs, no report, no ML screen, no export, no viewer, no
figure.
EOF
mkdir -p "$OUT/menu"
cp structure.fdf calc.fdf "$OUT/menu/"
(cd "$OUT/menu" && printf '7.1.1\n\n\n2\n16\n0.015 0.03 0.045 0.06\n-3 -1.5 1.5 3 4.5\n9\nn\nn\nmenu_runs\nn\nn\n\nn\nn\n\n0\n' \
    | stb-suite > menu.log 2>&1) || true
n_same=0; n_diff=0
for d in mlff_runs/configs_round_01/mlff_config_* mlff_runs/mlff_test/test_*; do
    m="$OUT/menu/menu_runs/${d#mlff_runs/}"
    if cmp -s "$d/structure.fdf" "$m/structure.fdf"; then n_same=$((n_same + 1)); else n_diff=$((n_diff + 1)); fi
done
echo "  structure.fdf identical in $n_same folder(s), different in $n_diff"
pause


banner "Done"
cat <<EOF
Everything is in $OUT:
  run/mlff_runs/         Stage 1 of the graphene 3x3 run, labelled with the reference SIESTA runs --
                         ready for 'stb-mlffAnalysis' if you want to train the real committee
                         yourself (~25 min on 2 cores; README section 8.4 has the reference numbers)
  primitive/mlff_runs/   the one-minute loop of Case 8 (round_01/, reports, figures)
  menu/                  Case 10's menu run

The full reports of README section 8 are in reference_results/graphene_3x3/mlff_runs/.
EOF
