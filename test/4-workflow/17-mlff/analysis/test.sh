#!/bin/bash

# --- Setup ---
# Smoke test for stb-mlffAnalysis (Custom ML Force Field Train, item 4.17.2)
#
# Fixtures (mlff_config_001..008) are 8 REAL finished SIESTA calculations:
# an O2-in-vacuum reference (SZ basis) rattled 8 ways by a real stb-mlff run,
# then actually run through the system's installed SIESTA (v5.4.2) to
# produce genuine calc.fdf/.out/.FA/.XV in each folder -- this whole
# pipeline (aggregation -> training_set.xyz -> real mace_run_train
# fine-tuning -> loadable .model) was verified live end-to-end during
# development. Needs the optional 'ml' extra (mace-torch) -- skips
# gracefully (not a failure) if mace_run_train isn't on PATH.
#
# The fixtures sit in the older flat layout (mlff_config_* in the run folder);
# section 6 copies them into the per-round layout (configs_round_01/).
export MPLBACKEND=Agg
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$FIXTURE_DIR/test_files"
FAST="--epochs 2 --batch-size 4 --device cpu --foundation-model small --no-stability"

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
NC='\033[0m'

PASS=0
FAIL=0

check_success() {
    if [ -s "$1" ]; then
        echo -e " ... ${GREEN}OK${NC} (file '$1' created)"
        PASS=$((PASS+1))
    else
        echo -e " ... ${RED}FAIL${NC} (file '$1' was not created)"
        FAIL=$((FAIL+1))
    fi
}

check_contains() {
    if grep -q -- "$1" "$2" 2>/dev/null; then
        echo -e "   -> ${GREEN}Verified:${NC} '$1' found in '$2'"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} '$1' NOT found in '$2'"
        FAIL=$((FAIL+1))
    fi
}

check_not_contains() {
    if grep -q -- "$1" "$2" 2>/dev/null; then
        echo -e "   -> ${RED}Failed:${NC} '$1' unexpectedly found in '$2'"
        FAIL=$((FAIL+1))
    else
        echo -e "   -> ${GREEN}Verified:${NC} '$1' absent from '$2'"
        PASS=$((PASS+1))
    fi
}

check_exit_code() {
    if [ "$1" -eq "$2" ]; then
        echo -e "   -> ${GREEN}Verified:${NC} exit code $1 (expected $2)"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} exit code $1 (expected $2)"
        FAIL=$((FAIL+1))
    fi
}

check_python() {   # $1 = description, $2 = python code exiting 0 on success
    python3 -c "$2"
    if [ $? -eq 0 ]; then
        echo -e "   -> ${GREEN}Verified:${NC} $1"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} $1"
        FAIL=$((FAIL+1))
    fi
}


echo "--- Starting tester for stb-mlffAnalysis (item 4.17.2) ---"

if ! command -v mace_run_train &> /dev/null; then
    echo -e "${YELLOW}SKIPPED${NC}: 'mace_run_train' not on PATH (needs 'pip install stb_suite[ml]')."
    exit 0
fi

rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
for d in "$FIXTURE_DIR"/mlff_config_*; do
    cp -r "$d" "$TEST_DIR/"
done
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 1. Basic fine-tuning run (fast: few epochs, small model) ---
echo -e "\n--- Testing a basic fine-tuning run (--epochs 5, MD stability check on) ---"
stb-mlffAnalysis --path "mlff_config_*" --committee 1 -o . --epochs 5 --batch-size 4 --device cpu \
    --foundation-model small --stability-steps 20 --name test_model --save-data --show \
    --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "Folders           : 8 listed, 8 usable, 0 skipped" log_basic.txt
check_contains "Member 1/1 -- seed 123" log_basic.txt
check_contains "-> ./test_model_compiled.model" log_basic.txt
check_success test_model_compiled.model
check_success train_set.xyz
check_success valid_set.xyz
check_success mace_run_train.log
check_success mace_run_train.cmd
check_python "train_set.xyz + valid_set.xyz together have 8 frames with REF_energy/REF_forces" "
import sys
from ase.io import read
frames = read('train_set.xyz', index=':') + read('valid_set.xyz', index=':')
sys.exit(0 if len(frames) == 8 and all('REF_energy' in a.info and 'REF_forces' in a.arrays for a in frames) else 1)
"

echo "Testing: the report sections, from the data to the warnings"
for sec in "\[0\] RUN METADATA" "\[1\] TRAINING DATA" "\[2\] TRAINING SETUP" "\[3\] TRAINING" \
           "\[4\] ERRORS ON THE DATA" "\[6\] COMPARISON WITH EARLIER ROUNDS" "\[7\] MODEL CHECKS" \
           "\[8\] SUMMARY & FILES" "\[9\] NEXT STEPS" "\[10\] WARNINGS"; do
    check_contains "$sec" log_basic.txt
done
check_contains "\[OK\] Same level of theory in all 8 folders" log_basic.txt
check_contains "\[OK\] SCF converged in every config" log_basic.txt
check_contains "Library messages  :" log_basic.txt
check_not_contains "WARNING:root" log_basic.txt
check_not_contains "cuequivariance" log_basic.txt
check_python "every warning is listed after [10] WARNINGS (none in the middle of the report)" "
import sys, re
t = open('log_basic.txt').read()
head, tail = t.split('[10] WARNINGS', 1)
sys.exit(0 if '[WARNING]' not in head and re.search(r' 1\. \w', tail) else 1)
"

echo "Testing: errors vs the raw foundation model (its absolute E is a scale offset, not an error)"
check_contains "foundation (small) | all" log_basic.txt
check_contains "(scale offset)" log_basic.txt
check_contains "Fine-tuning gain  :" log_basic.txt
check_contains "Worst-predicted configs" log_basic.txt
check_success test_model_parity.png
check_success mlff_training.png
check_success mlff_training_errors.dat
check_contains "--show needs a display" log_basic.txt

echo "Testing: MD stability and speed of the model"
check_contains "MD cell           :" log_basic.txt
check_contains "Speed             :" log_basic.txt
check_python "stability and speed stored in the round" "
import json, sys
r = json.load(open('mlff_workflow.json'))['rounds'][0]
st = r['model_checks']['stability'][0]
sys.exit(0 if st['status'] in ('OK', 'UNSTABLE') and r['model_checks']['speed'] and r['n_usable'] == 8 else 1)
"

echo "Testing: stress is actually in the loss (--loss stress), lr 0.001, E0s from the foundation residual"
check_contains "stress 10000 (--loss stress)" log_basic.txt
check_contains "Learning rate     : 0.001" log_basic.txt
check_contains "\-\-loss stress" mace_run_train.cmd
check_contains "\-\-stress_weight 10000.0" mace_run_train.cmd
check_contains "\-\-lr 0.001" mace_run_train.cmd
check_contains "\-\-E0s {8: " mace_run_train.cmd
check_contains "WeightedEnergyForcesStressLoss" mace_run_train.log
check_contains "E0s             : O " log_basic.txt
check_contains "('average'" log_basic.txt
check_contains "start E RMSE" log_basic.txt

echo "Testing: --e0s average still passes mace_run_train's own fit"
stb-mlffAnalysis --path "mlff_config_*" --committee 1 -o e0avg $FAST --name test_model_avg --e0s average \
    --skip-foundation-comparison --workflow wf_e0avg.json --no-intro > log_e0avg.txt 2>&1
check_exit_code $? 0
check_contains "\-\-E0s average" e0avg/mace_run_train.cmd
rm -rf e0avg wf_e0avg.json


# --- 1b. --export-lammps ---
echo -e "\n--- Testing --export-lammps ---"
if command -v mace_create_lammps_model &> /dev/null; then
    stb-mlffAnalysis --path "mlff_config_*" --committee 1 -o lammps $FAST --name test_model_lammps \
        --export-lammps --skip-foundation-comparison --workflow wf_lammps.json --no-intro > log_lammps.txt 2>&1
    check_exit_code $? 0
    check_success lammps/test_model_lammps_compiled.model-lammps.pt
    check_contains "test_model_lammps_compiled.model-lammps.pt" log_lammps.txt
    rm -rf lammps wf_lammps.json
else
    echo -e "   -> ${YELLOW}SKIPPED${NC}: 'mace_create_lammps_model' not on PATH."
fi


# --- 2. Too few usable configs (expects failure) ---
echo -e "\n--- Testing a folder with too few usable configs (expects failure) ---"
mkdir -p too_few/mlff_config_001
cp mlff_config_001/calc.fdf mlff_config_001/calc.out too_few/mlff_config_001/ 2>/dev/null
(cd too_few && stb-mlffAnalysis --path "mlff_config_*" --committee 1 -o . --no-intro > ../log_toofew.txt 2>&1)
check_exit_code $? 1
check_contains "Only" log_toofew.txt
rm -rf too_few


# --- 3. No matching folders ---
echo -e "\n--- Testing --path matching nothing (expects failure) ---"
stb-mlffAnalysis --path "no_such_folder_*" --no-intro > log_nomatch.txt 2>&1
check_exit_code $? 1
check_contains "No config folders found (--path 'no_such_folder_\*')" log_nomatch.txt


# --- 3b. Committee + rounds through the workflow file ---
echo -e "\n--- Testing --committee 2 and rounds recorded in a workflow file ---"
rm -f wf_test.json
stb-mlffAnalysis --path "mlff_config_*" --committee 2 $FAST --workflow wf_test.json -o c_round1 \
    --skip-foundation-comparison --no-intro > log_committee1.txt 2>&1
check_exit_code $? 0
check_contains "Committee         : 2 model(s), seed(s) 123, 124" log_committee1.txt
check_contains "Member 2/2 -- seed 124" log_committee1.txt
check_contains "\[5\] COMMITTEE" log_committee1.txt
check_contains "Force spread      :" log_committee1.txt
check_success c_round1/m_s123/train_set.xyz
check_success c_round1/m_s124/valid_set.xyz
check_success wf_test.json
echo "Testing: the second run takes the configs from the workflow file, goes to round_02/, continues from round 1"
stb-mlffAnalysis --committee 2 $FAST --workflow wf_test.json --init-from previous \
    --skip-foundation-comparison --no-intro > log_committee2.txt 2>&1
check_exit_code $? 0
check_contains "Configs from      : wf_test.json (training configs)" log_committee2.txt
check_contains "Round directory   : round_02" log_committee2.txt
check_contains "Starting weights  : c_round1/m_s123/mlff_model.model, c_round1/m_s124/mlff_model.model" log_committee2.txt
check_contains "2 (this)" log_committee2.txt
check_python "2 rounds, 2 models each, paths/metrics/history/init_from recorded" "
import json, os, sys
d = json.load(open('wf_test.json'))
r1, r2 = d['rounds']
ok = (d['format'] == 'stb-mlff-workflow' and [r['round'] for r in d['rounds']] == [1, 2]
      and r1['dir'] == 'c_round1' and r2['dir'] == 'round_02'
      and len(r2['committee']) == 2 and all(os.path.isfile(m['model']) for m in r1['committee'] + r2['committee'])
      and [m['seed'] for m in r2['committee']] == [123, 124]
      and all('F_rmse_meV' in m['valid'] for m in r2['committee'])
      and r1['init_from'] == 'foundation' and r2['init_from'] == 'previous'
      and {'out_of_sample', 'committee'} <= set(r2['training_errors'])
      and len(r2['configs']) == 8 and [h['stage'] for h in d['history']] == ['train', 'train'])
sys.exit(0 if ok else 1)
"
echo "Testing: --committee 0 is rejected"
stb-mlffAnalysis --committee 0 --workflow wf_test.json --no-intro > log_committee0.txt 2>&1
check_exit_code $? 1
check_contains "--committee must be >= 1" log_committee0.txt
rm -rf c_round1 round_02 wf_test.json


# --- 3c. --cv K: every config gets an out-of-sample prediction ---
echo -e "\n--- Testing --cv 4 ---"
stb-mlffAnalysis --path "mlff_config_*" --cv 4 $FAST -o cv_round --workflow wf_cv.json --no-intro > log_cv.txt 2>&1
check_exit_code $? 0
check_contains "4-fold cross-validation" log_cv.txt
check_contains "Committee         : 4 model(s)" log_cv.txt
check_contains "Out-of-sample pool: 8/8 configs" log_cv.txt
check_python "the round records the 4 folds" "
import json, sys
r = json.load(open('wf_cv.json'))['rounds'][0]
sys.exit(0 if r['training_errors']['cv_folds'] == 4 and len(r['committee']) == 4 else 1)
"
rm -rf cv_round wf_cv.json


# --- 4. --save-report ---
echo -e "\n--- Testing --save-report (written in the round directory) ---"
stb-mlffAnalysis --path "mlff_config_*" --committee 1 -o rep $FAST --name test_model2 --save-report \
    --skip-foundation-comparison --workflow wf_rep.json --no-intro > log_savereport.txt 2>&1
check_success rep/stb_mlffAnalysis_report.txt
check_contains "STB-MLFFANALYSIS REPORT" rep/stb_mlffAnalysis_report.txt
check_contains "\[10\] WARNINGS" rep/stb_mlffAnalysis_report.txt
rm -rf rep wf_rep.json


# --- 5. Robustness / errors ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: --version"
stb-mlffAnalysis --version > log_version.txt 2>&1
check_contains "stb-mlffAnalysis" log_version.txt

echo "Testing: --help documents the options"
stb-mlffAnalysis --help > log_help.txt 2>&1
for opt in foundation-model epochs valid-fraction export-lammps skip-foundation-comparison stress-weight \
           no-stress-training foundation-residual committee workflow "\-\-cv" init-from no-stability \
           view-worst save-data skip-unlabelled; do
    check_contains "$opt" log_help.txt
done

echo "Testing: --skip-foundation-comparison skips the foundation-model evaluation"
stb-mlffAnalysis --path "mlff_config_*" --committee 1 -o nofound $FAST --name test_model_nofound \
    --skip-foundation-comparison --workflow wf_nofound.json --no-intro > log_nofound.txt 2>&1
check_exit_code $? 0
check_contains "seed 123 *| train" log_nofound.txt
check_not_contains "foundation (small) | all" log_nofound.txt
rm -rf nofound wf_nofound.json

echo "Testing: --no-stress-training keeps the energy+forces loss"
stb-mlffAnalysis --path "mlff_config_*" --committee 1 -o nostress $FAST --name test_model_nostress \
    --no-stress-training --skip-foundation-comparison --workflow wf_nostress.json --no-intro > log_nostress.txt 2>&1
check_exit_code $? 0
check_contains "no-stress-training although every config has a stress" log_nostress.txt
check_not_contains "\-\-loss stress" nostress/mace_run_train.cmd
rm -rf nostress wf_nostress.json

echo "Testing: a level-of-theory mismatch between folders is reported"
rm -rf theory && mkdir theory && cp -r mlff_config_00[1-5] theory/
sed -i 's/^MeshCutoff .*/MeshCutoff          200 Ry/' theory/mlff_config_005/calc.fdf
(cd theory && stb-mlffAnalysis --path "mlff_config_*" --committee 1 -o . $FAST --skip-foundation-comparison \
    --no-intro > ../log_theory.txt 2>&1)
check_contains "level of theory differs in 1 folder(s) for Mesh.Cutoff" log_theory.txt
rm -rf theory


# --- 6. Per-round layout: configs_round_NN/, unlabelled folders, the Stage-4 guard ---
echo -e "\n--- Testing the per-round layout (configs_round_01/, configs_round_02/) ---"
rm -rf layout && mkdir -p layout/configs_round_01 && cp -r mlff_config_* layout/configs_round_01/
mkdir -p layout/configs_round_01/mlff_config_009 && cp mlff_config_001/calc.fdf layout/configs_round_01/mlff_config_009/
pushd layout > /dev/null
stb-mlffAnalysis --committee 1 $FAST --skip-foundation-comparison --no-intro > log_layout1.txt 2>&1
check_exit_code $? 0
check_contains "Configs from      : '\[configs_round_\*/\]mlff_config_\*'" log_layout1.txt
check_contains "Folders           : 9 listed, 8 usable, 1 skipped" log_layout1.txt
check_contains "By folder         : configs_round_01 8/9 labelled" log_layout1.txt
check_contains "SIESTA: 1 of 9 folder(s) in configs_round_01/ have no finished run" log_layout1.txt
check_success mlff_workflow.json
echo "Testing: a batch meant for round 2 with no finished SIESTA run stops Stage 2 (and says where to run it)"
mkdir -p configs_round_02/mlff_config_010 && cp ../mlff_config_002/calc.fdf configs_round_02/mlff_config_010/
python3 - <<'PY'
import json
d = json.load(open('mlff_workflow.json'))
d['configs'].append({'folder': 'configs_round_02/mlff_config_010', 'set': 'train', 'kind': 'al-rattle',
                     'origin': 'test', 'stage': 'active-learning', 'round': 1, 'n_atoms': 2})
json.dump(d, open('mlff_workflow.json', 'w'), indent=1)
PY
stb-mlffAnalysis --committee 1 $FAST --skip-foundation-comparison --no-intro > log_layout2.txt 2>&1
check_exit_code $? 1
check_contains "None of the 1 configs meant for round 2 (configs_round_02/)" log_layout2.txt
check_contains "for d in configs_round_02/mlff_config_\*" log_layout2.txt
stb-mlffAnalysis --committee 1 $FAST --skip-foundation-comparison --skip-unlabelled --no-intro > log_layout3.txt 2>&1
check_exit_code $? 0
check_contains "SIESTA: 1 of 1 folder(s) in configs_round_02/ have no finished run" log_layout3.txt
check_contains "Round directory   : round_02" log_layout3.txt
popd > /dev/null
rm -rf layout


# --- 7. Interactive path (stb-suite, shortcut 4.17.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.17.2) ---"
# mlff_workflow.json (from section 1) is here: no path / work-dir prompts. Answers:
# small model, no CV, committee 1, not from the previous round, 3 epochs, batch 4,
# default stress weight, cpu, name, no stability test, no worst configs, no LAMMPS,
# no --show, no data file, no report; Enter; 0 = quit.
printf '4.17.2\n1\n\n1\nn\n3\n4\n\ncpu\ntest_model3\nn\n\nn\nn\nn\nn\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_contains "mlff_workflow.json: 8 training" log_menu.txt
check_contains "Configs from      : mlff_workflow.json" log_menu.txt
check_contains "Round directory   : round_02" log_menu.txt
check_success round_02/test_model3_compiled.model


popd > /dev/null

# --- 7b. Level of theory: k x a measured in reference cells ---
echo "Testing: strained cells and divided supercell grids with the reference's k-points do not warn"
check_python "same grid on a +-2 % cell and a 3x supercell with a third of the grid match; other grids do not" "
import sys
from stb.mlff_analysis import _kxa_in_reference_cell
L = 2.483
def sig(k, l):
    return {'k x a (Ang)': 'raw', '_kgrid': (k, l)}
s = [sig((15, 15, 1), (L, L, 20.0)), sig((15, 15, 1), (L * 0.98, L * 0.98, 20.0)),
     sig((15, 15, 1), (L * 1.02, L * 1.02, 20.0)), sig((8, 8, 1), (2 * L, 2 * L, 20.0)),
     sig((18, 18, 1), (L, L, 20.0)), sig((5, 5, 1), (3 * L, 3 * L, 20.0))]
_kxa_in_reference_cell(s)
v = [x['k x a (Ang)'] for x in s]
sys.exit(0 if v[0] == v[1] == v[2] == v[5] and v[3] != v[4] and v[4] != v[0] and all('_kgrid' not in x for x in s)
         and v[0].endswith('G') else 1)
"

# --- 8. Summary ---
echo -e "\n--- Tests Complete ---"
echo -e "${GREEN}Passed: $PASS${NC}   ${RED}Failed: $FAIL${NC}"

read -p "Remove the '$TEST_DIR' directory and all test files? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "Cleaning up test files..."
    rm -rf "$TEST_DIR"
    echo -e "${GREEN}Cleanup complete.${NC}"
else
    echo "Test files were kept in '$TEST_DIR/' for inspection."
fi

[ "$FAIL" -eq 0 ]
