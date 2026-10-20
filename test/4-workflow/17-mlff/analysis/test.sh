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
export MPLBACKEND=Agg
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$FIXTURE_DIR/test_files"

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

check_exit_code() {
    if [ "$1" -eq "$2" ]; then
        echo -e "   -> ${GREEN}Verified:${NC} exit code $1 (expected $2)"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} exit code $1 (expected $2)"
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
echo -e "\n--- Testing a basic fine-tuning run (--epochs 5 --device cpu) ---"
stb-mlffAnalysis --path "mlff_config_*" --committee 1 -o . --epochs 5 --batch-size 4 \
    --device cpu --name test_model --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "8 usable configuration(s) (0 skipped)" log_basic.txt
check_contains "mace_run_train finished" log_basic.txt
check_contains "Fine-tuned model" log_basic.txt
check_contains "Model loads and runs" log_basic.txt
check_success train_set.xyz
check_success valid_set.xyz
check_success mace_run_train.log

echo "Testing: train_set.xyz + valid_set.xyz together have 8 frames with REF_energy/REF_forces"
python3 -c "
import sys
from ase.io import read
frames = read('train_set.xyz', index=':') + read('valid_set.xyz', index=':')
sys.exit(0 if len(frames) == 8 and all(
    'REF_energy' in a.info and 'REF_forces' in a.arrays for a in frames
) else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} train_set.xyz + valid_set.xyz have 8 frames total with REF_energy/REF_forces"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} train/valid sets missing expected frames/fields"
    FAIL=$((FAIL+1))
fi

echo "Testing: parity plot + comparison against the raw foundation model"
check_success test_model_parity.png
check_contains "Fine-tuned:" log_basic.txt
check_contains "Foundation (small, no fine-tuning):" log_basic.txt
check_contains "E RMSE (offset removed)" log_basic.txt
check_contains "Stress RMSE" log_basic.txt

echo "Testing: the foundation block has no absolute E RMSE (it is the SIESTA/MACE scale offset)"
python3 -c "
import sys
text = open('log_basic.txt').read()
block = text.split('Foundation (small, no fine-tuning):', 1)[1].split('validation config', 1)[0]
sys.exit(0 if 'absolute' not in block else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} no absolute E RMSE reported for the foundation model"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} absolute E RMSE reported for the foundation model"
    FAIL=$((FAIL+1))
fi

echo "Testing: stress is actually in the loss (--loss stress), not only computed"
check_contains "trained on, with weight 10000" log_basic.txt
check_contains "--loss stress --stress_weight 10000.0" log_basic.txt
check_contains "--lr 0.001" log_basic.txt
check_contains "WeightedEnergyForcesStressLoss" mace_run_train.log

echo "Testing: E0s default to the foundation residual (explicit {Z: E0} on the SIESTA scale)"
check_contains "--E0s {8: " log_basic.txt
check_contains "Starting E RMSE" log_basic.txt
check_contains "\"average\" would give" log_basic.txt

echo "Testing: --e0s average still passes mace_run_train's own fit"
stb-mlffAnalysis --path "mlff_config_*" --committee 1 -o . --epochs 2 --batch-size 4 --device cpu \
    --name test_model_avg --e0s average --skip-foundation-comparison --no-intro > log_e0avg.txt 2>&1
check_exit_code $? 0
check_contains "--E0s average" log_e0avg.txt


# --- 1b. --export-lammps ---
echo -e "\n--- Testing --export-lammps ---"
if command -v mace_create_lammps_model &> /dev/null; then
    stb-mlffAnalysis --path "mlff_config_*" --committee 1 -o . --epochs 3 --batch-size 4 --device cpu \
        --name test_model_lammps --export-lammps --no-intro > log_lammps.txt 2>&1
    check_exit_code $? 0
    check_contains "Exported LAMMPS model" log_lammps.txt
    check_success test_model_lammps_compiled.model-lammps.pt
else
    echo -e "   -> ${YELLOW}SKIPPED${NC}: 'mace_create_lammps_model' not on PATH."
fi


# --- 2. Missing mace_run_train dependency check (simulated via empty PATH) ---
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
stb-mlffAnalysis --path "mlff_config_*" --committee 2 --epochs 2 --batch-size 4 --device cpu \
    --workflow wf_test.json -o c_round1 --skip-foundation-comparison --no-intro > log_committee1.txt 2>&1
check_exit_code $? 0
check_contains "Committee         : 2 model(s), seed(s) 123, 124" log_committee1.txt
check_contains "MEMBER 2/2 (seed 124)" log_committee1.txt
check_contains "COMMITTEE SUMMARY" log_committee1.txt
check_success c_round1/m_s123/train_set.xyz
check_success c_round1/m_s124/valid_set.xyz
check_success wf_test.json
echo "Testing: the second run takes the configs from the workflow file and goes to round_02/"
stb-mlffAnalysis --committee 2 --epochs 2 --batch-size 4 --device cpu \
    --workflow wf_test.json --skip-foundation-comparison --no-intro > log_committee2.txt 2>&1
check_exit_code $? 0
check_contains "Configs from      : wf_test.json (training configs)" log_committee2.txt
check_contains "Round directory   : round_02" log_committee2.txt
python3 -c "
import json, os, sys
d = json.load(open('wf_test.json'))
r1, r2 = d['rounds']
ok = (d['format'] == 'stb-mlff-workflow' and [r['round'] for r in d['rounds']] == [1, 2]
      and r1['dir'] == 'c_round1' and r2['dir'] == 'round_02'
      and len(r2['committee']) == 2 and all(os.path.isfile(m['model']) for m in r1['committee'] + r2['committee'])
      and [m['seed'] for m in r2['committee']] == [123, 124]
      and all('F_rmse_meV' in m['valid'] for m in r2['committee'])
      and len(r2['configs']) == 8 and [h['stage'] for h in d['history']] == ['train', 'train'])
sys.exit(0 if ok else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} 2 rounds, 2 models each, paths/metrics/history recorded"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} workflow file does not hold the expected rounds"
    FAIL=$((FAIL+1))
fi
echo "Testing: --committee 0 is rejected"
stb-mlffAnalysis --committee 0 --workflow wf_test.json --no-intro > log_committee0.txt 2>&1
check_exit_code $? 1
check_contains "--committee must be >= 1" log_committee0.txt
rm -rf c_round1 round_02 wf_test.json


# --- 4. --save-report ---
echo -e "\n--- Testing --save-report ---"
rm -f stb_mlffAnalysis_report.txt
stb-mlffAnalysis --path "mlff_config_*" --committee 1 -o . --epochs 2 --batch-size 4 --device cpu \
    --name test_model2 --save-report --no-intro > log_savereport.txt 2>&1
check_success stb_mlffAnalysis_report.txt
check_contains "STB-MLFFANALYSIS REPORT" stb_mlffAnalysis_report.txt
rm -f stb_mlffAnalysis_report.txt


# --- 5. Robustness / errors ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: --version"
stb-mlffAnalysis --version > log_version.txt 2>&1
check_contains "stb-mlffAnalysis" log_version.txt

echo "Testing: --help documents --foundation-model, --epochs, --valid-fraction, --export-lammps, --skip-foundation-comparison"
stb-mlffAnalysis --help > log_help.txt 2>&1
check_contains "foundation-model" log_help.txt
check_contains "epochs" log_help.txt
check_contains "valid-fraction" log_help.txt
check_contains "export-lammps" log_help.txt
check_contains "skip-foundation-comparison" log_help.txt
check_contains "stress-weight" log_help.txt
check_contains "no-stress-training" log_help.txt
check_contains "foundation-residual" log_help.txt
check_contains "committee" log_help.txt
check_contains "workflow" log_help.txt

echo "Testing: --skip-foundation-comparison skips the foundation-model evaluation"
stb-mlffAnalysis --path "mlff_config_*" --committee 1 -o . --epochs 2 --batch-size 4 --device cpu \
    --name test_model_nofound --skip-foundation-comparison --no-intro > log_nofound.txt 2>&1
check_exit_code $? 0
check_contains "Fine-tuned:" log_nofound.txt
grep -q "Foundation (" log_nofound.txt
if [ $? -ne 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} foundation metrics correctly absent with --skip-foundation-comparison"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} foundation metrics unexpectedly present"
    FAIL=$((FAIL+1))
fi

echo "Testing: --no-stress-training keeps the energy+forces loss"
stb-mlffAnalysis --path "mlff_config_*" --committee 1 -o . --epochs 2 --batch-size 4 --device cpu \
    --name test_model_nostress --no-stress-training --skip-foundation-comparison \
    --no-intro > log_nostress.txt 2>&1
check_exit_code $? 0
check_contains "no-stress-training: training on energy+forces only" log_nostress.txt
grep -q -- "--loss stress" log_nostress.txt
if [ $? -ne 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} no '--loss stress' with --no-stress-training"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} '--loss stress' passed despite --no-stress-training"
    FAIL=$((FAIL+1))
fi


# --- 6. Interactive path (stb-suite, shortcut 4.17.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.17.2) ---"
# mlff_workflow.json exists here (written by the runs above): no path / work-dir prompts.
printf '4.17.2\n1\n1\n3\n4\n\ncpu\ntest_model3\nn\nn\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_contains "mlff_workflow.json: 8 training" log_menu.txt
check_contains "Fine-tuned model" log_menu.txt
check_contains "Configs from      : mlff_workflow.json" log_menu.txt


popd > /dev/null

# --- 7. Summary ---
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
