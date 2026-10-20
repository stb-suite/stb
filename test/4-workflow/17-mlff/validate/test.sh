#!/bin/bash

# --- Setup ---
# Smoke test for stb-mlffValidate (Custom ML Force Field Validate, item 4.17.3)
#
# Reuses the 8 real finished-SIESTA fixtures from ../analysis/: 6 become the
# training set, 2 the held-out test set (mlff_test/test_00N, listed in
# mlff_workflow.json as Stage 1 --n-test would). Two quick training rounds
# (2-model committees, --epochs 2) are scored round by round. The phonon
# comparison (--phonon-ref) needs a finished stb-phononsCreate run of the same
# system and is not exercised here. Needs the optional 'ml' extra -- skips
# gracefully (not a failure) if mace_run_train isn't on PATH.
export MPLBACKEND=Agg
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MLFF_FIXTURES="$FIXTURE_DIR/../analysis"
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


echo "--- Starting tester for stb-mlffValidate (item 4.17.3) ---"

if ! command -v mace_run_train &> /dev/null; then
    echo -e "${YELLOW}SKIPPED${NC}: 'mace_run_train' not on PATH (needs 'pip install stb_suite[ml]')."
    exit 0
fi
if [ ! -d "$MLFF_FIXTURES/mlff_config_001" ]; then
    echo -e "${YELLOW}SKIPPED${NC}: ../analysis fixtures not found."
    exit 0
fi

rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR/mlff_test"
for i in 1 2 3 4 5 6; do cp -r "$MLFF_FIXTURES/mlff_config_00$i" "$TEST_DIR/"; done
cp -r "$MLFF_FIXTURES/mlff_config_007" "$TEST_DIR/mlff_test/test_001"
cp -r "$MLFF_FIXTURES/mlff_config_008" "$TEST_DIR/mlff_test/test_002"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null

# A workflow file as Stage 1 would write it: 6 training + 2 test configs.
python3 -c "
from stb.core import mlff_workflow as wf
s = wf.new_state('mlff_config_001/calc.fdf', 2, ['O'], '.')
wf.add_configs(s, [{'folder': f'mlff_config_00{i}', 'set': 'train', 'kind': 'rattle', 'origin': 'fixture',
                    'stage': 'prep', 'round': 0} for i in range(1, 7)], '.')
wf.add_configs(s, [{'folder': 'mlff_test/test_001', 'set': 'test', 'kind': 'rattle', 'origin': 'fixture', 'stage': 'prep', 'round': 0},
                   {'folder': 'mlff_test/test_002', 'set': 'test', 'kind': 'strain', 'origin': 'fixture', 'stage': 'prep', 'round': 0}], '.')
wf.save(s, 'mlff_workflow.json')
"
check_success mlff_workflow.json


# --- 1. Errors before any round ---
echo -e "\n--- Testing: no training round yet ---"
stb-mlffValidate --no-intro > log_noround.txt 2>&1
check_exit_code $? 1
check_contains "No training rounds" log_noround.txt


# --- 2. Round 1 ---
echo -e "\n--- Training round 1 (2-model committee) and validating it ---"
stb-mlffAnalysis --committee 2 --epochs 2 --batch-size 2 --device cpu \
    --skip-foundation-comparison --no-intro > log_train1.txt 2>&1
check_exit_code $? 0
check_contains "Round directory   : round_01" log_train1.txt

stb-mlffValidate --save-report --no-intro > log_val1.txt 2>&1
check_exit_code $? 0
check_contains "Test set from     : mlff_workflow.json (test set)" log_val1.txt
check_contains "2 finished test config(s) (0 skipped): rattle 1, strain 1" log_val1.txt
check_contains "Round 1 (2 model(s), committee mean)" log_val1.txt
check_contains "F RMSE / MAE" log_val1.txt
check_contains "Members F RMSE" log_val1.txt
check_contains "Committee spread" log_val1.txt
check_contains "Foundation (small, no fine-tuning)" log_val1.txt
check_contains "\[CONTINUE\] First scored round" log_val1.txt
check_contains "Fewer than 5 test configs" log_val1.txt
check_success mlff_validation.png
check_success stb_mlffValidate_report.txt
check_contains "STB-MLFFVALIDATE REPORT" stb_mlffValidate_report.txt
check_python "round 1 validation stored in mlff_workflow.json" "
import json, sys
d = json.load(open('mlff_workflow.json'))
v = d['rounds'][0]['validation']
sys.exit(0 if v['n_test'] == 2 and v['overall']['F_rmse_meV'] > 0 and len(v['members']) == 2
         and set(v['groups']) == {'rattle', 'strain'} and 'calibration' in v and v['verdict'] == 'continue'
         and d['history'][-1]['stage'] == 'validate' else 1)
"


# --- 3. Round 2: two rounds compared, verdicts ---
echo -e "\n--- Training round 2 and comparing the rounds ---"
stb-mlffAnalysis --committee 2 --epochs 3 --batch-size 2 --device cpu \
    --skip-foundation-comparison --no-intro > log_train2.txt 2>&1
check_exit_code $? 0
stb-mlffValidate --skip-foundation-comparison --no-intro > log_val2.txt 2>&1
check_exit_code $? 0
check_contains "Round 2 (2 model(s), committee mean)" log_val2.txt
check_python "verdict table lists rounds 1 and 2, then a CONTINUE or PLATEAU verdict" "
import re, sys
t = open('log_val2.txt').read().split('[3] VERDICT', 1)[1]
rows = re.findall(r'^\s+(\d+)\s+[\d.]+\s+[\d.]+', t, re.M)
sys.exit(0 if rows[:2] == ['1', '2'] and re.search(r'\[(CONTINUE|PLATEAU)\] Round 2', t) else 1)
"

echo "Testing: --target-force reached -> STOP"
stb-mlffValidate --target-force 1e9 --skip-foundation-comparison --no-intro > log_val_stop.txt 2>&1
check_contains "\[STOP\] Target met" log_val_stop.txt

echo "Testing: --plateau 1000 -> any change counts as a plateau"
stb-mlffValidate --plateau 1000 --skip-foundation-comparison --no-intro > log_val_plateau.txt 2>&1
check_contains "\[PLATEAU\] Round 2" log_val_plateau.txt

echo "Testing: --rounds 1 scores only round 1"
stb-mlffValidate --rounds 1 --skip-foundation-comparison --no-intro > log_val_r1.txt 2>&1
check_contains "Scoring           : round 1 (2 model(s))" log_val_r1.txt
grep -q "Round 2 (" log_val_r1.txt
if [ $? -ne 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} round 2 not scored"; PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} round 2 scored despite --rounds 1"; FAIL=$((FAIL+1))
fi
stb-mlffValidate --rounds 7 --no-intro > log_val_r7.txt 2>&1
check_exit_code $? 1
check_contains "Round 7 is not in" log_val_r7.txt


# --- 4. Explicit models / test folders, no workflow ---
echo -e "\n--- Testing --model + --test-path (no workflow file) ---"
M1=$(python3 -c "import json; print(json.load(open('mlff_workflow.json'))['rounds'][0]['committee'][0]['model'])")
stb-mlffValidate --workflow none.json --model "$M1" --test-path 'mlff_test/test_*' \
    --skip-foundation-comparison --no-intro > log_val_explicit.txt 2>&1
check_exit_code $? 0
check_contains "1 model(s) from --model" log_val_explicit.txt
check_contains "2 finished test config(s) (0 skipped): test 2" log_val_explicit.txt
check_contains "no round history to compare" log_val_explicit.txt

echo "Testing: test folders that were also trained on are flagged"
stb-mlffValidate --test-path 'mlff_config_00[1-2]' --rounds 1 --skip-foundation-comparison --no-intro > log_val_overlap.txt 2>&1
check_contains "were also trained on" log_val_overlap.txt

echo "Testing: no test set -> clear error"
stb-mlffValidate --test-path 'no_such_*' --no-intro > log_val_notest.txt 2>&1
check_exit_code $? 1
check_contains "No test configs found" log_val_notest.txt

echo "Testing: --phonon-ref without phonopy_disp.yaml -> clear error"
stb-mlffValidate --phonon-ref round_01 --skip-foundation-comparison --no-intro > log_val_phref.txt 2>&1
check_exit_code $? 1
check_contains "has no phonopy_disp.yaml" log_val_phref.txt


# --- 5. Help / version ---
echo -e "\n--- Testing --help / --version ---"
stb-mlffValidate --version > log_version.txt 2>&1
check_contains "stb-mlffValidate" log_version.txt
stb-mlffValidate --help > log_help.txt 2>&1
for opt in workflow rounds model test-path phonon-ref phonon-mesh target-force plateau; do
    check_contains "$opt" log_help.txt
done


# --- 6. Interactive path (stb-suite, shortcut 4.17.3) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.17.3) ---"
printf '4.17.3\n\n\n\nn\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_contains "mlff_workflow.json: 6 training + 2 test config(s), 2 round(s)" log_menu.txt
check_contains "VERDICT" log_menu.txt


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
