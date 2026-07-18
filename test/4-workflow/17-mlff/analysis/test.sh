#!/bin/bash

# --- Setup ---
# Smoke test for stb-mlffAnalysis (Custom ML Force Field Analysis, item 4.17.2)
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
stb-mlffAnalysis --path "mlff_config_*" --epochs 5 --batch-size 4 \
    --device cpu --name test_model --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "8 usable configuration(s) (0 skipped)" log_basic.txt
check_contains "mace_run_train finished" log_basic.txt
check_contains "Fine-tuned model" log_basic.txt
check_contains "Model loads and runs" log_basic.txt
check_success training_set.xyz
check_success mace_run_train.log

echo "Testing: training_set.xyz has the right per-atom/info fields"
python3 -c "
import sys
from ase.io import read
frames = read('training_set.xyz', index=':')
sys.exit(0 if len(frames) == 8 and all(
    'REF_energy' in a.info and 'REF_forces' in a.arrays for a in frames
) else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} training_set.xyz has 8 frames with REF_energy/REF_forces"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} training_set.xyz missing expected frames/fields"
    FAIL=$((FAIL+1))
fi

echo "Testing: E0s uses 'average' (SIESTA-scale), never 'foundation' (wrong absolute scale)"
check_contains "E0s average" log_basic.txt


# --- 1b. --export-lammps ---
echo -e "\n--- Testing --export-lammps ---"
if command -v mace_create_lammps_model &> /dev/null; then
    stb-mlffAnalysis --path "mlff_config_*" --epochs 3 --batch-size 4 --device cpu \
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
(cd too_few && stb-mlffAnalysis --path "mlff_config_*" --no-intro > ../log_toofew.txt 2>&1)
check_exit_code $? 1
check_contains "Only" log_toofew.txt
rm -rf too_few


# --- 3. No matching folders ---
echo -e "\n--- Testing --path matching nothing (expects failure) ---"
stb-mlffAnalysis --path "no_such_folder_*" --no-intro > log_nomatch.txt 2>&1
check_exit_code $? 1
check_contains "No directories matched" log_nomatch.txt


# --- 4. --save-report ---
echo -e "\n--- Testing --save-report ---"
rm -f stb_mlffAnalysis_report.txt
stb-mlffAnalysis --path "mlff_config_*" --epochs 2 --batch-size 4 --device cpu \
    --name test_model2 --save-report --no-intro > log_savereport.txt 2>&1
check_success stb_mlffAnalysis_report.txt
check_contains "STB-MLFFANALYSIS REPORT" stb_mlffAnalysis_report.txt
rm -f stb_mlffAnalysis_report.txt


# --- 5. Robustness / errors ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: --version"
stb-mlffAnalysis --version > log_version.txt 2>&1
check_contains "stb-mlffAnalysis" log_version.txt

echo "Testing: --help documents --foundation-model, --epochs, --valid-fraction, --export-lammps"
stb-mlffAnalysis --help > log_help.txt 2>&1
check_contains "foundation-model" log_help.txt
check_contains "epochs" log_help.txt
check_contains "valid-fraction" log_help.txt
check_contains "export-lammps" log_help.txt


# --- 6. Interactive path (stb-suite, shortcut 4.17.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.17.2) ---"
printf '4.17.2\nmlff_config_*\n1\n3\n4\ncpu\ntest_model3\n.\nn\nn\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_contains "Fine-tuned model" log_menu.txt


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
