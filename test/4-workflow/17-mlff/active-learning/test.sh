#!/bin/bash

# --- Setup ---
# Smoke test for stb-mlffActiveLearning (Custom ML Force Field Active
# Learning, item 4.17.3)
#
# Reuses the 8 real finished-SIESTA fixtures from ../analysis/ to quickly
# fine-tune a tiny model on the fly (no large binary .model committed to
# the repo) via stb-mlffAnalysis, then feeds that model into
# stb-mlffActiveLearning to confirm the whole screening/ranking/writing
# pipeline works against a genuine model. Needs the optional 'ml' extra --
# skips gracefully (not a failure) if mace_run_train isn't on PATH.
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


echo "--- Starting tester for stb-mlffActiveLearning (item 4.17.3) ---"

if ! command -v mace_run_train &> /dev/null; then
    echo -e "${YELLOW}SKIPPED${NC}: 'mace_run_train' not on PATH (needs 'pip install stb_suite[ml]')."
    exit 0
fi
if [ ! -d "$MLFF_FIXTURES" ]; then
    echo -e "${YELLOW}SKIPPED${NC}: ../analysis fixtures not found."
    exit 0
fi

rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp -r "$MLFF_FIXTURES"/mlff_config_* "$TEST_DIR/"
# The analysis fixtures deliberately have no pseudopotentials (stb-mlffAnalysis
# never needs them) -- add a placeholder one for stb-mlffActiveLearning's
# --pseudo-dir (default "."), content doesn't matter (only copied, never read).
echo "dummy pseudopotential content" > "$TEST_DIR/O.psf"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null

echo -e "\n--- Fine-tuning a quick model to screen with (--epochs 2) ---"
stb-mlffAnalysis --path "mlff_config_*" --epochs 2 --batch-size 4 --device cpu \
    --name quick_model --skip-foundation-comparison --no-intro > log_quicktrain.txt 2>&1
QUICK_MODEL=$(find . -maxdepth 1 -name "quick_model*.model" | sort | tail -1)
if [ -z "$QUICK_MODEL" ]; then
    echo -e "${RED}FAIL${NC}: could not produce a quick fine-tuned model to test against."
    popd > /dev/null
    exit 1
fi
echo "Using model: $QUICK_MODEL"


# --- 1. Basic screening run ---
echo -e "\n--- Testing a basic screening run (--n-candidates 15 --top-k 3) ---"
stb-mlffActiveLearning --model "$QUICK_MODEL" --file mlff_config_001/calc.fdf --pseudo-dir . \
    --n-candidates 15 --stdev 0.2 --top-k 3 --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "Screened 15/15 candidate(s)" log_basic.txt
check_success mlff_config_009/calc.fdf
check_success mlff_config_010/calc.fdf
check_success mlff_config_011/calc.fdf
check_success mlff_config_009/O.psf

echo "Testing: numbering continues past the 8 already-present mlff_config_* folders"
check_contains "mlff_config_009" log_basic.txt

echo "Testing: suggested configs are ranked by decreasing predicted max|F|"
python3 -c "
import re, sys
text = open('log_basic.txt').read()
scores = [float(x) for x in re.findall(r'max\|F\| = ([\d.]+) eV/Ang\)', text)]
sys.exit(0 if scores == sorted(scores, reverse=True) and len(scores) == 3 else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} 3 suggested configs, ranked by decreasing predicted force"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} suggested configs not correctly ranked"
    FAIL=$((FAIL+1))
fi

echo "Testing: each suggested config keeps the reference basis/SCF settings"
check_contains "PAO.BasisSize       SZ" mlff_config_009/calc.fdf


# --- 2. --save-report ---
echo -e "\n--- Testing --save-report ---"
rm -f stb_mlffActiveLearning_report.txt
stb-mlffActiveLearning --model "$QUICK_MODEL" --file mlff_config_001/calc.fdf --pseudo-dir . \
    --n-candidates 5 --top-k 1 --output-dir report_test --save-report --no-intro > log_savereport.txt 2>&1
check_success stb_mlffActiveLearning_report.txt
check_contains "STB-MLFFACTIVELEARNING REPORT" stb_mlffActiveLearning_report.txt
rm -f stb_mlffActiveLearning_report.txt
rm -rf report_test


# --- 3. Robustness / errors ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: missing --model / --file"
stb-mlffActiveLearning --no-intro > log_missing_args.txt 2>&1
check_exit_code $? 2

echo "Testing: nonexistent model file"
stb-mlffActiveLearning --model does_not_exist.model --file mlff_config_001/calc.fdf --no-intro > log_missing_model.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_model.txt

echo "Testing: --version"
stb-mlffActiveLearning --version > log_version.txt 2>&1
check_contains "stb-mlffActiveLearning" log_version.txt

echo "Testing: --help documents --n-candidates, --top-k, --stdev"
stb-mlffActiveLearning --help > log_help.txt 2>&1
check_contains "n-candidates" log_help.txt
check_contains "top-k" log_help.txt
check_contains "stdev" log_help.txt


# --- 4. Interactive path (stb-suite, shortcut 4.17.3) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.17.3) ---"
rm -rf mlff_config_012 mlff_config_013
printf '4.17.3\n%s\nmlff_config_001/calc.fdf\n\n5\n0.2\n1\n.\nn\n\n0\n' "$QUICK_MODEL" | stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_contains "Suggested configs" log_menu.txt


popd > /dev/null

# --- 5. Summary ---
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
