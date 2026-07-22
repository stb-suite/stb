#!/bin/bash

# --- Setup ---
# Smoke test for stb-mlconvergence (ML Model-Size Convergence, item 5.6)
#
# Needs the optional 'ml' extra (pip install stb_suite[ml] -- PyTorch +
# mace-torch). The whole file is skipped with a clear message if `mace`
# isn't importable, same gating pattern as the other test/5-mlsimulations/*
# tools.
#
# si8.fdf: 8-atom bulk Si (diamond cubic, a=5.43 Ang), same fixture as
# test/5-mlsimulations/3-mlelastic. Verified live: small vs. medium MACE-MP-0
# differ by ~29 meV/atom and ~0.6% in cell volume for this structure --
# correctly flagged NOT CONVERGED at the tight default tolerances (5 meV/
# atom, 1%), and correctly flagged OK at looser ones.
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


echo "--- Starting tester for stb-mlconvergence (item 5.6) ---"

if ! python3 -c "import mace" 2>/dev/null; then
    echo -e "${YELLOW}Skipped entirely:${NC} the optional 'ml' extra is not installed."
    echo "Install with: pip install stb_suite[ml]  (then re-run this test)"
    exit 0
fi

rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/si8.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 1. small vs medium, tight tolerances (expects NOT CONVERGED) ---
echo -e "\n--- Testing small vs medium at tight (default) tolerances ---"
stb-mlconvergence --file si8.fdf --models small medium \
    --save-data --save-report --no-intro -o conv_tight > log_tight.txt 2>&1
check_exit_code $? 0
check_contains "PER-MODEL CALCULATION" log_tight.txt
check_contains "CONVERGENCE COMPARISON" log_tight.txt
check_contains "small -> medium" log_tight.txt
check_contains "NOT CONVERGED" log_tight.txt
check_success conv_tight/convergence.png
check_success conv_tight/convergence.dat
check_success conv_tight/stb_mlconvergence_report.txt
check_contains "STB-MLCONVERGENCE REPORT" conv_tight/stb_mlconvergence_report.txt
rm -rf conv_tight


# --- 2. small vs medium, loose tolerances (expects OK + a recommendation) ---
echo -e "\n--- Testing small vs medium at loose tolerances ---"
stb-mlconvergence --file si8.fdf --models small medium --energy-tolerance 50 \
    --lattice-tolerance 2 --no-intro -o conv_loose > log_loose.txt 2>&1
check_exit_code $? 0
check_contains "All consecutive pairs agree within tolerance" log_loose.txt
check_contains "Recommended       : small" log_loose.txt
rm -rf conv_loose


# --- 3. Convergence data is internally consistent (energies/volumes present for both models) ---
echo -e "\n--- Testing convergence.dat has one row per model ---"
stb-mlconvergence --file si8.fdf --models small medium large --save-data \
    --no-intro -o conv_all > log_all.txt 2>&1
check_exit_code $? 0
python3 -c "
import sys
with open('conv_all/convergence.dat') as f:
    lines = [l for l in f if not l.startswith('#')]
sys.exit(0 if len(lines) == 3 else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} convergence.dat has one row per compared model (3)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} convergence.dat does not have 3 rows"
    FAIL=$((FAIL+1))
fi
rm -rf conv_all


# --- 4. Error and robustness cases ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: missing --file"
stb-mlconvergence --no-intro > log_missing_file.txt 2>&1
check_exit_code $? 2

echo "Testing: nonexistent input file"
stb-mlconvergence --file nope.fdf --no-intro > log_missing_input.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_input.txt

echo "Testing: --custom-models with a nonexistent file"
stb-mlconvergence --file si8.fdf --custom-models does_not_exist.model --no-intro > log_custommodel.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_custommodel.txt

echo "Testing: --version"
stb-mlconvergence --version > log_version.txt 2>&1
check_contains "stb-mlconvergence" log_version.txt

echo "Testing: --help documents --models, --custom-models, --energy-tolerance, --lattice-tolerance"
stb-mlconvergence --help > log_help.txt 2>&1
check_contains "\-\-models" log_help.txt
check_contains "custom-models" log_help.txt
check_contains "energy-tolerance" log_help.txt
check_contains "lattice-tolerance" log_help.txt
check_contains "no-relax" log_help.txt


# --- 5. Interactive path (stb-suite, shortcut 5.6) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 5.6) ---"
rm -rf interactive_out
printf '5.6\nsi8.fdf\nsmall medium\n\n\n\ninteractive_out\nn\nn\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success interactive_out/convergence.png


popd > /dev/null

# --- 6. Summary ---
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
