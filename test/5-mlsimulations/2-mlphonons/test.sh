#!/bin/bash

# --- Setup ---
# Smoke test for stb-mlphonons (ML Phonons, item 5.2)
#
# Needs the optional 'ml' extra (pip install stb_suite[ml] -- PyTorch +
# mace-torch). The whole file is skipped with a clear message if `mace`
# isn't importable, same gating pattern as test/2-structures/12-mlrelax and
# test/5-mlsimulations/1-mlmd. si_bulk.fdf is the same 64-atom Ge:Si
# substitutional-defect fixture used there (bulk, 3D periodic). --dim 1 1 1
# and a small --mesh/--band-points keep this a fast smoke test, not a
# converged production phonon calculation.
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


echo "--- Starting tester for stb-mlphonons (item 5.2) ---"

if ! python3 -c "import mace" 2>/dev/null; then
    echo -e "${YELLOW}Skipped entirely:${NC} the optional 'ml' extra is not installed."
    echo "Install with: pip install stb_suite[ml]  (then re-run this test)"
    exit 0
fi

rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/si_bulk.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 1. Basic run ---
echo -e "\n--- Testing a basic run (--dim 1 1 1, small mesh/band-points) ---"
stb-mlphonons --file si_bulk.fdf --dim 1 1 1 --mesh 8 8 8 --band-points 11 \
    --save-data --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "Displacements needed" log_basic.txt
check_contains "BAND STRUCTURE" log_basic.txt
check_contains "DENSITY OF STATES" log_basic.txt
check_contains "THERMAL PROPERTIES" log_basic.txt
check_success mlphonons_out/phonopy_disp.yaml
check_success mlphonons_out/bands.png
check_success mlphonons_out/bands.dat
check_success mlphonons_out/dos.png
check_success mlphonons_out/dos.dat
check_success mlphonons_out/thermal.png
check_success mlphonons_out/thermal.dat

echo "Testing: thermal properties data is finite (no NaN/Inf) and entropy starts near 0 (3rd law)"
python3 -c "
import numpy as np, sys
data = np.loadtxt('mlphonons_out/thermal.dat')
temps, free_e, entropy, cv = data[:,0], data[:,1], data[:,2], data[:,3]
finite = np.all(np.isfinite(data))
sys.exit(0 if finite and entropy[0] < 1.0 else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} thermal properties finite, entropy(T~0) near 0"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} thermal properties data missing/invalid"
    FAIL=$((FAIL+1))
fi

rm -rf mlphonons_out


# --- 2. --custom-model with a nonexistent file (expects failure) ---
echo -e "\n--- Testing --custom-model with a nonexistent file (expects failure) ---"
stb-mlphonons --file si_bulk.fdf --custom-model does_not_exist.model --no-intro > log_custommodel.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_custommodel.txt


# --- 3. --save-report ---
echo -e "\n--- Testing --save-report ---"
stb-mlphonons --file si_bulk.fdf --dim 1 1 1 --mesh 8 8 8 --band-points 11 \
    -o report_test --save-report --no-intro > log_savereport.txt 2>&1
check_exit_code $? 0
check_success report_test/stb_mlphonons_report.txt
check_contains "STB-MLPHONONS REPORT" report_test/stb_mlphonons_report.txt
rm -rf report_test


# --- 4. Robustness / errors ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: missing --file"
stb-mlphonons --no-intro > log_missing_file.txt 2>&1
check_exit_code $? 2

echo "Testing: nonexistent input file"
stb-mlphonons --file nope.fdf --no-intro > log_missing_input.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_input.txt

echo "Testing: --version"
stb-mlphonons --version > log_version.txt 2>&1
check_contains "stb-mlphonons" log_version.txt

echo "Testing: --help documents --dim, --custom-model, --mesh, --no-relax"
stb-mlphonons --help > log_help.txt 2>&1
check_contains "dim" log_help.txt
check_contains "custom-model" log_help.txt
check_contains "mesh" log_help.txt
check_contains "no-relax" log_help.txt


# --- 5. Interactive path (stb-suite, shortcut 5.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 5.2) ---"
rm -rf interactive_out
printf '5.2\nsi_bulk.fdf\n\nsmall\n1 1 1\ninteractive_out\nn\nn\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success interactive_out/bands.png
check_success interactive_out/dos.png
check_success interactive_out/thermal.png


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
