#!/bin/bash

# --- Setup ---
# Smoke test for stb-mlmelting (ML Melting Point, item 5.5)
#
# Needs the optional 'ml' extra (pip install stb_suite[ml] -- PyTorch +
# mace-torch). The whole file is skipped with a clear message if `mace`
# isn't importable, same gating pattern as the other test/5-mlsimulations/*
# tools.
#
# al4.fdf: 4-atom bulk Al (fcc conventional cubic cell, a=4.05 Ang). Very
# short --equilibration-steps/--production-steps keep this a fast smoke
# test, not a converged melting-point calculation -- verified live with
# fuller step counts (200/500) that the Lindemann index rises monotonically
# with temperature and correctly crosses the default 0.10 threshold
# (estimated ~3095 K, vastly overestimating real Al's ~933 K -- the well-
# known "one-phase" superheating artifact for a tiny, defect-free periodic
# cell with no free surface to nucleate melting from, not a code bug).
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


echo "--- Starting tester for stb-mlmelting (item 5.5) ---"

if ! python3 -c "import mace" 2>/dev/null; then
    echo -e "${YELLOW}Skipped entirely:${NC} the optional 'ml' extra is not installed."
    echo "Install with: pip install stb_suite[ml]  (then re-run this test)"
    exit 0
fi

rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/al4.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 1. Basic scan (fast: short equilibration/production) ---
echo -e "\n--- Testing a basic temperature scan (4 points, short steps) ---"
stb-mlmelting --file al4.fdf --temperatures 300 1500 3000 4500 \
    --equilibration-steps 30 --production-steps 50 --stride 5 --seed 3 \
    --save-data --save-report --no-intro -o al_melt > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "TEMPERATURE SCAN" log_basic.txt
check_contains "MELTING-POINT ESTIMATE" log_basic.txt
check_success al_melt/melting_curve.png
check_success al_melt/melting_curve.dat
check_success al_melt/stb_mlmelting_report.txt
check_contains "STB-MLMELTING REPORT" al_melt/stb_mlmelting_report.txt

echo "Testing: Lindemann index increases monotonically with temperature"
python3 -c "
import numpy as np, sys
data = np.loadtxt('al_melt/melting_curve.dat')
lindemann = data[:, 1]
sys.exit(0 if np.all(np.diff(lindemann) > 0) else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} Lindemann index rises monotonically with T"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} Lindemann index did not rise monotonically with T"
    FAIL=$((FAIL+1))
fi

echo "Testing: diffusion coefficient is never negative (clipped)"
python3 -c "
import numpy as np, sys
data = np.loadtxt('al_melt/melting_curve.dat')
d = data[:, 2]
sys.exit(0 if np.all(d >= 0.0) else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} diffusion coefficient is never negative"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} diffusion coefficient is negative (should be clipped to 0)"
    FAIL=$((FAIL+1))
fi
rm -rf al_melt


# --- 2. Melting point estimate (lower threshold guarantees a crossing) ---
echo -e "\n--- Testing melting-point crossing detection (--lindemann-threshold 0.06) ---"
stb-mlmelting --file al4.fdf --temperatures 300 1500 3000 4500 \
    --equilibration-steps 30 --production-steps 50 --stride 5 --seed 3 \
    --lindemann-threshold 0.06 --no-intro -o al_melt_cross > log_cross.txt 2>&1
check_exit_code $? 0
check_contains "Estimated melting point" log_cross.txt
check_contains "SUPERHEAT" log_cross.txt
rm -rf al_melt_cross


# --- 3. Bulk-only check (rejects vacuum-padded structures) ---
echo -e "\n--- Testing bulk-only rejection (graphene, vacuum along c) ---"
cp "$FIXTURE_DIR/../3-mlelastic/graphene.fdf" .
stb-mlmelting --file graphene.fdf --no-intro > log_vacuum.txt 2>&1
check_exit_code $? 1
check_contains "bulk (3D periodic) only" log_vacuum.txt
rm -f graphene.fdf


# --- 4. Error and robustness cases ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: missing --file"
stb-mlmelting --no-intro > log_missing_file.txt 2>&1
check_exit_code $? 2

echo "Testing: nonexistent input file"
stb-mlmelting --file nope.fdf --no-intro > log_missing_input.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_input.txt

echo "Testing: --custom-model with a nonexistent file"
stb-mlmelting --file al4.fdf --custom-model does_not_exist.model --no-intro > log_custommodel.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_custommodel.txt

echo "Testing: a single temperature is rejected (needs at least 2)"
stb-mlmelting --file al4.fdf --temperatures 500 --no-intro > log_singletemp.txt 2>&1
check_exit_code $? 2
check_contains "at least 2 temperatures" log_singletemp.txt

echo "Testing: --version"
stb-mlmelting --version > log_version.txt 2>&1
check_contains "stb-mlmelting" log_version.txt

echo "Testing: --help documents --temp-min, --temp-max, --temperatures, --lindemann-threshold"
stb-mlmelting --help > log_help.txt 2>&1
check_contains "temp-min" log_help.txt
check_contains "temp-max" log_help.txt
check_contains "temperatures" log_help.txt
check_contains "lindemann-threshold" log_help.txt
check_contains "custom-model" log_help.txt
check_contains "save-trajectories" log_help.txt


# --- 5. Interactive path (stb-suite, shortcut 5.5) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 5.5) ---"
rm -rf interactive_out
printf '5.5\nal4.fdf\n\nsmall\n1\n300\n1800\n1500\n30\n50\ninteractive_out\nn\nn\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success interactive_out/melting_curve.png


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
