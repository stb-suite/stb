#!/bin/bash

# --- Setup ---
# Smoke test for stb-eosInputs (Equation of State Prep, item 4.18.1)
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


# --- 1. Preparation ---
echo "--- Starting tester for stb-eosInputs (item 4.18.1) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/si_cubic.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/structure.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Basic scan (bulk cubic Si, default settings) ---
echo -e "\n--- Testing a basic scan (bulk Si, default settings) ---"
stb-eosInputs --file si_cubic.fdf --no-intro -o eos_runs > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "Detected dimensionality: 3D" log_basic.txt
check_success eos_runs/vol_0.00/si_cubic.fdf
check_success eos_runs/vol_5.00/si_cubic.fdf
check_success eos_runs/vol_m5.00/si_cubic.fdf

echo "Testing: vol_0.00's raw lattice block is unchanged (zero strain, factor=1)"
check_contains "1.000000000000        0.000000000000        0.000000000000" eos_runs/vol_0.00/si_cubic.fdf

echo "Testing: vol_5.00's lattice is scaled by (1.05)^(1/3)"
check_contains "1.016396" eos_runs/vol_5.00/si_cubic.fdf

echo "Testing: atomic positions (fractional) are preserved verbatim"
check_contains "0.250000000   0.250000000   0.250000000   1" eos_runs/vol_5.00/si_cubic.fdf

echo "Testing: LatticeConstant line is untouched"
check_contains "LatticeConstant 5.43 Ang" eos_runs/vol_5.00/si_cubic.fdf
rm -rf eos_runs


# --- 3. Custom --n-volumes/--strain-range ---
echo -e "\n--- Testing --n-volumes 5 --strain-range 10.0 ---"
stb-eosInputs --file si_cubic.fdf --n-volumes 5 --strain-range 10.0 --no-intro \
    -o eos_custom > log_custom.txt 2>&1
check_exit_code $? 0
check_success eos_custom/vol_10.00/si_cubic.fdf
check_success eos_custom/vol_5.00/si_cubic.fdf
check_success eos_custom/vol_0.00/si_cubic.fdf
check_success eos_custom/vol_m5.00/si_cubic.fdf
check_success eos_custom/vol_m10.00/si_cubic.fdf
if [ -d eos_custom/vol_2.00 ]; then
    echo -e "   -> ${RED}Failed:${NC} unexpected vol_2.00 folder for a 5-point (-10..10) scan"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} exactly 5 points scanned (no off-grid vol_2.00 folder)"
    PASS=$((PASS+1))
fi
rm -rf eos_custom


# --- 4. Physics: refuse a vacuum-padded structure ---
echo -e "\n--- Testing vacuum-padded (2D) structure rejection ---"
stb-eosInputs --file structure.fdf --no-intro -o eos_2d > log_2d.txt 2>&1
check_exit_code $? 1
check_contains "bulk (3D periodic) only" log_2d.txt


# --- 5. --n-volumes below the minimum is rejected ---
echo -e "\n--- Testing --n-volumes 3 is rejected (needs >= 5) ---"
stb-eosInputs --file si_cubic.fdf --n-volumes 3 --no-intro -o eos_toofew > log_toofew.txt 2>&1
check_exit_code $? 2


# --- 6. Error and robustness cases ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: missing --file"
stb-eosInputs --no-intro > log_missing_args.txt 2>&1
check_exit_code $? 2

echo "Testing: nonexistent input file"
stb-eosInputs --file nope.fdf --no-intro > log_missing_input.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_input.txt

echo "Testing: --version"
stb-eosInputs --version > log_version.txt 2>&1
check_contains "stb-eosInputs" log_version.txt

echo "Testing: --help documents --file, --n-volumes, --strain-range, --output-dir, --vacuum-gap"
stb-eosInputs --help > log_help.txt 2>&1
check_contains "\-\-file" log_help.txt
check_contains "n-volumes" log_help.txt
check_contains "strain-range" log_help.txt
check_contains "output-dir" log_help.txt
check_contains "vacuum-gap" log_help.txt


# --- 7. Interactive path (stb-suite, shortcut 4.18.1) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.18.1) ---"
rm -rf eos_interactive
printf '4.18.1\nsi_cubic.fdf\n5\n10.0\neos_interactive\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success eos_interactive/vol_10.00/si_cubic.fdf


popd > /dev/null

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
