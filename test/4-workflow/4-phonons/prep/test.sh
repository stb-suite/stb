#!/bin/bash

# --- Setup ---
# Smoke test for stb-phononsCreate (Phonons Prep, item 4.4.1)
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$FIXTURE_DIR/test_files"

# Output colors
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

PASS=0
FAIL=0

# --- Check helpers ---

# Checks that file $1 exists and is not empty
check_success() {
    if [ -s "$1" ]; then
        echo -e " ... ${GREEN}OK${NC} (file '$1' created)"
        PASS=$((PASS+1))
    else
        echo -e " ... ${RED}FAIL${NC} (file '$1' was not created)"
        FAIL=$((FAIL+1))
    fi
}

# Checks that file $2 contains (grep -q) pattern $1
check_contains() {
    if grep -q "$1" "$2" 2>/dev/null; then
        echo -e "   -> ${GREEN}Verified:${NC} '$1' found in '$2'"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} '$1' NOT found in '$2'"
        FAIL=$((FAIL+1))
    fi
}

# Checks that $1 (actual exit code) equals $2 (expected exit code)
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
echo "--- Starting tester for STB-PhononsCreate prep (item 4.4.1) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/structure.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/calc.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/Sn.psf" "$TEST_DIR/"
cp "$FIXTURE_DIR/O.psf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Basic displacement generation (-dim 1 1 1 -- same 14-atom unit cell,
#     kept small so the ML-preview test below stays fast; dim isn't what's
#     under test here) ---
echo -e "\n--- Testing basic displacement generation ---"
stb-phononsCreate -s structure.fdf -c calc.fdf -p . -dim 1 1 1 --no-intro > log_create.txt 2>&1
check_contains "Displacement folders : 84" log_create.txt
check_success phonon_runs/phonopy_disp.yaml
check_success phonon_runs/disp-001/structure.fdf
check_success phonon_runs/disp-001/calc.fdf
check_success phonon_runs/disp-001/Sn.psf
check_success phonon_runs/disp-001/O.psf
check_success phonon_runs/disp-084/structure.fdf
check_success phonon_runs/phonon_prep_properties.txt
check_contains "RUN METADATA" phonon_runs/phonon_prep_properties.txt
check_contains "Displacements needed : 84" phonon_runs/phonon_prep_properties.txt


# --- 3. --ml-prerelax pre-flight check (MACE-MP-0 already cached locally --
#     fast, no network needed) ---
echo -e "\n--- Testing --ml-prerelax ---"
rm -rf phonon_runs
stb-phononsCreate -s structure.fdf -c calc.fdf -p . -dim 1 1 1 --ml-prerelax --no-intro > log_ml_prerelax.txt 2>&1
check_contains "ML PRE-FLIGHT CHECK" log_ml_prerelax.txt
check_contains "Max residual force on input structure" log_ml_prerelax.txt
check_success phonon_runs/phonopy_disp.yaml


# --- 4. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: missing structure file"
rm -rf phonon_runs
stb-phononsCreate -s does_not_exist.fdf -c calc.fdf -p . -dim 1 1 1 --no-intro > log_missing_struct.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_struct.txt

echo "Testing: missing calc file"
stb-phononsCreate -s structure.fdf -c does_not_exist.fdf -p . -dim 1 1 1 --no-intro > log_missing_calc.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_calc.txt

echo "Testing: missing pseudopotentials"
mkdir -p empty_pseudo_dir
stb-phononsCreate -s structure.fdf -c calc.fdf -p empty_pseudo_dir -dim 1 1 1 --no-intro > log_missing_pseudo.txt 2>&1
check_exit_code $? 1
check_contains "Missing pseudopotentials" log_missing_pseudo.txt

echo "Testing: re-running on top of an existing phonon_runs/ refuses to overwrite"
stb-phononsCreate -s structure.fdf -c calc.fdf -p . -dim 1 1 1 --no-intro > log_regen_1.txt 2>&1
stb-phononsCreate -s structure.fdf -c calc.fdf -p . -dim 1 1 1 --no-intro > log_regen_2.txt 2>&1
check_exit_code $? 1
check_contains "already contains" log_regen_2.txt
rm -rf phonon_runs

echo "Testing: --version"
stb-phononsCreate --version > log_version.txt 2>&1
check_contains "stb-phononsCreate" log_version.txt

echo "Testing: --help documents -dim/-p"
stb-phononsCreate --help > log_help.txt 2>&1
check_contains "dim" log_help.txt
check_contains "pseudo" log_help.txt


# --- 5. Interactive path (stb-suite, shortcut 4.4.1) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.4.1) ---"

echo "Testing: navigate 4.4.1 -> defaults except dim 1 1 1, skip ML -> quit"
rm -rf phonon_runs
printf '4.4.1\n\n\n1 1 1\n\n\nn\nn\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_contains "Displacement folders : 84" log_menu.txt
check_success phonon_runs/phonopy_disp.yaml


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
