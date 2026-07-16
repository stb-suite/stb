#!/bin/bash

# --- Setup ---
# Smoke test for stb-raman (Raman Spectrum, Stage 1: Phonon Displacements, item 4.11.1)
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
echo "--- Starting tester for STB-Raman stage 1: phonon displacements (item 4.11.1) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/structure.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/calc.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/O.psf" "$TEST_DIR/"
cp "$FIXTURE_DIR/Sn.psf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Default run (dim 1 1 1, small enough to be fast) ---
echo -e "\n--- Testing default generation (-dim 1 1 1) ---"
stb-raman -s structure.fdf -c calc.fdf -p . -dim 1 1 1 --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_success raman_study/phonon_disp/disp-001/structure.fdf
check_success raman_study/phonon_disp/disp-001/calc.fdf
check_success raman_study/phonon_disp/disp-001/O.psf
check_success raman_study/phonon_disp/disp-001/Sn.psf
check_success raman_study/phonon_disp/phonopy_disp.yaml
check_success raman_study/raman_stage1.txt
check_contains "\[0\] RUN METADATA" raman_study/raman_stage1.txt
check_contains "\[1\] PSEUDOPOTENTIALS" raman_study/raman_stage1.txt
check_contains "\[2\] PHONON DISPLACEMENTS" raman_study/raman_stage1.txt
check_contains "\[3\] SUMMARY & NEXT STEPS" raman_study/raman_stage1.txt
check_contains "Displacement folders : 84" log_basic.txt
check_contains "stb-ramanModes --directory raman_study" raman_study/raman_stage1.txt


# --- 3. Refuses to overwrite an existing phonon_disp/ ---
echo -e "\n--- Testing the re-run guard (existing phonon_disp/) ---"
stb-raman -s structure.fdf -c calc.fdf -p . -dim 1 1 1 --no-intro > log_regen.txt 2>&1
check_exit_code $? 1
check_contains "already contains" log_regen.txt


# --- 4. Error cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: missing structure file"
rm -rf raman_study
stb-raman -s does_not_exist.fdf -c calc.fdf -p . --no-intro > log_missing_struct.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_struct.txt

echo "Testing: missing calc file"
stb-raman -s structure.fdf -c does_not_exist.fdf -p . --no-intro > log_missing_calc.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_calc.txt

echo "Testing: missing pseudopotentials"
mkdir -p empty_pp
stb-raman -s structure.fdf -c calc.fdf -p empty_pp --no-intro > log_missing_pseudo.txt 2>&1
check_exit_code $? 1
check_contains "Missing pseudopotentials" log_missing_pseudo.txt

echo "Testing: --version"
stb-raman --version > log_version.txt 2>&1
check_contains "stb-raman" log_version.txt

echo "Testing: --help documents -dim/-d/-p/-O"
stb-raman --help > log_help.txt 2>&1
check_contains "dim" log_help.txt
check_contains "distance" log_help.txt
check_contains "pseudo" log_help.txt
check_contains "output-dir" log_help.txt


# --- 5. Interactive path (stb-suite, shortcut 4.11.1) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.11.1) ---"

echo "Testing: navigate 4.11.1 -> defaults except dim 1 1 1 -> quit"
rm -rf raman_study
{
  echo "4.11.1"
  echo "structure.fdf"  # structure_file
  echo "calc.fdf"       # calc_file
  echo "1 1 1"          # dim
  echo ""               # distance (default 0.02)
  echo ""               # pseudo_dir (skip -> defaults to ".", fixtures are in cwd)
  echo ""               # output_dir (default raman_study)
  echo ""               # show_advanced (default -> skip)
  echo ""               # press enter to continue
  echo "0"              # quit stage submenu
} | stb-suite > log_menu.txt 2>&1
check_contains "Displacement folders : 84" log_menu.txt
check_success raman_study/phonon_disp/disp-001/structure.fdf


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
