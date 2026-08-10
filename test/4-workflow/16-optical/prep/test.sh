#!/bin/bash

# --- Setup ---
# Smoke test for stb-optical (Optical Properties Stage 1: Direction
# Folders, item 4.16.1). 2 fixtures: cubic.fdf (3D bulk) and mos2.fdf
# (2D, exercises the vacuum-dilution report note).
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

check_not_contains() {
    if grep -q -- "$1" "$2" 2>/dev/null; then
        echo -e "   -> ${RED}Failed:${NC} '$1' unexpectedly found in '$2'"
        FAIL=$((FAIL+1))
    else
        echo -e "   -> ${GREEN}Verified:${NC} '$1' NOT found in '$2' (as expected)"
        PASS=$((PASS+1))
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
echo "--- Starting tester for STB-OPTICAL stage 1: direction folders (item 4.16.1) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/cubic.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/mos2.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/calc.fdf" "$TEST_DIR/"
for sym in Na Mo S; do
    echo "# placeholder pseudopotential" > "$TEST_DIR/$sym.psf"
done
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Default run (3D fixture, --directions x y z default) ---
echo -e "\n--- Testing default run (cubic.fdf, all 3 directions) ---"
rm -rf optical_study_cubic
stb-optical -f cubic.fdf -c calc.fdf -p . -O optical_study_cubic --no-intro \
    > log_cubic.txt 2>&1
check_exit_code $? 0
check_success optical_study_cubic/dir_x/structure.fdf
check_success optical_study_cubic/dir_y/structure.fdf
check_success optical_study_cubic/dir_z/structure.fdf
check_success optical_study_cubic/optical_stage1.txt

echo "Testing: every folder forces MD.TypeOfRun CG + MD.Steps 0 (single point)"
check_contains "MD.TypeOfRun          CG" optical_study_cubic/dir_x/calc.fdf
check_contains "MD.Steps              0" optical_study_cubic/dir_x/calc.fdf

echo "Testing: OpticalCalculation/Optical.Broaden/Optical.Mesh/Optical.Vector blocks are correct per axis"
check_contains "OpticalCalculation T" optical_study_cubic/dir_x/calc.fdf
check_contains "Optical.Broaden        0.2 eV" optical_study_cubic/dir_x/calc.fdf
check_contains "10  10  10" optical_study_cubic/dir_x/calc.fdf
check_contains "1.0000  0.0000  0.0000" optical_study_cubic/dir_x/calc.fdf
check_contains "0.0000  1.0000  0.0000" optical_study_cubic/dir_y/calc.fdf
check_contains "0.0000  0.0000  1.0000" optical_study_cubic/dir_z/calc.fdf

echo "Testing: kgrid.MonkhorstPack left untouched, exactly as the template had it (no recompute)"
check_contains "kgrid.MonkhorstPack   \[1  1  1\]" optical_study_cubic/dir_x/calc.fdf

echo "Testing: 3D fixture report has no [KNOWN LIMITATION] note"
check_not_contains "KNOWN LIMITATION" optical_study_cubic/optical_stage1.txt


# --- 3. --optical-mesh/--optical-broaden/--optical-nbands propagate ---
echo -e "\n--- Testing --optical-mesh/--optical-broaden/--optical-nbands ---"
rm -rf optical_study_params
stb-optical -f cubic.fdf -c calc.fdf -p . -O optical_study_params --directions x \
    --optical-mesh 20 20 20 --optical-broaden 0.05 --optical-nbands 40 --no-intro \
    > log_params.txt 2>&1
check_exit_code $? 0
check_contains "20  20  20" optical_study_params/dir_x/calc.fdf
check_contains "Optical.Broaden        0.05 eV" optical_study_params/dir_x/calc.fdf
check_contains "Optical.NumberOfBands  40" optical_study_params/dir_x/calc.fdf


# --- 4. --directions restricts which folders are written ---
echo -e "\n--- Testing --directions z (only out-of-plane) ---"
rm -rf optical_study_z
stb-optical -f cubic.fdf -c calc.fdf -p . -O optical_study_z --directions z --no-intro \
    > log_z_only.txt 2>&1
check_exit_code $? 0
check_success optical_study_z/dir_z/structure.fdf
python3 -c "
import os
assert not os.path.isdir('optical_study_z/dir_x'), 'dir_x should not exist'
assert not os.path.isdir('optical_study_z/dir_y'), 'dir_y should not exist'
print('OK')
" > log_z_only_check.txt 2>&1
check_contains "OK" log_z_only_check.txt


# --- 5. 2D fixture: KNOWN LIMITATION note ---
echo -e "\n--- Testing 2D fixture (mos2.fdf) -- vacuum-dilution note ---"
rm -rf optical_study_mos2
stb-optical -f mos2.fdf -c calc.fdf -p . -O optical_study_mos2 --directions x z --no-intro \
    > log_mos2.txt 2>&1
check_exit_code $? 0
check_success optical_study_mos2/dir_x/structure.fdf
check_success optical_study_mos2/dir_z/structure.fdf
check_contains "KNOWN LIMITATION" optical_study_mos2/optical_stage1.txt
check_contains "2D (e.g., a slab or surface)" optical_study_mos2/optical_stage1.txt


# --- 6. Error cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: --version"
stb-optical --version > log_version.txt 2>&1
check_contains "stb-optical" log_version.txt

echo "Testing: --help documents --directions/--optical-mesh/--optical-broaden/--optical-nbands"
stb-optical --help > log_help.txt 2>&1
check_contains "directions" log_help.txt
check_contains "optical-mesh" log_help.txt
check_contains "optical-broaden" log_help.txt
check_contains "optical-nbands" log_help.txt

echo "Testing: missing structure file is rejected"
stb-optical -f nonexistent.fdf -c calc.fdf -p . -O optical_study_bad --no-intro \
    > log_missing_structure.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_structure.txt


# --- 7. Interactive path (stb-suite, shortcut 4.16.1) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.16.1) ---"

echo "Testing: navigate 4.16.1 -> cubic fixture -> defaults -> quit"
rm -rf optical_study_menu
{
  echo "4.16.1"
  echo "cubic.fdf"        # structure
  echo "calc.fdf"           # calc template
  echo "3"                   # pseudo_dir -> option 3 = Custom path
  echo "."                    # custom pseudo path
  echo ""                      # directions (default x y z)
  echo "optical_study_menu"     # output dir
  echo ""                        # advanced settings (default N)
  echo ""                         # press enter to continue
  echo "0"                         # quit stage submenu
} | stb-suite > log_menu.txt 2>&1
check_success optical_study_menu/dir_x/structure.fdf
check_success optical_study_menu/dir_y/structure.fdf
check_success optical_study_menu/dir_z/structure.fdf


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
