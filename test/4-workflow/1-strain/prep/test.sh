#!/bin/bash

# --- Setup ---
# Smoke test for stb-strain (Stress-Strain Prep, item 4.1.1)
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
echo "--- Starting tester for STB-Strain prep (item 4.1.1) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR/subdir"
cp "$FIXTURE_DIR/structure.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/structure.fdf" "$TEST_DIR/subdir/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Uniaxial strain (x), written inside the strain_runs/ wrapper dir ---
echo -e "\n--- Testing --stdir x (uniaxial) ---"
stb-strain --file structure.fdf --stdir x --stmin 0 --stmax 2 --step 2 --no-intro > log_uniaxial.txt 2>&1
check_contains "Detected strain type: uniaxial" log_uniaxial.txt
check_contains "Detected dimensionality: 2D" log_uniaxial.txt
check_success strain_runs/strain_x_0.00/structure.fdf
check_success strain_runs/strain_x_2.00/structure.fdf
check_contains "6.026159992248" strain_runs/strain_x_2.00/structure.fdf
check_contains "5.116478079000" strain_runs/strain_x_2.00/structure.fdf


# --- 3. Biaxial strain (xy) ---
echo -e "\n--- Testing --stdir yx (biaxial, normalized to xy) ---"
stb-strain --file structure.fdf --stdir yx --stmin 5 --stmax 5 --step 1 --no-intro > log_biaxial.txt 2>&1
check_contains "Detected strain type: biaxial" log_biaxial.txt
check_contains "Strain direction: xy" log_biaxial.txt
check_success strain_runs/strain_xy_5.00/structure.fdf
check_contains "6.203399992020" strain_runs/strain_xy_5.00/structure.fdf
check_contains "5.372301982950" strain_runs/strain_xy_5.00/structure.fdf


# --- 4. Physics: refuse straining a vacuum-padded axis ---
echo -e "\n--- Testing vacuum-axis rejection (--stdir z, the vacuum axis of this 2D fixture) ---"
stb-strain --file structure.fdf --stdir z --stmin 0 --stmax 1 --no-intro > log_vacuum_z.txt 2>&1
check_exit_code $? 1
check_contains "vacuum-padded axis" log_vacuum_z.txt
check_contains "Periodic axis/axes available for this structure: x, y" log_vacuum_z.txt

echo "Testing: biaxial direction that includes the vacuum axis (xz)"
stb-strain --file structure.fdf --stdir xz --stmin 0 --stmax 1 --no-intro > log_vacuum_xz.txt 2>&1
check_exit_code $? 1
check_contains "vacuum-padded axis" log_vacuum_xz.txt


# --- 5. Bug fix: --file with a directory component ---
echo -e "\n--- Testing --file with a subdirectory in the path ---"
rm -rf strain_runs
stb-strain --file subdir/structure.fdf --stdir x --stmin 0 --stmax 2 --step 2 --no-intro > log_subdir.txt 2>&1
check_exit_code $? 0
check_success strain_runs/strain_x_2.00/structure.fdf


# --- 6. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: --stmin > --stmax"
stb-strain --file structure.fdf --stdir x --stmin 10 --stmax 5 --no-intro > log_bad_range.txt 2>&1
check_exit_code $? 1
check_contains "\[ERROR\] Minimum strain cannot be greater than maximum strain" log_bad_range.txt

echo "Testing: invalid --stdir"
stb-strain --file structure.fdf --stdir w --no-intro > log_bad_dir.txt 2>&1
check_exit_code $? 1
check_contains "\[ERROR\] Invalid direction" log_bad_dir.txt

echo "Testing: --step 0 (must not crash with a raw ZeroDivisionError traceback)"
stb-strain --file structure.fdf --stdir x --step 0 --no-intro > log_zero_step.txt 2>&1
check_exit_code $? 1
check_contains "\[ERROR\] Step cannot be zero" log_zero_step.txt

echo "Testing: missing required args"
stb-strain --file structure.fdf --no-intro > log_missing_args.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-strain --version > log_version.txt 2>&1
check_contains "stb-strain" log_version.txt

echo "Testing: --help documents uniaxial/biaxial/output-dir/vacuum-gap"
stb-strain --help > log_help.txt 2>&1
check_contains "biaxial" log_help.txt
check_contains "output-dir" log_help.txt
check_contains "vacuum-gap" log_help.txt


# --- 7. Interactive path (stb-suite, shortcut 4.1.1) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.1.1) ---"

echo "Testing: navigate 4.1.1 -> invalid file then valid -> x -> 0/2/2 -> quit"
rm -rf strain_runs
printf '4.1.1\ndoes_not_exist.fdf\nstructure.fdf\nx\n0\n2\n2\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_contains "File not found" log_menu.txt
check_success strain_runs/strain_x_2.00/structure.fdf


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
