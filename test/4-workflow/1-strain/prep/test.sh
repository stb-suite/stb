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
cp "$FIXTURE_DIR/si_cubic.fdf" "$TEST_DIR/"
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


# --- 5. Symmetry advisory (uniaxial only, informational, never blocks) ---
echo -e "\n--- Testing the symmetry-equivalence advisory ---"

echo "Testing: cubic structure (m-3m) -- x/y/z are truly equivalent, advisory must fire"
rm -rf strain_runs
stb-strain --file si_cubic.fdf --stdir x --stmin 0 --stmax 1 --no-intro > log_sym_cubic.txt 2>&1
check_exit_code $? 0
check_contains "AXIS SYMMETRY (uniaxial direction 'x')" log_sym_cubic.txt
check_contains "point group m-3m" log_sym_cubic.txt
check_contains "48 operation" log_sym_cubic.txt
check_contains "y.*EQUIVALENT.*x" log_sym_cubic.txt
check_contains "z.*EQUIVALENT.*x" log_sym_cubic.txt
check_contains "y, z are equivalent to 'x' by symmetry" log_sym_cubic.txt

echo "Testing: hexagonal 2D fixture -- x and y are NOT equivalent (known zigzag/armchair"
echo "         asymmetry), advisory must NOT fire"
rm -rf strain_runs
stb-strain --file structure.fdf --stdir x --stmin 0 --stmax 1 --no-intro > log_sym_hex.txt 2>&1
check_exit_code $? 0
if grep -q "equivalent to 'x' by symmetry" log_sym_hex.txt 2>/dev/null; then
    echo -e "   -> ${RED}Failed:${NC} unexpected symmetry advisory for the hexagonal 2D fixture"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no symmetry advisory for x/y on the hexagonal 2D fixture"
    PASS=$((PASS+1))
fi

echo "Testing: biaxial direction -- out of scope for v1, advisory must NOT fire"
rm -rf strain_runs
stb-strain --file si_cubic.fdf --stdir xy --stmin 0 --stmax 1 --no-intro > log_sym_biaxial.txt 2>&1
check_exit_code $? 0
if grep -q "equivalent to" log_sym_biaxial.txt 2>/dev/null; then
    echo -e "   -> ${RED}Failed:${NC} unexpected symmetry advisory for a biaxial direction"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no symmetry advisory for biaxial directions (v1 scope)"
    PASS=$((PASS+1))
fi

echo "Testing: --help documents --symprec/--angle-tolerance"
stb-strain --help > log_help_sym.txt 2>&1
check_contains "symprec" log_help_sym.txt
check_contains "angle-tolerance" log_help_sym.txt


# --- 6. Bug fix: --file with a directory component ---
echo -e "\n--- Testing --file with a subdirectory in the path ---"
rm -rf strain_runs
stb-strain --file subdir/structure.fdf --stdir x --stmin 0 --stmax 2 --step 2 --no-intro > log_subdir.txt 2>&1
check_exit_code $? 0
check_success strain_runs/strain_x_2.00/structure.fdf


# --- 7. Error and robustness cases ---
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


# --- 8. Interactive path (stb-suite, shortcut 4.1.1) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.1.1) ---"

echo "Testing: navigate 4.1.1 -> invalid file then valid -> skip advanced -> x (no equivalent on this 2D fixture) -> 0/2/2 -> quit"
rm -rf strain_runs
# Prompts in order: file (retry), advanced settings (n -> skip), direction
# (x -- no symmetry-equivalent axis on this hexagonal 2D fixture, so no
# "generate equivalent(s) too?" follow-up), stmin, stmax, step, then the
# "Press Enter to continue" pause, then quit.
printf '4.1.1\ndoes_not_exist.fdf\nstructure.fdf\nn\nx\n0\n2\n2\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_contains "File not found" log_menu.txt
check_contains "Detected dimensionality: 2D" log_menu.txt
check_contains "AXIS SYMMETRY (uniaxial direction 'x')" log_menu.txt
check_contains "CONFIGURATION SUMMARY" log_menu.txt
check_contains "Dimensionality.*2D" log_menu.txt
check_success strain_runs/strain_x_2.00/structure.fdf

echo "Testing: navigate 4.1.1 -> biaxial direction -> dimensionality still detected/shown, no axis-symmetry table"
rm -rf strain_runs
printf '4.1.1\nstructure.fdf\nn\nxy\n0\n2\n2\n\n0\n' | stb-suite > log_menu_biaxial.txt 2>&1
check_contains "Detected dimensionality: 2D" log_menu_biaxial.txt
check_contains "Dimensionality.*2D" log_menu_biaxial.txt
if grep -q "AXIS SYMMETRY" log_menu_biaxial.txt 2>/dev/null; then
    echo -e "   -> ${RED}Failed:${NC} axis-symmetry table shown for a biaxial direction (out of scope)"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no axis-symmetry table for a biaxial direction"
    PASS=$((PASS+1))
fi
check_success strain_runs/strain_xy_2.00/structure.fdf

echo "Testing: navigate 4.1.1 -> cubic fixture -> x -> accept generating the equivalent y,z folders too"
rm -rf strain_runs
# Prompts in order: file, advanced settings (n), direction (x -- y,z ARE
# equivalent on this cubic fixture), accept the follow-up (y), stmin, stmax,
# step, then one "Press Enter to continue" pause (not three -- run_tool's
# pause is suppressed for every direction but the last), then quit.
printf '4.1.1\nsi_cubic.fdf\nn\nx\ny\n0\n2\n2\n\n0\n' | stb-suite > log_menu_multi.txt 2>&1
check_contains "y, z are equivalent to 'x' by symmetry" log_menu_multi.txt
check_contains "Direction(s) to generate.*x, y, z" log_menu_multi.txt
check_success strain_runs/strain_x_2.00/si_cubic.fdf
check_success strain_runs/strain_y_2.00/si_cubic.fdf
check_success strain_runs/strain_z_2.00/si_cubic.fdf
if [ "$(grep -c "Press Enter to continue" log_menu_multi.txt)" -eq 1 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} only one 'Press Enter to continue' pause for all 3 directions"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} expected exactly 1 'Press Enter to continue' pause, "\
"got $(grep -c "Press Enter to continue" log_menu_multi.txt)"
    FAIL=$((FAIL+1))
fi


popd > /dev/null

# --- 9. Summary ---
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
