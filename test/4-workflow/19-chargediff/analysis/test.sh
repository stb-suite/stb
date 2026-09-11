#!/bin/bash

# --- Setup ---
# Smoke test for stb-chargediffAnalysis (Charge Density Difference, Stage 2 - Analysis, item 4.19.2)
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
# Tiny synthetic .RHO/.fdf fixtures are generated on the fly (via sisl) rather
# than committed to the repo -- see make_synthetic_rho.py's own docstring for
# why (this test only needs to exercise the subtraction/shape-check/cube
# -write plumbing, not real physics).
echo "--- Starting tester for stb-chargediffAnalysis (item 4.19.2) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
python3 "$FIXTURE_DIR/make_synthetic_rho.py"
mv "$FIXTURE_DIR"/*.RHO "$FIXTURE_DIR"/combined.fdf "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Basic 3-file subtraction (combined - frag1 - frag2), default slice mode ---
echo -e "\n--- Testing a basic 2-fragment Delta rho (default slice mode) ---"
stb-chargediffAnalysis --combined-label combined --fragment-labels frag1 frag2 \
    --no-intro --output-dir out_basic --save-report > log_basic.txt 2>&1
check_exit_code $? 0
check_success out_basic/combined_chargediff.dat
check_contains "Integrated Delta Charge Density" out_basic/chargediff_analysis_report.txt


# --- 3. --profile mode ---
echo -e "\n--- Testing --profile mode ---"
stb-chargediffAnalysis --combined-label combined --fragment-labels frag1 frag2 \
    --profile --axis 2 --no-intro --output-dir out_profile > log_profile.txt 2>&1
check_exit_code $? 0
check_success out_profile/combined_chargediff.dat


# --- 4. --cube ---
echo -e "\n--- Testing --cube (needs combined.fdf for geometry) ---"
stb-chargediffAnalysis --combined-label combined --fragment-labels frag1 frag2 \
    --cube --no-intro --output-dir out_cube > log_cube.txt 2>&1
check_exit_code $? 0
check_success out_cube/combined_chargediff.cube


# --- 5. Physics: refuse a shape-mismatched fragment ---
echo -e "\n--- Testing a shape-mismatched fragment is rejected ---"
stb-chargediffAnalysis --combined-label combined --fragment-labels frag1 badshape \
    --no-intro --output-dir out_bad_shape > log_bad_shape.txt 2>&1
check_exit_code $? 1
check_contains "Grid shape mismatch" log_bad_shape.txt


# --- 6. Physics: refuse fewer than 2 fragment labels ---
echo -e "\n--- Testing --fragment-labels with only 1 label is rejected ---"
stb-chargediffAnalysis --combined-label combined --fragment-labels frag1 \
    --no-intro --output-dir out_one_frag > log_one_frag.txt 2>&1
check_exit_code $? 1
check_contains "at least 2 labels" log_one_frag.txt


# --- 7. Error and robustness cases ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: missing --combined-label/--fragment-labels"
stb-chargediffAnalysis --no-intro > log_missing_args.txt 2>&1
check_exit_code $? 2

echo "Testing: nonexistent combined label"
stb-chargediffAnalysis --combined-label nope --fragment-labels frag1 frag2 \
    --no-intro --output-dir out_missing > log_missing_input.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_input.txt

echo "Testing: --version"
stb-chargediffAnalysis --version > log_version.txt 2>&1
check_contains "stb-chargediffAnalysis" log_version.txt

echo "Testing: --help documents --combined-label, --fragment-labels, --profile, --cube"
stb-chargediffAnalysis --help > log_help.txt 2>&1
check_contains "combined-label" log_help.txt
check_contains "fragment-labels" log_help.txt
check_contains "profile" log_help.txt
check_contains "cube" log_help.txt


# --- 8. Interactive path (stb-suite, shortcut 4.19.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.19.2) ---"
rm -rf out_interactive
printf '4.19.2\ncombined\nfrag1 frag2\n1\nn\nn\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success combined_chargediff.dat


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
