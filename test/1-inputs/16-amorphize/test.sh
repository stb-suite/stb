#!/bin/bash

# --- Setup ---
# Smoke test for stb-amorphize (Amorphous Structure Generator, item 1.16)
#
# Needs the optional 'ml' extra (pip install stb_suite[ml] -- PyTorch +
# mace-torch). The whole file is skipped with a clear message if `mace`
# isn't importable, same pattern as 8-defect/15-mlrelax's tests.
#
# Uses a tiny 8-atom cell with deliberately tiny --melt-steps/
# --quench-steps so this smoke test finishes quickly -- it only checks
# that the tool runs correctly end-to-end and produces valid output, not
# that the result is a high-quality amorphous structure (that physics was
# validated separately during planning with a realistic 64-atom cell and
# step counts; see CLAUDE.md's stb-amorphize entry).
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


echo "--- Starting tester for STB-Amorphize (item 1.16) ---"

if ! python3 -c "import mace" 2>/dev/null; then
    echo -e "${YELLOW}Skipped entirely:${NC} the optional 'ml' extra is not installed."
    echo "Install with: pip install stb_suite[ml]  (then re-run this test)"
    exit 0
fi


# --- 1. Preparation ---
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/si8.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Happy path: quick melt-quench on si8.fdf ---
echo -e "\n--- Testing default protocol (tiny step counts for speed) ---"
rm -f amorphous.fdf
stb-amorphize -f si8.fdf --melt-steps 20 --quench-steps 30 --no-intro \
    > log_default.txt 2>&1
check_exit_code $? 0
check_contains "Input bond-angle mean/std" log_default.txt
check_contains "Post-quench bond-angle mean/std" log_default.txt
check_contains "Steps used" log_default.txt
check_success amorphous.fdf


# --- 3. --no-final-relax skips the cleanup relax ---
echo -e "\n--- Testing --no-final-relax ---"
rm -f noreelax.fdf
stb-amorphize -f si8.fdf --melt-steps 10 --quench-steps 10 --no-final-relax \
    -o noreelax.fdf --no-intro > log_norelax.txt 2>&1
check_exit_code $? 0
if grep -q "Final static relax" log_norelax.txt; then
    echo -e "   -> ${RED}Failed:${NC} 'Final static relax' should NOT appear with --no-final-relax"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} 'Final static relax' correctly absent"
    PASS=$((PASS+1))
fi
check_success noreelax.fdf


# --- 4. Rejected on a vacuum-padded slab (bulk-only tool) ---
echo -e "\n--- Testing rejection on a vacuum-padded slab ---"
stb-slab -f si8.fdf --hkl 1 0 0 -o slab_for_test.fdf --no-intro > /dev/null 2>&1
stb-amorphize -f slab_for_test.fdf --no-intro > log_vacuum.txt 2>&1
check_exit_code $? 1
check_contains "only supports bulk" log_vacuum.txt
rm -f slab_for_test.fdf


# --- 5. --help / --version ---
echo -e "\n--- Testing --help / --version ---"
stb-amorphize --version > log_version.txt 2>&1
check_contains "stb-amorphize" log_version.txt

stb-amorphize --help > log_help.txt 2>&1
check_contains "melt-temp" log_help.txt
check_contains "quench-steps" log_help.txt
check_contains "no-final-relax" log_help.txt


# --- 6. Interactive path (stb-suite, shortcut 1.16) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 1.16) ---"

echo "Testing: navigate 1.16 -> invalid file then valid -> tiny steps -> final relax yes -> default output -> quit"
rm -f amorphous.fdf
printf '1.16\ndoes_not_exist.fdf\nsi8.fdf\n3000\n15\n300\n20\n\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_contains "File not found" log_menu.txt
check_contains "Post-quench bond-angle mean/std" log_menu.txt
check_success amorphous.fdf


popd > /dev/null

# --- 7. Summary ---
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
