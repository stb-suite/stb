#!/bin/bash

# --- Setup ---
# Smoke test for stb-cohesive (Cohesive Energy Prep, item 4.3.1)
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
    if grep -q -- "$1" "$2" 2>/dev/null; then
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
echo "--- Starting tester for STB-Cohesive prep (item 4.3.1) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR/pp"
cp "$FIXTURE_DIR/structure.fdf" "$TEST_DIR/"
# Minimal placeholder pseudopotentials -- only existence/format matters here,
# not real pseudopotential content.
echo "dummy pseudopotential content" > "$TEST_DIR/pp/C.psf"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Default run: vacuum-aware k-grid (10-atom graphene-like slab, 20 Ang
#     out-of-plane vacuum -> c-axis forced to a single division) ---
echo -e "\n--- Testing default run (vacuum-aware k-grid) ---"
rm -rf structure atoms
stb-cohesive -s structure.fdf --no-intro > log_default.txt 2>&1
check_exit_code $? 0
check_success structure/calc.fdf
check_success structure/structure.fdf
check_success atoms/C/calc.fdf
check_success atoms/C/structure.fdf

echo "Verifying the in-plane axes get a real k-grid but the vacuum axis is forced to 1"
check_contains "Calculated K-grid for full structure: 7 7 1" log_default.txt
check_contains "kgrid.MonkhorstPack   \[7  7  1\]" structure/calc.fdf

echo "Verifying the isolated atom is always Gamma-only"
check_contains "kgrid.MonkhorstPack   \[1  1  1\]" atoms/C/calc.fdf

echo "Verifying spin: full structure defaults to non-polarized, isolated atom is always polarized"
check_contains "Spin                non-polarized" structure/calc.fdf
check_contains "Spin                polarized" atoms/C/calc.fdf


# --- 3. --spin: only the full structure calc changes ---
echo -e "\n--- Testing --spin ---"
rm -rf structure atoms
stb-cohesive -s structure.fdf --spin --no-intro > log_spin.txt 2>&1
check_exit_code $? 0
check_contains "Spin                polarized" structure/calc.fdf
check_contains "Spin polarization ENABLED" log_spin.txt


# --- 4. --pp-path: .psf pseudopotentials are linked (regression test -- this
#     tool used to only recognize .psml and silently drop .psf, leaving the
#     generated directories without any pseudopotential at all) ---
echo -e "\n--- Testing --pp-path with a .psf pseudopotential ---"
rm -rf structure atoms
stb-cohesive -s structure.fdf -p pp --no-intro > log_psf.txt 2>&1
check_exit_code $? 0
check_success structure/C.psf
check_success atoms/C/C.psf
if grep -qi "not found" log_psf.txt; then
    echo -e "   -> ${RED}Failed:${NC} a pseudopotential 'not found' warning was printed despite pp/C.psf existing"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no 'not found' warning printed"
    PASS=$((PASS+1))
fi


# --- 5. --pp-path: .psml takes priority over .psf when both exist ---
echo -e "\n--- Testing .psml priority over .psf ---"
rm -rf structure atoms
echo "dummy pseudopotential content" > pp/C.psml
stb-cohesive -s structure.fdf -p pp --no-intro > log_psml.txt 2>&1
check_exit_code $? 0
check_success structure/C.psml
if [ -e structure/C.psf ]; then
    echo -e "   -> ${RED}Failed:${NC} C.psf was also linked even though C.psml exists"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} only C.psml was linked, not C.psf"
    PASS=$((PASS+1))
fi


# --- 6. --pp-path accepts a bundled bank name ---
echo -e "\n--- Testing --pp-path with a bundled bank name (dojo) ---"
rm -rf structure atoms
stb-cohesive -s structure.fdf -p dojo --no-intro > log_bank.txt 2>&1
check_exit_code $? 0
check_success structure/C.psml
rm -f pp/C.psml

echo "Testing: using a bundled bank prints its citation/origin"
check_contains "van Setten" log_bank.txt
check_contains "pseudo-dojo.org" log_bank.txt


# --- 6. --vacuum: isolated-atom box size is configurable ---
echo -e "\n--- Testing --vacuum ---"
rm -rf structure atoms
stb-cohesive -s structure.fdf --vacuum 30 --no-intro > log_vacuum.txt 2>&1
check_exit_code $? 0
check_contains "30.000000   0.000000   0.000000" atoms/C/structure.fdf

echo "Testing: --vacuum 0 is rejected"
stb-cohesive -s structure.fdf --vacuum 0 --no-intro > log_vacuum_zero.txt 2>&1
check_exit_code $? 2


# --- 7. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: missing structure file"
stb-cohesive -s does_not_exist.fdf --no-intro > log_missing.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing.txt

echo "Testing: missing required -s"
stb-cohesive --no-intro > log_missing_arg.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-cohesive --version > log_version.txt 2>&1
check_contains "stb-cohesive" log_version.txt


# --- 8. Interactive path (stb-suite, shortcut 3.3.1) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 3.3.1) ---"

echo "Testing: navigate 3.3.1 -> structure.fdf -> defaults -> no pp -> no spin -> default vacuum"
rm -rf structure atoms
printf '4.3.1\nstructure.fdf\n\n\nn\n\n' | stb-suite > log_menu.txt 2>&1
check_success structure/calc.fdf
check_success atoms/C/calc.fdf


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
