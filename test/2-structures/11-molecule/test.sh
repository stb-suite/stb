#!/bin/bash

# --- Setup ---
# Smoke test for stb-molecule (Reference Molecule Builder, item 2.11)
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
echo "--- Starting tester for STB-Molecule (item 2.11) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Happy path: H2O ---
echo -e "\n--- Testing --name H2O ---"
rm -f h2o.fdf
stb-molecule --name H2O -o h2o.fdf --no-intro > log_h2o.txt 2>&1
check_contains "Formula:.*H2O" log_h2o.txt
check_contains "Atoms:.*3" log_h2o.txt
check_success h2o.fdf
check_contains "NumberofAtoms      3" h2o.fdf
check_contains "NumberOfSpecies    2" h2o.fdf


# --- 3. Multi-element molecule: CH3OH ---
echo -e "\n--- Testing --name CH3OH (3 species: C, O, H) ---"
rm -f ch3oh.fdf
stb-molecule --name CH3OH -o ch3oh.fdf --no-intro > log_ch3oh.txt 2>&1
check_contains "Atoms:.*6" log_ch3oh.txt
check_success ch3oh.fdf
check_contains "NumberOfSpecies    3" ch3oh.fdf
check_contains "NumberofAtoms      6" ch3oh.fdf


# --- 4. --vacuum changes the output cell size ---
echo -e "\n--- Testing --vacuum override ---"
rm -f h2o_small.fdf h2o_big.fdf
stb-molecule --name H2O --vacuum 5 -o h2o_small.fdf --no-intro > /dev/null 2>&1
stb-molecule --name H2O --vacuum 20 -o h2o_big.fdf --no-intro > /dev/null 2>&1
small_a=$(sed -n '/LatticeVectors/{n;p;q}' h2o_small.fdf | awk '{print $1}')
big_a=$(sed -n '/LatticeVectors/{n;p;q}' h2o_big.fdf | awk '{print $1}')
if awk -v a="$big_a" -v b="$small_a" 'BEGIN { exit !(a > b) }'; then
    echo -e "   -> ${GREEN}Verified:${NC} --vacuum 20 cell ($big_a) > --vacuum 5 cell ($small_a)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} --vacuum did not scale the cell ($small_a vs $big_a)"
    FAIL=$((FAIL+1))
fi


# --- 5. --list ---
echo -e "\n--- Testing --list ---"
stb-molecule --list --no-intro > log_list.txt 2>&1
check_exit_code $? 0
check_contains "Available molecules (162)" log_list.txt
check_contains "H2O" log_list.txt
check_contains "C6H6" log_list.txt


# --- 6. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: unknown molecule name"
stb-molecule --name NotAMolecule --no-intro > log_bad_name.txt 2>&1
check_exit_code $? 1
check_contains "is not a known molecule name" log_bad_name.txt
check_contains "\-\-list" log_bad_name.txt

echo "Testing: wrong-case name gives a 'did you mean' pointing at the exact name"
stb-molecule --name h2o --no-intro > log_wrong_case.txt 2>&1
check_exit_code $? 1
check_contains "Did you mean 'H2O'" log_wrong_case.txt

echo "Testing: --name omitted (and no --list)"
stb-molecule --no-intro > log_missing_name.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-molecule --version > log_version.txt 2>&1
check_contains "stb-molecule" log_version.txt

echo "Testing: --help documents --name/--list/--vacuum"
stb-molecule --help > log_help.txt 2>&1
check_contains "\-\-name" log_help.txt
check_contains "\-\-list" log_help.txt
check_contains "vacuum" log_help.txt


# --- 7. Interactive path (stb-suite, shortcut 1.14) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 1.14) ---"

echo "Testing: navigate 1.14 -> H2O -> default vacuum -> default output -> quit"
rm -f molecule.fdf
printf '2.11\nH2O\n\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_contains "Formula:.*H2O" log_menu.txt
check_success molecule.fdf


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
