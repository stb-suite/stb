#!/bin/bash

# --- Setup ---
# Smoke test for stb-crystalbuilder (Crystal Builder, item 1.11)
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
echo "--- Starting tester for STB-Crystalbuilder (item 1.11) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. FCC Ni (single site, no free parameter) ---
echo -e "\n--- Testing FCC Ni (--spacegroup Fm-3m) ---"

echo "Testing: symbol form of --spacegroup"
rm -f fcc_ni.fdf
stb-crystalbuilder --spacegroup Fm-3m --a 3.52 --site Ni 0 0 0 -o fcc_ni.fdf --no-intro > log_fcc_symbol.txt 2>&1
check_contains "Output atoms:.*4" log_fcc_symbol.txt
check_contains "Detected space group:.*Fm-3m (No. 225)" log_fcc_symbol.txt
check_success fcc_ni.fdf
check_contains "NumberofAtoms      4" fcc_ni.fdf

echo "Testing: int form of --spacegroup (should behave the same)"
rm -f fcc_ni2.fdf
stb-crystalbuilder --spacegroup 225 --a 3.52 --site Ni 0 0 0 -o fcc_ni2.fdf --no-intro > log_fcc_int.txt 2>&1
check_contains "Output atoms:.*4" log_fcc_int.txt
check_success fcc_ni2.fdf
if diff -q fcc_ni.fdf fcc_ni2.fdf > /dev/null; then
    echo -e "   -> ${GREEN}Verified:${NC} symbol and int spacegroup forms produce identical output"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} symbol and int spacegroup forms differ"
    FAIL=$((FAIL+1))
fi


# --- 3. Magnetite (multi-site, real oxide, cross-checked with stb-unitcell) ---
echo -e "\n--- Testing magnetite (--spacegroup Fd-3m, 2 Fe sites + 1 O site) ---"

rm -f magnetite.fdf
stb-crystalbuilder --spacegroup Fd-3m --a 8.396 \
    --site Fe 0.125 0.125 0.125 --site Fe 0.5 0.5 0.5 --site O 0.379 0.379 0.379 \
    -o magnetite.fdf --no-intro > log_magnetite.txt 2>&1
check_contains "Output formula:.*Fe3O4" log_magnetite.txt
check_contains "Output atoms:.*56" log_magnetite.txt
check_contains "Detected space group:.*Fd-3m (No. 227)" log_magnetite.txt
if grep -q "unusually close" log_magnetite.txt; then
    echo -e "   -> ${RED}Failed:${NC} unexpected overlap warning for physically valid magnetite"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no spurious overlap warning"
    PASS=$((PASS+1))
fi
check_success magnetite.fdf

echo "Testing: cross-check with stb-unitcell (56-atom conventional -> 14-atom primitive)"
stb-unitcell -f magnetite.fdf --mode primitive --no-intro -o magnetite_prim.fdf > log_crosscheck.txt 2>&1
check_contains "Output atoms:.*14" log_crosscheck.txt
check_contains "Space group:.*Fd-3m" log_crosscheck.txt
check_success magnetite_prim.fdf
rm -f unitcell_report.txt


# --- 4. Overlap sanity check ---
echo -e "\n--- Testing overlap sanity check (two --site at the same coordinate) ---"
stb-crystalbuilder --spacegroup Fm-3m --a 3.52 --site Ni 0 0 0 --site Fe 0 0 0 \
    -o overlap.fdf --no-intro > log_overlap.txt 2>&1
check_exit_code $? 0
check_contains "unusually close" log_overlap.txt
check_success overlap.fdf


# --- 5. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: lattice incompatible with requested space group"
stb-crystalbuilder --spacegroup Fm-3m --a 3 --b 4 --site Ni 0 0 0 --no-intro > log_bad_lattice.txt 2>&1
check_exit_code $? 1
check_contains "incompatible" log_bad_lattice.txt

echo "Testing: unrecognized space group symbol"
stb-crystalbuilder --spacegroup NotASpaceGroup --a 3.52 --site Ni 0 0 0 --no-intro > log_bad_sg.txt 2>&1
check_exit_code $? 1

echo "Testing: invalid element symbol in --site"
stb-crystalbuilder --spacegroup Fm-3m --a 3.52 --site Xx 0 0 0 --no-intro > log_bad_element.txt 2>&1
check_exit_code $? 1
check_contains "not a valid Element" log_bad_element.txt

echo "Testing: non-numeric --site coordinate"
stb-crystalbuilder --spacegroup Fm-3m --a 3.52 --site Ni a b c --no-intro > log_bad_coords.txt 2>&1
check_exit_code $? 1
check_contains "must be numbers" log_bad_coords.txt

echo "Testing: --site omitted entirely (required)"
stb-crystalbuilder --spacegroup Fm-3m --a 3.52 --no-intro > log_missing_site.txt 2>&1
check_exit_code $? 2

echo "Testing: --spacegroup omitted (required)"
stb-crystalbuilder --a 3.52 --site Ni 0 0 0 --no-intro > log_missing_sg.txt 2>&1
check_exit_code $? 2

echo "Testing: --a omitted (required)"
stb-crystalbuilder --spacegroup Fm-3m --site Ni 0 0 0 --no-intro > log_missing_a.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-crystalbuilder --version > log_version.txt 2>&1
check_contains "stb-crystalbuilder" log_version.txt

echo "Testing: --help documents --site and --spacegroup"
stb-crystalbuilder --help > log_help.txt 2>&1
check_contains "site" log_help.txt
check_contains "spacegroup" log_help.txt


# --- 6. Interactive path (stb-suite, shortcut 1.11) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 1.11) ---"

echo "Testing: navigate 1.11 -> Fm-3m -> a=3.52, defaults for b/c/angles -> site 'Ni 0 0 0' -> blank to finish -> default output -> quit"
rm -f crystal.fdf
printf '1.11\nFm-3m\n3.52\n\n\n\n\n\nNi 0 0 0\n\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_contains "Output atoms:.*4" log_menu.txt
check_success crystal.fdf


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
