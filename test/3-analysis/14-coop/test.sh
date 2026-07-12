#!/bin/bash

# --- Setup ---
# Smoke test for stb-coop (COOP/COHP Bonding Analysis, item 3.14)
#
# Fixture note: Sn3O4.selected.WFSX + Sn3O4.HSX are a real SIESTA run
# (restarted from the converged Sn3O4.DM in test/5-utils/4-wantibexos/, but
# with a small explicit 2x2x2-equivalent k-mesh via WriteWaveFunctions T +
# %block WaveFuncKPoints instead of a band path -- COOP/COHP needs a real
# BZ sampling, not a band-path WFSX like stb-fatbands uses). Small (~6 MB)
# since this is only 8 k-points. calc.fdf's real fdf filename does not
# match SystemLabel (Sn3O4) -- same pattern as test/3-analysis/10-fatbands/,
# so tests use --wfsx/--hsx-file explicitly where the --label path can't
# find <label>.fdf (not needed here since coop.py never reads .fdf directly).
export MPLBACKEND=Agg
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
echo "--- Starting tester for stb-coop (item 3.14) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR"/Sn3O4.selected.WFSX "$FIXTURE_DIR"/Sn3O4.HSX "$FIXTURE_DIR"/Sn3O4.XV "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null

# Keep runtimes reasonable: --npoints modest, since sisl's COOP/COHP is
# memory/CPU-heavy per its own docstring warning.
COMMON="--erange -10 5 --npoints 25 --sigma 300"


# --- 2. Basic run: --quantity coop, --pair (label auto-detect via --label) ---
echo -e "\n--- Testing a basic COOP run (--label Sn3O4 --pair 0 1) ---"

timeout 120 stb-coop --label Sn3O4 --quantity coop --pair 0 1 $COMMON --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "Using 'Sn3O4.HSX'" log_basic.txt
check_contains "Selected pairs: \['0-1'\]" log_basic.txt
check_success coop.dat


# --- 3. --quantity cohp + --pair-species ---
echo -e "\n--- Testing --quantity cohp with --pair-species ---"

timeout 120 stb-coop --label Sn3O4 --quantity cohp --pair-species Sn O $COMMON -o cohp_out --no-intro > log_cohp.txt 2>&1
check_exit_code $? 0
check_contains "Selected pairs: \['Sn-O'\]" log_cohp.txt
check_success cohp_out/cohp.dat


# --- 4. Multiple --pair given together ---
echo -e "\n--- Testing multiple --pair (repeatable) ---"

timeout 120 stb-coop --label Sn3O4 --quantity coop --pair 0 1 --pair 6 7 $COMMON -o multi_out --no-intro > log_multi.txt 2>&1
check_exit_code $? 0
check_success multi_out/coop.dat
echo "Testing: both pair columns present in the header"
check_contains "0-1" multi_out/coop.dat
check_contains "6-7" multi_out/coop.dat


# --- 5. --bond-order complementary sanity check ---
echo -e "\n--- Testing --bond-order ---"

timeout 120 stb-coop --label Sn3O4 --quantity coop --pair 0 1 $COMMON --bond-order -o bo_out --no-intro > log_bo.txt 2>&1
check_exit_code $? 0
check_contains "Complementary bond_order" log_bo.txt


# --- 6. Explicit --wfsx/--hsx-file (no --label) ---
echo -e "\n--- Testing explicit --wfsx/--hsx-file (no --label) ---"

timeout 120 stb-coop --wfsx Sn3O4.selected.WFSX --hsx-file Sn3O4.HSX --quantity coop --pair 0 1 $COMMON -o explicit_out --no-intro > log_explicit.txt 2>&1
check_exit_code $? 0
check_success explicit_out/coop.dat


# --- 7. --shift fermi ---
echo -e "\n--- Testing --shift fermi ---"

timeout 120 stb-coop --label Sn3O4 --quantity coop --pair 0 1 --erange -3 3 --npoints 20 --sigma 300 \
    --shift fermi --fermi -3.200055 -o fermi_out --no-intro > log_fermi.txt 2>&1
check_exit_code $? 0
check_success fermi_out/coop.dat


# --- 8. Robustness / errors ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: --wfsx without --hsx-file"
stb-coop --wfsx Sn3O4.selected.WFSX --quantity coop --pair 0 1 $COMMON --no-intro > log_missing_hsx.txt 2>&1
check_exit_code $? 2

echo "Testing: missing --pair/--pair-species"
stb-coop --label Sn3O4 --quantity coop $COMMON --no-intro > log_missing_pair.txt 2>&1
check_exit_code $? 2

echo "Testing: --pair index out of range"
stb-coop --label Sn3O4 --quantity coop --pair 0 999 $COMMON --no-intro > log_bad_pair.txt 2>&1
check_exit_code $? 2
check_contains "out of range" log_bad_pair.txt

echo "Testing: missing --quantity"
stb-coop --label Sn3O4 --pair 0 1 $COMMON --no-intro > log_missing_quantity.txt 2>&1
check_exit_code $? 2

echo "Testing: missing --label/--wfsx"
stb-coop --quantity coop --pair 0 1 $COMMON --no-intro > log_missing_label.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-coop --version > log_version.txt 2>&1
check_contains "stb-coop" log_version.txt

echo "Testing: --help documents --quantity, --pair, --pair-species, --erange, --sigma, --bond-order"
stb-coop --help > log_help.txt 2>&1
check_contains "quantity" log_help.txt
check_contains "pair-species" log_help.txt
check_contains "erange" log_help.txt
check_contains "sigma" log_help.txt
check_contains "bond-order" log_help.txt


# --- 9. Interactive path (stb-suite, shortcut 3.14) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 3.14) ---"

echo "Testing: navigate 3.14 -> label Sn3O4 -> quantity 1 (coop) -> pair mode 1 (index) -> 0 1 -> erange -10 5 -> sigma 300 -> default output"
rm -f coop.dat
printf '3.14\nSn3O4\n1\n1\n0\n1\n-10\n5\n300\n\n\n0\n' | timeout 120 stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success coop.dat


popd > /dev/null

# --- 10. Summary ---
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
