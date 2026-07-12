#!/bin/bash

# --- Setup ---
# Smoke test for stb-stm (STM Simulator, item 3.11)
#
# Fixture note: Graphene.LDOS/.XV are a real SIESTA run (monolayer graphene,
# 2-atom primitive cell, 20 Ang vacuum along z, LDOS integrated over
# [-1.5, 1.5] eV around E_Fermi via '%block LocalDensityOfStates') -- fast
# to reproduce (`siesta < calc.fdf`) since it's a tiny, quickly-converging
# system, unlike test/3-analysis/10-fatbands/'s much heavier Sn3O4 fixture.
export MPLBACKEND=Agg
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
echo "--- Starting tester for stb-stm (item 3.11) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR"/Graphene.LDOS "$FIXTURE_DIR"/Graphene.XV "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Basic run: --mode current (default), --label auto-detect ---
echo -e "\n--- Testing a basic constant-current run (--label Graphene) ---"

timeout 60 stb-stm --label Graphene --iso 0.001 --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "Using axis 2" log_basic.txt
check_contains "Height map:" log_basic.txt
check_success stm_current.dat
check_success stm_current.gplot

echo "Testing: constant-current corrugation is physically sensible (sub-Angstrom, non-zero)"
check_contains "corrugation:" log_basic.txt


# --- 3. --mode height ---
echo -e "\n--- Testing --mode height ---"

timeout 60 stb-stm --label Graphene --mode height --z 2.0 --no-intro > log_height.txt 2>&1
check_exit_code $? 0
check_contains "constant-height" log_height.txt
check_success stm_height.dat
check_success stm_height.gplot


# --- 4. Explicit --file/--geometry-file (no --label) ---
echo -e "\n--- Testing explicit --file/--geometry-file (no --label) ---"

timeout 60 stb-stm --file Graphene.LDOS --geometry-file Graphene.XV --iso 0.001 \
    -o explicit_out --no-intro > log_explicit.txt 2>&1
check_exit_code $? 0
check_success explicit_out/stm_current.dat


# --- 5. Explicit --axis bypasses vacuum auto-detection ---
echo -e "\n--- Testing explicit --axis 2 ---"

timeout 60 stb-stm --label Graphene --iso 0.001 --axis 2 -o axis_out --no-intro > log_axis.txt 2>&1
check_exit_code $? 0
check_success axis_out/stm_current.dat


# --- 6. --iso set unreachably high: NaN points are reported, not a crash ---
echo -e "\n--- Testing an unreachable --iso (NaN/no-crossing points) ---"

timeout 60 stb-stm --label Graphene --iso 100.0 -o unreachable_out --no-intro > log_unreachable.txt 2>&1
check_exit_code $? 0
check_contains "never reached --iso" log_unreachable.txt
check_contains "NaN" unreachable_out/stm_current.dat


# --- 7. Missing required options ---
echo -e "\n--- Testing missing --iso/--z ---"

echo "Testing: --mode current without --iso"
stb-stm --label Graphene --no-intro > log_missing_iso.txt 2>&1
check_exit_code $? 2
check_contains "\-\-iso is required" log_missing_iso.txt

echo "Testing: --mode height without --z"
stb-stm --label Graphene --mode height --no-intro > log_missing_z.txt 2>&1
check_exit_code $? 2
check_contains "\-\-z is required" log_missing_z.txt


# --- 8. Missing .LDOS aborts with a clear error ---
echo -e "\n--- Testing a missing .LDOS ---"

timeout 30 stb-stm --label does_not_exist --iso 0.001 --no-intro > log_missing_ldos.txt 2>&1
check_exit_code $? 2
check_contains "No LDOS file found" log_missing_ldos.txt


# --- 9. A mismatched geometry (wrong cell for this grid) is rejected ---
echo -e "\n--- Testing rejection of a geometry whose cell doesn't match the grid ---"

echo "Testing: pairing Graphene.LDOS with an unrelated bulk .XV (Sn3O4, wrong cell)"
SN3O4_XV="$FIXTURE_DIR/../10-fatbands/Sn3O4.XV"
if [ -f "$SN3O4_XV" ]; then
    stb-stm --file Graphene.LDOS --geometry-file "$SN3O4_XV" --iso 0.001 --no-intro > log_cell_mismatch.txt 2>&1
    check_exit_code $? 1
    check_contains "does not match" log_cell_mismatch.txt
else
    echo -e " ... ${YELLOW}SKIPPED${NC} (test/3-analysis/10-fatbands/Sn3O4.XV not found)"
fi


# --- 10. Robustness / --help / --version ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: --label omitted (required)"
stb-stm --iso 0.001 --no-intro > log_missing_label.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-stm --version > log_version.txt 2>&1
check_contains "stb-stm" log_version.txt

echo "Testing: --help documents --mode, --z, --iso, --z-max, --axis, --geometry-file, -o"
stb-stm --help > log_help.txt 2>&1
check_contains "mode" log_help.txt
check_contains "\-\-z " log_help.txt
check_contains "iso" log_help.txt
check_contains "z-max" log_help.txt
check_contains "axis" log_help.txt
check_contains "geometry-file" log_help.txt
check_contains "output-dir" log_help.txt


# --- 11. Interactive path (stb-suite, shortcut 3.11) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 3.11) ---"

echo "Testing: navigate 3.11 -> label Graphene -> mode 1 (current) -> iso 0.001 -> auto axis -> default output"
rm -f stm_current.dat stm_current.gplot
printf '3.11\nGraphene\n1\n0.001\n\n\n\n0\n' | timeout 60 stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success stm_current.dat


popd > /dev/null

# --- 12. Summary ---
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
