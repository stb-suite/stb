#!/bin/bash

# --- Setup ---
# Smoke test for stb-nativecharges (Native Charges reader, item 2.19).
#
# Fixtures:
#   both/    -- a real spin-polarized SIESTA .out with Charge.Hirshfeld end
#               and Charge.Voronoi end both set (Te/NO2 adsorption system,
#               84 atoms).
#               Atom 73 is the Te atom the proposal itself calls out: Hirshfeld
#               charge +0.393899 e, Voronoi +0.400221 e.
#   neither/ -- the same singlepoint fixture used by 6-status/6-archive
#               (Sn3O4), which has neither native charges section, to
#               exercise the graceful [WARNING] path.
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$FIXTURE_DIR/test_files"

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

PASS=0
FAIL=0

check_contains() {
    if grep -q "$1" "$2" 2>/dev/null; then
        echo -e "   -> ${GREEN}Verified:${NC} '$1' found in '$2'"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} '$1' NOT found in '$2'"
        FAIL=$((FAIL+1))
    fi
}

check_success() {
    if [ -s "$1" ]; then
        echo -e " ... ${GREEN}OK${NC} (file '$1' created)"
        PASS=$((PASS+1))
    else
        echo -e " ... ${RED}FAIL${NC} (file '$1' was not created)"
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
echo "--- Starting tester for STB-NATIVECHARGES (item 2.19) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR/both" "$TEST_DIR/neither"
cp "$FIXTURE_DIR/both/calc.out" "$FIXTURE_DIR/both/calc.fdf" "$TEST_DIR/both/"
cp "$FIXTURE_DIR/neither/calc.out" "$FIXTURE_DIR/neither/calc.fdf" "$TEST_DIR/neither/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Both sections present (auto-detects SystemLabel=siesta) ---
echo -e "\n--- Testing folder with both Hirshfeld and Voronoi sections ---"
rm -f stb_nativecharges_report*.txt stb_nativecharges_references*.bib
stb-nativecharges --no-intro --path both > log_both.txt 2>&1
check_exit_code $? 0
check_contains "Log               : both/calc.out (SIESTA 5.4.2)" log_both.txt
check_contains "Hirshfeld         : found" log_both.txt
check_contains "Voronoi           : found" log_both.txt
check_contains "Mulliken          : not in the log (add WriteMullikenPop 1)" log_both.txt
echo "Testing: [2] CHARGES PER ATOM puts every method side by side (Te 73: Hirshfeld +0.3939, Voronoi +0.4002)"
check_contains "\[2\] CHARGES PER ATOM" log_both.txt
check_contains "73 | Te | 6     | +0.3939     | +0.4002   | 0.0063" log_both.txt
check_contains "\[3\] PER SPECIES" log_both.txt
check_contains "Te      | 9     | +0.2918 +-0.0361 | +0.2907 +-0.0387" log_both.txt
echo "Testing: [4] CHECKS (sum of q vs -NetCharge, population + q = Z_val, sum of Sz)"
check_contains "\[4\] CHECKS" log_both.txt
check_contains "Hirshfeld | -0.000084    | 4.4e-16" log_both.txt
check_contains "Method dependence : up to 0.027 e (N82)" log_both.txt
check_contains "Hirshfeld         : in 1 of 1 folder(s)" log_both.txt
echo "Testing: the report is always written, in the current folder, with its references file"
check_success stb_nativecharges_report.txt
check_contains "STB-NATIVECHARGES REPORT" stb_nativecharges_report.txt
check_success stb_nativecharges_references.bib


# --- 3. Neither section present (graceful degradation) ---
echo -e "\n--- Testing folder with neither section (graceful [WARNING]) ---"
stb-nativecharges --no-intro --path neither > log_neither.txt 2>&1
check_exit_code $? 0
check_contains "Hirshfeld         : not in the log (add Charge.Hirshfeld end)" log_neither.txt
check_contains "no Mulliken, Hirshfeld or Voronoi block in the log" log_neither.txt
check_contains "Hirshfeld         : in 0 of 1 folder(s)" log_neither.txt
echo "Testing: a second run never overwrites the first report (numbered _1)"
check_success stb_nativecharges_report_1.txt
check_contains "Target            : neither" stb_nativecharges_report_1.txt


# --- 4. --label override ---
echo -e "\n--- Testing --label override ---"
stb-nativecharges --no-intro --path both --label siesta > log_labelover.txt 2>&1
check_exit_code $? 0
check_contains "73 | Te | 6     | +0.3939" log_labelover.txt


# --- 5. Batch mode (glob matching both folders) ---
echo -e "\n--- Testing batch mode (--path matching both folders) ---"
stb-nativecharges --no-intro --path "*" > log_batch.txt 2>&1
check_contains "Target            : \* (2 folders)" log_batch.txt
check_contains "\[1\] CHARGES PER SPECIES" log_batch.txt
check_contains "both    | DZP   | Te      | 9     | -               | +0.2918          | +0.2907" log_batch.txt
check_contains "Hirshfeld         : in 1 of 2 folder(s)" log_batch.txt


# --- 6. --save-report is a hidden, ignored option (old scripts keep working) ---
echo -e "\n--- Testing the old --save-report flag ---"
stb-nativecharges --no-intro --path both --save-report > log_savereport.txt 2>&1
check_exit_code $? 0
check_contains "Native charges report complete" log_savereport.txt


# --- 7. Missing directory ---
echo -e "\n--- Testing missing directory ---"
stb-nativecharges --no-intro --path does_not_exist > log_missing.txt 2>&1
check_exit_code $? 2
check_contains "is not a directory" log_missing.txt


# --- 8. Interactive path (stb-suite, shortcut 2.19) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 2.19) ---"
printf '2.19\nboth\n\n\n\n0\n' | stb-suite > log_interactive.txt 2>&1
check_contains "\[2\] CHARGES PER ATOM" log_interactive.txt
check_contains "Native charges report complete" log_interactive.txt


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
