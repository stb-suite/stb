#!/bin/bash

# --- Setup ---
# Smoke test for stb-status (Status Checker, item 5.6).
#
# Fixtures are 3 real, distinctly-typed SIESTA runs:
#   singlepoint/  -- Sn3O4.out (calc.out/calc.fdf), single-point, shared with
#                    the other 5-utils fixtures.
#   aimd/         -- the same 5-step Verlet AIMD run used by stb-ani2traj's
#                    own fixtures (5 steps, final Temp_ion = 9917.825 K).
#   relax/        -- a real 3-step CG relaxation (isolated O2, converges to
#                    max force 0.044366 eV/Ang, below MD.MaxForceTol 0.05).
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
echo "--- Starting tester for STB-STATUS (item 5.6) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR/singlepoint" "$TEST_DIR/aimd" "$TEST_DIR/relax"
cp "$FIXTURE_DIR/singlepoint/calc.out" "$FIXTURE_DIR/singlepoint/calc.fdf" "$TEST_DIR/singlepoint/"
cp "$FIXTURE_DIR/aimd/aimd.out" "$FIXTURE_DIR/aimd/aimd.fdf" "$FIXTURE_DIR/aimd/aimd.ANI" "$FIXTURE_DIR/aimd/aimd.XV" "$TEST_DIR/aimd/"
cp "$FIXTURE_DIR/relax/relax.out" "$FIXTURE_DIR/relax/relax.fdf" "$TEST_DIR/relax/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Single-point folder ---
echo -e "\n--- Testing single-point folder (auto-detects SystemLabel=Sn3O4) ---"
stb-status --no-intro --path singlepoint > log_singlepoint.txt 2>&1
check_contains "SystemLabel       : Sn3O4" log_singlepoint.txt
check_contains "Run type          : single-point" log_singlepoint.txt
check_contains "SCF converged     : yes" log_singlepoint.txt
check_contains "Free energy (eV)  : -4054.427104" log_singlepoint.txt


# --- 3. AIMD folder ---
echo -e "\n--- Testing AIMD folder (5 steps, final Temp_ion) ---"
stb-status --no-intro --path aimd > log_aimd.txt 2>&1
check_contains "SystemLabel       : aimd" log_aimd.txt
check_contains "Run type          : aimd" log_aimd.txt
check_contains "Steps             : 5" log_aimd.txt
check_contains "Final temperature : 9917.825 K" log_aimd.txt
check_contains "AIMD trajectory     : .*yes" log_aimd.txt


# --- 4. Relaxation folder ---
echo -e "\n--- Testing relaxation folder (3 CG opt. moves) ---"
stb-status --no-intro --path relax > log_relax.txt 2>&1
check_contains "SystemLabel       : relax" log_relax.txt
check_contains "Run type          : relaxation" log_relax.txt
check_contains "Steps             : 3" log_relax.txt
check_contains "Max force (eV/Ang): 0.044366" log_relax.txt


# --- 5. --label override ---
echo -e "\n--- Testing --label override ---"
stb-status --no-intro --path singlepoint --label Sn3O4 > log_labelover.txt 2>&1
check_contains "SystemLabel       : Sn3O4" log_labelover.txt


# --- 6. Batch mode (glob matching all 3 folders) ---
echo -e "\n--- Testing batch mode (--path matching all 3 folders) ---"
stb-status --no-intro --path "*" > log_batch.txt 2>&1
check_contains "3 director" log_batch.txt
check_contains "singlepoint" log_batch.txt
check_contains "aimd" log_batch.txt
check_contains "relax" log_batch.txt
check_contains "Folders inspected : 3" log_batch.txt


# --- 7. --save-report ---
echo -e "\n--- Testing --save-report ---"
rm -f stb_status_report.txt
stb-status --no-intro --path aimd --save-report > log_savereport.txt 2>&1
check_success stb_status_report.txt
check_contains "STB-STATUS REPORT" stb_status_report.txt
rm -f stb_status_report.txt


# --- 8. Missing directory ---
echo -e "\n--- Testing missing directory ---"
stb-status --no-intro --path does_not_exist > log_missing.txt 2>&1
check_exit_code $? 2
check_contains "is not a directory" log_missing.txt


# --- 9. Interactive path (stb-suite, shortcut 5.6) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 5.6) ---"
printf '5.6\naimd\n\n\n0\n' | stb-suite > log_interactive.txt 2>&1
check_contains "Run type          : aimd" log_interactive.txt
check_contains "Status check complete" log_interactive.txt


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
