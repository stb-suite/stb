#!/bin/bash

# --- Setup ---
# Smoke test for stb-status (Status Checker, item 3.6).
#
# Fixtures are 3 real, distinctly-typed SIESTA runs:
#   singlepoint/  -- Sn3O4.out (calc.out/calc.fdf), single-point, shared with
#                    the other 3-utils fixtures.
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

check_not_contains() {
    if grep -q -- "$1" "$2" 2>/dev/null; then
        echo -e "   -> ${RED}Failed:${NC} '$1' found in '$2' but should not be"
        FAIL=$((FAIL+1))
    else
        echo -e "   -> ${GREEN}Verified:${NC} '$1' not found in '$2' (as expected)"
        PASS=$((PASS+1))
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
echo "--- Starting tester for STB-STATUS (item 3.6) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR/singlepoint" "$TEST_DIR/aimd" "$TEST_DIR/relax"
cp "$FIXTURE_DIR/singlepoint/calc.out" "$FIXTURE_DIR/singlepoint/calc.fdf" "$TEST_DIR/singlepoint/"
cp "$FIXTURE_DIR/aimd/aimd.out" "$FIXTURE_DIR/aimd/aimd.fdf" "$FIXTURE_DIR/aimd/aimd.ANI" "$FIXTURE_DIR/aimd/aimd.XV" "$TEST_DIR/aimd/"
cp "$FIXTURE_DIR/relax/relax.out" "$FIXTURE_DIR/relax/relax.fdf" "$TEST_DIR/relax/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Single-point folder ---
echo -e "\n--- Testing single-point folder (auto-detects SystemLabel=Sn3O4) ---"
rm -f stb_status_report*.txt stb_status_references*.bib
stb-status --no-intro --path singlepoint > log_singlepoint.txt 2>&1
check_exit_code $? 0
check_contains "SystemLabel       : Sn3O4" log_singlepoint.txt
check_contains "State             : .*FINISHED" log_singlepoint.txt
check_contains "Run type          : single point" log_singlepoint.txt
check_contains "Converged steps   : 1 of 1" log_singlepoint.txt
check_contains "Total energy      : -4054.427104 eV" log_singlepoint.txt
check_contains "Mesh cutoff       : 320 Ry required, 325.4 Ry used" log_singlepoint.txt
check_contains "Verdict           : FINISHED -- result usable" log_singlepoint.txt
echo "Testing: the report is always written, in the current folder, with its references file"
check_success stb_status_report.txt
check_contains "STB-STATUS REPORT" stb_status_report.txt
check_success stb_status_references.bib


# --- 3. AIMD folder ---
echo -e "\n--- Testing AIMD folder (5 MD steps, per-step table) ---"
stb-status --no-intro --path aimd > log_aimd.txt 2>&1
check_contains "SystemLabel       : aimd" log_aimd.txt
check_contains "Run type          : MD (Verlet MD run)" log_aimd.txt
check_contains "MD steps          : 5 of 5" log_aimd.txt
check_contains "5    | 6       | -857.5875 | 3.5559" log_aimd.txt
check_contains "MD/relaxation trajectory | aimd.ANI | stb-ani2traj (3.5)" log_aimd.txt
echo "Testing: a second run never overwrites the first report (numbered _1)"
check_success stb_status_report_1.txt


# --- 4. Relaxation folder ---
echo -e "\n--- Testing relaxation folder (2 CG moves after the input geometry) ---"
stb-status --no-intro --path relax > log_relax.txt 2>&1
check_contains "SystemLabel       : relax" log_relax.txt
check_contains "Run type          : relaxation (CG coord. optimization)" log_relax.txt
check_contains "Moves             : 2 of at most 10 (plus the input geometry)" log_relax.txt
check_contains "Relaxation        : converged ('Relaxed atomic coordinates')" log_relax.txt
check_contains "Max force (last)  : 0.0444 eV/Ang vs tolerance 0.05 (ok)" log_relax.txt
check_contains "Max force         : 0.044366 eV/Ang" log_relax.txt


# --- 5. --label override ---
echo -e "\n--- Testing --label override ---"
stb-status --no-intro --path singlepoint --label Sn3O4 > log_labelover.txt 2>&1
check_contains "SystemLabel       : Sn3O4" log_labelover.txt


# --- 6. Batch mode (glob matching all 3 folders) ---
echo -e "\n--- Testing batch mode (--path matching all 3 folders) ---"
stb-status --no-intro --path "*" > log_batch.txt 2>&1
check_contains "Target            : \* -- 3 folder(s)" log_batch.txt
check_contains "\[1\] FOLDERS" log_batch.txt
check_contains "aimd        | FINISHED | MD " log_batch.txt
check_contains "relax       | FINISHED | relaxation (relaxed) | 2 " log_batch.txt
check_contains "singlepoint | FINISHED | single point " log_batch.txt
check_contains "States            : 3 FINISHED" log_batch.txt
check_contains "Folders inspected : 3" log_batch.txt


# --- 7. --save-report is a hidden, ignored option (old scripts keep working) ---
echo -e "\n--- Testing the old --save-report flag ---"
stb-status --no-intro --path aimd --save-report > log_savereport.txt 2>&1
check_exit_code $? 0
check_contains "Status check complete" log_savereport.txt


# --- 8. Missing directory ---
echo -e "\n--- Testing missing directory ---"
stb-status --no-intro --path does_not_exist > log_missing.txt 2>&1
check_exit_code $? 2
check_contains "is not a directory" log_missing.txt


# --- 9. Interactive path (stb-suite, shortcut 3.6) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 3.6) ---"
# Prompts in order: folder, recursive, label, CSV file
printf '3.6\naimd\nn\n\n\n\n0\n' | stb-suite > log_interactive.txt 2>&1
check_contains "Run type          : MD (Verlet MD run)" log_interactive.txt
check_contains "Status check complete" log_interactive.txt
check_not_contains "Traceback" log_interactive.txt


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
