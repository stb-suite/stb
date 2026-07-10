#!/bin/bash

# --- Setup ---
# Smoke test for stb-convergenceAnalysis (Convergence Test Analysis, item 4.5.2)
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
echo "--- Starting tester for STB-ConvergenceAnalysis (item 4.5.2) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR/convergence_runs"
for d in 150.0000 200.0000 250.0000 300.0000; do
    cp -r "$FIXTURE_DIR/convergence_meshcutoff_$d" "$TEST_DIR/convergence_runs/"
done
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Convergence analysis (synthetic FreeEng-only calc.out fixtures, same
#     minimal style already used by 1-strain/2-elastic/3-cohesive) ---
echo -e "\n--- Testing analysis of a 4-point meshcutoff sweep ---"
stb-convergenceAnalysis --dir convergence_runs --tolerance 0.001 --no-intro > log_analysis.txt 2>&1
check_contains "Runs found:.*4" log_analysis.txt
check_contains "Converged at:.*meshcutoff = 300.0000" log_analysis.txt
check_success convergence_curve.dat
check_success convergence_report.txt
check_contains "Converged at: meshcutoff = 300.0000" convergence_report.txt


# --- 3. Tighter tolerance moves the converged point ---
echo -e "\n--- Testing a tolerance no point satisfies ---"
rm -f convergence_curve.dat convergence_report.txt
stb-convergenceAnalysis --dir convergence_runs --tolerance 0.0000001 --no-intro > log_tight.txt 2>&1
check_contains "No point met the" log_tight.txt
check_contains "Not converged within the tested range" convergence_report.txt


# --- 4. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: empty directory (no convergence_* folders)"
mkdir -p empty_runs
stb-convergenceAnalysis --dir empty_runs --no-intro > log_empty.txt 2>&1
check_exit_code $? 1
check_contains "No 'convergence_\*' folders found" log_empty.txt

echo "Testing: missing directory"
stb-convergenceAnalysis --dir does_not_exist --no-intro > log_missing_dir.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_dir.txt

echo "Testing: folder missing calc.out (no valid data)"
mkdir -p bad_runs/convergence_meshcutoff_100.0000
cp convergence_runs/convergence_meshcutoff_150.0000/structure.fdf bad_runs/convergence_meshcutoff_100.0000/
stb-convergenceAnalysis --dir bad_runs --no-intro > log_bad_runs.txt 2>&1
check_exit_code $? 1
check_contains "No valid data found" log_bad_runs.txt

echo "Testing: --version"
stb-convergenceAnalysis --version > log_version.txt 2>&1
check_contains "stb-convergenceAnalysis" log_version.txt

echo "Testing: --help documents tolerance/dir"
stb-convergenceAnalysis --help > log_help.txt 2>&1
check_contains "tolerance" log_help.txt
check_contains "dir" log_help.txt


# --- 5. Interactive path (stb-suite, shortcut 3.5.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 3.5.2) ---"

echo "Testing: navigate 3.5.2 -> defaults -> quit"
rm -f convergence_curve.dat convergence_report.txt
printf '4.5.2\nconvergence_runs\n\n0.001\n0\n' | stb-suite > log_menu.txt 2>&1
check_contains "Converged at:.*meshcutoff = 300.0000" log_menu.txt
check_success convergence_curve.dat


popd > /dev/null

# --- 6. Summary ---
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
