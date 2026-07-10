#!/bin/bash

# --- Setup ---
# Smoke test for stb-hubbarduAnalysis (Hubbard U Linear Response, Stage 3: Analysis, item 4.7.3)
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
# The fixture 'hubbardu_runs/' is committed with synthetic .out files whose
# occupation-vs-alpha response is an EXACT known analytic linear function
# (chi=0.30 1/eV screened, chi0=0.20 1/eV bare, N0=5.0), so the fitted U has a
# known ground truth: U = 1/chi0 - 1/chi = 1/0.20 - 1/0.30 = 1.666667 eV. The
# frozen branch has its OWN frozen_alpha_0.0000 point (not a reuse of
# reference), matching what stb-hubbarduAlphas now actually generates.
# scf_alpha_-0.1500 is deliberately marked SCF_NOT_CONV to exercise the
# per-run convergence-status reporting.
echo "--- Starting tester for STB-HubbarduAnalysis stage 3: analysis (item 4.7.3) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp -r "$FIXTURE_DIR/hubbardu_runs" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Analysis of the known-ground-truth fixture ---
echo -e "\n--- Testing analysis recovers the known U ---"
stb-hubbarduAnalysis --dir hubbardu_runs --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "0.300000 1/eV" log_basic.txt
check_contains "0.200000 1/eV" log_basic.txt
check_contains "(R^2=1.0000)" log_basic.txt
check_contains "1.6667 eV" log_basic.txt
check_success Mn_LDAU.fdf

echo "Testing: per-run SCF convergence status is reported"
check_contains "SCF: converged" log_basic.txt
check_contains "SCF: NOT CONVERGED, used anyway" log_basic.txt

echo "Testing: the final .fdf has no LDAU.PotentialShift (production block, not a perturbation run)"
if grep -q "PotentialShift" Mn_LDAU.fdf; then
    echo -e "   -> ${RED}Failed:${NC} final block should NOT contain PotentialShift"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} final block has no PotentialShift line"
    PASS=$((PASS+1))
fi
check_contains "Mn   1" Mn_LDAU.fdf
check_contains "n=3    2" Mn_LDAU.fdf
check_contains "1.667    0.000" Mn_LDAU.fdf

echo "Testing: literature reference comparison is printed for a known species (Mn)"
check_contains "Literature reference" log_basic.txt
check_contains "3.90 eV" log_basic.txt


# --- 3. --output override ---
echo -e "\n--- Testing -o/--output ---"
stb-hubbarduAnalysis --dir hubbardu_runs -o custom_LDAU.fdf --no-intro > log_output.txt 2>&1
check_exit_code $? 0
check_success custom_LDAU.fdf


# --- 4. R^2 warning path (noisy/non-linear synthetic data) ---
echo -e "\n--- Testing the low-R^2 warning ---"
mkdir -p noisy_runs
python3 -c "
import json, os
manifest = {'species': 'Mn', 'n': 3, 'l': 2, 'j': 0.0, 'label': 'siesta', 'runs': {
    'reference': {'kind': 'reference', 'alpha': 0.0},
    'scf_alpha_0.1000': {'kind': 'scf', 'alpha': 0.1},
    'scf_alpha_-0.1000': {'kind': 'scf', 'alpha': -0.1},
    'frozen_alpha_0.0000': {'kind': 'frozen', 'alpha': 0.0},
    'frozen_alpha_0.1000': {'kind': 'frozen', 'alpha': 0.1},
    'frozen_alpha_-0.1000': {'kind': 'frozen', 'alpha': -0.1},
}}
with open('noisy_runs/run_manifest.json', 'w') as f:
    json.dump(manifest, f)
data = {'reference': 5.0, 'scf_alpha_0.1000': 5.2, 'scf_alpha_-0.1000': 5.15,
        'frozen_alpha_0.0000': 5.0, 'frozen_alpha_0.1000': 5.02, 'frozen_alpha_-0.1000': 4.99}
for folder, occ in data.items():
    os.makedirs(f'noisy_runs/{folder}', exist_ok=True)
    with open(f'noisy_runs/{folder}/calc.out', 'w') as f:
        f.write('hubbard_term: Total projector shell\n')
        f.write(f'Occupations:   {occ:.6f}\n')
        f.write('SCF Convergence by DM+H criterion\n')
"
stb-hubbarduAnalysis --dir noisy_runs --no-intro > log_noisy.txt 2>&1
check_exit_code $? 0
check_contains "below the 0.98 tolerance" log_noisy.txt


# --- 5. Error and robustness cases ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: missing directory"
stb-hubbarduAnalysis --dir does_not_exist --no-intro > log_missing_dir.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_dir.txt

echo "Testing: directory without run_manifest.json"
mkdir -p not_a_hubbardu_dir
stb-hubbarduAnalysis --dir not_a_hubbardu_dir --no-intro > log_no_manifest.txt 2>&1
check_exit_code $? 1
check_contains "run_manifest.json' not found" log_no_manifest.txt

echo "Testing: folders with unparseable/missing .out (not enough valid points)"
mkdir -p bad_runs/reference bad_runs/scf_alpha_0.1000
python3 -c "
import json
manifest = {'species': 'Fe', 'n': 3, 'l': 2, 'j': 0.0, 'label': 'siesta', 'runs': {
    'reference': {'kind': 'reference', 'alpha': 0.0},
    'scf_alpha_0.1000': {'kind': 'scf', 'alpha': 0.1},
}}
with open('bad_runs/run_manifest.json', 'w') as f:
    json.dump(manifest, f)
"
echo "no occupation line here" > bad_runs/reference/calc.out
echo "no occupation line here" > bad_runs/scf_alpha_0.1000/calc.out
stb-hubbarduAnalysis --dir bad_runs --no-intro > log_bad_runs.txt 2>&1
check_exit_code $? 1
check_contains "Not enough valid runs" log_bad_runs.txt

echo "Testing: --version"
stb-hubbarduAnalysis --version > log_version.txt 2>&1
check_contains "stb-hubbarduAnalysis" log_version.txt

echo "Testing: --help documents dir/file/r2-tolerance"
stb-hubbarduAnalysis --help > log_help.txt 2>&1
check_contains "dir" log_help.txt
check_contains "r2-tolerance" log_help.txt


# --- 6. Interactive path (stb-suite, shortcut 4.7.3) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.7.3) ---"

echo "Testing: navigate 4.7.3 -> defaults -> quit"
rm -f Mn_LDAU.fdf
printf '4.7.3\nhubbardu_runs\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_contains "1.6667 eV" log_menu.txt
check_success Mn_LDAU.fdf


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
