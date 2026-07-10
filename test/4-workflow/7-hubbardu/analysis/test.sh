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
# per-run convergence-status reporting for a GENUINE convergence problem.
# Every frozen_alpha_* .out is marked SCF_NOT_CONV too, matching real SIESTA
# (MaxSCFIterations caps them before real self-consistency, by design) --
# this must be reported as expected/green, not as a warning.
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
check_contains "R^2=1.0000" log_basic.txt
check_contains "1.6667 eV" log_basic.txt
check_success Mn_LDAU.fdf

echo "Testing: per-run SCF convergence status is reported"
check_contains "SCF: converged" log_basic.txt
check_contains "SCF: NOT CONVERGED, used anyway" log_basic.txt

echo "Testing: frozen-density runs that hit SCF_NOT_CONV (expected, by design) are reported as OK, not a warning"
check_contains "SCF: frozen evaluation complete (stopped at MaxSCFIterations, as designed)" log_basic.txt
if grep -A1 "frozen_alpha_0.0000" log_basic.txt | grep -q "NOT CONVERGED, used anyway"; then
    echo -e "   -> ${RED}Failed:${NC} a frozen run's expected SCF_NOT_CONV was reported as a real problem"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} frozen run's SCF_NOT_CONV not flagged as a real problem"
    PASS=$((PASS+1))
fi

echo "Testing: a convergence summary is reported, split between scf/reference and frozen runs"
check_contains "Convergence summary:" log_basic.txt
check_contains "14/14 run(s) read" log_basic.txt
check_contains "of the 7 self-consistent (scf/reference) run(s), 6/7 fully converged, 1 NOT converged" log_basic.txt
check_contains "Frozen-density runs: 7/7 completed as expected" log_basic.txt

echo "Testing: fit intercepts are reported"
check_contains "intercept=5.000000" log_basic.txt

echo "Testing: a response plot (data + gplot) is generated"
check_success Mn_hubbardu_response.dat
check_success Mn_hubbardu_response.gplot
check_contains "index 0" Mn_hubbardu_response.gplot
check_contains "index 1" Mn_hubbardu_response.gplot
check_contains "f_scf(x) = 0.30000000" Mn_hubbardu_response.gplot
check_contains "f_frozen(x) = 0.20000000" Mn_hubbardu_response.gplot

echo "Testing: no spurious sign/negative-U warnings on this well-behaved fixture"
if grep -qE "OPPOSITE signs|Computed U is NEGATIVE" log_basic.txt; then
    echo -e "   -> ${RED}Failed:${NC} unexpected sign warning on well-behaved data"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no sign warnings on well-behaved data"
    PASS=$((PASS+1))
fi

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


# --- 4b. Opposite-sign chi/chi0 -> negative U warnings ---
echo -e "\n--- Testing the opposite-sign / negative-U warnings ---"
mkdir -p opposite_sign_runs
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
with open('opposite_sign_runs/run_manifest.json', 'w') as f:
    json.dump(manifest, f)
# scf increases with alpha (chi>0), frozen DECREASES with alpha (chi0<0) -- unphysical
data = {'reference': 5.0, 'scf_alpha_0.1000': 5.03, 'scf_alpha_-0.1000': 4.97,
        'frozen_alpha_0.0000': 5.0, 'frozen_alpha_0.1000': 4.98, 'frozen_alpha_-0.1000': 5.02}
for folder, occ in data.items():
    os.makedirs(f'opposite_sign_runs/{folder}', exist_ok=True)
    with open(f'opposite_sign_runs/{folder}/calc.out', 'w') as f:
        f.write('hubbard_term: Total projector shell\n')
        f.write(f'Occupations:   {occ:.6f}\n')
        f.write('SCF Convergence by DM+H criterion\n')
"
stb-hubbarduAnalysis --dir opposite_sign_runs --no-intro > log_opposite_sign.txt 2>&1
check_exit_code $? 0
check_contains "OPPOSITE signs" log_opposite_sign.txt
check_contains "Computed U is NEGATIVE" log_opposite_sign.txt


# --- 4c. Intercept-deviation warning (reference is an outlier vs. its own branch's line) ---
echo -e "\n--- Testing the intercept-deviation warning ---"
mkdir -p intercept_outlier_runs
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
with open('intercept_outlier_runs/run_manifest.json', 'w') as f:
    json.dump(manifest, f)
# reference sits well off the line the other two scf points define
data = {'reference': 5.10, 'scf_alpha_0.1000': 5.03, 'scf_alpha_-0.1000': 4.97,
        'frozen_alpha_0.0000': 5.0, 'frozen_alpha_0.1000': 5.02, 'frozen_alpha_-0.1000': 4.98}
for folder, occ in data.items():
    os.makedirs(f'intercept_outlier_runs/{folder}', exist_ok=True)
    with open(f'intercept_outlier_runs/{folder}/calc.out', 'w') as f:
        f.write('hubbard_term: Total projector shell\n')
        f.write(f'Occupations:   {occ:.6f}\n')
        f.write('SCF Convergence by DM+H criterion\n')
"
stb-hubbarduAnalysis --dir intercept_outlier_runs --no-intro > log_intercept.txt 2>&1
check_exit_code $? 0
check_contains "differs from the actually-measured alpha=0 occupation" log_intercept.txt


# --- 4c-bis. Exactly 2 points per branch warns that R^2 is trivially 1.0 ---
echo -e "\n--- Testing the 2-points-per-branch R^2 warning ---"
mkdir -p two_point_runs
python3 -c "
import json, os
manifest = {'species': 'Mn', 'n': 3, 'l': 2, 'j': 0.0, 'label': 'siesta', 'runs': {
    'reference': {'kind': 'reference', 'alpha': 0.0},
    'scf_alpha_0.1000': {'kind': 'scf', 'alpha': 0.1},
    'frozen_alpha_0.0000': {'kind': 'frozen', 'alpha': 0.0},
    'frozen_alpha_0.1000': {'kind': 'frozen', 'alpha': 0.1},
}}
with open('two_point_runs/run_manifest.json', 'w') as f:
    json.dump(manifest, f)
data = {'reference': 5.0, 'scf_alpha_0.1000': 5.03, 'frozen_alpha_0.0000': 5.0, 'frozen_alpha_0.1000': 4.98}
for folder, occ in data.items():
    os.makedirs(f'two_point_runs/{folder}', exist_ok=True)
    with open(f'two_point_runs/{folder}/calc.out', 'w') as f:
        f.write('hubbard_term: Total projector shell\n')
        f.write(f'Occupations:   {occ:.6f}\n')
        f.write('SCF Convergence by DM+H criterion\n')
"
stb-hubbarduAnalysis --dir two_point_runs --no-intro > log_two_point.txt 2>&1
check_exit_code $? 0
check_contains "only 2 points" log_two_point.txt


# --- 4d. Mismatched --frozen-iterations across the sweep is flagged ---
echo -e "\n--- Testing the frozen_iterations-mismatch warning ---"
mkdir -p mixed_iter_runs
python3 -c "
import json, os
manifest = {'species': 'Mn', 'n': 3, 'l': 2, 'j': 0.0, 'label': 'siesta', 'runs': {
    'reference': {'kind': 'reference', 'alpha': 0.0},
    'scf_alpha_0.1000': {'kind': 'scf', 'alpha': 0.1},
    'scf_alpha_-0.1000': {'kind': 'scf', 'alpha': -0.1},
    'frozen_alpha_0.0000': {'kind': 'frozen', 'alpha': 0.0, 'frozen_iterations': 2},
    'frozen_alpha_0.1000': {'kind': 'frozen', 'alpha': 0.1, 'frozen_iterations': 3},
    'frozen_alpha_-0.1000': {'kind': 'frozen', 'alpha': -0.1, 'frozen_iterations': 2},
}}
with open('mixed_iter_runs/run_manifest.json', 'w') as f:
    json.dump(manifest, f)
data = {'reference': 5.0, 'scf_alpha_0.1000': 5.03, 'scf_alpha_-0.1000': 4.97,
        'frozen_alpha_0.0000': 5.0, 'frozen_alpha_0.1000': 4.98, 'frozen_alpha_-0.1000': 5.02}
for folder, occ in data.items():
    os.makedirs(f'mixed_iter_runs/{folder}', exist_ok=True)
    with open(f'mixed_iter_runs/{folder}/calc.out', 'w') as f:
        f.write('hubbard_term: Total projector shell\n')
        f.write(f'Occupations:   {occ:.6f}\n')
        f.write('SCF Convergence by DM+H criterion\n')
"
stb-hubbarduAnalysis --dir mixed_iter_runs --no-intro > log_mixed_iter.txt 2>&1
check_exit_code $? 0
check_contains "used different --frozen-iterations values" log_mixed_iter.txt


# --- 4f. check_scf_converged uses the LAST marker, not "any marker present" ---
echo -e "\n--- Testing that a failed-then-retried-and-converged .out is read as converged ---"
mkdir -p retry_runs/reference
python3 -c "
import json
manifest = {'species': 'Mn', 'n': 3, 'l': 2, 'j': 0.0, 'label': 'siesta', 'runs': {
    'reference': {'kind': 'reference', 'alpha': 0.0},
    'scf_alpha_0.1000': {'kind': 'scf', 'alpha': 0.1},
    'frozen_alpha_0.0000': {'kind': 'frozen', 'alpha': 0.0},
    'frozen_alpha_0.1000': {'kind': 'frozen', 'alpha': 0.1},
}}
with open('retry_runs/run_manifest.json', 'w') as f:
    json.dump(manifest, f)
"
for d in scf_alpha_0.1000 frozen_alpha_0.0000 frozen_alpha_0.1000; do mkdir -p retry_runs/$d; done
for d in reference scf_alpha_0.1000 frozen_alpha_0.0000 frozen_alpha_0.1000; do
    cat > retry_runs/$d/calc.out << 'EOF'
hubbard_term: Total projector shell
Occupations:   4.500000
SCF_NOT_CONV: SCF did not converge in maximum number of steps (required).
hubbard_term: Total projector shell
Occupations:   5.000000
SCF Convergence by DM+H criterion
EOF
done
stb-hubbarduAnalysis --dir retry_runs --no-intro > log_retry.txt 2>&1
check_contains "SCF: converged" log_retry.txt
if grep -q "NOT CONVERGED" log_retry.txt; then
    echo -e "   -> ${RED}Failed:${NC} a converged-on-retry .out should not be reported as not-converged"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} the LAST marker (converged) wins over an earlier failed attempt"
    PASS=$((PASS+1))
fi


# --- 4g. Non-polarized (nspin=1) 'Occupations:' lines print TWO numbers, and
#     BOTH equal the per-spin occupation (SIESTA's sum(oc) degenerates to
#     oc(1) when there's only one spin channel) -- verified against
#     Src/dftu.F: occu() is computed as Dij*Sik*Sjk/(3-nspin), i.e. divided
#     by 2 for nspin=1, since Dscf already holds the full (both-spin)
#     density in its single channel. So the true total is 2x the printed
#     value, not the raw last number -- a real bug found analyzing a real
#     non-polarized Mn run that silently halved the occupation, doubling
#     the computed U. ---
echo -e "\n--- Testing that non-polarized (2-number) Occupations lines are doubled to get the true total ---"
mkdir -p nonpolarized_runs
python3 -c "
import os, json
manifest = {'species': 'Mn', 'n': 3, 'l': 2, 'j': 0.0, 'label': 'siesta', 'runs': {
    'reference': {'kind': 'reference', 'alpha': 0.0},
    'scf_alpha_0.1000': {'kind': 'scf', 'alpha': 0.1},
    'scf_alpha_-0.1000': {'kind': 'scf', 'alpha': -0.1},
    'frozen_alpha_0.0000': {'kind': 'frozen', 'alpha': 0.0},
    'frozen_alpha_0.1000': {'kind': 'frozen', 'alpha': 0.1},
    'frozen_alpha_-0.1000': {'kind': 'frozen', 'alpha': -0.1},
}}
with open('nonpolarized_runs/run_manifest.json', 'w') as f:
    json.dump(manifest, f)
# Per-spin (halved) occupations: true totals are 5.0/5.06/4.94.
data = {'reference': 2.5, 'scf_alpha_0.1000': 2.53, 'scf_alpha_-0.1000': 2.47,
        'frozen_alpha_0.0000': 2.5, 'frozen_alpha_0.1000': 2.51, 'frozen_alpha_-0.1000': 2.49}
for folder, occ in data.items():
    os.makedirs(f'nonpolarized_runs/{folder}', exist_ok=True)
    with open(f'nonpolarized_runs/{folder}/calc.out', 'w') as f:
        f.write('hubbard_term: Total projector shell\n')
        f.write(f'Occupations:   {occ:.6f}   {occ:.6f}\n')
        f.write('SCF Convergence by DM+H criterion\n')
"
stb-hubbarduAnalysis --dir nonpolarized_runs --no-intro > log_nonpolarized.txt 2>&1
check_exit_code $? 0
check_contains "occupation=5.000000" log_nonpolarized.txt
check_contains "occupation=5.060000" log_nonpolarized.txt
check_contains "occupation=4.940000" log_nonpolarized.txt
if grep -q "occupation=2\." log_nonpolarized.txt; then
    echo -e "   -> ${RED}Failed:${NC} a non-polarized run's occupation was NOT doubled (still using the per-spin value)"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} non-polarized occupations were correctly doubled to the true total"
    PASS=$((PASS+1))
fi


# --- 4h. A manifest missing a required field is a clean error, not a raw KeyError ---
echo -e "\n--- Testing the manifest-validation guard ---"
mkdir -p old_format_runs/reference
python3 -c "
import json
# Pre-3-stage-split manifest shape: no 'j' field
manifest = {'species': 'Mn', 'n': 3, 'l': 2, 'label': 'siesta', 'runs': {
    'reference': {'kind': 'reference', 'alpha': 0.0},
}}
with open('old_format_runs/run_manifest.json', 'w') as f:
    json.dump(manifest, f)
"
stb-hubbarduAnalysis --dir old_format_runs --no-intro > log_old_format.txt 2>&1
check_exit_code $? 1
check_contains "missing required field" log_old_format.txt


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
printf '4.7.3\nhubbardu_runs\n\n\n\n\n0\n' | stb-suite > log_menu.txt 2>&1
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
