#!/bin/bash

# --- Setup ---
# Smoke test for stb-hirshfeldAnalysis (Hirshfeld-I Stage 3: Analysis,
# item 4.20.3). Chains real stb-hirshfeldPrep + stb-hirshfeldIons (Stages
# 1-2, tested separately) runs on the ../prep fixture, fabricating each
# neutral/<species>/'s and ion/<species>/{cation,anion}/'s own .RHO via
# ../make_synthetic_rho.py (standing in for "SIESTA has already been run
# there") before running stb-hirshfeldAnalysis for real.
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREP_DIR="$(cd "$FIXTURE_DIR/../prep" && pwd)"
IONS_DIR="$(cd "$FIXTURE_DIR/../ions" && pwd)"
GEN_SCRIPT="$FIXTURE_DIR/../make_synthetic_rho.py"
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

# Runs stb-hirshfeldPrep (writes combined/ + neutral/<species>/ itself,
# nothing pre-computed needed) + fabricates combined/neutral .RHOs +
# stb-hirshfeldIons (writes BOTH cation/anion per species) + fabricates
# cation/anion .RHOs -- the full chain stb-hirshfeldAnalysis needs.
make_prep_ions_and_rho() {
    rm -rf hirshfeld_study
    stb-hirshfeldPrep -s structure.fdf --calc calc.fdf -p . -O hirshfeld_study --no-intro \
        > /dev/null 2>&1
    python3 "$GEN_SCRIPT" hirshfeld_study/combined/structure.fdf \
        hirshfeld_study/combined/siesta.RHO 24 "0.8,1.4" "1.4,1.0" > /dev/null 2>&1
    python3 "$GEN_SCRIPT" hirshfeld_study/neutral/C/structure.fdf \
        hirshfeld_study/neutral/C/siesta.RHO 18 "0.8" "1.4" > /dev/null 2>&1
    python3 "$GEN_SCRIPT" hirshfeld_study/neutral/O/structure.fdf \
        hirshfeld_study/neutral/O/siesta.RHO 18 "1.4" "1.0" > /dev/null 2>&1
    stb-hirshfeldIons -O hirshfeld_study -p . --no-intro > /dev/null 2>&1
    python3 "$GEN_SCRIPT" hirshfeld_study/ions/C/cation/structure.fdf \
        hirshfeld_study/ions/C/cation/siesta.RHO 18 "0.6" "1.3" > /dev/null 2>&1
    python3 "$GEN_SCRIPT" hirshfeld_study/ions/C/anion/structure.fdf \
        hirshfeld_study/ions/C/anion/siesta.RHO 18 "1.0" "1.5" > /dev/null 2>&1
    python3 "$GEN_SCRIPT" hirshfeld_study/ions/O/cation/structure.fdf \
        hirshfeld_study/ions/O/cation/siesta.RHO 18 "1.1" "0.9" > /dev/null 2>&1
    python3 "$GEN_SCRIPT" hirshfeld_study/ions/O/anion/structure.fdf \
        hirshfeld_study/ions/O/anion/siesta.RHO 18 "1.7" "1.1" > /dev/null 2>&1
}


# --- 1. Preparation ---
echo "--- Starting tester for STB-HIRSHFELDANALYSIS (item 4.20.3) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$PREP_DIR/structure.fdf" "$PREP_DIR/calc.fdf" "$PREP_DIR/C.psf" "$PREP_DIR/O.psf" "$TEST_DIR/"
cp "$IONS_DIR/structure_mixed.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Missing Stage-2 output guard ---
echo -e "\n--- Testing the missing-Stage-2-output guard ---"
rm -rf hirshfeld_study
stb-hirshfeldAnalysis -O hirshfeld_study --no-intro > log_no_stage2.txt 2>&1
check_exit_code $? 1
check_contains "run stb-hirshfeldPrep then stb-hirshfeldIons first" log_no_stage2.txt


# --- 3. Default run: live iteration output, convergence + per-atom report ---
echo -e "\n--- Testing default run (live iteration output, convergence, per-atom report) ---"
make_prep_ions_and_rho
stb-hirshfeldAnalysis -O hirshfeld_study --no-intro > log_default.txt 2>&1
check_exit_code $? 0

echo "Testing: each round prints a live max|Delta q| line naming the worst atom/species,"
echo "         its current charge/population + whether it's leaning cation/anion-like,"
echo "         elapsed time, and convergence status"
check_contains "ITERATING HIRSHFELD-I" log_default.txt
check_contains "Iteration  1/20: max|Delta q| = " log_default.txt
check_contains "-- atom #" log_default.txt
check_contains "-like, pop\." log_default.txt
check_contains "s elapsed --" log_default.txt

check_contains "CONVERGENCE HISTORY" log_default.txt
check_contains "Converged after" log_default.txt
check_contains "SIMPLE HIRSHFELD vs. HIRSHFELD-I" log_default.txt
check_contains "PER-SPECIES SUMMARY (HIRSHFELD-I)" log_default.txt
check_contains "Cation-like | Anion-like" log_default.txt
check_contains "Iterations run  :" log_default.txt
check_contains "Iterating time  :" log_default.txt


# --- 4. --tol / --max-iter flags ---
echo -e "\n--- Testing --tol/--max-iter flags ---"
stb-hirshfeldAnalysis -O hirshfeld_study --tol 0.05 --max-iter 3 --no-intro > log_tol.txt 2>&1
check_exit_code $? 0
check_contains "Tolerance      : 0.05 e-" log_tol.txt
check_contains "Max iterations : 3" log_tol.txt


# --- 4b. --ref (explicit reference file, e.g. not named with a .out extension) ---
echo -e "\n--- Testing --ref (Z_val detection from an explicitly-named, non-.out reference file) ---"
cat > hirshfeld_study/combined/run_log.txt << 'ZVALEOF'
atom: Called for C(Z=6)
Vna: chval, zval:    4.00000   4.00000
atom: Called for O(Z=8)
Vna: chval, zval:    6.00000   6.00000
ZVALEOF
stb-hirshfeldAnalysis -O hirshfeld_study --ref hirshfeld_study/combined/run_log.txt \
    --no-intro > log_ref.txt 2>&1
check_exit_code $? 0
check_contains "Total charge" log_ref.txt


# --- 5. --save-report ---
echo -e "\n--- Testing --save-report ---"
rm -f hirshfeld_analysis_report.txt
stb-hirshfeldAnalysis -O hirshfeld_study --save-report --no-intro > log_savereport.txt 2>&1
check_success hirshfeld_analysis_report.txt
check_contains "STB-HIRSHFELDANALYSIS REPORT" hirshfeld_analysis_report.txt
rm -f hirshfeld_analysis_report.txt


# --- 6. End-to-end proof: atoms of the SAME species converge to OPPOSITE ion
# states (the actual payoff of writing both cation and anion per species) ---
echo -e "\n--- Testing the full pipeline on a genuinely mixed-sign species end-to-end ---"
rm -rf hirshfeld_study_mixed mixed_combined
mkdir -p mixed_combined
cp structure_mixed.fdf mixed_combined/structure.fdf
cp calc.fdf mixed_combined/
stb-hirshfeldPrep -s mixed_combined/structure.fdf --calc mixed_combined/calc.fdf -p . \
    -O hirshfeld_study_mixed --mesh-cutoff 400 --vacuum 18 --no-intro > /dev/null 2>&1
python3 "$GEN_SCRIPT" hirshfeld_study_mixed/combined/structure.fdf \
    hirshfeld_study_mixed/combined/siesta.RHO 40 "0.02,3.0,1.0" "0.5,0.6,1.0" > /dev/null 2>&1
python3 "$GEN_SCRIPT" hirshfeld_study_mixed/neutral/C/structure.fdf \
    hirshfeld_study_mixed/neutral/C/siesta.RHO 24 "0.5" "0.8" > /dev/null 2>&1
python3 "$GEN_SCRIPT" hirshfeld_study_mixed/neutral/O/structure.fdf \
    hirshfeld_study_mixed/neutral/O/siesta.RHO 24 "1.0" "0.8" > /dev/null 2>&1
stb-hirshfeldIons -O hirshfeld_study_mixed -p . --no-intro > /dev/null 2>&1
python3 "$GEN_SCRIPT" hirshfeld_study_mixed/ions/C/cation/structure.fdf \
    hirshfeld_study_mixed/ions/C/cation/siesta.RHO 24 "0.3" "0.7" > /dev/null 2>&1
python3 "$GEN_SCRIPT" hirshfeld_study_mixed/ions/C/anion/structure.fdf \
    hirshfeld_study_mixed/ions/C/anion/siesta.RHO 24 "0.9" "1.0" > /dev/null 2>&1
python3 "$GEN_SCRIPT" hirshfeld_study_mixed/ions/O/cation/structure.fdf \
    hirshfeld_study_mixed/ions/O/cation/siesta.RHO 24 "0.7" "0.8" > /dev/null 2>&1
python3 "$GEN_SCRIPT" hirshfeld_study_mixed/ions/O/anion/structure.fdf \
    hirshfeld_study_mixed/ions/O/anion/siesta.RHO 24 "1.3" "0.9" > /dev/null 2>&1
stb-hirshfeldAnalysis -O hirshfeld_study_mixed --no-intro > log_mixed_analysis.txt 2>&1
check_exit_code $? 0
echo "Testing: the converged per-species table shows species 'C' with BOTH a cation-like"
echo "         and an anion-like atom, highlighted -- exactly what a single-sign-per-species"
echo "         approach could never represent correctly"
check_contains "C    | 2 | 1           | 1          " log_mixed_analysis.txt
check_contains "atoms genuinely leaning toward BOTH ion states" log_mixed_analysis.txt


# --- 7. Interactive path (stb-suite, shortcut 4.20.3) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.20.3) ---"
printf '4.20.3\nhirshfeld_study\n\n\n\nn\n\n0\n' | stb-suite > log_interactive.txt 2>&1
check_contains "Hirshfeld-I Analysis (Stage 3) complete" log_interactive.txt


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
