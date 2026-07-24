#!/bin/bash

# --- Setup ---
# Smoke test for stb-sqs (Special Quasirandom Structure Generator, item 2.6)
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
echo "--- Starting tester for STB-SQS (item 2.6) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/fcc_ni.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. NiFe SQS (Monte Carlo, capped steps for speed) ---
echo -e "\n--- Testing NiFe SQS (monte_carlo, --mc-steps 1000) ---"
rm -f sqs.fdf
stb-sqs -f fcc_ni.fdf --sublattice Ni --composition Ni:0.5,Fe:0.5 --scaling 4 --mc-steps 1000 -o sqs.fdf --no-intro > log_sqs.txt 2>&1
check_contains "\[3\] DISORDERED SUBLATTICE" log_sqs.txt
check_contains "Species            : Ni" log_sqs.txt
check_contains "Sites              : 4 (of 4 atoms" log_sqs.txt
check_contains "\[4\] SUPERCELL SIZING (icet)" log_sqs.txt
check_contains "Minimal valid scaling  : 2" log_sqs.txt
check_contains "\[5\] SQS SEARCH" log_sqs.txt
check_contains "Output atoms       : 4" log_sqs.txt
check_contains "Objective function :" log_sqs.txt
check_contains "\[8\] SYMMETRY ANALYSIS (BEFORE / AFTER)" log_sqs.txt
check_contains "Property.*Before.*After" log_sqs.txt
check_success sqs.fdf
check_contains "SQS structure built by stb-sqs" sqs.fdf
check_contains "Disordered sublattice: Ni (4 site" sqs.fdf
check_contains "NumberOfSpecies    2" sqs.fdf
check_contains " 2   26   Fe" sqs.fdf

echo "Testing: stb_sqs_report.txt is NOT written unless --save-report is given"
[ ! -f stb_sqs_report.txt ] && { echo -e "   -> ${GREEN}Verified:${NC} stb_sqs_report.txt absent by default"; PASS=$((PASS+1)); } \
    || { echo -e "   -> ${RED}Failed:${NC} stb_sqs_report.txt was written without --save-report"; FAIL=$((FAIL+1)); }


# --- 3. Round-trip: output is a valid, readable fdf ---
echo -e "\n--- Verifying the generated SQS is valid input for another tool ---"
stb-kgrid --file sqs.fdf --density 0.3 --no-intro > log_roundtrip.txt 2>&1
check_contains "Suggested Monkhorst-Pack grid" log_roundtrip.txt


# --- 4. --save-report ---
echo -e "\n--- Testing --save-report ---"
rm -f sqs.fdf stb_sqs_report.txt
stb-sqs -f fcc_ni.fdf --sublattice Ni --composition Ni:0.5,Fe:0.5 --scaling 4 --mc-steps 1000 -o sqs.fdf --save-report --no-intro > log_report.txt 2>&1
check_exit_code $? 0
check_success stb_sqs_report.txt
check_contains "\[8\] SYMMETRY ANALYSIS (BEFORE / AFTER)" stb_sqs_report.txt
check_contains "Report         : stb_sqs_report.txt" log_report.txt
rm -f stb_sqs_report.txt sqs.fdf


# --- 5. --ml-relax (needs the optional 'ml' extra) ---
echo -e "\n--- Testing --ml-relax (MACE pre-optimization), if the 'ml' extra is available ---"

if python3 -c "import mace" 2>/dev/null; then
    rm -f sqs.fdf stb_sqs_report.txt
    stb-sqs -f fcc_ni.fdf --sublattice Ni --composition Ni:0.5,Fe:0.5 --scaling 4 --mc-steps 500 \
        -o sqs.fdf --ml-relax --save-report --no-intro > log_ml_relax.txt 2>&1
    check_exit_code $? 0
    check_contains "\[6\] ML PRE-RELAXATION (MACE)" stb_sqs_report.txt
    check_success sqs.fdf
    check_contains "ML pre-relaxed with MACE-MP-0" sqs.fdf
    rm -f sqs.fdf stb_sqs_report.txt

    echo "Testing: --ml-relax-cell without --ml-relax is rejected"
    stb-sqs -f fcc_ni.fdf --sublattice Ni --composition Ni:0.5,Fe:0.5 --scaling 4 --ml-relax-cell --no-intro > log_ml_relax_cell_only.txt 2>&1
    check_exit_code $? 2
else
    echo -e "   -> ${YELLOW}SKIPPED${NC}: the optional 'ml' extra is not installed "
    echo "      (pip install stb_suite[ml] to also exercise --ml-relax)."
fi


# --- 6. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: --sublattice not present in structure"
stb-sqs -f fcc_ni.fdf --sublattice Xe --composition Ni:0.5,Fe:0.5 --scaling 4 --no-intro > log_bad_sublattice.txt 2>&1
check_exit_code $? 1
check_contains "not found in the structure" log_bad_sublattice.txt

echo "Testing: --composition not summing to 1.0"
stb-sqs -f fcc_ni.fdf --sublattice Ni --composition Ni:0.5,Fe:0.3 --scaling 4 --no-intro > log_bad_composition.txt 2>&1
check_exit_code $? 1
check_contains "must sum to 1.0" log_bad_composition.txt

echo "Testing: --scaling incompatible with composition (caught by minimal_scaling(), before reaching icet)"
timeout 30 stb-sqs -f fcc_ni.fdf --sublattice Ni --composition Ni:0.3,Fe:0.7 --scaling 4 --mc-steps 500 --no-intro > log_infeasible.txt 2>&1
check_exit_code $? 1
check_contains "is incompatible with composition" log_infeasible.txt
check_contains "multiple of 10" log_infeasible.txt

echo "Testing: --scaling 0"
stb-sqs -f fcc_ni.fdf --sublattice Ni --composition Ni:0.5,Fe:0.5 --scaling 0 --no-intro > log_zero_scaling.txt 2>&1
check_exit_code $? 2

echo "Testing: malformed --cluster-cutoffs"
stb-sqs -f fcc_ni.fdf --sublattice Ni --composition Ni:0.5,Fe:0.5 --scaling 4 --cluster-cutoffs "bad" --no-intro > log_bad_cutoffs.txt 2>&1
check_exit_code $? 1
check_contains "size:shell" log_bad_cutoffs.txt

echo "Testing: missing structure file"
stb-sqs -f does_not_exist.fdf --sublattice Ni --composition Ni:0.5,Fe:0.5 --scaling 4 --no-intro > log_missing.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing.txt

echo "Testing: missing required args"
stb-sqs -f fcc_ni.fdf --no-intro > log_missing_args.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-sqs --version > log_version.txt 2>&1
check_contains "stb-sqs" log_version.txt

echo "Testing: --help documents sublattice/composition/scaling"
stb-sqs --help > log_help.txt 2>&1
check_contains "sublattice" log_help.txt
check_contains "composition" log_help.txt


# --- 7. Interactive path (stb-suite, menu 2.6) ---
echo -e "\n--- Testing the interactive path via stb-suite (menu 2.6) ---"

echo "Testing: navigate 2.6 -> invalid file then valid -> Ni -> Ni:0.5,Fe:0.5 -> scaling 2 -> Monte Carlo ->"
echo "         default temperature -> default output -> no ml-relax -> no save-report -> no view -> quit"
rm -f sqs.fdf
printf '2.6\ndoes_not_exist.fdf\nfcc_ni.fdf\nNi\nNi:0.5,Fe:0.5\n2\n1\n\n\n\n\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_contains "File not found" log_menu.txt
check_contains "Structure written to" log_menu.txt
check_success sqs.fdf


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
