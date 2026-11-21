#!/bin/bash

# --- Setup ---
# Smoke test for stb-siesta2wtb (WanTiBEXOS interface, item 3.4).
#
# Fixtures: a real, finished SIESTA 5.2.2 single point of bulk Sn3O4 with SaveHS
# (calc.fdf %includes structure.fdf; Sn3O4.HSX, Sn3O4.EIG, calc.out). The HSX is
# version 1: it stores E_F and sisl shifts H by -E_F S on reading, so the written
# H is referred to E_F while the header carries -o's E_F -- the report says so
# and warns; with -f 0 the two match.
#
# Checks the report around the writers (system, basis, Hamiltonian range, Fermi
# level, band edges from .EIG), --dry-run, -f 0, the full write with its
# read-back check (H(k), S(k) rebuilt from tb-NP.ham vs sisl, within 5 eV of E_F)
# and the interactive menu path. The full write takes ~1.5 min and a ~200 MB
# tb-NP.ham, removed right after its checks.
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$FIXTURE_DIR/test_files"

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

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

check_absent() {
    if [ ! -e "$1" ]; then
        echo -e "   -> ${GREEN}Verified:${NC} '$1' does not exist (as expected)"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} '$1' exists but should not"
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

check_not_contains() {
    if grep -q -- "$1" "$2" 2>/dev/null; then
        echo -e "   -> ${RED}Failed:${NC} '$1' found in '$2' but should not be"
        FAIL=$((FAIL+1))
    else
        echo -e "   -> ${GREEN}Verified:${NC} '$1' not found in '$2' (as expected)"
        PASS=$((PASS+1))
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
echo "--- Starting tester for STB-SIESTA2WTB (item 3.4) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR"/calc.fdf "$FIXTURE_DIR"/structure.fdf "$FIXTURE_DIR"/calc.out \
   "$FIXTURE_DIR"/Sn3O4.HSX "$FIXTURE_DIR"/Sn3O4.EIG "$FIXTURE_DIR"/Sn3O4.XV \
   "$FIXTURE_DIR"/O.psf "$FIXTURE_DIR"/Sn.psf "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Dry run: the whole report, nothing written ---
echo -e "\n--- Testing --dry-run (report only) ---"
stb-siesta2wtb -i calc.fdf -o calc.out --dry-run --no-intro -q > log_dry.txt 2>&1
check_exit_code $? 0
check_contains "Input             : calc.fdf (SystemLabel Sn3O4)" log_dry.txt
check_contains "SaveHS in input   : yes" log_dry.txt
check_contains "Atoms             : 14 (Sn6O8)" log_dry.txt
check_contains "Sn      | 50 | 6     | 13       | 5s(2z) 5p(2z) 5d(1z)" log_dry.txt
check_contains "Orbitals          : 182 in the cell" log_dry.txt
check_contains "Spin{unpolarized} -> WanTiBEXOS type NP, 182 x 182 matrices" log_dry.txt
check_contains "nsc = 5 x 5 x 3 = 75 cells" log_dry.txt
check_contains "Fermi level       : -3.203856 eV (from calc.out)" log_dry.txt
check_contains "VBM -4.4747, CBM -2.3722 eV -> gap 2.1025 eV (indirect" log_dry.txt
check_contains "NBANDSV can be at most 36, NBANDSC at most 146" log_dry.txt
check_contains "Would write       : basis_set-NP.basis (182 orbital lines)" log_dry.txt
check_contains "Skipped           : --dry-run" log_dry.txt
echo "Testing: the energy zero of H (HSX version 1 stores E_F) is reported and warned about"
check_contains "Energy zero of H  : the run's Fermi level" log_dry.txt
check_contains "a reader that subtracts the header's E_F from H moves every band by +3.2039 eV" log_dry.txt
check_absent tb-NP.ham
check_absent basis_set-NP.basis
echo "Testing: the report is always written, in the current folder, with its references file"
check_success stb_siesta2wtb_report.txt
check_contains "STB-SIESTA2WTB REPORT" stb_siesta2wtb_report.txt
check_success stb_siesta2wtb_references.bib


# --- 3. -f 0: the header matches H, no energy-zero warning ---
echo -e "\n--- Testing -f 0 (header E_F matches H) ---"
stb-siesta2wtb -i calc.fdf -f 0 --dry-run --no-intro -q > log_f0.txt 2>&1
check_exit_code $? 0
check_contains "Fermi level       : 0.000000 eV (given with -f)" log_f0.txt
check_contains "header E_F = 0.0 matches H" log_f0.txt
check_not_contains "a reader that subtracts the header's E_F" log_f0.txt
echo "Testing: a second run never overwrites the first report (numbered _1)"
check_success stb_siesta2wtb_report_1.txt


# --- 4. Full write + read-back check ---
echo -e "\n--- Testing the full write and the read-back check (~1.5 min) ---"
stb-siesta2wtb -i calc.fdf -o calc.out --no-intro -q > log_write.txt 2>&1
check_exit_code $? 0
check_success tb-NP.ham
check_success basis_set-NP.basis
check_contains "Header            : type NP, scissor 0.0, E_F -3.203856 eV" log_write.txt
check_contains "File read         : 2484300 lines" log_write.txt
check_contains "Verdict           : ok -- within 5 eV of E_F the file reproduces the SIESTA bands" log_write.txt
check_contains "Written           : tb-NP.ham, basis_set-NP.basis" log_write.txt
echo "Testing: the .ham header and the .basis file"
check_contains "-3.203856" tb-NP.ham
basis_lines=$(grep -c . basis_set-NP.basis)
if [ "$basis_lines" -ge 182 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} basis_set-NP.basis has $basis_lines lines (>= 182 orbitals)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} basis_set-NP.basis has $basis_lines lines (expected >= 182)"
    FAIL=$((FAIL+1))
fi
rm -f tb-NP.ham basis_set-NP.basis


# --- 5. Missing input ---
echo -e "\n--- Testing a missing input file ---"
stb-siesta2wtb -i does_not_exist.fdf --no-intro -q > log_missing.txt 2>&1
rc=$?
if [ "$rc" -ne 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} non-zero exit code ($rc)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} exit code 0 for a missing input"
    FAIL=$((FAIL+1))
fi
check_contains "Input file 'does_not_exist.fdf' not found" log_missing.txt


# --- 6. Interactive path (stb-suite, shortcut 3.4) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 3.4) ---"
# Prompts in order: input .fdf, SIESTA output, Fermi level, dry run
printf '3.4\ncalc.fdf\ncalc.out\n\ny\n\n0\n' | stb-suite > log_interactive.txt 2>&1
check_contains "Dry run           : nothing written" log_interactive.txt
check_not_contains "Traceback" log_interactive.txt
check_absent tb-NP.ham


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
