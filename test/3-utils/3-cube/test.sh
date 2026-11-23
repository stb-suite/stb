#!/bin/bash

# --- Setup ---
# Smoke test for stb-cube (Grid to Cube Converter, item 3.3).
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
    if grep -q "$1" "$2" 2>/dev/null; then
        echo -e "   -> ${RED}Failed:${NC} '$1' unexpectedly found in '$2'"
        FAIL=$((FAIL+1))
    else
        echo -e "   -> ${GREEN}Verified:${NC} '$1' NOT found in '$2'"
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
# Sn3O4.RHO/.XV: real, non-spin-polarized SIESTA output (shared with other
# 3-utils fixtures). o2.RHO/.XV/.HSX: a real spin-polarized SIESTA run
# generated specifically for the --spin option and .HSX-based auto-detection
# (isolated O2 molecule, triplet ground state -- see o2.fdf here for the
# exact input; SIESTA itself reports spin moment |S| = 2.0, which is what
# the --spin diff checks below verify against, and o2.HSX correctly reads
# back as Spin{polarized} via sisl).
echo "--- Starting tester for STB-Cube (item 3.3) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/Sn3O4.RHO" "$TEST_DIR/"
cp "$FIXTURE_DIR/Sn3O4.XV" "$TEST_DIR/"
cp "$FIXTURE_DIR/o2.RHO" "$TEST_DIR/"
cp "$FIXTURE_DIR/o2.XV" "$TEST_DIR/"
cp "$FIXTURE_DIR/o2.HSX" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. RHO -> cube, default output name, unit-corrected density + XV mapping ---
echo -e "\n--- Testing RHO -> cube (default output name) ---"
stb-cube --no-intro --label Sn3O4 --type RHO > log_rho.txt 2>&1
check_success Sn3O4_RHO.cube
check_contains "Integral" log_rho.txt
echo "Testing: the cube is written by the tool (Gaussian-style header ASE parses, Bohr, e/Bohr^3)"
check_contains "lengths in Bohr, values in e/Bohr^3" Sn3O4_RHO.cube
check_contains "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z" Sn3O4_RHO.cube
check_contains "Integral (cube)   : 71.99" log_rho.txt
check_not_contains "lattice mismatch" log_rho.txt
check_success stb_cube_report.txt
check_success stb_cube_references.bib


# --- 3. Custom -o output filename ---
echo -e "\n--- Testing custom -o output filename ---"
rm -f custom_name.cube
stb-cube --no-intro --label Sn3O4 --type RHO -o custom_name.cube > log_custom.txt 2>&1
check_success custom_name.cube


# --- 4. Missing input file ---
echo -e "\n--- Testing missing input file (label with no matching grid file) ---"
stb-cube --no-intro --label DoesNotExist --type VT > log_missing.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing.txt


# --- 5. Missing .XV: still converts, warns, no mismatch check possible ---
echo -e "\n--- Testing missing .XV (atoms will be absent, no crash) ---"
mv Sn3O4.XV Sn3O4.XV.bak
rm -f Sn3O4_RHO.cube
stb-cube --no-intro --label Sn3O4 --type RHO > log_noxv.txt 2>&1
check_contains "no Sn3O4.XV: the cube has no atoms" log_noxv.txt
check_success Sn3O4_RHO.cube
mv Sn3O4.XV.bak Sn3O4.XV


# --- 6. Mismatched .XV: warns instead of silently mismapping atoms ---
echo -e "\n--- Testing mismatched .XV (different lattice) ---"
cat > Sn3O4.XV << 'EOF'
1.0 0.0 0.0    0.0 0.0 0.0
0.0 1.0 0.0    0.0 0.0 0.0
0.0 0.0 1.0    0.0 0.0 0.0
1
1 14   0.0 0.0 0.0   0.0 0.0 0.0
EOF
stb-cube --no-intro --label Sn3O4 --type RHO > log_mismatch.txt 2>&1
check_contains "lattice mismatch between the RHO grid and 'Sn3O4.XV'" log_mismatch.txt
cp "$FIXTURE_DIR/Sn3O4.XV" .


# --- 7. --spin channels, verified against a real spin-polarized run ---
# SIESTA itself reports "spin moment: {S} ... |S| = 2.00000" for this O2
# calculation (see o2.out from the original run, or re-run o2.fdf) -- total
# must equal up+down, and diff (up-down) must equal that same 2.0 exactly.
# o2.HSX makes the auto-detection ("Spin configuration (from o2.HSX): ...")
# kick in too, verified against a real sisl.read_hamiltonian() read.
echo -e "\n--- Testing --spin (total/up/down/diff) against a real spin-polarized O2 run ---"
stb-cube --no-intro --label o2 --type RHO --spin total > log_spin_total.txt 2>&1
check_contains "Spin configuration (from o2.HSX): polarized" log_spin_total.txt
check_contains "12.0000 e" log_spin_total.txt
check_success o2_RHO.cube

stb-cube --no-intro --label o2 --type RHO --spin up > log_spin_up.txt 2>&1
check_contains "up                | 7.0000       |        | <- written" log_spin_up.txt
check_contains "Integral          : 6.99999" log_spin_up.txt
check_success o2_RHO_up.cube

stb-cube --no-intro --label o2 --type RHO --spin down > log_spin_down.txt 2>&1
check_contains "down              | 5.0000       |        | <- written" log_spin_down.txt
check_success o2_RHO_down.cube

stb-cube --no-intro --label o2 --type RHO --spin diff > log_spin_diff.txt 2>&1
check_contains "Converting        : the 'diff' channel" log_spin_diff.txt
check_contains "Integral          : 2.000000 e" log_spin_diff.txt
check_success o2_RHO_diff.cube

echo "Testing: --spin on a non-spin-polarized file falls back to total, with a warning"
stb-cube --no-intro --label Sn3O4 --type RHO --spin up -o fallback_check.cube > log_spin_fallback.txt 2>&1
check_contains "not spin-polarized (1 component) -- ignoring --spin up" log_spin_fallback.txt
check_contains "Integral          : 71.99" log_spin_fallback.txt
check_success fallback_check.cube

echo "Testing: .HSX/grid component-count mismatch warning (o2.HSX says polarized, paired with Sn3O4's 1-component RHO)"
cp o2.HSX Sn3O4.HSX
stb-cube --no-intro --label Sn3O4 --type RHO > log_hsx_mismatch.txt 2>&1
check_contains "Sn3O4.HSX reports polarized (2 component(s) expected) but 'Sn3O4.RHO' has 1" log_hsx_mismatch.txt
rm -f Sn3O4.HSX


# --- 8. Reports: always written, never overwritten; --save-report is a hidden, ignored option ---
echo -e "\n--- Testing the report files (and the old --save-report flag) ---"
rm -f stb_cube_report*.txt stb_cube_references*.bib
stb-cube --no-intro --label o2 --type RHO --save-report > log_savereport.txt 2>&1
check_exit_code $? 0
check_success stb_cube_report.txt
check_contains "STB-CUBE REPORT" stb_cube_report.txt
stb-cube --no-intro --label o2 --type RHO > log_savereport2.txt 2>&1
check_success stb_cube_report_1.txt
check_success stb_cube_references_1.bib


# --- 9. Interactive path (stb-suite, shortcut 3.3) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 3.3) ---"
rm -f Sn3O4_RHO.cube interactive_out.cube
# Prompts in order: type choice, label, spin choice, output filename, dry run
printf '3.3\n1\nSn3O4\n1\ninteractive_out.cube\nn\n\n0\n' | stb-suite > log_interactive.txt 2>&1
check_success interactive_out.cube
check_contains "Conversion complete" log_interactive.txt
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
