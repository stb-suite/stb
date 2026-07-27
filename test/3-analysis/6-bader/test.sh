#!/bin/bash

# --- Setup ---
# Smoke test for stb-bader (Bader Charge Analysis, item 3.6)
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

# Checks that file $1 does NOT exist
check_absent() {
    if [ ! -e "$1" ]; then
        echo -e " ... ${GREEN}OK${NC} (file '$1' absent, as expected)"
        PASS=$((PASS+1))
    else
        echo -e " ... ${RED}FAIL${NC} (file '$1' should not exist)"
        FAIL=$((FAIL+1))
    fi
}

# Checks that file $2 contains (grep -q) pattern $1
check_contains() {
    if grep -q -- "$1" "$2" 2>/dev/null; then
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
echo "--- Starting tester for STB-Bader (item 3.6) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/Sn3O4.out" "$FIXTURE_DIR/Sn3O4.RHO" "$FIXTURE_DIR/Sn3O4.XV" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Basic run ---
echo -e "\n--- Testing a basic run (--label Sn3O4) ---"

timeout 120 stb-bader --label Sn3O4 --no-intro --output-dir out1 --save-report > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "Z_vals detected" log_basic.txt
check_success out1/stb_bader_report.txt
check_contains "===== STB-BADER REPORT =====" out1/stb_bader_report.txt
check_contains "\[0\] RUN METADATA" out1/stb_bader_report.txt
check_contains "\[1\] VALENCE (Z_val) SETUP" out1/stb_bader_report.txt
check_contains "\[2\] PYBADER CONFIGURATION" out1/stb_bader_report.txt
check_contains "\[3\] PER-ATOM BADER POPULATIONS" out1/stb_bader_report.txt
check_contains "\[4\] PER-SPECIES SUMMARY" out1/stb_bader_report.txt
check_contains "\[5\] DIAGNOSTICS & WARNINGS" out1/stb_bader_report.txt
check_contains "\[6\] REFERENCES" out1/stb_bader_report.txt
check_contains "\[7\] SUMMARY & FILES" out1/stb_bader_report.txt
check_contains "Refine method    | neargrid" out1/stb_bader_report.txt
check_contains "Total Integrated:" out1/stb_bader_report.txt
check_success out1/Sn3O4.cube
check_success out1/references.bib
check_contains "@article{Soler2002," out1/references.bib

echo "Testing: essentially-zero-population atoms are flagged (this fixture genuinely has some)"
check_contains "WARN\] Atom(s).*essentially zero Bader population" log_basic.txt
check_contains "WARNING\] Atom(s).*essentially zero Bader population" out1/stb_bader_report.txt

echo "Testing: Volume/SurfDist columns and per-species summary are reported"
check_contains "Vol(Å³)" log_basic.txt
check_contains "SurfD(Å)" log_basic.txt
check_contains "Vol(Å³)" out1/stb_bader_report.txt
check_contains "Elem | N | Mean(e-) | Std(e-)" out1/stb_bader_report.txt


# --- 3. --speed fast uses PyBader's own on-grid profile ---
echo -e "\n--- Testing --speed fast ---"

timeout 120 stb-bader --label Sn3O4 --no-intro --speed fast --output-dir out2 --save-report > log_fast.txt 2>&1
check_exit_code $? 0
check_contains "CAUTION.*--speed fast trades" log_fast.txt
check_contains "Method           | ongrid" out2/stb_bader_report.txt


# --- 4. --threads ---
echo -e "\n--- Testing --threads 2 ---"

timeout 120 stb-bader --label Sn3O4 --no-intro --threads 2 --output-dir out3 --save-report > log_threads.txt 2>&1
check_exit_code $? 0
check_contains "Threads          | 2" out3/stb_bader_report.txt

echo "Testing: --threads 0 is rejected"
stb-bader --label Sn3O4 --no-intro --threads 0 > log_threads0.txt 2>&1
check_exit_code $? 2
check_contains "must be at least 1" log_threads0.txt


# --- 5. --vacuum-tol ---
echo -e "\n--- Testing --vacuum-tol ---"

timeout 120 stb-bader --label Sn3O4 --no-intro --vacuum-tol 0.001 --output-dir out4 --save-report > log_vac.txt 2>&1
check_exit_code $? 0
check_contains "Vacuum tolerance | 0.001" out4/stb_bader_report.txt


# --- 6. --no-cube removes the intermediate .cube file(s) ---
echo -e "\n--- Testing --no-cube ---"

timeout 120 stb-bader --label Sn3O4 --no-intro --no-cube --output-dir out5 --save-report > log_nocube.txt 2>&1
check_exit_code $? 0
check_absent out5/Sn3O4.cube
check_contains "Cube file(s)   : deleted (--no-cube)" out5/stb_bader_report.txt


# --- 6b. --export-volumes writes one .cube per atom inside --output-dir,
#     and never leaks PyBader's own bader.p ---
echo -e "\n--- Testing --export-volumes ---"

timeout 120 stb-bader --label Sn3O4 --no-intro --export-volumes --output-dir out6 --save-report > log_export.txt 2>&1
check_exit_code $? 0
check_contains "Exported 14 per-atom Bader volume" log_export.txt
check_success out6/Bader-atoms-0.cube
check_success out6/Bader-atoms-13.cube
check_absent out6/bader.p
check_absent bader.p
check_absent Bader-atoms-0.cube


# --- 7. --ref pointing at a missing file falls back to hardcoded defaults ---
echo -e "\n--- Testing --ref with a missing file (falls back, doesn't abort) ---"

timeout 120 stb-bader --label Sn3O4 --no-intro --ref does_not_exist.out --output-dir out7 --save-report > log_fallback.txt 2>&1
check_exit_code $? 0
check_contains "USING DEFAULT VALUES" log_fallback.txt
check_contains "Hardcoded Defaults (FALLBACK)" out7/stb_bader_report.txt


# --- 8. Missing .RHO aborts with a non-zero exit code ---
echo -e "\n--- Testing a missing .RHO file ---"

timeout 30 stb-bader --label does_not_exist --no-intro > log_missing_rho.txt 2>&1
check_exit_code $? 1
check_contains "Grid file 'does_not_exist.RHO' not found" log_missing_rho.txt


# --- 9. --output-dir / --save-report: a plain run writes references.bib
#     only; the old <label>_BADER.txt is no longer written unconditionally
#     (and a stale one from an older version is cleaned up automatically) ---
echo -e "\n--- Testing default output (no --save-report) and stale <label>_BADER.txt cleanup ---"

rm -rf out8
echo "stale content from an older version" > Sn3O4_BADER.txt
timeout 120 stb-bader --label Sn3O4 --no-intro --output-dir out8 > log_default.txt 2>&1
check_exit_code $? 0
check_success out8/references.bib
if [ -e out8/stb_bader_report.txt ] || [ -e Sn3O4_BADER.txt ]; then
    echo -e "   -> ${RED}Failed:${NC} a report file exists without --save-report, or the stale file wasn't cleaned up"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no stb_bader_report.txt without --save-report, and the stale Sn3O4_BADER.txt is gone"
    PASS=$((PASS+1))
fi


# --- 10. Robustness / --help / --version ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: --label omitted (required)"
stb-bader --no-intro > log_missing_label.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-bader --version > log_version.txt 2>&1
check_contains "stb-bader" log_version.txt

echo "Testing: --help documents --label, --output-dir, --save-report, --speed, --threads, --vacuum-tol, --no-cube, --export-volumes"
stb-bader --help > log_help.txt 2>&1
check_contains "label" log_help.txt
check_contains "output-dir" log_help.txt
check_contains "save-report" log_help.txt
check_contains "speed" log_help.txt
check_contains "threads" log_help.txt
check_contains "vacuum-tol" log_help.txt
check_contains "no-cube" log_help.txt
check_contains "export-volumes" log_help.txt


# --- 11. Interactive path (stb-suite, shortcut 3.6) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 3.6) ---"

echo "Testing: navigate 3.6 -> label Sn3O4 -> default output dir -> default ref -> speed 1 -> no vacuum-tol -> no threads -> keep cube -> no export -> save-report=y -> quit"
printf '3.6\nSn3O4\n\n\n1\n\n\ny\nn\ny\n\n0\n' | timeout 120 stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success stb_bader_report.txt
check_success references.bib


popd > /dev/null

# --- 12. Summary ---
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
