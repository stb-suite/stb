#!/bin/bash

# --- Setup ---
# Smoke test for stb-xrd (Powder XRD Simulator, item 3.9)
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$FIXTURE_DIR/test_files"

# Non-interactive matplotlib backend: --view calls plt.show(), which would
# otherwise block forever waiting for a GUI event loop that doesn't exist here.
# Under Agg, plt.show() just warns and returns immediately -- verified directly
# before writing this test (a plain plt.show() under MPLBACKEND=Agg returned in
# ~0.006s) -- so, unlike some other tools' test.sh in this suite, --view itself
# is exercised directly below rather than only documented via --help.
export MPLBACKEND=Agg

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

# Checks that file $1 does NOT exist (or is empty)
check_absent() {
    if [ -s "$1" ]; then
        echo -e " ... ${RED}FAIL${NC} (file '$1' exists but shouldn't have been written)"
        FAIL=$((FAIL+1))
    else
        echo -e " ... ${GREEN}OK${NC} (file '$1' correctly not written)"
        PASS=$((PASS+1))
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
echo "--- Starting tester for STB-XRD (item 3.9) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/siesta.STRUCT_OUT" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null

# Build a small NaCl structure with stb-crystalcast to use as a known fixture
stb-crystalcast --group 225 --species Na Cl --num-ions 4 4 \
    --lattice 5.64 5.64 5.64 90 90 90 -o nacl.fdf --no-intro > /dev/null 2>&1


# --- 2. Basic run: numbered report, no unconditional data file ---
echo -e "\n--- Testing a basic run (--file nacl.fdf --format fdf) ---"

rm -f xrd_pattern.dat xrd_pattern.gplot references.bib
stb-xrd --file nacl.fdf --format fdf --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "\[0\] RUN METADATA" log_basic.txt
check_contains "\[1\] STRUCTURE" log_basic.txt
check_contains "\[2\] DIFFRACTION PATTERN" log_basic.txt
check_contains "\[4\] OUTPUT DATA & PLOTS" log_basic.txt
check_contains "\[5\] REFERENCES" log_basic.txt
check_contains "\[6\] SUMMARY & FILES" log_basic.txt
check_contains "Formula.*NaCl" log_basic.txt
check_contains "Wavelength.*CuKa" log_basic.txt
check_contains "Space Group" log_basic.txt
check_contains "Crystal System" log_basic.txt
check_contains "Peaks found in range" log_basic.txt
check_contains "Strongest peak" log_basic.txt

echo "Testing: without --save-gnuplot, no .dat/.gplot is written (used to be unconditional)"
check_absent xrd_pattern.dat
check_absent xrd_pattern.gplot

echo "Testing: references.bib is always written (SIESTA + pyxtal)"
check_success references.bib
check_contains "Soler2002" references.bib
check_contains "Garcia2020" references.bib
check_contains "Fredericks2021pyxtal" references.bib

echo "Testing: no text report without --save-report"
check_absent stb_xrd_report.txt


# --- 3. --save-report / --save-gnuplot ---
echo -e "\n--- Testing --save-report / --save-gnuplot ---"

rm -f xrd_pattern.dat xrd_pattern.gplot stb_xrd_report.txt
stb-xrd --file nacl.fdf --format fdf --save-report --save-gnuplot --no-intro > log_saved.txt 2>&1
check_exit_code $? 0
check_success stb_xrd_report.txt
check_contains "\[0\] RUN METADATA" stb_xrd_report.txt
check_success xrd_pattern.dat
check_contains "wavelength/thetas" xrd_pattern.dat
check_success xrd_pattern.gplot
check_contains "using 1:6 with impulses" xrd_pattern.gplot
check_contains "\[OK\] Data written to" log_saved.txt
check_contains "\[OK\] Gnuplot script written to" log_saved.txt


# --- 4. -o/--output-dir ---
echo -e "\n--- Testing -o/--output-dir ---"

rm -rf out_subdir
stb-xrd --file nacl.fdf --format fdf --output-dir out_subdir --save-report --save-gnuplot \
    --no-intro > log_outdir.txt 2>&1
check_exit_code $? 0
check_success out_subdir/references.bib
check_success out_subdir/stb_xrd_report.txt
check_success out_subdir/xrd_pattern.dat
check_success out_subdir/xrd_pattern.gplot


# --- 5. --wavelength (named source and numeric) ---
echo -e "\n--- Testing --wavelength ---"

stb-xrd --file nacl.fdf --format fdf --wavelength MoKa --no-intro > log_wl_named.txt 2>&1
check_exit_code $? 0
check_contains "Wavelength.*MoKa" log_wl_named.txt

stb-xrd --file nacl.fdf --format fdf --wavelength 1.5406 --no-intro > log_wl_numeric.txt 2>&1
check_exit_code $? 0
check_contains "1.54060 Ang" log_wl_numeric.txt

echo "Testing: unknown --wavelength name"
stb-xrd --file nacl.fdf --format fdf --wavelength NotASource --no-intro > log_wl_bad.txt 2>&1
check_exit_code $? 1
check_contains "not a known source name" log_wl_bad.txt


# --- 6. --two-theta-range ---
echo -e "\n--- Testing --two-theta-range ---"

stb-xrd --file nacl.fdf --format fdf --two-theta-range 20 60 --no-intro > log_range.txt 2>&1
check_exit_code $? 0
check_contains "20.0 - 60.0 deg" log_range.txt

echo "Testing: --two-theta-range MIN >= MAX is rejected"
stb-xrd --file nacl.fdf --format fdf --two-theta-range 60 20 --no-intro > log_range_bad.txt 2>&1
check_exit_code $? 1
check_contains "MIN must be less than MAX" log_range_bad.txt

echo "Testing: an empty --two-theta-range window (no peaks) is reported cleanly"
stb-xrd --file nacl.fdf --format fdf --two-theta-range 0.1 0.5 --no-intro > log_range_empty.txt 2>&1
check_exit_code $? 1
check_contains "no peaks found" log_range_empty.txt


# --- 7. --top ---
echo -e "\n--- Testing --top ---"

stb-xrd --file nacl.fdf --format fdf --top 3 --no-intro > log_top.txt 2>&1
check_exit_code $? 0
check_contains "Showing top 3 of" log_top.txt


# --- 8. struct_out format ---
echo -e "\n--- Testing --format struct_out ---"

stb-xrd --file siesta.STRUCT_OUT --format struct_out --no-intro > log_struct_out.txt 2>&1
check_exit_code $? 0
check_contains "\[1\] STRUCTURE" log_struct_out.txt


# --- 9. --view (renamed from --plot; MPLBACKEND=Agg above makes this safe/non-blocking) ---
echo -e "\n--- Testing --view (headless via MPLBACKEND=Agg, only checking exit code) ---"

timeout 20 stb-xrd --file nacl.fdf --format fdf --view --no-intro > log_view.txt 2>&1
check_exit_code $? 0


# --- 10. --compare-to ---
echo -e "\n--- Testing --compare-to (pyxtal's Similarity() is genuinely slow, ~20-30s -- generous timeouts below) ---"

rm -f xrd_pattern.dat
stb-xrd --file nacl.fdf --format fdf --save-gnuplot --no-intro > /dev/null 2>&1

# Build a synthetic "experimental" file from this structure's own simulated
# peaks (2theta + intensity columns only) -- not a perfect broadened profile,
# but enough to exercise the full --compare-to code path end to end.
awk 'NR>1 {print $1, $6}' xrd_pattern.dat > experimental.dat

timeout 60 stb-xrd --file nacl.fdf --format fdf --compare-to experimental.dat \
    --no-intro > log_compare.txt 2>&1
check_exit_code $? 0
check_contains "\[3\] EXPERIMENTAL COMPARISON" log_compare.txt
check_contains "Experimental file.*experimental.dat" log_compare.txt
check_contains "Similarity score" log_compare.txt

echo "Testing: --compare-to combined with --view (overlay branch, headless)"
timeout 60 stb-xrd --file nacl.fdf --format fdf --compare-to experimental.dat --view \
    --no-intro > log_compare_view.txt 2>&1
check_exit_code $? 0

echo "Testing: --compare-to with a nonexistent file"
stb-xrd --file nacl.fdf --format fdf --compare-to does_not_exist.dat --no-intro > log_compare_missing.txt 2>&1
check_exit_code $? 1
check_contains "No such file or directory" log_compare_missing.txt

echo "Testing: --compare-to with a malformed file"
printf "not a number here\n1.0 2.0\n" > bad_experimental.dat
stb-xrd --file nacl.fdf --format fdf --compare-to bad_experimental.dat --no-intro > log_compare_bad.txt 2>&1
check_exit_code $? 1
check_contains "expected two numbers" log_compare_bad.txt

echo "Testing: --compare-to with too few points (pyxtal's cubic interpolation needs >= 4)"
printf "10.0 1.0\n20.0 2.0\n30.0 1.5\n" > sparse_experimental.dat
stb-xrd --file nacl.fdf --format fdf --compare-to sparse_experimental.dat --no-intro > log_compare_sparse.txt 2>&1
check_exit_code $? 1
check_contains "at least 4 are needed" log_compare_sparse.txt


# --- 11. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: nonexistent structure file"
stb-xrd --file does_not_exist.fdf --format fdf --no-intro > log_missing_file.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_file.txt

echo "Testing: --file omitted (required)"
stb-xrd --format fdf --no-intro > log_missing_arg_file.txt 2>&1
check_exit_code $? 2

echo "Testing: --format omitted (required)"
stb-xrd --file nacl.fdf --no-intro > log_missing_arg_format.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-xrd --version > log_version.txt 2>&1
check_contains "stb-xrd" log_version.txt

echo "Testing: --help documents --file, --format, --wavelength, --two-theta-range, --top, --compare-to, --output-dir, --save-report, --save-gnuplot, --view"
stb-xrd --help > log_help.txt 2>&1
check_contains "file" log_help.txt
check_contains "format" log_help.txt
check_contains "wavelength" log_help.txt
check_contains "two-theta-range" log_help.txt
check_contains "top" log_help.txt
check_contains "compare-to" log_help.txt
check_contains "output-dir" log_help.txt
check_contains "save-report" log_help.txt
check_contains "save-gnuplot" log_help.txt
check_contains "\-\-view" log_help.txt


# --- 12. Interactive path (stb-suite, shortcut 3.9) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 3.9) ---"

echo "Testing: navigate 3.9 -> nacl.fdf -> format fdf (1) -> default wavelength -> default range -> no top -> no compare-to -> no save-report -> no save-gnuplot -> no view -> quit"
rm -rf menu_test && mkdir menu_test && cp nacl.fdf menu_test/
pushd menu_test > /dev/null
printf '3.9\nnacl.fdf\n1\n\n\n\n\n\nn\nn\nn\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_contains "Formula.*NaCl" log_menu.txt
check_success references.bib
check_absent xrd_pattern.dat
popd > /dev/null


popd > /dev/null

# --- 13. Summary ---
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
