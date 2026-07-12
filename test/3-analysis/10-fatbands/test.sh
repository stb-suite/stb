#!/bin/bash

# --- Setup ---
# Smoke test for stb-fatbands (Fatbands Analyzer, item 3.10)
#
# Fixture note: Sn3O4.bands/.bands.WFSX here are a SHORT band path
# (Gamma-X-Y-Gamma, 16 k-points) re-run from the same converged Sn3O4.DM
# used by test/5-utils/4-wantibexos/'s full run, purely so the .WFSX stays
# a few MB instead of the >100MB a full ~500-k-point path produces (SIESTA
# writes full wavefunction coefficients per k, not just eigenvalues).
# Sn3O4.HSX is reused as-is (its size doesn't depend on the band path).
export MPLBACKEND=Agg
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
echo "--- Starting tester for stb-fatbands (item 3.10) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR"/Sn3O4.bands "$FIXTURE_DIR"/Sn3O4.bands.WFSX "$FIXTURE_DIR"/Sn3O4.HSX \
   "$FIXTURE_DIR"/Sn3O4.XV "$FIXTURE_DIR"/calc.fdf "$FIXTURE_DIR"/structure.fdf \
   "$FIXTURE_DIR"/Sn.ion "$FIXTURE_DIR"/Sn.ion.xml "$FIXTURE_DIR"/O.ion "$FIXTURE_DIR"/O.ion.xml \
   "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Basic run: --label auto-detect (.bands.WFSX + .HSX), --projection l (default) ---
echo -e "\n--- Testing a basic run (--label Sn3O4, default --projection l) ---"

timeout 120 stb-fatbands --label Sn3O4 --shift fermi --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "Using 'Sn3O4.HSX' for overlap-aware orbital weights" log_basic.txt
check_contains "bands / .WFSX correspondence check passed" log_basic.txt
check_success fatbands_s.dat
check_success fatbands_p.dat
check_success fatbands_d.dat
check_success fatbands.gplot

echo "Testing: --projection l category files carry (k, energy, weight) columns"
check_contains "k_position" fatbands_s.dat


# --- 3. --projection species ---
echo -e "\n--- Testing --projection species ---"

timeout 120 stb-fatbands --label Sn3O4 --shift fermi --projection species --no-intro > log_species.txt 2>&1
check_exit_code $? 0
check_success fatbands_Sn.dat
check_success fatbands_O.dat


# --- 4. --projection ml ---
echo -e "\n--- Testing --projection ml ---"

timeout 120 stb-fatbands --label Sn3O4 --shift fermi --projection ml --no-intro > log_ml.txt 2>&1
check_exit_code $? 0
check_success fatbands_px.dat
check_success fatbands_dxy.dat


# --- 5. --projection atom + --category filtering ---
echo -e "\n--- Testing --projection atom with --category filtering ---"

rm -f fatbands_0.dat fatbands_1.dat
timeout 120 stb-fatbands --label Sn3O4 --shift fermi --projection atom --category 0 1 --no-intro > log_atom.txt 2>&1
check_exit_code $? 0
check_success fatbands_0.dat
check_success fatbands_1.dat
echo "Testing: an atom index NOT in --category was not written"
if [ -e fatbands_2.dat ]; then
    echo -e " ... ${RED}FAIL${NC} (fatbands_2.dat should not have been written)"
    FAIL=$((FAIL+1))
else
    echo -e " ... ${GREEN}OK${NC} (fatbands_2.dat correctly absent)"
    PASS=$((PASS+1))
fi
rm -f fatbands_2.dat fatbands_3.dat fatbands_4.dat fatbands_5.dat fatbands_6.dat \
      fatbands_7.dat fatbands_8.dat fatbands_9.dat fatbands_10.dat fatbands_11.dat \
      fatbands_12.dat fatbands_13.dat


# --- 6. Explicit --file/--wfsx/--hsx-file path (no --label) ---
echo -e "\n--- Testing explicit --file/--wfsx/--hsx-file (no --label) ---"

timeout 120 stb-fatbands --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX \
    --shift fermi --projection species -o explicit_out --no-intro > log_explicit.txt 2>&1
check_exit_code $? 0
check_success explicit_out/fatbands_Sn.dat


# --- 7. No .HSX given via --label: no <label>.fdf exists either (the real
#        fdf here is named calc.fdf, not Sn3O4.fdf), so this must fail
#        cleanly rather than silently produce a wrong/empty plot.
echo -e "\n--- Testing --label with no .HSX and no <label>.fdf (clean failure) ---"

mv Sn3O4.HSX Sn3O4.HSX.bak
timeout 120 stb-fatbands --label Sn3O4 --shift fermi --projection l -o approx_out --no-intro > log_approx.txt 2>&1
check_exit_code $? 2
check_contains "No geometry/Hamiltonian could be resolved" log_approx.txt

# --- 7b. No .HSX, but an explicit --geometry-file (.fdf + its .ion/.ion.xml
#         files) is given: the approximate orthogonal-basis fallback should
#         actually work and print its accuracy warning.
echo -e "\n--- Testing the --geometry-file fallback (no .HSX, approximate weights) ---"

timeout 120 stb-fatbands --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --geometry-file calc.fdf \
    --shift fermi --projection l -o approx_out2 --no-intro > log_approx2.txt 2>&1
check_exit_code $? 0
check_contains "No .HSX found" log_approx2.txt
check_success approx_out2/fatbands_s.dat
mv Sn3O4.HSX.bak Sn3O4.HSX


# --- 8. k-point count mismatch is a hard error, not a silent wrong plot ---
echo -e "\n--- Testing the k-point/.WFSX count-mismatch guard ---"

WANTIBEXOS_BANDS="$FIXTURE_DIR/../../5-utils/4-wantibexos/Sn3O4.bands"
if [ -f "$WANTIBEXOS_BANDS" ]; then
    timeout 120 stb-fatbands --file "$WANTIBEXOS_BANDS" --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX \
        --shift fermi --no-intro > log_mismatch.txt 2>&1
    check_exit_code $? 2
    check_contains "k-point count mismatch" log_mismatch.txt
else
    echo -e " ... ${YELLOW}SKIPPED${NC} (test/5-utils/4-wantibexos/Sn3O4.bands not found)"
fi


# --- 9. Missing .WFSX aborts with a clear error ---
echo -e "\n--- Testing a missing .WFSX ---"

timeout 30 stb-fatbands --label does_not_exist --shift fermi --no-intro > log_missing_wfsx.txt 2>&1
check_exit_code $? 2
check_contains "not found" log_missing_wfsx.txt


# --- 10. Robustness / --help / --version ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: --label omitted (required)"
stb-fatbands --shift fermi --no-intro > log_missing_label.txt 2>&1
check_exit_code $? 2

echo "Testing: --shift omitted (required)"
stb-fatbands --label Sn3O4 --no-intro > log_missing_shift.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-fatbands --version > log_version.txt 2>&1
check_contains "stb-fatbands" log_version.txt

echo "Testing: --help documents --wfsx, --hsx-file, --geometry-file, --projection, --category, -o"
stb-fatbands --help > log_help.txt 2>&1
check_contains "wfsx" log_help.txt
check_contains "hsx-file" log_help.txt
check_contains "geometry-file" log_help.txt
check_contains "projection" log_help.txt
check_contains "category" log_help.txt
check_contains "output-dir" log_help.txt


# --- 11. Interactive path (stb-suite, shortcut 3.10) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 3.10) ---"

echo "Testing: navigate 3.10 -> label Sn3O4 -> shift fermi -> projection l -> default output"
rm -f fatbands_s.dat fatbands_p.dat fatbands_d.dat
printf '3.10\nSn3O4\n3\n1\n\n\n0\n' | timeout 120 stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success fatbands_s.dat


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
