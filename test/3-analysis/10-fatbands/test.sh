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


# --- 2. Basic run: --label auto-detect (.bands.WFSX + .HSX), default --projection (species_l) ---
echo -e "\n--- Testing a basic run (--label Sn3O4, default --projection species_l) ---"

timeout 120 stb-fatbands --label Sn3O4 --shift fermi --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "Using 'Sn3O4.HSX' for overlap-aware orbital weights" log_basic.txt
check_contains "bands / .WFSX correspondence check passed" log_basic.txt
check_contains "\[0\] RUN METADATA" log_basic.txt
check_contains "\[1\] INPUT DATA" log_basic.txt
check_contains "\[2\] BAND GAP ANALYSIS" log_basic.txt
check_contains "\[3\] ORBITAL PROJECTION" log_basic.txt
check_contains "\[4\] WRITING OUTPUT FILES" log_basic.txt
check_contains "\[5\] REFERENCES" log_basic.txt
check_contains "\[6\] SUMMARY & FILES" log_basic.txt
check_contains "VBM" log_basic.txt
echo "Testing: default projection (no --projection given) is species_l, not plain l"
check_contains "Projection      : species_l" log_basic.txt
check_contains "Categories found : 6 (O-s, O-p, O-d, Sn-s, Sn-p, Sn-d)" log_basic.txt

echo "Testing: without --save-gnuplot, no .dat/.gplot is written (used to be unconditional)"
if [ -e fatbands_Sn-s.dat ] || [ -e fatbands.gplot ]; then
    echo -e " ... ${RED}FAIL${NC} (fatbands_Sn-s.dat/fatbands.gplot should not have been written)"
    FAIL=$((FAIL+1))
else
    echo -e " ... ${GREEN}OK${NC} (fatbands_Sn-s.dat/fatbands.gplot correctly absent)"
    PASS=$((PASS+1))
fi
check_contains "Not written (off by default" log_basic.txt

echo "Testing: no text report without --save-report"
if [ -e stb_fatbands_report.txt ]; then
    echo -e " ... ${RED}FAIL${NC} (stb_fatbands_report.txt should not have been written)"
    FAIL=$((FAIL+1))
else
    echo -e " ... ${GREEN}OK${NC} (stb_fatbands_report.txt correctly absent)"
    PASS=$((PASS+1))
fi

echo "Testing: references.bib is always written (SIESTA)"
check_success references.bib


# --- 2b. --save-report / --save-gnuplot (explicit --projection l, the pre-2.1.0 default) ---
echo -e "\n--- Testing --save-report / --save-gnuplot (--projection l) ---"

timeout 120 stb-fatbands --label Sn3O4 --shift fermi --projection l \
    --save-report --save-gnuplot --no-intro > log_saved.txt 2>&1
check_exit_code $? 0
check_success stb_fatbands_report.txt
check_contains "\[0\] RUN METADATA" stb_fatbands_report.txt
check_success fatbands_s.dat
check_success fatbands_p.dat
check_success fatbands_d.dat
check_success fatbands.gplot
check_contains "\[OK\] Data written" log_saved.txt
check_contains "\[OK\] Gnuplot script written to" log_saved.txt

echo "Testing: --projection l category files carry (k, energy, weight) columns"
check_contains "k_position" fatbands_s.dat


# --- 3. --projection species ---
echo -e "\n--- Testing --projection species ---"

rm -f fatbands_Sn.dat fatbands_O.dat
timeout 120 stb-fatbands --label Sn3O4 --shift fermi --projection species --save-gnuplot --no-intro > log_species.txt 2>&1
check_exit_code $? 0
check_success fatbands_Sn.dat
check_success fatbands_O.dat


# --- 3b. --projection species_l (species AND s/p/d/f combined) ---
echo -e "\n--- Testing --projection species_l ---"

timeout 120 stb-fatbands --label Sn3O4 --shift fermi --projection species_l --save-gnuplot --no-intro > log_species_l.txt 2>&1
check_exit_code $? 0
check_contains "Categories found : 6 (O-s, O-p, O-d, Sn-s, Sn-p, Sn-d)" log_species_l.txt
check_success fatbands_Sn-s.dat
check_success fatbands_Sn-p.dat
check_success fatbands_Sn-d.dat
check_success fatbands_O-s.dat
check_success fatbands_O-p.dat
check_success fatbands_O-d.dat

echo "Testing: --projection species_l with --category filtering"
rm -f fatbands_Sn-s.dat fatbands_Sn-p.dat fatbands_Sn-d.dat \
      fatbands_O-s.dat fatbands_O-p.dat fatbands_O-d.dat
timeout 120 stb-fatbands --label Sn3O4 --shift fermi --projection species_l --category Sn-s O-p \
    --save-gnuplot --no-intro > log_species_l_cat.txt 2>&1
check_exit_code $? 0
check_success fatbands_Sn-s.dat
check_success fatbands_O-p.dat
if [ -e fatbands_Sn-p.dat ]; then
    echo -e " ... ${RED}FAIL${NC} (fatbands_Sn-p.dat should not have been written -- not in --category)"
    FAIL=$((FAIL+1))
else
    echo -e " ... ${GREEN}OK${NC} (fatbands_Sn-p.dat correctly absent)"
    PASS=$((PASS+1))
fi


# --- 3c. Spin-polarized (nspin=2): categories must split per spin ---
# Fixture note: Ospin.bands/.bands.WFSX/.HSX are a REAL SIESTA run (not
# synthetic) -- a single O atom in a large vacuum box, Spin polarized,
# seeded via %block DM.InitSpin to converge to its physical 2 Bohr-magneton
# triplet ground state (verified in Ospin.out: "spin moment ... |S| = 2.0").
# This is the fixture that caught a real bug: the original weight loop
# merged both spin channels into one category with no spin label at all,
# silently combining two very different band sets (here, spin-up/spin-down
# CBM differ by ~29 eV) into one indistinguishable series.
echo -e "\n--- Testing --projection species_l on a real spin-polarized (nspin=2) calculation ---"

rm -rf spin_test && mkdir spin_test
cp "$FIXTURE_DIR"/spin/Ospin.fdf "$FIXTURE_DIR"/spin/Ospin.bands "$FIXTURE_DIR"/spin/Ospin.bands.WFSX \
   "$FIXTURE_DIR"/spin/Ospin.HSX "$FIXTURE_DIR"/spin/O.ion "$FIXTURE_DIR"/spin/O.ion.xml spin_test/
pushd spin_test > /dev/null

timeout 120 stb-fatbands --label Ospin --shift fermi --projection species_l \
    --save-gnuplot --no-intro > log_spin.txt 2>&1
check_exit_code $? 0
check_contains "Spin channels" log_spin.txt
check_contains "2 (polarized)" log_spin.txt
check_contains "Spin up:" log_spin.txt
check_contains "Spin down:" log_spin.txt
echo "Testing: each category is split into <category>_up/<category>_down for nspin=2"
check_contains "Spin-resolved    : yes -- each category split into <category>_up/<category>_down (6 series total)." log_spin.txt
check_success fatbands_O-s_up.dat
check_success fatbands_O-s_down.dat
check_success fatbands_O-p_up.dat
check_success fatbands_O-p_down.dat
check_success fatbands_O-d_up.dat
check_success fatbands_O-d_down.dat

echo "Testing: spin-up and spin-down data for the same category cover different energy ranges"
echo "(a real, physically-expected spin splitting -- not a bug if these numbers differ)"
UP_LINES=$(grep -vc "^#" fatbands_O-s_up.dat)
DOWN_LINES=$(grep -vc "^#" fatbands_O-s_down.dat)
if [ "$UP_LINES" -eq "$DOWN_LINES" ] && [ "$UP_LINES" -gt 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} fatbands_O-s_up.dat and fatbands_O-s_down.dat both have "
    echo -e "      $UP_LINES data lines each (13 bands x 6 k-points, one spin channel each --"
    echo -e "      not merged into one 2x-length file anymore)."
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} expected equal, nonzero line counts (got up=$UP_LINES, down=$DOWN_LINES)"
    FAIL=$((FAIL+1))
fi

echo "Testing: --view on a spin-polarized run (multi-series plot path, headless)"
timeout 120 stb-fatbands --label Ospin --shift fermi --projection species_l --view \
    --no-intro > log_spin_view.txt 2>&1
check_exit_code $? 0

popd > /dev/null


# --- 4. --projection ml ---
echo -e "\n--- Testing --projection ml ---"

timeout 120 stb-fatbands --label Sn3O4 --shift fermi --projection ml --save-gnuplot --no-intro > log_ml.txt 2>&1
check_exit_code $? 0
check_success fatbands_px.dat
check_success fatbands_dxy.dat


# --- 5. --projection atom + --category filtering ---
echo -e "\n--- Testing --projection atom with --category filtering ---"

rm -f fatbands_0.dat fatbands_1.dat
timeout 120 stb-fatbands --label Sn3O4 --shift fermi --projection atom --category 0 1 \
    --save-gnuplot --no-intro > log_atom.txt 2>&1
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


# --- 5b. --view (headless via MPLBACKEND=Agg) ---
echo -e "\n--- Testing --view (headless via MPLBACKEND=Agg, only checking exit code) ---"

timeout 120 stb-fatbands --label Sn3O4 --shift fermi --view --no-intro > log_view.txt 2>&1
check_exit_code $? 0

echo "Testing: --view with a multi-category projection (plot_multi_series_on_bands path)"
timeout 120 stb-fatbands --label Sn3O4 --shift fermi --projection species --view --no-intro > log_view_multi.txt 2>&1
check_exit_code $? 0


# --- 6. Explicit --file/--wfsx/--hsx-file path (no --label) ---
echo -e "\n--- Testing explicit --file/--wfsx/--hsx-file (no --label) ---"

timeout 120 stb-fatbands --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX \
    --shift fermi --projection species -o explicit_out --save-gnuplot --no-intro > log_explicit.txt 2>&1
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
    --shift fermi --projection l -o approx_out2 --save-gnuplot --no-intro > log_approx2.txt 2>&1
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
check_contains "species_l" log_help.txt
check_contains "category" log_help.txt
check_contains "output-dir" log_help.txt
check_contains "save-report" log_help.txt
check_contains "save-gnuplot" log_help.txt
check_contains "\-\-view" log_help.txt


# --- 11. Interactive path (stb-suite, shortcut 3.10) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 3.10) ---"

echo "Testing: navigate 3.10 -> label Sn3O4 -> shift left blank (defaults to fermi) -> projection"
echo "left blank (defaults to species_l) -> default output -> no save-report -> save-gnuplot=y -> no view"
rm -f fatbands_Sn-s.dat fatbands_Sn-p.dat fatbands_Sn-d.dat \
      fatbands_O-s.dat fatbands_O-p.dat fatbands_O-d.dat fatbands.gplot stb_fatbands_report.txt
printf '3.10\nSn3O4\n\n\n\nn\ny\nn\n\n0\n' | timeout 120 stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_contains "Shift mode      : fermi" log_menu.txt
check_contains "Projection      : species_l" log_menu.txt
check_success fatbands_Sn-s.dat
check_success fatbands_O-p.dat
check_success fatbands.gplot
if [ -e stb_fatbands_report.txt ]; then
    echo -e " ... ${RED}FAIL${NC} (stb_fatbands_report.txt should not have been written -- save-report was 'n')"
    FAIL=$((FAIL+1))
else
    echo -e " ... ${GREEN}OK${NC} (stb_fatbands_report.txt correctly absent)"
    PASS=$((PASS+1))
fi

echo "Testing: navigate 3.10 -> projection l (option 1, non-default choice)"
rm -f fatbands_s.dat
printf '3.10\nSn3O4\n1\n1\n\nn\ny\nn\n\n0\n' | timeout 120 stb-suite > log_menu_l.txt 2>&1
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
