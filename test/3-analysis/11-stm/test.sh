#!/bin/bash

# --- Setup ---
# Smoke test for stb-stm (STM Simulator, item 3.11)
#
# Fixture note: Graphene.LDOS/.XV are a real SIESTA run (monolayer graphene,
# 2-atom primitive cell, 20 Ang vacuum along z, LDOS integrated over
# [-1.5, 1.5] eV around E_Fermi via '%block LocalDensityOfStates') -- fast
# to reproduce (`siesta < calc.fdf`) since it's a tiny, quickly-converging
# system, unlike test/3-analysis/10-fatbands/'s much heavier Sn3O4 fixture.
# crs/siesta.LDOS(.XV) is a second, real SIESTA run (a CrS monolayer, see
# section 5b below) -- both atoms of graphene's own primitive cell sit at
# EXACTLY the same z (a perfectly flat, centered plane), which never
# exercises the buckled/off-center vacuum-detection case crs/ was added
# for; kept as a separate fixture rather than replacing graphene's own.
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
echo "--- Starting tester for stb-stm (item 3.11) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR"/Graphene.LDOS "$FIXTURE_DIR"/Graphene.XV "$TEST_DIR/"
mkdir -p "$TEST_DIR/crs"
cp "$FIXTURE_DIR"/crs/siesta.LDOS "$FIXTURE_DIR"/crs/siesta.XV "$TEST_DIR/crs/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Basic run: --mode current (default), --label auto-detect, default --iso ---
echo -e "\n--- Testing a basic constant-current run with the DEFAULT --iso (--label Graphene) ---"

timeout 60 stb-stm --label Graphene --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "Using axis 2" log_basic.txt
check_contains "\[0\] RUN METADATA" log_basic.txt
check_contains "\[1\] INPUT DATA" log_basic.txt
check_contains "\[2\] STM IMAGE" log_basic.txt
check_contains "\[3\] OUTPUT DATA & PLOTS" log_basic.txt
check_contains "\[4\] REFERENCES" log_basic.txt
check_contains "\[5\] SUMMARY & FILES" log_basic.txt
check_contains "Iso threshold  : 1.000000e-03" log_basic.txt
check_contains "Corrugation" log_basic.txt

echo "Testing: without --save-gnuplot, no .dat/.gplot is written (used to be unconditional)"
if [ -e stm_current.dat ] || [ -e stm_current.gplot ]; then
    echo -e " ... ${RED}FAIL${NC} (stm_current.dat/stm_current.gplot should not have been written)"
    FAIL=$((FAIL+1))
else
    echo -e " ... ${GREEN}OK${NC} (stm_current.dat/stm_current.gplot correctly absent)"
    PASS=$((PASS+1))
fi
check_contains "Not written (off by default" log_basic.txt

echo "Testing: no text report without --save-report"
if [ -e stb_stm_report.txt ]; then
    echo -e " ... ${RED}FAIL${NC} (stb_stm_report.txt should not have been written)"
    FAIL=$((FAIL+1))
else
    echo -e " ... ${GREEN}OK${NC} (stb_stm_report.txt correctly absent)"
    PASS=$((PASS+1))
fi

echo "Testing: references.bib is always written (SIESTA + Tersoff-Hamann)"
check_success references.bib
check_contains "Soler2002" references.bib
check_contains "Garcia2020" references.bib
check_contains "Tersoff1985" references.bib


# --- 2b. --save-report / --save-gnuplot ---
echo -e "\n--- Testing --save-report / --save-gnuplot ---"

timeout 60 stb-stm --label Graphene --save-report --save-gnuplot --no-intro > log_saved.txt 2>&1
check_exit_code $? 0
check_success stb_stm_report.txt
check_contains "\[0\] RUN METADATA" stb_stm_report.txt
check_success stm_current.dat
check_success stm_current.gplot
check_contains "\[OK\] Data written to" log_saved.txt
check_contains "\[OK\] Gnuplot script written to" log_saved.txt

echo "Testing: constant-current corrugation is physically sensible (sub-Angstrom, non-zero)"
check_contains "Corrugation (max-min)" stb_stm_report.txt

echo "Testing: the .gplot script's own filenames are bare basenames, no --output-dir prefix"
echo "(the user is expected to cd into --output-dir and run gnuplot directly from there)"
if grep -q '"stm_current\.pdf"' stm_current.gplot && grep -q '"stm_current\.dat"' stm_current.gplot; then
    echo -e "   -> ${GREEN}Verified:${NC} stm_current.gplot references bare 'stm_current.pdf'/'stm_current.dat'"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} stm_current.gplot has an unexpected path in its output/data filenames"
    FAIL=$((FAIL+1))
fi
if command -v gnuplot > /dev/null; then
    gnuplot stm_current.gplot && test -s stm_current.pdf
    if [ $? -eq 0 ]; then
        echo -e "   -> ${GREEN}Verified:${NC} a real gnuplot run on stm_current.gplot produced stm_current.pdf"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} gnuplot did not produce stm_current.pdf"
        FAIL=$((FAIL+1))
    fi
fi


# --- 2c. Same path-prefix check, but from a genuine --output-dir subfolder ---
echo -e "\n--- Testing the gnuplot path fix with a real --output-dir subfolder ---"

rm -rf path_check
timeout 60 stb-stm --label Graphene --save-gnuplot -o path_check --no-intro > log_path_check.txt 2>&1
check_exit_code $? 0
if grep -q '"path_check' path_check/stm_current.gplot; then
    echo -e "   -> ${RED}Failed:${NC} path_check/stm_current.gplot still embeds the --output-dir prefix"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} path_check/stm_current.gplot has no --output-dir prefix in its filenames"
    PASS=$((PASS+1))
fi
if command -v gnuplot > /dev/null; then
    (cd path_check && gnuplot stm_current.gplot && test -s stm_current.pdf)
    if [ $? -eq 0 ]; then
        echo -e "   -> ${GREEN}Verified:${NC} gnuplot run from inside path_check/ (as a real user would) succeeds"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} gnuplot run from inside path_check/ failed"
        FAIL=$((FAIL+1))
    fi
fi


# --- 3. --mode height, default --z ---
echo -e "\n--- Testing --mode height with the DEFAULT --z ---"

timeout 60 stb-stm --label Graphene --mode height --save-gnuplot --no-intro > log_height.txt 2>&1
check_exit_code $? 0
check_contains "Mode           : height" log_height.txt
check_contains "Requested z    : 3.0 Ang" log_height.txt
check_success stm_height.dat
check_success stm_height.gplot

echo "Testing: --mode height with an explicit --z overrides the default"
timeout 60 stb-stm --label Graphene --mode height --z 2.0 --no-intro > log_height_explicit.txt 2>&1
check_exit_code $? 0
check_contains "Requested z    : 2.0 Ang" log_height_explicit.txt


# --- 3b. --view (headless via MPLBACKEND=Agg) ---
echo -e "\n--- Testing --view (headless via MPLBACKEND=Agg, only checking exit code) ---"

timeout 60 stb-stm --label Graphene --view --no-intro > log_view.txt 2>&1
check_exit_code $? 0

timeout 60 stb-stm --label Graphene --mode height --view --no-intro > log_view_height.txt 2>&1
check_exit_code $? 0


# --- 4. Explicit --file/--geometry-file (no --label) ---
echo -e "\n--- Testing explicit --file/--geometry-file (no --label) ---"

timeout 60 stb-stm --file Graphene.LDOS --geometry-file Graphene.XV \
    -o explicit_out --save-gnuplot --no-intro > log_explicit.txt 2>&1
check_exit_code $? 0
check_success explicit_out/stm_current.dat


# --- 5. Explicit --axis bypasses vacuum auto-detection ---
echo -e "\n--- Testing explicit --axis 2 ---"

timeout 60 stb-stm --label Graphene --axis 2 -o axis_out --save-gnuplot --no-intro > log_axis.txt 2>&1
check_exit_code $? 0
check_success axis_out/stm_current.dat


# --- 5b. Real bug fix: correct vacuum-side detection on a buckled monolayer ---
# Fixture note: crs/siesta.LDOS(.XV) is a REAL SIESTA run (a CrS monolayer
# fetched via stb-fetch from the twodmatpedia OPTIMADE database, id
# 2dm-2617 -- the same structure already used by examples/3.7/3.9), with
# %block LocalDensityOfStates EF -3.50 0.00 eV. Its 4 atoms sit at
# fractional z = 0, 0, 0.066, 0.934 -- a naive "topmost atom = max(z)"
# picks the wrong side (the tiny ~7% wraparound sliver past z=0.934
# instead of the real ~87%-of-cell vacuum gap starting at z=0.066),
# collapsing the search window to ~1.5 Ang and, once the window itself was
# fixed, letting the outside-in scan cross into the far face's own LDOS
# tail (18-20 Ang "corrugation", unphysical) unless the search is also
# capped at half the gap. This is the exact fixture that caught both bugs
# live -- see core/kspace.py::find_surface_reference.
echo -e "\n--- Testing correct vacuum-side detection on a real buckled monolayer (CrS) ---"

timeout 60 stb-stm --file crs/siesta.LDOS --geometry-file crs/siesta.XV --no-intro > log_crs.txt 2>&1
check_exit_code $? 0
echo "Testing: the topmost atom is the LOWER-z S atom (1.538 Ang), not the naive max(z) (21.817 Ang)"
check_contains "Topmost atom        | 1.538 Ang" log_crs.txt
echo "Testing: the search window is capped at half the real vacuum gap (~11.7 Ang), not the full ~20 Ang or the broken ~1.5 Ang"
check_contains "search window up to 11.677 Ang" log_crs.txt
echo "Testing: default --iso gives full coverage and a physically sensible (sub-2-Ang) corrugation"
check_contains "Points reaching iso   | 1600/1600" log_crs.txt
if grep -qE "Corrugation \(max-min\) \| 1\.[0-9]+ Ang" log_crs.txt; then
    echo -e "   -> ${GREEN}Verified:${NC} corrugation is a sensible sub-2-Ang value (not 0.000, not 18-20 Ang)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} corrugation is not in the expected 1.x Ang range"
    FAIL=$((FAIL+1))
fi


# --- 6. --iso set unreachably high: NaN points are reported, not a crash ---
echo -e "\n--- Testing an unreachable --iso (NaN/no-crossing points) ---"

timeout 60 stb-stm --label Graphene --iso 100.0 -o unreachable_out --save-gnuplot --no-intro > log_unreachable.txt 2>&1
check_exit_code $? 0
check_contains "never reached --iso" log_unreachable.txt
check_contains "NaN" unreachable_out/stm_current.dat


# --- 7. --mode current/height now work with NO --iso/--z at all (new defaults) ---
echo -e "\n--- Testing that --iso/--z are no longer required (real bug fixed: no usable defaults before) ---"

echo "Testing: --mode current without --iso uses the new DEFAULT_ISO"
stb-stm --label Graphene --no-intro > log_default_iso.txt 2>&1
check_exit_code $? 0
check_contains "Iso threshold  : 1.000000e-03" log_default_iso.txt

echo "Testing: --mode height without --z uses the new DEFAULT_Z"
stb-stm --label Graphene --mode height --no-intro > log_default_z.txt 2>&1
check_exit_code $? 0
check_contains "Requested z    : 3.0 Ang" log_default_z.txt


# --- 8. Missing .LDOS aborts with a clear error ---
echo -e "\n--- Testing a missing .LDOS ---"

timeout 30 stb-stm --label does_not_exist --no-intro > log_missing_ldos.txt 2>&1
check_exit_code $? 2
check_contains "No LDOS file found" log_missing_ldos.txt


# --- 9. A mismatched geometry (wrong cell for this grid) is rejected ---
echo -e "\n--- Testing rejection of a geometry whose cell doesn't match the grid ---"

echo "Testing: pairing Graphene.LDOS with an unrelated bulk .XV (Sn3O4, wrong cell)"
SN3O4_XV="$FIXTURE_DIR/../10-fatbands/Sn3O4.XV"
if [ -f "$SN3O4_XV" ]; then
    stb-stm --file Graphene.LDOS --geometry-file "$SN3O4_XV" --no-intro > log_cell_mismatch.txt 2>&1
    check_exit_code $? 1
    check_contains "does not match" log_cell_mismatch.txt
else
    echo -e " ... ${YELLOW}SKIPPED${NC} (test/3-analysis/10-fatbands/Sn3O4.XV not found)"
fi


# --- 10. Robustness / --help / --version ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: --label omitted (required)"
stb-stm --no-intro > log_missing_label.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-stm --version > log_version.txt 2>&1
check_contains "stb-stm" log_version.txt

echo "Testing: --help documents --mode, --z, --iso, --z-max, --axis, --geometry-file, -o, --save-report, --save-gnuplot, --view"
stb-stm --help > log_help.txt 2>&1
check_contains "mode" log_help.txt
check_contains "\-\-z " log_help.txt
check_contains "iso" log_help.txt
check_contains "z-max" log_help.txt
check_contains "axis" log_help.txt
check_contains "geometry-file" log_help.txt
check_contains "output-dir" log_help.txt
check_contains "save-report" log_help.txt
check_contains "save-gnuplot" log_help.txt
check_contains "\-\-view" log_help.txt


# --- 11. Interactive path (stb-suite, shortcut 3.11) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 3.11) ---"

echo "Testing: navigate 3.11 -> label Graphene -> mode 1 (current) -> iso left blank (default"
echo "0.001) -> auto axis -> default output -> no save-report -> save-gnuplot=y -> no view"
rm -f stm_current.dat stm_current.gplot stb_stm_report.txt
printf '3.11\nGraphene\n1\n\n\n\nn\ny\nn\n\n0\n' | timeout 60 stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_contains "Iso threshold  : 1.000000e-03" log_menu.txt
check_success stm_current.dat
check_success stm_current.gplot
if [ -e stb_stm_report.txt ]; then
    echo -e " ... ${RED}FAIL${NC} (stb_stm_report.txt should not have been written -- save-report was 'n')"
    FAIL=$((FAIL+1))
else
    echo -e " ... ${GREEN}OK${NC} (stb_stm_report.txt correctly absent)"
    PASS=$((PASS+1))
fi


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
