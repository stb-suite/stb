#!/bin/bash

# --- Setup ---
# Smoke test for stb-sts (STS Spectroscopy Simulator, item 3.13)
#
# Fixture note: Graphene.selected.WFSX is a real SIESTA run (monolayer
# graphene, 9 explicit k-points via %block WaveFuncKPoints + WriteWaveFunctions T
# -- SIESTA's own naming convention for this mechanism, confirmed live: bare
# WriteWaveFunctions T with no WaveFuncKPoints block silently writes nothing).
# calc.fdf's real fdf filename does not match SystemLabel (Graphene) -- same
# pattern as test/3-analysis/10-fatbands/, so tests use --wfsx/--geometry-file
# explicitly rather than --label.
#
# Important physics note verified live: SIESTA's confined PAO basis has a
# finite cutoff radius (~2.6 Ang for this DZP carbon basis) -- a tip height
# much beyond that gives an exactly-zero (not just small) curve. The basic
# test below uses --height 1.5 (well within the cutoff) to get a real,
# non-trivial spectrum; a separate test exercises the "too far" warning path.
export MPLBACKEND=Agg
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
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


# --- 1. Preparation ---
echo "--- Starting tester for stb-sts (item 3.13) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR"/Graphene.selected.WFSX "$FIXTURE_DIR"/Graphene.XV "$FIXTURE_DIR"/calc.fdf \
   "$FIXTURE_DIR"/structure.fdf "$FIXTURE_DIR"/C.ion "$FIXTURE_DIR"/C.ion.xml "$TEST_DIR/"
# Same-name fdf copy for the interactive-menu test (--geometry-file default there).
cp "$FIXTURE_DIR"/calc.fdf "$TEST_DIR/Graphene.fdf"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Basic run: --xy/--height tip position within the PAO cutoff ---
echo -e "\n--- Testing a basic run (--xy 0 0 --height 1.5) ---"

timeout 60 stb-sts --wfsx Graphene.selected.WFSX --geometry-file calc.fdf \
    --xy 0 0 --height 1.5 --erange -3 3 --sigma 100 --shift fermi --fermi -4.461866 \
    --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "Using axis 2 as the surface normal" log_basic.txt
check_contains "Accumulated 234" log_basic.txt
check_success sts.dat

echo "Testing: the curve has at least one non-trivial (non-zero) point"
python3 -c "
import numpy as np, sys
data = np.loadtxt('sts.dat')
sys.exit(0 if (data[:,1] > 1e-6).any() else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} sts.dat has non-trivial signal"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} sts.dat is all-zero"
    FAIL=$((FAIL+1))
fi


# --- 3. Tip point beyond the PAO cutoff radius: warned, not silently wrong ---
echo -e "\n--- Testing a tip point beyond the confined-basis cutoff (--height 3.0) ---"

timeout 60 stb-sts --wfsx Graphene.selected.WFSX --geometry-file calc.fdf \
    --xy 0 0 --height 3.0 --erange -3 3 --sigma 100 -o farout_out --no-intro > log_farout.txt 2>&1
check_exit_code $? 0
check_contains "exactly zero for every state" log_farout.txt


# --- 4. --point (absolute Cartesian) instead of --xy/--height ---
echo -e "\n--- Testing --point ---"

timeout 60 stb-sts --wfsx Graphene.selected.WFSX --geometry-file calc.fdf \
    --point 0 0 11.5 --erange -3 3 --fwhm 200 -o point_out --no-intro > log_point.txt 2>&1
check_exit_code $? 0
check_success point_out/sts.dat


# --- 5. --view (replaces the old, on-by-default --no-plot) ---
echo -e "\n--- Testing --view (off by default, opt-in matplotlib preview) ---"

timeout 60 stb-sts --wfsx Graphene.selected.WFSX --geometry-file calc.fdf \
    --xy 0 0 --height 1.5 --erange -3 3 --sigma 100 --view -o view_out --no-intro > log_view.txt 2>&1
check_exit_code $? 0
check_success view_out/sts.dat

echo "Testing: a plain run with no --view/--save-report/--save-gnuplot only writes sts.dat + references.bib"
rm -rf plain_out
timeout 60 stb-sts --wfsx Graphene.selected.WFSX --geometry-file calc.fdf \
    --xy 0 0 --height 1.5 --erange -3 3 --sigma 100 -o plain_out --no-intro > log_plain.txt 2>&1
check_exit_code $? 0
check_success plain_out/sts.dat
check_success plain_out/references.bib
if [ ! -f plain_out/sts.gplot ] && [ ! -f plain_out/stb_sts_report.txt ]; then
    echo -e "   -> ${GREEN}Verified:${NC} no .gplot/report written without --save-gnuplot/--save-report"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} .gplot or report file was written despite being off by default"
    FAIL=$((FAIL+1))
fi


# --- 6. --save-report / --save-gnuplot (numbered [0]...[6] report, real .gplot script) ---
echo -e "\n--- Testing --save-report / --save-gnuplot ---"

timeout 60 stb-sts --wfsx Graphene.selected.WFSX --geometry-file calc.fdf \
    --xy 0 0 --height 1.5 --erange -3 3 --sigma 100 --save-report --save-gnuplot \
    -o report_out --no-intro > log_report.txt 2>&1
check_exit_code $? 0
check_success report_out/sts.dat
check_success report_out/sts.gplot
check_success report_out/stb_sts_report.txt
check_success report_out/references.bib
check_contains "\[0\] RUN METADATA" report_out/stb_sts_report.txt
check_contains "\[1\] INPUT DATA" report_out/stb_sts_report.txt
check_contains "\[2\] TIP POSITION" report_out/stb_sts_report.txt
check_contains "\[3\] STS CURVE" report_out/stb_sts_report.txt
check_contains "\[4\] OUTPUT DATA & PLOTS" report_out/stb_sts_report.txt
check_contains "\[5\] REFERENCES" report_out/stb_sts_report.txt
check_contains "\[6\] SUMMARY & FILES" report_out/stb_sts_report.txt

echo "Testing: sts.gplot uses a bare basename (not the full --output-dir path) in its plot/output lines"
if grep -q '"sts.dat"' report_out/sts.gplot && ! grep -q "report_out" report_out/sts.gplot; then
    echo -e "   -> ${GREEN}Verified:${NC} sts.gplot references bare basenames, no leaked output-dir path"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} sts.gplot leaked the --output-dir path (same bug class fixed in stm.py/grid_export.py)"
    FAIL=$((FAIL+1))
fi

echo "Testing: the saved .gplot script actually renders with a real gnuplot"
if command -v gnuplot >/dev/null 2>&1; then
    (cd report_out && gnuplot sts.gplot)
    check_success report_out/sts.pdf
else
    echo -e "   -> ${YELLOW}Skipped:${NC} gnuplot not installed in this environment"
fi


# --- 7. --shift fermi: --fermi, --bands-file, --fermi-file, and .out auto-detection ---
echo -e "\n--- Testing --shift fermi (--fermi explicit value) ---"

timeout 60 stb-sts --wfsx Graphene.selected.WFSX --geometry-file calc.fdf \
    --xy 0 0 --height 1.5 --erange -3 3 --sigma 100 --shift fermi --fermi -4.461866 \
    -o fermi_out --no-intro > log_fermi.txt 2>&1
check_exit_code $? 0
check_contains "Fermi source: --fermi (explicit value)" log_fermi.txt

echo -e "\n--- Testing --shift fermi with --fermi-file (an arbitrarily-named .out log, decoupled from --label) ---"
cat > my_weird_name.out << 'EOF'
Some SIESTA log noise before the summary
siesta:         Fermi = -4.461866
More noise after
EOF
timeout 60 stb-sts --wfsx Graphene.selected.WFSX --geometry-file calc.fdf \
    --xy 0 0 --height 1.5 --erange -3 3 --sigma 100 --shift fermi --fermi-file my_weird_name.out \
    -o fermifile_out --no-intro > log_fermifile.txt 2>&1
check_exit_code $? 0
check_contains "Fermi source: 'my_weird_name.out' (--fermi-file)" log_fermifile.txt

echo -e "\n--- Testing --shift fermi auto-detecting a generic .out (no --label match needed) ---"
timeout 60 stb-sts --wfsx Graphene.selected.WFSX --geometry-file calc.fdf \
    --xy 0 0 --height 1.5 --erange -3 3 --sigma 100 --shift fermi \
    -o fermiauto_out --no-intro > log_fermiauto.txt 2>&1
check_exit_code $? 0
check_contains "auto-detected .out" log_fermiauto.txt
rm -f my_weird_name.out

echo "Testing: --shift fermi with no Fermi source at all is a clean error"
stb-sts --wfsx Graphene.selected.WFSX --geometry-file calc.fdf --xy 0 0 --height 1.5 --erange -3 3 --sigma 50 --shift fermi --no-intro > log_missing_fermi.txt 2>&1
check_exit_code $? 2
check_contains "needs a Fermi energy" log_missing_fermi.txt


# --- 8. --label + --geometry-file (real bug: used to be rejected outright) ---
echo -e "\n--- Testing --label + --geometry-file together (SystemLabel != real fdf filename) ---"

timeout 60 stb-sts --label Graphene --geometry-file calc.fdf \
    --xy 0 0 --height 1.5 --erange -3 3 --sigma 100 -o labelgeom_out --no-intro > log_labelgeom.txt 2>&1
check_exit_code $? 0
check_contains "Geometry source: calc.fdf" log_labelgeom.txt

echo "Testing: --label + --wfsx together is still rejected (ambiguous)"
stb-sts --label Graphene --wfsx Graphene.selected.WFSX --xy 0 0 --height 1.5 --erange -3 3 --sigma 50 --no-intro > log_labelwfsx.txt 2>&1
check_exit_code $? 2
check_contains "cannot be combined with --wfsx" log_labelwfsx.txt


# --- 9. Robustness / errors ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: --xy without --height"
stb-sts --wfsx Graphene.selected.WFSX --geometry-file calc.fdf --xy 0 0 --erange -3 3 --sigma 50 --no-intro > log_missing_height.txt 2>&1
check_exit_code $? 2

echo "Testing: missing --wfsx/--label"
stb-sts --xy 0 0 --height 1.5 --erange -3 3 --sigma 50 --no-intro > log_missing_label.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-sts --version > log_version.txt 2>&1
check_contains "stb-sts" log_version.txt

echo "Testing: --help documents --xy, --point, --height, --erange, --sigma, --fwhm, --shift, --save-report, --save-gnuplot, --view"
stb-sts --help > log_help.txt 2>&1
check_contains "xy" log_help.txt
check_contains "point" log_help.txt
check_contains "height" log_help.txt
check_contains "erange" log_help.txt
check_contains "sigma" log_help.txt
check_contains "fwhm" log_help.txt
check_contains "shift" log_help.txt
check_contains "save-report" log_help.txt
check_contains "save-gnuplot" log_help.txt
check_contains "view" log_help.txt


# --- 10. Interactive path (stb-suite, shortcut 3.13) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 3.13) ---"

echo "Testing: navigate 3.13 -> label Graphene -> default fdf -> tip mode 1 (xy+height) -> 0 0 1.5 -> erange -3 3 -> sigma 100 -> no shift -> default output -> no report/gnuplot/view"
rm -f sts.dat
printf '3.13\nGraphene\n\n1\n0\n0\n1.5\n-3\n3\n100\n\n\nn\nn\nn\n\n0\n' | timeout 60 stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success sts.dat

echo "Testing: interactive menu's Fermi-source submenu (option 4 -> auto-detect a .out)"
cat > Graphene.out << 'EOF'
siesta:         Fermi = -4.461866
EOF
rm -f sts.dat
printf '3.13\nGraphene\n\n1\n0\n0\n1.5\n-3\n3\n100\n2\n4\n\nn\nn\nn\n\n0\n' | timeout 60 stb-suite > log_menu_fermi.txt 2>&1
check_exit_code $? 0
check_contains "auto-detected .out" log_menu_fermi.txt
rm -f Graphene.out


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
