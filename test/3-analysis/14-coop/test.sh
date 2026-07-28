#!/bin/bash

# --- Setup ---
# Smoke test for stb-coop (COOP/COHP Bonding Analysis, item 3.14)
#
# Fixture note: Sn3O4.selected.WFSX + Sn3O4.HSX are a real SIESTA run
# (restarted from the converged Sn3O4.DM in test/6-utils/4-wantibexos/, but
# with a small explicit 2x2x2-equivalent k-mesh via WriteWaveFunctions T +
# %block WaveFuncKPoints instead of a band path -- COOP/COHP needs a real
# BZ sampling, not a band-path WFSX like stb-fatbands uses). Small (~6 MB)
# since this is only 8 k-points. calc.fdf's real fdf filename does not
# match SystemLabel (Sn3O4) -- same pattern as test/3-analysis/10-fatbands/,
# so tests use --wfsx/--hsx-file explicitly where the --label path can't
# find <label>.fdf (not needed here since coop.py never reads .fdf directly).
# Sn3O4.DM (same test/6-utils/4-wantibexos/ calc) is the real density matrix
# needed for --bond-order -- see the module docstring for the real bug this
# fixture caught (H.bond_order() called on the Hamiltonian instead of a real
# density matrix, giving numbers off by orders of magnitude).
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
echo "--- Starting tester for stb-coop (item 3.14) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR"/Sn3O4.selected.WFSX "$FIXTURE_DIR"/Sn3O4.HSX "$FIXTURE_DIR"/Sn3O4.XV \
   "$FIXTURE_DIR"/Sn3O4.DM "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null

# Keep runtimes reasonable: --npoints modest, since sisl's COOP/COHP is
# memory/CPU-heavy per its own docstring warning.
COMMON="--erange -10 5 --npoints 25 --sigma 300"
# 0-6 is a real Sn-O bond (~2.08 Ang); 0-1 is Sn-Sn 3.8 Ang apart (not
# bonded) -- see the module docstring's live verification of both the
# COOP sign convention and the --bond-order fix.
# -25 to the real Fermi energy (-3.203856 eV, from this calc's own .out)
# spans (approximately) the occupied states, needed for the bond-order
# sign cross-check to be a fair comparison (see report section [4]'s own
# in-tool caveat about --erange coverage).
OCCUPIED="--erange -25 -3.203856 --npoints 40 --sigma 300"


# --- 2. Basic run: --quantity coop, --pair (label auto-detect via --label) ---
echo -e "\n--- Testing a basic COOP run (--label Sn3O4 --pair 0 1), numbered report ---"

timeout 120 stb-coop --label Sn3O4 --quantity coop --pair 0 1 $COMMON --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "Using 'Sn3O4.HSX'" log_basic.txt
check_contains "Selected pairs: \['0-1'\]" log_basic.txt
check_success coop.dat
check_contains "\[0\] RUN METADATA" log_basic.txt
check_contains "\[1\] INPUT DATA" log_basic.txt
check_contains "\[2\] PAIR SELECTION" log_basic.txt
check_contains "\[3\] COOP CURVE" log_basic.txt
check_contains "\[4\] BOND ORDER" log_basic.txt
check_contains "\[5\] OUTPUT DATA & PLOTS" log_basic.txt
check_contains "\[6\] REFERENCES" log_basic.txt
check_contains "\[7\] SUMMARY & FILES" log_basic.txt
check_contains "Not requested (pass --bond-order" log_basic.txt


# --- 3. --quantity cohp + --pair-species ---
echo -e "\n--- Testing --quantity cohp with --pair-species ---"

timeout 120 stb-coop --label Sn3O4 --quantity cohp --pair-species Sn O $COMMON -o cohp_out --no-intro > log_cohp.txt 2>&1
check_exit_code $? 0
check_contains "Selected pairs: \['Sn-O'\]" log_cohp.txt
check_success cohp_out/cohp.dat


# --- 4. Multiple --pair given together ---
echo -e "\n--- Testing multiple --pair (repeatable) ---"

timeout 120 stb-coop --label Sn3O4 --quantity coop --pair 0 1 --pair 6 7 $COMMON -o multi_out --no-intro > log_multi.txt 2>&1
check_exit_code $? 0
check_success multi_out/coop.dat
echo "Testing: both pair columns present in the header"
check_contains "0-1" multi_out/coop.dat
check_contains "6-7" multi_out/coop.dat


# --- 5. --bond-order: the real bug fix, verified against a genuinely bonded
# and a genuinely non-bonded pair ---
echo -e "\n--- Testing --bond-order (fixed: reads a real .DM, not H.bond_order()) ---"

timeout 120 stb-coop --label Sn3O4 --quantity coop --pair 0 6 --pair 0 1 $OCCUPIED \
    --bond-order -o bo_out --no-intro > log_bo.txt 2>&1
check_exit_code $? 0
check_contains "Using 'Sn3O4.DM' for --bond-order" log_bo.txt
check_contains "\[4\] BOND ORDER" log_bo.txt
check_contains "signs agree" log_bo.txt

echo "Testing: the bonded pair's bond order is a physically sensible O(0.1-1) value, not O(10-100)"
python3 -c "
import re, sys
text = open('log_bo.txt').read()
section = text.split('BOND ORDER')[1]
m = re.search(r'0-6\s+\|\s+(-?[\d.]+)', section)
assert m, 'bond order for 0-6 not found in the BOND ORDER report section'
value = float(m.group(1))
sys.exit(0 if 0.1 < value < 2.0 else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} bonded pair's Mulliken bond order is in the physically sensible 0.1-2.0 range"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} bonded pair's Mulliken bond order is NOT in the sensible range"
    FAIL=$((FAIL+1))
fi

echo "Testing: --bond-order without any usable .DM is a clean error, not a silently wrong number"
mv Sn3O4.DM Sn3O4.DM.hidden
stb-coop --label Sn3O4 --quantity coop --pair 0 1 $COMMON --bond-order --no-intro > log_bo_missing.txt 2>&1
check_exit_code $? 2
check_contains "needs a real density matrix" log_bo_missing.txt
mv Sn3O4.DM.hidden Sn3O4.DM

echo "Testing: explicit --dm-file (renamed .DM, decoupled from --label)"
cp Sn3O4.DM my_density.DM
stb-coop --label Sn3O4 --quantity coop --pair 0 6 $OCCUPIED --bond-order --dm-file my_density.DM \
    -o dmfile_out --no-intro > log_dmfile.txt 2>&1
check_exit_code $? 0
check_contains "Using 'my_density.DM' for --bond-order" log_dmfile.txt
rm -f my_density.DM


# --- 6. Explicit --wfsx/--hsx-file (no --label) ---
echo -e "\n--- Testing explicit --wfsx/--hsx-file (no --label) ---"

timeout 120 stb-coop --wfsx Sn3O4.selected.WFSX --hsx-file Sn3O4.HSX --quantity coop --pair 0 1 $COMMON -o explicit_out --no-intro > log_explicit.txt 2>&1
check_exit_code $? 0
check_success explicit_out/coop.dat


# --- 7. --label + --hsx-file together (real bug: used to be rejected outright) ---
echo -e "\n--- Testing --label + --hsx-file together (a renamed .HSX, decoupled from --label) ---"

cp Sn3O4.HSX calc.HSX
timeout 120 stb-coop --label Sn3O4 --hsx-file calc.HSX --quantity coop --pair 0 1 $COMMON \
    -o labelhsx_out --no-intro > log_labelhsx.txt 2>&1
check_exit_code $? 0
check_contains "Hamiltonian src: calc.HSX" log_labelhsx.txt
rm -f calc.HSX

echo "Testing: --label + --wfsx together is still rejected (ambiguous)"
stb-coop --label Sn3O4 --wfsx Sn3O4.selected.WFSX --quantity coop --pair 0 1 $COMMON --no-intro > log_labelwfsx.txt 2>&1
check_exit_code $? 2
check_contains "cannot be combined with --wfsx" log_labelwfsx.txt


# --- 8. --shift fermi: --fermi, --bands-file, --fermi-file, and .out auto-detection ---
echo -e "\n--- Testing --shift fermi (--fermi explicit value) ---"

timeout 120 stb-coop --label Sn3O4 --quantity coop --pair 0 1 --erange -3 3 --npoints 20 --sigma 300 \
    --shift fermi --fermi -3.203856 -o fermi_out --no-intro > log_fermi.txt 2>&1
check_exit_code $? 0
check_success fermi_out/coop.dat
check_contains "Fermi source: --fermi (explicit value)" log_fermi.txt

echo -e "\n--- Testing --shift fermi with --fermi-file (an arbitrarily-named .out log) ---"
cat > my_weird_name.out << 'EOF'
Some SIESTA log noise before the summary
siesta:         Fermi = -3.203856
More noise after
EOF
timeout 120 stb-coop --label Sn3O4 --quantity coop --pair 0 1 --erange -3 3 --npoints 20 --sigma 300 \
    --shift fermi --fermi-file my_weird_name.out -o fermifile_out --no-intro > log_fermifile.txt 2>&1
check_exit_code $? 0
check_contains "Fermi source: 'my_weird_name.out' (--fermi-file)" log_fermifile.txt

echo -e "\n--- Testing --shift fermi auto-detecting a generic .out ---"
timeout 120 stb-coop --label Sn3O4 --quantity coop --pair 0 1 --erange -3 3 --npoints 20 --sigma 300 \
    --shift fermi -o fermiauto_out --no-intro > log_fermiauto.txt 2>&1
check_exit_code $? 0
check_contains "auto-detected .out" log_fermiauto.txt
rm -f my_weird_name.out

echo "Testing: --shift fermi with no Fermi source at all is a clean error"
stb-coop --label Sn3O4 --quantity coop --pair 0 1 $COMMON --shift fermi --no-intro > log_missing_fermi.txt 2>&1
check_exit_code $? 2
check_contains "needs a Fermi energy" log_missing_fermi.txt


# --- 9. --save-report / --save-gnuplot (off by default) ---
echo -e "\n--- Testing --save-report / --save-gnuplot ---"

timeout 120 stb-coop --label Sn3O4 --quantity coop --pair 0 1 $COMMON \
    --save-report --save-gnuplot -o report_out --no-intro > log_report.txt 2>&1
check_exit_code $? 0
check_success report_out/coop.dat
check_success report_out/coop.gplot
check_success report_out/stb_coop_report.txt
check_success report_out/references.bib

echo "Testing: coop.gplot uses a bare basename (not the full --output-dir path)"
if grep -q '"coop.dat"' report_out/coop.gplot && ! grep -q "report_out" report_out/coop.gplot; then
    echo -e "   -> ${GREEN}Verified:${NC} coop.gplot references bare basenames, no leaked output-dir path"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} coop.gplot leaked the --output-dir path"
    FAIL=$((FAIL+1))
fi

echo "Testing: the saved .gplot script actually renders with a real gnuplot"
if command -v gnuplot >/dev/null 2>&1; then
    (cd report_out && gnuplot coop.gplot)
    check_success report_out/coop.pdf
else
    echo -e "   -> ${YELLOW}Skipped:${NC} gnuplot not installed in this environment"
fi

echo "Testing: a plain run with no --view/--save-report/--save-gnuplot only writes coop.dat + references.bib"
rm -rf plain_out
timeout 120 stb-coop --label Sn3O4 --quantity coop --pair 0 1 $COMMON -o plain_out --no-intro > log_plain.txt 2>&1
check_exit_code $? 0
check_success plain_out/coop.dat
check_success plain_out/references.bib
if [ ! -f plain_out/coop.gplot ] && [ ! -f plain_out/stb_coop_report.txt ]; then
    echo -e "   -> ${GREEN}Verified:${NC} no .gplot/report written without --save-gnuplot/--save-report"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} .gplot or report file was written despite being off by default"
    FAIL=$((FAIL+1))
fi


# --- 10. --view (replaces the old, on-by-default --no-plot) ---
echo -e "\n--- Testing --view (off by default, opt-in matplotlib preview) ---"

timeout 120 stb-coop --label Sn3O4 --quantity coop --pair 0 1 $COMMON --view -o view_out --no-intro > log_view.txt 2>&1
check_exit_code $? 0
check_success view_out/coop.dat


# --- 11. Robustness / errors ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: --wfsx without --hsx-file"
stb-coop --wfsx Sn3O4.selected.WFSX --quantity coop --pair 0 1 $COMMON --no-intro > log_missing_hsx.txt 2>&1
check_exit_code $? 2

echo "Testing: missing --pair/--pair-species"
stb-coop --label Sn3O4 --quantity coop $COMMON --no-intro > log_missing_pair.txt 2>&1
check_exit_code $? 2

echo "Testing: --pair index out of range"
stb-coop --label Sn3O4 --quantity coop --pair 0 999 $COMMON --no-intro > log_bad_pair.txt 2>&1
check_exit_code $? 2
check_contains "out of range" log_bad_pair.txt

echo "Testing: missing --quantity"
stb-coop --label Sn3O4 --pair 0 1 $COMMON --no-intro > log_missing_quantity.txt 2>&1
check_exit_code $? 2

echo "Testing: missing --label/--wfsx"
stb-coop --quantity coop --pair 0 1 $COMMON --no-intro > log_missing_label.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-coop --version > log_version.txt 2>&1
check_contains "stb-coop" log_version.txt

echo "Testing: --help documents --quantity, --pair, --pair-species, --erange, --sigma, --bond-order, --dm-file, --save-report, --save-gnuplot, --view"
stb-coop --help > log_help.txt 2>&1
check_contains "quantity" log_help.txt
check_contains "pair-species" log_help.txt
check_contains "erange" log_help.txt
check_contains "sigma" log_help.txt
check_contains "bond-order" log_help.txt
check_contains "dm-file" log_help.txt
check_contains "save-report" log_help.txt
check_contains "save-gnuplot" log_help.txt
check_contains "view" log_help.txt


# --- 12. Interactive path (stb-suite, shortcut 3.14) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 3.14) ---"

echo "Testing: navigate 3.14 -> label Sn3O4 -> default HSX -> quantity 1 (coop) -> pair mode 1 (index) -> 0 1 -> erange -10 5 -> sigma 300 -> no shift -> no bond-order -> default output -> no report/gnuplot/view"
rm -f coop.dat
printf '3.14\nSn3O4\n\n1\n1\n0\n1\n-10\n5\n300\n\nn\n\nn\nn\nn\n\n0\n' | timeout 120 stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success coop.dat

echo "Testing: interactive menu's --bond-order + default .DM prompt"
rm -f coop.dat
printf '3.14\nSn3O4\n\n1\n1\n0\n6\n-25\n-3.203856\n300\n\ny\n\n\nn\nn\nn\n\n0\n' | timeout 120 stb-suite > log_menu_bo.txt 2>&1
check_exit_code $? 0
check_contains "signs agree" log_menu_bo.txt


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
