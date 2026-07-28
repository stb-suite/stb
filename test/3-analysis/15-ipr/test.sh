#!/bin/bash

# --- Setup ---
# Smoke test for stb-ipr (IPR Analyzer, item 3.15)
#
# Reuses test/3-analysis/10-fatbands/'s Sn3O4 fixture (calc.fdf + .ion files
# + Sn3O4.bands.WFSX + Sn3O4.HSX) -- IPR only needs orbital expansion
# coefficients (norm2/hadamard), the exact same input stb-fatbands already
# uses, so no new SIESTA run is needed. Also reuses 10-fatbands/spin/'s real
# spin-polarized isolated-O-atom fixture (Ospin.*) for the nspin=2
# regression case (same fixture stb-fatbands' own rewrite verified live).
export MPLBACKEND=Agg
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="$FIXTURE_DIR/../10-fatbands"
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
echo "--- Starting tester for stb-ipr (item 3.15) ---"
if [ ! -f "$SOURCE_DIR/Sn3O4.bands.WFSX" ]; then
    echo -e "${RED}FATAL: $SOURCE_DIR/Sn3O4.bands.WFSX not found -- run test/3-analysis/10-fatbands/test.sh's fixture setup first.${NC}"
    exit 1
fi
if [ ! -f "$SOURCE_DIR/spin/Ospin.bands.WFSX" ]; then
    echo -e "${RED}FATAL: $SOURCE_DIR/spin/Ospin.bands.WFSX not found.${NC}"
    exit 1
fi
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$SOURCE_DIR"/Sn3O4.bands "$SOURCE_DIR"/Sn3O4.bands.WFSX "$SOURCE_DIR"/Sn3O4.HSX \
   "$SOURCE_DIR"/calc.fdf "$SOURCE_DIR"/structure.fdf "$SOURCE_DIR"/Sn.ion "$SOURCE_DIR"/Sn.ion.xml \
   "$SOURCE_DIR"/O.ion "$SOURCE_DIR"/O.ion.xml "$TEST_DIR/"
# Same-name fdf copy for the interactive-menu test (default .HSX/.fdf
# suggestion by exact name).
cp "$SOURCE_DIR"/calc.fdf "$TEST_DIR/Sn3O4.fdf"
mkdir -p "$TEST_DIR/spin"
cp "$SOURCE_DIR"/spin/Ospin.bands "$SOURCE_DIR"/spin/Ospin.bands.WFSX "$SOURCE_DIR"/spin/Ospin.HSX \
   "$SOURCE_DIR"/spin/Ospin.fdf "$SOURCE_DIR"/spin/O.ion "$SOURCE_DIR"/spin/O.ion.xml "$TEST_DIR/spin/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Basic run: --file/--wfsx/--hsx-file explicit (accurate, overlap-aware), numbered report ---
echo -e "\n--- Testing a basic run (--hsx-file, accurate path), numbered report ---"

timeout 60 stb-ipr --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX \
    --shift fermi --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "Using 'Sn3O4.HSX' for overlap-aware IPR" log_basic.txt
check_contains "correspondence check passed" log_basic.txt
check_contains "\[0\] RUN METADATA" log_basic.txt
check_contains "\[1\] INPUT DATA" log_basic.txt
check_contains "\[2\] BAND GAP ANALYSIS" log_basic.txt
check_contains "\[3\] IPR ANALYSIS" log_basic.txt
check_contains "\[4\] WRITING OUTPUT FILES" log_basic.txt
check_contains "\[5\] REFERENCES" log_basic.txt
check_contains "\[6\] SUMMARY & FILES" log_basic.txt
check_contains "no universal absolute scale" log_basic.txt

echo "Testing: a plain run with no --view/--save-report/--save-gnuplot only writes references.bib"
if [ ! -f ipr.dat ] && [ ! -f ipr.gplot ] && [ ! -f stb_ipr_report.txt ] && [ -f references.bib ]; then
    echo -e "   -> ${GREEN}Verified:${NC} no ipr.dat/.gplot/report written without --save-gnuplot/--save-report"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} ipr.dat/.gplot/report was written despite being off by default"
    FAIL=$((FAIL+1))
fi
rm -f references.bib

echo "Testing: deep (core-like) states have higher IPR (more localized) than near-Fermi states"
timeout 60 stb-ipr --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX \
    --shift fermi --save-gnuplot --no-intro > log_basic_gnuplot.txt 2>&1
check_exit_code $? 0
check_success ipr.dat
check_success ipr.gplot
python3 -c "
import numpy as np, sys
data = np.loadtxt('ipr.dat')
energy, ipr = data[:,1], data[:,2]
deep = ipr[energy < -15]
near_fermi = ipr[(energy > -2) & (energy < 2)]
sys.exit(0 if deep.mean() > near_fermi.mean() else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} deep-state IPR > near-Fermi IPR (physically sensible)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} deep-state IPR did not exceed near-Fermi IPR"
    FAIL=$((FAIL+1))
fi

echo "Testing: ipr.gplot uses a bare basename (not the full --output-dir path)"
if grep -q '"ipr.dat"' ipr.gplot; then
    echo -e "   -> ${GREEN}Verified:${NC} ipr.gplot references bare basenames"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} ipr.gplot did not reference a bare basename"
    FAIL=$((FAIL+1))
fi

echo "Testing: the saved .gplot script actually renders with a real gnuplot"
if command -v gnuplot >/dev/null 2>&1; then
    gnuplot ipr.gplot
    check_success ipr.pdf
else
    echo -e "   -> ${YELLOW}Skipped:${NC} gnuplot not installed in this environment"
fi
rm -f ipr.dat ipr.gplot ipr.pdf references.bib


# --- 3. --save-report / --view ---
echo -e "\n--- Testing --save-report / --view ---"

timeout 60 stb-ipr --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX \
    --shift fermi --save-report -o savereport_out --no-intro > log_savereport.txt 2>&1
check_exit_code $? 0
check_success savereport_out/stb_ipr_report.txt

timeout 60 stb-ipr --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX \
    --shift fermi --view -o view_out --no-intro > log_view.txt 2>&1
check_exit_code $? 0
check_success view_out/references.bib


# --- 4. --geometry-file fallback (no .HSX) ---
echo -e "\n--- Testing the --geometry-file fallback (no .HSX) ---"

timeout 60 stb-ipr --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --geometry-file calc.fdf \
    --shift fermi -o approx_out --no-intro > log_approx.txt 2>&1
check_exit_code $? 0
check_contains "implicit-orthogonal-basis approximation" log_approx.txt


# --- 5. --q (custom IPR order + the new q>=2 validation) ---
echo -e "\n--- Testing --q 4 ---"

timeout 60 stb-ipr --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX \
    --shift fermi --q 4 -o q4_out --no-intro > log_q4.txt 2>&1
check_exit_code $? 0
check_contains "IPR order (q)  : 4" log_q4.txt

echo "Testing: --q 1 is a clean error, not a raw sisl AssertionError traceback"
stb-ipr --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX --shift fermi --q 1 --no-intro > log_q1.txt 2>&1
check_exit_code $? 2
check_contains "\-\-q must be >= 2" log_q1.txt
if ! grep -q "Traceback" log_q1.txt; then
    echo -e "   -> ${GREEN}Verified:${NC} no raw Python traceback leaked to the user"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} a raw traceback leaked instead of a clean error"
    FAIL=$((FAIL+1))
fi


# --- 6. --shift manual ---
echo -e "\n--- Testing --shift manual ---"

timeout 60 stb-ipr --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX \
    --shift manual --manual-value -3.2 -o manual_out --no-intro > log_manual.txt 2>&1
check_exit_code $? 0
check_contains "manual (-3.2 eV)" log_manual.txt


# --- 7. --label + --hsx-file (real bug: used to be rejected outright) ---
echo -e "\n--- Testing --label + --hsx-file together (a renamed .HSX, decoupled from --label) ---"

cp Sn3O4.HSX calc.HSX
timeout 60 stb-ipr --label Sn3O4 --hsx-file calc.HSX --shift fermi -o labelhsx_out --no-intro > log_labelhsx.txt 2>&1
check_exit_code $? 0
check_contains "Geometry/Hamiltonian source | calc.HSX" log_labelhsx.txt
rm -f calc.HSX

echo "Testing: --label + --wfsx together is still rejected (ambiguous)"
stb-ipr --label Sn3O4 --wfsx Sn3O4.bands.WFSX --shift fermi --no-intro > log_labelwfsx.txt 2>&1
check_exit_code $? 2
check_contains "cannot be combined with --file/--wfsx" log_labelwfsx.txt


# --- 8. nspin=2 fix: spin channels no longer silently merged ---
echo -e "\n--- Testing nspin=2 (real spin-polarized fixture): ipr_up/ipr_down written separately ---"

pushd spin > /dev/null
timeout 60 stb-ipr --label Ospin --shift fermi --save-gnuplot -o spin_out --no-intro > log_spin.txt 2>&1
check_exit_code $? 0
check_contains "Spin-resolved  : yes" log_spin.txt
check_success spin_out/ipr_up.dat
check_success spin_out/ipr_down.dat
if [ ! -f spin_out/ipr.dat ]; then
    echo -e "   -> ${GREEN}Verified:${NC} no merged ipr.dat written for a spin-polarized run"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} a merged ipr.dat was written instead of separate spin series"
    FAIL=$((FAIL+1))
fi

echo "Testing: spin-up and spin-down mean IPR are genuinely different (not silently averaged)"
python3 -c "
import numpy as np, sys
up = np.loadtxt('spin_out/ipr_up.dat')[:,2]
down = np.loadtxt('spin_out/ipr_down.dat')[:,2]
sys.exit(0 if abs(up.mean() - down.mean()) > 0.01 else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} spin-up/spin-down mean IPR genuinely differ"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} spin-up/spin-down mean IPR are suspiciously identical"
    FAIL=$((FAIL+1))
fi
popd > /dev/null


# --- 9. Robustness / errors ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: --file without --wfsx"
stb-ipr --file Sn3O4.bands --shift fermi --no-intro > log_missing_wfsx.txt 2>&1
check_exit_code $? 2

echo "Testing: missing --shift"
stb-ipr --label Sn3O4 --no-intro > log_missing_shift.txt 2>&1
check_exit_code $? 2

echo "Testing: --shift manual without --manual-value"
stb-ipr --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX --shift manual --no-intro > log_missing_manual.txt 2>&1
check_exit_code $? 2

echo "Testing: missing --label/--file"
stb-ipr --shift fermi --no-intro > log_missing_label.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-ipr --version > log_version.txt 2>&1
check_contains "stb-ipr" log_version.txt

echo "Testing: --help documents --q, --hsx-file, --geometry-file, --k-tol, --save-report, --save-gnuplot, --view"
stb-ipr --help > log_help.txt 2>&1
check_contains "\-\-q " log_help.txt
check_contains "hsx-file" log_help.txt
check_contains "geometry-file" log_help.txt
check_contains "k-tol" log_help.txt
check_contains "save-report" log_help.txt
check_contains "save-gnuplot" log_help.txt
check_contains "view" log_help.txt


# --- 10. Interactive path (stb-suite, shortcut 3.15) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 3.15) ---"

echo "Testing: navigate 3.15 -> label Sn3O4 -> default HSX -> shift fermi -> q default -> default output -> no report/gnuplot/view"
rm -f references.bib
printf '3.15\nSn3O4\n\n3\n\n\nn\nn\nn\n\n0\n' | timeout 60 stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success references.bib

echo "Testing: interactive menu's --save-gnuplot prompt"
rm -f ipr.dat ipr.gplot
printf '3.15\nSn3O4\n\n3\n\n\nn\ny\nn\n\n0\n' | timeout 60 stb-suite > log_menu_gnuplot.txt 2>&1
check_exit_code $? 0
check_success ipr.dat
check_success ipr.gplot


popd > /dev/null

# --- 11. Summary ---
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
