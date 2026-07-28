#!/bin/bash

# --- Setup ---
# Smoke test for stb-spintexture (Spin Texture Analyzer, item 3.17)
#
# Fixture note: IsolatedO.HSX/.selected.WFSX are a real SIESTA run (a
# single O atom in a large vacuum box, Spin non-collinear, Gamma-point
# only). WITHOUT an explicit initial spin canting (%block DM.InitSpin --
# its exact syntax could not be confirmed reliably during development),
# this converges with a physically sensible, strongly non-zero Sz (~+-1.0,
# expected for an open-shell atom) and Sx/Sy ~0 (numerical noise) -- this
# validates the numerical pipeline, but does not demonstrate a genuinely
# "textured" (canted) spin texture. See the module docstring for details.
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
echo "--- Starting tester for stb-spintexture (item 3.17) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR"/IsolatedO.selected.WFSX "$FIXTURE_DIR"/IsolatedO.HSX "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Basic run: --label auto-detect, no shift, numbered report ---
echo -e "\n--- Testing a basic run (--label IsolatedO), numbered report ---"

timeout 60 stb-spintexture --label IsolatedO --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "nspin=4 (non-collinear)" log_basic.txt
check_contains "Using 'IsolatedO.HSX'" log_basic.txt
check_contains "\[0\] RUN METADATA" log_basic.txt
check_contains "\[1\] INPUT DATA" log_basic.txt
check_contains "\[2\] ENERGY REFERENCE" log_basic.txt
check_contains "\[3\] SPIN TEXTURE ANALYSIS" log_basic.txt
check_contains "\[4\] OUTPUT DATA & PLOTS" log_basic.txt
check_contains "\[5\] REFERENCES" log_basic.txt
check_contains "\[6\] SUMMARY & FILES" log_basic.txt

echo "Testing: a plain run with no --view/--save-report/--save-gnuplot only writes references.bib"
if [ ! -f spintexture_Sx.dat ] && [ ! -f spintexture.gplot ] && [ ! -f stb_spintexture_report.txt ] && [ -f references.bib ]; then
    echo -e "   -> ${GREEN}Verified:${NC} no data/gplot/report written without --save-gnuplot/--save-report"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} data/gplot/report was written despite being off by default"
    FAIL=$((FAIL+1))
fi
rm -f references.bib

echo "Testing: Sz is strongly non-zero (~+-1.0), Sx/Sy are ~0 (documented fixture limitation)"
timeout 60 stb-spintexture --label IsolatedO --save-gnuplot --no-intro > log_basic_gnuplot.txt 2>&1
python3 -c "
import numpy as np, sys
sz = np.loadtxt('spintexture_Sz.dat')[:,2]
sx = np.loadtxt('spintexture_Sx.dat')[:,2]
sy = np.loadtxt('spintexture_Sy.dat')[:,2]
ok = (np.abs(sz) > 0.9).any() and (np.abs(sx) < 1e-6).all() and (np.abs(sy) < 1e-6).all()
sys.exit(0 if ok else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} Sz ~+-1.0, Sx/Sy ~0 (matches documented fixture behavior)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} spin moments did not match expected pattern"
    FAIL=$((FAIL+1))
fi

echo "Testing: the |S| <= 1 normalization check is reported and passes"
check_contains "stays within the physical bound of 1" log_basic_gnuplot.txt

echo "Testing: spintexture.gplot uses bare basenames (not the full --output-dir path)"
if grep -q '"spintexture_Sx.dat"' spintexture.gplot; then
    echo -e "   -> ${GREEN}Verified:${NC} spintexture.gplot references bare basenames"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} spintexture.gplot did not reference a bare basename"
    FAIL=$((FAIL+1))
fi

echo "Testing: the saved .gplot script actually renders with a real gnuplot"
if command -v gnuplot >/dev/null 2>&1; then
    gnuplot spintexture.gplot
    check_success spintexture.pdf
else
    echo -e "   -> ${YELLOW}Skipped:${NC} gnuplot not installed in this environment"
fi
rm -f spintexture_S*.dat spintexture.gplot spintexture.pdf references.bib


# --- 3. --save-report / --view ---
echo -e "\n--- Testing --save-report / --view ---"

timeout 60 stb-spintexture --label IsolatedO --save-report -o savereport_out --no-intro > log_savereport.txt 2>&1
check_exit_code $? 0
check_success savereport_out/stb_spintexture_report.txt

timeout 60 stb-spintexture --label IsolatedO --view -o view_out --no-intro > log_view.txt 2>&1
check_exit_code $? 0
check_success view_out/references.bib


# --- 4. --shift fermi/vbm/cbm: the Fermi-source hierarchy (same as stb-effmass) ---
echo -e "\n--- Testing --shift fermi (--fermi explicit value) ---"

timeout 60 stb-spintexture --label IsolatedO --shift fermi --fermi -12.0 -o fermi_out --no-intro > log_fermi.txt 2>&1
check_exit_code $? 0
check_success fermi_out/references.bib
check_contains "Fermi source: --fermi (explicit value)" log_fermi.txt

echo -e "\n--- Testing --shift fermi with --fermi-file (an arbitrarily-named .out log) ---"
cat > my_weird_name.out << 'EOF'
Some SIESTA log noise before the summary
siesta:         Fermi = -12.0
More noise after
EOF
timeout 60 stb-spintexture --label IsolatedO --shift fermi --fermi-file my_weird_name.out \
    -o fermifile_out --no-intro > log_fermifile.txt 2>&1
check_exit_code $? 0
check_contains "Fermi source: 'my_weird_name.out' (--fermi-file)" log_fermifile.txt

echo -e "\n--- Testing --shift vbm auto-detecting a generic .out ---"
timeout 60 stb-spintexture --label IsolatedO --shift vbm -o vbmauto_out --no-intro > log_vbmauto.txt 2>&1
check_exit_code $? 0
check_contains "auto-detected .out" log_vbmauto.txt
check_contains "VBM = " log_vbmauto.txt
rm -f my_weird_name.out

echo "Testing: --shift fermi with no Fermi source at all is a clean error"
stb-spintexture --label IsolatedO --shift fermi --no-intro > log_missing_fermi2.txt 2>&1
check_exit_code $? 2
check_contains "needs a Fermi energy" log_missing_fermi2.txt


# --- 5. --label + --hsx-file (real bug: used to be rejected outright) ---
echo -e "\n--- Testing --label + --hsx-file together (a renamed .HSX, decoupled from --label) ---"

cp IsolatedO.HSX calc.HSX
timeout 60 stb-spintexture --label IsolatedO --hsx-file calc.HSX -o labelhsx_out --no-intro > log_labelhsx.txt 2>&1
check_exit_code $? 0
check_contains "Hamiltonian src: calc.HSX" log_labelhsx.txt
rm -f calc.HSX

echo "Testing: --label + --wfsx together is still rejected (ambiguous)"
stb-spintexture --label IsolatedO --wfsx IsolatedO.selected.WFSX --no-intro > log_labelwfsx.txt 2>&1
check_exit_code $? 2
check_contains "cannot be combined with --wfsx" log_labelwfsx.txt


# --- 6. Explicit --wfsx/--hsx-file (no --label) ---
echo -e "\n--- Testing explicit --wfsx/--hsx-file (no --label) ---"

timeout 60 stb-spintexture --wfsx IsolatedO.selected.WFSX --hsx-file IsolatedO.HSX \
    -o explicit_out --no-intro > log_explicit.txt 2>&1
check_exit_code $? 0
check_success explicit_out/references.bib


# --- 7. nspin validation: reject a non-nc/soc WFSX ---
echo -e "\n--- Testing rejection of a collinear (nspin=1) WFSX ---"

BULK_WFSX="$FIXTURE_DIR/../14-coop/Sn3O4.selected.WFSX"
BULK_HSX="$FIXTURE_DIR/../14-coop/Sn3O4.HSX"
if [ -f "$BULK_WFSX" ] && [ -f "$BULK_HSX" ]; then
    stb-spintexture --wfsx "$BULK_WFSX" --hsx-file "$BULK_HSX" --no-intro > log_nspin.txt 2>&1
    check_exit_code $? 2
    check_contains "only makes physical sense for a non-collinear" log_nspin.txt
else
    echo -e " ... ${YELLOW}SKIPPED${NC} (test/3-analysis/14-coop fixture not found)"
fi


# --- 8. Robustness / errors ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: --shift manual without --manual-value"
stb-spintexture --label IsolatedO --shift manual --no-intro > log_missing_manual.txt 2>&1
check_exit_code $? 2

echo "Testing: missing --label/--wfsx"
stb-spintexture --no-intro > log_missing_label.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-spintexture --version > log_version.txt 2>&1
check_contains "stb-spintexture" log_version.txt

echo "Testing: --help documents --shift, --fermi, --hsx-file, --geometry-file, --save-report, --save-gnuplot, --view"
stb-spintexture --help > log_help.txt 2>&1
check_contains "shift" log_help.txt
check_contains "fermi" log_help.txt
check_contains "hsx-file" log_help.txt
check_contains "geometry-file" log_help.txt
check_contains "save-report" log_help.txt
check_contains "save-gnuplot" log_help.txt
check_contains "view" log_help.txt


# --- 9. Interactive path (stb-suite, shortcut 3.17) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 3.17) ---"

echo "Testing: navigate 3.17 -> label IsolatedO -> default HSX -> shift 1 (none) -> default output -> no report/gnuplot/view"
rm -f references.bib
printf '3.17\nIsolatedO\n\n1\n\nn\nn\nn\n\n0\n' | timeout 60 stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success references.bib

echo "Testing: interactive menu's Fermi-source submenu (shift 4=fermi, option 1=explicit)"
rm -f references.bib
printf '3.17\nIsolatedO\n\n4\n1\n-12.0\n\nn\nn\nn\n\n0\n' | timeout 60 stb-suite > log_menu_fermi.txt 2>&1
check_exit_code $? 0
check_contains "Fermi source: --fermi (explicit value)" log_menu_fermi.txt


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
