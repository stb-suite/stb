#!/bin/bash

# --- Setup ---
# Smoke test for stb-effmass (Effective Mass / Velocity, item 3.16)
#
# Reuses test/3-analysis/14-coop/'s Sn3O4 fixture (Sn3O4.HSX +
# Sn3O4.selected.WFSX, a real full-BZ 2x2x2-equivalent k-mesh run) -- no
# new SIESTA run needed, since effective_mass()/velocity() just need a
# real Hamiltonian + WFSX coefficients, same input stb-coop already uses.
# Also reuses test/3-analysis/10-fatbands/spin/'s real spin-polarized
# isolated-O-atom fixture (Ospin.*) for the nspin=2 regression (the exact
# bug a real user hit: sisl's ddPk() has no 'spin' parameter at all), and
# test/3-analysis/17-spintexture/'s non-collinear fixture for nspin=4.
export MPLBACKEND=Agg
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="$FIXTURE_DIR/../14-coop"
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
echo "--- Starting tester for stb-effmass (item 3.16) ---"
if [ ! -f "$SOURCE_DIR/Sn3O4.selected.WFSX" ]; then
    echo -e "${RED}FATAL: $SOURCE_DIR/Sn3O4.selected.WFSX not found -- run test/3-analysis/14-coop/test.sh's fixture setup first.${NC}"
    exit 1
fi
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$SOURCE_DIR"/Sn3O4.selected.WFSX "$SOURCE_DIR"/Sn3O4.HSX "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Basic run: --k-index 0 (Gamma) --band 1, numbered report ---
echo -e "\n--- Testing a basic run (--k-index 0 --band 1, at Gamma), numbered report ---"

timeout 60 stb-effmass --wfsx Sn3O4.selected.WFSX --hsx-file Sn3O4.HSX \
    --k-index 0 --band 1 --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "\[0\] RUN METADATA" log_basic.txt
check_contains "\[1\] INPUT DATA" log_basic.txt
check_contains "\[2\] STATE SELECTION" log_basic.txt
check_contains "\[3\] EFFECTIVE MASS (per-axis Voigt" log_basic.txt
check_contains "\[4\] EFFECTIVE MASS (principal, full tensor)" log_basic.txt
check_contains "\[5\] BAND VELOCITY" log_basic.txt
check_contains "\[6\] OUTPUT DATA & PLOTS" log_basic.txt
check_contains "\[7\] REFERENCES" log_basic.txt
check_contains "\[8\] SUMMARY & FILES" log_basic.txt

echo "Testing: velocity is ~0 at Gamma (time-reversal-invariant point, physically expected)"
python3 -c "
import re, sys
text = open('log_basic.txt').read()
vals = [float(x) for x in re.findall(r'v_[xyz] *\| *(-?[\d.]+)', text)]
sys.exit(0 if all(abs(v) < 0.01 for v in vals) else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} velocity ~0 at Gamma"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} velocity not ~0 at Gamma"
    FAIL=$((FAIL+1))
fi

echo "Testing: xx effective mass matches the independently-verified reference value (0.4445 m0)"
python3 -c "
import re, sys
text = open('log_basic.txt').read()
m = re.search(r'm\*_xx *\| *(-?[\d.]+)', text)
val = float(m.group(1))
sys.exit(0 if abs(val - 0.4445) < 0.001 else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} m*_xx matches expected reference value"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} m*_xx did not match expected reference value"
    FAIL=$((FAIL+1))
fi

echo "Testing: a plain run with no --view/--save-report/--save-gnuplot only writes references.bib"
if [ ! -f effmass.dat ] && [ ! -f velocity.dat ] && [ ! -f effmass.gplot ] && [ ! -f stb_effmass_report.txt ] && [ -f references.bib ]; then
    echo -e "   -> ${GREEN}Verified:${NC} no data/gplot/report written without --save-gnuplot/--save-report"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} data/gplot/report was written despite being off by default"
    FAIL=$((FAIL+1))
fi
rm -f references.bib


# --- 3. --band vbm: negative effective mass, and the real off-diagonal/principal-mass finding ---
echo -e "\n--- Testing --band vbm (expect negative effective mass, hole-like; principal-mass diagonalization) ---"

timeout 60 stb-effmass --wfsx Sn3O4.selected.WFSX --hsx-file Sn3O4.HSX \
    --band vbm --fermi -3.200055 -o vbm_out --no-intro --save-report --save-gnuplot > log_vbm.txt 2>&1
check_exit_code $? 0
check_contains "VBM found at k-index" log_vbm.txt
check_contains "Fermi source: --fermi (explicit value)" log_vbm.txt

echo "Testing: at least one per-axis Voigt component is negative at the VBM"
python3 -c "
import re, sys
text = open('vbm_out/stb_effmass_report.txt').read()
section = text.split('[3] EFFECTIVE MASS')[1].split('[4] EFFECTIVE MASS')[0]
vals = [float(x) for x in re.findall(r'm\*_\w+ *\| *(-?[\d.]+)', section)]
sys.exit(0 if any(v < 0 for v in vals) else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} negative per-axis effective mass found at VBM (physically expected)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} no negative per-axis effective mass at VBM"
    FAIL=$((FAIL+1))
fi

echo "Testing: the principal (diagonalized) effective masses differ from sisl's own naive per-axis values"
echo "(a real, verified physics finding: off-diagonal curvature is substantial for this state)"
check_contains "significant fraction of (or exceeds) the diagonal curvature" vbm_out/stb_effmass_report.txt
python3 -c "
import re, sys
text = open('vbm_out/stb_effmass_report.txt').read()
section = text.split('[4] EFFECTIVE MASS')[1].split('[5] BAND VELOCITY')[0]
principal = sorted(float(x) for x in re.findall(r'(-?[\d.]+) *\|', section))
naive_diag = sorted(float(x) for x in re.findall(r'm\*_(?:xx|yy|zz) *\| *(-?[\d.]+)', text))
diffs = [abs(p - n) for p, n in zip(principal, naive_diag)]
sys.exit(0 if max(diffs) > 0.1 else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} principal effective masses meaningfully differ from the naive per-axis diagonal values"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} principal masses unexpectedly matched the naive diagonal values"
    FAIL=$((FAIL+1))
fi

check_success vbm_out/effmass.dat
check_success vbm_out/velocity.dat
check_success vbm_out/effmass.gplot

echo "Testing: effmass.gplot uses bare basenames (not the full --output-dir path)"
if grep -q '"effmass.dat"' vbm_out/effmass.gplot && ! grep -q "vbm_out" vbm_out/effmass.gplot; then
    echo -e "   -> ${GREEN}Verified:${NC} effmass.gplot references bare basenames"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} effmass.gplot leaked the --output-dir path"
    FAIL=$((FAIL+1))
fi

echo "Testing: the saved .gplot script actually renders with a real gnuplot"
if command -v gnuplot >/dev/null 2>&1; then
    (cd vbm_out && gnuplot effmass.gplot)
    check_success vbm_out/effmass.pdf
else
    echo -e "   -> ${YELLOW}Skipped:${NC} gnuplot not installed in this environment"
fi


# --- 4. --view ---
echo -e "\n--- Testing --view (opt-in matplotlib preview) ---"

timeout 60 stb-effmass --wfsx Sn3O4.selected.WFSX --hsx-file Sn3O4.HSX \
    --band vbm --fermi -3.200055 --view -o view_out --no-intro > log_view.txt 2>&1
check_exit_code $? 0
check_success view_out/references.bib


# --- 5. --k-point (match by vector) ---
echo -e "\n--- Testing --k-point ---"

timeout 60 stb-effmass --wfsx Sn3O4.selected.WFSX --hsx-file Sn3O4.HSX \
    --k-point 0 0 0 --band 2 -o kpoint_out --no-intro > log_kpoint.txt 2>&1
check_exit_code $? 0
check_success kpoint_out/references.bib


# --- 6. --label + --hsx-file (real bug: used to be rejected outright) ---
echo -e "\n--- Testing --label + --hsx-file together (a renamed .HSX, decoupled from --label) ---"

cp Sn3O4.HSX calc.HSX
timeout 60 stb-effmass --label Sn3O4 --hsx-file calc.HSX --k-index 0 --band 1 \
    -o labelhsx_out --no-intro > log_labelhsx.txt 2>&1
check_exit_code $? 0
check_contains "Hamiltonian src: calc.HSX" log_labelhsx.txt
rm -f calc.HSX

echo "Testing: --label + --wfsx together is still rejected (ambiguous)"
stb-effmass --label Sn3O4 --wfsx Sn3O4.selected.WFSX --band 1 --no-intro > log_labelwfsx.txt 2>&1
check_exit_code $? 2
check_contains "cannot be combined with --wfsx" log_labelwfsx.txt


# --- 7. --shift fermi's Fermi-source hierarchy: --fermi-file, auto-detected .out ---
echo -e "\n--- Testing --band vbm with --fermi-file (an arbitrarily-named .out log) ---"

cat > my_weird_name.out << 'EOF'
Some SIESTA log noise before the summary
siesta:         Fermi = -3.200055
More noise after
EOF
timeout 60 stb-effmass --wfsx Sn3O4.selected.WFSX --hsx-file Sn3O4.HSX --band vbm \
    --fermi-file my_weird_name.out -o fermifile_out --no-intro > log_fermifile.txt 2>&1
check_exit_code $? 0
check_contains "Fermi source: 'my_weird_name.out' (--fermi-file)" log_fermifile.txt

echo -e "\n--- Testing --band vbm auto-detecting a generic .out ---"
timeout 60 stb-effmass --wfsx Sn3O4.selected.WFSX --hsx-file Sn3O4.HSX --band vbm \
    -o fermiauto_out --no-intro > log_fermiauto.txt 2>&1
check_exit_code $? 0
check_contains "auto-detected .out" log_fermiauto.txt
rm -f my_weird_name.out

echo "Testing: --band vbm/cbm with no Fermi source at all is a clean error"
stb-effmass --wfsx Sn3O4.selected.WFSX --hsx-file Sn3O4.HSX --band vbm --no-intro > log_missing_fermi2.txt 2>&1
check_exit_code $? 2
check_contains "needs a Fermi energy" log_missing_fermi2.txt


# --- 8. THE REAL BUG FIX: nspin=2 (real spin-polarized) no longer crashes ---
echo -e "\n--- Testing graceful degradation on a REAL spin-polarized (nspin=2) Hamiltonian ---"
echo "(regression test for a real user-reported crash: sisl's ddPk() has no 'spin' parameter"
echo " at all in this version, so effective_mass() used to crash with a raw TypeError for"
echo " ANY spin-polarized calculation, not just non-collinear/SOC as previously assumed)"

SPIN_WFSX="$FIXTURE_DIR/../10-fatbands/spin/Ospin.bands.WFSX"
SPIN_HSX="$FIXTURE_DIR/../10-fatbands/spin/Ospin.HSX"
if [ -f "$SPIN_WFSX" ] && [ -f "$SPIN_HSX" ]; then
    timeout 60 stb-effmass --wfsx "$SPIN_WFSX" --hsx-file "$SPIN_HSX" --k-index 0 --band 3 \
        -o spin_out --no-intro > log_spin.txt 2>&1
    check_exit_code $? 0
    check_contains "not supported by sisl for spin-resolved Hamiltonians (nspin=2)" log_spin.txt
    if ! grep -q "Traceback" log_spin.txt; then
        echo -e "   -> ${GREEN}Verified:${NC} no raw Python traceback leaked (the original bug report)"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} a raw traceback leaked instead of a clean warning"
        FAIL=$((FAIL+1))
    fi
    check_contains "BAND VELOCITY" log_spin.txt
else
    echo -e " ... ${YELLOW}SKIPPED${NC} (test/3-analysis/10-fatbands/spin fixture not found)"
fi


# --- 9. Non-collinear/SOC Hamiltonian: still degrades gracefully (pre-existing case) ---
echo -e "\n--- Testing graceful degradation on a non-collinear (nspin=4) Hamiltonian ---"

NC_WFSX="$FIXTURE_DIR/../17-spintexture/IsolatedO.selected.WFSX"
NC_HSX="$FIXTURE_DIR/../17-spintexture/IsolatedO.HSX"
if [ -f "$NC_WFSX" ] && [ -f "$NC_HSX" ]; then
    timeout 60 stb-effmass --wfsx "$NC_WFSX" --hsx-file "$NC_HSX" --k-index 0 --band 3 \
        -o nc_out --no-intro > log_nc.txt 2>&1
    check_exit_code $? 0
    check_contains "not supported by sisl for spin-resolved Hamiltonians (nspin=4)" log_nc.txt
    check_contains "BAND VELOCITY" log_nc.txt
else
    echo -e " ... ${YELLOW}SKIPPED${NC} (test/3-analysis/17-spintexture fixture not found)"
fi


# --- 10. Robustness / errors ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: --wfsx without --hsx-file"
stb-effmass --wfsx Sn3O4.selected.WFSX --band 1 --no-intro > log_missing_hsx.txt 2>&1
check_exit_code $? 2

echo "Testing: invalid --band"
stb-effmass --wfsx Sn3O4.selected.WFSX --hsx-file Sn3O4.HSX --band notanumber --no-intro > log_bad_band.txt 2>&1
check_exit_code $? 2

echo "Testing: missing --label/--wfsx"
stb-effmass --band 1 --no-intro > log_missing_label.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-effmass --version > log_version.txt 2>&1
check_contains "stb-effmass" log_version.txt

echo "Testing: --help documents --k-index, --k-point, --band, --fermi, --hsx-file, --save-report, --save-gnuplot, --view"
stb-effmass --help > log_help.txt 2>&1
check_contains "k-index" log_help.txt
check_contains "k-point" log_help.txt
check_contains "band" log_help.txt
check_contains "fermi" log_help.txt
check_contains "hsx-file" log_help.txt
check_contains "save-report" log_help.txt
check_contains "save-gnuplot" log_help.txt
check_contains "view" log_help.txt


# --- 11. Interactive path (stb-suite, shortcut 3.16) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 3.16) ---"

echo "Testing: navigate 3.16 -> label Sn3O4 -> default HSX -> band selection 1 -> k-index 0 -> band 1 -> default output -> no report/gnuplot/view"
rm -f references.bib
printf '3.16\nSn3O4\n\n1\n0\n1\n\nn\nn\nn\n\n0\n' | timeout 60 stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success references.bib

echo "Testing: interactive menu's VBM + Fermi-source submenu (option 1, explicit value)"
rm -f references.bib
printf '3.16\nSn3O4\n\n2\n1\n-3.200055\n\nn\nn\nn\n\n0\n' | timeout 60 stb-suite > log_menu_vbm.txt 2>&1
check_exit_code $? 0
check_contains "VBM found at k-index" log_menu_vbm.txt


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
