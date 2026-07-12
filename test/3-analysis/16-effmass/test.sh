#!/bin/bash

# --- Setup ---
# Smoke test for stb-effmass (Effective Mass / Velocity, item 3.16)
#
# Reuses test/3-analysis/14-coop/'s Sn3O4 fixture (Sn3O4.HSX +
# Sn3O4.selected.WFSX, a real full-BZ 2x2x2-equivalent k-mesh run) -- no
# new SIESTA run needed, since effective_mass()/velocity() just need a
# real Hamiltonian + WFSX coefficients, same input stb-coop already uses.
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


# --- 2. Basic run: --k-index 0 (Gamma) --band 1 ---
echo -e "\n--- Testing a basic run (--k-index 0 --band 1, at Gamma) ---"

timeout 60 stb-effmass --wfsx Sn3O4.selected.WFSX --hsx-file Sn3O4.HSX \
    --k-index 0 --band 1 --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "Effective mass tensor" log_basic.txt
check_contains "Band velocity" log_basic.txt
check_success effmass.txt

echo "Testing: velocity is ~0 at Gamma (time-reversal-invariant point, physically expected)"
python3 -c "
import re, sys
text = open('effmass.txt').read()
vals = [float(x) for x in re.findall(r'v_[xyz] = *(-?[\d.]+) km/s', text)]
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
text = open('effmass.txt').read()
m = re.search(r'm\*_xx *= *(-?[\d.]+) m0', text)
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


# --- 3. --band vbm: expect NEGATIVE effective mass (valence band curvature) ---
echo -e "\n--- Testing --band vbm (expect negative effective mass, hole-like) ---"

timeout 60 stb-effmass --wfsx Sn3O4.selected.WFSX --hsx-file Sn3O4.HSX \
    --band vbm --fermi -3.200055 -o vbm_out --no-intro > log_vbm.txt 2>&1
check_exit_code $? 0
check_contains "VBM found at k-index" log_vbm.txt
echo "Testing: at least one Voigt component is negative at the VBM"
python3 -c "
import re, sys
text = open('vbm_out/effmass.txt').read()
vals = [float(x) for x in re.findall(r'm\*_\w+ *= *(-?[\d.]+) m0', text)]
sys.exit(0 if any(v < 0 for v in vals) else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} negative effective mass found at VBM (physically expected)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} no negative effective mass at VBM"
    FAIL=$((FAIL+1))
fi


# --- 4. --k-point (match by vector) ---
echo -e "\n--- Testing --k-point ---"

timeout 60 stb-effmass --wfsx Sn3O4.selected.WFSX --hsx-file Sn3O4.HSX \
    --k-point 0 0 0 --band 2 -o kpoint_out --no-intro > log_kpoint.txt 2>&1
check_exit_code $? 0
check_success kpoint_out/effmass.txt


# --- 5. Non-collinear/SOC Hamiltonian: effective mass degrades gracefully
#         to velocity-only (no crash) -- regression test for a real bug
#         found and fixed during manual review: sisl's effective_mass()
#         raises a bare TypeError for nspin=4/8 (ddSk not implemented
#         there), which this tool now detects and reports cleanly instead.
echo -e "\n--- Testing graceful degradation on a non-collinear (nspin=4) Hamiltonian ---"

NC_WFSX="$FIXTURE_DIR/../17-spintexture/IsolatedO.selected.WFSX"
NC_HSX="$FIXTURE_DIR/../17-spintexture/IsolatedO.HSX"
if [ -f "$NC_WFSX" ] && [ -f "$NC_HSX" ]; then
    timeout 60 stb-effmass --wfsx "$NC_WFSX" --hsx-file "$NC_HSX" --k-index 0 --band 3 \
        -o nc_out --no-intro > log_nc.txt 2>&1
    check_exit_code $? 0
    check_contains "not supported by sisl for non-collinear/spin-orbit" log_nc.txt
    check_contains "Effective mass tensor: N/A" nc_out/effmass.txt
    check_contains "Band velocity" nc_out/effmass.txt
else
    echo -e " ... ${YELLOW}SKIPPED${NC} (test/3-analysis/17-spintexture fixture not found)"
fi


# --- 6. Robustness / errors ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: --wfsx without --hsx-file"
stb-effmass --wfsx Sn3O4.selected.WFSX --band 1 --no-intro > log_missing_hsx.txt 2>&1
check_exit_code $? 2

echo "Testing: --band vbm without --fermi"
stb-effmass --wfsx Sn3O4.selected.WFSX --hsx-file Sn3O4.HSX --band vbm --no-intro > log_missing_fermi.txt 2>&1
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

echo "Testing: --help documents --k-index, --k-point, --band, --fermi, --hsx-file"
stb-effmass --help > log_help.txt 2>&1
check_contains "k-index" log_help.txt
check_contains "k-point" log_help.txt
check_contains "band" log_help.txt
check_contains "fermi" log_help.txt
check_contains "hsx-file" log_help.txt


# --- 7. Interactive path (stb-suite, shortcut 3.16) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 3.16) ---"

echo "Testing: navigate 3.16 -> label Sn3O4 -> band selection 1 -> k-index 0 -> band 1 -> default output"
rm -f effmass.txt
printf '3.16\nSn3O4\n1\n0\n1\n\n\n0\n' | timeout 60 stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success effmass.txt


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
