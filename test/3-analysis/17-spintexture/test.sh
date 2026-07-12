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


# --- 2. Basic run: --label auto-detect, no shift ---
echo -e "\n--- Testing a basic run (--label IsolatedO) ---"

timeout 60 stb-spintexture --label IsolatedO --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "nspin=4 (non-collinear)" log_basic.txt
check_contains "Using 'IsolatedO.HSX'" log_basic.txt
check_success spintexture_Sx.dat
check_success spintexture_Sy.dat
check_success spintexture_Sz.dat
check_success spintexture.gplot

echo "Testing: Sz is strongly non-zero (~+-1.0), Sx/Sy are ~0 (documented fixture limitation)"
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


# --- 3. --shift fermi ---
echo -e "\n--- Testing --shift fermi ---"

timeout 60 stb-spintexture --label IsolatedO --shift fermi --fermi -12.0 -o fermi_out --no-intro > log_fermi.txt 2>&1
check_exit_code $? 0
check_success fermi_out/spintexture_Sz.dat


# --- 4. Explicit --wfsx/--hsx-file (no --label) ---
echo -e "\n--- Testing explicit --wfsx/--hsx-file (no --label) ---"

timeout 60 stb-spintexture --wfsx IsolatedO.selected.WFSX --hsx-file IsolatedO.HSX \
    -o explicit_out --no-intro > log_explicit.txt 2>&1
check_exit_code $? 0
check_success explicit_out/spintexture_Sz.dat


# --- 5. nspin validation: reject a non-nc/soc WFSX ---
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


# --- 6. Robustness / errors ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: --shift manual without --manual-value"
stb-spintexture --label IsolatedO --shift manual --no-intro > log_missing_manual.txt 2>&1
check_exit_code $? 2

echo "Testing: --shift fermi without --fermi"
stb-spintexture --label IsolatedO --shift fermi --no-intro > log_missing_fermi.txt 2>&1
check_exit_code $? 2

echo "Testing: missing --label/--wfsx"
stb-spintexture --no-intro > log_missing_label.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-spintexture --version > log_version.txt 2>&1
check_contains "stb-spintexture" log_version.txt

echo "Testing: --help documents --shift, --fermi, --hsx-file, --geometry-file"
stb-spintexture --help > log_help.txt 2>&1
check_contains "shift" log_help.txt
check_contains "fermi" log_help.txt
check_contains "hsx-file" log_help.txt
check_contains "geometry-file" log_help.txt


# --- 7. Interactive path (stb-suite, shortcut 3.17) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 3.17) ---"

echo "Testing: navigate 3.17 -> label IsolatedO -> shift 1 (none) -> default output"
rm -f spintexture_Sx.dat spintexture_Sy.dat spintexture_Sz.dat
printf '3.17\nIsolatedO\n1\n\n\n0\n' | timeout 60 stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success spintexture_Sz.dat


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
