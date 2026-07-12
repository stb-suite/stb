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
# Same-name fdf copy for the interactive-menu test (no --geometry-file prompt there).
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


# --- 5. --no-plot ---
echo -e "\n--- Testing --no-plot ---"

timeout 60 stb-sts --wfsx Graphene.selected.WFSX --geometry-file calc.fdf \
    --xy 0 0 --height 1.5 --erange -3 3 --sigma 100 --no-plot -o noplot_out --no-intro > log_noplot.txt 2>&1
check_exit_code $? 0
check_success noplot_out/sts.dat


# --- 6. Robustness / errors ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: --xy without --height"
stb-sts --wfsx Graphene.selected.WFSX --geometry-file calc.fdf --xy 0 0 --erange -3 3 --sigma 50 --no-intro > log_missing_height.txt 2>&1
check_exit_code $? 2

echo "Testing: --shift fermi without --fermi"
stb-sts --wfsx Graphene.selected.WFSX --geometry-file calc.fdf --xy 0 0 --height 1.5 --erange -3 3 --sigma 50 --shift fermi --no-intro > log_missing_fermi.txt 2>&1
check_exit_code $? 2

echo "Testing: missing --wfsx/--label"
stb-sts --xy 0 0 --height 1.5 --erange -3 3 --sigma 50 --no-intro > log_missing_label.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-sts --version > log_version.txt 2>&1
check_contains "stb-sts" log_version.txt

echo "Testing: --help documents --xy, --point, --height, --erange, --sigma, --fwhm, --shift"
stb-sts --help > log_help.txt 2>&1
check_contains "xy" log_help.txt
check_contains "point" log_help.txt
check_contains "height" log_help.txt
check_contains "erange" log_help.txt
check_contains "sigma" log_help.txt
check_contains "fwhm" log_help.txt
check_contains "shift" log_help.txt


# --- 7. Interactive path (stb-suite, shortcut 3.13) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 3.13) ---"

echo "Testing: navigate 3.13 -> label Graphene -> tip mode 1 (xy+height) -> 0 0 1.5 -> erange -3 3 -> sigma 100 -> default output"
rm -f sts.dat
printf '3.13\nGraphene\n1\n0\n0\n1.5\n-3\n3\n100\n\n\n0\n' | timeout 60 stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success sts.dat


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
