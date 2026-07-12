#!/bin/bash

# --- Setup ---
# Smoke test for stb-ipr (IPR Analyzer, item 3.15)
#
# Reuses test/3-analysis/10-fatbands/'s Sn3O4 fixture (calc.fdf + .ion files
# + Sn3O4.bands.WFSX + Sn3O4.HSX) -- IPR only needs orbital expansion
# coefficients (norm2/hadamard), the exact same input stb-fatbands already
# uses, so no new SIESTA run is needed.
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
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$SOURCE_DIR"/Sn3O4.bands "$SOURCE_DIR"/Sn3O4.bands.WFSX "$SOURCE_DIR"/Sn3O4.HSX \
   "$SOURCE_DIR"/calc.fdf "$SOURCE_DIR"/structure.fdf "$SOURCE_DIR"/Sn.ion "$SOURCE_DIR"/Sn.ion.xml \
   "$SOURCE_DIR"/O.ion "$SOURCE_DIR"/O.ion.xml "$TEST_DIR/"
# Same-name fdf copy for the interactive-menu test (no --geometry-file
# prompt there, only auto-detected <label>.HSX/.fdf by exact name).
cp "$SOURCE_DIR"/calc.fdf "$TEST_DIR/Sn3O4.fdf"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Basic run: --file/--wfsx/--hsx-file explicit (accurate, overlap-aware) ---
echo -e "\n--- Testing a basic run (--hsx-file, accurate path) ---"

timeout 60 stb-ipr --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX \
    --shift fermi --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "Using 'Sn3O4.HSX' for overlap-aware IPR" log_basic.txt
check_contains "correspondence check passed" log_basic.txt
check_success ipr.dat
check_success ipr.gplot

echo "Testing: deep (core-like) states have higher IPR (more localized) than near-Fermi states"
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


# --- 3. --geometry-file fallback (no .HSX, approximate) ---
echo -e "\n--- Testing the --geometry-file fallback (no .HSX) ---"

timeout 60 stb-ipr --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --geometry-file calc.fdf \
    --shift fermi -o approx_out --no-intro > log_approx.txt 2>&1
check_exit_code $? 0
check_contains "implicit-orthogonal-basis approximation" log_approx.txt
check_success approx_out/ipr.dat


# --- 4. --q (custom IPR order) ---
echo -e "\n--- Testing --q 4 ---"

timeout 60 stb-ipr --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX \
    --shift fermi --q 4 -o q4_out --no-intro > log_q4.txt 2>&1
check_exit_code $? 0
check_success q4_out/ipr.dat


# --- 5. --shift vbm/manual ---
echo -e "\n--- Testing --shift manual ---"

timeout 60 stb-ipr --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX \
    --shift manual --manual-value -3.2 -o manual_out --no-intro > log_manual.txt 2>&1
check_exit_code $? 0
check_success manual_out/ipr.dat


# --- 6. Robustness / errors ---
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

echo "Testing: --help documents --q, --hsx-file, --geometry-file, --k-tol"
stb-ipr --help > log_help.txt 2>&1
check_contains "\-\-q " log_help.txt
check_contains "hsx-file" log_help.txt
check_contains "geometry-file" log_help.txt
check_contains "k-tol" log_help.txt


# --- 7. Interactive path (stb-suite, shortcut 3.15) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 3.15) ---"

echo "Testing: navigate 3.15 -> label Sn3O4 -> shift fermi -> q default -> default output"
rm -f ipr.dat ipr.gplot
printf '3.15\nSn3O4\n3\n\n\n\n0\n' | timeout 60 stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success ipr.dat


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
