#!/bin/bash

# --- Setup ---
# Smoke test for stb-wfdensity (Wavefunction Density, item 3.12)
#
# Reuses test/3-analysis/10-fatbands/'s Sn3O4 fixture (calc.fdf + .ion files
# + Sn3O4.bands.WFSX) -- no new SIESTA run needed, since a band-path WFSX is
# perfectly valid input for picking one (k-index, band) state (unlike
# stb-sts/stb-coop, which need a real full-BZ mesh).
export MPLBACKEND=Agg
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="$FIXTURE_DIR/../10-fatbands"
TEST_DIR="$FIXTURE_DIR/test_files"

# Output colors
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

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
echo "--- Starting tester for stb-wfdensity (item 3.12) ---"
if [ ! -f "$SOURCE_DIR/Sn3O4.bands.WFSX" ]; then
    echo -e "${RED}FATAL: $SOURCE_DIR/Sn3O4.bands.WFSX not found -- run test/3-analysis/10-fatbands/test.sh's fixture setup first.${NC}"
    exit 1
fi
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$SOURCE_DIR"/Sn3O4.bands.WFSX "$SOURCE_DIR"/calc.fdf "$SOURCE_DIR"/structure.fdf \
   "$SOURCE_DIR"/Sn.ion "$SOURCE_DIR"/Sn.ion.xml "$SOURCE_DIR"/O.ion "$SOURCE_DIR"/O.ion.xml "$TEST_DIR/"
# The interactive stb-suite menu (unlike the CLI) has no --geometry-file
# prompt -- it only auto-detects <label>.fdf by exact name. The real fdf
# here is named calc.fdf (a common real-world mismatch with SystemLabel),
# so give the interactive-path test (section 8) a same-named copy too.
cp "$SOURCE_DIR"/calc.fdf "$TEST_DIR/Sn3O4.fdf"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Basic run: --wfsx/--geometry-file explicit, --k-index 0 --band 1 ---
echo -e "\n--- Testing a basic run (--k-index 0 --band 1) ---"

timeout 60 stb-wfdensity --wfsx Sn3O4.bands.WFSX --geometry-file calc.fdf \
    --k-index 0 --band 1 --spacing 0.4 --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "Normalization check" log_basic.txt
check_success wfdensity_k0_b1.cube
check_success wfdensity_k0_b1_slice.dat


# --- 3. --cube-only skips the slice/profile output ---
echo -e "\n--- Testing --cube-only ---"

timeout 60 stb-wfdensity --wfsx Sn3O4.bands.WFSX --geometry-file calc.fdf \
    --k-index 0 --band 2 --spacing 0.5 --cube-only -o cubeonly_out --no-intro > log_cubeonly.txt 2>&1
check_exit_code $? 0
check_success cubeonly_out/wfdensity_k0_b2.cube
if [ -e cubeonly_out/wfdensity_k0_b2_slice.dat ]; then
    echo -e " ... ${RED}FAIL${NC} (slice file should not exist with --cube-only)"
    FAIL=$((FAIL+1))
else
    echo -e " ... ${GREEN}OK${NC} (slice file correctly absent)"
    PASS=$((PASS+1))
fi


# --- 4. --profile ---
echo -e "\n--- Testing --profile ---"

timeout 60 stb-wfdensity --wfsx Sn3O4.bands.WFSX --geometry-file calc.fdf \
    --k-index 0 --band 1 --spacing 0.5 --profile -o profile_out --no-intro > log_profile.txt 2>&1
check_exit_code $? 0
check_success profile_out/wfdensity_k0_b1_profile.dat


# --- 5. --band vbm/cbm with --fermi ---
echo -e "\n--- Testing --band vbm ---"

timeout 60 stb-wfdensity --wfsx Sn3O4.bands.WFSX --geometry-file calc.fdf \
    --band vbm --fermi -3.2038563052897451 --spacing 0.6 --cube-only -o vbm_out --no-intro > log_vbm.txt 2>&1
check_exit_code $? 0
check_contains "VBM found at k-index" log_vbm.txt


# --- 6. --k-point matches by vector instead of index ---
echo -e "\n--- Testing --k-point ---"

timeout 60 stb-wfdensity --wfsx Sn3O4.bands.WFSX --geometry-file calc.fdf \
    --k-point 0 0 0 --band 2 --spacing 0.6 --cube-only -o kpoint_out --no-intro > log_kpoint.txt 2>&1
check_exit_code $? 0
check_success kpoint_out/wfdensity_k0_b2.cube


# --- 7. Errors ---
echo -e "\n--- Testing error cases ---"

echo "Testing: --wfsx without --geometry-file"
stb-wfdensity --wfsx Sn3O4.bands.WFSX --band 1 --no-intro > log_missing_geo.txt 2>&1
check_exit_code $? 2

echo "Testing: invalid --band"
stb-wfdensity --wfsx Sn3O4.bands.WFSX --geometry-file calc.fdf --band notanumber --no-intro > log_bad_band.txt 2>&1
check_exit_code $? 2

echo "Testing: --k-index out of range"
stb-wfdensity --wfsx Sn3O4.bands.WFSX --geometry-file calc.fdf --k-index 999 --band 1 --no-intro > log_bad_k.txt 2>&1
check_exit_code $? 2
check_contains "out of range" log_bad_k.txt

echo "Testing: missing --label/--wfsx"
stb-wfdensity --band 1 --no-intro > log_missing_label.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-wfdensity --version > log_version.txt 2>&1
check_contains "stb-wfdensity" log_version.txt

echo "Testing: --help documents --k-index, --k-point, --band, --spacing, --cube-only, --geometry-file"
stb-wfdensity --help > log_help.txt 2>&1
check_contains "k-index" log_help.txt
check_contains "k-point" log_help.txt
check_contains "band" log_help.txt
check_contains "spacing" log_help.txt
check_contains "cube-only" log_help.txt
check_contains "geometry-file" log_help.txt


# --- 8. Interactive path (stb-suite, shortcut 3.12) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 3.12) ---"

echo "Testing: navigate 3.12 -> label Sn3O4 -> band selection 1 -> k-index 0 -> band 1 -> spacing 0.5 -> default output"
rm -f wfdensity_k0_b1.cube
printf '3.12\nSn3O4\n1\n0\n1\n0.5\n\n\n0\n' | timeout 60 stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success wfdensity_k0_b1.cube


popd > /dev/null

# --- 9. Summary ---
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
