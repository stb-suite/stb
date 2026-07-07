#!/bin/bash

# --- Setup ---
# Smoke test for stb-bands (Bands Analyzer, item 2.1)
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$FIXTURE_DIR/test_files"

# Non-interactive matplotlib backend: main() calls plt.show(), which would
# otherwise block waiting for a window to close.
export MPLBACKEND=Agg

# Output colors
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

PASS=0
FAIL=0

# --- Check helpers ---

# Checks that file $1 exists and is not empty
check_success() {
    if [ -s "$1" ]; then
        echo -e " ... ${GREEN}OK${NC} (file '$1' created)"
        PASS=$((PASS+1))
    else
        echo -e " ... ${RED}FAIL${NC} (file '$1' was not created)"
        FAIL=$((FAIL+1))
    fi
}

# Checks that file $2 contains (grep -q) pattern $1
check_contains() {
    if grep -q "$1" "$2" 2>/dev/null; then
        echo -e "   -> ${GREEN}Verified:${NC} '$1' found in '$2'"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} '$1' NOT found in '$2'"
        FAIL=$((FAIL+1))
    fi
}

# Checks that $2 does NOT contain (grep -q) pattern $1
check_not_contains() {
    if grep -q "$1" "$2" 2>/dev/null; then
        echo -e "   -> ${RED}Failed:${NC} '$1' found in '$2' (should not be there)"
        FAIL=$((FAIL+1))
    else
        echo -e "   -> ${GREEN}Verified:${NC} '$1' absent from '$2'"
        PASS=$((PASS+1))
    fi
}

# Checks that $1 (actual exit code) equals $2 (expected exit code)
check_exit_code() {
    if [ "$1" -eq "$2" ]; then
        echo -e "   -> ${GREEN}Verified:${NC} exit code $1 (expected $2)"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} exit code $1 (expected $2)"
        FAIL=$((FAIL+1))
    fi
}

# Counts the non-blank lines in a whitespace-separated .dat block belonging
# to a single k-point (bands_gnuplot.dat has one row per band, blank line
# between k-points) and compares against an expected count.
check_column_count() {
    local file="$1" expected="$2" line
    line=$(head -n 1 "$file")
    local actual
    actual=$(echo "$line" | wc -w)
    if [ "$actual" -eq "$expected" ]; then
        echo -e "   -> ${GREEN}Verified:${NC} '$file' first row has $actual columns (expected $expected)"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} '$file' first row has $actual columns (expected $expected)"
        FAIL=$((FAIL+1))
    fi
}


# --- 1. Preparation ---
echo "--- Starting tester for STB-Bands (item 2.1) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/siesta.bands" "$TEST_DIR/"
cp "$FIXTURE_DIR/siesta_edge_nbands21.bands" "$TEST_DIR/"
cp "$FIXTURE_DIR/siesta_spin2.bands" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Real fixture: --shift fermi ---
echo -e "\n--- Testing --shift fermi (real fixture, 182 bands x 501 k-points, nspin=1) ---"
rm -f bands_gnuplot.dat bands.gplot bands_analysis.txt
stb-bands --file siesta.bands --shift fermi --no-intro > log_fermi.txt 2>&1
check_exit_code $? 0
check_success bands_gnuplot.dat
check_success bands.gplot
check_success bands_analysis.txt

echo "Verifying bands_analysis.txt content"
check_contains "Fermi energy" bands_analysis.txt
check_contains "VBM" bands_analysis.txt
check_contains "CBM" bands_analysis.txt
check_contains "Band gap" bands_analysis.txt
check_contains "Gap type" bands_analysis.txt
check_contains "Indirect" bands_analysis.txt

echo "Verifying high-symmetry labels are not quoted in the gnuplot script"
check_contains "set xtics" bands.gplot
check_not_contains "'GAMMA'" bands.gplot
check_not_contains "''" bands.gplot
check_contains '{/Symbol G}' bands.gplot

echo "Verifying bands_gnuplot.dat row shape (k, energy) pairs, one band per block"
check_column_count bands_gnuplot.dat 2


# --- 3. Real fixture: other --shift modes ---
echo -e "\n--- Testing --shift vbm/cbm/manual ---"
stb-bands --file siesta.bands --shift vbm --no-intro > log_vbm.txt 2>&1
check_exit_code $? 0
stb-bands --file siesta.bands --shift cbm --no-intro > log_cbm.txt 2>&1
check_exit_code $? 0
stb-bands --file siesta.bands --shift manual --manual-value 0.5 --no-intro > log_manual.txt 2>&1
check_exit_code $? 0
stb-bands --file siesta.bands --shift manual --no-intro > log_manual_missing.txt 2>&1
check_exit_code $? 2


# --- 4. Edge-case fixture: nbands mod 10 == 1 (old parser dropped a value here) ---
echo -e "\n--- Testing continuation-line edge case (21 bands, trailing 1-value line) ---"
rm -f bands_gnuplot.dat bands.gplot bands_analysis.txt
stb-bands --file siesta_edge_nbands21.bands --shift fermi --no-intro > log_edge.txt 2>&1
check_exit_code $? 0
check_success bands_gnuplot.dat

echo "Verifying row shape (k, energy) pairs"
check_column_count bands_gnuplot.dat 2

echo "Verifying no energy value was dropped: 21 bands x (4 k-points + 1 blank separator) = 105 lines"
lines_total=$(wc -l < bands_gnuplot.dat)
if [ "$lines_total" -eq 105 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} bands_gnuplot.dat has 105 lines (would be 100 if one value/k-point were dropped)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} bands_gnuplot.dat has $lines_total lines (expected 105)"
    FAIL=$((FAIL+1))
fi

echo "Verifying repeated Gamma labels (start and end of path) both render clean"
check_contains "Gap type" bands_analysis.txt
check_contains '{/Symbol G}' bands.gplot
check_not_contains "'GAMMA'" bands.gplot


# --- 4b. Synthetic nspin=2 fixture: half-metallic (spin-up metallic, spin-down 2 eV gap) ---
echo -e "\n--- Testing spin-polarized bands (nspin=2, half-metallic synthetic fixture) ---"
rm -f bands_gnuplot.dat bands.gplot bands_analysis.txt
stb-bands --file siesta_spin2.bands --shift fermi --no-intro > log_spin2.txt 2>&1
check_exit_code $? 0
check_success bands_gnuplot.dat
check_success bands.gplot
check_success bands_analysis.txt

echo "Verifying bands_gnuplot.dat has 3 columns (k, spin-up, spin-down)"
check_column_count bands_gnuplot.dat 3

echo "Verifying bands.gplot has one plot command per spin channel"
check_contains 'using 1:2 with lines ls 1 title "Spin Up"' bands.gplot
check_contains 'using 1:3 with lines ls 2 title "Spin Down"' bands.gplot

echo "Verifying bands_analysis.txt reports per-spin channel + half-metallic character"
check_contains "Per-spin channel" bands_analysis.txt
check_contains "Spin Up" bands_analysis.txt
check_contains "Spin Down" bands_analysis.txt
check_contains "Half-metallic: Yes" bands_analysis.txt

echo "Testing: --gap-tol tightened below the synthetic up-channel gap (no longer half-metallic)"
stb-bands --file siesta_spin2.bands --shift fermi --gap-tol 0.0001 --no-intro > log_spin2_tol.txt 2>&1
check_exit_code $? 0
check_contains "Half-metallic: No" bands_analysis.txt


# --- 5. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: missing --file (required)"
stb-bands --shift fermi --no-intro > log_missing_file_arg.txt 2>&1
check_exit_code $? 2

echo "Testing: missing --shift (required)"
stb-bands --file siesta.bands --no-intro > log_missing_shift.txt 2>&1
check_exit_code $? 2

echo "Testing: nonexistent input file"
stb-bands --file does_not_exist.bands --shift fermi --no-intro > log_missing_input.txt 2>&1
check_exit_code $? 1

echo "Testing: --version"
stb-bands --version > log_version.txt 2>&1
check_contains "stb-bands" log_version.txt


# --- 6. Interactive path (stb-suite, shortcut 2.1) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 2.1) ---"

echo "Testing: navigate 2.1 -> siesta.bands -> shift by Fermi level"
rm -f bands_gnuplot.dat bands.gplot bands_analysis.txt
printf '2.1\nsiesta.bands\n3\n' | stb-suite > log_menu.txt 2>&1
check_success bands_gnuplot.dat
check_success bands_analysis.txt


popd > /dev/null

# --- 7. Summary ---
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
