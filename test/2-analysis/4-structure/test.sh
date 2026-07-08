#!/bin/bash

# --- Setup ---
# Smoke test for stb-structural (Structure Analyzer / ECN, item 2.4)
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$FIXTURE_DIR/test_files"

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
    if grep -q -- "$1" "$2" 2>/dev/null; then
        echo -e "   -> ${GREEN}Verified:${NC} '$1' found in '$2'"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} '$1' NOT found in '$2'"
        FAIL=$((FAIL+1))
    fi
}

# Checks that $2 does NOT contain (grep -q) pattern $1
check_not_contains() {
    if grep -q -- "$1" "$2" 2>/dev/null; then
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

# Checks that file $1 has exactly $2 lines
check_line_count() {
    local actual
    actual=$(wc -l < "$1")
    if [ "$actual" -eq "$2" ]; then
        echo -e "   -> ${GREEN}Verified:${NC} '$1' has $actual lines (expected $2)"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} '$1' has $actual lines (expected $2)"
        FAIL=$((FAIL+1))
    fi
}


# --- 1. Preparation ---
echo "--- Starting tester for STB-Structural (item 2.4) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/siesta.STRUCT_OUT" "$TEST_DIR/"
cp "$FIXTURE_DIR/structure.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. --format fdf, mode mean (14-atom Sn3O4 structure) ---
echo -e "\n--- Testing --format fdf --mode mean ---"
rm -f structural_information.dat warnings.log
stb-structural --file structure.fdf --format fdf --mode mean --no-intro > log_fdf_mean.txt 2>&1
check_exit_code $? 0
check_success structural_information.dat

echo "Verifying lattice parameters"
check_contains "a = 4.883 Å" structural_information.dat
check_contains "b = 5.907 Å" structural_information.dat
check_contains "c = 8.238 Å" structural_information.dat

echo "Verifying effective (weighted) CN is broken down per species, not lumped into one"
check_contains "Effective Coordination Number (weighted), per species:" structural_information.dat
check_contains " O (8 atoms):" structural_information.dat
check_contains " Sn (6 atoms):" structural_information.dat
check_contains "CrystalNN      : 2.956" structural_information.dat
check_contains "CrystalNN      : 4.004" structural_information.dat

echo "Verifying bond distance is broken down per species pair, not a single lumped average"
check_contains "Average bond distance, per species pair:" structural_information.dat
check_contains "O-Sn      : 2.1109 Å  (n=48)" structural_information.dat


# --- 3. --format struct_out, mode mean -- same physical structure (post-
#     relaxation), so lattice/CN/bond-distance summaries should match the
#     .fdf run to the precision printed (atomic positions differ in the
#     last digit or two, which isn't checked here). ---
echo -e "\n--- Testing --format struct_out --mode mean (same structure, different source format) ---"
rm -f structural_information.dat warnings.log
stb-structural --file siesta.STRUCT_OUT --format struct_out --mode mean --no-intro > log_struct_mean.txt 2>&1
check_exit_code $? 0
check_contains "a = 4.883 Å" structural_information.dat
check_contains "CrystalNN      : 2.956" structural_information.dat
check_contains "O-Sn      : 2.1109 Å  (n=48)" structural_information.dat


# --- 4. --mode list: per-atom effective CN for specific 1-based indices ---
echo -e "\n--- Testing --mode list ---"
rm -f structural_information.dat warnings.log
stb-structural --file structure.fdf --format fdf --mode list --list 1,7 --no-intro > log_list.txt 2>&1
check_exit_code $? 0
check_contains "Effective Coordination Number (weighted) for specified atoms" structural_information.dat
check_contains "Element: Sn     Cartesian Position:" structural_information.dat
check_contains "Element: O     Cartesian Position:" structural_information.dat

echo "Verifying CrystalNN no longer errors under use_weights=True (needs weighted_cn=True at construction)"
check_not_contains "CrystalNN failed" log_list.txt


# --- 5. --output-dir: writes into (and creates) a chosen directory ---
echo -e "\n--- Testing --output-dir ---"
rm -rf out_dir
stb-structural --file structure.fdf --format fdf --mode mean --output-dir out_dir --no-intro > log_outdir.txt 2>&1
check_exit_code $? 0
check_success out_dir/structural_information.dat
check_success out_dir/warnings.log


# --- 6. warnings.log doesn't grow across repeated runs (logging.basicConfig
#     defaults to append mode, which used to let it accumulate stale
#     warnings from unrelated previous runs forever) ---
echo -e "\n--- Testing warnings.log is reset each run, not appended to ---"
rm -f warnings.log
stb-structural --file structure.fdf --format fdf --mode mean --no-intro > /dev/null 2>&1
lines_first=$(wc -l < warnings.log)
stb-structural --file structure.fdf --format fdf --mode mean --no-intro > /dev/null 2>&1
lines_second=$(wc -l < warnings.log)
if [ "$lines_first" -eq "$lines_second" ]; then
    echo -e "   -> ${GREEN}Verified:${NC} warnings.log stayed at $lines_first lines across two runs (not appended)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} warnings.log grew from $lines_first to $lines_second lines across two runs"
    FAIL=$((FAIL+1))
fi


# --- 7. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: atom index out of range"
stb-structural --file structure.fdf --format fdf --mode list --list 1,99 --no-intro > log_oob.txt 2>&1
check_exit_code $? 1
check_contains "out of range" log_oob.txt

echo "Testing: nonexistent input file"
stb-structural --file does_not_exist.fdf --format fdf --mode mean --no-intro > log_missing.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing.txt

echo "Testing: --mode list without --list"
stb-structural --file structure.fdf --format fdf --mode list --no-intro > log_missing_list.txt 2>&1
check_exit_code $? 2

echo "Testing: invalid --format (cif/poscar no longer accepted)"
stb-structural --file structure.fdf --format cif --mode mean --no-intro > log_badfmt.txt 2>&1
check_exit_code $? 2

echo "Testing: missing required arguments"
stb-structural --file structure.fdf --no-intro > log_missing_args.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-structural --version > log_version.txt 2>&1
check_contains "stb-structural" log_version.txt


# --- 8. Interactive path (stb-suite, shortcut 2.4) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 2.4) ---"

echo "Testing: navigate 2.4 -> structure.fdf -> fdf -> mean -> default output dir"
rm -f structural_information.dat warnings.log
printf '2.4\nstructure.fdf\nfdf\nmean\n\n' | stb-suite > log_menu.txt 2>&1
check_success structural_information.dat


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
