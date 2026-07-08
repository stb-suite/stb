#!/bin/bash

# --- Setup ---
# Smoke test for stb-symmetry (Crystal Symmetry Analyzer, item 2.5)
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


# --- 1. Preparation ---
echo "--- Starting tester for STB-Symmetry (item 2.5) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/structure.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/siesta.STRUCT_OUT" "$TEST_DIR/"
cp "$FIXTURE_DIR/nacl.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. --format fdf on a P1 (no symmetry) structure ---
echo -e "\n--- Testing --format fdf on the P1 Sn3O4 structure ---"
rm -f symmetry.dat
stb-symmetry --file structure.fdf --format fdf --no-intro > log_p1.txt 2>&1
check_exit_code $? 0
check_success symmetry.dat

echo "Verifying the report has a title block and generation timestamp"
check_contains "CRYSTAL SYMMETRY REPORT - STB Suite" symmetry.dat
check_contains "Generated        :" symmetry.dat
check_contains "Source file      : structure.fdf  (format: fdf)" symmetry.dat

echo "Verifying space group / point group / crystal system for a P1 (no symmetry) structure"
check_contains "Space group      : P1 (No. 1)" symmetry.dat
check_contains "Hall symbol      : P 1" symmetry.dat
check_contains "Point group      : 1" symmetry.dat
check_contains "Crystal system   : triclinic" symmetry.dat
check_contains "Pearson symbol   : aP14" symmetry.dat
check_contains "Symmetry operations: 1" symmetry.dat

echo "Verifying lattice parameters and volume"
check_contains "a = 4.883 Å   b = 5.907 Å   c = 8.238 Å" symmetry.dat
check_contains "Volume = 236.821 Å³" symmetry.dat

echo "Verifying every atom is its own Wyckoff orbit (14 distinct sites for 14 atoms, no symmetry)"
check_contains "SYMMETRICALLY DISTINCT SITES: 14" symmetry.dat
check_contains "   1  Sn   1a  " symmetry.dat


# --- 3. --format struct_out -- same physical structure, same symmetry result ---
echo -e "\n--- Testing --format struct_out (same P1 structure, different source format) ---"
rm -f symmetry.dat
stb-symmetry --file siesta.STRUCT_OUT --format struct_out --no-intro > log_struct_out.txt 2>&1
check_exit_code $? 0
check_contains "Space group      : P1 (No. 1)" symmetry.dat
check_contains "Volume = 236.821 Å³" symmetry.dat


# --- 4. A real, high-symmetry structure: NaCl rock salt (Fm-3m, No. 225) --
#     exercises Wyckoff multiplicity (4a/4b) and a large symmetry-operation
#     count, neither of which the P1 fixture can. ---
echo -e "\n--- Testing on NaCl rock salt (space group Fm-3m, No. 225) ---"
rm -f symmetry.dat
stb-symmetry --file nacl.fdf --format fdf --no-intro > log_nacl.txt 2>&1
check_exit_code $? 0

echo "Verifying space group / Pearson symbol / lattice type for a textbook Fm-3m structure"
check_contains "Space group      : Fm-3m (No. 225)" symmetry.dat
check_contains "Hall symbol      : -F 4 2 3" symmetry.dat
check_contains "Point group      : m-3m" symmetry.dat
check_contains "Crystal system   : cubic" symmetry.dat
check_contains "Lattice type     : cubic" symmetry.dat
check_contains "Pearson symbol   : cF8" symmetry.dat
check_contains "Symmetry operations: 192" symmetry.dat

echo "Verifying the two Wyckoff orbits (4a for Na, 4b for Cl), not 8 separate 1-atom sites"
check_contains "SYMMETRICALLY DISTINCT SITES: 2" symmetry.dat
check_contains "4a        Na        4         1" symmetry.dat
check_contains "4b        Cl        4         5" symmetry.dat
check_contains "   1  Na   4a  " symmetry.dat
check_contains "   5  Cl   4b  " symmetry.dat

echo "Verifying all 192 symmetry operations are listed, in compact x,y,z notation"
check_contains "SYMMETRY OPERATIONS (192), in x,y,z notation" symmetry.dat
check_contains "   1: x, y, z" symmetry.dat
check_contains "  49: x, y+1/2, z+1/2" symmetry.dat
check_contains " 192: -x+1/2, z+1/2, -y" symmetry.dat

n_op_lines=$(grep -c "^ *[0-9]*: " symmetry.dat)
if [ "$n_op_lines" -eq 192 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} exactly 192 operation lines found"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} found $n_op_lines operation lines, expected 192"
    FAIL=$((FAIL+1))
fi


# --- 5. --output-dir: writes into (and creates) a chosen directory ---
echo -e "\n--- Testing --output-dir ---"
rm -rf out_dir
stb-symmetry --file structure.fdf --format fdf --output-dir out_dir --no-intro > log_outdir.txt 2>&1
check_exit_code $? 0
check_success out_dir/symmetry.dat


# --- 6. --symprec / --angle-tolerance are configurable and validated ---
echo -e "\n--- Testing --symprec and --angle-tolerance ---"
stb-symmetry --file structure.fdf --format fdf --symprec 1e-4 --angle-tolerance 2.0 --no-intro > log_symprec.txt 2>&1
check_exit_code $? 0
check_contains "symprec=0.0001, angle_tolerance=2°" symmetry.dat

echo "Testing: --symprec 0 is rejected"
stb-symmetry --file structure.fdf --format fdf --symprec 0 --no-intro > log_symprec_zero.txt 2>&1
check_exit_code $? 2

echo "Testing: --angle-tolerance 0 is rejected"
stb-symmetry --file structure.fdf --format fdf --angle-tolerance 0 --no-intro > log_angletol_zero.txt 2>&1
check_exit_code $? 2


# --- 7. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: nonexistent input file"
stb-symmetry --file does_not_exist.fdf --format fdf --no-intro > log_missing.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing.txt
check_not_contains "Traceback" log_missing.txt

echo "Testing: invalid --format (cif/poscar no longer accepted)"
stb-symmetry --file structure.fdf --format cif --no-intro > log_badfmt.txt 2>&1
check_exit_code $? 2

echo "Testing: missing required arguments"
stb-symmetry --file structure.fdf --no-intro > log_missing_args.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-symmetry --version > log_version.txt 2>&1
check_contains "stb-symmetry" log_version.txt


# --- 8. Interactive path (stb-suite, shortcut 2.5) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 2.5) ---"

echo "Testing: navigate 2.5 -> nacl.fdf -> format=1(fdf) -> default output dir"
rm -f symmetry.dat
printf '2.5\nnacl.fdf\n1\n\n\n' | stb-suite > log_menu.txt 2>&1
check_success symmetry.dat
check_contains "Select input file format:" log_menu.txt
check_contains "Fm-3m" symmetry.dat


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
