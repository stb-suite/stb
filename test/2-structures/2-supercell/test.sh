#!/bin/bash

# --- Setup ---
# Smoke test for stb-supercell (Supercell Builder, item 2.2)
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
    if grep -q "$1" "$2" 2>/dev/null; then
        echo -e "   -> ${GREEN}Verified:${NC} '$1' found in '$2'"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} '$1' NOT found in '$2'"
        FAIL=$((FAIL+1))
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
echo "--- Starting tester for STB-Supercell (item 2.2) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/structure.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/graphene.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null

cat > malformed.fdf << 'EOF'
NumberOfSpecies    1
NumberofAtoms      1
LatticeConstant 1.0 Ang
AtomicCoordinatesFormat  Fractional
%block LatticeVectors
 5.0 0.0 0.0
 0.0 5.0 0.0
 0.0 0.0 5.0
%endblock LatticeVectors
%block AtomicCoordinatesAndAtomicSpecies
 0.0 0.0 0.0 1
%endblock AtomicCoordinatesAndAtomicSpecies
EOF


# --- 2. Diagonal supercells ---
echo -e "\n--- Testing diagonal supercells ---"

echo "Testing: 2x2x2 diagonal (2 atoms -> 16 atoms, lattice doubled to 10.86 Ang)"
rm -f supercell.fdf
stb-supercell -f structure.fdf -d 2 2 2 --no-intro > log_2x2x2.txt 2>&1
check_contains "Output atoms:.*16" log_2x2x2.txt
check_contains "Determinant (cell multiplication factor):.*8" log_2x2x2.txt
check_success supercell.fdf
check_contains "NumberofAtoms      16" supercell.fdf
check_contains "10.86000000" supercell.fdf
mv supercell.fdf supercell_2x2x2.fdf

echo "Testing: non-uniform diagonal 1x2x3 (2 atoms -> 12 atoms, determinant 6)"
rm -f supercell.fdf
stb-supercell -f structure.fdf -d 1 2 3 --no-intro > log_1x2x3.txt 2>&1
check_contains "Output atoms:.*12" log_1x2x3.txt
check_contains "Determinant (cell multiplication factor):.*6" log_1x2x3.txt
check_success supercell.fdf
rm -f supercell.fdf


# --- 3. Full 3x3 matrix (non-orthogonal lattice) ---
echo -e "\n--- Testing the full row-major 3x3 matrix path ---"

echo "Testing: sqrt(3)xsqrt(3) reconstruction on graphene (2 atoms -> 6 atoms, determinant 3)"
stb-supercell -f graphene.fdf -d 1 1 0 -1 2 0 0 0 1 --no-intro > log_sqrt3.txt 2>&1
check_contains "Output atoms:.*6" log_sqrt3.txt
check_contains "Determinant (cell multiplication factor):.*3" log_sqrt3.txt
check_success supercell.fdf
mv supercell.fdf supercell_sqrt3.fdf


# --- 4. Round-trip: output is valid, readable fdf ---
echo -e "\n--- Verifying the generated supercell is valid input for another tool ---"
stb-kgrid --file supercell_2x2x2.fdf --type fdf --density 0.2 --no-intro > log_roundtrip.txt 2>&1
check_contains "Suggested Monkhorst-Pack grid: 3 3 3" log_roundtrip.txt


# --- 5. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: singular matrix (determinant zero)"
stb-supercell -f structure.fdf -d 1 0 0 0 1 0 0 0 0 --no-intro > log_singular.txt 2>&1
check_exit_code $? 1
check_contains "singular" log_singular.txt

echo "Testing: non-integer entries"
stb-supercell -f structure.fdf -d 1.5 1 1 --no-intro > log_noninteger.txt 2>&1
check_exit_code $? 1
check_contains "must have integer entries" log_noninteger.txt

echo "Testing: wrong count of numbers (4, neither 3 nor 9)"
stb-supercell -f structure.fdf -d 1 2 3 4 --no-intro > log_wrong_count.txt 2>&1
check_exit_code $? 1
check_contains "Expected 3 (diagonal) or 9 (full matrix)" log_wrong_count.txt

echo "Testing: negative determinant (mirrored cell) -> warning, still succeeds"
rm -f supercell.fdf
stb-supercell -f structure.fdf -d -1 0 0 0 1 0 0 0 1 --no-intro > log_negative_det.txt 2>&1
check_exit_code $? 0
check_contains "negative determinant" log_negative_det.txt
check_success supercell.fdf
rm -f supercell.fdf

echo "Testing: missing structure file"
stb-supercell -f does_not_exist.fdf -d 2 2 2 --no-intro > log_missing.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing.txt

echo "Testing: malformed fdf (missing ChemicalSpeciesLabel)"
stb-supercell -f malformed.fdf -d 2 2 2 --no-intro > log_malformed.txt 2>&1
check_exit_code $? 1
check_contains "ChemicalSpeciesLabel not found" log_malformed.txt

echo "Testing: -f/--file missing (required)"
stb-supercell -d 2 2 2 --no-intro > log_missing_file_arg.txt 2>&1
check_exit_code $? 2

echo "Testing: -d/--dim missing (required)"
stb-supercell -f structure.fdf --no-intro > log_missing_dim_arg.txt 2>&1
check_exit_code $? 2

echo "Testing: -o/--output custom filename"
stb-supercell -f structure.fdf -d 2 2 2 -o my_supercell.fdf --no-intro > log_custom_output.txt 2>&1
check_success my_supercell.fdf
rm -f my_supercell.fdf

echo "Testing: --version"
stb-supercell --version > log_version.txt 2>&1
check_contains "stb-supercell" log_version.txt

echo "Testing: --help documents the dim convention"
stb-supercell --help > log_help.txt 2>&1
check_contains "diagonal supercell" log_help.txt
check_contains "row-major" log_help.txt


# --- 6. Interactive path (stb-suite, shortcut 1.5) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 1.5) ---"

echo "Testing: navigate 1.5 -> invalid file then valid -> dims '2 2 2' -> default output -> quit"
rm -f supercell.fdf
printf '2.2\ndoes_not_exist.fdf\nstructure.fdf\n2 2 2\n\n0\n' | stb-suite > log_interactive.txt 2>&1
check_contains "File not found" log_interactive.txt
check_contains "Output atoms:.*16" log_interactive.txt
check_success supercell.fdf


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
