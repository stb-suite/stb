#!/bin/bash

# --- Setup ---
# Smoke test for stb-cohesiveAnalysis (Cohesive Energy Analysis, item 4.3.2)
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
echo "--- Starting tester for STB-CohesiveAnalysis (item 4.3.2) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp -r "$FIXTURE_DIR/structure" "$TEST_DIR/"
cp -r "$FIXTURE_DIR/atoms" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Real fixture: single-atom-in-cell synthetic FreeEng values ---
echo -e "\n--- Testing default run (synthetic 1-atom fixture) ---"
rm -f cohesive_results.dat
stb-cohesiveAnalysis -o calc.out --no-intro > log_default.txt 2>&1
check_exit_code $? 0
check_success cohesive_results.dat

echo "Verifying the cohesive energy arithmetic (E_bulk - E_isolated) / N_atoms"
check_contains "Total atoms in cell: 1" log_default.txt
check_contains "Total Cohesive Energy:      -450.1111 eV" log_default.txt
check_contains "Cohesive Energy per Atom:   -450.1111 eV/atom" log_default.txt
check_contains "Total Cohesive Energy:          -450.1111 eV" cohesive_results.dat
check_contains "Cohesive Energy per Atom:       -450.1111 eV/atom" cohesive_results.dat


# --- 3. Zero-atom structure.fdf: used to crash with an uncaught
#     ZeroDivisionError (raw traceback) instead of a clean [ERROR] message --
#     regression test for the shared-parser migration. ---
echo -e "\n--- Testing a structure.fdf with a declared species but zero atoms ---"
mkdir -p zero_atoms/structure zero_atoms/atoms/C
cat > zero_atoms/structure/structure.fdf << 'EOF'
%block ChemicalSpeciesLabel
 1   6   C
%endblock ChemicalSpeciesLabel

LatticeConstant 1.0 Ang

AtomicCoordinatesFormat Fractional

%block LatticeVectors
 20.0 0.0 0.0
 0.0 20.0 0.0
 0.0 0.0 20.0
%endblock LatticeVectors

%block AtomicCoordinatesAndAtomicSpecies
%endblock AtomicCoordinatesAndAtomicSpecies
EOF
cp structure/calc.out zero_atoms/structure/
cp atoms/C/calc.out zero_atoms/atoms/C/
stb-cohesiveAnalysis -o calc.out -d zero_atoms --no-intro > log_zero_atoms.txt 2>&1
check_exit_code $? 1
check_contains "not found (or empty)" log_zero_atoms.txt
check_not_contains "Traceback" log_zero_atoms.txt


# --- 4. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: missing 'structure'/'atoms' directories"
stb-cohesiveAnalysis -o calc.out -d /does/not/exist --no-intro > log_missing_dirs.txt 2>&1
check_exit_code $? 1
check_contains "Required directories 'structure' and/or 'atoms' not found" log_missing_dirs.txt

echo "Testing: missing calc.out for the isolated atom"
mkdir -p missing_atom_out/structure missing_atom_out/atoms/C
cp structure/structure.fdf structure/calc.out missing_atom_out/structure/
stb-cohesiveAnalysis -o calc.out -d missing_atom_out --no-intro > log_missing_atom_out.txt 2>&1
check_exit_code $? 1
check_contains "Could not find 'calc.out' or results for isolated atom: C" log_missing_atom_out.txt
check_contains "Cannot calculate cohesive energy because some calculations are missing" log_missing_atom_out.txt

echo "Testing: missing required -o"
stb-cohesiveAnalysis --no-intro > log_missing_arg.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-cohesiveAnalysis --version > log_version.txt 2>&1
check_contains "stb-cohesiveAnalysis" log_version.txt

echo "Testing: --help documents -o/-d"
stb-cohesiveAnalysis --help > log_help.txt 2>&1
check_contains "\-\-out" log_help.txt
check_contains "\-\-dir" log_help.txt


# --- 5. Interactive path (stb-suite, shortcut 3.3.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 3.3.2) ---"

echo "Testing: navigate 3.3.2 -> calc.out -> default dir"
rm -f cohesive_results.dat
printf '4.3.2\ncalc.out\n\n' | stb-suite > log_menu.txt 2>&1
check_contains "Total Cohesive Energy:      -450.1111 eV" log_menu.txt
check_success cohesive_results.dat


popd > /dev/null

# --- 6. Summary ---
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
