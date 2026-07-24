#!/bin/bash

# --- Setup ---
# Smoke test for stb-unitcell (Unit Cell Finder, item 2.7)
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
echo "--- Starting tester for STB-Unitcell (item 2.7) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/fcc_ni.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/graphene.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/fcc_ni_noisy.fdf" "$TEST_DIR/"
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


# --- 2. Primitive mode: real reduction (FCC Ni, 4 -> 1 atom) ---
echo -e "\n--- Testing primitive mode (reduction case) ---"

echo "Testing: default mode (no --mode) behaves as primitive"
rm -f unitcell.fdf
stb-unitcell -f fcc_ni.fdf --no-intro > log_default.txt 2>&1
check_contains "Output atoms       : 1" log_default.txt
check_contains "Space group    : Fm-3m" log_default.txt
check_success unitcell.fdf
check_contains "NumberofAtoms      1" unitcell.fdf
check_contains "Unit cell reduced by stb-unitcell" unitcell.fdf
mv unitcell.fdf unitcell_primitive.fdf

echo "Testing: explicit --mode primitive"
rm -f unitcell.fdf
stb-unitcell -f fcc_ni.fdf --mode primitive --no-intro > log_primitive.txt 2>&1
check_contains "Output atoms       : 1" log_primitive.txt
check_success unitcell.fdf
rm -f unitcell.fdf

echo "Testing: stb_unitcell_report.txt is NOT written unless --save-report is given"
[ ! -f stb_unitcell_report.txt ] && { echo -e "   -> ${GREEN}Verified:${NC} stb_unitcell_report.txt absent by default"; PASS=$((PASS+1)); } \
    || { echo -e "   -> ${RED}Failed:${NC} stb_unitcell_report.txt was written without --save-report"; FAIL=$((FAIL+1)); }


# --- 2b. Refined mode: conventional-sized cell, positions snapped to symmetry ---
echo -e "\n--- Testing refined mode ---"

echo "Testing: --mode refined on fcc_ni.fdf (conventional-sized, 4 -> 4, no false 'already' note)"
rm -f unitcell.fdf
stb-unitcell -f fcc_ni.fdf --mode refined --no-intro > log_refined.txt 2>&1
check_contains "Output atoms       : 4" log_refined.txt
if grep -q "already the refined cell" log_refined.txt; then
    echo -e "   -> ${RED}Failed:${NC} unexpected 'already the refined cell' note for --mode refined"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no spurious 'already the refined cell' note"
    PASS=$((PASS+1))
fi
check_success unitcell.fdf
check_contains "NumberofAtoms      4" unitcell.fdf
rm -f unitcell.fdf

echo "Testing: --mode refined on fcc_ni_noisy.fdf (visibly snaps noisy positions to 0/0.5)"
rm -f unitcell.fdf
stb-unitcell -f fcc_ni_noisy.fdf --mode refined --no-intro > log_refined_noisy.txt 2>&1
check_contains "Space group    : Fm-3m" log_refined_noisy.txt
check_contains "Output atoms       : 4" log_refined_noisy.txt
check_success unitcell.fdf
if grep -qE "0\.9999|0\.0000[1-9]|0\.5000[1-9]|0\.4999" unitcell.fdf; then
    echo -e "   -> ${RED}Failed:${NC} refined output still contains noisy (non-exact) coordinates"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} noisy input coordinates were snapped to exact values"
    PASS=$((PASS+1))
fi
rm -f unitcell.fdf


# --- 3. Conventional mode: no-op on an already-conventional cell ---
echo -e "\n--- Testing conventional mode (no-op case) ---"

echo "Testing: --mode conventional on fcc_ni.fdf (already conventional, 4 -> 4)"
stb-unitcell -f fcc_ni.fdf --mode conventional --no-intro > log_conventional.txt 2>&1
check_contains "Output atoms       : 4" log_conventional.txt
check_contains "already the conventional cell" log_conventional.txt
check_success unitcell.fdf
rm -f unitcell.fdf


# --- 4. Already-primitive input (graphene, 2 -> 2) ---
echo -e "\n--- Testing already-primitive input ---"

echo "Testing: primitive mode on graphene.fdf (already primitive, 2 -> 2, vacuum axis preserved)"
stb-unitcell -f graphene.fdf --no-intro > log_graphene.txt 2>&1
check_contains "Output atoms       : 2" log_graphene.txt
check_contains "already the primitive cell" log_graphene.txt
check_success unitcell.fdf
check_contains "20.00000000" unitcell.fdf
rm -f unitcell.fdf


# --- 5. --save-report ---
echo -e "\n--- Testing --save-report ---"

stb-unitcell -f fcc_ni.fdf --save-report --no-intro > log_report.txt 2>&1
check_exit_code $? 0
check_success stb_unitcell_report.txt
check_contains "Space group    : Fm-3m" stb_unitcell_report.txt
check_contains "Reduction factor   : 4x" stb_unitcell_report.txt
check_contains "\[7\] SYMMETRY ANALYSIS (BEFORE / AFTER)" stb_unitcell_report.txt
check_contains "Before/After match exactly, as expected" stb_unitcell_report.txt
check_contains "Report         : stb_unitcell_report.txt" log_report.txt
rm -f unitcell.fdf stb_unitcell_report.txt


# --- 6. --ml-relax (needs the optional 'ml' extra) ---
echo -e "\n--- Testing --ml-relax (MACE pre-optimization), if the 'ml' extra is available ---"

if python3 -c "import mace" 2>/dev/null; then
    rm -f unitcell.fdf stb_unitcell_report.txt
    stb-unitcell -f fcc_ni.fdf --ml-relax --save-report --no-intro > log_ml_relax.txt 2>&1
    check_exit_code $? 0
    check_contains "\[5\] ML PRE-RELAXATION (MACE)" stb_unitcell_report.txt
    check_success unitcell.fdf
    check_contains "ML pre-relaxed with MACE-MP-0" unitcell.fdf
    rm -f unitcell.fdf stb_unitcell_report.txt

    echo "Testing: --ml-relax-cell without --ml-relax is rejected"
    stb-unitcell -f fcc_ni.fdf --ml-relax-cell --no-intro > log_ml_relax_cell_only.txt 2>&1
    check_exit_code $? 2
else
    echo -e "   -> ${YELLOW}SKIPPED${NC}: the optional 'ml' extra is not installed "
    echo "      (pip install stb_suite[ml] to also exercise --ml-relax)."
fi


# --- 7. Round-trip: output is valid, readable fdf ---
echo -e "\n--- Verifying the generated unit cell is valid input for another tool ---"
stb-unitcell -f fcc_ni.fdf --no-intro -o prim.fdf > /dev/null 2>&1
stb-kgrid --file prim.fdf --density 0.2 --no-intro > log_roundtrip.txt 2>&1
check_contains "Suggested Monkhorst-Pack grid" log_roundtrip.txt
rm -f prim.fdf


# --- 8. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: missing structure file"
stb-unitcell -f does_not_exist.fdf --no-intro > log_missing.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing.txt

echo "Testing: malformed fdf (missing ChemicalSpeciesLabel)"
stb-unitcell -f malformed.fdf --no-intro > log_malformed.txt 2>&1
check_exit_code $? 1
check_contains "ChemicalSpeciesLabel not found" log_malformed.txt

echo "Testing: -f/--file missing (required)"
stb-unitcell --no-intro > log_missing_file_arg.txt 2>&1
check_exit_code $? 2

echo "Testing: invalid --mode value"
stb-unitcell -f fcc_ni.fdf --mode wrong --no-intro > log_bad_mode.txt 2>&1
check_exit_code $? 2

echo "Testing: -o/--output custom filename"
stb-unitcell -f fcc_ni.fdf --no-intro -o my_unitcell.fdf > log_custom_output.txt 2>&1
check_success my_unitcell.fdf
rm -f my_unitcell.fdf

echo "Testing: --version"
stb-unitcell --version > log_version.txt 2>&1
check_contains "stb-unitcell" log_version.txt

echo "Testing: --help documents --mode"
stb-unitcell --help > log_help.txt 2>&1
check_contains "primitive" log_help.txt
check_contains "conventional" log_help.txt
check_contains "refined" log_help.txt


# --- 9. Interactive path (stb-suite, menu 2.7) ---
echo -e "\n--- Testing the interactive path via stb-suite (menu 2.7) ---"

echo "Testing: navigate 2.7 -> invalid file then valid -> default mode -> default symprec -> default output ->"
echo "         no ml-relax -> no save-report -> no view -> quit"
rm -f unitcell.fdf
printf '2.7\ndoes_not_exist.fdf\nfcc_ni.fdf\n\n\n\n\n\n\n0\n' | stb-suite > log_interactive.txt 2>&1
check_contains "File not found" log_interactive.txt
check_contains "Output atoms       : 1" log_interactive.txt
check_success unitcell.fdf


popd > /dev/null

# --- 10. Summary ---
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
