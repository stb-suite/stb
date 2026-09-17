#!/bin/bash

# --- Setup ---
# Smoke test for stb-hirshfeldPrep (Hirshfeld-I Stage 1: Prep, item
# 4.20.1). Fixtures: a synthetic 2-species (C, O) combined structure +
# calc.fdf + placeholder pseudopotentials. Stage 1 does NOT require
# anything pre-computed -- it writes both combined/ (the full system,
# forced to single-point + SaveRho) and neutral/<species>/ itself, so this
# test needs no synthetic .RHO at all (unlike ions/analysis's own tests).
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
echo "--- Starting tester for STB-HIRSHFELDPREP (item 4.20.1) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/structure.fdf" "$FIXTURE_DIR/calc.fdf" "$FIXTURE_DIR/C.psf" "$FIXTURE_DIR/O.psf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Default run (no pre-computed .RHO anywhere) ---
echo -e "\n--- Testing default run (writes combined/ + neutral/<species>/, nothing pre-computed) ---"
rm -rf hirshfeld_study
stb-hirshfeldPrep -s structure.fdf --calc calc.fdf -p . -O hirshfeld_study --no-intro \
    > log_default.txt 2>&1
check_exit_code $? 0
check_success hirshfeld_study/hirshfeld_manifest.json
check_success hirshfeld_study/combined/structure.fdf
check_success hirshfeld_study/combined/calc.fdf
check_success hirshfeld_study/combined/config_extra.fdf
check_success hirshfeld_study/combined/C.psf
check_success hirshfeld_study/combined/O.psf
check_success hirshfeld_study/neutral/C/structure.fdf
check_success hirshfeld_study/neutral/C/calc.fdf
check_success hirshfeld_study/neutral/C/config_extra.fdf
check_success hirshfeld_study/neutral/C/C.psf
check_success hirshfeld_study/neutral/O/structure.fdf
check_success hirshfeld_study/neutral/O/O.psf

echo "Testing: manifest records both species + the combined/ folder just written (absolute paths)"
check_contains '"symbol": "C"' hirshfeld_study/hirshfeld_manifest.json
check_contains '"symbol": "O"' hirshfeld_study/hirshfeld_manifest.json
check_contains '"neutral_folder": "neutral/C"' hirshfeld_study/hirshfeld_manifest.json
check_contains "hirshfeld_study/combined" hirshfeld_study/hirshfeld_manifest.json

echo "Testing: combined/ has the full 2-atom system, fixed cell + single-point + SaveRho"
echo "         + mandatory native Hirshfeld/Voronoi charges, but NO forced k-grid/Spin"
echo "         (left as whatever calc.fdf already has)"
check_contains "NumberofAtoms      2" hirshfeld_study/combined/structure.fdf
check_contains "MD.VariableCell false" hirshfeld_study/combined/config_extra.fdf
check_contains "MD.Steps              0" hirshfeld_study/combined/config_extra.fdf
check_contains "MeshCutoff          400 Ry" hirshfeld_study/combined/config_extra.fdf
check_contains "SaveRho             true" hirshfeld_study/combined/config_extra.fdf
check_contains "Charge.Hirshfeld end" hirshfeld_study/combined/config_extra.fdf
check_contains "Charge.Voronoi end" hirshfeld_study/combined/config_extra.fdf
check_contains "%include config_extra.fdf" hirshfeld_study/combined/calc.fdf

echo "Testing: neutral/C is a 1-atom, single-point, Gamma-only, spin-polarized cell"
check_contains "NumberofAtoms      1" hirshfeld_study/neutral/C/structure.fdf
check_contains "MD.Steps              0" hirshfeld_study/neutral/C/config_extra.fdf
check_contains "kgrid.MonkhorstPack   \[1  1  1\]" hirshfeld_study/neutral/C/config_extra.fdf
check_contains "Spin                polarized" hirshfeld_study/neutral/C/config_extra.fdf
check_contains "MeshCutoff          400 Ry" hirshfeld_study/neutral/C/config_extra.fdf
check_contains "SaveRho             true" hirshfeld_study/neutral/C/config_extra.fdf
check_contains "%include config_extra.fdf" hirshfeld_study/neutral/C/calc.fdf

echo "Testing: Spin polarized is ALSO forced directly into the literal calc.fdf body"
echo "         (not just via config_extra.fdf) -- the fixture's own calc.fdf has"
echo "         'Spin          non-polarized', which must be overwritten, not just shadowed"
check_contains "Spin                polarized" hirshfeld_study/neutral/C/calc.fdf
grep -q "Spin          non-polarized" hirshfeld_study/neutral/C/calc.fdf
if [ $? -ne 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} original 'Spin          non-polarized' line no longer present"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} original 'Spin          non-polarized' line still present"
    FAIL=$((FAIL+1))
fi


# --- 3. --mesh-cutoff / --vacuum overrides ---
echo -e "\n--- Testing --mesh-cutoff/--vacuum overrides ---"
rm -rf hirshfeld_study_override
stb-hirshfeldPrep -s structure.fdf --calc calc.fdf -p . -O hirshfeld_study_override \
    --mesh-cutoff 250 --vacuum 12 --no-intro > log_override.txt 2>&1
check_exit_code $? 0
check_contains "MeshCutoff          250 Ry" hirshfeld_study_override/combined/config_extra.fdf
check_contains "MeshCutoff          250 Ry" hirshfeld_study_override/neutral/O/config_extra.fdf
check_contains " 12.000000   0.000000   0.000000" hirshfeld_study_override/neutral/O/structure.fdf


# --- 4. Missing structure file ---
echo -e "\n--- Testing missing structure file ---"
stb-hirshfeldPrep -s does_not_exist.fdf --calc calc.fdf -p . -O hirshfeld_study_missing \
    --no-intro > log_missingstruct.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missingstruct.txt


# --- 5. Interactive path (stb-suite, shortcut 4.20.1) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.20.1) ---"
rm -rf hirshfeld_study_interactive
printf '4.20.1\nstructure.fdf\ncalc.fdf\n3\n.\n400\n20\nhirshfeld_study_interactive\nn\n\n0\n' \
    | stb-suite > log_interactive.txt 2>&1
check_contains "Hirshfeld-I Prep (Stage 1) complete" log_interactive.txt
check_success hirshfeld_study_interactive/hirshfeld_manifest.json


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
