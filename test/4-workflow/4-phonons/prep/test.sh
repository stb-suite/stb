#!/bin/bash

# --- Setup ---
# Smoke test for stb-phononsCreate (Phonons Prep, item 4.4.1)
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
echo "--- Starting tester for STB-PhononsCreate prep (item 4.4.1) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/structure.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/calc.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/Sn.psf" "$TEST_DIR/"
cp "$FIXTURE_DIR/O.psf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Basic displacement generation (-dim 1 1 1 -- same 14-atom unit cell,
#     kept small so the ML-preview test below stays fast; dim isn't what's
#     under test here) ---
echo -e "\n--- Testing basic displacement generation ---"
stb-phononsCreate -s structure.fdf -c calc.fdf -p . -dim 1 1 1 --save-report --no-intro > log_create.txt 2>&1
check_contains "Displacement folders : 42" log_create.txt
check_success phonon_runs/phonopy_disp.yaml
check_success phonon_runs/disp-001/structure.fdf
check_success phonon_runs/disp-001/calc.fdf
check_success phonon_runs/disp-001/Sn.psf
check_success phonon_runs/disp-001/O.psf
check_success phonon_runs/disp-042/structure.fdf
check_success phonon_runs/phonon_prep_properties.txt
check_contains "RUN METADATA" phonon_runs/phonon_prep_properties.txt
check_contains "Displacements needed" phonon_runs/phonon_prep_properties.txt
check_contains "42 (of 84 without symmetry reduction)" phonon_runs/phonon_prep_properties.txt

echo "Verifying the symmetry section: default symprec, highlighted space group, print_table format"
check_contains "Symmetry precision: 0.01 Ang" log_create.txt
check_contains ">>> Detected space group: P-1 (2)" log_create.txt
check_contains "Symmetry precision (symprec) | 0.01 Ang" log_create.txt
check_contains "Point group" log_create.txt
check_contains "Reduction" log_create.txt

echo "Verifying [6] LIBRARY WARNINGS section always appears (clean run -> no warnings)"
check_contains "\[6\] LIBRARY WARNINGS" log_create.txt
check_contains "No library warnings." log_create.txt

echo "Verifying the supercell k-grid is auto-suggested (default density 0.2) and written into config_extra.fdf"
check_contains "Supercell k-grid  : 7 6 4 (auto-suggested, density=0.2 1/Ang)" log_create.txt
check_contains "kgrid.MonkhorstPack   \[7  6  4\]" phonon_runs/disp-001/config_extra.fdf


echo "Verifying single-point SCF is forced (calc.fdf has MD.VariableCell .true. -- a real relaxation setting)"
check_contains "SINGLE-POINT SCF ENFORCEMENT" log_create.txt
check_contains "MD.VariableCell=.true." log_create.txt
check_success phonon_runs/disp-001/config_extra.fdf
kgrid_line=$(grep -n "kgrid.MonkhorstPack" phonon_runs/disp-001/config_extra.fdf | head -1 | cut -d: -f1)
md_line=$(grep -n "MD.TypeOfRun" phonon_runs/disp-001/config_extra.fdf | head -1 | cut -d: -f1)
if [ -n "$kgrid_line" ] && [ -n "$md_line" ] && [ "$kgrid_line" -lt "$md_line" ]; then
    echo -e "   -> ${GREEN}Verified:${NC} kgrid.MonkhorstPack (line $kgrid_line) is defined before MD.TypeOfRun (line $md_line) in config_extra.fdf"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} kgrid.MonkhorstPack does not come before MD.TypeOfRun in config_extra.fdf"
    FAIL=$((FAIL+1))
fi
check_contains "MD.VariableCell    false" phonon_runs/disp-001/config_extra.fdf
first_line=$(head -1 phonon_runs/disp-001/calc.fdf)
if [ "$first_line" = "%include config_extra.fdf" ]; then
    echo -e "   -> ${GREEN}Verified:${NC} '%include config_extra.fdf' is the FIRST line of the generated calc.fdf"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} first line is '$first_line', expected '%include config_extra.fdf'"
    FAIL=$((FAIL+1))
fi
echo "Verifying the original MD.VariableCell directive is still present further down (harmlessly overridden, not stripped)"
check_contains "MD.VariableCell  .true." phonon_runs/disp-001/calc.fdf


# --- 2b. --symprec is configurable ---
echo -e "\n--- Testing --symprec ---"
rm -rf phonon_runs
stb-phononsCreate -s structure.fdf -c calc.fdf -p . -dim 1 1 1 --symprec 1e-5 --no-intro > log_symprec.txt 2>&1
check_contains "Symmetry precision: 1e-05 Ang" log_symprec.txt
check_contains ">>> Detected space group: P1 (1)" log_symprec.txt
check_contains "Displacement folders : 84" log_symprec.txt
rm -rf phonon_runs


# --- 2c. --kgrid (explicit) / --kgrid-density (auto-suggestion knob) ---
echo -e "\n--- Testing --kgrid and --kgrid-density ---"
stb-phononsCreate -s structure.fdf -c calc.fdf -p . -dim 1 1 1 --kgrid 3 3 3 --no-intro > log_kgrid.txt 2>&1
check_contains "Supercell k-grid  : 3 3 3 (explicit --kgrid)" log_kgrid.txt
check_contains "kgrid.MonkhorstPack   \[3  3  3\]" phonon_runs/disp-001/config_extra.fdf
rm -rf phonon_runs

stb-phononsCreate -s structure.fdf -c calc.fdf -p . -dim 1 1 1 --kgrid-density 0.4 --no-intro > log_kgrid_density.txt 2>&1
check_contains "auto-suggested, density=0.4 1/Ang" log_kgrid_density.txt
rm -rf phonon_runs


# --- 3. --ml-prerelax pre-flight check (MACE-MP-0 already cached locally --
#     fast, no network needed) ---
echo -e "\n--- Testing --ml-prerelax ---"
rm -rf phonon_runs
stb-phononsCreate -s structure.fdf -c calc.fdf -p . -dim 1 1 1 --ml-prerelax --no-intro > log_ml_prerelax.txt 2>&1
check_contains "ML PRE-FLIGHT CHECK" log_ml_prerelax.txt
check_contains "Max residual force on input structure" log_ml_prerelax.txt
check_contains "Verdict:" log_ml_prerelax.txt
check_success phonon_runs/phonopy_disp.yaml

echo "Verifying MACE's own library noise (torch/cuequivariance messages) was captured into [6], not left interleaved in [0b]"
sed -n '/\[0b\] ML PRE-FLIGHT CHECK/,/\[1\] SINGLE-POINT/p' log_ml_prerelax.txt > ml_section_only.txt
check_not_contains "cuequivariance" ml_section_only.txt
check_contains "\[6\] LIBRARY WARNINGS" log_ml_prerelax.txt
check_contains "\[MACE import\]" log_ml_prerelax.txt
check_contains "cuequivariance" log_ml_prerelax.txt


# --- 4. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: missing structure file"
rm -rf phonon_runs
stb-phononsCreate -s does_not_exist.fdf -c calc.fdf -p . -dim 1 1 1 --no-intro > log_missing_struct.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_struct.txt

echo "Testing: missing calc file"
stb-phononsCreate -s structure.fdf -c does_not_exist.fdf -p . -dim 1 1 1 --no-intro > log_missing_calc.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_calc.txt

echo "Testing: missing pseudopotentials"
mkdir -p empty_pseudo_dir
stb-phononsCreate -s structure.fdf -c calc.fdf -p empty_pseudo_dir -dim 1 1 1 --no-intro > log_missing_pseudo.txt 2>&1
check_exit_code $? 1
check_contains "Missing pseudopotentials" log_missing_pseudo.txt

echo "Testing: re-running on top of an existing phonon_runs/ refuses to overwrite"
stb-phononsCreate -s structure.fdf -c calc.fdf -p . -dim 1 1 1 --no-intro > log_regen_1.txt 2>&1
stb-phononsCreate -s structure.fdf -c calc.fdf -p . -dim 1 1 1 --no-intro > log_regen_2.txt 2>&1
check_exit_code $? 1
check_contains "already contains" log_regen_2.txt
rm -rf phonon_runs

echo "Testing: --version"
stb-phononsCreate --version > log_version.txt 2>&1
check_contains "stb-phononsCreate" log_version.txt

echo "Testing: --help documents -dim/-p/--symprec/--kgrid"
stb-phononsCreate --help > log_help.txt 2>&1
check_contains "dim" log_help.txt
check_contains "pseudo" log_help.txt
check_contains "symprec" log_help.txt
check_contains "kgrid" log_help.txt


# --- 5. Interactive path (stb-suite, shortcut 4.4.1) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.4.1) ---"

echo "Testing: navigate 4.4.1 -> defaults except dim 1 1 1, accept previewed k-grid, skip ML, save report (y) -> quit"
rm -rf phonon_runs
# Prompts in order: structure file (blank -> default), calc file (blank ->
# default), dim ('1 1 1'), k-grid (blank -> accepts the previewed
# '[X Y Z]' default shown in brackets), distance (blank -> default),
# pseudo source (blank -> skip, uses '.'), ML pre-flight check (n),
# advanced settings (n -> skip), save report (y), then the "Press Enter to
# continue" pause, then quit.
printf '4.4.1\n\n\n1 1 1\n\n\n\nn\nn\ny\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_contains "Supercell k-grid \[7 6 4\]:" log_menu.txt
check_contains "Displacement folders : 42" log_menu.txt
check_contains "Supercell k-grid  : 7 6 4 (explicit --kgrid)" log_menu.txt
check_success phonon_runs/phonopy_disp.yaml
check_success phonon_runs/phonon_prep_properties.txt

echo "Testing: navigate 4.4.1 -> explicit k-grid '3 3 3' -> advanced settings (y) -> custom symprec 1e-5 -> save report (n) -> quit"
rm -rf phonon_runs
# Prompts in order: structure file (blank), calc file (blank), dim ('1 1
# 1'), k-grid ('3 3 3'), distance (blank), pseudo source (blank), ML
# pre-flight check (n), advanced settings (y), vacuum-gap (blank ->
# default), symprec ('1e-5'), save report (n), then "Press Enter to
# continue", then quit.
printf '4.4.1\n\n\n1 1 1\n3 3 3\n\n\nn\ny\n\n1e-5\nn\n\n0\n' | stb-suite > log_menu_symprec.txt 2>&1
check_contains "Symmetry precision: 1e-05 Ang" log_menu_symprec.txt
check_contains "Displacement folders : 84" log_menu_symprec.txt
check_contains "Supercell k-grid  : 3 3 3 (explicit --kgrid)" log_menu_symprec.txt

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
