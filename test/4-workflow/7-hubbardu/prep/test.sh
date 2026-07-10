#!/bin/bash

# --- Setup ---
# Smoke test for stb-hubbardu (Hubbard U Linear Response, Stage 1: Reference Prep, item 4.7.1)
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
echo "--- Starting tester for STB-Hubbardu stage 1: reference prep (item 4.7.1) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/structure.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/calc.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Default reference generation (species with a known default shell) ---
echo -e "\n--- Testing reference generation for Mn (auto-detected 3d shell) ---"
stb-hubbardu -s structure.fdf -c calc.fdf --species Mn --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "Mn (3d: n=3, l=2)" log_basic.txt
check_contains "Generated 'reference/'" log_basic.txt
check_success hubbardu_runs/run_manifest.json
check_success hubbardu_runs/reference/calc.fdf
check_success hubbardu_runs/reference/structure.fdf
check_success hubbardu_runs/_template_structure.fdf
check_success hubbardu_runs/_template_calc.fdf

echo "Testing: only the reference folder is created (no scf/frozen folders yet)"
if [ -d hubbardu_runs/scf_alpha_0.1000 ] || [ -d hubbardu_runs/frozen_alpha_0.1000 ]; then
    echo -e "   -> ${RED}Failed:${NC} stage 1 should not create scf_alpha_*/frozen_alpha_* folders"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no scf_alpha_*/frozen_alpha_* folders created by stage 1"
    PASS=$((PASS+1))
fi

echo "Testing: the generated LDAU.proj block matches the verified SIESTA syntax (U=0, J=0 default)"
check_contains "LDAU.PotentialShift T" hubbardu_runs/reference/calc.fdf
check_contains "Mn   1" hubbardu_runs/reference/calc.fdf
check_contains "n=3    2" hubbardu_runs/reference/calc.fdf
check_contains "0.000    0.000" hubbardu_runs/reference/calc.fdf

echo "Testing: run_manifest.json records species/shell/j and the reference run"
check_contains '"species": "Mn"' hubbardu_runs/run_manifest.json
check_contains '"j": 0.0' hubbardu_runs/run_manifest.json
check_contains '"kind": "reference"' hubbardu_runs/run_manifest.json


# --- 3. Explicit --shell and --j override ---
echo -e "\n--- Testing --shell and --j overrides ---"
rm -rf hubbardu_runs
stb-hubbardu -s structure.fdf -c calc.fdf --species Mn --shell 4f --j 0.5 --no-intro > log_override.txt 2>&1
check_exit_code $? 0
check_contains "Mn (4f: n=4, l=3)" log_override.txt
check_contains "0.500" log_override.txt
check_contains "n=4    3" hubbardu_runs/reference/calc.fdf
check_contains "0.000    0.500" hubbardu_runs/reference/calc.fdf
check_contains '"j": 0.5' hubbardu_runs/run_manifest.json


# --- 4. Error and robustness cases ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: species not present in the structure"
stb-hubbardu -s structure.fdf -c calc.fdf --species Fe --no-intro > log_missing_species.txt 2>&1
check_exit_code $? 1
check_contains "not found in" log_missing_species.txt

echo "Testing: species with no default shell and no --shell given"
cat > structure_si.fdf << 'EOF'
NumberOfSpecies    2
NumberofAtoms      2
%block ChemicalSpeciesLabel
 1   25   Mn
 2   14   Si
%endblock ChemicalSpeciesLabel
LatticeConstant 1.0 Ang
AtomicCoordinatesFormat  Fractional
%block LatticeVectors
 2.8900000000   0.0000000000   0.0000000000
 0.0000000000   2.8900000000   0.0000000000
 0.0000000000   0.0000000000   2.8900000000
%endblock LatticeVectors
%block AtomicCoordinatesAndAtomicSpecies
 0.000000000   0.000000000   0.000000000   1
 0.500000000   0.500000000   0.500000000   2
%endblock AtomicCoordinatesAndAtomicSpecies
EOF
stb-hubbardu -s structure_si.fdf -c calc.fdf --species Si --no-intro > log_no_shell.txt 2>&1
check_exit_code $? 1
check_contains "No default correlated shell known" log_no_shell.txt

echo "Testing: --shell explicitly supplied for a species with no default works"
stb-hubbardu -s structure_si.fdf -c calc.fdf --species Si --shell 3d -o si_runs --no-intro > log_shell_override.txt 2>&1
check_exit_code $? 0

echo "Testing: a template that already has an LDAU.proj block/PotentialShift is rejected"
cp hubbardu_runs/reference/calc.fdf already_has_ldau.fdf
stb-hubbardu -s structure.fdf -c already_has_ldau.fdf --species Mn --no-intro > log_guard.txt 2>&1
check_exit_code $? 1
check_contains "already contains an LDAU.proj block" log_guard.txt

echo "Testing: missing structure file"
stb-hubbardu -s does_not_exist.fdf -c calc.fdf --species Mn --no-intro > log_missing.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing.txt

echo "Testing: missing required args"
stb-hubbardu -s structure.fdf -c calc.fdf --no-intro > log_missing_args.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-hubbardu --version > log_version.txt 2>&1
check_contains "stb-hubbardu" log_version.txt

echo "Testing: --help documents species/shell/j and SCF convergence tips"
stb-hubbardu --help > log_help.txt 2>&1
check_contains "species" log_help.txt
check_contains "shell" log_help.txt
check_contains "\-\-j J" log_help.txt
check_contains "DM.MixingWeight" log_help.txt
check_contains "SCF.H.Tolerance" log_help.txt


# --- 5. Interactive path (stb-suite, shortcut 4.7.1) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.7.1) ---"

echo "Testing: navigate 4.7.1 -> invalid file then valid -> calc.fdf -> Mn -> auto shell -> default J -> default output -> quit"
rm -rf hubbardu_runs
printf '4.7.1\ndoes_not_exist.fdf\nstructure.fdf\ncalc.fdf\nMn\n\n\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_contains "File not found" log_menu.txt
check_contains "Generated 'reference/'" log_menu.txt
check_success hubbardu_runs/run_manifest.json


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
