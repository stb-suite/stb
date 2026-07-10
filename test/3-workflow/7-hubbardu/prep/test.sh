#!/bin/bash

# --- Setup ---
# Smoke test for stb-hubbardu (Hubbard U Linear Response, Prep, item 4.7.1)
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
echo "--- Starting tester for STB-Hubbardu prep (item 4.7.1) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/structure.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/calc.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Default sweep (species with a known default shell) ---
echo -e "\n--- Testing a default sweep for Mn (auto-detected 3d shell) ---"
stb-hubbardu -s structure.fdf -c calc.fdf --species Mn --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "Mn (3d: n=3, l=2)" log_basic.txt
check_contains "Generated 13 run(s)" log_basic.txt
check_success hubbardu_runs/run_manifest.json
check_success hubbardu_runs/reference/calc.fdf
check_success hubbardu_runs/scf_alpha_-0.1500/calc.fdf
check_success hubbardu_runs/frozen_alpha_0.1500/calc.fdf

echo "Testing: the generated LDAU.proj block matches the verified SIESTA syntax"
check_contains "LDAU.PotentialShift T" hubbardu_runs/scf_alpha_0.1000/calc.fdf
check_contains "Mn   1" hubbardu_runs/scf_alpha_0.1000/calc.fdf
check_contains "n=3    2" hubbardu_runs/scf_alpha_0.1000/calc.fdf
check_contains "0.100    0.000" hubbardu_runs/scf_alpha_0.1000/calc.fdf

echo "Testing: the frozen-density folders additionally get MaxSCFIterations/DM.UseSaveDM"
check_contains "MaxSCFIterations 1" hubbardu_runs/frozen_alpha_0.1000/calc.fdf
check_contains "DM.UseSaveDM T" hubbardu_runs/frozen_alpha_0.1000/calc.fdf
if grep -q "MaxSCFIterations" hubbardu_runs/scf_alpha_0.1000/calc.fdf; then
    echo -e "   -> ${RED}Failed:${NC} scf_alpha_0.1000 should NOT have MaxSCFIterations"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} scf_alpha_0.1000 has no MaxSCFIterations override"
    PASS=$((PASS+1))
fi

echo "Testing: reference folder also carries LDAU.PotentialShift (needed for the occupation printout)"
check_contains "LDAU.PotentialShift T" hubbardu_runs/reference/calc.fdf
check_contains "0.000    0.000" hubbardu_runs/reference/calc.fdf

echo "Testing: run_manifest.json records species/shell and every run's kind+alpha"
check_contains '"species": "Mn"' hubbardu_runs/run_manifest.json
check_contains '"kind": "reference"' hubbardu_runs/run_manifest.json
check_contains '"kind": "frozen"' hubbardu_runs/run_manifest.json


# --- 3. Explicit --shell and --alphas override ---
echo -e "\n--- Testing --shell and --alphas overrides ---"
rm -rf hubbardu_runs
stb-hubbardu -s structure.fdf -c calc.fdf --species Mn --shell 4f --alphas -0.2 0.2 --no-intro > log_override.txt 2>&1
check_exit_code $? 0
check_contains "Mn (4f: n=4, l=3)" log_override.txt
check_contains "Generated 5 run(s)" log_override.txt
check_contains "n=4    3" hubbardu_runs/scf_alpha_0.2000/calc.fdf


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
stb-hubbardu -s structure_si.fdf -c calc.fdf --species Si --shell 3d --alphas 0.1 --no-intro > log_shell_override.txt 2>&1
check_exit_code $? 0

echo "Testing: alpha=0.0 is rejected (already covered by 'reference')"
stb-hubbardu -s structure.fdf -c calc.fdf --species Mn --alphas 0.0 0.1 --no-intro > log_zero_alpha.txt 2>&1
check_exit_code $? 1
check_contains "is not a valid perturbation strength" log_zero_alpha.txt

echo "Testing: duplicate alphas are rejected"
stb-hubbardu -s structure.fdf -c calc.fdf --species Mn --alphas 0.1 0.1 --no-intro > log_dup_alpha.txt 2>&1
check_exit_code $? 1
check_contains "duplicate values" log_dup_alpha.txt

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

echo "Testing: --help documents species/shell/alphas"
stb-hubbardu --help > log_help.txt 2>&1
check_contains "species" log_help.txt
check_contains "shell" log_help.txt
check_contains "alphas" log_help.txt


# --- 5. Interactive path (stb-suite, shortcut 4.7.1) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.7.1) ---"

echo "Testing: navigate 4.7.1 -> invalid file then valid -> calc.fdf -> Mn -> auto shell -> default alphas -> default output -> quit"
rm -rf hubbardu_runs
printf '4.7.1\ndoes_not_exist.fdf\nstructure.fdf\ncalc.fdf\nMn\n\n\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_contains "File not found" log_menu.txt
check_contains "Generated 13 run(s)" log_menu.txt
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
