#!/bin/bash

# --- Setup ---
# Smoke test for stb-convergence (Convergence Test Prep, item 4.5.1)
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
echo "--- Starting tester for STB-Convergence prep (item 4.5.1) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/structure.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/calc.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Mesh.CutOff sweep ---
echo -e "\n--- Testing --parameter meshcutoff ---"
rm -rf convergence_runs
stb-convergence -s structure.fdf -c calc.fdf --parameter meshcutoff --min 150 --max 300 --step 50 --no-intro > log_meshcutoff.txt 2>&1
check_contains "Generated 4 run(s)" log_meshcutoff.txt
check_success convergence_runs/convergence_meshcutoff_250.0000/calc.fdf
check_contains "Mesh.CutOff           250.0000  Ry" convergence_runs/convergence_meshcutoff_250.0000/calc.fdf
check_success convergence_runs/convergence_meshcutoff_250.0000/structure.fdf


# --- 3. PAO.EnergyShift sweep ---
echo -e "\n--- Testing --parameter energyshift ---"
rm -rf convergence_runs
stb-convergence -s structure.fdf -c calc.fdf --parameter energyshift --min 0.01 --max 0.05 --step 0.02 --no-intro > log_energyshift.txt 2>&1
check_contains "Generated 3 run(s)" log_energyshift.txt
check_contains "PAO.EnergyShift     0.0300 Ry" convergence_runs/convergence_energyshift_0.0300/calc.fdf


# --- 4. K-grid density sweep (cross-checked against stb-kgrid directly) ---
echo -e "\n--- Testing --parameter kgrid ---"
rm -rf convergence_runs
stb-convergence -s structure.fdf -c calc.fdf --parameter kgrid --min 0.10 --max 0.20 --step 0.05 --no-intro > log_kgrid.txt 2>&1
check_contains "Generated 3 run(s)" log_kgrid.txt
check_contains "kgrid.MonkhorstPack   \[8  8  8\]" convergence_runs/convergence_kgrid_0.1500/calc.fdf

stb-kgrid --file structure.fdf --type fdf --density 0.15 --no-intro > log_kgrid_crosscheck.txt 2>&1
check_contains "Suggested Monkhorst-Pack grid: 8 8 8" log_kgrid_crosscheck.txt


# --- 5. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: calc.fdf missing the swept tag"
cat > calc_notag.fdf << 'EOF'
SystemLabel siesta
%include structure.fdf
EOF
stb-convergence -s structure.fdf -c calc_notag.fdf --parameter meshcutoff --min 150 --max 300 --step 50 --no-intro > log_notag.txt 2>&1
check_exit_code $? 1
check_contains "Could not find a 'meshcutoff' tag" log_notag.txt

echo "Testing: --step 0"
stb-convergence -s structure.fdf -c calc.fdf --parameter meshcutoff --min 150 --max 300 --step 0 --no-intro > log_zero_step.txt 2>&1
check_exit_code $? 2

echo "Testing: --max <= --min"
stb-convergence -s structure.fdf -c calc.fdf --parameter meshcutoff --min 300 --max 150 --step 50 --no-intro > log_bad_range.txt 2>&1
check_exit_code $? 2

echo "Testing: missing structure file"
stb-convergence -s does_not_exist.fdf -c calc.fdf --parameter meshcutoff --min 150 --max 300 --step 50 --no-intro > log_missing.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing.txt

echo "Testing: missing required args"
stb-convergence -s structure.fdf -c calc.fdf --no-intro > log_missing_args.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-convergence --version > log_version.txt 2>&1
check_contains "stb-convergence" log_version.txt

echo "Testing: --help documents the 3 parameters"
stb-convergence --help > log_help.txt 2>&1
check_contains "meshcutoff" log_help.txt
check_contains "kgrid" log_help.txt


# --- 6. Interactive path (stb-suite, shortcut 3.5.1) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 3.5.1) ---"

echo "Testing: navigate 3.5.1 -> invalid file then valid -> calc.fdf -> meshcutoff -> 150/300/50 -> default output -> quit"
rm -rf convergence_runs
printf '4.5.1\ndoes_not_exist.fdf\nstructure.fdf\ncalc.fdf\n1\n150\n300\n50\n\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_contains "File not found" log_menu.txt
check_contains "Generated 4 run(s)" log_menu.txt
check_success convergence_runs/convergence_meshcutoff_250.0000/calc.fdf


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
