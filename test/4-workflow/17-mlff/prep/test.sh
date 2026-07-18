#!/bin/bash

# --- Setup ---
# Smoke test for stb-mlff (Custom ML Force Field Prep, item 4.17.1)
#
# calc.fdf is a real O2-in-vacuum SIESTA input (SZ basis, same reference
# used to verify this whole workflow live: 8 rattled configs -> real SIESTA
# -> stb-mlffAnalysis -> a genuine fine-tuned MACE model). Pseudopotential
# content doesn't matter here (prep only copies the file, never reads it),
# so a placeholder is enough -- same convention as test/4-workflow/3-cohesive/prep.
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
echo "--- Starting tester for stb-mlff (item 4.17.1) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/calc.fdf" "$TEST_DIR/"
echo "dummy pseudopotential content" > "$TEST_DIR/O.psf"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Basic rattling run ---
echo -e "\n--- Testing a basic rattling run (--n-configs 6 --stdev 0.05 0.1) ---"
stb-mlff --file calc.fdf --n-configs 6 --stdev 0.05 0.1 --pseudo-dir . --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "Generated 6 rattled configuration(s)" log_basic.txt
for i in 001 002 003 004 005 006; do
    check_success "mlff_config_$i/calc.fdf"
    check_success "mlff_config_$i/O.psf"
done

echo "Testing: each config has the same species/basis/SCF settings but different positions"
check_contains "PAO.BasisSize       SZ" mlff_config_001/calc.fdf
python3 -c "
import numpy as np, sys, re
def read_positions(path):
    text = open(path).read()
    block = re.search(r'%block AtomicCoordinatesAndAtomicSpecies\n(.*?)%endblock', text, re.S).group(1)
    return np.array([[float(x) for x in line.split()[:3]] for line in block.strip().splitlines()])
p1 = read_positions('mlff_config_001/calc.fdf')
p2 = read_positions('mlff_config_002/calc.fdf')
ref = read_positions('calc.fdf')
sys.exit(0 if not np.allclose(p1, p2) and not np.allclose(p1, ref) else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} rattled configs have distinct, displaced positions"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} rattled configs are not distinct from each other/the reference"
    FAIL=$((FAIL+1))
fi


# --- 3. --n-configs 0 with no --from-aimd (expects failure) ---
echo -e "\n--- Testing --n-configs 0 with no --from-aimd (expects failure) ---"
stb-mlff --file calc.fdf --n-configs 0 --pseudo-dir . --no-intro > log_empty.txt 2>&1
check_exit_code $? 1
check_contains "No configurations were generated" log_empty.txt


# --- 4. Missing pseudopotential (expects failure) ---
echo -e "\n--- Testing a missing pseudopotential (expects failure) ---"
mkdir -p empty_pp_dir
stb-mlff --file calc.fdf --n-configs 2 --pseudo-dir empty_pp_dir --no-intro > log_missingpp.txt 2>&1
check_exit_code $? 1
check_contains "Missing pseudopotential" log_missingpp.txt


# --- 5. --save-report ---
echo -e "\n--- Testing --save-report ---"
rm -f stb_mlff_report.txt
stb-mlff --file calc.fdf --n-configs 2 --pseudo-dir . --save-report --no-intro > log_savereport.txt 2>&1
check_success stb_mlff_report.txt
check_contains "STB-MLFF REPORT" stb_mlff_report.txt
rm -f stb_mlff_report.txt


# --- 6. Robustness / errors ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: missing --file"
stb-mlff --no-intro > log_missing_file.txt 2>&1
check_exit_code $? 2

echo "Testing: nonexistent reference file"
stb-mlff --file nope.fdf --no-intro > log_missing_ref.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_ref.txt

echo "Testing: --version"
stb-mlff --version > log_version.txt 2>&1
check_contains "stb-mlff" log_version.txt

echo "Testing: --help documents --stdev, --from-aimd, --n-samples"
stb-mlff --help > log_help.txt 2>&1
check_contains "stdev" log_help.txt
check_contains "from-aimd" log_help.txt
check_contains "n-samples" log_help.txt


# --- 7. Interactive path (stb-suite, shortcut 4.17.1) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.17.1) ---"
rm -rf mlff_config_001 mlff_config_002 mlff_config_003
printf '4.17.1\ncalc.fdf\n\n3\n0.05\nn\n\nn\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success mlff_config_001/calc.fdf
check_success mlff_config_003/calc.fdf


popd > /dev/null

# --- 8. Summary ---
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
