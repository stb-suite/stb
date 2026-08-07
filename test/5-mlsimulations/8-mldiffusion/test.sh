#!/bin/bash

# --- Setup ---
# Smoke test for stb-mldiffusion (ML Vacancy Diffusion, item 5.8)
#
# Needs the optional 'ml' extra (pip install stb_suite[ml] -- PyTorch +
# mace-torch). The whole file is skipped with a clear message if `mace`
# isn't importable, same gating pattern as the other test/5-mlsimulations/*
# tools.
#
# si8.fdf: 8-atom bulk Si (diamond cubic, a=5.43 Ang), same fixture as
# test/5-mlsimulations/3-mlelastic. Verified live: removing atom index 1
# (a site with exactly 4 tetrahedral nearest neighbors, the correct Si
# coordination number) auto-detects all 4 as candidate hops, and ALL FOUR
# give the exact same ~0.59 eV barrier -- expected by cubic symmetry, since
# all 4 nearest-neighbor directions in diamond Si are crystallographically
# equivalent. A strong physical correctness signal for both the neighbor
# auto-detection and the barrier computation.
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


echo "--- Starting tester for stb-mldiffusion (item 5.8) ---"

if ! python3 -c "import mace" 2>/dev/null; then
    echo -e "${YELLOW}Skipped entirely:${NC} the optional 'ml' extra is not installed."
    echo "Install with: pip install stb_suite[ml]  (then re-run this test)"
    exit 0
fi

rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/si8.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 1. Basic run: vacancy at index 1, 4 tetrahedral neighbors ---
echo -e "\n--- Testing a basic scan (vacancy at index 1, --n-images 5) ---"
stb-mldiffusion --file si8.fdf --vacancy-index 1 --n-images 5 \
    --save-report --no-intro -o vac_scan > log_scan.txt 2>&1
check_exit_code $? 0
check_contains "NEIGHBOR SHELL DETECTION" log_scan.txt
check_contains "MIGRATION-BARRIER SCAN" log_scan.txt
check_contains "Lowest barrier" log_scan.txt
check_success vac_scan/migration_barriers.png
check_success vac_scan/lowest_barrier_preview.png
check_success vac_scan/lowest_barrier_path.xyz
check_success vac_scan/stb_mldiffusion_report.txt
check_contains "STB-MLDIFFUSION REPORT" vac_scan/stb_mldiffusion_report.txt

echo "Testing: exactly 4 tetrahedral neighbors auto-detected (correct Si coordination number)"
python3 -c "
import re, sys
text = open('vac_scan/stb_mldiffusion_report.txt').read()
m = re.search(r'Candidate hop target\(s\).*: \[([\d, ]+)\]', text)
if not m:
    sys.exit(1)
targets = [int(x) for x in m.group(1).split(',')]
sys.exit(0 if len(targets) == 4 else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} exactly 4 neighbors detected (correct Si coordination number)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} did not detect exactly 4 neighbors"
    FAIL=$((FAIL+1))
fi

echo "Testing: all 4 symmetric-equivalent hops give the same barrier (cubic symmetry)"
python3 -c "
import re, sys
text = open('vac_scan/stb_mldiffusion_report.txt').read()
barriers = [float(x) for x in re.findall(r'^\s*\d+\s+([\d.]+)\s+', text, re.M)]
if len(barriers) < 4:
    sys.exit(1)
sys.exit(0 if max(barriers) - min(barriers) < 0.01 else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} all symmetric-equivalent hops give the same barrier"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} symmetric-equivalent hops gave different barriers"
    FAIL=$((FAIL+1))
fi

echo "Testing: reported barrier is positive and physically sensible (0.05-3.0 eV)"
python3 -c "
import re, sys
text = open('vac_scan/stb_mldiffusion_report.txt').read()
m = re.search(r'Lowest barrier\s*:\s*([\d.]+) eV', text)
if not m:
    sys.exit(1)
barrier = float(m.group(1))
sys.exit(0 if 0.05 < barrier < 3.0 else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} lowest barrier is positive and physically sensible"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} lowest barrier is outside the expected physical range"
    FAIL=$((FAIL+1))
fi
rm -rf vac_scan


# --- 2. Error and robustness cases ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: --vacancy-index out of range"
stb-mldiffusion --file si8.fdf --vacancy-index 99 --no-intro > log_outofrange.txt 2>&1
check_exit_code $? 1
check_contains "out of range" log_outofrange.txt

echo "Testing: missing --file"
stb-mldiffusion --vacancy-index 0 --no-intro > log_missing_file.txt 2>&1
check_exit_code $? 2

echo "Testing: nonexistent input file"
stb-mldiffusion --file nope.fdf --vacancy-index 0 --no-intro > log_missing_input.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_input.txt

echo "Testing: --custom-model with a nonexistent file"
stb-mldiffusion --file si8.fdf --vacancy-index 0 --custom-model does_not_exist.model \
    --no-intro > log_custommodel.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_custommodel.txt

echo "Testing: --version"
stb-mldiffusion --version > log_version.txt 2>&1
check_contains "stb-mldiffusion" log_version.txt

echo "Testing: --help documents --vacancy-index, --shell-tolerance, --custom-model"
stb-mldiffusion --help > log_help.txt 2>&1
check_contains "vacancy-index" log_help.txt
check_contains "shell-tolerance" log_help.txt
check_contains "custom-model" log_help.txt
check_contains "freeze-substrate" log_help.txt


# --- 3. Interactive path (stb-suite, shortcut 5.8) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 5.8) ---"
rm -rf interactive_out
printf '5.8\nsi8.fdf\n1\n\nsmall\n\n5\ninteractive_out\nn\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success interactive_out/migration_barriers.png


popd > /dev/null

# --- 4. Summary ---
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
