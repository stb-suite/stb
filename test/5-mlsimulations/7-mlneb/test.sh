#!/bin/bash

# --- Setup ---
# Smoke test for stb-mlneb (ML NEB, item 5.7)
#
# Needs the optional 'ml' extra (pip install stb_suite[ml] -- PyTorch +
# mace-torch). The whole file is skipped with a clear message if `mace`
# isn't importable, same gating pattern as the other test/5-mlsimulations/*
# tools.
#
# vac_initial.fdf/vac_final.fdf: a real single-atom vacancy-migration hop
# in 7-atom bulk Si (see each fixture's own header comment for the exact
# site-to-site move). Verified live: the tool reports a ~0.59 eV barrier
# (physically the right order of magnitude for Si vacancy migration,
# commonly cited ~0.2-0.5 eV in the DFT literature) and a reaction energy
# of essentially 0 eV (expected by symmetry -- both endpoints are the same
# vacancy-in-bulk-Si situation, just at two symmetric-equivalent sites).
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


echo "--- Starting tester for stb-mlneb (item 5.7) ---"

if ! python3 -c "import mace" 2>/dev/null; then
    echo -e "${YELLOW}Skipped entirely:${NC} the optional 'ml' extra is not installed."
    echo "Install with: pip install stb_suite[ml]  (then re-run this test)"
    exit 0
fi

rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/vac_initial.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/vac_final.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 1. Basic run: vacancy migration barrier ---
echo -e "\n--- Testing a basic NEB run (vacancy migration in Si) ---"
stb-mlneb --initial vac_initial.fdf --final vac_final.fdf --n-images 7 \
    --save-images --save-report --no-intro -o vac_neb > log_neb.txt 2>&1
check_exit_code $? 0
check_contains "NEB PIPELINE (MACE)" log_neb.txt
check_contains "Barrier (forward)" log_neb.txt
check_contains "converged" log_neb.txt
check_success vac_neb/neb_preview.png
check_success vac_neb/neb_path.xyz
check_success vac_neb/image_00.fdf
check_success vac_neb/image_06.fdf
check_success vac_neb/stb_mlneb_report.txt
check_contains "STB-MLNEB REPORT" vac_neb/stb_mlneb_report.txt

echo "Testing: barrier is positive and in a physically sensible range (0.05-3.0 eV)"
python3 -c "
import re, sys
text = open('vac_neb/stb_mlneb_report.txt').read()
m = re.search(r'Barrier\s*:\s*([\d.]+) eV', text)
if not m:
    sys.exit(1)
barrier = float(m.group(1))
sys.exit(0 if 0.05 < barrier < 3.0 else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} barrier is positive and physically sensible"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} barrier is outside the expected physical range"
    FAIL=$((FAIL+1))
fi

echo "Testing: reaction energy is near 0 (symmetric endpoints, expected by construction)"
python3 -c "
import re, sys
text = open('vac_neb/stb_mlneb_report.txt').read()
m = re.search(r'Reaction energy\s*:\s*(-?[\d.]+) eV', text)
if not m:
    sys.exit(1)
de = float(m.group(1))
sys.exit(0 if abs(de) < 0.05 else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} reaction energy is near 0, as expected for symmetric endpoints"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} reaction energy is unexpectedly far from 0"
    FAIL=$((FAIL+1))
fi
rm -rf vac_neb


# --- 2. --freeze-substrate runs without crashing ---
echo -e "\n--- Testing --freeze-substrate ---"
stb-mlneb --initial vac_initial.fdf --final vac_final.fdf --n-images 5 \
    --freeze-substrate --no-intro -o vac_frozen > log_frozen.txt 2>&1
check_exit_code $? 0
check_contains "Freezing" log_frozen.txt
rm -rf vac_frozen


# --- 3. Composition mismatch is rejected ---
echo -e "\n--- Testing composition mismatch rejection ---"
stb-mlneb --initial vac_initial.fdf --final ../../3-mlelastic/si8.fdf --no-intro > log_mismatch.txt 2>&1
check_exit_code $? 1
check_contains "different composition" log_mismatch.txt


# --- 4. Error and robustness cases ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: missing --initial"
stb-mlneb --final vac_final.fdf --no-intro > log_missing_initial_arg.txt 2>&1
check_exit_code $? 2

echo "Testing: nonexistent initial file"
stb-mlneb --initial nope.fdf --final vac_final.fdf --no-intro > log_missing_initial.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_initial.txt

echo "Testing: --custom-model with a nonexistent file"
stb-mlneb --initial vac_initial.fdf --final vac_final.fdf --custom-model does_not_exist.model \
    --no-intro > log_custommodel.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_custommodel.txt

echo "Testing: --version"
stb-mlneb --version > log_version.txt 2>&1
check_contains "stb-mlneb" log_version.txt

echo "Testing: --help documents --initial, --final, --n-images, --custom-model, --freeze-substrate"
stb-mlneb --help > log_help.txt 2>&1
check_contains "initial" log_help.txt
check_contains "final" log_help.txt
check_contains "n-images" log_help.txt
check_contains "custom-model" log_help.txt
check_contains "freeze-substrate" log_help.txt
check_contains "idpp" log_help.txt


# --- 5. Interactive path (stb-suite, shortcut 5.7) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 5.7) ---"
rm -rf interactive_out
printf '5.7\nvac_initial.fdf\nvac_final.fdf\n\nsmall\n\n5\nn\ninteractive_out\nn\nn\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success interactive_out/neb_preview.png
check_success interactive_out/neb_path.xyz


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
