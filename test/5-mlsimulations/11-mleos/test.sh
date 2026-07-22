#!/bin/bash

# --- Setup ---
# Smoke test for stb-mleos (ML Equation of State, item 5.11)
#
# Needs the optional 'ml' extra (pip install stb_suite[ml] -- PyTorch +
# mace-torch). The whole file is skipped with a clear message if `mace`
# isn't importable, same gating pattern as the other test/5-mlsimulations/*
# tools.
#
# si8.fdf: 8-atom bulk Si (diamond cubic, a=5.43 Ang), same fixture as
# test/5-mlsimulations/3-mlelastic and .../6-mlconvergence. Verified live:
# default settings give V0=163.19 Ang^3, E0=-42.967 eV (matches
# test/5-mlsimulations/4-mlsearch's independently-found energy minimum for
# the same structure, ~-42.966 eV), B0=71.22 GPa, B0'=4.32, R^2=0.999993 --
# and B0 agrees with stb-mlelastic's own stress-strain-derived Bulk Modulus
# on the exact same structure (70.41 GPa) to ~1%, an independent-method
# cross-check on the MACE potential itself.
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


echo "--- Starting tester for stb-mleos (item 5.11) ---"

if ! python3 -c "import mace" 2>/dev/null; then
    echo -e "${YELLOW}Skipped entirely:${NC} the optional 'ml' extra is not installed."
    echo "Install with: pip install stb_suite[ml]  (then re-run this test)"
    exit 0
fi

rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/si8.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/../3-mlelastic/graphene.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 1. Basic run (foundation model, default settings) ---
echo -e "\n--- Testing a basic volume scan + EOS fit (bulk Si, default settings) ---"
stb-mleos --file si8.fdf --save-data --save-report --no-intro -o si_basic > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "READING STRUCTURE" log_basic.txt
check_contains "REFERENCE PRE-RELAX" log_basic.txt
check_contains "VOLUME SCAN (MACE)" log_basic.txt
check_contains "EQUATION-OF-STATE FIT" log_basic.txt
check_contains "Bulk modulus" log_basic.txt
check_success si_basic/eos_curve.png
check_success si_basic/stb_mleos_report.txt
check_contains "STB-MLEOS REPORT" si_basic/stb_mleos_report.txt

echo "Testing: bulk modulus is a physically sensible value for Si (10-500 GPa)"
python3 -c "
import re, sys
text = re.sub(r'\x1b\[[0-9;]*m', '', open('log_basic.txt').read())
m = re.search(r'Bulk modulus\s*\(B0\) : ([\-0-9.]+) GPa', text)
sys.exit(0 if m and 10.0 < float(m.group(1)) < 500.0 else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} bulk modulus is a physically sensible Si value"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} bulk modulus was not physically sensible"
    FAIL=$((FAIL+1))
fi

echo "Testing: fit quality R^2 is close to 1 (good EOS fit)"
python3 -c "
import re, sys
text = re.sub(r'\x1b\[[0-9;]*m', '', open('log_basic.txt').read())
m = re.search(r'Fit quality\s*\(R\^2\): ([0-9.]+)', text)
sys.exit(0 if m and float(m.group(1)) > 0.99 else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} EOS fit quality (R^2) is close to 1"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} EOS fit quality was not close to 1"
    FAIL=$((FAIL+1))
fi
rm -rf si_basic


# --- 2. Alternate EOS form (--eos vinet) ---
echo -e "\n--- Testing --eos vinet ---"
stb-mleos --file si8.fdf --eos vinet --no-intro -o si_vinet > log_vinet.txt 2>&1
check_exit_code $? 0
check_contains "EOS form          : vinet" log_vinet.txt
check_contains "Pressure derivative" log_vinet.txt
rm -rf si_vinet


# --- 3. --eos birch_murnaghan alias normalizes to ASE's own spelling ---
echo -e "\n--- Testing --eos birch_murnaghan (alias for birchmurnaghan) ---"
stb-mleos --file si8.fdf --eos birch_murnaghan --no-intro -o si_alias > log_alias.txt 2>&1
check_exit_code $? 0
check_contains "EOS form          : birchmurnaghan" log_alias.txt
rm -rf si_alias


# --- 4. --n-volumes below the minimum is rejected ---
echo -e "\n--- Testing --n-volumes 3 is rejected (needs >= 5) ---"
stb-mleos --file si8.fdf --n-volumes 3 --no-intro -o si_toofew > log_toofew.txt 2>&1
check_exit_code $? 2
rm -rf si_toofew


# --- 5. Vacuum-padded structure is rejected (bulk-only) ---
echo -e "\n--- Testing a vacuum-padded (2D) structure is rejected ---"
stb-mleos --file graphene.fdf --no-intro -o si_graphene > log_graphene.txt 2>&1
check_exit_code $? 1
check_contains "bulk (3D periodic) only" log_graphene.txt
rm -rf si_graphene


# --- 6. Fine-tuning a quick Si model to test the foundation-model comparison with ---
echo -e "\n--- Fine-tuning a quick Si model to test the foundation-model comparison with ---"
PHONONS_FIXTURES="$FIXTURE_DIR/../2-mlphonons/si_mlff_fixtures"
if [ -d "$PHONONS_FIXTURES" ]; then
    cp -r "$PHONONS_FIXTURES" mlff_config_src
    mv mlff_config_src/mlff_config_* .
    rmdir mlff_config_src
    stb-mlffAnalysis --path "mlff_config_*" --epochs 3 --batch-size 4 --device cpu \
        --name si_quick_model --skip-foundation-comparison --no-intro > log_quicktrain.txt 2>&1
    QUICK_MODEL=$(find . -maxdepth 1 -name "si_quick_model*.model" | sort | tail -1)
    rm -rf mlff_config_001 mlff_config_002 mlff_config_003 mlff_config_004 \
           mlff_config_005 mlff_config_006 mlff_config_007 mlff_config_008
else
    QUICK_MODEL=""
fi

if [ -n "$QUICK_MODEL" ]; then
    echo "Using model: $QUICK_MODEL"

    echo -e "\n--- Testing --custom-model with foundation comparison (default on) ---"
    stb-mleos --file si8.fdf --custom-model "$QUICK_MODEL" --n-volumes 5 \
        --no-intro -o si_compare > log_compare.txt 2>&1
    check_exit_code $? 0
    check_contains "Comparison        : also fitting with foundation model" log_compare.txt
    check_contains "FOUNDATION MODEL COMPARISON" log_compare.txt
    check_contains "(foundation)" log_compare.txt
    check_success si_compare/eos_curve.png
    rm -rf si_compare

    echo "Testing: --skip-foundation-comparison disables the second run"
    stb-mleos --file si8.fdf --custom-model "$QUICK_MODEL" --n-volumes 5 \
        --skip-foundation-comparison --no-intro -o si_nocompare > log_nocompare.txt 2>&1
    check_exit_code $? 0
    if grep -q "FOUNDATION MODEL COMPARISON" log_nocompare.txt; then
        echo -e "   -> ${RED}Failed:${NC} foundation comparison ran despite --skip-foundation-comparison"
        FAIL=$((FAIL+1))
    else
        echo -e "   -> ${GREEN}Verified:${NC} foundation comparison correctly skipped"
        PASS=$((PASS+1))
    fi
    rm -rf si_nocompare
else
    echo -e "${YELLOW}SKIPPED${NC}: foundation-comparison tests (no quick model available)."
fi


# --- 7. Error and robustness cases ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: missing --file"
stb-mleos --no-intro > log_missing_file.txt 2>&1
check_exit_code $? 2

echo "Testing: nonexistent input file"
stb-mleos --file nope.fdf --no-intro > log_missing_input.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_input.txt

echo "Testing: --custom-model with a nonexistent file"
stb-mleos --file si8.fdf --custom-model does_not_exist.model --no-intro > log_custommodel.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_custommodel.txt

echo "Testing: --version"
stb-mleos --version > log_version.txt 2>&1
check_contains "stb-mleos" log_version.txt

echo "Testing: --help documents --file, --eos, --n-volumes, --strain-range, --custom-model"
stb-mleos --help > log_help.txt 2>&1
check_contains "\-\-file" log_help.txt
check_contains "eos" log_help.txt
check_contains "n-volumes" log_help.txt
check_contains "strain-range" log_help.txt
check_contains "custom-model" log_help.txt
check_contains "skip-foundation-comparison" log_help.txt


# --- 8. Interactive path (stb-suite, shortcut 5.11) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 5.11) ---"
rm -rf interactive_out
printf '5.11\nsi8.fdf\n\n\n\n\n\ninteractive_out\nn\nn\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success interactive_out/eos_curve.png


popd > /dev/null

# --- 9. Summary ---
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
