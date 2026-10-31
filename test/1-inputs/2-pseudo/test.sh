#!/bin/bash

# --- Setup ---
# Smoke test for stb-pseudo (Pseudopotential Resolver, item 1.2)
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

check_missing() {
    if [ -e "$1" ]; then
        echo -e " ... ${RED}FAIL${NC} (file '$1' should NOT have been created)"
        FAIL=$((FAIL+1))
    else
        echo -e " ... ${GREEN}OK${NC} (file '$1' correctly not created)"
        PASS=$((PASS+1))
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
echo "--- Starting tester for STB-Pseudo (item 1.2) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null

cat > sic.fdf <<'EOF'
SystemName sic-test
SystemLabel sic
LatticeConstant 4.36 Ang
%block LatticeVectors
0.0 0.5 0.5
0.5 0.0 0.5
0.5 0.5 0.0
%endblock LatticeVectors
NumberOfSpecies 2
%block ChemicalSpeciesLabel
1 14 Si
2  6 C
%endblock ChemicalSpeciesLabel
NumberOfAtoms 2
AtomicCoordinatesFormat Fractional
%block AtomicCoordinatesAndAtomicSpecies
0.00 0.00 0.00 1
0.25 0.25 0.25 2
%endblock AtomicCoordinatesAndAtomicSpecies
EOF


# --- 2. --list-elements: browse a bundled bank, no structure needed ---
echo -e "\n--- Testing --list-elements (no structure/pp-path required) ---"
stb-pseudo --list-elements dojo --no-intro > log_list_elements.txt 2>&1
check_exit_code $? 0
check_contains "72 element(s) available" log_list_elements.txt
check_contains "Si" log_list_elements.txt

echo "Testing: --list-elements rejects an unknown bank"
stb-pseudo --list-elements not_a_bank --no-intro > log_list_bad.txt 2>&1
check_exit_code $? 2


# --- 3. -f/--file: resolve + copy from a real multi-species structure ---
echo -e "\n--- Testing -f/--file against a real structure (both species in 'dojo') ---"
rm -rf out_sic
stb-pseudo -f sic.fdf -p dojo -o out_sic --no-intro > log_sic.txt 2>&1
check_exit_code $? 0
check_contains "Species found in structure : Si, C" log_sic.txt
check_contains "Resolved : 2/2" log_sic.txt
check_success out_sic/Si.psml
check_success out_sic/C.psml


# --- 4. --species + --fallback-dir: fill a gap the primary bank can't cover ---
echo -e "\n--- Testing --species with --fallback-dir (At not in 'dojo', is in 'virtual_vault') ---"
rm -rf out_fallback
stb-pseudo --species Si At -p dojo --fallback-dir virtual_vault -o out_fallback --no-intro > log_fallback.txt 2>&1
check_exit_code $? 0
check_contains "Si  FOUND    (primary, .psml)" log_fallback.txt
check_contains "At  FOUND    (fallback, .psf)" log_fallback.txt
check_success out_fallback/Si.psml
check_success out_fallback/At.psf

echo "Testing: without --fallback-dir, the same request is reported (and exits) as incomplete"
stb-pseudo --species Si At -p dojo --dry-run --no-intro > log_no_fallback.txt 2>&1
check_exit_code $? 1
check_contains "At  MISSING" log_no_fallback.txt
check_contains "INCOMPLETE" log_no_fallback.txt


# --- 5. --dry-run: report only, nothing copied ---
echo -e "\n--- Testing --dry-run copies nothing ---"
rm -rf out_dryrun
stb-pseudo -f sic.fdf -p dojo -o out_dryrun --dry-run --no-intro > log_dryrun.txt 2>&1
check_exit_code $? 0
check_contains "dry-run: no files copied" log_dryrun.txt
check_missing out_dryrun/Si.psml


# --- 6. Error and robustness cases ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: neither -f nor --species given"
stb-pseudo -p dojo --no-intro > log_no_input.txt 2>&1
check_exit_code $? 1
check_contains "Either -f/--file or --species is required" log_no_input.txt

echo "Testing: both -f and --species given"
stb-pseudo -f sic.fdf --species Si -p dojo --no-intro > log_both_input.txt 2>&1
check_exit_code $? 1
check_contains "not both" log_both_input.txt

echo "Testing: -p/--pp-path missing"
stb-pseudo -f sic.fdf --no-intro > log_no_pp.txt 2>&1
check_exit_code $? 1
check_contains "\-p/--pp-path is required" log_no_pp.txt

echo "Testing: unknown bank/path given to -p/--pp-path"
stb-pseudo --species Si -p not_a_bank_or_path --no-intro > log_bad_pp.txt 2>&1
check_exit_code $? 1
check_contains "is not a recognized pseudopotential bank" log_bad_pp.txt

echo "Testing: nonexistent structure file is a clean error"
stb-pseudo -f does_not_exist.fdf -p dojo --no-intro > log_bad_fdf.txt 2>&1
check_exit_code $? 1

echo "Testing: --version"
stb-pseudo --version > log_version.txt 2>&1
check_contains "stb-pseudo" log_version.txt

echo "Testing: --help documents file/species/pp-path/fallback-dir/list-elements/dry-run/save-report"
stb-pseudo --help > log_help.txt 2>&1
check_contains "pp-path" log_help.txt
check_contains "fallback-dir" log_help.txt
check_contains "list-elements" log_help.txt
check_contains "dry-run" log_help.txt
check_contains "save-report" log_help.txt


# --- 7. --save-report and references.bib ---
echo -e "\n--- Testing --save-report and references.bib ---"

echo "Testing: --save-report writes stb_pseudo_report.txt matching the console output"
rm -f stb_pseudo_report.txt
stb-pseudo -f sic.fdf -p dojo -o out_sic --save-report --no-intro > log_save_report_console.txt 2>&1
check_success stb_pseudo_report.txt
check_contains "Resolved : 2/2" stb_pseudo_report.txt
check_contains "Report" log_save_report_console.txt

echo "Testing: references.bib always written, with SIESTA + bank citations"
rm -f references.bib
stb-pseudo --species Si -p dojo -o out_bib --no-intro > /dev/null 2>&1
check_success references.bib
check_contains "@article{Soler2002," references.bib
check_contains "@article{Garcia2020," references.bib
check_contains "@article{vanSetten2018," references.bib

echo "Testing: fallback bank's citation is also written when a fallback resolves an element"
rm -f references.bib
stb-pseudo --species At -p dojo --fallback-dir virtual_vault -o out_bib2 --no-intro > /dev/null 2>&1
check_contains "@misc{VirtualVault," references.bib


# --- 8. Interactive path (stb-suite, shortcut 1.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 1.2) ---"

echo "Testing: navigate 1.2 -> option 1 (structure file) -> invalid then valid file -> dojo -> defaults -> save report -> quit"
rm -f stb_pseudo_report.txt
printf '1.2\n1\ndoes_not_exist.fdf\nsic.fdf\ndojo\n\n\n\ny\n\n0\n' | timeout 20 stb-suite > log_menu_file.txt 2>&1
check_exit_code $? 0
check_contains "Species found in structure : Si, C" log_menu_file.txt
check_success stb_pseudo_report.txt

echo "Testing: navigate 1.2 -> option 2 (explicit species) -> dry-run -> no report -> quit"
printf '1.2\n2\nSi C\ndojo\n\n\nn\nn\n\n0\n' | timeout 20 stb-suite > log_menu_species.txt 2>&1
check_exit_code $? 0
check_contains "Species requested : Si, C" log_menu_species.txt

echo "Testing: navigate 1.2 -> option 3 (browse bank) -> dojo -> quit"
printf '1.2\n3\ndojo\n\n0\n' | timeout 20 stb-suite > log_menu_bank.txt 2>&1
check_exit_code $? 0
check_contains "72 element(s) available" log_menu_bank.txt


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
