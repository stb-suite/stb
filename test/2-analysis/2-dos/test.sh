#!/bin/bash

# --- Setup ---
# Smoke test for stb-dos (PDOS XML Parser, item 2.2)
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
    if grep -q -- "$1" "$2" 2>/dev/null; then
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
echo "--- Starting tester for STB-DOS (item 2.2) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/siesta.PDOS.xml" "$TEST_DIR/"
cp "$FIXTURE_DIR/spin2.PDOS.xml" "$TEST_DIR/"
cp "$FIXTURE_DIR/lmax4.PDOS.xml" "$TEST_DIR/"
cp "$FIXTURE_DIR/nspin_missing.PDOS.xml" "$TEST_DIR/"
cp "$FIXTURE_DIR/no_energy_values.PDOS.xml" "$TEST_DIR/"
cp "$FIXTURE_DIR/no_orbitals.PDOS.xml" "$TEST_DIR/"
cp "$FIXTURE_DIR/malformed.PDOS.xml" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Real fixture: default run (nspin=1, 182 orbitals, 14 atoms) ---
echo -e "\n--- Testing default run (real fixture, nspin=1, 182 orbitals) ---"
rm -rf dos_total.dat dos_per_atom dos_per_species
stb-dos siesta.PDOS.xml --no-intro > log_default.txt 2>&1
check_exit_code $? 0
check_success dos_total.dat
check_success dos_per_atom/Sn_1.dat
check_success dos_per_species/dos_Sn.dat
check_success dos_per_species/dos_O.dat

echo "Verifying dos_total.dat columns (nspin=1: no spin suffix)"
check_contains "#Energy(eV)" dos_total.dat
check_contains "s" dos_total.dat
check_contains "p" dos_total.dat
check_contains "d" dos_total.dat
check_not_contains "_up" dos_total.dat
check_not_contains "_down" dos_total.dat


# --- 3. --type filtering ---
echo -e "\n--- Testing --type filtering ---"
rm -rf dos_total.dat dos_per_atom dos_per_species
stb-dos siesta.PDOS.xml --type total --no-intro > log_type_total.txt 2>&1
check_exit_code $? 0
check_success dos_total.dat
if [ -d dos_per_atom ] || [ -d dos_per_species ]; then
    echo -e "   -> ${RED}Failed:${NC} --type total also created dos_per_atom/dos_per_species"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} --type total did not create dos_per_atom/dos_per_species"
    PASS=$((PASS+1))
fi


# --- 4. --shift modes ---
echo -e "\n--- Testing --shift fermi/manual/invalid ---"
stb-dos siesta.PDOS.xml --shift fermi --no-intro > log_shift_fermi.txt 2>&1
check_exit_code $? 0
check_contains "Using automatic Fermi energy shift" log_shift_fermi.txt

stb-dos siesta.PDOS.xml --shift 0.0 --no-intro > log_shift_zero.txt 2>&1
check_exit_code $? 0
check_contains "Using manual energy shift: 0.0 eV" log_shift_zero.txt

stb-dos siesta.PDOS.xml --shift notanumber --no-intro > log_shift_bad.txt 2>&1
check_exit_code $? 1
check_contains "Invalid shift value" log_shift_bad.txt


# --- 4b. --output-dir: writes into (and creates) a nested directory ---
echo -e "\n--- Testing --output-dir (nested, auto-created) ---"
rm -rf my_out dos_total.dat dos_per_atom dos_per_species
stb-dos siesta.PDOS.xml --type total --output-dir my_out/nested --no-intro > log_outdir.txt 2>&1
check_exit_code $? 0
check_success my_out/nested/dos_total.dat
check_contains "Saved Total DOS to my_out/nested/dos_total.dat" log_outdir.txt
if [ -e dos_total.dat ]; then
    echo -e "   -> ${RED}Failed:${NC} dos_total.dat was also written to CWD instead of only --output-dir"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} nothing was written to CWD, only to --output-dir"
    PASS=$((PASS+1))
fi
rm -rf my_out


# --- 5. --projection ml (detailed m-resolved columns) ---
echo -e "\n--- Testing --projection ml ---"
rm -rf dos_total.dat dos_per_atom dos_per_species
stb-dos siesta.PDOS.xml --type total --projection ml --no-intro > log_ml.txt 2>&1
check_exit_code $? 0
check_contains "px" dos_total.dat
check_contains "py" dos_total.dat
check_contains "pz" dos_total.dat


# --- 6. nspin=2: spin-polarized data must survive (regression test for the
# silent-data-loss bug: <nspin> was never read, so every orbital's data
# length mismatched num_energy_points and got silently dropped -- output
# used to contain only the Energy(eV) column). ---
echo -e "\n--- Testing nspin=2 (spin-polarized PDOS data is not silently dropped) ---"
rm -rf dos_total.dat dos_per_atom dos_per_species
stb-dos spin2.PDOS.xml --shift 0.0 --no-intro > log_spin2.txt 2>&1
check_exit_code $? 0
check_contains "Detected nspin=2" log_spin2.txt

echo "Verifying dos_total.dat has spin-resolved s_up/s_down columns with the right values"
check_contains "s_up" dos_total.dat
check_contains "s_down" dos_total.dat
check_contains "-1.000000E+00	  1.000000E+00	  5.000000E-01" dos_total.dat
check_contains "0.000000E+00	  2.000000E+00	  6.000000E-01" dos_total.dat
check_contains "1.000000E+00	  3.000000E+00	  7.000000E-01" dos_total.dat


# --- 7. <nspin> tag missing: falls back to nspin=1 instead of erroring ---
echo -e "\n--- Testing missing <nspin> tag (falls back to nspin=1) ---"
rm -rf dos_total.dat dos_per_atom dos_per_species
stb-dos nspin_missing.PDOS.xml --shift 0.0 --no-intro > log_nspin_missing.txt 2>&1
check_exit_code $? 0
check_contains "<nspin> tag not found" log_nspin_missing.txt
check_success dos_total.dat
check_not_contains "_up" dos_total.dat


# --- 8. l > 3 (g-orbitals) are counted and reported, not silently dropped ---
echo -e "\n--- Testing l>3 orbitals are excluded but reported ---"
rm -rf dos_total.dat dos_per_atom dos_per_species
stb-dos lmax4.PDOS.xml --shift 0.0 --no-intro > log_lmax4.txt 2>&1
check_exit_code $? 0
check_contains "skipped 1 orbital" log_lmax4.txt
check_success dos_total.dat
check_contains "#Energy(eV)" dos_total.dat


# --- 9. Error and robustness cases (all must now exit non-zero, not 0) ---
echo -e "\n--- Testing error cases ---"

echo "Testing: nonexistent input file"
stb-dos does_not_exist.xml --no-intro > log_missing_input.txt 2>&1
check_exit_code $? 1
check_contains "File not found" log_missing_input.txt

echo "Testing: malformed XML"
stb-dos malformed.PDOS.xml --no-intro > log_malformed.txt 2>&1
check_exit_code $? 1
check_contains "Error parsing XML file" log_malformed.txt

echo "Testing: no <energy_values> tag"
stb-dos no_energy_values.PDOS.xml --no-intro > log_no_energy.txt 2>&1
check_exit_code $? 1
check_contains "energy_values" log_no_energy.txt

echo "Testing: no <orbital> tags"
stb-dos no_orbitals.PDOS.xml --no-intro > log_no_orbitals.txt 2>&1
check_exit_code $? 1
check_contains "No <orbital> tags found" log_no_orbitals.txt

echo "Testing: --version"
stb-dos --version > log_version.txt 2>&1
check_contains "stb-dos" log_version.txt


# --- 10. Interactive path (stb-suite, shortcut 2.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 2.2) ---"

echo "Testing: navigate 2.2 -> siesta.PDOS.xml -> defaults"
rm -rf dos_total.dat dos_per_atom dos_per_species
printf '2.2\nsiesta.PDOS.xml\n\n\n\n\n' | stb-suite > log_menu.txt 2>&1
check_success dos_total.dat

echo "Testing: interactive path surfaces a tool failure (malformed XML -> non-zero exit -> run_tool reports it)"
printf '2.2\nmalformed.PDOS.xml\n\n\n\n\n' | stb-suite > log_menu_fail.txt 2>&1
check_contains "Error running stb-dos" log_menu_fail.txt

echo "Testing: interactive path forwards a custom output directory"
rm -rf menu_out
printf '2.2\nsiesta.PDOS.xml\n\n\n\nmenu_out\n' | stb-suite > log_menu_outdir.txt 2>&1
check_success menu_out/dos_total.dat


popd > /dev/null

# --- 11. Summary ---
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
