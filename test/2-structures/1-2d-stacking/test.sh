#!/bin/bash

# --- Setup ---
# Smoke test for stb-2Dstacking (2D Monolayer Stacker, item 2.1)
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
echo "--- Starting tester for STB-2Dstacking (item 2.1) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/graphene.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/hbn.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null

cat > malformed.fdf << 'EOF'
NumberOfSpecies    1
NumberofAtoms      1
LatticeConstant 1.0 Ang
AtomicCoordinatesFormat  Fractional
%block LatticeVectors
 5.0 0.0 0.0
 0.0 5.0 0.0
 0.0 0.0 5.0
%endblock LatticeVectors
%block AtomicCoordinatesAndAtomicSpecies
 0.0 0.0 0.0 1
%endblock AtomicCoordinatesAndAtomicSpecies
EOF


# --- 2. Basic stacking (graphene + h-BN, ~1.8% lattice mismatch) ---
echo -e "\n--- Testing basic stacking (mismatched pair: graphene + h-BN) ---"

rm -f stacked_structure.fdf symmetry_report.txt
stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf --no-intro > log_basic.txt 2>&1
check_contains "Selected Match ID 0" log_basic.txt
check_contains "BC2N" log_basic.txt
check_success stacked_structure.fdf
check_success symmetry_report.txt
check_contains "%block LatticeVectors" stacked_structure.fdf
check_contains "%block AtomicCoordinatesAndAtomicSpecies" stacked_structure.fdf
mv stacked_structure.fdf stacked_basic.fdf


# --- 3. Identical-layer bilayer (perfect lattice match, 0% strain) ---
echo -e "\n--- Testing identical-layer bilayer (expected 0% strain) ---"

rm -f stacked_structure.fdf symmetry_report.txt
stb-2Dstacking -l1 graphene.fdf -l2 graphene.fdf --no-intro > log_identical.txt 2>&1
check_contains "Max Applied Linear Strain" log_identical.txt
check_contains "0.00%" log_identical.txt
check_success stacked_structure.fdf
mv stacked_structure.fdf stacked_identical.fdf


# --- 4. --batch_sym (hexagonal high-symmetry stackings) ---
echo -e "\n--- Testing --batch_sym (hexagonal detection -> 6 configurations) ---"

rm -f stacked_structure_*.fdf symmetry_report_*.txt
stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf --batch_sym --no-intro > log_batch_sym.txt 2>&1
check_contains "Hexagonal lattice detected" log_batch_sym.txt
for name in AA Bridge_X Bridge_Y Center Hex_AB Hex_BA; do
    check_success "stacked_structure_${name}.fdf"
    check_success "symmetry_report_${name}.txt"
done

echo "Testing: --batch_sym ignores custom shifts with a warning"
stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf --batch_sym -tx 0.1 --no-intro > log_batch_sym_shift.txt 2>&1
check_contains "Custom shifts .* are IGNORED" log_batch_sym_shift.txt
rm -f stacked_structure_*.fdf symmetry_report_*.txt


# --- 5. --gap_range (energy-curve mode) ---
echo -e "\n--- Testing --gap_range (energy-curve mode, 3 points) ---"

rm -f stacked_structure_*.fdf symmetry_report_*.txt
stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf --gap_range 3.0 4.0 3 --no-intro > log_gap_range.txt 2>&1
check_contains "Generating 3 distance points" log_gap_range.txt
check_success stacked_structure_Custom_gap_3.00.fdf
check_success stacked_structure_Custom_gap_3.50.fdf
check_success stacked_structure_Custom_gap_4.00.fdf
rm -f stacked_structure_*.fdf symmetry_report_*.txt


# --- 6. Geometric controls (twist, shifts, strain mode, vacuum, custom output) ---
echo -e "\n--- Testing geometric controls ---"

echo "Testing: -t/--twist"
rm -f stacked_structure.fdf
stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf -t 5.0 --no-intro > log_twist.txt 2>&1
check_contains "Applying initial twist angle of 5.0" log_twist.txt
rm -f stacked_structure.fdf symmetry_report.txt

echo "Testing: -sm bottom / -sm sym change the applied strain"
stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf -sm top --no-intro > log_sm_top.txt 2>&1
stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf -sm bottom --no-intro > log_sm_bottom.txt 2>&1
stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf -sm sym --no-intro > log_sm_sym.txt 2>&1
STRAIN_TOP=$(grep "Max Applied Linear Strain" log_sm_top.txt)
STRAIN_BOTTOM=$(grep "Max Applied Linear Strain" log_sm_bottom.txt)
STRAIN_SYM=$(grep "Max Applied Linear Strain" log_sm_sym.txt)
if [ "$STRAIN_TOP" != "$STRAIN_BOTTOM" ] && [ "$STRAIN_TOP" != "$STRAIN_SYM" ]; then
    echo -e "   -> ${GREEN}Verified:${NC} top/bottom/sym strain modes give distinct strain values"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} strain modes did not produce distinct values"
    FAIL=$((FAIL+1))
fi
rm -f stacked_structure.fdf symmetry_report.txt

echo "Testing: --vacuum custom value"
stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf --vacuum 15.0 --no-intro > log_vacuum.txt 2>&1
check_contains "Vacuum Corrected" log_vacuum.txt
rm -f stacked_structure.fdf symmetry_report.txt

echo "Testing: -o/--output and --sym_out custom filenames"
stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf -o my_hetero.fdf --sym_out my_symmetry.txt --no-intro > log_custom_names.txt 2>&1
check_success my_hetero.fdf
check_success my_symmetry.txt
rm -f my_hetero.fdf my_symmetry.txt

echo "Testing: -sp/--symprec custom tolerance"
stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf -sp 0.5 --no-intro > log_symprec.txt 2>&1
check_contains "Tolerance: 0.5" log_symprec.txt
rm -f stacked_structure.fdf symmetry_report.txt


# --- 7. Match selection (-id, -i) ---
echo -e "\n--- Testing match selection ---"

echo "Testing: -id/--match_id direct selection"
stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf -id 0 --no-intro > log_match_id.txt 2>&1
check_contains "Selected Match ID 0" log_match_id.txt
rm -f stacked_structure.fdf symmetry_report.txt

echo "Testing: -id out of range (genuinely, beyond the ~4481 ZSL matches for this pair)"
stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf -id 999999 --no-intro > log_match_id_bad.txt 2>&1
check_exit_code $? 1
check_contains "is out of range" log_match_id_bad.txt

echo "Testing: -i/--interactive with piped selection (single-process, no nesting)"
rm -f stacked_structure.fdf symmetry_report.txt
printf '0\n' | stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf -i --no-intro > log_interactive_cli.txt 2>&1
check_contains "ZSL Commensurate Supercells Found" log_interactive_cli.txt
check_success stacked_structure.fdf
rm -f stacked_structure.fdf symmetry_report.txt

echo "Testing: -i/--interactive quitting with 'q'"
printf 'q\n' | stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf -i --no-intro > log_interactive_quit.txt 2>&1
check_exit_code $? 0
if [ ! -s stacked_structure.fdf ]; then
    echo -e "   -> ${GREEN}Verified:${NC} no file created after quitting with 'q'"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} a file was created despite quitting with 'q'"
    FAIL=$((FAIL+1))
fi


# --- 8. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: missing layer file"
stb-2Dstacking -l1 does_not_exist.fdf -l2 hbn.fdf --no-intro > log_missing.txt 2>&1
check_exit_code $? 1
check_contains "Error parsing" log_missing.txt

echo "Testing: malformed fdf (missing ChemicalSpeciesLabel)"
stb-2Dstacking -l1 malformed.fdf -l2 hbn.fdf --no-intro > log_malformed.txt 2>&1
check_exit_code $? 1
check_contains "Error parsing" log_malformed.txt

echo "Testing: -i and -id together (mutually exclusive)"
stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf -i -id 0 --no-intro > log_mutex_match.txt 2>&1
check_exit_code $? 2
check_contains "not allowed with argument" log_mutex_match.txt

echo "Testing: -g and --gap_range together (mutually exclusive)"
stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf -g 3.0 --gap_range 3.0 4.0 3 --no-intro > log_mutex_gap.txt 2>&1
check_exit_code $? 2
check_contains "not allowed with argument" log_mutex_gap.txt

echo "Testing: no lattice match found (max_area/max_strain too tight)"
stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf -a 1 -s 0.0001 --no-intro > log_no_match.txt 2>&1
check_exit_code $? 1
check_contains "No lattice match found" log_no_match.txt

echo "Testing: -l1/-l2 missing (required)"
stb-2Dstacking --no-intro > log_missing_required.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-2Dstacking --version > log_version.txt 2>&1
check_contains "stb-2Dstacking" log_version.txt

echo "Testing: --help documents the main options"
stb-2Dstacking --help > log_help.txt 2>&1
check_contains "batch_sym" log_help.txt
check_contains "gap_range" log_help.txt


# --- 9. Interactive path (stb-suite, shortcut 1.4) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 1.4) ---"
# Note: run_2d_stacker() always passes '-i' to stb-2Dstacking, so it launches
# a SECOND, nested interactive prompt inside the subprocess. Piping enough
# stdin for both layers of prompts does not work reliably here: Python's
# input() in the parent process buffers ahead and can drain the whole piped
# stream before handing the (now EOF) file descriptor to the subprocess --
# a well-known stdio artifact of chaining input() across a subprocess
# boundary over a non-tty pipe, not a bug in the tool itself (confirmed by
# the standalone `-i` test above working fine as a single process). So this
# check only verifies navigation, file validation, and that the wrapper
# builds/launches the command correctly -- not the nested prompt itself.

echo "Testing: navigate 1.4 -> invalid layer1 then valid -> valid layer2 -> defaults -> subprocess launches"
printf '2.1\ndoes_not_exist.fdf\ngraphene.fdf\nhbn.fdf\n\n\n\n\n\n0\n\n0\n' | stb-suite > log_interactive.txt 2>&1
check_contains "File not found" log_interactive.txt
check_contains "Searching for commensurate supercells" log_interactive.txt
rm -f stacked_structure.fdf symmetry_report.txt


popd > /dev/null

# --- 10. Summary ---
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
