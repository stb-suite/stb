#!/bin/bash

# --- Setup ---
# Smoke test for stb-elasticInputs (Elastic Constants Prep, item 4.2.1)
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
    if grep -q "$1" "$2" 2>/dev/null; then
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
echo "--- Starting tester for STB-ElasticInputs prep (item 4.2.1) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/structure.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/si_cubic.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/triclinic.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Physics: vacuum-axis filtering on the 2D hexagonal fixture ---
echo -e "\n--- Testing --dirs all on a 2D fixture (vacuum in c) ---"
stb-elasticInputs --file structure.fdf --dirs all --no-intro > log_2d_all.txt 2>&1
check_exit_code $? 0
check_contains "Detected dimensionality: 2D" log_2d_all.txt
check_contains "Skipping direction 'zz'" log_2d_all.txt
check_contains "Skipping direction 'xz'" log_2d_all.txt
check_contains "Skipping direction 'yz'" log_2d_all.txt
check_contains "Modes: \['xx', 'yy', 'xy'\]" log_2d_all.txt
check_success strain_xx_2.00/structure.fdf
check_success strain_yy_2.00/structure.fdf
check_success strain_xy_2.00/structure.fdf
if [ -d strain_zz_2.00 ] || [ -d strain_xz_2.00 ] || [ -d strain_yz_2.00 ]; then
    echo -e "   -> ${RED}Failed:${NC} a vacuum-touching direction was generated anyway"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no vacuum-touching direction folders were created"
    PASS=$((PASS+1))
fi

echo "Testing: explicit --dirs zz on the same 2D fixture -> nothing to generate"
rm -rf strain_*
stb-elasticInputs --file structure.fdf --dirs zz --no-intro > log_2d_zz.txt 2>&1
check_exit_code $? 1
check_contains "None of the requested directions are physically valid" log_2d_zz.txt
if [ -d strain_zz_2.00 ]; then
    echo -e "   -> ${RED}Failed:${NC} strain_zz_2.00 was created despite being vacuum-only"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} nothing generated"
    PASS=$((PASS+1))
fi


# --- 3. Physics: cubic 3D fixture -- 'all' auto-reduces to 1 direction per
#     symmetry-equivalent group (real DFT-cost cut); an explicit full list
#     is never reduced. No vacuum warnings either way (no vacuum here). ---
echo -e "\n--- Testing --dirs all on a cubic 3D fixture (auto-reduction) ---"
rm -rf strain_* reference_structure.fdf
stb-elasticInputs --file si_cubic.fdf --dirs all --no-intro > log_cubic.txt 2>&1
check_exit_code $? 0
check_contains "Detected dimensionality: 3D" log_cubic.txt
check_contains "DEFORMATION DIRECTIONS (symmetry-method basic)" log_cubic.txt
check_contains "point group m-3m" log_cubic.txt
check_contains "48 operation" log_cubic.txt
check_contains "yy.*SUPPRESSED.*xx" log_cubic.txt
check_contains "zz.*SUPPRESSED.*xx" log_cubic.txt
check_contains "xz.*SUPPRESSED.*xy" log_cubic.txt
check_contains "yz.*SUPPRESSED.*xy" log_cubic.txt
check_contains "Modes: \['xx', 'xy'\]" log_cubic.txt
check_success strain_xx_2.00/structure.fdf
check_success strain_xy_2.00/structure.fdf
check_success reference_structure.fdf
for d in yy zz xz yz; do
    if [ -d "strain_${d}_2.00" ]; then
        echo -e "   -> ${RED}Failed:${NC} redundant direction '$d' was generated anyway"
        FAIL=$((FAIL+1))
    else
        echo -e "   -> ${GREEN}Verified:${NC} redundant direction '$d' correctly skipped"
        PASS=$((PASS+1))
    fi
done
if grep -qi "WARNING" log_cubic.txt 2>/dev/null; then
    echo -e "   -> ${RED}Failed:${NC} unexpected warning for a high-symmetry 3D structure"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no vacuum/symmetry warnings for the cubic fixture"
    PASS=$((PASS+1))
fi

echo "Testing: explicit full --dirs list is never auto-reduced"
rm -rf strain_* reference_structure.fdf
stb-elasticInputs --file si_cubic.fdf --dirs xx yy zz xy xz yz --no-intro > log_cubic_explicit.txt 2>&1
check_exit_code $? 0
if grep -q "Skipping" log_cubic_explicit.txt 2>/dev/null; then
    echo -e "   -> ${RED}Failed:${NC} an explicit --dirs list was reduced (should never be)"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} explicit --dirs list generated in full, no reduction"
    PASS=$((PASS+1))
fi
for d in xx yy zz xy xz yz; do
    check_success "strain_${d}_2.00/structure.fdf"
done


# --- 4. Physics: triclinic fixture triggers the symmetry warning ---
echo -e "\n--- Testing the triclinic/monoclinic symmetry warning ---"
rm -rf strain_*
stb-elasticInputs --file triclinic.fdf --dirs xx --no-intro > log_triclinic.txt 2>&1
check_exit_code $? 0
check_contains "Detected triclinic symmetry" log_triclinic.txt
check_contains "keep all 6 default directions" log_triclinic.txt


# --- 5. Physics: --symmetry-method full -- the general tensor-invariance
#     method, which (unlike 'basic') also catches hexagonal/trigonal
#     reductions since it constrains the FULL elastic tensor's symmetry
#     rather than matching individual strain-tensor pairs. ---
echo -e "\n--- Testing --symmetry-method full (general tensor-invariance reduction) ---"

echo "Testing: cubic fixture -- full method agrees with basic (2 of 6 directions)"
rm -rf strain_* reference_structure.fdf
stb-elasticInputs --file si_cubic.fdf --dirs all --symmetry-method full --no-intro > log_full_cubic.txt 2>&1
check_exit_code $? 0
check_contains "Point group m-3m has 3 independent elastic constant" log_full_cubic.txt
check_contains "Modes: \['xx', 'xy'\]" log_full_cubic.txt
check_success strain_xx_2.00/structure.fdf
check_success strain_xy_2.00/structure.fdf
for d in yy zz xz yz; do
    if [ -d "strain_${d}_2.00" ]; then
        echo -e "   -> ${RED}Failed:${NC} redundant direction '$d' was generated anyway (full method)"
        FAIL=$((FAIL+1))
    else
        echo -e "   -> ${GREEN}Verified:${NC} redundant direction '$d' correctly skipped (full method)"
        PASS=$((PASS+1))
    fi
done

echo "Testing: 2D hexagonal fixture -- full method reduces to 1 direction (basic needs all 6)"
rm -rf strain_* reference_structure.fdf
stb-elasticInputs --file structure.fdf --dirs all --symmetry-method full --no-intro > log_full_hex.txt 2>&1
check_exit_code $? 0
check_contains "Point group -6m2 has 5 independent elastic constant" log_full_hex.txt
check_contains "Modes: \['xx'\]" log_full_hex.txt
check_success strain_xx_2.00/structure.fdf
for d in yy xy; do
    if [ -d "strain_${d}_2.00" ]; then
        echo -e "   -> ${RED}Failed:${NC} redundant direction '$d' was generated anyway (full method)"
        FAIL=$((FAIL+1))
    else
        echo -e "   -> ${GREEN}Verified:${NC} redundant direction '$d' correctly skipped (full method)"
        PASS=$((PASS+1))
    fi
done

echo "Testing: triclinic fixture -- full method needs all 6 directions (no reduction possible)"
rm -rf strain_* reference_structure.fdf
stb-elasticInputs --file triclinic.fdf --dirs all --symmetry-method full --no-intro > log_full_triclinic.txt 2>&1
check_exit_code $? 0
check_contains "Modes: \['xx', 'yy', 'zz', 'xy', 'xz', 'yz'\]" log_full_triclinic.txt
for d in xx yy zz xy xz yz; do
    check_success "strain_${d}_2.00/structure.fdf"
done

echo "Testing: explicit --dirs list is never reduced, even with --symmetry-method full"
rm -rf strain_* reference_structure.fdf
stb-elasticInputs --file si_cubic.fdf --dirs xx yy zz xy xz yz --symmetry-method full --no-intro > log_full_explicit.txt 2>&1
check_exit_code $? 0
if grep -qi "Skipping" log_full_explicit.txt 2>/dev/null; then
    echo -e "   -> ${RED}Failed:${NC} an explicit --dirs list was reduced under --symmetry-method full"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} explicit --dirs list generated in full, no reduction"
    PASS=$((PASS+1))
fi
for d in xx yy zz xy xz yz; do
    check_success "strain_${d}_2.00/structure.fdf"
done


# --- 6. Physics: --method energy -- energy-strain (parabolic fit) patterns.
#     Unlike --method stress, a pure single-component strain's energy
#     curvature only ever determines its OWN diagonal constant (C_mm);
#     off-diagonal constants need a genuinely COMBINED pattern (two Voigt
#     components strained at once, e.g. 'xx+yy') -- so this method always
#     needs at least as many patterns as independent constants (never fewer
#     DFT calculations than --method stress for the same crystal system). ---
echo -e "\n--- Testing --method energy (energy-strain pattern selection) ---"

echo "Testing: --method energy rejects an explicit --dirs list"
rm -rf strain_* reference_structure.fdf
stb-elasticInputs --file si_cubic.fdf --dirs xx --method energy --no-intro > log_energy_baddirs.txt 2>&1
check_exit_code $? 1
check_contains "always auto-selects its own strain patterns" log_energy_baddirs.txt

echo "Testing: cubic fixture -- needs 3 patterns (2 pure + 1 combined) for 3 independent constants"
rm -rf strain_* reference_structure.fdf
stb-elasticInputs --file si_cubic.fdf --method energy --no-intro > log_energy_cubic.txt 2>&1
check_exit_code $? 0
check_contains "Point group m-3m has 3 independent elastic constant" log_energy_cubic.txt
check_contains "Modes: \['xx', 'yz', 'xx+yy'\]" log_energy_cubic.txt
check_contains "DEFORMATION DIRECTIONS (symmetry-method energy)" log_energy_cubic.txt
check_contains "yy.*SUPPRESSED.*symmetry-allowed fit" log_energy_cubic.txt
check_contains "3 of 21 run -- 18 suppressed" log_energy_cubic.txt
check_success strain_xx_2.00/structure.fdf
check_success strain_yz_2.00/structure.fdf
check_success "strain_xx+yy_2.00/structure.fdf"

echo "Testing: 2D hexagonal fixture -- only in-plane patterns are candidates, needs 2 (xx, xy)"
rm -rf strain_* reference_structure.fdf
stb-elasticInputs --file structure.fdf --method energy --no-intro > log_energy_hex.txt 2>&1
check_exit_code $? 0
check_contains "Point group -6m2 has 5 independent elastic constant" log_energy_hex.txt
check_contains "Modes: \['xx', 'xy'\]" log_energy_hex.txt
check_contains "DEFORMATION DIRECTIONS (symmetry-method energy)" log_energy_hex.txt
check_contains "point group -6m2" log_energy_hex.txt
check_contains "2 of 6 run -- 4 suppressed" log_energy_hex.txt
check_success strain_xx_2.00/structure.fdf
check_success strain_xy_2.00/structure.fdf
if [ -d "strain_zz_2.00" ] || [ -d "strain_xx+zz_2.00" ]; then
    echo -e "   -> ${RED}Failed:${NC} an out-of-plane (vacuum-touching) pattern was generated anyway"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no out-of-plane pattern was generated"
    PASS=$((PASS+1))
fi

echo "Testing: triclinic fixture -- needs all 21 patterns (no reduction possible)"
rm -rf strain_* reference_structure.fdf
stb-elasticInputs --file triclinic.fdf --method energy --no-intro > log_energy_triclinic.txt 2>&1
check_exit_code $? 0
check_contains "Point group -1 has 21 independent elastic constant" log_energy_triclinic.txt
check_contains "All 21 direction(s) run -- no symmetry reduction" log_energy_triclinic.txt
if [ "$(ls -d strain_*/ 2>/dev/null | wc -l)" -eq 84 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} 21 patterns x 4 strain steps = 84 folders generated"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} expected 84 strain folders (21 patterns x 4 steps)"
    FAIL=$((FAIL+1))
fi


# --- 7. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: missing required args"
stb-elasticInputs --dirs xx --no-intro > log_missing_args.txt 2>&1
check_exit_code $? 2

echo "Testing: nonexistent input file"
stb-elasticInputs --file does_not_exist.fdf --dirs xx --no-intro > log_missing_file.txt 2>&1
check_exit_code $? 1

echo "Testing: --version"
stb-elasticInputs --version > log_version.txt 2>&1
check_contains "stb-elasticInputs" log_version.txt

echo "Testing: --help documents --dirs/--vacuum-gap/--symprec/--angle-tolerance"
stb-elasticInputs --help > log_help.txt 2>&1
check_contains "dirs" log_help.txt
check_contains "vacuum-gap" log_help.txt
check_contains "symprec" log_help.txt
check_contains "angle-tolerance" log_help.txt
check_contains "symmetry-method" log_help.txt
check_contains "method" log_help.txt


# --- 8. Interactive path (stb-suite, shortcut 4.2.1) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.2.1) ---"

echo "Testing: navigate 4.2.1 -> invalid file then valid -> defaults -> stress method -> basic symmetry reduction -> skip advanced -> quit"
rm -rf strain_*
# Prompts in order: file (retry), max strain, steps, method (blank -> stress),
# mode (1 -> basic, the only "which direction subset" question left -- all 6
# canonical directions are always requested, auto-filtered by dimensionality,
# and symmetry-reduced), advanced settings (n -> skip), then the "Press Enter
# to continue" pause, then quit.
printf '4.2.1\ndoes_not_exist.fdf\nstructure.fdf\n\n\n\n1\nn\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_contains "File not found" log_menu.txt
check_contains "Modes: \['xx', 'yy', 'xy'\]" log_menu.txt
check_contains "CONFIGURATION SUMMARY" log_menu.txt
check_contains "DEFORMATION DIRECTIONS (symmetry-method basic)" log_menu.txt
check_success strain_xx_2.00/structure.fdf
if grep -q "Traceback" log_menu.txt 2>/dev/null; then
    echo -e "   -> ${RED}Failed:${NC} unexpected traceback in interactive session"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} interactive session exited cleanly"
    PASS=$((PASS+1))
fi

echo "Testing: navigate 4.2.1 -> energy method never asks a symmetry-method question (elastic_inputs.py never reads it for --method energy)"
rm -rf strain_*
# Prompts in order: file, max strain, steps, method (2 -> energy, no
# direction/symmetry-method menu at all), advanced settings (n -> skip),
# then the "Press Enter to continue" pause, then quit.
printf '4.2.1\nstructure.fdf\n\n\n2\nn\n\n0\n' | stb-suite > log_energy_menu.txt 2>&1
check_contains "auto-selects its own pure+combined strain patterns" log_energy_menu.txt
check_contains "CONFIGURATION SUMMARY" log_energy_menu.txt
if grep -q "Symmetry reduction" log_energy_menu.txt 2>/dev/null; then
    echo -e "   -> ${RED}Failed:${NC} symmetry-method question appeared for --method energy (it's never read there)"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no symmetry-method question for --method energy"
    PASS=$((PASS+1))
fi
if grep -q "Traceback" log_energy_menu.txt 2>/dev/null; then
    echo -e "   -> ${RED}Failed:${NC} unexpected traceback in interactive session"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} interactive session exited cleanly"
    PASS=$((PASS+1))
fi


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
