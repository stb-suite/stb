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
cp "$FIXTURE_DIR/calc.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/calc_cubic.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/calc_triclinic.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/calc_relaxation.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/cartesian_coords.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Physics: vacuum-axis filtering on the 2D hexagonal fixture ---
echo -e "\n--- Testing --dirs all on a 2D fixture (vacuum in c) ---"
rm -rf elastic_runs
stb-elasticInputs -s structure.fdf -c calc.fdf --dirs all --no-intro > log_2d_all.txt 2>&1
check_exit_code $? 0
check_contains "Detected dimensionality: 2D" log_2d_all.txt
check_contains "Skipping direction 'zz'" log_2d_all.txt
check_contains "Skipping direction 'xz'" log_2d_all.txt
check_contains "Skipping direction 'yz'" log_2d_all.txt
check_contains "Modes: \['xx', 'yy', 'xy'\]" log_2d_all.txt
check_success elastic_runs/xx/strain_xx_2.00/structure.fdf
check_success elastic_runs/yy/strain_yy_2.00/structure.fdf
check_success elastic_runs/xy/strain_xy_2.00/structure.fdf
if [ -d elastic_runs/zz ] || [ -d elastic_runs/xz ] || [ -d elastic_runs/yz ]; then
    echo -e "   -> ${RED}Failed:${NC} a vacuum-touching direction was generated anyway"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no vacuum-touching direction folders were created"
    PASS=$((PASS+1))
fi

echo "Testing: explicit --dirs zz on the same 2D fixture -> nothing to generate"
rm -rf elastic_runs
stb-elasticInputs -s structure.fdf -c calc.fdf --dirs zz --no-intro > log_2d_zz.txt 2>&1
check_exit_code $? 1
check_contains "None of the requested directions are physically valid" log_2d_zz.txt
if [ -d elastic_runs/zz ]; then
    echo -e "   -> ${RED}Failed:${NC} elastic_runs/zz was created despite being vacuum-only"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} nothing generated"
    PASS=$((PASS+1))
fi


# --- 3. Physics: cubic 3D fixture -- 'all' auto-reduces to 1 direction per
#     symmetry-equivalent group (real DFT-cost cut); an explicit full list
#     is never reduced. No vacuum warnings either way (no vacuum here). ---
echo -e "\n--- Testing --dirs all on a cubic 3D fixture (auto-reduction) ---"
rm -rf elastic_runs
stb-elasticInputs -s si_cubic.fdf -c calc_cubic.fdf --dirs all --no-intro > log_cubic.txt 2>&1
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
check_success elastic_runs/xx/strain_xx_2.00/si_cubic.fdf
check_success elastic_runs/xy/strain_xy_2.00/si_cubic.fdf
check_success elastic_runs/reference_structure.fdf
for d in yy zz xz yz; do
    if [ -d "elastic_runs/${d}" ]; then
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
rm -rf elastic_runs
stb-elasticInputs -s si_cubic.fdf -c calc_cubic.fdf --dirs xx yy zz xy xz yz --no-intro > log_cubic_explicit.txt 2>&1
check_exit_code $? 0
if grep -q "Skipping" log_cubic_explicit.txt 2>/dev/null; then
    echo -e "   -> ${RED}Failed:${NC} an explicit --dirs list was reduced (should never be)"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} explicit --dirs list generated in full, no reduction"
    PASS=$((PASS+1))
fi
for d in xx yy zz xy xz yz; do
    check_success "elastic_runs/${d}/strain_${d}_2.00/si_cubic.fdf"
done


# --- 4. Physics: triclinic fixture triggers the symmetry warning ---
echo -e "\n--- Testing the triclinic/monoclinic symmetry warning ---"
rm -rf elastic_runs
stb-elasticInputs -s triclinic.fdf -c calc_triclinic.fdf --dirs xx --no-intro > log_triclinic.txt 2>&1
check_exit_code $? 0
check_contains "Detected triclinic symmetry" log_triclinic.txt
check_contains "keep all 6 default directions" log_triclinic.txt


# --- 5. Physics: --symmetry-method full -- the general tensor-invariance
#     method, which (unlike 'basic') also catches hexagonal/trigonal
#     reductions since it constrains the FULL elastic tensor's symmetry
#     rather than matching individual strain-tensor pairs. ---
echo -e "\n--- Testing --symmetry-method full (general tensor-invariance reduction) ---"

echo "Testing: cubic fixture -- full method agrees with basic (2 of 6 directions)"
rm -rf elastic_runs
stb-elasticInputs -s si_cubic.fdf -c calc_cubic.fdf --dirs all --symmetry-method full --no-intro > log_full_cubic.txt 2>&1
check_exit_code $? 0
check_contains "Point group m-3m has 3 independent elastic constant" log_full_cubic.txt
check_contains "Modes: \['xx', 'xy'\]" log_full_cubic.txt
check_success elastic_runs/xx/strain_xx_2.00/si_cubic.fdf
check_success elastic_runs/xy/strain_xy_2.00/si_cubic.fdf
for d in yy zz xz yz; do
    if [ -d "elastic_runs/${d}" ]; then
        echo -e "   -> ${RED}Failed:${NC} redundant direction '$d' was generated anyway (full method)"
        FAIL=$((FAIL+1))
    else
        echo -e "   -> ${GREEN}Verified:${NC} redundant direction '$d' correctly skipped (full method)"
        PASS=$((PASS+1))
    fi
done

echo "Testing: 2D hexagonal fixture -- full method reduces to 1 direction (basic needs all 6)"
rm -rf elastic_runs
stb-elasticInputs -s structure.fdf -c calc.fdf --dirs all --symmetry-method full --no-intro > log_full_hex.txt 2>&1
check_exit_code $? 0
check_contains "Point group -6m2 has 5 independent elastic constant" log_full_hex.txt
check_contains "Modes: \['xx'\]" log_full_hex.txt
check_success elastic_runs/xx/strain_xx_2.00/structure.fdf
for d in yy xy; do
    if [ -d "elastic_runs/${d}" ]; then
        echo -e "   -> ${RED}Failed:${NC} redundant direction '$d' was generated anyway (full method)"
        FAIL=$((FAIL+1))
    else
        echo -e "   -> ${GREEN}Verified:${NC} redundant direction '$d' correctly skipped (full method)"
        PASS=$((PASS+1))
    fi
done

echo "Testing: triclinic fixture -- full method needs all 6 directions (no reduction possible)"
rm -rf elastic_runs
stb-elasticInputs -s triclinic.fdf -c calc_triclinic.fdf --dirs all --symmetry-method full --no-intro > log_full_triclinic.txt 2>&1
check_exit_code $? 0
check_contains "Modes: \['xx', 'yy', 'zz', 'xy', 'xz', 'yz'\]" log_full_triclinic.txt
for d in xx yy zz xy xz yz; do
    check_success "elastic_runs/${d}/strain_${d}_2.00/triclinic.fdf"
done

echo "Testing: explicit --dirs list is never reduced, even with --symmetry-method full"
rm -rf elastic_runs
stb-elasticInputs -s si_cubic.fdf -c calc_cubic.fdf --dirs xx yy zz xy xz yz --symmetry-method full --no-intro > log_full_explicit.txt 2>&1
check_exit_code $? 0
if grep -qi "Skipping" log_full_explicit.txt 2>/dev/null; then
    echo -e "   -> ${RED}Failed:${NC} an explicit --dirs list was reduced under --symmetry-method full"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} explicit --dirs list generated in full, no reduction"
    PASS=$((PASS+1))
fi
for d in xx yy zz xy xz yz; do
    check_success "elastic_runs/${d}/strain_${d}_2.00/si_cubic.fdf"
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
rm -rf elastic_runs
stb-elasticInputs -s si_cubic.fdf -c calc_cubic.fdf --dirs xx --method energy --no-intro > log_energy_baddirs.txt 2>&1
check_exit_code $? 1
check_contains "always auto-selects its own strain patterns" log_energy_baddirs.txt

echo "Testing: cubic fixture -- needs 3 patterns (2 pure + 1 combined) for 3 independent constants"
rm -rf elastic_runs
stb-elasticInputs -s si_cubic.fdf -c calc_cubic.fdf --method energy --no-intro > log_energy_cubic.txt 2>&1
check_exit_code $? 0
check_contains "Point group m-3m has 3 independent elastic constant" log_energy_cubic.txt
check_contains "Modes: \['xx', 'yz', 'xx+yy'\]" log_energy_cubic.txt
check_contains "DEFORMATION DIRECTIONS (symmetry-method energy)" log_energy_cubic.txt
check_contains "yy.*SUPPRESSED.*symmetry-allowed fit" log_energy_cubic.txt
check_contains "3 of 21 run -- 18 suppressed" log_energy_cubic.txt
check_success elastic_runs/xx/strain_xx_2.00/si_cubic.fdf
check_success elastic_runs/yz/strain_yz_2.00/si_cubic.fdf
check_success "elastic_runs/xx+yy/strain_xx+yy_2.00/si_cubic.fdf"

echo "Testing: 2D hexagonal fixture -- only in-plane patterns are candidates, needs 2 (xx, xy)"
rm -rf elastic_runs
stb-elasticInputs -s structure.fdf -c calc.fdf --method energy --no-intro > log_energy_hex.txt 2>&1
check_exit_code $? 0
check_contains "Point group -6m2 has 5 independent elastic constant" log_energy_hex.txt
check_contains "Modes: \['xx', 'xy'\]" log_energy_hex.txt
check_contains "DEFORMATION DIRECTIONS (symmetry-method energy)" log_energy_hex.txt
check_contains "point group -6m2" log_energy_hex.txt
check_contains "2 of 6 run -- 4 suppressed" log_energy_hex.txt
check_success elastic_runs/xx/strain_xx_2.00/structure.fdf
check_success elastic_runs/xy/strain_xy_2.00/structure.fdf
if [ -d "elastic_runs/zz" ] || [ -d "elastic_runs/xx+zz" ]; then
    echo -e "   -> ${RED}Failed:${NC} an out-of-plane (vacuum-touching) pattern was generated anyway"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no out-of-plane pattern was generated"
    PASS=$((PASS+1))
fi

echo "Testing: triclinic fixture -- needs all 21 patterns (no reduction possible)"
rm -rf elastic_runs
stb-elasticInputs -s triclinic.fdf -c calc_triclinic.fdf --method energy --no-intro > log_energy_triclinic.txt 2>&1
check_exit_code $? 0
check_contains "Point group -1 has 21 independent elastic constant" log_energy_triclinic.txt
check_contains "All 21 direction(s) run -- no symmetry reduction" log_energy_triclinic.txt
if [ "$(find elastic_runs -mindepth 2 -maxdepth 2 -type d -name 'strain_*' 2>/dev/null | wc -l)" -eq 84 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} 21 patterns x 4 strain steps = 84 folders generated"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} expected 84 strain folders (21 patterns x 4 steps)"
    FAIL=$((FAIL+1))
fi


# --- 7. Physics: Cartesian-coordinate input is rejected ---
echo -e "\n--- Testing rejection of Cartesian-coordinate input ---"
rm -rf elastic_runs
stb-elasticInputs -s cartesian_coords.fdf -c calc_cubic.fdf --dirs xx --no-intro > log_cartesian.txt 2>&1
check_exit_code $? 1
check_contains "Atomic coordinates in this file are Cartesian" log_cartesian.txt


# --- 8. New feature: calc.fdf copy + forced single-point SCF, verifying the
#     actual precedence guarantee (config_extra.fdf's override must be
#     PREPENDED, and physically read before the base template's own MD
#     directives, or SIESTA's first-occurrence-wins fdf parser would
#     silently keep the base template's real relaxation settings instead). ---
echo -e "\n--- Testing calc.fdf copy + config_extra.fdf + forced single-point precedence ---"
rm -rf elastic_runs
stb-elasticInputs -s si_cubic.fdf -c calc_relaxation.fdf --dirs xx --max 1.0 --steps 3 --no-intro > log_singlepoint.txt 2>&1
check_exit_code $? 0
check_contains "SINGLE-POINT SCF ENFORCEMENT" log_singlepoint.txt
check_contains "MD.TypeOfRun=CG  Steps: MD.Steps=100  MD.VariableCell=true" log_singlepoint.txt
check_success elastic_runs/xx/strain_xx_1.00/si_cubic.fdf
check_success elastic_runs/xx/strain_xx_1.00/calc_relaxation.fdf
check_success elastic_runs/xx/strain_xx_1.00/config_extra.fdf
check_contains "MD.Steps           0" elastic_runs/xx/strain_xx_1.00/config_extra.fdf
check_contains "MD.VariableCell    false" elastic_runs/xx/strain_xx_1.00/config_extra.fdf

echo "Testing: '%include config_extra.fdf' is the very FIRST line (prepended, not appended)"
first_line=$(head -1 elastic_runs/xx/strain_xx_1.00/calc_relaxation.fdf)
if [ "$first_line" = "%include config_extra.fdf" ]; then
    echo -e "   -> ${GREEN}Verified:${NC} config_extra.fdf is %included as the very first line"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} expected '%include config_extra.fdf' as the first line, got '$first_line'"
    FAIL=$((FAIL+1))
fi

echo "Testing: the base template's own (real, active) MD.Steps=100/MD.VariableCell=true"
echo "         lines are preserved verbatim further down (only the %include is added)"
check_contains "MD.Steps                100" elastic_runs/xx/strain_xx_1.00/calc_relaxation.fdf
check_contains "MD.VariableCell         true" elastic_runs/xx/strain_xx_1.00/calc_relaxation.fdf


# --- 9. New feature: pseudopotential linking ---
echo -e "\n--- Testing pseudopotential linking (-p dojo) ---"
rm -rf elastic_runs
stb-elasticInputs -s si_cubic.fdf -c calc_cubic.fdf -p dojo --dirs xx --max 1.0 --steps 3 --no-intro > log_pseudo.txt 2>&1
check_exit_code $? 0
check_contains "Using bundled pseudopotential bank 'dojo'" log_pseudo.txt
check_contains "\[OK\] All required pseudopotentials found" log_pseudo.txt
check_success elastic_runs/xx/strain_xx_1.00/Si.psml
if [ -L elastic_runs/xx/strain_xx_1.00/Si.psml ]; then
    echo -e "   -> ${GREEN}Verified:${NC} Si.psml is a symlink (not a copy)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} Si.psml is not a symlink"
    FAIL=$((FAIL+1))
fi

echo "Testing: missing-pseudopotential warning (custom empty folder)"
mkdir -p empty_bank
rm -rf elastic_runs
stb-elasticInputs -s si_cubic.fdf -c calc_cubic.fdf -p empty_bank --dirs xx --no-intro > log_missing_pp.txt 2>&1
check_exit_code $? 0
check_contains "Si.*MISSING" log_missing_pp.txt
check_contains "\[WARNING\] Missing pseudopotential(s) for: Si" log_missing_pp.txt
if ls elastic_runs/xx/strain_xx_2.00/*.psml elastic_runs/xx/strain_xx_2.00/*.psf > /dev/null 2>&1; then
    echo -e "   -> ${RED}Failed:${NC} a pseudopotential file was linked despite none being available"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no pseudopotential file linked when unavailable"
    PASS=$((PASS+1))
fi

echo "Testing: nonexistent pseudo bank/path is rejected cleanly"
stb-elasticInputs -s si_cubic.fdf -c calc_cubic.fdf -p does_not_exist_bank --dirs xx --no-intro > log_bad_pseudo.txt 2>&1
check_exit_code $? 1
check_contains "not a recognized pseudopotential bank" log_bad_pseudo.txt


# --- 10. New feature: --output-dir + --save-report ---
echo -e "\n--- Testing --save-report + a custom --output-dir ---"
rm -rf my_elastic_out
stb-elasticInputs -s si_cubic.fdf -c calc_cubic.fdf --dirs xx --output-dir my_elastic_out --save-report --no-intro > log_savereport.txt 2>&1
check_exit_code $? 0
check_success my_elastic_out/elastic_stage1.txt
check_success my_elastic_out/reference_structure.fdf
check_success my_elastic_out/xx/strain_xx_2.00/si_cubic.fdf
rm -rf my_elastic_out


# --- 11. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: missing required args"
stb-elasticInputs --dirs xx --no-intro > log_missing_args.txt 2>&1
check_exit_code $? 2

echo "Testing: nonexistent input file"
stb-elasticInputs -s does_not_exist.fdf -c calc.fdf --dirs xx --no-intro > log_missing_file.txt 2>&1
check_exit_code $? 1
check_contains "Structure file 'does_not_exist.fdf' not found" log_missing_file.txt

echo "Testing: nonexistent calc.fdf"
stb-elasticInputs -s si_cubic.fdf -c does_not_exist_calc.fdf --dirs xx --no-intro > log_missing_calc.txt 2>&1
check_exit_code $? 1
check_contains "Calc file 'does_not_exist_calc.fdf' not found" log_missing_calc.txt

echo "Testing: --version"
stb-elasticInputs --version > log_version.txt 2>&1
check_contains "stb-elasticInputs" log_version.txt

echo "Testing: --help documents --structure/--calc/--pseudo-dir/--dirs/--vacuum-gap/--symprec/--angle-tolerance/--output-dir/--save-report"
stb-elasticInputs --help > log_help.txt 2>&1
check_contains "\-\-structure" log_help.txt
check_contains "\-\-calc" log_help.txt
check_contains "pseudo-dir" log_help.txt
check_contains "dirs" log_help.txt
check_contains "vacuum-gap" log_help.txt
check_contains "symprec" log_help.txt
check_contains "angle-tolerance" log_help.txt
check_contains "symmetry-method" log_help.txt
check_contains "method" log_help.txt
check_contains "output-dir" log_help.txt
check_contains "save-report" log_help.txt


# --- 12. Interactive path (stb-suite, shortcut 4.2.1) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.2.1) ---"

echo "Testing: navigate 4.2.1 -> invalid file then valid -> default calc.fdf -> skip pseudo ->"
echo "         defaults -> stress method -> basic symmetry reduction -> skip advanced -> quit"
rm -rf elastic_runs
# Prompts in order: structure file (retry), calc.fdf (blank -> default 'calc.fdf'),
# pseudo source (blank -> skip), max strain (blank), steps (blank), method (blank
# -> stress), mode (1 -> basic), advanced settings (n -> skip), save report (blank
# -> no), "Press Enter to continue", quit.
printf '4.2.1\ndoes_not_exist.fdf\nstructure.fdf\n\n\n\n\n\n1\nn\n\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_contains "File not found" log_menu.txt
check_contains "Modes: \['xx', 'yy', 'xy'\]" log_menu.txt
check_contains "CONFIGURATION SUMMARY" log_menu.txt
check_contains "DEFORMATION DIRECTIONS (symmetry-method basic)" log_menu.txt
check_success elastic_runs/xx/strain_xx_2.00/structure.fdf
if grep -q "Traceback" log_menu.txt 2>/dev/null; then
    echo -e "   -> ${RED}Failed:${NC} unexpected traceback in interactive session"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} interactive session exited cleanly"
    PASS=$((PASS+1))
fi

echo "Testing: navigate 4.2.1 -> energy method never asks a symmetry-method question"
rm -rf elastic_runs
# Prompts in order: structure file, calc.fdf (blank -> default), pseudo (blank ->
# skip), max strain (blank), steps (blank), method (2 -> energy, no direction/
# symmetry-method menu at all), advanced settings (n -> skip), save report (blank
# -> no), "Press Enter to continue", quit.
printf '4.2.1\nstructure.fdf\n\n\n\n\n2\nn\n\n\n0\n' | stb-suite > log_energy_menu.txt 2>&1
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

# --- 13. Summary ---
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
