#!/bin/bash

# --- Setup ---
# Smoke test for stb-strain (Stress-Strain Prep, item 4.1.1)
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

# Checks that file $1 does NOT exist (or is empty)
check_absent() {
    if [ ! -s "$1" ]; then
        echo -e " ... ${GREEN}OK${NC} (file '$1' correctly not created)"
        PASS=$((PASS+1))
    else
        echo -e " ... ${RED}FAIL${NC} (file '$1' was unexpectedly created)"
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

# Checks that file $2 does NOT contain (grep -q) pattern $1
check_not_contains() {
    if grep -q "$1" "$2" 2>/dev/null; then
        echo -e "   -> ${RED}Failed:${NC} '$1' unexpectedly found in '$2'"
        FAIL=$((FAIL+1))
    else
        echo -e "   -> ${GREEN}Verified:${NC} '$1' correctly absent from '$2'"
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
echo "--- Starting tester for STB-Strain prep (item 4.1.1) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR/subdir"
cp "$FIXTURE_DIR/structure.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/structure.fdf" "$TEST_DIR/subdir/"
cp "$FIXTURE_DIR/calc.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/calc.fdf" "$TEST_DIR/subdir/"
cp "$FIXTURE_DIR/si_cubic.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/calc_cubic.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/calc_relaxed_cubic.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/calc_zero_steps.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/cartesian_coords.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Uniaxial strain (x), stress-constrained mode, written inside the strain_runs/ wrapper dir ---
echo -e "\n--- Testing --stdir x (uniaxial), --relax-mode stress-constrained, numbered report ---"
stb-strain -s structure.fdf -c calc.fdf --relax-mode stress-constrained --stdir x --stmin 0 --stmax 2 --step 2 --no-intro > log_uniaxial.txt 2>&1
check_contains "\[0\] RUN METADATA" log_uniaxial.txt
check_contains "\[1\] INPUT STRUCTURE" log_uniaxial.txt
check_contains "\[2\] DIMENSIONALITY & DIRECTION VALIDATION" log_uniaxial.txt
check_contains "\[3\] AXIS SYMMETRY ADVISORY" log_uniaxial.txt
check_contains "\[4\] RELAXATION MODE & CELL CONSTRAINTS" log_uniaxial.txt
check_contains "\[5\] PSEUDOPOTENTIALS" log_uniaxial.txt
check_contains "\[6\] GENERATED STRAIN FOLDERS" log_uniaxial.txt
check_contains "\[7\] SUMMARY & NEXT STEPS" log_uniaxial.txt
check_contains "Structure file    : structure.fdf" log_uniaxial.txt
check_contains "Calc template     : calc.fdf" log_uniaxial.txt
check_contains "Relax mode        : stress-constrained" log_uniaxial.txt
check_contains "Direction         : x (uniaxial)" log_uniaxial.txt
check_contains "Detected dimensionality: 2D" log_uniaxial.txt
check_success strain_runs/x/strain_x_0.00/structure.fdf
check_success strain_runs/x/strain_x_0.00/calc.fdf
check_success strain_runs/x/strain_x_2.00/structure.fdf
check_contains "6.026159992248" strain_runs/x/strain_x_2.00/structure.fdf
check_contains "5.116478079000" strain_runs/x/strain_x_2.00/structure.fdf

echo "Testing: calc.fdf and config_extra.fdf are also copied into every generated folder, identical across folders"
check_contains "%include structure.fdf" strain_runs/x/strain_x_2.00/calc.fdf
check_contains "%include config_extra.fdf" strain_runs/x/strain_x_2.00/calc.fdf
check_success strain_runs/x/strain_x_2.00/config_extra.fdf
if diff -q strain_runs/x/strain_x_0.00/calc.fdf strain_runs/x/strain_x_2.00/calc.fdf > /dev/null 2>&1; then
    echo -e "   -> ${GREEN}Verified:${NC} calc.fdf is byte-identical across strain steps"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} calc.fdf differs between strain steps (should be identical)"
    FAIL=$((FAIL+1))
fi
if diff -q strain_runs/x/strain_x_0.00/config_extra.fdf strain_runs/x/strain_x_2.00/config_extra.fdf > /dev/null 2>&1; then
    echo -e "   -> ${GREEN}Verified:${NC} config_extra.fdf is byte-identical across strain steps"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} config_extra.fdf differs between strain steps (should be identical)"
    FAIL=$((FAIL+1))
fi

echo "Testing: stress-constrained mode on the 2D fixture (vacuum along z) -- x is imposed, z is vacuum"
echo "         (protected), y is periodic and free, XY shear is free, XZ/YZ shear fixed (touch vacuum)"
check_contains "stress 1  # Fixes XX" strain_runs/x/strain_x_2.00/config_extra.fdf
check_not_contains "stress 2  # Fixes YY" strain_runs/x/strain_x_2.00/config_extra.fdf
check_contains "stress 3  # Fixes ZZ" strain_runs/x/strain_x_2.00/config_extra.fdf
check_contains "stress 4  # Fixes YZ" strain_runs/x/strain_x_2.00/config_extra.fdf
check_contains "stress 5  # Fixes XZ" strain_runs/x/strain_x_2.00/config_extra.fdf
check_not_contains "stress 6  # Fixes XY" strain_runs/x/strain_x_2.00/config_extra.fdf

echo "Testing: no -p/--pseudo-dir given -- no pseudopotential files linked, note printed"
check_contains "Not given (pass -p/--pseudo-dir" log_uniaxial.txt
if ls strain_runs/x/strain_x_2.00/*.psml strain_runs/x/strain_x_2.00/*.psf > /dev/null 2>&1; then
    echo -e "   -> ${RED}Failed:${NC} a pseudopotential file was linked despite no -p given"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no pseudopotential file linked without -p"
    PASS=$((PASS+1))
fi

echo "Testing: no --save-report by default -- no report file written"
check_absent strain_runs/x/strain_stage1.txt


# --- 3. Biaxial strain (xy) ---
echo -e "\n--- Testing --stdir yx (biaxial, normalized to xy) ---"
stb-strain -s structure.fdf -c calc.fdf --relax-mode cell-fixed --stdir yx --stmin 5 --stmax 5 --step 1 --no-intro > log_biaxial.txt 2>&1
check_contains "Direction         : xy (biaxial)" log_biaxial.txt
check_success strain_runs/xy/strain_xy_5.00/structure.fdf
check_contains "6.203399992020" strain_runs/xy/strain_xy_5.00/structure.fdf
check_contains "5.372301982950" strain_runs/xy/strain_xy_5.00/structure.fdf

echo "Testing: cell-fixed mode fixes all 6 Voigt components regardless of direction"
for i in 1 2 3 4 5 6; do
    check_contains "stress $i  # Fixes" strain_runs/xy/strain_xy_5.00/config_extra.fdf
done


# --- 4. Physics: refuse straining a vacuum-padded axis ---
echo -e "\n--- Testing vacuum-axis rejection (--stdir z, the vacuum axis of this 2D fixture) ---"
stb-strain -s structure.fdf -c calc.fdf --relax-mode stress-constrained --stdir z --stmin 0 --stmax 1 --no-intro > log_vacuum_z.txt 2>&1
check_exit_code $? 1
check_contains "vacuum-padded axis" log_vacuum_z.txt
check_contains "Periodic axis/axes available for this structure: x, y" log_vacuum_z.txt

echo "Testing: biaxial direction that includes the vacuum axis (xz)"
stb-strain -s structure.fdf -c calc.fdf --relax-mode stress-constrained --stdir xz --stmin 0 --stmax 1 --no-intro > log_vacuum_xz.txt 2>&1
check_exit_code $? 1
check_contains "vacuum-padded axis" log_vacuum_xz.txt


# --- 5. Symmetry advisory (uniaxial only, informational, never blocks) ---
echo -e "\n--- Testing the symmetry-equivalence advisory ---"

echo "Testing: cubic structure (m-3m) -- x/y/z are truly equivalent, advisory must fire"
rm -rf strain_runs
stb-strain -s si_cubic.fdf -c calc_cubic.fdf --relax-mode cell-fixed --stdir x --stmin 0 --stmax 1 --no-intro > log_sym_cubic.txt 2>&1
check_exit_code $? 0
check_contains "AXIS SYMMETRY (uniaxial direction 'x')" log_sym_cubic.txt
check_contains "point group m-3m" log_sym_cubic.txt
check_contains "48 operation" log_sym_cubic.txt
check_contains "y.*EQUIVALENT.*x" log_sym_cubic.txt
check_contains "z.*EQUIVALENT.*x" log_sym_cubic.txt
check_contains "y, z are equivalent to 'x' by symmetry" log_sym_cubic.txt

echo "Testing: hexagonal 2D fixture -- x and y are NOT equivalent (known zigzag/armchair"
echo "         asymmetry), advisory must NOT fire"
rm -rf strain_runs
stb-strain -s structure.fdf -c calc.fdf --relax-mode stress-constrained --stdir x --stmin 0 --stmax 1 --no-intro > log_sym_hex.txt 2>&1
check_exit_code $? 0
if grep -q "equivalent to 'x' by symmetry" log_sym_hex.txt 2>/dev/null; then
    echo -e "   -> ${RED}Failed:${NC} unexpected symmetry advisory for the hexagonal 2D fixture"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no symmetry advisory for x/y on the hexagonal 2D fixture"
    PASS=$((PASS+1))
fi

echo "Testing: biaxial direction -- out of scope for v1, advisory must NOT fire (section always"
echo "         prints, but only says 'Skipped')"
rm -rf strain_runs
stb-strain -s si_cubic.fdf -c calc_cubic.fdf --relax-mode cell-fixed --stdir xy --stmin 0 --stmax 1 --no-intro > log_sym_biaxial.txt 2>&1
check_exit_code $? 0
check_contains "Skipped -- symmetry-equivalence check is uniaxial-only" log_sym_biaxial.txt
if grep -q "equivalent to" log_sym_biaxial.txt 2>/dev/null; then
    echo -e "   -> ${RED}Failed:${NC} unexpected symmetry advisory for a biaxial direction"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no symmetry advisory for biaxial directions (v1 scope)"
    PASS=$((PASS+1))
fi

echo "Testing: --help documents --symprec/--angle-tolerance"
stb-strain --help > log_help_sym.txt 2>&1
check_contains "symprec" log_help_sym.txt
check_contains "angle-tolerance" log_help_sym.txt


# --- 6. Bug fix (pre-existing): -s/--structure with a directory component ---
echo -e "\n--- Testing --structure with a subdirectory in the path ---"
rm -rf strain_runs
stb-strain -s subdir/structure.fdf -c calc.fdf --relax-mode cell-fixed --stdir x --stmin 0 --stmax 2 --step 2 --no-intro > log_subdir.txt 2>&1
check_exit_code $? 0
check_success strain_runs/x/strain_x_2.00/structure.fdf


# --- 7. Fractional-coordinate hard check ---
echo -e "\n--- Testing rejection of Cartesian-coordinate input ---"
rm -rf strain_runs
stb-strain -s cartesian_coords.fdf -c calc_cubic.fdf --relax-mode cell-fixed --stdir x --stmin 0 --stmax 1 --no-intro > log_cartesian.txt 2>&1
check_exit_code $? 1
check_contains "Atomic coordinates in this file are Cartesian" log_cartesian.txt
check_contains "Coordinate format        | cartesian" log_cartesian.txt


# --- 8. cell-fixed mode on a calc.fdf with real relaxation directives (MD.VariableCell already T) ---
echo -e "\n--- Testing cell-fixed mode + pseudopotential linking (bundled dojo bank) ---"
rm -rf strain_runs
stb-strain -s si_cubic.fdf -c calc_relaxed_cubic.fdf --relax-mode cell-fixed -p dojo --stdir x --stmin 0 --stmax 2 --step 2 \
    --save-report --no-intro > log_md.txt 2>&1
check_exit_code $? 0
check_contains "\[4\] RELAXATION MODE & CELL CONSTRAINTS" log_md.txt
check_contains "MD.TypeOfRun=CG  Steps: MD.NumCGsteps=100  MD.VariableCell=T" log_md.txt
check_contains "Using bundled pseudopotential bank 'dojo'" log_md.txt
check_contains "\[OK\] All required pseudopotentials found" log_md.txt

echo "Testing: no VariableCell/steps warnings fire (MD.VariableCell already T, 100 steps > 0)"
check_not_contains "\[WARNING\] MD.VariableCell is not enabled" log_md.txt
check_not_contains "\[WARNING\] Relaxation step count is" log_md.txt

echo "Testing: the generated folder's calc.fdf preserves the source's own MD.NumCGsteps/MD.VariableCell"
echo "         untouched -- only the %include config_extra.fdf line is added"
check_success strain_runs/x/strain_x_2.00/calc_relaxed_cubic.fdf
check_contains "MD.NumCGsteps         100" strain_runs/x/strain_x_2.00/calc_relaxed_cubic.fdf
check_contains "MD.VariableCell       T" strain_runs/x/strain_x_2.00/calc_relaxed_cubic.fdf
check_contains "%include config_extra.fdf" strain_runs/x/strain_x_2.00/calc_relaxed_cubic.fdf

echo "Testing: cell-fixed mode's config_extra.fdf fixes all 6 Voigt components"
for i in 1 2 3 4 5 6; do
    check_contains "stress $i  # Fixes" strain_runs/x/strain_x_2.00/config_extra.fdf
done

echo "Testing: the Si pseudopotential was linked into every generated folder"
check_success strain_runs/x/strain_x_0.00/Si.psml
check_success strain_runs/x/strain_x_2.00/Si.psml
if [ -L strain_runs/x/strain_x_2.00/Si.psml ]; then
    echo -e "   -> ${GREEN}Verified:${NC} Si.psml is a symlink (not a copy)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} Si.psml is not a symlink"
    FAIL=$((FAIL+1))
fi

echo "Testing: --save-report wrote strain_stage1.txt (under the direction's own subfolder) with"
echo "         the full numbered report"
check_success strain_runs/x/strain_stage1.txt
check_contains "\[0\] RUN METADATA" strain_runs/x/strain_stage1.txt
check_contains "\[7\] SUMMARY & NEXT STEPS" strain_runs/x/strain_stage1.txt


# --- 9. stress-constrained mode on a 3D bulk fixture (no vacuum) -- only the imposed direction fixed ---
echo -e "\n--- Testing stress-constrained mode on a 3D bulk fixture (si_cubic.fdf, no vacuum) ---"
rm -rf strain_runs
stb-strain -s si_cubic.fdf -c calc_relaxed_cubic.fdf --relax-mode stress-constrained --stdir x --stmin 0 --stmax 2 --step 2 --no-intro > log_stress3d.txt 2>&1
check_exit_code $? 0
check_contains "stress 1  # Fixes XX" strain_runs/x/strain_x_2.00/config_extra.fdf
check_not_contains "stress 2" strain_runs/x/strain_x_2.00/config_extra.fdf
check_not_contains "stress 3" strain_runs/x/strain_x_2.00/config_extra.fdf
check_not_contains "stress 4" strain_runs/x/strain_x_2.00/config_extra.fdf
check_not_contains "stress 5" strain_runs/x/strain_x_2.00/config_extra.fdf
check_not_contains "stress 6" strain_runs/x/strain_x_2.00/config_extra.fdf


# --- 10. MD.VariableCell / step-count advisory warnings ---
echo -e "\n--- Testing MD.VariableCell-not-enabled and zero-relaxation-steps warnings ---"

echo "Testing: calc.fdf with no MD block at all -- both warnings fire"
rm -rf strain_runs
stb-strain -s si_cubic.fdf -c calc_cubic.fdf --relax-mode cell-fixed --stdir x --stmin 0 --stmax 1 --no-intro > log_warn_absent.txt 2>&1
check_exit_code $? 0
check_contains "\[WARNING\] MD.VariableCell is not enabled" log_warn_absent.txt
check_contains "\[WARNING\] Relaxation step count is (absent)" log_warn_absent.txt

echo "Testing: calc.fdf using the real 'MD.Steps' spelling (not MD.NumCGsteps), set to 0 -- steps"
echo "         warning fires (VariableCell is already true, so that warning does NOT fire)"
rm -rf strain_runs
stb-strain -s si_cubic.fdf -c calc_zero_steps.fdf --relax-mode cell-fixed --stdir x --stmin 0 --stmax 1 --no-intro > log_warn_zero.txt 2>&1
check_exit_code $? 0
check_contains "Steps: MD.Steps=0" log_warn_zero.txt
check_contains "\[WARNING\] Relaxation step count is MD.Steps=0" log_warn_zero.txt
check_not_contains "\[WARNING\] MD.VariableCell is not enabled" log_warn_zero.txt


# --- 11. Missing pseudopotential (custom empty folder) ---
echo -e "\n--- Testing missing-pseudopotential warning (custom folder without Si) ---"
mkdir -p empty_bank
rm -rf strain_runs
stb-strain -s si_cubic.fdf -c calc_cubic.fdf --relax-mode cell-fixed -p empty_bank --stdir x --stmin 0 --stmax 1 --no-intro > log_missing_pp.txt 2>&1
check_exit_code $? 0
check_contains "Si.*MISSING" log_missing_pp.txt
check_contains "\[WARNING\] Missing pseudopotential(s) for: Si" log_missing_pp.txt
if ls strain_runs/x/strain_x_1.00/*.psml strain_runs/x/strain_x_1.00/*.psf > /dev/null 2>&1; then
    echo -e "   -> ${RED}Failed:${NC} a pseudopotential file was linked despite none being available"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no pseudopotential file linked when unavailable"
    PASS=$((PASS+1))
fi

echo "Testing: nonexistent pseudo bank/path is rejected cleanly"
stb-strain -s si_cubic.fdf -c calc_cubic.fdf --relax-mode cell-fixed -p does_not_exist_bank --stdir x --stmin 0 --stmax 1 --no-intro \
    > log_bad_pseudo.txt 2>&1
check_exit_code $? 1
check_contains "not a recognized pseudopotential bank" log_bad_pseudo.txt


# --- 12. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: missing --structure file"
stb-strain -s does_not_exist.fdf -c calc.fdf --relax-mode cell-fixed --stdir x --no-intro > log_missing_structure.txt 2>&1
check_exit_code $? 1
check_contains "\[ERROR\] Structure file 'does_not_exist.fdf' not found" log_missing_structure.txt

echo "Testing: missing --calc file"
stb-strain -s structure.fdf -c does_not_exist.fdf --relax-mode cell-fixed --stdir x --no-intro > log_missing_calc.txt 2>&1
check_exit_code $? 1
check_contains "\[ERROR\] Calc file 'does_not_exist.fdf' not found" log_missing_calc.txt

echo "Testing: --stmin > --stmax"
stb-strain -s structure.fdf -c calc.fdf --relax-mode cell-fixed --stdir x --stmin 10 --stmax 5 --no-intro > log_bad_range.txt 2>&1
check_exit_code $? 1
check_contains "\[ERROR\] Minimum strain cannot be greater than maximum strain" log_bad_range.txt

echo "Testing: invalid --stdir"
stb-strain -s structure.fdf -c calc.fdf --relax-mode cell-fixed --stdir w --no-intro > log_bad_dir.txt 2>&1
check_exit_code $? 1
check_contains "\[FAIL\] Invalid direction" log_bad_dir.txt

echo "Testing: --step 0 (must not crash with a raw ZeroDivisionError traceback)"
stb-strain -s structure.fdf -c calc.fdf --relax-mode cell-fixed --stdir x --step 0 --no-intro > log_zero_step.txt 2>&1
check_exit_code $? 1
check_contains "\[ERROR\] Step cannot be zero" log_zero_step.txt

echo "Testing: missing required args (both --stdir and --relax-mode omitted)"
stb-strain -s structure.fdf -c calc.fdf --no-intro > log_missing_args.txt 2>&1
check_exit_code $? 2

echo "Testing: --relax-mode omitted (--stdir given) -- still rejected by argparse"
stb-strain -s structure.fdf -c calc.fdf --stdir x --no-intro > log_missing_relaxmode.txt 2>&1
check_exit_code $? 2

echo "Testing: invalid --relax-mode value"
stb-strain -s structure.fdf -c calc.fdf --relax-mode bogus --stdir x --no-intro > log_bad_relaxmode.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-strain --version > log_version.txt 2>&1
check_contains "stb-strain" log_version.txt

echo "Testing: --help documents structure/calc/relax-mode/pseudo-dir/uniaxial/biaxial/output-dir/vacuum-gap/save-report"
stb-strain --help > log_help.txt 2>&1
check_contains "\-\-structure" log_help.txt
check_contains "\-\-calc" log_help.txt
check_contains "relax-mode" log_help.txt
check_contains "cell-fixed" log_help.txt
check_contains "stress-constrained" log_help.txt
check_contains "pseudo-dir" log_help.txt
check_contains "biaxial" log_help.txt
check_contains "output-dir" log_help.txt
check_contains "vacuum-gap" log_help.txt
check_contains "save-report" log_help.txt


# --- 13. Interactive path (stb-suite, shortcut 4.1.1) ---
# NOTE: vacuum-padded axes and symmetry-equivalent directions (uniaxial AND
# biaxial) are now detected up front and offered as a NUMBERED menu (see
# strain.py::print_direction_selection_table) -- the direction is chosen by
# picking a number, not typing 'x'/'xy' as free text, and there's no more
# separate "also generate the equivalent axis/axes too?" y/n follow-up.
# Individual independent directions are listed first, followed by grouped
# "run several at once" choices: ALL UNIAXIAL (only shown if there's at least
# one independent uniaxial direction), ALL BIAXIAL (same, biaxial), and
# finally ALL (uniaxial + biaxial together), always last. On structure.fdf
# (2D, vacuum along z, point group -6m2, x/y NOT symmetry-equivalent -- see
# section 5 above), the independent/non-vacuum menu is: 1) x  2) y  3) xy
# 4) ALL UNIAXIAL (x, y)  5) ALL BIAXIAL (xy)  6) ALL (x, y, xy). On
# si_cubic.fdf (m-3m, all 3 axes equivalent), it collapses to: 1) x  2) xy
# 3) ALL UNIAXIAL (x)  4) ALL BIAXIAL (xy)  5) ALL (x, xy).
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.1.1) ---"

echo "Testing: navigate 4.1.1 -> structure -> calc -> relax-mode (1=cell-fixed) -> skip pseudo -> skip"
echo "         advanced -> direction menu, pick 1 (x) -> 0/2/2 -> no report -> quit"
rm -rf strain_runs
# Prompts in order: structure file (default suggestion, blank -> structure.fdf),
# calc file (default suggestion, blank -> calc.fdf), relax mode (1 -> cell-fixed),
# pseudo source (blank -> skip), advanced settings (n -> skip), direction menu
# (1 -> x), stmin, stmax, step, save-report (n), then the "Press Enter to
# continue" pause, then quit.
printf '4.1.1\n\n\n1\n\nn\n1\n0\n2\n2\nn\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_contains "Detected dimensionality: 2D" log_menu.txt
check_contains "STRAIN DIRECTION SELECTION" log_menu.txt
check_contains "CONFIGURATION SUMMARY" log_menu.txt
check_contains "Dimensionality.*2D" log_menu.txt
check_contains "Relax mode.*cell-fixed" log_menu.txt
check_contains "Pseudo source.*not given" log_menu.txt
check_contains "Save report.*no" log_menu.txt
check_success strain_runs/x/strain_x_2.00/structure.fdf
check_success strain_runs/x/strain_x_2.00/calc.fdf
check_success strain_runs/x/strain_x_2.00/config_extra.fdf
check_absent strain_runs/x/strain_stage1.txt

echo "Testing: navigate 4.1.1 -> invalid structure file then valid"
rm -rf strain_runs
printf '4.1.1\ndoes_not_exist.fdf\nstructure.fdf\ncalc.fdf\n1\n\nn\n1\n0\n2\n2\nn\n\n0\n' | stb-suite > log_menu_retry.txt 2>&1
check_contains "File not found" log_menu_retry.txt
check_success strain_runs/x/strain_x_2.00/structure.fdf

echo "Testing: navigate 4.1.1 -> direction menu, pick the biaxial 'xy' entry -> generated correctly"
rm -rf strain_runs
printf '4.1.1\nstructure.fdf\ncalc.fdf\n1\n\nn\n3\n0\n2\n2\nn\n\n0\n' | stb-suite > log_menu_biaxial.txt 2>&1
check_contains "Detected dimensionality: 2D" log_menu_biaxial.txt
check_contains "Dimensionality.*2D" log_menu_biaxial.txt
check_contains "xy.*biaxial.*INDEPENDENT" log_menu_biaxial.txt
check_contains "z.*uniaxial.*VACUUM" log_menu_biaxial.txt
check_success strain_runs/xy/strain_xy_2.00/structure.fdf

echo "Testing: navigate 4.1.1 -> invalid relax-mode choice then valid"
rm -rf strain_runs
printf '4.1.1\nstructure.fdf\ncalc.fdf\n9\n1\n\nn\n1\n0\n2\n2\nn\n\n0\n' | stb-suite > log_menu_badmode.txt 2>&1
check_contains "Invalid choice" log_menu_badmode.txt
check_success strain_runs/x/strain_x_2.00/structure.fdf

echo "Testing: navigate 4.1.1 -> invalid direction-menu choice then valid"
rm -rf strain_runs
printf '4.1.1\nstructure.fdf\ncalc.fdf\n1\n\nn\n9\n1\n0\n2\n2\nn\n\n0\n' | stb-suite > log_menu_baddir.txt 2>&1
check_contains "Invalid choice! Enter a number from 1 to 6" log_menu_baddir.txt
check_success strain_runs/x/strain_x_2.00/structure.fdf

echo "Testing: navigate 4.1.1 -> direction menu, pick ALL UNIAXIAL (x, y) -> only x/y generated, not xy"
rm -rf strain_runs
printf '4.1.1\nstructure.fdf\ncalc.fdf\n1\n\nn\n4\n0\n2\n2\nn\n\n0\n' | stb-suite > log_menu_alluni.txt 2>&1
check_contains "Direction(s) to generate.*x, y" log_menu_alluni.txt
check_success strain_runs/x/strain_x_2.00/structure.fdf
check_success strain_runs/y/strain_y_2.00/structure.fdf
check_absent strain_runs/xy/strain_xy_2.00/structure.fdf

echo "Testing: navigate 4.1.1 -> direction menu, pick ALL BIAXIAL (xy) -> only xy generated, not x/y"
rm -rf strain_runs
printf '4.1.1\nstructure.fdf\ncalc.fdf\n1\n\nn\n5\n0\n2\n2\nn\n\n0\n' | stb-suite > log_menu_allbi.txt 2>&1
check_contains "Direction(s) to generate.*xy" log_menu_allbi.txt
check_success strain_runs/xy/strain_xy_2.00/structure.fdf
check_absent strain_runs/x/strain_x_2.00/structure.fdf
check_absent strain_runs/y/strain_y_2.00/structure.fdf

echo "Testing: navigate 4.1.1 -> cubic fixture -> relax-mode (2=stress-constrained) -> pseudo bank 1 (dojo)"
echo "         -> direction menu, pick ALL (independent directions only: x, xy -- NOT y, z, xz, yz,"
echo "         which are symmetry-redundant) -> save report"
rm -rf strain_runs
# Prompts in order: structure file, calc file, relax mode (2 -- stress
# -constrained), pseudo source (1 -- bundled dojo), advanced settings (n),
# direction menu (5 -- 'ALL' entry, the last of 5 entries: x, xy, ALL
# UNIAXIAL, ALL BIAXIAL, ALL, on this cubic fixture), stmin, stmax, step,
# save-report (y), then one "Press Enter to continue" pause (not two --
# run_tool's pause is suppressed for every direction but the last), then quit.
printf '4.1.1\nsi_cubic.fdf\ncalc_cubic.fdf\n2\n1\nn\n5\n0\n2\n2\ny\n\n0\n' | stb-suite > log_menu_multi.txt 2>&1
check_contains "2 independent, non-vacuum direction(s): x, xy" log_menu_multi.txt
check_contains "Direction(s) to generate.*x, xy" log_menu_multi.txt
check_contains "Relax mode.*stress-constrained" log_menu_multi.txt
check_success strain_runs/x/strain_x_2.00/si_cubic.fdf
check_success strain_runs/xy/strain_xy_2.00/si_cubic.fdf
check_success strain_runs/x/strain_x_2.00/Si.psml
echo "Testing: each direction's own strain_stage1.txt was saved separately (no clobbering)"
check_success strain_runs/x/strain_stage1.txt
check_success strain_runs/xy/strain_stage1.txt
if [ "$(grep -c "Press Enter to continue" log_menu_multi.txt)" -eq 1 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} only one 'Press Enter to continue' pause for all directions"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} expected exactly 1 'Press Enter to continue' pause, "\
"got $(grep -c "Press Enter to continue" log_menu_multi.txt)"
    FAIL=$((FAIL+1))
fi


popd > /dev/null

# --- 14. Summary ---
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
