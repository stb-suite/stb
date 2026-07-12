#!/bin/bash

# --- Setup ---
# Smoke test for stb-strainAnalysis (Stress-Strain Analysis, item 4.1.2)
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
echo "--- Starting tester for STB-StrainAnalysis (item 4.1.2) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
for d in strain_xx_1.00 strain_xx_2.00 strain_xx_m1.00; do
    mkdir -p "$TEST_DIR/$d"
    cp "$FIXTURE_DIR/$d/calc.out" "$TEST_DIR/$d/"
done
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Uniaxial analysis (existing checked-in fixtures, --dir . since they
#     sit directly here rather than under a strain_runs/ wrapper) ---
echo -e "\n--- Testing analysis of a uniaxial (xx) strain sweep ---"
stb-strainAnalysis --file calc.out --dir . -o mech.dat --no-intro > log_uniaxial.txt 2>&1
check_success mech.dat
check_success mechanical_report.txt
check_contains "Initial Slope" mechanical_report.txt
check_contains "1000.0000 GPa" mechanical_report.txt
check_contains "Peak Stress" mechanical_report.txt
check_contains "20.0000 GPa" mechanical_report.txt
check_contains "Critical Strain: 2.0000 %" mechanical_report.txt
# Fixture only goes up to 2% strain, so the peak sits at the edge of the
# tested range -- the boundary warning must fire.
check_contains "WARNING: Peak occurred at the edge of the tested strain range" mechanical_report.txt
# Yield must NOT appear by default (opt-in via --yield now).
if grep -q "Yield" mechanical_report.txt 2>/dev/null; then
    echo -e "   -> ${RED}Failed:${NC} Yield line present without --yield"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no Yield line without --yield"
    PASS=$((PASS+1))
fi
# Transverse-stress diagnostic must be present for a uniaxial 3D direction.
check_contains "Transverse Stress @ peak" mechanical_report.txt
# Uniaxial 'xx' must NOT trigger the biaxial warning (regression check for the
# len(direction)==2 vs "same axis repeated" distinction).
if grep -q "Biaxial direction" log_uniaxial.txt 2>/dev/null; then
    echo -e "   -> ${RED}Failed:${NC} unexpected biaxial WARNING for uniaxial 'xx' direction"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no biaxial warning for uniaxial 'xx' direction"
    PASS=$((PASS+1))
fi

echo -e "\n--- Testing --yield opt-in ---"
stb-strainAnalysis --file calc.out --dir . -o mech.dat --yield --no-intro > log_yield.txt 2>&1
check_contains "Yield Stress (0.2%)" log_yield.txt
check_contains "macroscopic-plasticity concept" log_yield.txt

echo -e "\n--- Testing the companion .gplot/.dat plot output ---"
check_success mech.gplot
check_contains 'set ylabel "Stress (GPa)"' mech.gplot
check_contains 'index 1 using 2:3' mech.gplot
if command -v gnuplot > /dev/null 2>&1; then
    gnuplot mech.gplot > log_gnuplot.txt 2>&1
    check_exit_code $? 0
    check_success mech.pdf
else
    echo -e "   -> ${YELLOW}Skipped:${NC} gnuplot not installed, skipping actual render"
fi


# --- 3. Biaxial analysis (synthetic fixtures, same minimal style already
#     used by 1-strain/2-elastic/3-cohesive: only the parsed line matters).
#     Run in its own subdirectory -- strain_analysis.py scans ALL 'strain_*'
#     folders inside --dir, so this must not see the uniaxial 'strain_xx_*'
#     dirs from step 2. ---
echo -e "\n--- Testing analysis of a biaxial (xy) strain sweep ---"
mkdir -p biaxial_run
pushd biaxial_run > /dev/null
mkdir -p strain_xy_m1.00 strain_xy_0.00 strain_xy_1.00 strain_xy_2.00
for pair in "strain_xy_m1.00:-100" "strain_xy_0.00:0" "strain_xy_1.00:100" "strain_xy_2.00:200"; do
    folder="${pair%%:*}"
    xx="${pair##*:}"
    cat > "$folder/calc.out" <<EOF
outcell: Unit cell vectors (Ang):
        5.000000    0.000000    0.000000
        0.000000    5.000000    0.000000
        0.000000    0.000000   20.000000

siesta: Stress tensor Voigt (kbar):     ${xx}.00      ${xx}.00      0.00      0.00      0.00      0.00
EOF
done
stb-strainAnalysis --file calc.out --dir . -o mech_biaxial.dat --no-intro > log_biaxial.txt 2>&1
check_contains "WARNING" log_biaxial.txt
check_contains "Biaxial direction 'xy' detected" log_biaxial.txt
check_contains "first axis ('xx')" log_biaxial.txt
check_success mech_biaxial.dat
# Transverse-stress diagnostic is skipped for a genuine biaxial direction.
if grep -q "Transverse Stress" mechanical_report.txt 2>/dev/null; then
    echo -e "   -> ${RED}Failed:${NC} transverse-stress diagnostic unexpectedly shown for biaxial direction"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no transverse-stress diagnostic for biaxial direction"
    PASS=$((PASS+1))
fi
popd > /dev/null


# --- 4. Multi-direction handling (--compare), synthetic fixtures with 2
#     different strain directions ('strain_x_*' + 'strain_y_*') in the same
#     directory -- own subdirectory, same reason as step 3. ---
echo -e "\n--- Testing multi-direction detection and --compare ---"
mkdir -p compare_run
pushd compare_run > /dev/null
for pair in "strain_x_m1.00:-80" "strain_x_0.00:0" "strain_x_1.00:80" "strain_x_2.00:160"; do
    folder="${pair%%:*}"
    val="${pair##*:}"
    mkdir -p "$folder"
    cat > "$folder/calc.out" <<EOF
outcell: Unit cell vectors (Ang):
        5.000000    0.000000    0.000000
        0.000000    5.000000    0.000000
        0.000000    0.000000   20.000000

siesta: Stress tensor Voigt (kbar):     ${val}.00      0.00      0.00      0.00      0.00      0.00
EOF
done
for pair in "strain_y_m1.00:-60" "strain_y_0.00:0" "strain_y_1.00:60" "strain_y_2.00:120"; do
    folder="${pair%%:*}"
    val="${pair##*:}"
    mkdir -p "$folder"
    cat > "$folder/calc.out" <<EOF
outcell: Unit cell vectors (Ang):
        5.000000    0.000000    0.000000
        0.000000    5.000000    0.000000
        0.000000    0.000000   20.000000

siesta: Stress tensor Voigt (kbar):     0.00      ${val}.00      0.00      0.00      0.00      0.00
EOF
done

echo "Testing: multiple directions without --compare -> clean error"
stb-strainAnalysis --file calc.out --dir . --no-intro > log_nocompare.txt 2>&1
check_exit_code $? 1
check_contains "Multiple strain directions found (x, y)" log_nocompare.txt

echo "Testing: multiple directions with --compare -> comparison table"
stb-strainAnalysis --file calc.out --dir . --compare -o mech.dat --no-intro > log_compare.txt 2>&1
check_exit_code $? 0
check_contains "MECHANICAL STRESS COMPARISON" log_compare.txt
check_success mechanical_comparison.txt
check_contains "800.00" mechanical_comparison.txt
check_contains "600.00" mechanical_comparison.txt
check_success mech_comparison.dat
check_success mech_comparison.gplot
if command -v gnuplot > /dev/null 2>&1; then
    gnuplot mech_comparison.gplot > log_gnuplot_compare.txt 2>&1
    check_exit_code $? 0
    check_success mech_comparison.pdf
else
    echo -e "   -> ${YELLOW}Skipped:${NC} gnuplot not installed, skipping actual render"
fi
popd > /dev/null


# --- 5. 1D mode (--1d / --cross-section), synthetic nanotube-like fixture:
#     large in-plane vacuum (30x30 Ang), single periodic axis (z). ---
echo -e "\n--- Testing --1d (Force in nN, auto cross-section) ---"
mkdir -p oned_run
pushd oned_run > /dev/null
for pair in "strain_z_m1.00:-50" "strain_z_0.00:0" "strain_z_1.00:50" "strain_z_2.00:100"; do
    folder="${pair%%:*}"
    zz="${pair##*:}"
    mkdir -p "$folder"
    cat > "$folder/calc.out" <<EOF
outcell: Unit cell vectors (Ang):
       30.000000    0.000000    0.000000
        0.000000   30.000000    0.000000
        0.000000    0.000000   10.000000

siesta: Stress tensor Voigt (kbar):     0.00      0.00      ${zz}.00      0.00      0.00      0.00
EOF
done

stb-strainAnalysis --file calc.out --dir . --1d -o mech1d.dat --no-intro > log_1d.txt 2>&1
check_exit_code $? 0
# Cell cross-section = |30x0x0 . 0x30x0| = 900 Ang^2; 100 kBar * 900 * 1e-3 = 90 nN.
check_contains "Peak Force" log_1d.txt
check_contains "90.00 nN" log_1d.txt
check_contains "MECHANICAL FORCE REPORT" log_1d.txt

echo "Testing: --1d --cross-section (matches cell area -> recovers plain 10 GPa)"
stb-strainAnalysis --file calc.out --dir . --1d --cross-section 900 -o mech1d.dat --no-intro > log_1d_xsec.txt 2>&1
check_contains "Peak Stress (conventional, cross-section=900.00" log_1d_xsec.txt
check_contains "10.00 GPa" log_1d_xsec.txt

echo "Testing: --1d gnuplot output"
if command -v gnuplot > /dev/null 2>&1; then
    gnuplot mech1d.gplot > log_gnuplot_1d.txt 2>&1
    check_exit_code $? 0
    check_success mech1d.pdf
else
    echo -e "   -> ${YELLOW}Skipped:${NC} gnuplot not installed, skipping actual render"
fi

echo "Testing: --1d and --2d are mutually exclusive"
stb-strainAnalysis --file calc.out --dir . --1d --2d --no-intro > log_1d2d.txt 2>&1
check_exit_code $? 1
check_contains "mutually exclusive" log_1d2d.txt

echo "Testing: --cross-section without --1d is rejected"
stb-strainAnalysis --file calc.out --dir . --cross-section 100 --no-intro > log_xsec_only.txt 2>&1
check_exit_code $? 1
check_contains "only applies together with --1d" log_xsec_only.txt
popd > /dev/null


# --- 6. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: --dir not found"
stb-strainAnalysis --file calc.out --dir does_not_exist --no-intro > log_missing_dir.txt 2>&1
check_exit_code $? 1
check_contains "Directory 'does_not_exist' not found" log_missing_dir.txt

echo "Testing: no strain_* folders"
mkdir -p empty_runs
stb-strainAnalysis --file calc.out --dir empty_runs --no-intro > log_empty.txt 2>&1
check_exit_code $? 1
check_contains "No 'strain_\*' folders found" log_empty.txt

echo "Testing: strain_* folders present but no matching output file"
mkdir -p missing_out/strain_x_1.00
stb-strainAnalysis --file does_not_exist.out --dir missing_out --no-intro > log_missing_out.txt 2>&1
check_exit_code $? 1
check_contains "No stress data found" log_missing_out.txt

echo "Testing: only 1 valid strain step -- must not silently report NaN"
mkdir -p single_point/strain_x_1.00
cat > single_point/strain_x_1.00/calc.out <<EOF
outcell: Unit cell vectors (Ang):
        5.000000    0.000000    0.000000
        0.000000    5.000000    0.000000
        0.000000    0.000000   20.000000

siesta: Stress tensor Voigt (kbar):     100.00      0.00      0.00      0.00      0.00      0.00
EOF
stb-strainAnalysis --file calc.out --dir single_point --no-intro > log_single_point.txt 2>&1
check_exit_code $? 1
check_contains "at least 2 are needed to fit a slope" log_single_point.txt

echo "Testing: unrecognized direction folder -- must warn, not silently default to xx"
mkdir -p weird_dir/strain_ab_1.00 weird_dir/strain_ab_2.00
for s in 1.00:50 2.00:100; do
    val="${s%%:*}"; kbar="${s##*:}"
    cat > "weird_dir/strain_ab_$val/calc.out" <<EOF
outcell: Unit cell vectors (Ang):
        5.000000    0.000000    0.000000
        0.000000    5.000000    0.000000
        0.000000    0.000000   20.000000

siesta: Stress tensor Voigt (kbar):     ${kbar}.00      0.00      0.00      0.00      0.00      0.00
EOF
done
stb-strainAnalysis --file calc.out --dir weird_dir --no-intro > log_weird_dir.txt 2>&1
check_exit_code $? 0
check_contains "Unrecognized direction 'ab'" log_weird_dir.txt

echo "Testing: 'xx' and 'x' folder-naming conventions are treated as the SAME direction"
mkdir -p canon_check/strain_xx_1.00 canon_check/strain_x_2.00
cat > canon_check/strain_xx_1.00/calc.out <<EOF
outcell: Unit cell vectors (Ang):
        5.000000    0.000000    0.000000
        0.000000    5.000000    0.000000
        0.000000    0.000000   20.000000

siesta: Stress tensor Voigt (kbar):     50.00      0.00      0.00      0.00      0.00      0.00
EOF
cat > canon_check/strain_x_2.00/calc.out <<EOF
outcell: Unit cell vectors (Ang):
        5.000000    0.000000    0.000000
        0.000000    5.000000    0.000000
        0.000000    0.000000   20.000000

siesta: Stress tensor Voigt (kbar):     100.00      0.00      0.00      0.00      0.00      0.00
EOF
stb-strainAnalysis --file calc.out --dir canon_check --no-intro > log_canon.txt 2>&1
check_exit_code $? 0
if grep -q "Multiple strain directions found" log_canon.txt 2>/dev/null; then
    echo -e "   -> ${RED}Failed:${NC} 'xx' and 'x' were treated as 2 different directions"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} 'xx' and 'x' folders were grouped as the same direction"
    PASS=$((PASS+1))
fi

echo "Testing: --version"
stb-strainAnalysis --version > log_version.txt 2>&1
check_contains "stb-strainAnalysis" log_version.txt

echo "Testing: --help documents --2d/--1d/--cross-section/--yield/--dir"
stb-strainAnalysis --help > log_help.txt 2>&1
check_contains "2D" log_help.txt
check_contains "1D" log_help.txt
check_contains "cross-section" log_help.txt
check_contains "yield" log_help.txt
check_contains "strain_runs" log_help.txt


# --- 7. Interactive path (stb-suite, shortcut 4.1.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.1.2) ---"

echo "Testing: navigate 4.1.2 -> dir '.' -> calc.out -> 3D (default) -> quit"
rm -f mechanical_curve.dat mechanical_report.txt
printf '4.1.2\n.\ncalc.out\n1\nn\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_contains "Initial Slope" log_menu.txt
check_success mechanical_curve.dat
if grep -q "Traceback" log_menu.txt 2>/dev/null; then
    echo -e "   -> ${RED}Failed:${NC} unexpected Python traceback in interactive session (log_menu.txt)"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} interactive session exited cleanly, no traceback"
    PASS=$((PASS+1))
fi


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
