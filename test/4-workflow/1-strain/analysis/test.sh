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
#     sit directly here rather than under a strain_runs/ wrapper). No .fdf
#     structure file is present here -- exercises the "no structure file
#     found, fall back to 3D" auto-detection path (the common/expected case
#     for a hand-built/legacy folder). ---
echo -e "\n--- Testing analysis of a uniaxial (xx) strain sweep ---"
stb-strainAnalysis --file calc.out --dir . --save-report --save-gnuplot --no-intro > log_uniaxial.txt 2>&1
check_success stb_strainAnalysis_report.txt
check_contains "No parseable .fdf structure file found" stb_strainAnalysis_report.txt
check_contains "Dimensionality    : 3D (auto-detected)" stb_strainAnalysis_report.txt
check_contains "Initial Slope" stb_strainAnalysis_report.txt
check_contains "1000.0000 GPa" stb_strainAnalysis_report.txt
check_contains "Peak Stress" stb_strainAnalysis_report.txt
check_contains "20.0000 GPa" stb_strainAnalysis_report.txt
check_contains "Critical Strain : 2.0000 %" stb_strainAnalysis_report.txt
# Fixture only goes up to 2% strain, so the peak sits at the edge of the
# tested range -- the boundary warning must fire.
check_contains "Peak occurred at the edge of the tested strain range" stb_strainAnalysis_report.txt
# Yield must NOT appear by default (opt-in via --yield now).
if grep -q "Yield" stb_strainAnalysis_report.txt 2>/dev/null; then
    echo -e "   -> ${RED}Failed:${NC} Yield line present without --yield"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no Yield line without --yield"
    PASS=$((PASS+1))
fi
# Transverse-stress diagnostic must be present for a uniaxial 3D direction.
check_contains "Transverse Stress @ peak" stb_strainAnalysis_report.txt
# Uniaxial 'xx' must NOT trigger the biaxial warning (regression check for the
# len(direction)==2 vs "same axis repeated" distinction).
if grep -q "Biaxial direction" log_uniaxial.txt 2>/dev/null; then
    echo -e "   -> ${RED}Failed:${NC} unexpected biaxial WARNING for uniaxial 'xx' direction"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no biaxial warning for uniaxial 'xx' direction"
    PASS=$((PASS+1))
fi

echo -e "\n--- Testing the initial-slope fit window (coarse-step / large-strain sweep) ---"
echo "Testing: only 1 point falls within +-2% (10% step size) -- must widen to the nearest"
echo "         few points, NOT silently fit the entire (nonlinear/plateauing) sweep"
mkdir -p coarse_run
pushd coarse_run > /dev/null
# strain(%) : stress(kbar, xx) -- perfectly linear 0->20%, then a plateau
# (yielded/plastic) from 20% on. A fit using ONLY the nearest 3 points
# (0/10/20%) is exact (slope=2000 kbar/frac=200 GPa, R^2=1.0); a fit using
# the ENTIRE sweep (the old, buggy fallback) is pulled down by the plateau
# points to a materially different, much worse-fitting line -- this is
# exactly the real bug found live on a real 0-40% SIESTA sweep (R^2 as low
# as ~0.01 using the old fallback).
for pair in "0.00:0" "10.00:200" "20.00:400" "30.00:400" "40.00:400"; do
    val="${pair%%:*}"
    kbar="${pair##*:}"
    mkdir -p "strain_x_$val"
    cat > "strain_x_$val/calc.out" <<EOF
outcell: Unit cell vectors (Ang):
        5.000000    0.000000    0.000000
        0.000000    5.000000    0.000000
        0.000000    0.000000   20.000000

siesta: Stress tensor Voigt (kbar):     ${kbar}.00      0.00      0.00      0.00      0.00      0.00
EOF
done
stb-strainAnalysis --file calc.out --dir . --no-intro > log_coarse.txt 2>&1
check_contains "200.0000 GPa" log_coarse.txt
check_contains "R²=1.0000" log_coarse.txt
check_contains "Fit window    : 3 point(s), strain 0.0000% to 20.0000%" log_coarse.txt
check_contains "Fewer than 3 strain steps fall within +-2%" log_coarse.txt
if grep -q "100.0000 GPa" log_coarse.txt 2>/dev/null; then
    echo -e "   -> ${RED}Failed:${NC} modulus matches the OLD buggy full-range fit (100.0000 GPa) -- fit window fallback regressed"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} modulus does not match the old buggy full-range fit"
    PASS=$((PASS+1))
fi
popd > /dev/null


echo -e "\n--- Testing --yield opt-in ---"
stb-strainAnalysis --file calc.out --dir . --yield --no-intro > log_yield.txt 2>&1
check_contains "Yield Stress (0.2%)" log_yield.txt
check_contains "macroscopic-plasticity concept" log_yield.txt

echo -e "\n--- Testing the companion .gplot/.dat plot output (--save-gnuplot) ---"
check_success x_curve.dat
check_success x_curve.gplot
check_contains 'set ylabel "Stress (GPa)"' x_curve.gplot
check_contains 'index 1 using 2:3' x_curve.gplot
if command -v gnuplot > /dev/null 2>&1; then
    gnuplot x_curve.gplot > log_gnuplot.txt 2>&1
    check_exit_code $? 0
    check_success x_curve.pdf
else
    echo -e "   -> ${YELLOW}Skipped:${NC} gnuplot not installed, skipping actual render"
fi

echo -e "\n--- Testing --save-gnuplot with a relative -o/--output-dir (not '.') ---"
echo "Testing: the .gplot script must reference its data/output files by bare filename,"
echo "         NOT prefixed with --output-dir -- gnuplot resolves them relative to ITS OWN"
echo "         cwd when run, and the natural usage is 'cd <output-dir> && gnuplot <script>'"
rm -rf gnuplot_subdir_out
stb-strainAnalysis --file calc.out --dir . -o gnuplot_subdir_out --save-gnuplot --no-intro > log_gnuplot_outdir.txt 2>&1
check_success gnuplot_subdir_out/x_curve.dat
check_success gnuplot_subdir_out/x_curve.gplot
if grep -q "gnuplot_subdir_out" gnuplot_subdir_out/x_curve.gplot 2>/dev/null; then
    echo -e "   -> ${RED}Failed:${NC} x_curve.gplot references its own --output-dir path -- will"
    echo -e "      fail when run from inside that same directory (doubly-nested path)"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} x_curve.gplot references bare filenames only"
    PASS=$((PASS+1))
fi
check_contains 'set output "x_curve.pdf"' gnuplot_subdir_out/x_curve.gplot
check_contains 'plot "x_curve.dat"' gnuplot_subdir_out/x_curve.gplot
check_contains "Render with: cd gnuplot_subdir_out" log_gnuplot_outdir.txt
if command -v gnuplot > /dev/null 2>&1; then
    (cd gnuplot_subdir_out && gnuplot x_curve.gplot > gnuplot_render.log 2>&1)
    check_exit_code $? 0
    check_success gnuplot_subdir_out/x_curve.pdf
else
    echo -e "   -> ${YELLOW}Skipped:${NC} gnuplot not installed, skipping actual render"
fi

echo "Testing: without --save-report/--save-gnuplot, nothing is written"
rm -f x_curve.dat x_curve.gplot stb_strainAnalysis_report.txt
stb-strainAnalysis --file calc.out --dir . --no-intro > log_nosave.txt 2>&1
check_contains "Not written (off by default" log_nosave.txt
if [ -s x_curve.dat ] || [ -s stb_strainAnalysis_report.txt ]; then
    echo -e "   -> ${RED}Failed:${NC} a file was written despite no --save-* flag"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no file written without --save-report/--save-gnuplot"
    PASS=$((PASS+1))
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
stb-strainAnalysis --file calc.out --dir . --save-report --save-gnuplot --no-intro > log_biaxial.txt 2>&1
check_contains "WARNING" log_biaxial.txt
check_contains "Biaxial direction 'xy' detected" log_biaxial.txt
check_contains "first axis ('xx')" log_biaxial.txt
check_success xy_curve.dat
# Transverse-stress diagnostic is skipped for a genuine biaxial direction --
# checked against THIS run's own report file (own subdirectory, so no risk
# of reading a leftover report from a different run).
check_success stb_strainAnalysis_report.txt
if grep -q "Transverse Stress" stb_strainAnalysis_report.txt 2>/dev/null; then
    echo -e "   -> ${RED}Failed:${NC} transverse-stress diagnostic unexpectedly shown for biaxial direction"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no transverse-stress diagnostic for biaxial direction"
    PASS=$((PASS+1))
fi
popd > /dev/null


# --- 4. Multi-direction handling, synthetic fixtures with 2 different
#     strain directions ('strain_x_*' + 'strain_y_*') in the same directory
#     -- own subdirectory, same reason as step 3. Comparison now happens
#     AUTOMATICALLY whenever more than 1 direction is found -- no --compare
#     flag needed. ---
echo -e "\n--- Testing automatic multi-direction comparison ---"
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

echo "Testing: 2 directions found under --dir -> compared automatically, no flag needed"
stb-strainAnalysis --file calc.out --dir . --save-report --save-gnuplot --no-intro > log_compare.txt 2>&1
check_exit_code $? 0
check_contains "2 direction(s) found -- comparing automatically" log_compare.txt
check_success stb_strainAnalysis_report.txt
check_contains "800.0000" stb_strainAnalysis_report.txt
check_contains "600.0000" stb_strainAnalysis_report.txt
check_success comparison_curve.dat
check_success comparison_curve.gplot
if command -v gnuplot > /dev/null 2>&1; then
    gnuplot comparison_curve.gplot > log_gnuplot_compare.txt 2>&1
    check_exit_code $? 0
    check_success comparison_curve.pdf
else
    echo -e "   -> ${YELLOW}Skipped:${NC} gnuplot not installed, skipping actual render"
fi
popd > /dev/null


# --- 4b. Nested per-direction layout (stb-strain's own default output layout:
#     one <direction>/ subfolder per strain direction, each holding that
#     direction's own strain_<direction>_<pct> folders) -- own subdirectory,
#     same reason as step 3/4. ---
echo -e "\n--- Testing the nested <direction>/strain_<direction>_<pct> layout ---"
mkdir -p nested_run/x nested_run/y
pushd nested_run > /dev/null
for pair in "x/strain_x_0.00:0" "x/strain_x_1.00:80" "x/strain_x_2.00:160"; do
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
for pair in "y/strain_y_0.00:0" "y/strain_y_1.00:60" "y/strain_y_2.00:120"; do
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

echo "Testing: --dir pointed at a single direction's own subfolder -- that direction alone"
stb-strainAnalysis --file calc.out --dir x --no-intro > log_nested_single.txt 2>&1
check_exit_code $? 0
check_contains "MECHANICAL PROPERTIES" log_nested_single.txt
check_contains "Direction       : X" log_nested_single.txt

echo "Testing: --dir at the top level -> every nested direction compared automatically"
stb-strainAnalysis --file calc.out --dir . --no-intro > log_nested_compare.txt 2>&1
check_exit_code $? 0
check_contains "2 direction(s) found -- comparing automatically" log_nested_compare.txt
check_contains "800.0000" log_nested_compare.txt
check_contains "600.0000" log_nested_compare.txt
popd > /dev/null


# --- 4c. Auto-detected dimensionality from a REAL .fdf structure file (the
#     genuine Stage-1-generated case -- steps 2-4b above all lack a .fdf,
#     so they only ever exercise the "no structure file found, fall back to
#     3D" path; these fixtures exercise real vacuum-axis detection). ---
echo -e "\n--- Testing auto dimensionality detection from a real .fdf structure file ---"

echo "Testing: 2D (vacuum along z only) -- no --dimensionality flag given"
mkdir -p auto2d_run/x/strain_x_0.00 auto2d_run/x/strain_x_1.00 auto2d_run/x/strain_x_2.00
for pair in "0.00:0" "1.00:40" "2.00:80"; do
    val="${pair%%:*}"
    kbar="${pair##*:}"
    folder="auto2d_run/x/strain_x_$val"
    cat > "$folder/calc.out" <<EOF
outcell: Unit cell vectors (Ang):
        5.900000    0.000000    0.000000
        0.000000    5.900000    0.000000
        0.000000    0.000000   20.000000

siesta: Stress tensor Voigt (kbar):     ${kbar}.00      0.00      0.00      0.00      0.00      0.00
EOF
    # Real structure.fdf, adapted from test/4-workflow/1-strain/prep/structure.fdf
    # -- a real 10-atom sheet with a genuine vacuum gap along z only.
    cat > "$folder/structure.fdf" <<'EOF'
NumberOfSpecies    1
NumberofAtoms      10

%block ChemicalSpeciesLabel
 1   6   C
%endblock ChemicalSpeciesLabel

LatticeConstant 1.0 Ang

AtomicCoordinatesFormat  Fractional

%block LatticeVectors
 5.9079999924   0.0000000000   0.0000000000
 -2.9539999962   5.1164780790   0.0000000000
 0.0000000000   0.0000000000   20.0000000000
%endblock LatticeVectors

%block AtomicCoordinatesAndAtomicSpecies
  0.776380002   0.288109988   0.500000000   1
  0.711880028   0.488260001   0.500000000   1
  0.511740029   0.223619998   0.500000000   1
  0.288109988   0.776380002   0.500000000   1
  0.488260001   0.711880028   0.500000000   1
  0.223619998   0.511739969   0.500000000   1
  0.000000000   0.000000000   0.500000000   1
  0.000000000   0.251800001   0.500000000   1
  0.748190045   0.748189986   0.500000000   1
  0.251809984   0.000000000   0.500000000   1
%endblock AtomicCoordinatesAndAtomicSpecies
EOF
done
stb-strainAnalysis --file calc.out --dir auto2d_run --no-intro > log_auto2d.txt 2>&1
check_exit_code $? 0
check_contains "Auto-detected from a real structure file" log_auto2d.txt
check_contains "Dimensionality    : 2D (auto-detected)" log_auto2d.txt
check_contains "N/m" log_auto2d.txt

echo "Testing: 1D (vacuum along x and y, wire along z) -- no --dimensionality flag given"
mkdir -p auto1d_run/z/strain_z_0.00 auto1d_run/z/strain_z_1.00 auto1d_run/z/strain_z_2.00
for pair in "0.00:0" "1.00:50" "2.00:100"; do
    val="${pair%%:*}"
    kbar="${pair##*:}"
    folder="auto1d_run/z/strain_z_$val"
    cat > "$folder/calc.out" <<EOF
outcell: Unit cell vectors (Ang):
       30.000000    0.000000    0.000000
        0.000000   30.000000    0.000000
        0.000000    0.000000   10.000000

siesta: Stress tensor Voigt (kbar):     0.00      0.00      ${kbar}.00      0.00      0.00      0.00
EOF
    # Minimal 2-atom wire along z, centered in the (vacuum-padded) x/y
    # cross-section -- same 30x30x10 Ang cell as the existing --dimensionality
    # 1d fixture below, but with a genuine structure file this time.
    cat > "$folder/wire.fdf" <<'EOF'
NumberOfSpecies    1
NumberofAtoms      2

%block ChemicalSpeciesLabel
 1  14   Si
%endblock ChemicalSpeciesLabel

LatticeConstant 1.0 Ang

AtomicCoordinatesFormat  Fractional

%block LatticeVectors
 30.000000000000    0.000000000000    0.000000000000
  0.000000000000   30.000000000000    0.000000000000
  0.000000000000    0.000000000000   10.000000000000
%endblock LatticeVectors

%block AtomicCoordinatesAndAtomicSpecies
  0.500000000   0.500000000   0.000000000   1
  0.500000000   0.500000000   0.500000000   1
%endblock AtomicCoordinatesAndAtomicSpecies
EOF
done
stb-strainAnalysis --file calc.out --dir auto1d_run --no-intro > log_auto1d.txt 2>&1
check_exit_code $? 0
check_contains "Auto-detected from a real structure file" log_auto1d.txt
check_contains "Dimensionality    : 1D (auto-detected)" log_auto1d.txt
check_contains "Peak Force" log_auto1d.txt


# --- 5. 1D mode (--dimensionality 1d / --cross-section), synthetic
#     nanotube-like fixture: large in-plane vacuum (30x30 Ang), single
#     periodic axis (z). No .fdf here -- exercises the explicit-override
#     path (auto-detection is exercised separately in step 4c above). ---
echo -e "\n--- Testing --dimensionality 1d (Force in nN, auto cross-section) ---"
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

stb-strainAnalysis --file calc.out --dir . --dimensionality 1d --save-gnuplot --no-intro > log_1d.txt 2>&1
check_exit_code $? 0
# Cell cross-section = |30x0x0 . 0x30x0| = 900 Ang^2; 100 kBar * 900 * 1e-3 = 90 nN.
check_contains "Peak Force" log_1d.txt
check_contains "90.0000 nN" log_1d.txt

echo "Testing: --dimensionality 1d --cross-section (matches cell area -> recovers plain 10 GPa)"
stb-strainAnalysis --file calc.out --dir . --dimensionality 1d --cross-section 900 --no-intro > log_1d_xsec.txt 2>&1
check_contains "Peak Stress (conventional, cross-section=900.0000" log_1d_xsec.txt
check_contains "10.0000 GPa" log_1d_xsec.txt

echo "Testing: --save-gnuplot output for a 1D run"
check_success z_curve.dat
check_success z_curve.gplot
if command -v gnuplot > /dev/null 2>&1; then
    gnuplot z_curve.gplot > log_gnuplot_1d.txt 2>&1
    check_exit_code $? 0
    check_success z_curve.pdf
else
    echo -e "   -> ${YELLOW}Skipped:${NC} gnuplot not installed, skipping actual render"
fi

echo "Testing: invalid --dimensionality value is rejected by argparse"
stb-strainAnalysis --file calc.out --dir . --dimensionality bogus --no-intro > log_bad_dim.txt 2>&1
check_exit_code $? 2
check_contains "invalid choice" log_bad_dim.txt

echo "Testing: --cross-section without --dimensionality 1d/2d is rejected immediately (static check)"
stb-strainAnalysis --file calc.out --dir . --dimensionality 3d --cross-section 100 --no-intro > log_xsec_static.txt 2>&1
check_exit_code $? 2
check_contains "only applies together with --dimensionality 1d" log_xsec_static.txt

echo "Testing: --cross-section with --dimensionality auto resolving to 3D is rejected after"
echo "         resolution (deferred check -- this fixture has no .fdf, so auto resolves 3D)"
stb-strainAnalysis --file calc.out --dir . --dimensionality auto --cross-section 100 --no-intro > log_xsec_deferred.txt 2>&1
check_exit_code $? 1
check_contains "only applies when the run is 1D" log_xsec_deferred.txt
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
if grep -q "comparing automatically" log_canon.txt 2>/dev/null; then
    echo -e "   -> ${RED}Failed:${NC} 'xx' and 'x' were treated as 2 different directions"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} 'xx' and 'x' folders were grouped as the same direction"
    PASS=$((PASS+1))
fi

echo "Testing: --version"
stb-strainAnalysis --version > log_version.txt 2>&1
check_contains "stb-strainAnalysis" log_version.txt

echo "Testing: --help documents --dimensionality/--cross-section/--yield/--dir/--save-report/--save-gnuplot/--view"
stb-strainAnalysis --help > log_help.txt 2>&1
check_contains "\-\-dimensionality" log_help.txt
check_contains "auto,3d,2d,1d" log_help.txt
check_contains "cross-section" log_help.txt
check_contains "yield" log_help.txt
check_contains "strain_runs" log_help.txt
check_contains "\-\-save-report" log_help.txt
check_contains "\-\-save-gnuplot" log_help.txt
check_contains "\-\-view" log_help.txt


# --- 7. Interactive path (stb-suite, shortcut 4.1.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.1.2) ---"

echo "Testing: navigate 4.1.2 -> dir '.' -> calc.out -> dimensionality auto -> skip advanced ->"
echo "         no yield -> default output-dir -> no report -> save gnuplot -> no view -> quit"
rm -f x_curve.dat x_curve.gplot stb_strainAnalysis_report.txt
# Prompts in order: dir ('.'), file ('calc.out'), dimensionality (blank ->
# auto), advanced settings (n -> skip), yield (n), output-dir (blank ->
# '.'), save-report (n), save-gnuplot (y -- so a real file is produced for
# check_success below), view (n), "Press Enter to continue", quit.
printf '4.1.2\n.\ncalc.out\n\nn\nn\n\nn\ny\nn\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_contains "Initial Slope" log_menu.txt
check_success x_curve.dat
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
