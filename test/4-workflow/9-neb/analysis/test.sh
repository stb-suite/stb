#!/bin/bash

# --- Setup ---
# Smoke test for stb-nebAnalysis (NEB Analysis, item 4.9.2). Uses stb-neb
# itself (real tool, no SIESTA needed for prep) to build a 5-image band
# against the same fixture as ../prep/ (linear interpolation, no ML --
# cheap and deterministic), then fabricates "siesta: FreeEng"/SCF/Max-force
# lines per image (synthetic, same style already used by
# 8-adsorption/analysis/test.sh) to exercise the analysis side without
# needing a real SIESTA run.
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREP_DIR="$(cd "$FIXTURE_DIR/../prep" && pwd)"
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
echo "--- Starting tester for STB-NebAnalysis (item 4.9.2) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$PREP_DIR/initial.fdf" "$TEST_DIR/"
cp "$PREP_DIR/final.fdf" "$TEST_DIR/"
cp "$PREP_DIR/calc.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Build a real 5-image band via stb-neb, then fabricate FreeEng/SCF/Max ---
echo -e "\n--- Testing analysis of a 5-image NEB band ---"
stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 5 --no-intro > log_prep.txt 2>&1
check_exit_code $? 0
n_images=$(find . -maxdepth 1 -type d -name 'image_*' | wc -l)
if [ "$n_images" -eq 5 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} prep produced 5 image folders"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} expected 5 image folders, found $n_images"
    FAIL=$((FAIL+1))
fi

# Asymmetric profile so forward != backward barrier is actually exercised;
# image_02 is deliberately unconverged / high residual force (verifies
# check_scf_and_force, freshly moved into core/siesta_log.py this session,
# still works correctly for a second consumer).
printf 'siesta: FreeEng =    -300.000000\nSCF cycle converged after 12 iterations\nsiesta: Atomic forces (eV/Ang):\n   Max    0.010000\n' > image_00/calc.out
printf 'siesta: FreeEng =    -299.500000\nSCF cycle converged after 12 iterations\nsiesta: Atomic forces (eV/Ang):\n   Max    0.010000\n' > image_01/calc.out
printf 'siesta: FreeEng =    -298.700000\nsiesta: Atomic forces (eV/Ang):\n   Max    0.350000\n' > image_02/calc.out
printf 'siesta: FreeEng =    -299.100000\nSCF cycle converged after 12 iterations\nsiesta: Atomic forces (eV/Ang):\n   Max    0.010000\n' > image_03/calc.out
printf 'siesta: FreeEng =    -299.900000\nSCF cycle converged after 12 iterations\nsiesta: Atomic forces (eV/Ang):\n   Max    0.010000\n' > image_04/calc.out

stb-nebAnalysis --dir . --no-intro > log_analysis.txt 2>&1
check_exit_code $? 0
check_contains "Highest-energy image (approx. TS) : image_02" log_analysis.txt
check_contains "Forward barrier  (TS - initial)   : 1.300000 eV" log_analysis.txt
check_contains "Backward barrier (TS - final)     : 1.200000 eV" log_analysis.txt
check_contains "Reaction energy  (final - initial) : 0.100000 eV, endothermic" log_analysis.txt
check_contains "Spline-fitted barrier (smoothed)" log_analysis.txt
check_contains "Approximate, interpolated-path estimate" log_analysis.txt
check_contains "never confirmed SCF convergence.*image_02" log_analysis.txt
check_contains "residual force above --force-tolerance.*image_02" log_analysis.txt
check_success neb_curve.dat
check_success neb_curve.gplot
check_success neb_report.txt
check_contains "\[0\] RUN METADATA" neb_report.txt
check_contains "\[1\] IMAGE ENERGIES" neb_report.txt
check_contains "\[2\] BARRIER ANALYSIS" neb_report.txt
check_contains "\[3\] SUMMARY" neb_report.txt
check_contains "0.615000  -299.500000" neb_curve.dat


# --- 2a2. Too few images (3, the minimum) for a spline fit -- NOTE, not an error ---
echo -e "\n--- Testing spline-fit NOTE with only 3 images ---"
mkdir -p tiny
cp initial.fdf final.fdf calc.fdf tiny/
(
    cd tiny
    stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 3 --no-intro > log_tiny_prep.txt 2>&1
    printf 'siesta: FreeEng =    -300.000000\nSCF cycle converged after 12 iterations\nsiesta: Atomic forces (eV/Ang):\n   Max    0.010000\n' > image_00/calc.out
    printf 'siesta: FreeEng =    -298.900000\nSCF cycle converged after 12 iterations\nsiesta: Atomic forces (eV/Ang):\n   Max    0.010000\n' > image_01/calc.out
    printf 'siesta: FreeEng =    -299.900000\nSCF cycle converged after 12 iterations\nsiesta: Atomic forces (eV/Ang):\n   Max    0.010000\n' > image_02/calc.out
    stb-nebAnalysis --dir . --no-intro > log_tiny_analysis.txt 2>&1
)
check_contains "Not enough images .need >= 4. for a spline-fitted barrier" tiny/log_tiny_analysis.txt


# --- 2b. --ml-neb-sourced band gets the stronger caveat wording ---
echo -e "\n--- Testing --ml-neb caveat wording ---"
mkdir -p mlneb
cp initial.fdf final.fdf calc.fdf mlneb/
(
    cd mlneb
    stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 5 --ml-neb --ml-max-steps 30 --no-intro \
        > log_mlneb_prep.txt 2>&1
    printf 'siesta: FreeEng =    -300.000000\nSCF cycle converged after 12 iterations\nsiesta: Atomic forces (eV/Ang):\n   Max    0.010000\n' > image_00/calc.out
    printf 'siesta: FreeEng =    -299.600000\nSCF cycle converged after 12 iterations\nsiesta: Atomic forces (eV/Ang):\n   Max    0.010000\n' > image_01/calc.out
    printf 'siesta: FreeEng =    -298.900000\nSCF cycle converged after 12 iterations\nsiesta: Atomic forces (eV/Ang):\n   Max    0.010000\n' > image_02/calc.out
    printf 'siesta: FreeEng =    -299.200000\nSCF cycle converged after 12 iterations\nsiesta: Atomic forces (eV/Ang):\n   Max    0.010000\n' > image_03/calc.out
    printf 'siesta: FreeEng =    -299.900000\nSCF cycle converged after 12 iterations\nsiesta: Atomic forces (eV/Ang):\n   Max    0.010000\n' > image_04/calc.out
    stb-nebAnalysis --dir . --no-intro > log_mlneb_analysis.txt 2>&1
)
check_contains "already converged a real climbing-image band on MACE-MP-0" mlneb/log_mlneb_analysis.txt


# --- 2c. --apply copies the highest-energy image's structure.fdf ---
echo -e "\n--- Testing --apply ---"
stb-nebAnalysis --dir . --apply ts_guess.fdf --no-intro > log_apply.txt 2>&1
check_contains "Applied.*image_02" log_apply.txt
check_success ts_guess.fdf
check_contains "NumberofAtoms      9" ts_guess.fdf


# --- 3. A folder missing calc.out is skipped, not fatal ---
echo -e "\n--- Testing that an image missing calc.out is skipped ---"
mv image_03/calc.out image_03/calc.out.bak
stb-nebAnalysis --dir . --no-intro > log_partial.txt 2>&1
check_contains "SKIP" log_partial.txt
check_contains "skipped: 1" log_partial.txt
mv image_03/calc.out.bak image_03/calc.out


# --- 4. neb_setup.txt missing entirely (fallback path) ---
echo -e "\n--- Testing fallback when neb_setup.txt is missing ---"
mv neb_setup.txt neb_setup.txt.bak
stb-nebAnalysis --dir . --no-intro > log_fallback.txt 2>&1
check_exit_code $? 0
check_contains "No 'neb_setup.txt' image table found" log_fallback.txt
check_contains "Highest-energy image (approx. TS) : image_02" log_fallback.txt
mv neb_setup.txt.bak neb_setup.txt


# --- 5. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: missing directory entirely"
stb-nebAnalysis --dir does_not_exist --no-intro > log_no_dir.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_no_dir.txt

echo "Testing: directory with no image_* folders"
mkdir -p empty_dir
stb-nebAnalysis --dir empty_dir --no-intro > log_no_images.txt 2>&1
check_exit_code $? 1
check_contains "Did you run stb-neb" log_no_images.txt

echo "Testing: missing initial/final endpoint energy"
mkdir -p missing_endpoint/image_00 missing_endpoint/image_01
echo "siesta: FreeEng =    -1.0" > missing_endpoint/image_01/calc.out
stb-nebAnalysis --dir missing_endpoint --no-intro > log_missing_endpoint.txt 2>&1
check_exit_code $? 1
check_contains "Could not read the initial and/or final endpoint" log_missing_endpoint.txt

echo "Testing: --force-tolerance"
stb-nebAnalysis --dir . --force-tolerance 1.0 --no-intro > log_force_tol.txt 2>&1
if grep -q "residual force above --force-tolerance.*image_02" log_force_tol.txt; then
    echo -e "   -> ${RED}Failed:${NC} unexpected force warning with a loose --force-tolerance"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no force warning with --force-tolerance 1.0 (image_02's 0.35 eV/Ang is under it)"
    PASS=$((PASS+1))
fi

echo "Testing: --version"
stb-nebAnalysis --version > log_version.txt 2>&1
check_contains "stb-nebAnalysis" log_version.txt

echo "Testing: --help documents --dir/--file/--apply"
stb-nebAnalysis --help > log_help.txt 2>&1
check_contains "dir" log_help.txt
check_contains "file" log_help.txt
check_contains "apply" log_help.txt


# --- 6. Interactive path (stb-suite, shortcut 4.9.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.9.2) ---"

echo "Testing: navigate 4.9.2 -> defaults -> quit"
rm -f neb_curve.dat neb_report.txt
# 4.9.2 (menu code) / . (dir) / "" (out_file default) / "" (force-tolerance
# default) / "" (apply_target: skip) / "" (Press Enter to continue) / 0 (quit)
printf '4.9.2\n.\n\n\n\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_contains "Highest-energy image (approx. TS) : image_02" log_menu.txt
check_success neb_curve.dat


popd > /dev/null

# --- 7. Summary ---
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
