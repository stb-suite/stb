#!/bin/bash

# --- Setup ---
# Smoke test for stb-nebAnalysis (NEB Analysis, item 4.9.2). Uses stb-neb
# itself (real tool) to build bands against the same fixture as ../prep/,
# then fabricates "siesta: FreeEng"/SCF/Max-force lines per image
# (synthetic, same style already used by 8-adsorption/analysis/test.sh) to
# exercise the analysis side without needing a real SIESTA run. Now covers
# all 4 modes stb-neb can produce (see CLAUDE.md's mode table): the main
# block below uses --mode 4 (no MACE, fastest, cycle_00/ nesting) for the
# core aggregation-logic tests; separate blocks cover mode 2 (MACE +
# single-point, image_NN/ directly), mode 1 (JSON only), and mode 3's
# multi-cycle auto-discovery + convergence-history plot.
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


# --- 2. Build a real 5-image band via stb-neb --mode 4 (no MACE, cycle_00/), then fabricate FreeEng/SCF/Max ---
echo -e "\n--- Testing analysis of a 5-image NEB band (--mode 4) ---"
stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 5 --mode 4 --no-intro > log_prep.txt 2>&1
check_exit_code $? 0
n_images=$(find cycle_00 -maxdepth 1 -type d -name 'image_*' | wc -l)
if [ "$n_images" -eq 5 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} prep produced 5 image folders"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} expected 5 image folders, found $n_images"
    FAIL=$((FAIL+1))
fi

# Asymmetric profile so forward != backward barrier is actually exercised;
# image_02 is deliberately unconverged / high residual force (verifies
# check_scf_and_force still works correctly for this consumer too).
printf 'siesta: FreeEng =    -300.000000\nSCF cycle converged after 12 iterations\nsiesta: Atomic forces (eV/Ang):\n   Max    0.010000\n' > cycle_00/image_00/calc.out
printf 'siesta: FreeEng =    -299.500000\nSCF cycle converged after 12 iterations\nsiesta: Atomic forces (eV/Ang):\n   Max    0.010000\n' > cycle_00/image_01/calc.out
printf 'siesta: FreeEng =    -298.700000\nsiesta: Atomic forces (eV/Ang):\n   Max    0.350000\n' > cycle_00/image_02/calc.out
printf 'siesta: FreeEng =    -299.100000\nSCF cycle converged after 12 iterations\nsiesta: Atomic forces (eV/Ang):\n   Max    0.010000\n' > cycle_00/image_03/calc.out
printf 'siesta: FreeEng =    -299.900000\nSCF cycle converged after 12 iterations\nsiesta: Atomic forces (eV/Ang):\n   Max    0.010000\n' > cycle_00/image_04/calc.out

stb-nebAnalysis --dir . --save-report --no-intro > log_analysis.txt 2>&1
check_exit_code $? 0
check_contains "\[NOTE\] Analyzing 'cycle_00' -- NOT YET CONVERGED" log_analysis.txt
check_contains "Highest-energy image (approx. TS) : image_02" log_analysis.txt
check_contains "Forward barrier  (TS - initial)   : 1.300000 eV" log_analysis.txt
check_contains "Backward barrier (TS - final)     : 1.200000 eV" log_analysis.txt
check_contains "Reaction energy  (final - initial) : 0.100000 eV, endothermic" log_analysis.txt
check_contains "Spline-fitted barrier (smoothed)" log_analysis.txt
check_contains "gone through real-DFT NEB refinement" log_analysis.txt
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

echo "Testing: only 1 complete cycle so far -- no convergence-history plot yet"
check_contains "Not enough complete cycles yet" neb_report.txt


# --- 2a2. Too few images (3, the minimum) for a spline fit -- NOTE, not an error ---
echo -e "\n--- Testing spline-fit NOTE with only 3 images ---"
mkdir -p tiny
cp initial.fdf final.fdf calc.fdf tiny/
(
    cd tiny
    stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 3 --mode 4 --no-intro > log_tiny_prep.txt 2>&1
    printf 'siesta: FreeEng =    -300.000000\nSCF cycle converged after 12 iterations\nsiesta: Atomic forces (eV/Ang):\n   Max    0.010000\n' > cycle_00/image_00/calc.out
    printf 'siesta: FreeEng =    -298.900000\nSCF cycle converged after 12 iterations\nsiesta: Atomic forces (eV/Ang):\n   Max    0.010000\n' > cycle_00/image_01/calc.out
    printf 'siesta: FreeEng =    -299.900000\nSCF cycle converged after 12 iterations\nsiesta: Atomic forces (eV/Ang):\n   Max    0.010000\n' > cycle_00/image_02/calc.out
    stb-nebAnalysis --dir . --no-intro > log_tiny_analysis.txt 2>&1
)
check_contains "Not enough images .need >= 4. for a spline-fitted barrier" tiny/log_tiny_analysis.txt


# --- 2b. Mode 2 (MACE-MP-0) -sourced band gets the stronger caveat wording ---
echo -e "\n--- Testing --mode 2 caveat wording ---"
mkdir -p mode2
cp initial.fdf final.fdf calc.fdf mode2/
(
    cd mode2
    stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 5 --mode 2 --ml-max-steps 30 --no-intro \
        > log_mode2_prep.txt 2>&1
    printf 'siesta: FreeEng =    -300.000000\nSCF cycle converged after 12 iterations\nsiesta: Atomic forces (eV/Ang):\n   Max    0.010000\n' > image_00/calc.out
    printf 'siesta: FreeEng =    -299.600000\nSCF cycle converged after 12 iterations\nsiesta: Atomic forces (eV/Ang):\n   Max    0.010000\n' > image_01/calc.out
    printf 'siesta: FreeEng =    -298.900000\nSCF cycle converged after 12 iterations\nsiesta: Atomic forces (eV/Ang):\n   Max    0.010000\n' > image_02/calc.out
    printf 'siesta: FreeEng =    -299.200000\nSCF cycle converged after 12 iterations\nsiesta: Atomic forces (eV/Ang):\n   Max    0.010000\n' > image_03/calc.out
    printf 'siesta: FreeEng =    -299.900000\nSCF cycle converged after 12 iterations\nsiesta: Atomic forces (eV/Ang):\n   Max    0.010000\n' > image_04/calc.out
    stb-nebAnalysis --dir . --no-intro > log_mode2_analysis.txt 2>&1
)
check_contains "already converged a real climbing-image band on MACE-MP-0" mode2/log_mode2_analysis.txt
echo "Testing: mode 2's images live directly under the root (no cycle_NN/ nesting)"
check_success mode2/image_02/calc.out


# --- 2c. Mode 1 (100% MACE, JSON only) analysis ---
echo -e "\n--- Testing --mode 1 analysis (JSON only, no calc.out anywhere) ---"
mkdir -p mode1
cp initial.fdf final.fdf calc.fdf mode1/
(
    cd mode1
    stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 5 --mode 1 --ml-max-steps 30 --no-intro \
        > log_mode1_prep.txt 2>&1
    stb-nebAnalysis --dir . --save-report --no-intro > log_mode1_analysis.txt 2>&1
)
check_exit_code $? 0
check_contains "MODE 1, MACE-MP-0" mode1/log_mode1_analysis.txt
check_contains "Mode 1 never touches SIESTA" mode1/log_mode1_analysis.txt
check_contains "MACE-MP-0 barrier (forward)" mode1/log_mode1_analysis.txt
check_success mode1/neb_curve.dat
check_success mode1/neb_report.txt
echo "Testing: --apply is refused (cleanly) for mode 1 -- no per-image structure.fdf exists"
(cd mode1 && stb-nebAnalysis --dir . --apply ts.fdf --no-intro > log_mode1_apply.txt 2>&1)
check_contains "Mode 1 has no per-image structure.fdf" mode1/log_mode1_apply.txt
if [ -f mode1/ts.fdf ]; then
    echo -e "   -> ${RED}Failed:${NC} ts.fdf unexpectedly written for mode 1"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no structure.fdf written for mode 1 --apply"
    PASS=$((PASS+1))
fi


# --- 2d. Mode 3: multi-cycle auto-discovery + convergence-history plot ---
echo -e "\n--- Testing --mode 3 multi-cycle analysis ---"
mkdir -p mode3
cp initial.fdf final.fdf calc.fdf mode3/
(
    cd mode3
    stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 5 --mode 3 --ml-max-steps 30 --no-intro \
        > log_mode3_prep.txt 2>&1
    printf 'siesta: FreeEng =    -300.000000\nSCF cycle converged after 12 iterations\n' > cycle_00/image_00/calc.out
    printf 'siesta: FreeEng =    -299.600000\nSCF cycle converged after 12 iterations\n' > cycle_00/image_01/calc.out
    printf 'siesta: FreeEng =    -298.950000\nSCF cycle converged after 12 iterations\n' > cycle_00/image_02/calc.out
    printf 'siesta: FreeEng =    -299.300000\nSCF cycle converged after 12 iterations\n' > cycle_00/image_03/calc.out
    printf 'siesta: FreeEng =    -299.900000\nSCF cycle converged after 12 iterations\n' > cycle_00/image_04/calc.out
)
echo "Testing: only cycle_00 exists -> analyzed, correctly flagged NOT YET CONVERGED"
(cd mode3 && stb-nebAnalysis --dir . --no-intro > log_mode3_analysis1.txt 2>&1)
check_contains "\[NOTE\] Analyzing 'cycle_00' -- NOT YET CONVERGED" mode3/log_mode3_analysis1.txt

echo "Testing: a second cycle (cycle_01, a refined re-evaluation) makes the analysis pick it instead"
cp -r mode3/cycle_00 mode3/cycle_01
printf 'siesta: FreeEng =    -300.000000\nSCF cycle converged after 12 iterations\n' > mode3/cycle_01/image_00/calc.out
printf 'siesta: FreeEng =    -299.650000\nSCF cycle converged after 12 iterations\n' > mode3/cycle_01/image_01/calc.out
printf 'siesta: FreeEng =    -298.900000\nSCF cycle converged after 12 iterations\n' > mode3/cycle_01/image_02/calc.out
printf 'siesta: FreeEng =    -299.350000\nSCF cycle converged after 12 iterations\n' > mode3/cycle_01/image_03/calc.out
printf 'siesta: FreeEng =    -299.900000\nSCF cycle converged after 12 iterations\n' > mode3/cycle_01/image_04/calc.out
(cd mode3 && stb-nebAnalysis --dir . --save-report --no-intro > log_mode3_analysis2.txt 2>&1)
check_contains "\[NOTE\] Analyzing 'cycle_01' -- NOT YET CONVERGED" mode3/log_mode3_analysis2.txt
check_contains "Highest-energy image (approx. TS) : image_02  (E = -298.900000 eV)" mode3/log_mode3_analysis2.txt

echo "Testing: barrier-vs-cycle convergence-history plot, 2 complete cycles"
check_contains "\[2b\] BARRIER VS. CYCLE" mode3/neb_report.txt
check_success mode3/neb_cycle_convergence.png
check_contains "cycle_00: forward barrier" mode3/neb_report.txt
check_contains "cycle_01: forward barrier" mode3/neb_report.txt

echo "Testing: NEB_CONVERGED sentinel changes the [NOTE] wording"
touch mode3/NEB_CONVERGED
(cd mode3 && stb-nebAnalysis --dir . --no-intro > log_mode3_analysis3.txt 2>&1)
check_contains "\[NOTE\] Analyzing 'cycle_01' -- NEB_CONVERGED." mode3/log_mode3_analysis3.txt

echo "Testing: --apply on a mode-3 cycle copies from the auto-discovered cycle_01, not cycle_00"
(cd mode3 && stb-nebAnalysis --dir . --apply ts_guess.fdf --no-intro > log_mode3_apply.txt 2>&1)
check_contains "Applied.*image_02" mode3/log_mode3_apply.txt
check_success mode3/ts_guess.fdf
python3 -c "
import sys
from stb.core import structure_io
import numpy as np
applied = structure_io.to_pymatgen(structure_io.read_fdf('mode3/ts_guess.fdf'))
cycle1 = structure_io.to_pymatgen(structure_io.read_fdf('mode3/cycle_01/image_02/structure.fdf'))
sys.exit(0 if np.allclose(applied.cart_coords, cycle1.cart_coords) else 1)
"
check_exit_code $? 0


# --- 2e. --apply copies the highest-energy image's structure.fdf (mode 4, root case) ---
echo -e "\n--- Testing --apply (mode 4) ---"
stb-nebAnalysis --dir . --apply ts_guess.fdf --no-intro > log_apply.txt 2>&1
check_contains "Applied.*image_02" log_apply.txt
check_success ts_guess.fdf
check_contains "NumberofAtoms      9" ts_guess.fdf


# --- 3. A folder missing calc.out is skipped, not fatal ---
echo -e "\n--- Testing that an image missing calc.out is skipped ---"
mv cycle_00/image_03/calc.out cycle_00/image_03/calc.out.bak
stb-nebAnalysis --dir . --no-intro > log_partial.txt 2>&1
check_contains "SKIP" log_partial.txt
check_contains "skipped: 1" log_partial.txt
mv cycle_00/image_03/calc.out.bak cycle_00/image_03/calc.out


# --- 4. neb_setup.txt missing entirely (fallback path, mode defaults to 2) ---
echo -e "\n--- Testing fallback when neb_setup.txt is missing ---"
mkdir -p no_setup
cp -r cycle_00/image_00 cycle_00/image_01 cycle_00/image_02 cycle_00/image_03 cycle_00/image_04 no_setup/
(cd no_setup && stb-nebAnalysis --dir . --no-intro > log_fallback.txt 2>&1)
check_exit_code $? 0
check_contains "No 'neb_setup.txt' image table found" no_setup/log_fallback.txt
check_contains "Highest-energy image (approx. TS) : image_02" no_setup/log_fallback.txt


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

echo "Testing: mode 3/4 directory with no cycle_NN folders at all"
mkdir -p empty_cycle_dir
echo "# MODE: 3" > empty_cycle_dir/neb_setup.txt
stb-nebAnalysis --dir empty_cycle_dir --no-intro > log_no_cycles.txt 2>&1
check_exit_code $? 1
check_contains "No 'cycle_NN' folder found" log_no_cycles.txt

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
