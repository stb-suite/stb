#!/bin/bash

# --- Setup ---
# Smoke test for stb-herRefs (HER Stage 2: References & ZPE Prep, item
# 4.13.2). Chains a real stb-her (Stage 1, tested separately) run on the
# graphene fixture, then fabricates a relaxed calc.out (known FreeEng +
# a deliberately-shifted "outcoor: Relaxed atomic coordinates" block, to
# confirm Stage 2 uses the RELAXED geometry, not Stage 1's pre-relaxation
# input) for the single candidate site.
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
    if grep -q -- "$1" "$2" 2>/dev/null; then
        echo -e "   -> ${GREEN}Verified:${NC} '$1' found in '$2'"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} '$1' NOT found in '$2'"
        FAIL=$((FAIL+1))
    fi
}

check_not_contains() {
    if grep -q -- "$1" "$2" 2>/dev/null; then
        echo -e "   -> ${RED}Failed:${NC} '$1' unexpectedly found in '$2'"
        FAIL=$((FAIL+1))
    else
        echo -e "   -> ${GREEN}Verified:${NC} '$1' NOT found in '$2' (as expected)"
        PASS=$((PASS+1))
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

# Runs stb-her (Stage 1, tested separately) on the graphene fixture with
# a single ontop site, then fabricates a relaxed calc.out for it: known
# FreeEng, and a relaxed-coordinates block where H has moved slightly
# from its Stage-1 input position -- Stage 2 must use THIS relaxed
# position, not the original.
make_site_disp() {
    rm -rf her_study
    stb-her -s structure.fdf -c calc.fdf -p . --site-type ontop --height 1.2 --no-intro \
        > /dev/null 2>&1
    python3 - <<'PYEOF'
with open("her_study/sites/site_1_ontop/calc.out", "w") as f:
    f.write("siesta: SCF Convergence by DM criterion\n")
    f.write("SCF cycle converged after 15 iterations\n")
    f.write("siesta: FreeEng =        -499.950000\n")
    f.write("outcoor: Relaxed atomic coordinates (fractional)\n")
    f.write("    0.00000000    0.00000000    0.50000000   1\n")
    f.write("    0.33333333    0.66666667    0.50000000   1\n")
    f.write("    0.00100000    0.00100000    0.58000000   2\n")
    f.write("siesta: Atomic forces (eV/Ang):\n")
    f.write("     1    0.001    0.001    0.001\n")
    f.write("     2    0.001    0.001    0.001\n")
    f.write("     3    0.001    0.001    0.001\n")
    f.write("     Max    0.001\n")
PYEOF
}


# --- 1. Preparation ---
echo "--- Starting tester for STB-HERRefs stage 2: references & ZPE prep (item 4.13.2) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$PREP_DIR/structure.fdf" "$TEST_DIR/"
cp "$PREP_DIR/calc.fdf" "$TEST_DIR/"
echo "# placeholder pseudopotential" > "$TEST_DIR/C.psf"
echo "# placeholder pseudopotential" > "$TEST_DIR/H.psf"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Missing Stage-1 output guard ---
echo -e "\n--- Testing the missing-Stage-1-output guard ---"
rm -rf her_study
stb-herRefs --directory her_study -p . --no-intro > log_no_stage1.txt 2>&1
check_exit_code $? 1
check_contains "run stb-her" log_no_stage1.txt


# --- 3. Default run (--zpe-mode local) ---
echo -e "\n--- Testing default generation (--zpe-mode local) ---"
make_site_disp
stb-herRefs --directory her_study -p . --zpe-mode local --no-intro > log_local.txt 2>&1
check_exit_code $? 0
check_contains "Winning site    : site_1_ontop" log_local.txt
check_success her_study/00_clean_slab/structure.fdf
check_success her_study/02_h2_molecule/structure.fdf
check_success her_study/03_slab_deformed/structure.fdf
check_success her_study/04_slab_ghost/structure.fdf
check_success her_study/06_h_ghost_slab/structure.fdf
check_success her_study/07_h_isolated/structure.fdf
check_success her_study/05_zpe_calc/disp_001/structure.fdf
check_success her_study/05_zpe_calc/disp_006/structure.fdf
check_success her_study/05_zpe_calc/zpe_local_meta.json
check_success her_study/her_stage2.txt

echo "Testing: 03_slab_deformed has no H (2 atoms, C only)"
check_contains "NumberofAtoms      2" her_study/03_slab_deformed/structure.fdf
check_not_contains "   1   H" her_study/03_slab_deformed/structure.fdf

echo "Testing: 04_slab_ghost has H_ghost (negative Z), 06_h_ghost_slab has C_ghost"
check_contains "H_ghost" her_study/04_slab_ghost/structure.fdf
check_contains "C_ghost" her_study/06_h_ghost_slab/structure.fdf

echo "Testing: 07_h_isolated is a single H atom at the RELAXED z=0.58 (not Stage 1's input z=0.56)"
check_contains "NumberofAtoms      1" her_study/07_h_isolated/structure.fdf
check_contains "0.58000000" her_study/07_h_isolated/structure.fdf

echo "Testing: H2 molecule relaxes (not single-point), Gamma-only + spin-polarized"
check_contains "MD.TypeOfRun          CG" her_study/02_h2_molecule/calc.fdf
check_not_contains "MD.NumCGsteps         0" her_study/02_h2_molecule/calc.fdf
check_contains "kgrid.MonkhorstPack   \[1  1  1\]" her_study/02_h2_molecule/calc.fdf
check_contains "Spin.*polarized" her_study/02_h2_molecule/calc.fdf

echo "Testing: every other reference folder is forced single-point"
check_contains "MD.NumCGsteps         0" her_study/00_clean_slab/calc.fdf
check_contains "MD.NumCGsteps         0" her_study/03_slab_deformed/calc.fdf
check_contains "MD.NumCGsteps         0" her_study/04_slab_ghost/calc.fdf


# --- 4. --zpe-mode standard (nothing generated) ---
echo -e "\n--- Testing --zpe-mode standard (no ZPE folders) ---"
make_site_disp
stb-herRefs --directory her_study -p . --zpe-mode standard --no-intro > log_standard.txt 2>&1
check_exit_code $? 0
check_contains "Standard mode -- no folders generated" log_standard.txt
python3 -c "
import os
assert not os.path.isdir('her_study/05_zpe_calc'), '05_zpe_calc should not exist in standard mode'
print('OK')
" > log_standard_check.txt 2>&1
check_contains "OK" log_standard_check.txt


# --- 5. --zpe-mode full (two phonon displacement sets) ---
echo -e "\n--- Testing --zpe-mode full (site + clean slab phonon displacement sets) ---"
make_site_disp
stb-herRefs --directory her_study -p . --zpe-mode full --supercell 1 1 1 --no-intro \
    > log_full.txt 2>&1
check_exit_code $? 0
check_contains "needs a full phonon calculation of BOTH" log_full.txt
check_success her_study/05_zpe_calc_site/disp-001/structure.fdf
check_success her_study/05_zpe_calc_site/phonopy_disp.yaml
check_success her_study/05_zpe_calc_clean/disp-001/structure.fdf
check_success her_study/05_zpe_calc_clean/phonopy_disp.yaml
check_contains "SystemLabel her_zpe_site" her_study/05_zpe_calc_site/disp-001/calc.fdf
check_contains "SystemLabel her_zpe_clean" her_study/05_zpe_calc_clean/disp-001/calc.fdf


# --- 6. Error cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: --version"
stb-herRefs --version > log_version.txt 2>&1
check_contains "stb-herRefs" log_version.txt

echo "Testing: --help documents --zpe-mode/--displacement/--supercell"
stb-herRefs --help > log_help.txt 2>&1
check_contains "zpe-mode" log_help.txt
check_contains "displacement" log_help.txt
check_contains "supercell" log_help.txt


# --- 7. Interactive path (stb-suite, shortcut 4.13.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.13.2) ---"

echo "Testing: navigate 4.13.2 -> zpe-mode local -> quit"
make_site_disp
{
  echo "4.13.2"
  echo ""               # run_dir (default her_study)
  echo ""               # output_filename (default calc.out)
  echo "3"              # pseudo_dir -> option 3 = Custom path
  echo "."              # custom pseudo path
  echo "local"          # zpe_mode
  echo ""               # show_advanced (default N)
  echo ""               # press enter to continue
  echo "0"              # quit stage submenu
} | stb-suite > log_menu.txt 2>&1
check_success her_study/00_clean_slab/structure.fdf
check_success her_study/05_zpe_calc/disp_001/structure.fdf


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
