#!/bin/bash

# --- Setup ---
# Smoke test for stb-her (HER Stage 1: Adsorption Sites, item 4.13.1).
# Reuses the same graphene fixture as test/4-workflow/8-adsorption/prep/
# (free-standing 2D, so it also exercises --both-sides) -- stb-her is a
# leaner, H-only rewrite of the same pymatgen AdsorbateSiteFinder pattern.
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
    if grep -q -- "$1" "$2" 2>/dev/null; then
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
echo "--- Starting tester for STB-HER stage 1: adsorption sites (item 4.13.1) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/structure.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/calc.fdf" "$TEST_DIR/"
echo "# placeholder pseudopotential" > "$TEST_DIR/C.psf"
echo "# placeholder pseudopotential" > "$TEST_DIR/H.psf"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Default run (--site-type all) ---
echo -e "\n--- Testing default generation (--site-type all) ---"
stb-her -s structure.fdf -c calc.fdf -p . --height 1.2 --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_success her_study/clean_slab_source/structure.fdf
check_success her_study/sites/site_1_ontop/structure.fdf
check_success her_study/sites/site_1_ontop/calc.fdf
check_success her_study/sites/site_1_ontop/C.psf
check_success her_study/sites/site_1_ontop/H.psf
check_success her_study/her_stage1.txt
check_contains "\[0\] RUN METADATA" her_study/her_stage1.txt
check_contains "\[1\] CLEAN SLAB REFERENCE" her_study/her_stage1.txt
check_contains "\[2\] ADSORPTION SITES" her_study/her_stage1.txt
check_contains "\[3\] SUMMARY & NEXT STEPS" her_study/her_stage1.txt
check_contains "stb-herRefs --directory her_study" her_study/her_stage1.txt

echo "Testing: H is appended as the LAST atom (3rd), Slab.DipoleCorrection forced on"
check_contains "NumberofAtoms      3" her_study/sites/site_1_ontop/structure.fdf
check_contains "Slab.DipoleCorrection   T" her_study/sites/site_1_ontop/calc.fdf

echo "Testing: SITE_TABLE written, only 2 columns (label, dir -- no adsorbate/height column, always H)"
check_contains "# SITE_TABLE" her_study/her_stage1.txt
check_contains "site_1_ontop" her_study/her_stage1.txt


# --- 3. --site-type ontop, single site, explicit height ---
echo -e "\n--- Testing --site-type ontop (single site type on graphene's 1 distinct C site) ---"
rm -rf her_study
stb-her -s structure.fdf -c calc.fdf -p . --site-type ontop --height 1.5 --no-intro \
    > log_ontop.txt 2>&1
check_exit_code $? 0
python3 -c "
import glob
dirs = sorted(glob.glob('her_study/sites/site_*'))
assert all('ontop' in d for d in dirs), f'expected only ontop sites, got {dirs}'
print('OK')
" > log_ontop_check.txt 2>&1
check_contains "OK" log_ontop_check.txt


# --- 4. --both-sides (free-standing 2D graphene) ---
echo -e "\n--- Testing --both-sides ---"
rm -rf her_study
stb-her -s structure.fdf -c calc.fdf -p . --site-type ontop --both-sides --no-intro \
    > log_bothsides.txt 2>&1
check_exit_code $? 0
check_contains "bothsides" log_bothsides.txt

echo "Testing: --both-sides requires a concrete --site-type"
rm -rf her_study
stb-her -s structure.fdf -c calc.fdf -p . --site-type all --both-sides --no-intro \
    > log_bothsides_all.txt 2>&1
check_exit_code $? 1
check_contains "concrete --site-type" log_bothsides_all.txt


# --- 5. Overlap warning (unrealistically small height) ---
echo -e "\n--- Testing the too-close warning (height too small) ---"
rm -rf her_study
stb-her -s structure.fdf -c calc.fdf -p . --site-type ontop --height 0.1 --no-intro \
    > log_overlap.txt 2>&1
check_exit_code $? 0
check_contains "likely overlapping atoms" log_overlap.txt


# --- 6. Error cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: missing structure file"
rm -rf her_study
stb-her -s does_not_exist.fdf -c calc.fdf -p . --no-intro > log_missing_struct.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_struct.txt

echo "Testing: bulk 3D structure rejected (no single vacuum axis)"
cat > bulk3d.fdf << 'FDFEOF'
NumberOfSpecies    1
NumberofAtoms      1

%block ChemicalSpeciesLabel
 1   6   C
%endblock ChemicalSpeciesLabel

LatticeConstant 1.0 Ang
AtomicCoordinatesFormat  Fractional

%block LatticeVectors
  3.0   0.0   0.0
  0.0   3.0   0.0
  0.0   0.0   3.0
%endblock LatticeVectors

%block AtomicCoordinatesAndAtomicSpecies
  0.0   0.0   0.0   1
%endblock AtomicCoordinatesAndAtomicSpecies
FDFEOF
stb-her -s bulk3d.fdf -c calc.fdf -p . --no-intro > log_bulk3d.txt 2>&1
check_exit_code $? 1
check_contains "vacuum along exactly one axis" log_bulk3d.txt

echo "Testing: --version"
stb-her --version > log_version.txt 2>&1
check_contains "stb-her" log_version.txt

echo "Testing: --help documents --site-type/--height/--both-sides/--symprec"
stb-her --help > log_help.txt 2>&1
check_contains "site-type" log_help.txt
check_contains "height" log_help.txt
check_contains "both-sides" log_help.txt
check_contains "symprec" log_help.txt


# --- 7. Interactive path (stb-suite, shortcut 4.13.1) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.13.1) ---"

echo "Testing: navigate 4.13.1 -> defaults -> quit"
rm -rf her_study
{
  echo "4.13.1"
  echo "structure.fdf"  # structure_file
  echo "calc.fdf"       # calc_file
  echo "3"              # pseudo_dir -> option 3 = Custom path (dojo=1, virtual_vault=2, custom=3)
  echo "."              # custom pseudo path
  echo ""               # site_type (default all)
  echo ""               # height (default 1.5)
  echo ""               # both_sides_choice (default N)
  echo ""               # output_dir (default her_study)
  echo ""               # show_advanced (default N)
  echo ""               # press enter to continue
  echo "0"              # quit stage submenu
} | stb-suite > log_menu.txt 2>&1
check_success her_study/sites/site_1_ontop/structure.fdf


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
