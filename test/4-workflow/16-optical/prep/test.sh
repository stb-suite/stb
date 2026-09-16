#!/bin/bash

# --- Setup ---
# Smoke test for stb-optical (Optical Properties Stage 1: Direction
# Folders, item 4.16.1). 2 fixtures: bulk.fdf (3D, simple-cubic NaCl --
# exercises the diagonal/biaxial direction folders and the config_extra.fdf
# override mechanism) and slab.fdf (2D graphene monolayer -- exercises the
# vacuum-padded [KNOWN LIMITATION] dimensionality note).
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$FIXTURE_DIR/test_files"

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
NC='\033[0m'

PASS=0
FAIL=0

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

# --- 1. Preparation ---
echo "--- Starting tester for STB-OPTICAL stage 1: direction folders (item 4.16.1) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/bulk.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/slab.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/calc.fdf" "$TEST_DIR/"
for sym in Na Cl C; do
    echo "# placeholder pseudopotential" > "$TEST_DIR/$sym.psf"
done
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Default behavior: -c/--directions both default, config_extra.fdf mechanism ---
echo -e "\n--- Testing default behavior (bulk.fdf, no -c, no --directions) ---"
rm -rf study_default
stb-optical -f bulk.fdf -O study_default --no-intro > log_default.txt 2>&1
check_exit_code $? 0
check_success study_default/dir_xx/structure.fdf
check_success study_default/dir_yy/structure.fdf
check_success study_default/dir_zz/structure.fdf
check_success study_default/dir_xx/calc.fdf
check_success study_default/dir_xx/config_extra.fdf

echo "Testing: -c defaults to calc.fdf in the current directory (no error despite no -c given)"
check_not_contains "\[ERROR\]" log_default.txt

echo "Testing: calc.fdf is the untouched user template, %include'ing config_extra.fdf at the top"
check_contains "%include config_extra.fdf" study_default/dir_xx/calc.fdf
check_contains "MD.Steps 150" study_default/dir_xx/calc.fdf

echo "Testing: config_extra.fdf carries the forced single-point + Optical.* block, not calc.fdf"
check_contains "MD.TypeOfRun          CG" study_default/dir_xx/config_extra.fdf
check_contains "MD.Steps              0" study_default/dir_xx/config_extra.fdf
check_contains "OpticalCalculation T" study_default/dir_xx/config_extra.fdf
check_contains "1.0000  0.0000  0.0000" study_default/dir_xx/config_extra.fdf
check_contains "0.0000  1.0000  0.0000" study_default/dir_yy/config_extra.fdf
check_contains "0.0000  0.0000  1.0000" study_default/dir_zz/config_extra.fdf

echo "Testing: calc.fdf is byte-identical across direction folders (shared template)"
diff -q study_default/dir_xx/calc.fdf study_default/dir_yy/calc.fdf > /dev/null 2>&1
check_exit_code $? 0

echo "Testing: config_extra.fdf DIFFERS across direction folders (direction-specific Optical.Vector)"
diff -q study_default/dir_xx/config_extra.fdf study_default/dir_yy/config_extra.fdf > /dev/null 2>&1
[ $? -ne 0 ] && echo -e "   -> ${GREEN}Verified:${NC} config_extra.fdf differs between dir_xx and dir_yy" && PASS=$((PASS+1)) \
             || { echo -e "   -> ${RED}Failed:${NC} config_extra.fdf unexpectedly identical between dir_xx and dir_yy"; FAIL=$((FAIL+1)); }

echo "Testing: report documents formula/cell volume (regression guard: LatticeConstant must be applied)"
check_contains "Formula         : Na1 Cl1" log_default.txt
check_contains "Cell volume     : 64.0000 Ang\^3" log_default.txt

echo "Testing: 3D bulk input -- no 2D known-limitation warning, dimensionality reported correctly"
check_contains "Detected : 3D (bulk material)" log_default.txt
check_not_contains "\[KNOWN LIMITATION\] Vacuum-padded" log_default.txt

echo "Testing: summary line counts diagonal/biaxial directions correctly"
check_contains "3 direction folder(s) written under 'study_default' (3 diagonal, 0 biaxial)" log_default.txt


# --- 3. Biaxial direction + missing-diagonal-pair warning ---
echo -e "\n--- Testing biaxial direction xy without its yy pair (missing-pair warning) ---"
rm -rf study_xy
stb-optical -f bulk.fdf -c calc.fdf -O study_xy --directions xx xy --no-intro > log_xy.txt 2>&1
check_exit_code $? 0
check_success study_xy/dir_xy/structure.fdf
check_success study_xy/dir_xy/config_extra.fdf
check_contains "0.7071  0.7071  0.0000" study_xy/dir_xy/config_extra.fdf
check_contains "missing yy in this run" log_xy.txt
check_contains "(1 diagonal, 1 biaxial)" log_xy.txt


# --- 4. Full 6-direction run (diagonal + biaxial complete set) ---
echo -e "\n--- Testing full 6-direction run (xx yy zz xy xz yz) ---"
rm -rf study_full
stb-optical -f bulk.fdf -c calc.fdf -O study_full --directions xx yy zz xy xz yz --no-intro \
    > log_full.txt 2>&1
check_exit_code $? 0
check_contains "(3 diagonal, 3 biaxial)" log_full.txt
check_success study_full/dir_xz/config_extra.fdf
check_contains "0.7071  0.0000  0.7071" study_full/dir_xz/config_extra.fdf
check_success study_full/dir_yz/config_extra.fdf
check_contains "0.0000  0.7071  0.7071" study_full/dir_yz/config_extra.fdf

echo "Testing: dielectric-tensor mapping table documents direct vs. needs-reconstruction"
check_contains "eps_xx (direct)" log_full.txt
check_contains "eps_xy (needs reconstruction)" log_full.txt

echo "Testing: with all pairs present, no missing-pair warning and the automation note is shown"
check_not_contains "missing" log_full.txt
check_contains "stb-opticalAnalysis applies this automatically" log_full.txt


# --- 5. 2D input: known-limitation dimensionality note ---
echo -e "\n--- Testing 2D fixture (slab.fdf, vacuum-padded) ---"
rm -rf study_slab
stb-optical -f slab.fdf -c calc.fdf -O study_slab --no-intro > log_slab.txt 2>&1
check_exit_code $? 0
check_contains "Detected : 2D (e.g., a slab or surface)" log_slab.txt
check_contains "\[KNOWN LIMITATION\] Vacuum-padded (2D/slab) input" log_slab.txt
echo "Testing: out-of-plane (zz) direction is still written for a 2D input, not skipped"
check_success study_slab/dir_zz/structure.fdf


# --- 6. Pseudopotentials copied into every direction folder ---
echo -e "\n--- Testing pseudopotential copying (-p .) ---"
rm -rf study_pseudo
stb-optical -f bulk.fdf -c calc.fdf -p . -O study_pseudo --no-intro > log_pseudo.txt 2>&1
check_exit_code $? 0
check_success study_pseudo/dir_xx/Na.psf
check_success study_pseudo/dir_xx/Cl.psf


# --- 7. Error cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: --version"
stb-optical --version > log_version.txt 2>&1
check_contains "stb-optical" log_version.txt

echo "Testing: --help documents the diagonal/biaxial naming and config_extra.fdf convention"
stb-optical --help > log_help.txt 2>&1
check_contains "xx yy zz" log_help.txt
check_contains "BIAXIAL" log_help.txt
check_contains "config_extra.fdf" log_help.txt

echo "Testing: missing structure file is rejected"
stb-optical -f does_not_exist.fdf -c calc.fdf --no-intro > log_bad_structure.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_bad_structure.txt

echo "Testing: missing --calc file is rejected"
stb-optical -f bulk.fdf -c does_not_exist.fdf --no-intro > log_bad_calc.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_bad_calc.txt

echo "Testing: default -c (calc.fdf) missing in the current directory is rejected the same way"
mkdir -p no_calc_here
( cd no_calc_here && cp ../bulk.fdf . && stb-optical -f bulk.fdf --no-intro > ../log_no_calc.txt 2>&1 )
check_exit_code $? 1
check_contains "not found" log_no_calc.txt

echo "Testing: an unrecognized --directions value is rejected by argparse"
stb-optical -f bulk.fdf -c calc.fdf --directions bogus --no-intro > log_bad_direction.txt 2>&1
check_exit_code $? 2
check_contains "invalid choice" log_bad_direction.txt


# --- 8. Interactive path (stb-suite, shortcut 4.16.1) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.16.1) ---"

echo "Testing: navigate 4.16.1 -> bulk.fdf -> defaults -> quit"
rm -rf study_menu
{
  echo "4.16.1"
  echo "bulk.fdf"        # structure file
  echo ""                 # calc.fdf (default)
  echo ""                  # pseudo source (skip)
  echo ""                   # directions (default xx yy zz)
  echo "study_menu"          # output dir
  echo "n"                    # advanced settings
  echo ""                      # press enter to continue
  echo "0"                      # quit stage submenu
} | stb-suite > log_menu.txt 2>&1
check_success study_menu/dir_xx/structure.fdf
check_success study_menu/dir_yy/structure.fdf
check_success study_menu/dir_zz/structure.fdf

echo "Testing: interactive-menu and direct-CLI results agree (same geometry)"
python3 -c "
import sys
from stb.core import structure_io
import numpy as np
a = structure_io.read_fdf('study_menu/dir_xx/structure.fdf')
b = structure_io.read_fdf('study_default/dir_xx/structure.fdf')
a_frac = sorted([tuple(np.round(pos, 6)) for _s, pos in a.atoms])
b_frac = sorted([tuple(np.round(pos, 6)) for _s, pos in b.atoms])
sys.exit(0 if a_frac == b_frac else 1)
" > log_menu_check.txt 2>&1
check_exit_code $? 0


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
