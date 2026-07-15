#!/bin/bash

# --- Setup ---
# Smoke test for stb-neb (NEB Prep, item 4.9.1)
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
echo "--- Starting tester for STB-Neb prep (item 4.9.1) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/initial.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/final.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/final_bad_composition.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/final_bad_lattice.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/wrap_initial.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/wrap_final.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/calc.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Default linear interpolation ---
echo -e "\n--- Testing default linear interpolation, -n 5 ---"
stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 5 --no-intro > log_linear.txt 2>&1
check_exit_code $? 0
check_success image_00/structure.fdf
check_success image_00/calc.fdf
check_success image_04/structure.fdf
check_contains "MD.TypeOfRun.*CG" image_02/calc.fdf
check_contains "MD.NumCGsteps.*0" image_02/calc.fdf
check_success neb_setup.txt
check_contains "\[0\] RUN METADATA" neb_setup.txt
check_contains "\[1\] INTERPOLATION" neb_setup.txt
check_contains "\[2\] IMAGE FOLDERS" neb_setup.txt
check_contains "\[3\] SUMMARY" neb_setup.txt
check_contains "# ML_NEB_USED: no" neb_setup.txt
check_contains "IMAGE_TABLE" neb_setup.txt
check_success neb_path.xyz
if grep -q "WARNING" log_linear.txt; then
    echo -e "   -> ${RED}Failed:${NC} unexpected path-quality WARNING for a well-behaved 5-image band"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no path-quality warnings for a well-behaved band"
    PASS=$((PASS+1))
fi

echo "Testing: reaction_coord is non-decreasing and image_00 is 0.0"
python3 -c "
import re
rows = []
in_table = False
with open('neb_setup.txt') as f:
    for line in f:
        if line.startswith('# IMAGE_TABLE'):
            in_table = True
            continue
        if not in_table:
            continue
        if line.startswith('#') or not line.strip():
            continue
        parts = line.split()
        rows.append((parts[0], int(parts[1]), float(parts[2])))
rows.sort(key=lambda r: r[1])
assert rows[0][2] == 0.0, f'image_00 reaction_coord should be 0.0, got {rows[0][2]}'
coords = [r[2] for r in rows]
assert coords == sorted(coords), f'reaction_coord not non-decreasing: {coords}'
print('OK')
" > log_coord_check.txt 2>&1
check_contains "OK" log_coord_check.txt


# --- 3. Composition mismatch (hard error) ---
echo -e "\n--- Testing composition mismatch (hard error) ---"
stb-neb -i initial.fdf -f final_bad_composition.fdf -c calc.fdf --no-intro > log_bad_composition.txt 2>&1
check_exit_code $? 1
check_contains "different composition" log_bad_composition.txt


# --- 4. Lattice mismatch (warning, lattice override) ---
echo -e "\n--- Testing lattice mismatch (warning, not error) ---"
rm -rf image_* neb_setup.txt
stb-neb -i initial.fdf -f final_bad_lattice.fdf -c calc.fdf -n 5 --no-intro > log_bad_lattice.txt 2>&1
check_exit_code $? 0
check_contains "WARNING.*different lattices" log_bad_lattice.txt
python3 -c "
from stb.core import structure_io
initial = structure_io.read_fdf('initial.fdf')
final_image = structure_io.read_fdf('image_04/structure.fdf')
import numpy as np
assert np.allclose(initial.lattice, final_image.lattice), 'image_04 lattice should match initial.fdf, not final_bad_lattice.fdf'
print('OK')
" > log_lattice_check.txt 2>&1
check_contains "OK" log_lattice_check.txt


# --- 4b. Coordinate-wrapping bug regression: two endpoints describing the
#     SAME physical position via different fractional-wrapping conventions
#     (frac z=0.95 vs z=-0.05 in a 10 Ang cell) must interpolate to a
#     near-static path, not a spurious ~10 Ang round trip through the cell ---
echo -e "\n--- Testing coordinate-wrapping regression (same point, different wrapping) ---"
rm -rf image_* neb_setup.txt
stb-neb -i wrap_initial.fdf -f wrap_final.fdf -c calc.fdf -n 5 --no-intro > log_wrap.txt 2>&1
check_exit_code $? 0
python3 -c "
from stb.core import structure_io
for i in range(5):
    s = structure_io.read_fdf(f'image_0{i}/structure.fdf')
    z = s.atoms[0][1][2]
    assert abs(z - 0.95) < 1e-6, f'image_0{i}: expected frac z ~0.95 (same wrapped point), got {z}'
print('OK')
" > log_wrap_check.txt 2>&1
check_contains "OK" log_wrap_check.txt


# --- 5. -n minimum ---
echo -e "\n--- Testing -n minimum (n=2) ---"
stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 2 --no-intro > log_nmin.txt 2>&1
check_exit_code $? 2


# --- 6. --idpp ---
echo -e "\n--- Testing --idpp ---"
rm -rf image_* neb_setup.txt
stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 5 --idpp --no-intro > log_idpp.txt 2>&1
check_exit_code $? 0
check_contains "IDPP" log_idpp.txt
check_success image_02/structure.fdf


# --- 7. --ml-neb (MACE-MP-0, cached locally) ---
echo -e "\n--- Testing --ml-neb ---"
rm -rf image_* neb_setup.txt
stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 5 --ml-neb --ml-max-steps 30 --no-intro > log_mlneb.txt 2>&1
check_exit_code $? 0
check_contains "ML-NEB" log_mlneb.txt
check_contains "barrier" log_mlneb.txt
check_success neb_ml_preview.png
check_success neb_path.xyz
check_contains "# ML_NEB_USED: yes" neb_setup.txt
check_contains "ML-NEB freeze substrate: yes" neb_setup.txt
check_contains "Freezing 8/9 atom" log_mlneb.txt


# --- 7b. --no-ml-freeze-substrate (freeze-substrate is default ON) ---
echo -e "\n--- Testing --no-ml-freeze-substrate ---"
rm -rf image_* neb_setup.txt neb_ml_preview.png neb_path.xyz
stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 5 --ml-neb --ml-max-steps 10 \
    --no-ml-freeze-substrate --no-intro > log_nofreeze.txt 2>&1
check_exit_code $? 0
check_contains "ML-NEB freeze substrate: no" log_nofreeze.txt
if grep -q "Freezing" log_nofreeze.txt; then
    echo -e "   -> ${RED}Failed:${NC} unexpected 'Freezing' message with --no-ml-freeze-substrate"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no atoms frozen with --no-ml-freeze-substrate"
    PASS=$((PASS+1))
fi


# --- 7c. Path-quality warning (deliberately dense -n triggers "nearly identical") ---
echo -e "\n--- Testing path-quality warning (dense -n) ---"
rm -rf image_* neb_setup.txt neb_ml_preview.png neb_path.xyz
stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 60 --no-intro > log_dense.txt 2>&1
check_exit_code $? 0
check_contains "WARNING.*nearly identical" log_dense.txt


# --- 8. --ml-prerelax-endpoints without --ml-neb (independence check) ---
echo -e "\n--- Testing --ml-prerelax-endpoints without --ml-neb ---"
rm -rf image_* neb_setup.txt neb_ml_preview.png
stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 5 --ml-prerelax-endpoints --no-intro > log_prerelax_endpoints.txt 2>&1
check_exit_code $? 0
check_contains "ML pre-relax" log_prerelax_endpoints.txt
check_contains "ML-NEB (MACE-MP-0): no" log_prerelax_endpoints.txt
if grep -q "running a real climbing-image NEB" log_prerelax_endpoints.txt; then
    echo -e "   -> ${RED}Failed:${NC} unexpected climbing-image NEB run without --ml-neb"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no climbing-image NEB run without --ml-neb (flags are independent)"
    PASS=$((PASS+1))
fi


# --- 9. Inert ML/IDPP flags without their master flag ---
echo -e "\n--- Testing inert ML/IDPP flags without --ml-neb/--idpp ---"
rm -rf image_* neb_setup.txt
stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 5 --ml-k 0.2 --ml-max-steps 50 \
    --idpp-fmax 0.2 --no-intro > log_inert_flags.txt 2>&1
check_exit_code $? 0


# --- 10. Missing files ---
echo -e "\n--- Testing missing files ---"
stb-neb -i does_not_exist.fdf -f final.fdf -c calc.fdf --no-intro > log_missing.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing.txt


# --- 11. --version / --help ---
echo -e "\n--- Testing --version / --help ---"
stb-neb --version > log_version.txt 2>&1
check_contains "stb-neb" log_version.txt
stb-neb --help > log_help.txt 2>&1
check_contains "n-images" log_help.txt
check_contains "ml-neb" log_help.txt
check_contains "ml-freeze-substrate" log_help.txt
check_contains "ml-freeze-threshold" log_help.txt


# --- 12. Interactive path (stb-suite, shortcut 4.9.1) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.9.1) ---"

echo "Testing: navigate 4.9.1 -> defaults -> quit"
rm -rf image_* neb_setup.txt
{
  echo "4.9.1"
  echo "initial.fdf"   # initial_file
  echo "final.fdf"     # final_file
  echo "calc.fdf"      # calc_file
  echo ""              # pp_path (skip)
  echo "5"             # n_images
  echo ""              # idpp_choice (default N)
  echo ""              # ml_neb_choice (default N)
  echo ""              # ml_prerelax_choice (default N)
  echo ""              # show_advanced (default -> skip)
  echo ""              # press enter to continue
  echo "0"             # quit stage submenu
} | stb-suite > log_menu.txt 2>&1
check_contains "Success:.*5 image folder" log_menu.txt
check_success image_00/structure.fdf
check_success image_04/structure.fdf

echo "Testing: navigate 4.9.1 -> --idpp on -> quit"
rm -rf image_* neb_setup.txt
{
  echo "4.9.1"
  echo "initial.fdf"
  echo "final.fdf"
  echo "calc.fdf"
  echo ""
  echo "5"
  echo "y"              # idpp_choice -> Y
  echo ""               # ml_neb_choice
  echo ""               # ml_prerelax_choice
  echo ""               # show_advanced
  echo ""               # press enter
  echo "0"
} | stb-suite > log_menu_idpp.txt 2>&1
check_contains "Success:.*5 image folder" log_menu_idpp.txt
check_success image_02/structure.fdf

echo "Testing: navigate 4.9.1 -> --ml-neb on, freeze-substrate default -> quit"
rm -rf image_* neb_setup.txt neb_ml_preview.png neb_path.xyz
{
  echo "4.9.1"
  echo "initial.fdf"
  echo "final.fdf"
  echo "calc.fdf"
  echo ""
  echo "5"
  echo ""               # idpp_choice
  echo "y"               # ml_neb_choice -> Y
  echo "0.1"             # ml_k
  echo "10"              # ml_max_steps (small, keep the test fast)
  echo ""                # freeze_choice (default Y)
  echo ""                # ml_freeze_threshold (default 0.3)
  echo ""                # ml_prerelax_choice
  echo ""                # show_advanced
  echo ""                # press enter
  echo "0"
} | stb-suite > log_menu_mlneb.txt 2>&1
check_contains "Success:.*5 image folder" log_menu_mlneb.txt
check_success neb_ml_preview.png
check_success neb_path.xyz
check_contains "ML-NEB freeze substrate.*ON" log_menu_mlneb.txt


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
