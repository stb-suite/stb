#!/bin/bash

# --- Setup ---
# Smoke test for stb-ramanModes (Raman Spectrum, Stage 2: Modes & Optical
# Displacements, item 4.11.2). Chains a real stb-raman (Stage 1) run to get
# a genuine raman_study/phonon_disp/ tree, then fabricates disp-*/Sn3O4.FA
# force files (same 14-atom Sn3O4 system/format as
# test/4-workflow/4-phonons/analysis/Sn3O4.FA -- physically meaningless but
# format-correct, enough to exercise FORCE_SETS extraction and the
# Gamma-mode/Optical-displacement wiring).
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

check_exit_code() {
    if [ "$1" -eq "$2" ]; then
        echo -e "   -> ${GREEN}Verified:${NC} exit code $1 (expected $2)"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} exit code $1 (expected $2)"
        FAIL=$((FAIL+1))
    fi
}

# Runs stb-raman (Stage 1, tested separately) then drops a fabricated,
# per-folder-scaled Sn3O4.FA into every disp-*/ -- same "each folder gets a
# slightly different scale factor" trick as test/4-workflow/4-phonons/
# analysis/test.sh, for the same reason: identical forces everywhere make
# the dynamical matrix perfectly degenerate and crash phonopy's own code.
make_phonon_disp() {
    rm -rf raman_study
    stb-raman -s structure.fdf -c calc.fdf -p . -dim 1 1 1 --no-intro > /dev/null 2>&1
    python3 - <<'PYEOF'
import glob
with open("Sn3O4.FA") as f:
    lines = f.readlines()
header, rows = lines[0], lines[1:]
for i, d in enumerate(sorted(glob.glob("raman_study/phonon_disp/disp-*")), start=1):
    scale = 1.0 + 0.01 * (i % 17)
    with open(f"{d}/Sn3O4.FA", "w") as out:
        out.write(header)
        for row in rows:
            parts = row.split()
            idx = parts[0]
            vals = [float(v) * scale for v in parts[1:]]
            out.write(f"{idx:>6}  {vals[0]: .9E}  {vals[1]: .9E}  {vals[2]: .9E}\n")
PYEOF
}


# --- 1. Preparation ---
echo "--- Starting tester for STB-RamanModes stage 2: modes & optical displacements (item 4.11.2) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$PREP_DIR/structure.fdf" "$TEST_DIR/"
cp "$PREP_DIR/calc.fdf" "$TEST_DIR/"
cp "$PREP_DIR/O.psf" "$TEST_DIR/"
cp "$PREP_DIR/Sn.psf" "$TEST_DIR/"
cp "$FIXTURE_DIR/Sn3O4.FA" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Missing Stage-1 output guard ---
echo -e "\n--- Testing the missing-Stage-1-output guard ---"
rm -rf raman_study
stb-ramanModes --directory raman_study --calc calc.fdf --no-intro > log_no_stage1.txt 2>&1
check_exit_code $? 1
check_contains "run stb-raman" log_no_stage1.txt


# --- 3. Normal generation: 3 selected modes ---
echo -e "\n--- Testing default generation (--modes 1 2 3) ---"
make_phonon_disp
stb-ramanModes --directory raman_study --calc calc.fdf --modes 1 2 3 --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "Force files found : 84" log_basic.txt
check_contains "FORCE_SETS       : generated successfully" log_basic.txt
check_contains "Non-acoustic Gamma modes : 39" log_basic.txt
check_contains "Selected modes    : 3" log_basic.txt
check_contains "Folders to write  : 3 modes x 2 signs x 3 axes = 18" log_basic.txt
check_success raman_study/optical_disp/mode_01_plus_x/structure.fdf
check_success raman_study/optical_disp/mode_01_plus_x/calc.fdf
check_success raman_study/optical_disp/mode_01_minus_z/structure.fdf
check_success raman_study/optical_disp/mode_01_plus_x/O.psf
check_success raman_study/optical_disp/mode_01_plus_x/Sn.psf
check_contains "OpticalCalculation T" raman_study/optical_disp/mode_01_plus_x/calc.fdf
check_contains "Optical.Vector" raman_study/optical_disp/mode_01_plus_x/calc.fdf
check_contains "MD.TypeOfRun.*CG" raman_study/optical_disp/mode_01_plus_x/calc.fdf
check_contains "MD.NumCGsteps.*0" raman_study/optical_disp/mode_01_plus_x/calc.fdf
check_contains "MODE_TABLE" raman_study/raman_stage2.txt
check_contains "\[2\] PSEUDOPOTENTIALS" raman_study/raman_stage2.txt
check_contains "Elements needed   : O, Sn" raman_study/raman_stage2.txt
check_contains "reused from Stage 1's disp-\* folders" raman_study/raman_stage2.txt
check_contains "Found all required" raman_study/raman_stage2.txt
check_contains "O.psf" raman_study/raman_stage2.txt
check_contains "Sn.psf" raman_study/raman_stage2.txt

echo "Testing: mode_01_plus_x and mode_01_plus_y share the SAME displaced structure (only Optical.Vector differs)"
python3 -c "
from stb.core import structure_io
import numpy as np
sx = structure_io.read_fdf('raman_study/optical_disp/mode_01_plus_x/structure.fdf')
sy = structure_io.read_fdf('raman_study/optical_disp/mode_01_plus_y/structure.fdf')
px = np.array([p for _, p in sx.atoms])
py = np.array([p for _, p in sy.atoms])
assert np.allclose(px, py), 'plus_x and plus_y structures should be identical'
print('OK')
" > log_structure_check.txt 2>&1
check_contains "OK" log_structure_check.txt

echo "Testing: mode_01_plus_x and mode_01_minus_x differ (opposite displacement sign), same cell"
python3 -c "
from stb.core import structure_io
import numpy as np
sp = structure_io.read_fdf('raman_study/optical_disp/mode_01_plus_x/structure.fdf')
sm = structure_io.read_fdf('raman_study/optical_disp/mode_01_minus_x/structure.fdf')
pp = np.array([p for _, p in sp.atoms])
pm = np.array([p for _, p in sm.atoms])
assert not np.allclose(pp, pm), 'plus_x and minus_x should differ'
assert np.allclose(sp.lattice, sm.lattice), 'cell should be unchanged by the mode displacement'
print('OK')
" > log_sign_check.txt 2>&1
check_contains "OK" log_sign_check.txt


# --- 4. --freq-min/--freq-max filtering ---
echo -e "\n--- Testing --freq-min/--freq-max filtering ---"
make_phonon_disp
stb-ramanModes --directory raman_study --calc calc.fdf --freq-min 1.0 --no-intro > log_freqfilter.txt 2>&1
check_exit_code $? 0
check_contains "Selected modes    : 1" log_freqfilter.txt


# --- 4b. Missing pseudopotentials in Stage 1's disp-* folders -> hard error, not silent ---
echo -e "\n--- Testing the missing-pseudopotential guard (no silent 0-pseudo folders) ---"
make_phonon_disp
rm -f raman_study/phonon_disp/disp-*/O.psf raman_study/phonon_disp/disp-*/Sn.psf
stb-ramanModes --directory raman_study --calc calc.fdf --modes 1 --no-intro > log_missing_pseudo.txt 2>&1
check_exit_code $? 1
check_contains "Missing pseudopotential" log_missing_pseudo.txt
check_contains "O, Sn" log_missing_pseudo.txt

echo "Testing: -p/--pseudo-dir override recovers from that same missing-pseudo situation"
stb-ramanModes --directory raman_study --calc calc.fdf --modes 1 -p . --no-intro > log_pseudo_override.txt 2>&1
check_exit_code $? 0
check_contains "-p/--pseudo-dir override" log_pseudo_override.txt
check_success raman_study/optical_disp/mode_01_plus_x/O.psf
check_success raman_study/optical_disp/mode_01_plus_x/Sn.psf


# --- 5. --full-tensor: full symmetric Raman tensor (6 axes instead of 3) ---
echo -e "\n--- Testing --full-tensor ---"
make_phonon_disp
stb-ramanModes --directory raman_study --calc calc.fdf --modes 1 --full-tensor --no-intro > log_fulltensor.txt 2>&1
check_exit_code $? 0
check_contains "Tensor scope      : full symmetric tensor" log_fulltensor.txt
check_contains "Folders to write  : 1 modes x 2 signs x 6 axes = 12" log_fulltensor.txt
check_success raman_study/optical_disp/mode_01_plus_xy/structure.fdf
check_success raman_study/optical_disp/mode_01_plus_xz/structure.fdf
check_success raman_study/optical_disp/mode_01_plus_yz/structure.fdf
check_success raman_study/optical_disp/mode_01_minus_xy/structure.fdf

echo "Testing: mode_01_plus_xy uses the (x+y)/sqrt(2) Optical.Vector direction"
check_contains "0.7071  0.7071  0.0000" raman_study/optical_disp/mode_01_plus_xy/calc.fdf


# --- 6. Error cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: missing calc file"
stb-ramanModes --directory raman_study --calc does_not_exist.fdf --no-intro > log_missing_calc.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_calc.txt

echo "Testing: --version"
stb-ramanModes --version > log_version.txt 2>&1
check_contains "stb-ramanModes" log_version.txt

echo "Testing: --help documents --modes/--optical-mesh/--optical-broaden/--full-tensor"
stb-ramanModes --help > log_help.txt 2>&1
check_contains "modes" log_help.txt
check_contains "optical-mesh" log_help.txt
check_contains "optical-broaden" log_help.txt
check_contains "full-tensor" log_help.txt


# --- 7. Interactive path (stb-suite, shortcut 4.11.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.11.2) ---"

echo "Testing: navigate 4.11.2 -> modes 1 2, skip full-tensor -> quit"
make_phonon_disp
{
  echo "4.11.2"
  echo "raman_study"    # run_dir
  echo "calc.fdf"       # calc_file
  echo "1 2"            # modes
  echo ""               # full_tensor_choice (default N)
  echo ""               # displacement (default 0.02)
  echo ""               # show_advanced (default -> skip)
  echo ""               # press enter to continue
  echo "0"              # quit stage submenu
} | stb-suite > log_menu.txt 2>&1
check_contains "Folders to write  : 2 modes x 2 signs x 3 axes = 12" log_menu.txt
check_success raman_study/optical_disp/mode_01_plus_x/structure.fdf

echo "Testing: navigate 4.11.2 -> modes 1, full-tensor ON -> quit"
make_phonon_disp
{
  echo "4.11.2"
  echo "raman_study"    # run_dir
  echo "calc.fdf"       # calc_file
  echo "1"              # modes
  echo "y"              # full_tensor_choice -> Y
  echo ""               # displacement (default 0.02)
  echo ""               # show_advanced (default -> skip)
  echo ""               # press enter to continue
  echo "0"              # quit stage submenu
} | stb-suite > log_menu_fulltensor.txt 2>&1
check_contains "Folders to write  : 1 modes x 2 signs x 6 axes = 12" log_menu_fulltensor.txt
check_success raman_study/optical_disp/mode_01_plus_xy/structure.fdf


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
