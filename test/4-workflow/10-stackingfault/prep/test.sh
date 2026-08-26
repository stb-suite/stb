#!/bin/bash

# --- Setup ---
# Smoke test for stb-stackingfault (Stacking Fault Prep, item 4.10.1)
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

check_missing() {
    if [ ! -e "$1" ]; then
        echo -e " ... ${GREEN}OK${NC} (file '$1' correctly NOT created)"
        PASS=$((PASS+1))
    else
        echo -e " ... ${RED}FAIL${NC} (file '$1' unexpectedly exists)"
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
echo "--- Starting tester for STB-Stackingfault prep (item 4.10.1) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/graphene.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/hbn.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/calc.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Default homobilayer grid (graphene/graphene, -nx 3 -ny 3), everything
#     nested under sf_run/, manifest always written, report only with
#     --save-report -- explicit --mode 3 (fixed gap, plain single-point) since
#     this block asserts single-point config_extra.fdf content, independent
#     of --mode's own default (see the dedicated default-mode test below) ---
echo -e "\n--- Testing homobilayer grid (graphene/graphene, -nx 3 -ny 3, --save-report) ---"
stb-stackingfault -l1 graphene.fdf -l2 graphene.fdf -c calc.fdf -nx 3 -ny 3 --mode 3 \
    --save-report --no-intro > log_homo.txt 2>&1
check_exit_code $? 0
check_success sf_run/shift_00_00/structure.fdf
check_success sf_run/shift_00_00/calc.fdf
check_success sf_run/shift_02_02/structure.fdf
check_success sf_run/shift_01_01/config_extra.fdf
check_contains "MD.TypeOfRun.*CG" sf_run/shift_01_01/config_extra.fdf
check_contains "MD.Steps.*0" sf_run/shift_01_01/config_extra.fdf
check_contains "%include config_extra.fdf" sf_run/shift_01_01/calc.fdf
echo "Testing: %include config_extra.fdf is the first line of calc.fdf"
first_line=$(head -n1 sf_run/shift_01_01/calc.fdf)
if [[ "$first_line" == *"%include config_extra.fdf"* ]]; then
    echo -e "   -> ${GREEN}Verified:${NC} config_extra.fdf is included first (overrides the template)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} first line was '$first_line'"
    FAIL=$((FAIL+1))
fi
echo "Testing: D3 is ON by default (no --d3/--no-d3 passed) -- config_extra.fdf has DFTD3"
check_contains "D3 dispersion   : yes" log_homo.txt
if grep -q "DFTD3" sf_run/shift_01_01/config_extra.fdf; then
    echo -e "   -> ${GREEN}Verified:${NC} DFTD3 tag present in config_extra.fdf by default"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} expected DFTD3 in config_extra.fdf by default"
    FAIL=$((FAIL+1))
fi

check_success sf_run/stackingfault_setup.txt
check_contains "\[0\] RUN METADATA" sf_run/stackingfault_setup.txt
check_contains "\[1\] ZSL MATCH" sf_run/stackingfault_setup.txt
check_contains "Each grid point: 4 atoms (2 layer 1 + 2 layer 2)" sf_run/stackingfault_setup.txt
check_contains "9 grid point(s) total" sf_run/stackingfault_setup.txt
check_contains "\[2\] GRID FOLDERS" sf_run/stackingfault_setup.txt
check_contains "\[3\] SUMMARY & NEXT STEPS" sf_run/stackingfault_setup.txt
check_contains "\[4\] LIBRARY WARNINGS" sf_run/stackingfault_setup.txt
check_contains "No library warnings\." sf_run/stackingfault_setup.txt
check_success sf_run/sf_manifest.json

echo "Testing: sf_manifest.json has exactly 9 rows (3x3), correct grid_nx/grid_ny, shift AA registry (0,0) present"
python3 -c "
import json
with open('sf_run/sf_manifest.json') as f:
    manifest = json.load(f)
assert manifest['grid_nx'] == 3 and manifest['grid_ny'] == 3, manifest
rows = manifest['rows']
assert len(rows) == 9, f'expected 9 grid rows, got {len(rows)}'
labels = {r['label'] for r in rows}
assert 'shift_00_00' in labels
print('OK')
" > log_manifest_check.txt 2>&1
check_contains "OK" log_manifest_check.txt

echo "Testing: every grid folder has the same atom count (only layer-2 position changes)"
python3 -c "
from stb.core import structure_io
import glob
counts = set()
for d in sorted(glob.glob('sf_run/shift_*')):
    s = structure_io.read_fdf(f'{d}/structure.fdf')
    counts.add(len(s.atoms))
assert len(counts) == 1, f'expected identical atom count across all grid points, got {counts}'
print('OK')
" > log_natoms_check.txt 2>&1
check_contains "OK" log_natoms_check.txt

echo "Testing: shift_01_01 layer-2 atoms are shifted by exactly (1/3, 1/3) relative to shift_00_00"
python3 -c "
from stb.core import structure_io
import numpy as np
s00 = structure_io.read_fdf('sf_run/shift_00_00/structure.fdf')
s11 = structure_io.read_fdf('sf_run/shift_01_01/structure.fdf')
# layer 2 = last 2 atoms (4-atom cell, 2 per layer)
l2_00 = np.array([pos for _, pos in s00.atoms[2:]])
l2_11 = np.array([pos for _, pos in s11.atoms[2:]])
delta = (l2_11 - l2_00) % 1.0
expected = np.array([1/3, 1/3, 0.0])
assert np.allclose(delta, expected, atol=1e-4), f'expected shift {expected}, got {delta}'
print('OK')
" > log_shift_check.txt 2>&1
check_contains "OK" log_shift_check.txt

echo "Testing: --no-d3 opts out of the DFT-D3 dispersion correction"
rm -rf sf_run
stb-stackingfault -l1 graphene.fdf -l2 graphene.fdf -c calc.fdf -nx 2 -ny 2 --mode 3 --no-d3 \
    --no-intro > log_nod3.txt 2>&1
check_exit_code $? 0
check_contains "D3 dispersion   : no" log_nod3.txt
if grep -q "DFTD3" sf_run/shift_00_00/config_extra.fdf; then
    echo -e "   -> ${RED}Failed:${NC} unexpected DFTD3 in config_extra.fdf with --no-d3"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no DFTD3 tag in config_extra.fdf with --no-d3"
    PASS=$((PASS+1))
fi


# --- 2b. Without --save-report, the manifest is still written but the
#     narrative report is not ---
echo -e "\n--- Testing that the manifest is always written, the report only with --save-report ---"
rm -rf sf_run
stb-stackingfault -l1 graphene.fdf -l2 graphene.fdf -c calc.fdf -nx 3 -ny 3 --no-intro \
    > log_noreport.txt 2>&1
check_exit_code $? 0
check_success sf_run/sf_manifest.json
check_missing sf_run/stackingfault_setup.txt
check_contains "Success:.*9 grid folder" log_noreport.txt


# --- 2c. Asymmetric grid (-nx 4 -ny 2) ---
echo -e "\n--- Testing an asymmetric grid (-nx 4 -ny 2) ---"
rm -rf sf_run
stb-stackingfault -l1 graphene.fdf -l2 graphene.fdf -c calc.fdf -nx 4 -ny 2 --save-report \
    --no-intro > log_asym.txt 2>&1
check_exit_code $? 0
check_contains "Grid            : 4 x 2 (asymmetric)" sf_run/stackingfault_setup.txt
check_contains "8 grid point(s) total" sf_run/stackingfault_setup.txt
check_success sf_run/shift_03_01/structure.fdf
check_missing sf_run/shift_00_02
python3 -c "
import json
with open('sf_run/sf_manifest.json') as f:
    manifest = json.load(f)
assert manifest['grid_nx'] == 4 and manifest['grid_ny'] == 2, manifest
assert len(manifest['rows']) == 8, manifest['rows']
print('OK')
" > log_asym_manifest_check.txt 2>&1
check_contains "OK" log_asym_manifest_check.txt


# --- 2d. --d3 forces DFTD3 in config_extra.fdf, alongside the single-point block ---
echo -e "\n--- Testing --d3 (Grimme DFT-D3 dispersion correction) ---"
rm -rf sf_run
stb-stackingfault -l1 graphene.fdf -l2 graphene.fdf -c calc.fdf -nx 3 -ny 3 --d3 --mode 3 \
    --save-report --no-intro > log_d3.txt 2>&1
check_exit_code $? 0
check_contains "D3 dispersion   : yes" log_d3.txt
check_contains "MD.TypeOfRun.*CG" sf_run/shift_00_00/config_extra.fdf
check_contains "MD.Steps.*0" sf_run/shift_00_00/config_extra.fdf
check_contains "DFTD3.*\.true\." sf_run/shift_00_00/config_extra.fdf
check_contains "%include config_extra.fdf" sf_run/shift_00_00/calc.fdf
if grep -q "DFTD3" sf_run/shift_00_00/calc.fdf; then
    echo -e "   -> ${RED}Failed:${NC} DFTD3 should live in config_extra.fdf, not calc.fdf directly"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} calc.fdf itself carries no DFTD3 tag (it's in config_extra.fdf)"
    PASS=$((PASS+1))
fi


# --- 3. Heterostructure grid (graphene/hBN) ---
echo -e "\n--- Testing heterostructure grid (graphene/hBN) ---"
rm -rf sf_run
stb-stackingfault -l1 graphene.fdf -l2 hbn.fdf -c calc.fdf -nx 3 -ny 3 --no-intro > log_hetero.txt 2>&1
check_exit_code $? 0
check_success sf_run/shift_00_00/structure.fdf
check_contains " 1   6   C" sf_run/shift_00_00/structure.fdf
check_contains " 2   5   B" sf_run/shift_00_00/structure.fdf
check_contains " 3   7   N" sf_run/shift_00_00/structure.fdf


# --- 3b. --ml-prerelax-layers (MACE-MP-0, cached locally) ---
echo -e "\n--- Testing --ml-prerelax-layers ---"
rm -rf sf_run
stb-stackingfault -l1 graphene.fdf -l2 graphene.fdf -c calc.fdf -nx 3 -ny 3 --ml-prerelax-layers \
    --save-report --no-intro > log_prerelax.txt 2>&1
check_exit_code $? 0
check_contains "ML pre-relax" log_prerelax.txt
check_contains "Converged (layer 1)" log_prerelax.txt
check_contains "Converged (layer 2)" log_prerelax.txt
check_success sf_run/shift_00_00/structure.fdf

# --- 3c. --ml-preview (MACE-MP-0, cached locally) ---
echo -e "\n--- Testing --ml-preview ---"
rm -rf sf_run
stb-stackingfault -l1 graphene.fdf -l2 graphene.fdf -c calc.fdf -nx 3 -ny 3 --ml-preview \
    --save-report --no-intro > log_preview.txt 2>&1
check_exit_code $? 0
check_success sf_run/stackingfault_ml_preview.png
check_contains "ML preview: predicted equilibrium at shift_" log_preview.txt
check_contains "predicted corrugation" log_preview.txt

# --- 3d. --mode 2 (MACE-MP-0 relaxes z, then SIESTA single-point; MACE-MP-0
#     cached locally) ---
echo -e "\n--- Testing --mode 2 (MACE relaxes z, then single-point SIESTA) ---"
rm -rf sf_run
stb-stackingfault -l1 graphene.fdf -l2 graphene.fdf -c calc.fdf -nx 3 -ny 3 --mode 2 --no-d3 \
    --save-report --no-intro > log_mode2.txt 2>&1
check_exit_code $? 0
check_contains "Mode            : 2" log_mode2.txt
check_contains "MACE-relaxed gap ranged" log_mode2.txt
check_contains "gap:.*Ang (ML-relaxed)" log_mode2.txt
check_contains "MD.TypeOfRun.*CG" sf_run/shift_00_00/config_extra.fdf
check_contains "MD.Steps.*0" sf_run/shift_00_00/config_extra.fdf
check_contains "ML model.*MACE-MP-0.*(no D3 dispersion)" log_mode2.txt

echo "Testing: [4] LIBRARY WARNINGS actually captures MACE/torch noise (not just 'No library warnings.')"
check_contains "\[4\] LIBRARY WARNINGS" log_mode2.txt
if grep -qE "MACE (calculator setup|z-relax)" log_mode2.txt; then
    echo -e "   -> ${GREEN}Verified:${NC} MACE-related entries present in the library warnings section"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} no MACE-related entries found in the library warnings section"
    FAIL=$((FAIL+1))
fi

echo "Testing: the written structure.fdf's actual interlayer gap matches the reported MACE-relaxed value"
python3 -c "
import re
from stb.core import structure_io
with open('log_mode2.txt') as f:
    text = f.read()
m = re.search(r'shift_00_00.*?gap: ([\d.]+) Ang \(ML-relaxed\)', text)
assert m, 'could not find shift_00_00 reported gap in log'
reported_gap = float(m.group(1))
s = structure_io.read_fdf('sf_run/shift_00_00/structure.fdf')
zs = sorted(set(round(p[2], 6) for _, p in s.atoms))
actual_gap = (zs[1] - zs[0]) * s.lattice[2][2]
assert abs(actual_gap - reported_gap) < 0.01, f'reported {reported_gap} vs actual {actual_gap}'
print('OK')
" > log_gap_check.txt 2>&1
check_contains "OK" log_gap_check.txt

echo "Testing: --mode 2 --d3 matches the MACE calculator's own dispersion to --d3"
rm -rf sf_run
stb-stackingfault -l1 graphene.fdf -l2 graphene.fdf -c calc.fdf -nx 2 -ny 2 --mode 2 --d3 \
    --no-intro > log_mode2_d3.txt 2>&1
check_exit_code $? 0
check_contains "ML model.*MACE-MP-0.*D3(BJ) dispersion (matches --d3)" log_mode2_d3.txt

# --- 3e. --mode 1 (SIESTA relaxes z for real -- restricted CG, x/y frozen) ---
echo -e "\n--- Testing --mode 1 (restricted SIESTA z-relaxation) ---"
rm -rf sf_run
stb-stackingfault -l1 graphene.fdf -l2 graphene.fdf -c calc.fdf -nx 2 -ny 2 --mode 1 \
    --relax-z-steps 50 --save-report --no-intro > log_mode1.txt 2>&1
check_exit_code $? 0
check_contains "Mode            : 1" log_mode1.txt
check_contains "NOT verified against a real SIESTA run" log_mode1.txt
check_contains "MD.TypeOfRun.*CG" sf_run/shift_00_00/config_extra.fdf
check_contains "MD.Steps.*50" sf_run/shift_00_00/config_extra.fdf
check_contains "MD.VariableCell        false" sf_run/shift_00_00/config_extra.fdf
check_contains "position from 1 to .* 1.0 0.0 0.0" sf_run/shift_00_00/config_extra.fdf
check_contains "position from 1 to .* 0.0 1.0 0.0" sf_run/shift_00_00/config_extra.fdf
if grep -qE "MD\.Steps\s+0$" sf_run/shift_00_00/config_extra.fdf; then
    echo -e "   -> ${RED}Failed:${NC} --mode 1 must NOT force MD.Steps to 0 (SIESTA needs to relax)"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} MD.Steps is NOT forced to 0 under --mode 1"
    PASS=$((PASS+1))
fi

# --- 3f. Regression: omitting --mode entirely defaults to mode 1 (SIESTA
#     relaxes z), not the old mode-3 (fixed gap/single-point) default ---
echo -e "\n--- Testing that --mode defaults to 1 when omitted ---"
rm -rf sf_run
stb-stackingfault -l1 graphene.fdf -l2 graphene.fdf -c calc.fdf -nx 2 -ny 2 --no-intro \
    > log_defaultmode.txt 2>&1
check_exit_code $? 0
check_contains "Mode            : 1" log_defaultmode.txt
check_contains "MD.TypeOfRun.*CG" sf_run/shift_00_00/config_extra.fdf
check_contains "MD.Steps.*200" sf_run/shift_00_00/config_extra.fdf
check_contains "position from 1 to .* 1.0 0.0 0.0" sf_run/shift_00_00/config_extra.fdf


# --- 4. Grid resolution minimum ---
echo -e "\n--- Testing -nx/-ny minimum (nx=1) ---"
stb-stackingfault -l1 graphene.fdf -l2 graphene.fdf -c calc.fdf -nx 1 -ny 3 --no-intro > log_nmin.txt 2>&1
check_exit_code $? 2

echo -e "\n--- Testing -nx/-ny minimum (ny=1) ---"
stb-stackingfault -l1 graphene.fdf -l2 graphene.fdf -c calc.fdf -nx 3 -ny 1 --no-intro > log_nmin2.txt 2>&1
check_exit_code $? 2


# --- 5. Error cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: missing layer1 file"
stb-stackingfault -l1 does_not_exist.fdf -l2 graphene.fdf -c calc.fdf --no-intro > log_missing.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing.txt

echo "Testing: --ml-custom-model with a nonexistent file"
stb-stackingfault -l1 graphene.fdf -l2 graphene.fdf -c calc.fdf --mode 2 \
    --ml-custom-model does_not_exist.model --no-intro > log_custommodel.txt 2>&1
check_exit_code $? 2
check_contains "\-\-ml-custom-model file not found" log_custommodel.txt

echo "Testing: --version"
stb-stackingfault --version > log_version.txt 2>&1
check_contains "stb-stackingfault" log_version.txt

echo "Testing: --help documents --grid-nx/--grid-ny/--gap/--save-report"
stb-stackingfault --help > log_help.txt 2>&1
check_contains "grid-nx" log_help.txt
check_contains "grid-ny" log_help.txt
check_contains "gap" log_help.txt
check_contains "\-\-d3" log_help.txt
check_contains "\-\-mode" log_help.txt
check_contains "relax-z-steps" log_help.txt
check_contains "save-report" log_help.txt
check_contains "ml-prerelax-layers" log_help.txt
check_contains "ml-preview" log_help.txt


# --- 6. Interactive path (stb-suite, shortcut 4.10.1) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.10.1) ---"

echo "Testing: navigate 4.10.1 -> homobilayer defaults, asymmetric 4x2 grid, save report -> quit"
rm -rf sf_run
{
  echo "4.10.1"
  echo "graphene.fdf"  # layer1_file
  echo "graphene.fdf"  # layer2_file
  echo "calc.fdf"      # calc_file
  echo ""              # mode (default 1)
  echo ""              # relax_z_steps (default 200, mode 1 asks this)
  echo ""              # pp_path (skip)
  echo "4"             # grid_nx
  echo "2"             # grid_ny
  echo ""              # gap (default 3.2)
  echo ""              # d3_choice (default Y now)
  echo ""              # ml_prerelax_choice (default N)
  echo ""              # ml_preview_choice (default N)
  echo ""              # show_advanced (default -> skip)
  echo "y"             # save_report -> Y
  echo ""              # press enter to continue
  echo "0"             # quit stage submenu
} | stb-suite > log_menu.txt 2>&1
check_contains "Success:.*8 grid folder" log_menu.txt
check_success sf_run/shift_00_00/structure.fdf
check_success sf_run/shift_03_01/structure.fdf
check_success sf_run/stackingfault_setup.txt

echo "Testing: navigate 4.10.1 -> ML preview on -> quit"
rm -rf sf_run
{
  echo "4.10.1"
  echo "graphene.fdf"
  echo "graphene.fdf"
  echo "calc.fdf"
  echo ""              # mode (default 1)
  echo ""              # relax_z_steps (default 200, mode 1 asks this)
  echo ""
  echo "3"              # grid_nx
  echo "3"              # grid_ny
  echo ""
  echo ""              # d3_choice
  echo ""              # ml_prerelax_choice
  echo "y"              # ml_preview_choice -> Y
  echo ""              # show_advanced
  echo "n"              # save_report -> N
  echo ""              # press enter
  echo "0"
} | stb-suite > log_menu_preview.txt 2>&1
check_contains "Success:.*9 grid folder" log_menu_preview.txt
check_success sf_run/stackingfault_ml_preview.png

echo "Testing: navigate 4.10.1 -> mode 1 (SIESTA relaxes z), custom relax-z-steps -> quit"
rm -rf sf_run
{
  echo "4.10.1"
  echo "graphene.fdf"
  echo "graphene.fdf"
  echo "calc.fdf"
  echo "1"                # mode -> 1
  echo "40"               # relax_z_steps
  echo ""
  echo "2"              # grid_nx
  echo "2"              # grid_ny
  echo ""                # gap
  echo ""                # d3_choice
  echo ""                # ml_prerelax_choice
  echo ""                # ml_preview_choice
  echo ""                # show_advanced
  echo "n"                # save_report -> N
  echo ""
  echo "0"
} | stb-suite > log_menu_mode1.txt 2>&1
check_contains "Success:.*4 grid folder" log_menu_mode1.txt
check_contains "MD.Steps.*40" sf_run/shift_00_00/config_extra.fdf


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
