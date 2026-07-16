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


# --- 2. Default homobilayer grid (graphene/graphene, -n 3) ---
echo -e "\n--- Testing homobilayer grid (graphene/graphene, -n 3) ---"
stb-stackingfault -l1 graphene.fdf -l2 graphene.fdf -c calc.fdf -n 3 --no-intro > log_homo.txt 2>&1
check_exit_code $? 0
check_success shift_00_00/structure.fdf
check_success shift_00_00/calc.fdf
check_success shift_02_02/structure.fdf
check_contains "MD.TypeOfRun.*CG" shift_01_01/calc.fdf
check_contains "MD.NumCGsteps.*0" shift_01_01/calc.fdf
check_success stackingfault_setup.txt
check_contains "\[0\] RUN METADATA" stackingfault_setup.txt
check_contains "\[1\] ZSL MATCH" stackingfault_setup.txt
check_contains "Each grid point: 4 atoms (2 layer 1 + 2 layer 2)" stackingfault_setup.txt
check_contains "9 grid point(s) total" stackingfault_setup.txt
check_contains "\[2\] GRID FOLDERS" stackingfault_setup.txt
check_contains "\[3\] SUMMARY" stackingfault_setup.txt
check_contains "GRID_TABLE" stackingfault_setup.txt

echo "Testing: GRID_TABLE has exactly 9 rows (3x3) and shift AA registry (0,0) is present"
python3 -c "
rows = []
in_table = False
with open('stackingfault_setup.txt') as f:
    for line in f:
        if line.startswith('# GRID_TABLE'):
            in_table = True
            continue
        if not in_table:
            continue
        if line.startswith('#') or not line.strip():
            continue
        rows.append(line.split())
assert len(rows) == 9, f'expected 9 grid rows, got {len(rows)}'
labels = {r[0] for r in rows}
assert 'shift_00_00' in labels
print('OK')
" > log_table_check.txt 2>&1
check_contains "OK" log_table_check.txt

echo "Testing: every grid folder has the same atom count (only layer-2 position changes)"
python3 -c "
from stb.core import structure_io
import glob
counts = set()
for d in sorted(glob.glob('shift_*')):
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
s00 = structure_io.read_fdf('shift_00_00/structure.fdf')
s11 = structure_io.read_fdf('shift_01_01/structure.fdf')
# layer 2 = last 2 atoms (4-atom cell, 2 per layer)
l2_00 = np.array([pos for _, pos in s00.atoms[2:]])
l2_11 = np.array([pos for _, pos in s11.atoms[2:]])
delta = (l2_11 - l2_00) % 1.0
expected = np.array([1/3, 1/3, 0.0])
assert np.allclose(delta, expected, atol=1e-4), f'expected shift {expected}, got {delta}'
print('OK')
" > log_shift_check.txt 2>&1
check_contains "OK" log_shift_check.txt


# --- 3. Heterostructure grid (graphene/hBN) ---
echo -e "\n--- Testing heterostructure grid (graphene/hBN) ---"
rm -rf shift_* stackingfault_setup.txt
stb-stackingfault -l1 graphene.fdf -l2 hbn.fdf -c calc.fdf -n 3 --no-intro > log_hetero.txt 2>&1
check_exit_code $? 0
check_success shift_00_00/structure.fdf
check_contains " 1   6   C" shift_00_00/structure.fdf
check_contains " 2   5   B" shift_00_00/structure.fdf
check_contains " 3   7   N" shift_00_00/structure.fdf


# --- 4. Grid resolution minimum ---
echo -e "\n--- Testing -n minimum (n=1) ---"
stb-stackingfault -l1 graphene.fdf -l2 graphene.fdf -c calc.fdf -n 1 --no-intro > log_nmin.txt 2>&1
check_exit_code $? 2


# --- 5. Error cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: missing layer1 file"
stb-stackingfault -l1 does_not_exist.fdf -l2 graphene.fdf -c calc.fdf --no-intro > log_missing.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing.txt

echo "Testing: --version"
stb-stackingfault --version > log_version.txt 2>&1
check_contains "stb-stackingfault" log_version.txt

echo "Testing: --help documents --grid-n/--gap"
stb-stackingfault --help > log_help.txt 2>&1
check_contains "grid-n" log_help.txt
check_contains "gap" log_help.txt


# --- 6. Interactive path (stb-suite, shortcut 4.10.1) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.10.1) ---"

echo "Testing: navigate 4.10.1 -> homobilayer defaults -> quit"
rm -rf shift_* stackingfault_setup.txt
{
  echo "4.10.1"
  echo "graphene.fdf"  # layer1_file
  echo "graphene.fdf"  # layer2_file
  echo "calc.fdf"      # calc_file
  echo ""              # pp_path (skip)
  echo "3"             # grid_n
  echo ""              # gap (default 3.2)
  echo ""              # show_advanced (default -> skip)
  echo ""              # press enter to continue
  echo "0"             # quit stage submenu
} | stb-suite > log_menu.txt 2>&1
check_contains "Success:.*9 grid folder" log_menu.txt
check_success shift_00_00/structure.fdf
check_success shift_02_02/structure.fdf


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
