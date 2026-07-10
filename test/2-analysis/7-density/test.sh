#!/bin/bash

# --- Setup ---
# Smoke test for stb-density (Density Plotter, item 3.8)
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
echo "--- Starting tester for STB-Density (item 3.8) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/Sn3O4.RHO" "$TEST_DIR/"
cp "$FIXTURE_DIR/Sn3O4.RHO" "$TEST_DIR/Sn3O4_copy.RHO"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Basic 2D slice (default output name derived from --label) ---
echo -e "\n--- Testing a basic 2D slice ---"

timeout 60 stb-density -l Sn3O4 --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "Mode: 2D Map (Slice)" log_basic.txt
check_contains "Integrated Charge Density:" log_basic.txt
check_success Sn3O4_density.dat
check_success Sn3O4_density.gplot
check_contains "Charge Density (e/Ang^3)" Sn3O4_density.gplot
check_contains "White -> Yellow -> Red" Sn3O4_density.gplot
echo "Testing: a non-negative quantity gets a cbrange anchored at 0"
check_contains "set cbrange \[0:" Sn3O4_density.gplot


# --- 3. --pos and explicit --axis ---
echo -e "\n--- Testing --pos and --axis ---"

timeout 60 stb-density -l Sn3O4 --no-intro --axis 1 --pos 2.5 -o axis1.dat > log_axis1.txt 2>&1
check_exit_code $? 0
check_contains "Mapping plane perpendicular to Y axis" log_axis1.txt
check_success axis1.dat

echo "Testing: a skewed cut plane is flagged (axis 0 here isn't ~90deg between its in-plane vectors)"
timeout 60 stb-density -l Sn3O4 --no-intro --axis 0 -o axis0.dat > log_axis0.txt 2>&1
check_exit_code $? 0
check_contains "cut plane is skewed" log_axis0.txt
check_contains "cut plane is skewed" axis0.gplot

echo "Testing: --pos out of bounds is rejected"
stb-density -l Sn3O4 --no-intro --pos 9999 -o oob.dat > log_oob.txt 2>&1
check_exit_code $? 1
check_contains "Position out of bounds" log_oob.txt


# --- 4. 3D volume export + --iso-min filtering ---
echo -e "\n--- Testing --3d ---"

timeout 60 stb-density -l Sn3O4 --no-intro --3d -o full3d.dat > log_3d.txt 2>&1
check_exit_code $? 0
check_contains "Mode: Full 3D Volume Export" log_3d.txt
check_success full3d.dat
FULL_LINES=$(wc -l < full3d.dat)

echo "Testing: --iso-min reduces the point count"
timeout 60 stb-density -l Sn3O4 --no-intro --3d --iso-min 0.1 -o filtered3d.dat > log_iso.txt 2>&1
check_exit_code $? 0
check_contains "--iso-min 0.1: kept" log_iso.txt
check_success filtered3d.dat
FILTERED_LINES=$(wc -l < filtered3d.dat)
if [ "$FILTERED_LINES" -lt "$FULL_LINES" ]; then
    echo -e "   -> ${GREEN}Verified:${NC} filtered3d.dat ($FILTERED_LINES lines) is smaller than full3d.dat ($FULL_LINES lines)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} --iso-min did not reduce the point count"
    FAIL=$((FAIL+1))
fi

echo "Testing: --iso-min above every value in the grid is a clean error"
stb-density -l Sn3O4 --no-intro --3d --iso-min 999 -o empty3d.dat > log_iso_empty.txt 2>&1
check_exit_code $? 1
check_contains "No points survived --iso-min" log_iso_empty.txt


# --- 5. --spin ---
echo -e "\n--- Testing --spin on a non-spin-polarized .RHO (clean error) ---"

stb-density -l Sn3O4 --no-intro --spin -o spin_test.dat > log_spin.txt 2>&1
check_exit_code $? 1
check_contains "no spin component" log_spin.txt


# --- 6. --rho2 (charge density difference) ---
echo -e "\n--- Testing --rho2 (difference mode) ---"

timeout 60 stb-density -l Sn3O4 --no-intro --rho2 Sn3O4_copy.RHO -o diff.dat > log_diff.txt 2>&1
check_exit_code $? 0
check_contains "Delta Charge Density" diff.gplot
check_contains "Blue -> White -> Red" diff.gplot
python3 -c "
import numpy as np
import sys
d = np.loadtxt('diff.dat')
sys.exit(0 if np.abs(d[:, 3]).max() < 1e-9 else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} identical files subtract to exactly zero everywhere"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} difference of identical files was not zero"
    FAIL=$((FAIL+1))
fi

echo "Testing: --rho2 with a mismatched grid shape is rejected"
python3 -c "
import sisl
g = sisl.get_sile('Sn3O4.RHO').read_grid()
small = sisl.Grid([4, 4, 4], lattice=g.lattice)
small.grid[:] = 0.0
small.write('mismatched.RHO')
"
stb-density -l Sn3O4 --no-intro --rho2 mismatched.RHO -o baddiff.dat > log_shapemismatch.txt 2>&1
check_exit_code $? 1
check_contains "shape mismatch" log_shapemismatch.txt

echo "Testing: a signed quantity (--rho2 diff) gets a colorbar range symmetric around zero"
python3 -c "
import sisl
sile = sisl.get_sile('Sn3O4.RHO')
g = sile.read_grid()
g.grid = g.grid * 0.5
g.write('Sn3O4_scaled.RHO')
"
timeout 60 stb-density -l Sn3O4 --no-intro --rho2 Sn3O4_scaled.RHO -o diffscale.dat > log_diffscale.txt 2>&1
check_exit_code $? 0
python3 -c "
import re, sys
txt = open('diffscale.gplot').read()
m = re.search(r'set cbrange \[(-?[0-9.eE+-]+):(-?[0-9.eE+-]+)\]', txt)
sys.exit(0 if m and abs(float(m.group(1)) + float(m.group(2))) < 1e-9 else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} cbrange is symmetric around zero for a signed quantity"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} cbrange for a signed quantity is not symmetric around zero"
    FAIL=$((FAIL+1))
fi


# --- 7. --profile (planar-averaged 1D profile) ---
echo -e "\n--- Testing --profile ---"

timeout 60 stb-density -l Sn3O4 --no-intro --profile --axis 2 -o profile.dat > log_profile.txt 2>&1
check_exit_code $? 0
check_contains "Mode: Planar-Averaged Profile" log_profile.txt
check_success profile.dat
check_success profile.gplot
check_contains "Planar-averaged Charge Density along Z" profile.gplot
check_contains "with lines" profile.gplot

echo "Testing: --pos has no effect (and warns) in --profile mode"
stb-density -l Sn3O4 --no-intro --profile --pos 1.0 -o profilepos.dat > log_profilepos.txt 2>&1
check_exit_code $? 0
check_contains "pos has no effect in --profile mode" log_profilepos.txt

echo "Testing: --3d and --profile are mutually exclusive"
stb-density -l Sn3O4 --no-intro --3d --profile -o bad.dat > log_mutex.txt 2>&1
check_exit_code $? 2
check_contains "not allowed with argument" log_mutex.txt


# --- 8. --cube (Gaussian cube export) ---
echo -e "\n--- Testing --cube ---"

echo "Testing: --cube without a geometry file (.XV/.fdf) is a clean error"
stb-density -l Sn3O4 --no-intro --cube -o nogeo.dat > log_cube_nogeo.txt 2>&1
check_exit_code $? 1
check_contains "requires a geometry file" log_cube_nogeo.txt

echo "Testing: --cube with a matching geometry file succeeds"
python3 -c "
import sisl
g = sisl.get_sile('Sn3O4.RHO').read_grid()
geom = sisl.Geometry([[0, 0, 0], [1, 1, 1]], atoms=[sisl.Atom(50), sisl.Atom(8)], lattice=g.lattice.cell)
geom.write('Sn3O4.fdf')
"
timeout 60 stb-density -l Sn3O4 --no-intro --cube -o withcube.dat > log_cube_ok.txt 2>&1
check_exit_code $? 0
check_contains "Cube file saved" log_cube_ok.txt
check_success withcube.cube

echo "Testing: --cube with a mismatched geometry cell is rejected"
cp Sn3O4.RHO Mismatch.RHO
python3 -c "
import sisl
geom = sisl.Geometry([[0, 0, 0]], atoms=[sisl.Atom(14)], lattice=[[5.0, 0, 0], [0, 5.0, 0], [0, 0, 5.0]])
geom.write('Mismatch.fdf')
"
stb-density -l Mismatch --no-intro --cube -o mismatch_cube.dat > log_cube_mismatch.txt 2>&1
check_exit_code $? 1
check_contains "does not match the density grid" log_cube_mismatch.txt


# --- 9. --vmin/--vmax (manual colorbar range) ---
echo -e "\n--- Testing --vmin/--vmax ---"

timeout 60 stb-density -l Sn3O4 --no-intro --vmin -0.5 --vmax 2.0 -o vrange.dat > log_vrange.txt 2>&1
check_exit_code $? 0
check_contains "set cbrange \[-0.5:2\]" vrange.gplot


# --- 10. --contour ---
echo -e "\n--- Testing --contour ---"

timeout 60 stb-density -l Sn3O4 --no-intro --contour -o contourtest.dat > log_contour.txt 2>&1
check_exit_code $? 0
check_contains "set contour base" contourtest.gplot
check_contains "set cntrparam levels" contourtest.gplot

echo "Testing: --contour is ignored (with a warning) outside slice mode"
stb-density -l Sn3O4 --no-intro --3d --contour -o contour3d.dat > log_contour3d.txt 2>&1
check_exit_code $? 0
check_contains "contour only applies to the default 2D slice mode" log_contour3d.txt


# --- 11. Error and robustness cases ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: missing .RHO file"
stb-density -l does_not_exist --no-intro > log_missing.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing.txt

echo "Testing: --label omitted (required)"
stb-density --no-intro > log_missing_label.txt 2>&1
check_exit_code $? 2

echo "Testing: --axis outside 0-2 is rejected"
stb-density -l Sn3O4 --axis 5 --no-intro > log_badaxis.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-density --version > log_version.txt 2>&1
check_contains "stb-density" log_version.txt

echo "Testing: --help documents --label, --3d, --profile, --axis, --spin, --rho2, --iso-min, --cube, --vmin, --vmax, --contour"
stb-density --help > log_help.txt 2>&1
check_contains "label" log_help.txt
check_contains "3d" log_help.txt
check_contains "profile" log_help.txt
check_contains "axis" log_help.txt
check_contains "spin" log_help.txt
check_contains "rho2" log_help.txt
check_contains "iso-min" log_help.txt
check_contains "cube" log_help.txt
check_contains "vmin" log_help.txt
check_contains "vmax" log_help.txt
check_contains "contour" log_help.txt


# --- 12. Interactive path (stb-suite, shortcut 2.8) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 2.8) ---"

echo "Testing: navigate 2.8 -> label Sn3O4 -> mode 1 (2D) -> axis 2 -> center -> no contour -> no spin -> no rho2 -> no manual range -> Enter -> quit"
printf '3.8\nSn3O4\n1\n2\n\nn\nn\n\nn\n\n0\n' | timeout 60 stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_contains "Mode: 2D Map (Slice)" log_menu.txt
check_success Sn3O4_density.dat

echo "Testing: navigate 2.8 -> profile mode (3) -> axis 2 -> no spin -> no rho2 -> Enter -> quit"
printf '3.8\nSn3O4\n3\n2\nn\n\n\n0\n' | timeout 60 stb-suite > log_menu_profile.txt 2>&1
check_exit_code $? 0
check_contains "Mode: Planar-Averaged Profile" log_menu_profile.txt


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

[ "$FAIL" -eq 0 ]
