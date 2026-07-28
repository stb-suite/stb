#!/bin/bash

# --- Setup ---
# Smoke test for stb-density (Density Plotter, item 3.8)
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$FIXTURE_DIR/test_files"

# Non-interactive matplotlib backend: --view calls plt.show(), which would
# otherwise block waiting for a window to close.
export MPLBACKEND=Agg

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

check_absent() {
    if [ ! -e "$1" ]; then
        echo -e " ... ${GREEN}OK${NC} (file '$1' absent, as expected)"
        PASS=$((PASS+1))
    else
        echo -e " ... ${RED}FAIL${NC} (file '$1' should not exist)"
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

echo "Generating a synthetic spin-polarized .RHO fixture (nspin=2: total, spin)..."
python3 << 'PYEOF'
import sisl
import numpy as np

# A 2-component (nspin=2) .RHO stores the raw UP and DOWN spin channels
# directly (SIESTA's own collinear-polarized convention, confirmed against
# a real spin-polarized O2 calculation in test/6-utils/3-cube/o2.RHO) --
# NOT "total" and "spin" as separate components. up > down everywhere here,
# by construction, so the derived total (up+down) and net spin (up-down)
# are both well-defined, non-degenerate, physically sensible quantities.
lattice = sisl.Lattice([[4.0, 0, 0], [0, 4.0, 0], [0, 0, 10.0]])
nx, ny, nz = 8, 8, 20
rng = np.random.default_rng(1)
up = np.abs(rng.normal(0.6, 0.15, (nx, ny, nz)))
down = np.abs(rng.normal(0.4, 0.1, (nx, ny, nz)))

g_up = sisl.Grid([nx, ny, nz], lattice=lattice)
g_up.grid[:] = up
g_down = sisl.Grid([nx, ny, nz], lattice=lattice)
g_down.grid[:] = down

sisl.get_sile("spintest.RHO", mode="w").write_grid(g_up, g_down)
print("spintest.RHO written.")
PYEOF


# --- 2. Basic 2D slice ---
echo -e "\n--- Testing a basic 2D slice ---"

timeout 60 stb-density -l Sn3O4 --no-intro --output-dir out_basic --save-gnuplot > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "\[0\] RUN METADATA" log_basic.txt
check_contains "\[1\] CHARGE DENSITY" log_basic.txt
check_contains "\[4\] REFERENCES" log_basic.txt
check_contains "\[5\] SUMMARY & FILES" log_basic.txt
check_contains "Spin-polarized : no" log_basic.txt
check_contains "Mode: 2D Map (Slice)" log_basic.txt
check_contains "Integrated Charge Density:" log_basic.txt
check_success out_basic/Sn3O4_density.dat
check_success out_basic/Sn3O4_density.gplot
check_success out_basic/references.bib
check_contains "@article{Soler2002," out_basic/references.bib
check_contains "Charge Density (e/Ang^3)" out_basic/Sn3O4_density.gplot
check_contains "White -> Yellow -> Red" out_basic/Sn3O4_density.gplot
echo "Testing: a non-negative quantity gets a cbrange anchored at 0"
check_contains "set cbrange \[0:" out_basic/Sn3O4_density.gplot

echo "Testing: a plain run (no --save-gnuplot) writes only the .dat file + references.bib"
timeout 60 stb-density -l Sn3O4 --no-intro --output-dir out_plain > log_plain.txt 2>&1
check_exit_code $? 0
check_success out_plain/Sn3O4_density.dat
check_success out_plain/references.bib
check_absent out_plain/Sn3O4_density.gplot


# --- 3. --pos and explicit --axis ---
echo -e "\n--- Testing --pos and --axis ---"

timeout 60 stb-density -l Sn3O4 --no-intro --axis 1 --pos 2.5 --output-dir out_axis1 > log_axis1.txt 2>&1
check_exit_code $? 0
check_contains "Mapping plane perpendicular to Y axis" log_axis1.txt
check_success out_axis1/Sn3O4_density.dat

echo "Testing: a skewed cut plane is flagged (axis 0 here isn't ~90deg between its in-plane vectors)"
timeout 60 stb-density -l Sn3O4 --no-intro --axis 0 --output-dir out_axis0 --save-gnuplot > log_axis0.txt 2>&1
check_exit_code $? 0
check_contains "cut plane is skewed" log_axis0.txt
check_contains "cut plane is skewed" out_axis0/Sn3O4_density.gplot

echo "Testing: --pos out of bounds is rejected"
stb-density -l Sn3O4 --no-intro --pos 9999 --output-dir out_oob > log_oob.txt 2>&1
check_exit_code $? 1
check_contains "Position out of bounds" log_oob.txt


# --- 4. 3D volume export + --iso-min filtering ---
echo -e "\n--- Testing --3d ---"

timeout 60 stb-density -l Sn3O4 --no-intro --3d --output-dir out_full3d > log_3d.txt 2>&1
check_exit_code $? 0
check_contains "Mode: Full 3D Volume Export" log_3d.txt
check_success out_full3d/Sn3O4_density.dat
FULL_LINES=$(wc -l < out_full3d/Sn3O4_density.dat)

echo "Testing: --iso-min reduces the point count"
timeout 60 stb-density -l Sn3O4 --no-intro --3d --iso-min 0.1 --output-dir out_filtered3d > log_iso.txt 2>&1
check_exit_code $? 0
check_contains "--iso-min 0.1: kept" log_iso.txt
check_success out_filtered3d/Sn3O4_density.dat
FILTERED_LINES=$(wc -l < out_filtered3d/Sn3O4_density.dat)
if [ "$FILTERED_LINES" -lt "$FULL_LINES" ]; then
    echo -e "   -> ${GREEN}Verified:${NC} filtered (${FILTERED_LINES} lines) is smaller than full (${FULL_LINES} lines)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} --iso-min did not reduce the point count"
    FAIL=$((FAIL+1))
fi

echo "Testing: --iso-min above every value in the grid is a clean error"
stb-density -l Sn3O4 --no-intro --3d --iso-min 999 --output-dir out_empty3d > log_iso_empty.txt 2>&1
check_exit_code $? 1
check_contains "No points survived --iso-min" log_iso_empty.txt


# --- 5. --spin ---
echo -e "\n--- Testing --spin on a non-spin-polarized .RHO (clean error) ---"

stb-density -l Sn3O4 --no-intro --spin --output-dir out_spinerr > log_spin.txt 2>&1
check_exit_code $? 1
check_contains "no spin component" log_spin.txt

echo "Testing: a spin-polarized .RHO is auto-detected -- no --spin flag needed to see it"
timeout 60 stb-density -l spintest --no-intro --output-dir out_spinauto --save-gnuplot > log_spinauto.txt 2>&1
check_exit_code $? 0
check_contains "Spin-polarized : yes (auto-detected)" log_spinauto.txt
check_contains "Spin-polarized .RHO detected -- the net spin (magnetization) density will be" log_spinauto.txt
check_contains "\[1\] CHARGE DENSITY" log_spinauto.txt
check_contains "\[2\] SPIN DENSITY" log_spinauto.txt
check_success out_spinauto/spintest_density.dat
check_success out_spinauto/spintest_density.gplot
check_success out_spinauto/spintest_density_spin.dat
check_success out_spinauto/spintest_density_spin.gplot

echo "Testing: --spin on a spin-polarized file processes ONLY spin, skipping charge"
timeout 60 stb-density -l spintest --no-intro --spin --output-dir out_spinonly > log_spinonly.txt 2>&1
check_exit_code $? 0
check_contains "Skipped (--spin: only the spin density is processed this run)." log_spinonly.txt
check_contains "\[2\] SPIN DENSITY" log_spinonly.txt


# --- 6. --rho2 (charge density difference) ---
echo -e "\n--- Testing --rho2 (difference mode) ---"

timeout 60 stb-density -l Sn3O4 --no-intro --rho2 Sn3O4_copy.RHO --output-dir out_diff --save-gnuplot > log_diff.txt 2>&1
check_exit_code $? 0
check_contains "Delta Charge Density" out_diff/Sn3O4_density.gplot
check_contains "Blue -> White -> Red" out_diff/Sn3O4_density.gplot
python3 -c "
import numpy as np
import sys
d = np.loadtxt('out_diff/Sn3O4_density.dat')
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
stb-density -l Sn3O4 --no-intro --rho2 mismatched.RHO --output-dir out_baddiff > log_shapemismatch.txt 2>&1
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
timeout 60 stb-density -l Sn3O4 --no-intro --rho2 Sn3O4_scaled.RHO --output-dir out_diffscale --save-gnuplot > log_diffscale.txt 2>&1
check_exit_code $? 0
python3 -c "
import re, sys
txt = open('out_diffscale/Sn3O4_density.gplot').read()
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

timeout 60 stb-density -l Sn3O4 --no-intro --profile --axis 2 --output-dir out_profile --save-gnuplot > log_profile.txt 2>&1
check_exit_code $? 0
check_contains "Mode: Planar-Averaged Profile" log_profile.txt
check_success out_profile/Sn3O4_density.dat
check_success out_profile/Sn3O4_density.gplot
check_contains "Planar-averaged Charge Density along Z" out_profile/Sn3O4_density.gplot
check_contains "with lines" out_profile/Sn3O4_density.gplot

echo "Testing: --pos has no effect (and warns) in --profile mode"
stb-density -l Sn3O4 --no-intro --profile --pos 1.0 --output-dir out_profilepos > log_profilepos.txt 2>&1
check_exit_code $? 0
check_contains "pos has no effect in --profile mode" log_profilepos.txt

echo "Testing: --3d and --profile are mutually exclusive"
stb-density -l Sn3O4 --no-intro --3d --profile --output-dir out_bad > log_mutex.txt 2>&1
check_exit_code $? 2
check_contains "not allowed with argument" log_mutex.txt


# --- 8. --cube (Gaussian cube export) ---
echo -e "\n--- Testing --cube ---"

echo "Testing: --cube without a geometry file (.XV/.fdf) degrades gracefully (no cube, run still succeeds)"
stb-density -l Sn3O4 --no-intro --cube --output-dir out_nogeo > log_cube_nogeo.txt 2>&1
check_exit_code $? 0
check_contains "requires a geometry file" log_cube_nogeo.txt
check_absent out_nogeo/Sn3O4_density.cube

echo "Testing: --cube with a matching geometry file succeeds"
python3 -c "
import sisl
g = sisl.get_sile('Sn3O4.RHO').read_grid()
geom = sisl.Geometry([[0, 0, 0], [1, 1, 1]], atoms=[sisl.Atom(50), sisl.Atom(8)], lattice=g.lattice.cell)
geom.write('Sn3O4.fdf')
"
timeout 60 stb-density -l Sn3O4 --no-intro --cube --output-dir out_withcube > log_cube_ok.txt 2>&1
check_exit_code $? 0
check_contains "Cube file saved" log_cube_ok.txt
check_success out_withcube/Sn3O4_density.cube

echo "Testing: --cube with a mismatched geometry cell degrades gracefully (no cube, run still succeeds)"
cp Sn3O4.RHO Mismatch.RHO
python3 -c "
import sisl
geom = sisl.Geometry([[0, 0, 0]], atoms=[sisl.Atom(14)], lattice=[[5.0, 0, 0], [0, 5.0, 0], [0, 0, 5.0]])
geom.write('Mismatch.fdf')
"
stb-density -l Mismatch --no-intro --cube --output-dir out_cubemismatch > log_cube_mismatch.txt 2>&1
check_exit_code $? 0
check_contains "does not match the density grid" log_cube_mismatch.txt
check_absent out_cubemismatch/Mismatch_density.cube


# --- 9. --vmin/--vmax (manual colorbar range) ---
echo -e "\n--- Testing --vmin/--vmax ---"

timeout 60 stb-density -l Sn3O4 --no-intro --vmin -0.5 --vmax 2.0 --output-dir out_vrange --save-gnuplot > log_vrange.txt 2>&1
check_exit_code $? 0
check_contains "set cbrange \[-0.5:2\]" out_vrange/Sn3O4_density.gplot


# --- 10. --contour ---
echo -e "\n--- Testing --contour ---"

timeout 60 stb-density -l Sn3O4 --no-intro --contour --output-dir out_contour --save-gnuplot > log_contour.txt 2>&1
check_exit_code $? 0
check_contains "set contour base" out_contour/Sn3O4_density.gplot
check_contains "set cntrparam levels" out_contour/Sn3O4_density.gplot

echo "Testing: --contour is ignored (with a warning) outside slice mode"
stb-density -l Sn3O4 --no-intro --3d --contour --output-dir out_contour3d > log_contour3d.txt 2>&1
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

echo "Testing: --help documents --label, --output-dir, --3d, --profile, --axis, --spin, --rho2, --iso-min, --cube, --vmin, --vmax, --contour, --save-report, --save-gnuplot, --view"
stb-density --help > log_help.txt 2>&1
check_contains "label" log_help.txt
check_contains "output-dir" log_help.txt
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
check_contains "save-report" log_help.txt
check_contains "save-gnuplot" log_help.txt
check_contains "\-\-view" log_help.txt


# --- 12. Interactive path (stb-suite, shortcut 3.8) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 3.8) ---"

echo "Testing: navigate 3.8 -> label Sn3O4 -> mode 1 (2D) -> axis 2 -> center -> no contour -> no spin-only -> no rho2 -> no manual range -> default output dir -> no save-report/gnuplot/view -> quit"
printf '3.8\nSn3O4\n1\n2\n\nn\nn\n\nn\n\nn\nn\nn\n\n0\n' | timeout 60 stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_contains "Mode: 2D Map (Slice)" log_menu.txt
check_success Sn3O4_density.dat

echo "Testing: navigate 3.8 -> profile mode (3) -> axis 2 -> no spin-only -> no rho2 -> no manual range n/a -> default output dir -> no save-report/gnuplot/view -> quit"
rm -f Sn3O4_density.dat
printf '3.8\nSn3O4\n3\n2\nn\n\n\nn\nn\nn\n\n0\n' | timeout 60 stb-suite > log_menu_profile.txt 2>&1
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
