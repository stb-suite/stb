#!/bin/bash

# --- Setup ---
# Smoke test for stb-workfunction (Work Function Calculator, item 2.7)
# NOTE: no real SIESTA .VT/.XV/.out fixtures exist for this tool, so this test
# generates small synthetic-but-realistic ones itself (same approach
# test_translate/test.sh uses for structure files) instead of committing
# binary grid fixtures to the repo.
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$FIXTURE_DIR/test_files"

# Non-interactive matplotlib backend: the interactive path (no --no-plot toggle exposed
# there) calls plt.show(), which would otherwise block waiting for a window to close.
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
echo "--- Starting tester for STB-WorkFunction (item 2.7) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null

echo "Generating synthetic .VT/.XV/.out fixtures..."
python3 << 'PYEOF'
import sisl
import numpy as np

lattice = sisl.Lattice([[4.0, 0, 0], [0, 4.0, 0], [0, 0, 30.0]])
Nx, Ny, Nz = 10, 10, 300
z = np.linspace(0, 30.0, Nz, endpoint=False)

def make_vt(fname, region_fn):
    data = np.zeros((Nx, Ny, Nz))
    for k in range(Nz):
        data[:, :, k] = region_fn(z[k])
    grid = sisl.Grid([Nx, Ny, Nz], lattice=lattice)
    grid.grid[:] = data
    grid.write(fname)

# Symmetric slab: one flat vacuum level (5.0 eV) wrapping around the boundary.
make_vt("symmetric.VT", lambda zz: 2.0 + 3.0*np.sin(2*np.pi*10*(zz-12.0)/6.0) if 12.0 <= zz <= 18.0 else 5.0)

# Asymmetric slab: two distinct vacuum levels (5.0 below, 6.5 above).
def asym(zz):
    if 13.0 <= zz <= 17.0:
        return 2.0 + 3.0*np.sin(2*np.pi*10*(zz-13.0)/4.0)
    return 5.0 if zz < 13.0 else 6.5
make_vt("asymmetric.VT", asym)

# Bulk: no real vacuum anywhere along z (pure oscillation).
rng = np.random.default_rng(0)
make_vt("bulk.VT", lambda zz: 2.0 + 3.0*np.sin(2*np.pi*10*zz/30.0) + 0.01*rng.standard_normal())

# Matching geometry for the symmetric case, atoms only in z=[12,18] -- used to
# test --axis auto-detection (z is the only vacuum-padded direction here).
atoms = [sisl.Atom(6)] * 6
xyz = np.array([[2.0, 2.0, zz] for zz in np.linspace(12.5, 17.5, 6)])
sisl.Geometry(xyz, atoms=atoms, lattice=lattice).write("symmetric.fdf")

with open("symmetric.out", "w") as f:
    f.write("siesta:         Fermi energy:    -4.500000 eV\n")

# A geometry file with a totally different cell, paired (same label) with the real
# symmetric.VT -- used to test the geometry/grid cell cross-check.
mismatch_lattice = sisl.Lattice([[8.0, 0, 0], [0, 3.0, 0], [0, 0, 3.0]])
mismatch_atoms = [sisl.Atom(14)] * 3
mismatch_xyz = np.array([[x, 1.5, 1.5] for x in [1.0, 4.0, 7.0]])
sisl.Geometry(mismatch_xyz, atoms=mismatch_atoms, lattice=mismatch_lattice).write("mismatch.fdf")
make_vt("mismatch.VT", lambda zz: 2.0 + 3.0*np.sin(2*np.pi*10*(zz-12.0)/6.0) if 12.0 <= zz <= 18.0 else 5.0)

# Same cell/vacuum as the symmetric case, but with an extra "adsorbate" atom sitting
# inside the vacuum plateau (z=25) -- used to test the atom-in-vacuum-plateau warning.
withads_atoms = [sisl.Atom(6)] * 6 + [sisl.Atom(1)]
withads_xyz = np.array([[2.0, 2.0, zz] for zz in np.linspace(12.5, 17.5, 6)] + [[2.0, 2.0, 25.0]])
sisl.Geometry(withads_xyz, atoms=withads_atoms, lattice=lattice).write("withads.fdf")
make_vt("withads.VT", lambda zz: 2.0 + 3.0*np.sin(2*np.pi*10*(zz-12.0)/6.0) if 12.0 <= zz <= 18.0 else 5.0)

print("Fixtures generated.")
PYEOF


# --- 2. Basic run: single vacuum plateau ---
echo -e "\n--- Testing a basic run (symmetric slab, single vacuum plateau) ---"

timeout 60 stb-workfunction -l symmetric --grid symmetric.VT --fermi 0.0 --no-plot --no-intro \
    > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "Vacuum Level   =   5.0000 eV" log_basic.txt
check_contains "Work Function  =   5.0000 eV" log_basic.txt
check_success workfunction_data.dat
check_success workfunction.gplot
check_contains "Vacuum region 1:" workfunction_data.dat


# --- 3. Fermi energy auto-detected from .out ---
echo -e "\n--- Testing Fermi energy auto-detected from symmetric.out ---"

timeout 60 stb-workfunction -l symmetric --grid symmetric.VT --no-plot --no-intro \
    > log_fermi.txt 2>&1
check_exit_code $? 0
check_contains "Fermi Energy (Ef): -4.500000 eV" log_fermi.txt


# --- 4. --axis auto-detection from symmetric.fdf ---
echo -e "\n--- Testing --axis auto-detection ---"

timeout 60 stb-workfunction -l symmetric --grid symmetric.VT --fermi 0.0 --no-plot --no-intro \
    > log_autoaxis.txt 2>&1
check_exit_code $? 0
check_contains "Auto-detected vacuum axis: 2" log_autoaxis.txt

echo "Testing: explicit wrong --axis (0) warns instead of silently proceeding"
timeout 60 stb-workfunction -l symmetric --grid symmetric.VT --fermi 0.0 --axis 0 --no-plot --no-intro \
    > log_wrongaxis.txt 2>&1
check_exit_code $? 0
check_contains "doesn't look vacuum-padded" log_wrongaxis.txt

echo "Testing: a geometry file with an unrelated cell is ignored, not trusted"
timeout 60 stb-workfunction -l mismatch --grid mismatch.VT --fermi 0.0 --no-plot --no-intro \
    > log_mismatch.txt 2>&1
check_exit_code $? 0
check_contains "cell doesn't match" log_mismatch.txt

echo "Testing: an atom projecting into the vacuum plateau is flagged"
timeout 60 stb-workfunction -l withads --grid withads.VT --fermi 0.0 --no-plot --no-intro \
    > log_withads.txt 2>&1
check_exit_code $? 0
check_contains "atom(s) project into the vacuum region" log_withads.txt


# --- 5. Asymmetric slab: two distinct work functions ---
echo -e "\n--- Testing an asymmetric slab (two distinct vacuum levels) ---"

timeout 60 stb-workfunction -l asymmetric --grid asymmetric.VT --fermi 0.0 --no-plot --no-intro \
    > log_asym.txt 2>&1
check_exit_code $? 0
check_contains "2 distinct vacuum levels detected" log_asym.txt
check_contains "Vacuum Level 1 (near z=" log_asym.txt
check_contains "Work Function 1 " log_asym.txt
check_contains "Vacuum Level 2 (near z=" log_asym.txt
check_contains "Work Function 2 " log_asym.txt
check_contains "Average Work Function" log_asym.txt

echo "Testing: --asymmetric-tol wide enough merges the two levels into one"
timeout 60 stb-workfunction -l asymmetric --grid asymmetric.VT --fermi 0.0 --asymmetric-tol 5.0 \
    --no-plot --no-intro > log_wide_tol.txt 2>&1
check_exit_code $? 0
check_contains "Vacuum Level   =" log_wide_tol.txt

echo "Testing: --asymmetric-tol must be positive"
stb-workfunction -l asymmetric --asymmetric-tol -1 --no-intro > log_bad_tol.txt 2>&1
check_exit_code $? 2
check_contains "must be positive" log_bad_tol.txt


# --- 6. Bulk (no vacuum): clean failure instead of a bogus number ---
echo -e "\n--- Testing a bulk grid (no vacuum along the chosen axis) ---"

timeout 60 stb-workfunction -l bulk --grid bulk.VT --fermi 0.0 --no-plot --no-intro \
    > log_bulk.txt 2>&1
check_exit_code $? 1
check_contains "No vacuum plateau found" log_bulk.txt


# --- 7. Error and robustness cases ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: --axis outside 0-2 is rejected"
stb-workfunction -l symmetric --axis 5 --no-intro > log_badaxis.txt 2>&1
check_exit_code $? 2

echo "Testing: missing grid file"
stb-workfunction -l does_not_exist --fermi 0.0 --no-plot --no-intro > log_missing_grid.txt 2>&1
check_exit_code $? 1
check_contains "Grid file" log_missing_grid.txt

echo "Testing: no Fermi energy available (no .out, no --fermi)"
stb-workfunction -l no_fermi_label --no-plot --no-intro > log_missing_fermi.txt 2>&1
check_exit_code $? 1
check_contains "Could not determine Fermi Energy" log_missing_fermi.txt

echo "Testing: --label omitted (required)"
stb-workfunction --no-intro > log_missing_label.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-workfunction --version > log_version.txt 2>&1
check_contains "stb-workfunction" log_version.txt

echo "Testing: --help documents --label, --axis, --fermi, --no-plot, --asymmetric-tol"
stb-workfunction --help > log_help.txt 2>&1
check_contains "label" log_help.txt
check_contains "axis" log_help.txt
check_contains "fermi" log_help.txt
check_contains "no-plot" log_help.txt
check_contains "asymmetric-tol" log_help.txt


# --- 8. Interactive path (stb-suite, shortcut 2.7) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 2.7) ---"

echo "Testing: navigate 2.7 -> label symmetric -> default grid -> default fermi file -> blank axis (auto) -> quit"
printf '2.7\nsymmetric\n\n\n\n\n0\n' | timeout 60 stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_contains "Auto-detected vacuum axis: 2" log_menu.txt


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
