#!/bin/bash

# --- Setup ---
# Smoke test for stb-adsorbGibbs (Adsorption Gibbs Free Energy, item 4.8.4).
# Builds a real site + isolated-adsorbate reference via stb-adsorb, fabricates
# relaxed siesta.XV files (same sisl-write recipe already used by
# ../bsse/test.sh), runs stb-adsorbAnalysis --compute-gibbs to generate the
# 'gibbs/' Hessian displacement folders (this tool's own prep step), then
# fabricates .FA force files for those folders to exercise stb-adsorbGibbs
# itself without needing a real SIESTA run.
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


# --- 1. Preparation: a real relaxed site + isolated-adsorbate reference,
#     then stb-adsorbAnalysis --compute-gibbs writes the Hessian folders
#     this tool itself needs. ---
echo "--- Starting tester for STB-AdsorbGibbs (item 4.8.4) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$PREP_DIR/structure.fdf" "$PREP_DIR/calc.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null

stb-adsorb -s structure.fdf -c calc.fdf --adsorbate O --site-type ontop -O . --no-intro \
    > log_prep.txt 2>&1
check_success sites/site_1_ontop/structure.fdf

printf 'siesta: FreeEng =    -428.500000\nSCF cycle converged after 14 iterations\nsiesta: Atomic forces (eV/Ang):\n   Max    0.020000\n' \
    > sites/site_1_ontop/calc.out
printf 'siesta: FreeEng =    -213.900000\nSCF cycle converged after 10 iterations\nsiesta: Atomic forces (eV/Ang):\n   Max    0.010000\n' \
    > clean_slab/calc.out
printf 'siesta: FreeEng =    -214.500000\nSCF cycle converged after 8 iterations\nsiesta: Atomic forces (eV/Ang):\n   Max    0.005000\n' \
    > adsorbate/calc.out
# E_ads (raw) = -428.5 - (-213.9) - (-214.5) = -0.100000 eV -- checked below
# against stb-adsorbGibbs's own [1] ELECTRONIC ENERGIES section.

python3 -c "
import sisl
from stb.core import structure_io

for site_dir in ('sites/site_1_ontop', 'adsorbate'):
    fdf = structure_io.read_fdf(f'{site_dir}/structure.fdf')
    pmg = structure_io.to_pymatgen(fdf)
    cart = pmg.cart_coords.copy()
    if 'sites' in site_dir:
        cart[-1, 2] -= 0.4  # relax the O atom closer to the substrate
    atoms = [sisl.Atom(str(s.specie)) for s in pmg]
    geom = sisl.Geometry(cart, atoms=atoms, lattice=sisl.Lattice(pmg.lattice.matrix))
    sisl.get_sile(f'{site_dir}/siesta.XV', mode='w').write_geometry(geom)
"
check_success sites/site_1_ontop/siesta.XV
check_success adsorbate/siesta.XV

stb-adsorbAnalysis --dir . --compute-gibbs --zpe-mode local --no-intro > log_analysis.txt 2>&1
check_success gibbs/site_1_ontop/gibbs_local_meta.json
check_success gibbs/O_isolated/gibbs_local_meta.json


# --- 2. Fabricate .FA force files: an IDENTICAL isotropic harmonic
#     "spring" (k=20 eV/Ang^2) on BOTH the site's O atom and the isolated
#     reference's O atom -- same mass, same force constant, so their
#     vibrational frequencies/ZPE/entropy are mathematically IDENTICAL and
#     must cancel exactly: DZPE = DTS = 0, so DG must come out EXACTLY
#     equal to the electronic-only E_ads (-0.100000 eV) -- a strong,
#     exact-value cross-check of the whole Hessian/ZPE/entropy/DG pipeline,
#     not just "did it run". ---
echo -e "\n--- Testing DG with a symmetric synthetic Hessian (DZPE = DTS = 0 by construction) ---"
python3 -c "
k = 20.0  # eV/Ang^2
d = 0.015
order = [(0,1.0),(0,-1.0),(1,1.0),(1,-1.0),(2,1.0),(2,-1.0)]

def write_fa(path, natoms, moved_1based, force):
    with open(path, 'w') as f:
        f.write(f'{natoms}\n')
        for i in range(1, natoms + 1):
            if i == moved_1based:
                f.write(f'{i} {force[0]:.8f} {force[1]:.8f} {force[2]:.8f}\n')
            else:
                f.write(f'{i} 0.0 0.0 0.0\n')

for i, (axis, sign) in enumerate(order, start=1):
    force = [0.0, 0.0, 0.0]
    force[axis] = -sign * k * d
    write_fa(f'gibbs/site_1_ontop/disp_{i:03d}/gibbs_site.FA', 3, 3, force)
    write_fa(f'gibbs/O_isolated/disp_{i:03d}/gibbs_isolated.FA', 1, 1, force)
"
stb-adsorbGibbs --dir . --tmin 300 --tmax 300 --tstep 25 --save-report --no-intro \
    > log_gibbs.txt 2>&1
check_exit_code $? 0
check_contains "\[0\] RUN METADATA" log_gibbs.txt
check_contains "\[1\] ELECTRONIC ENERGIES" log_gibbs.txt
check_contains "E_ads (raw)  =    -0.100000 eV" log_gibbs.txt
check_contains "E_ads used for DG below: raw" log_gibbs.txt
check_contains "\[2\] VIBRATIONAL/THERMAL TERMS" log_gibbs.txt
check_contains "LIMITATION" log_gibbs.txt
check_contains "\[3\] GIBBS FREE ENERGY vs. TEMPERATURE" log_gibbs.txt
check_contains "DZPE = +0.0000 eV" log_gibbs.txt
check_contains "DTS = +0.0000 eV" log_gibbs.txt
check_contains "DG = -0.1000 eV" log_gibbs.txt
check_contains "FINAL RESULT" log_gibbs.txt
check_contains "DG(adsorption) = -0.1000 eV" log_gibbs.txt
check_success adsorption_gibbs.png
check_success adsorption_gibbs_report.txt


# --- 3. Physical correctness: a truly isolated single atom (zero force at
#     every displacement -- nothing to restore it) must give EXACTLY
#     zero ZPE/entropy once its 3 pure-translation modes are excluded --
#     not close to zero, exactly zero. Verified by calling the module's own
#     function directly (not through the CLI, to get the raw float back). ---
echo -e "\n--- Testing that a genuinely free atom gives exactly ZPE=0 ---"
python3 -c "
def write_fa(path, natoms, moved_1based, force):
    with open(path, 'w') as f:
        f.write(f'{natoms}\n')
        for i in range(1, natoms + 1):
            f.write(f'{i} 0.0 0.0 0.0\n')

for i in range(1, 7):
    write_fa(f'gibbs/O_isolated/disp_{i:03d}/gibbs_isolated.FA', 1, 1, [0.0, 0.0, 0.0])

from stb.adsorb_gibbs import compute_local_hessian_thermo
zpe, ts = compute_local_hessian_thermo('gibbs/O_isolated', 300.0, None, 'test')
assert zpe == 0.0, f'expected exactly 0.0 ZPE for a free atom, got {zpe}'
assert ts == 0.0, f'expected exactly 0.0 TS for a free atom, got {ts}'
print('OK')
" > log_free_atom_check.txt 2>&1
check_contains "OK" log_free_atom_check.txt
# restore the symmetric-Hessian fixture for the sections below
python3 -c "
k = 20.0
d = 0.015
order = [(0,1.0),(0,-1.0),(1,1.0),(1,-1.0),(2,1.0),(2,-1.0)]
for i, (axis, sign) in enumerate(order, start=1):
    force = [0.0, 0.0, 0.0]
    force[axis] = -sign * k * d
    with open(f'gibbs/O_isolated/disp_{i:03d}/gibbs_isolated.FA', 'w') as f:
        f.write('1\n')
        f.write(f'1 {force[0]:.8f} {force[1]:.8f} {force[2]:.8f}\n')
"


# --- 4. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: missing 'gibbs/' directory entirely"
mkdir -p no_gibbs_dir
stb-adsorbGibbs --dir no_gibbs_dir --no-intro > log_no_gibbs.txt 2>&1
check_exit_code $? 1
check_contains "run stb-adsorbAnalysis --compute-gibbs" log_no_gibbs.txt

echo "Testing: ambiguous 'gibbs/' layout (two site-like folders, ni isolated ref missing)"
mkdir -p ambiguous_dir/gibbs/site_a ambiguous_dir/gibbs/site_b
stb-adsorbGibbs --dir ambiguous_dir --no-intro > log_ambiguous.txt 2>&1
check_exit_code $? 1
check_contains "Expected exactly one" log_ambiguous.txt

echo "Testing: a missing .FA file is a clean error, not a crash"
mv gibbs/O_isolated/disp_001/gibbs_isolated.FA gibbs/O_isolated/disp_001/gibbs_isolated.FA.bak
stb-adsorbGibbs --dir . --tmin 300 --tmax 300 --no-intro > log_missing_fa.txt 2>&1
check_exit_code $? 1
check_contains "Could not read forces" log_missing_fa.txt
mv gibbs/O_isolated/disp_001/gibbs_isolated.FA.bak gibbs/O_isolated/disp_001/gibbs_isolated.FA

echo "Testing: --version"
stb-adsorbGibbs --version > log_version.txt 2>&1
check_contains "stb-adsorbGibbs" log_version.txt

echo "Testing: --help documents --dir/--tmin/--tmax"
stb-adsorbGibbs --help > log_help.txt 2>&1
check_contains "dir" log_help.txt
check_contains "tmin" log_help.txt
check_contains "tmax" log_help.txt


# --- 5. Interactive path (stb-suite, shortcut 4.8.4) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.8.4) ---"

echo "Testing: navigate 4.8.4 -> defaults -> quit"
rm -f adsorption_gibbs_report.txt adsorption_gibbs.png
# 4.8.4 (menu code) / . (dir) / "" (out_file default) / "" (tmin default) /
# "" (tmax default) / "" (tstep default) / "" (save_report: N) /
# "" (view_plots: N) / "" (Press Enter to continue) / 0 (quit)
printf '4.8.4\n.\n\n\n\n\n\n\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_contains "FINAL RESULT" log_menu.txt
check_success adsorption_gibbs.png


popd > /dev/null

# --- 6. Summary ---
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
