#!/bin/bash

# --- Setup ---
# Smoke test for stb-oerIntermediates (OER Stage 2: O*/OOH* Intermediates,
# item 4.14.2). Chains a real stb-oer (Stage 1, tested separately) run on
# the graphene fixture, then fabricates a relaxed calc.out (known FreeEng
# + a shifted "outcoor: Relaxed atomic coordinates" block) for the single
# candidate OH* site.
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREP_DIR="$(cd "$FIXTURE_DIR/../prep" && pwd)"
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

# Runs stb-oer (Stage 1, tested separately) on the graphene fixture with a
# single ontop site, then fabricates a relaxed calc.out for it: known
# FreeEng, and a relaxed-coordinates block where O/H have moved slightly
# from their Stage-1 input positions.
make_site_disp() {
    rm -rf oer_study
    stb-oer -s structure.fdf -c calc.fdf -p . --site-type ontop --height 1.6 --no-intro \
        > /dev/null 2>&1
    python3 - <<'PYEOF'
with open("oer_study/sites/site_1_ontop/calc.out", "w") as f:
    f.write("siesta: SCF Convergence by DM criterion\n")
    f.write("SCF cycle converged after 15 iterations\n")
    f.write("siesta: FreeEng =        -499.950000\n")
    f.write("outcoor: Relaxed atomic coordinates (fractional)\n")
    f.write("    0.00000000    0.00000000    0.50000000   1\n")
    f.write("    0.33333333    0.66666667    0.50000000   1\n")
    f.write("    0.00100000    0.00100000    0.58000000   2\n")
    f.write("    0.00100000    0.00100000    0.62000000   3\n")
    f.write("siesta: Atomic forces (eV/Ang):\n")
    f.write("     1    0.001    0.001    0.001\n")
    f.write("     2    0.001    0.001    0.001\n")
    f.write("     3    0.001    0.001    0.001\n")
    f.write("     4    0.001    0.001    0.001\n")
    f.write("     Max    0.001\n")
PYEOF
}


# --- 1. Preparation ---
echo "--- Starting tester for STB-OERIntermediates stage 2: O*/OOH* intermediates (item 4.14.2) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$PREP_DIR/structure.fdf" "$TEST_DIR/"
cp "$PREP_DIR/calc.fdf" "$TEST_DIR/"
echo "# placeholder pseudopotential" > "$TEST_DIR/C.psf"
echo "# placeholder pseudopotential" > "$TEST_DIR/O.psf"
echo "# placeholder pseudopotential" > "$TEST_DIR/H.psf"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Missing Stage-1 output guard ---
echo -e "\n--- Testing the missing-Stage-1-output guard ---"
rm -rf oer_study
stb-oerIntermediates --directory oer_study -p . --no-intro > log_no_stage1.txt 2>&1
check_exit_code $? 1
check_contains "run stb-oer" log_no_stage1.txt


# --- 3. Default run (O*/OOH* always derived from the winning OH* site) ---
echo -e "\n--- Testing default generation (O*/OOH* derived from the winning OH* site) ---"
make_site_disp
stb-oerIntermediates --directory oer_study -p . --no-intro > log_derived.txt 2>&1
check_exit_code $? 0
check_contains "Winning OH\* site : site_1_ontop" log_derived.txt
check_success oer_study/intermediates/o_star/structure.fdf
check_success oer_study/intermediates/o_star/calc.fdf
check_success oer_study/intermediates/ooh_star/structure.fdf
check_success oer_study/intermediates/ooh_star/calc.fdf
check_success oer_study/oer_stage2.txt

echo "Testing: O* has 3 atoms (2C + O), no H"
check_contains "NumberofAtoms      3" oer_study/intermediates/o_star/structure.fdf

echo "Testing: OOH* has 5 atoms (2C + O + O + H)"
check_contains "NumberofAtoms      5" oer_study/intermediates/ooh_star/structure.fdf

echo "Testing: every derived folder forces MD.TypeOfRun CG (real relaxation, not single-point)"
check_contains "MD.TypeOfRun          CG" oer_study/intermediates/o_star/calc.fdf
check_contains "MD.TypeOfRun          CG" oer_study/intermediates/ooh_star/calc.fdf
check_not_contains "MD.Steps              0" oer_study/intermediates/o_star/calc.fdf

echo "Testing: fixed cell + Slab.DipoleCorrection + Spin polarized + DFTD3 (all mandatory)"
echo "         forced via config_extra.fdf (4.8/4.11/4.12/4.13 model), %include'd on top of"
echo "         the untouched calc.fdf, not edited in place"
check_contains "%include config_extra.fdf" oer_study/intermediates/o_star/calc.fdf
check_contains "%include config_extra.fdf" oer_study/intermediates/ooh_star/calc.fdf
check_success oer_study/intermediates/o_star/config_extra.fdf
check_success oer_study/intermediates/ooh_star/config_extra.fdf
check_contains "MD.VariableCell false" oer_study/intermediates/o_star/config_extra.fdf
check_contains "Slab.DipoleCorrection      .true." oer_study/intermediates/o_star/config_extra.fdf
check_contains "Spin                polarized" oer_study/intermediates/o_star/config_extra.fdf
check_contains "DFTD3                   .true." oer_study/intermediates/o_star/config_extra.fdf
check_contains "Slab.DipoleCorrection      .true." oer_study/intermediates/ooh_star/config_extra.fdf

echo "Testing: derived O*/OOH* geometries also saved as single-frame extended-XYZ (OVITO/VMD),"
echo "         with bare (real) element symbols, not the fragment-labeled O_ads/H_ads"
check_success oer_study/intermediates/o_trajectory.xyz
check_success oer_study/intermediates/ooh_trajectory.xyz
check_contains "site_label=o_star" oer_study/intermediates/o_trajectory.xyz
check_contains "site_label=ooh_star" oer_study/intermediates/ooh_trajectory.xyz
check_not_contains "O_ads" oer_study/intermediates/o_trajectory.xyz
python3 -c "
import ase.io
o = ase.io.read('oer_study/intermediates/o_trajectory.xyz', index=':')
ooh = ase.io.read('oer_study/intermediates/ooh_trajectory.xyz', index=':')
assert len(o) == 1 and len(ooh) == 1, (len(o), len(ooh))
assert sorted(o[0].get_chemical_symbols()) == ['C', 'C', 'O'], o[0].get_chemical_symbols()
assert sorted(ooh[0].get_chemical_symbols()) == ['C', 'C', 'H', 'O', 'O'], ooh[0].get_chemical_symbols()
print('OK')
" > log_derived_trajectory_check.txt 2>&1
check_contains "OK" log_derived_trajectory_check.txt


# --- 3b. --ml-prerelax (MACE-MP-0 already cached locally) ---
echo -e "\n--- Testing --ml-prerelax (derived O*/OOH*, substrate fixed) ---"
make_site_disp
stb-oerIntermediates --directory oer_study -p . --ml-prerelax --no-intro > log_mlprerelax.txt 2>&1
check_exit_code $? 0
check_contains "ML pre-relax" log_mlprerelax.txt
check_success oer_study/intermediates/o_star/structure.fdf
check_success oer_study/intermediates/ooh_star/structure.fdf
check_contains "NumberofAtoms      3" oer_study/intermediates/o_star/structure.fdf
check_contains "NumberofAtoms      5" oer_study/intermediates/ooh_star/structure.fdf

echo "Testing: substrate (first 2 atoms) unchanged by --ml-prerelax (FixAtoms)"
python3 -c "
import re

def frac_rows(path):
    with open(path) as f:
        text = f.read()
    m = re.search(r'%block AtomicCoordinatesAndAtomicSpecies\n(.*?)%endblock', text, re.DOTALL)
    return [[float(x) for x in r.split()[:3]] for r in m.group(1).strip().split(chr(10))]

def close(a, b, tol=1e-6):
    return all(abs(x - y) < tol for x, y in zip(a, b))

rows = frac_rows('oer_study/intermediates/o_star/structure.fdf')
# substrate = first 2 rows (the graphene C atoms), must be untouched (FixAtoms)
assert close(rows[0], [0.0, 0.0, 0.5]), f'substrate atom 1 moved: {rows[0]}'
assert close(rows[1], [0.33333333, 0.66666667, 0.5]), f'substrate atom 2 moved: {rows[1]}'
print('OK')
" > log_mlprerelax_substrate_check.txt 2>&1
check_contains "OK" log_mlprerelax_substrate_check.txt


# --- 3c. OOH* orientation sampling, unscreened (no --ml-prerelax) ---
echo -e "\n--- Testing --ooh-n-orientations-polar/-azimuthal WITHOUT --ml-prerelax (unscreened) ---"
make_site_disp
stb-oerIntermediates --directory oer_study -p . --ooh-n-orientations-polar 2 \
    --ooh-n-orientations-azimuthal 2 --no-intro > log_orient_unscreened.txt 2>&1
check_exit_code $? 0
check_contains "unscreened" log_orient_unscreened.txt
check_success oer_study/intermediates/ooh_star_orient1/structure.fdf
check_success oer_study/intermediates/ooh_star_orient2/structure.fdf
check_success oer_study/intermediates/ooh_star_orient3/structure.fdf
check_success oer_study/intermediates/ooh_star_orient4/structure.fdf

echo "Testing: O* is unaffected (still a single ooh_star-style folder, no orientation DOF)"
check_success oer_study/intermediates/o_star/structure.fdf
python3 -c "
import os
assert not os.path.isdir('oer_study/intermediates/ooh_star'), \
    'single ooh_star/ should not exist once orientation sampling is requested'
print('OK')
" > log_orient_no_single_check.txt 2>&1
check_contains "OK" log_orient_no_single_check.txt

echo "Testing: WITHOUT --ml-prerelax nothing relaxes O1 before writing, so every orientation's"
echo "         O1 sits at the EXACT SAME Cartesian position (the winning OH* site's own relaxed"
echo "         oxygen) -- orientation sampling never moves the anchor/site (confirms 'mesmo sitio')"
python3 -c "
import glob
import numpy as np
from stb.core import structure_io

dirs = sorted(glob.glob('oer_study/intermediates/ooh_star_orient*'))
assert len(dirs) == 4
o1_positions = []
for d in dirs:
    s = structure_io.read_fdf(f'{d}/structure.fdf')
    sym, frac = s.atoms[2]  # 2 substrate C + O1, then O2/Hnew appended last
    assert sym.startswith('O'), f'expected O1 at index 2, got {sym} in {d}'
    o1_positions.append(np.array(frac) @ s.lattice)
ref = o1_positions[0]
for i, pos in enumerate(o1_positions[1:], start=2):
    assert np.allclose(pos, ref, atol=1e-9), \
        f'O1 moved between unscreened orientation candidates: {ref} vs {pos} (candidate {i})'
print('OK')
" > log_orient_same_site_unscreened_check.txt 2>&1
check_contains "OK" log_orient_same_site_unscreened_check.txt


# --- 3d. OOH* orientation sampling WITH --ml-prerelax + --orientation-top-k ---
echo -e "\n--- Testing --ooh-n-orientations-polar/-azimuthal WITH --ml-prerelax --orientation-top-k 2 ---"
make_site_disp
stb-oerIntermediates --directory oer_study -p . --ooh-n-orientations-polar 2 \
    --ooh-n-orientations-azimuthal 2 --ml-prerelax --orientation-top-k 2 --no-intro \
    > log_orient_ranked.txt 2>&1
check_exit_code $? 0
check_contains "E_MACE" log_orient_ranked.txt
check_contains "unique kept" log_orient_ranked.txt
check_success oer_study/intermediates/ooh_star_orient1/structure.fdf

echo "Testing: at most 2 orientation folders kept (--orientation-top-k 2)"
python3 -c "
import glob
dirs = sorted(glob.glob('oer_study/intermediates/ooh_star_orient*'))
assert 1 <= len(dirs) <= 2, f'expected 1-2 kept orientations, got {len(dirs)}: {dirs}'
print('OK')
" > log_orient_topk_check.txt 2>&1
check_contains "OK" log_orient_topk_check.txt

echo "Testing: WITH --ml-prerelax, O1 is free to relax too (adsorbate = O1+O2+H, same"
echo "         convention as the single-orientation --ml-prerelax path -- O1 is NOT frozen"
echo "         during MACE screening), so exact equality isn't expected here -- but every kept"
echo "         orientation's O1 must still stay CLOSE to the winning OH* site's own oxygen"
echo "         (a small MACE relaxation drift, not a jump to a whole different part of the"
echo "         surface the way the old, now-removed per-intermediate site search could)"
python3 -c "
import glob
import numpy as np
from stb.core import structure_io

dirs = sorted(glob.glob('oer_study/intermediates/ooh_star_orient*'))
assert len(dirs) >= 1
oh_structure = structure_io.read_fdf('oer_study/sites/site_1_ontop/structure.fdf')
sym, frac = oh_structure.atoms[2]  # OH*'s own O (2 substrate C + O + H)
assert sym.startswith('O'), sym
oh_o_cart = np.array(frac) @ oh_structure.lattice
for d in dirs:
    s = structure_io.read_fdf(f'{d}/structure.fdf')
    sym, frac = s.atoms[2]
    assert sym.startswith('O'), f'expected O1 at index 2, got {sym} in {d}'
    cart = np.array(frac) @ s.lattice
    dist = np.linalg.norm(cart - oh_o_cart)
    # Tolerance is loose (this tiny 2-atom-substrate graphene fixture has much less steric
    # constraint than a real slab, so MACE can legitimately relax O1 several Ang here) --
    # the invariant under test is same site (a few Ang), not zero movement; a genuinely
    # different SITE (the old, removed search mode's failure) would show a much larger
    # separation, comparable to the distance between distinct adsorption sites.
    assert dist < 3.0, f'{d}: O1 drifted {dist:.3f} Ang from the winning OH* site -- too far ' \
        'for a same-site MACE relaxation, looks like a different site'
print('OK')
" > log_orient_same_site_check.txt 2>&1
check_contains "OK" log_orient_same_site_check.txt

echo "Testing: --help documents the OOH* orientation-sampling flags"
stb-oerIntermediates --help > log_help_orient.txt 2>&1
check_contains "ooh-n-orientations-polar" log_help_orient.txt
check_contains "ooh-n-orientations-azimuthal" log_help_orient.txt
check_contains "orientation-top-k" log_help_orient.txt
check_contains "orientation-rmsd-tol" log_help_orient.txt


# --- 4. Error cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: --version"
stb-oerIntermediates --version > log_version.txt 2>&1
check_contains "stb-oerIntermediates" log_version.txt

echo "Testing: --help documents --oo-bond-length/--ml-prerelax, and no longer offers"
echo "         a search strategy (removed -- O*/OOH* are always derived from the same site)"
stb-oerIntermediates --help > log_help.txt 2>&1
check_contains "oo-bond-length" log_help.txt
check_contains "ml-prerelax" log_help.txt
check_contains "ml-model" log_help.txt
check_not_contains "site search" log_help.txt

echo "Testing: no --structure/--calc/--o-strategy/--ooh-strategy flags exist anymore (removed"
echo "         entirely, not just deprecated -- search mode/standalone mode no longer exist)"
stb-oerIntermediates --o-strategy derived --no-intro > log_removed_flag.txt 2>&1
check_exit_code $? 2
stb-oerIntermediates --structure structure.fdf --calc calc.fdf --no-intro > log_removed_flag2.txt 2>&1
check_exit_code $? 2


# --- 5. Interactive path (stb-suite, shortcut 4.14.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.14.2) ---"

echo "Testing: navigate 4.14.2 -> defaults -> quit"
make_site_disp
{
  echo "4.14.2"
  echo ""               # run_dir (default oer_study)
  echo ""               # output_filename (default calc.out)
  echo "3"              # pseudo_dir -> option 3 = Custom path
  echo "."              # custom pseudo path
  echo ""               # ml_prerelax_choice (default N)
  echo ""               # ooh_orient_grid (blank = skip orientation sampling)
  echo ""               # show_advanced (default N)
  echo ""               # press enter to continue
  echo "0"              # quit stage submenu
} | stb-suite > log_menu.txt 2>&1
check_success oer_study/intermediates/o_star/structure.fdf
check_success oer_study/intermediates/ooh_star/structure.fdf


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
