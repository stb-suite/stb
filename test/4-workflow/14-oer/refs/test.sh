#!/bin/bash

# --- Setup ---
# Smoke test for stb-oerRefs (OER Stage 3: References, BSSE & ZPE Prep,
# item 4.14.3). Chains real stb-oer + stb-oerIntermediates (both --strategy
# derived, tested separately) runs, then fabricates "relaxed" calc.out
# files (known FreeEng + relaxed-coordinates blocks) for the winning OH*
# site and for intermediates/o_star, intermediates/ooh_star.
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

check_exit_code() {
    if [ "$1" -eq "$2" ]; then
        echo -e "   -> ${GREEN}Verified:${NC} exit code $1 (expected $2)"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} exit code $1 (expected $2)"
        FAIL=$((FAIL+1))
    fi
}

# Full chain: Stage1 (1 ontop site) -> fabricate OH* relaxed calc.out ->
# Stage2 (derived/derived) -> fabricate relaxed calc.out for o_star and
# ooh_star.
make_full_chain() {
    rm -rf oer_study
    stb-oer -s structure.fdf -c calc.fdf -p . --site-type ontop --height 1.6 --no-intro \
        > /dev/null 2>&1
    python3 - <<'PYEOF'
with open("oer_study/sites/site_1_ontop/calc.out", "w") as f:
    f.write("siesta: FreeEng =        -499.950000\n")
    f.write("outcoor: Relaxed atomic coordinates (fractional)\n")
    f.write("    0.00000000    0.00000000    0.50000000   1\n")
    f.write("    0.33333333    0.66666667    0.50000000   1\n")
    f.write("    0.00100000    0.00100000    0.58000000   2\n")
    f.write("    0.00100000    0.00100000    0.62000000   3\n")
    f.write("siesta: Atomic forces (eV/Ang):\n")
    f.write("     Max    0.001\n")
PYEOF
    stb-oerIntermediates --directory oer_study -p . --no-intro > /dev/null 2>&1
    python3 - <<'PYEOF'
def write_relaxed(path, freeeng, coords):
    with open(path, "w") as f:
        f.write(f"siesta: FreeEng =        {freeeng:.6f}\n")
        f.write("outcoor: Relaxed atomic coordinates (fractional)\n")
        for row in coords:
            f.write(f"    {row[0]:.8f}    {row[1]:.8f}    {row[2]:.8f}   {row[3]}\n")
        f.write("siesta: Atomic forces (eV/Ang):\n")
        f.write("     Max    0.001\n")

# O* = 3 atoms (2C + O)
write_relaxed("oer_study/intermediates/o_star/calc.out", -480.000000, [
    (0.0, 0.0, 0.5, 1), (0.333333333, 0.666666667, 0.5, 1), (0.001, 0.001, 0.58, 2),
])
# OOH* = 5 atoms (2C + O + O + H)
write_relaxed("oer_study/intermediates/ooh_star/calc.out", -520.000000, [
    (0.0, 0.0, 0.5, 1), (0.333333333, 0.666666667, 0.5, 1), (0.001, 0.001, 0.58, 2),
    (0.001, 0.001, 0.66, 2), (0.001, 0.001, 0.70, 3),
])
PYEOF
}


# --- 1. Preparation ---
echo "--- Starting tester for STB-OERRefs stage 3: references, BSSE & ZPE prep (item 4.14.3) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$PREP_DIR/structure.fdf" "$TEST_DIR/"
cp "$PREP_DIR/calc.fdf" "$TEST_DIR/"
echo "# placeholder pseudopotential" > "$TEST_DIR/C.psf"
echo "# placeholder pseudopotential" > "$TEST_DIR/O.psf"
echo "# placeholder pseudopotential" > "$TEST_DIR/H.psf"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Missing Stage-2 output guard ---
echo -e "\n--- Testing the missing-Stage-2-output guard ---"
rm -rf oer_study
stb-oerRefs --directory oer_study -p . --no-intro > log_no_stage2.txt 2>&1
check_exit_code $? 1
check_contains "run stb-oerIntermediates" log_no_stage2.txt


# --- 3. Default run (--bsse-mode shared --zpe-mode local) ---
echo -e "\n--- Testing default generation (shared BSSE, local ZPE) ---"
make_full_chain
stb-oerRefs --directory oer_study -p . --no-intro > log_local.txt 2>&1
check_exit_code $? 0
check_success oer_study/00_clean_slab/structure.fdf
check_success oer_study/02_h2_molecule/structure.fdf
check_success oer_study/03_h2o_molecule/structure.fdf
check_success oer_study/04_slab_deformed/structure.fdf
check_success oer_study/05_bsse_slab_only/structure.fdf
check_success oer_study/05_bsse_slab_ghost/structure.fdf
check_success oer_study/05_bsse_adsorbate_ghost_slab/structure.fdf
check_success oer_study/05_bsse_isolated/structure.fdf
check_success oer_study/08_zpe_calc_OH/disp_001/structure.fdf
check_success oer_study/08_zpe_calc_OH/disp_012/structure.fdf
check_success oer_study/08_zpe_calc_O/disp_001/structure.fdf
check_success oer_study/08_zpe_calc_O/disp_006/structure.fdf
check_success oer_study/08_zpe_calc_OOH/disp_001/structure.fdf
check_success oer_study/08_zpe_calc_OOH/disp_018/structure.fdf
check_success oer_study/08_zpe_calc_H2O/disp_001/structure.fdf
check_success oer_study/08_zpe_calc_H2O/disp_018/structure.fdf
check_success oer_study/oer_stage3.txt

echo "Testing: O*'s local ZPE folder count (6, disp_007 must NOT exist) matches HER's own single-H-atom case exactly"
python3 -c "
import os
assert os.path.isdir('oer_study/08_zpe_calc_O/disp_006')
assert not os.path.isdir('oer_study/08_zpe_calc_O/disp_007'), 'O* should have exactly 6 local displacement folders'
print('OK')
" > log_o_zpe_count.txt 2>&1
check_contains "OK" log_o_zpe_count.txt

echo "Testing: 04_slab_deformed has no O/H (2 atoms, C only -- diagnostic geometry, OH* based)"
check_contains "NumberofAtoms      2" oer_study/04_slab_deformed/structure.fdf

echo "Testing: BSSE ghost species present (O_ghost/H_ghost on slab_ghost, C_ghost on adsorbate_ghost_slab)"
check_contains "O_ghost" oer_study/05_bsse_slab_ghost/structure.fdf
check_contains "C_ghost" oer_study/05_bsse_adsorbate_ghost_slab/structure.fdf

echo "Testing: 05_bsse_isolated has exactly 3 atoms (OOH* is the shared triad's geometry)"
check_contains "NumberofAtoms      3" oer_study/05_bsse_isolated/structure.fdf

echo "Testing: H2/H2O molecules relax (not single-point), Gamma-only + spin-polarized"
check_contains "MD.TypeOfRun          CG" oer_study/02_h2_molecule/calc.fdf
check_contains "MD.TypeOfRun          CG" oer_study/03_h2o_molecule/calc.fdf
check_contains "kgrid.MonkhorstPack   \[1  1  1\]" oer_study/02_h2_molecule/calc.fdf
check_contains "Spin.*polarized" oer_study/03_h2o_molecule/calc.fdf

echo "Testing: 00_clean_slab/04_slab_deformed/BSSE folders are forced single-point"
check_contains "MD.NumCGsteps         0" oer_study/00_clean_slab/calc.fdf
check_contains "MD.NumCGsteps         0" oer_study/04_slab_deformed/calc.fdf
check_contains "MD.NumCGsteps         0" oer_study/05_bsse_slab_only/calc.fdf


# --- 4. --bsse-mode full ---
echo -e "\n--- Testing --bsse-mode full (3 separate triads) ---"
make_full_chain
stb-oerRefs --directory oer_study -p . --bsse-mode full --no-intro > log_bsse_full.txt 2>&1
check_exit_code $? 0
check_success oer_study/05_bsse_OH_slab_only/structure.fdf
check_success oer_study/05_bsse_O_slab_only/structure.fdf
check_success oer_study/05_bsse_OOH_slab_only/structure.fdf
python3 -c "
import os
assert not os.path.isdir('oer_study/05_bsse_slab_only'), 'shared triad should not exist in full mode'
print('OK')
" > log_bsse_full_check.txt 2>&1
check_contains "OK" log_bsse_full_check.txt


# --- 5. --zpe-mode full ---
echo -e "\n--- Testing --zpe-mode full (4 full phonon displacement sets: clean+OH+O+OOH+H2O) ---"
make_full_chain
stb-oerRefs --directory oer_study -p . --zpe-mode full --supercell 1 1 1 --no-intro \
    > log_zpe_full.txt 2>&1
check_exit_code $? 0
check_success oer_study/09_zpe_calc_clean/disp-001/structure.fdf
check_success oer_study/09_zpe_calc_OH/disp-001/structure.fdf
check_success oer_study/09_zpe_calc_O/disp-001/structure.fdf
check_success oer_study/09_zpe_calc_OOH/disp-001/structure.fdf
check_success oer_study/09_zpe_calc_H2O/disp-001/structure.fdf
check_contains "SystemLabel oer_zpe_clean" oer_study/09_zpe_calc_clean/disp-001/calc.fdf
check_contains "SystemLabel oer_zpe_ooh" oer_study/09_zpe_calc_OOH/disp-001/calc.fdf


# --- 6. Error cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: --version"
stb-oerRefs --version > log_version.txt 2>&1
check_contains "stb-oerRefs" log_version.txt

echo "Testing: --help documents --bsse-mode/--zpe-mode/--displacement/--supercell"
stb-oerRefs --help > log_help.txt 2>&1
check_contains "bsse-mode" log_help.txt
check_contains "zpe-mode" log_help.txt
check_contains "displacement" log_help.txt
check_contains "supercell" log_help.txt

echo "Testing: --zpe-mode standard is rejected (not offered in v1)"
stb-oerRefs --directory oer_study --zpe-mode standard --no-intro > log_standard_reject.txt 2>&1
check_exit_code $? 2
check_contains "invalid choice" log_standard_reject.txt


# --- 7. Interactive path (stb-suite, shortcut 4.14.3) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.14.3) ---"

echo "Testing: navigate 4.14.3 -> shared/local -> quit"
make_full_chain
{
  echo "4.14.3"
  echo ""               # run_dir (default oer_study)
  echo ""               # output_filename (default calc.out)
  echo "3"              # pseudo_dir -> option 3 = Custom path
  echo "."              # custom pseudo path
  echo ""               # bsse_mode (default shared)
  echo ""               # zpe_mode (default local)
  echo ""               # show_advanced (default N)
  echo ""               # press enter to continue
  echo "0"              # quit stage submenu
} | stb-suite > log_menu.txt 2>&1
check_success oer_study/00_clean_slab/structure.fdf
check_success oer_study/08_zpe_calc_OH/disp_001/structure.fdf


popd > /dev/null

# --- 8. Summary ---
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
