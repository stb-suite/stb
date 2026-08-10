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


# --- 3. Default run (--o-strategy derived --ooh-strategy derived) ---
echo -e "\n--- Testing default generation (both derived) ---"
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

echo "Testing: Slab.DipoleCorrection forced on both derived folders"
check_contains "Slab.DipoleCorrection   T" oer_study/intermediates/o_star/calc.fdf
check_contains "Slab.DipoleCorrection   T" oer_study/intermediates/ooh_star/calc.fdf


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


# --- 4. --o-strategy search ---
echo -e "\n--- Testing --o-strategy search (independent site search for O*) ---"
make_site_disp
stb-oerIntermediates --directory oer_study -p . --o-strategy search --o-site-type ontop \
    --ooh-strategy derived --no-intro > log_o_search.txt 2>&1
check_exit_code $? 0
check_success oer_study/intermediates/o_candidates/site_1_ontop/structure.fdf
check_contains "NumberofAtoms      3" oer_study/intermediates/o_candidates/site_1_ontop/structure.fdf
check_contains "MD.TypeOfRun          CG" oer_study/intermediates/o_candidates/site_1_ontop/calc.fdf
python3 -c "
import os
assert not os.path.isdir('oer_study/intermediates/o_star'), 'o_star should not exist in search mode'
print('OK')
" > log_o_search_check.txt 2>&1
check_contains "OK" log_o_search_check.txt


# --- 5. --ooh-strategy search ---
echo -e "\n--- Testing --ooh-strategy search (independent site search for OOH*) ---"
make_site_disp
stb-oerIntermediates --directory oer_study -p . --o-strategy derived \
    --ooh-strategy search --ooh-site-type ontop --no-intro > log_ooh_search.txt 2>&1
check_exit_code $? 0
check_success oer_study/intermediates/ooh_candidates/site_1_ontop/structure.fdf
check_contains "NumberofAtoms      5" oer_study/intermediates/ooh_candidates/site_1_ontop/structure.fdf


# --- 6. Error cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: --version"
stb-oerIntermediates --version > log_version.txt 2>&1
check_contains "stb-oerIntermediates" log_version.txt

echo "Testing: --help documents --o-strategy/--ooh-strategy/--oo-bond-length/--ml-prerelax"
stb-oerIntermediates --help > log_help.txt 2>&1
check_contains "o-strategy" log_help.txt
check_contains "ooh-strategy" log_help.txt
check_contains "oo-bond-length" log_help.txt
check_contains "ml-prerelax" log_help.txt
check_contains "ml-model" log_help.txt


# --- 7. Interactive path (stb-suite, shortcut 4.14.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.14.2) ---"

echo "Testing: navigate 4.14.2 -> defaults (derived/derived) -> quit"
make_site_disp
{
  echo "4.14.2"
  echo ""               # run_dir (default oer_study)
  echo ""               # output_filename (default calc.out)
  echo "3"              # pseudo_dir -> option 3 = Custom path
  echo "."              # custom pseudo path
  echo ""               # o_strategy (default derived)
  echo ""               # ooh_strategy (default derived)
  echo ""               # ml_prerelax_choice (default N)
  echo ""               # show_advanced (default N)
  echo ""               # press enter to continue
  echo "0"              # quit stage submenu
} | stb-suite > log_menu.txt 2>&1
check_success oer_study/intermediates/o_star/structure.fdf
check_success oer_study/intermediates/ooh_star/structure.fdf


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
