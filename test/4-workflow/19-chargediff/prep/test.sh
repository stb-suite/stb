#!/bin/bash

# --- Setup ---
# Smoke test for stb-chargediffPrep (Charge Density Difference, Stage 1 - Prep, item 4.19.1)
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
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


# --- 1. Preparation ---
echo "--- Starting tester for stb-chargediffPrep (item 4.19.1) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp -r "$FIXTURE_DIR/combined" "$TEST_DIR/"
cp -r "$FIXTURE_DIR/manifest_fixture" "$TEST_DIR/"
cp "$FIXTURE_DIR/calc_template.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/calc_no_mesh.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. --fragment-sizes: a 3-fragment split ---
echo -e "\n--- Testing --fragment-sizes \"2,2,2\" (3 fragments) ---"
stb-chargediffPrep --structure combined/structure.fdf --fragment-sizes "2,2,2" \
    --calc calc_template.fdf --no-intro --output-dir cdd_sizes > log_sizes.txt 2>&1
check_exit_code $? 0
check_success cdd_sizes/fragment1/structure.fdf
check_success cdd_sizes/fragment2/structure.fdf
check_success cdd_sizes/fragment3/structure.fdf
check_success cdd_sizes/chargediff_manifest.json

echo "Testing: fragment1 keeps its own real species and ghosts the other two"
check_contains "C_frag1" cdd_sizes/fragment1/structure.fdf
check_contains "B_frag2_ghost" cdd_sizes/fragment1/structure.fdf
check_contains "N_frag3_ghost" cdd_sizes/fragment1/structure.fdf

echo "Testing: config_extra.fdf mandatorily forces MeshCutoff (default 400 Ry, NOT the 300 Ry "\
"calc_template.fdf happens to declare)/SaveRho/fixed cell/single-point (no relaxation)"
check_contains "MeshCutoff          400 Ry" cdd_sizes/fragment1/config_extra.fdf
check_contains "SaveRho             true" cdd_sizes/fragment1/config_extra.fdf
check_contains "MD.VariableCell false" cdd_sizes/fragment1/config_extra.fdf
check_contains "MD.TypeOfRun          CG" cdd_sizes/fragment1/config_extra.fdf
check_contains "MD.Steps              0" cdd_sizes/fragment1/config_extra.fdf

echo "Testing: calc.fdf includes config_extra.fdf"
check_contains "%include config_extra.fdf" cdd_sizes/fragment1/calc.fdf


# --- 3. --fragment-species: a 2-fragment split by element ---
echo -e "\n--- Testing --fragment-species \"C;B,N\" (2 fragments) ---"
stb-chargediffPrep --structure combined/structure.fdf --fragment-species "C;B,N" \
    --calc calc_template.fdf --no-intro --output-dir cdd_species > log_species.txt 2>&1
check_exit_code $? 0
check_success cdd_species/fragment1/structure.fdf
check_success cdd_species/fragment2/structure.fdf
check_contains "C_frag1" cdd_species/fragment1/structure.fdf
check_contains "B_frag2_ghost" cdd_species/fragment1/structure.fdf
check_contains "N_frag2_ghost" cdd_species/fragment1/structure.fdf


# --- 4. --fragment-manifest: reuse an stb-adsorb-style split ---
echo -e "\n--- Testing --fragment-manifest (reusing a slab/adsorbate split) ---"
stb-chargediffPrep --structure manifest_fixture/structure.fdf \
    --fragment-manifest manifest_fixture/fragment_manifest.json --calc calc_template.fdf \
    --no-intro --output-dir cdd_manifest > log_manifest.txt 2>&1
check_exit_code $? 0
check_success cdd_manifest/slab/structure.fdf
check_success cdd_manifest/adsorbate/structure.fdf
check_contains "C_slab" cdd_manifest/slab/structure.fdf
check_contains "O_ads_ghost" cdd_manifest/slab/structure.fdf
check_contains "O_ads" cdd_manifest/adsorbate/structure.fdf
check_contains "C_slab_ghost" cdd_manifest/adsorbate/structure.fdf


# --- 5. Physics: refuse a --fragment-sizes that doesn't sum to the atom count ---
echo -e "\n--- Testing --fragment-sizes sum mismatch is rejected ---"
stb-chargediffPrep --structure combined/structure.fdf --fragment-sizes "2,2" \
    --calc calc_template.fdf --no-intro --output-dir cdd_bad_sizes > log_bad_sizes.txt 2>&1
check_exit_code $? 1
check_contains "must add up exactly" log_bad_sizes.txt


# --- 6. MeshCutoff: default (400 Ry, never auto-detected from --calc) and explicit override ---
echo -e "\n--- Testing the default MeshCutoff (400 Ry) is applied when --mesh-cutoff is omitted ---"
mkdir -p combined_nomesh
cp combined/structure.fdf combined_nomesh/
stb-chargediffPrep --structure combined_nomesh/structure.fdf --fragment-sizes "2,2,2" \
    --calc calc_no_mesh.fdf --no-intro --output-dir cdd_defaultmesh > log_defaultmesh.txt 2>&1
check_exit_code $? 0
check_contains "MeshCutoff          400 Ry" cdd_defaultmesh/fragment1/config_extra.fdf
check_contains "MeshCutoff          400 Ry" cdd_defaultmesh/fragment2/config_extra.fdf
check_contains "MeshCutoff          400 Ry" cdd_defaultmesh/fragment3/config_extra.fdf

echo "Testing: --mesh-cutoff explicit override is accepted, even though calc_template.fdf has "\
"its own (different) MeshCutoff line -- never auto-detected"
stb-chargediffPrep --structure combined/structure.fdf --fragment-sizes "2,2,2" \
    --calc calc_template.fdf --mesh-cutoff 250 --no-intro --output-dir cdd_explicit_mesh \
    > log_explicit_mesh.txt 2>&1
check_exit_code $? 0
check_contains "MeshCutoff          250 Ry" cdd_explicit_mesh/fragment1/config_extra.fdf
check_contains "MeshCutoff          250 Ry" cdd_explicit_mesh/fragment2/config_extra.fdf


# --- 6b. Defaults: --structure 'structure.fdf' and --calc 'calc.fdf' when omitted ---
echo -e "\n--- Testing --structure/--calc defaults ('structure.fdf' / 'calc.fdf') ---"
mkdir -p cdd_defaults_src
cp combined/structure.fdf cdd_defaults_src/
cp calc_template.fdf cdd_defaults_src/calc.fdf
( cd cdd_defaults_src && stb-chargediffPrep --fragment-sizes "2,2,2" \
    --no-intro --output-dir ../cdd_defaults_out > ../log_defaults.txt 2>&1 )
check_exit_code $? 0
check_success cdd_defaults_out/fragment1/structure.fdf
check_contains "Structure file  : structure.fdf" log_defaults.txt
check_contains "Calc template   : calc.fdf" log_defaults.txt


# --- 6c. -p/--pseudo-dir: bundled bank selection ---
echo -e "\n--- Testing -p dojo copies real pseudopotential files (not just a warning) ---"
stb-chargediffPrep --structure combined/structure.fdf --fragment-species "C;B,N" \
    --calc calc_template.fdf -p dojo --no-intro --output-dir cdd_pseudo_bank \
    > log_pseudo_bank.txt 2>&1
check_exit_code $? 0
check_success cdd_pseudo_bank/fragment1/C_frag1.psml
check_contains "Pseudo source" log_pseudo_bank.txt


# --- 6d. --view: writes fragments then opens ASE (bounded by timeout, safe headless or not) ---
echo -e "\n--- Testing --view writes fragments then attempts to open ASE without crashing ---"
rm -rf cdd_view
timeout 8 stb-chargediffPrep --structure combined/structure.fdf --fragment-species "C;B,N" \
    --calc calc_template.fdf --no-intro --output-dir cdd_view --view > log_view.txt 2>&1
view_exit=$?
if [ "$view_exit" -eq 0 ] || [ "$view_exit" -eq 124 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} exit code $view_exit (0 = no display/ase gui failed fast and was "\
"caught, 124 = timeout killed an actually-opened viewer window -- both are OK, neither is a crash)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} unexpected exit code $view_exit"
    FAIL=$((FAIL+1))
fi
check_success cdd_view/fragment1/structure.fdf
check_success cdd_view/fragment2/structure.fdf


# --- 7. Error and robustness cases ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: missing fragment-definition group (still required; --structure/--calc now default)"
stb-chargediffPrep --no-intro > log_missing_args.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-chargediffPrep --version > log_version.txt 2>&1
check_contains "stb-chargediffPrep" log_version.txt

echo "Testing: --help documents the fragment-definition flags and the new options"
stb-chargediffPrep --help > log_help.txt 2>&1
check_contains "fragment-sizes" log_help.txt
check_contains "fragment-species" log_help.txt
check_contains "fragment-manifest" log_help.txt
check_contains "mesh-cutoff" log_help.txt
check_contains "pseudo-dir" log_help.txt
check_contains "\-\-view" log_help.txt
check_contains "\-\-structure" log_help.txt


# --- 8. Interactive path (stb-suite, shortcut 4.19.1) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.19.1) ---"
rm -rf cdd_interactive
printf '4.19.1\ncombined/structure.fdf\n1\n2,2,2\ncalc_template.fdf\n\n\ncdd_interactive\nn\nn\n\n0\n' \
    | stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success cdd_interactive/fragment1/structure.fdf


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
