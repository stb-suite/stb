#!/bin/bash

# --- Setup ---
# Smoke test for stb-crystalcast (Random Crystal Generator, item 1.17)
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

# Checks that file $1 exists and is not empty
check_success() {
    if [ -s "$1" ]; then
        echo -e " ... ${GREEN}OK${NC} (file '$1' created)"
        PASS=$((PASS+1))
    else
        echo -e " ... ${RED}FAIL${NC} (file '$1' was not created)"
        FAIL=$((FAIL+1))
    fi
}

# Checks that file $2 contains (grep -q) pattern $1
check_contains() {
    if grep -q -- "$1" "$2" 2>/dev/null; then
        echo -e "   -> ${GREEN}Verified:${NC} '$1' found in '$2'"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} '$1' NOT found in '$2'"
        FAIL=$((FAIL+1))
    fi
}

# Checks that $1 (actual exit code) equals $2 (expected exit code)
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
echo "--- Starting tester for STB-Crystalcast (item 1.17) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Single random structure, dim=3 default (rock-salt NiO, group by number) ---
echo -e "\n--- Testing a single random dim=3 structure (--group 225) ---"

rm -f nio.fdf
stb-crystalcast --group 225 --species Ni O --num-ions 4 8 -o nio.fdf --no-intro > log_single.txt 2>&1
check_exit_code $? 0
check_contains "space group Fm-3m (No. 225)" log_single.txt
check_success nio.fdf
check_contains "NumberofAtoms      12" nio.fdf


# --- 3. Symbol form of --group (should behave the same as the number form) ---
echo -e "\n--- Testing symbol form of --group (Fm-3m) ---"

stb-crystalcast --group Fm-3m --species Ni O --num-ions 4 8 --seed 42 -o nio_symbol.fdf --no-intro > log_symbol.txt 2>&1
stb-crystalcast --group 225 --species Ni O --num-ions 4 8 --seed 42 -o nio_number.fdf --no-intro > log_number.txt 2>&1
if diff -q nio_symbol.fdf nio_number.fdf > /dev/null; then
    echo -e "   -> ${GREEN}Verified:${NC} symbol and int group forms produce identical output for the same seed"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} symbol and int group forms differ"
    FAIL=$((FAIL+1))
fi


# --- 4. --seed reproducibility ---
echo -e "\n--- Testing --seed reproducibility ---"

stb-crystalcast --group 225 --species Ni O --num-ions 4 8 --seed 7 -o repro.fdf --no-intro > /dev/null 2>&1
cp repro.fdf repro_run1.fdf
stb-crystalcast --group 225 --species Ni O --num-ions 4 8 --seed 7 -o repro.fdf --no-intro > /dev/null 2>&1
if diff -q repro.fdf repro_run1.fdf > /dev/null; then
    echo -e "   -> ${GREEN}Verified:${NC} same --seed reproduces the same structure"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} same --seed produced different structures"
    FAIL=$((FAIL+1))
fi


# --- 5. --count > 1 (numbered outputs, all distinct) ---
echo -e "\n--- Testing --count 3 (numbered outputs) ---"

rm -f batch_*.fdf
stb-crystalcast --group 225 --species Ni O --num-ions 4 8 --count 3 --seed 1 -o batch.fdf --no-intro > log_batch.txt 2>&1
check_exit_code $? 0
check_contains "3 of 3 structure(s) written" log_batch.txt
check_success batch_1.fdf
check_success batch_2.fdf
check_success batch_3.fdf
if [ ! -f batch.fdf ]; then
    echo -e "   -> ${GREEN}Verified:${NC} no unnumbered 'batch.fdf' left behind with --count > 1"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} unexpected 'batch.fdf' present alongside numbered outputs"
    FAIL=$((FAIL+1))
fi
if diff -q batch_1.fdf batch_2.fdf > /dev/null; then
    echo -e "   -> ${RED}Failed:${NC} batch_1.fdf and batch_2.fdf are identical (expected distinct random structures)"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} batch_1.fdf and batch_2.fdf differ"
    PASS=$((PASS+1))
fi


# --- 6. --dim 2 (layer group, periodic in a/b with vacuum along c) ---
echo -e "\n--- Testing --dim 2 (--group 65, layer group) ---"

rm -f layer.fdf
stb-crystalcast --dim 2 --group 65 --species C --num-ions 6 --thickness 3.4 -o layer.fdf --no-intro > log_dim2.txt 2>&1
check_exit_code $? 0
check_contains "layer group p3 (No. 65)" log_dim2.txt
check_success layer.fdf
check_contains "NumberofAtoms      6" layer.fdf


# --- 7. --dim 1 (rod group, periodic along c with vacuum in a/b) ---
echo -e "\n--- Testing --dim 1 (--group 10, rod group) ---"

rm -f rod.fdf
stb-crystalcast --dim 1 --group 10 --species C --num-ions 4 -o rod.fdf --no-intro > log_dim1.txt 2>&1
check_exit_code $? 0
check_contains "rod group p11m (No. 10)" log_dim1.txt
check_success rod.fdf


# --- 8. --dim 0 (point group, isolated cluster in a vacuum box) ---
echo -e "\n--- Testing --dim 0 (--group D3d, point group, symbol form) ---"

rm -f cluster.fdf
stb-crystalcast --dim 0 --group D3d --species C --num-ions 6 --vacuum 12 -o cluster.fdf --no-intro > log_dim0.txt 2>&1
check_exit_code $? 0
check_contains "point group D3d (No. 20)" log_dim0.txt
check_success cluster.fdf
check_contains "AtomicCoordinatesFormat  Ang" cluster.fdf

echo "Testing: --dim 0 with point group given by number (should behave the same)"
stb-crystalcast --dim 0 --group 20 --species C --num-ions 6 --vacuum 12 --seed 3 -o cluster_num.fdf --no-intro > /dev/null 2>&1
stb-crystalcast --dim 0 --group D3d --species C --num-ions 6 --vacuum 12 --seed 3 -o cluster_sym.fdf --no-intro > /dev/null 2>&1
if diff -q cluster_num.fdf cluster_sym.fdf > /dev/null; then
    echo -e "   -> ${GREEN}Verified:${NC} numeric and Schoenflies point-group forms produce identical output for the same seed"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} numeric and Schoenflies point-group forms differ"
    FAIL=$((FAIL+1))
fi


# --- 9. --molecular (packs whole rigid molecules instead of bare atoms) ---
echo -e "\n--- Testing --molecular (--group 19, --species H2O, dim=3) ---"

rm -f ice.fdf
stb-crystalcast --molecular --group 19 --species H2O --num-ions 4 -o ice.fdf --no-intro > log_molecular.txt 2>&1
check_exit_code $? 0
check_contains "Molecule:.*H2O x 4" log_molecular.txt
check_contains "space group P2_12_12_1 (No. 19)" log_molecular.txt
check_success ice.fdf
check_contains "NumberofAtoms      12" ice.fdf

echo "Testing: --molecular with --dim 2 (molecular layer)"
rm -f ice_layer.fdf
stb-crystalcast --molecular --dim 2 --group 65 --species H2O --num-ions 3 -o ice_layer.fdf --no-intro > log_molecular_dim2.txt 2>&1
check_exit_code $? 0
check_contains "layer group p3 (No. 65)" log_molecular_dim2.txt
check_success ice_layer.fdf

echo "Testing: --molecular with --dim 0 is rejected (upstream pyxtal limitation)"
stb-crystalcast --molecular --dim 0 --group D3d --species H2O --num-ions 3 --no-intro > log_molecular_dim0.txt 2>&1
check_exit_code $? 2
check_contains "does not support --dim 0" log_molecular_dim0.txt

echo "Testing: --list-molecules"
stb-crystalcast --list-molecules --no-intro > log_list_molecules.txt 2>&1
check_exit_code $? 0
check_contains "H2O" log_list_molecules.txt
check_contains "aspirin" log_list_molecules.txt

echo "Testing: unknown molecule name in --species (with a suggestion)"
stb-crystalcast --molecular --group 19 --species H2Oo --num-ions 4 --no-intro > log_bad_molecule.txt 2>&1
check_exit_code $? 1
check_contains "not a known molecule name" log_bad_molecule.txt
check_contains "Did you mean: H2O" log_bad_molecule.txt

echo "Testing: mix of a valid and an unknown molecule name is still rejected"
echo "  (regression test: pyxtal's own molecule lookup is a process-wide singleton that"
echo "   silently reuses the last successfully resolved molecule instead of raising --"
echo "   this must be caught by stb-crystalcast's own pre-validation, not left to pyxtal)"
stb-crystalcast --molecular --group 14 --species H2O BadMolecule --num-ions 4 4 --no-intro > log_bad_molecule_mix.txt 2>&1
check_exit_code $? 1
check_contains "BadMolecule.*not a known molecule name" log_bad_molecule_mix.txt

echo "Testing: nonexistent .xyz file path in --species"
stb-crystalcast --molecular --group 19 --species does_not_exist.xyz --num-ions 4 --no-intro > log_bad_molecule_file.txt 2>&1
check_exit_code $? 1
check_contains "file not found" log_bad_molecule_file.txt

echo "Testing: unsupported file extension in --species"
stb-crystalcast --molecular --group 19 --species molecule.pdb --num-ions 4 --no-intro > log_bad_molecule_ext.txt 2>&1
check_exit_code $? 1
check_contains "unsupported file extension" log_bad_molecule_ext.txt


# --- 10. --analyze (Wyckoff decomposition of an existing structure) ---
echo -e "\n--- Testing --analyze on the dim=3 NiO structure from step 2 ---"

stb-crystalcast --analyze -f nio.fdf --no-intro > log_analyze.txt 2>&1
check_exit_code $? 0
check_contains "Space group:.*Fm-3m (No. 225)" log_analyze.txt
check_contains "--site Ni" log_analyze.txt
check_contains "--site O" log_analyze.txt
check_contains "Wyckoff" log_analyze.txt

echo "Testing: --analyze without -f/--file"
stb-crystalcast --analyze --no-intro > log_analyze_nofile.txt 2>&1
check_exit_code $? 1
check_contains "requires -f/--file" log_analyze_nofile.txt

echo "Testing: --analyze with a nonexistent file"
stb-crystalcast --analyze -f does_not_exist.fdf --no-intro > log_analyze_missing.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_analyze_missing.txt


# --- 11. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: composition incompatible with the requested space group"
stb-crystalcast --group 225 --species Ni --num-ions 5 --no-intro > log_incompatible.txt 2>&1
check_exit_code $? 1
check_contains "not compatible" log_incompatible.txt

echo "Testing: --species/--num-ions length mismatch"
stb-crystalcast --group 225 --species Ni O --num-ions 4 --no-intro > log_mismatch.txt 2>&1
check_exit_code $? 1
check_contains "must match one-to-one" log_mismatch.txt

echo "Testing: invalid element symbol in --species"
stb-crystalcast --group 225 --species Xx --num-ions 4 --no-intro > log_bad_element.txt 2>&1
check_exit_code $? 1
check_contains "not a valid Element" log_bad_element.txt

echo "Testing: unrecognized space group symbol"
stb-crystalcast --group NotASpaceGroup --species Ni O --num-ions 4 8 --no-intro > log_bad_sg.txt 2>&1
check_exit_code $? 1

echo "Testing: --count 0 is rejected"
stb-crystalcast --group 225 --species Ni O --num-ions 4 8 --count 0 --no-intro > log_bad_count.txt 2>&1
check_exit_code $? 1
check_contains "count must be at least 1" log_bad_count.txt

echo "Testing: --species omitted (required unless --analyze)"
stb-crystalcast --group 225 --num-ions 4 --no-intro > log_missing_species.txt 2>&1
check_exit_code $? 2

echo "Testing: --group omitted (required unless --analyze)"
stb-crystalcast --species Ni --num-ions 4 --no-intro > log_missing_group.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-crystalcast --version > log_version.txt 2>&1
check_contains "stb-crystalcast" log_version.txt

echo "Testing: --help documents --species, --group, --dim, --analyze and --molecular"
stb-crystalcast --help > log_help.txt 2>&1
check_contains "species" log_help.txt
check_contains "group" log_help.txt
check_contains "dim" log_help.txt
check_contains "analyze" log_help.txt
check_contains "molecular" log_help.txt


# --- 12. Interactive path (stb-suite, shortcut 1.17) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 1.17) ---"

echo "Testing: navigate 1.17 -> generate (default) -> dim 3 (default) -> molecular? n (default) -> group 225 -> species 'Ni 4', 'O 8', blank to finish -> defaults -> default output -> quit"
rm -f crystalcast.fdf
printf '1.17\n\n\n\n225\nNi 4\nO 8\n\n\n\n\n\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_contains "space group Fm-3m (No. 225)" log_menu.txt
check_success crystalcast.fdf

echo "Testing: navigate 1.17 -> generate (default) -> dim 3 (default) -> molecular? y -> group 19 -> species 'H2O 4', blank to finish -> defaults -> default output -> quit"
rm -f crystalcast.fdf
printf '1.17\n\n\ny\n19\nH2O 4\n\n\n\n\n\n\n0\n' | stb-suite > log_menu_molecular.txt 2>&1
check_exit_code $? 0
check_contains "space group P2_12_12_1 (No. 19)" log_menu_molecular.txt
check_success crystalcast.fdf

echo "Testing: navigate 1.17 -> analyze -> nio.fdf -> quit"
printf '1.17\nanalyze\nnio.fdf\n\n0\n' | stb-suite > log_menu_analyze.txt 2>&1
check_exit_code $? 0
check_contains "--site Ni" log_menu_analyze.txt


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
