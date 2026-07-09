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


# --- 11. --lattice (fix the cell instead of estimating it) ---
echo -e "\n--- Testing --lattice ---"

rm -f nacl_lattice.fdf
stb-crystalcast --group 225 --species Na Cl --num-ions 4 4 \
    --lattice 5.64 5.64 5.64 90 90 90 -o nacl_lattice.fdf --no-intro > log_lattice.txt 2>&1
check_exit_code $? 0
check_contains "Fixed lattice:.*a=5.64" log_lattice.txt
check_success nacl_lattice.fdf
check_contains "5.64000000   0.00000000   0.00000000" nacl_lattice.fdf

echo "Testing: --lattice rejected with --dim 0"
stb-crystalcast --dim 0 --group D3d --species C --num-ions 6 \
    --lattice 20 20 20 90 90 90 --no-intro > log_lattice_dim0.txt 2>&1
check_exit_code $? 2
check_contains "not valid with --dim 0" log_lattice_dim0.txt


# --- 12. --sites (pre-assign Wyckoff positions) ---
echo -e "\n--- Testing --sites ---"

rm -f nacl_sites.fdf
stb-crystalcast --group 225 --species Na Cl --num-ions 4 4 \
    --sites 4a 4b -o nacl_sites.fdf --no-intro > log_sites.txt 2>&1
check_exit_code $? 0
check_contains "Na -> 4a" log_sites.txt
check_contains "Cl -> 4b" log_sites.txt
check_success nacl_sites.fdf

echo "Testing: --sites entry count mismatch with --species"
stb-crystalcast --group 225 --species Na Cl --num-ions 4 4 \
    --sites 4a --no-intro > log_sites_mismatch.txt 2>&1
check_exit_code $? 1
check_contains "one --sites entry is needed per --species" log_sites_mismatch.txt


# --- 13. --substitute (element substitution preserving symmetry) ---
echo -e "\n--- Testing --substitute ---"

rm -f naf.fdf
stb-crystalcast --substitute Cl:F -f nacl_lattice.fdf -o naf.fdf --no-intro > log_substitute.txt 2>&1
check_exit_code $? 0
check_contains "Substituting:.*Cl -> F" log_substitute.txt
check_contains "Output formula:.*NaF" log_substitute.txt
check_success naf.fdf
check_contains "9   F" naf.fdf

echo "Testing: --substitute with an invalid element symbol"
stb-crystalcast --substitute Xx:F -f nacl_lattice.fdf --no-intro > log_substitute_badelem.txt 2>&1
check_exit_code $? 1
check_contains "not a valid Element" log_substitute_badelem.txt

echo "Testing: --substitute with an element not present in the structure"
stb-crystalcast --substitute Br:F -f nacl_lattice.fdf --no-intro > log_substitute_missing.txt 2>&1
check_exit_code $? 1
check_contains "not present in" log_substitute_missing.txt

echo "Testing: --substitute with a malformed OLD:NEW pair"
stb-crystalcast --substitute ClF -f nacl_lattice.fdf --no-intro > log_substitute_malformed.txt 2>&1
check_exit_code $? 1
check_contains "expected OLD:NEW" log_substitute_malformed.txt

echo "Testing: --substitute without -f/--file"
stb-crystalcast --substitute Cl:F --no-intro > log_substitute_nofile.txt 2>&1
check_exit_code $? 1
check_contains "requires -f/--file" log_substitute_nofile.txt


# --- 14. --subgroup (lower-symmetry distortion search) ---
echo -e "\n--- Testing --subgroup ---"

rm -f distorted.fdf distorted_*.fdf
stb-crystalcast --subgroup -f nacl_lattice.fdf --count 3 -o distorted.fdf --no-intro > log_subgroup.txt 2>&1
check_exit_code $? 0
check_contains "Subgroup candidates found:" log_subgroup.txt
check_success distorted_1.fdf
check_success distorted_2.fdf
check_success distorted_3.fdf
check_contains "more subgroup candidates found but not written" log_subgroup.txt

echo "Testing: --subgroup with --target-group filter"
rm -f targeted.fdf targeted_*.fdf
stb-crystalcast --subgroup -f nacl_lattice.fdf --target-group 139 --count 2 \
    -o targeted.fdf --no-intro > log_subgroup_target.txt 2>&1
check_exit_code $? 0
check_contains "space group I4/mmm (No. 139)" log_subgroup_target.txt
check_success targeted_1.fdf

echo "Testing: --subgroup with no reachable candidates for an unrelated --target-group"
stb-crystalcast --subgroup -f nacl_lattice.fdf --target-group 2 --no-intro > log_subgroup_none.txt 2>&1
check_exit_code $? 1
check_contains "no subgroup candidates found" log_subgroup_none.txt

echo "Testing: --subgroup without -f/--file"
stb-crystalcast --subgroup --no-intro > log_subgroup_nofile.txt 2>&1
check_exit_code $? 1
check_contains "requires -f/--file" log_subgroup_nofile.txt


# --- 15. --supergroup (higher-symmetry parent search) ---
echo -e "\n--- Testing --supergroup ---"

echo "Testing: --supergroup requires --target-group"
stb-crystalcast --supergroup -f distorted_1.fdf --no-intro > log_supergroup_notarget.txt 2>&1
check_exit_code $? 2
check_contains "requires --target-group" log_supergroup_notarget.txt

echo "Testing: --supergroup with no candidates found (clean error, not a crash)"
stb-crystalcast --supergroup --target-group 225 -f distorted_1.fdf \
    -o parent.fdf --no-intro > log_supergroup_none.txt 2>&1
check_exit_code $? 1
check_contains "no supergroup structure found" log_supergroup_none.txt

echo "Testing: --supergroup without -f/--file"
stb-crystalcast --supergroup --target-group 225 --no-intro > log_supergroup_nofile.txt 2>&1
check_exit_code $? 1
check_contains "requires -f/--file" log_supergroup_nofile.txt


# --- 16. --ml-rank ---
echo -e "\n--- Testing --ml-rank ---"

echo "Testing: --ml-rank without the optional 'ml' extra installed (clean error, not a crash)"
stb-crystalcast --group 225 --species Ni O --num-ions 4 8 --ml-rank \
    -o mltest.fdf --no-intro > log_mlrank.txt 2>&1
check_exit_code $? 1
check_contains "mace-torch" log_mlrank.txt
if [ ! -f mltest.fdf ]; then
    echo -e "   -> ${GREEN}Verified:${NC} fails before generating anything (no mltest.fdf written)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} mltest.fdf was written despite the missing dependency"
    FAIL=$((FAIL+1))
fi

echo "Testing: --ml-rank rejected outside generation mode"
stb-crystalcast --analyze -f nio.fdf --ml-rank --no-intro > log_mlrank_wrongmode.txt 2>&1
check_exit_code $? 2
check_contains "only valid in generation mode" log_mlrank_wrongmode.txt


# --- 17. Mode mutual exclusivity ---
echo -e "\n--- Testing mode mutual exclusivity ---"

echo "Testing: --analyze and --substitute together"
stb-crystalcast --analyze --substitute Cl:F -f nio.fdf --no-intro > log_mode_clash.txt 2>&1
check_exit_code $? 2
check_contains "mutually exclusive" log_mode_clash.txt

echo "Testing: --subgroup and --supergroup together"
stb-crystalcast --subgroup --supergroup --target-group 225 -f nio.fdf --no-intro > log_mode_clash2.txt 2>&1
check_exit_code $? 2
check_contains "mutually exclusive" log_mode_clash2.txt


# --- 18. Error and robustness cases ---
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

echo "Testing: --help documents --species, --group, --dim, --analyze, --molecular, --lattice, --sites, --substitute, --subgroup, --supergroup and --ml-rank"
stb-crystalcast --help > log_help.txt 2>&1
check_contains "species" log_help.txt
check_contains "group" log_help.txt
check_contains "dim" log_help.txt
check_contains "analyze" log_help.txt
check_contains "molecular" log_help.txt
check_contains "lattice" log_help.txt
check_contains "sites" log_help.txt
check_contains "substitute" log_help.txt
check_contains "subgroup" log_help.txt
check_contains "supergroup" log_help.txt
check_contains "ml-rank" log_help.txt


# --- 19. Interactive path (stb-suite, shortcut 1.17) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 1.17) ---"

echo "Testing: navigate 1.17 -> generate (default) -> dim 3 (default) -> molecular? n -> group 225 -> species 'Ni 4', 'O 8', blank to finish -> no lattice/sites -> defaults -> ml-rank? n -> default output -> quit"
rm -f crystalcast.fdf
printf '1.17\n\n\nn\n225\nNi 4\nO 8\n\n\n\n\n\n\nn\n\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_contains "space group Fm-3m (No. 225)" log_menu.txt
check_success crystalcast.fdf

echo "Testing: navigate 1.17 -> generate (default) -> dim 3 (default) -> molecular? y -> group 19 -> species 'H2O 4', blank to finish -> no lattice/sites -> defaults -> ml-rank? n -> default output -> quit"
rm -f crystalcast.fdf
printf '1.17\n\n\ny\n19\nH2O 4\n\n\n\n\n\n\nn\n\n\n0\n' | stb-suite > log_menu_molecular.txt 2>&1
check_exit_code $? 0
check_contains "space group P2_12_12_1 (No. 19)" log_menu_molecular.txt
check_success crystalcast.fdf

echo "Testing: navigate 1.17 -> generate -> dim 3 -> molecular? n -> group 225 -> species 'Na 4', 'Cl 4' -> lattice '5.64 5.64 5.64 90 90 90' -> sites '4a 4b' -> defaults -> ml-rank? n -> output -> quit"
rm -f menu_lattice.fdf
printf '1.17\n\n\nn\n225\nNa 4\nCl 4\n\n5.64 5.64 5.64 90 90 90\n4a 4b\n\n\n\nn\nmenu_lattice.fdf\n\n0\n' | stb-suite > log_menu_lattice.txt 2>&1
check_exit_code $? 0
check_contains "Fixed lattice:.*a=5.64" log_menu_lattice.txt
check_contains "Na -> 4a" log_menu_lattice.txt
check_success menu_lattice.fdf

echo "Testing: navigate 1.17 -> analyze -> nio.fdf -> quit"
printf '1.17\nanalyze\nnio.fdf\n\n0\n' | stb-suite > log_menu_analyze.txt 2>&1
check_exit_code $? 0
check_contains "--site Ni" log_menu_analyze.txt

echo "Testing: navigate 1.17 -> substitute -> nio.fdf -> 'Ni:Co', blank to finish -> output -> quit"
printf '1.17\nsubstitute\nnio.fdf\nNi:Co\n\nmenu_sub.fdf\n\n0\n' | stb-suite > log_menu_substitute.txt 2>&1
check_exit_code $? 0
check_contains "Substituting:.*Ni -> Co" log_menu_substitute.txt
check_success menu_sub.fdf

echo "Testing: navigate 1.17 -> subgroup -> nacl_lattice.fdf -> auto target -> eps default -> count 2 -> output -> quit"
rm -f menu_subgroup.fdf menu_subgroup_*.fdf
printf '1.17\nsubgroup\nnacl_lattice.fdf\n\n\n2\nmenu_subgroup.fdf\n\n0\n' | stb-suite > log_menu_subgroup.txt 2>&1
check_exit_code $? 0
check_contains "Subgroup candidates found:" log_menu_subgroup.txt
check_success menu_subgroup_1.fdf

echo "Testing: navigate 1.17 -> supergroup -> distorted_1.fdf -> target 225 -> d-tol default -> count 1 -> output -> quit (no candidates is still a clean exit from the menu wrapper)"
printf '1.17\nsupergroup\ndistorted_1.fdf\n225\n\n1\nmenu_super.fdf\n\n0\n' | stb-suite > log_menu_supergroup.txt 2>&1
check_exit_code $? 0
check_contains "Target space group:.*225" log_menu_supergroup.txt


# --- 20. Regression tests (second review round) ---
echo -e "\n--- Testing second-round review fixes ---"

echo "Testing: plain (non-molecular) generation preserves --species order in ChemicalSpeciesLabel"
rm -f order.fdf
stb-crystalcast --group 225 --species O Ni --num-ions 8 4 --seed 1 -o order.fdf --no-intro > log_order.txt 2>&1
check_exit_code $? 0
check_success order.fdf
if grep -A3 "ChemicalSpeciesLabel" order.fdf | grep -q "1   8   O"; then
    echo -e "   -> ${GREEN}Verified:${NC} O (first --species) got id 1, matching --species order"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} species id order does not match --species order"
    FAIL=$((FAIL+1))
fi

echo "Testing: --eps nan is rejected instead of silently writing NaN coordinates"
stb-crystalcast --subgroup -f nacl_lattice.fdf --eps nan --no-intro > log_eps_nan.txt 2>&1
check_exit_code $? 1
check_contains "eps must be a positive, finite number" log_eps_nan.txt

echo "Testing: --lattice incompatible with the group's crystal system is rejected"
stb-crystalcast --group 225 --species Na Cl --num-ions 4 4 \
    --lattice 5 6 7 90 90 90 --no-intro > log_lattice_incompatible.txt 2>&1
check_exit_code $? 1
check_contains "not a valid cubic cell" log_lattice_incompatible.txt

echo "Testing: --substitute that collides with an existing species warns about the cell collapsing"
stb-crystalcast --substitute Na:Cl -f nacl_lattice.fdf -o collapsed.fdf --no-intro > log_substitute_collapse.txt 2>&1
check_exit_code $? 0
check_contains "atom count changed" log_substitute_collapse.txt

echo "Testing: --count 0 is rejected for --subgroup"
stb-crystalcast --subgroup -f nacl_lattice.fdf --count 0 --no-intro > log_subgroup_count0.txt 2>&1
check_exit_code $? 1
check_contains "count must be at least 1" log_subgroup_count0.txt

echo "Testing: --subgroup reports partial success when fewer candidates exist than --count"
rm -f fewer.fdf fewer_*.fdf
stb-crystalcast --subgroup -f nacl_lattice.fdf --target-group 139 --count 5 \
    -o fewer.fdf --no-intro > log_subgroup_partial.txt 2>&1
check_exit_code $? 1
check_contains "Partial success:" log_subgroup_partial.txt

echo "Testing: --eps/--group-type/--d-tol/--target-group warn when given outside their mode"
stb-crystalcast --group 225 --species Na Cl --num-ions 4 4 --eps 0.5 --group-type k \
    --d-tol 2.0 --target-group 99 -o silent.fdf --seed 2 --no-intro > log_ignored_flags.txt 2>&1
check_exit_code $? 0
check_contains "eps is ignored" log_ignored_flags.txt
check_contains "group-type is ignored" log_ignored_flags.txt
check_contains "d-tol is ignored" log_ignored_flags.txt
check_contains "target-group is ignored" log_ignored_flags.txt

echo "Testing: --sites pinning only zero-dof Wyckoff positions warns that --count > 1 is pointless"
stb-crystalcast --group 225 --species Na Cl --num-ions 4 4 --sites 4a 4b --count 2 \
    -o dof.fdf --no-intro > log_sites_dof.txt 2>&1
check_exit_code $? 0
check_contains "zero free parameters" log_sites_dof.txt


popd > /dev/null

# --- 22. Summary ---
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
