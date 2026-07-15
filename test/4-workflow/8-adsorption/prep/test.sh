#!/bin/bash

# --- Setup ---
# Smoke test for stb-adsorb (Adsorption Prep, item 4.8.1)
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
    if grep -q "$1" "$2" 2>/dev/null; then
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
echo "--- Starting tester for STB-Adsorb prep (item 4.8.1) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/structure.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/calc.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Single atom, ontop site (default) -- BSSE correction is ON by
#     default, so this also covers the standard (uncorrected + BSSE) case ---
echo -e "\n--- Testing a single-atom adsorbate (O), default ontop site, default BSSE ---"
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate O --no-intro > log_atom.txt 2>&1
check_contains "Success:.*1 site folder" log_atom.txt
check_success clean_slab/structure.fdf
check_success clean_slab/calc.fdf
check_success adsorbate/structure.fdf
check_contains "kgrid.MonkhorstPack   \[1  1  1\]" adsorbate/calc.fdf
check_contains "Spin.*polarized" adsorbate/calc.fdf
check_success sites/site_1_ontop/structure.fdf
check_contains "NumberofAtoms      3" sites/site_1_ontop/structure.fdf
check_success sites/adsorption_sites.txt
check_success sites/adsorption_sites.png
check_success sites/site_1_ontop/bsse_slab/structure.fdf
check_success sites/site_1_ontop/bsse_adsorbate/structure.fdf
check_contains "O_ghost" sites/site_1_ontop/bsse_slab/structure.fdf
check_contains "C_ghost" sites/site_1_ontop/bsse_adsorbate/structure.fdf
# bsse_slab keeps the real slab atoms (C) at the same positions as the site
check_contains "NumberofAtoms      3" sites/site_1_ontop/bsse_slab/structure.fdf
check_contains "NumberofAtoms      3" sites/site_1_ontop/bsse_adsorbate/structure.fdf
# persistent numbered report: [0]-[3] sections + full run metadata
check_contains "\[0\] RUN METADATA" sites/adsorption_sites.txt
check_contains "\[1\] REFERENCE FOLDERS" sites/adsorption_sites.txt
check_contains "\[2\] ADSORPTION SITES" sites/adsorption_sites.txt
check_contains "\[3\] SUMMARY" sites/adsorption_sites.txt
check_contains "Calc template   : calc.fdf" sites/adsorption_sites.txt
check_contains "BSSE correction : ON" sites/adsorption_sites.txt
check_contains "SITE_TABLE" sites/adsorption_sites.txt


# --- 2b. --no-bsse-correction skips the ghost-fragment folders entirely ---
echo -e "\n--- Testing --no-bsse-correction ---"
rm -rf clean_slab adsorbate sites
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate O --no-bsse-correction --no-intro \
    > log_no_bsse.txt 2>&1
check_success sites/site_1_ontop/structure.fdf
if [ -d sites/site_1_ontop/bsse_slab ]; then
    echo -e "   -> ${RED}Failed:${NC} 'sites/site_1_ontop/bsse_slab' unexpectedly created"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} 'sites/site_1_ontop/bsse_slab' NOT created (as expected)"
    PASS=$((PASS+1))
fi


# --- 3. Molecule adsorbate (G2), all site types, --all-sites ---
echo -e "\n--- Testing a molecule adsorbate (H2O), --site-type all --all-sites ---"
rm -rf clean_slab adsorbate sites
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate H2O --site-type all --all-sites --no-intro > log_molecule.txt 2>&1
check_contains "3 atom(s)" log_molecule.txt
n_sites=$(find sites -maxdepth 1 -type d -name 'site_*' | wc -l)
if [ "$n_sites" -ge 2 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} multiple site folders written ($n_sites)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} expected >= 2 site folders, found $n_sites"
    FAIL=$((FAIL+1))
fi


# --- 4. --ml-rank --top-k (MACE-MP-0 already cached locally) ---
echo -e "\n--- Testing --ml-rank --top-k 2 ---"
rm -rf clean_slab adsorbate sites
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate H --site-type all --all-sites \
    --ml-rank --top-k 2 --no-intro > log_mlrank.txt 2>&1
check_contains "ML pre-screen" log_mlrank.txt
check_contains "Rank" log_mlrank.txt
n_kept=$(find sites -maxdepth 1 -type d -name 'site_*' | wc -l)
if [ "$n_kept" -eq 2 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} exactly 2 site folders kept (--top-k 2)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} expected exactly 2 site folders, found $n_kept"
    FAIL=$((FAIL+1))
fi


# --- 5. --both-sides (free-standing 2D material, vacuum on both sides) ---
echo -e "\n--- Testing --both-sides ---"
rm -rf clean_slab adsorbate sites
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate H --site-type ontop --both-sides --no-intro \
    > log_bothsides.txt 2>&1
check_success sites/site_1_ontop_bothsides/structure.fdf
check_contains "NumberofAtoms      4" sites/site_1_ontop_bothsides/structure.fdf


# --- 6. Manual --position override ---
echo -e "\n--- Testing --position (manual override) ---"
rm -rf clean_slab adsorbate sites
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate O --position 0.0 0.0 --height 2.5 --no-intro \
    > log_position.txt 2>&1
check_success sites/site_1_manual/structure.fdf


# --- 6b. Multiple adsorbates in one call ---
echo -e "\n--- Testing --adsorbate O,N (multiple adsorbates) ---"
rm -rf clean_slab adsorbate adsorbate_O adsorbate_N sites
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate O,N --site-type ontop --no-bsse-correction \
    --no-intro > log_multi_ads.txt 2>&1
check_success adsorbate_O/structure.fdf
check_success adsorbate_N/structure.fdf
check_success sites/site_1_ontop_O/structure.fdf
check_success sites/site_1_ontop_N/structure.fdf
check_contains "Adsorbate: O" sites/adsorption_sites.txt
check_contains "Adsorbate: N" sites/adsorption_sites.txt
check_contains "single atom: its isolated reference is forced spin-polarized" log_multi_ads.txt

echo "Testing: an invalid name inside a comma-separated list is reported clearly"
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate O,Xx --no-intro > log_multi_bad.txt 2>&1
check_exit_code $? 1
check_contains "'Xx' is not a recognized" log_multi_bad.txt


# --- 6c. Height sweep (approach curve) ---
echo -e "\n--- Testing --height-sweep ---"
rm -rf clean_slab adsorbate sites
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate O --site-type ontop \
    --height-sweep 1.5 3.0 0.5 --no-bsse-correction --no-intro > log_height_sweep.txt 2>&1
for h in 1.50 2.00 2.50 3.00; do
    check_success "sites/site_1_ontop_h${h}/structure.fdf"
done
python3 -c "
import sys
from stb.core import structure_io
for h in (1.5, 2.0, 2.5, 3.0):
    s = structure_io.to_pymatgen(structure_io.read_fdf(f'sites/site_1_ontop_h{h:.2f}/structure.fdf'))
    z_slab = max(site.coords[2] for site in s if site.specie.symbol == 'C')
    z_ads = [site.coords[2] for site in s if site.specie.symbol == 'O'][0]
    actual = round(z_ads - z_slab, 4)
    if abs(actual - h) > 1e-6:
        print(f'height {h}: got {actual}')
        sys.exit(1)
sys.exit(0)
"
check_exit_code $? 0

echo "Testing: --height-sweep with a non-positive step is rejected"
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate O --height-sweep 1.5 3.0 0 --no-intro \
    > log_height_sweep_bad.txt 2>&1
check_exit_code $? 2


# --- 6d. --ml-prerelax (isolated adsorbate) ---
echo -e "\n--- Testing --ml-prerelax ---"
rm -rf clean_slab adsorbate sites
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate H2O --site-type ontop --ml-prerelax \
    --no-bsse-correction --no-intro > log_prerelax.txt 2>&1
check_contains "ML pre-relax" log_prerelax.txt
check_contains "Converged" log_prerelax.txt
check_success adsorbate/structure.fdf


# --- 6e. Overlap/clash warning: --height too small places the adsorbate
#     right on top of a slab atom ---
echo -e "\n--- Testing overlap warning (--position with a too-small height) ---"
rm -rf clean_slab adsorbate sites
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate O --position 0.0 0.0 --height 0.2 --no-intro \
    > log_overlap.txt 2>&1
check_contains "WARNING.*closest slab-adsorbate distance" log_overlap.txt
check_contains "WARNING.*closest slab-adsorbate distance" sites/adsorption_sites.txt


# --- 6f. Vacuum-box-vs-molecule-size warning: a molecule too large for a
#     tiny isolated-reference box may self-interact with its own images ---
echo -e "\n--- Testing --vacuum-box sanity warning (H2O in a 3.0 Ang box) ---"
rm -rf clean_slab adsorbate sites
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate H2O --vacuum-box 3.0 --no-intro \
    > log_smallbox.txt 2>&1
check_contains "WARNING.*spans.*more than half of --vacuum-box" log_smallbox.txt


# --- 7. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: invalid adsorbate name"
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate Xx --no-intro > log_bad_adsorbate.txt 2>&1
check_exit_code $? 1
check_contains "not a recognized element symbol" log_bad_adsorbate.txt

echo "Testing: vacuum axis not along c is auto-relabeled (cyclic permutation, handedness preserved)"
python3 -c "
from stb.core import structure_io
import numpy as np
s = structure_io.read_fdf('structure.fdf')
s.lattice = s.lattice[[2, 1, 0]]
s.atoms = [(sym, np.array([pos[2], pos[1], pos[0]])) for sym, pos in s.atoms]
structure_io.write_fdf(s, 'wrong_vacuum.fdf')
"
rm -rf clean_slab adsorbate sites
stb-adsorb -s wrong_vacuum.fdf -c calc.fdf --adsorbate O --no-intro > log_wrong_vacuum.txt 2>&1
check_exit_code $? 0
check_contains "vacuum axis was 'a', not c -- relabeled" log_wrong_vacuum.txt
check_success sites/site_1_ontop/structure.fdf
python3 -c "
import sys
import numpy as np
from stb.core import structure_io
orig = structure_io.to_pymatgen(structure_io.read_fdf('structure.fdf'))
relabeled = structure_io.to_pymatgen(structure_io.read_fdf('clean_slab/structure.fdf'))
orig_coords = np.array(sorted(orig.cart_coords.tolist()))
relabeled_coords = np.array(sorted(relabeled.cart_coords.tolist()))
sys.exit(0 if np.allclose(orig_coords, relabeled_coords, atol=1e-6) else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} relabeled clean_slab recovers the same Cartesian atomic positions"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} relabeled clean_slab does NOT match the original structure's Cartesian positions"
    FAIL=$((FAIL+1))
fi

echo "Testing: vacuum axis not along c is auto-relabeled (handedness-flipping permutation, corrected)"
python3 -c "
from stb.core import structure_io
import numpy as np
s = structure_io.read_fdf('structure.fdf')
s.lattice = s.lattice[[0, 2, 1]]
s.atoms = [(sym, np.array([pos[0], pos[2], pos[1]])) for sym, pos in s.atoms]
structure_io.write_fdf(s, 'vacuum_on_b.fdf')
"
rm -rf clean_slab adsorbate sites
stb-adsorb -s vacuum_on_b.fdf -c calc.fdf --adsorbate O --no-intro > log_vacuum_on_b.txt 2>&1
check_exit_code $? 0
check_contains "vacuum axis was 'b', not c -- relabeled" log_vacuum_on_b.txt
python3 -c "
import sys
from stb.core import structure_io
s = structure_io.read_fdf('clean_slab/structure.fdf')
import numpy as np
sys.exit(0 if np.linalg.det(s.lattice) > 0 else 1)
"
check_exit_code $? 0
python3 -c "
import sys
import numpy as np
from stb.core import structure_io
orig = structure_io.to_pymatgen(structure_io.read_fdf('structure.fdf'))
relabeled = structure_io.to_pymatgen(structure_io.read_fdf('clean_slab/structure.fdf'))
orig_coords = np.array(sorted(orig.cart_coords.tolist()))
relabeled_coords = np.array(sorted(relabeled.cart_coords.tolist()))
sys.exit(0 if np.allclose(orig_coords, relabeled_coords, atol=1e-6) else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} relabeled clean_slab recovers the same Cartesian atomic positions"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} relabeled clean_slab does NOT match the original structure's Cartesian positions"
    FAIL=$((FAIL+1))
fi

echo "Testing: no vacuum axis at all (bulk 3D) cannot be auto-fixed"
cat > bulk_si.fdf << 'EOF'
NumberOfSpecies    1
NumberofAtoms      2
%block ChemicalSpeciesLabel
 1   14   Si
%endblock ChemicalSpeciesLabel
LatticeConstant 1.0 Ang
AtomicCoordinatesFormat  Fractional
%block LatticeVectors
 5.43000000   0.00000000   0.00000000
 0.00000000   5.43000000   0.00000000
 0.00000000   0.00000000   5.43000000
%endblock LatticeVectors
%block AtomicCoordinatesAndAtomicSpecies
  0.000000000   0.000000000   0.000000000   1
  0.250000000   0.250000000   0.250000000   1
%endblock AtomicCoordinatesAndAtomicSpecies
EOF
stb-adsorb -s bulk_si.fdf -c calc.fdf --adsorbate O --no-intro > log_bulk.txt 2>&1
check_exit_code $? 1
check_contains "single well-defined surface" log_bulk.txt

echo "Testing: --both-sides with --site-type all"
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate O --site-type all --both-sides --no-intro \
    > log_bothsides_all.txt 2>&1
check_exit_code $? 2

echo "Testing: --top-k without --ml-rank"
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate O --all-sites --top-k 1 --no-intro \
    > log_topk_no_mlrank.txt 2>&1
check_exit_code $? 2

echo "Testing: missing structure file"
stb-adsorb -s does_not_exist.fdf -c calc.fdf --adsorbate O --no-intro > log_missing.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing.txt

echo "Testing: --version"
stb-adsorb --version > log_version.txt 2>&1
check_contains "stb-adsorb" log_version.txt

echo "Testing: --help documents --adsorbate/--site-type"
stb-adsorb --help > log_help.txt 2>&1
check_contains "adsorbate" log_help.txt
check_contains "site-type" log_help.txt

echo "Testing: --list"
stb-adsorb --list --no-intro > log_list.txt 2>&1
check_contains "H2O" log_list.txt


# --- 8. Interactive path (stb-suite, shortcut 4.8.1) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.8.1) ---"

echo "Testing: navigate 4.8.1 -> single ontop site, default BSSE (on) -> quit"
rm -rf clean_slab adsorbate sites
{
  echo "4.8.1"
  echo "structure.fdf"  # struct_file
  echo "calc.fdf"       # calc_file
  echo ""               # pp_path (skip)
  echo "O"               # adsorbate(s)
  echo ""               # ml_prerelax_choice (default N)
  echo "1"               # site_choice (ontop)
  echo ""               # height_sweep_choice (default N)
  echo ""               # height (default)
  echo "n"               # all_sites_choice -> single site
  echo "0"               # site_index
  echo "n"               # both_sides_choice
  echo ""               # bsse_choice (default -> ON)
  echo ""               # show_advanced (default -> skip)
  echo ""               # press enter to continue
  echo "0"               # quit stage submenu
} | stb-suite > log_menu.txt 2>&1
check_contains "Success:.*1 site folder" log_menu.txt
check_success sites/site_1_ontop/structure.fdf
check_success sites/site_1_ontop/bsse_slab/structure.fdf

echo "Testing: navigate 4.8.1 -> single ontop site, BSSE off -> quit"
rm -rf clean_slab adsorbate sites
{
  echo "4.8.1"
  echo "structure.fdf"
  echo "calc.fdf"
  echo ""
  echo "O"
  echo ""               # ml_prerelax_choice
  echo "1"
  echo ""               # height_sweep_choice
  echo ""
  echo "n"
  echo "0"
  echo "n"
  echo "n"               # bsse_choice -> OFF
  echo ""
  echo ""
  echo "0"
} | stb-suite > log_menu_no_bsse.txt 2>&1
check_success sites/site_1_ontop/structure.fdf
if [ -d sites/site_1_ontop/bsse_slab ]; then
    echo -e "   -> ${RED}Failed:${NC} 'sites/site_1_ontop/bsse_slab' unexpectedly created (menu BSSE=off)"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} 'sites/site_1_ontop/bsse_slab' NOT created (menu BSSE=off, as expected)"
    PASS=$((PASS+1))
fi

echo "Testing: navigate 4.8.1 -> two adsorbates + height sweep -> quit"
rm -rf clean_slab adsorbate adsorbate_* sites
{
  echo "4.8.1"
  echo "structure.fdf"
  echo "calc.fdf"
  echo ""               # pp_path
  echo "O,N"             # adsorbate(s)
  echo ""               # ml_prerelax_choice
  echo "1"               # site_choice (ontop)
  echo "y"               # height_sweep_choice -> Y
  echo "1.5"              # h_min
  echo "2.5"              # h_max
  echo "0.5"              # h_step
  echo "n"               # all_sites_choice -> single site
  echo "0"               # site_index
  echo ""               # bsse_choice (default -> ON)
  echo ""               # show_advanced
  echo ""               # press enter
  echo "0"               # quit
} | stb-suite > log_menu_sweep.txt 2>&1
check_contains "Success:.*6 site folder" log_menu_sweep.txt
check_success sites/site_1_ontop_O_h1.50/structure.fdf
check_success sites/site_1_ontop_N_h2.50/structure.fdf
check_success adsorbate_O/structure.fdf
check_success adsorbate_N/structure.fdf


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
