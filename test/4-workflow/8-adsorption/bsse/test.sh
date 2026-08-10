#!/bin/bash

# --- Setup ---
# Smoke test for stb-adsorbBsse (Adsorption BSSE Prep, item 4.8.2). Builds a
# real site via stb-adsorb, then fabricates a "relaxed" siesta.XV (written via
# sisl) with a geometry DIFFERENT from stb-adsorb's initial guess -- the exact
# scenario stb-adsorbBsse exists to fix: BSSE must be evaluated at the site's
# real RELAXED geometry, not the pre-relaxation guess (see examples/
# 4.8-adsorption/README.md and this stage's own [WHY THIS STAGE EXISTS] note).
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

check_absent() {
    if [ ! -e "$1" ]; then
        echo -e "   -> ${GREEN}Verified:${NC} '$1' absent, as expected"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} '$1' should not exist"
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
echo "--- Starting tester for STB-AdsorbBsse (item 4.8.2) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$PREP_DIR/structure.fdf" "$TEST_DIR/"
cp "$PREP_DIR/calc.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Build 2 sites (O ready/relaxed, N not yet relaxed) via stb-adsorb ---
echo -e "\n--- Testing BSSE generation at the site's RELAXED geometry ---"
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate O,N --site-type ontop -O . --no-intro \
    > log_prep.txt 2>&1
check_success sites/site_1_ontop_O/structure.fdf
check_success sites/site_1_ontop_N/structure.fdf

# Fake pseudopotential files, standing in for what a real stb-adsorb -p run
# would have copied into each site's own folder -- stb-adsorbBsse reuses
# these directly (no -p flag of its own), so the ghost-species dest_label
# copies (e.g. O_ghost.psml) can be verified below.
echo "fake C pseudopotential" > sites/site_1_ontop_O/C.psml
echo "fake O pseudopotential" > sites/site_1_ontop_O/O.psml
echo "fake C pseudopotential" > sites/site_1_ontop_N/C.psml
echo "fake N pseudopotential" > sites/site_1_ontop_N/N.psml

printf 'siesta: FreeEng =    -214.100000\nSCF cycle converged after 14 iterations\nsiesta: Atomic forces (eV/Ang):\n   Max    0.020000\n' \
    > sites/site_1_ontop_O/calc.out

# Fabricate a RELAXED siesta.XV for site O: moves the adsorbate (last atom,
# appended after the substrate -- see make_ghost_variant's own docstring)
# noticeably closer to the substrate than stb-adsorb's initial 2.0 Ang guess
# -- the real, physically-expected direction of relaxation for a genuine
# chemisorption bond (this fix was motivated by a live SIESTA run where an
# O-Si bond shrank from a 2.0 Ang guess to ~1.62 Ang after real relaxation).
python3 -c "
import sisl
from stb.core import structure_io

fdf = structure_io.read_fdf('sites/site_1_ontop_O/structure.fdf')
pmg = structure_io.to_pymatgen(fdf)
cart = pmg.cart_coords.copy()
cart[-1, 2] -= 0.5  # relax the adsorbate ~0.5 Ang closer to the substrate
atoms = [sisl.Atom(str(s.specie)) for s in pmg]
geom = sisl.Geometry(cart, atoms=atoms, lattice=sisl.Lattice(pmg.lattice.matrix))
sisl.get_sile('sites/site_1_ontop_O/siesta.XV', mode='w').write_geometry(geom)
with open('relaxed_z.txt', 'w') as f:
    f.write(f'{cart[-1, 2]:.6f}\n')
"
check_success sites/site_1_ontop_O/siesta.XV
check_success relaxed_z.txt

# site_1_ontop_N is deliberately left with NO siesta.XV -- simulates SIESTA
# not having finished relaxing it yet.

stb-adsorbBsse --dir . --save-report --no-intro > log_bsse.txt 2>&1
check_exit_code $? 0

echo "Testing: numbered report sections"
check_contains "\[0\] RUN METADATA" log_bsse.txt
check_contains "WHY THIS STAGE EXISTS" log_bsse.txt
check_contains "\[1\] SITE SCAN" log_bsse.txt
check_contains "\[2\] WRITING BSSE FOLDERS" log_bsse.txt
check_contains "\[3\] SUMMARY & NEXT STEPS" log_bsse.txt
check_contains "\[4\] LIBRARY WARNINGS" log_bsse.txt
check_success adsorption_bsse_report.txt

echo "Testing: site_1_ontop_O is ready, site_1_ontop_N is reported as not yet relaxed"
check_contains "site_1_ontop_O.*ready" log_bsse.txt
check_contains "site_1_ontop_N.*SKIP (no .XV yet -- not relaxed)" log_bsse.txt
check_contains "Ready (relaxed)    : 1" log_bsse.txt
check_contains "Not yet relaxed    : 1 (site_1_ontop_N)" log_bsse.txt
check_contains "BSSE folder pair(s) written : 1 of 2 site(s)" log_bsse.txt

echo "Testing: BSSE folders written only for the ready site"
check_success bsse/site_1_ontop_O/bsse_slab/structure.fdf
check_success bsse/site_1_ontop_O/bsse_adsorbate/structure.fdf
check_absent bsse/site_1_ontop_N

echo "Testing: ghost species labels (Boys-Bernardi counterpoise convention)"
check_contains "O_ghost" bsse/site_1_ontop_O/bsse_slab/structure.fdf
check_contains "C_ghost" bsse/site_1_ontop_O/bsse_adsorbate/structure.fdf

echo "Testing: BSSE config_extra.fdf forces single-point SCF on a fixed cell"
check_contains "MD.VariableCell false" bsse/site_1_ontop_O/bsse_slab/config_extra.fdf
check_contains "MD.TypeOfRun       CG" bsse/site_1_ontop_O/bsse_slab/config_extra.fdf
check_contains "MD.Steps           0" bsse/site_1_ontop_O/bsse_slab/config_extra.fdf
check_contains "%include config_extra.fdf" bsse/site_1_ontop_O/bsse_slab/calc.fdf

echo "Testing: BSSE ghost fragments inherit the site's own level of theory"
# stb-adsorb defaults to --force-spin (site) and an unconditional dipole
# correction (see write_reference_folder's force_spin/force_dipole), so
# site_1_ontop_O/config_extra.fdf carries both -- the BSSE ghost fragments
# (bsse_slab/bsse_adsorbate) must inherit BOTH, or the counterpoise
# correction mixes a spin-restriction/missing-dipole-correction penalty
# into what's supposed to be a small basis-set effect (a real, reported bug
# fixed by read_site_theory_flags()).
check_contains "Spin" sites/site_1_ontop_O/config_extra.fdf
check_contains "Slab.DipoleCorrection" sites/site_1_ontop_O/config_extra.fdf
check_contains "DFTD3" sites/site_1_ontop_O/config_extra.fdf
check_contains "Spin                polarized" bsse/site_1_ontop_O/bsse_slab/config_extra.fdf
check_contains "Slab.DipoleCorrection      .true." bsse/site_1_ontop_O/bsse_slab/config_extra.fdf
check_contains "DFTD3                   .true." bsse/site_1_ontop_O/bsse_slab/config_extra.fdf
check_contains "Spin                polarized" bsse/site_1_ontop_O/bsse_adsorbate/config_extra.fdf
check_contains "Slab.DipoleCorrection      .true." bsse/site_1_ontop_O/bsse_adsorbate/config_extra.fdf
check_contains "DFTD3                   .true." bsse/site_1_ontop_O/bsse_adsorbate/config_extra.fdf
check_contains "spin: yes, dipole: yes, vdw: yes" log_bsse.txt

echo "Testing: real + ghost pseudopotentials copied directly from the site's own folder"
check_success bsse/site_1_ontop_O/bsse_slab/C.psml
check_success bsse/site_1_ontop_O/bsse_slab/O_ghost.psml
check_success bsse/site_1_ontop_O/bsse_adsorbate/O.psml
check_success bsse/site_1_ontop_O/bsse_adsorbate/C_ghost.psml

echo "Testing: the written BSSE geometry matches the RELAXED position, not the original guess"
python3 -c "
import sys
import numpy as np
from stb.core import structure_io

# Read atoms/lattice directly (not via to_pymatgen -- ghost species labels
# like 'O_ghost' aren't valid pymatgen element symbols) and convert the
# adsorbate's own fractional position to Cartesian by hand.
relaxed_z = float(open('relaxed_z.txt').read().strip())
for sub in ('bsse_slab', 'bsse_adsorbate'):
    fdf = structure_io.read_fdf(f'bsse/site_1_ontop_O/{sub}/structure.fdf')
    frac = fdf.atoms[-1][1]
    written_z = (np.asarray(frac) @ fdf.lattice)[2]
    if abs(written_z - relaxed_z) > 1e-4:
        print(f'{sub}: expected z={relaxed_z:.6f}, got {written_z:.6f}')
        sys.exit(1)
sys.exit(0)
"
check_exit_code $? 0


# --- 2b. Negative case: a site built with --no-force-spin must NOT get a
# Spin line propagated into its BSSE ghost fragments either -- confirms
# read_site_theory_flags() actually reads the site's own config_extra.fdf
# instead of unconditionally forcing spin. ---
echo -e "\n--- Testing that --no-force-spin sites don't get a Spin line in BSSE either ---"
mkdir -p nospin && pushd nospin > /dev/null
cp "$PREP_DIR/structure.fdf" .
cp "$PREP_DIR/calc.fdf" .
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate O --site-type ontop --no-force-spin -O . --no-intro \
    > log_prep.txt 2>&1
check_success sites/site_1_ontop/structure.fdf
echo "fake C pseudopotential" > sites/site_1_ontop/C.psml
echo "fake O pseudopotential" > sites/site_1_ontop/O.psml
printf 'siesta: FreeEng =    -214.100000\nSCF cycle converged after 14 iterations\nsiesta: Atomic forces (eV/Ang):\n   Max    0.020000\n' \
    > sites/site_1_ontop/calc.out
python3 -c "
import sisl
from stb.core import structure_io

fdf = structure_io.read_fdf('sites/site_1_ontop/structure.fdf')
pmg = structure_io.to_pymatgen(fdf)
cart = pmg.cart_coords.copy()
cart[-1, 2] -= 0.5
atoms = [sisl.Atom(str(s.specie)) for s in pmg]
geom = sisl.Geometry(cart, atoms=atoms, lattice=sisl.Lattice(pmg.lattice.matrix))
sisl.get_sile('sites/site_1_ontop/siesta.XV', mode='w').write_geometry(geom)
"
check_success sites/site_1_ontop/siesta.XV
stb-adsorbBsse --dir . --save-report --no-intro > log_bsse.txt 2>&1
check_exit_code $? 0
check_contains "spin: no, dipole: yes, vdw: yes" log_bsse.txt
check_contains "DFTD3                   .true." bsse/site_1_ontop/bsse_slab/config_extra.fdf
if grep -q "^Spin" bsse/site_1_ontop/bsse_slab/config_extra.fdf 2>/dev/null; then
    echo -e "   -> ${RED}Failed:${NC} unexpected 'Spin' line in bsse_slab/config_extra.fdf (--no-force-spin site)"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no 'Spin' line in bsse_slab/config_extra.fdf, as expected"
    PASS=$((PASS+1))
fi
popd > /dev/null


# --- 3. Mismatched/stale folder (atom count doesn't match structure.fdf) ---
echo -e "\n--- Testing the atom-count-mismatch guard ---"
mkdir -p sites/site_9_fake
python3 -c "
import numpy as np
import sisl
from stb.core import structure_io

# A minimal, valid 2-atom structure.fdf (independent of the real fixture) --
# only its atom COUNT matters for this guard, not real physics.
lattice = np.eye(3) * 10.0
species_meta = {'C': {'id': '1', 'Z': 6}}
atoms = [('C', np.array([0.0, 0.0, 0.0])), ('C', np.array([0.5, 0.5, 0.5]))]
fdf = structure_io.FdfStructure(lattice=lattice, lattice_constant=1.0, species=['C'],
                                 species_meta=species_meta, atoms=atoms,
                                 coord_format='fractional', raw_lines=[])
structure_io.write_fdf(fdf, 'sites/site_9_fake/structure.fdf')

# A 3-atom .XV -- deliberately one atom too many.
geom = sisl.Geometry(np.array([[0., 0., 0.], [1., 1., 1.], [2., 2., 2.]]),
                      atoms=[sisl.Atom('C')] * 3, lattice=sisl.Lattice(lattice))
sisl.get_sile('sites/site_9_fake/siesta.XV', mode='w').write_geometry(geom)
"
cp calc.fdf sites/site_9_fake/calc.fdf
stb-adsorbBsse --dir . --save-report --no-intro > log_mismatch.txt 2>&1
check_exit_code $? 0
check_contains "mismatched/stale folder" log_mismatch.txt
check_contains "Mismatched/stale   : 1 (site_9_fake)" log_mismatch.txt
rm -rf sites/site_9_fake bsse/site_9_fake


# --- 4. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: missing 'sites' directory entirely"
mkdir -p no_sites_dir
stb-adsorbBsse --dir no_sites_dir --no-intro > log_no_sites.txt 2>&1
check_exit_code $? 1
check_contains "Did you run stb-adsorb" log_no_sites.txt

echo "Testing: missing clean_slab/structure.fdf"
mkdir -p missing_clean_slab/sites/site_1_ontop
stb-adsorbBsse --dir missing_clean_slab --no-intro > log_missing_clean_slab.txt 2>&1
check_exit_code $? 1
check_contains "not found. Did you run stb-adsorb" log_missing_clean_slab.txt

echo "Testing: no site ready at all (none relaxed yet) is a clear error, not a silent no-op"
mkdir -p none_ready/clean_slab
mkdir -p none_ready/sites/site_1_ontop
cp structure.fdf none_ready/clean_slab/structure.fdf
cp sites/site_1_ontop_N/structure.fdf none_ready/sites/site_1_ontop/structure.fdf
stb-adsorbBsse --dir none_ready --no-intro > log_none_ready.txt 2>&1
check_exit_code $? 1
check_contains "No site is ready yet" log_none_ready.txt

echo "Testing: --version"
stb-adsorbBsse --version > log_version.txt 2>&1
check_contains "stb-adsorbBsse" log_version.txt

echo "Testing: --help documents --dir/--file"
stb-adsorbBsse --help > log_help.txt 2>&1
check_contains "dir" log_help.txt
check_contains "file" log_help.txt


# --- 5. Interactive path (stb-suite, shortcut 4.8.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.8.2) ---"

echo "Testing: navigate 4.8.2 -> defaults -> quit"
rm -rf bsse
# 4.8.2 (menu code) / . (dir) / "" (out_file default) / "" (force-tolerance
# default) / "" (save_report: N) / "" (Press Enter to continue) / 0 (quit)
printf '4.8.2\n.\n\n\n\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_contains "BSSE folder pair(s) written : 1 of 2 site(s)" log_menu.txt
check_success bsse/site_1_ontop_O/bsse_slab/structure.fdf


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
