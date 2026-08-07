#!/bin/bash

# --- Setup ---
# Smoke test for stb-nebSites (NEB Site Selection, item 4.9.1). Reuses the
# bare 2-atom graphene slab fixture from ../../8-adsorption/prep/ (no
# adsorbate placed yet -- exactly what stb-nebSites needs, unlike
# ../prep/initial.fdf/final.fdf, which already have H placed at one
# specific position each).
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SLAB_DIR="$(cd "$FIXTURE_DIR/../../8-adsorption/prep" && pwd)"
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


# --- 1. Preparation ---
echo "--- Starting tester for STB-NebSites (item 4.9.1) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$SLAB_DIR/structure.fdf" "$TEST_DIR/"
cp "$SLAB_DIR/calc.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/structure_near_boundary.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null

RUN=nebsites_run


# --- 2. No --site-a/--site-b: prints the candidate table and stops ---
echo -e "\n--- Testing the ambiguous-selection path (no --site-a/--site-b) ---"
stb-nebSites -s structure.fdf -c calc.fdf --adsorbate H --no-intro > log_ambiguous.txt 2>&1
check_exit_code $? 1
check_contains "\[0\] RUN METADATA" log_ambiguous.txt
check_contains "\[1\] REFERENCE SLAB" log_ambiguous.txt
check_contains "\[2\] CANDIDATE SITES: FINDING & COUNT" log_ambiguous.txt
check_contains "\[3\] SITE SELECTION" log_ambiguous.txt
check_contains "Index | Type" log_ambiguous.txt
check_contains "Suggested.*closest candidate pair" log_ambiguous.txt
check_contains "STOP.*Pass --site-a/--site-b" log_ambiguous.txt
if [ -d "$RUN/site_A" ]; then
    echo -e "   -> ${RED}Failed:${NC} site_A unexpectedly written before a selection was made"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no site_A/site_B written yet (ambiguous path only lists candidates)"
    PASS=$((PASS+1))
fi
check_success "$RUN/candidate_sites.png"

echo "Testing: this bare 2-atom cell is laterally too small (same Section 3.2 check stb-adsorb uses)"
check_contains "WARNING.*own periodic image in the ab-plane" log_ambiguous.txt


# --- 3. Explicit --site-a/--site-b: writes site_A/site_B ---
echo -e "\n--- Testing --site-a 0 --site-b 3 (explicit selection) ---"
rm -rf $RUN
stb-nebSites -s structure.fdf -c calc.fdf --adsorbate H --site-a 0 --site-b 3 --save-report \
    --no-intro > log_selected.txt 2>&1
check_exit_code $? 0
check_contains "\[4\] WRITING ENDPOINT FOLDERS" log_selected.txt
check_contains "\[5\] SUMMARY & NEXT STEPS" log_selected.txt
check_contains "Chosen: site 0 (ontop) <-> site 3 (bridge)" log_selected.txt
check_success $RUN/site_A/structure.fdf
check_success $RUN/site_A/calc.fdf
check_success $RUN/site_A/config_extra.fdf
check_contains "MD.VariableCell false" $RUN/site_A/config_extra.fdf
check_contains "%include config_extra.fdf" $RUN/site_A/calc.fdf
check_success $RUN/site_B/structure.fdf
check_contains "NumberofAtoms      3" $RUN/site_A/structure.fdf
check_success $RUN/neb_sites_prep_report.txt
check_contains "\[0\] RUN METADATA" $RUN/neb_sites_prep_report.txt

echo "Testing: this fixture's slab is already centered in the vacuum gap -- no spurious recentering"
check_contains "Slab already centered in the vacuum gap" log_selected.txt

echo "Testing: --force-spin is ON by default -- Spin polarized forced via config_extra.fdf,"
echo "  overriding calc.fdf's own Spin tag (first-occurrence-wins, config_extra.fdf included first)"
check_contains "Force spin      : yes -> Spin polarized" log_selected.txt
check_contains "Spin                polarized" $RUN/site_A/config_extra.fdf
check_contains "Spin                non-polarized" calc.fdf
python3 -c "
with open('$RUN/site_A/calc.fdf') as f:
    lines = [l for l in f if l.strip()]
assert lines[0].startswith('%include config_extra.fdf'), lines[0]
print('OK')
" > log_spin_order_check.txt 2>&1
check_contains "OK" log_spin_order_check.txt

echo "Testing: --no-force-spin leaves calc.fdf's own Spin tag untouched"
rm -rf $RUN
stb-nebSites -s structure.fdf -c calc.fdf --adsorbate H --site-a 0 --site-b 3 --no-force-spin \
    --no-intro > log_nospin.txt 2>&1
check_exit_code $? 0
check_contains "Force spin      : no" log_nospin.txt
if grep -q "Spin" $RUN/site_A/config_extra.fdf; then
    echo -e "   -> ${RED}Failed:${NC} unexpected Spin tag in config_extra.fdf with --no-force-spin"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no Spin tag added to config_extra.fdf"
    PASS=$((PASS+1))
fi
rm -rf $RUN
stb-nebSites -s structure.fdf -c calc.fdf --adsorbate H --site-a 0 --site-b 3 --save-report \
    --no-intro > log_selected.txt 2>&1

echo "Testing: neb_manifest.json is written into BOTH site_A/ and site_B/ (atom-order proof for stb-neb)"
check_success $RUN/site_A/neb_manifest.json
check_success $RUN/site_B/neb_manifest.json
check_contains "atom-order proof for stb-neb" log_selected.txt
python3 -c "
import sys
from stb.core import neb_manifest
a = neb_manifest.load_manifest('$RUN/site_A/neb_manifest.json')
b = neb_manifest.load_manifest('$RUN/site_B/neb_manifest.json')
missing = [f for f in neb_manifest.REQUIRED_MANIFEST_FIELDS if f not in a or f not in b]
assert not missing, f'missing fields: {missing}'
assert a['pair_id'] == b['pair_id'], 'pair_id should match between site_A/site_B'
assert a['species_sequence'] == b['species_sequence'], 'species_sequence should match exactly'
assert len(a['species_sequence']) == a['n_atoms'] == a['n_substrate'] + a['n_adsorbate']
assert a['site_label'] == 'site_A' and b['site_label'] == 'site_B'
print('OK')
" > log_manifest_check.txt 2>&1
check_contains "OK" log_manifest_check.txt

echo "Testing: chosen_sites.png marks exactly the 2 chosen endpoints (distinct from candidate_sites.png)"
check_success $RUN/chosen_sites.png
check_contains "\[Saved\].*chosen_sites.png" log_selected.txt

echo "Testing: --view opens ASE's interactive viewer with clean_slab/site_A/site_B (headless-safe:"
echo "  MPLBACKEND=Agg forces a non-interactive matplotlib backend, DISPLAY= makes ASE's own viewer"
echo "  fail fast with a graceful [FAIL] instead of hanging -- same convention already used by"
echo "  ../../8-adsorption/prep/test.sh's own --view/--view-plots test)"
rm -rf $RUN
MPLBACKEND=Agg DISPLAY= stb-nebSites -s structure.fdf -c calc.fdf --adsorbate H --site-a 0 --site-b 3 \
    --view --view-plots --no-intro > log_view.txt 2>&1
check_contains "opening 3 frame(s)" log_view.txt
check_contains "0 = clean_slab" log_view.txt
check_contains "1 = site_A (A)" log_view.txt
check_contains "2 = site_B (B)" log_view.txt
check_contains "Could not open the interactive 3D viewer" log_view.txt
check_success $RUN/chosen_sites.png

echo "Testing: site_A and site_B share the exact same 2 substrate atoms (trivial atom "
echo "  correspondence for the NEB step that follows -- only the appended H differs)"
python3 -c "
import sys
from stb.core import structure_io
a = structure_io.read_fdf('$RUN/site_A/structure.fdf')
b = structure_io.read_fdf('$RUN/site_B/structure.fdf')
import numpy as np
sub_a = [pos for sym, pos in a.atoms[:2]]
sub_b = [pos for sym, pos in b.atoms[:2]]
sys.exit(0 if np.allclose(sub_a, sub_b, atol=1e-9) else 1)
"
check_exit_code $? 0


# --- 4. Chaining into stb-neb: --initial/--final accept these folders directly ---
echo -e "\n--- Testing the chain into stb-neb (folder input, not yet relaxed) ---"
rm -rf image_run
stb-neb --initial $RUN/site_A --final $RUN/site_B -c calc.fdf -n 5 --mode 2 -O image_run \
    --no-intro > log_neb_chain.txt 2>&1
check_exit_code $? 0
check_contains "NOT YET RELAXED" log_neb_chain.txt
check_success image_run/image_00/structure.fdf
check_success image_run/image_04/structure.fdf


# --- 4b. Slab-centering regression: a real bug fix -- a slab left near a cell
#     boundary along c (frac z=0.02 here, vacuum stacked entirely above) risks an
#     atom drifting slightly negative during SIESTA relaxation and wrapping to
#     frac z~1 (the TOP of the cell), which reads as an enormous, spurious
#     displacement when stb-neb later compares site_A against site_B.
#     center_slab_in_vacuum() must detect this and rigidly shift the whole slab
#     (all atoms, same relative geometry) to the middle of the vacuum gap before
#     any site is enumerated/written. ---
echo -e "\n--- Testing slab recentering for a near-boundary input structure ---"
rm -rf $RUN
stb-nebSites -s structure_near_boundary.fdf -c calc.fdf --adsorbate H --site-a 0 --site-b 2 \
    --no-intro > log_recenter.txt 2>&1
check_exit_code $? 0
check_contains "\[INFO\] Slab was near a cell boundary along c" log_recenter.txt
check_contains "shifted by +9.6" log_recenter.txt
check_success $RUN/site_A/structure.fdf

echo "Testing: the written site_A slab is now centered (frac z ~ 0.5), not near the boundary"
python3 -c "
from stb.core import structure_io
s = structure_io.read_fdf('$RUN/site_A/structure.fdf')
slab_z = [pos[2] for sym, pos in s.atoms[:2]]
assert all(abs(z - 0.5) < 1e-6 for z in slab_z), f'slab not centered: {slab_z}'
print('OK')
" > log_recenter_check.txt 2>&1
check_contains "OK" log_recenter_check.txt

echo "Testing: the adsorbate is still placed ABOVE the (now recentered) slab, at the requested height"
python3 -c "
from stb.core import structure_io
s = structure_io.read_fdf('$RUN/site_A/structure.fdf')
ads_z = s.atoms[-1][1][2]
# slab at frac z=0.5, cell height 20 Ang, default --height 2.0 Ang -> frac 0.1 above slab
assert abs(ads_z - 0.60) < 1e-6, f'expected adsorbate at frac z~0.60, got {ads_z}'
print('OK')
" > log_recenter_ads_check.txt 2>&1
check_contains "OK" log_recenter_ads_check.txt
rm -rf $RUN


# --- 4c. cluster_candidate_coords / nearest_candidate_pair: symmetrically-equivalent
#     candidates that land in different periodic images (pymatgen's own symm_reduce/
#     put_coord_inside wrap each candidate independently) must be re-wrapped into the
#     same periodic image, and nearest_candidate_pair must find a genuinely close pair
#     even through the periodic boundary (minimum-image, not plain Cartesian norm). ---
echo -e "\n--- Testing candidate clustering + minimum-image nearest-pair (real bug fix) ---"
python3 -c "
import numpy as np
from stb.core.adsorption_sites import cluster_candidate_coords
from stb.neb_sites import nearest_candidate_pair

lattice = np.eye(3) * 10.0
c_near = np.array([0.001, 0.5, 0.5]) @ lattice
c_far = np.array([0.999, 0.5, 0.5]) @ lattice
wrapped = cluster_candidate_coords([c_near, c_far], lattice)
w_far_frac = wrapped[1] @ np.linalg.inv(lattice)
assert abs(w_far_frac[0] - (-0.001)) < 1e-9, f'expected ~-0.001, got {w_far_frac[0]}'

# nearest_candidate_pair: a genuinely far-apart candidate (unclustered, plain Cartesian
# distance ~9.98 Ang) and a truly-close-through-the-boundary pair (~0.02 Ang via PBC) --
# the minimum-image distance must correctly pick the close pair, not the far one.
candidates = [
    (0, 'ontop', c_near),
    (1, 'ontop', c_far),               # ~0.02 Ang from c_near, THROUGH the boundary
    (2, 'ontop', np.array([0.5, 0.5, 0.5]) @ lattice),   # genuinely far from both
]
i, j, dist = nearest_candidate_pair(candidates, lattice)
assert {i, j} == {0, 1}, f'expected the near-through-boundary pair (0,1), got ({i},{j})'
assert dist < 0.1, f'expected a small minimum-image distance, got {dist}'
print('OK')
" > log_cluster_check.txt 2>&1
check_contains "OK" log_cluster_check.txt


# --- 5. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: out-of-range --site-a"
rm -rf $RUN
stb-nebSites -s structure.fdf -c calc.fdf --adsorbate H --site-a 99 --site-b 0 --no-intro \
    > log_oor.txt 2>&1
check_exit_code $? 1
check_contains "out of range" log_oor.txt

echo "Testing: --site-a equal to --site-b"
stb-nebSites -s structure.fdf -c calc.fdf --adsorbate H --site-a 1 --site-b 1 --no-intro \
    > log_same.txt 2>&1
check_exit_code $? 2

echo "Testing: only one of --site-a/--site-b given"
stb-nebSites -s structure.fdf -c calc.fdf --adsorbate H --site-a 0 --no-intro > log_onlyone.txt 2>&1
check_exit_code $? 2

echo "Testing: --site-type with too few candidates (ontop reduces to 1 site on this cell)"
rm -rf $RUN
stb-nebSites -s structure.fdf -c calc.fdf --adsorbate H --site-type ontop --no-intro \
    > log_toofew.txt 2>&1
check_exit_code $? 1
check_contains "at least 2" log_toofew.txt
echo "Testing: candidate_sites.png is still written for context, even on this error path"
check_success "$RUN/candidate_sites.png"

echo "Testing: invalid adsorbate name"
stb-nebSites -s structure.fdf -c calc.fdf --adsorbate Xx --no-intro > log_bad_ads.txt 2>&1
check_exit_code $? 1
check_contains "not a recognized" log_bad_ads.txt

echo "Testing: missing structure file"
stb-nebSites -s does_not_exist.fdf -c calc.fdf --adsorbate H --no-intro > log_missing.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing.txt

echo "Testing: --version"
stb-nebSites --version > log_version.txt 2>&1
check_contains "stb-nebSites" log_version.txt

echo "Testing: --help documents --site-a/--site-b/--adsorbate"
stb-nebSites --help > log_help.txt 2>&1
check_contains "site-a" log_help.txt
check_contains "site-b" log_help.txt
check_contains "adsorbate" log_help.txt

echo "Testing: --list"
stb-nebSites --list --no-intro > log_list.txt 2>&1
check_contains "H2O" log_list.txt


# --- 6. Interactive path (stb-suite, shortcut 4.9.1) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.9.1) ---"

echo "Testing: navigate 4.9.1 -> list candidates -> pick 0/3 -> quit"
rm -rf $RUN
{
  echo "4.9.1"
  echo "structure.fdf"  # struct_file
  echo "calc.fdf"       # calc_file
  echo ""               # pp_path (skip)
  echo "H"               # adsorbate
  echo ""               # site_choice (default: all)
  echo ""               # height (default 2.0)
  echo ""               # force_spin (default: Y)
  echo "n"               # show_advanced
  echo "0"               # site_a
  echo "3"               # site_b
  echo "y"               # save_report
  echo "n"               # view_plots
  echo ""               # press enter to continue
  echo "0"               # quit stage submenu
} | stb-suite > log_menu.txt 2>&1
check_contains "Chosen: site 0 (ontop) <-> site 3 (bridge)" log_menu.txt
check_success $RUN/site_A/structure.fdf
check_success $RUN/neb_sites_prep_report.txt


popd > /dev/null

# --- 7. Summary ---
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
