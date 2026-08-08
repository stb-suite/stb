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

# Every case below writes under the tool's own default --output-dir
# (adsorption_run/ -- everything for one study lives under one folder, same
# convention as stb-hubbardu's hubbardu_runs), so "rm -rf adsorption_run"
# alone resets state between cases.
RUN=adsorption_run


# --- 2. Single atom, ontop site (default) ---
echo -e "\n--- Testing a single-atom adsorbate (O), default ontop site ---"
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate O --no-intro > log_atom.txt 2>&1
check_contains "Success:.*1 site folder" log_atom.txt
check_success $RUN/clean_slab/structure.fdf
check_success $RUN/clean_slab/calc.fdf
check_success $RUN/adsorbate/structure.fdf
check_contains "kgrid.MonkhorstPack   \[1  1  1\]" $RUN/adsorbate/calc.fdf
check_contains "Spin.*polarized" $RUN/adsorbate/calc.fdf
check_success $RUN/sites/site_1_ontop/structure.fdf
check_contains "NumberofAtoms      3" $RUN/sites/site_1_ontop/structure.fdf
check_success $RUN/sites/adsorption_sites.txt
check_success $RUN/sites/adsorption_sites.png

echo "Testing: dipole correction is forced (via config_extra.fdf, not calc.fdf directly) on"
echo "  clean_slab/ and sites/site_*/ (real slabs), but NOT on adsorbate/ (isolated molecule"
echo "  in an all-around vacuum box, not a slab)"
check_contains "Slab.DipoleCorrection      .true." $RUN/clean_slab/config_extra.fdf
check_contains "Slab.DipoleCorrection      .true." $RUN/sites/site_1_ontop/config_extra.fdf
if grep -q "Slab.DipoleCorrection" $RUN/clean_slab/calc.fdf $RUN/sites/site_1_ontop/calc.fdf; then
    echo -e "   -> ${RED}Failed:${NC} Slab.DipoleCorrection should live in config_extra.fdf, not calc.fdf directly"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} calc.fdf itself carries no dipole correction tag (it's in config_extra.fdf)"
    PASS=$((PASS+1))
fi
if grep -q "Slab.DipoleCorrection" $RUN/adsorbate/config_extra.fdf $RUN/adsorbate/calc.fdf; then
    echo -e "   -> ${RED}Failed:${NC} unexpected Slab.DipoleCorrection in the isolated-adsorbate folder"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no dipole correction tag in the isolated-adsorbate reference"
    PASS=$((PASS+1))
fi

echo "Testing: DFT-D3 (van der Waals) is forced ON, unconditionally, in EVERY folder --"
echo "  clean_slab/, adsorbate/ (isolated reference), and sites/site_*/ alike (unlike the"
echo "  dipole correction, vdW is physically meaningful for a boxed molecule too)"
check_contains "DFTD3                   .true." $RUN/clean_slab/config_extra.fdf
check_contains "DFTD3                   .true." $RUN/adsorbate/config_extra.fdf
check_contains "DFTD3                   .true." $RUN/sites/site_1_ontop/config_extra.fdf
check_contains "van der Waals   : yes -> DFTD3 T" log_atom.txt
if grep -q "DFTD3" $RUN/clean_slab/calc.fdf $RUN/sites/site_1_ontop/calc.fdf $RUN/adsorbate/calc.fdf; then
    echo -e "   -> ${RED}Failed:${NC} DFTD3 should live in config_extra.fdf, not calc.fdf directly"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} calc.fdf itself carries no DFTD3 tag (it's in config_extra.fdf)"
    PASS=$((PASS+1))
fi

echo "Testing: --force-spin is ON by default -- sites/site_*/ get Spin polarized via"
echo "  config_extra.fdf (overriding calc.fdf's own 'Spin non-polarized'), clean_slab/ does not"
check_contains "Force spin      : yes -> Spin polarized" log_atom.txt
check_contains "Dipole correction: yes" log_atom.txt
check_contains "Spin                polarized" $RUN/sites/site_1_ontop/config_extra.fdf
python3 -c "
with open('$RUN/sites/site_1_ontop/calc.fdf') as f:
    lines = [l for l in f if l.strip()]
assert lines[0].startswith('%include config_extra.fdf'), lines[0]
print('OK')
" > log_spin_order_check.txt 2>&1
check_contains "OK" log_spin_order_check.txt
if grep -q "Spin" $RUN/clean_slab/config_extra.fdf; then
    echo -e "   -> ${RED}Failed:${NC} unexpected Spin tag in clean_slab/config_extra.fdf (no adsorbate there)"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} clean_slab/config_extra.fdf has no Spin override"
    PASS=$((PASS+1))
fi

echo "Testing: --no-force-spin leaves sites/site_*/ calc.fdf's own Spin tag untouched"
rm -rf $RUN
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate O --no-force-spin --no-intro > log_nospin.txt 2>&1
check_contains "Force spin      : isolated ref only" log_nospin.txt
if grep -q "Spin" $RUN/sites/site_1_ontop/config_extra.fdf; then
    echo -e "   -> ${RED}Failed:${NC} unexpected Spin tag in config_extra.fdf with --no-force-spin"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no Spin tag added to config_extra.fdf with --no-force-spin"
    PASS=$((PASS+1))
fi
check_contains "Slab.DipoleCorrection      .true." $RUN/sites/site_1_ontop/config_extra.fdf
check_contains "\[NOTE\].*'O':.*means the COMBINED" log_nospin.txt
rm -rf $RUN
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate O --no-intro > log_atom.txt 2>&1

echo "Testing: stb-adsorb no longer writes 'bsse/' at prep time -- that's stb-adsorbBsse's job, "
echo "  once sites have actually relaxed (see ../bsse/test.sh)"
if [ -d $RUN/bsse ]; then
    echo -e "   -> ${RED}Failed:${NC} 'bsse/' unexpectedly created by stb-adsorb"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} 'bsse/' NOT created by stb-adsorb (as expected)"
    PASS=$((PASS+1))
fi

# numbered report ([0]-[6], print_section style): printed to console (captured
# in log_atom.txt) always; sites/adsorption_sites.txt is now a lean machine-
# readable file (just the SITE_TABLE), NOT gated behind --save-report -- see
# the --save-report sub-case further below for the full persisted narrative.
check_contains "\[0\] RUN METADATA" log_atom.txt
check_contains "\[1\] REFERENCE FOLDERS" log_atom.txt
check_contains "\[2\] ADSORPTION SITES: FINDING & COUNT" log_atom.txt
check_contains "\[3\] ML PRE-SCREENING" log_atom.txt
check_contains "\[4\] WRITING SITE FOLDERS" log_atom.txt
check_contains "\[5\] SUMMARY & NEXT STEPS" log_atom.txt
check_contains "\[6\] LIBRARY WARNINGS" log_atom.txt
check_contains "Calc template   : calc.fdf" log_atom.txt
check_contains "Not requested" log_atom.txt
check_contains "stb-adsorbBsse next" log_atom.txt
check_contains "SITE_TABLE" $RUN/sites/adsorption_sites.txt

echo "Testing: everything lands under the default adsorption_run/ folder"
check_contains "Output dir      : adsorption_run" log_atom.txt

echo "Testing: config_extra.fdf enforces a fixed cell everywhere (single-point BSSE config_extra "
echo "  is now written by stb-adsorbBsse instead, see ../bsse/test.sh)"
check_success $RUN/clean_slab/config_extra.fdf
check_contains "MD.VariableCell false" $RUN/clean_slab/config_extra.fdf
check_contains "%include config_extra.fdf" $RUN/clean_slab/calc.fdf
check_success $RUN/sites/site_1_ontop/config_extra.fdf
check_contains "MD.VariableCell false" $RUN/sites/site_1_ontop/config_extra.fdf
if grep -q "MD.TypeOfRun" $RUN/sites/site_1_ontop/config_extra.fdf; then
    echo -e "   -> ${RED}Failed:${NC} the site's OWN config_extra.fdf should not force single-point"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} the site's own config_extra.fdf only fixes the cell (positions may relax)"
    PASS=$((PASS+1))
fi

echo "Testing: [1] warns that this small graphene primitive cell is laterally too small (adsorbate too close to its own periodic image)"
check_contains "WARNING.*own periodic image in the ab-plane" log_atom.txt

echo "Testing: --view-plots shows the generated plots (headless-safe via MPLBACKEND=Agg) without crashing"
rm -rf $RUN
MPLBACKEND=Agg stb-adsorb -s structure.fdf -c calc.fdf --adsorbate O --site-type ontop \
    --view-plots --no-intro > log_view_plots.txt 2>&1
check_exit_code $? 0
check_success $RUN/sites/adsorption_sites.png
check_contains "\[6\] LIBRARY WARNINGS" log_view_plots.txt

echo "Testing: without --save-report, sites/adsorption_sites.txt has no report narrative"
if grep -q "RUN METADATA" $RUN/sites/adsorption_sites.txt; then
    echo -e "   -> ${RED}Failed:${NC} adsorption_sites.txt unexpectedly contains the full narrative"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} adsorption_sites.txt has no report narrative (SITE_TABLE only)"
    PASS=$((PASS+1))
fi

echo "Testing: --save-report writes adsorption_run/adsorption_prep_report.txt with the full narrative"
rm -rf $RUN
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate O --save-report --no-intro > log_save_report.txt 2>&1
check_success $RUN/adsorption_prep_report.txt
check_contains "\[0\] RUN METADATA" $RUN/adsorption_prep_report.txt
check_contains "\[6\] LIBRARY WARNINGS" $RUN/adsorption_prep_report.txt
check_contains "Report          : " log_save_report.txt

echo "Testing: -s/--structure and -c/--calc default to structure.fdf/calc.fdf"
rm -rf $RUN
stb-adsorb --adsorbate O --no-intro > log_defaults.txt 2>&1
check_exit_code $? 0
check_contains "Structure       : structure.fdf" log_defaults.txt
check_contains "Calc template   : calc.fdf" log_defaults.txt
check_success $RUN/sites/site_1_ontop/structure.fdf

echo "Testing: -O/--output-dir still overrides the default (single-folder convention is opt-out, not forced)"
rm -rf $RUN custom_out
stb-adsorb --adsorbate O -O custom_out --no-intro > log_custom_out.txt 2>&1
check_success custom_out/clean_slab/structure.fdf
check_success custom_out/sites/site_1_ontop/structure.fdf
rm -rf custom_out



# --- 3. Molecule adsorbate (G2), all site types, --all-sites ---
echo -e "\n--- Testing a molecule adsorbate (H2O), --site-type all --all-sites ---"
rm -rf $RUN
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate H2O --site-type all --all-sites --no-intro > log_molecule.txt 2>&1
check_contains "3 atom(s)" log_molecule.txt
n_sites=$(find $RUN/sites -maxdepth 1 -type d -name 'site_*' | wc -l)
if [ "$n_sites" -ge 2 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} multiple site folders written ($n_sites)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} expected >= 2 site folders, found $n_sites"
    FAIL=$((FAIL+1))
fi

echo "Testing: [2] reports raw-vs-symmetry-reduced candidate counts per site type"
check_contains "Raw candidates" log_molecule.txt
check_contains "After symm. reduction" log_molecule.txt
check_contains "TOTAL" log_molecule.txt
check_contains "Configuration count :" log_molecule.txt


# --- 4. --ml-rank --top-k (MACE-MP-0 already cached locally) ---
echo -e "\n--- Testing --ml-rank --top-k 2 ---"
rm -rf $RUN
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate H --site-type all --all-sites \
    --ml-rank --top-k 2 --no-intro > log_mlrank.txt 2>&1
check_contains "Relaxing candidates for" log_mlrank.txt
check_contains "Rank" log_mlrank.txt
check_contains "\[1/4\]" log_mlrank.txt
check_contains "converged\|hit step cap" log_mlrank.txt
check_contains "top-k 2: keeping every orientation of the 2 best-ranked" log_mlrank.txt
check_success $RUN/sites/ml_rank_ranking.png

echo "Testing: the relaxed adsorbate-slab distance (not just the pre-relax height guess) is"
echo "  reported per candidate, and it's explicit that this relaxed geometry -- not the"
echo "  initial guess -- is what gets written to the site folders below"
check_contains "relaxed distance = " log_mlrank.txt
check_contains "Relaxed dist (Ang)" log_mlrank.txt
check_contains "Init\. h" log_mlrank.txt
check_contains "this relaxed geometry, not the initial guess, is what gets written" log_mlrank.txt
check_contains "\[4\] WRITING SITE FOLDERS" log_mlrank.txt
check_contains "every site folder below is written from the MACE-MP-0 -relaxed geometry" log_mlrank.txt

echo "Testing: the written structure.fdf reflects the RELAXED distance, not the 2.0 Ang initial guess"
python3 -c "
import glob
from stb.core import structure_io
site_dir = sorted(glob.glob('$RUN/sites/site_*'))[0]
s = structure_io.read_fdf(f'{site_dir}/structure.fdf')
ads_z = s.atoms[-1][1][2]
slab_z = s.atoms[0][1][2]
delta_frac = abs(ads_z - slab_z)
c_length = s.lattice[2][2]
delta_ang = delta_frac * c_length
# the initial guess was 2.0 Ang; a real MACE relax on this fixture should have
# moved it measurably away from that value (verified live: ~1.2-3.1 Ang range)
assert abs(delta_ang - 2.0) > 0.1, f'adsorbate height looks unrelaxed: {delta_ang:.3f} Ang'
print('OK')
" > log_mlrank_relaxed_check.txt 2>&1
check_contains "OK" log_mlrank_relaxed_check.txt
n_kept=$(find $RUN/sites -maxdepth 1 -type d -name 'site_*' | wc -l)
if [ "$n_kept" -eq 2 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} exactly 2 site folders kept (--top-k 2)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} expected exactly 2 site folders, found $n_kept"
    FAIL=$((FAIL+1))
fi

echo "Testing: every written site's fractional coordinates are wrapped into [0, 1)"
echo "  (AdsorbateSiteFinder places the adsorbate via a raw Cartesian offset with no"
echo "  to_unit_cell=True, so an unwrapped coordinate can otherwise land outside the cell)"
python3 -c "
import glob
from stb.core import structure_io
# Tolerance matches write_fdf's own .8f formatting: a truly-wrapped coordinate
# a hair below 1.0 (e.g. 0.9999999997, from MACE relax's own floating-point
# noise) rounds to the printed string '1.00000000' -- not a real out-of-cell
# bug, just 8-decimal display rounding at the exact 0/1 boundary.
bad = []
for site_dir in sorted(glob.glob('$RUN/sites/site_*')):
    s = structure_io.read_fdf(f'{site_dir}/structure.fdf')
    for symbol, frac in s.atoms:
        for c in frac:
            if not (-1e-6 <= c <= 1.0 + 1e-6):
                bad.append(f'{site_dir}: {symbol} has out-of-range fractional coord {c}')
if bad:
    print('\n'.join(bad))
else:
    print('OK')
" > log_frac_wrap_check.txt 2>&1
check_contains "OK" log_frac_wrap_check.txt

echo "Testing: a MACE relax that legitimately drifts the adsorbate outside the cell"
echo "  (e.g. an unstable orientation sliding across a periodic in-plane boundary) is"
echo "  wrapped back into [0, 1) right after relaxing -- BEFORE the site folder is"
echo "  written AND before orientation_trajectory.xyz (opened directly in VMD/ase view,"
echo "  not just the final SIESTA input) is exported"
python3 -c "
import sys, os, json
sys.argv = ['stb-adsorb', '-s', 'structure.fdf', '-c', 'calc.fdf', '--adsorbate', 'H2O',
            '--site-type', 'ontop', '--all-sites', '--height', '2.0', '--ml-rank',
            '--n-orientations-polar', '2', '--n-orientations-azimuthal', '1',
            '--no-intro', '-O', 'wrap_after_relax_run']

from stb.core import mace_relax
_orig_relax = mace_relax.relax
pre_shift_frac = {}  # orient_j (1-based call order) -> unwrapped frac coords, before the
                      # synthetic shift below -- the real MACE-relaxed position, whatever
                      # it physically ended up being (not assumed to be any specific value).
_call_counter = [0]

def _patched_relax(atoms, calc, cell_mask=None, optimizer='FIRE', fmax=0.05, max_steps=200,
                    step_history=None):
    converged, steps = _orig_relax(atoms, calc, cell_mask=cell_mask, optimizer=optimizer,
                                    fmax=fmax, max_steps=max_steps, step_history=step_history)
    _call_counter[0] += 1
    pre_shift_frac[_call_counter[0]] = atoms.get_scaled_positions(wrap=False)[2:].copy()
    # Simulate the real-world case this fix targets: the optimizer's own
    # trajectory legitimately carries the free adsorbate atom(s) a full+half
    # lattice vector past the cell boundary along both in-plane axes -- e.g.
    # an unstable starting orientation sliding across x=0/x=1 while relaxing.
    # Shifts EVERY free (non-substrate) atom, not just one -- H2O has 3.
    cell = atoms.get_cell()
    atoms.positions[2:] += 1.5 * cell[0] + 1.5 * cell[1]
    return converged, steps

mace_relax.relax = _patched_relax
from stb.adsorb import main
main()
with open('wrap_after_relax_run/pre_shift_frac.json', 'w') as f:
    json.dump({k: v.tolist() for k, v in pre_shift_frac.items()}, f)
" > log_wrap_after_relax.txt 2>&1
check_exit_code $? 0

echo "Testing: orientation_trajectory.xyz frames are wrapped correctly -- the same"
echo "  physical position modulo the injected 1.5-lattice-vector shift (opened directly"
echo "  in VMD/ase view, not just the final SIESTA input written in [4])"
python3 -c "
import json
import numpy as np
import ase.io
with open('wrap_after_relax_run/pre_shift_frac.json') as f:
    pre_shift_frac = {int(k): np.array(v) for k, v in json.load(f).items()}
frames = ase.io.read('wrap_after_relax_run/sites/orientation_trajectory.xyz', index=':')
assert len(frames) == len(pre_shift_frac), \
    f'expected {len(pre_shift_frac)} sampled orientation frame(s), found {len(frames)}'
for atoms in frames:
    orient_j = atoms.info['orientation']
    wrapped = atoms.get_scaled_positions(wrap=False)[2:]
    assert np.all((wrapped >= -1e-6) & (wrapped <= 1.0 + 1e-6)), \
        f'orientation {orient_j}: out-of-range fractional coord(s): {wrapped}'
    # Wrapped value must equal (pre-shift + 1.5 along a/b only, z untouched)
    # mod 1, whatever the real relaxed position actually was -- proves the
    # SAME fractional part survived (only the whole-lattice-vector part was
    # removed), not silently clamped/discarded/corrupted by the wrap.
    expected = pre_shift_frac[orient_j].copy()
    expected[:, :2] += 1.5
    expected %= 1.0
    delta = np.abs(wrapped - expected)
    delta = np.minimum(delta, 1.0 - delta)  # tolerate a wrap-boundary flip (e.g. -1e-9 vs 1-1e-9)
    assert np.all(delta < 1e-4), \
        f'orientation {orient_j}: wrapped {wrapped} != expected {expected} (pre-shift {pre_shift_frac[orient_j]})'
print('OK')
" > log_wrap_xyz_check.txt 2>&1
check_contains "OK" log_wrap_xyz_check.txt

echo "Testing: the same wrapped positions are what's actually written to the site folders"
python3 -c "
import glob
from stb.core import structure_io
bad = []
for site_dir in sorted(glob.glob('wrap_after_relax_run/sites/site_*')):
    s = structure_io.read_fdf(f'{site_dir}/structure.fdf')
    for symbol, frac in s.atoms:
        if not all(-1e-6 <= c <= 1.0 + 1e-6 for c in frac):
            bad.append(f'{site_dir}: {symbol} has out-of-range fractional coord {frac}')
assert not bad, '\n'.join(bad)
print('OK')
" > log_wrap_after_relax_check.txt 2>&1
check_contains "OK" log_wrap_after_relax_check.txt
rm -rf wrap_after_relax_run


# --- 4a. generate_systematic_orientations / deduplicate_orientations: pure
# unit tests (no MACE, no SIESTA) on the two new core/adsorption_sites.py
# functions, before exercising them live through --ml-rank below. ---
echo -e "\n--- Testing generate_systematic_orientations / deduplicate_orientations (unit) ---"
python3 -c "
import numpy as np
from pymatgen.core import Molecule
from ase import Atoms
from stb.core.adsorption_sites import generate_systematic_orientations, deduplicate_orientations

mol = Molecule(['O', 'H', 'H'], [[0, 0, 0], [0.76, 0.59, 0], [-0.76, 0.59, 0]])

# default (1x1) is a no-op -- exact same molecule, unchanged
out = generate_systematic_orientations(mol, 1, 1)
assert len(out) == 1
assert np.allclose(out[0].cart_coords, mol.cart_coords)

# 6x4 = 24 distinct orientations, all rigid rotations (bond lengths/COM preserved)
out = generate_systematic_orientations(mol, 6, 4)
assert len(out) == 24
coords = [o.cart_coords for o in out]
for i in range(len(coords)):
    for j in range(i + 1, len(coords)):
        assert not np.allclose(coords[i], coords[j], atol=1e-6), f'orientations {i},{j} identical'
com0 = mol.center_of_mass
d0 = np.linalg.norm(mol.cart_coords[0] - mol.cart_coords[1])
for o in out:
    assert np.allclose(o.center_of_mass, com0, atol=1e-6)
    assert abs(np.linalg.norm(o.cart_coords[0] - o.cart_coords[1]) - d0) < 1e-6

# deduplicate_orientations: needs BOTH close RMSD and close energy to collapse
n_sub = 2
def entry(ads_pos, e):
    return e, Atoms(['C', 'C'] + ['O'] * len(ads_pos), positions=[[0, 0, 0], [1, 0, 0]] + ads_pos)

near_geom_near_e = [entry([[0, 0, 2.0]], -10.000), entry([[0.01, 0, 2.0]], -10.001)]
assert deduplicate_orientations(near_geom_near_e, n_sub, rmsd_tol=0.3, energy_tol=0.01) == [0]

near_geom_far_e = [entry([[0, 0, 2.0]], -10.000), entry([[0.02, 0, 2.0]], -8.000)]
assert deduplicate_orientations(near_geom_far_e, n_sub, rmsd_tol=0.3, energy_tol=0.01) == [0, 1]

far_geom_near_e = [entry([[0, 0, 2.0]], -10.000), entry([[3.0, 3.0, 2.0]], -10.0005)]
assert deduplicate_orientations(far_geom_near_e, n_sub, rmsd_tol=0.3, energy_tol=0.01) == [0, 1]

print('OK')
" > log_orientation_unit_check.txt 2>&1
check_contains "OK" log_orientation_unit_check.txt


# --- 4a2. --n-orientations-polar/--n-orientations-azimuthal live, through
# --ml-rank, with a real (multi-atom) H2O adsorbate on the ontop site (MACE
# -MP-0 already cached locally from section 4 above). Kept small (2x2=4) for
# test runtime. ---
echo -e "\n--- Testing --n-orientations-polar/--n-orientations-azimuthal (live, H2O) ---"
rm -rf $RUN
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate H2O --site-type ontop --all-sites \
    --ml-rank --n-orientations-polar 2 --n-orientations-azimuthal 2 --no-intro \
    > log_orientations.txt 2>&1
check_exit_code $? 0
check_contains "4 orientation(s) sampled" log_orientations.txt
check_contains "Orientation" log_orientations.txt
n_site_dirs=$(find $RUN/sites -maxdepth 1 -type d -name 'site_1_ontop_orient*' | wc -l)
if [ "$n_site_dirs" -ge 1 ] && [ "$n_site_dirs" -le 4 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} 1-4 unique orientation folder(s) written ($n_site_dirs)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} expected 1-4 orientation folders, found $n_site_dirs"
    FAIL=$((FAIL+1))
fi
check_success $RUN/sites/site_1_ontop_orient1/structure.fdf
check_contains "NumberofAtoms      5" $RUN/sites/site_1_ontop_orient1/structure.fdf

echo "Testing: succinct per-orientation progress counter is printed"
check_contains "orientation 1/4 relaxed (E = " log_orientations.txt
check_contains "orientation 4/4 relaxed (E = " log_orientations.txt

echo "Testing: orientation_trajectory.xyz carries every generated orientation (not just kept ones)"
check_success $RUN/sites/orientation_trajectory.xyz
n_frames=$(grep -c "^Lattice=" $RUN/sites/orientation_trajectory.xyz 2>/dev/null)
if [ "$n_frames" -eq 4 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} orientation_trajectory.xyz has all 4 generated frames"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} expected 4 frames in orientation_trajectory.xyz, found $n_frames"
    FAIL=$((FAIL+1))
fi
check_contains "site=1 site_type=ontop" $RUN/sites/orientation_trajectory.xyz
check_contains "energy_eV=" $RUN/sites/orientation_trajectory.xyz

echo "Testing: --ml-rank via orientation sampling on a single --site-index (no --all-sites)"
rm -rf $RUN
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate H2O --site-type ontop --site-index 0 \
    --ml-rank --n-orientations-polar 2 --n-orientations-azimuthal 2 --no-intro \
    > log_orientations_singlesite.txt 2>&1
check_exit_code $? 0
check_contains "site 1/1" log_orientations_singlesite.txt
check_success $RUN/sites/site_1_ontop_orient1/structure.fdf
check_success $RUN/sites/orientation_trajectory.xyz

echo "Testing: --ml-rank without --all-sites and without orientation sampling is still rejected"
rm -rf $RUN
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate H2O --site-type ontop --site-index 0 \
    --ml-rank --no-intro > log_mlrank_singlesite_rejected.txt 2>&1
check_exit_code $? 2
check_contains "UNLESS orientation sampling" log_mlrank_singlesite_rejected.txt

echo "Testing: --orientation-top-k caps the number of kept orientations"
rm -rf $RUN
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate H2O --site-type ontop --all-sites \
    --ml-rank --n-orientations-polar 2 --n-orientations-azimuthal 2 --orientation-top-k 1 \
    --no-intro > log_orientations_topk.txt 2>&1
check_exit_code $? 0
n_kept=$(find $RUN/sites -maxdepth 1 -type d -name 'site_1_ontop_orient*' | wc -l)
if [ "$n_kept" -eq 1 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} --orientation-top-k 1 kept exactly 1 orientation"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} expected exactly 1 orientation folder, found $n_kept"
    FAIL=$((FAIL+1))
fi

echo "Testing: orientation flags without --ml-rank are rejected"
rm -rf $RUN
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate H2O --site-type ontop \
    --n-orientations-polar 2 --no-intro > log_orientations_no_mlrank.txt 2>&1
check_exit_code $? 2
check_contains "only valid with --ml-rank" log_orientations_no_mlrank.txt


# --- 4b. --ml-device cuda without a usable GPU: warn and fall back to CPU,
# not abort (this machine has no GPU, so this exercises the real fallback
# path rather than a mocked one). ---
echo -e "\n--- Testing --ml-device cuda falls back to CPU with a warning (no GPU here) ---"
rm -rf $RUN
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate H --site-type ontop --all-sites \
    --ml-rank --top-k 1 --ml-device cuda --no-intro > log_ml_device_fallback.txt 2>&1
check_exit_code $? 0
check_contains "\[WARNING\] --ml-device cuda requested, but .* -- falling back to CPU\." log_ml_device_fallback.txt
check_success $RUN/sites/site_1_ontop/structure.fdf


# --- 4c. A --height large enough to push the adsorbate's fractional z
# decisively past 1.0 (not just floating-point noise at the 0/1 boundary --
# see the [4] wrap check above): AdsorbateSiteFinder.add_adsorbate places
# the adsorbate via a raw Cartesian offset with no to_unit_cell=True, so
# this would land at frac z=1.25 (25 Ang / 20 Ang cell) without the
# wrap_into_cell() fix in adsorb.py's site-writing loop. ---
echo -e "\n--- Testing a --height that pushes the adsorbate decisively past frac z=1.0 ---"
rm -rf $RUN
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate H --site-type ontop --height 15 --no-intro \
    > log_height_wrap.txt 2>&1
check_exit_code $? 0
python3 -c "
from stb.core import structure_io
s = structure_io.read_fdf('$RUN/sites/site_1_ontop/structure.fdf')
symbol, frac = s.atoms[-1]
z = frac[2]
# unwrapped would be 25/20 = 1.25; wrapped, it's 0.25
assert 0.0 <= z < 1.0, f'{symbol} z={z} is not in [0, 1)'
assert abs(z - 0.25) < 1e-6, f'{symbol} z={z}, expected ~0.25 (1.25 wrapped)'
print('OK')
" > log_height_wrap_check.txt 2>&1
check_contains "OK" log_height_wrap_check.txt


# --- 5. --both-sides (free-standing 2D material, vacuum on both sides) ---
echo -e "\n--- Testing --both-sides ---"
rm -rf $RUN
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate H --site-type ontop --both-sides --no-intro \
    > log_bothsides.txt 2>&1
check_success $RUN/sites/site_1_ontop_bothsides/structure.fdf
check_contains "NumberofAtoms      4" $RUN/sites/site_1_ontop_bothsides/structure.fdf


# --- 6. Manual --position override ---
echo -e "\n--- Testing --position (manual override) ---"
rm -rf $RUN
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate O --position 0.0 0.0 --height 2.5 --no-intro \
    > log_position.txt 2>&1
check_success $RUN/sites/site_1_manual/structure.fdf


# --- 6b. Multiple adsorbates in one call ---
echo -e "\n--- Testing --adsorbate O,N (multiple adsorbates) ---"
rm -rf $RUN
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate O,N --site-type ontop \
    --no-intro > log_multi_ads.txt 2>&1
check_success $RUN/adsorbate_O/structure.fdf
check_success $RUN/adsorbate_N/structure.fdf
check_success $RUN/sites/site_1_ontop_O/structure.fdf
check_success $RUN/sites/site_1_ontop_N/structure.fdf
# "Adsorbate:" itself is ANSI-color-wrapped in the raw console log (unlike
# the old adsorption_sites.txt, print_dual doesn't strip color codes going to
# stdout), so match the plain text on either side of the color reset instead.
check_contains "O (1 atom(s))" log_multi_ads.txt
check_contains "N (1 atom(s))" log_multi_ads.txt
check_contains "'O': its isolated reference is forced spin-polarized" log_multi_ads.txt
check_contains "'N': its isolated reference is forced spin-polarized" log_multi_ads.txt

echo "Testing: an invalid name inside a comma-separated list is reported clearly"
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate O,Xx --no-intro > log_multi_bad.txt 2>&1
check_exit_code $? 1
check_contains "'Xx' is not a recognized" log_multi_bad.txt


# --- 6c. Height sweep (approach curve) ---
echo -e "\n--- Testing --height-sweep ---"
rm -rf $RUN
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate O --site-type ontop \
    --height-sweep 1.5 3.0 0.5 --no-intro > log_height_sweep.txt 2>&1
for h in 1.50 2.00 2.50 3.00; do
    check_success "$RUN/sites/site_1_ontop_h${h}/structure.fdf"
done
python3 -c "
import sys
from stb.core import structure_io
for h in (1.5, 2.0, 2.5, 3.0):
    s = structure_io.to_pymatgen(structure_io.read_fdf(f'adsorption_run/sites/site_1_ontop_h{h:.2f}/structure.fdf'))
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
rm -rf $RUN
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate H2O --site-type ontop --ml-prerelax \
    --no-intro > log_prerelax.txt 2>&1
check_contains "ML pre-relax" log_prerelax.txt
check_contains "Converged" log_prerelax.txt
check_success $RUN/adsorbate/structure.fdf


# --- 6e. Overlap/clash warning: --height too small places the adsorbate
#     right on top of a slab atom ---
echo -e "\n--- Testing overlap warning (--position with a too-small height) ---"
rm -rf $RUN
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate O --position 0.0 0.0 --height 0.2 --no-intro \
    > log_overlap.txt 2>&1
check_contains "WARNING.*closest slab-adsorbate distance" log_overlap.txt


# --- 6f. Vacuum-box-vs-molecule-size warning: a molecule too large for a
#     tiny isolated-reference box may self-interact with its own images ---
echo -e "\n--- Testing --vacuum-box sanity warning (H2O in a 3.0 Ang box) ---"
rm -rf $RUN
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate H2O --vacuum-box 3.0 --no-intro \
    > log_smallbox.txt 2>&1
check_contains "WARNING.*spans.*more than half of --vacuum-box" log_smallbox.txt


# --- 6g. cluster_candidate_coords: symmetrically-equivalent candidates that land in
#     different periodic images (pymatgen's own AdsorbateSiteFinder.symm_reduce/
#     put_coord_inside wrap each candidate independently, no shared reference point)
#     must be re-wrapped into the SAME periodic image -- a real bug reported live: a
#     user's adsorption_sites.png showed most candidates clustered together but 1-2
#     appearing far away, even though they were the same physical site type. ---
echo -e "\n--- Testing cluster_candidate_coords (symmetrically-equivalent site clustering) ---"
python3 -c "
import numpy as np
from stb.core.adsorption_sites import cluster_candidate_coords

lattice = np.eye(3) * 10.0
# two physically-equivalent points just on either side of a cell boundary
c_near = np.array([0.001, 0.5, 0.5]) @ lattice
c_far = np.array([0.999, 0.5, 0.5]) @ lattice   # pymatgen's own floor-wrap would leave this here
wrapped = cluster_candidate_coords([c_near, c_far], lattice)
w_far_frac = wrapped[1] @ np.linalg.inv(lattice)
assert abs(w_far_frac[0] - (-0.001)) < 1e-9, f'expected ~-0.001, got {w_far_frac[0]}'
assert abs(w_far_frac[1] - 0.5) < 1e-9 and abs(w_far_frac[2] - 0.5) < 1e-9

# an explicit reference point (not coords[0]) is honored
wrapped2 = cluster_candidate_coords([c_far], lattice, reference=c_near)
w2_frac = wrapped2[0] @ np.linalg.inv(lattice)
assert abs(w2_frac[0] - (-0.001)) < 1e-9

# empty input -> unchanged (no crash)
assert cluster_candidate_coords([], lattice) == []
print('OK')
" > log_cluster_check.txt 2>&1
check_contains "OK" log_cluster_check.txt


# --- 6h. wrap_markers_into_cell: adsorption_sites.png must mark ONLY the
#     candidates' representative position INSIDE the highlighted unit cell
#     (per explicit request), one marker per candidate, no periodic clutter
#     outside it. ---
echo -e "\n--- Testing wrap_markers_into_cell (every candidate marked exactly once, inside the cell) ---"
python3 -c "
import numpy as np
from pymatgen.analysis.adsorption import get_rot
from pymatgen.core import Structure, Lattice
from stb.core.adsorption_sites import wrap_markers_into_cell

lattice = Lattice.from_parameters(2.46, 2.46, 20.0, 90, 90, 120)
structure = Structure(lattice, ['C', 'C'], [[0, 0, 0.5], [1/3, 2/3, 0.5]])
symm_op = get_rot(structure)
a2d = np.asarray(symm_op.operate(lattice.matrix[0])[:2])
b2d = np.asarray(symm_op.operate(lattice.matrix[1])[:2])
inv_ab = np.linalg.inv(np.column_stack([a2d, b2d]))

def frac_ab(xy):
    return inv_ab @ np.asarray(xy)

# a candidate deliberately outside the [0,1) cell in fractional terms
# (frac x=-0.3, physically the same as x=0.7) -- must come back wrapped
# into [0,1) x [0,1), one point, no extra copies.
cart_outside = np.array([-0.3, 0.5, 0.5]) @ lattice.matrix
# a candidate already inside [0,1) -- must come back unchanged (up to fp noise)
cart_inside = np.array([0.2, 0.4, 0.5]) @ lattice.matrix

wrapped = wrap_markers_into_cell([cart_outside, cart_inside], lattice.matrix, symm_op)
assert len(wrapped) == 2, f'expected exactly 1 marker per candidate, got {len(wrapped)}'
for xy in wrapped:
    alpha, beta = frac_ab(xy)
    assert -1e-9 <= alpha < 1.0 + 1e-9 and -1e-9 <= beta < 1.0 + 1e-9, \
        f'marker not inside the cell: alpha={alpha}, beta={beta}'

alpha0, beta0 = frac_ab(wrapped[0])
assert abs(alpha0 - 0.7) < 1e-6, f'expected alpha~0.7 (wrapped from -0.3), got {alpha0}'
alpha1, beta1 = frac_ab(wrapped[1])
assert abs(alpha1 - 0.2) < 1e-6 and abs(beta1 - 0.4) < 1e-6

# empty input -> empty output (no crash)
assert wrap_markers_into_cell([], lattice.matrix, symm_op) == []
print('OK')
" > log_periodic_marker_check.txt 2>&1
check_contains "OK" log_periodic_marker_check.txt


# --- 6i. center_slab_in_vacuum: a slab pinned near a cell boundary along c
#     (uneven vacuum split) must be recentered before any site is written --
#     user-requested, same fix already shipped for stb-nebSites, now shared
#     via core/adsorption_sites.py once stb-adsorb became a second consumer.
#     Deliberately build a structure with frac z=0.02 in a 20 Ang cell (same
#     scenario stb-nebSites' own live verification used). ---
echo -e "\n--- Testing center_slab_in_vacuum (slab near a cell boundary gets recentered) ---"
python3 -c "
import numpy as np
from stb.core import structure_io

lattice = np.array([[2.46, 0.0, 0.0], [-1.23, 2.130422, 0.0], [0.0, 0.0, 20.0]])
species_meta = {'C': {'id': '1', 'Z': 6}}
atoms = [('C', np.array([0.0, 0.0, 0.02])), ('C', np.array([1 / 3, 2 / 3, 0.02]))]
s = structure_io.FdfStructure(lattice=lattice, lattice_constant=1.0, species=['C'],
                               species_meta=species_meta, atoms=atoms,
                               coord_format='fractional', raw_lines=[])
structure_io.write_fdf(s, 'structure_offcenter.fdf')
"
check_success structure_offcenter.fdf
rm -rf $RUN
stb-adsorb -s structure_offcenter.fdf -c calc.fdf --adsorbate H --site-type ontop --no-intro \
    > log_offcenter.txt 2>&1
check_exit_code $? 0
check_contains "\[INFO\] Slab was near a cell boundary along c" log_offcenter.txt
check_contains "shifted by +9.600 Ang to center it" log_offcenter.txt
python3 -c "
from stb.core import structure_io
s = structure_io.read_fdf('$RUN/clean_slab/structure.fdf')
for symbol, frac in s.atoms:
    assert abs(frac[2] - 0.5) < 1e-6, f'{symbol} z={frac[2]}, expected 0.5 (recentered)'
print('OK')
" > log_offcenter_check.txt 2>&1
check_contains "OK" log_offcenter_check.txt

echo "Testing: an already-centered slab is left unchanged (no gratuitous rewrite)"
rm -rf $RUN
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate H --site-type ontop --no-intro \
    > log_centered.txt 2>&1
check_contains "Slab already centered in the vacuum gap" log_centered.txt


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
rm -rf $RUN
stb-adsorb -s wrong_vacuum.fdf -c calc.fdf --adsorbate O --no-intro > log_wrong_vacuum.txt 2>&1
check_exit_code $? 0
check_contains "vacuum axis was 'a', not c -- relabeled" log_wrong_vacuum.txt
check_success $RUN/sites/site_1_ontop/structure.fdf
python3 -c "
import sys
import numpy as np
from stb.core import structure_io
orig = structure_io.to_pymatgen(structure_io.read_fdf('structure.fdf'))
relabeled = structure_io.to_pymatgen(structure_io.read_fdf('adsorption_run/clean_slab/structure.fdf'))
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
rm -rf $RUN
stb-adsorb -s vacuum_on_b.fdf -c calc.fdf --adsorbate O --no-intro > log_vacuum_on_b.txt 2>&1
check_exit_code $? 0
check_contains "vacuum axis was 'b', not c -- relabeled" log_vacuum_on_b.txt
python3 -c "
import sys
from stb.core import structure_io
s = structure_io.read_fdf('adsorption_run/clean_slab/structure.fdf')
import numpy as np
sys.exit(0 if np.linalg.det(s.lattice) > 0 else 1)
"
check_exit_code $? 0
python3 -c "
import sys
import numpy as np
from stb.core import structure_io
orig = structure_io.to_pymatgen(structure_io.read_fdf('structure.fdf'))
relabeled = structure_io.to_pymatgen(structure_io.read_fdf('adsorption_run/clean_slab/structure.fdf'))
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
check_contains "view-plots" log_help.txt
check_contains "adsorption_run" log_help.txt

echo "Testing: --list"
stb-adsorb --list --no-intro > log_list.txt 2>&1
check_contains "H2O" log_list.txt


# --- 8. Interactive path (stb-suite, shortcut 4.8.1) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.8.1) ---"

echo "Testing: navigate 4.8.1 -> single ontop site -> quit"
rm -rf $RUN
{
  echo "4.8.1"
  echo "structure.fdf"  # struct_file
  echo "calc.fdf"       # calc_file
  echo ""               # pp_path (skip)
  echo "O"               # adsorbate(s)
  echo ""               # ml_prerelax (default N)
  echo "1"               # site_choice (ontop)
  echo ""               # height (blank = default 2.0, no sweep)
  echo "n"               # all_sites_choice -> single site
  echo "0"               # site_index
  echo ""               # force_spin (default Y)
  echo ""               # orient_grid (blank -> skip orientation sampling)
  echo "n"               # both_sides_choice
  echo ""               # show_advanced (default -> skip, so output_dir defaults to adsorption_run)
  echo "y"               # save_report -> yes
  echo "n"               # view_choice -> no (headless)
  echo "n"               # view_plots_choice -> no (headless)
  echo ""               # press enter to continue
  echo "0"               # quit stage submenu
} | stb-suite > log_menu.txt 2>&1
check_contains "Success:.*1 site folder" log_menu.txt
check_success $RUN/sites/site_1_ontop/structure.fdf
check_success $RUN/adsorption_prep_report.txt
check_contains "\[0\] RUN METADATA" $RUN/adsorption_prep_report.txt

echo "Testing: navigate 4.8.1 -> two adsorbates + height sweep -> quit"
rm -rf $RUN
{
  echo "4.8.1"
  echo "structure.fdf"
  echo "calc.fdf"
  echo ""               # pp_path
  echo "O,N"             # adsorbate(s)
  echo ""               # ml_prerelax
  echo "1"               # site_choice (ontop)
  echo "sweep"           # height -> trigger the approach-curve scan
  echo "1.5"              # h_min
  echo "2.5"              # h_max
  echo "0.5"              # h_step
  echo "n"               # all_sites_choice -> single site
  echo "0"               # site_index
  echo ""               # force_spin (default Y)
  echo ""               # orient_grid (blank -> skip orientation sampling)
  echo ""               # show_advanced
  echo "n"               # save_report -> no
  echo "n"               # view_choice -> no
  echo "n"               # view_plots_choice -> no
  echo ""               # press enter
  echo "0"               # quit
} | stb-suite > log_menu_sweep.txt 2>&1
check_contains "Success:.*6 site folder" log_menu_sweep.txt
check_success $RUN/sites/site_1_ontop_O_h1.50/structure.fdf
check_success $RUN/sites/site_1_ontop_N_h2.50/structure.fdf
check_success $RUN/adsorbate_O/structure.fdf
check_success $RUN/adsorbate_N/structure.fdf

echo "Testing: navigate 4.8.1 -> --ml-rank + orientation sampling (H2O) -> quit"
rm -rf $RUN
{
  echo "4.8.1"
  echo "structure.fdf"
  echo "calc.fdf"
  echo ""               # pp_path
  echo "H2O"             # adsorbate(s)
  echo ""               # ml_prerelax
  echo "1"               # site_choice (ontop)
  echo ""               # height (default, no sweep)
  echo "y"               # all_sites_choice -> Y
  echo ""               # force_spin (default Y)
  echo "y"               # ml_rank_choice -> Y
  echo ""               # top_k (blank = keep all)
  echo "2x2"             # orientation grid (polar x azimuthal), merged prompt
  echo "1"               # orient_top_k -> keep only 1
  echo ""               # show_advanced (orient_rmsd_tol default applies unasked)
  echo ""               # save_report
  echo ""               # view_choice
  echo ""               # view_plots_choice
  echo ""               # press enter
  echo "0"               # quit
} | stb-suite > log_menu_orient.txt 2>&1
check_contains "Orientation sampling" log_menu_orient.txt
check_contains "2x2, top 1" log_menu_orient.txt
check_contains "4 orientation(s) sampled -> 1 unique kept" log_menu_orient.txt
check_success $RUN/sites/site_1_ontop_orient1/structure.fdf
check_contains "DFTD3                   .true." $RUN/sites/site_1_ontop_orient1/config_extra.fdf

echo "Testing: navigate 4.8.1 -> single site + orientation sampling (H2O) -> quit"
rm -rf $RUN
{
  echo "4.8.1"
  echo "structure.fdf"
  echo "calc.fdf"
  echo ""               # pp_path
  echo "H2O"             # adsorbate(s)
  echo ""               # ml_prerelax
  echo "1"               # site_choice (ontop)
  echo ""               # height (default, no sweep)
  echo "n"               # all_sites_choice -> single site
  echo "0"               # site_index
  echo ""               # force_spin (default Y)
  echo "2x2"             # orientation grid -> also turns --ml-rank on here
  echo "1"               # orient_top_k -> keep only 1
  echo ""               # show_advanced
  echo ""               # save_report
  echo ""               # view_choice
  echo ""               # view_plots_choice
  echo ""               # press enter
  echo "0"               # quit
} | stb-suite > log_menu_orient_singlesite.txt 2>&1
check_contains "ML pre-screen" log_menu_orient_singlesite.txt
check_contains "ON (single site)" log_menu_orient_singlesite.txt
check_contains "Orientation sampling" log_menu_orient_singlesite.txt
check_contains "2x2, top 1" log_menu_orient_singlesite.txt
check_success $RUN/sites/site_1_ontop_orient1/structure.fdf
check_success $RUN/sites/orientation_trajectory.xyz


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
