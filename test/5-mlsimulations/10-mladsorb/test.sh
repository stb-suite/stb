#!/bin/bash

# --- Setup ---
# Smoke test for stb-mladsorb (ML Adsorption-Site Screening, item 5.10)
#
# Needs the optional 'ml' extra (pip install stb_suite[ml] -- PyTorch +
# mace-torch). The whole file is skipped with a clear message if `mace`
# isn't importable, same gating pattern as the other test/5-mlsimulations/*
# tools.
#
# graphene_slab.fdf: 2-atom graphene monolayer (vacuum along c), copied
# from test/4-workflow/8-adsorption/prep/structure.fdf (the real fixture
# stb-adsorb's own SIESTA workflow test already uses). Verified live: H
# chemisorption on graphene gives an ML adsorption energy of ~-0.7 to
# -1.0 eV depending on site (bridge most stable here), the right physical
# magnitude for a real H-carbon covalent bond, matching typical DFT
# literature values within the expected MACE-vs-DFT heuristic spread.
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


echo "--- Starting tester for stb-mladsorb (item 5.10) ---"

if ! python3 -c "import mace" 2>/dev/null; then
    echo -e "${YELLOW}Skipped entirely:${NC} the optional 'ml' extra is not installed."
    echo "Install with: pip install stb_suite[ml]  (then re-run this test)"
    exit 0
fi

rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/graphene_slab.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 1. Basic run: all site types, single adsorbate ---
echo -e "\n--- Testing a basic screen (H on graphene, all site types) ---"
stb-mladsorb --slab graphene_slab.fdf --adsorbate H --site-type all --height 2.0 \
    --save-data --save-report --no-intro -o h_screen > log_screen.txt 2>&1
check_exit_code $? 0
check_contains "REFERENCE ENERGIES" log_screen.txt
check_contains "SITE SCREENING" log_screen.txt
check_contains "RANKING" log_screen.txt
check_contains "Best site" log_screen.txt
check_success h_screen/adsorption_sites.png
check_success h_screen/site_ranking.png
check_success h_screen/site_ranking.dat
check_success h_screen/stb_mladsorb_report.txt
check_contains "STB-MLADSORB REPORT" h_screen/stb_mladsorb_report.txt

echo "Testing: 4 candidate sites found (2-atom graphene: 2 ontop + 2 bridge, no hollow for a 2-atom basis)"
python3 -c "
import sys
with open('h_screen/site_ranking.dat') as f:
    lines = [l for l in f if not l.startswith('#')]
sys.exit(0 if len(lines) == 4 else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} 4 candidate sites found and ranked"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} did not find/rank exactly 4 candidate sites"
    FAIL=$((FAIL+1))
fi

echo "Testing: best (most negative) adsorption energy is physically sensible chemisorption (-3.0 to -0.1 eV)"
python3 -c "
import numpy as np, sys
data = np.genfromtxt('h_screen/site_ranking.dat', comments='#', dtype=None, encoding='utf-8')
e_ads = [float(row[5]) for row in np.atleast_1d(data)]
best = min(e_ads)
sys.exit(0 if -3.0 < best < -0.1 else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} best adsorption energy is physically sensible chemisorption"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} best adsorption energy is outside the expected physical range"
    FAIL=$((FAIL+1))
fi

echo "Testing: ranking is sorted ascending (rank 1 = most negative/stable)"
python3 -c "
import numpy as np, sys
data = np.genfromtxt('h_screen/site_ranking.dat', comments='#', dtype=None, encoding='utf-8')
e_ads = [float(row[5]) for row in np.atleast_1d(data)]
sys.exit(0 if all(e_ads[i] <= e_ads[i+1] for i in range(len(e_ads)-1)) else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} sites are ranked ascending by adsorption energy"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} ranking is not sorted ascending"
    FAIL=$((FAIL+1))
fi
rm -rf h_screen


# --- 2. --site-type filter + --top-k ---
echo -e "\n--- Testing --site-type ontop --top-k 1 ---"
stb-mladsorb --slab graphene_slab.fdf --adsorbate H --site-type ontop --top-k 1 \
    --no-intro -o h_ontop > log_ontop.txt 2>&1
check_exit_code $? 0
check_success h_ontop/rank01_site1_ontop_H.fdf
check_contains "ontop" log_ontop.txt
echo "Testing: --site-type ontop excludes bridge sites entirely"
check_not_contains_local() {
    if grep -q -- "bridge" "$1" 2>/dev/null; then
        echo -e "   -> ${RED}Failed:${NC} 'bridge' unexpectedly found in '$1'"
        FAIL=$((FAIL+1))
    else
        echo -e "   -> ${GREEN}Verified:${NC} 'bridge' NOT found in '$1' (as expected)"
        PASS=$((PASS+1))
    fi
}
check_not_contains_local log_ontop.txt
rm -rf h_ontop


# --- 2b. --height-sweep produces a binding-energy-vs-height curve ---
echo -e "\n--- Testing --height-sweep (binding-energy curve) ---"
stb-mladsorb --slab graphene_slab.fdf --adsorbate H --site-type ontop \
    --height-sweep 1.0 2.5 0.5 --save-data --no-intro -o h_sweep > log_sweep.txt 2>&1
check_exit_code $? 0
check_success h_sweep/height_curve.png
python3 -c "
import sys
with open('h_sweep/site_ranking.dat') as f:
    lines = [l for l in f if not l.startswith('#')]
sys.exit(0 if len(lines) == 4 else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} 4 heights screened for the single ontop site"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} did not screen the expected 4 heights"
    FAIL=$((FAIL+1))
fi
rm -rf h_sweep


# --- 2c. --n-orientations improves (or matches) a poorly-oriented molecular adsorbate ---
echo -e "\n--- Testing --n-orientations (H2O, orientation sampling) ---"
stb-mladsorb --slab graphene_slab.fdf --adsorbate H2O --site-type ontop --height 2.5 \
    --n-orientations 1 --seed 1 --no-intro -o h2o_noorient > log_h2o_1.txt 2>&1
check_exit_code $? 0
stb-mladsorb --slab graphene_slab.fdf --adsorbate H2O --site-type ontop --height 2.5 \
    --n-orientations 8 --seed 1 --no-intro -o h2o_orient > log_h2o_8.txt 2>&1
check_exit_code $? 0
echo "Testing: sampling 8 orientations finds an energy at least as good as the default one"
python3 -c "
import re, sys
def best_e(path):
    text = open(path).read()
    m = re.search(r'Best site.*E_ads = (-?[\d.]+) eV', text)
    return float(m.group(1))
e1 = best_e('log_h2o_1.txt')
e8 = best_e('log_h2o_8.txt')
sys.exit(0 if e8 <= e1 + 1e-6 else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} orientation sampling finds an energy at least as good as the default"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} orientation sampling gave a worse result than the default orientation"
    FAIL=$((FAIL+1))
fi
rm -rf h2o_noorient h2o_orient


# --- 2d. --diffusion-barrier: NEB between the 2 best sites, cross-checked
#     against the independently-computed site-ranking energies ---
echo -e "\n--- Testing --diffusion-barrier (cross-checked against site ranking) ---"
stb-mladsorb --slab graphene_slab.fdf --adsorbate H --site-type all --height 2.0 \
    --diffusion-barrier --diffusion-n-images 5 --save-data --save-report \
    --no-intro -o h_diff > log_diff.txt 2>&1
check_exit_code $? 0
check_contains "SURFACE DIFFUSION BARRIER" log_diff.txt
check_contains "Diffusion path" log_diff.txt
check_contains "Diffusion barrier" log_diff.txt

echo "Testing: NEB reaction energy matches the independently-computed site-ranking energy difference"
python3 -c "
import re, sys
import numpy as np
text = open('log_diff.txt').read()
m_path = re.search(r'Diffusion path\s*:\s*site (\d+) .* site (\d+)', text)
m_de = re.search(r'Reaction energy\s*:\s*(-?[\d.]+) eV', text)
if not (m_path and m_de):
    sys.exit(1)
site_a, site_b = int(m_path.group(1)), int(m_path.group(2))
de_neb = float(m_de.group(1))
data = np.genfromtxt('h_diff/site_ranking.dat', comments='#', dtype=None, encoding='utf-8')
e_by_site = {int(row[1]): float(row[5]) for row in np.atleast_1d(data)}
de_ranking = e_by_site[site_b] - e_by_site[site_a]
sys.exit(0 if abs(de_neb - de_ranking) < 0.01 else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} NEB reaction energy matches the site-ranking energy difference"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} NEB reaction energy does not match the site-ranking energy difference"
    FAIL=$((FAIL+1))
fi
rm -rf h_diff


# --- 2e. Foundation-model comparison (quick-fine-tuned Si model, reusing
#     test/5-mlsimulations/2-mlphonons's real SIESTA fixtures, on a
#     synthetic single-atom Si 'slab' + Si adsorbate so every element
#     involved is one the quick model was actually trained on -- there's
#     no committed carbon+hydrogen SIESTA dataset to fine-tune against,
#     so this section deliberately doesn't reuse graphene_slab.fdf/H) ---
echo -e "\n--- Fine-tuning a quick Si model to test the foundation-model comparison with ---"
PHONONS_FIXTURES="$FIXTURE_DIR/../2-mlphonons/si_mlff_fixtures"
if [ -d "$PHONONS_FIXTURES" ]; then
    cp -r "$PHONONS_FIXTURES" mlff_config_src
    mv mlff_config_src/mlff_config_* .
    rmdir mlff_config_src
    stb-mlffAnalysis --path "mlff_config_*" --epochs 3 --batch-size 4 --device cpu \
        --name si_quick_model --skip-foundation-comparison --no-intro > log_quicktrain.txt 2>&1
    QUICK_MODEL=$(find . -maxdepth 1 -name "si_quick_model*.model" | sort | tail -1)
    rm -rf mlff_config_001 mlff_config_002 mlff_config_003 mlff_config_004 \
           mlff_config_005 mlff_config_006 mlff_config_007 mlff_config_008
else
    QUICK_MODEL=""
fi

if [ -n "$QUICK_MODEL" ]; then
    echo "Using model: $QUICK_MODEL"
    python3 -c "
from stb.core import structure_io
from ase import Atoms
from pymatgen.io.ase import AseAtomsAdaptor
atoms = Atoms('Si', positions=[[1.5, 1.5, 10.0]], cell=[3.0, 3.0, 20.0], pbc=True)
pmg = AseAtomsAdaptor.get_structure(atoms)
fdf_s = structure_io.from_pymatgen(pmg, coord_format='fractional')
structure_io.write_fdf(fdf_s, 'si_slab.fdf')
"
    echo -e "\n--- Testing --custom-model with foundation comparison (default on) ---"
    stb-mladsorb --slab si_slab.fdf --adsorbate Si --site-type ontop --height 2.0 \
        --custom-model "$QUICK_MODEL" --no-intro -o si_compare > log_compare.txt 2>&1
    check_exit_code $? 0
    check_contains "Comparison        : also evaluating with foundation model" log_compare.txt
    check_contains "FOUNDATION MODEL COMPARISON" log_compare.txt
    check_contains "Best site (fine-tuned)" log_compare.txt
    check_contains "Best site (foundation)" log_compare.txt
    check_success si_compare/site_ranking.png
    rm -rf si_compare

    echo "Testing: --skip-foundation-comparison disables the second run"
    stb-mladsorb --slab si_slab.fdf --adsorbate Si --site-type ontop --height 2.0 \
        --custom-model "$QUICK_MODEL" --skip-foundation-comparison --no-intro \
        -o si_nocompare > log_nocompare.txt 2>&1
    check_exit_code $? 0
    if grep -q "FOUNDATION MODEL COMPARISON" log_nocompare.txt; then
        echo -e "   -> ${RED}Failed:${NC} foundation comparison ran despite --skip-foundation-comparison"
        FAIL=$((FAIL+1))
    else
        echo -e "   -> ${GREEN}Verified:${NC} foundation comparison correctly skipped"
        PASS=$((PASS+1))
    fi
    rm -rf si_nocompare si_slab.fdf
else
    echo -e "${YELLOW}SKIPPED${NC}: foundation-comparison tests (no quick model available)."
fi


# --- 3. Error and robustness cases ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: unknown adsorbate"
stb-mladsorb --slab graphene_slab.fdf --adsorbate Xx --no-intro > log_badads.txt 2>&1
check_exit_code $? 1
check_contains "not a recognized element symbol" log_badads.txt

echo "Testing: missing --slab"
stb-mladsorb --adsorbate H --no-intro > log_missing_file.txt 2>&1
check_exit_code $? 2

echo "Testing: nonexistent slab file"
stb-mladsorb --slab nope.fdf --adsorbate H --no-intro > log_missing_slab.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_slab.txt

echo "Testing: --custom-model with a nonexistent file"
stb-mladsorb --slab graphene_slab.fdf --adsorbate H --custom-model does_not_exist.model \
    --no-intro > log_custommodel.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_custommodel.txt

echo "Testing: --version"
stb-mladsorb --version > log_version.txt 2>&1
check_contains "stb-mladsorb" log_version.txt

echo "Testing: --help documents --slab, --adsorbate, --site-type, --custom-model, --top-k"
stb-mladsorb --help > log_help.txt 2>&1
check_contains "slab" log_help.txt
check_contains "adsorbate" log_help.txt
check_contains "site-type" log_help.txt
check_contains "custom-model" log_help.txt
check_contains "top-k" log_help.txt
check_contains "height-sweep" log_help.txt
check_contains "skip-foundation-comparison" log_help.txt
check_contains "n-orientations" log_help.txt
check_contains "diffusion-barrier" log_help.txt
check_contains "diffusion-adsorbate" log_help.txt
check_contains "diffusion-n-images" log_help.txt


# --- 4. Interactive path (stb-suite, shortcut 5.10) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 5.10) ---"
rm -rf interactive_out
printf '5.10\ngraphene_slab.fdf\nH\n\nsmall\nontop\n2.0\n1\n\nn\ninteractive_out\nn\nn\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success interactive_out/site_ranking.png


popd > /dev/null

# --- 5. Summary ---
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
