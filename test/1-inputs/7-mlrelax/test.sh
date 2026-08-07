#!/bin/bash

# --- Setup ---
# Smoke test for stb-mlrelax (ML Pre-Relaxation, item 1.7)
#
# Needs the optional 'ml' extra (pip install stb_suite[ml] -- PyTorch +
# mace-torch). The whole file is skipped with a clear message if `mace`
# isn't importable, so a plain (non-ml) install of the suite doesn't fail
# this test -- same gating pattern as stb-fetch's PMG_MAPI_KEY-gated
# Materials Project tests. The MACE-MP-0 "small" model (~31 MB) is
# downloaded on first use and cached under ~/.cache/mace/.
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


echo "--- Starting tester for STB-Mlrelax (item 1.7) ---"

if ! python3 -c "import mace" 2>/dev/null; then
    echo -e "${YELLOW}Skipped entirely:${NC} the optional 'ml' extra is not installed."
    echo "Install with: pip install stb_suite[ml]  (then re-run this test)"
    exit 0
fi


# --- 1. Preparation ---
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/si_wrong_a.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/si_ge_defect.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. --relax-cell: bulk Si with a wrong lattice constant ---
echo -e "\n--- Testing --relax-cell (si_wrong_a.fdf, 5.00 -> ~5.46 Ang) ---"
rm -f si_relaxed.fdf
stb-mlrelax -f si_wrong_a.fdf --relax-cell -o si_relaxed.fdf --no-intro > log_cell.txt 2>&1
check_contains "Lattice a,b,c (Ang)" log_cell.txt
check_success si_relaxed.fdf
python3 - <<'EOF'
from stb.core import structure_io
s = structure_io.to_pymatgen(structure_io.read_fdf("si_relaxed.fdf"))
a = s.lattice.abc[0]
assert 5.40 <= a <= 5.50, f"lattice constant {a} out of expected [5.40, 5.50] range"
print(f"OK: relaxed lattice constant a={a:.4f} Ang")
EOF
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} relaxed lattice constant lands in the expected [5.40, 5.50] Ang range"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} relaxed lattice constant out of range"
    FAIL=$((FAIL+1))
fi


# --- 3. Default (positions-only): Ge:Si substitutional defect ---
echo -e "\n--- Testing default positions-only mode (si_ge_defect.fdf) ---"
rm -f ge_relaxed.fdf
stb-mlrelax -f si_ge_defect.fdf -o ge_relaxed.fdf --no-intro > log_positions.txt 2>&1
check_contains "Mode.*positions only" log_positions.txt
check_contains "Displacement.*mean" log_positions.txt
check_success ge_relaxed.fdf
python3 - <<'EOF'
from stb.core import structure_io
s = structure_io.to_pymatgen(structure_io.read_fdf("ge_relaxed.fdf"))
ge = next(site for site in s if site.specie.symbol == "Ge")
si_sites = [site for site in s if site.specie.symbol == "Si"]
nearest = min(ge.distance(si) for si in si_sites)
assert 2.38 <= nearest <= 2.46, f"Ge-Si distance {nearest} out of expected [2.38, 2.46] range"
print(f"OK: relaxed Ge-Si distance = {nearest:.4f} Ang")
EOF
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} relaxed Ge-Si distance lands in the expected [2.38, 2.46] Ang range"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} relaxed Ge-Si distance out of range"
    FAIL=$((FAIL+1))
fi


# --- 4. --relax-cell on a vacuum-padded slab: relax in-plane, keep vacuum fixed ---
echo -e "\n--- Testing --relax-cell on a 2D slab (in-plane relaxes, vacuum axis held exactly fixed) ---"
stb-slab -f si_ge_defect.fdf --hkl 1 0 0 -o slab_for_test.fdf --no-intro > /dev/null 2>&1
stb-mlrelax -f slab_for_test.fdf --relax-cell -o slab_relaxed.fdf --no-intro > log_2d_cell.txt 2>&1
check_exit_code $? 0
check_contains "Vacuum axis/axes held fixed:.*c" log_2d_cell.txt
python3 - <<'EOF'
from stb.core import structure_io
before = structure_io.to_pymatgen(structure_io.read_fdf("slab_for_test.fdf"))
after = structure_io.to_pymatgen(structure_io.read_fdf("slab_relaxed.fdf"))
c_before, c_after = before.lattice.abc[2], after.lattice.abc[2]
assert abs(c_before - c_after) < 1e-6, f"vacuum axis c changed: {c_before} -> {c_after}"
print(f"OK: vacuum axis c held exactly fixed ({c_before:.4f} Ang)")
EOF
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} vacuum axis (c) unchanged by --relax-cell on a slab"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} vacuum axis (c) was altered by --relax-cell"
    FAIL=$((FAIL+1))
fi
rm -f slab_for_test.fdf slab_relaxed.fdf


# --- 4b. --relax-cell on a fully isolated molecule: no-op, falls back gracefully ---
echo -e "\n--- Testing --relax-cell on a molecule (0D, no periodic direction) ---"
stb-molecule --name H2O -o h2o_for_test.fdf --no-intro > /dev/null 2>&1
stb-mlrelax -f h2o_for_test.fdf --relax-cell -o h2o_relaxed.fdf --no-intro > log_0d.txt 2>&1
check_exit_code $? 0
check_contains "no periodic direction to relax" log_0d.txt
check_success h2o_relaxed.fdf
rm -f h2o_for_test.fdf h2o_relaxed.fdf


# --- 4c. --custom-model: relax with a fine-tuned MACE model instead of the foundation one ---
echo -e "\n--- Testing --custom-model (a real model quickly fine-tuned from stb-mlffAnalysis fixtures) ---"
MLFF_FIXTURES="$FIXTURE_DIR/../../4-workflow/17-mlff/analysis"
if command -v mace_run_train &> /dev/null && [ -d "$MLFF_FIXTURES" ]; then
    mkdir -p custom_model_src
    cp -r "$MLFF_FIXTURES"/mlff_config_* custom_model_src/
    (cd custom_model_src && stb-mlffAnalysis --path "mlff_config_*" --epochs 2 --batch-size 4 \
        --device cpu --name quick_model --skip-foundation-comparison --no-intro > ../log_quicktrain.txt 2>&1)
    QUICK_MODEL=$(find custom_model_src -maxdepth 1 -name "quick_model*.model" | sort | tail -1)
    if [ -n "$QUICK_MODEL" ]; then
        # Relax the SAME O2 reference the model was fine-tuned on (not
        # si_ge_defect.fdf) -- a model trained only on O never learned Si/Ge,
        # so it would (correctly) fail on elements outside its training set.
        rm -f custom_relaxed.fdf
        stb-mlrelax -f custom_model_src/mlff_config_001/calc.fdf --custom-model "$QUICK_MODEL" \
            -o custom_relaxed.fdf --no-intro > log_custommodel.txt 2>&1
        check_exit_code $? 0
        check_contains "Model.*custom" log_custommodel.txt
        check_success custom_relaxed.fdf
    else
        echo -e "   -> ${RED}Failed:${NC} quick fine-tuned model was not produced"
        FAIL=$((FAIL+1))
    fi
    rm -rf custom_model_src custom_relaxed.fdf
else
    echo -e "   -> ${YELLOW}SKIPPED${NC}: 'mace_run_train' not on PATH or mlff fixtures not found."
fi

echo "Testing: --custom-model with a nonexistent file (expects failure)"
stb-mlrelax -f si_ge_defect.fdf --custom-model does_not_exist.model --no-intro > log_custommodel_missing.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_custommodel_missing.txt


# --- 5. --help / --version ---
echo -e "\n--- Testing --help / --version ---"
stb-mlrelax --version > log_version.txt 2>&1
check_contains "stb-mlrelax" log_version.txt

stb-mlrelax --help > log_help.txt 2>&1
check_contains "relax-cell" log_help.txt
check_contains "fmax" log_help.txt
check_contains "model" log_help.txt
check_contains "custom-model" log_help.txt
check_contains "save-report" log_help.txt
check_contains "save-data" log_help.txt
check_contains "view" log_help.txt


# --- 5b. --save-report / --save-data / references.bib / structure validation ---
echo -e "\n--- Testing --save-report, --save-data, references.bib, and structure validation sections ---"
rm -rf report_test && mkdir report_test
(cd report_test && stb-mlrelax -f ../si_ge_defect.fdf --save-report --save-data --no-intro > console.log 2>&1)
check_success report_test/stb_mlrelax_report.txt
check_contains "STRUCTURE VALIDATION (pre-relaxation)" report_test/stb_mlrelax_report.txt
check_contains "STRUCTURE VALIDATION (post-relaxation)" report_test/stb_mlrelax_report.txt
check_contains "Atom proximity" report_test/stb_mlrelax_report.txt
check_contains "Lattice handedness" report_test/stb_mlrelax_report.txt
check_contains "\[OK\]" report_test/stb_mlrelax_report.txt
check_contains "BEFORE / AFTER COMPARISON" report_test/stb_mlrelax_report.txt
check_contains "Quantity" report_test/stb_mlrelax_report.txt
check_contains "Energy (eV)" report_test/stb_mlrelax_report.txt
check_success report_test/references.bib
check_contains "Soler2002" report_test/references.bib
check_contains "Batatia" report_test/references.bib
check_success report_test/relaxed_relax_convergence.dat
rm -rf report_test


# --- 6. Interactive path (stb-suite, shortcut 1.7) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 1.7) ---"

echo "Testing: navigate 1.7 -> invalid file then valid -> no cell relax -> defaults -> default output -> no report -> no view -> quit"
rm -f relaxed.fdf
printf '1.7\ndoes_not_exist.fdf\nsi_ge_defect.fdf\nn\n\n\n\n\n\n\n\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_contains "File not found" log_menu.txt
check_contains "Displacement" log_menu.txt
check_contains "Optimizer      : FIRE" log_menu.txt
check_success relaxed.fdf

echo "Testing: navigate 1.7 -> choosing BFGS as the optimizer"
rm -f relaxed.fdf
printf '1.7\nsi_ge_defect.fdf\nn\n\n\n\n\nBFGS\n\nn\nn\n\n0\n' | stb-suite > log_menu_bfgs.txt 2>&1
check_contains "Optimizer      : BFGS" log_menu_bfgs.txt
check_success relaxed.fdf


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
