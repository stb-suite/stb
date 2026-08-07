#!/bin/bash

# --- Setup ---
# Smoke test for stb-mlphonons (ML Phonons, item 5.2)
#
# Needs the optional 'ml' extra (pip install stb_suite[ml] -- PyTorch +
# mace-torch). The whole file is skipped with a clear message if `mace`
# isn't importable, same gating pattern as test/1-inputs/7-mlrelax and
# test/5-mlsimulations/1-mlmd. si_bulk.fdf is the same 64-atom Ge:Si
# substitutional-defect fixture used there (bulk, 3D periodic). --dim 1 1 1
# and a small --mesh/--band-points keep this a fast smoke test, not a
# converged production phonon calculation.
#
# si8.fdf + si_mlff_fixtures/ (8 real finished SIESTA single-point
# calculations, 8-atom bulk Si, real siesta binary) are used to quickly
# fine-tune a Si-compatible model on the fly (stb-mlffAnalysis), to exercise
# --custom-model / the foundation-model comparison -- si_bulk.fdf (Ge:Si)
# can't be used for that since no fixture model there is trained on Ge.
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


echo "--- Starting tester for stb-mlphonons (item 5.2) ---"

if ! python3 -c "import mace" 2>/dev/null; then
    echo -e "${YELLOW}Skipped entirely:${NC} the optional 'ml' extra is not installed."
    echo "Install with: pip install stb_suite[ml]  (then re-run this test)"
    exit 0
fi

rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/si_bulk.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 1. Basic run ---
echo -e "\n--- Testing a basic run (--dim 1 1 1, small mesh/band-points) ---"
stb-mlphonons --file si_bulk.fdf --dim 1 1 1 --mesh 8 8 8 --band-points 11 \
    --save-data --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "Displacements needed" log_basic.txt
check_contains "Detected symmetry" log_basic.txt
check_contains "Reduction" log_basic.txt
check_contains "BAND STRUCTURE" log_basic.txt
check_contains "DENSITY OF STATES" log_basic.txt
check_contains "THERMAL PROPERTIES" log_basic.txt
check_contains "Species (PDOS)    : Ge, Si" log_basic.txt
check_success mlphonons_out/phonopy_disp.yaml
check_success mlphonons_out/bands.png
check_success mlphonons_out/bands.dat
check_success mlphonons_out/dos.png
check_success mlphonons_out/dos.dat
check_success mlphonons_out/thermal.png
check_success mlphonons_out/thermal.dat

echo "Testing: dos.dat has per-species PDOS columns (total + Ge + Si)"
python3 -c "
import numpy as np, sys
data = np.loadtxt('mlphonons_out/dos.dat')
sys.exit(0 if data.shape[1] == 4 and np.all(np.isfinite(data)) else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} dos.dat has 4 columns (freq, total, Ge, Si)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} dos.dat missing expected PDOS columns"
    FAIL=$((FAIL+1))
fi

echo "Testing: acoustic sum rule correction -- no spurious imaginary-mode warning"
grep -q "No imaginary modes found" log_basic.txt
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} band path reports no imaginary modes (ASR-corrected)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} expected 'No imaginary modes found' after ASR correction"
    FAIL=$((FAIL+1))
fi

echo "Testing: thermal properties data is finite (no NaN/Inf) and entropy starts near 0 (3rd law)"
python3 -c "
import numpy as np, sys
data = np.loadtxt('mlphonons_out/thermal.dat')
temps, free_e, entropy, cv = data[:,0], data[:,1], data[:,2], data[:,3]
finite = np.all(np.isfinite(data))
sys.exit(0 if finite and entropy[0] < 1.0 else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} thermal properties finite, entropy(T~0) near 0"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} thermal properties data missing/invalid"
    FAIL=$((FAIL+1))
fi

rm -rf mlphonons_out


# --- 1b. Fine-tune a quick Si-compatible model, for --custom-model tests ---
echo -e "\n--- Fine-tuning a quick Si model to test --custom-model with (--epochs 3) ---"
cp "$FIXTURE_DIR/si8.fdf" .
cp -r "$FIXTURE_DIR/si_mlff_fixtures" mlff_config_src
mv mlff_config_src/mlff_config_* .
rmdir mlff_config_src
stb-mlffAnalysis --path "mlff_config_*" --epochs 3 --batch-size 4 --device cpu \
    --name si_quick_model --skip-foundation-comparison --no-intro > log_quicktrain.txt 2>&1
QUICK_MODEL=$(find . -maxdepth 1 -name "si_quick_model*.model" | sort | tail -1)
if [ -z "$QUICK_MODEL" ]; then
    echo -e "${RED}FAIL${NC}: could not produce a quick fine-tuned Si model to test against."
else
    echo "Using model: $QUICK_MODEL"
fi
rm -rf mlff_config_001 mlff_config_002 mlff_config_003 mlff_config_004 \
       mlff_config_005 mlff_config_006 mlff_config_007 mlff_config_008


# --- 1c. --custom-model + foundation-model comparison (default on) ---
if [ -n "$QUICK_MODEL" ]; then
    echo -e "\n--- Testing --custom-model with foundation comparison (default on) ---"
    stb-mlphonons --file si8.fdf --custom-model "$QUICK_MODEL" --dim 1 1 1 --mesh 8 8 8 \
        --band-points 11 --no-intro -o mlphonons_compare > log_compare.txt 2>&1
    check_exit_code $? 0
    check_contains "Comparison        : also evaluating with foundation model" log_compare.txt
    check_contains "\[Fine-tuned\] Minimum band-path frequency" log_compare.txt
    check_contains "\[Foundation (small)\] Minimum band-path frequency" log_compare.txt
    check_success mlphonons_compare/bands.png
    check_success mlphonons_compare/dos.png
    check_success mlphonons_compare/thermal.png
    rm -rf mlphonons_compare

    echo "Testing: --skip-foundation-comparison disables the second run"
    stb-mlphonons --file si8.fdf --custom-model "$QUICK_MODEL" --dim 1 1 1 --mesh 8 8 8 \
        --band-points 11 --skip-foundation-comparison --no-intro -o mlphonons_nocompare > log_nocompare.txt 2>&1
    check_exit_code $? 0
    grep -q "also evaluating with foundation model" log_nocompare.txt
    if [ $? -ne 0 ]; then
        echo -e "   -> ${GREEN}Verified:${NC} foundation comparison correctly skipped"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} foundation comparison ran despite --skip-foundation-comparison"
        FAIL=$((FAIL+1))
    fi
    rm -rf mlphonons_nocompare
else
    echo -e "${YELLOW}SKIPPED${NC}: comparison tests (no quick model available)."
fi


# --- 1d. --qha ---
echo -e "\n--- Testing --qha (small volume scan) ---"
stb-mlphonons --file si8.fdf --dim 1 1 1 --mesh 8 8 8 --qha --n-volumes 4 --strain-range 3.0 \
    --save-data --no-intro -o mlphonons_qha > log_qha.txt 2>&1
check_exit_code $? 0
check_contains "VOLUME SCAN" log_qha.txt
check_contains "EQUILIBRIUM FIT" log_qha.txt
check_contains "THERMAL EXPANSION SUMMARY" log_qha.txt
check_contains "B0 (bulk modulus)" log_qha.txt
check_success mlphonons_qha/qha.png
check_success mlphonons_qha/qha_thermal_expansion.dat

echo "Testing: --qha requires --n-volumes >= 4 (expects failure)"
stb-mlphonons --file si8.fdf --qha --n-volumes 3 --no-intro > log_qha_toofew.txt 2>&1
check_exit_code $? 2
check_contains "n-volumes must be at least 4" log_qha_toofew.txt

rm -rf mlphonons_qha


# --- 1e. --check-convergence ---
echo -e "\n--- Testing --check-convergence (2 supercell sizes) ---"
stb-mlphonons --file si8.fdf --check-convergence --convergence-dims 1 1 1 2 2 2 \
    --mesh 4 4 4 --band-points 5 --save-data --no-intro -o mlphonons_conv > log_conv.txt 2>&1
check_exit_code $? 0
check_contains "SUPERCELL CONVERGENCE CHECK" log_conv.txt
check_contains "1x1x1" log_conv.txt
check_contains "2x2x2" log_conv.txt
check_success mlphonons_conv/convergence.png
check_success mlphonons_conv/convergence.dat
rm -rf mlphonons_conv

echo "Testing: --convergence-dims must be a multiple of 3 (expects failure)"
stb-mlphonons --file si8.fdf --check-convergence --convergence-dims 1 1 --no-intro > log_conv_bad.txt 2>&1
check_exit_code $? 2
check_contains "multiple of 3" log_conv_bad.txt

echo "Testing: --check-convergence + --qha are mutually exclusive (expects failure)"
stb-mlphonons --file si8.fdf --check-convergence --qha --no-intro > log_conv_qha.txt 2>&1
check_exit_code $? 2
check_contains "not supported together" log_conv_qha.txt


# --- 1f. --freeze-unstable-mode (needs a genuinely unstable structure --
#     compress si8.fdf's cell to force a real imaginary mode) ---
echo -e "\n--- Testing --freeze-unstable-mode ---"
python3 -c "
from stb.core import structure_io
from pymatgen.io.ase import AseAtomsAdaptor
s = structure_io.read_fdf('si8.fdf')
pmg = structure_io.to_pymatgen(s)
atoms = AseAtomsAdaptor.get_atoms(pmg)
atoms.set_cell(atoms.get_cell() * 0.80, scale_atoms=True)
pmg2 = AseAtomsAdaptor.get_structure(atoms)
new_s = structure_io.from_pymatgen(pmg2, species_meta=s.species_meta)
structure_io.write_fdf(new_s, 'si8_compressed.fdf')
"
stb-mlphonons --file si8_compressed.fdf --dim 1 1 1 --mesh 6 6 6 --band-points 5 \
    --no-relax --freeze-unstable-mode --freeze-amplitude 0.2 --no-intro \
    -o mlphonons_freeze > log_freeze.txt 2>&1
check_exit_code $? 0
check_contains "MODE FREEZE" log_freeze.txt
check_contains "Negative (imaginary) frequency found" log_freeze.txt
check_contains "Softest mode" log_freeze.txt
check_success mlphonons_freeze/frozen_mode.fdf

echo "Testing: --freeze-unstable-mode calibrates the requested Angstrom amplitude correctly"
python3 -c "
import re, sys
text = open('log_freeze.txt').read()
m = re.search(r'requested ([0-9.]+) Ang, achieved ([0-9.]+) Ang', text)
if not m:
    sys.exit(1)
requested, achieved = float(m.group(1)), float(m.group(2))
sys.exit(0 if abs(requested - achieved) < 1e-3 else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} achieved displacement matches requested --freeze-amplitude"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} achieved displacement does not match requested --freeze-amplitude"
    FAIL=$((FAIL+1))
fi
rm -rf mlphonons_freeze

echo "Testing: no imaginary mode -- nothing to freeze"
stb-mlphonons --file si8.fdf --dim 1 1 1 --mesh 6 6 6 --band-points 5 \
    --freeze-unstable-mode --no-intro -o mlphonons_freeze_stable > log_freeze_stable.txt 2>&1
check_exit_code $? 0
check_contains "nothing to freeze" log_freeze_stable.txt
rm -rf mlphonons_freeze_stable


# --- 1g. --animate-mode ---
echo -e "\n--- Testing --animate-mode ---"
stb-mlphonons --file si8.fdf --dim 1 1 1 --mesh 4 4 4 --band-points 5 \
    --animate-mode 5 --animate-frames 6 --animate-amplitude 0.3 --animate-format xsf \
    --no-intro -o mlphonons_animate > log_animate.txt 2>&1
check_exit_code $? 0
check_contains "MODE ANIMATION" log_animate.txt
check_contains "Band index        : 5" log_animate.txt
check_success mlphonons_animate/mode_5.xsf
rm -rf mlphonons_animate

echo "Testing: --animate-mode out of range is flagged instead of crashing"
stb-mlphonons --file si8.fdf --dim 1 1 1 --mesh 4 4 4 --band-points 5 \
    --animate-mode 99 --no-intro -o mlphonons_animate_bad > log_animate_bad.txt 2>&1
check_exit_code $? 0
check_contains "out of range" log_animate_bad.txt
rm -rf mlphonons_animate_bad


# --- 2. --custom-model with a nonexistent file (expects failure) ---
echo -e "\n--- Testing --custom-model with a nonexistent file (expects failure) ---"
stb-mlphonons --file si_bulk.fdf --custom-model does_not_exist.model --no-intro > log_custommodel.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_custommodel.txt


# --- 3. --save-report ---
echo -e "\n--- Testing --save-report ---"
stb-mlphonons --file si_bulk.fdf --dim 1 1 1 --mesh 8 8 8 --band-points 11 \
    -o report_test --save-report --no-intro > log_savereport.txt 2>&1
check_exit_code $? 0
check_success report_test/stb_mlphonons_report.txt
check_contains "STB-MLPHONONS REPORT" report_test/stb_mlphonons_report.txt
rm -rf report_test


# --- 4. Robustness / errors ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: missing --file"
stb-mlphonons --no-intro > log_missing_file.txt 2>&1
check_exit_code $? 2

echo "Testing: nonexistent input file"
stb-mlphonons --file nope.fdf --no-intro > log_missing_input.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_input.txt

echo "Testing: --version"
stb-mlphonons --version > log_version.txt 2>&1
check_contains "stb-mlphonons" log_version.txt

echo "Testing: --help documents --dim, --custom-model, --mesh, --no-relax, --qha, --skip-foundation-comparison"
stb-mlphonons --help > log_help.txt 2>&1
check_contains "dim" log_help.txt
check_contains "custom-model" log_help.txt
check_contains "mesh" log_help.txt
check_contains "no-relax" log_help.txt
check_contains "\-\-qha" log_help.txt
check_contains "skip-foundation-comparison" log_help.txt
check_contains "n-volumes" log_help.txt
check_contains "strain-range" log_help.txt
check_contains "eos" log_help.txt
check_contains "check-convergence" log_help.txt
check_contains "convergence-dims" log_help.txt
check_contains "freeze-unstable-mode" log_help.txt
check_contains "freeze-amplitude" log_help.txt
check_contains "animate-mode" log_help.txt
check_contains "animate-qpoint" log_help.txt
check_contains "animate-dim" log_help.txt
check_contains "animate-amplitude" log_help.txt
check_contains "animate-frames" log_help.txt
check_contains "animate-format" log_help.txt


# --- 5. Interactive path (stb-suite, shortcut 5.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 5.2) ---"
rm -rf interactive_out
printf '5.2\nsi_bulk.fdf\n\nsmall\n\n1 1 1\ninteractive_out\n\nn\nn\nn\nn\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success interactive_out/bands.png
check_success interactive_out/dos.png
check_success interactive_out/thermal.png


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
