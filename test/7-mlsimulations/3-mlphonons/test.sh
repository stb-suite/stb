#!/bin/bash

# --- Setup ---
# Smoke test for stb-mlphonons (ML Phonons, item 7.3)
#
# Needs the optional 'ml' extra (pip install stb_suite[ml] -- PyTorch +
# mace-torch). The whole file is skipped with a clear message if `mace`
# isn't importable, same gating pattern as test/1-inputs/7-mlrelax and
# test/7-mlsimulations/2-mlmd. si_bulk.fdf is the same 64-atom Ge:Si
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

check_not_contains() {
    if grep -q -- "$1" "$2" 2>/dev/null; then
        echo -e "   -> ${RED}Failed:${NC} '$1' found in '$2'"
        FAIL=$((FAIL+1))
    else
        echo -e "   -> ${GREEN}Verified:${NC} '$1' not in '$2'"
        PASS=$((PASS+1))
    fi
}

check_python() {   # $1 = description, $2 = python code exiting 0 on success
    python3 -c "$2"
    if [ $? -eq 0 ]; then
        echo -e "   -> ${GREEN}Verified:${NC} $1"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} $1"
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


echo "--- Starting tester for stb-mlphonons (item 7.3) ---"

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
check_contains "Displacements     :" log_basic.txt
check_contains "Residual forces   : max" log_basic.txt
check_contains "subtracted from every displacement" log_basic.txt
check_contains "Symmetry          :" log_basic.txt
check_contains "fewer evaluations" log_basic.txt
check_contains "Acoustic sum rule :" log_basic.txt
for sec in "\[0\] RUN METADATA" "\[1\] STRUCTURE" "\[2\] MODEL" "\[3\] DISPLACEMENTS & FORCE CONSTANTS" \
           "\[4\] FREQUENCIES AT HIGH-SYMMETRY POINTS" "\[5\] BAND STRUCTURE & STABILITY" "\[6\] DENSITY OF STATES" \
           "\[7\] THERMAL PROPERTIES" "\[8\] SUMMARY & FILES" "\[9\] WARNINGS"; do
    check_contains "$sec" log_basic.txt
done
check_contains "Species (PDOS)    : Ge, Si" log_basic.txt
check_contains "Sound velocities" log_basic.txt
check_contains "Debye temperature" log_basic.txt
check_contains "Zero-point energy :" log_basic.txt
check_contains "Gamma-point modes" log_basic.txt
check_contains "per 64-atom primitive cell" log_basic.txt
check_contains "more) ..." log_basic.txt
echo "Testing: library messages (MACE/torch/phonopy) are counted, not printed in the report"
check_not_contains "UserWarning" log_basic.txt
check_not_contains "PrimitiveMatrixAutoDefaultWarning" log_basic.txt
check_contains "Library messages  :" log_basic.txt
echo "Testing: plots are written, gnuplot/--show are opt-in"
check_python "no plot/ folder without --save-gnuplot" "
import os, sys
sys.exit(0 if not os.path.exists('mlphonons_out/plot') else 1)
"
check_contains "Also available    : --save-gnuplot" log_basic.txt
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
stb-mlffAnalysis --path "mlff_config_*" --committee 1 -o . --epochs 3 --batch-size 4 --device cpu \
    --name si_quick_model --skip-foundation-comparison --no-intro > log_quicktrain.txt 2>&1
QUICK_MODEL=$(find . -maxdepth 1 -name "si_quick_model*.model" | sort | tail -1)
QUICK_CARD="$(dirname "$QUICK_MODEL")/mlff_model_card.json"
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
    stb-mlphonons --file si8.fdf --custom-model "$QUICK_MODEL" --custom-model-card "$QUICK_CARD" --dim 1 1 1 --mesh 8 8 8 \
        --band-points 11 --no-intro -o mlphonons_compare > log_compare.txt 2>&1
    check_exit_code $? 0
    check_contains "Comparison        : also evaluating with foundation model" log_compare.txt
    check_contains "\[Fine-tuned\] Band path" log_compare.txt
    check_contains "\[Foundation (medium)\] Band path" log_compare.txt
    check_contains "FINE-TUNED vs FOUNDATION" log_compare.txt
    check_contains "Model card -- identity" log_compare.txt
    check_success mlphonons_compare/bands.png
    check_success mlphonons_compare/dos.png
    check_success mlphonons_compare/thermal.png
    rm -rf mlphonons_compare

    echo "Testing: --skip-foundation-comparison disables the second run"
    stb-mlphonons --file si8.fdf --custom-model "$QUICK_MODEL" --custom-model-card "$QUICK_CARD" --dim 1 1 1 --mesh 8 8 8 \
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
check_contains "MODE GRUENEISEN PARAMETERS" log_qha.txt
check_contains "Cell factor       : 4 (input cell 8 atoms / primitive 2)" log_qha.txt
check_success mlphonons_qha/qha.png
check_success mlphonons_qha/qha_thermal_expansion.dat
echo "Testing: the static B0 is printed in GPa (it was eV/Ang^3 labelled GPa: 0.58 for silicon)"
check_python "B0 of silicon between 40 and 150 GPa" "
import re, sys
t = re.sub(r'\x1b\[[0-9;]*m', '', open('log_qha.txt').read())
m = re.search(r'B0 \(bulk modulus\) : ([\d.]+) GPa', t)
sys.exit(0 if m and 40 < float(m.group(1)) < 150 else 1)
"
echo "Testing: the QHA free energy is per INPUT cell (8 atoms), not per 2-atom primitive cell"
check_python "F_vib in vol_*/phonopy_disp.yaml x cell factor reproduces the report's F_vib(300K) per atom" "
import glob, re, sys, warnings, numpy as np, phonopy
warnings.filterwarnings('ignore')
p = phonopy.load(sorted(glob.glob('mlphonons_qha/vol_*'))[0] + '/phonopy_disp.yaml', log_level=0)
p.run_mesh([8, 8, 8]); p.run_thermal_properties(temperatures=[300])
f_prim = p.get_thermal_properties_dict()['free_energy'][0]          # kJ/mol per 2-atom cell
expect = f_prim * 4 / 8 * 10.3642688                                 # meV/atom
t = re.sub(r'\x1b\[[0-9;]*m', '', open('log_qha.txt').read())
row = [l for l in t.splitlines() if re.match(r'\s+-3\.00\s', l)][0]
got = float(row.split()[8])
sys.exit(0 if abs(got - expect) < 0.05 else 1)
"

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
check_contains "RMSE THz" log_conv.txt
check_contains "Frequencies at the special points" log_conv.txt
check_success mlphonons_conv/convergence.png
check_success mlphonons_conv/convergence.dat
echo "Testing: convergence is judged on the band path, not on the Gamma maximum (identical for every supercell)"
check_python "1x1x1 differs from the 2x2x2 reference along the band path (RMSE > 0.01 THz), the reference itself has RMSE 0" "
import sys, numpy as np
d = np.genfromtxt('mlphonons_conv/convergence.dat', dtype=None, encoding=None)
sys.exit(0 if d[0][6] > 0.01 and d[1][6] == 0 else 1)
"
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

echo "Testing: --animate-mode with no value animates the first optical mode (band 3)"
stb-mlphonons --file si8.fdf --dim 1 1 1 --mesh 4 4 4 --band-points 5 \
    --animate-mode --animate-frames 6 --no-intro -o mlphonons_animate_opt > log_animate_opt.txt 2>&1
check_exit_code $? 0
check_contains "first optical mode" log_animate_opt.txt
check_contains "Band index        : 3" log_animate_opt.txt
check_success mlphonons_animate_opt/mode_3.xsf
rm -rf mlphonons_animate_opt

echo "Testing: --animate-mode out of range is flagged instead of crashing"
stb-mlphonons --file si8.fdf --dim 1 1 1 --mesh 4 4 4 --band-points 5 \
    --animate-mode 99 --no-intro -o mlphonons_animate_bad > log_animate_bad.txt 2>&1
check_exit_code $? 0
check_contains "out of range" log_animate_bad.txt
rm -rf mlphonons_animate_bad


# --- 1h. --save-gnuplot writes plot/ and runs gnuplot; --show without a display is a warning ---
echo -e "\n--- Testing --save-gnuplot and --show ---"
MPLBACKEND=Agg stb-mlphonons --file si8.fdf --dim 1 1 1 --mesh 6 6 6 --band-points 11 --save-gnuplot --show \
    --no-intro -o mlphonons_gp > log_gp.txt 2>&1
check_exit_code $? 0
for f in bands.gplot bands_0.dat dos.gplot dos_0.dat thermal.gplot thermal_0.dat; do
    check_success mlphonons_gp/plot/$f
done
if command -v gnuplot &> /dev/null; then
    for f in bands.pdf dos.pdf thermal.pdf; do check_success mlphonons_gp/plot/$f; done
    check_contains "gnuplot figure" log_gp.txt
fi
check_contains "show needs a display" log_gp.txt
MPLBACKEND=Agg stb-mlphonons --file si8.fdf --qha --n-volumes 4 --dim 1 1 1 --mesh 6 6 6 --save-gnuplot \
    --no-intro -o mlphonons_gp_qha > log_gp_qha.txt 2>&1
check_exit_code $? 0
check_success mlphonons_gp_qha/plot/qha.gplot
check_success mlphonons_gp_qha/plot/qha_fv.dat
if command -v gnuplot &> /dev/null; then check_success mlphonons_gp_qha/plot/qha.pdf; fi
rm -rf mlphonons_gp mlphonons_gp_qha


# --- 2. --custom-model with a nonexistent file (expects failure) ---
echo -e "\n--- Testing --custom-model with a nonexistent file (expects failure) ---"
stb-mlphonons --file si_bulk.fdf --custom-model does_not_exist.model --custom-model-card does_not_exist.json --no-intro > log_custommodel.txt 2>&1
check_exit_code $? 2
check_contains "not found" log_custommodel.txt

echo "Testing: a custom model without its model card (--custom-model-card) is rejected"
stb-mlphonons --file si_bulk.fdf --custom-model does_not_exist.model --no-intro > log_custommodel_nocard.txt 2>&1
check_exit_code $? 2
check_contains "custom-model-card" log_custommodel_nocard.txt


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
check_contains "custom-model-card" log_help.txt
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
check_contains "save-gnuplot" log_help.txt
check_contains "\-\-show" log_help.txt
check_contains "convergence-tol" log_help.txt


# --- 5. Interactive path (stb-suite, shortcut 7.3) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 7.3) ---"
rm -rf interactive_out
printf '7.3\nsi_bulk.fdf\n\nsmall\n\n1 1 1\ninteractive_out\n\nn\nn\nn\ny\nn\nn\n\n0\n' | MPLBACKEND=Agg stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success interactive_out/bands.png
check_success interactive_out/dos.png
check_success interactive_out/thermal.png
check_success interactive_out/plot/bands.gplot

# the device question must also come up with a custom model (it used to be skipped: always cpu)
if [ -n "$QUICK_MODEL" ]; then
    echo "Testing: the 7.3 menu asks for the device with a custom model too"
    rm -rf interactive_custom_out
    printf '7.3\nsi8.fdf\n%s\n%s\ncpu\n1 1 1\ninteractive_custom_out\n\nn\nn\nn\nn\nn\nn\nn\n\n0\n' \
        "$QUICK_MODEL" "$QUICK_CARD" | MPLBACKEND=Agg stb-suite > log_menu_custom.txt 2>&1
    check_exit_code $? 0
    check_contains "Device \[cpu/cuda" log_menu_custom.txt
    check_contains "Device            : cpu" log_menu_custom.txt
    check_success interactive_custom_out/bands.png
fi


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
