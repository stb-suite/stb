#!/bin/bash

# --- Setup ---
# Smoke test for stb-mlmd (ML Molecular Dynamics, item 5.1)
#
# Needs the optional 'ml' extra (pip install stb_suite[ml]). The whole file
# is skipped with a clear message if `mace` isn't importable, same gating
# pattern as test/1-inputs/6-mlrelax/test.sh. si_bulk.fdf is the same
# 64-atom Ge:Si substitutional-defect fixture used there (bulk, 3D
# periodic -- suitable for NVE/NVT/NPT alike).
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


echo "--- Starting tester for stb-mlmd (item 5.1) ---"

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


# --- 1. NVT (default ensemble) ---
echo -e "\n--- Testing NVT (default), 100 steps, stride 20 ---"
stb-mlmd --file si_bulk.fdf --ensemble nvt --temperature 300 --n-steps 100 --stride 20 \
    -o nvt_traj.xsf --output-structure nvt_final.fdf --save-data --no-intro > log_nvt.txt 2>&1
check_exit_code $? 0
check_contains "Ensemble          : NVT" log_nvt.txt
check_contains "Collected 5 frame(s)" log_nvt.txt
check_contains "Temperature       : mean" log_nvt.txt
check_success nvt_traj.xsf
check_success nvt_final.fdf
check_success si_bulk_md_diagnostics.png
check_success si_bulk_md_diagnostics.dat


# --- 1b. --equilibration-steps discards the initial transient ---
echo -e "\n--- Testing --equilibration-steps ---"
stb-mlmd --file si_bulk.fdf --ensemble nvt --equilibration-steps 20 --n-steps 40 --stride 10 \
    --no-intro > log_equil.txt 2>&1
check_exit_code $? 0
check_contains "EQUILIBRATION (discarded)" log_equil.txt
check_contains "Equilibration done" log_equil.txt
check_contains "20 equilibration" log_equil.txt


# --- 2. NVE (checks total-energy conservation reporting) ---
echo -e "\n--- Testing NVE (energy conservation report) ---"
stb-mlmd --file si_bulk.fdf --ensemble nve --n-steps 200 --stride 20 --seed 1 \
    -o nve_traj.xsf --output-structure nve_final.fdf --no-intro > log_nve.txt 2>&1
check_exit_code $? 0
check_contains "Ensemble          : NVE" log_nve.txt
check_contains "Relative drift" log_nve.txt
check_success nve_traj.xsf

echo "Testing: NVE total-energy drift is small (well-conserved integrator)"
python3 -c "
import re, sys
text = open('log_nve.txt').read()
m = re.search(r'Relative drift    : ([\d.eE+-]+)%', text)
sys.exit(0 if m and abs(float(m.group(1))) < 1.0 else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} NVE relative energy drift is well under 1%"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} NVE relative energy drift missing or too large"
    FAIL=$((FAIL+1))
fi


# --- 3. NPT on bulk (allowed) ---
echo -e "\n--- Testing NPT on a bulk structure (allowed) ---"
stb-mlmd --file si_bulk.fdf --ensemble npt --temperature 300 --pressure 0.0 --n-steps 40 --stride 10 \
    -o npt_traj.xsf --output-structure npt_final.fdf --no-intro > log_npt.txt 2>&1
check_exit_code $? 0
check_contains "Ensemble          : NPT" log_npt.txt
check_success npt_traj.xsf


# --- 4. NPT on a vacuum-padded slab (must fail) ---
echo -e "\n--- Testing NPT on a 2D slab (expects failure) ---"
stb-slab -f si_bulk.fdf --hkl 1 0 0 -o slab_for_test.fdf --no-intro > /dev/null 2>&1
stb-mlmd --file slab_for_test.fdf --ensemble npt --n-steps 10 --no-intro > log_npt_slab.txt 2>&1
check_exit_code $? 1
check_contains "requires a bulk" log_npt_slab.txt
rm -f slab_for_test.fdf


# --- 5. --custom-model with a nonexistent file (expects failure) ---
echo -e "\n--- Testing --custom-model with a nonexistent file (expects failure) ---"
stb-mlmd --file si_bulk.fdf --custom-model does_not_exist.model --n-steps 10 --no-intro > log_custommodel.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_custommodel.txt


# --- 6. --out-format xyz embeds per-frame Epot/Temp/Time ---
echo -e "\n--- Testing --out-format xyz (per-frame Epot/Temp/Time embedding) ---"
stb-mlmd --file si_bulk.fdf --ensemble nvt --n-steps 40 --stride 10 --out-format xyz \
    -o nvt_traj.xyz --output-structure xyz_final.fdf --no-intro > log_xyz.txt 2>&1
check_exit_code $? 0
python3 -c "
import sys
from ase.io import read
frames = read('nvt_traj.xyz', index=':')
sys.exit(0 if len(frames) == 4 and all(
    'Epot' in a.info and 'Temp' in a.info and 'Time' in a.info for a in frames
) else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} xyz trajectory has 4 frames with Epot/Temp/Time embedded"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} xyz trajectory missing expected frames/fields"
    FAIL=$((FAIL+1))
fi


# --- 7. --save-report ---
echo -e "\n--- Testing --save-report ---"
rm -f stb_mlmd_report.txt
stb-mlmd --file si_bulk.fdf --n-steps 10 --save-report --no-intro > log_savereport.txt 2>&1
check_success stb_mlmd_report.txt
check_contains "STB-MLMD REPORT" stb_mlmd_report.txt
rm -f stb_mlmd_report.txt


# --- 8. Robustness / errors ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: missing --file"
stb-mlmd --no-intro > log_missing_file.txt 2>&1
check_exit_code $? 2

echo "Testing: nonexistent input file"
stb-mlmd --file nope.fdf --no-intro > log_missing_input.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_input.txt

echo "Testing: --version"
stb-mlmd --version > log_version.txt 2>&1
check_contains "stb-mlmd" log_version.txt

echo "Testing: --help documents --ensemble, --custom-model, --friction, --pressure, --equilibration-steps, --taut/--taup/--compressibility"
stb-mlmd --help > log_help.txt 2>&1
check_contains "ensemble" log_help.txt
check_contains "custom-model" log_help.txt
check_contains "friction" log_help.txt
check_contains "pressure" log_help.txt
check_contains "equilibration-steps" log_help.txt
check_contains "taut" log_help.txt
check_contains "taup" log_help.txt
check_contains "compressibility" log_help.txt

echo "Testing: --taut/--taup/--compressibility are accepted on an NPT run"
stb-mlmd --file si_bulk.fdf --ensemble npt --n-steps 20 --stride 10 \
    --taut 50 --taup 300 --compressibility 1.0e-6 --no-intro > log_nptparams.txt 2>&1
check_exit_code $? 0


# --- 9. Interactive path (stb-suite, shortcut 5.1) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 5.1) ---"
rm -f interactive_traj.xsf md_final.fdf
printf '5.1\nsi_bulk.fdf\n\nsmall\n1\n300\n0\n30\n1.0\nn\ninteractive_traj.xsf\nn\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success interactive_traj.xsf
check_success md_final.fdf


popd > /dev/null

# --- 10. Summary ---
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
