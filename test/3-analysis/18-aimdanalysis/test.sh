#!/bin/bash

# --- Setup ---
# Smoke test for stb-aimdAnalysis (AIMD Trajectory Analysis, item 3.18).
#
# Fixtures are a real SIESTA AIMD run (5-step Verlet, O2 dimer in vacuum,
# same aimd.ANI/.XV/.fdf/.out used by test/6-utils/5-ani2traj) -- too short
# for MSD/VDOS to carry real physical meaning, but enough to confirm every
# plot/data file is generated without error and the numbers aren't NaN.
export MPLBACKEND=Agg
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


# --- 1. Preparation ---
echo "--- Starting tester for stb-aimdAnalysis (item 3.18) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR"/aimd.ANI "$FIXTURE_DIR"/aimd.fdf "$FIXTURE_DIR"/aimd.out "$FIXTURE_DIR"/aimd.XV "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Basic run: all 3 analyses, --save-data ---
echo -e "\n--- Testing a basic run (--label aimd --save-data) ---"

stb-aimdAnalysis --label aimd --save-data --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "Lattice read per-frame from aimd.out" log_basic.txt
check_contains "Composition       : O2" log_basic.txt
check_contains "RDF computed over 5 frame(s)" log_basic.txt
check_contains "Diffusion coefficient" log_basic.txt
check_contains "VDOS computed" log_basic.txt
check_success aimd_rdf.png
check_success aimd_rdf.dat
check_success aimd_msd.png
check_success aimd_msd.dat
check_success aimd_vdos.png
check_success aimd_vdos.dat

echo "Testing: RDF first peak lands near the real O-O bond length (~1.0-1.3 Ang)"
python3 -c "
import numpy as np, sys
data = np.loadtxt('aimd_rdf.dat')
r, g = data[:,0], data[:,1]
peak_r = r[np.nanargmax(g)]
sys.exit(0 if 1.0 <= peak_r <= 1.3 else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} RDF first peak in the expected O-O bond-length range"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} RDF first peak outside the expected range"
    FAIL=$((FAIL+1))
fi

echo "Testing: MSD/VDOS data are finite (no NaN/Inf)"
python3 -c "
import numpy as np, sys
msd = np.loadtxt('aimd_msd.dat')[:,1]
vdos = np.loadtxt('aimd_vdos.dat')[:,1]
sys.exit(0 if np.all(np.isfinite(msd)) and np.all(np.isfinite(vdos)) else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} MSD/VDOS data are all finite"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} MSD/VDOS data contain NaN/Inf"
    FAIL=$((FAIL+1))
fi


# --- 3. --pair (same-species RDF) ---
echo -e "\n--- Testing --pair O-O ---"
stb-aimdAnalysis --label aimd --pair O-O --no-intro > log_pair.txt 2>&1
check_exit_code $? 0
check_contains "RDF pair          : O-O" log_pair.txt
check_success aimd_rdf_OO.png


# --- 4. --pair with a species not present ---
echo -e "\n--- Testing --pair with a nonexistent species (expects failure) ---"
stb-aimdAnalysis --label aimd --pair O-H --no-intro > log_badpair.txt 2>&1
check_exit_code $? 1
check_contains "species not found" log_badpair.txt


# --- 5. --stride / --skip ---
echo -e "\n--- Testing --stride 2 --skip 1 ---"
stb-aimdAnalysis --label aimd --stride 2 --skip 1 --no-intro > log_stride.txt 2>&1
check_exit_code $? 0
check_contains "stride 2, skip 1" log_stride.txt


# --- 6. --skip too large ---
echo -e "\n--- Testing --skip larger than the available frames (expects failure) ---"
stb-aimdAnalysis --label aimd --skip 999 --no-intro > log_badskip.txt 2>&1
check_exit_code $? 1
check_contains "discards all" log_badskip.txt


# --- 7. --save-report ---
echo -e "\n--- Testing --save-report ---"
rm -f stb_aimdAnalysis_report.txt
stb-aimdAnalysis --label aimd --save-report --no-intro > log_savereport.txt 2>&1
check_success stb_aimdAnalysis_report.txt
check_contains "STB-AIMDANALYSIS REPORT" stb_aimdAnalysis_report.txt
rm -f stb_aimdAnalysis_report.txt


# --- 7b. --trajectory: generic ASE-readable input (e.g. stb-mlmd's own output) ---
echo -e "\n--- Testing --trajectory (generic ASE trajectory, independent of SIESTA/.ANI) ---"
python3 -c "
from ase import Atoms
from ase.io import write
import numpy as np
rng = np.random.default_rng(0)
cell = [[10.0, 0, 0], [0, 10.0, 0], [0, 0, 10.0]]
frames = []
for i in range(6):
    pos = np.array([[5.0, 5.0, 4.395], [5.0, 5.0, 5.605]]) + rng.normal(scale=0.05, size=(2, 3))
    a = Atoms(symbols=['O', 'O'], positions=pos, cell=cell, pbc=True)
    a.info['Time'] = i * 15.0
    frames.append(a)
write('synthetic_traj.xyz', frames, format='extxyz')
"
check_success synthetic_traj.xyz

echo "Testing: auto-detected dt from embedded 'Time' info (no --dt needed)"
stb-aimdAnalysis --trajectory synthetic_traj.xyz --no-intro > log_trajectory.txt 2>&1
check_exit_code $? 0
check_contains "Auto-detected dt = 15.0000 fs" log_trajectory.txt
check_success synthetic_traj_rdf.png
check_success synthetic_traj_msd.png
check_success synthetic_traj_vdos.png

echo "Testing: --trajectory + --label are mutually exclusive"
stb-aimdAnalysis --trajectory synthetic_traj.xyz --label aimd --no-intro > log_mutex.txt 2>&1
check_exit_code $? 2
check_contains "mutually exclusive" log_mutex.txt

echo "Testing: --dt overrides/is required for a format without embedded Time"
python3 -c "
from ase.io import read, write
frames = read('synthetic_traj.xyz', index=':')
write('synthetic_traj.xsf', frames, format='xsf')
"
stb-aimdAnalysis --trajectory synthetic_traj.xsf --no-intro > log_missing_dt.txt 2>&1
check_exit_code $? 1
check_contains "\-\-dt is required" log_missing_dt.txt
stb-aimdAnalysis --trajectory synthetic_traj.xsf --dt 15 --no-intro > log_with_dt.txt 2>&1
check_exit_code $? 0


# --- 8. Robustness / errors ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: missing --label and --trajectory"
stb-aimdAnalysis --no-intro > log_missing_label.txt 2>&1
check_exit_code $? 2
check_contains "one of --label or --trajectory is required" log_missing_label.txt

echo "Testing: missing .ANI file"
stb-aimdAnalysis --label nope --no-intro > log_missing_ani.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_ani.txt

echo "Testing: --version"
stb-aimdAnalysis --version > log_version.txt 2>&1
check_contains "stb-aimdAnalysis" log_version.txt

echo "Testing: --help documents --pair, --skip, --fit-start, --save-data, --trajectory, --dt"
stb-aimdAnalysis --help > log_help.txt 2>&1
check_contains "pair" log_help.txt
check_contains "skip" log_help.txt
check_contains "fit-start" log_help.txt
check_contains "save-data" log_help.txt
check_contains "trajectory" log_help.txt
check_contains "\-\-dt" log_help.txt


# --- 9. Interactive path (stb-suite, shortcut 3.18) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 3.18) ---"
rm -f aimd_rdf.png aimd_msd.png aimd_vdos.png
printf '3.18\n\naimd\n1\n0\n\nn\nn\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success aimd_rdf.png
check_success aimd_msd.png
check_success aimd_vdos.png


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
