#!/bin/bash

# --- Setup ---
# Smoke test for stb-ani2traj (AIMD Trajectory Converter, item 6.5).
#
# Fixtures: aimd.ANI/.XV/.fdf/.out are from a real SIESTA AIMD run (5-step
# Verlet MD on an isolated O2 molecule in a 10x10x10 Ang box -- see aimd.fdf
# for the exact input). aimd.ANI itself has no lattice information (SIESTA's
# format doesn't carry one). aimd.out has a real "Begin MD step = N" /
# "outcell:" / "siesta: E_KS(eV) =" / "siesta: Temp_ion =" block for each of
# the 5 steps (SIESTA prints outcell every step regardless of whether the
# cell actually varies) -- stb-ani2traj prefers this per-frame source over
# the single static lattice from aimd.XV/aimd.fdf when available.
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$FIXTURE_DIR/test_files"

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

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
    if grep -q "$1" "$2" 2>/dev/null; then
        echo -e "   -> ${GREEN}Verified:${NC} '$1' found in '$2'"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} '$1' NOT found in '$2'"
        FAIL=$((FAIL+1))
    fi
}

check_not_contains() {
    if grep -q "$1" "$2" 2>/dev/null; then
        echo -e "   -> ${RED}Failed:${NC} '$1' unexpectedly found in '$2'"
        FAIL=$((FAIL+1))
    else
        echo -e "   -> ${GREEN}Verified:${NC} '$1' NOT found in '$2'"
        PASS=$((PASS+1))
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
echo "--- Starting tester for STB-ANI2TRAJ (item 6.5) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/aimd.ANI" "$TEST_DIR/"
cp "$FIXTURE_DIR/aimd.XV" "$TEST_DIR/"
cp "$FIXTURE_DIR/aimd.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/aimd.out" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. XSF (default format), per-frame lattice from .out (preferred over .XV/.fdf) ---
echo -e "\n--- Testing --out-format xsf (default), per-frame lattice from .out ---"
stb-ani2traj --no-intro --label aimd > log_xsf.txt 2>&1
check_contains "Lattice read per-frame from aimd.out (5 MD step(s) found)" log_xsf.txt
check_contains "Read 5 frame(s)" log_xsf.txt
check_success aimd_traj.xsf
check_contains "ANIMSTEPS 5" aimd_traj.xsf
check_contains "PRIMVEC" aimd_traj.xsf
echo "Testing: energy/temperature enrichment is skipped (not embeddable) for xsf"
check_contains "only embedded for --out-format xyz" log_xsf.txt


# --- 3. PDB ---
echo -e "\n--- Testing --out-format pdb ---"
stb-ani2traj --no-intro --label aimd --out-format pdb > log_pdb.txt 2>&1
check_success aimd_traj.pdb
check_contains "CRYST1" aimd_traj.pdb
check_contains "MODEL" aimd_traj.pdb


# --- 4. Extended XYZ, with real per-frame E_KS/Temp_ion/Time enrichment ---
# SIESTA itself reports (aimd.out): step 1 E_KS=-856.3116 eV, Temp_ion=300.000 K;
# step 5 E_KS=-857.5875 eV, Temp_ion=9917.825 K. MD.LengthTimeStep is 1.0 fs
# (aimd.fdf), so step N's Time should be N*1.0 fs.
echo -e "\n--- Testing --out-format xyz with real per-frame E_KS/Temp_ion/Time ---"
stb-ani2traj --no-intro --label aimd --out-format xyz > log_xyz.txt 2>&1
check_success aimd_traj.xyz
check_contains "Lattice=" aimd_traj.xyz
check_contains "E_KS=-856.3116 Temp_ion=300.0 Time=1.0" aimd_traj.xyz
check_contains "E_KS=-857.5875 Temp_ion=9917.825 Time=5.0" aimd_traj.xyz


# --- 5. --stride, staying aligned with the per-frame .out data ---
echo -e "\n--- Testing --stride 2 (5 frames -> 3, keeping steps 1/3/5) ---"
stb-ani2traj --no-intro --label aimd --out-format xyz --stride 2 -o strided.xyz > log_stride.txt 2>&1
check_contains "Read 3 frame(s)" log_stride.txt
check_contains "E_KS=-856.3116" strided.xyz
check_contains "E_KS=-856.8561" strided.xyz
check_contains "E_KS=-857.5875" strided.xyz
check_not_contains "E_KS=-856.4574" strided.xyz


# --- 6. --unwrap ---
echo -e "\n--- Testing --unwrap ---"
stb-ani2traj --no-intro --label aimd --out-format xyz --unwrap -o unwrapped.xyz > log_unwrap.txt 2>&1
check_contains "Unwrapped periodic-boundary jumps" log_unwrap.txt
check_success unwrapped.xyz


# --- 7. Custom -o ---
echo -e "\n--- Testing custom -o output filename ---"
stb-ani2traj --no-intro --label aimd --out-format pdb -o custom_name.pdb > log_custom.txt 2>&1
check_success custom_name.pdb


# --- 8. Lattice source priority: .out missing -> falls back to .XV ---
echo -e "\n--- Testing fallback to .XV when .out is absent ---"
mv aimd.out aimd.out.bak
stb-ani2traj --no-intro --label aimd --out-format xsf -o xv_fallback.xsf > log_xv_fallback.txt 2>&1
check_contains "Lattice read from aimd.XV" log_xv_fallback.txt
check_success xv_fallback.xsf


# --- 9. Lattice source priority: .out and .XV missing -> falls back to .fdf ---
echo -e "\n--- Testing fallback to .fdf when .out and .XV are both absent ---"
mv aimd.XV aimd.XV.bak
stb-ani2traj --no-intro --label aimd --out-format xsf -o fdf_fallback.xsf > log_fdf_fallback.txt 2>&1
check_contains "Lattice read from aimd.fdf" log_fdf_fallback.txt
check_success fdf_fallback.xsf
mv aimd.XV.bak aimd.XV
mv aimd.out.bak aimd.out


# --- 10. Missing .ANI ---
echo -e "\n--- Testing missing input file ---"
stb-ani2traj --no-intro --label doesnotexist --out-format xsf > log_missing.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing.txt


# --- 11. Missing lattice source entirely ---
echo -e "\n--- Testing missing .out/.XV/.fdf (no lattice source at all) ---"
mkdir -p no_lattice
cp aimd.ANI no_lattice/
(cd no_lattice && stb-ani2traj --no-intro --label aimd --out-format xsf > ../log_nolattice.txt 2>&1)
check_exit_code $? 1
check_contains "a lattice is required" log_nolattice.txt
rm -rf no_lattice


# --- 12. Interactive path (stb-suite, shortcut 6.5) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 6.5) ---"
rm -f aimd_traj.xsf interactive_traj.xsf
# Prompts in order: label, format choice, stride, output filename
printf '6.5\naimd\n1\n1\ninteractive_traj.xsf\n\n0\n' | stb-suite > log_interactive.txt 2>&1
check_success interactive_traj.xsf
check_contains "Converting AIMD trajectory" log_interactive.txt


popd > /dev/null

# --- 13. Summary ---
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
