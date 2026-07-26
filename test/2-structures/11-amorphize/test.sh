#!/bin/bash

# --- Setup ---
# Smoke test for stb-amorphize (Amorphous Structure Generator, item 2.11)
#
# Needs the optional 'ml' extra (pip install stb_suite[ml] -- PyTorch +
# mace-torch). The whole file is skipped with a clear message if `mace`
# isn't importable, same pattern as 8-defect/15-mlrelax's tests.
#
# Uses a tiny 8-atom cell with deliberately tiny --melt-steps/
# --quench-steps so this smoke test finishes quickly -- it only checks
# that the tool runs correctly end-to-end and produces valid output, not
# that the result is a high-quality amorphous structure (that physics was
# validated separately during planning with a realistic 64-atom cell and
# step counts; see CLAUDE.md's stb-amorphize entry).
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


echo "--- Starting tester for STB-Amorphize (item 2.11) ---"

if ! python3 -c "import mace" 2>/dev/null; then
    echo -e "${YELLOW}Skipped entirely:${NC} the optional 'ml' extra is not installed."
    echo "Install with: pip install stb_suite[ml]  (then re-run this test)"
    exit 0
fi


# --- 1. Preparation ---
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/si8.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Happy path: quick melt-quench on si8.fdf ---
echo -e "\n--- Testing default protocol (tiny step counts for speed) ---"
rm -f amorphous.fdf
stb-amorphize -f si8.fdf --melt-steps 20 --quench-steps 30 --no-intro \
    > log_default.txt 2>&1
check_exit_code $? 0
check_contains "\[1\] INPUT STRUCTURE" log_default.txt
check_contains "Bond-angle mean/std : " log_default.txt
check_contains "Steps used" log_default.txt
check_success amorphous.fdf


# --- 3. --no-final-relax skips the cleanup relax ---
echo -e "\n--- Testing --no-final-relax ---"
rm -f noreelax.fdf
stb-amorphize -f si8.fdf --melt-steps 10 --quench-steps 10 --no-final-relax \
    -o noreelax.fdf --no-intro > log_norelax.txt 2>&1
check_exit_code $? 0
if grep -q "FINAL STATIC RELAX" log_norelax.txt; then
    echo -e "   -> ${RED}Failed:${NC} '[4] FINAL STATIC RELAX' should NOT appear with --no-final-relax"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} '[4] FINAL STATIC RELAX' correctly absent"
    PASS=$((PASS+1))
fi
check_success noreelax.fdf


# --- 4. Rejected on a vacuum-padded slab (bulk-only tool) ---
echo -e "\n--- Testing rejection on a vacuum-padded slab ---"
stb-slab -f si8.fdf --hkl 1 0 0 -o slab_for_test.fdf --no-intro > /dev/null 2>&1
stb-amorphize -f slab_for_test.fdf --no-intro > log_vacuum.txt 2>&1
check_exit_code $? 1
check_contains "only supports bulk" log_vacuum.txt
rm -f slab_for_test.fdf


# --- 5. --help / --version ---
echo -e "\n--- Testing --help / --version ---"
stb-amorphize --version > log_version.txt 2>&1
check_contains "stb-amorphize" log_version.txt

stb-amorphize --help > log_help.txt 2>&1
check_contains "melt-temp" log_help.txt
check_contains "quench-steps" log_help.txt
check_contains "no-final-relax" log_help.txt
check_contains "taut" log_help.txt
check_contains "taup" log_help.txt
check_contains "compressibility" log_help.txt
check_contains "Bulk (3D periodic) structures only" log_help.txt
check_contains "save-data" log_help.txt
check_contains "save-traj" log_help.txt
check_contains "traj-format" log_help.txt
check_contains "stride" log_help.txt


# --- 6. stb-standard report: full symmetry, --save-report, provenance header, taut/taup fix ---
echo -e "\n--- Testing the standardized report (full symmetry, --save-report, provenance header) ---"

echo "Testing: --save-report writes the full report to file, including the (expected-to-differ) symmetry table"
rm -f stb_amorphize_report.txt
stb-amorphize -f si8.fdf --melt-steps 10 --quench-steps 10 --save-report -o report_test.fdf --no-intro \
    > log_save_report.txt 2>&1
check_exit_code $? 0
check_success stb_amorphize_report.txt
check_contains "\[6\] SYMMETRY ANALYSIS (BEFORE / AFTER)" stb_amorphize_report.txt
check_contains "EXPECTED to differ" stb_amorphize_report.txt
check_contains "Point Group" stb_amorphize_report.txt

echo "Testing: provenance header written into the output .fdf"
check_contains "# Amorphous structure generated by stb-amorphize" report_test.fdf
check_contains "# Protocol: melt" report_test.fdf

echo "Testing: --taut/--taup/--compressibility overrides are accepted and recorded"
stb-amorphize -f si8.fdf --melt-steps 5 --quench-steps 5 --taut 30 --taup 150 \
    --compressibility 1e-5 -o custom_barostat.fdf --no-intro > log_barostat.txt 2>&1
check_exit_code $? 0
check_contains "taut/taup      : 30.0 fs / 150.0 fs" log_barostat.txt
check_contains "Compressibility: 1e-05 eV/Ang\^3" log_barostat.txt
check_success custom_barostat.fdf


# --- 6b. --save-data / --save-traj: MD diagnostics (.dat + gnuplot) and OVITO/VMD trajectory ---
echo -e "\n--- Testing --save-data / --save-traj ---"

echo "Testing: --save-data writes the diagnostics .dat + .gplot"
rm -f savedata.fdf savedata_md_diagnostics.dat savedata_md_diagnostics.gplot savedata_md_diagnostics.pdf
stb-amorphize -f si8.fdf --melt-steps 20 --quench-steps 20 --stride 5 --save-data \
    -o savedata.fdf --no-intro > log_savedata.txt 2>&1
check_exit_code $? 0
check_success savedata_md_diagnostics.dat
check_success savedata_md_diagnostics.gplot
check_contains "step  time_fs  temp_K  epot_eV  ekin_eV  etot_eV  volume_Ang3" savedata_md_diagnostics.dat
check_contains "index 0: melt" savedata_md_diagnostics.dat
check_contains "index 1: quench" savedata_md_diagnostics.dat
check_contains "MD diagnostics" log_savedata.txt

echo "Testing: rendering the .gplot with gnuplot produces a PDF"
if command -v gnuplot > /dev/null 2>&1; then
    gnuplot savedata_md_diagnostics.gplot > log_gnuplot.txt 2>&1
    check_exit_code $? 0
    check_success savedata_md_diagnostics.pdf
else
    echo -e "   -> ${YELLOW}Skipped:${NC} gnuplot not installed, skipping actual render"
fi

echo "Testing: --save-traj writes a multi-frame trajectory (default xsf format)"
rm -f savetraj.fdf savetraj_md_traj.xsf
stb-amorphize -f si8.fdf --melt-steps 10 --quench-steps 10 --stride 5 --save-traj \
    -o savetraj.fdf --no-intro > log_savetraj.txt 2>&1
check_exit_code $? 0
check_success savetraj_md_traj.xsf
check_contains "MD trajectory" log_savetraj.txt

echo "Testing: --traj-format pdb/xyz"
rm -f trajpdb.fdf trajpdb_md_traj.pdb trajxyz.fdf trajxyz_md_traj.xyz
stb-amorphize -f si8.fdf --melt-steps 10 --quench-steps 10 --stride 5 --save-traj \
    --traj-format pdb -o trajpdb.fdf --no-intro > log_trajpdb.txt 2>&1
check_exit_code $? 0
check_success trajpdb_md_traj.pdb
stb-amorphize -f si8.fdf --melt-steps 10 --quench-steps 10 --stride 5 --save-traj \
    --traj-format xyz -o trajxyz.fdf --no-intro > log_trajxyz.txt 2>&1
check_exit_code $? 0
check_success trajxyz_md_traj.xyz


# --- 7. Interactive path (stb-suite, shortcut 2.11) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 2.11) ---"

echo "Testing: navigate 2.11 -> invalid file then valid -> tiny steps -> no save-data/save-traj -> final relax yes -> default model -> no save-report/view -> default output -> quit"
rm -f amorphous.fdf
menu_inputs=(
    "2.11" "does_not_exist.fdf" "si8.fdf" "3000" "15" "300" "20" ""
    "" "" "" "" "" "" ""
    "0"
)
printf '%s\n' "${menu_inputs[@]}" | stb-suite > log_menu.txt 2>&1
check_contains "File not found" log_menu.txt
check_contains "Bond-angle mean/std : " log_menu.txt
check_success amorphous.fdf

echo "Testing: navigate 2.11 -> tiny steps -> save-data yes -> save-traj yes -> xyz format -> stride 5 -> final relax yes -> default model -> no save-report/view -> custom output -> quit"
rm -f interactive_traj.fdf interactive_traj_md_diagnostics.dat interactive_traj_md_traj.xyz
menu_inputs=(
    "2.11" "si8.fdf" "3000" "10" "300" "10" "y" "y" "xyz" "5"
    "" "" "" "" "" "interactive_traj.fdf"
    "0"
)
printf '%s\n' "${menu_inputs[@]}" | stb-suite > log_menu_traj.txt 2>&1
check_contains "Bond-angle mean/std : " log_menu_traj.txt
check_success interactive_traj.fdf
check_success interactive_traj_md_diagnostics.dat
check_success interactive_traj_md_traj.xyz


popd > /dev/null

# --- 8. Summary ---
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
