#!/bin/bash

# --- Setup ---
# Smoke test for stb-convergenceAnalysis (Convergence Test Analysis, item 4.5.2)
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
echo "--- Starting tester for STB-ConvergenceAnalysis (item 4.5.2) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR/convergence_runs"
for d in 150.0000 200.0000 250.0000 300.0000; do
    cp -r "$FIXTURE_DIR/convergence_meshcutoff_$d" "$TEST_DIR/convergence_runs/"
done
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Convergence analysis (synthetic FreeEng+SCF-converged calc.out
#     fixtures -- richer than the FreeEng-only style used by
#     1-strain/2-elastic/3-cohesive, since this tool also gates on SCF
#     convergence) ---
echo -e "\n--- Testing analysis of a 4-point meshcutoff sweep ---"
stb-convergenceAnalysis --dir convergence_runs --tolerance 0.001 --no-intro > log_analysis.txt 2>&1
check_contains "Runs found:.*4" log_analysis.txt
check_contains "Converged at:.*meshcutoff = 300.0000" log_analysis.txt
check_success convergence_curve.dat
check_success convergence_curve.gplot
check_success convergence_report.txt
check_contains "Converged at: meshcutoff = 300.0000" convergence_report.txt


# --- 3. Tighter tolerance moves the converged point ---
echo -e "\n--- Testing a tolerance no point satisfies ---"
rm -f convergence_curve.dat convergence_report.txt
stb-convergenceAnalysis --dir convergence_runs --tolerance 0.0000001 --no-intro > log_tight.txt 2>&1
check_contains "No point met the" log_tight.txt
check_contains "No point met the" convergence_report.txt


# --- 3b. A folder from a different sweep in the same --dir is skipped, not
#     silently merged into the curve (both stb-convergence invocations
#     default to the same convergence_runs/ output dir) ---
echo -e "\n--- Testing that a mismatched-parameter folder is skipped ---"
mkdir -p convergence_runs/convergence_kgrid_0.2000
cp convergence_runs/convergence_meshcutoff_150.0000/structure.fdf convergence_runs/convergence_kgrid_0.2000/
echo "siesta: FreeEng =    -999.000000" > convergence_runs/convergence_kgrid_0.2000/calc.out
stb-convergenceAnalysis --dir convergence_runs --tolerance 0.001 --no-intro > log_mixed.txt 2>&1
check_contains "different parameter: 'kgrid', expected 'meshcutoff'" log_mixed.txt
check_contains "Runs found:.*4" log_mixed.txt
rm -rf convergence_runs/convergence_kgrid_0.2000


# --- 3c. A run whose SCF cycle never confirmed convergence is flagged and
#     cannot be the reported converged point ---
echo -e "\n--- Testing that an unconfirmed-SCF run is excluded from the result ---"
mkdir -p scf_runs/convergence_meshcutoff_100.0000 scf_runs/convergence_meshcutoff_200.0000
cp convergence_runs/convergence_meshcutoff_150.0000/structure.fdf scf_runs/convergence_meshcutoff_100.0000/
cp convergence_runs/convergence_meshcutoff_150.0000/structure.fdf scf_runs/convergence_meshcutoff_200.0000/
echo "siesta: FreeEng =    -500.850000" > scf_runs/convergence_meshcutoff_100.0000/calc.out
echo "siesta: FreeEng =    -500.850001" > scf_runs/convergence_meshcutoff_200.0000/calc.out
stb-convergenceAnalysis --dir scf_runs --tolerance 0.001 --no-intro > log_scf.txt 2>&1
check_contains "never confirmed SCF convergence" log_scf.txt
check_contains "No point met the" log_scf.txt


# --- 3d. --apply writes the converged value back into a target calc.fdf ---
echo -e "\n--- Testing --apply patches a calc.fdf in place ---"
printf 'Mesh.CutOff   150.0000 Ry\nPAO.EnergyShift 0.0200 Ry\n' > apply_target.fdf
stb-convergenceAnalysis --dir convergence_runs --tolerance 0.001 --apply apply_target.fdf --no-intro > log_apply.txt 2>&1
check_contains "Applied.*meshcutoff = 300.0000" log_apply.txt
check_contains "Mesh.CutOff   300.0000 Ry" apply_target.fdf
check_contains "PAO.EnergyShift 0.0200 Ry" apply_target.fdf


# --- 3e. --apply on a kgrid sweep (never exercised before this review --
#     goes through substitute_kgrid_tag instead of substitute_numeric_tag,
#     which needs a structure to compute Monkhorst-Pack divisions from) ---
echo -e "\n--- Testing --apply on a kgrid sweep ---"
mkdir -p kgrid_runs/convergence_kgrid_0.1000 kgrid_runs/convergence_kgrid_0.2000
cp convergence_runs/convergence_meshcutoff_150.0000/structure.fdf kgrid_runs/convergence_kgrid_0.1000/
cp convergence_runs/convergence_meshcutoff_150.0000/structure.fdf kgrid_runs/convergence_kgrid_0.2000/
printf 'SCF cycle converged after 10 iterations\nsiesta: FreeEng =    -500.850000\n' > kgrid_runs/convergence_kgrid_0.1000/calc.out
printf 'SCF cycle converged after 9 iterations\nsiesta: FreeEng =    -500.850050\n' > kgrid_runs/convergence_kgrid_0.2000/calc.out
printf 'kgrid.MonkhorstPack   [1  1  1]\n' > kgrid_apply_target.fdf
stb-convergenceAnalysis --dir kgrid_runs --tolerance 0.001 --apply kgrid_apply_target.fdf --no-intro > log_kgrid_apply.txt 2>&1
check_contains "Converged at:.*kgrid = 0.1000" log_kgrid_apply.txt
check_contains "Applied.*kgrid = 0.1000" log_kgrid_apply.txt
# Bulk Si, 5.43 Ang cubic cell, no vacuum: reciprocal length ~1.157/Ang,
# density 0.1000 -> ceil(1.157/0.1) = 12 divisions per axis.
check_contains "kgrid.MonkhorstPack   \[12  12  12\]" kgrid_apply_target.fdf


# --- 3f. A tied parameter count (e.g. 1 meshcutoff + 1 kgrid folder) is
#     flagged instead of silently picked by alphabetical folder order ---
echo -e "\n--- Testing that a parameter-count tie is flagged ---"
mkdir -p tie_runs/convergence_meshcutoff_150.0000
cp convergence_runs/convergence_meshcutoff_150.0000/structure.fdf tie_runs/convergence_meshcutoff_150.0000/
cp convergence_runs/convergence_meshcutoff_150.0000/calc.out tie_runs/convergence_meshcutoff_150.0000/
cp -r kgrid_runs/convergence_kgrid_0.1000 tie_runs/
stb-convergenceAnalysis --dir tie_runs --tolerance 0.001 --no-intro > log_tie.txt 2>&1
check_contains "Ambiguous parameter mix: kgrid, meshcutoff" log_tie.txt


# --- 3g. --apply against a read-only target fails cleanly (no traceback) --
#     calc.fdf templates are sometimes chmod'd read-only on purpose ---
echo -e "\n--- Testing --apply against a read-only target ---"
printf 'Mesh.CutOff   150.0000 Ry\n' > readonly_target.fdf
chmod 444 readonly_target.fdf
stb-convergenceAnalysis --dir convergence_runs --tolerance 0.001 --apply readonly_target.fdf --no-intro > log_apply_readonly.txt 2>&1
check_contains "\[ERROR\] Could not write" log_apply_readonly.txt
chmod 644 readonly_target.fdf


# --- 4. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: empty directory (no convergence_* folders)"
mkdir -p empty_runs
stb-convergenceAnalysis --dir empty_runs --no-intro > log_empty.txt 2>&1
check_exit_code $? 1
check_contains "No 'convergence_\*' folders found" log_empty.txt

echo "Testing: missing directory"
stb-convergenceAnalysis --dir does_not_exist --no-intro > log_missing_dir.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_dir.txt

echo "Testing: folder missing calc.out (no valid data)"
mkdir -p bad_runs/convergence_meshcutoff_100.0000
cp convergence_runs/convergence_meshcutoff_150.0000/structure.fdf bad_runs/convergence_meshcutoff_100.0000/
stb-convergenceAnalysis --dir bad_runs --no-intro > log_bad_runs.txt 2>&1
check_exit_code $? 1
check_contains "No valid data found" log_bad_runs.txt

echo "Testing: --version"
stb-convergenceAnalysis --version > log_version.txt 2>&1
check_contains "stb-convergenceAnalysis" log_version.txt

echo "Testing: --help documents tolerance/dir"
stb-convergenceAnalysis --help > log_help.txt 2>&1
check_contains "tolerance" log_help.txt
check_contains "dir" log_help.txt


# --- 5. Interactive path (stb-suite, shortcut 3.5.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 3.5.2) ---"

echo "Testing: navigate 3.5.2 -> defaults -> quit"
rm -f convergence_curve.dat convergence_report.txt
printf '4.5.2\nconvergence_runs\n\n0.001\n\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_contains "Converged at:.*meshcutoff = 300.0000" log_menu.txt
check_success convergence_curve.dat

echo "Testing: navigate 3.5.2 -> custom curve name + --apply -> quit"
rm -f convergence_curve.dat convergence_report.txt my_curve.dat
printf 'Mesh.CutOff   150.0000 Ry\n' > menu_apply_target.fdf
printf '4.5.2\nconvergence_runs\n\n0.001\nmy_curve.dat\nmenu_apply_target.fdf\n\n0\n' | stb-suite > log_menu_apply.txt 2>&1
check_contains "Applied.*meshcutoff = 300.0000" log_menu_apply.txt
check_success my_curve.dat
check_contains "Mesh.CutOff   300.0000 Ry" menu_apply_target.fdf


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
