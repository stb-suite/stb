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
# The 4 committed 'convergence_meshcutoff_*' fixtures carry both a
# 'siesta: FreeEng' line AND a real 'outcell: Unit cell vectors' block, with
# energy and relaxed-cell-volume trends that DELIBERATELY disagree: energy
# looks converged by meshcutoff=250, but the relaxed cell volume only
# stabilizes at 300 -- the exact physical scenario stb-convergence's own
# module docstring warns about (energy converging before the relaxed
# geometry does), used below to exercise the dual energy+structure
# convergence criteria and their cross-check note.
echo "--- Starting tester for STB-ConvergenceAnalysis (item 4.5.2) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR/convergence_runs/meshcutoff"
for d in 150.0000 200.0000 250.0000 300.0000; do
    cp -r "$FIXTURE_DIR/convergence_meshcutoff_$d" "$TEST_DIR/convergence_runs/meshcutoff/"
done
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null
STRUCT="convergence_runs/meshcutoff/convergence_meshcutoff_150.0000/structure.fdf"


# --- 2. Basic nested-layout analysis: energy vs. structure disagreement ---
echo -e "\n--- Testing nested-layout analysis of a 4-point meshcutoff sweep ---"
stb-convergenceAnalysis --dir convergence_runs --no-intro > log_basic.txt 2>&1
check_contains "Folders found: 4" log_basic.txt
check_contains "Energy converged at:.*meshcutoff = 250.0000" log_basic.txt
check_contains "Structure converged at:.*meshcutoff = 300.0000" log_basic.txt
check_contains "Energy and relaxed structure disagree" log_basic.txt
if grep -q "NOT CONVERGED" log_basic.txt; then
    echo -e "   -> ${RED}Failed:${NC} unexpected NOT CONVERGED for a resolved sweep"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no NOT CONVERGED (recommended value was found)"
    PASS=$((PASS+1))
fi


# --- 3. Tolerance overrides ---
echo -e "\n--- Testing a very tight energy tolerance (energy WARNING, structure still resolves) ---"
stb-convergenceAnalysis --dir convergence_runs --tolerance 0.0000001 --no-intro > log_tight_energy.txt 2>&1
check_contains "No point met the 1e-07 eV/atom energy tolerance" log_tight_energy.txt
check_contains "Structure converged at:.*meshcutoff = 300.0000" log_tight_energy.txt

echo -e "\n--- Testing a very tight volume tolerance (structure WARNING) ---"
stb-convergenceAnalysis --dir convergence_runs --volume-tolerance 0.0000001 --no-intro > log_tight_volume.txt 2>&1
check_contains "No point met the 1e-07% relaxed-volume tolerance" log_tight_volume.txt

echo -e "\n--- Testing a volume tolerance wide enough that energy and structure AGREE ---"
stb-convergenceAnalysis --dir convergence_runs --volume-tolerance 2.0 --no-intro > log_agree.txt 2>&1
check_contains "Structure converged at:.*meshcutoff = 250.0000" log_agree.txt
check_contains "Energy and structure agree on the same converged value." log_agree.txt


# --- 4. Data-quality flags (SCF not confirmed, residual force, step-limited) ---
echo -e "\n--- Testing data-quality flags: SCF?, F>tol, steps ---"
mkdir -p quality_runs/meshcutoff/convergence_meshcutoff_100.0000
cp "$STRUCT" quality_runs/meshcutoff/convergence_meshcutoff_100.0000/
cat > quality_runs/meshcutoff/convergence_meshcutoff_100.0000/calc.out <<'EOF'
siesta: FreeEng =    -500.900000

outcell: Unit cell vectors (Ang):
        5.35    0.000000    0.000000
        0.000000    5.35    0.000000
        0.000000    0.000000    5.35
EOF

mkdir -p quality_runs/meshcutoff/convergence_meshcutoff_200.0000
cp "$STRUCT" quality_runs/meshcutoff/convergence_meshcutoff_200.0000/
cat > quality_runs/meshcutoff/convergence_meshcutoff_200.0000/calc.out <<'EOF'
SCF cycle converged after 10 iterations
siesta: FreeEng =    -500.850000

outcell: Unit cell vectors (Ang):
        5.42    0.000000    0.000000
        0.000000    5.42    0.000000
        0.000000    0.000000    5.42

siesta: Atomic forces (eV/Ang):
     1    0.500000    0.500000    0.500000
   Max    0.500000
EOF

mkdir -p quality_runs/meshcutoff/convergence_meshcutoff_300.0000
cp "$STRUCT" quality_runs/meshcutoff/convergence_meshcutoff_300.0000/
cat > quality_runs/meshcutoff/convergence_meshcutoff_300.0000/calc.out <<'EOF'
SCF cycle converged after 10 iterations
siesta: FreeEng =    -500.845000

outcell: Unit cell vectors (Ang):
        5.45    0.000000    0.000000
        0.000000    5.45    0.000000
        0.000000    0.000000    5.45

                        Begin CG opt. move =      0
                        Begin CG opt. move =      1
                        Begin CG opt. move =      2
EOF
printf 'MD.NumCGsteps      3\n' > quality_runs/meshcutoff/convergence_meshcutoff_300.0000/config_extra.fdf

stb-convergenceAnalysis --dir quality_runs --no-intro > log_quality.txt 2>&1
check_contains "SCF?" log_quality.txt
check_contains "F>tol" log_quality.txt
check_contains "steps" log_quality.txt
check_contains "3/3 run(s) flagged" log_quality.txt

echo -e "\n--- Testing --force-tolerance override clears the F>tol flag ---"
stb-convergenceAnalysis --dir quality_runs --no-intro --force-tolerance 10 > log_force_tol.txt 2>&1
check_contains "200.0000 | -250.425000 | 0.025000       | 159.2201      | 3.9768         | OK" log_force_tol.txt


# --- 5. Flat layout (pointing --dir directly at one parameter's own subfolder) ---
echo -e "\n--- Testing flat layout (--dir pointed directly at the meshcutoff subfolder) ---"
stb-convergenceAnalysis --dir convergence_runs/meshcutoff --no-intro > log_flat.txt 2>&1
check_contains "Energy converged at:.*meshcutoff = 250.0000" log_flat.txt


# --- 6. Multi-parameter nested discovery (all 3 swept parameters in one --dir) ---
echo -e "\n--- Testing multi-parameter discovery (meshcutoff + energyshift + kgrid together) ---"
mkdir -p convergence_runs/energyshift/convergence_energyshift_0.0050
cp "$STRUCT" convergence_runs/energyshift/convergence_energyshift_0.0050/
printf 'SCF cycle converged after 10 iterations\nsiesta: FreeEng =    -500.850000\n' \
    > convergence_runs/energyshift/convergence_energyshift_0.0050/calc.out
mkdir -p convergence_runs/energyshift/convergence_energyshift_0.0200
cp "$STRUCT" convergence_runs/energyshift/convergence_energyshift_0.0200/
printf 'SCF cycle converged after 9 iterations\nsiesta: FreeEng =    -500.850050\n' \
    > convergence_runs/energyshift/convergence_energyshift_0.0200/calc.out

mkdir -p convergence_runs/kgrid/convergence_kgrid_0.1000
cp "$STRUCT" convergence_runs/kgrid/convergence_kgrid_0.1000/
printf 'SCF cycle converged after 10 iterations\nsiesta: FreeEng =    -500.850000\n' \
    > convergence_runs/kgrid/convergence_kgrid_0.1000/calc.out
mkdir -p convergence_runs/kgrid/convergence_kgrid_0.2000
cp "$STRUCT" convergence_runs/kgrid/convergence_kgrid_0.2000/
printf 'SCF cycle converged after 9 iterations\nsiesta: FreeEng =    -500.850050\n' \
    > convergence_runs/kgrid/convergence_kgrid_0.2000/calc.out

stb-convergenceAnalysis --dir convergence_runs --no-intro > log_multi.txt 2>&1
check_contains "\[2.1\] MESHCUTOFF SWEEP" log_multi.txt
check_contains "\[2.2\] ENERGYSHIFT SWEEP" log_multi.txt
check_contains "\[2.3\] KGRID SWEEP" log_multi.txt
check_contains "Energy converged at:.*energyshift = 0.0050" log_multi.txt
check_contains "Energy converged at:.*kgrid = 0.1000" log_multi.txt


# --- 7. --save-report / --save-gnuplot / --output-dir ---
echo -e "\n--- Testing --save-report/--save-gnuplot/--output-dir ---"
rm -rf saved_out
stb-convergenceAnalysis --dir convergence_runs --no-intro --save-report --save-gnuplot \
    --output-dir saved_out > log_saved.txt 2>&1
check_success saved_out/convergence_report.txt
check_success saved_out/meshcutoff_convergence.dat
check_success saved_out/meshcutoff_convergence.gplot
check_success saved_out/energyshift_convergence.dat
check_success saved_out/kgrid_convergence.dat
check_contains "Energy converged at: meshcutoff = 250.0000" saved_out/convergence_report.txt

echo "Testing: the .gplot script actually renders with gnuplot"
if command -v gnuplot > /dev/null 2>&1; then
    ( cd saved_out && gnuplot meshcutoff_convergence.gplot )
    check_success saved_out/meshcutoff_convergence.pdf
else
    echo -e "   -> ${YELLOW}Skipped:${NC} gnuplot not installed"
fi


# --- 8. --apply (multi-parameter nested tree, kgrid via substitute_kgrid_tag) ---
echo -e "\n--- Testing --apply writes every converged parameter into one calc.fdf ---"
printf 'Mesh.CutOff   150.0000 Ry\nPAO.EnergyShift 0.0500 Ry\nkgrid.MonkhorstPack   [1  1  1]\n' \
    > apply_target.fdf
stb-convergenceAnalysis --dir convergence_runs --apply apply_target.fdf --no-intro > log_apply.txt 2>&1
check_contains "\[Applied\].*meshcutoff = 300.0000" log_apply.txt
check_contains "\[Applied\].*energyshift = 0.0050" log_apply.txt
check_contains "\[Applied\].*kgrid = 0.1000" log_apply.txt
check_contains "Mesh.CutOff   300.0000 Ry" apply_target.fdf
check_contains "PAO.EnergyShift 0.0050 Ry" apply_target.fdf
# Bulk Si, 5.43 Ang cubic cell (from convergence_meshcutoff_150.0000/structure.fdf),
# no vacuum: density 0.1000 -> ceil(reciprocal length / 0.1) = 12 divisions/axis.
check_contains "kgrid.MonkhorstPack   \[12  12  12\]" apply_target.fdf

echo -e "\n--- Testing --apply against a read-only target fails cleanly (no traceback) ---"
printf 'Mesh.CutOff   150.0000 Ry\n' > readonly_target.fdf
chmod 444 readonly_target.fdf
stb-convergenceAnalysis --dir convergence_runs --apply readonly_target.fdf --no-intro > log_apply_readonly.txt 2>&1
check_contains "\[ERROR\] Could not write" log_apply_readonly.txt
chmod 644 readonly_target.fdf


# --- 9. Error and robustness cases ---
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

echo "Testing: folder missing calc.out is skipped, not a hard failure (partial sweep)"
mkdir -p bad_runs/meshcutoff/convergence_meshcutoff_100.0000
cp "$STRUCT" bad_runs/meshcutoff/convergence_meshcutoff_100.0000/
stb-convergenceAnalysis --dir bad_runs --no-intro > log_bad_runs.txt 2>&1
check_exit_code $? 0
check_contains "No valid data found for this parameter" log_bad_runs.txt
check_contains "NOT CONVERGED" log_bad_runs.txt

echo "Testing: --version"
stb-convergenceAnalysis --version > log_version.txt 2>&1
check_contains "stb-convergenceAnalysis" log_version.txt

echo "Testing: --help documents tolerance/volume-tolerance/force-tolerance/save-gnuplot/view/apply"
stb-convergenceAnalysis --help > log_help.txt 2>&1
check_contains "tolerance" log_help.txt
check_contains "volume-tolerance" log_help.txt
check_contains "force-tolerance" log_help.txt
check_contains "save-gnuplot" log_help.txt
check_contains "\-\-view" log_help.txt
check_contains "apply" log_help.txt


# --- 10. Interactive path (stb-suite, shortcut 4.5.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.5.2) ---"

echo "Testing: navigate 4.5.2 -> all defaults (no advanced settings, no save/view/apply) -> quit"
printf '4.5.2\n\n\n\nn\n\nn\nn\nn\n\n\n0\n' | stb-suite > log_menu_defaults.txt 2>&1
check_contains "Energy converged at:.*meshcutoff = 250.0000" log_menu_defaults.txt

echo "Testing: navigate 4.5.2 -> advanced settings + save-report/gnuplot + apply -> quit"
rm -rf menu_out
printf 'Mesh.CutOff   150.0000 Ry\n' > menu_apply_target.fdf
printf '4.5.2\n\n\n\ny\n2.0\n0.05\nmenu_out\ny\ny\nn\nmenu_apply_target.fdf\n\n\n0\n' \
    | stb-suite > log_menu_advanced.txt 2>&1
check_contains "Volume tolerance  : 2.0 %" log_menu_advanced.txt
check_contains "Energy and structure agree on the same converged value." log_menu_advanced.txt
check_success menu_out/convergence_report.txt
check_success menu_out/meshcutoff_convergence.dat
check_contains "Applied.*meshcutoff = 250.0000" log_menu_advanced.txt
check_contains "Mesh.CutOff   250.0000 Ry" menu_apply_target.fdf


popd > /dev/null

# --- 11. Summary ---
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
