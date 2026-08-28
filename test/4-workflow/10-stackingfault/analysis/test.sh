#!/bin/bash

# --- Setup ---
# Smoke test for stb-stackingfaultAnalysis (Stacking Fault Analysis, item
# 4.10.2). Uses stb-stackingfault itself (real tool, no SIESTA needed for
# prep) to build a 3x3 grid against the same fixture as ../prep/ (homobilayer
# graphene/graphene, cheap and deterministic), then fabricates
# "siesta: FreeEng"/SCF/Max-force lines per grid point (synthetic, same
# printf style already used by 8-adsorption/9-neb's analysis test.sh) to
# exercise the analysis side without needing a real SIESTA run. Since Stage 1
# nests everything under 'sf_run/' and every grid folder under its own
# 'positions/' subfolder, every fabricated calc.out below goes through
# 'sf_run/positions/<label>/'.
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREP_DIR="$(cd "$FIXTURE_DIR/../prep" && pwd)"
TEST_DIR="$FIXTURE_DIR/test_files"

# Output colors
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

PASS=0
FAIL=0

# --- Check helpers ---

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
echo "--- Starting tester for STB-StackingfaultAnalysis (item 4.10.2) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$PREP_DIR/graphene.fdf" "$TEST_DIR/"
cp "$PREP_DIR/calc.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Build a real 3x3 grid via stb-stackingfault, fabricate FreeEng/SCF/Max ---
echo -e "\n--- Testing analysis of a 3x3 stacking-fault grid ---"
stb-stackingfault -l1 graphene.fdf -l2 graphene.fdf -c calc.fdf -nx 3 -ny 3 --mode 3 --no-intro > log_prep.txt 2>&1
check_exit_code $? 0
n_grid=$(find sf_run/positions -maxdepth 1 -type d -name 'shift_*' | wc -l)
if [ "$n_grid" -eq 9 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} prep produced 9 grid folders under sf_run/positions/"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} expected 9 grid folders, found $n_grid"
    FAIL=$((FAIL+1))
fi
check_success sf_run/sf_manifest.json

# Non-trivial energy landscape: AA (shift_00_00, eclipsed) is highest, a
# specific off-symmetry point (shift_01_02) is the minimum, one point
# (shift_02_01) is deliberately unconverged / high residual force.
# Corrugation = 1.6 eV over the fixture's 5.2408 Ang^2 in-plane cell area ->
# 305.295 meV/Ang^2 (hand-verified against a real run of this exact fixture).
declare -A energies=(
  [shift_00_00]=-198.900000 [shift_00_01]=-199.700000 [shift_00_02]=-199.750000
  [shift_01_00]=-199.720000 [shift_01_01]=-199.900000 [shift_01_02]=-200.500000
  [shift_02_00]=-199.680000 [shift_02_01]=-199.100000 [shift_02_02]=-199.730000
)
for label in "${!energies[@]}"; do
    e="${energies[$label]}"
    if [ "$label" = "shift_02_01" ]; then
        printf "siesta: FreeEng =    %s\nsiesta: Atomic forces (eV/Ang):\n   Max    0.350000\n" "$e" > "sf_run/positions/$label/calc.out"
    else
        printf "siesta: FreeEng =    %s\nSCF cycle converged after 10 iterations\nsiesta: Atomic forces (eV/Ang):\n   Max    0.010000\n" "$e" > "sf_run/positions/$label/calc.out"
    fi
done

stb-stackingfaultAnalysis --dir sf_run --save-gnuplot --no-intro > log_analysis.txt 2>&1
check_exit_code $? 0
check_contains "Equilibrium stacking (min)    : shift_01_02" log_analysis.txt
check_contains "Highest-energy registry (max) : shift_00_00" log_analysis.txt
check_contains "Corrugation (stacking-fault) energy : 305.295 meV/Ang\^2" log_analysis.txt
check_contains "Interface area (in-plane): 5.2408 Ang\^2" log_analysis.txt
check_contains "never confirmed SCF convergence.*shift_02_01" log_analysis.txt
check_contains "residual force above --force-tolerance.*shift_02_01" log_analysis.txt
check_contains "dE(meV/A2)" log_analysis.txt
check_success sf_run/plot/stackingfault_surface.dat
check_success sf_run/plot/stackingfault_surface.gplot
check_success sf_run/stackingfault_animation.xyz
check_contains "\[0\] RUN METADATA" log_analysis.txt
check_contains "\[1\] GRID ENERGIES" log_analysis.txt
check_contains "\[2\] STACKING FAULT ANALYSIS" log_analysis.txt
check_contains "\[3\] SUMMARY" log_analysis.txt
check_contains "Grid       : 3 x 3" log_analysis.txt
check_contains "set pm3d map" sf_run/plot/stackingfault_surface.gplot
check_contains "stb-stackingfaultAnalysis" sf_run/plot/stackingfault_surface.gplot
check_contains "Stacking Fault Energy (meV/A\^2)" sf_run/plot/stackingfault_surface.gplot
check_contains "dE(meV/Ang\^2)" sf_run/plot/stackingfault_surface.dat
if grep -q "plotdensity.py" sf_run/plot/stackingfault_surface.gplot; then
    echo -e "   -> ${RED}Failed:${NC} .gplot incorrectly says 'Generated by plotdensity.py'"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} .gplot does not say 'Generated by plotdensity.py'"
    PASS=$((PASS+1))
fi
# Hexagonal lattice (120 deg) -> the skewed-plane warning should fire
check_contains "cut plane is skewed" sf_run/plot/stackingfault_surface.gplot

echo "Testing: extended-XYZ animation is genuine extended XYZ (Lattice + per-frame .info), 9 frames"
python3 - <<'PYEOF'
import sys
from ase.io import read
frames = read("sf_run/stackingfault_animation.xyz", index=":")
ok = (len(frames) == 9
      and all(f.cell.volume > 0 for f in frames)
      and "dE_meV_per_A2" in frames[0].info
      and "label" in frames[0].info)
sys.exit(0 if ok else 1)
PYEOF
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} animation has 9 frames with Lattice + per-frame .info"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} animation is missing frames, Lattice, or .info metadata"
    FAIL=$((FAIL+1))
fi

echo "Testing: --save-report writes the report under <dir>/, not cwd"
rm -f sf_run/stackingfault_report.txt
stb-stackingfaultAnalysis --dir sf_run --save-report --no-intro > log_savereport.txt 2>&1
check_success sf_run/stackingfault_report.txt

echo "Testing: without --save-gnuplot, no plot/ files are (re)written but the animation still is"
rm -rf sf_run/plot sf_run/stackingfault_animation.xyz
stb-stackingfaultAnalysis --dir sf_run --no-intro > log_nognuplot.txt 2>&1
check_contains "Gnuplot data+script : not written" log_nognuplot.txt
check_success sf_run/stackingfault_animation.xyz
if [ -d sf_run/plot ]; then
    echo -e "   -> ${RED}Failed:${NC} 'sf_run/plot/' was created despite no --save-gnuplot"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no 'sf_run/plot/' created without --save-gnuplot"
    PASS=$((PASS+1))
fi
# Restore the plot/ files for the tests below that expect them present.
stb-stackingfaultAnalysis --dir sf_run --save-gnuplot --no-intro > /dev/null 2>&1

echo "Testing: --view (headless-safe, Agg backend) does not crash and writes no PNG"
MPLBACKEND=Agg stb-stackingfaultAnalysis --dir sf_run --view --no-intro > log_view.txt 2>&1
check_exit_code $? 0
if find sf_run -maxdepth 1 -name '*.png' | grep -q .; then
    echo -e "   -> ${RED}Failed:${NC} --view wrote a PNG (should be on-screen-only, per WORKFLOW_TOOLS convention)"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} --view wrote no PNG"
    PASS=$((PASS+1))
fi

echo "Testing: --view-animation without a display fails gracefully (no crash)"
DISPLAY= stb-stackingfaultAnalysis --dir sf_run --view-animation --no-intro > log_viewanim.txt 2>&1
check_exit_code $? 0
check_contains "view-animation" log_viewanim.txt


# --- 3. --apply copies the equilibrium (min-energy) point's structure.fdf ---
echo -e "\n--- Testing --apply ---"
stb-stackingfaultAnalysis --dir sf_run --apply equilibrium.fdf --no-intro > log_apply.txt 2>&1
check_contains "Applied.*shift_01_02" log_apply.txt
check_success equilibrium.fdf
check_contains "NumberofAtoms      4" equilibrium.fdf


# --- 4. A grid point missing calc.out is skipped, not fatal ---
echo -e "\n--- Testing that a grid point missing calc.out is skipped ---"
mv sf_run/positions/shift_02_02/calc.out sf_run/positions/shift_02_02/calc.out.bak
# clear the plot/ files left by the earlier complete-grid run (test 2)
# BEFORE the incomplete-grid run below, so we can tell whether THIS run
# (re)creates them, not whether they were already there from before.
rm -rf sf_run/plot
stb-stackingfaultAnalysis --dir sf_run --save-gnuplot --no-intro > log_partial.txt 2>&1
check_contains "SKIP" log_partial.txt
check_contains "skipped: 1" log_partial.txt

# --- 4b. Regression: an INCOMPLETE grid must skip the pm3d surface plot
#     entirely (not write a broken/blank one) -- gnuplot's pm3d map renders
#     the WHOLE plot blank if even a single grid point is missing (verified
#     live against real gnuplot 6.0, not just assumed), so the fix must
#     actually suppress the plot files, not just print a warning alongside
#     a broken one. The always-on animation is unaffected (fewer frames,
#     never blocked).
echo -e "\n--- Testing that an incomplete grid skips the surface plot (not a broken one) ---"
check_contains "Skipping the gamma-surface plot" log_partial.txt
if [ -d sf_run/plot ]; then
    echo -e "   -> ${RED}Failed:${NC} 'sf_run/plot/' was written despite an incomplete grid"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no 'sf_run/plot/' written for an incomplete grid"
    PASS=$((PASS+1))
fi
check_contains "Animation     -> sf_run/stackingfault_animation.xyz (8 frame(s))" log_partial.txt
mv sf_run/positions/shift_02_02/calc.out.bak sf_run/positions/shift_02_02/calc.out


# --- 5. sf_manifest.json missing entirely (fallback path) ---
echo -e "\n--- Testing fallback when sf_manifest.json is missing ---"
mv sf_run/sf_manifest.json sf_run/sf_manifest.json.bak
stb-stackingfaultAnalysis --dir sf_run --no-intro > log_fallback.txt 2>&1
check_exit_code $? 0
check_contains "No 'sf_manifest.json' found" log_fallback.txt
check_contains "Equilibrium stacking (min)    : shift_01_02" log_fallback.txt
mv sf_run/sf_manifest.json.bak sf_run/sf_manifest.json


# --- 5b. --dir smart default: 'sf_run' from the output root, '.' from
#     inside sf_run/ itself ---
echo -e "\n--- Testing --dir smart default ---"
stb-stackingfaultAnalysis --no-intro > log_smartdefault_root.txt 2>&1
check_contains "Directory  : sf_run" log_smartdefault_root.txt
check_contains "Equilibrium stacking (min)    : shift_01_02" log_smartdefault_root.txt
pushd sf_run > /dev/null
stb-stackingfaultAnalysis --no-intro > ../log_smartdefault_inside.txt 2>&1
popd > /dev/null
check_contains "Directory  : \." log_smartdefault_inside.txt
check_contains "Equilibrium stacking (min)    : shift_01_02" log_smartdefault_inside.txt


# --- 6. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: missing directory entirely"
stb-stackingfaultAnalysis --dir does_not_exist --no-intro > log_no_dir.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_no_dir.txt

echo "Testing: directory with no shift_* folders"
mkdir -p empty_dir
stb-stackingfaultAnalysis --dir empty_dir --no-intro > log_no_grid.txt 2>&1
check_exit_code $? 1
check_contains "Did you run stb-stackingfault" log_no_grid.txt

echo "Testing: --force-tolerance"
stb-stackingfaultAnalysis --dir sf_run --force-tolerance 1.0 --no-intro > log_force_tol.txt 2>&1
if grep -q "residual force above --force-tolerance.*shift_02_01" log_force_tol.txt; then
    echo -e "   -> ${RED}Failed:${NC} unexpected force warning with a loose --force-tolerance"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no force warning with --force-tolerance 1.0 (shift_02_01's 0.35 eV/Ang is under it)"
    PASS=$((PASS+1))
fi

echo "Testing: --version"
stb-stackingfaultAnalysis --version > log_version.txt 2>&1
check_contains "stb-stackingfaultAnalysis" log_version.txt

echo "Testing: --help documents --dir/--file/--apply/--save-report/--save-gnuplot/--view/--view-animation"
stb-stackingfaultAnalysis --help > log_help.txt 2>&1
check_contains "dir" log_help.txt
check_contains "file" log_help.txt
check_contains "apply" log_help.txt
check_contains "save-report" log_help.txt
check_contains "save-gnuplot" log_help.txt
check_contains "\-\-view" log_help.txt
check_contains "view-animation" log_help.txt


# --- 6b. --mode 1 grid: manifest carries mode/n_layer1_atoms, analysis
#     reports them and adds the GapFinal(A) column (falls back to the
#     pre-relax structure.fdf when no '*.XV' is present -- no real SIESTA
#     binary/sisl fixture available in this environment, so this only
#     exercises the wiring, not a genuinely-relaxed gap value) ---
echo -e "\n--- Testing --mode 1 grid is recognized/reported by the analysis stage ---"
rm -rf mode1_root
stb-stackingfault -l1 graphene.fdf -l2 graphene.fdf -c calc.fdf -nx 2 -ny 2 --mode 1 \
    -O mode1_root --no-intro > log_prep_mode1.txt 2>&1
check_exit_code $? 0
for label in shift_00_00 shift_00_01 shift_01_00 shift_01_01; do
    printf "siesta: FreeEng =    -199.000000\nSCF cycle converged after 10 iterations\nsiesta: Atomic forces (eV/Ang):\n   Max    0.010000\n" \
        > "mode1_root/sf_run/positions/$label/calc.out"
done
stb-stackingfaultAnalysis --dir mode1_root/sf_run --no-intro > log_analysis_mode1.txt 2>&1
check_exit_code $? 0
check_contains "Mode       : 1" log_analysis_mode1.txt
check_contains "GapFinal(A)" log_analysis_mode1.txt


# --- 7. Interactive path (stb-suite, shortcut 4.10.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.10.2) ---"

echo "Testing: navigate 4.10.2 -> defaults (sf_run, everything off) -> quit"
rm -rf sf_run/plot sf_run/stackingfault_animation.xyz sf_run/stackingfault_report.txt
# 4.10.2 (menu code) / "" (dir default -> sf_run) / "" (out_file default) /
# "" (force-tolerance default) / "" (apply_target: skip) / "" (save_report:
# N) / "" (save_gnuplot: N) / "" (view: N) / "" (view_animation: N) / ""
# (Press Enter to continue) / 0 (quit)
printf '4.10.2\n\n\n\n\n\n\n\n\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_contains "Equilibrium stacking (min)    : shift_01_02" log_menu.txt
check_success sf_run/stackingfault_animation.xyz
if [ -d sf_run/plot ]; then
    echo -e "   -> ${RED}Failed:${NC} 'sf_run/plot/' was created despite answering N to save_gnuplot"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} 'sf_run/plot/' not created when save_gnuplot answered N"
    PASS=$((PASS+1))
fi


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
