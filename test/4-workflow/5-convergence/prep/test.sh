#!/bin/bash

# --- Setup ---
# Smoke test for stb-convergence (Convergence Test Prep, item 4.5.1)
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
echo "--- Starting tester for STB-Convergence prep (item 4.5.1) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/structure.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/calc.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Single-parameter sweep, explicit range, pseudopotentials, full relaxation forcing ---
echo -e "\n--- Testing --parameter meshcutoff (--meshcutoff-range, -p dojo, custom --relax-steps) ---"
rm -rf convergence_runs
stb-convergence -s structure.fdf -c calc.fdf -p dojo --parameter meshcutoff --meshcutoff-range 150 300 50 --relax-steps 200 --no-intro > log_meshcutoff.txt 2>&1
check_contains "4 folder(s) written under 'convergence_runs' (1 parameter(s))" log_meshcutoff.txt
check_contains "RELAXATION & PARAMETER ENFORCEMENT" log_meshcutoff.txt
check_contains "PSEUDOPOTENTIALS" log_meshcutoff.txt
check_success convergence_runs/meshcutoff/convergence_meshcutoff_250.0000/calc.fdf
check_contains "%include config_extra.fdf" convergence_runs/meshcutoff/convergence_meshcutoff_250.0000/calc.fdf
check_success convergence_runs/meshcutoff/convergence_meshcutoff_250.0000/structure.fdf
check_success convergence_runs/meshcutoff/convergence_meshcutoff_250.0000/config_extra.fdf
check_contains "Mesh.CutOff           250.0000  Ry" convergence_runs/meshcutoff/convergence_meshcutoff_250.0000/config_extra.fdf
check_contains "MD.TypeOfRun       CG" convergence_runs/meshcutoff/convergence_meshcutoff_250.0000/config_extra.fdf
check_contains "MD.Steps.*200" convergence_runs/meshcutoff/convergence_meshcutoff_250.0000/config_extra.fdf
check_contains "MD.VariableCell    true" convergence_runs/meshcutoff/convergence_meshcutoff_250.0000/config_extra.fdf
if [ -L convergence_runs/meshcutoff/convergence_meshcutoff_250.0000/Si.psml ]; then
    echo -e "   -> ${GREEN}Verified:${NC} Si.psml linked into the generated folder"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} Si.psml was not linked into the generated folder"
    FAIL=$((FAIL+1))
fi


# --- 3. All 3 parameters in ONE invocation, only meshcutoff customized -- the others
#     must still show up in the SAME report with their own suggested default ---
echo -e "\n--- Testing --parameter all with only --meshcutoff-range customized (energyshift/kgrid still shown, suggested defaults) ---"
rm -rf convergence_runs
stb-convergence -s structure.fdf -c calc.fdf -p dojo --parameter all --meshcutoff-range 150 300 50 --no-intro > log_mixed.txt 2>&1
check_contains "Parameter(s)      : meshcutoff, energyshift, kgrid" log_mixed.txt
check_contains "Relax steps       : 100" log_mixed.txt
check_contains "meshcutoff  | custom (--meshcutoff-range)" log_mixed.txt
check_contains "energyshift | suggested default" log_mixed.txt
check_contains "kgrid       | suggested default" log_mixed.txt
check_contains "16 folder(s) written under 'convergence_runs' (3 parameter(s))" log_mixed.txt
check_success convergence_runs/meshcutoff/convergence_meshcutoff_250.0000/config_extra.fdf
check_success convergence_runs/energyshift/convergence_energyshift_0.0010/config_extra.fdf
check_success convergence_runs/kgrid/convergence_kgrid_0.0500/config_extra.fdf

echo "Testing: [5]'s Folder column is just the leaf name, not the full path"
check_contains "convergence_meshcutoff_250.0000 | meshcutoff" log_mixed.txt
if grep -q "convergence_runs/meshcutoff/convergence_meshcutoff_250.0000 " log_mixed.txt; then
    echo -e "   -> ${RED}Failed:${NC} report table still shows the full path"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} report table does not repeat the full path"
    PASS=$((PASS+1))
fi


# --- 4. Suggested default range (no --<parameter>-range), default relax-steps ---
echo -e "\n--- Testing --parameter energyshift (suggested default range, default --relax-steps) ---"
rm -rf convergence_runs
stb-convergence -s structure.fdf -c calc.fdf -p dojo --parameter energyshift --no-intro > log_energyshift.txt 2>&1
check_contains "energyshift | suggested default | 0.001 to 0.05, step 0.01 (Ry)" log_energyshift.txt
check_contains "6 folder(s) written under 'convergence_runs' (1 parameter(s))" log_energyshift.txt
check_success convergence_runs/energyshift/convergence_energyshift_0.0010/calc.fdf
check_contains "PAO.EnergyShift     0.0010 Ry" convergence_runs/energyshift/convergence_energyshift_0.0010/config_extra.fdf
check_contains "MD.Steps.*100" convergence_runs/energyshift/convergence_energyshift_0.0010/config_extra.fdf


# --- 5. K-grid density sweep (cross-checked against stb-kgrid directly) ---
echo -e "\n--- Testing --parameter kgrid ---"
rm -rf convergence_runs
stb-convergence -s structure.fdf -c calc.fdf -p dojo --parameter kgrid --kgrid-range 0.10 0.20 0.05 --no-intro > log_kgrid.txt 2>&1
check_contains "3 folder(s) written under 'convergence_runs' (1 parameter(s))" log_kgrid.txt
check_contains "kgrid.MonkhorstPack   \[8  8  8\]" convergence_runs/kgrid/convergence_kgrid_0.1500/config_extra.fdf

stb-kgrid -f structure.fdf -d 0.15 --no-intro > log_kgrid_crosscheck.txt 2>&1
check_contains "Suggested Monkhorst-Pack grid" log_kgrid_crosscheck.txt
check_contains "8 8 8" log_kgrid_crosscheck.txt


# --- 6. --parameter all -- 3 independent per-parameter subfolders, all default ranges ---
echo -e "\n--- Testing --parameter all (default ranges, all 3 parameters) ---"
rm -rf convergence_runs
stb-convergence -s structure.fdf -c calc.fdf -p dojo --parameter all --save-report --no-intro > log_all.txt 2>&1
check_contains "Parameter(s)      : meshcutoff, energyshift, kgrid" log_all.txt
check_contains "19 folder(s) written under 'convergence_runs' (3 parameter(s))" log_all.txt
check_success convergence_runs/meshcutoff/convergence_meshcutoff_100.0000/calc.fdf
check_success convergence_runs/energyshift/convergence_energyshift_0.0010/calc.fdf
check_success convergence_runs/kgrid/convergence_kgrid_0.0500/calc.fdf
check_success convergence_runs/convergence_stage1.txt
check_contains "stb-convergenceAnalysis (Stage 2) doesn't yet support this per-parameter layout" log_all.txt


# --- 7. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: calc.fdf WITHOUT the swept tag still works (forced via config_extra.fdf, no longer needs to pre-exist)"
cat > calc_notag.fdf << 'EOF'
SystemLabel siesta
%include structure.fdf
EOF
stb-convergence -s structure.fdf -c calc_notag.fdf -p dojo --parameter meshcutoff --meshcutoff-range 150 300 50 --no-intro > log_notag.txt 2>&1
check_exit_code $? 0
check_contains "Mesh.CutOff           250.0000  Ry" convergence_runs/meshcutoff/convergence_meshcutoff_250.0000/config_extra.fdf

echo "Testing: --meshcutoff-range with step 0"
stb-convergence -s structure.fdf -c calc.fdf --parameter meshcutoff --meshcutoff-range 150 300 0 --no-intro > log_zero_step.txt 2>&1
check_exit_code $? 2

echo "Testing: --meshcutoff-range with max <= min"
stb-convergence -s structure.fdf -c calc.fdf --parameter meshcutoff --meshcutoff-range 300 150 50 --no-intro > log_bad_range.txt 2>&1
check_exit_code $? 2

echo "Testing: --relax-steps 0"
stb-convergence -s structure.fdf -c calc.fdf --parameter meshcutoff --relax-steps 0 --no-intro > log_zero_relax.txt 2>&1
check_exit_code $? 2
check_contains "relax-steps must be" log_zero_relax.txt

echo "Testing: --kgrid-range given but 'kgrid' not in --parameter (mismatched range flag)"
stb-convergence -s structure.fdf -c calc.fdf --parameter meshcutoff --kgrid-range 0.1 0.2 0.05 --no-intro > log_mismatched.txt 2>&1
check_exit_code $? 2
check_contains "kgrid-range was given but 'kgrid' is not in --parameter" log_mismatched.txt

echo "Testing: unrecognized --parameter value"
stb-convergence -s structure.fdf -c calc.fdf --parameter bogus --no-intro > log_bad_param.txt 2>&1
check_exit_code $? 2
check_contains "unrecognized value 'bogus'" log_bad_param.txt

echo "Testing: missing structure file"
stb-convergence -s does_not_exist.fdf -c calc.fdf --parameter meshcutoff --no-intro > log_missing.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing.txt

echo "Testing: missing required --parameter"
stb-convergence -s structure.fdf -c calc.fdf --no-intro > log_missing_args.txt 2>&1
check_exit_code $? 2

echo "Testing: -s/-c default to structure.fdf/calc.fdf when omitted"
rm -rf convergence_runs
stb-convergence --parameter meshcutoff --no-intro > log_default_files.txt 2>&1
check_exit_code $? 0
check_contains "Structure file    : structure.fdf" log_default_files.txt
check_contains "Calc template     : calc.fdf" log_default_files.txt

echo "Testing: --version"
stb-convergence --version > log_version.txt 2>&1
check_contains "stb-convergence" log_version.txt

echo "Testing: --help documents the 3 parameters, per-parameter range flags, pseudo-dir/save-report/relax-steps, -s/-c defaults"
stb-convergence --help > log_help.txt 2>&1
check_contains "meshcutoff" log_help.txt
check_contains "kgrid" log_help.txt
check_contains "meshcutoff-range" log_help.txt
check_contains "energyshift-range" log_help.txt
check_contains "kgrid-range" log_help.txt
check_contains "pseudo-dir" log_help.txt
check_contains "save-report" log_help.txt
check_contains "relax-steps" log_help.txt
check_contains "default: structure.fdf" log_help.txt
check_contains "default: calc.fdf" log_help.txt


# --- 8. Interactive path (stb-suite, shortcut 4.5.1) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.5.1) ---"

echo "Testing: navigate 4.5.1 -> dojo bank -> 'all' shortcut -> no customization for any parameter -> default relax-steps -> single combined run -> save report (y) -> quit"
rm -rf convergence_runs
printf '4.5.1\nstructure.fdf\ncalc.fdf\n1\nall\nn\nn\nn\n\n\n\ny\n\n0\n' | stb-suite > log_menu_all.txt 2>&1
check_contains "Selected: meshcutoff, energyshift, kgrid" log_menu_all.txt
check_contains "Relax steps       : 100" log_menu_all.txt
check_contains "19 folder(s) written under 'convergence_runs' (3 parameter(s))" log_menu_all.txt
check_success convergence_runs/convergence_stage1.txt

echo "Testing: navigate 4.5.1 -> single parameter (meshcutoff) -> customize range -> default relax-steps -> quit"
rm -rf convergence_runs
printf '4.5.1\nstructure.fdf\ncalc.fdf\n1\n1\ny\n150\n300\n50\n\n\n\n\n0\n' | stb-suite > log_menu_single.txt 2>&1
check_contains "Selected: meshcutoff" log_menu_single.txt
check_contains "4 folder(s) written under 'convergence_runs' (1 parameter(s))" log_menu_single.txt
check_contains "custom (--meshcutoff-range)" log_menu_single.txt

echo "Testing: navigate 4.5.1 -> 'all' shortcut -> customize ONLY meshcutoff, defaults for the other two -> ONE combined run, ALL 3 shown together -> quit"
rm -rf convergence_runs
printf '4.5.1\nstructure.fdf\ncalc.fdf\n1\nall\ny\n120\n280\n40\nn\nn\n\n\n\n\n\n0\n' | stb-suite > log_menu_mixed.txt 2>&1
check_contains "Parameter(s)      : meshcutoff, energyshift, kgrid" log_menu_mixed.txt
check_contains "custom (--meshcutoff-range)" log_menu_mixed.txt
check_contains "energyshift | suggested default" log_menu_mixed.txt
check_contains "kgrid       | suggested default" log_menu_mixed.txt
check_contains "17 folder(s) written under 'convergence_runs' (3 parameter(s))" log_menu_mixed.txt
check_success convergence_runs/meshcutoff/convergence_meshcutoff_120.0000/calc.fdf
check_success convergence_runs/energyshift/convergence_energyshift_0.0010/calc.fdf
check_success convergence_runs/kgrid/convergence_kgrid_0.0500/calc.fdf

echo "Testing: blank structure/calc default to structure.fdf/calc.fdf; an unrecognized token in the parameter list is flagged, falls back to meshcutoff"
rm -rf convergence_runs
printf '4.5.1\n\n\n\nbogus\nn\n\n\n\n\n0\n' | stb-suite > log_menu_typo.txt 2>&1
check_contains "Structure file    : structure.fdf" log_menu_typo.txt
check_contains "Calc template     : calc.fdf" log_menu_typo.txt
check_contains "Ignored unrecognized option(s): bogus" log_menu_typo.txt
check_contains "Selected: meshcutoff" log_menu_typo.txt


popd > /dev/null

# --- 9. Summary ---
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
