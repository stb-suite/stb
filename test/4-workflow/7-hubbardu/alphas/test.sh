#!/bin/bash

# --- Setup ---
# Smoke test for stb-hubbarduAlphas (Hubbard U Linear Response, Stage 2: Perturbation Prep, item 4.7.2)
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

check_absent() {
    if [ ! -e "$1" ]; then
        echo -e "   -> ${GREEN}Verified:${NC} '$1' absent, as expected"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} '$1' should not exist"
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

# Simulates the user having run SIESTA in reference/ to convergence: stage 1
# (stb-hubbardu, tested separately) already produces a correct reference/
# folder -- we only need to fake the .DM file it would have written.
make_reference() {
    local outdir="$1"
    rm -rf "$outdir"
    stb-hubbardu -s structure.fdf -c calc.fdf --species Mn --output-dir "$outdir" --no-intro > /dev/null 2>&1
    echo "fake density matrix placeholder" > "$outdir/reference/siesta.DM"
}


# --- 1. Preparation ---
echo "--- Starting tester for STB-HubbarduAlphas stage 2: perturbation prep (item 4.7.2) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/structure.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/structure_2mn.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/calc.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Missing-.DM guard (reference/ generated but SIESTA never run) ---
echo -e "\n--- Testing the missing-.DM guard ---"
rm -rf hubbardu_runs
stb-hubbardu -s structure.fdf -c calc.fdf --species Mn --no-intro > /dev/null 2>&1
stb-hubbarduAlphas --dir hubbardu_runs --no-intro > log_no_dm.txt 2>&1
check_exit_code $? 1
check_contains "run SIESTA in 'hubbardu_runs/reference/' first" log_no_dm.txt


# --- 3. Normal generation (default alphas) ---
echo -e "\n--- Testing default perturbation-sweep generation ---"
make_reference hubbardu_runs
stb-hubbarduAlphas --dir hubbardu_runs --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "Generated 9 run(s)" log_basic.txt
check_success hubbardu_runs/scf_alpha_0.1000/calc.fdf
check_success hubbardu_runs/frozen_alpha_0.1000/calc.fdf

echo "Testing: calc.fdf just includes config_extra.fdf now (DFT+U/MaxSCFIterations/etc. live there)"
check_contains "%include config_extra.fdf" hubbardu_runs/scf_alpha_0.1000/calc.fdf
check_success hubbardu_runs/scf_alpha_0.1000/config_extra.fdf

echo "Testing: the frozen branch gets its OWN alpha=0.0000 folder (not a reuse of reference)"
check_success hubbardu_runs/frozen_alpha_0.0000/calc.fdf
check_contains "0.000    0.000" hubbardu_runs/frozen_alpha_0.0000/config_extra.fdf

echo "Testing: frozen folders use MaxSCFIterations 2 (not the old buggy 1)"
check_contains "MaxSCFIterations 2" hubbardu_runs/frozen_alpha_0.1000/config_extra.fdf
if grep -q "MaxSCFIterations 1$" hubbardu_runs/frozen_alpha_0.1000/config_extra.fdf; then
    echo -e "   -> ${RED}Failed:${NC} should not use MaxSCFIterations 1"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no stray MaxSCFIterations 1"
    PASS=$((PASS+1))
fi

echo "Testing: the reference .DM was auto-copied into every scf/frozen folder (no manual step)"
check_success hubbardu_runs/scf_alpha_0.1000/siesta.DM
check_success hubbardu_runs/frozen_alpha_0.1000/siesta.DM
check_success hubbardu_runs/frozen_alpha_0.0000/siesta.DM

echo "Testing: scf folders get DM.UseSaveDM but no MaxSCFIterations cap"
check_contains "DM.UseSaveDM T" hubbardu_runs/scf_alpha_0.1000/config_extra.fdf
if grep -q "MaxSCFIterations" hubbardu_runs/scf_alpha_0.1000/config_extra.fdf; then
    echo -e "   -> ${RED}Failed:${NC} scf_alpha_0.1000 should NOT have MaxSCFIterations"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} scf_alpha_0.1000 has no MaxSCFIterations override"
    PASS=$((PASS+1))
fi

echo "Testing: run_manifest.json now has reference + 4 scf + 5 frozen = 10 entries"
N_RUNS=$(python3 -c "import json; print(len(json.load(open('hubbardu_runs/run_manifest.json'))['runs']))")
if [ "$N_RUNS" -eq 10 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} 10 total run entries in manifest"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} expected 10 run entries, got $N_RUNS"
    FAIL=$((FAIL+1))
fi


# --- 4. Explicit --alphas override ---
echo -e "\n--- Testing --alphas override ---"
make_reference hubbardu_runs2
stb-hubbarduAlphas --dir hubbardu_runs2 --alphas -0.2 0.2 --no-intro > log_override.txt 2>&1
check_exit_code $? 0
check_contains "Generated 5 run(s)" log_override.txt
check_success hubbardu_runs2/scf_alpha_0.2000/calc.fdf
check_success hubbardu_runs2/frozen_alpha_0.0000/calc.fdf


# --- 4b. --frozen-iterations override ---
echo -e "\n--- Testing --frozen-iterations ---"
make_reference hubbardu_runs5
stb-hubbarduAlphas --dir hubbardu_runs5 --alphas 0.1 --frozen-iterations 3 --no-intro > log_frozen_iter.txt 2>&1
check_exit_code $? 0
check_contains "MaxSCFIterations 3" hubbardu_runs5/frozen_alpha_0.1000/config_extra.fdf

echo "Testing: --frozen-iterations below the empirical minimum warns but still runs"
rm -rf hubbardu_runs6
make_reference hubbardu_runs6
stb-hubbarduAlphas --dir hubbardu_runs6 --alphas 0.1 --frozen-iterations 1 --no-intro > log_frozen_iter_low.txt 2>&1
check_exit_code $? 0
check_contains "below the empirically-needed minimum" log_frozen_iter_low.txt
check_contains "MaxSCFIterations 1" hubbardu_runs6/frozen_alpha_0.1000/config_extra.fdf

echo "Testing: --frozen-iterations 0 is rejected"
stb-hubbarduAlphas --dir hubbardu_runs6 --alphas 0.1 --frozen-iterations 0 --no-intro > log_frozen_iter_zero.txt 2>&1
check_exit_code $? 1
check_contains "must be >= 1" log_frozen_iter_zero.txt


# --- 4c. Rounding-collision guard (--alphas 0.00001 rounds to the same folder as the implicit zero point) ---
echo -e "\n--- Testing the alpha rounding-collision guard ---"
make_reference hubbardu_runs7
stb-hubbarduAlphas --dir hubbardu_runs7 --alphas 0.00001 0.1 --no-intro > log_round_zero.txt 2>&1
check_exit_code $? 1
check_contains "an alpha rounds to 0.0000 eV" log_round_zero.txt

stb-hubbarduAlphas --dir hubbardu_runs7 --alphas 0.10001 0.1 --no-intro > log_round_collide.txt 2>&1
check_exit_code $? 1
check_contains "round to the same folder name" log_round_collide.txt


# --- 4d. Pseudopotentials saved by stage 1 are auto-copied into every folder ---
echo -e "\n--- Testing pseudopotential auto-copy ---"
mkdir -p pseudos
echo "fake Mn pseudopotential" > pseudos/Mn.psml
rm -rf hubbardu_runs8
stb-hubbardu -s structure.fdf -c calc.fdf --species Mn --pseudo-dir pseudos --output-dir hubbardu_runs8 --no-intro > /dev/null 2>&1
echo "fake density matrix placeholder" > hubbardu_runs8/reference/siesta.DM
stb-hubbarduAlphas --dir hubbardu_runs8 --alphas 0.1 --no-intro > log_pseudo.txt 2>&1
check_exit_code $? 0
check_success hubbardu_runs8/scf_alpha_0.1000/Mn.psml
check_success hubbardu_runs8/frozen_alpha_0.0000/Mn.psml
check_success hubbardu_runs8/frozen_alpha_0.1000/Mn.psml

echo "Testing: no pseudopotentials saved -> a clear reminder note is printed"
stb-hubbarduAlphas --dir hubbardu_runs2 --alphas 0.05 --no-intro > log_no_pseudo.txt 2>&1
check_contains "No pseudopotentials saved by stage 1" log_no_pseudo.txt


# --- 4e. Stage 2 also rejects a hand-edited template that already has DFT+U set up ---
echo -e "\n--- Testing the DFT+U guard is also enforced in stage 2 ---"
make_reference hubbardu_runs9
# reference/calc.fdf is now just '%include config_extra.fdf' (the DFT+U block itself lives in
# the sibling config_extra.fdf, not inline) -- simulates a user hand-editing the template
# snapshot back to a previously-generated run folder's own calc.fdf by mistake.
cat "hubbardu_runs9/reference/calc.fdf" > "hubbardu_runs9/_template_calc.fdf"
stb-hubbarduAlphas --dir hubbardu_runs9 --alphas 0.1 --no-intro > log_stage2_guard.txt 2>&1
check_exit_code $? 1
check_contains "already '%include config_extra.fdf'" log_stage2_guard.txt


# --- 4f. Aliased (multi-atom-species) reference propagates through stage 2 ---
echo -e "\n--- Testing stage 2 with an aliased perturbed_species (from stb-hubbardu --atom-index) ---"
rm -rf hubbardu_runs10
mkdir -p pseudos_2mn
echo "fake Mn pseudopotential" > pseudos_2mn/Mn.psml
stb-hubbardu -s structure_2mn.fdf -c calc.fdf --species Mn --atom-index 1 --pseudo-dir pseudos_2mn \
    --output-dir hubbardu_runs10 --no-intro > /dev/null 2>&1
echo "fake density matrix placeholder" > hubbardu_runs10/reference/siesta.DM
stb-hubbarduAlphas --dir hubbardu_runs10 --alphas 0.1 --no-intro > log_alias.txt 2>&1
check_exit_code $? 0
check_contains "isolated via alias species 'Mn_pert'" log_alias.txt

echo "Testing: generated scf/frozen folders perturb 'Mn_pert', not 'Mn'"
check_contains "Mn_pert   1" hubbardu_runs10/scf_alpha_0.1000/config_extra.fdf
check_contains "Mn_pert   1" hubbardu_runs10/frozen_alpha_0.1000/config_extra.fdf
check_contains "Mn_pert   1" hubbardu_runs10/frozen_alpha_0.0000/config_extra.fdf
if grep -q "^Mn   1" hubbardu_runs10/scf_alpha_0.1000/config_extra.fdf; then
    echo -e "   -> ${RED}Failed:${NC} the unperturbed Mn should not appear in the LDAU.proj block"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} only the aliased atom is in the LDAU.proj block"
    PASS=$((PASS+1))
fi

echo "Testing: the aliased structure.fdf (with Mn_pert) propagated into every generated folder"
check_contains "Mn_pert" hubbardu_runs10/scf_alpha_0.1000/structure.fdf
check_contains "Mn_pert" hubbardu_runs10/frozen_alpha_0.1000/structure.fdf

echo "Testing: the aliased pseudopotential file propagated without any extra stage-2 filtering"
check_success hubbardu_runs10/scf_alpha_0.1000/Mn_pert.psml
check_success hubbardu_runs10/frozen_alpha_0.1000/Mn_pert.psml


# --- 5. Error and robustness cases ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: missing directory"
stb-hubbarduAlphas --dir does_not_exist --no-intro > log_missing_dir.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_dir.txt

echo "Testing: directory without run_manifest.json"
mkdir -p not_a_hubbardu_dir
stb-hubbarduAlphas --dir not_a_hubbardu_dir --no-intro > log_no_manifest.txt 2>&1
check_exit_code $? 1
check_contains "run stb-hubbardu" log_no_manifest.txt

echo "Testing: alpha=0.0 is rejected (frozen already gets its own zero point automatically)"
make_reference hubbardu_runs3
stb-hubbarduAlphas --dir hubbardu_runs3 --alphas 0.0 0.1 --no-intro > log_zero_alpha.txt 2>&1
check_exit_code $? 1
check_contains "an alpha rounds to 0.0000 eV" log_zero_alpha.txt

echo "Testing: duplicate alphas are rejected"
stb-hubbarduAlphas --dir hubbardu_runs3 --alphas 0.1 0.1 --no-intro > log_dup_alpha.txt 2>&1
check_exit_code $? 1
check_contains "duplicate values" log_dup_alpha.txt

echo "Testing: --version"
stb-hubbarduAlphas --version > log_version.txt 2>&1
check_contains "stb-hubbarduAlphas" log_version.txt

echo "Testing: --help documents dir/alphas/frozen-iterations"
stb-hubbarduAlphas --help > log_help.txt 2>&1
check_contains "dir" log_help.txt
check_contains "alphas" log_help.txt
check_contains "frozen-iterations" log_help.txt


# --- 6. Interactive path (stb-suite, shortcut 4.7.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.7.2) ---"

echo "Testing: navigate 4.7.2 -> default dir -> default alphas -> default frozen-iterations ->"
echo "  save-report=N -> quit"
make_reference hubbardu_runs4
printf '4.7.2\nhubbardu_runs4\n\n\nn\n\n0\n' | timeout 20 stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_contains "Generated 9 run(s)" log_menu.txt
check_success hubbardu_runs4/frozen_alpha_0.0000/calc.fdf
check_absent hubbardu_runs4/hubbardu_stage2.txt

echo "Testing: navigate 4.7.2 again -> save-report=Y -> quit"
rm -rf hubbardu_runs4b
make_reference hubbardu_runs4b
printf '4.7.2\nhubbardu_runs4b\n\n\ny\n\n0\n' | timeout 20 stb-suite > log_menu_report.txt 2>&1
check_exit_code $? 0
check_contains "Generated 9 run(s)" log_menu_report.txt
check_success hubbardu_runs4b/hubbardu_stage2.txt


popd > /dev/null

# --- 7. Summary ---
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
