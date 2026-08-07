#!/bin/bash

# --- Setup ---
# Smoke test for stb-hubbardu (Hubbard U Linear Response, Stage 1: Reference Prep, item 4.7.1)
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


# --- 1. Preparation ---
echo "--- Starting tester for STB-Hubbardu stage 1: reference prep (item 4.7.1) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/structure.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/structure_2mn.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/calc.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Default reference generation (species with a known default shell) ---
echo -e "\n--- Testing reference generation for Mn (auto-detected 3d shell) ---"
stb-hubbardu -s structure.fdf -c calc.fdf --species Mn --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "Mn (3d: n=3, l=2)" log_basic.txt
check_contains "Generated 'reference/'" log_basic.txt
check_success hubbardu_runs/run_manifest.json
check_success hubbardu_runs/reference/calc.fdf
check_success hubbardu_runs/reference/structure.fdf
check_success hubbardu_runs/_template_structure.fdf
check_success hubbardu_runs/_template_calc.fdf

echo "Testing: only the reference folder is created (no scf/frozen folders yet)"
if [ -d hubbardu_runs/scf_alpha_0.1000 ] || [ -d hubbardu_runs/frozen_alpha_0.1000 ]; then
    echo -e "   -> ${RED}Failed:${NC} stage 1 should not create scf_alpha_*/frozen_alpha_* folders"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no scf_alpha_*/frozen_alpha_* folders created by stage 1"
    PASS=$((PASS+1))
fi

echo "Testing: calc.fdf now just includes config_extra.fdf, instead of embedding the DFT+U block inline"
check_contains "%include config_extra.fdf" hubbardu_runs/reference/calc.fdf

echo "Testing: the generated LDAU.proj block (in config_extra.fdf) matches the verified SIESTA syntax (U=0, J=0 default)"
check_contains "LDAU.PotentialShift T" hubbardu_runs/reference/config_extra.fdf
check_contains "Mn   1" hubbardu_runs/reference/config_extra.fdf
check_contains "n=3    2" hubbardu_runs/reference/config_extra.fdf
check_contains "0.000    0.000" hubbardu_runs/reference/config_extra.fdf

echo "Testing: config_extra.fdf also forces single-point SCF (no ionic/cell relaxation)"
check_contains "MD.TypeOfRun       CG" hubbardu_runs/reference/config_extra.fdf
check_contains "MD.Steps           0" hubbardu_runs/reference/config_extra.fdf
check_contains "MD.VariableCell    false" hubbardu_runs/reference/config_extra.fdf

echo "Testing: run_manifest.json records species/shell/j and the reference run"
check_contains '"species": "Mn"' hubbardu_runs/run_manifest.json
check_contains '"j": 0.0' hubbardu_runs/run_manifest.json
check_contains '"kind": "reference"' hubbardu_runs/run_manifest.json


# --- 3. Explicit --shell and --j override ---
echo -e "\n--- Testing --shell and --j overrides ---"
rm -rf hubbardu_runs
stb-hubbardu -s structure.fdf -c calc.fdf --species Mn --shell 4f --j 0.5 --no-intro > log_override.txt 2>&1
check_exit_code $? 0
check_contains "Mn (4f: n=4, l=3)" log_override.txt
check_contains "0.500" log_override.txt
check_contains "n=4    3" hubbardu_runs/reference/config_extra.fdf
check_contains "0.000    0.500" hubbardu_runs/reference/config_extra.fdf
check_contains '"j": 0.5' hubbardu_runs/run_manifest.json


# --- 3b. --pseudo-dir, including '~' expansion (argv never goes through a shell
#     when the interactive menu calls subprocess.run(), so os.path.expanduser
#     must do this -- 'is not a directory' otherwise even for a real path) ---
echo -e "\n--- Testing --pseudo-dir (plain path and '~' expansion) ---"
mkdir -p pseudos
echo "fake Mn pseudopotential" > pseudos/Mn.psml
rm -rf hubbardu_runs
stb-hubbardu -s structure.fdf -c calc.fdf --species Mn --pseudo-dir pseudos --no-intro > log_pseudo.txt 2>&1
check_exit_code $? 0
check_contains "Copied 1 pseudopotential file(s)" log_pseudo.txt
check_success hubbardu_runs/reference/Mn.psml
check_success hubbardu_runs/_pseudopotentials/Mn.psml

echo "Testing: no --pseudo-dir given -> a reminder is printed"
rm -rf hubbardu_runs
stb-hubbardu -s structure.fdf -c calc.fdf --species Mn --no-intro > log_no_pseudo.txt 2>&1
check_contains "No --pseudo-dir given" log_no_pseudo.txt

echo "Testing: --pseudo-dir with a literal '~' is expanded to \$HOME"
mkdir -p ~/stb_hubbardu_pseudo_test
echo "fake" > ~/stb_hubbardu_pseudo_test/Mn.psml
rm -rf hubbardu_runs
stb-hubbardu -s structure.fdf -c calc.fdf --species Mn --pseudo-dir '~/stb_hubbardu_pseudo_test' --no-intro > log_tilde.txt 2>&1
check_exit_code $? 0
check_success hubbardu_runs/reference/Mn.psml
rm -rf ~/stb_hubbardu_pseudo_test

echo "Testing: a --pseudo-dir that isn't a bundled bank or a real directory is a clean error"
stb-hubbardu -s structure.fdf -c calc.fdf --species Mn --pseudo-dir does_not_exist_dir --no-intro > log_bad_pseudo.txt 2>&1
check_exit_code $? 1
check_contains "nor an existing directory" log_bad_pseudo.txt

echo "Testing: --pseudo-dir accepts a bundled bank name (dojo)"
rm -rf hubbardu_runs
stb-hubbardu -s structure.fdf -c calc.fdf --species Mn --pseudo-dir dojo --no-intro > log_bank.txt 2>&1
check_exit_code $? 0
check_success hubbardu_runs/reference/Mn.psml

echo "Testing: using a bundled bank prints its citation/origin"
check_contains "van Setten" log_bank.txt
check_contains "pseudo-dojo.org" log_bank.txt

echo "Testing: using a custom --pseudo-dir folder (not a bundled bank) prints NO citation"
rm -rf hubbardu_runs_custom_pp
stb-hubbardu -s structure.fdf -c calc.fdf --species Mn --pseudo-dir pseudos \
    --output-dir hubbardu_runs_custom_pp --no-intro > log_custom_pp.txt 2>&1
check_exit_code $? 0
if grep -q "\[INFO\] Using bundled pseudopotential bank" log_custom_pp.txt; then
    echo -e "   -> ${RED}Failed:${NC} a custom folder should never print a bundled-bank citation"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no citation printed for a custom folder"
    PASS=$((PASS+1))
fi

echo "Testing: --pseudo-dir accepts the other bundled bank (virtual_vault) with its own citation"
rm -rf hubbardu_runs_vv
stb-hubbardu -s structure.fdf -c calc.fdf --species Mn --pseudo-dir virtual_vault \
    --output-dir hubbardu_runs_vv --no-intro > log_bank_vv.txt 2>&1
check_exit_code $? 0
check_contains "NNIN/C" log_bank_vv.txt
check_contains "nninc.cnf.cornell.edu" log_bank_vv.txt

echo "Testing: a bank has ~70-80 elements, but only the structure's own species get copied"
PP_COUNT=$(ls hubbardu_runs/reference/*.psml 2>/dev/null | wc -l)
if [ "$PP_COUNT" -eq 1 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} exactly 1 pseudopotential copied from the bank (found $PP_COUNT)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} expected exactly 1 pseudopotential copied from the bank, found $PP_COUNT"
    FAIL=$((FAIL+1))
fi
check_absent hubbardu_runs/reference/Fe.psml
check_absent hubbardu_runs/_pseudopotentials/Fe.psml


# --- 3c. Multi-atom species -- perturbation must be isolated to one atom ---
echo -e "\n--- Testing --atom-index (species appearing more than once) ---"

echo "Testing: --species matching 2 atoms, no --atom-index -> clean error listing candidates"
stb-hubbardu -s structure_2mn.fdf -c calc.fdf --species Mn --no-intro > log_multi_no_index.txt 2>&1
check_exit_code $? 1
check_contains "appears 2 times" log_multi_no_index.txt
check_contains "Pass --atom-index" log_multi_no_index.txt
check_contains "atom-index | Wyckoff | Note" log_multi_no_index.txt
check_contains "symmetry-equivalent atom(s) -- any one of them gives the same result" log_multi_no_index.txt

echo "Testing: --species matching 1 atom (existing fixture) needs no --atom-index -- unaffected"
rm -rf hubbardu_runs
stb-hubbardu -s structure.fdf -c calc.fdf --species Mn --no-intro > log_single_unaffected.txt 2>&1
check_exit_code $? 0
check_contains "Mn   1" hubbardu_runs/reference/config_extra.fdf
if grep -q "Mn_pert" hubbardu_runs/reference/config_extra.fdf; then
    echo -e "   -> ${RED}Failed:${NC} single-atom species should never be aliased"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no aliasing for a single-atom species"
    PASS=$((PASS+1))
fi

echo "Testing: --atom-index pointing at the wrong species is a clean error"
stb-hubbardu -s structure_2mn.fdf -c calc.fdf --species Mn --atom-index 3 --no-intro > log_multi_wrong_species.txt 2>&1
check_exit_code $? 1
check_contains "is species 'O', not 'Mn'" log_multi_wrong_species.txt

echo "Testing: --atom-index out of range is a clean error"
stb-hubbardu -s structure_2mn.fdf -c calc.fdf --species Mn --atom-index 99 --no-intro > log_multi_out_of_range.txt 2>&1
check_exit_code $? 1
check_contains "out of range (1 to 4)" log_multi_out_of_range.txt

echo "Testing: valid --atom-index isolates the perturbation via an aliased species"
rm -rf hubbardu_runs
stb-hubbardu -s structure_2mn.fdf -c calc.fdf --species Mn --atom-index 1 --no-intro > log_multi_ok.txt 2>&1
check_exit_code $? 0
check_contains "Perturbed atom" log_multi_ok.txt
check_contains "Mn_pert" log_multi_ok.txt
check_contains "Mn_pert   1" hubbardu_runs/reference/config_extra.fdf
check_contains "perturbed_species" hubbardu_runs/run_manifest.json
check_contains "Mn_pert" hubbardu_runs/run_manifest.json
check_contains "\"atom_index\": 1" hubbardu_runs/run_manifest.json

echo "Testing: aliased structure.fdf has 3 species (Mn, O, Mn_pert) and NumberOfSpecies bumped to 3"
check_contains "NumberOfSpecies    3" hubbardu_runs/reference/structure.fdf
check_contains "Mn_pert" hubbardu_runs/reference/structure.fdf
check_contains "Mn_pert" hubbardu_runs/_template_structure.fdf

echo "Testing: only ONE atom's species index was changed -- the other Mn line still points at species 1"
MN_LINE_COUNT=$(grep -c "   1$" hubbardu_runs/reference/structure.fdf)
if [ "$MN_LINE_COUNT" -ge 1 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} the other Mn atom still uses the original species index"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} expected the unperturbed Mn atom to keep species index 1"
    FAIL=$((FAIL+1))
fi

echo "Testing: alias species collision is a clean error"
cat > structure_2mn_collision.fdf << 'EOF'
NumberOfSpecies    3
NumberofAtoms      5
%block ChemicalSpeciesLabel
 1   25   Mn
 2   8   O
 3   25   Mn_pert
%endblock ChemicalSpeciesLabel
LatticeConstant 1.0 Ang
AtomicCoordinatesFormat  Fractional
%block LatticeVectors
 0.00000000   4.40800000   4.40800000
 2.20400000   0.00000000   2.20400000
 2.20400000   2.20400000   0.00000000
%endblock LatticeVectors
%block AtomicCoordinatesAndAtomicSpecies
  0.00000000   0.00000000   0.00000000   1
  0.50000000   0.00000000   0.00000000   1
  0.25000000   0.50000000   0.50000000   2
  0.75000000   0.50000000   0.50000000   2
  0.10000000   0.10000000   0.10000000   3
%endblock AtomicCoordinatesAndAtomicSpecies
EOF
stb-hubbardu -s structure_2mn_collision.fdf -c calc.fdf --species Mn --atom-index 1 --no-intro > log_alias_collision.txt 2>&1
check_exit_code $? 1
check_contains "already exists" log_alias_collision.txt

echo "Testing: --pseudo-dir also copies an aliased pseudopotential file"
mkdir -p pseudos_2mn
echo "fake Mn pseudopotential" > pseudos_2mn/Mn.psml
rm -rf hubbardu_runs
stb-hubbardu -s structure_2mn.fdf -c calc.fdf --species Mn --atom-index 1 --pseudo-dir pseudos_2mn --no-intro > log_multi_pseudo.txt 2>&1
check_exit_code $? 0
check_success hubbardu_runs/reference/Mn_pert.psml
check_success hubbardu_runs/_pseudopotentials/Mn_pert.psml


# --- 4. Error and robustness cases ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: species not present in the structure"
stb-hubbardu -s structure.fdf -c calc.fdf --species Fe --no-intro > log_missing_species.txt 2>&1
check_exit_code $? 1
check_contains "not found in" log_missing_species.txt

echo "Testing: species with no default shell and no --shell given"
cat > structure_si.fdf << 'EOF'
NumberOfSpecies    2
NumberofAtoms      2
%block ChemicalSpeciesLabel
 1   25   Mn
 2   14   Si
%endblock ChemicalSpeciesLabel
LatticeConstant 1.0 Ang
AtomicCoordinatesFormat  Fractional
%block LatticeVectors
 2.8900000000   0.0000000000   0.0000000000
 0.0000000000   2.8900000000   0.0000000000
 0.0000000000   0.0000000000   2.8900000000
%endblock LatticeVectors
%block AtomicCoordinatesAndAtomicSpecies
 0.000000000   0.000000000   0.000000000   1
 0.500000000   0.500000000   0.500000000   2
%endblock AtomicCoordinatesAndAtomicSpecies
EOF
stb-hubbardu -s structure_si.fdf -c calc.fdf --species Si --no-intro > log_no_shell.txt 2>&1
check_exit_code $? 1
check_contains "No default correlated shell known" log_no_shell.txt

echo "Testing: --shell explicitly supplied for a species with no default works"
stb-hubbardu -s structure_si.fdf -c calc.fdf --species Si --shell 3d -o si_runs --no-intro > log_shell_override.txt 2>&1
check_exit_code $? 0

echo "Testing: a template with an embedded LDAU.proj block/PotentialShift is rejected"
cat > already_has_ldau_inline.fdf << 'EOF'
SystemLabel siesta
%include structure.fdf
LDAU.PotentialShift T
%block LDAU.proj
Mn 1
n=3 2
0.000 0.000
0.000 0.000
%endblock LDAU.proj
EOF
stb-hubbardu -s structure.fdf -c already_has_ldau_inline.fdf --species Mn --no-intro > log_guard.txt 2>&1
check_exit_code $? 1
check_contains "already contains an LDAU.proj" log_guard.txt

echo "Testing: re-pointing --calc straight at a PREVIOUSLY-GENERATED reference/calc.fdf is rejected"
echo "  (its DFT+U block now lives in a sibling config_extra.fdf, not inline -- the raw LDAU regex"
echo "  alone can no longer catch this, so check_no_existing_ldau also flags the bare %include line)"
cp hubbardu_runs/reference/calc.fdf already_has_ldau.fdf
stb-hubbardu -s structure.fdf -c already_has_ldau.fdf --species Mn --no-intro > log_guard_include.txt 2>&1
check_exit_code $? 1
check_contains "already '%include config_extra.fdf'" log_guard_include.txt

echo "Testing: the DFTU.proj/DFTU.PotentialShift spelling is also rejected (SIESTA accepts both)"
cat > already_has_dftu.fdf << 'EOF'
SystemLabel siesta
%include structure.fdf
DFTU.PotentialShift T
%block DFTU.proj
Mn 1
n=3 2
0.000 0.000
0.000 0.000
%endblock DFTU.proj
EOF
stb-hubbardu -s structure.fdf -c already_has_dftu.fdf --species Mn --no-intro > log_guard_dftu.txt 2>&1
check_exit_code $? 1
check_contains "already contains an LDAU.proj" log_guard_dftu.txt

echo "Testing: missing structure file"
stb-hubbardu -s does_not_exist.fdf -c calc.fdf --species Mn --no-intro > log_missing.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing.txt

echo "Testing: missing required args"
stb-hubbardu -s structure.fdf -c calc.fdf --no-intro > log_missing_args.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-hubbardu --version > log_version.txt 2>&1
check_contains "stb-hubbardu" log_version.txt

echo "Testing: --help documents species/shell/j and SCF convergence tips"
stb-hubbardu --help > log_help.txt 2>&1
check_contains "species" log_help.txt
check_contains "shell" log_help.txt
check_contains "\-\-j J" log_help.txt
check_contains "DM.MixingWeight" log_help.txt
check_contains "SCF.H.Tolerance" log_help.txt


# --- 5. Interactive path (stb-suite, shortcut 4.7.1) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.7.1) ---"

echo "Testing: navigate 4.7.1 -> blank structure/calc (defaults to structure.fdf/calc.fdf in cwd) ->"
echo "  Mn -> no atom-index -> auto shell -> default J -> no pseudo -> default output -> save-report=N -> quit"
rm -rf hubbardu_runs
printf '4.7.1\n\n\nMn\n\n\n\n\n\nn\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_contains "Generated 'reference/'" log_menu.txt
check_success hubbardu_runs/run_manifest.json
check_absent hubbardu_runs/hubbardu_stage1.txt

echo "Testing: navigate 4.7.1 -> structure_2mn.fdf -> calc.fdf -> Mn -> atom-index 1 -> auto shell ->"
echo "  default J -> no pseudo -> default output -> save-report=Y -> quit"
rm -rf hubbardu_runs
printf '4.7.1\nstructure_2mn.fdf\ncalc.fdf\nMn\n1\n\n\n\n\ny\n\n0\n' | stb-suite > log_menu_atom_index.txt 2>&1
check_contains "Perturbed atom" log_menu_atom_index.txt
check_contains "Generated 'reference/'" log_menu_atom_index.txt
check_contains "Mn_pert   1" hubbardu_runs/reference/config_extra.fdf
check_success hubbardu_runs/hubbardu_stage1.txt

echo "Testing: navigate 4.7.1 with a genuinely missing structure file -> the CLI itself (not the"
echo "  interactive wrapper) reports the error, since the wrapper no longer pre-validates paths"
rm -rf hubbardu_runs
printf '4.7.1\ndoes_not_exist.fdf\ncalc.fdf\nMn\n\n\n\n\n\nn\n\n0\n' | stb-suite > log_menu_missing.txt 2>&1
check_contains "does_not_exist.fdf' not found" log_menu_missing.txt


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
