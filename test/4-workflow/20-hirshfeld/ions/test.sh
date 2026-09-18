#!/bin/bash

# --- Setup ---
# Smoke test for stb-hirshfeldIons (Hirshfeld-I Stage 2: Ion Prep, item
# 4.20.2). Chains a real stb-hirshfeldPrep (Stage 1, tested separately) run
# on the ../prep fixture, then fabricates each neutral/<species>/'s own
# .RHO via ../make_synthetic_rho.py (standing in for "SIESTA has already
# been run there") before running stb-hirshfeldIons for real. Writes BOTH
# a cation and an anion reference folder per species, unconditionally --
# no more per-species sign decision (see hirshfeld_ions.py's own module
# docstring for why: individual atoms of the same species can converge
# toward opposite ion states, per the literal Bultinck et al. formulation).
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREP_DIR="$(cd "$FIXTURE_DIR/../prep" && pwd)"
GEN_SCRIPT="$FIXTURE_DIR/../make_synthetic_rho.py"
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

# Runs stb-hirshfeldPrep on the shared fixture (which itself writes
# combined/ + neutral/<species>/, nothing pre-computed needed), then
# fabricates a synthetic .RHO in combined/ and every neutral/<species>/ it
# wrote -- purely for PLUMBING coverage (folder writing, manifest
# chaining), not a physically meaningful charge.
make_prep_and_neutral_rho() {
    rm -rf hirshfeld_study
    stb-hirshfeldPrep -s structure.fdf --calc calc.fdf -p . -O hirshfeld_study --no-intro \
        > /dev/null 2>&1
    python3 "$GEN_SCRIPT" hirshfeld_study/combined/structure.fdf \
        hirshfeld_study/combined/siesta.RHO 24 "0.8,1.4" "1.4,1.0" > /dev/null 2>&1
    python3 "$GEN_SCRIPT" hirshfeld_study/neutral/C/structure.fdf \
        hirshfeld_study/neutral/C/siesta.RHO 18 "0.8" "1.4" > /dev/null 2>&1
    python3 "$GEN_SCRIPT" hirshfeld_study/neutral/O/structure.fdf \
        hirshfeld_study/neutral/O/siesta.RHO 18 "1.4" "1.0" > /dev/null 2>&1
}


# --- 1. Preparation ---
echo "--- Starting tester for STB-HIRSHFELDIONS (item 4.20.2) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$PREP_DIR/structure.fdf" "$PREP_DIR/calc.fdf" "$PREP_DIR/C.psf" "$PREP_DIR/O.psf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Missing Stage-1 output guard ---
echo -e "\n--- Testing the missing-Stage-1-output guard ---"
rm -rf hirshfeld_study
stb-hirshfeldIons -O hirshfeld_study -p . --no-intro > log_no_stage1.txt 2>&1
check_exit_code $? 1
check_contains "run stb-hirshfeldPrep first" log_no_stage1.txt


# --- 3. Default run: BOTH cation and anion written for every species ---
echo -e "\n--- Testing default run (cation + anion folders, every species) ---"
make_prep_and_neutral_rho
stb-hirshfeldIons -O hirshfeld_study -p . --no-intro > log_default.txt 2>&1
check_exit_code $? 0
check_success hirshfeld_study/hirshfeld_ions_manifest.json
check_success hirshfeld_study/ions/C/cation/structure.fdf
check_success hirshfeld_study/ions/C/cation/calc.fdf
check_success hirshfeld_study/ions/C/cation/config_extra.fdf
check_success hirshfeld_study/ions/C/anion/config_extra.fdf
check_success hirshfeld_study/ions/O/cation/config_extra.fdf
check_success hirshfeld_study/ions/O/anion/config_extra.fdf

echo "Testing: Stage 1 recap, native-charges cross-check, and pass-0 preview"
echo "         all land in the console report -- Stage 2's report is self-contained"
check_contains "STAGE 1 RECAP" log_default.txt
check_contains "Combined structure :" log_default.txt
check_contains "Species written by Stage 1:" log_default.txt
check_contains "NATIVE CHARGES CROSS-CHECK" log_default.txt
check_contains "PASS-0 (SIMPLE HIRSHFELD) CHARGES -- PREVIEW ONLY" log_default.txt
check_contains "N+ / N- = how many atoms" log_default.txt
check_contains "WRITING ION REFERENCES (cation + anion, every species)" log_default.txt

echo "Testing: the species table distinguishes Z (atomic number, for the pseudopotential)"
echo "         from Z_val (valence charge, for the charge formula) -- C: Z=6, Z_val=4"
check_contains "C       | 6 | 4     | hardcoded fallback" log_default.txt
check_contains "Z_val | Population(e-) | Charge(e-)" log_default.txt

echo "Testing: each ion folder carries a NetCharge directive (+1.0/-1.0 respectively) and"
echo "         is otherwise single-point/Gamma-only/spin-polarized like Stage 1 --"
echo "         Spin polarized forced BOTH via config_extra.fdf AND directly in calc.fdf"
check_contains "NetCharge            +1.0" hirshfeld_study/ions/C/cation/config_extra.fdf
check_contains "NetCharge            -1.0" hirshfeld_study/ions/C/anion/config_extra.fdf
check_contains "NetCharge            +1.0" hirshfeld_study/ions/O/cation/config_extra.fdf
check_contains "NetCharge            -1.0" hirshfeld_study/ions/O/anion/config_extra.fdf
check_contains "MD.Steps              0" hirshfeld_study/ions/C/cation/config_extra.fdf
check_contains "kgrid.MonkhorstPack   \[1  1  1\]" hirshfeld_study/ions/O/anion/config_extra.fdf
check_contains "%include config_extra.fdf" hirshfeld_study/ions/C/cation/calc.fdf
check_contains "Spin                polarized" hirshfeld_study/ions/C/cation/calc.fdf

echo "Testing: the ions manifest records both folders + pass0_charges per species (no more"
echo "         single net_charge/sign_source -- there's no single sign to decide anymore)"
check_contains '"symbol": "C"' hirshfeld_study/hirshfeld_ions_manifest.json
check_contains '"cation_folder": "ions/C/cation"' hirshfeld_study/hirshfeld_ions_manifest.json
check_contains '"anion_folder": "ions/C/anion"' hirshfeld_study/hirshfeld_ions_manifest.json
check_contains '"pass0_charges"' hirshfeld_study/hirshfeld_ions_manifest.json


# --- 4. Native-charges cross-check, once combined/ has a real .out ---
echo -e "\n--- Testing the native-charges cross-check (fabricated combined/siesta.out) ---"
cat > hirshfeld_study/combined/siesta.out << 'RHOEOF'
Hirshfeld Atomic Populations:
Atom #   charge [q] valence [e]      Sz [e]  Species
     1    -0.049444    4.049444   -0.005202  C
     2     0.049444    7.950556   -0.005202  O
-------------------------------------------
 Total     0.000000                0.000000

Voronoi Atomic Populations:
Atom #   charge [q] valence [e]      Sz [e]  Species
     1    -0.047137    4.047137   -0.011423  C
     2     0.047137    7.952863   -0.011423  O
-------------------------------------------
 Total     0.000000                0.000000
RHOEOF
stb-hirshfeldIons -O hirshfeld_study -p . --no-intro > log_nativecheck.txt 2>&1
check_exit_code $? 0
check_contains "Hirshfeld(e-) | Voronoi(e-)" log_nativecheck.txt
check_contains "\-0.0494" log_nativecheck.txt
check_contains "Native Hirshfeld(e-)" log_nativecheck.txt
rm -f hirshfeld_study/combined/siesta.out


# --- 5. Missing neutral .RHO (combined/ ran, neutral/<species>/ didn't) guard ---
echo -e "\n--- Testing the missing-neutral-.RHO guard ---"
rm -rf hirshfeld_study
stb-hirshfeldPrep -s structure.fdf --calc calc.fdf -p . -O hirshfeld_study --no-intro \
    > /dev/null 2>&1
python3 "$GEN_SCRIPT" hirshfeld_study/combined/structure.fdf \
    hirshfeld_study/combined/siesta.RHO 24 "0.8,1.4" "1.4,1.0" > /dev/null 2>&1
stb-hirshfeldIons -O hirshfeld_study -p . --no-intro > log_norho.txt 2>&1
check_exit_code $? 1
check_contains "has SIESTA been run there yet" log_norho.txt


# --- 6. Genuinely mixed-sign species -- BOTH ion folders written regardless ---
echo -e "\n--- Testing a genuinely mixed-sign species (structure_mixed.fdf) ---"
rm -rf hirshfeld_study_mixed mixed_combined
mkdir -p mixed_combined
cp "$FIXTURE_DIR/structure_mixed.fdf" mixed_combined/structure.fdf
cp calc.fdf mixed_combined/
stb-hirshfeldPrep -s mixed_combined/structure.fdf --calc mixed_combined/calc.fdf -p . \
    -O hirshfeld_study_mixed --mesh-cutoff 400 --vacuum 18 --no-intro > /dev/null 2>&1
python3 "$GEN_SCRIPT" hirshfeld_study_mixed/combined/structure.fdf \
    hirshfeld_study_mixed/combined/siesta.RHO 40 "0.02,3.0,1.0" "0.5,0.6,1.0" > /dev/null 2>&1
python3 "$GEN_SCRIPT" hirshfeld_study_mixed/neutral/C/structure.fdf \
    hirshfeld_study_mixed/neutral/C/siesta.RHO 24 "0.5" "0.8" > /dev/null 2>&1
python3 "$GEN_SCRIPT" hirshfeld_study_mixed/neutral/O/structure.fdf \
    hirshfeld_study_mixed/neutral/O/siesta.RHO 24 "1.0" "0.8" > /dev/null 2>&1
stb-hirshfeldIons -O hirshfeld_study_mixed -p . --no-intro > log_mixed.txt 2>&1
check_exit_code $? 0
echo "Testing: C reads genuinely mixed sign at pass-0 (1 cation-like, 1 anion-like) and the"
echo "         table is highlighted -- but BOTH ions/C/cation and ions/C/anion still get"
echo "         written (no fallback/vote needed anymore)"
check_contains "C    | 2 | 1  | 1  " log_mixed.txt
check_success hirshfeld_study_mixed/ions/C/cation/config_extra.fdf
check_success hirshfeld_study_mixed/ions/C/anion/config_extra.fdf
check_contains "NetCharge            +1.0" hirshfeld_study_mixed/ions/C/cation/config_extra.fdf
check_contains "NetCharge            -1.0" hirshfeld_study_mixed/ions/C/anion/config_extra.fdf


# --- 7. --ref (explicit reference file, e.g. not named with a .out extension) ---
echo -e "\n--- Testing --ref (Z_val detection from an explicitly-named, non-.out reference file) ---"
make_prep_and_neutral_rho
cat > hirshfeld_study/combined/run_log.txt << 'ZVALEOF'
atom: Called for C(Z=6)
Vna: chval, zval:    4.00000   4.00000
atom: Called for O(Z=8)
Vna: chval, zval:    6.00000   6.00000
ZVALEOF
stb-hirshfeldIons -O hirshfeld_study -p . --ref hirshfeld_study/combined/run_log.txt \
    --no-intro > log_ref.txt 2>&1
check_exit_code $? 0
check_contains "detected (SIESTA log)" log_ref.txt


# --- 8. Interactive path (stb-suite, shortcut 4.20.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.20.2) ---"
make_prep_and_neutral_rho
printf '4.20.2\nhirshfeld_study\n3\n.\n\nn\n\n0\n' | stb-suite > log_interactive.txt 2>&1
check_contains "Hirshfeld-I Ions (Stage 2) complete" log_interactive.txt
check_success hirshfeld_study/hirshfeld_ions_manifest.json


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
