#!/bin/bash

# --- Setup ---
# Smoke test for stb-eosInputs (Equation of State Prep, item 4.18.1)
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
echo "--- Starting tester for stb-eosInputs (item 4.18.1) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/si_cubic.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/structure.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/calc.fdf" "$TEST_DIR/calc_tmpl.fdf"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null

# a template that already has a k-grid, one that never %includes the structure, one with its own lattice
sed 's/^MD.VariableCell.*/MD.VariableCell true\nkgrid.MonkhorstPack [5 5 5]/' calc_tmpl.fdf > calc_kgrid.fdf
grep -v "^%include" calc_tmpl.fdf > calc_noinclude.fdf
sed 's/^%include si_cubic.fdf/%include si_cubic.fdf\nLatticeConstant 5.0 Ang/' calc_tmpl.fdf > calc_ownlattice.fdf
# a self-contained input (structure + settings in one file)
cat si_cubic.fdf > si_full.fdf
printf 'XC.Functional GGA\nXC.Authors PBE\nPAO.BasisSize DZP\nMesh.Cutoff 200 Ry\nkgrid.MonkhorstPack [4 4 4]\n' >> si_full.fdf
# a hexagonal crystal (2-atom hcp Mg) for the non-cubic note
cat > mg_hcp.fdf <<'EOF'
NumberOfSpecies 1
NumberOfAtoms 2
%block ChemicalSpeciesLabel
 1 12 Mg
%endblock ChemicalSpeciesLabel
LatticeConstant 1.0 Ang
AtomicCoordinatesFormat Fractional
%block LatticeVectors
 3.21 0.0 0.0
 -1.605 2.77994 0.0
 0.0 0.0 5.21
%endblock LatticeVectors
%block AtomicCoordinatesAndAtomicSpecies
 0.333333333 0.666666667 0.25 1
 0.666666667 0.333333333 0.75 1
%endblock AtomicCoordinatesAndAtomicSpecies
EOF


# --- 2. Basic scan, no calc.fdf in the folder: structure-only mode ---
echo -e "\n--- Testing a basic scan (bulk Si, no calc template: scaled structures only) ---"
stb-eosInputs --file si_cubic.fdf --no-intro -o eos_runs > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "Detected        : 3D" log_basic.txt
check_contains "Calc template   : (none" log_basic.txt
check_contains "\[3\] SCAN DEFINITION" log_basic.txt
check_success eos_runs/vol_0.00/si_cubic.fdf
check_success eos_runs/vol_5.00/si_cubic.fdf
check_success eos_runs/vol_m5.00/si_cubic.fdf
check_success eos_runs/eos_scan.json

echo "Testing: structure-only mode writes nothing but the structure"
if [ -e eos_runs/vol_0.00/calc.fdf ] || [ -e eos_runs/vol_0.00/config_extra.fdf ]; then
    echo -e "   -> ${RED}Failed:${NC} calc.fdf/config_extra.fdf written without a template"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} only the scaled structure is in the folder"
    PASS=$((PASS+1))
fi

echo "Testing: vol_0.00's raw lattice block is unchanged (zero strain, factor=1)"
check_contains "1.000000000000        0.000000000000        0.000000000000" eos_runs/vol_0.00/si_cubic.fdf

echo "Testing: vol_5.00's lattice is scaled by (1.05)^(1/3)"
check_contains "1.016396" eos_runs/vol_5.00/si_cubic.fdf

echo "Testing: atomic positions (fractional) are preserved verbatim"
check_contains "0.250000000   0.250000000   0.250000000   1" eos_runs/vol_5.00/si_cubic.fdf

echo "Testing: LatticeConstant line is untouched"
check_contains "LatticeConstant 5.43 Ang" eos_runs/vol_5.00/si_cubic.fdf

echo "Testing: the scan record (eos_scan.json) matches the plan"
python3 -c "
import json, sys
m = json.load(open('eos_runs/eos_scan.json'))
pts = m['scan']['points']
v0 = m['structure']['volume']
ok = (m['format'] == 'stb-eos-scan' and len(pts) == 9 and m['structure']['n_atoms'] == 8
      and m['structure']['species'] == ['Si'] and abs(v0 - 5.43 ** 3) < 1e-3
      and abs(pts[-1]['volume'] / v0 - 1.05) < 1e-9 and abs(pts[0]['strain_pct'] + 5.0) < 1e-12
      and m['calc']['mode'] == 'structure-only' and abs(m['structure']['latticeconstant'] - 5.43) < 1e-9)
sys.exit(0 if ok else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} 9 planned volumes from -5 % to +5 %, reference V and LatticeConstant recorded"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} eos_scan.json does not match the plan"
    FAIL=$((FAIL+1))
fi
rm -rf eos_runs


# --- 3. With a template and pseudopotentials: config_extra.fdf ---
echo -e "\n--- Testing a scan with a calc template (-c) and a pseudopotential bank (-p dojo) ---"
stb-eosInputs --file si_cubic.fdf -c calc_tmpl.fdf -p dojo --n-volumes 5 --strain-range 4.0 --no-intro \
    -o eos_tmpl --save-report > log_tmpl.txt 2>&1
check_exit_code $? 0
for n in 0 1 2 3 4 5 6; do
    check_contains "\[$n\] " log_tmpl.txt
done
check_success eos_tmpl/vol_0.00/calc.fdf
check_success eos_tmpl/vol_0.00/config_extra.fdf
check_success eos_tmpl/vol_0.00/si_cubic.fdf
check_success eos_tmpl/vol_0.00/Si.psml
check_success eos_tmpl/stb_eosInputs_report.txt
check_contains "STB-EOSINPUTS REPORT" eos_tmpl/stb_eosInputs_report.txt
check_contains "SCF k-grid forced : " log_tmpl.txt
check_contains "the template defines no k-grid" log_tmpl.txt

echo "Testing: calc.fdf = '%include config_extra.fdf' on top + the template, untouched"
head -1 eos_tmpl/vol_0.00/calc.fdf | grep -q "^%include config_extra.fdf$"
if [ $? -eq 0 ] && diff <(tail -n +3 eos_tmpl/vol_0.00/calc.fdf) calc_tmpl.fdf > /dev/null; then
    echo -e "   -> ${GREEN}Verified:${NC} the template is copied byte for byte below the %include"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} calc.fdf is not 'include + untouched template'"
    FAIL=$((FAIL+1))
fi

echo "Testing: config_extra.fdf forces a fixed-cell single point and the missing k-grid"
check_contains "^MD.VariableCell       false" eos_tmpl/vol_0.00/config_extra.fdf
check_contains "^MD.Steps              0" eos_tmpl/vol_0.00/config_extra.fdf
check_contains "^kgrid.MonkhorstPack" eos_tmpl/vol_0.00/config_extra.fdf
check_contains "template had MD.VariableCell true\|The template had MD.VariableCell true" log_tmpl.txt
for v in vol_m4.00 vol_m2.00 vol_0.00 vol_2.00 vol_4.00; do
    check_success eos_tmpl/$v/config_extra.fdf
done
rm -rf eos_tmpl


# --- 4. K-grid handling ---
echo -e "\n--- Testing k-grid handling ---"
stb-eosInputs --file si_cubic.fdf -c calc_kgrid.fdf --n-volumes 5 --no-intro -o eos_kg > log_kg.txt 2>&1
check_exit_code $? 0
check_contains "the template's own grid is kept" log_kg.txt
if grep -q "kgrid" eos_kg/vol_0.00/config_extra.fdf; then
    echo -e "   -> ${RED}Failed:${NC} a k-grid was forced although the template has one"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} the template's grid is kept (nothing forced)"
    PASS=$((PASS+1))
fi
stb-eosInputs --file si_cubic.fdf -c calc_kgrid.fdf --n-volumes 5 --kgrid 7 7 7 --no-intro -o eos_kg7 > log_kg7.txt 2>&1
check_exit_code $? 0
check_contains "kgrid.MonkhorstPack   \[7  7  7\]" eos_kg7/vol_0.00/config_extra.fdf
stb-eosInputs --file si_cubic.fdf -c calc_kgrid.fdf --n-volumes 5 --kgrid-density 0.1 --no-intro -o eos_kgd > log_kgd.txt 2>&1
check_exit_code $? 0
check_contains "explicit\|--kgrid-density 0.1" log_kgd.txt
check_contains "^kgrid.MonkhorstPack" eos_kgd/vol_0.00/config_extra.fdf
rm -rf eos_kg eos_kg7 eos_kgd


# --- 5. --relax-positions ---
echo -e "\n--- Testing --relax-positions / --relax-steps ---"
stb-eosInputs --file si_cubic.fdf -c calc_tmpl.fdf --n-volumes 5 --relax-positions --relax-steps 25 \
    --no-intro -o eos_relax > log_relax.txt 2>&1
check_exit_code $? 0
check_contains "^MD.Steps              25" eos_relax/vol_0.00/config_extra.fdf
check_contains "^MD.MaxForceTol        0.01 eV/Ang" eos_relax/vol_0.00/config_extra.fdf
check_contains "^MD.VariableCell       false" eos_relax/vol_0.00/config_extra.fdf
check_contains "relax at every volume" log_relax.txt
rm -rf eos_relax


# --- 6. Self-contained input (--file is also the template) ---
echo -e "\n--- Testing a self-contained input (-c equal to --file) ---"
stb-eosInputs --file si_full.fdf -c si_full.fdf --n-volumes 5 --no-intro -o eos_sc > log_sc.txt 2>&1
check_exit_code $? 0
check_contains "the structure file is itself the input" log_sc.txt
check_success eos_sc/vol_0.00/si_full.fdf
check_success eos_sc/vol_0.00/config_extra.fdf
head -1 eos_sc/vol_0.00/si_full.fdf | grep -q "^%include config_extra.fdf$"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} the scaled self-contained file starts with %include config_extra.fdf"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} the scaled self-contained file does not start with the %include"
    FAIL=$((FAIL+1))
fi
check_contains "1.016396" eos_sc/vol_5.00/si_full.fdf
rm -rf eos_sc


# --- 7. Template problems are warned about ---
echo -e "\n--- Testing template advisories ---"
stb-eosInputs --file si_cubic.fdf -c calc_noinclude.fdf --n-volumes 5 --no-intro -o eos_w1 > log_w1.txt 2>&1
check_exit_code $? 0
check_contains "does not %include 'si_cubic.fdf'" log_w1.txt
stb-eosInputs --file si_cubic.fdf -c calc_ownlattice.fdf --n-volumes 5 --no-intro -o eos_w2 > log_w2.txt 2>&1
check_exit_code $? 0
check_contains "defines a lattice of its own" log_w2.txt
stb-eosInputs --file si_cubic.fdf -c calc_tmpl.fdf --n-volumes 6 --no-intro -o eos_w3 > log_w3.txt 2>&1
check_exit_code $? 0
check_contains "even number of volumes" log_w3.txt
stb-eosInputs --file si_cubic.fdf -c calc_tmpl.fdf --n-volumes 5 --strain-range 1.0 --no-intro -o eos_w4 > log_w4.txt 2>&1
check_contains "very narrow" log_w4.txt
rm -rf eos_w1 eos_w2 eos_w3 eos_w4


# --- 8. Non-cubic crystal: the constrained-path note ---
echo -e "\n--- Testing the non-cubic note (hcp Mg) ---"
stb-eosInputs --file mg_hcp.fdf --n-volumes 5 --no-intro -o eos_hcp > log_hcp.txt 2>&1
check_exit_code $? 0
check_contains "hexagonal" log_hcp.txt
check_contains "keeps its shape" log_hcp.txt
python3 -c "
import json, sys
m = json.load(open('eos_hcp/eos_scan.json'))
sys.exit(0 if m['structure']['crystal_system'] == 'hexagonal' else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} crystal system recorded as hexagonal"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} crystal system not recorded"
    FAIL=$((FAIL+1))
fi
rm -rf eos_hcp


# --- 9. Custom --n-volumes/--strain-range ---
echo -e "\n--- Testing --n-volumes 5 --strain-range 10.0 ---"
stb-eosInputs --file si_cubic.fdf --n-volumes 5 --strain-range 10.0 --no-intro \
    -o eos_custom > log_custom.txt 2>&1
check_exit_code $? 0
check_success eos_custom/vol_10.00/si_cubic.fdf
check_success eos_custom/vol_5.00/si_cubic.fdf
check_success eos_custom/vol_0.00/si_cubic.fdf
check_success eos_custom/vol_m5.00/si_cubic.fdf
check_success eos_custom/vol_m10.00/si_cubic.fdf
if [ -d eos_custom/vol_2.00 ]; then
    echo -e "   -> ${RED}Failed:${NC} unexpected vol_2.00 folder for a 5-point (-10..10) scan"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} exactly 5 points scanned (no off-grid vol_2.00 folder)"
    PASS=$((PASS+1))
fi
rm -rf eos_custom


# --- 10. Physics: refuse a vacuum-padded structure ---
echo -e "\n--- Testing vacuum-padded (2D) structure rejection ---"
stb-eosInputs --file structure.fdf --no-intro -o eos_2d > log_2d.txt 2>&1
check_exit_code $? 1
check_contains "bulk (3D periodic) only" log_2d.txt


# --- 11. Argument validation ---
echo -e "\n--- Testing argument validation ---"
echo "Testing: --n-volumes 3 is rejected (needs >= 5)"
stb-eosInputs --file si_cubic.fdf --n-volumes 3 --no-intro -o eos_toofew > log_toofew.txt 2>&1
check_exit_code $? 2
echo "Testing: --strain-range 0 is rejected"
stb-eosInputs --file si_cubic.fdf --strain-range 0 --no-intro -o eos_zero > log_zero.txt 2>&1
check_exit_code $? 2
echo "Testing: missing --file"
stb-eosInputs --no-intro > log_missing_args.txt 2>&1
check_exit_code $? 2
echo "Testing: nonexistent input file"
stb-eosInputs --file nope.fdf --no-intro > log_missing_input.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_input.txt
echo "Testing: an explicit -c that does not exist is an error (the default calc.fdf is not)"
stb-eosInputs --file si_cubic.fdf -c nope_calc.fdf --no-intro -o eos_nocalc > log_nocalc.txt 2>&1
check_exit_code $? 1
check_contains "Calc file 'nope_calc.fdf' not found" log_nocalc.txt
echo "Testing: an unknown pseudopotential source"
stb-eosInputs --file si_cubic.fdf -c calc_tmpl.fdf -p no_such_bank_or_dir --no-intro -o eos_nopp > log_nopp.txt 2>&1
check_exit_code $? 1
check_contains "not a recognized pseudopotential bank" log_nopp.txt

echo "Testing: --version"
stb-eosInputs --version > log_version.txt 2>&1
check_contains "stb-eosInputs" log_version.txt

echo "Testing: --help documents the options"
stb-eosInputs --help > log_help.txt 2>&1
for opt in "\-\-file" "\-\-calc" "pseudo-dir" "n-volumes" "strain-range" "relax-positions" "relax-steps" \
           "\-\-kgrid" "kgrid-density" "output-dir" "vacuum-gap" "save-report"; do
    check_contains "$opt" log_help.txt
done


# --- 12. Interactive path (stb-suite, shortcut 4.18.1) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.18.1) ---"
rm -rf eos_interactive
printf '4.18.1\nsi_cubic.fdf\ncalc_tmpl.fdf\n\n5\n10.0\neos_interactive\n\n\n\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success eos_interactive/vol_10.00/si_cubic.fdf
check_success eos_interactive/vol_10.00/config_extra.fdf
check_success eos_interactive/vol_10.00/calc.fdf


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
