#!/bin/bash

# Test for stb-chargediffAnalysis (Charge Density Difference, Stage 2 - Analysis, item 5.16.2).
# A finished stb-chargediffPrep run with a KNOWN answer is built (make_synthetic_run.py: Gaussian atoms, a planted
# 0.05 + 0.02 e transfer from the C atoms to the N atoms, fixed energies), then every output, mode and failure path of
# the analysis is checked against that answer.
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$FIXTURE_DIR/test_files"

GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'
PASS=0
FAIL=0

ok()   { echo -e "   -> ${GREEN}Verified:${NC} $1"; PASS=$((PASS+1)); }
bad()  { echo -e "   -> ${RED}Failed:${NC} $1"; FAIL=$((FAIL+1)); }

check_success() { if [ -s "$1" ]; then ok "file '$1' created"; else bad "file '$1' was not created"; fi; }
check_absent()  { if [ ! -e "$1" ]; then ok "file '$1' not written (opt-in output)"; else bad "file '$1' written without being asked"; fi; }
check_contains() { if grep -q -- "$1" "$2" 2>/dev/null; then ok "'$1' found in '$2'"; else bad "'$1' NOT found in '$2'"; fi; }
check_exit_code() { if [ "$1" -eq "$2" ]; then ok "exit code $1 (expected $2)"; else bad "exit code $1 (expected $2)"; fi; }
check_nonzero() { if [ "$1" -ne 0 ]; then ok "exit code $1 (a failure, as expected)"; else bad "exit code 0 (a failure was expected)"; fi; }

# check_json <results.json> <python expression on d> <expected> <tolerance> <description>
check_json() {
    local got
    got=$(python3 -c "import json,sys; d=json.load(open('$1')); print($2)" 2>/dev/null)
    if python3 -c "import sys; sys.exit(0 if abs(float('$got')-($3)) <= $4 else 1)" 2>/dev/null; then
        ok "$5 = $got (expected $3 +/- $4)"
    else
        bad "$5 = '$got' (expected $3 +/- $4)"
    fi
}

echo "--- Starting tester for stb-chargediffAnalysis (item 5.16.2) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/structure.fdf" "$TEST_DIR/structure.fdf"
cp "$FIXTURE_DIR/../prep/calc_template.fdf" "$TEST_DIR/calc.fdf"
pushd "$TEST_DIR" > /dev/null

# C and N atoms are interleaved in the combined structure (C N C N), the fragment folders list them grouped by label
stb-chargediffPrep --structure structure.fdf --fragment-species "C;N" --calc calc.fdf --no-intro \
    --output-dir run > log_prep.txt 2>&1
check_exit_code $? 0
python3 "$FIXTURE_DIR/make_synthetic_run.py" run > log_synth.txt 2>&1
check_exit_code $? 0
check_success run/combined/siesta.RHO
check_success run/fragment2/calc.out


# --- 1. Default run: slice, report, only the opt-out outputs ---
echo -e "\n--- Default run (slice) ---"
stb-chargediffAnalysis --dir run -o out_default --no-intro > log_default.txt 2>&1
check_exit_code $? 0
for s in "\[0\] RUN METADATA" "\[1\] SYSTEM" "\[2\] CONSISTENCY CHECKS" "\[3\] COUNTERPOISE" "\[4\] DELTA CHARGE DENSITY" \
         "\[5\] CHARGE CHANGE BY FRAGMENT" "\[6\] MAP AND FILES" "\[7\] REFERENCES" "\[8\] SUMMARY" "\[9\] WARNINGS"; do
    check_contains "$s" log_default.txt
done
check_contains "checks 8 ok, 0 warn" log_default.txt
check_contains "E_int (CP)     : -0.1500 eV" log_default.txt
check_success out_default/run_chargediff_slice.dat
check_success out_default/run_chargediff.cube
check_success out_default/references.bib
echo "Testing: plots and data tables are opt-in"
check_absent out_default/run_chargediff_summary.png
check_absent out_default/run_chargediff_slice.gplot
check_absent out_default/run_chargediff_results.json
check_absent out_default/run_chargediff_atoms.csv


# --- 2. The numbers against the planted answer ---
echo -e "\n--- --png --save-data: results against the known answer ---"
stb-chargediffAnalysis --dir run -o out_full --png --save-data --no-intro > log_full.txt 2>&1
check_exit_code $? 0
check_success out_full/run_chargediff_summary.png
check_success out_full/run_chargediff_atoms.csv
R=out_full/run_chargediff_results.json
check_json $R "d['integral_e']" 0 1e-3 "integral of Delta rho (e)"
check_json $R "d['charge_moved_e']" 0.07 2e-3 "charge moved (e)"
check_json $R "d['e_int_cp_ev']" -0.150 1e-6 "counterpoise E_int (eV)"
check_json $R "d['per_fragment_e']['fragment1']" -0.07 2e-3 "C fragment (e)"
check_json $R "d['per_fragment_e']['fragment2']" 0.07 2e-3 "N fragment (e)"
echo "Testing: per-atom charges follow the atoms although the folders list them in another order"
check_contains "1,C,C_frag1,fragment1,-0.04" out_full/run_chargediff_atoms.csv
check_contains "2,N,N_frag2,fragment2,0.04" out_full/run_chargediff_atoms.csv
check_contains "3,C,C_frag1,fragment1,-0.01" out_full/run_chargediff_atoms.csv
check_contains "4,N,N_frag2,fragment2,0.01" out_full/run_chargediff_atoms.csv


# --- 3. The cube is in e/bohr^3 (its positive part, over its own voxel volume in bohr^3, is the charge moved) ---
echo -e "\n--- Cube units ---"
python3 - out_full/run_chargediff.cube > log_cube.txt 2>&1 <<'PY'
import sys
import numpy as np
L = open(sys.argv[1]).read().split("\n")
na = int(L[2].split()[0])
n = [int(L[i].split()[0]) for i in (3, 4, 5)]
vec = np.array([[float(x) for x in L[i].split()[1:]] for i in (3, 4, 5)])      # bohr per voxel
d = np.array(" ".join(L[6 + na:]).split(), float)
dv = abs(np.linalg.det(vec))
pos = d[d > 0].sum() * dv
print(f"cube positive part: {pos:.5f} e")
sys.exit(0 if abs(pos - 0.07) < 2e-3 else 1)
PY
check_exit_code $? 0
cat log_cube.txt


# --- 4. Modes and options ---
echo -e "\n--- --profile ---"
stb-chargediffAnalysis --dir run -o out_profile --profile --save-gnuplot --no-intro > log_profile.txt 2>&1
check_exit_code $? 0
check_success out_profile/run_chargediff_profile.dat
check_success out_profile/run_chargediff_profile.gplot

echo -e "\n--- --3d ---"
stb-chargediffAnalysis --dir run -o out_3d --3d --iso-min 0.005 --save-gnuplot --no-intro > log_3d.txt 2>&1
check_exit_code $? 0
check_success out_3d/run_chargediff_3d.dat
check_success out_3d/run_chargediff_3d.gplot

echo -e "\n--- --axis/--pos/--contour with --no-cube ---"
stb-chargediffAnalysis --dir run -o out_axis --axis 0 --pos 3.0 --contour --no-cube --save-gnuplot --no-intro > log_axis.txt 2>&1
check_exit_code $? 0
check_contains "plane normal to x" log_axis.txt
check_success out_axis/run_chargediff_slice.gplot
check_absent out_axis/run_chargediff.cube

echo -e "\n--- --save-report ---"
stb-chargediffAnalysis --dir run -o out_report --save-report --no-intro > log_report.txt 2>&1
check_exit_code $? 0
check_success out_report/chargediff_analysis_report.txt
check_contains "COUNTERPOISE" out_report/chargediff_analysis_report.txt

echo -e "\n--- --show without a display falls back to saving ---"
env -u DISPLAY -u WAYLAND_DISPLAY stb-chargediffAnalysis --dir run -o out_show --show --no-intro > log_show.txt 2>&1
check_exit_code $? 0
check_contains "needs a display" log_show.txt
check_success out_show/run_chargediff_summary.png


# --- 5. Warnings and failures ---
echo -e "\n--- SCF not confirmed in one fragment: a warning, not a failure ---"
cp -r run run_scf
sed -i '/SCF cycle converged/d' run_scf/fragment1/calc.out
stb-chargediffAnalysis --dir run_scf -o out_scf --no-intro > log_scf.txt 2>&1
check_exit_code $? 0
check_contains "not confirmed: fragment1" log_scf.txt
check_contains "1 warning(s); checks 7 ok, 1 warn" log_scf.txt

echo -e "\n--- Different mesh actually used: warning ---"
cp -r run run_mesh
sed -i 's/= *400.000 *400.000/= 400.000 450.000/' run_mesh/fragment2/calc.out
stb-chargediffAnalysis --dir run_mesh -o out_mesh --no-intro > log_mesh.txt 2>&1
check_exit_code $? 0
check_contains "Same mesh actually used by SIESTA" log_mesh.txt
check_contains "fragment2: 450.0" log_mesh.txt

echo -e "\n--- Different level of theory between folders: warning ---"
cp -r run run_theory
echo "XC.functional GGA" >> run_theory/fragment2/config_extra.fdf
echo "XC.authors BLYP" >> run_theory/fragment2/config_extra.fdf
stb-chargediffAnalysis --dir run_theory -o out_theory --no-intro > log_theory.txt 2>&1
check_exit_code $? 0
check_contains "WARN" log_theory.txt

echo -e "\n--- Grid of a fragment on another mesh: error ---"
cp -r run run_grid
python3 - > /dev/null 2>&1 <<'PY'
import sisl
g = sisl.get_sile('run_grid/fragment2/siesta.RHO').read_grid()
g2 = sisl.Grid((16, 16, 36), geometry=g.geometry)
g2.grid[...] = 1
sisl.get_sile('run_grid/fragment2/siesta.RHO', mode='w').write_grid(g2)
PY
stb-chargediffAnalysis --dir run_grid -o out_grid --no-intro > log_grid.txt 2>&1
check_nonzero $?
check_contains "Grid shape mismatch" log_grid.txt

echo -e "\n--- Folder without a prep manifest: error ---"
stb-chargediffAnalysis --dir run_missing -o out_missing --no-intro > log_missing.txt 2>&1
check_nonzero $?
check_contains "chargediff_manifest.json" log_missing.txt

echo -e "\n--- --spin on a run that is not spin-polarized: error ---"
stb-chargediffAnalysis --dir run --spin -o out_spin --no-intro > log_spin.txt 2>&1
check_nonzero $?
check_contains "no spin component" log_spin.txt

popd > /dev/null
echo -e "\n--- Tests Complete ---"
echo -e "Passed: ${GREEN}${PASS}${NC}   Failed: ${RED}${FAIL}${NC}"
echo "Test files were kept in '$TEST_DIR' for inspection."
[ "$FAIL" -eq 0 ]
