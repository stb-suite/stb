#!/bin/bash

# --- Setup ---
# Smoke test for stb-eosAnalysis (Equation of State Analysis, item 5.15.2)
#
# Since no real SIESTA run is available here, this test generates a real
# 'vol_*' sweep via stb-eosInputs (so the scanned volumes are genuine, not
# hand-picked) and then injects SYNTHETIC calc.out files whose energies AND
# stress tensors follow an EXACT, known Birch-Murnaghan curve (E0=-300 eV,
# B0=90 GPa, B0'=4.3, V0=160.103 Ang^3, matching si_cubic.fdf's real 5.43 Ang
# cubic-Si volume): the energy from ase.eos.birchmurnaghan, the pressure from
# the ANALYTIC Birch-Murnaghan P(V) formula (independent of the numerical
# derivative stb-eosAnalysis uses). The fit, the pressure cross-check and the
# stress-derived bulk modulus can then all be checked against an exact ground
# truth, not just "did it run without crashing". (The real-SIESTA case is in
# examples/5.15-eos/.)
export MPLBACKEND=Agg
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

check_absent() {
    if [ -e "$1" ]; then
        echo -e " ... ${RED}FAIL${NC} (file '$1' exists but should not)"
        FAIL=$((FAIL+1))
    else
        echo -e " ... ${GREEN}OK${NC} (file '$1' correctly absent)"
        PASS=$((PASS+1))
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

# check_python "description" "code exiting 0 on success"
check_python() {
    python3 -c "$2"
    if [ $? -eq 0 ]; then
        echo -e "   -> ${GREEN}Verified:${NC} $1"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} $1"
        FAIL=$((FAIL+1))
    fi
}

echo "--- Starting tester for stb-eosAnalysis (item 5.15.2) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/../prep/si_cubic.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 1. Generate a real 'vol_*' sweep, then inject synthetic (known-EOS) calc.out ---
echo -e "\n--- Generating a real vol_* sweep and injecting a known synthetic EOS ---"
stb-eosInputs --file si_cubic.fdf --n-volumes 7 --strain-range 6.0 --no-intro -o eos_runs \
    > log_gen.txt 2>&1
check_exit_code $? 0

cat > inject.py << 'INJECT_EOF'
import sys
import numpy as np
import glob, os
from ase.eos import birchmurnaghan
from ase.units import kJ
from stb.core import structure_io

E0_true, B0_true_gpa, Bp_true, V0_true = -300.0, 90.0, 4.3, 160.1030
B0_true_evA3 = B0_true_gpa / (1.0 / kJ * 1.0e24)
root = sys.argv[1]


def pressure_evA3(v):
    x = V0_true / v          # analytic third-order Birch-Murnaghan pressure
    return 1.5 * B0_true_evA3 * (x ** (7 / 3) - x ** (5 / 3)) * (1 + 0.75 * (Bp_true - 4) * (x ** (2 / 3) - 1))


for folder in sorted(glob.glob(os.path.join(root, "vol_*"))):
    structure = structure_io.read_fdf(os.path.join(folder, "si_cubic.fdf"))
    lattice = structure.lattice
    volume = abs(np.linalg.det(lattice))
    energy = float(birchmurnaghan(volume, E0_true, B0_true_evA3, Bp_true, V0_true))
    sigma = -pressure_evA3(volume)          # SIESTA stress: positive = tensile, so sigma = -P
    with open(os.path.join(folder, "calc.out"), "w") as f:
        f.write("Version         : 5.4.2\n")
        f.write("outcell: Unit cell vectors (Ang):\n")
        for row in lattice:
            f.write(f"        {row[0]:.6f}    {row[1]:.6f}    {row[2]:.6f}\n")
        f.write("\n")
        f.write("siesta: SCF cycle converged after 12 iterations\n")
        f.write(f"siesta: FreeEng =      {energy:.6f}\n")
        f.write("siesta: Atomic forces (eV/Ang):\n")
        f.write("   Max    0.000100\n")
        f.write("siesta: Stress tensor (static) (eV/Ang**3):\n")
        for i in range(3):
            f.write("siesta: " + " ".join(f"{(sigma if i == j else 0.0):14.8f}" for j in range(3)) + "\n")
INJECT_EOF
python3 inject.py eos_runs
echo "Synthetic calc.out files written (known Birch-Murnaghan ground truth: "
echo "E0=-300 eV, B0=90 GPa, B0'=4.3, V0=160.103 Ang^3, analytic stress tensor)."


# --- 2. The report; the fit recovers the known ground truth exactly ---
echo -e "\n--- Testing the report and that the EOS fit recovers the known synthetic ground truth ---"
stb-eosAnalysis --dir eos_runs --no-intro --save-report > log_fit.txt 2>&1
check_exit_code $? 0
for sec in "\[0\] RUN METADATA" "\[1\] SCAN & LEVEL OF THEORY" "\[2\] READING FOLDERS" \
           "\[3\] VOLUME-ENERGY-PRESSURE TABLE" "\[4\] EQUATION-OF-STATE FIT" \
           "\[5\] PRESSURE CROSS-CHECK" "\[6\] TARGET PRESSURE" "\[7\] CHECKS" "\[8\] SUMMARY & FILES" \
           "\[9\] NEXT STEPS" "\[10\] WARNINGS"; do
    check_contains "$sec" log_fit.txt
done
check_success stb_eosAnalysis_report.txt
check_contains "STB-EOSANALYSIS REPORT" stb_eosAnalysis_report.txt
check_contains "Scan record     : eos_scan.json" log_fit.txt
check_contains "Uncertainty (leave-one-out)" log_fit.txt
check_contains "Window sensitivity" log_fit.txt
check_contains "None\." log_fit.txt

check_python "V0/E0/B0/B0'/R^2 all match the known synthetic ground truth" "
import re, sys
text = re.sub(r'\x1b\[[0-9;]*m', '', open('log_fit.txt').read())
v0 = float(re.search(r'Equilibrium volume\s+V0 : ([\-0-9.]+)', text).group(1))
e0 = float(re.search(r'Equilibrium energy\s+E0 : ([\-0-9.]+)', text).group(1))
b0 = float(re.search(r'Bulk modulus\s+B0 : ([\-0-9.]+)', text).group(1))
bp = float(re.search(r\"Pressure derivative\s+B0': ([\-0-9.]+)\", text).group(1))
r2 = float(re.search(r'Fit quality\s*: R\^2 ([0-9.]+)', text).group(1))
ok = (abs(v0 - 160.1030) < 0.01 and abs(e0 - (-300.0)) < 0.01
      and abs(b0 - 90.0) < 0.5 and abs(bp - 4.3) < 0.05 and r2 > 0.9999)
sys.exit(0 if ok else 1)
"

echo "Testing: the pressure cross-check against the ANALYTIC Birch-Murnaghan pressure"
check_python "P_SIESTA (analytic) agrees with -dE/dV of the fit, B from the stress ~ 90 GPa, P(V0) ~ 0" "
import re, sys
text = re.sub(r'\x1b\[[0-9;]*m', '', open('log_fit.txt').read())
dp = float(re.search(r'Max \|P_SIESTA - P_fit\|\s*: ([0-9.]+) GPa', text).group(1))
bs = float(re.search(r'B0 from the stress alone\s*: ([0-9.]+) GPa', text).group(1))
p0 = float(re.search(r'SIESTA pressure at V0\s*: ([\-+0-9.]+) GPa', text).group(1))
sys.exit(0 if dp < 0.02 and abs(bs - 90.0) < 2.0 and abs(p0) < 0.02 else 1)
"
echo "Testing: the checks table has no warning on a perfect synthetic curve"
python3 - << 'CHK_EOF'
import re, sys
text = re.sub(r'\x1b\[[0-9;]*m', '', open('log_fit.txt').read())
sec = text[text.index('[7] CHECKS'):text.index('[8] SUMMARY')]
print("   checks table rows:", len(re.findall(r'\| (OK|WARN)\s+\|', sec)), "| WARN:", len(re.findall(r'\| WARN\s+\|', sec)))
sys.exit(0 if len(re.findall(r'\| OK\s+\|', sec)) >= 12 and not re.findall(r'\| WARN\s+\|', sec) else 1)
CHK_EOF
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} every check in [7] is OK"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} a check of [7] is WARN (or the table is incomplete)"
    FAIL=$((FAIL+1))
fi
check_contains "Suggested LatticeConstant  : 5.43000 -> 5.43" log_fit.txt

echo "Testing: no plot or data file is written unless asked for"
check_absent eos_curve.dat
check_absent eos_curve.gplot
check_absent eos_curve_summary.png
check_absent eos_curve_results.json


# --- 3. Alternate EOS form (--eos vinet) still runs and fits well ---
echo -e "\n--- Testing --eos vinet ---"
stb-eosAnalysis --dir eos_runs --eos vinet --no-intro > log_vinet.txt 2>&1
check_exit_code $? 0
check_contains "EOS form        : vinet" log_vinet.txt
check_contains "Pressure derivative" log_vinet.txt


# --- 3b. --eos all fits every form on the same data and they all agree ---
echo -e "\n--- Testing --eos all ---"
stb-eosAnalysis --dir eos_runs --eos all --no-intro > log_all.txt 2>&1
check_exit_code $? 0
check_contains "EQUATION-OF-STATE FIT (ALL FORMS)" log_all.txt
check_contains "Spread between EOS forms" log_all.txt
check_contains "EOS forms agree on B0" log_all.txt
check_python "all 4 EOS forms agree on B0 ~= 90 GPa (known ground truth)" "
import re, sys
text = re.sub(r'\x1b\[[0-9;]*m', '', open('log_all.txt').read())
rows = re.findall(r'^(birchmurnaghan|vinet|murnaghan|sj)\s*\*?\s*\|\s*([\-0-9.]+)\s*\|\s*([\-0-9.]+)\s*\|\s*([\-0-9.]+)\s*\|\s*([\-0-9.]+)', text, re.M)
b0s = {r[0]: float(r[4]) for r in rows}
sys.exit(0 if len(b0s) == 4 and all(abs(b0 - 90.0) < 1.0 for b0 in b0s.values()) else 1)
"


# --- 3c. --target-pressure inverts the fit correctly ---
echo -e "\n--- Testing --target-pressure ---"
stb-eosAnalysis --dir eos_runs --target-pressure 0 5 -2 --no-intro > log_press.txt 2>&1
check_exit_code $? 0
check_contains "TARGET PRESSURE" log_press.txt
check_python "P=0 gives V0, P=5 GPa the analytic Birch-Murnaghan volume" "
import re, sys
import numpy as np
from scipy.optimize import brentq
text = re.sub(r'\x1b\[[0-9;]*m', '', open('log_press.txt').read())
rows = re.findall(r'^\s*([\-0-9.]+)\s*\|\s*([0-9.]+)\s*\|\s*[0-9.]+\s*\|\s*[0-9.]+\s*\|\s*within the scan', text, re.M)
got = {round(float(p), 3): float(v) for p, v in rows}
B0, Bp, V0 = 90.0 / 160.21766208, 4.3, 160.1030
def P(v):
    x = V0 / v
    return 1.5 * B0 * (x ** (7 / 3) - x ** (5 / 3)) * (1 + 0.75 * (Bp - 4) * (x ** (2 / 3) - 1)) * 160.21766208
ok = len(got) == 3
for p, v in got.items():
    ok = ok and abs(v - brentq(lambda w: P(w) - p, 100, 220)) < 0.05
sys.exit(0 if ok else 1)
"
echo "Testing: a target pressure outside the scanned range is flagged"
stb-eosAnalysis --dir eos_runs --target-pressure 50 --no-intro > log_extrap.txt 2>&1
check_exit_code $? 0
check_contains "EXTRAPOLATED" log_extrap.txt
check_contains "outside the range the scan achieved" log_extrap.txt
echo "Testing: --eos all with --target-pressure is rejected"
stb-eosAnalysis --dir eos_runs --eos all --target-pressure 5 --no-intro > log_mutex.txt 2>&1
check_exit_code $? 2
check_contains "needs one specific --eos form" log_mutex.txt


# --- 4. Outputs on request: gnuplot, PNG dashboard, data files ---
echo -e "\n--- Testing --save-gnuplot / --png / --save-data / --show ---"
stb-eosAnalysis --dir eos_runs --no-intro -o eos_out --save-gnuplot --png --save-data > log_out.txt 2>&1
check_exit_code $? 0
check_success eos_out.dat
check_success eos_out.gplot
check_success eos_out_summary.png
check_success eos_out_table.csv
check_success eos_out_results.json
check_contains "set multiplot layout 1,2" eos_out.gplot
check_contains "block 0: SIESTA points" eos_out.dat
check_python "the JSON has the fit, the checks and every point; the CSV one row per volume" "
import json, csv, sys
r = json.load(open('eos_out_results.json'))
rows = list(csv.DictReader(open('eos_out_table.csv')))
ok = (abs(r['fits']['birchmurnaghan']['b0_gpa'] - 90.0) < 0.5 and len(r['points']) == 7 and len(rows) == 7
      and all(c['ok'] for c in r['checks']) and abs(r['b0_from_stress_gpa'] - 90.0) < 2.0
      and r['jackknife']['b0_sigma'] < 0.5)
sys.exit(0 if ok else 1)
"
echo "Testing: --show without a display warns and still saves the figure"
env -u DISPLAY -u WAYLAND_DISPLAY stb-eosAnalysis --dir eos_runs --no-intro -o eos_show --show --png > log_show.txt 2>&1
check_exit_code $? 0
check_contains "needs a display" log_show.txt
check_success eos_show_summary.png


# --- 5. Warnings: planned vs actual, level of theory, unconverged SCF ---
echo -e "\n--- Testing the consistency warnings ---"
rm -rf eos_bad && cp -r eos_runs eos_bad
for d in eos_bad/vol_*; do
    printf 'Mesh.Cutoff 300 Ry\nPAO.BasisSize DZP\n' > "$d/calc.fdf"
done
printf 'Mesh.Cutoff 400 Ry\nPAO.BasisSize DZP\n' > eos_bad/vol_2.00/calc.fdf
sed -i '/SCF cycle converged/d' eos_bad/vol_m2.00/calc.out
python3 - << 'EDIT_EOF'
import re
p = "eos_bad/vol_0.00/calc.out"
s = open(p).read()
lines = s.split("\n")
out, n = [], 0
for line in lines:
    m = re.match(r"^(\s+)([\-0-9.]+)(\s+)([\-0-9.]+)(\s+)([\-0-9.]+)$", line)
    if m and n < 3:
        vals = [float(m.group(i)) * 1.02 for i in (2, 4, 6)]
        line = f"        {vals[0]:.6f}    {vals[1]:.6f}    {vals[2]:.6f}"
        n += 1
    out.append(line)
open(p, "w").write("\n".join(out))
EDIT_EOF
stb-eosAnalysis --dir eos_bad --no-intro > log_bad.txt 2>&1
check_exit_code $? 0
check_contains "level of theory differs between folders for meshcutoff" log_bad.txt
check_contains "SCF convergence not confirmed in vol_m2.00" log_bad.txt
check_contains "vol_0.00: SIESTA ran a cell of" log_bad.txt
check_contains "SCF converged in every run *| WARN" log_bad.txt
rm -rf eos_bad


# --- 6. A folder missing calc.out is skipped, not fatal ---
echo -e "\n--- Testing a folder missing calc.out is skipped (not fatal) ---"
mkdir -p eos_runs/vol_99.00
stb-eosAnalysis --dir eos_runs --no-intro > log_partial.txt 2>&1
check_exit_code $? 0
check_contains "vol_99.00.*SKIP" log_partial.txt
check_contains "1 skipped" log_partial.txt
check_contains "All planned runs available *| OK" log_partial.txt
rm -rf eos_runs/vol_99.00

echo "Testing: a scan that does not bracket the minimum warns"
rm -rf eos_edge && mkdir eos_edge
for d in vol_0.00 vol_2.00 vol_4.00 vol_6.00; do cp -r "eos_runs/$d" "eos_edge/$d"; done
stb-eosAnalysis --dir eos_edge --no-intro > log_edge.txt 2>&1
check_exit_code $? 0
check_contains "Points on both sides of V0 (>= 2 each) *| WARN" log_edge.txt
rm -rf eos_edge


# --- 7. Error and robustness cases ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: nonexistent --dir"
stb-eosAnalysis --dir nope_dir --no-intro > log_missing_dir.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_dir.txt

echo "Testing: --dir with no vol_* folders"
mkdir -p empty_dir
stb-eosAnalysis --dir empty_dir --no-intro > log_empty_dir.txt 2>&1
check_exit_code $? 1
check_contains "No 'vol_\*' folders found" log_empty_dir.txt

echo "Testing: fewer than 4 valid volumes (needs >= 4 for a 4-parameter fit)"
mkdir -p few_runs
for d in vol_0.00 vol_2.00 vol_4.00; do
    cp -r "eos_runs/$d" "few_runs/$d"
done
stb-eosAnalysis --dir few_runs --no-intro > log_toofew.txt 2>&1
check_exit_code $? 1
check_contains "at least 4" log_toofew.txt

echo "Testing: --version"
stb-eosAnalysis --version > log_version.txt 2>&1
check_contains "stb-eosAnalysis" log_version.txt

echo "Testing: --help documents the options"
stb-eosAnalysis --help > log_help.txt 2>&1
for opt in "\-\-dir" "\-\-file" "\-\-eos" "target-pressure" "\-\-output" "\-\-show" "\-\-png" "save-gnuplot" \
           "save-data" "save-report"; do
    check_contains "$opt" log_help.txt
done


# --- 8. Interactive path (stb-suite, shortcut 5.15.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 5.15.2) ---"
printf '5.15.2\neos_runs\ncalc.out\nbirchmurnaghan\neos_interactive\n\n\n\ny\n\n\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success eos_interactive.dat
check_success eos_interactive.gplot
check_absent eos_interactive_summary.png


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
