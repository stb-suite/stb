#!/bin/bash

# --- Setup ---
# Smoke test for stb-eosAnalysis (Equation of State Analysis, item 4.18.2)
#
# Since no real SIESTA run is available in this dev environment, this test
# generates a real 'vol_*' sweep via stb-eosInputs (so the scanned volumes
# are genuine, not hand-picked) and then injects SYNTHETIC calc.out files
# whose energies follow an EXACT, known Birch-Murnaghan curve (E0=-300 eV,
# B0=90 GPa, B0'=4.3, V0=160.103 Ang^3, matching si_cubic.fdf's real
# 5.43 Ang cubic-Si volume) -- built with ase.eos.birchmurnaghan itself, the
# same function stb-eosAnalysis's own fit uses internally. This lets the fit
# be checked against an exact ground truth (same rigor as stb-mleos's own
# synthetic-data verification), not just "did it run without crashing".
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


echo "--- Starting tester for stb-eosAnalysis (item 4.18.2) ---"
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

python3 << 'PYEOF'
import numpy as np
import glob, os
from ase.eos import birchmurnaghan
from ase.units import kJ
from stb.core import structure_io

E0_true, B0_true_gpa, Bp_true, V0_true = -300.0, 90.0, 4.3, 160.1030
B0_true_evA3 = B0_true_gpa / (1.0 / kJ * 1.0e24)

for folder in sorted(glob.glob("eos_runs/vol_*")):
    fdf = os.path.join(folder, "si_cubic.fdf")
    structure = structure_io.read_fdf(fdf)
    lattice = structure.lattice
    volume = abs(np.linalg.det(lattice))
    energy = float(birchmurnaghan(volume, E0_true, B0_true_evA3, Bp_true, V0_true))
    with open(os.path.join(folder, "calc.out"), "w") as f:
        f.write("outcell: Unit cell vectors (Ang):\n")
        for row in lattice:
            f.write(f"        {row[0]:.6f}    {row[1]:.6f}    {row[2]:.6f}\n")
        f.write("\n")
        f.write("siesta: SCF cycle converged after 12 iterations\n")
        f.write(f"siesta: FreeEng =      {energy:.6f}\n")
PYEOF
echo "Synthetic calc.out files written (known Birch-Murnaghan ground truth: "
echo "E0=-300 eV, B0=90 GPa, B0'=4.3, V0=160.103 Ang^3)."


# --- 2. Fit recovers the known ground truth exactly ---
echo -e "\n--- Testing the EOS fit recovers the known synthetic ground truth ---"
stb-eosAnalysis --dir eos_runs --no-intro -o eos_curve > log_fit.txt 2>&1
check_exit_code $? 0
check_contains "READING FOLDERS" log_fit.txt
check_contains "VOLUME-ENERGY TABLE" log_fit.txt
check_contains "EQUATION-OF-STATE FIT" log_fit.txt
check_success eos_curve.dat
check_success eos_curve.gplot
check_success stb_eosAnalysis_report.txt
check_contains "STB-EOSANALYSIS REPORT" stb_eosAnalysis_report.txt

python3 -c "
import re, sys
text = re.sub(r'\x1b\[[0-9;]*m', '', open('log_fit.txt').read())
v0 = float(re.search(r'Equilibrium volume \(V0\): ([\-0-9.]+)', text).group(1))
e0 = float(re.search(r'Equilibrium energy \(E0\): ([\-0-9.]+)', text).group(1))
b0 = float(re.search(r'Bulk modulus\s*\(B0\): ([\-0-9.]+)', text).group(1))
bp = float(re.search(r\"Pressure derivative \(B0'\): ([\-0-9.]+)\", text).group(1))
r2 = float(re.search(r'Fit quality\s*\(R\^2\): ([0-9.]+)', text).group(1))
ok = (abs(v0 - 160.1030) < 0.01 and abs(e0 - (-300.0)) < 0.01
      and abs(b0 - 90.0) < 0.5 and abs(bp - 4.3) < 0.05 and r2 > 0.9999)
sys.exit(0 if ok else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} V0/E0/B0/B0'/R^2 all match the known synthetic ground truth"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} fit did not recover the known synthetic ground truth"
    FAIL=$((FAIL+1))
fi


# --- 3. Alternate EOS form (--eos vinet) still runs and fits well ---
echo -e "\n--- Testing --eos vinet ---"
stb-eosAnalysis --dir eos_runs --eos vinet --no-intro -o eos_vinet > log_vinet.txt 2>&1
check_exit_code $? 0
check_contains "EOS form   : vinet" log_vinet.txt
check_contains "Pressure derivative" log_vinet.txt


# --- 3b. --eos all fits every form on the same data and they all agree ---
echo -e "\n--- Testing --eos all ---"
stb-eosAnalysis --dir eos_runs --eos all --no-intro -o eos_all > log_all.txt 2>&1
check_exit_code $? 0
check_contains "EQUATION-OF-STATE FIT (ALL FORMS)" log_all.txt
check_contains "^birchmurnaghan " log_all.txt
check_contains "^vinet " log_all.txt
check_contains "^murnaghan " log_all.txt
check_contains "^sj " log_all.txt

python3 -c "
import re, sys
text = re.sub(r'\x1b\[[0-9;]*m', '', open('log_all.txt').read())
rows = re.findall(r'^(birchmurnaghan|vinet|murnaghan|sj)\s+([\-0-9.]+)\s+([\-0-9.]+)\s+([\-0-9.]+)', text, re.M)
b0s = {name: float(b0) for name, v0, e0, b0 in rows}
sys.exit(0 if len(b0s) == 4 and all(abs(b0 - 90.0) < 1.0 for b0 in b0s.values()) else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} all 4 EOS forms agree on B0 ~= 90 GPa (known ground truth)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} EOS forms did not agree on the known B0"
    FAIL=$((FAIL+1))
fi


# --- 3c. --target-pressure inverts the fit correctly ---
echo -e "\n--- Testing --target-pressure ---"
stb-eosAnalysis --dir eos_runs --target-pressure 0 5 -2 --no-intro -o eos_press > log_press.txt 2>&1
check_exit_code $? 0
check_contains "TARGET PRESSURE" log_press.txt

python3 -c "
import re, sys
text = re.sub(r'\x1b\[[0-9;]*m', '', open('log_press.txt').read())
rows = re.findall(r'^\s*([\-0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s*\$', text, re.M)
by_p = {float(p): float(v) for p, v, scale in rows}
# P=0 should return essentially V0 (160.103); P=5 (compression) should give a
# smaller volume than P=-2 (tension) -- checked against the analytic
# Birch-Murnaghan pressure formula at the returned volumes, same B0/B0'/V0
# ground truth the synthetic calc.out files were built from.
def bm_pressure(V, V0=160.1030, B0=90.0, Bp=4.3):
    x = (V0 / V) ** (1.0 / 3.0)
    return 1.5 * B0 * (x**7 - x**5) * (1 + 0.75 * (Bp - 4) * (x**2 - 1))
ok = (abs(by_p[0.0] - 160.103) < 0.05
      and abs(bm_pressure(by_p[5.0]) - 5.0) < 0.01
      and abs(bm_pressure(by_p[-2.0]) - (-2.0)) < 0.01)
sys.exit(0 if ok else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} predicted volumes invert the known analytic BM pressure correctly"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} --target-pressure predictions did not match the analytic BM formula"
    FAIL=$((FAIL+1))
fi

echo "Testing: a wildly out-of-range pressure is flagged [EXTRAPOLATED]"
stb-eosAnalysis --dir eos_runs --target-pressure 50 --no-intro -o eos_extrap > log_extrap.txt 2>&1
check_exit_code $? 0
check_contains "EXTRAPOLATED" log_extrap.txt
check_contains "fall outside the range actually achieved" log_extrap.txt

echo "Testing: --eos all + --target-pressure together is rejected"
stb-eosAnalysis --dir eos_runs --eos all --target-pressure 5 --no-intro > log_mutex.txt 2>&1
check_exit_code $? 2
check_contains "needs one specific --eos form" log_mutex.txt


# --- 4. A folder missing calc.out is skipped, not fatal ---
echo -e "\n--- Testing a folder missing calc.out is skipped (not fatal) ---"
mkdir -p eos_runs/vol_99.00
stb-eosAnalysis --dir eos_runs --no-intro -o eos_partial > log_partial.txt 2>&1
check_exit_code $? 0
check_contains "vol_99.00.*SKIP" log_partial.txt
check_contains "skipped: 1" log_partial.txt
rm -rf eos_runs/vol_99.00


# --- 5. Error and robustness cases ---
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
check_contains "at least 4 are needed" log_toofew.txt

echo "Testing: --version"
stb-eosAnalysis --version > log_version.txt 2>&1
check_contains "stb-eosAnalysis" log_version.txt

echo "Testing: --help documents --dir, --file, --eos, --target-pressure, --output"
stb-eosAnalysis --help > log_help.txt 2>&1
check_contains "\-\-dir" log_help.txt
check_contains "\-\-file" log_help.txt
check_contains "eos" log_help.txt
check_contains "target-pressure" log_help.txt
check_contains "\-\-output" log_help.txt


# --- 6. Interactive path (stb-suite, shortcut 4.18.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.18.2) ---"
printf '4.18.2\neos_runs\ncalc.out\nbirchmurnaghan\neos_interactive\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success eos_interactive.dat


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
