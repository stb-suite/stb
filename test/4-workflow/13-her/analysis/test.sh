#!/bin/bash

# --- Setup ---
# Smoke test for stb-herAnalysis (HER Stage 3: Analysis, item 4.13.3).
# Builds a minimal her_study/ tree BY HAND (fabricated FreeEng/.FA data
# with known, hand-computed closed-form answers) rather than chaining a
# real stb-her/stb-herRefs run for every scenario -- Stage 1/2's own CLI
# wiring (including the real relaxed-geometry/ghost-atom/H2 machinery) is
# already covered end to end in test/.../13-her/{prep,refs}/test.sh. Same
# "only Stage 3's reading side is under test here" discipline as
# test/.../11-raman/analysis/test.sh's own "4c" reconstruction test.
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

# Writes the 6 common reference-energy calc.out files shared by every
# scenario below (only the ZPE-mode-specific pieces differ).
write_common_energies() {
    python3 - <<'PYEOF'
import os

def write_calc_out(path, freeeng):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write("siesta: SCF Convergence by DM criterion\n")
        f.write("SCF cycle converged after 12 iterations\n")
        f.write(f"siesta: FreeEng =        {freeeng:.6f}\n")
        f.write("siesta: Atomic forces (eV/Ang):\n")
        f.write("     Max    0.001\n")

os.makedirs("her_study/sites/site_1_ontop", exist_ok=True)
write_calc_out("her_study/sites/site_1_ontop/calc.out", -499.950000)
write_calc_out("her_study/00_clean_slab/calc.out", -400.000000)
write_calc_out("her_study/02_h2_molecule/calc.out", -31.500000)
write_calc_out("her_study/03_slab_deformed/calc.out", -399.900000)
write_calc_out("her_study/04_slab_ghost/calc.out", -399.950000)
write_calc_out("her_study/06_h_ghost_slab/calc.out", -100.050000)
write_calc_out("her_study/07_h_isolated/calc.out", -100.000000)

with open("her_study/her_stage2.txt", "w") as f:
    f.write("Winning site   : sites/site_1_ontop\n")
PYEOF
}


# --- 1. Preparation ---
echo "--- Starting tester for STB-HERAnalysis stage 3: analysis (item 4.13.3) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Missing Stage-2 output guard ---
echo -e "\n--- Testing the missing-Stage-2-output guard ---"
rm -rf her_study
stb-herAnalysis --directory her_study --no-intro > log_no_stage2.txt 2>&1
check_exit_code $? 1
check_contains "run stb-herRefs" log_no_stage2.txt


# --- 3. --zpe-mode standard ---
echo -e "\n--- Testing --zpe-mode standard (fixed Norskov offset) ---"
rm -rf her_study
write_common_energies
python3 -c "
with open('her_study/her_stage2.txt', 'a') as f:
    f.write('ZPE mode       : standard\n')
"
stb-herAnalysis --directory her_study --no-intro > log_standard.txt 2>&1
check_exit_code $? 0
check_contains "Standard Norskov offset" log_standard.txt
check_contains "+0.2400 eV" log_standard.txt

echo "Testing: BSSE correction and Delta-E_H match hand-computed values"
python3 -c "
import re
with open('her_study/her_stage3.txt') as f:
    text = f.read()
# E_clean=-400, E_slab_h=-499.95, E_h2=-31.5, E_deformed=-399.9,
# E_ghost=-399.95, E_h_ghost=-100.05, E_h_iso=-100.0
bsse_slab = -399.9 - (-399.95)      # = 0.05
bsse_ads = -100.0 - (-100.05)        # = 0.05
bsse_total = bsse_slab + bsse_ads    # = 0.10
delta_e_raw = -499.95 - (-400.0) - 0.5*(-31.5)  # = -84.2
delta_e_corrected = delta_e_raw + bsse_total     # = -84.1
delta_g_expected = delta_e_corrected + 0.24      # = -83.86

m = re.search(r'Delta-G_H\* = ([\-+\d.]+) eV', text)
assert m, 'could not find Delta-G_H* line'
got = float(m.group(1))
assert abs(got - delta_g_expected) < 1e-6, f'got {got} vs expected {delta_g_expected}'
print('OK')
" > log_standard_check.txt 2>&1
check_contains "OK" log_standard_check.txt


# --- 4. --zpe-mode local ---
echo -e "\n--- Testing --zpe-mode local (partial Hessian, known spring constant) ---"
rm -rf her_study
write_common_energies
python3 -c "
with open('her_study/her_stage2.txt', 'a') as f:
    f.write('ZPE mode       : local\n')
"
python3 - <<'PYEOF'
import os, json

zpe_dir = "her_study/05_zpe_calc"
os.makedirs(zpe_dir, exist_ok=True)

order = [{"axis": 0, "sign": 1.0}, {"axis": 0, "sign": -1.0},
         {"axis": 1, "sign": 1.0}, {"axis": 1, "sign": -1.0},
         {"axis": 2, "sign": 1.0}, {"axis": 2, "sign": -1.0}]
d = 0.015
with open(os.path.join(zpe_dir, "zpe_local_meta.json"), "w") as f:
    json.dump({"h_index": 2, "displacement_ang": d, "order": order}, f)

# Isotropic harmonic H* with a KNOWN spring constant k (eV/Ang^2).
k = 5.0
for i, entry in enumerate(order, start=1):
    axis, sign = entry["axis"], entry["sign"]
    disp = sign * d
    force = -k * disp
    fxyz = [0.0, 0.0, 0.0]
    fxyz[axis] = force
    disp_dir = os.path.join(zpe_dir, f"disp_{i:03d}")
    os.makedirs(disp_dir, exist_ok=True)
    with open(os.path.join(disp_dir, "her_zpe.FA"), "w") as f:
        f.write("3\n")
        f.write("1   0.000000000E+00   0.000000000E+00   0.000000000E+00\n")
        f.write("2   0.000000000E+00   0.000000000E+00   0.000000000E+00\n")
        f.write(f"3   {fxyz[0]: .9E}   {fxyz[1]: .9E}   {fxyz[2]: .9E}\n")
print("fixture ready")
PYEOF
stb-herAnalysis --directory her_study --temp 298.15 --no-intro > log_local.txt 2>&1
check_exit_code $? 0
check_contains "Local (partial-Hessian) H\* vibrational modes" log_local.txt

echo "Testing: ZPE(H*, local) matches the analytic isotropic-oscillator result"
python3 -c "
import re
import numpy as np

k = 5.0
H_MASS_AMU = 1.00794
FREQ_CONV = 15.633302
EV_PER_THZ = 0.00413566733
d_eig = k / H_MASS_AMU
freq_thz = FREQ_CONV * np.sqrt(d_eig)
e_mode = freq_thz * EV_PER_THZ
zpe_expected = 3 * 0.5 * e_mode  # 3 degenerate modes

with open('her_study/her_stage3.txt') as f:
    text = f.read()
m = re.search(r'ZPE\(H\*, local\)\s*=\s*([\-\d.]+) eV', text)
assert m, 'could not find ZPE(H*, local) line'
got = float(m.group(1))
assert abs(got - zpe_expected) < 1e-3, f'got {got} vs expected {zpe_expected}'
print('OK')
" > log_local_check.txt 2>&1
check_contains "OK" log_local_check.txt


# --- 5. Error cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: --version"
stb-herAnalysis --version > log_version.txt 2>&1
check_contains "stb-herAnalysis" log_version.txt

echo "Testing: --help documents --temp/--force-tolerance"
stb-herAnalysis --help > log_help.txt 2>&1
check_contains "temp" log_help.txt
check_contains "force-tolerance" log_help.txt
check_contains "clean slab's own" log_help.txt


# --- 6. Interactive path (stb-suite, shortcut 4.13.3) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.13.3) ---"

echo "Testing: navigate 4.13.3 -> defaults -> quit"
rm -rf her_study
write_common_energies
python3 -c "
with open('her_study/her_stage2.txt', 'a') as f:
    f.write('ZPE mode       : standard\n')
"
{
  echo "4.13.3"
  echo ""               # run_dir (default her_study)
  echo ""               # output_filename (default calc.out)
  echo ""               # temperature (default 298.15)
  echo ""               # force_tolerance (default 0.05)
  echo ""               # press enter to continue
  echo "0"              # quit stage submenu
} | stb-suite > log_menu.txt 2>&1
check_contains "Delta-G_H\*" log_menu.txt
check_success her_study/HER_report.txt


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
