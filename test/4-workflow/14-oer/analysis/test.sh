#!/bin/bash

# --- Setup ---
# Smoke test for stb-oerAnalysis (OER Stage 4: Analysis, item 4.14.4).
# Builds a minimal oer_study/ tree BY HAND (fabricated FreeEng/.FA data
# with known, hand-computed closed-form answers) rather than chaining a
# real stb-oer/stb-oerIntermediates/stb-oerRefs run for every scenario --
# Stages 1-3's own CLI wiring is already covered end to end in
# test/.../14-oer/{prep,intermediates,refs}/test.sh. Same "only Stage 4's
# reading side is under test here" discipline as test/.../13-her/
# analysis/test.sh.
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

# Writes every FreeEng-bearing calc.out this Stage needs, all with known,
# round, hand-picked values, plus the Stage2/Stage3 report lines Stage 4
# reads back (strategy/BSSE mode/ZPE mode).
write_common_energies() {
    python3 - <<'PYEOF'
import os

def write_calc_out(path, freeeng):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write("siesta: SCF Convergence by DM criterion\n")
        f.write(f"siesta: FreeEng =        {freeeng:.6f}\n")
        f.write("siesta: Atomic forces (eV/Ang):\n")
        f.write("     Max    0.001\n")

os.makedirs("oer_study/sites/site_1_ontop", exist_ok=True)
os.makedirs("oer_study/intermediates/o_star", exist_ok=True)
os.makedirs("oer_study/intermediates/ooh_star", exist_ok=True)
write_calc_out("oer_study/sites/site_1_ontop/calc.out", -499.950000)     # E_OH
write_calc_out("oer_study/intermediates/o_star/calc.out", -480.000000)   # E_O
write_calc_out("oer_study/intermediates/ooh_star/calc.out", -520.000000) # E_OOH
write_calc_out("oer_study/00_clean_slab/calc.out", -400.000000)          # E_clean
write_calc_out("oer_study/02_h2_molecule/calc.out", -31.500000)          # E_H2
write_calc_out("oer_study/03_h2o_molecule/calc.out", -470.000000)        # E_H2O
write_calc_out("oer_study/04_slab_deformed/calc.out", -399.900000)       # diagnostic only

# Shared BSSE triad, small known round shifts (same style as HER's own
# bsse_slab=0.05/bsse_ads=0.05 fixture).
write_calc_out("oer_study/05_bsse_slab_only/calc.out", -520.050000)
write_calc_out("oer_study/05_bsse_slab_ghost/calc.out", -520.000000)
write_calc_out("oer_study/05_bsse_adsorbate_ghost_slab/calc.out", -120.050000)
write_calc_out("oer_study/05_bsse_isolated/calc.out", -120.000000)

with open("oer_study/oer_stage2.txt", "w") as f:
    f.write("O* strategy     : derived\n")
    f.write("OOH* strategy   : derived\n")

with open("oer_study/oer_stage3.txt", "w") as f:
    f.write("BSSE mode       : shared\n")
    f.write("ZPE mode        : local\n")
PYEOF
}

# Writes an isotropic, diagonal, decoupled local-mode Hessian fixture for
# one intermediate: k_O eV/Ang^2 on every O atom, k_H eV/Ang^2 on every H
# atom (k_H reuses HER's own k=5.0 value on purpose, so the two test
# suites cross-check each other's mass-weighting/conversion constants).
write_local_zpe_fixture() {
    python3 - <<PYEOF
import os, json

zpe_dir = "oer_study/08_zpe_calc_$1"
local_indices = [$2]
local_symbols = [$3]
system_label = "$4"
k_by_symbol = {"O": $5, "H": $6}

os.makedirs(zpe_dir, exist_ok=True)
order = []
for atom_index in local_indices:
    for axis in range(3):
        for sign in (1.0, -1.0):
            order.append({"atom_index": atom_index, "axis": axis, "sign": sign})

d = 0.015
with open(os.path.join(zpe_dir, "zpe_local_meta.json"), "w") as f:
    json.dump({"local_indices": local_indices, "local_symbols": local_symbols,
               "displacement_ang": d, "order": order, "system_label": system_label}, f)

n_atoms = max(local_indices) + 1
for i, entry in enumerate(order, start=1):
    atom_index, axis, sign = entry["atom_index"], entry["axis"], entry["sign"]
    moved_symbol = local_symbols[local_indices.index(atom_index)]
    k = k_by_symbol[moved_symbol]
    disp = sign * d
    force = -k * disp
    disp_dir = os.path.join(zpe_dir, f"disp_{i:03d}")
    os.makedirs(disp_dir, exist_ok=True)
    with open(os.path.join(disp_dir, f"{system_label}.FA"), "w") as f:
        f.write(f"{n_atoms}\n")
        for j in range(n_atoms):
            if j == atom_index:
                fxyz = [0.0, 0.0, 0.0]
                fxyz[axis] = force
            else:
                fxyz = [0.0, 0.0, 0.0]
            f.write(f"{j+1}   {fxyz[0]: .9E}   {fxyz[1]: .9E}   {fxyz[2]: .9E}\n")
PYEOF
}


# --- 1. Preparation ---
echo "--- Starting tester for STB-OERAnalysis stage 4: analysis (item 4.14.4) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Missing Stage-3 output guard ---
echo -e "\n--- Testing the missing-Stage-3-output guard ---"
rm -rf oer_study
stb-oerAnalysis --directory oer_study --no-intro > log_no_stage3.txt 2>&1
check_exit_code $? 1
check_contains "run stb-oerRefs" log_no_stage3.txt


# --- 3. Full default run, local ZPE mode ---
echo -e "\n--- Testing default generation (--zpe-mode local, isotropic known spring constants) ---"
rm -rf oer_study
write_common_energies
# OH* = C,C,O,H (local_indices 2,3 = O,H); O* = C,C,O (local_indices 2 = O);
# OOH* = C,C,O,O,H (local_indices 2,3,4 = O,O,H); H2O = O,H,H (local_indices 0,1,2 = O,H,H)
write_local_zpe_fixture "OH" "2, 3" "'O', 'H'" "oer_zpe_oh" 8.0 5.0
write_local_zpe_fixture "O" "2" "'O'" "oer_zpe_o" 8.0 5.0
write_local_zpe_fixture "OOH" "2, 3, 4" "'O', 'O', 'H'" "oer_zpe_ooh" 8.0 5.0
write_local_zpe_fixture "H2O" "0, 1, 2" "'O', 'H', 'H'" "oer_zpe_h2o" 8.0 5.0

stb-oerAnalysis --directory oer_study --temp 298.15 --no-intro > log_local.txt 2>&1
check_exit_code $? 0
check_contains "Local (partial-Hessian, 1-atom) vibrational modes" log_local.txt
check_contains "Local (partial-Hessian, 2-atom) vibrational modes" log_local.txt
check_contains "Local (partial-Hessian, 3-atom) vibrational modes" log_local.txt
check_success oer_study/oer_stage4.txt
check_success oer_study/OER_report.txt

echo "Testing: ZPE(O*, local) matches the analytic isotropic single-atom result (k_O=8.0 eV/Ang^2)"
python3 -c "
import re
import numpy as np

k_O = 8.0
O_MASS_AMU = 15.999
FREQ_CONV = 15.633302
EV_PER_THZ = 0.00413566733
freq_thz = FREQ_CONV * np.sqrt(k_O / O_MASS_AMU)
e_mode = freq_thz * EV_PER_THZ
zpe_expected = 3 * 0.5 * e_mode  # 3 degenerate modes, single isotropic O atom

with open('oer_study/oer_stage4.txt') as f:
    text = f.read()
m = re.search(r'ZPE\(O\*\) = ([\-\d.]+) eV', text)
assert m, 'could not find ZPE(O*) line'
got = float(m.group(1))
assert abs(got - zpe_expected) < 1e-3, f'got {got} vs expected {zpe_expected}'
print('OK')
" > log_o_zpe_check.txt 2>&1
check_contains "OK" log_o_zpe_check.txt

echo "Testing: ZPE(OH*, local) matches the analytic 2-atom (1 O + 1 H, decoupled/diagonal) result"
python3 -c "
import re
import numpy as np

k_O, k_H = 8.0, 5.0
O_MASS_AMU, H_MASS_AMU = 15.999, 1.00794
FREQ_CONV = 15.633302
EV_PER_THZ = 0.00413566733

def zpe_for(k, m, n_modes=3):
    freq_thz = FREQ_CONV * np.sqrt(k / m)
    e_mode = freq_thz * EV_PER_THZ
    return n_modes * 0.5 * e_mode

zpe_expected = zpe_for(k_O, O_MASS_AMU) + zpe_for(k_H, H_MASS_AMU)

with open('oer_study/oer_stage4.txt') as f:
    text = f.read()
m = re.search(r'ZPE\(OH\*\) = ([\-\d.]+) eV', text)
assert m, 'could not find ZPE(OH*) line'
got = float(m.group(1))
assert abs(got - zpe_expected) < 1e-3, f'got {got} vs expected {zpe_expected}'
print('OK')
" > log_oh_zpe_check.txt 2>&1
check_contains "OK" log_oh_zpe_check.txt

echo "Testing: G(O2)/4.92 eV identity holds exactly (dG1+dG2+dG3+dG4 == 4.92), independent of input energies"
python3 -c "
import re
with open('oer_study/oer_stage4.txt') as f:
    text = f.read()
m = re.search(r'Sum dG1\+dG2\+dG3\+dG4 = ([\-\d.]+) eV', text)
assert m, 'could not find the internal sanity-check sum line'
total = float(m.group(1))
assert abs(total - 4.920) < 1e-6, f'got {total}, expected 4.920 (identity by construction of G(O2))'
print('OK')
" > log_identity_check.txt 2>&1
check_contains "OK" log_identity_check.txt

echo "Testing: eta reported matches max(dG1..dG4) - 1.23, and PDS matches the argmax"
python3 -c "
import re
with open('oer_study/oer_stage4.txt') as f:
    text = f.read()
rows = re.findall(r'Step (\d):.*?dG = ([\-+\d.]+) eV', text)
dg_map = {int(k): float(v) for k, v in rows}
assert len(dg_map) == 4, f'expected 4 steps, got {dg_map}'
best_step = max(dg_map, key=dg_map.get)
eta_expected = dg_map[best_step] - 1.23

m_eta = re.search(r'overpotential eta = ([\-+\d.]+) V', text)
assert m_eta, 'could not find eta line'
eta_got = float(m_eta.group(1))
assert abs(eta_got - eta_expected) < 1e-6, f'eta: got {eta_got} vs expected {eta_expected}'

m_pds = re.search(r'PDS\) = Step (\d)', text)
assert m_pds, 'could not find PDS line'
pds_got = int(m_pds.group(1))
assert pds_got == best_step, f'PDS: got {pds_got} vs expected {best_step}'
print('OK')
" > log_eta_pds_check.txt 2>&1
check_contains "OK" log_eta_pds_check.txt


# --- 4. Error cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: --version"
stb-oerAnalysis --version > log_version.txt 2>&1
check_contains "stb-oerAnalysis" log_version.txt

echo "Testing: --help documents --temp/--force-tolerance and the PDS/eta terminology"
stb-oerAnalysis --help > log_help.txt 2>&1
check_contains "temp" log_help.txt
check_contains "force-tolerance" log_help.txt
check_contains "potential-determining step" log_help.txt
check_contains "4.92 eV" log_help.txt


# --- 5. Interactive path (stb-suite, shortcut 4.14.4) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.14.4) ---"

echo "Testing: navigate 4.14.4 -> defaults -> quit"
rm -rf oer_study
write_common_energies
write_local_zpe_fixture "OH" "2, 3" "'O', 'H'" "oer_zpe_oh" 8.0 5.0
write_local_zpe_fixture "O" "2" "'O'" "oer_zpe_o" 8.0 5.0
write_local_zpe_fixture "OOH" "2, 3, 4" "'O', 'O', 'H'" "oer_zpe_ooh" 8.0 5.0
write_local_zpe_fixture "H2O" "0, 1, 2" "'O', 'H', 'H'" "oer_zpe_h2o" 8.0 5.0
{
  echo "4.14.4"
  echo ""               # run_dir (default oer_study)
  echo ""               # output_filename (default calc.out)
  echo ""               # temperature (default 298.15)
  echo ""               # force_tolerance (default 0.05)
  echo ""               # press enter to continue
  echo "0"              # quit stage submenu
} | stb-suite > log_menu.txt 2>&1
check_contains "eta" log_menu.txt
check_success oer_study/OER_report.txt


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
