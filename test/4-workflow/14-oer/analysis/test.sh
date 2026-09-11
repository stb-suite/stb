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

check_not_contains() {
    if grep -q -- "$1" "$2" 2>/dev/null; then
        echo -e "   -> ${RED}Failed:${NC} '$1' unexpectedly found in '$2'"
        FAIL=$((FAIL+1))
    else
        echo -e "   -> ${GREEN}Verified:${NC} '$1' NOT found in '$2' (as expected)"
        PASS=$((PASS+1))
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
# reads back (ZPE mode) or merely checks exist (Stage 2's own report --
# O*/OOH* are always derived from the winning OH* site now, no strategy
# line left to parse; Stage 3's own BSSE mode line is gone too -- BSSE is
# always 3 separate per-intermediate triads now, no shared/full choice).
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

# 3 separate per-intermediate BSSE triads, small known round shifts (same
# style as HER's own bsse_slab=0.05/bsse_ads=0.05 fixture) -- slab/ads
# SPLIT deliberately differs per intermediate (0.05/0.05, 0.08/0.08,
# 0.03/0.03) to prove 3 genuinely distinct triads are read (not one
# reused triad, the old --bsse-mode shared behavior, now removed), while
# each triad's TOTAL still cancels to exactly 0.0 so every downstream
# dG/eta assertion below (independent hand-computation etc.) needs no
# change -- E_corrected == E_raw for all 3 intermediates either way.
write_calc_out("oer_study/05_bsse_OH_slab_only/calc.out", -520.050000)
write_calc_out("oer_study/05_bsse_OH_slab_ghost/calc.out", -520.000000)
write_calc_out("oer_study/05_bsse_OH_adsorbate_ghost_slab/calc.out", -120.050000)
write_calc_out("oer_study/05_bsse_OH_isolated/calc.out", -120.000000)
write_calc_out("oer_study/05_bsse_O_slab_only/calc.out", -520.080000)
write_calc_out("oer_study/05_bsse_O_slab_ghost/calc.out", -520.000000)
write_calc_out("oer_study/05_bsse_O_adsorbate_ghost_slab/calc.out", -120.080000)
write_calc_out("oer_study/05_bsse_O_isolated/calc.out", -120.000000)
write_calc_out("oer_study/05_bsse_OOH_slab_only/calc.out", -520.030000)
write_calc_out("oer_study/05_bsse_OOH_slab_ghost/calc.out", -520.000000)
write_calc_out("oer_study/05_bsse_OOH_adsorbate_ghost_slab/calc.out", -120.030000)
write_calc_out("oer_study/05_bsse_OOH_isolated/calc.out", -120.000000)

with open("oer_study/oer_stage2.txt", "w") as f:
    f.write("===== OER STAGE 2 REPORT (O*/OOH* INTERMEDIATES) =====\n")

with open("oer_study/oer_stage3.txt", "w") as f:
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

stb-oerAnalysis --directory oer_study --temp 298.15 --save-report --no-intro > log_local.txt 2>&1
check_exit_code $? 0
check_contains "Local (partial-Hessian, 1-atom) vibrational modes" log_local.txt
check_contains "Local (partial-Hessian, 2-atom) vibrational modes" log_local.txt
check_contains "Local (partial-Hessian, 3-atom) vibrational modes" log_local.txt
check_success oer_study/oer_stage4.txt
check_success oer_study/OER_report.txt

echo "Testing: consolidated Markdown report (OER_report.md) is always written, with a real"
echo "         embedded PNG (fig.savefig, headless-safe) and every section's data"
check_contains "\[Saved\].*OER_report.md" log_local.txt
check_success oer_study/OER_report.md
check_success oer_study/plot/oer_free_energy_diagram.png
check_contains "# OER Study Report" oer_study/OER_report.md
check_contains "## 1. Overview" oer_study/OER_report.md
check_contains "## 2. Electronic Energies (FreeEng)" oer_study/OER_report.md
check_contains "## 3. BSSE Correction" oer_study/OER_report.md
check_contains "## 4. Thermal Correction" oer_study/OER_report.md
check_contains "## 5. Gas-Phase O2 Reference" oer_study/OER_report.md
check_contains "## 6. Reaction Steps & Potential-Determining Step" oer_study/OER_report.md
check_contains "## 7. Free-Energy Diagram" oer_study/OER_report.md
check_contains "## 8. Final Result" oer_study/OER_report.md
check_contains "## 9. Files" oer_study/OER_report.md
check_contains "!\[OER free-energy diagram\](plot/oer_free_energy_diagram.png)" oer_study/OER_report.md
check_contains "\*\*<-- PDS\*\*" oer_study/OER_report.md

echo "Testing: OER_report.md's eta/PDS match the plain-text report exactly (same run, same numbers)"
python3 -c "
import re
txt = open('oer_study/oer_stage4.txt').read()
md = open('oer_study/OER_report.md').read()
m_txt = re.search(r'overpotential eta = ([\-+\d.]+) V', txt)
m_md = re.search(r'Theoretical overpotential eta.*?\*\*([\-+\d.]+) V\*\*', md)
assert m_txt and m_md, 'could not find eta in one of the two reports'
assert abs(float(m_txt.group(1)) - float(m_md.group(1))) < 1e-9, \
    f'eta mismatch: txt={m_txt.group(1)} vs md={m_md.group(1)}'
pds_txt = re.search(r'Potential-determining step \(PDS\) = Step (\d+)', txt)
pds_md = re.search(r'Potential-determining step.*?Step (\d+)', md)
assert pds_txt and pds_md, 'could not find PDS in one of the two reports'
assert pds_txt.group(1) == pds_md.group(1), \
    f'PDS mismatch: txt=Step {pds_txt.group(1)} vs md=Step {pds_md.group(1)}'
print('OK')
" > log_md_consistency_check.txt 2>&1
check_contains "OK" log_md_consistency_check.txt

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

echo "Testing: dG1..dG4/eta match a FULLY independent hand-computation (correct T*S, i.e."
echo "         NOT multiplied by --temp again -- this is the exact regression check for a"
echo "         real bug found 2026-09-07: every g_* line in oer_analysis.py's main() used to"
echo "         do '- args.temp * d_ts_*', double-counting temperature into thermal_term()'s"
echo "         own already-T-integrated ts_ev return (same convention as HER's her_analysis.py,"
echo "         which never re-multiplies) -- inflated eta by a factor of ~298 at room temp."
python3 -c "
import re
import numpy as np

BOLTZMANN_EV_K = 8.617333262e-5
FREQ_CONV = 15.633302
EV_PER_THZ = 0.00413566733
T = 298.15
MASS = {'O': 15.999, 'H': 1.00794}
K = {'O': 8.0, 'H': 5.0}

def mode_zpe_ts(symbol, temperature_k=T):
    freq_thz = FREQ_CONV * np.sqrt(K[symbol] / MASS[symbol])
    e_mode = freq_thz * EV_PER_THZ
    zpe = 0.5 * e_mode
    x = e_mode / (BOLTZMANN_EV_K * temperature_k)
    S = BOLTZMANN_EV_K * (x / np.expm1(x) - np.log1p(-np.exp(-x)))
    ts = temperature_k * S
    return zpe, ts

def zpe_ts_for(symbols):
    # local mode: 3 isotropic (degenerate) modes per local atom, decoupled
    # (diagonal fixture) -- matches oer_analysis.py's compute_local_zpe_entropy
    # applied to this test's own diagonal/isotropic k_O=8.0/k_H=5.0 fixture.
    zpe = ts = 0.0
    for sym in symbols:
        z, t = mode_zpe_ts(sym)
        zpe += 3 * z
        ts += 3 * t
    return zpe, ts

zpe_oh, ts_oh = zpe_ts_for(['O', 'H'])
zpe_o, ts_o = zpe_ts_for(['O'])
zpe_ooh, ts_ooh = zpe_ts_for(['O', 'O', 'H'])
zpe_h2o, ts_h2o = zpe_ts_for(['O', 'H', 'H'])

H2_ZPE_EV, H2_TS_298K_EV = 0.270, 0.400  # fixed literature constants, oer_analysis.py's own

# Known fixture energies (write_common_energies above) + each intermediate's
# own BSSE triad (slab/ads split differs per intermediate, but every triad's
# TOTAL cancels to exactly 0.0 by construction -- see write_common_energies).
e_clean, e_oh, e_o, e_ooh = -400.000000, -499.950000, -480.000000, -520.000000
e_h2, e_h2o = -31.500000, -470.000000
bsse_total = 0.0
e_oh_corr, e_o_corr, e_ooh_corr = e_oh + bsse_total, e_o + bsse_total, e_ooh + bsse_total

# G = E(corrected) + ZPE - T*S -- ts_* above is ALREADY T*S (eV), no further
# multiplication by T -- this is the exact step the real bug got wrong.
g_clean = e_clean
g_h2 = e_h2 + H2_ZPE_EV - H2_TS_298K_EV
g_h2o = e_h2o + zpe_h2o - ts_h2o
g_oh = e_oh_corr + zpe_oh - ts_oh
g_o = e_o_corr + zpe_o - ts_o
g_ooh = e_ooh_corr + zpe_ooh - ts_ooh
g_o2 = 2.0 * g_h2o - 2.0 * g_h2 + 4.920

dG1 = g_oh + 0.5 * g_h2 - g_clean - g_h2o
dG2 = g_o + 0.5 * g_h2 - g_oh
dG3 = g_ooh + 0.5 * g_h2 - g_o - g_h2o
dG4 = g_clean + g_o2 - g_ooh + 0.5 * g_h2
eta_expected = max(dG1, dG2, dG3, dG4) - 1.23

with open('oer_study/oer_stage4.txt') as f:
    text = f.read()
rows = re.findall(r'Step (\d):.*?dG = ([\-+\d.]+) eV', text)
dg_got = {int(k): float(v) for k, v in rows}
expected = {1: dG1, 2: dG2, 3: dG3, 4: dG4}
for step in (1, 2, 3, 4):
    assert abs(dg_got[step] - expected[step]) < 1e-3, \
        f'dG{step}: got {dg_got[step]} vs independently hand-computed {expected[step]}'

m_eta = re.search(r'overpotential eta = ([\-+\d.]+) V', text)
eta_got = float(m_eta.group(1))
assert abs(eta_got - eta_expected) < 1e-3, \
    f'eta: got {eta_got} vs independently hand-computed {eta_expected} (temperature ' \
    f'double-counting regression: with the old bug this would be off by a factor of ~{T:.0f}x)'
print('OK')
" > log_independent_eta_check.txt 2>&1
check_contains "OK" log_independent_eta_check.txt


# --- 3b. --plot (gnuplot free-energy diagram) ---
echo -e "\n--- Testing --plot (gnuplot .dat/.gplot free-energy diagram) ---"
stb-oerAnalysis --directory oer_study --plot --no-show --no-intro > log_plot.txt 2>&1
check_exit_code $? 0
check_contains "\[Saved\]" log_plot.txt
check_success oer_study/plot/oer_free_energy_diagram.dat
check_success oer_study/plot/oer_free_energy_diagram.gplot

echo "Testing: .dat holds the plateau/CG-step form -- 4 index blocks (U0 plateaus, U0"
echo "         connectors, Ueq plateaus, Ueq connectors), and .gplot's xtics carry the"
echo "         5 state labels (moved out of the .dat itself in the plateau redesign)"
check_contains '"\* + H2O"' oer_study/plot/oer_free_energy_diagram.gplot
check_contains '"\* + O2"' oer_study/plot/oer_free_energy_diagram.gplot
check_contains '"OH\*"' oer_study/plot/oer_free_energy_diagram.gplot
check_contains '"OOH\*"' oer_study/plot/oer_free_energy_diagram.gplot

python3 -c "
with open('oer_study/plot/oer_free_energy_diagram.dat') as f:
    text = f.read()
# strip comment/header lines, keep blank-line structure (gnuplot: 1 blank line
# breaks a trace, 2+ blank lines start a new 'index' block)
body = '\n'.join(l for l in text.splitlines() if not l.startswith('#'))
blocks = [b for b in body.split('\n\n\n') if b.strip()]
assert len(blocks) == 4, f'expected 4 index blocks (plateaus/connectors x U0/Ueq), got {len(blocks)}'

def segments(block):
    return [seg for seg in block.strip('\n').split('\n\n') if seg.strip()]

u0_plateaus, u0_connectors, ueq_plateaus, ueq_connectors = [segments(b) for b in blocks]
assert len(u0_plateaus) == 5, f'expected 5 U0 plateau segments (1/state), got {len(u0_plateaus)}'
assert len(u0_connectors) == 4, f'expected 4 U0 connector segments (between consecutive states), got {len(u0_connectors)}'
assert len(ueq_plateaus) == 5 and len(ueq_connectors) == 4

def plateau_y(seg):
    p1, p2 = seg.strip().split('\n')
    y1, y2 = float(p1.split()[1]), float(p2.split()[1])
    assert abs(y1 - y2) < 1e-9, f'plateau segment should be flat (same y at both x): {seg}'
    return y1

g_u0 = [plateau_y(s) for s in u0_plateaus]
g_ueq = [plateau_y(s) for s in ueq_plateaus]
assert abs(g_u0[0]) < 1e-9, 'first state (* + H2O) should be G=0 at U=0'
assert abs(g_u0[-1] - 4.92) < 1e-3, f'last state (* + O2) should be G=4.92 eV at U=0, got {g_u0[-1]}'
assert abs(g_ueq[-1]) < 1e-3, f'last state should return to G=0 at U=1.23V (equilibrium), got {g_ueq[-1]}'
for n in range(5):
    assert abs(g_ueq[n] - (g_u0[n] - n * 1.23)) < 1e-9, f'state {n}: G_Ueq must equal G_U0 - n*1.23'
print('OK')
" > log_plot_check.txt 2>&1
check_contains "OK" log_plot_check.txt

echo "Testing: .gplot shades the PDS step, draws 4 index blocks (solid plateaus + dashed"
echo "         connectors, 2 colors), explicit x (no histogram style)"
check_contains "set object 1 rect from" oer_study/plot/oer_free_energy_diagram.gplot
check_contains "index 0 using 1:2 with lines" oer_study/plot/oer_free_energy_diagram.gplot
check_contains "index 1 using 1:2 with lines lw 1.5 dt 2" oer_study/plot/oer_free_energy_diagram.gplot
check_contains "index 2 using 1:2 with lines" oer_study/plot/oer_free_energy_diagram.gplot
check_contains "index 3 using 1:2 with lines lw 1.5 dt 2" oer_study/plot/oer_free_energy_diagram.gplot
check_not_contains "with linespoints" oer_study/plot/oer_free_energy_diagram.gplot
check_contains "set yrange \[" oer_study/plot/oer_free_energy_diagram.gplot

if command -v gnuplot > /dev/null 2>&1; then
    echo "Testing: gnuplot actually renders the .gplot into a PDF (real render, not just a syntax check)"
    ( cd oer_study/plot && gnuplot oer_free_energy_diagram.gplot ) > log_gnuplot_render.txt 2>&1
    check_exit_code $? 0
    check_success oer_study/plot/oer_free_energy_diagram.pdf
else
    echo "gnuplot not installed -- skipping the real-render check (files-written checks above still apply)."
fi

echo "Testing: --plot ALSO saves the BSSE-correction chart (.dat/.gplot), same invocation"
check_success oer_study/plot/oer_bsse_correction.dat
check_success oer_study/plot/oer_bsse_correction.gplot
check_contains "OH\* (slab)" oer_study/plot/oer_bsse_correction.dat
check_contains "OH\* (total)" oer_study/plot/oer_bsse_correction.dat
check_contains "OOH\* (total)" oer_study/plot/oer_bsse_correction.dat
check_contains "using 1:2:4:xtic(3) with boxes lc rgb variable" oer_study/plot/oer_bsse_correction.gplot

echo "Testing: BSSE .dat's 9 rows (3 intermediates x slab/ads./total) match [2]'s own printed"
echo "         BSSE table exactly -- and the 3 SLAB components genuinely differ (0.05/0.08/0.03"
echo "         eV per the fixture), proving 3 separate triads are read, not one shared/reused"
echo "         triad (the old --bsse-mode shared behavior, now removed)"
python3 -c "
with open('oer_study/plot/oer_bsse_correction.dat') as f:
    rows = [l.split(None, 3) for l in f if l.strip() and not l.startswith('#')]
assert len(rows) == 9, f'expected 9 bars (3 intermediates x slab/ads/total), got {len(rows)}'
slabs = [float(rows[i][1]) for i in (0, 3, 6)]  # OH*/O*/OOH* slab rows, 0-indexed
assert len(set(round(s, 6) for s in slabs)) == 3, \
    f'expected 3 DISTINCT slab components (separate triads), got {slabs}'
totals = [float(rows[i][1]) for i in (2, 5, 8)]  # OH*/O*/OOH* total rows
assert all(abs(t - 0.0) < 1e-9 for t in totals), \
    f'this fixture is built so every triad total cancels to 0.0: {totals}'
for i in (0, 3, 6):  # slab rows
    slab, ads, total = float(rows[i][1]), float(rows[i+1][1]), float(rows[i+2][1])
    assert abs((slab + ads) - total) < 1e-9, f'slab+ads should equal total at row {i}'
print('OK')
" > log_bsse_plot_check.txt 2>&1
check_contains "OK" log_bsse_plot_check.txt

if command -v gnuplot > /dev/null 2>&1; then
    echo "Testing: gnuplot actually renders the BSSE .gplot into a PDF (real render)"
    ( cd oer_study/plot && gnuplot oer_bsse_correction.gplot ) > log_gnuplot_bsse_render.txt 2>&1
    check_exit_code $? 0
    check_success oer_study/plot/oer_bsse_correction.pdf
fi

echo "Testing: --plot ALSO saves the raw-vs-corrected comparison chart (3 small-multiple"
echo "         panels, one per intermediate, each own y-scale -- see its own docstring for"
echo "         why a shared/zero-based axis can't show a ~0.4 eV shift on a ~-3700 eV base)"
check_success oer_study/plot/oer_bsse_raw_corrected.dat
check_success oer_study/plot/oer_bsse_raw_corrected.gplot
check_contains "set multiplot layout 1,3" oer_study/plot/oer_bsse_raw_corrected.gplot
check_contains "\"OH\*\"" oer_study/plot/oer_bsse_raw_corrected.dat
check_contains "\"OOH\*\"" oer_study/plot/oer_bsse_raw_corrected.dat

echo "Testing: raw-vs-corrected panels are drawn as bars (with boxes), not points"
check_contains "with boxes lc rgb \"#2255cc\"" oer_study/plot/oer_bsse_raw_corrected.gplot
check_contains "with boxes lc rgb \"#22aa55\"" oer_study/plot/oer_bsse_raw_corrected.gplot
check_contains "set boxwidth" oer_study/plot/oer_bsse_raw_corrected.gplot
check_not_contains "with points" oer_study/plot/oer_bsse_raw_corrected.gplot

python3 -c "
with open('oer_study/plot/oer_bsse_raw_corrected.dat') as f:
    rows = [l.split() for l in f if l.strip() and not l.startswith('#')]
assert len(rows) == 3, f'expected 1 row per intermediate (OH*/O*/OOH*), got {len(rows)}'
for label, e_raw, bsse, e_corr in rows:
    e_raw, bsse, e_corr = float(e_raw), float(bsse), float(e_corr)
    assert abs((e_raw + bsse) - e_corr) < 1e-6, f'{label}: raw + BSSE should equal corrected'
print('OK')
" > log_bsse_rc_check.txt 2>&1
check_contains "OK" log_bsse_rc_check.txt

if command -v gnuplot > /dev/null 2>&1; then
    echo "Testing: gnuplot actually renders the raw-vs-corrected multiplot into a PDF"
    ( cd oer_study/plot && gnuplot oer_bsse_raw_corrected.gplot ) > log_gnuplot_bsse_rc_render.txt 2>&1
    check_exit_code $? 0
    check_success oer_study/plot/oer_bsse_raw_corrected.pdf
fi

echo "Testing: the Markdown report always embeds all 3 chart PNGs (free energy, BSSE"
echo "         components, BSSE raw-vs-corrected), unconditional, independent of --plot/--show"
check_success oer_study/plot/oer_bsse_correction.png
check_success oer_study/plot/oer_bsse_raw_corrected.png
check_contains "OER BSSE correction: raw vs. corrected" oer_study/OER_report.md
check_contains "OER BSSE correction: slab/adsorbate components" oer_study/OER_report.md
check_contains "plot/oer_bsse_raw_corrected.png" oer_study/OER_report.md
check_contains "plot/oer_bsse_correction.png" oer_study/OER_report.md

echo "Testing: --no-plot --no-show skip both, no interactive prompt hang"
stb-oerAnalysis --directory oer_study --no-plot --no-show --no-intro > log_noplot.txt 2>&1
check_exit_code $? 0
check_contains "No plot generated" log_noplot.txt

echo "Testing: --help documents --plot/--no-plot/--show/--no-show"
stb-oerAnalysis --help > log_help_plot.txt 2>&1
check_contains "no-plot" log_help_plot.txt
check_contains "no-show" log_help_plot.txt


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
check_contains "lattice oxygen evolution mechanism" log_help.txt

echo "Testing: final report includes the AEM/LOER limitation note"
check_contains "LIMITATION" oer_study/oer_stage4.txt


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
