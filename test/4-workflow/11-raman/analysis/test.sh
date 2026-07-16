#!/bin/bash

# --- Setup ---
# Smoke test for stb-ramanAnalysis (Raman Spectrum, Stage 3: Analysis, item
# 4.11.3). Chains real stb-raman + stb-ramanModes runs (Stages 1-2, tested
# separately) to get a genuine raman_study/ tree with a # MODE_TABLE, then
# fabricates SystemLabel.EPSIMG + calc.out in every optical_disp/mode_*/
# folder from a known analytic Lorentz oscillator (eps2 in closed form),
# shifted slightly per mode/sign/axis so the finite-difference Raman tensor
# comes out non-degenerate -- physically arbitrary, but format-correct and
# enough to exercise the Kramers-Kronig static-dielectric readout, the
# central finite difference, and the spectrum/report wiring end to end.
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREP_DIR="$(cd "$FIXTURE_DIR/../prep" && pwd)"
MODES_DIR="$(cd "$FIXTURE_DIR/../modes" && pwd)"
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

# Runs stb-raman + stb-ramanModes (Stages 1-2, tested separately) with the
# given selected modes (default: 1 2 3, all imaginary in this fixture --
# pass e.g. "1 39" to also include mode 39, the one genuinely real/positive
# -frequency vibration in this fixture, when a test needs one), then
# fabricates SystemLabel.EPSIMG + calc.out in every resulting
# optical_disp/mode_*/ folder.
make_optical_disp() {
    local modes="${1:-1 2 3}"
    rm -rf raman_study
    stb-raman -s structure.fdf -c calc.fdf -p . -dim 1 1 1 --no-intro > /dev/null 2>&1
    python3 - <<'PYEOF'
import glob
with open("Sn3O4.FA") as f:
    lines = f.readlines()
header, rows = lines[0], lines[1:]
for i, d in enumerate(sorted(glob.glob("raman_study/phonon_disp/disp-*")), start=1):
    scale = 1.0 + 0.01 * (i % 17)
    with open(f"{d}/Sn3O4.FA", "w") as out:
        out.write(header)
        for row in rows:
            parts = row.split()
            idx = parts[0]
            vals = [float(v) * scale for v in parts[1:]]
            out.write(f"{idx:>6}  {vals[0]: .9E}  {vals[1]: .9E}  {vals[2]: .9E}\n")
PYEOF
    stb-ramanModes --directory raman_study --calc calc.fdf --modes $modes --no-intro > /dev/null 2>&1
    python3 - <<'PYEOF'
import glob, os
import numpy as np

w0_by_axis = {"x": 3.0, "y": 3.5, "z": 4.0}
wp, gamma = 5.0, 0.3
w = np.linspace(0.05, 15, 400)

for d in sorted(glob.glob("raman_study/optical_disp/mode_*")):
    name = os.path.basename(d)
    parts = name.split("_")
    axis, sign, mode_idx = parts[-1], parts[-2], int(parts[1])
    delta_shift = 0.01 if sign == "plus" else -0.01
    w0 = w0_by_axis[axis] + delta_shift * mode_idx
    denom = (w0**2 - w**2)**2 + (gamma*w)**2
    eps2 = wp**2 * gamma * w / denom
    with open(os.path.join(d, "Sn3O4.EPSIMG"), "w") as f:
        f.write("## Minimum and maximum energy in eV\n")
        f.write(f"##  {w[0]:.6f}  {w[-1]:.6f}\n")
        f.write("## Number of spin components\n")
        f.write("## 1\n")
        for wi, e2i in zip(w, eps2):
            f.write(f"{wi:.6f}  {e2i:.6f}\n")
    with open(os.path.join(d, "calc.out"), "w") as f:
        f.write("siesta: SCF Convergence by DM criterion\n")
        f.write("SCF cycle converged after 12 iterations\n")
PYEOF
}

# Runs stb-raman + stb-ramanModes --full-tensor (1 mode -> 12 folders), then
# fabricates SystemLabel.EPSIMG so that its Kramers-Kronig static value
# encodes a KNOWN ground-truth symmetric matrix M (with off-diagonal terms)
# along each probed direction: target_static = BASE + sign*delta*(n^T M n),
# solved into a Lorentz oscillator (BASE=5.0 keeps eps1 safely positive for
# both signs; it cancels out of the +/-delta central difference, so the
# recovered Rxx..Ryz only depend on M, not BASE). Lets the test assert the
# recovered tensor against the exact known M, not just "some plausible
# numbers", for the new off-diagonal-recovery math specifically.
make_optical_disp_full() {
    rm -rf raman_study
    stb-raman -s structure.fdf -c calc.fdf -p . -dim 1 1 1 --no-intro > /dev/null 2>&1
    python3 - <<'PYEOF'
import glob
with open("Sn3O4.FA") as f:
    lines = f.readlines()
header, rows = lines[0], lines[1:]
for i, d in enumerate(sorted(glob.glob("raman_study/phonon_disp/disp-*")), start=1):
    scale = 1.0 + 0.01 * (i % 17)
    with open(f"{d}/Sn3O4.FA", "w") as out:
        out.write(header)
        for row in rows:
            parts = row.split()
            idx = parts[0]
            vals = [float(v) * scale for v in parts[1:]]
            out.write(f"{idx:>6}  {vals[0]: .9E}  {vals[1]: .9E}  {vals[2]: .9E}\n")
PYEOF
    stb-ramanModes --directory raman_study --calc calc.fdf --modes 1 --full-tensor --no-intro > /dev/null 2>&1
    python3 - <<'PYEOF'
import glob, os
import numpy as np

M = np.array([
    [2.0, 0.3, -0.5],
    [0.3, 1.5,  0.2],
    [-0.5, 0.2, 1.0],
])
np.save("ground_truth_M.npy", M)
delta = 0.02
BASE = 5.0
inv_sqrt2 = 0.70710678118654752440
directions = {
    "x": np.array([1, 0, 0]), "y": np.array([0, 1, 0]), "z": np.array([0, 0, 1]),
    "xy": np.array([inv_sqrt2, inv_sqrt2, 0]),
    "xz": np.array([inv_sqrt2, 0, inv_sqrt2]),
    "yz": np.array([0, inv_sqrt2, inv_sqrt2]),
}
wp, gamma = 5.0, 0.3
w = np.linspace(0.05, 15, 400)

for d in sorted(glob.glob("raman_study/optical_disp/mode_*")):
    name = os.path.basename(d)
    parts = name.split("_")
    axis, sign = parts[-1], parts[-2]
    n = directions[axis]
    quad = n @ M @ n
    sign_factor = 1.0 if sign == "plus" else -1.0
    target_static = BASE + sign_factor * delta * quad
    w0 = (wp**2 / (target_static - 1.0)) ** 0.5
    denom = (w0**2 - w**2)**2 + (gamma*w)**2
    eps2 = wp**2 * gamma * w / denom
    with open(os.path.join(d, "Sn3O4.EPSIMG"), "w") as f:
        f.write("## Minimum and maximum energy in eV\n")
        f.write(f"##  {w[0]:.6f}  {w[-1]:.6f}\n")
        f.write("## Number of spin components\n")
        f.write("## 1\n")
        for wi, e2i in zip(w, eps2):
            f.write(f"{wi:.6f}  {e2i:.6f}\n")
    with open(os.path.join(d, "calc.out"), "w") as f:
        f.write("siesta: SCF Convergence by DM criterion\n")
        f.write("SCF cycle converged after 12 iterations\n")
PYEOF
}


# --- 1. Preparation ---
echo "--- Starting tester for STB-RamanAnalysis stage 3: analysis (item 4.11.3) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$PREP_DIR/structure.fdf" "$TEST_DIR/"
cp "$PREP_DIR/calc.fdf" "$TEST_DIR/"
cp "$PREP_DIR/O.psf" "$TEST_DIR/"
cp "$PREP_DIR/Sn.psf" "$TEST_DIR/"
cp "$MODES_DIR/Sn3O4.FA" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Missing Stage-2 output guard ---
echo -e "\n--- Testing the missing-Stage-2-output guard ---"
rm -rf raman_study
stb-ramanAnalysis --directory raman_study --no-intro > log_no_stage2.txt 2>&1
check_exit_code $? 1
check_contains "run stb-ramanModes" log_no_stage2.txt


# --- 3. Normal analysis ---
echo -e "\n--- Testing default analysis ---"
make_optical_disp
stb-ramanAnalysis --directory raman_study --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "Modes in table    : 3" log_basic.txt
check_contains "Modes analyzed      : 3/3" log_basic.txt
check_contains "RAMAN-ACTIVE MODES SUMMARY" log_basic.txt
check_contains "diagonal Raman tensor" log_basic.txt
check_success raman_study/raman_stage3.txt
check_success raman_study/raman_spectrum.dat
check_success raman_study/raman_spectrum.gplot
check_contains "\[0\] RUN METADATA" raman_study/raman_stage3.txt
check_contains "\[1\] READING RUNS" raman_study/raman_stage3.txt
check_contains "\[2\] RAMAN-ACTIVE MODES SUMMARY" raman_study/raman_stage3.txt
check_contains "\[3\] SPECTRUM" raman_study/raman_stage3.txt
check_contains "\[4\] SUMMARY & FILES" raman_study/raman_stage3.txt

echo "Testing: Rxx/Ryy/Rzz differ per mode (finite-difference across genuinely different fabricated dielectric data)"
python3 -c "
with open('raman_study/raman_stage3.txt') as f:
    text = f.read()
assert 'Rxx=' in text
print('OK')
" > log_tensor_check.txt 2>&1
check_contains "OK" log_tensor_check.txt

echo "Testing: spectrum .dat has the expected 2-column header"
check_contains "frequency.cm" raman_study/raman_spectrum.dat

echo "Testing: depolarization ratio (rho) reported with a polarized/depolarized tag"
check_contains "rho~" raman_study/raman_stage3.txt
python3 -c "
import re
with open('raman_study/raman_stage3.txt') as f:
    text = f.read()
matches = re.findall(r'rho~([\d.]+) \((polarized|depolarized)\)', text)
assert matches, 'no rho values found'
for rho_str, tag in matches:
    rho = float(rho_str)
    assert 0.0 <= rho <= 0.75 + 1e-6, f'rho {rho} out of physical range [0, 0.75]'
    expected_tag = 'polarized' if rho < 0.75 - 1e-6 else 'depolarized'
    assert tag == expected_tag, f'rho={rho} tagged {tag}, expected {expected_tag}'
print('OK')
" > log_rho_check.txt 2>&1
check_contains "OK" log_rho_check.txt


# --- 3b. --temperature: Bose-Einstein weighting changes the spectrum ---
# Modes "1 39": mode 1 is imaginary (must be left unweighted, see the
# [WARNING] below), mode 39 is the one genuinely real/positive-frequency
# vibration in this fixture (2.7882 THz) -- needed so the weighting
# actually has something real to change.
echo -e "\n--- Testing --temperature (Bose-Einstein thermal weighting) ---"
make_optical_disp "1 39"
stb-ramanAnalysis --directory raman_study --temperature 300 -o raman_spectrum_300k --no-intro \
    > log_temperature.txt 2>&1
check_exit_code $? 0
check_contains "Thermal weighting : Bose-Einstein Stokes factor (n+1) at 300.0 K" log_temperature.txt
check_contains "imaginary-frequency mode(s) left unweighted" log_temperature.txt
check_success raman_study/raman_spectrum_300k.dat

echo "Testing: weighted spectrum differs from the unweighted default (T=300K changes intensities)"
stb-ramanAnalysis --directory raman_study -o raman_spectrum_unweighted --no-intro > /dev/null 2>&1
python3 -c "
import numpy as np
weighted = np.loadtxt('raman_study/raman_spectrum_300k.dat')
unweighted = np.loadtxt('raman_study/raman_spectrum_unweighted.dat')
assert weighted.shape == unweighted.shape
assert not np.allclose(weighted[:, 1], unweighted[:, 1]), 'thermal weighting had no effect'
# Bose weight (n+1) is always >= 1, so every weighted intensity must be >= the unweighted one.
assert (weighted[:, 1] >= unweighted[:, 1] - 1e-9).all(), 'weighted intensity fell below unweighted'
print('OK')
" > log_temperature_check.txt 2>&1
check_contains "OK" log_temperature_check.txt


# --- 3c. --experimental: overlay + nearest-neighbor peak match ---
# Mode 39 alone -> exactly one simulated peak, at 2.7882 THz * 33.35641 =
# ~92.99 cm^-1 (known from this fixture, see make_optical_disp's comment).
# Fabricates a synthetic experimental spectrum with ONE clean Gaussian
# peak at a nearby, deliberately different position (95.0 cm^-1) so the
# match has a known, checkable non-zero delta (~2.0 cm^-1).
echo -e "\n--- Testing --experimental (peak overlay + matching) ---"
make_optical_disp "39"
python3 -c "
import numpy as np
freq = np.linspace(50, 140, 500)
intensity = 100.0 * np.exp(-0.5 * ((freq - 95.0) / 3.0) ** 2)
with open('fake_experimental.dat', 'w') as f:
    f.write('# fake experimental spectrum for testing\n')
    for fr, it in zip(freq, intensity):
        f.write(f'{fr:.4f} {it:.6f}\n')
"
stb-ramanAnalysis --directory raman_study --experimental fake_experimental.dat --no-intro \
    > log_experimental.txt 2>&1
check_exit_code $? 0
check_contains "\[3b\] EXPERIMENTAL COMPARISON" log_experimental.txt
check_contains "Peaks found       : 1 simulated, 1 experimental" log_experimental.txt
check_success raman_study/raman_spectrum_experimental.dat

echo "Testing: matched peak delta is close to the known ~2.0 cm^-1 offset (95.0 exp vs ~92.99 sim)"
python3 -c "
import re
with open('raman_study/raman_stage3.txt') as f:
    text = f.read()
section = text.split('EXPERIMENTAL COMPARISON')[1]
m = re.search(r'(\d+\.\d+)\s+(\d+\.\d+)\s+([+-]\d+\.\d+)', section)
assert m, 'could not find the peak-match row'
exp_f, sim_f, delta = float(m.group(1)), float(m.group(2)), float(m.group(3))
assert abs(exp_f - 95.0) < 0.5, f'expected exp peak near 95.0, got {exp_f}'
assert abs(sim_f - 92.99) < 1.0, f'expected sim peak near 92.99, got {sim_f}'
assert abs(delta - (sim_f - exp_f)) < 1e-6
print('OK')
" > log_experimental_check.txt 2>&1
check_contains "OK" log_experimental_check.txt

echo "Testing: missing --experimental file is a clear error, not a crash"
stb-ramanAnalysis --directory raman_study --experimental does_not_exist.dat --no-intro \
    > log_experimental_missing.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_experimental_missing.txt


# --- 4. Missing .EPSIMG in one folder -> that mode is skipped, others still analyzed ---
echo -e "\n--- Testing graceful skip of an incomplete mode (missing .EPSIMG) ---"
make_optical_disp
rm -f raman_study/optical_disp/mode_02_plus_x/Sn3O4.EPSIMG
stb-ramanAnalysis --directory raman_study --no-intro > log_partial.txt 2>&1
check_exit_code $? 0
check_contains "SKIP (incomplete" log_partial.txt
check_contains "Modes analyzed      : 2/3" log_partial.txt


# --- 4b. --full-tensor analysis: recovered tensor must match the known ground truth M ---
echo -e "\n--- Testing full-tensor analysis against a known ground-truth matrix ---"
make_optical_disp_full
stb-ramanAnalysis --directory raman_study --no-intro > log_fulltensor.txt 2>&1
check_exit_code $? 0
check_contains "Rxy=" log_fulltensor.txt
check_contains "FULL Raman tensor" log_fulltensor.txt

echo "Testing: recovered Rxx/Ryy/Rzz/Rxy/Rxz/Ryz match the known ground-truth M within 1%"
python3 -c "
import re
import numpy as np
M = np.load('ground_truth_M.npy')
with open('raman_study/raman_stage3.txt') as f:
    text = f.read()
m = re.search(
    r'Rxx=([\-\d.]+)\s+Ryy=([\-\d.]+)\s+Rzz=([\-\d.]+)\s+Rxy=([\-\d.]+)\s+Rxz=([\-\d.]+)\s+Ryz=([\-\d.]+)',
    text)
assert m, 'could not find the full tensor line in the report'
rxx, ryy, rzz, rxy, rxz, ryz = (float(v) for v in m.groups())
recovered = np.array([[rxx, rxy, rxz], [rxy, ryy, ryz], [rxz, ryz, rzz]])
assert np.allclose(recovered, M, atol=0.02), f'recovered {recovered} vs ground truth {M}'
print('OK')
" > log_fulltensor_check.txt 2>&1
check_contains "OK" log_fulltensor_check.txt


# --- 4c. Reconstruction from a stb-ramanModes --reduce-components sidecar ---
# Hand-crafts a mode_XX_tensor_form.json (a KNOWN, purely off-diagonal T0 --
# same shape as water's real B2 mode, verified against core.raman_symmetry.
# tensor_form directly in test/.../modes/test.sh) and fabricates EPSIMG data
# for ONLY the single probe-axis pair it names, encoding a KNOWN scale
# factor c=3.0 -- lets the test assert the reconstructed Rxx..Ryz match
# c*T0 EXACTLY, not just "some plausible numbers", for the new single-
# measurement reconstruction math specifically. Builds a minimal
# raman_study/ tree by hand (doesn't need real Stage 1/2 output -- only
# Stage 3's reading side is under test here).
echo -e "\n--- Testing reconstruction from a --reduce-components sidecar (known ground truth) ---"
rm -rf raman_study
python3 - <<'PYEOF'
import json, os
import numpy as np

os.makedirs("raman_study/optical_disp/mode_05_plus_yz", exist_ok=True)
os.makedirs("raman_study/optical_disp/mode_05_minus_yz", exist_ok=True)

T0 = np.array([
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 1.0],
    [0.0, 1.0, 0.0],
])
np.save("reduced_ground_truth_T0.npy", T0)
c_true = 3.0
delta = 0.02
BASE = 5.0
inv_sqrt2 = 0.70710678118654752440
n_probe = np.array([0.0, inv_sqrt2, inv_sqrt2])
quad = n_probe @ T0 @ n_probe  # = 1.0 for this T0

wp, gamma = 5.0, 0.3
w = np.linspace(0.05, 15, 400)
for sign, folder in (("plus", "raman_study/optical_disp/mode_05_plus_yz"),
                      ("minus", "raman_study/optical_disp/mode_05_minus_yz")):
    sign_factor = 1.0 if sign == "plus" else -1.0
    target_static = BASE + sign_factor * delta * c_true * quad
    w0 = (wp**2 / (target_static - 1.0)) ** 0.5
    denom = (w0**2 - w**2)**2 + (gamma*w)**2
    eps2 = wp**2 * gamma * w / denom
    with open(os.path.join(folder, "Sn3O4.EPSIMG"), "w") as f:
        f.write("## Minimum and maximum energy in eV\n")
        f.write(f"##  {w[0]:.6f}  {w[-1]:.6f}\n")
        f.write("## Number of spin components\n")
        f.write("## 1\n")
        for wi, e2i in zip(w, eps2):
            f.write(f"{wi:.6f}  {e2i:.6f}\n")
    with open(os.path.join(folder, "calc.out"), "w") as f:
        f.write("siesta: SCF Convergence by DM criterion\n")
        f.write("SCF cycle converged after 12 iterations\n")
    with open(os.path.join(folder, "calc.fdf"), "w") as f:
        f.write("SystemLabel Sn3O4\n")

with open("raman_study/optical_disp/mode_05_tensor_form.json", "w") as f:
    json.dump({"T0": T0.tolist(), "probe_axis": "yz"}, f)

with open("raman_study/raman_stage2.txt", "w") as f:
    f.write("Displacement      : 0.02 Ang\n")
    f.write("\n# MODE_TABLE -- parsed by stb-ramanAnalysis, do not reorder the first 6 columns\n")
    f.write(f"# {'label':<22}{'mode_index':<12}{'band_index':<12}{'frequency_thz':<16}{'sign':<8}{'axis':<6}{'dir'}\n")
    f.write(f"{'mode_05_plus_yz':<24}{5:<12}{7:<12}{10.0:<16.6f}{'plus':<8}{'yz':<6}raman_study/optical_disp/mode_05_plus_yz\n")
    f.write(f"{'mode_05_minus_yz':<24}{5:<12}{7:<12}{10.0:<16.6f}{'minus':<8}{'yz':<6}raman_study/optical_disp/mode_05_minus_yz\n")
print("fixture ready")
PYEOF
stb-ramanAnalysis --directory raman_study --no-intro > log_reconstruct.txt 2>&1
check_exit_code $? 0
check_contains "Component reduction" log_reconstruct.txt

echo "Testing: reconstructed Rxx..Ryz match the known ground truth c*T0 exactly (c=3.0)"
python3 -c "
import re
import numpy as np
T0 = np.load('reduced_ground_truth_T0.npy')
expected = 3.0 * T0
with open('raman_study/raman_stage3.txt') as f:
    text = f.read()
m = re.search(
    r'Rxx=([\-\d.]+)\s+Ryy=([\-\d.]+)\s+Rzz=([\-\d.]+)\s+Rxy=([\-\d.]+)\s+Rxz=([\-\d.]+)\s+Ryz=([\-\d.]+)',
    text)
assert m, 'could not find the reconstructed tensor line in the report'
rxx, ryy, rzz, rxy, rxz, ryz = (float(v) for v in m.groups())
recovered = np.array([[rxx, rxy, rxz], [rxy, ryy, ryz], [rxz, ryz, rzz]])
assert np.allclose(recovered, expected, atol=0.02), f'recovered {recovered} vs expected {expected}'
print('OK')
" > log_reconstruct_check.txt 2>&1
check_contains "OK" log_reconstruct_check.txt


# --- 4d. --skip-degenerate: a DERIVED mode reuses its representative's tensor exactly ---
# Hand-crafts a minimal raman_study/ tree (doesn't need a real degenerate
# structure -- only Stage 3's reading side is under test here, same
# "build the MODE_TABLE by hand" approach as 4c above): mode 5 gets 6 real
# folders (a KNOWN diagonal Raman tensor, same Lorentz-oscillator trick as
# every other fixture in this file) and mode 6 is written as a bare
# DERIVED row pointing at mode 5 (same format stb-ramanModes --skip-
# degenerate itself writes -- see raman_modes.py). Real CLI-level
# correctness of the degenerate-GROUP DETECTION (an actual doubly-
# degenerate mode collapsing to 1 representative) is covered end to end in
# test/.../modes/test.sh's graphene E2g gate; this test is Stage 3's half:
# given a DERIVED row, does it reuse the exact right numbers.
echo -e "\n--- Testing --skip-degenerate reuse (a DERIVED mode reports its representative's exact tensor) ---"
rm -rf raman_study
python3 - <<'PYEOF'
import os
import numpy as np

for axis in ("x", "y", "z"):
    for sign in ("plus", "minus"):
        os.makedirs(f"raman_study/optical_disp/mode_05_{sign}_{axis}", exist_ok=True)

R_true = {"x": 1.0, "y": 1.0, "z": -2.0}
delta = 0.02
BASE = 5.0
wp, gamma = 5.0, 0.3
w = np.linspace(0.05, 15, 400)

for axis in ("x", "y", "z"):
    for sign in ("plus", "minus"):
        folder = f"raman_study/optical_disp/mode_05_{sign}_{axis}"
        sign_factor = 1.0 if sign == "plus" else -1.0
        target_static = BASE + sign_factor * delta * R_true[axis]
        w0 = (wp**2 / (target_static - 1.0)) ** 0.5
        denom = (w0**2 - w**2)**2 + (gamma*w)**2
        eps2 = wp**2 * gamma * w / denom
        with open(os.path.join(folder, "Sn3O4.EPSIMG"), "w") as f:
            f.write("## Minimum and maximum energy in eV\n")
            f.write(f"##  {w[0]:.6f}  {w[-1]:.6f}\n")
            f.write("## Number of spin components\n")
            f.write("## 1\n")
            for wi, e2i in zip(w, eps2):
                f.write(f"{wi:.6f}  {e2i:.6f}\n")
        with open(os.path.join(folder, "calc.out"), "w") as f:
            f.write("siesta: SCF Convergence by DM criterion\n")
            f.write("SCF cycle converged after 12 iterations\n")
        with open(os.path.join(folder, "calc.fdf"), "w") as f:
            f.write("SystemLabel Sn3O4\n")

with open("raman_study/raman_stage2.txt", "w") as f:
    f.write("Displacement      : 0.02 Ang\n")
    f.write("\n# MODE_TABLE -- parsed by stb-ramanAnalysis, do not reorder the first 6 columns\n")
    f.write(f"# {'label':<22}{'mode_index':<12}{'band_index':<12}{'frequency_thz':<16}"
            f"{'sign':<8}{'axis':<6}{'dir':<24}{'derived_from'}\n")
    for axis in ("x", "y", "z"):
        for sign in ("plus", "minus"):
            label = f"mode_05_{sign}_{axis}"
            f.write(f"{label:<24}{5:<12}{7:<12}{10.0:<16.6f}{sign:<8}{axis:<6}"
                    f"raman_study/optical_disp/{label}  -\n")
    f.write(f"{'mode_06_derived':<24}{6:<12}{8:<12}{10.0:<16.6f}{'DERIVED':<8}{'-':<6}-  5\n")
print("fixture ready")
PYEOF
stb-ramanAnalysis --directory raman_study --no-intro > log_skipdeg_reuse.txt 2>&1
check_exit_code $? 0
check_contains "\[Degenerate partner\]" log_skipdeg_reuse.txt

echo "Testing: mode 6's activity/rho exactly match mode 5's (the representative)"
python3 -c "
import re
with open('raman_study/raman_stage3.txt') as f:
    text = f.read()

def block(mode_num):
    idx = text.index(f'Mode {mode_num} (')
    nxt = text.find('\n  Mode ', idx + 1)
    if nxt == -1:
        nxt = text.find('\n\n', idx)
    return text[idx:nxt]

b5, b6 = block(5), block(6)
assert '[derived]' not in b5, b5
assert '[derived]' in b6, b6
a5 = float(re.search(r'activity~([\-\d.]+)', b5).group(1))
a6 = float(re.search(r'activity~([\-\d.]+)', b6).group(1))
assert a5 == a6, f'activity mismatch: mode 5={a5} vs mode 6={a6}'
r5 = float(re.search(r'rho~([\-\d.]+)', b5).group(1))
r6 = float(re.search(r'rho~([\-\d.]+)', b6).group(1))
assert r5 == r6, f'rho mismatch: mode 5={r5} vs mode 6={r6}'
print('OK')
" > log_skipdeg_reuse_check.txt 2>&1
check_contains "OK" log_skipdeg_reuse_check.txt

echo "Testing: both mode 5 and derived mode 6 count toward Modes analyzed"
check_contains "Modes analyzed      : 2/2" log_skipdeg_reuse.txt


# --- 5. Error cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: --version"
stb-ramanAnalysis --version > log_version.txt 2>&1
check_contains "stb-ramanAnalysis" log_version.txt

echo "Testing: --help documents --linewidth/--file/--temperature/--experimental"
stb-ramanAnalysis --help > log_help.txt 2>&1
check_contains "linewidth" log_help.txt
check_contains "file" log_help.txt
check_contains "temperature" log_help.txt
check_contains "experimental" log_help.txt
check_contains "peak-prominence" log_help.txt


# --- 6. Interactive path (stb-suite, shortcut 4.11.3) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.11.3) ---"

echo "Testing: navigate 4.11.3 -> defaults -> quit"
make_optical_disp
{
  echo "4.11.3"
  echo "raman_study"    # run_dir
  echo ""               # output_filename (default calc.out)
  echo ""               # linewidth (default 10.0)
  echo ""               # temperature (optional, skip)
  echo ""               # experimental_file (optional, skip)
  echo ""               # press enter to continue
  echo "0"              # quit stage submenu
} | stb-suite > log_menu.txt 2>&1
check_contains "Modes analyzed      : 3/3" log_menu.txt
check_success raman_study/raman_spectrum.dat


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
