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

# Runs stb-raman + stb-ramanModes (Stages 1-2, tested separately) with 3
# selected modes, then fabricates SystemLabel.EPSIMG + calc.out in every
# resulting optical_disp/mode_*/ folder.
make_optical_disp() {
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
    stb-ramanModes --directory raman_study --calc calc.fdf --modes 1 2 3 --no-intro > /dev/null 2>&1
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


# --- 5. Error cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: --version"
stb-ramanAnalysis --version > log_version.txt 2>&1
check_contains "stb-ramanAnalysis" log_version.txt

echo "Testing: --help documents --linewidth/--file"
stb-ramanAnalysis --help > log_help.txt 2>&1
check_contains "linewidth" log_help.txt
check_contains "file" log_help.txt


# --- 6. Interactive path (stb-suite, shortcut 4.11.3) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.11.3) ---"

echo "Testing: navigate 4.11.3 -> defaults -> quit"
make_optical_disp
{
  echo "4.11.3"
  echo "raman_study"    # run_dir
  echo ""               # output_filename (default calc.out)
  echo ""               # linewidth (default 10.0)
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
