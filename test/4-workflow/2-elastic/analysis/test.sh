#!/bin/bash

# --- Setup ---
# Smoke test for stb-elasticAnalysis (Elastic Constants Analysis, item 4.2.2)
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
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
    if grep -q "$1" "$2" 2>/dev/null; then
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
echo "--- Starting tester for STB-ElasticAnalysis (item 4.2.2) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
for d in strain_xx_1.00 strain_xx_2.00 strain_xx_m1.00; do
    mkdir -p "$TEST_DIR/$d"
    cp "$FIXTURE_DIR/$d/calc.out" "$TEST_DIR/$d/"
done
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Baseline: existing xx-only fixtures (regression, still works) ---
echo -e "\n--- Testing baseline analysis (xx-only, existing checked-in fixtures) ---"
stb-elasticAnalysis --no-intro > log_baseline.txt 2>&1
check_exit_code $? 0
check_success mechanical_properties.txt
check_contains "C11:    16.0 GPa" mechanical_properties.txt
check_contains "coupling" mechanical_properties.txt


# --- 3. Physics: corrected 2D Poisson's ratio on an anisotropic fixture
#     (C11=200, C22=100, C12=30 GPa by construction) -- own subdirectory,
#     since stb-elasticAnalysis scans ALL 'strain_*' folders in cwd. ---
echo -e "\n--- Testing the corrected 2D Poisson's ratio (anisotropic C11 != C22) ---"
mkdir -p aniso_run
pushd aniso_run > /dev/null
# CONV_EVA3_TO_GPA = 160.21766 -- values below are eV/Ang^3 chosen so the
# fitted GPa slopes come out to exactly C11=200, C22=100, C12=30.
python3 - <<'PYEOF'
import os
CONV = 160.21766
C11, C22, C12 = 200.0, 100.0, 30.0
# Lz = 10.0 Ang makes the --2d conversion factor (Lz * 16.021766) equal the
# plain 3D one (160.21766), so the numbers below map 1:1 to N/m.
OUTCELL = ("outcell: Unit cell vectors (Ang):\n"
           "        5.000000    0.000000    0.000000\n"
           "        0.000000    5.000000    0.000000\n"
           "        0.000000    0.000000   10.000000\n\n")
for pct, tag in [(1.0, "1.00"), (2.0, "2.00"), (-1.0, "m1.00")]:
    eps = pct / 100.0
    sxx, syy = C11 * eps / CONV, C12 * eps / CONV
    folder = f"strain_xx_{tag}"
    os.makedirs(folder, exist_ok=True)
    with open(f"{folder}/calc.out", "w") as f:
        f.write(OUTCELL)
        f.write("siesta: Stress tensor (static) (eV/Ang**3):\n")
        f.write(f"siesta:  {sxx:.8f}   0.00000000   0.00000000\n")
        f.write(f"siesta:   0.00000000  {syy:.8f}   0.00000000\n")
        f.write("siesta:   0.00000000   0.00000000   0.00000000\n")
for pct, tag in [(1.0, "1.00"), (2.0, "2.00"), (-1.0, "m1.00")]:
    eps = pct / 100.0
    syy, sxx = C22 * eps / CONV, C12 * eps / CONV
    folder = f"strain_yy_{tag}"
    os.makedirs(folder, exist_ok=True)
    with open(f"{folder}/calc.out", "w") as f:
        f.write(OUTCELL)
        f.write("siesta: Stress tensor (static) (eV/Ang**3):\n")
        f.write(f"siesta:  {sxx:.8f}   0.00000000   0.00000000\n")
        f.write(f"siesta:   0.00000000  {syy:.8f}   0.00000000\n")
        f.write("siesta:   0.00000000   0.00000000   0.00000000\n")
PYEOF
stb-elasticAnalysis --2d --no-intro > log_aniso.txt 2>&1
check_exit_code $? 0
check_contains "C11.*200.00" log_aniso.txt
check_contains "C22.*100.00" log_aniso.txt
check_contains "C12.*3.00\|C12.*30.00" log_aniso.txt
# Correct values (derived by hand + cross-checked via the Maxwell-Betti
# reciprocal relation nu_yx/Ex == nu_xy/Ey): nu_yx=C12/C22, nu_xy=C12/C11.
# With this fixture's C11/C22/C12 ratios, nu_yx=0.300, nu_xy=0.150 --
# the OLD (buggy) code would have printed 0.150 as "v_21".
check_contains "Poisson's Ratio (v_yx).*0.300" log_aniso.txt
check_contains "Poisson's Ratio (v_xy).*0.150" log_aniso.txt
popd > /dev/null


# --- 4. Physics: 'x'/'xx' folder-naming merge (was silently discarding
#     one series before) -- own subdirectory. ---
echo -e "\n--- Testing 'x' and 'xx' folders merge into the same series ---"
mkdir -p merge_run
pushd merge_run > /dev/null
python3 - <<'PYEOF'
import os
CONV = 160.21766
C11 = 200.0
for folder, pct in [("strain_xx_1.00", 1.0), ("strain_x_2.00", 2.0)]:
    eps = pct / 100.0
    sxx = C11 * eps / CONV
    os.makedirs(folder, exist_ok=True)
    with open(f"{folder}/calc.out", "w") as f:
        f.write("siesta: Stress tensor (static) (eV/Ang**3):\n")
        f.write(f"siesta:  {sxx:.8f}   0.00000000   0.00000000\n")
        f.write("siesta:   0.00000000   0.00000000   0.00000000\n")
        f.write("siesta:   0.00000000   0.00000000   0.00000000\n")
PYEOF
stb-elasticAnalysis --no-intro > log_merge.txt 2>&1
check_exit_code $? 0
check_contains "Loaded 2 calculations from 2 folders" log_merge.txt
# If the two folders hadn't merged, only 1 point would survive for 'xx' and
# calculate_slope's <2-points guard would report C11=0.0, not 200.0.
check_contains "C11:   200.0 GPa" log_merge.txt
popd > /dev/null


# --- 5. Physics: full 21-constant (P1/triclinic) tensor recovery. All 6
#     canonical directions with a FULL 3x3 stress response each (not just
#     the "obvious" component) -- verifies the normal-shear coupling
#     constants (C14-C16/C24-C26/C34-C36/C45/C46/C56) are correctly
#     recovered from data already collected by the standard 6-direction
#     sweep, no new strain patterns needed. Reference C matrix + expected
#     values cross-checked by hand and by an independent numpy script
#     before implementing (max error 3e-14 recovering this exact matrix). ---
echo -e "\n--- Testing full 21-constant tensor recovery (P1/triclinic-like) ---"
mkdir -p full_tensor_run
pushd full_tensor_run > /dev/null
python3 - <<'PYEOF'
import os
import numpy as np

CONV = 160.21766  # eV/Ang^3 -> GPa

C_true = np.zeros((6, 6))
C_true[0,0], C_true[1,1], C_true[2,2] = 200.0, 150.0, 180.0
C_true[0,1] = C_true[1,0] = 40.0
C_true[0,2] = C_true[2,0] = 35.0
C_true[1,2] = C_true[2,1] = 30.0
C_true[3,3], C_true[4,4], C_true[5,5] = 60.0, 55.0, 65.0
C_true[0,3] = C_true[3,0] = 5.0    # C14
C_true[0,4] = C_true[4,0] = 8.0    # C15
C_true[0,5] = C_true[5,0] = -3.0   # C16
C_true[1,3] = C_true[3,1] = 2.0    # C24
C_true[1,5] = C_true[5,1] = 4.0    # C26
C_true[3,4] = C_true[4,3] = 1.5    # C45
C_true[4,5] = C_true[5,4] = -2.0   # C56

def strain_tensor(direction, delta):
    eps = np.zeros((3, 3))
    if direction == 'xx': eps[0,0] = delta
    elif direction == 'yy': eps[1,1] = delta
    elif direction == 'zz': eps[2,2] = delta
    elif direction == 'yz': eps[1,2] = eps[2,1] = delta
    elif direction == 'xz': eps[0,2] = eps[2,0] = delta
    elif direction == 'xy': eps[0,1] = eps[1,0] = delta
    return eps

def voigt_strain(eps_t):
    return np.array([eps_t[0,0], eps_t[1,1], eps_t[2,2],
                      2*eps_t[1,2], 2*eps_t[0,2], 2*eps_t[0,1]])

def voigt_to_tensor(sig_v):
    s = np.zeros((3, 3))
    s[0,0], s[1,1], s[2,2], s[1,2], s[0,2], s[0,1] = sig_v
    s[2,1], s[2,0], s[1,0] = s[1,2], s[0,2], s[0,1]
    return s

for direction in ['xx', 'yy', 'zz', 'yz', 'xz', 'xy']:
    for pct, tag in [(-2.0, "m2.00"), (-0.67, "m0.67"), (0.67, "0.67"), (2.0, "2.00")]:
        delta = pct / 100.0
        sigma_v_gpa = C_true @ voigt_strain(strain_tensor(direction, delta))
        sigma_t = voigt_to_tensor(sigma_v_gpa / CONV)
        folder = f"strain_{direction}_{tag}"
        os.makedirs(folder, exist_ok=True)
        with open(f"{folder}/calc.out", "w") as f:
            f.write("siesta: Stress tensor (static) (eV/Ang**3):\n")
            for row in range(3):
                f.write("siesta: " + "   ".join(f"{sigma_t[row, c]:.10f}" for c in range(3)) + "\n")
PYEOF
stb-elasticAnalysis --no-intro > log_full_tensor.txt 2>&1
check_exit_code $? 0
check_contains "C11:   200.0 GPa" log_full_tensor.txt
check_contains "C22:   150.0 GPa" log_full_tensor.txt
check_contains "C33:   180.0 GPa" log_full_tensor.txt
check_contains "COUPLING CONSTANTS" log_full_tensor.txt
check_contains "C14:     5.0 GPa" log_full_tensor.txt
check_contains "C15:     8.0 GPa" log_full_tensor.txt
check_contains "C16:    -3.0 GPa" log_full_tensor.txt
check_contains "C24:     2.0 GPa" log_full_tensor.txt
check_contains "C26:     4.0 GPa" log_full_tensor.txt
check_contains "C45:     1.5 GPa" log_full_tensor.txt
check_contains "C56:    -2.0 GPa" log_full_tensor.txt
check_contains "VERDICT:" log_full_tensor.txt
check_contains "STABLE" log_full_tensor.txt
popd > /dev/null


# --- 6. Physics: high-symmetry (cubic-like) fixture must NOT show a
#     coupling-constants section (all near 0), and a partial direction set
#     must report exactly which directions/columns are missing. ---
echo -e "\n--- Testing high-symmetry (no coupling) vs. partial-data (missing directions) ---"
mkdir -p cubic_run
pushd cubic_run > /dev/null
python3 - <<'PYEOF'
import os
import numpy as np
CONV = 160.21766
C_true = np.zeros((6, 6))
C_true[0,0] = C_true[1,1] = C_true[2,2] = 166.0
C_true[0,1] = C_true[1,0] = C_true[0,2] = C_true[2,0] = C_true[1,2] = C_true[2,1] = 64.0
C_true[3,3] = C_true[4,4] = C_true[5,5] = 80.0  # cubic: zero coupling

def strain_tensor(direction, delta):
    eps = np.zeros((3, 3))
    if direction == 'xx': eps[0,0] = delta
    elif direction == 'yy': eps[1,1] = delta
    elif direction == 'zz': eps[2,2] = delta
    elif direction == 'yz': eps[1,2] = eps[2,1] = delta
    elif direction == 'xz': eps[0,2] = eps[2,0] = delta
    elif direction == 'xy': eps[0,1] = eps[1,0] = delta
    return eps

def voigt_strain(eps_t):
    return np.array([eps_t[0,0], eps_t[1,1], eps_t[2,2],
                      2*eps_t[1,2], 2*eps_t[0,2], 2*eps_t[0,1]])

def voigt_to_tensor(sig_v):
    s = np.zeros((3, 3))
    s[0,0], s[1,1], s[2,2], s[1,2], s[0,2], s[0,1] = sig_v
    s[2,1], s[2,0], s[1,0] = s[1,2], s[0,2], s[0,1]
    return s

for direction in ['xx', 'yy', 'zz', 'yz', 'xz', 'xy']:
    for pct, tag in [(-2.0, "m2.00"), (-0.67, "m0.67"), (0.67, "0.67"), (2.0, "2.00")]:
        delta = pct / 100.0
        sigma_v_gpa = C_true @ voigt_strain(strain_tensor(direction, delta))
        sigma_t = voigt_to_tensor(sigma_v_gpa / CONV)
        folder = f"strain_{direction}_{tag}"
        os.makedirs(folder, exist_ok=True)
        with open(f"{folder}/calc.out", "w") as f:
            f.write("siesta: Stress tensor (static) (eV/Ang**3):\n")
            for row in range(3):
                f.write("siesta: " + "   ".join(f"{sigma_t[row, c]:.10f}" for c in range(3)) + "\n")
PYEOF

echo "Testing: all 6 directions present, cubic symmetry -> no coupling section"
stb-elasticAnalysis --no-intro > log_cubic_full.txt 2>&1
check_exit_code $? 0
check_contains "C11:   166.0 GPa" log_cubic_full.txt
if grep -q "COUPLING CONSTANTS" log_cubic_full.txt 2>/dev/null; then
    echo -e "   -> ${RED}Failed:${NC} unexpected coupling-constants section for a cubic-symmetric fixture"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no coupling-constants section for near-zero coupling"
    PASS=$((PASS+1))
fi

echo "Testing: only xx/yy/zz present -> missing-direction note names yz, xz, xy"
mkdir -p partial_run
pushd partial_run > /dev/null
cp -r ../strain_xx_* ../strain_yy_* ../strain_zz_* .
stb-elasticAnalysis --no-intro > log_partial.txt 2>&1
check_exit_code $? 0
check_contains "No data for direction(s) yz, xz, xy" log_partial.txt
check_contains "0 by missing data, not by symmetry" log_partial.txt
popd > /dev/null
popd > /dev/null


# --- 7. Physics: symmetry-based DFT-cost reduction round-trip. Runs the
#     REAL stb-elasticInputs (not a hand-built fixture) so it auto-reduces
#     'all' to just xx+xy for a cubic structure, then confirms
#     stb-elasticAnalysis reconstructs the other 4 columns exactly via the
#     real point-group operation (not a guessed permutation table). ---
echo -e "\n--- Testing the full prep+analysis symmetry-reduction round-trip (cubic) ---"
mkdir -p reduction_run
pushd reduction_run > /dev/null
cp "$FIXTURE_DIR/../prep/si_cubic.fdf" .
stb-elasticInputs --file si_cubic.fdf --dirs all --no-intro > log_prep.txt 2>&1
check_contains "Modes: \['xx', 'xy'\]" log_prep.txt
check_success reference_structure.fdf
python3 - <<'PYEOF'
import os, re
import numpy as np

CONV = 160.21766
C_true = np.zeros((6, 6))
C_true[0,0] = C_true[1,1] = C_true[2,2] = 166.0
C_true[0,1] = C_true[1,0] = C_true[0,2] = C_true[2,0] = C_true[1,2] = C_true[2,1] = 64.0
C_true[3,3] = C_true[4,4] = C_true[5,5] = 80.0

def strain_tensor(mode, delta):
    eps = np.zeros((3, 3))
    if mode == 'xx': eps[0,0] = delta
    elif mode == 'xy': eps[0,1] = eps[1,0] = delta
    return eps

def voigt_strain(t):
    return np.array([t[0,0], t[1,1], t[2,2], 2*t[1,2], 2*t[0,2], 2*t[0,1]])

def voigt_to_tensor(v):
    s = np.zeros((3, 3))
    s[0,0], s[1,1], s[2,2], s[1,2], s[0,2], s[0,1] = v
    s[2,1], s[2,0], s[1,0] = s[1,2], s[0,2], s[0,1]
    return s

for folder in os.listdir('.'):
    if not folder.startswith('strain_'):
        continue
    m = re.match(r"strain_(xx|xy)_(m?)(\d+\.\d+)", folder)
    mode, sign, val = m.group(1), m.group(2), float(m.group(3))
    pct = -val if sign == 'm' else val
    delta = pct / 100.0
    sigma_v = C_true @ voigt_strain(strain_tensor(mode, delta))
    sigma_t = voigt_to_tensor(sigma_v / CONV)
    with open(f"{folder}/calc.out", "w") as f:
        f.write("siesta: Stress tensor (static) (eV/Ang**3):\n")
        for row in range(3):
            f.write("siesta: " + "   ".join(f"{sigma_t[row, c]:.10f}" for c in range(3)) + "\n")
PYEOF
stb-elasticAnalysis --no-intro > log_reduction.txt 2>&1
check_exit_code $? 0
check_contains "Detected point group m-3m" log_reduction.txt
check_contains "C11:   166.0 GPa" log_reduction.txt
check_contains "C22:   166.0 GPa" log_reduction.txt
check_contains "C33:   166.0 GPa" log_reduction.txt
check_contains "C44:    80.0 GPa" log_reduction.txt
check_contains "C55:    80.0 GPa" log_reduction.txt
check_contains "C66:    80.0 GPa" log_reduction.txt
check_contains "C12:    64.0 GPa" log_reduction.txt
check_contains "C13:    64.0 GPa" log_reduction.txt
check_contains "C23:    64.0 GPa" log_reduction.txt
check_contains "Column(s) yy, zz, yz, xz were not computed directly" log_reduction.txt
popd > /dev/null


# --- 8. Physics: pooling/diagnostic across an explicit (non-reduced) full
#     6-direction run, with a different measured value per equivalent
#     direction (synthetic DFT noise) -- exercises the sign-aware stress
#     transform for ALL 3 pairings within the shear group {yz,xz,xy}, which
#     is exactly where a real bug was caught during development: naively
#     negating R (instead of the final transformed stress) for a sign-flip
#     mapping left C44/C55/C66 silently wrong (1/3 of the true value) while
#     the normal-strain group's C11/C22/C33 pooling looked fine -- this
#     case must never regress. ---
echo -e "\n--- Testing pooled fitting across explicit equivalent directions (cubic, with noise) ---"
mkdir -p pooling_run
pushd pooling_run > /dev/null
cp "$FIXTURE_DIR/../prep/si_cubic.fdf" .
stb-elasticInputs --file si_cubic.fdf --dirs xx yy zz xy xz yz --no-intro > log_prep.txt 2>&1
check_contains "Modes: \['xx', 'yy', 'zz', 'xy', 'xz', 'yz'\]" log_prep.txt
python3 - <<'PYEOF'
import os, re
import numpy as np

CONV = 160.21766
C_by_mode = {'xx': 165.8, 'yy': 166.3, 'zz': 166.0}
C_shear = 80.0

def voigt_to_tensor(v):
    s = np.zeros((3, 3))
    s[0,0], s[1,1], s[2,2], s[1,2], s[0,2], s[0,1] = v
    s[2,1], s[2,0], s[1,0] = s[1,2], s[0,2], s[0,1]
    return s

for folder in os.listdir('.'):
    if not folder.startswith('strain_'):
        continue
    m = re.match(r"strain_(xx|yy|zz|xy|xz|yz)_(m?)(\d+\.\d+)", folder)
    mode, sign, val = m.group(1), m.group(2), float(m.group(3))
    pct = -val if sign == 'm' else val
    delta = pct / 100.0
    sigma_v = np.zeros(6)
    if mode in C_by_mode:
        diag = C_by_mode[mode]
        for k, mm in enumerate(['xx', 'yy', 'zz']):
            sigma_v[k] = diag * delta if mm == mode else 64.0 * delta
    else:
        col = {'yz': 3, 'xz': 4, 'xy': 5}[mode]
        sigma_v[col] = C_shear * (2 * delta)
    sigma_t = voigt_to_tensor(sigma_v / CONV)
    with open(f"{folder}/calc.out", "w") as f:
        f.write("siesta: Stress tensor (static) (eV/Ang**3):\n")
        for row in range(3):
            f.write("siesta: " + "   ".join(f"{sigma_t[row, c]:.10f}" for c in range(3)) + "\n")
PYEOF
stb-elasticAnalysis --no-intro > log_pooling.txt 2>&1
check_exit_code $? 0
check_contains "\['xx', 'yy', 'zz'\] are symmetry-equivalent" log_pooling.txt
check_contains "xx: 165.80 GPa (independent, pre-pooling)" log_pooling.txt
check_contains "yy: 166.30 GPa (independent, pre-pooling)" log_pooling.txt
check_contains "zz: 166.00 GPa (independent, pre-pooling)" log_pooling.txt
check_contains "xx combined (12 pooled points): 166.03 GPa" log_pooling.txt
# The regression this guards: shear-group pooling must recover the true
# C44=C55=C66=80.0 GPa, not 26.7 GPa (=80/3, the sign-flip-bug symptom).
check_contains "C44:    80.0 GPa" log_pooling.txt
check_contains "C55:    80.0 GPa" log_pooling.txt
check_contains "C66:    80.0 GPa" log_pooling.txt
popd > /dev/null


# --- 9. Physics: --symmetry-method full -- the general tensor-invariance
#     method. Unlike 'basic', it constrains the FULL elastic tensor's
#     symmetry (not a pairwise strain-tensor match), so it also catches
#     hexagonal/trigonal reductions that 'basic' structurally cannot (see
#     core/symmetry.py's module docstring). ---
echo -e "\n--- Testing --symmetry-method full (general tensor-invariance reconstruction) ---"

echo "Testing: cubic fixture, full method agrees exactly with the basic-method result"
mkdir -p full_cubic_run
pushd full_cubic_run > /dev/null
cp "$FIXTURE_DIR/../prep/si_cubic.fdf" .
stb-elasticInputs --file si_cubic.fdf --dirs all --symmetry-method full --no-intro > log_prep.txt 2>&1
check_contains "Modes: \['xx', 'xy'\]" log_prep.txt
python3 - <<'PYEOF'
import os, re
import numpy as np

CONV = 160.21766
C_true = np.zeros((6, 6))
C_true[0,0] = C_true[1,1] = C_true[2,2] = 166.0
C_true[0,1] = C_true[1,0] = C_true[0,2] = C_true[2,0] = C_true[1,2] = C_true[2,1] = 64.0
C_true[3,3] = C_true[4,4] = C_true[5,5] = 80.0

def strain_tensor(mode, delta):
    eps = np.zeros((3, 3))
    if mode == 'xx': eps[0,0] = delta
    elif mode == 'xy': eps[0,1] = eps[1,0] = delta
    return eps

def voigt_strain(t):
    return np.array([t[0,0], t[1,1], t[2,2], 2*t[1,2], 2*t[0,2], 2*t[0,1]])

def voigt_to_tensor(v):
    s = np.zeros((3, 3))
    s[0,0], s[1,1], s[2,2], s[1,2], s[0,2], s[0,1] = v
    s[2,1], s[2,0], s[1,0] = s[1,2], s[0,2], s[0,1]
    return s

for folder in os.listdir('.'):
    if not folder.startswith('strain_'):
        continue
    m = re.match(r"strain_(xx|xy)_(m?)(\d+\.\d+)", folder)
    mode, sign, val = m.group(1), m.group(2), float(m.group(3))
    pct = -val if sign == 'm' else val
    delta = pct / 100.0
    sigma_v = C_true @ voigt_strain(strain_tensor(mode, delta))
    sigma_t = voigt_to_tensor(sigma_v / CONV)
    with open(f"{folder}/calc.out", "w") as f:
        f.write("siesta: Stress tensor (static) (eV/Ang**3):\n")
        for row in range(3):
            f.write("siesta: " + "   ".join(f"{sigma_t[row, c]:.10f}" for c in range(3)) + "\n")
PYEOF
stb-elasticAnalysis --symmetry-method full --no-intro > log_full_cubic.txt 2>&1
check_exit_code $? 0
check_contains "Point group m-3m has 3 independent elastic constant" log_full_cubic.txt
check_contains "C11:   166.0 GPa" log_full_cubic.txt
check_contains "C22:   166.0 GPa" log_full_cubic.txt
check_contains "C44:    80.0 GPa" log_full_cubic.txt
check_contains "C12:    64.0 GPa" log_full_cubic.txt
popd > /dev/null

echo "Testing: 2D hexagonal fixture -- the real target case (basic method finds no equivalence at all)"
mkdir -p full_hex_run
pushd full_hex_run > /dev/null
cp "$FIXTURE_DIR/../prep/structure.fdf" .
stb-elasticInputs --file structure.fdf --dirs all --symmetry-method full --no-intro > log_prep.txt 2>&1
check_contains "Point group -6m2 has 5 independent elastic constant" log_prep.txt
check_contains "Modes: \['xx'\]" log_prep.txt
python3 - <<'PYEOF'
import os, re
import numpy as np

CONV = 160.21766
C11, C12 = 200.0, 50.0
C13, C33, C44 = 40.0, 300.0, 60.0
C66 = (C11 - C12) / 2
C_true = np.zeros((6, 6))
C_true[0,0] = C_true[1,1] = C11
C_true[0,1] = C_true[1,0] = C12
C_true[0,2] = C_true[2,0] = C_true[1,2] = C_true[2,1] = C13
C_true[2,2] = C33
C_true[3,3] = C_true[4,4] = C44
C_true[5,5] = C66

def strain_tensor(mode, delta):
    eps = np.zeros((3, 3))
    if mode == 'xx': eps[0,0] = delta
    return eps

def voigt_strain(t):
    return np.array([t[0,0], t[1,1], t[2,2], 2*t[1,2], 2*t[0,2], 2*t[0,1]])

def voigt_to_tensor(v):
    s = np.zeros((3, 3))
    s[0,0], s[1,1], s[2,2], s[1,2], s[0,2], s[0,1] = v
    s[2,1], s[2,0], s[1,0] = s[1,2], s[0,2], s[0,1]
    return s

for folder in os.listdir('.'):
    if not folder.startswith('strain_'):
        continue
    m = re.match(r"strain_(xx)_(m?)(\d+\.\d+)", folder)
    mode, sign, val = m.group(1), m.group(2), float(m.group(3))
    pct = -val if sign == 'm' else val
    delta = pct / 100.0
    sigma_v = C_true @ voigt_strain(strain_tensor(mode, delta))
    sigma_t = voigt_to_tensor(sigma_v / CONV)
    with open(f"{folder}/calc.out", "w") as f:
        f.write("siesta: Stress tensor (static) (eV/Ang**3):\n")
        for row in range(3):
            f.write("siesta: " + "   ".join(f"{sigma_t[row, c]:.10f}" for c in range(3)) + "\n")
PYEOF
# Only 'xx' was ever run, yet the 2D report must show C22=C11 and
# C66=(C11-C12)/2 exactly -- the hexagonal reduction 'basic' cannot see.
stb-elasticAnalysis --2d --symmetry-method full --no-intro > log_full_hex.txt 2>&1
check_exit_code $? 0
check_contains "Point group -6m2 has 5 independent elastic constant" log_full_hex.txt
check_contains "fit used 1 of 6 canonical direction(s): xx" log_full_hex.txt
check_contains "C11.*20.00" log_full_hex.txt
check_contains "C22.*20.00" log_full_hex.txt
check_contains "C12.*5.00" log_full_hex.txt
check_contains "C66.*7.50" log_full_hex.txt
check_contains "Direction(s) yy, zz, yz, xz, xy were never run directly" log_full_hex.txt
popd > /dev/null

echo "Testing: triclinic fixture, full method requires --reference-structure (else a clear error)"
mkdir -p full_noref_run
pushd full_noref_run > /dev/null
mkdir -p strain_xx_2.00
cp "$FIXTURE_DIR/strain_xx_2.00/calc.out" strain_xx_2.00/ 2>/dev/null || true
stb-elasticAnalysis --symmetry-method full --no-intro > log_noref.txt 2>&1
check_exit_code $? 1
check_contains "requires a readable reference_structure.fdf" log_noref.txt
popd > /dev/null


# --- 10. Physics: --method energy -- energy-strain (parabolic fit) method.
#     Fits total energy vs. strain (E = E0 + (V0/2) gamma^T C gamma)
#     instead of stress vs. strain -- a physically independent cross-check.
#     Runs the REAL stb-elasticInputs --method energy so its own pattern
#     selection (pure + combined Voigt patterns) is exercised end to end,
#     not hand-picked. ---
echo -e "\n--- Testing --method energy (energy-strain parabolic fit) ---"

echo "Testing: cubic fixture -- 3 patterns (2 pure + 1 combined), exact recovery"
mkdir -p energy_cubic_run
pushd energy_cubic_run > /dev/null
cp "$FIXTURE_DIR/../prep/si_cubic.fdf" .
stb-elasticInputs --file si_cubic.fdf --method energy --no-intro > log_prep.txt 2>&1
check_contains "Modes: \['xx', 'yz', 'xx+yy'\]" log_prep.txt
python3 - <<'PYEOF'
import os, re
import numpy as np
from stb.core import structure_io, symmetry

CONV = 160.21766
structure = structure_io.read_fdf("reference_structure.fdf")
pmg = structure_io.to_pymatgen(structure)
V0 = pmg.volume

C_true = np.zeros((6, 6))
C_true[0,0] = C_true[1,1] = C_true[2,2] = 166.0
C_true[0,1] = C_true[1,0] = C_true[0,2] = C_true[2,0] = C_true[1,2] = C_true[2,1] = 64.0
C_true[3,3] = C_true[4,4] = C_true[5,5] = 80.0
C_raw = C_true / CONV

for folder in os.listdir('.'):
    if not folder.startswith('strain_'):
        continue
    m = re.match(r"strain_([a-z+]+)_(m?)(\d+\.\d+)", folder)
    pattern, sign, val = m.group(1), m.group(2), float(m.group(3))
    pct = -val if sign == 'm' else val
    delta = pct / 100.0
    eps = sum((symmetry.strain_tensor(p, delta) for p in pattern.split('+')), np.zeros((3, 3)))
    gamma = symmetry.voigt_strain_engineering(eps)
    dE = 0.5 * V0 * (gamma @ C_raw @ gamma)
    with open(f"{folder}/calc.out", "w") as f:
        f.write(f"siesta: FreeEng = {dE:.10f}\n")
PYEOF
stb-elasticAnalysis --method energy --no-intro > log_energy_cubic.txt 2>&1
check_exit_code $? 0
check_contains "Point group m-3m has 3 independent elastic constant" log_energy_cubic.txt
check_contains "C11:   166.0 GPa" log_energy_cubic.txt
check_contains "C22:   166.0 GPa" log_energy_cubic.txt
check_contains "C44:    80.0 GPa" log_energy_cubic.txt
check_contains "C12:    64.0 GPa" log_energy_cubic.txt
popd > /dev/null

echo "Testing: 2D hexagonal fixture -- only 2 in-plane patterns needed, exact recovery"
mkdir -p energy_hex_run
pushd energy_hex_run > /dev/null
cp "$FIXTURE_DIR/../prep/structure.fdf" .
stb-elasticInputs --file structure.fdf --method energy --no-intro > log_prep.txt 2>&1
check_contains "Modes: \['xx', 'xy'\]" log_prep.txt
python3 - <<'PYEOF'
import os, re
import numpy as np
from stb.core import structure_io, symmetry

CONV = 16.021766  # eV/A^2 -> N/m
structure = structure_io.read_fdf("reference_structure.fdf")
pmg = structure_io.to_pymatgen(structure)
lat = pmg.lattice.matrix
Area0 = np.linalg.norm(np.cross(lat[0], lat[1]))

C11, C12 = 200.0, 50.0
C13, C33, C44 = 40.0, 300.0, 60.0
C66 = (C11 - C12) / 2
C_true = np.zeros((6, 6))
C_true[0,0] = C_true[1,1] = C11
C_true[0,1] = C_true[1,0] = C12
C_true[0,2] = C_true[2,0] = C_true[1,2] = C_true[2,1] = C13
C_true[2,2] = C33
C_true[3,3] = C_true[4,4] = C44
C_true[5,5] = C66
C_raw = C_true / CONV

for folder in os.listdir('.'):
    if not folder.startswith('strain_'):
        continue
    m = re.match(r"strain_([a-z+]+)_(m?)(\d+\.\d+)", folder)
    pattern, sign, val = m.group(1), m.group(2), float(m.group(3))
    pct = -val if sign == 'm' else val
    delta = pct / 100.0
    eps = sum((symmetry.strain_tensor(p, delta) for p in pattern.split('+')), np.zeros((3, 3)))
    gamma = symmetry.voigt_strain_engineering(eps)
    dE = 0.5 * Area0 * (gamma @ C_raw @ gamma)
    with open(f"{folder}/calc.out", "w") as f:
        f.write(f"siesta: FreeEng = {dE:.10f}\n")
PYEOF
stb-elasticAnalysis --2d --method energy --no-intro > log_energy_hex.txt 2>&1
check_exit_code $? 0
check_contains "Point group -6m2 has 5 independent elastic constant" log_energy_hex.txt
check_contains "fit used 2 strain pattern(s): \['xx', 'xy'\]" log_energy_hex.txt
check_contains "C11.*200.00" log_energy_hex.txt
check_contains "C22.*200.00" log_energy_hex.txt
check_contains "C12.*50.00" log_energy_hex.txt
check_contains "C66.*75.00" log_energy_hex.txt
popd > /dev/null

echo "Testing: diagnostic -- non-negligible linear energy term warns about a poorly relaxed reference"
mkdir -p energy_linear_run
pushd energy_linear_run > /dev/null
cp "$FIXTURE_DIR/../prep/si_cubic.fdf" .
stb-elasticInputs --file si_cubic.fdf --method energy --no-intro > log_prep.txt 2>&1
python3 - <<'PYEOF'
import os, re
import numpy as np
from stb.core import structure_io, symmetry

CONV = 160.21766
structure = structure_io.read_fdf("reference_structure.fdf")
pmg = structure_io.to_pymatgen(structure)
V0 = pmg.volume

C_true = np.zeros((6, 6))
C_true[0,0] = C_true[1,1] = C_true[2,2] = 166.0
C_true[0,1] = C_true[1,0] = C_true[0,2] = C_true[2,0] = C_true[1,2] = C_true[2,1] = 64.0
C_true[3,3] = C_true[4,4] = C_true[5,5] = 80.0
C_raw = C_true / CONV

for folder in os.listdir('.'):
    if not folder.startswith('strain_'):
        continue
    m = re.match(r"strain_([a-z+]+)_(m?)(\d+\.\d+)", folder)
    pattern, sign, val = m.group(1), m.group(2), float(m.group(3))
    pct = -val if sign == 'm' else val
    delta = pct / 100.0
    eps = sum((symmetry.strain_tensor(p, delta) for p in pattern.split('+')), np.zeros((3, 3)))
    gamma = symmetry.voigt_strain_engineering(eps)
    dE = 0.5 * V0 * (gamma @ C_raw @ gamma)
    if pattern == 'xx':
        dE += 0.5 * delta  # deliberate un-relaxed-reference residual force
    with open(f"{folder}/calc.out", "w") as f:
        f.write(f"siesta: FreeEng = {dE:.10f}\n")
PYEOF
stb-elasticAnalysis --method energy --no-intro > log_energy_linear.txt 2>&1
check_contains "linear energy term.*not.*negligible" log_energy_linear.txt
popd > /dev/null

echo "Testing: --method energy requires --reference-structure (else a clear error)"
mkdir -p energy_noref_run
pushd energy_noref_run > /dev/null
mkdir -p strain_xx_2.00
cp "$FIXTURE_DIR/strain_xx_2.00/calc.out" strain_xx_2.00/ 2>/dev/null || true
stb-elasticAnalysis --method energy --no-intro > log_noref.txt 2>&1
check_exit_code $? 1
check_contains "requires a readable reference_structure.fdf" log_noref.txt
popd > /dev/null


# --- 11. Physics: [3] NUMERICAL QUALITY DIAGNOSTICS -- 3 mandatory,
#     zero-new-DFT-calculation checks whenever --method stress (the
#     default) is used: eggbox cross-check (stress-linear-fit vs.
#     energy-parabolic-fit of the SAME calc.out files -- see below), per
#     -direction fit quality (R^2 + residual stress at zero strain, reusing
#     scipy.stats.linregress's own r_value/intercept, previously discarded
#     by calculate_slope), and tensor symmetry (C vs C^T, only meaningful
#     for --symmetry-method basic, where C[i,j]/C[j,i] can come from two
#     independently-run directions). ---
echo -e "\n--- Testing the mandatory eggbox cross-check (--method stress) ---"

echo "Testing: consistent stress+energy data -> section [3] shows 0% diff, no warning"
mkdir -p eggbox_consistent_run
pushd eggbox_consistent_run > /dev/null
cp "$FIXTURE_DIR/../prep/si_cubic.fdf" .
stb-elasticInputs --file si_cubic.fdf --dirs all --symmetry-method full --no-intro > log_prep.txt 2>&1
python3 - <<'PYEOF'
import os, re
import numpy as np
from stb.core import structure_io, symmetry

CONV = 160.21766
structure = structure_io.read_fdf("reference_structure.fdf")
pmg = structure_io.to_pymatgen(structure)
V0 = pmg.volume

C_true = np.zeros((6, 6))
C_true[0,0] = C_true[1,1] = C_true[2,2] = 166.0
C_true[0,1] = C_true[1,0] = C_true[0,2] = C_true[2,0] = C_true[1,2] = C_true[2,1] = 64.0
C_true[3,3] = C_true[4,4] = C_true[5,5] = 80.0
C_raw = C_true / CONV

def voigt_to_tensor(v):
    s = np.zeros((3, 3))
    s[0,0], s[1,1], s[2,2], s[1,2], s[0,2], s[0,1] = v
    s[2,1], s[2,0], s[1,0] = s[1,2], s[0,2], s[0,1]
    return s

for folder in os.listdir('.'):
    if not folder.startswith('strain_'):
        continue
    m = re.match(r"strain_([a-z+]+)_(m?)(\d+\.\d+)", folder)
    pattern, sign, val = m.group(1), m.group(2), float(m.group(3))
    pct = -val if sign == 'm' else val
    delta = pct / 100.0
    eps = symmetry.strain_tensor(pattern, delta)
    gamma = symmetry.voigt_strain_engineering(eps)
    sigma_t = voigt_to_tensor(C_raw @ gamma)
    dE = 0.5 * V0 * (gamma @ C_raw @ gamma)
    with open(f"{folder}/calc.out", "w") as f:
        f.write("siesta: Stress tensor (static) (eV/Ang**3):\n")
        for row in range(3):
            f.write("siesta: " + "   ".join(f"{sigma_t[row, c]:.10f}" for c in range(3)) + "\n")
        f.write(f"siesta: FreeEng = {dE:.10f}\n")
PYEOF
stb-elasticAnalysis --no-intro > log_consistent.txt 2>&1
check_exit_code $? 0
check_contains "\[3\] NUMERICAL QUALITY DIAGNOSTICS" log_consistent.txt
check_contains "Eggbox cross-check (stress vs. energy, same data)" log_consistent.txt
check_contains "xx: stress =   166.00 GPa | energy =   166.00 GPa | diff =  0.00%" log_consistent.txt
check_contains "xy: stress =    80.00 GPa | energy =    80.00 GPa | diff =  0.00%" log_consistent.txt
check_contains "All checked direction(s) agree within tolerance" log_consistent.txt
popd > /dev/null

echo "Testing: contaminated energy on one direction only -> flagged, other direction unaffected"
mkdir -p eggbox_contaminated_run
pushd eggbox_contaminated_run > /dev/null
cp "$FIXTURE_DIR/../prep/si_cubic.fdf" .
stb-elasticInputs --file si_cubic.fdf --dirs all --symmetry-method full --no-intro > log_prep.txt 2>&1
python3 - <<'PYEOF'
import os, re
import numpy as np
from stb.core import structure_io, symmetry

CONV = 160.21766
structure = structure_io.read_fdf("reference_structure.fdf")
pmg = structure_io.to_pymatgen(structure)
V0 = pmg.volume

C_true = np.zeros((6, 6))
C_true[0,0] = C_true[1,1] = C_true[2,2] = 166.0
C_true[0,1] = C_true[1,0] = C_true[0,2] = C_true[2,0] = C_true[1,2] = C_true[2,1] = 64.0
C_true[3,3] = C_true[4,4] = C_true[5,5] = 80.0
C_raw = C_true / CONV

def voigt_to_tensor(v):
    s = np.zeros((3, 3))
    s[0,0], s[1,1], s[2,2], s[1,2], s[0,2], s[0,1] = v
    s[2,1], s[2,0], s[1,0] = s[1,2], s[0,2], s[0,1]
    return s

for folder in os.listdir('.'):
    if not folder.startswith('strain_'):
        continue
    m = re.match(r"strain_([a-z+]+)_(m?)(\d+\.\d+)", folder)
    pattern, sign, val = m.group(1), m.group(2), float(m.group(3))
    pct = -val if sign == 'm' else val
    delta = pct / 100.0
    eps = symmetry.strain_tensor(pattern, delta)
    gamma = symmetry.voigt_strain_engineering(eps)
    sigma_t = voigt_to_tensor(C_raw @ gamma)  # stress stays exact
    dE = 0.5 * V0 * (gamma @ C_raw @ gamma)
    if pattern == 'xx':
        dE *= 1.25  # simulated eggbox-corrupted energy curvature, xx only
    with open(f"{folder}/calc.out", "w") as f:
        f.write("siesta: Stress tensor (static) (eV/Ang**3):\n")
        for row in range(3):
            f.write("siesta: " + "   ".join(f"{sigma_t[row, c]:.10f}" for c in range(3)) + "\n")
        f.write(f"siesta: FreeEng = {dE:.10f}\n")
PYEOF
stb-elasticAnalysis --no-intro > log_contaminated.txt 2>&1
check_exit_code $? 0
check_contains "xx: stress =   166.00 GPa | energy =   207.50 GPa | diff = 25.00%.*WARNING" log_contaminated.txt
check_contains "xy: stress =    80.00 GPa | energy =    80.00 GPa | diff =  0.00%.*OK" log_contaminated.txt
check_contains "possible eggbox" log_contaminated.txt
popd > /dev/null

echo "Testing: no reference structure -> [WARNING] skip, rest of the report unaffected"
mkdir -p eggbox_noref_run
pushd eggbox_noref_run > /dev/null
for d in strain_xx_1.00 strain_xx_2.00 strain_xx_m1.00; do
    mkdir -p "$d"
    cp "$FIXTURE_DIR/$d/calc.out" "$d/"
done
stb-elasticAnalysis --no-intro > log_noref_eggbox.txt 2>&1
check_exit_code $? 0
check_contains "Eggbox cross-check skipped: no readable reference_structure.fdf" log_noref_eggbox.txt
check_contains "STABILITY AND PROPERTIES" log_noref_eggbox.txt
popd > /dev/null

echo "Testing: --method energy never runs the eggbox check"
mkdir -p eggbox_energy_run
pushd eggbox_energy_run > /dev/null
cp "$FIXTURE_DIR/../prep/si_cubic.fdf" .
stb-elasticInputs --file si_cubic.fdf --method energy --no-intro > log_prep.txt 2>&1
python3 - <<'PYEOF'
import os, re
import numpy as np
from stb.core import structure_io, symmetry

CONV = 160.21766
structure = structure_io.read_fdf("reference_structure.fdf")
pmg = structure_io.to_pymatgen(structure)
V0 = pmg.volume
C_true = np.zeros((6, 6))
C_true[0,0] = C_true[1,1] = C_true[2,2] = 166.0
C_true[0,1] = C_true[1,0] = C_true[0,2] = C_true[2,0] = C_true[1,2] = C_true[2,1] = 64.0
C_true[3,3] = C_true[4,4] = C_true[5,5] = 80.0
C_raw = C_true / CONV
for folder in os.listdir('.'):
    if not folder.startswith('strain_'):
        continue
    m = re.match(r"strain_([a-z+]+)_(m?)(\d+\.\d+)", folder)
    pattern, sign, val = m.group(1), m.group(2), float(m.group(3))
    pct = -val if sign == 'm' else val
    delta = pct / 100.0
    eps = sum((symmetry.strain_tensor(p, delta) for p in pattern.split('+')), np.zeros((3, 3)))
    gamma = symmetry.voigt_strain_engineering(eps)
    dE = 0.5 * V0 * (gamma @ C_raw @ gamma)
    with open(f"{folder}/calc.out", "w") as f:
        f.write(f"siesta: FreeEng = {dE:.10f}\n")
PYEOF
stb-elasticAnalysis --method energy --no-intro > log_energy_noeggbox.txt 2>&1
check_exit_code $? 0
if grep -qi "EGGBOX" log_energy_noeggbox.txt 2>/dev/null; then
    echo -e "   -> ${RED}Failed:${NC} --method energy unexpectedly ran the eggbox check"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no eggbox section for --method energy"
    PASS=$((PASS+1))
fi
popd > /dev/null

echo "Testing: per-direction fit quality -- residual stress and low R^2 flagged independently"
mkdir -p quality_run
pushd quality_run > /dev/null
cp "$FIXTURE_DIR/../prep/si_cubic.fdf" .
stb-elasticInputs --file si_cubic.fdf --dirs xx yy zz xy xz yz --no-intro > log_prep.txt 2>&1
python3 - <<'PYEOF'
import os, re
import numpy as np
from stb.core import structure_io, symmetry

CONV = 160.21766
structure = structure_io.read_fdf("reference_structure.fdf")
pmg = structure_io.to_pymatgen(structure)
V0 = pmg.volume

C_true = np.zeros((6, 6))
C_true[0,0] = C_true[1,1] = C_true[2,2] = 166.0
C_true[0,1] = C_true[1,0] = C_true[0,2] = C_true[2,0] = C_true[1,2] = C_true[2,1] = 64.0
C_true[3,3] = C_true[4,4] = C_true[5,5] = 80.0
C_raw = C_true / CONV

def voigt_to_tensor(v):
    s = np.zeros((3, 3))
    s[0,0], s[1,1], s[2,2], s[1,2], s[0,2], s[0,1] = v
    s[2,1], s[2,0], s[1,0] = s[1,2], s[0,2], s[0,1]
    return s

rng = np.random.default_rng(42)
for folder in os.listdir('.'):
    if not folder.startswith('strain_'):
        continue
    m = re.match(r"strain_([a-z+]+)_(m?)(\d+\.\d+)", folder)
    pattern, sign, val = m.group(1), m.group(2), float(m.group(3))
    pct = -val if sign == 'm' else val
    delta = pct / 100.0
    eps = symmetry.strain_tensor(pattern, delta)
    gamma = symmetry.voigt_strain_engineering(eps)
    sigma_v = C_raw @ gamma
    if pattern == 'xx':
        sigma_v[0] += 2.0 / CONV  # constant residual stress (unrelaxed reference), well above threshold
    if pattern == 'zz':
        sigma_v[2] += rng.normal(0, 1.0) / CONV  # scattered data, should tank R^2
    sigma_t = voigt_to_tensor(sigma_v)
    dE = 0.5 * V0 * (gamma @ C_raw @ gamma)
    with open(f"{folder}/calc.out", "w") as f:
        f.write("siesta: Stress tensor (static) (eV/Ang**3):\n")
        for row in range(3):
            f.write("siesta: " + "   ".join(f"{sigma_t[row, c]:.10f}" for c in range(3)) + "\n")
        f.write(f"siesta: FreeEng = {dE:.10f}\n")
PYEOF
stb-elasticAnalysis --no-intro > log_quality.txt 2>&1
check_exit_code $? 0
check_contains "xx: R\^2 = 1.000000 | residual stress =    2.00 GPa.*WARNING" log_quality.txt
check_contains "yy: R\^2 = 1.000000 | residual stress =    0.00 GPa.*OK" log_quality.txt
check_contains "zz: R\^2 = 0..* | residual stress =.*WARNING" log_quality.txt
check_contains "R\^2 below --fit-quality-tolerance" log_quality.txt
popd > /dev/null

echo "Testing: tensor symmetry (C vs C^T) -- flagged for --symmetry-method basic on independent (unpooled) directions"
mkdir -p symcheck_run
pushd symcheck_run > /dev/null
cp "$FIXTURE_DIR/../prep/triclinic.fdf" .
stb-elasticInputs --file triclinic.fdf --dirs xx yy zz xy xz yz --no-intro > log_prep.txt 2>&1
python3 - <<'PYEOF'
import os, re
import numpy as np
from stb.core import structure_io, symmetry

CONV = 160.21766
structure = structure_io.read_fdf("reference_structure.fdf")
pmg = structure_io.to_pymatgen(structure)
V0 = pmg.volume

C_true = np.array([
    [200, 50, 45, 5, 8, -3],
    [50, 150, 40, 2, 0, 4],
    [45, 40, 180, 0, 0, 0],
    [5, 2, 0, 70, 1.5, 0],
    [8, 0, 0, 1.5, 65, -2],
    [-3, 4, 0, 0, -2, 75],
], dtype=float)
C_raw = C_true / CONV

def voigt_to_tensor(v):
    s = np.zeros((3, 3))
    s[0,0], s[1,1], s[2,2], s[1,2], s[0,2], s[0,1] = v
    s[2,1], s[2,0], s[1,0] = s[1,2], s[0,2], s[0,1]
    return s

for folder in os.listdir('.'):
    if not folder.startswith('strain_'):
        continue
    m = re.match(r"strain_([a-z+]+)_(m?)(\d+\.\d+)", folder)
    pattern, sign, val = m.group(1), m.group(2), float(m.group(3))
    pct = -val if sign == 'm' else val
    delta = pct / 100.0
    eps = symmetry.strain_tensor(pattern, delta)
    gamma = symmetry.voigt_strain_engineering(eps)
    sigma_v = C_raw @ gamma
    if pattern == 'yy':
        # inflates the SLOPE of sigma_xx vs strain for the yy series only,
        # so C[0,1] (from yy's own column) disagrees with C[1,0] (xx's own
        # column, untouched) by exactly 15 GPa -- triclinic (point group
        # -1) has no symmetry-equivalent directions, so nothing pools these
        # two together and masks the disagreement (verified live: the same
        # perturbation on the cubic fixture gets silently averaged away by
        # --symmetry-method basic's pooling across {xx,yy,zz}).
        sigma_v[0] += (15.0 / CONV) * delta
    sigma_t = voigt_to_tensor(sigma_v)
    dE = 0.5 * V0 * (gamma @ C_raw @ gamma)
    with open(f"{folder}/calc.out", "w") as f:
        f.write("siesta: Stress tensor (static) (eV/Ang**3):\n")
        for row in range(3):
            f.write("siesta: " + "   ".join(f"{sigma_t[row, c]:.10f}" for c in range(3)) + "\n")
        f.write(f"siesta: FreeEng = {dE:.10f}\n")
PYEOF
stb-elasticAnalysis --no-intro > log_symcheck.txt 2>&1
check_exit_code $? 0
check_contains "Tensor symmetry (C vs C\^T): max |C_ij - C_ji| = 15.000.*i=1, j=2.*WARNING" log_symcheck.txt
check_contains "C is not close to symmetric" log_symcheck.txt

echo "Testing: --symmetry-method full never runs the tensor symmetry check (always exact by construction)"
stb-elasticAnalysis --symmetry-method full --no-intro > log_symcheck_full.txt 2>&1
check_exit_code $? 0
if grep -q "Tensor symmetry" log_symcheck_full.txt 2>/dev/null; then
    echo -e "   -> ${RED}Failed:${NC} --symmetry-method full unexpectedly ran the tensor symmetry check"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no tensor symmetry line for --symmetry-method full"
    PASS=$((PASS+1))
fi
popd > /dev/null


# --- 12. Physics: universal anisotropy index (A^U) and compliance matrix
#     (S = C^-1) -- both reuse numbers check_stability_and_report/the
#     stiffness-matrix printer already compute internally (B_v/G_v/B_r/G_r
#     for the Hill average; S was only ever used for the Reuss bounds), now
#     also surfaced directly. Verified against real silicon (the si_cubic
#     fixture's C11=166/C12=64/C44=80 GPa are literal textbook Si values):
#     A^U should land in silicon's known ~0.2-0.3 range, and S11 should
#     equal 1/E[100] (~130 GPa for Si), S44 should equal 1/C44 exactly
#     (cubic symmetry). ---
echo -e "\n--- Testing the universal anisotropy index and compliance matrix ---"

echo "Testing: cubic Si fixture -- A^U in the known ~0.2-0.3 range, S11 = 1/E[100], S44 = 1/C44"
mkdir -p compliance_cubic_run
pushd compliance_cubic_run > /dev/null
cp "$FIXTURE_DIR/../prep/si_cubic.fdf" .
stb-elasticInputs --file si_cubic.fdf --dirs all --symmetry-method full --no-intro > log_prep.txt 2>&1
python3 - <<'PYEOF'
import os, re
import numpy as np
from stb.core import structure_io, symmetry

CONV = 160.21766
structure = structure_io.read_fdf("reference_structure.fdf")
pmg = structure_io.to_pymatgen(structure)
V0 = pmg.volume

C_true = np.zeros((6, 6))
C_true[0,0] = C_true[1,1] = C_true[2,2] = 166.0
C_true[0,1] = C_true[1,0] = C_true[0,2] = C_true[2,0] = C_true[1,2] = C_true[2,1] = 64.0
C_true[3,3] = C_true[4,4] = C_true[5,5] = 80.0
C_raw = C_true / CONV

def voigt_to_tensor(v):
    s = np.zeros((3, 3))
    s[0,0], s[1,1], s[2,2], s[1,2], s[0,2], s[0,1] = v
    s[2,1], s[2,0], s[1,0] = s[1,2], s[0,2], s[0,1]
    return s

for folder in os.listdir('.'):
    if not folder.startswith('strain_'):
        continue
    m = re.match(r"strain_([a-z+]+)_(m?)(\d+\.\d+)", folder)
    pattern, sign, val = m.group(1), m.group(2), float(m.group(3))
    pct = -val if sign == 'm' else val
    delta = pct / 100.0
    eps = symmetry.strain_tensor(pattern, delta)
    gamma = symmetry.voigt_strain_engineering(eps)
    sigma_t = voigt_to_tensor(C_raw @ gamma)
    dE = 0.5 * V0 * (gamma @ C_raw @ gamma)
    with open(f"{folder}/calc.out", "w") as f:
        f.write("siesta: Stress tensor (static) (eV/Ang**3):\n")
        for row in range(3):
            f.write("siesta: " + "   ".join(f"{sigma_t[row, c]:.10f}" for c in range(3)) + "\n")
        f.write(f"siesta: FreeEng = {dE:.10f}\n")
PYEOF
stb-elasticAnalysis --no-intro > log_compliance_cubic.txt 2>&1
check_exit_code $? 0
check_contains "\[1b\] 3D COMPLIANCE MATRIX (1/GPa)" log_compliance_cubic.txt
check_contains "7.670e-03" log_compliance_cubic.txt
check_contains "1.250e-02" log_compliance_cubic.txt
check_contains "Universal Anisotropy Index (A\^U).*0.247" log_compliance_cubic.txt
popd > /dev/null

echo "Testing: 2D hexagonal fixture -- 2D compliance matrix, S11 = 1/Ex"
mkdir -p compliance_hex_run
pushd compliance_hex_run > /dev/null
cp "$FIXTURE_DIR/../prep/structure.fdf" .
stb-elasticInputs --file structure.fdf --dirs all --symmetry-method full --no-intro > log_prep.txt 2>&1
python3 - <<'PYEOF'
import os, re
import numpy as np
from stb.core import structure_io, symmetry

CONV = 16.021766
structure = structure_io.read_fdf("reference_structure.fdf")
pmg = structure_io.to_pymatgen(structure)
lat = pmg.lattice.matrix
Area0 = np.linalg.norm(np.cross(lat[0], lat[1]))

C11, C12 = 200.0, 50.0
C13, C33, C44 = 40.0, 300.0, 60.0
C66 = (C11 - C12) / 2
C_true = np.zeros((6, 6))
C_true[0,0] = C_true[1,1] = C11
C_true[0,1] = C_true[1,0] = C12
C_true[0,2] = C_true[2,0] = C_true[1,2] = C_true[2,1] = C13
C_true[2,2] = C33
C_true[3,3] = C_true[4,4] = C44
C_true[5,5] = C66
C_raw = C_true / CONV

def voigt_to_tensor(v):
    s = np.zeros((3, 3))
    s[0,0], s[1,1], s[2,2], s[1,2], s[0,2], s[0,1] = v
    s[2,1], s[2,0], s[1,0] = s[1,2], s[0,2], s[0,1]
    return s

for folder in os.listdir('.'):
    if not folder.startswith('strain_'):
        continue
    m = re.match(r"strain_([a-z+]+)_(m?)(\d+\.\d+)", folder)
    pattern, sign, val = m.group(1), m.group(2), float(m.group(3))
    pct = -val if sign == 'm' else val
    delta = pct / 100.0
    eps = symmetry.strain_tensor(pattern, delta)
    gamma = symmetry.voigt_strain_engineering(eps)
    sigma_t = voigt_to_tensor(C_raw @ gamma)
    dE = 0.5 * Area0 * (gamma @ C_raw @ gamma)
    with open(f"{folder}/calc.out", "w") as f:
        f.write("siesta: Stress tensor (static) (eV/Ang**3):\n")
        for row in range(3):
            f.write("siesta: " + "   ".join(f"{sigma_t[row, c]:.10f}" for c in range(3)) + "\n")
        f.write(f"siesta: FreeEng = {dE:.10f}\n")
PYEOF
stb-elasticAnalysis --2d --symmetry-method full --no-intro > log_compliance_hex.txt 2>&1
check_exit_code $? 0
check_contains "\[1b\] 2D COMPLIANCE MATRIX (m/N)" log_compliance_hex.txt
check_contains "5.333e-03" log_compliance_hex.txt
popd > /dev/null

echo "Testing: singular stiffness matrix -> compliance matrix warns instead of crashing"
mkdir -p compliance_singular_run
pushd compliance_singular_run > /dev/null
for d in strain_xx_1.00 strain_xx_2.00 strain_xx_m1.00; do
    mkdir -p "$d"
    cp "$FIXTURE_DIR/$d/calc.out" "$d/"
done
stb-elasticAnalysis --no-intro > log_singular.txt 2>&1
check_exit_code $? 0
check_contains "Could not invert the 3D stiffness matrix" log_singular.txt
popd > /dev/null


# --- 13. Physics: richer report ([0] RUN METADATA / [5] SUMMARY & FILES)
#     and the gnuplot data folder (--plot-dir, default elastic_plots/) --
#     one .dat/.gplot pair per direction/pattern + a combined overview,
#     written from data/fits already computed for the report (no new DFT
#     calculations, no re-fitting). If gnuplot itself is installed, actually
#     render one script to a PDF as an end-to-end sanity check. ---
echo -e "\n--- Testing report metadata/summary sections and the gnuplot data folder ---"

echo "Testing: --method stress -- elastic_plots/ has per-direction + combined dat/gplot pairs"
mkdir -p plots_stress_run
pushd plots_stress_run > /dev/null
cp "$FIXTURE_DIR/../prep/si_cubic.fdf" .
stb-elasticInputs --file si_cubic.fdf --dirs all --symmetry-method full --no-intro > log_prep.txt 2>&1
python3 - <<'PYEOF'
import os, re
import numpy as np
from stb.core import structure_io, symmetry

CONV = 160.21766
structure = structure_io.read_fdf("reference_structure.fdf")
pmg = structure_io.to_pymatgen(structure)
V0 = pmg.volume

C_true = np.zeros((6, 6))
C_true[0,0] = C_true[1,1] = C_true[2,2] = 166.0
C_true[0,1] = C_true[1,0] = C_true[0,2] = C_true[2,0] = C_true[1,2] = C_true[2,1] = 64.0
C_true[3,3] = C_true[4,4] = C_true[5,5] = 80.0
C_raw = C_true / CONV

def voigt_to_tensor(v):
    s = np.zeros((3, 3))
    s[0,0], s[1,1], s[2,2], s[1,2], s[0,2], s[0,1] = v
    s[2,1], s[2,0], s[1,0] = s[1,2], s[0,2], s[0,1]
    return s

for folder in os.listdir('.'):
    if not folder.startswith('strain_'):
        continue
    m = re.match(r"strain_([a-z+]+)_(m?)(\d+\.\d+)", folder)
    pattern, sign, val = m.group(1), m.group(2), float(m.group(3))
    pct = -val if sign == 'm' else val
    delta = pct / 100.0
    eps = symmetry.strain_tensor(pattern, delta)
    gamma = symmetry.voigt_strain_engineering(eps)
    sigma_t = voigt_to_tensor(C_raw @ gamma)
    dE = 0.5 * V0 * (gamma @ C_raw @ gamma)
    with open(f"{folder}/calc.out", "w") as f:
        f.write("siesta: Stress tensor (static) (eV/Ang**3):\n")
        for row in range(3):
            f.write("siesta: " + "   ".join(f"{sigma_t[row, c]:.10f}" for c in range(3)) + "\n")
        f.write(f"siesta: FreeEng = {dE:.10f}\n")
PYEOF
stb-elasticAnalysis --no-intro > log_report.txt 2>&1
check_exit_code $? 0
check_contains "\[0\] RUN METADATA" log_report.txt
check_contains "Method               : stress (--symmetry-method basic)" log_report.txt
check_contains "Tolerances           : eggbox=5.0%" log_report.txt
check_success elastic_plots/xx_fit.dat
check_success elastic_plots/xx_fit.gplot
check_success elastic_plots/xy_fit.dat
check_success elastic_plots/xy_fit.gplot
check_success elastic_plots/all_directions_fit.dat
check_success elastic_plots/all_directions_fit.gplot
check_contains "C_xx = 166.0000 GPa" elastic_plots/xx_fit.dat
check_contains "f(x) = 166" elastic_plots/xx_fit.gplot
check_contains "\[5\] SUMMARY & FILES" log_report.txt
check_contains "E=163.27 GPa.*A\^U=0.247" log_report.txt
check_contains "Plot data   : elastic_plots/ (3 file pair(s)" log_report.txt
if command -v gnuplot > /dev/null 2>&1; then
    ( cd elastic_plots && gnuplot xx_fit.gplot > /dev/null 2>gnuplot_err.txt )
    check_success elastic_plots/xx_fit.pdf
    if [ -s elastic_plots/gnuplot_err.txt ]; then
        echo -e "   -> ${RED}Failed:${NC} gnuplot reported errors rendering xx_fit.gplot"
        cat elastic_plots/gnuplot_err.txt
        FAIL=$((FAIL+1))
    else
        echo -e "   -> ${GREEN}Verified:${NC} gnuplot rendered xx_fit.gplot to a PDF with no errors"
        PASS=$((PASS+1))
    fi
else
    echo "   (gnuplot not installed -- skipping actual render check)"
fi
popd > /dev/null

echo "Testing: --method energy -- elastic_plots/ has pattern-named dat/gplot pairs (incl. combined 'xx+yy')"
mkdir -p plots_energy_run
pushd plots_energy_run > /dev/null
cp "$FIXTURE_DIR/../prep/si_cubic.fdf" .
stb-elasticInputs --file si_cubic.fdf --method energy --no-intro > log_prep.txt 2>&1
python3 - <<'PYEOF'
import os, re
import numpy as np
from stb.core import structure_io, symmetry

CONV = 160.21766
structure = structure_io.read_fdf("reference_structure.fdf")
pmg = structure_io.to_pymatgen(structure)
V0 = pmg.volume
C_true = np.zeros((6, 6))
C_true[0,0] = C_true[1,1] = C_true[2,2] = 166.0
C_true[0,1] = C_true[1,0] = C_true[0,2] = C_true[2,0] = C_true[1,2] = C_true[2,1] = 64.0
C_true[3,3] = C_true[4,4] = C_true[5,5] = 80.0
C_raw = C_true / CONV
for folder in os.listdir('.'):
    if not folder.startswith('strain_'):
        continue
    m = re.match(r"strain_([a-z+]+)_(m?)(\d+\.\d+)", folder)
    pattern, sign, val = m.group(1), m.group(2), float(m.group(3))
    pct = -val if sign == 'm' else val
    delta = pct / 100.0
    eps = sum((symmetry.strain_tensor(p, delta) for p in pattern.split('+')), np.zeros((3, 3)))
    gamma = symmetry.voigt_strain_engineering(eps)
    dE = 0.5 * V0 * (gamma @ C_raw @ gamma)
    with open(f"{folder}/calc.out", "w") as f:
        f.write(f"siesta: FreeEng = {dE:.10f}\n")
PYEOF
stb-elasticAnalysis --method energy --no-intro > log_report.txt 2>&1
check_exit_code $? 0
check_contains "\[0\] RUN METADATA" log_report.txt
check_contains "Method               : energy" log_report.txt
check_success "elastic_plots/xx+yy_fit.dat"
check_success "elastic_plots/xx+yy_fit.gplot"
check_contains "f(x) = " "elastic_plots/xx+yy_fit.gplot"
check_contains "x\*\*2" "elastic_plots/xx+yy_fit.gplot"
check_success elastic_plots/all_patterns_fit.dat
check_contains "\[5\] SUMMARY & FILES" log_report.txt
popd > /dev/null


# --- 14. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: no strain_* folders"
mkdir -p empty_run
pushd empty_run > /dev/null
stb-elasticAnalysis --no-intro > ../log_empty.txt 2>&1
check_exit_code $? 1
popd > /dev/null
check_contains "No valid data found in strain folders" log_empty.txt

echo "Testing: --version"
stb-elasticAnalysis --version > log_version.txt 2>&1
check_contains "stb-elasticAnalysis" log_version.txt

echo "Testing: --help documents --2d/--file"
stb-elasticAnalysis --help > log_help.txt 2>&1
check_contains "2D" log_help.txt
check_contains "file" log_help.txt
check_contains "symmetry-method" log_help.txt
check_contains "method" log_help.txt
check_contains "eggbox-tolerance" log_help.txt
check_contains "fit-quality-tolerance" log_help.txt
check_contains "plot-dir" log_help.txt


# --- 15. Interactive path (stb-suite, shortcut 4.2.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.2.2) ---"

echo "Testing: navigate 4.2.2 -> calc.out -> not 2D -> quit"
rm -f mechanical_properties.txt
printf '4.2.2\ncalc.out\nn\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_contains "STABILITY AND PROPERTIES" log_menu.txt
check_success mechanical_properties.txt
if grep -q "Traceback" log_menu.txt 2>/dev/null; then
    echo -e "   -> ${RED}Failed:${NC} unexpected traceback in interactive session"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} interactive session exited cleanly"
    PASS=$((PASS+1))
fi


popd > /dev/null

# --- 16. Summary ---
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
