#!/bin/bash

# --- Setup ---
# Smoke test for stb-irAnalysis (IR Spectrum, Stage 3: Analysis, item
# 4.12.3). Builds minimal ir_study/ trees BY HAND (writing ir_stage2.txt's
# # MODE_TABLE + metadata lines directly, same "only Stage 3's reading side
# is under test" discipline as test/.../11-raman/analysis/test.sh's own
# "4c" reconstruction test) -- avoids a real phonon/MACE calculation for
# every fixture, since Stage 2's own CLI wiring (both paths, including a
# real degenerate-mode collapse) is already covered end to end in
# test/.../12-ir/modes/test.sh.
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


# --- 1. Preparation ---
echo "--- Starting tester for STB-IRAnalysis stage 3: analysis (item 4.12.3) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$PREP_DIR/structure.fdf" "$TEST_DIR/"
cp "$PREP_DIR/calc.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Missing Stage-2 output guard ---
echo -e "\n--- Testing the missing-Stage-2-output guard ---"
rm -rf ir_study
stb-irAnalysis --directory ir_study --no-intro > log_no_stage2.txt 2>&1
check_exit_code $? 1
check_contains "run stb-irModes" log_no_stage2.txt


# --- 2b. core.born_charges.read_born_charges: both real .BC layouts ---
echo -e "\n--- Testing read_born_charges against both the unlabeled (real, confirmed) and "
echo "    x/y/z-labeled (SIESTA tutorial material) .BC layouts ---"
python3 -c "
import numpy as np
from stb.core.born_charges import read_born_charges

# Real confirmed layout: bare header, N*3 unlabeled rows.
with open('/tmp/unlabeled.BC', 'w') as f:
    f.write('BC matrix\n')
    f.write('   1.0000000   0.0000000   0.0000000\n')
    f.write('   0.0000000   1.0000000   0.0000000\n')
    f.write('   0.0000000   0.0000000   1.0000000\n')
labels, Z = read_born_charges('/tmp/unlabeled.BC')
assert Z.shape == (1, 3, 3), Z.shape
assert np.allclose(Z[0], np.eye(3)), Z

# x/y/z-labeled layout (never seen live, but keep supporting it).
with open('/tmp/labeled.BC', 'w') as f:
    f.write('BC matrix\n')
    f.write('            x              y              z\n')
    f.write('x    2.0000000   0.0000000   0.0000000\n')
    f.write('y    0.0000000   2.0000000   0.0000000\n')
    f.write('z    0.0000000   0.0000000   2.0000000\n')
labels, Z = read_born_charges('/tmp/labeled.BC')
assert Z.shape == (1, 3, 3), Z.shape
assert np.allclose(Z[0], np.eye(3) * 2), Z
print('OK')
" > log_bc_parser_check.txt 2>&1
check_contains "OK" log_bc_parser_check.txt


# --- 3. Bulk path: known Born charges + eigendisplacements -> exact dmu/dQ ---
echo -e "\n--- Testing bulk-path analysis against known ground truth (fabricated Z* + eigendisplacement) ---"
rm -rf ir_study
python3 - <<'PYEOF'
import json, os
import numpy as np

os.makedirs("ir_study/born_charge_disp/equilibrium", exist_ok=True)

# 3 atoms, simple KNOWN Z* tensors (diagonal, scaled by atom index) and a
# KNOWN eigendisplacement per mode -- independently recomputed below via
# plain numpy, same "recompute the expected answer independently from the
# same fabricated closed-form data" discipline as raman_analysis's own
# EPSIMG-based tests.
Z_star = np.array([np.eye(3) * (1.0 + i) for i in range(3)])
np.save("ground_truth_Z.npy", Z_star)

eigendisp_mode1 = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
eigendisp_mode2 = np.array([[0.5, 0.5, 0.0], [-0.5, 0.5, 0.0], [0.0, 0.0, 1.0]])
np.save("ground_truth_e1.npy", eigendisp_mode1)
np.save("ground_truth_e2.npy", eigendisp_mode2)

with open("ir_study/born_charge_disp/mode_01_eigendisplacement.json", "w") as f:
    json.dump({"eigendisplacement": eigendisp_mode1.tolist()}, f)
with open("ir_study/born_charge_disp/mode_02_eigendisplacement.json", "w") as f:
    json.dump({"eigendisplacement": eigendisp_mode2.tolist()}, f)

with open("ir_study/born_charge_disp/equilibrium/Sn3O4.BC", "w") as f:
    # Real confirmed SIESTA .BC layout (MaX-1.2/PSML build, 2026-07): a
    # bare "BC matrix" header, then N*3 UNLABELED rows (no per-row x/y/z
    # prefix, no column header) -- NOT the x/y/z-labeled layout originally
    # guessed from SIESTA's tutorial material, which this exact fixture
    # used to write (and which passed even though it never exercised the
    # real unlabeled format -- caught live against a real .BC file, see
    # core/born_charges.py's own docstring).
    f.write("BC matrix\n")
    for Z in Z_star:
        for row in Z:
            f.write(f"   {row[0]: .7f}   {row[1]: .7f}   {row[2]: .7f}\n")

with open("ir_study/born_charge_disp/equilibrium/calc.fdf", "w") as f:
    f.write("SystemLabel Sn3O4\n")
with open("ir_study/born_charge_disp/equilibrium/calc.out", "w") as f:
    f.write("siesta: SCF Convergence by DM criterion\n")
    f.write("SCF cycle converged after 12 iterations\n")

with open("ir_study/ir_stage2.txt", "w") as f:
    f.write("Displacement      : 0.02 Ang\n")
    f.write("Vacuum axes (abc) : False False False\n")
    f.write("Lattice vectors (Ang, for Stage 3's dipole-axis alignment check):\n")
    f.write("  a  5.000000  0.000000  0.000000\n")
    f.write("  b  0.000000  5.000000  0.000000\n")
    f.write("  c  0.000000  0.000000  5.000000\n")
    f.write("\n# MODE_TABLE -- parsed by stb-irAnalysis, do not reorder the first 7 columns\n")
    f.write(f"# {'label':<22}{'mode_index':<12}{'band_index':<12}{'frequency_thz':<16}"
            f"{'path':<8}{'sign':<8}{'dir':<24}{'derived_from'}\n")
    eq_dir = "ir_study/born_charge_disp/equilibrium"
    f.write(f"{'mode_01_born_charge':<24}{1:<12}{3:<12}{10.0:<16.6f}{'BULK':<8}{'-':<8}{eq_dir}  -\n")
    f.write(f"{'mode_02_born_charge':<24}{2:<12}{4:<12}{12.0:<16.6f}{'BULK':<8}{'-':<8}{eq_dir}  -\n")
print("fixture ready")
PYEOF
stb-irAnalysis --directory ir_study --no-intro > log_bulk_analysis.txt 2>&1
check_exit_code $? 0
check_contains "Modes analyzed      : 2/2" log_bulk_analysis.txt

echo "Testing: reported dmu/dQ for mode 1 and 2 exactly match sum_atom Z*_atom @ e_atom(mode)"
python3 -c "
import re
import numpy as np
Z = np.load('ground_truth_Z.npy')
e1 = np.load('ground_truth_e1.npy')
e2 = np.load('ground_truth_e2.npy')
expected1 = np.einsum('atb,at->b', Z, e1)
expected2 = np.einsum('atb,at->b', Z, e2)
with open('ir_study/ir_stage3.txt') as f:
    text = f.read()
matches = re.findall(r'dmu/dQ=\(([\-\d.]+), ([\-\d.]+), ([\-\d.]+)\)', text)
assert len(matches) == 2, f'expected 2 dmu/dQ lines, got {len(matches)}'
got1 = np.array([float(v) for v in matches[0]])
got2 = np.array([float(v) for v in matches[1]])
assert np.allclose(got1, expected1, atol=1e-5), f'mode 1: got {got1} vs expected {expected1}'
assert np.allclose(got2, expected2, atol=1e-5), f'mode 2: got {got2} vs expected {expected2}'
print('OK')
" > log_bulk_check.txt 2>&1
check_contains "OK" log_bulk_check.txt


# --- 4. Non-bulk (0D) path: known linear dipole data -> exact central difference ---
echo -e "\n--- Testing non-bulk-path analysis against known ground truth (fabricated dipole data) ---"
rm -rf ir_study
python3 - <<'PYEOF'
import os
import numpy as np

delta = 0.02
np.save("ground_truth_slope.npy", np.array([1.0, 0.5, -0.3]))
slope = np.array([1.0, 0.5, -0.3])
base = np.array([0.01, 0.02, -0.01])

for sign, factor in (("plus", 1.0), ("minus", -1.0)):
    folder = f"ir_study/dipole_disp/mode_01_{sign}"
    os.makedirs(folder, exist_ok=True)
    mu = base + factor * delta * slope
    with open(os.path.join(folder, "calc.fdf"), "w") as f:
        f.write("SystemLabel H2O\n")
    with open(os.path.join(folder, "calc.out"), "w") as f:
        f.write("siesta: SCF Convergence by DM criterion\n")
        f.write("SCF cycle converged after 10 iterations\n")
        f.write(f"siesta: Electric dipole (a.u.) =    {mu[0]:.6f}    {mu[1]:.6f}    {mu[2]:.6f}\n")

with open("ir_study/ir_stage2.txt", "w") as f:
    f.write("Displacement      : 0.02 Ang\n")
    f.write("Vacuum axes (abc) : True True True\n")
    f.write("Lattice vectors (Ang, for Stage 3's dipole-axis alignment check):\n")
    f.write("  a  16.000000  0.000000  0.000000\n")
    f.write("  b  0.000000  16.000000  0.000000\n")
    f.write("  c  0.000000  0.000000  16.000000\n")
    f.write("\n# MODE_TABLE -- parsed by stb-irAnalysis, do not reorder the first 7 columns\n")
    f.write(f"# {'label':<22}{'mode_index':<12}{'band_index':<12}{'frequency_thz':<16}"
            f"{'path':<8}{'sign':<8}{'dir':<24}{'derived_from'}\n")
    f.write(f"{'mode_01_plus':<24}{1:<12}{6:<12}{44.5:<16.6f}{'NONBULK':<8}{'plus':<8}"
            f"ir_study/dipole_disp/mode_01_plus  -\n")
    f.write(f"{'mode_01_minus':<24}{1:<12}{6:<12}{44.5:<16.6f}{'NONBULK':<8}{'minus':<8}"
            f"ir_study/dipole_disp/mode_01_minus  -\n")
print("fixture ready")
PYEOF
stb-irAnalysis --directory ir_study --no-intro > log_nonbulk_analysis.txt 2>&1
check_exit_code $? 0
check_contains "Modes analyzed      : 1/1" log_nonbulk_analysis.txt

echo "Testing: reported dmu/dQ exactly matches the known slope (0D -- no axis masking applies)"
python3 -c "
import re
import numpy as np
slope = np.load('ground_truth_slope.npy')
with open('ir_study/ir_stage3.txt') as f:
    text = f.read()
m = re.search(r'dmu/dQ=\(([\-\d.]+), ([\-\d.]+), ([\-\d.]+)\)', text)
assert m, 'could not find dmu/dQ line'
got = np.array([float(v) for v in m.groups()])
assert np.allclose(got, slope, atol=1e-6), f'got {got} vs expected {slope}'
print('OK')
" > log_nonbulk_check.txt 2>&1
check_contains "OK" log_nonbulk_check.txt


# --- 5. Non-bulk (2D-like) path: dipole-axis masking ---
echo -e "\n--- Testing dipole-axis masking (periodic in-plane component discarded) ---"
rm -rf ir_study
python3 - <<'PYEOF'
import os
import numpy as np

delta = 0.02
# in-plane (x, y -- periodic, "a"/"b" vacuum=False) components are
# DELIBERATELY large; only the out-of-plane (z -- vacuum=True, "c") slope
# is the real physical signal Stage 3 must keep.
slope = np.array([5.0, -4.0, 0.7])
np.save("ground_truth_slope_2d.npy", slope)
base = np.array([0.1, -0.2, 0.01])

for sign, factor in (("plus", 1.0), ("minus", -1.0)):
    folder = f"ir_study/dipole_disp/mode_01_{sign}"
    os.makedirs(folder, exist_ok=True)
    mu = base + factor * delta * slope
    with open(os.path.join(folder, "calc.fdf"), "w") as f:
        f.write("SystemLabel graphene\n")
    with open(os.path.join(folder, "calc.out"), "w") as f:
        f.write("siesta: SCF Convergence by DM criterion\n")
        f.write("SCF cycle converged after 10 iterations\n")
        f.write(f"siesta: Electric dipole (a.u.) =    {mu[0]:.6f}    {mu[1]:.6f}    {mu[2]:.6f}\n")

with open("ir_study/ir_stage2.txt", "w") as f:
    f.write("Displacement      : 0.02 Ang\n")
    # a, b periodic (in-plane); c vacuum-padded (out-of-plane) -- graphene-like.
    f.write("Vacuum axes (abc) : False False True\n")
    f.write("Lattice vectors (Ang, for Stage 3's dipole-axis alignment check):\n")
    f.write("  a  2.460000  0.000000  0.000000\n")
    f.write("  b  -1.230000  2.130842  0.000000\n")
    f.write("  c  0.000000  0.000000  20.000000\n")
    f.write("\n# MODE_TABLE -- parsed by stb-irAnalysis, do not reorder the first 7 columns\n")
    f.write(f"# {'label':<22}{'mode_index':<12}{'band_index':<12}{'frequency_thz':<16}"
            f"{'path':<8}{'sign':<8}{'dir':<24}{'derived_from'}\n")
    f.write(f"{'mode_01_plus':<24}{1:<12}{3:<12}{25.0:<16.6f}{'NONBULK':<8}{'plus':<8}"
            f"ir_study/dipole_disp/mode_01_plus  -\n")
    f.write(f"{'mode_01_minus':<24}{1:<12}{3:<12}{25.0:<16.6f}{'NONBULK':<8}{'minus':<8}"
            f"ir_study/dipole_disp/mode_01_minus  -\n")
print("fixture ready")
PYEOF
stb-irAnalysis --directory ir_study --no-intro > log_mask_analysis.txt 2>&1
check_exit_code $? 0
check_contains "Non-negligible dipole-derivative component along a periodic direction" log_mask_analysis.txt

echo "Testing: only the out-of-plane (z) component survives masking, in-plane (x,y) zeroed"
python3 -c "
import re
import numpy as np
slope = np.load('ground_truth_slope_2d.npy')
with open('ir_study/ir_stage3.txt') as f:
    text = f.read()
m = re.search(r'dmu/dQ=\(([\-\d.]+), ([\-\d.]+), ([\-\d.]+)\)', text)
assert m, 'could not find dmu/dQ line'
got = np.array([float(v) for v in m.groups()])
expected = np.array([0.0, 0.0, slope[2]])
assert np.allclose(got, expected, atol=1e-6), f'got {got} vs expected {expected} (masked)'
print('OK')
" > log_mask_check.txt 2>&1
check_contains "OK" log_mask_check.txt


# --- 6. --temperature: Bose-Einstein weighting changes the spectrum ---
echo -e "\n--- Testing --temperature (Bose-Einstein thermal weighting) ---"
stb-irAnalysis --directory ir_study --temperature 300 -o ir_spectrum_300k --no-intro \
    > log_temperature.txt 2>&1
check_exit_code $? 0
check_contains "Thermal weighting : Bose-Einstein Stokes factor (n+1) at 300.0 K" log_temperature.txt
check_success ir_study/ir_spectrum_300k.dat


# --- 7. --experimental: overlay + nearest-neighbor peak match ---
echo -e "\n--- Testing --experimental (peak overlay + matching) ---"
python3 -c "
import numpy as np
freq = np.linspace(700, 900, 200)
inten = 100 * np.exp(-((freq - 800.0) ** 2) / (2 * 5.0 ** 2))
with open('fake_experimental.dat', 'w') as f:
    f.write('# fake experimental spectrum for testing\n')
    for fr, it in zip(freq, inten):
        f.write(f'{fr:.4f} {it:.6f}\n')
"
stb-irAnalysis --directory ir_study --experimental fake_experimental.dat --no-intro \
    > log_experimental.txt 2>&1
check_exit_code $? 0
check_contains "\[3b\] EXPERIMENTAL COMPARISON" log_experimental.txt
check_success ir_study/ir_spectrum_experimental.dat

echo "Testing: missing --experimental file is a clear error, not a crash"
stb-irAnalysis --directory ir_study --experimental does_not_exist.dat --no-intro \
    > log_experimental_missing.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_experimental_missing.txt


# --- 8. Error cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: --version"
stb-irAnalysis --version > log_version.txt 2>&1
check_contains "stb-irAnalysis" log_version.txt

echo "Testing: --help documents --linewidth/--file/--temperature/--experimental"
stb-irAnalysis --help > log_help.txt 2>&1
check_contains "linewidth" log_help.txt
check_contains "file" log_help.txt
check_contains "temperature" log_help.txt
check_contains "experimental" log_help.txt
check_contains "peak-prominence" log_help.txt


# --- 9. Interactive path (stb-suite, shortcut 4.12.3) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.12.3) ---"

echo "Testing: navigate 4.12.3 -> defaults -> quit"
{
  echo "4.12.3"
  echo ""               # run_dir (default ir_study)
  echo ""               # output_filename (default calc.out)
  echo ""               # linewidth (default 10.0)
  echo ""               # temperature (optional, skip)
  echo ""               # experimental_file (optional, skip)
  echo ""               # press enter to continue
  echo "0"              # quit stage submenu
} | stb-suite > log_menu.txt 2>&1
check_contains "Modes analyzed      : 1/1" log_menu.txt
check_success ir_study/ir_spectrum.dat


popd > /dev/null

# --- 10. Summary ---
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
