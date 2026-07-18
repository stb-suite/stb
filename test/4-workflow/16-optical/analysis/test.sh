#!/bin/bash

# --- Setup ---
# Smoke test for stb-opticalAnalysis (Optical Properties Stage 2:
# Analysis, item 4.16.2). Chains a real stb-optical (Stage 1, tested
# separately) run on the cubic.fdf fixture, then fabricates a Lorentz-
# oscillator SystemLabel.EPSIMG (closed-form eps1(w)+eps2(w), distinct
# w0 per axis) in each dir_{x,y,z} folder plus a converged calc.out --
# same "fabricate a closed-form Lorentz oscillator" trick
# core/dielectric.py's own docstring says was used to validate
# kramers_kronig (~0.2% error).
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREP_DIR="$(cd "$FIXTURE_DIR/../prep" && pwd)"
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

# --- 1. Preparation ---
echo "--- Starting tester for STB-OPTICALANALYSIS stage 2: analysis (item 4.16.2) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$PREP_DIR/cubic.fdf" "$TEST_DIR/"
cp "$PREP_DIR/calc.fdf" "$TEST_DIR/"
echo "# placeholder pseudopotential" > "$TEST_DIR/Na.psf"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null

echo -e "\n--- Running stb-optical (Stage 1) once on cubic.fdf to build direction folders ---"
rm -rf optical_study
stb-optical -f cubic.fdf -c calc.fdf -p . -O optical_study --no-intro > /dev/null 2>&1
check_success optical_study/dir_x/structure.fdf
check_success optical_study/dir_y/structure.fdf
check_success optical_study/dir_z/structure.fdf

echo -e "\n--- Fabricating a closed-form Lorentz oscillator EPSIMG per direction ---"
python3 - <<'PYEOF'
import numpy as np

# Distinct (w0, wp, gamma) per axis -- lets the analysis gate below
# confirm the 3 directions never get mixed up.
params = {"x": (5.0, 3.0, 0.3), "y": (7.0, 2.5, 0.4), "z": (3.0, 4.0, 0.2)}
w = np.linspace(0.05, 15.0, 1500)

for axis, (w0, wp, gamma) in params.items():
    denom = (w0**2 - w**2)**2 + (gamma*w)**2
    eps2 = wp**2 * gamma * w / denom
    with open(f"optical_study/dir_{axis}/siesta.EPSIMG", "w") as f:
        f.write("# energy(eV)  eps2\n")
        for wi, e2i in zip(w, eps2):
            f.write(f"{wi:.6f}  {e2i:.8f}\n")
    with open(f"optical_study/dir_{axis}/calc.out", "w") as f:
        f.write("siesta: SCF Convergence by DM criterion\n")
        f.write("SCF cycle converged after 10 iterations\n")
print("fixtures written")
PYEOF


# --- 2. Missing-results guard ---
echo -e "\n--- Testing the missing-results guard ---"
rm -rf optical_study_empty
mkdir -p optical_study_empty
stb-opticalAnalysis --directory optical_study_empty --no-intro > log_no_results.txt 2>&1
check_exit_code $? 1
check_contains "run stb-optical" log_no_results.txt


# --- 3. Aggregated run: independently recompute eps1/n/k/alpha/R/L/sigma1 ---
echo -e "\n--- Testing aggregated analysis against the closed-form Lorentz oscillator ---"
stb-opticalAnalysis --directory optical_study --no-intro > log_analysis.txt 2>&1
check_exit_code $? 0
check_success optical_study/optical_results.csv
check_success optical_study/optical_results_x.dat
check_success optical_study/optical_results_y.dat
check_success optical_study/optical_results_z.dat
check_success optical_study/optical_results.gplot
check_contains "Folders found   : 3" log_analysis.txt

python3 -c "
import csv
import numpy as np
from scipy import constants

params = {'x': (5.0, 3.0, 0.3), 'y': (7.0, 2.5, 0.4), 'z': (3.0, 4.0, 0.2)}
HC = constants.h * constants.c / constants.e * 100.0
E_test = 6.0

with open('optical_study/optical_results.csv') as f:
    rows = list(csv.DictReader(f))

def analytic(w, w0, wp, gamma):
    denom = (w0**2 - w**2)**2 + (gamma*w)**2
    eps1 = 1.0 + wp**2 * (w0**2 - w**2) / denom
    eps2 = wp**2 * gamma * w / denom
    return eps1, eps2

def relerr(a, b):
    return abs(a - b) / abs(b) if b != 0 else abs(a - b)

for axis, (w0, wp, gamma) in params.items():
    eps1_e, eps2_e = analytic(E_test, w0, wp, gamma)
    eps_mag = np.sqrt(eps1_e**2 + eps2_e**2)
    n_e = np.sqrt((eps_mag + eps1_e) / 2)
    k_e = np.sqrt((eps_mag - eps1_e) / 2)
    alpha_e = 4 * np.pi * k_e * E_test / HC
    R_e = ((n_e - 1)**2 + k_e**2) / ((n_e + 1)**2 + k_e**2)
    L_e = eps2_e / (eps1_e**2 + eps2_e**2)
    omega_rad = E_test * constants.e / constants.hbar
    sigma1_e = constants.epsilon_0 * omega_rad * eps2_e

    axis_rows = [r for r in rows if r['direction'] == axis]
    assert axis_rows, f'no rows found for direction {axis}'
    closest = min(axis_rows, key=lambda r: abs(float(r['E_eV']) - E_test))

    assert relerr(float(closest['eps1']), eps1_e) < 0.03, (axis, 'eps1', closest['eps1'], eps1_e)
    assert relerr(float(closest['eps2']), eps2_e) < 0.03, (axis, 'eps2', closest['eps2'], eps2_e)
    assert relerr(float(closest['n']), n_e) < 0.03, (axis, 'n', closest['n'], n_e)
    assert relerr(float(closest['k']), k_e) < 0.03, (axis, 'k', closest['k'], k_e)
    assert relerr(float(closest['alpha_cm-1']), alpha_e) < 0.03, (axis, 'alpha')
    assert relerr(float(closest['R']), R_e) < 0.03, (axis, 'R', closest['R'], R_e)
    assert relerr(float(closest['L']), L_e) < 0.05, (axis, 'L', closest['L'], L_e)
    assert relerr(float(closest['sigma1_Spm']), sigma1_e) < 0.03, (axis, 'sigma1')

print('OK')
" > log_analysis_check.txt 2>&1
check_contains "OK" log_analysis_check.txt


# --- 3b. Isotropic average (all 3 directions present) ---
echo -e "\n--- Testing the isotropic ('avg') average when x/y/z are all present ---"
check_contains "Direction avg" log_analysis.txt

python3 -c "
import csv
with open('optical_study/optical_results.csv') as f:
    rows = list(csv.DictReader(f))
by_axis = {}
for r in rows:
    by_axis.setdefault(r['direction'], []).append(r)
assert 'avg' in by_axis, 'no avg direction in CSV'
# eps1(E->0) for avg must equal the arithmetic mean of x/y/z's own eps1(E->0)
def eps1_at_min_e(axis):
    axis_rows = by_axis[axis]
    return min(axis_rows, key=lambda r: float(r['E_eV']))
e_x = float(eps1_at_min_e('x')['eps1'])
e_y = float(eps1_at_min_e('y')['eps1'])
e_z = float(eps1_at_min_e('z')['eps1'])
e_avg = float(eps1_at_min_e('avg')['eps1'])
expected = (e_x + e_y + e_z) / 3.0
assert abs(e_avg - expected) < 1e-6, (e_avg, expected)
print('OK')
" > log_avg_check.txt 2>&1
check_contains "OK" log_avg_check.txt


# --- 3c. --plot-quantity all ---
echo -e "\n--- Testing --plot-quantity all (writes one .gplot per quantity) ---"
stb-opticalAnalysis --directory optical_study --plot-quantity all -o all_quantities --no-intro \
    > log_all_quantities.txt 2>&1
check_exit_code $? 0
check_success optical_study/all_quantities_eps.gplot
check_success optical_study/all_quantities_n_k.gplot
check_success optical_study/all_quantities_alpha.gplot
check_success optical_study/all_quantities_R.gplot
check_success optical_study/all_quantities_L.gplot
check_success optical_study/all_quantities_sigma1.gplot


# --- 3d. --experimental (peak overlay + matching) ---
echo -e "\n--- Testing --experimental (peak overlay + matching against the isotropic average) ---"
python3 -c "
import numpy as np
w = np.linspace(1, 12, 200)
# 2 Gaussian peaks, deliberately close to the fabricated Lorentz oscillators'
# own alpha peaks (around 3 and 5 eV for this parameter set) so the nearest-
# neighbor peak match below has something meaningful to find.
alpha_fake = 5e5*np.exp(-((w-5.2)/0.5)**2) + 3e5*np.exp(-((w-3.1)/0.4)**2)
with open('uvvis_fake.dat', 'w') as f:
    f.write('# energy(eV) absorbance(arb)\n')
    for wi, ai in zip(w, alpha_fake):
        f.write(f'{wi:.4f} {ai:.4f}\n')
"
stb-opticalAnalysis --directory optical_study --plot-quantity alpha --experimental uvvis_fake.dat \
    -o exp_results --no-intro > log_experimental.txt 2>&1
check_exit_code $? 0
check_contains "\[3b\] EXPERIMENTAL COMPARISON" log_experimental.txt
check_contains "Simulated reference : direction 'avg'" log_experimental.txt
check_success optical_study/exp_results_experimental.dat
check_contains "axes x1y2" optical_study/exp_results.gplot

echo "Testing: --experimental + a paired --plot-quantity (eps) is rejected"
stb-opticalAnalysis --directory optical_study --plot-quantity eps --experimental uvvis_fake.dat \
    --no-intro > log_experimental_paired.txt 2>&1
check_exit_code $? 1
check_contains "single-valued" log_experimental_paired.txt

echo "Testing: missing --experimental file is a clear error, not a crash"
stb-opticalAnalysis --directory optical_study --plot-quantity alpha --experimental nonexistent.dat \
    --no-intro > log_experimental_missing.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_experimental_missing.txt


# --- 4. Only one direction present (partial results, still runs OK) ---
echo -e "\n--- Testing partial results (only dir_x present) ---"
rm -rf optical_study_partial
cp -r optical_study optical_study_partial
rm -rf optical_study_partial/dir_y optical_study_partial/dir_z
stb-opticalAnalysis --directory optical_study_partial -o partial_results --no-intro \
    > log_partial.txt 2>&1
check_exit_code $? 0
check_contains "Folders found   : 1" log_partial.txt
check_success optical_study_partial/partial_results_x.dat


# --- 5. Direction detected from Optical.Vector block, not folder name ---
echo -e "\n--- Testing direction detection is content-based, not folder-name-based ---"
rm -rf optical_study_renamed
cp -r optical_study optical_study_renamed
mv optical_study_renamed/dir_x optical_study_renamed/some_random_folder_name
stb-opticalAnalysis --directory optical_study_renamed -o renamed_results --no-intro \
    > log_renamed.txt 2>&1
check_exit_code $? 0
check_contains "Direction x (optical_study_renamed/some_random_folder_name)" log_renamed.txt
check_success optical_study_renamed/renamed_results_x.dat


# --- 6. Unconverged SCF -> warning, does not block ---
echo -e "\n--- Testing unconverged-SCF warning ---"
rm -rf optical_study_unconv
cp -r optical_study optical_study_unconv
echo "siesta: did not converge" > optical_study_unconv/dir_x/calc.out
stb-opticalAnalysis --directory optical_study_unconv -o unconv_results --no-intro \
    > log_unconv.txt 2>&1
check_exit_code $? 0
check_contains "Could not confirm SCF convergence" log_unconv.txt
check_success optical_study_unconv/unconv_results_x.dat


# --- 7. Error cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: --version"
stb-opticalAnalysis --version > log_version.txt 2>&1
check_contains "stb-opticalAnalysis" log_version.txt

echo "Testing: --help documents --plot-quantity/--energy-min/--energy-max/--output/--experimental/--peak-prominence"
stb-opticalAnalysis --help > log_help.txt 2>&1
check_contains "plot-quantity" log_help.txt
check_contains "energy-min" log_help.txt
check_contains "energy-max" log_help.txt
check_contains "output" log_help.txt
check_contains "experimental" log_help.txt
check_contains "peak-prominence" log_help.txt


# --- 8. Interactive path (stb-suite, shortcut 4.16.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.16.2) ---"

echo "Testing: navigate 4.16.2 -> optical_study -> alpha -> quit"
{
  echo "4.16.2"
  echo "optical_study"    # directory
  echo ""                   # output filename (default calc.out)
  echo "3"                   # alpha
  echo ""                     # advanced settings (default N)
  echo ""                      # press enter to continue
  echo "0"                      # quit stage submenu
} | stb-suite > log_menu.txt 2>&1
check_success optical_study/optical_results.csv
check_success optical_study/optical_stage2.txt


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
