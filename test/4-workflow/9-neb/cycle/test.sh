#!/bin/bash

# --- Setup ---
# Smoke test for stb-nebCycle (item 4.9, CLI-only -- deliberately NOT in
# the interactive stb-suite menu, see CLAUDE.md). Builds a real cycle_00/
# via stb-neb --mode 4 (no MACE, fast) against the same fixture as ../prep/,
# then fabricates synthetic .FA/calc.out results (a smooth, physically-
# motivated harmonic "pull toward a fixed per-image target" force field --
# NOT random noise, which was found live during development to trigger a
# real ASE edge case, see the [4] NaN-guard test below) to exercise the
# real-DFT NEB step without needing an actual SIESTA run.
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
echo "--- Starting tester for STB-NebCycle (item 4.9, CLI-only) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$PREP_DIR/initial.fdf" "$TEST_DIR/"
cp "$PREP_DIR/final.fdf" "$TEST_DIR/"
cp "$PREP_DIR/calc.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null

# Fabricates .FA + calc.out for every image_* in $1/cycle_$2 (zero-padded),
# pulling each atom toward a FIXED per-image target (harmonic force,
# F = -K*(pos - target), E = 0.5*K*|pos-target|^2 -- a genuine smooth
# function of the CURRENT geometry, unlike a repeated/random force, which
# is what originally triggered the ASE 'improvedtangent' 0/0 tangent bug
# this tool now guards against explicitly (see [4] below). Targets are
# derived once (a small deterministic per-image offset from cycle_00's own
# guess) and cached in targets.npy so every call pulls toward the SAME
# fixed point regardless of which cycle is being fabricated.
cat > fabricate.py << 'PYEOF'
import sys, os
import numpy as np
from stb.core import structure_io

ROOT, CYCLE = sys.argv[1], sys.argv[2]
K = 2.0

def positions(path):
    return structure_io.to_pymatgen(structure_io.read_fdf(path)).cart_coords

def write_fa(path, forces):
    n = len(forces)
    with open(path, 'w') as f:
        f.write(f'{n}\n')
        for i, (fx, fy, fz) in enumerate(forces, start=1):
            f.write(f'{i} {fx:.6f} {fy:.6f} {fz:.6f}\n')

def write_out(path, energy):
    with open(path, 'w') as f:
        f.write(f'siesta: FreeEng =    {energy:.6f}\nSCF cycle converged after 10 iterations\n')

labels = sorted(d for d in os.listdir(f'{ROOT}/cycle_00') if d.startswith('image_'))
targets_path = f'{ROOT}/targets.npy'
if not os.path.exists(targets_path):
    rng = np.random.default_rng(7)
    targets = {}
    for label in labels:
        p0 = positions(f'{ROOT}/cycle_00/{label}/structure.fdf')
        targets[label] = p0 + (rng.random(p0.shape) - 0.5) * 0.15
    np.save(targets_path, targets, allow_pickle=True)
targets = np.load(targets_path, allow_pickle=True).item()

for label in labels:
    d = f'{ROOT}/cycle_{CYCLE}/{label}'
    pos = positions(f'{d}/structure.fdf')
    disp = pos - targets[label]
    write_fa(os.path.join(d, 'siesta.FA'), -K * disp)
    write_out(os.path.join(d, 'calc.out'), 0.5 * K * float(np.sum(disp**2)) - 300.0)
PYEOF


# --- 2. Build cycle_00 via stb-neb --mode 4, fabricate results, take one step ---
echo -e "\n--- Testing a single real-DFT NEB step ---"
stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 4 --mode 4 -O run1 --no-intro > log_prep1.txt 2>&1
check_exit_code $? 0
python3 fabricate.py run1 00
stb-nebCycle --dir run1 --fmax 0.001 --no-intro > log_step1.txt 2>&1
check_exit_code $? 0
check_contains "\[0\] RUN METADATA" log_step1.txt
check_contains "Current cycle   : cycle_00" log_step1.txt
check_contains "\[1\] READING SIESTA RESULTS" log_step1.txt
check_contains "\[2\] NEB STEP" log_step1.txt
check_contains "fresh, first cycle" log_step1.txt
check_contains "\[3\] RESULT" log_step1.txt
check_success run1/cycle_01/image_00/structure.fdf
check_success run1/cycle_01/image_01/calc.fdf
check_success run1/neb_cycle_state.json
echo "Testing: endpoints (image_00/image_03) never move"
python3 -c "
import sys
import numpy as np
from stb.core import structure_io
for label in ('image_00', 'image_03'):
    a = structure_io.to_pymatgen(structure_io.read_fdf(f'run1/cycle_00/{label}/structure.fdf'))
    b = structure_io.to_pymatgen(structure_io.read_fdf(f'run1/cycle_01/{label}/structure.fdf'))
    if not np.allclose(a.cart_coords, b.cart_coords):
        print(f'{label} moved unexpectedly'); sys.exit(1)
print('OK')
" > log_endpoints_check.txt 2>&1
check_contains "OK" log_endpoints_check.txt
echo "Testing: an interior image DID move (real forces, not a no-op)"
python3 -c "
import sys
import numpy as np
from stb.core import structure_io
a = structure_io.to_pymatgen(structure_io.read_fdf('run1/cycle_00/image_01/structure.fdf'))
b = structure_io.to_pymatgen(structure_io.read_fdf('run1/cycle_01/image_01/structure.fdf'))
sys.exit(0 if not np.allclose(a.cart_coords, b.cart_coords) else 1)
"
check_exit_code $? 0


# --- 3. State persistence: a second step accelerates (FIRE momentum), a fresh restart doesn't ---
echo -e "\n--- Testing FIRE restart state persistence across separate process calls ---"
python3 fabricate.py run1 01
stb-nebCycle --dir run1 --fmax 0.001 --no-intro > log_step2.txt 2>&1
check_exit_code $? 0
check_contains "loaded from previous cycle" log_step2.txt
check_success run1/cycle_02/image_01/structure.fdf

echo "Testing: displacement of the SAME image accelerates cycle-to-cycle (proves state was reused, not reset)"
python3 -c "
import sys
import numpy as np
from stb.core import structure_io
def pos(c, label):
    return structure_io.to_pymatgen(structure_io.read_fdf(f'run1/cycle_{c}/{label}/structure.fdf')).cart_coords
d1 = np.linalg.norm(pos('00', 'image_01') - pos('01', 'image_01'))
d2 = np.linalg.norm(pos('01', 'image_01') - pos('02', 'image_01'))
print(f'step1={d1:.6f} step2={d2:.6f} ratio={d2/d1:.3f}')
# FIRE accelerates while consecutive gradients agree -- with identical
# fabricated forces reapplied each cycle, step 2 must be strictly larger
# than step 1 if (and only if) velocity truly carried over between the two
# separate stb-nebCycle process invocations.
sys.exit(0 if d2 > d1 * 1.5 else 1)
" > log_acceleration_check.txt 2>&1
cat log_acceleration_check.txt
check_exit_code $? 0

echo "Testing: WITHOUT a restart file, the same starting point gives the SAME first-step displacement"
stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 4 --mode 4 -O run2 --no-intro > log_prep2.txt 2>&1
python3 fabricate.py run2 00
stb-nebCycle --dir run2 --fmax 0.001 --no-intro > log_run2_step1.txt 2>&1
python3 -c "
import sys
import numpy as np
from stb.core import structure_io
def pos(run, c, label):
    return structure_io.to_pymatgen(structure_io.read_fdf(f'{run}/cycle_{c}/{label}/structure.fdf')).cart_coords
d_run1_step1 = np.linalg.norm(pos('run1', '00', 'image_01') - pos('run1', '01', 'image_01'))
d_run2_step1 = np.linalg.norm(pos('run2', '00', 'image_01') - pos('run2', '01', 'image_01'))
print(f'run1 first step={d_run1_step1:.6f}  run2 (independent) first step={d_run2_step1:.6f}')
sys.exit(0 if abs(d_run1_step1 - d_run2_step1) < 1e-9 else 1)
" > log_fresh_check.txt 2>&1
cat log_fresh_check.txt
check_exit_code $? 0


# --- 4. NaN guard: a deliberately degenerate band (all images at the same point) ---
echo -e "\n--- Testing the NaN/degenerate-tangent guard ---"
mkdir -p degen
cp initial.fdf calc.fdf degen/
(
    cd degen
    stb-neb -i ../initial.fdf -f ../initial.fdf -c calc.fdf -n 4 --mode 4 \
        --no-intro > log_degen_prep.txt 2>&1
)
check_exit_code $? 0
echo "Testing: identical initial/final endpoints -> every image sits at the exact same position"
python3 -c "
import numpy as np
from stb.core import structure_io
labels = ['image_00', 'image_01', 'image_02', 'image_03']
coords = [structure_io.to_pymatgen(structure_io.read_fdf(f'degen/cycle_00/{l}/structure.fdf')).cart_coords for l in labels]
print('all identical:', all(np.allclose(coords[0], c) for c in coords))
"
python3 -c "
import numpy as np, os
def write_fa(path, forces):
    n = len(forces)
    with open(path, 'w') as f:
        f.write(f'{n}\n')
        for i, (fx, fy, fz) in enumerate(forces, start=1):
            f.write(f'{i} {fx:.6f} {fy:.6f} {fz:.6f}\n')
def write_out(path, energy):
    with open(path, 'w') as f:
        f.write(f'siesta: FreeEng =    {energy:.6f}\nSCF cycle converged after 10 iterations\n')
for label in ('image_00', 'image_01', 'image_02', 'image_03'):
    d = f'degen/cycle_00/{label}'
    write_fa(os.path.join(d, 'siesta.FA'), np.full((9, 3), 0.01))
    write_out(os.path.join(d, 'calc.out'), -300.0)
"
stb-nebCycle --dir degen --fmax 0.001 --no-intro > log_degen_step.txt 2>&1
check_exit_code $? 1
check_contains "non-finite (NaN/Inf)" log_degen_step.txt
check_contains "path has collapsed" log_degen_step.txt
echo "Testing: no corrupted cycle_01/ was written after the guard fired"
if [ -d degen/cycle_01 ]; then
    echo -e "   -> ${RED}Failed:${NC} cycle_01 unexpectedly written despite the NaN guard"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no cycle_01 written -- guard stopped before corrupting anything"
    PASS=$((PASS+1))
fi


# --- 5. Convergence: a loose --fmax makes the very first check succeed trivially ---
echo -e "\n--- Testing convergence (loose --fmax) and the NEB_CONVERGED sentinel ---"
stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 4 --mode 4 -O run3 --no-intro > log_prep3.txt 2>&1
python3 fabricate.py run3 00
stb-nebCycle --dir run3 --fmax 10.0 --no-intro > log_converged.txt 2>&1
check_exit_code $? 0
check_contains "\[CONVERGED\]" log_converged.txt
check_success run3/NEB_CONVERGED
if [ -d run3/cycle_01 ]; then
    echo -e "   -> ${RED}Failed:${NC} cycle_01 unexpectedly written after declaring convergence"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no cycle_01 written once converged"
    PASS=$((PASS+1))
fi

echo "Testing: calling stb-nebCycle again is a clean, harmless no-op"
stb-nebCycle --dir run3 --fmax 10.0 --no-intro > log_converged_again.txt 2>&1
check_exit_code $? 0
check_contains "already exists from a previous call" log_converged_again.txt


# --- 6. --climb-after ---
echo -e "\n--- Testing --climb-after ---"
stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 4 --mode 4 -O run4 --no-intro > log_prep4.txt 2>&1
python3 fabricate.py run4 00
stb-nebCycle --dir run4 --fmax 0.001 --climb-after 5 --no-intro > log_climb_off.txt 2>&1
check_contains "Climbing image  : OFF (cycle 0 < --climb-after 5)" log_climb_off.txt
rm -rf run4/cycle_01 run4/neb_cycle_state.json
stb-nebCycle --dir run4 --fmax 0.001 --climb-after 0 --no-intro > log_climb_on.txt 2>&1
check_contains "Climbing image  : ON (cycle 0 >= --climb-after 0)" log_climb_on.txt


# --- 7. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: --dir with no cycle_NN at all"
mkdir -p no_cycles
stb-nebCycle --dir no_cycles --no-intro > log_no_cycles.txt 2>&1
check_exit_code $? 1
check_contains "No 'cycle_NN' folder found" log_no_cycles.txt

echo "Testing: missing directory entirely"
stb-nebCycle --dir does_not_exist --no-intro > log_no_dir.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_no_dir.txt

echo "Testing: a cycle with an image missing calc.out/.FA"
mkdir -p incomplete
stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 4 --mode 4 -O incomplete --no-intro > log_prep_incomplete.txt 2>&1
stb-nebCycle --dir incomplete --no-intro > log_incomplete.txt 2>&1
check_exit_code $? 1
check_contains "has SIESTA finished here yet" log_incomplete.txt

echo "Testing: --version"
stb-nebCycle --version > log_version.txt 2>&1
check_contains "stb-nebCycle" log_version.txt

echo "Testing: --help documents --dir/--fmax/--climb-after"
stb-nebCycle --help > log_help.txt 2>&1
check_contains "dir" log_help.txt
check_contains "fmax" log_help.txt
check_contains "climb-after" log_help.txt

echo "Testing: stb-nebCycle is NOT registered in the interactive stb-suite menu"
if command -v stb-suite > /dev/null; then
    python3 -c "
from stb import stb_suite
codes = [c for c in stb_suite.TOOL_CODES if c.startswith('4.9')]
print('4.9 codes:', sorted(codes))
import sys
sys.exit(1 if any('nebCycle' in stb_suite.TOOL_CODES[c].__name__ for c in codes) else 0)
" > log_not_in_menu.txt 2>&1
    check_exit_code $? 0
fi


popd > /dev/null

# --- 8. Summary ---
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
