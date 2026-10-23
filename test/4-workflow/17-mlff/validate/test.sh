#!/bin/bash

# --- Setup ---
# Smoke test for stb-mlffValidate (Custom ML Force Field Validate, item 4.17.3)
#
# Reuses the 8 real finished-SIESTA fixtures from ../analysis/: 6 become the
# training set (configs_round_01/, the per-round layout), 2 the held-out test
# set (mlff_test/test_00N, listed in mlff_workflow.json as Stage 1 --n-test
# would). Two quick training rounds (2-model committees, --epochs 2) are scored
# round by round. The phonon comparison (--phonon-ref) runs against a synthetic
# stb-phononsCreate folder whose ".FA" forces come from a trained model plus a
# constant residual -- the residual must be subtracted on both sides. Needs the
# optional 'ml' extra -- skips gracefully (not a failure) if mace_run_train
# isn't on PATH.
export MPLBACKEND=Agg
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MLFF_FIXTURES="$FIXTURE_DIR/../analysis"
TEST_DIR="$FIXTURE_DIR/test_files"
R1=configs_round_01
TRAIN="--epochs 2 --batch-size 2 --device cpu --foundation-model small --no-stability --skip-foundation-comparison"

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
        echo -e "   -> ${GREEN}Verified:${NC} '$1' absent from '$2'"
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

check_python() {   # $1 = description, $2 = python code exiting 0 on success
    python3 -c "$2"
    if [ $? -eq 0 ]; then
        echo -e "   -> ${GREEN}Verified:${NC} $1"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} $1"
        FAIL=$((FAIL+1))
    fi
}


echo "--- Starting tester for stb-mlffValidate (item 4.17.3) ---"

if ! command -v mace_run_train &> /dev/null; then
    echo -e "${YELLOW}SKIPPED${NC}: 'mace_run_train' not on PATH (needs 'pip install stb_suite[ml]')."
    exit 0
fi
if [ ! -d "$MLFF_FIXTURES/mlff_config_001" ]; then
    echo -e "${YELLOW}SKIPPED${NC}: ../analysis fixtures not found."
    exit 0
fi

rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR/mlff_test" "$TEST_DIR/$R1"
for i in 1 2 3 4 5 6; do cp -r "$MLFF_FIXTURES/mlff_config_00$i" "$TEST_DIR/$R1/"; done
cp -r "$MLFF_FIXTURES/mlff_config_007" "$TEST_DIR/mlff_test/test_001"
cp -r "$MLFF_FIXTURES/mlff_config_008" "$TEST_DIR/mlff_test/test_002"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null

# A workflow file as Stage 1 would write it: 6 training + 2 test configs.
python3 -c "
from stb.core import mlff_workflow as wf
s = wf.new_state('$R1/mlff_config_001/calc.fdf', 2, ['O'], '.')
wf.add_configs(s, [{'folder': f'$R1/mlff_config_00{i}', 'set': 'train', 'kind': 'rattle', 'origin': 'fixture',
                    'stage': 'prep', 'round': 0} for i in range(1, 7)], '.')
wf.add_configs(s, [{'folder': 'mlff_test/test_001', 'set': 'test', 'kind': 'rattle', 'origin': 'fixture', 'stage': 'prep', 'round': 0},
                   {'folder': 'mlff_test/test_002', 'set': 'test', 'kind': 'strain', 'origin': 'fixture', 'stage': 'prep', 'round': 0}], '.')
wf.save(s, 'mlff_workflow.json')
"
check_success mlff_workflow.json


# --- 1. Errors before any round ---
echo -e "\n--- Testing: no training round yet ---"
stb-mlffValidate --no-intro > log_noround.txt 2>&1
check_exit_code $? 1
check_contains "No training rounds" log_noround.txt


# --- 2. Round 1 ---
echo -e "\n--- Training round 1 (2-model committee) and validating it ---"
stb-mlffAnalysis --committee 2 $TRAIN --no-intro > log_train1.txt 2>&1
check_exit_code $? 0
check_contains "Round directory   : round_01" log_train1.txt

stb-mlffValidate --save-report --save-data --foundation-model small --no-intro > log_val1.txt 2>&1
check_exit_code $? 0
for sec in "\[0\] RUN METADATA" "\[1\] TEST SET" "\[2\] MODELS SCORED" "\[3\] TEST ERRORS" "\[4\] ROUND 1 IN DETAIL" \
           "\[5\] COMMITTEE (round 1)" "\[6\] PHONONS" "\[7\] EQUATION OF STATE" "\[8\] MD CHECK" \
           "\[9\] ROUND COMPARISON" "\[10\] ROUND HISTORY, LEARNING CURVE & VERDICT" "\[11\] SUMMARY & FILES" \
           "\[12\] NEXT STEPS" "\[13\] WARNINGS"; do
    check_contains "$sec" log_val1.txt
done
echo "Testing: the best model and its card are copied to best_model/"
check_success best_model/mlff_model_compiled.model
check_success best_model/mlff_model_card.json
check_contains "Best model        : round" log_val1.txt
check_python "best_model/ holds the model with a matching hash, a selection section and the 'best' entry of the workflow file" "
import json, hashlib, sys
c = json.load(open('best_model/mlff_model_card.json'))
b = json.load(open('mlff_workflow.json'))['best']
f = {x['file']: x['sha256'] for x in c['model']['files']}
ok = (c['selection']['round'] == b['round'] == 1 and c['selection']['seed'] == b['seed'] and b['folder'] == 'best_model'
      and f['mlff_model_compiled.model'] == hashlib.sha256(open('best_model/mlff_model_compiled.model', 'rb').read()).hexdigest()
      and len(c['selection']['rounds_compared']) == 1 and c['validation'].get('overall'))
sys.exit(0 if ok else 1)
"
check_contains "Test set from     : mlff_workflow.json (test set)" log_val1.txt
check_contains "Configs           : 2 finished of 2 folder(s), 4 atoms, rattle 1, strain 1" log_val1.txt
check_contains "\[OK\] Same level of theory in all 8 folders" log_val1.txt
check_contains "Training range    :" log_val1.txt
check_contains "foundation |" log_val1.txt
check_contains "Per element" log_val1.txt
check_contains "Per config (worst first)" log_val1.txt
check_contains "Averaging         :" log_val1.txt
check_contains "\[CONTINUE\] First scored round" log_val1.txt
check_contains "Test set: only 2 config(s)" log_val1.txt
check_contains "Learning curve    : needs rounds with 3+ different training-set sizes" log_val1.txt
check_not_contains "WARNING:root" log_val1.txt
check_not_contains "cuequivariance" log_val1.txt
check_success mlff_validation.png
check_success mlff_validation_errors.dat
check_success stb_mlffValidate_report.txt
check_contains "STB-MLFFVALIDATE REPORT" stb_mlffValidate_report.txt
check_python "round 1 validation stored in mlff_workflow.json" "
import json, sys
d = json.load(open('mlff_workflow.json'))
v = d['rounds'][0]['validation']
sys.exit(0 if v['n_test'] == 2 and v['overall']['F_rmse_meV'] > 0 and len(v['members']) == 2
         and set(v['groups']) == {'rattle', 'strain'} and 'calibration' in v and v['verdict'] == 'continue'
         and 'O' in v['elements'] and len(v['per_config']) == 2
         and v['overall'].get('F_rmse_sigma_meV') is not None
         and d['history'][-1]['stage'] == 'validate' else 1)
"


# --- 3. Round 2: two rounds compared, verdicts ---
echo -e "\n--- Training round 2 and comparing the rounds ---"
stb-mlffAnalysis --committee 2 $TRAIN --no-intro > log_train2.txt 2>&1
check_exit_code $? 0
stb-mlffValidate --skip-foundation-comparison --compare 1 2 --no-intro > log_val2.txt 2>&1
check_exit_code $? 0
check_contains "\[4\] ROUND 2 IN DETAIL" log_val2.txt
check_python "history table lists rounds 1 and 2, then a CONTINUE, PLATEAU or REGRESSED verdict" "
import re, sys
t = open('log_val2.txt').read().split('[10] ROUND HISTORY', 1)[1]
rows = re.findall(r'^(?:\x1b\[\d+m)?(\d+)\s+\|', t, re.M)
sys.exit(0 if rows[:2] == ['1', '2'] and re.search(r'\[(CONTINUE|PLATEAU|REGRESSED)\] Round 2', t) else 1)
"
check_contains "Overall           : F " log_val2.txt
check_contains "better," log_val2.txt
check_contains "Per config (largest regression first)" log_val2.txt

echo "Testing: --target-force reached -> STOP"
stb-mlffValidate --target-force 1e9 --skip-foundation-comparison --no-intro > log_val_stop.txt 2>&1
check_contains "\[STOP\] Target met" log_val_stop.txt

echo "Testing: --plateau 1000 -> any change counts as a plateau"
stb-mlffValidate --plateau 1000 --skip-foundation-comparison --no-intro > log_val_plateau.txt 2>&1
check_contains "\[PLATEAU\] Round 2" log_val_plateau.txt

echo "Testing: --plateau 0 and a worse round 2 -> REGRESSED (else CONTINUE), naming the best round"
stb-mlffValidate --plateau 0 --skip-foundation-comparison --no-intro > log_val_regr.txt 2>&1
check_python "the verdict follows the rounds' test errors" "
import json, re, sys
v = [r['validation']['overall']['F_rmse_meV'] for r in json.load(open('mlff_workflow.json'))['rounds']]
t = open('log_val_regr.txt').read()
ok = (re.search(r'\[REGRESSED\] Round 2 .* than round 1', t) if v[1] > v[0] else re.search(r'\[CONTINUE\] Round 2', t))
sys.exit(0 if ok else 1)
"

echo "Testing: --rounds 1 scores only round 1; --detail-round"
stb-mlffValidate --rounds 1 --skip-foundation-comparison --no-intro > log_val_r1.txt 2>&1
check_contains "Scoring           : round 1 (2 model(s))" log_val_r1.txt
check_not_contains "round 2 (2 model(s))" log_val_r1.txt
check_contains "was chosen among the rounds scored in this run only" log_val_r1.txt
stb-mlffValidate --detail-round 1 --skip-foundation-comparison --no-intro > log_val_d1.txt 2>&1
check_contains "\[4\] ROUND 1 IN DETAIL" log_val_d1.txt
stb-mlffValidate --rounds 7 --no-intro > log_val_r7.txt 2>&1
check_exit_code $? 1
check_contains "Round 7 is not in .* (it has round(s) 1, 2)" log_val_r7.txt
stb-mlffValidate --compare 1 9 --no-intro > log_val_c9.txt 2>&1
check_exit_code $? 1

echo "Testing: --view-worst writes the worst test configs with per-atom errors (no display here)"
env -u DISPLAY -u WAYLAND_DISPLAY stb-mlffValidate --view-worst 1 --skip-foundation-comparison --no-intro > log_val_worst.txt 2>&1
check_success mlff_test_worst.xyz
check_python "worst-config file carries the per-atom F_error" "
import sys
from ase.io import read
a = read('mlff_test_worst.xyz', ':')
sys.exit(0 if len(a) == 1 and 'F_error' in a[0].arrays else 1)
"


# --- 4. Explicit models / test folders, no workflow ---
echo -e "\n--- Testing --model + --test-path (no workflow file) ---"
M1=$(python3 -c "import json; print(json.load(open('mlff_workflow.json'))['rounds'][0]['committee'][0]['model'])")
stb-mlffValidate --workflow none.json --model "$M1" --test-path 'mlff_test/test_*' \
    --skip-foundation-comparison --target-force 1e9 --no-intro > log_val_explicit.txt 2>&1
check_exit_code $? 0
check_contains "1 model(s) from --model" log_val_explicit.txt
check_contains "Configs           : 2 finished of 2 folder(s), 4 atoms, test 2" log_val_explicit.txt
check_contains "no round history to compare" log_val_explicit.txt
check_contains "\[TARGET MET\]" log_val_explicit.txt

echo "Testing: test folders that were also trained on are flagged"
stb-mlffValidate --test-path "$R1/mlff_config_00[1-2]" --rounds 1 --skip-foundation-comparison --no-intro > log_val_overlap.txt 2>&1
check_contains "were also trained on" log_val_overlap.txt

echo "Testing: no test set -> clear error"
stb-mlffValidate --test-path 'no_such_*' --no-intro > log_val_notest.txt 2>&1
check_exit_code $? 1
check_contains "No test configs found" log_val_notest.txt

echo "Testing: --phonon-ref without phonopy_disp.yaml -> clear error"
stb-mlffValidate --phonon-ref round_01 --skip-foundation-comparison --no-intro > log_val_phref.txt 2>&1
check_exit_code $? 1
check_contains "has no phonopy_disp.yaml" log_val_phref.txt


# --- 4b. --phonon-ref against a synthetic stb-phononsCreate folder ---
# A symmetric 2D O lattice (P6/mmm: ONE '+' displacement, as for graphene -- a low-symmetry
# cell gets +/- pairs, which cancel a constant residual by themselves). Its ".FA" forces are
# round 1's first model plus a staggered 0.05 eV/Ang residual: the same model on both sides,
# so the phonon RMSE is exactly 0 only if the residual is subtracted on both.
echo -e "\n--- Testing --phonon-ref (model forces + a residual, residual subtraction) ---"
mkdir -p phref && pushd phref > /dev/null
echo "dummy" > O.psf
cat > structure.fdf <<'FDF'
NumberOfSpecies 1
NumberOfAtoms 2
%block ChemicalSpeciesLabel
 1 8 O
%endblock ChemicalSpeciesLabel
LatticeConstant 1.0 Ang
%block LatticeVectors
 1.25 -2.16506351 0.0
 1.25  2.16506351 0.0
 0.0 0.0 15.0
%endblock LatticeVectors
AtomicCoordinatesFormat Fractional
%block AtomicCoordinatesAndAtomicSpecies
 0.33333333 0.66666667 0.5 1
 0.66666667 0.33333333 0.5 1
%endblock AtomicCoordinatesAndAtomicSpecies
FDF
printf 'SystemLabel calc\n%%include structure.fdf\nPAO.BasisSize SZ\nMeshCutoff 150 Ry\n' > calc.fdf
stb-phononsCreate --no-intro -dim 2 2 1 -s structure.fdf -c calc.fdf -p . > log_create.txt 2>&1
check_contains "Displacements needed         | 1 " log_create.txt
check_success phonon_runs/perfect_supercell/structure.fdf
python3 - "../$M1" <<'PY'
import glob, sys, warnings
import numpy as np
warnings.filterwarnings("ignore")
from pymatgen.io.ase import AseAtomsAdaptor
from mace.calculators import MACECalculator
from stb.core import structure_io as io
calc = MACECalculator(model_paths=[sys.argv[1]], device="cpu", default_dtype="float64")
for d in sorted(glob.glob("phonon_runs/disp-*")) + ["phonon_runs/perfect_supercell"]:
    a = AseAtomsAdaptor.get_atoms(io.to_pymatgen(io.read_fdf(f"{d}/structure.fdf")))   # written in Bohr
    a.calc = calc
    resid = np.zeros((len(a), 3))
    resid[:, 1] = 0.05 * (-1) ** np.arange(len(a))
    f = a.get_forces() + resid
    with open(f"{d}/calc.FA", "w") as fh:
        fh.write(f"{len(f)}\n" + "".join(f"{i + 1:6d} {x:14.8f} {y:14.8f} {z:14.8f}\n" for i, (x, y, z) in enumerate(f)))
PY
popd > /dev/null
stb-mlffValidate --workflow none.json --model "$M1" --test-path 'mlff_test/test_*' --phonon-ref phref/phonon_runs \
    --phonon-mesh 2 2 1 --skip-foundation-comparison --no-intro > log_val_phonon.txt 2>&1
check_exit_code $? 0
check_contains "\[6\] PHONONS vs SIESTA" log_val_phonon.txt
check_contains "Residual forces  : subtracted (undisplaced supercell, max 0.0500" log_val_phonon.txt
check_python "the same model on both sides, residual subtracted on both: phonon RMSE 0.00 THz" "
import re, sys
t = open('log_val_phonon.txt').read().split('[6] PHONONS vs SIESTA', 1)[1]
m = re.search(r'--model \| ([\d.]+) ', t)
sys.exit(0 if m and float(m.group(1)) == 0.0 else 1)
"
echo "Testing: without the undisplaced supercell's forces the residual stays in (and the report says so)"
rm phref/phonon_runs/perfect_supercell/calc.FA
stb-mlffValidate --workflow none.json --model "$M1" --test-path 'mlff_test/test_*' --phonon-ref phref/phonon_runs \
    --phonon-mesh 2 2 1 --skip-foundation-comparison --no-intro > log_val_phonon_raw.txt 2>&1
check_contains "Residual forces  : NOT subtracted" log_val_phonon_raw.txt
check_contains "perfect_supercell/ has no finished SIESTA run" log_val_phonon_raw.txt
check_python "phonon RMSE well above 0 once the residual is left in (several THz)" "
import re, sys
t = open('log_val_phonon_raw.txt').read().split('[6] PHONONS vs SIESTA', 1)[1]
m = re.search(r'--model \| ([\d.]+) ', t)
sys.exit(0 if m and float(m.group(1)) > 1.0 else 1)
"
check_python "the phonon loader uses symprec 0.01, like stb-phononsCreate" "
import sys
from stb import mlff_validate
sys.exit(0 if mlff_validate.PHONON_SYMPREC == 0.01 else 1)
"


# --- 5. Help / version ---
echo -e "\n--- Testing --help / --version ---"
stb-mlffValidate --version > log_version.txt 2>&1
check_contains "stb-mlffValidate" log_version.txt
stb-mlffValidate --help > log_help.txt 2>&1
for opt in workflow rounds model test-path phonon-ref phonon-mesh target-force plateau detail-round compare \
           md-check md-check-to-train show save-data view-worst; do
    check_contains "$opt" log_help.txt
done


# --- 6. Interactive path (stb-suite, shortcut 4.17.3) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.17.3) ---"
# target, phonon ref, device, show, data file, worst configs, compare rounds, MD check, report
printf '4.17.3\n\n\n\nn\nn\n0\n1 2\nn\nn\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_contains "mlff_workflow.json: 6 training + 2 test config(s), 2 round(s)" log_menu.txt
check_contains "two of 1, 2, e.g. '1 2'" log_menu.txt
check_contains "\[9\] ROUND COMPARISON" log_menu.txt
check_contains "VERDICT" log_menu.txt


# --- 7. Test folders on disk vs the workflow file ---
echo -e "\n--- Testing a test folder missing from the workflow file and one without SIESTA ---"
cp -r mlff_test/test_002 mlff_test/test_003                 # finished, not listed
mkdir -p mlff_test/test_004 && cp mlff_test/test_001/calc.fdf mlff_test/test_004/   # listed, never run
python3 - <<'PY'
import json
d = json.load(open('mlff_workflow.json'))
d['configs'].append({'folder': 'mlff_test/test_004', 'set': 'test', 'kind': 'rattle', 'origin': 'x',
                     'stage': 'prep', 'round': 0})
json.dump(d, open('mlff_workflow.json', 'w'), indent=1)
PY
stb-mlffValidate --skip-foundation-comparison --no-intro > log_val_disk.txt 2>&1
check_exit_code $? 0
check_contains "1 folder(s) in mlff_test/ are not in the workflow file and were not scored (test_003)" log_val_disk.txt
check_contains "SIESTA: 1 of 3 test folder(s) in mlff_test/ have no finished run" log_val_disk.txt


popd > /dev/null

# --- 7b. Learning-curve sign, configs far above the test set in energy ---
echo "Testing: a growing learning curve is written N^+b; far_training_configs flags high-energy training configs"
check_python "power law printed with one sign; far configs = more than one test spread above the top" "
import sys
from ase import Atoms
from stb.mlff_validate import learning_curve, _power_law, far_training_configs
lc = learning_curve([(14, 21.3), (19, 31.1), (24, 37.5)])
up = _power_law(lc)
down = _power_law(learning_curve([(10, 40.0), (20, 28.0), (40, 20.0)]))
def cfg(name, e):
    a = Atoms('C2', positions=[[0, 0, 0], [1.4, 0, 0]]); a.info.update(REF_energy=2 * e, folder='/x/' + name); return a
test = [cfg('t1', 0.0), cfg('t2', 0.3)]
far = far_training_configs([cfg('a', 0.1), cfg('b', 0.5), cfg('c', 1.5)], test)
sys.exit(0 if lc['b'] < 0 and '^+' in up and '--' not in up and '^-' in down
         and [n for n, _ in far] == ['c'] and abs(far[0][1] - 1200) < 1e-6 else 1)
"

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
