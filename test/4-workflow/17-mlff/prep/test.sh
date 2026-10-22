#!/bin/bash

# --- Setup ---
# Smoke test for stb-mlff (Custom ML Force Field Prep, item 4.17.1)
#
# calc.fdf is a real O2-in-vacuum SIESTA input (SZ basis, same reference
# used to verify this whole workflow live: 8 rattled configs -> real SIESTA
# -> stb-mlffAnalysis -> a genuine fine-tuned MACE model). Pseudopotential
# content doesn't matter here (prep only copies the file, never reads it),
# so a placeholder is enough -- same convention as test/4-workflow/3-cohesive/prep.
#
# Layout checked here: the training configs go to <run>/configs_round_01/
# (the round they train), the test set to <run>/mlff_test/.
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$FIXTURE_DIR/test_files"
export MPLBACKEND=Agg
R1=configs_round_01

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

# --- 1. Preparation ---
echo "--- Starting tester for stb-mlff (item 4.17.1) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/calc.fdf" "$TEST_DIR/"
echo "dummy pseudopotential content" > "$TEST_DIR/O.psf"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Basic rattling run ---
echo -e "\n--- Testing a basic rattling run (--n-configs 6 --stdev 0.05 0.1) ---"
stb-mlff -o . --file calc.fdf --n-configs 6 --stdev 0.05 0.1 --pseudo-dir . --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "Generated 6 rattled configuration(s)" log_basic.txt
for i in 001 002 003 004 005 006; do
    check_success "$R1/mlff_config_$i/calc.fdf"
    check_success "$R1/mlff_config_$i/O.psf"
done

echo "Testing: each config has the same species/basis/SCF settings but different positions"
check_contains "PAO.BasisSize       SZ" $R1/mlff_config_001/calc.fdf
check_python "rattled configs have distinct, displaced positions" "
import numpy as np, sys, re
def read_positions(path):
    text = open(path).read()
    block = re.search(r'%block AtomicCoordinatesAndAtomicSpecies\n(.*?)%endblock', text, re.S).group(1)
    return np.array([[float(x) for x in line.split()[:3]] for line in block.strip().splitlines()])
p1 = read_positions('$R1/mlff_config_001/calc.fdf')
p2 = read_positions('$R1/mlff_config_002/calc.fdf')
ref = read_positions('calc.fdf')
sys.exit(0 if not np.allclose(p1, p2) and not np.allclose(p1, ref) else 1)
"

echo "Testing: every folder forces a single-point run through config_extra.fdf"
check_success $R1/mlff_config_001/config_extra.fdf
check_contains "MD.Steps              0" $R1/mlff_config_001/config_extra.fdf
check_contains "MD.VariableCell       false" $R1/mlff_config_001/config_extra.fdf
if [ "$(head -1 $R1/mlff_config_001/calc.fdf)" = "%include config_extra.fdf" ]; then
    echo -e "   -> ${GREEN}Verified:${NC} '%include config_extra.fdf' is the first line of calc.fdf"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} '%include config_extra.fdf' is not the first line of calc.fdf"
    FAIL=$((FAIL+1))
fi

echo "Testing: the report -- numbered sections, warnings last, where to run SIESTA"
check_contains "\[5\] CONFIGURATIONS" log_basic.txt
check_contains "\[6\] COVERAGE & CHECKS" log_basic.txt
check_contains "\[9\] WARNINGS" log_basic.txt
check_contains "Run SIESTA in     : ./$R1/ (6 folders)" log_basic.txt
check_contains "SIESTA: run it in ./$R1/ (6 folders) before Stage 2" log_basic.txt
check_contains "for d in ./$R1/mlff_config_\*" log_basic.txt
check_contains "Status            : OK (" log_basic.txt
check_python "the SIESTA: notice is listed in the last section, after [9] WARNINGS" "
import sys
t = open('log_basic.txt').read()
sys.exit(0 if t.rfind('SIESTA: run it in') > t.find('[9] WARNINGS') > 0 else 1)
"
check_success coverage/mlff_coverage.png
check_success coverage/coverage_configs.gplot
check_contains "Coverage figure   : ./coverage/mlff_coverage.png" log_basic.txt


# --- 2b. --volume-strains ---
echo -e "\n--- Testing --volume-strains (rattled + 2 strained configs) ---"
rm -rf $R1
stb-mlff -o . --file calc.fdf --n-configs 2 --stdev 0.05 --volume-strains -3 3 --pseudo-dir . \
    --no-intro > log_strain.txt 2>&1
check_exit_code $? 0
check_contains "Generated 2 strained configuration(s)" log_strain.txt
check_success $R1/mlff_config_004/calc.fdf
echo "Testing: strained cells have volume (1 +/- 0.03) x reference (a molecule box: every axis is scaled)"
check_python "strained volumes are 0.97 / 1.03 x reference" "
import numpy as np, sys
from stb.core import structure_io
v = lambda p: abs(np.linalg.det(structure_io.read_fdf(p).lattice))
v0 = v('calc.fdf')
ok = (np.isclose(v('$R1/mlff_config_001/calc.fdf'), v0) and np.isclose(v('$R1/mlff_config_002/calc.fdf'), v0)
      and np.isclose(v('$R1/mlff_config_003/calc.fdf'), 0.97 * v0) and np.isclose(v('$R1/mlff_config_004/calc.fdf'), 1.03 * v0))
sys.exit(0 if ok else 1)
"
stb-mlff -o . --file calc.fdf --n-configs 1 --volume-strains -100 --pseudo-dir . --no-intro > log_strain_bad.txt 2>&1
check_exit_code $? 1
check_contains "must be > -100" log_strain_bad.txt
rm -rf $R1


# --- 2c. --n-test + the workflow file ---
echo -e "\n--- Testing --n-test (held-out test set) and mlff_workflow.json ---"
rm -rf $R1 mlff_test mlff_workflow.json*
stb-mlff -o . --file calc.fdf --n-configs 4 --stdev 0.05 0.1 --volume-strains -3 3 --n-test 5 \
    --pseudo-dir . --no-intro > log_ntest.txt 2>&1
check_exit_code $? 0
check_contains "Generated 5 held-out test configuration(s)" log_ntest.txt
check_contains "Configurations    : 6 training + 5 test" log_ntest.txt
check_contains "Run SIESTA in     : ./$R1/ (6 folders) and ./mlff_test/ (5 folders)" log_ntest.txt
check_contains "Test-set E(V)     : no -- Stage 3 needs 5+ --volume-strains" log_ntest.txt
check_success mlff_test/test_001/calc.fdf
check_success mlff_test/test_005/calc.fdf
check_success mlff_test/test_005/config_extra.fdf
check_success mlff_workflow.json
echo "Testing: workflow file lists train/test configs with kinds; test kinds cycle rattle/strain; tests differ from training"
check_python "6 train + 5 test entries, kinds/volumes right, test geometries independent, per-config geometry stored" "
import json, sys, numpy as np
from stb.core import structure_io
d = json.load(open('mlff_workflow.json'))
train = [c for c in d['configs'] if c['set'] == 'train']
test = [c for c in d['configs'] if c['set'] == 'test']
kinds_test = [c['kind'] for c in test]
pos = lambda p: np.array([a[1] for a in structure_io.read_fdf(p + '/calc.fdf').atoms])
v = lambda p: abs(np.linalg.det(structure_io.read_fdf(p + '/calc.fdf').lattice))
v0 = abs(np.linalg.det(structure_io.read_fdf('calc.fdf').lattice))
distinct = all(not np.allclose(pos(t['folder']), pos(r['folder'])) for t in test for r in train)
ok = (d['format'] == 'stb-mlff-workflow' and d['reference']['file'] == 'calc.fdf'
      and d['reference']['pseudo_dir'] == '.' and len(train) == 6 and len(test) == 5
      and all(c['folder'].startswith('$R1/mlff_config_') for c in train)
      and [c['kind'] for c in train] == ['rattle'] * 4 + ['strain'] * 2
      and kinds_test == ['rattle', 'rattle', 'strain', 'strain', 'rattle']
      and np.isclose(v('mlff_test/test_003'), 0.97 * v0) and np.isclose(v('mlff_test/test_004'), 1.03 * v0)
      and distinct and d['rounds'] == [] and d['history'][0]['stage'] == 'prep'
      and all({'volume', 'dV_pct', 'u_rms', 'd_min'} <= set(c['geometry']) for c in d['configs']))
sys.exit(0 if ok else 1)
"
echo "Testing: a second run keeps the old workflow file as .bak and flags the folders it did not rewrite"
stb-mlff -o . --file calc.fdf --n-configs 1 --pseudo-dir . --no-intro > log_ntest2.txt 2>&1
check_contains "already existed -- copied to" log_ntest2.txt
check_success mlff_workflow.json.bak
check_contains "were not rewritten and are not in the new workflow file" log_ntest2.txt
echo "Testing: --n-test negative is rejected"
stb-mlff -o . --file calc.fdf --n-configs 1 --n-test -1 --pseudo-dir . --no-intro > log_ntest_bad.txt 2>&1
check_exit_code $? 1
check_contains "--n-test must be >= 0" log_ntest_bad.txt
rm -rf $R1 mlff_test mlff_workflow.json*

echo "Testing: 5 volume strains and a large enough --n-test give Stage 3 a test-set E(V)"
stb-mlff -o eos_plan --file calc.fdf --n-configs 2 --stdev 0.05 --volume-strains -4 -2 2 4 6 --n-test 6 \
    --pseudo-dir . --no-plots --no-intro > log_eos_plan.txt 2>&1
check_contains "Test-set E(V)     : yes" log_eos_plan.txt
check_not_contains "Coverage figure" log_eos_plan.txt
rm -rf eos_plan


# --- 2c'. --d3: dispersion in every folder, recorded for the later stages ---
echo "Testing: --d3 writes DFTD3 .true. into every config_extra.fdf and records it in the workflow file"
stb-mlff -o d3_run --file calc.fdf --n-configs 2 --stdev 0.05 --n-test 1 --pseudo-dir . --d3 --no-plots \
    --no-intro > log_d3.txt 2>&1
check_exit_code $? 0
check_contains "DFT-D3 forced in every folder" log_d3.txt
check_contains "and DFT-D3 (--d3)" log_d3.txt
check_contains "DFTD3                 .true." d3_run/$R1/mlff_config_001/config_extra.fdf
check_contains "DFTD3                 .true." d3_run/mlff_test/test_001/config_extra.fdf
check_python "reference.dftd3 = on, forced, source --d3" "
import json, sys
r = json.load(open('d3_run/mlff_workflow.json'))['reference']['dftd3']
sys.exit(0 if r['on'] and r['forced'] and r['source'].startswith('--d3') else 1)
"
stb-mlff -o nod3_run --file calc.fdf --n-configs 1 --pseudo-dir . --no-plots --no-intro > log_nod3.txt 2>&1
check_contains "Dispersion           : none" log_nod3.txt
check_not_contains "DFTD3" nod3_run/$R1/mlff_config_001/config_extra.fdf
rm -rf d3_run nod3_run

# --- 2d. Extra outputs: --export-traj, near-duplicates, --show without a display ---
echo -e "\n--- Testing --export-traj, --dup-threshold and --show (no display) ---"
stb-mlff -o extras --file calc.fdf --n-configs 3 --stdev 0.05 --n-test 1 --pseudo-dir . --export-traj \
    --dup-threshold 1000 --show --no-intro > log_extras.txt 2>&1
check_exit_code $? 0
check_success extras/mlff_configs.xsf
check_python "the trajectory holds the reference + every config (1 + 3 + 1 frames)" "
import sys
from ase.io import read
sys.exit(0 if len(read('extras/mlff_configs.xsf', ':')) == 5 else 1)
"
check_contains "near-duplicate config(s)" log_extras.txt
check_contains "--show needs a display" log_extras.txt
rm -rf extras


# --- 2e. Foundation-model screen (one shared 'unphysical' threshold) ---
echo -e "\n--- Testing --ml-screen (MACE-MP small) ---"
stb-mlff -o mlscreen --file calc.fdf --n-configs 2 --stdev 0.05 0.6 --pseudo-dir . --ml-screen --ml-model small \
    --no-plots --no-intro > log_mlscreen.txt 2>&1
check_exit_code $? 0
check_contains "dE_ML" log_mlscreen.txt
check_contains "Fmax_ML" log_mlscreen.txt
check_not_contains "cuequivariance" log_mlscreen.txt
check_contains "Library messages  :" log_mlscreen.txt
check_python "ml_screen numbers stored in the workflow file" "
import json, sys
d = json.load(open('mlscreen/mlff_workflow.json'))
sys.exit(0 if all('f_max' in c['ml_screen'] and 'dE_meV_atom' in c['ml_screen'] for c in d['configs']) else 1)
"
check_python "Stage 1's unphysical threshold is the workflow's shared one (15 eV/Ang, 1000 meV/atom)" "
import sys
from stb import mlff
from stb.core import mlff_workflow as w
sys.exit(0 if (mlff.ML_FORCE_WARN, mlff.ML_ENERGY_WARN) == (w.UNPHYSICAL_FMAX, w.UNPHYSICAL_DE) == (15.0, 1000.0) else 1)
"
rm -rf mlscreen


# --- 3. --n-configs 0 with no --from-aimd (expects failure) ---
echo -e "\n--- Testing --n-configs 0 with no --from-aimd (expects failure) ---"
stb-mlff -o . --file calc.fdf --n-configs 0 --pseudo-dir . --no-intro > log_empty.txt 2>&1
check_exit_code $? 1
check_contains "No configurations were generated" log_empty.txt


# --- 4. Missing pseudopotential (expects failure) ---
echo -e "\n--- Testing a missing pseudopotential (expects failure) ---"
mkdir -p empty_pp_dir
stb-mlff -o . --file calc.fdf --n-configs 2 --pseudo-dir empty_pp_dir --no-intro > log_missingpp.txt 2>&1
check_exit_code $? 1
check_contains "Missing pseudopotential" log_missingpp.txt


# --- 5. --save-report (written in the run folder) ---
echo -e "\n--- Testing --save-report ---"
rm -f stb_mlff_report.txt
stb-mlff -o report_run --file calc.fdf --n-configs 2 --pseudo-dir . --save-report --no-intro > log_savereport.txt 2>&1
check_success report_run/stb_mlff_report.txt
check_contains "STB-MLFF REPORT" report_run/stb_mlff_report.txt
check_contains "\[9\] WARNINGS" report_run/stb_mlff_report.txt
rm -rf report_run $R1


# --- 6. Robustness / errors ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: no input -> default structure.fdf missing"
stb-mlff --no-intro > log_missing_file.txt 2>&1
check_exit_code $? 1
check_contains "Structure file 'structure.fdf' not found" log_missing_file.txt

echo "Testing: nonexistent reference file"
stb-mlff --file nope.fdf --no-intro > log_missing_ref.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_ref.txt

echo "Testing: --version"
stb-mlff --version > log_version.txt 2>&1
check_contains "stb-mlff" log_version.txt

echo "Testing: --help documents the options"
stb-mlff --help > log_help.txt 2>&1
for opt in stdev from-aimd n-samples volume-strains n-test workflow ml-screen export-traj dup-threshold show configs_round_01; do
    check_contains "$opt" log_help.txt
done


# --- 6b. structure.fdf + calc.fdf pair (-s / -c, the suite's usual inputs) ---
echo -e "\n--- Testing -s structure.fdf -c calc.fdf ---"
mkdir -p pair && cp O.psf pair/
# Split the self-contained fixture: structure blocks -> structure.fdf, the rest
# -> calc.fdf with '%include structure.fdf'.
python3 - <<'PY'
lines = open('calc.fdf').read().splitlines(keepends=True)
start = next(i for i, l in enumerate(lines) if l.startswith('NumberOfSpecies'))
end = next(i for i, l in enumerate(lines) if l.startswith('%endblock AtomicCoordinates')) + 1
struct = [l for l in lines[start:end] if not l.startswith(('PAO.', 'MeshCutoff'))]
calc = lines[:start] + ['%include structure.fdf\n'] + [l for l in lines[start:end] if l.startswith(('PAO.', 'MeshCutoff'))] + lines[end:]
open('pair/structure.fdf', 'w').writelines(struct)
open('pair/calc.fdf', 'w').writelines(calc)
open('pair/calc_noinclude.fdf', 'w').writelines(lines[:start] + lines[end:])
PY
pushd pair > /dev/null
stb-mlff -o . -s structure.fdf -c calc.fdf --n-configs 2 --volume-strains 3 --n-test 2 --pseudo-dir . \
    --no-intro > log_pair.txt 2>&1
check_exit_code $? 0
check_contains "Structure source  : structure.fdf" log_pair.txt
check_contains "Calc template     : calc.fdf (%include structure.fdf)" log_pair.txt
check_success $R1/mlff_config_001/structure.fdf
check_success $R1/mlff_config_001/calc.fdf
check_success mlff_test/test_002/structure.fdf
check_python "calc.fdf copied with config_extra on top and its own %include kept; structure.fdf has new positions" "
import sys
calc = open('$R1/mlff_config_001/calc.fdf').read().splitlines()
ref = open('calc.fdf').read().splitlines()
st = open('$R1/mlff_config_001/structure.fdf').read()
sys.exit(0 if calc[0] == '%include config_extra.fdf' and calc[2:] == ref
         and '%include' not in st and st != open('structure.fdf').read() else 1)
"
check_python "strained config rewrites the lattice in structure.fdf, calc.fdf untouched" "
import sys
from stb.core import structure_io as s
a = s.read_fdf('$R1/mlff_config_003/structure.fdf').lattice[0][0]
sys.exit(0 if abs(a - 10.0 * 1.03 ** (1 / 3)) < 1e-6 and
         open('$R1/mlff_config_003/calc.fdf').read().split('\n', 2)[2] == open('calc.fdf').read() else 1)
"
check_python "workflow records the pair (file = structure, calc = template)" "
import json, sys
r = json.load(open('mlff_workflow.json'))['reference']
sys.exit(0 if r['file'] == 'structure.fdf' and r['calc'] == 'calc.fdf' else 1)
"
echo "Testing: calc.fdf without %include structure.fdf -> clear error"
stb-mlff -o . -s structure.fdf -c calc_noinclude.fdf --n-configs 1 --pseudo-dir . --no-intro > log_pair_noinc.txt 2>&1
check_exit_code $? 1
check_contains "has no '%include structure.fdf'" log_pair_noinc.txt
echo "Testing: bundled bank --pseudo-dir virtual_vault (resolved, copied, recorded by name; no banner mid-report)"
stb-mlff -o bank_test -s structure.fdf -c calc.fdf --n-configs 1 --pseudo-dir virtual_vault --no-intro > log_pair_bank.txt 2>&1
check_exit_code $? 0
check_contains "bundled bank 'virtual_vault' -> O.psf" log_pair_bank.txt
check_not_contains "\[INFO\] Using bundled pseudopotential bank" log_pair_bank.txt
check_python "O.psf copied from the bank, bank name stored in the workflow" "
import json, sys, os
r = json.load(open('bank_test/mlff_workflow.json'))['reference']
sys.exit(0 if r['pseudo_dir'] == 'virtual_vault' and os.path.getsize('bank_test/$R1/mlff_config_001/O.psf') > 1000 else 1)
"
rm -rf bank_test
echo "Testing: unknown pseudo source -> clear error"
stb-mlff -o . -s structure.fdf -c calc.fdf --n-configs 1 --pseudo-dir no_such_bank --no-intro > log_pair_badbank.txt 2>&1
check_exit_code $? 1
check_contains "not a recognized pseudopotential bank" log_pair_badbank.txt
echo "Testing: missing calc.fdf -> clear error"
stb-mlff -o . -s structure.fdf -c nope.fdf --n-configs 1 --pseudo-dir . --no-intro > log_pair_nocalc.txt 2>&1
check_exit_code $? 1
check_contains "Calc file 'nope.fdf' not found" log_pair_nocalc.txt
echo "Testing: defaults (structure.fdf + calc.fdf here, everything written to mlff_runs/)"
rm -rf $R1 mlff_test mlff_workflow.json
stb-mlff --n-configs 1 --n-test 1 --pseudo-dir . --no-intro > log_pair_default.txt 2>&1
check_exit_code $? 0
check_success mlff_runs/$R1/mlff_config_001/structure.fdf
check_success mlff_runs/$R1/mlff_config_001/calc.fdf
check_success mlff_runs/mlff_test/test_001/calc.fdf
check_success mlff_runs/mlff_workflow.json
check_contains "Run SIESTA in each mlff_runs/$R1/mlff_config_NNN/" log_pair_default.txt
check_python "workflow inside mlff_runs/ points back to ../structure.fdf, ../calc.fdf and the pseudos in .." "
import json, sys
d = json.load(open('mlff_runs/mlff_workflow.json'))
r = d['reference']
sys.exit(0 if r['file'] == '../structure.fdf' and r['calc'] == '../calc.fdf' and r['pseudo_dir'] == '..'
         and [c['folder'] for c in d['configs']] == ['$R1/mlff_config_001', 'mlff_test/test_001'] else 1)
"
check_python "the later stages find mlff_runs/mlff_workflow.json by default, and ./mlff_workflow.json from inside it" "
import os, sys
from stb.core import mlff_workflow as wf
a = wf.default_path()
os.chdir('mlff_runs'); b = wf.default_path()
sys.exit(0 if a == os.path.join('mlff_runs', 'mlff_workflow.json') and b == 'mlff_workflow.json' else 1)
"
check_python "the per-round helpers: folder, project-wide numbering, SIESTA status" "
import sys
from stb.core import mlff_workflow as w
ok = (w.configs_dir('mlff_runs', 2) == 'mlff_runs/configs_round_02'
      and w.find_config_dirs('mlff_runs') == ['mlff_runs/$R1/mlff_config_001']
      and w.next_config_index('mlff_runs') == 2
      and w.pending_by_dir(w.find_config_dirs('mlff_runs')) == [('mlff_runs/$R1', 1, 1)])
sys.exit(0 if ok else 1)
"
popd > /dev/null


# --- 6c. A 2D cell: volume strains change the periodic axes only ---
echo -e "\n--- Testing --volume-strains on a 2D slab (vacuum axis kept) ---"
mkdir -p slab && cp O.psf slab/C.psf
cat > slab/slab.fdf <<'FDF'
SystemLabel graphene
NumberOfSpecies 1
NumberOfAtoms 2
%block ChemicalSpeciesLabel
 1 6 C
%endblock ChemicalSpeciesLabel
LatticeConstant 1.0 Ang
%block LatticeVectors
 1.24172751 -2.15069131 0.0
 1.24172751  2.15069131 0.0
 0.0 0.0 20.0
%endblock LatticeVectors
AtomicCoordinatesFormat Fractional
%block AtomicCoordinatesAndAtomicSpecies
 0.33333333 0.66666667 0.5 1
 0.66666667 0.33333333 0.5 1
%endblock AtomicCoordinatesAndAtomicSpecies
PAO.BasisSize SZ
FDF
pushd slab > /dev/null
stb-mlff -o . --file slab.fdf --n-configs 1 --stdev 0.02 --volume-strains 4 --pseudo-dir . --no-plots \
    --no-intro > log_slab.txt 2>&1
check_exit_code $? 0
check_contains "Volume strains    : periodic axes only" log_slab.txt
check_not_contains "also scales the vacuum" log_slab.txt
check_python "the strained slab keeps c = 20 Ang and changes the area by +4 %" "
import numpy as np, sys
from stb.core import structure_io as s
c0 = s.read_fdf('slab.fdf').lattice; c1 = s.read_fdf('$R1/mlff_config_002/slab.fdf').lattice
area = lambda c: np.linalg.norm(np.cross(c[0], c[1]))
sys.exit(0 if np.isclose(c1[2][2], 20.0) and np.isclose(area(c1) / area(c0), 1.04) else 1)
"
popd > /dev/null


# --- 7. Interactive path (stb-suite, shortcut 4.17.1) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.17.1) ---"
rm -rf mlff_runs
# Structure 'calc.fdf', calc '-' = one self-contained .fdf (--file); then pseudo dir,
# 3 configs, stdev, no strains, 2 tests, no D3, no AIMD, default run folder, no report,
# no ML screen, no export, no viewer, no --show; Enter; 0 = quit.
printf '4.17.1\ncalc.fdf\n-\n\n3\n0.05\n\n2\nn\nn\n\nn\nn\n\nn\nn\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success mlff_runs/$R1/mlff_config_001/calc.fdf
check_success mlff_runs/$R1/mlff_config_003/calc.fdf
check_success mlff_runs/mlff_test/test_002/calc.fdf
check_success mlff_runs/mlff_workflow.json
check_contains "Display the coverage figure on screen" log_menu.txt

echo "Testing: menu with the default structure.fdf + calc.fdf pair and a custom run folder"
pushd pair > /dev/null
rm -rf mlff_runs
printf '4.17.1\n\n\n\n2\n0.05\n\n0\nn\nn\nmy_mlff\nn\nn\n\nn\nn\n\n0\n' | stb-suite > log_menu_pair.txt 2>&1
check_exit_code $? 0
check_success my_mlff/$R1/mlff_config_002/structure.fdf
check_success my_mlff/$R1/mlff_config_002/calc.fdf
check_success my_mlff/mlff_workflow.json
popd > /dev/null


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
