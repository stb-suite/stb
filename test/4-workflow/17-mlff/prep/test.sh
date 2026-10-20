#!/bin/bash

# --- Setup ---
# Smoke test for stb-mlff (Custom ML Force Field Prep, item 4.17.1)
#
# calc.fdf is a real O2-in-vacuum SIESTA input (SZ basis, same reference
# used to verify this whole workflow live: 8 rattled configs -> real SIESTA
# -> stb-mlffAnalysis -> a genuine fine-tuned MACE model). Pseudopotential
# content doesn't matter here (prep only copies the file, never reads it),
# so a placeholder is enough -- same convention as test/4-workflow/3-cohesive/prep.
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
stb-mlff --file calc.fdf --n-configs 6 --stdev 0.05 0.1 --pseudo-dir . --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "Generated 6 rattled configuration(s)" log_basic.txt
for i in 001 002 003 004 005 006; do
    check_success "mlff_config_$i/calc.fdf"
    check_success "mlff_config_$i/O.psf"
done

echo "Testing: each config has the same species/basis/SCF settings but different positions"
check_contains "PAO.BasisSize       SZ" mlff_config_001/calc.fdf
python3 -c "
import numpy as np, sys, re
def read_positions(path):
    text = open(path).read()
    block = re.search(r'%block AtomicCoordinatesAndAtomicSpecies\n(.*?)%endblock', text, re.S).group(1)
    return np.array([[float(x) for x in line.split()[:3]] for line in block.strip().splitlines()])
p1 = read_positions('mlff_config_001/calc.fdf')
p2 = read_positions('mlff_config_002/calc.fdf')
ref = read_positions('calc.fdf')
sys.exit(0 if not np.allclose(p1, p2) and not np.allclose(p1, ref) else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} rattled configs have distinct, displaced positions"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} rattled configs are not distinct from each other/the reference"
    FAIL=$((FAIL+1))
fi


echo "Testing: every folder forces a single-point run through config_extra.fdf"
check_success mlff_config_001/config_extra.fdf
check_contains "MD.Steps              0" mlff_config_001/config_extra.fdf
check_contains "MD.VariableCell       false" mlff_config_001/config_extra.fdf
if [ "$(head -1 mlff_config_001/calc.fdf)" = "%include config_extra.fdf" ]; then
    echo -e "   -> ${GREEN}Verified:${NC} '%include config_extra.fdf' is the first line of calc.fdf"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} '%include config_extra.fdf' is not the first line of calc.fdf"
    FAIL=$((FAIL+1))
fi


# --- 2b. --volume-strains ---
echo -e "\n--- Testing --volume-strains (rattled + 2 strained configs) ---"
rm -rf mlff_config_*
stb-mlff --file calc.fdf --n-configs 2 --stdev 0.05 --volume-strains -3 3 --pseudo-dir . \
    --no-intro > log_strain.txt 2>&1
check_exit_code $? 0
check_contains "Generated 2 strained configuration(s)" log_strain.txt
check_success mlff_config_004/calc.fdf
echo "Testing: strained cells have volume (1 +/- 0.03) x reference, rattled configs keep it"
python3 -c "
import numpy as np, sys
from stb.core import structure_io
v = lambda p: abs(np.linalg.det(structure_io.read_fdf(p).lattice))
v0 = v('calc.fdf')
ok = (np.isclose(v('mlff_config_001/calc.fdf'), v0) and np.isclose(v('mlff_config_002/calc.fdf'), v0)
      and np.isclose(v('mlff_config_003/calc.fdf'), 0.97 * v0) and np.isclose(v('mlff_config_004/calc.fdf'), 1.03 * v0))
sys.exit(0 if ok else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} strained volumes are 0.97 / 1.03 x reference"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} strained volumes are wrong"
    FAIL=$((FAIL+1))
fi
stb-mlff --file calc.fdf --n-configs 1 --volume-strains -100 --pseudo-dir . --no-intro > log_strain_bad.txt 2>&1
check_exit_code $? 1
check_contains "must be > -100" log_strain_bad.txt
rm -rf mlff_config_*


# --- 2c. --n-test + the workflow file ---
echo -e "\n--- Testing --n-test (held-out test set) and mlff_workflow.json ---"
rm -rf mlff_config_* mlff_test mlff_workflow.json*
stb-mlff --file calc.fdf --n-configs 4 --stdev 0.05 0.1 --volume-strains -3 3 --n-test 5 \
    --pseudo-dir . --no-intro > log_ntest.txt 2>&1
check_exit_code $? 0
check_contains "Generated 5 held-out test configuration(s)" log_ntest.txt
check_contains "Configurations    : 6 training + 5 test" log_ntest.txt
check_success mlff_test/test_001/calc.fdf
check_success mlff_test/test_005/calc.fdf
check_success mlff_test/test_005/config_extra.fdf
check_success mlff_workflow.json
echo "Testing: workflow file lists train/test configs with kinds; test kinds cycle rattle/strain; tests differ from training"
python3 -c "
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
      and [c['kind'] for c in train] == ['rattle'] * 4 + ['strain'] * 2
      and kinds_test == ['rattle', 'rattle', 'strain', 'strain', 'rattle']
      and np.isclose(v('mlff_test/test_003'), 0.97 * v0) and np.isclose(v('mlff_test/test_004'), 1.03 * v0)
      and distinct and d['rounds'] == [] and d['history'][0]['stage'] == 'prep')
sys.exit(0 if ok else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} 6 train + 5 test entries, kinds/volumes right, test geometries independent"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} workflow file / test set not as expected"
    FAIL=$((FAIL+1))
fi
echo "Testing: a second run keeps the old workflow file as .bak"
stb-mlff --file calc.fdf --n-configs 1 --pseudo-dir . --no-intro > log_ntest2.txt 2>&1
check_contains "already existed -- copied to" log_ntest2.txt
check_success mlff_workflow.json.bak
echo "Testing: --n-test negative is rejected"
stb-mlff --file calc.fdf --n-configs 1 --n-test -1 --pseudo-dir . --no-intro > log_ntest_bad.txt 2>&1
check_exit_code $? 1
check_contains "--n-test must be >= 0" log_ntest_bad.txt
rm -rf mlff_config_* mlff_test mlff_workflow.json*


# --- 3. --n-configs 0 with no --from-aimd (expects failure) ---
echo -e "\n--- Testing --n-configs 0 with no --from-aimd (expects failure) ---"
stb-mlff --file calc.fdf --n-configs 0 --pseudo-dir . --no-intro > log_empty.txt 2>&1
check_exit_code $? 1
check_contains "No configurations were generated" log_empty.txt


# --- 4. Missing pseudopotential (expects failure) ---
echo -e "\n--- Testing a missing pseudopotential (expects failure) ---"
mkdir -p empty_pp_dir
stb-mlff --file calc.fdf --n-configs 2 --pseudo-dir empty_pp_dir --no-intro > log_missingpp.txt 2>&1
check_exit_code $? 1
check_contains "Missing pseudopotential" log_missingpp.txt


# --- 5. --save-report ---
echo -e "\n--- Testing --save-report ---"
rm -f stb_mlff_report.txt
stb-mlff --file calc.fdf --n-configs 2 --pseudo-dir . --save-report --no-intro > log_savereport.txt 2>&1
check_success stb_mlff_report.txt
check_contains "STB-MLFF REPORT" stb_mlff_report.txt
rm -f stb_mlff_report.txt


# --- 6. Robustness / errors ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: missing --file"
stb-mlff --no-intro > log_missing_file.txt 2>&1
check_exit_code $? 2

echo "Testing: nonexistent reference file"
stb-mlff --file nope.fdf --no-intro > log_missing_ref.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_ref.txt

echo "Testing: --version"
stb-mlff --version > log_version.txt 2>&1
check_contains "stb-mlff" log_version.txt

echo "Testing: --help documents --stdev, --from-aimd, --n-samples"
stb-mlff --help > log_help.txt 2>&1
check_contains "stdev" log_help.txt
check_contains "from-aimd" log_help.txt
check_contains "n-samples" log_help.txt
check_contains "volume-strains" log_help.txt
check_contains "n-test" log_help.txt
check_contains "workflow" log_help.txt


# --- 7. Interactive path (stb-suite, shortcut 4.17.1) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.17.1) ---"
rm -rf mlff_config_001 mlff_config_002 mlff_config_003
printf '4.17.1\ncalc.fdf\n\n3\n0.05\n\n2\nn\n\nn\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success mlff_config_001/calc.fdf
check_success mlff_config_003/calc.fdf
check_success mlff_test/test_002/calc.fdf
check_success mlff_workflow.json


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
