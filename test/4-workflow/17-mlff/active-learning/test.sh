#!/bin/bash

# --- Setup ---
# Smoke test for stb-mlffActiveLearning (Custom ML Force Field Active
# Learning, item 4.17.4)
#
# Reuses the 8 real finished-SIESTA fixtures from ../analysis/ to quickly
# fine-tune a tiny model on the fly (no large binary .model committed to
# the repo) via stb-mlffAnalysis, then feeds that model into
# stb-mlffActiveLearning to confirm the whole screening/ranking/writing
# pipeline works against a genuine model. Needs the optional 'ml' extra --
# skips gracefully (not a failure) if mace_run_train isn't on PATH.
export MPLBACKEND=Agg
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MLFF_FIXTURES="$FIXTURE_DIR/../analysis"
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


echo "--- Starting tester for stb-mlffActiveLearning (item 4.17.4) ---"

if ! command -v mace_run_train &> /dev/null; then
    echo -e "${YELLOW}SKIPPED${NC}: 'mace_run_train' not on PATH (needs 'pip install stb_suite[ml]')."
    exit 0
fi
if [ ! -d "$MLFF_FIXTURES" ]; then
    echo -e "${YELLOW}SKIPPED${NC}: ../analysis fixtures not found."
    exit 0
fi

rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp -r "$MLFF_FIXTURES"/mlff_config_* "$TEST_DIR/"
# The analysis fixtures deliberately have no pseudopotentials (stb-mlffAnalysis
# never needs them) -- add a placeholder one for stb-mlffActiveLearning's
# --pseudo-dir (default "."), content doesn't matter (only copied, never read).
echo "dummy pseudopotential content" > "$TEST_DIR/O.psf"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null

echo -e "\n--- Fine-tuning a quick model to screen with (--epochs 2) ---"
stb-mlffAnalysis --path "mlff_config_*" --epochs 2 --batch-size 4 --device cpu --committee 1 -o . \
    --name quick_model --skip-foundation-comparison --no-intro > log_quicktrain.txt 2>&1
QUICK_MODEL=$(find . -maxdepth 1 -name "quick_model*.model" | sort | tail -1)
if [ -z "$QUICK_MODEL" ]; then
    echo -e "${RED}FAIL${NC}: could not produce a quick fine-tuned model to test against."
    popd > /dev/null
    exit 1
fi
echo "Using model: $QUICK_MODEL"


# --- 1. Basic screening run ---
echo -e "\n--- Testing a basic screening run (--n-candidates 15 --top-k 3) ---"
stb-mlffActiveLearning --model "$QUICK_MODEL" --file mlff_config_001/calc.fdf --pseudo-dir . \
    --n-candidates 15 --stdev 0.2 --top-k 3 --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "Screened 15/15 candidate(s)" log_basic.txt
check_success mlff_config_009/calc.fdf
check_success mlff_config_010/calc.fdf
check_success mlff_config_011/calc.fdf
check_success mlff_config_009/O.psf

echo "Testing: numbering continues past the 8 already-present mlff_config_* folders"
check_contains "mlff_config_009" log_basic.txt

echo "Testing: single model -> score is the largest predicted force, with a committee hint"
check_contains "largest predicted force (single model" log_basic.txt
check_contains "train 2-3 more models" log_basic.txt
check_contains "Novelty reference: 8 config(s)" log_basic.txt

echo "Testing: suggested configs are ranked by decreasing predicted max|F|"
python3 -c "
import re, sys
text = open('log_basic.txt').read()
scores = [float(x) for x in re.findall(r'\(score ([\d.]+) eV/Ang', text)]
sys.exit(0 if scores == sorted(scores, reverse=True) and len(scores) == 3 else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} 3 suggested configs, ranked by decreasing predicted force"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} suggested configs not correctly ranked"
    FAIL=$((FAIL+1))
fi

echo "Testing: each suggested config keeps the reference basis/SCF settings"
check_contains "PAO.BasisSize       SZ" mlff_config_009/calc.fdf
check_contains "MD.VariableCell       false" mlff_config_009/config_extra.fdf


# --- 1b. --sampling-method md: MD-driven candidate generation (stb-mlmd's own dynamics) ---
echo -e "\n--- Testing --sampling-method md ---"
rm -rf mlff_config_009 mlff_config_010 mlff_config_011
stb-mlffActiveLearning --model "$QUICK_MODEL" --file mlff_config_001/calc.fdf --pseudo-dir . \
    --sampling-method md --md-steps 100 --md-temperature 800 --stride 20 --top-k 3 \
    --no-intro > log_md.txt 2>&1
check_exit_code $? 0
check_contains "Sampling method   : md" log_md.txt
check_contains "Screened 5/5 candidate(s)" log_md.txt
check_success mlff_config_009/calc.fdf
check_success mlff_config_010/calc.fdf
check_success mlff_config_011/calc.fdf

echo "Testing: MD-sampled configs also keep the reference basis/SCF settings"
check_contains "PAO.BasisSize       SZ" mlff_config_009/calc.fdf

echo "Testing: MD-sampled configs are also ranked by decreasing predicted max|F|"
python3 -c "
import re, sys
text = open('log_md.txt').read()
scores = [float(x) for x in re.findall(r'\(score ([\d.]+) eV/Ang', text)]
sys.exit(0 if scores == sorted(scores, reverse=True) and len(scores) == 3 else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} MD-sampled configs ranked by decreasing predicted force"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} MD-sampled configs not correctly ranked"
    FAIL=$((FAIL+1))
fi
rm -rf mlff_config_009 mlff_config_010 mlff_config_011


# --- 1c. Committee (2 models), novelty/force filters, several MD temperatures ---
echo -e "\n--- Fine-tuning a second quick model (different --seed) for a committee ---"
stb-mlffAnalysis --path "mlff_config_*" --epochs 2 --batch-size 4 --device cpu --seed 7 --committee 1 \
    --name quick_model_b --work-dir committee_b --skip-foundation-comparison --no-intro > log_quicktrain_b.txt 2>&1
QUICK_MODEL_B=$(find committee_b -maxdepth 1 -name "quick_model_b*.model" | sort | tail -1)
if [ -z "$QUICK_MODEL_B" ]; then
    echo -e "${RED}FAIL${NC}: could not produce the second committee model."
    FAIL=$((FAIL+1))
else
    echo -e "\n--- Testing a 2-model committee ---"
    stb-mlffActiveLearning --model "$QUICK_MODEL" "$QUICK_MODEL_B" --file mlff_config_001/calc.fdf \
        --pseudo-dir . --n-candidates 15 --stdev 0.2 --top-k 3 --no-intro > log_committee.txt 2>&1
    check_exit_code $? 0
    check_contains "committee force disagreement (2 models" log_committee.txt
    check_contains "Committee spread range" log_committee.txt
    python3 -c "
import re, sys
text = open('log_committee.txt').read()
scores = [float(x) for x in re.findall(r'\(score ([\d.]+) eV/Ang', text)]
spread = re.search(r'Committee spread range: ([\d.]+) - ([\d.]+)', text)
lo, hi = float(spread.group(1)), float(spread.group(2))
ok = len(scores) == 3 and scores == sorted(scores, reverse=True) and all(lo - 1e-3 <= s <= hi + 1e-3 for s in scores)
sys.exit(0 if ok else 1)
"
    if [ $? -eq 0 ]; then
        echo -e "   -> ${GREEN}Verified:${NC} committee picks ranked by decreasing spread, within the reported range"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} committee scores not ranked / outside the spread range"
        FAIL=$((FAIL+1))
    fi
    rm -rf mlff_config_009 mlff_config_010 mlff_config_011
fi

echo "Testing: --min-novelty huge -> every candidate skipped, clear warning, no folder"
stb-mlffActiveLearning --model "$QUICK_MODEL" --file mlff_config_001/calc.fdf --pseudo-dir . \
    --n-candidates 6 --top-k 2 --min-novelty 1000 --no-intro > log_novelty.txt 2>&1
check_exit_code $? 0
check_contains "Skipped (novelty < 1000)      : 6" log_novelty.txt
check_contains "Only 0 of the requested 2" log_novelty.txt
[ ! -d mlff_config_009 ] && { echo -e "   -> ${GREEN}Verified:${NC} no folder written"; PASS=$((PASS+1)); } \
    || { echo -e "   -> ${RED}Failed:${NC} a folder was written"; FAIL=$((FAIL+1)); rm -rf mlff_config_009; }

echo "Testing: --min-novelty 0 disables the filter"
stb-mlffActiveLearning --model "$QUICK_MODEL" --file mlff_config_001/calc.fdf --pseudo-dir . \
    --n-candidates 4 --top-k 1 --min-novelty 0 --output-dir novelty_off --no-intro > log_novelty_off.txt 2>&1
check_contains "Novelty filter    : off" log_novelty_off.txt
check_success novelty_off/mlff_config_001/calc.fdf
rm -rf novelty_off

echo "Testing: --max-force skips candidates above the limit"
stb-mlffActiveLearning --model "$QUICK_MODEL" --file mlff_config_001/calc.fdf --pseudo-dir . \
    --n-candidates 5 --top-k 1 --max-force 0.000001 --no-intro > log_maxforce.txt 2>&1
check_contains "Skipped (max|F| > 1e-06 eV/Ang) : 5" log_maxforce.txt

echo "Testing: several MD temperatures, one trajectory each"
stb-mlffActiveLearning --model "$QUICK_MODEL" --file mlff_config_001/calc.fdf --pseudo-dir . \
    --sampling-method md --md-steps 40 --md-temperature 300 900 --stride 20 --top-k 2 \
    --min-novelty 0 --output-dir md_temps --no-intro > log_md_temps.txt 2>&1
check_exit_code $? 0
check_contains "300, 900 K" log_md_temps.txt
check_contains "Screened 4/4 candidate(s)" log_md_temps.txt
rm -rf md_temps


# --- 1d. --supercell ---
echo -e "\n--- Testing --supercell 2 1 1 ---"
stb-mlffActiveLearning --model "$QUICK_MODEL" --file mlff_config_001/calc.fdf --pseudo-dir . \
    --n-candidates 4 --stdev 0.1 --top-k 1 --supercell 2 1 1 --min-novelty 0 \
    --output-dir supercell_test --no-intro > log_supercell.txt 2>&1
check_exit_code $? 0
check_contains "Supercell 2x1x1: 4 atoms (reference 2)" log_supercell.txt
check_contains "template has no k-grid (Gamma only)" log_supercell.txt
check_success supercell_test/mlff_config_001/calc.fdf
check_contains "NumberOfAtoms    4" supercell_test/mlff_config_001/calc.fdf
check_contains "PAO.BasisSize       SZ" supercell_test/mlff_config_001/calc.fdf
check_contains "MD.VariableCell       false" supercell_test/mlff_config_001/config_extra.fdf
echo "Testing: supercell .fdf = ase repeat() of the reference (cell, atom order), rattled ~0.1 Ang"
python3 -c "
import sys, numpy as np
from stb.core import structure_io as sio
from pymatgen.io.ase import AseAtomsAdaptor
get = lambda p: AseAtomsAdaptor.get_atoms(sio.to_pymatgen(sio.read_fdf(p)))
ref = get('mlff_config_001/calc.fdf').repeat((2, 1, 1))
new = get('supercell_test/mlff_config_001/calc.fdf')
d = np.linalg.norm(new.positions - ref.positions, axis=1)
ok = len(new) == 4 and np.allclose(new.cell[:], ref.cell[:]) and d.max() < 0.6 and new.get_chemical_symbols() == ref.get_chemical_symbols()
sys.exit(0 if ok else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} doubled cell, same atom order, displacements within rattle range"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} supercell geometry does not match the reference repeated"
    FAIL=$((FAIL+1))
fi
echo "Testing: replicate_fdf_cell with a k-grid template -> k-grid divided in config_extra.fdf"
python3 -c "
import sys
from stb.mlff_active_learning import supercell_kgrid_block
block, note = supercell_kgrid_block('kgrid.MonkhorstPack [4 4 4]\n', [2, 2, 1])
block2, note2 = supercell_kgrid_block('kgrid.MonkhorstPack [5 4 4]\n', [2, 1, 1])
ok = 'kgrid.MonkhorstPack [2 2 4]' in block and 'equivalent' in note and '[3 4 4]' in block2 and 'not divisible' in note2
sys.exit(0 if ok else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} 4x4x4 / 2x2x1 -> 2x2x4; 5x4x4 / 2x1x1 -> 3x4x4 with a warning"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} supercell k-grid not divided correctly"
    FAIL=$((FAIL+1))
fi
rm -rf supercell_test

echo "Testing: invalid --supercell"
stb-mlffActiveLearning --model "$QUICK_MODEL" --file mlff_config_001/calc.fdf --pseudo-dir . \
    --supercell 0 1 1 --no-intro > log_supercell_bad.txt 2>&1
check_exit_code $? 1
check_contains "dimensions must be >= 1" log_supercell_bad.txt


# --- 1e. Driven by mlff_workflow.json (committee + reference + training set from the file) ---
echo -e "\n--- Testing a run driven by mlff_workflow.json ---"
rm -rf mlff_config_009 mlff_config_01*
python3 -c "
import json
d = json.load(open('mlff_workflow.json'))
d['reference'] = {'file': 'mlff_config_001/calc.fdf', 'n_atoms': 2, 'species': ['O'], 'pseudo_dir': '.'}
json.dump(d, open('mlff_workflow.json', 'w'), indent=1)
"
LAST_ROUND=$(python3 -c "import json; print(max(r['round'] for r in json.load(open('mlff_workflow.json'))['rounds']))")
stb-mlffActiveLearning --n-candidates 6 --stdev 0.2 --top-k 2 --no-intro > log_workflow.txt 2>&1
check_exit_code $? 0
check_contains "Committee from    : mlff_workflow.json, round $LAST_ROUND" log_workflow.txt
check_contains "2 config(s) added to the training set" log_workflow.txt
check_success mlff_config_009/calc.fdf
python3 -c "
import json, sys
d = json.load(open('mlff_workflow.json'))
picked = [c for c in d['configs'] if c.get('stage') == 'active-learning' and c['folder'] in ('mlff_config_009', 'mlff_config_010')]
r = max(d['rounds'], key=lambda r: r['round'])
ok = (len(picked) == 2 and all(c['set'] == 'train' and c['kind'] == 'al-rattle' and c['round'] == r['round'] for c in picked)
      and r['active_learning'][-1]['picked'] == ['mlff_config_009', 'mlff_config_010']
      and d['history'][-1]['stage'] == 'active-learning')
sys.exit(0 if ok else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} picks recorded as training configs of the proposing round"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} picks not recorded in mlff_workflow.json"
    FAIL=$((FAIL+1))
fi
echo "Testing: the next run's novelty reference includes the recorded picks (8 + 2 configs)"
stb-mlffActiveLearning --n-candidates 4 --top-k 1 --output-dir wf_next --no-intro > log_workflow2.txt 2>&1
check_contains "Novelty reference: 10 config(s) from 'mlff_workflow.json (training configs)'" log_workflow2.txt
rm -rf wf_next mlff_config_009 mlff_config_010


# --- 2. --save-report ---
echo -e "\n--- Testing --save-report ---"
rm -f stb_mlffActiveLearning_report.txt
stb-mlffActiveLearning --model "$QUICK_MODEL" --file mlff_config_001/calc.fdf --pseudo-dir . \
    --n-candidates 5 --top-k 1 --output-dir report_test --save-report --no-intro > log_savereport.txt 2>&1
check_success stb_mlffActiveLearning_report.txt
check_contains "STB-MLFFACTIVELEARNING REPORT" stb_mlffActiveLearning_report.txt
rm -f stb_mlffActiveLearning_report.txt
rm -rf report_test


# --- 3. Robustness / errors ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: no --model and no workflow file -> clear error"
stb-mlffActiveLearning --workflow no_such_workflow.json --no-intro > log_missing_args.txt 2>&1
check_exit_code $? 1
check_contains "pass --model and --file" log_missing_args.txt

echo "Testing: workflow without a recorded reference .fdf -> asks for --file"
python3 -c "
import json
d = json.load(open('mlff_workflow.json')); d['reference'] = {}
json.dump(d, open('wf_noref.json', 'w'))
"
stb-mlffActiveLearning --workflow wf_noref.json --no-intro > log_missing_file.txt 2>&1
check_exit_code $? 1
check_contains "records no reference .fdf -- pass --file" log_missing_file.txt
rm -f wf_noref.json

echo "Testing: nonexistent model file"
stb-mlffActiveLearning --model does_not_exist.model --file mlff_config_001/calc.fdf --no-intro > log_missing_model.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_model.txt

echo "Testing: --version"
stb-mlffActiveLearning --version > log_version.txt 2>&1
check_contains "stb-mlffActiveLearning" log_version.txt

echo "Testing: --help documents the sampling, committee, supercell and filter options"
stb-mlffActiveLearning --help > log_help.txt 2>&1
check_contains "n-candidates" log_help.txt
check_contains "top-k" log_help.txt
check_contains "stdev" log_help.txt
check_contains "sampling-method" log_help.txt
check_contains "md-steps" log_help.txt
check_contains "md-temperature" log_help.txt
check_contains "supercell" log_help.txt
check_contains "min-novelty" log_help.txt
check_contains "max-force" log_help.txt
check_contains "train-path" log_help.txt
check_contains "committee" log_help.txt
check_contains "workflow" log_help.txt
check_contains "round" log_help.txt


# --- 4. Interactive path (stb-suite, shortcut 4.17.4) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.17.4) ---"
rm -rf mlff_config_012 mlff_config_013
printf '4.17.4\nn\n%s\nmlff_config_001/calc.fdf\n\n\n\n5\n0.2\n1\n.\n\nn\n\n0\n' "$QUICK_MODEL" | stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_contains "Suggested configs" log_menu.txt
rm -rf mlff_config_009

echo "Testing: menu with a 2-model committee and a supercell"
printf '4.17.4\nn\n%s %s\nmlff_config_001/calc.fdf\n\n\n2 1 1\n4\n0.1\n1\n.\n\nn\n\n0\n' "$QUICK_MODEL" "$QUICK_MODEL_B" | stb-suite > log_menu_committee.txt 2>&1
check_exit_code $? 0
check_contains "committee force disagreement (2 models" log_menu_committee.txt
check_contains "Supercell         : 2 x 1 x 1" log_menu_committee.txt
check_contains "Suggested configs" log_menu_committee.txt

echo "Testing: menu driven by mlff_workflow.json (Enter = use the latest round)"
printf '4.17.4\n\n\n\n4\n0.1\n1\n\nn\n\n0\n' | stb-suite > log_menu_workflow.txt 2>&1
check_exit_code $? 0
check_contains "Committee from    : mlff_workflow.json" log_menu_workflow.txt
check_contains "Suggested configs : 1" log_menu_workflow.txt


popd > /dev/null

# --- 5. Summary ---
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
