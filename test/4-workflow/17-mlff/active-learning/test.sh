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
#
# The quick training writes ./mlff_workflow.json (round 1), so picks go to
# ./configs_round_02/ (the round they will train). O2 rattled by 0.2 Ang has
# forces far above the default --max-force (15 eV/Ang): the runs that check the
# ranking pass --max-force 0 --max-energy 0 --frontier 0 (no force or energy limit) and --selection greedy (score order).
export MPLBACKEND=Agg
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MLFF_FIXTURES="$FIXTURE_DIR/../analysis"
TEST_DIR="$FIXTURE_DIR/test_files"
R2=configs_round_02
REF="--file mlff_config_001/calc.fdf --pseudo-dir ."

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

check_absent() {
    if [ -e "$1" ]; then
        echo -e " ... ${RED}FAIL${NC} ('$1' exists but should not)"
        FAIL=$((FAIL+1))
    else
        echo -e " ... ${GREEN}OK${NC} ('$1' not written)"
        PASS=$((PASS+1))
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

# Scores of the picks table ([4]), in pick order.
pick_scores='
import re, sys
def scores(path):
    t = open(path).read().split("[4] PICKED CONFIGURATIONS", 1)[1].split("[5]", 1)[0]
    return [float(m) for m in re.findall(r"^(?:\x1b\[\d+m)?\d+ \| [^|]+\| ([\d.]+) ", t, re.M)]
'


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
    --name quick_model --foundation-model small --no-stability --skip-foundation-comparison \
    --no-intro > log_quicktrain.txt 2>&1
QUICK_MODEL=$(find . -maxdepth 1 -name "quick_model_compiled.model" | sort | tail -1)
if [ -z "$QUICK_MODEL" ]; then
    echo -e "${RED}FAIL${NC}: could not produce a quick fine-tuned model to test against."
    popd > /dev/null
    exit 1
fi
echo "Using model: $QUICK_MODEL"


# --- 1. Basic screening run (workflow file present: picks -> configs_round_02/) ---
echo -e "\n--- Testing a basic screening run (--n-candidates 15 --top-k 3) ---"
stb-mlffActiveLearning --model "$QUICK_MODEL" $REF --n-candidates 15 --stdev 0.2 --top-k 3 \
    --max-force 0 --max-energy 0 --frontier 0 --selection greedy --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
for sec in "\[0\] RUN METADATA" "\[1\] REFERENCE, COMMITTEE & TRAINING SET" "\[2\] CANDIDATE POOL" "\[3\] SELECTION" \
           "\[4\] PICKED CONFIGURATIONS" "\[5\] IS ANOTHER ROUND WORTH IT?" "\[6\] ACTIVE-LEARNING HISTORY" \
           "\[7\] SIESTA COST ESTIMATE" "\[8\] SUMMARY & FILES" "\[9\] NEXT STEPS" "\[10\] WARNINGS"; do
    check_contains "$sec" log_basic.txt
done
check_contains "Screened          : 15 of 15 candidate(s)" log_basic.txt
check_contains "Output            : ./$R2/ (the configs that will train round 2)" log_basic.txt
check_contains "seed 2042 (--seed 42 + 1000 x round 2)" log_basic.txt
check_contains "Committee from    : --model (1 model(s)); picks recorded as proposed by round 1" log_basic.txt
check_success $R2/mlff_config_009/calc.fdf
check_success $R2/mlff_config_010/calc.fdf
check_success $R2/mlff_config_011/calc.fdf
check_success $R2/mlff_config_009/O.psf
check_success mlff_al_r01_1.png

echo "Testing: numbering continues past the 8 already-present mlff_config_* folders (older flat layout)"
check_contains "Pick 1            : ./$R2/mlff_config_009" log_basic.txt

echo "Testing: single model -> score is the largest predicted force, with a committee hint"
check_contains "Score             : largest predicted force (single model)" log_basic.txt
check_contains "Score: single model -- the score is its largest predicted force" log_basic.txt
check_contains "Training set      : 8 config(s) from 'mlff_workflow.json (training configs)', 8 labelled by SIESTA" log_basic.txt
check_contains "Training range    : SIESTA max|F|" log_basic.txt
check_contains "New territory     :" log_basic.txt
check_contains "Estimate          :" log_basic.txt

echo "Testing: where to run SIESTA, and the warning count in [8] matches the list in [10]"
check_contains "Run SIESTA in     : $R2/ (3 new folder(s)" log_basic.txt
check_contains "SIESTA: run it in $R2/ (3 new folder(s))" log_basic.txt
check_python "status line counts every warning listed in [10]" "
import re, sys
t = open('log_basic.txt').read()
n = int(re.search(r'\((\d+) warning\(s\), see \[10\]\)', t).group(1))
listed = re.findall(r' (\d+)\. ', t.split('[10] WARNINGS', 1)[1])
sys.exit(0 if n == len(listed) and n > 0 else 1)
"

echo "Testing: suggested configs are ranked by decreasing predicted max|F| (greedy selection)"
check_python "3 suggested configs, ranked by decreasing predicted force" "$pick_scores
s = scores('log_basic.txt')
sys.exit(0 if len(s) == 3 and s == sorted(s, reverse=True) else 1)
"

echo "Testing: each suggested config keeps the reference basis/SCF settings and is wrapped into the cell"
check_contains "PAO.BasisSize       SZ" $R2/mlff_config_009/calc.fdf
check_contains "MD.VariableCell       false" $R2/mlff_config_009/config_extra.fdf
check_python "picks recorded in the workflow file (proposed by round 1, with d_min), the run with its seed and pool" "
import json, sys
d = json.load(open('mlff_workflow.json'))
picked = [c for c in d['configs'] if c.get('stage') == 'active-learning']
al = d['rounds'][0]['active_learning'][-1]
ok = (len(picked) == 3 and all(c['folder'].startswith('$R2/') and c['round'] == 1 and 'd_min' in c for c in picked)
      and al['seed'] == 2042 and al['seed_base'] == 42 and al['pool']['n'] == 15 and al['selection'] == 'greedy'
      and al['picked'] == ['$R2/mlff_config_009', '$R2/mlff_config_010', '$R2/mlff_config_011'])
sys.exit(0 if ok else 1)
"


# --- 1b. The default --max-force (the workflow's 'unphysical' line) and --dry-run ---
echo -e "\n--- Testing the default --max-force and --dry-run ---"
cp mlff_workflow.json wf_before.json
stb-mlffActiveLearning --model "$QUICK_MODEL" $REF --n-candidates 15 --stdev 0.2 --top-k 3 --dry-run \
    --no-intro > log_dry.txt 2>&1
check_exit_code $? 0
check_contains "max|F| <= 15 eV/Ang" log_dry.txt
check_contains "DRY RUN, nothing is written" log_dry.txt
check_python "candidates above 15 eV/Ang are skipped by default" "
import re, sys
m = re.search(r'Skipped: max\|F\| > --max-force +\| (\d+)', open('log_dry.txt').read())
sys.exit(0 if m and int(m.group(1)) > 0 else 1)
"
check_absent $R2/mlff_config_012
cmp -s mlff_workflow.json wf_before.json && { echo -e "   -> ${GREEN}Verified:${NC} workflow file untouched"; PASS=$((PASS+1)); } \
    || { echo -e "   -> ${RED}Failed:${NC} the dry run changed the workflow file"; FAIL=$((FAIL+1)); }
check_python "the default --max-force is the workflow's shared threshold" "
import sys, subprocess
from stb.core import mlff_workflow as w
h = subprocess.run(['stb-mlffActiveLearning', '--help'], capture_output=True, text=True).stdout
sys.exit(0 if w.UNPHYSICAL_FMAX == 15.0 and 'default: 15' in h else 1)
"


# --- 1c. --sampling-method md ---
echo -e "\n--- Testing --sampling-method md ---"
stb-mlffActiveLearning --model "$QUICK_MODEL" $REF --sampling-method md --md-steps 100 --md-temperature 800 \
    --stride 20 --top-k 3 --max-force 0 --max-energy 0 --frontier 0 --min-novelty 0 --selection greedy --dry-run --no-intro > log_md.txt 2>&1
check_exit_code $? 0
check_contains "Sampling          : md, 100 NVT steps at 800 K" log_md.txt
check_contains "Screened          : 5 of 5 candidate(s)" log_md.txt
check_contains "T target (K)" log_md.txt
check_contains "md 800 K" log_md.txt
check_python "MD-sampled configs ranked by decreasing predicted force" "$pick_scores
s = scores('log_md.txt')
sys.exit(0 if len(s) == 3 and s == sorted(s, reverse=True) else 1)
"

echo "Testing: several MD temperatures, one trajectory each"
stb-mlffActiveLearning --model "$QUICK_MODEL" $REF --sampling-method md --md-steps 40 --md-temperature 300 900 \
    --stride 20 --top-k 2 --min-novelty 0 --max-force 0 --max-energy 0 --frontier 0 --dry-run --no-intro > log_md_temps.txt 2>&1
check_exit_code $? 0
check_contains "at 300, 900 K" log_md_temps.txt
check_contains "Screened          : 4 of 4 candidate(s)" log_md_temps.txt

echo "Testing: --md-ensemble npt needs a bulk cell (this molecule box has vacuum)"
stb-mlffActiveLearning --model "$QUICK_MODEL" $REF --sampling-method md --md-ensemble npt --dry-run \
    --no-intro > log_npt.txt 2>&1
check_exit_code $? 1
check_contains "md-ensemble npt needs a bulk cell" log_npt.txt


# --- 1d. Committee (2 models): energy score, diverse vs greedy selection, target verdicts ---
echo -e "\n--- Fine-tuning a second quick model (different --seed) for a committee ---"
stb-mlffAnalysis --path "mlff_config_*" --epochs 2 --batch-size 4 --device cpu --seed 7 --committee 1 \
    --foundation-model small --no-stability --name quick_model_b --work-dir committee_b \
    --workflow wf_b.json --skip-foundation-comparison --no-intro > log_quicktrain_b.txt 2>&1
QUICK_MODEL_B=$(find committee_b -maxdepth 1 -name "quick_model_b_compiled.model" | sort | tail -1)
if [ -z "$QUICK_MODEL_B" ]; then
    echo -e "${RED}FAIL${NC}: could not produce the second committee model."
    FAIL=$((FAIL+1))
else
    echo -e "\n--- Testing a 2-model committee ---"
    stb-mlffActiveLearning --model "$QUICK_MODEL" "$QUICK_MODEL_B" $REF --n-candidates 15 --stdev 0.2 --top-k 3 \
        --max-force 0 --max-energy 0 --frontier 0 --selection greedy --dry-run --save-data --no-intro > log_committee.txt 2>&1
    check_exit_code $? 0
    check_contains "Score             : force disagreement" log_committee.txt
    check_contains "F spread med/max" log_committee.txt
    check_contains "Energy spread     :" log_committee.txt
    check_contains "Spread vs force   : corr" log_committee.txt
    check_contains "Committee quality : --model given" log_committee.txt
    check_python "committee picks ranked by decreasing spread, none above the pool maximum" "$pick_scores
s = scores('log_committee.txt')
hi = float(re.search(r'max ([\d.]+) eV/Ang', open('log_committee.txt').read().split('Score (force)', 1)[1]).group(1))
sys.exit(0 if len(s) == 3 and s == sorted(s, reverse=True) and max(s) <= hi + 1e-3 else 1)
"
    check_success mlff_al_r01_2_pool.dat
    check_python "the pool data file has one line per candidate" "
import sys
rows = [l for l in open('mlff_al_r01_2_pool.dat') if not l.startswith('#')]
sys.exit(0 if len(rows) == 15 else 1)
"

    echo "Testing: --score energy (spread of the energy per atom, members' constant offsets removed)"
    stb-mlffActiveLearning --model "$QUICK_MODEL" "$QUICK_MODEL_B" $REF --n-candidates 10 --stdev 0.1 --top-k 2 \
        --max-force 0 --max-energy 0 --frontier 0 --score energy --dry-run --no-intro > log_energy.txt 2>&1
    check_exit_code $? 0
    check_contains "Score (energy)    :" log_energy.txt
    check_contains "eV/atom" log_energy.txt
    check_python "the energy spread is offset-free: a constant shift of one member changes nothing" "
import sys, numpy as np
from stb.mlff_active_learning import center_energy_spread
rng = np.random.default_rng(0)
e = rng.normal(size=(6, 2)) * 0.01
evs = [{'e_models': list(r), 'sigma_e': 0.0} for r in e]
center_energy_spread(evs); a = [v['sigma_e'] for v in evs]
evs2 = [{'e_models': [r[0] + 5.0, r[1]], 'sigma_e': 0.0} for r in e]
center_energy_spread(evs2); b = [v['sigma_e'] for v in evs2]
sys.exit(0 if np.allclose(a, b) and max(a) < 0.05 else 1)
"

    echo "Testing: diverse selection (default) is compared with the greedy one"
    stb-mlffActiveLearning --model "$QUICK_MODEL" "$QUICK_MODEL_B" $REF --sampling-method md --md-steps 200 \
        --md-temperature 600 --stride 5 --top-k 3 --min-novelty 0 --max-force 0 --max-energy 0 --frontier 0 --dry-run --no-intro > log_diverse.txt 2>&1
    check_exit_code $? 0
    check_contains "Selection         : diverse (diversity 1)" log_diverse.txt
    check_contains "vs greedy         :" log_diverse.txt

    echo "Testing: --target-force verdicts (no Stage-3 calibration for --model: the spread is a lower bound)"
    stb-mlffActiveLearning --model "$QUICK_MODEL" "$QUICK_MODEL_B" $REF --n-candidates 8 --stdev 0.1 --top-k 1 \
        --max-force 0 --max-energy 0 --frontier 0 --target-force 1e9 --dry-run --no-intro > log_target_hi.txt 2>&1
    check_contains "\[UNSURE\]" log_target_hi.txt
    check_contains "LOWER bound" log_target_hi.txt
    stb-mlffActiveLearning --model "$QUICK_MODEL" "$QUICK_MODEL_B" $REF --n-candidates 8 --stdev 0.1 --top-k 1 \
        --max-force 0 --max-energy 0 --frontier 0 --target-force 1e-6 --dry-run --no-intro > log_target_lo.txt 2>&1
    check_contains "\[CONTINUE\]" log_target_lo.txt

    echo "Testing: --view-picks writes the picks with the per-atom spread (no display here)"
    env -u DISPLAY -u WAYLAND_DISPLAY stb-mlffActiveLearning --model "$QUICK_MODEL" "$QUICK_MODEL_B" $REF \
        --n-candidates 8 --stdev 0.1 --top-k 2 --max-force 0 --max-energy 0 --frontier 0 --dry-run --view-picks --no-intro > log_view.txt 2>&1
    check_python "picks file carries the per-atom committee spread" "
import glob, sys
from ase.io import read
f = sorted(glob.glob('mlff_al_r01_*_picks.xyz'))
a = read(f[-1], ':') if f else []
sys.exit(0 if len(a) == 2 and 'spread' in a[0].arrays else 1)
"
fi

echo "Testing: --min-novelty huge -> every candidate skipped, clear warning, no folder"
stb-mlffActiveLearning --model "$QUICK_MODEL" $REF --workflow none.json -o nov_huge --train-path "mlff_config_*" \
    --n-candidates 6 --top-k 2 --min-novelty 1000 --max-force 0 --max-energy 0 --frontier 0 --selection greedy --no-intro > log_novelty.txt 2>&1
check_exit_code $? 0
check_contains "Skipped: novelty < --min-novelty | 6" log_novelty.txt
check_contains "Selection: only 0 of the requested 2" log_novelty.txt
check_absent nov_huge/mlff_config_001

echo "Testing: --min-novelty 0 disables the filter (no workflow file: picks in the output dir itself)"
stb-mlffActiveLearning --model "$QUICK_MODEL" $REF --workflow none.json -o novelty_off --n-candidates 4 --top-k 1 \
    --min-novelty 0 --max-force 0 --max-energy 0 --frontier 0 --no-intro > log_novelty_off.txt 2>&1
check_contains "novelty off" log_novelty_off.txt
check_success novelty_off/mlff_config_001/calc.fdf
rm -rf novelty_off nov_huge

echo "Testing: --max-force skips candidates above the limit"
stb-mlffActiveLearning --model "$QUICK_MODEL" $REF --n-candidates 5 --top-k 1 --max-force 0.000001 --dry-run \
    --no-intro > log_maxforce.txt 2>&1
check_contains "Skipped: max|F| > --max-force    | 5" log_maxforce.txt

echo "Testing: --auto-sampling scans the stdev and picks the largest still-physical one"
stb-mlffActiveLearning --model "$QUICK_MODEL" $REF --auto-sampling --scan-values 0.02 0.05 0.1 --scan-n 4 \
    --n-candidates 6 --top-k 1 --max-force 0 --max-energy 0 --frontier 0 --dry-run --no-intro > log_auto.txt 2>&1
check_exit_code $? 0
check_contains "Auto-sampling scan (4 candidates per setting" log_auto.txt
check_contains "Chosen            :" log_auto.txt
check_contains "E med (meV/atom)" log_auto.txt

echo "Testing: the energy cap (--max-energy, --frontier x the training set's 90th percentile)"
stb-mlffActiveLearning --model "$QUICK_MODEL" $REF --n-candidates 6 --stdev 0.2 --top-k 2 --max-force 0 \
    --max-energy 1e-6 --frontier 0 --dry-run --no-intro > log_ecap.txt 2>&1
check_exit_code $? 0
check_contains "Energy cap        : picks at most 0 meV/atom above the reference" log_ecap.txt
check_contains "above the cap (0)" log_ecap.txt
check_contains "Selection: the energy cap" log_ecap.txt
check_python "choose_scanned: the frontier stops at 1.5 x the training energy; 3/2 kT conversion" "
import sys
from stb.mlff_active_learning import choose_scanned, equipartition_t
rows = [dict(value=v, n=4, stopped=False, close_pct=0, maxf_med=1.0, de_med=de, score_med=sc)
        for v, de, sc in [(0.05, 60, 0.03), (0.1, 250, 0.05), (0.15, 560, 0.12), (0.2, 980, 0.24)]]
a, _ = choose_scanned(rows, 15, max_energy=1000, train_de=370, frontier=1.5)   # limit 555 -> 0.1
b, _ = choose_scanned(rows, 15, max_energy=1000, train_de=None)                # old rule -> 0.2
c, why = choose_scanned(rows, 15, max_energy=1000, train_de=20, frontier=1.5)   # all beyond -> mildest
ok = a == 0.1 and b == 0.2 and c == 0.05 and 'mildest' in why and abs(equipartition_t(1.5 * 0.08617333 * 300) - 300) < 1e-6
sys.exit(0 if ok else 1)
"

echo "Testing: rattle candidates do not repeat Stage 1's displacement directions, and each round draws new ones"
check_python "cos(Stage-1 config i, candidate i) far from 1; rounds 2 and 3 give different pools" "
import sys, numpy as np
from ase.build import bulk
from stb.mlff import build_rattled_configs
from stb.mlff_active_learning import rattle_candidates, ROUND_SEED_STEP
a = bulk('Si', 'diamond', a=5.43, cubic=True); p0 = a.get_positions()
prep = build_rattled_configs(a, 'cartesian', 6, [0.05], 42)
cand = [c[0].get_positions() for c in rattle_candidates(a, 6, [0.15], 42 + ROUND_SEED_STEP * 2)]
cand3 = [c[0].get_positions() for c in rattle_candidates(a, 6, [0.15], 42 + ROUND_SEED_STEP * 3)]
cos = lambda x, y: abs(np.dot((x - p0).ravel(), (y - p0).ravel())) / np.linalg.norm(x - p0) / np.linalg.norm(y - p0)
sys.exit(0 if max(cos(prep[i], cand[i]) for i in range(6)) < 0.9 and not np.allclose(cand[0], cand3[0]) else 1)
"


# --- 1e. --supercell ---
echo -e "\n--- Testing --supercell 2 1 1 ---"
stb-mlffActiveLearning --model "$QUICK_MODEL" $REF --workflow none.json --n-candidates 4 --stdev 0.1 --top-k 1 \
    --supercell 2 1 1 --min-novelty 0 --max-force 0 --max-energy 0 --frontier 0 --output-dir supercell_test --no-intro > log_supercell.txt 2>&1
check_exit_code $? 0
check_contains "Supercell         : 2x1x1 = 4 atoms; template has no k-grid (Gamma only)" log_supercell.txt
check_contains "Supercell: template has no k-grid" log_supercell.txt
check_success supercell_test/mlff_config_001/calc.fdf
check_contains "NumberOfAtoms    4" supercell_test/mlff_config_001/calc.fdf
check_contains "PAO.BasisSize       SZ" supercell_test/mlff_config_001/calc.fdf
check_contains "MD.VariableCell       false" supercell_test/mlff_config_001/config_extra.fdf
check_python "doubled cell, same atom order, displacements within the rattle range" "
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
check_python "4x4x4 / 2x2x1 -> 2x2x4; 5x4x4 / 2x1x1 -> 3x4x4 with a warning" "
import sys
from stb.mlff_active_learning import supercell_kgrid_block
block, note = supercell_kgrid_block('kgrid.MonkhorstPack [4 4 4]\n', [2, 2, 1])
block2, note2 = supercell_kgrid_block('kgrid.MonkhorstPack [5 4 4]\n', [2, 1, 1])
ok = 'kgrid.MonkhorstPack [2 2 4]' in block and 'equivalent' in note and '[3 4 4]' in block2 and 'not divisible' in note2
sys.exit(0 if ok else 1)
"
rm -rf supercell_test

echo "Testing: invalid --supercell"
stb-mlffActiveLearning --model "$QUICK_MODEL" $REF --supercell 0 1 1 --no-intro > log_supercell_bad.txt 2>&1
check_exit_code $? 1
check_contains "dimensions must be >= 1" log_supercell_bad.txt


# --- 1f. --sampling-method strain on a 2D cell (the vacuum axis is not strained) ---
echo -e "\n--- Testing --sampling-method strain on a 2D slab ---"
mkdir -p slab && cp O.psf slab/
cat > slab/slab.fdf <<'FDF'
SystemLabel slab
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
PAO.BasisSize SZ
FDF
stb-mlffActiveLearning --model "$QUICK_MODEL" --file slab/slab.fdf --pseudo-dir slab --workflow none.json \
    -o strain_test --sampling-method strain --strain-max 4 --stdev 0.01 --n-candidates 6 --top-k 2 \
    --min-novelty 0 --max-force 0 --max-energy 0 --frontier 0 --no-intro > log_strain.txt 2>&1
check_exit_code $? 0
check_contains "Sampling          : strain, 6 candidates, strain up to 4 %" log_strain.txt
check_contains "with their own cell" log_strain.txt
check_python "each pick has its own in-plane cell, the vacuum axis stays 15 Ang" "
import sys, glob, numpy as np
from stb.core import structure_io as s
c0 = s.read_fdf('slab/slab.fdf').lattice
cells = [s.read_fdf(f).lattice for f in sorted(glob.glob('strain_test/mlff_config_*/slab.fdf'))]
ok = len(cells) == 2 and all(np.isclose(c[2][2], 15.0) and not np.allclose(c[:2], c0[:2]) for c in cells)
sys.exit(0 if ok else 1)
"
rm -rf strain_test


# --- 1g. structure.fdf + calc.fdf pair (-s / -c) ---
echo -e "\n--- Testing -s structure.fdf -c calc.fdf (explicit, then from the workflow) ---"
mkdir -p pair
# NumberOfAtoms stays in calc.fdf on purpose: a supercell must override it.
python3 - <<'PY'
lines = open('mlff_config_001/calc.fdf').read().splitlines(keepends=True)
start = next(i for i, l in enumerate(lines) if l.startswith('NumberOfSpecies'))
end = next(i for i, l in enumerate(lines) if l.startswith('%endblock AtomicCoordinates')) + 1
keep = lambda l: l.startswith(('PAO.', 'MeshCutoff'))
open('pair/structure.fdf', 'w').writelines(l for l in lines[start:end] if not keep(l))
open('pair/calc.fdf', 'w').writelines(lines[:start] + ['%include structure.fdf\n'] +
                                      [l for l in lines[start:end] if keep(l)] + lines[end:])
PY
stb-mlffActiveLearning --model "$QUICK_MODEL" -s pair/structure.fdf -c pair/calc.fdf --pseudo-dir . --workflow none.json \
    --n-candidates 4 --stdev 0.1 --top-k 1 --supercell 2 1 1 --min-novelty 0 --max-force 0 --max-energy 0 --frontier 0 \
    --output-dir pair_test --no-intro > log_pair.txt 2>&1
check_exit_code $? 0
check_contains "Reference         : pair/structure.fdf + pair/calc.fdf" log_pair.txt
check_success pair_test/mlff_config_001/structure.fdf
check_success pair_test/mlff_config_001/calc.fdf
check_contains "NumberOfAtoms 4" pair_test/mlff_config_001/config_extra.fdf
check_contains "%include structure.fdf" pair_test/mlff_config_001/calc.fdf
check_python "supercell written to structure.fdf (4 atoms, a doubled)" "
import sys
from stb.core import structure_io as sio
s = sio.read_fdf('pair_test/mlff_config_001/structure.fdf')
sys.exit(0 if len(s.atoms) == 4 and abs(s.lattice[0][0] - 20.0) < 1e-6 else 1)
"
echo "Testing: -s/-c together with --file -> clear error"
stb-mlffActiveLearning --model "$QUICK_MODEL" -s pair/structure.fdf --file mlff_config_001/calc.fdf \
    --no-intro > log_pair_both.txt 2>&1
check_exit_code $? 1
check_contains "not both" log_pair_both.txt
rm -rf pair_test


# --- 1h. Driven by mlff_workflow.json (committee + reference + training set from the file) ---
echo -e "\n--- Testing a run driven by mlff_workflow.json ---"
python3 -c "
import json
d = json.load(open('mlff_workflow.json'))
d['reference'] = {'file': 'mlff_config_001/calc.fdf', 'n_atoms': 2, 'species': ['O'], 'pseudo_dir': '.'}
json.dump(d, open('mlff_workflow.json', 'w'), indent=1)
"
stb-mlffActiveLearning --n-candidates 6 --stdev 0.2 --top-k 2 --max-force 0 --max-energy 0 --frontier 0 --no-intro > log_workflow.txt 2>&1
check_exit_code $? 0
check_contains "Committee from    : round 1 (1 model(s): seed 123)" log_workflow.txt
check_contains "2 config(s) added to the training set, proposed by round 1" log_workflow.txt
check_contains "Pick 1            : ./$R2/mlff_config_012" log_workflow.txt
check_contains "Round.run" log_workflow.txt
check_python "picks recorded as training configs of the proposing round, second AL entry of round 1" "
import json, sys
d = json.load(open('mlff_workflow.json'))
r = d['rounds'][0]
ok = (len(r['active_learning']) == 2 and r['active_learning'][-1]['picked'] == ['$R2/mlff_config_012', '$R2/mlff_config_013']
      and d['history'][-1]['stage'] == 'active-learning')
sys.exit(0 if ok else 1)
"
echo "Testing: the next run's training set includes the recorded picks (8 + 3 + 2 configs)"
stb-mlffActiveLearning --n-candidates 4 --top-k 1 --dry-run --no-intro > log_workflow2.txt 2>&1
check_contains "Training set      : 13 config(s) from 'mlff_workflow.json (training configs)'" log_workflow2.txt
echo "Testing: workflow recording a structure.fdf + calc.fdf pair"
cp mlff_workflow.json wf_single.json
python3 -c "
import json
d = json.load(open('mlff_workflow.json'))
d['reference'].update({'file': 'pair/structure.fdf', 'calc': 'pair/calc.fdf'})
json.dump(d, open('mlff_workflow.json', 'w'), indent=1)
"
stb-mlffActiveLearning --n-candidates 4 --stdev 0.2 --top-k 1 --max-force 0 --max-energy 0 --frontier 0 --min-novelty 0 --output-dir wf_pair \
    --no-intro > log_workflow_pair.txt 2>&1
check_exit_code $? 0
check_contains "Reference         : pair/structure.fdf + pair/calc.fdf" log_workflow_pair.txt
check_success wf_pair/$R2/mlff_config_001/structure.fdf
check_success wf_pair/$R2/mlff_config_001/calc.fdf
mv wf_single.json mlff_workflow.json
rm -rf wf_pair


# --- 2. --save-report (written in the output directory) ---
echo -e "\n--- Testing --save-report ---"
stb-mlffActiveLearning --model "$QUICK_MODEL" $REF --workflow none.json --n-candidates 5 --top-k 1 --max-force 0 --max-energy 0 --frontier 0 \
    --output-dir report_test --save-report --no-intro > log_savereport.txt 2>&1
check_success report_test/stb_mlffActiveLearning_report.txt
check_contains "STB-MLFFACTIVELEARNING REPORT" report_test/stb_mlffActiveLearning_report.txt
check_contains "\[10\] WARNINGS" report_test/stb_mlffActiveLearning_report.txt
rm -rf report_test


# --- 2b. Dispersion: picks keep Stage 1's DFT-D3 ---
echo "Testing: --d3 (and a workflow file recording Stage 1 --d3) puts DFTD3 into every pick's config_extra.fdf"
stb-mlffActiveLearning --model "$QUICK_MODEL" $REF --workflow none.json --n-candidates 4 --top-k 1 --max-force 0 \
    --max-energy 0 --frontier 0 --d3 --output-dir d3_picks --no-intro > log_d3.txt 2>&1
check_exit_code $? 0
check_contains "Dispersion        : DFT-D3 in every pick's config_extra.fdf (--d3)" log_d3.txt
check_contains "DFTD3                 .true." d3_picks/mlff_config_001/config_extra.fdf
check_python "dftd3_block follows the workflow's record; fdf_true reads fdf logicals" "
import sys
from stb.core import mlff_workflow as w
on = {'reference': {'dftd3': {'on': True, 'forced': True, 'source': '--d3 (config_extra.fdf)'}}}
tmpl = {'reference': {'dftd3': {'on': True, 'forced': False, 'source': 'template'}}}
ok = ('DFTD3' in w.dftd3_block(on) and w.dftd3_block(tmpl) == '' and w.dftd3_block({}) == ''
      and w.fdf_true('.true.') and w.fdf_true('T') and not w.fdf_true('.false.') and not w.fdf_true('default'))
sys.exit(0 if ok else 1)
"
rm -rf d3_picks


# --- 3. Robustness / errors ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: no --model and no workflow file -> clear error"
stb-mlffActiveLearning --workflow no_such_workflow.json --no-intro > log_missing_args.txt 2>&1
check_exit_code $? 1
check_contains "pass --model and -s/-c (or --file)" log_missing_args.txt

echo "Testing: workflow without a recorded reference .fdf -> asks for -s/-c or --file"
python3 -c "
import json
d = json.load(open('mlff_workflow.json')); d['reference'] = {}
json.dump(d, open('wf_noref.json', 'w'))
"
stb-mlffActiveLearning --workflow wf_noref.json --no-intro > log_missing_file.txt 2>&1
check_exit_code $? 1
check_contains "records no reference .fdf -- pass -s/-c or --file" log_missing_file.txt
rm -f wf_noref.json

echo "Testing: --round K that does not exist lists the rounds"
stb-mlffActiveLearning --round 9 --no-intro > log_round9.txt 2>&1
check_exit_code $? 1
check_contains "(it has: 1)" log_round9.txt

echo "Testing: nonexistent model file"
stb-mlffActiveLearning --model does_not_exist.model $REF --no-intro > log_missing_model.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_model.txt

echo "Testing: --version"
stb-mlffActiveLearning --version > log_version.txt 2>&1
check_contains "stb-mlffActiveLearning" log_version.txt

echo "Testing: --help documents the options"
stb-mlffActiveLearning --help > log_help.txt 2>&1
for opt in n-candidates top-k stdev sampling-method md-steps md-temperature supercell min-novelty max-force train-path \
           committee workflow round strain-max md-ensemble "\-\-score" selection diversity auto-sampling scan-values \
           target-force dry-run "\-\-show" save-data view-picks configs_round_NN; do
    check_contains "$opt" log_help.txt
done


# --- 4. Interactive path (stb-suite, shortcut 4.17.4) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.17.4) ---"
# mlff_workflow.json is here: 'n' = give the models; then model, structure, calc '-' (one file), pseudos,
# sampling 1, no supercell, no auto-sampling, 5 candidates, stdev 0.02, force score, diverse, no target,
# top-k 1, output dir, cpu, no dry run / show / data / viewer / report; Enter; 0 = quit.
printf '4.17.4\nn\n%s\nmlff_config_001/calc.fdf\n-\n\n\n\n\n5\n0.02\n\n\n\n1\n.\n\nn\nn\nn\nn\nn\n\n0\n' "$QUICK_MODEL" \
    | stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_contains "Pick 1            :" log_menu.txt
check_contains "Let a quick scan choose the sampling" log_menu.txt

if [ -n "$QUICK_MODEL_B" ]; then
echo "Testing: menu with a 2-model committee and a supercell (dry run)"
printf '4.17.4\nn\n%s %s\nmlff_config_001/calc.fdf\n-\n\n\n2 1 1\n\n4\n0.02\n\n\n\n1\n.\n\ny\nn\nn\nn\nn\n\n0\n' "$QUICK_MODEL" "$QUICK_MODEL_B" \
    | stb-suite > log_menu_committee.txt 2>&1
check_exit_code $? 0
check_contains "Score             : force disagreement" log_menu_committee.txt
check_contains "Supercell         : 2x1x1 = 4 atoms" log_menu_committee.txt
check_contains "DRY RUN" log_menu_committee.txt
fi

echo "Testing: menu driven by mlff_workflow.json (Enter = use it; one round, so no round question)"
printf '4.17.4\n\n\n\n\n4\n0.02\n\n\n\n1\n\nn\nn\nn\nn\nn\n\n0\n' | stb-suite > log_menu_workflow.txt 2>&1
check_exit_code $? 0
check_contains "Committee from    : round 1" log_menu_workflow.txt
check_contains "Status            : OK" log_menu_workflow.txt


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
