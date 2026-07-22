#!/bin/bash

# --- Setup ---
# Smoke test for stb-mlsearch (ML Structure Search, item 5.4)
#
# Needs the optional 'ml' extra (pip install stb_suite[ml] -- PyTorch +
# mace-torch). The whole file is skipped with a clear message if `mace`
# isn't importable, same gating pattern as the other test/5-mlsimulations/*
# tools.
#
# si8.fdf: 8-atom bulk Si (diamond cubic, a=5.43 Ang), same fixture as
# test/5-mlsimulations/3-mlelastic. A heavily rattled copy (stdev 0.3 Ang,
# generated inline below) gives the search something genuine to improve on
# -- verified live that both --algorithm basin-hopping (--steps 30, --seed
# 7) and --algorithm simulated-annealing recover almost exactly the same
# low energy a full local relax alone finds for this simple, single-basin
# system (diamond Si has one dominant global minimum at this perturbation
# scale), which is itself the correct, expected result, not a sign the
# search did nothing. The exponential cooling schedule's exact values were
# also verified by hand against the closed-form geometric-decay formula.
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


echo "--- Starting tester for stb-mlsearch (item 5.4) ---"

if ! python3 -c "import mace" 2>/dev/null; then
    echo -e "${YELLOW}Skipped entirely:${NC} the optional 'ml' extra is not installed."
    echo "Install with: pip install stb_suite[ml]  (then re-run this test)"
    exit 0
fi

rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/si8.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null

python3 -c "
from stb.core import structure_io
from pymatgen.io.ase import AseAtomsAdaptor
import numpy as np
s = structure_io.read_fdf('si8.fdf')
pmg = structure_io.to_pymatgen(s)
atoms = AseAtomsAdaptor.get_atoms(pmg)
rng = np.random.default_rng(42)
atoms.set_positions(atoms.get_positions() + rng.normal(scale=0.3, size=atoms.get_positions().shape))
pmg2 = AseAtomsAdaptor.get_structure(atoms)
new_s = structure_io.from_pymatgen(pmg2, species_meta=s.species_meta)
structure_io.write_fdf(new_s, 'si8_rattled.fdf')
"


# --- 1. Basic run on a heavily rattled structure ---
echo -e "\n--- Testing a basic search on a rattled structure (--steps 30) ---"
stb-mlsearch --file si8_rattled.fdf --steps 30 --dr 0.3 --temperature 1500 --seed 7 \
    --save-data --save-report --no-intro -o search_out > log_search.txt 2>&1
check_exit_code $? 0
check_contains "BASIN-HOPPING SEARCH" log_search.txt
check_contains "Best energy found" log_search.txt
check_contains "FINAL POLISH" log_search.txt
check_success search_out/best_structure.fdf
check_success search_out/search_history.png
check_success search_out/search_history.dat
check_success search_out/search_trajectory.xsf
check_success search_out/stb_mlsearch_report.txt
check_contains "STB-MLSEARCH REPORT" search_out/stb_mlsearch_report.txt

echo "Testing: search energy is no worse than the input (single-relax) energy"
python3 -c "
import re, sys
text = open('search_out/stb_mlsearch_report.txt').read()
m_in = re.search(r'Input energy\s*:\s*(-?[\d.]+) eV', text)
m_best = re.search(r'Best energy\s*:\s*(-?[\d.]+) eV', text)
if not (m_in and m_best):
    sys.exit(1)
e_in, e_best = float(m_in.group(1)), float(m_best.group(1))
sys.exit(0 if e_best <= e_in + 1e-6 else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} best energy found is no worse than the input energy"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} best energy found is worse than the input energy (should never happen)"
    FAIL=$((FAIL+1))
fi
rm -rf search_out


# --- 1b. Simulated-annealing algorithm ---
echo -e "\n--- Testing --algorithm simulated-annealing (--steps 100) ---"
stb-mlsearch --file si8_rattled.fdf --algorithm simulated-annealing --steps 100 \
    --temp-start 1500 --temp-end 50 --snapshot-interval 20 --seed 7 \
    --save-data --no-intro -o sa_out > log_sa.txt 2>&1
check_exit_code $? 0
check_contains "SIMULATED-ANNEALING SEARCH" log_sa.txt
check_contains "cooling 1500 K -> 50 K" log_sa.txt
check_contains "Best energy found" log_sa.txt
check_success sa_out/best_structure.fdf
check_success sa_out/search_history.png
check_success sa_out/search_history.dat
check_success sa_out/search_trajectory.xsf

echo "Testing: cooling schedule matches the exact exponential-decay formula"
python3 -c "
import numpy as np, sys
data = np.loadtxt('sa_out/search_history.dat')
steps, temps = data[:, 0], data[:, 1]
n = len(steps)
expected = 1500.0 * (50.0 / 1500.0) ** (np.arange(n) / (n - 1))
sys.exit(0 if np.allclose(temps, expected, rtol=1e-6) else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} cooling schedule matches the exact exponential formula"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} cooling schedule does not match the expected exponential decay"
    FAIL=$((FAIL+1))
fi

echo "Testing: --schedule linear gives an evenly-spaced cooling schedule"
stb-mlsearch --file si8_rattled.fdf --algorithm simulated-annealing --steps 100 \
    --temp-start 1500 --temp-end 50 --snapshot-interval 20 --schedule linear --seed 7 \
    --save-data --no-intro -o sa_linear > log_sa_linear.txt 2>&1
check_exit_code $? 0
python3 -c "
import numpy as np, sys
data = np.loadtxt('sa_linear/search_history.dat')
temps = data[:, 1]
diffs = np.diff(temps)
sys.exit(0 if np.allclose(diffs, diffs[0], rtol=1e-6) else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} --schedule linear gives evenly-spaced temperature steps"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} --schedule linear did not give evenly-spaced steps"
    FAIL=$((FAIL+1))
fi
rm -rf sa_out sa_linear


# --- 2. --seed reproducibility ---
echo -e "\n--- Testing --seed gives reproducible results ---"
stb-mlsearch --file si8_rattled.fdf --steps 10 --seed 99 --no-intro -o run1 > log_run1.txt 2>&1
stb-mlsearch --file si8_rattled.fdf --steps 10 --seed 99 --no-intro -o run2 > log_run2.txt 2>&1
E1=$(grep "Best energy found" log_run1.txt)
E2=$(grep "Best energy found" log_run2.txt)
if [ -n "$E1" ] && [ "$E1" == "$E2" ]; then
    echo -e "   -> ${GREEN}Verified:${NC} same --seed gives identical best energy across runs"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} same --seed gave different results ('$E1' vs '$E2')"
    FAIL=$((FAIL+1))
fi
rm -rf run1 run2


# --- 3. --no-pre-relax / --no-final-cell-relax ---
echo -e "\n--- Testing --no-pre-relax and --no-final-cell-relax ---"
stb-mlsearch --file si8_rattled.fdf --steps 5 --no-pre-relax --no-final-cell-relax \
    --no-intro -o norelax > log_norelax.txt 2>&1
check_exit_code $? 0
check_contains "no-pre-relax: starting the search from the structure as given" log_norelax.txt
check_contains "Skipped (--no-final-cell-relax)" log_norelax.txt
rm -rf norelax


# --- 4. Error and robustness cases ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: missing --file"
stb-mlsearch --no-intro > log_missing_file.txt 2>&1
check_exit_code $? 2

echo "Testing: nonexistent input file"
stb-mlsearch --file nope.fdf --no-intro > log_missing_input.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_input.txt

echo "Testing: --custom-model with a nonexistent file"
stb-mlsearch --file si8.fdf --custom-model does_not_exist.model --no-intro > log_custommodel.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_custommodel.txt

echo "Testing: --version"
stb-mlsearch --version > log_version.txt 2>&1
check_contains "stb-mlsearch" log_version.txt

echo "Testing: --help documents --steps, --dr, --temperature, --seed, --custom-model"
stb-mlsearch --help > log_help.txt 2>&1
check_contains "steps" log_help.txt
check_contains "\-\-dr" log_help.txt
check_contains "temperature" log_help.txt
check_contains "seed" log_help.txt
check_contains "custom-model" log_help.txt
check_contains "no-pre-relax" log_help.txt
check_contains "no-final-cell-relax" log_help.txt
check_contains "trajectory-format" log_help.txt
check_contains "algorithm" log_help.txt
check_contains "temp-start" log_help.txt
check_contains "temp-end" log_help.txt
check_contains "schedule" log_help.txt
check_contains "snapshot-interval" log_help.txt


# --- 5. Interactive path (stb-suite, shortcut 5.4) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 5.4) ---"
rm -rf interactive_out
printf '5.4\nsi8.fdf\n\nsmall\n1\n5\n0.3\n1000\n1\ninteractive_out\nn\nn\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success interactive_out/best_structure.fdf
check_success interactive_out/search_history.png


popd > /dev/null

# --- 6. Summary ---
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
