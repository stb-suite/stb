#!/bin/bash

# --- Setup ---
# Smoke test for stb-mlgcmc (ML Monte Carlo / GCMC, item 5.9)
#
# Needs the optional 'ml' extra (pip install stb_suite[ml] -- PyTorch +
# mace-torch). The whole file is skipped with a clear message if `mace`
# isn't importable, same gating pattern as the other test/5-mlsimulations/*
# tools.
#
# graphene.fdf: 2-atom primitive graphene monolayer (vacuum along c), same
# fixture as test/5-mlsimulations/3-mlelastic -- a simple, fast physisorption
# host for Ar. Verified live: --ensemble canonical (fixed N=2) keeps the
# adsorbate count exactly fixed with a sensible (~50-80%) move acceptance
# rate; the insertion/deletion Metropolis formulas themselves were verified
# BEFORE ever touching MACE, against the exact analytic ideal-gas
# grand-canonical result (<N> = V*z, variance(N) ~= mean(N), the Poisson
# signature) using a toy zero-interaction test.
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


echo "--- Starting tester for stb-mlgcmc (item 5.9) ---"

if ! python3 -c "import mace" 2>/dev/null; then
    echo -e "${YELLOW}Skipped entirely:${NC} the optional 'ml' extra is not installed."
    echo "Install with: pip install stb_suite[ml]  (then re-run this test)"
    exit 0
fi

rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/graphene.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 1. Canonical ensemble: fixed N ---
echo -e "\n--- Testing --ensemble canonical (fixed N=2 Ar atoms) ---"
stb-mlgcmc --host graphene.fdf --adsorbate Ar --ensemble canonical --n-initial 2 \
    --steps 200 --equilibration-steps 50 --dr 0.5 --seed 1 \
    --save-data --save-report --no-intro -o canon_out > log_canon.txt 2>&1
check_exit_code $? 0
check_contains "MONTE CARLO" log_canon.txt
check_contains "move moves" log_canon.txt
check_success canon_out/mc_history.png
check_success canon_out/mc_history.dat
check_success canon_out/final_config.xsf
check_success canon_out/stb_mlgcmc_report.txt
check_contains "STB-MLGCMC REPORT" canon_out/stb_mlgcmc_report.txt

echo "Testing: N stays exactly fixed at 2 throughout the canonical run"
python3 -c "
import numpy as np, sys
data = np.loadtxt('canon_out/mc_history.dat')
n_col = data[:, 0]
sys.exit(0 if np.all(n_col == 2) else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} N stays exactly fixed at 2 (no insert/delete in canonical)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} N changed during a canonical (fixed-N) run"
    FAIL=$((FAIL+1))
fi

echo "Testing: move acceptance rate is nonzero and sensible (10-95%)"
python3 -c "
import re, sys
text = open('canon_out/stb_mlgcmc_report.txt').read()
m = re.search(r'move moves: (\d+)/(\d+) accepted \(([\d.]+)%\)', text)
if not m:
    sys.exit(1)
rate = float(m.group(3))
sys.exit(0 if 10.0 < rate < 95.0 else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} move acceptance rate is sensible"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} move acceptance rate is outside the sensible range"
    FAIL=$((FAIL+1))
fi
rm -rf canon_out


# --- 2. Grand-canonical: single mu ---
# --mu -0.30 is deliberately calibrated (not -0.5, which gives ~zero
# loading for this box -- see the module docstring's live-verified
# example): mu is extremely sensitive (~50x change in loading per 0.1 eV
# at 300 K for a light adsorbate), so an arbitrary mu mostly just gives
# "no adsorption" or "runaway overcrowding" for a box this small. -0.30
# was found by solving for an O(1) ideal-gas reference loading (V*z) for
# this exact 105 Ang^3 graphene cell/Ar/300K combination.
echo -e "\n--- Testing --ensemble grand-canonical (single --mu, calibrated) ---"
stb-mlgcmc --host graphene.fdf --adsorbate Ar --mu -0.30 \
    --steps 150 --equilibration-steps 50 --insertion-zmin 0.55 --insertion-zmax 0.95 \
    --seed 3 --no-intro -o gc_out > log_gc.txt 2>&1
check_exit_code $? 0
check_contains "Ideal-gas reference loading" log_gc.txt
check_contains "insert moves" log_gc.txt
check_contains "delete moves" log_gc.txt
check_contains "<N> =" log_gc.txt
check_success gc_out/mc_history.png
check_success gc_out/final_config.xsf
rm -rf gc_out


# --- 2b. --mu-scan: loading increases with mu (2 well-separated,
#     calibrated points -- more robust than a fine-grained scan, which is
#     noisy at smoke-test step counts; see module docstring for how these
#     were chosen) ---
echo -e "\n--- Testing --mu-scan: loading increases with mu (2-point isotherm) ---"
stb-mlgcmc --host graphene.fdf --adsorbate Ar --mu-scan -1.0 -0.28 \
    --steps 100 --equilibration-steps 20 --insertion-zmin 0.55 --insertion-zmax 0.95 \
    --seed 9 --save-data --no-intro -o twopoint > log_twopoint.txt 2>&1
check_exit_code $? 0
check_contains "far from the O(0.1-10) range" log_twopoint.txt
check_success twopoint/isotherm.png
check_success twopoint/isotherm.dat

echo "Testing: mean loading is (near) zero at mu=-1.0 and clearly nonzero at mu=-0.28"
python3 -c "
import numpy as np, sys
data = np.loadtxt('twopoint/isotherm.dat')
mus, means = data[:, 0], data[:, 1]
lo = means[np.argmin(mus)]
hi = means[np.argmax(mus)]
sys.exit(0 if (lo < 0.5 and hi > lo) else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} loading is near-zero at the very negative mu and clearly higher at the calibrated one"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} loading did not increase with mu as expected"
    FAIL=$((FAIL+1))
fi
rm -rf twopoint


# --- 3. Error and robustness cases ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: grand-canonical without --mu/--mu-scan is rejected"
stb-mlgcmc --host graphene.fdf --adsorbate Ar --no-intro > log_nomu.txt 2>&1
check_exit_code $? 2
check_contains "needs --mu or --mu-scan" log_nomu.txt

echo "Testing: canonical without --n-initial > 0 is rejected"
stb-mlgcmc --host graphene.fdf --adsorbate Ar --ensemble canonical --no-intro > log_non0.txt 2>&1
check_exit_code $? 2
check_contains "needs --n-initial > 0" log_non0.txt

echo "Testing: --equilibration-steps >= --steps is rejected (would silently collect no statistics)"
stb-mlgcmc --host graphene.fdf --adsorbate Ar --ensemble canonical --n-initial 2 --steps 150 \
    --no-intro > log_eqcheck.txt 2>&1
check_exit_code $? 2
check_contains "equilibration-steps must be less than" log_eqcheck.txt

echo "Testing: unknown adsorbate element symbol"
stb-mlgcmc --host graphene.fdf --adsorbate Xx --mu -0.5 --no-intro > log_badelem.txt 2>&1
check_exit_code $? 1
check_contains "Unknown element symbol" log_badelem.txt

echo "Testing: missing --host"
stb-mlgcmc --adsorbate Ar --mu -0.5 --no-intro > log_missing_file.txt 2>&1
check_exit_code $? 2

echo "Testing: nonexistent host file"
stb-mlgcmc --host nope.fdf --adsorbate Ar --mu -0.5 --no-intro > log_missing_host.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_host.txt

echo "Testing: --custom-model with a nonexistent file"
stb-mlgcmc --host graphene.fdf --adsorbate Ar --mu -0.5 --custom-model does_not_exist.model \
    --no-intro > log_custommodel.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_custommodel.txt

echo "Testing: --version"
stb-mlgcmc --version > log_version.txt 2>&1
check_contains "stb-mlgcmc" log_version.txt

echo "Testing: --help documents --ensemble, --mu, --mu-scan, --adsorbate, --custom-model"
stb-mlgcmc --help > log_help.txt 2>&1
check_contains "ensemble" log_help.txt
check_contains "\-\-mu" log_help.txt
check_contains "mu-scan" log_help.txt
check_contains "adsorbate" log_help.txt
check_contains "custom-model" log_help.txt
check_contains "insertion-zmin" log_help.txt


# --- 4. Interactive path (stb-suite, shortcut 5.9) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 5.9) ---"
rm -rf interactive_out
printf '5.9\ngraphene.fdf\nAr\n2\n2\n300\n150\n30\n\ninteractive_out\nn\nn\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success interactive_out/mc_history.png


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
