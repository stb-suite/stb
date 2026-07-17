#!/bin/bash

# --- Setup ---
# Smoke test for stb-gqcaAnalysis (GQCA Stage 2: Mass-Action Analysis,
# item 4.15.2). Chains a real stb-gqca (Stage 1, tested separately) run
# on the cubic.fdf fixture, then fabricates calc.out files with known
# FreeEng values in each cluster_n{0,1,2} folder to exercise 3 analytic
# gates recomputed independently in Python (same discipline as HER/OER's
# own analysis test.sh files):
#   1. Pure end-member limits (x=0 -> [1,0,0], x=1 -> [0,0,1]).
#   2. Equal cluster energies reduce EXACTLY to the random/binomial
#      mixing closed form (Delta_j == 0 for all j regardless of any
#      ambiguity in the mixing-energy formula itself).
#   3. An ordering-favoring energy set drives x_1 above the random value
#      and SRO away from 0, with the expected (ordering) sign -- a
#      directional check, not an exact-value one, since the SRO formula
#      is itself flagged [UNVERIFIED].
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

# Writes calc.out with a known FreeEng into cluster_n$1 (1 sublattice
# site: n0/n2) or cluster_n1 (2 sublattice sites), a converged/low-force
# SCF block so report_quality_diagnostics stays silent.
write_energy() {
    local dir="$1"
    local free_eng="$2"
    cat > "gqca_study/cluster_n${dir}/calc.out" <<EOF
siesta: SCF Convergence by DM criterion
SCF cycle converged after 12 iterations
siesta: FreeEng =        ${free_eng}
siesta: Atomic forces (eV/Ang):
     1    0.001    0.001    0.001
     Max    0.001
EOF
}


# --- 1. Preparation ---
echo "--- Starting tester for STB-GQCAAnalysis stage 2: mass-action analysis (item 4.15.2) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$PREP_DIR/cubic.fdf" "$TEST_DIR/"
cp "$PREP_DIR/calc.fdf" "$TEST_DIR/"
for sym in Na K; do
    echo "# placeholder pseudopotential" > "$TEST_DIR/$sym.psf"
done
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null

echo -e "\n--- Running stb-gqca (Stage 1) once on cubic.fdf to build cluster folders ---"
rm -rf gqca_study
stb-gqca -f cubic.fdf --sublattice Na --species-b K -c calc.fdf -p . -O gqca_study --no-intro \
    > /dev/null 2>&1
check_success gqca_study/cluster_n0/structure.fdf
check_success gqca_study/cluster_n1/structure.fdf
check_success gqca_study/cluster_n2/structure.fdf


# --- 2. Missing Stage-1-results guard ---
echo -e "\n--- Testing the missing-results guard ---"
rm -rf gqca_study_empty
mkdir -p gqca_study_empty
stb-gqcaAnalysis --directory gqca_study_empty --temp 500 --no-intro > log_no_results.txt 2>&1
check_exit_code $? 1
check_contains "run stb-gqca" log_no_results.txt


# --- 3. Gate 1 + Gate 2: equal cluster energies (Delta_j == 0 for all j) ---
echo -e "\n--- Testing Gate 1 (pure end-member limits) + Gate 2 (random/binomial mixing) ---"
# cluster_n0/n2 have 1 sublattice site each; cluster_n1 has 2 -- FreeEng
# chosen so every folder normalizes to the SAME per-site energy (-50.0
# eV/site), making Delta_j == 0 for all j exactly.
write_energy 0 "-50.000000"
write_energy 2 "-50.000000"
write_energy 1 "-100.000000"

stb-gqcaAnalysis --directory gqca_study --temp 500 --x-min 0.0 --x-max 1.0 --x-points 5 \
    -o gate12 --no-intro > log_gate12.txt 2>&1
check_exit_code $? 0
check_success gqca_study/gate12.csv
check_contains "x=0 pure-A limit (x_j = \[1,0,0\]) : " log_gate12.txt
check_contains "x=1 pure-B limit (x_j = \[0,0,1\]) : " log_gate12.txt

echo "Testing: miscibility-gap section runs and documents its own structural-convexity limitation"
check_contains "MISCIBILITY GAP" log_gate12.txt
check_contains "STRUCTURAL LIMITATION" log_gate12.txt
check_contains "fully miscible" log_gate12.txt

python3 -c "
import csv, math

kB = 8.617333262e-5

with open('gqca_study/gate12.csv') as f:
    rows = list(csv.DictReader(f))

def close(a, b, tol):
    return abs(a - b) < tol

row0 = next(r for r in rows if close(float(r['x']), 0.0, 1e-9))
assert close(float(row0['x_0']), 1.0, 1e-6), row0
assert close(float(row0['x_1']), 0.0, 1e-6), row0
assert close(float(row0['x_2']), 0.0, 1e-6), row0

row1 = next(r for r in rows if close(float(r['x']), 1.0, 1e-9))
assert close(float(row1['x_2']), 1.0, 1e-6), row1
assert close(float(row1['x_0']), 0.0, 1e-6), row1
assert close(float(row1['x_1']), 0.0, 1e-6), row1

# Gate 2: at x=0.5, equal cluster energies (Delta_j == 0) must reduce
# exactly to the binomial(N=2, x) random-mixing distribution:
# eta = x/(1-x) = 1, x_j = [0.25, 0.5, 0.25], H_mix = 0, SRO = 0. S_mix
# is recomputed independently here from the same closed form (not
# imported from stb's own solver) as the gate.
row_half = next(r for r in rows if close(float(r['x']), 0.5, 1e-9))
x_j = [float(row_half['x_0']), float(row_half['x_1']), float(row_half['x_2'])]
assert close(x_j[0], 0.25, 1e-6), x_j
assert close(x_j[1], 0.50, 1e-6), x_j
assert close(x_j[2], 0.25, 1e-6), x_j
assert close(float(row_half['H_mix_eV']), 0.0, 1e-6), row_half
assert close(float(row_half['SRO']), 0.0, 1e-6), row_half

g = [1, 2, 1]
s_expected = -kB * sum(xj * math.log(xj / gj) for xj, gj in zip(x_j, g))
assert close(float(row_half['S_mix_eV_per_K']), s_expected, 1e-9), \
    (row_half['S_mix_eV_per_K'], s_expected)

print('OK')
" > log_gate12_check.txt 2>&1
check_contains "OK" log_gate12_check.txt


# --- 4. Gate 3: ordering-favoring energy set ---
echo -e "\n--- Testing Gate 3 (ordering-favoring energy: x_1 > random value, SRO != 0) ---"
# e_1 (AB, per-site) pushed well below the e_0/e_2 linear interpolation
# -- Delta_1 = -2.0 eV/site, strongly ordering-favoring at 500 K
# (kB*T ~ 0.043 eV, so a 2 eV bias essentially fully saturates x_1 -> 1).
write_energy 0 "-50.000000"
write_energy 2 "-50.000000"
write_energy 1 "-104.000000"

stb-gqcaAnalysis --directory gqca_study --temp 500 --x-min 0.0 --x-max 1.0 --x-points 5 \
    -o gate3 --no-intro > log_gate3.txt 2>&1
check_exit_code $? 0
check_success gqca_study/gate3.csv

python3 -c "
import csv

with open('gqca_study/gate3.csv') as f:
    rows = list(csv.DictReader(f))

row_half = next(r for r in rows if abs(float(r['x']) - 0.5) < 1e-9)
x1 = float(row_half['x_1'])
sro = float(row_half['SRO'])
# Ordering-favoring Delta_1 < 0 must push x_1 (the AB/mixed-cluster
# fraction) ABOVE its random-mixing value of 0.5 at x=0.5 (gate 2's own
# result), and SRO must be non-zero with the ordering sign (negative,
# per sro_parameter's own documented convention).
assert x1 > 0.5, f'x_1={x1} not above the random value 0.5'
assert sro < -1e-6, f'SRO={sro} not negative (ordering-favoring) as expected'
print('OK')
" > log_gate3_check.txt 2>&1
check_contains "OK" log_gate3_check.txt


# --- 5. Temperature sweep ---
echo -e "\n--- Testing --temp-min/--temp-max/--temp-points (sweep) ---"
stb-gqcaAnalysis --directory gqca_study --temp-min 300 --temp-max 1200 --temp-points 4 \
    --x-points 3 -o sweep --no-intro > log_sweep.txt 2>&1
check_exit_code $? 0
check_success gqca_study/sweep.csv
python3 -c "
import csv
with open('gqca_study/sweep.csv') as f:
    rows = list(csv.DictReader(f))
temps = sorted(set(float(r['T_K']) for r in rows))
assert len(temps) == 4, temps
assert abs(temps[0] - 300.0) < 1e-6 and abs(temps[-1] - 1200.0) < 1e-6, temps
assert len(rows) == 4 * 3, len(rows)
print('OK')
" > log_sweep_check.txt 2>&1
check_contains "OK" log_sweep_check.txt


# --- 6. Error cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: --temp and --temp-min/--temp-max are mutually exclusive"
stb-gqcaAnalysis --directory gqca_study --temp 500 --temp-min 300 --temp-max 800 \
    --no-intro > log_mutex.txt 2>&1
check_exit_code $? 1
check_contains "mutually exclusive" log_mutex.txt

echo "Testing: --version"
stb-gqcaAnalysis --version > log_version.txt 2>&1
check_contains "stb-gqcaAnalysis" log_version.txt

echo "Testing: --help documents --temp/--x-min/--output and the [UNVERIFIED] entropy note"
stb-gqcaAnalysis --help > log_help.txt 2>&1
check_contains "temp-min" log_help.txt
check_contains "x-min" log_help.txt
check_contains "output" log_help.txt
check_contains "UNVERIFIED" log_help.txt


# --- 7. Interactive path (stb-suite, shortcut 4.15.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.15.2) ---"

echo "Testing: navigate 4.15.2 -> gqca_study -> fixed T -> quit"
{
  echo "4.15.2"
  echo "gqca_study"     # directory
  echo ""                # output filename (default calc.out)
  echo "1"                 # single fixed temperature
  echo "500"                 # T
  echo ""                     # advanced settings (default N)
  echo ""                      # press enter to continue
  echo "0"                      # quit stage submenu
} | stb-suite > log_menu.txt 2>&1
check_success gqca_study/gqca_results.csv
check_success gqca_study/gqca_stage2.txt


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
