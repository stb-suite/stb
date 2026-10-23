#!/bin/bash

# --- Setup ---
# Smoke test for stb-hirshfeldAnalysis (Hirshfeld-I Stage 3: Analysis,
# item 4.20.3). Chains real stb-hirshfeldPrep + stb-hirshfeldIons (Stages
# 1-2, tested separately) runs on the ../prep fixture, fabricating each
# neutral/<species>/'s and ion/<species>/{cation,anion}/'s own .RHO via
# ../make_synthetic_rho.py (standing in for "SIESTA has already been run
# there") before running stb-hirshfeldAnalysis for real.
#
# Two kinds of synthetic density are used:
#  - raw Gaussians (sections 2-6b): exercise the PLUMBING only -- they do not
#    hold Z_val electrons, so the charges are meaningless and every reference
#    check fires (which section 3 relies on);
#  - a "physical" study (section 6c, make_physical_study): every Gaussian is
#    normalized to the electrons a real SIESTA run would hold (Z_val, Z_val-+1),
#    and each ion has the SAME shape as its neutral atom, so the Hirshfeld-I
#    reference at charge q is exactly the atom of the combined density. With
#    3.6 e- on C and 6.4 e- on O, Hirshfeld-I has an EXACT answer, q = +-0.4,
#    and spherical atoms make the point-charge dipole exact too (-0.6 e*Ang).
# Section 6d checks core.hirshfeld's numerics directly (radial profile,
# remaining-error estimate, unwrapping for the dipole).
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREP_DIR="$(cd "$FIXTURE_DIR/../prep" && pwd)"
IONS_DIR="$(cd "$FIXTURE_DIR/../ions" && pwd)"
GEN_SCRIPT="$FIXTURE_DIR/../make_synthetic_rho.py"
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

# Runs stb-hirshfeldPrep (writes combined/ + neutral/<species>/ itself,
# nothing pre-computed needed) + fabricates combined/neutral .RHOs +
# stb-hirshfeldIons (writes BOTH cation/anion per species) + fabricates
# cation/anion .RHOs -- the full chain stb-hirshfeldAnalysis needs.
make_prep_ions_and_rho() {
    rm -rf hirshfeld_study
    stb-hirshfeldPrep -s structure.fdf --calc calc.fdf -p . -O hirshfeld_study --no-intro \
        > /dev/null 2>&1
    python3 "$GEN_SCRIPT" hirshfeld_study/combined/structure.fdf \
        hirshfeld_study/combined/siesta.RHO 24 "0.8,1.4" "1.4,1.0" > /dev/null 2>&1
    python3 "$GEN_SCRIPT" hirshfeld_study/neutral/C/structure.fdf \
        hirshfeld_study/neutral/C/siesta.RHO 18 "0.8" "1.4" > /dev/null 2>&1
    python3 "$GEN_SCRIPT" hirshfeld_study/neutral/O/structure.fdf \
        hirshfeld_study/neutral/O/siesta.RHO 18 "1.4" "1.0" > /dev/null 2>&1
    stb-hirshfeldIons -O hirshfeld_study -p . --no-intro > /dev/null 2>&1
    python3 "$GEN_SCRIPT" hirshfeld_study/ions/C/cation/structure.fdf \
        hirshfeld_study/ions/C/cation/siesta.RHO 18 "0.6" "1.3" > /dev/null 2>&1
    python3 "$GEN_SCRIPT" hirshfeld_study/ions/C/anion/structure.fdf \
        hirshfeld_study/ions/C/anion/siesta.RHO 18 "1.0" "1.5" > /dev/null 2>&1
    python3 "$GEN_SCRIPT" hirshfeld_study/ions/O/cation/structure.fdf \
        hirshfeld_study/ions/O/cation/siesta.RHO 18 "1.1" "0.9" > /dev/null 2>&1
    python3 "$GEN_SCRIPT" hirshfeld_study/ions/O/anion/structure.fdf \
        hirshfeld_study/ions/O/anion/siesta.RHO 18 "1.7" "1.1" > /dev/null 2>&1
}


# Normalized "physical" study (see header): CO molecule, C sigma $1, O sigma $2
# (Ang), combined density 3.6 e- on C and 6.4 e- on O, ions with the neutral
# atom's shape. Written into $3 with a fresh Stage 1
# (cp -r of a finished study is safe: section 6e).
make_physical_study() {
    local sc="$1" so="$2" d="$3" st
    rm -rf "$d"
    stb-hirshfeldPrep -s structure.fdf --calc calc.fdf -p . -O "$d" --vacuum 12 --no-intro \
        > /dev/null 2>&1
    python3 "$GEN_SCRIPT" "$d/combined/structure.fdf" "$d/combined/siesta.RHO" 80 \
        "3.6,6.4" "$sc,$so" --electrons > /dev/null 2>&1
    for st in "neutral/C 4 $sc" "neutral/O 6 $so"; do
        set -- $st
        python3 "$GEN_SCRIPT" "$d/$1/structure.fdf" "$d/$1/siesta.RHO" 64 "$2" "$3" \
            --electrons > /dev/null 2>&1
    done
    stb-hirshfeldIons -O "$d" -p . --no-intro > /dev/null 2>&1
    for st in "ions/C/cation 3 $sc" "ions/C/anion 5 $sc" "ions/O/cation 5 $so" "ions/O/anion 7 $so"; do
        set -- $st
        python3 "$GEN_SCRIPT" "$d/$1/structure.fdf" "$d/$1/siesta.RHO" 64 "$2" "$3" \
            --electrons > /dev/null 2>&1
    done
}

# Prints field $2 (1-based) of atom $3's row in a hirshfeld_charges.dat ($1).
dat_field() {
    awk -v col="$2" -v idx="$3" '!/^#/ && $1 == idx {print $col}' "$1"
}

# check_close <label> <value> <expected> <tolerance>
check_close() {
    if python3 -c "import sys; sys.exit(0 if abs(float('$2') - float('$3')) <= float('$4') else 1)" \
        2>/dev/null; then
        echo -e "   -> ${GREEN}Verified:${NC} $1 = $2 (expected $3 +- $4)"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} $1 = '$2' (expected $3 +- $4)"
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

# --- 1. Preparation ---
echo "--- Starting tester for STB-HIRSHFELDANALYSIS (item 4.20.3) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$PREP_DIR/structure.fdf" "$PREP_DIR/calc.fdf" "$PREP_DIR/C.psf" "$PREP_DIR/O.psf" "$TEST_DIR/"
cp "$IONS_DIR/structure_mixed.fdf" "$IONS_DIR/structure_hydrogen.fdf" "$IONS_DIR/H.psf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Missing Stage-2 output guard ---
echo -e "\n--- Testing the missing-Stage-2-output guard ---"
rm -rf hirshfeld_study
stb-hirshfeldAnalysis -O hirshfeld_study --no-intro > log_no_stage2.txt 2>&1
check_exit_code $? 1
check_contains "run stb-hirshfeldPrep then stb-hirshfeldIons first" log_no_stage2.txt


# --- 3. Default run: live iteration output, convergence + per-atom report ---
echo -e "\n--- Testing default run (live iteration output, convergence, per-atom report) ---"
make_prep_ions_and_rho
stb-hirshfeldAnalysis -O hirshfeld_study --no-intro > log_default.txt 2>&1
check_exit_code $? 0

echo "Testing: each round prints a live max|Delta q| line naming the worst atom/species,"
echo "         its current charge/population + whether it's leaning cation/anion-like,"
echo "         elapsed time, and convergence status"
check_contains "ITERATING HIRSHFELD-I" log_default.txt
check_contains "Iteration  1/20: max|Delta q| = " log_default.txt
check_contains "-- atom #" log_default.txt
check_contains "-like, pop\." log_default.txt
check_contains "s elapsed --" log_default.txt

check_contains "CONVERGENCE HISTORY" log_default.txt
check_contains "Converged after" log_default.txt
check_contains "SIMPLE HIRSHFELD vs. HIRSHFELD-I" log_default.txt
check_contains "PER-SPECIES SUMMARY (HIRSHFELD-I)" log_default.txt
check_contains "Cation-like | Anion-like" log_default.txt
check_contains "Iterations run  :" log_default.txt

echo "Testing: table [4] has Idx/position(x,y,z)/Elem/Voronoi/Simple Hirshfeld/Hirshfeld-I columns"
check_contains "VORONOI vs. SIMPLE HIRSHFELD vs. HIRSHFELD-I" log_default.txt
check_contains "X(Ang)" log_default.txt
check_contains "Voronoi(e-)" log_default.txt
check_contains "Idx | Elem | X(Ang)" log_default.txt
echo "Testing: this synthetic fixture's combined/ has no real .out with a native Voronoi"
echo "         block, so the Voronoi column gracefully degrades to N/A with an [INFO] note"
check_contains "Native Voronoi charges not found" log_default.txt
check_contains "Iterating time  :" log_default.txt

echo "Testing: the fuller report -- combined-system electron count, every reference atom"
echo "         checked (electrons on the grid and in the radial profile, radius vs box),"
echo "         estimated remaining error, per-atom reference, charge balance, dipole check"
echo "         (the fixture is 3D periodic: skipped) and ALL warnings in a last section"
check_contains "\[1\] COMBINED SYSTEM" log_default.txt
check_contains "Electrons      :" log_default.txt
check_contains "\[2\] REFERENCE ATOMS" log_default.txt
check_contains "N on grid | N radial" log_default.txt
check_contains "Ratio to previous" log_default.txt
check_contains "Reference (neutral + frac ion)" log_default.txt
check_contains "\[7\] CHARGE BALANCE AND CROSS-CHECKS" log_default.txt
check_contains "Assigned to atoms (HI)" log_default.txt
check_contains "\[8\] DIPOLE CHECK" log_default.txt
check_contains "\[10\] WARNINGS" log_default.txt
echo "Testing: the synthetic Gaussians do not hold Z_val electrons -> flagged, in WARNINGS"
check_contains "\[WARNING 1\]" log_default.txt
check_contains "Files written   : none" log_default.txt


# --- 4. --tol / --max-iter flags ---
echo -e "\n--- Testing --tol/--max-iter flags ---"
stb-hirshfeldAnalysis -O hirshfeld_study --tol 0.05 --max-iter 3 --no-intro > log_tol.txt 2>&1
check_exit_code $? 0
check_contains "Tolerance      : 0.05 e-" log_tol.txt
check_contains "Max iterations : 3" log_tol.txt


# --- 4b. --ref (explicit reference file, e.g. not named with a .out extension) ---
# Also carries a real 'Voronoi Atomic Populations:' block (same file --ref
# points at doubles as the Z_val AND native-Voronoi source, exactly like a
# real combined/*.out would) -- proves table [4]'s Voronoi column's happy
# path, not just its N/A fallback (already exercised by the default run
# above, whose synthetic combined/ has no real .out at all).
echo -e "\n--- Testing --ref (Z_val detection from an explicitly-named, non-.out reference file) ---"
cat > hirshfeld_study/combined/run_log.txt << 'ZVALEOF'
atom: Called for C(Z=6)
Vna: chval, zval:    4.00000   4.00000
atom: Called for O(Z=8)
Vna: chval, zval:    6.00000   6.00000
Voronoi Atomic Populations:
Atom #   charge [q] valence [e]      Sz [e]  Species
     1    -0.050000    4.050000    0.000000  C
     2     0.050000    5.950000    0.000000  O
-------------------------------------------
 Total     0.000000                0.000000
ZVALEOF
stb-hirshfeldAnalysis -O hirshfeld_study --ref hirshfeld_study/combined/run_log.txt \
    --no-intro > log_ref.txt 2>&1
check_exit_code $? 0
check_contains "Total charge" log_ref.txt

echo "Testing: table [4]'s Voronoi column happy path -- real values, not N/A, read"
echo "         straight from the same --ref file's native 'Voronoi Atomic Populations:' block"
check_contains "-0.0500" log_ref.txt


# --- 4c. Optional outputs: --save-data and --png (off unless asked) ---
echo -e "\n--- Testing --save-data / --png ---"
rm -f hirshfeld_charges.dat hirshfeld_iterations.dat hirshfeld_charges.xyz hirshfeld_analysis.png
stb-hirshfeldAnalysis -O hirshfeld_study --no-intro > log_nodata.txt 2>&1
if [ ! -e hirshfeld_charges.dat ] && [ ! -e hirshfeld_analysis.png ]; then
    echo -e "   -> ${GREEN}Verified:${NC} no data/figure written without the flags"; PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} data/figure written without being asked"; FAIL=$((FAIL+1))
fi
stb-hirshfeldAnalysis -O hirshfeld_study --save-data --png --no-intro > log_data.txt 2>&1
check_exit_code $? 0
check_success hirshfeld_charges.dat
check_success hirshfeld_iterations.dat
check_success hirshfeld_charges.xyz
check_success hirshfeld_analysis.png
check_contains "hirshfeld_i" hirshfeld_charges.xyz
check_contains "# idx elem" hirshfeld_charges.dat
check_contains "Files written   : hirshfeld_charges.dat" log_data.txt


# --- 5. --save-report ---
echo -e "\n--- Testing --save-report ---"
rm -f hirshfeld_analysis_report.txt
stb-hirshfeldAnalysis -O hirshfeld_study --save-report --no-intro > log_savereport.txt 2>&1
check_success hirshfeld_analysis_report.txt
check_contains "STB-HIRSHFELDANALYSIS REPORT" hirshfeld_analysis_report.txt
rm -f hirshfeld_analysis_report.txt


# --- 6. End-to-end proof: atoms of the SAME species converge to OPPOSITE ion
# states (the actual payoff of writing both cation and anion per species) ---
echo -e "\n--- Testing the full pipeline on a genuinely mixed-sign species end-to-end ---"
rm -rf hirshfeld_study_mixed mixed_combined
mkdir -p mixed_combined
cp structure_mixed.fdf mixed_combined/structure.fdf
cp calc.fdf mixed_combined/
stb-hirshfeldPrep -s mixed_combined/structure.fdf --calc mixed_combined/calc.fdf -p . \
    -O hirshfeld_study_mixed --mesh-cutoff 400 --vacuum 18 --no-intro > /dev/null 2>&1
python3 "$GEN_SCRIPT" hirshfeld_study_mixed/combined/structure.fdf \
    hirshfeld_study_mixed/combined/siesta.RHO 40 "0.02,3.0,1.0" "0.5,0.6,1.0" > /dev/null 2>&1
python3 "$GEN_SCRIPT" hirshfeld_study_mixed/neutral/C/structure.fdf \
    hirshfeld_study_mixed/neutral/C/siesta.RHO 24 "0.5" "0.8" > /dev/null 2>&1
python3 "$GEN_SCRIPT" hirshfeld_study_mixed/neutral/O/structure.fdf \
    hirshfeld_study_mixed/neutral/O/siesta.RHO 24 "1.0" "0.8" > /dev/null 2>&1
stb-hirshfeldIons -O hirshfeld_study_mixed -p . --no-intro > /dev/null 2>&1
python3 "$GEN_SCRIPT" hirshfeld_study_mixed/ions/C/cation/structure.fdf \
    hirshfeld_study_mixed/ions/C/cation/siesta.RHO 24 "0.3" "0.7" > /dev/null 2>&1
python3 "$GEN_SCRIPT" hirshfeld_study_mixed/ions/C/anion/structure.fdf \
    hirshfeld_study_mixed/ions/C/anion/siesta.RHO 24 "0.9" "1.0" > /dev/null 2>&1
python3 "$GEN_SCRIPT" hirshfeld_study_mixed/ions/O/cation/structure.fdf \
    hirshfeld_study_mixed/ions/O/cation/siesta.RHO 24 "0.7" "0.8" > /dev/null 2>&1
python3 "$GEN_SCRIPT" hirshfeld_study_mixed/ions/O/anion/structure.fdf \
    hirshfeld_study_mixed/ions/O/anion/siesta.RHO 24 "1.3" "0.9" > /dev/null 2>&1
stb-hirshfeldAnalysis -O hirshfeld_study_mixed --no-intro > log_mixed_analysis.txt 2>&1
check_exit_code $? 0
echo "Testing: the converged per-species table shows species 'C' with BOTH a cation-like"
echo "         and an anion-like atom, highlighted -- exactly what a single-sign-per-species"
echo "         approach could never represent correctly"
check_contains "C    | 2 | 1           | 1          " log_mixed_analysis.txt
check_contains "atoms genuinely leaning toward BOTH ion states" log_mixed_analysis.txt


# --- 6b. End-to-end proof: hydrogen (Z_val=1) -- its cation is a bare proton, so
# Stage 2 writes NO ions/H/cation/ folder (NetCharge +1.0 leaves 0 electrons, an
# SCF SIESTA cannot converge) and Stage 3 substitutes the analytically exact
# zero-density profile instead of aborting on a missing .RHO ---
echo -e "\n--- Testing the full pipeline with hydrogen (bare-proton cation, zero density) ---"
rm -rf hirshfeld_study_h h_combined
mkdir -p h_combined
cp structure_hydrogen.fdf h_combined/structure.fdf
cp calc.fdf h_combined/
stb-hirshfeldPrep -s h_combined/structure.fdf --calc h_combined/calc.fdf -p . \
    -O hirshfeld_study_h --no-intro > /dev/null 2>&1
python3 "$GEN_SCRIPT" hirshfeld_study_h/combined/structure.fdf \
    hirshfeld_study_h/combined/siesta.RHO 24 "0.8,0.5" "1.4,0.9" > /dev/null 2>&1
python3 "$GEN_SCRIPT" hirshfeld_study_h/neutral/C/structure.fdf \
    hirshfeld_study_h/neutral/C/siesta.RHO 18 "0.8" "1.4" > /dev/null 2>&1
python3 "$GEN_SCRIPT" hirshfeld_study_h/neutral/H/structure.fdf \
    hirshfeld_study_h/neutral/H/siesta.RHO 18 "0.5" "0.9" > /dev/null 2>&1
stb-hirshfeldIons -O hirshfeld_study_h -p . --no-intro > /dev/null 2>&1
# Deliberately NO cation .RHO for H -- there is no ions/H/cation/ folder to put one in.
python3 "$GEN_SCRIPT" hirshfeld_study_h/ions/C/cation/structure.fdf \
    hirshfeld_study_h/ions/C/cation/siesta.RHO 18 "0.6" "1.3" > /dev/null 2>&1
python3 "$GEN_SCRIPT" hirshfeld_study_h/ions/C/anion/structure.fdf \
    hirshfeld_study_h/ions/C/anion/siesta.RHO 18 "1.0" "1.5" > /dev/null 2>&1
python3 "$GEN_SCRIPT" hirshfeld_study_h/ions/H/anion/structure.fdf \
    hirshfeld_study_h/ions/H/anion/siesta.RHO 18 "0.8" "1.2" > /dev/null 2>&1
stb-hirshfeldAnalysis -O hirshfeld_study_h --no-intro > log_hydrogen_analysis.txt 2>&1
check_exit_code $? 0
echo "Testing: H's cation is announced as the analytical zero-density reference (not read"
echo "         from a folder, not a silent omission), C's is still read from its folder,"
echo "         and the iteration converges to a full report"
check_contains "H (cation)  : analytical zero density" log_hydrogen_analysis.txt
check_contains "C (cation)  : read from" log_hydrogen_analysis.txt
check_contains "H (anion)   : read from" log_hydrogen_analysis.txt
check_contains "Converged after" log_hydrogen_analysis.txt
check_contains "Total charge" log_hydrogen_analysis.txt


echo "Testing: an atom past |q| = 1 is flagged -- its reference is capped at the +-1 ion"
check_contains "the reference is capped at the +-1 ion" log_mixed_analysis.txt
check_contains "(capped)" log_mixed_analysis.txt


# --- 6c. Physical study: Hirshfeld-I against its EXACT answer ---
echo -e "\n--- Testing the normalized physical study (exact answer q = +-0.4) ---"
make_physical_study 1.1 1.0 hirshfeld_study_phys
# A SIESTA log for combined/: Z_val, SCF, native Hirshfeld/Voronoi and the
# dipole of this density, -0.6 e*Ang along z = -1.133836 a.u. = -2.881923 D.
cat > hirshfeld_study_phys/combined/calc.out << 'LOGEOF'
atom: Called for C(Z=6)
Vna: chval, zval:    4.00000   4.00000
atom: Called for O(Z=8)
Vna: chval, zval:    6.00000   6.00000
SCF cycle converged after 7 iterations
Hirshfeld Atomic Populations:
Atom #   charge [q] valence [e]      Sz [e]  Species
     1     0.150000    3.850000    0.000000  C
     2    -0.150000    6.150000    0.000000  O
-------------------------------------------
 Total     0.000000                0.000000

Voronoi Atomic Populations:
Atom #   charge [q] valence [e]      Sz [e]  Species
     1     0.120000    3.880000    0.000000  C
     2    -0.120000    6.120000    0.000000  O
-------------------------------------------
 Total     0.000000                0.000000
siesta: Electric dipole (a.u.)  =    0.000000    0.000000   -1.133836
siesta: Electric dipole (Debye) =    0.000000    0.000000   -2.881923
LOGEOF

echo "Testing: default --tol -- this density converges slowly (contraction ratio ~0.64),"
echo "         so the last step is under --tol while the charges are not yet: the report"
echo "         estimates the remaining error and warns"
rm -f hirshfeld_charges.dat
stb-hirshfeldAnalysis -O hirshfeld_study_phys --save-data --no-intro > log_phys_default.txt 2>&1
check_exit_code $? 0
check_contains "Converged after" log_phys_default.txt
check_contains "Contraction ratio      : 0\.6" log_phys_default.txt
check_contains "Est. remaining error" log_phys_default.txt
check_contains "from the converged Hirshfeld-I values" log_phys_default.txt
check_contains "rerun with --tol" log_phys_default.txt
q_default=$(dat_field hirshfeld_charges.dat 10 1)
est_remaining=$(grep -o "Est. remaining error   : [0-9.eE+-]*" log_phys_default.txt | awk '{print $NF}')

echo "Testing: tight --tol -- the exact answer, a clean report"
rm -f hirshfeld_charges.dat
stb-hirshfeldAnalysis -O hirshfeld_study_phys --tol 1e-5 --max-iter 100 --save-data \
    --no-intro > log_phys_tight.txt 2>&1
check_exit_code $? 0
q_c=$(dat_field hirshfeld_charges.dat 10 1)
q_o=$(dat_field hirshfeld_charges.dat 10 2)
q0_c=$(dat_field hirshfeld_charges.dat 9 1)
check_close "Hirshfeld-I q(C)" "$q_c" 0.4 0.005
check_close "Hirshfeld-I q(O)" "$q_o" -0.4 0.005
echo "Testing: simple Hirshfeld (neutral references) understates the transfer"
check_close "simple Hirshfeld q(C), well below 0.4" "$q0_c" 0.17 0.12

echo "Testing: the estimated remaining error at the default --tol matches the real one"
echo "         (|q_default - q_converged|) within 40 %"
actual_remaining=$(python3 -c "print('%.4g' % abs(float('$q_default') - float('$q_c')))")
check_close "remaining error estimated (vs actual)" "$est_remaining" "$actual_remaining" \
    "$(python3 -c "print('%.2g' % (0.4 * float('$actual_remaining')))")"

echo "Testing: every reference holds Z_val / Z_val -+ 1 electrons; Z_val, SCF and the native"
echo "         SIESTA columns come from the log -- and there is nothing to warn about"
check_contains "C    | neutral | 4        | 4.0000" log_phys_tight.txt
check_contains "O    | anion   | 7        | 7.0000" log_phys_tight.txt
check_contains "Electrons      : 10.0000 e-" log_phys_tight.txt
check_contains "detected from" log_phys_tight.txt
check_contains "SCF            : converged after 7 iterations" log_phys_tight.txt
check_contains "+0.1500" log_phys_tight.txt
check_contains "+0.1200" log_phys_tight.txt
check_contains "Simple Hirshfeld here vs SIESTA's own" log_phys_tight.txt
check_contains "unassigned: [+-]0.0000 e-" log_phys_tight.txt
check_contains "Periodicity    : 0D" log_phys_tight.txt
check_contains "0.40 cation" log_phys_tight.txt
check_contains "0.40 anion" log_phys_tight.txt
check_not_contains "WARNING 1" log_phys_tight.txt
check_contains "  none" log_phys_tight.txt

echo "Testing: dipole -- spherical atoms make the Hirshfeld-I point charges reproduce the"
echo "         density's dipole (2.882 D); SIESTA's value is read from the log"
dip_line=$(sed 's/\x1b\[[0-9;]*m//g' log_phys_tight.txt | grep "^|mu| ")
dip_hi=$(echo "$dip_line" | awk -F'|' '{print $(NF-1)}' | tr -d ' ')
dip_siesta=$(echo "$dip_line" | awk -F'|' '{print $NF}' | tr -d ' ')
check_close "|mu| SIESTA (from the log)" "$dip_siesta" 2.882 0.001
check_close "|mu| Hirshfeld-I point charges" "$dip_hi" 2.882 0.03
check_contains "Component (Debye) | Voronoi (SIESTA) | Hirshfeld (SIESTA)" log_phys_tight.txt

echo "Testing: --max-iter too small -> not converged, said in WARNINGS"
stb-hirshfeldAnalysis -O hirshfeld_study_phys --max-iter 1 --no-intro > log_phys_maxiter.txt 2>&1
check_exit_code $? 0
check_contains "NOT converged within --max-iter" log_phys_maxiter.txt
check_contains "Did not confirm convergence within --max-iter 1" log_phys_maxiter.txt

echo "Testing: --show with no display falls back to saving the figure"
rm -f hirshfeld_analysis.png
env -u DISPLAY -u WAYLAND_DISPLAY MPLBACKEND=Agg \
    stb-hirshfeldAnalysis -O hirshfeld_study_phys --show --no-intro > log_phys_show.txt 2>&1
check_exit_code $? 0
check_contains "needs a display" log_phys_show.txt
check_success hirshfeld_analysis.png

echo "Testing: an anion whose density reaches its box edge is flagged (O anion widened to"
echo "         sigma 1.8 Ang in its 12 Ang box)"
python3 "$GEN_SCRIPT" hirshfeld_study_phys/ions/O/anion/structure.fdf \
    hirshfeld_study_phys/ions/O/anion/siesta.RHO 64 "7" "1.8" --electrons > /dev/null 2>&1
stb-hirshfeldAnalysis -O hirshfeld_study_phys --no-intro > log_phys_box.txt 2>&1
check_exit_code $? 0
check_contains "O anion: 99.9 % of its density lies within" log_phys_box.txt
check_contains "rerun Stage 1 with a larger --vacuum" log_phys_box.txt


# --- 6e. A copied study reads ITS OWN combined/, never the original's ---
# The manifests used to hold absolute paths: a `cp -r` of a study then read
# the ORIGINAL combined/ density with the copy's references, silently (in a
# real case +1.25 e- instead of the right +0.4). The copy here gets a density
# with twice the electrons (20 e-), so any mix-up shows in the electron count.
echo -e "\n--- Testing a copied study (relative paths; old absolute manifests) ---"
make_physical_study 1.1 1.0 study_orig
rm -rf study_copy
cp -r study_orig study_copy
python3 "$GEN_SCRIPT" study_copy/combined/structure.fdf study_copy/combined/siesta.RHO 80 \
    "7.2,12.8" "1.1,1.0" --electrons > /dev/null 2>&1
stb-hirshfeldAnalysis -O study_orig --no-intro > log_copy_orig.txt 2>&1
stb-hirshfeldAnalysis -O study_copy --no-intro > log_copy_new.txt 2>&1
check_exit_code $? 0
echo "Testing: relative manifests -- each study reads its own density"
check_contains "Electrons      : 10.0000 e-" log_copy_orig.txt
check_contains "Electrons      : 20.0000 e-" log_copy_new.txt
check_contains "Structure      : study_copy/combined/structure.fdf" log_copy_new.txt

echo "Testing: an OLD manifest (absolute paths, written before the fix) copied along with"
echo "         its study: the copy's own combined/ wins, and the report says so"
python3 - "$PWD/study_orig" << 'ABSEOF'
import json, os, sys
orig = sys.argv[1]
m = json.load(open("study_copy/hirshfeld_manifest.json"))
m["combined_dir"] = os.path.join(orig, "combined")
m["combined_structure"] = os.path.join(orig, "combined", "structure.fdf")
json.dump(m, open("study_copy/hirshfeld_manifest.json", "w"), indent=2)
im = json.load(open("study_copy/hirshfeld_ions_manifest.json"))
im["hirshfeld_manifest"] = os.path.join(orig, "hirshfeld_manifest.json")
json.dump(im, open("study_copy/hirshfeld_ions_manifest.json", "w"), indent=2)
ABSEOF
stb-hirshfeldAnalysis -O study_copy --no-intro > log_copy_old.txt 2>&1
check_exit_code $? 0
check_contains "Electrons      : 20.0000 e-" log_copy_old.txt
check_contains "outside this study (the study was copied or moved)" log_copy_old.txt
stb-hirshfeldIons -O study_copy -p . --no-intro > log_copy_old_ions.txt 2>&1
check_contains "outside this study (the study was copied or moved)" log_copy_old_ions.txt
check_contains "Combined structure : study_copy/combined/structure.fdf" log_copy_old_ions.txt

echo "Testing: the old manifest's foreign combined/ is used only when the copy has none --"
echo "         and then it is a WARNING"
python3 - "$PWD/study_orig" << 'ABSEOF'
import json, os, sys
orig = sys.argv[1]
m = json.load(open("study_copy/hirshfeld_manifest.json"))
m["combined_dir"] = os.path.join(orig, "combined")
m["combined_structure"] = os.path.join(orig, "combined", "structure.fdf")
json.dump(m, open("study_copy/hirshfeld_manifest.json", "w"), indent=2)
ABSEOF
rm -rf study_copy/combined
stb-hirshfeldAnalysis -O study_copy --no-intro > log_copy_nocombined.txt 2>&1
check_exit_code $? 0
check_contains "Electrons      : 10.0000 e-" log_copy_nocombined.txt
check_contains "in ANOTHER study folder -- is read instead" log_copy_nocombined.txt
rm -rf study_orig study_copy


# --- 6d. core.hirshfeld numerics, called directly ---
echo -e "\n--- Testing core.hirshfeld numerics (radial profile, remaining error, unwrap) ---"
python3 "$FIXTURE_DIR/check_core_numerics.py" > log_core.txt 2>&1
check_exit_code $? 0
check_contains "OK profile electron count (default bins)" log_core.txt
check_contains "OK no empty bins (fine bins)" log_core.txt
check_contains "OK profile radius" log_core.txt
check_contains "OK remaining-error estimate" log_core.txt
check_contains "OK unwrap across the boundary" log_core.txt


# --- 7. Interactive path (stb-suite, shortcut 4.20.3) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.20.3) ---"
rm -f hirshfeld_analysis.png
printf '4.20.3\nhirshfeld_study\n\n\n\nn\nn\ny\nn\n\n0\n' | stb-suite > log_interactive.txt 2>&1
check_contains "Hirshfeld-I Analysis (Stage 3) complete" log_interactive.txt
check_contains "Optional outputs:" log_interactive.txt
check_success hirshfeld_analysis.png


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
