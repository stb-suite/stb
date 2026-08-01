#!/bin/bash

# --- Setup ---
# Smoke test for stb-cohesiveAnalysis (Cohesive Energy Analysis, item 4.3.2)
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$FIXTURE_DIR/test_files"

# Output colors
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

PASS=0
FAIL=0

# --- Check helpers ---

# Checks that file $1 exists and is not empty
check_success() {
    if [ -s "$1" ]; then
        echo -e " ... ${GREEN}OK${NC} (file '$1' created)"
        PASS=$((PASS+1))
    else
        echo -e " ... ${RED}FAIL${NC} (file '$1' was not created)"
        FAIL=$((FAIL+1))
    fi
}

# Checks that file $2 contains (grep -q) pattern $1
check_contains() {
    if grep -q -- "$1" "$2" 2>/dev/null; then
        echo -e "   -> ${GREEN}Verified:${NC} '$1' found in '$2'"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} '$1' NOT found in '$2'"
        FAIL=$((FAIL+1))
    fi
}

# Checks that $2 does NOT contain (grep -q) pattern $1
check_not_contains() {
    if grep -q -- "$1" "$2" 2>/dev/null; then
        echo -e "   -> ${RED}Failed:${NC} '$1' found in '$2' (should not be there)"
        FAIL=$((FAIL+1))
    else
        echo -e "   -> ${GREEN}Verified:${NC} '$1' absent from '$2'"
        PASS=$((PASS+1))
    fi
}

# Checks that $1 (actual exit code) equals $2 (expected exit code)
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
echo "--- Starting tester for STB-CohesiveAnalysis (item 4.3.2) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp -r "$FIXTURE_DIR/structure" "$TEST_DIR/"
cp -r "$FIXTURE_DIR/atoms" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Real fixture: single-atom-in-cell synthetic FreeEng values (no
#     atoms_bsse/ here -- uncorrected-only path) ---
echo -e "\n--- Testing default run (synthetic 1-atom fixture, no BSSE data) ---"
rm -f cohesive_results.dat cohesive_correction.png
stb-cohesiveAnalysis -o calc.out -d . --save-report --no-intro > log_default.txt 2>&1
check_exit_code $? 0
check_success cohesive_results.dat

echo "Verifying the cohesive energy arithmetic (E_bulk - E_isolated) / N_atoms"
check_contains "Total atoms | 1" log_default.txt
check_contains "Uncorrected.*-450.1111.*-450.1111" log_default.txt
check_contains "Uncorrected.*-450.1111.*-450.1111" cohesive_results.dat

echo "Verifying the BSSE caveat note appears when atoms_bsse/ isn't present"
check_contains "NOTE. This result is NOT corrected for BSSE" log_default.txt
check_not_contains "atoms_bsse" log_default.txt

echo "Verifying [4] CORRECTION PLOT is skipped (no BSSE data) and no PNG is written"
check_contains "INFO. Skipped -- no complete BSSE-corrected data" log_default.txt
if [ -f cohesive_correction.png ]; then
    echo -e "   -> ${RED}Failed:${NC} 'cohesive_correction.png' created (should be skipped, no BSSE data)"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} 'cohesive_correction.png' absent (no BSSE data)"
    PASS=$((PASS+1))
fi

echo "Verifying the SCF-convergence diagnostic: this minimal synthetic calc.out has no"
echo "  'SCF cycle converged' line at all, so convergence can't be confirmed (advisory only)"
check_contains "WARNING. Could not confirm SCF convergence for the full structure" log_default.txt


# --- 2b. Numerical-quality diagnostics: SCF convergence confirmed + a large
#     residual force flagged. Synthetic calc.out mirroring the exact format
#     stb-cohesiveAnalysis's parsers look for (verified against a real
#     SIESTA run's calc.out during manual review of this feature). ---
echo -e "\n--- Testing SCF/force diagnostics (synthetic, real-format calc.out) ---"
mkdir -p quality_run/structure quality_run/atoms/C
cp structure/structure.fdf quality_run/structure/structure.fdf
cat > quality_run/structure/calc.out << 'EOF'
SCF cycle converged after 2 iterations

siesta: Atomic forces (eV/Ang):
----------------------------------------
   Tot    0.000534   -0.000910    0.000000
----------------------------------------
   Max    0.689325
   Res    0.317001    sqrt( Sum f_i^2 / 3N )
----------------------------------------
   Max    0.689325    constrained

siesta: FreeEng =    -500.123456
EOF
cp atoms/C/calc.out quality_run/atoms/C/calc.out
stb-cohesiveAnalysis -o calc.out -d quality_run --no-intro > log_quality.txt 2>&1
check_contains "SCF convergence.*converged after 2 iteration.s." log_quality.txt
check_contains "WARNING. Residual force on the full structure .0.6893 eV/Ang. exceeds" log_quality.txt


# --- 2c. --force-tolerance is configurable ---
echo -e "\n--- Testing --force-tolerance ---"
stb-cohesiveAnalysis -o calc.out -d quality_run --force-tolerance 1.0 --no-intro > log_force_tol.txt 2>&1
check_not_contains "exceeds --force-tolerance" log_force_tol.txt
check_contains "Residual force.*0.6893 eV/Ang .within --force-tolerance 1.0" log_force_tol.txt


# --- 3. atoms_bsse/ present and complete: both cohesive energies reported ---
echo -e "\n--- Testing with atoms_bsse/ present and complete ---"
mkdir -p atoms_bsse/C
echo "siesta: FreeEng =     -50.512345" > atoms_bsse/C/calc.out
rm -f cohesive_results.dat cohesive_correction.png
stb-cohesiveAnalysis -o calc.out -d . --save-report --no-intro > log_bsse.txt 2>&1
check_exit_code $? 0
check_contains "BSSE reference.*found" log_bsse.txt
check_contains "C.*-50.012345.*-50.512345" log_bsse.txt

echo "Verifying the BSSE-corrected cohesive energy (E_bulk - E_atom_BSSE) is less negative than uncorrected"
check_contains "Uncorrected.*-450.1111.*-450.1111" log_bsse.txt
check_contains "BSSE-corrected.*-449.6111.*-449.6111" log_bsse.txt
check_contains "BSSE correction per atom: +0.5000 eV/atom" log_bsse.txt
check_contains "BSSE-corrected" cohesive_results.dat

echo "Verifying [4] CORRECTION PLOT is built and saved when BSSE data is complete"
check_contains "OK. Saved:.*cohesive_correction.png" log_bsse.txt
check_success cohesive_correction.png

echo "Testing: --view shows an interactive matplotlib plot without crashing (non-interactive backend)"
MPLBACKEND=Agg stb-cohesiveAnalysis -o calc.out -d . --view --no-intro > log_view.txt 2>&1
check_exit_code $? 0
check_not_contains "Traceback" log_view.txt
check_contains "shown interactively via --view" log_view.txt


# --- 4. atoms_bsse/ present but incomplete: falls back to uncorrected-only,
#     with a warning -- never blocks the (still valid) uncorrected result. ---
echo -e "\n--- Testing with atoms_bsse/ present but incomplete ---"
rm -f atoms_bsse/C/calc.out
rm -f cohesive_results.dat
stb-cohesiveAnalysis -o calc.out -d . --no-intro > log_bsse_incomplete.txt 2>&1
check_exit_code $? 0
check_contains "WARNING. Could not find 'calc.out' or results for the BSSE-corrected isolated atom: C" log_bsse_incomplete.txt
check_contains "Uncorrected.*-450.1111.*-450.1111" log_bsse_incomplete.txt
check_not_contains "BSSE-corrected |" log_bsse_incomplete.txt
rm -rf atoms_bsse


# --- 4b. Multi-site BSSE layout: atoms_bsse/<sym>/site_<wyckoff>_x<mult>/ is
#     read as a multiplicity-weighted average, not just the first site found
#     (verified by hand: (6*-151.0 + 1*-152.0 + 3*-150.5)/10 = -150.95). ---
echo -e "\n--- Testing multi-site BSSE weighted average ---"
mkdir -p atoms_bsse/C/site_k_x6 atoms_bsse/C/site_b_x1 atoms_bsse/C/site_g_x3
echo "siesta: FreeEng =     -151.000000" > atoms_bsse/C/site_k_x6/calc.out
echo "siesta: FreeEng =     -152.000000" > atoms_bsse/C/site_b_x1/calc.out
echo "siesta: FreeEng =     -150.500000" > atoms_bsse/C/site_g_x3/calc.out
stb-cohesiveAnalysis -o calc.out -d . --no-intro > log_multisite.txt 2>&1
check_exit_code $? 0
check_contains "C.*-50.012345.*-150.950000" log_multisite.txt
check_contains "BSSE-corrected.*-349.1735.*-349.1735" log_multisite.txt

echo "Testing: one incomplete site drops the whole weighted average (never a partial/wrong average)"
rm -f atoms_bsse/C/site_b_x1/calc.out
stb-cohesiveAnalysis -o calc.out -d . --no-intro > log_multisite_incomplete.txt 2>&1
check_exit_code $? 0
check_contains "WARNING. Could not find 'calc.out' or results for the BSSE-corrected isolated atom: C" log_multisite_incomplete.txt
check_not_contains "BSSE-corrected |" log_multisite_incomplete.txt


# --- 4c. atoms_bsse_check/ present: reports a convergence-check delta on
#     top of the (single-site, for simplicity) BSSE-corrected result. ---
echo -e "\n--- Testing atoms_bsse_check/ convergence-check reporting ---"
rm -rf atoms_bsse
mkdir -p atoms_bsse/C atoms_bsse_check/C
echo "siesta: FreeEng =     -50.512345" > atoms_bsse/C/calc.out
echo "siesta: FreeEng =     -50.612345" > atoms_bsse_check/C/calc.out
stb-cohesiveAnalysis -o calc.out -d . --no-intro > log_conv_check.txt 2>&1
check_exit_code $? 0
check_contains "BSSE convergence check.*found" log_conv_check.txt
check_contains "BSSE check .larger cutoff..*-449.5111.*-449.5111" log_conv_check.txt
check_contains "BSSE cutoff convergence shift: .0.1000 eV/atom" log_conv_check.txt
rm -rf atoms_bsse atoms_bsse_check


# --- 4d. Multi-cutoff BSSE convergence SCAN (atoms_bsse_check_<cutoff>/, 2+
#     points -- from stb-cohesive's --bsse-convergence-increment given
#     multiple values) -- distinct from the single atoms_bsse_check/ above.
#     Base cutoff (4.0 Ang) recovered from atoms_bsse/C/bsse_cutoff.txt. ---
echo -e "\n--- Testing multi-cutoff BSSE convergence scan (atoms_bsse_check_<cutoff>/) ---"
rm -rf atoms_bsse atoms_bsse_check
mkdir -p atoms_bsse/C atoms_bsse_check_6.0/C atoms_bsse_check_8.0/C
echo "siesta: FreeEng =     -50.512345" > atoms_bsse/C/calc.out
echo "4.0" > atoms_bsse/C/bsse_cutoff.txt
echo "siesta: FreeEng =     -50.612345" > atoms_bsse_check_6.0/C/calc.out
echo "siesta: FreeEng =     -50.650000" > atoms_bsse_check_8.0/C/calc.out
MPLBACKEND=Agg stb-cohesiveAnalysis -o calc.out -d . --no-intro > log_scan.txt 2>&1
check_exit_code $? 0
check_contains "BSSE cutoff scan        . found, 2 point.s. .6.0, 8.0 Ang." log_scan.txt
check_contains "\[3b\] BSSE CUTOFF CONVERGENCE SCAN" log_scan.txt
check_contains "4.0 Ang .base.*-449.6111.*reference" log_scan.txt
check_contains "6.0 Ang.*-449.5111.*.0.1000" log_scan.txt
check_contains "8.0 Ang.*-449.4735.*.0.1377" log_scan.txt
check_contains "Last step change: .0.0377 eV/atom .cutoff 6.0 -> 8.0 Ang." log_scan.txt
check_success cohesive_correction.png

echo "Testing: an incomplete scan point (missing calc.out) is skipped with a warning, others unaffected"
rm -f atoms_bsse_check_8.0/C/calc.out
stb-cohesiveAnalysis -o calc.out -d . --no-intro > log_scan_incomplete.txt 2>&1
check_exit_code $? 0
check_contains "WARNING. Could not find 'calc.out' or results for the BSSE convergence-scan reference at cutoff 8.0 Ang for: C" log_scan_incomplete.txt
check_contains "6.0 Ang.*-449.5111" log_scan_incomplete.txt
check_not_contains "8.0 Ang.*-449.4735" log_scan_incomplete.txt

echo "Testing: base cutoff unrecoverable (no bsse_cutoff.txt, e.g. an older run) degrades gracefully"
rm -f atoms_bsse/C/bsse_cutoff.txt
echo "siesta: FreeEng =     -50.650000" > atoms_bsse_check_8.0/C/calc.out
stb-cohesiveAnalysis -o calc.out -d . --no-intro > log_scan_nocutoff.txt 2>&1
check_exit_code $? 0
check_contains "INFO. Base --bsse-cutoff value not recoverable" log_scan_nocutoff.txt
check_not_contains "Ang .base." log_scan_nocutoff.txt
rm -rf atoms_bsse atoms_bsse_check_6.0 atoms_bsse_check_8.0 cohesive_correction.png


# --- 5. Phantom species: a species declared in ChemicalSpeciesLabel but with
#     0 real atoms alongside a real, complete species must NOT block the
#     analysis -- distinct from the "totally empty structure" case below
#     (a real bug found by code review: this used to require a calc.out for
#     a species that never appears in the actual structure). ---
echo -e "\n--- Testing a real structure with an extra declared-but-unused species ---"
mkdir -p phantom_species/structure phantom_species/atoms/C
cat > phantom_species/structure/structure.fdf << 'EOF'
%block ChemicalSpeciesLabel
 1   6   C
 2   14  Si
%endblock ChemicalSpeciesLabel

LatticeConstant 1.0 Ang

AtomicCoordinatesFormat Fractional

%block LatticeVectors
 20.0 0.0 0.0
 0.0 20.0 0.0
 0.0 0.0 20.0
%endblock LatticeVectors

%block AtomicCoordinatesAndAtomicSpecies
 0.5 0.5 0.5 1
%endblock AtomicCoordinatesAndAtomicSpecies
EOF
cp structure/calc.out phantom_species/structure/
cp atoms/C/calc.out phantom_species/atoms/C/
stb-cohesiveAnalysis -o calc.out -d phantom_species --no-intro > log_phantom.txt 2>&1
check_exit_code $? 0
check_contains "Declared but never placed.*Si" log_phantom.txt
check_contains "Uncorrected.*-450.1111.*-450.1111" log_phantom.txt


# --- 6. Zero-atom structure.fdf: used to crash with an uncaught
#     ZeroDivisionError (raw traceback) instead of a clean [ERROR] message --
#     regression test for the shared-parser migration. This is a DIFFERENT
#     scenario from #5 above: here NO species has any real atom at all,
#     which read_fdf() itself rejects as malformed input. ---
echo -e "\n--- Testing a structure.fdf with a declared species but zero atoms ---"
mkdir -p zero_atoms/structure zero_atoms/atoms/C
cat > zero_atoms/structure/structure.fdf << 'EOF'
%block ChemicalSpeciesLabel
 1   6   C
%endblock ChemicalSpeciesLabel

LatticeConstant 1.0 Ang

AtomicCoordinatesFormat Fractional

%block LatticeVectors
 20.0 0.0 0.0
 0.0 20.0 0.0
 0.0 0.0 20.0
%endblock LatticeVectors

%block AtomicCoordinatesAndAtomicSpecies
%endblock AtomicCoordinatesAndAtomicSpecies
EOF
cp structure/calc.out zero_atoms/structure/
cp atoms/C/calc.out zero_atoms/atoms/C/
stb-cohesiveAnalysis -o calc.out -d zero_atoms --no-intro > log_zero_atoms.txt 2>&1
check_exit_code $? 1
check_contains "not found (or empty)" log_zero_atoms.txt
check_not_contains "Traceback" log_zero_atoms.txt


# --- 7. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: -d points at a directory that doesn't exist at all"
stb-cohesiveAnalysis -o calc.out -d /does/not/exist --no-intro > log_missing_dirs.txt 2>&1
check_exit_code $? 1
check_contains "FAIL. Directory '/does/not/exist' not found" log_missing_dirs.txt

echo "Testing: -d points at a real directory missing 'structure'/'atoms' subfolders"
mkdir -p empty_run
stb-cohesiveAnalysis -o calc.out -d empty_run --no-intro > log_missing_subdirs.txt 2>&1
check_exit_code $? 1
check_contains "Required directories 'structure' and/or 'atoms' not found" log_missing_subdirs.txt

echo "Testing: missing calc.out for the isolated atom"
mkdir -p missing_atom_out/structure missing_atom_out/atoms/C
cp structure/structure.fdf structure/calc.out missing_atom_out/structure/
stb-cohesiveAnalysis -o calc.out -d missing_atom_out --no-intro > log_missing_atom_out.txt 2>&1
check_exit_code $? 1
check_contains "Could not find 'calc.out' or results for isolated atom: C" log_missing_atom_out.txt
check_contains "Cannot calculate cohesive energy because some calculations are missing" log_missing_atom_out.txt

echo "Testing: missing required -o"
stb-cohesiveAnalysis --no-intro > log_missing_arg.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-cohesiveAnalysis --version > log_version.txt 2>&1
check_contains "stb-cohesiveAnalysis" log_version.txt

echo "Testing: --help documents -o/-d/--force-tolerance"
stb-cohesiveAnalysis --help > log_help.txt 2>&1
check_contains "\-\-out" log_help.txt
check_contains "\-\-dir" log_help.txt
check_contains "force-tolerance" log_help.txt


# --- 8. Interactive path (stb-suite, shortcut 4.3.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.3.2) ---"

echo "Testing: navigate 4.3.2 -> calc.out -> dir '.' -> default force-tolerance -> save report (y) -> view (n)"
rm -f cohesive_results.dat
# Prompts in order: out_file, dir_path ('.', since the CLI's own new default
# 'cohesive_runs' doesn't exist in this flat fixture), force-tolerance (blank
# -> default 0.05), save report (y), view (n), then the "Press Enter to
# continue" pause, then quit.
printf '4.3.2\ncalc.out\n.\n\ny\nn\n\n0\n' | MPLBACKEND=Agg stb-suite > log_menu.txt 2>&1
check_contains "CONFIGURATION SUMMARY" log_menu.txt
check_contains "Force tolerance.*0.05 eV/Ang" log_menu.txt
check_contains "Save report.*yes" log_menu.txt
check_contains "Uncorrected.*-450.1111.*-450.1111" log_menu.txt
check_success cohesive_results.dat


popd > /dev/null

# --- 9. Summary ---
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
