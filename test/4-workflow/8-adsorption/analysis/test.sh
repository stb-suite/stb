#!/bin/bash

# --- Setup ---
# Smoke test for stb-adsorbAnalysis (Adsorption Analysis, item 4.8.2). Uses
# stb-adsorb itself (real tool, no SIESTA needed for prep) to build the
# clean_slab/adsorbate/sites folder layout against the same fixture as
# ../prep/, then fabricates a "siesta: FreeEng" line per folder (synthetic,
# same style already used by 3-cohesive/5-convergence/4-phonons) to exercise
# the analysis side without needing a real SIESTA run.
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREP_DIR="$(cd "$FIXTURE_DIR/../prep" && pwd)"
TEST_DIR="$FIXTURE_DIR/test_files"

# Output colors
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

PASS=0
FAIL=0

# --- Check helpers ---

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
    if grep -q "$1" "$2" 2>/dev/null; then
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
echo "--- Starting tester for STB-AdsorbAnalysis (item 4.8.2) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$PREP_DIR/structure.fdf" "$TEST_DIR/"
cp "$PREP_DIR/calc.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Build 3 candidate sites via stb-adsorb, then fabricate FreeEng ---
echo -e "\n--- Testing analysis of a 3-site adsorption study ---"
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate H --site-type all --all-sites --no-intro \
    > log_prep.txt 2>&1
n_sites=$(find sites -maxdepth 1 -type d -name 'site_*' | wc -l)
if [ "$n_sites" -ge 3 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} prep produced >= 3 candidate sites ($n_sites)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} expected >= 3 candidate sites, found $n_sites"
    FAIL=$((FAIL+1))
fi

printf 'siesta: FreeEng =    -200.000000\nSCF cycle converged after 12 iterations\nsiesta: Atomic forces (eV/Ang):\n   Max    0.010000\n' > clean_slab/calc.out
printf 'siesta: FreeEng =     -13.500000\nSCF cycle converged after 9 iterations\nsiesta: Atomic forces (eV/Ang):\n   Max    0.005000\n' > adsorbate/calc.out
site_dirs=(sites/site_*/)
printf 'siesta: FreeEng =    -214.100000\nSCF cycle converged after 14 iterations\nsiesta: Atomic forces (eV/Ang):\n   Max    0.020000\n' > "${site_dirs[0]}calc.out"
printf 'siesta: FreeEng =    -214.350000\nSCF cycle converged after 15 iterations\nsiesta: Atomic forces (eV/Ang):\n   Max    0.015000\n' > "${site_dirs[1]}calc.out"
# site 2 is deliberately left unconverged AND with a large residual force, to
# verify report_quality_diagnostics/the per-site table flag it instead of
# silently trusting E_ads
printf 'siesta: FreeEng =    -213.900000\nsiesta: Atomic forces (eV/Ang):\n   Max    0.350000\n' > "${site_dirs[2]}calc.out"
most_stable=$(basename "${site_dirs[1]}")
worst_quality=$(basename "${site_dirs[2]}")

# bsse_slab/bsse_adsorbate folders already exist (stb-adsorb's --bsse-
# correction defaults ON) but have no calc.out yet -- every site's BSSE
# correction is "incomplete", so this run should report it as entirely
# unavailable, not partially available.
stb-adsorbAnalysis --dir . --no-intro > log_analysis.txt 2>&1
check_contains "E_clean_slab : -200.000000 eV" log_analysis.txt
check_contains "E_adsorbate (H) : -13.500000 eV" log_analysis.txt
check_contains "\-0.850000" log_analysis.txt
check_contains "Most stable site (uncorrected):.*$most_stable" log_analysis.txt
check_contains "exothermic" log_analysis.txt
check_contains "No BSSE-corrected results found" log_analysis.txt
check_success adsorption_curve.dat
check_success adsorption_curve.gplot
check_success adsorption_report.txt
check_contains "Most stable site (uncorrected):.*$most_stable" adsorption_report.txt
# SCF-convergence / residual-force diagnostics: clean_slab and adsorbate both
# converged cleanly above, so no warning should be emitted for either
if grep -q "Could not confirm SCF convergence for clean_slab" log_analysis.txt; then
    echo -e "   -> ${RED}Failed:${NC} unexpected SCF warning for clean_slab (fixture has a converged line)"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no spurious SCF warning for clean_slab"
    PASS=$((PASS+1))
fi
check_contains "never confirmed SCF convergence.*$worst_quality" log_analysis.txt
check_contains "0.3500" log_analysis.txt
check_contains "residual force above --force-tolerance" log_analysis.txt


# --- 2b. BSSE-corrected analysis: complete for all 3 sites ---
echo -e "\n--- Testing BSSE-corrected analysis (complete for all sites) ---"
# Deliberately make site 0's BSSE-corrected energy the most stable
# (E_ads_bsse = -0.700), different from the uncorrected ranking (site 1,
# E_ads_bsse = -0.300) -- demonstrates the correction can actually change
# the answer, not just shift every site uniformly.
echo "siesta: FreeEng =    -200.020000" > "${site_dirs[0]}bsse_slab/calc.out"
echo "siesta: FreeEng =     -13.380000" > "${site_dirs[0]}bsse_adsorbate/calc.out"
echo "siesta: FreeEng =    -200.500000" > "${site_dirs[1]}bsse_slab/calc.out"
echo "siesta: FreeEng =     -13.550000" > "${site_dirs[1]}bsse_adsorbate/calc.out"
echo "siesta: FreeEng =    -200.200000" > "${site_dirs[2]}bsse_slab/calc.out"
echo "siesta: FreeEng =     -13.600000" > "${site_dirs[2]}bsse_adsorbate/calc.out"
bsse_most_stable=$(basename "${site_dirs[0]}")

stb-adsorbAnalysis --dir . --no-intro > log_bsse.txt 2>&1
check_contains "E_ads_BSSE" log_bsse.txt
check_contains "Most stable site (uncorrected):.*$most_stable" log_bsse.txt
check_contains "Most stable site (BSSE-corrected):.*$bsse_most_stable" log_bsse.txt
check_contains "BSSE correction at that site:" log_bsse.txt
check_contains "5:Label" adsorption_curve.dat
check_contains "(BSSE)" adsorption_curve.gplot


# --- 2c. BSSE-corrected analysis: incomplete for one site ---
echo -e "\n--- Testing BSSE-corrected analysis (incomplete for one site) ---"
rm -f "${site_dirs[2]}bsse_adsorbate/calc.out"
stb-adsorbAnalysis --dir . --no-intro > log_bsse_partial.txt 2>&1
check_contains "incomplete/unreadable results" log_bsse_partial.txt
check_contains "available for only 2/3 site" log_bsse_partial.txt
echo "siesta: FreeEng =     -13.600000" > "${site_dirs[2]}bsse_adsorbate/calc.out"


# --- 2d. Multi-adsorbate + height-sweep round trip: real stb-adsorb prep,
#     fabricated energies, stb-adsorbAnalysis reads them back correctly ---
echo -e "\n--- Testing multi-adsorbate + height-sweep round trip ---"
mkdir -p multi
cp structure.fdf calc.fdf multi/
(
    cd multi
    stb-adsorb -s structure.fdf -c calc.fdf --adsorbate O,N --site-type ontop \
        --height-sweep 1.5 2.5 0.5 --no-bsse-correction --no-intro > log_multi_prep.txt 2>&1
    echo "siesta: FreeEng =    -200.000000" > clean_slab/calc.out
    echo "siesta: FreeEng =      -6.500000" > adsorbate_O/calc.out
    echo "siesta: FreeEng =      -8.200000" > adsorbate_N/calc.out
    # O approach curve: weakest at short/long, minimum around h=2.0
    echo "siesta: FreeEng =    -206.300000" > sites/site_1_ontop_O_h1.50/calc.out
    echo "siesta: FreeEng =    -206.900000" > sites/site_1_ontop_O_h2.00/calc.out
    echo "siesta: FreeEng =    -206.600000" > sites/site_1_ontop_O_h2.50/calc.out
    # N approach curve: most stable at the shortest height tested
    echo "siesta: FreeEng =    -209.100000" > sites/site_1_ontop_N_h1.50/calc.out
    echo "siesta: FreeEng =    -208.700000" > sites/site_1_ontop_N_h2.00/calc.out
    echo "siesta: FreeEng =    -208.300000" > sites/site_1_ontop_N_h2.50/calc.out
    stb-adsorbAnalysis --dir . --no-intro > log_multi_analysis.txt 2>&1
)
check_contains "E_adsorbate (N) : -8.200000 eV" multi/log_multi_analysis.txt
check_contains "E_adsorbate (O) : -6.500000 eV" multi/log_multi_analysis.txt
check_contains "Most stable site (uncorrected):.*site_1_ontop_N_h1.50" multi/log_multi_analysis.txt
check_contains "Best site per adsorbate" multi/log_multi_analysis.txt
check_contains "N: site_1_ontop_N_h1.50" multi/log_multi_analysis.txt
check_contains "O: site_1_ontop_O_h2.00" multi/log_multi_analysis.txt
check_success multi/height_curve_site_1_ontop_N.dat
check_success multi/height_curve_site_1_ontop_O.dat
check_contains "1.5000  -0.900000" multi/height_curve_site_1_ontop_N.dat


# --- 2e. --apply copies the most stable site's structure.fdf ---
echo -e "\n--- Testing --apply ---"
(
    cd multi
    stb-adsorbAnalysis --dir . --apply best_prod.fdf --no-intro > log_apply.txt 2>&1
)
check_contains "Applied.*site_1_ontop_N_h1.50" multi/log_apply.txt
check_success multi/best_prod.fdf
check_contains "NumberofAtoms      3" multi/best_prod.fdf


# --- 3. A folder missing calc.out is skipped, not fatal ---
echo -e "\n--- Testing that a site missing calc.out is skipped ---"
rm -f "${site_dirs[2]}calc.out"
stb-adsorbAnalysis --dir . --no-intro > log_partial.txt 2>&1
check_contains "SKIP" log_partial.txt
check_contains "skipped: 2" log_partial.txt
echo "siesta: FreeEng =    -213.900000" > "${site_dirs[2]}calc.out"


# --- 4. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: missing 'sites' directory entirely"
mkdir -p no_sites_dir
stb-adsorbAnalysis --dir no_sites_dir --no-intro > log_no_sites.txt 2>&1
check_exit_code $? 1
check_contains "Did you run stb-adsorb" log_no_sites.txt

echo "Testing: missing clean_slab reference energy"
mkdir -p missing_ref/sites/site_1_ontop
echo "siesta: FreeEng =    -1.0" > missing_ref/sites/site_1_ontop/calc.out
stb-adsorbAnalysis --dir missing_ref --no-intro > log_missing_ref.txt 2>&1
check_exit_code $? 1
check_contains "Could not read energy" log_missing_ref.txt

echo "Testing: --version"
stb-adsorbAnalysis --version > log_version.txt 2>&1
check_contains "stb-adsorbAnalysis" log_version.txt

echo "Testing: --help documents --dir/--file"
stb-adsorbAnalysis --help > log_help.txt 2>&1
check_contains "dir" log_help.txt
check_contains "file" log_help.txt


# --- 5. Interactive path (stb-suite, shortcut 4.8.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.8.2) ---"

echo "Testing: navigate 4.8.2 -> defaults -> quit"
rm -f adsorption_curve.dat adsorption_report.txt
# 4.8.2 (menu code) / . (dir) / "" (out_file default) / "" (force-tolerance
# default) / "" (apply_target: skip) / "" (Press Enter to continue) / 0 (quit)
printf '4.8.2\n.\n\n\n\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_contains "Most stable site (uncorrected):.*$most_stable" log_menu.txt
check_contains "Most stable site (BSSE-corrected):.*$bsse_most_stable" log_menu.txt
check_contains "Force tolerance" log_menu.txt
check_success adsorption_curve.dat


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
