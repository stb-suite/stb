#!/bin/bash

# --- Setup ---
# Smoke test for stb-cohesive (Cohesive Energy Prep, item 4.3.1)
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
echo "--- Starting tester for STB-Cohesive prep (item 4.3.1) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR/pp"
cp "$FIXTURE_DIR/structure.fdf" "$TEST_DIR/"
# Minimal placeholder pseudopotentials -- only existence/format matters here,
# not real pseudopotential content.
echo "dummy pseudopotential content" > "$TEST_DIR/pp/C.psf"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Default run: vacuum-aware k-grid (10-atom graphene-like slab, 20 Ang
#     out-of-plane vacuum -> c-axis forced to a single division). BSSE is ON
#     by default -- disabled here (--no-bsse-correction) to isolate the
#     pre-existing k-grid/spin behavior; BSSE itself is tested in its own
#     section below. ---
echo -e "\n--- Testing default run (vacuum-aware k-grid) ---"
rm -rf structure atoms atoms_bsse
stb-cohesive -s structure.fdf --no-bsse-correction --no-intro > log_default.txt 2>&1
check_exit_code $? 0
check_success structure/calc.fdf
check_success structure/structure.fdf
check_success atoms/C/calc.fdf
check_success atoms/C/structure.fdf

echo "Verifying dimensionality is now reported"
check_contains "Detected dimensionality: 2D" log_default.txt

echo "Verifying the in-plane axes get a real k-grid but the vacuum axis is forced to 1"
check_contains "Calculated K-grid for full structure: 7 7 1" log_default.txt
check_contains "kgrid.MonkhorstPack   \[7  7  1\]" structure/calc.fdf

echo "Verifying the isolated atom is always Gamma-only"
check_contains "kgrid.MonkhorstPack   \[1  1  1\]" atoms/C/calc.fdf

echo "Verifying spin: full structure defaults to non-polarized, isolated atom is always polarized"
check_contains "Spin                non-polarized" structure/calc.fdf
check_contains "Spin                polarized" atoms/C/calc.fdf

echo "Verifying --no-bsse-correction skips atoms_bsse/ and prints the caveat note"
if [ -d atoms_bsse ]; then
    echo -e "   -> ${RED}Failed:${NC} atoms_bsse/ was created despite --no-bsse-correction"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} atoms_bsse/ not created"
    PASS=$((PASS+1))
fi
check_contains "NOTE. --bsse-correction is OFF" log_default.txt


# --- 3. --spin: only the full structure calc changes ---
echo -e "\n--- Testing --spin ---"
rm -rf structure atoms atoms_bsse
stb-cohesive -s structure.fdf --spin --no-bsse-correction --no-intro > log_spin.txt 2>&1
check_exit_code $? 0
check_contains "Spin                polarized" structure/calc.fdf
check_contains "Spin polarization ENABLED" log_spin.txt


# --- 4. --pp-path: .psf pseudopotentials are linked (regression test -- this
#     tool used to only recognize .psml and silently drop .psf, leaving the
#     generated directories without any pseudopotential at all) ---
echo -e "\n--- Testing --pp-path with a .psf pseudopotential ---"
rm -rf structure atoms atoms_bsse
stb-cohesive -s structure.fdf -p pp --no-bsse-correction --no-intro > log_psf.txt 2>&1
check_exit_code $? 0
check_success structure/C.psf
check_success atoms/C/C.psf
if grep -qi "not found" log_psf.txt; then
    echo -e "   -> ${RED}Failed:${NC} a pseudopotential 'not found' warning was printed despite pp/C.psf existing"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no 'not found' warning printed"
    PASS=$((PASS+1))
fi


# --- 5. --pp-path: .psml takes priority over .psf when both exist ---
echo -e "\n--- Testing .psml priority over .psf ---"
rm -rf structure atoms atoms_bsse
echo "dummy pseudopotential content" > pp/C.psml
stb-cohesive -s structure.fdf -p pp --no-bsse-correction --no-intro > log_psml.txt 2>&1
check_exit_code $? 0
check_success structure/C.psml
if [ -e structure/C.psf ]; then
    echo -e "   -> ${RED}Failed:${NC} C.psf was also linked even though C.psml exists"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} only C.psml was linked, not C.psf"
    PASS=$((PASS+1))
fi


# --- 6. --pp-path accepts a bundled bank name ---
echo -e "\n--- Testing --pp-path with a bundled bank name (dojo) ---"
rm -rf structure atoms atoms_bsse
stb-cohesive -s structure.fdf -p dojo --no-bsse-correction --no-intro > log_bank.txt 2>&1
check_exit_code $? 0
check_success structure/C.psml
rm -f pp/C.psml

echo "Testing: using a bundled bank prints its citation/origin"
check_contains "van Setten" log_bank.txt
check_contains "pseudo-dojo.org" log_bank.txt


# --- 7. --vacuum: isolated-atom box size is configurable ---
echo -e "\n--- Testing --vacuum ---"
rm -rf structure atoms atoms_bsse
stb-cohesive -s structure.fdf --vacuum 30 --no-bsse-correction --no-intro > log_vacuum.txt 2>&1
check_exit_code $? 0
check_contains "30.000000   0.000000   0.000000" atoms/C/structure.fdf

echo "Testing: --vacuum 0 is rejected"
stb-cohesive -s structure.fdf --vacuum 0 --no-intro > log_vacuum_zero.txt 2>&1
check_exit_code $? 2


# --- 8. --vacuum-gap: dimensionality detection threshold is configurable ---
echo -e "\n--- Testing --vacuum-gap ---"
rm -rf structure atoms atoms_bsse
echo "Testing: a --vacuum-gap larger than the fixture's 20 Ang out-of-plane gap reports 3D instead"
stb-cohesive -s structure.fdf --vacuum-gap 25 --no-bsse-correction --no-intro > log_vacuum_gap.txt 2>&1
check_exit_code $? 0
check_contains "Detected dimensionality: 3D" log_vacuum_gap.txt
if grep -q "Calculated K-grid for full structure: 7 7 [^1]" log_vacuum_gap.txt; then
    echo -e "   -> ${GREEN}Verified:${NC} c-axis is no longer forced to 1 division"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} c-axis still forced to 1 with a larger --vacuum-gap"
    FAIL=$((FAIL+1))
fi


# --- 9. -O/--output-dir: structure/atoms isolated under a chosen root ---
echo -e "\n--- Testing --output-dir ---"
rm -rf structure atoms atoms_bsse run1
stb-cohesive -s structure.fdf --output-dir run1 --no-bsse-correction --no-intro > log_outdir.txt 2>&1
check_exit_code $? 0
check_success run1/structure/calc.fdf
check_success run1/atoms/C/calc.fdf
if [ -e structure/calc.fdf ]; then
    echo -e "   -> ${RED}Failed:${NC} structure/calc.fdf was also created outside --output-dir"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} nothing created outside run1/"
    PASS=$((PASS+1))
fi
rm -rf run1


# --- 10. Phantom species: a species declared in ChemicalSpeciesLabel but
#     never placed in the coordinates block must NOT get an isolated-atom
#     (or BSSE) folder -- it would contribute exactly 0 to the cohesive-
#     energy sum regardless of its own energy, so requiring one only wastes
#     DFT time and (on the analysis side) incorrectly blocks the whole
#     result. Regression test for a real bug found by code review. ---
echo -e "\n--- Testing a declared-but-unused ('phantom') species ---"
rm -rf structure atoms atoms_bsse
cat > phantom.fdf << 'EOF'
NumberOfSpecies    2
NumberofAtoms      2

%block ChemicalSpeciesLabel
 1   6   C
 2   14  Si
%endblock ChemicalSpeciesLabel

LatticeConstant 1.0 Ang

AtomicCoordinatesFormat Fractional

%block LatticeVectors
 5.43 0.0 0.0
 0.0 5.43 0.0
 0.0 0.0 5.43
%endblock LatticeVectors

%block AtomicCoordinatesAndAtomicSpecies
 0.0 0.0 0.0 1
 0.25 0.25 0.25 1
%endblock AtomicCoordinatesAndAtomicSpecies
EOF
stb-cohesive -s phantom.fdf --no-bsse-correction --no-intro > log_phantom.txt 2>&1
check_exit_code $? 0
check_contains "Detected species: C" log_phantom.txt
check_contains "Declared but never placed.*Si" log_phantom.txt
check_success atoms/C/calc.fdf
if [ -d atoms/Si ]; then
    echo -e "   -> ${RED}Failed:${NC} atoms/Si/ was created for a species with 0 real atoms"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no atoms/Si/ folder created"
    PASS=$((PASS+1))
fi


# --- 11. --bsse-correction (the default): ghost-cluster generation ---
echo -e "\n--- Testing --bsse-correction (default ON) ---"
rm -rf structure atoms atoms_bsse
stb-cohesive -s phantom.fdf --no-intro > log_bsse.txt 2>&1
check_exit_code $? 0
check_contains "BSSE (COUNTERPOISE) CORRECTION" log_bsse.txt
check_success atoms_bsse/C/calc.fdf
check_success atoms_bsse/C/structure.fdf
if [ -d atoms_bsse/Si ]; then
    echo -e "   -> ${RED}Failed:${NC} atoms_bsse/Si/ was created for a phantom species"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no atoms_bsse/Si/ folder created"
    PASS=$((PASS+1))
fi

echo "Verifying the ghost cluster has a real anchor (Z=6) plus a ghost neighbor (Z=-6)"
check_contains "1   6   C" atoms_bsse/C/structure.fdf
check_contains "2   -6   C_ghost" atoms_bsse/C/structure.fdf
check_contains "NumberOfSpecies    2" atoms_bsse/C/structure.fdf

echo "Verifying the ghost pseudopotential is symlinked from the real element's file"
stb-cohesive -s phantom.fdf -p pp --no-intro > log_bsse_pp.txt 2>&1
check_success atoms_bsse/C/C_ghost.psf
if [ "$(readlink -f atoms_bsse/C/C_ghost.psf)" = "$(readlink -f atoms_bsse/C/C.psf)" ]; then
    echo -e "   -> ${GREEN}Verified:${NC} C_ghost.psf resolves to the same file as the real C.psf"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} C_ghost.psf does not resolve to the real C.psf"
    FAIL=$((FAIL+1))
fi


# --- 12. --bsse-cutoff: configurable, and a --vacuum too small for it warns ---
echo -e "\n--- Testing --bsse-cutoff safety check ---"
rm -rf structure atoms atoms_bsse
stb-cohesive -s phantom.fdf --vacuum 10 --bsse-cutoff 4 --no-intro > log_bsse_small_box.txt 2>&1
check_exit_code $? 0
check_contains "WARNING. --vacuum 10" log_bsse_small_box.txt


# --- 12b. Multi-site BSSE (default ON): the graphene-like fixture has 3
# symmetrically distinct C sites (multiplicities 6, 1, 3 = 10 atoms total)
# -- one ghost cluster per site, in nested atoms_bsse/<sym>/site_<x>_<n>/. ---
echo -e "\n--- Testing multi-site BSSE (default, 3 distinct sites) ---"
rm -rf structure atoms atoms_bsse
stb-cohesive -s structure.fdf --no-intro > log_multisite.txt 2>&1
check_exit_code $? 0
check_contains "3 symmetrically distinct site.s. detected for C" log_multisite.txt
check_success atoms_bsse/C/site_k_x6/structure.fdf
check_success atoms_bsse/C/site_b_x1/structure.fdf
check_success atoms_bsse/C/site_g_x3/structure.fdf
if [ -f atoms_bsse/C/structure.fdf ]; then
    echo -e "   -> ${RED}Failed:${NC} a flat atoms_bsse/C/structure.fdf was also created (should be nested-only)"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no flat atoms_bsse/C/structure.fdf (nested layout only)"
    PASS=$((PASS+1))
fi

echo "Testing: --no-bsse-multi-site collapses to a single (first) site, with a warning naming the skipped ones"
rm -rf structure atoms atoms_bsse
stb-cohesive -s structure.fdf --no-bsse-multi-site --no-intro > log_single_site.txt 2>&1
check_exit_code $? 0
check_contains "WARNING. --no-bsse-multi-site: using only the first site" log_single_site.txt
check_contains "skipping b .x1., g .x3." log_single_site.txt
check_success atoms_bsse/C/structure.fdf
if [ -d atoms_bsse/C/site_k_x6 ]; then
    echo -e "   -> ${RED}Failed:${NC} nested per-site folders were created despite --no-bsse-multi-site"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} flat single-site layout used"
    PASS=$((PASS+1))
fi


# --- 12c. --bsse-convergence-check: a second, larger-cutoff reference per
# site/species, mirroring the same layout under atoms_bsse_check/ ---
echo -e "\n--- Testing --bsse-convergence-check ---"
rm -rf structure atoms atoms_bsse atoms_bsse_check
stb-cohesive -s phantom.fdf --bsse-convergence-check --no-intro > log_conv_check.txt 2>&1
check_exit_code $? 0
check_contains "Also generated .atoms_bsse_check/. at cutoff 6.0 Ang" log_conv_check.txt
check_success atoms_bsse_check/C/structure.fdf
check_success atoms_bsse_check/C/calc.fdf

echo "Testing: --bsse-convergence-check requires --bsse-correction"
stb-cohesive -s phantom.fdf --no-bsse-correction --bsse-convergence-check --no-intro > log_conv_check_noreq.txt 2>&1
check_exit_code $? 2


# --- 12d. --dispersion: DFT-D3 enabled uniformly across every calc.fdf ---
echo -e "\n--- Testing --dispersion (DFT-D3) ---"
rm -rf structure atoms atoms_bsse
stb-cohesive -s phantom.fdf --dispersion --no-intro > log_dispersion.txt 2>&1
check_exit_code $? 0
check_contains "DFT-D3 dispersion correction ENABLED" log_dispersion.txt
check_contains "DFTD3                   .true." structure/calc.fdf
check_contains "DFTD3                   .true." atoms/C/calc.fdf
check_contains "DFTD3                   .true." atoms_bsse/C/calc.fdf

echo "Testing: --dispersion OFF by default (DFTD3 .false. everywhere)"
rm -rf structure atoms atoms_bsse
stb-cohesive -s phantom.fdf --no-intro > log_no_dispersion.txt 2>&1
check_contains "DFTD3                   .false." structure/calc.fdf
check_contains "DFTD3                   .false." atoms/C/calc.fdf


# --- 13. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: missing structure file"
stb-cohesive -s does_not_exist.fdf --no-intro > log_missing.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing.txt

echo "Testing: missing required -s"
stb-cohesive --no-intro > log_missing_arg.txt 2>&1
check_exit_code $? 2

echo "Testing: --bsse-cutoff 0 is rejected"
stb-cohesive -s structure.fdf --bsse-cutoff 0 --no-intro > log_bsse_cutoff_zero.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-cohesive --version > log_version.txt 2>&1
check_contains "stb-cohesive" log_version.txt

echo "Testing: --help documents --vacuum-gap/--output-dir/--bsse-correction/--bsse-cutoff/--bsse-multi-site/--bsse-convergence-check/--dispersion"
stb-cohesive --help > log_help.txt 2>&1
check_contains "vacuum-gap" log_help.txt
check_contains "output-dir" log_help.txt
check_contains "bsse-correction" log_help.txt
check_contains "bsse-cutoff" log_help.txt
check_contains "bsse-multi-site" log_help.txt
check_contains "bsse-convergence-check" log_help.txt
check_contains "dispersion" log_help.txt
check_contains "symprec" log_help.txt


# --- 14. Interactive path (stb-suite, shortcut 4.3.1) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.3.1) ---"

echo "Testing: navigate 4.3.1 -> structure.fdf -> defaults (dispersion N, BSSE on, multi-site on, conv-check N) -> skip advanced -> quit"
rm -rf structure atoms atoms_bsse
# Prompts in order: file, k-density (blank -> default), pp-path (blank ->
# skip), spin (n), dispersion (n), vacuum (blank -> default), BSSE choice
# (blank -> Y, default), multi-site choice (blank -> Y, default),
# convergence-check choice (n), advanced settings (n -> skip), then the
# "Press Enter to continue" pause, then quit.
printf '4.3.1\nstructure.fdf\n\n\nn\nn\n\n\n\nn\n\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_contains "CONFIGURATION SUMMARY" log_menu.txt
check_success structure/calc.fdf
check_success atoms/C/calc.fdf
# This fixture has 3 symmetrically distinct C sites (multi-site is ON by
# default), so the BSSE layout is nested rather than a flat atoms_bsse/C/.
check_success atoms_bsse/C/site_k_x6/calc.fdf

echo "Testing: navigate 4.3.1 -> decline BSSE -> no atoms_bsse/ created (multi-site/conv-check prompts skipped entirely)"
rm -rf structure atoms atoms_bsse
printf '4.3.1\nstructure.fdf\n\n\nn\nn\n\nn\n\n\n0\n' | stb-suite > log_menu_nobsse.txt 2>&1
check_contains "BSSE correction.*OFF" log_menu_nobsse.txt
check_success structure/calc.fdf
if [ -d atoms_bsse ]; then
    echo -e "   -> ${RED}Failed:${NC} atoms_bsse/ was created despite declining BSSE in the menu"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no atoms_bsse/ created"
    PASS=$((PASS+1))
fi


popd > /dev/null

# --- 15. Summary ---
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
