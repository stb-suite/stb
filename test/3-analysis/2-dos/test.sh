#!/bin/bash

# --- Setup ---
# Smoke test for stb-dos (PDOS XML Parser, item 3.2)
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$FIXTURE_DIR/test_files"

# stb-dos now imports matplotlib for --view; Agg keeps that a no-op
# instead of blocking on a GUI window, same convention test/3-analysis/1-bands/test.sh already uses.
export MPLBACKEND=Agg

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
echo "--- Starting tester for STB-DOS (item 3.2) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/siesta.PDOS.xml" "$TEST_DIR/"
cp "$FIXTURE_DIR/spin2.PDOS.xml" "$TEST_DIR/"
cp "$FIXTURE_DIR/lmax4.PDOS.xml" "$TEST_DIR/"
cp "$FIXTURE_DIR/nspin_missing.PDOS.xml" "$TEST_DIR/"
cp "$FIXTURE_DIR/no_energy_values.PDOS.xml" "$TEST_DIR/"
cp "$FIXTURE_DIR/no_orbitals.PDOS.xml" "$TEST_DIR/"
cp "$FIXTURE_DIR/malformed.PDOS.xml" "$TEST_DIR/"
cp "$FIXTURE_DIR/siesta.bands" "$TEST_DIR/"
cp "$FIXTURE_DIR/siesta.EIG" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Real fixture: default run (nspin=1, 182 orbitals, 14 atoms) ---
echo -e "\n--- Testing default run (real fixture, nspin=1, 182 orbitals) ---"
rm -rf dos_total.dat dos_per_atom dos_per_species
stb-dos siesta.PDOS.xml --no-intro > log_default.txt 2>&1
check_exit_code $? 0
check_success dos_total.dat
check_success dos_per_atom/Sn_1.dat
check_success dos_per_species/dos_Sn.dat
check_success dos_per_species/dos_O.dat

echo "Verifying dos_total.dat columns (nspin=1: no spin suffix)"
check_contains "#Energy(eV)" dos_total.dat
check_contains "s" dos_total.dat
check_contains "p" dos_total.dat
check_contains "d" dos_total.dat
check_not_contains "_up" dos_total.dat
check_not_contains "_down" dos_total.dat


# --- 3. --type filtering ---
echo -e "\n--- Testing --type filtering ---"
rm -rf dos_total.dat dos_per_atom dos_per_species
stb-dos siesta.PDOS.xml --type total --no-intro > log_type_total.txt 2>&1
check_exit_code $? 0
check_success dos_total.dat
if [ -d dos_per_atom ] || [ -d dos_per_species ]; then
    echo -e "   -> ${RED}Failed:${NC} --type total also created dos_per_atom/dos_per_species"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} --type total did not create dos_per_atom/dos_per_species"
    PASS=$((PASS+1))
fi


# --- 4. --shift modes ---
echo -e "\n--- Testing --shift fermi/manual/invalid ---"
stb-dos siesta.PDOS.xml --shift fermi --no-intro > log_shift_fermi.txt 2>&1
check_exit_code $? 0
check_contains "Using automatic Fermi energy shift" log_shift_fermi.txt

stb-dos siesta.PDOS.xml --shift 0.0 --no-intro > log_shift_zero.txt 2>&1
check_exit_code $? 0
check_contains "Using manual energy shift: 0.0 eV" log_shift_zero.txt

stb-dos siesta.PDOS.xml --shift notanumber --no-intro > log_shift_bad.txt 2>&1
check_exit_code $? 1
check_contains "Invalid shift value" log_shift_bad.txt


# --- 4a. --shift vbm/cbm hierarchy: .bands -> .EIG -> DOS estimate.
# siesta.bands/siesta.EIG are the SAME real SIESTA calculation as
# siesta.PDOS.xml (matching Fermi energy), so the VBM/CBM reported here
# must match exactly what stb-bands itself reports for these files. ---
echo -e "\n--- Testing --shift vbm/cbm hierarchy (.bands -> .EIG -> DOS estimate) ---"

echo "Testing: both siesta.bands and siesta.EIG present -> .bands wins"
stb-dos siesta.PDOS.xml --shift vbm --no-intro > log_shift_vbm_bands.txt 2>&1
check_exit_code $? 0
check_contains "Using VBM shift from './siesta.bands' (k-path .bands): -4.471500 eV" log_shift_vbm_bands.txt

stb-dos siesta.PDOS.xml --shift cbm --no-intro > log_shift_cbm_bands.txt 2>&1
check_exit_code $? 0
check_contains "Using CBM shift from './siesta.bands' (k-path .bands): -2.333400 eV" log_shift_cbm_bands.txt

echo "Testing: only siesta.EIG present (no .bands) -> falls back to .EIG"
mv siesta.bands siesta.bands.bak
stb-dos siesta.PDOS.xml --shift vbm --no-intro > log_shift_vbm_eig.txt 2>&1
check_exit_code $? 0
check_contains "Using VBM shift from './siesta.EIG' (k-mesh .EIG): -4.474676 eV" log_shift_vbm_eig.txt

echo "Testing: neither .bands nor .EIG present, no --estimate-from-dos -> clear error"
mv siesta.EIG siesta.EIG.bak
stb-dos siesta.PDOS.xml --shift vbm --no-intro > log_shift_vbm_none.txt 2>&1
check_exit_code $? 1
check_contains "no 'siesta.bands' or 'siesta.EIG' found" log_shift_vbm_none.txt

echo "Testing: neither present, with --estimate-from-dos -> approximation + warning"
stb-dos siesta.PDOS.xml --shift vbm --estimate-from-dos --no-intro > log_shift_vbm_estimate.txt 2>&1
check_exit_code $? 0
check_contains "estimating VBM directly from the DOS" log_shift_vbm_estimate.txt
check_contains "This is an approximation" log_shift_vbm_estimate.txt

# Restore for later tests / cleanliness.
mv siesta.bands.bak siesta.bands
mv siesta.EIG.bak siesta.EIG

echo "Testing: --label alone (no positional filename) resolves <label>.PDOS.xml"
rm -rf dos_total.dat dos_per_atom dos_per_species
stb-dos --label siesta --no-intro > log_label_only.txt 2>&1
check_exit_code $? 0
check_success dos_total.dat

echo "Testing: neither filename nor --label given -> argparse error"
stb-dos --no-intro > log_no_filename_no_label.txt 2>&1
check_exit_code $? 2
check_contains "one of filename or --label is required" log_no_filename_no_label.txt

echo "Testing: explicit filename (no --label) still auto-derives the label for .bands/.EIG lookup"
stb-dos siesta.PDOS.xml --shift vbm --no-intro > log_derived_label.txt 2>&1
check_exit_code $? 0
check_contains "Using VBM shift from './siesta.bands' (k-path .bands): -4.471500 eV" log_derived_label.txt


# --- 4b. --output-dir: writes into (and creates) a nested directory ---
echo -e "\n--- Testing --output-dir (nested, auto-created) ---"
rm -rf my_out dos_total.dat dos_per_atom dos_per_species
stb-dos siesta.PDOS.xml --type total --output-dir my_out/nested --no-intro > log_outdir.txt 2>&1
check_exit_code $? 0
check_success my_out/nested/dos_total.dat
check_contains "Saved Total DOS to my_out/nested/dos_total.dat" log_outdir.txt
if [ -e dos_total.dat ]; then
    echo -e "   -> ${RED}Failed:${NC} dos_total.dat was also written to CWD instead of only --output-dir"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} nothing was written to CWD, only to --output-dir"
    PASS=$((PASS+1))
fi
rm -rf my_out


# --- 5. --projection ml (detailed m-resolved columns) ---
echo -e "\n--- Testing --projection ml ---"
rm -rf dos_total.dat dos_per_atom dos_per_species
stb-dos siesta.PDOS.xml --type total --projection ml --no-intro > log_ml.txt 2>&1
check_exit_code $? 0
check_contains "px" dos_total.dat
check_contains "py" dos_total.dat
check_contains "pz" dos_total.dat


# --- 6. nspin=2: spin-polarized data must survive (regression test for the
# silent-data-loss bug: <nspin> was never read, so every orbital's data
# length mismatched num_energy_points and got silently dropped -- output
# used to contain only the Energy(eV) column). ---
echo -e "\n--- Testing nspin=2 (spin-polarized PDOS data is not silently dropped) ---"
rm -rf dos_total.dat dos_per_atom dos_per_species
stb-dos spin2.PDOS.xml --shift 0.0 --no-intro > log_spin2.txt 2>&1
check_exit_code $? 0
check_contains "Detected nspin=2" log_spin2.txt

echo "Verifying dos_total.dat has spin-resolved s_up/s_down columns with the right values"
check_contains "s_up" dos_total.dat
check_contains "s_down" dos_total.dat
check_contains "-1.000000E+00	  1.000000E+00	  5.000000E-01" dos_total.dat
check_contains "0.000000E+00	  2.000000E+00	  6.000000E-01" dos_total.dat
check_contains "1.000000E+00	  3.000000E+00	  7.000000E-01" dos_total.dat


# --- 7. <nspin> tag missing: falls back to nspin=1 instead of erroring ---
echo -e "\n--- Testing missing <nspin> tag (falls back to nspin=1) ---"
rm -rf dos_total.dat dos_per_atom dos_per_species
stb-dos nspin_missing.PDOS.xml --shift 0.0 --no-intro > log_nspin_missing.txt 2>&1
check_exit_code $? 0
check_contains "<nspin> tag not found" log_nspin_missing.txt
check_success dos_total.dat
check_not_contains "_up" dos_total.dat


# --- 8. l > 3 (g-orbitals) are counted and reported, not silently dropped ---
echo -e "\n--- Testing l>3 orbitals are excluded but reported ---"
rm -rf dos_total.dat dos_per_atom dos_per_species
stb-dos lmax4.PDOS.xml --shift 0.0 --no-intro > log_lmax4.txt 2>&1
check_exit_code $? 0
check_contains "skipped 1 orbital" log_lmax4.txt
check_success dos_total.dat
check_contains "#Energy(eV)" dos_total.dat


# --- 8b. Standardized report: --save-report/--save-gnuplot/--view,
# all opt-in and off by default ---
echo -e "\n--- Testing --save-report/--save-gnuplot/--view (all off by default) ---"

echo "Testing: default run -- no report file, no .gplot files"
rm -rf dos_total.dat dos_per_atom dos_per_species stb_dos_report.txt references.bib
stb-dos siesta.PDOS.xml --shift fermi --no-intro > log_default_report.txt 2>&1
check_exit_code $? 0
check_contains "\[0\] RUN METADATA" log_default_report.txt
check_contains "\[5\] SUMMARY & FILES" log_default_report.txt
if [ -e stb_dos_report.txt ] || [ -e dos_total.gplot ]; then
    echo -e "   -> ${RED}Failed:${NC} stb_dos_report.txt/dos_total.gplot exists without the opt-in flags"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no stb_dos_report.txt, no dos_total.gplot without the opt-in flags"
    PASS=$((PASS+1))
fi

echo "Testing: --save-report writes the full numbered report to disk"
stb-dos siesta.PDOS.xml --shift fermi --save-report --no-intro > log_save_report.txt 2>&1
check_exit_code $? 0
check_success stb_dos_report.txt
check_contains "\[0\] RUN METADATA" stb_dos_report.txt
check_contains "\[1\] INPUT DATA" stb_dos_report.txt
check_contains "\[2\] ENERGY SHIFT" stb_dos_report.txt
check_contains "\[3\] WRITING OUTPUT FILES" stb_dos_report.txt
check_contains "\[4\] REFERENCES" stb_dos_report.txt
check_contains "\[5\] SUMMARY & FILES" stb_dos_report.txt
check_success references.bib

echo "Testing: --save-gnuplot writes one .gplot for 'total' (per-orbital) and one per category (one curve per atom/orbital, not summed)"
rm -rf dos_total.dat dos_per_atom dos_per_species
stb-dos siesta.PDOS.xml --shift fermi --save-gnuplot --no-intro > log_save_gnuplot.txt 2>&1
check_exit_code $? 0
check_success dos_total.gplot
check_success dos_per_atom/dos_per_atom.gplot
check_success dos_per_species/dos_per_species.gplot
if [ -e dos_per_atom/Sn_1.gplot ] || [ -e dos_per_species/dos_Sn.gplot ]; then
    echo -e "   -> ${RED}Failed:${NC} a per-file .gplot was written inside dos_per_atom/dos_per_species (should be one per category, not one per file)"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no per-file .gplot inside dos_per_atom/dos_per_species, only the one per-category script"
    PASS=$((PASS+1))
fi
check_contains 'plot for \[i=2:' dos_total.gplot
check_contains 'columnheader(i)' dos_total.gplot
check_not_contains 'sum \[i=2:' dos_per_atom/dos_per_atom.gplot
check_not_contains 'sum \[i=2:' dos_per_species/dos_per_species.gplot
check_contains 'Sn_1: s' dos_per_atom/dos_per_atom.gplot
check_contains 'Sn_1: p' dos_per_atom/dos_per_atom.gplot
check_contains 'Sn: s' dos_per_species/dos_per_species.gplot
if command -v gnuplot > /dev/null 2>&1; then
    echo "Testing: the saved .gplot files actually render with a real gnuplot"
    (gnuplot dos_total.gplot > /dev/null 2>&1)
    check_success dos_total.pdf
    (cd dos_per_atom && gnuplot dos_per_atom.gplot > /dev/null 2>&1)
    check_success dos_per_atom/dos_per_atom.pdf
    (cd dos_per_species && gnuplot dos_per_species.gplot > /dev/null 2>&1)
    check_success dos_per_species/dos_per_species.pdf
else
    echo -e "   -> ${YELLOW}Skipped:${NC} gnuplot not installed, cannot render the .gplot files"
fi

echo "Testing: --view (matplotlib preview) runs to completion under MPLBACKEND=Agg"
rm -rf dos_total.dat dos_per_atom dos_per_species
stb-dos siesta.PDOS.xml --shift fermi --view --no-intro > log_view.txt 2>&1
check_exit_code $? 0
check_success dos_total.dat


# --- 9. Error and robustness cases (all must now exit non-zero, not 0) ---
echo -e "\n--- Testing error cases ---"

echo "Testing: nonexistent input file"
stb-dos does_not_exist.xml --no-intro > log_missing_input.txt 2>&1
check_exit_code $? 1
check_contains "File not found" log_missing_input.txt

echo "Testing: malformed XML"
stb-dos malformed.PDOS.xml --no-intro > log_malformed.txt 2>&1
check_exit_code $? 1
check_contains "Error parsing XML file" log_malformed.txt

echo "Testing: no <energy_values> tag"
stb-dos no_energy_values.PDOS.xml --no-intro > log_no_energy.txt 2>&1
check_exit_code $? 1
check_contains "energy_values" log_no_energy.txt

echo "Testing: no <orbital> tags"
stb-dos no_orbitals.PDOS.xml --no-intro > log_no_orbitals.txt 2>&1
check_exit_code $? 1
check_contains "No <orbital> tags found" log_no_orbitals.txt

echo "Testing: --version"
stb-dos --version > log_version.txt 2>&1
check_contains "stb-dos" log_version.txt


# --- 10. Interactive path (stb-suite, shortcut 2.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 2.2) ---"

echo "Testing: navigate 2.2 -> siesta.PDOS.xml -> defaults (incl. save-report/save-gnuplot/view prompts)"
rm -rf dos_total.dat dos_per_atom dos_per_species
printf '3.2\nsiesta.PDOS.xml\n\n\n\n\n\n\n\n' | stb-suite > log_menu.txt 2>&1
check_success dos_total.dat

echo "Testing: interactive path surfaces a tool failure (malformed XML -> non-zero exit -> run_tool reports it)"
printf '3.2\nmalformed.PDOS.xml\n\n\n\n\n\n\n\n' | stb-suite > log_menu_fail.txt 2>&1
check_contains "Error running stb-dos" log_menu_fail.txt

echo "Testing: interactive path forwards a custom output directory"
rm -rf menu_out
printf '3.2\nsiesta.PDOS.xml\n\n\n\nmenu_out\n\n\n\n' | stb-suite > log_menu_outdir.txt 2>&1
check_success menu_out/dos_total.dat

echo "Testing: interactive path answering y/y/n saves a report and per-category gnuplot scripts, skips --view"
rm -rf dos_total.dat dos_per_atom dos_per_species stb_dos_report.txt
printf '3.2\nsiesta.PDOS.xml\n\n\n\n\ny\ny\nn\n' | stb-suite > log_menu_saveflags.txt 2>&1
check_success stb_dos_report.txt
check_success dos_total.gplot

echo "Testing: interactive path, numbered shift menu choice '4' (custom value)"
rm -rf dos_total.dat dos_per_atom dos_per_species dos_total.gplot stb_dos_report.txt
printf '3.2\nsiesta.PDOS.xml\n\n4\n-1.5\n\n\n\n\n\n' | stb-suite > log_menu_manual.txt 2>&1
check_contains "Using manual energy shift: -1.5 eV" log_menu_manual.txt

echo "Testing: interactive path, invalid shift menu choice is rejected and re-prompted"
rm -rf dos_total.dat dos_per_atom dos_per_species dos_total.gplot stb_dos_report.txt
printf '3.2\nsiesta.PDOS.xml\n\n9\n2\n\n\n\n\n\n' | stb-suite > log_menu_invalid_choice.txt 2>&1
check_contains "Invalid choice! Please select 1-4." log_menu_invalid_choice.txt
check_contains "Using VBM shift from './siesta.bands' (k-path .bands): -4.471500 eV" log_menu_invalid_choice.txt

echo "Testing: interactive path answering y to the --view prompt runs to completion (MPLBACKEND=Agg)"
rm -rf dos_total.dat dos_per_atom dos_per_species dos_total.gplot stb_dos_report.txt
printf '3.2\nsiesta.PDOS.xml\n\n\n\n\n\n\ny\n' | stb-suite > log_menu_view.txt 2>&1
check_success dos_total.dat
check_not_contains "Error running stb-dos" log_menu_view.txt

echo "Testing: interactive path with just a label, numbered shift menu choice '2' (VBM) resolved via .bands"
rm -rf dos_total.dat dos_per_atom dos_per_species dos_total.gplot stb_dos_report.txt
printf '3.2\nsiesta\n\n2\n\n\n\n\n\n\n0\n' | stb-suite > log_menu_vbm.txt 2>&1
check_contains "Will use 'siesta.bands' to compute VBM" log_menu_vbm.txt
check_contains "Using VBM shift from './siesta.bands' (k-path .bands): -4.471500 eV" log_menu_vbm.txt

echo "Testing: interactive path, numbered shift choice '2' (VBM) with no .bands/.EIG -> falls back to fermi with a warning"
mv siesta.bands siesta.bands.bak
mv siesta.EIG siesta.EIG.bak
printf '3.2\nsiesta\n\n2\nn\n\n\n\n\n\n0\n' | stb-suite > log_menu_vbm_none.txt 2>&1
check_contains "No '.bands'/'.EIG' found. Estimate VBM/CBM from the DOS itself instead?" log_menu_vbm_none.txt
check_contains "Falling back to --shift fermi" log_menu_vbm_none.txt
mv siesta.bands.bak siesta.bands
mv siesta.EIG.bak siesta.EIG


popd > /dev/null

# --- 11. Summary ---
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
