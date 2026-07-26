#!/bin/bash

# --- Setup ---
# Smoke test for stb-structural (Structure Analyzer / ECN, item 3.4)
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
echo "--- Starting tester for STB-Structural (item 3.4) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/siesta.STRUCT_OUT" "$TEST_DIR/"
cp "$FIXTURE_DIR/structure.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. --format fdf, mode mean (14-atom Sn3O4 structure) ---
echo -e "\n--- Testing --format fdf --mode mean ---"
rm -f stb_structural_report.txt structural_information.dat warnings.log rdf.dat references.bib
stb-structural --file structure.fdf --format fdf --mode mean --no-intro > log_fdf_mean.txt 2>&1
check_exit_code $? 0

echo "Verifying the old structural_information.dat is no longer written unconditionally, and no report file without --save-report"
if [ -e structural_information.dat ] || [ -e stb_structural_report.txt ]; then
    echo -e "   -> ${RED}Failed:${NC} a report file exists without --save-report"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no structural_information.dat, no stb_structural_report.txt"
    PASS=$((PASS+1))
fi
check_success references.bib

echo "Verifying no warnings.log is written, and CrystalNN's expected, non-actionable 'no oxidation states' warnings don't print either"
if [ -e warnings.log ]; then
    echo -e "   -> ${RED}Failed:${NC} warnings.log was written (should not be -- warnings print to console only)"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} warnings.log absent"
    PASS=$((PASS+1))
fi
check_not_contains "\[WARNING\]" log_fdf_mean.txt

echo "Verifying lattice parameters, cell volume, and density"
check_contains "^a        | 4.883 Å" log_fdf_mean.txt
check_contains "^b        | 5.907 Å" log_fdf_mean.txt
check_contains "^c        | 8.238 Å" log_fdf_mean.txt
check_contains "^Volume   | 236.821 Å³" log_fdf_mean.txt
check_contains "^Density  | 5.892 g/cm³" log_fdf_mean.txt

echo "Verifying the report has the numbered [0]-[11] section headers"
check_contains "\[0\] RUN METADATA" log_fdf_mean.txt
check_contains "\[1\] LATTICE" log_fdf_mean.txt
check_contains "\[2\] EFFECTIVE COORDINATION NUMBER" log_fdf_mean.txt
check_contains "\[3\] BOND DISTANCES" log_fdf_mean.txt
check_contains "\[4\] BOND ANGLES" log_fdf_mean.txt
check_contains "\[5\] COORDINATION POLYHEDRON DISTORTION" log_fdf_mean.txt
check_contains "\[6\] SAME-SPECIES MINIMUM DISTANCE" log_fdf_mean.txt
check_contains "\[7\] COORDINATION POLYHEDRON CONNECTIVITY" log_fdf_mean.txt
check_contains "\[8\] RADIAL DISTRIBUTION FUNCTION" log_fdf_mean.txt
check_contains "\[9\] ATOMIC POSITIONS" log_fdf_mean.txt
check_contains "\[10\] REFERENCES" log_fdf_mean.txt
check_contains "\[11\] SUMMARY & FILES" log_fdf_mean.txt
check_contains "Input file     : structure.fdf (format: fdf)" log_fdf_mean.txt
check_contains "Analysis mode  : mean" log_fdf_mean.txt

echo "Verifying effective (weighted) CN is broken down per species, not lumped into one"
check_contains "^O (8 atoms):" log_fdf_mean.txt
check_contains "^Sn (6 atoms):" log_fdf_mean.txt
check_contains "CrystalNN | 2.956" log_fdf_mean.txt
check_contains "CrystalNN | 4.004" log_fdf_mean.txt

echo "Verifying bond distance is broken down per species pair, not a single lumped average"
check_contains "O-Sn    | 2.1109           | 48" log_fdf_mean.txt

echo "Verifying bond angles are broken down per ligand-center-ligand triplet"
check_contains "O-Sn-O  | 102.210       | 42" log_fdf_mean.txt
check_contains "Sn-O-Sn | 119.808       | 24" log_fdf_mean.txt

echo "Verifying coordination polyhedron distortion (BLD, BAV, volume) per species"
check_contains "O       | 8       | 0.790   | 239.172    | N/A" log_fdf_mean.txt
check_contains "Sn      | 6       | 0.695   | 537.811    | 12.053" log_fdf_mean.txt

echo "Verifying same-species minimum distance (whole structure)"
check_contains "O-O   | 2.6218" log_fdf_mean.txt
check_contains "Sn-Sn | 3.2799" log_fdf_mean.txt

echo "Verifying coordination polyhedron connectivity (corner/edge/face-sharing) per species"
check_contains "O-O          | 16             | 1            | 0" log_fdf_mean.txt
check_contains "Sn-Sn        | 6              | 2            | 0" log_fdf_mean.txt

echo "Verifying the radial distribution function g(r) is written and integrates to the known O-Sn coordination number"
check_success rdf.dat
check_contains "Full curve written to ./rdf.dat" log_fdf_mean.txt
check_contains "g_O-Sn(r)" rdf.dat
python3 -c "
import numpy as np
data = np.loadtxt('rdf.dat')
r, g_osn = data[:,0], data[:,3]
dr = r[1]-r[0]
shell = 4*np.pi*r**2*dr
rho_Sn = 6/236.82072939038554
mask = r < 2.6
n = np.sum(g_osn[mask]*shell[mask]*rho_Sn)
assert abs(n - 3.0) < 0.01, f'integrated O-Sn coordination {n} != 3.0'
print('OK')
" > rdf_check.txt 2>&1
check_contains "OK" rdf_check.txt

echo "Verifying atomic positions are cleanly formatted (fixed-width columns, not numpy's raw array repr)"
check_contains "^1    | Sn      | 4.574068  | 3.168518 | 4.142994" log_fdf_mean.txt
check_not_contains "\[4.57406764" log_fdf_mean.txt

echo "Verifying references.bib has the SIESTA citation"
check_contains "@article{Soler2002," references.bib


# --- 2b. --save-report persists the exact same report to disk ---
echo -e "\n--- Testing --save-report ---"
rm -f stb_structural_report.txt
stb-structural --file structure.fdf --format fdf --mode mean --save-report --no-intro > log_saved.txt 2>&1
check_exit_code $? 0
check_success stb_structural_report.txt
check_contains "\[0\] RUN METADATA" stb_structural_report.txt
check_contains "\[11\] SUMMARY & FILES" stb_structural_report.txt
check_contains "Report         : ./stb_structural_report.txt" stb_structural_report.txt
check_contains "CrystalNN | 2.956" stb_structural_report.txt


# --- 3. --format struct_out, mode mean -- same physical structure (post-
#     relaxation), so lattice/CN/bond-distance summaries should match the
#     .fdf run to the precision printed (atomic positions differ in the
#     last digit or two, which isn't checked here). ---
echo -e "\n--- Testing --format struct_out --mode mean (same structure, different source format) ---"
stb-structural --file siesta.STRUCT_OUT --format struct_out --mode mean --no-intro > log_struct_mean.txt 2>&1
check_exit_code $? 0
check_contains "^a        | 4.883 Å" log_struct_mean.txt
check_contains "CrystalNN | 2.956" log_struct_mean.txt
check_contains "O-Sn    | 2.1109           | 48" log_struct_mean.txt


# --- 4. --mode list: per-atom effective CN for specific 1-based indices ---
echo -e "\n--- Testing --mode list ---"
stb-structural --file structure.fdf --format fdf --mode list --list 1,7 --no-intro > log_list.txt 2>&1
check_exit_code $? 0
check_contains "\[2\] EFFECTIVE COORDINATION NUMBER" log_list.txt
check_contains "Analysis mode  : list (atoms: 1, 7)" log_list.txt
check_contains "Atom 1 (Sn), position: 4.574068  3.168518  4.142994" log_list.txt
check_contains "Atom 7 (O), position:" log_list.txt

echo "Verifying CrystalNN no longer errors under use_weights=True (needs weighted_cn=True at construction)"
check_not_contains "CrystalNN failed" log_list.txt

echo "Verifying per-atom CN values are formatted to 3 decimals, not printed as raw ~15-digit floats"
check_contains "MinDistNN | 5.936" log_list.txt
check_not_contains "5.935876979917421" log_list.txt

echo "Verifying JmolNN was dropped (consistently the worst-agreement outlier of the five methods)"
check_not_contains "JmolNN" log_list.txt

echo "Verifying per-atom distortion (BAV independently hand-verified for atom 1: variance of its 15 pairwise angles = 1429.05 deg²)"
check_contains "1    | Sn      | 0.910   | 1429.050" log_list.txt


# --- 5. --output-dir: writes into (and creates) a chosen directory ---
echo -e "\n--- Testing --output-dir ---"
rm -rf out_dir
stb-structural --file structure.fdf --format fdf --mode mean --output-dir out_dir --save-report --no-intro > log_outdir.txt 2>&1
check_exit_code $? 0
check_success out_dir/stb_structural_report.txt
check_success out_dir/references.bib
if [ -e out_dir/warnings.log ]; then
    echo -e "   -> ${RED}Failed:${NC} out_dir/warnings.log was written (should not be)"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} out_dir/warnings.log absent"
    PASS=$((PASS+1))
fi


# --- 5b. --no-rdf: skips rdf.dat, and cleans up a stale one from a
#     previous run in the same directory ---
echo -e "\n--- Testing --no-rdf ---"
rm -f rdf.dat
stb-structural --file structure.fdf --format fdf --mode mean --no-intro > /dev/null 2>&1
check_success rdf.dat
stb-structural --file structure.fdf --format fdf --mode mean --no-rdf --no-intro > log_nordf.txt 2>&1
check_exit_code $? 0
if [ -e rdf.dat ]; then
    echo -e "   -> ${RED}Failed:${NC} rdf.dat still exists after --no-rdf (stale file not cleaned up)"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} rdf.dat absent after --no-rdf (and stale one removed)"
    PASS=$((PASS+1))
fi
check_not_contains "RADIAL DISTRIBUTION FUNCTION" log_nordf.txt

echo "Testing: --rdf-rmax changes the cutoff used"
stb-structural --file structure.fdf --format fdf --mode mean --rdf-rmax 5 --no-intro > log_rmax.txt 2>&1
check_exit_code $? 0
check_contains "r_max = 5.0 Å" log_rmax.txt

echo "Testing: --rdf-rmax 0 is rejected"
stb-structural --file structure.fdf --format fdf --mode mean --rdf-rmax 0 --no-intro > log_rmax_zero.txt 2>&1
check_exit_code $? 2


# --- 5c. A stale structural_information.dat left by an older version of
#     this tool (which always wrote one) is cleaned up rather than left
#     to linger and be mistaken for current output. ---
echo -e "\n--- Testing a stale structural_information.dat from a previous (older) version is removed ---"
echo "stale content from an older run" > structural_information.dat
stb-structural --file structure.fdf --format fdf --mode mean --no-intro > /dev/null 2>&1
if [ -e structural_information.dat ]; then
    echo -e "   -> ${RED}Failed:${NC} stale structural_information.dat was not cleaned up"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} stale structural_information.dat was removed"
    PASS=$((PASS+1))
fi


# --- 6. A stale warnings.log from an older version of this tool (which
#     always wrote one) is cleaned up rather than left to linger. ---
echo -e "\n--- Testing a stale warnings.log from a previous run is removed ---"
echo "stale content from an older run" > warnings.log
stb-structural --file structure.fdf --format fdf --mode mean --no-intro > /dev/null 2>&1
if [ -e warnings.log ]; then
    echo -e "   -> ${RED}Failed:${NC} stale warnings.log was not cleaned up"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} stale warnings.log was removed"
    PASS=$((PASS+1))
fi


# --- 7. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: atom index out of range"
stb-structural --file structure.fdf --format fdf --mode list --list 1,99 --no-intro > log_oob.txt 2>&1
check_exit_code $? 1
check_contains "out of range" log_oob.txt

echo "Testing: nonexistent input file"
stb-structural --file does_not_exist.fdf --format fdf --mode mean --no-intro > log_missing.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing.txt

echo "Testing: --mode list without --list"
stb-structural --file structure.fdf --format fdf --mode list --no-intro > log_missing_list.txt 2>&1
check_exit_code $? 2

echo "Testing: invalid --format (cif/poscar no longer accepted)"
stb-structural --file structure.fdf --format cif --mode mean --no-intro > log_badfmt.txt 2>&1
check_exit_code $? 2

echo "Testing: missing required arguments"
stb-structural --file structure.fdf --no-intro > log_missing_args.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-structural --version > log_version.txt 2>&1
check_contains "stb-structural" log_version.txt


# --- 8. Interactive path (stb-suite, shortcut 3.4) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 3.4) ---"

echo "Testing: navigate 3.4 -> structure.fdf -> format=1(fdf) -> mode=1(mean) -> default output dir -> default RDF -> save-report=N (default)"
rm -f stb_structural_report.txt structural_information.dat warnings.log rdf.dat
printf '3.4\nstructure.fdf\n1\n1\n\n\n\n\n' | stb-suite > log_menu.txt 2>&1
check_success rdf.dat
if [ -e stb_structural_report.txt ] || [ -e structural_information.dat ]; then
    echo -e "   -> ${RED}Failed:${NC} a report file was created despite the default (N) save-report answer"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no report file created (default save-report answer is N)"
    PASS=$((PASS+1))
fi
check_contains "Select input file format:" log_menu.txt
check_contains "Select analysis mode:" log_menu.txt
check_contains "Also save a text report to file?" log_menu.txt

echo "Testing: navigate 3.4 -> siesta.STRUCT_OUT -> format=2(struct_out) -> mode=2(list) -> atoms 1,7 -> RDF=n -> save-report=y"
rm -f stb_structural_report.txt structural_information.dat warnings.log rdf.dat
printf '3.4\nsiesta.STRUCT_OUT\n2\n2\n\n1,7\nn\ny\n' | stb-suite > log_menu_list.txt 2>&1
check_success stb_structural_report.txt
if [ -e rdf.dat ]; then
    echo -e "   -> ${RED}Failed:${NC} rdf.dat was created despite answering 'n' to the RDF prompt"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} rdf.dat absent (RDF prompt answered 'n')"
    PASS=$((PASS+1))
fi
check_contains "\[2\] EFFECTIVE COORDINATION NUMBER" stb_structural_report.txt
check_contains "atoms: 1, 7" stb_structural_report.txt


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
