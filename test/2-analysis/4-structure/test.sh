#!/bin/bash

# --- Setup ---
# Smoke test for stb-structural (Structure Analyzer / ECN, item 2.4)
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

# Checks that file $1 has exactly $2 lines
check_line_count() {
    local actual
    actual=$(wc -l < "$1")
    if [ "$actual" -eq "$2" ]; then
        echo -e "   -> ${GREEN}Verified:${NC} '$1' has $actual lines (expected $2)"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} '$1' has $actual lines (expected $2)"
        FAIL=$((FAIL+1))
    fi
}


# --- 1. Preparation ---
echo "--- Starting tester for STB-Structural (item 2.4) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/siesta.STRUCT_OUT" "$TEST_DIR/"
cp "$FIXTURE_DIR/structure.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. --format fdf, mode mean (14-atom Sn3O4 structure) ---
echo -e "\n--- Testing --format fdf --mode mean ---"
rm -f structural_information.dat warnings.log
stb-structural --file structure.fdf --format fdf --mode mean --no-intro > log_fdf_mean.txt 2>&1
check_exit_code $? 0
check_success structural_information.dat

echo "Verifying lattice parameters, cell volume, and density"
check_contains "a = 4.883 Å" structural_information.dat
check_contains "b = 5.907 Å" structural_information.dat
check_contains "c = 8.238 Å" structural_information.dat
check_contains "Volume = 236.821 Å³   Density = 5.892 g/cm³" structural_information.dat

echo "Verifying the report has a title block and generation timestamp"
check_contains "STRUCTURAL PROPERTIES REPORT - STB Suite" structural_information.dat
check_contains "Generated        :" structural_information.dat
check_contains "Source file      : structure.fdf  (format: fdf)" structural_information.dat

echo "Verifying effective (weighted) CN is broken down per species, not lumped into one"
check_contains "EFFECTIVE COORDINATION NUMBER (weighted), PER SPECIES" structural_information.dat
check_contains "O (8 atoms):" structural_information.dat
check_contains "Sn (6 atoms):" structural_information.dat
check_contains "  CrystalNN      :   2.956" structural_information.dat
check_contains "  CrystalNN      :   4.004" structural_information.dat

echo "Verifying bond distance is broken down per species pair, not a single lumped average"
check_contains "AVERAGE BOND DISTANCE, PER SPECIES PAIR" structural_information.dat
check_contains "O-Sn      : 2.1109 Å  (n=48)" structural_information.dat

echo "Verifying bond angles are broken down per ligand-center-ligand triplet"
check_contains "AVERAGE BOND ANGLE, PER LIGAND-CENTER-LIGAND TRIPLET" structural_information.dat
check_contains "O-Sn-O      : 102.210°  (n=42)" structural_information.dat
check_contains "Sn-O-Sn     : 119.808°  (n=24)" structural_information.dat

echo "Verifying coordination polyhedron distortion (BLD, BAV, volume) per species"
check_contains "COORDINATION POLYHEDRON DISTORTION" structural_information.dat
check_contains "  BLD   :   0.790 %" structural_information.dat
check_contains "  BAV   : 239.172 deg²" structural_information.dat
check_contains "  Volume:  12.053 Å³" structural_information.dat

echo "Verifying a 3-coordinate site (O) correctly reports polyhedron volume as N/A (needs >=4 neighbors for a 3D hull)"
check_contains "  Volume:     N/A Å³" structural_information.dat

echo "Verifying same-species minimum distance (whole structure, both species aligned to the same field width)"
check_contains "SAME-SPECIES MINIMUM DISTANCE" structural_information.dat
check_contains "  O-O       : 2.6218 Å" structural_information.dat
check_contains "  Sn-Sn     : 3.2799 Å" structural_information.dat

echo "Verifying coordination polyhedron connectivity (corner/edge/face-sharing) per species"
check_contains "COORDINATION POLYHEDRON CONNECTIVITY" structural_information.dat
check_contains "corner-sharing: 16 pair(s)" structural_information.dat
check_contains "edge-sharing  : 1 pair(s)" structural_information.dat
check_contains "corner-sharing: 6 pair(s)" structural_information.dat
check_contains "edge-sharing  : 2 pair(s)" structural_information.dat

echo "Verifying the radial distribution function g(r) is written and integrates to the known O-Sn coordination number"
check_success rdf.dat
check_contains "RADIAL DISTRIBUTION FUNCTION g(r)" structural_information.dat
check_contains "Full curve written to ./rdf.dat" structural_information.dat
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
check_contains "ATOMIC POSITIONS (Cartesian, Å)" structural_information.dat
check_contains "   1  Sn       4.574068      3.168518      4.142994" structural_information.dat
check_not_contains "\[4.57406764" structural_information.dat


# --- 3. --format struct_out, mode mean -- same physical structure (post-
#     relaxation), so lattice/CN/bond-distance summaries should match the
#     .fdf run to the precision printed (atomic positions differ in the
#     last digit or two, which isn't checked here). ---
echo -e "\n--- Testing --format struct_out --mode mean (same structure, different source format) ---"
rm -f structural_information.dat warnings.log
stb-structural --file siesta.STRUCT_OUT --format struct_out --mode mean --no-intro > log_struct_mean.txt 2>&1
check_exit_code $? 0
check_contains "a = 4.883 Å" structural_information.dat
check_contains "  CrystalNN      :   2.956" structural_information.dat
check_contains "O-Sn      : 2.1109 Å  (n=48)" structural_information.dat


# --- 4. --mode list: per-atom effective CN for specific 1-based indices ---
echo -e "\n--- Testing --mode list ---"
rm -f structural_information.dat warnings.log
stb-structural --file structure.fdf --format fdf --mode list --list 1,7 --no-intro > log_list.txt 2>&1
check_exit_code $? 0
check_contains "EFFECTIVE COORDINATION NUMBER (weighted), SPECIFIED ATOMS" structural_information.dat
check_contains "Atom 1 (Sn), position: 4.574068  3.168518  4.142994" structural_information.dat
check_contains "Atom 7 (O), position:" structural_information.dat

echo "Verifying CrystalNN no longer errors under use_weights=True (needs weighted_cn=True at construction)"
check_not_contains "CrystalNN failed" log_list.txt

echo "Verifying per-atom CN values are formatted to 3 decimals, not printed as raw ~15-digit floats"
check_contains "  JmolNN         :   8.970" structural_information.dat
check_not_contains "8.970272100975317" structural_information.dat

echo "Verifying per-atom distortion (BAV independently hand-verified for atom 1: variance of its 15 pairwise angles = 1429.05 deg²)"
check_contains "Atom 1 (Sn):" structural_information.dat
check_contains "  BAV   : 1429.050 deg²" structural_information.dat


# --- 5. --output-dir: writes into (and creates) a chosen directory ---
echo -e "\n--- Testing --output-dir ---"
rm -rf out_dir
stb-structural --file structure.fdf --format fdf --mode mean --output-dir out_dir --no-intro > log_outdir.txt 2>&1
check_exit_code $? 0
check_success out_dir/structural_information.dat
check_success out_dir/warnings.log


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
check_not_contains "RADIAL DISTRIBUTION FUNCTION" structural_information.dat

echo "Testing: --rdf-rmax changes the cutoff used"
stb-structural --file structure.fdf --format fdf --mode mean --rdf-rmax 5 --no-intro > log_rmax.txt 2>&1
check_exit_code $? 0
check_contains "r_max = 5.0 Å" structural_information.dat

echo "Testing: --rdf-rmax 0 is rejected"
stb-structural --file structure.fdf --format fdf --mode mean --rdf-rmax 0 --no-intro > log_rmax_zero.txt 2>&1
check_exit_code $? 2


# --- 6. warnings.log doesn't grow across repeated runs (logging.basicConfig
#     defaults to append mode, which used to let it accumulate stale
#     warnings from unrelated previous runs forever) ---
echo -e "\n--- Testing warnings.log is reset each run, not appended to ---"
rm -f warnings.log
stb-structural --file structure.fdf --format fdf --mode mean --no-intro > /dev/null 2>&1
lines_first=$(wc -l < warnings.log)
stb-structural --file structure.fdf --format fdf --mode mean --no-intro > /dev/null 2>&1
lines_second=$(wc -l < warnings.log)
if [ "$lines_first" -eq "$lines_second" ]; then
    echo -e "   -> ${GREEN}Verified:${NC} warnings.log stayed at $lines_first lines across two runs (not appended)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} warnings.log grew from $lines_first to $lines_second lines across two runs"
    FAIL=$((FAIL+1))
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


# --- 8. Interactive path (stb-suite, shortcut 2.4) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 2.4) ---"

echo "Testing: navigate 2.4 -> structure.fdf -> fdf -> mean -> default output dir"
rm -f structural_information.dat warnings.log
printf '2.4\nstructure.fdf\nfdf\nmean\n\n' | stb-suite > log_menu.txt 2>&1
check_success structural_information.dat


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
