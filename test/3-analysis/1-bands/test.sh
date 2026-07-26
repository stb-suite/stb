#!/bin/bash

# --- Setup ---
# Smoke test for stb-bands (Bands Analyzer, item 3.1)
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$FIXTURE_DIR/test_files"

# Non-interactive matplotlib backend: main() calls plt.show(), which would
# otherwise block waiting for a window to close.
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

# Checks that file $1 does NOT exist
check_absent() {
    if [ ! -e "$1" ]; then
        echo -e " ... ${GREEN}OK${NC} (file '$1' correctly not created)"
        PASS=$((PASS+1))
    else
        echo -e " ... ${RED}FAIL${NC} (file '$1' should not have been created)"
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

# Checks that pattern $1 occurs exactly $3 times (grep -c) in file $2
check_occurrences() {
    local pattern="$1" file="$2" expected="$3" actual
    actual=$(grep -c -- "$pattern" "$file" 2>/dev/null)
    if [ "$actual" -eq "$expected" ]; then
        echo -e "   -> ${GREEN}Verified:${NC} '$pattern' occurs $actual time(s) in '$file' (expected $expected)"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} '$pattern' occurs $actual time(s) in '$file' (expected $expected)"
        FAIL=$((FAIL+1))
    fi
}

# Counts the non-blank lines in a whitespace-separated .dat block belonging
# to a single k-point (bands_gnuplot.dat has one row per band, blank line
# between k-points) and compares against an expected count.
check_column_count() {
    local file="$1" expected="$2" line
    line=$(head -n 1 "$file")
    local actual
    actual=$(echo "$line" | wc -w)
    if [ "$actual" -eq "$expected" ]; then
        echo -e "   -> ${GREEN}Verified:${NC} '$file' first row has $actual columns (expected $expected)"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} '$file' first row has $actual columns (expected $expected)"
        FAIL=$((FAIL+1))
    fi
}


# --- 1. Preparation ---
echo "--- Starting tester for STB-Bands (item 3.1) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/siesta.bands" "$TEST_DIR/"
cp "$FIXTURE_DIR/siesta_edge_nbands21.bands" "$TEST_DIR/"
cp "$FIXTURE_DIR/siesta_spin2.bands" "$TEST_DIR/"
cp "$FIXTURE_DIR/Sn3O4.EIG" "$TEST_DIR/"
cp "$FIXTURE_DIR/siesta_spin2_mesh.EIG" "$TEST_DIR/"
cp "$FIXTURE_DIR/siesta_bad_fermi.bands" "$TEST_DIR/"
cp "$FIXTURE_DIR/mesh_kxkyz.EIG" "$TEST_DIR/"
cp "$FIXTURE_DIR/mesh_kxkyz.KP" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Real fixture: --shift fermi, default (no report file, no --save-report) ---
echo -e "\n--- Testing --shift fermi (real fixture, 182 bands x 501 k-points, nspin=1) ---"
rm -f bands_gnuplot.dat bands.gplot bands_analysis.txt stb_bands_report.txt references.bib
stb-bands --file siesta.bands --shift fermi --no-intro > log_fermi.txt 2>&1
check_exit_code $? 0
check_success bands_gnuplot.dat
check_success bands.gplot
check_success references.bib

echo "Testing: no report file is written without --save-report (neither the old nor the new name)"
check_absent bands_analysis.txt
check_absent stb_bands_report.txt

echo "Verifying references.bib has the SIESTA citation"
check_contains "@article{Soler2002," references.bib

echo "Verifying high-symmetry labels are not quoted in the gnuplot script"
check_contains "set xtics" bands.gplot
check_not_contains "'GAMMA'" bands.gplot
check_not_contains "''" bands.gplot
check_contains '{/Symbol G}' bands.gplot

echo "Verifying bands_gnuplot.dat row shape (k, energy) pairs, one band per block"
check_column_count bands_gnuplot.dat 2


# --- 2b. --save-report: numbered [0]...[7] sections, stb_bands_report.txt ---
echo -e "\n--- Testing --save-report (numbered report, stb_bands_report.txt) ---"
rm -f bands_gnuplot.dat bands.gplot stb_bands_report.txt references.bib
stb-bands --file siesta.bands --shift fermi --save-report --no-intro > log_save_report.txt 2>&1
check_exit_code $? 0
check_success stb_bands_report.txt
check_absent bands_analysis.txt

echo "Verifying stb_bands_report.txt has the numbered sections"
check_contains "\[0\] RUN METADATA" stb_bands_report.txt
check_contains "\[1\] INPUT DATA" stb_bands_report.txt
check_contains "\[2\] BAND GAP ANALYSIS" stb_bands_report.txt
check_contains "\[4\] MESH VS LINE COMPARISON" stb_bands_report.txt
check_contains "\[5\] WRITING OUTPUT FILES" stb_bands_report.txt
check_contains "\[6\] REFERENCES" stb_bands_report.txt
check_contains "\[7\] SUMMARY & FILES" stb_bands_report.txt

echo "Verifying the physics content itself (VBM/CBM/gap) is still there"
check_contains "Fermi energy" stb_bands_report.txt
check_contains "VBM" stb_bands_report.txt
check_contains "CBM" stb_bands_report.txt
check_contains "Indirect gap" stb_bands_report.txt
check_contains "Direct gap" stb_bands_report.txt
check_contains "Gap type" stb_bands_report.txt
check_contains "Indirect" stb_bands_report.txt


# --- 3. Real fixture: other --shift modes ---
echo -e "\n--- Testing --shift vbm/cbm/manual ---"
stb-bands --file siesta.bands --shift vbm --no-intro > log_vbm.txt 2>&1
check_exit_code $? 0
stb-bands --file siesta.bands --shift cbm --no-intro > log_cbm.txt 2>&1
check_exit_code $? 0
stb-bands --file siesta.bands --shift manual --manual-value 0.5 --no-intro > log_manual.txt 2>&1
check_exit_code $? 0
stb-bands --file siesta.bands --shift manual --no-intro > log_manual_missing.txt 2>&1
check_exit_code $? 2


# --- 4. Edge-case fixture: nbands mod 10 == 1 (old parser dropped a value here) ---
echo -e "\n--- Testing continuation-line edge case (21 bands, trailing 1-value line) ---"
rm -f bands_gnuplot.dat bands.gplot stb_bands_report.txt references.bib
stb-bands --file siesta_edge_nbands21.bands --shift fermi --save-report --no-intro > log_edge.txt 2>&1
check_exit_code $? 0
check_success bands_gnuplot.dat

echo "Verifying row shape (k, energy) pairs"
check_column_count bands_gnuplot.dat 2

echo "Verifying no energy value was dropped: 21 bands x (4 k-points + 1 blank separator) = 105 lines"
lines_total=$(wc -l < bands_gnuplot.dat)
if [ "$lines_total" -eq 105 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} bands_gnuplot.dat has 105 lines (would be 100 if one value/k-point were dropped)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} bands_gnuplot.dat has $lines_total lines (expected 105)"
    FAIL=$((FAIL+1))
fi

echo "Verifying repeated Gamma labels (start and end of path) both render clean"
check_contains "Gap type" stb_bands_report.txt
check_contains '{/Symbol G}' bands.gplot
check_not_contains "'GAMMA'" bands.gplot


# --- 4b. Synthetic nspin=2 fixture: half-metallic (spin-up metallic, spin-down 2 eV gap) ---
echo -e "\n--- Testing spin-polarized bands (nspin=2, half-metallic synthetic fixture) ---"
rm -f bands_gnuplot.dat bands.gplot stb_bands_report.txt references.bib
stb-bands --file siesta_spin2.bands --shift fermi --save-report --no-intro > log_spin2.txt 2>&1
check_exit_code $? 0
check_success bands_gnuplot.dat
check_success bands.gplot
check_success stb_bands_report.txt

echo "Verifying bands_gnuplot.dat has 3 columns (k, spin-up, spin-down)"
check_column_count bands_gnuplot.dat 3

echo "Verifying bands.gplot has one plot command per spin channel"
check_contains 'using 1:2 with lines ls 1 title "Spin Up"' bands.gplot
check_contains 'using 1:3 with lines ls 2 title "Spin Down"' bands.gplot

echo "Verifying stb_bands_report.txt reports per-spin channel + half-metallic character"
check_contains "\[3\] SPIN-POLARIZED ANALYSIS" stb_bands_report.txt
check_contains "Spin up:" stb_bands_report.txt
check_contains "Spin down:" stb_bands_report.txt
check_contains "Half-metallic : Yes" stb_bands_report.txt
check_contains "Half-metallic character detected" stb_bands_report.txt

echo "Testing: --gap-tol tightened below the synthetic up-channel gap (no longer half-metallic)"
stb-bands --file siesta_spin2.bands --shift fermi --gap-tol 0.0001 --save-report --no-intro > log_spin2_tol.txt 2>&1
check_exit_code $? 0
check_contains "Half-metallic : No" stb_bands_report.txt


# --- 4c. --eig-file: mesh (k-grid) vs line (k-path) gap comparison ---
# Sn3O4.EIG is the same SIESTA calculation as siesta.bands (same Fermi
# energy, same nbands/nspin) but samples the full SCF k-mesh (89 points)
# instead of just the high-symmetry path.
echo -e "\n--- Testing --eig-file (mesh vs line gap comparison) ---"
rm -f bands_gnuplot.dat bands.gplot stb_bands_report.txt references.bib
stb-bands --file siesta.bands --eig-file Sn3O4.EIG --shift fermi --save-report --no-intro > log_eig.txt 2>&1
check_exit_code $? 0
check_success stb_bands_report.txt

echo "Verifying stb_bands_report.txt has the [2]/[4]/[7] sections with direct + indirect gaps"
check_contains "\[2\] BAND GAP ANALYSIS" stb_bands_report.txt
check_contains "\[4\] MESH VS LINE COMPARISON" stb_bands_report.txt
check_contains "\[7\] SUMMARY & FILES" stb_bands_report.txt
check_contains "VBM          | -4.471500 eV (k = 0.744142)" stb_bands_report.txt
check_contains "CBM          | -2.333400 eV (k = 1.343953)" stb_bands_report.txt
check_contains "Indirect gap | 2.138100 eV" stb_bands_report.txt
check_contains "Direct gap   | 2.444800 eV (same-k minimum, k = 1.370772)" stb_bands_report.txt
check_contains "VBM          | -4.474676 eV (k = index 42)" stb_bands_report.txt
check_contains "CBM          | -2.372184 eV (k = index 74)" stb_bands_report.txt
check_contains "Indirect gap | 2.102492 eV" stb_bands_report.txt
check_contains "Direct gap   | 2.448347 eV (same-k minimum, k = index 6)" stb_bands_report.txt
check_contains "Mesh gap is smaller than the line gap" stb_bands_report.txt
check_contains "Best fundamental (indirect) gap: 2.102492 eV (Indirect), from the k-mesh" stb_bands_report.txt

echo "Testing: --file without --eig-file still shows a MESH section, explicitly marked unavailable"
rm -f stb_bands_report.txt
stb-bands --file siesta.bands --shift fermi --save-report --no-intro > log_no_eig.txt 2>&1
check_contains "\[4\] MESH VS LINE COMPARISON" stb_bands_report.txt
check_contains "Not available -- no k-mesh eigenvalues were given" stb_bands_report.txt
check_not_contains "Indirect gap diff" stb_bands_report.txt

echo "Testing: --eig-file with a spin-polarized mesh reports the mesh per-spin section too"
rm -f bands_gnuplot.dat bands.gplot stb_bands_report.txt references.bib
stb-bands --file siesta_spin2.bands --eig-file siesta_spin2_mesh.EIG --shift fermi --save-report --no-intro > log_eig_spin.txt 2>&1
check_exit_code $? 0
check_occurrences "Spin up:" stb_bands_report.txt 2
check_occurrences "Spin down:" stb_bands_report.txt 2
check_contains "Half-metallic : Yes" stb_bands_report.txt


# --- 4c-bis. --label: auto-detect <label>.bands + <label>.EIG from a shared label ---
echo -e "\n--- Testing --label (auto-detects .bands + .EIG together) ---"
rm -f bands_gnuplot.dat bands.gplot stb_bands_report.txt references.bib
cp Sn3O4.EIG siesta.EIG
stb-bands --label siesta --shift fermi --save-report --no-intro > log_label.txt 2>&1
check_exit_code $? 0
check_success stb_bands_report.txt

echo "Verifying --label auto-detected siesta.EIG for the mesh comparison"
check_contains "Indirect gap diff" stb_bands_report.txt
rm -f siesta.EIG

echo "Testing: --label without a matching .EIG still computes the bands, mesh gap just unavailable"
rm -f bands_gnuplot.dat bands.gplot stb_bands_report.txt references.bib
stb-bands --label siesta_edge_nbands21 --shift fermi --save-report --no-intro > log_label_no_eig.txt 2>&1
check_exit_code $? 0
check_success stb_bands_report.txt
check_contains "Not available -- no k-mesh eigenvalues were given" stb_bands_report.txt
check_not_contains "Indirect gap diff" stb_bands_report.txt
check_contains "No 'siesta_edge_nbands21.EIG' found" log_label_no_eig.txt

echo "Testing: --label combined with --file/--eig-file is rejected"
stb-bands --label siesta --file siesta.bands --shift fermi --no-intro > log_label_conflict.txt 2>&1
check_exit_code $? 2


# --- 4c-ter. --kp-file: mesh VBM/CBM/direct-gap k-points reported as (kx, ky, kz) ---
# mesh_kxkyz.{EIG,KP} is a small synthetic 1-band-pair, 3-k-point mesh with
# known-by-hand extrema: VBM/direct gap at k-index 1 (0.1,0.2,0.3 1/Ang),
# CBM at k-index 2 (0.4,0.1,0.0 1/Ang); indirect gap 0.9 eV, direct gap 1.1 eV.
echo -e "\n--- Testing --kp-file (mesh k-points reported as kx,ky,kz instead of a bare index) ---"
rm -f bands_gnuplot.dat bands.gplot stb_bands_report.txt references.bib
stb-bands --file siesta.bands --eig-file mesh_kxkyz.EIG --kp-file mesh_kxkyz.KP --shift fermi --save-report --no-intro > log_kp.txt 2>&1
check_exit_code $? 0
check_success stb_bands_report.txt

echo "Verifying mesh VBM/CBM/direct-gap k-points are reported as (kx, ky, kz)"
check_contains "Mesh k-points  : 3" stb_bands_report.txt
check_contains "VBM          | -0.500000 eV (k = kx=0.100000, ky=0.200000, kz=0.300000 1/Ang (index 1))" stb_bands_report.txt
check_contains "CBM          | 0.400000 eV (k = kx=0.400000, ky=0.100000, kz=0.000000 1/Ang (index 2))" stb_bands_report.txt
check_contains "Indirect gap | 0.900000 eV" stb_bands_report.txt
check_contains "Direct gap   | 1.100000 eV (same-k minimum, k = kx=0.100000, ky=0.200000, kz=0.300000 1/Ang (index 1))" stb_bands_report.txt

echo "Testing: --kp-file without --eig-file is rejected"
stb-bands --file siesta.bands --kp-file mesh_kxkyz.KP --shift fermi --no-intro > log_kp_only.txt 2>&1
check_exit_code $? 2


# --- 4d. Robustness: Fermi energy outside the eigenvalue range ---
echo -e "\n--- Testing Fermi energy outside the eigenvalue range (should error, not report 'inf') ---"
stb-bands --file siesta_bad_fermi.bands --shift fermi --no-intro > log_bad_fermi.txt 2>&1
check_exit_code $? 1
check_contains "No occupied/empty states found" log_bad_fermi.txt


# --- 5. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: missing --file (required)"
stb-bands --shift fermi --no-intro > log_missing_file_arg.txt 2>&1
check_exit_code $? 2

echo "Testing: missing --shift (required)"
stb-bands --file siesta.bands --no-intro > log_missing_shift.txt 2>&1
check_exit_code $? 2

echo "Testing: nonexistent input file"
stb-bands --file does_not_exist.bands --shift fermi --no-intro > log_missing_input.txt 2>&1
check_exit_code $? 1

echo "Testing: --version"
stb-bands --version > log_version.txt 2>&1
check_contains "stb-bands" log_version.txt

echo "Testing: --help documents --save-report"
stb-bands --help > log_help.txt 2>&1
check_contains "save-report" log_help.txt
check_contains "stb_bands_report.txt" log_help.txt


# --- 6. Interactive path (stb-suite, shortcut 3.1) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 3.1) ---"

echo "Testing: navigate 3.1 -> label 'siesta' -> shift by Fermi level -> output dir '.' -> no save-report -> quit"
rm -f bands_gnuplot.dat bands.gplot stb_bands_report.txt references.bib
printf '3.1\nsiesta\n3\n.\nn\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_success bands_gnuplot.dat
check_absent stb_bands_report.txt

echo "Testing: navigate 3.1 -> label 'siesta' -> shift by Fermi level -> output dir '.' -> save-report yes -> quit"
rm -f bands_gnuplot.dat bands.gplot stb_bands_report.txt references.bib
printf '3.1\nsiesta\n3\n.\ny\n\n0\n' | stb-suite > log_menu_report.txt 2>&1
check_success bands_gnuplot.dat
check_success stb_bands_report.txt


popd > /dev/null

# --- 7. Summary ---
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
