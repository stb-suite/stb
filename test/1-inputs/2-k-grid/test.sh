#!/bin/bash

# --- Setup ---
# Smoke test for stb-kgrid (K-Grid Generator, item 1.2)
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$FIXTURE_DIR/test_files"
STRUCT="$FIXTURE_DIR/structure.fdf"

# Output colors
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

PASS=0
FAIL=0

# --- Check helpers ---

# Checks that file $2 contains (grep -q) pattern $1
check_contains() {
    if grep -q "$1" "$2" 2>/dev/null; then
        echo -e "   -> ${GREEN}Verified:${NC} '$1' found in '$2'"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} '$1' NOT found in '$2'"
        FAIL=$((FAIL+1))
    fi
}

# Checks that file $2 does NOT contain (grep -q) pattern $1
check_not_contains() {
    if ! grep -q "$1" "$2" 2>/dev/null; then
        echo -e "   -> ${GREEN}Verified:${NC} '$1' absent from '$2' (as expected)"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} '$1' found in '$2' (should not be there)"
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
echo "--- Starting tester for STB-KGrid (item 1.2) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$STRUCT" "$TEST_DIR/structure.fdf"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null

# --- 2. Fixtures ---
# A cubic 5.43 Ang lattice, expressed in all 4 supported formats, so the
# expected grid at density=0.2 is identical (6 6 6) and independently
# verifiable by hand: 2*pi/5.43 = 1.157 1/Ang; ceil(1.157/0.2) = 6.
echo -e "\n--- Generating fixtures ---"

cat > cubic.poscar << 'EOF'
Cubic Silicon Cell
1.0
5.430 0.000 0.000
0.000 5.430 0.000
0.000 0.000 5.430
Si
2
Direct
0.00 0.00 0.00
0.25 0.25 0.25
EOF

cat > cubic.cif << 'EOF'
data_silicon
_cell_length_a     5.43
_cell_length_b     5.43
_cell_length_c     5.43
_cell_angle_alpha  90.0
_cell_angle_beta   90.0
_cell_angle_gamma  90.0
_symmetry_space_group_name_H-M 'P 1'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Si1 Si 0.00 0.00 0.00
Si2 Si 0.25 0.25 0.25
EOF

cat > cubic.fhi << 'EOF'
lattice_vector 5.430 0.000 0.000
lattice_vector 0.000 5.430 0.000
lattice_vector 0.000 0.000 5.430
atom_frac 0.00 0.00 0.00 Si
atom_frac 0.25 0.25 0.25 Si
EOF

cat > cubic.fdf << 'EOF'
SystemName          Cubic Silicon

LatticeConstant     1.0 Ang

%block LatticeVectors
5.430 0.000 0.000
0.000 5.430 0.000
0.000 0.000 5.430
%endblock LatticeVectors

AtomicCoordinatesFormat  Fractional

%block ChemicalSpeciesLabel
1  14  Si
%endblock ChemicalSpeciesLabel

%block AtomicCoordinatesAndAtomicSpecies
0.00 0.00 0.00 1
0.25 0.25 0.25 1
%endblock AtomicCoordinatesAndAtomicSpecies
EOF

# Orthogonal boxes for the dimensionality heuristic (vacuum sides use 35 Ang,
# large enough that 2*pi/35 = 0.1795 1/Ang < density=0.2, so ceil(...) = 1).
cat > vac_0d.poscar << 'EOF'
0D Molecule Box
1.0
40.000 0.000 0.000
0.000 40.000 0.000
0.000 0.000 40.000
Si
1
Direct
0.00 0.00 0.00
EOF

cat > vac_1d.poscar << 'EOF'
1D Chain Box
1.0
5.430 0.000 0.000
0.000 35.000 0.000
0.000 0.000 35.000
Si
1
Direct
0.00 0.00 0.00
EOF

cat > vac_2d.poscar << 'EOF'
2D Slab Box
1.0
5.430 0.000 0.000
0.000 5.430 0.000
0.000 0.000 35.000
Si
1
Direct
0.00 0.00 0.00
EOF

cat > gap_boundary.poscar << 'EOF'
Boundary Gap Cell (2 atoms along c, 12 Ang gap either way)
1.0
5.430 0.000 0.000
0.000 5.430 0.000
0.000 0.000 24.000
Si
2
Direct
0.00 0.00 0.00
0.00 0.00 0.50
EOF

cat > malformed.poscar << 'EOF'
Malformed Cell
not_a_number
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
Si
1
Direct
0.0 0.0 0.0
EOF

cat > malformed.in << 'EOF'
atom_frac 0.0 0.0 0.0 Si
EOF

echo "Fixtures generated."


# --- 3. Correctness across the 4 readers (cubic lattice, density=0.2 -> 6 6 6) ---
echo -e "\n--- Testing the 4 structure readers (expected grid: 6 6 6) ---"

echo "Testing: POSCAR reader"
stb-kgrid --file cubic.poscar --type poscar --density 0.2 --no-intro > log_reader_poscar.txt 2>&1
check_contains "Suggested Monkhorst-Pack grid: 6 6 6" log_reader_poscar.txt

echo "Testing: CIF reader (ASE)"
stb-kgrid --file cubic.cif --type cif --density 0.2 --no-intro > log_reader_cif.txt 2>&1
check_contains "Suggested Monkhorst-Pack grid: 6 6 6" log_reader_cif.txt

echo "Testing: FHI reader"
stb-kgrid --file cubic.fhi --type fhi --density 0.2 --no-intro > log_reader_fhi.txt 2>&1
check_contains "Suggested Monkhorst-Pack grid: 6 6 6" log_reader_fhi.txt

echo "Testing: FDF reader"
stb-kgrid --file cubic.fdf --type fdf --density 0.2 --no-intro > log_reader_fdf.txt 2>&1
check_contains "Suggested Monkhorst-Pack grid: 6 6 6" log_reader_fdf.txt


# --- 4. Dimensionality heuristic ---
echo -e "\n--- Testing the dimensionality heuristic ---"

echo "Testing: 0D (isolated molecule, expected grid 1 1 1)"
stb-kgrid --file vac_0d.poscar --type poscar --density 0.2 --no-intro > log_dim_0d.txt 2>&1
check_contains "Suggested Monkhorst-Pack grid: 1 1 1" log_dim_0d.txt
check_contains "System appears to be 0D" log_dim_0d.txt

echo "Testing: 1D (chain, expected grid 6 1 1)"
stb-kgrid --file vac_1d.poscar --type poscar --density 0.2 --no-intro > log_dim_1d.txt 2>&1
check_contains "Suggested Monkhorst-Pack grid: 6 1 1" log_dim_1d.txt
check_contains "System appears to be 1D" log_dim_1d.txt

echo "Testing: 2D (slab, expected grid 6 6 1)"
stb-kgrid --file vac_2d.poscar --type poscar --density 0.2 --no-intro > log_dim_2d.txt 2>&1
check_contains "Suggested Monkhorst-Pack grid: 6 6 1" log_dim_2d.txt
check_contains "System appears to be 2D" log_dim_2d.txt

echo "Testing: 3D (cubic bulk, expected grid 6 6 6)"
stb-kgrid --file cubic.poscar --type poscar --density 0.2 --no-intro > log_dim_3d.txt 2>&1
check_contains "System appears to be 3D" log_dim_3d.txt

echo "Testing: real-world fixture structure.fdf (all atoms at z=0.5 -> c is vacuum, expected grid 7 7 1, 2D)"
stb-kgrid --file structure.fdf --type fdf --density 0.2 --no-intro > log_dim_real.txt 2>&1
check_contains "Suggested Monkhorst-Pack grid: 7 7 1" log_dim_real.txt
check_contains "System appears to be 2D" log_dim_real.txt


# --- 4b. --vacuum-gap threshold: same fixture, two different thresholds ---
echo -e "\n--- Testing the --vacuum-gap threshold override ---"
# gap_boundary.poscar: a=b=5.43 Ang (periodic), c=24 Ang with 2 atoms 12 Ang apart
# either way around the ring -> the gap is exactly 12 Ang.

echo "Testing: default --vacuum-gap (10.0) -> 12 Ang gap counts as vacuum (grid 6 6 1, 2D)"
stb-kgrid --file gap_boundary.poscar --type poscar --density 0.2 --no-intro > log_vacuum_gap_default.txt 2>&1
check_contains "Suggested Monkhorst-Pack grid: 6 6 1" log_vacuum_gap_default.txt
check_contains "System appears to be 2D" log_vacuum_gap_default.txt

echo "Testing: --vacuum-gap 15 -> the same 12 Ang gap no longer counts as vacuum (grid 6 6 2, 3D)"
stb-kgrid --file gap_boundary.poscar --type poscar --density 0.2 --vacuum-gap 15 --no-intro > log_vacuum_gap_15.txt 2>&1
check_contains "Suggested Monkhorst-Pack grid: 6 6 2" log_vacuum_gap_15.txt
check_contains "System appears to be 3D" log_vacuum_gap_15.txt


# --- 5. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: missing poscar/fdf/fhi file (returns cleanly, exit 0)"
stb-kgrid --file does_not_exist.poscar --type poscar --density 0.2 --no-intro > log_missing_poscar.txt 2>&1
check_exit_code $? 0
check_contains "not found" log_missing_poscar.txt

stb-kgrid --file does_not_exist.fdf --type fdf --density 0.2 --no-intro > log_missing_fdf.txt 2>&1
check_exit_code $? 0
check_contains "not found" log_missing_fdf.txt

stb-kgrid --file does_not_exist.in --type fhi --density 0.2 --no-intro > log_missing_fhi.txt 2>&1
check_exit_code $? 0
check_contains "not found" log_missing_fhi.txt

echo "Testing: missing cif file (sys.exit(1), different message)"
stb-kgrid --file does_not_exist.cif --type cif --density 0.2 --no-intro > log_missing_cif.txt 2>&1
check_exit_code $? 1
check_contains "Error reading CIF file" log_missing_cif.txt

echo "Testing: malformed poscar (generic parse error, returns cleanly, exit 0)"
stb-kgrid --file malformed.poscar --type poscar --density 0.2 --no-intro > log_malformed_poscar.txt 2>&1
check_exit_code $? 0
check_contains "unexpected error" log_malformed_poscar.txt

echo "Testing: malformed fhi (no lattice_vector lines, sys.exit(1))"
stb-kgrid --file malformed.in --type fhi --density 0.2 --no-intro > log_malformed_fhi.txt 2>&1
check_exit_code $? 1
check_contains "Could not find 3 'lattice_vector' lines" log_malformed_fhi.txt

echo "Testing: density = 0 (ValueError from kspace, sys.exit(1))"
stb-kgrid --file cubic.poscar --type poscar --density 0 --no-intro > log_density_zero.txt 2>&1
check_exit_code $? 1
check_contains "k_density must be positive" log_density_zero.txt

echo "Testing: negative density (same ValueError)"
stb-kgrid --file cubic.poscar --type poscar --density -1 --no-intro > log_density_negative.txt 2>&1
check_exit_code $? 1
check_contains "k_density must be positive" log_density_negative.txt

echo "Testing: invalid --type choice (argparse rejects, exit 2)"
stb-kgrid --file cubic.poscar --type xyz --density 0.2 --no-intro > log_bad_type.txt 2>&1
check_exit_code $? 2
check_contains "invalid choice" log_bad_type.txt

echo "Testing: missing required arguments (--density, --file, --type)"
stb-kgrid --file cubic.poscar --type poscar --no-intro > log_missing_density_arg.txt 2>&1
check_exit_code $? 2
stb-kgrid --type poscar --density 0.2 --no-intro > log_missing_file_arg.txt 2>&1
check_exit_code $? 2
stb-kgrid --file cubic.poscar --density 0.2 --no-intro > log_missing_type_arg.txt 2>&1
check_exit_code $? 2

echo "Testing: -v/--version"
stb-kgrid --version > log_version.txt 2>&1
check_contains "stb-kgrid" log_version.txt

echo "Testing: --help documents --vacuum-gap with its default"
stb-kgrid --help > log_help.txt 2>&1
check_contains "vacuum-gap" log_help.txt
check_contains "Default: 10.0" log_help.txt


# --- 6. Interactive path (stb-suite, shortcut 1.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 1.2) ---"

echo "Testing: navigate 1.2 -> invalid file then valid -> type picked from a numbered list (out of range then valid=1/fdf) -> invalid density then valid (0.2) -> quit"
printf '1.2\ndoes_not_exist.fdf\nstructure.fdf\n9\n1\n-1\n0.2\n\n0\n' | stb-suite > log_interactive.txt 2>&1
check_contains "File not found" log_interactive.txt
check_contains "Invalid choice! Please select between 1 and 4" log_interactive.txt
check_contains "as type 'fdf'" log_interactive.txt
check_contains "Density must be a positive number" log_interactive.txt
check_contains "Suggested Monkhorst-Pack grid: 7 7 1" log_interactive.txt

echo "Testing: navigate 1.2 -> valid file -> press Enter for both type and density (defaults: fdf, 0.2) -> quit"
printf '1.2\nstructure.fdf\n\n\n0\n' | stb-suite > log_interactive_defaults.txt 2>&1
check_contains "as type 'fdf'" log_interactive_defaults.txt
check_contains "Suggested Monkhorst-Pack grid: 7 7 1" log_interactive_defaults.txt


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
