#!/bin/bash

# --- Setup ---
# Smoke test for stb-kgrid (K-Grid Generator, item 1.3)
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

# Checks that file $2 does NOT contain (grep -q) pattern $1
check_not_contains() {
    if ! grep -q -- "$1" "$2" 2>/dev/null; then
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
echo "--- Starting tester for STB-KGrid (item 1.3) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$STRUCT" "$TEST_DIR/structure.fdf"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null

# --- 2. Fixtures ---
# stb-kgrid is fdf-only now (poscar/cif/fhi support removed) -- every
# fixture below is a plain SIESTA .fdf.
echo -e "\n--- Generating fixtures ---"

# Cubic 5.43 Ang Si lattice: 2*pi/5.43 = 1.157 1/Ang; ceil(1.157/0.2) = 6.
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
cat > vac_0d.fdf << 'EOF'
NumberOfSpecies    1
NumberofAtoms      1

%block ChemicalSpeciesLabel
 1  14  Si
%endblock ChemicalSpeciesLabel

LatticeConstant 1.0 Ang

AtomicCoordinatesFormat  Fractional

%block LatticeVectors
 40.000   0.000   0.000
 0.000   40.000   0.000
 0.000   0.000   40.000
%endblock LatticeVectors

%block AtomicCoordinatesAndAtomicSpecies
 0.00 0.00 0.00 1
%endblock AtomicCoordinatesAndAtomicSpecies
EOF

cat > vac_1d.fdf << 'EOF'
NumberOfSpecies    1
NumberofAtoms      1

%block ChemicalSpeciesLabel
 1  14  Si
%endblock ChemicalSpeciesLabel

LatticeConstant 1.0 Ang

AtomicCoordinatesFormat  Fractional

%block LatticeVectors
 5.430   0.000   0.000
 0.000   35.000   0.000
 0.000   0.000   35.000
%endblock LatticeVectors

%block AtomicCoordinatesAndAtomicSpecies
 0.00 0.00 0.00 1
%endblock AtomicCoordinatesAndAtomicSpecies
EOF

cat > vac_2d.fdf << 'EOF'
NumberOfSpecies    1
NumberofAtoms      1

%block ChemicalSpeciesLabel
 1  14  Si
%endblock ChemicalSpeciesLabel

LatticeConstant 1.0 Ang

AtomicCoordinatesFormat  Fractional

%block LatticeVectors
 5.430   0.000   0.000
 0.000   5.430   0.000
 0.000   0.000   35.000
%endblock LatticeVectors

%block AtomicCoordinatesAndAtomicSpecies
 0.00 0.00 0.00 1
%endblock AtomicCoordinatesAndAtomicSpecies
EOF

cat > gap_boundary.fdf << 'EOF'
NumberOfSpecies    1
NumberofAtoms      2

%block ChemicalSpeciesLabel
 1  14  Si
%endblock ChemicalSpeciesLabel

LatticeConstant 1.0 Ang

AtomicCoordinatesFormat  Fractional

%block LatticeVectors
 5.430   0.000   0.000
 0.000   5.430   0.000
 0.000   0.000   24.000
%endblock LatticeVectors

%block AtomicCoordinatesAndAtomicSpecies
 0.00 0.00 0.00 1
 0.00 0.00 0.50 1
%endblock AtomicCoordinatesAndAtomicSpecies
EOF

# Degenerate lattice (two identical vectors -> zero volume), same pattern
# already used by test/1-inputs/1-input_file/test.sh's own zero-volume case.
cat > malformed.fdf << 'EOF'
NumberOfSpecies    1
NumberofAtoms      1

%block ChemicalSpeciesLabel
 1  14  Si
%endblock ChemicalSpeciesLabel

LatticeConstant 1.0 Ang

AtomicCoordinatesFormat  Fractional

%block LatticeVectors
 5.0   0.0   0.0
 5.0   0.0   0.0
 0.0   0.0   5.0
%endblock LatticeVectors

%block AtomicCoordinatesAndAtomicSpecies
 0.0   0.0   0.0   1
%endblock AtomicCoordinatesAndAtomicSpecies
EOF

echo "Fixtures generated."


# --- 3. Base correctness (cubic lattice, density=0.2 -> 6 6 6) ---
echo -e "\n--- Testing base correctness ---"

echo "Testing: cubic.fdf at density 0.2"
stb-kgrid --file cubic.fdf --density 0.2 --no-intro > log_reader_fdf.txt 2>&1
check_contains "Suggested Monkhorst-Pack grid : 6 6 6" log_reader_fdf.txt
check_contains "Dimensionality : 3D" log_reader_fdf.txt


# --- 4. Dimensionality heuristic ---
echo -e "\n--- Testing the dimensionality heuristic ---"

echo "Testing: 0D (isolated molecule, expected grid 1 1 1)"
stb-kgrid --file vac_0d.fdf --density 0.2 --no-intro > log_dim_0d.txt 2>&1
check_contains "Suggested Monkhorst-Pack grid : 1 1 1" log_dim_0d.txt
check_contains "Dimensionality : 0D" log_dim_0d.txt

echo "Testing: 1D (chain, expected grid 6 1 1)"
stb-kgrid --file vac_1d.fdf --density 0.2 --no-intro > log_dim_1d.txt 2>&1
check_contains "Suggested Monkhorst-Pack grid : 6 1 1" log_dim_1d.txt
check_contains "Dimensionality : 1D" log_dim_1d.txt

echo "Testing: 2D (slab, expected grid 6 6 1)"
stb-kgrid --file vac_2d.fdf --density 0.2 --no-intro > log_dim_2d.txt 2>&1
check_contains "Suggested Monkhorst-Pack grid : 6 6 1" log_dim_2d.txt
check_contains "Dimensionality : 2D" log_dim_2d.txt

echo "Testing: 3D (cubic bulk, expected grid 6 6 6)"
stb-kgrid --file cubic.fdf --density 0.2 --no-intro > log_dim_3d.txt 2>&1
check_contains "Dimensionality : 3D" log_dim_3d.txt

echo "Testing: real-world fixture structure.fdf (all atoms at z=0.5 -> c is vacuum, expected grid 7 7 1, 2D)"
stb-kgrid --file structure.fdf --density 0.2 --no-intro > log_dim_real.txt 2>&1
check_contains "Suggested Monkhorst-Pack grid : 7 7 1" log_dim_real.txt
check_contains "Dimensionality : 2D" log_dim_real.txt


# --- 4b. --vacuum-gap threshold: same fixture, two different thresholds ---
echo -e "\n--- Testing the --vacuum-gap threshold override ---"
# gap_boundary.fdf: a=b=5.43 Ang (periodic), c=24 Ang with 2 atoms 12 Ang apart
# either way around the ring -> the gap is exactly 12 Ang.

echo "Testing: default --vacuum-gap (10.0) -> 12 Ang gap counts as vacuum (grid 6 6 1, 2D)"
stb-kgrid --file gap_boundary.fdf --density 0.2 --no-intro > log_vacuum_gap_default.txt 2>&1
check_contains "Suggested Monkhorst-Pack grid : 6 6 1" log_vacuum_gap_default.txt
check_contains "Dimensionality : 2D" log_vacuum_gap_default.txt

echo "Testing: --vacuum-gap 15 -> the same 12 Ang gap no longer counts as vacuum (grid 6 6 2, 3D)"
stb-kgrid --file gap_boundary.fdf --density 0.2 --vacuum-gap 15 --no-intro > log_vacuum_gap_15.txt 2>&1
check_contains "Suggested Monkhorst-Pack grid : 6 6 2" log_vacuum_gap_15.txt
check_contains "Dimensionality : 3D" log_vacuum_gap_15.txt


# --- 5. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: missing structure file (exit 1)"
stb-kgrid --file does_not_exist.fdf --density 0.2 --no-intro > log_missing_fdf.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_fdf.txt

echo "Testing: malformed fdf (degenerate/zero-volume lattice, exit 1)"
stb-kgrid --file malformed.fdf --density 0.2 --no-intro > log_malformed.txt 2>&1
check_exit_code $? 1
check_contains "Cell volume is zero" log_malformed.txt

echo "Testing: density = 0 (ValueError from kspace, exit 1)"
stb-kgrid --file cubic.fdf --density 0 --no-intro > log_density_zero.txt 2>&1
check_exit_code $? 1
check_contains "k_density must be positive" log_density_zero.txt

echo "Testing: negative density (same ValueError)"
stb-kgrid --file cubic.fdf --density -1 --no-intro > log_density_negative.txt 2>&1
check_exit_code $? 1
check_contains "k_density must be positive" log_density_negative.txt

echo "Testing: missing required arguments (--density, --file)"
stb-kgrid --file cubic.fdf --no-intro > log_missing_density_arg.txt 2>&1
check_exit_code $? 2
stb-kgrid --density 0.2 --no-intro > log_missing_file_arg.txt 2>&1
check_exit_code $? 2

echo "Testing: -v/--version"
stb-kgrid --version > log_version.txt 2>&1
check_contains "stb-kgrid" log_version.txt

echo "Testing: --help no longer mentions --type (removed) but still documents --vacuum-gap"
stb-kgrid --help > log_help.txt 2>&1
check_not_contains "--type" log_help.txt
check_contains "vacuum-gap" log_help.txt
check_contains "Default: 10.0" log_help.txt


# --- 5b. --save-report and references.bib ---
echo -e "\n--- Testing --save-report and references.bib ---"

echo "Testing: --save-report writes stb_kgrid_report.txt matching the console output"
rm -f stb_kgrid_report.txt
stb-kgrid --file cubic.fdf --density 0.2 --save-report --no-intro > log_save_report_console.txt 2>&1
check_success stb_kgrid_report.txt
check_contains "Suggested Monkhorst-Pack grid : 6 6 6" stb_kgrid_report.txt
check_contains "Report" log_save_report_console.txt

echo "Testing: references.bib always written, with SIESTA + Monkhorst-Pack citations"
rm -f references.bib
stb-kgrid --file cubic.fdf --density 0.2 --no-intro > /dev/null 2>&1
check_success references.bib
check_contains "@article{Soler2002," references.bib
check_contains "@article{Garcia2020," references.bib
check_contains "@article{MonkhorstPack1976," references.bib


# --- 6. Interactive path (stb-suite, shortcut 1.3) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 1.3) ---"

echo "Testing: navigate 1.3 -> invalid file then valid -> invalid density then valid (0.2) -> quit"
printf '1.3\ndoes_not_exist.fdf\nstructure.fdf\n-1\n0.2\n\n0\n' | stb-suite > log_interactive.txt 2>&1
check_contains "File not found" log_interactive.txt
check_contains "K-Point Density Recommendation Guide" log_interactive.txt
check_contains "Density must be a positive number" log_interactive.txt
check_contains "Suggested Monkhorst-Pack grid : 7 7 1" log_interactive.txt

echo "Testing: density guide appears BEFORE the density prompt (not after)"
GUIDE_LINE=$(grep -n "K-Point Density Recommendation Guide" log_interactive.txt | head -1 | cut -d: -f1)
PROMPT_LINE=$(grep -n "K-point density (e.g., 0.2)" log_interactive.txt | head -1 | cut -d: -f1)
if [ -n "$GUIDE_LINE" ] && [ -n "$PROMPT_LINE" ] && [ "$GUIDE_LINE" -lt "$PROMPT_LINE" ]; then
    echo -e "   -> ${GREEN}Verified:${NC} guide (line $GUIDE_LINE) appears before the density prompt (line $PROMPT_LINE)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} guide/prompt ordering wrong (guide=$GUIDE_LINE, prompt=$PROMPT_LINE)"
    FAIL=$((FAIL+1))
fi

echo "Testing: navigate 1.3 -> valid file -> press Enter for density (default: 0.2) -> quit"
printf '1.3\nstructure.fdf\n\n\n0\n' | stb-suite > log_interactive_defaults.txt 2>&1
check_contains "Suggested Monkhorst-Pack grid : 7 7 1" log_interactive_defaults.txt


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
