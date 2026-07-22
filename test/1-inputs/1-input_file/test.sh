#!/bin/bash

# --- Setup ---
# Smoke test for stb-inputfile (Input File Generator, item 1.1)
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

# Checks that the last command failed (exit code != 0)
check_failure_exit() {
    if [ "$1" -ne 0 ]; then
        echo -e "   -> ${GREEN}Verified:${NC} command failed as expected (exit $1)"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} command should have failed, but returned 0"
        FAIL=$((FAIL+1))
    fi
}

# --- 1. Preparation ---
echo "--- Starting tester for STB-InputFile (item 1.1) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$STRUCT" "$TEST_DIR/structure.fdf"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null

# --- 2. Matrix of the 4 calculation modes x DFT-D3 on/off ---
# DFT-D3 is a standalone flag (-d3/--d3), independent of -t/--type -- not a
# "<mode>+d3" suffix on the mode string anymore.
echo -e "\n--- Testing the 4 calculation modes, with and without -d3 ---"

MODES=(total_energy relax aimd bands)

for mode in "${MODES[@]}"; do
    for d3_flag in "" "-d3"; do
        label="$mode"
        [ -n "$d3_flag" ] && label="${mode}_d3"
        echo -n "Testing mode: $mode (d3=${d3_flag:-off})"
        rm -f calc.fdf
        stb-inputfile structure.fdf -t "$mode" $d3_flag --no-intro > "log_${label}.txt" 2>&1
        check_success calc.fdf
        mv calc.fdf "calc_${label}.fdf"
    done
done

echo -e "\n--- Verifying generated content ---"

echo "Checking %include (should use only the basename):"
check_contains "%include structure.fdf" calc_total_energy.fdf

echo "Checking DFTD3 (.false. for modes without +d3):"
check_contains "DFTD3                   .false." calc_relax.fdf
check_contains "DFTD3                   .false." calc_total_energy.fdf
check_contains "DFTD3                   .false." calc_aimd.fdf
check_contains "DFTD3                   .false." calc_bands.fdf

echo "Checking DFTD3 (.true. for +d3 modes):"
check_contains "DFTD3                   .true." calc_relax_d3.fdf
check_contains "DFTD3                   .true." calc_total_energy_d3.fdf
check_contains "DFTD3                   .true." calc_aimd_d3.fdf
check_contains "DFTD3                   .true." calc_bands_d3.fdf

echo "Checking that AIMD's k-grid is preserved from the template ([1  1  1]):"
check_contains "kgrid.MonkhorstPack   \[1  1  1\]" calc_aimd.fdf

echo "Checking that the k-grid was recalculated (should not be the template default) for relax:"
check_not_contains "kgrid.MonkhorstPack   \[6  6  6\]" calc_relax.fdf

echo "Checking that the k-grid was recalculated for total_energy:"
check_not_contains "kgrid.MonkhorstPack   \[4  13  1\]" calc_total_energy.fdf

echo "Checking the 'relax' template's unique marker (MD.TypeOfRun CG):"
check_contains "MD.TypeOfRun            CG" calc_relax.fdf

echo "Checking the 'aimd' template's unique markers (MD.TypeOfRun Nose + active Supercell block):"
check_contains "MD.TypeOfRun  Nose" calc_aimd.fdf
check_contains "^%block Supercell" calc_aimd.fdf

echo "Checking the 'bands' template's unique markers (kpath_bs.fdf + PDOS):"
check_contains "%include kpath_bs.fdf" calc_bands.fdf
check_contains "%block  ProjectedDensityOfStates" calc_bands.fdf

echo "Checking that 'total_energy' has NONE of the relax/aimd/bands active markers:"
check_not_contains "MD.TypeOfRun            CG" calc_total_energy.fdf
check_not_contains "^%block Supercell" calc_total_energy.fdf
check_not_contains "%include kpath_bs.fdf" calc_total_energy.fdf

echo "Checking that Spin defaults to non-polarized (no -s given):"
check_contains "Spin                non-polarized" calc_relax.fdf


# --- 2b. Spin polarization (-s / --spin-polarized), independent of -t and -d3 ---
echo -e "\n--- Testing -s / --spin-polarized ---"

echo -n "Testing: -s enables spin polarization"
rm -f calc.fdf
stb-inputfile structure.fdf -t relax -s --no-intro > log_spin.txt 2>&1
check_success calc.fdf
check_contains "Spin                polarized" calc.fdf
mv calc.fdf calc_spin.fdf

echo -n "Testing: -s combined with -d3 (both independent flags at once)"
rm -f calc.fdf
stb-inputfile structure.fdf -t relax -s -d3 --no-intro > log_spin_d3.txt 2>&1
check_success calc.fdf
check_contains "Spin                polarized" calc.fdf
check_contains "DFTD3                   .true." calc.fdf
mv calc.fdf calc_spin_d3.fdf


# --- 2c. Structure validation (symmetry + malformation warnings), section [2] ---
echo -e "\n--- Testing [2] STRUCTURE VALIDATION ---"

echo -n "Testing: normal structure -> symmetry reported, explicit per-check [OK] table"
rm -f calc.fdf
stb-inputfile structure.fdf -t total_energy --no-intro > log_validation_ok.txt 2>&1
check_contains "Space group :" log_validation_ok.txt
check_contains "Crystal system :" log_validation_ok.txt
check_contains "Atom proximity" log_validation_ok.txt
check_contains "Lattice handedness" log_validation_ok.txt
check_contains "Atomic density" log_validation_ok.txt
check_contains "\[OK\]" log_validation_ok.txt
check_not_contains "Unusual atomic density" log_validation_ok.txt
rm -f calc.fdf

echo -n "Testing: two atoms closer than 0.5 Ang -> too-close-atoms warning"
cat > structure_close_atoms.fdf << 'EOF'
NumberOfSpecies    1
NumberofAtoms      2

%block ChemicalSpeciesLabel
 1   6   C
%endblock ChemicalSpeciesLabel

LatticeConstant 1.0 Ang

AtomicCoordinatesFormat  Ang

%block LatticeVectors
 10.0   0.0   0.0
 0.0   10.0   0.0
 0.0   0.0   10.0
%endblock LatticeVectors

%block AtomicCoordinatesAndAtomicSpecies
 5.0   5.0   5.0   1
 5.1   5.0   5.0   1
%endblock AtomicCoordinatesAndAtomicSpecies
EOF
stb-inputfile structure_close_atoms.fdf -t total_energy --no-intro > log_validation_close.txt 2>&1
check_contains "unusually close" log_validation_close.txt
rm -f calc.fdf

echo -n "Testing: left-handed lattice (negative determinant) -> warning"
cat > structure_left_handed.fdf << 'EOF'
NumberOfSpecies    1
NumberofAtoms      1

%block ChemicalSpeciesLabel
 1   6   C
%endblock ChemicalSpeciesLabel

LatticeConstant 1.0 Ang

AtomicCoordinatesFormat  Fractional

%block LatticeVectors
 5.0   0.0   0.0
 0.0   5.0   0.0
 0.0   0.0  -5.0
%endblock LatticeVectors

%block AtomicCoordinatesAndAtomicSpecies
 0.0   0.0   0.0   1
%endblock AtomicCoordinatesAndAtomicSpecies
EOF
stb-inputfile structure_left_handed.fdf -t total_energy --no-intro > log_validation_lefthanded.txt 2>&1
check_contains "left-handed" log_validation_lefthanded.txt
rm -f calc.fdf

echo -n "Testing: implausibly low density in a genuine 3D bulk cell -> warning"
cat > structure_dilute.fdf << 'EOF'
NumberOfSpecies    1
NumberofAtoms      8

%block ChemicalSpeciesLabel
 1  14  Si
%endblock ChemicalSpeciesLabel

LatticeConstant 1.0 Ang

AtomicCoordinatesFormat  Fractional

%block LatticeVectors
 13.575   0.00000   0.00000
 0.00000   13.575   0.00000
 0.00000   0.00000   13.575
%endblock LatticeVectors

%block AtomicCoordinatesAndAtomicSpecies
  0.50000000   0.00000000   0.50000000   1
  0.25000000   0.25000000   0.25000000   1
  0.00000000   0.00000000   0.00000000   1
  0.25000000   0.75000000   0.75000000   1
  0.75000000   0.25000000   0.75000000   1
  0.00000000   0.50000000   0.50000000   1
  0.50000000   0.50000000   0.00000000   1
  0.75000000   0.75000000   0.25000000   1
%endblock AtomicCoordinatesAndAtomicSpecies
EOF
stb-inputfile structure_dilute.fdf -t total_energy --no-intro > log_validation_dilute.txt 2>&1
check_contains "Unusual atomic density" log_validation_dilute.txt
rm -f calc.fdf

echo -n "Testing: density check skipped for a vacuum-padded structure (would otherwise trigger)"
cat > structure_vacuum_dilute.fdf << 'EOF'
NumberOfSpecies    1
NumberofAtoms      1

%block ChemicalSpeciesLabel
 1   6   C
%endblock ChemicalSpeciesLabel

LatticeConstant 1.0 Ang

AtomicCoordinatesFormat  Ang

%block LatticeVectors
 3.0   0.0   0.0
 0.0   3.0   0.0
 0.0   0.0   20.0
%endblock LatticeVectors

%block AtomicCoordinatesAndAtomicSpecies
 1.5   1.5   10.0   1
%endblock AtomicCoordinatesAndAtomicSpecies
EOF
stb-inputfile structure_vacuum_dilute.fdf -t total_energy --no-intro > log_validation_vacuum_dilute.txt 2>&1
check_not_contains "Unusual atomic density" log_validation_vacuum_dilute.txt
rm -f calc.fdf

echo -n "Testing: NumberOfAtoms header mismatch -> warning"
cat > structure_bad_header.fdf << 'EOF'
NumberOfSpecies    1
NumberofAtoms      5

%block ChemicalSpeciesLabel
 1   6   C
%endblock ChemicalSpeciesLabel

LatticeConstant 1.0 Ang

AtomicCoordinatesFormat  Fractional

%block LatticeVectors
 5.0   0.0   0.0
 0.0   5.0   0.0
 0.0   0.0   5.0
%endblock LatticeVectors

%block AtomicCoordinatesAndAtomicSpecies
 0.0   0.0   0.0   1
%endblock AtomicCoordinatesAndAtomicSpecies
EOF
stb-inputfile structure_bad_header.fdf -t total_energy --no-intro > log_validation_badheader.txt 2>&1
check_contains "Declared NumberOfAtoms" log_validation_badheader.txt
rm -f calc.fdf

echo -n "Testing: --view is accepted and does not crash (no real display forced via DISPLAY=)"
rm -f calc.fdf
DISPLAY= stb-inputfile structure.fdf -t relax --view --no-intro > log_view.txt 2>&1
EXIT_CODE=$?
if [ "$EXIT_CODE" -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} --view did not crash the run (exit 0)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} --view run exited with $EXIT_CODE"
    FAIL=$((FAIL+1))
fi
check_success calc.fdf
mv calc.fdf calc_view.fdf


# --- 3. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo -n "Testing: nonexistent structure_file"
rm -f calc.fdf
stb-inputfile does_not_exist.fdf -t relax --no-intro > log_missing_struct.txt 2>&1
check_contains "not found" log_missing_struct.txt
if [ ! -s calc.fdf ]; then
    echo -e "   -> ${GREEN}Verified:${NC} calc.fdf was not created"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} calc.fdf was created despite a missing structure_file"
    FAIL=$((FAIL+1))
fi

echo -n "Testing: -t with an invalid value"
stb-inputfile structure.fdf -t invalid_mode --no-intro > log_bad_mode.txt 2>&1
check_failure_exit $?

echo -n "Testing: old '<mode>+d3' suffix syntax is now rejected (-d3 is its own flag)"
stb-inputfile structure.fdf -t relax+d3 --no-intro > log_old_suffix.txt 2>&1
check_failure_exit $?
check_contains "invalid choice" log_old_suffix.txt

echo -n "Testing: -t missing (required)"
stb-inputfile structure.fdf --no-intro > log_missing_type.txt 2>&1
check_failure_exit $?
check_contains "required" log_missing_type.txt

echo -n "Testing: -v/--version"
stb-inputfile --version > log_version.txt 2>&1
check_contains "stb-inputfile" log_version.txt

echo -n "Testing: structure_file inside a subdirectory (only the basename should reach %include)"
mkdir -p subdir
cp structure.fdf subdir/structure.fdf
rm -f calc.fdf
stb-inputfile subdir/structure.fdf -t relax --no-intro > log_subdir.txt 2>&1
check_success calc.fdf
check_contains "%include structure.fdf" calc.fdf
check_not_contains "%include subdir/structure.fdf" calc.fdf
mv calc.fdf calc_subdir.fdf

echo -n "Testing: degenerate lattice (zero volume) -> handled error, no crash"
cat > structure_zero_volume.fdf << 'EOF'
NumberOfSpecies    1
NumberofAtoms      1

%block ChemicalSpeciesLabel
 1   6   C
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
rm -f calc.fdf
stb-inputfile structure_zero_volume.fdf -t relax --no-intro > log_zero_volume.txt 2>&1
EXIT_CODE=$?
if [ "$EXIT_CODE" -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} process did not crash (exit 0, error handled internally)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} process exited with $EXIT_CODE (expected 0 with a handled error)"
    FAIL=$((FAIL+1))
fi
check_contains "ERROR" log_zero_volume.txt
if [ ! -s calc.fdf ]; then
    echo -e "   -> ${GREEN}Verified:${NC} calc.fdf was not created after the error"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} calc.fdf was created despite the zero-volume error"
    FAIL=$((FAIL+1))
fi


# --- 4. Pseudopotential copying (-p / --pp-path) ---
echo -e "\n--- Testing pseudopotential copying (-p) ---"

echo -n "Testing: -p pointing to a nonexistent directory"
rm -f calc.fdf
stb-inputfile structure.fdf -t relax -p ./pp_does_not_exist --no-intro > log_pp_missing_dir.txt 2>&1
check_contains "not a valid directory" log_pp_missing_dir.txt
check_success calc.fdf
mv calc.fdf calc_pp_missing_dir.fdf

echo -n "Testing: -p valid but missing the species' .psml/.psf (C)"
mkdir -p pp_empty
rm -f calc.fdf
stb-inputfile structure.fdf -t relax -p ./pp_empty --no-intro > log_pp_empty.txt 2>&1
check_contains "WARNING" log_pp_empty.txt
check_success calc.fdf
mv calc.fdf calc_pp_empty.fdf

echo -n "Testing: -p valid with C.psf present -> should be copied"
mkdir -p pp_with_psf
echo "# dummy psf" > pp_with_psf/C.psf
rm -f calc.fdf C.psf C.psml
stb-inputfile structure.fdf -t relax -p ./pp_with_psf --no-intro > log_pp_psf.txt 2>&1
check_success calc.fdf
check_success C.psf
mv calc.fdf calc_pp_psf.fdf

echo -n "Testing: -p with both C.psml AND C.psf -> .psml takes priority"
mkdir -p pp_with_both
echo "# dummy psml" > pp_with_both/C.psml
echo "# dummy psf" > pp_with_both/C.psf
rm -f calc.fdf C.psf C.psml
stb-inputfile structure.fdf -t relax -p ./pp_with_both --no-intro > log_pp_both.txt 2>&1
check_success calc.fdf
check_success C.psml

echo -n "Testing: -p accepts a bundled bank name (dojo)"
rm -f calc.fdf C.psf C.psml
stb-inputfile structure.fdf -t relax -p dojo --no-intro > log_pp_bank.txt 2>&1
check_success calc.fdf
check_success C.psml
if [ ! -s C.psf ]; then
    echo -e "   -> ${GREEN}Verified:${NC} C.psf was not copied (psml took priority)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} C.psf was copied even though C.psml already existed"
    FAIL=$((FAIL+1))
fi
mv calc.fdf calc_pp_both.fdf
rm -f C.psf C.psml


# --- 5. Interactive path (stb-suite, shortcut 1.1) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 1.1) ---"

echo -n "Testing: navigate 1.1 -> invalid file then valid -> mode out of range then valid (2=relax) -> d3=y -> spin=y -> skip PP -> view=n -> quit"
rm -f calc.fdf
printf '1.1\ndoes_not_exist.fdf\nstructure.fdf\n99\n2\ny\ny\n\nn\nn\n\n0\n' | stb-suite > log_interactive.txt 2>&1
check_success calc.fdf
check_contains "File not found" log_interactive.txt
check_contains "Invalid choice" log_interactive.txt
check_contains "DFTD3                   .true." calc.fdf
check_contains "Spin                polarized" calc.fdf
mv calc.fdf calc_interactive.fdf


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
