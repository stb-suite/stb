#!/bin/bash

# --- Setup ---
# Smoke test for stb-kpath (K-Path Generator, item 1.3)
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
    if grep -q "$1" "$2" 2>/dev/null; then
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
echo "--- Starting tester for STB-KPath (item 1.3) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/structure.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/example_0D_molecule.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/example_1D_chain.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/example_2D_graphene.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/example_3D_silicon.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null

cat > malformed.fdf << 'EOF'
NumberOfSpecies    1
NumberofAtoms      1
LatticeConstant 1.0 Ang
AtomicCoordinatesFormat  Fractional
%block LatticeVectors
 5.0 0.0 0.0
 0.0 5.0 0.0
 0.0 0.0 5.0
%endblock LatticeVectors
%block AtomicCoordinatesAndAtomicSpecies
 0.0 0.0 0.0 1
%endblock AtomicCoordinatesAndAtomicSpecies
EOF


# --- 2. Dimensionality matrix (0D/1D/2D/3D) ---
echo -e "\n--- Testing dimensionality detection ---"

echo "Testing: 0D (isolated CH4 molecule) -> no periodic direction, no file written"
rm -f kpath_bs.fdf
stb-kpath --file example_0D_molecule.fdf --no-intro > log_0d.txt 2>&1
check_contains "0D (isolated molecule)" log_0d.txt
check_contains "Skipping file generation" log_0d.txt
if [ ! -s kpath_bs.fdf ]; then
    echo -e "   -> ${GREEN}Verified:${NC} kpath_bs.fdf was not created"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} kpath_bs.fdf was created for a 0D system"
    FAIL=$((FAIL+1))
fi

echo "Testing: 1D (linear carbon chain) -> expected path Gamma-X"
rm -f kpath_bs.fdf
stb-kpath --file example_1D_chain.fdf --no-intro > log_1d.txt 2>&1
check_contains "Dimensionality:.*1D" log_1d.txt
check_contains "LINE" log_1d.txt
check_contains "Path:.*Gamma-X" log_1d.txt
check_success kpath_bs.fdf
mv kpath_bs.fdf kpath_1d.fdf

echo "Testing: 2D (primitive graphene, 2 atoms) -> expected path Gamma-M-K-Gamma"
rm -f kpath_bs.fdf
stb-kpath --file example_2D_graphene.fdf --no-intro > log_2d.txt 2>&1
check_contains "Dimensionality:.*2D" log_2d.txt
check_contains "HEX2D" log_2d.txt
check_contains "Path:.*Gamma-M-K-.Gamma" log_2d.txt
check_contains "not a 3D space group" log_2d.txt
check_success kpath_bs.fdf
mv kpath_bs.fdf kpath_2d.fdf

echo "Testing: 2D real-world fixture structure.fdf (10-atom supercell) -> same Gamma-M-K-Gamma, no primitive-cell warning"
rm -f kpath_bs.fdf
stb-kpath --file structure.fdf --no-intro > log_2d_real.txt 2>&1
check_contains "Dimensionality:.*2D" log_2d_real.txt
check_contains "HEX2D" log_2d_real.txt
check_contains "Path:.*Gamma-M-K-.Gamma" log_2d_real.txt
check_success kpath_bs.fdf
mv kpath_bs.fdf kpath_2d_real.fdf

echo "Testing: 3D (simple-cubic Si fixture) -> expected path Gamma-X-M-Gamma-R-X | M-R"
rm -f kpath_bs.fdf
stb-kpath --file example_3D_silicon.fdf --no-intro > log_3d.txt 2>&1
check_contains "Dimensionality:.*3D" log_3d.txt
check_contains "CUB" log_3d.txt
check_contains "Path:.*Gamma-X-M-.Gamma-R-X | M-R" log_3d.txt
check_contains "Disjointed path detected" log_3d.txt
check_success kpath_bs.fdf
mv kpath_bs.fdf kpath_3d.fdf


# --- 3. Content of the generated kpath_bs.fdf ---
echo -e "\n--- Verifying generated kpath_bs.fdf content ---"
check_contains "%block BandLines" kpath_2d.fdf
check_contains "%endblock BandLines" kpath_2d.fdf
check_contains "Gamma" kpath_2d.fdf
check_contains "  M$" kpath_2d.fdf


# --- 4. --vacuum-gap threshold override ---
echo -e "\n--- Testing the --vacuum-gap threshold override ---"

echo "Testing: raising --vacuum-gap above the 20 Ang vacuum -> falls back to full 3D (old behavior)"
rm -f kpath_bs.fdf
stb-kpath --file example_2D_graphene.fdf --vacuum-gap 25 --no-intro > log_vacuum_gap_override.txt 2>&1
check_contains "Dimensionality:.*3D" log_vacuum_gap_override.txt
check_success kpath_bs.fdf
rm -f kpath_bs.fdf


# --- 5. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: missing file"
stb-kpath --file does_not_exist.fdf --no-intro > log_missing.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing.txt

echo "Testing: malformed fdf (missing ChemicalSpeciesLabel) -> handled, no crash"
rm -f kpath_bs.fdf
stb-kpath --file malformed.fdf --no-intro > log_malformed.txt 2>&1
check_exit_code $? 0
check_contains "ChemicalSpeciesLabel not found" log_malformed.txt
if [ ! -s kpath_bs.fdf ]; then
    echo -e "   -> ${GREEN}Verified:${NC} kpath_bs.fdf was not created after a parse error"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} kpath_bs.fdf was created despite the parse error"
    FAIL=$((FAIL+1))
fi

echo "Testing: -f/--file missing (required)"
stb-kpath --no-intro > log_missing_file_arg.txt 2>&1
check_exit_code $? 2

echo "Testing: -v/--version"
stb-kpath --version > log_version.txt 2>&1
check_contains "stb-kpath" log_version.txt

echo "Testing: --help documents --vacuum-gap and -p/--prec defaults"
stb-kpath --help > log_help.txt 2>&1
check_contains "vacuum-gap" log_help.txt
check_contains "Default: 10.0" log_help.txt
check_contains "Default: 0.0002" log_help.txt


# --- 6. Interactive path (stb-suite, shortcut 1.3) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 1.3) ---"

echo "Testing: navigate 1.3 -> invalid file then valid -> press Enter for precision (default 0.0002) -> quit"
rm -f kpath_bs.fdf
printf '1.3\ndoes_not_exist.fdf\nexample_2D_graphene.fdf\n\n0\n' | stb-suite > log_interactive.txt 2>&1
check_contains "File not found" log_interactive.txt
check_contains "Dimensionality:.*2D" log_interactive.txt
check_contains "Path:.*Gamma-M-K-.Gamma" log_interactive.txt
check_success kpath_bs.fdf


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
