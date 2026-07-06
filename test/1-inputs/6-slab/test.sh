#!/bin/bash

# --- Setup ---
# Smoke test for stb-slab (Slab Builder, item 1.6)
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
echo "--- Starting tester for STB-Slab (item 1.6) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/structure.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/rocksalt.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Default path: single-element cubic fixture, (1,0,0) ---
echo -e "\n--- Testing default termination selection (structure.fdf, hkl 1 0 0) ---"
rm -f slab.fdf
stb-slab -f structure.fdf --hkl 1 0 0 --no-intro > log_default.txt 2>&1
check_contains "Terminations found:.*2" log_default.txt
check_contains "Slab Thickness.*9.50 Ang" log_default.txt
check_contains "Vacuum Thickness.*17.65 Ang" log_default.txt
check_success slab.fdf
check_contains "NumberofAtoms      4" slab.fdf


# --- 2b. --passivate (cut + passivate in one step, reuses core/passivation.py) ---
echo -e "\n--- Testing --passivate (structure.fdf, hkl 1 1 1) ---"
rm -f slab_passivated.fdf
stb-slab -f structure.fdf --hkl 1 1 1 --passivate -o slab_passivated.fdf --no-intro \
    > log_passivate.txt 2>&1
check_contains "Dangling bonds found:" log_passivate.txt
check_contains "Auto-passivated:.* with H" log_passivate.txt
check_success slab_passivated.fdf
check_contains " 2   1   H" slab_passivated.fdf

echo "Testing: --passivant rejected without --passivate"
stb-slab -f structure.fdf --hkl 1 0 0 --passivant F --no-intro > log_passivant_bad.txt 2>&1
check_exit_code $? 2
check_contains "only valid with --passivate" log_passivant_bad.txt


# --- 3. --all writes every termination ---
echo -e "\n--- Testing --all (writes one file per termination) ---"
rm -f slab_term0.fdf slab_term1.fdf
stb-slab -f structure.fdf --hkl 1 0 0 --all --no-intro > log_all.txt 2>&1
check_contains "Termination 0" log_all.txt
check_contains "Termination 1" log_all.txt
check_success slab_term0.fdf
check_success slab_term1.fdf


# --- 4. Direct --termination selection ---
echo -e "\n--- Testing direct --termination selection ---"
rm -f slab.fdf
stb-slab -f structure.fdf --hkl 1 0 0 --termination 1 --no-intro > log_term1.txt 2>&1
check_contains "Termination 1" log_term1.txt
check_success slab.fdf


# --- 5. Round-trip: output is valid, 2D-detected input for another tool ---
echo -e "\n--- Verifying the generated slab is valid, vacuum-padded input for stb-kgrid ---"
stb-kgrid --file slab.fdf --type fdf --density 0.2 --no-intro > log_roundtrip.txt 2>&1
check_contains "System appears to be 2D" log_roundtrip.txt


# --- 6. Two-species fixture: --symmetrize changes the termination set ---
echo -e "\n--- Testing rocksalt.fdf (hkl 1 1 1), with and without --symmetrize ---"
rm -f slab.fdf sym_term0.fdf sym_term1.fdf
stb-slab -f rocksalt.fdf --hkl 1 1 1 --no-intro > log_rocksalt.txt 2>&1
check_contains "Terminations found:.*1" log_rocksalt.txt
check_success slab.fdf

stb-slab -f rocksalt.fdf --hkl 1 1 1 --symmetrize --all -o sym.fdf --no-intro > log_rocksalt_sym.txt 2>&1
check_contains "Terminations found:.*2" log_rocksalt_sym.txt
check_contains "Symmetric.*Yes" log_rocksalt_sym.txt
check_success sym_term0.fdf
check_success sym_term1.fdf


# --- 7. Interactive termination selection (-i) ---
echo -e "\n--- Testing -i/--interactive termination selection ---"
rm -f slab.fdf
printf '0\n' | stb-slab -f structure.fdf --hkl 1 0 0 -i --no-intro > log_interactive_sel.txt 2>&1
check_contains "Select a termination ID" log_interactive_sel.txt
check_success slab.fdf


# --- 8. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: missing structure file"
stb-slab -f does_not_exist.fdf --hkl 1 0 0 --no-intro > log_missing.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing.txt

echo "Testing: invalid Miller index (0, 0, 0)"
stb-slab -f structure.fdf --hkl 0 0 0 --no-intro > log_zero_hkl.txt 2>&1
check_exit_code $? 1
check_contains "not valid" log_zero_hkl.txt

echo "Testing: --termination out of range"
stb-slab -f structure.fdf --hkl 1 0 0 --termination 99 --no-intro > log_oob.txt 2>&1
check_exit_code $? 1
check_contains "out of range" log_oob.txt

echo "Testing: -f/--file missing (required)"
stb-slab --hkl 1 0 0 --no-intro > log_missing_file_arg.txt 2>&1
check_exit_code $? 2

echo "Testing: --hkl missing (required)"
stb-slab -f structure.fdf --no-intro > log_missing_hkl_arg.txt 2>&1
check_exit_code $? 2

echo "Testing: -i and --termination are mutually exclusive"
stb-slab -f structure.fdf --hkl 1 0 0 -i --termination 0 --no-intro > log_mutex.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-slab --version > log_version.txt 2>&1
check_contains "stb-slab" log_version.txt

echo "Testing: --help documents the Miller index, vacuum, and passivate options"
stb-slab --help > log_help.txt 2>&1
check_contains "Miller index" log_help.txt
check_contains "vacuum" log_help.txt
check_contains "passivate" log_help.txt


# --- 9. Interactive path (stb-suite, shortcut 1.6) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 1.6) ---"

echo "Testing: navigate 1.6 -> invalid file then valid -> hkl '1 0 0' -> defaults -> mode 1 -> default output -> quit"
rm -f slab.fdf
printf '1.6\ndoes_not_exist.fdf\nstructure.fdf\n1 0 0\n\n\nn\nn\nn\n1\nn\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_contains "File not found" log_menu.txt
check_contains "Slab Summary" log_menu.txt
check_success slab.fdf


popd > /dev/null

# --- 10. Summary ---
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
