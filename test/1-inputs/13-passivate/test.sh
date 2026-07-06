#!/bin/bash

# --- Setup ---
# Smoke test for stb-passivate (Surface Passivator, item 1.13)
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
echo "--- Starting tester for STB-Passivate (item 1.13) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/si111_slab.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/si111_slab_bad_termination.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Happy path: single-missing-bond termination ---
echo -e "\n--- Testing default passivation (si111_slab.fdf, 8 deficit-1 sites) ---"
rm -f passivated.fdf
stb-passivate -f si111_slab.fdf -o passivated.fdf --no-intro > log_happy.txt 2>&1
check_contains "Dangling bonds found:.*8" log_happy.txt
check_contains "Auto-passivated:.*8" log_happy.txt
check_contains "Output atoms:.*40" log_happy.txt
check_success passivated.fdf
check_contains "NumberofAtoms      40" passivated.fdf
check_contains " 2   1   H" passivated.fdf

echo "Testing: 'needs manual attention' warning is NOT printed on the good termination"
if grep -q "manual attention" log_happy.txt; then
    echo -e "   -> ${RED}Failed:${NC} unexpected manual-attention warning on the deficit-1 fixture"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no spurious manual-attention warning"
    PASS=$((PASS+1))
fi


# --- 3. Bad termination: multi-bond deficit, unresolved warning ---
echo -e "\n--- Testing the deficit-3 termination (needs-manual-attention path) ---"
rm -f passivated_bad.fdf
stb-passivate -f si111_slab_bad_termination.fdf -o passivated_bad.fdf --no-intro \
    > log_bad_termination.txt 2>&1
check_contains "Dangling bonds found:.*8" log_bad_termination.txt
check_contains "Auto-passivated:.*0" log_bad_termination.txt
check_contains "missing 2+ bonds" log_bad_termination.txt
check_contains "Output atoms:.*32" log_bad_termination.txt
check_success passivated_bad.fdf
check_contains "NumberofAtoms      32" passivated_bad.fdf


# --- 4. --bond-length / --cutoff overrides ---
echo -e "\n--- Testing --bond-length override ---"
rm -f custom.fdf
stb-passivate -f si111_slab.fdf --bond-length 1.6 -o custom.fdf --no-intro > log_bondlength.txt 2>&1
check_success custom.fdf
python3 - <<'EOF'
from stb.core import structure_io
s = structure_io.to_pymatgen(structure_io.read_fdf("custom.fdf"))
h = next(site for site in s if site.specie.symbol == "H")
si_sites = [site for site in s if site.specie.symbol == "Si"]
nearest = min(h.distance(si) for si in si_sites)
assert abs(nearest - 1.6) < 1e-3, f"expected ~1.6 Ang, got {nearest}"
print("OK: nearest Si-H distance is", round(nearest, 4))
EOF
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} --bond-length 1.6 honored in the output geometry"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} --bond-length override not reflected in output geometry"
    FAIL=$((FAIL+1))
fi


# --- 5. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: invalid --passivant element"
stb-passivate -f si111_slab.fdf --passivant Xx --no-intro > log_bad_element.txt 2>&1
check_exit_code $? 2
check_contains "not a valid Element" log_bad_element.txt

echo "Testing: missing structure file"
stb-passivate -f does_not_exist.fdf --no-intro > log_missing.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing.txt

echo "Testing: -f/--file missing (required)"
stb-passivate --no-intro > log_missing_file_arg.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-passivate --version > log_version.txt 2>&1
check_contains "stb-passivate" log_version.txt

echo "Testing: --help documents --passivant/--cutoff/--bond-length"
stb-passivate --help > log_help.txt 2>&1
check_contains "passivant" log_help.txt
check_contains "cutoff" log_help.txt
check_contains "bond-length" log_help.txt


# --- 6. Interactive path (stb-suite, shortcut 1.13) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 1.13) ---"

echo "Testing: navigate 1.13 -> invalid file then valid -> default passivant/cutoff/bond-length -> default output -> quit"
rm -f passivated.fdf
printf '1.13\ndoes_not_exist.fdf\nsi111_slab.fdf\n\n\n\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_contains "File not found" log_menu.txt
check_contains "Auto-passivated:.*8" log_menu.txt
check_success passivated.fdf


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
