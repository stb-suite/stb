#!/bin/bash

# --- Setup ---
# Smoke test for stb-nanotube (Nanotube/Nanoribbon Builder, item 2.4)
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
echo "--- Starting tester for STB-Nanotube (item 2.4) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/graphene.fdf" "$TEST_DIR/"
cp "$(cd "$FIXTURE_DIR/../3-k-path" && pwd)/example_3D_silicon.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Tube mode: (6, 0) zigzag-style tube from graphene ---
echo -e "\n--- Testing tube mode (graphene.fdf, chirality 6 0) ---"
rm -f nanotube.fdf
stb-nanotube -f graphene.fdf --chirality 6 0 --no-intro > log_tube.txt 2>&1
check_contains "Translation vector (t1, t2):.*(0, 1)" log_tube.txt
check_contains "Cells in periodic unit (N_cells):.*6" log_tube.txt
check_contains "Tube diameter (2R):.*4.6983 Ang" log_tube.txt
check_contains "Output atoms:.*12" log_tube.txt
check_success nanotube.fdf
check_contains "NumberofAtoms      12" nanotube.fdf


# --- 3. Round-trip: tube is a valid, 1D-detected input ---
echo -e "\n--- Verifying the generated tube is valid, vacuum-padded (1D) input for stb-kgrid ---"
stb-kgrid --file nanotube.fdf --type fdf --density 0.2 --no-intro > log_tube_roundtrip.txt 2>&1
check_contains "System appears to be 1D" log_tube_roundtrip.txt
check_contains "Suggested Monkhorst-Pack grid: 1 1 13" log_tube_roundtrip.txt


# --- 4. Ribbon mode: finite width scales with --repeats ---
echo -e "\n--- Testing ribbon mode (graphene.fdf, chirality 6 0, repeats 4) ---"
rm -f ribbon.fdf
stb-nanotube -f graphene.fdf --chirality 6 0 --mode ribbon --repeats 4 -o ribbon.fdf --no-intro > log_ribbon.txt 2>&1
check_contains "Ribbon width:.*49.7099 Ang" log_ribbon.txt
check_contains "Output atoms:.*48" log_ribbon.txt
check_success ribbon.fdf

echo -e "\n--- Verifying the generated ribbon is also a valid 1D input for stb-kgrid ---"
stb-kgrid --file ribbon.fdf --type fdf --density 0.2 --no-intro > log_ribbon_roundtrip.txt 2>&1
check_contains "System appears to be 1D" log_ribbon_roundtrip.txt


# --- 5. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: chirality (0, 0)"
stb-nanotube -f graphene.fdf --chirality 0 0 --no-intro > log_zero_chir.txt 2>&1
check_exit_code $? 1
check_contains "not valid" log_zero_chir.txt

echo "Testing: --repeats 0"
stb-nanotube -f graphene.fdf --chirality 6 0 --repeats 0 --no-intro > log_zero_repeats.txt 2>&1
check_exit_code $? 1
check_contains "repeats must be" log_zero_repeats.txt

echo "Testing: non-2D input (bulk silicon)"
stb-nanotube -f example_3D_silicon.fdf --chirality 6 0 --no-intro > log_non2d.txt 2>&1
check_exit_code $? 1
check_contains "must be a 2D monolayer" log_non2d.txt

echo "Testing: missing structure file"
stb-nanotube -f does_not_exist.fdf --chirality 6 0 --no-intro > log_missing.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing.txt

echo "Testing: -f/--file missing (required)"
stb-nanotube --chirality 6 0 --no-intro > log_missing_file_arg.txt 2>&1
check_exit_code $? 2

echo "Testing: --chirality missing (required)"
stb-nanotube -f graphene.fdf --no-intro > log_missing_chir_arg.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-nanotube --version > log_version.txt 2>&1
check_contains "stb-nanotube" log_version.txt

echo "Testing: --help documents chirality and vacuum options"
stb-nanotube --help > log_help.txt 2>&1
check_contains "Chirality" log_help.txt
check_contains "vacuum" log_help.txt


# --- 6. Interactive path (stb-suite, shortcut 1.7) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 1.7) ---"

echo "Testing: navigate 1.7 -> invalid file then valid -> chirality '6 0' -> mode 1 (tube) -> defaults -> quit"
rm -f nanotube.fdf
printf '2.4\ndoes_not_exist.fdf\ngraphene.fdf\n6 0\n1\n\n\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_contains "File not found" log_menu.txt
check_contains "Success" log_menu.txt
check_success nanotube.fdf


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
