#!/bin/bash

# --- Setup ---
# Smoke test for stb-dftu (DFT+U / Hubbard Block Generator, item 1.4)
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

check_success() {
    if [ -s "$1" ]; then
        echo -e " ... ${GREEN}OK${NC} (file '$1' created)"
        PASS=$((PASS+1))
    else
        echo -e " ... ${RED}FAIL${NC} (file '$1' was not created)"
        FAIL=$((FAIL+1))
    fi
}

check_contains() {
    if grep -q -- "$1" "$2" 2>/dev/null; then
        echo -e "   -> ${GREEN}Verified:${NC} '$1' found in '$2'"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} '$1' NOT found in '$2'"
        FAIL=$((FAIL+1))
    fi
}

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
echo "--- Starting tester for STB-DFTU (item 1.4) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Single species (auto-detected shell) ---
echo -e "\n--- Testing a single species with auto-detected shell ---"
stb-dftu --species Mn --u 3.9 --no-intro > log_single.txt 2>&1
check_exit_code $? 0
check_contains "shell=3d (n=3, l=2)" log_single.txt
check_contains "U=3.900 eV, J=0.000 eV" log_single.txt
check_contains "%block LDAU.proj" log_single.txt
check_contains "Mn   1" log_single.txt
check_contains "n=3    2" log_single.txt
check_contains "3.900    0.000" log_single.txt
check_contains "Literature reference for Mn: 3.90 eV" log_single.txt


# --- 3. Multiple species, explicit shell and J, saved to file ---
echo -e "\n--- Testing multiple species with explicit --shell/--j and -o ---"
stb-dftu --species Fe Co --u 5.3 3.32 --j 0.5 0.5 --shell 3d 3d -o multi_LDAU.fdf --no-intro > log_multi.txt 2>&1
check_exit_code $? 0
check_success multi_LDAU.fdf
check_contains "Fe   1" multi_LDAU.fdf
check_contains "Co   1" multi_LDAU.fdf
check_contains "5.300    0.500" multi_LDAU.fdf
check_contains "3.320    0.500" multi_LDAU.fdf

echo "Testing: exactly one %block/%endblock wraps both species stanzas"
BLOCK_COUNT=$(grep -c "%block LDAU.proj" multi_LDAU.fdf)
ENDBLOCK_COUNT=$(grep -c "%endblock LDAU.proj" multi_LDAU.fdf)
if [ "$BLOCK_COUNT" -eq 1 ] && [ "$ENDBLOCK_COUNT" -eq 1 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} single %block/%endblock pair wraps both stanzas"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} expected exactly one %block/%endblock pair"
    FAIL=$((FAIL+1))
fi


# --- 4. --list-reference and --suggest ---
echo -e "\n--- Testing --list-reference ---"
stb-dftu --list-reference --no-intro > log_list.txt 2>&1
check_exit_code $? 0
check_contains "Wang/Maxisch/Ceder" log_list.txt
check_contains "Ni  6.20 eV" log_list.txt
check_contains "starting-point sanity check" log_list.txt

echo "Testing: --list-reference does not require --species/--u"
if grep -q "%block" log_list.txt; then
    echo -e "   -> ${RED}Failed:${NC} --list-reference should not print a block"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} --list-reference prints no block"
    PASS=$((PASS+1))
fi

echo -e "\n--- Testing --suggest ---"
stb-dftu --suggest Ni --no-intro > log_suggest.txt 2>&1
check_exit_code $? 0
check_contains "6.20 eV" log_suggest.txt
check_contains "not a validation" log_suggest.txt

echo "Testing: --suggest for an element with no reference value is a clean error"
stb-dftu --suggest Xx --no-intro > log_suggest_bad.txt 2>&1
check_exit_code $? 1
check_contains "No reference value tabulated" log_suggest_bad.txt


# --- 5. Error and robustness cases ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: --u never used silently -- --species alone is rejected"
stb-dftu --species Mn --no-intro > log_no_u.txt 2>&1
check_exit_code $? 2

echo "Testing: mismatched --species/--u length"
stb-dftu --species Mn Fe --u 3.9 --no-intro > log_mismatch_u.txt 2>&1
check_exit_code $? 1
check_contains "must match one-to-one" log_mismatch_u.txt

echo "Testing: mismatched --species/--j length"
stb-dftu --species Mn Fe --u 3.9 5.3 --j 0.5 --no-intro > log_mismatch_j.txt 2>&1
check_exit_code $? 1
check_contains "must match one-to-one" log_mismatch_j.txt

echo "Testing: mismatched --species/--shell length"
stb-dftu --species Mn Fe --u 3.9 5.3 --shell 3d --no-intro > log_mismatch_shell.txt 2>&1
check_exit_code $? 1
check_contains "must match one-to-one" log_mismatch_shell.txt

echo "Testing: species with no default shell and no --shell given"
stb-dftu --species Si --u 2.0 --no-intro > log_no_shell.txt 2>&1
check_exit_code $? 1
check_contains "No default correlated shell known" log_no_shell.txt

echo "Testing: --shell explicitly supplied for a species with no default works"
stb-dftu --species Si --u 2.0 --shell 3d --no-intro > log_shell_override.txt 2>&1
check_exit_code $? 0

echo "Testing: --version"
stb-dftu --version > log_version.txt 2>&1
check_contains "stb-dftu" log_version.txt

echo "Testing: --help documents species/u/j/shell/list-reference/suggest"
stb-dftu --help > log_help.txt 2>&1
check_contains "species" log_help.txt
check_contains "list-reference" log_help.txt
check_contains "suggest" log_help.txt
check_contains "linear-response" log_help.txt


# --- 6. Interactive path (stb-suite, shortcut 1.4) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 1.4) ---"

echo "Testing: navigate 1.4 -> generate (not list) -> Mn -> U=3.9 -> J default -> auto shell -> finish -> no save -> quit"
printf '1.4\nn\nMn\n3.9\n0.0\n\n\n\n\n0\n' | timeout 20 stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_contains "shell=3d (n=3, l=2)" log_menu.txt
check_contains "%block LDAU.proj" log_menu.txt

echo "Testing: navigate 1.4 -> list-reference path"
printf '1.4\ny\n\n0\n' | timeout 20 stb-suite > log_menu_list.txt 2>&1
check_exit_code $? 0
check_contains "Wang/Maxisch/Ceder" log_menu_list.txt


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
