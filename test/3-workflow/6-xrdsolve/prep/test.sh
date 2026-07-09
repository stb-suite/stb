#!/bin/bash

# --- Setup ---
# Smoke test for stb-xrdsearch (Structure Solution XRD, Prep, item 3.6.1)
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
echo "--- Starting tester for STB-XRDSearch (item 3.6.1) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Basic run across a few space groups ---
echo -e "\n--- Testing a basic run (--species Na Cl --num-ions 4 4 --groups 225,221,62) ---"

stb-xrdsearch --species Na Cl --num-ions 4 4 --groups 225,221,62 --seed 1 \
    -o candidates --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "Space groups:.*225, 221, 62" log_basic.txt
check_contains "3 candidate folder(s) written" log_basic.txt
check_success candidates/group_225/structure.fdf
check_success candidates/group_221/structure.fdf
check_success candidates/group_62/structure.fdf
check_contains "Next:.*stb-xrdrank" log_basic.txt


# --- 3. --count-per-group ---
echo -e "\n--- Testing --count-per-group 2 ---"

stb-xrdsearch --species Na Cl --num-ions 4 4 --groups 225 --count-per-group 2 \
    -o multi --no-intro > log_multi.txt 2>&1
check_exit_code $? 0
check_success multi/group_225_1/structure.fdf
check_success multi/group_225_2/structure.fdf


# --- 4. Space-separated --groups ---
echo -e "\n--- Testing space-separated --groups ---"

stb-xrdsearch --species Na Cl --num-ions 4 4 --groups "225 62" -o spaced --no-intro > log_spaced.txt 2>&1
check_exit_code $? 0
check_success spaced/group_225/structure.fdf
check_success spaced/group_62/structure.fdf

echo "Testing: duplicate --groups entries are deduped (not silently overwritten)"
stb-xrdsearch --species Na Cl --num-ions 4 4 --groups 225,227,225 -o dup --no-intro > log_dup.txt 2>&1
check_contains "Warning:.*--groups had duplicate entries" log_dup.txt
check_contains "Space groups:.*225, 227" log_dup.txt


# --- 5. --ml-rank / --keep-top-per-group ---
echo -e "\n--- Testing --ml-rank --keep-top-per-group 1 (slow: MACE relaxes each candidate) ---"

timeout 180 stb-xrdsearch --species Na Cl --num-ions 4 4 --groups 225 --count-per-group 3 \
    --ml-rank --keep-top-per-group 1 --seed 7 -o ranked --no-intro > log_mlrank.txt 2>&1
check_exit_code $? 0
check_contains "1 folder(s) (2 of 3 skipped/dropped) (group_225_rank1)" log_mlrank.txt
check_success ranked/group_225_rank1/structure.fdf

echo "Testing: --keep-top-per-group without --ml-rank is rejected"
stb-xrdsearch --species Na Cl --num-ions 4 4 --groups 225 --keep-top-per-group 1 \
    -o err_keeptop --no-intro > log_keeptop_noml.txt 2>&1
check_exit_code $? 2
check_contains "requires --ml-rank" log_keeptop_noml.txt


# --- 6. Refuses to reuse a non-empty output folder ---
echo -e "\n--- Testing refusal to overwrite an existing candidate folder ---"

stb-xrdsearch --species Na Cl --num-ions 4 4 --groups 225 -o candidates --no-intro > log_reuse.txt 2>&1
check_exit_code $? 1
check_contains "already has content" log_reuse.txt


# --- 7. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: --species/--num-ions length mismatch"
stb-xrdsearch --species Na Cl --num-ions 4 --groups 225 -o err1 --no-intro > log_mismatch.txt 2>&1
check_exit_code $? 2
check_contains "must match one-to-one" log_mismatch.txt

echo "Testing: --count-per-group 0 is rejected"
stb-xrdsearch --species Na Cl --num-ions 4 4 --groups 225 --count-per-group 0 -o err2 --no-intro > log_zero.txt 2>&1
check_exit_code $? 2
check_contains "must be at least 1" log_zero.txt

echo "Testing: non-numeric --groups"
stb-xrdsearch --species Na Cl --num-ions 4 4 --groups "not-a-number" -o err3 --no-intro > log_badgroups.txt 2>&1
check_exit_code $? 2
check_contains "must be a list of space group numbers" log_badgroups.txt

echo "Testing: every requested group fails (invalid space group number)"
stb-xrdsearch --species Na Cl --num-ions 4 4 --groups 99999 -o err4 --no-intro > log_allfail.txt 2>&1
check_exit_code $? 1
check_contains "no candidates were generated" log_allfail.txt

echo "Testing: --species omitted (required)"
stb-xrdsearch --num-ions 4 4 --groups 225 --no-intro > log_missing_species.txt 2>&1
check_exit_code $? 2

echo "Testing: --groups omitted (required)"
stb-xrdsearch --species Na Cl --num-ions 4 4 --no-intro > log_missing_groups.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-xrdsearch --version > log_version.txt 2>&1
check_contains "stb-xrdsearch" log_version.txt

echo "Testing: --help documents --species, --groups, --count-per-group"
stb-xrdsearch --help > log_help.txt 2>&1
check_contains "species" log_help.txt
check_contains "groups" log_help.txt
check_contains "count-per-group" log_help.txt


# --- 8. Interactive path (stb-suite, shortcut 3.6.1) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 3.6.1) ---"

echo "Testing: navigate 3.6.1 -> species 'Na 4', 'Cl 4', blank to finish -> groups '225,139' -> count 1 -> output -> quit"
rm -rf menu_candidates
printf '3.6.1\nNa 4\nCl 4\n\n225,139\n\nmenu_candidates\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_contains "2 candidate folder(s) written" log_menu.txt
check_success menu_candidates/group_225/structure.fdf
check_success menu_candidates/group_139/structure.fdf


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
