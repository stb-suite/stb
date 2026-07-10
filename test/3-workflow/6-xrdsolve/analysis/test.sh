#!/bin/bash

# --- Setup ---
# Smoke test for stb-xrdrank (Structure Solution XRD, Analysis, item 4.6.2)
# NOTE: pyxtal's Similarity() has no fast mode (~10-30s per candidate), so this
# test intentionally uses only 2-3 tiny candidates per --input-dir to keep the
# total runtime reasonable. Generous timeouts below (not a hang, just slow).
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
echo "--- Starting tester for STB-XRDRank (item 4.6.2) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/experimental.dat" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null

# Generate a tiny candidate set with stb-xrdsearch itself (fast: no Similarity()
# calls happen until stb-xrdrank runs below).
stb-xrdsearch --species Na Cl --num-ions 4 4 --groups 225,62 --seed 3 \
    -o candidates --no-intro > /dev/null 2>&1


# --- 2. Basic run (this is the slow part: 2 candidates x ~10-30s each) ---
echo -e "\n--- Testing a basic run (--input-dir candidates --experimental experimental.dat) ---"

timeout 120 stb-xrdrank --input-dir candidates --experimental experimental.dat \
    --no-intro -o rank.txt > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "Candidates found:.*2 in 'candidates'" log_basic.txt
check_contains "Ranking (2 of 2 shown, best match first):" log_basic.txt
check_contains "group_225" log_basic.txt
check_contains "group_62" log_basic.txt
check_success rank.txt
check_contains "group_225" rank.txt


# --- 3. --top ---
echo -e "\n--- Testing --top 1 ---"

timeout 120 stb-xrdrank --input-dir candidates --experimental experimental.dat \
    --top 1 --no-intro -o rank_top.txt > log_top.txt 2>&1
check_exit_code $? 0
check_contains "Ranking (1 of 2 shown, best match first):" log_top.txt
check_success rank_top.txt


# --- 4. Prefers a relaxed .STRUCT_OUT over the raw structure.fdf ---
echo -e "\n--- Testing .STRUCT_OUT preference (simulating a candidate SIESTA has already relaxed) ---"

cp "$FIXTURE_DIR/../../../2-analysis/8-xrd/siesta.STRUCT_OUT" candidates/group_225/
timeout 120 stb-xrdrank --input-dir candidates --experimental experimental.dat \
    --no-intro -o rank_relaxed.txt > log_relaxed.txt 2>&1
check_exit_code $? 0
check_contains "group_225: similarity.*relaxed" log_relaxed.txt
check_contains "group_62: similarity.*raw" log_relaxed.txt
check_contains "group_225 .*relaxed" rank_relaxed.txt
check_contains "group_62 .*raw" rank_relaxed.txt
rm candidates/group_225/siesta.STRUCT_OUT


# --- 5. Skips a flat .fdf that collides with a subfolder candidate's name ---
echo -e "\n--- Testing name-collision skip (flat .fdf with the same name as a subfolder) ---"

mkdir -p dup_candidates/group_225
cp candidates/group_225/structure.fdf dup_candidates/group_225/structure.fdf
cp candidates/group_225/structure.fdf dup_candidates/group_225.fdf
timeout 60 stb-xrdrank --input-dir dup_candidates --experimental experimental.dat \
    --no-intro -o rank_dup.txt > log_dup.txt 2>&1
check_exit_code $? 0
check_contains "Warning: skipping.*group_225.*collides" log_dup.txt
check_contains "Candidates found:.*1 in 'dup_candidates'" log_dup.txt


# --- 6. Error and robustness cases (all fast: fail before any Similarity() call) ---
echo -e "\n--- Testing error cases ---"

echo "Testing: nonexistent --input-dir"
stb-xrdrank --input-dir does_not_exist --experimental experimental.dat --no-intro > log_missing_dir.txt 2>&1
check_exit_code $? 1
check_contains "is not a directory" log_missing_dir.txt

echo "Testing: --input-dir with no candidate structures"
mkdir -p empty_dir
stb-xrdrank --input-dir empty_dir --experimental experimental.dat --no-intro > log_empty_dir.txt 2>&1
check_exit_code $? 1
check_contains "no candidate structures found" log_empty_dir.txt

echo "Testing: nonexistent --experimental file"
stb-xrdrank --input-dir candidates --experimental does_not_exist.dat --no-intro > log_missing_exp.txt 2>&1
check_exit_code $? 1

echo "Testing: --experimental with too few points"
printf "10.0 1.0\n20.0 2.0\n" > tiny_exp.dat
stb-xrdrank --input-dir candidates --experimental tiny_exp.dat --no-intro > log_tiny_exp.txt 2>&1
check_exit_code $? 1
check_contains "at least 4 are needed" log_tiny_exp.txt

echo "Testing: --input-dir omitted (required)"
stb-xrdrank --experimental experimental.dat --no-intro > log_missing_arg_dir.txt 2>&1
check_exit_code $? 2

echo "Testing: --experimental omitted (required)"
stb-xrdrank --input-dir candidates --no-intro > log_missing_arg_exp.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-xrdrank --version > log_version.txt 2>&1
check_contains "stb-xrdrank" log_version.txt

echo "Testing: --help documents --input-dir, --experimental, --top"
stb-xrdrank --help > log_help.txt 2>&1
check_contains "input-dir" log_help.txt
check_contains "experimental" log_help.txt
check_contains "top" log_help.txt


# --- 7. Interactive path (stb-suite, shortcut 3.6.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 3.6.2) (slow: real Similarity() calls) ---"

echo "Testing: navigate 3.6.2 -> candidates dir -> experimental.dat -> default wavelength -> no top -> output -> quit"
printf '4.6.2\ncandidates\nexperimental.dat\n\n\nmenu_rank.txt\n\n0\n' | timeout 120 stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_contains "Ranking (2 of 2 shown, best match first):" log_menu.txt
check_success menu_rank.txt


popd > /dev/null

# --- 8. Summary ---
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
