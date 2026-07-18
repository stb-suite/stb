#!/bin/bash

# --- Setup ---
# Smoke test for stb-clean (Clean Files Tool, item 5.2).
#
# Uses its own disposable test_files/ subdir with small dummy files -- NOT
# the real ~12MB completed-SIESTA-run fixture already sitting next to this
# script (0_NORMAL_EXIT, Sn3O4.*, etc.), which is kept around as a realistic
# "messy calc folder" for manual inspection, not to be deleted by an
# automated test run.
#
# stb-clean's purpose (confirmed with the tool's author): prepare a SIESTA
# calc folder to restart a fresh run after a previous one that may or may not
# have finished OK -- delete every stale output/intermediate file, keep only
# what's needed to resubmit (pseudopotentials, .fdf, the submission script).
# This is deliberately the OPPOSITE of an analysis-results archiver, so the
# checks below assert that outputs (.out, .XV, .DM, ...) ARE removed by
# default, not preserved. Every removal is archived (tar.gz) first by
# default -- there's no --dry-run, the archive is the safety net.
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$FIXTURE_DIR/test_files"

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

PASS=0
FAIL=0

# Checks that file $1 exists.
check_exists() {
    if [ -f "$1" ]; then
        echo -e " ... ${GREEN}OK${NC} ('$1' still present)"
        PASS=$((PASS+1))
    else
        echo -e " ... ${RED}FAIL${NC} ('$1' missing, expected it to be kept)"
        FAIL=$((FAIL+1))
    fi
}

# Checks that file $1 does NOT exist.
check_removed() {
    if [ ! -f "$1" ]; then
        echo -e " ... ${GREEN}OK${NC} ('$1' removed)"
        PASS=$((PASS+1))
    else
        echo -e " ... ${RED}FAIL${NC} ('$1' still present, expected it to be removed)"
        FAIL=$((FAIL+1))
    fi
}

# Checks that $2 (a file) contains (grep -q) pattern $1.
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

# Checks that exactly one *.tar.gz backup archive exists directly under $1.
check_archive_created() {
    if compgen -G "$1"/*.tar.gz > /dev/null; then
        echo -e " ... ${GREEN}OK${NC} (backup archive created in '$1')"
        PASS=$((PASS+1))
    else
        echo -e " ... ${RED}FAIL${NC} (no backup archive found in '$1')"
        FAIL=$((FAIL+1))
    fi
}

# Populates $TEST_DIR (or $1, if given) with a small stand-in for a completed
# (or crashed) SIESTA run: a few "keep" files (input/submission) and a few
# "delete" ones (stale outputs/intermediates) -- named after real SIESTA
# file kinds, but tiny, throwaway content. rm -f's the *contents*, never
# rm -rf's the directory itself -- some calls below run while already
# pushd'd inside it, and removing the cwd out from under the shell breaks
# every command after it.
populate_fixture() {
    local dir="${1:-$TEST_DIR}"
    mkdir -p "$dir"
    rm -f "$dir"/*
    echo "SystemName test" > "$dir/calc.fdf"
    echo "dummy psf" > "$dir/Si.psf"
    echo "#!/bin/bash" > "$dir/run.sh"
    echo "dummy out" > "$dir/calc.out"
    echo "dummy xv" > "$dir/calc.XV"
    echo "dummy dm" > "$dir/calc.DM"
    echo "normal exit" > "$dir/0_NORMAL_EXIT"
}


# --- 1. Preparation ---
echo "--- Starting tester for STB-Clean (item 5.2) ---"
populate_fixture
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. --no-confirm: default --keep removes outputs, keeps inputs, backs up first ---
echo -e "\n--- Testing --no-confirm (default --keep .psml .psf .fdf .sh, archive on by default) ---"
stb-clean --no-intro --no-confirm > "$FIXTURE_DIR/log_default.txt" 2>&1
check_exists calc.fdf
check_exists Si.psf
check_exists run.sh
check_removed calc.out
check_removed calc.XV
check_removed calc.DM
check_removed 0_NORMAL_EXIT
check_archive_created "$FIXTURE_DIR"
rm -f "$FIXTURE_DIR"/*.tar.gz


# --- 3. Custom --keep + --no-archive ---
echo -e "\n--- Testing custom --keep (only .fdf) with --no-archive ---"
populate_fixture
stb-clean --no-intro --no-confirm --no-archive --keep .fdf > "$FIXTURE_DIR/log_custom.txt" 2>&1
check_exists calc.fdf
check_removed Si.psf
check_removed run.sh
check_removed calc.out
if compgen -G "$FIXTURE_DIR"/*.tar.gz > /dev/null; then
    echo -e "   -> ${RED}Failed:${NC} a backup archive was created despite --no-archive"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no backup archive created (--no-archive)"
    PASS=$((PASS+1))
fi


# --- 4. --keep validation (bare extension without a leading dot) ---
echo -e "\n--- Testing --keep validation (rejects 'fdf' without a leading dot) ---"
stb-clean --no-intro --path . --keep fdf > "$FIXTURE_DIR/log_badkeep.txt" 2>&1
check_exit_code $? 2
check_contains "must start with" "$FIXTURE_DIR/log_badkeep.txt"


# --- 5. --recursive ---
echo -e "\n--- Testing --recursive (cleans subdirectories too) ---"
populate_fixture
mkdir -p sub
echo "nested dm" > sub/calc.DM
echo "Testing: non-recursive leaves the subdirectory's file alone"
stb-clean --no-intro --no-confirm --no-archive > "$FIXTURE_DIR/log_norecursive.txt" 2>&1
check_exists sub/calc.DM
echo "Testing: --recursive removes it"
stb-clean --no-intro --no-confirm --no-archive --recursive > "$FIXTURE_DIR/log_recursive.txt" 2>&1
check_removed sub/calc.DM
rmdir sub 2>/dev/null


# --- 6. Batch clean via a glob pattern matching several directories ---
echo -e "\n--- Testing batch mode (--path matching multiple directories) ---"
popd > /dev/null
rm -rf "$FIXTURE_DIR/batch_a" "$FIXTURE_DIR/batch_b"
populate_fixture "$FIXTURE_DIR/batch_a"
populate_fixture "$FIXTURE_DIR/batch_b"
pushd "$FIXTURE_DIR" > /dev/null
stb-clean --no-intro --no-confirm --no-archive --path "batch_*" > "$FIXTURE_DIR/log_batch.txt" 2>&1
check_contains "2 director" "$FIXTURE_DIR/log_batch.txt"
check_exists batch_a/calc.fdf
check_removed batch_a/calc.out
check_exists batch_b/calc.fdf
check_removed batch_b/calc.out
rm -rf "$FIXTURE_DIR/batch_a" "$FIXTURE_DIR/batch_b"
popd > /dev/null
pushd "$TEST_DIR" > /dev/null


# --- 7. --save-report ---
echo -e "\n--- Testing --save-report ---"
populate_fixture
rm -f stb_clean_report.txt
stb-clean --no-intro --no-confirm --no-archive --save-report > "$FIXTURE_DIR/log_savereport.txt" 2>&1
check_exists stb_clean_report.txt
check_contains "STB-CLEAN REPORT" stb_clean_report.txt
rm -f stb_clean_report.txt


# --- 8. Interactive path (stb-suite, shortcut 5.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 5.2) ---"
populate_fixture
# Prompts in order: path, extensions (default), recursive?, archive?, no-confirm?, save-report?
printf '5.2\n%s\n\nn\nn\ny\nn\n\n0\n' "$TEST_DIR" | stb-suite > "$FIXTURE_DIR/log_interactive.txt" 2>&1
check_contains "Cleanup complete" "$FIXTURE_DIR/log_interactive.txt"
check_exists calc.fdf
check_removed calc.out


popd > /dev/null

# Test-runner logs, not stb-clean output to inspect -- always removed,
# regardless of whether test_files/ itself is kept below. Written outside
# test_files/ in the first place (see the runs above) since anything inside
# it is fair game for stb-clean itself to delete (a .txt extension isn't in
# any of the --keep lists exercised here).
rm -f "$FIXTURE_DIR"/log_*.txt "$FIXTURE_DIR"/*.tar.gz

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
