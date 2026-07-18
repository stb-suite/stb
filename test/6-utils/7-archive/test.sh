#!/bin/bash

# --- Setup ---
# Smoke test for stb-archive (Archive Packager, item 6.7).
#
# Fixtures are real files from a completed single-point Sn3O4 calculation
# (shared with the other 5-utils fixtures): calc.fdf/calc.out (SystemLabel
# Sn3O4, even though the log isn't named Sn3O4.out -- exercises
# core.siesta_log.find_out_file's generic-filename fallback), O.psf/Sn.psf,
# Sn3O4.bands/.DOS/.XV (expected to be archived by default), plus Sn3O4.RHO
# and MESSAGES (expected to be EXCLUDED by default -- a large regenerable
# intermediate and a no-extension log file, respectively).
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$FIXTURE_DIR/test_files"

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

PASS=0
FAIL=0

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
    if grep -q "$1" "$2" 2>/dev/null; then
        echo -e "   -> ${GREEN}Verified:${NC} '$1' found in '$2'"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} '$1' NOT found in '$2'"
        FAIL=$((FAIL+1))
    fi
}

check_archive_contains() {
    if tar tzf "$1" 2>/dev/null | grep -qx "$2"; then
        echo -e "   -> ${GREEN}Verified:${NC} '$2' present in '$1'"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} '$2' NOT present in '$1'"
        FAIL=$((FAIL+1))
    fi
}

check_archive_not_contains() {
    if tar tzf "$1" 2>/dev/null | grep -qx "$2"; then
        echo -e "   -> ${RED}Failed:${NC} '$2' unexpectedly present in '$1'"
        FAIL=$((FAIL+1))
    else
        echo -e "   -> ${GREEN}Verified:${NC} '$2' NOT present in '$1'"
        PASS=$((PASS+1))
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
echo "--- Starting tester for STB-ARCHIVE (item 6.7) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR"/calc.fdf "$FIXTURE_DIR"/calc.out "$FIXTURE_DIR"/O.psf "$FIXTURE_DIR"/Sn.psf \
   "$FIXTURE_DIR"/Sn3O4.bands "$FIXTURE_DIR"/Sn3O4.DOS "$FIXTURE_DIR"/Sn3O4.RHO \
   "$FIXTURE_DIR"/Sn3O4.XV "$FIXTURE_DIR"/MESSAGES "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Default archive: includes inputs+essentials, excludes .RHO/MESSAGES ---
echo -e "\n--- Testing default --include set ---"
stb-archive --no-intro --path . -o default.tar.gz > log_default.txt 2>&1
check_success default.tar.gz
check_archive_contains default.tar.gz calc.fdf
check_archive_contains default.tar.gz calc.out
check_archive_contains default.tar.gz O.psf
check_archive_contains default.tar.gz Sn3O4.bands
check_archive_contains default.tar.gz Sn3O4.XV
check_archive_contains default.tar.gz MANIFEST.txt
check_archive_not_contains default.tar.gz Sn3O4.RHO
check_archive_not_contains default.tar.gz MESSAGES
echo "Testing: manifest picked up SystemLabel/run type/energy despite calc.out's generic name"
tar xzOf default.tar.gz MANIFEST.txt > manifest_default.txt
check_contains "SystemLabel       : Sn3O4" manifest_default.txt
check_contains "Run type          : single-point" manifest_default.txt
check_contains "Free energy (eV)  : -4054.427104" manifest_default.txt


# --- 3. --include adds .RHO back ---
echo -e "\n--- Testing --include .RHO ---"
stb-archive --no-intro --path . --include .RHO -o with_rho.tar.gz > log_include.txt 2>&1
check_archive_contains with_rho.tar.gz Sn3O4.RHO


# --- 4. --exclude removes .bands ---
echo -e "\n--- Testing --exclude .bands ---"
stb-archive --no-intro --path . --exclude .bands -o no_bands.tar.gz > log_exclude.txt 2>&1
check_archive_not_contains no_bands.tar.gz Sn3O4.bands
check_archive_contains no_bands.tar.gz Sn3O4.XV


# --- 5. --include/--exclude validation (must start with '.') ---
echo -e "\n--- Testing --include validation (rejects 'RHO' without a leading dot) ---"
stb-archive --no-intro --path . --include RHO > log_badinclude.txt 2>&1
check_exit_code $? 2
check_contains "must start with" log_badinclude.txt


# --- 6. Batch mode ---
echo -e "\n--- Testing batch mode (--path matching multiple directories) ---"
mkdir -p batchA batchB
cp calc.fdf calc.out batchA/
cp calc.fdf calc.out batchB/
stb-archive --no-intro --path "batch*" > log_batch.txt 2>&1
check_contains "2 director" log_batch.txt
check_success batchA_archive.tar.gz
check_success batchB_archive.tar.gz
rm -rf batchA batchB batchA_archive.tar.gz batchB_archive.tar.gz


# --- 7. --save-report ---
echo -e "\n--- Testing --save-report ---"
rm -f stb_archive_report.txt
stb-archive --no-intro --path . -o report_test.tar.gz --save-report > log_savereport.txt 2>&1
check_success stb_archive_report.txt
check_contains "STB-ARCHIVE REPORT" stb_archive_report.txt
rm -f stb_archive_report.txt report_test.tar.gz


# --- 8. Nothing to archive (empty directory) ---
echo -e "\n--- Testing an empty directory (nothing matches --include) ---"
mkdir -p empty_dir
stb-archive --no-intro --path empty_dir > log_empty.txt 2>&1
check_exit_code $? 1
check_contains "nothing to archive" log_empty.txt
rm -rf empty_dir


# --- 9. Interactive path (stb-suite, shortcut 6.7) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 6.7) ---"
rm -f "$TEST_DIR"_archive.tar.gz
printf '6.7\n%s\n\n\n\n0\n' "$TEST_DIR" | stb-suite > log_interactive.txt 2>&1
check_success "$TEST_DIR"_archive.tar.gz
check_contains "Archived" log_interactive.txt
rm -f "$TEST_DIR"_archive.tar.gz


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
