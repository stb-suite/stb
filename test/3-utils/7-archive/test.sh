#!/bin/bash

# --- Setup ---
# Smoke test for stb-archive (Archive Packager, item 3.7).
#
# Fixtures are real files from a completed single-point Sn3O4 calculation
# (shared with the other 3-utils fixtures): calc.fdf/calc.out (SystemLabel
# Sn3O4, even though the log isn't named Sn3O4.out -- exercises
# core.siesta_log.find_out_file's generic-filename fallback), O.psf/Sn.psf,
# Sn3O4.bands/.DOS/.XV and MESSAGES (SIESTA's own log; expected to be archived by
# default), plus Sn3O4.RHO (EXCLUDED by default: a real-space grid). calc.fdf
# %includes a structure.fdf that is not there, so the report must say the archive
# cannot rerun the calculation (true for this fixture). An archive holds ONE top
# folder (the archived folder's name) with MANIFEST.txt and SHA256SUMS.
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

check_equal() {
    if [ "$1" = "$2" ]; then
        echo -e "   -> ${GREEN}Verified:${NC} $3 ($1)"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} $3 (got '$1', expected '$2')"
        FAIL=$((FAIL+1))
    fi
}

check_archive_contains() {
    if tar tzf "$1" 2>/dev/null | grep -qx "[^/]*/$2"; then
        echo -e "   -> ${GREEN}Verified:${NC} '$2' present in '$1'"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} '$2' NOT present in '$1'"
        FAIL=$((FAIL+1))
    fi
}

check_archive_not_contains() {
    if tar tzf "$1" 2>/dev/null | grep -qx "[^/]*/$2"; then
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
echo "--- Starting tester for STB-ARCHIVE (item 3.7) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR"/calc.fdf "$FIXTURE_DIR"/calc.out "$FIXTURE_DIR"/O.psf "$FIXTURE_DIR"/Sn.psf \
   "$FIXTURE_DIR"/Sn3O4.bands "$FIXTURE_DIR"/Sn3O4.DOS "$FIXTURE_DIR"/Sn3O4.RHO \
   "$FIXTURE_DIR"/Sn3O4.XV "$FIXTURE_DIR"/MESSAGES "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Default archive (results mode): inputs + logs + small outputs, no grid ---
echo -e "\n--- Testing the default (results) mode ---"
rm -f stb_archive_report*.txt stb_archive_references*.bib
stb-archive --no-intro --path . -o default.tar.gz > log_default.txt 2>&1
check_exit_code $? 0
check_success default.tar.gz
echo "Testing: one top folder in the archive, named after the archived folder"
top=$(tar tzf default.tar.gz | cut -d/ -f1 | sort -u)
check_equal "$top" "test_files" "top folders in default.tar.gz"
check_archive_contains default.tar.gz calc.fdf
check_archive_contains default.tar.gz calc.out
check_archive_contains default.tar.gz O.psf
check_archive_contains default.tar.gz Sn.psf
check_archive_contains default.tar.gz Sn3O4.bands
check_archive_contains default.tar.gz Sn3O4.XV
check_archive_contains default.tar.gz MESSAGES
check_archive_contains default.tar.gz MANIFEST.txt
check_archive_contains default.tar.gz SHA256SUMS
check_archive_not_contains default.tar.gz Sn3O4.RHO
check_contains "Verification      : read back, every file's size and SHA-256 compared -- ok" log_default.txt
check_contains "Add back with     : --include .RHO" log_default.txt
echo "Testing: the fixture's calc.fdf reads a structure.fdf that is not there -> not rerunnable, warned"
check_contains "Rerunnable        : NO" log_default.txt
check_contains "an .fdf reads structure.fdf, which is not in the folder" log_default.txt
echo "Testing: SHA256SUMS verifies with sha256sum -c after extraction"
mkdir -p extract && tar xzf default.tar.gz -C extract
(cd extract/test_files && sha256sum -c --quiet SHA256SUMS) > log_sha.txt 2>&1
check_exit_code $? 0
echo "Testing: the manifest picked up the run (via stb-status) despite calc.out's generic name"
cp extract/test_files/MANIFEST.txt manifest_default.txt
rm -rf extract
check_contains "SystemLabel       : Sn3O4" manifest_default.txt
check_contains "State             : FINISHED" manifest_default.txt
check_contains "Total energy      : -4054.427104 eV" manifest_default.txt
check_contains "FILES IN THIS ARCHIVE" manifest_default.txt
echo "Testing: the report is always written, in the current folder, with its references file"
check_success stb_archive_report.txt
check_contains "STB-ARCHIVE REPORT" stb_archive_report.txt
check_success stb_archive_references.bib


# --- 3. --include adds .RHO back; the earlier archive (named with -o) stays out ---
echo -e "\n--- Testing --include .RHO ---"
stb-archive --no-intro --path . --include .RHO -o with_rho.tar.gz > log_include.txt 2>&1
check_archive_contains with_rho.tar.gz Sn3O4.RHO
echo "Testing: an earlier archive of any name (default.tar.gz, written with -o) never goes in"
check_archive_not_contains with_rho.tar.gz default.tar.gz


# --- 4. --exclude removes .bands ---
echo -e "\n--- Testing --exclude .bands ---"
stb-archive --no-intro --path . --exclude .bands -o no_bands.tar.gz > log_exclude.txt 2>&1
check_archive_not_contains no_bands.tar.gz Sn3O4.bands
check_archive_contains no_bands.tar.gz Sn3O4.XV
check_archive_not_contains no_bands.tar.gz with_rho.tar.gz


# --- 4b. --mode inputs: only what the run needs; inputs are in every mode ---
echo -e "\n--- Testing --mode inputs ---"
stb-archive --no-intro --path . --mode inputs -o inputs.tar.gz > log_inputs.txt 2>&1
check_archive_contains inputs.tar.gz calc.fdf
check_archive_contains inputs.tar.gz O.psf
check_archive_not_contains inputs.tar.gz calc.out
check_archive_not_contains inputs.tar.gz Sn3O4.bands


# --- 4c. Never overwrites: an existing -o needs --force ---
echo -e "\n--- Testing that an existing archive is not overwritten ---"
stb-archive --no-intro --path . -o default.tar.gz > log_exists.txt 2>&1
check_exit_code $? 1
check_contains "exists (pass --force to overwrite)" log_exists.txt


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
check_contains "batch\* (2 folders, one archive each)" log_batch.txt
check_contains "Archived          : 2 of 2 folder(s)" log_batch.txt
check_success batchA_archive.tar.gz
check_success batchB_archive.tar.gz
echo "Testing: a second run writes _1 next to the first archive (never overwritten)"
stb-archive --no-intro --path batchA > log_batch2.txt 2>&1
check_success batchA_archive_1.tar.gz
rm -rf batchA batchB batchA_archive*.tar.gz batchB_archive*.tar.gz


# --- 7. --save-report is a hidden, ignored option (old scripts keep working) ---
echo -e "\n--- Testing the old --save-report flag ---"
stb-archive --no-intro --path . -o report_test.tar.gz --save-report > log_savereport.txt 2>&1
check_exit_code $? 0
check_success report_test.tar.gz
rm -f report_test.tar.gz


# --- 8. Nothing to archive (empty directory) ---
echo -e "\n--- Testing an empty directory ---"
mkdir -p empty_dir
stb-archive --no-intro --path empty_dir > log_empty.txt 2>&1
check_exit_code $? 1
check_contains "no file to archive" log_empty.txt
check_contains "Archived          : 0 of 1 folder(s)" log_empty.txt
rm -rf empty_dir


# --- 9. Interactive path (stb-suite, shortcut 3.7) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 3.7) ---"
rm -f "$TEST_DIR"_archive*.tar.gz
# Prompts in order: folder, mode, recursive, extra suffixes, suffixes to leave out, dry run
printf '3.7\n%s\n1\nn\n\n\nn\n\n0\n' "$TEST_DIR" | stb-suite > log_interactive.txt 2>&1
check_success "$TEST_DIR"_archive.tar.gz
check_contains "Archived          : 1 of 1 folder(s)" log_interactive.txt
rm -f "$TEST_DIR"_archive*.tar.gz


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
