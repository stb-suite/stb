#!/bin/bash

# --- Setup ---
# Smoke test for the suite's single version (stb/_version.py): `1.9.<commits>`,
# shared by every tool. Runs check_versions.py, which checks that stb.__version__
# is the number of commits, that every command prints it for --version, that
# no tool defines a version of its own, that the banners carry one year, and
# that building/installing the package carries the number through (on a
# throw-away git repository, see the script's docstring).
#
# Needs the stb-suite package with its 'ml' extra installed (the ML tools refuse
# to start without MACE), plus git and setuptools; skipped with a message if not.
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$FIXTURE_DIR/test_files"

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
NC='\033[0m'

PASS=0
FAIL=0

check_exit_code() {
    if [ "$1" -eq "$2" ]; then
        echo -e "   -> ${GREEN}Verified:${NC} exit code $1 (expected $2)"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} exit code $1 (expected $2)"
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


echo "--- Starting tester for the suite's version (stb/_version.py) ---"

if ! python3 -c "import stb, mace" 2>/dev/null || ! command -v git > /dev/null; then
    echo -e "${YELLOW}Skipped entirely:${NC} needs the stb-suite package with its 'ml' extra, and git."
    echo "Install with: pip install -e \"stb-suite[ml]\"  (then re-run this test)"
    exit 0
fi

rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"

echo -e "\n--- Testing the version, every command, the banners and the build ---"
python3 "$FIXTURE_DIR/check_versions.py" "$TEST_DIR/work" > "$TEST_DIR/log_versions.txt" 2>&1
check_exit_code $? 0
sed 's/^/      /' "$TEST_DIR/log_versions.txt"
check_contains "stb.__version__ = 1.9." "$TEST_DIR/log_versions.txt"
check_contains "commands print '... 1.9." "$TEST_DIR/log_versions.txt"
check_contains "VERSION is 1.9." "$TEST_DIR/log_versions.txt"
check_contains "no tool defines its own VERSION literal" "$TEST_DIR/log_versions.txt"
check_contains "the banners carry one year" "$TEST_DIR/log_versions.txt"
check_contains "a build in a 2-commit checkout is version 1.9.2" "$TEST_DIR/log_versions.txt"
check_contains "an installed copy (no checkout) reports its metadata version" "$TEST_DIR/log_versions.txt"
check_contains "STB_VERSION sets the built version" "$TEST_DIR/log_versions.txt"
check_contains "0 failure(s)" "$TEST_DIR/log_versions.txt"


# --- Summary ---
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
