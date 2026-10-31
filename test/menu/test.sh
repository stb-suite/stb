#!/bin/bash

# --- Setup ---
# Smoke test for the menu's plugin hook (stb_suite._load_plugins): extension
# packages add entries to the interactive menu through the "stb.plugins" entry-point
# group. check_plugins.py exercises it with fake plugins (no package installed):
# merging, numeric order, a taken number, a plugin that fails to load, and the
# stb.core.menu API extensions import. Needs the stb-suite package; skipped if absent.
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

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

echo "--- Starting tester for the menu's plugin hook ---"
if ! python3 -c "import stb.stb_suite" 2>/dev/null; then
    echo -e "${YELLOW}Skipped entirely:${NC} the stb-suite package is not installed."
    exit 0
fi

python3 "$FIXTURE_DIR/check_plugins.py"
check_exit_code $? 0

echo -e "\n--- Tests Complete ---"
echo -e "${GREEN}Passed: $PASS${NC}   ${RED}Failed: $FAIL${NC}"
[ "$FAIL" -eq 0 ]
