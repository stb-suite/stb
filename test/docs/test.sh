#!/bin/bash

# --- Setup ---
# Smoke test for the documentation site (mkdocs.yml + docs/gen_pages.py).
# Builds the whole site with `mkdocs build --strict` into test_files/site and
# checks what came out: every example became a guide page, the sidebar is
# generated in menu order, links between guides were rewritten, build tooling
# was not published, and each README keeps the list/code structure GitHub gives
# it (check_markdown_structure.py).
#
# Needs the documentation dependencies (not part of the stb-suite package):
#   pip install -r docs/requirements.txt markdown-it-py
# The whole file is skipped with a clear message if they are missing.
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$FIXTURE_DIR/../.." && pwd)"
TEST_DIR="$FIXTURE_DIR/test_files"
SITE_DIR="$TEST_DIR/site"

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
NC='\033[0m'

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

check_absent() {
    if [ ! -e "$1" ]; then
        echo -e "   -> ${GREEN}Verified:${NC} '$1' does not exist (as expected)"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} '$1' exists but should not"
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


echo "--- Starting tester for the documentation site ---"

if ! python3 -c "import mkdocs, mkdocs_gen_files, mkdocs_literate_nav, material, markdown_it" 2>/dev/null; then
    echo -e "${YELLOW}Skipped entirely:${NC} the documentation dependencies are not installed."
    echo "Install with: pip install -r docs/requirements.txt markdown-it-py  (then re-run this test)"
    exit 0
fi

rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
echo "Test directory '$TEST_DIR' prepared."


# --- 1. Strict build ---
echo -e "\n--- Testing the strict build (warnings are errors) ---"
(cd "$REPO_DIR" && python3 -m mkdocs build --strict -f mkdocs.yml -d "$SITE_DIR") > "$TEST_DIR/log_build.txt" 2>&1
check_exit_code $? 0
check_success "$SITE_DIR/index.html"
check_success "$SITE_DIR/guides/index.html"
check_success "$SITE_DIR/assets/stb-icon.png"


# --- 2. Every example became a guide page ---
echo -e "\n--- Testing that every example folder became a guide page ---"
expected=$(ls -d "$REPO_DIR"/examples/[0-9]*/ | wc -l)
built=$(find "$SITE_DIR/guides" -mindepth 3 -name index.html | wc -l)
check_equal "$built" "$expected" "guide pages built vs. example folders"
check_success "$SITE_DIR/guides/1-inputs/1.3-stb-kgrid/index.html"
check_success "$SITE_DIR/guides/3-analysis/3.11-stm/index.html"
check_success "$SITE_DIR/guides/4-workflows/4.8-adsorption/index.html"


# --- 3. Sidebar: menu categories, numeric (not alphabetical) order, short labels ---
echo -e "\n--- Testing the generated sidebar ---"
page="$SITE_DIR/guides/1-inputs/1.3-stb-kgrid/index.html"
check_contains "1 · Inputs" "$page"
check_contains "4 · Workflows" "$page"
check_contains "1.3 K-Grid Generator" "$page"
check_contains "4.7 Hubbard U (Linear Response)" "$SITE_DIR/guides/4-workflows/4.7-hubbardu/index.html"
echo "Testing: 2.10 sorts after 2.9 (numeric), not after 2.1 (alphabetical)"
python3 - "$SITE_DIR/guides/2-structures/2.1-stb-2Dstacking/index.html" << 'PYEOF'
import re, sys
html = open(sys.argv[1], encoding="utf-8").read()
# the open page also shows up in its own table of contents, so de-duplicate
# (keeping first-appearance order) before comparing
labels = list(dict.fromkeys(re.findall(r'<span class="md-ellipsis">\s*(2\.\d+) ', html)))
sys.exit(0 if labels == [f"2.{i}" for i in range(1, 13)] else 1)
PYEOF
check_exit_code $? 0


# --- 4. Links between guides and back to GitHub ---
echo -e "\n--- Testing link rewriting and the GitHub footer ---"
check_contains 'href="../4.3-cohesive/"' "$SITE_DIR/guides/4-workflows/4.1-strain/index.html"
check_contains 'href="4-workflows/4.9-neb/"' "$SITE_DIR/guides/index.html"
check_contains "github.com/stb-suite/stb/tree/main/examples/1.3-stb-kgrid" "$page"
check_contains "example_1.3.sh" "$page"


# --- 5. Build tooling is not published ---
echo -e "\n--- Testing that build tooling and the nav file are not published ---"
check_absent "$SITE_DIR/gen_pages.py"
check_absent "$SITE_DIR/hooks.py"
check_absent "$SITE_DIR/requirements.txt"
check_absent "$SITE_DIR/__pycache__"
check_absent "$SITE_DIR/guides/SUMMARY"
if grep -q "SUMMARY" "$SITE_DIR/sitemap.xml"; then
    echo -e "   -> ${RED}Failed:${NC} SUMMARY listed in sitemap.xml"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} SUMMARY not listed in sitemap.xml"
    PASS=$((PASS+1))
fi


# --- 6. GitHub-flavoured lists survive Python-Markdown ---
echo -e "\n--- Testing list/code structure against GitHub's (CommonMark) rendering ---"
python3 "$FIXTURE_DIR/check_markdown_structure.py" > "$TEST_DIR/log_structure.txt" 2>&1
check_exit_code $? 0
check_contains "documents keep GitHub's list/code structure" "$TEST_DIR/log_structure.txt"


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
