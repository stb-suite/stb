#!/bin/bash

# --- Setup ---
# Smoke test for the documentation site (mkdocs.yml + docs/gen_pages.py).
# Builds the whole site with `mkdocs build --strict` into test_files/site and
# checks what came out: every menu item has a guide page (a real one from
# examples/, or a generated stub), every console command has a reference page,
# the sidebars follow menu order, links between guides and to the reference
# were written, build tooling was not published, each README keeps the
# list/table/code structure GitHub gives it (check_markdown_structure.py), and
# the committed --help snapshot is up to date (only where the stb-suite package
# and its ml extra are installed).
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

check_not_contains() {
    if grep -q -- "$1" "$2" 2>/dev/null; then
        echo -e "   -> ${RED}Failed:${NC} '$1' found in '$2' but should not be"
        FAIL=$((FAIL+1))
    else
        echo -e "   -> ${GREEN}Verified:${NC} '$1' not found in '$2' (as expected)"
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

echo "Testing: the home page carries the large 'still under development' notice (and only the home page)"
check_contains "This version is still under development" "$SITE_DIR/index.html"
check_contains "Verify every result" "$SITE_DIR/index.html"
check_contains "dev-notice" "$SITE_DIR/index.html"
check_success "$SITE_DIR/assets/notice.css"
check_not_contains "still under development" "$SITE_DIR/guides/index.html"


# --- 2. Every menu item has a guide page (real or stub) ---
echo -e "\n--- Testing that every menu item has a guide page ---"
menu_items=$(python3 -c "import sys; sys.path.insert(0, '$REPO_DIR/docs'); from stbdocs import catalog; print(len(catalog.load_menu()))")
built=$(find "$SITE_DIR/guides" -mindepth 3 -name index.html | wc -l)
check_equal "$built" "$menu_items" "guide pages built vs. top-level menu items"
check_success "$SITE_DIR/guides/1-inputs/1.3-stb-kgrid/index.html"
check_success "$SITE_DIR/guides/2-analysis/2.11-stm/index.html"
check_success "$SITE_DIR/guides/5-workflows/5.8-adsorption/index.html"

echo "Testing: an item with no example gets a generated stub that says so and lists its commands"
stub="$SITE_DIR/guides/3-utils/3.3-grid-to-cube-converter/index.html"
check_success "$stub"
check_contains "no hands-on guide for this tool yet" "$stub"
check_contains 'href="../../../reference/stb-cube/"' "$stub"
echo "Testing: 3.2-stb-clean (examples/3.2-stb-clean/) is a real guide, with its figures"
real62="$SITE_DIR/guides/3-utils/3.2-stb-clean/index.html"
check_success "$real62"
check_absent "$SITE_DIR/guides/3-utils/3.2-clean-file-tools"
check_contains "UseSaveDM" "$real62"
check_success "$SITE_DIR/guides/3-utils/3.2-stb-clean/figures/rerun.png"
check_not_contains "no hands-on guide for this tool yet" "$real62"
echo "Testing: 3.1-stb-translate (examples/3.1-stb-translate/) is a real guide, with its figures"
real61="$SITE_DIR/guides/3-utils/3.1-stb-translate/index.html"
check_success "$real61"
check_absent "$SITE_DIR/guides/3-utils/3.1-file-translator"
check_contains "minimum image" "$real61"
check_success "$SITE_DIR/guides/3-utils/3.1-stb-translate/figures/tolerance.png"
check_not_contains "no hands-on guide for this tool yet" "$real61"
echo "Testing: 7.2-mlmd (examples/7.2-mlmd/) is a real guide, with its figures copied next to the page"
real51="$SITE_DIR/guides/7-ml-simulations/7.2-mlmd/index.html"
check_success "$real51"
check_contains "velocity Verlet" "$real51"
check_success "$SITE_DIR/guides/7-ml-simulations/7.2-mlmd/figures/verlet.png"
echo "Testing: 5.5-convergence (examples/5.5-convergence/ added this session) is a real guide now, not the old stub"
real45="$SITE_DIR/guides/5-workflows/5.5-convergence/index.html"
check_success "$real45"
check_absent "$SITE_DIR/guides/5-workflows/5.5-convergence-tests"
check_contains "5.5 Convergence Tests" "$real45"
check_contains 'href="../../../reference/stb-convergence/"' "$real45"
check_contains 'href="../../../reference/stb-convergenceAnalysis/"' "$real45"
check_not_contains "no hands-on guide for this tool yet" "$real45"
check_contains "Tools without a hands-on guide yet" "$SITE_DIR/guides/index.html"

echo "Testing: 5.13-gqca shows a figure -- the image is copied next to the page and the link works"
gqca="$SITE_DIR/guides/5-workflows/5.13-gqca/index.html"
check_success "$gqca"
check_success "$SITE_DIR/guides/5-workflows/5.13-gqca/figures/ingaas_dashboard.png"
check_contains 'src="figures/ingaas_dashboard.png"' "$gqca"
check_contains 'href="../../../reference/stb-gqcaAnalysis/"' "$gqca"

echo "Testing: 5.14-optical shows two figures -- copied next to the page and linked"
opt16="$SITE_DIR/guides/5-workflows/5.14-optical/index.html"
check_success "$opt16"
check_success "$SITE_DIR/guides/5-workflows/5.14-optical/figures/si_dashboard.png"
check_success "$SITE_DIR/guides/5-workflows/5.14-optical/figures/graphene_dashboard.png"
check_contains 'src="figures/si_dashboard.png"' "$opt16"
check_contains 'href="../../../reference/stb-opticalAnalysis/"' "$opt16"

echo "Testing: 5.16-chargediff shows its four figures -- copied next to the page and linked"
cdd19="$SITE_DIR/guides/5-workflows/5.16-chargediff/index.html"
check_success "$cdd19"
for fig in water_dashboard water_studies nh3bh3_dashboard li_graphene_dashboard; do
    check_success "$SITE_DIR/guides/5-workflows/5.16-chargediff/figures/$fig.png"
    check_contains "src=\"figures/$fig.png\"" "$cdd19"
done
check_contains 'href="../../../reference/stb-chargediffAnalysis/"' "$cdd19"

echo "Testing: 5.17-hirshfeld shows its six figures and links to the three stages' commands"
hi20="$SITE_DIR/guides/5-workflows/5.17-hirshfeld/index.html"
check_success "$hi20"
check_contains "5.17 Hirshfeld-I Charges" "$hi20"
for fig in h2o_dashboard stockholder_weights h2o_studies n2o_dashboard nacl_dashboard mgo_dashboard; do
    check_success "$SITE_DIR/guides/5-workflows/5.17-hirshfeld/figures/$fig.png"
    check_contains "src=\"figures/$fig.png\"" "$hi20"
done
for cmd in stb-hirshfeldPrep stb-hirshfeldIons stb-hirshfeldAnalysis; do
    check_contains "href=\"../../../reference/$cmd/\"" "$hi20"
done
check_not_contains "no hands-on guide for this tool yet" "$hi20"


# --- 3. Sidebar: menu categories, numeric (not alphabetical) order, short labels ---
echo -e "\n--- Testing the generated sidebar ---"
page="$SITE_DIR/guides/1-inputs/1.3-stb-kgrid/index.html"
check_contains "1 · Inputs" "$page"
check_contains "5 · Workflows" "$page"
check_contains "6 · Catalysis" "$page"
check_contains "7 · ML Simulations" "$page"
echo "Testing: the Catalysis section is always shown; with no extension installed it explains where its tools come from"
check_success "$SITE_DIR/guides/6-catalysis/index.html"
check_contains "extension" "$SITE_DIR/guides/6-catalysis/index.html"
echo "Testing: the MLFF workflow is the first item of ML Simulations (7.1, four stages)"
check_success "$SITE_DIR/guides/7-ml-simulations/7.1-mlff/index.html"
check_contains 'href="../../../reference/stb-mlffActiveLearning/"' "$SITE_DIR/guides/7-ml-simulations/7.1-mlff/index.html"
check_contains "1.3 K-Grid Generator" "$page"
check_contains "5.7 Hubbard U (Linear Response)" "$SITE_DIR/guides/5-workflows/5.7-hubbardu/index.html"
echo "Testing: 4.10 sorts after 4.9 (numeric), not after 4.1 (alphabetical)"
python3 - "$SITE_DIR/guides/4-structures/4.1-stb-2Dstacking/index.html" << 'PYEOF'
import re, sys
html = open(sys.argv[1], encoding="utf-8").read()
# the open page also shows up in its own table of contents, so de-duplicate
# (keeping first-appearance order) before comparing
labels = list(dict.fromkeys(re.findall(r'<span class="md-ellipsis">\s*(4\.\d+) ', html)))
sys.exit(0 if labels == [f"4.{i}" for i in range(1, 13)] else 1)
PYEOF
check_exit_code $? 0


# --- 4. Links between guides and back to GitHub ---
echo -e "\n--- Testing link rewriting and the GitHub footer ---"
check_contains 'href="../5.3-cohesive/"' "$SITE_DIR/guides/5-workflows/5.1-strain/index.html"
check_contains 'href="5-workflows/5.9-neb/"' "$SITE_DIR/guides/index.html"
check_contains "github.com/stb-suite/stb/tree/main/examples/1.3-stb-kgrid" "$page"
check_contains "example_1.3.sh" "$page"
check_contains 'href="../../../reference/stb-kgrid/"' "$page"
echo "Testing: a workflow guide links to the reference page of every stage's command"
check_contains 'href="../../../reference/stb-strain/"' "$SITE_DIR/guides/5-workflows/5.1-strain/index.html"
check_contains 'href="../../../reference/stb-strainAnalysis/"' "$SITE_DIR/guides/5-workflows/5.1-strain/index.html"


# --- 4b. Reference: one page per console command ---
echo -e "\n--- Testing the command reference ---"
commands=$(python3 -c "import sys; sys.path.insert(0, '$REPO_DIR/docs'); from stbdocs import catalog; print(len(catalog.load_scripts()))")
built=$(find "$SITE_DIR/reference" -mindepth 2 -name index.html | wc -l)
check_equal "$built" "$commands" "reference pages built vs. commands in [project.scripts]"
check_success "$SITE_DIR/reference/index.html"
ref="$SITE_DIR/reference/stb-kgrid/index.html"
check_success "$ref"
check_contains "usage: stb-kgrid" "$ref"
check_contains "menu code" "$ref"
check_contains 'href="../../guides/1-inputs/1.3-stb-kgrid/"' "$ref"
check_contains "1 · Inputs" "$ref"
check_contains "Not in the menu" "$ref"
echo "Testing: command pages wrap long --help lines with their own stylesheet; guides do not load it"
check_success "$SITE_DIR/assets/reference.css"
check_contains "reference.css" "$ref"
check_not_contains "reference.css" "$page"
check_not_contains "reference.css" "$SITE_DIR/reference/index.html"
echo "Testing: a workflow stage's page names its stage and links to the workflow's guide"
check_contains "Stage 1 - Prep" "$SITE_DIR/reference/stb-strain/index.html"
check_contains 'href="../../guides/5-workflows/5.1-strain/"' "$SITE_DIR/reference/stb-strain/index.html"
echo "Testing: a command outside the menu says so; the menu map lists every code"
check_contains "Not in the" "$SITE_DIR/reference/stb-nebCycle/index.html"
check_contains "5.1.2" "$SITE_DIR/reference/stb-suite/index.html"
check_contains "7.1.4" "$SITE_DIR/reference/stb-suite/index.html"
check_contains "<code>5.5.2</code>" "$SITE_DIR/reference/stb-suite/index.html"


# --- 5. Build tooling is not published ---
echo -e "\n--- Testing that build tooling and the nav file are not published ---"
check_absent "$SITE_DIR/gen_pages.py"
check_absent "$SITE_DIR/hooks.py"
check_absent "$SITE_DIR/dump_help.py"
check_absent "$SITE_DIR/reference_help.json"
check_absent "$SITE_DIR/stbdocs"
check_absent "$SITE_DIR/overrides"
check_absent "$SITE_DIR/requirements.txt"
check_absent "$SITE_DIR/__pycache__"
check_absent "$SITE_DIR/guides/SUMMARY"
check_absent "$SITE_DIR/reference/SUMMARY"
if grep -q "SUMMARY" "$SITE_DIR/sitemap.xml"; then
    echo -e "   -> ${RED}Failed:${NC} SUMMARY listed in sitemap.xml"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} SUMMARY not listed in sitemap.xml"
    PASS=$((PASS+1))
fi


# --- 5b. Extension repositories as extra documentation sources (STB_DOCS_EXTRA_ROOTS) ---
# fixture_plugin/ is a fictional extension (workflow 5.99, command stb-demoext): it proves that
# guides, menu entries and command references from another repository join the site, without
# the site tooling knowing anything about what that repository contains.
echo -e "\n--- Testing extra documentation sources (STB_DOCS_EXTRA_ROOTS) ---"
EXT_SITE="$TEST_DIR/site_ext"
(cd "$REPO_DIR" && STB_DOCS_EXTRA_ROOTS="$FIXTURE_DIR/fixture_plugin" python3 -m mkdocs build --strict -f mkdocs.yml -d "$EXT_SITE") > "$TEST_DIR/log_build_ext.txt" 2>&1
check_exit_code $? 0
demo="$EXT_SITE/guides/5-workflows/5.99-demo/index.html"
check_success "$demo"
check_contains "Demo Extension" "$demo"
check_contains "5.99 Demo Extension" "$demo"
check_contains 'href="../../../reference/stb-demoext/"' "$demo"
check_contains "example.invalid/demo-extension/tree/main/examples/5.99-demo" "$demo"
check_contains "example.invalid/demo-extension/blob/main/examples/5.99-demo/example_5.99.sh" "$demo"
echo "Testing: an extension page has no 'edit this page' link (its source is not in this repository)"
check_not_contains "/edit/main/" "$demo"
check_contains "usage: stb-demoext" "$EXT_SITE/reference/stb-demoext/index.html"
check_contains 'href="../../guides/5-workflows/5.99-demo/"' "$EXT_SITE/reference/stb-demoext/index.html"
check_contains "5.99.1" "$EXT_SITE/reference/stb-suite/index.html"
echo "Testing: the regular pages are still there, and the main site (built without the variable) has none of the extension's"
check_success "$EXT_SITE/guides/1-inputs/1.3-stb-kgrid/index.html"
check_absent "$SITE_DIR/guides/5-workflows/5.99-demo"
check_absent "$SITE_DIR/reference/stb-demoext"
echo "Testing: a root without docs_plugin.toml is rejected with a clear message"
(cd "$REPO_DIR" && STB_DOCS_EXTRA_ROOTS="$FIXTURE_DIR" python3 -m mkdocs build --strict -f mkdocs.yml -d "$TEST_DIR/site_bad") > "$TEST_DIR/log_build_bad.txt" 2>&1
check_exit_code $? 1
check_contains "docs_plugin.toml" "$TEST_DIR/log_build_bad.txt"


# --- 6. GitHub-flavoured Markdown survives Python-Markdown ---
echo -e "\n--- Testing list/table/code structure against GitHub's (CommonMark) rendering ---"
python3 "$FIXTURE_DIR/check_markdown_structure.py" > "$TEST_DIR/log_structure.txt" 2>&1
check_exit_code $? 0
check_contains "0 difference(s) from GitHub's rendering" "$TEST_DIR/log_structure.txt"


# --- 7. The committed --help snapshot matches the tools ---
echo -e "\n--- Testing that docs/reference_help.json is up to date ---"
if python3 -c "import stb, mace" 2>/dev/null; then
    python3 "$REPO_DIR/docs/dump_help.py" --check > "$TEST_DIR/log_help_snapshot.txt" 2>&1
    check_exit_code $? 0
    check_contains "is up to date" "$TEST_DIR/log_help_snapshot.txt"
else
    echo -e "   -> ${YELLOW}Skipped:${NC} needs the stb-suite package and its 'ml' extra installed"
    echo "      (pip install -e \"stb-suite[ml]\"); then: python docs/dump_help.py --check"
fi


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
