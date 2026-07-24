#!/bin/bash

# --- Setup ---
# Smoke test for stb-defect (Point Defect Generator, item 2.5)
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

# Checks that file $1 exists and is not empty
check_success() {
    if [ -s "$1" ]; then
        echo -e " ... ${GREEN}OK${NC} (file '$1' created)"
        PASS=$((PASS+1))
    else
        echo -e " ... ${RED}FAIL${NC} (file '$1' was not created)"
        FAIL=$((FAIL+1))
    fi
}

# Checks that file $2 contains (grep -q) pattern $1
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


# --- 1. Preparation ---
echo "--- Starting tester for STB-Defect (item 2.5) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/structure.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/graphene.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/magnetite.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Vacancy by index ---
echo -e "\n--- Testing vacancy by --index ---"
rm -f vac.fdf
stb-defect -f structure.fdf --type vacancy --index 1 -o vac.fdf --no-intro > log_vacancy.txt 2>&1
check_contains "Selected site.*#1 (Si)" log_vacancy.txt
check_contains "Output atoms.*1" log_vacancy.txt
check_success vac.fdf
check_contains "NumberofAtoms      1" vac.fdf
check_contains "Vacancy defect introduced by stb-defect" vac.fdf


# --- 3. Substitution by index ---
echo -e "\n--- Testing substitution by --index (Si -> Ge) ---"
rm -f sub.fdf
stb-defect -f structure.fdf --type substitution --index 2 --new-species Ge -o sub.fdf --no-intro > log_substitution.txt 2>&1
check_contains "Substitution.*-> Ge" log_substitution.txt
check_contains "Output formula.*SiGe" log_substitution.txt
check_success sub.fdf
check_contains " 2   32   Ge" sub.fdf


# --- 4. Interstitial ---
echo -e "\n--- Testing interstitial ---"
rm -f inter.fdf
stb-defect -f structure.fdf --type interstitial --position 0.5 0.5 0.5 --species N -o inter.fdf --no-intro > log_interstitial.txt 2>&1
check_contains "Output atoms.*3" log_interstitial.txt
check_success inter.fdf
check_contains " 2   7   N" inter.fdf
check_contains "0.50000000   0.50000000   0.50000000   2" inter.fdf


# --- 5. --nearest with periodic wraparound ---
echo -e "\n--- Testing --nearest (graphene.fdf, target just outside the cell) ---"
rm -f vac_near.fdf
stb-defect -f graphene.fdf --type vacancy --nearest 1.05 0.02 0.5 -o vac_near.fdf --no-intro > log_nearest.txt 2>&1
check_contains "Selected site.*#1 (C)" log_nearest.txt
check_success vac_near.fdf

echo -e "\n--- Testing --nearest with --filter-species (no match found) ---"
stb-defect -f graphene.fdf --type vacancy --nearest 0 0 0.5 --filter-species N --no-intro > log_nearest_nomatch.txt 2>&1
check_exit_code $? 1
check_contains "No atoms of species" log_nearest_nomatch.txt


# --- 5b. --all-inequivalent-sites (magnetite.fdf: 2 distinct Fe sites, 1 O site) ---
echo -e "\n--- Testing --all-inequivalent-sites ---"

echo "Testing: vacancy, --all-inequivalent-sites --filter-species Fe (expect 2 sites, 2 files)"
rm -f magvac_site*.fdf
stb-defect -f magnetite.fdf --type vacancy --all-inequivalent-sites --filter-species Fe \
    -o magvac.fdf --no-intro > log_all_sites_fe.txt 2>&1
check_contains "Space group.*Fd-3m" log_all_sites_fe.txt
check_contains "Symmetrically distinct sites found: 2" log_all_sites_fe.txt
check_contains "2 structure(s) written" log_all_sites_fe.txt
check_success magvac_site1.fdf
check_success magvac_site5.fdf
check_contains "NumberofAtoms      13" magvac_site1.fdf
check_contains "Vacancy defect at site #1 (Fe, Wyckoff" magvac_site1.fdf

echo "Testing: substitution, --all-inequivalent-sites, no filter (expect 3 sites: 2 Fe + 1 O)"
rm -f magsub_site*.fdf
stb-defect -f magnetite.fdf --type substitution --new-species Ni --all-inequivalent-sites \
    -o magsub.fdf --no-intro > log_all_sites_all.txt 2>&1
check_contains "Symmetrically distinct sites found: 3" log_all_sites_all.txt
check_contains "3 structure(s) written" log_all_sites_all.txt
check_success magsub_site1.fdf
check_success magsub_site5.fdf
check_success magsub_site7.fdf
check_contains " 3   28   Ni" magsub_site7.fdf

echo "Testing: --all-inequivalent-sites with --filter-species not present (error)"
stb-defect -f magnetite.fdf --type vacancy --all-inequivalent-sites --filter-species Au --no-intro \
    > log_all_sites_nomatch.txt 2>&1
check_exit_code $? 1
check_contains "no atoms of species" log_all_sites_nomatch.txt

echo "Testing: --all-inequivalent-sites and --index together (mutually exclusive)"
stb-defect -f magnetite.fdf --type vacancy --index 1 --all-inequivalent-sites --no-intro \
    > log_all_sites_mutex.txt 2>&1
check_exit_code $? 2

echo "Testing: --all-inequivalent-sites with --type interstitial (invalid combination)"
stb-defect -f magnetite.fdf --type interstitial --position 0 0 0 --species N \
    --all-inequivalent-sites --no-intro > log_all_sites_interstitial.txt 2>&1
check_exit_code $? 2
check_contains "only valid with --type vacancy/substitution" log_all_sites_interstitial.txt


# --- 5c. --ml-rank (needs the optional 'ml' extra; skipped otherwise) ---
echo -e "\n--- Testing --ml-rank ---"

echo "Testing: --ml-rank without --all-inequivalent-sites (invalid, doesn't need mace)"
stb-defect -f magnetite.fdf --type vacancy --index 1 --ml-rank --no-intro > log_ml_rank_mutex.txt 2>&1
check_exit_code $? 2
check_contains "only valid with --all-inequivalent-sites" log_ml_rank_mutex.txt

echo "Testing: --ml-rank and --ml-relax together (mutually exclusive, doesn't need mace)"
stb-defect -f magnetite.fdf --type vacancy --all-inequivalent-sites --ml-rank --ml-relax --no-intro \
    > log_ml_rank_relax_mutex.txt 2>&1
check_exit_code $? 2
check_contains "mutually exclusive" log_ml_rank_relax_mutex.txt

if python3 -c "import mace" 2>/dev/null; then
    echo "Testing: vacancy, --all-inequivalent-sites --filter-species Fe --ml-rank (site #1 ranked above #5)"
    rm -f magrank_site*.fdf
    stb-defect -f magnetite.fdf --type vacancy --all-inequivalent-sites --filter-species Fe --ml-rank \
        -o magrank.fdf --no-intro > log_ml_rank.txt 2>&1
    check_exit_code $? 0
    check_contains "ML-ranked sites" log_ml_rank.txt
    check_contains "1    #1     Fe       d" log_ml_rank.txt
    check_contains "2    #5     Fe       a" log_ml_rank.txt
    check_success magrank_site1.fdf
    check_success magrank_site5.fdf
else
    echo -e "${YELLOW}Skipped:${NC} the optional 'ml' extra is not installed (pip install stb_suite[ml])."
fi


# --- 5d. Numbered report, before/after symmetry table, --save-report ---
echo -e "\n--- Testing the numbered report, before/after symmetry table, and --save-report ---"
rm -f stb_defect_report.txt vac.fdf
stb-defect -f structure.fdf --type vacancy --index 1 --save-report --no-intro > log_report.txt 2>&1
check_success stb_defect_report.txt
check_contains "===== STB-DEFECT REPORT =====" stb_defect_report.txt
check_contains "\[5\] STRUCTURE VALIDATION, SYMMETRY & WRITING OUTPUT FILE" stb_defect_report.txt
check_contains "Property.*Before.*After" stb_defect_report.txt
check_contains "Report         : stb_defect_report.txt" log_report.txt
rm -f stb_defect_report.txt vac.fdf


# --- 5e. --ml-relax (needs the optional 'ml' extra): single defect and multi-site paths ---
echo -e "\n--- Testing --ml-relax (MACE pre-optimization), if the 'ml' extra is available ---"

if python3 -c "import mace" 2>/dev/null; then
    echo "Testing: single defect, --ml-relax --ml-relax-cell"
    rm -f vacrelax.fdf stb_defect_report.txt
    stb-defect -f structure.fdf --type vacancy --index 1 --ml-relax --ml-relax-cell --save-report \
        -o vacrelax.fdf --no-intro > log_ml_relax.txt 2>&1
    check_exit_code $? 0
    check_contains "\[4\] ML PRE-RELAXATION (MACE)" stb_defect_report.txt
    check_success vacrelax.fdf
    check_contains "ML pre-relaxed with MACE-MP-0" vacrelax.fdf
    rm -f vacrelax.fdf stb_defect_report.txt

    echo "Testing: --all-inequivalent-sites, --ml-relax (no ranking table)"
    rm -f magrelax_site*.fdf
    stb-defect -f magnetite.fdf --type vacancy --all-inequivalent-sites --filter-species Fe --ml-relax \
        -o magrelax.fdf --no-intro > log_ml_relax_multi.txt 2>&1
    check_exit_code $? 0
    check_contains "\[4\] ML PRE-RELAXATION (MACE)" log_ml_relax_multi.txt
    check_success magrelax_site1.fdf
    check_contains "ML pre-relaxed with MACE-MP-0" magrelax_site1.fdf

    echo "Testing: --ml-relax-cell without --ml-relax is rejected"
    stb-defect -f structure.fdf --type vacancy --index 1 --ml-relax-cell --no-intro \
        > log_ml_relax_cell_only.txt 2>&1
    check_exit_code $? 2
else
    echo -e "   -> ${YELLOW}SKIPPED${NC}: the optional 'ml' extra is not installed "
    echo "      (pip install stb_suite[ml] to also exercise --ml-relax)."
fi


# --- 6. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: index out of range"
stb-defect -f structure.fdf --type vacancy --index 5 --no-intro > log_oob.txt 2>&1
check_exit_code $? 1
check_contains "out of range" log_oob.txt

echo "Testing: duplicate index"
stb-defect -f structure.fdf --type vacancy --index 1,1 --no-intro > log_dup.txt 2>&1
check_exit_code $? 1
check_contains "duplicate" log_dup.txt

echo "Testing: --index and --nearest together (mutually exclusive)"
stb-defect -f structure.fdf --type vacancy --index 1 --nearest 0 0 0 --no-intro > log_mutex.txt 2>&1
check_exit_code $? 2

echo "Testing: substitution missing --new-species"
stb-defect -f structure.fdf --type substitution --index 1 --no-intro > log_missing_newsp.txt 2>&1
check_exit_code $? 2

echo "Testing: interstitial missing --position/--species"
stb-defect -f structure.fdf --type interstitial --no-intro > log_missing_inter.txt 2>&1
check_exit_code $? 2

echo "Testing: invalid element symbol"
stb-defect -f structure.fdf --type substitution --index 1 --new-species Xx --no-intro > log_bad_element.txt 2>&1
check_exit_code $? 1
check_contains "not a valid Element" log_bad_element.txt

echo "Testing: missing structure file"
stb-defect -f does_not_exist.fdf --type vacancy --index 1 --no-intro > log_missing_file.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_file.txt

echo "Testing: -f/--file missing (required)"
stb-defect --type vacancy --index 1 --no-intro > log_missing_file_arg.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-defect --version > log_version.txt 2>&1
check_contains "stb-defect" log_version.txt

echo "Testing: --help documents index/nearest/position/all-inequivalent-sites/ml-rank"
stb-defect --help > log_help.txt 2>&1
check_contains "index" log_help.txt
check_contains "nearest" log_help.txt
check_contains "all-inequivalent-sites" log_help.txt
check_contains "ml-rank" log_help.txt


# --- 7. Interactive path (stb-suite, menu 2.5) ---
echo -e "\n--- Testing the interactive path via stb-suite (menu 2.5) ---"

echo "Testing: navigate 2.5 -> invalid file then valid -> vacancy -> by index '1' -> default output ->"
echo "         no ml-relax -> no save-report -> no view -> quit"
rm -f defect.fdf
printf '2.5\ndoes_not_exist.fdf\nstructure.fdf\n1\n1\n1\n\n\n\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_contains "File not found" log_menu.txt
check_contains "Structure written to" log_menu.txt
check_success defect.fdf

echo "Testing: navigate 2.5 -> magnetite.fdf -> vacancy -> every inequivalent site -> filter Fe ->"
echo "         no ml-rank -> no ml-relax -> no save-report -> no view -> default output -> quit"
rm -f defect_site1.fdf defect_site5.fdf
printf '2.5\nmagnetite.fdf\n1\n3\nFe\nn\n\n\n\n\n0\n' | stb-suite > log_menu_all_sites.txt 2>&1
check_contains "Symmetrically distinct sites found: 2" log_menu_all_sites.txt
check_success defect_site1.fdf
check_success defect_site5.fdf

if python3 -c "import mace" 2>/dev/null; then
    echo "Testing: navigate 2.5 -> magnetite.fdf -> vacancy -> every inequivalent site -> filter Fe ->"
    echo "         ml-rank yes -> default output -> no save-report -> no view -> quit"
    rm -f defect_site1.fdf defect_site5.fdf
    printf '2.5\nmagnetite.fdf\n1\n3\nFe\ny\n\n\n\n0\n' | stb-suite > log_menu_ml_rank.txt 2>&1
    check_contains "ML-ranked sites" log_menu_ml_rank.txt
    check_success defect_site1.fdf
    check_success defect_site5.fdf
else
    echo -e "${YELLOW}Skipped:${NC} the optional 'ml' extra is not installed (pip install stb_suite[ml])."
fi


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
