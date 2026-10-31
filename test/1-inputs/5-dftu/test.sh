#!/bin/bash

# --- Setup ---
# Smoke test for stb-dftu (DFT+U / Hubbard Block Generator, item 1.5)
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
echo "--- Starting tester for STB-DFTU (item 1.5) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Single species (auto-detected shell) ---
echo -e "\n--- Testing a single species with auto-detected shell ---"
stb-dftu --species Mn --u 3.9 --no-intro > log_single.txt 2>&1
check_exit_code $? 0
check_contains "shell=3d (n=3, l=2)" log_single.txt
check_contains "U=3.900 eV, J=0.000 eV" log_single.txt
check_contains "%block LDAU.proj" log_single.txt
check_contains "Mn   1" log_single.txt
check_contains "n=3    2" log_single.txt
check_contains "3.900    0.000" log_single.txt
check_contains "Literature reference for Mn: 3.90 eV" log_single.txt


# --- 3. Multiple species, explicit shell and J, saved to file ---
echo -e "\n--- Testing multiple species with explicit --shell/--j and -o ---"
stb-dftu --species Fe Co --u 5.3 3.32 --j 0.5 0.5 --shell 3d 3d -o multi_LDAU.fdf --no-intro > log_multi.txt 2>&1
check_exit_code $? 0
check_success multi_LDAU.fdf
check_contains "Fe   1" multi_LDAU.fdf
check_contains "Co   1" multi_LDAU.fdf
check_contains "5.300    0.500" multi_LDAU.fdf
check_contains "3.320    0.500" multi_LDAU.fdf

echo "Testing: exactly one %block/%endblock wraps both species stanzas"
BLOCK_COUNT=$(grep -c "%block LDAU.proj" multi_LDAU.fdf)
ENDBLOCK_COUNT=$(grep -c "%endblock LDAU.proj" multi_LDAU.fdf)
if [ "$BLOCK_COUNT" -eq 1 ] && [ "$ENDBLOCK_COUNT" -eq 1 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} single %block/%endblock pair wraps both stanzas"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} expected exactly one %block/%endblock pair"
    FAIL=$((FAIL+1))
fi


# --- 3b. --fdf (verify a structure and use tabulated Materials Project U values) ---
echo -e "\n--- Testing --fdf + --use-reference (auto species detection + tabulated U) ---"
# Sc: metal, no tabulated value (-> skipped with a warning). Fe: metal, tabulated
# (-> included). O: non-metal (-> skipped silently, no warning).
cat > sc_fe_o.fdf <<'EOF'
NumberOfSpecies    3
NumberofAtoms      3

%block ChemicalSpeciesLabel
 1   21   Sc
 2   26   Fe
 3    8   O
%endblock ChemicalSpeciesLabel

LatticeConstant 1.0 Ang

AtomicCoordinatesFormat  Fractional

%block LatticeVectors
 4.0000000000   0.0000000000   0.0000000000
 0.0000000000   4.0000000000   0.0000000000
 0.0000000000   0.0000000000   4.0000000000
%endblock LatticeVectors

%block AtomicCoordinatesAndAtomicSpecies
  0.000000000   0.000000000   0.000000000   1
  0.500000000   0.500000000   0.500000000   2
  0.250000000   0.250000000   0.250000000   3
%endblock AtomicCoordinatesAndAtomicSpecies
EOF

echo "Testing: --fdf detects species; --species Fe restricts to it; --use-reference fills U"
stb-dftu --fdf sc_fe_o.fdf --species Fe --use-reference --no-intro > log_fdf_fe.txt 2>&1
check_exit_code $? 0
check_contains "Species found in 'sc_fe_o.fdf' : Sc, Fe, O" log_fdf_fe.txt
check_contains "Using tabulated Materials Project reference U values" log_fdf_fe.txt
check_contains "Fe   1" log_fdf_fe.txt
check_contains "5.300    0.000" log_fdf_fe.txt

echo "Testing: --fdf + --use-reference over all species -- untabulated metal (Sc) warns and is"
echo "  skipped, non-metal (O) is silently skipped, tabulated metal (Fe) is used"
stb-dftu --fdf sc_fe_o.fdf --use-reference --no-intro > log_fdf_mixed.txt 2>&1
check_exit_code $? 0
check_contains "No tabulated reference U for metal species" log_fdf_mixed.txt
check_contains "Sc" log_fdf_mixed.txt
check_contains "Fe   1" log_fdf_mixed.txt

echo "Testing: O (non-metal) never appears in the generated block"
if grep -q "^O   1" log_fdf_mixed.txt; then
    echo -e "   -> ${RED}Failed:${NC} non-metal species O should never appear in the block"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} non-metal species O is absent from the block"
    PASS=$((PASS+1))
fi

echo "Testing: --fdf + --use-reference where nothing qualifies is a clean error"
stb-dftu --fdf sc_fe_o.fdf --species Sc --use-reference --no-intro > log_fdf_none.txt 2>&1
check_exit_code $? 1
check_contains "nothing to generate" log_fdf_none.txt

echo "Testing: --fdf + --species not present in the file is a clean error"
stb-dftu --fdf sc_fe_o.fdf --species Mn --use-reference --no-intro > log_fdf_unknown.txt 2>&1
check_exit_code $? 1
check_contains "not found in 'sc_fe_o.fdf'" log_fdf_unknown.txt

echo "Testing: --fdf on a nonexistent file is a clean error"
stb-dftu --fdf does_not_exist.fdf --use-reference --no-intro > log_fdf_notfound.txt 2>&1
check_exit_code $? 1
check_contains "\[ERROR\] Reading 'does_not_exist.fdf'" log_fdf_notfound.txt


# --- 4. --list-reference and --suggest ---
echo -e "\n--- Testing --list-reference ---"
stb-dftu --list-reference --no-intro > log_list.txt 2>&1
check_exit_code $? 0
check_contains "Wang/Maxisch/Ceder" log_list.txt
check_contains "Ni  6.20 eV" log_list.txt
check_contains "starting-point sanity check" log_list.txt

echo "Testing: --list-reference does not require --species/--u"
if grep -q "%block" log_list.txt; then
    echo -e "   -> ${RED}Failed:${NC} --list-reference should not print a block"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} --list-reference prints no block"
    PASS=$((PASS+1))
fi

echo -e "\n--- Testing --suggest ---"
stb-dftu --suggest Ni --no-intro > log_suggest.txt 2>&1
check_exit_code $? 0
check_contains "6.20 eV" log_suggest.txt
check_contains "not a validation" log_suggest.txt

echo "Testing: --suggest for an element with no reference value is a clean error"
stb-dftu --suggest Xx --no-intro > log_suggest_bad.txt 2>&1
check_exit_code $? 1
check_contains "No reference value tabulated" log_suggest_bad.txt


# --- 5. Error and robustness cases ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: --u never used silently -- --species alone is rejected"
stb-dftu --species Mn --no-intro > log_no_u.txt 2>&1
check_exit_code $? 1

echo "Testing: mismatched --species/--u length"
stb-dftu --species Mn Fe --u 3.9 --no-intro > log_mismatch_u.txt 2>&1
check_exit_code $? 1
check_contains "must match one-to-one" log_mismatch_u.txt

echo "Testing: mismatched --species/--j length"
stb-dftu --species Mn Fe --u 3.9 5.3 --j 0.5 --no-intro > log_mismatch_j.txt 2>&1
check_exit_code $? 1
check_contains "must match one-to-one" log_mismatch_j.txt

echo "Testing: mismatched --species/--shell length"
stb-dftu --species Mn Fe --u 3.9 5.3 --shell 3d --no-intro > log_mismatch_shell.txt 2>&1
check_exit_code $? 1
check_contains "must match one-to-one" log_mismatch_shell.txt

echo "Testing: species with no default shell and no --shell given"
stb-dftu --species Si --u 2.0 --no-intro > log_no_shell.txt 2>&1
check_exit_code $? 1
check_contains "No default correlated shell known" log_no_shell.txt

echo "Testing: --shell explicitly supplied for a species with no default works"
stb-dftu --species Si --u 2.0 --shell 3d --no-intro > log_shell_override.txt 2>&1
check_exit_code $? 0

echo "Testing: --version"
stb-dftu --version > log_version.txt 2>&1
check_contains "stb-dftu" log_version.txt

echo "Testing: --help documents species/u/j/shell/list-reference/suggest/save-report"
stb-dftu --help > log_help.txt 2>&1
check_contains "species" log_help.txt
check_contains "list-reference" log_help.txt
check_contains "suggest" log_help.txt
check_contains "linear-response" log_help.txt
check_contains "save-report" log_help.txt


# --- 5b. --save-report and references.bib ---
echo -e "\n--- Testing --save-report and references.bib ---"

echo "Testing: --save-report writes stb_dftu_report.txt matching the console output"
rm -f stb_dftu_report.txt
stb-dftu --species Mn --u 3.9 --save-report --no-intro > log_save_report_console.txt 2>&1
check_success stb_dftu_report.txt
check_contains "U=3.900 eV, J=0.000 eV" stb_dftu_report.txt
check_contains "Report" log_save_report_console.txt

echo "Testing: references.bib always written, with SIESTA citations"
rm -f references.bib
stb-dftu --species Si --u 2.0 --shell 3d --no-intro > /dev/null 2>&1
check_success references.bib
check_contains "@article{Soler2002," references.bib
check_contains "@article{Garcia2020," references.bib

echo "Testing: references.bib also gets the Wang/Maxisch/Ceder citation when the reference table is touched (Mn is tabulated)"
rm -f references.bib
stb-dftu --species Mn --u 3.9 --no-intro > /dev/null 2>&1
check_contains "@article{WangMaxischCeder2006," references.bib

echo "Testing: --list-reference and --suggest also write references.bib (Wang/Maxisch/Ceder only)"
rm -f references.bib
stb-dftu --list-reference --no-intro > /dev/null 2>&1
check_success references.bib
check_contains "@article{WangMaxischCeder2006," references.bib
if grep -q "@article{Soler2002," references.bib; then
    echo -e "   -> ${RED}Failed:${NC} SIESTA citation should not appear for a pure reference-table lookup"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} SIESTA citation absent (no SIESTA-related output was generated)"
    PASS=$((PASS+1))
fi


# --- 6. Interactive path (stb-suite, shortcut 1.5) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 1.5) ---"

cat > fe_only.fdf <<'EOF'
NumberOfSpecies    1
NumberofAtoms      1

%block ChemicalSpeciesLabel
 1   26   Fe
%endblock ChemicalSpeciesLabel

LatticeConstant 1.0 Ang

AtomicCoordinatesFormat  Fractional

%block LatticeVectors
 4.0000000000   0.0000000000   0.0000000000
 0.0000000000   4.0000000000   0.0000000000
 0.0000000000   0.0000000000   4.0000000000
%endblock LatticeVectors

%block AtomicCoordinatesAndAtomicSpecies
  0.000000000   0.000000000   0.000000000   1
%endblock AtomicCoordinatesAndAtomicSpecies
EOF

echo "Testing: navigate 1.5 -> verify .fdf + tabulated U path (option 1)"
printf '1.5\n1\nfe_only.fdf\n\n\n\n0\n' | timeout 20 stb-suite > log_menu_fdf.txt 2>&1
check_exit_code $? 0
check_contains "Species found in 'fe_only.fdf' : Fe" log_menu_fdf.txt
check_contains "Using tabulated Materials Project reference U values" log_menu_fdf.txt
check_contains "Fe   1" log_menu_fdf.txt

echo "Testing: navigate 1.5 -> blank mode choice defaults to option 1 (Materials Project table)"
printf '1.5\n\nfe_only.fdf\n\n\n\n0\n' | timeout 20 stb-suite > log_menu_default.txt 2>&1
check_exit_code $? 0
check_contains "Species found in 'fe_only.fdf' : Fe" log_menu_default.txt
check_contains "Using tabulated Materials Project reference U values" log_menu_default.txt

echo "Testing: navigate 1.5 -> generate by hand (option 2) -> Mn -> U=3.9 -> J default -> auto shell -> finish -> no save -> quit"
printf '1.5\n2\nMn\n3.9\n0.0\n\n\n\n\n\n0\n' | timeout 20 stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_contains "shell=3d (n=3, l=2)" log_menu.txt
check_contains "%block LDAU.proj" log_menu.txt


popd > /dev/null

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
