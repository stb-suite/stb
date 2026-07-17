#!/bin/bash

# --- Setup ---
# Smoke test for stb-gqca (GQCA Stage 1: Cluster Structure Generation,
# item 4.15.1). 2 fixtures: cubic.fdf (3D, bipartite/simple-cubic
# sublattice -- exact -1.0 AB ordering reachable) and mos2.fdf (2D,
# triangular/frustrated Mo sublattice -- exercises the 2D cell-treatment
# fallback and the vacuum-aware AB enumeration filter).
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$FIXTURE_DIR/test_files"

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
        echo -e "   -> ${RED}Failed:${NC} '$1' unexpectedly found in '$2'"
        FAIL=$((FAIL+1))
    else
        echo -e "   -> ${GREEN}Verified:${NC} '$1' NOT found in '$2' (as expected)"
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
echo "--- Starting tester for STB-GQCA stage 1: cluster structure generation (item 4.15.1) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/cubic.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/mos2.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/calc.fdf" "$TEST_DIR/"
for sym in Na K Mo S W; do
    echo "# placeholder pseudopotential" > "$TEST_DIR/$sym.psf"
done
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. 3D fixture (simple-cubic, bipartite) ---
echo -e "\n--- Testing 3D fixture (cubic.fdf, Na sublattice -> K) ---"
rm -rf gqca_study_cubic
stb-gqca -f cubic.fdf --sublattice Na --species-b K -c calc.fdf -p . \
    -O gqca_study_cubic --no-intro > log_cubic.txt 2>&1
check_exit_code $? 0
check_success gqca_study_cubic/cluster_n0/structure.fdf
check_success gqca_study_cubic/cluster_n1/structure.fdf
check_success gqca_study_cubic/cluster_n2/structure.fdf
check_success gqca_study_cubic/gqca_stage1.txt

echo "Testing: every folder forces MD.TypeOfRun CG (real relaxation, not single-point)"
check_contains "MD.TypeOfRun          CG" gqca_study_cubic/cluster_n0/calc.fdf
check_contains "MD.TypeOfRun          CG" gqca_study_cubic/cluster_n1/calc.fdf
check_contains "MD.TypeOfRun          CG" gqca_study_cubic/cluster_n2/calc.fdf

echo "Testing: 3D input gets MD.VariableCell T (full cell relaxes)"
check_contains "MD.VariableCell       T" gqca_study_cubic/cluster_n0/calc.fdf
check_contains "MD.VariableCell       T" gqca_study_cubic/cluster_n1/calc.fdf

echo "Testing: bipartite (simple-cubic) sublattice reaches the ideal -1.0 AB ordering exactly"
check_not_contains "deviates by" gqca_study_cubic/gqca_stage1.txt

echo "Testing: report documents the pair orbit and Sublattice/Species B footer"
check_contains "PAIR ORBIT" gqca_study_cubic/gqca_stage1.txt
check_contains "Sublattice     : Na" gqca_study_cubic/gqca_stage1.txt
check_contains "Species B      : K" gqca_study_cubic/gqca_stage1.txt


# --- 3. 2D fixture (MoS2, frustrated/triangular) ---
echo -e "\n--- Testing 2D fixture (mos2.fdf, Mo sublattice -> W) ---"
rm -rf gqca_study_mos2
stb-gqca -f mos2.fdf --sublattice Mo --species-b W -c calc.fdf -p . \
    -O gqca_study_mos2 --no-intro > log_mos2.txt 2>&1
check_exit_code $? 0
check_success gqca_study_mos2/cluster_n0/structure.fdf
check_success gqca_study_mos2/cluster_n1/structure.fdf
check_success gqca_study_mos2/cluster_n2/structure.fdf
check_success gqca_study_mos2/gqca_stage1.txt

echo "Testing: 2D input is detected and cell relaxation left OFF (positions-only, KNOWN LIMITATION)"
check_contains "KNOWN LIMITATION" gqca_study_mos2/gqca_stage1.txt
check_not_contains "MD.VariableCell       T" gqca_study_mos2/cluster_n0/calc.fdf
check_not_contains "MD.VariableCell       T" gqca_study_mos2/cluster_n1/calc.fdf
check_contains "MD.TypeOfRun          CG" gqca_study_mos2/cluster_n1/calc.fdf

echo "Testing: frustrated (triangular) Mo sublattice does NOT reach the ideal -1.0 exactly"
check_contains "deviates by" gqca_study_mos2/gqca_stage1.txt

echo "Testing: frustrated-lattice Delta_1 sensitivity check ran (larger-cell probe)"
check_contains "sensitivity check" gqca_study_mos2/gqca_stage1.txt

echo "Testing: AB structure's vacuum axis was preserved (not expanded into by the enumeration)"
python3 -c "
import re
with open('gqca_study_mos2/cluster_n1/structure.fdf') as f:
    text = f.read()
m = re.search(r'%block LatticeVectors\n(.*?)%endblock', text, re.DOTALL)
rows = [[float(x) for x in r.split()] for r in m.group(1).strip().split(chr(10))]
lens = [sum(v*v for v in row) ** 0.5 for row in rows]
assert any(abs(l - 33.13) < 0.5 for l in lens), f'no lattice vector close to the 33.13 Ang vacuum length: {lens}'
print('OK')
" > log_mos2_vacuum_check.txt 2>&1
check_contains "OK" log_mos2_vacuum_check.txt


# --- 4. Error cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: --version"
stb-gqca --version > log_version.txt 2>&1
check_contains "stb-gqca" log_version.txt

echo "Testing: --help documents --sublattice/--species-b/--kdensity/--vacuum-gap"
stb-gqca --help > log_help.txt 2>&1
check_contains "sublattice" log_help.txt
check_contains "species-b" log_help.txt
check_contains "kdensity" log_help.txt
check_contains "vacuum-gap" log_help.txt

echo "Testing: unknown --sublattice species is rejected"
stb-gqca -f cubic.fdf --sublattice Xx --species-b K -c calc.fdf -p . \
    -O gqca_study_bad --no-intro > log_bad_sublattice.txt 2>&1
check_exit_code $? 1
check_contains "not found in the structure" log_bad_sublattice.txt

echo "Testing: --sublattice == --species-b is rejected"
stb-gqca -f cubic.fdf --sublattice Na --species-b Na -c calc.fdf -p . \
    -O gqca_study_bad2 --no-intro > log_same_species.txt 2>&1
check_exit_code $? 1
check_contains "must be different species" log_same_species.txt


# --- 5. Interactive path (stb-suite, shortcut 4.15.1) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.15.1) ---"

echo "Testing: navigate 4.15.1 -> cubic fixture -> defaults -> quit"
rm -rf gqca_study_menu
{
  echo "4.15.1"
  echo "cubic.fdf"        # structure
  echo "Na"                # sublattice
  echo "K"                  # species-b
  echo "calc.fdf"           # calc template
  echo "3"                   # pseudo_dir -> option 3 = Custom path
  echo "."                    # custom pseudo path
  echo "gqca_study_menu"       # output dir
  echo ""                       # advanced settings (default N)
  echo ""                        # press enter to continue
  echo "0"                        # quit stage submenu
} | stb-suite > log_menu.txt 2>&1
check_success gqca_study_menu/cluster_n0/structure.fdf
check_success gqca_study_menu/cluster_n1/structure.fdf
check_success gqca_study_menu/cluster_n2/structure.fdf


popd > /dev/null

# --- 6. Summary ---
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
