#!/bin/bash

# --- Setup ---
# Smoke test for stb-opticalAnalysis (Optical Properties Stage 2: Analysis,
# item 4.16.2). Chains a real stb-optical (Stage 1, tested separately) run
# on 3 fixtures (bulk.fdf, slab.fdf, mol.fdf -- all from ../prep/), then
# fabricates SystemLabel.EPSIMG + calc.out files with HAND-CHOSEN, constant
# eps2(E) values in each dir_* folder, so the off-diagonal eps_ij
# reconstruction and the 2D/0D dimensionality corrections can be checked
# against an analytically known answer (same "worked backward from a known
# answer" discipline as the HER/OER/GQCA analysis test.sh files):
#
#   eps2_xx = 2.0, eps2_yy = 4.0, eps2_zz = 6.0 (bulk.fdf, all constant in E)
#   eps2_xy (mixed/bisector direction) = (eps2_xx+eps2_yy)/2 + 0.5 = 3.5
#     -> reconstruction must recover eps2_xy = 0.5000 EXACTLY:
#        eps_xy = eps_(bisector) - (eps_xx+eps_yy)/2 = 3.5 - 3.0 = 0.5
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREP_DIR="$(cd "$FIXTURE_DIR/../prep" && pwd)"
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

# Writes a constant-eps2(E) SystemLabel.EPSIMG (20 points, 0.5-5.0 eV) into
# one direction folder.
write_epsimg() {
    local folder="$1" eps2_const="$2"
    python3 -c "
import numpy as np
omega = np.linspace(0.5, 5.0, 20)
with open('$folder/siesta.EPSIMG', 'w') as f:
    f.write('# synthetic eps2(E), constant = $eps2_const\n')
    for e in omega:
        f.write(f'{e:.6f}  {$eps2_const:.6f}\n')
"
}

# Writes a minimal calc.out with a clean SCF-convergence line into one
# direction folder.
write_calc_out() {
    local folder="$1"
    printf 'siesta: SCF Convergence by DM criterion\nSCF cycle converged after 10 iterations\n' \
        > "$folder/calc.out"
}


# --- 1. Preparation ---
echo "--- Starting tester for STB-OPTICALANALYSIS stage 2: analysis (item 4.16.2) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$PREP_DIR/bulk.fdf" "$TEST_DIR/"
cp "$PREP_DIR/slab.fdf" "$TEST_DIR/"
cp "$PREP_DIR/mol.fdf" "$TEST_DIR/"
cp "$PREP_DIR/calc.fdf" "$TEST_DIR/"
for sym in Na Cl C; do
    echo "# placeholder pseudopotential" > "$TEST_DIR/$sym.psf"
done
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null

echo -e "\n--- Running stb-optical (Stage 1) on bulk.fdf: xx, yy, zz, xy ---"
rm -rf bulk_study
stb-optical -f bulk.fdf -c calc.fdf -O bulk_study --directions xx yy zz xy --no-intro \
    > /dev/null 2>&1
check_success bulk_study/dir_xx/structure.fdf
check_success bulk_study/dir_yy/structure.fdf
check_success bulk_study/dir_zz/structure.fdf
check_success bulk_study/dir_xy/structure.fdf


# --- 2. Missing-results guard (folders exist, but no .EPSIMG yet) ---
echo -e "\n--- Testing the missing-.EPSIMG guard (before fabricating any data) ---"
stb-opticalAnalysis --directory bulk_study --no-intro > log_no_epsimg.txt 2>&1
check_exit_code $? 1
check_contains "Could not read '.EPSIMG'" log_no_epsimg.txt
check_contains "No usable direction results" log_no_epsimg.txt

echo "Testing: a directory with no direction folders at all is rejected"
rm -rf empty_study
mkdir -p empty_study
stb-opticalAnalysis --directory empty_study --no-intro > log_empty.txt 2>&1
check_exit_code $? 1
check_contains "run stb-optical first" log_empty.txt


# --- 3. Fabricate known eps2(E) values, then run the real analysis ---
echo -e "\n--- Fabricating known eps2(E) values (xx=2.0, yy=4.0, zz=6.0, xy_mixed=3.5) ---"
write_epsimg bulk_study/dir_xx 2.0
write_epsimg bulk_study/dir_yy 4.0
write_epsimg bulk_study/dir_zz 6.0
write_epsimg bulk_study/dir_xy 3.5
for d in dir_xx dir_yy dir_zz dir_xy; do
    write_calc_out "bulk_study/$d"
done

echo -e "\n--- Gate 1: full analysis (isotropic average + off-diagonal reconstruction) ---"
stb-opticalAnalysis --directory bulk_study --no-intro > log_gate1.txt 2>&1
check_exit_code $? 0
check_contains "Direction xx (bulk_study/dir_xx):" log_gate1.txt
check_contains "Direction xy (bulk_study/dir_xy):" log_gate1.txt

echo "Testing: isotropic average ('avg') is computed once xx/yy/zz are all present"
check_contains "Direction avg ((isotropic average of xx, yy, zz)):" log_gate1.txt

echo "Testing: off-diagonal eps_xy is reconstructed to EXACTLY the known value (0.5000)"
check_contains "\[2b\] OFF-DIAGONAL DIELECTRIC TENSOR" log_gate1.txt
check_contains "eps_xy (from dir_xy, dir_xx, dir_yy):" log_gate1.txt
check_contains "eps2_xy(E->0)  : 0.5000" log_gate1.txt

echo "Testing: eps2_xy = 0.5 independently, reading the .dat file's own column (not just the log)"
python3 -c "
with open('bulk_study/optical_results_xy_offdiag.dat') as f:
    lines = [l for l in f if not l.startswith('#')]
cols = lines[0].split()
e, eps1_xy, eps2_xy = float(cols[0]), float(cols[1]), float(cols[2])
assert abs(eps2_xy - 0.5) < 1e-6, f'eps2_xy={eps2_xy}, expected 0.5'
print('OK')
" > log_offdiag_check.txt 2>&1
check_contains "OK" log_offdiag_check.txt

echo "Testing: all expected output files were written"
check_success bulk_study/optical_results.csv
check_success bulk_study/optical_results_xx.dat
check_success bulk_study/optical_results_avg.dat
check_success bulk_study/optical_results_offdiagonal.csv
check_success bulk_study/optical_results_xy_offdiag.dat
check_success bulk_study/optical_results_offdiagonal.gplot
check_success bulk_study/optical_results.gplot


# --- 4. Gate 2: a biaxial direction WITHOUT its matching diagonal pair ---
echo -e "\n--- Gate 2: xy present but yy missing -- no reconstruction possible ---"
rm -rf partial_study
stb-optical -f bulk.fdf -c calc.fdf -O partial_study --directions xx xy --no-intro \
    > /dev/null 2>&1
write_epsimg partial_study/dir_xx 2.0
write_epsimg partial_study/dir_xy 3.5
write_calc_out partial_study/dir_xx
write_calc_out partial_study/dir_xy

stb-opticalAnalysis --directory partial_study --no-intro > log_gate2.txt 2>&1
check_exit_code $? 0
check_not_contains "\[2b\] OFF-DIAGONAL DIELECTRIC TENSOR" log_gate2.txt
if [ -f partial_study/optical_results_offdiagonal.csv ]; then
    echo -e "   -> ${RED}Failed:${NC} optical_results_offdiagonal.csv unexpectedly written"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} optical_results_offdiagonal.csv NOT written (as expected)"
    PASS=$((PASS+1))
fi


# --- 5. Dimensionality correction: 3D bulk (no-op) ---
echo -e "\n--- Testing --dimensionality-correction on a 3D (bulk) input (no-op) ---"
stb-opticalAnalysis --directory bulk_study --dimensionality-correction --no-intro \
    > log_3d_correction.txt 2>&1
check_exit_code $? 0
check_contains "3D (bulk) input -- nothing to correct." log_3d_correction.txt


# --- 6. Dimensionality correction: 2D slab (perpendicular vs. parallel) ---
echo -e "\n--- Running stb-optical (Stage 1) on slab.fdf: xx, yy, zz ---"
rm -rf slab_study
stb-optical -f slab.fdf -c calc.fdf -O slab_study --directions xx yy zz --no-intro \
    > /dev/null 2>&1
write_epsimg slab_study/dir_xx 3.0
write_epsimg slab_study/dir_yy 3.0
write_epsimg slab_study/dir_zz 1.0
for d in dir_xx dir_yy dir_zz; do
    write_calc_out "slab_study/$d"
done

echo -e "\n--- Testing --dimensionality-correction on a 2D input (vacuum along z) ---"
stb-opticalAnalysis --directory slab_study --dimensionality-correction --thickness 3.35 \
    --no-intro > log_2d_correction.txt 2>&1
check_exit_code $? 0
check_contains "Direction xx (parallel)" log_2d_correction.txt
check_contains "Direction yy (parallel)" log_2d_correction.txt
check_contains "Direction zz (perpendicular)" log_2d_correction.txt
check_success slab_study/optical_results_xx_2Dcorrected.dat
check_success slab_study/optical_results_2Dcorrected.csv

echo "Testing: --dimensionality-correction on a 2D input without --thickness is rejected"
stb-opticalAnalysis --directory slab_study --dimensionality-correction --no-intro \
    > log_2d_nothickness.txt 2>&1
check_exit_code $? 1
check_contains "requires --thickness" log_2d_nothickness.txt

echo "Testing: --thickness larger than the vacuum axis length is rejected"
stb-opticalAnalysis --directory slab_study --dimensionality-correction --thickness 25 \
    --no-intro > log_2d_badthickness.txt 2>&1
check_exit_code $? 1
check_contains "must be smaller than" log_2d_badthickness.txt


# --- 7. Dimensionality correction: 0D isolated atom (molecular polarizability) ---
echo -e "\n--- Running stb-optical (Stage 1) on mol.fdf: xx, yy, zz ---"
rm -rf mol_study
stb-optical -f mol.fdf -c calc.fdf -O mol_study --directions xx yy zz --no-intro \
    > /dev/null 2>&1
write_epsimg mol_study/dir_xx 1.0
write_epsimg mol_study/dir_yy 1.0
write_epsimg mol_study/dir_zz 1.0
for d in dir_xx dir_yy dir_zz; do
    write_calc_out "mol_study/$d"
done

echo -e "\n--- Testing --dimensionality-correction on a 0D input (molecular polarizability) ---"
stb-opticalAnalysis --directory mol_study --dimensionality-correction --no-intro \
    > log_0d_correction.txt 2>&1
check_exit_code $? 0
check_contains "Extracting molecular polarizability" log_0d_correction.txt
check_contains "Direction xx -- alpha(E->0)" log_0d_correction.txt
check_success mol_study/optical_results_polarizability.csv
check_success mol_study/optical_results_xx_polarizability.dat


# --- 8. Energy-range trimming ---
echo -e "\n--- Testing --energy-min/--energy-max trimming ---"
stb-opticalAnalysis --directory bulk_study --energy-min 2.0 --energy-max 4.0 \
    -o trimmed --no-intro > log_trim.txt 2>&1
check_exit_code $? 0
python3 -c "
with open('bulk_study/trimmed_xx.dat') as f:
    rows = [l for l in f if not l.startswith('#')]
energies = [float(r.split()[0]) for r in rows]
assert all(2.0 - 1e-9 <= e <= 4.0 + 1e-9 for e in energies), energies
assert len(energies) > 0
print('OK')
" > log_trim_check.txt 2>&1
check_contains "OK" log_trim_check.txt


# --- 9. --plot-quantity all ---
echo -e "\n--- Testing --plot-quantity all ---"
stb-opticalAnalysis --directory bulk_study --plot-quantity all -o allq --no-intro \
    > log_allq.txt 2>&1
check_exit_code $? 0
check_success bulk_study/allq_eps.gplot
check_success bulk_study/allq_alpha.gplot
check_success bulk_study/allq_sigma1.gplot


# --- 10. --experimental overlay ---
echo -e "\n--- Testing --experimental overlay ---"
cat > fake_experimental.dat << 'EOF'
1.0  0.10
2.0  0.55
3.0  0.90
4.0  0.40
EOF
stb-opticalAnalysis --directory bulk_study --plot-quantity alpha --experimental fake_experimental.dat \
    -o exp --no-intro > log_experimental.txt 2>&1
check_exit_code $? 0
check_contains "\[3b\] EXPERIMENTAL COMPARISON" log_experimental.txt
check_success bulk_study/exp_experimental.dat

echo "Testing: --experimental with a paired quantity (eps) is rejected"
stb-opticalAnalysis --directory bulk_study --plot-quantity eps --experimental fake_experimental.dat \
    --no-intro > log_experimental_bad.txt 2>&1
check_exit_code $? 1
check_contains "\[ERROR\]" log_experimental_bad.txt


# --- 11. Error cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: --version"
stb-opticalAnalysis --version > log_version.txt 2>&1
check_contains "stb-opticalAnalysis" log_version.txt

echo "Testing: --help documents the diagonal/biaxial naming and off-diagonal reconstruction"
stb-opticalAnalysis --help > log_help.txt 2>&1
check_contains "xx" log_help.txt
check_contains "biaxial" log_help.txt
check_contains "OFF-DIAGONAL RECONSTRUCTION" log_help.txt

echo "Testing: a nonexistent --directory is rejected"
stb-opticalAnalysis --directory does_not_exist --no-intro > log_bad_dir.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_bad_dir.txt


# --- 12. Interactive path (stb-suite, shortcut 4.16.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.16.2) ---"

echo "Testing: navigate 4.16.2 -> bulk_study -> defaults -> quit"
{
  echo "4.16.2"
  echo "bulk_study"    # directory
  echo ""               # output filename (default calc.out)
  echo ""                # plot quantity (default eps)
  echo ""                 # dimensionality correction (default N)
  echo ""                  # advanced settings (default N)
  echo ""                   # press enter to continue
  echo "0"                   # quit stage submenu
} | stb-suite > log_menu.txt 2>&1
check_contains "eps2_xy(E->0)  : 0.5000" log_menu.txt


popd > /dev/null

# --- 13. Summary ---
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
