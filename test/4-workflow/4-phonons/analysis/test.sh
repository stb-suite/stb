#!/bin/bash

# --- Setup ---
# Smoke test for stb-phononsPos (Phonons Analysis, item 4.4.2), covering both
# ways it can be fed: real SIESTA .FA forces (fabricated here -- same 14-atom
# Sn3O4 system and format as test/5-utils/2-clean/Sn3O4.FA, physically
# meaningless but format-correct, enough to exercise FORCE_SETS extraction),
# and an ML-computed phonopy_disp.yaml (stb-phononsML, no SIESTA needed).
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREP_DIR="$(cd "$FIXTURE_DIR/../prep" && pwd)"
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

# Checks that file $2 does NOT contain (grep -q) pattern $1
check_not_contains() {
    if grep -q "$1" "$2" 2>/dev/null; then
        echo -e "   -> ${RED}Failed:${NC} '$1' unexpectedly found in '$2'"
        FAIL=$((FAIL+1))
    else
        echo -e "   -> ${GREEN}Verified:${NC} '$1' NOT found in '$2' (as expected)"
        PASS=$((PASS+1))
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
echo "--- Starting tester for STB-PhononsPos analysis (item 4.4.2) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$PREP_DIR/structure.fdf" "$TEST_DIR/"
cp "$PREP_DIR/calc.fdf" "$TEST_DIR/"
cp "$PREP_DIR/Sn.psf" "$TEST_DIR/"
cp "$PREP_DIR/O.psf" "$TEST_DIR/"
cp "$FIXTURE_DIR/Sn3O4.FA" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Real-SIESTA-forces path: build displacement folders, drop a
#     fabricated Sn3O4.FA into every disp-*/ (format-correct, not physically
#     meaningful -- fine for a wiring smoke test). Each folder gets a
#     slightly different scale factor rather than the exact same forces
#     everywhere -- identical forces on every displacement make the
#     resulting dynamical matrix perfectly degenerate, which crashes
#     phonopy's thermal-displacement-matrix code with an unrelated assertion
#     error (reproduced once while writing this test; not a stb bug, just an
#     unrealistic fixture) ---
echo -e "\n--- Testing the real-SIESTA-forces path (--bands --dos --pdos --thermal-displacements --freeze-unstable-mode) ---"
rm -rf phonon_runs
stb-phononsCreate -s structure.fdf -c calc.fdf -p . -dim 1 1 1 --no-intro > log_create.txt 2>&1
python3 - <<'PYEOF'
import glob

with open("Sn3O4.FA") as f:
    lines = f.readlines()
header, rows = lines[0], lines[1:]

for i, d in enumerate(sorted(glob.glob("phonon_runs/disp-*")), start=1):
    scale = 1.0 + 0.01 * (i % 17)
    with open(f"{d}/Sn3O4.FA", "w") as out:
        out.write(header)
        for row in rows:
            parts = row.split()
            idx = parts[0]
            vals = [float(v) * scale for v in parts[1:]]
            out.write(f"{idx:>6}  {vals[0]: .9E}  {vals[1]: .9E}  {vals[2]: .9E}\n")
PYEOF
stb-phononsPos -dir phonon_runs -m 4 4 4 --bands --dos --pdos --thermal-displacements \
    --freeze-unstable-mode --no-intro > log_pos_siesta.txt 2>&1
check_contains "Force files found : 84" log_pos_siesta.txt
check_contains "FORCE_SETS       : generated successfully" log_pos_siesta.txt
check_contains "THERMAL PROPERTIES SUMMARY" log_pos_siesta.txt
check_contains "BAND STRUCTURE" log_pos_siesta.txt
check_contains "Bravais lattice" log_pos_siesta.txt
check_contains "DENSITY OF STATES" log_pos_siesta.txt
check_contains "MODE FREEZE" log_pos_siesta.txt
check_success phonon_runs/phonon_properties.txt
check_success phonon_runs/phonon_bands.png
check_success phonon_runs/band.yaml
check_success phonon_runs/thermal_properties.dat
check_success phonon_runs/phonon_plots/bands.dat
check_success phonon_runs/phonon_plots/bands.gplot
check_success phonon_runs/phonon_plots/dos.dat


# --- 3. ML-computed force constants path (no SIESTA .FA needed) -- this is
#     the exact pipeline that surfaced the "bands didn't do anything" menu
#     bug (typing "y" at the old free-text prompt silently dropped --bands);
#     covered again in the interactive section below ---
echo -e "\n--- Testing the ML-computed (stb-phononsML) path with --bands ---"
rm -rf phonon_ml_runs
stb-phononsML -s structure.fdf -dim 1 1 1 -o phonon_ml_runs --no-relax --no-intro > log_ml.txt 2>&1
check_success phonon_ml_runs/phonopy_disp.yaml

stb-phononsPos -dir phonon_ml_runs -m 4 4 4 --bands --no-intro > log_pos_ml.txt 2>&1
check_contains "Force constants already embedded" log_pos_ml.txt
check_contains "BAND STRUCTURE" log_pos_ml.txt
check_success phonon_ml_runs/phonon_bands.png
check_success phonon_ml_runs/band.yaml


# --- 4. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: missing directory"
stb-phononsPos -dir does_not_exist --no-intro > log_missing_dir.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_dir.txt

echo "Testing: directory without phonopy_disp.yaml"
mkdir -p empty_phonon_dir
stb-phononsPos -dir empty_phonon_dir --no-intro > log_no_yaml.txt 2>&1
check_exit_code $? 1
check_contains "Did you run the displacement generator" log_no_yaml.txt

echo "Testing: --version"
stb-phononsPos --version > log_version.txt 2>&1
check_contains "stb-phononsPos" log_version.txt

echo "Testing: --help documents --bands/--dos"
stb-phononsPos --help > log_help.txt 2>&1
check_contains "bands" log_help.txt
check_contains "dos" log_help.txt


# --- 5. Interactive path (stb-suite, shortcut 4.4.2) -- exercises the
#     numbered "Additional analyses" list added to fix the silent-no-op bug ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.4.2) ---"

echo "Testing: select '1' from the numbered list -> bands actually runs"
rm -f phonon_ml_runs/phonon_bands.png phonon_ml_runs/band.yaml
printf '4.4.2\nphonon_ml_runs\n\n4 4 4\n\n\n\n1\nn\n\n0\n' | stb-suite > log_menu_bands.txt 2>&1
check_contains "Selected: bands" log_menu_bands.txt
check_contains "BAND STRUCTURE" log_menu_bands.txt
check_success phonon_ml_runs/phonon_bands.png

echo "Testing: an unrecognized token (old y/n habit) is now flagged instead of silently ignored"
rm -f phonon_ml_runs/phonon_bands.png phonon_ml_runs/band.yaml
printf '4.4.2\nphonon_ml_runs\n\n4 4 4\n\n\n\ny\nn\n\n0\n' | stb-suite > log_menu_typo.txt 2>&1
check_contains "Ignored unrecognized option(s): y" log_menu_typo.txt
check_not_contains "BAND STRUCTURE" log_menu_typo.txt


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
