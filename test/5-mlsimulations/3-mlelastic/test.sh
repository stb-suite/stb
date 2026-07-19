#!/bin/bash

# --- Setup ---
# Smoke test for stb-mlelastic (ML Elastic Constants, item 5.3)
#
# Needs the optional 'ml' extra (pip install stb_suite[ml] -- PyTorch +
# mace-torch). The whole file is skipped with a clear message if `mace`
# isn't importable, same gating pattern as test/5-mlsimulations/1-mlmd and
# test/5-mlsimulations/2-mlphonons.
#
# si8.fdf: 8-atom bulk Si (diamond cubic, a=5.43 Ang) -- exercises the 3D
# path (all 6 canonical directions, cubic C11=C22=C33/C12=C13=C23/
# C44=C55=C66 symmetry).
# graphene.fdf: 2-atom primitive graphene monolayer (vacuum along c, same
# fixture as test/2-structures/7-unitcell) -- exercises the 2D path
# (vacuum-blocked zz/yz/xz directions dropped, Lz dilution correction,
# layer modulus/Poisson ratio report). Both verified live to give physically
# sensible values (Si: cubic stiffness matrix, Born-stable; graphene: layer
# modulus ~260-320 N/m, in the right ballpark vs. the ~340 N/m experimental
# value, Born-stable) before this test was written.
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


echo "--- Starting tester for stb-mlelastic (item 5.3) ---"

if ! python3 -c "import mace" 2>/dev/null; then
    echo -e "${YELLOW}Skipped entirely:${NC} the optional 'ml' extra is not installed."
    echo "Install with: pip install stb_suite[ml]  (then re-run this test)"
    exit 0
fi

rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/si8.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/graphene.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 1. 3D bulk Si: full 6-direction cubic stiffness matrix ---
echo -e "\n--- Testing a basic 3D run (bulk Si, all 6 directions) ---"
stb-mlelastic --file si8.fdf --max 1.0 --steps 4 --save-data --save-report \
    --no-intro -o si_elastic > log_si.txt 2>&1
check_exit_code $? 0
check_contains "3D STIFFNESS MATRIX" log_si.txt
check_contains "STABILITY AND PROPERTIES" log_si.txt
check_contains "Young's Modulus" log_si.txt
check_contains "Universal Anisotropy Index" log_si.txt
check_success si_elastic/stress_strain.png
check_success si_elastic/stress_strain.dat
check_success si_elastic/stb_mlelastic_report.txt
check_contains "STB-MLELASTIC REPORT" si_elastic/stb_mlelastic_report.txt
check_contains "Full report" si_elastic/stb_mlelastic_report.txt
check_contains "SOUND VELOCITIES" log_si.txt
check_contains "Debye temperature" log_si.txt
check_contains "DIRECTIONAL ANISOTROPY" log_si.txt
check_contains "Anisotropy ratio" log_si.txt
check_success si_elastic/anisotropy_surface.png

echo "Testing: Debye temperature and sound velocities are physically sensible for Si (100-1200 K, Vs 1000-6000 m/s)"
python3 -c "
import re, sys
# stdout (log_si.txt) still carries ANSI color codes around the numeric
# value -- read the saved report file instead, which print_dual wrote with
# color codes already stripped.
text = open('si_elastic/stb_mlelastic_report.txt').read()
m = re.search(r'Debye temperature\s*:\s*([\d.]+) K', text)
m2 = re.search(r'Shear \(Vs\)\s*:\s*([\d.]+) m/s', text)
if not (m and m2):
    sys.exit(1)
theta_d, vs = float(m.group(1)), float(m2.group(1))
sys.exit(0 if (100.0 < theta_d < 1200.0 and 1000.0 < vs < 6000.0) else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} Debye temperature and shear sound velocity in a physically sensible range"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} Debye temperature/sound velocity outside the expected physical range"
    FAIL=$((FAIL+1))
fi

echo "Testing: --skip-anisotropy-surface omits section [8]"
stb-mlelastic --file si8.fdf --max 1.0 --steps 4 --skip-anisotropy-surface \
    --no-intro -o si_noaniso > log_noaniso.txt 2>&1
check_exit_code $? 0
check_not_contains "DIRECTIONAL ANISOTROPY" log_noaniso.txt
rm -rf si_noaniso

echo "Testing: cubic symmetry recovered (C11~=C22~=C33, C12~=C13~=C23, C44~=C55~=C66)"
python3 -c "
import re, sys
text = open('log_si.txt').read()
m = re.search(r'C11:\s*([\d.]+).*C22:\s*([\d.]+).*C33:\s*([\d.]+)', text, re.S)
m2 = re.search(r'C44:\s*([\d.]+).*C55:\s*([\d.]+).*C66:\s*([\d.]+)', text, re.S)
m3 = re.search(r'C12:\s*([\d.]+).*C13:\s*([\d.]+).*C23:\s*([\d.]+)', text, re.S)
if not (m and m2 and m3):
    sys.exit(1)
c11, c22, c33 = (float(x) for x in m.groups())
c44, c55, c66 = (float(x) for x in m2.groups())
c12, c13, c23 = (float(x) for x in m3.groups())
ok = (abs(c11 - c22) < 2 and abs(c22 - c33) < 2 and
      abs(c44 - c55) < 2 and abs(c55 - c66) < 2 and
      abs(c12 - c13) < 2 and abs(c13 - c23) < 2 and
      c11 > 0 and c44 > 0)
sys.exit(0 if ok else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} cubic symmetry recovered from independent per-direction fits"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} cubic symmetry not recovered (unexpected for bulk Si)"
    FAIL=$((FAIL+1))
fi

echo "Testing: Born stability criteria pass for a real, stable crystal"
check_contains "VERDICT:.*STABLE" si_elastic/stb_mlelastic_report.txt
rm -rf si_elastic


# --- 2. 2D graphene: vacuum-blocked directions dropped, layer modulus ---
echo -e "\n--- Testing a 2D run (graphene monolayer) ---"
stb-mlelastic --file graphene.fdf --max 0.5 --steps 4 --no-intro \
    -o graphene_elastic > log_graphene.txt 2>&1
check_exit_code $? 0
check_contains "dimensionality: 2D" log_graphene.txt
check_contains "2D STIFFNESS MATRIX" log_graphene.txt
check_contains "In-Plane Stiffness" log_graphene.txt
check_contains "2D Stability Criteria" log_graphene.txt
check_contains "No data for direction(s) zz, yz, xz" log_graphene.txt
check_success graphene_elastic/stress_strain.png

echo "Testing: graphene layer modulus is in a physically sensible range (150-450 N/m)"
python3 -c "
import re, sys
text = open('log_graphene.txt').read()
m = re.search(r'Ex=([\d.]+) N/m', text)
if not m:
    sys.exit(1)
ex = float(m.group(1))
sys.exit(0 if 150.0 < ex < 450.0 else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} layer modulus in a physically sensible range"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} layer modulus outside the expected physical range"
    FAIL=$((FAIL+1))
fi
rm -rf graphene_elastic


# --- 3. --dimensionality manual override ---
echo -e "\n--- Testing --dimensionality manual override ---"
stb-mlelastic --file graphene.fdf --dimensionality 3d --max 0.5 --steps 4 \
    --no-intro -o graphene_3d_override > log_override.txt 2>&1
check_exit_code $? 0
check_contains "3D (manual override)" log_override.txt
rm -rf graphene_3d_override


# --- 3b. Fine-tune a quick Si model (reusing test/5-mlsimulations/2-mlphonons's
#     8 real finished SIESTA fixtures) to test --custom-model / foundation
#     comparison with. Not expected to be physically accurate (3 epochs on 8
#     rattled configs, no strain data) -- only exercises the comparison
#     wiring/report sections, same scope as test/5-mlsimulations/2-mlphonons's
#     own --custom-model test. ---
echo -e "\n--- Fine-tuning a quick Si model to test --custom-model with (--epochs 3) ---"
PHONONS_FIXTURES="$FIXTURE_DIR/../2-mlphonons/si_mlff_fixtures"
if [ -d "$PHONONS_FIXTURES" ]; then
    cp -r "$PHONONS_FIXTURES" mlff_config_src
    mv mlff_config_src/mlff_config_* .
    rmdir mlff_config_src
    stb-mlffAnalysis --path "mlff_config_*" --epochs 3 --batch-size 4 --device cpu \
        --name si_quick_model --skip-foundation-comparison --no-intro > log_quicktrain.txt 2>&1
    QUICK_MODEL=$(find . -maxdepth 1 -name "si_quick_model*.model" | sort | tail -1)
    rm -rf mlff_config_001 mlff_config_002 mlff_config_003 mlff_config_004 \
           mlff_config_005 mlff_config_006 mlff_config_007 mlff_config_008
else
    QUICK_MODEL=""
fi
if [ -z "$QUICK_MODEL" ]; then
    echo -e "${YELLOW}SKIPPED${NC}: could not produce a quick fine-tuned Si model (fixtures not found)."
else
    echo "Using model: $QUICK_MODEL"
fi


# --- 3c. --custom-model + foundation-model comparison (default on) ---
if [ -n "$QUICK_MODEL" ]; then
    echo -e "\n--- Testing --custom-model with foundation comparison (default on) ---"
    stb-mlelastic --file si8.fdf --custom-model "$QUICK_MODEL" --max 1.0 --steps 4 \
        --no-intro -o si_compare > log_compare.txt 2>&1
    check_exit_code $? 0
    check_contains "Comparison        : also evaluating with foundation model" log_compare.txt
    check_contains "FOUNDATION MODEL COMPARISON" log_compare.txt
    check_contains "Fine-tuned" log_compare.txt
    check_contains "Foundation (small)" log_compare.txt
    check_success si_compare/stress_strain.png
    rm -rf si_compare

    echo "Testing: --skip-foundation-comparison disables the second run"
    stb-mlelastic --file si8.fdf --custom-model "$QUICK_MODEL" --skip-foundation-comparison \
        --max 1.0 --steps 4 --no-intro -o si_nocompare > log_nocompare.txt 2>&1
    check_exit_code $? 0
    check_not_contains "FOUNDATION MODEL COMPARISON" log_nocompare.txt
    rm -rf si_nocompare
else
    echo -e "${YELLOW}SKIPPED${NC}: comparison tests (no quick model available)."
fi


# --- 3d. --check-convergence ---
echo -e "\n--- Testing --check-convergence (strain-amplitude scan) ---"
stb-mlelastic --file si8.fdf --check-convergence --convergence-strains 0.5 1.0 2.0 \
    --save-data --no-intro -o si_conv > log_conv.txt 2>&1
check_exit_code $? 0
check_contains "STRAIN-AMPLITUDE CONVERGENCE CHECK" log_conv.txt
check_contains "MaxStrain%" log_conv.txt
check_success si_conv/strain_convergence.png
check_success si_conv/strain_convergence.dat
rm -rf si_conv

echo "Testing: --convergence-strains needs at least 2 values (expects failure)"
stb-mlelastic --file si8.fdf --check-convergence --convergence-strains 1.0 \
    --no-intro > log_conv_bad.txt 2>&1
check_exit_code $? 2
check_contains "at least 2 values" log_conv_bad.txt


# --- 4. --no-relax / --no-pre-relax / --custom-model errors ---
echo -e "\n--- Testing robustness / error cases ---"

echo "Testing: --no-pre-relax and --no-relax run without crashing"
stb-mlelastic --file si8.fdf --no-pre-relax --no-relax --max 1.0 --steps 3 \
    --no-intro -o si_norelax > log_norelax.txt 2>&1
check_exit_code $? 0
check_contains "no-pre-relax: straining the structure as given" log_norelax.txt
rm -rf si_norelax

echo "Testing: missing --file"
stb-mlelastic --no-intro > log_missing_file.txt 2>&1
check_exit_code $? 2

echo "Testing: nonexistent input file"
stb-mlelastic --file nope.fdf --no-intro > log_missing_input.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_input.txt

echo "Testing: --custom-model with a nonexistent file"
stb-mlelastic --file si8.fdf --custom-model does_not_exist.model --no-intro > log_custommodel.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_custommodel.txt

echo "Testing: --version"
stb-mlelastic --version > log_version.txt 2>&1
check_contains "stb-mlelastic" log_version.txt

echo "Testing: --help documents --max, --steps, --dimensionality, --custom-model, --no-relax"
stb-mlelastic --help > log_help.txt 2>&1
check_contains "\-\-max" log_help.txt
check_contains "steps" log_help.txt
check_contains "dimensionality" log_help.txt
check_contains "custom-model" log_help.txt
check_contains "no-relax" log_help.txt
check_contains "no-pre-relax" log_help.txt
check_contains "fit-quality-tolerance" log_help.txt
check_contains "skip-foundation-comparison" log_help.txt
check_contains "check-convergence" log_help.txt
check_contains "convergence-strains" log_help.txt
check_contains "skip-anisotropy-surface" log_help.txt


# --- 5. Interactive path (stb-suite, shortcut 5.3) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 5.3) ---"
rm -rf interactive_out
printf '5.3\nsi8.fdf\n\nsmall\n1.0\n4\n\ninteractive_out\nn\nn\nn\nn\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success interactive_out/stress_strain.png


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
