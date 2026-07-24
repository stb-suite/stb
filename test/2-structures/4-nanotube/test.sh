#!/bin/bash

# --- Setup ---
# Smoke test for stb-nanotube (Nanotube/Nanoribbon Builder, item 2.4)
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
echo "--- Starting tester for STB-Nanotube (item 2.4) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/graphene.fdf" "$TEST_DIR/"
cp "$(cd "$FIXTURE_DIR/../../1-inputs/3-k-path" && pwd)/example_3D_silicon.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Tube mode: (6, 0) zigzag-style tube from graphene ---
echo -e "\n--- Testing tube mode (graphene.fdf, chirality 6 0) ---"
rm -f nanotube.fdf
stb-nanotube -f graphene.fdf --chirality 6 0 --no-intro > log_tube.txt 2>&1
check_contains "Translation vector (t1, t2).*(1, 2)" log_tube.txt
check_contains "Cells in periodic unit (N_cells).*12" log_tube.txt
check_contains "Tube diameter (2R): 4.6983 Ang" log_tube.txt
check_contains "CNT type         : zigzag" log_tube.txt
check_contains "Electronic hint  : metallic" log_tube.txt
check_contains "Output atoms.*24" log_tube.txt
check_success nanotube.fdf
check_contains "NumberofAtoms      24" nanotube.fdf
check_contains "Tube built by stb-nanotube" nanotube.fdf

echo "Testing: physical correctness -- every atom is 3-fold coordinated at the real"
echo "         graphene C-C bond length (~1.41-1.42 Ang), matching ase.build.nanotube"
python3 - > log_tube_physics.txt <<'PYEOF'
from stb.core import structure_io
import numpy as np
s = structure_io.to_pymatgen(structure_io.read_fdf("nanotube.fdf"))
dm = s.distance_matrix
np.fill_diagonal(dm, np.inf)
coord = (dm < 1.6).sum(axis=1)
nn = dm.min(axis=1)
ok = bool((coord == 3).all() and (nn > 1.35).all() and (nn < 1.45).all())
print("PHYSICS_CHECK_OK" if ok else f"PHYSICS_CHECK_FAIL coord={coord.tolist()} nn_range={nn.min():.4f}-{nn.max():.4f}")
PYEOF
check_contains "PHYSICS_CHECK_OK" log_tube_physics.txt


# --- 2b. Generic chiral index with gcd(n,m) > 1: still fully, correctly bonded ---
echo -e "\n--- Testing a generic chiral tube with gcd(n,m) > 1 (graphene.fdf, chirality 8 2) ---"
rm -f chiral82.fdf
stb-nanotube -f graphene.fdf --chirality 8 2 -o chiral82.fdf --no-intro > log_chiral82.txt 2>&1
check_contains "CNT type         : chiral" log_chiral82.txt
check_contains "Note: gcd(n, m) > 1" log_chiral82.txt
check_success chiral82.fdf
python3 - > log_chiral82_physics.txt <<'PYEOF'
from stb.core import structure_io
import numpy as np
s = structure_io.to_pymatgen(structure_io.read_fdf("chiral82.fdf"))
dm = s.distance_matrix
np.fill_diagonal(dm, np.inf)
coord = (dm < 1.6).sum(axis=1)
nn = dm.min(axis=1)
ok = bool((coord == 3).all() and (nn > 1.35).all() and (nn < 1.45).all())
print("PHYSICS_CHECK_OK" if ok else f"PHYSICS_CHECK_FAIL coord={coord.tolist()} nn_range={nn.min():.4f}-{nn.max():.4f}")
PYEOF
check_contains "PHYSICS_CHECK_OK" log_chiral82_physics.txt


# --- 3. Round-trip: tube is a valid, 1D-detected input ---
echo -e "\n--- Verifying the generated tube is valid, vacuum-padded (1D) input for stb-kgrid ---"
stb-kgrid --file nanotube.fdf --density 0.2 --no-intro > log_tube_roundtrip.txt 2>&1
check_contains "Dimensionality : 1D" log_tube_roundtrip.txt
check_contains "Suggested Monkhorst-Pack grid.*1 1 8" log_tube_roundtrip.txt


# --- 4. Ribbon mode: finite width scales with --repeats ---
echo -e "\n--- Testing ribbon mode (graphene.fdf, chirality 6 0, repeats 4) ---"
rm -f ribbon.fdf
stb-nanotube -f graphene.fdf --chirality 6 0 --mode ribbon --repeats 4 -o ribbon.fdf --no-intro > log_ribbon.txt 2>&1
check_contains "Ribbon width.*57.8100 Ang" log_ribbon.txt
check_contains "Output atoms.*96" log_ribbon.txt
check_success ribbon.fdf

echo "Testing: ribbon physical correctness -- exact 1.4203 Ang bond length everywhere"
echo "         (flat tiling is exact regardless of chirality), edge atoms show a real"
echo "         (unpassivated) dangling bond (coordination 2), bulk atoms coordination 3"
python3 - > log_ribbon_physics.txt <<'PYEOF'
from stb.core import structure_io
import numpy as np
s = structure_io.to_pymatgen(structure_io.read_fdf("ribbon.fdf"))
dm = s.distance_matrix
np.fill_diagonal(dm, np.inf)
coord = (dm < 1.6).sum(axis=1)
nn = dm.min(axis=1)
ok = bool(abs(nn.max() - 1.4203) < 1e-3 and set(coord.tolist()) <= {2, 3} and (coord == 2).sum() > 0)
print("RIBBON_PHYSICS_CHECK_OK" if ok else f"RIBBON_PHYSICS_CHECK_FAIL coord_set={set(coord.tolist())} nn_max={nn.max():.4f}")
PYEOF
check_contains "RIBBON_PHYSICS_CHECK_OK" log_ribbon_physics.txt

echo -e "\n--- Verifying the generated ribbon is also a valid 1D input for stb-kgrid ---"
stb-kgrid --file ribbon.fdf --density 0.2 --no-intro > log_ribbon_roundtrip.txt 2>&1
check_contains "Dimensionality : 1D" log_ribbon_roundtrip.txt


# --- 4b. --passivate: a tube has no dangling bonds; a ribbon's edges genuinely do ---
echo -e "\n--- Testing --passivate: a closed tube has 0 dangling bonds ---"
rm -f tube_pass.fdf
stb-nanotube -f graphene.fdf --chirality 6 0 --passivate -o tube_pass.fdf --no-intro > log_tube_pass.txt 2>&1
check_contains "Dangling bonds found : 0" log_tube_pass.txt
check_contains "Auto-passivated      : 0 with H" log_tube_pass.txt
check_contains "Note: no dangling bonds found" log_tube_pass.txt
check_success tube_pass.fdf
check_contains "NumberOfSpecies    1" tube_pass.fdf

echo -e "\n--- Testing --passivate: a ribbon's two edges are genuinely passivated ---"
rm -f ribbon_pass.fdf
stb-nanotube -f graphene.fdf --chirality 6 0 --mode ribbon --repeats 4 --passivate -o ribbon_pass.fdf --no-intro \
    > log_ribbon_pass.txt 2>&1
check_contains "Dangling bonds found : 4" log_ribbon_pass.txt
check_contains "Auto-passivated      : 4 with H" log_ribbon_pass.txt
check_success ribbon_pass.fdf
check_contains "Passivated 4 dangling bond" ribbon_pass.fdf
check_contains " 2   1   H" ribbon_pass.fdf
python3 - > log_ribbon_pass_physics.txt <<'PYEOF'
from stb.core import structure_io
import numpy as np
s = structure_io.to_pymatgen(structure_io.read_fdf("ribbon_pass.fdf"))
dm = s.distance_matrix
np.fill_diagonal(dm, np.inf)
coord = (dm < 1.6).sum(axis=1)
symbols = [str(sp) for sp in s.species]
h_coord = [c for c, sym in zip(coord, symbols) if sym == "H"]
c_coord = [c for c, sym in zip(coord, symbols) if sym == "C"]
ok = bool(len(h_coord) == 4 and all(c == 1 for c in h_coord) and all(c == 3 for c in c_coord))
print("PASSIVATION_PHYSICS_OK" if ok else f"PASSIVATION_PHYSICS_FAIL h={h_coord} c_set={set(c_coord)}")
PYEOF
check_contains "PASSIVATION_PHYSICS_OK" log_ribbon_pass_physics.txt

echo "Testing: --passivant rejected without --passivate"
stb-nanotube -f graphene.fdf --chirality 6 0 --passivant F --no-intro > log_passivant_bad.txt 2>&1
check_exit_code $? 2
check_contains "only valid with --passivate" log_passivant_bad.txt


# --- 5. Numbered report, before/after symmetry table, --save-report ---
echo -e "\n--- Testing the numbered report, before/after symmetry table, and --save-report ---"
rm -f stb_nanotube_report.txt nanotube.fdf
stb-nanotube -f graphene.fdf --chirality 6 0 --save-report --no-intro > log_report.txt 2>&1
check_success stb_nanotube_report.txt
check_contains "===== STB-NANOTUBE REPORT =====" stb_nanotube_report.txt
check_contains "\[6\] SYMMETRY ANALYSIS (BEFORE / AFTER)" stb_nanotube_report.txt
check_contains "Property.*Monolayer (before).*Tube/Ribbon (after)" stb_nanotube_report.txt
check_contains "Layer Group.*p6/mmm" stb_nanotube_report.txt
check_contains "Layer Group.*N/A (not 2D-periodic)" stb_nanotube_report.txt
check_contains "Report         : stb_nanotube_report.txt" log_report.txt
rm -f stb_nanotube_report.txt nanotube.fdf


# --- 6. --ml-relax (needs the optional 'ml' extra) ---
echo -e "\n--- Testing --ml-relax (MACE pre-optimization), if the 'ml' extra is available ---"

if python3 -c "import mace" 2>/dev/null; then
    rm -f nanotube.fdf stb_nanotube_report.txt
    stb-nanotube -f graphene.fdf --chirality 6 0 --ml-relax --ml-relax-cell --save-report --no-intro \
        > log_ml_relax.txt 2>&1
    check_exit_code $? 0
    check_contains "\[4\] ML PRE-RELAXATION (MACE)" stb_nanotube_report.txt
    check_success nanotube.fdf
    check_contains "ML pre-relaxed with MACE-MP-0" nanotube.fdf
    rm -f nanotube.fdf stb_nanotube_report.txt

    echo "Testing: --ml-relax-cell without --ml-relax is rejected"
    stb-nanotube -f graphene.fdf --chirality 6 0 --ml-relax-cell --no-intro > log_ml_relax_cell_only.txt 2>&1
    check_exit_code $? 2
else
    echo -e "   -> ${YELLOW}SKIPPED${NC}: the optional 'ml' extra is not installed "
    echo "      (pip install stb_suite[ml] to also exercise --ml-relax)."
fi


# --- 7. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: chirality (0, 0)"
stb-nanotube -f graphene.fdf --chirality 0 0 --no-intro > log_zero_chir.txt 2>&1
check_exit_code $? 1
check_contains "not valid" log_zero_chir.txt

echo "Testing: --repeats 0"
stb-nanotube -f graphene.fdf --chirality 6 0 --repeats 0 --no-intro > log_zero_repeats.txt 2>&1
check_exit_code $? 1
check_contains "repeats must be" log_zero_repeats.txt

echo "Testing: non-2D input (bulk silicon)"
stb-nanotube -f example_3D_silicon.fdf --chirality 6 0 --no-intro > log_non2d.txt 2>&1
check_exit_code $? 1
check_contains "must be a 2D monolayer" log_non2d.txt

echo "Testing: missing structure file"
stb-nanotube -f does_not_exist.fdf --chirality 6 0 --no-intro > log_missing.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing.txt

echo "Testing: -f/--file missing (required)"
stb-nanotube --chirality 6 0 --no-intro > log_missing_file_arg.txt 2>&1
check_exit_code $? 2

echo "Testing: --chirality missing (required)"
stb-nanotube -f graphene.fdf --no-intro > log_missing_chir_arg.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-nanotube --version > log_version.txt 2>&1
check_contains "stb-nanotube" log_version.txt

echo "Testing: --help documents chirality, vacuum, and passivate options"
stb-nanotube --help > log_help.txt 2>&1
check_contains "Chirality" log_help.txt
check_contains "vacuum" log_help.txt
check_contains "passivate" log_help.txt


# --- 8. Interactive path (stb-suite, menu 2.4) ---
echo -e "\n--- Testing the interactive path via stb-suite (menu 2.4) ---"

echo "Testing: navigate 2.4 -> invalid file then valid -> chirality '6 0' -> mode 1 (tube) -> defaults ->"
echo "         no passivate -> no ml-relax -> no save-report -> no view -> quit"
rm -f nanotube.fdf
printf '2.4\ndoes_not_exist.fdf\ngraphene.fdf\n6 0\n1\n\n\n\n\n\n\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_contains "File not found" log_menu.txt
check_contains "Structure written to" log_menu.txt
check_success nanotube.fdf

echo "Testing: menu defaults the output filename to nanoribbon.fdf when mode=2 (ribbon)"
rm -f nanoribbon.fdf
printf '2.4\ngraphene.fdf\n6 0\n2\n\n\n\n\n\n\n\n0\n' | stb-suite > log_menu_ribbon.txt 2>&1
check_contains "Output file name \[default: nanoribbon.fdf\]" log_menu_ribbon.txt
check_success nanoribbon.fdf


popd > /dev/null

# --- 9. Summary ---
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
