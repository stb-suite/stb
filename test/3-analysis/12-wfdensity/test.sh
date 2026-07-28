#!/bin/bash

# --- Setup ---
# Smoke test for stb-wfdensity (Wavefunction Density, item 3.12)
#
# Reuses test/3-analysis/10-fatbands/'s Sn3O4 fixture (calc.fdf + .ion files
# + Sn3O4.bands.WFSX) -- no new SIESTA run needed, since a band-path WFSX is
# perfectly valid input for picking one (k-index, band) state (unlike
# stb-sts/stb-coop, which need a real full-BZ mesh).
export MPLBACKEND=Agg
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="$FIXTURE_DIR/../10-fatbands"
TEST_DIR="$FIXTURE_DIR/test_files"

# Output colors
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

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


# --- 1. Preparation ---
echo "--- Starting tester for stb-wfdensity (item 3.12) ---"
if [ ! -f "$SOURCE_DIR/Sn3O4.bands.WFSX" ]; then
    echo -e "${RED}FATAL: $SOURCE_DIR/Sn3O4.bands.WFSX not found -- run test/3-analysis/10-fatbands/test.sh's fixture setup first.${NC}"
    exit 1
fi
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$SOURCE_DIR"/Sn3O4.bands.WFSX "$SOURCE_DIR"/calc.fdf "$SOURCE_DIR"/structure.fdf \
   "$SOURCE_DIR"/Sn.ion "$SOURCE_DIR"/Sn.ion.xml "$SOURCE_DIR"/O.ion "$SOURCE_DIR"/O.ion.xml "$TEST_DIR/"
# The interactive stb-suite menu now asks for the .fdf path separately from
# the label (a real fix -- it used to only auto-detect <label>.fdf by exact
# name, with no way to point it at a differently-named real input file like
# this fixture's own calc.fdf). Give the interactive-path tests a same-named
# copy too, so their own default-suggestion prompt has something to accept.
cp "$SOURCE_DIR"/calc.fdf "$TEST_DIR/Sn3O4.fdf"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Basic run: --wfsx/--geometry-file explicit, --k-index 0 --band 1 ---
echo -e "\n--- Testing a basic run (--k-index 0 --band 1) ---"

timeout 60 stb-wfdensity --wfsx Sn3O4.bands.WFSX --geometry-file calc.fdf \
    --k-index 0 --band 1 --spacing 0.4 --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "\[0\] RUN METADATA" log_basic.txt
check_contains "\[1\] INPUT DATA" log_basic.txt
check_contains "\[2\] STATE SELECTION" log_basic.txt
check_contains "\[3\] WAVEFUNCTION DENSITY" log_basic.txt
check_contains "\[4\] OUTPUT DATA & PLOTS" log_basic.txt
check_contains "\[5\] REFERENCES" log_basic.txt
check_contains "\[6\] SUMMARY & FILES" log_basic.txt
check_contains "Normalization   : integral |psi|^2 dV" log_basic.txt
check_success wfdensity_k0_b1.cube

echo "Testing: without --save-gnuplot, no slice/profile .dat/.gplot is written (used to be unconditional, opted out via the now-removed --cube-only)"
if [ -e wfdensity_k0_b1_slice.dat ]; then
    echo -e " ... ${RED}FAIL${NC} (wfdensity_k0_b1_slice.dat should not have been written)"
    FAIL=$((FAIL+1))
else
    echo -e " ... ${GREEN}OK${NC} (wfdensity_k0_b1_slice.dat correctly absent)"
    PASS=$((PASS+1))
fi
check_contains "Slice/profile data not written (off by default" log_basic.txt

echo "Testing: no text report without --save-report"
if [ -e stb_wfdensity_report.txt ]; then
    echo -e " ... ${RED}FAIL${NC} (stb_wfdensity_report.txt should not have been written)"
    FAIL=$((FAIL+1))
else
    echo -e " ... ${GREEN}OK${NC} (stb_wfdensity_report.txt correctly absent)"
    PASS=$((PASS+1))
fi

echo "Testing: references.bib is always written (SIESTA)"
check_success references.bib


# --- 2b. --save-report / --save-gnuplot, and the gnuplot-path fix ---
echo -e "\n--- Testing --save-report / --save-gnuplot (a real, verified path bug fixed here too) ---"

timeout 60 stb-wfdensity --wfsx Sn3O4.bands.WFSX --geometry-file calc.fdf \
    --k-index 0 --band 1 --spacing 0.4 --save-report --save-gnuplot --no-intro > log_saved.txt 2>&1
check_exit_code $? 0
check_success stb_wfdensity_report.txt
check_contains "\[0\] RUN METADATA" stb_wfdensity_report.txt
check_success wfdensity_k0_b1_slice.dat
check_success wfdensity_k0_b1_slice.gplot
check_contains "\[OK\] Data written to" log_saved.txt
check_contains "\[OK\] Gnuplot script written to" log_saved.txt

echo "Testing: the .gplot script's own filenames are bare basenames, no --output-dir prefix"
echo "(this exact bug -- shared core/grid_export.py::write_gnuplot_script -- was also found and"
echo " fixed in stb-density, which uses the same function)"
if grep -q '"wfdensity_k0_b1_slice\.pdf"' wfdensity_k0_b1_slice.gplot && \
   grep -q '"wfdensity_k0_b1_slice\.dat"' wfdensity_k0_b1_slice.gplot; then
    echo -e "   -> ${GREEN}Verified:${NC} wfdensity_k0_b1_slice.gplot references bare basenames"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} wfdensity_k0_b1_slice.gplot has an unexpected path"
    FAIL=$((FAIL+1))
fi
if command -v gnuplot > /dev/null; then
    gnuplot wfdensity_k0_b1_slice.gplot && test -s wfdensity_k0_b1_slice.pdf
    if [ $? -eq 0 ]; then
        echo -e "   -> ${GREEN}Verified:${NC} a real gnuplot run produced wfdensity_k0_b1_slice.pdf"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} gnuplot did not produce wfdensity_k0_b1_slice.pdf"
        FAIL=$((FAIL+1))
    fi
fi

echo "Testing: same check from a genuine --output-dir subfolder"
rm -rf path_check
timeout 60 stb-wfdensity --wfsx Sn3O4.bands.WFSX --geometry-file calc.fdf \
    --k-index 0 --band 1 --spacing 0.5 --save-gnuplot -o path_check --no-intro > log_path_check.txt 2>&1
check_exit_code $? 0
if grep -q '"path_check' path_check/wfdensity_k0_b1_slice.gplot; then
    echo -e "   -> ${RED}Failed:${NC} path_check/wfdensity_k0_b1_slice.gplot still embeds the --output-dir prefix"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no --output-dir prefix embedded"
    PASS=$((PASS+1))
fi


# --- 2c. --view (headless via MPLBACKEND=Agg) ---
echo -e "\n--- Testing --view (headless via MPLBACKEND=Agg, only checking exit code) ---"

timeout 60 stb-wfdensity --wfsx Sn3O4.bands.WFSX --geometry-file calc.fdf \
    --k-index 0 --band 1 --spacing 0.5 --view --no-intro > log_view.txt 2>&1
check_exit_code $? 0

echo "Testing: --view combined with --profile (plot_matplotlib_profile path)"
timeout 60 stb-wfdensity --wfsx Sn3O4.bands.WFSX --geometry-file calc.fdf \
    --k-index 0 --band 1 --spacing 0.5 --profile --view --no-intro > log_view_profile.txt 2>&1
check_exit_code $? 0

echo "Testing: --view alone (no --save-gnuplot) does not falsely report a Data/Gnuplot file in SUMMARY"
if grep -q "^Data " log_view.txt || grep -q "^Gnuplot script " log_view.txt; then
    echo -e "   -> ${RED}Failed:${NC} SUMMARY reported a Data/Gnuplot file that was never written"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} SUMMARY correctly omits Data/Gnuplot lines"
    PASS=$((PASS+1))
fi


# --- 2d. --label + --geometry-file together (a real bug fix): the .WFSX is
# auto-detected from --label, but the geometry comes from an explicitly
# -named fdf that does NOT match the label -- found live from a real user
# report (SystemLabel "siesta", real input file called calc.fdf, no
# siesta.fdf anywhere) that this combination used to be rejected outright.
echo -e "\n--- Testing --label + --geometry-file together (real bug: used to be rejected outright) ---"

timeout 60 stb-wfdensity --label Sn3O4 --geometry-file calc.fdf --k-index 0 --band 1 \
    --spacing 0.5 --no-intro > log_label_plus_geom.txt 2>&1
check_exit_code $? 0
check_contains "WFSX file      : Sn3O4.bands.WFSX" log_label_plus_geom.txt
check_contains "Geometry source: calc.fdf" log_label_plus_geom.txt

echo "Testing: --label + --wfsx together is still rejected (--wfsx IS what --label auto-detects)"
stb-wfdensity --label Sn3O4 --wfsx Sn3O4.bands.WFSX --band 1 --no-intro > log_label_plus_wfsx.txt 2>&1
check_exit_code $? 2
check_contains "cannot be combined with --wfsx" log_label_plus_wfsx.txt


# --- 2e. --mode slice position: auto-detected |psi|^2 peak (default) vs --pos ---
# Real bug fixed here: the slice used to always cut at the geometric center
# of the cell (dim_size // 2), landing in empty vacuum for a structure whose
# atoms aren't centered along --axis. Sn3O4 here is a bulk structure (not the
# vacuum-padded CrS monolayer that originally exposed the bug), so this only
# checks the new mechanism (auto-detect + --pos + bounds) runs and reports
# correctly, not a specific "slice must land on an atom" assertion.
echo -e "\n--- Testing --mode slice position (auto-detected |psi|^2 peak, --pos, bounds) ---"

timeout 60 stb-wfdensity --wfsx Sn3O4.bands.WFSX --geometry-file calc.fdf \
    --k-index 0 --band 1 --spacing 0.5 --save-gnuplot -o slicepos_default_out --no-intro > log_slicepos_default.txt 2>&1
check_exit_code $? 0
check_contains "auto-detected |psi|^2 peak" log_slicepos_default.txt

echo "Testing: --pos overrides the auto-detected default"
timeout 60 stb-wfdensity --wfsx Sn3O4.bands.WFSX --geometry-file calc.fdf \
    --k-index 0 --band 1 --spacing 0.5 --pos 2.0 --save-gnuplot -o slicepos_manual_out --no-intro > log_slicepos_manual.txt 2>&1
check_exit_code $? 0
check_contains "position 2.000 Ang (--pos)" log_slicepos_manual.txt

echo "Testing: --pos out of bounds is a clean error, not a crash"
timeout 60 stb-wfdensity --wfsx Sn3O4.bands.WFSX --geometry-file calc.fdf \
    --k-index 0 --band 1 --spacing 0.5 --pos 999 --save-gnuplot -o slicepos_bad_out --no-intro > log_slicepos_bad.txt 2>&1
check_exit_code $? 1
check_contains "out of bounds" log_slicepos_bad.txt

echo "Testing: --pos is ignored (no effect / no error) in --profile mode"
timeout 60 stb-wfdensity --wfsx Sn3O4.bands.WFSX --geometry-file calc.fdf \
    --k-index 0 --band 1 --spacing 0.5 --pos 2.0 --profile --save-gnuplot -o slicepos_profile_out --no-intro > log_slicepos_profile.txt 2>&1
check_exit_code $? 0
check_success slicepos_profile_out/wfdensity_k0_b1_profile.dat


# --- 3. --profile ---
echo -e "\n--- Testing --profile ---"

timeout 60 stb-wfdensity --wfsx Sn3O4.bands.WFSX --geometry-file calc.fdf \
    --k-index 0 --band 1 --spacing 0.5 --profile --save-gnuplot -o profile_out --no-intro > log_profile.txt 2>&1
check_exit_code $? 0
check_success profile_out/wfdensity_k0_b1_profile.dat
check_success profile_out/wfdensity_k0_b1_profile.gplot


# --- 4. --band vbm/cbm: --fermi, --fermi-file, and .out auto-detection ---
echo -e "\n--- Testing --band vbm (--fermi) ---"

timeout 60 stb-wfdensity --wfsx Sn3O4.bands.WFSX --geometry-file calc.fdf \
    --band vbm --fermi -3.2038563052897451 --spacing 0.6 -o vbm_out --no-intro > log_vbm.txt 2>&1
check_exit_code $? 0
check_contains "VBM found at k-index" log_vbm.txt
check_contains "Fermi source: --fermi (explicit value)" log_vbm.txt

echo -e "\n--- Testing --band vbm with --fermi-file (an arbitrarily-named .out log, decoupled from --label) ---"
cat > my_weird_name.out << 'EOF'
Some SIESTA log noise before the summary
siesta:         Fermi = -3.203872
More noise after
EOF
timeout 60 stb-wfdensity --wfsx Sn3O4.bands.WFSX --geometry-file calc.fdf \
    --band vbm --fermi-file my_weird_name.out --spacing 0.6 -o vbm_fermifile_out --no-intro > log_vbm_fermifile.txt 2>&1
check_exit_code $? 0
check_contains "Fermi source: 'my_weird_name.out' (--fermi-file)" log_vbm_fermifile.txt

echo -e "\n--- Testing --band vbm auto-detecting a generic .out (no --label, decoupled from any SystemLabel) ---"
timeout 60 stb-wfdensity --wfsx Sn3O4.bands.WFSX --geometry-file calc.fdf \
    --band vbm --spacing 0.6 -o vbm_autoout_out --no-intro > log_vbm_autoout.txt 2>&1
check_exit_code $? 0
check_contains "auto-detected .out" log_vbm_autoout.txt
rm -f my_weird_name.out

echo "Testing: --band vbm/cbm with no Fermi source at all is a clean error"
stb-wfdensity --wfsx Sn3O4.bands.WFSX --geometry-file calc.fdf --band vbm --no-intro > log_vbm_missing.txt 2>&1
check_exit_code $? 2
check_contains "needs a Fermi energy" log_vbm_missing.txt


# --- 5. --k-point matches by vector instead of index ---
echo -e "\n--- Testing --k-point ---"

timeout 60 stb-wfdensity --wfsx Sn3O4.bands.WFSX --geometry-file calc.fdf \
    --k-point 0 0 0 --band 2 --spacing 0.6 -o kpoint_out --no-intro > log_kpoint.txt 2>&1
check_exit_code $? 0
check_success kpoint_out/wfdensity_k0_b2.cube


# --- 6. Errors ---
echo -e "\n--- Testing error cases ---"

echo "Testing: --wfsx without --geometry-file"
stb-wfdensity --wfsx Sn3O4.bands.WFSX --band 1 --no-intro > log_missing_geo.txt 2>&1
check_exit_code $? 2

echo "Testing: invalid --band"
stb-wfdensity --wfsx Sn3O4.bands.WFSX --geometry-file calc.fdf --band notanumber --no-intro > log_bad_band.txt 2>&1
check_exit_code $? 2

echo "Testing: --k-index out of range"
stb-wfdensity --wfsx Sn3O4.bands.WFSX --geometry-file calc.fdf --k-index 999 --band 1 --no-intro > log_bad_k.txt 2>&1
check_exit_code $? 2
check_contains "out of range" log_bad_k.txt

echo "Testing: missing --label/--wfsx"
stb-wfdensity --band 1 --no-intro > log_missing_label.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-wfdensity --version > log_version.txt 2>&1
check_contains "stb-wfdensity" log_version.txt

echo "Testing: --help documents --k-index, --k-point, --band, --spacing, --pos, --fermi-file, --geometry-file, --save-report, --save-gnuplot, --view"
stb-wfdensity --help > log_help.txt 2>&1
check_contains "k-index" log_help.txt
check_contains "k-point" log_help.txt
check_contains "band" log_help.txt
check_contains "spacing" log_help.txt
check_contains "\-\-pos" log_help.txt
check_contains "fermi-file" log_help.txt
check_contains "geometry-file" log_help.txt
check_contains "save-report" log_help.txt
check_contains "save-gnuplot" log_help.txt
check_contains "\-\-view" log_help.txt


# --- 7. Interactive path (stb-suite, shortcut 3.12) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 3.12) ---"

echo "Testing: navigate 3.12 -> label Sn3O4 -> band selection 1 -> k-index 0 -> band 1 -> spacing 0.5"
echo "-> axis/mode/pos defaults -> default output -> no save-report -> save-gnuplot=y -> no view"
rm -f wfdensity_k0_b1.cube wfdensity_k0_b1_slice.dat wfdensity_k0_b1_slice.gplot stb_wfdensity_report.txt
printf '3.12\nSn3O4\n\n1\n0\n1\n0.5\n\n\n\n\nn\ny\nn\n\n0\n' | timeout 60 stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success wfdensity_k0_b1.cube
check_success wfdensity_k0_b1_slice.dat
check_success wfdensity_k0_b1_slice.gplot
if [ -e stb_wfdensity_report.txt ]; then
    echo -e " ... ${RED}FAIL${NC} (stb_wfdensity_report.txt should not have been written -- save-report was 'n')"
    FAIL=$((FAIL+1))
else
    echo -e " ... ${GREEN}OK${NC} (stb_wfdensity_report.txt correctly absent)"
    PASS=$((PASS+1))
fi

echo "Testing: navigate 3.12 -> band selection 3 (CBM) -> Fermi source submenu"
echo "(previously: no way in the interactive menu to point at a .out log at all)"
cat > Sn3O4.out << 'EOF'
noise
siesta:         Fermi = -3.203872
noise
EOF
echo "  -> option 3 (explicit .out path)"
printf '3.12\nSn3O4\n\n3\n3\nSn3O4.out\n\n\n\n\n\nn\nn\nn\n\n0\n' | timeout 60 stb-suite > log_menu_fermi3.txt 2>&1
check_exit_code $? 0
check_contains "Fermi source: 'Sn3O4.out' (--fermi-file)" log_menu_fermi3.txt
echo "  -> option 4 (auto-detect, default)"
printf '3.12\nSn3O4\n\n3\n\n\n\n\n\n\nn\nn\nn\n\n0\n' | timeout 60 stb-suite > log_menu_fermi4.txt 2>&1
check_exit_code $? 0
check_contains "auto-detected .out" log_menu_fermi4.txt
rm -f Sn3O4.out


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
