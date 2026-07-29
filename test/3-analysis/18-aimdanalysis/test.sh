#!/bin/bash

# --- Setup ---
# Smoke test for stb-aimdAnalysis (AIMD Trajectory Analysis, item 3.18).
#
# Fixtures are a real SIESTA AIMD run (5-step Verlet, O2 dimer in vacuum,
# same aimd.ANI/.XV/.fdf/.out used by test/6-utils/5-ani2traj) -- too short
# for MSD/VDOS to carry real physical meaning, but enough to confirm every
# report section/data file is generated without error and the numbers
# aren't NaN. The 2-atom O2 composition also makes --track-atom 0 /
# --track-pair 0-1 a physically meaningful demo (the O-O bond stretch).
export MPLBACKEND=Agg
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
echo "--- Starting tester for stb-aimdAnalysis (item 3.18) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR"/aimd.ANI "$FIXTURE_DIR"/aimd.fdf "$FIXTURE_DIR"/aimd.out "$FIXTURE_DIR"/aimd.XV "$TEST_DIR/"
cp "$FIXTURE_DIR"/siesta.ANI "$FIXTURE_DIR"/siesta.XV "$FIXTURE_DIR"/siesta.MDE \
   "$FIXTURE_DIR"/calc.fdf "$FIXTURE_DIR"/structure.fdf "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Basic run: numbered report, no flags -> only references.bib ---
echo -e "\n--- Testing a basic run (--label aimd), numbered report ---"

stb-aimdAnalysis --label aimd --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "Lattice read per-frame from aimd.out" log_basic.txt
check_contains "Composition" log_basic.txt
check_contains "O2" log_basic.txt
check_contains "RDF computed over 5 frame(s)" log_basic.txt
check_contains "Diffusion coefficient\|D (cm" log_basic.txt
check_contains "VDOS computed" log_basic.txt
check_contains "\[0\] RUN METADATA" log_basic.txt
check_contains "\[1\] INPUT DATA" log_basic.txt
check_contains "\[2\] RADIAL DISTRIBUTION FUNCTION" log_basic.txt
check_contains "\[3\] MEAN-SQUARED DISPLACEMENT" log_basic.txt
check_contains "\[4\] VELOCITY AUTOCORRELATION" log_basic.txt
check_contains "\[5\] SINGLE-ATOM DISPLACEMENT" log_basic.txt
check_contains "\[6\] ATOM-PAIR RELATIVE DISTANCE" log_basic.txt
check_contains "\[7\] THERMODYNAMIC TIME SERIES" log_basic.txt
check_contains "\[8\] WRITING OUTPUT FILES" log_basic.txt
check_contains "\[9\] REFERENCES" log_basic.txt
check_contains "\[10\] SUMMARY & FILES" log_basic.txt

echo "Testing: [7] reports Volume always (from the cell), even with no .MDE present"
check_contains "Volume (Ang^3)" log_basic.txt

echo "Testing: a plain run with no --save-gnuplot/--save-report/--track-* only writes references.bib"
if [ ! -f aimd_rdf.dat ] && [ ! -f aimd_rdf.gplot ] && [ ! -f stb_aimdAnalysis_report.txt ] && [ -f references.bib ]; then
    echo -e "   -> ${GREEN}Verified:${NC} no data/gplot/report written without --save-gnuplot/--save-report"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} data/gplot/report was written despite being off by default"
    FAIL=$((FAIL+1))
fi
rm -f references.bib


# --- 2b. --geometry-file (real bug fix: <label>.fdf is never auto-found in practice) ---
echo -e "\n--- Testing --geometry-file (the real .fdf is NOT named <label>.fdf) ---"
mkdir -p geomfile_test
cp aimd.ANI aimd.out aimd.XV geomfile_test/
cp aimd.fdf geomfile_test/calc.fdf
(cd geomfile_test && stb-aimdAnalysis --label aimd --no-intro > ../log_geom_missing.txt 2>&1)
check_contains "MD.LengthTimeStep not found" log_geom_missing.txt
check_contains "pass --geometry-file" log_geom_missing.txt

(cd geomfile_test && stb-aimdAnalysis --label aimd --geometry-file calc.fdf --no-intro > ../log_geom_fixed.txt 2>&1)
check_exit_code $? 0
if grep -q "MD.LengthTimeStep not found" log_geom_fixed.txt; then
    echo -e "   -> ${RED}Failed:${NC} --geometry-file calc.fdf did not resolve the MD timestep warning"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} --geometry-file calc.fdf resolves the MD timestep (no more warning)"
    PASS=$((PASS+1))
fi
check_contains "Geometry file     : calc.fdf" log_geom_fixed.txt
rm -rf geomfile_test

echo "Testing: --geometry-file rejected together with --trajectory"
stb-aimdAnalysis --trajectory aimd.ANI --geometry-file aimd.fdf --no-intro > log_geom_traj.txt 2>&1
check_exit_code $? 2
check_contains "only applies to --label mode" log_geom_traj.txt

echo "Testing: --geometry-file pointing to a nonexistent file"
stb-aimdAnalysis --label aimd --geometry-file nope.fdf --no-intro > log_geom_missing2.txt 2>&1
check_exit_code $? 2
check_contains "not found" log_geom_missing2.txt


# --- 2c. Real 500-step SIESTA Nose (NVT) AIMD run: --geometry-file + [7] THERMODYNAMIC TIME SERIES ---
# siesta.ANI/.XV/.MDE + calc.fdf/structure.fdf are a REAL SIESTA run (8-atom
# SiC supercell, Nose thermostat, target 500 K, 500 MD steps) -- SystemLabel
# 'siesta' but the real input is 'calc.fdf', reproducing the exact --geometry-
# file bug scenario live, plus a real siesta.MDE for the new [7] section.
echo -e "\n--- Testing a real 500-step SIESTA NVT run: --geometry-file + [7] THERMODYNAMIC TIME SERIES ---"
mkdir -p real_nvt
cp siesta.ANI siesta.XV siesta.MDE calc.fdf structure.fdf real_nvt/

(cd real_nvt && stb-aimdAnalysis --label siesta --geometry-file calc.fdf --save-gnuplot \
    --no-intro > ../log_real_nvt.txt 2>&1)
check_exit_code $? 0
check_contains "\[7\] THERMODYNAMIC TIME SERIES" log_real_nvt.txt
check_contains "Geometry file     : calc.fdf" log_real_nvt.txt

echo "Testing: Volume is exactly constant (a fixed-cell NVT run)"
python3 -c "
import numpy as np, sys
v = np.loadtxt('real_nvt/siesta_volume.dat')[:, 1]
# tolerance, not exact 0.0 -- numpy's std() of bit-identical values can still
# report a tiny ~1e-13 nonzero result due to internal floating-point rounding
sys.exit(0 if v.std() < 1e-6 and abs(v.mean() - 691.4382) < 1e-3 else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} volume is exactly constant at 691.4382 Ang^3"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} volume is not constant as expected for a fixed-cell NVT run"
    FAIL=$((FAIL+1))
fi

echo "Testing: energy per atom is also reported"
check_contains "Energy, total (eV/atom)" log_real_nvt.txt

echo "Testing: siesta_energy.dat's PLOTTED columns (2/3) are eV/atom, not the raw absolute total"
python3 -c "
import numpy as np, sys
with open('real_nvt/siesta_energy.dat') as f:
    header = f.readline()
d = np.loadtxt('real_nvt/siesta_energy.dat')
e_tot_per_atom, e_pot_per_atom = d[:, 1], d[:, 2]
e_tot_abs, e_pot_abs = d[:, 3], d[:, 4]
ok = ('eV_per_atom' in header and header.index('Energy_total(eV_per_atom)') < header.index('Energy_total(eV)')
      and np.allclose(e_tot_per_atom * 8, e_tot_abs) and np.allclose(e_pot_per_atom * 8, e_pot_abs))
sys.exit(0 if ok else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} per-atom columns come first and match abs/natoms (8 atoms)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} siesta_energy.dat column layout/values are wrong"
    FAIL=$((FAIL+1))
fi

echo "Testing: siesta_energy.gplot and siesta_thermo.gplot both plot in eV/atom"
check_contains 'set ylabel "Energy (eV/atom)"' real_nvt/siesta_energy.gplot
grep -A2 '"Energy -- siesta"' real_nvt/siesta_thermo.gplot | grep -q 'using 1:2'
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} siesta_thermo.gplot's Energy panel plots column 2 (eV/atom)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} siesta_thermo.gplot's Energy panel doesn't plot column 2"
    FAIL=$((FAIL+1))
fi

echo "Testing: Temperature fluctuates around the 500 K Nose target"
python3 -c "
import numpy as np, sys
t = np.loadtxt('real_nvt/siesta_temperature.dat')[:, 1]
sys.exit(0 if 400 <= t.mean() <= 600 else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} mean temperature lands near the 500 K Nose target"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} mean temperature far from the 500 K target"
    FAIL=$((FAIL+1))
fi

echo "Testing: E_tot fluctuates far less than E_KS (E_tot is the physically conserved-ish quantity)"
python3 -c "
import numpy as np, sys
e_tot, e_ks = np.loadtxt('real_nvt/siesta_energy.dat', usecols=(1, 2), unpack=True)
sys.exit(0 if e_tot.std() < e_ks.std() / 10 else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} E_tot standard deviation is at least 10x smaller than E_KS's"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} E_tot is not markedly more stable than E_KS"
    FAIL=$((FAIL+1))
fi

check_success real_nvt/siesta_volume.dat
check_success real_nvt/siesta_energy.dat
check_success real_nvt/siesta_energy.gplot
check_success real_nvt/siesta_temperature.dat
check_success real_nvt/siesta_pressure.dat
check_success real_nvt/siesta_thermo.gplot

echo "Testing: the combined 4-panel gnuplot script actually renders with a real gnuplot"
if command -v gnuplot >/dev/null 2>&1; then
    (cd real_nvt && gnuplot siesta_thermo.gplot 2>/dev/null)
    check_success real_nvt/siesta_thermo.pdf
else
    echo -e "   -> ${YELLOW}Skipped:${NC} gnuplot not installed in this environment"
fi
rm -rf real_nvt


# --- 3. --track-atom / --track-pair (new features) ---
echo -e "\n--- Testing --track-atom 0 --track-pair 0-1 (new: displacement + relative distance) ---"

stb-aimdAnalysis --label aimd --track-atom 0 --track-pair 0-1 --save-gnuplot --no-intro > log_track.txt 2>&1
check_exit_code $? 0
check_contains "Tracked atom      : index 0 (species O)" log_track.txt
check_contains "Tracked pair      : atom 0 (O) -- atom 1 (O)" log_track.txt
check_success aimd_disp_atom0.dat
check_success aimd_disp_atom0.gplot
check_success aimd_dist_0_1.dat
check_success aimd_dist_0_1.gplot

echo "Testing: the tracked pair's minimum-image distance matches the RDF's own first peak (O-O bond, ~1.0-1.3 Ang)"
python3 -c "
import numpy as np, sys
dist = np.loadtxt('aimd_dist_0_1.dat')[:, 4]
sys.exit(0 if np.all((dist >= 1.0) & (dist <= 1.3)) else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} atom-pair distance stays in the expected O-O bond-length range"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} atom-pair distance outside the expected range"
    FAIL=$((FAIL+1))
fi

echo "Testing: single-atom displacement magnitude is finite and starts at 0"
python3 -c "
import numpy as np, sys
data = np.loadtxt('aimd_disp_atom0.dat')
disp_mag = data[:, 4]
sys.exit(0 if disp_mag[0] == 0.0 and np.all(np.isfinite(disp_mag)) else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} displacement magnitude starts at 0 and stays finite"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} displacement magnitude data looks wrong"
    FAIL=$((FAIL+1))
fi

echo "Testing: the saved .gplot scripts actually render with a real gnuplot"
if command -v gnuplot >/dev/null 2>&1; then
    gnuplot aimd_disp_atom0.gplot && check_success aimd_disp_atom0.pdf
    gnuplot aimd_dist_0_1.gplot && check_success aimd_dist_0_1.pdf
else
    echo -e "   -> ${YELLOW}Skipped:${NC} gnuplot not installed in this environment"
fi
rm -f aimd_*.dat aimd_*.gplot aimd_*.pdf references.bib


# --- 4. --track-atom / --track-pair validation ---
echo -e "\n--- Testing --track-atom out of range (expects failure) ---"
stb-aimdAnalysis --label aimd --track-atom 5 --no-intro > log_badatom.txt 2>&1
check_exit_code $? 1
check_contains "out of range" log_badatom.txt

echo -e "\n--- Testing --track-pair with the same index twice (expects failure) ---"
stb-aimdAnalysis --label aimd --track-pair 0-0 --no-intro > log_badpair2.txt 2>&1
check_exit_code $? 2
check_contains "DIFFERENT atom indices" log_badpair2.txt

echo -e "\n--- Testing --track-pair malformed (expects failure) ---"
stb-aimdAnalysis --label aimd --track-pair 0 --no-intro > log_malformedpair.txt 2>&1
check_exit_code $? 2
check_contains "must be 'I-J'" log_malformedpair.txt


# --- 5. --pair (same-species RDF) ---
echo -e "\n--- Testing --pair O-O ---"
stb-aimdAnalysis --label aimd --pair O-O --save-gnuplot --no-intro > log_pair.txt 2>&1
check_exit_code $? 0
check_contains "RDF pair          : O-O" log_pair.txt
check_success aimd_rdf_OO.dat
rm -f aimd_*.dat aimd_*.gplot references.bib


# --- 6. --pair with a species not present ---
echo -e "\n--- Testing --pair with a nonexistent species (expects failure) ---"
stb-aimdAnalysis --label aimd --pair O-H --no-intro > log_badpair.txt 2>&1
check_exit_code $? 1
check_contains "species not found" log_badpair.txt


# --- 7. --stride / --skip ---
echo -e "\n--- Testing --stride 2 --skip 1 ---"
stb-aimdAnalysis --label aimd --stride 2 --skip 1 --no-intro > log_stride.txt 2>&1
check_exit_code $? 0
check_contains "stride 2, skip 1" log_stride.txt


# --- 8. --skip too large ---
echo -e "\n--- Testing --skip larger than the available frames (expects failure) ---"
stb-aimdAnalysis --label aimd --skip 999 --no-intro > log_badskip.txt 2>&1
check_exit_code $? 1
check_contains "discards all" log_badskip.txt


# --- 9. --save-report / -o output-dir ---
echo -e "\n--- Testing --save-report -o custom_out ---"
stb-aimdAnalysis --label aimd --save-report -o custom_out --no-intro > log_savereport.txt 2>&1
check_exit_code $? 0
check_success custom_out/stb_aimdAnalysis_report.txt
check_contains "STB-AIMDANALYSIS REPORT" custom_out/stb_aimdAnalysis_report.txt
rm -rf custom_out


# --- 10. --view (should not block with MPLBACKEND=Agg) ---
echo -e "\n--- Testing --view (matplotlib preview, Agg backend) ---"
timeout 30 stb-aimdAnalysis --label aimd --track-atom 0 --track-pair 0-1 --view --no-intro > log_view.txt 2>&1
check_exit_code $? 0
rm -f references.bib


# --- 11. --trajectory: generic ASE-readable input (e.g. stb-mlmd's own output) ---
echo -e "\n--- Testing --trajectory (generic ASE trajectory, independent of SIESTA/.ANI) ---"
python3 -c "
from ase import Atoms
from ase.io import write
import numpy as np
rng = np.random.default_rng(0)
cell = [[10.0, 0, 0], [0, 10.0, 0], [0, 0, 10.0]]
frames = []
for i in range(6):
    pos = np.array([[5.0, 5.0, 4.395], [5.0, 5.0, 5.605]]) + rng.normal(scale=0.05, size=(2, 3))
    a = Atoms(symbols=['O', 'O'], positions=pos, cell=cell, pbc=True)
    a.info['Time'] = i * 15.0
    frames.append(a)
write('synthetic_traj.xyz', frames, format='extxyz')
"
check_success synthetic_traj.xyz

echo "Testing: auto-detected dt from embedded 'Time' info (no --dt needed)"
stb-aimdAnalysis --trajectory synthetic_traj.xyz --track-atom 0 --track-pair 0-1 --no-intro > log_trajectory.txt 2>&1
check_exit_code $? 0
check_contains "Auto-detected dt = 15.0000 fs" log_trajectory.txt
check_contains "No SIESTA-specific references for a generic --trajectory input" log_trajectory.txt

echo "Testing: --trajectory + --label are mutually exclusive"
stb-aimdAnalysis --trajectory synthetic_traj.xyz --label aimd --no-intro > log_mutex.txt 2>&1
check_exit_code $? 2
check_contains "mutually exclusive" log_mutex.txt

echo "Testing: --dt overrides/is required for a format without embedded Time"
python3 -c "
from ase.io import read, write
frames = read('synthetic_traj.xyz', index=':')
write('synthetic_traj.xsf', frames, format='xsf')
"
stb-aimdAnalysis --trajectory synthetic_traj.xsf --no-intro > log_missing_dt.txt 2>&1
check_exit_code $? 1
check_contains "\-\-dt is required" log_missing_dt.txt
stb-aimdAnalysis --trajectory synthetic_traj.xsf --dt 15 --no-intro > log_with_dt.txt 2>&1
check_exit_code $? 0


# --- 12. Robustness / errors ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: missing --label and --trajectory"
stb-aimdAnalysis --no-intro > log_missing_label.txt 2>&1
check_exit_code $? 2
check_contains "one of --label or --trajectory is required" log_missing_label.txt

echo "Testing: missing .ANI file"
stb-aimdAnalysis --label nope --no-intro > log_missing_ani.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_ani.txt

echo "Testing: --version"
stb-aimdAnalysis --version > log_version.txt 2>&1
check_contains "stb-aimdAnalysis" log_version.txt

echo "Testing: --help documents --pair, --skip, --fit-start, --track-atom, --track-pair, --geometry-file, --save-gnuplot, --view, --trajectory, --dt"
stb-aimdAnalysis --help > log_help.txt 2>&1
check_contains "pair" log_help.txt
check_contains "skip" log_help.txt
check_contains "fit-start" log_help.txt
check_contains "track-atom" log_help.txt
check_contains "track-pair" log_help.txt
check_contains "list-atoms" log_help.txt
check_contains "geometry-file" log_help.txt
check_contains "save-gnuplot" log_help.txt
check_contains "view" log_help.txt
check_contains "trajectory" log_help.txt
check_contains "\-\-dt" log_help.txt


# --- 12b. --list-atoms (index/species/coordinates preview, early exit) ---
echo -e "\n--- Testing --list-atoms (early-exit atom index/species/coordinates preview) ---"
stb-aimdAnalysis --label aimd --list-atoms --no-intro > log_list_atoms.txt 2>&1
check_exit_code $? 0
check_contains "2 atom(s), coordinates from the first frame" log_list_atoms.txt
check_contains "Index | Species | X" log_list_atoms.txt
check_contains "0     | O" log_list_atoms.txt
check_contains "1     | O" log_list_atoms.txt

echo "Testing: --list-atoms exits before running RDF/MSD/VDOS/anything else (no report sections)"
if grep -q "RADIAL DISTRIBUTION\|MEAN-SQUARED\|VELOCITY AUTOCORRELATION" log_list_atoms.txt; then
    echo -e "   -> ${RED}Failed:${NC} --list-atoms ran the full analysis instead of exiting early"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} --list-atoms exited before any RDF/MSD/VDOS section"
    PASS=$((PASS+1))
fi

echo "Testing: --list-atoms works for --trajectory input too"
stb-aimdAnalysis --trajectory synthetic_traj.xyz --list-atoms --no-intro > log_list_atoms_traj.txt 2>&1
check_exit_code $? 0
check_contains "2 atom(s), coordinates from the first frame" log_list_atoms_traj.txt


# --- 13. Interactive path (stb-suite, shortcut 3.18) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 3.18) ---"
rm -f references.bib aimd_disp_atom0.dat aimd_dist_0_1.dat
printf '3.18\n\naimd\n\n1\n0\n\nn\n0\n0-1\n\nn\ny\nn\n\n0\n' | timeout 60 stb-suite > log_menu.txt 2>&1
check_contains "Tracked atom      : index 0" log_menu.txt
check_contains "Tracked pair      : atom 0 (O) -- atom 1 (O)" log_menu.txt
check_success aimd_disp_atom0.dat
check_success aimd_dist_0_1.dat

echo "Testing: interactive menu asks before listing atoms -- 'y' shows the table, 'n' skips it"
rm -f references.bib aimd_disp_atom0.dat aimd_dist_0_1.dat
printf '3.18\n\naimd\n\n1\n0\n\ny\n0\n0-1\n\nn\nn\nn\n\n0\n' | timeout 60 stb-suite > log_menu_list.txt 2>&1
check_contains "List every atom's index/species/coordinates" log_menu_list.txt
check_contains "2 atom(s), coordinates from the first frame" log_menu_list.txt
check_contains "0     | O" log_menu_list.txt


popd > /dev/null

# --- 14. Summary ---
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
