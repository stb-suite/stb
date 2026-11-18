#!/bin/bash

# --- Setup ---
# Smoke test for stb-convergenceAnalysis (Convergence Test Analysis, item 5.5.2)
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
# The 4 committed 'convergence_meshcutoff_*' fixtures carry both a
# 'siesta: FreeEng' line AND a real 'outcell: Unit cell vectors' block, with
# energy and relaxed-cell-volume trends that DELIBERATELY disagree: energy
# looks converged by meshcutoff=250, but the relaxed cell volume only
# stabilizes at 300 -- the exact physical scenario stb-convergence's own
# module docstring warns about (energy converging before the relaxed
# geometry does), used below to exercise the dual energy+structure
# convergence criteria and their cross-check note.
echo "--- Starting tester for STB-ConvergenceAnalysis (item 5.5.2) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR/convergence_runs/meshcutoff"
for d in 150.0000 200.0000 250.0000 300.0000; do
    cp -r "$FIXTURE_DIR/convergence_meshcutoff_$d" "$TEST_DIR/convergence_runs/meshcutoff/"
done
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null
STRUCT="convergence_runs/meshcutoff/convergence_meshcutoff_150.0000/structure.fdf"


# --- 2. Basic nested-layout analysis: energy vs. structure disagreement ---
echo -e "\n--- Testing nested-layout analysis of a 4-point meshcutoff sweep ---"
stb-convergenceAnalysis --dir convergence_runs --no-intro > log_basic.txt 2>&1
check_contains "Folders found: 4" log_basic.txt
check_contains "Energy converged at:.*meshcutoff = 250.0000" log_basic.txt
check_contains "Structure converged at:.*meshcutoff = 300.0000" log_basic.txt
check_contains "Energy and relaxed structure disagree" log_basic.txt
if grep -q "NOT CONVERGED" log_basic.txt; then
    echo -e "   -> ${RED}Failed:${NC} unexpected NOT CONVERGED for a resolved sweep"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no NOT CONVERGED (recommended value was found)"
    PASS=$((PASS+1))
fi


# --- 3. Tolerance overrides ---
echo -e "\n--- Testing a very tight energy tolerance (energy WARNING, structure still resolves) ---"
stb-convergenceAnalysis --dir convergence_runs --tolerance 0.0000001 --no-intro > log_tight_energy.txt 2>&1
check_contains "No point met the 1e-07 eV/atom energy tolerance" log_tight_energy.txt
check_contains "Structure converged at:.*meshcutoff = 300.0000" log_tight_energy.txt

echo -e "\n--- Testing a very tight volume tolerance (structure WARNING) ---"
stb-convergenceAnalysis --dir convergence_runs --volume-tolerance 0.0000001 --no-intro > log_tight_volume.txt 2>&1
check_contains "No point met the 1e-07% relaxed-volume tolerance" log_tight_volume.txt

echo -e "\n--- Testing a volume tolerance wide enough that energy and structure AGREE ---"
stb-convergenceAnalysis --dir convergence_runs --volume-tolerance 2.0 --no-intro > log_agree.txt 2>&1
check_contains "Structure converged at:.*meshcutoff = 250.0000" log_agree.txt
check_contains "Energy and structure agree on the same converged value." log_agree.txt


# --- 4. Data-quality flags (SCF not confirmed, residual force, step-limited) ---
echo -e "\n--- Testing data-quality flags: SCF?, F>tol, steps ---"
mkdir -p quality_runs/meshcutoff/convergence_meshcutoff_100.0000
cp "$STRUCT" quality_runs/meshcutoff/convergence_meshcutoff_100.0000/
cat > quality_runs/meshcutoff/convergence_meshcutoff_100.0000/calc.out <<'EOF'
siesta: FreeEng =    -500.900000

outcell: Unit cell vectors (Ang):
        5.35    0.000000    0.000000
        0.000000    5.35    0.000000
        0.000000    0.000000    5.35
EOF

mkdir -p quality_runs/meshcutoff/convergence_meshcutoff_200.0000
cp "$STRUCT" quality_runs/meshcutoff/convergence_meshcutoff_200.0000/
cat > quality_runs/meshcutoff/convergence_meshcutoff_200.0000/calc.out <<'EOF'
SCF cycle converged after 10 iterations
siesta: FreeEng =    -500.850000

outcell: Unit cell vectors (Ang):
        5.42    0.000000    0.000000
        0.000000    5.42    0.000000
        0.000000    0.000000    5.42

siesta: Atomic forces (eV/Ang):
     1    0.500000    0.500000    0.500000
   Max    0.500000
EOF

mkdir -p quality_runs/meshcutoff/convergence_meshcutoff_300.0000
cp "$STRUCT" quality_runs/meshcutoff/convergence_meshcutoff_300.0000/
cat > quality_runs/meshcutoff/convergence_meshcutoff_300.0000/calc.out <<'EOF'
SCF cycle converged after 10 iterations
siesta: FreeEng =    -500.845000

outcell: Unit cell vectors (Ang):
        5.45    0.000000    0.000000
        0.000000    5.45    0.000000
        0.000000    0.000000    5.45

                        Begin CG opt. move =      0
                        Begin CG opt. move =      1
                        Begin CG opt. move =      2
EOF
printf 'MD.NumCGsteps      3\n' > quality_runs/meshcutoff/convergence_meshcutoff_300.0000/config_extra.fdf

stb-convergenceAnalysis --dir quality_runs --no-intro > log_quality.txt 2>&1
check_contains "SCF?" log_quality.txt
check_contains "F>tol" log_quality.txt
check_contains "steps" log_quality.txt
check_contains "3/3 run(s) flagged" log_quality.txt

echo -e "\n--- Testing --force-tolerance override clears the F>tol flag ---"
stb-convergenceAnalysis --dir quality_runs --no-intro --force-tolerance 10 > log_force_tol.txt 2>&1
check_contains "200.0000 | -250.425000 | 0.025000       | 159.2201      | 3.9768         | OK" log_force_tol.txt


# --- 5. Flat layout (pointing --dir directly at one parameter's own subfolder) ---
echo -e "\n--- Testing flat layout (--dir pointed directly at the meshcutoff subfolder) ---"
stb-convergenceAnalysis --dir convergence_runs/meshcutoff --no-intro > log_flat.txt 2>&1
check_contains "Energy converged at:.*meshcutoff = 250.0000" log_flat.txt


# --- 6. Multi-parameter nested discovery (all 3 swept parameters in one --dir) ---
echo -e "\n--- Testing multi-parameter discovery (meshcutoff + energyshift + kgrid together) ---"
mkdir -p convergence_runs/energyshift/convergence_energyshift_0.0050
cp "$STRUCT" convergence_runs/energyshift/convergence_energyshift_0.0050/
printf 'SCF cycle converged after 10 iterations\nsiesta: FreeEng =    -500.850000\n' \
    > convergence_runs/energyshift/convergence_energyshift_0.0050/calc.out
mkdir -p convergence_runs/energyshift/convergence_energyshift_0.0200
cp "$STRUCT" convergence_runs/energyshift/convergence_energyshift_0.0200/
printf 'SCF cycle converged after 9 iterations\nsiesta: FreeEng =    -500.850050\n' \
    > convergence_runs/energyshift/convergence_energyshift_0.0200/calc.out

mkdir -p convergence_runs/kgrid/convergence_kgrid_0.1000
cp "$STRUCT" convergence_runs/kgrid/convergence_kgrid_0.1000/
printf 'SCF cycle converged after 10 iterations\nsiesta: FreeEng =    -500.850000\n' \
    > convergence_runs/kgrid/convergence_kgrid_0.1000/calc.out
mkdir -p convergence_runs/kgrid/convergence_kgrid_0.2000
cp "$STRUCT" convergence_runs/kgrid/convergence_kgrid_0.2000/
printf 'SCF cycle converged after 9 iterations\nsiesta: FreeEng =    -500.850050\n' \
    > convergence_runs/kgrid/convergence_kgrid_0.2000/calc.out

stb-convergenceAnalysis --dir convergence_runs --no-intro > log_multi.txt 2>&1
check_contains "\[2.1\] MESHCUTOFF SWEEP" log_multi.txt
check_contains "\[2.2\] ENERGYSHIFT SWEEP" log_multi.txt
check_contains "\[2.3\] KGRID SWEEP" log_multi.txt
# EnergyShift is judged by the relaxed volume only (no plateau in E); these logs carry no cell
check_contains "PAO.EnergyShift is NOT judged by the total energy" log_multi.txt
check_contains "Energy converged at:.*kgrid = 0.1000" log_multi.txt


# --- 7. --save-report / --save-gnuplot / --output-dir ---
echo -e "\n--- Testing --save-report/--save-gnuplot/--output-dir ---"
rm -rf saved_out
stb-convergenceAnalysis --dir convergence_runs --no-intro --save-report --save-gnuplot \
    --output-dir saved_out > log_saved.txt 2>&1
check_success saved_out/convergence_report.txt
check_success saved_out/meshcutoff_convergence.dat
check_success saved_out/meshcutoff_convergence.gplot
check_success saved_out/energyshift_convergence.dat
check_success saved_out/kgrid_convergence.dat
check_contains "Energy converged at: meshcutoff = 250.0000" saved_out/convergence_report.txt

echo "Testing: the .gplot script actually renders with gnuplot"
if command -v gnuplot > /dev/null 2>&1; then
    ( cd saved_out && gnuplot meshcutoff_convergence.gplot )
    check_success saved_out/meshcutoff_convergence.pdf
else
    echo -e "   -> ${YELLOW}Skipped:${NC} gnuplot not installed"
fi


# --- 8. --apply (multi-parameter nested tree, kgrid via substitute_kgrid_tag) ---
echo -e "\n--- Testing --apply writes every converged parameter into one calc.fdf ---"
printf 'Mesh.CutOff   150.0000 Ry\nPAO.EnergyShift 0.0500 Ry\nkgrid.MonkhorstPack   [1  1  1]\n' \
    > apply_target.fdf
stb-convergenceAnalysis --dir convergence_runs --apply apply_target.fdf --no-intro > log_apply.txt 2>&1
check_contains "\[Applied\].*meshcutoff = 300.0000" log_apply.txt
check_contains "Skipped energyshift: no converged value to apply" log_apply.txt
check_contains "\[Applied\].*kgrid = 0.1000" log_apply.txt
check_contains "Mesh.CutOff   300.0000 Ry" apply_target.fdf
check_contains "PAO.EnergyShift 0.0500 Ry" apply_target.fdf
# Bulk Si, 5.43 Ang cubic cell (from convergence_meshcutoff_150.0000/structure.fdf),
# no vacuum: density 0.1000 -> ceil(reciprocal length / 0.1) = 12 divisions/axis.
check_contains "kgrid.MonkhorstPack   \[12  12  12\]" apply_target.fdf

echo -e "\n--- Testing --apply against a read-only target fails cleanly (no traceback) ---"
printf 'Mesh.CutOff   150.0000 Ry\n' > readonly_target.fdf
chmod 444 readonly_target.fdf
stb-convergenceAnalysis --dir convergence_runs --apply readonly_target.fdf --no-intro > log_apply_readonly.txt 2>&1
check_contains "\[ERROR\] Could not write" log_apply_readonly.txt
chmod 644 readonly_target.fdf


# --- 9. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: empty directory (no convergence_* folders)"
mkdir -p empty_runs
stb-convergenceAnalysis --dir empty_runs --no-intro > log_empty.txt 2>&1
check_exit_code $? 1
check_contains "No 'convergence_\*' folders found" log_empty.txt

echo "Testing: missing directory"
stb-convergenceAnalysis --dir does_not_exist --no-intro > log_missing_dir.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_dir.txt

echo "Testing: folder missing calc.out is skipped, not a hard failure (partial sweep)"
mkdir -p bad_runs/meshcutoff/convergence_meshcutoff_100.0000
cp "$STRUCT" bad_runs/meshcutoff/convergence_meshcutoff_100.0000/
stb-convergenceAnalysis --dir bad_runs --no-intro > log_bad_runs.txt 2>&1
check_exit_code $? 0
check_contains "No valid data found for this parameter" log_bad_runs.txt
check_contains "NOT CONVERGED" log_bad_runs.txt

echo "Testing: --version"
stb-convergenceAnalysis --version > log_version.txt 2>&1
check_contains "stb-convergenceAnalysis" log_version.txt

echo "Testing: --help documents tolerance/volume-tolerance/force-tolerance/save-gnuplot/view/apply"
stb-convergenceAnalysis --help > log_help.txt 2>&1
check_contains "tolerance" log_help.txt
check_contains "volume-tolerance" log_help.txt
check_contains "force-tolerance" log_help.txt
check_contains "save-gnuplot" log_help.txt
check_contains "\-\-view" log_help.txt
check_contains "apply" log_help.txt


# --- 10. Interactive path (stb-suite, shortcut 5.5.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 5.5.2) ---"

echo "Testing: navigate 5.5.2 -> all defaults (no advanced settings, no save/view/apply) -> quit"
printf '5.5.2\n\n\n\nn\n\nn\nn\nn\n\n\n0\n' | stb-suite > log_menu_defaults.txt 2>&1
check_contains "Energy converged at:.*meshcutoff = 250.0000" log_menu_defaults.txt

echo "Testing: navigate 5.5.2 -> advanced settings + save-report/gnuplot + apply -> quit"
rm -rf menu_out
printf 'Mesh.CutOff   150.0000 Ry\n' > menu_apply_target.fdf
printf '5.5.2\n\n\n\ny\n2.0\n0.05\nmenu_out\ny\ny\nn\nmenu_apply_target.fdf\n\n\n0\n' \
    | stb-suite > log_menu_advanced.txt 2>&1
check_contains "Volume tolerance  : 2.0 %" log_menu_advanced.txt
check_contains "Energy and structure agree on the same converged value." log_menu_advanced.txt
check_success menu_out/convergence_report.txt
check_success menu_out/meshcutoff_convergence.dat
check_contains "Applied.*meshcutoff = 250.0000" log_menu_advanced.txt
check_contains "Mesh.CutOff   250.0000 Ry" menu_apply_target.fdf


# --- 10b. Relaxations that never moved the cell (SIESTA's 1 GPa default stress tolerance) ---
echo -e "\n--- Testing: relaxation stopped at the input cell is flagged (P tolerance 1 GPa, V=input) ---"
rm -rf es_relax
CELL=$(awk '/%block LatticeVectors/{f=1;next} /%endblock LatticeVectors/{f=0} f' "$FIXTURE_DIR/convergence_meshcutoff_150.0000/structure.fdf")
for pair in 0.0200:-500.90 0.0100:-500.95 0.0050:-500.97; do
    v=${pair%%:*}; e=${pair##*:}
    d=es_relax/energyshift/convergence_energyshift_$v
    mkdir -p $d
    cp "$FIXTURE_DIR/convergence_meshcutoff_150.0000/structure.fdf" $d/
    printf 'PAO.EnergyShift %s Ry\nMD.TypeOfRun CG\nMD.Steps 100\nMD.VariableCell true\n' $v > $d/config_extra.fdf
    { echo "Begin CG opt. move =      0"; echo "SCF cycle converged after 10 iterations"; echo "siesta: FreeEng =    $e";
      echo "outcell: Unit cell vectors (Ang):"; echo "$CELL";
      echo "siesta: Stress tensor (static) (eV/Ang**3):";
      echo "siesta:     -0.005990    0.000000    0.000000"; echo "siesta:      0.000000   -0.005990    0.000000";
      echo "siesta:      0.000000    0.000000   -0.005990"; } > $d/calc.out
done
stb-convergenceAnalysis --dir es_relax --no-intro > log_es_relax.txt 2>&1
check_contains "Ptol=1GPa" log_es_relax.txt
check_contains "V=input" log_es_relax.txt
check_contains "relaxation(s) ended with exactly the input cell" log_es_relax.txt

# --- 10b2. A variable-cell relaxation of an fcc cell that ends at 59.91 deg (stress noise):
#     R-3m at 0.001 Ang, still Fm-3m at 0.01 -> warning; a cell that keeps 60 deg -> none. ---
echo -e "\n--- Testing: relaxed cell that lost the input symmetry to noise is flagged ---"
rm -rf fcc_relax
python3 - <<'PYEOF'
import os
from ase.cell import Cell
def fdf_cell(c): return "\n".join("   ".join(f"{x:12.8f}" for x in r) for r in c)
ideal = Cell.fromcellpar([2.55, 2.55, 2.55, 60, 60, 60])[:]
for value, alpha, e in ((300, 59.91, -1500.90), (400, 59.92, -1500.95), (500, 59.92, -1500.96)):
    d = f"fcc_relax/meshcutoff/convergence_meshcutoff_{value:.4f}"
    os.makedirs(d, exist_ok=True)
    open(f"{d}/structure.fdf", "w").write(
        "NumberOfAtoms 1\nNumberOfSpecies 1\n%block ChemicalSpeciesLabel\n 1 29 Cu\n%endblock ChemicalSpeciesLabel\n"
        "LatticeConstant 1.0 Ang\n%block LatticeVectors\n" + fdf_cell(ideal) + "\n%endblock LatticeVectors\n"
        "AtomicCoordinatesFormat Fractional\n%block AtomicCoordinatesAndAtomicSpecies\n 0 0 0 1\n"
        "%endblock AtomicCoordinatesAndAtomicSpecies\n")
    open(f"{d}/config_extra.fdf", "w").write(f"Mesh.Cutoff {value} Ry\nMD.TypeOfRun CG\nMD.Steps 100\nMD.VariableCell true\n")
    relaxed = Cell.fromcellpar([2.55, 2.55, 2.55, alpha, alpha, alpha])[:]
    open(f"{d}/calc.out", "w").write(
        "Begin CG opt. move =      3\nSCF cycle converged after 10 iterations\n"
        f"siesta: FreeEng =    {e}\noutcell: Unit cell vectors (Ang):\n{fdf_cell(relaxed)}\n")
PYEOF
stb-convergenceAnalysis --dir fcc_relax --no-intro > log_fcc_relax.txt 2>&1
check_contains "relaxed cell(s) lost the input symmetry to numerical noise: Fm-3m (No. 225) in, R-3m (No. 166) out" log_fcc_relax.txt
echo "Testing: the earlier relaxations (cell unchanged) give no such warning"
if grep -q "lost the input symmetry" log_es_relax.txt; then
    echo -e "   -> ${RED}Failed:${NC} symmetry-loss warning for an unchanged cell"; FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no symmetry-loss warning for an unchanged cell"; PASS=$((PASS+1))
fi

# --- 10c. --parameter basis, 2D: area / a0 / N/m (N7), 'every tighter value' closing rule (N13),
#     SIESTA 5.4.2 filter lines (N5), -p in the next-step command (N11). Synthetic E(A) parabolas. ---
echo -e "\n--- Testing the basis analysis of a 2D range stage (synthetic, false plateau at 0.02-0.005 Ry) ---"
rm -rf basis2d
python3 - <<'PYEOF'
import json, os, numpy as np
c = 20.0
def cell(a): return np.array([[a, 0, 0], [-a / 2, a * np.sqrt(3) / 2, 0], [0, 0, c]])
A_ref = np.linalg.norm(np.cross(cell(2.46)[0], cell(2.46)[1]))
k = 200.0 / 16.0218                                   # 200 N/m in eV/Ang^2
# area factor per EnergyShift: 0.02-0.005 agree with each other (a false plateau), the tighter ones +3 %
factors = {0.02: 1.000, 0.01: 1.004, 0.005: 1.002, 0.002: 1.030, 0.001: 1.031, 0.0005: 1.0295}
subs = []
for shift, f in factors.items():
    name = f"range_{shift:g}"
    subs.append({"name": name, "value": shift, "mesh_cutoff": 500, "energyshift": shift, "size": "DZP"})
    A0 = A_ref * f
    for st in np.linspace(-8, 8, 9):
        A = A_ref * (1 + st / 100); a = 2.46 * np.sqrt(A / A_ref)
        E = -310.0 + 0.5 * k / A0 * (A - A0) ** 2
        d = f"basis2d/basis/{name}/vol_{'m' if st < 0 else ''}{abs(st):.2f}"
        os.makedirs(d, exist_ok=True)
        lines = ["Filter: Requested cutoff:  150.000000 Ry", "PAO Filter: Cutoff used:   73.50 Ry",
                 "outcell: Unit cell vectors (Ang):"] + ["   ".join(f"{x:12.6f}" for x in r) for r in cell(a)]
        lines += ["", "SCF cycle converged after 10 iterations", f"siesta: Etot    =    {E:.8f}",
                  f"siesta: E_KS(eV) =    {E:.8f}"]
        open(f"{d}/calc.out", "w").write("\n".join(lines) + "\n")
json.dump({"format": "stb-basis-scan", "version": 1, "tool": "test", "stage": "range", "structure": "structure.fdf",
           "calc": "calc.fdf", "pseudo": "dojo", "n_atoms": 2, "volume": A_ref * c, "n_periodic": 2,
           "periodic_axes": [0, 1], "measure": "area", "scan_range_pct": 8.0, "scan_points": 9, "sub_scans": subs},
          open("basis2d/basis/basis_scan_range.json", "w"), indent=1)
PYEOF
stb-convergenceAnalysis --dir basis2d --file calc.out --no-intro > log_basis2d.txt 2>&1
check_contains "A0 (Ang^2) | a0 (Ang)" log_basis2d.txt
check_contains "stiffness (N/m)" log_basis2d.txt
check_contains "0.0005 Ry  *| 9/9    | 5.3940     | 2.4957" log_basis2d.txt
check_contains "VERDICT: PAO.EnergyShift = 0.002 Ry" log_basis2d.txt
check_contains "whose area is within 0.3 % of every tighter setting" log_basis2d.txt
check_contains "the logs show the orbital filter" log_basis2d.txt
check_contains "Next: stb-convergence --parameter basis --basis-stage size --energyshift 0.002 --mesh-cutoff 500 -s structure.fdf -c calc.fdf -p dojo" log_basis2d.txt
if grep -q "0.0005 Ry .*| 20[0-9]\.[0-9]" log_basis2d.txt; then
    echo -e "   -> ${GREEN}Verified:${NC} in-plane stiffness ~200 N/m"; PASS=$((PASS+1))
else echo -e "   -> ${RED}Failed:${NC} in-plane stiffness not ~200 N/m"; FAIL=$((FAIL+1)); fi

popd > /dev/null

# --- 11. Summary ---
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
