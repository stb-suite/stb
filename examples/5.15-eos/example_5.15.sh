#!/bin/bash
# Guided example: Workflow 5.15 -- Equation of State (both stages: stb-eosInputs, code 5.15.1, and
# stb-eosAnalysis, code 5.15.2, in the stb-suite menu).
#
# Not an automated test (see test/5-workflow/15-eos/{prep,analysis}/test.sh for that) -- a commented
# walk-through: it runs real commands, one case at a time, into its own output/<case>/ folder, and
# shows the piece of output that proves what just happened. Pauses between sections so you can read
# before moving on. Safe to re-run any time -- it always starts by wiping its own output/.
#
# Stage 1 (stb-eosInputs) runs for real in every case below (cheap -- it only writes .fdf files, no
# SCF). Stage 2 (stb-eosAnalysis) runs on REAL SIESTA results: reference_results/si/ holds the finished
# runs of bulk silicon -- the 11-volume scan of the base setup and 13 variants that change one thing at
# a time (k-grid, Mesh.Cutoff, basis size, PAO.EnergyShift) -- computed with the Stage 1 commands of
# this script and the inputs in this folder (SIESTA 5.4.2, mpirun -np 2; see README section 6 and
# make_reference.sh to regenerate them). Nothing in Stage 2 is fabricated; this script never runs SIESTA.

set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

# --show needs a display; MPLBACKEND=Agg keeps this script headless-safe (the figure is still saved
# with --png).
export MPLBACKEND=Agg

OUT="$DIR/output"
rm -rf "$OUT"
mkdir -p "$OUT"
REF="$DIR/reference_results/si"

pause() {
    echo
    read -p "  [Press Enter to continue] " -r
    echo
}

banner() {
    echo "=================================================================="
    echo " $1"
    echo "=================================================================="
}

strip_colors() { sed 's/\x1b\[[0-9;]*m//g'; }

# section N of a report: from "[N] TITLE" up to the next "[M] " header
section() {   # $1 = file, $2 = section number
    strip_colors < "$1" | awk -v n="[$2]" 'index($0, n) == 1 {on = 1; print; next}
        on && /^\[[0-9]+\] / {exit} on {print}'
}


banner "Welcome: the equation of state of a crystal, from a handful of SIESTA single points"
cat <<'EOF'
The energy of a crystal as a function of its volume, E(V), holds its equilibrium volume V0, its
cohesive-energy minimum E0, its bulk modulus B0 = V d2E/dV2 and the pressure derivative B0' -- and,
by differentiating, its pressure at any volume. SIESTA gives one point E(V) per run; this workflow
scales the cell isotropically, runs the scan, and fits an equation of state (Birch-Murnaghan, Vinet,
Murnaghan or stabilized jellium).

Stage 1 (stb-eosInputs) prepares one folder per volume: the scaled cell, YOUR calc.fdf untouched, a
config_extra.fdf with the directives the scan needs (fixed cell, single point), the
pseudopotentials. Stage 2 (stb-eosAnalysis) reads the finished runs and reports the fit AND whether
to trust it: SIESTA's own pressure (from the stress) against -dE/dV, a second B0 from the stress,
uncertainties, residuals, level-of-theory consistency, and a dashboard.

This tour: Stage 1 on silicon, a template with problems, a hexagonal crystal and a slab (Cases 1-5),
then Stage 2 on real SIESTA results for silicon (Cases 6-11): the fit, the pressure cross-check, the
outputs, the convergence study (k-points, mesh, basis) and the comparison with an all-electron
reference, which shows why B0 alone is not a good figure of merit for a basis set. See README.md for the
theory (section 1) and the worked example (section 5).
EOF
pause


banner "Case 1: Stage 1 -- bulk silicon, what stb-eosInputs reports"
cat <<'EOF'
structure.fdf is the 2-atom silicon cell (a = 5.469 Ang, fractional coordinates); calc.fdf is the
template: PBE, DZP, a k-grid, and %include structure.fdf. Eleven volumes, -8 % .. +8 %. -p dojo copies
the PseudoDojo pseudopotential into every folder. Watch [1] (structure, space group, valence
electrons), [2] (the level of theory read from the template) and [3] (the scan: volume, cell
lengths and shortest bond at every point).
EOF
mkdir -p "$OUT/case1-si"
cp structure.fdf calc.fdf "$OUT/case1-si/"
echo "\$ stb-eosInputs -f structure.fdf -p dojo --n-volumes 11 --strain-range 8.0 --save-report --no-intro"
(cd "$OUT/case1-si" && stb-eosInputs -f structure.fdf -p dojo --n-volumes 11 --strain-range 8.0 \
    --save-report --no-intro > stage1.log 2>&1)
section "$OUT/case1-si/stage1.log" 1 | cut -c1-150
section "$OUT/case1-si/stage1.log" 2 | cut -c1-150
section "$OUT/case1-si/stage1.log" 3 | cut -c1-150
pause


banner "Case 2: config_extra.fdf -- everything mandatory goes there, your template is not edited"
cat <<'EOF'
Each folder gets a config_extra.fdf, and calc.fdf starts with '%include config_extra.fdf' followed by
your template byte for byte. SIESTA's fdf reader is first-occurrence-wins, so the forced directives
beat anything the template says: the cell is the scan variable and must never relax, and every folder
is a single point. The scaled structure keeps your LatticeConstant and the fractional positions;
only the LatticeVectors block changes.
EOF
F="$OUT/case1-si/eos_runs/vol_8.00"
ls "$OUT/case1-si/eos_runs" | tr '\n' ' '; echo; echo "--- vol_8.00:"; ls "$F"
echo
echo "\$ cat eos_runs/vol_8.00/config_extra.fdf"
cat "$F/config_extra.fdf"
echo
echo "\$ head -3 eos_runs/vol_8.00/calc.fdf; diff <(tail -n +3 eos_runs/vol_8.00/calc.fdf) calc.fdf"
head -3 "$F/calc.fdf"
diff <(tail -n +3 "$F/calc.fdf") calc.fdf && echo "-> the template is copied byte for byte below the %include"
echo
echo "\$ grep -A4 'block LatticeVectors' structure.fdf eos_runs/vol_8.00/structure.fdf   # factor (1.08)^(1/3) = 1.0260"
grep -A4 "block LatticeVectors" structure.fdf | head -5
echo "--- scaled:"
grep -A4 "block LatticeVectors" "$F/structure.fdf" | head -5
echo
echo "\$ python3 -c ...  # eos_scan.json, the plan Stage 2 compares the runs with"
python3 - "$OUT/case1-si/eos_runs/eos_scan.json" <<'EOF2'
import json, sys
m = json.load(open(sys.argv[1]))
pts = m["scan"]["points"]
print(f"  {m['structure']['formula']}, V = {m['structure']['volume']:.4f} Ang^3, {m['structure']['space_group']}, "
      f"LatticeConstant {m['structure']['latticeconstant']} Ang")
print(f"  {len(pts)} planned volumes: {pts[0]['folder']} ({pts[0]['volume']:.3f}) ... {pts[-1]['folder']} ({pts[-1]['volume']:.3f}) Ang^3")
EOF2
pause


banner "Case 3: the knobs -- k-grid, ionic relaxation, a different scan"
cat <<'EOF'
The template's own k-grid is kept; --kgrid / --kgrid-density override it, and a template with NO grid
gets one suggested (SIESTA would otherwise sample Gamma only) -- all through config_extra.fdf.
--relax-positions lets the ions move at every volume with the cell fixed (needed when the crystal has
free internal coordinates, e.g. wurtzite); the default is a single point.
EOF
mkdir -p "$OUT/case3-knobs"
cp structure.fdf calc.fdf "$OUT/case3-knobs/"
grep -v "^kgrid" calc.fdf > "$OUT/case3-knobs/calc_nogrid.fdf"
cd "$OUT/case3-knobs"
echo "--- the template has a grid (10x10x10): nothing is forced"
stb-eosInputs -f structure.fdf --n-volumes 5 -o keep --no-intro 2>&1 | strip_colors | grep -E "SCF k-grid forced" | cut -c1-150
echo "--- --kgrid-density 0.1 (finer than the template's)"
stb-eosInputs -f structure.fdf --n-volumes 5 -o dens --kgrid-density 0.1 --no-intro 2>&1 | strip_colors | grep -E "SCF k-grid forced" | cut -c1-170
echo "--- a template without a grid"
stb-eosInputs -f structure.fdf -c calc_nogrid.fdf --n-volumes 5 -o nogrid --no-intro 2>&1 | strip_colors | grep -E "SCF k-grid forced" | cut -c1-170
echo
echo "--- --relax-positions --relax-steps 40: config_extra.fdf of a folder"
stb-eosInputs -f structure.fdf --n-volumes 5 -o relax --relax-positions --relax-steps 40 --no-intro > /dev/null 2>&1
grep -v "^#" relax/vol_0.00/config_extra.fdf | grep -v "^$"
cd "$DIR"
pause


banner "Case 4: the advisories -- what Stage 1 warns about BEFORE you spend CPU time"
cat <<'EOF'
A template that never %includes the structure (SIESTA would NOT see the scaled cell: the whole scan
would silently be the same calculation 11 times), one that defines a lattice of its own, a scan so
narrow that the curvature drowns in SCF noise, and an even number of volumes (no point at V0 of the
input). Each is printed where it applies and repeated at the end.
EOF
mkdir -p "$OUT/case4-advisories"
cp structure.fdf calc.fdf "$OUT/case4-advisories/"
cd "$OUT/case4-advisories"
grep -v "^%include" calc.fdf > calc_noinclude.fdf
sed 's/^%include structure.fdf/%include structure.fdf\nLatticeConstant 5.0 Ang/' calc.fdf > calc_ownlattice.fdf
echo "\$ stb-eosInputs -f structure.fdf -c calc_noinclude.fdf -p dojo"
stb-eosInputs -f structure.fdf -c calc_noinclude.fdf -p dojo --n-volumes 5 -o w1 --no-intro 2>&1 | strip_colors | grep -E "^  - " | cut -c1-210
echo "\$ stb-eosInputs -f structure.fdf -c calc_ownlattice.fdf -p dojo"
stb-eosInputs -f structure.fdf -c calc_ownlattice.fdf -p dojo --n-volumes 5 -o w2 --no-intro 2>&1 | strip_colors | grep -E "^  - " | cut -c1-210
echo "\$ stb-eosInputs -f structure.fdf --n-volumes 6 --strain-range 1"
stb-eosInputs -f structure.fdf --n-volumes 6 --strain-range 1 -o w3 --no-intro 2>&1 | strip_colors | grep -E "even number|very narrow" | cut -c1-210
cd "$DIR"
pause


banner "Case 5: crystal system and dimensionality"
cat <<'EOF'
A hexagonal crystal (hcp Mg) is scaled isotropically, so the scan keeps c/a: the fit then gives the
bulk modulus along that constrained path, which is not the fully relaxed solid's -- the report says
so. A vacuum-padded cell (a graphene sheet) has no bounded volume at all: the tool refuses it.
EOF
mkdir -p "$OUT/case5-symmetry"
cat > "$OUT/case5-symmetry/mg_hcp.fdf" <<'EOF2'
NumberOfSpecies 1
NumberOfAtoms 2
%block ChemicalSpeciesLabel
 1 12 Mg
%endblock ChemicalSpeciesLabel
LatticeConstant 1.0 Ang
AtomicCoordinatesFormat Fractional
%block LatticeVectors
 3.21 0.0 0.0
 -1.605 2.77994 0.0
 0.0 0.0 5.21
%endblock LatticeVectors
%block AtomicCoordinatesAndAtomicSpecies
 0.333333333 0.666666667 0.25 1
 0.666666667 0.333333333 0.75 1
%endblock AtomicCoordinatesAndAtomicSpecies
EOF2
cat > "$OUT/case5-symmetry/graphene.fdf" <<'EOF2'
NumberOfSpecies 1
NumberOfAtoms 2
%block ChemicalSpeciesLabel
 1 6 C
%endblock ChemicalSpeciesLabel
LatticeConstant 1.0 Ang
AtomicCoordinatesFormat Fractional
%block LatticeVectors
 2.46 0.0 0.0
 -1.23 2.130422 0.0
 0.0 0.0 20.0
%endblock LatticeVectors
%block AtomicCoordinatesAndAtomicSpecies
 0.333333333 0.666666667 0.5 1
 0.666666667 0.333333333 0.5 1
%endblock AtomicCoordinatesAndAtomicSpecies
EOF2
cd "$OUT/case5-symmetry"
echo "\$ stb-eosInputs -f mg_hcp.fdf --n-volumes 5"
stb-eosInputs -f mg_hcp.fdf --n-volumes 5 -o hcp --no-intro 2>&1 | strip_colors | grep -E "Space group|^  - " | cut -c1-210
echo "\$ stb-eosInputs -f graphene.fdf          # exit code follows"
set +e
stb-eosInputs -f graphene.fdf -o sheet --no-intro 2>&1 | strip_colors | grep -E "ERROR" | cut -c1-210
echo "exit code: ${PIPESTATUS[0]}"
set -e
cd "$DIR"
pause


banner "Case 6: Stage 2 -- the fit of the real silicon scan, every EOS form"
cat <<'EOF'
reference_results/si/base/ is the real scan of Case 1's commands: 11 finished SIESTA runs (DZP,
Mesh.Cutoff 320 Ry, 10x10x10 k-points). --eos all fits Birch-Murnaghan, Vinet, Murnaghan and
stabilized jellium on the same data. Read [1] (is every folder at the same level of theory? did SIESTA
run the planned cells?), [3] (the data, with the pressure from each run's own stress tensor) and [4]
(the fit, its uncertainty and its sensitivity to the fitted window).
EOF
mkdir -p "$OUT/case6-fit"
echo "\$ stb-eosAnalysis --dir reference_results/si/base --eos all --save-report --no-intro"
(cd "$OUT/case6-fit" && stb-eosAnalysis --dir "$REF/base" --eos all --save-report --no-intro > stage2.log 2>&1)
section "$OUT/case6-fit/stage2.log" 1 | cut -c1-150
section "$OUT/case6-fit/stage2.log" 3 | cut -c1-175
section "$OUT/case6-fit/stage2.log" 4 | cut -c1-175
pause


banner "Case 7: the pressure cross-check and the checks table"
cat <<'EOF'
SIESTA prints a stress tensor in every run. Its trace gives a pressure P_SIESTA that must equal -dE/dV
of the fitted curve -- two independent routes to the same physics (the atom-centred orbitals keep their
shape and range when the cell is scaled, so the derivative is exact up to the real-space mesh). The
stress alone also gives a second B0 = -V dP/dV. [7] then collects every test the report makes on the
fit -- bracketing, residuals, B0', agreement between forms, jackknife and window stability, pressure
agreement, residual forces, same level of theory -- as OK/WARN.
EOF
section "$OUT/case6-fit/stage2.log" 5 | cut -c1-175
section "$OUT/case6-fit/stage2.log" 7 | cut -c1-175
pause


banner "Case 8: plots and data on request -- dashboard, gnuplot, CSV/JSON, target pressure"
cat <<'EOF'
--show opens a matplotlib dashboard (6 panels: E(V) with every fit, residuals, P(V) from the stress and
from the fit, B along the curve, B0 by EOS form and by method with the leave-one-out band, and a
summary); --png saves it; --save-gnuplot writes the .dat/.gplot pair; --save-data a CSV table and a
JSON with every number. Nothing is written unless asked for. This script runs headless, so --show is
replaced by --png (try --show on a desktop). --target-pressure inverts the fit -- it needs one specific
EOS form, so it is a separate command: the volume (and the linear scale of the cell) at 0, 5 and 10 GPa.
Notice the last row: the scan only reached +8.65 GPa (at -8 % volume), so 10 GPa is an extrapolation,
and the report says so.
EOF
mkdir -p "$OUT/case8-outputs"
echo "\$ stb-eosAnalysis --dir reference_results/si/base --eos all --png --save-gnuplot --save-data -o eos --no-intro"
(cd "$OUT/case8-outputs" && stb-eosAnalysis --dir "$REF/base" --eos all --png --save-gnuplot --save-data \
    -o eos --no-intro > stage2.log 2>&1)
ls -1 "$OUT/case8-outputs" | grep -v "\.log$" | sed 's/^/  /'
cp "$OUT/case8-outputs/eos_summary.png" "$DIR/figures/si_dashboard.png"
echo "  (figures/si_dashboard.png is this dashboard)"
python3 - "$OUT/case8-outputs/eos_results.json" <<'EOF2'
import json, sys
r = json.load(open(sys.argv[1]))
f = r["fits"]["birchmurnaghan"]
print(f"  eos_results.json: V0 {f['v0']:.4f} Ang^3, B0 {f['b0_gpa']:.2f} +/- {r['jackknife']['b0_sigma']:.2f} GPa, "
      f"B0' {f['bprime']:.3f}; B0 from the stress {r['b0_from_stress_gpa']:.2f} GPa; "
      f"{sum(c['ok'] for c in r['checks'])} of {len(r['checks'])} checks OK; suggested LatticeConstant {r['suggested_latticeconstant']:.5f} Ang")
EOF2
echo
echo "\$ stb-eosAnalysis --dir reference_results/si/base --target-pressure 0 5 10 --no-intro"
(cd "$OUT/case8-outputs" && stb-eosAnalysis --dir "$REF/base" --target-pressure 0 5 10 --no-intro > target.log 2>&1)
section "$OUT/case8-outputs/target.log" 6 | cut -c1-175
pause


banner "Case 9: convergence -- which knob moves which number"
cat <<'EOF'
Thirteen more real scans, each changing exactly one thing from the base: the k-grid (4 .. 14), the
Mesh.Cutoff (150 .. 600 Ry) and the basis (DZP/TZP with PAO.EnergyShift 0.02, 0.01, 0.005 Ry).
convergence_table.py fits each with the same code as Stage 2. E0 shift is the total-energy change
relative to the densest setting; max dP is the worst disagreement between the stress and -dE/dV.
EOF
echo "\$ python3 convergence_table.py --png"
python3 convergence_table.py --png 2>&1 | grep -v Warn
pause


banner "Case 10: compared with an all-electron calculation -- why B0 alone misleads"
cat <<'EOF'
The tables add the all-electron PBE reference of the Delta test (WIEN2k: V0 = 20.453 Ang^3/atom,
B0 = 88.545 GPa, B0' = 4.31) and Delta, the root-mean-square energy difference between the two
equation-of-state curves (meV/atom, curves aligned at their minima). 'B(V_AE)' is the bulk modulus of
OUR curve at the all-electron volume. Read the DZP rows: lengthening the orbitals (smaller
PAO.EnergyShift) brings B0 from 94.3 to the reference 88.5 GPa -- but V0 grows from 20.33 to 20.70
(the reference is 20.45), at the reference volume the modulus is 91.9 vs 93.0 GPa, and Delta DOUBLES.
And the smallest Delta of all the k-point runs (0.60) belongs to the 4x4x4 grid, whose total energy
is 93 meV/atom off: Delta and B0 can be right for the wrong reason. Read V0, B at a fixed volume,
Delta AND the convergence of the total energy together.
EOF
python3 convergence_table.py 2>&1 | grep -v Warn | grep -E "^(run|k04|k14|dzp|tzp|All)|^base|^ " | awk '!seen[$0]++'
pause


banner "Case 11: guard rails -- what Stage 2 does with a damaged scan"
cat <<'EOF'
Three damaged copies of the base scan: one run missing (skipped, not fatal), one folder edited to a
different mesh cutoff (the level of theory differs: E(V) is no longer one curve), and a scan that does
not bracket the minimum (only the compressed side: V0 is outside it). Each is reported where it
applies, in the checks table, and again in the warnings.
EOF
mkdir -p "$OUT/case11-damaged"
cd "$OUT/case11-damaged"
cp -r "$REF/base" a_missing; rm a_missing/vol_3.20/calc.out
cp -r "$REF/base" b_theory; sed -i 's/^Mesh.Cutoff.*/Mesh.Cutoff 200 Ry/' b_theory/vol_1.60/calc.fdf
mkdir c_edge; for v in vol_m8.00 vol_m6.40 vol_m4.80 vol_m3.20 vol_m1.60; do cp -r "$REF/base/$v" c_edge/; done
cp "$REF/base/eos_scan.json" c_edge/
for d in a_missing b_theory c_edge; do
    echo "------------------------------------------------------------------"
    echo "\$ stb-eosAnalysis --dir $d --no-intro      # warnings only"
    stb-eosAnalysis --dir $d --no-intro 2>&1 | strip_colors | sed -n '/^\[10\] WARNINGS/,/^\[INFO\]/p' | grep -E "^ *[0-9]+\. " | cut -c1-210
done
cd "$DIR"

echo
echo "Done. Every case wrote into $OUT/; the figures of the README are in figures/."
