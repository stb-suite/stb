#!/bin/bash
# Guided example: Workflow 5.13 -- GQCA, Generalized Quasi-Chemical
# Approximation (both stages: stb-gqca, code 5.13.1, and stb-gqcaAnalysis,
# code 5.13.2, in the stb-suite menu).
#
# Not an automated test (see test/5-workflow/13-gqca/{prep,analysis}/test.sh
# for that) -- a commented walk-through: it runs real commands, one case at a
# time, into its own output/<case>/ folder, and shows you the piece of output
# that proves what just happened. Pauses between sections so you can read
# before moving on. Safe to re-run any time -- it always starts by wiping its
# own output/.
#
# structure.fdf/calc.fdf are GaAs (zincblende, DZP/GGA-PBE) -- the host of the
# two alloys studied here: Ga(1-x)In(x)As (size-mismatched, ~7%) and
# Ga(1-x)Al(x)As (nearly lattice-matched, ~0.1%).
#
# Stage 1 (stb-gqca) runs for real in every case below (cheap -- it only
# writes .fdf files, no SCF). Stage 2 (stb-gqcaAnalysis) runs on REAL SIESTA
# results: reference_results/{ingaas,algaas}/ hold the 3 finished cluster
# relaxations of each alloy (calc.out, siesta.XV, the exact input files),
# computed with the stb-gqca output of Case 1 and the calc.fdf in this folder
# (SIESTA 5.4.2, mpirun -np 2, ~15 min per alloy). Nothing in Stage 2 is
# fabricated -- every number the tool prints is genuine output on genuine DFT
# data. See the README (section 6) for how to regenerate them yourself.

set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

# Stage 2's --show / --view need a display; MPLBACKEND=Agg keeps this script
# headless-safe (the figure is still saved with --png).
export MPLBACKEND=Agg

OUT="$DIR/output"
rm -rf "$OUT"
mkdir -p "$OUT"

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


banner "Welcome: from 3 DFT relaxations to a whole alloy phase behaviour"
cat <<'EOF'
A substitutional alloy A(1-x)Bx has astronomically many atomic arrangements,
and you cannot run DFT on all of them. The pair-cluster (N=2) Generalized
Quasi-Chemical Approximation (Sher, van Schilfgaarde, Chen, Chen, Phys. Rev.
B 36, 4279, 1987) sidesteps that: it treats the alloy as a gas of independent
nearest-neighbour PAIRS, each AA, AB or BB, in a mass-action equilibrium
    2 AB  <->  AA + BB
driven by only THREE numbers -- the energies of the three ordered reference
structures. Stage 1 (stb-gqca) builds and prepares them; you run SIESTA on
them ONCE; Stage 2 (stb-gqcaAnalysis) turns the three energies into the
mixing enthalpy, entropy, free energy and short-range order at ANY
composition x and ANY temperature T.

This tour: Stage 1 on GaAs -> InGaAs (Cases 1-6), then Stage 2 on REAL SIESTA
results for InGaAs and AlGaAs (Cases 7-11). See README.md for the theory
(section 1) and the worked, literature-checked example (section 5).
EOF
pause


banner "Case 1: Stage 1 -- GaAs -> InGaAs, what stb-gqca reports"
cat <<'EOF'
--sublattice Ga marks the cation sites as THE alloy sublattice, --species-b In
is what replaces Ga on it. No -c flag: calc.fdf in the current folder is the
default template. -p dojo copies the PseudoDojo pseudopotentials into each
folder. Watch [1] (what was read), [2] (the nearest-neighbour pair the whole
method is built on) and [4] (the 3 structures written).
EOF
mkdir -p "$OUT/case1-stage1"
cp structure.fdf calc.fdf "$OUT/case1-stage1/"
echo "\$ stb-gqca -f structure.fdf --sublattice Ga --species-b In -p dojo --no-intro"
(cd "$OUT/case1-stage1" && stb-gqca -f structure.fdf --sublattice Ga --species-b In \
    -p dojo --no-intro > stage1.log 2>&1)
strip_colors < "$OUT/case1-stage1/stage1.log" | sed -n '/\[1\] INPUT STRUCTURE/,/\[3\] CONFIG_EXTRA/p' | sed '$d'
echo
strip_colors < "$OUT/case1-stage1/stage1.log" | sed -n '/Written structures/,/Folder contents/p' | sed '$d'
echo
echo "On disk:"
(cd "$OUT/case1-stage1" && ls gqca_study gqca_study/cluster_n1)
pause


banner "Case 2: config_extra.fdf -- how the relaxation is forced WITHOUT touching your calc.fdf"
cat <<'EOF'
The 3 clusters must be REALLY relaxed (positions AND cell: In and Ga have
different bond lengths, and a fixed cell would put spurious strain into the
mixing energy -- README section 1.5). stb-gqca forces that through a
config_extra.fdf written next to each calc.fdf, %included at the very TOP
of it. SIESTA's fdf reader is first-occurrence-wins, so the forced values win
even over directives already in your template -- and your template is copied
byte for byte.
EOF
F="$OUT/case1-stage1/gqca_study/cluster_n1"
echo "\$ cat gqca_study/cluster_n1/config_extra.fdf"
cat "$F/config_extra.fdf"
echo
echo "\$ head -4 gqca_study/cluster_n1/calc.fdf     # the %include is line 1"
head -4 "$F/calc.fdf"
echo
echo "The rest of calc.fdf is your template, unchanged:"
if diff <(tail -n +3 "$F/calc.fdf") calc.fdf > /dev/null; then
    echo "  confirmed -- calc.fdf (minus the 2 include lines) is identical to ./calc.fdf"
fi
echo
echo "--- A hostile template: kgrid 1x1x1, MD.NumCGsteps 0 (= no relaxation!), Broyden ---"
mkdir -p "$OUT/case2-hostile"
cp structure.fdf "$OUT/case2-hostile/"
{ cat calc.fdf; printf 'kgrid.MonkhorstPack [1 1 1]\nMD.TypeOfRun Broyden\n'; } > "$OUT/case2-hostile/calc_hostile.fdf"
sed -i 's/^MD.NumCGsteps 150/MD.NumCGsteps 0/' "$OUT/case2-hostile/calc_hostile.fdf"
echo "\$ stb-gqca -f structure.fdf --sublattice Ga --species-b In -c calc_hostile.fdf --no-intro"
(cd "$OUT/case2-hostile" && stb-gqca -f structure.fdf --sublattice Ga --species-b In \
    -c calc_hostile.fdf --no-intro 2>&1 | strip_colors | sed -n '/\[3\] CONFIG_EXTRA/,/^Why:/p')
echo
echo "Resulting config_extra.fdf (steps raised, k-grid + run type overridden):"
cat "$OUT/case2-hostile/gqca_study/cluster_n0/config_extra.fdf"
pause


banner "Case 3: the k-point density -- --kdensity (default 0.2 1/Ang)"
cat <<'EOF'
The k-grid is recomputed PER FOLDER from a density (same convention as
stb-kgrid and the other workflows), because the 3 cells differ in shape and
size -- the same density gives each its own proper grid. 0.2 is the suite's
usual starting point for semiconductors; tighten it (0.1, 0.05) for metals or
when you want a convergence check of the mixing energy itself (README 5.4).
EOF
for d in 0.5 0.2 0.1; do
    mkdir -p "$OUT/case3-k$d"
    cp structure.fdf calc.fdf "$OUT/case3-k$d/"
    echo "\$ stb-gqca ... --kdensity $d"
    (cd "$OUT/case3-k$d" && stb-gqca -f structure.fdf --sublattice Ga --species-b In \
        --kdensity "$d" --no-intro 2>&1 | strip_colors | grep -E "^(cluster_n[012]|Folder) +\\|" | cut -d'|' -f1,10)
    echo
done
pause


banner "Case 4: three lattices, three behaviours -- bipartite vs frustrated vs 2D"
cat <<'EOF'
The AB structure should make every nearest-neighbour pair UNLIKE. Whether
that is geometrically possible depends on the lattice:
  * simple cubic (cubic.fdf, Na -> K)   -- bipartite, a checkerboard works
  * FCC cation sublattice (GaAs)        -- triangles, can't be perfect
  * triangular Mo sublattice of a 2D MoS2 monolayer (mos2.fdf, Mo -> W)
    -- frustrated AND vacuum-padded, so positions-only relaxation
Compare the "Pair correlation" line (ideal AB = -1) and the cell treatment.
EOF
mkdir -p "$OUT/case4-lattices"
cp cubic.fdf mos2.fdf calc.fdf "$OUT/case4-lattices/"
for spec in "cubic.fdf Na K" "structure.fdf Ga In" "mos2.fdf Mo W"; do
    set -- $spec
    [ "$1" = "structure.fdf" ] && cp structure.fdf "$OUT/case4-lattices/"
    echo "------------------------------------------------------------------"
    echo "\$ stb-gqca -f $1 --sublattice $2 --species-b $3 -O study_$2$3"
    (cd "$OUT/case4-lattices" && stb-gqca -f "$1" --sublattice "$2" --species-b "$3" \
        -O "study_$2$3" --no-intro 2>&1 | strip_colors > "log_$2$3.txt")
    grep -E "^Detected|^Cell treatment|Pair correlation|Deviation from ideal" "$OUT/case4-lattices/log_$2$3.txt"
done
echo
echo "A 2D input also prints a [KNOWN LIMITATION] note (cell relaxation OFF, positions only):"
grep -c "KNOWN LIMITATION" "$OUT/case4-lattices/log_MoW.txt" | sed 's/^/  occurrences in the MoS2 report: /'
pause


banner "Case 5: validation -- what stb-gqca correctly refuses"
cat <<'EOF'
Four real mistakes, each exiting with a clear [ERROR] and a non-zero status:
an unknown sublattice species, sublattice == species B, a missing structure,
and a missing calc.fdf (the default template that is not there).
EOF
mkdir -p "$OUT/case5-validation"
cp structure.fdf "$OUT/case5-validation/"     # note: NO calc.fdf copied
cd "$OUT/case5-validation"
echo "\$ stb-gqca -f structure.fdf --sublattice Xx --species-b In -c $DIR/calc.fdf"
stb-gqca -f structure.fdf --sublattice Xx --species-b In -c "$DIR/calc.fdf" --no-intro 2>&1 | strip_colors | grep ERROR || true
echo
echo "\$ stb-gqca -f structure.fdf --sublattice Ga --species-b Ga -c $DIR/calc.fdf"
stb-gqca -f structure.fdf --sublattice Ga --species-b Ga -c "$DIR/calc.fdf" --no-intro 2>&1 | strip_colors | grep ERROR || true
echo
echo "\$ stb-gqca -f nothere.fdf --sublattice Ga --species-b In -c $DIR/calc.fdf"
stb-gqca -f nothere.fdf --sublattice Ga --species-b In -c "$DIR/calc.fdf" --no-intro 2>&1 | strip_colors | grep ERROR || true
echo
echo "\$ stb-gqca -f structure.fdf --sublattice Ga --species-b In      # no ./calc.fdf here"
stb-gqca -f structure.fdf --sublattice Ga --species-b In --no-intro 2>&1 | strip_colors | grep ERROR || true
cd "$DIR"
pause


banner "Case 6: CLI vs. the interactive stb-suite menu (5.13.1) -- same result"
cat <<'EOF'
The menu asks the same questions instead of flags. Reproducing Case 1 through
stb-suite -> 5.13.1 (structure, sublattice Ga, species B In, calc.fdf default,
pseudopotential bank 1 = dojo, output name, no advanced settings) and
comparing the generated folders with Case 1's proves the two paths equal.
EOF
mkdir -p "$OUT/case6-interactive"
cp structure.fdf calc.fdf "$OUT/case6-interactive/"
(cd "$OUT/case6-interactive" && \
    printf '5.13.1\nstructure.fdf\nGa\nIn\n\n1\n\n\n\n0\n' | stb-suite > menu.log 2>&1)
if diff -rq -x '*.psml' "$OUT/case1-stage1/gqca_study" "$OUT/case6-interactive/gqca_study" 2>&1 \
        | grep -v "gqca_stage1.txt" | grep -q .; then
    echo "UNEXPECTED: CLI and interactive menu folders differ -- see"
    echo "  diff -rq -x '*.psml' '$OUT/case1-stage1/gqca_study' '$OUT/case6-interactive/gqca_study'"
else
    echo "CLI (Case 1) and interactive menu (5.13.1) produced identical cluster folders"
    echo "(the Stage 1 report differs only by its own timestamp)."
fi
pause


banner "Case 7: Stage 2 on REAL SIESTA data -- Ga(1-x)In(x)As"
cat <<'EOF'
reference_results/ingaas/ holds the 3 finished SIESTA relaxations of this
alloy. Stage 2 reads each folder's calc.out, normalizes the energies per
sublattice site (the 3 cells have different sizes!), and solves the pair
mass-action equilibrium. Default: --directory gqca_study (here we point it at
the reference), calc.out files, and a 300-1500 K sweep. Four sections of the
report, in turn: level of theory, convergence/geometry, energies, tables.
EOF
cp -r reference_results/ingaas "$OUT/ingaas"
echo "\$ stb-gqcaAnalysis --directory output/ingaas --x-target 0.25 0.5 --png --no-intro"
stb-gqcaAnalysis --directory "$OUT/ingaas" --x-target 0.25 0.5 --png \
    --no-intro 2>&1 | strip_colors > "$OUT/ingaas_stage2.log"
echo
echo "----- [1] Did the 3 folders use the SAME level of theory? -----"
sed -n '/\[1\] LEVEL OF THEORY/,/^All 3 folders/p' "$OUT/ingaas_stage2.log" | tail -n +2
pause
echo "----- [2] Are they converged, and what did they relax to? -----"
sed -n '/^Convergence:/,/AB volume vs Vegard/p' "$OUT/ingaas_stage2.log"
pause
echo "----- [3] The three energies that carry the whole method -----"
sed -n '/\[3\] MIXING ENERGIES/,/H_mix(x=0.5/p' "$OUT/ingaas_stage2.log" | tail -n +2
sed -n '/T_c = Omega/p' "$OUT/ingaas_stage2.log"
pause
echo "----- [4] What happens at x = 0.5 as the temperature rises -----"
sed -n '/x = 0.500/,/^$/p' "$OUT/ingaas_stage2.log"
echo "(random-alloy x_1 would be 0.50; SRO > 0 = like pairs clustered, -> 0 as T grows)"
pause
echo "----- [5]-[6] Sanity checks and the plain-language interpretation -----"
sed -n '/\[5\] SANITY/,/\[7\] MISCIBILITY/p' "$OUT/ingaas_stage2.log" | sed '$d' | tail -n +2
pause


banner "Case 8: the contrast -- REAL Ga(1-x)Al(x)As (and why you must not trust one basis)"
cat <<'EOF'
Same workflow, same level of theory, on the NEARLY LATTICE-MATCHED alloy: AlAs
and GaAs differ by ~0.1% experimentally, so there is almost no strain to pay
for and the mixing enthalpy should be tiny. Put the two alloys side by side
(both with the DZP basis of calc.fdf). Then read the surprise: Delta_1 of
AlGaAs is NOT near zero. Case 9 investigates.
EOF
cp -r reference_results/algaas "$OUT/algaas"
stb-gqcaAnalysis --directory "$OUT/algaas" --x-target 0.5 --no-intro 2>&1 | strip_colors > "$OUT/algaas_stage2.log"
echo "\$ stb-gqcaAnalysis --directory output/algaas --x-target 0.5 --no-intro      # [3] only"
sed -n '/\[3\] MIXING ENERGIES/,/T_c = Omega/p' "$OUT/algaas_stage2.log" | tail -n +2 | cut -c1-110
echo
python3 "$DIR/compare_runs.py" "$OUT" ingaas algaas
pause


banner "Case 9: convergence of the mixing energy -- basis (DZP/TZP) and mesh cut-off"
cat <<'EOF'
Delta_1 is a small difference (tens of meV) of three total energies of
thousands of eV: it inherits every incompleteness of the numerical set-up.
reference_results/ also holds the same 3 relaxations of each alloy re-run with
PAO.BasisSize TZP, and of AlGaAs with Mesh.Cutoff 450 Ry. The tool reads the
basis and mesh back from each folder (watch the first two rows). Watch also
the GaAs lattice constant (experiment 5.653 A; plane-wave PBE ~5.67 A).
EOF
for r in ingaas_tzp algaas_tzp algaas_mesh450; do
    cp -r reference_results/$r "$OUT/$r"
    stb-gqcaAnalysis --directory "$OUT/$r" --x-target 0.5 --no-intro 2>&1 | strip_colors > "$OUT/${r}_stage2.log"
done
python3 "$DIR/compare_runs.py" "$OUT" ingaas ingaas_tzp algaas algaas_tzp algaas_mesh450
echo
echo "Reading it: TZP moves Delta_1 by about -7 meV (InGaAs) and -11 meV (AlGaAs)"
echo "and brings the GaAs lattice constant down 0.05 A; the denser mesh does NOT"
echo "lower Delta_1 (+3 meV). Delta_1 of AlGaAs stays far from the ~few-meV value"
echo "an almost ideal alloy should have: DO NOT quote a number from one set-up."
echo "README section 5.3 lists what was ruled out and what remains open."
pause


banner "Case 10: the optional extras -- dashboard, CSV, gnuplot, relaxed structures"
cat <<'EOF'
Besides the report, Stage 2 writes machine-readable and visual outputs, and
can show them on screen:
  --png / --show        6-panel dashboard (saved / opened with matplotlib)
  --view                the 3 relaxed structures in ASE's 3D viewer
  --export-relaxed FMT  the 3 relaxed structures (cif|xsf|xyz|vasp) to files
  --save-report         the text report as gqca_stage2.txt
Here the headless-safe ones (--show/--view need a display, see the README):
EOF
echo "\$ stb-gqcaAnalysis --directory output/ingaas --png --export-relaxed cif --save-report --no-intro"
stb-gqcaAnalysis --directory "$OUT/ingaas" --png --export-relaxed cif --save-report \
    --no-intro 2>&1 | strip_colors | sed -n '/\[8\] SUMMARY/,/UNVERIFIED/p' | sed '$d' | sed "s#$DIR/##g"
echo
echo "gqca_results_clusters.csv (one row per cluster):"
cat "$OUT/ingaas/gqca_results_clusters.csv"
echo
echo "gqca_results.csv, first rows (x, T, eta, x_0, x_1, x_2, H, S, G, SRO):"
head -3 "$OUT/ingaas/gqca_results.csv"
echo
echo "Open the dashboard yourself:  xdg-open output/ingaas/gqca_results_summary.png"
echo "On a desktop you can also run:  stb-gqcaAnalysis --directory output/ingaas --show --view"
pause


banner "Case 11: Stage 2 guard rails -- what it flags when the data is not trustworthy"
cat <<'EOF'
Three realistic problems, all detected from the files alone:
  (a) results missing         -> clean error, nothing computed
  (b) one folder computed with a different level of theory (e.g. a 200 Ry
      mesh in cluster_n1) -> WARNING in [1]: mixing energies are only
      meaningful when all 3 clusters are computed identically
  (c) a calc.out that never confirmed SCF convergence -> WARNING in [2]
EOF
echo "--- (a) an empty folder"
mkdir -p "$OUT/guard_empty"
stb-gqcaAnalysis --directory "$OUT/guard_empty" --temp 500 --no-intro 2>&1 | strip_colors | grep ERROR | sed "s#$DIR/##" || true
echo
echo "--- (b) different mesh cutoff in cluster_n1"
cp -r reference_results/ingaas "$OUT/guard_theory"
sed -i 's/^Mesh.Cutoff 250 Ry/Mesh.Cutoff 200 Ry/' "$OUT/guard_theory/cluster_n1/calc.fdf"
stb-gqcaAnalysis --directory "$OUT/guard_theory" --temp 500 --no-intro 2>&1 | strip_colors \
    | grep -E "^Mesh cutoff|WARNING.*level of theory" | cut -c1-175
echo
echo "--- (c) SCF convergence line removed from cluster_n2's calc.out"
cp -r reference_results/ingaas "$OUT/guard_scf"
sed -i '/SCF cycle converged after/d' "$OUT/guard_scf/cluster_n2/calc.out"
stb-gqcaAnalysis --directory "$OUT/guard_scf" --temp 500 --no-intro 2>&1 | strip_colors \
    | grep -E "WARNING.*SCF|all clusters SCF" | sed "s#$DIR/##" | cut -c1-175
pause


banner "Case 12: CLI vs. the interactive stb-suite menu (5.13.2) -- same result"
cat <<'EOF'
Stage 2 through stb-suite -> 5.13.2: directory, default output name, option 1
(single temperature) = 500 K, default target composition, none of the
optional extras except saving the report, no advanced settings.
EOF
cp -r reference_results/ingaas "$OUT/menu_study"
(cd "$OUT" && printf '5.13.2\nmenu_study\n\n1\n500\n\n\n\n\n\ny\n\n\n0\n' | stb-suite > menu_stage2.log 2>&1)
stb-gqcaAnalysis --directory "$OUT/ingaas" --temp 500 -o cli_check --no-intro > /dev/null 2>&1
if diff -q "$OUT/menu_study/gqca_results.csv" "$OUT/ingaas/cli_check.csv" > /dev/null; then
    echo "CLI and interactive menu (5.13.2) produced an identical results CSV."
else
    echo "UNEXPECTED: the two results CSVs differ."
fi
pause


banner "Workflow 5.13 complete"
cat <<EOF
Stage 1 (stb-gqca, 5.13.1): 3 relaxation-ready folders (cluster_n0/n1/n2):
structure.fdf, your untouched calc.fdf, a config_extra.fdf that forces the
relaxation + the per-cell k-grid, and the pseudopotentials.
Stage 2 (stb-gqcaAnalysis, 5.13.2): level-of-theory check, convergence and
relaxed geometry, Delta_j in several units, regular-solution Omega and T_c,
tables at your target compositions, sanity checks, an interpretation, CSV /
gnuplot / PNG dashboard / relaxed structures.

On YOUR alloy (any lattice, 2D or 3D):
  stb-gqca -f structure.fdf --sublattice <A> --species-b <B> -p dojo
  for c in n0 n2 n1; do (cd gqca_study/cluster_\$c && mpirun -np 4 siesta calc.fdf > calc.out); done
  stb-gqcaAnalysis --x-target 0.25 0.5 0.75 --png
Read [WARNING]/[NOTE] lines in both reports first -- README section 5.5
explains how far to trust the numbers.
EOF
