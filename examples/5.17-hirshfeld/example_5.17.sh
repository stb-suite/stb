#!/bin/bash
# Guided example: Workflow 5.17 -- Hirshfeld-I charges (all three stages: stb-hirshfeldPrep, code 5.17.1,
# stb-hirshfeldIons, code 5.17.2, and stb-hirshfeldAnalysis, code 5.17.3, in the stb-suite menu).
#
# Not an automated test (see test/5-workflow/17-hirshfeld/{prep,ions,analysis}/test.sh for that) -- a commented
# walk-through: it runs real commands, one case at a time, into its own output/<case>/ folder, and shows the
# piece of output that proves what just happened. Pauses between sections so you can read before moving on.
# Safe to re-run any time -- it always starts by wiping its own output/.
#
# Stages 1 and 2 only write folders (cheap, no SCF) and run for real here. Stage 3 runs on REAL SIESTA
# results: reference_results/ holds the finished runs of H2O, N2O, NaCl and MgO (every combined/,
# neutral/<El>/ and ions/<El>/<cation|anion>/ folder) plus the H2O studies (mesh, box, basis), computed with
# the commands of this script and the inputs in this folder (SIESTA 5.4.2; see README section 6 and
# make_reference.sh to regenerate them). The density grids are stored compressed (siesta.RHO.xz) and unpacked
# into output/ when a case needs them. Nothing in Stage 3 is fabricated; this script never runs SIESTA.

set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

# --show needs a display; MPLBACKEND=Agg keeps this script headless-safe (the figure is still saved with --png).
export MPLBACKEND=Agg

OUT="$DIR/output"
rm -rf "$OUT"
mkdir -p "$OUT"
REF="$DIR/reference_results"

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

# a finished reference study, copied into output/ with every density grid unpacked
unpack() {   # $1 = study name, $2 = destination
    cp -r "$REF/$1" "$2"
    find "$2" -name "siesta.RHO.xz" -exec xz -d {} \;
    return 0
}


banner "Welcome: how many electrons does an atom in a molecule own?"
cat <<'EOF'
The density of a molecule or a crystal is an observable; the charge of one of its atoms is not -- it needs a
rule that cuts the density into atoms. The Hirshfeld (stockholder) rule shares the density at every point in
proportion to free reference atoms placed on the nuclei. With NEUTRAL references the charges come out too
small; Hirshfeld-I makes every reference atom hold as many electrons as the atom it defines, by iteration,
mixing in the cation or the anion of each element:

    reference(q) = (1 - |q|) neutral + |q| ion,        iterate until q stops changing

Stage 1 (stb-hirshfeldPrep) writes the combined system (single point, SaveRho) and one isolated NEUTRAL atom
per element. Stage 2 (stb-hirshfeldIons) reads the finished neutral atoms and writes the cation and the anion
of every element. Stage 3 (stb-hirshfeldAnalysis) builds the reference profiles, iterates, and reports the
charges next to Voronoi, SIESTA's Hirshfeld and simple Hirshfeld, with checks of every reference, how far
the iteration still is from its fixed point, and a dipole check.

This tour: the three stages on water (Cases 1-4), how converged "converged" is (Case 5), N2O with two
nitrogens of opposite sign (Case 6), the NaCl crystal and its periodic images (Case 7), MgO beyond +-1
(Case 8), what the charges depend on (Case 9) and the guard rails (Case 10). See README.md for the theory
(section 1) and the worked example (section 5).
EOF
pause


banner "Case 1: Stage 1 -- the combined system and the neutral atoms"
cat <<'EOF'
reference_results/relax/h2o/relaxed.fdf is water relaxed with calc.fdf (PBE, DZP, PAO.EnergyShift 0.01 Ry)
in a 14 Ang box. Stage 1 writes combined/ (the molecule itself, forced to a single point with SaveRho and
SIESTA's own Hirshfeld and Voronoi charges) and neutral/O/, neutral/H/ (one atom in a 20 Ang box, Gamma only,
spin-polarized), all with the same MeshCutoff. -p dojo copies the PseudoDojo pseudopotentials.
EOF
mkdir -p "$OUT/case1-water"
cp "$REF/relax/h2o/relaxed.fdf" "$OUT/case1-water/h2o.fdf"
cp calc.fdf "$OUT/case1-water/"
echo "\$ stb-hirshfeldPrep -s h2o.fdf --calc calc.fdf -p dojo --mesh-cutoff 400 --save-report --no-intro"
(cd "$OUT/case1-water" && stb-hirshfeldPrep -s h2o.fdf --calc calc.fdf -p dojo --mesh-cutoff 400 --save-report \
    --no-intro > stage1.log 2>&1)
section "$OUT/case1-water/stage1.log" 2 | cut -c1-150
section "$OUT/case1-water/stage1.log" 3 | cut -c1-200
S="$OUT/case1-water/hirshfeld_study"
echo "\$ ls hirshfeld_study hirshfeld_study/neutral"
ls "$S"; ls "$S/neutral"
echo
echo "\$ grep -v '^#' hirshfeld_study/neutral/O/config_extra.fdf"
grep -v "^#" "$S/neutral/O/config_extra.fdf" | grep -v "^$" | sed 's/^/  /'
echo
echo "\$ cat hirshfeld_study/hirshfeld_manifest.json     # paths relative to the study: it can be copied or moved"
python3 -c "import json,sys; m=json.load(open(sys.argv[1])); print('  ', {k: m[k] for k in ('combined_dir','combined_structure') if k in m})" \
    "$S/hirshfeld_manifest.json"
pause


banner "Case 2: Stage 2 -- after SIESTA ran in combined/ and neutral/*, the ions"
cat <<'EOF'
reference_results/h2o/ is the same study after SIESTA ran in every folder. Here its ion folders are removed
and Stage 2 writes them again: it reads Z_val from the SIESTA log, previews the simple-Hirshfeld charges (pass
0, what Stage 3 will start from) and writes BOTH ions of every element, NetCharge +1 and -1, because each atom
picks its own ion in Stage 3. H has Z_val = 1: its cation is a bare proton with no electrons, exactly zero
density, so no folder is written for it.
EOF
mkdir -p "$OUT/case2-ions"
unpack h2o "$OUT/case2-ions/h2o"
rm -rf "$OUT/case2-ions/h2o/ions" "$OUT/case2-ions/h2o/hirshfeld_ions_manifest.json"
echo "\$ stb-hirshfeldIons -O h2o -p dojo --save-report --no-intro"
(cd "$OUT/case2-ions" && stb-hirshfeldIons -O h2o -p dojo --save-report --no-intro > stage2.log 2>&1)
section "$OUT/case2-ions/stage2.log" 1 | sed -n '/^Species/,$p' | cut -c1-150
section "$OUT/case2-ions/stage2.log" 4 | sed -n '/^Idx/,/^$/p' | cut -c1-150
section "$OUT/case2-ions/stage2.log" 5 | cut -c1-220
echo "\$ grep -v '^#' h2o/ions/O/anion/config_extra.fdf | grep -E 'NetCharge|Spin'"
grep -v "^#" "$OUT/case2-ions/h2o/ions/O/anion/config_extra.fdf" | grep -E "NetCharge|Spin" | sed 's/^/  /'
pause


banner "Case 3: Stage 3 -- the reference atoms and the iteration"
cat <<'EOF'
With SIESTA finished in the ion folders too, Stage 3 reads every reference, averages it over angles and
checks it ([2]: electrons on the grid and in the radial profile against Z_val, Z_val - 1, Z_val + 1, the
radius that holds 99.9 % of it against its box, SCF). Then it iterates ([3], [4]): the step max|Delta q|
shrinks by a constant ratio, about 0.67 for water. Note the warning: the last step is below --tol, the
charges are not -- at that ratio they are still ~2x the last step from the fixed point.
EOF
mkdir -p "$OUT/case3-water"
unpack h2o "$OUT/case3-water/h2o"
echo "\$ stb-hirshfeldAnalysis -O h2o --save-data --png --save-report --no-intro"
(cd "$OUT/case3-water" && stb-hirshfeldAnalysis -O h2o --save-data --png --save-report --no-intro > stage3.log 2>&1)
section "$OUT/case3-water/stage3.log" 2 | sed -n '/^Elem/,$p' | cut -c1-150
section "$OUT/case3-water/stage3.log" 4 | cut -c1-150
section "$OUT/case3-water/stage3.log" 10 | cut -c1-220
pause


banner "Case 4: the charges, the cross-checks, the dipole, the files"
cat <<'EOF'
[5]: per atom, Voronoi and Hirshfeld as SIESTA computes them, simple Hirshfeld (pass 0) and Hirshfeld-I --
2.7x larger here -- and the reference each atom ended with. [7]: the charges add up, every electron is
assigned. [8]: the dipole of the point charges of each method against SIESTA's dipole of the density.
EOF
section "$OUT/case3-water/stage3.log" 5 | cut -c1-175
section "$OUT/case3-water/stage3.log" 7 | cut -c1-175
section "$OUT/case3-water/stage3.log" 8 | cut -c1-150
echo "\$ ls    # --save-data, --png, --save-report"
ls -1 "$OUT/case3-water" | grep -v "^h2o$\|stage3.log" | sed 's/^/  /'
cp "$OUT/case3-water/hirshfeld_analysis.png" "$DIR/figures/h2o_dashboard.png"
echo "  (figures/h2o_dashboard.png is this dashboard; hirshfeld_charges.xyz opens in 'ase gui' coloured by charge)"
pause


banner "Case 5: how converged is 'converged'? Stage 3 again, nothing else"
cat <<'EOF'
The warning of Case 3 proposed --tol 0.0024. The iteration is cheap (seconds) and needs no new SIESTA run.
Below: the suggested --tol and a tight 1e-6; the estimate printed at the default tolerance (0.0079 e) is
compared with what the converged charges show. study_table.py [C] lists step, estimate and true error round
by round for all four systems.
EOF
for t in 0.0024 1e-6; do
    mkdir -p "$OUT/case5-tol/t$t"
    (cd "$OUT/case5-tol/t$t" && stb-hirshfeldAnalysis -O ../../case3-water/h2o --tol $t --max-iter 200 --save-data \
        --no-intro > log 2>&1)
    printf "  --tol %-7s %s   q(O) = %s\n" "$t" "$(strip_colors < "$OUT/case5-tol/t$t/log" | grep -o 'Converged after [0-9]* iteration(s)')" \
        "$(awk '!/^#/ && $2=="O" {print $10}' "$OUT/case5-tol/t$t/hirshfeld_charges.dat")"
done
echo "  --tol 0.005   (Case 3)                       q(O) = $(awk '!/^#/ && $2=="O" {print $10}' "$OUT/case3-water/hirshfeld_charges.dat")"
echo
python3 study_table.py 2>&1 | sed -n '/^\[C\] H2O/,/^$/p' | head -16
pause


banner "Case 6: N2O -- two atoms of one element, opposite signs"
cat <<'EOF'
In N=N=O the central nitrogen is positive and the terminal one negative. Each atom chooses its own ion every
round, so the terminal N mixes in N- and the central N mixes in N+: Stage 2 wrote both for this reason ([6]
marks such a species). The dipole of N2O is small (0.27 D here, 0.16 D measured) and the point charges
of every method miss it, even its sign: the atomic dipoles of Eq. 11 dominate.
EOF
mkdir -p "$OUT/case6-n2o"
unpack n2o "$OUT/case6-n2o/n2o"
(cd "$OUT/case6-n2o" && stb-hirshfeldAnalysis -O n2o --png --no-intro > stage3.log 2>&1)
section "$OUT/case6-n2o/stage3.log" 5 | cut -c1-175
section "$OUT/case6-n2o/stage3.log" 6 | cut -c1-140
section "$OUT/case6-n2o/stage3.log" 8 | cut -c1-150
cp "$OUT/case6-n2o/hirshfeld_analysis.png" "$DIR/figures/n2o_dashboard.png"
pause


banner "Case 7: NaCl -- a crystal, every periodic image, and Bader for comparison"
cat <<'EOF'
The primitive cell of rock salt has 2 atoms and is 4 Ang across; every point sees several images of each
atom (6 Na around each Cl), and the promolecule sums all of them within reach of each reference. The Na
pseudopotential has 9 valence electrons (2s2 2p6 3s1), so Na+ keeps 8. Simple Hirshfeld agrees with SIESTA's
own Hirshfeld; Hirshfeld-I is 3.4x larger; Bader (stb-bader on the same density) is larger still.
EOF
mkdir -p "$OUT/case7-nacl"
unpack nacl "$OUT/case7-nacl/nacl"
(cd "$OUT/case7-nacl" && stb-hirshfeldAnalysis -O nacl --png --no-intro > stage3.log 2>&1)
section "$OUT/case7-nacl/stage3.log" 1 | grep -E "Periodicity|Density grid|Valence" | cut -c1-150
section "$OUT/case7-nacl/stage3.log" 5 | cut -c1-175
echo "--- Bader charges of the same combined/ density (reference_results/nacl/bader/stb_bader_report.txt):"
strip_colors < "$REF/nacl/bader/stb_bader_report.txt" | sed -n '/^Idx/,/^Total/p' | cut -c1-90
cp "$OUT/case7-nacl/hirshfeld_analysis.png" "$DIR/figures/nacl_dashboard.png"
pause


banner "Case 8: MgO -- charges beyond +-1"
cat <<'EOF'
Hirshfeld-I with references up to +-1 cannot follow an atom past |q| = 1: the reference stops at Mg+ / O-.
MgO goes there (+-1.09): Stage 3 marks the atoms '(capped)' and warns. The fully self-consistent value would
need Mg2+ and O2- -- and O2- is not bound in vacuum, the known ambiguity of Hirshfeld-I for oxides.
EOF
mkdir -p "$OUT/case8-mgo"
unpack mgo "$OUT/case8-mgo/mgo"
(cd "$OUT/case8-mgo" && stb-hirshfeldAnalysis -O mgo --png --no-intro > stage3.log 2>&1)
section "$OUT/case8-mgo/stage3.log" 4 | grep -E "^Converged|Contraction" | cut -c1-150
section "$OUT/case8-mgo/stage3.log" 5 | cut -c1-175
section "$OUT/case8-mgo/stage3.log" 10 | cut -c1-220
cp "$OUT/case8-mgo/hirshfeld_analysis.png" "$DIR/figures/mgo_dashboard.png"
pause


banner "Case 9: what the charges depend on -- mesh, box, basis (water)"
cat <<'EOF'
Eleven more complete studies of water: MeshCutoff 200-600 Ry, the isolated-atom box 10-28 Ang, and six basis
sets. The mesh and the box hardly matter (SIESTA's atoms vanish beyond their orbital radii, well inside any of
these boxes). The basis matters: it shapes both the density and the reference atoms.
EOF
echo "\$ python3 study_table.py --png"
python3 study_table.py --png 2>&1 | sed -n '/^\[D\]/,$p'
pause


banner "Case 10: guard rails -- a missing reference, a short iteration, a copied study"
cat <<'EOF'
A reference whose SIESTA run is missing stops Stage 3 with a clear error (the reference cannot be guessed).
--max-iter too small is reported, not hidden. And a study copied with 'cp -r' reads its OWN combined/ --
the manifests store paths relative to the study folder.
EOF
mkdir -p "$OUT/case10-guards"
cd "$OUT/case10-guards"
unpack h2o no_anion; rm -f no_anion/ions/O/anion/siesta.RHO
set +e
echo "\$ stb-hirshfeldAnalysis -O no_anion          # ions/O/anion has no .RHO"
stb-hirshfeldAnalysis -O no_anion --no-intro 2>&1 | strip_colors | grep -E "ERROR" | cut -c1-200
echo "\$ stb-hirshfeldAnalysis -O h2o --max-iter 3"
unpack h2o h2o
stb-hirshfeldAnalysis -O h2o --max-iter 3 --no-intro 2>&1 | strip_colors | sed -n '/^\[10\] WARNINGS/,$p' | grep WARNING | cut -c1-200
echo "\$ cp -r h2o h2o_copy && stb-hirshfeldAnalysis -O h2o_copy"
cp -r h2o h2o_copy
stb-hirshfeldAnalysis -O h2o_copy --no-intro 2>&1 | strip_colors | grep -E "^Structure|INFO|Per species" | cut -c1-160
set -e
cd "$DIR"

echo
echo "Done. Every case wrote into $OUT/; the figures of the README are in figures/."
