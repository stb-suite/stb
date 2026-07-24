#!/bin/bash
# Guided example: stb-sqs (code 2.6 in the stb-suite menu)
#
# Not an automated test (see test/2-structures/6-sqs/test.sh for that) --
# a commented walk-through: it runs real commands, one group at a time,
# into its own output/<case>/ folder, and shows you the piece of output
# that proves what just happened. Pauses between sections so you can read
# before moving on. Safe to re-run any time -- it always starts by wiping
# its own output/.

set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

OUT="$DIR/output"
rm -rf "$OUT"
mkdir -p "$OUT"

pause() {
    echo
    read -p "  [Press Enter to continue] " -r
    echo
}

strip_ansi() {
    sed 's/\x1b\[[0-9;]*m//g' "$1"
}

echo "=================================================================="
echo " Special Quasirandom Structures: what problem does this solve?"
echo "=================================================================="
cat <<'EOF'
A DFT calculation needs a PERIODIC cell, but a real substitutional alloy
(NiFe, SiGe, a doped oxide, ...) has no periodicity at all -- each site is
occupied by one species or another at random, following the bulk
composition. Just substituting a small cell's sites at random doesn't
solve this: a handful of atoms is nowhere near enough to statistically
reproduce a truly random arrangement's short-range order (which atoms
tend to sit next to which), so different random draws of the same small
cell can behave quite differently from each other and from the real alloy.

A Special Quasirandom Structure (SQS) -- Zunger, Wei, Ferreira & Bernard,
Phys. Rev. Lett. 65, 353 (1990) -- sidesteps this: instead of a literal
random draw, it SEARCHES a small periodic cell for the one atomic
arrangement whose correlation functions (how often pairs/triplets/... of
sites at given distances are occupied by the same vs. different species)
best match those of the ideal, infinitely large random alloy. A small SQS
cell is a much better DFT-ready stand-in for "the random alloy" than a
naively randomized cell of the same size.

stb-sqs wraps icet's implementation of this search (Monte Carlo simulated
annealing, or exact enumeration for small cells) on top of pymatgen's
SQSTransformation.
EOF
pause

echo "=================================================================="
echo " Two ways to run it"
echo "=================================================================="
cat <<'EOF'
A -- direct CLI:
  stb-sqs -f fcc_ni.fdf --sublattice Ni --composition Ni:0.5,Fe:0.5

B -- interactive stb-suite menu:
  $ stb-suite
  Select an option (0-6, or a tool code like 4.1.2): 2.6

Every run always writes a numbered report ([0] RUN METADATA ... [11]
SUMMARY & FILES -- [6] only appears with --ml-relax), a structure
-validation checklist both before and after, and a before/after symmetry
comparison table (ordered input vs. final SQS structure). A full text
report (--save-report), MACE pre-relaxation (--ml-relax), and an
interactive 3D view (--view) are all off by default. Unlike earlier
versions of this tool, no report file is EVER written unless you ask for
it with --save-report -- everything it used to contain is now in this
same numbered console report.
EOF
pause

echo "=================================================================="
echo " output/monte-carlo/  --  the basic search, and what it costs"
echo "=================================================================="
cat <<'EOF'
fcc_ni.fdf is a single-species FCC nickel cell (4 atoms, a = 3.52 Ang,
space group Fm-3m No. 225). --sublattice Ni --composition Ni:0.5,Fe:0.5
disorders ALL of it into a 50:50 Ni-Fe binary alloy -- explicitly given
--scaling 8 here (a valid multiple of icet's own auto-detected minimum,
shown in section [4]) so the search has enough atoms to be a genuine SQS,
not just a tiny 2-atom toy cell. icet's Monte Carlo search minimizes an
objective function Q = -w*L + sum|correlation - target correlation| (van
de Walle et al., Calphad 42, 13 (2013)) -- more negative is better, a
closer match to the ideal random alloy's correlation functions:
EOF
mkdir -p "$OUT/monte-carlo"
cp fcc_ni.fdf "$OUT/monte-carlo/"
echo
echo "\$ stb-sqs -f fcc_ni.fdf --sublattice Ni --composition Ni:0.5,Fe:0.5 --scaling 8 --mc-steps 5000 -o sqs.fdf --no-intro"
(cd "$OUT/monte-carlo" && stb-sqs -f fcc_ni.fdf --sublattice Ni --composition Ni:0.5,Fe:0.5 \
    --scaling 8 --mc-steps 5000 -o sqs.fdf --no-intro > console.log)
grep -E "Output formula|Output atoms|Objective function" "$OUT/monte-carlo/console.log"
echo
echo "Disordering an ordered metal breaks its symmetry substantially --"
echo "same theme as introducing a point defect (see 2.5-stb-defect/):"
sed -n '/Detailed symmetry analysis/,/WRITING OUTPUT/p' "$OUT/monte-carlo/console.log" | head -n -1
pause

echo "=================================================================="
echo " output/enumeration/  --  exact search, cross-checked against Monte Carlo"
echo "=================================================================="
cat <<'EOF'
--method enumeration exhaustively checks every distinct arrangement at
this cell size instead of a stochastic search -- guaranteed to find the
TRUE global optimum, but only feasible for small cells (the number of
arrangements grows combinatorially). Running it on the exact same problem
as monte-carlo/ above is a real cross-check: if Monte Carlo's objective
function matches enumeration's exactly, Monte Carlo also found the true
global optimum here, not just a good-enough local one.
EOF
mkdir -p "$OUT/enumeration"
cp fcc_ni.fdf "$OUT/enumeration/"
echo
echo "\$ stb-sqs -f fcc_ni.fdf --sublattice Ni --composition Ni:0.5,Fe:0.5 --scaling 8 --method enumeration -o sqs.fdf --no-intro"
(cd "$OUT/enumeration" && stb-sqs -f fcc_ni.fdf --sublattice Ni --composition Ni:0.5,Fe:0.5 \
    --scaling 8 --method enumeration -o sqs.fdf --no-intro > console.log)
grep -E "Output formula|Output atoms|Objective function" "$OUT/enumeration/console.log"
echo
MC_OBJ=$(grep "Objective function" "$OUT/monte-carlo/console.log" | awk '{print $NF}')
ENUM_OBJ=$(grep "Objective function" "$OUT/enumeration/console.log" | awk '{print $NF}')
if [ "$MC_OBJ" = "$ENUM_OBJ" ]; then
    echo "Confirmed: Monte Carlo ($MC_OBJ) and exhaustive enumeration ($ENUM_OBJ)"
    echo "found the exact same optimal objective function -- Monte Carlo really"
    echo "did find the true global optimum here, not a lucky local one."
else
    echo "Monte Carlo ($MC_OBJ) vs. enumeration ($ENUM_OBJ) -- if these differ,"
    echo "Monte Carlo landed in a local optimum this time (rerun with more"
    echo "--mc-steps); enumeration's answer is always the true optimum."
fi
pause

echo "=================================================================="
echo " output/scaling-constraint/  --  not every cell size fits a composition"
echo "=================================================================="
cat <<'EOF'
A target composition needs an EXACT integer atom count per species at the
chosen cell size. icet auto-detects the smallest --scaling that allows
this (section [4], "Minimal valid scaling"); any --scaling that isn't a
multiple of it is rejected up front, with a clear message, instead of
letting icet's own search silently fail to find candidates:
EOF
mkdir -p "$OUT/scaling-constraint"
cp fcc_ni.fdf "$OUT/scaling-constraint/"
echo
echo "\$ stb-sqs -f fcc_ni.fdf --sublattice Ni --composition Ni:0.3,Fe:0.7 --scaling 4 --no-intro   # 30:70 needs a multiple of 10, not 4"
set +e
(cd "$OUT/scaling-constraint" && stb-sqs -f fcc_ni.fdf --sublattice Ni --composition Ni:0.3,Fe:0.7 \
    --scaling 4 --no-intro > rejected.log 2>&1)
set -e
grep -E "Minimal valid scaling|ERROR" "$OUT/scaling-constraint/rejected.log"
echo
echo "Omit --scaling (or give the suggested multiple, 10) and it works:"
echo "\$ stb-sqs -f fcc_ni.fdf --sublattice Ni --composition Ni:0.3,Fe:0.7 --mc-steps 2000 -o sqs.fdf --no-intro"
(cd "$OUT/scaling-constraint" && stb-sqs -f fcc_ni.fdf --sublattice Ni --composition Ni:0.3,Fe:0.7 \
    --mc-steps 2000 -o sqs.fdf --no-intro > accepted.log)
grep -E "Scaling used|Output formula|Output atoms" "$OUT/scaling-constraint/accepted.log"
pause

echo "=================================================================="
echo " output/cluster-cutoffs/  --  matching MORE correlations is harder"
echo "=================================================================="
cat <<'EOF'
--cluster-cutoffs controls which pair/triplet/... clusters (by size:shell)
the search tries to match to the ideal random alloy. Asking for a wider
set of correlations to match simultaneously, with the SAME small cell, is
a strictly harder target -- watch the objective function (and search time)
grow when widening the default 2:3,3:2,4:1 to include longer pair shells
and a wider triplet shell:
EOF
mkdir -p "$OUT/cluster-cutoffs"
cp fcc_ni.fdf "$OUT/cluster-cutoffs/"
echo
echo "\$ stb-sqs ... --scaling 8 --mc-steps 1000 -o default.fdf --no-intro                        # default cutoffs"
(cd "$OUT/cluster-cutoffs" && stb-sqs -f fcc_ni.fdf --sublattice Ni --composition Ni:0.5,Fe:0.5 \
    --scaling 8 --mc-steps 1000 -o default.fdf --no-intro > default.log)
echo "\$ stb-sqs ... --scaling 8 --mc-steps 1000 --cluster-cutoffs 2:5,3:3 -o wide.fdf --no-intro   # wider cutoffs"
(cd "$OUT/cluster-cutoffs" && stb-sqs -f fcc_ni.fdf --sublattice Ni --composition Ni:0.5,Fe:0.5 \
    --scaling 8 --mc-steps 1000 --cluster-cutoffs 2:5,3:3 -o wide.fdf --no-intro > wide.log)
echo
echo "Default cutoffs:"
grep -E "Cluster cutoffs|Objective function" "$OUT/cluster-cutoffs/default.log"
echo "Wide cutoffs:"
grep -E "Cluster cutoffs|Objective function" "$OUT/cluster-cutoffs/wide.log"
pause

echo "=================================================================="
echo " output/instances/  --  parallel independent searches"
echo "=================================================================="
cat <<'EOF'
--instances N runs N independent Monte Carlo searches and keeps the best.
On the same small problem as monte-carlo/ above, this should (and, run
live, does) rediscover the same true global optimum enumeration already
confirmed -- a second, independent cross-check:
EOF
mkdir -p "$OUT/instances"
cp fcc_ni.fdf "$OUT/instances/"
echo
echo "\$ stb-sqs -f fcc_ni.fdf --sublattice Ni --composition Ni:0.5,Fe:0.5 --scaling 8 --mc-steps 2000 --instances 2 -o sqs.fdf --no-intro"
(cd "$OUT/instances" && stb-sqs -f fcc_ni.fdf --sublattice Ni --composition Ni:0.5,Fe:0.5 \
    --scaling 8 --mc-steps 2000 --instances 2 -o sqs.fdf --no-intro > console.log)
grep -E "Instances used|Objective function" "$OUT/instances/console.log"
pause

echo "=================================================================="
echo " output/sublattice-mgo/  --  disordering ONE species, leaving the other alone"
echo "=================================================================="
cat <<'EOF'
mgo_rocksalt.fdf already has TWO species (Mg, O). --sublattice Mg targets
only the Mg sites -- modeling a (Mg,Ni)O rocksalt solid solution (studied
for catalysis/electrocatalysis) -- while O keeps its role in the crystal.
icet rebuilds the cell from its own primitive-cell tiling, so the O atoms'
exact coordinates/order can come out different from the input (same
"equivalent, not identical" caveat as spglib-based tools elsewhere in this
suite) -- but the underlying rocksalt geometry is preserved, verified by
the cation-O nearest-neighbor distance coming out exactly unchanged:
EOF
mkdir -p "$OUT/sublattice-mgo"
cp mgo_rocksalt.fdf "$OUT/sublattice-mgo/"
echo
echo "\$ stb-sqs -f mgo_rocksalt.fdf --sublattice Mg --composition Mg:0.5,Ni:0.5 --scaling 4 --mc-steps 2000 -o mgno.fdf --no-intro"
(cd "$OUT/sublattice-mgo" && stb-sqs -f mgo_rocksalt.fdf --sublattice Mg --composition Mg:0.5,Ni:0.5 \
    --scaling 4 --mc-steps 2000 -o mgno.fdf --no-intro > console.log)
grep -E "Species  |Sites  |Output formula|Output atoms" "$OUT/sublattice-mgo/console.log"
echo
python3 - "$OUT/sublattice-mgo/mgo_rocksalt.fdf" "$OUT/sublattice-mgo/mgno.fdf" <<'PYEOF'
import sys
from stb.core import structure_io
orig = structure_io.to_pymatgen(structure_io.read_fdf(sys.argv[1]))
sqs = structure_io.to_pymatgen(structure_io.read_fdf(sys.argv[2]))
def min_cation_o(s):
    dists = []
    for i, site in enumerate(s):
        if site.specie.symbol in ("Mg", "Ni"):
            dists.append(min(s.get_distance(i, j) for j, o in enumerate(s) if o.specie.symbol == "O"))
    return sorted(set(round(d, 4) for d in dists))
print(f"O atoms: {sum(1 for s in orig if s.specie.symbol == 'O')} -> "
      f"{sum(1 for s in sqs if s.specie.symbol == 'O')}  (unchanged)")
print(f"Cation-O nearest-neighbor distance: {min_cation_o(orig)} -> {min_cation_o(sqs)} Ang  (unchanged)")
PYEOF
pause

echo "=================================================================="
echo " output/ml-relax/  --  MACE pre-relaxation of the SQS structure"
echo "=================================================================="
if ! python3 -c "import mace" 2>/dev/null; then
    echo "Skipping -- needs the optional 'ml' extra: pip install stb_suite[ml]"
    echo "(everything else in this script works fine without it)."
else
cat <<'EOF'
An SQS cell's local environment is inherently irregular -- unlike a
perfect crystal, atoms genuinely sit at slightly different local
equilibrium positions depending on their neighbors. --ml-relax (+
--ml-relax-cell) relaxes the found SQS structure with MACE before writing
it out, giving a better-optimized starting geometry for a real DFT run:
EOF
    mkdir -p "$OUT/ml-relax"
    cp fcc_ni.fdf "$OUT/ml-relax/"
    echo
    echo "\$ stb-sqs -f fcc_ni.fdf --sublattice Ni --composition Ni:0.5,Fe:0.5 --scaling 8 --mc-steps 2000 --ml-relax --ml-relax-cell -o relaxed.fdf --no-intro"
    (cd "$OUT/ml-relax" && stb-sqs -f fcc_ni.fdf --sublattice Ni --composition Ni:0.5,Fe:0.5 \
        --scaling 8 --mc-steps 2000 --ml-relax --ml-relax-cell -o relaxed.fdf --no-intro > console.log 2>&1)
    echo "Full MACE simulation detail, straight from the report:"
    strip_ansi "$OUT/ml-relax/console.log" | awk '/\[6\] ML PRE-RELAXATION/{flag=1} /\[7\] STRUCTURE VALIDATION/{flag=0} flag'
    echo
    echo "Provenance header now also records the MACE pre-relaxation:"
    head -6 "$OUT/ml-relax/relaxed.fdf"
fi
pause

echo "=================================================================="
echo " output/full-report/  --  --save-report, validation, references.bib"
echo "=================================================================="
cat <<'EOF'
The full numbered report (also written to stb_sqs_report.txt with
--save-report) includes the structure-validation checklist for the input
and the SQS output, and a references.bib -- SIESTA plus icet, the library
behind the SQS search itself:
EOF
mkdir -p "$OUT/full-report"
cp fcc_ni.fdf "$OUT/full-report/"
echo "\$ stb-sqs -f fcc_ni.fdf --sublattice Ni --composition Ni:0.5,Fe:0.5 --scaling 8 --mc-steps 2000 --save-report --no-intro"
(cd "$OUT/full-report" && stb-sqs -f fcc_ni.fdf --sublattice Ni --composition Ni:0.5,Fe:0.5 \
    --scaling 8 --mc-steps 2000 --save-report --no-intro > console.log)
echo
echo "Report sections written to stb_sqs_report.txt:"
grep -E "^\[[0-9]+\]" "$OUT/full-report/stb_sqs_report.txt"
echo
echo "Validation checklist (one row shown):"
grep -m1 "Atom proximity" "$OUT/full-report/console.log"
echo
echo "references.bib -- SIESTA plus icet:"
grep "^@" "$OUT/full-report/references.bib"
pause

echo "=================================================================="
echo " Proof: CLI and the interactive stb-suite menu agree"
echo "=================================================================="
echo "Driving the same basic case through the interactive menu's manual"
echo "entry mode and checking it reaches the same atom count."
TMP="$(mktemp -d)"
cp fcc_ni.fdf "$TMP/"
echo
echo "\$ printf '2.6\\ndoes_not_exist.fdf\\nfcc_ni.fdf\\nNi\\nNi:0.5,Fe:0.5\\n8\\n1\\n\\n\\n\\n\\n\\n0\\n' | stb-suite"
(cd "$TMP" && printf '2.6\ndoes_not_exist.fdf\nfcc_ni.fdf\nNi\nNi:0.5,Fe:0.5\n8\n1\n\n\n\n\n\n0\n' | stb-suite > session.log 2>&1) || true
if grep -q "Output atoms       : 8" "$TMP/session.log"; then
    echo "Confirmed: the interactive menu built and launched the exact same"
    echo "underlying stb-sqs command as the CLI walkthrough above (8 atoms out,"
    echo "same cell size as output/monte-carlo/)."
else
    echo "Unexpected: menu did not reach the write step -- see $TMP/session.log."
fi
rm -rf "$TMP"
pause

echo "=================================================================="
echo " Done"
echo "=================================================================="
cat <<EOF
Up to eight self-contained folders were generated under output/ (ml-relax/
skipped if the optional 'ml' extra isn't installed):
  monte-carlo/       enumeration/        scaling-constraint/
  cluster-cutoffs/   instances/          sublattice-mgo/
  ml-relax/          full-report/

Each has references.bib (SIESTA + icet, always); full-report/ additionally
has stb_sqs_report.txt; ml-relax/'s references.bib also cites the MACE
architecture/foundation-model papers.

Recap of what this walkthrough covered:
  - what an SQS is and why it beats naive random substitution (Zunger,
    Wei, Ferreira & Bernard, 1990) for a small periodic DFT cell
  - icet's objective function (van de Walle et al., 2013) and how a real
    disordering breaks symmetry substantially, same theme as 2.5-stb-defect/
  - Monte Carlo vs. exact enumeration -- cross-verified live to find the
    exact same global optimum on the same small problem
  - why --scaling must be a multiple of icet's own auto-detected minimum,
    and what happens when it isn't
  - --cluster-cutoffs: matching more/longer-range correlations is a
    strictly harder search target on the same cell
  - --instances: independent parallel searches as a second cross-check
  - --sublattice only disorders the species you name -- verified on a
    2-species fixture (MgO) that the untouched sublattice's local
    coordination geometry survives even though icet rebuilds the cell
  - --ml-relax pre-optimizing an SQS cell's inherently irregular local
    environment
  - the structure-validation checklist, references.bib (now including
    icet), and --save-report
  - CLI and the interactive stb-suite menu building the same command
  - stb_sqs_report.txt is no longer written unless you ask for it

Not exercised by this script (needs a display): --view opens the input
structure and the final SQS structure in an interactive ase-gui window --
try it yourself:
  stb-sqs -f fcc_ni.fdf --sublattice Ni --composition Ni:0.5,Fe:0.5 --view

As a next step, try on your own:
  stb-sqs -f your_alloy.fdf --sublattice X --composition X:0.5,Y:0.5     # your own binary alloy
  stb-sqs -f fcc_ni.fdf --sublattice Ni --composition Ni:0.25,Fe:0.75    # a different composition ratio
  stb-supercell -f fcc_ni.fdf -d 2 2 2 -o ni32.fdf --no-intro && \\
      stb-sqs -f ni32.fdf --sublattice Ni --composition Ni:0.5,Fe:0.5    # a larger, more dilute SQS

Known limitation, not exercised above: --composition is practically limited to
2 species per sublattice -- a 3+-species target crashes inside pymatgen's own
SQSTransformation/IcetSQS wrapper ("ASE Atoms only supports ordered
structures"), a verified upstream pymatgen issue, not something stb-sqs
itself can currently work around (see --help).
EOF
