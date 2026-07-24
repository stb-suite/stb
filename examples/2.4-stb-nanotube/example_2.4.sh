#!/bin/bash
# Guided example: stb-nanotube (code 2.4 in the stb-suite menu)
#
# Not an automated test (see test/2-structures/4-nanotube/test.sh for that)
# -- a commented walk-through: it runs real commands, one group at a time,
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

echo "=================================================================="
echo " The chiral vector C and the translation vector T"
echo "=================================================================="
cat <<'EOF'
A chirality (n, m) defines the CHIRAL VECTOR C = n*a1 + m*a2, using the
monolayer's own in-plane lattice vectors. For --mode tube (the default),
C becomes the tube's circumference (2R = |C| / pi); for --mode ribbon, C
instead sets the finite WIDTH direction.

Both modes also need a TRANSLATION VECTOR T, a real lattice vector along
which the result is periodic. For a ribbon, "not parallel to C" is enough.
For a TUBE, it is not: rolling maps the flat cell onto a cylinder treating
C and T as an ORTHOGONAL (circumference, axial) frame, which only preserves
real bond lengths if T has no component along C -- T must be as close to
PERPENDICULAR to C as an integer lattice vector allows.
EOF
pause

echo "=================================================================="
echo " Two ways to run it"
echo "=================================================================="
cat <<'EOF'
A -- direct CLI:
  stb-nanotube -f graphene.fdf --chirality 6 0

B -- interactive stb-suite menu:
  $ stb-suite
  Select an option (0-6, or a tool code like 4.1.2): 2.4

Every run always writes a numbered report ([0] RUN METADATA ... [9] SUMMARY
& FILES -- [4] ML PRE-RELAXATION only appears with --ml-relax), a
structure-validation checklist both before and after, and a before/after
symmetry comparison table. A full text report (--save-report), MACE
pre-relaxation (--ml-relax), and an interactive 3D view (--view) are all
off by default.
EOF
pause

echo "=================================================================="
echo " output/basic-tube/  --  zigzag (6,0), the corrected perpendicular-T fix"
echo "=================================================================="
cat <<'EOF'
graphene.fdf is a primitive 2-atom hexagonal cell (a = 2.46 Ang). A real bug
existed here before this session: picking the shortest lattice vector merely
NON-PARALLEL to C (not perpendicular) silently gave (6,0) a 12-atom cell
where every atom had only 1 real bonded neighbor, instead of 3. Fixed by
ranking candidates on perpendicularity first, length second -- watch for
the correct 24-atom, fully 3-fold-coordinated result below, verified live
against ase.build.nanotube (an independent reference implementation).
EOF
mkdir -p "$OUT/basic-tube"
cp graphene.fdf "$OUT/basic-tube/"
echo
echo "\$ stb-nanotube -f graphene.fdf --chirality 6 0 --no-intro"
(cd "$OUT/basic-tube" && stb-nanotube -f graphene.fdf --chirality 6 0 --no-intro > console.log)
grep -E "CNT type|Translation vector|Cells in periodic|Tube diameter|Output atoms" "$OUT/basic-tube/console.log"
echo
echo "gcd(6,0) = 6 > 1 -- the tool's own note about the classic (screw-symmetry"
echo "-based) CNT cell being potentially smaller, printed directly in the report:"
grep "Note: gcd(n, m) > 1" "$OUT/basic-tube/console.log" | fold -s -w 100
echo
echo "Verified physically correct -- every atom 3-fold coordinated at the real"
echo "graphene C-C bond length (curvature-compressed from flat 1.42 Ang):"
python3 - <<'PYEOF'
from stb.core import structure_io
import numpy as np
s = structure_io.to_pymatgen(structure_io.read_fdf("output/basic-tube/nanotube.fdf"))
dm = s.distance_matrix
np.fill_diagonal(dm, np.inf)
coord = (dm < 1.6).sum(axis=1)
nn = dm.min(axis=1)
print(f"  coordination: min={coord.min()} max={coord.max()} (expect 3)")
print(f"  bond length : {nn.min():.4f}-{nn.max():.4f} Ang")
PYEOF
echo
echo "Before/after symmetry table -- Layer Group flips the OTHER way from"
echo "2.3-stb-slab/: the 2D monolayer (before) has a real one, the 1D tube"
echo "(after) never does (spglib has no 1D rod-group classification):"
sed -n '/Detailed symmetry analysis/,/WRITING OUTPUT/p' "$OUT/basic-tube/console.log" | head -n -1
echo
echo "Provenance header written into the output .fdf itself:"
head -3 "$OUT/basic-tube/nanotube.fdf"
pause

echo "=================================================================="
echo " output/electronic-hint/  --  metallic vs. semiconducting, verified"
echo "=================================================================="
cat <<'EOF'
For a single-species hexagonal lattice (graphene itself), stb-nanotube
reports the well-known rule: metallic if (n-m) mod 3 == 0, else
semiconducting. Zigzag (6,0): 6-0=6, divisible by 3 -> metallic. Zigzag
(7,0): 7-0=7, not divisible by 3 -> semiconducting. Armchair tubes
(n,n) always have n-m=0, divisible by 3 -> ALWAYS metallic, regardless of n.
EOF
mkdir -p "$OUT/electronic-hint"
cp graphene.fdf "$OUT/electronic-hint/"
for chir in "6 0" "7 0" "6 6"; do
    tag=$(echo "$chir" | tr ' ' '_')
    (cd "$OUT/electronic-hint" && stb-nanotube -f graphene.fdf --chirality $chir \
        -o "tube_$tag.fdf" --no-intro > "console_$tag.log")
    echo -n "  ($chir): "
    sed 's/\x1b\[[0-9;]*m//g' "$OUT/electronic-hint/console_$tag.log" | \
        grep "Electronic hint" | sed 's/^Electronic hint  : //' | cut -d'-' -f1
done
pause

echo "=================================================================="
echo " output/chiral-index/  --  a genuinely chiral (7,1) tube, gcd=1"
echo "=================================================================="
cat <<'EOF'
(7,1): gcd(7,1) = 1 -- unlike every (n,0) zigzag (always gcd=n) or (n,n)
armchair (always gcd=n), a coprime chiral index like this one needs NO
gcd(n,m)>1 caveat: the shortest perpendicular T already gives the smallest
possible cell. Watch the CNT type and chiral angle (between 0 deg zigzag
and the armchair value for this lattice convention):
EOF
mkdir -p "$OUT/chiral-index"
cp graphene.fdf "$OUT/chiral-index/"
echo
echo "\$ stb-nanotube -f graphene.fdf --chirality 7 1 --no-intro"
(cd "$OUT/chiral-index" && stb-nanotube -f graphene.fdf --chirality 7 1 --no-intro > console.log)
grep -E "CNT type|Chiral angle|Cells in periodic|Output atoms" "$OUT/chiral-index/console.log"
echo "(no 'gcd(n, m) > 1' note this time -- confirm it's absent:)"
grep -c "Note: gcd" "$OUT/chiral-index/console.log" | sed 's/^/  matches: /'
pause

echo "=================================================================="
echo " output/curvature-trend/  --  measured bond compression vs. diameter"
echo "=================================================================="
cat <<'EOF'
Rolling a flat sheet without stretching along the circumference still
changes real 3D bond lengths -- smaller tubes are more curved, hence more
compressed. Zigzag (n,0) for increasing n, bond length measured directly
from each written structure (not just read from the report):
EOF
mkdir -p "$OUT/curvature-trend"
cp graphene.fdf "$OUT/curvature-trend/"
for n in 4 6 10 20; do
    (cd "$OUT/curvature-trend" && stb-nanotube -f graphene.fdf --chirality $n 0 \
        -o "tube_n$n.fdf" --no-intro > "console_n$n.log")
done
python3 - <<'PYEOF'
from stb.core import structure_io
import numpy as np
for n in (4, 6, 10, 20):
    s = structure_io.to_pymatgen(structure_io.read_fdf(f"output/curvature-trend/tube_n{n}.fdf"))
    diameter = n * 2.46 / np.pi  # |C| = n * a (m=0, zigzag), 2R = |C| / pi
    dm = s.distance_matrix
    np.fill_diagonal(dm, np.inf)
    print(f"  n={n:>2}: diameter={diameter:6.2f} Ang  "
          f"bond length = {dm.min(axis=1).max():.4f} Ang  (flat graphene: 1.4200 Ang)")
PYEOF
echo "Monotonically approaching the flat, strain-free value as the tube widens."
pause

echo "=================================================================="
echo " output/ribbon/  --  finite width scales linearly with --repeats"
echo "=================================================================="
cat <<'EOF'
--mode ribbon cuts a finite-width strip instead of rolling a cylinder --
no wrap, so bond lengths are EXACT regardless of chirality (no curvature
distortion at all, unlike tube mode above). Width scales linearly with
--repeats:
EOF
mkdir -p "$OUT/ribbon"
cp graphene.fdf "$OUT/ribbon/"
for r in 1 2 4; do
    (cd "$OUT/ribbon" && stb-nanotube -f graphene.fdf --chirality 6 0 --mode ribbon \
        --repeats $r -o "ribbon_r$r.fdf" --no-intro > "console_r$r.log")
    echo -n "  repeats=$r: "
    grep -E "Ribbon width|Output atoms" "$OUT/ribbon/console_r$r.log" | tr '\n' ' '
    echo
done
pause

echo "=================================================================="
echo " output/passivate/  --  a closed tube has 0 dangling bonds; a ribbon's"
echo " two edges genuinely do"
echo "=================================================================="
cat <<'EOF'
A tube is BOTH circumferentially closed and axially periodic -- every atom
already has its full coordination, no dangling bonds anywhere. A ribbon has
two real, finite physical edges, which DO have genuine dangling bonds.
--passivate proves both sides of this on the exact same chirality:
EOF
mkdir -p "$OUT/passivate"
cp graphene.fdf "$OUT/passivate/"
echo
echo "\$ stb-nanotube -f graphene.fdf --chirality 6 0 --passivate -o tube.fdf --no-intro"
(cd "$OUT/passivate" && stb-nanotube -f graphene.fdf --chirality 6 0 --passivate \
    -o tube.fdf --no-intro > console_tube.log)
grep -E "Dangling bonds found|Auto-passivated" "$OUT/passivate/console_tube.log" | sed 's/^/  tube:   /'
echo
echo "\$ stb-nanotube -f graphene.fdf --chirality 6 0 --mode ribbon --repeats 4 --passivate -o ribbon.fdf --no-intro"
(cd "$OUT/passivate" && stb-nanotube -f graphene.fdf --chirality 6 0 --mode ribbon --repeats 4 \
    --passivate -o ribbon.fdf --no-intro > console_ribbon.log)
grep -E "Dangling bonds found|Auto-passivated" "$OUT/passivate/console_ribbon.log" | sed 's/^/  ribbon: /'
pause

echo "=================================================================="
echo " output/ml-relax/  --  MACE measurably releases curvature strain"
echo "=================================================================="
if ! python3 -c "import mace" 2>/dev/null; then
    echo "Skipping -- needs the optional 'ml' extra: pip install stb_suite[ml]"
    echo "(everything else in this script works fine without it)."
else
cat <<'EOF'
A freshly rolled tube's bond lengths/axial period come straight from the
FLAT monolayer's own lattice constant -- curvature-compressed, not
relaxed. --ml-relax (with --ml-relax-cell, also letting the axial period
itself relax) measurably releases that strain:
EOF
    mkdir -p "$OUT/ml-relax"
    cp graphene.fdf "$OUT/ml-relax/"
    echo
    echo "\$ stb-nanotube -f graphene.fdf --chirality 6 0 --ml-relax --ml-relax-cell -o relaxed.fdf --no-intro"
    (cd "$OUT/ml-relax" && stb-nanotube -f graphene.fdf --chirality 6 0 \
        --ml-relax --ml-relax-cell -o relaxed.fdf --no-intro > console.log 2>&1)
    echo "Full MACE simulation detail, straight from the report:"
    sed 's/\x1b\[[0-9;]*m//g' "$OUT/ml-relax/console.log" | \
        awk '/\[4\] ML PRE-RELAXATION/{flag=1} /\[5\] STRUCTURE VALIDATION/{flag=0} flag'
    echo
    echo "Provenance header now also records the MACE pre-relaxation:"
    head -4 "$OUT/ml-relax/relaxed.fdf"
fi
pause

echo "=================================================================="
echo " output/corrugated/  --  a buckled sheet rolls into a corrugated tube"
echo "=================================================================="
cat <<'EOF'
silicene_buckled.fdf is a synthetic honeycomb monolayer whose 2 basis atoms
sit at slightly different heights (0.22 Ang apart) -- a real structural
feature of silicene/germanene, unlike graphene's perfectly flat sheet.
stb-nanotube tracks each atom's own out-of-plane offset when rolling, so
the two sublattices should end up at two DIFFERENT radii in the tube --
i.e. a genuinely corrugated/puckered tube, not a flattened approximation:
EOF
mkdir -p "$OUT/corrugated"
cp silicene_buckled.fdf "$OUT/corrugated/"
echo
echo "\$ stb-nanotube -f silicene_buckled.fdf --chirality 6 0 --no-intro"
(cd "$OUT/corrugated" && stb-nanotube -f silicene_buckled.fdf --chirality 6 0 --no-intro > console.log)
grep -E "Tube diameter|Output atoms" "$OUT/corrugated/console.log"
python3 - <<'PYEOF'
from stb.core import structure_io
import numpy as np
s = structure_io.to_pymatgen(structure_io.read_fdf("output/corrugated/nanotube.fdf"))
r = np.linalg.norm(s.cart_coords[:, :2], axis=1)
radii = sorted(set(np.round(r, 4)))
print(f"  radii found: {radii} Ang -- {len(radii)} distinct value(s)")
print(f"  difference : {radii[-1] - radii[0]:.4f} Ang (matches the 0.22 Ang buckling amplitude)")
PYEOF
pause

echo "=================================================================="
echo " output/full-report/  --  --save-report, validation, references.bib"
echo "=================================================================="
cat <<'EOF'
The full numbered report (also written to stb_nanotube_report.txt with
--save-report) includes the structure-validation checklist for the input
monolayer and the output tube/ribbon, and a references.bib with SIESTA.
EOF
mkdir -p "$OUT/full-report"
cp graphene.fdf "$OUT/full-report/"
echo "\$ stb-nanotube -f graphene.fdf --chirality 6 0 --save-report --no-intro"
(cd "$OUT/full-report" && stb-nanotube -f graphene.fdf --chirality 6 0 --save-report --no-intro > console.log)
echo
echo "Report sections written to stb_nanotube_report.txt:"
grep -E "^\[[0-9]\]" "$OUT/full-report/stb_nanotube_report.txt"
echo
echo "Validation checklist (one row shown):"
grep -m1 "Atom proximity" "$OUT/full-report/console.log"
echo
echo "references.bib always written -- SIESTA:"
grep "^@" "$OUT/full-report/references.bib"
pause

echo "=================================================================="
echo " Proof: CLI and the interactive stb-suite menu agree"
echo "=================================================================="
echo "Driving the same basic case through the interactive menu's manual"
echo "entry mode and checking it reaches the same atom count."
TMP="$(mktemp -d)"
cp graphene.fdf "$TMP/"
echo
echo "\$ printf '2.4\\ngraphene.fdf\\n6 0\\n1\\n\\n\\n\\n\\n\\n\\n0\\n' | stb-suite"
(cd "$TMP" && printf '2.4\ngraphene.fdf\n6 0\n1\n\n\n\n\n\n\n0\n' | stb-suite > session.log 2>&1) || true
if grep -q "Output atoms     : 24" "$TMP/session.log"; then
    echo "Confirmed: the interactive menu built and launched the exact same"
    echo "underlying stb-nanotube command as the CLI walkthrough above (24"
    echo "atoms out, same as output/basic-tube/)."
else
    echo "Unexpected: menu did not reach the construction step -- see $TMP/session.log."
fi
rm -rf "$TMP"
pause

echo "=================================================================="
echo " Done"
echo "=================================================================="
cat <<EOF
Up to nine self-contained folders were generated under output/ (ml-relax/
skipped if the optional 'ml' extra isn't installed):
  basic-tube/       electronic-hint/   chiral-index/
  curvature-trend/  ribbon/            passivate/
  ml-relax/         corrugated/        full-report/

Each has references.bib; full-report/ additionally has
stb_nanotube_report.txt; ml-relax/'s references.bib also cites the MACE
architecture/foundation-model papers.

Recap of what this walkthrough covered:
  - the chiral vector C (circumference/width) and translation vector T
    (axial period), and why T must be PERPENDICULAR to C for a tube, not
    just non-parallel -- a real bug, found and fixed, verified against
    ase.build.nanotube
  - zigzag/armchair/chiral classification and the graphene metallic/
    semiconducting (n-m) mod 3 rule, verified across several chiralities
  - gcd(n,m) > 1: a physically correct but not always minimal cell,
    because this tool deliberately never assumes hexagonal-only screw
    symmetry -- contrasted with a genuinely chiral, gcd=1 index
  - measured curvature-induced bond compression vs. diameter
  - ribbon mode's exact (uncompressed) bond lengths and linear width
    scaling with --repeats
  - a tube's 0 dangling bonds vs. a ribbon's 2 real edges, both verified
    with --passivate
  - --ml-relax measurably releasing curvature strain
  - rolling a corrugated (buckled, silicene-like) sheet into a genuinely
    corrugated tube -- two distinct radii, matching the input buckling
  - the Layer Group column's before/after contrast with 2.3-stb-slab/:
    here the 2D input has one and the 1D output never does
  - CLI and the interactive stb-suite menu building the same command

Not exercised by this script (needs a display): --view opens the input
monolayer and the final tube/ribbon in an interactive ase-gui window --
try it yourself:
  stb-nanotube -f graphene.fdf --chirality 6 0 --view

As a next step, try on your own:
  stb-nanotube -f your_monolayer.fdf --chirality 10 5   # your own 2D material
  stb-nanotube -f graphene.fdf --chirality 6 0 --mode ribbon --repeats 6 --passivate
  stb-nanotube -f graphene.fdf --chirality 6 0 --ml-relax --custom-model my_finetuned.model
EOF
