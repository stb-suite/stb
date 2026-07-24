#!/bin/bash
# Guided example: stb-unitcell (code 2.7 in the stb-suite menu)
#
# Not an automated test (see test/2-structures/7-unitcell/test.sh for that)
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

strip_ansi() {
    sed 's/\x1b\[[0-9;]*m//g' "$1"
}

echo "=================================================================="
echo " What is a unit cell, and why does the choice matter?"
echo "=================================================================="
cat <<'EOF'
A crystal's periodicity can be described by more than one repeating cell --
all of them describe the exact same infinite structure, just tiled
differently:

PRIMITIVE    -- the smallest possible repeating unit. Fewer atoms means a
                cheaper DFT calculation (fewer basis functions, smaller
                Hamiltonian, faster SCF) -- the natural choice for a real
                calculation, especially phonon/band-structure work.
CONVENTIONAL -- a larger, standardized cell (the International Tables /
                IUCr convention). Easier for humans to read and compare
                against a database entry (most databases store this one),
                at the cost of more atoms than strictly necessary.
REFINED      -- conventional-sized, but with atomic positions snapped
                exactly onto the detected symmetry -- cleans up numerical
                noise from a relaxation or a hand-built/CIF structure
                WITHOUT changing which atoms are present.

stb-unitcell uses pymatgen's symmetry analyzer (spglib) to detect a
structure's true symmetry and rebuild whichever of these three cells you
ask for.
EOF
pause

echo "=================================================================="
echo " Two ways to run it"
echo "=================================================================="
cat <<'EOF'
A -- direct CLI:
  stb-unitcell -f fcc_ni.fdf --mode primitive

B -- interactive stb-suite menu:
  $ stb-suite
  Select an option (0-6, or a tool code like 4.1.2): 2.7

Every run always writes a numbered report ([0] RUN METADATA ... [10]
SUMMARY & FILES -- [5] only appears with --ml-relax), a structure
-validation checklist both before and after, and a before/after symmetry
comparison table. A full text report (--save-report), MACE pre-relaxation
(--ml-relax), and an interactive 3D view (--view) are all off by default.
Unlike earlier versions of this tool, no report file is EVER written
unless you ask for it with --save-report.
EOF
pause

echo "=================================================================="
echo " output/primitive-fcc/  --  the classic case: 4 atoms -> 1"
echo "=================================================================="
cat <<'EOF'
fcc_ni.fdf is FCC nickel's conventional cubic cell -- 4 atoms, because the
conventional cell is chosen for human readability, not minimal size. The
TRUE primitive cell has just 1 atom -- watch the reduction, and notice the
primitive cell's own lattice vectors: 60-degree angles, NOT looking cubic
at all, even though the crystal SYSTEM is still reported as cubic (crystal
system comes from the space group's full symmetry, not from how any one
particular cell choice happens to look):
EOF
mkdir -p "$OUT/primitive-fcc"
cp fcc_ni.fdf "$OUT/primitive-fcc/"
echo
echo "\$ stb-unitcell -f fcc_ni.fdf --mode primitive -o prim.fdf --no-intro"
(cd "$OUT/primitive-fcc" && stb-unitcell -f fcc_ni.fdf --mode primitive -o prim.fdf --no-intro > console.log)
grep -E "Output formula|Output atoms|Reduction factor" "$OUT/primitive-fcc/console.log"
echo
echo "The primitive cell's own lattice vectors (a rhombohedral 60-degree cell):"
python3 - "$OUT/primitive-fcc/prim.fdf" <<'PYEOF'
import sys
from stb.core import structure_io
s = structure_io.to_pymatgen(structure_io.read_fdf(sys.argv[1]))
a, b, c, al, be, ga = s.lattice.parameters
print(f"  a=b=c = {a:.4f} Ang   alpha=beta=gamma = {al:.1f} deg")
PYEOF
pause

echo "=================================================================="
echo " output/primitive-rocksalt/  --  BOTH species' ratio is preserved"
echo "=================================================================="
cat <<'EOF'
nacl_rocksalt.fdf has TWO species (Na, Cl) in its conventional 8-atom
cell. Reduction to the primitive cell keeps the exact 1:1 Na:Cl ratio --
watch it go from 4 Na + 4 Cl straight to 1 Na + 1 Cl, same physical
crystal, same space group Fm-3m as fcc_ni.fdf above (same number,
different chemistry -- a real coincidence worth noticing, not a bug):
EOF
mkdir -p "$OUT/primitive-rocksalt"
cp nacl_rocksalt.fdf "$OUT/primitive-rocksalt/"
echo
echo "\$ stb-unitcell -f nacl_rocksalt.fdf --mode primitive -o prim.fdf --no-intro"
(cd "$OUT/primitive-rocksalt" && stb-unitcell -f nacl_rocksalt.fdf --mode primitive -o prim.fdf --no-intro > console.log)
grep -E "Space group    :|Output formula|Output atoms|Reduction factor" "$OUT/primitive-rocksalt/console.log" | head -4
pause

echo "=================================================================="
echo " output/conventional-noop/  --  when there's nothing to reduce"
echo "=================================================================="
cat <<'EOF'
fcc_ni.fdf is ALREADY the conventional cell -- asking for --mode
conventional on it is a legitimate no-op, clearly flagged rather than
silently doing nothing:
EOF
mkdir -p "$OUT/conventional-noop"
cp fcc_ni.fdf "$OUT/conventional-noop/"
echo
echo "\$ stb-unitcell -f fcc_ni.fdf --mode conventional -o conv.fdf --no-intro"
(cd "$OUT/conventional-noop" && stb-unitcell -f fcc_ni.fdf --mode conventional -o conv.fdf --no-intro > console.log)
grep -E "Output atoms|NOTE" "$OUT/conventional-noop/console.log"
pause

echo "=================================================================="
echo " output/refined-noise/  --  cleaning up numerical noise"
echo "=================================================================="
cat <<'EOF'
nacl_noisy.fdf is the same NaCl crystal, but with a tiny (~2e-5 Ang)
asymmetric perturbation on every position -- the kind of residual noise a
real DFT relaxation or a hand-typed/CIF structure often has. --mode
refined keeps all 8 atoms (same cell size) but snaps every position back
to its EXACT symmetry-consistent value:
EOF
mkdir -p "$OUT/refined-noise"
cp nacl_noisy.fdf "$OUT/refined-noise/"
echo
echo "\$ stb-unitcell -f nacl_noisy.fdf --mode refined -o refined.fdf --no-intro"
(cd "$OUT/refined-noise" && stb-unitcell -f nacl_noisy.fdf --mode refined -o refined.fdf --no-intro > console.log)
grep -E "Space group    :|Output atoms" "$OUT/refined-noise/console.log" | head -2
echo
echo "Input positions (noisy):"
sed -n '/AtomicCoordinatesAndAtomicSpecies/,/endblock/p' "$OUT/refined-noise/nacl_noisy.fdf" | head -3
echo "Refined positions (exact):"
sed -n '/AtomicCoordinatesAndAtomicSpecies/,/endblock/p' "$OUT/refined-noise/refined.fdf" | head -3
pause

echo "=================================================================="
echo " output/symprec-sensitivity/  --  a real, dramatic consequence"
echo "=================================================================="
cat <<'EOF'
--symprec is the tolerance spglib uses to decide whether two positions
count as symmetry-equivalent. Run the SAME noisy structure through TWO
different tolerances: the default (1e-3, looser than the ~2e-5 noise) sees
straight through it to the true Fm-3m symmetry; a tolerance TIGHTER than
the noise itself (1e-8) sees only genuine asymmetry -- P1, no symmetry at
all -- and correctly refuses to reduce the cell (there is nothing it can
prove is redundant at that tolerance):
EOF
mkdir -p "$OUT/symprec-sensitivity"
cp nacl_noisy.fdf "$OUT/symprec-sensitivity/"
echo
echo "\$ stb-unitcell -f nacl_noisy.fdf --mode primitive -o loose.fdf --no-intro                   # default symprec"
(cd "$OUT/symprec-sensitivity" && stb-unitcell -f nacl_noisy.fdf --mode primitive -o loose.fdf --no-intro > loose.log)
echo "\$ stb-unitcell -f nacl_noisy.fdf --mode primitive --symprec 1e-8 -o tight.fdf --no-intro     # tolerance tighter than the noise"
(cd "$OUT/symprec-sensitivity" && stb-unitcell -f nacl_noisy.fdf --mode primitive --symprec 1e-8 -o tight.fdf --no-intro > tight.log)
echo
echo "Default symprec (1e-3):"
grep -E "^Space group    :|Output atoms" "$OUT/symprec-sensitivity/loose.log" | head -2
echo "Tight symprec (1e-8):"
grep -E "^Space group    :|Output atoms" "$OUT/symprec-sensitivity/tight.log" | head -2
echo
echo "Every section of the report agrees on which tolerance was used -- a"
echo "real inconsistency (found and fixed while building this example): the"
echo "validation sections used to silently default to 1e-3 regardless of"
echo "--symprec, so a tight-tolerance run would misleadingly show 'Fm-3m'"
echo "in one section and 'P1' in another for the exact same run."
pause

echo "=================================================================="
echo " output/vacuum-warning/  --  2D materials need a spot-check"
echo "=================================================================="
cat <<'EOF'
graphene.fdf is a 2D monolayer (vacuum along c). stb-unitcell operates on
the literal 3D periodic cell as given -- it has no special vacuum-axis
handling like stb-slab/stb-supercell do, so it flags this explicitly. The
vacuum thickness itself is preserved, but watch the IN-PLANE lattice
vectors come back re-expressed differently from the input (still the
exact same physical lattice -- spglib is free to pick any symmetry
-equivalent cell choice/origin, the same documented caveat noted in
core/symmetry.py elsewhere in this suite):
EOF
mkdir -p "$OUT/vacuum-warning"
cp graphene.fdf "$OUT/vacuum-warning/"
echo
echo "\$ stb-unitcell -f graphene.fdf --mode primitive -o prim.fdf --no-intro"
(cd "$OUT/vacuum-warning" && stb-unitcell -f graphene.fdf --mode primitive -o prim.fdf --no-intro > console.log)
grep -E "WARNING|Output atoms|NOTE" "$OUT/vacuum-warning/console.log"
echo
echo "Input in-plane lattice vectors:"
sed -n '/LatticeVectors/,/endblock LatticeVectors/p' "$OUT/vacuum-warning/graphene.fdf" | head -3
echo "Output in-plane lattice vectors (vacuum axis c = 20 Ang, unchanged):"
sed -n '/LatticeVectors/,/endblock LatticeVectors/p' "$OUT/vacuum-warning/prim.fdf" | head -4
pause

echo "=================================================================="
echo " output/ml-relax/  --  MACE pre-relaxation of the reduced cell"
echo "=================================================================="
if ! python3 -c "import mace" 2>/dev/null; then
    echo "Skipping -- needs the optional 'ml' extra: pip install stb_suite[ml]"
    echo "(everything else in this script works fine without it)."
else
cat <<'EOF'
--ml-relax (+ --ml-relax-cell) pre-relaxes the reduced cell with MACE
before writing it out -- useful e.g. after switching to the (cheaper)
primitive cell, to get a consistent, relaxed starting geometry for it:
EOF
    mkdir -p "$OUT/ml-relax"
    cp fcc_ni.fdf "$OUT/ml-relax/"
    echo
    echo "\$ stb-unitcell -f fcc_ni.fdf --ml-relax --ml-relax-cell -o relaxed.fdf --no-intro"
    (cd "$OUT/ml-relax" && stb-unitcell -f fcc_ni.fdf --ml-relax --ml-relax-cell \
        -o relaxed.fdf --no-intro > console.log 2>&1)
    echo "Full MACE simulation detail, straight from the report:"
    strip_ansi "$OUT/ml-relax/console.log" | awk '/\[5\] ML PRE-RELAXATION/{flag=1} /\[6\] STRUCTURE VALIDATION/{flag=0} flag'
    echo
    echo "Provenance header now also records the MACE pre-relaxation:"
    head -5 "$OUT/ml-relax/relaxed.fdf"
fi
pause

echo "=================================================================="
echo " output/full-report/  --  --save-report, validation, references.bib"
echo "=================================================================="
cat <<'EOF'
The full numbered report (also written to stb_unitcell_report.txt with
--save-report) includes the structure-validation checklist for the input
and the reduced cell, and a references.bib with SIESTA:
EOF
mkdir -p "$OUT/full-report"
cp fcc_ni.fdf "$OUT/full-report/"
echo "\$ stb-unitcell -f fcc_ni.fdf --save-report --no-intro"
(cd "$OUT/full-report" && stb-unitcell -f fcc_ni.fdf --save-report --no-intro > console.log)
echo
echo "Report sections written to stb_unitcell_report.txt:"
grep -E "^\[[0-9]+\]" "$OUT/full-report/stb_unitcell_report.txt"
echo
echo "Validation checklist (one row shown):"
grep -m1 "Atom proximity" "$OUT/full-report/console.log"
echo
echo "references.bib -- SIESTA:"
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
echo "\$ printf '2.7\\ndoes_not_exist.fdf\\nfcc_ni.fdf\\n\\n\\n\\n\\n\\n\\n0\\n' | stb-suite"
(cd "$TMP" && printf '2.7\ndoes_not_exist.fdf\nfcc_ni.fdf\n\n\n\n\n\n\n0\n' | stb-suite > session.log 2>&1) || true
if grep -q "Output atoms       : 1" "$TMP/session.log"; then
    echo "Confirmed: the interactive menu built and launched the exact same"
    echo "underlying stb-unitcell command as the CLI walkthrough above (1"
    echo "atom out, same as output/primitive-fcc/)."
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
  primitive-fcc/       primitive-rocksalt/    conventional-noop/
  refined-noise/       symprec-sensitivity/   vacuum-warning/
  ml-relax/            full-report/

Each has references.bib (SIESTA, always); full-report/ additionally has
stb_unitcell_report.txt; ml-relax/'s references.bib also cites the MACE
architecture/foundation-model papers.

Recap of what this walkthrough covered:
  - primitive vs. conventional vs. refined, and why the primitive cell is
    the cheap choice for a real DFT calculation
  - a real reduction (FCC Ni, 4 -> 1 atom) and its rhombohedral primitive
    lattice vectors, despite the crystal system staying cubic
  - a 2-species compound (NaCl) keeping its exact stoichiometric ratio
    through the same reduction
  - the --mode conventional/primitive no-op case, clearly flagged instead
    of silently doing nothing
  - --mode refined snapping noisy positions back to their exact symmetry
    -consistent values
  - --symprec's real consequence: a tolerance tighter than a structure's
    own position noise sees P1 (no symmetry) instead of the true Fm-3m --
    plus a real report inconsistency this found and fixed (validation
    sections silently ignoring --symprec)
  - the vacuum-axis caveat for 2D materials, and the "same crystal,
    different origin/lattice-vector choice" caveat this suite documents
    elsewhere (core/symmetry.py)
  - --ml-relax pre-optimizing the reduced cell
  - the structure-validation checklist, references.bib, and --save-report
  - CLI and the interactive stb-suite menu building the same command
  - stb_unitcell_report.txt is no longer written unless you ask for it

Not exercised by this script (needs a display): --view opens the input
structure and the final reduced cell in an interactive ase-gui window --
try it yourself:
  stb-unitcell -f fcc_ni.fdf --view

As a next step, try on your own:
  stb-unitcell -f your_structure.fdf --mode primitive           # smallest cell for a cheap DFT calc
  stb-unitcell -f your_relaxed_structure.fdf --mode refined      # clean up post-relaxation noise
  stb-unitcell -f fcc_ni.fdf --ml-relax --ml-relax-cell --custom-model my_finetuned.model
EOF
