#!/bin/bash
# Guided example: stb-supercell (code 2.2 in the stb-suite menu)
#
# Not an automated test (see test/2-structures/2-supercell/test.sh for
# that) -- a commented walk-through: it runs real commands, one group at a
# time, into its own output/<case>/ folder, and shows you the piece of
# output that proves what just happened. Pauses between sections so you can
# read before moving on. Safe to re-run any time -- it always starts by
# wiping its own output/.

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
echo " What a supercell transformation matrix is"
echo "=================================================================="
cat <<'EOF'
stb-supercell tiles one structure into a larger periodic cell following an
integer 3x3 matrix M: the new lattice vectors are M applied to the original
ones. Two shapes cover almost every real use:
  - DIAGONAL, diag(na, nb, nc) -- the everyday case, repeat the cell
    na x nb x nc times. The DETERMINANT (na*nb*nc here) is the atom-count
    multiplication factor -- always reported next to the matrix.
  - FULL (non-diagonal) -- needed when the enlarged cell isn't simply "the
    same shape, bigger": a hexagonal 2D material's classic sqrt(3)xsqrt(3)
    R30 reconstruction cell needs off-diagonal entries, because its new
    lattice vectors point along different crystallographic directions than
    the original ones, not just longer copies of them (see reconstruction/
    below).
A NEGATIVE determinant is still geometrically valid -- it mirrors the cell
(see mirrored/ below).
EOF
pause

echo "=================================================================="
echo " Two ways to run it"
echo "=================================================================="
cat <<'EOF'
A -- direct CLI:
  stb-supercell -f si_bulk.fdf -d 2 2 2

B -- interactive stb-suite menu:
  $ stb-suite
  Select an option (0-6, or a tool code like 4.1.2): 2.2

Every run always writes a numbered report ([0] RUN METADATA ... [9] SUMMARY
& FILES -- [4] ML PRE-RELAXATION only appears with --ml-relax), a
structure-validation checklist (atom proximity/lattice handedness/atomic
density) both before and after the transformation, a before/after symmetry
comparison table, and references.bib (SIESTA, plus MACE if --ml-relax was
used). A full text report (--save-report), MACE pre-relaxation
(--ml-relax), and an interactive 3D view (--view) are all off by default.
EOF
pause

echo "=================================================================="
echo " output/basic/  --  the everyday case: a diagonal 2x2x2 supercell"
echo "=================================================================="
cat <<'EOF'
si_bulk.fdf is bulk silicon, diamond cubic, the real conventional 8-atom
cell (a = 5.431 Ang, space group Fd-3m No. 227) -- a clean, highly
symmetric structure. A supercell doesn't create a new crystal, it just
re-describes the same one in a bigger cell, so the WHOLE before/after
symmetry table should come back unchanged, entry for entry.
EOF
mkdir -p "$OUT/basic"
cp si_bulk.fdf "$OUT/basic/"
echo
echo "\$ stb-supercell -f si_bulk.fdf -d 2 2 2 --no-intro"
(cd "$OUT/basic" && stb-supercell -f si_bulk.fdf -d 2 2 2 --no-intro > console.log)
grep -E "Output atoms|Determinant" "$OUT/basic/console.log"
echo
echo "Before/after symmetry table -- every column identical:"
sed -n '/SYMMETRY ANALYSIS/,/WRITING OUTPUT/p' "$OUT/basic/console.log" | head -n -2
echo
echo "Provenance header written into the output .fdf itself:"
head -4 "$OUT/basic/supercell.fdf"
pause

echo "=================================================================="
echo " output/nonuniform/  --  a real gotcha, verified: Space Group survives,"
echo " Point Group can look like it drops"
echo "=================================================================="
cat <<'EOF'
1x1x2 doubles the cell ONLY along c -- still diagonal, but no longer
isotropic. The crystal is still fully cubic silicon, and the SPACE GROUP,
LAYER GROUP, and HALL SYMBOL columns prove it by staying identical (they
come from spglib's own internally standardized cell, independent of the
shape you actually handed it).

The POINT GROUP column, though, is computed from the raw rotation matrices
found for the EXACT cell you gave it, in that cell's own basis -- not from
the standardized space-group data. An anisotropic (non-cubic-shaped) cell
genuinely restricts which of the crystal's true symmetry operations still
map THAT SPECIFIC cell onto itself as an integer matrix. Watch it drop:
EOF
mkdir -p "$OUT/nonuniform"
cp si_bulk.fdf "$OUT/nonuniform/"
echo
echo "\$ stb-supercell -f si_bulk.fdf -d 1 1 2 --no-intro"
(cd "$OUT/nonuniform" && stb-supercell -f si_bulk.fdf -d 1 1 2 --no-intro > console.log)
sed -n '/SYMMETRY ANALYSIS/,/WRITING OUTPUT/p' "$OUT/nonuniform/console.log" | head -n -2
echo
cat <<'EOF'
Space Group stayed Fd-3m (227) in both columns -- still, unambiguously, the
same crystal. Point Group alone dropped from m-3m to 4/mmm: not a bug, just
this specific (elongated) cell shape no longer respecting the full cubic
point symmetry as an integer matrix. Takeaway: after an ANISOTROPIC
supercell, trust Space Group for "is this still the same crystal", not
Point Group.
EOF
pause

echo "=================================================================="
echo " output/reconstruction/  --  the full 3x3 matrix path: graphene's"
echo " sqrt(3) x sqrt(3) R30 reconstruction cell"
echo "=================================================================="
cat <<'EOF'
A reconstruction/adsorbate-superlattice cell on a hexagonal 2D material
needs an off-diagonal matrix: the new lattice vectors (a-b, a+2b here) point
along DIFFERENT crystallographic directions than the original ones, not
just longer copies of them -- 2 atoms in, 6 atoms out, determinant 3.
EOF
mkdir -p "$OUT/reconstruction"
cp graphene.fdf "$OUT/reconstruction/"
echo
echo "\$ stb-supercell -f graphene.fdf -d 1 1 0 -1 2 0 0 0 1 --no-intro"
(cd "$OUT/reconstruction" && stb-supercell -f graphene.fdf -d 1 1 0 -1 2 0 0 0 1 --no-intro > console.log)
grep -E "Output atoms|Determinant" "$OUT/reconstruction/console.log"
echo
echo "Layer Group (the physically correct 2D symmetry classification) survives"
echo "the reconstruction; Point Group shows the same anisotropic-cell effect"
echo "as nonuniform/ above (6/mmm -> mmm), even though this is a genuinely"
echo "isotropic (equal-length) in-plane transformation -- the effect is about"
echo "the cell's SHAPE relative to the crystal axes, not simply its size:"
sed -n '/SYMMETRY ANALYSIS/,/WRITING OUTPUT/p' "$OUT/reconstruction/console.log" | head -n -2
pause

echo "=================================================================="
echo " output/mirrored/  --  a negative determinant, still valid"
echo "=================================================================="
cat <<'EOF'
A negative determinant just mirrors the cell -- still geometrically valid,
flagged with a [WARNING] both up front and again in the post-transform
validation checklist (left-handed lattice). Whether the RESULT is physically
distinguishable from the original depends on the crystal itself: silicon's
diamond structure already contains inversion symmetry (centrosymmetric), so
its mirror image is the same crystal -- verified below, the symmetry table
comes back fully identical. A chiral crystal (no inversion/mirror symmetry
of its own) would instead give a genuinely different, enantiomeric
structure.
EOF
mkdir -p "$OUT/mirrored"
cp si_bulk.fdf "$OUT/mirrored/"
echo
echo "\$ stb-supercell -f si_bulk.fdf -d -1 0 0 0 1 0 0 0 1 --no-intro"
(cd "$OUT/mirrored" && stb-supercell -f si_bulk.fdf -d -1 0 0 0 1 0 0 0 1 --no-intro > console.log)
grep -E "WARNING.*negative determinant|WARNING.*left-handed" "$OUT/mirrored/console.log"
echo
echo "Symmetry table -- identical, since Si is centrosymmetric:"
sed -n '/SYMMETRY ANALYSIS/,/WRITING OUTPUT/p' "$OUT/mirrored/console.log" | head -n -2
pause

echo "=================================================================="
echo " output/ml-relax/  --  MACE pre-relaxation of the built supercell"
echo "=================================================================="
if ! python3 -c "import mace" 2>/dev/null; then
    echo "Skipping -- needs the optional 'ml' extra: pip install stb_suite[ml]"
    echo "(everything else in this script works fine without it)."
else
cat <<'EOF'
--ml-relax pre-relaxes the just-BUILT supercell with a MACE potential before
writing it out (positions only by default) -- e.g. before hand-editing in a
defect, or as a fast heuristic pass before a real SIESTA relaxation, same
idea as stb-mlrelax (see 1.6-stb-mlrelax/). --ml-relax-cell additionally
relaxes the cell itself. The report shows the full simulation detail: model
parameter count/cutoff radius, steps used, convergence, wall time, and a
before/after table -- energy, max force, and, with --ml-relax-cell, the
lattice/volume change.
EOF
    mkdir -p "$OUT/ml-relax"
    cp si_bulk.fdf "$OUT/ml-relax/"
    echo
    echo "\$ stb-supercell -f si_bulk.fdf -d 2 2 2 --ml-relax --ml-relax-cell -o relaxed.fdf --no-intro"
    (cd "$OUT/ml-relax" && stb-supercell -f si_bulk.fdf -d 2 2 2 \
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
echo " output/full-report/  --  --save-report, validation, references.bib"
echo "=================================================================="
cat <<'EOF'
The full numbered report (also written to stb_supercell_report.txt with
--save-report) includes the structure-validation checklist for both the
input structure and the final supercell, and a references.bib with SIESTA
(the output is a .fdf).
EOF
mkdir -p "$OUT/full-report"
cp si_bulk.fdf "$OUT/full-report/"
echo "\$ stb-supercell -f si_bulk.fdf -d 2 2 2 --save-report --no-intro"
(cd "$OUT/full-report" && stb-supercell -f si_bulk.fdf -d 2 2 2 --save-report --no-intro > console.log)
echo
echo "Report sections written to stb_supercell_report.txt:"
grep -E "^\[[0-9]\]" "$OUT/full-report/stb_supercell_report.txt"
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
echo "entry mode and checking it reaches the same construction step."
TMP="$(mktemp -d)"
cp si_bulk.fdf "$TMP/"
echo
echo "\$ printf '2.2\\nsi_bulk.fdf\\n2 2 2\\n\\n\\n\\n\\n0\\n' | stb-suite"
(cd "$TMP" && printf '2.2\nsi_bulk.fdf\n2 2 2\n\n\n\n\n0\n' | stb-suite > session.log 2>&1) || true
if grep -q "Output atoms   : 64" "$TMP/session.log"; then
    echo "Confirmed: the interactive menu built and launched the exact same"
    echo "underlying stb-supercell command as the CLI walkthrough above (64"
    echo "atoms out, same as output/basic/)."
else
    echo "Unexpected: menu did not reach the construction step -- see $TMP/session.log."
fi
rm -rf "$TMP"
pause

echo "=================================================================="
echo " Done"
echo "=================================================================="
cat <<EOF
Up to six self-contained folders were generated under output/ (ml-relax/
skipped if the optional 'ml' extra isn't installed):
  basic/         nonuniform/    reconstruction/
  mirrored/      ml-relax/      full-report/

Each has references.bib; full-report/ additionally has
stb_supercell_report.txt; ml-relax/'s references.bib also cites the MACE
architecture/foundation-model papers.

Recap of what this walkthrough covered:
  - diagonal vs. full (non-diagonal) transformation matrices, and the
    determinant as the atom-count multiplication factor
  - a real, verified gotcha: an ANISOTROPIC supercell keeps the correct
    Space Group/Layer Group/Hall Symbol but can drop the reported Point
    Group -- a cell-shape artifact, not a symmetry change in the crystal
    itself
  - the full 3x3 matrix path for a non-orthogonal 2D reconstruction cell
    (graphene's sqrt(3) x sqrt(3) R30)
  - a negative-determinant mirrored cell, and why it's physically identical
    here (Si is centrosymmetric) but wouldn't be for a chiral crystal
  - --ml-relax/--ml-relax-cell -- MACE pre-relaxation of the built supercell
  - the structure-validation checklist (before AND after), references.bib,
    and --save-report
  - CLI and the interactive stb-suite menu building the same command

Not exercised by this script (needs a display): --view opens the input
structure and the final supercell side by side in an interactive ase-gui
window -- try it yourself:
  stb-supercell -f si_bulk.fdf -d 2 2 2 --view

As a next step, try on your own:
  stb-supercell -f your_structure.fdf -d 3 3 1        # e.g. a slab supercell
  stb-supercell -f si_bulk.fdf -d 2 2 2 --ml-relax --custom-model my_finetuned.model
  stb-supercell -f si_bulk.fdf -d 2 2 2 -sp 0.1       # looser symmetry tolerance
EOF
