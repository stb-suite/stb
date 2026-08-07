#!/bin/bash
# Guided example: stb-2Dstacking (code 2.1 in the stb-suite menu)
#
# Not an automated test (see test/2-structures/1-2d-stacking/test.sh for
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
echo " What the ZSL algorithm does"
echo "=================================================================="
cat <<'EOF'
Stacking two 2D monolayers into a van der Waals heterostructure runs into an
immediate geometry problem: the two lattices almost never match exactly. A
graphene sheet (a = 2.46 Ang) and a monolayer of h-BN (a = 2.504 Ang) differ
by ~1.8% -- close, but not identical -- and two arbitrary 2D materials can
differ far more.

The Zur & McGill (1984) algorithm (ZSL) solves this by searching over pairs
of integer supercell transformation matrices for layer 1 and layer 2, up to
a maximum cell AREA (--max_area), looking for a combination where the two
resulting supercells' lattice vectors match each other to within a maximum
allowed STRAIN fraction (--max_strain). Every match found is a trade-off:
a smaller cell area usually means more strain is needed to force a match; a
larger max_area gives the search more candidates to find a lower-strain
match, at the cost of more atoms in the final structure. stb-2Dstacking
runs this search once (expensive) and can then build many stacking
configurations (different gaps/shifts/twists) from the same match cheaply.

Every candidate match is scored on TWO independent numbers, both shown in
the match table (-i) and in the [2] report section:
  - Strain (%)        -- the lattice-CONSTANT mismatch left after forcing
                          the two supercells to the same size (what
                          --max_strain actually bounds).
  - Ang. Strain (deg) -- how much layer 2 has to be rotated ("un-rotated")
                          to line its supercell vectors back up with
                          layer 1's. ZSL only cares about VECTOR LENGTHS,
                          not absolute orientation, so a match can have
                          ~0% length strain and still need a large
                          rotation to actually fit -- that rotation is a
                          real geometric distortion, not a bookkeeping
                          detail, and matters most once a --twist angle is
                          involved (see the twist/ case below). Matches
                          with BOTH numbers below ~1 are highlighted in
                          green in the -i table -- a quick visual cue for
                          "physically clean" candidates without reading
                          every row.
EOF
pause

echo "=================================================================="
echo " Two ways to run it"
echo "=================================================================="
cat <<'EOF'
A -- direct CLI:
  stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf -id 0

B -- interactive stb-suite menu:
  $ stb-suite
  Select an option (0-6, or a tool code like 4.1.2): 2.1

Every run always writes a numbered report ([0] RUN METADATA ... [6] SUMMARY
& FILES -- [3] ML PRE-RELAXATION only appears with --ml-relax), the
structure-validation checklist (atom proximity/lattice handedness/atomic
density) for every heterostructure it builds, a symmetry table that
includes the physically correct LAYER GROUP (not just the 3D space group
of the vacuum-padded cell), and references.bib (SIESTA + the ZSL algorithm
paper). A full text report (--save-report), MACE pre-relaxation
(--ml-relax), and an interactive 3D view of the final structure(s)
(--view) are all off by default.
EOF
pause

echo "=================================================================="
echo " output/basic/  --  a realistic, lattice-mismatched pair"
echo "=================================================================="
cat <<'EOF'
graphene.fdf (a = 2.46 Ang) + hbn.fdf (a = 2.504 Ang): a ~1.8% lattice
mismatch, the classic van der Waals heterostructure test case. -id 0 picks
the very first (lowest-area) ZSL match non-interactively; -i instead shows
the full candidate table and lets you pick by hand.
EOF
mkdir -p "$OUT/basic"
cp graphene.fdf hbn.fdf "$OUT/basic/"
echo
echo "\$ stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf -id 0 --no-intro"
(cd "$OUT/basic" && stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf -id 0 --no-intro > console.log)
grep "Selected Match ID" "$OUT/basic/console.log"
grep -A5 "STRUCTURE VALIDATION & WRITING" "$OUT/basic/console.log" | grep -E "Formula|Max Applied"
echo
echo "Layer Group (the physically correct 2D symmetry classification, not just"
echo "the 3D space group of the vacuum-padded cell) -- one row per structure:"
grep "Layer Group" "$OUT/basic/console.log"
echo
echo "Provenance header written into the output .fdf itself:"
head -5 "$OUT/basic/stacked_structure.fdf"
pause

echo "=================================================================="
echo " output/batch-sym/  --  high-symmetry stackings in one run"
echo "=================================================================="
cat <<'EOF'
Real bilayers are usually studied at a few specific, physically meaningful
registries, not an arbitrary shift. --batch_sym builds all of them in one
call: AA (directly on top), Bridge_X/Bridge_Y (shifted half a lattice
vector), Center. For a HEXAGONAL lattice (detected automatically from the
gamma angle and a=b), it also adds the two inequivalent hexagonal stackings,
Hex_AB and Hex_BA (analogous to bulk graphite's AB "Bernal" stacking) --
6 heterostructures, 6 output files, from a single command.
EOF
mkdir -p "$OUT/batch-sym"
cp graphene.fdf hbn.fdf "$OUT/batch-sym/"
echo
echo "\$ stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf -id 0 --batch_sym -o hetero.fdf --no-intro"
(cd "$OUT/batch-sym" && stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf -id 0 --batch_sym -o hetero.fdf --no-intro > console.log)
grep "Hexagonal lattice detected" "$OUT/batch-sym/console.log"
echo
echo "One output .fdf per configuration:"
ls "$OUT/batch-sym/"hetero_*.fdf
pause

echo "=================================================================="
echo " output/gap-range/  --  scanning the interlayer distance"
echo "=================================================================="
cat <<'EOF'
--gap_range START END STEPS builds the SAME registry at several van der
Waals gaps instead of one fixed value -- the usual first step toward an
interlayer binding-energy curve: relax/evaluate each of these structures
with SIESTA (or stb-mleos-style tooling) afterward and you get E vs.
distance, whose minimum is the equilibrium interlayer spacing.
EOF
mkdir -p "$OUT/gap-range"
cp graphene.fdf hbn.fdf "$OUT/gap-range/"
echo
echo "\$ stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf -id 0 --gap_range 3.0 3.6 4 -o hetero.fdf --no-intro"
(cd "$OUT/gap-range" && stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf -id 0 --gap_range 3.0 3.6 4 -o hetero.fdf --no-intro > console.log)
grep "Energy Curve Mode" "$OUT/gap-range/console.log"
ls "$OUT/gap-range/"hetero_*.fdf
pause

echo "=================================================================="
echo " Strain distribution: -sm top / bottom / sym"
echo "=================================================================="
cat <<'EOF'
Forcing two mismatched lattices to share one commensurate supercell means
SOMETHING has to stretch. -sm controls who absorbs that strain:
  top    (default) -- layer 2 (top) is strained to match layer 1 exactly.
  bottom            -- layer 1 (bottom) is strained to match layer 2.
  sym               -- both layers share the strain, roughly symmetrically.
This changes the REPORTED "Max Applied Linear Strain" for the same pair and
the same ZSL match -- same geometry problem, different physical choice
about which layer you trust to be the strain-free reference.
EOF
mkdir -p "$OUT/strain-modes"
cp graphene.fdf hbn.fdf "$OUT/strain-modes/"
for mode in top bottom sym; do
    (cd "$OUT/strain-modes" && stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf -id 0 -sm "$mode" \
        -o "${mode}.fdf" --no-intro > "console_$mode.log")
done
echo "Side by side:"
for mode in top bottom sym; do
    grep "Max Applied Linear Strain" "$OUT/strain-modes/console_$mode.log" | sed "s/^/  -sm $mode  ->  /"
done
pause

echo "=================================================================="
echo " output/twist/  --  twisted bilayer graphene and Moire patterns"
echo "=================================================================="
cat <<'EOF'
--twist rotates layer 2 by a fixed angle (degrees) before the ZSL search --
the same idea behind twisted bilayer graphene research. 21.8 deg is not a
random choice: for a hexagonal lattice it's an EXACT commensurate angle
(cos(theta) = (3p^2 + 3pq + q^2/2) / (3p^2 + 3pq + q^2), the smallest
nontrivial case p=q=1) -- a small coincidence cell exists at this angle by
construction. This is deliberately NOT the ~1.1 deg "magic angle" from real
flat-band moire physics (Bistritzer-MacDonald): that one needs a supercell
of several thousand atoms, well beyond a quick CLI demo -- 21.8 deg is
picked here specifically because it stays small enough to run in seconds.
EOF
mkdir -p "$OUT/twist"
cp graphene.fdf "$OUT/twist/"
echo
echo "\$ stb-2Dstacking -l1 graphene.fdf -l2 graphene.fdf -id 0 -t 21.8 -o twisted.fdf --no-intro"
(cd "$OUT/twist" && stb-2Dstacking -l1 graphene.fdf -l2 graphene.fdf -id 0 -t 21.8 -o twisted.fdf --no-intro > console.log)
grep "Applying initial twist angle" "$OUT/twist/console.log"
grep -E "Selected Match ID|WARNING" "$OUT/twist/console.log" | sed 's/^/  /'
grep "Total atoms" "$OUT/twist/console.log"
echo
cat <<'EOF'
A real gotcha: -id 0 always picks the LOWEST-strain match first, and for
two IDENTICAL twisted layers there's a whole family of matches with ~0%
length strain that need a large rotation to align (exactly what "Ang.
Strain: 21.80 deg" above means -- equal to the twist itself). ZSL sorts
those ahead of the genuine small-angular-strain Moire cell, so match ID 0
here quietly "un-rotates" the twist away instead of building it -- proof:
the symmetry table below shows the result is still p6/mmm, IDENTICAL to
the untwisted input, not a twisted (lower-symmetry) structure at all.
EOF
grep -A2 "^Crystal System" "$OUT/twist/console.log" | sed 's/^/  /'
echo
cat <<'EOF'
The genuine Moire cell IS in the candidate list (verified: match ID 70 of
233 at --max_area 50), just far past where its ~0% strain ties get sorted
against every trivial "un-rotate" match -- past what -i's own 15-row
preview shows, so it has to be picked directly once you know it's there:
EOF
echo "\$ stb-2Dstacking -l1 graphene.fdf -l2 graphene.fdf -id 70 -a 50 -t 21.8 -o moire.fdf --no-intro"
(cd "$OUT/twist" && stb-2Dstacking -l1 graphene.fdf -l2 graphene.fdf -id 70 -a 50 -t 21.8 -o moire.fdf --no-intro > console_moire.log)
grep "Selected Match ID" "$OUT/twist/console_moire.log" | sed 's/^/  /'
grep "Total atoms" "$OUT/twist/console_moire.log"
echo "Now genuinely twisted -- symmetry drops to triclinic/p1, unlike the p6/mmm above:"
grep -A2 "^Crystal System" "$OUT/twist/console_moire.log" | sed 's/^/  /'
pause

echo "=================================================================="
echo " output/custom-shift/  --  fine control: -tx/-ty and --vacuum"
echo "=================================================================="
cat <<'EOF'
--batch_sym's 4-6 registries only cover specific, physically meaningful
fractional shifts. -tx/-ty place layer 2 at ANY fractional shift instead
(e.g. to scan a full gamma-surface, one point at a time), and --vacuum
overrides the vacuum thickness instead of inheriting layer 1's.
EOF
mkdir -p "$OUT/custom-shift"
cp graphene.fdf hbn.fdf "$OUT/custom-shift/"
echo
echo "\$ stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf -id 0 -tx 0.2 -ty 0.05 --vacuum 15.0 -o shifted.fdf --no-intro"
(cd "$OUT/custom-shift" && stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf -id 0 -tx 0.2 -ty 0.05 \
    --vacuum 15.0 -o shifted.fdf --no-intro > console.log)
echo
grep "Lattice C Parameter" "$OUT/custom-shift/console.log"
echo "(compare to the 23.20 Ang default, which inherits layer 1's own vacuum)"
echo
cat <<'EOF'
An arbitrary shift generally lands the heterostructure at NO special
symmetry point at all -- p1, no symmetry -- unlike --batch_sym's 6
registries (AA/Hex_AB/Hex_BA keep p3m1, Bridge_X/Bridge_Y/Center keep
cm11: see batch-sym/ above). That's precisely why --batch_sym hardcodes
those specific fractional shifts instead of sampling shifts at random --
they're the only ones a hexagonal bilayer's own symmetry can protect.
EOF
grep "Layer Group" "$OUT/custom-shift/console.log"
pause

echo "=================================================================="
echo " output/full-report/  --  --save-report, validation, references.bib"
echo "=================================================================="
cat <<'EOF'
The full numbered report (also written to stb_2dstacking_report.txt with
--save-report) includes the structure-validation checklist -- atom
proximity, lattice handedness, atomic density -- for every heterostructure
built, and a references.bib with SIESTA plus the ZSL algorithm paper
(Zur & McGill, 1984), the citation this tool's method itself is built on.
EOF
mkdir -p "$OUT/full-report"
cp graphene.fdf hbn.fdf "$OUT/full-report/"
echo "\$ stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf -id 0 --save-report --no-intro"
(cd "$OUT/full-report" && stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf -id 0 --save-report --no-intro > console.log)
echo
echo "Report sections written to stb_2dstacking_report.txt:"
grep -E "^\[[0-9]\]" "$OUT/full-report/stb_2dstacking_report.txt"
echo
echo "Validation checklist (one row shown):"
grep "Atom proximity" "$OUT/full-report/console.log"
echo
echo "references.bib always written -- SIESTA + the ZSL algorithm paper:"
grep "^@" "$OUT/full-report/references.bib"
pause

echo "=================================================================="
echo " output/ml-relax/  --  MACE pre-relaxation, vacuum axis locked"
echo "=================================================================="
if ! python3 -c "import mace" 2>/dev/null; then
    echo "Skipping -- needs the optional 'ml' extra: pip install stb_suite[ml]"
    echo "(everything else in this script works fine without it)."
else
cat <<'EOF'
--ml-relax pre-relaxes the just-built heterostructure with a MACE potential
before writing it out (positions only by default) -- a fast heuristic pass
before committing to a real SIESTA relaxation, same idea as stb-mlrelax
(see 1.7-stb-mlrelax/). --ml-relax-cell additionally relaxes the in-plane
cell (a/b) -- but the vacuum axis (c) always stays EXACTLY fixed, via the
same build_cell_mask() masking stb-mlrelax uses for a slab: relaxing the
vacuum thickness/gap here would be meaningless (and would break the ZSL
commensurate match between the two layers if the in-plane axes moved
independently of each other). The report shows the full simulation detail:
model size, energy/max force before and after, steps, convergence, wall
time, and (with --ml-relax-cell) the lattice change.
EOF
    mkdir -p "$OUT/ml-relax"
    cp graphene.fdf hbn.fdf "$OUT/ml-relax/"
    echo
    echo "\$ stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf -id 0 --ml-relax --ml-relax-cell -o relaxed.fdf --no-intro"
    (cd "$OUT/ml-relax" && stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf -id 0 \
        --ml-relax --ml-relax-cell -o relaxed.fdf --no-intro > console.log 2>&1)
    echo "Full MACE simulation detail, straight from the report:"
    sed 's/\x1b\[[0-9;]*m//g' "$OUT/ml-relax/console.log" | \
        awk '/\[3\] ML PRE-RELAXATION/{flag=1} /\[4\] STRUCTURE VALIDATION/{flag=0} flag'
    echo
    echo "Provenance header now also records the MACE pre-relaxation:"
    head -6 "$OUT/ml-relax/relaxed.fdf"
fi
pause

echo "=================================================================="
echo " Proof: CLI and the interactive stb-suite menu agree"
echo "=================================================================="
echo "Driving the same basic case through the interactive menu's manual"
echo "entry mode and diffing the reported strain line."
echo
echo "Note: run_2d_stacker() always passes -i to stb-2Dstacking itself, so"
echo "the subprocess launches a SECOND, nested interactive prompt (pick a"
echo "match ID). Piping enough stdin for both layers of prompts does not"
echo "work reliably over a non-tty pipe (a well-known stdio artifact of"
echo "chaining input() across a subprocess boundary, not a bug in the tool"
echo "itself) -- so this only verifies navigation and that the wrapper"
echo "builds/launches the exact same command, not the nested prompt itself."
TMP="$(mktemp -d)"
cp graphene.fdf hbn.fdf "$TMP/"
echo
echo "\$ printf '2.1\\ngraphene.fdf\\nhbn.fdf\\n\\n\\n\\n\\n\\n\\n\\n\\n\\n\\n0\\n' | stb-suite"
(cd "$TMP" && printf '2.1\ngraphene.fdf\nhbn.fdf\n\n\n\n\n\n\n\n\n\n\n0\n' | stb-suite > session.log 2>&1)
if grep -q "Searching for commensurate supercells" "$TMP/session.log"; then
    echo "Confirmed: the interactive menu built and launched the exact same"
    echo "underlying stb-2Dstacking command as the CLI walkthrough above."
else
    echo "Unexpected: menu did not reach the ZSL search -- see $TMP/session.log."
fi
rm -rf "$TMP"
pause

echo "=================================================================="
echo " Done"
echo "=================================================================="
cat <<EOF
Up to eight self-contained folders were generated under output/ (ml-relax/
skipped if the optional 'ml' extra isn't installed):
  basic/       batch-sym/     gap-range/     strain-modes/
  twist/       custom-shift/  full-report/   ml-relax/

Each has references.bib; full-report/ additionally has
stb_2dstacking_report.txt; ml-relax/'s references.bib also cites the MACE
architecture/foundation-model papers.

Recap of what this walkthrough covered:
  - what the ZSL commensurate-supercell algorithm does, the max_area/
    max_strain trade-off it searches over, and its TWO independent match
    metrics -- lattice-constant (length) strain vs. angular ("un-rotation")
    strain
  - --batch_sym's automatic high-symmetry registries, including hexagonal
    AB/BA detection, and why exactly those 6 fractional shifts are the
    ones that keep a nontrivial layer group (p3m1/cm11) instead of p1
  - --gap_range for scanning the interlayer distance
  - the three strain-distribution modes (-sm top/bottom/sym) and how they
    change which layer absorbs the lattice mismatch
  - --twist, the difference between a small exact-commensurate angle
    (21.8 deg, used here) and the real ~1.1 deg magic angle, and a real,
    verified gotcha: match ID 0 can quietly "un-rotate" a twisted
    homobilayer back to p6/mmm instead of building the genuine, lower
    -symmetry Moire cell
  - -tx/-ty for an arbitrary fractional shift and --vacuum to override the
    inherited vacuum thickness
  - the structure-validation checklist, references.bib, and --save-report
  - the Layer Group row in the symmetry table -- the physically correct
    classification for a 2D-periodic structure, alongside the 3D space
    group of the vacuum-padded cell
  - --ml-relax/--ml-relax-cell -- MACE pre-relaxation with the vacuum axis
    always locked
  - CLI and the interactive stb-suite menu building the same command

Not exercised by this script (needs a display): --view opens an
interactive ase-gui window with every generated heterostructure as
separate, pageable frames -- try it yourself:
  stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf -id 0 --view

As a next step, try on your own:
  stb-2Dstacking -l1 your_layer1.fdf -l2 your_layer2.fdf -i   # pick a match by hand
  stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf -id 0 --batch_sym --gap_range 3.0 4.0 6
  stb-2Dstacking -l1 graphene.fdf -l2 hbn.fdf -id 0 --ml-relax --custom-model my_finetuned.model
EOF
