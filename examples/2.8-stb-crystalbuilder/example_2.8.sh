#!/bin/bash
# Guided example: stb-crystalbuilder (code 2.8 in the stb-suite menu)
#
# Not an automated test (see test/2-structures/8-crystalbuilder/test.sh for
# that) -- a commented walk-through: it runs real commands, one group at a
# time, into its own output/<case>/ folder, and shows you the piece of
# output that proves what just happened. Pauses between sections so you
# can read before moving on. Safe to re-run any time -- it always starts
# by wiping its own output/. Unlike most 2-structures examples, no input
# fixture files are needed: stb-crystalbuilder builds everything from
# --spacegroup/--site flags.

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
echo " What is a space group, and why build FROM one?"
echo "=================================================================="
cat <<'EOF'
Every 3D crystal's full symmetry is one of exactly 230 international SPACE
GROUPS -- the complete set of translations/rotations/mirrors/screw-glide
operations that map the infinite lattice onto itself. Rather than typing
out every atom's position by hand, you only need to give the minimal,
symmetrically-distinct set of atoms (their WYCKOFF POSITIONS) -- the space
group's own operations generate every other atom automatically. This is
exactly how real crystal structures are published (a paper or database
entry lists a space group + a handful of Wyckoff sites, not a giant atom
list) -- stb-crystalbuilder's --site flag mirrors that convention: one
entry per Wyckoff letter, not every atom.

stb-crystalbuilder is the INVERSE of stb-unitcell (2.7): that tool takes a
full structure and reduces it down to a smaller cell; this one takes the
minimal symmetry description and builds the full structure up.
EOF
pause

echo "=================================================================="
echo " Is it only for 3D systems? YES -- and this is the first thing to know"
echo "=================================================================="
cat <<'EOF'
--spacegroup only accepts the 230 INTERNATIONAL SPACE GROUPS -- the
complete classification for 3D-PERIODIC crystals. There is no 2D (layer
group) or 1D (rod group) equivalent here: pymatgen's Structure.from_
spacegroup (what this tool wraps) is inherently a 3D-crystal builder.

Every run's own report now says so explicitly, in [2] STRUCTURE
VALIDATION:
  Dimensionality : 3D (bulk material) -- every space group is inherently
  3D-periodic; use stb-slab/stb-nanotube/stb-2Dstacking on the output to
  cut a 2D/1D structure from it.

If you actually want a slab, monolayer, nanotube, or wire: build the 3D
bulk crystal here FIRST, then cut it down with the tool that understands
vacuum/lower dimensionality (stb-slab, stb-nanotube, stb-2Dstacking). The
LAST case in this walkthrough (bulk-graphite-to-slab/) proves this
concretely with a real material.

Two more places this is made explicit:
  - The [7] SYMMETRY ANALYSIS (BEFORE/AFTER) table never shows a "Layer
    Group" row (unlike stb-unitcell/stb-slab) -- it would always read
    "N/A (not 2D-periodic)" here, never real information, so it's left
    out entirely instead of printing a permanently-useless row.
  - The interactive stb-suite menu (code 2.8) prints this same bulk-only
    warning immediately, before even asking for a space group -- see the
    proof near the end of this script.
EOF
pause

echo "=================================================================="
echo " Two ways to run it"
echo "=================================================================="
cat <<'EOF'
A -- direct CLI:
  stb-crystalbuilder --spacegroup Fm-3m --a 3.52 --site Ni 0 0 0

B -- interactive stb-suite menu:
  $ stb-suite
  Select an option (0-6, or a tool code like 4.1.2): 2.8

Every run always writes a numbered report ([0] RUN METADATA ... [10]
SUMMARY & FILES -- [4] only with --reduce, [5] only with --ml-relax), a
structure-validation checklist before and after, full symmetry detection,
and a before/after symmetry comparison table. A full text report
(--save-report), cell reduction (--reduce), MACE pre-relaxation
(--ml-relax), and an interactive 3D view (--view) are all off by default.
EOF
pause

echo "=================================================================="
echo " output/fcc-nickel/  --  the basic case: 1 site -> 4 atoms"
echo "=================================================================="
cat <<'EOF'
FCC nickel: space group Fm-3m (No. 225), a single Wyckoff site (Ni at the
origin, Wyckoff letter 4a). Watch the symmetry operations expand this ONE
given site into 4 atoms in the conventional cubic cell -- and confirm the
symbol form (Fm-3m) and the equivalent international-number form (225)
build the exact same structure:
EOF
mkdir -p "$OUT/fcc-nickel"
echo
echo "\$ stb-crystalbuilder --spacegroup Fm-3m --a 3.52 --site Ni 0 0 0 -o fcc_symbol.fdf --no-intro"
(cd "$OUT/fcc-nickel" && stb-crystalbuilder --spacegroup Fm-3m --a 3.52 --site Ni 0 0 0 \
    -o fcc_symbol.fdf --no-intro > symbol.log)
echo "\$ stb-crystalbuilder --spacegroup 225   --a 3.52 --site Ni 0 0 0 -o fcc_int.fdf    --no-intro"
(cd "$OUT/fcc-nickel" && stb-crystalbuilder --spacegroup 225 --a 3.52 --site Ni 0 0 0 \
    -o fcc_int.fdf --no-intro > int.log)
grep -E "Output formula|Output atoms|Space group    :" "$OUT/fcc-nickel/symbol.log" | tail -3
echo
if diff -q <(grep -v '^#' "$OUT/fcc-nickel/fcc_symbol.fdf") <(grep -v '^#' "$OUT/fcc-nickel/fcc_int.fdf") > /dev/null; then
    echo "Confirmed: symbol ('Fm-3m') and number ('225') forms produce the"
    echo "identical structure (only the provenance header's own \"Requested"
    echo "space group\" line differs, since that literally records what you typed)."
fi
pause

echo "=================================================================="
echo " output/magnetite-spinel/  --  a real oxide, 2 distinct Fe sites"
echo "=================================================================="
cat <<'EOF'
Magnetite (Fe3O4), space group Fd-3m (No. 227): TWO symmetrically-distinct
iron sites (tetrahedral and octahedral coordination) plus one oxygen site
-- exactly the 3 Wyckoff positions a real crystallography table lists for
this mineral. Watch 3 --site entries expand into the full 56-atom
conventional spinel cell:
EOF
mkdir -p "$OUT/magnetite-spinel"
echo
echo "\$ stb-crystalbuilder --spacegroup Fd-3m --a 8.396 \\"
echo "    --site Fe 0.125 0.125 0.125 --site Fe 0.5 0.5 0.5 --site O 0.379 0.379 0.379 \\"
echo "    -o magnetite.fdf --no-intro"
(cd "$OUT/magnetite-spinel" && stb-crystalbuilder --spacegroup Fd-3m --a 8.396 \
    --site Fe 0.125 0.125 0.125 --site Fe 0.5 0.5 0.5 --site O 0.379 0.379 0.379 \
    -o magnetite.fdf --no-intro > console.log)
grep -E "Output formula|Output atoms|Point group|Crystal system" "$OUT/magnetite-spinel/console.log"
pause

echo "=================================================================="
echo " output/special-position-polonium/  --  you can ask for LESS symmetry than you get"
echo "=================================================================="
cat <<'EOF'
A single atom placed at a lattice's origin can accidentally sit on MORE
symmetry than the space group you requested. Ask for the tetragonal group
P4/mmm (No. 123) but give a perfectly cubic lattice (a=b=c, all angles 90)
with one atom at the origin -- alpha-polonium's real structure, one of the
only elements that crystallizes simple-cubic. The actual structure snaps
to the FULL cubic symmetry (Pm-3m, No. 221), and the tool reports the
mismatch explicitly instead of silently building something you didn't
ask for:
EOF
mkdir -p "$OUT/special-position-polonium"
echo
echo "\$ stb-crystalbuilder --spacegroup P4/mmm --a 3.52 --site Po 0 0 0 -o po.fdf --no-intro"
(cd "$OUT/special-position-polonium" && stb-crystalbuilder --spacegroup P4/mmm --a 3.52 \
    --site Po 0 0 0 -o po.fdf --no-intro > console.log)
grep -E "Requested space group :|^Space group    :|NOTE" "$OUT/special-position-polonium/console.log" | sed -n '2,4p'
pause

echo "=================================================================="
echo " output/reduce-to-primitive/  --  build AND reduce in one command"
echo "=================================================================="
cat <<'EOF'
Space-group construction naturally gives the CONVENTIONAL cell (the one a
database/paper actually lists) -- often bigger than a real DFT calculation
needs. --reduce (new in this version) folds in stb-unitcell's own
reduction directly: magnetite's 56-atom conventional cell becomes its
14-atom primitive cell, no separate stb-unitcell call needed (see
2.7-stb-unitcell/ for the full primitive/conventional/refined theory):
EOF
mkdir -p "$OUT/reduce-to-primitive"
echo
echo "\$ stb-crystalbuilder --spacegroup Fd-3m --a 8.396 \\"
echo "    --site Fe 0.125 0.125 0.125 --site Fe 0.5 0.5 0.5 --site O 0.379 0.379 0.379 \\"
echo "    --reduce primitive -o magnetite_prim.fdf --no-intro"
(cd "$OUT/reduce-to-primitive" && stb-crystalbuilder --spacegroup Fd-3m --a 8.396 \
    --site Fe 0.125 0.125 0.125 --site Fe 0.5 0.5 0.5 --site O 0.379 0.379 0.379 \
    --reduce primitive -o magnetite_prim.fdf --no-intro > console.log)
strip_ansi "$OUT/reduce-to-primitive/console.log" | grep -A5 "UNIT CELL REDUCTION"
pause

echo "=================================================================="
echo " output/out-of-range-error/  --  a fixed bug: a clear error now"
echo "=================================================================="
cat <<'EOF'
Space group numbers only go from 1 to 230. Requesting an out-of-range
number used to fail with a confusing, pymatgen-internal message ("Bad
international symbol '255'") -- caused by pymatgen's own from_spacegroup()
swallowing its clearer "must be between 1 and 230" ValueError and
falling back to a symbol lookup instead. Fixed: the range is now checked
directly, with a clear error:
EOF
mkdir -p "$OUT/out-of-range-error"
echo
echo "\$ stb-crystalbuilder --spacegroup 255 --a 3.52 --site Ni 0 0 0 --no-intro"
set +e
(cd "$OUT/out-of-range-error" && stb-crystalbuilder --spacegroup 255 --a 3.52 \
    --site Ni 0 0 0 --no-intro > console.log 2>&1)
EXIT_CODE=$?
set -e
grep "ERROR" "$OUT/out-of-range-error/console.log"
echo "(exit code: $EXIT_CODE)"
pause

echo "=================================================================="
echo " output/ml-relax/  --  MACE pre-relaxation of the built + reduced cell"
echo "=================================================================="
if ! python3 -c "import mace" 2>/dev/null; then
    echo "Skipping -- needs the optional 'ml' extra: pip install stb_suite[ml]"
    echo "(everything else in this script works fine without it)."
else
cat <<'EOF'
--ml-relax (+ --ml-relax-cell) pre-relaxes the final structure with MACE
before writing it out. Since every space-group-built structure is
inherently 3D-periodic, --ml-relax-cell here always relaxes all 3 cell
axes -- no vacuum-axis masking needed (unlike stb-unitcell, which must
also handle arbitrary user-supplied slabs/wires):
EOF
    mkdir -p "$OUT/ml-relax"
    echo
    echo "\$ stb-crystalbuilder --spacegroup Fm-3m --a 3.52 --site Ni 0 0 0 --reduce primitive \\"
    echo "    --ml-relax --ml-relax-cell -o relaxed.fdf --no-intro"
    (cd "$OUT/ml-relax" && stb-crystalbuilder --spacegroup Fm-3m --a 3.52 --site Ni 0 0 0 \
        --reduce primitive --ml-relax --ml-relax-cell -o relaxed.fdf --no-intro > console.log 2>&1)
    echo "Full MACE simulation detail, straight from the report:"
    strip_ansi "$OUT/ml-relax/console.log" | awk '/\[5\] ML PRE-RELAXATION/{flag=1} /\[6\] STRUCTURE VALIDATION/{flag=0} flag'
    echo
    echo "Provenance header now records the full history: space group, reduce mode, and MACE relax:"
    head -6 "$OUT/ml-relax/relaxed.fdf"
fi
pause

echo "=================================================================="
echo " output/full-report/  --  --save-report, validation, references.bib"
echo "=================================================================="
cat <<'EOF'
The full numbered report (also written to stb_crystalbuilder_report.txt
with --save-report) includes the structure-validation checklist and a
references.bib with SIESTA:
EOF
mkdir -p "$OUT/full-report"
echo "\$ stb-crystalbuilder --spacegroup Fm-3m --a 3.52 --site Ni 0 0 0 --save-report --no-intro"
(cd "$OUT/full-report" && stb-crystalbuilder --spacegroup Fm-3m --a 3.52 --site Ni 0 0 0 \
    --save-report --no-intro > console.log)
echo
echo "Report sections written to stb_crystalbuilder_report.txt:"
grep -E "^\[[0-9]+\]" "$OUT/full-report/stb_crystalbuilder_report.txt"
echo
echo "references.bib -- SIESTA:"
grep "^@" "$OUT/full-report/references.bib"
pause

echo "=================================================================="
echo " output/bulk-graphite-to-slab/  --  THE 3D-only answer, proven"
echo "=================================================================="
cat <<'EOF'
Real bulk graphite: space group P6_3/mmc (No. 194), hexagonal, TWO Wyckoff
carbon sites (the two sublattices of each graphene layer, AB-stacked) --
watch the built structure's own C-C bond length come out at the real,
correct value (~1.42 Ang) purely from the space group + lattice constants
you give it. stb-crystalbuilder itself has NO concept of a vacuum gap or
a monolayer -- it only ever builds the 3D bulk. To get an isolated 2D
structure out of it, hand the result to stb-slab (2.3):
EOF
mkdir -p "$OUT/bulk-graphite-to-slab"
echo
echo "\$ stb-crystalbuilder --spacegroup \"P6_3/mmc\" --a 2.464 --c 6.711 --gamma 120 \\"
echo "    --site C 0 0 0.25 --site C 0.333333333 0.666666667 0.25 \\"
echo "    -o graphite.fdf --no-intro"
(cd "$OUT/bulk-graphite-to-slab" && stb-crystalbuilder --spacegroup "P6_3/mmc" --a 2.464 \
    --c 6.711 --gamma 120 --site C 0 0 0.25 --site C 0.333333333 0.666666667 0.25 \
    -o graphite.fdf --no-intro > build.log)
echo "Bulk graphite, straight from stb-crystalbuilder:"
grep -E "Dimensionality|Atom proximity|Output formula|Output atoms|^Space group    :" "$OUT/bulk-graphite-to-slab/build.log"
echo
echo "\$ stb-slab -f graphite.fdf --hkl 0 0 1 --min-slab-size 3.0 --min-vacuum-size 15 -o graphene.fdf --no-intro"
(cd "$OUT/bulk-graphite-to-slab" && stb-slab -f graphite.fdf --hkl 0 0 1 --min-slab-size 3.0 \
    --min-vacuum-size 15 -o graphene.fdf --no-intro > slab.log)
echo
echo "The slab's before/after table: Layer Group is N/A for the 3D bulk"
echo "(exactly the limitation this whole example is about), but becomes a"
echo "REAL, detected 2D layer group for the cut slab -- direct proof that"
echo "stb-slab is what turns a 3D-only crystal-builder output into a"
echo "genuinely 2D-periodic structure:"
strip_ansi "$OUT/bulk-graphite-to-slab/slab.log" | grep -E "Layer Group|Crystal System|Space Group    \|"
pause

echo "=================================================================="
echo " Proof: CLI and the interactive stb-suite menu agree"
echo "=================================================================="
echo "Driving the same fcc-nickel case through the interactive menu's"
echo "manual entry mode and checking it reaches the same atom count."
TMP="$(mktemp -d)"
echo
echo "\$ printf '2.8\\nFm-3m\\n3.52\\n\\n\\n\\n\\n\\nNi 0 0 0\\n\\n\\n\\n\\n\\n\\n0\\n' | stb-suite"
(cd "$TMP" && printf '2.8\nFm-3m\n3.52\n\n\n\n\n\nNi 0 0 0\n\n\n\n\n\n\n0\n' | stb-suite > session.log 2>&1) || true
if grep -q "Output atoms   : 4" "$TMP/session.log"; then
    echo "Confirmed: the interactive menu built and launched the exact same"
    echo "underlying stb-crystalbuilder command as the CLI walkthrough above"
    echo "(4 atoms out, same as output/fcc-nickel/)."
else
    echo "Unexpected: menu did not reach the write step -- see $TMP/session.log."
fi
echo
echo "The menu also prints the bulk-only warning immediately, before even"
echo "asking for a space group -- the first thing a new user sees, not just"
echo "a line buried in --help:"
strip_ansi "$TMP/session.log" | grep "Bulk (3D periodic) structures only"
rm -rf "$TMP"
pause

echo "=================================================================="
echo " Done"
echo "=================================================================="
cat <<EOF
Up to eight self-contained folders were generated under output/ (ml-relax/
skipped if the optional 'ml' extra isn't installed):
  fcc-nickel/                magnetite-spinel/       special-position-polonium/
  reduce-to-primitive/       out-of-range-error/     ml-relax/
  full-report/               bulk-graphite-to-slab/

Each has references.bib (SIESTA, always); full-report/ additionally has
stb_crystalbuilder_report.txt; ml-relax/'s references.bib also cites the
MACE architecture/foundation-model papers.

Recap of what this walkthrough covered:
  - space groups and Wyckoff positions: give the minimal symmetrically
    -distinct sites, the space group's own operations build the rest
  - stb-crystalbuilder is 3D (bulk) ONLY -- every space group is
    inherently a 3D-periodic group, there is no 2D/1D analog here
  - the fix: build the 3D bulk crystal here, then cut it down with
    stb-slab/stb-nanotube/stb-2Dstacking for anything lower-dimensional
    -- proven directly with real bulk graphite -> a 2D slab
  - a real reduction (FCC Ni, symbol vs. number spacegroup forms agree)
  - a real oxide (magnetite) with 2 distinct Fe sites + 1 O site
  - special positions: asking for less symmetry than your sites actually
    have (real alpha-polonium, simple-cubic) -- and the tool's explicit
    [NOTE] about it
  - --reduce: build AND reduce (primitive/conventional/refined) in one
    command, reusing stb-unitcell's own reduction
  - the fixed "Bad international symbol" bug: out-of-range space group
    numbers now give a clear, direct error message
  - --ml-relax pre-optimizing the final structure
  - the structure-validation checklist, references.bib, and --save-report
  - CLI and the interactive stb-suite menu building the same command

Not exercised by this script (needs a display): --view opens the built
and final structures in an interactive ase-gui window -- try it yourself:
  stb-crystalbuilder --spacegroup Fm-3m --a 3.52 --site Ni 0 0 0 --view

As a next step, try on your own:
  stb-crystalbuilder --spacegroup <symbol-or-number> --a <a> --site <SYMBOL> <x> <y> <z> ...
  stb-crystalbuilder --spacegroup Fm-3m --a 3.52 --site Ni 0 0 0 --reduce primitive
  stb-crystalbuilder --spacegroup "P6_3/mmc" --a 2.464 --c 6.711 --gamma 120 \\
      --site C 0 0 0.25 --site C 0.333333 0.666667 0.25 -o graphite.fdf
  stb-slab -f graphite.fdf --hkl 0 0 1   # then cut it down to 2D
EOF
