#!/bin/bash
# Guided example: stb-inputfile (code 1.1 in the stb-suite menu)
#
# Not an automated test (see test/1-inputs/1-input_file/test.sh for that) --
# a commented walk-through: it runs real commands, one group at a time, and
# shows you the piece of output that proves what just happened. It pauses
# between sections so you can read before moving on.
#
# Every case below is written into its own folder under output/, and always
# includes pseudopotentials (-p dojo) -- so if you want to actually hand any
# of them to SIESTA, the folder is already self-contained and ready to run,
# no extra assembly needed.

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
echo " What stb-inputfile does"
echo "=================================================================="
cat <<'EOF'
stb-inputfile takes a "bare" structure (only lattice + positions + species)
and generates a complete calc.fdf, ready to run in SIESTA: a basis set, a
Monkhorst-Pack k-point grid computed from the cell's own density, and the
right SIESTA blocks for whichever kind of calculation you're preparing.

Basic usage:
  stb-inputfile <structure.fdf> -t <mode> [-p <pseudo-source>]

The four calculation modes, and what each is for:
  total_energy   single-point energy -- structure is already final
  relax          optimize atomic positions (and cell, depending on setup)
  aimd           ab initio molecular dynamics
  bands          total_energy + a k-path and PDOS block, for band structure

Any mode also accepts a "+d3" suffix (e.g. "relax+d3") to switch on the
DFT-D3 dispersion correction. "-p dojo" pulls pseudopotentials from a bank
bundled with stb_suite instead of a local folder -- used throughout this
example so every generated folder is immediately runnable.
EOF
pause

echo "=================================================================="
echo " Two ways to run it"
echo "=================================================================="
cat <<'EOF'
A -- direct CLI (used for every case below):
  stb-inputfile structure.fdf -t relax -p dojo

B -- interactive stb-suite menu, same questions as guided prompts:
  $ stb-suite
  Select an option (0-6, or a tool code like 4.1.2): 1.1

"1.1" jumps straight to this tool. It asks, one at a time: which structure
file, which mode (numbered list), an optional pseudopotential source,
whether to save a report -- then runs the exact same stb-inputfile command
underneath. Proven identical to the CLI further below.
EOF
pause

echo "=================================================================="
echo " output/relax/  --  -t relax"
echo "=================================================================="
echo "Before any real analysis, you almost always relax the structure first."
mkdir -p "$OUT/relax"
cp structure.fdf "$OUT/relax/"
echo
echo "\$ stb-inputfile structure.fdf -t relax -p dojo --no-intro"
(cd "$OUT/relax" && stb-inputfile structure.fdf -t relax -p dojo --no-intro > /dev/null)
echo "What makes 'relax' different: it turns on structural relaxation (CG):"
grep "MD.TypeOfRun" "$OUT/relax/calc.fdf"
echo
echo "Folder contents (self-contained, ready for SIESTA):"
ls "$OUT/relax"
pause

echo "=================================================================="
echo " output/bands/  --  -t bands"
echo "=================================================================="
echo "Computing bands needs a high-symmetry k-point path in the input --"
echo "stb-inputfile references one (kpath_bs.fdf); stb-kpath (code 1.3) is"
echo "what actually generates that file, so we run it too to complete the folder."
mkdir -p "$OUT/bands"
cp structure.fdf "$OUT/bands/"
echo
echo "\$ stb-inputfile structure.fdf -t bands -p dojo --no-intro"
(cd "$OUT/bands" && stb-inputfile structure.fdf -t bands -p dojo --no-intro > /dev/null)
echo "'bands' mode's fingerprints: the k-path %include and the PDOS block:"
grep -E "kpath_bs.fdf|ProjectedDensityOfStates" "$OUT/bands/calc.fdf"
echo
echo "\$ stb-kpath -f structure.fdf --no-intro"
(cd "$OUT/bands" && stb-kpath -f structure.fdf --no-intro > /dev/null)
echo
echo "Folder contents (now genuinely complete -- kpath_bs.fdf included):"
ls "$OUT/bands"
pause

echo "=================================================================="
echo " Proof: CLI and the interactive stb-suite menu agree"
echo "=================================================================="
echo "Driving the same 'relax' prompts non-interactively (via a piped printf)"
echo "and diffing the result against output/relax/calc.fdf -- no separate"
echo "output folder is kept for this, it's just a check."
TMP="$(mktemp -d)"
cp structure.fdf "$TMP/"
echo
echo "\$ printf '1.1\\nstructure.fdf\\n3\\n\\nn\\n\\n0\\n' | stb-suite"
(cd "$TMP" && printf '1.1\nstructure.fdf\n3\n\nn\n\n0\n' | stb-suite > /dev/null 2>&1)
if diff -q "$OUT/relax/calc.fdf" "$TMP/calc.fdf" > /dev/null; then
    echo "Confirmed: identical to output/relax/calc.fdf."
else
    echo "Unexpected: the two files differ -- see 'diff \"$OUT/relax/calc.fdf\" \"$TMP/calc.fdf\"'."
fi
rm -rf "$TMP"
pause

echo "=================================================================="
echo " Why dimensionality matters"
echo "=================================================================="
cat <<'EOF'
A k-point grid only makes physical sense along a genuinely periodic
direction. An axis padded with a large empty gap (a vacuum layer used to
isolate a molecule, a wire, or a slab from its own periodic images) should
get a single k-point (1), not a dense grid.

stb-inputfile detects this automatically: it measures the largest empty gap
along each lattice vector and, above a threshold (10 Ang), treats that axis
as vacuum instead of periodic when computing the k-grid -- no flag needed,
it just reads the geometry. The next 4 folders run the exact same command
against structures that only differ in how many axes are periodic:

  0D  structure_0D.fdf   CH4 molecule       -- 0 periodic axes
  1D  structure_1D.fdf   carbon chain       -- 1 periodic axis
  2D  structure_2D.fdf   graphene monolayer -- 2 periodic axes
  3D  structure.fdf      bulk silicon       -- 3 periodic axes
EOF
pause

for dim in 0D 1D 2D 3D; do
    case "$dim" in
        0D) struct="structure_0D.fdf"; label="0D -- isolated molecule (CH4)" ;;
        1D) struct="structure_1D.fdf"; label="1D -- carbon chain" ;;
        2D) struct="structure_2D.fdf"; label="2D -- graphene monolayer" ;;
        3D) struct="structure.fdf";    label="3D -- bulk silicon" ;;
    esac

    echo "=================================================================="
    echo " output/$dim/  --  $label"
    echo "=================================================================="
    mkdir -p "$OUT/$dim"
    cp "$struct" "$OUT/$dim/"
    echo "\$ stb-inputfile $struct -t total_energy -p dojo --no-intro"
    (cd "$OUT/$dim" && stb-inputfile "$struct" -t total_energy -p dojo --no-intro > /dev/null)
    grep "kgrid.MonkhorstPack" "$OUT/$dim/calc.fdf"
    ls "$OUT/$dim"
    pause
done

echo "=================================================================="
echo " Side by side"
echo "=================================================================="
for dim in 0D 1D 2D 3D; do
    printf "  %-3s " "$dim"
    grep "kgrid.MonkhorstPack" "$OUT/$dim/calc.fdf"
done
cat <<'EOF'

  0D: 1 1 1     -- no periodic axis        -> single k-point everywhere
  1D: N 1 1     -- periodic along a only   -> real grid on a, 1 elsewhere
  2D: N N 1     -- periodic along a and b  -> real grid in-plane, 1 out-of-plane
  3D: N N N     -- fully periodic          -> real grid on all 3 axes
EOF

echo
echo "=================================================================="
echo " Done"
echo "=================================================================="
cat <<EOF
Six self-contained, SIESTA-ready folders were generated under output/:
  relax/   bands/   0D/   1D/   2D/   3D/

Each one already has its structure file, calc.fdf, and pseudopotentials
(plus kpath_bs.fdf for bands/) -- hand any of them to SIESTA as-is.

As a next step, try on your own:
  stb-inputfile structure.fdf -t total_energy -p dojo
  stb-inputfile structure.fdf -t aimd -p dojo
  stb-inputfile structure.fdf -t relax+d3 -p dojo   # any mode + "+d3" adds DFT-D3
EOF
