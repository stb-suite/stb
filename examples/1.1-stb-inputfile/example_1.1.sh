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

DFT-D3 (dispersion correction) and spin polarization are separate,
independent flags: -d3/--d3 and -s/--spin-polarized, each combinable with
any of the four modes above (e.g. "-t relax -d3 -s"). Spin defaults to
non-polarized if -s isn't given. "-p dojo" pulls pseudopotentials from a
bank bundled with stb_suite instead of a local folder -- used throughout
this example so every generated folder is immediately runnable.

Before writing calc.fdf, it also always runs a structure validation pass
(symmetry + malformation checks -- atoms too close, a left-handed cell, an
implausible density, a stale atom-count header) and reports whatever it
finds as [WARNING] lines -- never fatal, just a heads-up. Watch for the
"[2] STRUCTURE VALIDATION" section in every run below. --view (not used in
this unattended script, but try it yourself) opens the structure in ASE's
interactive 3D viewer right before the tool exits.
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
file, which mode (numbered list), whether to enable DFT-D3, whether to
enable spin polarization, an optional pseudopotential source, whether to
save a report -- then runs the exact same stb-inputfile command underneath.
Proven identical to the CLI further below.
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
echo " output/relax_d3/  --  -t relax -d3"
echo "=================================================================="
echo "-d3 is independent of the mode: same 'relax' template, DFT-D3 switched on."
mkdir -p "$OUT/relax_d3"
cp structure.fdf "$OUT/relax_d3/"
echo
echo "\$ stb-inputfile structure.fdf -t relax -d3 -p dojo --no-intro"
(cd "$OUT/relax_d3" && stb-inputfile structure.fdf -t relax -d3 -p dojo --no-intro > /dev/null)
echo "Only the DFTD3 line differs from output/relax/calc.fdf:"
diff <(grep -v "^#" "$OUT/relax/calc.fdf") <(grep -v "^#" "$OUT/relax_d3/calc.fdf") || true
pause

echo "=================================================================="
echo " output/relax_spin/  --  -t relax -s"
echo "=================================================================="
echo "-s is also independent of the mode -- and of -d3. Non-polarized (the"
echo "default) is what every other folder in this example uses."
mkdir -p "$OUT/relax_spin"
cp structure.fdf "$OUT/relax_spin/"
echo
echo "\$ stb-inputfile structure.fdf -t relax -s -p dojo --no-intro"
(cd "$OUT/relax_spin" && stb-inputfile structure.fdf -t relax -s -p dojo --no-intro > /dev/null)
echo "Only the Spin line differs from output/relax/calc.fdf:"
diff <(grep -v "^#" "$OUT/relax/calc.fdf") <(grep -v "^#" "$OUT/relax_spin/calc.fdf") || true
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
echo "\$ printf '1.1\\nstructure.fdf\\n2\\nn\\nn\\n\\nn\\nn\\n\\n0\\n' | stb-suite"
(cd "$TMP" && printf '1.1\nstructure.fdf\n2\nn\nn\n\nn\nn\n\n0\n' | stb-suite > /dev/null 2>&1)
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

  molecule/  structure_molecule.fdf   CH4, isolated       -- 0D, 0 periodic axes
  chain/     structure_chain.fdf      carbon chain        -- 1D, 1 periodic axis
  graphene/  structure_graphene.fdf   graphene monolayer  -- 2D, 2 periodic axes
  silicon/   structure.fdf            bulk silicon        -- 3D, 3 periodic axes
EOF
pause

for case in molecule chain graphene silicon; do
    case "$case" in
        molecule) struct="structure_molecule.fdf"; label="molecule/ -- CH4, isolated (0D)" ;;
        chain)    struct="structure_chain.fdf";    label="chain/ -- carbon chain (1D)" ;;
        graphene) struct="structure_graphene.fdf"; label="graphene/ -- graphene monolayer (2D)" ;;
        silicon)  struct="structure.fdf";          label="silicon/ -- bulk silicon (3D)" ;;
    esac

    echo "=================================================================="
    echo " output/$case  --  $label"
    echo "=================================================================="
    mkdir -p "$OUT/$case"
    cp "$struct" "$OUT/$case/"
    echo "\$ stb-inputfile $struct -t total_energy -p dojo --no-intro"
    (cd "$OUT/$case" && stb-inputfile "$struct" -t total_energy -p dojo --no-intro > /dev/null)
    grep "kgrid.MonkhorstPack" "$OUT/$case/calc.fdf"
    ls "$OUT/$case"
    pause
done

echo "=================================================================="
echo " Side by side"
echo "=================================================================="
for case in molecule chain graphene silicon; do
    printf "  %-10s " "$case"
    grep "kgrid.MonkhorstPack" "$OUT/$case/calc.fdf"
done
cat <<'EOF'

  molecule (0D): 1 1 1     -- no periodic axis        -> single k-point everywhere
  chain    (1D): N 1 1     -- periodic along a only   -> real grid on a, 1 elsewhere
  graphene (2D): N N 1     -- periodic along a and b  -> real grid in-plane, 1 out-of-plane
  silicon  (3D): N N N     -- fully periodic          -> real grid on all 3 axes
EOF

echo
echo "=================================================================="
echo " Done"
echo "=================================================================="
cat <<EOF
Eight self-contained, SIESTA-ready folders were generated under output/:
  relax/   relax_d3/   relax_spin/   bands/
  molecule/   chain/   graphene/   silicon/

Each one already has its structure file, calc.fdf, and pseudopotentials
(plus kpath_bs.fdf for bands/) -- hand any of them to SIESTA as-is.

As a next step, try on your own:
  stb-inputfile structure.fdf -t total_energy -p dojo
  stb-inputfile structure.fdf -t aimd -p dojo
  stb-inputfile structure.fdf -t bands -d3 -s -p dojo   # -d3 and -s combine with any mode
EOF
