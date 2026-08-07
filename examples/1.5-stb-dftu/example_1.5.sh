#!/bin/bash
# Guided example: stb-dftu (code 1.5 in the stb-suite menu)
#
# Not an automated test (see test/1-inputs/5-dftu/test.sh for that) -- a
# commented walk-through: it runs real commands, one group at a time, and
# shows you the piece of output that proves what just happened. It pauses
# between sections so you can read before moving on.

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
echo " Why plain DFT gets d/f-electron materials wrong"
echo "=================================================================="
cat <<'EOF'
Ordinary DFT (LDA/GGA) has a self-interaction error: an electron partly
interacts with itself. For most materials this barely matters, but for
LOCALIZED d or f electrons (transition-metal, lanthanide, and actinide
compounds) it artificially favors delocalization -- the classic failure is
NiO, a real antiferromagnetic insulator that plain GGA predicts as a
metal.

DFT+U (the Hubbard correction) patches this with a mean-field term that
penalizes FRACTIONAL occupation of the correlated orbital, pushing it
toward an integer occupation (0 or 1) -- closer to how a real localized
electron behaves. It's a cheap, physically-motivated patch, not a
first-principles fix.

Two parameters describe it:
  U  -- on-site Coulomb repulsion: the energy cost of putting two
        electrons in the same correlated orbital.
  J  -- Hund's-rule exchange energy (usually smaller than U).

The SIESTA %block LDAU.proj this tool generates also needs n/l -- the
correlated shell's principal and angular quantum number (l=2 for a d
shell, l=3 for f) -- auto-detected by periodic-table block (transition
metals -> nd, lanthanides -> 4f, actinides -> 5f).
EOF
pause

echo "=================================================================="
echo " Why stb-dftu never guesses U for you"
echo "=================================================================="
cat <<'EOF'
Unlike n/l (uncontroversial chemistry), U is NOT a universal constant --
it depends on the material, the functional, the pseudopotential, and the
basis set (screening, hybridization with neighboring ligands all shift
it). stb-dftu bundles one literature table (Wang, Maxisch & Ceder, Phys.
Rev. B 73, 195107 (2006) -- the Materials Project's own GGA+U oxide
calibration, and the exact reference this tool writes to references.bib)
purely as a starting-point SANITY CHECK, never a validated value -- you
always have to pass --u explicitly, or opt in with --use-reference.

For a real, first-principles U computed from YOUR OWN system (not a
generic oxide table), use stb-hubbardu / stb-hubbarduAnalysis (Workflow
menu) instead -- the Cococcioni & de Gironcoli linear-response method.
stb-dftu itself already points there in its own --help text.
EOF
pause

echo "=================================================================="
echo " Two ways to run it"
echo "=================================================================="
cat <<'EOF'
A -- direct CLI:
  stb-dftu --species Mn --u 3.9

B -- interactive stb-suite menu:
  $ stb-suite
  Select an option (0-6, or a tool code like 4.1.2): 1.5

"1.5" itself offers two sub-modes: [1] read species off a structure file
and auto-fill U from the reference table, or [2] enter species one at a
time by hand. Both end with the same save-block / save-report prompts.
EOF
pause

echo "=================================================================="
echo " output/single/  --  one species, explicit U"
echo "=================================================================="
mkdir -p "$OUT/single"
echo "\$ stb-dftu --species Mn --u 3.9 --save-report --no-intro"
(cd "$OUT/single" && stb-dftu --species Mn --u 3.9 --save-report --no-intro > console.log)
grep -E "shell=3d|Literature reference" "$OUT/single/console.log"
echo
echo "Folder contents:"
ls "$OUT/single"
pause

echo "=================================================================="
echo " output/multi/  --  two species, explicit U/J/shell"
echo "=================================================================="
echo "\$ stb-dftu --species Fe Co --u 5.3 3.32 --j 0.5 0.5 --shell 3d 3d -o ldau_block.fdf --save-report --no-intro"
mkdir -p "$OUT/multi"
(cd "$OUT/multi" && stb-dftu --species Fe Co --u 5.3 3.32 --j 0.5 0.5 --shell 3d 3d \
    -o ldau_block.fdf --save-report --no-intro > console.log)
cat "$OUT/multi/ldau_block.fdf"
echo
echo "One %block/%endblock pair wraps both stanzas -- confirmed:"
BLOCK_COUNT=$(grep -c "%block LDAU.proj" "$OUT/multi/ldau_block.fdf")
ENDBLOCK_COUNT=$(grep -c "%endblock LDAU.proj" "$OUT/multi/ldau_block.fdf")
echo "  %block occurrences: $BLOCK_COUNT, %endblock occurrences: $ENDBLOCK_COUNT"
pause

echo "=================================================================="
echo " output/from-structure/  --  species read from a structure, U auto-filled"
echo "=================================================================="
cat <<'EOF'
sc_fe_o.fdf has 3 species: Sc (a metal with NO tabulated reference), Fe (a
metal WITH one), and O (a non-metal). Watch --use-reference handle all
three differently:
EOF
mkdir -p "$OUT/from-structure"
cp sc_fe_o.fdf "$OUT/from-structure/"
echo
echo "\$ stb-dftu --fdf sc_fe_o.fdf --use-reference --save-report --no-intro"
(cd "$OUT/from-structure" && stb-dftu --fdf sc_fe_o.fdf --use-reference --save-report --no-intro > console.log)
grep -E "Species found|WARNING|^  Fe:" "$OUT/from-structure/console.log"
echo
if grep -q "^Fe   1" "$OUT/from-structure/console.log" || grep -q "Fe   1" "$OUT/from-structure/console.log"; then
    echo "Confirmed: Fe (tabulated) is IN the block."
fi
if ! grep -q "^O   1" "$OUT/from-structure/console.log"; then
    echo "Confirmed: O (non-metal) never appears in the block, no warning either."
fi
pause

echo "=================================================================="
echo " output/list-reference/ and output/suggest/  --  pure lookups"
echo "=================================================================="
echo "Neither generates a block -- both still write references.bib (just the"
echo "reference-table citation, no SIESTA entries, since nothing SIESTA-bound"
echo "was produced)."
mkdir -p "$OUT/list-reference" "$OUT/suggest"
echo
echo "\$ stb-dftu --list-reference --save-report --no-intro"
(cd "$OUT/list-reference" && stb-dftu --list-reference --save-report --no-intro > console.log)
grep "eV$" "$OUT/list-reference/console.log" | head -3
echo "..."
echo
echo "\$ stb-dftu --suggest Ni --save-report --no-intro"
(cd "$OUT/suggest" && stb-dftu --suggest Ni --save-report --no-intro > console.log)
grep "Ni:" "$OUT/suggest/console.log"
echo
for case in list-reference suggest; do
    if grep -q "@article{Soler2002," "$OUT/$case/references.bib"; then
        echo "Unexpected: $case/references.bib has a SIESTA entry it shouldn't."
    else
        echo "Confirmed: $case/references.bib has no SIESTA entry (nothing SIESTA-bound was generated)."
    fi
done
pause

echo "=================================================================="
echo " stb-dftu never guesses -- an unrecognized species is a hard error"
echo "=================================================================="
cat <<'EOF'
Si has no default correlated shell (it's not a transition metal/
lanthanide/actinide) and no tabulated reference U. Passing --u alone
isn't enough -- stb-dftu still needs to know which shell to correlate:
EOF
echo
echo "\$ stb-dftu --species Si --u 2.0 --no-intro"
set +e
SI_OUTPUT="$(stb-dftu --species Si --u 2.0 --no-intro 2>&1)"
SI_EXIT=$?
set -e
echo "$SI_OUTPUT" | grep "ERROR"
if [ "$SI_EXIT" -eq 1 ]; then
    echo "Confirmed: exit code 1, a real error, not a silent guess."
fi
pause

echo "=================================================================="
echo " Proof: CLI and the interactive stb-suite menu agree"
echo "=================================================================="
echo "Driving the same single-species (Mn) case through the interactive"
echo "menu's manual entry mode and diffing the block against output/single/."
TMP="$(mktemp -d)"
echo
echo "\$ printf '1.5\\n2\\nMn\\n3.9\\n0.0\\n\\n\\n\\n\\n\\n0\\n' | stb-suite"
(cd "$TMP" && printf '1.5\n2\nMn\n3.9\n0.0\n\n\n\n\n\n0\n' | stb-suite > session.log 2>&1)
CLI_BLOCK=$(sed -n '/^%block LDAU\.proj$/,/^%endblock LDAU\.proj$/p' "$OUT/single/console.log")
MENU_BLOCK=$(sed -n '/^%block LDAU\.proj$/,/^%endblock LDAU\.proj$/p' "$TMP/session.log")
if [ "$CLI_BLOCK" = "$MENU_BLOCK" ]; then
    echo "Confirmed: identical block from the CLI and the interactive menu."
else
    echo "Unexpected: blocks differ -- see $TMP/session.log."
fi
rm -rf "$TMP"
pause

echo "=================================================================="
echo " Done"
echo "=================================================================="
cat <<EOF
Five self-contained folders were generated under output/:
  single/   multi/   from-structure/   list-reference/   suggest/

Each has stb_dftu_report.txt and references.bib (SIESTA citations for the
3 block-generating cases; just the reference-table citation for the 2
lookup-only ones).

As a next step, try on your own:
  stb-dftu --species Ni --use-reference     # single species, auto-filled U
  stb-dftu --species Cr --u 3.7 --shell 3d

For a real U (not a generic literature table), see stb-hubbardu /
stb-hubbarduAnalysis in the Workflow menu.
EOF
