#!/bin/bash
# Guided example: stb-stm (code 3.11 in the stb-suite menu)
#
# Not an automated test (see test/3-analysis/11-stm/test.sh for that, which
# uses its own separate graphene fixture with a conventional, centered
# vacuum axis) -- a commented walk-through: it runs real commands, one
# group at a time, into its own output/<case>/ folder, and shows you the
# piece of output that proves what just happened. Pauses between sections
# so you can read before moving on. Safe to re-run any time -- it always
# starts by wiping its own output/.
#
# siesta.LDOS/siesta.XV/calc.fdf/structure.fdf are a REAL SIESTA calculation
# of a CrS monolayer (fetched via stb-fetch from the twodmatpedia OPTIMADE
# database, id 2dm-2617 -- the same structure already used in
# examples/3.7-stb-workfunction/ and examples/3.9-stb-xrd/), with
# %block LocalDensityOfStates EF -3.50 0.00 eV (occupied states only --
# see the README for what that means and why the fdf syntax must be on one
# line). This exact fixture is what caught a real vacuum-side bug while
# building this walkthrough -- see below and the README's own theory
# section.

set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

# --view calls plt.show() -- MPLBACKEND=Agg makes that a no-op instead of
# blocking on a GUI window, same convention test.sh itself already uses.
export MPLBACKEND=Agg

OUT="$DIR/output"
rm -rf "$OUT"
mkdir -p "$OUT"

pause() {
    echo
    read -p "  [Press Enter to continue] " -r
    echo
}

echo "=================================================================="
echo " What stb-stm computes, and the fdf setting it needs"
echo "=================================================================="
cat <<'EOF'
The Tersoff-Hamann approximation: for a structureless, s-wave STM tip,
tunneling current is PROPORTIONAL to the sample's own local density of
states (LDOS) at the tip position -- so an STM image can be simulated
directly from a finished DFT calculation, no separate tip needed.

SIESTA writes that LDOS as a real-space grid via

    %block LocalDensityOfStates
    EF -3.50 0.00  eV
    %endblock LocalDensityOfStates

-- "EF" as the leading token means the two energies are measured RELATIVE
TO the Fermi level: -3.50 to 0.00 integrates the density of OCCUPIED
states only (3.50 eV below E_F up to E_F itself), the STM analog of a
negative sample bias. This is exactly the block used to generate this
folder's own siesta.LDOS, for a real CrS monolayer (the same structure
already used in examples/3.7-stb-workfunction/ and examples/3.9-stb-xrd/).
Verified live while building this example: SIESTA's parser needs
"EF <Emin> <Emax> <unit>" all on ONE line -- putting the unit on its own
line raises a hard parse error and aborts the run.

Without this block, SIESTA never writes a .LDOS file at all -- stb-stm
refuses to run (a clear error) rather than doing nothing silently. The
bias window itself is baked into the file; stb-stm cannot change it after
the fact, only image whatever window you (or whoever ran the calculation)
chose.

Every run prints a numbered [0]...[5] report. --save-report persists it;
--save-gnuplot writes stm_<mode>.dat/.gplot (off by default -- this tool
used to write both unconditionally); --view shows a matplotlib preview
(this tool previously had none at all).
EOF
pause

echo "=================================================================="
echo " A real bug this exact fixture caught: which atom is 'the surface'?"
echo "=================================================================="
cat <<'EOF'
This CrS monolayer's 4 atoms sit at fractional z = 0, 0, 0.066, 0.934 (2 Cr
at the cell's stored origin, 2 S buckled slightly to either side). A naive
"topmost atom = max(z)" picks the z=0.934 S atom and searches UPWARD from
there -- straight into a tiny ~7% sliver before wrapping back to the Cr
atoms, instead of the genuine ~20 Ang vacuum gap that actually separates
this thin monolayer's two faces (the REAL gap is the ~87% span between the
z=0.066 and z=0.934 atoms, not the ~7% wraparound remainder past the naive
max). Fixed in core/kspace.py::find_surface_reference, which reuses the
same gap-finding logic that already decides an axis IS vacuum-padded to
also find WHERE the vacuum actually starts.

A second subtlety this same structure exposed: a periodic stacking axis is
topologically a RING, so any compact atomic region surrounded by vacuum
has TWO faces exposed to that SAME gap (a "top" and a "bottom"), not one.
Searching the full ~20 Ang gap gave wildly unphysical "corrugation" values
(18-20 Ang -- real STM corrugation is sub-Angstrom to a few Ang) at several
--iso values tried live, because the search eventually crossed into the
FAR face's own LDOS tail. Capping the search at HALF the identified gap
(matching this tool's own pre-existing documented "only images one face"
limitation) fixed it -- watch [1] INPUT DATA below report a sensible
~11.7 Ang search window, not the full ~20 Ang gap nor the broken ~1.5 Ang
sliver:
EOF
mkdir -p "$OUT/constant-current"
cp siesta.LDOS siesta.XV "$OUT/constant-current/"
echo
echo "\$ stb-stm --file siesta.LDOS --geometry-file siesta.XV --save-gnuplot --no-intro"
(cd "$OUT/constant-current" && stb-stm --file siesta.LDOS --geometry-file siesta.XV \
    --save-gnuplot --no-intro > console.log 2>&1)
awk '/\[1\] INPUT DATA/{flag=1} /\[2\]/{flag=0} flag' "$OUT/constant-current/console.log"
pause

echo "=================================================================="
echo " output/constant-current/  --  the real STM feedback-loop image"
echo "=================================================================="
cat <<'EOF'
--mode current (the default) mimics a real STM's feedback loop: at every
(x, y), the tip retreats from far above the surface until the LDOS first
reaches --iso (default 0.001 e/Bohr^3) -- searched OUTSIDE-IN so the first
crossing found is where a real retreating tip would actually stop. Watch
[2] STM IMAGE report full coverage and a sub-2-Angstrom corrugation -- a
physically sensible height variation for a real 2D monolayer, only visible
once the surface-reference bug above was fixed:
EOF
awk '/\[2\] STM IMAGE/{flag=1} /\[3\]/{flag=0} flag' "$OUT/constant-current/console.log"
pause

echo "=================================================================="
echo " output/constant-height/  --  the simpler (and riskier) alternative"
echo "=================================================================="
cat <<'EOF'
--mode height just reads the LDOS grid at ONE fixed z everywhere (default
3.0 Ang above the topmost atom) -- much cheaper, but not what a real
constant-current STM does: a flat height that's safe over one atom could,
on a more corrugated surface, mean the tip has crashed into a neighboring
protrusion. [2] STM IMAGE reports the actual grid height used (nearest
available grid point to the requested one) and the plain LDOS statistics
at that height -- no "corrugation" concept here, since height itself is
fixed by construction.

Watch what happens at the DEFAULT 3.0 Ang here specifically: this occupied
-states LDOS has already decayed to essentially zero by that height (the
constant-current run above found its own contour around ~1-2 Ang) -- a
flat, uninformative image, silently. This is exactly the failure mode
--mode current's feedback loop avoids by construction: it always finds
SOME height with real signal (or clearly reports when it can't, see
iso-sensitivity/ next), while --mode height trusts YOU to have picked a
sensible z. A closer z (1.5 Ang) recovers real contrast on the same grid:
EOF
mkdir -p "$OUT/constant-height"
cp siesta.LDOS siesta.XV "$OUT/constant-height/"
echo
echo "\$ stb-stm --file siesta.LDOS --geometry-file siesta.XV --mode height --save-gnuplot --no-intro   # default z=3.0"
(cd "$OUT/constant-height" && stb-stm --file siesta.LDOS --geometry-file siesta.XV \
    --mode height --save-gnuplot --no-intro > console_z3.log 2>&1)
awk '/\[2\] STM IMAGE/{flag=1} /\[3\]/{flag=0} flag' "$OUT/constant-height/console_z3.log"
echo
echo "\$ stb-stm --file siesta.LDOS --geometry-file siesta.XV --mode height --z 1.5 --no-intro"
(cd "$OUT/constant-height" && stb-stm --file siesta.LDOS --geometry-file siesta.XV \
    --mode height --z 1.5 --no-intro > console_z1.5.log 2>&1)
awk '/\[2\] STM IMAGE/{flag=1} /\[3\]/{flag=0} flag' "$OUT/constant-height/console_z1.5.log"
pause

echo "=================================================================="
echo " output/iso-sensitivity/  --  --iso is relative, not an absolute current"
echo "=================================================================="
cat <<'EOF'
--iso is a threshold in the .LDOS file's OWN units (e/Bohr^3) -- a
proportionality, per Tersoff-Hamann, not a calibrated tunneling current in
Amperes. The same structure, three thresholds:
  - too high (0.5): almost no point in the (correctly capped) search
    window ever reaches it -- [2] reports a genuine [WARNING], NaN written
    for those points, not a silent wrong answer.
  - the default (0.001): full coverage, sensible sub-2-Angstrom corrugation.
  - very low (0.00001): still full coverage, but the tip now stops
    farther out (deeper in the vacuum tail) -- height shifts up, the
    corrugation itself changes slightly too.
EOF
mkdir -p "$OUT/iso-sensitivity"
cp siesta.LDOS siesta.XV "$OUT/iso-sensitivity/"
for iso in 0.5 0.001 0.00001; do
    echo
    echo "\$ stb-stm --file siesta.LDOS --geometry-file siesta.XV --iso $iso --no-intro"
    (cd "$OUT/iso-sensitivity" && stb-stm --file siesta.LDOS --geometry-file siesta.XV \
        --iso "$iso" --no-intro > "console_$iso.log" 2>&1)
    grep -E "Points reaching iso|never reaching iso|Corrugation" "$OUT/iso-sensitivity/console_$iso.log"
done
pause

echo "=================================================================="
echo " output/full-report/  --  --save-report / --save-gnuplot (off by default)"
echo "=================================================================="
cat <<'EOF'
Every run always prints the numbered [0]...[5] report to the console.
--save-report additionally persists it to stb_stm_report.txt, and
--save-gnuplot writes stm_<mode>.dat/stm_<mode>.gplot -- both off by
default, so a plain run only ever writes references.bib:
EOF
mkdir -p "$OUT/full-report"
cp siesta.LDOS siesta.XV "$OUT/full-report/"
echo
echo "\$ stb-stm --file siesta.LDOS --geometry-file siesta.XV --no-intro   # default: no report, no data/gnuplot"
(cd "$OUT/full-report" && stb-stm --file siesta.LDOS --geometry-file siesta.XV \
    --no-intro > console_default.log 2>&1)
ls "$OUT/full-report/" | grep -v "^siesta\|console"
echo "(no stm_*.dat/.gplot, no stb_stm_report.txt -- only references.bib)"
echo
echo "\$ stb-stm --file siesta.LDOS --geometry-file siesta.XV --save-report --save-gnuplot --no-intro"
(cd "$OUT/full-report" && stb-stm --file siesta.LDOS --geometry-file siesta.XV \
    --save-report --save-gnuplot --no-intro > console_saved.log 2>&1)
echo "Report sections written to stb_stm_report.txt:"
grep -E "^\[[0-9]+\]" "$OUT/full-report/stb_stm_report.txt"
echo
echo "references.bib -- SIESTA + Tersoff-Hamann (the STM simulation method itself):"
grep "^@" "$OUT/full-report/references.bib"
echo
echo "The .gplot script's own filenames are bare basenames -- no --output-dir"
echo "prefix (a real, verified bug fixed while rewriting this tool):"
grep -E "set output|splot" "$OUT/full-report/stm_current.gplot"
if command -v gnuplot > /dev/null; then
    (cd "$OUT/full-report" && gnuplot stm_current.gplot)
    echo "(rendered stm_current.pdf with the real, installed gnuplot, run from inside its own folder)"
fi
pause

echo "=================================================================="
echo " Two ways to run it"
echo "=================================================================="
cat <<'EOF'
A -- direct CLI:
  stb-stm --file siesta.LDOS --geometry-file siesta.XV

B -- interactive stb-suite menu:
  $ stb-suite
  Select an option (0-6, or a tool code like 4.1.2): 3.11

Both paths call the exact same underlying tool -- proven directly below.
The menu defaults --iso/--z the same way the CLI does (press Enter
through both).
EOF
TMP="$(mktemp -d)"
cp siesta.LDOS siesta.XV "$TMP/"
echo
echo "\$ printf '3.11\\nsiesta\\n1\\n\\n\\n\\nn\\nn\\nn\\n' | stb-suite     # mode 1 = current, defaults"
(cd "$TMP" && printf '3.11\nsiesta\n1\n\n\n\nn\nn\nn\n' | stb-suite > session.log 2>&1) || true
CLI_LINE=$(grep "Corrugation" "$OUT/constant-current/console.log")
MENU_LINE=$(grep "Corrugation" "$TMP/session.log")
if [ "$CLI_LINE" = "$MENU_LINE" ]; then
    echo "Confirmed: identical corrugation result from the CLI and the interactive menu."
    echo "  $CLI_LINE"
else
    echo "Unexpected: results differ -- see $TMP/session.log."
fi
rm -rf "$TMP"
pause

echo "=================================================================="
echo " Done"
echo "=================================================================="
cat <<EOF
Four self-contained folders were generated under output/:
  constant-current/   constant-height/   iso-sensitivity/   full-report/

Each has references.bib; the --save-gnuplot runs additionally have
stm_<mode>.dat/stm_<mode>.gplot, and full-report/ has stb_stm_report.txt
(only from its --save-report run).

Recap of what this walkthrough covered:
  - the Tersoff-Hamann approximation (tunneling current proportional to
    LDOS at the tip position) and why that makes --iso a RELATIVE
    threshold, not a calibrated absolute current
  - %block LocalDensityOfStates EF -3.50 0.00 eV -- required for a .LDOS
    file to exist at all, occupied-states-only via the EF-relative
    syntax, and the real fdf parser gotcha (unit must be on the same line)
  - a real, verified bug this exact fixture caught: a naive "topmost atom"
    search picked the wrong side of this monolayer's buckled atoms,
    collapsing the search window to a nonsensical ~1.5 Ang; fixed with a
    new gap-aware core/kspace.py helper, capped at half the vacuum gap to
    avoid conflating the slab's two opposite-facing surfaces
  - --mode current (the real feedback-loop image, searched outside-in)
    vs. --mode height (cheaper, fixed-height, no feedback) -- and a real
    demonstration of --mode height silently going flat/uninformative at
    a poorly-chosen z on this exact grid, fixed by picking a closer one
  - --iso sensitivity: too high leaves points as NaN (a real, reported
    warning), the default gives full coverage and sensible corrugation,
    very low shifts the whole image farther into vacuum
  - the numbered [0]...[5] report, --save-report, --save-gnuplot,
    references.bib (now including the Tersoff-Hamann citation)
  - the fixed .gplot path bug: set output/splot always use bare
    basenames now, verified with a real gnuplot run from inside the
    output folder
  - CLI and the interactive stb-suite menu building the same command

As a next step, try on your own with a real SIESTA slab calculation
(remember %block LocalDensityOfStates in the fdf):
  stb-stm --label my_slab --save-report --save-gnuplot
  stb-stm --label my_slab --mode height --z 4.0 --view
EOF
