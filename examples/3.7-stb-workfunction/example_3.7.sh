#!/bin/bash
# Guided example: stb-workfunction (code 3.7 in the stb-suite menu)
#
# Not an automated test (see test/3-analysis/7-workfunction/test.sh for that
# -- it uses small synthetic grids instead, to exercise every code path:
# asymmetric slabs, a sloped/uncorrected dipole, mismatched geometry, etc.,
# without committing large binary fixtures to the repo). This walkthrough
# uses a real, converged SIESTA calculation instead, so every number below
# is genuine DFT output, not synthetic data.
#
# calc.out/siesta.VT/siesta.XV are a real, finished SIESTA calculation
# (copied from test/3-analysis/7-workfunction/) on a CrS monolayer
# (SystemLabel "siesta") -- a genuine 2D material (fetched from the
# twodmatpedia OPTIMADE database), with a large vacuum gap along c.

set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

OUT="$DIR/output"
rm -rf "$OUT"
mkdir -p "$OUT"

# Non-interactive matplotlib backend for the "Two ways to run it" proof below
# (which navigates the interactive menu) -- --view isn't exercised in this
# script (see its own section), but this keeps any accidental plot from
# blocking the walkthrough.
export MPLBACKEND=Agg

pause() {
    echo
    read -p "  [Press Enter to continue] " -r
    echo
}

echo "=================================================================="
echo " What stb-workfunction computes, and why"
echo "=================================================================="
cat <<'EOF'
stb-workfunction reads a finished SIESTA calculation's total electrostatic
potential (.VT), planar-averages it along the direction normal to the
surface, and finds the flat vacuum plateau(s) in that profile to compute
Phi = E_vacuum - E_F -- the minimum energy to remove an electron from the
material and place it at rest in vacuum just outside the surface.

This only makes sense for a structure with a real vacuum gap (a slab, a
2D monolayer, a wire) -- a genuinely 3D-periodic bulk crystal has no
vacuum anywhere in its cell, so E_vacuum simply isn't defined there. See
the README for the full theory: why the planar average "flattens" in
vacuum, what an asymmetric slab's two distinct work functions and their
difference (dV) mean physically, and why a systematic slope across an
otherwise-flat vacuum region is the classic symptom of a missing
SlabDipoleCorrection.

Every run prints a numbered [0]...[7] report, the same style every newer
tool in this suite uses. --save-report additionally persists that report
to stb_workfunction_report.txt -- off by default. --save-gnuplot writes
workfunction_data.dat + an annotated workfunction.gplot -- also off by
default now (this tool used to write both files UNCONDITIONALLY on every
run; that's no longer the case). --view shows a matplotlib preview,
likewise off by default (replacing the old --no-plot, which showed a
plot by default and only let you turn it off).
EOF
pause

echo "=================================================================="
echo " output/basic/  --  a real, converged DFT calculation (CrS monolayer)"
echo "=================================================================="
cat <<'EOF'
Watch [1] read the real Fermi energy straight out of calc.out, [2] confirm
the vacuum axis (c) was auto-detected from siesta.XV, [3] find the one
genuinely flat vacuum plateau, and [4] report the real work function --
note the vacuum size (~3.4 Ang, the genuinely FLAT part of the profile)
is much narrower than the real geometric vacuum gap (~20 Ang) -- the
README explains exactly why these are different, both-correct numbers:
EOF
mkdir -p "$OUT/basic"
cp calc.out siesta.VT siesta.XV "$OUT/basic/"
echo
echo "\$ stb-workfunction --label siesta --file calc.out --no-intro"
(cd "$OUT/basic" && stb-workfunction --label siesta --file calc.out \
    --no-intro > console.log 2>&1)
awk '/\[1\] FERMI/{flag=1} /\[5\] OUTPUT/{flag=0} flag' "$OUT/basic/console.log"
pause

echo "=================================================================="
echo " output/save-gnuplot/  --  --save-gnuplot (off by default)"
echo "=================================================================="
cat <<'EOF'
--save-gnuplot writes workfunction_data.dat (the planar-averaged
potential profile, with the detected Fermi/vacuum levels noted in its own
header) and an annotated workfunction.gplot -- Fermi level and each vacuum
level drawn as dashed reference lines, with the work function itself
marked between them, ready to render with a real gnuplot install:
EOF
mkdir -p "$OUT/save-gnuplot"
cp calc.out siesta.VT siesta.XV "$OUT/save-gnuplot/"
echo
echo "\$ stb-workfunction --label siesta --file calc.out --save-gnuplot --no-intro"
(cd "$OUT/save-gnuplot" && stb-workfunction --label siesta --file calc.out \
    --save-gnuplot --no-intro > console.log 2>&1)
echo "workfunction_data.dat header:"
head -5 "$OUT/save-gnuplot/workfunction_data.dat"
echo
echo "workfunction.gplot (annotated with the real Fermi/vacuum/work-function values):"
grep "^set label\|^plot" "$OUT/save-gnuplot/workfunction.gplot"
pause

echo "=================================================================="
echo " output/full-report/  --  --save-report (off by default)"
echo "=================================================================="
cat <<'EOF'
Every run always prints the numbered [0]...[7] report to the console.
--save-report additionally persists it to stb_workfunction_report.txt --
off by default, so a plain run only ever writes references.bib, no text
report file (and no data/gnuplot files either, without --save-gnuplot):
EOF
mkdir -p "$OUT/full-report"
cp calc.out siesta.VT siesta.XV "$OUT/full-report/"
echo
echo "\$ stb-workfunction --label siesta --file calc.out --no-intro   # default: no report, no data/gnuplot"
(cd "$OUT/full-report" && stb-workfunction --label siesta --file calc.out \
    --no-intro > console_default.log 2>&1)
ls "$OUT/full-report/" | grep -v "^calc.out$\|^siesta.VT$\|^siesta.XV$\|console"
echo "(no workfunction_data.dat/.gplot, no stb_workfunction_report.txt -- only references.bib)"
echo
echo "\$ stb-workfunction --label siesta --file calc.out --save-report --no-intro"
(cd "$OUT/full-report" && stb-workfunction --label siesta --file calc.out \
    --save-report --no-intro > console_saved.log 2>&1)
echo "Report sections written to stb_workfunction_report.txt:"
grep -E "^\[[0-9]+\]" "$OUT/full-report/stb_workfunction_report.txt"
echo
echo "references.bib -- SIESTA (every stb-workfunction run analyzes a finished SIESTA calculation):"
grep "^@" "$OUT/full-report/references.bib"
pause

echo "=================================================================="
echo " Two ways to run it"
echo "=================================================================="
cat <<'EOF'
A -- direct CLI:
  stb-workfunction --label siesta --file calc.out

B -- interactive stb-suite menu:
  $ stb-suite
  Select an option (0-6, or a tool code like 4.1.2): 3.7

Both paths call the exact same underlying tool -- proven directly below.
EOF
TMP="$(mktemp -d)"
cp calc.out siesta.VT siesta.XV "$TMP/"
echo
echo "\$ printf '3.7\\nsiesta\\n\\ncalc.out\\n\\nn\\nn\\nn\\n' | stb-suite     # default grid, --file calc.out, auto axis, no save-report/gnuplot/view"
(cd "$TMP" && printf '3.7\nsiesta\n\ncalc.out\n\nn\nn\nn\n' | stb-suite > session.log 2>&1) || true
CLI_LINE=$(grep "Work Function |" "$OUT/basic/console.log" | head -1)
MENU_LINE=$(grep "Work Function |" "$TMP/session.log" | head -1)
if [ "$CLI_LINE" = "$MENU_LINE" ]; then
    echo "Confirmed: identical work-function line from the CLI and the interactive menu."
    echo "  $CLI_LINE"
else
    echo "Unexpected: results differ -- see $TMP/session.log."
fi
rm -rf "$TMP"
pause

echo "=================================================================="
echo " --view (needs a display)"
echo "=================================================================="
cat <<'EOF'
Not exercised by this script (needs a display): --view shows an
interactive matplotlib preview of the planar-averaged potential profile
(Fermi level and vacuum level(s) marked) before exiting:

  stb-workfunction --label siesta --file calc.out --view
EOF
pause

echo "=================================================================="
echo " Done"
echo "=================================================================="
cat <<EOF
Three self-contained folders were generated under output/:
  basic/          save-gnuplot/   full-report/

Each has references.bib; save-gnuplot/ additionally has
workfunction_data.dat/workfunction.gplot, and full-report/ has
stb_workfunction_report.txt (only from its --save-report run).

Recap of what this walkthrough covered:
  - the work function Phi = E_vacuum - E_F, and why it needs a real
    vacuum gap (no meaning for a 3D-periodic bulk crystal)
  - the planar average, why it flattens in vacuum, and why the reported
    "vacuum size" (the flat plateau) is narrower than the real geometric
    vacuum gap
  - --save-gnuplot's annotated gnuplot script vs. the default (nothing
    written unconditionally anymore)
  - the numbered [0]...[7] report, --save-report, references.bib
  - CLI and the interactive stb-suite menu building the same command

As a next step, try on your own right after a slab calculation:
  stb-workfunction --label my_slab --save-report --save-gnuplot
EOF
