#!/bin/bash
# Guided example: stb-bader (code 3.6 in the stb-suite menu)
#
# Not an automated test (see test/3-analysis/6-bader/test.sh for that)
# -- a commented walk-through: it runs real commands, one group at a time,
# into its own output/<case>/ folder, and shows you the piece of output
# that proves what just happened. Pauses between sections so you can read
# before moving on. Safe to re-run any time -- it always starts by wiping
# its own output/. Each run partitions a real, fairly dense charge-density
# grid, so this takes a little while (tens of seconds per run) -- that's
# PyBader working, not a hang.
#
# Sn3O4.out/.RHO/.XV are a real, finished SIESTA calculation on a 14-atom
# Sn3O4 structure, copied from test/3-analysis/6-bader/ -- chosen because
# it genuinely exercises several diagnostics at once: donor Sn vs. acceptor
# O, a real near-zero-population case on 2 of the 6 Sn atoms, and enough
# symmetry for the equivalent-atom cross-check to run.

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
echo " What stb-bader computes, and why"
echo "=================================================================="
cat <<'EOF'
stb-bader reads a finished SIESTA calculation's charge density (.RHO)
plus its geometry (.XV or .fdf), converts it to a .cube file via sisl,
and partitions it into atomic (Bader) basins with PyBader -- reporting,
per atom:
  - Bader population (electrons integrated in its own basin)
  - net charge (Z_val - population) and donor/acceptor state
  - Bader volume and minimum surface distance
  - net spin, if the .RHO is spin-polarized
Plus a per-species mean/std summary and a cross-check that
symmetry-equivalent atoms agree on charge.

There is no quantum-mechanical operator for "how much charge belongs to
this atom" -- Bader's Atoms-in-Molecules theory divides the density into
regions bounded by zero-flux surfaces (where the density gradient has no
component crossing the surface) instead. See the README for the full
theory, the grid-based algorithm's own approximation, and every
documented limitation (non-nuclear attractors, Z_val accuracy, the
near-zero-population red flag, --speed fast's tradeoff).

Every run prints a numbered [0]...[7] report, the same style every newer
tool in this suite uses. --save-report additionally persists that report
to stb_bader_report.txt -- off by default, so a plain run only ever
writes the .cube file(s) and references.bib, no text report file. The
tool used to always write <label>_BADER.txt unconditionally; that's gone
now (a stale one from an older run is cleaned up automatically instead).
EOF
pause

echo "=================================================================="
echo " output/basic/  --  Z_val detection, per-atom populations, per-species summary"
echo "=================================================================="
cat <<'EOF'
Watch [1] detect Z_val straight from Sn3O4.out (Sn=4, O=6 valence
electrons -- the actual pseudopotentials used, not a generic guess),
[3] report every atom's population/net charge/donor-or-acceptor state,
and [5] flag atoms #1/#2 (Sn) with an essentially-zero Bader population
-- a real, genuine case in this fixture, not a synthetic example: PyBader
found almost no density of its own to anchor a basin on there, so their
whole region got folded into a neighboring basin instead (see the README
for why this happens and why it's a red flag, not a precise result):
EOF
mkdir -p "$OUT/basic"
cp Sn3O4.out Sn3O4.RHO Sn3O4.XV "$OUT/basic/"
echo
echo "\$ stb-bader --label Sn3O4 --no-intro"
(cd "$OUT/basic" && timeout 120 stb-bader --label Sn3O4 --no-intro > console.log 2>&1)
awk '/\[1\] VALENCE/{flag=1} /\[2\] PYBADER/{flag=0} flag' "$OUT/basic/console.log"
awk '/\[3\] PER-ATOM/{flag=1} /\[4\] PER-SPECIES/{flag=0} flag' "$OUT/basic/console.log"
awk '/\[5\] DIAGNOSTICS/{flag=1} /\[6\] REFERENCES/{flag=0} flag' "$OUT/basic/console.log"
pause

echo "=================================================================="
echo " output/speed-fast/  --  --speed fast: PyBader's on-grid method, less precise edges"
echo "=================================================================="
cat <<'EOF'
--speed fast switches PyBader from the default near-grid method (this
tool's normal mode -- corrects a known on-grid bias where basin edges can
spuriously align with the grid axes) to the plain, faster on-grid method.
Useful for a quick look; not recommended for a number you intend to quote:
EOF
mkdir -p "$OUT/speed-fast"
cp Sn3O4.out Sn3O4.RHO Sn3O4.XV "$OUT/speed-fast/"
echo
echo "\$ stb-bader --label Sn3O4 --speed fast --no-intro"
(cd "$OUT/speed-fast" && timeout 120 stb-bader --label Sn3O4 --speed fast --no-intro > console.log 2>&1)
grep "CAUTION" "$OUT/speed-fast/console.log"
grep "Method " "$OUT/speed-fast/console.log"
pause

echo "=================================================================="
echo " output/export-volumes/  --  one .cube per atom, for VESTA/VMD"
echo "=================================================================="
cat <<'EOF'
--export-volumes writes each atom's own Bader volume as a separate
Bader-atoms-<N>.cube file, inside --output-dir -- for visual inspection
of the actual basin shapes in VESTA/VMD (one file per atom, each as large
as the main grid, so this is opt-in):
EOF
mkdir -p "$OUT/export-volumes"
cp Sn3O4.out Sn3O4.RHO Sn3O4.XV "$OUT/export-volumes/"
echo
echo "\$ stb-bader --label Sn3O4 --export-volumes --no-intro"
(cd "$OUT/export-volumes" && timeout 120 stb-bader --label Sn3O4 --export-volumes --no-intro > console.log 2>&1)
grep "Exported" "$OUT/export-volumes/console.log"
echo "Files written:"
ls "$OUT/export-volumes/" | grep "Bader-atoms" | wc -l
echo "Bader-atoms-<N>.cube files (one per atom)."
pause

echo "=================================================================="
echo " output/full-report/  --  --save-report (off by default)"
echo "=================================================================="
cat <<'EOF'
Every run always prints the numbered [0]...[7] report to the console.
--save-report additionally persists it to stb_bader_report.txt -- off by
default, so a plain run only ever writes the .cube file and
references.bib, no text report file at all:
EOF
mkdir -p "$OUT/full-report"
cp Sn3O4.out Sn3O4.RHO Sn3O4.XV "$OUT/full-report/"
echo
echo "\$ stb-bader --label Sn3O4 --no-intro   # default: no report file"
(cd "$OUT/full-report" && timeout 120 stb-bader --label Sn3O4 --no-intro > console_default.log 2>&1)
ls "$OUT/full-report/" | grep -v "\.out$\|\.RHO$\|\.XV$\|\.cube$\|console"
echo "(no <label>_BADER.txt, no stb_bader_report.txt -- only references.bib)"
echo
echo "\$ stb-bader --label Sn3O4 --save-report --no-intro"
(cd "$OUT/full-report" && timeout 120 stb-bader --label Sn3O4 --save-report --no-intro > console_saved.log 2>&1)
echo "Report sections written to stb_bader_report.txt:"
grep -E "^\[[0-9]+\]" "$OUT/full-report/stb_bader_report.txt"
echo
echo "references.bib -- SIESTA (every stb-bader run analyzes a finished SIESTA calculation):"
grep "^@" "$OUT/full-report/references.bib"
pause

echo "=================================================================="
echo " Two ways to run it"
echo "=================================================================="
cat <<'EOF'
A -- direct CLI:
  stb-bader --label Sn3O4

B -- interactive stb-suite menu:
  $ stb-suite
  Select an option (0-6, or a tool code like 4.1.2): 3.6

Both paths call the exact same underlying tool -- proven directly below.
EOF
TMP="$(mktemp -d)"
cp Sn3O4.out Sn3O4.RHO Sn3O4.XV "$TMP/"
echo
echo "\$ printf '3.6\\nSn3O4\\n\\n\\n1\\n\\n\\ny\\nn\\nn\\n\\n0\\n' | stb-suite     # default output dir, default ref, speed=normal, no vacuum-tol/threads, keep cube, no export, no report"
(cd "$TMP" && printf '3.6\nSn3O4\n\n\n1\n\n\ny\nn\nn\n\n0\n' | timeout 120 stb-suite > session.log 2>&1) || true
CLI_LINE=$(grep "Total Integrated" "$OUT/basic/console.log" | head -1)
MENU_LINE=$(grep "Total Integrated" "$TMP/session.log" | head -1)
if [ "$CLI_LINE" = "$MENU_LINE" ]; then
    echo "Confirmed: identical Total Integrated line from the CLI and the interactive menu."
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
  basic/            speed-fast/
  export-volumes/   full-report/

Each has references.bib and the .cube file(s); full-report/ additionally
has stb_bader_report.txt (only from its --save-report run).

Recap of what this walkthrough covered:
  - Bader's Atoms-in-Molecules theory and the grid-based algorithm's own
    approximation (see the README for the full theory/limitations)
  - Z_val auto-detection from the real .out pseudopotential log
  - the per-atom population/net-charge/donor-or-acceptor table
  - a real near-zero-population red flag on this fixture's own Sn #1/#2
  - --speed fast's on-grid vs. the default near-grid method
  - --export-volumes for per-atom VESTA/VMD visualization
  - the numbered [0]...[7] report, --save-report, references.bib
  - CLI and the interactive stb-suite menu building the same command

As a next step, try on your own with a finished SIESTA calculation:
  stb-bader --label my_calc --save-report
  stb-bader --label my_calc --ref relax/my_calc.out --export-volumes
EOF
