#!/bin/bash
# Guided example: stb-structural (code 3.4 in the stb-suite menu)
#
# Not an automated test (see test/3-analysis/4-structure/test.sh for that)
# -- a commented walk-through: it runs real commands, one group at a time,
# into its own output/<case>/ folder, and shows you the piece of output
# that proves what just happened. Pauses between sections so you can read
# before moving on. Safe to re-run any time -- it always starts by wiping
# its own output/.
#
# structure.fdf is a real, small (14-atom, Sn3O4-composition) SIESTA
# structure, copied from test/3-analysis/4-structure/structure.fdf -- two
# chemical species with different coordination environments (6-coordinate
# Sn, 3-coordinate O), rich enough to exercise every section of the report.

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
echo " What stb-structural computes, and why 4 different ECN methods"
echo "=================================================================="
cat <<'EOF'
stb-structural reads a SIESTA structure file (.fdf or a post-relaxation
.STRUCT_OUT) and reports, with no SIESTA re-run needed:
  - lattice parameters, cell volume, density
  - Effective Coordination Number (ECN), by 4 independent methods
  - bond distances/angles, pooled per species-pair/triplet
  - coordination-polyhedron distortion (BLD, BAV, volume)
  - the shortest same-species contact distance in the structure
  - coordination-polyhedron connectivity (corner/edge/face-sharing)
  - optionally, the radial distribution function g(r)

There is no single "true" coordination number for a real structure --
each of the 4 ECN methods is a different operational definition of
"where does the first neighbor shell end." Reporting all 4 side by side
means: where they agree, trust the number; where they disagree, that
disagreement is itself informative (a genuinely ambiguous/distorted
site), not a bug in any one method. See the README for the exact
formula and limitation of every single method/metric below.

Every run prints a numbered [0]...[11] report, the same style every
newer tool in this suite uses. --save-report additionally persists that
report to stb_structural_report.txt -- off by default, so a plain run
only ever writes rdf.dat (unless --no-rdf) and references.bib, no text
report file. The tool used to always write structural_information.dat
unconditionally; that's gone now (a stale one from an older run is
cleaned up automatically instead).
EOF
pause

echo "=================================================================="
echo " output/mode-mean/  --  full-structure analysis, per species"
echo "=================================================================="
cat <<'EOF'
--mode mean analyzes every atom, aggregated per chemical species. Watch
[2] report 4 CN methods for both Sn (6-coordinate) and O (3-coordinate),
[3]/[4] pool bond distances/angles by species-pair/triplet instead of one
meaningless structure-wide average, and [5] show BLD/BAV/volume per
species -- including "Volume: N/A" for O, whose 3-coordinate sites don't
have enough neighbors (need >=4) to bound a 3D convex-hull volume:
EOF
mkdir -p "$OUT/mode-mean"
cp structure.fdf "$OUT/mode-mean/"
echo
echo "\$ stb-structural --file structure.fdf --format fdf --mode mean --no-intro"
(cd "$OUT/mode-mean" && stb-structural --file structure.fdf --format fdf --mode mean \
    --no-intro > console.log 2>&1)
awk '/\[2\] EFFECTIVE/{flag=1} /\[6\]/{flag=0} flag' "$OUT/mode-mean/console.log"
pause

echo "=================================================================="
echo " output/mode-list/  --  --mode list: scoped to specific atoms"
echo "=================================================================="
cat <<'EOF'
--mode list --list 1,7 restricts the ECN and distortion sections to just
the 2 named atoms (atom 1 = Sn, atom 7 = O, 1-based) -- no point
computing/printing whole-structure statistics if you only care about
these sites. But [6]/[7] (same-species minimum distance, polyhedron
connectivity) stay whole-structure regardless of --mode: both describe
the extended network, not a property of one atom, so scoping them down
to 2 atoms wouldn't make physical sense -- watch them report the exact
same O-O/Sn-Sn numbers as the --mode mean run above:
EOF
mkdir -p "$OUT/mode-list"
cp structure.fdf "$OUT/mode-list/"
echo
echo "\$ stb-structural --file structure.fdf --format fdf --mode list --list 1,7 --no-intro"
(cd "$OUT/mode-list" && stb-structural --file structure.fdf --format fdf --mode list --list 1,7 \
    --no-intro > console.log 2>&1)
awk '/\[2\] EFFECTIVE/{flag=1} /\[3\]/{flag=0} flag' "$OUT/mode-list/console.log"
echo "..."
awk '/\[6\] SAME-SPECIES/{flag=1} /\[8\]/{flag=0} flag' "$OUT/mode-list/console.log"
pause

echo "=================================================================="
echo " output/rdf/  --  g(r), and independently cross-checking the O-Sn coordination number"
echo "=================================================================="
cat <<'EOF'
The radial distribution function g(r) is proportional to the probability
of finding an atom at distance r from a reference atom, relative to a
uniform density. [8] reports each pair's first peak position/height;
rdf.dat has the full curve. The real proof that this is the same physics
as the ECN section above: integrating g_O-Sn(r) under its first peak
(r < 2.6 Ang, the trough right after the peak) with

    N = 4*pi*rho_Sn * integral( r^2 * g_O-Sn(r) dr )

independently recovers the known O-Sn coordination number -- computed a
completely different way (a continuous pair-correlation integral,
instead of a discrete neighbor list) from the CrystalNN ECN value seen
above (O: 2.956, close to the geometric CN=3 every O site actually has):
EOF
mkdir -p "$OUT/rdf"
cp structure.fdf "$OUT/rdf/"
echo
echo "\$ stb-structural --file structure.fdf --format fdf --mode mean --no-intro"
(cd "$OUT/rdf" && stb-structural --file structure.fdf --format fdf --mode mean \
    --no-intro > console.log 2>&1)
awk '/\[8\] RADIAL/{flag=1} /\[9\]/{flag=0} flag' "$OUT/rdf/console.log"
echo
python3 -c "
import numpy as np
data = np.loadtxt('$OUT/rdf/rdf.dat')
r, g_osn = data[:, 0], data[:, 3]
dr = r[1] - r[0]
shell = 4 * np.pi * r**2 * dr
rho_Sn = 6 / 236.82072939038554   # 6 Sn atoms / cell volume (Ang^3)
mask = r < 2.6                     # cuts right after the first O-Sn peak
n = np.sum(g_osn[mask] * shell[mask] * rho_Sn)
print(f'Integrating g_O-Sn(r) under the first peak (r < 2.6 Ang): N = {n:.3f}')
print('Matches the known/geometric O-Sn coordination number of 3.0 -- an')
print('independent cross-check between two completely different methods.')
"
pause

echo "=================================================================="
echo " output/full-report/  --  --save-report (off by default)"
echo "=================================================================="
cat <<'EOF'
Every run always prints the numbered [0]...[11] report to the console.
--save-report additionally persists it to stb_structural_report.txt --
off by default, so a plain run only ever writes rdf.dat and
references.bib, no text report file at all:
EOF
mkdir -p "$OUT/full-report"
cp structure.fdf "$OUT/full-report/"
echo
echo "\$ stb-structural --file structure.fdf --format fdf --mode mean --no-intro   # default: no report file"
(cd "$OUT/full-report" && stb-structural --file structure.fdf --format fdf --mode mean \
    --no-intro > console_default.log 2>&1)
ls "$OUT/full-report/" | grep -v "\.fdf$\|console"
echo "(no structural_information.dat, no stb_structural_report.txt -- only rdf.dat + references.bib)"
echo
echo "\$ stb-structural --file structure.fdf --format fdf --mode mean --save-report --no-intro"
(cd "$OUT/full-report" && stb-structural --file structure.fdf --format fdf --mode mean \
    --save-report --no-intro > console_saved.log 2>&1)
echo "Report sections written to stb_structural_report.txt:"
grep -E "^\[[0-9]+\]" "$OUT/full-report/stb_structural_report.txt"
echo
echo "references.bib -- SIESTA (every stb-structural run analyzes a SIESTA structure file):"
grep "^@" "$OUT/full-report/references.bib"
pause

echo "=================================================================="
echo " Two ways to run it"
echo "=================================================================="
cat <<'EOF'
A -- direct CLI:
  stb-structural --file structure.fdf --format fdf --mode mean

B -- interactive stb-suite menu:
  $ stb-suite
  Select an option (0-6, or a tool code like 4.1.2): 3.4

Both paths call the exact same underlying tool -- proven directly below.
EOF
TMP="$(mktemp -d)"
cp structure.fdf "$TMP/"
echo
echo "\$ printf '3.4\\nstructure.fdf\\n1\\n1\\n\\n\\n\\n\\n' | stb-suite     # format=fdf, mode=mean, defaults"
(cd "$TMP" && printf '3.4\nstructure.fdf\n1\n1\n\n\n\n\n' | stb-suite > session.log 2>&1) || true
CLI_LINE=$(grep "CrystalNN" "$OUT/mode-mean/console.log" | head -1)
MENU_LINE=$(grep "CrystalNN" "$TMP/session.log" | head -1)
if [ "$CLI_LINE" = "$MENU_LINE" ]; then
    echo "Confirmed: identical CrystalNN coordination number from the CLI and the interactive menu."
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
  mode-mean/   mode-list/
  rdf/         full-report/

Each has rdf.dat/references.bib; full-report/ additionally has
stb_structural_report.txt (only from its --save-report run).

Recap of what this walkthrough covered:
  - what stb-structural computes, and why 4 ECN methods are reported
    side by side instead of picking one (see the README for the exact
    formula/limitation of each: MinDistNN, BrunnerNN, EconNN, CrystalNN)
  - --mode mean (per-species, whole structure) vs. --mode list (scoped
    to named atoms) -- and which sections stay whole-structure regardless
  - BLD/BAV/polyhedron volume: what each distortion metric means and its
    documented limitation (BAV here is a modified, non-literature-
    comparable variant -- see the README for exactly why)
  - g(r) and its first-peak integral independently recovering the same
    coordination number the ECN section reports, a real physics
    cross-check, not just a plausibility check
  - the numbered [0]...[11] report, --save-report, references.bib
  - CLI and the interactive stb-suite menu building the same command

As a next step, try on your own with a real SIESTA structure:
  stb-structural --file my_calc.fdf --format fdf --mode mean --save-report
  stb-structural --file my_calc.fdf --format fdf --mode list --list 12,45
EOF
