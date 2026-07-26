#!/bin/bash
# Guided example: stb-bands (code 3.1 in the stb-suite menu)
#
# Not an automated test (see test/3-analysis/1-bands/test.sh for that) --
# a commented walk-through: it runs real commands, one group at a time,
# into its own output/<case>/ folder, and shows you the piece of output
# that proves what just happened. Pauses between sections so you can read
# before moving on. Safe to re-run any time -- it always starts by wiping
# its own output/.
#
# First example for the Analysis category (3.x) -- fixtures here are small
# synthetic .bands/.EIG/.KP files (copied from test/3-analysis/1-bands/,
# renamed for clarity), not a real SIESTA calculation's heavy output: a
# real .bands file for a decent-sized cell is easily 1+ MB (hundreds of
# bands x hundreds of k-points), far too heavy for a lightweight example.
# These tiny files were purpose-built to demonstrate exactly the physics
# below (a real, if small, indirect-gap band structure; a synthetic
# half-metal; a synthetic 3-point k-mesh) without the weight.

set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

# main() calls plt.show() at the very end (an interactive matplotlib
# preview) -- MPLBACKEND=Agg makes that a no-op instead of blocking on a
# GUI window, same convention test.sh itself already uses.
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
echo " What a SIESTA .bands file is, and what stb-bands does with it"
echo "=================================================================="
cat <<'EOF'
A SIESTA .bands file holds the eigenvalues (energies) of every electronic
band, at every k-point along a chosen high-symmetry PATH through the
Brillouin zone (e.g. Gamma -> X -> M -> Gamma) -- exactly what you need to
PLOT a band structure. Its layout: line 1 is the Fermi energy; line 4 is
"nbands nspin nk"; then one block per k-point (k value + nbands*nspin
energies, wrapped onto continuation lines if there are many bands); a
footer lists which k-points are the high-symmetry ones, with labels.

stb-bands parses this directly (no SIESTA re-run needed) and reports:
  - the Valence Band Maximum (VBM) and Conduction Band Minimum (CBM),
    and whether the material is a semiconductor/insulator or a metal
  - whether the gap is DIRECT (VBM and CBM at the same k-point) or
    INDIRECT (different k-points)
  - (optionally) how that k-PATH-based gap compares to the true gap from
    the full SCF k-MESH (a .EIG file) -- the path might miss the real
    extrema entirely
  - (for spin-polarized calculations) each spin channel separately, and
    whether the material is HALF-METALLIC

This walkthrough's manual code+physics review found no bugs in any of
this (see the session notes) -- the numbers were already correct. What
changed here is how the report is presented: a numbered [0]...[7] report
like every 1.x/2.x tool now has, and an opt-in --save-report instead of
always writing bands_analysis.txt (removed) whether you wanted a file or
not.
EOF
pause

echo "=================================================================="
echo " VBM/CBM and direct vs. indirect gaps -- a real (if tiny) example"
echo "=================================================================="
cat <<'EOF'
semiconductor.bands: a real 21-band structure along Gamma -> X -> M ->
Gamma (4 k-points only, kept tiny on purpose). Watch [2] BAND GAP ANALYSIS
report the VBM/CBM at their own k-points, and both the fundamental
(indirect) gap and the same-k (direct) gap -- these usually differ, and
the smaller of the two is never the direct one (the indirect gap is a
minimum over ALL k-pairs; the direct gap is a minimum only over PAIRS AT
THE SAME k, a more restrictive search, so indirect <= direct always):
EOF
mkdir -p "$OUT/gap-analysis"
cp semiconductor.bands "$OUT/gap-analysis/"
echo
echo "\$ stb-bands --file semiconductor.bands --shift vbm --no-intro"
(cd "$OUT/gap-analysis" && stb-bands --file semiconductor.bands --shift vbm \
    --no-intro > console.log 2>&1)
awk '/\[2\] BAND GAP ANALYSIS/{flag=1} /\[3\]|\[4\]/{flag=0} flag' "$OUT/gap-analysis/console.log"
pause

echo "=================================================================="
echo " output/metallic-threshold/  --  --gap-tol decides Metallic vs. Indirect"
echo "=================================================================="
cat <<'EOF'
A "gap" can never be negative (bands crossing the Fermi level just means
CBM < VBM at the global level, clamped to exactly 0.0 -- see the review
notes), but a small POSITIVE gap and a true metal look identical without
a threshold. --gap-tol (default 0.01 eV) is that threshold: below it, the
material is classified "Metallic" instead of "Indirect"/"Direct". Same
file, only --gap-tol changes -- watch the classification flip live:
EOF
mkdir -p "$OUT/metallic-threshold"
cp semiconductor.bands "$OUT/metallic-threshold/"
echo
echo "\$ stb-bands --file semiconductor.bands --shift vbm --gap-tol 0.0001 --no-intro   # tight"
(cd "$OUT/metallic-threshold" && stb-bands --file semiconductor.bands --shift vbm \
    --gap-tol 0.0001 --no-intro > console_tight.log 2>&1)
grep "Gap type" "$OUT/metallic-threshold/console_tight.log"
echo
echo "\$ stb-bands --file semiconductor.bands --shift vbm --gap-tol 5.0 --no-intro      # very loose"
(cd "$OUT/metallic-threshold" && stb-bands --file semiconductor.bands --shift vbm \
    --gap-tol 5.0 --no-intro > console_loose.log 2>&1)
grep "Gap type" "$OUT/metallic-threshold/console_loose.log"
echo
echo "The real gap (~1.25 eV) never changed -- only which side of --gap-tol"
echo "it falls on. A --gap-tol of 5.0 eV would swallow almost any real"
echo "semiconductor gap into 'Metallic' -- --gap-tol needs to be small"
echo "compared to the material's expected gap, not a generic default."
pause

echo "=================================================================="
echo " output/half-metallic/  --  one spin channel with a gap, the other without"
echo "=================================================================="
cat <<'EOF'
half_metal.bands: a synthetic spin-polarized (nspin=2) structure built so
the spin-up channel is (nearly) metallic and the spin-down channel has a
clean 2 eV gap -- a textbook HALF-METAL, a material that conducts current
of only one spin orientation. [3] SPIN-POLARIZED ANALYSIS reports each
channel independently plus the half-metallic flag:
EOF
mkdir -p "$OUT/half-metallic"
cp half_metal.bands "$OUT/half-metallic/"
echo
echo "\$ stb-bands --file half_metal.bands --shift fermi --no-intro"
(cd "$OUT/half-metallic" && stb-bands --file half_metal.bands --shift fermi \
    --no-intro > console.log 2>&1)
awk '/\[3\] SPIN-POLARIZED ANALYSIS/{flag=1} /\[4\]/{flag=0} flag' "$OUT/half-metallic/console.log"
pause

echo "=================================================================="
echo " output/mesh-vs-line/  --  a k-path can miss the true VBM/CBM entirely"
echo "=================================================================="
cat <<'EOF'
A .bands file only samples a 1-D PATH through the Brillouin zone -- the
true VBM/CBM might sit at some other k-point entirely, off that path. A
.EIG file (the full SCF k-MESH, denser and unbiased) can catch this.
mesh.EIG/mesh.KP is a tiny synthetic 3-point mesh with known-by-hand
extrema. The key physical fact --eig-file's comparison relies on: a
denser mesh can only find a fundamental gap SMALLER THAN OR EQUAL TO a
1-D subset of it, never larger -- so a mesh gap smaller than the line gap
means "trust the mesh, the path missed the real extrema"; a mesh gap
LARGER than the line gap would be a red flag (mesh too coarse), never
expected from a genuinely denser sampling:
EOF
mkdir -p "$OUT/mesh-vs-line"
cp semiconductor.bands mesh.EIG mesh.KP "$OUT/mesh-vs-line/"
echo
echo "\$ stb-bands --file semiconductor.bands --eig-file mesh.EIG --kp-file mesh.KP \\"
echo "      --shift fermi --no-intro"
(cd "$OUT/mesh-vs-line" && stb-bands --file semiconductor.bands --eig-file mesh.EIG \
    --kp-file mesh.KP --shift fermi --no-intro > console.log 2>&1)
awk '/\[4\] MESH VS LINE COMPARISON/{flag=1} /\[5\]/{flag=0} flag' "$OUT/mesh-vs-line/console.log"
pause

echo "=================================================================="
echo " Two ways to run it"
echo "=================================================================="
cat <<'EOF'
A -- direct CLI:
  stb-bands --file semiconductor.bands --shift fermi

B -- interactive stb-suite menu:
  $ stb-suite
  Select an option (0-6, or a tool code like 4.1.2): 3.1

Every run always writes bands_gnuplot.dat + bands.gplot (the actual
plot data -- run gnuplot on the .gplot for a PDF) and references.bib
(SIESTA). --save-report additionally persists the full numbered console
report to stb_bands_report.txt -- off by default, and (this session's
change) bands_analysis.txt is never written anymore, under any flag.
EOF
mkdir -p "$OUT/full-report"
cp semiconductor.bands "$OUT/full-report/"
echo
echo "\$ stb-bands --file semiconductor.bands --shift fermi --no-intro   # default: no report file"
(cd "$OUT/full-report" && stb-bands --file semiconductor.bands --shift fermi \
    --no-intro > console_default.log 2>&1)
ls "$OUT/full-report/" | grep -v "\.bands$\|console"
echo "(no bands_analysis.txt, no stb_bands_report.txt -- only the plot data + references.bib)"
echo
echo "\$ stb-bands --file semiconductor.bands --shift fermi --save-report --no-intro"
(cd "$OUT/full-report" && stb-bands --file semiconductor.bands --shift fermi \
    --save-report --no-intro > console_saved.log 2>&1)
echo "Report sections written to stb_bands_report.txt:"
grep -E "^\[[0-9]+\]" "$OUT/full-report/stb_bands_report.txt"
echo
echo "references.bib -- SIESTA (every stb-bands run is post-processing a SIESTA output):"
grep "^@" "$OUT/full-report/references.bib"
pause

echo "=================================================================="
echo " Proof: CLI and the interactive stb-suite menu agree"
echo "=================================================================="
echo "Driving the same case through the interactive menu's manual entry"
echo "mode and checking it reaches the same gap classification."
TMP="$(mktemp -d)"
cp semiconductor.bands "$TMP/siesta.bands"
echo
echo "\$ printf '3.1\\nsiesta\\n1\\n.\\nn\\n\\n0\\n' | stb-suite     # '1' = shift by VBM"
(cd "$TMP" && printf '3.1\nsiesta\n1\n.\nn\n\n0\n' | stb-suite > session.log 2>&1) || true
CLI_LINE=$(grep "Gap type" "$OUT/gap-analysis/console.log")
MENU_LINE=$(grep "Gap type" "$TMP/session.log")
if [ "$CLI_LINE" = "$MENU_LINE" ]; then
    echo "Confirmed: identical gap classification from the CLI and the interactive menu."
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
Five self-contained folders were generated under output/:
  gap-analysis/       metallic-threshold/   half-metallic/
  mesh-vs-line/       full-report/

Each has bands_gnuplot.dat/bands.gplot/references.bib; full-report/
additionally has stb_bands_report.txt (only from its --save-report run).

Recap of what this walkthrough covered:
  - what a .bands file is, and VBM/CBM/direct-vs-indirect gap physics on
    a real (if tiny) 21-band indirect-gap structure
  - --gap-tol: the threshold between "Metallic" and a genuine small gap,
    demonstrated flipping live on the exact same data
  - half-metallic character: one spin channel with a gap, the other
    without, on a purpose-built synthetic fixture
  - --eig-file/--kp-file: why a k-mesh can catch a true VBM/CBM a 1-D
    k-path misses, and why the comparison can only go one direction
  - bands_gnuplot.dat/bands.gplot/references.bib always written;
    --save-report (opt-in) for the full numbered text report; the old
    always-on bands_analysis.txt is gone
  - CLI and the interactive stb-suite menu building the same command

Not exercised by this script (needs a display): the interactive
matplotlib preview (main() calls plt.show() at the very end) -- try it
yourself without MPLBACKEND=Agg:
  stb-bands --file semiconductor.bands --shift fermi

As a next step, try on your own with a real SIESTA .bands file:
  stb-bands --label <your_siesta_label> --shift fermi --save-report
  stb-bands --file your.bands --eig-file your.EIG --shift fermi
EOF
