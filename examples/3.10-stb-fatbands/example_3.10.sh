#!/bin/bash
# Guided example: stb-fatbands (code 3.10 in the stb-suite menu)
#
# Not an automated test (see test/3-analysis/10-fatbands/test.sh for that) --
# a commented walk-through: it runs real commands, one group at a time,
# into its own output/<case>/ folder, and shows you the piece of output
# that proves what just happened. Pauses between sections so you can read
# before moving on. Safe to re-run any time -- it always starts by wiping
# its own output/.
#
# Two real (not synthetic) SIESTA fixtures are used, both copied/derived
# from test/3-analysis/10-fatbands/ -- a .WFSX/.HSX pair cannot be
# meaningfully hand-crafted the way a small .bands/.PDOS.xml can (see
# examples/3.1-stb-bands, examples/3.2-stb-dos):
#   - Sn3O4.bands(.WFSX)/Sn3O4.HSX: a real, converged bulk Sn3O4 cell, a
#     deliberately SHORT band path (16 k-points) so the .WFSX stays a few
#     MB instead of the 100+ MB a finely-sampled path would produce.
#   - spin/Ospin.*: a real, freshly-run spin-polarized isolated O atom,
#     converged to its physical 2 Bohr-magneton triplet moment -- see the
#     README for why this exists and what it demonstrates.

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
echo " What a fatbands plot shows, and why WFS.Write.For.Bands is required"
echo "=================================================================="
cat <<'EOF'
SIESTA is an LCAO code: every band's wavefunction is a linear combination
of localized atomic orbitals,

    psi_n(k, r) = sum_mu  c_mu,n(k) * phi_mu(r)

The complex coefficients c_mu,n(k) -- one per (orbital, band, k-point) --
are exactly what a .WFSX file stores. A plain .bands file only has the
eigenvalues (WHERE a state sits in energy); stb-fatbands additionally
needs the .WFSX to know WHAT that state is made of, so it can color/size
each point by orbital, atom, or species character.

SIESTA only writes this file if explicitly asked -- it can be large (one
complex number per orbital/band/k-point), so it is never saved by
default. The fdf setting that matters:

    WFS.Write.For.Bands T
    BandLinesScale ReciprocalLatticeVectors
    %block BandLines
      1   0.000  0.000  0.000  GAMMA
      5   0.500  0.000  0.000  X
    %endblock BandLines

Without it, SIESTA still writes <label>.bands but NO <label>.bands.WFSX --
stb-fatbands refuses to run at all in that case (a clear error), rather
than silently doing nothing. See the real spin/Ospin.fdf in this folder
for this exact block in a complete, minimal, working input:
EOF
grep -A4 "WFS.Write.For.Bands" spin/Ospin.fdf
pause

echo "=================================================================="
echo " output/species-l-default/  --  the default projection, accurate weights"
echo "=================================================================="
cat <<'EOF'
The orbital weight itself is a Mulliken-like population: not the naive
|c_mu|^2, but Re[c*_mu * sum_nu S_mu,nu * c_nu] -- SIESTA's numerical
atomic orbitals are NOT orthogonal to each other, so the overlap matrix S
(from a companion .HSX file, SaveHS T) matters. Watch [1] INPUT DATA
report which one was actually used ("accurate" here, since Sn3O4.HSX is
right next to the .bands/.WFSX), and [3] ORBITAL PROJECTION show the
default species_l categories -- species AND angular momentum combined
(e.g. 'Sn-s', 'O-p'), the most informative single view for chemical
bonding character:
EOF
mkdir -p "$OUT/species-l-default"
cp Sn3O4.bands Sn3O4.bands.WFSX Sn3O4.HSX "$OUT/species-l-default/"
echo
echo "\$ stb-fatbands --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX --shift fermi --save-gnuplot --no-intro"
(cd "$OUT/species-l-default" && stb-fatbands --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX \
    --hsx-file Sn3O4.HSX --shift fermi --save-gnuplot --no-intro > console.log 2>&1)
awk '/\[1\] INPUT DATA/{flag=1} /\[2\]/{flag=0} flag' "$OUT/species-l-default/console.log"
awk '/\[3\] ORBITAL PROJECTION/{flag=1} /\[4\]/{flag=0} flag' "$OUT/species-l-default/console.log"
pause

echo "=================================================================="
echo " output/projection-modes/  --  l vs. species vs. species_l, same data"
echo "=================================================================="
cat <<'EOF'
Same structure, three --projection modes: l (orbital character only), atom
species (which element only), species_l (both combined). Watch how many
categories each finds, and that species_l alone tells you BOTH which
element AND which orbital shell dominate a band -- something neither l nor
species can show on its own:
EOF
mkdir -p "$OUT/projection-modes"
cp Sn3O4.bands Sn3O4.bands.WFSX Sn3O4.HSX "$OUT/projection-modes/"
for proj in l species species_l; do
    echo
    echo "\$ stb-fatbands --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX \\"
    echo "      --shift fermi --projection $proj --no-intro"
    (cd "$OUT/projection-modes" && stb-fatbands --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX \
        --hsx-file Sn3O4.HSX --shift fermi --projection "$proj" --no-intro > "console_$proj.log" 2>&1)
    grep "Categories found" "$OUT/projection-modes/console_$proj.log"
done
pause

echo "=================================================================="
echo " output/category-filter/  --  --category restricts which series are plotted/saved"
echo "=================================================================="
cat <<'EOF'
--category keeps only the requested categories (any --projection mode) --
here, only the Sn-s and O-p characters, out of species_l's full 6:
EOF
mkdir -p "$OUT/category-filter"
cp Sn3O4.bands Sn3O4.bands.WFSX Sn3O4.HSX "$OUT/category-filter/"
echo
echo "\$ stb-fatbands --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX \\"
echo "      --shift fermi --projection species_l --category Sn-s O-p --save-gnuplot --no-intro"
(cd "$OUT/category-filter" && stb-fatbands --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX \
    --hsx-file Sn3O4.HSX --shift fermi --projection species_l --category Sn-s O-p \
    --save-gnuplot --no-intro > console.log 2>&1)
grep "Categories used" "$OUT/category-filter/console.log"
echo "Files written:"
ls "$OUT/category-filter/"*.dat
pause

echo "=================================================================="
echo " output/accuracy-fallback/  --  no .HSX: the approximate |c|^2 weights"
echo "=================================================================="
cat <<'EOF'
Without a .HSX (no real overlap matrix available), stb-fatbands falls back
to an implicit-orthogonal-basis approximation instead of refusing to run --
same order of magnitude, not exact, and it always tells you so. Here we
give --geometry-file calc.fdf (the real SIESTA input, which also needs its
own structure.fdf + Sn/O .ion files alongside it) instead of a .HSX:
EOF
mkdir -p "$OUT/accuracy-fallback"
cp Sn3O4.bands Sn3O4.bands.WFSX calc.fdf structure.fdf Sn.ion Sn.ion.xml O.ion O.ion.xml \
   "$OUT/accuracy-fallback/"
echo
echo "\$ stb-fatbands --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --geometry-file calc.fdf \\"
echo "      --shift fermi --no-intro"
(cd "$OUT/accuracy-fallback" && stb-fatbands --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX \
    --geometry-file calc.fdf --shift fermi --no-intro > console.log 2>&1)
grep -A1 "WARNING.*No .HSX found" "$OUT/accuracy-fallback/console.log"
grep "Weight accuracy" "$OUT/accuracy-fallback/console.log"
pause

echo "=================================================================="
echo " output/spin-polarized/  --  a REAL spin-polarized calculation"
echo "=================================================================="
cat <<'EOF'
spin/Ospin.*: a single O atom in a large vacuum box, Spin polarized, seeded
via %block DM.InitSpin to converge to its physical ground state -- this
calculation's own log confirms a converged spin moment |S| = 2.00000 (the
textbook triplet oxygen atom). [1] INPUT DATA reports 2 spin channels;
[2] BAND GAP ANALYSIS reports spin-up and spin-down VBM/CBM SEPARATELY
(they differ enormously here -- a real, physical consequence of the
spin-split occupation, not a bug); and [3] ORBITAL PROJECTION splits every
species_l category into _up/_down. This exact check caught a real bug
while building this tool: the original code merged both spin channels
into ONE series per category, silently averaging together two band sets
whose energies differ by tens of eV:
EOF
mkdir -p "$OUT/spin-polarized"
cp spin/Ospin.fdf spin/Ospin.bands spin/Ospin.bands.WFSX spin/Ospin.HSX spin/O.ion spin/O.ion.xml \
   "$OUT/spin-polarized/"
echo
echo "\$ stb-fatbands --label Ospin --shift fermi --save-gnuplot --no-intro"
(cd "$OUT/spin-polarized" && stb-fatbands --label Ospin --shift fermi \
    --save-gnuplot --no-intro > console.log 2>&1)
awk '/\[1\] INPUT DATA/{flag=1} /\[2\]/{flag=0} flag' "$OUT/spin-polarized/console.log"
awk '/\[2\] BAND GAP ANALYSIS/{flag=1} /\[3\]/{flag=0} flag' "$OUT/spin-polarized/console.log"
awk '/\[3\] ORBITAL PROJECTION/{flag=1} /\[4\]/{flag=0} flag' "$OUT/spin-polarized/console.log"
echo
echo "Data files -- one per (category, spin), never merged:"
ls "$OUT/spin-polarized/"fatbands_O-s_*.dat "$OUT/spin-polarized/"fatbands_O-p_*.dat
pause

echo "=================================================================="
echo " output/full-report/  --  --save-report / --save-gnuplot (off by default)"
echo "=================================================================="
cat <<'EOF'
Every run always prints the numbered [0]...[6] report to the console.
--save-report additionally persists it to stb_fatbands_report.txt, and
--save-gnuplot writes fatbands_<category>.dat/fatbands.gplot -- both off
by default, so a plain run only ever writes references.bib:
EOF
mkdir -p "$OUT/full-report"
cp Sn3O4.bands Sn3O4.bands.WFSX Sn3O4.HSX "$OUT/full-report/"
echo
echo "\$ stb-fatbands --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX \\"
echo "      --shift fermi --no-intro   # default: no report, no data/gnuplot"
(cd "$OUT/full-report" && stb-fatbands --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX \
    --hsx-file Sn3O4.HSX --shift fermi --no-intro > console_default.log 2>&1)
ls "$OUT/full-report/" | grep -v "^Sn3O4\|console"
echo "(no fatbands_*.dat/.gplot, no stb_fatbands_report.txt -- only references.bib)"
echo
echo "\$ stb-fatbands --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX \\"
echo "      --shift fermi --save-report --save-gnuplot --no-intro"
(cd "$OUT/full-report" && stb-fatbands --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX \
    --hsx-file Sn3O4.HSX --shift fermi --save-report --save-gnuplot --no-intro \
    > console_saved.log 2>&1)
echo "Report sections written to stb_fatbands_report.txt:"
grep -E "^\[[0-9]+\]" "$OUT/full-report/stb_fatbands_report.txt"
echo
echo "references.bib -- SIESTA:"
grep "^@" "$OUT/full-report/references.bib"
pause

echo "=================================================================="
echo " Two ways to run it"
echo "=================================================================="
cat <<'EOF'
A -- direct CLI:
  stb-fatbands --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX --shift fermi

B -- interactive stb-suite menu:
  $ stb-suite
  Select an option (0-6, or a tool code like 4.1.2): 3.10

Both paths call the exact same underlying tool -- proven directly below.
The menu defaults --shift to Fermi level and --projection to species_l
(press Enter through both).
EOF
TMP="$(mktemp -d)"
cp Sn3O4.bands Sn3O4.bands.WFSX Sn3O4.HSX "$TMP/"
echo
echo "\$ printf '3.10\\nSn3O4\\n\\n\\n\\nn\\nn\\nn\\n' | stb-suite     # defaults: fermi, species_l, no save/view"
(cd "$TMP" && printf '3.10\nSn3O4\n\n\n\nn\nn\nn\n' | stb-suite > session.log 2>&1) || true
CLI_LINE=$(grep "Best (indirect) gap" "$OUT/species-l-default/console.log")
MENU_LINE=$(grep "Best (indirect) gap" "$TMP/session.log")
if [ "$CLI_LINE" = "$MENU_LINE" ]; then
    echo "Confirmed: identical gap result from the CLI and the interactive menu."
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
Not exercised interactively by this script (needs a display, though it IS
covered headlessly by test/3-analysis/10-fatbands/test.sh under
MPLBACKEND=Agg). --view shows the fatbands scatter plot (single colorbar
-encoded series for one category, or one discrete-colored series per
category with a shared legend for several) and blocks until you close it:

  stb-fatbands --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX --shift fermi --view
  stb-fatbands --label Ospin --shift fermi --view   # spin-split series, from spin/
EOF
pause

echo "=================================================================="
echo " Done"
echo "=================================================================="
cat <<EOF
Six self-contained folders were generated under output/:
  species-l-default/   projection-modes/   category-filter/
  accuracy-fallback/    spin-polarized/      full-report/

Each has references.bib; the --save-gnuplot runs additionally have
fatbands_<category>.dat/fatbands.gplot, and full-report/ has
stb_fatbands_report.txt (only from its --save-report run).

Recap of what this walkthrough covered:
  - LCAO wavefunctions, the .WFSX coefficients, and why
    WFS.Write.For.Bands T is required in the SIESTA calculation -- no
    .WFSX means stb-fatbands cannot run at all
  - the Mulliken-like orbital weight (needs the .HSX overlap matrix for
    exact values; an approximate |c|^2 fallback otherwise, always labeled)
  - the five --projection modes, and why species_l (species + angular
    momentum combined) is the default -- more informative than either
    species or l alone
  - --category restricting the output to specific characters
  - a REAL spin-polarized calculation: [2] BAND GAP ANALYSIS reporting
    spin-up/spin-down separately, and [3] ORBITAL PROJECTION splitting
    every category into _up/_down -- the exact check that caught a real
    bug (both spin channels silently merged into one series) while
    building this tool
  - the numbered [0]...[6] report, --save-report, --save-gnuplot,
    references.bib
  - CLI and the interactive stb-suite menu building the same command

As a next step, try on your own with a real SIESTA calculation (remember
WFS.Write.For.Bands T + SaveHS T in the fdf):
  stb-fatbands --label my_calc --shift fermi --save-report --save-gnuplot
  stb-fatbands --label my_magnetic_calc --shift fermi --view
EOF
