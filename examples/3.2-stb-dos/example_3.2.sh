#!/bin/bash
# Guided example: stb-dos (code 3.2 in the stb-suite menu)
#
# Not an automated test (see test/3-analysis/2-dos/test.sh for that) --
# a commented walk-through: it runs real commands, one group at a time,
# into its own output/<case>/ folder, and shows you the piece of output
# that proves what just happened. Pauses between sections so you can read
# before moving on. Safe to re-run any time -- it always starts by wiping
# its own output/.
#
# Fixtures here are small, purpose-built synthetic files (example.PDOS.xml
# + a matching example.bands/example.EIG "same calculation" trio), not a
# real SIESTA PDOS.xml: a real one is easily tens of MB (broadened DOS on
# a fine energy grid, one <orbital> block per basis function per atom),
# far too heavy for a lightweight example. These were built by hand with
# known VBM/CBM values so every number below can be checked directly.

set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

# --view opens an interactive matplotlib preview at the very end of a run
# -- MPLBACKEND=Agg makes that a no-op instead of blocking on a GUI
# window, same convention test.sh and example_3.1.sh already use.
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
echo " What a PDOS.xml file is, and what stb-dos does with it"
echo "=================================================================="
cat <<'EOF'
A SIESTA <label>.PDOS.xml file holds the projected density of states: for
every orbital of every atom, the DOS contribution (already broadened --
smoothed with a Gaussian/Lorentzian, SIESTA's own PDOS.info setting) on a
shared energy grid, plus the calculation's Fermi energy and nspin.

stb-dos reads this and writes:
  - dos_total.dat        -- the DOS summed over every atom/orbital
  - dos_per_atom/*.dat    -- one file per atom
  - dos_per_species/*.dat -- one file per chemical species
each with one column per orbital angular momentum (s, p, d, f -- or, with
--projection ml, the individual m sub-orbitals: px, py, pz, ...), split
into spin-up/spin-down columns for a spin-polarized (nspin=2) calculation.

--shift controls which energy sits at 0 eV in the written files:
  fermi  -- the Fermi energy (default)
  vbm/cbm -- the Valence Band Maximum / Conduction Band Minimum (see below)
  a number -- any manual value
EOF
pause

echo "=================================================================="
echo " output/basic-run/  --  the three output types, --type/--projection"
echo "=================================================================="
cat <<'EOF'
example.PDOS.xml: 1 atom (Si), 1 orbital (s), 13 energy points from -3 to
+3 eV, Fermi energy 0.0 eV. Deliberately built with a DOS "gap" shape:
5.0 (peak) far from E_F, dropping to 0.01 (near-zero) between -1 and +1 eV
-- used again below for the --estimate-from-dos case.
EOF
mkdir -p "$OUT/basic-run"
cp example.PDOS.xml "$OUT/basic-run/"
echo
echo "\$ stb-dos example.PDOS.xml --shift fermi --no-intro"
(cd "$OUT/basic-run" && stb-dos example.PDOS.xml --shift fermi --no-intro > console.log 2>&1)
cat "$OUT/basic-run/console.log"
echo
ls "$OUT/basic-run" | grep -v console
pause

echo "=================================================================="
echo " VBM/CBM and direct vs. indirect gaps -- why the DOS alone isn't enough"
echo "=================================================================="
cat <<'EOF'
A PDOS.xml is a broadened, continuous curve -- there is no discrete list
of eigenvalues in it, so there is no exact VBM/CBM to read off directly
(the same way stb-bands reads one from discrete .bands/.EIG eigenvalues).
--shift vbm/cbm instead reuses stb-bands' own VBM/CBM machinery
(core/siesta_bands.py) on a companion file next to the PDOS.xml, found
from a single --label (or auto-derived from the PDOS.xml's own filename):

    1. <label>.bands  (a k-PATH, the same file stb-bands itself defaults
       to)                                                       -- tried first
    2. <label>.EIG    (the full SCF k-MESH, sisl-read)             -- fallback
    3. an approximation from the DOS itself (--estimate-from-dos)  -- last resort

example.bands and example.EIG below are a purpose-built "same calculation"
trio alongside example.PDOS.xml (same Fermi energy, 0.0 eV) with
known-by-hand VBM/CBM at each level of the hierarchy -- watch the exact
same number stb-bands would report come out of stb-dos here.
EOF
mkdir -p "$OUT/vbm-hierarchy"
cp example.PDOS.xml example.bands example.EIG "$OUT/vbm-hierarchy/"
echo
echo "\$ stb-dos example.PDOS.xml --shift vbm --no-intro   # example.bands present -> used first"
(cd "$OUT/vbm-hierarchy" && stb-dos example.PDOS.xml --shift vbm --no-intro > console_bands.log 2>&1)
grep "Using VBM shift" "$OUT/vbm-hierarchy/console_bands.log"
echo
echo "Removing example.bands -- falls back to the k-mesh example.EIG:"
mv "$OUT/vbm-hierarchy/example.bands" "$OUT/vbm-hierarchy/example.bands.aside"
echo "\$ stb-dos example.PDOS.xml --shift vbm --no-intro   # no .bands now -> .EIG used"
(cd "$OUT/vbm-hierarchy" && stb-dos example.PDOS.xml --shift vbm --no-intro > console_eig.log 2>&1)
grep "Using VBM shift" "$OUT/vbm-hierarchy/console_eig.log"
echo
echo "Note the different value (-1.000000 eV vs -0.500000 eV): the .bands"
echo "k-PATH and the .EIG k-MESH are, in general, different samplings of"
echo "the same Brillouin zone -- exactly the same distinction stb-bands'"
echo "own --eig-file mesh-vs-line comparison is built around."
pause

echo "=================================================================="
echo " output/vbm-hierarchy/  --  neither .bands nor .EIG: error, or --estimate-from-dos"
echo "=================================================================="
cat <<'EOF'
With neither a .bands nor a .EIG file available, an exact VBM/CBM simply
cannot be computed -- stb-dos refuses by default rather than silently
guessing a number with different physics behind it:
EOF
mv "$OUT/vbm-hierarchy/example.EIG" "$OUT/vbm-hierarchy/example.EIG.aside"
echo
echo "\$ stb-dos example.PDOS.xml --shift vbm --no-intro   # neither file present"
set +e
(cd "$OUT/vbm-hierarchy" && stb-dos example.PDOS.xml --shift vbm --no-intro > console_none.log 2>&1)
EXIT_NONE=$?
set -e
echo "Exit code: $EXIT_NONE"
grep "^Error:" "$OUT/vbm-hierarchy/console_none.log"
echo
cat <<'EOF'
--estimate-from-dos opts into a heuristic fallback instead: walk outward
from the Fermi energy on each side until the (summed) total DOS first
rises above --dos-threshold-frac * its own peak value -- the point where
the curve leaves the near-zero gap region. This is explicitly an
APPROXIMATION, blurred by the PDOS's own energy broadening, printed as an
explicit warning rather than presented as an exact result:
EOF
echo
echo "\$ stb-dos example.PDOS.xml --shift vbm --estimate-from-dos --no-intro"
(cd "$OUT/vbm-hierarchy" && stb-dos example.PDOS.xml --shift vbm --estimate-from-dos \
    --no-intro > console_estimate.log 2>&1)
grep "\[WARNING\]" "$OUT/vbm-hierarchy/console_estimate.log"
echo
echo "Compare: -1.500000 eV (DOS estimate) vs. -1.000000 eV (real .bands VBM"
echo "above) vs. -0.500000 eV (real .EIG VBM) -- all in the same right"
echo "ballpark (the DOS genuinely does drop to near-zero around there), but"
echo "not identical, which is exactly why this fallback needs an explicit"
echo "opt-in flag instead of being used silently."
mv "$OUT/vbm-hierarchy/example.bands.aside" "$OUT/vbm-hierarchy/example.bands"
mv "$OUT/vbm-hierarchy/example.EIG.aside" "$OUT/vbm-hierarchy/example.EIG"
pause

echo "=================================================================="
echo " --label: one label covers the PDOS.xml AND its .bands/.EIG"
echo "=================================================================="
cat <<'EOF'
--label lets you skip typing the PDOS.xml path at all: it resolves to
<label>.PDOS.xml, and (for --shift vbm/cbm) the exact same <label>.bands/
<label>.EIG lookup above uses that same label -- one name covers all
three files, since they all follow SIESTA's own <label>.<ext> convention.
EOF
mkdir -p "$OUT/label-shorthand"
cp example.PDOS.xml "$OUT/label-shorthand/example.PDOS.xml"
cp example.bands "$OUT/label-shorthand/example.bands"
echo
echo "\$ stb-dos --label example --shift vbm --no-intro"
(cd "$OUT/label-shorthand" && stb-dos --label example --shift vbm --no-intro > console.log 2>&1)
grep "Using VBM shift" "$OUT/label-shorthand/console.log"
ls "$OUT/label-shorthand" | grep -v console
pause

echo "=================================================================="
echo " output/full-report/  --  --save-report/--save-gnuplot (both off by default)"
echo "=================================================================="
cat <<'EOF'
Every run always prints the numbered [0]...[5] report to the console.
--save-report additionally persists it to stb_dos_report.txt, and
--save-gnuplot writes ONE .gplot script per output category (total,
atom, species) -- one overlay plot per category instead of one file per
individual atom/species -- both off by default, so a plain run only ever
writes the .dat files and references.bib:
EOF
mkdir -p "$OUT/full-report"
cp example.PDOS.xml "$OUT/full-report/"
echo
echo "\$ stb-dos example.PDOS.xml --shift fermi --no-intro   # default: no report file, no .gplot"
(cd "$OUT/full-report" && stb-dos example.PDOS.xml --shift fermi --no-intro > console_default.log 2>&1)
ls "$OUT/full-report/" | grep -v "\.PDOS\.xml$\|console"
echo "(no stb_dos_report.txt, no .gplot files -- only the .dat files + references.bib)"
echo
echo "\$ stb-dos example.PDOS.xml --shift fermi --save-report --save-gnuplot --no-intro"
(cd "$OUT/full-report" && stb-dos example.PDOS.xml --shift fermi --save-report --save-gnuplot \
    --no-intro > console_saved.log 2>&1)
echo "Report sections written to stb_dos_report.txt:"
grep -E "^\[[0-9]+\]" "$OUT/full-report/stb_dos_report.txt"
echo
echo "Gnuplot scripts, one per output category (dos_total.gplot uses columnheader() per"
echo "orbital column; dos_per_atom/species.gplot overlay one summed curve per atom/species):"
find "$OUT/full-report" -name "*.gplot" | sed "s#$OUT/##" | sort
echo
echo "references.bib -- SIESTA (every stb-dos run is post-processing a SIESTA output):"
grep "^@" "$OUT/full-report/references.bib"
pause

echo "=================================================================="
echo " Two ways to run it"
echo "=================================================================="
cat <<'EOF'
A -- direct CLI:
  stb-dos --label example --shift vbm

B -- interactive stb-suite menu:
  $ stb-suite
  Select an option (0-6, or a tool code like 4.1.2): 3.2

Both paths call the exact same underlying tool -- proven directly below.
EOF
TMP="$(mktemp -d)"
cp example.PDOS.xml example.bands example.EIG "$TMP/"
echo
echo "\$ printf '3.2\\nexample\\n\\n2\\n\\n\\n\\n\\n\\n0\\n' | stb-suite     # label 'example', shift menu choice 2 = VBM"
(cd "$TMP" && printf '3.2\nexample\n\n2\n\n\n\n\n\n0\n' | stb-suite > session.log 2>&1) || true
CLI_LINE=$(grep "Using VBM shift" "$OUT/vbm-hierarchy/console_bands.log")
MENU_LINE=$(grep "Using VBM shift" "$TMP/session.log")
if [ "$CLI_LINE" = "$MENU_LINE" ]; then
    echo "Confirmed: identical VBM shift from the CLI and the interactive menu."
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
  basic-run/         vbm-hierarchy/
  label-shorthand/   full-report/

Recap of what this walkthrough covered:
  - what a PDOS.xml is, and the three output types stb-dos writes
    (dos_total.dat, dos_per_atom/, dos_per_species/)
  - --shift vbm/cbm's .bands -> .EIG -> DOS-estimate hierarchy, on a
    purpose-built "same calculation" trio with known-by-hand values
  - why the DOS-estimate fallback needs an explicit --estimate-from-dos
    opt-in, and prints an explicit approximation warning
  - --label as shorthand for the PDOS.xml AND its .bands/.EIG lookup
  - the numbered [0]...[5] report, --save-report, --save-gnuplot (one
    .gplot per output category: total/atom/species), and references.bib
  - CLI and the interactive stb-suite menu building the same command
    (--view is offered in the menu too, as its own y/N prompt)

Not exercised by this script (needs a display): the interactive
matplotlib preview (--view, one figure per category) -- try it yourself
without MPLBACKEND=Agg:
  stb-dos example.PDOS.xml --shift fermi --view

As a next step, try on your own with a real SIESTA calculation:
  stb-dos --label my_calc --shift vbm --save-report
  stb-dos --label my_calc --shift cbm --estimate-from-dos   # no .bands/.EIG
EOF
