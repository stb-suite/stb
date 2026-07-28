#!/bin/bash
# Guided example: stb-spintexture (code 3.17 in the stb-suite menu)
#
# Not an automated test (see test/3-analysis/17-spintexture/test.sh for
# that) -- a commented walk-through: it runs real commands, one group at
# a time, into its own output/<case>/ folder, and shows you the piece of
# output that proves what just happened. Pauses between sections so you
# can read before moving on. Safe to re-run any time -- it always starts
# by wiping its own output/.
#
# IsolatedO.HSX/IsolatedO.selected.WFSX are a REAL, non-collinear SIESTA
# calculation (a single oxygen atom in a large vacuum box, Gamma-point
# only). calc.fdf is the same structure's geometry-only input, used below
# to demonstrate a real, live-verified accuracy gap.

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
echo " What a spin texture is"
echo "=================================================================="
cat <<'EOF'
A non-collinear/SOC calculation represents each electron as a full
2-component spinor, whose spin direction can point anywhere in 3D --
described by the Pauli expectation values <Sx>, <Sy>, <Sz>. For a
properly normalized spinor:

  |S| = sqrt(<Sx>^2 + <Sy>^2 + <Sz>^2)  <=  1   (always)

Needs the real overlap matrix (spin_moment() calls Sk()) -- the same
.HSX-accurate/.fdf-approximate tradeoff stb-fatbands/stb-ipr have. Unlike
those two, no companion .bands file is needed at all: the x-axis is a
plain k-INDEX, not a physical k-path arc length.

Every run prints a numbered [0]...[6] report. --save-report persists it;
--save-gnuplot writes spintexture_S{x,y,z}.dat + a real .gplot script
(off by default -- this tool used to write both unconditionally); --view
shows the same as a matplotlib preview (off by default -- this tool
previously always showed it, with no way to skip at all).
EOF
pause


echo "=================================================================="
echo " output/basic-run/  --  the numbered report, the |S|<=1 normalization check"
echo "=================================================================="
mkdir -p "$OUT/basic-run"
cp IsolatedO.HSX IsolatedO.selected.WFSX "$OUT/basic-run/"
cat <<'EOF'
This fixture converges (without an explicit initial spin canting) to a
physically sensible, strongly non-zero Sz (~+-1.0, expected for this
open-shell atom) with Sx/Sy ~0 (numerical noise) -- validating the
numerical pipeline, though not demonstrating a genuinely CANTED texture
(that needs real spin-orbit coupling on a heavier element). Watch [3]
report the |S| normalization check:
EOF
echo "\$ stb-spintexture --wfsx IsolatedO.selected.WFSX --hsx-file IsolatedO.HSX --no-intro"
(cd "$OUT/basic-run" && stb-spintexture --wfsx IsolatedO.selected.WFSX --hsx-file IsolatedO.HSX \
    --no-intro | tee console.log | sed -n '/\[3\] SPIN TEXTURE/,/\[4\] OUTPUT/p')
pause


echo "=================================================================="
echo " A real, live-verified reason the .HSX accuracy warning matters"
echo "=================================================================="
cat <<'EOF'
Without a real .HSX, spin_moment() falls back to an implicit-orthogonal
approximation -- watch this go from merely "less accurate" to outright
UNPHYSICAL (|S| exceeding its hard 1.0 bound by 18x!), and the new
normalization check catching it automatically:
EOF
mkdir -p "$OUT/accuracy"
cp IsolatedO.selected.WFSX calc.fdf "$OUT/accuracy/"
echo "\$ stb-spintexture --wfsx IsolatedO.selected.WFSX --geometry-file calc.fdf --no-intro"
(cd "$OUT/accuracy" && stb-spintexture --wfsx IsolatedO.selected.WFSX --geometry-file calc.fdf \
    --no-intro | sed -n '/\[3\] SPIN TEXTURE/,/\[4\] OUTPUT/p')
pause


echo "=================================================================="
echo " A real bug this tool inherited from stb-wfdensity/stb-sts/stb-coop/stb-effmass: --label + --hsx-file"
echo "=================================================================="
cat <<'EOF'
This combination used to be rejected outright. Simulating a renamed .HSX
(a common real-world mismatch, e.g. SystemLabel not matching the real
input filename):
EOF
mkdir -p "$OUT/label-plus-hsx"
cp IsolatedO.selected.WFSX "$OUT/label-plus-hsx/IsolatedO.selected.WFSX"
cp IsolatedO.HSX "$OUT/label-plus-hsx/calc.HSX"
echo "\$ stb-spintexture --label IsolatedO --hsx-file calc.HSX --no-intro"
(cd "$OUT/label-plus-hsx" && stb-spintexture --label IsolatedO --hsx-file calc.HSX \
    --no-intro > console.log 2>&1)
grep "WFSX file\|Hamiltonian src" "$OUT/label-plus-hsx/console.log"
pause


echo "=================================================================="
echo " output/fermi-shift/  --  --shift fermi/vbm/cbm's decoupled-from-label hierarchy"
echo "=================================================================="
cat <<'EOF'
The same priority-ordered hierarchy stb-effmass has: --fermi (explicit) >
--bands-file > --fermi-file > an auto-detected .out log in the current
directory. (No real .out was saved alongside this fixture's .WFSX/.HSX,
so calc.out below is a representative demonstration value.)
EOF
mkdir -p "$OUT/fermi-shift"
cp IsolatedO.HSX IsolatedO.selected.WFSX "$OUT/fermi-shift/"
cat > "$OUT/fermi-shift/calc.out" << 'EOF'
Some SIESTA log noise before the summary
siesta:         Fermi = -12.0
More noise after
EOF
echo "\$ stb-spintexture --wfsx IsolatedO.selected.WFSX --hsx-file IsolatedO.HSX --shift fermi --no-intro   # auto-detects calc.out"
(cd "$OUT/fermi-shift" && stb-spintexture --wfsx IsolatedO.selected.WFSX --hsx-file IsolatedO.HSX \
    --shift fermi --no-intro > console.log 2>&1)
grep "Shift mode" "$OUT/fermi-shift/console.log"
pause


echo "=================================================================="
echo " output/full-report/  --  --save-report / --save-gnuplot (off by default)"
echo "=================================================================="
cat <<'EOF'
Every run always prints the numbered [0]...[6] report to the console.
--save-report additionally persists it to stb_spintexture_report.txt,
and --save-gnuplot writes spintexture_S{x,y,z}.dat + a real .gplot
script -- both off by default, so a plain run only ever writes
references.bib:
EOF
mkdir -p "$OUT/full-report"
cp IsolatedO.HSX IsolatedO.selected.WFSX "$OUT/full-report/"
echo "\$ stb-spintexture --wfsx IsolatedO.selected.WFSX --hsx-file IsolatedO.HSX --no-intro   # default"
(cd "$OUT/full-report" && stb-spintexture --wfsx IsolatedO.selected.WFSX --hsx-file IsolatedO.HSX \
    --no-intro > console_default.log 2>&1)
(cd "$OUT/full-report" && ls)
echo "(only references.bib -- no spintexture_S{x,y,z}.dat/.gplot, no stb_spintexture_report.txt)"
echo
echo "\$ stb-spintexture --wfsx IsolatedO.selected.WFSX --hsx-file IsolatedO.HSX \\"
echo "      --save-report --save-gnuplot --no-intro"
(cd "$OUT/full-report" && stb-spintexture --wfsx IsolatedO.selected.WFSX --hsx-file IsolatedO.HSX \
    --save-report --save-gnuplot --no-intro > console_saved.log 2>&1)
echo "Report sections written to stb_spintexture_report.txt:"
grep -o "^\[[0-9]\] [A-Za-z& ]*" "$OUT/full-report/stb_spintexture_report.txt" | sort -u
echo
if command -v gnuplot >/dev/null 2>&1; then
    (cd "$OUT/full-report" && gnuplot spintexture.gplot 2>/dev/null)
    echo "(rendered spintexture.pdf with the real, installed gnuplot)"
fi
pause


echo "=================================================================="
echo " Two ways to run it"
echo "=================================================================="
cat <<EOF
A -- direct CLI:
  stb-spintexture --label IsolatedO --shift fermi --fermi -12.0

B -- interactive stb-suite menu:
  \$ stb-suite
  Select an option (0-6, or a tool code like 4.1.2): 3.17

Both paths call the exact same underlying tool -- proven directly below.
The menu asks for the label and the .HSX path SEPARATELY (default
suggested if <label>.HSX exists), gained a Fermi-source submenu for
fermi/vbm/cbm, and save-report/save-gnuplot/view prompts.
EOF
TMP="$(mktemp -d)"
cp IsolatedO.HSX IsolatedO.selected.WFSX "$TMP/"
echo
echo "\$ stb-spintexture --label IsolatedO --shift fermi --fermi -12.0 --no-intro"
(cd "$TMP" && stb-spintexture --label IsolatedO --shift fermi --fermi -12.0 \
    --no-intro > cli.log 2>&1)
echo
echo "\$ printf '3.17\\nIsolatedO\\n\\n4\\n1\\n-12.0\\n\\nn\\nn\\nn\\n\\n0\\n' | stb-suite"
(cd "$TMP" && printf '3.17\nIsolatedO\n\n4\n1\n-12.0\n\nn\nn\nn\n\n0\n' | timeout 60 stb-suite > session.log 2>&1) || true
CLI_LINE=$(grep "Most spin-polarized state" "$TMP/cli.log")
MENU_LINE=$(grep "Most spin-polarized state" "$TMP/session.log")
if [ "$CLI_LINE" = "$MENU_LINE" ]; then
    echo "Confirmed: identical result from the CLI and the interactive menu."
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
  basic-run/       accuracy/       label-plus-hsx/
  fermi-shift/     full-report/

Each has references.bib; full-report/ additionally has
spintexture_S{x,y,z}.dat/.gplot/.pdf and stb_spintexture_report.txt.

Recap of what this walkthrough covered:
  - spin texture theory: the Pauli expectation values <Sx>,<Sy>,<Sz>,
    and the |S| <= 1 physical bound for a normalized spinor
  - a real, live-verified finding: without a real .HSX, spin_moment()
    doesn't just get slightly less accurate -- it can give values that
    violate the hard |S| <= 1 bound by 18x, now caught by a new
    normalization check
  - a real bug: --label + --hsx-file used to be rejected outright; fixed
    the same way as stb-wfdensity/stb-sts/stb-coop/stb-effmass
  - --shift fermi/vbm/cbm's Fermi-source hierarchy, shared with
    stb-effmass
  - the numbered [0]...[6] report, --save-report, --save-gnuplot (now
    opt-in, previously unconditional), --view (now opt-in, previously
    always shown with no way to skip)
  - CLI and the interactive stb-suite menu building the same command

As a next step, try on your own with a real SIESTA calculation:
  stb-spintexture --label my_calc --shift fermi --fermi-file my_calc.out --save-gnuplot --view
EOF
