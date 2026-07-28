#!/bin/bash
# Guided example: stb-sts (code 3.13 in the stb-suite menu)
#
# Not an automated test (see test/3-analysis/13-sts/test.sh for that) -- a
# commented walk-through: it runs real commands, one group at a time, into
# its own output/<case>/ folder, and shows you the piece of output that
# proves what just happened. Pauses between sections so you can read
# before moving on. Safe to re-run any time -- it always starts by wiping
# its own output/.
#
# calc.fdf/structure.fdf/C.ion(.xml)/Graphene.selected.WFSX/Graphene.XV are
# a REAL, non-polarized SIESTA calculation of monolayer graphene (2-atom
# primitive hexagonal cell, 20 Ang vacuum along Z), with a genuine full-BZ
# WFSX (9 explicit k-points via WriteWaveFunctions T + %block
# WaveFuncKPoints -- NOT a band-path WFSX like stb-fatbands needs).
# SystemLabel is "Graphene", but there is no Graphene.fdf anywhere -- the
# real input file is calc.fdf, the same mismatch already exercised in the
# 3.12-wfdensity example.

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
echo " What stb-sts computes, and why it needs a full-BZ .WFSX"
echo "=================================================================="
cat <<'EOF'
In the Tersoff-Hamann tunneling approximation, dI/dV at a fixed tip point
r and bias energy E is proportional to the local density of states there:

  dI/dV(r, E) ~ LDOS(r, E) = sum_(n,k) |psi_n(k, r)|^2 * delta(E - eps_n(k))

stb-stm already uses this exact proxy for a SPATIAL image (fix the energy
window, scan r). stb-sts does the opposite reduction: fix r at one point,
scan E -- "PDOS projected onto a single real-space point instead of onto
atomic orbitals". Every (k, band, spin) state is Gaussian-broadened around
its own eigenvalue and weighted by its own Brillouin-zone weight, then
summed into one continuous curve.

This needs a full-BZ-sampled .WFSX (WriteWaveFunctions T + an explicit
%block WaveFuncKPoints listing a real k-mesh -- SIESTA names this
<label>.selected.WFSX) -- NOT the band-path .WFSX stb-fatbands uses, since
this is a DOS-like, Brillouin-zone-integrated quantity, and needs the real
per-orbital radial basis (.fdf + .ion/.ion.xml files), same requirement as
stb-wfdensity -- never a bare .XV/.HSX.

Every run prints a numbered [0]...[6] report. --save-report persists it;
--save-gnuplot writes sts.dat + a real .gplot script (off by default --
this tool used to write the .dat unconditionally but never actually wrote
a .gplot at all); --view shows a matplotlib preview (replaces the old
--no-plot, which used to be on by default).
EOF
pause


echo "=================================================================="
echo " output/basic-run/  --  --xy/--height, the numbered report"
echo "=================================================================="
mkdir -p "$OUT/basic-run"
cp calc.fdf structure.fdf C.ion C.ion.xml Graphene.selected.WFSX "$OUT/basic-run/"
cat <<'EOF'
A tip at (0, 0), 1.5 Ang above the graphene sheet (well within the carbon
DZP basis's own ~2.58 Ang cutoff radius) -- watch [2] TIP POSITION report
where the topmost surface atom actually is, and [3] STS CURVE report the
peak signal found:
EOF
echo "\$ stb-sts --label Graphene --geometry-file calc.fdf --xy 0 0 --height 1.5 --erange -5 5 --sigma 50 --no-intro"
(cd "$OUT/basic-run" && stb-sts --label Graphene --geometry-file calc.fdf \
    --xy 0 0 --height 1.5 --erange -5 5 --sigma 50 --no-intro | tee console.log | \
    sed -n '/\[2\] TIP POSITION/,/\[4\] OUTPUT/p')
pause


echo "=================================================================="
echo " A real bug this tool inherited from stb-wfdensity: --label + --geometry-file"
echo "=================================================================="
cat <<'EOF'
SystemLabel is "Graphene" here, but the real input file is calc.fdf -- no
Graphene.fdf exists anywhere. Watch --label alone fail cleanly, then
succeed once --geometry-file is added -- --label still auto-detects the
.WFSX either way, only the geometry source changes. This combination used
to be rejected outright ("cannot be combined") -- fixed the same way it
was fixed for stb-wfdensity.
EOF
mkdir -p "$OUT/label-plus-geometry"
cp calc.fdf structure.fdf C.ion C.ion.xml Graphene.selected.WFSX Graphene.XV "$OUT/label-plus-geometry/"
echo "\$ stb-sts --label Graphene --xy 0 0 --height 1.5 --erange -5 5 --sigma 50 --no-intro   # no --geometry-file: fails"
(cd "$OUT/label-plus-geometry" && stb-sts --label Graphene \
    --xy 0 0 --height 1.5 --erange -5 5 --sigma 50 --no-intro > console_fail.log 2>&1) || true
grep -A1 "Resolving geometry" "$OUT/label-plus-geometry/console_fail.log"
echo
echo "\$ stb-sts --label Graphene --geometry-file calc.fdf --xy 0 0 --height 1.5 --erange -5 5 --sigma 50 --no-intro"
(cd "$OUT/label-plus-geometry" && stb-sts --label Graphene --geometry-file calc.fdf \
    --xy 0 0 --height 1.5 --erange -5 5 --sigma 50 --no-intro > console_ok.log 2>&1)
grep "WFSX file\|Geometry source" "$OUT/label-plus-geometry/console_ok.log"
pause


echo "=================================================================="
echo " output/cutoff-radius/  --  SIESTA's confined PAO basis is exactly zero past its cutoff"
echo "=================================================================="
cat <<'EOF'
Unlike a plane-wave or Gaussian basis, SIESTA's numerical atomic orbitals
have a FINITE, hard cutoff radius -- for this fixture's carbon DZP basis,
verified directly on the geometry: max(orbital.R) = 2.576 Ang. A tip only
0.4 Ang beyond that gives an EXACTLY zero (not just small) signal, with no
error -- only a warning. Watch the peak signal collapse to 0.0 between
--height 1.5 and --height 3.0:
EOF
mkdir -p "$OUT/cutoff-radius"
cp calc.fdf structure.fdf C.ion C.ion.xml Graphene.selected.WFSX "$OUT/cutoff-radius/"
for h in 1.5 3.0; do
    echo
    echo "\$ stb-sts --label Graphene --geometry-file calc.fdf --xy 0 0 --height $h --erange -5 5 --sigma 50 --no-intro"
    (cd "$OUT/cutoff-radius" && stb-sts --label Graphene --geometry-file calc.fdf \
        --xy 0 0 --height "$h" --erange -5 5 --sigma 50 --no-intro > "console_h${h}.log" 2>&1) || true
    grep "Peak dI/dV\|WARNING" "$OUT/cutoff-radius/console_h${h}.log"
done
pause


echo "=================================================================="
echo " A real, verified physics finding: the vacuum axis is still periodic"
echo "=================================================================="
cat <<'EOF'
This cell's Z length is 20 Ang, atoms sit at Z=10 Ang -- so the periodic
IMAGE of the same atomic plane recurs at Z=10+20=30 Ang, i.e. --height=20.
A --height comparable to or beyond roughly half the vacuum thickness does
NOT probe "farther into empty space" -- it silently wraps toward the NEXT
periodic copy of the same surface. Watch the signal reappear at --height
19-20 after being genuinely zero (true vacuum, beyond the PAO cutoff
either way) at --height 8-15:
EOF
mkdir -p "$OUT/periodic-wraparound"
cp calc.fdf structure.fdf C.ion C.ion.xml Graphene.selected.WFSX "$OUT/periodic-wraparound/"
for h in 10.0 19.0 20.0; do
    echo
    echo "\$ stb-sts --label Graphene --geometry-file calc.fdf --xy 0 0 --height $h --erange -5 5 --sigma 50 --no-intro"
    (cd "$OUT/periodic-wraparound" && stb-sts --label Graphene --geometry-file calc.fdf \
        --xy 0 0 --height "$h" --erange -5 5 --sigma 50 --no-intro > "console_h${h}.log" 2>&1) || true
    grep "Absolute position\|Peak dI/dV" "$OUT/periodic-wraparound/console_h${h}.log"
done
echo
echo "(--height 20 lands almost exactly back on the periodic image of the atomic plane --"
echo " both the real tip position and its far-away 'copy' give strong signal, for the same reason.)"
pause


echo "=================================================================="
echo " output/fermi-shift/  --  --shift fermi's decoupled-from-label hierarchy"
echo "=================================================================="
cat <<'EOF'
--shift fermi used to accept ONLY an explicit --fermi value. It now has
the same priority-ordered hierarchy stb-wfdensity's --band vbm/cbm has:
--fermi (explicit) > --bands-file > --fermi-file > an auto-detected .out
log in the current directory -- decoupled from --label, since many real
SIESTA jobs redirect stdout to a generic name instead of <label>.out.
(No real .out was saved alongside this fixture's .WFSX, so calc.out below
is a representative demonstration value, not one re-derived from this
exact run.)
EOF
mkdir -p "$OUT/fermi-shift"
cp calc.fdf structure.fdf C.ion C.ion.xml Graphene.selected.WFSX "$OUT/fermi-shift/"
cat > "$OUT/fermi-shift/calc.out" << 'EOF'
Some SIESTA log noise before the summary
siesta:         Fermi = -4.461866
More noise after
EOF
echo "\$ stb-sts --label Graphene --geometry-file calc.fdf --xy 0 0 --height 1.5 --erange -5 5 --sigma 50 \\"
echo "      --shift fermi --fermi -4.461866 --no-intro   # explicit value"
(cd "$OUT/fermi-shift" && stb-sts --label Graphene --geometry-file calc.fdf \
    --xy 0 0 --height 1.5 --erange -5 5 --sigma 50 --shift fermi --fermi -4.461866 \
    --no-intro > console_explicit.log 2>&1)
grep "Energy shift" "$OUT/fermi-shift/console_explicit.log"
echo
echo "\$ stb-sts --label Graphene --geometry-file calc.fdf --xy 0 0 --height 1.5 --erange -5 5 --sigma 50 \\"
echo "      --shift fermi --no-intro   # no --fermi/--fermi-file: auto-detects calc.out"
(cd "$OUT/fermi-shift" && stb-sts --label Graphene --geometry-file calc.fdf \
    --xy 0 0 --height 1.5 --erange -5 5 --sigma 50 --shift fermi \
    --no-intro > console_auto.log 2>&1)
grep "Energy shift" "$OUT/fermi-shift/console_auto.log"
pause


echo "=================================================================="
echo " output/full-report/  --  --save-report / --save-gnuplot (off by default)"
echo "=================================================================="
cat <<'EOF'
Every run always prints the numbered [0]...[6] report to the console.
--save-report additionally persists it to stb_sts_report.txt, and
--save-gnuplot writes sts.dat + a real .gplot script -- both off by
default, so a plain run only ever writes sts.dat + references.bib:
EOF
mkdir -p "$OUT/full-report"
cp calc.fdf structure.fdf C.ion C.ion.xml Graphene.selected.WFSX "$OUT/full-report/"
echo "\$ stb-sts --label Graphene --geometry-file calc.fdf --xy 0 0 --height 1.5 --erange -5 5 --sigma 50 --no-intro   # default"
(cd "$OUT/full-report" && stb-sts --label Graphene --geometry-file calc.fdf \
    --xy 0 0 --height 1.5 --erange -5 5 --sigma 50 --no-intro > console_default.log 2>&1)
(cd "$OUT/full-report" && ls)
echo "(only sts.dat + references.bib -- no .gplot, no stb_sts_report.txt)"
echo
echo "\$ stb-sts --label Graphene --geometry-file calc.fdf --xy 0 0 --height 1.5 --erange -5 5 --sigma 50 \\"
echo "      --save-report --save-gnuplot --no-intro"
(cd "$OUT/full-report" && stb-sts --label Graphene --geometry-file calc.fdf \
    --xy 0 0 --height 1.5 --erange -5 5 --sigma 50 --save-report --save-gnuplot \
    --no-intro > console_saved.log 2>&1)
echo "Report sections written to stb_sts_report.txt:"
grep -o '\[[0-9]\] [A-Z& ]*' "$OUT/full-report/stb_sts_report.txt" | sort -u
echo
if command -v gnuplot >/dev/null 2>&1; then
    (cd "$OUT/full-report" && gnuplot sts.gplot)
    echo "(rendered sts.pdf with the real, installed gnuplot)"
fi
pause


echo "=================================================================="
echo " Two ways to run it"
echo "=================================================================="
cat <<EOF
A -- direct CLI:
  stb-sts --label Graphene --geometry-file calc.fdf --xy 0 0 --height 1.5 --erange -5 5 --sigma 50

B -- interactive stb-suite menu:
  \$ stb-suite
  Select an option (0-6, or a tool code like 4.1.2): 3.13

Both paths call the exact same underlying tool -- proven directly below.
The menu asks for the label and the .fdf path SEPARATELY, and offers an
energy-shift submenu (with the same Fermi-source options shown above) for
--shift fermi.
EOF
TMP="$(mktemp -d)"
cp calc.fdf structure.fdf C.ion C.ion.xml Graphene.selected.WFSX "$TMP/"
cp calc.fdf "$TMP/Graphene.fdf"
echo
echo "\$ printf '3.13\\nGraphene\\n\\n1\\n0\\n0\\n1.5\\n-5\\n5\\n50\\n\\n\\nn\\nn\\nn\\n\\n0\\n' | stb-suite"
(cd "$TMP" && printf '3.13\nGraphene\n\n1\n0\n0\n1.5\n-5\n5\n50\n\n\nn\nn\nn\n\n0\n' | timeout 60 stb-suite > session.log 2>&1) || true
CLI_LINE=$(grep "Peak dI/dV" "$OUT/basic-run/console.log")
MENU_LINE=$(grep "Peak dI/dV" "$TMP/session.log")
if [ "$CLI_LINE" = "$MENU_LINE" ]; then
    echo "Confirmed: identical peak dI/dV from the CLI and the interactive menu."
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
Six self-contained folders were generated under output/:
  basic-run/            label-plus-geometry/   cutoff-radius/
  periodic-wraparound/  fermi-shift/           full-report/

Each has references.bib and sts.dat; full-report/ additionally has
sts.gplot/sts.pdf and stb_sts_report.txt (only from its --save-report/
--save-gnuplot run).

Recap of what this walkthrough covered:
  - the Tersoff-Hamann LDOS proxy, extended from stb-stm's spatial map to
    a fixed-point energy curve, and where stb-sts sits among the suite's
    other real-space/energy tools
  - why a full-BZ .WFSX (WriteWaveFunctions T + WaveFuncKPoints) is
    required, not a band-path one like stb-fatbands uses
  - a real bug: --label + --geometry-file used to be rejected outright;
    fixed the same way it was for stb-wfdensity, with this exact fixture
    (SystemLabel != real fdf filename) as the case that caught it
  - SIESTA's confined PAO basis has an exact, finite cutoff radius
    (verified: 2.576 Ang for this carbon DZP basis) -- --height 3.0
    already gives an exactly-zero curve
  - a real, verified physics finding: the "vacuum" axis is still
    periodic, so a --height comparable to the vacuum thickness silently
    samples the NEXT periodic image of the surface, not true vacuum decay
  - --shift fermi's new priority-ordered Fermi-source hierarchy, shared
    with stb-wfdensity via core.siesta_bands.resolve_fermi_energy_hierarchy
  - the numbered [0]...[6] report, --save-report, --save-gnuplot (now
    with a REAL .gplot script, previously missing entirely), --view
    (replacing the old, on-by-default --no-plot)
  - CLI and the interactive stb-suite menu building the same command

As a next step, try on your own with a real SIESTA calculation:
  stb-sts --label my_calc --xy 0 0 --height 2.0 --erange -3 3 --sigma 50 --view
  stb-sts --label my_calc --point 1.2 0.7 12.0 --erange -5 5 --fwhm 100
EOF
