#!/bin/bash
# Guided example: stb-wfdensity (code 3.12 in the stb-suite menu)
#
# Not an automated test (see test/3-analysis/12-wfdensity/test.sh for that,
# which uses a separate, non-spin-polarized bulk Sn3O4 fixture) -- a
# commented walk-through: it runs real commands, one group at a time, into
# its own output/<case>/ folder, and shows you the piece of output that
# proves what just happened. Pauses between sections so you can read
# before moving on. Safe to re-run any time -- it always starts by wiping
# its own output/.
#
# siesta.bands.WFSX/siesta.XV/calc.fdf/structure.fdf/calc.out are a REAL,
# spin-polarized SIESTA band-structure calculation of a CrS monolayer
# (fetched via stb-fetch from the twodmatpedia OPTIMADE database, id
# 2dm-2617 -- the same structure already used in examples 3.7/3.9/3.11).
# SystemLabel is "siesta", but there is no siesta.fdf anywhere -- the real
# input file is calc.fdf, and this exact mismatch is what caught two real
# bugs while building this walkthrough (see the README's own theory
# section, and the two dedicated sections below).

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
echo " What stb-wfdensity computes, and why it needs a real .fdf"
echo "=================================================================="
cat <<'EOF'
SIESTA is an LCAO code: every state is a linear combination of localized
atomic orbitals, psi_n(k, r) = sum_mu c_mu,n(k) * phi_mu(r). stb-fatbands
uses this to compute how much WEIGHT each orbital contributes; stb
-wfdensity instead evaluates psi's actual VALUE at every point of a real
-space grid, then takes |psi(r)|^2 -- a real, plottable 3D density for
ONE chosen (k, band) state, written as a Gaussian .cube file.

This needs the orbitals' real numerical radial shape -- NOT available
from a bare .XV or .HSX. Verified live on this example's own structure:

  siesta.XV -> 4 orbitals total (1 per atom), cutoff radius = -1 (no shape)
  calc.fdf  -> 64 real orbitals, cutoff radius = 3.61 Ang

Only a .fdf, read together with its .ion/.ion.xml basis files, carries
the real radial functions. This is exactly why --geometry-file (or
--label's auto-detected <label>.fdf) is required, never optional.

Also worth knowing before interpreting any single-state density: |psi_n
(k,r)|^2 for ONE band at ONE k-point is NOT, by itself, a physical
ground-state observable -- the real charge density sums this over EVERY
occupied band and k-point (that's exactly what stb-density's own .RHO
-based total density already is). This tool is a diagnostic for one state
at a time, not a substitute for that.

Every run prints a numbered [0]...[6] report. --save-report persists it;
--save-gnuplot writes the slice/profile .dat + a real .gplot script (off
by default -- this tool used to write the .dat unconditionally but never
actually wrote a .gplot at all, despite its own docstring calling this
"gnuplot output"); --view shows a matplotlib preview (this tool
previously had none at all). The .cube file itself is always written --
it's this tool's primary output.
EOF
pause

echo "=================================================================="
echo " output/basic-run/  --  --k-index/--band, the .cube file, normalization"
echo "=================================================================="
cat <<'EOF'
A plain run: band 1 (a deep, core-like state -- eigenvalue ~-79 eV) at
k-index 0. Watch [3] WAVEFUNCTION DENSITY report the normalization check
(integral |psi|^2 dV, should be ~1) -- a real, built-in quality control,
not decoration:
EOF
mkdir -p "$OUT/basic-run"
cp siesta.bands.WFSX calc.fdf structure.fdf Cr.ion Cr.ion.xml S.ion S.ion.xml "$OUT/basic-run/"
echo
echo "\$ stb-wfdensity --wfsx siesta.bands.WFSX --geometry-file calc.fdf --k-index 0 --band 1 --no-intro"
(cd "$OUT/basic-run" && stb-wfdensity --wfsx siesta.bands.WFSX --geometry-file calc.fdf \
    --k-index 0 --band 1 --no-intro > console.log 2>&1)
awk '/\[1\] INPUT DATA/{flag=1} /\[2\]/{flag=0} flag' "$OUT/basic-run/console.log"
awk '/\[2\] STATE SELECTION/{flag=1} /\[3\]/{flag=0} flag' "$OUT/basic-run/console.log"
awk '/\[3\] WAVEFUNCTION DENSITY/{flag=1} /\[4\]/{flag=0} flag' "$OUT/basic-run/console.log"
pause

echo "=================================================================="
echo " A real bug this exact fixture caught: --label + --geometry-file"
echo "=================================================================="
cat <<'EOF'
SystemLabel is "siesta" here, but the real input file is calc.fdf -- there
is NO siesta.fdf anywhere (this folder actually has THREE .fdf files:
calc.fdf, structure.fdf, and kpath_bs.fdf, a k-path-only fragment that
doesn't even parse as a geometry -- none named after the label). Watch
--label alone fail cleanly, then succeed once --geometry-file is added --
--label still auto-detects the .WFSX either way, only the geometry source
changes. This combination used to be rejected outright as "cannot be
combined" -- fixed, since load_parent() already preferred an explicit
--geometry-file over <label>.fdf on its own:
EOF
mkdir -p "$OUT/label-plus-geometry"
cp siesta.bands.WFSX siesta.XV calc.fdf structure.fdf kpath_bs.fdf Cr.ion Cr.ion.xml S.ion S.ion.xml "$OUT/label-plus-geometry/"
echo
echo "\$ stb-wfdensity --label siesta --k-index 0 --band 1 --no-intro   # no --geometry-file: fails"
(cd "$OUT/label-plus-geometry" && stb-wfdensity --label siesta --k-index 0 --band 1 \
    --no-intro > console_fail.log 2>&1) || true
tail -2 "$OUT/label-plus-geometry/console_fail.log"
echo
echo "\$ stb-wfdensity --label siesta --geometry-file calc.fdf --k-index 0 --band 1 --no-intro"
(cd "$OUT/label-plus-geometry" && stb-wfdensity --label siesta --geometry-file calc.fdf \
    --k-index 0 --band 1 --no-intro > console_ok.log 2>&1)
grep -E "WFSX file|Geometry source" "$OUT/label-plus-geometry/console_ok.log"
pause

echo "=================================================================="
echo " output/vbm-cbm/  --  searching the whole k-mesh, Fermi from --fermi-file"
echo "=================================================================="
cat <<'EOF'
--band vbm/cbm searches the ENTIRE k-mesh in the .WFSX for the global
extremum -- the same search stb-fatbands/stb-bands use. It needs a Fermi
energy: --fermi (explicit) > --bands-file > --fermi-file > an
auto-detected .out in the current directory, in that priority order.
Here we use --fermi-file calc.out directly (a real SIESTA log, NOT
assumed to be named after the label -- the same decoupling fixed for the
same reason as the geometry file above). This is a genuinely spin
-polarized calculation (nspin=2) -- watch [2] STATE SELECTION report which
spin channel the extremum was actually found in:
EOF
mkdir -p "$OUT/vbm-cbm"
cp siesta.bands.WFSX calc.fdf structure.fdf calc.out Cr.ion Cr.ion.xml S.ion S.ion.xml "$OUT/vbm-cbm/"
echo
echo "\$ stb-wfdensity --label siesta --geometry-file calc.fdf --band vbm --fermi-file calc.out --spacing 0.2 --no-intro"
(cd "$OUT/vbm-cbm" && stb-wfdensity --label siesta --geometry-file calc.fdf --band vbm \
    --fermi-file calc.out --spacing 0.2 --no-intro > console_vbm.log 2>&1)
awk '/\[2\] STATE SELECTION/{flag=1} /\[3\]/{flag=0} flag' "$OUT/vbm-cbm/console_vbm.log"
echo
echo "\$ stb-wfdensity --label siesta --geometry-file calc.fdf --band cbm --fermi-file calc.out --spacing 0.2 --no-intro"
(cd "$OUT/vbm-cbm" && stb-wfdensity --label siesta --geometry-file calc.fdf --band cbm \
    --fermi-file calc.out --spacing 0.2 --no-intro > console_cbm.log 2>&1)
awk '/\[2\] STATE SELECTION/{flag=1} /\[3\]/{flag=0} flag' "$OUT/vbm-cbm/console_cbm.log"
pause

echo "=================================================================="
echo " output/slice-position/  --  a real bug: where does the slice cut?"
echo "=================================================================="
cat <<'EOF'
This structure's 4 atoms sit at fractional z = 0, 0, 0.066, 0.934 --
clustered near the cell boundary, with the real vacuum in the MIDDLE of
the cell. --mode slice output used to always cut at the geometric center
of the cell along --axis, with no way to change it -- for THIS structure,
that landed squarely in the empty vacuum, showing essentially nothing.
Fixed: the default is now the plane where the planar-averaged |psi|^2 is
largest. Watch the default land at index 0 (Z=0.00, right at the atoms),
versus --pos at the old geometric-center position (deep in the vacuum):
EOF
mkdir -p "$OUT/slice-position"
cp siesta.bands.WFSX calc.fdf structure.fdf calc.out Cr.ion Cr.ion.xml S.ion S.ion.xml "$OUT/slice-position/"
echo
echo "\$ stb-wfdensity --label siesta --geometry-file calc.fdf --band vbm --fermi-file calc.out \\"
echo "      --spacing 0.2 --save-gnuplot --no-intro   # default: auto-detected |psi|^2 peak"
(cd "$OUT/slice-position" && stb-wfdensity --label siesta --geometry-file calc.fdf --band vbm \
    --fermi-file calc.out --spacing 0.2 --save-gnuplot --no-intro > console_default.log 2>&1)
grep -E "Save-gnuplot mode|Mapping plane" "$OUT/slice-position/console_default.log"
echo
echo "\$ stb-wfdensity --label siesta --geometry-file calc.fdf --band vbm --fermi-file calc.out \\"
echo "      --spacing 0.2 --pos 11.68 --save-gnuplot -o center --no-intro   # old center-of-cell behavior"
mkdir -p "$OUT/slice-position/center"
(cd "$OUT/slice-position" && stb-wfdensity --label siesta --geometry-file calc.fdf --band vbm \
    --fermi-file calc.out --spacing 0.2 --pos 11.68 --save-gnuplot -o center --no-intro \
    > console_center.log 2>&1)
grep -E "Save-gnuplot mode|Mapping plane" "$OUT/slice-position/console_center.log"
echo
echo "Peak |psi|^2 value at each slice (from the saved .dat files, column 4):"
echo -n "  default (auto-peak) max: "
awk '!/^#/{print $4}' "$OUT/slice-position/wfdensity_k139_b22_slice.dat" | sort -g | tail -1
echo -n "  --pos 11.68 (old center) max: "
awk '!/^#/{print $4}' "$OUT/slice-position/center/wfdensity_k139_b22_slice.dat" | sort -g | tail -1
pause

echo "=================================================================="
echo " output/grid-resolution/  --  --spacing and the normalization warning"
echo "=================================================================="
cat <<'EOF'
The SAME VBM state, two grid spacings -- watch the normalization check
(integral |psi|^2 dV, should be ~1) fail at the coarser default and pass
once tightened. Not a bug: a coarser grid genuinely clips more of the
real orbitals' own several-Angstrom tails:
EOF
mkdir -p "$OUT/grid-resolution"
cp siesta.bands.WFSX calc.fdf structure.fdf calc.out Cr.ion Cr.ion.xml S.ion S.ion.xml "$OUT/grid-resolution/"
echo
echo "\$ stb-wfdensity --label siesta --geometry-file calc.fdf --band vbm --fermi-file calc.out --spacing 0.3 --no-intro"
(cd "$OUT/grid-resolution" && stb-wfdensity --label siesta --geometry-file calc.fdf --band vbm \
    --fermi-file calc.out --spacing 0.3 --no-intro > console_0.3.log 2>&1)
grep -A1 "Normalization" "$OUT/grid-resolution/console_0.3.log"
echo
echo "\$ stb-wfdensity --label siesta --geometry-file calc.fdf --band vbm --fermi-file calc.out --spacing 0.2 --no-intro"
(cd "$OUT/grid-resolution" && stb-wfdensity --label siesta --geometry-file calc.fdf --band vbm \
    --fermi-file calc.out --spacing 0.2 --no-intro > console_0.2.log 2>&1)
grep -A1 "Normalization" "$OUT/grid-resolution/console_0.2.log"
pause

echo "=================================================================="
echo " output/full-report/  --  --save-report / --save-gnuplot (off by default)"
echo "=================================================================="
cat <<'EOF'
Every run always prints the numbered [0]...[6] report to the console.
--save-report additionally persists it to stb_wfdensity_report.txt, and
--save-gnuplot writes the slice .dat + a real .gplot script -- both off
by default, so a plain run only ever writes the .cube file + references.bib:
EOF
mkdir -p "$OUT/full-report"
cp siesta.bands.WFSX calc.fdf structure.fdf calc.out Cr.ion Cr.ion.xml S.ion S.ion.xml "$OUT/full-report/"
echo
echo "\$ stb-wfdensity --label siesta --geometry-file calc.fdf --band vbm --fermi-file calc.out \\"
echo "      --spacing 0.2 --no-intro   # default: no report, no slice data/gnuplot"
(cd "$OUT/full-report" && stb-wfdensity --label siesta --geometry-file calc.fdf --band vbm \
    --fermi-file calc.out --spacing 0.2 --no-intro > console_default.log 2>&1)
ls "$OUT/full-report/" | grep -v "^siesta\|^calc\|^structure\|console"
echo "(only the .cube file + references.bib -- no slice .dat/.gplot, no stb_wfdensity_report.txt)"
echo
echo "\$ stb-wfdensity --label siesta --geometry-file calc.fdf --band vbm --fermi-file calc.out \\"
echo "      --spacing 0.2 --save-report --save-gnuplot --no-intro"
(cd "$OUT/full-report" && stb-wfdensity --label siesta --geometry-file calc.fdf --band vbm \
    --fermi-file calc.out --spacing 0.2 --save-report --save-gnuplot --no-intro \
    > console_saved.log 2>&1)
echo "Report sections written to stb_wfdensity_report.txt:"
grep -E "^\[[0-9]+\]" "$OUT/full-report/stb_wfdensity_report.txt"
echo
echo "references.bib -- SIESTA:"
grep "^@" "$OUT/full-report/references.bib"
if command -v gnuplot > /dev/null; then
    (cd "$OUT/full-report" && gnuplot wfdensity_k139_b22_slice.gplot)
    echo "(rendered wfdensity_k139_b22_slice.pdf with the real, installed gnuplot)"
fi
pause

echo "=================================================================="
echo " Two ways to run it"
echo "=================================================================="
cat <<'EOF'
A -- direct CLI:
  stb-wfdensity --label siesta --geometry-file calc.fdf --band vbm --fermi-file calc.out

B -- interactive stb-suite menu:
  $ stb-suite
  Select an option (0-6, or a tool code like 4.1.2): 3.12

Both paths call the exact same underlying tool -- proven directly below.
The menu asks for the label and the .fdf path SEPARATELY (a real fix --
it used to only ask for one label and assume <label>.fdf named
everything), and offers a small Fermi-source submenu for VBM/CBM.
EOF
TMP="$(mktemp -d)"
cp siesta.bands.WFSX calc.fdf structure.fdf calc.out Cr.ion Cr.ion.xml S.ion S.ion.xml "$TMP/"
echo
echo "\$ printf '3.12\\nsiesta\\ncalc.fdf\\n1\\n0\\n1\\n0.3\\n\\n\\n\\n\\nn\\nn\\nn\\n\\n0\\n' | stb-suite"
(cd "$TMP" && printf '3.12\nsiesta\ncalc.fdf\n1\n0\n1\n0.3\n\n\n\n\nn\nn\nn\n\n0\n' | timeout 60 stb-suite > session.log 2>&1) || true
CLI_LINE=$(grep "Eigenvalue" "$OUT/basic-run/console.log")
MENU_LINE=$(grep "Eigenvalue" "$TMP/session.log")
if [ "$CLI_LINE" = "$MENU_LINE" ]; then
    echo "Confirmed: identical eigenvalue from the CLI and the interactive menu."
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
  basic-run/            label-plus-geometry/   vbm-cbm/
  slice-position/       grid-resolution/       full-report/

Each has references.bib and a .cube file; the --save-gnuplot runs
additionally have the slice .dat/.gplot, and full-report/ has
stb_wfdensity_report.txt (only from its --save-report run).

Recap of what this walkthrough covered:
  - LCAO wavefunctions evaluated at real-space points, |psi|^2, and why
    only a real .fdf (+ .ion/.ion.xml) carries the radial shapes needed --
    verified live (4 placeholder orbitals from .XV vs. 64 real ones from
    calc.fdf)
  - the normalization check (integral |psi|^2 dV ~ 1) as a genuine grid
    -resolution diagnostic, not decoration
  - |psi_n(k,r)|^2 for one state is NOT a physical ground-state observable
    by itself (stb-density's own .RHO-based total is)
  - a real bug: --label + --geometry-file used to be rejected outright;
    fixed, with this exact fixture (3 .fdf files, none named after the
    label) as the case that caught it
  - --band vbm/cbm's Fermi-energy priority order, ending in an
    auto-detected .out log decoupled from the label
  - a real bug: the slice used to always cut at the geometric center of
    the cell, landing in empty vacuum for this structure's off-center
    atoms; fixed to auto-detect the |psi|^2 peak instead, with --pos for
    manual override
  - --spacing vs. the normalization warning, demonstrated flipping live
  - the numbered [0]...[6] report, --save-report, --save-gnuplot (now
    with a REAL .gplot script, previously missing entirely), references.bib
  - CLI and the interactive stb-suite menu building the same command

As a next step, try on your own with a real SIESTA calculation:
  stb-wfdensity --label my_calc --band vbm --fermi-file my_calc.out --view
  stb-wfdensity --label my_calc --k-point 0.0 0.0 0.0 --band 5 --spacing 0.2
EOF
