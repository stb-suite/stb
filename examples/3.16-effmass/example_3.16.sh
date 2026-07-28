#!/bin/bash
# Guided example: stb-effmass (code 3.16 in the stb-suite menu)
#
# Not an automated test (see test/3-analysis/16-effmass/test.sh for that) --
# a commented walk-through: it runs real commands, one group at a time,
# into its own output/<case>/ folder, and shows you the piece of output
# that proves what just happened. Pauses between sections so you can read
# before moving on. Safe to re-run any time -- it always starts by wiping
# its own output/.
#
# Sn3O4.selected.WFSX/Sn3O4.HSX are the exact same real SIESTA fixture
# 3.14-coop/ uses (a full-BZ 2x2x2-equivalent k-mesh calculation).
# spin/Ospin.* is a real spin-polarized isolated-oxygen-atom calculation
# (2 Bohr-magneton triplet ground state) -- the exact fixture that
# reproduces the real crash a user hit on their own calculation.
# nc/IsolatedO.* is a real non-collinear/spin-orbit calculation.

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
echo " Welcome: band velocity and effective mass, from first principles"
echo "=================================================================="
cat <<'EOF'
Near a band extremum, E_n(k) ~ E_n(k0) + hbar*v.(k-k0) + (hbar^2/2)
(k-k0)^T [M^-1] (k-k0) + ... -- v is the band VELOCITY (1st derivative),
M is the EFFECTIVE MASS TENSOR (from the 2nd derivative). Both come from
sisl's ANALYTIC k-derivatives of the Hamiltonian -- no finite-difference
step size to tune.

IMPORTANT (verified live): the 2nd-order curvature correction needs
energy differences between ALL bands at that k-point, so it must always
be computed on the FULL multi-band eigenstate, never after isolating one
band -- doing so silently gives values wrong by orders of magnitude.

Every run prints a numbered [0]...[8] report. --save-report persists it;
--save-gnuplot writes effmass.dat/velocity.dat + a real 2-panel bar-chart
.gplot script (off by default); --view shows the same as a matplotlib
preview (off by default).
EOF
pause


echo "=================================================================="
echo " output/basic-run/  --  a deep state at Gamma: per-axis vs. principal masses"
echo "=================================================================="
mkdir -p "$OUT/basic-run"
cp Sn3O4.selected.WFSX Sn3O4.HSX "$OUT/basic-run/"
cat <<'EOF'
sisl's own per-axis Voigt values (m*_xx/yy/zz/yz/xz/xy) are each an
INDEPENDENT element-wise reciprocal of the curvature tensor -- NOT a
proper matrix inversion. Watch [4] properly assemble the FULL tensor,
invert it, and diagonalize it into the true PRINCIPAL effective masses,
flagging a [WARNING] when the off-diagonal curvature is significant:
EOF
echo "\$ stb-effmass --wfsx Sn3O4.selected.WFSX --hsx-file Sn3O4.HSX --k-index 0 --band 1 --no-intro"
(cd "$OUT/basic-run" && stb-effmass --wfsx Sn3O4.selected.WFSX --hsx-file Sn3O4.HSX \
    --k-index 0 --band 1 --no-intro | tee console.log | sed -n '/\[3\] EFFECTIVE MASS/,/\[5\] BAND/p')
pause


echo "=================================================================="
echo " output/vbm/  --  the valence band maximum: negative (hole-like) mass"
echo "=================================================================="
cat <<'EOF'
At a true VBM, every properly-diagonalized principal effective mass
should be NEGATIVE (the band curves downward away from the maximum) --
a textbook-correct sign check the per-axis Voigt values alone don't
guarantee to show cleanly:
EOF
mkdir -p "$OUT/vbm"
cp Sn3O4.selected.WFSX Sn3O4.HSX "$OUT/vbm/"
echo "\$ stb-effmass --wfsx Sn3O4.selected.WFSX --hsx-file Sn3O4.HSX --band vbm --fermi -3.200055 --no-intro"
(cd "$OUT/vbm" && stb-effmass --wfsx Sn3O4.selected.WFSX --hsx-file Sn3O4.HSX \
    --band vbm --fermi -3.200055 --no-intro | sed -n '/\[4\] EFFECTIVE MASS (principal/,/\[5\] BAND/p')
pause


echo "=================================================================="
echo " output/cbm/  --  a caution about coarse k-meshes"
echo "=================================================================="
cat <<'EOF'
At the CBM, watch the principal masses come out MIXED-SIGN -- not the
all-positive signature a genuine local minimum should have. This
fixture's WFSX only samples a coarse 8-k-point mesh; the identified
"global CBM" may not be an exact true minimum in every direction of the
real, continuous band structure -- an honest limitation of ANY
mesh-identified extremum, not a bug:
EOF
mkdir -p "$OUT/cbm"
cp Sn3O4.selected.WFSX Sn3O4.HSX "$OUT/cbm/"
echo "\$ stb-effmass --wfsx Sn3O4.selected.WFSX --hsx-file Sn3O4.HSX --band cbm --fermi -3.200055 --no-intro"
(cd "$OUT/cbm" && stb-effmass --wfsx Sn3O4.selected.WFSX --hsx-file Sn3O4.HSX \
    --band cbm --fermi -3.200055 --no-intro | sed -n '/\[4\] EFFECTIVE MASS (principal/,/\[5\] BAND/p')
pause


echo "=================================================================="
echo " THE REAL BUG: a spin-polarized calculation used to crash"
echo "=================================================================="
cat <<'EOF'
A user ran this tool on their own real, spin-polarized SIESTA calculation
and hit: TypeError: SparseOrbitalBZ._ddPk() got an unexpected keyword
argument 'spin'. Verified directly against sisl's own source: ddPk()
(the 2nd-order derivative effective_mass() needs) has NO 'spin' parameter
at all in this sisl version, while dPk() (1st-order, used by velocity())
does -- this tool used to only guard non-collinear/SOC (nspin=4/8),
wrongly assuming plain spin-polarized (nspin=2) was safe. Reproducing the
EXACT crash on a real spin-polarized isolated-oxygen-atom fixture (2
Bohr-magneton triplet ground state) -- now fixed, no traceback:
EOF
mkdir -p "$OUT/nspin-fix"
cp spin/Ospin.bands.WFSX spin/Ospin.HSX "$OUT/nspin-fix/"
echo "\$ stb-effmass --wfsx Ospin.bands.WFSX --hsx-file Ospin.HSX --k-index 0 --band 3 --no-intro"
(cd "$OUT/nspin-fix" && stb-effmass --wfsx Ospin.bands.WFSX --hsx-file Ospin.HSX \
    --k-index 0 --band 3 --no-intro | sed -n '/\[1\] INPUT DATA/,/\[5\] BAND/p')
pause


echo "=================================================================="
echo " Same graceful degradation, non-collinear/SOC (pre-existing case)"
echo "=================================================================="
cat <<'EOF'
This exact treatment already existed for nspin=4/8 before the fix above --
the fix simply extends the SAME logic to nspin=2, which used to be
(incorrectly) assumed safe:
EOF
mkdir -p "$OUT/nspin4"
cp nc/IsolatedO.selected.WFSX nc/IsolatedO.HSX "$OUT/nspin4/"
echo "\$ stb-effmass --wfsx IsolatedO.selected.WFSX --hsx-file IsolatedO.HSX --k-index 0 --band 3 --no-intro"
(cd "$OUT/nspin4" && stb-effmass --wfsx IsolatedO.selected.WFSX --hsx-file IsolatedO.HSX \
    --k-index 0 --band 3 --no-intro | grep "Spin channels\|not supported by sisl")
pause


echo "=================================================================="
echo " A real bug this tool inherited from stb-wfdensity/stb-sts/stb-coop: --label + --hsx-file"
echo "=================================================================="
cat <<'EOF'
This combination used to be rejected outright. Simulating a renamed .HSX
(a common real-world mismatch, e.g. SystemLabel not matching the real
input filename):
EOF
mkdir -p "$OUT/label-plus-hsx"
cp Sn3O4.selected.WFSX "$OUT/label-plus-hsx/"
cp Sn3O4.HSX "$OUT/label-plus-hsx/calc.HSX"
echo "\$ stb-effmass --label Sn3O4 --k-index 0 --band 1 --no-intro   # no --hsx-file: fails (no Sn3O4.HSX here)"
(cd "$OUT/label-plus-hsx" && stb-effmass --label Sn3O4 --k-index 0 --band 1 \
    --no-intro > console_fail.log 2>&1) || true
grep -A1 "Resolving Hamiltonian" "$OUT/label-plus-hsx/console_fail.log"
echo
echo "\$ stb-effmass --label Sn3O4 --hsx-file calc.HSX --k-index 0 --band 1 --no-intro"
(cd "$OUT/label-plus-hsx" && stb-effmass --label Sn3O4 --hsx-file calc.HSX --k-index 0 --band 1 \
    --no-intro > console_ok.log 2>&1)
grep "WFSX file\|Hamiltonian src" "$OUT/label-plus-hsx/console_ok.log"
pause


echo "=================================================================="
echo " output/fermi-shift/  --  --band vbm/cbm's decoupled-from-label Fermi hierarchy"
echo "=================================================================="
cat <<'EOF'
--fermi (explicit) > --bands-file > --fermi-file > an auto-detected .out
log in the current directory -- the same hierarchy stb-wfdensity/stb-sts/
stb-coop already have. (No real .out was saved alongside this fixture's
.WFSX/.HSX, so calc.out below is a representative demonstration value.)
EOF
mkdir -p "$OUT/fermi-shift"
cp Sn3O4.selected.WFSX Sn3O4.HSX "$OUT/fermi-shift/"
cat > "$OUT/fermi-shift/calc.out" << 'EOF'
Some SIESTA log noise before the summary
siesta:         Fermi = -3.200055
More noise after
EOF
echo "\$ stb-effmass --wfsx Sn3O4.selected.WFSX --hsx-file Sn3O4.HSX --band vbm --no-intro   # auto-detects calc.out"
(cd "$OUT/fermi-shift" && stb-effmass --wfsx Sn3O4.selected.WFSX --hsx-file Sn3O4.HSX \
    --band vbm --no-intro > console.log 2>&1)
grep "Band selection" "$OUT/fermi-shift/console.log"
pause


echo "=================================================================="
echo " output/full-report/  --  --save-report / --save-gnuplot (off by default)"
echo "=================================================================="
cat <<'EOF'
Every run always prints the numbered [0]...[8] report to the console.
--save-report additionally persists it to stb_effmass_report.txt, and
--save-gnuplot writes effmass.dat/velocity.dat + a real 2-panel bar-chart
.gplot script -- both off by default:
EOF
mkdir -p "$OUT/full-report"
cp Sn3O4.selected.WFSX Sn3O4.HSX "$OUT/full-report/"
echo "\$ stb-effmass --wfsx Sn3O4.selected.WFSX --hsx-file Sn3O4.HSX --band vbm --fermi -3.200055 --no-intro   # default"
(cd "$OUT/full-report" && stb-effmass --wfsx Sn3O4.selected.WFSX --hsx-file Sn3O4.HSX \
    --band vbm --fermi -3.200055 --no-intro > console_default.log 2>&1)
(cd "$OUT/full-report" && ls)
echo "(only references.bib -- no effmass.dat/velocity.dat/.gplot, no stb_effmass_report.txt)"
echo
echo "\$ stb-effmass --wfsx Sn3O4.selected.WFSX --hsx-file Sn3O4.HSX --band vbm --fermi -3.200055 \\"
echo "      --save-report --save-gnuplot --no-intro"
(cd "$OUT/full-report" && stb-effmass --wfsx Sn3O4.selected.WFSX --hsx-file Sn3O4.HSX \
    --band vbm --fermi -3.200055 --save-report --save-gnuplot --no-intro > console_saved.log 2>&1)
echo "Report sections written to stb_effmass_report.txt:"
grep -o "^\[[0-9]\] [A-Za-z& (),.'-]*" "$OUT/full-report/stb_effmass_report.txt" | sort -u
echo
if command -v gnuplot >/dev/null 2>&1; then
    (cd "$OUT/full-report" && gnuplot effmass.gplot)
    echo "(rendered effmass.pdf with the real, installed gnuplot)"
fi
pause


echo "=================================================================="
echo " Two ways to run it"
echo "=================================================================="
cat <<EOF
A -- direct CLI:
  stb-effmass --wfsx Sn3O4.selected.WFSX --hsx-file Sn3O4.HSX --band vbm --fermi -3.200055

B -- interactive stb-suite menu:
  \$ stb-suite
  Select an option (0-6, or a tool code like 4.1.2): 3.16

Both paths call the exact same underlying tool -- proven directly below.
The menu asks for the label and the .HSX path SEPARATELY (default
suggested if <label>.HSX exists), gained a Fermi-source submenu for
VBM/CBM, and save-report/save-gnuplot/view prompts.
EOF
TMP="$(mktemp -d)"
cp Sn3O4.selected.WFSX Sn3O4.HSX "$TMP/"
echo
echo "\$ stb-effmass --wfsx Sn3O4.selected.WFSX --hsx-file Sn3O4.HSX --band vbm --fermi -3.200055 --no-intro"
(cd "$TMP" && stb-effmass --wfsx Sn3O4.selected.WFSX --hsx-file Sn3O4.HSX --band vbm \
    --fermi -3.200055 --no-intro > cli.log 2>&1)
echo
echo "\$ printf '3.16\\nSn3O4\\n\\n2\\n1\\n-3.200055\\n\\nn\\nn\\nn\\n\\n0\\n' | stb-suite"
(cd "$TMP" && printf '3.16\nSn3O4\n\n2\n1\n-3.200055\n\nn\nn\nn\n\n0\n' | timeout 60 stb-suite > session.log 2>&1) || true
CLI_LINE=$(grep "Eigenvalue" "$TMP/cli.log")
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
Eight self-contained folders were generated under output/:
  basic-run/       vbm/            cbm/             nspin-fix/
  nspin4/          label-plus-hsx/ fermi-shift/     full-report/

Each has references.bib; full-report/ additionally has effmass.dat/
velocity.dat/effmass.gplot/effmass.pdf and stb_effmass_report.txt.

Recap of what this walkthrough covered:
  - band velocity and effective mass as the 1st/2nd derivatives of E(k),
    computed analytically (no finite-difference step size)
  - THE REAL BUG: sisl's ddPk() has no 'spin' parameter at all, so
    effective_mass() used to crash for ANY spin-polarized calculation
    (not just non-collinear/SOC as previously assumed) -- reproduced and
    fixed live on a real spin-polarized isolated-oxygen-atom fixture
  - a real, deeper finding: sisl's own per-axis Voigt values are
    element-wise reciprocals, not a true tensor inversion -- properly
    diagonalizing the FULL curvature tensor gives genuinely different
    (and more physically meaningful) principal effective masses,
    verified with a negative-mass sign check at the VBM
  - a real, honest limitation: a coarse-mesh-identified CBM can give
    mixed-sign principal masses if it isn't an exact true minimum
  - a real bug: --label + --hsx-file used to be rejected outright; fixed
    the same way as stb-wfdensity/stb-sts/stb-coop
  - --band vbm/cbm's Fermi-source hierarchy, shared with those same tools
  - the numbered [0]...[8] report, --save-report, --save-gnuplot (a real
    2-panel bar-chart script, previously no plot output existed at all),
    --view

As a next step, try on your own with a real SIESTA calculation:
  stb-effmass --label my_calc --band vbm --fermi-file my_calc.out --save-gnuplot --view
EOF
