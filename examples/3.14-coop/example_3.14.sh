#!/bin/bash
# Guided example: stb-coop (code 3.14 in the stb-suite menu)
#
# Not an automated test (see test/3-analysis/14-coop/test.sh for that) -- a
# commented walk-through: it runs real commands, one group at a time, into
# its own output/<case>/ folder, and shows you the piece of output that
# proves what just happened. Pauses between sections so you can read
# before moving on. Safe to re-run any time -- it always starts by wiping
# its own output/.
#
# Sn3O4.HSX/Sn3O4.selected.WFSX/Sn3O4.DM are a REAL, non-polarized
# GGA-PBE(+DFT-D3) SIESTA calculation of the tin oxide Sn3O4 (14 atoms,
# 182 orbitals, DZP basis), sampled at 8 explicit k-points (a
# 2x2x2-equivalent full-BZ mesh via WriteWaveFunctions T + %block
# WaveFuncKPoints -- NOT a band-path WFSX like stb-fatbands needs).
# Atom 0 (Sn) has a real Sn-O bond to atom 6 (2.077 Ang) and sits 3.800
# Ang from atom 1 (another Sn, NOT bonded) -- these two pairs are used
# throughout this walkthrough as a real bonded-vs-non-bonded contrast.

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
echo " What COOP/COHP measure, and why the energy window matters"
echo "=================================================================="
cat <<'EOF'
For a crystal orbital psi_n(k) = sum_mu c_mu,n(k) phi_mu, the overlap
population between orbital mu (atom I) and nu (atom J) is
c_mu,n(k) * c_nu,n(k) * S_mu,nu(k) -- positive for an in-phase (bonding)
combination, negative for out-of-phase (antibonding). Summed over every
orbital pair, band, and k-point (Gaussian-broadened around each
eigenvalue), this gives a continuous COOP(E) curve. COHP replaces the
overlap matrix S with the Hamiltonian matrix H in the same formula.

This suite's own sign convention: COOP/COHP > 0 = bonding.

Needs a full-BZ .WFSX (WriteWaveFunctions T + %block WaveFuncKPoints,
NOT a band-path one like stb-fatbands uses) and a real .HSX -- unlike
stb-wfdensity/stb-sts, there is NO approximate fallback: COOP needs the
real overlap matrix, COHP needs the real Hamiltonian.

Every run prints a numbered [0]...[7] report. --save-report persists it;
--save-gnuplot writes coop.dat/cohp.dat + a real .gplot script (off by
default -- this tool used to write the .dat unconditionally but never a
.gplot at all); --view shows a matplotlib preview (replaces the old
--no-plot, which used to be on by default).
EOF
pause


echo "=================================================================="
echo " output/basic-run/  --  a misleading window (includes unoccupied states)"
echo "=================================================================="
mkdir -p "$OUT/basic-run"
cp Sn3O4.HSX Sn3O4.selected.WFSX "$OUT/basic-run/"
cat <<'EOF'
Atom 0 (Sn) has a REAL bond to atom 6 (O, 2.077 Ang) and is NOT bonded to
atom 1 (Sn, 3.800 Ang). Watch [3] COOP CURVE with a window that goes well
ABOVE the real Fermi energy (-3.204 eV) into unoccupied states -- the
result looks backwards, and that's the point of the next section:
EOF
echo "\$ stb-coop --label Sn3O4 --quantity coop --pair 0 6 --pair 0 1 --erange -10 5 --npoints 25 --sigma 300 --no-intro"
(cd "$OUT/basic-run" && stb-coop --label Sn3O4 --quantity coop --pair 0 6 --pair 0 1 \
    --erange -10 5 --npoints 25 --sigma 300 --no-intro | tee console.log | \
    sed -n '/\[3\] COOP CURVE/,/\[4\] BOND ORDER/p')
pause


echo "=================================================================="
echo " output/occupied-window/  --  restricting to the real occupied states"
echo "=================================================================="
cat <<'EOF'
The SAME two pairs, now integrated only over the occupied states
(up to the real Fermi energy, -3.204 eV) -- the physically meaningful
window. Watch the real bond (0-6) now correctly read BONDING, and the
non-bond (0-1) collapse to essentially zero:
EOF
mkdir -p "$OUT/occupied-window"
cp Sn3O4.HSX Sn3O4.selected.WFSX "$OUT/occupied-window/"
echo "\$ stb-coop --label Sn3O4 --quantity coop --pair 0 6 --pair 0 1 --erange -25 -3.204 --npoints 40 --sigma 300 --no-intro"
(cd "$OUT/occupied-window" && stb-coop --label Sn3O4 --quantity coop --pair 0 6 --pair 0 1 \
    --erange -25 -3.204 --npoints 40 --sigma 300 --no-intro | \
    sed -n '/\[3\] COOP CURVE/,/\[4\] BOND ORDER/p')
pause


echo "=================================================================="
echo " A real, serious bug this fixture caught: --bond-order"
echo "=================================================================="
cat <<'EOF'
--bond-order is a cheap, energy-INTEGRATED cross-check on the curve
above: the Mulliken bond order, computed from a REAL density matrix
(--dm-file, or an auto-detected <label>.DM). This tool used to call
H.bond_order() directly on the Hamiltonian instead -- silently using
Hamiltonian matrix elements (eV-scale) as if they were a density-matrix
population. Verified live before the fix: the on-site term came out as
-92.38, a real Sn-O bond as -38.32, and the NON-bonded Sn-Sn pair as
-29.04 -- nearly the SAME magnitude regardless of real bonding, the
smoking gun that this measured nothing physical at all.

Fixed by reading Sn3O4.DM and splicing in the Hamiltonian's own overlap
column (SIESTA's .DM format doesn't store overlap at all). Watch the
fixed numbers below, AND the automatic cross-check against the
occupied-window COOP curve from the previous section:
EOF
mkdir -p "$OUT/bond-order"
cp Sn3O4.HSX Sn3O4.selected.WFSX Sn3O4.DM "$OUT/bond-order/"
echo "\$ stb-coop --label Sn3O4 --quantity coop --pair 0 6 --pair 0 1 --erange -25 -3.204 --npoints 40 --sigma 300 --bond-order --no-intro"
(cd "$OUT/bond-order" && stb-coop --label Sn3O4 --quantity coop --pair 0 6 --pair 0 1 \
    --erange -25 -3.204 --npoints 40 --sigma 300 --bond-order --no-intro | \
    sed -n '/\[4\] BOND ORDER/,/\[5\] OUTPUT/p')
pause


echo "=================================================================="
echo " output/cohp/  --  --quantity cohp with --pair-species"
echo "=================================================================="
cat <<'EOF'
--pair-species aggregates ALL atom-index pairs between two species into
one combined curve -- here every Sn-O interaction in the cell at once,
using the Hamiltonian instead of the overlap matrix:
EOF
mkdir -p "$OUT/cohp"
cp Sn3O4.HSX Sn3O4.selected.WFSX "$OUT/cohp/"
echo "\$ stb-coop --label Sn3O4 --quantity cohp --pair-species Sn O --erange -10 5 --npoints 25 --sigma 300 --no-intro"
(cd "$OUT/cohp" && stb-coop --label Sn3O4 --quantity cohp --pair-species Sn O \
    --erange -10 5 --npoints 25 --sigma 300 --no-intro | \
    sed -n '/\[2\] PAIR/,/\[4\] BOND ORDER/p')
pause


echo "=================================================================="
echo " A real bug this tool inherited from stb-wfdensity/stb-sts: --label + --hsx-file"
echo "=================================================================="
cat <<'EOF'
SystemLabel is "Sn3O4" here (matching this fixture's own filenames), but
real calculations very often have a mismatched name (e.g. calc.HSX with
SystemLabel "siesta") -- this combination used to be rejected outright.
Simulating that exact scenario with a renamed copy:
EOF
mkdir -p "$OUT/label-plus-hsx"
cp Sn3O4.selected.WFSX "$OUT/label-plus-hsx/"
cp Sn3O4.HSX "$OUT/label-plus-hsx/calc.HSX"
echo "\$ stb-coop --label Sn3O4 --quantity coop --pair 0 6 --erange -10 5 --npoints 15 --sigma 300 --no-intro   # no --hsx-file: fails (no Sn3O4.HSX here)"
(cd "$OUT/label-plus-hsx" && stb-coop --label Sn3O4 --quantity coop --pair 0 6 \
    --erange -10 5 --npoints 15 --sigma 300 --no-intro > console_fail.log 2>&1) || true
grep -A1 "Resolving Hamiltonian" "$OUT/label-plus-hsx/console_fail.log"
echo
echo "\$ stb-coop --label Sn3O4 --hsx-file calc.HSX --quantity coop --pair 0 6 --erange -10 5 --npoints 15 --sigma 300 --no-intro"
(cd "$OUT/label-plus-hsx" && stb-coop --label Sn3O4 --hsx-file calc.HSX --quantity coop --pair 0 6 \
    --erange -10 5 --npoints 15 --sigma 300 --no-intro > console_ok.log 2>&1)
grep "WFSX file\|Hamiltonian src" "$OUT/label-plus-hsx/console_ok.log"
pause


echo "=================================================================="
echo " output/fermi-shift/  --  --shift fermi's decoupled-from-label hierarchy"
echo "=================================================================="
cat <<'EOF'
--shift fermi now has the same priority-ordered hierarchy stb-wfdensity/
stb-sts have: --fermi (explicit) > --bands-file > --fermi-file > an
auto-detected .out log in the current directory. (No real .out was saved
alongside this fixture's .WFSX/.HSX, so calc.out below is a
representative demonstration value, not one re-derived from this exact
run.)
EOF
mkdir -p "$OUT/fermi-shift"
cp Sn3O4.HSX Sn3O4.selected.WFSX "$OUT/fermi-shift/"
cat > "$OUT/fermi-shift/calc.out" << 'EOF'
Some SIESTA log noise before the summary
siesta:         Fermi = -3.203856
More noise after
EOF
echo "\$ stb-coop --label Sn3O4 --quantity coop --pair 0 6 --erange -3 3 --npoints 15 --sigma 300 \\"
echo "      --shift fermi --fermi -3.203856 --no-intro   # explicit value"
(cd "$OUT/fermi-shift" && stb-coop --label Sn3O4 --quantity coop --pair 0 6 \
    --erange -3 3 --npoints 15 --sigma 300 --shift fermi --fermi -3.203856 \
    --no-intro > console_explicit.log 2>&1)
grep "Energy shift" "$OUT/fermi-shift/console_explicit.log"
echo
echo "\$ stb-coop --label Sn3O4 --quantity coop --pair 0 6 --erange -3 3 --npoints 15 --sigma 300 \\"
echo "      --shift fermi --no-intro   # no --fermi/--fermi-file: auto-detects calc.out"
(cd "$OUT/fermi-shift" && stb-coop --label Sn3O4 --quantity coop --pair 0 6 \
    --erange -3 3 --npoints 15 --sigma 300 --shift fermi \
    --no-intro > console_auto.log 2>&1)
grep "Energy shift" "$OUT/fermi-shift/console_auto.log"
pause


echo "=================================================================="
echo " output/full-report/  --  --save-report / --save-gnuplot (off by default)"
echo "=================================================================="
cat <<'EOF'
Every run always prints the numbered [0]...[7] report to the console.
--save-report additionally persists it to stb_coop_report.txt, and
--save-gnuplot writes coop.dat + a real .gplot script -- both off by
default, so a plain run only ever writes coop.dat + references.bib:
EOF
mkdir -p "$OUT/full-report"
cp Sn3O4.HSX Sn3O4.selected.WFSX "$OUT/full-report/"
echo "\$ stb-coop --label Sn3O4 --quantity coop --pair 0 6 --pair 0 1 --erange -10 5 --npoints 25 --sigma 300 --no-intro   # default"
(cd "$OUT/full-report" && stb-coop --label Sn3O4 --quantity coop --pair 0 6 --pair 0 1 \
    --erange -10 5 --npoints 25 --sigma 300 --no-intro > console_default.log 2>&1)
(cd "$OUT/full-report" && ls)
echo "(only coop.dat + references.bib -- no .gplot, no stb_coop_report.txt)"
echo
echo "\$ stb-coop --label Sn3O4 --quantity coop --pair 0 6 --pair 0 1 --erange -10 5 --npoints 25 --sigma 300 \\"
echo "      --save-report --save-gnuplot --no-intro"
(cd "$OUT/full-report" && stb-coop --label Sn3O4 --quantity coop --pair 0 6 --pair 0 1 \
    --erange -10 5 --npoints 25 --sigma 300 --save-report --save-gnuplot \
    --no-intro > console_saved.log 2>&1)
echo "Report sections written to stb_coop_report.txt:"
grep -o '\[[0-9]\] [A-Z& ]*' "$OUT/full-report/stb_coop_report.txt" | sort -u
echo
if command -v gnuplot >/dev/null 2>&1; then
    (cd "$OUT/full-report" && gnuplot coop.gplot)
    echo "(rendered coop.pdf with the real, installed gnuplot)"
fi
pause


echo "=================================================================="
echo " Two ways to run it"
echo "=================================================================="
cat <<EOF
A -- direct CLI:
  stb-coop --label Sn3O4 --quantity coop --pair 0 6 --erange -10 5 --sigma 300

B -- interactive stb-suite menu:
  \$ stb-suite
  Select an option (0-6, or a tool code like 4.1.2): 3.14

Both paths call the exact same underlying tool -- proven directly below.
The menu asks for the label and the .HSX path SEPARATELY, offers an
energy-shift submenu (with the same Fermi-source options shown above),
and a bond-order/.DM prompt.
EOF
TMP="$(mktemp -d)"
cp Sn3O4.HSX Sn3O4.selected.WFSX "$TMP/"
echo
echo "\$ stb-coop --label Sn3O4 --quantity coop --pair 0 6 --erange -10 5 --npoints 25 --sigma 300 --no-intro"
(cd "$TMP" && stb-coop --label Sn3O4 --quantity coop --pair 0 6 \
    --erange -10 5 --npoints 25 --sigma 300 --no-intro > cli.log 2>&1)
echo
echo "\$ printf '3.14\\nSn3O4\\n\\n1\\n1\\n0\\n6\\n-10\\n5\\n300\\n\\nn\\n\\nn\\nn\\nn\\n\\n0\\n' | stb-suite"
(cd "$TMP" && printf '3.14\nSn3O4\n\n1\n1\n0\n6\n-10\n5\n300\n\nn\n\nn\nn\nn\n\n0\n' | timeout 120 stb-suite > session.log 2>&1) || true
CLI_LINE=$(grep "Selected pairs" "$TMP/cli.log")
MENU_LINE=$(grep "Selected pairs" "$TMP/session.log")
if [ "$CLI_LINE" = "$MENU_LINE" ]; then
    echo "Confirmed: identical pair selection from the CLI and the interactive menu."
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
Seven self-contained folders were generated under output/:
  basic-run/          occupied-window/   bond-order/
  cohp/                label-plus-hsx/    fermi-shift/
  full-report/

Each has references.bib and coop.dat/cohp.dat; full-report/ additionally
has coop.gplot/coop.pdf and stb_coop_report.txt (only from its
--save-report/--save-gnuplot run).

Recap of what this walkthrough covered:
  - COOP/COHP theory: crystal-orbital overlap/Hamilton populations,
    positive = bonding by this suite's (now cross-checked) sign
    convention
  - why the chosen --erange matters for interpretation: a window that
    includes unoccupied states can make a real bond look antibonding
  - a real, serious bug: --bond-order used to call H.bond_order() on the
    Hamiltonian instead of a real density matrix, giving numbers off by
    orders of magnitude (verified: -92/-38/-29 for on-site/bonded/
    non-bonded, all roughly the same size); fixed to read a real .DM and
    splice in the Hamiltonian's own overlap column
  - the first independent cross-check of this suite's own COOP sign
    convention: integrated COOP and the fixed bond order agree in sign
    for both a bonded and a non-bonded pair
  - --quantity cohp with --pair-species (aggregate curves)
  - a real bug: --label + --hsx-file used to be rejected outright; fixed
    the same way as stb-wfdensity/stb-sts
  - --shift fermi's new priority-ordered Fermi-source hierarchy, shared
    with stb-wfdensity/stb-sts
  - the numbered [0]...[7] report, --save-report, --save-gnuplot (now
    with a REAL .gplot script, previously missing entirely), --view
    (replacing the old, on-by-default --no-plot)
  - CLI and the interactive stb-suite menu building the same command

As a next step, try on your own with a real SIESTA calculation:
  stb-coop --label my_calc --quantity coop --pair 0 1 --erange -10 5 --sigma 200 \\
      --bond-order --save-report --save-gnuplot --view
EOF
