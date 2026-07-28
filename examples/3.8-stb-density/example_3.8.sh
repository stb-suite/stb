#!/bin/bash
# Guided example: stb-density (code 3.8 in the stb-suite menu)
#
# Not an automated test (see test/3-analysis/8-density/test.sh for that).
# Uses two real, finished SIESTA calculations:
#   Sn3O4.RHO           -- a non-spin-polarized 14-atom Sn3O4 structure
#                          (copied from test/3-analysis/8-density/)
#   o2.RHO/o2.XV/o2.fdf -- a spin-polarized isolated O2 molecule, a real
#                          textbook triplet ground state (copied from
#                          test/6-utils/3-cube/) -- this exact fixture is
#                          what verified the up/down-vs-total/spin fix
#                          documented in the README.

set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

OUT="$DIR/output"
rm -rf "$OUT"
mkdir -p "$OUT"

export MPLBACKEND=Agg

pause() {
    echo
    read -p "  [Press Enter to continue] " -r
    echo
}

echo "=================================================================="
echo " What stb-density computes, and why"
echo "=================================================================="
cat <<'EOF'
stb-density reads a SIESTA charge-density grid (.RHO) and exports it as a
2D slice, a 3D point cloud, or a planar-averaged 1D profile -- with a
matplotlib preview and/or gnuplot .dat+.gplot export. If the .RHO is
spin-polarized, the net spin (magnetization) density is now detected and
processed automatically alongside the total charge -- no flag needed.

See the README for the full theory: what the charge density itself shows
you, why the spin-polarized case needed a real, verified fix (see the
next section for a live demonstration), planar averaging, and the
signed-vs-non-negative colorbar convention.

Every run prints a numbered [0]...[5] report. --save-report persists it
to stb_density_report.txt; --save-gnuplot writes the .dat+.gplot pair
(off by default -- this tool used to write the gnuplot script
UNCONDITIONALLY; that's no longer the case); --view shows a matplotlib
preview (heatmap/line/3D scatter, matching the mode).
EOF
pause

echo "=================================================================="
echo " A real, verified physics bug: up/down channels vs. total/spin"
echo "=================================================================="
cat <<'EOF'
A 2-component (spin-polarized) .RHO stores the raw UP-spin and DOWN-spin
densities directly -- NOT "total" and "spin" as separately meaningful
quantities. o2.RHO is a real SIESTA calculation on an isolated O2
molecule (a textbook triplet ground state -- SIESTA's own log reports
|S| = 2.0). Watch the naive (WRONG) reading below vs. stb-density's
actual (correct) one:
EOF
python3 -c "
import sisl
sile = sisl.get_sile('o2.RHO')

# Naive (WRONG): raw components 0 and 1 directly.
naive_charge = sile.read_grid(index=0)
naive_spin = sile.read_grid(index=1)

# Correct: sisl's own combination convention (up+down, up-down).
real_total = sile.read_grid(index='total')
real_spin = sile.read_grid(index='z')

import numpy as np
def integrate(grid):
    cell_vol = abs(np.linalg.det(grid.lattice.cell))
    voxel = cell_vol / grid.grid.size
    return grid.grid.sum() * voxel

print(f'  Naive index=0 as \"charge\":  {integrate(naive_charge):.4f} e   (WRONG -- this is just the up channel)')
print(f'  Naive index=1 as \"spin\":    {integrate(naive_spin):.4f} e   (WRONG -- this is just the down channel)')
print(f'  Correct total (up+down):   {integrate(real_total):.4f} e   (2 O atoms x 6 valence e- each)')
print(f'  Correct net spin (up-down):{integrate(real_spin):.4f} e   (matches SIESTA\\'s own reported |S| = 2.0)')
"
echo
echo "stb-density now always uses the correct (total/net-spin) combination -- see below."
pause

echo "=================================================================="
echo " output/basic-slice/  --  Sn3O4 total charge density (non-magnetic)"
echo "=================================================================="
cat <<'EOF'
A plain, non-spin-polarized calculation -- [1] reports the total charge
only, integrating to 72 e (the full valence electron count for this
14-atom Sn3O4 cell):
EOF
mkdir -p "$OUT/basic-slice"
cp Sn3O4.RHO "$OUT/basic-slice/"
echo
echo "\$ stb-density --label Sn3O4 --save-gnuplot --no-intro"
(cd "$OUT/basic-slice" && stb-density --label Sn3O4 --save-gnuplot --no-intro > console.log 2>&1)
awk '/\[1\] CHARGE/{flag=1} /\[4\] REFERENCES/{flag=0} flag' "$OUT/basic-slice/console.log"
pause

echo "=================================================================="
echo " output/spin-auto/  --  O2's spin-polarization auto-detected, correctly"
echo "=================================================================="
cat <<'EOF'
No --spin flag needed: [0] shows "Spin-polarized: yes (auto-detected)",
[1] reports the correct total (12 e), and [2] reports the correct net
spin (2 e) -- matching the hand-computed "Correct" values above, and
SIESTA's own reported |S| = 2.0 for this real O2 triplet ground state:
EOF
mkdir -p "$OUT/spin-auto"
cp o2.RHO o2.XV o2.fdf "$OUT/spin-auto/"
echo
echo "\$ stb-density --label o2 --save-gnuplot --no-intro"
(cd "$OUT/spin-auto" && stb-density --label o2 --save-gnuplot --no-intro > console.log 2>&1)
awk '/\[0\] RUN METADATA/{flag=1} /\[1\] CHARGE/{flag=0} flag' "$OUT/spin-auto/console.log"
awk '/\[1\] CHARGE/{flag=1} /\[4\] REFERENCES/{flag=0} flag' "$OUT/spin-auto/console.log"
pause

echo "=================================================================="
echo " output/rho2-diff/  --  --rho2: a charge-transfer difference"
echo "=================================================================="
cat <<'EOF'
--rho2 subtracts a second .RHO (Delta rho = rho1 - rho2) -- the standard
way to visualize charge transfer/bonding. Subtracting Sn3O4.RHO from
itself is trivially zero everywhere, but it's enough to prove the real
mechanics: a signed quantity gets a colorbar/palette symmetric around
zero (Blue -> White -> Red), not anchored at 0 like the non-negative
total charge above:
EOF
mkdir -p "$OUT/rho2-diff"
cp Sn3O4.RHO "$OUT/rho2-diff/Sn3O4.RHO"
cp Sn3O4.RHO "$OUT/rho2-diff/Sn3O4_copy.RHO"
echo
echo "\$ stb-density --label Sn3O4 --rho2 Sn3O4_copy.RHO --save-gnuplot --no-intro"
(cd "$OUT/rho2-diff" && stb-density --label Sn3O4 --rho2 Sn3O4_copy.RHO --save-gnuplot --no-intro > console.log 2>&1)
grep "Integrated Delta\|Blue -> White" "$OUT/rho2-diff/console.log" "$OUT/rho2-diff/Sn3O4_density.gplot"
pause

echo "=================================================================="
echo " output/full-report/  --  --save-report / --save-gnuplot (both off by default)"
echo "=================================================================="
cat <<'EOF'
Every run always prints the numbered [0]...[5] report to the console.
Without --save-report/--save-gnuplot, a plain run only ever writes the
.dat file(s) and references.bib -- no .gplot script, no text report file:
EOF
mkdir -p "$OUT/full-report"
cp Sn3O4.RHO "$OUT/full-report/"
echo
echo "\$ stb-density --label Sn3O4 --no-intro   # default: no .gplot, no report"
(cd "$OUT/full-report" && stb-density --label Sn3O4 --no-intro > console_default.log 2>&1)
ls "$OUT/full-report/" | grep -v "\.RHO$\|console"
echo "(no .gplot, no stb_density_report.txt -- only Sn3O4_density.dat + references.bib)"
echo
echo "\$ stb-density --label Sn3O4 --save-report --save-gnuplot --no-intro"
(cd "$OUT/full-report" && stb-density --label Sn3O4 --save-report --save-gnuplot --no-intro > console_saved.log 2>&1)
echo "Report sections written to stb_density_report.txt:"
grep -E "^\[[0-9]+\]" "$OUT/full-report/stb_density_report.txt"
echo
echo "references.bib -- SIESTA (every stb-density run analyzes a SIESTA structure file):"
grep "^@" "$OUT/full-report/references.bib"
pause

echo "=================================================================="
echo " Two ways to run it"
echo "=================================================================="
cat <<'EOF'
A -- direct CLI:
  stb-density --label Sn3O4

B -- interactive stb-suite menu:
  $ stb-suite
  Select an option (0-6, or a tool code like 4.1.2): 3.8

Both paths call the exact same underlying tool -- proven directly below.
EOF
TMP="$(mktemp -d)"
cp Sn3O4.RHO "$TMP/"
echo
echo "\$ printf '3.8\\nSn3O4\\n1\\n2\\n\\nn\\nn\\n\\nn\\n\\nn\\nn\\nn\\n\\n0\\n' | stb-suite     # 2D slice, axis 2, center, defaults"
(cd "$TMP" && printf '3.8\nSn3O4\n1\n2\n\nn\nn\n\nn\n\nn\nn\nn\n\n0\n' | stb-suite > session.log 2>&1) || true
CLI_LINE=$(grep "Integrated Charge Density" "$OUT/basic-slice/console.log" | head -1)
MENU_LINE=$(grep "Integrated Charge Density" "$TMP/session.log" | head -1)
if [ "$CLI_LINE" = "$MENU_LINE" ]; then
    echo "Confirmed: identical integrated-charge line from the CLI and the interactive menu."
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
  basic-slice/   spin-auto/
  rho2-diff/     full-report/

Each has references.bib; full-report/ additionally has
stb_density_report.txt and a .gplot script (only from its --save-report/
--save-gnuplot run).

Recap of what this walkthrough covered:
  - a live, hand-verified proof of the up/down-vs-total/spin fix, using a
    real spin-polarized O2 calculation
  - auto-detection: no --spin flag needed to notice/process a
    spin-polarized .RHO correctly
  - the non-negative-vs-signed colorbar convention (--rho2 diff)
  - the numbered [0]...[5] report, --save-report, --save-gnuplot,
    references.bib
  - CLI and the interactive stb-suite menu building the same command

As a next step, try on your own with a finished SIESTA calculation:
  stb-density --label my_calc --save-gnuplot --view
  stb-density --label my_spin_calc --view   # spin detected automatically
EOF
