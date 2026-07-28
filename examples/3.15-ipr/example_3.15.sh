#!/bin/bash
# Guided example: stb-ipr (code 3.15 in the stb-suite menu)
#
# Not an automated test (see test/3-analysis/15-ipr/test.sh for that) -- a
# commented walk-through: it runs real commands, one group at a time, into
# its own output/<case>/ folder, and shows you the piece of output that
# proves what just happened. Pauses between sections so you can read
# before moving on. Safe to re-run any time -- it always starts by wiping
# its own output/.
#
# Sn3O4.bands/Sn3O4.bands.WFSX/Sn3O4.HSX are the exact same real SIESTA
# fixture 3.10-stb-fatbands/ uses -- stb-ipr shares its entire .bands/
# .WFSX loading and HSX-accurate/fdf-approximate machinery with
# stb-fatbands (core/band_scatter.py was extracted from fatbands.py once
# stb-ipr became its second consumer). spin/Ospin.* is a real
# spin-polarized isolated-oxygen-atom calculation (2 Bohr-magneton
# triplet ground state), used for the nspin=2 section below.

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
echo " What IPR measures: Anderson localization"
echo "=================================================================="
cat <<'EOF'
IPR_q(k, band) = sum_i |psi_i|^(2q) / [sum_i |psi_i|^2]^q -- large for a
state whose weight is concentrated on a few orbitals (localized), small
for one spread evenly across the whole cell (extended). q=2 is the
standard choice; sisl's own ipr() hard-requires q >= 2.

Uses the EXACT SAME per-orbital quantity (norm2, hadamard projection)
stb-fatbands already uses for its own orbital weights -- this tool is
really just a different reduction of the same numbers, which is why it
shares stb-fatbands' entire .bands/.WFSX loading, cross-check, and
.HSX-accurate/.fdf-approximate machinery.

Every run prints a numbered [0]...[6] report. --save-report persists it;
--save-gnuplot writes ipr.dat + a real .gplot script (off by default --
this tool used to write both unconditionally); --view shows a matplotlib
preview (off by default -- this tool previously ALWAYS showed it, with
no way to skip at all).
EOF
pause


echo "=================================================================="
echo " output/basic-run/  --  the numbered report, deep vs. near-Fermi localization"
echo "=================================================================="
mkdir -p "$OUT/basic-run"
cp Sn3O4.bands Sn3O4.bands.WFSX Sn3O4.HSX "$OUT/basic-run/"
cat <<'EOF'
Watch [3] IPR ANALYSIS report the mean/min/max IPR, and which exact state
is the most localized (max IPR) vs. most extended (min IPR):
EOF
echo "\$ stb-ipr --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX --shift fermi --no-intro"
(cd "$OUT/basic-run" && stb-ipr --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX \
    --shift fermi --no-intro | tee console.log | sed -n '/\[3\] IPR ANALYSIS/,/\[4\] WRITING/p')
echo
echo "Testing the physical claim directly: deep (core-like) states are MORE localized"
echo "than states near the Fermi level:"
(cd "$OUT/basic-run" && stb-ipr --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX \
    --shift fermi --save-gnuplot --no-intro > /dev/null 2>&1 && python3 -c "
import numpy as np
data = np.loadtxt('ipr.dat')
energy, ipr = data[:,1], data[:,2]
deep = ipr[energy < -15]
near_fermi = ipr[(energy > -2) & (energy < 2)]
print(f'  Deep states (E < -15 eV):        mean IPR = {deep.mean():.4f}')
print(f'  Near-Fermi states (-2 to 2 eV):  mean IPR = {near_fermi.mean():.4f}')
")
pause


echo "=================================================================="
echo " The --q order parameter: higher q sharpens localized/extended contrast"
echo "=================================================================="
cat <<'EOF'
As q increases, every IPR value shrinks -- but the most-extended state's
value collapses toward 0 much faster than the most-localized state's
does. Watch the contrast sharpen across q=2, 4, 6:
EOF
mkdir -p "$OUT/q-sweep"
cp Sn3O4.bands Sn3O4.bands.WFSX Sn3O4.HSX "$OUT/q-sweep/"
for q in 2 4 6; do
    echo
    echo "\$ stb-ipr --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX --shift fermi --q $q --no-intro"
    (cd "$OUT/q-sweep" && stb-ipr --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX \
        --shift fermi --q "$q" --no-intro | grep "^ipr ")
done
pause


echo "=================================================================="
echo " Accuracy: real overlap matrix (.HSX) vs. the orthogonal-basis approximation"
echo "=================================================================="
cat <<'EOF'
Without a real .HSX, IPR falls back to an implicit-orthogonal
approximation (|c|^2, no true overlap matrix) -- watch the SAME states
give a genuinely different mean IPR:
EOF
mkdir -p "$OUT/accuracy"
cp Sn3O4.bands Sn3O4.bands.WFSX Sn3O4.HSX calc.fdf structure.fdf Sn.ion Sn.ion.xml O.ion O.ion.xml "$OUT/accuracy/"
(cd "$OUT/accuracy" && \
    stb-ipr --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX --shift fermi \
        --save-gnuplot --no-intro > /dev/null 2>&1 && mv ipr.dat ipr_hsx.dat && \
    stb-ipr --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --geometry-file calc.fdf --shift fermi \
        --save-gnuplot --no-intro > /dev/null 2>&1 && mv ipr.dat ipr_fdf.dat && \
    python3 -c "
import numpy as np
hsx = np.loadtxt('ipr_hsx.dat')
fdf = np.loadtxt('ipr_fdf.dat')
print(f'  Sn3O4.HSX (accurate, overlap-aware): mean IPR = {hsx[:,2].mean():.5f}')
print(f'  calc.fdf only (approximate):         mean IPR = {fdf[:,2].mean():.5f}')
")
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
cp Sn3O4.bands Sn3O4.bands.WFSX "$OUT/label-plus-hsx/"
cp Sn3O4.HSX "$OUT/label-plus-hsx/calc.HSX"
cp calc.fdf "$OUT/label-plus-hsx/Sn3O4.fdf"
echo "\$ stb-ipr --label Sn3O4 --hsx-file calc.HSX --shift fermi --no-intro"
(cd "$OUT/label-plus-hsx" && stb-ipr --label Sn3O4 --hsx-file calc.HSX --shift fermi \
    --no-intro > console_ok.log 2>&1)
grep "WFSX file\|Geometry/Hamiltonian source" "$OUT/label-plus-hsx/console_ok.log"
pause


echo "=================================================================="
echo " A real bug fixed: --q now validates sisl's own q >= 2 requirement"
echo "=================================================================="
cat <<'EOF'
sisl's EigenstateElectron.ipr() hard-asserts q >= 2 internally -- --q 1
used to crash with a raw, unfriendly AssertionError traceback. Now a
clean, actionable error:
EOF
mkdir -p "$OUT/q-validation"
cp Sn3O4.bands Sn3O4.bands.WFSX Sn3O4.HSX "$OUT/q-validation/"
echo "\$ stb-ipr --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX --shift fermi --q 1 --no-intro"
(cd "$OUT/q-validation" && stb-ipr --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX \
    --shift fermi --q 1 --no-intro > console.log 2>&1) || true
tail -1 "$OUT/q-validation/console.log"
pause


echo "=================================================================="
echo " A real bug fixed: nspin=2 spin channels used to be silently merged"
echo "=================================================================="
cat <<'EOF'
The exact same bug stb-fatbands already found and fixed this session,
inherited here too: both spin channels used to be dumped into ONE
flat, spin-blind series. Now split into ipr_up/ipr_down, written and
reported separately -- watch the two means genuinely differ:
EOF
mkdir -p "$OUT/nspin-fix"
cp spin/Ospin.bands spin/Ospin.bands.WFSX spin/Ospin.HSX "$OUT/nspin-fix/"
echo "\$ stb-ipr --label Ospin --shift fermi --save-gnuplot --no-intro"
(cd "$OUT/nspin-fix" && stb-ipr --label Ospin --shift fermi --save-gnuplot \
    --no-intro | sed -n '/\[3\] IPR ANALYSIS/,/\[4\] WRITING/p')
echo
echo "Files written:"
(cd "$OUT/nspin-fix" && ls ipr_up.dat ipr_down.dat)
echo "(no merged ipr.dat -- each spin channel gets its own file)"
pause


echo "=================================================================="
echo " output/full-report/  --  --save-report / --save-gnuplot (off by default)"
echo "=================================================================="
cat <<'EOF'
Every run always prints the numbered [0]...[6] report to the console.
--save-report additionally persists it to stb_ipr_report.txt, and
--save-gnuplot writes ipr.dat + a real .gplot script -- both off by
default, so a plain run only ever writes references.bib:
EOF
mkdir -p "$OUT/full-report"
cp Sn3O4.bands Sn3O4.bands.WFSX Sn3O4.HSX "$OUT/full-report/"
echo "\$ stb-ipr --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX --shift fermi --no-intro   # default"
(cd "$OUT/full-report" && stb-ipr --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX \
    --shift fermi --no-intro > console_default.log 2>&1)
(cd "$OUT/full-report" && ls)
echo "(only references.bib -- no ipr.dat/.gplot, no stb_ipr_report.txt)"
echo
echo "\$ stb-ipr --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX --shift fermi \\"
echo "      --save-report --save-gnuplot --no-intro"
(cd "$OUT/full-report" && stb-ipr --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX \
    --shift fermi --save-report --save-gnuplot --no-intro > console_saved.log 2>&1)
echo "Report sections written to stb_ipr_report.txt:"
grep -o '\[[0-9]\] [A-Za-z& ()-]*' "$OUT/full-report/stb_ipr_report.txt" | sort -u
echo
if command -v gnuplot >/dev/null 2>&1; then
    (cd "$OUT/full-report" && gnuplot ipr.gplot)
    echo "(rendered ipr.pdf with the real, installed gnuplot)"
fi
pause


echo "=================================================================="
echo " Two ways to run it"
echo "=================================================================="
cat <<EOF
A -- direct CLI:
  stb-ipr --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX --shift fermi

B -- interactive stb-suite menu:
  \$ stb-suite
  Select an option (0-6, or a tool code like 4.1.2): 3.15

Both paths call the exact same underlying tool -- proven directly below.
The menu asks for the label and the .HSX path SEPARATELY (default
suggested if <label>.HSX exists), and gained save-report/save-gnuplot/
view prompts.
EOF
TMP="$(mktemp -d)"
cp Sn3O4.bands Sn3O4.bands.WFSX Sn3O4.HSX "$TMP/"
cp calc.fdf "$TMP/Sn3O4.fdf"
echo
echo "\$ stb-ipr --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX --shift fermi --no-intro"
(cd "$TMP" && stb-ipr --file Sn3O4.bands --wfsx Sn3O4.bands.WFSX --hsx-file Sn3O4.HSX --shift fermi \
    --no-intro > cli.log 2>&1)
echo
echo "\$ printf '3.15\\nSn3O4\\n\\n3\\n\\n\\nn\\nn\\nn\\n\\n0\\n' | stb-suite"
(cd "$TMP" && printf '3.15\nSn3O4\n\n3\n\n\nn\nn\nn\n\n0\n' | timeout 60 stb-suite > session.log 2>&1) || true
CLI_LINE=$(grep "Best (indirect) gap" "$TMP/cli.log")
MENU_LINE=$(grep "Best (indirect) gap" "$TMP/session.log")
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
Seven self-contained folders were generated under output/:
  basic-run/       q-sweep/          accuracy/       label-plus-hsx/
  q-validation/    nspin-fix/        full-report/

Each has references.bib; full-report/ additionally has ipr.dat/ipr.gplot/
ipr.pdf and stb_ipr_report.txt (only from its --save-report/
--save-gnuplot run).

Recap of what this walkthrough covered:
  - IPR theory: sum_i |psi_i|^(2q), the Anderson-localization
    interpretation, and why it's the same norm2(hadamard) quantity
    stb-fatbands already uses for its own orbital weights
  - a real physics result: deep, core-like states are more localized
    than states near the Fermi level (verified: mean IPR ~0.107 vs
    ~0.052)
  - the --q order parameter sharpens the localized/extended contrast as
    it increases
  - the real .HSX-vs-fdf-only accuracy tradeoff, verified with a
    genuine ~12% difference in mean IPR on the exact same states
  - a real bug: --label + --hsx-file used to be rejected outright; fixed
    the same way as stb-wfdensity/stb-sts/stb-coop
  - a real bug: --q 1/0/negative used to crash with a raw sisl
    AssertionError; now a clean, actionable error
  - a real bug: nspin=2 spin channels used to be silently merged into
    one series (the same bug already found and fixed in stb-fatbands),
    verified live on a real spin-polarized oxygen atom: ipr_up mean
    0.9515 vs ipr_down mean 1.0870
  - the numbered [0]...[6] report, --save-report, --save-gnuplot (now
    opt-in, previously unconditional), --view (now opt-in, previously
    always shown with no way to skip)
  - CLI and the interactive stb-suite menu building the same command

As a next step, try on your own with a real SIESTA calculation:
  stb-ipr --label my_calc --shift fermi --q 4 --save-report --save-gnuplot --view
EOF
