#!/bin/bash
# Guided example: Workflow 4.1 -- Stress-Strain, Stage 1 (stb-strain, code
# 4.1.1 in the stb-suite menu). Stage 2 (stb-strainAnalysis, code 4.1.2)
# will be added to this same script in a follow-up update.
#
# Not an automated test (see test/4-workflow/1-strain/prep/test.sh for
# that) -- a commented walk-through: it runs real commands, one case at a
# time, into its own output/<case>/ folder, and shows you the piece of
# output that proves what just happened. Pauses between sections so you
# can read before moving on. Safe to re-run any time -- it always starts
# by wiping its own output/.
#
# structure.fdf is bulk silicon, an 8-atom conventional cubic cell.
# calc.fdf is the (correct) calc.fdf for its original relaxation --
# MD.VariableCell true, MD.Steps 100. calc_missing_relax.fdf is the same
# file with both deliberately wrong, used only for the "gotcha" case
# below -- never a real template.

set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

OUT="$DIR/output"
rm -rf "$OUT"
mkdir -p "$OUT"

pause() {
    echo
    read -p "  [Press Enter to continue] " -r
    echo
}

echo "=================================================================="
echo " What Workflow 4.1 -- Stage 1 (stb-strain) does"
echo "=================================================================="
cat <<'EOF'
Applies a small Cartesian strain (uniaxial x/y/z, or biaxial xy/xz/yz) to
a relaxed structure's lattice vectors, for a range of strain values, and
writes one ready-to-run SIESTA folder per value: the deformed structure,
a copy of calc.fdf, a new config_extra.fdf (see below), and any linked
pseudopotentials. Doesn't run SIESTA itself -- you run each folder
yourself, then (once added here) hand the results to Stage 2
(stb-strainAnalysis) for the stress-strain fit.

--relax-mode selects how the cell responds to the imposed strain:
  cell-fixed          cell locked exactly at the imposed strain; only
                       ions relax
  stress-constrained   only the imposed direction's own stress fixed;
                       every other periodic direction relaxes freely
Both are expressed through the same mechanism: a new config_extra.fdf
file with a %block Geometry.Constraints, %include-d into calc.fdf.
EOF
pause


echo "=================================================================="
echo " Case 1: --relax-mode cell-fixed"
echo "=================================================================="
cat <<'EOF'
calc.fdf is the SAME file used for structure.fdf's original relaxation
(same basis/k-points/XC/SCF) -- stb-strain only strains the lattice and
adds one %include line, it never regenerates the rest. Watch the [4]
report section: it reads calc.fdf's current MD.TypeOfRun/MD.Steps/
MD.VariableCell (never forcing any of them) and shows exactly which Voigt
stress components config_extra.fdf will fix.
EOF
mkdir -p "$OUT/cell-fixed"
cp structure.fdf calc.fdf "$OUT/cell-fixed/"
echo "\$ stb-strain -s structure.fdf -c calc.fdf --relax-mode cell-fixed --stdir x --stmin 0 --stmax 2 --step 2 --no-intro"
(cd "$OUT/cell-fixed" && stb-strain -s structure.fdf -c calc.fdf --relax-mode cell-fixed \
    --stdir x --stmin 0 --stmax 2 --step 2 --no-intro | tee console.log \
    | sed -n '/\[4\] RELAXATION/,/\[5\] PSEUDOPOTENTIALS/p')
echo
echo "config_extra.fdf actually written (all 6 Voigt components fixed):"
cat "$OUT/cell-fixed/strain_runs/x/strain_x_2.00/config_extra.fdf"
echo
echo "calc.fdf's own MD block, preserved exactly as given (nothing forced):"
grep -E "^MD\." "$OUT/cell-fixed/strain_runs/x/strain_x_2.00/calc.fdf"
echo "...and the only line ADDED to it:"
grep "%include" "$OUT/cell-fixed/strain_runs/x/strain_x_2.00/calc.fdf"
pause


echo "=================================================================="
echo " Case 2: --relax-mode stress-constrained -- same inputs, different mode"
echo "=================================================================="
cat <<'EOF'
Same structure.fdf/calc.fdf, same direction -- only --relax-mode changes.
config_extra.fdf now fixes only the imposed direction's own component.
EOF
mkdir -p "$OUT/stress-constrained"
cp structure.fdf calc.fdf "$OUT/stress-constrained/"
echo "\$ stb-strain -s structure.fdf -c calc.fdf --relax-mode stress-constrained --stdir x --stmin 0 --stmax 2 --step 2 --no-intro"
(cd "$OUT/stress-constrained" && stb-strain -s structure.fdf -c calc.fdf --relax-mode stress-constrained \
    --stdir x --stmin 0 --stmax 2 --step 2 --no-intro > console.log \
    && sed -n '/Component | Label/,/^$/p' console.log)
echo
echo "config_extra.fdf actually written (only the imposed direction fixed):"
cat "$OUT/stress-constrained/strain_runs/x/strain_x_2.00/config_extra.fdf"
pause


echo "=================================================================="
echo " Case 3: the gotcha -- a misconfigured calc.fdf, caught live"
echo "=================================================================="
cat <<'EOF'
calc_missing_relax.fdf is the same template with MD.VariableCell removed
and MD.Steps set to 0 -- exactly the 2 conditions the [4] report section
checks for. stb-strain never forces these (they live entirely in YOUR
calc.fdf), so it warns instead of silently producing folders that would
never actually relax anything.
EOF
mkdir -p "$OUT/gotcha"
cp structure.fdf calc_missing_relax.fdf "$OUT/gotcha/"
echo "\$ stb-strain -s structure.fdf -c calc_missing_relax.fdf --relax-mode cell-fixed --stdir x --stmin 0 --stmax 1 --no-intro"
(cd "$OUT/gotcha" && stb-strain -s structure.fdf -c calc_missing_relax.fdf --relax-mode cell-fixed \
    --stdir x --stmin 0 --stmax 1 --no-intro 2>&1 | grep "WARNING\]")
pause


echo "=================================================================="
echo " Case 4: the axis-symmetry advisory -- cubic Si, x/y/z all equivalent"
echo "=================================================================="
cat <<'EOF'
structure.fdf is cubic (point group m-3m) -- one of the few point groups
where every Cartesian axis is truly mechanically equivalent. The [3]
report section detects this automatically, before you waste a real DFT
run straining y or z too.
EOF
mkdir -p "$OUT/symmetry"
cp structure.fdf calc.fdf "$OUT/symmetry/"
echo "\$ stb-strain -s structure.fdf -c calc.fdf --relax-mode cell-fixed --stdir x --stmin 0 --stmax 1 --no-intro"
(cd "$OUT/symmetry" && stb-strain -s structure.fdf -c calc.fdf --relax-mode cell-fixed \
    --stdir x --stmin 0 --stmax 1 --no-intro | sed -n '/\[3\] AXIS SYMMETRY/,/\[4\] RELAXATION/p')
pause


echo "=================================================================="
echo " Case 5: 2 directions, same --output-dir -- one subfolder each"
echo "=================================================================="
cat <<'EOF'
Run x, then y, into the same strain_runs/ -- each direction gets its own
<output-dir>/<direction>/ subfolder, so nothing collides. This is also
what lets Stage 2 (once added here) read every direction found under
strain_runs/ at once.
EOF
mkdir -p "$OUT/multi-direction"
cp structure.fdf calc.fdf "$OUT/multi-direction/"
(cd "$OUT/multi-direction" && \
    stb-strain -s structure.fdf -c calc.fdf --relax-mode cell-fixed --stdir x --stmin 0 --stmax 2 --step 2 --no-intro > /dev/null && \
    stb-strain -s structure.fdf -c calc.fdf --relax-mode cell-fixed --stdir y --stmin 0 --stmax 2 --step 2 --no-intro > /dev/null)
echo "\$ find strain_runs -maxdepth 2 -type d | sort"
(cd "$OUT/multi-direction" && find strain_runs -maxdepth 2 -type d | sort)
pause


echo "=================================================================="
echo " Case 6: CLI vs. the interactive stb-suite menu (4.1.1) -- same result"
echo "=================================================================="
cat <<'EOF'
The interactive menu asks the same questions instead of flags, then calls
the exact same stb-strain underneath. Reproducing Case 1 (cell-fixed,
direction x) through stb-suite -> 4.1.1 and diffing against Case 1's own
output proves the 2 paths are equivalent.
EOF
mkdir -p "$OUT/interactive"
cp structure.fdf calc.fdf "$OUT/interactive/"
# Prompts in order: structure file (blank -> default structure.fdf), calc
# file (blank -> default calc.fdf), relax mode (1 -> cell-fixed), pseudo
# source (blank -> skip), advanced settings (n -> skip), direction (x --
# y,z ARE equivalent on this cubic fixture, so the "generate them too?"
# follow-up appears -- answer n to match Case 1's single direction),
# stmin, stmax, step, save-report (n), "Press Enter to continue", quit.
(cd "$OUT/interactive" && \
    printf '4.1.1\n\n\n1\n\nn\nx\nn\n0\n2\n2\nn\n\n0\n' | stb-suite > menu.log 2>&1)
if diff -rq "$OUT/cell-fixed/strain_runs/x" "$OUT/interactive/strain_runs/x" > /dev/null 2>&1; then
    echo "CLI (Case 1) and interactive menu (4.1.1) produced byte-identical folders."
else
    echo "UNEXPECTED: CLI and interactive menu results differ -- see"
    echo "  diff -rq '$OUT/cell-fixed/strain_runs/x' '$OUT/interactive/strain_runs/x'"
fi
pause


echo "=================================================================="
echo " Stage 1 recap"
echo "=================================================================="
cat <<'EOF'
Stage 1 (stb-strain, 4.1.1) is done: strain_runs/<direction>/ folders are
ready, each with a deformed structure, calc.fdf (with %include
config_extra.fdf added), config_extra.fdf itself, and any linked
pseudopotentials. Run SIESTA in each one.

Stage 2 (stb-strainAnalysis, 4.1.2) -- reading the finished SIESTA runs
back and fitting the stress-strain curve -- will be added to this same
script and README.md in a follow-up update.
EOF
