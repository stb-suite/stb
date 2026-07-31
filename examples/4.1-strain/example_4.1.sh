#!/bin/bash
# Guided example: Workflow 4.1 -- Stress-Strain (both stages: stb-strain,
# code 4.1.1, and stb-strainAnalysis, code 4.1.2, in the stb-suite menu).
#
# Not an automated test (see test/4-workflow/1-strain/{prep,analysis}/
# test.sh for that) -- a commented walk-through: it runs real commands,
# one case at a time, into its own output/<case>/ folder, and shows you
# the piece of output that proves what just happened. Pauses between
# sections so you can read before moving on. Safe to re-run any time -- it
# always starts by wiping its own output/.
#
# structure.fdf is bulk silicon, an 8-atom conventional cubic cell.
# calc.fdf is the (correct) calc.fdf for its original relaxation --
# MD.VariableCell true, MD.Steps 100. calc_missing_relax.fdf is the same
# file with both deliberately wrong, used only for the "gotcha" case
# below -- never a real template.

set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

# Stage 2's --view calls plt.show() -- MPLBACKEND=Agg makes that a no-op
# instead of blocking on a GUI window, same convention test.sh itself uses.
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
echo " Welcome: stress, strain, and what this workflow measures"
echo "=================================================================="
cat <<'EOF'
A small crystal deformation is F = I + strain_tensor (the standard small
-strain/linear-elasticity approximation) -- Stage 1 (stb-strain) applies
this to the lattice vectors for a range of strain values and writes one
ready-to-run SIESTA folder per value; Stage 2 (stb-strainAnalysis) reads
the finished runs back and fits the resulting stress-strain curve.

That curve's shape carries 4 physical numbers (a 5th is opt-in):
  Initial Slope   stiffness -- stress needed for a small strain (Hooke's
                  law), NOT automatically a full elastic tensor C_ij (that
                  needs the sibling stb-elasticAnalysis workflow instead)
  Peak (UTS)      the largest stress/force reached -- ultimate strength
  Critical Strain the strain at which that peak occurs
  Toughness       area under the whole curve -- energy absorbed before
                  failure, a distinct concept from stiffness or strength
  Yield (0.2%)    opt-in; a macroscopic metallurgy concept, only loosely
                  meaningful for a defect-free periodic DFT crystal

Units depend on dimensionality, auto-detected by Stage 2: 3D -> GPa (a
real stress); 2D -> N/m (SIESTA's own stress un-diluted by the arbitrary
vacuum height, the same convention 2D-materials literature uses, e.g.
graphene's ~340 N/m in-plane stiffness); 1D -> nN (a raw axial force, no
arbitrary cross-section assumed).

--relax-mode selects how the cell responds to the imposed strain:
  cell-fixed          cell locked exactly at the imposed strain; only
                       ions relax
  stress-constrained   only the imposed direction's own stress fixed;
                       every other periodic direction relaxes freely
Both are expressed through the same mechanism: a new config_extra.fdf
file with a %block Geometry.Constraints, %include-d into calc.fdf. See
this folder's README.md for the full theory (section 1) behind all of
this.
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
# source (blank -> skip), advanced settings (n -> skip), direction menu
# (1 -- 'x', the first of 5 entries on this cubic structure: x, xy, ALL
# UNIAXIAL, ALL BIAXIAL, ALL -- see Case 7 below), stmin, stmax, step,
# save-report (n), "Press Enter to continue", quit.
(cd "$OUT/interactive" && \
    printf '4.1.1\n\n\n1\n\nn\n1\n0\n2\n2\nn\n\n0\n' | stb-suite > menu.log 2>&1)
if diff -rq "$OUT/cell-fixed/strain_runs/x" "$OUT/interactive/strain_runs/x" > /dev/null 2>&1; then
    echo "CLI (Case 1) and interactive menu (4.1.1) produced byte-identical folders."
else
    echo "UNEXPECTED: CLI and interactive menu results differ -- see"
    echo "  diff -rq '$OUT/cell-fixed/strain_runs/x' '$OUT/interactive/strain_runs/x'"
fi
pause


echo "=================================================================="
echo " Case 7: the interactive menu's direction-picker -- pick ALL, get x AND xy"
echo "=================================================================="
cat <<'EOF'
The interactive menu (4.1.1) detects vacuum axes and symmetry-equivalent
directions (uniaxial AND biaxial) up front, and only offers independent,
non-redundant ones instead of the free-text prompt the CLI itself has.
On this cubic structure.fdf that's just 2 out of the 6 possible
directions: x (y, z are equivalent) and xy (xz, yz are equivalent).
Picking the menu's "ALL" entry runs Stage 1 once per independent
direction in a single pass -- "ALL UNIAXIAL"/"ALL BIAXIAL" do the same
scoped to just one strain type.
EOF
mkdir -p "$OUT/all-directions"
cp structure.fdf calc.fdf "$OUT/all-directions/"
# Prompts in order: structure file (blank -> default), calc file (blank ->
# default), relax mode (1 -> cell-fixed), pseudo source (blank -> skip),
# advanced settings (n -> skip), direction menu (5 -- 'ALL independent
# directions -- uniaxial + biaxial', the last of 5 entries: x, xy, ALL
# UNIAXIAL, ALL BIAXIAL, ALL), stmin, stmax, step, save-report (n), "Press
# Enter to continue", quit.
(cd "$OUT/all-directions" && \
    printf '4.1.1\n\n\n1\n\nn\n5\n0\n2\n2\nn\n\n0\n' | stb-suite > menu.log 2>&1)
echo "\$ find strain_runs -maxdepth 2 -type d | sort"
(cd "$OUT/all-directions" && find strain_runs -maxdepth 2 -type d | sort)
pause


echo "=================================================================="
echo " Stage 1 recap"
echo "=================================================================="
cat <<'EOF'
Stage 1 (stb-strain, 4.1.1) is done: strain_runs/<direction>/ folders are
ready, each with a deformed structure, calc.fdf (with %include
config_extra.fdf added), config_extra.fdf itself, and any linked
pseudopotentials. Run SIESTA in each one, then hand the results to Stage 2
(stb-strainAnalysis, 4.1.2) below.
EOF
pause


echo "=================================================================="
echo " Case 8: Stage 2 -- reading one direction's results back (stb-strainAnalysis, 4.1.2)"
echo "=================================================================="
cat <<'EOF'
This script doesn't run real SIESTA, so Stage 2 is demonstrated here
against small, hand-built calc.out data instead (same idea as this
folder's own calc_missing_relax.fdf -- illustrative, not a real result).
Stage 2 reports in the same numbered [N] section style as Stage 1's own
report, and auto-detects dimensionality from a real structure.fdf if one
is found in the folder (none is here, so it falls back to 3D with a
[NOTE], not an error).

This fixture also deliberately uses a coarse 10%/step (linear only up to
~20%, then plateauing -- a stand-in for real yielding) so only 1 point
falls within the usual +-2% window, showing the fit-window WIDEN to the
nearest 3 points instead of silently averaging in the plateau -- watch
"Fit window" and the [NOTE] below.
EOF
mkdir -p "$OUT/stage2-single/strain_runs/x/strain_x_0.00" \
         "$OUT/stage2-single/strain_runs/x/strain_x_10.00" \
         "$OUT/stage2-single/strain_runs/x/strain_x_20.00" \
         "$OUT/stage2-single/strain_runs/x/strain_x_30.00"
for pair in "0.00:0" "10.00:15" "20.00:28" "30.00:29"; do
    pct="${pair%%:*}"; kbar="${pair##*:}"
    cat > "$OUT/stage2-single/strain_runs/x/strain_x_$pct/calc.out" <<EOF
outcell: Unit cell vectors (Ang):
        5.470000    0.000000    0.000000
        0.000000    5.470000    0.000000
        0.000000    0.000000    5.470000

siesta: Stress tensor Voigt (kbar):     ${kbar}.00      0.00      0.00      0.00      0.00      0.00
EOF
done
echo "\$ stb-strainAnalysis --file calc.out --dir strain_runs --save-report --save-gnuplot --no-intro"
(cd "$OUT/stage2-single" && stb-strainAnalysis --file calc.out --dir strain_runs \
    --save-report --save-gnuplot --no-intro | sed -n '/\[2\] DIMENSIONALITY/,/\[4\] OUTPUT FILES/p')
pause


echo "=================================================================="
echo " Case 9: Stage 2 -- multiple directions compared automatically, no flag needed"
echo "=================================================================="
cat <<'EOF'
Point --dir at the top-level output directory (holding one <direction>/
subfolder per direction, e.g. from Case 7's "ALL" run) and every direction
found is compared automatically -- no --compare-style flag exists or is
needed. --view shows the same comparison as an interactive matplotlib plot
instead of (or alongside) --save-gnuplot's .dat+.gplot pair.
EOF
mkdir -p "$OUT/stage2-compare/strain_runs/x/strain_x_0.00" \
         "$OUT/stage2-compare/strain_runs/x/strain_x_1.00" \
         "$OUT/stage2-compare/strain_runs/x/strain_x_2.00" \
         "$OUT/stage2-compare/strain_runs/xy/strain_xy_0.00" \
         "$OUT/stage2-compare/strain_runs/xy/strain_xy_1.00" \
         "$OUT/stage2-compare/strain_runs/xy/strain_xy_2.00"
for pair in "0.00:0" "1.00:15" "2.00:28"; do
    pct="${pair%%:*}"; kbar="${pair##*:}"
    cat > "$OUT/stage2-compare/strain_runs/x/strain_x_$pct/calc.out" <<EOF
outcell: Unit cell vectors (Ang):
        5.470000    0.000000    0.000000
        0.000000    5.470000    0.000000
        0.000000    0.000000    5.470000

siesta: Stress tensor Voigt (kbar):     ${kbar}.00      0.00      0.00      0.00      0.00      0.00
EOF
done
for pair in "0.00:0" "1.00:22" "2.00:40"; do
    pct="${pair%%:*}"; kbar="${pair##*:}"
    cat > "$OUT/stage2-compare/strain_runs/xy/strain_xy_$pct/calc.out" <<EOF
outcell: Unit cell vectors (Ang):
        5.470000    0.000000    0.000000
        0.000000    5.470000    0.000000
        0.000000    0.000000    5.470000

siesta: Stress tensor Voigt (kbar):     ${kbar}.00      ${kbar}.00      0.00      0.00      0.00      0.00
EOF
done
echo "\$ stb-strainAnalysis --file calc.out --dir strain_runs --no-intro --view"
(cd "$OUT/stage2-compare" && stb-strainAnalysis --file calc.out --dir strain_runs \
    --no-intro --view | sed -n '/\[3\] MECHANICAL/,/\[4\] OUTPUT FILES/p')
echo "(--view opened/closed a matplotlib comparison plot -- MPLBACKEND=Agg above"
echo " makes that a silent no-op in this unattended script; drop it yourself to"
echo " actually see the plot pop up)"
pause


echo "=================================================================="
echo " Workflow 4.1 complete"
echo "=================================================================="
cat <<'EOF'
Stage 1 (stb-strain, 4.1.1): strained structures generated, one folder per
direction/strain value. Stage 2 (stb-strainAnalysis, 4.1.2): reads real
SIESTA runs back (auto-detecting dimensionality, auto-comparing multiple
directions), reports mechanical properties in a numbered [0]-[5] section
style, and can plot with gnuplot (--save-gnuplot) and/or matplotlib
(--view).
EOF
