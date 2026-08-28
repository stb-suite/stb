#!/bin/bash
# Guided example: Workflow 4.10 -- 2D Stacking Fault / Gamma-Surface (both
# stages: stb-stackingfault, code 4.10.1, and stb-stackingfaultAnalysis,
# code 4.10.2, in the stb-suite menu).
#
# Not an automated test (see test/4-workflow/10-stackingfault/{prep,analysis}/
# test.sh for that) -- a commented walk-through: it runs real commands, one
# case at a time, into its own output/<case>/ folder, and shows you the
# piece of output that proves what just happened. Pauses between sections
# so you can read before moving on. Safe to re-run any time -- it always
# starts by wiping its own output/.
#
# graphene.fdf/hbn.fdf are the same 2-atom primitive-cell monolayers used by
# the automated test fixture (a=2.46 Ang / a=2.504 Ang, ~1.8% lattice
# mismatch -- realistic enough to exercise the ZSL commensurate-supercell
# search for real, Case 3 below).
#
# Stage 1 runs for real everywhere below (cheap -- it only writes .fdf
# files plus, for --mode 1/2, a real z-only relaxation: --mode 1 needs a
# real SIESTA binary, so those grid points are only DEMONSTRATED via their
# generated config_extra.fdf, not actually relaxed here; --mode 2 needs the
# optional 'ml' extra and IS actually run for real, MACE-MP-0 is cheap
# enough for a handful of 4-atom single points). Stage 2 is demonstrated
# with a real physical energy landscape: Case 5 fabricates calc.out using
# the EXACT SIESTA FreeEng values from a real bilayer-graphene 1D scan
# (30 points, a genuine DFT calculation, PAO.BasisSize DZP/GGA-PBE) --
# authentic numbers, just written onto this folder's own tiny fixture
# instead of re-running real SIESTA here.

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

HAS_ML=0
if python3 -c "import mace" 2>/dev/null; then HAS_ML=1; fi
# --mode 2 + --d3 (default ON) asks MACE-MP-0 for D3(BJ) dispersion, which
# itself needs the separate optional 'torch-dftd' package -- checked here
# so Case 2 degrades gracefully (same spirit as the HAS_ML check above)
# instead of a raw traceback on a machine that has 'ml' but not this.
HAS_TORCH_DFTD=0
if python3 -c "import torch_dftd" 2>/dev/null; then HAS_TORCH_DFTD=1; fi

# The 30-point 1D FreeEng progression below is REAL: siesta: FreeEng values
# from an actual bilayer-graphene DFT relaxation-based --scan x sweep
# (PAO.BasisSize DZP, XC GGA-PBE, real SIESTA 5.4.2 runs) -- not invented
# for this walkthrough. It's written onto THIS folder's own small
# graphene.fdf/graphene.fdf pair in Case 5 below (a different, unrelaxed
# a=2.46 Ang lattice than the real run's own relaxed one, so the meV/Ang^2
# figure Stage 2 reports here is this fixture's own, not the original
# run's -- the point is the SHAPE of a genuine gamma-surface energy
# progression, not reproducing one specific paper's number).
REAL_FREEENG=(
  -654.900928 -654.901136 -654.901574 -654.902258 -654.903136
  -654.904736 -654.906246 -654.907901 -654.909639 -654.911400
  -654.913058 -654.914659 -654.915805 -654.916873 -654.917482
  -654.917635 -654.917482 -654.916873 -654.915805 -654.914659
  -654.913058 -654.911400 -654.909639 -654.907901 -654.906246
  -654.904736 -654.903136 -654.902258 -654.901574 -654.901136
)


echo "=================================================================="
echo " Welcome: what a gamma-surface is, and why no BSSE stage is needed"
echo "=================================================================="
cat <<'EOF'
Slide one 2D layer of a van der Waals bilayer laterally over the other,
recompute the total energy at every offset, and you get a 2D energy
landscape -- the generalized stacking-fault energy surface, or "gamma
surface". Its two headline numbers: the EQUILIBRIUM stacking (the
lowest-energy offset -- e.g. graphene's AB registry) and the CORRUGATION
(max - min over the sampled grid -- the energy barrier resisting
interlayer sliding/shear, directly relevant to friction and stacking-
fault energetics).

Stage 1 (stb-stackingfault) rigidly slides layer 2 across a grid of
lateral offsets and writes one SIESTA folder per point. Stage 2
(stb-stackingfaultAnalysis) reads every finished energy back and reports
the equilibrium stacking, corrugation, and the full gamma-surface map.

This workflow used to have a 3rd stage (BSSE/counterpoise correction).
It was REMOVED: every point in a stacking-fault scan has the exact SAME
atoms and the exact SAME basis set (only positions move) -- unlike
adsorption or cohesive-energy calculations, where the compared systems
genuinely have different atom counts/basis sizes. Since each isolated
layer's own energy is a CONSTANT across the whole scan (it doesn't
depend on the other layer's lateral position at all), that constant
cancels exactly in any max-min comparison -- the RAW total energy's
corrugation is already, exactly, the correct physical answer. See the
README's Section 1 and 5 for the full derivation (it's short and worth
reading once).
EOF
pause


echo "=================================================================="
echo " Case 1: Stage 1, a real 2D surface scan (homobilayer graphene/graphene)"
echo "=================================================================="
cat <<'EOF'
Passing the SAME file as --layer1/--layer2 is the canonical use case: a
material sliding against itself (e.g. graphite ABA vs. ABC stacking).
--mode 3 (fixed gap, plain single-point) is the cheapest of the 3 gap
strategies -- Case 2 below compares all 3. A 5x5 grid is small enough to
run in seconds; --d3 stays ON (default) since interlayer binding is van
der Waals-dominated, which plain GGA misses.
EOF
mkdir -p "$OUT/case1-surface"
cp graphene.fdf calc.fdf "$OUT/case1-surface/"
echo "\$ stb-stackingfault -l1 graphene.fdf -l2 graphene.fdf -c calc.fdf -nx 5 -ny 5 --mode 3 --no-intro"
(cd "$OUT/case1-surface" && stb-stackingfault -l1 graphene.fdf -l2 graphene.fdf -c calc.fdf \
    -nx 5 -ny 5 --mode 3 --no-intro | sed -n '/\[1\] ZSL MATCH/,/\[3\] SUMMARY/p' | head -n -1)
echo
echo "25 single-point SIESTA folders, ready to run:"
ls "$OUT/case1-surface/sf_run/positions" | head -5
echo "... ($(ls "$OUT/case1-surface/sf_run/positions" | wc -l) total)"
pause


echo "=================================================================="
echo " Case 2: the 3 interlayer-gap strategies, side by side (--mode 1/2/3)"
echo "=================================================================="
cat <<'EOF'
Mode 1 (SIESTA relaxes z for real) writes a restricted CG relaxation
(x,y frozen for every atom, only z free) via config_extra.fdf -- most
accurate, most expensive. Mode 2 (MACE-MP-0 relaxes z first) is the
middle ground: a cheap ML z-relaxation, then a plain SIESTA single-point
at the result. Mode 3 (fixed gap) never relaxes anything -- cheapest,
today's Case 1 above. Same 2-point --scan x sweep, 3 different modes:
watch each one's own config_extra.fdf.
EOF
for mode in 1 2 3; do
    mkdir -p "$OUT/case2-modes/mode$mode"
    cp graphene.fdf calc.fdf "$OUT/case2-modes/mode$mode/"
done
echo "\$ stb-stackingfault -l1 graphene.fdf -l2 graphene.fdf -c calc.fdf --scan x -n 2 --mode 1 --no-intro"
(cd "$OUT/case2-modes/mode1" && stb-stackingfault -l1 graphene.fdf -l2 graphene.fdf -c calc.fdf \
    --scan x -n 2 --mode 1 --no-intro > prep.log 2>&1)
echo "--- mode 1's config_extra.fdf (restricted CG, x/y frozen) ---"
cat "$OUT/case2-modes/mode1/sf_run/positions/shift_00_00/config_extra.fdf"
if [ "$HAS_ML" -eq 1 ] && [ "$HAS_TORCH_DFTD" -eq 1 ]; then
    echo "\$ stb-stackingfault -l1 graphene.fdf -l2 graphene.fdf -c calc.fdf --scan x -n 2 --mode 2 --no-intro"
    (cd "$OUT/case2-modes/mode2" && stb-stackingfault -l1 graphene.fdf -l2 graphene.fdf -c calc.fdf \
        --scan x -n 2 --mode 2 --no-intro > prep.log 2>&1)
    echo "--- mode 2's config_extra.fdf (MACE already relaxed z -- single-point only) ---"
    cat "$OUT/case2-modes/mode2/sf_run/positions/shift_00_00/config_extra.fdf"
elif [ "$HAS_ML" -eq 1 ]; then
    echo "Optional 'torch-dftd' package not installed (needed for MACE-MP-0's own D3(BJ)"
    echo "dispersion, requested whenever --d3 is on, the default) -- running mode 2 with"
    echo "--no-d3 instead, just to show the mechanism (its config_extra.fdf below has no"
    echo "DFTD3 line, unlike modes 1/3's -- that's this workaround, not mode 2 itself)."
    echo "\$ stb-stackingfault -l1 graphene.fdf -l2 graphene.fdf -c calc.fdf --scan x -n 2 --mode 2 --no-d3 --no-intro"
    (cd "$OUT/case2-modes/mode2" && stb-stackingfault -l1 graphene.fdf -l2 graphene.fdf -c calc.fdf \
        --scan x -n 2 --mode 2 --no-d3 --no-intro > prep.log 2>&1)
    echo "--- mode 2's config_extra.fdf (MACE already relaxed z -- single-point only) ---"
    cat "$OUT/case2-modes/mode2/sf_run/positions/shift_00_00/config_extra.fdf"
else
    echo "Optional 'ml' extra not installed (pip install stb_suite[ml]) -- skipping mode 2's"
    echo "live demonstration. It writes the SAME single-point block as mode 3 below, just at a"
    echo "MACE-MP-0-relaxed z instead of the fixed --gap."
fi
echo "\$ stb-stackingfault -l1 graphene.fdf -l2 graphene.fdf -c calc.fdf --scan x -n 2 --mode 3 --no-intro"
(cd "$OUT/case2-modes/mode3" && stb-stackingfault -l1 graphene.fdf -l2 graphene.fdf -c calc.fdf \
    --scan x -n 2 --mode 3 --no-intro > prep.log 2>&1)
echo "--- mode 3's config_extra.fdf (fixed gap, plain single-point) ---"
cat "$OUT/case2-modes/mode3/sf_run/positions/shift_00_00/config_extra.fdf"
pause


echo "=================================================================="
echo " Case 3: a real heterostructure (graphene/hBN) -- the ZSL match, live"
echo "=================================================================="
cat <<'EOF'
Different --layer1/--layer2 files study an interlayer/heterostructure
sliding landscape instead of a material sliding against itself.
graphene.fdf (a=2.46 Ang) and hbn.fdf (a=2.504 Ang) are ~1.8% lattice
-mismatched -- stb-stackingfault finds the lowest-strain COMMENSURATE
supercell (Zur & McGill's ZSL algorithm) that fits both, unlike Case 1's
trivial 1x1 (identical layers, 0 strain by construction).
EOF
mkdir -p "$OUT/case3-heterostructure"
cp graphene.fdf hbn.fdf calc.fdf "$OUT/case3-heterostructure/"
echo "\$ stb-stackingfault -l1 graphene.fdf -l2 hbn.fdf -c calc.fdf -nx 3 -ny 3 --mode 3 --no-intro"
(cd "$OUT/case3-heterostructure" && stb-stackingfault -l1 graphene.fdf -l2 hbn.fdf -c calc.fdf \
    -nx 3 -ny 3 --mode 3 --no-intro | sed -n '/\[1\] ZSL MATCH/,/\[2\] GRID FOLDERS/p' | head -n -1)
pause


echo "=================================================================="
echo " Case 4: 1D scans (--scan x/y/xy) -- a cheap profile before the full grid"
echo "=================================================================="
cat <<'EOF'
--scan x/y/xy trade the full 2D gamma-surface for a cheaper 1D line
through it -- along shift_x, shift_y, or the shift_x=shift_y diagonal.
Same folder-naming/manifest scheme as --scan surface (shift_II_JJ, one
index held at 0) -- Stage 2 needs no changes to read either shape, it
detects which one it's looking at from sf_manifest.json.
EOF
for axis in x y xy; do
    mkdir -p "$OUT/case4-scans/$axis"
    cp graphene.fdf calc.fdf "$OUT/case4-scans/$axis/"
    echo "\$ stb-stackingfault -l1 graphene.fdf -l2 graphene.fdf -c calc.fdf --scan $axis -n 5 --mode 3 --no-intro"
    (cd "$OUT/case4-scans/$axis" && stb-stackingfault -l1 graphene.fdf -l2 graphene.fdf -c calc.fdf \
        --scan "$axis" -n 5 --mode 3 --no-intro > prep.log 2>&1)
    n=$(find "$OUT/case4-scans/$axis/sf_run/positions" -maxdepth 1 -type d -name 'shift_*' | wc -l)
    echo "  -> $n folder(s): $(ls "$OUT/case4-scans/$axis/sf_run/positions" | tr '\n' ' ')"
done
pause


echo "=================================================================="
echo " Case 5: Stage 2 with a REAL 30-point energy landscape (stb-stackingfaultAnalysis, 4.10.2)"
echo "=================================================================="
cat <<'EOF'
Fabricating calc.out for all 30 points of a --scan x sweep, using REAL
siesta: FreeEng values from an actual bilayer-graphene DFT calculation
(see this script's own REAL_FREEENG array, top of file) -- a genuine
gamma-surface energy progression, not invented numbers. Running the full
Stage 2 report, with --save-gnuplot (writes into 'plot/'), --view
(matplotlib preview, headless-safe here via MPLBACKEND=Agg), and
--view-animation (opens the always-saved extended-XYZ animation in ASE's
viewer -- also headless-safe, DISPLAY is unset in this environment so it
fails gracefully with a [FAIL], same as it would over a plain SSH
session with no X11).
EOF
mkdir -p "$OUT/case5-analysis"
cp graphene.fdf calc.fdf "$OUT/case5-analysis/"
echo "\$ stb-stackingfault -l1 graphene.fdf -l2 graphene.fdf -c calc.fdf --scan x -n 30 --mode 3 --no-intro"
(cd "$OUT/case5-analysis" && stb-stackingfault -l1 graphene.fdf -l2 graphene.fdf -c calc.fdf \
    --scan x -n 30 --mode 3 --no-intro > prep.log 2>&1)
i=0
for label in $(cd "$OUT/case5-analysis/sf_run/positions" && ls -d shift_*_00 | sort); do
    e="${REAL_FREEENG[$i]}"
    printf "siesta: FreeEng =    %s\nSCF cycle converged after 10 iterations\nsiesta: Atomic forces (eV/Ang):\n   Max    0.010000\n" \
        "$e" > "$OUT/case5-analysis/sf_run/positions/$label/calc.out"
    i=$((i+1))
done
echo "\$ stb-stackingfaultAnalysis --dir sf_run --save-gnuplot --view --view-animation --no-intro"
DISPLAY= stb-stackingfaultAnalysis --dir "$OUT/case5-analysis/sf_run" --save-gnuplot --view \
    --view-animation --no-intro > "$OUT/case5-analysis.log" 2>&1
sed -n '/\[2\] STACKING FAULT ANALYSIS/,/\[3\] SUMMARY/p' "$OUT/case5-analysis.log" | head -n -1
echo
echo "First/last few rows of the [1] GRID ENERGIES table (meV/Ang^2, relative to the minimum):"
sed -n '/\[1\] GRID ENERGIES/,/\[2\] STACKING/p' "$OUT/case5-analysis.log" | head -8
echo "  ... (30 rows total, symmetric around shift_x=0.5 -- see the full table in "
echo "      output/case5-analysis.log or the README's worked example)"
echo
echo "Files written:"
ls "$OUT/case5-analysis/sf_run/plot/"
ls "$OUT/case5-analysis/sf_run/stackingfault_animation.xyz"
pause


echo "=================================================================="
echo " Case 6: --apply -- promoting the equilibrium point to a production file"
echo "=================================================================="
cat <<'EOF'
Copies the lowest-energy grid point's structure.fdf out of sf_run/ to a
path of your choosing -- e.g. as the starting geometry for a follow-up
calculation at that specific registry.
EOF
echo "\$ stb-stackingfaultAnalysis --dir sf_run --apply equilibrium.fdf --no-intro"
(cd "$OUT/case5-analysis" && stb-stackingfaultAnalysis --dir sf_run --apply equilibrium.fdf \
    --no-intro | sed -n '/\[4\] APPLY/,$p')
ls -la "$OUT/case5-analysis/equilibrium.fdf"
pause


echo "=================================================================="
echo " Case 7: CLI vs. the interactive stb-suite menu (4.10.1) -- same result"
echo "=================================================================="
cat <<'EOF'
The interactive menu asks the same questions instead of flags, then
calls the exact same stb-stackingfault underneath. Reproducing Case 1's
own grid (5x5, --mode 3) through stb-suite -> 4.10.1 and diffing against
it proves the 2 paths are equivalent.
EOF
mkdir -p "$OUT/case7-interactive"
cp graphene.fdf calc.fdf "$OUT/case7-interactive/"
# Prompts in order: layer1 (graphene.fdf), layer2 (graphene.fdf), calc.fdf,
# mode (3), scan shape (blank -> 1=surface), pseudopotential source (blank
# -> skip), grid Nx (5), grid Ny (5), gap (blank -> 3.2), D3 (blank -> Y),
# ML pre-relax (blank -> N), ML preview (blank -> N), advanced settings
# (blank -> N), save report (blank -> N), "Press Enter to continue", quit.
echo "\$ printf '4.10.1\\ngraphene.fdf\\ngraphene.fdf\\ncalc.fdf\\n3\\n\\n\\n5\\n5\\n\\n\\n\\n\\n\\n\\n\\n0\\n' | stb-suite"
(cd "$OUT/case7-interactive" && \
    printf '4.10.1\ngraphene.fdf\ngraphene.fdf\ncalc.fdf\n3\n\n\n5\n5\n\n\n\n\n\n\n\n0\n' \
    | stb-suite > menu.log 2>&1)
if diff -rq "$OUT/case1-surface/sf_run/positions" "$OUT/case7-interactive/sf_run/positions" > /dev/null 2>&1; then
    echo "CLI (Case 1) and interactive menu (4.10.1) produced byte-identical grid folders."
else
    echo "UNEXPECTED: CLI and interactive menu results differ -- see"
    echo "  diff -rq '$OUT/case1-surface/sf_run/positions' '$OUT/case7-interactive/sf_run/positions'"
fi
pause


echo "=================================================================="
echo " Workflow 4.10 complete"
echo "=================================================================="
cat <<'EOF'
Stage 1 (stb-stackingfault, 4.10.1): one SIESTA folder per grid point
(--scan surface/x/y/xy), 3 interlayer-gap strategies, real ZSL matching
for heterostructures, D3 dispersion by default, sf_manifest.json written
for Stage 2. Stage 2 (stb-stackingfaultAnalysis, 4.10.2): equilibrium
stacking + corrugation in meV/Ang^2, optional gnuplot data (--save-gnuplot,
under 'plot/'), an on-screen matplotlib preview (--view), an always-saved
extended-XYZ animation of the whole grid + optional ASE 3D viewer
(--view-animation), and --apply to promote the equilibrium geometry.

See the README's worked-example section for the full 30-point table this
script's Case 5 generated, and its "Known limitations" section for what
to watch out for with --mode 1 (a loose MD.MaxForceTol can let the CG
relaxation stop early at SOME grid points but not others, introducing a
kink in the corrugation curve that isn't present in the true relaxed
landscape -- verified live in the conversation that produced this tool).
EOF
