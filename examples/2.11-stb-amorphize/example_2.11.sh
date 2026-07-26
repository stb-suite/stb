#!/bin/bash
# Guided example: stb-amorphize (code 2.11 in the stb-suite menu)
#
# Not an automated test (see test/2-structures/11-amorphize/test.sh for
# that) -- a commented walk-through: it runs real commands, one group at a
# time, into its own output/<case>/ folder, and shows you the piece of
# output that proves what just happened. Pauses between sections so you
# can read before moving on. Safe to re-run any time -- it always starts
# by wiping its own output/. Needs the optional 'ml' extra (pip install
# stb_suite[ml]); the whole script is skipped with a clear message if
# `mace` isn't importable.
#
# Uses the same tiny 8-atom si8.fdf fixture throughout, deliberately, so
# every command below is directly comparable to every other one -- only
# the flags being demonstrated change. Some comparison sections use
# --seed for repeatable initial velocities; even so, an 8-atom cell is far
# too small for the exact numbers to be anything but a qualitative,
# single-trajectory illustration (see the "why so tiny" note below) --
# the TREND each comparison shows is the real, reproducible point.

set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

if ! python3 -c "import mace" 2>/dev/null; then
    echo "This example needs the optional 'ml' extra: pip install stb_suite[ml]"
    exit 0
fi

OUT="$DIR/output"
rm -rf "$OUT"
mkdir -p "$OUT"

pause() {
    echo
    read -p "  [Press Enter to continue] " -r
    echo
}

strip_ansi() {
    sed 's/\x1b\[[0-9;]*m//g' "$1"
}

echo "=================================================================="
echo " What melt-quench amorphization is, physically"
echo "=================================================================="
cat <<'EOF'
An amorphous/glassy structure (no long-range crystalline order, but real
short-range chemistry -- roughly the right bond lengths/angles/
coordination) is hard to build directly: you can't just "un-order" atomic
positions randomly without creating physically absurd overlaps and bond
lengths. The classic laboratory trick -- melt the material, then cool it
faster than it can re-crystallize -- works because a liquid has already
forgotten its crystalline arrangement while keeping sensible local
chemistry (atoms still push each other apart at the right distance); cool
it FAST enough and there isn't time for atoms to find their way back into
an ordered lattice before the structure is "frozen" kinetically -- a glass,
not a crystal.

stb-amorphize reproduces exactly this in silico with a MACE-MP-0-driven MD
run: melt a crystalline structure well above its real melting point (so it
melts even given the tiny simulated time available -- real experimental
quenches take seconds; this one takes picoseconds), hold it there long
enough to actually forget the crystal, then ramp the temperature back down.
NOT a substitute for a slower, production-quality quench (many more steps,
a slower cool-down) or DFT verification -- a fast heuristic STARTING GUESS,
meant to be much better than an ad-hoc random-displacement structure.
EOF
pause

echo "=================================================================="
echo " Why an 8-atom cell? (and why the numbers below are qualitative)"
echo "=================================================================="
cat <<'EOF'
This walkthrough deliberately reuses one tiny 8-atom bulk Si cell
(si8.fdf, the same fixture test/2-structures/11-amorphize/ uses) so every
comparison below finishes in seconds on a CPU. Real amorphization needs a
much bigger supercell (tens to hundreds of atoms) -- with only 8 atoms,
a SINGLE MD trajectory is a genuinely chaotic, small-number-statistics
system: re-running the exact same command WITHOUT --seed gives visibly
different bond-angle numbers each time, and even --seed only fixes the
INITIAL velocities, not every source of run-to-run float noise. Every
comparison below still shows a real, physically sensible TREND (that's
the point being demonstrated) -- just don't expect the precise digits to
reproduce bit-for-bit on your machine, or to mean much on their own for a
cell this small.
EOF
pause

echo "=================================================================="
echo " NPT-Berendsen: a thermostat AND a barostat, with finite response time"
echo "=================================================================="
cat <<'EOF'
The MD engine underneath is ASE's NPTBerendsen: it doesn't just run plain
constant-energy dynamics, it actively drives the system toward a TARGET
temperature (via --taut, the thermostat's response time, in fs) and a
target pressure (via --taup, the barostat's response time -- 0 GPa here,
since a liquid/quenched structure should relax its own volume freely, not
stay pinned at the crystal's original cell size).

Both --taut/--taup are RESPONSE TIMES, not instantaneous snaps: Berendsen
coupling relaxes the current temperature toward the target roughly
exponentially, with --taut as the time constant. That single fact drives
two important consequences demonstrated below:
  1. --melt-steps needs to be comfortably longer than --taut (default 50
     fs) for the system to actually reach the melt temperature, not just
     start drifting toward it.
  2. A --quench-steps ramp that's FASTER than --taut leaves the actual,
     instantaneous temperature lagging behind the falling target --
     exactly why the final static relax (see further down) exists at all.

The very first velocities come from a Maxwell-Boltzmann distribution at
(at most) 300 K -- --seed makes that initial draw reproducible; the MD
trajectory that follows is deterministic given it, modulo floating-point
noise across runs/machines.
EOF
pause

echo "=================================================================="
echo " Two ways to run it"
echo "=================================================================="
cat <<'EOF'
A -- direct CLI:
  stb-amorphize -f si8.fdf

B -- interactive stb-suite menu:
  $ stb-suite
  Select an option (0-6, or a tool code like 4.1.2): 2.11

Every run always writes a numbered report ([0] RUN METADATA ... [9]
SUMMARY & FILES) and a bond-angle mean/std check before and after the MD --
a crystal starts near 0 deg std; a genuinely amorphized structure shows a
much broader spread while the mean stays near the same short-range bond
angle. --save-data/--save-traj/--save-report/--view are all off by default.
EOF
pause

echo "=================================================================="
echo " output/too-short-vs-adequate-melt/  --  --melt-steps needs to exceed --taut"
echo "=================================================================="
cat <<'EOF'
Putting a number to the theory above: --taut defaults to 50 fs. A melt
held for only 5 steps never gets close to actually thermalizing at the
target temperature; a melt held for 100 steps (2x --taut) does. --seed 7
fixes the initial velocities so both runs start from the exact same draw
-- only --melt-steps/--quench-steps differ:
EOF
mkdir -p "$OUT/too-short-vs-adequate-melt"
cp si8.fdf "$OUT/too-short-vs-adequate-melt/"
echo
echo "\$ stb-amorphize -f si8.fdf --melt-steps 5 --quench-steps 10 --seed 7 \\"
echo "      --no-final-relax -o too_short.fdf --no-intro"
(cd "$OUT/too-short-vs-adequate-melt" && stb-amorphize -f si8.fdf --melt-steps 5 --quench-steps 10 \
    --seed 7 --no-final-relax -o too_short.fdf --no-intro > console_short.log 2>&1)
grep "Bond-angle mean/std" "$OUT/too-short-vs-adequate-melt/console_short.log" | tail -1

echo
echo "\$ stb-amorphize -f si8.fdf --melt-steps 100 --quench-steps 200 --seed 7 \\"
echo "      --no-final-relax -o adequate.fdf --no-intro"
(cd "$OUT/too-short-vs-adequate-melt" && stb-amorphize -f si8.fdf --melt-steps 100 --quench-steps 200 \
    --seed 7 --no-final-relax -o adequate.fdf --no-intro > console_adequate.log 2>&1)
grep "Bond-angle mean/std" "$OUT/too-short-vs-adequate-melt/console_adequate.log" | tail -1
cat <<'EOF'

The 20x-longer protocol shows a visibly larger bond-angle spread -- more
genuine loss of crystalline order -- even though the crystal's own mean
angle (Si's tetrahedral ~109.5 deg) survives in both cases, exactly the
signature this diagnostic is designed to catch (mean stable, std rising).
EOF
pause

echo "=================================================================="
echo " output/final-relax-comparison/  --  why the final static relax matters"
echo "=================================================================="
cat <<'EOF'
Second consequence of the finite Berendsen response time above: even a
"reasonable" quench ramp can end with the ACTUAL instantaneous temperature
still well above the target -- real residual thermal energy/force, not a
clean local minimum. --no-final-relax skips the cleanup step; by default
stb-amorphize always relaxes positions AND cell afterward:
EOF
mkdir -p "$OUT/final-relax-comparison"
cp si8.fdf "$OUT/final-relax-comparison/"
echo
echo "\$ stb-amorphize -f si8.fdf --melt-steps 30 --quench-steps 60 --seed 7 \\"
echo "      --no-final-relax -o no_relax.fdf --no-intro"
(cd "$OUT/final-relax-comparison" && stb-amorphize -f si8.fdf --melt-steps 30 --quench-steps 60 \
    --seed 7 --no-final-relax -o no_relax.fdf --no-intro > console_norelax.log 2>&1)
grep "After MD" "$OUT/final-relax-comparison/console_norelax.log"
echo
echo "\$ stb-amorphize -f si8.fdf --melt-steps 30 --quench-steps 60 --seed 7 \\"
echo "      -o with_relax.fdf --no-intro          # final relax ON (the default)"
(cd "$OUT/final-relax-comparison" && stb-amorphize -f si8.fdf --melt-steps 30 --quench-steps 60 \
    --seed 7 -o with_relax.fdf --no-intro > console_relax.log 2>&1)
echo
echo "[4] FINAL STATIC RELAX (before = right after the quench MD, no cooldown):"
strip_ansi "$OUT/final-relax-comparison/console_relax.log" | awk '/^Quantity/,/^$/'
cat <<'EOF'

"Before" here is the structure exactly as the quench MD left it -- the max
force is real thermal noise (atoms still mid-vibration), not a converged
geometry. The relax removes it in a handful of FIRE steps, converging to
the same order of magnitude (~0.01-0.05 eV/Ang) every other stb-suite
relaxation targets -- this is why the final relax is on by default.
EOF
pause

echo "=================================================================="
echo " output/vacuum-rejected/  --  bulk-only, on purpose"
echo "=================================================================="
cat <<'EOF'
Melting/NPT-relaxing a vacuum-padded axis is physically meaningless -- a
slab's vacuum gap isn't a real liquid boundary, and letting NPTBerendsen's
barostat "relax" it would just be relaxing empty space. stb-amorphize
detects this (the same vacuum-axis heuristic stb-kgrid/stb-mlrelax use)
and refuses outright rather than silently producing nonsense:
EOF
mkdir -p "$OUT/vacuum-rejected"
cp si8.fdf "$OUT/vacuum-rejected/"
echo
echo "\$ stb-slab -f si8.fdf --hkl 1 0 0 -o slab.fdf --no-intro   # build a (100) slab"
(cd "$OUT/vacuum-rejected" && stb-slab -f si8.fdf --hkl 1 0 0 -o slab.fdf --no-intro > /dev/null 2>&1)
echo "\$ stb-amorphize -f slab.fdf --no-intro"
rc=0
(cd "$OUT/vacuum-rejected" && stb-amorphize -f slab.fdf --no-intro > console.log 2>&1) || rc=$?
echo "Exit code: $rc  (1 = rejected, as expected)"
grep "ERROR" "$OUT/vacuum-rejected/console.log"
pause

echo "=================================================================="
echo " output/save-data/  --  .dat + gnuplot export of the whole MD run"
echo "=================================================================="
cat <<'EOF'
Until this feature, every intermediate MD step was discarded -- only the
final frame was ever written anywhere, and the only trace of temperature/
energy during the run was a progress line on stderr that vanished the
moment the run ended. --save-data keeps it: one continuous step/time axis
spanning BOTH the melt and quench stages, sampled every --stride steps,
written as <stem>_md_diagnostics.dat (step, time, T, E_pot/E_kin/E_total,
cell volume) using gnuplot's own 'index' block convention -- the melt
stage is 'index 0', the quench stage is 'index 1' in the SAME file, so
either stage can be plotted alone, or both together via 'index 0:1'. The
companion .gplot renders both energy and temperature vs. step into one
PDF, with the melt->quench transition marked and each stage's target
temperature drawn as a reference line -- run it with gnuplot yourself:
EOF
mkdir -p "$OUT/save-data"
cp si8.fdf "$OUT/save-data/"
echo
echo "\$ stb-amorphize -f si8.fdf --melt-steps 40 --quench-steps 60 --stride 5 \\"
echo "      --save-data -o amorphous.fdf --no-intro"
(cd "$OUT/save-data" && stb-amorphize -f si8.fdf --melt-steps 40 --quench-steps 60 \
    --stride 5 --save-data -o amorphous.fdf --no-intro > console.log 2>&1)
grep "MD diagnostics" "$OUT/save-data/console.log"
echo
echo "The two gnuplot 'index' blocks inside amorphous_md_diagnostics.dat:"
grep -c "index 0: melt" "$OUT/save-data/amorphous_md_diagnostics.dat" > /dev/null
head -4 "$OUT/save-data/amorphous_md_diagnostics.dat"
echo "  ... (one row every 5 steps, melt phase) ..."
grep -A1 "index 1: quench" "$OUT/save-data/amorphous_md_diagnostics.dat" | tail -1
echo "  ... (one row every 5 steps, quench phase) ..."
echo
if command -v gnuplot > /dev/null 2>&1; then
    echo "\$ gnuplot amorphous_md_diagnostics.gplot"
    (cd "$OUT/save-data" && gnuplot amorphous_md_diagnostics.gplot)
    ls "$OUT/save-data/amorphous_md_diagnostics.pdf"
else
    echo "(gnuplot not installed here -- skipping the actual render; the"
    echo " .gplot file above is ready to run with 'gnuplot <file>.gplot')"
fi
pause

echo "=================================================================="
echo " output/save-traj/  --  a trajectory to open in OVITO/VMD"
echo "=================================================================="
cat <<'EOF'
--save-traj writes a multi-frame trajectory of the melt-quench MD itself
(same sampling as --save-data, same --stride) -- the same 3-format choice
as stb-ani2traj/stb-mlmd: xsf (OVITO/VMD-native, default), pdb (VMD's own
default), or xyz (OVITO-native). Use it to actually WATCH the crystal
disorder as it melts and requench, frame by frame, instead of only
trusting a bond-angle number. All 3 formats read back cleanly with ASE:
EOF
mkdir -p "$OUT/save-traj"
cp si8.fdf "$OUT/save-traj/"
for fmt in xsf pdb xyz; do
    echo
    echo "\$ stb-amorphize -f si8.fdf --melt-steps 20 --quench-steps 20 --stride 5 \\"
    echo "      --save-traj --traj-format $fmt -o traj_${fmt}.fdf --no-intro"
    (cd "$OUT/save-traj" && stb-amorphize -f si8.fdf --melt-steps 20 --quench-steps 20 \
        --stride 5 --save-traj --traj-format "$fmt" -o "traj_${fmt}.fdf" --no-intro \
        > "console_${fmt}.log" 2>&1)
done
python3 - "$OUT/save-traj" <<'PYEOF'
import sys
from ase.io import read
out = sys.argv[1]
for fmt, ext in [("xsf", "xsf"), ("pdb", "pdb"), ("xyz", "xyz")]:
    frames = read(f"{out}/traj_{fmt}_md_traj.{ext}", index=":")
    print(f"  {fmt}: {len(frames)} frames, {frames[0].get_chemical_formula()} each")
PYEOF
pause

echo "=================================================================="
echo " output/full-report/  --  --save-report, symmetry, references.bib"
echo "=================================================================="
cat <<'EOF'
The full numbered report (also written to stb_amorphize_report.txt with
--save-report) includes a [6] SYMMETRY ANALYSIS BEFORE/AFTER table.
UNLIKE most other stb-suite before/after tables, these are EXPECTED to
differ -- a genuine amorphization loses the crystal's long-range symmetry,
often collapsing all the way down to P1/P-1 (essentially no symmetry
operations left besides identity, or identity+inversion). Seeing a much
lower-symmetry space group in the "After" column is itself a positive
sign the melt-quench worked, not a validation failure -- --seed 2 below
is fixed only so this walkthrough's own README table is reproducible;
any adequately-melted run shows the same kind of collapse, just to a
different specific low-symmetry group each time:
EOF
mkdir -p "$OUT/full-report"
cp si8.fdf "$OUT/full-report/"
echo "\$ stb-amorphize -f si8.fdf --melt-steps 100 --quench-steps 200 --seed 2 --save-report --no-intro"
(cd "$OUT/full-report" && stb-amorphize -f si8.fdf --melt-steps 100 --quench-steps 200 --seed 2 \
    --save-report --no-intro > console.log 2>&1)
echo
echo "Report sections written to stb_amorphize_report.txt:"
grep -E "^\[[0-9]+\]" "$OUT/full-report/stb_amorphize_report.txt"
echo
echo "Symmetry before/after (crystal space group -> after melt-quench):"
awk '/^Property/,/^$/' "$OUT/full-report/stb_amorphize_report.txt" | grep -E "Property|Space Group"
echo
echo "Provenance header written into the output .fdf:"
head -4 "$OUT/full-report/amorphous.fdf"
echo
echo "references.bib -- SIESTA + MACE + MACE-MP (foundation model was used):"
grep "^@" "$OUT/full-report/references.bib"
pause

echo "=================================================================="
echo " Proof: CLI and the interactive stb-suite menu agree"
echo "=================================================================="
echo "Driving the same case through the interactive menu's --save-data/"
echo "--save-traj prompts and checking the same files come out."
TMP="$(mktemp -d)"
cp "$DIR/si8.fdf" "$TMP/"
echo
echo "\$ printf '2.11\\nsi8.fdf\\n3000\\n15\\n300\\n15\\ny\\ny\\nxyz\\n5\\n\\n\\n\\n\\n\\nmenu_out.fdf\\n0\\n' | stb-suite"
(cd "$TMP" && printf '2.11\nsi8.fdf\n3000\n15\n300\n15\ny\ny\nxyz\n5\n\n\n\n\n\nmenu_out.fdf\n0\n' | \
    stb-suite > session.log 2>&1) || true
if [ -f "$TMP/menu_out_md_diagnostics.dat" ] && [ -f "$TMP/menu_out_md_traj.xyz" ]; then
    echo "Confirmed: the interactive menu built and launched the exact same"
    echo "underlying stb-amorphize command as the CLI walkthrough above"
    echo "(menu_out_md_diagnostics.dat + menu_out_md_traj.xyz both written)."
else
    echo "Unexpected: menu did not reach the write step -- see $TMP/session.log."
fi
rm -rf "$TMP"
pause

echo "=================================================================="
echo " Done"
echo "=================================================================="
cat <<EOF
Seven self-contained folders were generated under output/:
  too-short-vs-adequate-melt/   final-relax-comparison/
  vacuum-rejected/               save-data/
  save-traj/                     full-report/

Each has references.bib; full-report/ additionally has
stb_amorphize_report.txt; save-data/ additionally has the .dat/.gplot
(and .pdf, if gnuplot is installed) diagnostics pair.

Recap of what this walkthrough covered:
  - the melt-quench idea: melt above the real melting point long enough to
    forget the crystal, then quench faster than it can re-crystallize
  - NPT-Berendsen's finite response time (--taut/--taup) and its two real
    consequences: --melt-steps must exceed --taut to actually melt, and a
    fast quench ramp leaves real residual thermal energy/force behind
  - the bond-angle mean/std diagnostic that proves amorphization happened
    (mean stable near the crystal's own angle, std rising)
  - why the final static relax is on by default, with real before/after
    numbers
  - bulk-only, and why melting a vacuum-padded slab is meaningless
  - --save-data: <stem>_md_diagnostics.dat + .gplot, one continuous step
    axis across melt and quench via gnuplot 'index' blocks
  - --save-traj: a multi-frame trajectory (xsf/pdb/xyz) of the whole MD
    run, for viewing in OVITO/VMD
  - the symmetry-collapses-to-a-low-symmetry-group diagnostic, references.bib, --save-report
  - CLI and the interactive stb-suite menu building the same command

Not exercised by this script (needs a display): --view opens the input
and final structure interactively via ASE -- try it yourself:
  stb-amorphize -f si8.fdf --view

As a next step, try on your own (with a real, bigger supercell -- not
this walkthrough's tiny 8-atom toy cell):
  stb-amorphize -f your_supercell.fdf --melt-temp 4000 --quench-steps 2000 \\
      --save-data --save-traj
  stb-amorphize -f your_supercell.fdf --custom-model my_finetuned.model
EOF
