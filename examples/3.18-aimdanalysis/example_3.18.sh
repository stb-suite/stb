#!/bin/bash
# Guided example: stb-aimdAnalysis (code 3.18 in the stb-suite menu)
#
# Not an automated test (see test/3-analysis/18-aimdanalysis/test.sh for
# that) -- a commented walk-through: it runs real commands, one group at
# a time, into its own output/<case>/ folder, and shows you the piece of
# output that proves what just happened. Pauses between sections so you
# can read before moving on. Safe to re-run any time -- it always starts
# by wiping its own output/.
#
# aimd.ANI/.XV/.fdf/.out are a REAL SIESTA AIMD run: a 5-step Verlet
# trajectory of an O2 dimer in vacuum. Too short for the MSD/VDOS numbers
# to carry real physical meaning, but its 2-atom composition makes
# --track-atom/--track-pair an easy, physically meaningful demo (the O-O
# bond stretching over the 5 steps).
#
# siesta.ANI/.XV/.MDE + calc.fdf/structure.fdf are a SECOND, real SIESTA
# run bundled here too: 500 MD steps, an 8-atom SiC supercell, Nose
# thermostat (NVT) at a 500 K target. SystemLabel 'siesta' but the real
# input file is 'calc.fdf' -- a genuine real-world --geometry-file
# mismatch, and a real siesta.MDE for the thermodynamic time series
# section further down.

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
echo " What stb-aimdAnalysis does"
echo "=================================================================="
cat <<'EOF'
Post-processes an AIMD/MD trajectory -- a SIESTA <label>.ANI (+ .out for
the per-step cell) or a generic ASE-readable trajectory (e.g. one written
by stb-mlmd) -- into five kinds of physics:

  [2] Radial distribution function g(r)         -- structure/bonding
  [3] Mean-squared displacement + diffusion D    -- how fast atoms move
  [4] Vibrational density of states (VDOS)       -- from the velocity
                                                     autocorrelation (VACF)
  [5]/[6] single-atom displacement tracking and atom-pair relative-
          distance tracking (--track-atom/--track-pair)
  [7] Thermodynamic time series: energy/temperature/volume/pressure,
      one 4-panel figure, straight from SIESTA's own <label>.MDE file

Every run prints a numbered [0]...[10] report. --save-report persists it;
--save-gnuplot writes a .dat + .gplot pair for every computed quantity
(off by default -- this tool used to write matplotlib PNGs unconditionally
on every run, with no gnuplot output at all); --view shows the same as an
interactive matplotlib preview (off by default -- this tool previously
always generated the PNGs, with no way to skip it).

Two more things worth knowing up front:
  --geometry-file <path>  the real SIESTA input .fdf is almost NEVER
                          named <label>.fdf in practice -- pass its real
                          path explicitly (needed for the true MD
                          timestep; a fixture below reproduces this bug
                          live)
  --list-atoms            prints every atom's index/species/coordinates
                          (first frame only, fast) then exits -- use this
                          to find valid --track-atom/--track-pair indices
                          instead of guessing
EOF
pause


echo "=================================================================="
echo " output/basic-run/  --  the numbered report, default output"
echo "=================================================================="
mkdir -p "$OUT/basic-run"
cp aimd.ANI aimd.fdf aimd.out aimd.XV "$OUT/basic-run/"
cat <<'EOF'
A plain run with no flags: the report prints to the console, but nothing
is written to disk except references.bib (this IS a real SIESTA run, so
the SIESTA citation is always written, unlike the --trajectory case
further down).
EOF
echo "\$ stb-aimdAnalysis --label aimd --no-intro"
(cd "$OUT/basic-run" && stb-aimdAnalysis --label aimd --no-intro \
    | tee console.log | sed -n '/\[1\] INPUT/,/\[5\] SINGLE-ATOM/p')
echo
echo "Files written:"
(cd "$OUT/basic-run" && ls)
echo "(only references.bib -- no .dat/.gplot/report without --save-gnuplot/--save-report)"
pause


echo "=================================================================="
echo " A real bug found via a user's real run: --label mode's .fdf lookup"
echo "=================================================================="
cat <<'EOF'
--label mode used to silently assume the real SIESTA input file is named
<label>.fdf -- needed for the true MD timestep (MD.LengthTimeStep). In
practice this almost never holds: SIESTA always names .XV/.ANI/.HSX/
.WFSX after SystemLabel, but the INPUT file's own name is chosen by the
user and is very often different (e.g. SystemLabel 'siesta' with the
real input called calc.fdf). This silently degraded to "assume 1 fs per
frame" with only an easy-to-miss warning. Fixed with --geometry-file,
the same --label-decoupled explicit path stb-sts/stb-coop/stb-ipr/
stb-effmass/stb-spintexture already use for their own .fdf/.HSX inputs:
EOF
mkdir -p "$OUT/geometry-file"
cp aimd.ANI aimd.out aimd.XV "$OUT/geometry-file/"
cp aimd.fdf "$OUT/geometry-file/calc.fdf"   # renamed -- simulates the real-world mismatch
echo "\$ stb-aimdAnalysis --label aimd --no-intro   # no aimd.fdf present -- only calc.fdf"
(cd "$OUT/geometry-file" && stb-aimdAnalysis --label aimd --no-intro 2>&1 \
    | grep "MD.LengthTimeStep")
echo
echo "\$ stb-aimdAnalysis --label aimd --geometry-file calc.fdf --no-intro"
(cd "$OUT/geometry-file" && stb-aimdAnalysis --label aimd --geometry-file calc.fdf --no-intro \
    2>&1 | grep "Geometry file\|Timestep")
pause


echo "=================================================================="
echo " NEW: --list-atoms -- find valid indices/coordinates before tracking anything"
echo "=================================================================="
cat <<'EOF'
--track-atom/--track-pair need 0-based atom indices -- --list-atoms prints
every atom's index, species, and Cartesian coordinates (first frame only,
so it's fast regardless of trajectory length), then exits immediately
(no RDF/MSD/VDOS/report at all). Off by default -- a real structure can
have hundreds of atoms, so this is opt-in, not part of every run's report.
The interactive stb-suite menu asks "List every atom's index/species/
coordinates?" (y/N) right before asking for --track-atom/--track-pair
(see section further down).
EOF
mkdir -p "$OUT/list-atoms"
cp aimd.ANI aimd.fdf aimd.out aimd.XV "$OUT/list-atoms/"
echo "\$ stb-aimdAnalysis --label aimd --list-atoms --no-intro"
(cd "$OUT/list-atoms" && stb-aimdAnalysis --label aimd --list-atoms --no-intro)
pause


echo "=================================================================="
echo " output/tracking/  --  NEW: --track-atom / --track-pair"
echo "=================================================================="
cat <<'EOF'
--track-atom N reports one atom's own Cartesian displacement from its
initial position (PBC-unwrapped, same convention as the MSD).

--track-pair I-J reports the relative distance between two SPECIFIC
atoms every frame, using the SAME minimum-image convention as the RDF --
deliberately NOT the unwrapped trajectory, because what matters for a
bond length is the atoms' TRUE instantaneous separation, not each atom's
independent drift path.

This fixture is a 2-atom O2 dimer, so atom 0 vs. atom 1 IS the O-O bond
-- watch the tracked distance below land in the same ~1.0-1.3 Ang range
as the RDF's own first peak (a real cross-check between two independent
code paths, same spirit as this session's other tools' verification):
EOF
mkdir -p "$OUT/tracking"
cp aimd.ANI aimd.fdf aimd.out aimd.XV "$OUT/tracking/"
echo "\$ stb-aimdAnalysis --label aimd --track-atom 0 --track-pair 0-1 --save-gnuplot --no-intro"
(cd "$OUT/tracking" && stb-aimdAnalysis --label aimd --track-atom 0 --track-pair 0-1 \
    --save-gnuplot --no-intro | tee console.log | sed -n '/\[5\] SINGLE-ATOM/,/\[7\] THERMODYNAMIC/p')
pause


echo "=================================================================="
echo " Cross-check: the tracked O-O distance and the RDF's first peak agree"
echo "=================================================================="
RDF_PEAK=$(grep "First peak" "$OUT/tracking/console.log" | grep -o 'r = [0-9.]*' | grep -o '[0-9.]*')
DIST_RANGE=$(python3 -c "
import numpy as np
d = np.loadtxt('$OUT/tracking/aimd_dist_0_1.dat')[:, 4]
print(f'{d.min():.3f} - {d.max():.3f}')
")
echo "RDF first peak       : $RDF_PEAK Ang"
echo "Tracked O-O distance : $DIST_RANGE Ang (min-max over the trajectory)"
echo "(both independently derived from the same trajectory -- the RDF from a"
echo " statistical histogram over all frames, the tracked pair from the exact"
echo " same two atoms' minimum-image separation frame by frame)"
if command -v gnuplot >/dev/null 2>&1; then
    (cd "$OUT/tracking" && gnuplot aimd_dist_0_1.gplot && gnuplot aimd_disp_atom0.gplot)
    echo "(rendered aimd_dist_0_1.pdf / aimd_disp_atom0.pdf with the real, installed gnuplot)"
fi
pause


echo "=================================================================="
echo " output/thermo/  --  NEW: [7] THERMODYNAMIC TIME SERIES, on a REAL 500-step NVT run"
echo "=================================================================="
cat <<'EOF'
siesta.ANI/.XV/.MDE + calc.fdf/structure.fdf are a SECOND, real SIESTA run
bundled with this example -- 500 MD steps, an 8-atom SiC supercell, Nose
thermostat (NVT) at a 500 K target. SystemLabel 'siesta' but the real
input is 'calc.fdf' -- the EXACT real-world mismatch --geometry-file
fixes (demonstrated above), hit live by a user running this tool for real.

[7] reads Energy/Temperature/Pressure straight from SIESTA's own small,
dedicated <label>.MDE file (Step/T/E_KS/E_tot/Vol/P) -- no .out log
-scraping needed for this part. Volume is ALWAYS available (from the
cell directly), regardless of input source.
EOF
mkdir -p "$OUT/thermo"
cp siesta.ANI siesta.XV siesta.MDE calc.fdf structure.fdf "$OUT/thermo/"
echo "\$ stb-aimdAnalysis --label siesta --geometry-file calc.fdf --save-gnuplot --no-intro"
(cd "$OUT/thermo" && stb-aimdAnalysis --label siesta --geometry-file calc.fdf --save-gnuplot \
    --no-intro | tee console.log | sed -n '/\[7\] THERMODYNAMIC/,/\[8\] WRITING/p')
pause


echo "=================================================================="
echo " A real physics lesson, live: E_tot vs. E_KS, and why it matters"
echo "=================================================================="
cat <<'EOF'
E_KS (SIESTA's own electronic/potential-like energy) trades with the
ions' kinetic energy by design -- NOT the quantity to judge stability
by. E_tot (kinetic+potential) is far more stable. Same lesson stb-mlmd
learned live this session for its own NVE energy tracking, now
confirmed here for a real SIESTA NVT run too. siesta_energy.dat's
columns 2/3 are the per-atom (eV/atom) values PLOTTED by both
siesta_energy.gplot and siesta_thermo.gplot -- energy is an EXTENSIVE
quantity (scales with atom count), so eV/atom (intensive, size-
independent) is the natural axis; columns 4/5 keep the raw absolute
eV totals for reference:
EOF
python3 -c "
import numpy as np
e_tot_pa, e_pot_pa, e_tot_abs, e_pot_abs = np.loadtxt(
    '$OUT/thermo/siesta_energy.dat', usecols=(1, 2, 3, 4), unpack=True)
t = np.loadtxt('$OUT/thermo/siesta_temperature.dat')[:, 1]
print(f'E_tot std : {e_tot_abs.std():.4f} eV  ({e_tot_pa.std():.6f} eV/atom)')
print(f'E_KS std  : {e_pot_abs.std():.4f} eV  ({e_pot_pa.std():.6f} eV/atom)')
print(f'  -> E_tot is {e_pot_abs.std()/e_tot_abs.std():.0f}x more stable than E_KS '
      '(same ratio either unit -- dividing both by the same atom count doesnt change it)')
print(f'Temperature: mean {t.mean():.1f} K (Nose target: 500 K)')
"
if command -v gnuplot >/dev/null 2>&1; then
    (cd "$OUT/thermo" && gnuplot siesta_thermo.gplot 2>/dev/null)
    echo "(rendered siesta_thermo.pdf -- the combined 4-panel energy/temperature/volume/pressure figure, Energy in eV/atom)"
fi
pause


echo "=================================================================="
echo " output/species-pair/  --  --pair restricts the RDF to one species pair"
echo "=================================================================="
mkdir -p "$OUT/species-pair"
cp aimd.ANI aimd.fdf aimd.out aimd.XV "$OUT/species-pair/"
echo "\$ stb-aimdAnalysis --label aimd --pair O-O --no-intro"
(cd "$OUT/species-pair" && stb-aimdAnalysis --label aimd --pair O-O --no-intro \
    | grep "RDF pair\|First peak")
pause


echo "=================================================================="
echo " output/full-report/  --  --save-report/--save-gnuplot/--view (all off by default)"
echo "=================================================================="
cat <<'EOF'
--save-report persists the full numbered report to disk;
--save-gnuplot writes a .dat + .gplot pair per computed quantity;
--view shows an interactive matplotlib preview (a no-op here, since
MPLBACKEND=Agg has no GUI to show it on -- try this yourself without
that env var set to actually see the plots pop up).
EOF
mkdir -p "$OUT/full-report"
cp aimd.ANI aimd.fdf aimd.out aimd.XV "$OUT/full-report/"
echo "\$ stb-aimdAnalysis --label aimd --track-atom 0 --track-pair 0-1 \\"
echo "      --save-report --save-gnuplot --view --no-intro"
(cd "$OUT/full-report" && stb-aimdAnalysis --label aimd --track-atom 0 --track-pair 0-1 \
    --save-report --save-gnuplot --view --no-intro > console.log 2>&1)
echo "Report sections written to stb_aimdAnalysis_report.txt:"
grep -o "^\[[0-9]\] [A-Za-z& ()-]*" "$OUT/full-report/stb_aimdAnalysis_report.txt" | sort -u
echo
echo "Files written:"
(cd "$OUT/full-report" && ls *.dat *.gplot *.txt)
pause


echo "=================================================================="
echo " output/trajectory/  --  --trajectory: generic ASE input (e.g. stb-mlmd's own output)"
echo "=================================================================="
cat <<'EOF'
Not every trajectory comes from SIESTA -- --trajectory reads any ASE
-readable multi-frame file (xsf/pdb/xyz) directly, e.g. one written by
stb-mlmd. An extended-xyz written with per-frame 'Time' info auto
-detects dt; otherwise --dt is required. Since this path is NOT
guaranteed to be a SIESTA run, the [8] REFERENCES section correctly
skips the SIESTA citation instead of assuming one:
EOF
mkdir -p "$OUT/trajectory"
python3 -c "
from ase import Atoms
from ase.io import write
import numpy as np
rng = np.random.default_rng(0)
cell = [[10.0, 0, 0], [0, 10.0, 0], [0, 0, 10.0]]
frames = []
for i in range(8):
    pos = np.array([[5.0, 5.0, 4.395], [5.0, 5.0, 5.605]]) + rng.normal(scale=0.05, size=(2, 3))
    a = Atoms(symbols=['O', 'O'], positions=pos, cell=cell, pbc=True)
    a.info['Time'] = i * 15.0
    frames.append(a)
write('$OUT/trajectory/synthetic_traj.xyz', frames, format='extxyz')
"
echo "\$ stb-aimdAnalysis --trajectory synthetic_traj.xyz --track-pair 0-1 --no-intro"
(cd "$OUT/trajectory" && stb-aimdAnalysis --trajectory synthetic_traj.xyz --track-pair 0-1 \
    --no-intro | grep "Auto-detected dt\|REFERENCES\|No SIESTA-specific")
pause


echo "=================================================================="
echo " Two ways to run it"
echo "=================================================================="
cat <<EOF
A -- direct CLI:
  stb-aimdAnalysis --label aimd --track-atom 0 --track-pair 0-1

B -- interactive stb-suite menu:
  \$ stb-suite
  Select an option (0-6, or a tool code like 4.1.2): 3.18

Both paths call the exact same underlying tool -- proven directly below.
The menu asks for the input source (label/generic trajectory), stride/
skip, the RDF species pair -- then asks whether to list every atom's
index/species/coordinates (y/N) right before the NEW --track-atom/
--track-pair prompts, so you can check valid indices before being asked
to pick one -- then the output directory, and save-report/save-gnuplot/
view prompts.
EOF
TMP="$(mktemp -d)"
cp aimd.ANI aimd.fdf aimd.out aimd.XV "$TMP/"
echo
echo "\$ stb-aimdAnalysis --label aimd --track-atom 0 --track-pair 0-1 --no-intro"
(cd "$TMP" && stb-aimdAnalysis --label aimd --track-atom 0 --track-pair 0-1 \
    --no-intro > cli.log 2>&1)
echo
echo "\$ printf '3.18\\n\\naimd\\n\\n\\n\\n\\nn\\n0\\n0-1\\n\\nn\\nn\\nn\\n\\n0\\n' | stb-suite"
(cd "$TMP" && printf '3.18\n\naimd\n\n\n\n\nn\n0\n0-1\n\nn\nn\nn\n\n0\n' | timeout 60 stb-suite > session.log 2>&1) || true
CLI_LINE=$(grep "Pair 0-1 distance" "$TMP/cli.log")
MENU_LINE=$(grep "Pair 0-1 distance" "$TMP/session.log")
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
Nine self-contained folders were generated under output/:
  basic-run/       geometry-file/   list-atoms/      tracking/
  thermo/          species-pair/    full-report/     trajectory/

Only full-report/ has the persisted .txt report + .dat/.gplot files;
the rest only ever write references.bib (SIESTA-input runs) or nothing
at all (the generic --trajectory run and --list-atoms).

Recap of what this walkthrough covered:
  - what stb-aimdAnalysis computes: RDF, MSD/diffusion, VACF-derived VDOS
  - two NEW features: --track-atom (one atom's own displacement) and
    --track-pair (minimum-image distance between two specific atoms),
    cross-checked live against the RDF's independently-computed peak
  - NEW --list-atoms: index/species/coordinates for every atom, fast
    (first frame only), then exits -- off by default (hundreds of atoms
    would make it unwieldy in every report); the interactive stb-suite
    menu now asks (y/N) before showing it, right before the
    --track-atom/--track-pair prompts
  - a real bug fix: --geometry-file, for when the real .fdf isn't named
    <label>.fdf (almost always the case in practice)
  - NEW [7] THERMODYNAMIC TIME SERIES: energy/temperature/volume/pressure
    (+ energy per atom) in one 4-panel figure, read from SIESTA's own
    <label>.MDE file -- verified live on a real 500-step Nose (NVT) run
    that E_tot is ~140x more stable than E_KS, the same E_pot-vs-E_total
    lesson stb-mlmd learned live this session
  - the numbered [0]...[10] report, --save-report, --save-gnuplot (now
    opt-in, previously unconditional PNGs with no gnuplot output at all),
    --view (now opt-in, previously always generated with no way to skip)
  - --trajectory: generic ASE input independent of SIESTA, and how the
    REFERENCES section correctly adapts to that (no assumed SIESTA cite)
  - CLI and the interactive stb-suite menu building the same command

As a next step, try on your own with a real SIESTA AIMD run:
  stb-aimdAnalysis --label my_calc --geometry-file my_calc.fdf \\
      --track-atom 12 --track-pair 12-47 \\
      --save-report --save-gnuplot --view
EOF
