#!/bin/bash
# Guided example: NEB / Reaction Path workflow
# (stb-neb / stb-nebCycle / stb-nebAnalysis, codes 4.9.1 / 4.9-CLI-only / 4.9.2)
#
# Not an automated test (see test/4-workflow/9-neb/{prep,cycle,analysis}/
# test.sh for that) -- a commented walk-through: it runs real commands, one
# group at a time, and shows you the piece of output that proves what just
# happened. Pauses between sections so you can read before moving on.
#
# stb-neb and stb-nebCycle are exercised for real (structure interpolation,
# a real MACE-MP-0 climbing-image relax for --mode 1, and a real ASE FIRE
# NEB step for --mode 2/--mode 3's refinement loop -- none of that needs a
# SIESTA binary). stb-nebAnalysis needs real SIESTA .out files to analyze,
# which this walkthrough doesn't have -- so its worked example (Section
# "output/mode3/" below) fabricates calc.out with a hand-chosen, exactly
# SYMMETRIC pair of energies (the H adatom hopping between two chemically
# -equivalent carbons), designed to make the forward/backward barrier and
# the spline-fitted transition state come out at clean, easy-to-check
# numbers. See the README's Section 8 for the full arithmetic.
#
# Needs the optional 'ml' extra (pip install stb_suite[ml]) for the
# --mode 1 section -- checked up front, skipped with a clear note if
# unavailable (the rest of the walkthrough does not need it).

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

# Fabricates a 'relaxed' siesta.XV for a directory holding a structure.fdf,
# at the SAME geometry (this walkthrough has no real SIESTA to relax it
# further with) -- purely to demonstrate the FILE mechanism stb-neb's
# read_relaxed_or_input uses (prefer a directory's own finished siesta.XV
# over its pre-relaxation structure.fdf guess), same recipe
# examples/4.8-adsorption/example_4.8.sh's own fabricate_relaxed_xv uses.
fabricate_relaxed_xv() {
    local dir="$1"
    python3 -c "
import sisl
from stb.core import structure_io
fdf = structure_io.read_fdf('$dir/structure.fdf')
pmg = structure_io.to_pymatgen(fdf)
atoms = [sisl.Atom(str(s.specie)) for s in pmg]
geom = sisl.Geometry(pmg.cart_coords, atoms=atoms, lattice=sisl.Lattice(pmg.lattice.matrix))
sisl.get_sile('$dir/siesta.XV', mode='w').write_geometry(geom)
"
}

# Writes a synthetic calc.out with a real 'siesta: FreeEng' line plus a
# clean SCF-convergence line and a residual-force line, so the quality
# diagnostics in stb-nebAnalysis read clean (or, deliberately, don't --
# image_02 below is given a large residual force on purpose).
write_freeeng() {
    local path="$1" energy="$2" maxforce="${3:-0.010000}"
    printf 'siesta: FreeEng =    %s\nSCF cycle converged after 12 iterations\nsiesta: Atomic forces (eV/Ang):\n   Max    %s\n' \
        "$energy" "$maxforce" > "$path"
}

# One real-DFT NEB step's worth of fabricated .FA/calc.out: a smooth,
# physically-motivated harmonic "pull toward a fixed per-image target"
# force field (F = -K*(pos-target), E = 0.5*K*|pos-target|^2) -- NOT
# random noise, which test/4-workflow/9-neb/cycle/test.sh found live
# triggers a real ASE 'improvedtangent' 0/0 tangent edge case. Targets are
# derived once from cycle_00's own images and cached in targets.npy, so
# every call pulls toward the SAME fixed point regardless of which cycle
# is being fabricated -- exactly the recipe that test.sh itself uses.
write_cycle_fabricator() {
cat > "$1/fabricate.py" << 'PYEOF'
import sys, os
import numpy as np
from stb.core import structure_io

ROOT, CYCLE = sys.argv[1], sys.argv[2]
K = 2.0

def positions(path):
    return structure_io.to_pymatgen(structure_io.read_fdf(path)).cart_coords

def write_fa(path, forces):
    n = len(forces)
    with open(path, 'w') as f:
        f.write(f'{n}\n')
        for i, (fx, fy, fz) in enumerate(forces, start=1):
            f.write(f'{i} {fx:.6f} {fy:.6f} {fz:.6f}\n')

def write_out(path, energy):
    with open(path, 'w') as f:
        f.write(f'siesta: FreeEng =    {energy:.6f}\nSCF cycle converged after 10 iterations\n')

labels = sorted(d for d in os.listdir(f'{ROOT}/cycle_00') if d.startswith('image_'))
targets_path = f'{ROOT}/targets.npy'
if not os.path.exists(targets_path):
    rng = np.random.default_rng(7)
    targets = {}
    for label in labels:
        p0 = positions(f'{ROOT}/cycle_00/{label}/structure.fdf')
        targets[label] = p0 + (rng.random(p0.shape) - 0.5) * 0.15
    np.save(targets_path, targets, allow_pickle=True)
targets = np.load(targets_path, allow_pickle=True).item()

for label in labels:
    d = f'{ROOT}/cycle_{CYCLE}/{label}'
    pos = positions(f'{d}/structure.fdf')
    disp = pos - targets[label]
    write_fa(os.path.join(d, 'siesta.FA'), -K * disp)
    write_out(os.path.join(d, 'calc.out'), 0.5 * K * float(np.sum(disp**2)) - 300.0)
PYEOF
}

HAS_ML=0
if python3 -c "import mace" 2>/dev/null; then HAS_ML=1; fi

echo "=================================================================="
echo " Why this workflow needs two stages, and what NEB actually computes"
echo "=================================================================="
cat <<'EOF'
A reaction/diffusion barrier is the energy cost of the HIGHEST point along
the path connecting two already-relaxed structures -- e.g. an adsorbate
hopping from one adsorption site to a neighboring one. The Nudged Elastic
Band method finds that path (and its saddle point) by placing a chain of
"images" between the two endpoints, connecting neighboring images with
virtual springs (so they can't bunch up or spread out unevenly), and
relaxing every image sideways (perpendicular to the path) while the spring
force keeps them evenly spaced along it. The "climbing image" refinement
then lets the single highest-energy image climb UPHILL along the path
direction too, converging it exactly onto the true saddle point instead of
just the highest of a handful of quantized sample points.

Stage 1 (stb-neb) turns two already-relaxed endpoint structures into a set
of image_NN/ SIESTA folders along that path -- interpolated, and in two of
its three modes, actually shaped into a real minimum-energy path first
(Section 4/5 below). You run SIESTA in every image_NN/ (always a single
-point evaluation -- Section 3 explains why). Stage 2 (stb-nebAnalysis)
reads every image's energy back and reports the barrier.
EOF
pause

echo "=================================================================="
echo " Section 2 live: endpoints must already be relaxed"
echo "=================================================================="
cat <<'EOF'
initial.fdf/final.fdf in this folder are given directly as bare .fdf files
-- stb-neb uses those exactly as given, no relaxation-status concept
applies (README Section 2.1). But --initial/--final also accept a
DIRECTORY, in which case stb-neb looks for that directory's own finished
siesta.XV, preferring it over the raw structure.fdf guess -- and WARNS if
it isn't there yet. This case builds that exact situation: two folders
holding only a structure.fdf, no siesta.XV.
EOF
mkdir -p "$OUT/endpoints_unrelaxed/init_dir" "$OUT/endpoints_unrelaxed/final_dir"
cp initial.fdf "$OUT/endpoints_unrelaxed/init_dir/structure.fdf"
cp final.fdf "$OUT/endpoints_unrelaxed/final_dir/structure.fdf"
echo "\$ stb-neb -i init_dir -f final_dir -c calc.fdf --mode 3"
stb-neb -i "$OUT/endpoints_unrelaxed/init_dir" -f "$OUT/endpoints_unrelaxed/final_dir" \
    -c calc.fdf --mode 3 -O "$OUT/endpoints_unrelaxed" --no-intro \
    > "$OUT/endpoints_unrelaxed.log"
grep "WARNING.*no finished 'siesta.XV'" "$OUT/endpoints_unrelaxed.log"
echo
echo "The band was still generated (useful to preview) -- but re-run once"
echo "SIESTA has actually relaxed each endpoint. Fabricating that now"
echo "(writing siesta.XV at the SAME geometry -- no real SIESTA available"
echo "here -- purely to demonstrate the file stb-neb looks for):"
fabricate_relaxed_xv "$OUT/endpoints_unrelaxed/init_dir"
fabricate_relaxed_xv "$OUT/endpoints_unrelaxed/final_dir"
stb-neb -i "$OUT/endpoints_unrelaxed/init_dir" -f "$OUT/endpoints_unrelaxed/final_dir" \
    -c calc.fdf --mode 3 -O "$OUT/endpoints_relaxed" --no-intro \
    > "$OUT/endpoints_relaxed.log"
if grep -q "WARNING.*no finished 'siesta.XV'" "$OUT/endpoints_relaxed.log"; then
    echo "Unexpected: the not-yet-relaxed WARNING is STILL present."
    exit 1
else
    echo "Confirmed: with siesta.XV present, the WARNING is GONE."
fi
pause

echo "=================================================================="
echo " Section 2 live: same composition, same lattice -- hard requirements"
echo "=================================================================="
echo "--- Composition mismatch (final.fdf missing the H atom) ---"
sed '/0.111111111   0.222222222   0.600000000   2/d; s/NumberofAtoms      19/NumberofAtoms      18/' \
    final.fdf > "$OUT/final_bad_composition.fdf"
mkdir -p "$OUT/composition_mismatch"
if stb-neb -i initial.fdf -f "$OUT/final_bad_composition.fdf" -c calc.fdf --mode 3 \
    -O "$OUT/composition_mismatch" --no-intro > "$OUT/composition_mismatch.log" 2>&1; then
    echo "Unexpected: stb-neb should have refused this."
    exit 1
fi
grep "ERROR.*different composition" "$OUT/composition_mismatch.log"
echo
echo "--- Lattice mismatch (final.fdf's cell scaled by ~1.6%) ---"
python3 -c "
with open('final.fdf') as f:
    text = f.read()
text = (text.replace('7.3800000000', '7.5000000000')
            .replace('-3.6900000000', '-3.7500000000')
            .replace('6.3912673200', '6.4948718700'))
with open('$OUT/final_bad_lattice.fdf', 'w') as f:
    f.write(text)
"
mkdir -p "$OUT/lattice_mismatch"
if stb-neb -i initial.fdf -f "$OUT/final_bad_lattice.fdf" -c calc.fdf --mode 3 \
    -O "$OUT/lattice_mismatch" --no-intro > "$OUT/lattice_mismatch.log" 2>&1; then
    echo "Unexpected: stb-neb should have refused this."
    exit 1
fi
grep "ERROR.*different lattices" "$OUT/lattice_mismatch.log"
echo
echo "Both are hard [ERROR]s, exit 1, nothing written -- stb-neb never"
echo "guesses or silently overrides one endpoint's cell with the other's"
echo "(README Section 2.3 explains exactly why: no downstream library"
echo "supports a per-image cell change anyway)."
pause

echo "=================================================================="
echo " output/mode3/  --  Stage 1 mode 3 (100% real-DFT, no MACE) + worked barrier"
echo "=================================================================="
echo "\$ stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 5 --mode 3"
stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 5 --mode 3 -O "$OUT/mode3" --no-intro \
    > "$OUT/mode3_prep.log"
echo
echo "5 images, real-DFT single-point folders directly (no MACE stage):"
ls "$OUT/mode3/neb_run/cycle_00"
echo
echo "Every image_NN/config_extra.fdf forces a single-point evaluation,"
echo "unconditionally, regardless of what calc.fdf's own MD.TypeOfRun says:"
cat "$OUT/mode3/neb_run/cycle_00/image_02/config_extra.fdf"
pause

cat <<'EOF'
Worked example: fabricated calc.out for all 5 images, hand-chosen so
image_00 and image_04 (physically equivalent carbons, by symmetry) share
the EXACT SAME energy -- a genuine symmetric hop, reaction energy = 0 by
construction. image_02 (the transition state) is also given a large
residual force on purpose, to show stb-nebAnalysis's quality warning.
EOF
write_freeeng "$OUT/mode3/neb_run/cycle_00/image_00/calc.out" "-300.000000" "0.010000"
write_freeeng "$OUT/mode3/neb_run/cycle_00/image_01/calc.out" "-299.700000" "0.015000"
write_freeeng "$OUT/mode3/neb_run/cycle_00/image_02/calc.out" "-299.500000" "0.400000"
write_freeeng "$OUT/mode3/neb_run/cycle_00/image_03/calc.out" "-299.700000" "0.015000"
write_freeeng "$OUT/mode3/neb_run/cycle_00/image_04/calc.out" "-300.000000" "0.010000"
echo "\$ stb-nebAnalysis --dir neb_run --save-report --save-gnuplot --save-path-xyz"
(cd "$OUT/mode3" && stb-nebAnalysis --dir neb_run --save-report --save-gnuplot \
    --save-path-xyz --no-intro > "$OUT/mode3_analysis.log")
sed -n '/\[1\] IMAGE ENERGIES/,/\[2\] CONSISTENCY/p' "$OUT/mode3_analysis.log" | head -n -1
sed -n '/\[3\] BARRIER ANALYSIS/,/\[3b\]/p' "$OUT/mode3_analysis.log" | head -n -1
echo
echo "Forward barrier = backward barrier = 0.5 eV, reaction energy = 0 eV,"
echo "TS sits at exactly 50% of the path -- exactly what a symmetric hop"
echo "between two equivalent carbons should give. The residual-force"
echo "WARNING on image_02 (0.40 eV/Ang, above --force-tolerance) also fired"
echo "as expected -- advisory only, doesn't block the barrier estimate."
pause

echo "=================================================================="
echo " output/mode1/  --  Stage 1 mode 1: MACE-MP-0 shapes the path first"
echo "=================================================================="
if [ "$HAS_ML" -eq 1 ]; then
cat <<'EOF'
Mode 1 relaxes a REAL climbing-image NEB on the MACE-MP-0 foundation model
before writing the single-point SIESTA folders -- a much better geometry
than plain linear interpolation, at zero SIESTA cost. The spring constant
k (eV/Ang^2) keeps neighboring images evenly spaced; two stages run in
sequence: stage 1 shapes the whole band (climb=False), stage 2 lets the
single highest-energy image climb uphill onto the true saddle point
(climb=False -> True, the same reasoning core/mace_relax.py documents).
--ml-max-steps is capped low here to keep this walkthrough fast (a real
run should raise it until [4] reports converged, not "hit step cap").
EOF
echo "\$ stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 5 --mode 1 --ml-max-steps 30"
stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 5 --mode 1 --ml-max-steps 30 \
    -O "$OUT/mode1" --no-intro > "$OUT/mode1.log"
sed -n '/\[4\] MACE-MP-0 PATH SHAPING/,/\[5\] PATH QUALITY/p' "$OUT/mode1.log" | head -n -1
echo
echo "image_NN/ folders were written directly under neb_run/ (no cycle_00/"
echo "nesting -- mode 1 is a single-shot, no real-DFT refinement loop):"
ls "$OUT/mode1/neb_run"
ls "$OUT/mode1/neb_run/image_02"
else
echo "Optional 'ml' extra not installed (pip install stb_suite[ml]) -- skipping"
echo "this section. --mode 2/--mode 3 below don't need it."
fi
pause

echo "=================================================================="
echo " output/mode2_cycle/  --  Stage 1 mode 2 + stb-nebCycle: real-DFT refinement"
echo "=================================================================="
cat <<'EOF'
Mode 2 writes cycle_00/image_NN/ (same as mode 3) but is meant to be
followed by real-DFT refinement cycles: run SIESTA in cycle_00/, then call
stb-nebCycle (a separate CLI-only tool, not in the interactive menu -- see
"CLI-only tools" in CLAUDE.md) to take ONE real-DFT NEB step using ASE's
FIRE optimizer (state persisted to neb_cycle_state.json between calls, so
each call is a fresh, short-lived process -- appropriate for a cluster job
array, not a long-running one). It either writes cycle_01/ (the next
geometry to run SIESTA on) or, once the max residual force drops below
--fmax, writes a NEB_CONVERGED sentinel and stops.
EOF
echo "\$ stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 5 --mode 2"
stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 5 --mode 2 -O "$OUT/mode2_cycle" --no-intro \
    > "$OUT/mode2_prep.log"
echo
echo "The printed cluster-submission snippet (run this yourself in a real"
echo "job-array loop -- --max-cycles/--siesta-exe/--mpirun-np/--conda-env"
echo "all customize it):"
sed -n '/\[7\] CLUSTER SUBMISSION/,/\[8\] SUMMARY/p' "$OUT/mode2_prep.log" | head -n -1
pause

echo "--- Fabricating one real-DFT NEB step (harmonic 'pull toward a fixed target' forces) ---"
write_cycle_fabricator "$OUT/mode2_cycle/neb_run"
(cd "$OUT/mode2_cycle/neb_run" && python3 fabricate.py . 00)
echo "\$ stb-nebCycle --dir neb_run --fmax 0.05 --climb-after 0"
(cd "$OUT/mode2_cycle" && stb-nebCycle --dir neb_run --fmax 0.05 --climb-after 0 --no-intro \
    > "$OUT/mode2_cycle1.log" 2>&1)
grep "Max residual force" "$OUT/mode2_cycle1.log"
echo
if [ -d "$OUT/mode2_cycle/neb_run/cycle_01" ]; then
    echo "Confirmed: cycle_01/ was written (the next geometry) -- residual force"
    echo "was still above --fmax."
else
    echo "Unexpected: cycle_01/ was not written."
    exit 1
fi
pause

echo "--- A second step, with a looser --fmax so it actually converges ---"
(cd "$OUT/mode2_cycle/neb_run" && python3 fabricate.py . 01)
echo "\$ stb-nebCycle --dir neb_run --fmax 0.5 --climb-after 0"
(cd "$OUT/mode2_cycle" && stb-nebCycle --dir neb_run --fmax 0.5 --climb-after 0 --no-intro \
    > "$OUT/mode2_cycle2.log" 2>&1)
grep "CONVERGED" "$OUT/mode2_cycle2.log"
if [ -f "$OUT/mode2_cycle/neb_run/NEB_CONVERGED" ]; then
    echo "Confirmed: NEB_CONVERGED written -- a real submission-script loop would"
    echo "check for this file and break out (see the printed snippet above)."
else
    echo "Unexpected: NEB_CONVERGED was not written."
    exit 1
fi
pause

echo "=================================================================="
echo " Proof: CLI and the interactive stb-suite menu agree (Stage 1, mode 3)"
echo "=================================================================="
echo "Driving 4.9.1 through the interactive menu (piped input) and comparing"
echo "against the same direct-CLI case from output/mode3/ above."
TMP="$(mktemp -d)"
cp initial.fdf final.fdf calc.fdf "$TMP/"
echo
echo "\$ printf '4.9.1\\ninitial.fdf\\nfinal.fdf\\ncalc.fdf\\n\\n\\n\\n\\n5\\nn\\n3\\n\\n\\n\\n\\n\\n\\n\\n\\n0\\n' | stb-suite"
(cd "$TMP" && printf '4.9.1\ninitial.fdf\nfinal.fdf\ncalc.fdf\n\n\n\n\n5\nn\n3\n\n\n\n\n\n\n\n\n0\n' \
    | stb-suite > menu1.log 2>&1)
MENU_IMG="$TMP/neb_run/cycle_00/image_02/structure.fdf"
CLI_IMG="$OUT/mode3/neb_run/cycle_00/image_02/structure.fdf"
if python3 -c "
import sys
from stb.core import structure_io
import numpy as np
a = structure_io.to_pymatgen(structure_io.read_fdf('$MENU_IMG'))
b = structure_io.to_pymatgen(structure_io.read_fdf('$CLI_IMG'))
sys.exit(0 if np.allclose(a.cart_coords, b.cart_coords, atol=1e-6) else 1)
"; then
    echo "Confirmed: the interactive-menu image_02 and the direct-CLI image_02 have"
    echo "the same geometry."
else
    echo "Unexpected: interactive-menu and direct-CLI results differ."
    exit 1
fi
rm -rf "$TMP"
pause

echo "=================================================================="
echo " Done"
echo "=================================================================="
cat <<'EOF'
Folders generated under output/:
  endpoints_unrelaxed/   endpoints_relaxed/   composition_mismatch/
  lattice_mismatch/      mode3/               mode1/ (needs the 'ml' extra)
  mode2_cycle/

output/mode3/neb_run/ has the full Stage 1+2 report pair (neb_setup.txt,
neb_report.txt), the energy-profile plot/gnuplot data, and
neb_path_current.xyz (viewable in VESTA/OVITO/ASE-GUI).

As a next step, on your OWN already-relaxed endpoints:
  stb-neb -i <initial> -f <final> -c <calc.fdf> --mode 2
  # run SIESTA in every neb_run/cycle_00/image_NN/ folder, then loop:
  stb-nebCycle --dir neb_run --fmax 0.05 --climb-after 0
  # run SIESTA in the newly-written cycle_NN/, repeat until NEB_CONVERGED
  stb-nebAnalysis --dir neb_run --save-report --save-gnuplot --save-path-xyz
EOF
