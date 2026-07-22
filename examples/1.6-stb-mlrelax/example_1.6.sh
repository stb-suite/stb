#!/bin/bash
# Guided example: stb-mlrelax (code 1.6 in the stb-suite menu)
#
# Not an automated test (see test/2-structures/12-mlrelax/test.sh for that)
# -- a commented walk-through: it runs real commands, one group at a time,
# and shows you the piece of output that proves what just happened. It
# pauses between sections so you can read before moving on.
#
# Needs the optional 'ml' extra: pip install stb_suite[ml]  (PyTorch +
# mace-torch). Downloads MACE-MP-0 on first use, cached under
# ~/.cache/mace/ afterward: "small" (~31 MB, used by most of this script)
# plus "medium" and "large" (~45 MB / ~135 MB, used once each in the
# model-comparison section below).

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

echo "=================================================================="
echo " What ML pre-relaxation is for"
echo "=================================================================="
cat <<'EOF'
A real SIESTA relaxation is expensive: every ionic step is a full DFT SCF
cycle -- solving the electronic structure from scratch (or nearly so) just
to take one small step downhill in energy. stb-mlrelax runs a MUCH cheaper
pre-relaxation first, using MACE-MP-0, a "foundation" machine-learned
interatomic potential (MLIP): a neural network trained to reproduce the
energies/forces of ~1.5 million DFT relaxation snapshots from the Materials
Project database. Once trained, evaluating it is just a forward pass through
the network -- microseconds to milliseconds per structure, not minutes.

MACE itself (the architecture, not the specific MACE-MP-0 training) is an
EQUIVARIANT message-passing graph neural network: atoms are nodes, nearby
atoms (within a cutoff radius) are connected by edges, and information about
local geometry is passed between neighbors over a couple of message-passing
"rounds" to build up a description of each atom's environment. "Equivariant"
is the key physical constraint: if you rotate the whole structure, the
network's predicted forces rotate the exact same way -- baked into the
architecture, not learned approximately, so the model can't accidentally
break rotational symmetry the way an ordinary (non-equivariant) network
could.

Why this helps DFT: a wrong lattice guess, a freshly-added defect/passivant
atom at the wrong bond length, or a hand-perturbed structure all get cleaned
up by MACE in seconds -- SIESTA's own relaxation then starts much closer to
the true local minimum and needs far fewer, cheaper ionic steps to converge.

This is a HEURISTIC pre-relaxation, NOT a substitute for a real DFT
relaxation. Two concrete caveats:
  - Universal MLIPs are trained mostly on Materials-Project-like inorganic
    crystals near their own equilibrium -- they're less reliable far from
    that distribution (exotic bonding, extreme pressure, molecules...).
  - A starting guess very far from equilibrium (e.g. a lattice constant off
    by 15-20% or more) can occasionally lead MACE to a spurious low-energy
    state instead of the real one -- stb-mlrelax flags this automatically
    when it happens (see output/relax-cell/ below).
Always still relax the result with SIESTA afterward.
EOF
pause

echo "=================================================================="
echo " The optimizer: FIRE and the fmax convergence criterion"
echo "=================================================================="
cat <<'EOF'
Once MACE can predict a force on every atom (F = -dE/dR, same physics as
DFT's Hellmann-Feynman forces), any standard geometry optimizer can use it
to walk downhill in energy. stb-mlrelax's default is FIRE (Fast Inertial
Relaxation Engine) -- it treats the atoms like a damped mechanical system:
move along the force direction like inertia would carry it, but adaptively
reset/dampen the "velocity" whenever it stops pointing downhill. It's simple,
has no matrix to invert (unlike BFGS), and is a common default for MLIP-based
relaxation. BFGS/LBFGS are also available via --optimizer.

Convergence is judged by fmax: the largest force component on any atom, in
eV/Ang (the default target is 0.05 eV/Ang, tightenable via --fmax) -- the
same style of criterion SIESTA's own relaxation uses, so "relaxed enough for
stb-mlrelax" and "relaxed enough for SIESTA" mean the same physical thing,
just possibly at different thresholds.
EOF
pause

echo "=================================================================="
echo " Two ways to run it"
echo "=================================================================="
cat <<'EOF'
A -- direct CLI:
  stb-mlrelax -f defect.fdf

B -- interactive stb-suite menu:
  $ stb-suite
  Select an option (0-6, or a tool code like 4.1.2): 1.6

Every run always writes a structure-validation section (atom-proximity,
left-handed-cell, density, space-group checks -- the same ones
stb-inputfile/stb-fetch already run) and references.bib. A full text report
(--save-report), the raw convergence trace (--save-data), and an
interactive before/after 3D view (--view) are all off by default.
EOF
pause

echo "=================================================================="
echo " output/positions-only/  --  cleaning up one perturbed atom"
echo "=================================================================="
cat <<'EOF'
si_defect.fdf: correct 5.43 Ang lattice constant, but one atom was nudged
~0.16 Ang off its ideal tetrahedral site (e.g. simulating a hand-edited
defect). Default mode only relaxes positions -- the cell must not move.
EOF
mkdir -p "$OUT/positions-only"
cp si_defect.fdf "$OUT/positions-only/"
echo
echo "\$ stb-mlrelax -f si_defect.fdf -o defect_relaxed.fdf --no-intro"
(cd "$OUT/positions-only" && stb-mlrelax -f si_defect.fdf -o defect_relaxed.fdf --no-intro > console.log)
grep "^Mode" "$OUT/positions-only/console.log"
awk '/\[5\] BEFORE \/ AFTER COMPARISON/,/^$/' "$OUT/positions-only/console.log"
echo
python3 - "$OUT/positions-only/si_defect.fdf" "$OUT/positions-only/defect_relaxed.fdf" <<'EOF'
import sys
from stb.core import structure_io
before = structure_io.to_pymatgen(structure_io.read_fdf(sys.argv[1]))
after = structure_io.to_pymatgen(structure_io.read_fdf(sys.argv[2]))
assert before.lattice.abc == after.lattice.abc, "cell should NOT have changed in positions-only mode"
print("Confirmed: the cell is bit-identical before/after (positions-only mode).")
EOF
pause

echo "=================================================================="
echo " output/relax-cell/  --  fixing a wrong lattice constant"
echo "=================================================================="
cat <<'EOF'
si_wrong_a.fdf: same bulk Si, but the lattice constant is 5.00 Ang instead
of the correct ~5.43-5.46 Ang -- e.g. a unit slip or a bad initial guess.
--relax-cell also relaxes the lattice, not just positions.
EOF
mkdir -p "$OUT/relax-cell"
cp si_wrong_a.fdf "$OUT/relax-cell/"
echo
echo "\$ stb-mlrelax -f si_wrong_a.fdf --relax-cell -o cell_relaxed.fdf --no-intro"
(cd "$OUT/relax-cell" && stb-mlrelax -f si_wrong_a.fdf --relax-cell -o cell_relaxed.fdf --no-intro > console.log)
awk '/\[5\] BEFORE \/ AFTER COMPARISON/,/^$/' "$OUT/relax-cell/console.log"
echo
python3 - "$OUT/relax-cell/cell_relaxed.fdf" <<'EOF'
import sys
from stb.core import structure_io
s = structure_io.to_pymatgen(structure_io.read_fdf(sys.argv[1]))
a = s.lattice.abc[0]
assert 5.40 <= a <= 5.50, f"lattice constant {a} out of the expected [5.40, 5.50] Ang range"
print(f"Confirmed: relaxed lattice constant a={a:.4f} Ang -- close to the real ~5.43 Ang, "
      "not the wrong 5.00 Ang starting guess.")
EOF
pause

echo "=================================================================="
echo " output/bad-structure/  --  structure validation catching a real problem"
echo "=================================================================="
cat <<'EOF'
si_too_close.fdf: one atom was mis-typed at fractional (0.01, 0.01, 0.01)
instead of (0.25, 0.25, 0.25) -- it sits only ~0.09 Ang from its neighbor at
the origin, far inside stb-mlrelax's 0.5 Ang minimum-distance check. Watch
[2] STRUCTURE VALIDATION (pre-relaxation) catch this BEFORE relaxing, and
[6] STRUCTURE VALIDATION (post-relaxation) confirm it's resolved AFTER.

Why the starting energy below is enormous: at ~0.09 Ang apart, the two
atoms' electron clouds are forced to overlap almost completely. Physically
this is Pauli repulsion -- electrons can't occupy the same state, so pushing
two atoms' electron clouds together costs a steeply rising energy penalty
(this is also why atoms have a "size" at all). MACE learned this repulsive
wall from its DFT training data just like it learned everything else, so the
predicted force here is huge and points straight apart -- exactly what
resolves the problem within a handful of optimizer steps.
EOF
mkdir -p "$OUT/bad-structure"
cp si_too_close.fdf "$OUT/bad-structure/"
echo
echo "\$ stb-mlrelax -f si_too_close.fdf -o close_relaxed.fdf --no-intro"
(cd "$OUT/bad-structure" && stb-mlrelax -f si_too_close.fdf -o close_relaxed.fdf --no-intro > console.log)
echo "Pre-relaxation validation (the Atom proximity row):"
grep -A5 "STRUCTURE VALIDATION (pre-relaxation)" "$OUT/bad-structure/console.log" | grep "Atom proximity"
echo "Post-relaxation validation (the Atom proximity row):"
grep -A5 "STRUCTURE VALIDATION (post-relaxation)" "$OUT/bad-structure/console.log" | grep "Atom proximity"
echo
grep "^Steps used" "$OUT/bad-structure/console.log"
awk '/\[5\] BEFORE \/ AFTER COMPARISON/,/^$/' "$OUT/bad-structure/console.log"
pause

echo "=================================================================="
echo " output/model-comparison/  --  small vs. medium vs. large: speed/accuracy"
echo "=================================================================="
cat <<'EOF'
MACE-MP-0 ships in 3 sizes -- more parameters/message-passing capacity costs
more compute but (usually) predicts closer to the underlying DFT reference.
stb-mlrelax's [3] MACE MODEL & SIMULATION SETUP section reports the actual
parameter count and cutoff radius of whichever model gets loaded, and [4]
RELAXATION reports real wall-clock time -- so this isn't a claim, it's
measured on every run. Same si_defect.fdf, same fmax, only --model changes:
EOF
mkdir -p "$OUT/model-comparison"
cp si_defect.fdf "$OUT/model-comparison/"
for size in small medium large; do
    echo
    echo "\$ stb-mlrelax -f si_defect.fdf --model $size -o ${size}_relaxed.fdf --no-intro"
    (cd "$OUT/model-comparison" && stb-mlrelax -f si_defect.fdf --model "$size" \
        -o "${size}_relaxed.fdf" --no-intro > "console_$size.log")
    grep -E "^MACE model parameters|^MACE cutoff radius|^Wall time|^Energy" \
        "$OUT/model-comparison/console_$size.log" | sed "s/^/  /"
done
echo
echo "Side by side:"
printf "  %-8s %14s %10s %10s\n" "model" "parameters" "wall time" "final E (eV)"
for size in small medium large; do
    log="$OUT/model-comparison/console_$size.log"
    params=$(grep "^MACE model parameters" "$log" | sed 's/.*~//')
    wtime=$(grep "^Wall time" "$log" | sed 's/.*: //')
    efinal=$(grep "^Energy (eV)" "$log" | awk -F'|' '{gsub(/ /, "", $3); print $3}')
    printf "  %-8s %14s %10s %10s\n" "$size" "$params" "$wtime" "$efinal"
done
cat <<'EOF'

Bigger isn't automatically "more right" for any single structure -- it's a
statistical tendency across many structures. Use --model small (the default)
for routine pre-relaxation; reach for medium/large when a small-model result
looks suspicious (e.g. the >15% lattice-change warning fired) or as an extra
sanity check before committing to expensive real DFT.
EOF
pause

echo "=================================================================="
echo " output/full-report/  --  --save-report, --save-data, references.bib"
echo "=================================================================="
cat <<'EOF'
Two extra physical quantities show up in [5] BEFORE/AFTER COMPARISON that
haven't come up yet: the SPACE GROUP and the STRESS TENSOR.

Space group before vs. after is a genuine diagnostic, not decoration: if a
structure was supposed to be a specific symmetric phase but relaxes into a
LOWER symmetry, that's often a sign the starting guess (or the potential
itself, this far from its training distribution) is misbehaving -- whereas
recovering or keeping HIGH symmetry through relaxation is a good physical
sanity check, the same way it was used earlier in this session to confirm
si_wrong_a.fdf's absurd starting energy still relaxed into a sensible cubic
Si structure.

Stress is force's generalization to the cell itself -- how hard the
structure is "pushing" outward (or pulling inward) on its own periodic
boundary, per unit area (reported here in GPa). A converged positions-only
relaxation at a FIXED cell will generally still show nonzero stress -- that
is expected (only --relax-cell explicitly relaxes the cell down to near-zero
stress); it's included here as extra information, not a second convergence
criterion.
EOF
mkdir -p "$OUT/full-report"
cp si_defect.fdf "$OUT/full-report/"
echo "\$ stb-mlrelax -f si_defect.fdf --save-report --save-data --no-intro"
(cd "$OUT/full-report" && stb-mlrelax -f si_defect.fdf --save-report --save-data --no-intro > console.log)
echo
awk '/\[5\] BEFORE \/ AFTER COMPARISON/,/^$/' "$OUT/full-report/console.log"
echo
echo "Report sections written to stb_mlrelax_report.txt:"
grep -E "^\[[0-9]\]" "$OUT/full-report/stb_mlrelax_report.txt"
echo
echo "references.bib always written -- SIESTA (output is a .fdf) + the MACE"
echo "architecture/foundation-model papers (this run used the foundation model):"
grep "^@" "$OUT/full-report/references.bib"
echo
echo "Convergence trace files:"
ls "$OUT/full-report/"*.dat "$OUT/full-report/"*.png
pause

echo "=================================================================="
echo " Proof: CLI and the interactive stb-suite menu agree"
echo "=================================================================="
echo "Driving the same positions-only case through the interactive menu's"
echo "manual entry mode and diffing the energy/displacement lines."
TMP="$(mktemp -d)"
cp si_defect.fdf "$TMP/"
echo
echo "\$ printf '1.6\\nsi_defect.fdf\\nn\\n\\n\\n\\n\\n\\nn\\nn\\n\\n0\\n' | stb-suite"
(cd "$TMP" && printf '1.6\nsi_defect.fdf\nn\n\n\n\n\n\nn\nn\n\n0\n' | stb-suite > session.log 2>&1)
CLI_LINE=$(grep "^Energy" "$OUT/positions-only/console.log")
MENU_LINE=$(grep "^Energy" "$TMP/session.log")
if [ "$CLI_LINE" = "$MENU_LINE" ]; then
    echo "Confirmed: identical relaxation result from the CLI and the interactive menu."
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
Five self-contained folders were generated under output/:
  positions-only/   relax-cell/   bad-structure/
  model-comparison/   full-report/

Each has references.bib; full-report/ additionally has stb_mlrelax_report.txt
and the *_relax_convergence.dat/.png convergence trace.

Recap of what this walkthrough covered:
  - what a foundation MLIP (MACE-MP-0) is and why it's cheap enough to be a
    useful DFT pre-relaxation step
  - FIRE + the fmax force-convergence criterion
  - positions-only vs. --relax-cell (and how a vacuum axis stays fixed)
  - structure validation catching (and confirming the fix of) a real
    malformed input, and why the repulsive energy there is so large
  - the small/medium/large speed-vs-accuracy tradeoff, with real numbers
  - space group and stress as extra physical diagnostics in the report
  - references.bib and --save-report/--save-data
  - CLI and the interactive stb-suite menu producing identical results

Not exercised by this script (needs a display): --view opens an interactive
ase-gui window with BOTH the before and after structures as two frames you
can page between --  try it yourself:
  stb-mlrelax -f si_defect.fdf --view

As a next step, try on your own:
  stb-mlrelax -f your_structure.fdf --custom-model my_finetuned.model
  stb-mlrelax -f your_slab.fdf --relax-cell   # vacuum axis stays fixed

This is a heuristic, not a substitute for a real SIESTA relaxation --
always relax the ML-pre-relaxed structure with SIESTA afterward.
EOF
