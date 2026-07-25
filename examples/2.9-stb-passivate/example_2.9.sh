#!/bin/bash
# Guided example: stb-passivate (code 2.9 in the stb-suite menu)
#
# Not an automated test (see test/2-structures/9-passivate/test.sh for that)
# -- a commented walk-through: it runs real commands, one group at a time,
# into its own output/<case>/ folder, and shows you the piece of output
# that proves what just happened. Pauses between sections so you can read
# before moving on. Safe to re-run any time -- it always starts by wiping
# its own output/. No committed input fixtures: like 2.8-stb-crystalbuilder's
# own bulk-graphite-to-slab/ case, this script builds the bulk crystal and
# cuts the slabs itself, live, at the start.

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

strip_ansi() {
    sed 's/\x1b\[[0-9;]*m//g' "$1"
}

echo "=================================================================="
echo " Why does cutting a crystal leave dangling bonds?"
echo "=================================================================="
cat <<'EOF'
A bulk crystal's atoms are all fully coordinated -- every bond they
"should" have is actually there. Cutting a slab out of that bulk (to model
a surface) necessarily severs some of those bonds: the atoms right at the
new surface end up with FEWER neighbors than the bulk value. In a real
material this gets satisfied somehow (reconstruction, adsorbates, an
oxide layer); a bare idealized DFT slab has none of that -- the dangling
bonds are unphysical, and introduce spurious mid-gap surface states that
can contaminate the very electronic structure (gap, DOS, work function)
you built the slab to study.

stb-passivate caps each dangling bond with a passivating atom (default H)
-- cheap, chemically simple, and enough to remove the spurious surface
states so the slab's interior behaves like the real bulk.
EOF
pause

echo "=================================================================="
echo " Setup: build bulk Si, then cut BOTH Si(111) terminations"
echo "=================================================================="
cat <<'EOF'
This example builds its own inputs live, the same 2.8 -> 2.3 pipeline a
real user would run first: bulk Si (space group Fd-3m, diamond cubic) via
stb-crystalbuilder, reduced to the 8-atom conventional cell, then cut
along (111) with stb-slab -- which finds TWO distinct terminations of the
exact same surface, cutting through two different atomic planes:
EOF
mkdir -p "$OUT/_setup"
echo
echo "\$ stb-crystalbuilder --spacegroup Fd-3m --a 5.43 --site Si 0 0 0 --reduce conventional -o si_bulk.fdf --no-intro"
(cd "$OUT/_setup" && stb-crystalbuilder --spacegroup Fd-3m --a 5.43 --site Si 0 0 0 \
    --reduce conventional -o si_bulk.fdf --no-intro > build.log)
echo "\$ stb-slab -f si_bulk.fdf --hkl 1 1 1 --all -o term.fdf --no-intro"
(cd "$OUT/_setup" && stb-slab -f si_bulk.fdf --hkl 1 1 1 --all -o term.fdf --no-intro > slab_all.log)
strip_ansi "$OUT/_setup/slab_all.log" | grep -E "Terminations found|Total atoms|Polar|Symmetric"
pause

echo "=================================================================="
echo " output/good-termination/  --  the clean case: 8 sites, 1 bond each"
echo "=================================================================="
cat <<'EOF'
Termination 0 cuts through the WEAK gap between Si bilayers -- each of the
8 surface atoms is missing exactly 1 bond (coordination 3 instead of the
bulk value 4). This is EXACTLY determined geometrically: watch every
capped H land at precisely 109.47 degrees from Si's existing bonds -- the
textbook ideal tetrahedral angle, not an approximation:
EOF
mkdir -p "$OUT/good-termination"
cp "$OUT/_setup/term_term0.fdf" "$OUT/good-termination/slab.fdf"
echo
echo "\$ stb-passivate -f slab.fdf -o passivated.fdf --no-intro"
(cd "$OUT/good-termination" && stb-passivate -f slab.fdf -o passivated.fdf --no-intro > console.log)
grep -E "Dangling bonds found|Auto-passivated|Output formula|Output atoms" "$OUT/good-termination/console.log"
echo
echo "Verifying the H-Si-Si bond angles in the output geometry:"
python3 - "$OUT/good-termination/passivated.fdf" <<'PYEOF'
import sys, itertools
import numpy as np
from stb.core import structure_io
s = structure_io.to_pymatgen(structure_io.read_fdf(sys.argv[1]))
h = next(site for site in s if site.specie.symbol == "H")
si_sites = [site for site in s if site.specie.symbol == "Si"]
host = min(si_sites, key=lambda si: si.distance(h))
neighbors = [n for n in s.get_neighbors(host, 2.6) if n.specie.symbol == "Si"]
vecs = [(n.coords - host.coords) for n in neighbors] + [h.coords - host.coords]
vecs = [v / np.linalg.norm(v) for v in vecs]
angles = [np.degrees(np.arccos(np.clip(np.dot(a, b), -1, 1))) for a, b in itertools.combinations(vecs, 2)]
print(f"  All {len(angles)} bond angles at the capped site: {[round(float(a), 2) for a in angles]} deg")
PYEOF
pause

echo "=================================================================="
echo " output/bad-termination/  --  same surface, one plane over: 3 bonds missing"
echo "=================================================================="
cat <<'EOF'
Termination 1 cuts through the STRONG interior of a Si bilayer instead --
the exact same (111) surface, just a different atomic plane. Now each of
the 8 surface atoms is missing 3 bonds, not 1. With only 1 remaining bond
vector, the missing directions are genuinely ambiguous (infinitely many
completions give the right coordination number) -- watch stb-passivate
correctly refuse to guess, reporting every site instead:
EOF
mkdir -p "$OUT/bad-termination"
cp "$OUT/_setup/term_term1.fdf" "$OUT/bad-termination/slab.fdf"
echo
echo "\$ stb-passivate -f slab.fdf -o passivated.fdf --no-intro"
(cd "$OUT/bad-termination" && stb-passivate -f slab.fdf -o passivated.fdf --no-intro > console.log)
grep -E "Dangling bonds found|Auto-passivated|WARNING|deficit|Output formula|Output atoms" "$OUT/bad-termination/console.log" | head -6
pause

echo "=================================================================="
echo " output/bond-length-override/  --  the default is an approximation"
echo "=================================================================="
cat <<'EOF'
The default bond length is the sum of the two species' atomic radii (Si+H
= 1.35 Ang) -- fast and generic, but not the real experimental Si-H value
(~1.48 Ang). --bond-length lets you dial in the real number:
EOF
mkdir -p "$OUT/bond-length-override"
cp "$OUT/_setup/term_term0.fdf" "$OUT/bond-length-override/slab.fdf"
echo
echo "\$ stb-passivate -f slab.fdf -o auto.fdf --no-intro                          # default (auto)"
(cd "$OUT/bond-length-override" && stb-passivate -f slab.fdf -o auto.fdf --no-intro > auto.log)
echo "\$ stb-passivate -f slab.fdf --bond-length 1.48 -o real.fdf --no-intro        # real experimental value"
(cd "$OUT/bond-length-override" && stb-passivate -f slab.fdf --bond-length 1.48 -o real.fdf --no-intro > real.log)
python3 - "$OUT/bond-length-override/auto.fdf" "$OUT/bond-length-override/real.fdf" <<'PYEOF'
import sys
from stb.core import structure_io
for label, path in zip(("Default (auto)", "--bond-length 1.48"), sys.argv[1:]):
    s = structure_io.to_pymatgen(structure_io.read_fdf(path))
    h = next(site for site in s if site.specie.symbol == "H")
    si_sites = [site for site in s if site.specie.symbol == "Si"]
    nearest = min(h.distance(si) for si in si_sites)
    print(f"  {label:<20}: Si-H = {nearest:.4f} Ang")
PYEOF
pause

echo "=================================================================="
echo " output/symmetry-preserved/  --  the tool's own before/after proof"
echo "=================================================================="
cat <<'EOF'
Since each capped bond is placed exactly along the direction the
structure's own local symmetry dictates, passivation (in the clean,
single-missing-bond case) never breaks the crystal's overall symmetry.
Every run's own [6] SYMMETRY ANALYSIS table proves this directly --
crystal system, 3D space group, LAYER group (meaningful here, since the
input is a genuinely 2D-periodic slab -- kept in this tool's own table,
unlike stb-crystalbuilder's bulk-only one), point group, and Hall symbol
all come back identical:
EOF
mkdir -p "$OUT/symmetry-preserved"
cp "$OUT/_setup/term_term0.fdf" "$OUT/symmetry-preserved/slab.fdf"
echo
echo "\$ stb-passivate -f slab.fdf --save-report -o passivated.fdf --no-intro"
(cd "$OUT/symmetry-preserved" && stb-passivate -f slab.fdf --save-report -o passivated.fdf --no-intro > console.log)
strip_ansi "$OUT/symmetry-preserved/stb_passivate_report.txt" | awk '/\[6\] SYMMETRY ANALYSIS/{flag=1} /\[7\] WRITING OUTPUT/{flag=0} flag'
pause

echo "=================================================================="
echo " output/ml-relax/  --  MACE pre-relaxation, vacuum axis fixed"
echo "=================================================================="
if ! python3 -c "import mace" 2>/dev/null; then
    echo "Skipping -- needs the optional 'ml' extra: pip install stb_suite[ml]"
    echo "(everything else in this script works fine without it)."
else
cat <<'EOF'
The capped H position is a geometric heuristic (an ideal, rigid bonding
direction) -- a real system relaxes slightly from there. --ml-relax
(+ --ml-relax-cell) corrects this cheaply with MACE before a real SIESTA
relaxation. Vacuum-aware: the vacuum axis (c) always stays EXACTLY fixed,
even with --ml-relax-cell -- verified directly below, not just assumed:
EOF
    mkdir -p "$OUT/ml-relax"
    cp "$OUT/_setup/term_term0.fdf" "$OUT/ml-relax/slab.fdf"
    echo
    echo "\$ stb-passivate -f slab.fdf --ml-relax --ml-relax-cell -o relaxed.fdf --no-intro"
    (cd "$OUT/ml-relax" && stb-passivate -f slab.fdf --ml-relax --ml-relax-cell \
        -o relaxed.fdf --no-intro > console.log 2>&1)
    echo "Full MACE simulation detail, straight from the report:"
    strip_ansi "$OUT/ml-relax/console.log" | awk '/\[4\] ML PRE-RELAXATION/{flag=1} /\[5\] STRUCTURE VALIDATION/{flag=0} flag'
    echo
    python3 - "$OUT/ml-relax/slab.fdf" "$OUT/ml-relax/relaxed.fdf" <<'PYEOF'
import sys
from stb.core import structure_io
before = structure_io.to_pymatgen(structure_io.read_fdf(sys.argv[1]))
after = structure_io.to_pymatgen(structure_io.read_fdf(sys.argv[2]))
print(f"Vacuum axis (c): {before.lattice.c:.4f} Ang -> {after.lattice.c:.4f} Ang "
      f"(unchanged: {abs(before.lattice.c - after.lattice.c) < 1e-6})")
PYEOF
    echo
    echo "Provenance header records the full history:"
    head -3 "$OUT/ml-relax/relaxed.fdf"
fi
pause

echo "=================================================================="
echo " output/full-report/  --  --save-report, validation, references.bib"
echo "=================================================================="
cat <<'EOF'
The full numbered report (also written to stb_passivate_report.txt with
--save-report) includes the structure-validation checklist and a
references.bib with SIESTA:
EOF
mkdir -p "$OUT/full-report"
cp "$OUT/_setup/term_term0.fdf" "$OUT/full-report/slab.fdf"
echo "\$ stb-passivate -f slab.fdf --save-report --no-intro"
(cd "$OUT/full-report" && stb-passivate -f slab.fdf --save-report --no-intro > console.log)
echo
echo "Report sections written to stb_passivate_report.txt:"
grep -E "^\[[0-9]+\]" "$OUT/full-report/stb_passivate_report.txt"
echo
echo "references.bib -- SIESTA:"
grep "^@" "$OUT/full-report/references.bib"
pause

echo "=================================================================="
echo " Proof: CLI and the interactive stb-suite menu agree"
echo "=================================================================="
echo "Driving the same good-termination case through the interactive"
echo "menu's manual entry mode and checking it reaches the same result."
TMP="$(mktemp -d)"
cp "$OUT/_setup/term_term0.fdf" "$TMP/slab.fdf"
echo
echo "\$ printf '2.9\\ndoes_not_exist.fdf\\nslab.fdf\\n\\n\\n\\n\\n\\n\\n\\n0\\n' | stb-suite"
(cd "$TMP" && printf '2.9\ndoes_not_exist.fdf\nslab.fdf\n\n\n\n\n\n\n\n0\n' | stb-suite > session.log 2>&1) || true
if grep -q "Auto-passivated      : 8 with H" "$TMP/session.log"; then
    echo "Confirmed: the interactive menu built and launched the exact same"
    echo "underlying stb-passivate command as the CLI walkthrough above"
    echo "(8 sites auto-passivated, same as output/good-termination/)."
else
    echo "Unexpected: menu did not reach the write step -- see $TMP/session.log."
fi
rm -rf "$TMP"
pause

echo "=================================================================="
echo " Done"
echo "=================================================================="
cat <<EOF
Up to six self-contained folders were generated under output/ (ml-relax/
skipped if the optional 'ml' extra isn't installed):
  good-termination/       bad-termination/       bond-length-override/
  symmetry-preserved/     ml-relax/              full-report/

Each has references.bib (SIESTA, always); full-report/ additionally has
stb_passivate_report.txt; ml-relax/'s references.bib also cites the MACE
architecture/foundation-model papers.

Recap of what this walkthrough covered:
  - why cutting a crystal leaves unphysical dangling bonds, and why that
    matters for a real DFT slab calculation
  - the exact geometry for a single missing bond (109.47 deg tetrahedral,
    verified directly) vs. the genuine ambiguity of 2+ missing bonds
  - the SAME physical surface, cut through 2 different atomic planes,
    giving completely different (and correctly handled) outcomes
  - --bond-length: the atomic-radii-sum default vs. dialing in a real
    experimental bond length
  - the tool's own before/after symmetry table proving passivation
    preserves the crystal's full symmetry in the clean case
  - --ml-relax correcting the geometric heuristic, with the vacuum axis
    provably unchanged
  - the structure-validation checklist, references.bib, and --save-report
  - CLI and the interactive stb-suite menu building the same command

Not exercised by this script (needs a display): --view opens the input
and passivated structures in an interactive ase-gui window -- try it
yourself:
  stb-passivate -f output/good-termination/slab.fdf --view

As a next step, try on your own:
  stb-crystalbuilder --spacegroup Fd-3m --a 5.43 --site Si 0 0 0 -o si_bulk.fdf
  stb-slab -f si_bulk.fdf --hkl 1 1 1 -o slab.fdf
  stb-passivate -f slab.fdf --bond-length 1.48 --ml-relax --ml-relax-cell
EOF
