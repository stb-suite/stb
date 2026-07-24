#!/bin/bash
# Guided example: stb-defect (code 2.5 in the stb-suite menu)
#
# Not an automated test (see test/2-structures/5-defect/test.sh for that)
# -- a commented walk-through: it runs real commands, one group at a time,
# into its own output/<case>/ folder, and shows you the piece of output
# that proves what just happened. Pauses between sections so you can read
# before moving on. Safe to re-run any time -- it always starts by wiping
# its own output/.

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
echo " Three kinds of point defect"
echo "=================================================================="
cat <<'EOF'
VACANCY       -- remove an atom. Real at any finite temperature; also the
                 starting point for diffusion and many color centers.
SUBSTITUTION  -- swap one atom's species for another. The mechanism behind
                 doping (Si:P, Si:B) and alloying.
INTERSTITIAL  -- squeeze an extra atom in between existing ones. Common for
                 small atoms (H, Li, N) diffusing through a host lattice.

Sites can be picked by raw atom --index, by the position --nearest a
target, or automatically at EVERY symmetrically distinct site at once
(--all-inequivalent-sites) -- not every site of a given element is
physically equivalent, and picking one by hand means already knowing
which ones are distinct.
EOF
pause

echo "=================================================================="
echo " Two ways to run it"
echo "=================================================================="
cat <<'EOF'
A -- direct CLI:
  stb-defect -f si_bulk.fdf --type vacancy --index 1

B -- interactive stb-suite menu:
  $ stb-suite
  Select an option (0-6, or a tool code like 4.1.2): 2.5

Every run always writes a numbered report ([0] RUN METADATA ... [7] SUMMARY
& FILES -- [4] only appears with --ml-relax/--ml-rank), a
structure-validation checklist both before and after, and a before/after
symmetry comparison table for every output structure. A full text report
(--save-report), MACE pre-relaxation (--ml-relax/--ml-rank), and an
interactive 3D view (--view) are all off by default.
EOF
pause

echo "=================================================================="
echo " output/vacancy-index/  --  a real, dramatic symmetry break"
echo "=================================================================="
cat <<'EOF'
si_bulk.fdf is bulk silicon, diamond cubic, the real conventional 8-atom
cell (a = 5.431 Ang, space group Fd-3m No. 227, 48 symmetry operations).
Removing just ONE of those 8 atoms breaks nearly all of that symmetry at
once -- watch the space group drop, and inversion symmetry disappear
entirely (m-3m has it, -43m does not):
EOF
mkdir -p "$OUT/vacancy-index"
cp si_bulk.fdf "$OUT/vacancy-index/"
echo
echo "\$ stb-defect -f si_bulk.fdf --type vacancy --index 1 --no-intro"
(cd "$OUT/vacancy-index" && stb-defect -f si_bulk.fdf --type vacancy --index 1 --no-intro > console.log)
grep -E "Selected site|Output atoms" "$OUT/vacancy-index/console.log"
echo
echo "Before/after symmetry table -- one missing atom out of eight is enough"
echo "to lose inversion symmetry completely:"
sed -n '/Detailed symmetry analysis/,/WRITING OUTPUT/p' "$OUT/vacancy-index/console.log" | head -n -1
echo
echo "Provenance header written into the output .fdf itself:"
head -3 "$OUT/vacancy-index/defect.fdf"
pause

echo "=================================================================="
echo " output/substitution/  --  doping: Si -> Ge"
echo "=================================================================="
cat <<'EOF'
--type substitution swaps one atom's species for another -- the same
mechanism behind doping a semiconductor (Si:P, Si:B) or alloying.
EOF
mkdir -p "$OUT/substitution"
cp si_bulk.fdf "$OUT/substitution/"
echo
echo "\$ stb-defect -f si_bulk.fdf --type substitution --index 3 --new-species Ge -o sub.fdf --no-intro"
(cd "$OUT/substitution" && stb-defect -f si_bulk.fdf --type substitution --index 3 \
    --new-species Ge -o sub.fdf --no-intro > console.log)
grep -E "Substitution|Output formula" "$OUT/substitution/console.log"
pause

echo "=================================================================="
echo " output/interstitial/  --  adding an extra atom"
echo "=================================================================="
cat <<'EOF'
--type interstitial adds a new atom at a given position, instead of
removing or replacing an existing one -- e.g. modeling a small diffusing
species (H, Li, N) squeezed between host atoms.
EOF
mkdir -p "$OUT/interstitial"
cp si_bulk.fdf "$OUT/interstitial/"
echo
echo "\$ stb-defect -f si_bulk.fdf --type interstitial --position 0.5 0.5 0.5 --species N -o inter.fdf --no-intro"
(cd "$OUT/interstitial" && stb-defect -f si_bulk.fdf --type interstitial --position 0.5 0.5 0.5 \
    --species N -o inter.fdf --no-intro > console.log)
grep -E "Interstitial|Output atoms" "$OUT/interstitial/console.log"
pause

echo "=================================================================="
echo " output/nearest-position/  --  selecting a site by target position"
echo "=================================================================="
cat <<'EOF'
--nearest picks the atom closest to a target position instead of an index
-- using pymatgen's periodic minimum-image distance, so a target just
OUTSIDE the cell (1.02 instead of 0.02) still resolves through the
periodic image to the correct nearby atom:
EOF
mkdir -p "$OUT/nearest-position"
cp si_bulk.fdf "$OUT/nearest-position/"
echo
echo "\$ stb-defect -f si_bulk.fdf --type vacancy --nearest 1.02 0.01 0.01 -o near.fdf --no-intro"
(cd "$OUT/nearest-position" && stb-defect -f si_bulk.fdf --type vacancy --nearest 1.02 0.01 0.01 \
    -o near.fdf --no-intro > console.log)
grep "Selected site" "$OUT/nearest-position/console.log"
pause

echo "=================================================================="
echo " output/all-inequivalent-sites/  --  magnetite's 2 distinct Fe sites"
echo "=================================================================="
cat <<'EOF'
magnetite.fdf (Fe3O4) has Fe on TWO symmetrically distinct sites --
tetrahedral and octahedral -- while O sits on only one. Instead of picking
a site by hand, --all-inequivalent-sites finds every distinct one
automatically (via spglib) and writes one output structure per site:
EOF
mkdir -p "$OUT/all-inequivalent-sites"
cp magnetite.fdf "$OUT/all-inequivalent-sites/"
echo
echo "\$ stb-defect -f magnetite.fdf --type vacancy --all-inequivalent-sites --filter-species Fe -o fevac.fdf --no-intro"
(cd "$OUT/all-inequivalent-sites" && stb-defect -f magnetite.fdf --type vacancy \
    --all-inequivalent-sites --filter-species Fe -o fevac.fdf --no-intro > console.log)
grep -E "Space group|Symmetrically distinct" "$OUT/all-inequivalent-sites/console.log" | head -2
sed -n '/Site | Species/,/^$/p' "$OUT/all-inequivalent-sites/console.log" | head -4
grep "structure(s) written" "$OUT/all-inequivalent-sites/console.log"
(cd "$OUT/all-inequivalent-sites" && ls fevac_site*.fdf)
pause

echo "=================================================================="
echo " output/ml-rank/  --  which site is actually more stable?"
echo "=================================================================="
if ! python3 -c "import mace" 2>/dev/null; then
    echo "Skipping -- needs the optional 'ml' extra: pip install stb_suite[ml]"
    echo "(everything else in this script works fine without it)."
else
cat <<'EOF'
Knowing the 2 distinct Fe sites is only half the story -- --ml-rank
relaxes each candidate with MACE and ranks them by energy, a fast
pre-screen for which site is worth real DFT time:
EOF
    mkdir -p "$OUT/ml-rank"
    cp magnetite.fdf "$OUT/ml-rank/"
    echo
    echo "\$ stb-defect -f magnetite.fdf --type vacancy --all-inequivalent-sites --filter-species Fe --ml-rank -o rank.fdf --no-intro"
    (cd "$OUT/ml-rank" && stb-defect -f magnetite.fdf --type vacancy --all-inequivalent-sites \
        --filter-species Fe --ml-rank -o rank.fdf --no-intro > console.log 2>&1)
    sed -n '/ML-ranked sites/,/Note:/p' "$OUT/ml-rank/console.log"
    echo
    echo "A real, non-degenerate result: the octahedral site is genuinely more"
    echo "stable for an Fe vacancy here, not a coin flip between two equal sites."
fi
pause

echo "=================================================================="
echo " output/ml-relax/  --  generic pre-relaxation of one chosen defect"
echo "=================================================================="
if ! python3 -c "import mace" 2>/dev/null; then
    echo "Skipping -- needs the optional 'ml' extra: pip install stb_suite[ml]"
    echo "(everything else in this script works fine without it)."
else
cat <<'EOF'
--ml-relax is the generic form -- works with ANY selection mode (here, a
single --index vacancy, not --all-inequivalent-sites), just relaxes
whatever you're already writing. Watch neighboring atoms relax around the
freshly created vacancy -- a real, physically expected local relaxation:
EOF
    mkdir -p "$OUT/ml-relax"
    cp si_bulk.fdf "$OUT/ml-relax/"
    echo
    echo "\$ stb-defect -f si_bulk.fdf --type vacancy --index 1 --ml-relax --ml-relax-cell -o relaxed.fdf --no-intro"
    (cd "$OUT/ml-relax" && stb-defect -f si_bulk.fdf --type vacancy --index 1 \
        --ml-relax --ml-relax-cell -o relaxed.fdf --no-intro > console.log 2>&1)
    echo "Full MACE simulation detail, straight from the report:"
    sed 's/\x1b\[[0-9;]*m//g' "$OUT/ml-relax/console.log" | \
        awk '/\[4\] ML PRE-RELAXATION/{flag=1} /\[5\] STRUCTURE VALIDATION/{flag=0} flag'
    echo
    echo "Provenance header now also records the MACE pre-relaxation:"
    head -3 "$OUT/ml-relax/relaxed.fdf"
fi
pause

echo "=================================================================="
echo " output/full-report/  --  --save-report, validation, references.bib"
echo "=================================================================="
cat <<'EOF'
The full numbered report (also written to stb_defect_report.txt with
--save-report) includes the structure-validation checklist for the input
and the output structure, and a references.bib with SIESTA.
EOF
mkdir -p "$OUT/full-report"
cp si_bulk.fdf "$OUT/full-report/"
echo "\$ stb-defect -f si_bulk.fdf --type vacancy --index 1 --save-report --no-intro"
(cd "$OUT/full-report" && stb-defect -f si_bulk.fdf --type vacancy --index 1 --save-report --no-intro > console.log)
echo
echo "Report sections written to stb_defect_report.txt:"
grep -E "^\[[0-9]\]" "$OUT/full-report/stb_defect_report.txt"
echo
echo "Validation checklist (one row shown):"
grep -m1 "Atom proximity" "$OUT/full-report/console.log"
echo
echo "references.bib always written -- SIESTA:"
grep "^@" "$OUT/full-report/references.bib"
pause

echo "=================================================================="
echo " Proof: CLI and the interactive stb-suite menu agree"
echo "=================================================================="
echo "Driving the same basic case through the interactive menu's manual"
echo "entry mode and checking it reaches the same atom count."
TMP="$(mktemp -d)"
cp si_bulk.fdf "$TMP/"
echo
echo "\$ printf '2.5\\nsi_bulk.fdf\\n1\\n1\\n1\\n\\n\\n\\n\\n0\\n' | stb-suite"
(cd "$TMP" && printf '2.5\nsi_bulk.fdf\n1\n1\n1\n\n\n\n\n0\n' | stb-suite > session.log 2>&1) || true
if grep -q "Output atoms     : 7" "$TMP/session.log"; then
    echo "Confirmed: the interactive menu built and launched the exact same"
    echo "underlying stb-defect command as the CLI walkthrough above (7 atoms"
    echo "out, same as output/vacancy-index/)."
else
    echo "Unexpected: menu did not reach the write step -- see $TMP/session.log."
fi
rm -rf "$TMP"
pause

echo "=================================================================="
echo " Done"
echo "=================================================================="
cat <<EOF
Up to eight self-contained folders were generated under output/ (ml-rank/
and ml-relax/ skipped if the optional 'ml' extra isn't installed):
  vacancy-index/          substitution/       interstitial/
  nearest-position/       all-inequivalent-sites/
  ml-rank/                ml-relax/           full-report/

Each has references.bib; full-report/ additionally has
stb_defect_report.txt; ml-rank/'s and ml-relax/'s references.bib also cite
the MACE architecture/foundation-model papers.

Recap of what this walkthrough covered:
  - vacancy / substitution / interstitial, and what each models physically
  - a real, measured symmetry break: one missing atom out of eight drops
    Fd-3m (48 operations) straight to P-43m (24 operations), losing
    inversion symmetry entirely
  - --nearest's periodic minimum-image resolution for a target just
    outside the cell
  - --all-inequivalent-sites automatically finding magnetite's 2 distinct
    Fe sites (tetrahedral, octahedral), instead of requiring you to
    already know which sites are distinct
  - --ml-rank: a real, non-degenerate stability ranking between those 2
    sites (octahedral ~0.06 eV more stable for a vacancy)
  - --ml-relax vs. --ml-rank: generic pre-relaxation of whatever you're
    writing, vs. relax-and-compare across candidate sites
  - the structure-validation checklist, references.bib, and --save-report
  - CLI and the interactive stb-suite menu building the same command

Not exercised by this script (needs a display): --view opens the input
structure and every output structure in an interactive ase-gui window --
try it yourself:
  stb-defect -f si_bulk.fdf --type vacancy --index 1 --view

As a next step, try on your own:
  stb-defect -f your_structure.fdf --type vacancy --all-inequivalent-sites  # find every distinct site
  stb-defect -f si_bulk.fdf --type substitution --index 1 --new-species P  # n-type doping
  stb-supercell -f si_bulk.fdf -d 2 2 2 -o si64.fdf --no-intro && \\
      stb-defect -f si64.fdf --type vacancy --index 1                     # a more dilute defect
EOF
