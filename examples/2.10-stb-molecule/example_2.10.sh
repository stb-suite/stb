#!/bin/bash
# Guided example: stb-molecule (code 2.10 in the stb-suite menu)
#
# Not an automated test (see test/2-structures/10-molecule/test.sh for
# that) -- a commented walk-through: it runs real commands, one group at a
# time, into its own output/<case>/ folder, and shows you the piece of
# output that proves what just happened. Pauses between sections so you
# can read before moving on. Safe to re-run any time -- it always starts
# by wiping its own output/. No committed input fixtures: like
# 2.8-stb-crystalbuilder, everything here is built from --name/--vacuum
# flags, not read from a structure file.

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
echo " Why does an 'isolated' molecule need a vacuum box at all?"
echo "=================================================================="
cat <<'EOF'
SIESTA (like almost every plane-wave/localized-basis DFT code) always
solves a PERIODIC crystal -- there is no native "just one molecule,
nowhere else" mode. To model an isolated molecule anyway, you place it in
a large, empty box and let THAT repeat periodically -- far enough apart
that the molecule's own periodic images don't meaningfully interact, so
the calculation behaves as if it were truly alone.

stb-molecule builds exactly this: a real, chemically sensible molecule
(from ASE's G2 database -- the "G2/97" computational-thermochemistry
benchmark set, not an arbitrary list) placed in a vacuum box of your
chosen size (--vacuum, default 10 Ang in every direction).
EOF
pause

echo "=================================================================="
echo " Two ways to run it"
echo "=================================================================="
cat <<'EOF'
A -- direct CLI:
  stb-molecule --name H2O

B -- interactive stb-suite menu:
  $ stb-suite
  Select an option (0-6, or a tool code like 4.1.2): 2.10

Every run always writes a numbered report ([0] RUN METADATA ... [8]
SUMMARY & FILES -- [3] only appears with --ml-relax) and a point-group
symmetry check before and after. A full text report (--save-report), MACE
pre-relaxation (--ml-relax), and an interactive 3D view (--view) are all
off by default.
EOF
pause

echo "=================================================================="
echo " output/point-group-gallery/  --  4 real, textbook point groups"
echo "=================================================================="
cat <<'EOF'
A crystal's symmetry is a SPACE group (2.8-stb-crystalbuilder); an
isolated molecule has no periodicity, so it's classified by its POINT
group instead (Schoenflies notation) -- a completely different detector
(pymatgen's PointGroupAnalyzer) from every periodic tool in this suite.
Four real molecules, four real point groups, verified directly from the
built geometry -- not just quoted from a textbook:
EOF
mkdir -p "$OUT/point-group-gallery"
for name in H2O CH4 C6H6 CO2; do
    echo
    echo "\$ stb-molecule --name $name -o ${name}.fdf --no-intro"
    (cd "$OUT/point-group-gallery" && stb-molecule --name "$name" -o "${name}.fdf" \
        --no-intro > "log_${name}.txt")
    point_group=$(grep -E "^Point group" "$OUT/point-group-gallery/log_${name}.txt" | head -1 | \
        sed 's/^Point group    : //')
    echo "  ${name}: ${point_group}"
done
echo
echo "(water bent C2v, methane tetrahedral Td, benzene hexagonal D6h,"
echo " CO2 linear D*h -- pymatgen's notation for the infinite-order D-inf-h)"
pause

echo "=================================================================="
echo " output/unknown-name-suggestions/  --  helpful errors, not just failures"
echo "=================================================================="
cat <<'EOF'
Names are case-sensitive and there's no fuzzy matching by default -- but
a wrong name still gets a helpful, specific suggestion instead of just
"not found":
EOF
mkdir -p "$OUT/unknown-name-suggestions"
echo
echo "\$ stb-molecule --name h2o --no-intro                # wrong case"
(cd "$OUT/unknown-name-suggestions" && stb-molecule --name h2o --no-intro > wrong_case.log 2>&1) || true
grep "ERROR" "$OUT/unknown-name-suggestions/wrong_case.log"
echo
echo "\$ stb-molecule --name CH3OJ --no-intro               # genuine typo"
(cd "$OUT/unknown-name-suggestions" && stb-molecule --name CH3OJ --no-intro > typo.log 2>&1) || true
grep "ERROR" "$OUT/unknown-name-suggestions/typo.log"
echo
echo "\$ stb-molecule --list --no-intro | head -2           # see all 162 names"
stb-molecule --list --no-intro | head -2
pause

echo "=================================================================="
echo " output/vacuum-and-dimensionality/  --  a real, verified pitfall"
echo "=================================================================="
cat <<'EOF'
stb-molecule always builds something physically 0D (fully isolated), but
the report's own "Dimensionality" line is a generic heuristic shared with
every periodic-aware tool in this suite: it checks whether the gap
between periodic images exceeds a fixed 10 Ang threshold on each axis.
The DEFAULT --vacuum 10.0 sits right at that threshold and reads
correctly. Turn it down, and watch the SAME physically-isolated molecule
get misclassified -- while the POINT GROUP (which never looks at the box
at all) stays exactly right either way:
EOF
mkdir -p "$OUT/vacuum-and-dimensionality"
echo
echo "\$ stb-molecule --name H2O -o default.fdf --no-intro                 # default vacuum=10"
(cd "$OUT/vacuum-and-dimensionality" && stb-molecule --name H2O -o default.fdf --no-intro > default.log)
echo "\$ stb-molecule --name H2O --vacuum 3 -o small.fdf --no-intro        # vacuum=3, below the threshold"
(cd "$OUT/vacuum-and-dimensionality" && stb-molecule --name H2O --vacuum 3 -o small.fdf --no-intro > small.log)
echo
echo "Default (--vacuum 10):"
grep -E "Dimensionality|Point group" "$OUT/vacuum-and-dimensionality/default.log" | head -2
echo "Turned down (--vacuum 3):"
grep -E "Dimensionality|Point group|Atomic density" "$OUT/vacuum-and-dimensionality/small.log" | head -3
pause

echo "=================================================================="
echo " output/ml-relax-symmetry-preserved/  --  MACE relax, Td unchanged"
echo "=================================================================="
if ! python3 -c "import mace" 2>/dev/null; then
    echo "Skipping -- needs the optional 'ml' extra: pip install stb_suite[ml]"
    echo "(everything else in this script works fine without it)."
else
cat <<'EOF'
The G2 database's own reference geometries are already close to
equilibrium, so relaxing with MACE is expected to correct bond lengths/
angles slightly WITHOUT changing the molecule's shape -- watch methane's
tetrahedral Td point group survive the relaxation exactly, proven by the
tool's own before/after table, not just assumed:
EOF
    mkdir -p "$OUT/ml-relax-symmetry-preserved"
    echo
    echo "\$ stb-molecule --name CH4 --ml-relax -o relaxed.fdf --no-intro"
    (cd "$OUT/ml-relax-symmetry-preserved" && stb-molecule --name CH4 --ml-relax \
        -o relaxed.fdf --no-intro > console.log 2>&1)
    echo "Full MACE simulation detail, straight from the report:"
    strip_ansi "$OUT/ml-relax-symmetry-preserved/console.log" | \
        awk '/\[3\] ML PRE-RELAXATION/{flag=1} /\[4\] STRUCTURE VALIDATION/{flag=0} flag'
    echo
    strip_ansi "$OUT/ml-relax-symmetry-preserved/console.log" | \
        awk '/\[5\] SYMMETRY ANALYSIS/{flag=1} /\[6\] WRITING OUTPUT/{flag=0} flag'
fi
pause

echo "=================================================================="
echo " output/full-report/  --  --save-report, validation, references.bib"
echo "=================================================================="
cat <<'EOF'
The full numbered report (also written to stb_molecule_report.txt with
--save-report) includes the structure-validation checklist and a
references.bib with SIESTA:
EOF
mkdir -p "$OUT/full-report"
echo "\$ stb-molecule --name H2O --save-report --no-intro"
(cd "$OUT/full-report" && stb-molecule --name H2O --save-report --no-intro > console.log)
echo
echo "Report sections written to stb_molecule_report.txt:"
grep -E "^\[[0-9]+\]" "$OUT/full-report/stb_molecule_report.txt"
echo
echo "Provenance header written into the output .fdf:"
head -3 "$OUT/full-report/molecule.fdf"
echo
echo "references.bib -- SIESTA:"
grep "^@" "$OUT/full-report/references.bib"
pause

echo "=================================================================="
echo " Proof: CLI and the interactive stb-suite menu agree"
echo "=================================================================="
echo "Driving the same H2O case through the interactive menu's manual"
echo "entry mode and checking it reaches the same point group."
TMP="$(mktemp -d)"
echo
echo "\$ printf '2.10\\nH2O\\n\\n\\n\\n\\n\\n\\n0\\n' | stb-suite"
(cd "$TMP" && printf '2.10\nH2O\n\n\n\n\n\n\n0\n' | stb-suite > session.log 2>&1) || true
if grep -q "Point group    : C2v" "$TMP/session.log"; then
    echo "Confirmed: the interactive menu built and launched the exact same"
    echo "underlying stb-molecule command as the CLI walkthrough above"
    echo "(point group C2v, same as output/point-group-gallery/log_H2O.txt)."
else
    echo "Unexpected: menu did not reach the write step -- see $TMP/session.log."
fi
rm -rf "$TMP"
pause

echo "=================================================================="
echo " Done"
echo "=================================================================="
cat <<EOF
Up to five self-contained folders were generated under output/ (ml-relax-
symmetry-preserved/ skipped if the optional 'ml' extra isn't installed):
  point-group-gallery/          unknown-name-suggestions/
  vacuum-and-dimensionality/     ml-relax-symmetry-preserved/
  full-report/

Each has references.bib (SIESTA, always); full-report/ additionally has
stb_molecule_report.txt; ml-relax-symmetry-preserved/'s references.bib
also cites the MACE architecture/foundation-model papers.

Recap of what this walkthrough covered:
  - why an "isolated" molecule in a periodic DFT code needs an artificial
    vacuum box, and what the G2 database actually is
  - point groups vs. space groups: a molecule's own symmetry, verified
    directly for 4 real molecules (C2v, Td, D6h, D*h)
  - helpful, specific error messages for a wrong molecule name
  - a real, verified pitfall: --vacuum below the 10 Ang detection
    threshold fools the report's dimensionality line (but never the point
    group, which stays correct regardless)
  - --ml-relax correcting the G2 geometry while providably preserving the
    molecule's point group
  - the structure-validation checklist, references.bib, and --save-report
  - CLI and the interactive stb-suite menu building the same command

Not exercised by this script (needs a display): --view opens the as-built
and final molecule in an interactive ase-gui window -- try it yourself:
  stb-molecule --name H2O --view

As a next step, try on your own:
  stb-molecule --list
  stb-molecule --name <any of the 162 names>
  stb-molecule --name CH3OH --ml-relax --custom-model my_finetuned.model
EOF
