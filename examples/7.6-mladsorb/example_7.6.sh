#!/bin/bash
# Guided example: 7.6 -- ML Adsorption-Site Screening (stb-mladsorb; code 7.6 in the stb-suite menu).
#
# Not an automated test (see test/7-mlsimulations/6-mladsorb/test.sh for that) -- a commented walk-through: it runs
# real commands, one case at a time, into its own output/<case>/ folder, and shows the piece of output that proves what
# just happened. Pauses between cases so you can read before moving on (`yes '' | bash example_7.6.sh` runs it straight
# through). Safe to re-run any time -- it always starts by wiping its own output/.
#
# The system is the one of README section 6: one silicon adatom on the ideal Si(111) surface, a 2x2 slab of three
# bilayers, with the round-3 silicon model of the 7.5 example (../7.5-mlneb/model_round3/). The live cases take about
# 1 minute each on a GPU (DEVICE=cuda) and a few minutes each on a CPU (the default) -- about 15 minutes in all on a
# CPU. Nothing here runs SIESTA: the SIESTA relaxations of README section 7 are in reference_results/siesta/, and
# siesta_compare.py reads them (Case 7). The full reference set of section 6 (relaxation depth, starting heights,
# binding curves, slab size, other models) is read from reference_results/ (Case 6; make_reference.sh regenerates it
# in about 30 minutes on a GPU).
#
# Needs the 'ml' extra: pip install stb_suite[ml]

set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

export MPLBACKEND=Agg           # headless-safe; the PNG figures are written either way
DEVICE="${DEVICE:-cpu}"
R3="$DIR/../7.5-mlneb/model_round3"
ROUND3="--custom-model $R3/mlff_model_compiled.model --custom-model-card $R3/mlff_model_card.json"
REF="$DIR/reference_results"

OUT="$DIR/output"
rm -rf "$OUT"
mkdir -p "$OUT"

pause() {
    echo
    read -p "  [Press Enter to continue] " -r || true
    echo
}

banner() {
    echo "=================================================================="
    echo " $1"
    echo "=================================================================="
}

strip_colors() { sed 's/\x1b\[[0-9;]*m//g'; }

# a report section by its title: from "[N] TITLE..." up to the next "[M] " header. By title, not number: which
# sections are present depends on the options, so their numbers move.
section() {   # $1 = file, $2 = start of the title
    strip_colors < "$1" | awk -v t="$2" '/^\[[0-9]+\] [A-Z]/ {if (on) exit; if (index($0, "] " t) > 0) on = 1}
        on {print}'
}

# ads <case folder> <flags...>: runs stb-mladsorb in output/<case>/ on the 2x2 slab, report saved, log kept
ads() {
    local case=$1; shift
    mkdir -p "$OUT/$case"
    local t0=$SECONDS
    (cd "$OUT/$case" && stb-mladsorb --slab "$DIR/si111_2x2.fdf" --adsorbate Si --device "$DEVICE" \
        -o . --save-report --no-intro "$@" < /dev/null > log.txt 2>&1)
    echo "  (stb-mladsorb took $((SECONDS - t0)) s on $DEVICE)"
}
report() { echo "$OUT/$1/stb_mladsorb_report.txt"; }


banner "Welcome: where does a silicon adatom sit on Si(111)?"
cat <<EOF
stb-mladsorb finds the candidate adsorption sites of a slab, relaxes the adsorbate on each with a MACE potential,
reports where every run ENDED, and groups the runs into distinct minima with their energies, multiplicities and
populations -- with the model card's trained range, the foundation model and the committee alongside.

The system: one Si adatom per 2x2 cell of the ideal Si(111) surface (1/4 monolayer). The answer is known: the adatom
of the 7x7 and 2x2 reconstructions sits on T4, above a second-layer atom (README section 1.5). The tour:
   1  the slab and its sites                          6  the reference set: depth, heights, size, models (stored)
   2  the main screening, top layers relaxed  [run]   7  the SIESTA check (stored)
   3  the files: data, gnuplot, extended XYZ          8  guard rails
   4  a fixed substrate, for comparison       [run]   9  the interactive menu (stb-suite, code 7.6)
   5  the foundation model on the same sites  [run]
Device: $DEVICE. Theory: README section 1. Numbers: README sections 6 and 7.
EOF
pause


banner "Case 1 -- The slab: Si(111) cut between bilayers, and its candidate sites"
cat <<'EOF'
make_slabs.py builds the three slabs of the example from the diamond structure at the 7.1 training lattice constant
(5.480238 Ang): a 2x2 surface cell of three bilayers (24 atoms, the main system), the same with four bilayers, and a 3x3
cell. It drops the single-atom layers diamond111 leaves at both faces, so every surface atom has ONE dangling bond.
EOF
mkdir -p "$OUT/case1"
(cd "$OUT/case1" && python3 "$DIR/make_slabs.py" 2>/dev/null | sed 's/^/    /')
if diff <(grep -v '^#' "$OUT/case1/si111_2x2.fdf") <(grep -v '^#' si111_2x2.fdf) > /dev/null; then
    echo "    si111_2x2.fdf: identical to the file shipped in this folder (used by every run below)"
fi
cat <<'EOF'

The sites pymatgen finds on it, as stb-mladsorb reports them (from the main run of the reference set). The 'what is
below' column is what names them: #1 ontop above L1 is T1, #2 ontop above L2 is T4, #3 bridge above L4 is H3, #5 the
B2 bridge between two surface atoms (README section 1.4, Figure 1):
EOF
section "$REF/main/stb_mladsorb_report.txt" "CANDIDATE SITES"
pause


banner "Case 2 -- The main screening: round-3 model, top 4 Ang relaxed, bulk reference   [run]"
cat <<'EOF'
  stb-mladsorb --slab si111_2x2.fdf --adsorbate Si --custom-model .../model_round3/mlff_model_compiled.model
               --custom-model-card .../mlff_model_card.json --relax-depth 4 --bulk-reference si_bulk.fdf
               --skip-foundation-comparison --no-committee --save-data --save-gnuplot --save-xyz --save-report
(The committee members of the model are not shipped -- 11 MB each; their numbers are in the reference report, Case 6.)
EOF
ads case2 $ROUND3 --relax-depth 4 --bulk-reference "$DIR/si_bulk.fdf" --skip-foundation-comparison --no-committee \
    --save-data --save-gnuplot --save-xyz
R=$(report case2)
section "$R" "CLEAN SLAB"
section "$R" "RELAXATIONS"
section "$R" "DISTINCT MINIMA"
cat <<'EOF'
  -> The slab cut from the bulk relaxes by -2.2 eV. Two of six runs end elsewhere than they started (the cyan rows):
     the bridges #4 and #6 are not minima. Four minima remain; T4 is first, as expected -- but read on: the order of
     the others and every number are what the SIESTA check of Case 7 puts to the test.
EOF
pause


banner "Case 3 -- The files: data tables, gnuplot, extended XYZ"
ls "$OUT/case2" | sed 's/^/    /'
echo
echo "  mladsorb_minima.dat (one line per minimum):"
head -4 "$OUT/case2/mladsorb_minima.dat" | sed 's/^/    /'
echo
echo "  mladsorb_minima.xyz, read back with ASE (each frame carries its energy, forces and the run's data):"
python3 - "$OUT/case2/mladsorb_minima.xyz" <<'EOF'
import sys
from ase.io import read
for a in read(sys.argv[1], ":"):
    i = a.info
    print(f"    rank {i['rank']}: {i['end_type']:<6} site {i['end_site']}  E_ads {i['E_ads']:.4f} eV  E_form "
          f"{i['E_form']:.4f} eV  {len(a)} atoms, adsorbate atoms {list(a.arrays['adsorbate']).count(1)}, "
          f"max|F| {abs(a.get_forces()).max():.3f} eV/Ang")
EOF
echo
if command -v gnuplot > /dev/null; then
    echo "  gnuplot was found: plot/ holds the .dat files, the .gplot scripts and their PDFs:"
else
    echo "  gnuplot is not installed: plot/ holds the .dat files and .gplot scripts (run 'gnuplot top_view.gplot' later):"
fi
ls "$OUT/case2/plot" | sed 's/^/    /'
pause


banner "Case 4 -- The same screening with the substrate fixed as given   [run]"
cat <<'EOF'
--relax-depth 0 (the default): only the adatom moves, the ideal surface stays exactly as cut.
EOF
ads case4 $ROUND3 --relax-depth 0 --bulk-reference "$DIR/si_bulk.fdf" --skip-foundation-comparison --no-committee
section "$(report case4)" "DISTINCT MINIMA" | head -8
section "$(report case4)" "WARNINGS"
cat <<'EOF'
  -> On an unrelaxed surface the adatom cannot pull its three neighbours towards it: T4 is still first, but bound
     about 1 eV less, H3 and B2 change, and the report warns that the forces on the top layer are not small.
EOF
pause


banner "Case 5 -- The foundation model on the same sites   [run]"
cat <<'EOF'
With a custom model the comparison is on by default: the whole screening again with MACE-MP-0 medium (its own clean
slab, references and relaxations). The table compares dE and E_form -- never E_ads, which carries each model's E0s.
EOF
ads case5 $ROUND3 --relax-depth 4 --bulk-reference "$DIR/si_bulk.fdf" --no-committee
section "$(report case5)" "FOUNDATION MODEL COMPARISON"
pause


banner "Case 6 -- The reference set of README section 6 (stored)"
cat <<'EOF'
make_reference.sh ran every case on a GPU: the main run with the committee, relaxation depths 0/1/4/7 Ang, four starting
heights, binding curves at fixed height, a thicker slab and a 3x3 cell, round 2 and the three foundation models.
results_figures.py reads them back (numpy + matplotlib, seconds) and prints the summary:
EOF
python3 results_figures.py | sed 's/^/  /'
echo
echo "  The committee on the main run's minima (reference report):"
section "$REF/main/stb_mladsorb_report.txt" "COMMITTEE" | head -12
pause


banner "Case 7 -- The SIESTA check (stored)"
cat <<'EOF'
Each minimum of the main run and the clean slab, relaxed with SIESTA at the training level of theory (PBE, DZP, 300 Ry,
k 4x4x1 = the training k x a), from the model's own geometry, the bottom bilayer fixed. siesta_compare.py reads the
relaxations (reference_results/siesta/): step 0 is SIESTA on the model's geometry, the last step SIESTA's own minimum.
EOF
python3 siesta_compare.py --no-model 2>/dev/null | sed 's/^/  /'
pause


banner "Case 8 -- Guard rails (seconds each)"
g() { echo "  \$ stb-mladsorb $*"; (cd "$OUT" && stb-mladsorb --no-intro --device "$DEVICE" "$@" < /dev/null 2>&1 \
        | strip_colors | grep -E 'error|FAIL|not a recognized' | head -2 | sed 's/^/      /') || true; }
g --slab "$DIR/nope.fdf" --adsorbate Si
g --slab "$DIR/si111_2x2.fdf" --adsorbate Xx
g --slab "$DIR/si111_2x2.fdf" --adsorbate Si --custom-model "$R3/mlff_model_compiled.model"
g --slab "$DIR/si111_2x2.fdf" --adsorbate Si --relax-depth -1
g --slab "$DIR/si111_2x2.fdf" --adsorbate Si --bulk-reference "$DIR/nope.fdf"
g --slab "$DIR/si_bulk.fdf" --adsorbate Si
cat <<'EOF'
  -> Every one stops before any relaxation: a missing slab, an unknown adsorbate, a custom model without its card, a
     negative depth, a missing bulk reference, and a "slab" with no vacuum.
EOF
pause


banner "Case 9 -- The interactive menu (stb-suite, code 7.6)"
cat <<'EOF'
The menu builds the same command line from its questions; here they are answered through a pipe: the 2x2 slab, Si,
the round-3 model and its card, no foundation comparison, the device, ontop sites only, default height, the top 4 Ang
relaxed, the bulk reference, all minima, no orientations, output folder, then data and report.
EOF
mkdir -p "$OUT/case10"
printf '7.6\n%s\nSi\n%s\n\nn\n%s\nontop\n\n4\n%s\n\n\nmenu_out\nn\ny\nn\nn\ny\n\n0\n' \
    "$DIR/si111_2x2.fdf" "$R3/mlff_model_compiled.model" "$DEVICE" "$DIR/si_bulk.fdf" \
    | (cd "$OUT/case10" && stb-suite > log.txt 2>&1) || true
R="$OUT/case10/menu_out/stb_mladsorb_report.txt"
if [ -f "$R" ]; then
    section "$R" "RUN METADATA"
    section "$R" "DISTINCT MINIMA" | head -6
    ls "$OUT/case10/menu_out" | sed 's/^/    /'
else
    echo "  The menu run did not produce a report -- see $OUT/case10/log.txt"
fi
pause


banner "Done"
cat <<EOF
Everything is in $OUT/. What to take from the tour:
  - A screening reports where each relaxation ENDED; the starting site is not a result. Read the sites by what is
    below them, not by pymatgen's names.
  - Relax the top of the slab (--relax-depth): a surface cut from the bulk is far from any model's minimum.
  - Measure energies from a reference the model knows (--bulk-reference): an isolated atom is outside a fine-tuned
    model, and its E_ads carries that error in full.
  - A model trained on the bulk has never seen a surface: the committee spread says so, and SIESTA says by how much
    (README section 7). Use the ranking to choose which sites to compute with SIESTA (stb-adsorb, 5.8).
EOF
