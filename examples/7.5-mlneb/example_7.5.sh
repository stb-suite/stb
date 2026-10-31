#!/bin/bash
# Guided example: 7.5 -- ML NEB (stb-mlneb; code 7.5 in the stb-suite menu).
#
# Not an automated test (see test/7-mlsimulations/5-mlneb/test.sh for that) -- a commented walk-through: it runs real
# commands, one case at a time, into its own output/<case>/ folder, and shows the piece of output that proves what just
# happened. Pauses between cases so you can read before moving on (`yes '' | bash example_7.5.sh` runs it straight
# through). Safe to re-run any time -- it always starts by wiping its own output/.
#
# The system is the one of README section 6: a vacancy hop in a 63-atom silicon supercell, with the fine-tuned models
# of the 7.1 silicon run (round 2: ../7.2-mlmd/model/; round 3: model_round3/). Three cases run a real NEB of that
# cell: about 1 minute each on a GPU (DEVICE=cuda), 9 minutes each on a CPU (the default) -- about 30 minutes in all on
# a CPU. Nothing here runs SIESTA: the SIESTA single points of the real loop are in reference_results/, and the script
# puts them into the folders stb-mlneb writes after checking that the geometry is the same. The full reference set of
# section 6 (foundation models, number of images, ...) is read from reference_results/ (Case 9; make_reference.sh
# regenerates it in about 8 minutes on a GPU).
#
# Needs the 'ml' extra: pip install stb_suite[ml]

set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

export MPLBACKEND=Agg           # headless-safe; the PNG figures are written either way
DEVICE="${DEVICE:-cpu}"
R2="$DIR/../7.2-mlmd/model"
R3="$DIR/model_round3"
ROUND2="--custom-model $R2/mlff_model_compiled.model --custom-model-card $R2/mlff_model_card.json"
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

# neb <case folder> <flags...>: runs stb-mlneb in output/<case>/ on the two endpoints, report saved, log kept
neb() {
    local case=$1; shift
    mkdir -p "$OUT/$case"
    local t0=$SECONDS
    (cd "$OUT/$case" && stb-mlneb --initial "$DIR/vac_initial.fdf" --final "$DIR/vac_final.fdf" --device "$DEVICE" \
        -o . --save-report --save-data --no-intro "$@" < /dev/null > log.txt 2>&1)
    echo "  (stb-mlneb took $((SECONDS - t0)) s on $DEVICE)"
}
report() { echo "$OUT/$1/stb_mlneb_report.txt"; }

# label <check folder> <finished folder>: copies SIESTA's output of a finished run into a folder stb-mlneb wrote,
# only when both were written for the same geometry (the stamp stb_mlneb_image.xyz, within 0.02 Ang)
label() {
    python3 - "$1" "$2" <<'EOF'
import os, shutil, sys
import numpy as np
from ase.io import read
new, old = sys.argv[1], sys.argv[2]
for name in sorted(os.listdir(new)):
    src = os.path.join(old, name)
    a = read(os.path.join(new, name, "stb_mlneb_image.xyz"))
    b = read(os.path.join(src, "stb_mlneb_image.xyz"))
    d = np.abs(a.get_positions() - b.get_positions()).max()
    if d > 0.02:
        print(f"    {name}: geometry differs by {d:.3f} Ang -- NOT labelled"); continue
    for f in ("siesta.out", "siesta.FA"):
        shutil.copy(os.path.join(src, f), os.path.join(new, name, f))
    print(f"    {name}: same geometry as the finished run (max difference {d:.1e} Ang) -> siesta.out + siesta.FA copied")
EOF
}


banner "Welcome: a migration barrier from a machine-learned potential"
cat <<EOF
stb-mlneb finds the minimum-energy path between two structures and its barrier with the climbing-image nudged elastic
band method, every energy and force from a MACE potential. With a fine-tuned model it also checks the band against
the model card, against the model's committee, and -- in two passes -- against SIESTA single points on the same
geometries; and it can feed the band back to the training set of the MLFF workflow (7.1).

The system: a vacancy in a 63-atom silicon supercell hops to a neighbouring site. The tour:
   1  the two endpoints                          6  the real loop: committee and retraining (stored)
   2  a workflow to export to (7.1 Stage 1)     7  round 3: the NEB again, checked     [NEB]
   3  round 2: NEB + SIESTA folders     [NEB]    8  guard rails
   4  SIESTA's results go into the folders       9  the reference set of README section 6 (stored)
   5  round 2 again: check + export     [NEB]   10  the interactive menu (GPU only)
Device: $DEVICE -- each [NEB] case takes about $([ "$DEVICE" = cuda ] && echo "1 minute" || echo "9 minutes").
Theory: README section 1. Numbers: README section 6.
EOF
pause


banner "Case 1 -- The endpoints: a vacancy before and after one hop"
cat <<'EOF'
make_endpoints.py builds both from structure.fdf (the 8-atom diamond cell at SIESTA's lattice constant, 5.48 Ang):
a 2x2x2 supercell (64 Si), the atom at the origin removed (initial), then its neighbour at (1/8,1/8,1/8) moved into
the hole (final). Both files list the same 63 atoms in the same order -- NEB interpolates atom j to atom j.
EOF
python3 make_endpoints.py "$OUT/case1"
echo
echo "  Same as the files shipped in this folder (used by every run below)?"
for f in vac_initial.fdf vac_final.fdf; do
    if diff <(grep -v '^#' "$OUT/case1/$f") <(grep -v '^#' "$f") > /dev/null; then echo "    $f: identical"
    else echo "    $f: DIFFERENT"; fi
done
echo
echo "  The only line that differs between the two endpoints (atom 1, the one that hops):"
diff vac_initial.fdf vac_final.fdf | grep '^[<>]' | grep -v '#' | sed 's/^/    /'
pause


banner "Case 2 -- A workflow to send the band to (7.1 Stage 1, seconds, no SIESTA)"
cat <<'EOF'
--siesta-check and --to-mlff take the level of theory from the MLFF workflow that trained the model: its calc.fdf,
its pseudopotentials, its dispersion and its k-grid. The real one (README section 6) is 2 GB; this case makes a small
one with the same reference (structure.fdf + calc.fdf: PBE, DZP, 300 Ry, k 6x6x6) and two rattled configs, so that
Case 5 has something to export to.
EOF
mkdir -p "$OUT/case2"
cp structure.fdf calc.fdf "$OUT/case2/"
(cd "$OUT/case2" && stb-mlff --no-intro -s structure.fdf -c calc.fdf --n-configs 2 --pseudo-dir dojo --no-plots \
    > log.txt 2>&1)
WF="$OUT/case2/mlff_runs/mlff_workflow.json"
echo "  $WF:"
python3 -c "
import json; d = json.load(open('$WF'))
print('    reference :', d['reference']['file'], '+', d['reference']['calc'], '| pseudopotentials', d['reference']['pseudo_dir'])
print('    configs   :', len(d['configs']), '(round 1, not labelled -- not needed here) | rounds trained:', len(d['rounds']))"
pause


banner "Case 3 -- Round 2: the NEB, and SIESTA folders for every image   [NEB]"
cat <<'EOF'
The model of the 7.2-7.4 examples: round 2 of the 7.1 silicon run, 57 SIESTA configs of the PERFECT crystal (rattled,
strained, MD), test force error 28.7 meV/Ang. --siesta-check --siesta-check-images all writes a single-point folder for
each of the 7 images, at the workflow's level of theory (--mlff-workflow). --no-committee: the committee members are
not shipped (Case 6 shows what they said); --skip-foundation-comparison: Case 9 has it.
EOF
neb case3 $ROUND2 --mlff-workflow "$WF" --no-committee --skip-foundation-comparison \
    --siesta-check --siesta-check-images all
R="$(report case3)"
section "$R" "ENDPOINTS"
section "$R" "NEB BAND"
section "$R" "BARRIER"
pause
section "$R" "SADDLE POINT"
section "$R" "TRAINED RANGE"
cat <<'EOF'
  -> The saddle: the hopping atom halfway, 2.78 Ang from the three atoms on each side (a 'split vacancy'). The band
     is INSIDE every trained range of the card -- volume, bond length, force, energy. Nothing looks wrong. 0.72 eV.
EOF
pause
section "$R" "SIESTA CHECK"
echo "  One of the folders, $OUT/case3/siesta_check/image_03/:"
ls "$OUT/case3/siesta_check/image_03" | sed 's/^/    /'
echo "  Its config_extra.fdf (%included on top of the workflow's calc.fdf -- first occurrence wins in SIESTA):"
sed 's/^/    /' "$OUT/case3/siesta_check/image_03/config_extra.fdf"
cat <<'EOF'
  -> k-grid 3x3x3: the training cell (5.48 Ang) used 6x6x6; the 2x2x2 supercell (10.96 Ang) keeps k x a = 32.9 Ang.
     stb_mlneb_image.xyz is the stamp: the geometry this folder was written for.
EOF
pause


banner "Case 4 -- SIESTA's results go into the folders"
cat <<'EOF'
Here you would run SIESTA in each folder (about 3 minutes per image on 6 cores, 63 atoms):
    for d in output/case3/siesta_check/image_*/; do (cd $d && mpirun -np 6 siesta < calc.fdf > siesta.out); done
They were run once, in the real loop (README section 6.4), and are kept in reference_results/loop_round2/. The band of
Case 3 is the same band (same model, same endpoints, same settings), so their output belongs in these folders -- the
script checks the stamps first, as stb-mlneb itself will:
EOF
label "$OUT/case3/siesta_check" "$REF/loop_round2/siesta_check"
echo
echo "  SIESTA's total energies (eV):"
for d in "$OUT"/case3/siesta_check/image_*; do
    printf "    %s  %s\n" "$(basename "$d")" "$(grep 'siesta: FreeEng' "$d/siesta.out" | tail -1 | awk '{print $4}')"
done
pause


banner "Case 5 -- Round 2 again: the check, and the export to training   [NEB]"
cat <<'EOF'
The SAME command as Case 3 (same -o), plus --to-mlff --to-mlff-rattle 1. The band is recomputed (the same, to the
digits); the finished folders whose stamp matches it are read; the images go to the workflow as training configs of
its next round -- the 7 checked ones copied WITH their SIESTA output, plus one rattled copy of each (0.05 Ang).
EOF
neb case3 $ROUND2 --mlff-workflow "$WF" --no-committee --skip-foundation-comparison \
    --siesta-check --siesta-check-images all --to-mlff --to-mlff-rattle 1
R="$(report case3)"
section "$R" "SIESTA CHECK"
cat <<'EOF'
  -> On its own path, SIESTA's barrier is 0.25 eV: the model is +186 % off, 0.47 eV at the saddle. The error grows
     toward the saddle (images 1, 2, 3: +0.04, +0.21, +0.47 eV) -- the farther from the crystal, the worse.
EOF
pause
section "$R" "EXPORT TO THE MLFF WORKFLOW"
python3 -c "
import json; d = json.load(open('$WF'))
neb = [c for c in d['configs'] if c.get('stage') == 'neb']
print(f'  The workflow now holds {len(neb)} NEB configs, e.g.:')
for c in neb[:2] + neb[-1:]:
    print('    ', {k: c[k] for k in ('folder', 'set', 'kind', 'round')}, '| E_rel model', c['neb']['E_rel_model_eV'], 'eV')
print('  Last history entry:', d['history'][-1]['stage'], d['history'][-1]['tool'])"
section "$R" "WARNINGS"
pause


banner "Case 6 -- The real loop: what the committee said, and round 3 (stored reports)"
cat <<'EOF'
In the real run (README section 6), the same NEB found round 2's committee next to the workflow. The three members'
barriers span 0.63-0.72 eV: they disagree by 0.09 eV, a fifth of their common error against SIESTA (0.47 eV). Seeds
fine-tuned from one foundation model are wrong together.
EOF
section "$REF/loop_round2/stb_mlneb_report.txt" "COMMITTEE UNCERTAINTY"
pause
cat <<'EOF'
Then: SIESTA on the 7 rattled copies (21 min on 6 cores), stb-mlffAnalysis --init-from previous (round 3, three
members, ~10 min each on the GPU), stb-mlffValidate. Round 3 is also the new best model on the perfect crystal:
EOF
strip_colors < "$REF/mlff_round3/stage3_r3.txt" | sed -n '/^\[3\] TEST ERRORS/,/^Best round/p' | head -9 | sed 's/^/  /'
strip_colors < "$REF/mlff_round3/stage3_r3.txt" | grep -E '^\[STOP\]' | sed 's/^/  /'
pause


banner "Case 7 -- Round 3: the NEB again, checked   [NEB]"
cat <<'EOF'
The round-3 model (model_round3/, its card says 'neb 7, neb-rattle 7' among the training kinds). The ends and the
saddle of its band were checked with SIESTA in the real loop; those folders are put in place BEFORE the run, so one
pass reads them -- stb-mlneb compares each stamp with its own band and would call a folder 'stale' if it differed.
EOF
mkdir -p "$OUT/case7"
cp -r "$REF/loop_round3/siesta_check" "$OUT/case7/"
neb case7 $ROUND3 --mlff-workflow "$WF" --no-committee --skip-foundation-comparison --siesta-check
R="$(report case7)"
section "$R" "MODEL" | grep -E 'Training set|Training kinds|test set'
section "$R" "BARRIER"
section "$R" "SIESTA CHECK"
cat <<'EOF'
  -> 0.31 eV against SIESTA's 0.25 on the same path (+23 %), from +186 %. SIESTA's force on the model's relaxed
     endpoints dropped from 0.23 to 0.04 eV/Ang: the model's vacancy is now close to SIESTA's minimum. The warning
     says what is left: one more round with these images, or SIESTA's number.
EOF
section "$R" "WARNINGS"
pause


banner "Case 8 -- Guard rails (seconds each)"
g() { echo "  \$ stb-mlneb $*"; (cd "$OUT" && stb-mlneb --no-intro --device "$DEVICE" "$@" < /dev/null 2>&1 \
        | strip_colors | grep -E 'error|FAIL' | head -2 | sed 's/^/      /') || true; }
mkdir -p "$OUT/case8"
sed 's/^\(  0.12500000   0.12500000   0.12500000\)   1/\1   2/' vac_initial.fdf \
    | sed 's/^NumberOfSpecies 1/NumberOfSpecies 2/; s/^ 1 14 Si/ 1 14 Si\n 2 32 Ge/' > "$OUT/case8/ge_initial.fdf"
g --initial "$DIR/vac_initial.fdf" --final "$DIR/vac_final.fdf" --n-images 2
g --initial "$DIR/vac_initial.fdf" --final "$DIR/vac_final.fdf" --custom-model "$R3/mlff_model_compiled.model"
g --initial "$DIR/vac_initial.fdf" --final "$DIR/vac_final.fdf" --to-mlff
g --initial "$DIR/vac_initial.fdf" --final "$DIR/vac_final.fdf" $ROUND3 --siesta-check
g --initial "$OUT/case8/ge_initial.fdf" --final "$DIR/vac_final.fdf"
cat <<'EOF'
  -> Every one stops before the NEB: an image count without an interior image, a custom model without its card,
     an export without a workflow, a SIESTA check without a level of theory, and endpoints that are not the same atoms.
EOF
pause


banner "Case 9 -- The reference set of README section 6 (stored)"
cat <<'EOF'
make_reference.sh ran every case on a GPU: three foundation models, rounds 2 and 3 with their SIESTA checks, the number
of images, IDPP and the endpoint relaxation. results_figures.py reads them back (numpy + ase + matplotlib, seconds):
EOF
python3 results_figures.py | sed 's/^/  /'
pause


banner "Case 10 -- The interactive menu (stb-suite, code 7.5)"
if [ "$DEVICE" != "cuda" ]; then
    cat <<'EOF'
  Skipped on the CPU (it runs one more 9-minute NEB). With DEVICE=cuda it answers the menu's questions through a pipe:
  the round-3 model and its card, cuda, 7 images, no foundation comparison, the workflow of Case 2, --siesta-check of
  the ends and the saddle, no export, output folder, then the optional outputs. README section 5 lists the questions.
EOF
else
    mkdir -p "$OUT/case10"
    printf '7.5\n%s\n%s\n%s\n\nmedium\ncuda\n7\nn\nn\n%s\ny\nn\nn\nmenu_out\nn\ny\ny\nn\ny\n\n0\n' \
        "$DIR/vac_initial.fdf" "$DIR/vac_final.fdf" "$R3/mlff_model_compiled.model" "$WF" \
        | (cd "$OUT/case10" && stb-suite > log.txt 2>&1) || true
    R="$OUT/case10/menu_out/stb_mlneb_report.txt"
    if [ -f "$R" ]; then
        section "$R" "RUN METADATA"
        section "$R" "SIESTA CHECK" | head -6
        ls "$OUT/case10/menu_out" | sed 's/^/    /'
    else
        echo "  The menu run did not produce a report -- see $OUT/case10/log.txt"
    fi
fi
pause


banner "Done"
cat <<EOF
Everything is in $OUT/. What to take from the tour:
  - A model validated on the perfect crystal can be 3x wrong on a barrier while its trained range, its committee and
    its test-set error all look fine. A transition state is not in a crystal's training set.
  - SIESTA single points on the model's own band (--siesta-check) expose it in a few minutes per image, and the same
    images, fed back to training (--to-mlff), fix most of it in one round.
  - For a number to publish: SIESTA's NEB (stb-neb, 5.9), starting from this band (--save-images), in a converged cell.
EOF
