#!/bin/bash
# Guided example: 7.4 -- ML Elastic Constants (stb-mlelastic; code 7.4 in the stb-suite menu).
#
# Not an automated test (see test/7-mlsimulations/4-mlelastic/test.sh for that) -- a commented walk-through: it runs
# real commands, one case at a time, into its own output/<case>/ folder, and shows the piece of output that proves
# what just happened. Pauses between sections so you can read before moving on. Safe to re-run any time -- it
# always starts by wiping its own output/.
#
# Every calculation is real and runs on the CPU in about 6 minutes in total (2- and 8-atom cells), with the
# fine-tuned model of the 7.2 example (../7.2-mlmd/model/, trained by the 7.1 workflow on SIESTA PBE/DZP data of
# bulk silicon) and the MACE-MP-0 medium foundation model (downloaded once on first use). Set DEVICE=cuda to run on
# a GPU. The full set of README section 5 is NOT rerun here: Case 12 reads it from reference_results/
# (make_reference.sh regenerates it in about 15 minutes on a GPU), and Case 13 reads the SIESTA runs at the model's
# own level of theory from reference_results/siesta_elastic/ (siesta_elastic.sh reruns them, ~6 min on 6 cores).
#
# Needs the 'ml' extra: pip install stb_suite[ml]

set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

export MPLBACKEND=Agg           # headless-safe; the PNG figures are written either way
DEVICE="${DEVICE:-cpu}"
MODEL="$DIR/../7.2-mlmd/model/mlff_model_compiled.model"
CARD="$DIR/../7.2-mlmd/model/mlff_model_card.json"
FT="--custom-model $MODEL --custom-model-card $CARD"

OUT="$DIR/output"
rm -rf "$OUT"
mkdir -p "$OUT"

pause() {
    echo
    read -p "  [Press Enter to continue] " -r
    echo
}

banner() {
    echo "=================================================================="
    echo " $1"
    echo "=================================================================="
}

strip_colors() { sed 's/\x1b\[[0-9;]*m//g'; }

# a report section by its title: from "[N] TITLE..." up to the next "[M] " header. By title, not number: the
# sections present depend on the mode and the dimensionality, so their numbers move. A header starts with a capital
# letter: the '[100] | ...' rows of the sound-velocity table are not headers.
section() {   # $1 = file, $2 = start of the title
    strip_colors < "$1" | awk -v t="$2" '/^\[[0-9]+\] [A-Z]/ {if (on) exit; if (index($0, "] " t) > 0) on = 1}
        on {print}'
}

# mlel <case folder> <structure> <flags...>: runs stb-mlelastic with output into that folder, report saved, log kept
mlel() {
    local case=$1 file=$2; shift 2
    mkdir -p "$OUT/$case"
    (cd "$OUT/$case" && stb-mlelastic --file "$file" --device "$DEVICE" -o . --save-report --no-intro "$@" \
        < /dev/null > log.txt 2>&1)
}
report() { echo "$OUT/$1/stb_mlelastic_report.txt"; }


banner "Welcome: elastic constants from a machine-learned potential"
cat <<'EOF'
stb-mlelastic measures the stiffness tensor C of a crystal by the stress-strain method: it relaxes the cell, strains
it along each of the six independent directions (xx, yy, zz, yz, xz, xy) at a few amplitudes up to +-1 %, relaxes
the atoms inside each strained cell, and fits stress against strain. The slope of each line is one column of C.
The potential is the MACE-MP-0 foundation model, or a model fine-tuned on YOUR SIESTA calculations (7.1), which
always comes with its model card.

From C the report derives: the constants allowed by the crystal's symmetry, mechanical stability, bulk/shear/Young
moduli and Poisson ratio, the effect of the internal relaxation, sound velocities and the Debye temperature, and
the directional Young's modulus -- and, with a fine-tuned model, the same with the foundation model.

This tour (CPU): a first run and its report (1-4), the force threshold (5), the strain amplitude (6), a hexagonal
crystal (7), a 2D sheet (8), a structure that is not at a minimum (9), the guard rails (10), the interactive menu
(11), and the full reference set of the README (12). Theory: README section 1. Numbers: README section 5.
EOF
pause


banner "Case 1: silicon, the fine-tuned model against the foundation model"
cat <<'EOF'
si8.fdf is the 8-atom conventional diamond cell at SIESTA's relaxed lattice constant (5.48024 Ang). Each model
first relaxes the cell to ITS OWN stress-free minimum: the fine-tuned model hardly moves it (it reproduces SIESTA),
the foundation model shrinks it by 1.35 % (its silicon is smaller). Then 6 directions x 5 strains each.
EOF
echo "\$ stb-mlelastic --file si8.fdf $FT --steps 5 --save-data --save-gnuplot --save-report"
mlel case1_si "$DIR/si8.fdf" $FT --steps 5 --save-data --save-gnuplot
R1=$(report case1_si)
section "$R1" "REFERENCE STATE"
pause
cat <<'EOF'
Stress-strain fits: one line per direction. R^2 ~ 0.99999, intercepts ~0 (the reference is stress-free). The last
columns are the internal relaxation: under the normal strains xx/yy/zz the atoms stay where the strain puts them
(diamond's symmetry pins them), under the shears they move by ~0.008 Ang per 1 % -- that is what lowers C44.
The trained-range table: every strained cell is inside the volumes and pressures the model was trained on.
EOF
section "$R1" "STRESS-STRAIN FITS" | sed -n '1,13p'
section "$R1" "TRAINED RANGE"
pause


banner "Case 2: the tensor, its symmetry, and stability"
cat <<'EOF'
The fitted 6x6 tensor is averaged over the 48 rotations of the point group m-3m (spglib on the relaxed cell): the
raw fit already HAS that symmetry (deviation 0.00 %), and only 3 constants are independent: C11, C12, C44. Born's
criteria for a cubic crystal: C11 - C12 > 0, C11 + 2 C12 > 0, C44 > 0 -- the general statement is that all
eigenvalues of C are positive.
EOF
section "$R1" "STIFFNESS TENSOR" | sed -n '1,20p'
section "$R1" "MECHANICAL STABILITY" | sed -n '1,9p'
pause


banner "Case 3: what the constants mean -- moduli, internal relaxation, sound, anisotropy"
cat <<'EOF'
Polycrystal moduli (Voigt/Reuss/Hill) -- experiment for silicon: B 97.8, G ~66 GPa (PBE is softer). Clamped vs
relaxed ions: C44 drops by a third when the atoms may relax inside the cell -- the measured constant is the relaxed
one. Christoffel velocities: the same numbers as the slopes of the acoustic phonons in 7.3. Directional Young's
modulus: silicon is stiffest along <111>, softest along <100>.
EOF
section "$R1" "ELASTIC MODULI" | sed -n '1,13p'
section "$R1" "INTERNAL RELAXATION" | sed -n '1,9p'
section "$R1" "SOUND VELOCITIES" | sed -n '1,12p'
section "$R1" "DIRECTIONAL" | sed -n '1,8p'
pause
cat <<'EOF'
The foundation model: same functional (PBE), same crystal -- and C44 52 % too low. Each model is strained about its
own minimum, so the comparison is fair.
EOF
section "$R1" "FINE-TUNED vs FOUNDATION"
pause


banner "Case 4: the outputs"
cat <<'EOF'
PNG figures are always written: stress_strain.png (every fit), directional_moduli.png (polar E and compressibility),
anisotropy_surface.png (3D E surface). --save-data writes the numbers, --save-gnuplot .dat + .gplot (+ PDF when
gnuplot is installed), --save-report the report, --show opens the figures. The last section lists every warning.
EOF
ls -1 "$OUT/case1_si" "$OUT/case1_si/plot" | sed 's/^/  /'
section "$R1" "SUMMARY" | sed -n '1,4p'
section "$R1" "WARNINGS"
pause


banner "Case 5: the force threshold of the internal relaxation"
cat <<'EOF'
Under shear the atoms move by only ~0.008 Ang per % strain. A relaxation stopped at a residual force f leaves them
~f/k off (k ~ 13 eV/Ang^2): f = 0.02 eV/Ang is 0.0015 Ang -- a fifth of the effect being measured. The default is
0.001; here 0.02 against it.
EOF
for f in 0.02 0.001; do
    mlel case5_fmax_$f "$DIR/si8.fdf" $FT --steps 5 --fmax $f --skip-foundation-comparison --no-clamped
    printf "  --fmax %-6s " "$f"
    section "$(report case5_fmax_$f)" "STRESS-STRAIN FITS" | grep -E "^yz " | awk -F'|' '{printf "C44 %s  R^2 %s  relax steps %s\n", $3, $4, $9}'
done
pause


banner "Case 6: the strain amplitude -- --check-convergence"
cat <<'EOF'
A finite strain also measures the next (third-order) terms of the energy. The same fit at 0.5, 1 and 2 %: the
constants drift smoothly -- C44 most -- and the report says by how much.
EOF
echo "\$ stb-mlelastic --file si8.fdf $FT --check-convergence --convergence-strains 0.5 1 2 --steps 5"
mlel case6_conv "$DIR/si8.fdf" $FT --check-convergence --convergence-strains 0.5 1 2 --steps 5
section "$(report case6_conv)" "STRAIN-AMPLITUDE"
pause


banner "Case 7: a hexagonal crystal -- hcp magnesium with the foundation model"
cat <<'EOF'
The 7.2 model knows only silicon; hcp Mg runs with MACE-MP-0 medium (no card). The point group 6/mmm leaves FIVE
independent constants: C11 = C22, C12, C13 = C23, C33, C44 = C55; C66 is not free -- it equals (C11 - C12)/2, a
linear relation the report recognizes and marks. Experiment: C11 59.4, C12 25.6, C13 21.4, C33 61.6, C44 16.4 GPa.
EOF
mlel case7_mg "$DIR/mg_hcp.fdf" --model medium --steps 5
section "$(report case7_mg)" "STIFFNESS TENSOR" | sed -n "1,22p"
section "$(report case7_mg)" "MECHANICAL STABILITY" | sed -n '1,9p'
pause


banner "Case 8: a 2D sheet -- graphene"
cat <<'EOF'
Vacuum along z: the directions zz, yz, xz would only strain empty space and are skipped. Stresses are multiplied by
the cell height, so the constants are per area, in N/m. Experiment (nanoindentation): E2D = 340 +- 50 N/m.
EOF
mlel case8_graphene "$DIR/graphene.fdf" --model medium --steps 5 --max 0.5
strip_colors < "$(report case8_graphene)" | grep -E "^Skipped|^Constants" | sed 's/^/  /'
section "$(report case8_graphene)" "IN-PLANE MODULI"
pause


banner "Case 9: a structure that is NOT at a minimum -- the pre-relax"
cat <<'EOF'
si8_rattled.fdf: every atom moved by ~0.03 Ang. Without the pre-relax the cell is strained around a point that is
not a minimum: forces of 0.85 eV/Ang, a residual stress, no symmetry -- 21 constants instead of 3. With it (the
default), the atoms return to their sites, spglib finds Fd-3m again, and the constants are the crystal's.
EOF
mlel case9_norelax "$DIR/si8_rattled.fdf" $FT --steps 5 --skip-foundation-comparison --no-clamped --no-pre-relax
mlel case9_relax "$DIR/si8_rattled.fdf" $FT --steps 5 --skip-foundation-comparison --no-clamped
for c in case9_norelax case9_relax; do
    echo "--- $c"
    section "$(report $c)" "REFERENCE STATE" | grep -E "^Pre-relax|^Stress"
    section "$(report $c)" "STIFFNESS TENSOR" | grep -E "^Crystal symmetry|^Independent|^Projection"
done
pause


banner "Case 10: guard rails"
try() {   # try <label> <args...>: runs and shows the error line
    echo "--- $1"
    shift
    (cd "$OUT" && stb-mlelastic --no-intro "$@" < /dev/null 2>&1 | strip_colors | grep -E "error|ERROR|FAIL" | head -2) || true
}
try "--custom-model without its card" --file "$DIR/si8.fdf" --custom-model "$MODEL"
try "--steps 1 (a fit needs two strains)" --file "$DIR/si8.fdf" --steps 1
try "--convergence-strains with one value" --file "$DIR/si8.fdf" --check-convergence --convergence-strains 1
try "a negative --max" --file "$DIR/si8.fdf" --max -1
try "an input file that does not exist" --file "$DIR/nothere.fdf"
pause


banner "Case 11: the same tool from the interactive menu (stb-suite, code 7.4)"
cat <<'EOF'
The menu asks the same questions and builds the same command line. Answers piped in: structure, the custom model
(its card next to it), cpu, 1 %, 5 points, auto dimensionality, output folder menu_out, no convergence check, no
foundation comparison, keep the anisotropy surface, no display, no gnuplot, no data, save the report.
EOF
mkdir -p "$OUT/case11_menu"
cd "$OUT/case11_menu"
printf '7.4\n%s\n%s\n\ncpu\n1.0\n5\n\nmenu_out\nn\nn\nn\nn\nn\nn\ny\n\n0\n' "$DIR/si8.fdf" "$MODEL" | \
    stb-suite > menu_log.txt 2>&1 || true
ls -1 menu_out 2>/dev/null | sed 's/^/  /'
strip_colors < menu_out/stb_mlelastic_report.txt 2>/dev/null | grep -E "^Device|^SHA-256|^Status|^Results" | sed 's/^/  /'
cd "$DIR"
pause


banner "Case 12: the full reference set of README section 5"
cat <<'EOF'
reference_results/ holds the GPU runs of make_reference.sh: the 7-point reference, the amplitude scan 0.25-3 %, four
force thresholds, three foundation models, a 216-atom supercell, the rattled cell both ways, hcp Mg and graphene.
results_figures.py turns them into the tables of the README (numpy + matplotlib only).
EOF
python3 "$DIR/results_figures.py" --no-plots | head -60
echo "   ... (python3 results_figures.py prints every table and redraws figures/)"
pause


banner "Case 13: SIESTA at the model's own level of theory"
cat <<'EOF'
The same constants computed with SIESTA and the training calc.fdf (PBE, DZP, 300 Ry, k 6x6x6), through the 5.2 workflow:
clamped and relaxed ions, plus an isotropic scan for B. siesta_compare.py reads the extracted runs and sets SIESTA,
the model and experiment side by side -- it separates the model's error from the error of PBE/DZP itself.
EOF
python3 "$DIR/siesta_compare.py" --no-figure | sed -n '/^BULK/,/^DERIVED/p'
echo "   ... (python3 siesta_compare.py also prints the derived ratios and redraws figures/siesta_vs_model.png)"
pause


banner "Done"
cat <<EOF
Every case wrote into $OUT/<case>/ (report, PNG figures, and whatever was asked for).
README.md explains the theory (section 1), the report (section 3) and every number of the examples (section 5).
EOF
