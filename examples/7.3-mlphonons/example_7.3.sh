#!/bin/bash
# Guided example: 7.3 -- ML Phonons (stb-mlphonons; code 7.3 in the stb-suite menu).
#
# Not an automated test (see test/7-mlsimulations/3-mlphonons/test.sh for that) -- a commented walk-through: it runs
# real commands, one case at a time, into its own output/<case>/ folder, and shows the piece of output that proves
# what just happened. Pauses between sections so you can read before moving on. Safe to re-run any time -- it
# always starts by wiping its own output/.
#
# Every calculation is real and runs on the CPU in about 3 minutes in total (bulk silicon, 2x2x2 supercells of 64
# atoms), with the fine-tuned model of the 7.2 example (../7.2-mlmd/model/, trained by the 7.1 workflow on SIESTA
# PBE/DZP data of bulk silicon) and, in Cases 1 and 9, the MACE-MP-0 medium foundation model (downloaded once on
# first use). Set DEVICE=cuda to run on a GPU. The full set of README section 5 (4x4x4 supercells, three foundation
# models, the SIESTA comparison) is NOT rerun here: Case 12 reads it from reference_results/ (make_reference.sh
# regenerates it in about 10 minutes on a CPU).
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
# optional sections (MODE FREEZE, MODE ANIMATION, FINE-TUNED vs FOUNDATION) shift the numbers after them.
section() {   # $1 = file, $2 = start of the title
    strip_colors < "$1" | awk -v t="$2" '/^\[[0-9]+\] / {if (on) exit; if (index($0, "] " t) > 0) on = 1}
        on {print}'
}

# mlph <case folder> <structure> <flags...>: runs stb-mlphonons with output into that folder, report saved, log kept
mlph() {
    local case=$1 file=$2; shift 2
    mkdir -p "$OUT/$case"
    (cd "$OUT/$case" && stb-mlphonons --file "$file" --device "$DEVICE" -o . --save-report --no-intro "$@" \
        < /dev/null > log.txt 2>&1)
}
report() { echo "$OUT/$1/stb_mlphonons_report.txt"; }


banner "Welcome: phonons from a machine-learned potential"
cat <<'EOF'
stb-mlphonons computes the phonons of a crystal by finite displacements (phonopy): it moves one atom of a supercell
by 0.01 Ang, asks a MACE potential for the forces on every atom, builds the force constants and diagonalizes the
dynamical matrix -- band structure, density of states, sound velocities, heat capacity, entropy, free energy. The
potential is the MACE-MP-0 foundation model, or a model fine-tuned on YOUR SIESTA calculations (7.1), which always
comes with its model card. Two more modes: --qha (thermal expansion from a volume scan) and --check-convergence
(the same calculation in several supercells).

This tour (bulk silicon, CPU): a first calculation and its report (1), the optional outputs (2), density of states
and thermal properties (3), a vibration as a movie (4), the supercell (5), the displacement size (6), a structure
that is not at a minimum (7), the quasi-harmonic approximation (8), the foundation model alone (9), the guard rails
(10), the interactive menu (11), and the full reference set of the README (12). Theory: README section 1. Worked
example with the numbers: README section 5.
EOF
pause


banner "Case 1: a first calculation -- silicon, 2x2x2 supercell, fine-tuned model vs foundation model"
cat <<'EOF'
si8.fdf is the 8-atom conventional diamond cell at SIESTA's relaxed lattice constant (5.48024 Ang). The model comes
with its card: --custom-model AND --custom-model-card, always both. With a custom model the same calculation is
repeated with the raw MACE-MP-0 medium model, to show what fine-tuning changed.

Report section [3]: the diamond structure has so much symmetry that ONE displaced supercell (64 atoms) gives every
force constant -- 1 calculation instead of 384. The forces on the undisplaced supercell are subtracted, and the
acoustic sum rule (a rigid translation costs no energy) is imposed.
EOF
echo "\$ stb-mlphonons --file si8.fdf $FT --dim 2 2 2 --save-data --save-gnuplot --save-report"
mlph case1_first "$DIR/si8.fdf" $FT --dim 2 2 2 --save-data --save-gnuplot
R1=$(report case1_first)
section "$R1" "DISPLACEMENTS"
pause
cat <<'EOF'
Section [4]: the frequencies at the special points of the FCC Brillouin zone. Silicon, by inelastic neutron
scattering (Nilsson & Nelin, 1972): Gamma 15.53; X 4.49 / 12.32 / 13.90; L 3.43 / 11.35 / 12.60 / 14.68 THz.
The fine-tuned model is within ~2 % (PBE is a little soft); the foundation model's optical modes are ~30 % too low.
At Gamma: three acoustic modes at zero (T1u) and the triply degenerate optical mode (T2g, Raman-active).
EOF
section "$R1" "FREQUENCIES"
pause
cat <<'EOF'
Section [5]: stability (no imaginary mode on the band path or the 20x20x20 mesh) and the sound velocities -- the
slopes of the acoustic branches at small q. rho v^2 is an elastic constant: C11 = rho v_LA[100]^2, C44 = rho
v_TA[100]^2 (README section 1.7).
EOF
section "$R1" "BAND STRUCTURE" | grep -v "^\[Foundation" | sed -n '1,14p'
pause


banner "Case 2: the optional outputs -- data files and gnuplot"
cat <<'EOF'
PNG figures (bands, DOS, thermal) are always written. --save-data adds the numbers behind them; --save-gnuplot
writes .dat + .gplot into plot/ and runs gnuplot there (one PDF per figure) when it is installed; --save-report
writes the report; --show opens the figures on screen.
EOF
ls -1 "$OUT/case1_first" "$OUT/case1_first/plot" | sed 's/^/  /'
echo
echo "--- head of bands.dat (path distance, then one column per branch) ---"
head -3 "$OUT/case1_first/bands.dat"
pause


banner "Case 3: density of states and thermal properties"
cat <<'EOF'
Section [8]: the DOS on the 20x20x20 q-mesh -- it must integrate to 3N states per primitive cell (6 for silicon's
2 atoms) -- and the Debye temperature from its second moment (calorimetric value for silicon: 645 K).
Section [9]: every mode is a quantum harmonic oscillator. Zero-point energy, free energy, entropy and heat capacity
per primitive cell AND per atom. At 300 K the real crystal has Cp = 20.0 J/K/mol = 2.41 kB per atom and
S = 18.8 J/K/mol = 2.26 kB per atom.
EOF
section "$R1" "DENSITY OF STATES" | sed -n '1,9p'
section "$R1" "THERMAL PROPERTIES" | sed -n '1,12p'
pause
cat <<'EOF'
Section [10]: fine-tuned against foundation on the same mesh. The foundation model's softer modes give it a larger
entropy and heat capacity and a lower free energy at every temperature: a phonon error is a free-energy error.
EOF
section "$R1" "FINE-TUNED vs FOUNDATION"
pause


banner "Case 4: a vibration as a movie -- the optical mode at Gamma"
cat <<'EOF'
--animate-mode with no value animates the first optical mode (band 3; bands 0-2 are the acoustic modes, at Gamma a
rigid translation of the whole crystal). 30 frames of one period, amplitude 0.5 Ang (exaggerated, to be seen).
Open mode_3.xsf in VESTA, XCrySDen, OVITO or 'ase gui': the two atoms of the primitive cell move against each other
along the bond -- the Raman mode of silicon.
EOF
echo "\$ stb-mlphonons --file si8.fdf $FT --skip-foundation-comparison --animate-mode"
mlph case4_animate "$DIR/si8.fdf" $FT --skip-foundation-comparison --animate-mode
section "$(report case4_animate)" "MODE ANIMATION"
pause


banner "Case 5: how large must the supercell be? --check-convergence"
cat <<'EOF'
A supercell of N cells gives exact frequencies only at the q-points commensurate with it; elsewhere it interpolates
(README section 1.4). The same calculation in 1x1x1, 2x2x2 and 3x3x3 supercells, compared with the largest. Gamma
and X are commensurate with the 8-atom cubic cell itself, so they never change; L and K tell the story.
EOF
echo "\$ stb-mlphonons --file si8.fdf $FT --check-convergence --convergence-dims 1 1 1 2 2 2 3 3 3"
mlph case5_conv "$DIR/si8.fdf" $FT --check-convergence --convergence-dims 1 1 1 2 2 2 3 3 3
R5=$(report case5_conv)
section "$R5" "SUPERCELL" | sed -n '1,9p'
echo
echo "Frequencies at L and K per supercell (THz):"
section "$R5" "SUPERCELL" | awk '/^(L|K) +\|/ {n = 3} n > 0 {print; n--}'
section "$R5" "SUPERCELL" | grep -E "^Converged|^Last step"
pause


banner "Case 6: the displacement size -- 0.001, 0.01 and 0.1 Ang"
cat <<'EOF'
A force constant is a finite difference: too large a displacement measures anharmonicity, too small one measures
the noise of the forces (README section 1.5). MACE in float64 has almost no noise, so 0.001 Ang still works here;
at 0.1 Ang the cubic anharmonicity of silicon starts to show.
EOF
for d in 0.001 0.01 0.1; do
    mlph case6_dist_$d "$DIR/si8.fdf" $FT --skip-foundation-comparison --distance $d
    printf "  %-6s Ang  " "$d"
    section "$(report case6_dist_$d)" "FREQUENCIES" | grep -E "^(Γ|X) " | awk -F'|' '{printf "%s %s   ", $1, $4}'
    echo
done
pause


banner "Case 7: a structure that is NOT at a minimum -- the pre-relax"
cat <<'EOF'
si8_rattled.fdf is si8.fdf with every atom moved at random by ~0.03 Ang. Harmonic phonons assume the atoms sit at a
minimum; at any other point the forces are not zero and the curvature is not that of the crystal. With --no-relax:
lower symmetry (many more displacements), spurious modes. With the default pre-relax (positions only, --fmax 0.05):
back to the crystal.
EOF
mlph case7_norelax "$DIR/si8_rattled.fdf" $FT --skip-foundation-comparison --no-relax
mlph case7_relax   "$DIR/si8_rattled.fdf" $FT --skip-foundation-comparison
for c in case7_norelax case7_relax; do
    echo "--- $c"
    section "$(report $c)" "DISPLACEMENTS" | grep -E "^Input forces|^Pre-relax|^Symmetry|^Displacements"
    section "$(report $c)" "BAND STRUCTURE" | grep -E "Band path :|Mesh 8000"
done
pause


banner "Case 8: the quasi-harmonic approximation -- thermal expansion"
cat <<'EOF'
--qha repeats the calculation at several volumes (here 5, over +-3 %), adds each volume's vibrational free energy to
its static energy, and finds the volume that minimizes F(V,T) at every temperature: V(T), the thermal expansion,
B(T), Cp, and the mode Grueneisen parameters gamma = -dln f / dln V. Silicon's transverse acoustic modes have
gamma < 0 -- they soften under compression -- and make the crystal SHRINK on heating below ~120 K.
EOF
echo "\$ stb-mlphonons --file si8.fdf $FT --qha --n-volumes 5 --strain-range 3"
mlph case8_qha "$DIR/si8.fdf" $FT --qha --n-volumes 5 --strain-range 3
section "$(report case8_qha)" "EQUILIBRIUM"
section "$(report case8_qha)" "THERMAL EXPANSION" | sed -n '1,16p'
section "$(report case8_qha)" "MODE GRUENEISEN" | sed -n '1,8p'
echo "Experiment (Okada & Tokumaru 1984): alpha_L(300 K) = 2.6e-6 /K; negative below 124 K."
pause


banner "Case 9: the foundation model alone"
cat <<'EOF'
Without --custom-model the tool uses MACE-MP-0 (--model small|medium|large) and has no card to check. It is
trained on the Materials Project's PBE relaxations, the same functional as the fine-tuned model -- yet its silicon
is soft: README section 5.5.
EOF
mlph case9_foundation "$DIR/si8.fdf" --model medium
section "$(report case9_foundation)" "MODEL"
section "$(report case9_foundation)" "FREQUENCIES" | grep -E "^(Γ|X|L) "
pause


banner "Case 10: guard rails"
try() {   # try <label> <args...>: runs and shows the error line
    echo "--- $1"
    shift
    (cd "$OUT" && stb-mlphonons --no-intro "$@" < /dev/null 2>&1 | strip_colors | grep -E "error|ERROR|FAIL" | head -2) || true
}
try "--custom-model without its card" --file "$DIR/si8.fdf" --custom-model "$MODEL"
try "--qha with 3 volumes (an EOS has 4 parameters)" --file "$DIR/si8.fdf" --qha --n-volumes 3
try "--qha with --animate-mode" --file "$DIR/si8.fdf" --qha --animate-mode
try "--convergence-dims not a multiple of 3" --file "$DIR/si8.fdf" --check-convergence --convergence-dims 1 1 1 2 2
try "--animate-mode with a word it does not know" --file "$DIR/si8.fdf" --animate-mode acoustic
try "an input file that does not exist" --file "$DIR/nothere.fdf"
pause


banner "Case 11: the same tool from the interactive menu (stb-suite, code 7.3)"
cat <<'EOF'
The menu asks the same questions and builds the same command line. Answers piped in: structure, the custom model
(Enter accepts the card next to it), cpu, 2x2x2, output folder menu_out, single volume, no foundation comparison, no
mode freeze, animate a mode (Enter: the first optical mode, at Gamma, default amplitude and frames, xsf), no
display, no gnuplot, no data, save the report.
EOF
mkdir -p "$OUT/case11_menu"
cd "$OUT/case11_menu"
printf '7.3\n%s\n%s\n\n\n\nmenu_out\n1\nn\nn\ny\n\n\n\n\n\n\nn\nn\nn\ny\n\n0\n' "$DIR/si8.fdf" "$MODEL" | \
    stb-suite > menu_log.txt 2>&1 || true
ls -1 menu_out 2>/dev/null | sed 's/^/  /'
strip_colors < menu_out/stb_mlphonons_report.txt 2>/dev/null | grep -E "^Model  |^SHA-256|^Mode  |^Band index|^Status|^Results" | sed 's/^/  /'
cd "$DIR"
pause


banner "Case 12: the full reference set of README section 5"
cat <<'EOF'
reference_results/ holds the reports and data of make_reference.sh: the 2x2x2 and 4x4x4 calculations, the three
foundation models, the 1x1x1..4x4x4 convergence, the QHA with both models, six displacement sizes, the rattled
structure three ways, the elastic constants of the same model (stb-mlelastic), and the SIESTA phonons of the 7.1
example against the model, with the force-noise experiment. results_figures.py turns them into the tables of the
README (numpy + matplotlib only).
EOF
python3 "$DIR/results_figures.py" --no-plots | head -60
echo "   ... (python3 results_figures.py prints every table and redraws figures/)"
pause


banner "Done"
cat <<EOF
Every case wrote into $OUT/<case>/ (report, PNG figures, and whatever was asked for).
README.md explains the theory (section 1), the report (section 3) and every number of the silicon example (5).
EOF
