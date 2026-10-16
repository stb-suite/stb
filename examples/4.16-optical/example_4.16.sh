#!/bin/bash
# Guided example: Workflow 4.16 -- Optical Properties (both stages: stb-optical, code 4.16.1,
# and stb-opticalAnalysis, code 4.16.2, in the stb-suite menu).
#
# Not an automated test (see test/4-workflow/16-optical/{prep,analysis}/test.sh for that) -- a
# commented walk-through: it runs real commands, one case at a time, into its own output/<case>/
# folder, and shows you the piece of output that proves what just happened. Pauses between
# sections so you can read before moving on. Safe to re-run any time -- it always starts by
# wiping its own output/.
#
# Stage 1 (stb-optical) runs for real in every case below (cheap -- it only writes .fdf files, no
# SCF). Stage 2 (stb-opticalAnalysis) runs on REAL SIESTA results: reference_results/ holds the
# finished runs of bulk silicon (the 6 direction folders and a convergence study) and of a graphene
# monolayer (xx and zz), computed with the Stage 1 commands of this script and the inputs in
# this folder (SIESTA 5.4.2, mpirun -np 2; see README section 6 to regenerate them). Nothing in
# Stage 2 is fabricated.

set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

# --show needs a display; MPLBACKEND=Agg keeps this script headless-safe (the figure is still
# saved with --png).
export MPLBACKEND=Agg

OUT="$DIR/output"
rm -rf "$OUT"
mkdir -p "$OUT"
REF="$DIR/reference_results"

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


banner "Welcome: how a crystal answers to light, from a ground-state DFT calculation"
cat <<'EOF'
A solid answers an oscillating electric field with a complex dielectric function
eps(E) = eps1 + i*eps2: eps2 is absorption, eps1 dispersion, and from the pair come the
refractive index, absorption coefficient, reflectivity, energy-loss function and optical
conductivity. SIESTA computes eps2 from interband transitions between the Kohn-Sham bands
(independent-particle approximation); Kramers-Kronig gives eps1; everything else follows.

Stage 1 (stb-optical) prepares one SIESTA folder per polarization direction of an ALREADY
RELAXED structure. Stage 2 (stb-opticalAnalysis) reads the finished runs and reports the
properties AND whether to trust them: sum rules, absorption edge, plasmon, anisotropy checked
against the crystal symmetry, a what-if scissor shift, and a dashboard.

This tour: Stage 1 on silicon, graphene and a water molecule (Cases 1-7), then Stage 2 on
REAL SIESTA results for silicon and graphene (Cases 8-14). See README.md for the theory (section 1)
and the worked, convergence-checked example (section 5).
EOF
pause


banner "Case 1: Stage 1 -- bulk silicon, what stb-optical reports"
cat <<'EOF'
structure.fdf is the already-relaxed 2-atom silicon cell; calc.fdf is the ground-state template
(PBE, TZP, small PAO.EnergyShift) with NO k-grid and NO MD directives. The tool needs only
-f; -c defaults to ./calc.fdf. Six directions: the 3 diagonal ones and the 3 biaxial ones that
give the off-diagonal components. --kgrid sets the SCF k-grid explicitly (Case 2 shows what
happens without it); the optical mesh, broadening and number of bands are the three optical knobs.
Watch [1] (structure and valence electrons), [3] (the optical mesh assessed) and [5] (the symmetry
analysis: for a cubic crystal 1 of the 6 directions is enough).
EOF
mkdir -p "$OUT/case1-si"
cp structure.fdf calc.fdf "$OUT/case1-si/"
echo "\$ stb-optical -f structure.fdf -p dojo --directions xx yy zz xy xz yz --kgrid 11 11 11 \\"
echo "      --optical-mesh 32 32 32 --optical-broaden 0.1 --optical-nbands 16 --no-intro"
(cd "$OUT/case1-si" && stb-optical -f structure.fdf -p dojo --directions xx yy zz xy xz yz --kgrid 11 11 11 \
    --optical-mesh 32 32 32 --optical-broaden 0.1 --optical-nbands 16 --no-intro > stage1.log 2>&1)
strip_colors < "$OUT/case1-si/stage1.log" | sed -n '/\[1\] INPUT STRUCTURE/,/\[2\] GROUND STATE/p' | sed '$d'
strip_colors < "$OUT/case1-si/stage1.log" | sed -n '/\[3\] OPTICAL PARAMETERS/,/\[4\] DIRECTIVES/p' | sed '$d' | cut -c1-150
strip_colors < "$OUT/case1-si/stage1.log" | sed -n '/What crystal symmetry says/,/\[6\] FOLDERS/p' | sed '$d' | cut -c1-175
pause


banner "Case 2: config_extra.fdf -- everything mandatory goes there, your template is not edited"
cat <<'EOF'
Each folder gets a config_extra.fdf, and calc.fdf starts with '%include config_extra.fdf' followed
by your template byte for byte. SIESTA's fdf reader is first-occurrence-wins, so the forced
directives beat anything the template says. Here: the SCF k-grid (required -- SIESTA would sample
Gamma only), the single-point guards, and the optical block with this folder's own polarization.
EOF
F="$OUT/case1-si/optical_study/dir_xy"
ls "$OUT/case1-si/optical_study" | tr '\n' ' '; echo
echo "\$ cat optical_study/dir_xy/config_extra.fdf          # the biaxial xy folder"
cat "$F/config_extra.fdf"
echo
echo "\$ head -3 optical_study/dir_xy/calc.fdf"
head -3 "$F/calc.fdf"
echo
echo "--- Without --kgrid: the template has none, so a grid is suggested (density 0.2 1/Ang) and forced"
mkdir -p "$OUT/case2-kgrid"
cp structure.fdf calc.fdf "$OUT/case2-kgrid/"
(cd "$OUT/case2-kgrid" && stb-optical -f structure.fdf -p dojo --directions xx --no-intro 2>&1 \
    | strip_colors | grep -E "SCF k-grid forced|kgrid.Monk" | cut -c1-175)
echo
echo "--- A template that already has a grid keeps it (nothing forced); --kgrid-density overrides it"
{ cat calc.fdf; printf 'kgrid.MonkhorstPack [6 6 6]\nMD.NumCGsteps 40\nMD.VariableCell T\n'; } > "$OUT/case2-kgrid/calc_own.fdf"
(cd "$OUT/case2-kgrid" && stb-optical -f structure.fdf -c calc_own.fdf --directions xx -O keep --no-intro 2>&1 \
    | strip_colors | grep -E "SCF k-grid forced|MD.VariableCell|MD.NumCGsteps" | cut -c1-150)
(cd "$OUT/case2-kgrid" && stb-optical -f structure.fdf -c calc_own.fdf --directions xx -O dens --kgrid-density 0.1 --no-intro 2>&1 \
    | strip_colors | grep -E "SCF k-grid forced" | cut -c1-150)
echo
echo "Even with MD.NumCGsteps 40 and MD.VariableCell T in the template, the folder runs one single point:"
grep -E "MD.Steps|MD.NumCGsteps|MD.VariableCell" "$OUT/case2-kgrid/keep/dir_xx/config_extra.fdf"
pause


banner "Case 3: the advisories -- what Stage 1 warns about BEFORE you spend CPU time"
cat <<'EOF'
Optical.Mesh too coarse (spiky eps2, unreliable peaks), Optical.NumberOfBands not exceeding the
occupied bands (no conduction states, so no absorption at all), and a broadening far below what the
mesh can resolve. Each is printed where it applies and repeated in a summary at the end.
EOF
mkdir -p "$OUT/case3-advisories"
cp structure.fdf calc.fdf "$OUT/case3-advisories/"
echo "\$ stb-optical -f structure.fdf --directions xx --optical-mesh 8 8 8 --optical-nbands 3 --optical-broaden 0.01"
(cd "$OUT/case3-advisories" && stb-optical -f structure.fdf --directions xx --optical-mesh 8 8 8 \
    --optical-nbands 3 --optical-broaden 0.01 --no-intro 2>&1 | strip_colors \
    | sed -n '/warning(s) above/,$p' | grep -E "warning|^  - " | cut -c1-170)
pause


banner "Case 4: symmetry decides how many directions you need"
cat <<'EOF'
The dielectric tensor must be invariant under the crystal's symmetry. For silicon (cubic) there is
1 independent component; for a hexagonal graphene monolayer 2 (in-plane xx = yy, and zz). The tool
finds this in the Cartesian frame of YOUR cell and tells you the minimal direction set.
EOF
mkdir -p "$OUT/case4-symmetry"
cp graphene.fdf calc_graphene.fdf "$OUT/case4-symmetry/"
for spec in "$DIR/structure.fdf $DIR/calc.fdf silicon" "graphene.fdf calc_graphene.fdf graphene"; do
    set -- $spec
    echo "------------------------------------------------------------------"
    echo "\$ stb-optical -f $(basename $1) -c $(basename $2) --directions xx yy zz --no-intro     # $3"
    (cd "$OUT/case4-symmetry" && stb-optical -f "$1" -c "$2" --directions xx yy zz -O "study_$3" --no-intro 2>&1 \
        | strip_colors | grep -E "Space group|Independent|Non-zero|Zero by|Minimal|NOTE\] " | cut -c1-175)
done
pause


banner "Case 5: dimensionality -- bulk, 2D slab and 0D molecule are told apart"
cat <<'EOF'
A vacuum-padded cell dilutes the dielectric response (README 1.9). Stage 1 detects vacuum axes and
prints the matching limitation, the single SCF k-point division along each vacuum axis, and what
Stage 2 can do about it (--dimensionality-correction).
EOF
mkdir -p "$OUT/case5-dim"
cp graphene.fdf h2o.fdf calc_graphene.fdf "$OUT/case5-dim/"
cat > "$OUT/case5-dim/calc_h2o.fdf" <<'EOF'
SystemLabel siesta
%include structure.fdf
PAO.BasisSize DZP
Mesh.Cutoff 200 Ry
EOF
for spec in "graphene.fdf calc_graphene.fdf" "h2o.fdf calc_h2o.fdf"; do
    set -- $spec
    echo "------------------------------------------------------------------"
    echo "\$ stb-optical -f $1 -c $2 --directions xx zz --no-intro"
    (cd "$OUT/case5-dim" && stb-optical -f "$1" -c "$2" --directions xx zz -O "study_${1%.fdf}" --no-intro 2>&1 \
        | strip_colors | grep -E "^Detected|KNOWN LIMITATION|SCF k-grid forced|more than 1 division" | cut -c1-210)
done
pause


banner "Case 6: validation -- what stb-optical correctly refuses"
cat <<'EOF'
A missing structure, a missing calc.fdf (the default template that is not there), and an
unknown direction -- each exits with a clear message and a non-zero status.
EOF
mkdir -p "$OUT/case6-validation"
cp structure.fdf "$OUT/case6-validation/"
cd "$OUT/case6-validation"
echo "\$ stb-optical -f nothere.fdf -c $DIR/calc.fdf"
stb-optical -f nothere.fdf -c "$DIR/calc.fdf" --no-intro 2>&1 | strip_colors | grep ERROR || true
echo
echo "\$ stb-optical -f structure.fdf          # no ./calc.fdf here"
stb-optical -f structure.fdf --no-intro 2>&1 | strip_colors | grep ERROR || true
echo
echo "\$ stb-optical -f structure.fdf -c $DIR/calc.fdf --directions bogus"
stb-optical -f structure.fdf -c "$DIR/calc.fdf" --directions bogus --no-intro 2>&1 | tail -1 | strip_colors || true
cd "$DIR"
pause


banner "Case 7: CLI vs. the interactive stb-suite menu (4.16.1) -- same result"
cat <<'EOF'
The menu asks the same questions as the flags. Reproducing Case 1's command through
stb-suite -> 4.16.1 (structure, calc.fdf default, dojo pseudos, six directions, output name,
advanced settings: SCF k-grid, optical mesh, broadening, bands) and diffing the generated
config_extra.fdf files against Case 1's proves the two paths equal.
EOF
mkdir -p "$OUT/case7-menu"
cp structure.fdf calc.fdf "$OUT/case7-menu/"
(cd "$OUT/case7-menu" && printf '4.16.1\nstructure.fdf\n\n1\nxx yy zz xy xz yz\noptical_study\ny\n11 11 11\n32 32 32\n0.1\n16\n\n0\n' \
    | stb-suite > menu.log 2>&1)
if diff -rq -x '*.psml' "$OUT/case1-si/optical_study" "$OUT/case7-menu/optical_study" > /dev/null 2>&1; then
    echo "CLI (Case 1) and interactive menu (4.16.1) produced identical direction folders."
else
    echo "UNEXPECTED: CLI and interactive menu folders differ -- see"
    echo "  diff -rq -x '*.psml' '$OUT/case1-si/optical_study' '$OUT/case7-menu/optical_study'"
fi
pause


banner "Case 8: Stage 2 on REAL SIESTA data -- bulk silicon"
cat <<'EOF'
reference_results/si/optical_study/ holds the 6 finished direction folders of Case 1 (SIESTA 5.4.2,
1-2 minutes each). Stage 2 detects each folder's direction from its own Optical.Vector,
reads the .EPSIMG, and does the rest. No argument is needed beyond --directory; here also
--png (dashboard) and --scissor (the what-if of Case 11). Report sections, in turn:
EOF
mkdir -p "$OUT/si"
cp -r "$REF/si/optical_study" "$OUT/si/optical_study"
echo "\$ stb-opticalAnalysis --directory output/si/optical_study --png --scissor 0.55 --no-intro"
stb-opticalAnalysis --directory "$OUT/si/optical_study" --png --scissor 0.55 --no-intro 2>&1 \
    | strip_colors > "$OUT/si_stage2.log"
echo
echo "----- system and calculation, as read back from the folders -----"
sed -n '/SYSTEM & CALCULATION/,/DATA QUALITY/p' "$OUT/si_stage2.log" | sed '$d' | sed '1,2d' | cut -c1-140
pause
echo "----- data quality: is the spectrum complete? (the f-sum rule) -----"
sed -n '/DATA QUALITY & CONVERGENCE/,/STATIC RESPONSE/p' "$OUT/si_stage2.log" | sed '$d' | sed '1,2d' | cut -c1-175
pause


banner "Case 9: static response, absorption edge, peaks and plasmon"
echo "----- static constants and the edge -----"
sed -n '/STATIC RESPONSE & ABSORPTION EDGE/,/SPECTRAL FEATURES/p' "$OUT/si_stage2.log" | sed '$d' | sed '1,2d' | cut -c1-175
echo "----- the peaks of eps2 and the plasmon -----"
sed -n '/SPECTRAL FEATURES/,/ANISOTROPY/p' "$OUT/si_stage2.log" | sed '$d' | sed '1,2d' | cut -c1-175
pause


banner "Case 10: anisotropy, symmetry check and the reconstructed off-diagonal tensor"
cat <<'EOF'
Silicon is cubic: eps_xx = eps_yy = eps_zz and eps_ij = 0. Stage 2 compares the finished spectra with
what the symmetry demands -- an end-to-end test of the whole calculation (all 6 directions
must agree to numerical noise, the bisector runs must reconstruct a vanishing eps_ij).
EOF
sed -n '/ANISOTROPY & SYMMETRY/,/OPTICAL CONSTANTS AT SELECTED/p' "$OUT/si_stage2.log" | sed '$d' | sed '1,2d' | cut -c1-175
pause


banner "Case 11: comparison with experiment, and the scissor what-if"
cat <<'EOF'
Independent-particle optics on a semilocal DFT ground state is red-shifted by about the gap error
(PBE gap of Si ~0.6 eV vs 1.17 eV) and overestimates eps1(0). The scissor what-if shifts eps2
rigidly by --scissor eV (eps1 re-derived by Kramers-Kronig) to show where the peaks would land.
Experiment (README 5.5): strongest peak (E2) ~4.2-4.3 eV, eps1(0) ~11.7.
EOF
sed -n '/SCISSOR WHAT-IF/,/INTERPRETATION/p' "$OUT/si_stage2.log" | sed '$d' | sed '1,2d' | cut -c1-175
echo "----- the plain-language interpretation -----"
sed -n '/INTERPRETATION/,/SUMMARY & FILES/p' "$OUT/si_stage2.log" | sed '$d' | sed '1,2d' | cut -c1-200
pause


banner "Case 12: convergence -- which knob moves which number (Optical.Mesh, basis, bands, broadening)"
cat <<'EOF'
reference_results/si/convergence/ holds 11 single-direction (xx) re-runs of the production calculation,
each changing ONE thing: the optical mesh (16/24/48 vs 32), PAO.EnergyShift (0.02/0.01/0.002 vs 0.005),
the basis (DZP vs TZP), the number of bands (26 = all), the broadening (0.05/0.2 vs 0.1) and the
mesh cut-off (450 vs 300 Ry). Their EPSIMG files are trimmed to 30 eV to keep the repository small,
so every row below is analysed with --energy-max 30 (the production row too).
EOF
mkdir -p "$OUT/conv"
for v in mesh16 mesh24 mesh48 shift020 shift010 shift002 dzp allbands broaden005 broaden020 cut450; do
    cp -r "$REF/si/convergence/$v" "$OUT/conv/$v"
done
mkdir -p "$OUT/conv/production/dir_xx"
cp "$REF/si/optical_study/dir_xx/"{structure.fdf,calc.fdf,config_extra.fdf,calc.out,siesta.EPSIMG} "$OUT/conv/production/dir_xx/"
for v in production mesh16 mesh24 mesh48 shift020 shift010 shift002 dzp allbands broaden005 broaden020 cut450; do
    stb-opticalAnalysis --directory "$OUT/conv/$v" --energy-max 30 --no-intro 2>&1 | strip_colors > "$OUT/conv/${v}_stage2.log"
done
python3 "$DIR/convergence_table.py" "$OUT/conv" production mesh16 mesh24 mesh48 shift020 shift010 shift002 dzp allbands broaden005 broaden020 cut450
cat <<'EOF'

Reading it: eps1(0) and the plasmon barely move with anything but the basis. The PEAK POSITION depends on
the optical mesh: 16 and 24 are not converged (3.66 and 3.83 eV), 32 and 48 agree (3.75, 3.74), and the peak
HEIGHT also depends on the broadening (0.05 -> 62, 0.2 -> 48). The basis (PAO.EnergyShift, DZP vs TZP) moves
eps1(0) by a few percent and the f-sum check, but at a converged mesh not the peak positions; the number of
bands and the mesh cut-off change essentially nothing. README 5.3-5.4 discusses each.
EOF
pause


banner "Case 13: a 2D material for real -- graphene and the dimensionality correction"
cat <<'EOF'
reference_results/graphene/ holds the xx and zz folders of a graphene monolayer (20 Ang cell along c,
~16.6 Ang of vacuum; GGA-PBE, TZP, virtual_vault carbon). The raw supercell response is diluted by the
vacuum and scales with the cell height; --dimensionality-correction restores the intrinsic 2D function
(README 1.9, equation 14) and needs --thickness: here 3.35 Ang, the graphite interlayer distance (a choice,
not a measurement -- see the README). Graphene is a zero-gap semimetal: no absorption edge, and in the
visible a free sheet absorbs the universal pi*alpha = 2.3 % of the light.
EOF
mkdir -p "$OUT/graphene"
cp -r "$REF/graphene/optical_study" "$OUT/graphene/optical_study"
echo "\$ stb-opticalAnalysis --directory output/graphene/optical_study --energy-max 15 --dimensionality-correction --thickness 3.35 --png"
stb-opticalAnalysis --directory "$OUT/graphene/optical_study" --energy-max 15 --dimensionality-correction \
    --thickness 3.35 --png --no-intro 2>&1 | strip_colors > "$OUT/graphene_stage2.log"
sed -n '/STATIC RESPONSE & ABSORPTION EDGE/,/SPECTRAL FEATURES/p' "$OUT/graphene_stage2.log" | sed '$d' | sed '1,2d' | head -8 | cut -c1-170
echo "..."
sed -n '/DIMENSIONALITY CORRECTION/,/INTERPRETATION/p' "$OUT/graphene_stage2.log" | sed '$d' | sed '1,2d' | cut -c1-170
echo
echo "The out-of-plane (zz) value depends violently on the thickness you choose -- it can even turn unphysical:"
for t in 3.35 6.7; do
    printf "  --thickness %-5s " "$t"
    stb-opticalAnalysis --directory "$OUT/graphene/optical_study" --energy-max 15 --dimensionality-correction \
        --thickness "$t" -o thick_"$t" --no-intro 2>&1 | strip_colors | grep "Direction zz (perpendicular)" | sed 's/^ *//'
done
echo
echo "Absorbance of the sheet, A(E) = E * t * eps2_2D(E) / (hbar c)  (thin-sheet limit), from the corrected in-plane spectrum:"
python3 - "$OUT/graphene/optical_study" <<'PYEOF'
import sys, numpy as np
d = np.loadtxt(f"{sys.argv[1]}/optical_results_xx_2Dcorrected.dat")
E, eps2 = d[:, 0], d[:, 2]
hbar_c_eV_nm = 197.3269804
t_nm = 0.335
for e in (1.0, 1.5, 2.0, 3.0):
    a = e * t_nm * np.interp(e, E, eps2) / hbar_c_eV_nm
    print(f"  E = {e:.1f} eV   A = {100 * a:.2f} %")
print("  universal value pi*alpha = 2.29 %  (alpha = 1/137.036)")
PYEOF
pause


banner "Case 14: guard rails and the interactive Stage 2 menu"
cat <<'EOF'
Stage 2 flags unreliable inputs from the files alone: a missing .EPSIMG (the folder is skipped),
no folder at all, an SCF that never confirmed convergence, folders run with different optical
settings, and a symmetry-forbidden component that is not zero.
EOF
echo "--- (a) folders with no EPSIMG yet"
mkdir -p "$OUT/guard_noeps/dir_xx"
cp "$REF/si/optical_study/dir_xx/"{structure.fdf,calc.fdf,config_extra.fdf} "$OUT/guard_noeps/dir_xx/"
stb-opticalAnalysis --directory "$OUT/guard_noeps" --no-intro 2>&1 | strip_colors | grep -E "Could not read|No usable" | sed "s#$DIR/##" | cut -c1-150 || true
echo
echo "--- (b) an SCF that never confirmed convergence (the line removed from dir_yy's calc.out)"
cp -r "$REF/si/optical_study" "$OUT/guard_scf"
sed -i '/SCF cycle converged after/d' "$OUT/guard_scf/dir_yy/calc.out"
stb-opticalAnalysis --directory "$OUT/guard_scf" --no-intro 2>&1 | strip_colors | grep -E "SCF convergence not confirmed" | cut -c1-150 | head -2
echo
echo "--- (c) dir_zz re-run with a different Optical.Mesh: the comparison is no longer legitimate"
cp -r "$REF/si/optical_study" "$OUT/guard_mesh"
sed -i 's/BZ mesh dimensions *= *32 *32 *32/BZ mesh dimensions =   16  16  16/' "$OUT/guard_mesh/dir_zz/calc.out"
stb-opticalAnalysis --directory "$OUT/guard_mesh" --no-intro 2>&1 | strip_colors | grep -E "NOT run with identical" | cut -c1-170
echo
echo "--- the interactive menu (4.16.2): directory, default log name, eps plot, no dimensionality correction,"
echo "    no screen dashboard, save the PNG, scissor 0.55 eV, default energies, no advanced settings"
cp -r "$REF/si/optical_study" "$OUT/menu_study"
(cd "$OUT" && printf '4.16.2\nmenu_study\n\n\n\n\ny\n0.55\n\n\n\n0\n' | stb-suite > menu2.log 2>&1)
if [ -f "$OUT/menu_study/optical_results_summary.png" ] && [ -f "$OUT/menu_study/optical_results_avg_scissor.dat" ]; then
    echo "The menu produced the dashboard PNG and the scissor spectrum, like the CLI run."
else
    echo "UNEXPECTED: the menu run did not write the dashboard/scissor files -- see output/menu2.log"
fi
pause


banner "Workflow 4.16 complete"
cat <<EOF
Stage 1 (stb-optical, 4.16.1): one folder per polarization direction -- structure.fdf, your untouched
calc.fdf, a config_extra.fdf with the SCF k-grid, the single-point guards and the optical block, and the
pseudopotentials -- and a report that checks the set-up (mesh, bands, broadening) and tells you the minimal
direction set from the crystal symmetry.
Stage 2 (stb-opticalAnalysis, 4.16.2): per-folder data quality (SCF, the f-sum rule), static constants,
absorption edge, peaks and plasmon, anisotropy and symmetry consistency, the reconstructed off-diagonal
tensor, optical constants at chosen energies, a scissor what-if, a dimensionality correction for 2D/0D,
and CSV, gnuplot and a 9-panel dashboard (--png / --show).

On YOUR crystal:
  # relax it first (stb-mlrelax or a SIESTA CG run), then
  stb-optical -f relaxed.fdf -c calc.fdf -p dojo --directions xx yy zz --optical-mesh 24 24 24
  for d in optical_study/dir_*; do (cd \$d && mpirun -np 4 siesta calc.fdf > calc.out); done
  stb-opticalAnalysis --png --energy-max 30
Read the [WARNING] lines of both reports first, and converge the optical mesh and the basis --
README section 5.3-5.4 shows how much they matter.
EOF
