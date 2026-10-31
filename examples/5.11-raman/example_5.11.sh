#!/bin/bash
# Guided example: Workflow 5.11 -- Raman Spectrum (stb-raman 5.11.1, stb-ramanModes 5.11.2, stb-ramanAnalysis 5.11.3).
#
# Not an automated test (see test/5-workflow/11-raman/{prep,modes,analysis}/test.sh for that) -- a commented
# walk-through: it runs real commands, one case at a time, inside output/, on a COPY of runs/ (finished SIESTA 5.4
# calculations, made by run_example.sh), and shows the piece of each report that proves what just happened.
# runs/ itself is never touched. No SIESTA is needed. Pauses between cases (`yes '' | bash example_5.11.sh` runs it
# straight through). Safe to re-run -- it always starts by wiping its own output/. README section 7 explains every
# number.

set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"
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

section() {   # $1 = file, $2 = the start of a section title ("[1] PHONON MODES")
    strip_colors < "$1" | awk -v t="$2" 'index($0, t) == 1 {on = 1; print; next} on && /^\[[0-9]+[a-z]?\] [A-Z]/ {exit} on {print}' \
        | grep -v "^-----"
}

copy_run() {  # copy_run <study> <dest>: a working copy of runs/<study>, the dielectric functions unpacked
    rm -rf "$2"; mkdir -p "$2"
    cp -a "$DIR/runs/$1/." "$2/"
    find "$2" -name "siesta.EPSIMG.gz" -exec gunzip -f {} \;
}

# run <log> <dir> <command...>: runs a stage inside <dir>, the report saved as <log>
run() {
    local log=$1 dir=$2; shift 2
    echo "  \$ $*"
    (cd "$dir" && "$@" --no-intro < /dev/null 2>&1 | strip_colors > "$log") || true
}


banner "Welcome: the Raman spectrum of a crystal from SIESTA, in three stages"
cat <<'EOF'
Stage 1 (stb-raman) displaces atoms to get the phonons; Stage 2 (stb-ramanModes) finds the Gamma modes and displaces
the cell along each of them, +/- delta, with an optical calculation in every folder; Stage 3 (stb-ramanAnalysis) turns
the change of the dielectric constant into a Raman tensor, an activity and a spectrum. SIESTA runs in between --
here it already ran (runs/).

The tour:
   1  Stage 1 on bulk Si                     5  how converged is the tensor (delta, Optical.Mesh)
   2  Stage 2: modes, symmetry, folders      6  graphene: selection rules, and a semimetal
   3  Stage 3 on the finished Si runs        7  the interactive menu
   4  the diagonal-only trap
EOF
pause


banner "Case 1 -- Stage 1 on bulk Si: the phonon displacements"
cat <<'EOF'
inputs/si holds a relaxed 2-atom Si cell (structure.fdf) and the level of theory (calc.fdf, PBE/TZP, 300 Ry, no
k-grid). A Gamma-point frequency needs no supercell (every cell moves in phase), so -dim 1 1 1; the k-grid of the
force runs is given explicitly. Symmetry leaves ONE displacement to compute, plus the undisplaced cell whose
residual forces Stage 2 subtracts.
EOF
mkdir -p "$OUT/case1"; cp inputs/si/structure.fdf inputs/si/calc.fdf "$OUT/case1/"
run "$OUT/case1/stage1.txt" "$OUT/case1" stb-raman -s structure.fdf -c calc.fdf -p dojo -dim 1 1 1 --kgrid 11 11 11
section "$OUT/case1/stage1.txt" "[4] PHONON DISPLACEMENTS"
ls "$OUT/case1/raman_study/phonon_disp" | sed 's/^/    /'
echo "  -> 1 disp-* folder + perfect_supercell/: two SIESTA runs, about 40 s on 2 cores."
pause


banner "Case 2 -- Stage 2: Gamma modes, their symmetry, and the optical folders"
cat <<'EOF'
On the finished Stage-1 runs (runs/si/raman_study/phonon_disp), Stage 2 builds FORCE_SETS, finds the three optical
modes -- one triply degenerate T2g at 512 cm^-1 -- and writes the +/- delta folders. --full-tensor probes six
directions (x, y, z and the face diagonals), --use-symmetry drops Raman-silent modes (none here), --skip-degenerate
computes one T2g partner for all three: 12 single points instead of 36.
EOF
mkdir -p "$OUT/case2/raman_study"; cp inputs/si/calc.fdf "$OUT/case2/"
cp -a runs/si/raman_study/phonon_disp "$OUT/case2/raman_study/"
run "$OUT/case2/stage2.txt" "$OUT/case2" stb-ramanModes --directory raman_study -c calc.fdf -p dojo --kgrid 11 11 11 \
    --full-tensor --use-symmetry --skip-degenerate --optical-mesh 32 32 32 --optical-broaden 0.1 --optical-nbands 16
section "$OUT/case2/stage2.txt" "[1] PHONON MODES" | grep -E "mode|Residual|Cells"
section "$OUT/case2/stage2.txt" "[1b] SYMMETRY" | grep -E "Space group|Point group|forbidden"
section "$OUT/case2/stage2.txt" "[4] OPTICAL" | grep -E "Degenerate|covers|Folders to write"
echo
echo "  Without --kgrid: the template has no k-grid, so one is suggested and forced, with a warning --"
echo "  SIESTA alone would sample Gamma only:"
rm -rf "$OUT/case2/raman_study/optical_disp"
run "$OUT/case2/stage2_nokgrid.txt" "$OUT/case2" stb-ramanModes --directory raman_study -c calc.fdf \
    --full-tensor --use-symmetry --skip-degenerate
grep -E "SCF k-grid forced|WARNING\] The template defines no" "$OUT/case2/stage2_nokgrid.txt" | fold -s -w 116 | sed 's/^/    /'
pause


banner "Case 3 -- Stage 3 on the finished Si runs: tensor, activity, spectrum"
cat <<'EOF'
The 12 optical runs of runs/si (Optical.Mesh 32^3, broadening 0.1 eV, 16 bands). For each folder, eps1(0) comes from
the Kramers-Kronig transform of eps2; R = [eps(+delta) - eps(-delta)] / 2 delta along every probe direction.
EOF
copy_run si "$OUT/case3"
run "$OUT/case3/stage3.txt" "$OUT/case3" stb-ramanAnalysis --directory raman_study --save-gnuplot
section "$OUT/case3/stage3.txt" "[2] RAMAN-ACTIVE" | grep -E "Mode|^  1 "
section "$OUT/case3/stage3.txt" "[3] SPECTRUM" | grep -E "Peaks|^  5"
cat <<'EOF'
  -> The diagonal is numerical noise (10^-3 1/Ang); all of the response is in Rxz and Rxy: the T2g tensor has only
     off-diagonal entries in the cubic axes. rho = 0.750 exactly (depolarized), one line at 512 cm^-1
     (experiment about 520 cm^-1: PBE is 1.5 % soft). The spectrum is in output/case3/raman_study/plot/.
EOF
pause


banner "Case 4 -- The diagonal-only trap"
cat <<'EOF'
The default of stb-ramanModes probes x, y and z only (6 folders per mode). For silicon that sees nothing: the same
6 calculations as in case 3 (their inputs are identical, so their outputs are reused here), analysed alone.
EOF
mkdir -p "$OUT/case4/raman_study"; cp inputs/si/calc.fdf "$OUT/case4/"
cp -a runs/si/raman_study/phonon_disp "$OUT/case4/raman_study/"
run "$OUT/case4/stage2.txt" "$OUT/case4" stb-ramanModes --directory raman_study -c calc.fdf -p dojo --kgrid 11 11 11 \
    --use-symmetry --skip-degenerate --optical-mesh 32 32 32 --optical-broaden 0.1 --optical-nbands 16
for f in "$OUT"/case4/raman_study/optical_disp/mode_01_*; do
    src="$OUT/case3/raman_study/optical_disp/$(basename "$f")"
    # same geometry and same optical settings (runs/ predates the forced Diag.ParallelOverK line)
    cmp -s "$f/structure.fdf" "$src/structure.fdf" \
        && diff -q <(grep -v ParallelOverK "$f/config_extra.fdf") <(grep -v ParallelOverK "$src/config_extra.fdf") > /dev/null \
        || { echo "  [unexpected] $(basename "$f") differs from case 3"; exit 1; }
    cp "$src/calc.out" "$src/siesta.EPSIMG" "$f/"
done
run "$OUT/case4/stage3.txt" "$OUT/case4" stb-ramanAnalysis --directory raman_study
grep -E "activity~" "$OUT/case4/stage3.txt" | head -1 | sed 's/^/  /'
grep -E "^\[WARNING\] Mode" "$OUT/case4/stage3.txt" | fold -s -w 116 | sed 's/^/  /'
echo "  -> Activity 2e-4 instead of 2891, and labelled 'polarized': without the warning it reads as 'no Raman line'."
pause


banner "Case 5 -- How converged is the tensor: displacement and Optical.Mesh"
cat <<'EOF'
The same Stage 2 with delta = 0.01 and 0.04 Ang (runs/si_d001, si_d004) and with Optical.Mesh 16^3 (si_mesh16).
EOF
printf "  %-30s %10s %10s %12s\n" "study" "|Rxz|" "|Rxy|" "activity"
for s in si_d001 si si_d004 si_mesh16; do
    copy_run "$s" "$OUT/case5/$s"
    run "$OUT/case5/$s/stage3.txt" "$OUT/case5/$s" stb-ramanAnalysis --directory raman_study
    awk -v s="$s" '/^  1 / {printf "  %-30s %10.3f %10.3f %12.1f\n", s, ($8<0?-$8:$8), ($7<0?-$7:$7), $10}' \
        "$OUT/case5/$s/stage3.txt"
done
echo "  -> README section 7.5 reads this table."
pause


banner "Case 6 -- Graphene: selection rules, and a semimetal"
cat <<'EOF'
The graphene monolayer of inputs/graphene (2x2x1 supercell for Stage 1; full tensor for ALL three optical modes, no
symmetry shortcut, so the rules can be checked rather than assumed): ZO (out of plane; B2g, B1g in Phonopy's axes) is Raman-silent, the two
E2g partners (the G line) must have the same activity.
EOF
copy_run graphene "$OUT/case6"
section "$OUT/case6/raman_study/raman_stage2.txt" "[1] PHONON MODES" | grep -E "mode" | sed 's/^/ /'
run "$OUT/case6/stage3.txt" "$OUT/case6" stb-ramanAnalysis --directory raman_study --save-gnuplot
section "$OUT/case6/stage3.txt" "[2] RAMAN-ACTIVE" | grep -E "Mode|^  [123] "
grep -E "^\[WARNING\] Degenerate" "$OUT/case6/stage3.txt" | fold -s -w 116 | sed 's/^/  /'
cat <<'EOT'
  -> What holds: ZO is silent (activity 7e-4) and the E2g tensors live in the plane (Rzz = Rxz = Ryz = 0).
     What does not: an E2g tensor must be traceless and both partners equal -- here they differ by 28x, because
     graphene has no gap and eps1(0) is not defined (README section 7.6).
EOT
pause


banner "Case 7 -- The interactive menu (stb-suite, code 5.11.3)"
cat <<'EOF'
The menu asks the same questions as the flags and runs the same command: here Stage 3 on the Si runs of case 3.
EOF
# directory, output file name, linewidth, temperature, experimental file, gnuplot, view, Enter, 0 = exit
(cd "$OUT/case3" && printf '5.11.3\nraman_study\n\n\n\n\nn\nn\n\n0\n' | stb-suite 2>&1 | strip_colors > "$OUT/case7_menu.txt") || true
grep -E "^  1 |Peaks detected" "$OUT/case7_menu.txt" | head -3 | sed 's/^/  /'
echo


banner " Done"
cat <<EOF
Everything is in $OUT/. What to take from the tour:
  - The Raman tensor is the derivative of the dielectric constant along a mode; its symmetry shape decides how many
    probe directions you need. For a cubic T2g mode, the default diagonal-only run sees nothing: use --full-tensor.
  - --use-symmetry and --skip-degenerate are exact (group theory); for Si they cut 36 runs to 12.
  - Frequencies come from the phonons (Stage 1); intensities from the optical runs, which must be converged in
    Optical.Mesh -- and are only qualitative for a zero-gap system like graphene.
EOF
