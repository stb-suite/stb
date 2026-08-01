#!/bin/bash
# Guided example: Workflow 4.3 -- Cohesive Energy (both stages: stb-cohesive,
# code 4.3.1, and stb-cohesiveAnalysis, code 4.3.2, in the stb-suite menu).
#
# Not an automated test (see test/4-workflow/3-cohesive/{prep,analysis}/
# test.sh for that) -- a commented walk-through: it runs real commands, one
# case at a time, into its own output/<case>/ folder, and shows you the
# piece of output that proves what just happened. Pauses between sections so
# you can read before moving on. Safe to re-run any time -- it always starts
# by wiping its own output/.
#
# structure.fdf is the REAL, relaxed GaAs structure this workflow was
# actually investigated against this session (README section 4). Unlike
# 4.1-strain/4.2-elastic, there is no separate calc.fdf here -- stb-cohesive
# always builds its own single-point SCF template internally (README
# section 2.2), so there's nothing to check in.
#
# Stage 1 (stb-cohesive) runs for real in every case below (cheap -- it only
# writes .fdf files, no SCF). Stage 2 (stb-cohesiveAnalysis) is demonstrated
# against small, hand-built calc.out data (same idea as 4.1-strain's and
# 4.2-elastic's own scripts) since this folder doesn't invoke real SIESTA --
# the target numbers are illustrative, NOT this material's real cohesive
# energy (see the README's own section 4 for that, with the real,
# physically-verified numbers this whole workflow was actually developed
# against).

set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

# Stage 2's --view calls plt.show() -- MPLBACKEND=Agg makes that a no-op
# instead of blocking on a GUI window, same convention test.sh itself uses.
export MPLBACKEND=Agg

OUT="$DIR/output"
rm -rf "$OUT"
mkdir -p "$OUT"

pause() {
    echo
    read -p "  [Press Enter to continue] " -r
    echo
}

# Writes a minimal, real-format SIESTA calc.out (just what
# core/siesta_log.py's get_free_energy/get_scf_convergence look for) into
# $1/calc.out, with FreeEng $2 and SCF iteration count $3 (default 5).
# Illustrative numbers only -- see the header comment above.
write_atom_out() {
    local folder="$1" free_eng="$2" iters="${3:-5}"
    mkdir -p "$folder"
    cat > "$folder/calc.out" <<EOF
SCF cycle converged after ${iters} iterations
siesta: FreeEng =    ${free_eng}
EOF
}


echo "=================================================================="
echo " Welcome: why isolated-atom LCAO cohesive energies need a correction"
echo "=================================================================="
cat <<'EOF'
Cohesive energy = (E_bulk - sum of isolated-atom energies) / N atoms --
simple arithmetic, EXCEPT that localized-basis (LCAO/PAO) DFT systematically
over-binds unless corrected: an atom computed alone only has its OWN basis
functions to relax into, while that same atom inside the solid also borrows
variational freedom from its neighbors' orbitals. stb-cohesive's BSSE
(counterpoise) correction fixes this by re-computing the isolated atom
surrounded by "ghost" copies of its real neighbors -- same basis, zero real
physics (no charge, no nuclear attraction).

Stage 1 (stb-cohesive) writes structure/, atoms/<species>/, and (by default)
atoms_bsse/<species>/ -- plus, optionally, one or more larger-cutoff
references to check the correction has actually converged. Stage 2
(stb-cohesiveAnalysis) reads the finished SIESTA runs back and reports the
cohesive energy, corrected and uncorrected, side by side.

See this folder's README.md for the full theory (section 1) and a real,
physically-verified worked example (section 4).
EOF
pause


echo "=================================================================="
echo " Case 1: Stage 1, default run (BSSE correction ON by default)"
echo "=================================================================="
cat <<'EOF'
structure.fdf is real, relaxed GaAs (2 atoms/cell). With no flags beyond
the pseudopotential source, stb-cohesive writes structure/, one atoms/<sp>/
per species, AND one atoms_bsse/<sp>/ per species (--bsse-correction is ON
by default) -- watch [4]/[5] below, then the folder layout on disk.
EOF
mkdir -p "$OUT/case1-default"
cp structure.fdf "$OUT/case1-default/"
echo "\$ stb-cohesive -s structure.fdf -p dojo --no-intro"
(cd "$OUT/case1-default" && stb-cohesive -s structure.fdf -p dojo --no-intro \
    | sed -n '/\[4\] ISOLATED ATOMS/,/\[6\] PSEUDOPOTENTIALS/p')
echo
echo "On disk:"
(cd "$OUT/case1-default" && find cohesive_runs -maxdepth 2 -type d | sort)
pause


echo "=================================================================="
echo " Case 2: the vacuum-buffer guarantee, live (README sections 1.6/2.3)"
echo "=================================================================="
cat <<'EOF'
A real bug lived exactly here: the ghost-cluster box used to be a FIXED
--vacuum-sized cube no matter how large --bsse-cutoff got, so a larger
convergence-check cutoff silently ate into the same fixed buffer. Fixed:
box side = 2*cutoff + --vacuum, so the buffer beyond the outermost ghost
atom stays constant. Proof, on this real structure: --bsse-cutoff 4 with a
+2 Ang convergence check (cutoff 6) should give boxes of 28.0 and 32.0 Ang
-- different sizes, SAME buffer.
EOF
mkdir -p "$OUT/case2-vacuum-buffer"
cp structure.fdf "$OUT/case2-vacuum-buffer/"
echo "\$ stb-cohesive -s structure.fdf -p dojo --bsse-cutoff 4 --bsse-convergence-check --bsse-convergence-increment 2 --no-intro"
(cd "$OUT/case2-vacuum-buffer" && stb-cohesive -s structure.fdf -p dojo --bsse-cutoff 4 \
    --bsse-convergence-check --bsse-convergence-increment 2 --no-intro \
    | sed -n '/Also generated/,/report the shift/p')
echo
echo "atoms_bsse/Ga/structure.fdf box (cutoff 4.0, expect 28.0 Ang):"
sed -n '/LatticeVectors/,/endblock/p' "$OUT/case2-vacuum-buffer/cohesive_runs/atoms_bsse/Ga/structure.fdf" | sed -n '2p'
echo "atoms_bsse_check/Ga/structure.fdf box (cutoff 6.0, expect 32.0 Ang):"
sed -n '/LatticeVectors/,/endblock/p' "$OUT/case2-vacuum-buffer/cohesive_runs/atoms_bsse_check/Ga/structure.fdf" | sed -n '2p'
pause


echo "=================================================================="
echo " Case 3: a full multi-cutoff convergence SCAN (README sections 1.5/2.4)"
echo "=================================================================="
cat <<'EOF'
--bsse-convergence-increment now accepts several values at once, generating
one atoms_bsse_check_<cutoff>/ per point instead of a single before/after
comparison -- the only way to actually SEE whether the correction has
plateaued. This is the feature this session's real GaAs investigation
(README section 4) directly motivated.
EOF
mkdir -p "$OUT/case3-scan"
cp structure.fdf "$OUT/case3-scan/"
echo "\$ stb-cohesive -s structure.fdf -p dojo --bsse-convergence-check --bsse-convergence-increment 2 4 6 --no-intro"
(cd "$OUT/case3-scan" && stb-cohesive -s structure.fdf -p dojo --bsse-convergence-check \
    --bsse-convergence-increment 2 4 6 --no-intro \
    | sed -n '/Also generated 3 BSSE/,/Space group/p')
echo
echo "On disk:"
(cd "$OUT/case3-scan" && find cohesive_runs -maxdepth 1 -type d -name "atoms_bsse_check*" | sort)
pause


echo "=================================================================="
echo " Case 4: multi-site BSSE weighting, live (README sections 1.7/2.5)"
echo "=================================================================="
cat <<'EOF'
structure_multisite.fdf is AB-stacked bilayer graphene: 2 real, physically
different carbon environments (eclipsed/dimer sites vs. non-eclipsed sites)
in the SAME structure. A single flat BSSE reference for "carbon" would
silently average over both -- --bsse-multi-site (default ON) instead builds
one ghost cluster PER site, weighted later by multiplicity.
EOF
mkdir -p "$OUT/case4-multisite"
cp structure_multisite.fdf "$OUT/case4-multisite/"
echo "\$ stb-cohesive -s structure_multisite.fdf -p dojo --no-intro"
(cd "$OUT/case4-multisite" && stb-cohesive -s structure_multisite.fdf -p dojo --no-intro \
    | sed -n '/\[5\] BSSE/,/Space group/p')
echo
echo "Nested per-site folders on disk:"
(cd "$OUT/case4-multisite" && find cohesive_runs/atoms_bsse -maxdepth 2 -type d | sort)
pause


echo "=================================================================="
echo " Case 5: a declared-but-unused ('phantom') species is correctly skipped"
echo "=================================================================="
cat <<'EOF'
A species can be declared in %block ChemicalSpeciesLabel but never actually
placed in the coordinates block (a template left over from a substitution,
say). Since it contributes exactly 0 to the cohesive-energy sum regardless
of its own energy, stb-cohesive correctly skips it -- no wasted DFT
calculation, no incorrect block on the analysis later.
EOF
mkdir -p "$OUT/case5-phantom"
cat > "$OUT/case5-phantom/structure.fdf" <<'EOF'
# GaAs with a phantom, never-placed In species (e.g. left over from a
# planned In-substitution that was never actually done to this geometry).
NumberOfSpecies    3
NumberofAtoms      2

%block ChemicalSpeciesLabel
 1   31   Ga
 2   33   As
 3   49   In
%endblock ChemicalSpeciesLabel

LatticeConstant 1.0 Ang

AtomicCoordinatesFormat  Fractional

%block LatticeVectors
 0.00395996   2.86331622   2.86331617
 2.86331622   0.00395996   2.86331617
 2.86331622   2.86331622   0.00395996
%endblock LatticeVectors

%block AtomicCoordinatesAndAtomicSpecies
  0.00014953   0.00014953   0.00014953   1
  0.24984602   0.24984602   0.24984602   2
%endblock AtomicCoordinatesAndAtomicSpecies
EOF
echo "\$ stb-cohesive -s structure.fdf -p dojo --no-intro"
(cd "$OUT/case5-phantom" && stb-cohesive -s structure.fdf -p dojo --no-intro \
    | sed -n '/\[1\] INPUT STRUCTURE/,/\[2\] DIMENSIONALITY/p')
if [ -d "$OUT/case5-phantom/cohesive_runs/atoms/In" ]; then
    echo "UNEXPECTED: atoms/In/ was generated despite In never being placed"
else
    echo "Confirmed: no atoms/In/ (or atoms_bsse/In/) generated."
fi
pause


echo "=================================================================="
echo " Case 6: CLI vs. the interactive stb-suite menu (4.3.1) -- same result"
echo "=================================================================="
cat <<'EOF'
The interactive menu asks the same questions instead of flags, then calls
the exact same stb-cohesive underneath. Reproducing Case 1's own generation
(default settings, dojo pseudopotentials) through stb-suite -> 4.3.1 and
diffing against it proves the 2 paths are equivalent.
EOF
mkdir -p "$OUT/case6-interactive"
cp structure.fdf "$OUT/case6-interactive/"
# Prompts in order: structure.fdf (no default, must type it), k-density
# (blank -> 0.2), pseudopotential source (1 -> dojo, the suite's 1st
# bundled bank), spin (n), dispersion (n), vacuum (blank -> 20.0), BSSE
# choice (blank -> Y, default), multi-site choice (blank -> Y, default),
# convergence-check choice (n), advanced settings (n -> skip), save report
# (n), then the "Press Enter to continue" pause, then quit.
(cd "$OUT/case6-interactive" && \
    printf '4.3.1\nstructure.fdf\n\n1\nn\nn\n\n\n\nn\nn\nn\n\n0\n' | stb-suite > menu.log 2>&1)
if diff -rq "$OUT/case1-default/cohesive_runs" "$OUT/case6-interactive/cohesive_runs" > /dev/null 2>&1; then
    echo "CLI (Case 1) and interactive menu (4.3.1) produced byte-identical folders."
else
    echo "UNEXPECTED: CLI and interactive menu results differ -- see"
    echo "  diff -rq '$OUT/case1-default/cohesive_runs' '$OUT/case6-interactive/cohesive_runs'"
fi
pause


echo "=================================================================="
echo " Case 7: Stage 2 -- a single BSSE convergence check (stb-cohesiveAnalysis, 4.3.2)"
echo "=================================================================="
cat <<'EOF'
This script doesn't run real SIESTA, so Stage 2 is demonstrated here
against small, hand-built calc.out data -- illustrative numbers only (see
the README's section 4 for this material's real, physically-verified
result). Watch [2]'s Delta columns and [3]'s BSSE correction per atom.
EOF
mkdir -p "$OUT/case7-stage2-single"
cp structure.fdf "$OUT/case7-stage2-single/"
(cd "$OUT/case7-stage2-single" && stb-cohesive -s structure.fdf -p dojo \
    --bsse-convergence-check --bsse-convergence-increment 2 --no-intro > /dev/null)
RUN7="$OUT/case7-stage2-single/cohesive_runs"
write_atom_out "$RUN7/structure" "-500.000000" 8
write_atom_out "$RUN7/atoms/Ga" "-50.000000"
write_atom_out "$RUN7/atoms/As" "-70.000000"
write_atom_out "$RUN7/atoms_bsse/Ga" "-50.150000"
write_atom_out "$RUN7/atoms_bsse/As" "-70.200000"
write_atom_out "$RUN7/atoms_bsse_check/Ga" "-50.220000"
write_atom_out "$RUN7/atoms_bsse_check/As" "-70.260000"
echo "\$ stb-cohesiveAnalysis -o calc.out --no-intro"
(cd "$RUN7" && stb-cohesiveAnalysis -o calc.out -d . --no-intro \
    | sed -n '/\[2\] ENERGY EXTRACTION/,/\[4\] CORRECTION PLOT/p')
[ -f "$RUN7/cohesive_correction.png" ] && echo "(cohesive_correction.png written, 2 panels)"
pause


echo "=================================================================="
echo " Case 8: Stage 2 -- reading a full multi-cutoff scan"
echo "=================================================================="
cat <<'EOF'
Same idea as Case 7, now against Case 3's 3-point scan layout
(atoms_bsse_check_6.0/, _8.0/, _10.0/) -- a fresh, deliberately converging
illustrative series this time, to show what [3b] BSSE CUTOFF CONVERGENCE
SCAN and the plot's 3rd panel look like when the correction IS settling
down (contrast with the real GaAs case in the README's own section 4,
where the shift was still large after the same kind of check).
EOF
mkdir -p "$OUT/case8-stage2-scan"
cp structure.fdf "$OUT/case8-stage2-scan/"
(cd "$OUT/case8-stage2-scan" && stb-cohesive -s structure.fdf -p dojo \
    --bsse-convergence-check --bsse-convergence-increment 2 4 6 --no-intro > /dev/null)
RUN8="$OUT/case8-stage2-scan/cohesive_runs"
write_atom_out "$RUN8/structure" "-500.000000" 8
write_atom_out "$RUN8/atoms/Ga" "-50.000000"
write_atom_out "$RUN8/atoms/As" "-70.000000"
write_atom_out "$RUN8/atoms_bsse/Ga" "-50.100000"
write_atom_out "$RUN8/atoms_bsse/As" "-70.140000"
write_atom_out "$RUN8/atoms_bsse_check_6.0/Ga" "-50.150000"
write_atom_out "$RUN8/atoms_bsse_check_6.0/As" "-70.180000"
write_atom_out "$RUN8/atoms_bsse_check_8.0/Ga" "-50.170000"
write_atom_out "$RUN8/atoms_bsse_check_8.0/As" "-70.195000"
write_atom_out "$RUN8/atoms_bsse_check_10.0/Ga" "-50.175000"
write_atom_out "$RUN8/atoms_bsse_check_10.0/As" "-70.198000"
echo "\$ stb-cohesiveAnalysis -o calc.out --no-intro"
(cd "$RUN8" && stb-cohesiveAnalysis -o calc.out -d . --no-intro \
    | sed -n '/\[3b\] BSSE CUTOFF/,/\[4\] CORRECTION PLOT/p')
[ -f "$RUN8/cohesive_correction.png" ] && echo "(cohesive_correction.png written, 3 panels -- the scan trend)"
pause


echo "=================================================================="
echo " Workflow 4.3 complete"
echo "=================================================================="
cat <<'EOF'
Stage 1 (stb-cohesive, 4.3.1): structure/, atoms/<species>/, and (by
default) atoms_bsse/<species>/[site_.../] written, with a guaranteed
--vacuum/2 buffer regardless of --bsse-cutoff, and an optional single check
or full multi-cutoff scan. Stage 2 (stb-cohesiveAnalysis, 4.3.2): reads real
SIESTA runs back (-d/--dir auto-finds them, default cohesive_runs, matching
Stage 1's own default -- no cd needed), reports uncorrected vs.
BSSE-corrected cohesive energy with SCF/force quality diagnostics, a
convergence-scan table when present, and a matplotlib plot (always saved,
--view for an interactive popup).

See the README's section 4 for this exact material's REAL, physically
-verified result -- including a real bug found and fixed while investigating
why the correction wasn't converging -- and section 6 for a step-by-step
guide to running this on your own structure.
EOF
