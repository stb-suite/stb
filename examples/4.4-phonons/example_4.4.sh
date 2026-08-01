#!/bin/bash
# Guided example: Workflow 4.4 -- Phonons (both stages: stb-phononsCreate,
# code 4.4.1, and stb-phononsPos, code 4.4.2, in the stb-suite menu).
#
# Not an automated test (see test/4-workflow/4-phonons/{prep,analysis}/
# test.sh for that) -- a commented walk-through: it runs real commands, one
# case at a time, into its own output/<case>/ folder, and shows you the
# piece of output that proves what just happened. Pauses between sections so
# you can read before moving on. Safe to re-run any time -- it always starts
# by wiping its own output/.
#
# structure.fdf/calc.fdf are the REAL, relaxed AlP structure/calc template
# this workflow was actually run against this session (README section 5).
#
# Stage 1 (stb-phononsCreate) runs for real in every case below (cheap --
# it only writes .fdf files, no SCF). Stage 2 (stb-phononsPos) needs real
# force CONSTANTS to say anything physically meaningful -- rather than
# hand-fabricate illustrative numbers (like 4.1-strain/4.2-elastic/
# 4.3-cohesive's own scripts do for their simpler energy/stress numbers),
# this script uses stb-mlphonons (ML Simulations, MACE-MP-0) once, up front,
# to compute REAL force constants for this same AlP structure with no
# SIESTA wait -- so every Stage 2 case below is genuinely physically
# meaningful, not illustrative. This is a tutorial-only shortcut: the real
# DFT-based pipeline this workflow documents never touches MACE (README
# section 2).

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


echo "=================================================================="
echo " Welcome: what a DFT phonon calculation from finite displacements is"
echo "=================================================================="
cat <<'EOF'
Stage 1 (stb-phononsCreate) displaces each symmetry-inequivalent atom a
tiny amount (0.01 Ang default) and writes one single-point SIESTA folder
per displacement -- the finite-difference method. Stage 2 (stb-phononsPos)
reads the resulting forces back, builds the force-constant matrix, and
reports dynamical stability, band structure, DOS/PDOS, and thermal
properties (free energy, entropy, heat capacity).

See this folder's README.md for the full theory (section 1), which
external libraries do the actual work (section 2), and a real, physically
-verified worked example on two real materials -- one clean (AlP), one
genuinely unstable (GaAs) -- in section 5.
EOF
pause


echo "=================================================================="
echo " Case 1: Stage 1, default run -- single-point SCF + symmetry reduction"
echo "=================================================================="
cat <<'EOF'
structure.fdf is real, relaxed AlP (zincblende, 2 atoms/cell). calc.fdf is
the SAME calc.fdf used for this structure's own relaxation (MD.TypeOfRun
CG, MD.VariableCell true) -- stb-phononsCreate forces every displacement
folder to pure single-point SCF regardless ([1]), and detects the correct
zincblende space group F-43m, reducing 96 possible displacements down to
just 2 by symmetry ([3]).
EOF
mkdir -p "$OUT/case1-stage1-default"
cp structure.fdf calc.fdf "$OUT/case1-stage1-default/"
echo "\$ stb-phononsCreate -s structure.fdf -c calc.fdf -p dojo -dim 2 2 2 --no-intro"
(cd "$OUT/case1-stage1-default" && stb-phononsCreate -s structure.fdf -c calc.fdf -p dojo -dim 2 2 2 --no-intro \
    | sed -n '/\[1\] SINGLE-POINT SCF ENFORCEMENT/,/\[4\] DISPLACEMENT FOLDERS/p')
echo
echo "On disk:"
(cd "$OUT/case1-stage1-default" && find phonon_runs -maxdepth 1 | sort)
pause


echo "=================================================================="
echo " Case 2: --symprec too tight reproduces the OLD misdetection bug"
echo "=================================================================="
cat <<'EOF'
This session fixed a real, suite-wide bug: several tools defaulted
--symprec to 1e-3 Ang, tighter than pymatgen's own real default (0.01) and
far tighter than typical DFT-relaxation noise. On this exact AlP structure,
a too-tight tolerance misdetects the true cubic F-43m as the lower
-symmetry trigonal R3m -- doubling the displacements needed AND silently
changing which physical symmetry the whole downstream force-constant
reconstruction relies on. Reproduced live below with an even tighter
--symprec 1e-5 (Phonopy's own raw internal default) for contrast.
EOF
mkdir -p "$OUT/case2-symprec-tight"
cp structure.fdf calc.fdf "$OUT/case2-symprec-tight/"
echo "\$ stb-phononsCreate -s structure.fdf -c calc.fdf -p dojo -dim 2 2 2 --symprec 1e-5 --no-intro"
(cd "$OUT/case2-symprec-tight" && stb-phononsCreate -s structure.fdf -c calc.fdf -p dojo -dim 2 2 2 --symprec 1e-5 --no-intro \
    | sed -n '/\[3\] SYMMETRY REDUCTION/,/\[4\] DISPLACEMENT FOLDERS/p')
echo
echo "(Case 1 above, at the correct default symprec=0.01, detected F-43m with only 2 displacements.)"
pause


echo "=================================================================="
echo " Case 3: the supercell k-grid, auto-suggested vs. explicit"
echo "=================================================================="
cat <<'EOF'
The supercell's own SCF needs a DIFFERENT k-grid than the small reference
unit cell's -- stb-phononsCreate auto-suggests one from a target k-point
DENSITY (--kgrid-density, default 0.2 1/Ang, same convention as
stb-kgrid/stb-strain) scaled to the actual supercell size, written first
into config_extra.fdf so it takes precedence. --kgrid overrides this with
an explicit grid instead.
EOF
mkdir -p "$OUT/case3-kgrid"
cp structure.fdf calc.fdf "$OUT/case3-kgrid/"
echo "\$ stb-phononsCreate -s structure.fdf -c calc.fdf -p dojo -dim 2 2 2 --kgrid-density 0.15 --no-intro"
(cd "$OUT/case3-kgrid" && stb-phononsCreate -s structure.fdf -c calc.fdf -p dojo -dim 2 2 2 --kgrid-density 0.15 --no-intro \
    | grep "Supercell k-grid")
echo "\$ stb-phononsCreate -s structure.fdf -c calc.fdf -p dojo -dim 2 2 2 --kgrid 7 7 7 --no-intro"
rm -rf "$OUT/case3-kgrid/phonon_runs"
(cd "$OUT/case3-kgrid" && stb-phononsCreate -s structure.fdf -c calc.fdf -p dojo -dim 2 2 2 --kgrid 7 7 7 --no-intro \
    | grep "Supercell k-grid")
pause


echo "=================================================================="
echo " Generating REAL force constants for Stage 2 (tutorial shortcut only)"
echo "=================================================================="
cat <<'EOF'
Stage 1 never runs SIESTA -- normally you would now cd into phonon_runs/
disp-001/, disp-002/ and run SIESTA yourself, then hand the results to
Stage 2. To demonstrate Stage 2's real, physically meaningful behavior
without that wait, this script uses stb-mlphonons (ML Simulations, MACE
-MP-0) ONCE here to compute real force constants for this same AlP
structure directly -- no SIESTA needed for this shortcut, but the DFT
-based Stage 1/2 pair itself never touches MACE.

Note: stb-mlphonons builds its own Phonopy object independently (a
different tool/category, untouched by this session's --symprec fix) using
Phonopy's raw, unexposed default tolerance -- so it detects R3m (4
displacements) here rather than the F-43m/2 displacements Stage 1 itself
correctly detects on the same structure above. Harmless for this tutorial
shortcut (more displacements than strictly needed is always safe, just
more MACE evaluations) -- not a discrepancy in the real DFT pipeline.
EOF
mkdir -p "$OUT/ml-force-constants"
cp structure.fdf "$OUT/ml-force-constants/"
echo "\$ stb-mlphonons --file structure.fdf --dim 2 2 2 --mesh 20 20 20 --output-dir phonon_runs --no-intro"
(cd "$OUT/ml-force-constants" && stb-mlphonons --file structure.fdf --dim 2 2 2 --mesh 20 20 20 \
    --output-dir phonon_runs --no-intro 2>/dev/null | sed -n '/\[2\] DISPLACEMENTS/,/\[OK\] Force constants/p')
PHONON_RUNS_ML="$OUT/ml-force-constants/phonon_runs"
pause


echo "=================================================================="
echo " Case 4: Stage 2, dynamical stability + the band-structure zero-line"
echo "=================================================================="
cat <<'EOF'
[2] checks the sampled mesh for negative (imaginary) frequencies -- soft
modes, an unrelaxed reference, or numerical noise. [2b] additionally
auto-detects a high-symmetry q-path (ASE's Cell.bandpath, README section
1.3) and checks it too (a mesh alone can miss a zone-boundary instability).
phonon_bands.png's thin dotted BLUE horizontal line at Frequency=0 is
Phonopy's OWN built-in reference guide (traced to
phonopy.phonon.band_structure.BandPlot.decorate() -- README section 4.3)
-- not something stb-suite draws, not a sign of a problem.
EOF
mkdir -p "$OUT/case4-bands"
cp -r "$PHONON_RUNS_ML" "$OUT/case4-bands/phonon_runs"
echo "\$ stb-phononsPos -dir phonon_runs -m 20 20 20 --bands --no-intro"
(cd "$OUT/case4-bands" && stb-phononsPos -dir phonon_runs -m 20 20 20 --bands --no-intro 2>/dev/null \
    | sed -n '/\[2\] DYNAMICAL STABILITY/,/Saved band data/p')
[ -f "$OUT/case4-bands/phonon_runs/phonon_bands.png" ] && echo "(phonon_bands.png written)"
pause


echo "=================================================================="
echo " Case 5: Stage 2 DOS + --view -- the Debye-fit BLUE curve, explained"
echo "=================================================================="
cat <<'EOF'
[3b] fits a Debye (free-particle, g(w) ~ w^2) model to the low-frequency
DOS to extract a Debye frequency/temperature. With --dos --view, Phonopy's
own phonon.plot_total_dos() overlays that SAME quadratic fit as a second,
BLUE curve on top of the real (red) DOS -- this is the OTHER "strange blue
line" investigated this session, and it IS physically meaningful (unlike
the band-structure one): it shows how well the real DOS matches ideal
Debye behavior at low frequency, up to the fitted Debye frequency where it
cuts off. Real, textbook-correct behavior: the fit tracks the real DOS well
only in the lowest few THz, then diverges once optical branches kick in --
exactly what you'd expect for a 2-atom-basis crystal, not any kind of bug.
--view is a live matplotlib popup (a safe no-op here under MPLBACKEND=Agg).
EOF
mkdir -p "$OUT/case5-dos"
cp -r "$PHONON_RUNS_ML" "$OUT/case5-dos/phonon_runs"
echo "\$ stb-phononsPos -dir phonon_runs -m 20 20 20 --dos --view --no-intro"
(cd "$OUT/case5-dos" && stb-phononsPos -dir phonon_runs -m 20 20 20 --dos --view --no-intro 2>/dev/null \
    | sed -n '/\[3b\] DENSITY OF STATES/,/Debye frequency/p')
pause


echo "=================================================================="
echo " Case 6: Stage 2 --pdos --thermal-displacements -- species/ADP detail"
echo "=================================================================="
cat <<'EOF'
--pdos splits the DOS per chemical species (heavier/more weakly-bonded
species shows a low-frequency bump in real materials). --thermal
-displacements computes each atom's isotropic mean-square displacement
(Debye-Waller U_iso) and writes an ADP CIF for visualization -- both
computed on a denser, symmetry-UNREDUCED mesh with eigenvectors, reusing
the same run rather than a second pass for each.
EOF
mkdir -p "$OUT/case6-pdos-adp"
cp -r "$PHONON_RUNS_ML" "$OUT/case6-pdos-adp/phonon_runs"
echo "\$ stb-phononsPos -dir phonon_runs -m 20 20 20 --pdos --thermal-displacements --no-intro"
(cd "$OUT/case6-pdos-adp" && stb-phononsPos -dir phonon_runs -m 20 20 20 --pdos --thermal-displacements --no-intro 2>/dev/null \
    | sed -n '/\[3b\] DENSITY OF STATES/,/  P U_iso/p')
pause


echo "=================================================================="
echo " Case 7: Stage 2 --symprec -- force-constant symmetrization + q-path"
echo "=================================================================="
cat <<'EOF'
--symprec at Stage 2 controls TWO things: the tolerance Phonopy uses to
symmetrize the loaded force constants (correcting small acoustic-sum-rule
violations), AND the tolerance ASE uses to classify the Bravais lattice for
the auto q-path. IMPORTANT, found and verified this session: a Stage-2
--symprec much TIGHTER than what Stage 1 used to reduce the displacement
set can fail to reconstruct the force constants at all (README section
4.6) -- loosen, never tighten blindly, if in doubt. This structure is
already very close to ideal (MACE-relaxed), so the loosened value below
doesn't change the detected Bravais lattice -- README section 5's real,
noisier SIESTA-relaxed GaAs structure is where this actually matters.
EOF
mkdir -p "$OUT/case7-symprec"
cp -r "$PHONON_RUNS_ML" "$OUT/case7-symprec/phonon_runs"
echo "\$ stb-phononsPos -dir phonon_runs -m 20 20 20 --bands --symprec 0.5 --no-intro"
(cd "$OUT/case7-symprec" && stb-phononsPos -dir phonon_runs -m 20 20 20 --bands --symprec 0.5 --no-intro 2>/dev/null \
    | grep -E "Symmetry precision|Bravais lattice")
pause


echo "=================================================================="
echo " Case 8: --save-gnuplot -- opt-in .dat/.gplot pairs (off by default)"
echo "=================================================================="
cat <<'EOF'
This session made gnuplot .dat/.gplot output OPT-IN (previously written
unconditionally on every run, matching this session's aimd_analysis.py/
effmass.py/coop.py convention) -- matplotlib PNGs (phonon_bands.png,
thermal_properties.png) are unaffected, still always saved.
EOF
mkdir -p "$OUT/case8-gnuplot"
cp -r "$PHONON_RUNS_ML" "$OUT/case8-gnuplot/phonon_runs"
echo "\$ stb-phononsPos -dir phonon_runs -m 20 20 20 --bands --dos --save-gnuplot --no-intro"
(cd "$OUT/case8-gnuplot" && stb-phononsPos -dir phonon_runs -m 20 20 20 --bands --dos --save-gnuplot --no-intro 2>/dev/null \
    | grep "Files")
echo
echo "phonon_plots/ contents:"
(cd "$OUT/case8-gnuplot" && ls phonon_runs/phonon_plots | sort)
pause


echo "=================================================================="
echo " Case 9: the interactive path (stb-suite -> 4.4.1 then 4.4.2)"
echo "=================================================================="
cat <<'EOF'
Both stages are reachable via stb-suite's interactive menu (dotted codes
4.4.1/4.4.2) with the exact same questions as the CLI flags above -- Stage
2's menu also has an 'all' shortcut for its additional-analyses list,
expanding to bands+dos+pdos+thermal (NOT mode-freeze, a different kind of
action -- it writes a structure file, not just a plot/report).
EOF
mkdir -p "$OUT/case9-interactive"
cp structure.fdf calc.fdf "$OUT/case9-interactive/"
echo "\$ stb-suite  (4.4.1, defaults, dojo bank)"
(cd "$OUT/case9-interactive" && printf '4.4.1\n\n\n\n\n\n1\nn\nn\nn\n\n0\n' | stb-suite 2>&1 \
    | sed -n '/\[3\] SYMMETRY REDUCTION/,/\[4\] DISPLACEMENT FOLDERS/p')
cp -r "$PHONON_RUNS_ML" "$OUT/case9-interactive/phonon_runs_ml"
echo "\$ stb-suite  (4.4.2, phonon_runs_ml, 'all' shortcut)"
(cd "$OUT/case9-interactive" && printf '4.4.2\nphonon_runs_ml\n\n20 20 20\n\n\n\nall\nn\nn\nn\nn\n\n0\n' | stb-suite 2>&1 \
    | grep "Selected:")
pause


echo "=================================================================="
echo " Workflow 4.4 complete"
echo "=================================================================="
cat <<'EOF'
Stage 1 (stb-phononsCreate, 4.4.1): forces single-point SCF on every
displacement (config_extra.fdf prepend, shared mechanism with
stb-elasticInputs), detects symmetry to minimize how many displacements are
actually needed, auto-suggests (or accepts an explicit) supercell k-grid.
Stage 2 (stb-phononsPos, 4.4.2): reads the finished force data back (real
SIESTA .FA files or an ML-embedded phonopy_disp.yaml), reports dynamical
stability with an explicit numerical tolerance (not a razor's-edge zero),
band structure, DOS/PDOS, thermal properties, and thermal displacements --
--save-gnuplot and --view are both opt-in.

See the README's section 5 for this exact material's real worked example on
TWO real materials (AlP: clean and stable; GaAs: a genuine, forensically
-traced instability), and section 7 for a step-by-step guide to running
this on your own structure.
EOF
