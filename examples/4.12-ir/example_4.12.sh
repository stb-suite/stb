#!/bin/bash
# Guided example: Workflow 4.12 -- IR Spectrum (all three stages: stb-ir,
# code 4.12.1, stb-irModes, code 4.12.2, and stb-irAnalysis, code 4.12.3,
# in the stb-suite menu).
#
# Not an automated test (see test/4-workflow/12-ir/{prep,modes,analysis}/
# test.sh for that) -- a commented walk-through: it runs real commands, one
# case at a time, into its own output/<case>/ folder, and shows you the
# piece of output that proves what just happened. Pauses between sections
# so you can read before moving on. Safe to re-run any time -- it always
# starts by wiping its own output/.
#
# structure.fdf/calc.fdf are the REAL, relaxed primitive-cell NaCl (rock
# salt, Fm-3m) structure and SIESTA calc template this workflow was
# actually run against (dojo pseudopotentials, DZP/GGA-PBE) -- see
# README.md Section 5 for the full worked-example writeup and literature
# comparison.
#
# Stage 1 (stb-ir) runs for real everywhere below -- cheap, it only writes
# displacement folders, no SCF. Stages 2 and 3 need REAL SIESTA-computed
# forces (phonon_disp/disp-*/siesta.FA) and a REAL Born-effective-charge
# run (born_charge_disp/equilibrium/siesta.BC + calc.out) to say anything
# physically meaningful -- rather than requiring a live SIESTA install (not
# guaranteed on the machine running this tutorial, and NaCl's own
# BornCharge run needs the `-np 1` workaround documented in README Section
# 4.4 for a real known SIESTA parallel-diagonalization bug), this script
# copies real_data/ -- genuine SIESTA output from this exact structure,
# same session -- into the folders Stage 1/2 just generated, verified
# byte-for-byte reproducible (same structure + same -dim/-d/--symprec ->
# identical displaced geometries, so the real forces genuinely apply).
# Anyone WITH a working SIESTA install can instead literally run it in
# each generated folder instead of this copy step -- the exact commands
# are printed at each step below.

set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

# Stage 3's --view calls plt.show() -- MPLBACKEND=Agg makes that a no-op
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
echo " Welcome: what an IR spectrum from finite displacements is"
echo "=================================================================="
cat <<'EOF'
stb-ir/stb-irModes/stb-irAnalysis compute the infrared-absorption spectrum
of a structure from first principles, in three stages:

  1. stb-ir        -- builds a supercell, symmetry-reduces how many atomic
                       displacements are actually needed, writes one SIESTA
                       folder per displacement (finite-difference force
                       constants, the same Phonopy machinery stb-raman and
                       item 4.4's stb-phononsCreate/stb-phononsPos share).
  2. stb-irModes    -- once you've run SIESTA in those folders, this loads
                       the force constants, finds the Gamma-point
                       vibrational modes, classifies which ones are
                       IR-active by symmetry, and writes whichever NEW
                       SIESTA folder(s) the intensity calculation needs
                       (auto-selected by dimensionality -- see below).
  3. stb-irAnalysis -- once THOSE are done, combines everything into
                       per-mode intensities and a synthetic spectrum.

IR intensity needs the DIPOLE MOMENT's derivative along each vibrational
mode (dmu/dQ) -- unlike Raman (a polarizability derivative, second order),
this is first order, but *how* you get dmu/dQ depends entirely on whether
the structure is periodic along a given direction:
  - 0D/1D (molecule/wire): SIESTA prints the total dipole moment directly
    -- a +/-delta displacement pair per mode, central-differenced.
  - 3D bulk (this example, NaCl): a naive dipole moment is ill-defined for
    a fully periodic crystal -- SIESTA's Born-effective-charge/Berry-phase
    machinery (BornCharge T) is used instead, ONE equilibrium run for
    every mode.
  - 2D slab (one vacuum axis): a HYBRID of both, one component per axis.

Full theory (why a periodic solid needs Berry phase at all, the harmonic
approximation, symmetry-reduced displacements, the T1u/degenerate-mode
story, and this workflow's own real, currently-unfixed physics limitation)
is in README.md Sections 1 and 6 -- read those alongside this script.
EOF
pause

echo "=================================================================="
echo " Case 1: Stage 1, default run -- supercell + symmetry reduction"
echo "=================================================================="
mkdir -p "$OUT/case1-stage1"
cp structure.fdf calc.fdf "$OUT/case1-stage1/"
(cd "$OUT/case1-stage1" && stb-ir -s structure.fdf -c calc.fdf -p dojo \
    -dim 3 3 3 --kgrid-density 0.15 --symprec 0.01 --no-intro) \
    | sed -n '/\[1\] INPUT STRUCTURE/,/\[2\] SINGLE-POINT/p'
echo
echo "\$ stb-ir -s structure.fdf -c calc.fdf -p dojo -dim 3 3 3 --kgrid-density 0.15 --symprec 0.01 --no-intro"
echo
echo "Only 2 of 3x3x3x2=54 atoms' worth of possible displacements are"
echo "actually needed (Fm-3m symmetry does the rest) -- on disk:"
find "$OUT/case1-stage1/ir_study/phonon_disp" -maxdepth 1 -mindepth 1 | sort
pause

echo "=================================================================="
echo " Case 2: --symprec -- same tolerance pitfall as stb-raman/stb-phononsCreate"
echo "=================================================================="
cat <<'EOF'
Phonopy's OWN raw symprec default (1e-5) is far tighter than any real
DFT-relaxed structure's residual numerical noise and can silently
misdetect a LOWER symmetry than the crystal actually has -- stb-ir uses
0.01 (pymatgen's own default) instead, same fix already applied to
stb-raman/stb-phononsCreate. This NaCl structure came from stb-crystal
builder + a real relaxation and is numerically very clean, so tightening
--symprec here does NOT change anything (still 2 displacements, still
Fm-3m) -- shown below for honesty, not because it demonstrates the bug.
For a structure where this DOES bite (a real misdetection, forensically
traced), see examples/4.4-phonons/README.md Section 5.2 (GaAs).
EOF
mkdir -p "$OUT/case2-symprec"
cp structure.fdf calc.fdf "$OUT/case2-symprec/"
echo "\$ stb-ir ... --symprec 1e-5   (vs. the default 0.01 used in Case 1)"
(cd "$OUT/case2-symprec" && stb-ir -s structure.fdf -c calc.fdf -p dojo \
    -dim 3 3 3 --kgrid-density 0.15 --symprec 1e-5 --no-intro) \
    | grep "Displacements needed"
pause

echo "=================================================================="
echo " Case 3: the supercell k-grid, auto-suggested vs. explicit"
echo "=================================================================="
cat <<'EOF'
--kgrid-density sets how finely the SUPERCELL's electronic structure is
sampled for the finite-displacement force calculation -- this is what the
phonon frequencies stb-irModes reports are actually built from. Too coarse
and the result can end up far from experiment (verified dramatically on
graphene's own Kohn-anomaly-sensitive G-band in stb-raman's example) --
NaCl is a wide-gap insulator, so it's far less sensitive than a
metal/semimetal, but the flag exists for exactly this reason.
EOF
mkdir -p "$OUT/case3a-auto" "$OUT/case3b-explicit"
cp structure.fdf calc.fdf "$OUT/case3a-auto/"
cp structure.fdf calc.fdf "$OUT/case3b-explicit/"
echo "\$ stb-ir ... --kgrid-density 0.15   (auto-suggests a grid from this target density)"
(cd "$OUT/case3a-auto" && stb-ir -s structure.fdf -c calc.fdf -p dojo \
    -dim 3 3 3 --kgrid-density 0.15 --symprec 0.01 --no-intro) \
    | grep "Supercell k-grid"
echo
echo "\$ stb-ir ... --kgrid 4 4 4   (explicit, overrides --kgrid-density)"
(cd "$OUT/case3b-explicit" && stb-ir -s structure.fdf -c calc.fdf -p dojo \
    -dim 3 3 3 --kgrid 4 4 4 --symprec 0.01 --no-intro) \
    | grep "Supercell k-grid"
pause

echo "=================================================================="
echo " Injecting REAL SIESTA data for Stages 2-3 (tutorial shortcut only)"
echo "=================================================================="
cat <<'EOF'
From here on, this script reuses Case 1's own ir_study/ (dojo, -dim 3 3 3,
--kgrid-density 0.15, --symprec 0.01) and copies real_data/'s two real
siesta.FA force files into disp-001/disp-002 -- genuine SIESTA output from
running this EXACT structure this same session (identical geometry,
verified reproducible above). If you have SIESTA installed, you would
instead literally run, in EACH disp-*/ folder:

  $ mpirun -np 6 siesta calc.fdf --out calc.out

(6 MPI ranks is just what was used originally -- any working SIESTA
install/rank count is fine here, this step doesn't hit the Cholesky bug
mentioned below.)
EOF
mkdir -p "$OUT/workflow"
cp structure.fdf calc.fdf "$OUT/workflow/"
cp -r "$OUT/case1-stage1/ir_study" "$OUT/workflow/"
cp real_data/phonon_disp/disp-001/siesta.FA "$OUT/workflow/ir_study/phonon_disp/disp-001/"
cp real_data/phonon_disp/disp-002/siesta.FA "$OUT/workflow/ir_study/phonon_disp/disp-002/"
echo "[OK] real siesta.FA copied into ir_study/phonon_disp/disp-001/ and disp-002/"
pause

echo "=================================================================="
echo " Case 4: Stage 2 -- Gamma modes, symmetry, and the BULK intensity path"
echo "=================================================================="
(cd "$OUT/workflow" && stb-irModes -dir ir_study -p dojo -c calc.fdf --symprec 0.01 --no-intro) \
    | sed -n '/\[1\] PHONON MODES AT GAMMA/,/\[3\] IR DISPLACEMENT/p'
echo
echo "\$ stb-irModes -dir ir_study -p dojo -c calc.fdf --symprec 0.01 --no-intro"
echo
echo "All 3 non-acoustic Gamma modes share ONE frequency and ONE irrep"
echo "(T1u) -- NaCl's rock-salt structure has only one independent optical"
echo "branch, triply degenerate by cubic symmetry. Note the [LIMITATION]"
echo "block above: this frequency has no LO-TO correction -- see README.md"
echo "Section 6 for the full story (root-caused, not fixed, this session)."
pause

echo "=================================================================="
echo " Injecting REAL Born-effective-charge data for Stage 3"
echo "=================================================================="
cat <<'EOF'
Same shortcut as before: born_charge_disp/equilibrium/ needs a real SCF +
BornCharge run (MD.TypeOfRun FC + BornCharge T + PolarizationGrids) to
produce siesta.BC (the Born effective charges) -- copied in from
real_data/ below. If you have SIESTA installed, you would run, inside
THAT one folder specifically:

  $ mpirun -np 1 siesta calc.fdf --out calc.out

-np 1 is NOT optional here -- SIESTA's polarization/Optical module hits a
real, reproducible "cdiag: Error in Cholesky factorisation" bug under
parallel (-np > 1) diagonalization specifically for BornCharge/Optical
runs (confirmed on this exact calculation, both -np 1 and -np 2 tested;
supercell force runs and dipole-difference runs elsewhere in this same
workflow are NOT affected and can use as many ranks as normal).
EOF
cp real_data/born_charge_equilibrium/siesta.BC "$OUT/workflow/ir_study/born_charge_disp/equilibrium/"
cp real_data/born_charge_equilibrium/calc.out "$OUT/workflow/ir_study/born_charge_disp/equilibrium/"
echo "[OK] real siesta.BC + calc.out copied into ir_study/born_charge_disp/equilibrium/"
pause

echo "=================================================================="
echo " Case 5: Stage 3 -- intensities, spectrum, degenerate-group combining"
echo "=================================================================="
(cd "$OUT/workflow" && stb-irAnalysis -dir ir_study --save-gnuplot --no-intro) \
    | sed -n '/\[2\] IR-ACTIVE MODES SUMMARY/,/\[3c\] MODE VIBRATIONS/p'
echo
echo "\$ stb-irAnalysis -dir ir_study --save-gnuplot --no-intro"
echo
echo "Individual modes 1/2/3 have DIFFERENT dmu/dQ directions (an arbitrary"
echo "Phonopy eigenvector-basis choice within the degenerate T1u subspace),"
echo "but the [Degenerate groups] COMBINED intensity (basis-independent) is"
echo "the physically meaningful number -- 3.317027, one clean peak."
echo
echo "plot/ + mode_animations/ contents:"
find "$OUT/workflow/ir_study/plot" "$OUT/workflow/ir_study/mode_animations" -maxdepth 1 -type f | sort
pause

echo "=================================================================="
echo " Case 6: comparing against the literature TO frequency"
echo "=================================================================="
cat <<'EOF'
NaCl is the textbook Lyddane-Sachs-Teller example -- inelastic-neutron-
scattering literature gives omega_TO ~ 164 cm^-1 / omega_LO ~ 264 cm^-1
(Raunio, Almqvist & Stedman, Phys. Rev. 178, 1496 (1969); cross-checked
via LST against eps(0)=5.9/eps(inf)=2.34, VASP Wiki). This script builds
the same comparison plot used in the real investigation.
EOF
mkdir -p "$OUT/case6-literature"
python3 - "$OUT/workflow/ir_study/plot/ir_spectrum.dat" "$OUT/case6-literature/ir_spectrum_vs_experimental.png" <<'PYEOF'
import sys
import numpy as np
import matplotlib.pyplot as plt

dat_path, out_path = sys.argv[1], sys.argv[2]
SIM_COLOR = "#2a78d6"
EXP_COLOR = "#eb6834"

x, y = np.loadtxt(dat_path, comments="#", unpack=True)
EXPERIMENTAL = [(164.0, "T1u (TO)")]

fig, ax = plt.subplots(figsize=(9, 5.5), dpi=200)
fig.patch.set_facecolor("#fcfcfb")
ax.set_facecolor("#fcfcfb")
ax.plot(x, y, color=SIM_COLOR, linewidth=2.2,
        label="Simulated (SIESTA/dojo, this example)",
        solid_capstyle="round", zorder=3)
for freq, label in EXPERIMENTAL:
    ax.axvline(freq, color=EXP_COLOR, linestyle="--", linewidth=1.6, zorder=2)
    peak_y = float(np.interp(freq, x, y))
    ax.annotate(f"{label}\n{freq:.0f} cm$^{{-1}}$", xy=(freq, peak_y), xytext=(0, 12),
                textcoords="offset points", ha="center", va="bottom",
                color=EXP_COLOR, fontsize=9)
ax.plot([], [], color=EXP_COLOR, linestyle="--", linewidth=1.6,
        label="Experimental (NaCl, Raunio et al. 1969)")
ax.set_xlabel("IR shift (cm$^{-1}$)", fontsize=11, color="#0b0b0b")
ax.set_ylabel("Intensity (arb. units)", fontsize=11, color="#0b0b0b")
ax.set_title("NaCl IR spectrum: simulated vs. experimental",
             fontsize=12, fontweight="bold", color="#0b0b0b", pad=14)
ax.grid(True, color="#e1e0d9", linewidth=0.8, zorder=0)
ax.set_axisbelow(True)
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)
for spine in ["left", "bottom"]:
    ax.spines[spine].set_color("#c3c2b7")
ax.tick_params(colors="#52514e", labelsize=9)
xmin = min(x.min(), min(f for f, _ in EXPERIMENTAL) - 10)
xmax = max(x.max(), max(f for f, _ in EXPERIMENTAL) + 10)
ax.set_xlim(xmin, xmax)
ax.set_ylim(0, max(y) * 1.35)
ax.legend(loc="upper left", frameon=False, fontsize=9.5, labelcolor="#0b0b0b")
fig.tight_layout()
fig.savefig(out_path, facecolor=fig.get_facecolor())
print("saved:", out_path)
PYEOF
echo
echo "Computed T1u : 185.6 cm^-1   |   Literature TO : ~164 cm^-1   (~13% high)"
echo "This ~13% gap is exactly README.md Section 6's LO-TO/NAC limitation,"
echo "not a k-grid/supercell convergence problem (both were checked)."
pause

echo "=================================================================="
echo " Case 7: the interactive path (stb-suite -> 4.12.1)"
echo "=================================================================="
cat <<'EOF'
Stage 1 is reachable via stb-suite's interactive menu (dotted code 4.12.1)
with the exact same questions as the CLI flags above. Stages 2/3 (4.12.2,
4.12.3) follow the identical pattern -- one guided question per CLI flag,
not scripted here since Stage 2/3's prompt sequence depends on how many
modes are found and is easy to explore live yourself: run 'stb-suite',
type 4.12.2, and answer the prompts against this script's own
output/workflow/ir_study/ folder.
EOF
mkdir -p "$OUT/case7-interactive"
cp structure.fdf calc.fdf "$OUT/case7-interactive/"
echo "\$ stb-suite  (4.12.1, dim=3 3 3, dojo bank, --kgrid-density 0.15)"
(cd "$OUT/case7-interactive" && printf '4.12.1\n\n\n3 3 3\n\n1\n\ny\n\n\n\n0.15\n\n0\n' | stb-suite 2>&1 \
    | grep "Displacement folders")
pause

echo "=================================================================="
echo " Workflow 4.12 complete"
echo "=================================================================="
cat <<'EOF'
Full theory (Berry phase, Born effective charges, the harmonic
approximation, symmetry-reduced displacements, degenerate mode groups,
and the LO-TO/non-analytic-correction limitation in detail) is in
README.md. All output from this run is under output/ (gitignored) --
inspect output/workflow/ir_study/plot/ and mode_animations/ directly, or
open output/case6-literature/ir_spectrum_vs_experimental.png.
EOF
