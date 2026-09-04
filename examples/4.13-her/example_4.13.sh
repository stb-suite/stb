#!/bin/bash
# Guided example: HER (Hydrogen Evolution Reaction) workflow
# (stb-her / stb-herRefs / stb-herAnalysis, codes 4.13.1/4.13.2/4.13.3)
#
# Not an automated test (see test/4-workflow/13-her/{prep,refs,analysis}/
# test.sh for that) -- a commented walk-through: it runs real commands, one
# group at a time, and shows you the piece of output that proves what just
# happened. It pauses between sections so you can read before moving on.
#
# stb-her and stb-herRefs are exercised for real (they only write input
# files, never run SIESTA themselves). stb-herAnalysis needs real SIESTA
# .out/.FA files to analyze, which this walkthrough doesn't have -- so the
# full-chain worked example (Section "output/workflow/" below) fabricates
# calc.out/.FA files with a HAND-CHOSEN, analytically known set of
# energies/forces, so the correct Delta-G_H* is known *before* running
# Stage 3, not just plausible after the fact. See the README's Section 6
# for the full arithmetic and why each number was picked.

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

# Writes a synthetic calc.out with a real 'siesta: FreeEng' line, a clean
# SCF-convergence line, and a small residual-force line, so the quality
# diagnostics in stb-herAnalysis read clean by default.
write_freeeng() {
    local path="$1" energy="$2" maxforce="${3:-0.005000}"
    printf 'siesta: SCF Convergence by DM criterion\nSCF cycle converged after 12 iterations\nsiesta: FreeEng =    %s\nsiesta: Atomic forces (eV/Ang):\n     Max    %s\n' \
        "$energy" "$maxforce" > "$path"
}

# Appends an 'outcoor: Relaxed atomic coordinates (fractional)' block to an
# already-written calc.out (see write_freeeng above), simulating a finished
# CG relaxation -- stb-herRefs reads coordinates from exactly this block,
# not from a separate .XV file. $frac_lines is the exact 'x y z species_id'
# text stb-herRefs expects, one atom per line, in the SAME order/count as
# the folder's own structure.fdf.
append_relaxed_coords() {
    local path="$1" frac_lines="$2"
    printf 'outcoor: Relaxed atomic coordinates (fractional):\n%s\n' "$frac_lines" >> "$path"
}

# Writes all 6 disp_NNN/her_zpe.FA files for the local-ZPE Hessian: an
# isotropic harmonic spring of stiffness k (eV/Ang^2) on the ADSORBED H
# (always the LAST atom, index natoms), zero force on every other atom --
# the exact same synthetic-Hessian trick test/4-workflow/13-her/analysis/
# test.sh uses to give ZPE(H*, local) a known, hand-verifiable closed form.
write_local_zpe_fa() {
    local zpe_dir="$1" natoms="$2" k="$3"
    python3 -c "
natoms = $natoms
k = $k
d = 0.015
order = [(0,1.0),(0,-1.0),(1,1.0),(1,-1.0),(2,1.0),(2,-1.0)]
for i, (axis, sign) in enumerate(order, start=1):
    force = [0.0, 0.0, 0.0]
    force[axis] = -sign * k * d
    with open(f'$zpe_dir/disp_{i:03d}/her_zpe.FA', 'w') as f:
        f.write(f'{natoms}\n')
        for atom in range(1, natoms + 1):
            if atom == natoms:
                f.write(f'{atom}   {force[0]: .9E}   {force[1]: .9E}   {force[2]: .9E}\n')
            else:
                f.write(f'{atom}   0.000000000E+00   0.000000000E+00   0.000000000E+00\n')
"
}

echo "=================================================================="
echo " Why this workflow needs three stages, and what Delta-G_H* means"
echo "=================================================================="
cat <<'EOF'
The computational hydrogen electrode (CHE, Norskov et al., J. Phys. Chem. B
108, 17886, 2004) reduces the whole electrochemical Hydrogen Evolution
Reaction to a single number: the Gibbs free energy of adsorbing one H atom
on your surface,

    Delta-G_H* = Delta-E_H (BSSE-corrected) + Delta-ZPE - T*Delta-S

A good HER catalyst needs Delta-G_H* close to ZERO (the Sabatier
principle): too positive and H won't stick in the first place (adsorption
itself is the bottleneck); too negative and H sticks too well to ever
leave as H2 gas (product release is the bottleneck). Getting one honest
number means combining SEVEN separate SIESTA energies (the combined
slab+H system, the bare slab, the H2 gas-phase reference, a BSSE
counterpoise triad, and a partial-Hessian ZPE/entropy correction) with the
right formula, the right sign conventions, and the right level of theory
applied consistently everywhere -- exactly what this 3-stage workflow
automates, and what the rest of this walkthrough proves numerically.

Stage 1 (stb-her) finds every symmetrically distinct H-adsorption site and
writes one relaxation folder per site. You run SIESTA in all of them.
Stage 2 (stb-herRefs) picks the lowest-energy (winning) site, then writes
the H2 gas-phase reference, the BSSE ghost triad, and the ZPE/entropy
calculation folder(s). You run SIESTA in those too. Stage 3
(stb-herAnalysis) combines everything into the final Delta-G_H*.
EOF
pause

echo "=================================================================="
echo " output/stage1/  --  Stage 1 mechanics: symmetry, site search, config_extra.fdf"
echo "=================================================================="
echo "structure.fdf is a bare 2-atom graphene primitive cell (free-standing,"
echo "vacuum along c, both faces equally exposed) -- small and fast, purely"
echo "to exercise the tools."
echo
echo "\$ stb-her -s structure.fdf -c calc.fdf -O output/stage1 --no-intro"
stb-her -s structure.fdf -c calc.fdf -O "$OUT/stage1" --no-intro > "$OUT/stage1_console.log"
echo
echo "[1] SLAB SYMMETRY reports the actual space/point group and layer group"
echo "of your slab -- P6/mmm for pristine graphene, one symmetrically"
echo "distinct carbon site (Wyckoff 'd', multiplicity 2):"
sed -n '/\[1\] SLAB SYMMETRY/,/\[2\] CLEAN SLAB/p' "$OUT/stage1_console.log" | head -n -1
echo
echo "[3] ADSORPTION SITES reports raw-vs-symmetry-reduced candidate counts"
echo "per site type, and (since --site-type defaults to 'all') a TOTAL row:"
sed -n '/Site type | Raw candidates/,/TOTAL/p' "$OUT/stage1_console.log"
echo
echo "4 site folders were written -- graphene's honeycomb lattice gives 1"
echo "distinct ontop site and 3 distinct bridge sites, 0 hollow (pymatgen's"
echo "site-finder never proposes a hollow candidate for this 2-atom basis):"
ls "$OUT/stage1/sites"
pause

echo "------------------------------------------------------------------"
echo " Fragment labels and config_extra.fdf, live"
echo "------------------------------------------------------------------"
echo "Every atom is labeled by WHICH FRAGMENT it came from -- '<symbol>_slab'"
echo "for the substrate, '<symbol>_ads' for the adsorbed H -- so a slab that"
echo "already has its own H (e.g. a passivated edge) can never be confused"
echo "with the freshly-adsorbed one after a write_fdf/read_fdf round trip:"
grep -A3 "ChemicalSpeciesLabel" "$OUT/stage1/sites/site_1_ontop/structure.fdf"
echo
echo "The forced directives (fixed cell, Slab.DipoleCorrection, Spin"
echo "polarized, DFTD3 -- all MANDATORY here, not opt-in) live in a"
echo "config_extra.fdf sidecar, %include'd on top of your OWN --calc"
echo "template rather than editing it in place:"
cat "$OUT/stage1/sites/site_1_ontop/config_extra.fdf"
echo
echo "Slab.DipoleCorrection is structurally required, not just advisory:"
echo "adsorbing H on only ONE face breaks whatever mirror symmetry the"
echo "clean slab had, giving the cell a net dipole along a PERIODIC"
echo "direction -- without the correction, that spurious field contaminates"
echo "every site's energy, and therefore the site ranking Stage 2 does."
pause

echo "=================================================================="
echo " output/stage1_both/  --  a free-standing 2D material has two faces"
echo "=================================================================="
echo "\$ stb-her --site-type ontop --both-sides"
stb-her -s structure.fdf -c calc.fdf --site-type ontop --both-sides \
    -O "$OUT/stage1_both" --no-intro > "$OUT/stage1_both.log"
grep "NumberofAtoms" "$OUT/stage1_both/sites/site_1_ontop_bothsides/structure.fdf"
echo "(2 C substrate atoms + 1 H per face = 4 atoms total)"
pause

echo "=================================================================="
echo " output/stage1_position/  --  manual site override with --position"
echo "=================================================================="
echo "For a site pymatgen's automatic finder doesn't propose (or just to"
echo "reproduce one exact point), --position X Y bypasses site-finding"
echo "entirely -- the height is measured along the TRUE surface normal, not"
echo "a naive Cartesian z-offset (matters for a tilted/non-orthogonal cell):"
echo
echo "\$ stb-her --position 0.0 0.0 --height 1.6"
stb-her -s structure.fdf -c calc.fdf --position 0.0 0.0 --height 1.6 \
    -O "$OUT/stage1_position" --no-intro > "$OUT/stage1_position.log"
grep "Manual override (--position)" "$OUT/stage1_position.log"
ls "$OUT/stage1_position/sites"
pause

echo "=================================================================="
echo " output/stage1_positionsfile/  --  round-tripping site_positions.dat"
echo "=================================================================="
echo "Every stb-her run -- automatic or manual -- always writes"
echo "sites/site_positions.dat: the exact fractional (a, b) positions used,"
echo "editable and feedable straight back via --positions-file. Reusing"
echo "output/stage1's own file (4 sites, --site-type all):"
cp "$OUT/stage1/sites/site_positions.dat" "$OUT/mypositions.dat"
echo
echo "\$ stb-her --positions-file mypositions.dat"
stb-her -s structure.fdf -c calc.fdf --positions-file "$OUT/mypositions.dat" \
    -O "$OUT/stage1_positionsfile" --no-intro > "$OUT/stage1_positionsfile.log"
ls "$OUT/stage1_positionsfile/sites"
echo "Same 4 sites, same labels -- a full round trip."
pause

echo "=================================================================="
echo " output/stage1_overlap/  --  the too-close-atoms safety check"
echo "=================================================================="
echo "\$ stb-her --site-type ontop --height 0.1"
stb-her -s structure.fdf -c calc.fdf --site-type ontop --height 0.1 \
    -O "$OUT/stage1_overlap" --no-intro > "$OUT/stage1_overlap.log"
grep "WARNING.*overlapping atoms" "$OUT/stage1_overlap.log"
echo "(0.1 Ang is far below any real X-H bond length -- this is a sanity"
echo "check on --height, not a real chemisorption geometry.)"
pause

echo "=================================================================="
echo " output/workflow/  --  the full 3-stage chain, a KNOWN Delta-G_H*"
echo "=================================================================="
cat <<'EOF'
Real SIESTA output isn't available inside this walkthrough (there's no
SIESTA binary to invoke), so -- exactly like stb-hubbarduAnalysis's own
linear-response verification (example 4.7) -- this section fabricates
calc.out/.FA files whose numbers were chosen BY HAND so the correct
Delta-G_H* is known in advance:

  E_clean (bare slab)        = -400.000000 eV
  E_slab+H (winning site)    = -416.250000 eV   ->  Delta-E_H (raw) = -0.5000 eV
  E_H2 (gas-phase reference) =  -31.500000 eV
  BSSE (slab side + H side)  =   +0.100000 eV   ->  Delta-E_H (corrected) = -0.4000 eV
  H* local Hessian           = isotropic harmonic spring, k = 5.0 eV/Ang^2

k = 5.0 eV/Ang^2 gives an EXACT, hand-verifiable ZPE/entropy via the
standard harmonic-oscillator formulas (Section 5 of the README) --
independently computed (not by calling stb-herAnalysis) as:

  Delta-ZPE  = +0.0810 eV
  Delta-TS   = -0.1981 eV
  Delta-G_H* = -0.4000 + 0.0810 - (-0.1981) = -0.1209 eV  (near-optimal)
EOF
pause

echo "--- Stage 1 (stb-her): a single ontop site, the winning-site candidate ---"
RUN="$OUT/workflow/her_study"
echo "\$ stb-her --site-type ontop -O output/workflow/her_study"
stb-her -s structure.fdf -c calc.fdf --site-type ontop -O "$RUN" --no-intro \
    > "$OUT/workflow_stage1.log"
ls "$RUN/sites"
echo
echo "Fabricating this site's 'finished, relaxed' calc.out -- FreeEng plus"
echo "an outcoor block with H pulled in slightly (0.575 -> 0.570 frac z,"
echo "i.e. height 1.5 -> 1.4 Ang, simulating a real CG relaxation):"
write_freeeng "$RUN/sites/site_1_ontop/calc.out" "-416.250000"
append_relaxed_coords "$RUN/sites/site_1_ontop/calc.out" "    0.00000000    0.00000000    0.50000000   1
    0.33333333    0.66666667    0.50000000   1
    0.00000000    0.00000000    0.57000000   2"
grep "FreeEng\|outcoor" "$RUN/sites/site_1_ontop/calc.out"
pause

echo "--- Stage 2 (stb-herRefs): winning-site scan, references, local ZPE prep ---"
echo "\$ stb-herRefs --directory output/workflow/her_study --zpe-mode local"
stb-herRefs --directory "$RUN" --zpe-mode local --no-intro > "$OUT/workflow_stage2.log"
echo
echo "[2] REFERENCE FOLDERS ends with a table of every folder written --"
echo "atom count, formula (ghost atoms called out explicitly), run type:"
sed -n '/Folder           | Atoms/,/07_h_isolated/p' "$OUT/workflow_stage2.log"
echo
echo "H2's own config_extra.fdf is the ONE folder forced spin-UNpolarized"
echo "instead of polarized (Section 4.4 of the README) -- H2's ground state"
echo "is an unambiguous closed-shell singlet, unlike the site itself:"
grep "Spin" "$RUN/02_h2_molecule/config_extra.fdf"
pause

echo "--- Fabricating every reference folder's calc.out/.FA ---"
write_freeeng "$RUN/00_clean_slab/calc.out"     "-400.000000"
write_freeeng "$RUN/02_h2_molecule/calc.out"    "-31.500000"
write_freeeng "$RUN/03_slab_deformed/calc.out"  "-399.900000"
write_freeeng "$RUN/04_slab_ghost/calc.out"     "-399.950000"
write_freeeng "$RUN/06_h_ghost_slab/calc.out"   "-100.050000"
write_freeeng "$RUN/07_h_isolated/calc.out"     "-100.000000"
write_local_zpe_fa "$RUN/05_zpe_calc" 3 5.0
echo "Done -- 6 reference calc.out + 6 disp_NNN/her_zpe.FA fabricated."
pause

echo "--- Stage 3 (stb-herAnalysis): does the math check out? ---"
echo "\$ stb-herAnalysis --directory output/workflow/her_study --save-report --no-plot --no-show"
stb-herAnalysis --directory "$RUN" --save-report --no-plot --no-show --no-intro \
    > "$OUT/workflow_stage3.log"
echo
echo "[1]/[2]: the energy table and the clearly-separated BSSE breakdown:"
sed -n '/Term       | Energy/,/E_H_iso/p' "$OUT/workflow_stage3.log"
sed -n '/\[2\] BSSE CORRECTION/,/Delta-E_H (corrected)/p' "$OUT/workflow_stage3.log" | tail -n +3
echo
echo "[4]: the final result, with every contribution itemized top to bottom:"
sed -n '/\[4\] FINAL RESULT/,/Delta-G_H\* (TOTAL)/p' "$OUT/workflow_stage3.log"

GOT=$(grep "Delta-G_H\* = " "$OUT/workflow_stage3.log" | grep -oE '[-+][0-9.]+' | head -1)
echo
if [ "$GOT" = "-0.1209" ]; then
    echo "Confirmed: Delta-G_H* = $GOT eV matches the hand-computed value exactly."
else
    echo "Unexpected: got Delta-G_H* = $GOT, expected -0.1209"
    exit 1
fi
pause

echo "------------------------------------------------------------------"
echo " A large residual force on 03_slab_deformed/04_slab_ghost: no warning"
echo "------------------------------------------------------------------"
cat <<'EOF'
03_slab_deformed and 04_slab_ghost deliberately evaluate the WINNING
SITE's own relaxed geometry with H removed/ghosted -- the atoms nearest
the former H position are never at their own equilibrium there BY
CONSTRUCTION, so a residual force well above --force-tolerance is EXPECTED,
not a red flag. stb-herAnalysis skips the force-quality check entirely for
exactly these two folders (never even a softened note) -- injecting an
enormous, obviously-wrong 2.5 eV/Ang force proves it stays silent:
EOF
sed -i 's/Max    0.005000/Max    2.500000/' "$RUN/03_slab_deformed/calc.out" "$RUN/04_slab_ghost/calc.out"
stb-herAnalysis --directory "$RUN" --no-plot --no-show --no-intro > "$OUT/workflow_largeforce.log" 2>&1
if grep -q "Residual force on E_deformed\|Residual force on E_ghost" "$OUT/workflow_largeforce.log"; then
    echo "Unexpected: a force warning/note was printed for 03/04 after all."
    exit 1
else
    echo "Confirmed: no [WARNING]/[NOTE] at all for E_deformed/E_ghost, even at 2.5 eV/Ang."
fi
echo
echo "The SAME check still fires normally for a folder that ISN'T a"
echo "deliberate non-equilibrium snapshot (E_clean):"
sed -i 's/Max    0.005000/Max    2.500000/' "$RUN/00_clean_slab/calc.out"
stb-herAnalysis --directory "$RUN" --no-plot --no-show --no-intro 2>&1 | grep "WARNING.*E_clean"
sed -i 's/Max    2.500000/Max    0.005000/' "$RUN/00_clean_slab/calc.out" \
    "$RUN/03_slab_deformed/calc.out" "$RUN/04_slab_ghost/calc.out"
pause

echo "=================================================================="
echo " output/workflow/her_study/plot/  --  the gnuplot energy-breakdown chart"
echo "=================================================================="
echo "\$ stb-herAnalysis --directory output/workflow/her_study --plot --no-show"
stb-herAnalysis --directory "$RUN" --plot --no-show --no-intro > "$OUT/workflow_plot.log" 2>&1
grep "\[Saved\]" "$OUT/workflow_plot.log"
if command -v gnuplot > /dev/null 2>&1; then
    ( cd "$RUN/plot" && gnuplot her_energy_breakdown.gplot )
    echo "Rendered: $RUN/plot/her_energy_breakdown.pdf"
else
    echo "gnuplot not installed -- the .dat/.gplot pair is still written, just not rendered here."
fi
pause

echo "=================================================================="
echo " Proof: CLI and the interactive stb-suite menu agree (Stage 1)"
echo "=================================================================="
echo "Driving 4.13.1 through the interactive menu (piped input, single ontop"
echo "site, pseudopotentials skipped) and comparing against the same"
echo "direct-CLI case from output/stage1/ above."
TMP="$(mktemp -d)"
cp structure.fdf calc.fdf "$TMP/"
echo
echo "\$ printf '4.13.1\\nstructure.fdf\\ncalc.fdf\\n\\nontop\\n\\n\\n\\n\\n\\n0\\n' | stb-suite"
(cd "$TMP" && printf '4.13.1\nstructure.fdf\ncalc.fdf\n\nontop\n\n\n\n\n\n0\n' \
    | stb-suite > menu1.log 2>&1)
python3 -c "
import sys
from stb.core import structure_io
import numpy as np
a = structure_io.read_fdf('$TMP/her_study/sites/site_1_ontop/structure.fdf')
b = structure_io.read_fdf('$OUT/stage1/sites/site_1_ontop/structure.fdf')
a_frac = sorted([tuple(np.round(pos, 6)) for _s, pos in a.atoms])
b_frac = sorted([tuple(np.round(pos, 6)) for _s, pos in b.atoms])
sys.exit(0 if a_frac == b_frac else 1)
"
if [ $? -eq 0 ]; then
    echo "Confirmed: the interactive-menu site and the direct-CLI site have the same geometry."
else
    echo "Unexpected: interactive-menu and direct-CLI results differ."
    exit 1
fi
rm -rf "$TMP"
pause

echo "=================================================================="
echo " Done"
echo "=================================================================="
cat <<EOF
Folders generated under output/:
  stage1/                  stage1_both/              stage1_position/
  stage1_positionsfile/    stage1_overlap/           workflow/

output/workflow/her_study/ has the full chain: clean_slab_source/, sites/,
00_clean_slab/, 02_h2_molecule/, 03_slab_deformed/, 04_slab_ghost/,
06_h_ghost_slab/, 07_h_isolated/, 05_zpe_calc/disp_00{1..6}/, her_stage2.txt,
her_stage3.txt, HER_report.txt, and plot/her_energy_breakdown.{dat,gplot,pdf}.

As a next step, on your OWN slab/2D structure:
  stb-her -s <structure.fdf> -c <calc.fdf> --site-type all -O her_study
  # run SIESTA in every her_study/sites/site_*/ folder, then:
  stb-herRefs --directory her_study --zpe-mode local
  # run SIESTA in every folder stb-herRefs just wrote, then:
  stb-herAnalysis --directory her_study --save-report
  # read [2]/[3] for any [WARNING] before trusting the final number,
  # and answer 'y' to the plot prompt for a saved energy-breakdown chart.
EOF
