#!/bin/bash
# Guided example: OER (Oxygen Evolution Reaction) workflow
# (stb-oer / stb-oerIntermediates / stb-oerRefs / stb-oerAnalysis,
# codes 4.14.1/4.14.2/4.14.3/4.14.4)
#
# Not an automated test (see test/4-workflow/14-oer/{prep,intermediates,
# refs,analysis}/test.sh for that) -- a commented walk-through: it runs
# real commands, one group at a time, and shows you the piece of output
# that proves what just happened. It pauses between sections so you can
# read before moving on.
#
# stb-oer, stb-oerIntermediates and stb-oerRefs are exercised for real
# (they only write input files, never run SIESTA themselves).
# stb-oerAnalysis needs real SIESTA .out/.FA files to analyze, which this
# walkthrough doesn't have -- so the full-chain worked example (Section
# "output/workflow/" below) fabricates calc.out/.FA files with a
# HAND-CHOSEN, analytically known set of energies/forces, so the correct
# eta and potential-determining step (PDS) are known *before* running
# Stage 4, not just plausible after the fact. See the README's Section 8
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
# diagnostics in stb-oerAnalysis read clean by default.
write_freeeng() {
    local path="$1" energy="$2" maxforce="${3:-0.005000}"
    printf 'siesta: SCF Convergence by DM criterion\nSCF cycle converged after 12 iterations\nsiesta: FreeEng =    %s\nsiesta: Atomic forces (eV/Ang):\n     Max    %s\n' \
        "$energy" "$maxforce" > "$path"
}

# Appends an 'outcoor: Relaxed atomic coordinates (fractional)' block to an
# already-written calc.out (see write_freeeng above), simulating a finished
# CG relaxation -- stb-oerIntermediates/stb-oerRefs read coordinates from
# exactly this block, not from a separate .XV file. $frac_lines is the
# exact 'x y z species_id' text expected, one atom per line, in the SAME
# order/count as the folder's own structure.fdf.
append_relaxed_coords() {
    local path="$1" frac_lines="$2"
    printf 'outcoor: Relaxed atomic coordinates (fractional):\n%s\n' "$frac_lines" >> "$path"
}

# Writes every disp_NNN/<system_label>.FA file for a LOCAL-mode ZPE folder:
# an isotropic harmonic spring of stiffness k_by_symbol[symbol] (eV/Ang^2)
# on EACH local atom (independent, decoupled oscillators -- see the
# README's Section 1.4), zero force on every other atom. The exact same
# synthetic-Hessian trick test/4-workflow/14-oer/analysis/test.sh's own
# write_local_zpe_fixture uses, so ZPE/TS come out to a known, hand
# -verifiable closed form for EVERY intermediate (O* alone, OH*, OOH*,
# and H2O), not just a single atom the way HER's own 1-atom H* case needs.
write_local_zpe_fa() {
    local zpe_dir="$1" system_label="$2" natoms="$3" local_indices_csv="$4" local_symbols_csv="$5"
    python3 -c "
import os, json

zpe_dir = '$zpe_dir'
system_label = '$system_label'
n_atoms = $natoms
local_indices = [$local_indices_csv]
local_symbols = [$local_symbols_csv]
k_by_symbol = {'O': 8.0, 'H': 5.0}
d = 0.015

order = []
for atom_index in local_indices:
    for axis in range(3):
        for sign in (1.0, -1.0):
            order.append({'atom_index': atom_index, 'axis': axis, 'sign': sign})

os.makedirs(zpe_dir, exist_ok=True)
with open(os.path.join(zpe_dir, 'zpe_local_meta.json'), 'w') as f:
    json.dump({'local_indices': local_indices, 'local_symbols': local_symbols,
               'displacement_ang': d, 'order': order, 'system_label': system_label}, f)

for i, entry in enumerate(order, start=1):
    atom_index, axis, sign = entry['atom_index'], entry['axis'], entry['sign']
    moved_symbol = local_symbols[local_indices.index(atom_index)]
    k = k_by_symbol[moved_symbol]
    force = -k * (sign * d)
    disp_dir = os.path.join(zpe_dir, f'disp_{i:03d}')
    os.makedirs(disp_dir, exist_ok=True)
    with open(os.path.join(disp_dir, f'{system_label}.FA'), 'w') as f:
        f.write(f'{n_atoms}\n')
        for j in range(n_atoms):
            fxyz = [0.0, 0.0, 0.0]
            if j == atom_index:
                fxyz[axis] = force
            f.write(f'{j+1}   {fxyz[0]: .9E}   {fxyz[1]: .9E}   {fxyz[2]: .9E}\n')
"
}

echo "=================================================================="
echo " Why this workflow needs four stages, and what eta/PDS mean"
echo "=================================================================="
cat <<'EOF'
The Oxygen Evolution Reaction (OER, 2H2O -> O2 + 4H+ + 4e-) is the anodic
half-reaction in water electrolysis, and typically the efficiency
bottleneck (its overpotential is almost always larger than HER's). The
4-electron computational hydrogen electrode (CHE, Rossmeisl et al. 2007;
Man et al. 2011) walks it through THREE adsorbed intermediates:

  * + H2O  -> OH*  + (H+ + e-)      dG1
  OH*      -> O*   + (H+ + e-)      dG2
  O* + H2O -> OOH*  + (H+ + e-)     dG3
  OOH*     -> *  + O2 + (H+ + e-)   dG4

eta = max(dG1..dG4) - 1.23 V; whichever step has the LARGEST dG is the
Potential-Determining Step (PDS) -- the step that needs the extra applied
voltage. A GOOD catalyst has all four steps close to 1.23 V/electron
(eta close to 0); a POOR one has one step that towers over the rest.

This needs, at minimum, SIX consistent SIESTA energies (OH*, O*, OOH*, the
clean slab, H2, H2O) plus a Boys-Bernardi BSSE triad and a ZPE/entropy
Hessian for EACH of the three intermediates AND H2O -- 4x more bookkeeping
than HER's single Delta-G_H*. Four stages split that up exactly the way
HER's three do, plus one: site search (Stage 1) -> derive O*/OOH* from
THAT SAME site (Stage 2, new relative to HER -- Section 1.2 of the README
explains why this has to be the SAME site, not independently re-searched)
-> reference/BSSE/ZPE prep (Stage 3) -> analysis (Stage 4).
EOF
pause

echo "=================================================================="
echo " output/stage1/  --  Stage 1 mechanics: symmetry, site search, config_extra.fdf"
echo "=================================================================="
echo "structure.fdf is the SAME free-standing 2-atom graphene primitive"
echo "cell 4.13-her's own example uses -- small and fast, purely to"
echo "exercise the tools (graphene itself is a chemically POOR OER"
echo "catalyst, see structure.fdf's own header comment)."
echo
echo "\$ stb-oer -s structure.fdf -c calc.fdf -O output/stage1 --no-intro"
stb-oer -s structure.fdf -c calc.fdf -O "$OUT/stage1" --no-intro > "$OUT/stage1_console.log"
echo
echo "[1] SLAB SYMMETRY -- the SAME space/point/layer group HER's own"
echo "example reports (same fixture): P6/mmm, one symmetrically distinct"
echo "carbon site (Wyckoff 'd', multiplicity 2):"
sed -n '/\[1\] SLAB SYMMETRY/,/Note: these are/p' "$OUT/stage1_console.log"
echo
echo "[3] ADSORPTION SITES: raw-vs-symmetry-reduced candidate counts, one"
echo "row per site type PLUS a TOTAL row (--site-type defaults to 'all'):"
sed -n '/Site type | Raw/,/TOTAL/p' "$OUT/stage1_console.log"
echo
echo "4 site folders were written -- 1 distinct ontop site, 3 distinct"
echo "bridge sites, 0 hollow (pymatgen's site-finder never proposes a"
echo "hollow candidate for this 2-atom honeycomb basis, at any cell size):"
ls "$OUT/stage1/sites/" | grep site_
pause

echo "------------------------------------------------------------------"
echo " Fragment labels and config_extra.fdf, live"
echo "------------------------------------------------------------------"
echo "Every atom is labeled by WHICH FRAGMENT it came from -- '<symbol>_slab'"
echo "for the substrate, 'O_ads'/'H_ads' for the adsorbed OH -- exactly"
echo "HER's own convention (Section 1.7), just with TWO adsorbate species"
echo "instead of one:"
grep -A4 "ChemicalSpeciesLabel" "$OUT/stage1/sites/site_1_ontop/structure.fdf"
echo
echo "The forced directives (fixed cell, Slab.DipoleCorrection, Spin"
echo "polarized, DFTD3 -- all MANDATORY here, not opt-in, same reasoning"
echo "as HER's own Section 1.5) live in a config_extra.fdf sidecar,"
echo "%include'd on top of your OWN --calc template rather than editing"
echo "it in place:"
cat "$OUT/stage1/sites/site_1_ontop/config_extra.fdf"
pause

echo "=================================================================="
echo " output/stage1_orient/  --  orientation sampling (new vs. HER)"
echo "=================================================================="
cat <<'EOF'
Unlike H* (a single atom, no orientation), OH* is a 2-atom rod that can
point at the surface many different ways. --n-orientations-polar/-azimuthal
samples a Fibonacci-sphere grid of starting orientations PER SITE (same
machinery stb-adsorb's own orientation sampling uses) -- optionally
MACE-MP-0 ranked via --ml-rank (skipped here: no ML dependency in this
walkthrough, see the README's Section 2), or, as below, every orientation
written UNSCREENED as its own folder:
EOF
echo "\$ stb-oer --site-type ontop --n-orientations-polar 2 --n-orientations-azimuthal 2"
stb-oer -s structure.fdf -c calc.fdf --site-type ontop --n-orientations-polar 2 \
    --n-orientations-azimuthal 2 -O "$OUT/stage1_orient" --no-intro \
    > "$OUT/stage1_orient.log"
echo
echo "2x2 = 4 orientations of the SAME site, 4 separate folders:"
ls "$OUT/stage1_orient/sites/" | grep site_
pause

echo "=================================================================="
echo " output/stage1_both/  --  a free-standing 2D material has two faces"
echo "=================================================================="
echo "\$ stb-oer --site-type ontop --both-sides"
stb-oer -s structure.fdf -c calc.fdf --site-type ontop --both-sides \
    -O "$OUT/stage1_both" --no-intro > "$OUT/stage1_both.log"
grep "NumberofAtoms" "$OUT/stage1_both/sites/site_1_ontop_bothsides/structure.fdf"
echo "(2 C substrate atoms + 1 OH per face = 6 atoms total)"
pause

echo "=================================================================="
echo " output/stage1_position/ and stage1_positionsfile/  --  manual overrides"
echo "=================================================================="
echo "\$ stb-oer --position 0.0 0.0 --height 1.6"
stb-oer -s structure.fdf -c calc.fdf --position 0.0 0.0 --height 1.6 \
    -O "$OUT/stage1_position" --no-intro > "$OUT/stage1_position.log"
ls "$OUT/stage1_position/sites/" | grep site_
echo
echo "Every stb-oer run -- automatic or manual -- always writes"
echo "sites/site_positions.dat, feedable straight back via --positions-file:"
cp "$OUT/stage1_position/sites/site_positions.dat" "$OUT/mypositions.dat"
echo "\$ stb-oer --positions-file mypositions.dat"
stb-oer -s structure.fdf -c calc.fdf --positions-file "$OUT/mypositions.dat" \
    -O "$OUT/stage1_positionsfile" --no-intro > "$OUT/stage1_positionsfile.log"
ls "$OUT/stage1_positionsfile/sites/" | grep site_
echo "Same site, same label -- a full round trip."
pause

echo "=================================================================="
echo " output/stage1_overlap/  --  the too-close-atoms safety check"
echo "=================================================================="
echo "\$ stb-oer --site-type ontop --height 0.1"
stb-oer -s structure.fdf -c calc.fdf --site-type ontop --height 0.1 \
    -O "$OUT/stage1_overlap" --no-intro > "$OUT/stage1_overlap.log"
grep "WARNING" "$OUT/stage1_overlap.log"
echo "(0.1 Ang is far below any real X-OH bond length -- this is a sanity"
echo "check on --height, not a real chemisorption geometry.)"
pause

echo "=================================================================="
echo " output/workflow/  --  the full 4-stage chain, a KNOWN eta/PDS"
echo "=================================================================="
cat <<'EOF'
Real SIESTA output isn't available inside this walkthrough (there's no
SIESTA binary to invoke), so -- exactly like 4.13-her's own worked example
-- this section fabricates calc.out/.FA files whose numbers were chosen
BY HAND (worked BACKWARD from clean, hand-checkable dG1..dG4 values) so
the correct eta/PDS are known in advance:

  dG1 = +2.20 eV   dG2 = +1.10 eV   dG3 = +0.70 eV   dG4 = +0.92 eV
  sum = 4.92 eV exactly (by construction, see the README's Section 1.1)
  eta = max(dG1..dG4) - 1.23 = 2.20 - 1.23 = +0.97 V
  PDS = Step 1 (* + H2O -> OH*)

Three DISTINCT BSSE corrections (Section 1.3 of the README explains why
this MUST be 3 separate triads, never one shared/reused triad):
  BSSE(OH*)  = +0.10 eV (slab +0.06, adsorbate +0.04)
  BSSE(O*)   = +0.15 eV (slab +0.09, adsorbate +0.06)
  BSSE(OOH*) = +0.05 eV (slab +0.02, adsorbate +0.03)

ZPE/TS via isotropic harmonic springs (Section 1.4), k_O = 8.0 eV/Ang^2,
k_H = 5.0 eV/Ang^2 -- an EXACT, hand-verifiable closed form for all three
intermediates AND H2O simultaneously (independently computed, not by
calling stb-oerAnalysis).
EOF
pause

echo "--- Stage 1 (stb-oer): a single ontop site, the winning-site candidate ---"
RUN="$OUT/workflow/oer_study"
echo "\$ stb-oer --site-type ontop -O output/workflow/oer_study"
stb-oer -s structure.fdf -c calc.fdf --site-type ontop -O "$RUN" --no-intro \
    > "$OUT/workflow_stage1.log"
ls "$RUN/sites"
echo
echo "Fabricating this site's 'finished, relaxed' calc.out -- FreeEng plus"
echo "an outcoor block with the OH group pulled in slightly, simulating a"
echo "real CG relaxation:"
write_freeeng "$RUN/sites/site_1_ontop/calc.out" "-851.870886"
append_relaxed_coords "$RUN/sites/site_1_ontop/calc.out" "    0.00000000    0.00000000    0.50000000   1
    0.33333333    0.66666667    0.50000000   1
    0.00100000    0.00100000    0.58000000   2
    0.00100000    0.00100000    0.62000000   3"
grep "FreeEng\|outcoor" "$RUN/sites/site_1_ontop/calc.out"
pause

echo "--- Stage 2 (stb-oerIntermediates): derive O*/OOH* from THAT SAME site ---"
echo "\$ stb-oerIntermediates --directory output/workflow/oer_study"
stb-oerIntermediates --directory "$RUN" --no-intro > "$OUT/workflow_stage2.log"
echo
echo "O* = OH* minus H (3 atoms); OOH* = OH* plus a second O-H group (5"
echo "atoms) -- BOTH derived from the SAME O1 the winning OH* site's own"
echo "relaxation converged to, never a different/independently-searched site:"
grep "atoms) (derived from the winning OH\* site" "$OUT/workflow_stage2.log" 2>/dev/null || true
sed -n '/\[2\] O\* GEOMETRY/,/\[4\] SUMMARY/p' "$OUT/workflow_stage2.log" | grep "\[OK\]"
pause

echo "--- Fabricating O*/OOH*'s own 'relaxed' calc.out (worked-backward values) ---"
python3 - "$RUN" <<'PYEOF'
import sys
run = sys.argv[1]

def write_relaxed(path, freeeng, coords):
    with open(path, "w") as f:
        f.write("siesta: SCF Convergence by DM criterion\n")
        f.write("SCF cycle converged after 12 iterations\n")
        f.write(f"siesta: FreeEng =    {freeeng:.6f}\n")
        f.write("outcoor: Relaxed atomic coordinates (fractional)\n")
        for row in coords:
            f.write(f"    {row[0]:.8f}    {row[1]:.8f}    {row[2]:.8f}   {row[3]}\n")
        f.write("siesta: Atomic forces (eV/Ang):\n")
        f.write("     Max    0.005000\n")

write_relaxed(f"{run}/intermediates/o_star/calc.out", -834.791772, [
    (0.0, 0.0, 0.5, 1), (0.333333333, 0.666666667, 0.5, 1), (0.001, 0.001, 0.58, 2),
])
write_relaxed(f"{run}/intermediates/ooh_star/calc.out", -1287.962658, [
    (0.0, 0.0, 0.5, 1), (0.333333333, 0.666666667, 0.5, 1), (0.001, 0.001, 0.58, 2),
    (0.001, 0.001, 0.66, 2), (0.001, 0.001, 0.70, 3),
])
PYEOF
echo "Done -- O*/OOH* both 'converged'."
pause

echo "--- Stage 3 (stb-oerRefs): references, 3 SEPARATE BSSE triads, local ZPE ---"
echo "\$ stb-oerRefs --directory output/workflow/oer_study"
stb-oerRefs --directory "$RUN" --no-intro > "$OUT/workflow_stage3.log"
echo
echo "[5] BSSE CORRECTION -- one triad PER INTERMEDIATE, each at ITS OWN"
echo "geometry (12 folders total: 3 x 4) -- see the README's Section 1.3"
echo "for why a single shared triad was removed from this suite entirely:"
sed -n '/Folder                           | Atoms/,/05_bsse_OOH_isolated/p' "$OUT/workflow_stage3.log"
echo
echo "[6] ZPE PREPARATION -- local-mode folder counts differ per intermediate"
echo "(O* has only 1 local atom -> 6 folders; OOH*/H2O have 3 -> 18 each):"
sed -n '/Intermediate | Local atom/,/H2O/p' "$OUT/workflow_stage3.log"
pause

echo "--- Fabricating every reference/BSSE/ZPE folder's calc.out/.FA ---"
write_freeeng "$RUN/00_clean_slab/calc.out"    "-400.000000"
write_freeeng "$RUN/02_h2_molecule/calc.out"   "-31.500000"
write_freeeng "$RUN/03_h2o_molecule/calc.out"  "-470.000000"
echo "04_slab_deformed gets a DELIBERATELY large residual force (1.5 eV/Ang)"
echo "-- see the next section for why."
write_freeeng "$RUN/04_slab_deformed/calc.out" "-399.900000" "1.500000"

write_freeeng "$RUN/05_bsse_OH_slab_only/calc.out"            "-400.000000"
write_freeeng "$RUN/05_bsse_OH_slab_ghost/calc.out"           "-400.060000"
write_freeeng "$RUN/05_bsse_OH_adsorbate_ghost_slab/calc.out" "-100.040000"
write_freeeng "$RUN/05_bsse_OH_isolated/calc.out"             "-100.000000"
write_freeeng "$RUN/05_bsse_O_slab_only/calc.out"              "-400.000000"
write_freeeng "$RUN/05_bsse_O_slab_ghost/calc.out"             "-400.090000"
write_freeeng "$RUN/05_bsse_O_adsorbate_ghost_slab/calc.out"   "-80.060000"
write_freeeng "$RUN/05_bsse_O_isolated/calc.out"               "-80.000000"
write_freeeng "$RUN/05_bsse_OOH_slab_only/calc.out"             "-400.000000"
write_freeeng "$RUN/05_bsse_OOH_slab_ghost/calc.out"            "-400.020000"
write_freeeng "$RUN/05_bsse_OOH_adsorbate_ghost_slab/calc.out"  "-120.030000"
write_freeeng "$RUN/05_bsse_OOH_isolated/calc.out"              "-120.000000"

write_local_zpe_fa "$RUN/08_zpe_calc_OH"  "oer_zpe_oh"  4 "2, 3"    "'O', 'H'"
write_local_zpe_fa "$RUN/08_zpe_calc_O"   "oer_zpe_o"   3 "2"       "'O'"
write_local_zpe_fa "$RUN/08_zpe_calc_OOH" "oer_zpe_ooh" 5 "2, 3, 4" "'O', 'O', 'H'"
write_local_zpe_fa "$RUN/08_zpe_calc_H2O" "oer_zpe_h2o" 3 "0, 1, 2" "'O', 'H', 'H'"
echo "Done -- 16 reference/BSSE calc.out + 54 disp_NNN/*.FA fabricated."
pause

echo "--- Stage 4 (stb-oerAnalysis): does the math check out? ---"
echo "\$ stb-oerAnalysis --directory output/workflow/oer_study --save-report --plot --no-show"
stb-oerAnalysis --directory "$RUN" --save-report --plot --no-show --no-intro \
    > "$OUT/workflow_stage4.log"
echo
echo "[2] BSSE CORRECTION: 3 GENUINELY DIFFERENT corrections, exactly the"
echo "hand-chosen values above:"
sed -n '/BSSE (05_bsse_OH):/p;/BSSE (05_bsse_O):/p;/BSSE (05_bsse_OOH):/p' "$OUT/workflow_stage4.log"
echo
echo "[5]/[6]: the four reaction steps and the final result:"
sed -n '/Step | Reaction/,/OOH\* -> \* + O2/p' "$OUT/workflow_stage4.log"
sed -n '/\[6\] FINAL RESULT/,/PDS)/p' "$OUT/workflow_stage4.log" | head -5

GOT_ETA=$(grep "overpotential eta = " "$OUT/workflow_stage4.log" | grep -oE '[-+][0-9.]+' | head -1)
GOT_PDS=$(grep "Potential-determining step (PDS) = Step" "$OUT/workflow_stage4.log" | grep -oE '[0-9]+$')
echo
if [ "$GOT_ETA" = "+0.9700" ] && [ "$GOT_PDS" = "1" ]; then
    echo "Confirmed: eta = $GOT_ETA V, PDS = Step $GOT_PDS -- both match the hand-computed values exactly."
else
    echo "Unexpected: got eta=$GOT_ETA PDS=$GOT_PDS, expected eta=+0.9700 PDS=1"
    exit 1
fi
pause

echo "------------------------------------------------------------------"
echo " A large residual force on 04_slab_deformed: NOT skipped"
echo " (a real difference from HER's own 03_slab_deformed/04_slab_ghost)"
echo "------------------------------------------------------------------"
cat <<'EOF'
HER's own 03_slab_deformed/04_slab_ghost are deliberately skipped by its
own quality check (Section 5.3 of ITS README) -- both are the winning
site's own relaxed geometry with H removed/ghosted, so a large residual
force there is EXPECTED, not a red flag, and printing it every single run
would be pure noise. OER's stb-oerAnalysis does NOT carry that same
skip-list (an honest limitation, not a bug this walkthrough papers over):
04_slab_deformed and every 05_bsse_*/ folder are EQUALLY non-equilibrium
snapshots by construction (all single-point evaluations of an already
-relaxed geometry with atoms removed/ghosted), but NONE of them are
exempted from the check here -- this walkthrough's 12 BSSE folders were
fabricated with small forces so ONLY 04_slab_deformed's deliberately large
one shows below, to isolate the point cleanly; in a REAL run, expect this
same [WARNING] on 05_bsse_*/ folders too, for the identical reason. Read
it with that context -- a large force on any of these specific folders is
expected, not a sign the calculation is broken:
EOF
grep "WARNING" "$OUT/workflow_stage4.log"
pause

echo "=================================================================="
echo " output/workflow/oer_study/plot/  --  all three gnuplot charts"
echo "=================================================================="
echo "stb-oerAnalysis --plot always saves THREE charts together (Section"
echo "5.5 of the README): the free-energy (CHE) diagram (plateau/step"
echo "form), the BSSE components breakdown, and BSSE raw-vs-corrected:"
grep "\[Saved\]" "$OUT/workflow_stage4.log"
if command -v gnuplot > /dev/null 2>&1; then
    ( cd "$RUN/plot" && gnuplot oer_free_energy_diagram.gplot \
                     && gnuplot oer_bsse_correction.gplot \
                     && gnuplot oer_bsse_raw_corrected.gplot )
    echo "Rendered: $RUN/plot/oer_{free_energy_diagram,bsse_correction,bsse_raw_corrected}.pdf"
else
    echo "gnuplot not installed -- the .dat/.gplot pairs are still written, just not rendered here."
fi
echo
echo "The Markdown report (OER_report.md) is ALWAYS written too, embedding"
echo "the same charts as PNGs -- unconditional, independent of --plot/--show:"
ls "$RUN/OER_report.md" "$RUN/plot/"*.png
pause

echo "=================================================================="
echo " Proof: CLI and the interactive stb-suite menu agree (Stage 1)"
echo "=================================================================="
echo "Driving 4.14.1 through the interactive menu (piped input, single"
echo "ontop site, pseudopotentials skipped) and comparing against the"
echo "same direct-CLI case from output/stage1/ above."
TMP="$(mktemp -d)"
cp structure.fdf calc.fdf "$TMP/"
echo
echo "\$ printf '4.14.1\\nstructure.fdf\\ncalc.fdf\\n\\nall\\n\\n\\n\\n\\n\\n\\n\\n0\\n' | stb-suite"
(cd "$TMP" && printf '4.14.1\nstructure.fdf\ncalc.fdf\n\nall\n\n\n\n\n\n\n\n0\n' \
    | stb-suite > menu1.log 2>&1)
python3 -c "
import sys
from stb.core import structure_io
import numpy as np
a = structure_io.read_fdf('$TMP/oer_study/sites/site_1_ontop/structure.fdf')
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
  stage1/                  stage1_orient/            stage1_both/
  stage1_position/         stage1_positionsfile/     stage1_overlap/
  workflow/

output/workflow/oer_study/ has the full chain: clean_slab_source/, sites/,
intermediates/{o_star,ooh_star}/, 00_clean_slab/, 02_h2_molecule/,
03_h2o_molecule/, 04_slab_deformed/, 05_bsse_{OH,O,OOH}_{slab_only,
slab_ghost,adsorbate_ghost_slab,isolated}/ (12 folders), 08_zpe_calc_
{OH,O,OOH,H2O}/disp_NNN/ (54 folders total), oer_stage2.txt, oer_stage3.txt,
oer_stage4.txt, OER_report.txt, OER_report.md, and plot/oer_{free_energy_
diagram,bsse_correction,bsse_raw_corrected}.{dat,gplot,pdf,png}.

As a next step, on your OWN slab/2D structure:
  stb-oer -s <structure.fdf> -c <calc.fdf> --site-type all -O oer_study
  # run SIESTA in every oer_study/sites/site_*/ folder, then:
  stb-oerIntermediates --directory oer_study
  # run SIESTA in intermediates/o_star/ and intermediates/ooh_star/, then:
  stb-oerRefs --directory oer_study
  # run SIESTA in every folder stb-oerRefs just wrote, then:
  stb-oerAnalysis --directory oer_study --save-report
  # read [1]-[3] for any [WARNING] before trusting eta -- see the README's
  # Section 7 for which ones are expected (04_slab_deformed, 05_bsse_*/)
  # and which aren't -- and answer 'y' to the plot prompt (or pass
  # --plot) for the three saved charts.
EOF
