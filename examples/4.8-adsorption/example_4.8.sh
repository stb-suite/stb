#!/bin/bash
# Guided example: Adsorption workflow
# (stb-adsorb / stb-adsorbBsse / stb-adsorbAnalysis, codes 4.8.1/4.8.2/4.8.3)
#
# Not an automated test (see test/4-workflow/8-adsorption/{prep,bsse,analysis}/
# test.sh for that) -- a commented walk-through: it runs real commands, one
# group at a time, and shows you the piece of output that proves what just
# happened. It pauses between sections so you can read before moving on.
#
# stb-adsorb and stb-adsorbBsse are exercised for real (they only write
# input files, never run SIESTA themselves). stb-adsorbAnalysis needs real
# SIESTA .out files to analyze, which this walkthrough doesn't have -- so
# the full-chain worked example (Section "output/workflow/" below)
# fabricates calc.out files with a hand-chosen set of energies designed to
# make one very concrete point: the BSSE correction can, and here does,
# CHANGE which site you'd conclude is the most stable one. See the
# README's Section 4 for the full arithmetic and why each number was
# picked.

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

# Fabricates a 'relaxed' siesta.XV for a site folder: reads its own
# structure.fdf, shifts the LAST atom (the adsorbate -- stb-adsorb always
# appends it after the substrate) 0.5 Ang closer to the substrate along z,
# and writes it via sisl -- the same recipe
# test/4-workflow/8-adsorption/bsse/test.sh uses to simulate a real SIESTA
# relaxation without actually running SIESTA.
fabricate_relaxed_xv() {
    local site_dir="$1"
    python3 -c "
import sisl
from stb.core import structure_io
fdf = structure_io.read_fdf('$site_dir/structure.fdf')
pmg = structure_io.to_pymatgen(fdf)
cart = pmg.cart_coords.copy()
cart[-1, 2] -= 0.5
atoms = [sisl.Atom(str(s.specie)) for s in pmg]
geom = sisl.Geometry(cart, atoms=atoms, lattice=sisl.Lattice(pmg.lattice.matrix))
sisl.get_sile('$site_dir/siesta.XV', mode='w').write_geometry(geom)
"
}

# Writes a synthetic calc.out with a real 'siesta: FreeEng' line plus a
# clean SCF-convergence line and a small residual-force line, so the
# quality diagnostics in stb-adsorbAnalysis/stb-adsorbBsse read clean.
write_freeeng() {
    local path="$1" energy="$2" maxforce="${3:-0.020000}"
    printf 'siesta: FreeEng =    %s\nSCF cycle converged after 12 iterations\nsiesta: Atomic forces (eV/Ang):\n   Max    %s\n' \
        "$energy" "$maxforce" > "$path"
}

echo "=================================================================="
echo " Why this workflow needs three stages, not two"
echo "=================================================================="
cat <<'EOF'
A real DFT adsorption energy compares three numbers -- the combined
slab+adsorbate system, the bare slab, and the isolated adsorbate:

    E_ads = E_site - E_clean_slab - E_adsorbate

Because SIESTA is a localized-basis (LCAO) code, that comparison is
contaminated by Basis Set Superposition Error (BSSE): in the COMBINED
calculation, each fragment "borrows" extra basis functions that belong to
the other fragment, artificially lowering the combined energy and making
every adsorption energy look more strongly bound than it really is. The
standard fix (Boys & Bernardi, 1970) re-evaluates each fragment WITH the
other fragment's basis functions present as chargeless "ghost" atoms -- but
only at the site's own RELAXED geometry, which doesn't exist until SIESTA
has actually finished relaxing it.

Stage 1 (stb-adsorb) prepares clean_slab/, adsorbate/, and every candidate
sites/site_*/. You run SIESTA in all of those. Stage 2 (stb-adsorbBsse)
then reads each site's finished, relaxed geometry and writes its ghost
-fragment BSSE folders. You run SIESTA in those too. Stage 3
(stb-adsorbAnalysis) combines everything into E_ads, with and without the
BSSE correction.
EOF
pause

echo "=================================================================="
echo " output/stage1/  --  Stage 1 mechanics, single O atom, default ontop site"
echo "=================================================================="
echo "structure.fdf is a bare 2-atom graphene primitive cell (free-standing,"
echo "vacuum along c) -- small and fast, purely to exercise the tools."
echo
echo "\$ stb-adsorb --adsorbate O -O output/stage1 --no-intro"
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate O -O "$OUT/stage1" --no-intro \
    > "$OUT/stage1_console.log"
echo
echo "The cell is forced fixed everywhere (config_extra.fdf, %include-d ahead"
echo "of your own calc.fdf template):"
cat "$OUT/stage1/clean_slab/config_extra.fdf"
echo
echo "clean_slab/, adsorbate/, and one candidate site were written:"
ls "$OUT/stage1"
ls "$OUT/stage1/sites"
echo
echo "Single-atom adsorbates are forced spin-polarized in their OWN isolated"
echo "reference (many, like O, have a non-zero ground-state spin), while the"
echo "combined slab+adsorbate calc.fdf is left exactly as you gave it:"
grep "single atom" "$OUT/stage1_console.log"
pause

echo "=================================================================="
echo " output/stage1_lateral_warning/  --  Section 3.2's in-plane self-interaction check"
echo "=================================================================="
echo "The SAME bare 2-atom cell is also deliberately too small laterally:"
echo "the adsorbate would sit only 2.46 Ang from its own periodic image in"
echo "neighboring cells -- far below the ~8 Ang rule of thumb for avoiding"
echo "lateral adsorbate-adsorbate interaction."
echo
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate O -O "$OUT/stage1_lateral_warning" --no-intro \
    > "$OUT/stage1_lateral_warning.log"
grep "WARNING.*periodic image in the ab-plane" "$OUT/stage1_lateral_warning.log"
pause

echo "=================================================================="
echo " output/stage1_fixed/  --  fixing it with stb-supercell"
echo "=================================================================="
echo "\$ stb-supercell -f structure.fdf -d 4 4 1 -o structure_super.fdf"
mkdir -p "$OUT/stage1_fixed"
(cd "$OUT/stage1_fixed" && stb-supercell -f "$DIR/structure.fdf" -d 4 4 1 -o structure_super.fdf --no-intro \
    > "$OUT/stage1_fixed_supercell.log")
grep "NumberofAtoms" "$OUT/stage1_fixed/structure_super.fdf"
echo
echo "\$ stb-adsorb -s structure_super.fdf --adsorbate O --site-type all --all-sites"
stb-adsorb -s "$OUT/stage1_fixed/structure_super.fdf" -c calc.fdf --adsorbate O \
    --site-type all --all-sites -O "$OUT/stage1_fixed/adsorption_run" --no-intro \
    > "$OUT/stage1_fixed_adsorb.log"
echo
if grep -q "WARNING.*periodic image in the ab-plane" "$OUT/stage1_fixed_adsorb.log"; then
    echo "Unexpected: the lateral-interaction warning is STILL present."
    exit 1
else
    echo "Confirmed: on the bigger 4x4 in-plane supercell, the warning is GONE."
fi
echo
echo "The candidate-site table on this bigger cell (Section [2] of the report):"
sed -n '/Site type | Raw candidates/,/TOTAL/p' "$OUT/stage1_fixed_adsorb.log"
echo "(this particular graphene basis never exposes a 'hollow' candidate to"
echo "pymatgen's site-finder, at any cell size -- ontop and bridge are what"
echo "this walkthrough's worked example below compares.)"
ls "$OUT/stage1_fixed/adsorption_run/sites"
pause

echo "=================================================================="
echo " output/stage1_bothsides/  --  a free-standing 2D material has two faces"
echo "=================================================================="
echo "\$ stb-adsorb --adsorbate H --site-type ontop --both-sides"
stb-adsorb -s structure.fdf -c calc.fdf --adsorbate H --site-type ontop --both-sides \
    -O "$OUT/stage1_bothsides" --no-intro > "$OUT/stage1_bothsides.log"
grep "NumberofAtoms" "$OUT/stage1_bothsides/sites/site_1_ontop_bothsides/structure.fdf"
echo "(2 C substrate atoms + 1 H per face = 4 atoms total)"
pause

echo "=================================================================="
echo " output/workflow/  --  the full 3-stage chain: watching BSSE flip the ranking"
echo "=================================================================="
cat <<'EOF'
This is the point of the whole workflow. Two candidate sites on the fixed
(4x4 supercell) structure -- ontop and bridge -- are given ENERGIES BY HAND
below (this walkthrough has no real SIESTA binary to call), chosen so that:

  * WITHOUT the BSSE correction, 'bridge' looks more stable.
  * WITH the BSSE correction, 'ontop' is actually more stable.

Exactly the failure mode Section 2 of the README warns about: trusting the
uncorrected number would make you carry the WRONG site forward to
production.
EOF
mkdir -p "$OUT/workflow"
RUN="$OUT/workflow/adsorption_run"

echo "--- Stage 1 (stb-adsorb): 4 candidate sites on the fixed supercell ---"
echo "\$ stb-adsorb -s structure_super.fdf --adsorbate O --site-type all --all-sites"
stb-adsorb -s "$OUT/stage1_fixed/structure_super.fdf" -c calc.fdf --adsorbate O \
    --site-type all --all-sites -O "$RUN" --no-intro > "$OUT/workflow_stage1.log"
ls "$RUN/sites"

echo
echo "Only 2 of the 4 candidate sites (site_1_ontop, site_2_bridge) will be"
echo "'finished' below -- site_3_bridge/site_4_bridge are left untouched on"
echo "purpose, to show Stage 2/3 skip an unfinished site instead of failing."
echo
echo "Placeholder pseudopotentials (stand in for what -p/--pseudo-dir would"
echo "have copied -- stb-adsorbBsse reuses these directly):"
for site in site_1_ontop site_2_bridge; do
    echo "fake C pseudopotential" > "$RUN/sites/$site/C.psml"
    echo "fake O pseudopotential" > "$RUN/sites/$site/O.psml"
done
echo "Fabricated 'relaxed' siesta.XV (adsorbate shifted 0.5 Ang closer to the"
echo "substrate than stb-adsorb's own initial guess -- same idea as a real"
echo "SIESTA relaxation would produce for a real chemisorption bond):"
fabricate_relaxed_xv "$RUN/sites/site_1_ontop"
fabricate_relaxed_xv "$RUN/sites/site_2_bridge"
ls "$RUN/sites/site_1_ontop/siesta.XV" "$RUN/sites/site_2_bridge/siesta.XV"

echo
echo "Fabricated calc.out for clean_slab/, adsorbate/, and the 2 finished sites:"
write_freeeng "$RUN/clean_slab/calc.out"           "-200.000000" "0.010000"
write_freeeng "$RUN/adsorbate/calc.out"            "-13.500000"  "0.005000"
write_freeeng "$RUN/sites/site_1_ontop/calc.out"   "-214.100000" "0.020000"
write_freeeng "$RUN/sites/site_2_bridge/calc.out"  "-214.150000" "0.018000"
python3 -c "
e_slab, e_ads, e_ontop, e_bridge = -200.0, -13.5, -214.100000, -214.150000
print(f'  E_ads(ontop,  uncorrected) = {e_ontop  - e_slab - e_ads:.6f} eV')
print(f'  E_ads(bridge, uncorrected) = {e_bridge - e_slab - e_ads:.6f} eV   <- more negative, looks MORE stable')
"
pause

echo "--- Stage 2 (stb-adsorbBsse): BSSE ghost-fragment folders, at the RELAXED geometry ---"
echo "\$ stb-adsorbBsse --dir . --save-report"
(cd "$RUN" && stb-adsorbBsse --dir . --save-report --no-intro > "$OUT/workflow_stage2.log")
sed -n '/\[1\] SITE SCAN/,/\[2\]/p' "$OUT/workflow_stage2.log" | head -n -1
echo
echo "BSSE folders were written ONLY for the 2 finished sites:"
ls "$RUN/bsse"
echo
echo "Fabricated calc.out for the 4 ghost-fragment folders:"
write_freeeng "$RUN/bsse/site_1_ontop/bsse_slab/calc.out"        "-200.050000"
write_freeeng "$RUN/bsse/site_1_ontop/bsse_adsorbate/calc.out"   "-13.520000"
write_freeeng "$RUN/bsse/site_2_bridge/bsse_slab/calc.out"       "-200.090000"
write_freeeng "$RUN/bsse/site_2_bridge/bsse_adsorbate/calc.out"  "-13.560000"
python3 -c "
e_ontop, e_bridge = -214.100000, -214.150000
e_bsse_slab_ontop, e_bsse_ads_ontop = -200.050000, -13.520000
e_bsse_slab_bridge, e_bsse_ads_bridge = -200.090000, -13.560000
print(f'  E_ads_BSSE(ontop)  = {e_ontop  - e_bsse_slab_ontop  - e_bsse_ads_ontop:.6f} eV')
print(f'  E_ads_BSSE(bridge) = {e_bridge - e_bsse_slab_bridge - e_bsse_ads_bridge:.6f} eV')
"
pause

echo "--- Stage 3 (stb-adsorbAnalysis): does the correction change the answer? ---"
echo "\$ stb-adsorbAnalysis --dir adsorption_run --save-report --apply best_production.fdf"
(cd "$RUN" && stb-adsorbAnalysis --dir . --save-report --apply best_production.fdf --no-intro \
    > "$OUT/workflow_stage3.log")
sed -n '/Site          | Adsorbate/,/\[3\]/p' "$OUT/workflow_stage3.log" | head -n -1
sed -n '/\[3\] SUMMARY & PLOT/,/\[4\]/p' "$OUT/workflow_stage3.log" | head -n -1

UNCORRECTED_BEST=$(grep "Most stable site (uncorrected):" "$OUT/workflow_stage3.log" | grep -oE 'site_[a-z0-9_]+' | head -1)
BSSE_BEST=$(grep "Most stable site (BSSE-corrected):" "$OUT/workflow_stage3.log" | grep -oE 'site_[a-z0-9_]+' | head -1)
echo
if [ "$UNCORRECTED_BEST" = "site_2_bridge" ] && [ "$BSSE_BEST" = "site_1_ontop" ]; then
    echo "Confirmed: the uncorrected ranking picks '$UNCORRECTED_BEST', but the"
    echo "BSSE-corrected ranking picks '$BSSE_BEST' instead -- the correction"
    echo "changed which site is the answer."
else
    echo "Unexpected: uncorrected='$UNCORRECTED_BEST' bsse='$BSSE_BEST' (expected bridge / ontop)"
    exit 1
fi
[ -s "$RUN/best_production.fdf" ] && echo "--apply wrote '$RUN/best_production.fdf' from the BSSE-corrected winner."
pause

echo "=================================================================="
echo " Proof: CLI and the interactive stb-suite menu agree (Stage 1)"
echo "=================================================================="
echo "Driving 4.8.1 through the interactive menu (piped input, single ontop"
echo "site on the bare primitive cell) and comparing against the same"
echo "direct-CLI case from output/stage1/ above."
TMP="$(mktemp -d)"
cp structure.fdf calc.fdf "$TMP/"
echo
echo "\$ printf '4.8.1\\nstructure.fdf\\ncalc.fdf\\n\\nO\\n\\n1\\n\\n\\nn\\n0\\nn\\n\\ny\\nn\\nn\\n' | stb-suite"
(cd "$TMP" && printf '4.8.1\nstructure.fdf\ncalc.fdf\n\nO\n\n1\n\n\nn\n0\nn\n\ny\nn\nn\n\n0\n' \
    | stb-suite > menu1.log 2>&1)
MENU_ONTOP="$TMP/adsorption_run/sites/site_1_ontop/structure.fdf"
CLI_ONTOP="$OUT/stage1/sites/site_1_ontop/structure.fdf"
if diff -q <(grep -v '^#' "$MENU_ONTOP") <(grep -v '^#' "$CLI_ONTOP") > /dev/null 2>&1 \
    || python3 -c "
import sys
from stb.core import structure_io
import numpy as np
a = structure_io.to_pymatgen(structure_io.read_fdf('$MENU_ONTOP'))
b = structure_io.to_pymatgen(structure_io.read_fdf('$CLI_ONTOP'))
sys.exit(0 if np.allclose(sorted(a.cart_coords.tolist()), sorted(b.cart_coords.tolist()), atol=1e-6) else 1)
"; then
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
cat <<'EOF'
Folders generated under output/:
  stage1/                 stage1_lateral_warning/   stage1_fixed/
  stage1_bothsides/        workflow/

output/workflow/adsorption_run/ has the full chain: clean_slab/, adsorbate/,
sites/site_*/, bsse/site_*/, adsorption_report.txt, adsorption_bsse_report.txt,
adsorption_ranking.png, and best_production.fdf (the BSSE-corrected winner).

As a next step, on your OWN slab/2D structure:
  stb-adsorb -s <structure.fdf> -c <calc.fdf> --adsorbate <El-or-molecule> \
      --site-type all --all-sites
  # run SIESTA in clean_slab/, every adsorbate*/, and every sites/site_*/, then:
  stb-adsorbBsse --dir adsorption_run
  # run SIESTA in every bsse/site_*/bsse_slab/ and bsse/site_*/bsse_adsorbate/, then:
  stb-adsorbAnalysis --dir adsorption_run --save-report --apply production.fdf
EOF
