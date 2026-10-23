#!/bin/bash
# Guided example: Workflow 4.19 -- Charge Density Difference (both stages: stb-chargediffPrep, code 4.19.1,
# and stb-chargediffAnalysis, code 4.19.2, in the stb-suite menu).
#
# Not an automated test (see test/4-workflow/19-chargediff/{prep,analysis}/test.sh for that) -- a commented
# walk-through: it runs real commands, one case at a time, into its own output/<case>/ folder, and shows the
# piece of output that proves what just happened. Pauses between sections so you can read before moving on.
# Safe to re-run any time -- it always starts by wiping its own output/.
#
# Stage 1 (stb-chargediffPrep) runs for real in every case below (cheap -- it only writes .fdf files, no SCF).
# Stage 2 (stb-chargediffAnalysis) runs on REAL SIESTA results: reference_results/ holds the finished runs of
# the water dimer, ammonia borane and Li on graphene, plus the water-dimer studies (basis, mesh, distance),
# computed with the Stage 1 commands of this script and the inputs in this folder (SIESTA 5.4.2,
# mpirun -np 2; see README section 6 and make_reference.sh to regenerate them). The density grids are
# stored compressed (siesta.RHO.xz) and unpacked into output/ when a case needs them. Nothing in Stage 2 is
# fabricated; this script never runs SIESTA.

set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

# --show needs a display; MPLBACKEND=Agg keeps this script headless-safe (the figure is still saved with --png).
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

# section N of a report: from "[N] TITLE" up to the next "[M] " header
section() {   # $1 = file, $2 = section number
    strip_colors < "$1" | awk -v n="[$2]" 'index($0, n) == 1 {on = 1; print; next}
        on && /^\[[0-9]+\] / {exit} on {print}'
}

# a finished reference set, copied into output/ with its density grids unpacked
unpack() {   # $1 = set name, $2 = destination
    cp -r "$REF/$1" "$2"
    for f in "$2"/*/siesta.RHO.xz; do [ -f "$f" ] && xz -d "$f"; done
    return 0
}


banner "Welcome: where do the electrons go when two things meet?"
cat <<'EOF'
Put two fragments together -- two molecules, a molecule and a surface, two layers -- and their electrons
rearrange. The charge density difference

    Delta rho(r) = rho_AB(r) - rho_A(r) - rho_B(r)

computed with every density on the SAME geometry, grid and basis, shows where charge piles up (bonds,
lone pairs pushed toward a partner) and where it is taken from. Its integrals give how much charge was
moved and how much crossed from one fragment to the other.

Stage 1 (stb-chargediffPrep) writes one SIESTA folder for the combined system and one per fragment. Each
fragment keeps every atom of the others as GHOST atoms (orbitals without nucleus or electrons), so all
runs share one basis set; config_extra.fdf forces the single point, the fixed cell, the same MeshCutoff
and SaveRho. Stage 2 (stb-chargediffAnalysis) subtracts the densities and reports the consistency of the
runs, the counterpoise-corrected interaction energy, the statistics of Delta rho, the charge change of
each fragment and atom, and a dashboard.

This tour: Stage 1 on the water dimer, its folders, the other ways to split a system (Cases 1-3); then
Stage 2 on real SIESTA runs: the water dimer (Cases 4-5), the basis-set superposition error, the mesh and
the distance (Cases 6-8), ammonia borane (Case 9), Li on graphene with the planar profile (Case 10), and
the guard rails (Case 11). See README.md for the theory (section 1) and the worked example (section 5).
EOF
pause


banner "Case 1: Stage 1 -- the water dimer, what stb-chargediffPrep reports"
cat <<'EOF'
water_dimer.fdf is the relaxed (H2O)2 in a 14 Ang box, the acceptor first (O H H), then the donor; calc.fdf
is the template (PBE, DZP, PAO.EnergyShift 0.01 Ry, %include structure.fdf). --fragment-sizes 3,3 splits
the atoms in file order. -p dojo copies the PseudoDojo pseudopotentials, renamed after each species label.
Watch [2] (the fragments, their valence electrons and the closest contact between them), [3] (the
template's level of theory and what config_extra.fdf forces) and [4] (real and ghost atoms per folder).
EOF
mkdir -p "$OUT/case1-water"
cp water_dimer.fdf calc.fdf "$OUT/case1-water/"
echo "\$ stb-chargediffPrep -s water_dimer.fdf --fragment-sizes 3,3 -p dojo --mesh-cutoff 320 --png --save-report --no-intro"
(cd "$OUT/case1-water" && stb-chargediffPrep -s water_dimer.fdf --fragment-sizes 3,3 -p dojo --mesh-cutoff 320 \
    --png --save-report --no-intro > stage1.log 2>&1)
section "$OUT/case1-water/stage1.log" 2 | cut -c1-150
section "$OUT/case1-water/stage1.log" 3 | cut -c1-150
section "$OUT/case1-water/stage1.log" 4 | cut -c1-150
pause


banner "Case 2: inside the folders -- labels, ghosts, config_extra.fdf"
cat <<'EOF'
Every species of every fragment gets its own label (O_frag1, H_frag2, ...) and its own pseudopotential
file, as in the adsorption workflow. In a fragment folder the atoms of the OTHER fragments are present as
<label>_ghost species with a NEGATIVE atomic number: SIESTA then places their orbitals but no nucleus and
no electrons. The calc.fdf of every folder is '%include config_extra.fdf' followed by your template byte for
byte; SIESTA's fdf reader is first-occurrence-wins, so the forced settings beat the template's.
EOF
R="$OUT/case1-water/chargediff_run"
ls "$R"; echo "--- fragment1/:"; ls "$R/fragment1"
echo
echo "\$ grep -A5 ChemicalSpeciesLabel chargediff_run/fragment1/structure.fdf"
grep -A5 "block ChemicalSpeciesLabel" "$R/fragment1/structure.fdf"
echo "\$ grep -A5 ChemicalSpeciesLabel chargediff_run/fragment2/structure.fdf"
grep -A5 "block ChemicalSpeciesLabel" "$R/fragment2/structure.fdf"
echo
echo "\$ grep -v '^#' chargediff_run/fragment1/config_extra.fdf"
grep -v "^#" "$R/fragment1/config_extra.fdf" | grep -v "^$"
echo
echo "\$ head -2 chargediff_run/fragment1/calc.fdf; diff <(tail -n +3 chargediff_run/fragment1/calc.fdf) calc.fdf"
head -2 "$R/fragment1/calc.fdf"
diff <(tail -n +3 "$R/fragment1/calc.fdf") calc.fdf && echo "-> the template is copied byte for byte below the %include"
echo
echo "\$ python3 -c ...   # chargediff_manifest.json, what Stage 2 reads"
python3 - "$R/chargediff_manifest.json" <<'EOF2'
import json, sys
m = json.load(open(sys.argv[1]))
print(f"  combined: {m['combined_folder']}/   fragments: {', '.join(m['fragment_folders'])}")
for f in m["fragments"]:
    print(f"  {f['folder']}: {f['n_atoms']} atoms, labels {', '.join(f['species_labels'])}")
EOF2
echo "  (figure of the split: output/case1-water/chargediff_run/chargediff_fragments.png)"
pause


banner "Case 3: other splits -- by species, by an adsorption manifest; spin and dipole"
cat <<'EOF'
--fragment-species groups atoms by element, whatever their order: 'C;Li' puts every carbon in one fragment
and the lithium in the other (it cannot separate two fragments that share an element: ammonia borane has
H in both halves, so it needs --fragment-sizes). --fragment-manifest reuses the slab/adsorbate split that
stb-adsorb (4.8) records next to each site. --force-spin / --force-dipole / --force-vdw put Spin polarized,
the slab dipole correction and DFT-D3 into EVERY folder (the default is to inherit them from the combined
folder's own config_extra.fdf). A split that does not add up is refused before anything is written.
EOF
mkdir -p "$OUT/case3-splits"
cp li_graphene.fdf calc_graphene.fdf nh3bh3.fdf calc.fdf "$OUT/case3-splits/"
cd "$OUT/case3-splits"
echo "\$ stb-chargediffPrep -s li_graphene.fdf --fragment-species 'C;Li' --calc calc_graphene.fdf --force-spin --force-dipole -p dojo"
stb-chargediffPrep -s li_graphene.fdf --fragment-species "C;Li" --calc calc_graphene.fdf --force-spin --force-dipole \
    -p dojo --mesh-cutoff 300 -o li_run --no-intro 2>&1 | strip_colors | sed -n '/^\[2\]/,/^\[3\]/p' | grep -E "Split|fragment[0-9]/" | cut -c1-150
grep -v "^#" li_run/fragment2/config_extra.fdf | grep -v "^$" | sed 's/^/  /'
echo
cat > slab_manifest.json <<'EOF2'
{"schema": "adsorb_fragment_manifest_v1", "site_label": "li_on_graphene", "slab_suffix": "_slab",
 "adsorbate_suffix": "_ads", "n_substrate": 18, "n_adsorbate": 1,
 "slab_species": {"C_slab": {"id": "1", "Z": 6, "count": 18}},
 "adsorbate_species": {"Li_ads": {"id": "2", "Z": 3, "count": 1}}, "adsorbate_formula": "Li"}
EOF2
# an adsorption site of stb-adsorb carries its own labels (C_slab, Li_ads) and the manifest of its split
sed 's/^ 1 6 C$/ 1 6 C_slab/; s/^ 2 3 Li$/ 2 3 Li_ads/' li_graphene.fdf > site.fdf
echo "\$ stb-chargediffPrep -s site.fdf --fragment-manifest slab_manifest.json ...   # an stb-adsorb-style site"
stb-chargediffPrep -s site.fdf --fragment-manifest slab_manifest.json --calc calc_graphene.fdf -p dojo \
    -o manifest_run --no-intro 2>&1 | strip_colors | grep -E "^(slab|adsorbate)/" | cut -c1-150
echo
echo "\$ stb-chargediffPrep -s nh3bh3.fdf --fragment-species 'N;B'      # H is in both halves: not expressible by species"
set +e
stb-chargediffPrep -s nh3bh3.fdf --fragment-species "N;B" -o x1 --no-intro 2>&1 | strip_colors | grep -E "ERROR" | cut -c1-200
echo "\$ stb-chargediffPrep -s nh3bh3.fdf --fragment-sizes 4,3           # 7 atoms of 8"
stb-chargediffPrep -s nh3bh3.fdf --fragment-sizes 4,3 -o x2 --no-intro 2>&1 | strip_colors | grep -E "ERROR" | cut -c1-200
set -e
cd "$DIR"
pause


banner "Case 4: Stage 2 -- the water dimer from real SIESTA runs"
cat <<'EOF'
reference_results/water/ is Case 1's folder after SIESTA ran in combined/, fragment1/ and fragment2/
(PBE, DZP, 320 Ry requested). Read [1] (each folder: atoms, real or ghost, the integral of its density and
the electron count SIESTA reports), [2] (is the subtraction legitimate? the same grid, cell, level of
theory, pseudopotentials, mesh SIESTA actually used, SCF convergence, electron count) and [3] (the
counterpoise-corrected interaction energy -- three numbers from the same runs, no extra calculation).
EOF
mkdir -p "$OUT/case4-water"
unpack water "$OUT/case4-water/water"
echo "\$ stb-chargediffAnalysis --dir water --save-report --no-intro"
(cd "$OUT/case4-water" && stb-chargediffAnalysis --dir water -o analysis --save-report --no-intro > stage2.log 2>&1)
section "$OUT/case4-water/stage2.log" 1 | cut -c1-150
section "$OUT/case4-water/stage2.log" 2 | cut -c1-150
section "$OUT/case4-water/stage2.log" 3 | cut -c1-200
pause


banner "Case 5: what moved, and where -- statistics, fragments, atoms, files"
cat <<'EOF'
[4]: the integral of Delta rho (zero: charge is only moved), the charge moved (half the integral of
|Delta rho|), the extremes and where they are, and a first isosurface level for VESTA in the unit of the
cube file (e/bohr^3). [5]: the grid split among the atoms by distance (nearest-atom partition) -- the net
charge each fragment and atom gained or lost. Then the outputs: nothing beyond the slice .dat and the
cube unless asked: --png the dashboard, --save-gnuplot the .gplot, --save-data a CSV and a JSON.
EOF
section "$OUT/case4-water/stage2.log" 4 | cut -c1-150
section "$OUT/case4-water/stage2.log" 5 | cut -c1-150
echo
echo "\$ stb-chargediffAnalysis --dir water --png --save-gnuplot --save-data -o files --no-intro"
(cd "$OUT/case4-water" && stb-chargediffAnalysis --dir water -o files --png --save-gnuplot --save-data --no-intro > files.log 2>&1)
ls -1 "$OUT/case4-water/files" | sed 's/^/  /'
cp "$OUT/case4-water/files/water_chargediff_summary.png" "$DIR/figures/water_dashboard.png"
echo "  (figures/water_dashboard.png is this dashboard)"
pause


banner "Case 6: the basis-set superposition error -- why the fragments carry ghost atoms"
cat <<'EOF'
Each reference set also ran every fragment WITHOUT the ghost atoms (fragmentN_alone/, its own basis only).
The difference is the basis-set superposition error: inside the complex each fragment borrows its
partner's orbitals and is artificially lowered; the counterpoise energy removes it. [A] gives the three
cases, with the deformation energy of the fragments (from the relaxed monomers) that turns an interaction
energy into a binding energy; [B] repeats the water dimer with five basis sets at the same geometry.
EOF
echo "\$ python3 study_table.py      # [A] and [B]"
python3 study_table.py 2>&1 | sed -n '/^\[A\]/,/^\[C\]/p' | sed '$d'
pause


banner "Case 7: the mesh -- the same grid for every density, and how fine"
cat <<'EOF'
Delta rho is a point-by-point subtraction: every density must live on the same grid, so the MeshCutoff is
forced into every folder. [C] repeats the dimer at five cutoffs. Note the column 'used': SIESTA rounds the
requested cutoff up to a grid that fits the cell, and reports the value it actually used.
EOF
python3 study_table.py 2>&1 | sed -n '/^\[C\]/,/^\[D\]/p' | sed '$d'
echo "--- and a set where one folder ran on another mesh (a copy of 'water' with fragment2's mesh edited):"
mkdir -p "$OUT/case7-mesh"
unpack water "$OUT/case7-mesh/mixed"
sed -i 's/Mesh cutoff (required, used) = *\([0-9.]*\) *\([0-9.]*\)/Mesh cutoff (required, used) =   450.000   487.111/' \
    "$OUT/case7-mesh/mixed/fragment2/calc.out"
(cd "$OUT/case7-mesh" && stb-chargediffAnalysis --dir mixed -o a --no-cube --no-intro 2>&1 | strip_colors | grep -E "^Same mesh|WARNING\] Same mesh" | cut -c1-170)
pause


banner "Case 8: the distance -- charge transfer is short-ranged, polarization is not"
cat <<'EOF'
The donor molecule moved rigidly along the O...O axis from 2.5 to 5 Ang (nine more real sets). [D] gives
the interaction energy with and without the counterpoise correction, the net charge the donor gains and
the charge moved. The net charge transfer dies off exponentially; the polarization (charge moved inside
each molecule) decays more slowly. figures/water_studies.png draws it.
EOF
python3 study_table.py --png 2>&1 | sed -n '/^\[D\]/,$p'
pause


banner "Case 9: ammonia borane -- a dative bond, an order of magnitude stronger"
cat <<'EOF'
H3N-BH3: the nitrogen lone pair is donated into the empty boron p orbital. Compare with the water dimer:
the interaction energy is ~11x larger, the charge crossing between the fragments ~7x, and the deformation energy (BH3 is planar
alone and pyramidal in the complex) is a large part of the story: here the interaction energy and the
binding energy are very different numbers.
EOF
mkdir -p "$OUT/case9-nh3bh3"
unpack nh3bh3 "$OUT/case9-nh3bh3/nh3bh3"
(cd "$OUT/case9-nh3bh3" && stb-chargediffAnalysis --dir nh3bh3 -o a --png --save-data --no-intro > stage2.log 2>&1)
section "$OUT/case9-nh3bh3/stage2.log" 3 | grep -E "E_int|E\(" | cut -c1-150
section "$OUT/case9-nh3bh3/stage2.log" 4 | cut -c1-150
section "$OUT/case9-nh3bh3/stage2.log" 5 | cut -c1-150
cp "$OUT/case9-nh3bh3/a/nh3bh3_chargediff_summary.png" "$DIR/figures/nh3bh3_dashboard.png"
echo "  (figures/nh3bh3_dashboard.png is its dashboard)"
pause


banner "Case 10: Li on graphene -- a periodic interface, the planar profile, spin"
cat <<'EOF'
A Li atom on a 3x3 graphene sheet: every run spin-polarized (the isolated Li has one unpaired electron)
and with the slab dipole correction, both forced by Stage 1. For an interface the natural view is the
PLANAR AVERAGE of Delta rho along the normal (--profile --axis 2) and its running integral: the charge
that crossed each plane. The report gives, without any partition of space, the charge on each side of the
plane halfway between Li and the sheet, and the dipole Delta rho induces -- checked against SIESTA's own
dipoles of the three runs -- with the work-function shift it causes.
EOF
mkdir -p "$OUT/case10-li"
unpack li_graphene "$OUT/case10-li/li_graphene"
(cd "$OUT/case10-li" && stb-chargediffAnalysis --dir li_graphene -o a --profile --axis 2 --png --save-data \
    --save-gnuplot --no-intro > stage2.log 2>&1)
section "$OUT/case10-li/stage2.log" 2 | cut -c1-150
section "$OUT/case10-li/stage2.log" 3 | grep -E "E_int|E\(" | cut -c1-150
section "$OUT/case10-li/stage2.log" 5 | grep -E "^19|^fragment[0-9] |Planar|Separating|Charge across|partition above|Induced|Work function" | cut -c1-170
cp "$OUT/case10-li/a/li_graphene_chargediff_summary.png" "$DIR/figures/li_graphene_dashboard.png"
echo "  (figures/li_graphene_dashboard.png is its dashboard)"
echo
echo "\$ stb-chargediffAnalysis --dir li_graphene --spin --profile --axis 2      # the spin-density difference"
echo "  (the isolated Li carries 1 mu_B; the combined single point ended with a small moment of its own -- README 5.6)"
for f in combined fragment1 fragment2; do
    echo "  $f: $(grep 'spin moment' "$OUT/case10-li/li_graphene/$f/calc.out" | tail -1 | awk '{print $NF}') mu_B"
done
(cd "$OUT/case10-li" && stb-chargediffAnalysis --dir li_graphene -o spin --spin --profile --axis 2 --no-cube \
    --save-data --no-intro > spin.log 2>&1)
section "$OUT/case10-li/spin.log" 4 | grep -E "Integral|Moment" | cut -c1-150
pause


banner "Case 11: guard rails -- what Stage 2 refuses or flags"
cat <<'EOF'
A fragment whose density sits on another grid (it cannot be subtracted point by point: refused), a run
whose SCF never converged (flagged), a folder with a different level of theory (flagged), and --spin on
runs that are not spin-polarized (refused). Each is reported where it applies and again at the end.
EOF
mkdir -p "$OUT/case11-guards"
cd "$OUT/case11-guards"
unpack water a_scf; sed -i '/SCF cycle converged/d' a_scf/fragment1/calc.out
unpack water b_theory; echo "XC.authors revPBE" >> b_theory/fragment2/config_extra.fdf
unpack water c_grid
python3 - <<'EOF2'
import sisl
g = sisl.get_sile("c_grid/fragment2/siesta.RHO").read_grid()
small = sisl.Grid([n // 2 for n in g.shape], geometry=g.geometry)
sisl.get_sile("c_grid/fragment2/siesta.RHO", mode="w").write_grid(small)
EOF2
set +e
for d in a_scf b_theory; do
    echo "\$ stb-chargediffAnalysis --dir $d      # warnings only"
    stb-chargediffAnalysis --dir $d -o out_$d --no-cube --no-intro 2>&1 | strip_colors | sed -n '/^\[9\] WARNINGS/,/^\[INFO\]/p' | grep -E "^ *[0-9]+\. " | cut -c1-200
done
echo "\$ stb-chargediffAnalysis --dir c_grid"
stb-chargediffAnalysis --dir c_grid -o out_c --no-intro 2>&1 | strip_colors | grep ERROR | cut -c1-200
echo "\$ stb-chargediffAnalysis --dir a_scf --spin"
stb-chargediffAnalysis --dir a_scf --spin -o out_s --no-intro 2>&1 | strip_colors | grep ERROR | cut -c1-200
set -e
cd "$DIR"

echo
echo "Done. Every case wrote into $OUT/; the figures of the README are in figures/."
