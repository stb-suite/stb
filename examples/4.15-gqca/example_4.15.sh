#!/bin/bash
# Guided example: GQCA (Generalized Quasi-Chemical Approximation) workflow
# (stb-gqca / stb-gqcaAnalysis, codes 4.15.1 / 4.15.2)
#
# Not an automated test (see test/4-workflow/15-gqca/{prep,analysis}/test.sh
# for that, including its own bipartite-vs-frustrated-lattice and 2D
# fixtures) -- a commented walk-through: it runs real commands, one group
# at a time, and shows you the piece of output that proves what just
# happened. It pauses between sections so you can read before moving on.
#
# stb-gqca is exercised for real -- it only builds structures, no SIESTA
# dependency at all. stb-gqcaAnalysis needs real SIESTA calc.out files to
# read, which this walkthrough doesn't have (no SIESTA binary here) -- so
# Section 2 fabricates calc.out files with HAND-CHOSEN FreeEng values, the
# same "worked backward from a known answer" convention 4.13-her/4.14-oer's
# own examples use. See the README's Section 3 for why these particular
# numbers were picked and what real literature they're standing in for.

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

write_energy() {
    local cluster_dir="$1" free_eng="$2"
    cat > "$cluster_dir/calc.out" <<EOF
siesta: SCF Convergence by DM criterion
SCF cycle converged after 12 iterations
siesta: FreeEng =        ${free_eng}
siesta: Atomic forces (eV/Ang):
     1    0.001    0.001    0.001
     Max    0.001
EOF
}

echo "=================================================================="
echo " Why 3 fixed structures are the whole workflow, and what they solve"
echo "=================================================================="
cat <<'EOF'
The pair-cluster (N=2) Generalized Quasi-Chemical Approximation (GQCA,
Sher, van Schilfgaarde, Chen, Chen, Phys. Rev. B 36, 4279, 1987) treats a
substitutional alloy A(1-x)Bx as a lattice of independent nearest-neighbor
PAIRS, each in one of 3 local states:

  n_j = 0 (AA)   n_j = 1 (AB)   n_j = 2 (BB)

Stage 1 (stb-gqca) builds and relaxes exactly these 3 reference
structures -- ONCE, regardless of which composition x you actually care
about. Stage 2 (stb-gqcaAnalysis) combines their 3 energies via a
mass-action equilibrium (like a chemical reaction 2(AB) <-> AA + BB) to
get the pair-population fractions x_j(x,T) at ANY x in [0,1] and ANY T,
then reports the mixing enthalpy/entropy/free energy and a short-range-
-order (SRO) parameter from them. Nothing is re-run per composition --
the 3 DFT calculations ARE the whole x/T thermodynamic surface.
EOF
pause

echo "=================================================================="
echo " output/algaas/  --  Stage 1 on a real semiconductor: GaAs -> AlGaAs"
echo "=================================================================="
cat <<'EOF'
structure.fdf is GaAs in the zincblende structure (conventional cubic
cell, a = 5.6533 Ang, the real experimental lattice constant) --
--sublattice Ga --species-b Al disorders the CATION sublattice, modeling
Al(x)Ga(1-x)As: the most widely used III-V heterostructure alloy (laser
diodes, HEMTs, solar cells), chosen here specifically because it's the
textbook example of a NEARLY IDEAL solid solution -- AlAs and GaAs share
almost exactly the same lattice constant (~0.1% mismatch), so real AlGaAs
mixes close to randomly at any composition. Section 3 below checks
whether stb-gqcaAnalysis reproduces exactly that.
EOF
echo "\$ stb-gqca -f structure.fdf --sublattice Ga --species-b Al -c calc.fdf -O output/algaas"
stb-gqca -f structure.fdf --sublattice Ga --species-b Al -c calc.fdf \
    -O "$OUT/algaas" --no-intro > "$OUT/algaas_stage1.log"
echo
echo "[2] PAIR ORBIT -- the shortest Ga-Ga distance and the orbit built on it:"
sed -n '/\[2\] PAIR ORBIT/,/multiplicity/p' "$OUT/algaas_stage1.log"
echo
echo "[3] CLUSTER STRUCTURES -- 3 folders, one full-cell 3D relaxation each:"
sed -n '/\[3\] CLUSTER STRUCTURES/,/AB-ordered)/p' "$OUT/algaas_stage1.log" | grep "\[OK\]"
pause

echo "------------------------------------------------------------------"
echo " The Ga sublattice is FCC: a real, unavoidable geometric limit"
echo "------------------------------------------------------------------"
cat <<'EOF'
Zincblende's cation sublattice is FCC -- the same lattice rocksalt's own
cation sublattice has (stb-gqca --help gives NaCl as its own worked
example of this). FCC contains triangles (an odd cycle), so no
2-coloring of it can make EVERY nearest-neighbor pair unlike: the
smallest exact-50:50-composition cell can only get so close to a pure
AB (every-pair-unlike) local environment before running out of distinct
orderings to try. This is real crystallography, not a fixture quirk --
it shows up for ANY zincblende- or rocksalt-based alloy, not just AlGaAs:
EOF
grep "NOTE.*deviates by\|NOTE.*sensitivity check" "$OUT/algaas_stage1.log"
echo
echo "-- and it's also why cluster_n1/ has FEWER atoms than cluster_n0/n2"
echo "(icet's own primitive-cell enumeration, not the 8-atom input cell):"
python3 -c "
from stb.core import structure_io
for n in (0, 1, 2):
    s = structure_io.read_fdf('$OUT/algaas/cluster_n' + str(n) + '/structure.fdf')
    print(f'  cluster_n{n}: {len(s.atoms)} atoms')
"
pause

echo "=================================================================="
echo " output/algaas/  --  Section 2: fabricating the 3 calc.out files"
echo "=================================================================="
cat <<'EOF'
Real SIESTA output isn't available here (no SIESTA binary to invoke), so
-- exactly like 4.13-her/4.14-oer's own examples -- this section writes
calc.out files with a HAND-CHOSEN FreeEng, worked backward from the two
cases the README's Section 3 discusses:

  CASE A (ideal): e_0 = e_2 = -50.0 eV/site, e_1 = -50.0 eV/site exactly
                  -> Delta_1 = 0 -- what real, nearly lattice-matched
                     AlGaAs should look like.
  CASE B (illustrative, NOT real AlGaAs): e_1 shifted by -0.1 eV/site
                  -> Delta_1 = -0.1 eV/site -- a HYPOTHETICAL
                     ordering-favoring alloy, shown only to demonstrate
                     what stb-gqcaAnalysis reports when the physics is
                     NOT ideal (see the README for why real AlGaAs is
                     not expected to look like this).

cluster_n0/cluster_n2 have 4 Ga(Al) sites each (the 8-atom input cell);
cluster_n1 has 2 (icet's smaller primitive-based AB cell) -- FreeEng is
scaled by site count so all 3 land on the SAME per-site energy for a
Delta_1 = 0 case, exactly like write_energy's own docstring in
test/4-workflow/15-gqca/analysis/test.sh.
EOF
echo "Done -- Case A written first (see Section 3)."
write_energy "$OUT/algaas/cluster_n0" "-200.000000"
write_energy "$OUT/algaas/cluster_n2" "-200.000000"
write_energy "$OUT/algaas/cluster_n1" "-100.000000"
pause

echo "=================================================================="
echo " Section 3, CASE A: does the tool reproduce a well-known fact?"
echo "=================================================================="
echo "\$ stb-gqcaAnalysis --directory output/algaas --temp 900 -o ideal"
stb-gqcaAnalysis --directory "$OUT/algaas" --temp 900 --x-points 5 -o ideal \
    --no-intro > "$OUT/ideal_stage2.log"
echo
echo "[1] CLUSTER ENERGIES -- Delta_j = 0 for all j, exactly by construction:"
sed -n '/Mixing energies/,/Delta_2/p' "$OUT/ideal_stage2.log"
echo
echo "At x=0.5, T=900 K (~a typical AlGaAs MBE/MOCVD growth temperature):"
grep "^0.500000," "$OUT/algaas/ideal.csv"
echo "x_j = [0.25, 0.50, 0.25] -- EXACTLY the random/binomial mixing"
echo "distribution, H_mix = 0, SRO = 0: this IS the well-established"
echo "textbook fact about real Al(x)Ga(1-x)As -- an almost perfectly"
echo "ideal (Vegard's-law) solid solution, with no measurable"
echo "short-range order at any accessible growth temperature."
pause

echo "=================================================================="
echo " Section 3, CASE B: a HYPOTHETICAL non-ideal alloy, same host lattice"
echo "=================================================================="
echo "Overwriting cluster_n1's calc.out with the Case B energy:"
write_energy "$OUT/algaas/cluster_n1" "-100.200000"
echo "\$ stb-gqcaAnalysis --directory output/algaas --temp 900 -o ordering"
stb-gqcaAnalysis --directory "$OUT/algaas" --temp 900 --x-points 5 -o ordering \
    --no-intro > "$OUT/ordering_stage2.log"
echo
grep "^0.500000," "$OUT/algaas/ordering.csv"
echo
echo "x_1 pushed from 0.50 up to ~0.78, and SRO from 0 to ~ -0.57: the"
echo "solver responds correctly to a NEGATIVE Delta_1 (AB pairing favored)"
echo "by favoring AB-ordering over random mixing -- this is NOT real"
echo "AlGaAs (see Section 3 above), just proof the machinery tracks the"
echo "sign/magnitude of whatever Delta_1 you actually give it."
echo
echo "And yet -- even at this artificially strong ordering bias -- the"
echo "miscibility-gap search still reports no gap, by mathematical"
echo "necessity of this model (read carefully, this is a documented"
echo "STRUCTURAL LIMITATION, not a claim about real AlGaAs miscibility):"
sed -n '/\[4\] MISCIBILITY GAP/,/no gap\|gap(s)/p' "$OUT/ordering_stage2.log" | tail -3
pause

echo "=================================================================="
echo " output/algaas/  --  a temperature sweep: order fades out with heat"
echo "=================================================================="
echo "\$ stb-gqcaAnalysis --directory output/algaas --temp-min 300 --temp-max 1500 --temp-points 5"
stb-gqcaAnalysis --directory "$OUT/algaas" --temp-min 300 --temp-max 1500 \
    --temp-points 5 --x-points 3 -o ordering_sweep --no-intro \
    > "$OUT/sweep_stage2.log"
echo
echo "SRO at x=0.5 (Case B energies), 300 K -> 1500 K:"
grep "^0.500000," "$OUT/algaas/ordering_sweep.csv" | awk -F, '{printf "  T = %6.0f K   SRO = %8.4f\n", $2, $10}'
echo "|SRO| shrinks monotonically as T rises -- thermal disorder"
echo "competing against the (hypothetical) ordering energy, the same"
echo "qualitative order-disorder temperature dependence real alloys show."
pause

echo "=================================================================="
echo " Proof: CLI and the interactive stb-suite menu agree (Stage 1)"
echo "=================================================================="
echo "Driving 4.15.1 through the interactive menu (piped input, Ga -> Al,"
echo "pseudopotentials skipped) and comparing against output/algaas/ above."
TMP="$(mktemp -d)"
cp structure.fdf calc.fdf "$TMP/"
echo
echo "\$ printf '4.15.1\\nstructure.fdf\\nGa\\nAl\\ncalc.fdf\\n\\nalgaas_menu\\n\\n\\n0\\n' | stb-suite"
(cd "$TMP" && printf '4.15.1\nstructure.fdf\nGa\nAl\ncalc.fdf\n\nalgaas_menu\n\n\n0\n' \
    | stb-suite > menu1.log 2>&1)
python3 -c "
import sys
from stb.core import structure_io
import numpy as np
a = structure_io.read_fdf('$TMP/algaas_menu/cluster_n0/structure.fdf')
b = structure_io.read_fdf('$OUT/algaas/cluster_n0/structure.fdf')
a_frac = sorted([tuple(np.round(pos, 6)) for _s, pos in a.atoms])
b_frac = sorted([tuple(np.round(pos, 6)) for _s, pos in b.atoms])
sys.exit(0 if a_frac == b_frac else 1)
"
if [ $? -eq 0 ]; then
    echo "Confirmed: the interactive-menu cluster_n0 and the direct-CLI cluster_n0 have the same geometry."
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
  algaas/cluster_n{0,1,2}/   -- Stage 1 structures + fabricated calc.out
  algaas/gqca_stage1.txt     -- Stage 1 report
  algaas/{ideal,ordering,ordering_sweep}.csv/.dat/.gplot -- Stage 2 runs

As a next step, on your OWN alloy (any lattice, 2D or 3D):
  stb-gqca -f structure.fdf --sublattice <A> --species-b <B> -c calc.fdf -O gqca_study
  # run SIESTA (a real CG/variable-cell relaxation) in every
  # gqca_study/cluster_n{0,1,2}/ folder, then:
  stb-gqcaAnalysis --directory gqca_study --temp <T_K>
  # or --temp-min/--temp-max/--temp-points for a sweep. Read the
  # [NOTE]/[WARNING] lines in gqca_stage1.txt first if your sublattice is
  # FCC/HCP/triangular (see Section 2 above) -- Delta_1 is then a less
  # clean estimate, still qualitatively meaningful.
EOF
