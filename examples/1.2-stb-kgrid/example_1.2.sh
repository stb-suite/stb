#!/bin/bash
# Guided example: stb-kgrid (code 1.2 in the stb-suite menu)
#
# Not an automated test (see test/1-inputs/2-k-grid/test.sh for that) -- a
# commented walk-through: it runs real commands, one group at a time, and
# shows you the piece of output that proves what just happened. It pauses
# between sections so you can read before moving on.

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

# Prints the total number of k-points in a "N N N" grid string.
total_kpoints() {
    echo "$1" | awk '{print $1*$2*$3}'
}

echo "=================================================================="
echo " Why a DFT calculation needs a k-grid at all"
echo "=================================================================="
cat <<'EOF'
A periodic crystal's electronic states are labeled by a wavevector k
(Bloch's theorem) that ranges continuously over the Brillouin zone (BZ) --
the reciprocal-space unit cell. Total energy, forces, and everything else
SIESTA computes are, in principle, an INTEGRAL over that whole zone. In
practice that integral is approximated by a weighted sum over a finite set
of sampled k-points -- too few of them and the sum is a poor approximation
of the true integral (the result changes noticeably if you add more);
enough of them and the result converges to a stable value.

"Monkhorst-Pack" (Monkhorst & Pack, 1976 -- the exact reference stb-kgrid
itself writes to references.bib below) is the standard convention for HOW
to lay out that finite set: a uniform mesh of N1 x N2 x N3 points along the
three reciprocal lattice directions, instead of picking points ad hoc.

stb-kgrid does NOT re-derive that original paper's own numerical-optimality
argument -- what it actually computes is simpler and more practical: given
a target k-point DENSITY (in 1/Ang), it picks a mesh size per axis,
  N_i = ceil( |b_i| / density )
where b_i is the length of reciprocal lattice vector i. A finer density
(smaller number) means more k-points, a more accurate (and slower)
calculation; a coarser one means fewer k-points, faster but less accurate.
Metals (sharp features at the Fermi surface) typically need a denser mesh
than semiconductors/insulators for the same accuracy.
EOF
pause

echo "=================================================================="
echo " Two ways to run it"
echo "=================================================================="
cat <<'EOF'
A -- direct CLI:
  stb-kgrid --file structure.fdf --density 0.2

B -- interactive stb-suite menu:
  $ stb-suite
  Select an option (0-6, or a tool code like 4.1.2): 1.2

"1.2" asks for the structure file, then -- as of this session -- shows the
density recommendation guide (below) BEFORE asking for a density, so the
choice is informed instead of guessed; then runs the same stb-kgrid command
underneath. Proven identical to the CLI further below.
EOF
pause

echo "=================================================================="
echo " output/silicon/  --  bulk silicon at density 0.2"
echo "=================================================================="
echo "The most common case: a real 3D bulk crystal."
mkdir -p "$OUT/silicon"
cp structure.fdf "$OUT/silicon/"
echo
echo "\$ stb-kgrid --file structure.fdf --density 0.2 --save-report --no-intro"
(cd "$OUT/silicon" && stb-kgrid --file structure.fdf --density 0.2 --save-report --no-intro > console.log)
grep -E "Dimensionality|Suggested Monkhorst-Pack" "$OUT/silicon/console.log"
echo
echo "Folder contents (report + citations, saved via --save-report):"
ls "$OUT/silicon"
pause

echo "=================================================================="
echo " The accuracy/cost tradeoff, made concrete"
echo "=================================================================="
echo "Same structure, three densities -- watch the k-point COUNT (the actual"
echo "cost driver: SCF cost scales with the number of sampled k-points)."
TMP="$(mktemp -d)"
cp structure.fdf "$TMP/"
for density in 0.1 0.2 0.3; do
    grid=$(cd "$TMP" && stb-kgrid --file structure.fdf --density "$density" --no-intro \
        | grep "Suggested Monkhorst-Pack" | sed 's/.*: //')
    n=$(total_kpoints "$grid")
    printf "  density=%-4s -> grid %-10s -> %s k-points\n" "$density" "$grid" "$n"
done
rm -rf "$TMP"
echo
echo "Going from 0.3 to 0.1 (3x finer) costs 27x more k-points -- the mesh"
echo "size scales with the CUBE of the linear density change for a 3D solid."
pause

echo "=================================================================="
echo " Why dimensionality matters here too"
echo "=================================================================="
cat <<'EOF'
Bloch's theorem, and the whole idea of a k-point, only applies along a
genuinely PERIODIC direction. stb-kgrid detects a vacuum-padded axis (the
same core/kspace.py::detect_vacuum_axes heuristic stb-inputfile also uses)
and forces a single k-point there regardless of density -- there is no
dispersion to sample along a direction with no real periodicity, so
"finer" sampling there wouldn't mean anything physically.
EOF
pause

for case in molecule chain graphene; do
    case "$case" in
        molecule) struct="structure_molecule.fdf"; label="molecule/ -- CH4, isolated (0D)" ;;
        chain)    struct="structure_chain.fdf";    label="chain/ -- carbon chain (1D)" ;;
        graphene) struct="structure_graphene.fdf"; label="graphene/ -- graphene monolayer (2D)" ;;
    esac

    echo "=================================================================="
    echo " output/$case  --  $label"
    echo "=================================================================="
    mkdir -p "$OUT/$case"
    cp "$struct" "$OUT/$case/"
    echo "\$ stb-kgrid --file $struct --density 0.2 --save-report --no-intro"
    (cd "$OUT/$case" && stb-kgrid --file "$struct" --density 0.2 --save-report --no-intro > console.log)
    grep -E "Dimensionality|Suggested Monkhorst-Pack" "$OUT/$case/console.log"
    pause
done

echo "=================================================================="
echo " Side by side"
echo "=================================================================="
for case in silicon molecule chain graphene; do
    printf "  %-9s " "$case"
    grep "Suggested Monkhorst-Pack" "$OUT/$case/console.log"
done
cat <<'EOF'

  silicon  (3D): N N N   -- fully periodic          -> real mesh on all 3 axes
  molecule (0D): 1 1 1   -- no periodic axis         -> single k-point (Gamma only)
  chain    (1D): N 1 1   -- periodic along a only    -> real mesh on a, 1 elsewhere
  graphene (2D): N N 1   -- periodic along a and b   -> real mesh in-plane, 1 out-of-plane
EOF
echo
echo "=================================================================="
echo " The --vacuum-gap threshold is a real, tunable choice"
echo "=================================================================="
cat <<'EOF'
Detecting "vacuum" isn't magic -- it's a threshold (default 10 Ang): the
largest empty gap along an axis, wrapped periodically, must exceed it to be
treated as non-periodic. A structure with a 12 Ang gap along c sits just
above the default threshold (-> 2D) but below a stricter 15 Ang one
(-> 3D, a real if sparse periodic direction) -- same gap, different
classification, different k-grid:
EOF
TMP="$(mktemp -d)"
cat > "$TMP/gap_boundary.fdf" << 'EOF'
NumberOfSpecies    1
NumberofAtoms      2

%block ChemicalSpeciesLabel
 1  14  Si
%endblock ChemicalSpeciesLabel

LatticeConstant 1.0 Ang

AtomicCoordinatesFormat  Fractional

%block LatticeVectors
 5.430   0.000   0.000
 0.000   5.430   0.000
 0.000   0.000   24.000
%endblock LatticeVectors

%block AtomicCoordinatesAndAtomicSpecies
 0.00 0.00 0.00 1
 0.00 0.00 0.50 1
%endblock AtomicCoordinatesAndAtomicSpecies
EOF
echo
echo "\$ stb-kgrid --file gap_boundary.fdf --density 0.2 --no-intro   (default --vacuum-gap 10)"
(cd "$TMP" && stb-kgrid --file gap_boundary.fdf --density 0.2 --no-intro) | grep -E "Dimensionality|Suggested Monkhorst-Pack"
echo
echo "\$ stb-kgrid --file gap_boundary.fdf --density 0.2 --vacuum-gap 15 --no-intro"
(cd "$TMP" && stb-kgrid --file gap_boundary.fdf --density 0.2 --vacuum-gap 15 --no-intro) | grep -E "Dimensionality|Suggested Monkhorst-Pack"
rm -rf "$TMP"
pause

echo "=================================================================="
echo " Proof: CLI and the interactive stb-suite menu agree"
echo "=================================================================="
echo "Driving the same silicon/density-0.2 case non-interactively (via a"
echo "piped printf) and diffing the grid line against output/silicon/ --"
echo "no separate output folder is kept for this, it's just a check."
TMP="$(mktemp -d)"
cp structure.fdf "$TMP/"
echo
echo "\$ printf '1.2\\nstructure.fdf\\n0.2\\n\\n0\\n' | stb-suite"
(cd "$TMP" && printf '1.2\nstructure.fdf\n0.2\n\n0\n' | stb-suite > session.log 2>&1)
GUIDE_LINE=$(grep -n "K-Point Density Recommendation Guide" "$TMP/session.log" | head -1 | cut -d: -f1)
PROMPT_LINE=$(grep -n "K-point density (e.g., 0.2)" "$TMP/session.log" | head -1 | cut -d: -f1)
if [ -n "$GUIDE_LINE" ] && [ -n "$PROMPT_LINE" ] && [ "$GUIDE_LINE" -lt "$PROMPT_LINE" ]; then
    echo "Confirmed: the density guide (line $GUIDE_LINE) prints before the density prompt (line $PROMPT_LINE)."
else
    echo "Unexpected: guide/prompt ordering looks wrong -- see $TMP/session.log."
fi
CLI_GRID=$(grep "Suggested Monkhorst-Pack" "$OUT/silicon/console.log")
MENU_GRID=$(grep "Suggested Monkhorst-Pack" "$TMP/session.log")
if [ "$CLI_GRID" = "$MENU_GRID" ]; then
    echo "Confirmed: identical grid ($MENU_GRID) from the CLI and the interactive menu."
else
    echo "Unexpected: grids differ -- CLI: $CLI_GRID / menu: $MENU_GRID"
fi
rm -rf "$TMP"
pause

echo "=================================================================="
echo " Done"
echo "=================================================================="
cat <<EOF
Four self-contained folders were generated under output/:
  silicon/   molecule/   chain/   graphene/

Each has stb_kgrid_report.txt (the exact console report) and references.bib
(SIESTA + Monkhorst-Pack citations) -- everything you need to justify the
k-grid choice in a paper or a lab notebook.

As a next step, try on your own:
  stb-kgrid --file structure.fdf --density 0.05   # high precision, expensive
  stb-kgrid --file structure.fdf --density 0.5    # low precision, cheap

In a real workflow you rarely need to run stb-kgrid by hand at all --
stb-inputfile (example 1.1) already calls the exact same core/kspace.py
logic internally and writes the resulting kgrid.MonkhorstPack line straight
into calc.fdf for you.
EOF
