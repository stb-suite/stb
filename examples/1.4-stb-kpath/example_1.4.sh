#!/bin/bash
# Guided example: stb-kpath (code 1.4 in the stb-suite menu)
#
# Not an automated test (see test/1-inputs/4-k-path/test.sh for that) -- a
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

echo "=================================================================="
echo " Why band structures need a PATH, not a grid"
echo "=================================================================="
cat <<'EOF'
stb-kgrid (example 1.2) builds a Monkhorst-Pack GRID -- a set of k-points
spread uniformly over the whole Brillouin zone (BZ), used to approximate an
INTEGRAL (total energy, charge density). A band structure plot is a
different job: it shows how the energy E(k) varies ALONG a route through
the BZ, so it's traced along a PATH connecting the BZ's high-symmetry
points instead -- E(k) is smooth, and the physically interesting features
(band extrema, gaps, crossings) are captured well by a handful of
well-chosen points and the segments between them.

Which points are "high-symmetry", and in what order to connect them, is
NOT arbitrary -- it depends on the shape of the Brillouin zone, which in
turn depends on the crystal's Bravais lattice. Setyawan & Curtarolo (2010)
-- the exact reference stb-kpath writes to references.bib -- standardized
this into one path convention per Bravais lattice type, now used across
essentially every DFT code for reproducible band-structure plots. ASE
implements that exact convention (Cell.bandpath) and extends it to 1D/2D
systems via a periodic-axes mask; stb-kpath is a thin, dimension-aware
wrapper around it -- it does not invent or choose paths itself.
EOF
pause

echo "=================================================================="
echo " Bravais lattice vs. space group -- two different questions"
echo "=================================================================="
cat <<'EOF'
stb-kpath reports BOTH, and they answer different questions:

  Bravais lattice  -- from the raw LATTICE VECTORS alone (ignoring the
                      atoms). Determines the BZ shape and which path
                      template applies. This is what the k-path is built
                      from.
  Space group      -- the FULL crystallographic symmetry, lattice + atomic
                      basis together (via the same core/symmetry.py
                      accessors stb-inputfile also uses).

These can legitimately disagree: a structure file can list atoms in a
non-primitive (conventional) cell, so the raw lattice vectors alone look
lower-symmetry than the real crystal actually is once its atomic
arrangement is taken into account. You'll see exactly this below with
bulk silicon.
EOF
pause

echo "=================================================================="
echo " Two ways to run it"
echo "=================================================================="
cat <<'EOF'
A -- direct CLI:
  stb-kpath --file structure.fdf

B -- interactive stb-suite menu:
  $ stb-suite
  Select an option (0-6, or a tool code like 4.1.2): 1.4

"1.4" asks for the structure file, then -- same "advanced settings" pattern
used elsewhere in the suite -- offers to configure vacuum-gap, the
Bravais-lattice tolerance, and the symmetry tolerances, all gated behind a
single [y/N] question so the essential flow stays short.
EOF
pause

echo "=================================================================="
echo " output/silicon/  --  bulk silicon (3D)"
echo "=================================================================="
mkdir -p "$OUT/silicon"
cp structure.fdf "$OUT/silicon/"
echo "\$ stb-kpath --file structure.fdf --save-report --no-intro"
(cd "$OUT/silicon" && stb-kpath --file structure.fdf --save-report --no-intro > console.log)
grep -E "Bravais lattice|Space group|Crystal system|^Path :" "$OUT/silicon/console.log"
echo
echo "Notice the mismatch: the Bravais lattice (CUB) comes from the plain"
echo "cubic lattice vectors written in structure.fdf; the space group"
echo "(Fd-3m, diamond structure) only emerges once the 8-atom basis is"
echo "considered. The k-path itself follows the Bravais lattice (CUB), since"
echo "that's what actually defines this Brillouin zone's shape."
pause

for case in chain graphene; do
    case "$case" in
        chain)    struct="structure_chain.fdf";    label="chain/ -- carbon chain (1D)" ;;
        graphene) struct="structure_graphene.fdf"; label="graphene/ -- graphene monolayer (2D)" ;;
    esac
    echo "=================================================================="
    echo " output/$case  --  $label"
    echo "=================================================================="
    mkdir -p "$OUT/$case"
    cp "$struct" "$OUT/$case/"
    echo "\$ stb-kpath --file $struct --save-report --no-intro"
    (cd "$OUT/$case" && stb-kpath --file "$struct" --save-report --no-intro > console.log)
    grep -E "Dimensionality|Bravais lattice|Space group|may not reflect|^Path :" "$OUT/$case/console.log"
    pause
done

echo "=================================================================="
echo " output/molecule/  --  CH4, isolated (0D) -- a DELIBERATE hard error"
echo "=================================================================="
cat <<'EOF'
Bloch's theorem (see example 1.2) only applies along a periodic direction.
An isolated molecule has none -- there is no k to speak of, so a k-path is
not just low-quality here, it's not physically meaningful at all.
stb-kpath treats this as a real error (exit code 1), not a soft warning:
EOF
mkdir -p "$OUT/molecule"
cp structure_molecule.fdf "$OUT/molecule/"
echo
echo "\$ stb-kpath --file structure_molecule.fdf --no-intro"
set +e
(cd "$OUT/molecule" && stb-kpath --file structure_molecule.fdf --no-intro > console.log 2>&1)
MOL_EXIT=$?
set -e
grep -E "Dimensionality|ERROR" "$OUT/molecule/console.log"
if [ "$MOL_EXIT" -eq 1 ] && [ ! -s "$OUT/molecule/kpath_bs.fdf" ]; then
    echo "Confirmed: exit code 1, and no kpath_bs.fdf was written."
else
    echo "Unexpected: exit code was $MOL_EXIT (expected 1), or a file was written anyway."
fi
pause

echo "=================================================================="
echo " The --vacuum-gap threshold changes the physics, not just a number"
echo "=================================================================="
cat <<'EOF'
Same idea as example 1.2: a structure with a 12 Ang gap along c sits just
above the default 10 Ang vacuum threshold (-> 2D, square BZ, path
M-Gamma-X-M) but below a stricter 15 Ang one (-> 3D, tetragonal BZ, a
completely different, longer path) -- same atoms, different threshold,
different Brillouin zone.
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
echo "\$ stb-kpath --file gap_boundary.fdf --no-intro   (default --vacuum-gap 10)"
(cd "$TMP" && stb-kpath --file gap_boundary.fdf --no-intro) | grep -E "Dimensionality|Bravais lattice|^Path :"
echo
echo "\$ stb-kpath --file gap_boundary.fdf --vacuum-gap 15 --no-intro"
(cd "$TMP" && stb-kpath --file gap_boundary.fdf --vacuum-gap 15 --no-intro) | grep -E "Dimensionality|Bravais lattice|^Path :"
rm -rf "$TMP"
pause

echo "=================================================================="
echo " Proof: CLI and the interactive stb-suite menu agree"
echo "=================================================================="
echo "Driving the same silicon case through the advanced-settings gate"
echo "(vacuum-gap/eps/symprec/angle all at their defaults, entered explicitly)"
echo "and diffing the path line against output/silicon/ -- no separate output"
echo "folder is kept for this, it's just a check."
TMP="$(mktemp -d)"
cp structure.fdf "$TMP/"
echo
echo "\$ printf '1.4\\nstructure.fdf\\ny\\n10.0\\n0.0002\\n0.001\\n5.0\\n\\n0\\n' | stb-suite"
(cd "$TMP" && printf '1.4\nstructure.fdf\ny\n10.0\n0.0002\n0.001\n5.0\n\n0\n' | stb-suite > session.log 2>&1)
CLI_PATH=$(grep "^Path :" "$OUT/silicon/console.log")
MENU_PATH=$(grep "^Path :" "$TMP/session.log")
if [ "$CLI_PATH" = "$MENU_PATH" ]; then
    echo "Confirmed: identical path ($MENU_PATH) from the CLI and the interactive menu."
else
    echo "Unexpected: paths differ -- CLI: $CLI_PATH / menu: $MENU_PATH"
fi
rm -rf "$TMP"
pause

echo "=================================================================="
echo " Done"
echo "=================================================================="
cat <<EOF
Three self-contained folders were generated under output/ (silicon/,
chain/, graphene/), each with kpath_bs.fdf, stb_kpath_report.txt, and
references.bib (SIESTA + Setyawan-Curtarolo citations). output/molecule/
deliberately has none of those -- 0D is a hard error.

As a next step, try on your own:
  stb-kpath --file structure.fdf --vacuum-gap 5    # a tighter vacuum threshold
  stb-kpath --file structure_chain.fdf -o my_path.fdf

In a real workflow, stb-inputfile's "bands" mode (example 1.1) already
calls stb-kpath itself to produce the kpath_bs.fdf its calc.fdf needs --
you rarely have to run this one by hand except to inspect or double-check
the path it chose.
EOF
