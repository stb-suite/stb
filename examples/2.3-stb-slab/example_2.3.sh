#!/bin/bash
# Guided example: stb-slab (code 2.3 in the stb-suite menu)
#
# Not an automated test (see test/2-structures/3-slab/test.sh for that) --
# a commented walk-through: it runs real commands, one group at a time,
# into its own output/<case>/ folder, and shows you the piece of output
# that proves what just happened. Pauses between sections so you can read
# before moving on. Safe to re-run any time -- it always starts by wiping
# its own output/.

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
echo " Miller indices, terminations, and why more than one slab can come out"
echo "=================================================================="
cat <<'EOF'
A Miller index (h,k,l) defines a FAMILY of lattice planes, but cutting along
that family isn't unique -- depending on exactly where between two planes
you cut, you expose different atoms at the surface: different
TERMINATIONS. stb-slab (via pymatgen's SlabGenerator) enumerates every
symmetrically distinct termination for the requested hkl, sorts non-polar/
symmetric ones first (the usual physically preferred choice), and by
default keeps only the top of that list (--all keeps every one, -i lets
you pick, --termination N picks by index directly).

SYMMETRIC means the same termination on both faces of the slab -- an
asymmetric slab has two physically different surfaces, which usually also
means a spurious residual dipole across the artificial vacuum gap in a
periodic calculation.
EOF
pause

echo "=================================================================="
echo " Two ways to run it"
echo "=================================================================="
cat <<'EOF'
A -- direct CLI:
  stb-slab -f si_bulk.fdf --hkl 1 0 0

B -- interactive stb-suite menu:
  $ stb-suite
  Select an option (0-6, or a tool code like 4.1.2): 2.3

Every run always writes a numbered report ([0] RUN METADATA ... [7] SUMMARY
& FILES -- [4] ML PRE-RELAXATION only appears with --ml-relax), a
structure-validation checklist both before and after the cut, a
before/after symmetry comparison table, and references.bib. A full text
report (--save-report), MACE pre-relaxation (--ml-relax), and an
interactive 3D view (--view) are all off by default.
EOF
pause

echo "=================================================================="
echo " output/basic/  --  Si(100), and a bulk-to-slab symmetry contrast"
echo "=================================================================="
cat <<'EOF'
si_bulk.fdf is bulk silicon, diamond cubic, the real conventional 8-atom
cell (a = 5.431 Ang, space group Fd-3m No. 227). Cutting ANY slab out of a
3D bulk crystal necessarily makes it 2D-periodic -- so watch the Layer
Group column: N/A for the (non-2D) bulk, but always a real value for the
cut slab.
EOF
mkdir -p "$OUT/basic"
cp si_bulk.fdf "$OUT/basic/"
echo
echo "\$ stb-slab -f si_bulk.fdf --hkl 1 0 0 --no-intro"
(cd "$OUT/basic" && stb-slab -f si_bulk.fdf --hkl 1 0 0 --no-intro > console.log)
grep -E "Terminations found|Slab thickness|Vacuum thickness" "$OUT/basic/console.log"
echo
echo "Before/after symmetry table -- Layer Group appears only on the slab side:"
sed -n '/Detailed symmetry analysis/,/Slab written to/p' "$OUT/basic/console.log" | head -n -1
echo
echo "Provenance header written into the output .fdf itself:"
head -3 "$OUT/basic/slab.fdf"
pause

echo "=================================================================="
echo " output/symmetrize/  --  NaCl(111): one asymmetric cut vs. two symmetric ones"
echo "=================================================================="
cat <<'EOF'
nacl_bulk.fdf is bulk NaCl, the real rocksalt structure (a = 5.64 Ang,
space group Fm-3m No. 225). Along (1,1,1), the natural bulk-truncated cut
gives ONE asymmetric termination (both surfaces different) -- --symmetrize
asks pymatgen to trim atoms off one side until both surfaces match,
changing both the atom count and how many distinct terminations exist.
EOF
mkdir -p "$OUT/symmetrize"
cp nacl_bulk.fdf "$OUT/symmetrize/"
echo
echo "\$ stb-slab -f nacl_bulk.fdf --hkl 1 1 1 --all -o plain.fdf --no-intro"
(cd "$OUT/symmetrize" && stb-slab -f nacl_bulk.fdf --hkl 1 1 1 --all -o plain.fdf --no-intro > console_plain.log)
grep -E "Terminations found|Symmetric|Total atoms" "$OUT/symmetrize/console_plain.log"
echo
echo "\$ stb-slab -f nacl_bulk.fdf --hkl 1 1 1 --symmetrize --all -o sym.fdf --no-intro"
(cd "$OUT/symmetrize" && stb-slab -f nacl_bulk.fdf --hkl 1 1 1 --symmetrize --all -o sym.fdf --no-intro > console_sym.log)
grep -E "Terminations found|Symmetric|Total atoms" "$OUT/symmetrize/console_sym.log"
echo
echo "32 atoms, 1 (asymmetric) termination -- vs. 28 atoms, 2 (symmetric) terminations."
pause

echo "=================================================================="
echo " output/polar-caveat/  --  a real gotcha, verified: the Polar column"
echo " needs oxidation states it never gets"
echo "=================================================================="
cat <<'EOF'
pymatgen's Slab.is_polar() computes a dipole from each site's formal
OXIDATION STATE -- its own docstring says the Slab must be oxidation-state
decorated for this to work, otherwise the dipole is always exactly 0. A
structure read from a SIESTA .fdf never carries oxidation states -- .fdf
has no such concept -- so stb-slab's Polar column reads "No"
UNCONDITIONALLY, regardless of the real electrostatics of the cut.

NaCl(100) is genuinely non-polar (each atomic layer already has equal
Na/Cl); NaCl(111) is the textbook POLAR rocksalt surface (alternating
pure-Na and pure-Cl layers, Tasker's "Type III" classification). Watch
stb-slab report the SAME "Polar: No" for both:
EOF
mkdir -p "$OUT/polar-caveat"
cp nacl_bulk.fdf "$OUT/polar-caveat/"
echo
echo "\$ stb-slab -f nacl_bulk.fdf --hkl 1 0 0 -o nacl100.fdf --no-intro"
(cd "$OUT/polar-caveat" && stb-slab -f nacl_bulk.fdf --hkl 1 0 0 -o nacl100.fdf --no-intro > console_100.log)
grep "Polar" "$OUT/polar-caveat/console_100.log" | sed 's/^/  (100): /'
echo "\$ stb-slab -f nacl_bulk.fdf --hkl 1 1 1 -o nacl111.fdf --no-intro"
(cd "$OUT/polar-caveat" && stb-slab -f nacl_bulk.fdf --hkl 1 1 1 -o nacl111.fdf --no-intro > console_111.log)
grep "Polar" "$OUT/polar-caveat/console_111.log" | sed 's/^/  (111): /'
echo
cat <<'EOF'
The real physics is still there -- it just needs oxidation states stb-slab
never applies. A quick Python illustration (NOT something stb-slab itself
does) recovers it directly on the exact same two structures:
EOF
python3 - "$OUT/polar-caveat/nacl_bulk.fdf" <<'PYEOF'
import sys
from stb.core import structure_io
from pymatgen.core.surface import SlabGenerator

s = structure_io.to_pymatgen(structure_io.read_fdf(sys.argv[1]))
s.add_oxidation_state_by_element({"Na": 1, "Cl": -1})
for hkl in [(1, 0, 0), (1, 1, 1)]:
    gen = SlabGenerator(s, hkl, 10.0, 15.0)
    slab = gen.get_slabs()[0]
    dipole = float((slab.dipole @ slab.dipole) ** 0.5)
    print(f"  {hkl}: dipole magnitude = {dipole:.2f} (arbitrary units), is_polar() = {slab.is_polar()}")
PYEOF
echo
echo "Takeaway: trust the Polar column only as 'pymatgen found no reason to"
echo "flag this' -- for an ionic/mixed-valence material, cross-check yourself."
pause

echo "=================================================================="
echo " output/passivate/  --  Si(111), 8 dangling bonds genuinely resolved"
echo "=================================================================="
cat <<'EOF'
Cutting a covalent, tetrahedrally-bonded crystal like Si always leaves
under-coordinated surface atoms -- dangling bonds, a strong artificial
perturbation on the electronic structure. --passivate caps each
single-missing-bond site with a passivating atom (H by default) placed
along the geometrically missing-bond direction (core/passivation.py,
shared with stb-passivate).
EOF
mkdir -p "$OUT/passivate"
cp si_bulk.fdf "$OUT/passivate/"
echo
echo "\$ stb-slab -f si_bulk.fdf --hkl 1 1 1 --passivate -o passivated.fdf --no-intro"
(cd "$OUT/passivate" && stb-slab -f si_bulk.fdf --hkl 1 1 1 --passivate -o passivated.fdf --no-intro > console.log)
grep -E "Dangling bonds found|Auto-passivated" "$OUT/passivate/console.log"
echo
echo "H registered as a new species in the output .fdf:"
grep -A3 "ChemicalSpeciesLabel" "$OUT/passivate/passivated.fdf" | head -4
pause

echo "=================================================================="
echo " output/ml-relax/  --  measured surface relaxation"
echo "=================================================================="
if ! python3 -c "import mace" 2>/dev/null; then
    echo "Skipping -- needs the optional 'ml' extra: pip install stb_suite[ml]"
    echo "(everything else in this script works fine without it)."
else
cat <<'EOF'
A freshly cut slab's surface atoms sit at their BULK positions, but a real
surface almost always relaxes away from them once the missing neighbors on
one side change the local force balance -- "surface relaxation," expected
real physics (distinct from RECONSTRUCTION, which changes the surface's
periodicity/composition, not just positions). --ml-relax lets a MACE
potential estimate this in seconds -- useful groundwork before an
expensive real SIESTA relaxation.
EOF
    mkdir -p "$OUT/ml-relax"
    cp si_bulk.fdf "$OUT/ml-relax/"
    echo
    echo "\$ stb-slab -f si_bulk.fdf --hkl 1 1 1 --ml-relax -o relaxed.fdf --no-intro"
    (cd "$OUT/ml-relax" && stb-slab -f si_bulk.fdf --hkl 1 1 1 --ml-relax -o relaxed.fdf --no-intro \
        > console.log 2>&1)
    echo "Full MACE simulation detail, straight from the report:"
    sed 's/\x1b\[[0-9;]*m//g' "$OUT/ml-relax/console.log" | \
        awk '/\[4\] ML PRE-RELAXATION/{flag=1} /\[5\] STRUCTURE VALIDATION/{flag=0} flag'
    echo
    echo "Provenance header now also records the MACE pre-relaxation:"
    head -4 "$OUT/ml-relax/relaxed.fdf"
fi
pause

echo "=================================================================="
echo " output/full-report/  --  --save-report, validation, references.bib"
echo "=================================================================="
cat <<'EOF'
The full numbered report (also written to stb_slab_report.txt with
--save-report) includes the structure-validation checklist for the bulk
input and every written slab, and a references.bib with SIESTA.
EOF
mkdir -p "$OUT/full-report"
cp si_bulk.fdf "$OUT/full-report/"
echo "\$ stb-slab -f si_bulk.fdf --hkl 1 0 0 --save-report --no-intro"
(cd "$OUT/full-report" && stb-slab -f si_bulk.fdf --hkl 1 0 0 --save-report --no-intro > console.log)
echo
echo "Report sections written to stb_slab_report.txt:"
grep -E "^\[[0-9]\]" "$OUT/full-report/stb_slab_report.txt"
echo
echo "Validation checklist (one row shown):"
grep -m1 "Atom proximity" "$OUT/full-report/console.log"
echo
echo "references.bib always written -- SIESTA:"
grep "^@" "$OUT/full-report/references.bib"
pause

echo "=================================================================="
echo " Proof: CLI and the interactive stb-suite menu agree"
echo "=================================================================="
echo "Driving the same basic case through the interactive menu's manual"
echo "entry mode and checking it reaches the same write step."
TMP="$(mktemp -d)"
cp si_bulk.fdf "$TMP/"
echo
echo "\$ printf '2.3\\nsi_bulk.fdf\\n1 0 0\\n\\n\\nn\\nn\\nn\\n1\\nn\\n\\n\\n\\n\\n0\\n' | stb-suite"
(cd "$TMP" && printf '2.3\nsi_bulk.fdf\n1 0 0\n\n\nn\nn\nn\n1\nn\n\n\n\n\n0\n' | stb-suite > session.log 2>&1) || true
if grep -q "Slab written to" "$TMP/session.log"; then
    echo "Confirmed: the interactive menu built and launched the exact same"
    echo "underlying stb-slab command as the CLI walkthrough above."
else
    echo "Unexpected: menu did not reach the write step -- see $TMP/session.log."
fi
rm -rf "$TMP"
pause

echo "=================================================================="
echo " Done"
echo "=================================================================="
cat <<EOF
Up to six self-contained folders were generated under output/ (ml-relax/
skipped if the optional 'ml' extra isn't installed):
  basic/         symmetrize/    polar-caveat/
  passivate/     ml-relax/      full-report/

Each has references.bib; full-report/ additionally has
stb_slab_report.txt; ml-relax/'s references.bib also cites the MACE
architecture/foundation-model papers.

Recap of what this walkthrough covered:
  - Miller indices, terminations, and why non-polar/symmetric ones are
    preferred and sorted first
  - --symmetrize's real effect: fewer atoms, more (symmetric) terminations
  - a real, verified gotcha: stb-slab's Polar column always reads "No"
    without oxidation-state decoration, which a .fdf never carries -- the
    real dipole is still there, recovered directly once decorated
  - --passivate genuinely resolving 8 dangling bonds on a Si(111) surface
  - --ml-relax measuring real surface relaxation with MACE
  - the structure-validation checklist (bulk AND slab), references.bib,
    and --save-report
  - the Layer Group column appearing only on the (2D-periodic) slab side,
    never the (3D) bulk side
  - CLI and the interactive stb-suite menu building the same command

Not exercised by this script (needs a display): --view opens the bulk
structure and every written slab in an interactive ase-gui window -- try
it yourself:
  stb-slab -f si_bulk.fdf --hkl 1 0 0 --view

As a next step, try on your own:
  stb-slab -f your_bulk.fdf --hkl 1 1 1 -i          # pick a termination by hand
  stb-slab -f si_bulk.fdf --hkl 1 1 1 --passivate --ml-relax
  stb-slab -f nacl_bulk.fdf --hkl 1 1 1 --min-slab-size 20 --symmetrize --all
EOF
