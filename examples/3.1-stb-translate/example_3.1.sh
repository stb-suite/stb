#!/bin/bash
# Guided example: 3.1 -- File Translator (stb-translate; code 3.1 in the stb-suite menu).
#
# Not an automated test (see test/3-utils/1-translator/test.sh for that) -- a commented walk-through: it runs real
# commands, one case at a time, into its own output/<case>/ folder, and shows the piece of output that proves what
# just happened. Pauses between cases so you can read before moving on (`yes '' | bash example_3.1.sh` runs it
# straight through). Safe to re-run any time -- it always starts by wiping its own output/. Every case takes seconds;
# nothing here needs a GPU, a model or SIESTA.
#
# The structures are real: bulk Si and a Si(111) slab at the lattice constant of the 7.1 SIESTA training set, and a
# SIESTA relaxation (siesta.STRUCT_OUT) of a Si adatom on that slab from the 7.6 example; make_inputs.py writes the
# hand-made ones (a CIF from a database, a cell in Bohr, a molecule, noisy and broken files). README section 6 explains
# every number.

set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"
export MPLBACKEND=Agg
OUT="$DIR/output"
rm -rf "$OUT"
mkdir -p "$OUT"

pause() {
    echo
    read -p "  [Press Enter to continue] " -r || true
    echo
}

banner() {
    echo "=================================================================="
    echo " $1"
    echo "=================================================================="
}

strip_colors() { sed 's/\x1b\[[0-9;]*m//g'; }

section() {   # $1 = file, $2 = start of a section title: from "[N] TITLE" to the next "[M] "
    strip_colors < "$1" | awk -v t="$2" '/^\[[0-9]+\] [A-Z]/ {if (on) exit; if (index($0, "] " t) > 0) on = 1} on {print}'
}

# st <case> <args...>: runs stb-translate in output/<case>/, the screen report saved as log.txt
st() {
    local case=$1; shift
    mkdir -p "$OUT/$case"
    echo "  \$ stb-translate $*"
    (cd "$OUT/$case" && stb-translate --no-intro "$@" < /dev/null 2>&1 | strip_colors > log.txt) || true
}
log() { echo "$OUT/$1/log.txt"; }


banner "Welcome: one structure, many file formats"
cat <<'EOF'
stb-translate converts a structure between the file formats around SIESTA -- fdf, STRUCT_OUT, POSCAR, CIF, XYZ, XSF,
FHI-aims geometry.in, DFTB+ .gen and LAMMPS data -- and reports what it read (composition, cell, density,
periodicity, shortest bonds, symmetry) and whether the written file reproduces it (it reads the file back).

The tour:
   1  the inputs                                   7  a CIF from a database: the asymmetric unit
   2  a SIESTA relaxation to CIF, formats auto     8  a molecule: XYZ has no cell
   3  one crystal, every format                    9  what the report catches: mix-ups, wrong scale
   4  fractional and Cartesian, by hand           10  batch conversion
   5  units: a lattice in Bohr                    11  guard rails
   6  symmetry and tolerance; CIF symmetrization  12  the interactive menu (stb-suite, code 3.1)
Theory: README section 1. Numbers: README section 6.
EOF
pause


banner "Case 1 -- The inputs"
python3 make_inputs.py 2>/dev/null | sed 's/^/  /'
ls -1 *.fdf *.vasp *.cif *.xyz *.STRUCT_OUT | sed 's/^/    /'
cat <<'EOF'
  si_bulk.fdf / si111_2x2.fdf: bulk Si and a Si(111) slab (SIESTA input format);
  si111_adatom_relaxed.STRUCT_OUT: what a SIESTA relaxation writes (cell + fractional positions, nothing else).
EOF
pause


banner "Case 2 -- A SIESTA relaxation to CIF: both formats from the file names"
st case2 -i "$DIR/si111_adatom_relaxed.STRUCT_OUT" -o si111_adatom.cif
cat "$(log case2)"
cat <<'EOF'
  -> No -if/-of: '.STRUCT_OUT' says SIESTA output, '.cif' says CIF. The report reads the structure (25 Si, a 2D slab
     with vacuum along c, P3m1), writes the CIF symmetrized at 0.01 Ang in the standard setting -- the cell is
     re-expressed, so the check is "the same crystal" (StructureMatcher) instead of atom by atom.
EOF
pause


banner "Case 3 -- One crystal, every output format"
for f in si.cif si.vasp si.xyz si.xsf geometry.in si.gen si.data si_copy.fdf; do
    st "case3_${f//./_}" -i "$DIR/si_bulk.fdf" -o "$f"
done
printf "  %-12s %-8s %-46s %s\n" "file" "format" "coordinates written" "read back"
for f in si.cif si.vasp si.xyz si.xsf geometry.in si.gen si.data si_copy.fdf; do
    L=$(log "case3_${f//./_}")
    fmt=$(grep -m1 "^Output            :" "$L" | sed 's/.* -- \([a-z]*\) .*/\1/')
    co=$(grep "^Coordinates       :" "$L" | tail -1 | cut -c21- | sed 's/ (as the input.*//; s/ *$//')
    dv=$(grep "largest atom displacement" "$L" | awk -F'|' '{print $3}' | tr -d ' ')
    printf "  %-12s %-8s %-46s %s\n" "$f" "$fmt" "$co" "max shift ${dv:-?} Ang"
done
cat <<'EOF'
  -> Eight formats, each read back to the same 8 atoms within 1e-15 Ang. What differs is what each format can say:
     CIF and XSF and XYZ fix the coordinate type, the others follow the input (-cf to choose).
EOF
pause


banner "Case 4 -- Fractional and Cartesian coordinates, by hand"
st case4_frac -i "$DIR/si111_2x2.fdf" -o slab_frac.fdf -cf direct
st case4_cart -i "$DIR/si111_2x2.fdf" -o slab_cart.fdf -cf cartesian
echo "  The 5th atom of the oblique slab cell (gamma = 60 deg), as each file stores it:"
echo "    fractional : $(sed -n '/AtomicCoordinatesAndAtomicSpecies/,/endblock/p' $OUT/case4_frac/slab_frac.fdf | sed -n 6p)"
echo "    Cartesian  : $(sed -n '/AtomicCoordinatesAndAtomicSpecies/,/endblock/p' $OUT/case4_cart/slab_cart.fdf | sed -n 6p)"
python3 - <<'EOF'
import numpy as np
A = np.array([[7.75022690, 0, 0], [3.87511345, 6.71189338, 0], [0, 0, 22.11903799]])   # rows: a, b, c
s = np.array([0.16666667, 0.16666667, 0.37483566])
print(f"    by hand    : r = u a + v b + w c = {s @ A}  (the rows of the cell matrix, weighted by u, v, w)")
print(f"    and back   : s = r A^-1 = {(s @ A) @ np.linalg.inv(A)}")
EOF
pause


banner "Case 5 -- Units: a lattice given in Bohr"
grep -A4 "LatticeConstant" si_bohr.fdf | sed 's/^/    /'
st case5 -i "$DIR/si_bohr.fdf" -o si_primitive.vasp
echo "  The POSCAR written from it (line 2 is VASP's scale factor):"
sed -n 1,9p "$OUT/case5/si_primitive.vasp" | sed 's/^/    /'
grep -E "^Cell   |^Volume|0.01          " "$(log case5)" | sed 's/^/  /'
cat <<'EOF'
  -> 'LatticeConstant 1.0 Bohr' becomes the POSCAR scale factor 0.529177 Ang: the vectors stay as written, the cell is
     the Si primitive cell (a = 3.875 Ang, 60 deg), 20.57 Ang^3 per atom -- the same as si_bulk.fdf's 8-atom cell.
EOF
pause


banner "Case 6 -- Symmetry and tolerance; how a CIF is symmetrized"
for a in 0.002 0.02 0.2; do
    st "case6_$a" -i "$DIR/si_noisy_$a.fdf" --dry-run
    echo "  si_noisy_$a.fdf (atoms displaced by up to $a Ang):"
    section "$(log case6_$a)" "SYMMETRY" | grep -E "^0\.|^Primitive" | sed 's/^/    /'
done
st case6_cif -i "$DIR/si_noisy_0.002.fdf" -o noisy_symmetrized.cif
st case6_p1 -i "$DIR/si_noisy_0.002.fdf" -o noisy_p1.cif --cif-symprec 0
echo "  Written as CIF (default --cif-symprec 0.01, then 0 = P1):"
grep -E "largest atom displacement|^Symmetrized" "$(log case6_cif)" | sed 's/^/    /'
grep -E "largest atom displacement" "$(log case6_p1)" | sed 's/^/    /'
cat <<'EOF'
  -> A tolerance decides what 'symmetric' means: noise of 0.002 Ang is Fd-3m at 0.01 Ang, noise of 0.02 Ang needs
     0.1 Ang, and 0.2 Ang is P1 at any reasonable tolerance (README Figure 3). A CIF written with a tolerance puts the
     atoms ON their symmetric positions (here they moved 0.0026 Ang); --cif-symprec 0 keeps them exactly as read.
EOF
pause


banner "Case 7 -- A CIF from a database: the asymmetric unit and its symmetry operations"
grep -E "_symmetry_space_group_name|_cell_length_a" nacl_asym.cif | sed 's/^/    /'
echo "    atom sites listed in the CIF:"
sed -n '/_atom_site_occupancy/,$p' nacl_asym.cif | tail -n +2 | sed 's/^/      /'
st case7 -i "$DIR/nacl_asym.cif" -o nacl.fdf
grep -E "^Formula|^Atoms|Cl-Na|^0.01  |^Primitive" "$(log case7)" | sed 's/^/  /'
cat <<'EOF'
  -> 2 sites in the CIF, 8 atoms in the fdf: the reader applied the 192 operations of Fm-3m (4 copies each, the
     'multiplicity' column). SIESTA needs every atom of the cell; the primitive cell has 2.
EOF
pause


banner "Case 8 -- A molecule: an XYZ file has no cell"
st case8_nobox -i "$DIR/water.xyz" -o water.fdf
grep -E "error|required" "$(log case8_nobox)" | sed 's/^/  /'
st case8 -i "$DIR/water.xyz" -o water.fdf --lattice-vectors 15 0 0 0 15 0 0 0 15
grep -E "^Formula|^Periodicity|H-O |^0.01  " "$(log case8)" | sed 's/^/  /'
cat <<'EOF'
  -> --lattice-vectors puts the molecule in a 15 Ang box; the report sees vacuum along a, b and c (a 0D system) and the
     point-group symmetry of H2O (mm2, written as the space group Pmm2 of the box).
EOF
pause


banner "Case 9 -- What the report catches"
st case9_mixup -i "$DIR/si_mixup.vasp" --dry-run
st case9_scale -i "$DIR/si_wrongscale.vasp" --dry-run
echo "  si_mixup.vasp (Cartesian numbers under 'Direct'):"
grep -E "^Coordinate range|^0.01  " "$(log case9_mixup)" | sed 's/^/    /'
grep -E "^  [0-9]+\." "$(log case9_mixup)" | sed 's/^/  /'
echo "  si_wrongscale.vasp (scale factor 0.5):"
grep -E "^Volume|Si-Si" "$(log case9_scale)" | sed 's/^/    /'
grep -E "^  [0-9]+\." "$(log case9_scale)" | sed 's/^/  /'
cat <<'EOF'
  -> Both files READ without error -- that is the danger. The mix-up becomes an R3m crystal with impossible fractional
     coordinates up to 4.1; the wrong scale a Si of 18 g/cm^3 with 1.19 Ang bonds. The report says so.
EOF
pause


banner "Case 10 -- Batch conversion"
st case10 -i "$DIR/si*.fdf" -of cif --out-dir cif
section "$(log case10)" "SUMMARY" | sed 's/^/  /'
pause


banner "Case 11 -- Guard rails"
g() { echo "  \$ stb-translate $*"; (cd "$OUT" && stb-translate --no-intro "$@" < /dev/null 2>&1 | strip_colors \
        | grep -E "error|FAIL" | head -2 | sed 's/^/      /') || true; }
g -i "$DIR/si_bulk.fdf" -o out.unknown
g -i "$DIR/si_bulk.fdf"
g -i "$DIR/nope.fdf" -o out.cif
g -i "$DIR/water.xyz" -o w.fdf
g -i "$DIR/si*.fdf" -o one.cif
cat <<'EOF'
  -> An unknown output extension, no output file, a missing input, an XYZ without a cell, and one output name for
     many inputs: each stops with a clear message.
EOF
pause


banner "Case 12 -- The interactive menu (stb-suite, code 3.1)"
mkdir -p "$OUT/case12"
printf '3.1\n%s\n0\n%s\n0\n3\nn\ny\n\n0\n' "$DIR/si111_adatom_relaxed.STRUCT_OUT" "$OUT/case12/menu.vasp" \
    | (cd "$OUT/case12" && stb-suite > log.txt 2>&1) || true
grep -E "Auto-detected|Output format:" "$OUT/case12/log.txt" | strip_colors | sed 's/^/  /'
ls "$OUT/case12" | sed 's/^/    /'
pause


banner "Done"
cat <<EOF
Everything is in $OUT/. What to take from the tour:
  - A structure file is a cell (three vectors x a lattice constant, in some unit), species and positions in fractional
    or Cartesian coordinates; every format stores a subset of that, in its own convention.
  - Let the names choose the formats, and read the report: the coordinate range, the density and the shortest bonds
    catch the classic mix-ups before SIESTA does.
  - Symmetry depends on the tolerance; a CIF written with one is symmetrized (--cif-symprec 0 keeps the data as is).
  - Every written file is read back: an unexpected shift means a lossy format or a bug (it found one in the .gen reader).
EOF
