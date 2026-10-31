#!/bin/bash
# Guided example: 2.19 -- Native Charges (stb-nativecharges; code 2.19 in the stb-suite menu).
#
# Not an automated test (see test/2-analysis/19-nativecharges/test.sh for that) -- a commented walk-through: it runs real
# commands, one case at a time, on a COPY of runs/ inside output/, and shows the piece of the report that proves what
# just happened. runs/ itself is never touched. Pauses between cases (`yes '' | bash example_2.19.sh` runs it straight
# through). Safe to re-run -- it always starts by wiping its own output/.
#
# runs/ holds real SIESTA 5.4 runs whose logs carry the Mulliken, Hirshfeld and Voronoi charges: water with four
# basis sizes, rock-salt NaCl and the O2 triplet. No SIESTA is needed. README section 6 explains every number.

set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"
OUT="$DIR/output"
rm -rf "$OUT"
mkdir -p "$OUT"
cp -a "$DIR/runs/." "$OUT/"

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

section() {   # $1 = file, $2 = a section title, exactly
    strip_colors < "$1" | awk -v t="$2" '/^\[[0-9]+\] [A-Z]/ {if (on) exit; sub(/^\[[0-9]+\] /, "", $0); if ($0 == t) {on = 1; next}} on {print}' \
        | grep -v "^-----"
}

# nc <log name> <args...>: runs stb-nativecharges inside output/, the report saved as output/<log name>
nc() {
    local log=$1; shift
    echo "  \$ stb-nativecharges $*"
    (cd "$OUT" && stb-nativecharges --no-intro "$@" < /dev/null 2>&1 | strip_colors > "$OUT/$log") || true
}
warnings() { section "$OUT/$1" "WARNINGS" | grep -E "^  [0-9]+\." | fold -s -w 116 | sed 's/^/  /'; }


banner "Welcome: the atomic charges SIESTA already printed, compared"
cat <<'EOF'
With WriteMullikenPop 1, Charge.Hirshfeld end and Charge.Voronoi end, SIESTA ends its log with three tables of atomic
charges, each from its own way of cutting one electron density into atoms. stb-nativecharges reads them, puts them
side by side, checks them, and shows how much they disagree -- which is the first thing to know about a partial
charge.

The tour:
   1  water (DZP): the whole report          4  NaCl: an ionic solid
   2  the dipole of point charges            5  the O2 triplet: spin
   3  four basis sets: Mulliken's sign flip  6  a log without the blocks; the menu
EOF
pause


banner "Case 1 -- Water, DZP: the whole report"
nc case1.txt --path h2o_basis/DZP
section "$OUT/case1.txt" "RUN" | sed 's/^/  /'
section "$OUT/case1.txt" "CHARGES PER ATOM" | sed 's/^/  /'
section "$OUT/case1.txt" "CHECKS" | sed 's/^/  /'
warnings case1.txt
cat <<'EOF'
  -> Hirshfeld and Voronoi make O negative (-0.24, -0.19), as chemistry expects; Mulliken makes it POSITIVE (+0.12).
     Every method adds up to zero and gives population + q = Z_val; they differ by 0.36 e on the same atom.
EOF
pause


banner "Case 2 -- The dipole of point charges"
section "$OUT/case1.txt" "DIPOLE" | sed 's/^/  /'
cat <<'EOF'
  -> Charges placed on the nuclei carry only 23-48 % of the dipole SIESTA computes from the density: the rest is in
     the shape of each atom's charge (the O lone pairs). Mulliken's even points the wrong way.
EOF
pause


banner "Case 3 -- Four basis sets: Mulliken's sign flip"
nc case3.txt --path 'h2o_basis/*' --csv h2o_basis_charges.csv
section "$OUT/case3.txt" "CHARGES PER SPECIES" | sed 's/^/  /'
cat <<'EOF'
  -> Hirshfeld barely moves with the basis (O: -0.29 to -0.24), Voronoi a little more; Mulliken goes from -0.65 (SZ)
     to +0.18 (TZP). Mulliken splits each overlap population half and half between two atoms, so diffuse and
     polarization orbitals on O that describe H's density count as O's (figures/basis_dependence.png).
EOF
pause


banner "Case 4 -- NaCl: an ionic solid"
nc case4.txt --path nacl
section "$OUT/case4.txt" "CHARGES PER ATOM" | sed 's/^/  /'
cat <<'EOF'
  -> A textbook ionic solid gets +-0.08 (Mulliken), +-0.14 (Hirshfeld), +-0.17 (Voronoi) -- far from +-1. These
     partitions share overlapping density between neighbours; Hirshfeld-I (5.17, its example) gives +-0.51.
     Compare charges only within one method and one basis.
EOF
pause


banner "Case 5 -- The O2 triplet: spin"
nc case5.txt --path o2
section "$OUT/case5.txt" "CHARGES PER ATOM" | sed 's/^/  /'
section "$OUT/case5.txt" "CHECKS" | grep -E "Hirsh|Expected" | sed 's/^/  /'
cat <<'EOF'
  -> A homonuclear molecule has no charge to split; the spin is shared exactly, 1 per atom, and the Sz of every
     method add up to the spin moment SIESTA prints, 2.
EOF
pause


banner "Case 6 -- A log without the blocks; the interactive menu"
mkdir -p "$OUT/nocharges"; cp "$OUT/nacl/nacl.XV" "$OUT/nocharges/"
sed '/^Charge\./d; /^WriteMullikenPop/d' "$OUT/nacl/calc.fdf" > "$OUT/nocharges/calc.fdf"
python3 - "$OUT/nacl/calc.out" "$OUT/nocharges/calc.out" <<'PY'
import sys, re
text = open(sys.argv[1]).read()
for h in ("Mulliken", "Hirshfeld", "Voronoi"):
    text = re.sub(rf"\n{h} Atomic Populations:\n.*?\n Total[^\n]*\n", "\n", text, flags=re.S)
open(sys.argv[2], "w").write(text)
PY
nc case6.txt --path nocharges
section "$OUT/case6.txt" "RUN" | grep -E "Mulliken|Hirshfeld|Voronoi" | sed 's/^/  /'
# path glob, no label, details n, no csv, no report, Enter, 0 = exit
(cd "$OUT" && printf '2.19\nh2o_basis/*\n\nn\n\nn\n\n0\n' | stb-suite 2>&1 | strip_colors > "$OUT/case6_menu.txt") || true
echo "  the menu (stb-suite, 2.19) on h2o_basis/*:"; grep -E "^h2o_basis/DZP .*O " "$OUT/case6_menu.txt" | sed 's/^/    /'
cat <<'EOF'
  -> Each missing table names the fdf flag that makes SIESTA print it. The menu runs the same command.
EOF
pause


banner "Done"
cat <<EOF
Everything is in $OUT/. What to take from the tour:
  - A partial charge depends on the partition: compare within one method (and one basis), never across.
  - Mulliken follows the basis, even in sign; Hirshfeld and Voronoi are stable but small for ionic solids.
  - Point charges miss the atoms' own dipoles: the molecular dipole comes from the density.
EOF
