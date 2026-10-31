#!/bin/bash
# Guided example: 3.3 -- Grid to Cube Converter (stb-cube; code 3.3 in the stb-suite menu).
#
# Not an automated test (see test/3-utils/3-cube/test.sh for that) -- a commented walk-through: it runs real
# commands, one case at a time, on COPIES of siesta_runs/ inside output/, and shows the piece of the report that
# proves what just happened. siesta_runs/ itself is never touched. Pauses between cases (`yes '' | bash
# example_3.3.sh` runs it straight through). Safe to re-run -- it always starts by wiping its own output/.
#
# siesta_runs/ holds three real, finished SIESTA 5.4 runs: bulk Si with every grid SIESTA can write, a water molecule
# in a 9 Ang box, and the spin-polarized O2 triplet. No SIESTA is needed. README section 6 explains every number.

set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"
OUT="$DIR/output"
rm -rf "$OUT"
mkdir -p "$OUT"
cp -a "$DIR/siesta_runs/si" "$DIR/siesta_runs/h2o" "$DIR/siesta_runs/o2" "$OUT/"

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

# sc <folder> <log name> <args...>: runs stb-cube inside output/<folder>/, the screen report saved as <log name>
sc() {
    local dir=$1 log=$2; shift 2
    echo "  \$ stb-cube $*"
    (cd "$OUT/$dir" && stb-cube --no-intro "$@" < /dev/null 2>&1 | strip_colors > "$OUT/$log") || true
}
warnings() { section "$OUT/$1" "WARNINGS" | grep -E "^  [0-9]+\." | fold -s -w 116 | sed 's/^/  /'; }


banner "Welcome: SIESTA's real-space grids as cube files"
cat <<'EOF'
SIESTA writes densities and potentials on its real-space mesh (<label>.RHO, .DRHO, .VT, .VH, ...) in a binary
format of its own. stb-cube writes them as Gaussian cube files, which VESTA, VMD, XCrySDen and the Bader codes read --
and its report says what the grid is: unit, mesh and cutoff, the geometry it belongs to, what it integrates to and
what it must integrate to, spin channels, the dipole of a molecule, and a read-back check of the cube.

The tour:
   1  bulk Si: the valence density, the whole report     5  spin: the O2 triplet
   2  units: what a naive conversion writes              6  what the report catches
   3  every grid of one run                              7  labels and folders
   4  a molecule in a box: dipole and vacuum level       8  the interactive menu
EOF
pause


banner "Case 1 -- Bulk Si: the valence density"
sc si case1.txt -l si -t RHO
section "$OUT/case1.txt" "INPUT GRID" | grep -E "Mesh|spacing|Voxel|Effective|angles" | sed 's/^/  /'
section "$OUT/case1.txt" "SIESTA RUN" | grep -E "Electrons|Mesh" | sed 's/^/  /'
section "$OUT/case1.txt" "VALUES" | sed 's/^/  /'
section "$OUT/case1.txt" "OUTPUT" | grep -E "File|Read-back|agree|Integral" | sed 's/^/  /'
echo "  the first lines of the cube:"; head -9 "$OUT/si/si_RHO.cube" | sed 's/^/    /'
cat <<'EOF'
  -> Mesh.Cutoff 200 Ry became 248.9 Ry: SIESTA rounds the mesh up to 30 points per vector (even, factors 2, 3, 5);
     the report computes it from the grid and agrees with the log. The grid integrates to the 8 valence electrons, and
     so does the cube with its own Bohr voxels. The minimum is AT the nucleus (a valence pseudo-density is small in
     the core); the maximum is on the Si-Si bond.
EOF
pause


banner "Case 2 -- Units: what a naive conversion writes"
(cd "$OUT/si" && python3 - <<'PY'
import sisl
from ase.io.cube import read_cube_data
g = sisl.get_sile("si.RHO").read_grid()
g.geometry = sisl.get_sile("si.XV").read_geometry()
g.write("naive.cube")                                     # sisl's own one-line conversion
for name in ("naive.cube", "si_RHO.cube"):
    data, atoms = read_cube_data(name)
    dv_bohr = atoms.get_volume() / data.size / 0.529177210903**3
    print(f"    {name:12s}  max {data.max():8.5f}   sum x voxel (Bohr^3) = {data.sum() * dv_bohr:8.4f} e")
PY
)
cat <<'EOF'
  -> A cube's lengths are in Bohr and its values carry no unit. sisl reads densities as e/Ang^3; written next to Bohr
     lengths they are 1/0.529177^3 = 6.748 times too large -- 54 electrons instead of 8. Every viewer shows the same
     picture at a higher level, so nothing looks wrong. stb-cube writes e/Bohr^3 and checks the integral.
EOF
pause


banner "Case 3 -- Every grid of one run"
printf "  %-6s %-14s %-48s\n" "grid" "integral (e)" "expected"
for t in RHO DRHO BADER TOCH IOCH; do
    sc si case3_$t.txt -l si -t $t --dry-run > /dev/null
    i=$(grep -m1 "^Integral " "$OUT/case3_$t.txt" | awk '{print $3}')
    e=$(grep -m1 "^Expected " "$OUT/case3_$t.txt" | sed 's/^Expected *: //; s/ -- difference.*//')
    printf "  %-6s %-14s %-48s\n" "$t" "$i" "$e"
done
grep "Charge moved" "$OUT/case3_DRHO.txt" | sed 's/^/  DRHO: /'
for t in VT VH VNA; do
    sc si case3_$t.txt -l si -t $t --dry-run > /dev/null
    echo "  $t: $(grep -E '^(minimum|maximum|mean) ' "$OUT/case3_$t.txt" | awk -F'|' '{gsub(/ /,"",$1); gsub(/ /,"",$2); printf "%s %s eV   ", $1, $2}')"
done
cat <<'EOF'
  -> Each grid must integrate to something known: RHO the 8 valence electrons, DRHO zero (bonding moved 0.46 e),
     BADER 28 = 8 + 20 core electrons (its maximum back on the nucleus, as a Bader analysis needs), TOCH zero and
     IOCH -8 up to the mesh's sampling of the ionic Gaussians. VH and VNA share their cell average; in a bulk cell
     that average is a convention, not a vacuum level.
EOF
pause


banner "Case 4 -- A molecule in a box: water"
sc h2o case4_rho.txt -l h2o -t RHO
sc h2o case4_drho.txt -l h2o -t DRHO
sc h2o case4_vh.txt -l h2o -t VH --potential-unit Ha
section "$OUT/case4_rho.txt" "GEOMETRY" | grep "Dimensionality" | sed 's/^/  /'
section "$OUT/case4_rho.txt" "VALUES" | grep -E "^Integral|^Dipole|^component|along|^\|p\|" | grep -v "^Vacuum" | sed 's/^/  /'
grep "Charge moved" "$OUT/case4_drho.txt" | sed 's/^/  DRHO: /'
section "$OUT/case4_vh.txt" "VALUES" | sed -n '/Vacuum planes/,$p' | grep -v "^---" | fold -s -w 116 | sed 's/^/  /'
cat <<'EOF'
  -> The dipole computed from the grid agrees with SIESTA's to 1e-5 D: grid, geometry and log belong together.
     The potential is flat in the vacuum along x and y, but slopes along z, the dipole axis: periodic images of a
     polar molecule set up a field across the box, so there is no single vacuum level (figures/h2o.png).
EOF
pause


banner "Case 5 -- Spin: the O2 triplet"
sc o2 case5.txt -l o2 -t RHO --spin diff
section "$OUT/case5.txt" "SPIN" | grep -v "^---" | sed 's/^/  /'
ls "$OUT/o2"/*.cube | sed "s#$OUT/#  written: output/#"
cat <<'EOF'
  -> Up 7, down 5, up - down 2 -- SIESTA's spin moment. The integral of |up - down| is a little larger (2.033):
     some density is polarized against the moment. --spin diff writes the magnetization density
     (figures/o2_spin.png); for a potential, 'total' is the average of up and down, not their sum.
EOF
pause


banner "Case 6 -- What the report catches"
sc si case6_missing.txt -l si -t LDOS
grep FAIL "$OUT/case6_missing.txt" | fold -s -w 116 | sed 's/^/  /'
mkdir -p "$OUT/stale" && cp "$OUT/si/si.RHO" "$OUT/si/calc.out" "$OUT/stale/" && cp "$OUT/o2/o2.XV" "$OUT/stale/si.XV"
sc stale case6_stale.txt -l si -t RHO --dry-run
echo "  a .XV from another run:"; warnings case6_stale.txt
sc si case6_spin.txt -l si -t RHO --spin diff --dry-run
echo "  --spin on an unpolarized grid:"; warnings case6_spin.txt
cat <<'EOF'
  -> A grid that was not saved stops with the list of grids that are there and the fdf flag that writes the one
     asked for. A geometry whose cell differs from the grid's is flagged (its atoms would land in the wrong places);
     a spin channel that does not exist falls back to the total, with a warning.
EOF
pause


banner "Case 7 -- Labels and folders"
sc . case7_path.txt -l si/si -t VT --potential-unit Ry --save-report
grep -E "^Input|^Log|^File  .*cube|^Report" "$OUT/case7_path.txt" | sed 's/^/  /'
sc h2o case7_auto.txt -t DRHO --dry-run
grep "^Input" "$OUT/case7_auto.txt" | sed 's/^/  /'
cat <<'EOF'
  -> -l may include a folder: the log, the .XV, the cube and the report are all taken from / written next to the
     grid. Without -l, the only file of that type in the current folder is used.
EOF
pause


banner "Case 8 -- The interactive menu (stb-suite, code 3.3)"
# type 10 = VH, label h2o, spin 1 = total, unit Ha, default output, no dry run, no report, Enter, 0 = exit
(cd "$OUT/h2o" && printf '3.3\n10\nh2o\n1\nHa\n\nn\nn\n\n0\n' | stb-suite 2>&1 | strip_colors > "$OUT/case8.txt") || true
grep -E "^Grid type  |^Potential unit|^Output|^Read-back" "$OUT/case8.txt" | sed 's/^/  /'
cat <<'EOF'
  -> The menu asks for the grid type (a numbered list with each type's fdf flag), the label, the spin channel, the
     potential unit, the output name, dry run and report, and runs the same command.
EOF
pause


banner "Done"
cat <<EOF
Everything is in $OUT/. What to take from the tour:
  - A cube carries no units: stb-cube writes lengths in Bohr, densities in e/Bohr^3, and says so in the file.
  - Read the integral line: each grid type integrates to something known (8, 0, 28, -8 for this Si).
  - For a molecule or slab, the dipole from the grid against SIESTA's proves that grid, geometry and log belong
    together; a sloped vacuum potential means there is no single vacuum level.
EOF
