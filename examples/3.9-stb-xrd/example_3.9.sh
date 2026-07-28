#!/bin/bash
# Guided example: stb-xrd (code 3.9 in the stb-suite menu)
#
# Not an automated test (see test/3-analysis/9-xrd/test.sh for that -- it
# uses a synthetic NaCl fixture instead, to exercise every code path
# cheaply). This walkthrough uses a real structure instead: the same CrS
# 2D monolayer geometry already used by examples/3.7-stb-workfunction (a
# real, converged SIESTA calculation, SystemLabel "siesta", originally
# fetched from the twodmatpedia OPTIMADE database via stb-fetch).
#
# siesta.XV is copied from examples/3.7-stb-workfunction/. stb-xrd only
# accepts .fdf/.STRUCT_OUT structure files, not .XV, so this script
# converts it first -- via sisl directly (NOT stb-translate's own "siesta"
# input reader: that one was found, while writing this example, to
# mis-handle this real file, treating its Cartesian Bohr positions as
# fractional and erroring out downstream -- a separate, pre-existing bug
# in translate.py, out of scope for this tool's own migration, reported
# instead of silently worked around). sisl itself (via
# core.deps.read_sisl_geometry_xv_or_fdf) already reads this exact file
# correctly -- it's the same mechanism stb-workfunction's own vacuum-axis
# detection already uses successfully on it.

set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

OUT="$DIR/output"
rm -rf "$OUT"
mkdir -p "$OUT"

# Non-interactive matplotlib backend for the "Two ways to run it" proof
# below (which navigates the interactive menu) -- --view isn't exercised
# in this script (see its own section), but this keeps any accidental
# plot from blocking the walkthrough.
export MPLBACKEND=Agg

pause() {
    echo
    read -p "  [Press Enter to continue] " -r
    echo
}

echo "=================================================================="
echo " Recovering the real CrS structure from siesta.XV"
echo "=================================================================="
cat <<'EOF'
stb-xrd needs a .fdf or .STRUCT_OUT structure file. siesta.XV (the CrS
geometry) is converted to crs_structure.fdf once, up front, via sisl
directly (see the comment at the top of this script for why not
stb-translate):
EOF
python3 - <<'PYEOF'
import sisl
from stb.core.structure_io import FdfStructure, write_fdf

geom = sisl.get_sile("siesta.XV").read_geometry()
z_to_symbol = {24: "Cr", 16: "S"}
species = []
for z in geom.atoms.Z:
    sym = z_to_symbol[int(z)]
    if sym not in species:
        species.append(sym)
atoms = [(z_to_symbol[int(z)], pos) for z, pos in zip(geom.atoms.Z, geom.fxyz)]
species_meta = {sym: {"id": str(i + 1), "Z": [z for z, s in z_to_symbol.items() if s == sym][0]}
                for i, sym in enumerate(species)}
structure = FdfStructure(lattice=geom.cell, lattice_constant=1.0, species=species,
                          species_meta=species_meta, atoms=atoms, coord_format="fractional")
write_fdf(structure, "crs_structure.fdf", header_comment=[
    "Real CrS 2D monolayer geometry, recovered from siesta.XV (a finished,",
    "converged SIESTA calculation originally fetched from the twodmatpedia",
    "OPTIMADE database via stb-fetch) via sisl -- the stb-xrd (3.9) example input.",
])
print(f"Wrote crs_structure.fdf: {len(atoms)} atoms, species {species}")
PYEOF
pause

echo "=================================================================="
echo " What stb-xrd computes, and why"
echo "=================================================================="
cat <<'EOF'
stb-xrd simulates a powder X-ray diffraction pattern from a structure:
every set of lattice planes (hkl) that satisfies Bragg's law,

    n*lambda = 2*d_hkl*sin(theta)

produces a diffraction peak at the corresponding 2*theta angle, with an
intensity set by the structure factor (how strongly that plane's atoms
scatter in phase). A space group's symmetry can force some structure
factors to vanish exactly -- "systematic absences" -- so the resulting
peak list is also a fingerprint of the crystal's symmetry, not just its
cell size.

See the README for the full theory (Bragg's law, structure factors and
extinction, wavelength choice, and what the experimental-comparison
similarity score measures) and an important limitation specific to THIS
structure: a 2D slab modeled with a large vacuum gap along one axis
still gets treated as fully 3D-periodic, so its huge, artificial c-axis
periodicity produces spurious low-angle (00l) peaks that don't correspond
to any real bulk diffraction pattern -- worth knowing before trusting a
low-angle peak from a vacuum-padded structure.

Every run prints a numbered [0]...[6] report -- the full peak table is
NOT part of it anymore (it used to be printed in full, which could run to
hundreds of lines for a low-symmetry structure over a wide 2-theta range);
[2] DIFFRACTION PATTERN instead prints just a compact summary (peak
count, strongest peak, resolution). --save-report persists the report to
stb_xrd_report.txt -- off by default. --save-gnuplot writes the actual
peak list to xrd_pattern.dat, now with a complete, human-readable header
(structure/formula/space group/wavelength/range, not just pyxtal's own
terse one-liner) plus a matching stick-pattern xrd_pattern.gplot -- also
off by default (this tool used to write the data file UNCONDITIONALLY on
every run; that's no longer the case). --top now controls how many peaks
land in that saved file (it used to just trim the console table, which
no longer exists). --view shows a matplotlib preview (renamed from the
old --plot, same off-by-default behavior) -- and, unlike before, actually
stays open until you close it: the old --plot delegated to a non-blocking
call (fig.show(), or pyxtal's own plot_pxrd(), which uses the same
non-blocking call internally) that made the window disappear as soon as
the script exited.
EOF
pause

echo "=================================================================="
echo " output/basic/  --  a real structure (CrS monolayer)"
echo "=================================================================="
cat <<'EOF'
Watch [1] report the real detected space group/crystal system/lattice
(P4/nmm, tetragonal -- CrS's actual symmetry) and [2] the pattern summary
-- peak count, strongest peak, and the resolution (smallest d-spacing)
reached over the scanned range (the full peak-by-peak list only ever
goes to xrd_pattern.dat, with --save-gnuplot -- see the next section):
EOF
mkdir -p "$OUT/basic"
cp crs_structure.fdf "$OUT/basic/"
echo
echo "\$ stb-xrd --file crs_structure.fdf --format fdf --no-intro"
(cd "$OUT/basic" && stb-xrd --file crs_structure.fdf --format fdf \
    --no-intro > console.log 2>&1)
awk '/\[1\] STRUCTURE/{flag=1} /\[4\] OUTPUT/{flag=0} flag' "$OUT/basic/console.log" | head -30
pause

echo "=================================================================="
echo " output/save-gnuplot/  --  --save-gnuplot (off by default)"
echo "=================================================================="
cat <<'EOF'
--save-gnuplot writes xrd_pattern.dat -- every peak (2theta, d, h, k, l,
intensity), with a complete header (structure/formula/space group/
wavelength/range) instead of just pyxtal's own terse one-liner -- and a
matching xrd_pattern.gplot stick-pattern plot script, ready to render
with a real gnuplot install:
EOF
mkdir -p "$OUT/save-gnuplot"
cp crs_structure.fdf "$OUT/save-gnuplot/"
echo
echo "\$ stb-xrd --file crs_structure.fdf --format fdf --save-gnuplot --no-intro"
(cd "$OUT/save-gnuplot" && stb-xrd --file crs_structure.fdf --format fdf \
    --save-gnuplot --no-intro > console.log 2>&1)
echo "xrd_pattern.dat header + first peaks:"
head -5 "$OUT/save-gnuplot/xrd_pattern.dat"
echo
echo "xrd_pattern.gplot:"
cat "$OUT/save-gnuplot/xrd_pattern.gplot"
if command -v gnuplot > /dev/null; then
    (cd "$OUT/save-gnuplot" && gnuplot xrd_pattern.gplot)
    echo "(rendered xrd_pattern.pdf with the real, installed gnuplot)"
fi
pause

echo "=================================================================="
echo " output/compare/  --  --compare-to an experimental pattern"
echo "=================================================================="
cat <<'EOF'
mock_experimental.dat is a SYNTHETIC "experimental" pattern -- built by
adding small angle jitter and intensity noise to this structure's own
simulated peaks (see the README), not independent lab data. Enough to
exercise the real comparison code path: pyxtal's Similarity() interpolates
both patterns onto a common grid and reports a cosine-weighted score
(0-1, higher is more similar):
EOF
mkdir -p "$OUT/compare"
cp crs_structure.fdf mock_experimental.dat "$OUT/compare/"
echo
echo "\$ stb-xrd --file crs_structure.fdf --format fdf --compare-to mock_experimental.dat --no-intro"
(cd "$OUT/compare" && stb-xrd --file crs_structure.fdf --format fdf \
    --compare-to mock_experimental.dat --no-intro > console.log 2>&1)
awk '/\[3\] EXPERIMENTAL/{flag=1} /\[4\] OUTPUT/{flag=0} flag' "$OUT/compare/console.log"
pause

echo "=================================================================="
echo " output/full-report/  --  --save-report (off by default)"
echo "=================================================================="
cat <<'EOF'
Every run always prints the numbered [0]...[6] report to the console.
--save-report additionally persists it to stb_xrd_report.txt -- off by
default, so a plain run only ever writes references.bib, no text report
file (and no data/gnuplot files either, without --save-gnuplot):
EOF
mkdir -p "$OUT/full-report"
cp crs_structure.fdf "$OUT/full-report/"
echo
echo "\$ stb-xrd --file crs_structure.fdf --format fdf --no-intro   # default: no report, no data/gnuplot"
(cd "$OUT/full-report" && stb-xrd --file crs_structure.fdf --format fdf \
    --no-intro > console_default.log 2>&1)
ls "$OUT/full-report/" | grep -v "^crs_structure.fdf$\|console"
echo "(no xrd_pattern.dat/.gplot, no stb_xrd_report.txt -- only references.bib)"
echo
echo "\$ stb-xrd --file crs_structure.fdf --format fdf --save-report --no-intro"
(cd "$OUT/full-report" && stb-xrd --file crs_structure.fdf --format fdf \
    --save-report --no-intro > console_saved.log 2>&1)
echo "Report sections written to stb_xrd_report.txt:"
grep -E "^\[[0-9]+\]" "$OUT/full-report/stb_xrd_report.txt"
echo
echo "references.bib -- SIESTA + pyxtal (the diffraction simulation engine):"
grep "^@" "$OUT/full-report/references.bib"
pause

echo "=================================================================="
echo " Two ways to run it"
echo "=================================================================="
cat <<'EOF'
A -- direct CLI:
  stb-xrd --file crs_structure.fdf --format fdf

B -- interactive stb-suite menu:
  $ stb-suite
  Select an option (0-6, or a tool code like 4.1.2): 3.9

Both paths call the exact same underlying tool -- proven directly below.
EOF
TMP="$(mktemp -d)"
cp crs_structure.fdf "$TMP/"
echo
echo "\$ printf '3.9\\ncrs_structure.fdf\\n1\\n\\n\\n\\n\\n\\nn\\nn\\nn\\n' | stb-suite     # format fdf, defaults, no top/compare-to/save-report/gnuplot/view"
(cd "$TMP" && printf '3.9\ncrs_structure.fdf\n1\n\n\n\n\n\nn\nn\nn\n' | stb-suite > session.log 2>&1) || true
CLI_LINE=$(grep "Space Group" "$OUT/basic/console.log" | head -1)
MENU_LINE=$(grep "Space Group" "$TMP/session.log" | head -1)
if [ "$CLI_LINE" = "$MENU_LINE" ]; then
    echo "Confirmed: identical space-group line from the CLI and the interactive menu."
    echo "  $CLI_LINE"
else
    echo "Unexpected: results differ -- see $TMP/session.log."
fi
rm -rf "$TMP"
pause

echo "=================================================================="
echo " --view (needs a display)"
echo "=================================================================="
cat <<'EOF'
Not exercised interactively by this script (needs a display, though it IS
covered headlessly by test/3-analysis/9-xrd/test.sh under MPLBACKEND=Agg,
which is exactly what caught the old fig.show()/plot_pxrd() bug not being
fixed would NOT have caught, since a non-blocking call "succeeds" either
way). --view shows the simulated stick pattern (or, with --compare-to, an
overlay of both patterns) and now genuinely blocks until you close it:

  stb-xrd --file crs_structure.fdf --format fdf --view
  stb-xrd --file crs_structure.fdf --format fdf --compare-to mock_experimental.dat --view
EOF
pause

echo "=================================================================="
echo " Done"
echo "=================================================================="
cat <<EOF
Four self-contained folders were generated under output/:
  basic/          save-gnuplot/   compare/        full-report/

Each has references.bib; save-gnuplot/ additionally has
xrd_pattern.dat/xrd_pattern.gplot, and full-report/ has
stb_xrd_report.txt (only from its --save-report run).

Recap of what this walkthrough covered:
  - Bragg's law and why symmetry can force some peaks to vanish
    (systematic absences) -- CrS's own P4/nmm space group
  - the new [1] STRUCTURE section: space group/crystal system/lattice
    info that this tool didn't report at all before
  - the peak table moved OUT of the console/report entirely -- [2] is now
    just a compact summary, and the full list only lives in
    xrd_pattern.dat (--save-gnuplot), now with a complete header instead
    of pyxtal's own terse one-liner
  - --top now trims that saved file, not a console table
  - --compare-to's similarity score against an experimental pattern
  - the fixed --view: a real, verified bug (the plot window disappearing
    immediately) traced to pyxtal's own plot_pxrd(), not just this tool's
    code, and fixed by building the preview independently
  - the numbered [0]...[6] report, --save-report, references.bib
    (now including pyxtal's own citation, not just SIESTA's)
  - CLI and the interactive stb-suite menu building the same command

As a next step, try on your own with any structure:
  stb-xrd --file my_structure.fdf --format fdf --save-gnuplot --view
  stb-xrd --file my_structure.fdf --format fdf --compare-to my_experimental.dat --view
EOF
