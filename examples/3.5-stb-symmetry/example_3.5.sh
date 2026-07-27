#!/bin/bash
# Guided example: stb-symmetry (code 3.5 in the stb-suite menu)
#
# Not an automated test (see test/3-analysis/5-symmetry/test.sh for that)
# -- a commented walk-through: it runs real commands, one group at a time,
# into its own output/<case>/ folder, and shows you the piece of output
# that proves what just happened. Pauses between sections so you can read
# before moving on. Safe to re-run any time -- it always starts by wiping
# its own output/.
#
# Fixtures used, all copied from test/3-analysis/5-symmetry/:
#   nacl.fdf          -- textbook rock-salt NaCl (space group Fm-3m, No. 225)
#   nacl_noisy.fdf     -- the same NaCl with a small random atomic displacement
#   graphene_slab.fdf  -- a 2-atom graphene layer, 20 Ang vacuum along c
#   molecule.fdf       -- an isolated water molecule, 10 Ang vacuum on a/b/c
#   wire.fdf           -- a 1D chain, 10 Ang vacuum on 2 axes (a, b)

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
echo " What stb-symmetry computes, and why"
echo "=================================================================="
cat <<'EOF'
stb-symmetry reads a SIESTA structure file (.fdf or a post-relaxation
.STRUCT_OUT) and reports, with no SIESTA re-run needed:
  - the 3D space group, Wyckoff positions/orbits, and symmetry operations
    (via spglib, wrapped by pymatgen's SpacegroupAnalyzer)
  - for a vacuum-padded 2D slab: the LAYER GROUP instead (spglib's own
    dedicated 2D detection -- the 3D space group is not physically
    meaningful once one axis is just empty space)
  - for a vacuum-padded isolated molecule: the POINT GROUP instead
    (pymatgen's separate, non-periodic detector)
  - a per-atom "distortion" (Å): how far each atom sits from where the
    detected symmetry operations would place it exactly
  - optionally, a tolerance-sensitivity scan (--scan-symprec), a
    comparison between two structures (--compare-to), and a
    symmetry-refined structure (--write-refined)

See the README for the underlying theory of each of these and their
documented limitations (no rod-group detection for a wire, the distortion
metric's dependence on --symprec, etc.).

Every run prints a numbered [0]...[11] report, the same style every
newer tool in this suite uses. --save-report additionally persists that
report to stb_symmetry_report.txt -- off by default, so a plain run only
ever writes references.bib, no text report file. The tool used to always
write symmetry.dat unconditionally; that's gone now (a stale one from an
older run is cleaned up automatically instead).
EOF
pause

echo "=================================================================="
echo " output/bulk-nacl/  --  space group, Wyckoff orbits, operations (NaCl)"
echo "=================================================================="
cat <<'EOF'
NaCl rock salt: a textbook Fm-3m (No. 225) structure. Watch [2] SPACE GROUP
report the full symmetry classification, [6] collapse the 8 atoms into
just 2 Wyckoff orbits (4a for Na, 4b for Cl -- symmetry-equivalent atoms
share one orbit instead of being listed as 8 unrelated sites), and [7]
show 0.0000 distortion for every atom (a perfectly symmetric structure):
EOF
mkdir -p "$OUT/bulk-nacl"
cp nacl.fdf "$OUT/bulk-nacl/"
echo
echo "\$ stb-symmetry --file nacl.fdf --format fdf --no-intro"
(cd "$OUT/bulk-nacl" && stb-symmetry --file nacl.fdf --format fdf \
    --no-intro > console.log 2>&1)
awk '/\[2\] SPACE GROUP/{flag=1} /\[4\] LATTICE/{flag=0} flag' "$OUT/bulk-nacl/console.log"
pause

echo "--write-refined writes the refined structure INSIDE --output-dir (a plain"
echo "filename, not a path of its own) with a header documenting how it was"
echo "generated -- self-describing even opened on its own, no console report needed:"
echo
echo "\$ stb-symmetry --file nacl.fdf --format fdf --write-refined refined.fdf --no-intro"
(cd "$OUT/bulk-nacl" && stb-symmetry --file nacl.fdf --format fdf \
    --write-refined refined.fdf --no-intro > console_refined.log 2>&1)
echo "output/bulk-nacl/refined.fdf's header:"
head -5 "$OUT/bulk-nacl/refined.fdf"
pause

echo "=================================================================="
echo " output/scan-symprec/  --  symmetry hidden by numerical noise"
echo "=================================================================="
cat <<'EOF'
nacl_noisy.fdf is the same NaCl with every atom nudged by a small random
displacement (simulating post-DFT-relaxation numerical drift). At the
default --symprec (1e-3 Ang) this is detected as P1 (no symmetry at all)
-- but --scan-symprec sweeps a range of tolerances and shows the real
Fm-3m symmetry is still there, just hidden below a looser tolerance:
EOF
mkdir -p "$OUT/scan-symprec"
cp nacl_noisy.fdf "$OUT/scan-symprec/"
echo
echo "\$ stb-symmetry --file nacl_noisy.fdf --format fdf --scan-symprec --no-intro"
(cd "$OUT/scan-symprec" && stb-symmetry --file nacl_noisy.fdf --format fdf \
    --scan-symprec --no-intro > console.log 2>&1)
awk '/\[5\] TOLERANCE/{flag=1} /\[6\]/{flag=0} flag' "$OUT/scan-symprec/console.log"
pause

echo "=================================================================="
echo " output/layer-group/  --  a 2D slab needs LAYER GROUP, not SPACE GROUP"
echo "=================================================================="
cat <<'EOF'
graphene_slab.fdf has 20 Ang of vacuum along c. spglib's 3D space-group
detection treats that vacuum as just an unusually tall periodic cell --
it reports a real group (P6/mmm, No. 191), but one that depends on the
arbitrary vacuum thickness, not the physical 2D symmetry of the sheet
itself. [3] LAYER GROUP is spglib's own dedicated 2D detection instead:
the correct, vacuum-thickness-independent classification (p6/mmm, No. 80,
Wyckoff 2b, site symmetry -6m2 -- textbook graphene crystallography):
EOF
mkdir -p "$OUT/layer-group"
cp graphene_slab.fdf "$OUT/layer-group/"
echo
echo "\$ stb-symmetry --file graphene_slab.fdf --format fdf --no-intro"
(cd "$OUT/layer-group" && stb-symmetry --file graphene_slab.fdf --format fdf \
    --no-intro > console.log 2>&1)
awk '/\[1\] DIMENSIONALITY/{flag=1} /\[4\] LATTICE/{flag=0} flag' "$OUT/layer-group/console.log"
pause

echo "=================================================================="
echo " output/point-group/  --  an isolated molecule needs POINT GROUP"
echo "=================================================================="
cat <<'EOF'
molecule.fdf is a water molecule in a vacuum box (10 Ang padding on all 3
axes -- a genuinely isolated, 0D system). Neither the 3D space group nor
spglib's layer-group detection apply here; [3] POINT GROUP uses pymatgen's
separate, non-periodic PointGroupAnalyzer instead, correctly finding C2v
and identifying the 2 H atoms as symmetry-equivalent:
EOF
mkdir -p "$OUT/point-group"
cp molecule.fdf "$OUT/point-group/"
echo
echo "\$ stb-symmetry --file molecule.fdf --format fdf --no-intro"
(cd "$OUT/point-group" && stb-symmetry --file molecule.fdf --format fdf \
    --no-intro > console.log 2>&1)
awk '/\[3\] POINT GROUP/{flag=1} /\[4\] LATTICE/{flag=0} flag' "$OUT/point-group/console.log"
pause

echo "=================================================================="
echo " output/limitation-wire/  --  a real, documented limitation: no rod groups"
echo "=================================================================="
cat <<'EOF'
wire.fdf has vacuum on 2 axes (a 1D-periodic chain). A slab (1 vacuum
axis) gets LAYER GROUP; an isolated molecule (3 vacuum axes) gets POINT
GROUP -- but spglib has no "rod group" detection at all (the periodic
analogue for a 1D system), so a wire gets neither section. This is
called out explicitly rather than silently falling back to the
physically-meaningless 3D space group:
EOF
mkdir -p "$OUT/limitation-wire"
cp wire.fdf "$OUT/limitation-wire/"
echo
echo "\$ stb-symmetry --file wire.fdf --format fdf --no-intro"
(cd "$OUT/limitation-wire" && stb-symmetry --file wire.fdf --format fdf \
    --no-intro > console.log 2>&1)
awk '/\[1\] DIMENSIONALITY/{flag=1} /\[2\] SPACE/{flag=0} flag' "$OUT/limitation-wire/console.log"
pause

echo "=================================================================="
echo " output/full-report/  --  --save-report (off by default)"
echo "=================================================================="
cat <<'EOF'
Every run always prints the numbered [0]...[11] report to the console.
--save-report additionally persists it to stb_symmetry_report.txt -- off
by default, so a plain run only ever writes references.bib, no text
report file at all:
EOF
mkdir -p "$OUT/full-report"
cp nacl.fdf "$OUT/full-report/"
echo
echo "\$ stb-symmetry --file nacl.fdf --format fdf --no-intro   # default: no report file"
(cd "$OUT/full-report" && stb-symmetry --file nacl.fdf --format fdf \
    --no-intro > console_default.log 2>&1)
ls "$OUT/full-report/" | grep -v "\.fdf$\|console"
echo "(no symmetry.dat, no stb_symmetry_report.txt -- only references.bib)"
echo
echo "\$ stb-symmetry --file nacl.fdf --format fdf --save-report --no-intro"
(cd "$OUT/full-report" && stb-symmetry --file nacl.fdf --format fdf \
    --save-report --no-intro > console_saved.log 2>&1)
echo "Report sections written to stb_symmetry_report.txt:"
grep -E "^\[[0-9]+\]" "$OUT/full-report/stb_symmetry_report.txt"
echo
echo "references.bib -- SIESTA (every stb-symmetry run analyzes a SIESTA structure file):"
grep "^@" "$OUT/full-report/references.bib"
pause

echo "=================================================================="
echo " Two ways to run it"
echo "=================================================================="
cat <<'EOF'
A -- direct CLI:
  stb-symmetry --file nacl.fdf --format fdf

B -- interactive stb-suite menu:
  $ stb-suite
  Select an option (0-6, or a tool code like 4.1.2): 3.5

Both paths call the exact same underlying tool -- proven directly below.
EOF
TMP="$(mktemp -d)"
cp nacl.fdf "$TMP/"
echo
echo "\$ printf '3.5\\nnacl.fdf\\n1\\n\\n\\n\\n\\n\\n\\n\\n\\n' | stb-suite     # format=fdf, all other prompts default/skip"
(cd "$TMP" && printf '3.5\nnacl.fdf\n1\n\n\n\n\n\n\n\n\n' | stb-suite > session.log 2>&1) || true
CLI_LINE=$(grep "Fm-3m" "$OUT/bulk-nacl/console.log" | head -1)
MENU_LINE=$(grep "Fm-3m" "$TMP/session.log" | head -1)
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
Not exercised by this script (needs a display): --view opens an
interactive 3D view (via ASE) of the analyzed structure -- or, with
--write-refined/--compare-to also given, both structures side by side,
paging through frames in ase-gui:

  stb-symmetry --file nacl.fdf --format fdf --write-refined refined.fdf --view
EOF
pause

echo "=================================================================="
echo " Done"
echo "=================================================================="
cat <<EOF
Six self-contained folders were generated under output/:
  bulk-nacl/       scan-symprec/    layer-group/
  point-group/     limitation-wire/ full-report/

Each has references.bib; full-report/ additionally has
stb_symmetry_report.txt (only from its --save-report run).

Recap of what this walkthrough covered:
  - space group, Wyckoff orbits, symmetry operations, and per-atom
    distortion on a textbook high-symmetry structure (NaCl)
  - --scan-symprec revealing symmetry hidden by small numerical noise
  - why a 2D slab needs LAYER GROUP and an isolated molecule needs POINT
    GROUP instead of trusting the plain 3D SPACE GROUP (see the README
    for the theory behind each)
  - the documented "no rod group" limitation for a 1D wire
  - the numbered [0]...[11] report, --save-report, references.bib
  - CLI and the interactive stb-suite menu building the same command

As a next step, try on your own with a real SIESTA structure:
  stb-symmetry --file my_calc.fdf --format fdf --scan-symprec --save-report
  stb-symmetry --file my_calc.fdf --format fdf --compare-to my_relaxed.STRUCT_OUT --compare-format struct_out
EOF
