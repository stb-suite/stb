#!/bin/bash
# Guided example: stb-fetch (code 1.6 in the stb-suite menu)
#
# Not an automated test -- a commented walk-through: it runs real commands,
# one group at a time, and shows you the piece of output that proves what
# just happened. It pauses between sections so you can read before moving
# on.
#
# Unlike every other example so far (1.1-1.4), this one needs INTERNET
# ACCESS -- stb-fetch is the suite's first network-dependent tool, querying
# the Crystallography Open Database (COD) and an OPTIMADE provider live. If
# a request fails (offline, or a remote endpoint is down), that section will
# error out -- re-run once connectivity is back.

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
echo " What stb-fetch does"
echo "=================================================================="
cat <<'EOF'
stb-fetch looks up a structure in an online database -- the Crystallography
Open Database (COD, no account needed), Materials Project (needs a free
PMG_MAPI_KEY), or any OPTIMADE-compliant database (AFLOW, JARVIS, 2D
materials databases, ...) -- and writes it straight out as a SIESTA .fdf,
so you never have to hand-build a structure file for a material that
already exists in the literature.

Real database entries are frequently disordered even for common materials
(e.g. magnetite's octahedral Fe site is often given as a Fe2+/Fe3+ mixed
-occupancy split). Same-element-different-oxidation-state disorder is
collapsed automatically to a single element; genuine multi-element disorder
(two different elements sharing a site) is refused with a clear error
instead of silently guessing.
EOF
pause

echo "=================================================================="
echo " Two ways to run it"
echo "=================================================================="
cat <<'EOF'
A -- direct CLI:
  stb-fetch --source cod --cod-id 1010369

B -- interactive stb-suite menu:
  $ stb-suite
  Select an option (0-6, or a tool code like 4.1.2): 1.6

"1.6" walks you through picking a source/provider, an exact id or a formula
search, an optional unit-cell reduction, and the save-report/--view prompts
-- the same choices the CLI flags below control directly.
EOF
pause

echo "=================================================================="
echo " output/cod-by-id/  --  exact id, bulk 3D material"
echo "=================================================================="
cat <<'EOF'
COD entry 1010369 is magnetite (Fe3O4) -- a real crystal with the
oxidation-state disorder mentioned above. Watch the report: the disorder
note in [2], and [3] STRUCTURE VALIDATION reporting a plain 3D space group
(no vacuum axis, so no layer/point group applies here).
EOF
mkdir -p "$OUT/cod-by-id"
echo "\$ stb-fetch --source cod --cod-id 1010369 --save-report --no-intro -o fe3o4.fdf"
(cd "$OUT/cod-by-id" && stb-fetch --source cod --cod-id 1010369 --save-report --no-intro -o fe3o4.fdf > console.log)
grep -E "disorder|Dimensionality|Space group|Layer group|Point group" "$OUT/cod-by-id/console.log"
pause

echo "=================================================================="
echo " output/cod-primitive/  --  same entry, reduced to its primitive cell"
echo "=================================================================="
mkdir -p "$OUT/cod-primitive"
echo "\$ stb-fetch --source cod --cod-id 1010369 --unitcell primitive --save-report --no-intro -o fe3o4_prim.fdf"
(cd "$OUT/cod-primitive" && stb-fetch --source cod --cod-id 1010369 --unitcell primitive --save-report --no-intro -o fe3o4_prim.fdf > console.log)
grep -E "^Atoms|Unit cell mode|Output atoms" "$OUT/cod-primitive/console.log"
echo
echo "56 conventional-cell atoms reduced down to the primitive cell above."
pause

echo "=================================================================="
echo " output/optimade-2d/  --  a genuinely 2D material: the layer group"
echo "=================================================================="
cat <<'EOF'
An ordinary 3D "space group" isn't the physically correct symmetry
classification for a slab: it treats the vacuum gap as just an unusually
tall periodic cell. For a structure with exactly one vacuum-padded axis,
stb-fetch now also reports the LAYER GROUP (same detection stb-symmetry,
code 3.5, uses) -- only the label here, not the full list of operations;
use stb-symmetry directly if you need those.

twodmatpedia (an OPTIMADE provider, no API key needed) entry 2dm-3150 is a
monolayer MoS2 -- exactly one vacuum-padded axis.
EOF
mkdir -p "$OUT/optimade-2d"
echo "\$ stb-fetch --source optimade --provider twodmatpedia --optimade-id 2dm-3150 --save-report --no-intro -o mos2.fdf"
(cd "$OUT/optimade-2d" && stb-fetch --source optimade --provider twodmatpedia --optimade-id 2dm-3150 --save-report --no-intro -o mos2.fdf > console.log)
grep -E "Dimensionality|Space group|Layer group" "$OUT/optimade-2d/console.log"
pause

echo "=================================================================="
echo " Multiple candidates: --formula without an exact id"
echo "=================================================================="
cat <<'EOF'
Searching by formula can match more than one entry. Non-interactively (as
in this script), stb-fetch never guesses which one you meant -- it prints
the candidate table and asks you to rerun with the exact id.
EOF
mkdir -p "$OUT/formula-search"
echo "\$ stb-fetch --source optimade --provider twodmatpedia --formula MoS2 --limit 5 --no-intro"
set +e
(cd "$OUT/formula-search" && stb-fetch --source optimade --provider twodmatpedia --formula MoS2 --limit 5 --no-intro > console.log 2>&1)
SEARCH_EXIT=$?
set -e
grep -E "Found|\[[0-9]+\]|ERROR" "$OUT/formula-search/console.log"
if [ "$SEARCH_EXIT" -eq 1 ]; then
    echo "Confirmed: exit code 1 -- an informative error, not a silent guess."
fi
pause

echo "=================================================================="
echo " references.bib merges across runs in the same folder"
echo "=================================================================="
cat <<'EOF'
references.bib is always written -- there's no flag for it. If you fetch
from two different sources into the SAME output folder, the citations
merge (by BibTeX key) instead of the second run erasing the first's.
EOF
mkdir -p "$OUT/merged-citations"
echo "\$ stb-fetch --source cod --cod-id 1010369 -o fe3o4.fdf --no-intro"
(cd "$OUT/merged-citations" && stb-fetch --source cod --cod-id 1010369 -o fe3o4.fdf --no-intro > /dev/null)
echo "\$ stb-fetch --source optimade --provider twodmatpedia --optimade-id 2dm-3150 -o mos2.fdf --no-intro"
(cd "$OUT/merged-citations" && stb-fetch --source optimade --provider twodmatpedia --optimade-id 2dm-3150 -o mos2.fdf --no-intro > /dev/null)
echo
echo "Citation keys present after both runs:"
grep "^@" "$OUT/merged-citations/references.bib"
if grep -q "Grazulis2009" "$OUT/merged-citations/references.bib" && grep -q "Andersen2021" "$OUT/merged-citations/references.bib"; then
    echo "Confirmed: both the COD and the OPTIMADE citation are present -- merged, not overwritten."
fi
pause

echo "=================================================================="
echo " Proof: CLI and the interactive stb-suite menu agree"
echo "=================================================================="
echo "Driving the same COD-by-id (1010369) case through the interactive"
echo "menu (main prompt -> '1.6' -> defaults -> save-report) and comparing"
echo "its report against output/cod-by-id/'s."
TMP="$(mktemp -d)"
echo
echo "\$ printf '1.6\\n\\n\\n1010369\\n\\nfetched_menu.fdf\\ny\\nn\\n\\n0\\n' | stb-suite"
(cd "$TMP" && printf '1.6\n\n\n1010369\n\nfetched_menu.fdf\ny\nn\n\n0\n' | stb-suite > session.log 2>&1)
CLI_LINES=$(grep -E "^Formula|^Space group" "$OUT/cod-by-id/console.log")
MENU_LINES=$(grep -E "^Formula|^Space group" "$TMP/stb_fetch_report.txt")
if [ "$CLI_LINES" = "$MENU_LINES" ]; then
    echo "Confirmed: identical Formula/Space group lines from the CLI and the interactive menu."
else
    echo "Note: menu output differs from the direct CLI run -- see $TMP/session.log."
    echo "$MENU_LINES"
fi
rm -rf "$TMP"
pause

echo "=================================================================="
echo " Done"
echo "=================================================================="
cat <<EOF
Folders generated under output/:
  cod-by-id/   cod-primitive/   optimade-2d/   formula-search/   merged-citations/

Each fetch also ran the same structure-validation checks stb-inputfile
(example 1.1) runs on a hand-built structure: atom-proximity, left-handed
-cell, and (bulk-only) density sanity checks, plus a dimension-aware
space/layer/point-group label.

As a next step, try on your own:
  stb-fetch --list-providers                         # see every OPTIMADE alias
  stb-fetch --source cod --formula Si --limit 5       # a formula search
  stb-fetch --source cod --cod-id 1010369 --view      # open it in ASE's viewer

For the FULL symmetry analysis (operations, Wyckoff sites) of a fetched
structure, run stb-symmetry (code 3.5) on the .fdf this tool wrote.
EOF
