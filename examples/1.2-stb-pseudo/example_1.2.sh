#!/bin/bash
# Guided example: stb-pseudo (code 1.2 in the stb-suite menu)
#
# Not an automated test (see test/1-inputs/2-pseudo/test.sh for that) -- a
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
echo " Why SIESTA needs a pseudopotential per species"
echo "=================================================================="
cat <<'EOF'
A pseudopotential replaces the true nucleus-plus-core-electrons Coulomb
potential with a smoother EFFECTIVE potential, tuned so that outside a
chosen core radius r_c the resulting pseudo-wavefunctions reproduce the
real (all-electron) valence wavefunctions almost exactly -- the core
electrons themselves are never resolved explicitly, which is what makes
the calculation cheap.

"Norm-conserving" pseudopotentials (both bundled banks below are this
kind) add one more constraint: pseudo- and all-electron wavefunctions
enclose the SAME CHARGE inside r_c. That single constraint is what makes
a pseudopotential TRANSFERABLE -- correct scattering behaviour across
different chemical environments, not just the one it was fitted on.

SIESTA resolves a species' pseudopotential purely by matching the FILE
NAME to the species LABEL declared in %block ChemicalSpeciesLabel, not
the atomic number -- so the file just needs to be named "<Label>.psf" or
"<Label>.psml" and sit in the run folder.
EOF
pause

echo "=================================================================="
echo " Browsing a bank without any structure at all"
echo "=================================================================="
echo "\$ stb-pseudo --list-elements dojo --no-intro"
stb-pseudo --list-elements dojo --no-intro
pause

echo "=================================================================="
echo " output/sic/  --  resolve + copy from a real structure"
echo "=================================================================="
echo "sic.fdf is a 2-atom SiC cell -- both Si and C are covered by 'dojo'."
mkdir -p "$OUT/sic"
echo
echo "\$ stb-pseudo -f sic.fdf -p dojo -o output/sic --save-report --no-intro"
stb-pseudo -f sic.fdf -p dojo -o "$OUT/sic" --save-report --no-intro > "$OUT/sic_console.log"
grep -E "Species found|Resolved" "$OUT/sic_console.log"
echo
echo "Copied into output/sic/:"
ls "$OUT/sic"
pause

echo "=================================================================="
echo " Bundled banks: what's actually inside dojo vs virtual_vault"
echo "=================================================================="
DOJO_N=$(stb-pseudo --list-elements dojo --no-intro | grep "element(s) available" | sed 's/ element.*//')
VV_N=$(stb-pseudo --list-elements virtual_vault --no-intro | grep "element(s) available" | sed 's/ element.*//')
cat <<EOF
  dojo           -- PseudoDojo v0.5, PBE, .psml format, $DOJO_N elements
  virtual_vault  -- SIESTA Virtual Vault, PBE, .psf format, $VV_N elements

Neither bank covers everything -- an overlapping but not identical set.
That's exactly why --fallback-dir exists below, instead of forcing you
to commit to one bank and live with its gaps.
EOF
pause

echo "=================================================================="
echo " output/fallback/  --  filling a gap with --fallback-dir"
echo "=================================================================="
echo "Astatine (At) is absent from 'dojo' but present in 'virtual_vault'."
echo
echo "Without a fallback, At is reported (and left) missing:"
echo "\$ stb-pseudo --species Si At -p dojo --dry-run --no-intro"
stb-pseudo --species Si At -p dojo --dry-run --no-intro > "$OUT/no_fallback.log" || true
grep -E "Si  FOUND|At  MISSING|Resolved" "$OUT/no_fallback.log"
echo
echo "With a fallback, the same At is resolved from the second source:"
mkdir -p "$OUT/fallback"
echo "\$ stb-pseudo --species Si At -p dojo --fallback-dir virtual_vault -o output/fallback --no-intro"
stb-pseudo --species Si At -p dojo --fallback-dir virtual_vault -o "$OUT/fallback" --no-intro > "$OUT/fallback_console.log"
grep -E "FOUND|Resolved" "$OUT/fallback_console.log"
echo
echo "Copied into output/fallback/:"
ls "$OUT/fallback"
pause

echo "=================================================================="
echo " output/dry_run/  --  report only, nothing copied"
echo "=================================================================="
mkdir -p "$OUT/dry_run"
echo "\$ stb-pseudo -f sic.fdf -p dojo -o output/dry_run --dry-run --no-intro"
stb-pseudo -f sic.fdf -p dojo -o "$OUT/dry_run" --dry-run --no-intro > "$OUT/dry_run_console.log"
grep "no files copied" "$OUT/dry_run_console.log"
echo
echo "output/dry_run/ is empty:"
ls -A "$OUT/dry_run" | wc -l | xargs echo "  file count:"
pause

echo "=================================================================="
echo " Proof: CLI and the interactive stb-suite menu agree"
echo "=================================================================="
echo "Driving the same sic.fdf/dojo case through the interactive menu"
echo "(non-interactively, via a piped printf) and comparing the"
echo "'Resolved : N/N' line against the direct CLI run above."
TMP="$(mktemp -d)"
cp sic.fdf "$TMP/"
echo
echo "\$ printf '1.2\\n1\\nsic.fdf\\ndojo\\n\\n\\n\\nn\\n\\n0\\n' | stb-suite"
(cd "$TMP" && printf '1.2\n1\nsic.fdf\ndojo\n\n\n\nn\n\n0\n' | stb-suite > session.log 2>&1)
CLI_RESOLVED=$(grep "Resolved :" "$OUT/sic_console.log" | head -1)
MENU_RESOLVED=$(grep "Resolved :" "$TMP/session.log" | head -1)
if [ "$CLI_RESOLVED" = "$MENU_RESOLVED" ]; then
    echo "Confirmed: identical result ($MENU_RESOLVED) from the CLI and the interactive menu."
else
    echo "Unexpected: results differ -- CLI: $CLI_RESOLVED / menu: $MENU_RESOLVED"
fi
rm -rf "$TMP"
pause

echo "=================================================================="
echo " Done"
echo "=================================================================="
cat <<'EOF'
Three self-contained folders were generated under output/:
  sic/   fallback/   dry_run/

output/sic/ and output/fallback/ each hold the actual pseudopotential
files, ready to reference from a calc.fdf; sic/ also has stb_pseudo_
report.txt (via --save-report) and this folder's own references.bib
(SIESTA + PseudoDojo citations).

As a next step, try on your own:
  stb-pseudo --species Fe O -p dojo -o out/
  stb-pseudo --list-elements virtual_vault

In a real workflow, point stb-inputfile's (example 1.1) -p/--pp-path
straight at whichever output/ folder stb-pseudo filled for you.
EOF
