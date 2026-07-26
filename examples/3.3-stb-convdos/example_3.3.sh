#!/bin/bash
# Guided example: stb-convdos (code 3.3 in the stb-suite menu)
#
# Not an automated test (see test/3-analysis/3-dos-convolution/test.sh for
# that) -- a commented walk-through: it runs real commands, one group at a
# time, into its own output/<case>/ folder, and shows you the piece of
# output that proves what just happened. Pauses between sections so you
# can read before moving on. Safe to re-run any time -- it always starts
# by wiping its own output/.
#
# dos_total.dat is a small, hand-built synthetic "stick spectrum" (a few
# sharp spikes, not a real SIESTA DOS) -- deliberately built so the
# broadening effect and the conservation check are easy to see and verify
# by hand. multi.PDOS.xml (the same fixture examples/3.2-stb-dos/ uses)
# is processed live via a real stb-dos run to demonstrate --dir on a
# genuine stb-dos output tree.

set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

# A run without --no-plot still calls plt.show() -- MPLBACKEND=Agg makes
# that a no-op instead of blocking on a GUI window, same convention
# test.sh and the other examples in this suite already use.
export MPLBACKEND=Agg

OUT="$DIR/output"
rm -rf "$OUT"
mkdir -p "$OUT"

pause() {
    echo
    read -p "  [Press Enter to continue] " -r
    echo
}

echo "=================================================================="
echo " What Gaussian broadening of a DOS is, and why it's needed"
echo "=================================================================="
cat <<'EOF'
A DFT calculation gives the DOS at a finite, discrete set of k-points --
especially with a coarse k-mesh, the "raw" DOS looks like a comb of
sharp spikes rather than the smooth curve a bulk crystal's DOS actually
is. stb-convdos convolves each DOS column with a normalized Gaussian
kernel to simulate finite instrumental/thermal broadening (or just to
turn a spiky, under-sampled DOS into something presentable).

stb-convdos --file <dos.dat> --sigma <meV> --out <filtered.dat>
  broadens one file (typically stb-dos's own output).
stb-convdos --dir <folder> --sigma <meV>
  recursively broadens EVERY .dat file under a folder (e.g. a whole
  stb-dos output tree) with the same width, mirroring the folder
  structure into --output-dir (default: <folder>_filtered).

Every run prints a numbered [0]...[6] report. --save-report additionally
persists it to stb_convdos_report.txt -- off by default.
EOF
pause

echo "=================================================================="
echo " output/basic-broadening/  --  sharp spikes -> smooth peaks"
echo "=================================================================="
cat <<'EOF'
dos_total.dat: 41 points (-2 to +2 eV, 0.1 eV spacing), 2 columns (s, p)
that are exactly zero except for 2 sharp spikes each -- a toy "stick
spectrum". Watch a spike of height 10.0 at E=-1.0 eV spread into a smooth
bump after broadening, and [3] CONSERVATION CHECK confirm the total area
under each column barely changes (a normalized kernel redistributes
weight in energy, it doesn't create or destroy states):
EOF
mkdir -p "$OUT/basic-broadening"
cp dos_total.dat "$OUT/basic-broadening/"
echo
echo "\$ stb-convdos --file dos_total.dat --sigma 100 --out dos_filtered.dat --no-plot --no-intro"
(cd "$OUT/basic-broadening" && stb-convdos --file dos_total.dat --sigma 100 \
    --out dos_filtered.dat --no-plot --no-intro > console.log 2>&1)
awk '/\[2\] BROADENING/{flag=1} /\[4\]/{flag=0} flag' "$OUT/basic-broadening/console.log"
echo
echo "The spike at E=-1.0 eV (s column, height 10.0) before vs. after broadening:"
grep "^-1.000000\|^-1.100000\|^-0.900000" "$OUT/basic-broadening/dos_total.dat"
echo "  -> becomes, in the filtered file:"
grep "^-1.000000\|^-1.100000\|^-0.900000" "$OUT/basic-broadening/dos_filtered.dat"
pause

echo "=================================================================="
echo " output/sigma-vs-fwhm/  --  --sigma and the equivalent --fwhm give the same result"
echo "=================================================================="
cat <<'EOF'
FWHM = 2.3548 * sigma (the standard Gaussian relation) -- --fwhm is
converted to the equivalent sigma internally before anything else
happens. --sigma 50 and --fwhm 117.741 (= 50 * 2.3548) should therefore
report the exact same internal sigma, and produce byte-identical output:
EOF
mkdir -p "$OUT/sigma-vs-fwhm"
cp dos_total.dat "$OUT/sigma-vs-fwhm/"
echo
echo "\$ stb-convdos --file dos_total.dat --sigma 50 --out via_sigma.dat --no-plot --no-intro"
(cd "$OUT/sigma-vs-fwhm" && stb-convdos --file dos_total.dat --sigma 50 \
    --out via_sigma.dat --no-plot --no-intro > console_sigma.log 2>&1)
grep "Sigma " "$OUT/sigma-vs-fwhm/console_sigma.log"
echo
echo "\$ stb-convdos --file dos_total.dat --fwhm 117.741 --out via_fwhm.dat --no-plot --no-intro"
(cd "$OUT/sigma-vs-fwhm" && stb-convdos --file dos_total.dat --fwhm 117.741 \
    --out via_fwhm.dat --no-plot --no-intro > console_fwhm.log 2>&1)
grep "Sigma " "$OUT/sigma-vs-fwhm/console_fwhm.log"
echo
diff -q "$OUT/sigma-vs-fwhm/via_sigma.dat" "$OUT/sigma-vs-fwhm/via_fwhm.dat" > /dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "Confirmed: --sigma 50 and --fwhm 117.741 write byte-identical output."
else
    echo "Unexpected: the two outputs differ."
fi
pause

echo "=================================================================="
echo " output/dir-mode/  --  broadening a whole stb-dos output tree at once"
echo "=================================================================="
cat <<'EOF'
multi.PDOS.xml (3 atoms/2 species) run through stb-dos gives a real
output tree: dos_total.dat + dos_per_atom/*.dat + dos_per_species/*.dat.
--dir broadens every .dat file found under it with the same width,
mirroring the folder structure (same subfolders, same filenames) into
--output-dir. Also watch: a deliberately broken .dat dropped into the
tree is skipped with a warning instead of aborting the whole batch, and
only the FIRST processed file's plot would ever be shown (not exercised
here, --no-plot is used, but see the README for why that's the default
--dir behavior instead of one plot window per file):
EOF
mkdir -p "$OUT/dir-mode/dos_output"
cp multi.PDOS.xml "$OUT/dir-mode/"
(cd "$OUT/dir-mode/dos_output" && stb-dos ../multi.PDOS.xml --shift fermi --no-intro > /dev/null 2>&1)
echo "Generated tree:"
(cd "$OUT/dir-mode" && find dos_output -name "*.dat" | sort)
echo
echo "# a deliberately invalid file dropped into the tree -- only 1 column, no energy grid"
echo "1.0" > "$OUT/dir-mode/dos_output/not_a_dos_file.dat"
echo "2.0" >> "$OUT/dir-mode/dos_output/not_a_dos_file.dat"
echo
echo "\$ stb-convdos --dir dos_output --sigma 50 --no-plot --no-intro"
(cd "$OUT/dir-mode" && stb-convdos --dir dos_output --sigma 50 --no-plot --no-intro > console.log 2>&1)
awk '/\[1\] INPUT DATA/{flag=1} /\[2\]/{flag=0} flag' "$OUT/dir-mode/console.log"
echo
echo "Mirrored output tree (dos_output_filtered/) -- same subfolders, same filenames:"
(cd "$OUT/dir-mode" && find dos_output_filtered -type f | sort)
pause

echo "=================================================================="
echo " output/full-report/  --  --save-report (off by default)"
echo "=================================================================="
cat <<'EOF'
Every run always prints the numbered [0]...[6] report to the console.
--save-report additionally persists it to stb_convdos_report.txt -- off
by default, so a plain run only ever writes the filtered .dat file(s)
and references.bib, no text report file:
EOF
mkdir -p "$OUT/full-report"
cp dos_total.dat "$OUT/full-report/"
echo
echo "\$ stb-convdos --file dos_total.dat --sigma 100 --out dos_filtered.dat --no-plot --no-intro"
(cd "$OUT/full-report" && stb-convdos --file dos_total.dat --sigma 100 \
    --out dos_filtered.dat --no-plot --no-intro > console_default.log 2>&1)
ls "$OUT/full-report/" | grep -v "\.dat$\|console"
echo "(no stb_convdos_report.txt -- only the filtered .dat + references.bib)"
echo
echo "\$ stb-convdos --file dos_total.dat --sigma 100 --out dos_filtered.dat --no-plot --save-report --no-intro"
(cd "$OUT/full-report" && stb-convdos --file dos_total.dat --sigma 100 --out dos_filtered.dat \
    --no-plot --save-report --no-intro > console_saved.log 2>&1)
echo "Report sections written to stb_convdos_report.txt:"
grep -E "^\[[0-9]+\]" "$OUT/full-report/stb_convdos_report.txt"
echo
echo "references.bib -- SIESTA (every DOS file here is ultimately SIESTA-derived):"
grep "^@" "$OUT/full-report/references.bib"
pause

echo "=================================================================="
echo " Two ways to run it"
echo "=================================================================="
cat <<'EOF'
A -- direct CLI:
  stb-convdos --file dos_total.dat --sigma 50 --out dos_filtered.dat

B -- interactive stb-suite menu:
  $ stb-suite
  Select an option (0-6, or a tool code like 4.1.2): 3.3

Both paths call the exact same underlying tool -- proven directly below.
The menu adds one extra first choice (single file vs. a whole folder)
before the usual --sigma/--fwhm/--size/plot/save-report prompts.
EOF
TMP="$(mktemp -d)"
cp dos_total.dat "$TMP/"
echo
echo "\$ printf '3.3\\n1\\ndos_total.dat\\n\\n\\n50\\n\\nn\\nn\\n' | stb-suite     # mode=1 (file), sigma 50 meV"
(cd "$TMP" && printf '3.3\n1\ndos_total.dat\n\n\n50\n\nn\nn\n' | stb-suite > session.log 2>&1) || true
CLI_LINE=$(grep "Sigma " "$OUT/sigma-vs-fwhm/console_sigma.log")
MENU_LINE=$(grep "Sigma " "$TMP/session.log")
if [ "$CLI_LINE" = "$MENU_LINE" ]; then
    echo "Confirmed: identical broadening (sigma/samples/kernel size) from the CLI and the interactive menu."
    echo "  $CLI_LINE"
else
    echo "Unexpected: results differ -- see $TMP/session.log."
fi
rm -rf "$TMP"
pause

echo "=================================================================="
echo " Done"
echo "=================================================================="
cat <<EOF
Four self-contained folders were generated under output/:
  basic-broadening/   sigma-vs-fwhm/
  dir-mode/            full-report/

Recap of what this walkthrough covered:
  - why a DOS needs Gaussian broadening at all (finite k-sampling ->
    spiky, non-physical "raw" DOS), shown on a real spike-to-peak example
  - sigma vs. FWHM (FWHM = 2.3548 * sigma), proven to give identical
    output either way
  - why --sigma/--fwhm are in meV (physical units), converted internally
    using each file's own median energy spacing
  - the automatic (3-sigma) kernel size, and why it must be odd
  - the conservation check, and what a large before/after drift means
  - --dir: broadening an entire stb-dos output tree in one command,
    mirroring its folder structure, skipping an invalid file without
    aborting the batch
  - the numbered [0]...[6] report, --save-report, references.bib
  - CLI and the interactive stb-suite menu building the same command

Not exercised by this script (needs a display): the interactive
before/after plot -- try it yourself without MPLBACKEND=Agg:
  stb-convdos --file dos_total.dat --sigma 100 --out dos_filtered.dat

As a next step, try on your own with a real stb-dos output:
  stb-convdos --dir my_calc_dos/ --sigma 50 --save-report
EOF
