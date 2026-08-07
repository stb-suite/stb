#!/bin/bash
# Guided example: Workflow 4.6 -- Structure Solution (XRD) (both stages:
# stb-xrdsearch, code 4.6.1, and stb-xrdrank, code 4.6.2, in the stb-suite
# menu).
#
# Not an automated test (see test/4-workflow/6-xrd/{prep,analysis}/test.sh
# for that) -- a commented walk-through: it runs real commands, one case at
# a time, into its own output/<case>/ folder, and shows you the piece of
# output that proves what just happened. Pauses between sections so you can
# read before moving on. Safe to re-run any time -- it always starts by
# wiping its own output/.
#
# Unlike 4.1-4.4's own scripts, there is no structure.fdf/calc.fdf fixture
# here -- this workflow starts from a bare COMPOSITION (README section
# intro), so Stage 1 itself builds everything from scratch, every time.
#
# experimental.dat / experimental_peaks.dat are REAL Si powder-XRD data
# (README section 8), not fabricated -- derived from the same real
# production file (/home/carlos/test/Si/xrd.int) README section 5 documents.
#
# Both stages run for real below. Stage 1's --mace-relax case takes
# roughly a minute (3 small candidates, MACE-MP-0 'small'); each Stage 2
# ranking call takes roughly a minute too (pyxtal.XRD.Similarity has no
# faster mode, ~20s/candidate x 3 candidates) -- expect this script to run
# for several minutes total, most of it real physics, not overhead.

set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

# Stage 2's --view calls plt.show() -- MPLBACKEND=Agg makes that a no-op
# instead of blocking on a GUI window, same convention example_4.4.sh uses.
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
echo " Welcome: casting candidates, then ranking them against real XRD data"
echo "=================================================================="
cat <<'EOF'
Stage 1 (stb-xrdsearch) casts random atomic arrangements for a composition
across one or more candidate space groups, and gets each one ready for
SIESTA (calc.fdf + pseudopotentials, both auto/default unless you say
otherwise). Stage 2 (stb-xrdrank) simulates each candidate's own powder
XRD pattern and ranks them by similarity to a real experimental pattern.

See this folder's README.md for the full theory (Bragg's law, the
structure factor, the similarity metric -- section 1), which external
libraries do the actual work (section 2), and a real production example
where this exact workflow returned nan/gave a wrong ranking before a real
bug fix, and correct numbers after (section 5).
EOF
pause


echo "=================================================================="
echo " Case 1: Stage 1, default run -- auto-generated calc.fdf"
echo "=================================================================="
cat <<'EOF'
No structure file needed -- just a composition (--species/--num-ions) and
a space group to try. By default calc.fdf is auto-generated: the exact
same standard relaxation template stb-inputfile -t relax produces.
EOF
mkdir -p "$OUT/case1-default"
echo "\$ stb-xrdsearch --species Si --num-ions 8 --groups 227 --no-intro"
(cd "$OUT/case1-default" && stb-xrdsearch --species Si --num-ions 8 --groups 227 --no-intro \
    | sed -n '/\[1\] CALC.FDF SETUP/,/\[2\] CANDIDATE CASTING/p')
echo
echo "On disk:"
(cd "$OUT/case1-default" && find xrd_search -maxdepth 2 | sort)
pause


echo "=================================================================="
echo " Case 2: Stage 1, --calc a custom template instead"
echo "=================================================================="
cat <<'EOF'
--calc switches to the second mode: your own template, copied verbatim
into every candidate folder -- always still named calc.fdf in either mode,
so a batch stays runnable the identical way regardless of which mode
produced it (README section 3.2).
EOF
mkdir -p "$OUT/case2-custom-calc"
cat > "$OUT/case2-custom-calc/my_template.fdf" <<'CALCEOF'
SystemLabel      siesta
%include structure.fdf
PAO.BasisType    split
PAO.BasisSize    DZP
PAO.EnergyShift  0.01 Ry
XC.functional    GGA
XC.authors       PBE
MeshCutoff       300 Ry
MD.TypeOfRun     CG
MD.VariableCell  true
MD.MaxForceTol   0.02 eV/Ang
MD.Steps         100
CALCEOF
echo "\$ stb-xrdsearch --species Si --num-ions 8 --groups 227 --calc my_template.fdf --no-intro"
(cd "$OUT/case2-custom-calc" && stb-xrdsearch --species Si --num-ions 8 --groups 227 \
    --calc my_template.fdf --no-intro | sed -n '/\[1\] CALC.FDF SETUP/,/\[2\] CANDIDATE CASTING/p')
pause


echo "=================================================================="
echo " Case 3: Stage 1, --mace-relax across 3 space groups (1 correct, 2 wrong)"
echo "=================================================================="
cat <<'EOF'
227 (Fd-3m) is the true diamond-cubic structure of silicon. 225 and 229
are physically WRONG space groups for elemental Si (single-site FCC/BCC
packing, not the tetrahedrally-bonded diamond motif) -- included on
purpose, to give Stage 2 something decisive to reject below.

--mace-relax pre-relaxes the FINAL candidates (positions + cell) with
MACE-MP-0 before calc.fdf is generated. Real, physically meaningful
result: 227 relaxes to -42.97 eV; 225/229 both collapse under their own
imposed symmetry into the SAME higher-energy -41.44 eV minimum -- already,
before any XRD comparison, MACE's own energy favors the true structure by
~1.5 eV (README section 3.4).

This is the candidate set every Stage 2 case below reuses.
EOF
mkdir -p "$OUT/xrd_search_shared"
echo "\$ stb-xrdsearch --species Si --num-ions 8 --groups 227,225,229 --mace-relax -p dojo --no-intro"
(cd "$OUT/xrd_search_shared" && stb-xrdsearch --species Si --num-ions 8 --groups 227,225,229 \
    --mace-relax -p dojo --no-intro 2>/dev/null | sed -n '/\[3\] MACE PRE-RELAXATION/,/\[4\] PSEUDOPOTENTIALS/p')
cp "$DIR/experimental.dat" "$DIR/experimental_peaks.dat" "$OUT/xrd_search_shared/"
pause


echo "=================================================================="
echo " Case 4: Stage 2, default -- a clean, decisive structure solution"
echo "=================================================================="
cat <<'EOF'
--input-dir defaults to 'xrd_search' (Stage 1's own default --output-dir)
-- no flag needed here. experimental.dat is the real Si pattern, the
canonical 2-column model (README section 1.5): '2theta intensity', nothing
else. --save-gnuplot + --view also requested here, together, to avoid a
second ~1-minute ranking pass just for plots (--view is a safe no-op under
MPLBACKEND=Agg).
EOF
echo "\$ stb-xrdrank --experimental experimental.dat --save-gnuplot --view --no-intro"
(cd "$OUT/xrd_search_shared" && stb-xrdrank --experimental experimental.dat \
    --save-gnuplot --view --no-intro 2>/dev/null | sed -n '/\[2\] EXPERIMENTAL PATTERN/,/\[5\] SUMMARY/p')
echo
echo "plot/ contents:"
(cd "$OUT/xrd_search_shared" && ls plot | sort)
pause


echo "=================================================================="
echo " Case 5: Stage 2, peak-list auto-broadening -- why it matters, live"
echo "=================================================================="
cat <<'EOF'
experimental_peaks.dat is the SAME real pattern, but as a sparse 9-point
peak list (the shape a literature table has) instead of a dense scan.
Auto-detected and Gaussian-broadened by default before comparing.
EOF
echo "\$ stb-xrdrank --experimental experimental_peaks.dat --no-intro"
(cd "$OUT/xrd_search_shared" && stb-xrdrank --experimental experimental_peaks.dat --no-intro 2>/dev/null \
    | sed -n '/\[2\] EXPERIMENTAL PATTERN/,/\[4\] OUTPUT/p')
cat <<'EOF'

Now the SAME file with --raw-experimental (skip broadening):
EOF
echo "\$ stb-xrdrank --experimental experimental_peaks.dat --raw-experimental --no-intro"
(cd "$OUT/xrd_search_shared" && stb-xrdrank --experimental experimental_peaks.dat --raw-experimental --no-intro 2>/dev/null \
    | sed -n '/\[2\] EXPERIMENTAL PATTERN/,/\[4\] OUTPUT/p')
cat <<'EOF'

Without broadening, the correct structure (227) LOSES to the two wrong
ones -- a complete inversion of Case 4's correct ranking (README section
4.4). This is a real demonstration, not a synthetic worst-case: comparing
a zero-width stick list directly against a broadened simulated profile
draws a spurious interpolated curve through the empty gaps between peaks.
EOF
pause


echo "=================================================================="
echo " Case 6: the real 3-column bug, reproduced and shown fixed"
echo "=================================================================="
cat <<'EOF'
Real diffractometer '.int' exports often have a 3rd column (e.g. an
unpopulated uncertainty/esd field, exactly zero on every line) -- this
USED TO silently get read as "intensity" (the old "always take the last
column" rule), giving similarity = nan with no visible cause beyond a bare
RuntimeWarning (README section 4.5/5.1, the exact real bug this session
found and fixed against /home/carlos/test/Si's own production data).
Reproduced here with the same real pattern plus a synthetic bogus 3rd
column mimicking that exact real-world mistake:
EOF
awk '$0 !~ /^#/ {print $1, $2, "0.00000"}' "$DIR/experimental.dat" > "$OUT/xrd_search_shared/experimental_3col.dat"
head -3 "$OUT/xrd_search_shared/experimental_3col.dat"
echo "\$ stb-xrdrank --experimental experimental_3col.dat --no-intro"
(cd "$OUT/xrd_search_shared" && stb-xrdrank --experimental experimental_3col.dat --no-intro 2>/dev/null \
    | sed -n '/\[3\] RANKING/,/\[4\] OUTPUT/p')
echo
echo "Identical to Case 4's clean 2-column result -- the fix correctly ignores the bogus constant column."
pause


echo "=================================================================="
echo " Case 7: a genuinely unusable file -- now a clear error, not a silent nan"
echo "=================================================================="
cat <<'EOF'
If EVERY column after 2-theta is constant (nothing usable at all as an
intensity signal), the fix raises a clear error instead of silently
proceeding to another nan (README section 4.5).
EOF
printf "1 0\n2 0\n3 0\n4 0\n5 0\n" > "$OUT/xrd_search_shared/broken.dat"
echo "\$ stb-xrdrank --experimental broken.dat --no-intro"
set +e
(cd "$OUT/xrd_search_shared" && stb-xrdrank --experimental broken.dat --no-intro 2>/dev/null \
    | grep "FAIL")
set -e
pause


echo "=================================================================="
echo " Case 8: [6] LIBRARY WARNINGS -- everything third-party, one place"
echo "=================================================================="
cat <<'EOF'
Every third-party call (pyxtal's pattern/similarity computation,
pymatgen/spglib's space-group lookup) is wrapped so its own stdout prints
and warnings.warn() calls land in ONE final section instead of leaking
into the terminal interleaved with the report (README section 4.7).
EOF
echo "\$ stb-xrdrank --experimental experimental.dat --no-intro"
(cd "$OUT/xrd_search_shared" && stb-xrdrank --experimental experimental.dat --no-intro 2>/dev/null \
    | sed -n '/\[6\] LIBRARY WARNINGS/,$p')
echo
echo "That line never appears anywhere else in the output -- not before the report banner, not interleaved with [3]'s own live progress lines."
pause


echo "=================================================================="
echo " Case 9: the interactive path (stb-suite -> 4.6.1 then 4.6.2)"
echo "=================================================================="
cat <<'EOF'
Both stages are reachable via stb-suite's interactive menu (dotted codes
4.6.1/4.6.2) with the exact same questions as the CLI flags above --
--input-dir defaults to 'xrd_search' there too (blank = accept).
EOF
mkdir -p "$OUT/case9-interactive"
echo "\$ stb-suite  (4.6.1: Si8, group 227, no MACE relax, no pseudo-dir)"
(cd "$OUT/case9-interactive" && printf '4.6.1\nSi 8\n\n227\n\nn\n\nn\nn\n\n\nn\nn\n\n0\n' | stb-suite 2>&1 \
    | sed -n '/\[2\] CANDIDATE CASTING/,/\[3\] MACE/p')
cp "$DIR/experimental.dat" "$OUT/case9-interactive/"
echo "\$ stb-suite  (4.6.2: default input folder, experimental.dat)"
(cd "$OUT/case9-interactive" && printf '4.6.2\n\nexperimental.dat\n\nn\n\n\nn\nn\nn\n\n0\n' | stb-suite 2>&1 \
    | sed -n '/\[3\] RANKING/,/Score gap/p')
pause


echo "=================================================================="
echo " Workflow 4.6 complete"
echo "=================================================================="
cat <<'EOF'
Stage 1 (stb-xrdsearch, 4.6.1): casts candidates from a bare composition +
space group(s), auto-generates calc.fdf (or copies your own template),
copies pseudopotentials, and optionally pre-relaxes with MACE -- every
candidate folder is ready for real SIESTA as-is.

Stage 2 (stb-xrdrank, 4.6.2): simulates each candidate's XRD pattern and
ranks by similarity to a real experimental pattern -- the canonical input
is 2 columns (2theta, intensity); a real 3-column-format bug that used to
silently give nan is now fixed and demonstrated live above (Cases 6-7);
[6] LIBRARY WARNINGS collects every third-party message in one place.

See the README's section 5 for the real /home/carlos/test/Si production
data this session actually debugged (the same nan -> real-numbers story,
on real SIESTA-relaxed candidates, plus a genuine relaxation-quality
finding the ranking alone can't resolve), and section 7 for a step-by-step
guide to running this on your own composition.
EOF
