#!/bin/bash
# Guided example: Workflow 4.5 -- Convergence Tests (both stages:
# stb-convergence, code 4.5.1, and stb-convergenceAnalysis, code 4.5.2, in
# the stb-suite menu).
#
# Not an automated test (see test/4-workflow/5-convergence/{prep,analysis}/
# test.sh for that) -- a commented walk-through: it runs real commands, one
# case at a time, into its own output/<case>/ folder, and shows you the
# piece of output that proves what just happened. Pauses between sections so
# you can read before moving on. Safe to re-run any time -- it always starts
# by wiping its own output/.
#
# structure.fdf/calc.fdf are bulk silicon (diamond cubic, DZP/GGA-PBE) -- the
# standard SIESTA convergence-test material, same fixture
# test/4-workflow/5-convergence/ itself uses.
#
# Stage 1 (stb-convergence) runs for real in every case below (cheap -- it
# only writes .fdf files, no SCF). Stage 2 (stb-convergenceAnalysis) is
# demonstrated against small, hand-built calc.out data (same idea as
# 4.1-strain's/4.2-elastic's/4.3-cohesive's own scripts) since this folder
# doesn't invoke real SIESTA -- the target numbers are realistic and
# physically plausible for bulk silicon, but NOT from a completed ab initio
# run (see the README's own Section 4 for the full, honest caveat -- every
# number the tool PRINTS below is still genuine, verbatim tool output).

set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

# Stage 2's --view calls plt.show() -- MPLBACKEND=Agg makes that a no-op
# instead of blocking on a GUI window, same convention test.sh itself uses.
export MPLBACKEND=Agg

OUT="$DIR/output"
rm -rf "$OUT"
mkdir -p "$OUT"

pause() {
    echo
    read -p "  [Press Enter to continue] " -r
    echo
}

# Writes a minimal, real-format SIESTA calc.out (just what
# core/siesta_log.py's get_free_energy/get_scf_convergence/get_max_force/
# get_outcell look for) into $1/calc.out, plus a structure.fdf copy (needed
# for the atom count). $2 = FreeEng (eV, whole cell), $3 = relaxed cubic cell
# side (Ang), $4 = SCF line (empty string omits it, simulating an
# unconverged SCF), $5 = extra "Max <force>" line (empty omits it).
# Illustrative-but-plausible numbers only -- see the header comment above.
write_run_out() {
    local folder="$1" free_eng="$2" side="$3" scf_line="${4-SCF cycle converged after 12 iterations}" max_force="${5:-}"
    mkdir -p "$folder"
    cp structure.fdf "$folder/structure.fdf"
    {
        [ -n "$scf_line" ] && echo "$scf_line"
        echo "siesta: FreeEng =    ${free_eng}"
        echo
        echo "outcell: Unit cell vectors (Ang):"
        echo "        ${side}    0.000000    0.000000"
        echo "        0.000000    ${side}    0.000000"
        echo "        0.000000    0.000000    ${side}"
        echo
        echo "outcell: Cell vector modules (Ang)   :    ${side}    ${side}    ${side}"
        if [ -n "$max_force" ]; then
            echo
            echo "siesta: Atomic forces (eV/Ang):"
            echo "   Max    ${max_force}"
        fi
    } > "$folder/calc.out"
}


echo "=================================================================="
echo " Welcome: why every SIESTA calculation needs a convergence test first"
echo "=================================================================="
cat <<'EOF'
Mesh.CutOff, PAO.EnergyShift, and the k-grid density all have real, physical
effects on the computed total energy AND the relaxed geometry -- effects
that have nothing to do with the actual physics you're trying to compute.
stb-convergence (Stage 1) sweeps one or more of these, fully relaxing
(positions + cell) at every single value -- never a fixed-geometry shortcut,
since a too-coarse setting can bias exactly where SIESTA relaxes atoms to
(the "eggbox effect", README section 1.3). stb-convergenceAnalysis (Stage 2)
reads the finished runs back and reports where BOTH the energy AND the
relaxed structure itself stop changing -- not just one or the other.

See this folder's README.md for the full theory (section 1) and a real,
code-verified worked example where the two criteria disagree (section 4).
EOF
pause


echo "=================================================================="
echo " Case 1: Stage 1, single-parameter default sweep (Mesh.CutOff)"
echo "=================================================================="
cat <<'EOF'
No range flags -- stb-convergence uses its own suggested default range
(100 to 400 Ry, step 50). Watch [2]'s range table, then [5]'s generated
folder list.
EOF
mkdir -p "$OUT/case1-default"
cp structure.fdf calc.fdf "$OUT/case1-default/"
echo "\$ stb-convergence -s structure.fdf -c calc.fdf -p dojo --parameter meshcutoff --no-intro"
(cd "$OUT/case1-default" && stb-convergence -s structure.fdf -c calc.fdf -p dojo \
    --parameter meshcutoff --no-intro \
    | sed -n '/\[2\] PARAMETER SELECTION/,/\[3\] RELAXATION/p')
echo
echo "On disk:"
(cd "$OUT/case1-default" && find convergence_runs -maxdepth 2 -type d | sort)
pause


echo "=================================================================="
echo " Case 2: Stage 1, --parameter all (3 independent sweeps in one call)"
echo "=================================================================="
cat <<'EOF'
--parameter all expands to meshcutoff + energyshift + kgrid -- each its OWN
independent sweep (README section 1.6), one subfolder tree per parameter.
kgrid sweeps a DENSITY (1/Ang), not a raw grid -- watch the resolved
kgrid.MonkhorstPack block for the densest point (0.05 1/Ang).
EOF
mkdir -p "$OUT/case2-all"
cp structure.fdf calc.fdf "$OUT/case2-all/"
echo "\$ stb-convergence -s structure.fdf -c calc.fdf -p dojo --parameter all --no-intro"
(cd "$OUT/case2-all" && stb-convergence -s structure.fdf -c calc.fdf -p dojo \
    --parameter all --no-intro \
    | sed -n '/\[2\] PARAMETER SELECTION/,/\[3\] RELAXATION/p')
echo
echo "3 independent subtrees on disk:"
(cd "$OUT/case2-all" && find convergence_runs -maxdepth 1 -type d | sort)
echo
echo "kgrid density 0.05 1/Ang resolved to (config_extra.fdf):"
grep "kgrid.MonkhorstPack" "$OUT/case2-all/convergence_runs/kgrid/convergence_kgrid_0.0500/config_extra.fdf"
pause


echo "=================================================================="
echo " Case 3: validation -- what stb-convergence correctly refuses"
echo "=================================================================="
cat <<'EOF'
Three real error cases: an unrecognized --parameter value, a range flag for
a parameter that isn't selected, and a non-positive range step. All exit
with a non-zero status rather than silently doing something unintended.
EOF
mkdir -p "$OUT/case3-validation"
cp structure.fdf calc.fdf "$OUT/case3-validation/"
cd "$OUT/case3-validation"
echo "\$ stb-convergence -s structure.fdf -c calc.fdf --parameter bogus --no-intro"
stb-convergence -s structure.fdf -c calc.fdf --parameter bogus --no-intro 2>&1 | tail -1 || true
echo
echo "\$ stb-convergence -s structure.fdf -c calc.fdf --parameter meshcutoff --energyshift-range 0.01 0.05 0.01 --no-intro"
stb-convergence -s structure.fdf -c calc.fdf --parameter meshcutoff \
    --energyshift-range 0.01 0.05 0.01 --no-intro 2>&1 | tail -1 || true
echo
echo "\$ stb-convergence -s structure.fdf -c calc.fdf --parameter meshcutoff --meshcutoff-range 100 400 0 --no-intro"
stb-convergence -s structure.fdf -c calc.fdf --parameter meshcutoff \
    --meshcutoff-range 100 400 0 --no-intro 2>&1 | tail -1 || true
cd "$DIR"
pause


echo "=================================================================="
echo " Case 4: CLI vs. the interactive stb-suite menu (4.5.1) -- same result"
echo "=================================================================="
cat <<'EOF'
The interactive menu asks the same questions instead of flags, then calls
the exact same stb-convergence underneath. Reproducing Case 1's own
generation (single meshcutoff sweep, dojo pseudopotentials, all defaults)
through stb-suite -> 4.5.1 and diffing against it proves the 2 paths are
equivalent.
EOF
mkdir -p "$OUT/case4-interactive"
cp structure.fdf calc.fdf "$OUT/case4-interactive/"
# Prompts in order: structure.fdf (blank -> default), calc.fdf (blank ->
# default), pseudo source (1 -> dojo), parameter select (1 -> meshcutoff),
# customize range (n -> suggested default), relax steps (blank -> 100),
# output dir (blank -> convergence_runs), save report (n), then the "Press
# Enter to continue" pause, then quit.
(cd "$OUT/case4-interactive" && \
    printf '4.5.1\n\n\n1\n1\nn\n\n\nn\n\n0\n' | stb-suite > menu.log 2>&1)
if diff -rq "$OUT/case1-default/convergence_runs" "$OUT/case4-interactive/convergence_runs" > /dev/null 2>&1; then
    echo "CLI (Case 1) and interactive menu (4.5.1) produced byte-identical folders."
else
    echo "UNEXPECTED: CLI and interactive menu results differ -- see"
    echo "  diff -rq '$OUT/case1-default/convergence_runs' '$OUT/case4-interactive/convergence_runs'"
fi
pause


echo "=================================================================="
echo " Case 5: Stage 2 -- energy convergence lies about geometry convergence"
echo "=================================================================="
cat <<'EOF'
The star of this workflow (README section 4): a Mesh.CutOff sweep where the
energy criterion and the relaxed-structure criterion disagree on where
convergence actually happens. Illustrative-but-plausible bulk-silicon
calc.out data (see the header comment above) -- watch [2.1]'s table, then
[4]'s Energy-converged vs. Structure-converged vs. Recommended columns.
EOF
RUN5="$OUT/case5-worked-example"
mkdir -p "$RUN5"
cp structure.fdf "$RUN5/"
write_run_out "$RUN5/convergence_meshcutoff_150.0000" "-500.900000" "5.35"
write_run_out "$RUN5/convergence_meshcutoff_200.0000" "-500.850000" "5.42"
write_run_out "$RUN5/convergence_meshcutoff_250.0000" "-500.845000" "5.45"
write_run_out "$RUN5/convergence_meshcutoff_300.0000" "-500.844000" "5.451"
echo "\$ stb-convergenceAnalysis --dir . --no-intro"
(cd "$RUN5" && stb-convergenceAnalysis --dir . --no-intro \
    | sed -n '/\[2\.1\]/,/\[3\] OUTPUT/p')
pause


echo "=================================================================="
echo " Case 6: Stage 2 -- the Quality column, live (SCF?, F>tol)"
echo "=================================================================="
cat <<'EOF'
Same sweep as Case 5, plus one deliberately poor 50 Ry point: no confirmed
SCF convergence and a large residual force. Watch its Quality flags, and
the [WARNING] summary line -- advisory only, it does NOT change the
Energy-converged/Structure-converged answer from Case 5 (README section
1.7).
EOF
RUN6="$OUT/case6-quality"
mkdir -p "$RUN6"
cp structure.fdf "$RUN6/"
write_run_out "$RUN6/convergence_meshcutoff_050.0000" "-498.200000" "5.10" \
    "" "0.612500"
write_run_out "$RUN6/convergence_meshcutoff_150.0000" "-500.900000" "5.35"
write_run_out "$RUN6/convergence_meshcutoff_200.0000" "-500.850000" "5.42"
write_run_out "$RUN6/convergence_meshcutoff_250.0000" "-500.845000" "5.45"
write_run_out "$RUN6/convergence_meshcutoff_300.0000" "-500.844000" "5.451"
echo "\$ stb-convergenceAnalysis --dir . --no-intro"
(cd "$RUN6" && stb-convergenceAnalysis --dir . --no-intro \
    | sed -n '/\[2\.1\]/,/Energy converged/p')
pause


echo "=================================================================="
echo " Case 7: Stage 2 -- --apply and --save-gnuplot"
echo "=================================================================="
cat <<'EOF'
--apply writes the recommended value straight into a real calc.fdf, in
place -- the same tag-substitution machinery Stage 1 itself uses.
--save-gnuplot writes a .dat + .gplot pair per parameter (2-panel: energy
per atom, relaxed cell volume -- the same 2 quantities the dual criterion
tracks).
EOF
RUN7="$OUT/case5-worked-example"
cp calc.fdf "$RUN7/my_production.fdf"
echo "\$ stb-convergenceAnalysis --dir . --no-intro --apply my_production.fdf"
(cd "$RUN7" && stb-convergenceAnalysis --dir . --no-intro --apply my_production.fdf \
    | sed -n '/\[5\] APPLY/,$p')
echo
echo "Applied Mesh.CutOff line:"
grep "Mesh.CutOff" "$RUN7/my_production.fdf"
echo
echo "\$ stb-convergenceAnalysis --dir . --no-intro --save-gnuplot -o ."
(cd "$RUN7" && stb-convergenceAnalysis --dir . --no-intro --save-gnuplot -o . \
    | sed -n '/\[3\] OUTPUT FILES/,/\[4\]/p')
pause


echo "=================================================================="
echo " Workflow 4.5 complete"
echo "=================================================================="
cat <<'EOF'
Stage 1 (stb-convergence, 4.5.1): one fully-relaxed folder per swept value
of meshcutoff/energyshift/kgrid (or all three, independently), under
convergence_runs/<parameter>/ -- never a fixed-geometry shortcut, since a
too-coarse setting can itself bias where the structure relaxes to. Stage 2
(stb-convergenceAnalysis, 4.5.2): reads the finished runs back, tracks BOTH
energy and relaxed-volume convergence per parameter, flags per-run data
-quality issues separately (advisory only), and can --apply the recommended
value straight into your own production calc.fdf.

See the README's section 4 for the real, code-verified worked example where
energy and relaxed structure disagree on where convergence happens -- and
section 6 for a step-by-step guide to running this on your own structure.
EOF
