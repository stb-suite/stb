#!/bin/bash
# Guided example: 3.6 -- Status Checker (stb-status; code 3.6 in the stb-suite menu).
#
# Not an automated test (see test/3-utils/6-status/test.sh for that) -- a commented walk-through: it runs real
# commands, one case at a time, on a COPY of project/ inside output/, and shows the piece of the report that proves
# what just happened. project/ itself is never touched. Pauses between cases (`yes '' | bash example_3.6.sh` runs it
# straight through). Safe to re-run -- it always starts by wiping its own output/.
#
# project/ holds real, finished SIESTA 5.4 runs (run_project.sh): a single point with bands, a relaxation that
# converges, one that runs out of moves, an SCF that fails, an MD run, a spin-polarized metal and a series of three
# lattice constants. README section 6 explains every number.

set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"
OUT="$DIR/output"
rm -rf "$OUT"
mkdir -p "$OUT"
cp -a "$DIR/project" "$OUT/project"

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

section() {   # $1 = file, $2 = a section title, exactly: from "[N] TITLE" to the next "[M] "
    strip_colors < "$1" | awk -v t="$2" '/^\[[0-9]+\] [A-Z]/ {if (on) exit; sub(/^\[[0-9]+\] /, "", $0); if ($0 == t) {on = 1; next}} on {print}' \
        | grep -v "^-----"
}

# st <log name> <args...>: runs stb-status inside output/, the report saved as output/<log name>
st() {
    local log=$1; shift
    echo "  \$ stb-status $*"
    (cd "$OUT" && stb-status --no-intro "$@" < /dev/null 2>&1 | strip_colors > "$OUT/$log") || true
}
warnings() { section "$OUT/$1" "WARNINGS" | grep -E "^  [0-9]+\." | fold -s -w 116 | sed 's/^/  /'; }
cut_run() {   # cut_run <dest>: project/relax stopped half-way -- log cut after move 4, no end of run anywhere
    rm -rf "$OUT/$1"; cp -a "$OUT/project/relax" "$OUT/$1"
    local n; n=$(grep -n "Begin CG opt. move =      4" "$OUT/$1/calc.out" | cut -d: -f1)
    head -n $((n + 30)) "$OUT/$1/calc.out" > "$OUT/$1/x" && mv "$OUT/$1/x" "$OUT/$1/calc.out"
    rm -f "$OUT/$1/MESSAGES" "$OUT/$1/0_NORMAL_EXIT"
    grep -v "End of run" "$OUT/$1/CLOCK" > "$OUT/$1/x" && mv "$OUT/$1/x" "$OUT/$1/CLOCK"
}


banner "Welcome: did it finish, did it converge, can I use it?"
cat <<'EOF'
After a batch of SIESTA jobs the questions are always the same: which finished, which failed and why, which are
still running, did the SCF and the relaxation converge, and what came out. stb-status answers them from the files
SIESTA leaves (log, MESSAGES, CLOCK, .EIG, ...) for one folder in full, or for many as one table.

The tour:
   1  a converged relaxation: the whole report    5  running, or stopped?
   2  a single point: results and next tools       6  a series: relative energies
   3  a relaxation that ran out of moves           7  the whole project (--recursive, --csv)
   4  a failed SCF                                 8  the interactive menu
EOF
pause


banner "Case 1 -- A converged relaxation: the whole report"
st case1.txt --path project/relax
section "$OUT/case1.txt" "RUN" | sed 's/^/  /'
section "$OUT/case1.txt" "PROGRESS" | sed 's/^/  /'
section "$OUT/case1.txt" "SCF" | sed 's/^/  /'
cat <<'EOF'
  -> SIESTA itself says whether a relaxation converged: 'Relaxed atomic coordinates' at the end, or 'Final
     (unrelaxed)'. The report shows that, with the last force and stress against MD.MaxForceTol/MD.MaxStressTol,
     and every move: SCF iterations, energy, force, pressure, and the net force (~0 when the mesh is fine).
EOF
pause


banner "Case 2 -- A single point: results and what to do next"
st case2.txt --path project/bulk_sp
section "$OUT/case2.txt" "RESULTS" | sed 's/^/  /'
section "$OUT/case2.txt" "FILES & NEXT STEPS" | sed 's/^/  /'
cat <<'EOF'
  -> Energies, Fermi level, the band gap from the .EIG, forces and pressure; then every analysis output in the
     folder with the suite tool that reads it -- the bands, the eigenvalues, the density.
EOF
pause


banner "Case 3 -- A relaxation that ran out of moves"
st case3.txt --path project/relax_short
section "$OUT/case3.txt" "PROGRESS" | grep -E "Moves|Relaxation|Max force|Max stress" | sed 's/^/  /'
section "$OUT/case3.txt" "SUMMARY & FILES" | grep Verdict | sed 's/^/  /'
warnings case3.txt
cat <<'EOF'
  -> The job ended normally -- SIESTA reached MD.NumCGsteps -- but the structure is not relaxed: the last move even
     overshot (0.93 eV/Ang, 41 meV/atom above the minimum; figures/relax_progress.png). 'FINISHED' is not 'usable'.
EOF
pause


banner "Case 4 -- A failed SCF"
st case4.txt --path project/scf_fail
section "$OUT/case4.txt" "RUN" | grep -E "State|SCF_NOT" | sed 's/^/  /'
section "$OUT/case4.txt" "SCF" | sed 's/^/  /'
grep -n "Tolerance\|MaxSCF" "$OUT/project/scf_fail/calc.fdf" | sed 's/^/  calc.fdf:/'
cat <<'EOF'
  -> The reason is SIESTA's own, from the MESSAGES file. Note the tolerance: calc.fdf asks for 1e-8 at its end, but
     the same label earlier says 1e-5, and SIESTA keeps the FIRST occurrence -- the report shows the value SIESTA
     used. The run failed on MaxSCFIterations 3 alone.
EOF
pause


banner "Case 5 -- Running, or stopped?"
cut_run running
st case5a.txt --path running
section "$OUT/case5a.txt" "RUN" | grep -A1 "State" | sed 's/^/  /'
section "$OUT/case5a.txt" "RESULTS" | sed 's/^/  /' | head -4
cut_run stopped; touch -d "3 hours ago" "$OUT/stopped/calc.out"
st case5b.txt --path stopped
section "$OUT/case5b.txt" "RUN" | grep -A1 "State" | sed 's/^/  /'
cat <<'EOF'
  -> The same log, cut after move 4: written seconds ago it is probably still RUNNING (results partial); written
     hours ago and never ended, it STOPPED -- killed, out of queue time, or crashed without a message.
EOF
pause


banner "Case 6 -- A series: relative energies"
st case6.txt --path 'project/strain/a_*'
section "$OUT/case6.txt" "FOLDERS" | sed 's/^/  /' | cut -c1-140
cat <<'EOF'
  -> With several folders the report is one table. dE is the energy per atom above the lowest folder of the same
     composition and level of theory (XC, spin, mesh, basis): a = 5.476 Ang is the lowest of the three; a parabola
     through them puts the minimum at 5.482 Ang (figures/strain_energy.png), the relaxation at 5.479 Ang.
EOF
pause


banner "Case 7 -- The whole project (--recursive, --csv)"
st case7.txt --path project --recursive --csv project_status.csv
section "$OUT/case7.txt" "FOLDERS" | sed 's/^/  /' | cut -c1-150
warnings case7.txt
echo "  project_status.csv:"; head -3 "$OUT/project_status.csv" | sed 's/^/    /'
cat <<'EOF'
  -> Every run folder below project/, one row each: the failed SCF in red, the unrelaxed relaxation in yellow, the
     MD left out of dE (its energy is not comparable), the spin-polarized Fe compared only with itself. --csv keeps
     the table for a spreadsheet or a script.
EOF
pause


banner "Case 8 -- The interactive menu (stb-suite, code 3.6)"
# path, recursive y, no label, no csv, no report, Enter, 0 = exit
(cd "$OUT" && printf '3.6\nproject\ny\n\n\nn\n\n0\n' | stb-suite 2>&1 | strip_colors > "$OUT/case8.txt") || true
grep -E "^Target|^States" "$OUT/case8.txt" | sed 's/^/  /'
cat <<'EOF'
  -> The menu asks for the folder or glob, recursive, a label, the CSV and the report, and runs the same command.
EOF
pause


banner "Done"
cat <<EOF
Everything is in $OUT/. What to take from the tour:
  - FINISHED is not the same as usable: read the relaxation verdict and the SCF line.
  - A failure carries SIESTA's own reason; an unfinished log is RUNNING? or STOPPED by its age.
  - Many folders: one table, with energies relative to the lowest of the same composition and level of theory.
EOF
