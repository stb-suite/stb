#!/bin/bash
# Guided example: 3.2 -- Clean File Tools (stb-clean; code 3.2 in the stb-suite menu).
#
# Not an automated test (see test/3-utils/2-clean/test.sh for that) -- a commented walk-through: it runs real
# commands, one case at a time, on COPIES of siesta_run/ inside output/<case>/, and shows the piece of output that
# proves what just happened. siesta_run/ itself is never touched. Pauses between cases (`yes '' | bash
# example_3.2.sh` runs it straight through). Safe to re-run -- it always starts by wiping its own output/.
#
# siesta_run/ is a real, finished SIESTA 5.4 run: bulk Si, variable-cell CG relaxation (13 steps), bands, DOS/PDOS
# and the density -- 48 files. Its pseudopotential (Si.psml) is not shipped: each copy gets it from the suite's
# bundled PseudoDojo bank. Case 4 reruns SIESTA when it is on PATH (~25 s); otherwise it shows the reference run.
# README section 6 explains every number.

set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"
OUT="$DIR/output"
rm -rf "$OUT"
mkdir -p "$OUT"
PSML="$(python3 -c 'from stb.core.pseudopotentials import bank_path; print(bank_path("dojo"))')/Si.psml"

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

copy_run() {  # copy_run <dest>: a fresh copy of the finished run + its pseudopotential
    rm -rf "$1"; mkdir -p "$(dirname "$1")"; cp -a "$DIR/siesta_run" "$1"; cp "$PSML" "$1/"
}

# sc <case> <args...>: runs stb-clean inside output/<case>/, the screen report saved as log.txt
sc() {
    local case=$1; shift
    mkdir -p "$OUT/$case"
    echo "  \$ stb-clean $*"
    (cd "$OUT/$case" && stb-clean --no-intro "$@" < /dev/null 2>&1 | strip_colors > log.txt) || true
}
log() { echo "$OUT/$1/log.txt"; }
warnings() { section "$(log "$1")" "WARNINGS" | grep -E "^  [0-9]+\." | fold -s -w 116 | sed 's/^/  /'; }


banner "Welcome: resetting a SIESTA folder for a fresh run"
cat <<'EOF'
A SIESTA run leaves its folder full of files: logs, restart files, bands, densities, basis tables. Before running
the same folder again -- a crashed job, a changed parameter, a convergence test -- you usually want only the inputs
back. stb-clean keeps the extensions in --keep (default .fdf .psf .psml .sh) and removes the rest, after putting
it into a checked tar.gz next to the folder. Its report says what the previous run was, what stays and why, what
goes by kind, and -- the part that matters for the next run -- what SIESTA will start from.

The tour:
   1  what a finished run leaves behind          6  what the report warns about
   2  --dry-run: look before deleting            7  several folders at once (a glob)
   3  the real clean: archive, check, restore    8  a project tree (--recursive)
   4  start fresh or continue (.DM / .XV)        9  the confirmation
   5  files a kept .fdf reads are protected     10  the interactive menu
EOF
pause


banner "Case 1 -- What a finished SIESTA run leaves behind"
copy_run "$OUT/case1/si_run"
sc case1 --path si_run --dry-run
section "$(log case1)" "FOLDER" | sed 's/^/  /'
section "$(log case1)" "REMOVED" | sed 's/^/  /'
cat <<'EOF'
  -> 48 files for a 2-atom cell. The inputs are 5 of them (86 KB); the rest is what SIESTA wrote, sorted by kind:
     restart files it can read back (.DM, .XV, .CG), logs, electronic structure (.HSX, .EIG, .bands, .PDOS),
     grids (.RHO, .VT), and its basis tables (.ion). figures/run_contents.png draws the same table.
EOF
pause


banner "Case 2 -- --dry-run: what stays, and why"
echo "  (the report of Case 1, the same --dry-run)"
section "$(log case1)" "KEPT" | sed 's/^/  /'
section "$(log case1)" "ARCHIVE" | sed 's/^/  /'
ls "$OUT/case1/si_run" | wc -l | sed 's/^/  files still in the folder: /'
cat <<'EOF'
  -> Nothing was touched. bandlines.txt stays although .txt is not in --keep: calc.fdf reads it
     ('%block BandLines < bandlines.txt'), and removing it would break the input (Case 5).
EOF
pause


banner "Case 3 -- The real clean: archive, check, remove; then restore"
copy_run "$OUT/case3/si_run"
sc case3 --path si_run -n --save-report
section "$(log case3)" "ARCHIVE" | sed "s#$OUT/case3/#output/case3/#g; s/^/  /"
echo "  left in the folder:"; ls -A "$OUT/case3/si_run" | sed 's/^/    /'
ARCH=$(ls "$OUT"/case3/si_run_backup_*.tar.gz)
echo "  in the archive: $(tar -tzf "$ARCH" | wc -l) files"
mkdir -p "$OUT/case3/restored" && cp -a "$OUT/case3/si_run/." "$OUT/case3/restored/" \
    && tar -xzf "$ARCH" -C "$OUT/case3/restored" && rm -f "$OUT/case3/restored/stb_clean_report.txt"
if diff -q <(cd "$OUT/case1/si_run" && ls -A | sort) <(cd "$OUT/case3/restored" && ls -A | sort) > /dev/null; then
    echo "  restored (kept + archive) = the original folder, file for file"
fi
cat <<'EOF'
  -> The archive is written next to the folder (not inside it), read back and compared file by file with what
     is about to go; only then is anything removed. Restoring is one tar command.
EOF
pause


banner "Case 4 -- Start fresh, or continue from the previous run"
copy_run "$OUT/case4/fresh"; copy_run "$OUT/case4/continue"
sc case4 --path fresh -n --no-archive
mv "$OUT/case4/log.txt" "$OUT/case4/log_fresh.txt"
sc case4 --path continue -n --no-archive --keep .fdf .psf .psml .sh .DM .XV
mv "$OUT/case4/log.txt" "$OUT/case4/log_continue.txt"
echo "  default --keep:"; section "$OUT/case4/log_fresh.txt" "RESTART" | sed 's/^/    /'
echo "  --keep .fdf .psf .psml .sh .DM .XV:"; section "$OUT/case4/log_continue.txt" "RESTART" | head -6 | sed 's/^/    /'
scf() { python3 - "$1" <<'EOF'
import re, sys
t = open(sys.argv[1]).read()
it = [len(re.findall(r"^\s+scf:\s+\d+", b, re.M)) for b in re.split(r"Begin CG opt\. move", t)[1:]]
e = [l.split("=")[1].strip() for l in t.splitlines() if "siesta:         Total =" in l][-1]
print(f"{len(it)} geometry step(s), {sum(it)} SCF iterations, final E = {e} eV")
EOF
}
if command -v siesta > /dev/null && command -v mpirun > /dev/null; then
    for d in fresh continue; do
        echo "  [siesta] rerun in $d/ ..."
        (cd "$OUT/case4/$d" && OMP_NUM_THREADS=1 mpirun -np 4 siesta calc.fdf < /dev/null > calc.out 2> siesta.err) || true
        echo "    $d: $(scf "$OUT/case4/$d/calc.out")"
    done
else
    echo "  (siesta not on PATH -- the reference reruns of siesta_rerun.sh:)"
    for d in fresh continue; do
        echo "    $d: $(scf "$DIR/reference_results/rerun_$d/calc.out"), wall $(awk '/End of run/ {print $NF}' \
             "$DIR/reference_results/rerun_$d/CLOCK") s"
    done
fi
cat <<'EOF'
  -> Same answer, very different cost. A default clean takes the run back to the .fdf geometry and atomic densities:
     the relaxation is done again (13 steps). Keeping si.XV and si.DM -- which this calc.fdf reads with
     MD.UseSaveXV / DM.UseSaveDM -- the rerun starts at the relaxed geometry and converged density: 1 step, 2 SCF
     iterations (figures/rerun.png). Restarting is a choice: keep them for "continue", drop them for "from scratch".
EOF
pause


banner "Case 5 -- Files a kept .fdf reads are protected"
grep -nE "include|<" "$DIR/siesta_run/calc.fdf" | sed 's/^/  calc.fdf:/'
section "$(log case1)" "KEPT" | grep -E "read by|Protected" | sed 's/^/  /'
copy_run "$OUT/case5/si_run"
echo "notes for myself" > "$OUT/case5/si_run/.notes"
echo "%include missing_options.fdf" >> "$OUT/case5/si_run/calc.fdf"
sc case5 --path si_run --dry-run
section "$(log case5)" "KEPT" | grep -E "notes|read by" | sed 's/^/  /'
warnings case5 | grep -i "missing" || true
cat <<'EOF'
  -> %include F, '%block X < F' and 'Label < F' are followed (also inside included files): what the input reads
     stays, whatever its extension. Hidden files stay too, and a kept .fdf that reads a file which is not there
     is reported -- the rerun would stop on it.
EOF
pause


banner "Case 6 -- What the report warns about"
copy_run "$OUT/case6/si_run"
sc case6 --path si_run --dry-run --keep .psml
echo "  --keep .psml (forgot the rest):"; warnings case6
mv "$OUT/case6/log.txt" "$OUT/case6/log_keep.txt"
copy_run "$OUT/case6/crashed"
head -c 60000 "$DIR/siesta_run/calc.out" > "$OUT/case6/crashed/calc.out"; rm -f "$OUT/case6/crashed/0_NORMAL_EXIT"
sc case6 --path crashed --dry-run
echo "  a calc.out that stops mid-run and was written seconds ago:"
section "$(log case6)" "FOLDER" | grep "Previous run" | sed 's/^/  /'
warnings case6
cat <<'EOF'
  -> No .fdf kept means nothing to rerun from; a file with #SBATCH/#PBS lines or a '#!' first line is flagged as a
     job script (submit.slurm is one: .slurm is not in the default --keep); an unfinished log written in the last
     two minutes may belong to a job that is still running.
EOF
pause


banner "Case 7 -- Several folders at once (a glob)"
for s in m1 0 p1; do copy_run "$OUT/case7/strain_$s"; done
sc case7 --path 'strain_*' -n --keep .fdf .psf .psml .sh .slurm
section "$(log case7)" "SUMMARY" | sed 's/^/  /'
cat <<'EOF'
  -> Each folder is cleaned on its own, with its own archive (strain_0_backup_<time>.tar.gz next to it), and the
     summary collects them. Adding .slurm to --keep silenced the job-script warning.
EOF
pause


banner "Case 8 -- A project tree (--recursive)"
mkdir -p "$OUT/case8/project/.git"
copy_run "$OUT/case8/project/relax"; copy_run "$OUT/case8/project/bands"
echo "ref: refs/heads/main" > "$OUT/case8/project/.git/HEAD"
(cd "$OUT/case8/project" && tar -czf relax_backup_20260101_120000.tar.gz relax/calc.out)
sc case8 --path project --recursive --dry-run
section "$(log case8)" "FOLDER" | sed 's/^/  /'
section "$(log case8)" "KEPT" | grep -E "backup|hidden|read by" | sed 's/^/  /'
cat <<'EOF'
  -> Subfolders are cleaned too; hidden folders (.git) are never entered, and an earlier stb-clean backup inside
     the tree is never removed.
EOF
pause


banner "Case 9 -- The confirmation"
copy_run "$OUT/case9/si_run"
echo "  \$ stb-clean --path si_run          (answering n)"
(cd "$OUT/case9" && printf 'n\n' | stb-clean --no-intro --path si_run 2>&1 | strip_colors > log.txt) || true
grep -E "will be removed|Delete these|Cancelled" "$(log case9)" | sed 's/^/  /'
echo "  files still there: $(ls "$OUT/case9/si_run" | wc -l)"
cat <<'EOF'
  -> One question per folder, after the report. -n skips it (for scripts); with no terminal to answer, the
     question reads end-of-file and cancels -- never deletes by default.
EOF
pause


banner "Case 10 -- The interactive menu (stb-suite, code 3.2)"
copy_run "$OUT/case10/si_run"
# path, extensions (Enter = default), recursive? n, dry run? y, save report? n, Enter, 0 = exit
printf '3.2\n%s\n\nn\ny\nn\n\n0\n' "$OUT/case10/si_run" | (cd "$OUT/case10" && stb-suite > log.txt 2>&1) || true
strip_colors < "$OUT/case10/log.txt" | grep -E "^Mode|^Dry run --|^Would remove" | sed 's/^/  /'
cat <<'EOF'
  -> The menu asks the same questions (folder or glob, extensions, recursive, dry run, archive, confirmation,
     report) and runs the same command; with a dry run it skips the archive and confirmation questions.
EOF
pause


banner "Done"
cat <<EOF
Everything is in $OUT/. What to take from the tour:
  - A SIESTA folder is a few inputs and many outputs; stb-clean keeps the inputs (and whatever they read), puts the
    rest in a checked tar.gz next to the folder, and removes it.
  - Read the RESTART section: removing .DM/.XV turns a rerun into a run from scratch; keeping them (with the
    UseSave flags on) continues where the last run stopped.
  - --dry-run first on anything you care about; the warnings catch a forgotten --keep, a job script, a missing
    include and a job that may still be running.
EOF
