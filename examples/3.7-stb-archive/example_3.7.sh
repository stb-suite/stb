#!/bin/bash
# Guided example: 3.7 -- Archive Packager (stb-archive; code 3.7 in the stb-suite menu).
#
# Not an automated test (see test/3-utils/7-archive/test.sh for that) -- a commented walk-through: it runs real
# commands, one case at a time, on COPIES of runs/ inside output/<case>/, and shows the piece of the report that proves
# what just happened. runs/ itself is never touched. Pauses between cases (`yes '' | bash example_3.7.sh` runs it
# straight through). Safe to re-run -- it always starts by wiping its own output/.
#
# runs/si_relax is a real SIESTA 5.4 run kept exactly as SIESTA left it (relaxation + bands + DOS/PDOS + density +
# H/S, a job script, an %include'd structure, a block read from bandlines.txt, an analysis/ subfolder); runs/strain/
# three single points. The pseudopotential (an input the archive must carry) is copied in from the suite's bundled
# bank. Case 4 reruns SIESTA from an archive when SIESTA is installed. README section 6 explains every number.

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

section() {   # $1 = file, $2 = a section title, exactly (or its start, for 'FOLDER <path>')
    strip_colors < "$1" | awk -v t="$2" '/^\[[0-9]+\] [A-Z]/ {if (on) exit; sub(/^\[[0-9]+\] /, "", $0); if ($0 == t || index($0, t " ") == 1) {on = 1; next}} on {print}' \
        | grep -v "^-----"
}

fresh() {     # fresh <dest> <run>: a copy of runs/<run> with its pseudopotential
    rm -rf "$1"; mkdir -p "$(dirname "$1")"; cp -a "$DIR/runs/$2" "$1"
    local f; for f in $(find "$1" -name calc.fdf); do cp "$PSML" "$(dirname "$f")/"; done
}

# ar <case> <args...>: runs stb-archive inside output/<case>/, the report saved as log.txt there
ar() {
    local case=$1; shift
    mkdir -p "$OUT/$case"
    echo "  \$ stb-archive $*"
    (cd "$OUT/$case" && stb-archive --no-intro "$@" < /dev/null 2>&1 | strip_colors > log.txt) || true
}
log() { echo "$OUT/$1/log.txt"; }
warnings() { section "$(log "$1")" "WARNINGS" | grep -E "^  [0-9]+\." | fold -s -w 116 | sed 's/^/  /'; }


banner "Welcome: packaging a calculation for sharing, or for the record"
cat <<'EOF'
A finished SIESTA folder holds a few kilobytes of inputs, a megabyte of results and several megabytes of files SIESTA
can regenerate. stb-archive sorts every file by kind, puts the inputs and the results into a .tar.gz, leaves the rest
out with a reason, adds a MANIFEST.txt (what the run was, every file with its SHA-256) and a SHA256SUMS, and reads the
archive back before saying ok. It is the complement of stb-clean (3.2), which removes outputs to prepare a rerun.

The tour:
   1  what would go in (--dry-run)                5  never overwritten
   2  the archive, the manifest, the checksums    6  what the report catches
   3  the three modes                             7  a series: one archive per folder
   4  an inputs-only archive reruns the result    8  the interactive menu
EOF
pause


banner "Case 1 -- What would go in: --dry-run"
fresh "$OUT/case1/si_relax" si_relax
ar case1 --path si_relax --dry-run
section "$(log case1)" "FOLDER" | sed 's/^/  /'
section "$(log case1)" "GOING IN" | sed 's/^/  /'
section "$(log case1)" "LEFT OUT" | sed 's/^/  /' | cut -c1-140
cat <<'EOF'
  -> Inputs first: both .fdf files, bandlines.txt (calc.fdf reads it), the pseudopotential and the job script. Then
     the results: logs, geometry, bands, DOS, PDOS, eigenvalues. Left out: the basis tables, the H/S binary, the
     density matrix and the grids -- SIESTA regenerates them, and they are two thirds of the folder. analysis/ is a
     subfolder: --recursive takes it.
EOF
pause


banner "Case 2 -- The archive, the manifest and the checksums"
fresh "$OUT/case2/si_relax" si_relax
ar case2 --path si_relax --recursive
section "$(log case2)" "ARCHIVE" | sed "s#$OUT/case2/##; s/^/  /"
mkdir -p "$OUT/case2/restore" && tar -xzf "$OUT/case2/si_relax_archive.tar.gz" -C "$OUT/case2/restore"
echo "  MANIFEST.txt (start):"; sed -n '/THE RUN/,/Wall time/p' "$OUT/case2/restore/si_relax/MANIFEST.txt" | sed 's/^/    /'
echo "  \$ cd restore/si_relax && sha256sum -c SHA256SUMS"
(cd "$OUT/case2/restore/si_relax" && LC_ALL=C sha256sum -c SHA256SUMS | tail -3 | sed 's/^/    /')
echo "tampered" >> "$OUT/case2/restore/si_relax/si.EIG"
echo "  ... after changing si.EIG:"
(cd "$OUT/case2/restore/si_relax" && LC_ALL=C sha256sum -c --quiet SHA256SUMS 2>&1 | sed 's/^/    /') || true
cat <<'EOF'
  -> 38 files, 1.0 MB, compress to 269 KB. The archive is checked file by file before the report says ok; after
     extraction, sha256sum -c checks it again -- a changed file is caught. The manifest says what the run was, so the
     archive explains itself years later.
EOF
pause


banner "Case 3 -- The three modes"
for m in inputs results full; do
    fresh "$OUT/case3_$m/si_relax" si_relax
    ar "case3_$m" --path si_relax --recursive --mode "$m" --dry-run > /dev/null
    echo "  --mode $m: $(grep '^Going in' "$(log "case3_$m")" | sed 's/^Going in *: //')"
done
cat <<'EOF'
  -> inputs: only what reruns the calculation (86 KB); results (default): what you would want to look at again;
     full: everything. --include .RHO adds the density to any of them (figures/contents.png).
EOF
pause


banner "Case 4 -- Does an inputs-only archive really reproduce the run?"
if command -v siesta > /dev/null && command -v mpirun > /dev/null; then
    bash "$DIR/rerun_from_archive.sh" 4 > /dev/null 2>&1 && sed 's/^/  /' "$DIR/reference_results/rerun/comparison.txt"
else
    echo "  (siesta not on PATH -- the reference rerun of rerun_from_archive.sh:)"
    sed 's/^/  /' "$DIR/reference_results/rerun/comparison.txt"
fi
cat <<'EOF'
  -> Unpacked in an empty folder and run again, the 86 KB inputs-only archive gives the same relaxation: 8 CG moves
     and the same energy to the last digit. 'Rerunnable: yes' in the report is checked, not assumed.
EOF
pause


banner "Case 5 -- Never overwritten"
ar case2 --path si_relax --recursive > /dev/null
grep "^Name" "$(log case2)" | sed 's/^/  /'
ar case2 --path si_relax -o si_relax_archive.tar.gz > /dev/null
grep "^Not written" "$(log case2)" | sed 's/^/  /'
ls "$OUT/case2" | grep tar.gz | sed 's/^/  in the folder: /'
cat <<'EOF'
  -> A second archive of the same folder gets _1, _2, ...; an -o that exists is refused unless --force.
EOF
pause


banner "Case 6 -- What the report catches"
fresh "$OUT/case6/nopsml" si_relax; rm "$OUT/case6/nopsml/Si.psml"
ar case6 --path nopsml --dry-run; echo "  without the pseudopotential:"; warnings case6
fresh "$OUT/case6/noref" si_relax; rm "$OUT/case6/noref/bandlines.txt"
ar case6 --path noref --dry-run; echo "  without the file calc.fdf reads:"; warnings case6
cat <<'EOF'
  -> 'Rerunnable' is checked: every species needs a pseudopotential, every file an .fdf reads must be there. A run
     that did not finish, or a relaxation that did not converge, is flagged too (stb-status's verdict).
EOF
pause


banner "Case 7 -- A series: one archive per folder"
for a in 5.40 5.476 5.55; do fresh "$OUT/case7/strain/a_$a" "strain/a_$a"; done
ar case7 --path 'strain/a_*'
section "$(log case7)" "SUMMARY & FILES" | head -6 | sed 's/^/  /'
cat <<'EOF'
  -> Quote the glob. Each folder gets its own archive next to it, checked, and the summary collects them.
EOF
pause


banner "Case 8 -- The interactive menu (stb-suite, code 3.7)"
fresh "$OUT/case8/si_relax" si_relax
# path, mode 2 = inputs, no recursive, no include, no exclude, dry run y, no report, Enter, 0 = exit
(cd "$OUT/case8" && printf '3.7\nsi_relax\n2\nn\n\n\ny\nn\n\n0\n' | stb-suite 2>&1 | strip_colors > log.txt) || true
grep -E "^Mode  |^Going in|^Rerunnable|^Dry run  " "$OUT/case8/log.txt" | sed 's/^/  /' | cut -c1-140
cat <<'EOF'
  -> The menu asks for the folder or glob, the mode, subfolders, extra suffixes in or out, dry run and report.
EOF
pause


banner "Done"
cat <<EOF
Everything is in $OUT/. What to take from the tour:
  - Inputs and results go in, regenerable bulk stays out, each file with its kind and reason.
  - The archive explains itself (MANIFEST.txt) and checks itself (SHA256SUMS); it is verified before 'ok'.
  - --mode inputs is enough to reproduce the run exactly; the report checks that it can.
EOF
