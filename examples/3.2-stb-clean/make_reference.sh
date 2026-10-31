#!/bin/bash
# =====================================================================================================
#  make_reference.sh -- regenerates reference_results/<case>/report.txt of the 3.2 example (stb-clean)
#
#  Every case works on a fresh copy of siesta_run/ (a real, finished SIESTA run) in a scratch folder, so
#  siesta_run/ itself is never touched; only the reports are kept (scratch paths written as <work>).
#  Seconds, no SIESTA. The SIESTA reruns of reference_results/rerun_*/ come from siesta_rerun.sh.
#
#      source ~/.anaconda3/etc/profile.d/conda.sh && conda activate stb    # any environment with stb_suite
#      bash make_reference.sh
# =====================================================================================================
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
REF="$HERE/reference_results"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
PSML="$(python3 -c 'from stb.core.pseudopotentials import bank_path; print(bank_path("dojo"))')/Si.psml"

fresh_copy() {   # fresh_copy <dest>: the finished run + its pseudopotential (from the suite's bundled bank)
    rm -rf "$1"; cp -a "$HERE/siesta_run" "$1"; cp "$PSML" "$1/"
}
report() {       # report <case> <command...>: runs in $WORK, stores the plain-text screen output
    local case="$1"; shift
    mkdir -p "$REF/$case"
    echo "  [run]  $case: $*"
    (cd "$WORK" && "$@" 2>&1) | sed 's/\x1b\[[0-9;]*m//g' | sed "s#$WORK#<work>#g" > "$REF/$case/report.txt" || true
}
listing() {      # listing <case> <folder>: what is left in a folder after the case
    (cd "$WORK/$2" && ls -A) > "$REF/$1/files_after.txt"
}

echo "== 1. What would go -- --dry-run"
fresh_copy "$WORK/si_run"
report dry_run stb-clean --no-intro --path si_run --dry-run

echo "== 2. The real clean: archive, check, remove; then restore from the archive"
report clean stb-clean --no-intro --path si_run -n --save-report
listing clean si_run
tar -tzf "$WORK"/si_run_backup_*.tar.gz | sort > "$REF/clean/archive_contents.txt"
mkdir "$WORK/restored" && cp -a "$WORK/si_run/." "$WORK/restored/" && tar -xzf "$WORK"/si_run_backup_*.tar.gz -C "$WORK/restored"
(cd "$WORK/restored" && rm -f stb_clean_report.txt && ls -A | sort) > "$REF/clean/restored_files.txt"
rm -f "$WORK"/si_run_backup_*.tar.gz

echo "== 3. Continue instead of starting fresh: keep .DM and .XV"
fresh_copy "$WORK/si_run"
report keep_restart stb-clean --no-intro --path si_run -n --no-archive --keep .fdf .psf .psml .sh .DM .XV
listing keep_restart si_run

echo "== 4. What the report warns about"
fresh_copy "$WORK/si_run"
report careless_keep stb-clean --no-intro --path si_run --dry-run --keep .psml
fresh_copy "$WORK/crashed"
head -c 60000 "$HERE/siesta_run/calc.out" > "$WORK/crashed/calc.out"; rm -f "$WORK/crashed/0_NORMAL_EXIT"
report running stb-clean --no-intro --path crashed --dry-run

echo "== 4b. Files a kept .fdf reads, hidden files, an include that is missing"
fresh_copy "$WORK/protected"
echo "notes for myself" > "$WORK/protected/.notes"
echo "%include missing_options.fdf" >> "$WORK/protected/calc.fdf"
report protected stb-clean --no-intro --path protected --dry-run

echo "== 5. Several folders at once (a glob)"
for s in m1 0 p1; do fresh_copy "$WORK/strain_$s"; done
report batch stb-clean --no-intro --path 'strain_*' -n --keep .fdf .psf .psml .sh .slurm

echo "== 6. A project tree (--recursive): hidden folders and earlier backups are left alone"
rm -rf "$WORK/project"; mkdir -p "$WORK/project/.git"
fresh_copy "$WORK/project/relax"; fresh_copy "$WORK/project/bands"
echo "ref: refs/heads/main" > "$WORK/project/.git/HEAD"
echo "notes on this project" > "$WORK/project/.notes"
(cd "$WORK/project" && tar -czf relax_backup_20260101_120000.tar.gz relax/calc.out)
report recursive stb-clean --no-intro --path project --recursive --dry-run

echo "== 7. Saying no at the confirmation"
fresh_copy "$WORK/si_run"
mkdir -p "$REF/cancel"
(cd "$WORK" && printf 'n\n' | stb-clean --no-intro --path si_run 2>&1) | sed 's/\x1b\[[0-9;]*m//g' \
    | sed "s#$WORK#<work>#g" > "$REF/cancel/report.txt" || true
listing cancel si_run

echo "done -- reference_results/<case>/report.txt"
