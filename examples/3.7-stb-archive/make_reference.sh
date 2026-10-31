#!/bin/bash
# =====================================================================================================
#  make_reference.sh -- regenerates reference_results/<case>/ of the 3.7 example (stb-archive)
#
#  Every case works on a copy of runs/ (real SIESTA 5.4 runs; run_example.sh) in a scratch folder, with the
#  pseudopotential copied in from the suite's bundled bank (an input the archive must carry). Only reports and
#  listings are kept; scratch paths become <work>, dates <date>, checksums of the archive itself <sha>. Seconds.
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
fresh() {        # fresh <dest> <run>: a copy of runs/<run> with its pseudopotential
    rm -rf "$WORK/$1"; mkdir -p "$(dirname "$WORK/$1")"; cp -a "$HERE/runs/$2" "$WORK/$1"
    local f; for f in $(find "$WORK/$1" -name calc.fdf); do cp "$PSML" "$(dirname "$f")/"; done
}
clean() { sed 's/\x1b\[[0-9;]*m//g' | sed "s#$WORK#<work>#g" \
          | sed -E 's/^(Date\/time *: ).*/\1<date>/; s/^(SHA-256 \(archive\) : ).*/\1<sha>/'; }
report() {       # report <case> <stb-archive args...>: runs in $WORK, stores the plain-text report
    local case="$1"; shift
    mkdir -p "$REF/$case"
    echo "  [run]  $case: stb-archive $*"
    (cd "$WORK" && stb-archive --no-intro "$@" 2>&1) | clean > "$REF/$case/report.txt" || true
}

echo "== 1. What would go in: --dry-run"
fresh si_relax si_relax
report dry_run --path si_relax --dry-run

echo "== 2. The archive, with its subfolder; then restore and check"
report archive --path si_relax --recursive
tar -tzf "$WORK/si_relax_archive.tar.gz" > "$REF/archive/members.txt"
mkdir -p "$WORK/restore" && tar -xzf "$WORK/si_relax_archive.tar.gz" -C "$WORK/restore"
clean < "$WORK/restore/si_relax/MANIFEST.txt" | sed -E 's/^(Created *: ).*/\1<date>/' > "$REF/archive/MANIFEST.txt"
(cd "$WORK/restore/si_relax" && sha256sum -c SHA256SUMS 2>&1 | LC_ALL=C sed 's/: SUCESSO$/: OK/') > "$REF/archive/sha256sum_check.txt"
echo "tampered" >> "$WORK/restore/si_relax/si.EIG"
(cd "$WORK/restore/si_relax" && LC_ALL=C sha256sum -c --quiet SHA256SUMS 2>&1 || true) > "$REF/archive/sha256sum_tampered.txt"

echo "== 3. The three modes, and --include"
for m in inputs results full; do report "mode_$m" --path si_relax --recursive --mode "$m" --dry-run; done
report include_rho --path si_relax --recursive --include .RHO --dry-run

echo "== 4. Never overwritten"
report second --path si_relax --recursive
report output_exists --path si_relax -o si_relax_archive.tar.gz

echo "== 5. What the report catches"
fresh nopsml si_relax; rm "$WORK/nopsml/Si.psml"
report guard_nopsml --path nopsml --dry-run
fresh noref si_relax; rm "$WORK/noref/bandlines.txt"
report guard_noref --path noref --dry-run
fresh unfinished si_relax
n=$(grep -n "Begin CG opt. move =      4" "$WORK/unfinished/calc.out" | cut -d: -f1)
head -n $((n + 30)) "$WORK/unfinished/calc.out" > "$WORK/x" && mv "$WORK/x" "$WORK/unfinished/calc.out"
rm -f "$WORK/unfinished/MESSAGES" "$WORK/unfinished/0_NORMAL_EXIT"
grep -v "End of run" "$WORK/unfinished/CLOCK" > "$WORK/x" && mv "$WORK/x" "$WORK/unfinished/CLOCK"
touch -d "3 hours ago" "$WORK/unfinished/calc.out"
report guard_unfinished --path unfinished --dry-run

echo "== 6. A series, one archive per folder"
for a in 5.40 5.476 5.55; do fresh "strain/a_$a" "strain/a_$a"; done
report series --path 'strain/a_*'

echo "== 7. The interactive menu (stb-suite, code 3.7)"
fresh menu/si_relax si_relax; mkdir -p "$REF/menu"
# path, mode 2 = inputs, no recursive, no include, no exclude, dry run y, no report, Enter, 0 = exit
(cd "$WORK/menu" && printf '3.7\nsi_relax\n2\nn\n\n\ny\nn\n\n0\n' | stb-suite 2>&1) | clean \
    | sed -n '/STB-ARCHIVE REPORT/,/science preserved/p' > "$REF/menu/report.txt"
echo "done -- reference_results/<case>/"
