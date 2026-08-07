#!/bin/bash
# Guided example: Hubbard U (Linear Response) workflow
# (stb-hubbardu / stb-hubbarduAlphas / stb-hubbarduAnalysis, codes 4.7.1/4.7.2/4.7.3)
#
# Not an automated test (see test/4-workflow/7-hubbardu/{prep,alphas,analysis}/
# test.sh for that) -- a commented walk-through: it runs real commands, one
# group at a time, and shows you the piece of output that proves what just
# happened. It pauses between sections so you can read before moving on.
#
# No SIESTA binary is invoked anywhere in this script (none of the three
# tools ever runs SIESTA themselves). Stages 1/2 are exercised for real
# (they only write input files). Stage 3 needs real SIESTA .out files to
# analyze, which this walkthrough doesn't have -- so the "full workflow"
# cases fabricate calc.out files whose occupation-vs-alpha response follows
# an EXACT, hand-chosen linear law (see the README's Section 6), proving
# the fitting/formula pipeline against a known answer instead of just
# printing a number and hoping it's right.

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

# Writes a synthetic calc.out into $1/$2 whose Occupations: line follows
# n(alpha) = $4 + $3*alpha exactly (a 3-number, already-total line -- no
# non-polarized doubling needed, see README Section 5.2).
write_occ() {
    local run_dir="$1" folder="$2" alpha="$3" chi="$4" n0="$5"
    local total half
    total=$(python3 -c "print(f'{$n0 + $chi*$alpha:.6f}')")
    half=$(python3 -c "print(f'{($n0 + $chi*$alpha)/2:.6f}')")
    cat > "$run_dir/$folder/calc.out" <<EOF
Siesta Version  : synthetic (tutorial fixture -- no real SIESTA run)
Occupations:   $half   $half   $total
SCF Convergence by DM+H criterion
EOF
}

# Populates every reference/scf_alpha_*/frozen_alpha_* folder under $1
# with a synthetic linear response: chi (screened) for scf+reference,
# chi0 (bare) for frozen, both around n0.
# Folder-name suffix for an alpha value -- must match stb-hubbarduAlphas'
# own f"{alpha:.4f}" naming exactly, hence python3 (not bash printf, whose
# %.4f is locale-sensitive -- e.g. uses a comma under a pt_BR locale).
alpha_suffix() { python3 -c "print(f'{$1:.4f}')"; }

populate_synthetic_response() {
    local run_dir="$1" chi="$2" chi0="$3" n0="$4"
    write_occ "$run_dir" "reference" 0.0 "$chi" "$n0"
    for a in -0.10 -0.05 0.05 0.10; do
        write_occ "$run_dir" "scf_alpha_$(alpha_suffix "$a")" "$a" "$chi" "$n0"
    done
    for a in 0.0 -0.10 -0.05 0.05 0.10; do
        write_occ "$run_dir" "frozen_alpha_$(alpha_suffix "$a")" "$a" "$chi0" "$n0"
    done
}

echo "=================================================================="
echo " Why this workflow needs three stages, not two"
echo "=================================================================="
cat <<'EOF'
DFT+U's own correction term needs a number, U -- the curvature of the
total energy with respect to a correlated shell's occupation. Cococcioni &
de Gironcoli (PRB 71, 035105, 2005) compute that curvature INDIRECTLY: by
measuring how the shell's SCF occupation responds to a small rigid
potential shift, alpha, applied only to that shell.

Two response slopes matter:
  chi  -- the SELF-CONSISTENT (screened) response: alpha is applied and
          the whole SCF is allowed to relax fully around it.
  chi0 -- the FROZEN-DENSITY (bare) response: alpha is applied but the SCF
          is capped at ~1 iteration, before the rest of the density has a
          chance to relax.

  U = 1/chi0 - 1/chi

Stage 1 (stb-hubbardu) prepares the alpha=0 reference point. Stage 2
(stb-hubbarduAlphas) prepares every other alpha, seeded from Stage 1's own
converged density. Stage 3 (stb-hubbarduAnalysis) fits chi/chi0 and
applies the formula.
EOF
pause

echo "=================================================================="
echo " output/stage1/  --  Stage 1 mechanics on a single-atom structure"
echo "=================================================================="
echo "structure.fdf is a simplified single-atom bcc-Mn cell."
echo
echo "\$ stb-hubbardu --species Mn -o output/stage1 --no-intro"
stb-hubbardu --species Mn -o "$OUT/stage1" --no-intro > "$OUT/stage1_console.log"
sed -n '/\[2\] DFT+U PERTURBATION SETUP/,/\[3\]/p' "$OUT/stage1_console.log" | head -n -1
echo
echo "Folder contents:"
ls "$OUT/stage1/reference"
pause

echo "=================================================================="
echo " output/stage1_2mn/  --  --atom-index and species aliasing"
echo "=================================================================="
echo "structure_2mn.fdf: conventional 2-atom bcc-Mn cell, both Mn atoms"
echo "symmetry-equivalent (space group Im-3m). No --atom-index given:"
echo
echo "\$ stb-hubbardu -s structure_2mn.fdf --species Mn -o output/stage1_2mn --no-intro"
set +e
stb-hubbardu -s structure_2mn.fdf --species Mn -o "$OUT/stage1_2mn" --no-intro > "$OUT/stage1_2mn_ambiguous.log" 2>&1
AMBIGUOUS_EXIT=$?
set -e
echo "(exit code: $AMBIGUOUS_EXIT)"
grep -E "FAIL|Space group|Wyckoff|symmetry-equivalent" "$OUT/stage1_2mn_ambiguous.log"
echo
echo "Now resolved with --atom-index 1:"
echo "\$ stb-hubbardu -s structure_2mn.fdf --species Mn --atom-index 1 -o output/stage1_2mn --no-intro"
stb-hubbardu -s structure_2mn.fdf --species Mn --atom-index 1 -o "$OUT/stage1_2mn" --no-intro > "$OUT/stage1_2mn_console.log"
grep "Perturbed atom" "$OUT/stage1_2mn_console.log"
echo
echo "The aliased structure.fdf -- one new same-Z species, one atom repointed:"
grep -A3 "ChemicalSpeciesLabel" "$OUT/stage1_2mn/reference/structure.fdf" | head -4
pause

echo "=================================================================="
echo " output/stage1_pp/  --  pseudopotentials from a bundled bank"
echo "=================================================================="
echo "\$ stb-hubbardu --species Mn --pseudo-dir dojo -o output/stage1_pp --no-intro"
stb-hubbardu --species Mn --pseudo-dir dojo -o "$OUT/stage1_pp" --no-intro > "$OUT/stage1_pp_console.log"
grep -A1 "PSEUDOPOTENTIALS" "$OUT/stage1_pp_console.log" | tail -1
ls "$OUT/stage1_pp/reference"
pause

echo "=================================================================="
echo " output/workflow/  --  the full 3-stage chain, a known-answer proof"
echo "=================================================================="
echo "Stage 1 (reference) + Stage 2 (perturbations), for real:"
mkdir -p "$OUT/workflow"
stb-hubbardu --species Mn -o "$OUT/workflow/hubbardu_runs" --no-intro > "$OUT/workflow_stage1.log"
echo "  \$ stb-hubbardu --species Mn -o hubbardu_runs --no-intro"
echo "This walkthrough never ran real SIESTA, so reference/siesta.DM is"
echo "simulated with a placeholder (a real run writes a real binary .DM):"
echo "siesta.DM placeholder -- simulates a converged SIESTA run" > "$OUT/workflow/hubbardu_runs/reference/siesta.DM"
stb-hubbarduAlphas --dir "$OUT/workflow/hubbardu_runs" --no-intro > "$OUT/workflow_stage2.log"
echo "  \$ stb-hubbarduAlphas --dir hubbardu_runs --no-intro"
grep "Success" "$OUT/workflow_stage2.log"
echo
echo "Now the synthetic ground truth: n0=5.0, chi=0.30 (screened), chi0=0.20"
echo "(bare) -- so the expected answer is known BEFORE running Stage 3:"
echo "  U = 1/chi0 - 1/chi = 1/0.20 - 1/0.30 = 1.666667 eV"
populate_synthetic_response "$OUT/workflow/hubbardu_runs" 0.30 0.20 5.0
echo
echo "\$ stb-hubbarduAnalysis --dir hubbardu_runs --no-intro"
stb-hubbarduAnalysis --dir "$OUT/workflow/hubbardu_runs" --no-intro > "$OUT/workflow_stage3.log"
sed -n '/\[2\] LINEAR RESPONSE FIT/,/\[4\]/p' "$OUT/workflow_stage3.log" | head -n -1
EXPECTED_U="1.666667"
GOT_U=$(grep "Computed U" "$OUT/workflow_stage3.log" | head -1 | grep -oE '[0-9]+\.[0-9]+')
echo
if python3 -c "import sys; sys.exit(0 if abs($GOT_U - $EXPECTED_U) < 0.001 else 1)"; then
    echo "Confirmed: computed U ($GOT_U eV) matches the hand-derived $EXPECTED_U eV."
else
    echo "Unexpected: computed U ($GOT_U eV) does NOT match the expected $EXPECTED_U eV."
fi
pause

echo "=================================================================="
echo " output/workflow_badU/  --  what a bad fit looks like: negative U"
echo "=================================================================="
echo "Same chain, chi and chi0 swapped (screened response now SMALLER than"
echo "the bare one) -- the formula is applied exactly the same way, and"
echo "the result is a physically impossible negative U:"
mkdir -p "$OUT/workflow_badU"
stb-hubbardu --species Mn -o "$OUT/workflow_badU/hubbardu_runs" --no-intro > /dev/null
echo "placeholder" > "$OUT/workflow_badU/hubbardu_runs/reference/siesta.DM"
stb-hubbarduAlphas --dir "$OUT/workflow_badU/hubbardu_runs" --no-intro > /dev/null
populate_synthetic_response "$OUT/workflow_badU/hubbardu_runs" 0.20 0.30 5.0
stb-hubbarduAnalysis --dir "$OUT/workflow_badU/hubbardu_runs" --no-intro > "$OUT/workflow_badU_stage3.log"
sed -n '/\[3\] COMPUTED U/,/\[4\]/p' "$OUT/workflow_badU_stage3.log" | head -n -1
pause

echo "=================================================================="
echo " output/workflow_noisy/  --  what a bad fit looks like: low R^2"
echo "=================================================================="
echo "Same clean dataset as output/workflow/, except ONE scf_alpha point"
echo "is replaced with an outlier -- everything else is untouched:"
mkdir -p "$OUT/workflow_noisy"
stb-hubbardu --species Mn -o "$OUT/workflow_noisy/hubbardu_runs" --no-intro > /dev/null
echo "placeholder" > "$OUT/workflow_noisy/hubbardu_runs/reference/siesta.DM"
stb-hubbarduAlphas --dir "$OUT/workflow_noisy/hubbardu_runs" --no-intro > /dev/null
populate_synthetic_response "$OUT/workflow_noisy/hubbardu_runs" 0.30 0.20 5.0
write_occ "$OUT/workflow_noisy/hubbardu_runs" "scf_alpha_-0.0500" -0.05 0.30 5.0  # placeholder, overwritten next
cat > "$OUT/workflow_noisy/hubbardu_runs/scf_alpha_-0.0500/calc.out" <<'EOF'
Occupations:   2.350000   2.350000   4.700000
SCF Convergence by DM+H criterion
EOF
stb-hubbarduAnalysis --dir "$OUT/workflow_noisy/hubbardu_runs" --no-intro > "$OUT/workflow_noisy_stage3.log" || true
sed -n '/\[2\] LINEAR RESPONSE FIT/,/\[3\]/p' "$OUT/workflow_noisy_stage3.log" | head -n -1
pause

echo "=================================================================="
echo " Proof: CLI and the interactive stb-suite menu agree"
echo "=================================================================="
echo "Driving the same single-atom Mn case through all three interactive"
echo "menu stages (non-interactively, via piped printf) and comparing the"
echo "final computed U against output/workflow/'s direct-CLI run."
TMP="$(mktemp -d)"
cp structure.fdf calc.fdf "$TMP/"
echo
echo "\$ printf '4.7.1\\nstructure.fdf\\ncalc.fdf\\nMn\\n\\n\\n\\n\\nhubbardu_runs\\nn\\n\\n0\\n' | stb-suite"
(cd "$TMP" && printf '4.7.1\nstructure.fdf\ncalc.fdf\nMn\n\n\n\n\nhubbardu_runs\nn\n\n0\n' | stb-suite > menu1.log 2>&1)
echo "placeholder" > "$TMP/hubbardu_runs/reference/siesta.DM"
echo "\$ printf '4.7.2\\nhubbardu_runs\\n\\n\\nn\\n\\n0\\n' | stb-suite"
(cd "$TMP" && printf '4.7.2\nhubbardu_runs\n\n\nn\n\n0\n' | stb-suite > menu2.log 2>&1)
populate_synthetic_response "$TMP/hubbardu_runs" 0.30 0.20 5.0
echo "\$ printf '4.7.3\\nhubbardu_runs\\n\\n\\n\\nn\\nn\\nn\\n\\n0\\n' | stb-suite"
(cd "$TMP" && printf '4.7.3\nhubbardu_runs\n\n\n\nn\nn\nn\n\n0\n' | stb-suite > menu3.log 2>&1)
CLI_U=$(grep "Computed U" "$OUT/workflow_stage3.log" | head -1)
MENU_U=$(grep "Computed U" "$TMP/menu3.log" | head -1)
if [ "$CLI_U" = "$MENU_U" ]; then
    echo "Confirmed: identical result ($MENU_U) from the CLI and the interactive menu."
else
    echo "Unexpected: results differ -- CLI: $CLI_U / menu: $MENU_U"
fi
rm -rf "$TMP"
rm -f Mn_LDAU.fdf
pause

echo "=================================================================="
echo " Done"
echo "=================================================================="
cat <<'EOF'
Folders generated under output/:
  stage1/          stage1_2mn/       stage1_pp/
  workflow/         workflow_badU/    workflow_noisy/

output/workflow/hubbardu_runs/ has the full chain: reference/,
scf_alpha_*/, frozen_alpha_*/, run_manifest.json, and (via Stage 3)
Mn_LDAU.fdf -- the production-ready %block LDAU.proj, no
LDAU.PotentialShift line.

As a next step, on your OWN structure:
  stb-hubbardu --species <El> --pseudo-dir dojo
  # run SIESTA in hubbardu_runs/reference/, then:
  stb-hubbarduAlphas --dir hubbardu_runs
  # run SIESTA in every scf_alpha_*/frozen_alpha_* folder, then:
  stb-hubbarduAnalysis --dir hubbardu_runs --save-gnuplot --view

Then paste the %block LDAU.proj from <species>_LDAU.fdf into your
production calc.fdf, or hand it to stb-inputfile (example 1.1).
EOF
