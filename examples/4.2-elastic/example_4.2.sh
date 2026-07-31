#!/bin/bash
# Guided example: Workflow 4.2 -- Elastic Constants (both stages:
# stb-elasticInputs, code 4.2.1, and stb-elasticAnalysis, code 4.2.2, in
# the stb-suite menu).
#
# Not an automated test (see test/4-workflow/2-elastic/{prep,analysis}/
# test.sh for that) -- a commented walk-through: it runs real commands,
# one case at a time, into its own output/<case>/ folder, and shows you
# the piece of output that proves what just happened. Pauses between
# sections so you can read before moving on. Safe to re-run any time -- it
# always starts by wiping its own output/.
#
# structure.fdf is a real, relaxed 2D SiC monolayer (2 atoms/cell, point
# group -6m2). calc.fdf is the real calc.fdf used for its relaxation --
# MD.VariableCell true, MD.Steps 200, kept exactly as-is on purpose (see
# Case 1 below: stb-elasticInputs forces single-point regardless).
#
# Stage 1 runs for real in every case below (cheap -- it only writes .fdf
# files, no SCF). Stage 2 is demonstrated against small, hand-built
# calc.out data (same idea as 4.1-strain's own script) since this folder
# doesn't invoke real SIESTA -- the target stiffness matrix (C11=C22=200,
# C12=40, C66=80 N/m) is a clean, hexagonal-consistent, illustrative
# number, not this material's real value (see the README's own section 4
# for the real, physically-verified result).

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

# Shared synthetic-data generator for the Stage 2 cases below -- writes a
# self-consistent calc.out (outcell + stress tensor + FreeEng) into every
# strain_* folder under $1/elastic_runs, targeting C11=C22=200 N/m,
# C12=40 N/m, C66=80 N/m (exactly hexagonal-consistent: C66=(C11-C12)/2).
write_synthetic_data() {
    local run_dir="$1"
    python3 - "$run_dir" <<'PYEOF'
import re, sys, glob
import numpy as np
from stb.core import structure_io

run_dir = sys.argv[1]
CONV = 16.021766  # eV/Ang^2 -> N/m
structure = structure_io.read_fdf(f"{run_dir}/elastic_runs/reference_structure.fdf")
pmg = structure_io.to_pymatgen(structure)
lat = pmg.lattice.matrix
Lz = np.linalg.norm(lat[2])
Area0 = np.linalg.norm(np.cross(lat[0], lat[1]))

C11, C12 = 200.0, 40.0
C66 = (C11 - C12) / 2.0
C11_raw, C12_raw, C66_raw = C11 / CONV, C12 / CONV, C66 / CONV

def voigt_to_tensor(v):
    s = np.zeros((3, 3))
    s[0,0], s[1,1], s[2,2], s[1,2], s[0,2], s[0,1] = v
    s[2,1], s[2,0], s[1,0] = s[1,2], s[0,2], s[0,1]
    return s

for folder in glob.glob(f"{run_dir}/elastic_runs/*/strain_*"):
    bname = folder.split('/')[-1]
    m = re.search(r"strain_(xx|yy|xy)_(m?)(\d+\.\d+)", bname)
    mode, sign, val = m.group(1), m.group(2), float(m.group(3))
    delta = (-val if sign == 'm' else val) / 100.0
    sigma_v = np.zeros(6)  # xx yy zz yz xz xy
    if mode == 'xx':
        sigma_v[0], sigma_v[1] = C11_raw * delta, C12_raw * delta
        dE = 0.5 * Area0 * C11_raw * delta**2
    elif mode == 'yy':
        sigma_v[0], sigma_v[1] = C12_raw * delta, C11_raw * delta
        dE = 0.5 * Area0 * C11_raw * delta**2
    else:  # xy shear -- engineering strain gamma = 2*delta
        sigma_v[5] = C66_raw * 2 * delta
        dE = 0.5 * Area0 * C66_raw * (2 * delta)**2
    sigma_v /= Lz  # dilute into SIESTA's own volumetric (eV/Ang^3) convention
    sigma_t = voigt_to_tensor(sigma_v)
    with open(f"{folder}/calc.out", "w") as f:
        f.write("outcell: Unit cell vectors (Ang):\n")
        for row in lat:
            f.write(f"        {row[0]:.6f}   {row[1]:.6f}   {row[2]:.6f}\n")
        f.write("\n")
        f.write("siesta: Stress tensor (static) (eV/Ang**3):\n")
        for row in range(3):
            f.write("siesta: " + "   ".join(f"{sigma_t[row, c]:.10f}" for c in range(3)) + "\n")
        f.write(f"siesta: FreeEng = {dE:.10f}\n")
PYEOF
}


echo "=================================================================="
echo " Welcome: from a single modulus to the full elastic tensor C_ij"
echo "=================================================================="
cat <<'EOF'
4.1-strain gives you one direction's own stiffness (a modulus, not a
tensor). This workflow gives you the real thing: the full symmetric 6x6
(3D) / 3x3 (2D, in-plane) stiffness matrix C_ij, via sigma = C . epsilon
(Hooke's law, generalized). Stage 1 (stb-elasticInputs) applies canonical
Voigt strains (xx/yy/zz/xy/xz/yz) to a relaxed structure and writes one
single-point SIESTA folder per (direction, strain value); Stage 2
(stb-elasticAnalysis) reads the finished runs back and fits C_ij.

Two independent physical routes to the SAME tensor:
  --method stress   sigma = C . epsilon, a LINEAR fit of the stress
                     tensor vs. strain (the default)
  --method energy    E = E0 + (V0/2) gamma^T C gamma, a PARABOLIC fit of
                     the total energy vs. strain -- needs combined
                     patterns (e.g. 'xx+yy') for off-diagonal constants,
                     since a single strain component's energy curvature
                     is structurally blind to cross terms

Point-group symmetry cuts the DFT cost of --method stress dramatically --
--symmetry-method controls HOW:
  basic   only pools directions related by a SINGLE point-group
          operation (misses e.g. hexagonal's C11=C22)
  full    fits the whole tensor at once inside the point group's
          symmetry-allowed subspace -- catches every crystal system's
          correct reduction, at the cost of losing basic's own
          independent tensor-symmetry cross-check (section 3 below)

See this folder's README.md for the full theory (section 1).
EOF
pause


echo "=================================================================="
echo " Case 1: Stage 1 + the single-point-SCF precedence proof, live"
echo "=================================================================="
cat <<'EOF'
calc.fdf is the SAME file used for structure.fdf's real relaxation --
MD.TypeOfRun CG, MD.VariableCell true, MD.Steps 200, all still active.
Unlike stb-strain, stb-elasticInputs needs no --relax-mode and never
trusts this block to already be single-point: it FORCES one, via a new
config_extra.fdf %include-d at the very TOP of the generated calc.fdf
(before your own directives) -- SIESTA's fdf reader is first-occurrence
-wins for duplicate labels, so this ordering guarantees the override
wins even though calc.fdf's own MD.Steps 200/MD.VariableCell true are
still sitting right there, untouched, further down. Watch the [4] report
section, then the proof itself.
EOF
mkdir -p "$OUT/case1-precedence"
cp structure.fdf calc.fdf "$OUT/case1-precedence/"
echo "\$ stb-elasticInputs -s structure.fdf -c calc.fdf -p virtual_vault --dirs all --symmetry-method full --steps 5 --no-intro"
(cd "$OUT/case1-precedence" && stb-elasticInputs -s structure.fdf -c calc.fdf -p virtual_vault \
    --dirs all --symmetry-method full --steps 5 --no-intro | tee console.log \
    | sed -n '/\[4\] SINGLE-POINT/,/\[5\] PSEUDOPOTENTIALS/p')
echo
echo "First line of the generated calc.fdf (the forced override, prepended):"
head -1 "$OUT/case1-precedence/elastic_runs/xx/strain_xx_2.00/calc.fdf"
echo
echo "...and calc.fdf's OWN MD block, still there, further down, silently overridden:"
grep -E "^MD\." "$OUT/case1-precedence/elastic_runs/xx/strain_xx_2.00/calc.fdf"
pause


echo "=================================================================="
echo " Case 2: --symmetry-method basic vs. full -- a real 3x DFT-cost gap"
echo "=================================================================="
cat <<'EOF'
Same structure.fdf/calc.fdf, same +/-2% x 5-step strain range -- only
--symmetry-method changes. This structure's point group (-6m2, hexagonal
-type) is exactly the case 'basic' cannot fully reduce: it needs all 3
in-plane directions (xx, yy, xy) to be sure it has enough independent
data, while 'full' proves 1 direction (xx alone) already determines the
whole in-plane tensor via the point group's own symmetry-allowed
subspace. Same physics, 3x the DFT cost for 'basic'.
EOF
mkdir -p "$OUT/case2-symmetry-cost"
cp structure.fdf calc.fdf "$OUT/case2-symmetry-cost/"
echo "\$ stb-elasticInputs -s structure.fdf -c calc.fdf -p virtual_vault --dirs all --symmetry-method full --steps 5 -o runs_full --no-intro"
(cd "$OUT/case2-symmetry-cost" && stb-elasticInputs -s structure.fdf -c calc.fdf -p virtual_vault \
    --dirs all --symmetry-method full --steps 5 -o runs_full --no-intro > /dev/null)
echo "\$ stb-elasticInputs -s structure.fdf -c calc.fdf -p virtual_vault --dirs all --symmetry-method basic --steps 5 -o runs_basic --no-intro"
(cd "$OUT/case2-symmetry-cost" && stb-elasticInputs -s structure.fdf -c calc.fdf -p virtual_vault \
    --dirs all --symmetry-method basic --steps 5 -o runs_basic --no-intro > /dev/null)
n_full=$(find "$OUT/case2-symmetry-cost/runs_full" -type d -name "strain_*" | wc -l)
n_basic=$(find "$OUT/case2-symmetry-cost/runs_basic" -type d -name "strain_*" | wc -l)
echo
echo "full:  $n_full SIESTA calculations (direction(s): $(find "$OUT/case2-symmetry-cost/runs_full" -maxdepth 1 -mindepth 1 -type d -printf '%f ' 2>/dev/null))"
echo "basic: $n_basic SIESTA calculations (direction(s): $(find "$OUT/case2-symmetry-cost/runs_basic" -maxdepth 1 -mindepth 1 -type d -printf '%f ' 2>/dev/null))"
echo "(the README's own worked example, section 4, shows both give the same real physics)"
pause


echo "=================================================================="
echo " Case 3: --method stress vs. energy -- the generation-time difference"
echo "=================================================================="
cat <<'EOF'
--method energy always auto-selects its own pure + combined Voigt
patterns from the point group (no direction menu, no --symmetry-method --
it always uses the general symmetry-allowed-subspace selection). On this
same structure it needs 2 patterns (xx, xy -- no combined pattern turns
out to be necessary for this specific point group's in-plane block) x 5
steps = 10 calculations, between --symmetry-method full's 5 and basic's
15 from Case 2.
EOF
mkdir -p "$OUT/case3-method"
cp structure.fdf calc.fdf "$OUT/case3-method/"
echo "\$ stb-elasticInputs -s structure.fdf -c calc.fdf -p virtual_vault --method energy --steps 5 --no-intro"
(cd "$OUT/case3-method" && stb-elasticInputs -s structure.fdf -c calc.fdf -p virtual_vault \
    --method energy --steps 5 --no-intro | sed -n '/\[3\] DEFORMATION/,/\[4\] SINGLE-POINT/p')
pause


echo "=================================================================="
echo " Stage 1 recap -- what reference_structure.fdf is, and why"
echo "=================================================================="
cat <<'EOF'
Every stb-elasticInputs run also writes ONE reference_structure.fdf at
the top of --output-dir (elastic_runs/ by default) -- a copy of the
UNDEFORMED input structure, not one from inside any strain_* folder.
Stage 2 (stb-elasticAnalysis) reads it back for 3 things:
  1. Dimensionality auto-detection (3D/2D/1D, from its vacuum padding)
  2. Symmetry/point-group detection (for --symmetry-method basic/full's
     own pooling/reconstruction, and --dimensionality auto's own logic)
  3. --method energy's reference volume/area/length (needed to convert
     the fitted energy curvature into a proper stiffness constant)
It must be undeformed: a single strained direction genuinely LOWERS the
symmetry (e.g. this -6m2 structure strained along xx alone is no longer
-6m2), so reading symmetry from inside a strain_* folder would silently
misdetect it. --dir (Stage 2, default elastic_runs, matching Stage 1's
own --output-dir default) points at the SAME directory for both the
strain_* folders and this file -- no separate flag needed in the common
case.
EOF
pause


echo "=================================================================="
echo " Case 4: Stage 2 -- a clean, consistent dataset (stb-elasticAnalysis, 4.2.2)"
echo "=================================================================="
cat <<'EOF'
This script doesn't run real SIESTA, so Stage 2 is demonstrated here
against small, hand-built calc.out data targeting a known, exactly
hexagonal-consistent stiffness matrix (C11=C22=200 N/m, C12=40 N/m,
C66=(C11-C12)/2=80 N/m) -- illustrative, not this material's real value
(see the README's section 4 for that). Both the stress-linear fit AND
the total energy are written self-consistently from the SAME target
matrix, so section [3]'s eggbox cross-check should show 0% disagreement
-- watch for it.
EOF
mkdir -p "$OUT/case4-stage2-clean"
cp structure.fdf calc.fdf "$OUT/case4-stage2-clean/"
(cd "$OUT/case4-stage2-clean" && stb-elasticInputs -s structure.fdf -c calc.fdf -p virtual_vault \
    --dirs xx yy xy --symmetry-method basic --steps 5 --no-intro > /dev/null)
write_synthetic_data "$OUT/case4-stage2-clean"
echo "\$ stb-elasticAnalysis --dir elastic_runs --symmetry-method basic --no-intro"
(cd "$OUT/case4-stage2-clean" && stb-elasticAnalysis --dir elastic_runs --symmetry-method basic \
    --no-intro | sed -n '/\[3\] NUMERICAL/,/\[4\] STABILITY/p')
pause


echo "=================================================================="
echo " Case 5: --symmetry-method full vs. basic at ANALYSIS time -- the tensor"
echo "         -symmetry check only 'basic' can run"
echo "=================================================================="
cat <<'EOF'
Same synthetic data as Case 4, analyzed both ways. 'basic' can genuinely
disagree with itself (C12 read off the xx column vs. the yy column are 2
INDEPENDENTLY measured numbers) -- that's exactly what the Tensor
symmetry (C vs C^T) line checks. 'full' constrains C12=C21 by
construction (it fits the whole tensor jointly inside the symmetry
-allowed subspace), so that same check would be structurally
uninformative -- it's skipped entirely, not just always-0.
EOF
echo "\$ stb-elasticAnalysis --dir elastic_runs --symmetry-method full --no-intro"
(cd "$OUT/case4-stage2-clean" && stb-elasticAnalysis --dir elastic_runs --symmetry-method full \
    --no-intro > full.log 2>&1)
if grep -q "Tensor symmetry" "$OUT/case4-stage2-clean/full.log"; then
    echo "UNEXPECTED: --symmetry-method full printed a Tensor symmetry line"
else
    echo "Confirmed: no 'Tensor symmetry' line for --symmetry-method full (see console.log"
    echo "from Case 4 above for 'basic''s own version of this same section)."
fi
pause


echo "=================================================================="
echo " Case 6: CLI vs. the interactive stb-suite menu (4.2.1) -- same result"
echo "=================================================================="
cat <<'EOF'
The interactive menu asks the same questions instead of flags, then
calls the exact same stb-elasticInputs underneath. Reproducing Case 1's
own generation (all directions, --symmetry-method full, 5 steps) through
stb-suite -> 4.2.1 and diffing against it proves the 2 paths are
equivalent.
EOF
mkdir -p "$OUT/case6-interactive"
cp structure.fdf calc.fdf "$OUT/case6-interactive/"
# Prompts in order: structure.fdf (no default, must type it), calc.fdf
# (blank -> default calc.fdf), pseudopotential source (2 -> virtual_vault,
# the suite's 2nd bundled bank), max strain (blank -> 2.0), steps (5),
# fitting method (blank -> stress), deformation mode (2 -> full), advanced
# settings (n -> skip, output-dir stays elastic_runs), save report (n),
# "Press Enter to continue", quit.
(cd "$OUT/case6-interactive" && \
    printf '4.2.1\nstructure.fdf\n\n2\n\n5\n\n2\nn\nn\n\n0\n' | stb-suite > menu.log 2>&1)
if diff -rq "$OUT/case1-precedence/elastic_runs/xx" "$OUT/case6-interactive/elastic_runs/xx" > /dev/null 2>&1; then
    echo "CLI (Case 1) and interactive menu (4.2.1) produced byte-identical folders."
else
    echo "UNEXPECTED: CLI and interactive menu results differ -- see"
    echo "  diff -rq '$OUT/case1-precedence/elastic_runs/xx' '$OUT/case6-interactive/elastic_runs/xx'"
fi
pause


echo "=================================================================="
echo " Workflow 4.2 complete"
echo "=================================================================="
cat <<'EOF'
Stage 1 (stb-elasticInputs, 4.2.1): single-point SIESTA folders generated,
symmetry-reduced by --symmetry-method, always forced single-point
regardless of your --calc's own MD settings, reference_structure.fdf
written alongside. Stage 2 (stb-elasticAnalysis, 4.2.2): reads real SIESTA
runs back (-d/--dir auto-finds them, default elastic_runs, matching Stage
1's own default -- no cd needed), fits the full stiffness matrix, runs 3
mandatory quality checks (eggbox cross-check, per-direction fit quality,
tensor symmetry -- the last one basic-only), reports Born mechanical
-stability criteria and a single glanceable [5] SUMMARY & FILES results
table, and can plot with gnuplot (--save-gnuplot) and/or matplotlib
(--view).

See the README's section 4 for this exact material's REAL, physically
-converged result (5 vs. 15 vs. 10 real SIESTA calculations compared
side by side) and section 6 for a step-by-step guide to running this on
your own structure.
EOF
