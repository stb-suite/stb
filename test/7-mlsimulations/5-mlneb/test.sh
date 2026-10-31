#!/bin/bash

# --- Setup ---
# Test for stb-mlneb (ML NEB, item 7.5)
#
# Needs the optional 'ml' extra (pip install stb_suite[ml] -- PyTorch + mace-torch); skipped with a clear message
# without it. Runs on the CPU with MACE-MP-0 small (downloaded once) and, for the custom-model cases, the fine-tuned
# silicon model of examples/7.2-mlmd/model/ (skipped when that folder is missing). About 10 minutes.
#
# vac_initial.fdf / vac_final.fdf: a single-atom vacancy hop in 7-atom bulk Si (the 8-atom diamond cell, a = 5.43 Ang,
# with one atom removed; one neighbour moves into the empty site, every other atom keeps its index and position). The
# two endpoints are symmetry-equivalent, so the forward and backward barriers must be equal, the reaction energy zero
# and the saddle at the middle of the path.
#
# Sections:
#   1  basic run with every output (sections, files, numbers of the band)
#   2  reproducibility
#   3  options: --n-images, --no-idpp, --freeze-substrate, --no-pre-relax-endpoints
#   4  endpoint checks: composition, lattice mismatch, identical endpoints
#   5  custom model with its card: card, trained range, foundation comparison, committee
#   6  --siesta-check with --calc (no workflow): template errors, --kgrid, pending folders, stale geometry
#   7  --siesta-check + --to-mlff with a workflow file: fake finished SIESTA runs, comparison, export
#   8  argument errors, --version, --help
#   9  interactive path (stb-suite 7.5): foundation model, then a custom model with a workflow
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$FIXTURE_DIR/test_files"
MODEL_DIR="$FIXTURE_DIR/../../../examples/7.2-mlmd/model"

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
NC='\033[0m'

PASS=0
FAIL=0

check_success() {
    if [ -s "$1" ]; then
        echo -e " ... ${GREEN}OK${NC} (file '$1' created)"
        PASS=$((PASS+1))
    else
        echo -e " ... ${RED}FAIL${NC} (file '$1' was not created)"
        FAIL=$((FAIL+1))
    fi
}

check_missing() {
    if [ -e "$1" ]; then
        echo -e " ... ${RED}FAIL${NC} ('$1' exists but should not)"
        FAIL=$((FAIL+1))
    else
        echo -e " ... ${GREEN}OK${NC} ('$1' not written, as expected)"
        PASS=$((PASS+1))
    fi
}

check_contains() {
    if grep -q -- "$1" "$2" 2>/dev/null; then
        echo -e "   -> ${GREEN}Verified:${NC} '$1' found in '$2'"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} '$1' NOT found in '$2'"
        FAIL=$((FAIL+1))
    fi
}

check_not_contains() {
    if grep -q -- "$1" "$2" 2>/dev/null; then
        echo -e "   -> ${RED}Failed:${NC} '$1' found in '$2' but should not be"
        FAIL=$((FAIL+1))
    else
        echo -e "   -> ${GREEN}Verified:${NC} '$1' not found in '$2' (as expected)"
        PASS=$((PASS+1))
    fi
}

check_exit_code() {
    if [ "$1" -eq "$2" ]; then
        echo -e "   -> ${GREEN}Verified:${NC} exit code $1 (expected $2)"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} exit code $1 (expected $2)"
        FAIL=$((FAIL+1))
    fi
}

# check_py "description" <<'PY' ... PY : a Python snippet that exits 0 when the check holds
check_py() {
    local desc="$1"
    if python3 - 2>/tmp/stb_mlneb_check_err.txt; then
        echo -e "   -> ${GREEN}Verified:${NC} $desc"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} $desc ($(tail -1 /tmp/stb_mlneb_check_err.txt))"
        FAIL=$((FAIL+1))
    fi
}

# plain text of a log (no colour codes)
plain() { sed 's/\x1b\[[0-9;]*m//g' "$1" > "$1.plain"; }


echo "--- Starting tester for stb-mlneb (item 7.5) ---"

if ! python3 -c "import mace" 2>/dev/null; then
    echo -e "${YELLOW}Skipped entirely:${NC} the optional 'ml' extra is not installed."
    echo "Install with: pip install stb_suite[ml]  (then re-run this test)"
    exit 0
fi

rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/vac_initial.fdf" "$FIXTURE_DIR/vac_final.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/../3-mlelastic/si8.fdf" "$TEST_DIR/si8.fdf"
echo "Test directory '$TEST_DIR' prepared."
pushd "$TEST_DIR" > /dev/null


# --- 1. Basic run with every output ---
echo -e "\n--- 1. Basic NEB run (MACE-MP-0 small, 7 images, every output) ---"
stb-mlneb --initial vac_initial.fdf --final vac_final.fdf --model small --n-images 7 \
    --save-images --save-report --save-data --save-gnuplot --no-intro -o vac_neb > log_neb.txt 2>&1
check_exit_code $? 0
plain log_neb.txt
for sec in "\[0\] RUN METADATA" "\[1\] STRUCTURES" "\[2\] MODEL" "\[3\] ENDPOINTS" "\[4\] STARTING PATH" \
           "\[5\] NEB BAND" "\[6\] BARRIER" "\[7\] SADDLE POINT" "\[8\] TRAINED RANGE" "\[9\] SUMMARY & FILES" \
           "\[10\] WARNINGS"; do
    check_contains "$sec" log_neb.txt.plain
done
check_contains "Device            : cpu" log_neb.txt.plain
check_contains "Atoms that move   : 1 by more than 0.1 Ang" log_neb.txt.plain
check_contains "climbing (saddle)" log_neb.txt.plain
check_contains "forward  (initial -> final)" log_neb.txt.plain
check_contains "No model card (foundation model)" log_neb.txt.plain
check_contains "Also available    : --show" log_neb.txt.plain
check_not_contains "FINE-TUNED vs FOUNDATION" log_neb.txt.plain
check_not_contains "COMMITTEE UNCERTAINTY" log_neb.txt.plain
check_not_contains "Traceback" log_neb.txt
check_success vac_neb/neb_profile.png
check_success vac_neb/neb_path.xyz
check_success vac_neb/neb_band.dat
check_success vac_neb/stb_mlneb_report.txt
check_success vac_neb/image_00.fdf
check_success vac_neb/image_06.fdf
check_missing vac_neb/image_07.fdf
check_success vac_neb/plot/neb_profile.gplot
check_success vac_neb/plot/band_0.dat
check_success vac_neb/plot/spline_0.dat
check_success vac_neb/plot/convergence_0.dat
if command -v gnuplot > /dev/null; then
    check_success vac_neb/plot/neb_profile.pdf
else
    check_contains "gnuplot is not installed" log_neb.txt.plain
fi
check_contains "STB-MLNEB REPORT" vac_neb/stb_mlneb_report.txt
check_not_contains $'\x1b\[' vac_neb/stb_mlneb_report.txt

check_py "the report's numbers: barrier 0.05-3 eV, forward = backward, reaction energy ~0, saddle at image 3 and 50 % of the path" <<'PY'
import re, sys
t = open("vac_neb/stb_mlneb_report.txt").read()
fwd = float(re.search(r"barrier ([\d.]+) eV forward", t).group(1))
bwd = float(re.search(r"([\d.]+) eV backward", t).group(1))
de = float(re.search(r"reaction energy ([-+]?[\d.]+) eV", t).group(1))
m = re.search(r"Saddle\s*: s = [\d.]+ of [\d.]+ Ang \(([\d.]+) % of the path\); climbing image (\d+)", t)
assert 0.05 < fwd < 3.0, fwd
assert abs(fwd - bwd) < 0.01, (fwd, bwd)
assert abs(de) < 0.01, de
assert abs(float(m.group(1)) - 50.0) < 2.0 and int(m.group(2)) == 3, m.groups()
PY

check_py "neb_band.dat: 7 images, E starts at 0, peaks at image 3 with the reported highest-image barrier, symmetric" <<'PY'
import re
import numpy as np
d = np.loadtxt("vac_neb/neb_band.dat")
t = open("vac_neb/stb_mlneb_report.txt").read()
img = float(re.search(r"forward  \(initial -> final\)\s*\|\s*[\d.]+\s*\|\s*([\d.]+)", t).group(1))
assert d.shape == (7, 6), d.shape
assert abs(d[0, 2]) < 1e-9 and int(np.argmax(d[:, 2])) == 3
assert abs(d[3, 2] - img) < 1e-3, (d[3, 2], img)
assert np.allclose(d[:, 2], d[::-1, 2], atol=0.01), d[:, 2]
assert np.all(np.diff(d[:, 1]) > 0), "reaction coordinate must increase"
PY

check_py "neb_path.xyz holds the RELAXED band (7 frames; its middle frame is image_03.fdf, not the interpolation)" <<'PY'
import numpy as np
from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor
from stb.core import structure_io
frames = read("vac_neb/neb_path.xyz", index=":")
assert len(frames) == 7, len(frames)
mid = AseAtomsAdaptor.get_atoms(structure_io.to_pymatgen(structure_io.read_fdf("vac_neb/image_03.fdf")))
assert np.abs(frames[3].get_positions() - mid.get_positions()).max() < 1e-4
lin = 0.5 * (frames[0].get_positions() + frames[6].get_positions())
assert np.abs(frames[3].get_positions() - lin).max() > 1e-3, "middle frame equals the linear interpolation"
PY


# --- 2. Reproducibility ---
echo -e "\n--- 2. Reproducibility: the same run gives the same band ---"
stb-mlneb --initial vac_initial.fdf --final vac_final.fdf --model small --n-images 7 --save-data \
    --no-intro -o vac_neb_again > log_again.txt 2>&1
check_exit_code $? 0
check_py "two identical runs give the same energies (1e-6 eV)" <<'PY'
import numpy as np
a, b = np.loadtxt("vac_neb/neb_band.dat"), np.loadtxt("vac_neb_again/neb_band.dat")
assert np.abs(a[:, 2] - b[:, 2]).max() < 1e-6
PY
rm -rf vac_neb_again


# --- 3. Options ---
echo -e "\n--- 3. Options: --n-images 5 --no-idpp, --freeze-substrate, --no-pre-relax-endpoints ---"
stb-mlneb --initial vac_initial.fdf --final vac_final.fdf --model small --n-images 5 --no-idpp \
    --save-images --save-data --no-intro -o vac_five > log_five.txt 2>&1
check_exit_code $? 0
plain log_five.txt
check_contains "Band              : 5 images (endpoints included), linear interpolation, atoms" log_five.txt.plain
check_not_contains "+ IDPP" log_five.txt.plain
check_success vac_five/image_04.fdf
check_missing vac_five/image_05.fdf
check_py "--n-images 5: 5 rows in neb_band.dat, saddle at image 2" <<'PY'
import numpy as np
d = np.loadtxt("vac_five/neb_band.dat")
assert d.shape[0] == 5 and int(np.argmax(d[:, 2])) == 2
PY
rm -rf vac_five

stb-mlneb --initial vac_initial.fdf --final vac_final.fdf --model small --n-images 5 \
    --freeze-substrate --no-intro -o vac_frozen > log_frozen.txt 2>&1
check_exit_code $? 0
plain log_frozen.txt
check_contains "substrate frozen below 0.3 Ang" log_frozen.txt.plain
check_contains "Frozen atoms      : 6/7" log_frozen.txt.plain
rm -rf vac_frozen

stb-mlneb --initial vac_initial.fdf --final vac_final.fdf --model small --n-images 5 \
    --no-pre-relax-endpoints --no-intro -o vac_norelax > log_norelax.txt 2>&1
check_exit_code $? 0
plain log_norelax.txt
check_contains "Endpoints relaxed : no (--no-pre-relax-endpoints)" log_norelax.txt.plain
check_contains "as given" log_norelax.txt.plain
check_contains "not relaxed with this model" log_norelax.txt.plain
rm -rf vac_norelax


# --- 4. Endpoint checks ---
echo -e "\n--- 4. Endpoint checks ---"
echo "Testing: a composition mismatch is rejected"
stb-mlneb --initial vac_initial.fdf --final si8.fdf --model small --no-intro > log_mismatch.txt 2>&1
check_exit_code $? 1
check_contains "different composition" log_mismatch.txt

echo "Testing: a lattice mismatch warns and adopts the INITIAL lattice"
sed 's/5\.43000000/5.45000000/g' vac_final.fdf > vac_final_stretched.fdf
stb-mlneb --initial vac_initial.fdf --final vac_final_stretched.fdf --model small --n-images 5 \
    --save-report --no-intro -o vac_stretched > log_lattice.txt 2>&1
check_exit_code $? 0
plain log_lattice.txt
check_contains "different lattices" log_lattice.txt.plain
check_contains "INITIAL structure's lattice is adopted" log_lattice.txt.plain
check_contains "different lattices" vac_stretched/stb_mlneb_report.txt
echo "Testing: matching lattices (section 1) do NOT warn"
check_not_contains "different lattices" log_neb.txt.plain
rm -rf vac_stretched vac_final_stretched.fdf

echo "Testing: identical endpoints -- no transition, no barrier, both said in WARNINGS"
stb-mlneb --initial vac_initial.fdf --final vac_initial.fdf --model small --n-images 5 \
    --no-intro -o vac_same > log_same.txt 2>&1
check_exit_code $? 0
plain log_same.txt
check_contains "No atom moves more than 0.1 Ang" log_same.txt.plain
check_contains "nearly identical" log_same.txt.plain
check_contains "no barrier along this band" log_same.txt.plain
rm -rf vac_same


# --- 5. Custom model with its card ---
if [ -f "$MODEL_DIR/mlff_model_compiled.model" ]; then
    M="$MODEL_DIR/mlff_model_compiled.model"
    C="$MODEL_DIR/mlff_model_card.json"
    echo -e "\n--- 5. Custom model (examples/7.2-mlmd) with its card and the foundation comparison ---"
    stb-mlneb --initial vac_initial.fdf --final vac_final.fdf --n-images 5 --model small \
        --custom-model "$M" --custom-model-card "$C" --save-data --no-intro -o vac_custom > log_custom.txt 2>&1
    check_exit_code $? 0
    plain log_custom.txt
    check_contains "Model card -- identity" log_custom.txt.plain
    check_contains "Training kinds    :" log_custom.txt.plain
    check_contains "against the model card's trained range" log_custom.txt.plain
    check_contains "FINE-TUNED vs FOUNDATION" log_custom.txt.plain
    check_contains "COMMITTEE UNCERTAINTY" log_custom.txt.plain
    check_contains "were not found next to it" log_custom.txt.plain
    check_success vac_custom/neb_band_foundation.dat
    check_success vac_custom/neb_path_foundation.xyz
    check_py "both models report a barrier; the comparison row matches the summary" <<'PY'
import re
t = open("log_custom.txt.plain").read()
rows = re.findall(r"^(fine-tuned|MACE-MP-0 small)\s*: barrier ([\d.]+) eV forward", t, re.M)
assert len(rows) == 2, rows
cmp_ = re.search(r"barrier, spline \(eV\)\s*\|\s*([\d.]+)\s*\|\s*([\d.]+)", t)
assert abs(float(cmp_.group(1)) - float(rows[0][1])) < 1e-3 and abs(float(cmp_.group(2)) - float(rows[1][1])) < 1e-3
PY
    rm -rf vac_custom

    echo "Testing: --committee-models (the same model twice: zero spread, equal member barriers)"
    stb-mlneb --initial vac_initial.fdf --final vac_final.fdf --n-images 5 --custom-model "$M" \
        --custom-model-card "$C" --skip-foundation-comparison --committee-models "$M" "$M" \
        --no-intro -o vac_committee > log_committee.txt 2>&1
    check_exit_code $? 0
    plain log_committee.txt
    check_contains "Committee         : 2 model(s) from --committee-models" log_committee.txt.plain
    check_contains "Calibration       : on the card's test set" log_committee.txt.plain
    check_not_contains "FINE-TUNED vs FOUNDATION" log_committee.txt.plain
    check_py "identical members: force spread 0.0 on every image, the same barrier twice" <<'PY'
import re
t = open("log_committee.txt.plain").read()
sec = t[t.index("COMMITTEE UNCERTAINTY"):t.index("SUMMARY & FILES")]
spreads = re.findall(r"^\d+\s*\|\s*[-+][\d.]+\s*\|\s*([\d.]+)\s*\|\s*([\d.]+)\s*\|", sec, re.M)
assert len(spreads) == 5 and all(float(e) == 0 and float(f) == 0 for e, f in spreads), spreads
b = re.search(r"Barrier per member: ([\d.]+), ([\d.]+) eV", sec)
assert b and b.group(1) == b.group(2), sec
PY
    rm -rf vac_committee

    echo "Testing: --no-committee"
    stb-mlneb --initial vac_initial.fdf --final vac_final.fdf --n-images 5 --custom-model "$M" \
        --custom-model-card "$C" --skip-foundation-comparison --no-committee --no-intro -o vac_nocom > log_nocom.txt 2>&1
    check_exit_code $? 0
    plain log_nocom.txt
    check_contains "Not evaluated     : --no-committee" log_nocom.txt.plain
    rm -rf vac_nocom
else
    echo -e "${YELLOW}Skipped section 5:${NC} examples/7.2-mlmd/model not found"
fi


# --- 6. --siesta-check with --calc (no workflow) ---
echo -e "\n--- 6. --siesta-check with --calc and --pseudo-dir ---"
cat > calc_noinclude.fdf <<'EOF'
SystemLabel siesta
PAO.BasisSize DZP
EOF
stb-mlneb --initial vac_initial.fdf --final vac_final.fdf --model small --n-images 5 --siesta-check \
    --calc calc_noinclude.fdf --no-intro -o vac_calc > log_noinclude.txt 2>&1
check_exit_code $? 1
check_contains "does not %include a structure file" log_noinclude.txt

stb-mlneb --initial vac_initial.fdf --final vac_final.fdf --model small --n-images 5 --siesta-check \
    --no-intro -o vac_calc > log_nocalc.txt 2>&1
check_exit_code $? 1
check_contains "needs --calc (or --mlff-workflow)" log_nocalc.txt
echo "Testing: a bad SIESTA set-up stops before the NEB, not after it"
check_not_contains "NEB BAND" log_noinclude.txt
check_not_contains "NEB BAND" log_nocalc.txt
mkdir -p empty_pseudo
stb-mlneb --initial vac_initial.fdf --final vac_final.fdf --model small --n-images 5 --siesta-check \
    --calc "$FIXTURE_DIR/../../../examples/7.5-mlneb/calc.fdf" --pseudo-dir empty_pseudo --no-intro -o vac_calc \
    > log_nopseudo.txt 2>&1
check_exit_code $? 1
check_contains "no pseudopotential for Si" log_nopseudo.txt
check_not_contains "NEB BAND" log_nopseudo.txt

cat > calc.fdf <<'EOF'
SystemLabel siesta
%include geometry.fdf
PAO.BasisSize DZP
Mesh.Cutoff 300 Ry
kgrid.MonkhorstPack [6 6 6]
EOF
stb-mlneb --initial vac_initial.fdf --final vac_final.fdf --model small --n-images 5 --siesta-check \
    --calc calc.fdf --pseudo-dir dojo --kgrid 2 2 2 --no-intro -o vac_calc > log_calc1.txt 2>&1
check_exit_code $? 0
plain log_calc1.txt
check_contains "SIESTA CHECK" log_calc1.txt.plain
check_contains "Template          : --calc calc.fdf" log_calc1.txt.plain
check_contains "Written           : 3 single-point folder(s)" log_calc1.txt.plain
check_contains "still need SIESTA" log_calc1.txt.plain
for i in 00 02 04; do
    check_success vac_calc/siesta_check/image_$i/geometry.fdf
done
check_missing vac_calc/siesta_check/image_01
check_contains "%include config_extra.fdf" vac_calc/siesta_check/image_02/calc.fdf
check_contains "kgrid.MonkhorstPack \[2 2 2\]" vac_calc/siesta_check/image_02/config_extra.fdf
check_success vac_calc/siesta_check/image_02/Si.psml
check_success vac_calc/siesta_check/image_02/stb_mlneb_image.xyz

echo "Testing: a second run without SIESTA keeps the folders pending and writes nothing new"
stb-mlneb --initial vac_initial.fdf --final vac_final.fdf --model small --n-images 5 --siesta-check \
    --calc calc.fdf --pseudo-dir dojo --kgrid 2 2 2 --no-intro -o vac_calc > log_calc2.txt 2>&1
check_exit_code $? 0
plain log_calc2.txt
check_not_contains "Written           :" log_calc2.txt.plain
check_contains "pending" log_calc2.txt.plain

echo "Testing: a folder written for another geometry is left alone and reported"
python3 - <<'PY'
from ase.io import read, write
a = read("vac_calc/siesta_check/image_00/stb_mlneb_image.xyz")
a.positions[0] += 0.3
write("vac_calc/siesta_check/image_00/stb_mlneb_image.xyz", a, format="extxyz")
PY
stb-mlneb --initial vac_initial.fdf --final vac_final.fdf --model small --n-images 5 --siesta-check \
    --calc calc.fdf --pseudo-dir dojo --kgrid 2 2 2 --no-intro -o vac_calc > log_calc3.txt 2>&1
check_exit_code $? 0
plain log_calc3.txt
check_contains "was written for another geometry" log_calc3.txt.plain
check_contains "stale" log_calc3.txt.plain
rm -rf vac_calc calc.fdf calc_noinclude.fdf


# --- 7. --siesta-check + --to-mlff with a workflow file ---
echo -e "\n--- 7. --siesta-check and --to-mlff with a small workflow file ---"
mkdir -p wf
cp si8.fdf wf/structure.fdf
cat > wf/calc.fdf <<'EOF'
SystemLabel siesta
%include structure.fdf
PAO.BasisSize DZP
Mesh.Cutoff 300 Ry
kgrid.MonkhorstPack [6 6 6]
EOF
cat > wf/mlff_workflow.json <<'EOF'
{"format": "stb-mlff-workflow", "version": 1, "reference": {"file": "structure.fdf", "calc": "calc.fdf",
 "pseudo_dir": "dojo"}, "configs": [], "rounds": [{"round": 1, "dir": "round_01"}], "history": []}
EOF
stb-mlneb --initial vac_initial.fdf --final vac_final.fdf --n-images 5 --model small \
    --mlff-workflow wf/mlff_workflow.json --siesta-check --no-intro -o vac_check > log_check1.txt 2>&1
check_exit_code $? 0
plain log_check1.txt
check_contains "Template          : the workflow's reference" log_check1.txt.plain
check_contains "Written           : 3 single-point folder(s)" log_check1.txt.plain
check_contains "NumberOfAtoms 7" vac_check/siesta_check/image_02/config_extra.fdf
check_contains "MD.Steps              0" vac_check/siesta_check/image_02/config_extra.fdf
echo "Testing: the k-grid keeps the reference's k x a (6 in 5.48 Ang -> 6 in 5.43 Ang)"
check_contains "kgrid.MonkhorstPack \[6 6 6\]" vac_check/siesta_check/image_02/config_extra.fdf

# fake finished SIESTA runs. First the saddle 10 eV above the endpoints: far from the model -> "not reliable".
fake_siesta() {   # $1 = saddle energy above the endpoints (eV)
    for i in 00 02 04; do
        python3 - "$i" "$1" <<'PY'
import sys
i, de = sys.argv[1], float(sys.argv[2])
e = -100.0 + (de if i == "02" else 0.0)
d = f"vac_check/siesta_check/image_{i}"
open(f"{d}/siesta.out", "w").write(f"siesta: FreeEng =   {e:.6f}\nEnd of run\n")
open(f"{d}/siesta.FA", "w").write("7\n" + "".join(f"{k + 1} 0.01 0.0 0.0\n" for k in range(7)))
PY
    done
}
fake_siesta 10.0
stb-mlneb --initial vac_initial.fdf --final vac_final.fdf --n-images 5 --model small \
    --mlff-workflow wf/mlff_workflow.json --siesta-check --no-intro -o vac_check > log_check2.txt 2>&1
check_exit_code $? 0
plain log_check2.txt
check_contains "SIESTA 10.0000" log_check2.txt.plain
check_contains "model - SIESTA" log_check2.txt.plain
check_contains "not reliable for this transition" log_check2.txt.plain
check_contains "SIESTA check      : barrier on the model's path 10.0000 eV" log_check2.txt.plain

echo "Testing: SIESTA equal to the model's barrier -> no SIESTA warning"
MB=$(python3 -c "
import re; t=open('log_check1.txt.plain').read()
print(re.search(r'forward  \(initial -> final\)\s*\|\s*[\d.]+\s*\|\s*([\d.]+)', t).group(1))")
fake_siesta "$MB"
stb-mlneb --initial vac_initial.fdf --final vac_final.fdf --n-images 5 --model small \
    --mlff-workflow wf/mlff_workflow.json --siesta-check --to-mlff --to-mlff-rattle 1 --save-report \
    --no-intro -o vac_check > log_check3.txt 2>&1
check_exit_code $? 0
plain log_check3.txt
check_contains "([+-]0 % model vs SIESTA)" log_check3.txt.plain
check_not_contains "SIESTA check: " log_check3.txt.plain
check_contains "EXPORT TO THE MLFF WORKFLOW" log_check3.txt.plain
check_contains "latest round 1; these configs train round 2" log_check3.txt.plain
check_contains "3 already labelled by --siesta-check" log_check3.txt.plain
check_contains "7 need SIESTA" log_check3.txt.plain
check_success wf/configs_round_02/neb_vac_check/image_02/siesta.out
check_success wf/configs_round_02/neb_vac_check/image_01/calc.fdf
check_missing wf/configs_round_02/neb_vac_check/image_01/siesta.out
check_success wf/configs_round_02/neb_vac_check/image_01_rattle_1/structure.fdf
check_success wf/mlff_workflow.json.bak
check_py "the workflow records 10 training configs (5 images + 5 rattled), kinds neb / neb-rattle, round 1" <<'PY'
import json
w = json.load(open("wf/mlff_workflow.json"))
c = [x for x in w["configs"] if x.get("stage") == "neb"]
assert len(c) == 10 and all(x["set"] == "train" and x["round"] == 1 and x["n_atoms"] == 7 for x in c)
assert sorted({x["kind"] for x in c}) == ["neb", "neb-rattle"]
assert all(x["folder"].startswith("configs_round_02/neb_vac_check/") for x in c), c[0]["folder"]
assert w["history"][-1]["tool"] == "stb-mlneb" and w["history"][-1]["stage"] == "neb"
PY
check_py "a rattled copy differs from its image by about --to-mlff-stdev (0.05 Ang)" <<'PY'
import numpy as np
from pymatgen.io.ase import AseAtomsAdaptor
from stb.core import structure_io
def atoms(p):
    return AseAtomsAdaptor.get_atoms(structure_io.to_pymatgen(structure_io.read_fdf(p)))
a = atoms("wf/configs_round_02/neb_vac_check/image_01/structure.fdf")
b = atoms("wf/configs_round_02/neb_vac_check/image_01_rattle_1/structure.fdf")
d = b.get_positions() - a.get_positions()
d -= np.round(d / a.cell.lengths()) * a.cell.lengths()
rms = float(np.sqrt((d ** 2).mean()))
assert 0.02 < rms < 0.1, rms
PY

echo "Testing: a second --to-mlff does not duplicate the configs"
stb-mlneb --initial vac_initial.fdf --final vac_final.fdf --n-images 5 --model small \
    --mlff-workflow wf/mlff_workflow.json --to-mlff --to-mlff-rattle 1 --no-intro -o vac_check > log_check4.txt 2>&1
check_exit_code $? 0
plain log_check4.txt
check_contains "10 already in the workflow" log_check4.txt.plain
check_py "still 10 NEB configs in the workflow" <<'PY'
import json
w = json.load(open("wf/mlff_workflow.json"))
assert sum(1 for x in w["configs"] if x.get("stage") == "neb") == 10
PY
rm -rf vac_check

echo "Testing: --to-mlff with -o . names the export folder after the real output folder"
mkdir -p dot_run
(cd dot_run && stb-mlneb --initial ../vac_initial.fdf --final ../vac_final.fdf --n-images 5 --model small \
    --mlff-workflow ../wf/mlff_workflow.json --to-mlff --no-intro -o . > ../log_dot.txt 2>&1)
check_exit_code $? 0
check_success wf/configs_round_02/neb_dot_run/image_02/calc.fdf
check_missing "wf/configs_round_02/neb_./image_02"
rm -rf dot_run


# --- 8. Argument errors, --version, --help ---
echo -e "\n--- 8. Argument errors, --version, --help ---"
stb-mlneb --initial vac_initial.fdf --final vac_final.fdf --n-images 2 --no-intro > log_nimages.txt 2>&1
check_exit_code $? 2
check_contains "at least 3" log_nimages.txt
stb-mlneb --initial vac_initial.fdf --final vac_final.fdf --to-mlff --no-intro > log_tomlff.txt 2>&1
check_exit_code $? 2
check_contains "needs --mlff-workflow" log_tomlff.txt
stb-mlneb --initial vac_initial.fdf --final vac_final.fdf --mlff-workflow nope.json --no-intro > log_nowf.txt 2>&1
check_exit_code $? 2
check_contains "not found" log_nowf.txt
stb-mlneb --final vac_final.fdf --no-intro > log_missing_initial_arg.txt 2>&1
check_exit_code $? 2
stb-mlneb --initial nope.fdf --final vac_final.fdf --model small --no-intro > log_missing_initial.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_initial.txt
stb-mlneb --initial vac_initial.fdf --final vac_final.fdf --custom-model does_not_exist.model \
    --custom-model-card does_not_exist.json --no-intro > log_custommodel.txt 2>&1
check_exit_code $? 2
check_contains "not found" log_custommodel.txt
stb-mlneb --initial vac_initial.fdf --final vac_final.fdf --custom-model does_not_exist.model \
    --no-intro > log_custommodel_nocard.txt 2>&1
check_exit_code $? 2
check_contains "custom-model-card" log_custommodel_nocard.txt
stb-mlneb --version > log_version.txt 2>&1
check_contains "stb-mlneb" log_version.txt
stb-mlneb --help > log_help.txt 2>&1
for opt in --initial --final --n-images --custom-model-card --skip-foundation-comparison --mlff-workflow \
           --no-committee --committee-models --siesta-check --siesta-check-images --calc --pseudo-dir --kgrid \
           --to-mlff --to-mlff-rattle --freeze-substrate --no-idpp --device --show --save-gnuplot --save-data; do
    check_contains "$opt" log_help.txt
done


# --- 9. Interactive path (stb-suite, shortcut 7.5) ---
echo -e "\n--- 9. Interactive path via stb-suite (shortcut 7.5) ---"
echo "Testing: foundation model, gnuplot files asked for"
# initial, final, custom (none), size, device, images, freeze, output, show, gnuplot, data, images, report, Enter, exit
printf '7.5\nvac_initial.fdf\nvac_final.fdf\n\nsmall\n\n5\nn\ninteractive_out\nn\ny\nn\nn\nn\n\n0\n' \
    | stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success interactive_out/neb_profile.png
check_success interactive_out/plot/neb_profile.gplot
check_success interactive_out/neb_path.xyz
check_missing interactive_out/stb_mlneb_report.txt

if [ -f "$MODEL_DIR/mlff_model_compiled.model" ]; then
    echo "Testing: custom model (card offered as default), device asked, workflow + --siesta-check from the menu"
    # initial, final, custom, card (default), size, device, images, freeze, compare (n), workflow, siesta-check (y),
    # all images (n), to-mlff (n), output, show, gnuplot, data, images, report (y), Enter, exit
    printf '7.5\nvac_initial.fdf\nvac_final.fdf\n%s\n\nsmall\ncpu\n5\nn\nn\nwf/mlff_workflow.json\ny\nn\nn\ninteractive_custom\nn\nn\nn\nn\ny\n\n0\n' \
        "$MODEL_DIR/mlff_model_compiled.model" | stb-suite > log_menu_custom.txt 2>&1
    check_exit_code $? 0
    check_success interactive_custom/stb_mlneb_report.txt
    check_contains "Device            : cpu" interactive_custom/stb_mlneb_report.txt
    check_contains "Model card -- identity" interactive_custom/stb_mlneb_report.txt
    check_contains "SIESTA CHECK" interactive_custom/stb_mlneb_report.txt
    check_not_contains "FINE-TUNED vs FOUNDATION" interactive_custom/stb_mlneb_report.txt
    check_success interactive_custom/siesta_check/image_02/calc.fdf
fi
rm -rf wf


popd > /dev/null

# --- Summary ---
echo -e "\n--- Tests Complete ---"
echo -e "${GREEN}Passed: $PASS${NC}   ${RED}Failed: $FAIL${NC}"

read -p "Remove the '$TEST_DIR' directory and all test files? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "Cleaning up test files..."
    rm -rf "$TEST_DIR"
    echo -e "${GREEN}Cleanup complete.${NC}"
else
    echo "Test files were kept in '$TEST_DIR/' for inspection."
fi

[ "$FAIL" -eq 0 ]
