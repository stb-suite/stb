#!/bin/bash

# --- Setup ---
# Smoke test for stb-mlmd (ML Molecular Dynamics, item 7.2)
#
# Needs the optional 'ml' extra (pip install stb_suite[ml]). The whole file
# is skipped with a clear message if `mace` isn't importable, same gating
# pattern as test/1-inputs/7-mlrelax/test.sh. si_bulk.fdf is the same
# 64-atom Ge:Si substitutional-defect fixture used there (bulk, 3D
# periodic -- suitable for NVE/NVT/NPT alike).
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$FIXTURE_DIR/test_files"

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
        echo -e "   -> ${RED}Failed:${NC} '$1' found in '$2'"
        FAIL=$((FAIL+1))
    else
        echo -e "   -> ${GREEN}Verified:${NC} '$1' not in '$2'"
        PASS=$((PASS+1))
    fi
}

check_python() {   # $1 = description, $2 = python code exiting 0 on success
    python3 -c "$2"
    if [ $? -eq 0 ]; then
        echo -e "   -> ${GREEN}Verified:${NC} $1"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} $1"
        FAIL=$((FAIL+1))
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


echo "--- Starting tester for stb-mlmd (item 7.2) ---"

if ! python3 -c "import mace" 2>/dev/null; then
    echo -e "${YELLOW}Skipped entirely:${NC} the optional 'ml' extra is not installed."
    echo "Install with: pip install stb_suite[ml]  (then re-run this test)"
    exit 0
fi

rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/si_bulk.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 1. NVT (default ensemble) ---
echo -e "\n--- Testing NVT (default), 100 steps, stride 20 ---"
stb-mlmd --file si_bulk.fdf --ensemble nvt --temperature 300 --n-steps 100 --stride 20 \
    -o nvt_traj.xsf --output-structure nvt_final.fdf --png --save-gnuplot --no-intro > log_nvt.txt 2>&1
check_exit_code $? 0
check_contains "Ensemble          : NVT" log_nvt.txt
check_contains "Collected 5 frame(s)" log_nvt.txt
check_contains "Temperature       : mean" log_nvt.txt
for sec in "\[0\] RUN METADATA" "\[1\] STRUCTURE" "\[2\] MODEL" "\[3\] MD SETUP" "\[4\] EQUILIBRATION" \
           "\[5\] PRODUCTION RUN" "\[6\] THERMODYNAMICS" "\[7\] ENSEMBLE CHECKS" \
           "\[8\] STRUCTURE ALONG THE RUN" "\[9\] SUMMARY & FILES" "\[10\] WARNINGS"; do
    check_contains "$sec" log_nvt.txt
done
check_contains "Speed             :" log_nvt.txt
check_contains "T fluctuation     :" log_nvt.txt
check_contains "g(r)              : first peak" log_nvt.txt
check_contains "Lindemann ratio" log_nvt.txt
check_contains "Model card        : none" log_nvt.txt
check_contains "coordination 4.00" log_nvt.txt
check_contains "Heat capacity     :" log_nvt.txt
check_contains "density" log_nvt.txt
echo "Testing: MACE/torch/ASE messages are counted, not printed in the report"
check_not_contains "UserWarning" log_nvt.txt
check_not_contains "cuequivariance" log_nvt.txt
check_not_contains "WARNING:root" log_nvt.txt
check_contains "Library messages  :" log_nvt.txt
echo "Testing: the WARNINGS section is the last one"
check_python "nothing but the library-messages line follows [10] WARNINGS" "
import re, sys
t = re.sub(r'\x1b\[[0-9;]*m', '', open('log_nvt.txt').read())
tail = t.split('[10] WARNINGS', 1)[1]
sys.exit(0 if not re.search(r'^\[\d+\] ', tail, re.M) else 1)
"
check_success nvt_traj.xsf
check_success nvt_final.fdf
check_success si_bulk_md_dashboard.png
check_success si_bulk_md_thermo.dat
check_success si_bulk_md_rdf.dat
check_success si_bulk_md_msd.dat
check_success si_bulk_md.gplot
check_python "thermo .dat: 12 columns, one row per frame (5), temperatures as in the report" "
import sys, numpy as np
a = np.loadtxt('si_bulk_md_thermo.dat')
ok = a.shape == (5, 12) and np.allclose(a[:, 0], [20, 40, 60, 80, 100]) and (a[:, 1] > 0).all() and (a[:, 6] > 15).all()
r = np.loadtxt('si_bulk_md_rdf.dat'); m = np.loadtxt('si_bulk_md_msd.dat')
ok = ok and r.shape[1] == 3 and abs(r[np.argmax(r[:, 1]), 0] - 2.37) < 0.15 and m.shape[1] >= 2
sys.exit(0 if ok else 1)
"
if command -v gnuplot &> /dev/null; then
    gnuplot si_bulk_md.gplot > log_gnuplot.txt 2>&1
    check_exit_code $? 0
    check_success si_bulk_md.pdf
fi

echo "Testing: plots and data files are opt-in"
rm -f si_bulk_md_dashboard.png si_bulk_md_thermo.dat si_bulk_md.gplot
stb-mlmd --file si_bulk.fdf --n-steps 20 --stride 10 -o plain.xsf --output-structure plain.fdf --no-intro > log_plain.txt 2>&1
check_exit_code $? 0
check_python "no PNG, .dat or .gplot written without --png/--save-data/--save-gnuplot" "
import os, sys
sys.exit(0 if not any(os.path.exists(f) for f in ('si_bulk_md_dashboard.png', 'si_bulk_md_thermo.dat', 'si_bulk_md.gplot')) else 1)
"
check_contains "Also available    : --png" log_plain.txt


# --- 1b. --equilibration-steps discards the initial transient ---
echo -e "\n--- Testing --equilibration-steps ---"
stb-mlmd --file si_bulk.fdf --ensemble nvt --equilibration-steps 20 --n-steps 40 --stride 10 \
    --no-intro > log_equil.txt 2>&1
check_exit_code $? 0
check_contains "EQUILIBRATION (discarded)" log_equil.txt
check_contains "Equilibration done" log_equil.txt
check_contains "20 equilibration" log_equil.txt


# --- 1c. --supercell repeats the input cell before the run ---
echo -e "\n--- Testing --supercell 1 1 2 (64 -> 128 atoms) ---"
stb-mlmd --file si_bulk.fdf --supercell 1 1 2 --n-steps 20 --stride 10 --seed 1 \
    -o sc_traj.xsf --output-structure sc_final.fdf --no-intro > log_supercell.txt 2>&1
check_exit_code $? 0
check_contains "Supercell         : 1 x 1 x 2" log_supercell.txt
check_contains "Supercell 1 x 1 x 2: 64 -> 128 atoms" log_supercell.txt
check_contains "Degrees of freedom: 381" log_supercell.txt
check_python "md_final.fdf and every trajectory frame hold 128 atoms, c doubled, a and b unchanged" "
import sys, numpy as np
from ase.io import read
from stb.core import structure_io as s
l0 = np.linalg.norm(s.read_fdf('si_bulk.fdf').lattice, axis=1)
f = s.read_fdf('sc_final.fdf'); l1 = np.linalg.norm(f.lattice, axis=1)
frames = read('sc_traj.xsf', index=':')
ok = len(f.atoms) == 128 and np.allclose(l1, l0 * [1, 1, 2], atol=1e-4) and all(len(a) == 128 for a in frames)
sys.exit(0 if ok else 1)
"
echo "Testing: the supercell's atoms are the input cell's, repeated (same order as ase.Atoms.repeat)"
check_python "first 64 atoms of frame 0 sit near the input positions; species follow the input order twice" "
import sys, numpy as np
from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor
from stb.core import structure_io as s
a0 = AseAtomsAdaptor.get_atoms(s.to_pymatgen(s.read_fdf('si_bulk.fdf')))
f0 = read('sc_traj.xsf', index=0)
d = f0.positions[:64] - a0.positions
d -= np.round(d @ np.linalg.inv(a0.cell)) @ a0.cell
ok = f0.get_chemical_symbols() == a0.get_chemical_symbols() * 2 and np.abs(d).max() < 0.5
sys.exit(0 if ok else 1)
"
echo "Testing: --supercell rejects values < 1"
stb-mlmd --file si_bulk.fdf --supercell 2 0 1 --no-intro > log_supercell_bad.txt 2>&1
check_exit_code $? 2
check_contains "three integers >= 1" log_supercell_bad.txt
echo "Testing: repeating a vacuum-padded axis is warned about"
stb-slab -f si_bulk.fdf --hkl 1 0 0 -o slab_sc.fdf --no-intro > /dev/null 2>&1
stb-mlmd --file slab_sc.fdf --supercell 1 1 2 --n-steps 10 --stride 10 -o slab_sc.xsf \
    --output-structure slab_sc_final.fdf --no-intro > log_supercell_slab.txt 2>&1
check_exit_code $? 0
check_contains "but it is vacuum-padded" log_supercell_slab.txt


# --- 2. NVE (checks total-energy conservation reporting) ---
echo -e "\n--- Testing NVE (energy conservation report) ---"
stb-mlmd --file si_bulk.fdf --ensemble nve --n-steps 200 --stride 20 --seed 1 \
    -o nve_traj.xsf --output-structure nve_final.fdf --no-intro > log_nve.txt 2>&1
check_exit_code $? 0
check_contains "Ensemble          : NVE" log_nve.txt
check_contains "Energy drift      :" log_nve.txt
check_contains "E_tot fluctuation :" log_nve.txt
check_success nve_traj.xsf

echo "Testing: NVE drift is reported in meV/atom/ps and is small (well-conserved integrator)"
check_python "NVE drift |d| < 1 meV/atom/ps (the old relative drift divided by the absolute E_tot and could never warn)" "
import re, sys
text = re.sub(r'\x1b\[[0-9;]*m', '', open('log_nve.txt').read())
m = re.search(r'Energy drift      : ([-+\d.eE]+) meV/atom/ps', text)
sys.exit(0 if m and abs(float(m.group(1))) < 1.0 else 1)
"
check_python "slope_per_ps: a linear E(t) gives its slope per ps" "
import sys, numpy as np
from stb.mlmd import slope_per_ps, block_se
t = np.arange(0, 1000, 10.0)
sys.exit(0 if abs(slope_per_ps(t, 2.5 * t / 1000 + 7) - 2.5) < 1e-9 and block_se(np.ones(50)) == 0 else 1)
"
echo "Testing: the model card warns once per KIND of message (a growing cell once gave 147 'volume per atom' lines)"
check_python "same message with new numbers -> one collected entry with a repeat count; a different kind -> its own entry" "
import sys
from stb.core import model_card as m
got = []
m.collect_warnings(got)
for v in (22.889, 22.928, 23.5):
    m._warn_once(f'volume per atom (Ang^3) is {v:.3f}, outside the trained range 18.722 .. 22.425')
m._warn_once('shortest interatomic distance 1.731 Ang is below the shortest trained one (1.852 Ang)')
m.collect_warnings(None)
sys.exit(0 if len(got) == 2 and 'first of 3 such warnings' in got[0] and '22.889' in got[0] else 1)
"


# --- 3. NPT on bulk (allowed) ---
echo -e "\n--- Testing NPT on a bulk structure (allowed) ---"
stb-mlmd --file si_bulk.fdf --ensemble npt --temperature 300 --pressure 0.0 --n-steps 40 --stride 10 \
    -o npt_traj.xsf --output-structure npt_final.fdf --no-intro > log_npt.txt 2>&1
check_exit_code $? 0
check_contains "Ensemble          : NPT" log_npt.txt
check_contains "Lattice (mean)    :" log_npt.txt
check_contains "Barostat          :" log_npt.txt
check_contains "placeholder (water)" log_npt.txt
check_success npt_traj.xsf


# --- 4. NPT on a vacuum-padded slab (must fail) ---
echo -e "\n--- Testing NPT on a 2D slab (expects failure) ---"
stb-slab -f si_bulk.fdf --hkl 1 0 0 -o slab_for_test.fdf --no-intro > /dev/null 2>&1
stb-mlmd --file slab_for_test.fdf --ensemble npt --n-steps 10 --no-intro > log_npt_slab.txt 2>&1
check_exit_code $? 1
check_contains "requires a bulk" log_npt_slab.txt
rm -f slab_for_test.fdf


# --- 5. --custom-model with a nonexistent file (expects failure) ---
echo -e "\n--- Testing --custom-model with a nonexistent file (expects failure) ---"
stb-mlmd --file si_bulk.fdf --custom-model does_not_exist.model --custom-model-card does_not_exist.json --n-steps 10 --no-intro > log_custommodel.txt 2>&1
check_exit_code $? 2
check_contains "not found" log_custommodel.txt

echo "Testing: a custom model without its model card (--custom-model-card) is rejected"
stb-mlmd --file si_bulk.fdf --custom-model does_not_exist.model --n-steps 10 --no-intro > log_custommodel_nocard.txt 2>&1
check_exit_code $? 2
check_contains "custom-model-card" log_custommodel_nocard.txt

echo "Testing: a card that belongs to another model is refused (SHA-256)"
python3 - <<'EOF'
import hashlib, json
open("dummy_a.model", "wb").write(b"model A")
json.dump({"format": "stb-mlff-model-card", "version": 1, "elements": ["Si"],
           "model": {"files": [{"file": "dummy_b.model", "sha256": hashlib.sha256(b"model B").hexdigest()}]}},
          open("dummy_card_b.json", "w"))
EOF
stb-mlmd --file si_bulk.fdf --custom-model dummy_a.model --custom-model-card dummy_card_b.json --n-steps 10 --no-intro > log_custommodel_sha.txt 2>&1
check_exit_code $? 2
check_contains "SHA-256" log_custommodel_sha.txt

echo "Testing: the report prints the model card and turns the card's own alarms into warnings"
check_python "print_model_card: identity, theory, trained range, accuracy; no validation and another member's note handled" "
import hashlib, io, sys, contextlib
from stb.mlmd import print_model_card
open('card_model.model', 'wb').write(b'weights')
sha = hashlib.sha256(b'weights').hexdigest()
card = {'format': 'stb-mlff-model-card', 'version': 1, 'elements': ['Si'],
        'model': {'name': 'round_02_seed_125', 'seed': 125, 'member': 3, 'files': [{'file': 'card_model.model', 'sha256': sha, 'kind': 'compiled'}],
                  'architecture': {'r_max': 6.0, 'num_interactions': 2, 'n_params': 1000}, 'foundation': {'name': 'mace-mp-0', 'size': 'medium'}},
        'committee': {'size': 3, 'round': 2},
        'level_of_theory': {'siesta_version': '5.4.2', 'settings': {'XC.authors': 'PBE', 'PAO.BasisSize': 'DZP'}, 'dftd3': {'on': False}},
        'training_data': {'n_configs': 57, 'kinds': {'rattle': 39}, 'labels': {'energy': True, 'forces': True}},
        'domain': {'elements': ['Si'], 'dimensionality': '3D', 'measure_name': 'volume per atom (Ang^3)', 'measure': {'min': 18.7, 'max': 22.4},
                   'd_min_ang': {'min': 1.85, 'max': 2.35}},
        'metrics': {'train': {'n': 51, 'F_rmse_meV': 45.8}, 'out_of_sample_committee': {'n': 16, 'F_rmse_meV': 56.0}},
        'checks': {'stability': {'status': 'OK', 'n_atoms': 64}},
        'provenance': {'training_warnings': ['Training: seed 124 kept epoch 98', 'Training: seed 125 overfits']}}
w = []
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    print_model_card(card, 'card.json', 'card_model.model', None, w.append)
out = buf.getvalue()
ok = all(k in out for k in ('Model card -- identity', '[verified]', 'level of theory', 'Trained range', '18.700 .. 22.400',
                            'Errors vs SIESTA', 'out-of-sample (committee)', 'member 3 of 3'))
ok = ok and any('no test-set validation' in x for x in w) and any('seed 125 overfits' in x for x in w) \
        and not any('seed 124' in x for x in w)
sys.exit(0 if ok else 1)
"

# --- 5b. A fine-tuned model with its card, end to end ---
echo -e "\n--- Testing a fine-tuned model with its model card, end to end (tiny O2 model from the 7.1 fixtures) ---"
MLFF_FIX="$FIXTURE_DIR/../../7-mlsimulations/1-mlff/analysis"
if command -v mace_run_train &> /dev/null && ls -d "$MLFF_FIX"/mlff_config_* > /dev/null 2>&1; then
    mkdir -p ft
    cp -r "$MLFF_FIX"/mlff_config_* ft/
    (cd ft && stb-mlffAnalysis --path "mlff_config_*" --committee 1 -o model --epochs 2 --batch-size 4 --device cpu \
        --foundation-model small --no-stability --name ft --no-intro > ../log_ft_train.txt 2>&1)
    check_exit_code $? 0
    check_success ft/model/ft_compiled.model
    check_success ft/model/mlff_model_card.json
    CM="--custom-model ft/model/ft_compiled.model --custom-model-card ft/model/mlff_model_card.json"

    echo "Testing: the report shows the card section by section and checks the run against the trained range"
    rm -f stb_mlmd_report.txt
    stb-mlmd --file ft/mlff_config_001/calc.fdf $CM --n-steps 60 --stride 10 --seed 1 --save-report \
        -o ft_traj.xsf --output-structure ft_final.fdf --no-intro > log_ft.txt 2>&1
    check_exit_code $? 0
    for sec in "Model card -- identity" "\[verified\]" "Model card -- level of theory" \
               "Model card -- training data and trained domain" "Trained range" "Model card -- recipe and accuracy" \
               "Errors vs SIESTA" "Model card -- stability, speed and usage" "Model domain along the run"; do
        check_contains "$sec" log_ft.txt
    done
    check_contains "has no test-set validation" log_ft.txt
    check_not_contains "Model card        : none" log_ft.txt
    echo "Testing: an O2 molecule in a 10 Ang box is not treated as a bulk crystal (card criterion: gaps >= 5 Ang)"
    check_contains "Empty space       : gaps >= 5 Ang along 3" log_ft.txt
    check_contains "(includes vacuum)" log_ft.txt
    check_not_contains "Lindemann ratio   :" log_ft.txt
    check_not_contains "Heat capacity     :" log_ft.txt
    check_not_contains "Pressure (GPa)" log_ft.txt
    check_python "the report's SHA-256 is the model file's, and the report file holds the card" "
import hashlib, re, sys
t = open('stb_mlmd_report.txt').read()
h = hashlib.sha256(open('ft/model/ft_compiled.model', 'rb').read()).hexdigest()
m = re.search(r'SHA-256           : ([0-9a-f]{16})\.\.\.', t)
sys.exit(0 if m and h.startswith(m.group(1)) and 'Model card -- identity' in t else 1)
"

    echo "Testing: a structure outside the trained range is flagged (O2 bond squeezed to 0.80 Ang)"
    python3 - <<'PYEOF'
import numpy as np
from stb.core import structure_io
s = structure_io.read_fdf('ft/mlff_config_001/calc.fdf')
p = structure_io.to_pymatgen(s)
c = p.cart_coords.mean(axis=0)
d = p.cart_coords[1] - p.cart_coords[0]
d /= np.linalg.norm(d)
for i, sign in ((0, -1), (1, 1)):
    p.replace(i, p[i].species_string, c + sign * 0.40 * d, coords_are_cartesian=True)
structure_io.write_fdf(structure_io.from_pymatgen(p, species_meta=s.species_meta), 'ft_short.fdf')
PYEOF
    stb-mlmd --file ft_short.fdf $CM --ensemble nve --timestep 0.1 --n-steps 20 --stride 5 --seed 1 \
        -o ft_short.xsf --output-structure ft_short_final.fdf --no-intro > log_ft_short.txt 2>&1
    check_exit_code $? 0
    check_contains "OUTSIDE" log_ft_short.txt
    check_contains "closer than any pair in the training" log_ft_short.txt
    check_contains "Model card check: shortest interatomic distance" log_ft_short.txt

    echo "Testing: elements the model never saw are refused with a clear message, not a traceback"
    stb-mlmd --file si_bulk.fdf $CM --n-steps 10 --no-intro > log_ft_elem.txt 2>&1
    check_exit_code $? 1
    check_contains "trained only on O and cannot evaluate Ge, Si" log_ft_elem.txt
    check_not_contains "Traceback" log_ft_elem.txt
else
    echo -e "${YELLOW}Skipped:${NC} needs mace_run_train and the 7.1 fixtures ($MLFF_FIX)."
fi


# --- 6. --out-format xyz embeds per-frame Epot/Temp/Time ---
echo -e "\n--- Testing --out-format xyz (per-frame Epot/Temp/Time embedding) ---"
stb-mlmd --file si_bulk.fdf --ensemble nvt --n-steps 40 --stride 10 --out-format xyz \
    -o nvt_traj.xyz --output-structure xyz_final.fdf --no-intro > log_xyz.txt 2>&1
check_exit_code $? 0
python3 -c "
import sys
from ase.io import read
frames = read('nvt_traj.xyz', index=':')
sys.exit(0 if len(frames) == 4 and all(
    'Epot' in a.info and 'Temp' in a.info and 'Time' in a.info for a in frames
) else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} xyz trajectory has 4 frames with Epot/Temp/Time embedded"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} xyz trajectory missing expected frames/fields"
    FAIL=$((FAIL+1))
fi


# --- 7. --save-report ---
echo -e "\n--- Testing --save-report ---"
rm -f stb_mlmd_report.txt
stb-mlmd --file si_bulk.fdf --n-steps 10 --save-report --no-intro > log_savereport.txt 2>&1
check_success stb_mlmd_report.txt
check_contains "STB-MLMD REPORT" stb_mlmd_report.txt
rm -f stb_mlmd_report.txt


# --- 8. Robustness / errors ---
echo -e "\n--- Testing error and robustness cases ---"

echo "Testing: missing --file"
stb-mlmd --no-intro > log_missing_file.txt 2>&1
check_exit_code $? 2

echo "Testing: nonexistent input file"
stb-mlmd --file nope.fdf --no-intro > log_missing_input.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_input.txt

echo "Testing: a long timestep and --show without a display end up in the WARNINGS section"
MPLBACKEND=Agg stb-mlmd --file si_bulk.fdf --n-steps 20 --stride 10 --timestep 3 --show -o warn.xsf \
    --output-structure warn.fdf --no-intro > log_warn.txt 2>&1
check_exit_code $? 0
check_python "both warnings are listed after [10] WARNINGS" "
import re, sys
t = re.sub(r'\x1b\[[0-9;]*m', '', open('log_warn.txt').read()).split('[10] WARNINGS', 1)[1]
sys.exit(0 if 'a 3 fs timestep is long' in t and '--show needs a display' in t else 1)
"

echo "Testing: --version"
stb-mlmd --version > log_version.txt 2>&1
check_contains "stb-mlmd" log_version.txt

echo "Testing: --help documents --ensemble, --custom-model, --friction, --pressure, --equilibration-steps, --taut/--taup/--compressibility, --supercell"
stb-mlmd --help > log_help.txt 2>&1
check_contains "ensemble" log_help.txt
check_contains "custom-model" log_help.txt
check_contains "custom-model-card" log_help.txt
check_contains "friction" log_help.txt
check_contains "pressure" log_help.txt
check_contains "equilibration-steps" log_help.txt
check_contains "taut" log_help.txt
check_contains "taup" log_help.txt
check_contains "compressibility" log_help.txt
check_contains "supercell" log_help.txt

echo "Testing: --taut/--taup/--compressibility are accepted on an NPT run"
stb-mlmd --file si_bulk.fdf --ensemble npt --n-steps 20 --stride 10 \
    --taut 50 --taup 300 --compressibility 1.0e-6 --no-intro > log_nptparams.txt 2>&1
check_exit_code $? 0

echo "Testing: --compressibility is in 1/bar and converted for ASE -- the barostat moves the cell"
check_contains "Compressibility, 1/bar" log_help.txt
python3 -c "
import sys
import numpy as np
from ase import units
from ase.build import bulk
from stb.mlmd import build_dynamics
dyn = build_dynamics(bulk('Si', 'diamond', 5.43), 'npt', 1.0, 300, 0.01, 0.0, 1, compressibility_bar=1e-6)
sys.exit(0 if np.isclose(dyn.compressibility, 1e-6 / units.bar) else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} compressibility handed to NPTBerendsen as 1e-6 / units.bar"; PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} compressibility not converted from 1/bar"; FAIL=$((FAIL+1))
fi
stb-mlmd --file si_bulk.fdf --ensemble npt --pressure 20 --temperature 300 --n-steps 300 --stride 50 \
    --compressibility 1.0e-6 -o npt_p20.xsf --output-structure npt_p20.fdf --no-intro > log_npt_p20.txt 2>&1
check_exit_code $? 0
python3 -c "
import sys, numpy as np
from stb.core import structure_io as s
v0 = abs(np.linalg.det(s.read_fdf('si_bulk.fdf').lattice))
v1 = abs(np.linalg.det(s.read_fdf('npt_p20.fdf').lattice))
print(f'   V(20 GPa)/V0 = {v1 / v0:.3f}')
sys.exit(0 if v1 / v0 < 0.95 else 1)
"
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} at 20 GPa the Si cell shrinks by more than 5 % (it stayed put when the unit was wrong)"; PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} the NPT cell did not shrink under 20 GPa"; FAIL=$((FAIL+1))
fi


# --- 9. Interactive path (stb-suite, shortcut 7.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 7.2) ---"
rm -f interactive_traj.xsf md_final.fdf
printf '7.2\nsi_bulk.fdf\n\nsmall\n1 1 2\n1\n300\n0\n30\n1.0\nn\ninteractive_traj.xsf\ny\nn\nn\nn\nn\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_exit_code $? 0
check_success interactive_traj.xsf
check_success md_final.fdf
check_success si_bulk_md_dashboard.png
check_python "the menu's supercell answer (1 1 2) reached stb-mlmd: md_final.fdf holds 128 atoms" "
import sys
from stb.core import structure_io as s
sys.exit(0 if len(s.read_fdf('md_final.fdf').atoms) == 128 else 1)
"


popd > /dev/null

# --- 10. Summary ---
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
