#!/bin/bash

# --- Setup ---
# Smoke test for stb-ramanModes (Raman Spectrum, Stage 2: Modes & Optical
# Displacements, item 4.11.2). Chains a real stb-raman (Stage 1) run to get
# a genuine raman_study/phonon_disp/ tree, then fabricates disp-*/Sn3O4.FA
# force files (same 14-atom Sn3O4 system/format as
# test/4-workflow/4-phonons/analysis/Sn3O4.FA -- physically meaningless but
# format-correct, enough to exercise FORCE_SETS extraction and the
# Gamma-mode/Optical-displacement wiring).
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREP_DIR="$(cd "$FIXTURE_DIR/../prep" && pwd)"
TEST_DIR="$FIXTURE_DIR/test_files"

# Output colors
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

PASS=0
FAIL=0

# --- Check helpers ---

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

check_exit_code() {
    if [ "$1" -eq "$2" ]; then
        echo -e "   -> ${GREEN}Verified:${NC} exit code $1 (expected $2)"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} exit code $1 (expected $2)"
        FAIL=$((FAIL+1))
    fi
}

# Runs stb-raman (Stage 1, tested separately) then drops a fabricated,
# per-folder-scaled Sn3O4.FA into every disp-*/ -- same "each folder gets a
# slightly different scale factor" trick as test/4-workflow/4-phonons/
# analysis/test.sh, for the same reason: identical forces everywhere make
# the dynamical matrix perfectly degenerate and crash phonopy's own code.
make_phonon_disp() {
    rm -rf raman_study
    stb-raman -s structure.fdf -c calc.fdf -p . -dim 1 1 1 --no-intro > /dev/null 2>&1
    python3 - <<'PYEOF'
import glob
with open("Sn3O4.FA") as f:
    lines = f.readlines()
header, rows = lines[0], lines[1:]
for i, d in enumerate(sorted(glob.glob("raman_study/phonon_disp/disp-*")), start=1):
    scale = 1.0 + 0.01 * (i % 17)
    with open(f"{d}/Sn3O4.FA", "w") as out:
        out.write(header)
        for row in rows:
            parts = row.split()
            idx = parts[0]
            vals = [float(v) * scale for v in parts[1:]]
            out.write(f"{idx:>6}  {vals[0]: .9E}  {vals[1]: .9E}  {vals[2]: .9E}\n")
PYEOF
}


# --- 1. Preparation ---
echo "--- Starting tester for STB-RamanModes stage 2: modes & optical displacements (item 4.11.2) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$PREP_DIR/structure.fdf" "$TEST_DIR/"
cp "$PREP_DIR/calc.fdf" "$TEST_DIR/"
cp "$PREP_DIR/O.psf" "$TEST_DIR/"
cp "$PREP_DIR/Sn.psf" "$TEST_DIR/"
cp "$FIXTURE_DIR/Sn3O4.FA" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Missing Stage-1 output guard ---
echo -e "\n--- Testing the missing-Stage-1-output guard ---"
rm -rf raman_study
stb-ramanModes --directory raman_study --calc calc.fdf --no-intro > log_no_stage1.txt 2>&1
check_exit_code $? 1
check_contains "run stb-raman" log_no_stage1.txt


# --- 3. Normal generation: 3 selected modes ---
echo -e "\n--- Testing default generation (--modes 1 2 3) ---"
make_phonon_disp
stb-ramanModes --directory raman_study --calc calc.fdf --modes 1 2 3 --no-intro > log_basic.txt 2>&1
check_exit_code $? 0
check_contains "Force files found : 84" log_basic.txt
check_contains "FORCE_SETS       : generated successfully" log_basic.txt
check_contains "Dimensionality    : 3D (bulk material)" log_basic.txt
check_contains "Non-acoustic Gamma modes : 39" log_basic.txt
check_contains "Selected modes    : 3" log_basic.txt
check_contains "Folders to write  : 3 modes x 2 signs x 3 axes = 18" log_basic.txt
check_success raman_study/optical_disp/mode_01_plus_x/structure.fdf
check_success raman_study/optical_disp/mode_01_plus_x/calc.fdf
check_success raman_study/optical_disp/mode_01_minus_z/structure.fdf
check_success raman_study/optical_disp/mode_01_plus_x/O.psf
check_success raman_study/optical_disp/mode_01_plus_x/Sn.psf
check_contains "OpticalCalculation T" raman_study/optical_disp/mode_01_plus_x/calc.fdf
check_contains "Optical.Vector" raman_study/optical_disp/mode_01_plus_x/calc.fdf
check_contains "MD.TypeOfRun.*CG" raman_study/optical_disp/mode_01_plus_x/calc.fdf
check_contains "MD.NumCGsteps.*0" raman_study/optical_disp/mode_01_plus_x/calc.fdf
check_contains "MODE_TABLE" raman_study/raman_stage2.txt
check_contains "\[2\] PSEUDOPOTENTIALS" raman_study/raman_stage2.txt
check_contains "Elements needed   : O, Sn" raman_study/raman_stage2.txt
check_contains "reused from Stage 1's disp-\* folders" raman_study/raman_stage2.txt
check_contains "Found all required" raman_study/raman_stage2.txt
check_contains "O.psf" raman_study/raman_stage2.txt
check_contains "Sn.psf" raman_study/raman_stage2.txt

echo "Testing: mode_01_plus_x and mode_01_plus_y share the SAME displaced structure (only Optical.Vector differs)"
python3 -c "
from stb.core import structure_io
import numpy as np
sx = structure_io.read_fdf('raman_study/optical_disp/mode_01_plus_x/structure.fdf')
sy = structure_io.read_fdf('raman_study/optical_disp/mode_01_plus_y/structure.fdf')
px = np.array([p for _, p in sx.atoms])
py = np.array([p for _, p in sy.atoms])
assert np.allclose(px, py), 'plus_x and plus_y structures should be identical'
print('OK')
" > log_structure_check.txt 2>&1
check_contains "OK" log_structure_check.txt

echo "Testing: mode_01_plus_x and mode_01_minus_x differ (opposite displacement sign), same cell"
python3 -c "
from stb.core import structure_io
import numpy as np
sp = structure_io.read_fdf('raman_study/optical_disp/mode_01_plus_x/structure.fdf')
sm = structure_io.read_fdf('raman_study/optical_disp/mode_01_minus_x/structure.fdf')
pp = np.array([p for _, p in sp.atoms])
pm = np.array([p for _, p in sm.atoms])
assert not np.allclose(pp, pm), 'plus_x and minus_x should differ'
assert np.allclose(sp.lattice, sm.lattice), 'cell should be unchanged by the mode displacement'
print('OK')
" > log_sign_check.txt 2>&1
check_contains "OK" log_sign_check.txt


# --- 4. --freq-min/--freq-max filtering ---
echo -e "\n--- Testing --freq-min/--freq-max filtering ---"
make_phonon_disp
stb-ramanModes --directory raman_study --calc calc.fdf --freq-min 1.0 --no-intro > log_freqfilter.txt 2>&1
check_exit_code $? 0
check_contains "Selected modes    : 1" log_freqfilter.txt


# --- 4b. Missing pseudopotentials in Stage 1's disp-* folders -> hard error, not silent ---
echo -e "\n--- Testing the missing-pseudopotential guard (no silent 0-pseudo folders) ---"
make_phonon_disp
rm -f raman_study/phonon_disp/disp-*/O.psf raman_study/phonon_disp/disp-*/Sn.psf
stb-ramanModes --directory raman_study --calc calc.fdf --modes 1 --no-intro > log_missing_pseudo.txt 2>&1
check_exit_code $? 1
check_contains "Missing pseudopotential" log_missing_pseudo.txt
check_contains "O, Sn" log_missing_pseudo.txt

echo "Testing: -p/--pseudo-dir override recovers from that same missing-pseudo situation"
stb-ramanModes --directory raman_study --calc calc.fdf --modes 1 -p . --no-intro > log_pseudo_override.txt 2>&1
check_exit_code $? 0
check_contains "-p/--pseudo-dir override" log_pseudo_override.txt
check_success raman_study/optical_disp/mode_01_plus_x/O.psf
check_success raman_study/optical_disp/mode_01_plus_x/Sn.psf


# --- 5. --full-tensor: full symmetric Raman tensor (6 axes instead of 3) ---
echo -e "\n--- Testing --full-tensor ---"
make_phonon_disp
stb-ramanModes --directory raman_study --calc calc.fdf --modes 1 --full-tensor --no-intro > log_fulltensor.txt 2>&1
check_exit_code $? 0
check_contains "Tensor scope      : full symmetric tensor" log_fulltensor.txt
check_contains "Folders to write  : 1 modes x 2 signs x 6 axes = 12" log_fulltensor.txt
check_success raman_study/optical_disp/mode_01_plus_xy/structure.fdf
check_success raman_study/optical_disp/mode_01_plus_xz/structure.fdf
check_success raman_study/optical_disp/mode_01_plus_yz/structure.fdf
check_success raman_study/optical_disp/mode_01_minus_xy/structure.fdf

echo "Testing: mode_01_plus_xy uses the (x+y)/sqrt(2) Optical.Vector direction"
check_contains "0.7071  0.7071  0.0000" raman_study/optical_disp/mode_01_plus_xy/calc.fdf


# --- 5a2. --export-animations: looping .axsf per selected mode ---
echo -e "\n--- Testing --export-animations ---"
make_phonon_disp
stb-ramanModes --directory raman_study --calc calc.fdf --modes 1 2 --export-animations \
    --animation-frames 8 --no-intro > log_animations.txt 2>&1
check_exit_code $? 0
check_contains "Mode animations   : 2 .axsf file(s), 8 frames each." log_animations.txt
check_success raman_study/optical_disp/mode_01_animation.axsf
check_success raman_study/optical_disp/mode_02_animation.axsf

echo "Testing: the .axsf is a real ASE-readable 8-frame animation, atom count matches the structure"
python3 -c "
import ase.io
frames = ase.io.read('raman_study/optical_disp/mode_01_animation.axsf', index=':')
assert len(frames) == 8, f'expected 8 frames, got {len(frames)}'
assert len(frames[0]) == 14, f'expected 14 atoms (Sn3O4), got {len(frames[0])}'
print('OK')
" > log_animations_check.txt 2>&1
check_contains "OK" log_animations_check.txt


# --- 5b. core.raman_symmetry correctness gate: silicon (point group m-3m/Oh) ---
# Independent of the Sn3O4 CLI fixture -- builds real phonons for bulk Si via
# MACE-MP-0 (same model already cached locally by earlier ML tests in this
# suite) and checks classify_modes() against the textbook-known result: the
# acoustic modes are T1u (translation, never Raman-active) and the triply
# degenerate optical mode is T2g (Si's actual, experimentally known Raman
# mode) -- this is the mandatory correctness gate for --use-symmetry's group
# theory (see core/raman_symmetry.py's own docstring for the derivation).
echo -e "\n--- Testing core.raman_symmetry.classify_modes() against known silicon (Oh) result ---"
python3 - <<'PYEOF' > log_symmetry_gate.txt 2>&1
import numpy as np
from ase.build import bulk
from ase import Atoms
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from stb.core import mace_relax
from stb.core.raman_symmetry import classify_modes

atoms = bulk('Si', 'diamond', a=5.43)
calc = mace_relax.get_calculator(model='small', device='cpu')
atoms.calc = calc
mace_relax.relax(atoms, calc, cell_mask=None, fmax=0.01, max_steps=200)

cell = PhonopyAtoms(symbols=atoms.get_chemical_symbols(), cell=atoms.get_cell(),
                     scaled_positions=atoms.get_scaled_positions())
phonon = Phonopy(cell, supercell_matrix=[[2, 0, 0], [0, 2, 0], [0, 0, 2]])
phonon.generate_displacements(distance=0.01)

force_sets = []
for scell in phonon.supercells_with_displacements:
    a2 = Atoms(symbols=scell.symbols, positions=scell.positions, cell=scell.cell, pbc=True)
    a2.calc = calc
    force_sets.append(a2.get_forces())
phonon.forces = force_sets
phonon.produce_force_constants()

mode_symmetries, point_group, error = classify_modes(phonon)
assert error is None, f"unexpected error: {error}"
assert point_group == "m-3m", f"expected m-3m, got {point_group}"

acoustic = [ms for ms in mode_symmetries if max(ms.band_indices) < 3]
optical = [ms for ms in mode_symmetries if min(ms.band_indices) >= 3]
assert len(acoustic) == 1 and acoustic[0].label == "T1u" and acoustic[0].is_raman_active is False, \
    f"acoustic modes must be T1u/inactive, got {acoustic}"
assert len(optical) == 1 and optical[0].label == "T2g" and optical[0].is_raman_active is True, \
    f"optical mode must be T2g/active, got {optical}"
print("SYMMETRY_GATE_OK")
PYEOF
check_contains "SYMMETRY_GATE_OK" log_symmetry_gate.txt


# --- 5b-2. core.raman_symmetry correctness gate: graphene (2D, point group D6h) ---
# A vacuum-padded 2D slab -- classify_modes uses the point group spglib
# detects on the padded cell (a standard, practical stand-in for the true
# "layer group" formalism). Checked against two well-known, widely-cited
# graphene Raman facts: the doubly-degenerate in-plane E2g mode (origin of
# graphene's G-band) is Raman-active, and the out-of-plane B1g mode is
# Raman-SILENT (not in D6h's quadratic basis) -- neither is a centrosymmetric
# g/u-only shortcut result, so this also exercises the general (non-g/u)
# path of the reduction formula on a real, physically meaningful case.
echo -e "\n--- Testing core.raman_symmetry.classify_modes() against known graphene (D6h) result ---"
python3 - <<'PYEOF' > log_symmetry_gate_2d.txt 2>&1
import numpy as np
from ase.build import graphene
from ase import Atoms
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from stb.core import mace_relax
from stb.core.raman_symmetry import classify_modes

atoms = graphene(formula='C2', a=2.46, vacuum=10.0)
calc = mace_relax.get_calculator(model='small', device='cpu')
atoms.calc = calc
mace_relax.relax(atoms, calc, cell_mask=[True, True, False, False, False, False],
                  fmax=0.01, max_steps=300)

cell = PhonopyAtoms(symbols=atoms.get_chemical_symbols(), cell=atoms.get_cell(),
                     scaled_positions=atoms.get_scaled_positions())
phonon = Phonopy(cell, supercell_matrix=[[2, 0, 0], [0, 2, 0], [0, 0, 1]])
phonon.generate_displacements(distance=0.01)

force_sets = []
for scell in phonon.supercells_with_displacements:
    a2 = Atoms(symbols=scell.symbols, positions=scell.positions, cell=scell.cell, pbc=True)
    a2.calc = calc
    force_sets.append(a2.get_forces())
phonon.forces = force_sets
phonon.produce_force_constants()

mode_symmetries, point_group, error = classify_modes(phonon)
assert error is None, f"unexpected error: {error}"
assert point_group == "6/mmm", f"expected 6/mmm (D6h), got {point_group}"

optical = [ms for ms in mode_symmetries if min(ms.band_indices) >= 3]
b1g = [ms for ms in optical if ms.label == "B1g"]
e2g = [ms for ms in optical if ms.label == "E2g"]
assert len(b1g) == 1 and b1g[0].is_raman_active is False, f"B1g (out-of-plane) must be silent, got {b1g}"
assert len(e2g) == 1 and e2g[0].is_raman_active is True, f"E2g (G-band) must be active, got {e2g}"
print("SYMMETRY_GATE_2D_OK")
PYEOF
check_contains "SYMMETRY_GATE_2D_OK" log_symmetry_gate_2d.txt


# --- 5b-3. get_gamma_modes 0D fix: H2O molecule (3 translations + 3 rotations) ---
# A genuinely isolated molecule has 6 trivial modes at Gamma (3 translation
# + 3 free rotation for a non-linear molecule like water), not 3 -- verified
# live via MACE-MP-0 during planning: bands 0-2 (translation) and 3-5
# (rotation) are both far below the genuine vibrational bending/stretching
# modes (~44.5/111.6/115.2 THz for water). get_gamma_modes must return
# exactly 3 non-acoustic modes here (the real vibrations), not 6.
echo -e "\n--- Testing get_gamma_modes 0D trivial-mode exclusion (H2O molecule) ---"
python3 - <<'PYEOF' > log_0d_gate.txt 2>&1
import numpy as np
from ase.build import molecule
from ase import Atoms
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from stb.core import mace_relax
from stb.core import kspace
from stb.core.phonon_workflow import get_gamma_modes

atoms = molecule('H2O')
atoms.center(vacuum=8.0)
atoms.pbc = True
calc = mace_relax.get_calculator(model='small', device='cpu')
atoms.calc = calc
mace_relax.relax(atoms, calc, cell_mask=None, fmax=0.01, max_steps=300)

frac = atoms.get_scaled_positions()
vac = kspace.detect_vacuum_axes(frac, np.array(atoms.get_cell()), 5.0)
assert all(vac), f"expected 0D (all axes vacuum-padded), got {vac}"

cell = PhonopyAtoms(symbols=atoms.get_chemical_symbols(), cell=atoms.get_cell(),
                     scaled_positions=atoms.get_scaled_positions())
phonon = Phonopy(cell, supercell_matrix=[[1, 0, 0], [0, 1, 0], [0, 0, 1]])
phonon.generate_displacements(distance=0.01)
force_sets = []
for scell in phonon.supercells_with_displacements:
    a2 = Atoms(symbols=scell.symbols, positions=scell.positions, cell=scell.cell, pbc=True)
    a2.calc = calc
    force_sets.append(a2.get_forces())
phonon.forces = force_sets
phonon.produce_force_constants()

# Without the 0D fix (extra_trivial_tol_thz=None): still 6 "non-acoustic"
# modes -- the 3 rotations wrongly included.
freqs_before, _, n_extra_before = get_gamma_modes(phonon, exclude_acoustic=True)
assert n_extra_before == 0
assert len(freqs_before) == 6, f"expected 6 without the 0D fix, got {len(freqs_before)}"

# With the fix (0D confirmed, as raman_modes.py itself now does): exactly
# the 3 real vibrations, all well above the rotational tolerance.
freqs_after, band_idx_after, n_extra_after = get_gamma_modes(
    phonon, exclude_acoustic=True, extra_trivial_tol_thz=5.0)
assert n_extra_after == 3, f"expected 3 extra trivial modes excluded, got {n_extra_after}"
assert len(freqs_after) == 3, f"expected 3 real vibrations, got {len(freqs_after)}"
assert all(f > 5.0 for f in freqs_after), f"leftover modes must be real vibrations, got {freqs_after}"
print("ZERO_D_FIX_OK")
PYEOF
check_contains "ZERO_D_FIX_OK" log_0d_gate.txt


# --- 5b-4. core.raman_symmetry.tensor_form() correctness gate: water (C2v/mm2) ---
# Reuses the same relaxed H2O/MACE setup as the 0D gate above. Water's 3
# real vibrational modes are known (textbook C2v triatomic result) to be
# 2xA1 + 1xB2 -- A1 has multiplicity 2 in the quadratic representation
# (xx+yy and zz are BOTH independently A1-invariant, so tensor_form must
# correctly decline to reduce it, at any --full-tensor scope), while B2's
# only invariant is the OFF-DIAGONAL xy-type combination (yz here, given
# this molecule's orientation) -- multiplicity exactly 1, so tensor_form
# must find a reduction, but ONLY when full_tensor=True (yz isn't among
# the diagonal-only x/y/z probe set, so there's nothing to reduce to
# without it -- correctly declining is just as important a check here as
# succeeding).
echo -e "\n--- Testing core.raman_symmetry.tensor_form() against known water (C2v) result ---"
python3 - <<'PYEOF' > log_tensor_form_gate.txt 2>&1
import numpy as np
from ase.build import molecule
from ase import Atoms
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from stb.core import mace_relax
from stb.core.raman_symmetry import classify_modes, tensor_form

atoms = molecule('H2O')
atoms.center(vacuum=8.0)
atoms.pbc = True
calc = mace_relax.get_calculator(model='small', device='cpu')
atoms.calc = calc
mace_relax.relax(atoms, calc, cell_mask=None, fmax=0.01, max_steps=300)

cell = PhonopyAtoms(symbols=atoms.get_chemical_symbols(), cell=atoms.get_cell(),
                     scaled_positions=atoms.get_scaled_positions())
phonon = Phonopy(cell, supercell_matrix=[[1, 0, 0], [0, 1, 0], [0, 0, 1]])
phonon.generate_displacements(distance=0.01)
force_sets = []
for scell in phonon.supercells_with_displacements:
    a2 = Atoms(symbols=scell.symbols, positions=scell.positions, cell=scell.cell, pbc=True)
    a2.calc = calc
    force_sets.append(a2.get_forces())
phonon.forces = force_sets
phonon.produce_force_constants()

mode_symmetries, point_group, error = classify_modes(phonon)
assert error is None and point_group == "mm2", f"expected mm2 (C2v), got {point_group}"

vibrational = [ms for ms in mode_symmetries if min(ms.band_indices) >= 6]
a1_bands = [ms.band_indices[0] for ms in vibrational if ms.label == "A1"]
b2_bands = [ms.band_indices[0] for ms in vibrational if ms.label == "B2"]
assert len(a1_bands) == 2 and len(b2_bands) == 1, f"expected 2xA1+1xB2, got A1={a1_bands} B2={b2_bands}"

for band_idx in a1_bands:
    assert tensor_form(phonon, band_idx, full_tensor=False) is None
    assert tensor_form(phonon, band_idx, full_tensor=True) is None, \
        "A1 has multiplicity 2 in the quadratic rep -- must never be reduced"

band_idx = b2_bands[0]
assert tensor_form(phonon, band_idx, full_tensor=False) is None, \
    "B2's only invariant is off-diagonal -- no reduction possible from x/y/z alone"
result = tensor_form(phonon, band_idx, full_tensor=True)
assert result is not None, "B2 has multiplicity 1 -- must reduce once xy/xz/yz are available"
assert result["probe_axis"] == "yz", f"expected probe_axis 'yz', got {result['probe_axis']}"
t0 = result["T0"]
assert abs(t0[0, 0]) < 1e-6 and abs(t0[1, 1]) < 1e-6 and abs(t0[2, 2]) < 1e-6, \
    f"B2 tensor form must have zero diagonal, got diag={np.diag(t0)}"
assert abs(t0[1, 2]) > 1e-3 and abs(t0[0, 1]) < 1e-6 and abs(t0[0, 2]) < 1e-6, \
    f"B2 tensor form must be pure yz, got {t0}"
print("TENSOR_FORM_GATE_OK")
PYEOF
check_contains "TENSOR_FORM_GATE_OK" log_tensor_form_gate.txt


# --- 5c. --use-symmetry CLI wiring (Sn3O4 fixture) ---
# The fabricated Sn3O4.FA data (see make_phonon_disp above) deliberately
# scales forces differently per disp-*/ folder to avoid a degenerate
# dynamical matrix -- as a side effect it has no real crystal symmetry
# (point group 1), so this only exercises the WIRING (new report section,
# no crash, correctly finds nothing to skip when there's no symmetry to
# exploit) -- the actual group-theory correctness is 5b's job, above.
echo -e "\n--- Testing --use-symmetry CLI wiring ---"
make_phonon_disp
stb-ramanModes --directory raman_study --calc calc.fdf --use-symmetry --no-intro > log_symmetry.txt 2>&1
check_exit_code $? 0
check_contains "\[1b\] SYMMETRY ANALYSIS" log_symmetry.txt
check_contains "Point group       : 1" log_symmetry.txt
check_contains "Symmetry-forbidden (Raman-inactive) : 0/39" log_symmetry.txt

echo "Testing: --use-symmetry + explicit --modes is informational only, never auto-skips"
stb-ramanModes --directory raman_study --calc calc.fdf --use-symmetry --modes 1 2 --no-intro \
    > log_symmetry_explicit.txt 2>&1
check_exit_code $? 0
check_contains "informational only" log_symmetry_explicit.txt
check_contains "Selected modes    : 2" log_symmetry_explicit.txt


# --- 5d. --reduce-components CLI wiring (Sn3O4 fixture, point group 1 -> no-op) ---
# Same caveat as the --use-symmetry wiring test above: this fabricated
# fixture has no real crystal symmetry (point group 1, where every
# quadratic function has multiplicity 6, never 1), so no mode actually
# qualifies for reduction here -- this only proves the flag doesn't crash
# and correctly finds nothing to reduce. The real group-theory correctness
# (including a genuine positive reduction) is 5b-4's job, above.
echo -e "\n--- Testing --reduce-components CLI wiring (no-op on this asymmetric fixture) ---"
stb-ramanModes --directory raman_study --calc calc.fdf --modes 1 2 --reduce-components --no-intro \
    > log_reduce.txt 2>&1
check_exit_code $? 0
check_contains "Folders to write  : 2 modes x 2 signs x 3 axes = 12" log_reduce.txt
check_success raman_study/optical_disp/mode_01_plus_x/structure.fdf

echo "Testing: no mode_XX_tensor_form.json sidecar is written when nothing was reduced"
python3 -c "
import glob
sidecars = glob.glob('raman_study/optical_disp/*_tensor_form.json')
assert sidecars == [], f'unexpected sidecar(s): {sidecars}'
print('OK')
" > log_reduce_check.txt 2>&1
check_contains "OK" log_reduce_check.txt


# --- 5e. --skip-degenerate CLI wiring (Sn3O4 fixture, point group 1 -> no-op) ---
# Same caveat as 5c/5d: no real degeneracy in this fabricated asymmetric
# fixture, so this only proves the flag doesn't crash and correctly finds
# no group to collapse. Real degenerate-group correctness (an actual
# doubly-degenerate mode getting collapsed to 1 representative) is 5f's
# job, below, on real graphene E2g phonons.
echo -e "\n--- Testing --skip-degenerate CLI wiring (no-op on this asymmetric fixture) ---"
stb-ramanModes --directory raman_study --calc calc.fdf --modes 1 2 --skip-degenerate --no-intro \
    > log_skipdeg_noop.txt 2>&1
check_exit_code $? 0
check_contains "Folders to write  : 2 modes x 2 signs x 3 axes = 12" log_skipdeg_noop.txt

echo "Testing: no DERIVED row is written when nothing was collapsed"
python3 -c "
from stb.raman_analysis import read_mode_table, group_by_mode
rows = read_mode_table('raman_study/raman_stage2.txt')
modes = group_by_mode(rows)
derived = {k: v['derived_from'] for k, v in modes.items() if v['derived_from'] is not None}
assert derived == {}, f'unexpected derived mode(s): {derived}'
print('OK')
" > log_skipdeg_noop_check.txt 2>&1
check_contains "OK" log_skipdeg_noop_check.txt


# --- 5f. --skip-degenerate correctness gate: graphene E2g (real MACE-MP-0 phonons) ---
# Real doubly-degenerate case, same graphene fixture/known result as 5b-2
# (in-plane E2g is Raman-active, out-of-plane B1g is silent) -- built
# directly as an ML-sourced (embedded force_constants) phonopy_disp.yaml,
# same has_embedded_fc path stb-mlphonons produces, so no SIESTA .FA/
# FORCE_SETS extraction or real stb-raman Stage 1 run is needed.
echo -e "\n--- Testing --skip-degenerate against a real doubly-degenerate mode (graphene E2g) ---"
rm -rf raman_study
python3 - <<'PYEOF' > log_skipdeg_fixture.txt 2>&1
import os
from ase.build import graphene
from ase import Atoms
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from stb.core import mace_relax

atoms = graphene(formula='C2', a=2.46, vacuum=10.0)
calc = mace_relax.get_calculator(model='small', device='cpu')
atoms.calc = calc
mace_relax.relax(atoms, calc, cell_mask=[True, True, False, False, False, False],
                  fmax=0.01, max_steps=300)

cell = PhonopyAtoms(symbols=atoms.get_chemical_symbols(), cell=atoms.get_cell(),
                     scaled_positions=atoms.get_scaled_positions())
phonon = Phonopy(cell, supercell_matrix=[[2, 0, 0], [0, 2, 0], [0, 0, 1]])
phonon.generate_displacements(distance=0.01)

force_sets = []
for scell in phonon.supercells_with_displacements:
    a2 = Atoms(symbols=scell.symbols, positions=scell.positions, cell=scell.cell, pbc=True)
    a2.calc = calc
    force_sets.append(a2.get_forces())
phonon.forces = force_sets
phonon.produce_force_constants()

os.makedirs("raman_study/phonon_disp", exist_ok=True)
phonon.save(filename="raman_study/phonon_disp/phonopy_disp.yaml", settings={'force_constants': True})
print("FIXTURE_OK")
PYEOF
check_contains "FIXTURE_OK" log_skipdeg_fixture.txt
touch C.psf

echo "Testing: --skip-degenerate alone (3 selected modes: B1g singlet + E2g pair -> 2 written)"
stb-ramanModes --directory raman_study --calc calc.fdf --pseudo-dir . --skip-degenerate --no-intro \
    > log_skipdeg.txt 2>&1
check_exit_code $? 0
check_contains "Degenerate groups : 1 mode(s) skipped" log_skipdeg.txt

python3 -c "
import glob
real_dirs = glob.glob('raman_study/optical_disp/mode_*_plus_*') + glob.glob('raman_study/optical_disp/mode_*_minus_*')
modes_with_folders = sorted(set(int(d.split('/')[-1].split('_')[1]) for d in real_dirs))
assert len(modes_with_folders) == 2, f'expected 2 modes with real folders, got {modes_with_folders}'
derived_dirs = glob.glob('raman_study/optical_disp/mode_*_derived')
assert derived_dirs == [], f'derived mode must not get a folder: {derived_dirs}'
print('OK')
" > log_skipdeg_folders.txt 2>&1
check_contains "OK" log_skipdeg_folders.txt

echo "Testing: MODE_TABLE records the derived mode against a fully-computed representative"
python3 -c "
from stb.raman_analysis import read_mode_table, group_by_mode
rows = read_mode_table('raman_study/raman_stage2.txt')
modes = group_by_mode(rows)
derived = {k: v['derived_from'] for k, v in modes.items() if v['derived_from'] is not None}
assert len(derived) == 1, f'expected exactly 1 derived mode, got {derived}'
rep_k = list(derived.values())[0]
assert modes[rep_k]['derived_from'] is None, 'representative must not itself be derived'
assert len(modes[rep_k]['folders']) == 6, f\"representative should have 6 folders (2 signs x 3 diagonal axes), got {len(modes[rep_k]['folders'])}\"
print('OK')
" > log_skipdeg_table.txt 2>&1
check_contains "OK" log_skipdeg_table.txt

echo "Testing: --skip-degenerate + --use-symmetry together (E2g pair survives, B1g filtered out)"
rm -rf raman_study/optical_disp raman_study/raman_stage2.txt
stb-ramanModes --directory raman_study --calc calc.fdf --pseudo-dir . --use-symmetry --skip-degenerate \
    --no-intro > log_skipdeg_sym.txt 2>&1
check_exit_code $? 0
check_contains "Selected modes    : 2" log_skipdeg_sym.txt
check_contains "E2g" log_skipdeg_sym.txt
check_contains "Degenerate groups : 1 mode(s) skipped" log_skipdeg_sym.txt


# --- 6. Error cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: missing calc file"
stb-ramanModes --directory raman_study --calc does_not_exist.fdf --no-intro > log_missing_calc.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_calc.txt

echo "Testing: --version"
stb-ramanModes --version > log_version.txt 2>&1
check_contains "stb-ramanModes" log_version.txt

echo "Testing: --help documents --modes/--optical-mesh/--optical-broaden/--full-tensor/--use-symmetry/--vacuum-gap"
stb-ramanModes --help > log_help.txt 2>&1
check_contains "modes" log_help.txt
check_contains "optical-mesh" log_help.txt
check_contains "optical-broaden" log_help.txt
check_contains "full-tensor" log_help.txt
check_contains "use-symmetry" log_help.txt
check_contains "vacuum-gap" log_help.txt
check_contains "rotational-mode-tol" log_help.txt
check_contains "export-animations" log_help.txt
check_contains "reduce-components" log_help.txt
check_contains "skip-degenerate" log_help.txt


# --- 7. Interactive path (stb-suite, shortcut 4.11.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.11.2) ---"

echo "Testing: navigate 4.11.2 -> modes 1 2, skip full-tensor, skip-degenerate ON, animations ON -> quit"
make_phonon_disp
{
  echo "4.11.2"
  echo "raman_study"    # run_dir
  echo "calc.fdf"       # calc_file
  echo "1 2"            # modes
  echo ""               # full_tensor_choice (default N)
  echo ""               # reduce_components_choice (default N)
  echo "y"              # skip_degenerate_choice -> Y (no-op on this asymmetric fixture)
  echo ""               # displacement (default 0.02)
  echo "y"              # animations_choice -> Y
  echo ""               # animation_frames (default 20)
  echo ""               # show_advanced (default -> skip)
  echo ""               # press enter to continue
  echo "0"              # quit stage submenu
} | stb-suite > log_menu.txt 2>&1
check_contains "Folders to write  : 2 modes x 2 signs x 3 axes = 12" log_menu.txt
check_success raman_study/optical_disp/mode_01_plus_x/structure.fdf
check_success raman_study/optical_disp/mode_01_animation.axsf

echo "Testing: navigate 4.11.2 -> modes 1, full-tensor ON -> quit"
make_phonon_disp
{
  echo "4.11.2"
  echo "raman_study"    # run_dir
  echo "calc.fdf"       # calc_file
  echo "1"              # modes
  echo "y"              # full_tensor_choice -> Y
  echo ""               # reduce_components_choice (default N)
  echo ""               # skip_degenerate_choice (default N)
  echo ""               # displacement (default 0.02)
  echo ""               # animations_choice (default N)
  echo ""               # show_advanced (default -> skip)
  echo ""               # press enter to continue
  echo "0"              # quit stage submenu
} | stb-suite > log_menu_fulltensor.txt 2>&1
check_contains "Folders to write  : 1 modes x 2 signs x 6 axes = 12" log_menu_fulltensor.txt
check_success raman_study/optical_disp/mode_01_plus_xy/structure.fdf

echo "Testing: navigate 4.11.2 -> default modes (empty), use-symmetry ON -> quit"
make_phonon_disp
{
  echo "4.11.2"
  echo "raman_study"    # run_dir
  echo "calc.fdf"       # calc_file
  echo ""               # modes (default -> every non-acoustic mode; unlocks the use-symmetry prompt)
  echo ""               # full_tensor_choice (default N)
  echo "y"              # use_symmetry_choice -> Y
  echo ""               # reduce_components_choice (default N)
  echo ""               # skip_degenerate_choice (default N)
  echo ""               # displacement (default 0.02)
  echo ""               # animations_choice (default N)
  echo ""               # show_advanced (default -> skip)
  echo ""               # press enter to continue
  echo "0"              # quit stage submenu
} | stb-suite > log_menu_symmetry.txt 2>&1
check_contains "\[1b\] SYMMETRY ANALYSIS" log_menu_symmetry.txt

echo "Testing: navigate 4.11.2 -> advanced settings ON (mesh/broaden/vacuum-gap/rotational-mode-tol)"
make_phonon_disp
{
  echo "4.11.2"
  echo "raman_study"    # run_dir
  echo "calc.fdf"       # calc_file
  echo "1"              # modes
  echo ""               # full_tensor_choice (default N)
  echo ""               # reduce_components_choice (default N)
  echo ""               # skip_degenerate_choice (default N)
  echo ""               # displacement (default 0.02)
  echo ""               # animations_choice (default N)
  echo "y"              # show_advanced -> Y
  echo "8 8 8"           # mesh_input
  echo "0.15"            # optical_broaden
  echo ""                # freq_min (skip)
  echo ""                # freq_max (skip)
  echo ""                # pseudo_dir prompt (skip)
  echo "12.5"             # vacuum_gap
  echo "1.5"              # rotational_mode_tol
  echo ""               # press enter to continue
  echo "0"              # quit stage submenu
} | stb-suite > log_menu_advanced.txt 2>&1
check_contains "Optical mesh      : 8 x 8 x 8" log_menu_advanced.txt
check_contains "Optical broaden   : 0.15 eV" log_menu_advanced.txt
check_success raman_study/optical_disp/mode_01_plus_x/structure.fdf


popd > /dev/null

# --- 7. Summary ---
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
