#!/bin/bash

# --- Setup ---
# Smoke test for stb-irModes (IR Spectrum, Stage 2: Modes & Displacements,
# item 4.12.2). Covers BOTH auto-selected paths:
#   - bulk (3D, Born-effective-charge): chains a real stb-ir (Stage 1) run
#     on the Sn3O4 fixture (genuinely bulk, no vacuum -- confirmed by its
#     own structure.fdf), fabricates disp-*/Sn3O4.FA force files (same
#     per-folder-scaled trick as test/.../11-raman/modes/test.sh).
#   - non-bulk (0D, dipole derivative): builds a real 0D phonon_disp/ via
#     ase.build.molecule('H2O') + MACE-MP-0 (embedded force constants,
#     bypassing SIESTA .FA entirely) -- same fixture recipe already
#     established and reused throughout the Raman modes test.
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

check_not_contains() {
    if grep -q -- "$1" "$2" 2>/dev/null; then
        echo -e "   -> ${RED}Failed:${NC} '$1' unexpectedly found in '$2'"
        FAIL=$((FAIL+1))
    else
        echo -e "   -> ${GREEN}Verified:${NC} '$1' NOT found in '$2' (as expected)"
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

# Runs stb-ir (Stage 1, tested separately) then drops a fabricated,
# per-folder-scaled Sn3O4.FA into every disp-*/ -- same trick as
# test/4-workflow/11-raman/modes/test.sh (identical forces everywhere
# make the dynamical matrix perfectly degenerate and crash phonopy).
make_phonon_disp() {
    rm -rf ir_study
    stb-ir -s structure.fdf -c calc.fdf -p . -dim 1 1 1 --no-intro > /dev/null 2>&1
    python3 - <<'PYEOF'
import glob
with open("Sn3O4.FA") as f:
    lines = f.readlines()
header, rows = lines[0], lines[1:]
for i, d in enumerate(sorted(glob.glob("ir_study/phonon_disp/disp-*")), start=1):
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
echo "--- Starting tester for STB-IRModes stage 2: modes & displacements (item 4.12.2) ---"
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
rm -rf ir_study
stb-irModes --directory ir_study --calc calc.fdf --no-intro > log_no_stage1.txt 2>&1
check_exit_code $? 1
check_contains "run stb-ir" log_no_stage1.txt


# --- 3. Bulk path: default generation (--modes 1 2) ---
echo -e "\n--- Testing bulk (3D) path: default generation (--modes 1 2) ---"
make_phonon_disp
stb-irModes --directory ir_study --calc calc.fdf --modes 1 2 --no-intro > log_bulk_basic.txt 2>&1
check_exit_code $? 0
check_contains "\[Path\] Bulk (3D periodic)" log_bulk_basic.txt
check_success ir_study/born_charge_disp/equilibrium/structure.fdf
check_success ir_study/born_charge_disp/equilibrium/calc.fdf
check_success ir_study/born_charge_disp/mode_01_eigendisplacement.json
check_success ir_study/born_charge_disp/mode_02_eigendisplacement.json
check_success ir_study/ir_stage2.txt
check_contains "MD.TypeOfRun          FC" ir_study/born_charge_disp/equilibrium/calc.fdf
check_contains "BornCharge            T" ir_study/born_charge_disp/equilibrium/calc.fdf
check_contains "%block PolarizationGrids" ir_study/born_charge_disp/equilibrium/calc.fdf

echo "Testing: exactly 1 folder written under born_charge_disp/, MODE_TABLE path column is BULK"
python3 -c "
import glob
folders = [d for d in glob.glob('ir_study/born_charge_disp/*') if '.' not in d.split('/')[-1]]
assert folders == ['ir_study/born_charge_disp/equilibrium'], f'expected exactly 1 folder, got {folders}'
from stb.ir_analysis import read_mode_table
rows = read_mode_table('ir_study/ir_stage2.txt')
assert len(rows) == 2, f'expected 2 MODE_TABLE rows, got {len(rows)}'
assert all(r['path'] == 'BULK' for r in rows), f'expected all BULK, got {[r[\"path\"] for r in rows]}'
print('OK')
" > log_bulk_check.txt 2>&1
check_contains "OK" log_bulk_check.txt


# --- 4. Bulk path ignores mode-selection flags for folder count ---
echo -e "\n--- Testing bulk path: --modes/--freq-min/--use-symmetry never change the folder count ---"
stb-irModes --directory ir_study --calc calc.fdf --modes 1 --no-intro > log_bulk_modes1.txt 2>&1
check_exit_code $? 0
check_contains "1 equilibrium Born-effective-charge SIESTA run" log_bulk_modes1.txt

python3 -c "
import glob
folders = [d for d in glob.glob('ir_study/born_charge_disp/*') if '.' not in d.split('/')[-1]]
assert folders == ['ir_study/born_charge_disp/equilibrium'], f'expected exactly 1 folder, got {folders}'
print('OK')
" > log_bulk_modes1_check.txt 2>&1
check_contains "OK" log_bulk_modes1_check.txt

echo "Testing: --skip-degenerate is a documented no-op on the bulk path"
stb-irModes --directory ir_study --calc calc.fdf --modes 1 2 --skip-degenerate --no-intro \
    > log_bulk_skipdeg.txt 2>&1
check_exit_code $? 0
check_contains "no effect on the bulk path" log_bulk_skipdeg.txt
check_success ir_study/born_charge_disp/mode_01_eigendisplacement.json
check_success ir_study/born_charge_disp/mode_02_eigendisplacement.json


# --- 5. Bulk path: custom --fc-displ/--polarization-grid ---
echo -e "\n--- Testing bulk path: --fc-displ/--polarization-grid ---"
stb-irModes --directory ir_study --calc calc.fdf --modes 1 --fc-displ 0.03 \
    --polarization-grid 8 3 3 3 8 3 3 3 8 --no-intro > log_bulk_custom.txt 2>&1
check_exit_code $? 0
check_contains "MD.FCDispl            0.03 bohr" ir_study/born_charge_disp/equilibrium/calc.fdf
check_contains "8  3  3" ir_study/born_charge_disp/equilibrium/calc.fdf


# --- 6. Missing pseudopotentials guard ---
echo -e "\n--- Testing the missing-pseudopotential guard ---"
make_phonon_disp
mkdir -p empty_pp
stb-irModes --directory ir_study --calc calc.fdf --modes 1 -p empty_pp --no-intro \
    > log_missing_pseudo.txt 2>&1
check_exit_code $? 1
check_contains "Missing pseudopotential" log_missing_pseudo.txt


# --- 7. IR symmetry correctness gates (real MACE-MP-0 phonons) ---
echo -e "\n--- Testing core.ir_symmetry.classify_ir_modes() against known results ---"

echo "Testing: silicon (m-3m) mutual-exclusion gate -- T2g is Raman-active AND IR-inactive"
python3 - <<'PYEOF' > log_ir_symmetry_si.txt 2>&1
from ase.build import bulk
from ase import Atoms
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from stb.core import mace_relax
from stb.core.raman_symmetry import classify_modes
from stb.core.ir_symmetry import classify_ir_modes

atoms = bulk('Si', 'diamond', a=5.43)
calc = mace_relax.get_calculator(model='small', device='cpu')
atoms.calc = calc
mace_relax.relax(atoms, calc, cell_mask=[True] * 6, fmax=0.01, max_steps=300)

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

ms_raman, pg_raman, err_raman = classify_modes(phonon)
ms_ir, pg_ir, err_ir = classify_ir_modes(phonon)
assert err_raman is None and err_ir is None, f"unexpected error: {err_raman} / {err_ir}"
assert pg_raman == "m-3m" and pg_ir == "m-3m", f"expected m-3m, got {pg_raman}/{pg_ir}"

optical_raman = [ms for ms in ms_raman if min(ms.band_indices) >= 3]
optical_ir = [ms for ms in ms_ir if min(ms.band_indices) >= 3]
assert len(optical_raman) == 1 and optical_raman[0].is_raman_active is True, \
    f"T2g must be Raman-active, got {[(m.label, m.is_raman_active) for m in optical_raman]}"
assert len(optical_ir) == 1 and optical_ir[0].is_ir_active is False, \
    f"T2g must be IR-inactive (mutual exclusion), got {[(m.label, m.is_ir_active) for m in optical_ir]}"
print("MUTUAL_EXCLUSION_GATE_OK")
PYEOF
check_contains "MUTUAL_EXCLUSION_GATE_OK" log_ir_symmetry_si.txt

echo "Testing: water (C2v/mm2) -- classify_ir_modes labels agree with classify_modes on the same bands"
python3 - <<'PYEOF' > log_ir_symmetry_h2o.txt 2>&1
from ase.build import molecule
from ase import Atoms
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from stb.core import mace_relax
from stb.core.raman_symmetry import classify_modes
from stb.core.ir_symmetry import classify_ir_modes

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

ms_raman, pg_raman, err_raman = classify_modes(phonon)
ms_ir, pg_ir, err_ir = classify_ir_modes(phonon)
assert err_raman is None and err_ir is None
assert pg_raman == "mm2" and pg_ir == "mm2", f"expected mm2 (C2v), got {pg_raman}/{pg_ir}"

# Vibrational bands only (>= 6, past the 3 translations + 3 rotations) --
# classify_ir_modes must agree with classify_modes' own (already-verified)
# labels for these same bands, since both read the exact same underlying
# phonopy IrReps labels via core.vibrational_symmetry.run_gamma_irreps.
# NOTE: this deliberately does NOT assert absolute is_ir_active TRUTH here
# (e.g. "both A1 modes must be IR-active") -- doing so requires the first
# 3 Gamma bands to be cleanly resolved pure translations, which turned out
# NOT to be reliable for this specific H2O/MACE-small/vacuum-box fixture
# (verified during planning: even at fmax=1e-4 the "acoustic" bands came
# out at ~1.5-1.9 THz, a real MACE-small numerical characteristic of this
# fixture, not a bug -- the exact same latent fragility already affects
# core.raman_symmetry's own chi_v extraction, which is why the existing
# Raman gate for this fixture (5b-4 in .../11-raman/modes/test.sh) also
# only checks labels/self-consistency, never absolute activity truth).
raman_labels = {tuple(ms.band_indices): ms.label for ms in ms_raman if min(ms.band_indices) >= 6}
ir_labels = {tuple(ms.band_indices): ms.label for ms in ms_ir if min(ms.band_indices) >= 6}
assert raman_labels == ir_labels, f"label mismatch: raman={raman_labels} ir={ir_labels}"
assert sorted(raman_labels.values()) == ["A1", "A1", "B2"], f"expected 2xA1+1xB2, got {raman_labels}"
print("LABEL_CONSISTENCY_GATE_OK")
PYEOF
check_contains "LABEL_CONSISTENCY_GATE_OK" log_ir_symmetry_h2o.txt


# --- 8. Non-bulk path: CLI wiring (H2O, real MACE-MP-0 phonons) ---
echo -e "\n--- Testing non-bulk (0D) path: CLI wiring (H2O) ---"
rm -rf ir_study
python3 - <<'PYEOF' > log_nonbulk_fixture.txt 2>&1
import os
from ase.build import molecule
from ase import Atoms
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from stb.core import mace_relax

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

os.makedirs("ir_study/phonon_disp", exist_ok=True)
phonon.save(filename="ir_study/phonon_disp/phonopy_disp.yaml", settings={'force_constants': True})
print("FIXTURE_OK")
PYEOF
check_contains "FIXTURE_OK" log_nonbulk_fixture.txt
touch H.psf

stb-irModes --directory ir_study --calc calc.fdf --pseudo-dir . --no-intro > log_nonbulk_basic.txt 2>&1
check_exit_code $? 0
check_contains "\[Path\] Non-bulk (0D/1D/2D)" log_nonbulk_basic.txt
check_contains "Folders to write  : 6 independent single-point SIESTA runs" log_nonbulk_basic.txt
check_success ir_study/dipole_disp/mode_01_plus/structure.fdf
check_success ir_study/dipole_disp/mode_01_minus/structure.fdf

echo "Testing: no axis suffix on non-bulk folder names (unlike Raman's mode_XX_plus_x/y/z)"
python3 -c "
import glob
dirs = sorted(d.split('/')[-1] for d in glob.glob('ir_study/dipole_disp/mode_*'))
expected = sorted(f'mode_{m:02d}_{s}' for m in (1, 2, 3) for s in ('plus', 'minus'))
assert dirs == expected, f'expected {expected}, got {dirs}'
print('OK')
" > log_nonbulk_dirs_check.txt 2>&1
check_contains "OK" log_nonbulk_dirs_check.txt

echo "Testing: no Optical.Vector/BornCharge leakage into a non-bulk calc.fdf"
check_contains "MD.TypeOfRun          CG" ir_study/dipole_disp/mode_01_plus/calc.fdf
check_not_contains "Optical" ir_study/dipole_disp/mode_01_plus/calc.fdf
check_not_contains "BornCharge" ir_study/dipole_disp/mode_01_plus/calc.fdf

echo "Testing: MODE_TABLE path column is NONBULK"
python3 -c "
from stb.ir_analysis import read_mode_table
rows = read_mode_table('ir_study/ir_stage2.txt')
assert all(r['path'] == 'NONBULK' for r in rows), f'expected all NONBULK, got {[r[\"path\"] for r in rows]}'
print('OK')
" > log_nonbulk_table_check.txt 2>&1
check_contains "OK" log_nonbulk_table_check.txt


# --- 9. --skip-degenerate on the non-bulk path ---
echo -e "\n--- Testing --skip-degenerate CLI wiring on H2O (non-bulk, no real degeneracy -> no-op) ---"
stb-irModes --directory ir_study --calc calc.fdf --pseudo-dir . --skip-degenerate --no-intro \
    > log_skipdeg_h2o.txt 2>&1
check_exit_code $? 0
check_contains "Folders to write  : 6 independent single-point SIESTA runs" log_skipdeg_h2o.txt

echo -e "\n--- Testing --skip-degenerate correctness on a real degenerate pair (graphene E2g, 2D non-bulk) ---"
rm -rf ir_study
python3 - <<'PYEOF' > log_graphene_fixture.txt 2>&1
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

os.makedirs("ir_study/phonon_disp", exist_ok=True)
phonon.save(filename="ir_study/phonon_disp/phonopy_disp.yaml", settings={'force_constants': True})
print("FIXTURE_OK")
PYEOF
check_contains "FIXTURE_OK" log_graphene_fixture.txt
touch C.psf

# 3 optical modes at Gamma for graphene: B1g singlet + E2g doubly-degenerate
# pair -- same fixture/known result already used in .../11-raman/modes/test.sh.
stb-irModes --directory ir_study --calc calc.fdf --pseudo-dir . --skip-degenerate --no-intro \
    > log_graphene_skipdeg.txt 2>&1
check_exit_code $? 0
check_contains "Degenerate groups : 1 mode(s) skipped" log_graphene_skipdeg.txt

python3 -c "
import glob
real_dirs = glob.glob('ir_study/dipole_disp/mode_*_plus') + glob.glob('ir_study/dipole_disp/mode_*_minus')
modes_with_folders = sorted(set(int(d.split('/')[-1].split('_')[1]) for d in real_dirs))
assert len(modes_with_folders) == 2, f'expected 2 modes with real folders (B1g + E2g representative), got {modes_with_folders}'
derived_dirs = glob.glob('ir_study/dipole_disp/mode_*_derived')
assert derived_dirs == [], f'derived mode must not get a folder: {derived_dirs}'
from stb.ir_analysis import read_mode_table, group_by_mode
rows = read_mode_table('ir_study/ir_stage2.txt')
modes = group_by_mode(rows)
derived = {k: v['derived_from'] for k, v in modes.items() if v['derived_from'] is not None}
assert len(derived) == 1, f'expected exactly 1 derived mode, got {derived}'
print('OK')
" > log_graphene_check.txt 2>&1
check_contains "OK" log_graphene_check.txt


# --- 10. Error cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: missing calc file"
stb-irModes --directory ir_study --calc does_not_exist.fdf --no-intro > log_missing_calc.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing_calc.txt

echo "Testing: --version"
stb-irModes --version > log_version.txt 2>&1
check_contains "stb-irModes" log_version.txt

echo "Testing: --help documents --modes/--use-symmetry/--skip-degenerate/--fc-displ/--polarization-grid"
stb-irModes --help > log_help.txt 2>&1
check_contains "modes" log_help.txt
check_contains "use-symmetry" log_help.txt
check_contains "skip-degenerate" log_help.txt
check_contains "fc-displ" log_help.txt
check_contains "polarization-grid" log_help.txt
check_contains "export-animations" log_help.txt


# --- 11. Interactive path (stb-suite, shortcut 4.12.2) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.12.2) ---"

echo "Testing: navigate 4.12.2 -> modes 1 2 (bulk path) -> quit"
make_phonon_disp
{
  echo "4.12.2"
  echo ""               # run_dir (default ir_study)
  echo "calc.fdf"       # calc_file
  echo "1 2"            # modes
  echo ""               # skip_degenerate_choice (default N)
  echo ""               # displacement (default 0.02)
  echo ""               # animations_choice (default N)
  echo ""               # show_advanced (default -> skip)
  echo ""               # press enter to continue
  echo "0"              # quit stage submenu
} | stb-suite > log_menu_bulk.txt 2>&1
check_contains "1 equilibrium Born-effective-charge SIESTA run" log_menu_bulk.txt
check_success ir_study/born_charge_disp/equilibrium/structure.fdf


popd > /dev/null

# --- 12. Summary ---
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
