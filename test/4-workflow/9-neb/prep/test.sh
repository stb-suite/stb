#!/bin/bash

# --- Setup ---
# Smoke test for stb-neb (NEB Prep, item 4.9.1 -- this is now Stage 1 of the
# NEB workflow: the old symmetry-site-enumeration Stage 1, stb-nebSites, was
# removed since restricting NEB endpoints to symmetry-derived candidate sites
# defeated the point of NEB -- stb-neb now accepts any pair of already
# -relaxed structures directly). Now a 4-mode tool (--mode 1-4, default 3)
# instead of a single --ml-neb on/off switch -- see examples/README.md /
# CLAUDE.md for the mode table. Tests that only exercise the common preamble
# (composition/lattice/interpolation mechanics, not mode-specific output)
# deliberately use --mode 4 (no MACE, fastest, no model load) unless they're
# specifically testing a MACE-touching mode.
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$FIXTURE_DIR/test_files"

# Output colors
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

PASS=0
FAIL=0

# --- Check helpers ---

# Checks that file $1 exists and is not empty
check_success() {
    if [ -s "$1" ]; then
        echo -e " ... ${GREEN}OK${NC} (file '$1' created)"
        PASS=$((PASS+1))
    else
        echo -e " ... ${RED}FAIL${NC} (file '$1' was not created)"
        FAIL=$((FAIL+1))
    fi
}

# Checks that file $2 contains (grep -q) pattern $1
check_contains() {
    if grep -q "$1" "$2" 2>/dev/null; then
        echo -e "   -> ${GREEN}Verified:${NC} '$1' found in '$2'"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} '$1' NOT found in '$2'"
        FAIL=$((FAIL+1))
    fi
}

# Checks that $1 (actual exit code) equals $2 (expected exit code)
check_exit_code() {
    if [ "$1" -eq "$2" ]; then
        echo -e "   -> ${GREEN}Verified:${NC} exit code $1 (expected $2)"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} exit code $1 (expected $2)"
        FAIL=$((FAIL+1))
    fi
}


# --- 1. Preparation ---
echo "--- Starting tester for STB-Neb prep (item 4.9.1) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/initial.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/final.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/final_bad_composition.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/final_bad_lattice.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/final_bad_autosort.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/final_many_moved.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/wrap_initial.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/wrap_final.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/calc.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Mode 4 (no MACE, cycle_00/ nesting) -- the common-preamble workhorse ---
echo -e "\n--- Testing --mode 4, -n 5 (interpolation mechanics, no MACE) ---"
stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 5 --mode 4 --no-intro > log_mode4.txt 2>&1
check_exit_code $? 0
check_success cycle_00/image_00/structure.fdf
check_success cycle_00/image_00/calc.fdf
check_success cycle_00/image_04/structure.fdf
check_contains "MD.TypeOfRun.*CG" cycle_00/image_02/calc.fdf
check_contains "MD.NumCGsteps.*0" cycle_00/image_02/calc.fdf
check_success neb_setup.txt
check_contains "\[0\] RUN METADATA" neb_setup.txt
check_contains "\[1\] INTERPOLATION" neb_setup.txt
check_contains "\[2\] IMAGE FOLDERS" neb_setup.txt
check_contains "\[3\] SUMMARY" neb_setup.txt
check_contains "\[4\] CLUSTER SUBMISSION" neb_setup.txt
check_contains "stb-nebCycle --dir . --fmax 0.05 --climb-after 5" neb_setup.txt
check_contains "# MODE: 4" neb_setup.txt
check_contains "MACE-MP-0 path shaping: no" neb_setup.txt
check_contains "IMAGE_TABLE" neb_setup.txt
check_success neb_path.xyz
if grep -q "WARNING" log_mode4.txt; then
    echo -e "   -> ${RED}Failed:${NC} unexpected path-quality WARNING for a well-behaved 5-image band"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no path-quality warnings for a well-behaved band"
    PASS=$((PASS+1))
fi

echo "Testing: reaction_coord is non-decreasing and image_00 is 0.0"
python3 -c "
import re
rows = []
in_table = False
with open('neb_setup.txt') as f:
    for line in f:
        if line.startswith('# IMAGE_TABLE'):
            in_table = True
            continue
        if not in_table:
            continue
        if line.startswith('#') or not line.strip():
            continue
        parts = line.split()
        rows.append((parts[0], int(parts[1]), float(parts[2])))
rows.sort(key=lambda r: r[1])
assert rows[0][2] == 0.0, f'image_00 reaction_coord should be 0.0, got {rows[0][2]}'
coords = [r[2] for r in rows]
assert coords == sorted(coords), f'reaction_coord not non-decreasing: {coords}'
print('OK')
" > log_coord_check.txt 2>&1
check_contains "OK" log_coord_check.txt


# --- 2b. config_extra.fdf / --force-spin / --force-vdw / --force-dipole
#     (all default ON) -- reuses the cycle_00/ written by Section 2 above,
#     since defaults are exactly what that run already used. ---
echo -e "\n--- Testing config_extra.fdf defaults (all three force flags ON by default) ---"
check_success cycle_00/image_02/config_extra.fdf
check_contains "Spin.*polarized" cycle_00/image_02/config_extra.fdf
check_contains "Slab.DipoleCorrection" cycle_00/image_02/config_extra.fdf
check_contains "DFTD3" cycle_00/image_02/config_extra.fdf
check_contains "%include config_extra.fdf" cycle_00/image_02/calc.fdf
echo "Testing: %include config_extra.fdf is the first line of calc.fdf"
first_line=$(head -1 cycle_00/image_02/calc.fdf)
if [[ "$first_line" == *"%include config_extra.fdf"* ]]; then
    echo -e "   -> ${GREEN}Verified:${NC} %include is the first line"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} %include is not the first line (got: '$first_line')"
    FAIL=$((FAIL+1))
fi
rm -rf cycle_00 neb_setup.txt

echo -e "\n--- Testing --no-force-spin --no-force-vdw --no-force-dipole ---"
stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 5 --mode 4 --no-force-spin --no-force-vdw \
    --no-force-dipole --no-intro > log_noforce.txt 2>&1
check_exit_code $? 0
echo "Testing: config_extra.fdf still exists (deliberately empty, not skipped)"
if [ -f cycle_00/image_02/config_extra.fdf ]; then
    echo -e "   -> ${GREEN}Verified:${NC} file 'cycle_00/image_02/config_extra.fdf' created (empty)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} file 'cycle_00/image_02/config_extra.fdf' was not created"
    FAIL=$((FAIL+1))
fi
if grep -qE "Spin.*polarized|Slab.DipoleCorrection|DFTD3" cycle_00/image_02/config_extra.fdf; then
    echo -e "   -> ${RED}Failed:${NC} config_extra.fdf should be empty of forced directives"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no forced directives in config_extra.fdf with all --no-force-*"
    PASS=$((PASS+1))
fi
check_contains "%include config_extra.fdf" cycle_00/image_02/calc.fdf
rm -rf cycle_00 neb_setup.txt

echo -e "\n--- Testing mixed flags (--force-spin --no-force-vdw --force-dipole) ---"
stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 5 --mode 4 --force-spin --no-force-vdw \
    --force-dipole --no-intro > log_mixed.txt 2>&1
check_exit_code $? 0
check_contains "Spin.*polarized" cycle_00/image_02/config_extra.fdf
check_contains "Slab.DipoleCorrection" cycle_00/image_02/config_extra.fdf
if grep -q "DFTD3" cycle_00/image_02/config_extra.fdf; then
    echo -e "   -> ${RED}Failed:${NC} DFTD3 should be absent with --no-force-vdw"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} DFTD3 correctly absent with --no-force-vdw"
    PASS=$((PASS+1))
fi
rm -rf cycle_00 neb_setup.txt

echo -e "\n--- Testing --mode 1 advisory note when a force flag is toggled off (no effect) ---"
stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 5 --mode 1 --no-force-spin --ml-max-steps 30 \
    --no-intro > log_mode1_noforce.txt 2>&1
check_exit_code $? 0
check_contains "\[NOTE\].*no effect with --mode 1" log_mode1_noforce.txt
rm -rf cycle_00 image_* neb_setup.txt neb_mace_result.json


# --- 3. Composition mismatch (hard error, before any mode-dependent code runs) ---
echo -e "\n--- Testing composition mismatch (hard error) ---"
stb-neb -i initial.fdf -f final_bad_composition.fdf -c calc.fdf --mode 4 --no-intro \
    > log_bad_composition.txt 2>&1
check_exit_code $? 1
check_contains "different composition" log_bad_composition.txt


# --- 4. Lattice mismatch (hard abort -- no more silent override) ---
echo -e "\n--- Testing lattice mismatch (hard abort) ---"
rm -rf cycle_00 neb_setup.txt
stb-neb -i initial.fdf -f final_bad_lattice.fdf -c calc.fdf -n 5 --mode 4 --no-intro \
    > log_bad_lattice.txt 2>&1
check_exit_code $? 1
check_contains "\[ERROR\].*different lattices" log_bad_lattice.txt
echo "Testing: no cycle_00/ written on lattice-mismatch abort"
if [ -d cycle_00 ]; then
    echo -e "   -> ${RED}Failed:${NC} unexpected cycle_00/ written despite lattice-mismatch abort"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no cycle_00/ written -- clean abort before any image is generated"
    PASS=$((PASS+1))
fi


# --- 4a2. autosort_tol mismatch regression: a real crash found live (a raw
#     pymatgen ValueError/traceback leaking straight to the user), and then
#     found AGAIN live even after the first fix turned it into a clean
#     [ERROR] -- a real user kept hitting the same dead end through the
#     interactive menu without knowing to set the advanced --autosort-tol 0
#     escape hatch by hand. Final fix: stb-neb now retries automatically
#     with --autosort-tol 0 (index-based, no distance matching -- always
#     correct once endpoints share a guaranteed matching atom order, e.g.
#     two structures built from the same base structure) and only prints a
#     [WARNING], not a hard [ERROR] -- never a raw traceback either way. ---
echo -e "\n--- Testing the autosort_tol mismatch -> automatic fallback (no more dead end) ---"
rm -rf cycle_00 neb_setup.txt
stb-neb -i initial.fdf -f final_bad_autosort.fdf -c calc.fdf -n 5 --mode 4 --autosort-tol 0.3 \
    --no-intro > log_bad_autosort.txt 2>&1
check_exit_code $? 0
check_contains "\[WARNING\] Could not match atoms at --autosort-tol 0.3" log_bad_autosort.txt
check_contains "Falling back to using --initial/--final's own atom order" log_bad_autosort.txt
check_success cycle_00/image_00/structure.fdf
if grep -q "Traceback (most recent call last)" log_bad_autosort.txt; then
    echo -e "   -> ${RED}Failed:${NC} a raw Python traceback leaked instead of a clean [WARNING]"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no raw traceback -- clean [WARNING] + automatic recovery instead"
    PASS=$((PASS+1))
fi
rm -rf cycle_00 neb_setup.txt

echo "Testing: --autosort-tol 0 directly (no mismatch to warn about, since it's already index-based)"
stb-neb -i initial.fdf -f final_bad_autosort.fdf -c calc.fdf -n 5 --mode 4 --autosort-tol 0 \
    --no-intro > log_autosort_zero.txt 2>&1
check_exit_code $? 0
check_success cycle_00/image_00/structure.fdf
if grep -q "\[WARNING\] Could not match atoms" log_autosort_zero.txt; then
    echo -e "   -> ${RED}Failed:${NC} unexpected fallback warning when already at --autosort-tol 0"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no fallback warning needed -- already index-based"
    PASS=$((PASS+1))
fi
rm -rf cycle_00 neb_setup.txt


# --- 4b. Coordinate-wrapping bug regression: two endpoints describing the
#     SAME physical position via different fractional-wrapping conventions
#     (frac z=0.95 vs z=-0.05 in a 10 Ang cell) must interpolate to a
#     near-static path, not a spurious ~10 Ang round trip through the cell ---
echo -e "\n--- Testing coordinate-wrapping regression (same point, different wrapping) ---"
rm -rf cycle_00 neb_setup.txt
stb-neb -i wrap_initial.fdf -f wrap_final.fdf -c calc.fdf -n 5 --mode 4 --no-intro > log_wrap.txt 2>&1
check_exit_code $? 0
python3 -c "
from stb.core import structure_io
for i in range(5):
    s = structure_io.read_fdf(f'cycle_00/image_0{i}/structure.fdf')
    z = s.atoms[0][1][2]
    assert abs(z - 0.95) < 1e-6, f'image_0{i}: expected frac z ~0.95 (same wrapped point), got {z}'
print('OK')
" > log_wrap_check.txt 2>&1
check_contains "OK" log_wrap_check.txt


# --- 5. -n minimum (fails before mode matters) ---
echo -e "\n--- Testing -n minimum (n=2) ---"
stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 2 --no-intro > log_nmin.txt 2>&1
check_exit_code $? 2


# --- 6. --idpp (orthogonal to mode) ---
echo -e "\n--- Testing --idpp ---"
rm -rf cycle_00 neb_setup.txt
stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 5 --idpp --mode 4 --no-intro > log_idpp.txt 2>&1
check_exit_code $? 0
check_contains "IDPP" log_idpp.txt
check_success cycle_00/image_02/structure.fdf


# --- 6b. Endpoint-displacement diagnostic: a broad multi-atom rearrangement
#     (not just the reacting atom) without --idpp should warn and suggest
#     --idpp; a normal single-atom reaction hop (final.fdf itself) should
#     never trigger it -- real physics finding from live user data (16/33
#     atoms > 0.3 Ang, no --idpp, ~47 eV unphysical MACE barrier). ---
echo -e "\n--- Testing the endpoint-displacement diagnostic ---"
rm -rf cycle_00 neb_setup.txt
stb-neb -i initial.fdf -f final_many_moved.fdf -c calc.fdf -n 5 --mode 4 --no-intro \
    > log_displacement.txt 2>&1
check_exit_code $? 0
check_contains "Endpoint displacement: max .* mean .*, 4/9 atom(s) moved" log_displacement.txt
check_contains "\[WARNING\] 4/9 atoms moved more than 0.3 Ang" log_displacement.txt
check_contains "Consider --idpp" log_displacement.txt
rm -rf cycle_00 neb_setup.txt

echo "Testing: --idpp suppresses the WARNING (info line stays, recommendation doesn't apply anymore)"
stb-neb -i initial.fdf -f final_many_moved.fdf -c calc.fdf -n 5 --idpp --mode 4 --no-intro \
    > log_displacement_idpp.txt 2>&1
check_exit_code $? 0
check_contains "Endpoint displacement:" log_displacement_idpp.txt
if grep -q "\[WARNING\].*atoms moved more than" log_displacement_idpp.txt; then
    echo -e "   -> ${RED}Failed:${NC} displacement WARNING should be suppressed once --idpp is used"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no displacement WARNING once --idpp is used"
    PASS=$((PASS+1))
fi
rm -rf cycle_00 neb_setup.txt

echo "Testing: a normal single-atom reaction hop (final.fdf, 1/9 atoms) never triggers this WARNING"
stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 5 --mode 4 --no-intro > log_displacement_normal.txt 2>&1
check_exit_code $? 0
check_contains "Endpoint displacement:" log_displacement_normal.txt
if grep -q "\[WARNING\].*atoms moved more than" log_displacement_normal.txt; then
    echo -e "   -> ${RED}Failed:${NC} unexpected displacement WARNING for a normal single-atom reaction hop"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no false-positive WARNING for a normal reaction hop"
    PASS=$((PASS+1))
fi
rm -rf cycle_00 neb_setup.txt


# --- 6c. Manifest-proven atom identity: the real fix for the live-reported
#     "estrutura completamente bagunçada" bug. read_relaxed_or_input's
#     --autosort-tol 0 fallback only trusts --initial/--final's own on-disk
#     atom order -- it was never itself PROOF that order is correct. Any
#     tool that builds two atom-index-corresponding structures can write a
#     neb_manifest.json into each (species_sequence, proven identical by
#     construction), and stb-neb uses it -- when present on BOTH sides --
#     to skip distance-based matching entirely instead of assuming/
#     guessing. No in-repo tool currently produces one (the symmetry-site
#     -enumeration tool that used to, stb-nebSites, was removed since it
#     over-restricted NEB endpoints), so this fabricates a pair by hand to
#     keep exercising resolve_manifest_pair/validate_manifest_pair. ---
echo -e "\n--- Testing manifest-proven atom correspondence (hand-built pair) ---"

build_manifest_pair() {
    # (Re)builds nebsites_run/site_A, nebsites_run/site_B (each carrying a
    # fresh neb_manifest.json) from initial.fdf/final.fdf (already sharing
    # composition/lattice) plus a synthetic "relaxed" siesta.XV in each
    # (same species/order as structure.fdf, only the adsorbate height
    # nudged) -- same sisl-write recipe as
    # ../../8-adsorption/bsse/test.sh's own fabricated-.XV fixture, so no
    # real SIESTA run is needed to exercise the manifest-consuming path.
    rm -rf nebsites_run
    mkdir -p nebsites_run/site_A nebsites_run/site_B
    cp initial.fdf nebsites_run/site_A/structure.fdf
    cp final.fdf nebsites_run/site_B/structure.fdf
    python3 -c "
import sisl
from stb.core import structure_io, neb_manifest
for label, dz in (('site_A', -0.1), ('site_B', -0.2)):
    fdf = structure_io.read_fdf(f'nebsites_run/{label}/structure.fdf')
    species_sequence = [sym for sym, _ in fdf.atoms]
    neb_manifest.write_manifest(f'nebsites_run/{label}', pair_id='synthetic-test-pair',
                                 site_label=label, species_sequence=species_sequence,
                                 n_substrate=len(species_sequence) - 1, n_adsorbate=1,
                                 adsorbate='H')
    pmg = structure_io.to_pymatgen(fdf)
    cart = pmg.cart_coords.copy()
    cart[-1, 2] += dz
    atoms = [sisl.Atom(str(s.specie)) for s in pmg]
    geom = sisl.Geometry(cart, atoms=atoms, lattice=sisl.Lattice(pmg.lattice.matrix))
    sisl.get_sile(f'nebsites_run/{label}/siesta.XV', mode='w').write_geometry(geom)
"
}

build_manifest_pair
check_success nebsites_run/site_A/neb_manifest.json
check_success nebsites_run/site_B/neb_manifest.json
check_success nebsites_run/site_A/siesta.XV
check_success nebsites_run/site_B/siesta.XV

rm -rf cycle_00 neb_setup.txt
stb-neb --initial nebsites_run/site_A --final nebsites_run/site_B -c calc.fdf -n 5 --mode 4 \
    --no-intro > log_manifest_proven.txt 2>&1
check_exit_code $? 0
check_contains "Atom correspondence: PROVEN via neb_manifest.json" log_manifest_proven.txt
check_contains "\[OK\] Atom order PROVEN via neb_manifest.json" log_manifest_proven.txt
check_success cycle_00/image_00/structure.fdf
if grep -q "Could not match atoms" log_manifest_proven.txt; then
    echo -e "   -> ${RED}Failed:${NC} distance-based matching text should never appear when the "
    echo "      manifest proves the correspondence"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} distance-based matching never attempted (manifest short-circuits it)"
    PASS=$((PASS+1))
fi

echo "Testing: no manifest on EITHER side -> falls back to today's distance-based behavior"
rm -rf cycle_00 neb_setup.txt
rm -f nebsites_run/site_A/neb_manifest.json nebsites_run/site_B/neb_manifest.json
stb-neb --initial nebsites_run/site_A --final nebsites_run/site_B -c calc.fdf -n 5 --mode 4 \
    --no-intro > log_no_manifest.txt 2>&1
check_exit_code $? 0
check_contains "Atom correspondence: distance-based matching" log_no_manifest.txt
if grep -q "PROVEN" log_no_manifest.txt; then
    echo -e "   -> ${RED}Failed:${NC} unexpected PROVEN messaging with no manifest present"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} clean fallback, no manifest messaging (zero regression)"
    PASS=$((PASS+1))
fi

echo "Testing: manifest present on only ONE side also falls back cleanly (not an error)"
build_manifest_pair
rm -f nebsites_run/site_B/neb_manifest.json
rm -rf cycle_00 neb_setup.txt
stb-neb --initial nebsites_run/site_A --final nebsites_run/site_B -c calc.fdf -n 5 --mode 4 \
    --no-intro > log_one_manifest.txt 2>&1
check_exit_code $? 0
check_contains "Atom correspondence: distance-based matching" log_one_manifest.txt

echo "Testing: manifest's species_sequence no longer matches its own folder's actual read-back"
echo "  structure (a stale/hand-edited manifest, e.g. left over after structure.fdf was"
echo "  regenerated) -> the self-consistency check catches this before ever reaching the"
echo "  cross-pair comparison"
build_manifest_pair
python3 -c "
import json
p = 'nebsites_run/site_A/neb_manifest.json'
m = json.load(open(p))
seq = m['species_sequence']
seq[0], seq[-1] = seq[-1], seq[0]
m['species_sequence'] = seq
json.dump(m, open(p, 'w'), indent=2)
"
rm -rf cycle_00 neb_setup.txt
stb-neb --initial nebsites_run/site_A --final nebsites_run/site_B -c calc.fdf -n 5 --mode 4 \
    --no-intro > log_stale_manifest.txt 2>&1
check_exit_code $? 1
check_contains "\[ERROR\].*does not match its own" log_stale_manifest.txt
if grep -q "Traceback (most recent call last)" log_stale_manifest.txt; then
    echo -e "   -> ${RED}Failed:${NC} a raw Python traceback leaked instead of a clean [ERROR]"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} clean [ERROR], no traceback"
    PASS=$((PASS+1))
fi

echo "Testing: manifest missing a required field -> clean [ERROR]"
build_manifest_pair
python3 -c "
import json
p = 'nebsites_run/site_A/neb_manifest.json'
m = json.load(open(p))
del m['pair_id']
json.dump(m, open(p, 'w'), indent=2)
"
rm -rf cycle_00 neb_setup.txt
stb-neb --initial nebsites_run/site_A --final nebsites_run/site_B -c calc.fdf -n 5 --mode 4 \
    --no-intro > log_incomplete_manifest.txt 2>&1
check_exit_code $? 1
check_contains "\[ERROR\].*missing required field" log_incomplete_manifest.txt

echo "Testing: two self-consistent manifests that genuinely disagree with EACH OTHER (both"
echo "  structure.fdf and its own manifest edited together for site_A, so it stays self"
echo "  -consistent -- isolates validate_manifest_pair's cross-check from the self"
echo "  -consistency check above; no fabricated .XV here, since read_relaxed_or_input's own"
echo "  .XV-vs-structure.fdf species check would otherwise catch this first for a different"
echo "  reason)"
build_manifest_pair
rm -f nebsites_run/site_A/siesta.XV nebsites_run/site_B/siesta.XV
python3 -c "
import json
# Edit the .fdf TEXT directly (swap the first and last coordinate LINES,
# species-index column included) rather than round-tripping through
# read_fdf/write_fdf -- write_fdf always re-groups atoms by species on
# write (see its own docstring/CLAUDE.md), which would silently UNDO an
# in-memory atoms-list reorder between two atoms of different species,
# defeating the point of this test.
path = 'nebsites_run/site_A/structure.fdf'
with open(path) as f:
    lines = f.readlines()
start = next(i for i, l in enumerate(lines) if 'AtomicCoordinatesAndAtomicSpecies' in l and 'block' in l and 'end' not in l) + 1
end = next(i for i, l in enumerate(lines) if 'endblock AtomicCoordinatesAndAtomicSpecies' in l)
coord_lines = lines[start:end]
coord_lines[0], coord_lines[-1] = coord_lines[-1], coord_lines[0]
lines[start:end] = coord_lines
with open(path, 'w') as f:
    f.writelines(lines)

p = 'nebsites_run/site_A/neb_manifest.json'
m = json.load(open(p))
seq = m['species_sequence']
seq[0], seq[-1] = seq[-1], seq[0]
m['species_sequence'] = seq
json.dump(m, open(p, 'w'), indent=2)
"
rm -rf cycle_00 neb_setup.txt
stb-neb --initial nebsites_run/site_A --final nebsites_run/site_B -c calc.fdf -n 5 --mode 4 \
    --no-intro > log_cross_pair_mismatch.txt 2>&1
check_exit_code $? 1
check_contains "\[ERROR\].*manifests disagree on species_sequence" log_cross_pair_mismatch.txt
if grep -q "Traceback (most recent call last)" log_cross_pair_mismatch.txt; then
    echo -e "   -> ${RED}Failed:${NC} a raw Python traceback leaked instead of a clean [ERROR]"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} clean [ERROR], no traceback"
    PASS=$((PASS+1))
fi
rm -rf nebsites_run cycle_00 neb_setup.txt


echo -e "\n--- Testing read_relaxed_or_input's per-index species hardening directly (no CLI) ---"
mkdir -p xv_species_swap_test
cp initial.fdf xv_species_swap_test/structure.fdf
python3 -c "
import sisl
import sys
from stb.core import structure_io
fdf = structure_io.read_fdf('xv_species_swap_test/structure.fdf')
pmg = structure_io.to_pymatgen(fdf)
cart = pmg.cart_coords.copy()
symbols = [str(s.specie) for s in pmg]
# swap species labels at indices 0 (C) and 8 (H) -- a real .XV would never
# do this (SIESTA never reorders atoms), but this is exactly the failure
# mode the count-only check (already existed before this fix) can't see:
# per-index species now disagree, while the total per-species COUNT is
# unchanged (still 8 C + 1 H overall).
symbols[0], symbols[8] = symbols[8], symbols[0]
atoms = [sisl.Atom(sym) for sym in symbols]
geom = sisl.Geometry(cart, atoms=atoms, lattice=sisl.Lattice(pmg.lattice.matrix))
sisl.get_sile('xv_species_swap_test/siesta.XV', mode='w').write_geometry(geom)
try:
    structure_io.read_relaxed_or_input('xv_species_swap_test')
    print('NO ERROR RAISED -- BUG')
except ValueError as e:
    assert 'disagree on atom species at index' in str(e), f'unexpected message: {e}'
    print('OK')
"  > log_xv_species_check.txt 2>&1
check_contains "OK" log_xv_species_check.txt
rm -rf xv_species_swap_test


# --- 7. Mode 2 (MACE-MP-0, cached locally, then 1 single-point per image) ---
echo -e "\n--- Testing --mode 2 ---"
rm -rf cycle_00 image_* neb_setup.txt
stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 5 --mode 2 --ml-max-steps 30 --no-intro \
    > log_mode2.txt 2>&1
check_exit_code $? 0
check_contains "MACE-MP-0" log_mode2.txt
check_contains "barrier" log_mode2.txt
check_success neb_ml_preview.png
check_success neb_path.xyz
check_success image_00/structure.fdf
check_success image_04/structure.fdf
check_contains "# MODE: 2" neb_setup.txt
check_contains "# ML_NEB_USED: yes" neb_setup.txt
check_contains "ML-NEB freeze substrate: yes" neb_setup.txt
check_contains "Freezing 8/9 atom" log_mode2.txt
echo "Testing: mode 2 writes directly under the output root, NOT under cycle_00/ (single-shot, no refinement loop)"
if [ -d cycle_00 ]; then
    echo -e "   -> ${RED}Failed:${NC} unexpected cycle_00/ written for --mode 2"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no cycle_00/ for --mode 2 (image_NN/ directly)"
    PASS=$((PASS+1))
fi


# --- 7b. --no-ml-freeze-substrate (freeze-substrate is default ON) ---
echo -e "\n--- Testing --no-ml-freeze-substrate ---"
rm -rf image_* neb_setup.txt neb_ml_preview.png neb_path.xyz
stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 5 --mode 2 --ml-max-steps 10 \
    --no-ml-freeze-substrate --no-intro > log_nofreeze.txt 2>&1
check_exit_code $? 0
check_contains "ML-NEB freeze substrate: no" log_nofreeze.txt
if grep -q "Freezing" log_nofreeze.txt; then
    echo -e "   -> ${RED}Failed:${NC} unexpected 'Freezing' message with --no-ml-freeze-substrate"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no atoms frozen with --no-ml-freeze-substrate"
    PASS=$((PASS+1))
fi


# --- 7c. Path-quality warning (deliberately dense -n triggers "nearly identical") ---
echo -e "\n--- Testing path-quality warning (dense -n) ---"
rm -rf cycle_00 image_* neb_setup.txt neb_ml_preview.png neb_path.xyz
stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 60 --mode 4 --no-intro > log_dense.txt 2>&1
check_exit_code $? 0
check_contains "WARNING.*nearly identical" log_dense.txt


# --- 8. --ml-prerelax-endpoints, independent of --mode ---
echo -e "\n--- Testing --ml-prerelax-endpoints with --mode 4 (no climbing-image NEB) ---"
rm -rf cycle_00 image_* neb_setup.txt neb_ml_preview.png
stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 5 --mode 4 --ml-prerelax-endpoints \
    --no-intro > log_prerelax_endpoints.txt 2>&1
check_exit_code $? 0
check_contains "ML pre-relax" log_prerelax_endpoints.txt
check_contains "MACE-MP-0 path shaping: no" log_prerelax_endpoints.txt
if grep -q "running a real climbing-image NEB" log_prerelax_endpoints.txt; then
    echo -e "   -> ${RED}Failed:${NC} unexpected climbing-image NEB run with --mode 4"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no climbing-image NEB run with --mode 4 (--ml-prerelax-endpoints is independent)"
    PASS=$((PASS+1))
fi


# --- 9. Inert ML/IDPP flags without their master mode/flag ---
echo -e "\n--- Testing inert ML flags with --mode 4 (no MACE stage to use them) ---"
rm -rf cycle_00 neb_setup.txt
stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 5 --mode 4 --ml-k 0.2 --ml-max-steps 50 \
    --idpp-fmax 0.2 --no-intro > log_inert_flags.txt 2>&1
check_exit_code $? 0


# --- 10. Mode 1 (100% MACE, JSON only, no SIESTA folders) ---
echo -e "\n--- Testing --mode 1 ---"
rm -rf cycle_00 image_* neb_setup.txt neb_mace_result.json
stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 5 --mode 1 --ml-max-steps 30 --no-intro \
    > log_mode1.txt 2>&1
check_exit_code $? 0
check_contains "\[2\] MACE RESULT (JSON)" log_mode1.txt
check_contains "# MODE: 1" neb_setup.txt
check_success neb_mace_result.json
python3 -c "
import json
with open('neb_mace_result.json') as f:
    d = json.load(f)
assert d['mode'] == 1
assert len(d['images']) == 5
assert 'barrier_forward_eV' in d and 'barrier_backward_eV' in d and 'reaction_energy_eV' in d
assert all({'index', 'label', 'reaction_coord', 'energy_eV', 'symbols', 'positions', 'cell'} <= set(img) for img in d['images'])
print('OK')
" > log_mode1_json_check.txt 2>&1
check_contains "OK" log_mode1_json_check.txt
echo "Testing: mode 1 writes no SIESTA folders at all"
if [ -d image_00 ] || [ -d cycle_00 ]; then
    echo -e "   -> ${RED}Failed:${NC} unexpected SIESTA folder(s) written for --mode 1"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} no SIESTA folders for --mode 1 (JSON only)"
    PASS=$((PASS+1))
fi


# --- 11. Mode 3 (MACE + cycle_00/, --climb-after 0 suggested) ---
echo -e "\n--- Testing --mode 3 (default) ---"
rm -rf cycle_00 neb_setup.txt neb_ml_preview.png neb_path.xyz
stb-neb -i initial.fdf -f final.fdf -c calc.fdf -n 5 --ml-max-steps 30 --no-intro \
    > log_mode3.txt 2>&1
check_exit_code $? 0
check_contains "# MODE: 3" neb_setup.txt
check_success cycle_00/image_00/structure.fdf
check_contains "\[4\] CLUSTER SUBMISSION" log_mode3.txt
check_contains "stb-nebCycle --dir . --fmax 0.05 --climb-after 0" log_mode3.txt
echo "Testing: --mode is not required (default is 3)"
check_exit_code $? 0


# --- 12. Missing files ---
echo -e "\n--- Testing missing files ---"
stb-neb -i does_not_exist.fdf -f final.fdf -c calc.fdf --no-intro > log_missing.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing.txt


# --- 13. --version / --help ---
echo -e "\n--- Testing --version / --help ---"
stb-neb --version > log_version.txt 2>&1
check_contains "stb-neb" log_version.txt
stb-neb --help > log_help.txt 2>&1
check_contains "n-images" log_help.txt
check_contains "mode" log_help.txt
check_contains "ml-freeze-substrate" log_help.txt
check_contains "ml-freeze-threshold" log_help.txt
check_contains "force-spin" log_help.txt
check_contains "force-vdw" log_help.txt
check_contains "force-dipole" log_help.txt


# --- 14. Interactive path (stb-suite, shortcut 4.9.1) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 4.9.1) ---"

echo "Testing: navigate 4.9.1 -> defaults (mode 3) -> quit"
rm -rf cycle_00 image_* neb_setup.txt neb_ml_preview.png neb_path.xyz
{
  echo "4.9.1"
  echo "initial.fdf"   # initial_file
  echo "final.fdf"     # final_file
  echo "calc.fdf"      # calc_file
  echo ""              # pp_path (skip)
  echo ""              # force_spin (default Y)
  echo ""              # force_vdw (default Y)
  echo ""              # force_dipole (default Y)
  echo "5"             # n_images
  echo ""              # idpp_choice (default N)
  echo ""              # mode_choice (default -> 3)
  echo "0.1"            # ml_k (mode 3 needs MACE tuning prompts)
  echo "30"             # ml_max_steps (small, keep the test fast)
  echo ""                # freeze_choice (default Y)
  echo ""                # ml_freeze_threshold (default 0.3)
  echo ""              # ml_prerelax_choice (default N)
  echo ""              # show_advanced (default -> skip)
  echo ""              # press enter to continue
  echo "0"             # quit stage submenu
} | stb-suite > log_menu.txt 2>&1
check_contains "Success:.*5 image folder" log_menu.txt
check_success cycle_00/image_00/structure.fdf
check_success cycle_00/image_04/structure.fdf

echo "Testing: navigate 4.9.1 -> --idpp on, mode 2 -> quit"
rm -rf cycle_00 image_* neb_setup.txt neb_ml_preview.png neb_path.xyz
{
  echo "4.9.1"
  echo "initial.fdf"
  echo "final.fdf"
  echo "calc.fdf"
  echo ""
  echo ""              # force_spin (default Y)
  echo ""              # force_vdw (default Y)
  echo ""              # force_dipole (default Y)
  echo "5"
  echo "y"              # idpp_choice -> Y
  echo "2"               # mode_choice -> 2
  echo "0.1"             # ml_k
  echo "30"              # ml_max_steps
  echo ""                # freeze_choice (default Y)
  echo ""                # ml_freeze_threshold
  echo ""               # ml_prerelax_choice
  echo ""               # show_advanced
  echo ""               # press enter
  echo "0"
} | stb-suite > log_menu_idpp.txt 2>&1
check_contains "Success:.*5 image folder" log_menu_idpp.txt
check_success image_02/structure.fdf

echo "Testing: navigate 4.9.1 -> mode 1 (JSON only) -> quit"
rm -rf cycle_00 image_* neb_setup.txt neb_ml_preview.png neb_path.xyz neb_mace_result.json
{
  echo "4.9.1"
  echo "initial.fdf"
  echo "final.fdf"
  echo "calc.fdf"
  echo ""
  echo ""              # force_spin (default Y)
  echo ""              # force_vdw (default Y)
  echo ""              # force_dipole (default Y)
  echo "5"
  echo ""               # idpp_choice
  echo "1"               # mode_choice -> 1
  echo "0.1"             # ml_k
  echo "30"              # ml_max_steps
  echo ""                # freeze_choice
  echo ""                # ml_freeze_threshold
  echo ""               # ml_prerelax_choice
  echo ""               # show_advanced
  echo ""               # press enter
  echo "0"
} | stb-suite > log_menu_mode1.txt 2>&1
check_success neb_mace_result.json
check_success neb_ml_preview.png
check_success neb_path.xyz
check_contains "MACE freeze substrate.*ON" log_menu_mode1.txt


popd > /dev/null

# --- 15. Summary ---
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
