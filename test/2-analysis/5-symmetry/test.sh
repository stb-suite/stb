#!/bin/bash

# --- Setup ---
# Smoke test for stb-symmetry (Crystal Symmetry Analyzer, item 2.5)
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
    if grep -q -- "$1" "$2" 2>/dev/null; then
        echo -e "   -> ${GREEN}Verified:${NC} '$1' found in '$2'"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} '$1' NOT found in '$2'"
        FAIL=$((FAIL+1))
    fi
}

# Checks that $2 does NOT contain (grep -q) pattern $1
check_not_contains() {
    if grep -q -- "$1" "$2" 2>/dev/null; then
        echo -e "   -> ${RED}Failed:${NC} '$1' found in '$2' (should not be there)"
        FAIL=$((FAIL+1))
    else
        echo -e "   -> ${GREEN}Verified:${NC} '$1' absent from '$2'"
        PASS=$((PASS+1))
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
echo "--- Starting tester for STB-Symmetry (item 2.5) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/structure.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/siesta.STRUCT_OUT" "$TEST_DIR/"
cp "$FIXTURE_DIR/nacl.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/nacl_noisy.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/rhombo.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/graphene_slab.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/graphene_distorted.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/molecule.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/wire.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/single_atom.fdf" "$TEST_DIR/"
cp "$FIXTURE_DIR/water_distorted.fdf" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. --format fdf on a P1 (no symmetry) structure ---
echo -e "\n--- Testing --format fdf on the P1 Sn3O4 structure ---"
rm -f symmetry.dat
stb-symmetry --file structure.fdf --format fdf --no-intro > log_p1.txt 2>&1
check_exit_code $? 0
check_success symmetry.dat

echo "Verifying the report has a title block and generation timestamp"
check_contains "CRYSTAL SYMMETRY REPORT - STB Suite" symmetry.dat
check_contains "Generated        :" symmetry.dat
check_contains "Source file      : structure.fdf  (format: fdf)" symmetry.dat

echo "Verifying space group / point group / crystal system for a P1 (no symmetry) structure"
check_contains "Space group      : P1 (No. 1)" symmetry.dat
check_contains "Hall symbol      : P 1 (Hall No. 1)" symmetry.dat
check_contains "Point group      : 1" symmetry.dat
check_contains "Crystal system   : triclinic" symmetry.dat
check_contains "Pearson symbol   : aP14" symmetry.dat
check_contains "Symmetry operations: 1" symmetry.dat

echo "Verifying site symmetry column shows trivial '1' symmetry for a P1 structure"
check_contains "1a        Sn        1         1               1" symmetry.dat

echo "Verifying lattice parameters and volume"
check_contains "a = 4.883 Å   b = 5.907 Å   c = 8.238 Å" symmetry.dat
check_contains "Volume = 236.821 Å³" symmetry.dat

echo "Verifying every atom is its own Wyckoff orbit (14 distinct sites for 14 atoms, no symmetry)"
check_contains "SYMMETRICALLY DISTINCT SITES: 14" symmetry.dat
check_contains "   1  Sn   1a  " symmetry.dat


# --- 3. --format struct_out -- same physical structure, same symmetry result ---
echo -e "\n--- Testing --format struct_out (same P1 structure, different source format) ---"
rm -f symmetry.dat
stb-symmetry --file siesta.STRUCT_OUT --format struct_out --no-intro > log_struct_out.txt 2>&1
check_exit_code $? 0
check_contains "Space group      : P1 (No. 1)" symmetry.dat
check_contains "Volume = 236.821 Å³" symmetry.dat


# --- 4. A real, high-symmetry structure: NaCl rock salt (Fm-3m, No. 225) --
#     exercises Wyckoff multiplicity (4a/4b) and a large symmetry-operation
#     count, neither of which the P1 fixture can. ---
echo -e "\n--- Testing on NaCl rock salt (space group Fm-3m, No. 225) ---"
rm -f symmetry.dat
stb-symmetry --file nacl.fdf --format fdf --no-intro > log_nacl.txt 2>&1
check_exit_code $? 0

echo "Verifying space group / Pearson symbol / lattice type for a textbook Fm-3m structure"
check_contains "Space group      : Fm-3m (No. 225)" symmetry.dat
check_contains "Hall symbol      : -F 4 2 3 (Hall No. 523)" symmetry.dat
check_contains "Point group      : m-3m" symmetry.dat
check_contains "Crystal system   : cubic" symmetry.dat
check_contains "Lattice type     : cubic" symmetry.dat
check_contains "Pearson symbol   : cF8" symmetry.dat
check_contains "Symmetry operations: 192" symmetry.dat
check_not_contains "Setting choice" symmetry.dat

echo "Verifying the reduced formula and Z (formula units in this cell)"
check_contains "Formula          : NaCl (Z = 4 formula units in this cell)" symmetry.dat

echo "Verifying no vacuum warning for a genuinely 3D bulk structure"
check_not_contains "WARNING: vacuum gap" symmetry.dat

echo "Verifying the two Wyckoff orbits (4a for Na, 4b for Cl), not 8 separate 1-atom sites, with site symmetry m-3m"
check_contains "SYMMETRICALLY DISTINCT SITES: 2" symmetry.dat
check_contains "4a        Na        4         m-3m            1" symmetry.dat
check_contains "4b        Cl        4         m-3m            5" symmetry.dat
check_contains "   1  Na   4a  " symmetry.dat
check_contains "   5  Cl   4b  " symmetry.dat

echo "Verifying the ATOMIC SITES header's 'Orbit'/'Distortion' columns line up with the data rows"
check_contains "Atom  Sp.  Wyckoff      Frac. x     Frac. y     Frac. z  Orbit Distortion(Å)" symmetry.dat
check_contains "   1  Na   4a          0.000000    0.000000    0.000000  1     0.0000" symmetry.dat

echo "Verifying every atom has zero distortion for a perfectly symmetric NaCl (each atom's position is exactly reproduced by the detected symmetry operations)"
n_zero_distortion=$(grep -cE "  [12]     0\.0000$" symmetry.dat)
if [ "$n_zero_distortion" -eq 8 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} all 8 atoms show 0.0000 distortion"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} only $n_zero_distortion/8 atoms show 0.0000 distortion"
    FAIL=$((FAIL+1))
fi

echo "Verifying all 192 symmetry operations are listed, in compact x,y,z notation"
check_contains "SYMMETRY OPERATIONS (192), in x,y,z notation" symmetry.dat
check_contains "   1: x, y, z" symmetry.dat
check_contains "  49: x, y+1/2, z+1/2" symmetry.dat
check_contains " 192: -x+1/2, z+1/2, -y" symmetry.dat

n_op_lines=$(grep -c "^ *[0-9]*: " symmetry.dat)
if [ "$n_op_lines" -eq 192 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} exactly 192 operation lines found"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} found $n_op_lines operation lines, expected 192"
    FAIL=$((FAIL+1))
fi


# --- 4b. --no-operations: skips the operations list entirely ---
echo -e "\n--- Testing --no-operations ---"
rm -f symmetry.dat
stb-symmetry --file nacl.fdf --format fdf --no-operations --no-intro > log_noops.txt 2>&1
check_exit_code $? 0
check_not_contains "SYMMETRY OPERATIONS" symmetry.dat
check_contains "SYMMETRICALLY DISTINCT SITES: 2" symmetry.dat


# --- 4c. --scan-symprec: sweeps a range of tolerances and reports where
#     the space group changes -- verified on a deliberately noisy NaCl
#     (small random displacement, simulating post-relaxation numerical
#     drift) that P1 is detected at tight tolerances and a higher symmetry
#     appears once symprec is loose enough. ---
echo -e "\n--- Testing --scan-symprec ---"
rm -f symmetry.dat
stb-symmetry --file nacl.fdf --format fdf --scan-symprec --no-intro > log_scan.txt 2>&1
check_exit_code $? 0
check_contains "SYMPREC SENSITIVITY SCAN" symmetry.dat
check_contains "1e-05       Fm-3m (No. 225)" symmetry.dat
check_contains "No change in detected space group across the scanned tolerance range." symmetry.dat

echo "Testing --scan-symprec on a deliberately noisy NaCl (small random displacement): P1 at tight tolerance, symmetry recovered at a looser one"
rm -f symmetry.dat
stb-symmetry --file nacl_noisy.fdf --format fdf --scan-symprec --no-intro > log_scan_noisy.txt 2>&1
check_exit_code $? 0
check_contains "1e-05       P1 (No. 1)" symmetry.dat
check_contains "0.1         Cm (No. 8)" symmetry.dat
check_contains "Symmetry changes at symprec >= 0.1: P1 (No. 1) -> Cm (No. 8)." symmetry.dat


# --- 4d. Per-atom Wyckoff distortion on the noisy NaCl: nonzero, matching
#     the values independently verified by hand (a direct residual of
#     applying each detected symmetry operation and measuring distance to
#     the nearest same-species atom -- not from
#     SpacegroupAnalyzer.get_refined_structure(), which was tried first and
#     empirically returns positions unchanged). ---
echo -e "\n--- Testing per-atom Wyckoff position distortion on noisy NaCl ---"
rm -f symmetry.dat
stb-symmetry --file nacl_noisy.fdf --format fdf --symprec 0.1 --no-intro > log_distortion.txt 2>&1
check_exit_code $? 0
check_contains "0.002765    0.006573    0.002644  1     0.0000" symmetry.dat
check_contains "0.489575    0.507243    0.003571  2     0.0382" symmetry.dat
check_contains "0.004791    0.000318    0.497660  5     0.0688" symmetry.dat


# --- 4e. --write-refined: writes the symmetry-refined structure as .fdf ---
echo -e "\n--- Testing --write-refined ---"
rm -f refined.fdf
stb-symmetry --file nacl.fdf --format fdf --write-refined refined.fdf --no-intro > log_refined.txt 2>&1
check_exit_code $? 0
check_success refined.fdf
check_contains "%block LatticeVectors" refined.fdf
check_contains "NumberofAtoms      8" refined.fdf


# --- 4f. --compare-to: same physical structure in two formats -> PRESERVED;
#     a clean vs. a deliberately displaced NaCl -> CHANGED. ---
echo -e "\n--- Testing --compare-to ---"
rm -f symmetry.dat
stb-symmetry --file structure.fdf --format fdf --compare-to siesta.STRUCT_OUT --compare-format struct_out --no-intro > log_compare_same.txt 2>&1
check_exit_code $? 0
check_contains "SYMMETRY COMPARISON" symmetry.dat
check_contains "Symmetry PRESERVED between the two structures: both P1 (No. 1, 1 ops)." symmetry.dat

rm -f symmetry.dat
stb-symmetry --file nacl.fdf --format fdf --compare-to nacl_noisy.fdf --compare-format fdf --no-intro > log_compare_diff.txt 2>&1
check_exit_code $? 0
check_contains "Space group         Fm-3m (No. 225)             P1 (No. 1)" symmetry.dat
check_contains "Symmetry CHANGED between the two structures: Fm-3m (No. 225, 192 ops) -> P1 (No. 1, 1 ops) -- lost symmetry operations." symmetry.dat

echo "Testing: --compare-to without --compare-format is rejected"
stb-symmetry --file nacl.fdf --format fdf --compare-to nacl_noisy.fdf --no-intro > log_compare_noformat.txt 2>&1
check_exit_code $? 2


# --- 4g. Setting choice (e.g. hexagonal vs rhombohedral axes for a trigonal
#     space group) is reported when spglib returns one, absent otherwise
#     (already implicitly covered: no fixture so far has triggered it). ---
echo -e "\n--- Testing setting choice on a rhombohedral structure (space group R3, No. 146) ---"
rm -f symmetry.dat
stb-symmetry --file rhombo.fdf --format fdf --no-intro > log_rhombo.txt 2>&1
check_exit_code $? 0
check_contains "Space group      : R3 (No. 146)" symmetry.dat
check_contains "Setting choice   : H" symmetry.dat


# --- 4h. Vacuum-axis note + real layer-group detection (spglib
#     get_layergroup(), added in spglib 2.1.0): a graphene-like slab with
#     20 Ang vacuum along c should get both the note pointing at LAYER
#     GROUP and a correct layer-group analysis (verified against known
#     graphene crystallography: layer group 80, p6/mmm, Wyckoff 2b, site
#     symmetry -6m2, 24 operations -- note the 3D SPACE GROUP section
#     above it independently labels the same 2 atoms "2d", since it's a
#     different group classification (space group 191) with its own
#     Wyckoff lettering, not the layer group's); the fixtures above (all
#     genuinely 3D bulk) must get neither section. ---
echo -e "\n--- Testing the vacuum-axis note and layer-group detection on a graphene slab (20 Ang vacuum along c) ---"
rm -f symmetry.dat
stb-symmetry --file graphene_slab.fdf --format fdf --no-intro > log_slab.txt 2>&1
check_exit_code $? 0
check_contains "NOTE: vacuum gap (>= 10 Å) detected along axis/axes c" symmetry.dat
check_contains "LAYER GROUP for that (spglib's dedicated 2D detection, not a heuristic)." symmetry.dat

echo "Verifying the layer-group section itself: correct group, Hall symbol, point group, site symmetry"
check_contains "LAYER GROUP (aperiodic axis: c)" symmetry.dat
check_contains "Layer group      : p6/mmm (No. 80)" symmetry.dat
check_contains "Hall symbol      : -p 6 2 (Hall No. -116)" symmetry.dat
check_contains "Point group      : 6/mmm" symmetry.dat
check_contains "Symmetrically distinct sites: 1" symmetry.dat
check_contains "2b        C         2         -6m2            1" symmetry.dat

echo "Verifying the layer-group atomic-sites table has its own Wyckoff distortion column (0.0000 for perfectly symmetric graphene)"
check_contains "Atom  Sp.  Wyckoff   Site symmetry   Orbit Distortion(Å)" symmetry.dat
check_contains "   1  C    2b        -6m2            1     0.0000" symmetry.dat

echo "Verifying --scan-symprec on a slab also scans the layer group (not just the 3D space group)"
rm -f symmetry.dat
stb-symmetry --file graphene_slab.fdf --format fdf --scan-symprec --no-intro > log_slab_scan.txt 2>&1
check_exit_code $? 0
check_contains "LAYER GROUP SYMPREC SENSITIVITY SCAN" symmetry.dat
check_contains "1e-05       p6/mmm (No. 80)" symmetry.dat
check_contains "No change in detected layer group across the scanned tolerance range." symmetry.dat

echo "Verifying all 24 layer-group symmetry operations are listed, in the same compact x,y,z notation as the 3D section"
n_lg_ops=$(grep -c "^ *[0-9]*: " symmetry.dat)
if [ "$n_lg_ops" -eq 48 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} 48 total operation lines (24 layer-group + 24 space-group)"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} found $n_lg_ops operation lines, expected 48 (24+24)"
    FAIL=$((FAIL+1))
fi

echo "Verifying a genuinely 3D bulk structure gets neither the vacuum note nor a layer-group section"
rm -f symmetry.dat
stb-symmetry --file nacl.fdf --format fdf --no-intro > log_nacl_novacuum.txt 2>&1
check_not_contains "NOTE: vacuum gap" symmetry.dat
check_not_contains "LAYER GROUP" symmetry.dat

echo "Verifying --no-operations suppresses both the space-group and layer-group operation lists, but keeps the layer-group site table"
rm -f symmetry.dat
stb-symmetry --file graphene_slab.fdf --format fdf --no-operations --no-intro > log_slab_noops.txt 2>&1
check_exit_code $? 0
check_contains "LAYER GROUP (aperiodic axis: c)" symmetry.dat
check_contains "Symmetrically distinct sites: 1" symmetry.dat
check_not_contains "in x,y,z notation" symmetry.dat


# --- 4i. --compare-to accounts for layer groups too when both structures
#     have one: preserved (two runs of the same slab) and changed
#     (graphene vs. a deliberately displaced-atom variant, p6/mmm No. 80
#     -> c2/m11 No. 18). ---
echo -e "\n--- Testing --compare-to's layer-group awareness ---"
rm -f symmetry.dat
stb-symmetry --file graphene_slab.fdf --format fdf --compare-to graphene_slab.fdf --compare-format fdf --no-intro > log_compare_slab_same.txt 2>&1
check_exit_code $? 0
check_contains "Layer group PRESERVED: both p6/mmm (No. 80, 24 ops)." symmetry.dat

rm -f symmetry.dat
stb-symmetry --file graphene_slab.fdf --format fdf --compare-to graphene_distorted.fdf --compare-format fdf --no-intro > log_compare_slab_diff.txt 2>&1
check_exit_code $? 0
check_contains "Layer group         p6/mmm (No. 80)             c2/m11 (No. 18)" symmetry.dat
check_contains "Layer group ops     24                          4" symmetry.dat
check_contains "Layer group CHANGED: p6/mmm (No. 80, 24 ops) -> c2/m11 (No. 18, 4 ops) -- lost symmetry operations." symmetry.dat

echo "Verifying the LAYER GROUP section reports its own setting choice (was silently dropped -- ds.choice was never even read), distinct from SPACE GROUP's"
rm -f symmetry.dat
stb-symmetry --file graphene_distorted.fdf --format fdf --no-intro > log_distorted_choice.txt 2>&1
check_exit_code $? 0
check_contains "Setting choice   : b1" symmetry.dat
check_contains "Setting choice   : a" symmetry.dat


# --- 4j. --write-refined on a vacuum-padded structure with exactly one
#     vacuum axis now refines against the LAYER group (std_lattice/
#     std_positions from spglib.get_layergroup()) instead of the
#     non-physical 3D space group -- verified the vacuum thickness (c=20
#     Ang) survives unchanged. A structure with 2+ vacuum axes (no
#     rod/point-group refinement available) is still blocked outright,
#     same reasoning as before. ---
echo -e "\n--- Testing --write-refined on a vacuum-padded structure uses the layer group ---"
rm -f refined_slab.fdf
stb-symmetry --file graphene_slab.fdf --format fdf --write-refined refined_slab.fdf --no-intro > log_refined_slab.txt 2>&1
check_exit_code $? 0
check_contains "\[INFO\] Writing layer-group-refined structure to refined_slab.fdf" log_refined_slab.txt
check_success refined_slab.fdf
check_contains " 0.00000000   0.00000000   20.00000000" refined_slab.fdf
check_contains "NumberofAtoms      2" refined_slab.fdf

echo "Verifying --write-refined still works normally for a genuinely 3D bulk structure"
rm -f refined_bulk.fdf
stb-symmetry --file nacl.fdf --format fdf --write-refined refined_bulk.fdf --no-intro > log_refined_bulk.txt 2>&1
check_exit_code $? 0
check_contains "\[INFO\] Writing symmetry-refined structure to refined_bulk.fdf" log_refined_bulk.txt
check_not_contains "ERROR" log_refined_bulk.txt
check_success refined_bulk.fdf

echo "Verifying --write-refined works for an isolated molecule (3 vacuum axes) via the point group, molecule re-centered back in the original box"
rm -f refined_mol.fdf
stb-symmetry --file molecule.fdf --format fdf --write-refined refined_mol.fdf --no-intro > log_refined_mol.txt 2>&1
check_exit_code $? 0
check_contains "\[INFO\] Writing point-group-refined structure to refined_mol.fdf" log_refined_mol.txt
check_success refined_mol.fdf
check_contains "0.50000000   0.50000000   0.50000000   1" refined_mol.fdf

echo "Verifying --write-refined is still blocked for a wire (2 vacuum axes -- no rod-group equivalent exists)"
rm -f refined_wire.fdf
stb-symmetry --file wire.fdf --format fdf --write-refined refined_wire.fdf --no-intro > log_refined_wire.txt 2>&1
check_exit_code $? 0
check_contains "\[ERROR\] --write-refined skipped" log_refined_wire.txt
check_contains "no rod-group equivalent is available" log_refined_wire.txt
if [ -e refined_wire.fdf ]; then
    echo -e "   -> ${RED}Failed:${NC} refined_wire.fdf was written despite 2 vacuum axes"
    FAIL=$((FAIL+1))
else
    echo -e "   -> ${GREEN}Verified:${NC} refined_wire.fdf was NOT written"
    PASS=$((PASS+1))
fi


# --- 4k. POINT GROUP report section for an isolated molecule (3 vacuum
#     axes): pymatgen's PointGroupAnalyzer, a separate non-periodic
#     detector from spglib. Verified against textbook water: C2v, the two
#     H atoms correctly identified as symmetry-equivalent. A wire (2
#     vacuum axes) must get neither this section nor a working detection
#     (spglib has no rod-group equivalent to fall back on). ---
echo -e "\n--- Testing the POINT GROUP section on an isolated water molecule (point group C2v) ---"
rm -f symmetry.dat
stb-symmetry --file molecule.fdf --format fdf --no-intro > log_molecule.txt 2>&1
check_exit_code $? 0
check_contains "NOTE: vacuum gap (>= 10 Å) detected along axis/axes a, b, c" symmetry.dat
check_contains "see POINT GROUP for that" symmetry.dat
check_contains "POINT GROUP (isolated molecule -- all 3 axes vacuum-padded)" symmetry.dat
check_contains "Point group      : C2v" symmetry.dat
check_contains "Rotational symmetry number: 2" symmetry.dat
check_contains "Symmetry operations: 4" symmetry.dat
check_contains "Symmetrically distinct atoms: 2" symmetry.dat
check_contains "O         1         1" symmetry.dat
check_contains "H         2         2, 3" symmetry.dat

echo "Verifying zero distortion for perfectly symmetric water (regression check: an earlier version operated point-group symmetry ops on uncentered Cartesian coordinates and got ~box-size-magnitude nonsense instead of ~0)"
check_contains "   1  O    0.0000" symmetry.dat
check_contains "   2  H    0.0000" symmetry.dat
check_contains "   3  H    0.0000" symmetry.dat

echo "Verifying nonzero, non-trivial distortion on a deliberately asymmetric water (one H displaced, breaking the C2 axis -> detected as Cs, not C2v)"
rm -f symmetry.dat
stb-symmetry --file water_distorted.fdf --format fdf --no-intro > log_water_distorted.txt 2>&1
check_exit_code $? 0
check_contains "Point group      : Cs" symmetry.dat
check_contains "   1  O    0.0079" symmetry.dat
check_contains "   2  H    0.0628" symmetry.dat
check_contains "   3  H    0.0628" symmetry.dat

echo "Verifying a wire (2 vacuum axes) gets neither LAYER GROUP nor POINT GROUP (no rod-group equivalent exists)"
rm -f symmetry.dat
stb-symmetry --file wire.fdf --format fdf --no-intro > log_wire.txt 2>&1
check_exit_code $? 0
check_contains "NOTE: vacuum gap (>= 10 Å) detected along axis/axes a, b" symmetry.dat
check_contains "no rod-group detection" symmetry.dat
check_not_contains "LAYER GROUP" symmetry.dat
check_not_contains "POINT GROUP" symmetry.dat


# --- 4l. A single isolated atom (3 vacuum axes, the same shape of
#     structure stb-cohesive generates as its isolated-atom reference) used
#     to crash PointGroupAnalyzer with an unrelated AttributeError deep
#     inside pymatgen, uncaught -- which aborted the ENTIRE run (exit 1,
#     "Symmetry detection failed") even though the primary 3D analysis
#     was fine. Now caught and reported as "detection... failed", same as
#     any other point-group failure, with everything else unaffected. ---
echo -e "\n--- Testing a single isolated atom doesn't crash the whole run ---"
rm -f symmetry.dat
stb-symmetry --file single_atom.fdf --format fdf --no-intro > log_single_atom.txt 2>&1
check_exit_code $? 0
check_not_contains "Traceback" log_single_atom.txt
check_not_contains "AttributeError" log_single_atom.txt
check_contains "Molecular point-group detection was attempted" symmetry.dat
check_not_contains "POINT GROUP (isolated molecule" symmetry.dat
check_contains "Space group      :" symmetry.dat


# --- 4m. --compare-to also compares the point group when both structures
#     have one (same reasoning as the layer-group comparison fix). ---
echo -e "\n--- Testing --compare-to's point-group awareness ---"
rm -f symmetry.dat
stb-symmetry --file molecule.fdf --format fdf --compare-to molecule.fdf --compare-format fdf --no-intro > log_compare_mol.txt 2>&1
check_exit_code $? 0
check_contains "Point group         C2v                         C2v" symmetry.dat
check_contains "Point group ops     4                           4" symmetry.dat
check_contains "Point group PRESERVED: both C2v (4 ops)." symmetry.dat


# --- 4n. --scan-symprec on an isolated molecule scans the point group's
#     own tolerance (not the crystallographic --symprec, a different
#     scale) instead of the layer group -- verified on water: stable C2v
#     across tight tolerances, breaking down to C1 once the tolerance (1
#     Å) is comparable to the O-H bond length itself. ---
echo -e "\n--- Testing --scan-symprec's point-group awareness ---"
rm -f symmetry.dat
stb-symmetry --file molecule.fdf --format fdf --scan-symprec --no-intro > log_scan_mol.txt 2>&1
check_exit_code $? 0
check_contains "POINT GROUP TOLERANCE SENSITIVITY SCAN" symmetry.dat
check_contains "0.01          C2v" symmetry.dat
check_contains "1             C1" symmetry.dat
check_contains "Point group changes at tolerance >= 1 Å: C2v -> C1." symmetry.dat
check_not_contains "LAYER GROUP SYMPREC SENSITIVITY SCAN" symmetry.dat


# --- 5. --output-dir: writes into (and creates) a chosen directory ---
echo -e "\n--- Testing --output-dir ---"
rm -rf out_dir
stb-symmetry --file structure.fdf --format fdf --output-dir out_dir --no-intro > log_outdir.txt 2>&1
check_exit_code $? 0
check_success out_dir/symmetry.dat


# --- 6. --symprec / --angle-tolerance are configurable and validated ---
echo -e "\n--- Testing --symprec and --angle-tolerance ---"
stb-symmetry --file structure.fdf --format fdf --symprec 1e-4 --angle-tolerance 2.0 --no-intro > log_symprec.txt 2>&1
check_exit_code $? 0
check_contains "symprec=0.0001, angle_tolerance=2°" symmetry.dat

echo "Testing: --symprec 0 is rejected"
stb-symmetry --file structure.fdf --format fdf --symprec 0 --no-intro > log_symprec_zero.txt 2>&1
check_exit_code $? 2

echo "Testing: --angle-tolerance 0 is rejected"
stb-symmetry --file structure.fdf --format fdf --angle-tolerance 0 --no-intro > log_angletol_zero.txt 2>&1
check_exit_code $? 2


# --- 7. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: nonexistent input file"
stb-symmetry --file does_not_exist.fdf --format fdf --no-intro > log_missing.txt 2>&1
check_exit_code $? 1
check_contains "not found" log_missing.txt
check_not_contains "Traceback" log_missing.txt

echo "Testing: invalid --format (cif/poscar no longer accepted)"
stb-symmetry --file structure.fdf --format cif --no-intro > log_badfmt.txt 2>&1
check_exit_code $? 2

echo "Testing: missing required arguments"
stb-symmetry --file structure.fdf --no-intro > log_missing_args.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-symmetry --version > log_version.txt 2>&1
check_contains "stb-symmetry" log_version.txt


# --- 8. Interactive path (stb-suite, shortcut 2.5) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 2.5) ---"

echo "Testing: navigate 2.5 -> nacl.fdf -> format=1(fdf) -> default output dir -> scan=y -> operations=n -> compare=skip -> write-refined=skip"
rm -f symmetry.dat
printf '2.5\nnacl.fdf\n1\n\ny\nn\n\n\n' | stb-suite > log_menu.txt 2>&1
check_success symmetry.dat
check_contains "Select input file format:" log_menu.txt
check_contains "Fm-3m" symmetry.dat
check_contains "SYMPREC SENSITIVITY SCAN" symmetry.dat
check_not_contains "SYMMETRY OPERATIONS" symmetry.dat

echo "Testing: navigate 2.5 -> structure.fdf -> format=1(fdf) -> default output dir -> scan=n -> operations=y -> compare=siesta.STRUCT_OUT(format 2) -> write-refined=refined_menu.fdf"
rm -f symmetry.dat refined_menu.fdf
printf '2.5\nstructure.fdf\n1\n\nn\ny\nsiesta.STRUCT_OUT\n2\nrefined_menu.fdf\n' | stb-suite > log_menu_compare.txt 2>&1
check_success symmetry.dat
check_contains "SYMMETRY COMPARISON" symmetry.dat
check_contains "Symmetry PRESERVED between the two structures: both P1 (No. 1, 1 ops)." symmetry.dat
check_success refined_menu.fdf


popd > /dev/null

# --- 9. Summary ---
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
