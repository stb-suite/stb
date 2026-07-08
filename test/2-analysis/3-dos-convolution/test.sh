#!/bin/bash

# --- Setup ---
# Smoke test for stb-convdos (DOS Convolution, item 2.3)
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$FIXTURE_DIR/test_files"

# Non-interactive matplotlib backend: a run without --no-plot still calls
# plt.show(), which would otherwise block waiting for a window to close.
export MPLBACKEND=Agg

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

# Checks that file $1 has exactly $2 lines
check_line_count() {
    local actual
    actual=$(wc -l < "$1")
    if [ "$actual" -eq "$2" ]; then
        echo -e "   -> ${GREEN}Verified:${NC} '$1' has $actual lines (expected $2)"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} '$1' has $actual lines (expected $2)"
        FAIL=$((FAIL+1))
    fi
}


# --- 1. Preparation ---
echo "--- Starting tester for STB-ConvDOS (item 2.3) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
cp "$FIXTURE_DIR/dos_total.dat" "$TEST_DIR/"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. Real fixture: default run (6000 points x 4 DOS columns from stb-dos,
#     energy grid spacing 0.01 eV) -- --sigma is in meV and --size is
#     auto-computed (not given here) ---
echo -e "\n--- Testing default run (real fixture, --sigma in meV, auto --size, --no-plot) ---"
rm -f filtered.dat
stb-convdos --file dos_total.dat --sigma 50 --out filtered.dat --no-plot --no-intro > log_default.txt 2>&1
check_exit_code $? 0
check_success filtered.dat
check_line_count filtered.dat 6002

echo "Verifying --sigma (meV) is converted using the file's own energy spacing, and --size is auto-sized"
check_contains "Energy grid spacing: 0.010000 eV -> sigma = 50.000 meV = 5.000 samples, kernel size = 31" log_default.txt

echo "Verifying the header lists every actual DOS column (used to be hardcoded 'Energy DOS_filtered' regardless of column count)"
check_contains "# Energy(eV) s_filtered p_filtered d_filtered f_filtered" filtered.dat

echo "Verifying the header's 2nd line records the broadening parameters actually used"
check_contains "# Gaussian broadening: sigma = 50.000 meV (5.000 samples), kernel size = 31" filtered.dat

echo "Verifying a known filtered data point (sigma=50 meV -> 5 samples, auto size=31)"
check_contains "13.188850 0.010066 0.788188 1.782186 0.000000" filtered.dat


# --- 2b. Explicit --size overrides the auto-computed default ---
echo -e "\n--- Testing an explicit --size override ---"
rm -f filtered_explicit.dat
stb-convdos --file dos_total.dat --sigma 50 --size 41 --out filtered_explicit.dat --no-plot --no-intro > log_explicit_size.txt 2>&1
check_exit_code $? 0
check_contains "kernel size = 41" log_explicit_size.txt


# --- 2c. --fwhm is an alternative to --sigma (FWHM = 2.3548 * sigma) ---
echo -e "\n--- Testing --fwhm as an alternative to --sigma ---"
rm -f filtered_fwhm.dat
stb-convdos --file dos_total.dat --fwhm 117.741 --out filtered_fwhm.dat --no-plot --no-intro > log_fwhm.txt 2>&1
check_exit_code $? 0
check_contains "sigma = 50.000 meV" log_fwhm.txt

echo "Testing: --sigma and --fwhm together are rejected"
stb-convdos --file dos_total.dat --sigma 50 --fwhm 100 --out f.dat --no-intro > log_both.txt 2>&1
check_exit_code $? 2

echo "Testing: neither --sigma nor --fwhm is rejected"
stb-convdos --file dos_total.dat --out f.dat --no-intro > log_neither.txt 2>&1
check_exit_code $? 2


# --- 2d. Kernel wider than the data errors out cleanly instead of crashing
#     with a raw traceback -- np.convolve(mode='same') returns a
#     max(M,N)-length array, not one matching the data, once the kernel
#     outgrows it, which used to reach the array-assignment line and raise
#     an uncaught ValueError. ---
echo -e "\n--- Testing the oversized-kernel error ---"
stb-convdos --file dos_total.dat --sigma 30000 --out f.dat --no-plot --no-intro > log_oversized.txt 2>&1
check_exit_code $? 1
check_contains "larger than the number of energy points" log_oversized.txt
check_not_contains "Traceback" log_oversized.txt


# --- 3. Kernel/energy-conservation sanity: total DOS is approximately
#     preserved by a normalized kernel (allowing small edge losses from
#     zero-padding where the real DOS is already ~0) ---
echo -e "\n--- Testing that Gaussian broadening approximately conserves total DOS ---"
python3 - << 'EOF'
import numpy as np
orig = np.loadtxt("dos_total.dat")
filt = np.loadtxt("filtered.dat")
ok = True
for col, name in zip(range(1, 5), ["s", "p", "d", "f"]):
    o, f = orig[:, col].sum(), filt[:, col].sum()
    if o > 0 and abs(f - o) / o > 0.01:
        print(f"FAIL: {name} sum drifted more than 1%: orig={o:.4f} filt={f:.4f}")
        ok = False
raise SystemExit(0 if ok else 1)
EOF
check_exit_code $? 0


# --- 3b. Many-column input (e.g. stb-dos --projection ml output) prints a
#     heads-up about the combined figure being tall instead of silently
#     opening/blocking on one window per column like it used to ---
echo -e "\n--- Testing the many-columns plotting heads-up ---"
python3 - << 'EOF'
import numpy as np
ncols = 16
energy = np.linspace(-10, 10, 200)
data = np.column_stack([energy, np.random.rand(200, ncols)])
header = "Energy(eV) " + " ".join(f"c{i}" for i in range(ncols))
np.savetxt("many_cols.dat", data, header=header)
EOF
stb-convdos --file many_cols.dat --sigma 0.5 --out many_cols_filtered.dat --no-intro > log_many_cols.txt 2>&1
check_exit_code $? 0
check_contains "16 columns to plot" log_many_cols.txt


# --- 4. Plotting path: a single figure (one plt.show() call), not one per
#     column -- used to require closing N blocking windows in a row. Under
#     MPLBACKEND=Agg this can't block, but it must still not crash and must
#     produce identical numeric output to the --no-plot run. ---
echo -e "\n--- Testing the plotting path (single combined figure) ---"
rm -f filtered_plot.dat
stb-convdos --file dos_total.dat --sigma 50 --out filtered_plot.dat --no-intro > log_plot.txt 2>&1
check_exit_code $? 0
diff -q filtered.dat filtered_plot.dat > /dev/null 2>&1
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verified:${NC} plotting does not change the numeric output"
    PASS=$((PASS+1))
else
    echo -e "   -> ${RED}Failed:${NC} filtered_plot.dat differs from filtered.dat"
    FAIL=$((FAIL+1))
fi


# --- 5. Validation: --size/--sigma reject values that used to silently
#     produce garbage (even size -> half-bin shift; size<=0 -> all-zero
#     output; sigma<=0 -> all-NaN output) ---
echo -e "\n--- Testing --size/--sigma validation ---"

echo "Testing: even --size is rejected"
stb-convdos --file dos_total.dat --sigma 50 --size 10 --out f.dat --no-intro > log_even.txt 2>&1
check_exit_code $? 2
check_contains "must be a positive odd number" log_even.txt

echo "Testing: --size 0 is rejected"
stb-convdos --file dos_total.dat --sigma 50 --size 0 --out f.dat --no-intro > log_size0.txt 2>&1
check_exit_code $? 2

echo "Testing: negative --size is rejected"
stb-convdos --file dos_total.dat --sigma 50 --size -5 --out f.dat --no-intro > log_sizeneg.txt 2>&1
check_exit_code $? 2

echo "Testing: --sigma 0 is rejected"
stb-convdos --file dos_total.dat --sigma 0 --out f.dat --no-intro > log_sigma0.txt 2>&1
check_exit_code $? 2
check_contains "must be positive" log_sigma0.txt

echo "Testing: negative --sigma is rejected"
stb-convdos --file dos_total.dat --sigma -1.0 --out f.dat --no-intro > log_sigmaneg.txt 2>&1
check_exit_code $? 2


# --- 6. Error and robustness cases ---
echo -e "\n--- Testing error cases ---"

echo "Testing: nonexistent input file"
stb-convdos --file does_not_exist.dat --sigma 50 --out f.dat --no-intro > log_missing.txt 2>&1
check_exit_code $? 1

echo "Testing: single-column input file (no DOS columns to filter)"
printf "1.0\n2.0\n3.0\n" > single_col.dat
stb-convdos --file single_col.dat --sigma 50 --out f.dat --no-intro > log_single_col.txt 2>&1
check_exit_code $? 1
check_contains "at least 2 columns" log_single_col.txt

echo "Testing: non-increasing energy column (can't determine grid spacing)"
printf "1.0 1.0\n1.0 2.0\n0.5 3.0\n" > unsorted.dat
stb-convdos --file unsorted.dat --sigma 50 --out f.dat --no-intro > log_unsorted.txt 2>&1
check_exit_code $? 1
check_contains "strictly increasing order" log_unsorted.txt

echo "Testing: missing required arguments"
stb-convdos --file dos_total.dat --no-intro > log_missing_args.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-convdos --version > log_version.txt 2>&1
check_contains "stb-convdos" log_version.txt


# --- 7. Interactive path (stb-suite, shortcut 2.3) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 2.3) ---"

echo "Testing: navigate 2.3 -> dos_total.dat -> default output -> sigma 50 meV -> auto size -> no plot"
rm -f dos_filtered.dat
printf '2.3\ndos_total.dat\n\n\n50\n\nn\n' | stb-suite > log_menu.txt 2>&1
check_success dos_filtered.dat
check_contains "# Energy(eV) s_filtered p_filtered d_filtered f_filtered" dos_filtered.dat

echo "Testing: navigate 2.3 -> dos_total.dat -> default output -> FWHM 117.741 meV -> auto size -> no plot"
rm -f dos_filtered.dat
printf '2.3\ndos_total.dat\n\n2\n117.741\n\nn\n' | stb-suite > log_menu_fwhm.txt 2>&1
check_success dos_filtered.dat
check_contains "sigma = 50.000 meV" log_menu_fwhm.txt


popd > /dev/null

# --- 8. Summary ---
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
