#!/bin/bash

# --- Setup ---
# Smoke test for stb-fetch (Structure Fetcher, item 1.12)
#
# COD and OPTIMADE (via the twodmatpedia provider) cases run for real (no
# API key needed, live network access confirmed reachable from this
# environment). Materials Project cases are limited to argument validation
# (no network) unless PMG_MAPI_KEY is set in the environment, in which case
# one live fetch is also exercised.
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
echo "--- Starting tester for STB-Fetch (item 1.12) ---"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR"
echo "Test directory '$TEST_DIR' prepared."

pushd "$TEST_DIR" > /dev/null


# --- 2. COD: direct id (real network call) ---
echo -e "\n--- Testing COD --cod-id (magnetite, COD 1010369) ---"
rm -f mag_cod.fdf
stb-fetch --source cod --cod-id 1010369 -o mag_cod.fdf --no-intro > log_cod_id.txt 2>&1
check_contains "Formula:.*Fe3O4" log_cod_id.txt
check_contains "Atoms:.*56" log_cod_id.txt
check_contains "Space group:.*Fd-3m" log_cod_id.txt
check_contains "oxidation-state-only disorder" log_cod_id.txt
check_success mag_cod.fdf
check_contains "NumberofAtoms      56" mag_cod.fdf


# --- 2b. --unitcell (reuses stb-unitcell's own logic) ---
echo -e "\n--- Testing --unitcell primitive/conventional on the fetched magnetite ---"

echo "Testing: --unitcell primitive (56 -> 14 atoms, matches stb-unitcell/stb-defect's fixture)"
rm -f mag_prim.fdf
stb-fetch --source cod --cod-id 1010369 --unitcell primitive -o mag_prim.fdf --no-intro > log_unitcell_prim.txt 2>&1
check_contains "Atoms:.*56" log_unitcell_prim.txt
check_contains "Unit cell mode:.*primitive" log_unitcell_prim.txt
check_contains "Output atoms:.*14" log_unitcell_prim.txt
check_success mag_prim.fdf
check_contains "NumberofAtoms      14" mag_prim.fdf

echo "Testing: --unitcell conventional (56 -> 56, no-op note)"
rm -f mag_conv.fdf
stb-fetch --source cod --cod-id 1010369 --unitcell conventional -o mag_conv.fdf --no-intro > log_unitcell_conv.txt 2>&1
check_contains "Output atoms:.*56" log_unitcell_conv.txt
check_contains "already the conventional cell" log_unitcell_conv.txt
check_success mag_conv.fdf


# --- 3. COD: formula search (real network call, multi-candidate) ---
echo -e "\n--- Testing COD --formula (Fe3O4, expect many candidates) ---"
stb-fetch --source cod --formula "Fe3O4" --no-intro > log_cod_formula.txt 2>&1
check_exit_code $? 1
check_contains "Found.*candidate" log_cod_formula.txt
check_contains "rerun with --cod-id" log_cod_formula.txt

echo "Testing: --formula without spaces still resolves to COD's Hill notation"
check_contains "COD 1010369" log_cod_formula.txt


# --- 4. COD: nonexistent id ---
echo -e "\n--- Testing COD nonexistent --cod-id ---"
stb-fetch --source cod --cod-id 999999999 --no-intro > log_cod_missing.txt 2>&1
check_exit_code $? 1
check_contains "network request failed" log_cod_missing.txt


# --- 4b. OPTIMADE: twodmatpedia (real network call, no credential needed) ---
echo -e "\n--- Testing OPTIMADE --provider twodmatpedia ---"

echo "Testing: --optimade-id exact fetch (2dm-2127, MoS monolayer)"
rm -f mos.fdf
stb-fetch --source optimade --provider twodmatpedia --optimade-id 2dm-2127 \
    -o mos.fdf --no-intro > log_optimade_id.txt 2>&1
check_contains "Fetched:.*optimade 2dm-2127" log_optimade_id.txt
check_contains "Formula:.*MoS" log_optimade_id.txt
check_contains "Atoms:.*4" log_optimade_id.txt
check_success mos.fdf
check_contains "NumberofAtoms      4" mos.fdf

echo "Testing: --formula search (MoS2 matches multiple stoichiometries, incl. exact MoS2)"
stb-fetch --source optimade --provider twodmatpedia --formula "MoS2" --limit 10 --no-intro \
    > log_optimade_formula.txt 2>&1
check_exit_code $? 1
check_contains "Found.*candidate" log_optimade_formula.txt
check_contains "rerun with --optimade-id" log_optimade_formula.txt
check_contains "2dm-3150" log_optimade_formula.txt

echo "Testing: --list-providers (works without --source)"
stb-fetch --list-providers > log_list_providers.txt 2>&1
check_exit_code $? 0
check_contains "twodmatpedia" log_list_providers.txt
check_contains "OptimadeRester.aliases" log_list_providers.txt

echo "Testing: nonexistent --optimade-id"
stb-fetch --source optimade --provider twodmatpedia --optimade-id does-not-exist-xyz --no-intro \
    > log_optimade_missing.txt 2>&1
check_exit_code $? 1


# --- 5. Argument validation (no network needed) ---
echo -e "\n--- Testing argument validation ---"

echo "Testing: --material-id rejected with --source cod"
stb-fetch --source cod --material-id mp-1 --no-intro > log_bad1.txt 2>&1
check_exit_code $? 2
check_contains "only valid with --source materials-project" log_bad1.txt

echo "Testing: --cod-id rejected with --source materials-project"
stb-fetch --source materials-project --cod-id 1 --no-intro > log_bad2.txt 2>&1
check_exit_code $? 2
check_contains "only valid with --source cod" log_bad2.txt

echo "Testing: --most-stable rejected with --source cod"
stb-fetch --source cod --formula Fe3O4 --most-stable --no-intro > log_bad4.txt 2>&1
check_exit_code $? 2

echo "Testing: --api-key rejected with --source cod"
stb-fetch --source cod --cod-id 1010369 --api-key xxx --no-intro > log_bad5.txt 2>&1
check_exit_code $? 2

echo "Testing: --provider rejected with --source cod"
stb-fetch --source cod --cod-id 1010369 --provider twodmatpedia --no-intro > log_bad6.txt 2>&1
check_exit_code $? 2
check_contains "only valid with --source optimade" log_bad6.txt

echo "Testing: --optimade-id rejected with --source materials-project"
stb-fetch --source materials-project --optimade-id x --no-intro > log_bad7.txt 2>&1
check_exit_code $? 2

echo "Testing: --source optimade without --provider"
stb-fetch --source optimade --optimade-id x --no-intro > log_bad8.txt 2>&1
check_exit_code $? 2
check_contains "requires --provider" log_bad8.txt

echo "Testing: --material-id rejected with --source optimade"
stb-fetch --source optimade --provider twodmatpedia --material-id mp-1 --no-intro > log_bad9.txt 2>&1
check_exit_code $? 2

echo "Testing: --most-stable rejected with --source optimade"
stb-fetch --source optimade --provider twodmatpedia --formula MoS2 --most-stable --no-intro \
    > log_bad10.txt 2>&1
check_exit_code $? 2

echo "Testing: --material-id/--cod-id/--optimade-id/--formula mutually exclusive"
stb-fetch --source cod --cod-id 1 --formula Fe3O4 --no-intro > log_mutex.txt 2>&1
check_exit_code $? 2

echo "Testing: no selection given (one of --material-id/--cod-id/--optimade-id/--formula required)"
stb-fetch --source cod --no-intro > log_missing_selection.txt 2>&1
check_exit_code $? 2

echo "Testing: --source missing (required)"
stb-fetch --cod-id 1 --no-intro > log_missing_source.txt 2>&1
check_exit_code $? 2

echo "Testing: --version"
stb-fetch --version > log_version.txt 2>&1
check_contains "stb-fetch" log_version.txt

echo "Testing: --help documents --source/--material-id/--cod-id/--optimade-id/--formula"
stb-fetch --help > log_help.txt 2>&1
check_contains "source" log_help.txt
check_contains "material-id" log_help.txt
check_contains "cod-id" log_help.txt
check_contains "optimade-id" log_help.txt
check_contains "provider" log_help.txt
check_contains "list-providers" log_help.txt
check_contains "formula" log_help.txt
check_contains "unitcell" log_help.txt


# --- 6. Materials Project: live smoke test, only if a key is configured ---
echo -e "\n--- Testing Materials Project (live fetch, skipped without a key) ---"
if [ -n "$PMG_MAPI_KEY" ]; then
    echo "PMG_MAPI_KEY is set -- running live materials-project fetches"
    rm -f mag_mp.fdf
    stb-fetch --source materials-project --material-id mp-19306 -o mag_mp.fdf --no-intro > log_mp.txt 2>&1
    check_contains "Formula:.*Fe3O4" log_mp.txt
    check_success mag_mp.fdf

    echo "Testing: --formula --most-stable ranks over ALL matches, not just the displayed --limit"
    rm -f mp_stable.fdf mp_stable_full.fdf
    stb-fetch --source materials-project --formula Fe3O4 --most-stable --limit 5 \
        -o mp_stable.fdf --no-intro > log_mp_stable.txt 2>&1
    check_contains "showing first 5" log_mp_stable.txt
    check_contains "most-stable selected:" log_mp_stable.txt
    check_success mp_stable.fdf
    # Same search without the low --limit must pick the identical material --
    # proving --most-stable isn't silently restricted to the displayed subset.
    # (strip ANSI color codes first -- a reset code sits between "selected:"
    # and the id, so a plain grep across that boundary would never match)
    stb-fetch --source materials-project --formula Fe3O4 --most-stable --limit 100 \
        -o mp_stable_full.fdf --no-intro > log_mp_stable_full.txt 2>&1
    picked_capped=$(sed -E 's/\x1b\[[0-9;]*m//g' log_mp_stable.txt | grep -oE 'selected: mp-[a-zA-Z0-9]+')
    picked_full=$(sed -E 's/\x1b\[[0-9;]*m//g' log_mp_stable_full.txt | grep -oE 'selected: mp-[a-zA-Z0-9]+')
    if [ -n "$picked_capped" ] && [ "$picked_capped" = "$picked_full" ]; then
        echo -e "   -> ${GREEN}Verified:${NC} --limit 5 and --limit 100 picked the same material ($picked_capped)"
        PASS=$((PASS+1))
    else
        echo -e "   -> ${RED}Failed:${NC} --most-stable picked differently under --limit 5 ($picked_capped) vs --limit 100 ($picked_full)"
        FAIL=$((FAIL+1))
    fi
else
    echo -e "   -> ${YELLOW}Skipped:${NC} PMG_MAPI_KEY not set in this environment"
fi


# --- 7. Interactive path (stb-suite, shortcut 1.12) ---
echo -e "\n--- Testing the interactive path via stb-suite (shortcut 1.12) ---"

echo "Testing: navigate 1.12 -> COD [default] -> by id [default] -> cod-id 1010369 -> default output -> quit"
rm -f fetched.fdf
printf '1.12\n\n\n1010369\n\n\n0\n' | stb-suite > log_menu.txt 2>&1
check_contains "Formula:.*Fe3O4" log_menu.txt
check_success fetched.fdf

echo "Testing: navigate 1.12 -> OPTIMADE -> twodmatpedia [default] -> by id -> 2dm-2127 -> default output -> quit"
rm -f fetched.fdf
printf '1.12\n3\n1\n1\n2dm-2127\n\n\n0\n' | stb-suite > log_menu_optimade.txt 2>&1
check_contains "Fetched:.*optimade 2dm-2127" log_menu_optimade.txt
check_contains "Formula:.*MoS" log_menu_optimade.txt
check_success fetched.fdf


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
