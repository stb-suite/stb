#!/bin/bash

# --- Setup ---
# Checks of the bundled PseudoDojo bank that need no SIESTA (check_bank.py): every folder against its
# manifest.json (files, md5, no ekb = 0 projectors), the projector sets of each table, the sr3plus
# valence, and the guards of core/pseudopotentials on small fixtures (ekb = 0 refused, fr files with no
# scalar part without spin-orbit, synthetic-atoms per folder index, Li/Be non-perturbative).
# The SIESTA regression of the bank (every element of every table, 1 and 8 MPI ranks) is
# xc_regression/; it needs SIESTA and MPI and is run by hand on a cluster, see its README.md.
FIXTURE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${PYTHON:-python}"

GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

if ! "$PYTHON" -c "import stb.core.pseudopotentials" 2>/dev/null; then
    echo -e "${RED}[SKIP]${NC} the stb package is not importable with '$PYTHON'"
    exit 0
fi

echo "--- Bundled PseudoDojo bank (no SIESTA) ---"
"$PYTHON" "$FIXTURE_DIR/check_bank.py"
status=$?
if [ $status -eq 0 ]; then
    echo -e "${GREEN}Bank checks passed${NC}"
else
    echo -e "${RED}Bank checks failed${NC}"
fi
exit $status
