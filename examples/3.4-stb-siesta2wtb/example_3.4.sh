#!/bin/bash
# Guided example: 3.4 -- Wantibexos Interface (stb-siesta2wtb; code 3.4 in the stb-suite menu).
#
# Not an automated test (see test/3-utils/4-wantibexos/ for the fixture) -- a commented walk-through: it runs real
# commands, one case at a time, on COPIES of siesta_runs/ inside output/<case>/, and shows the piece of the report
# that proves what just happened. siesta_runs/ itself is never touched. Pauses between cases (`yes '' | bash
# example_3.4.sh` runs it straight through). Safe to re-run -- it always starts by wiping its own output/.
#
# siesta_runs/ holds two real, finished SIESTA 5.4 runs with SaveHS true: an hBN monolayer (no spin) and
# ferromagnetic bcc Fe (collinear spin). No SIESTA and no WanTiBEXOS are needed. README section 6 explains every number.

set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"
OUT="$DIR/output"
rm -rf "$OUT"
mkdir -p "$OUT"

pause() {
    echo
    read -p "  [Press Enter to continue] " -r || true
    echo
}

banner() {
    echo "=================================================================="
    echo " $1"
    echo "=================================================================="
}

strip_colors() { sed 's/\x1b\[[0-9;]*m//g'; }

section() {   # $1 = file, $2 = start of a section title: from "[N] TITLE" to the next "[M] "
    strip_colors < "$1" | awk -v t="$2" '/^\[[0-9]+\] [A-Z]/ {if (on) exit; if (index($0, "] " t) > 0) on = 1} on {print}' \
        | grep -v "^-----"
}

# wtb <case> <run> <args...>: a fresh copy of siesta_runs/<run> in output/<case>/, stb-siesta2wtb run there
wtb() {
    local case=$1 run=$2; shift 2
    rm -rf "$OUT/$case"; cp -a "$DIR/siesta_runs/$run" "$OUT/$case"
    echo "  \$ stb-siesta2wtb $*"
    (cd "$OUT/$case" && stb-siesta2wtb --no-intro "$@" < /dev/null 2>&1 | strip_colors > log.txt) || true
}
log() { echo "$OUT/$1/log.txt"; }
warnings() { section "$(log "$1")" "WARNINGS" | grep -E "^  [0-9]+\." | fold -s -w 116 | sed 's/^/  /'; }


banner "Welcome: a SIESTA Hamiltonian as a WanTiBEXOS tight-binding model"
cat <<'EOF'
SIESTA works in a basis of localized atomic orbitals, so its Kohn-Sham Hamiltonian is a tight-binding model as it
stands: H(R) and S(R) between the orbitals of the cell and of a few neighbouring cells. With SaveHS true SIESTA writes
them to <label>.HSX; stb-siesta2wtb writes them as tb-<NP|SP|SOC>.ham + basis_set-<..>.basis for WanTiBEXOS, and its
report describes the model and checks the written file by rebuilding the bands from it.

The tour:
   1  hBN (no spin): the whole report        4  what 6 decimals do to the bands
   2  bcc Fe (spin): two blocks, 729 cells   5  what the report catches
   3  the Fermi level and the zero of H      6  the interactive menu
EOF
pause


banner "Case 1 -- hBN monolayer: the whole report"
wtb case1 hbn -i calc.fdf -o calc.out
section "$(log case1)" "SIESTA RUN" | grep -E "Hamiltonian|Fermi level|Energy zero" | fold -s -w 116 | sed 's/^/  /'
section "$(log case1)" "SYSTEM" | grep -E "Dimensionality|species|^B |^N |Orbitals" | sed 's/^/  /'
section "$(log case1)" "HAMILTONIAN" | sed 's/^/  /'
section "$(log case1)" "FERMI LEVEL" | grep -E "EIG|Occupied" | sed 's/^/  /'
section "$(log case1)" "READ-BACK" | sed 's/^/  /' | cut -c1-140
section "$(log case1)" "FOR WANTIBEXOS" | sed 's/^/  /'
echo "  the first lines of tb-NP.ham:"; head -n 11 "$OUT/case1/tb-NP.ham" | sed 's/^/    /'
cat <<'EOF'
  -> 13 DZP orbitals per atom, reaching 5.2 Ang: H couples the cell to its 24 nearest cells in the plane (nsc 5x5x1)
     and none across the vacuum. The bands rebuilt from the written file reproduce sisl's within 5 eV of E_F to
     ~0.01 meV at every k checked, with the 4 occupied sigma/pi bands below E_F (figures/hbn_bands.png: the whole
     Gamma-M-K-Gamma path against SIESTA's own hbn.bands).
EOF
pause


banner "Case 2 -- bcc Fe: two spin blocks"
wtb case2 fe -i calc.fdf -o calc.out
section "$(log case2)" "SYSTEM" | grep -E "^Fe |species" | sed 's/^/  /'
section "$(log case2)" "HAMILTONIAN" | grep -E "Spin|Neighbour|Stored|Longest|Layout" | sed 's/^/  /'
section "$(log case2)" "OUTPUT FILES" | grep -E "Hamiltonian|Body" | sed 's/^/  /'
section "$(log case2)" "READ-BACK" | grep -E "Non-zero|Verdict" | sed 's/^/  /'
echo "  the on-site 3s element of each spin block:"
awk 'NR > 9 && $1 == 0 && $2 == 0 && $3 == 0 && (($4 == 1 && $5 == 1) || ($4 == 20 && $5 == 20))' "$OUT/case2/tb-SP.ham" \
    | sed 's/^/    /'
cat <<'EOF'
  -> 19 orbitals per spin (with the 3s/3p semicore), diffuse 4s/4p reaching 8 Ang: 9x9x9 = 729 cells, one million
     lines (84 MB), 1.4 % of them non-zero. Up is orbitals 1-19, down 20-38. Even the deep 3s level is 2.5 eV lower
     for spin up -- the exchange field of the 3d moment (2.22 mu_B). figures/fe_bands.png checks both blocks.
EOF
pause


banner "Case 3 -- The Fermi level and the zero of H"
section "$(log case1)" "FERMI LEVEL" | grep -E "Fermi level|Reference" | fold -s -w 116 | sed 's/^/  /'
warnings case1 | grep -A3 "already referred" || true
wtb case3_f0 hbn -i calc.fdf -f 0 --no-check
echo "  with -f 0:"; head -n 3 "$OUT/case3_f0/tb-NP.ham" | sed 's/^/    /'
section "$(log case3_f0)" "FERMI LEVEL" | grep -E "Reference" | sed 's/^/  /'
cat <<'EOF'
  -> A version-2 .HSX stores E_F, and sisl shifts H by -E_F S on reading: the matrix elements in the file are already
     referred to E_F. With -o the header also carries E_F; a reader that subtracts it from H moves every band by
     -E_F (+4.02 eV here). With -f 0 the header matches H. Which one is right depends on how the reading program uses
     the header -- the report states both facts.
EOF
pause


banner "Case 4 -- What 6 decimals do to the bands"
section "$(log case1)" "HAMILTONIAN" | grep "S(Gamma)" | sed 's/^/  /'
section "$(log case1)" "READ-BACK" | grep -E "^k |^\(" | sed 's/^/  /' | cut -c1-140
cat <<'EOF'
  -> The file has 6 decimals. Where S(k) is nearly singular (smallest eigenvalue 6.6e-4 at Gamma) the rounding moves
     the highest, unphysical bands by tens of meV -- 70 to 170 eV above E_F. Within 5 eV of E_F the change is below
     0.01 meV (figures/precision.png: 12 random k-points, hBN and Fe).
EOF
pause


banner "Case 5 -- What the report catches"
rm -rf "$OUT/case5_nohsx"; cp -a "$DIR/siesta_runs/hbn" "$OUT/case5_nohsx"; rm -f "$OUT/case5_nohsx/hbn.HSX"
echo "  \$ stb-siesta2wtb -i calc.fdf -o calc.out          (hbn.HSX removed)"
(cd "$OUT/case5_nohsx" && stb-siesta2wtb --no-intro -i calc.fdf -o calc.out < /dev/null 2>&1 | strip_colors > log.txt) || true
echo "  a run without SaveHS (no hbn.HSX):"; grep FAIL "$(log case5_nohsx)" | sed 's/^/    /'
wtb case5_badlog hbn -i calc.fdf -o siesta.out
echo "  a log that is not there:"; grep FAIL "$(log case5_badlog)" | sed 's/^/    /'
wtb case5_dry fe -i calc.fdf -o calc.out --dry-run
echo "  --dry-run on Fe (nothing written):"; section "$(log case5_dry)" "OUTPUT FILES" | sed 's/^/    /'
ls "$OUT/case5_dry" | grep -c "tb-" | sed 's/^/    tb files in the folder: /' || true
cat <<'EOF'
  -> No .HSX means the run did not save H: rerun with SaveHS true. --dry-run says how many lines and roughly how many
     MB would be written before writing them.
EOF
pause


banner "Case 6 -- The interactive menu (stb-suite, code 3.4)"
rm -rf "$OUT/case6"; cp -a "$DIR/siesta_runs/hbn" "$OUT/case6"
# fdf, log, no manual E_F, no dry run, check yes (Enter), no report, Enter, 0 = exit
(cd "$OUT/case6" && printf '3.4\ncalc.fdf\ncalc.out\n\nn\n\nn\n\n0\n' | stb-suite 2>&1 | strip_colors > log.txt) || true
grep -E "^Input  |^Fermi level  |^Mode|^Verdict" "$OUT/case6/log.txt" | sed 's/^/  /'
cat <<'EOF'
  -> The menu asks for the .fdf, the log, a Fermi level by hand, dry run, the read-back check and the report, and
     runs the same command.
EOF
pause


banner "Done"
cat <<EOF
Everything is in $OUT/. What to take from the tour:
  - A SIESTA Hamiltonian is a tight-binding model as it stands; the files carry H and S for every neighbouring cell.
  - The read-back check rebuilds the bands from the written file: ~0.1 meV from SIESTA's within 5 eV of E_F.
  - With a modern .HSX the matrix elements are already referred to E_F: choose the header's Fermi level (-o or -f 0)
    for the program that reads it.
EOF
