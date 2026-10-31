#!/bin/bash
# Guided example: 3.5 -- AIMD Trajectory Converter (stb-ani2traj; code 3.5 in the stb-suite menu).
#
# Not an automated test (see test/3-utils/5-ani2traj/test.sh for that) -- a commented walk-through: it runs real
# commands, one case at a time, on COPIES of siesta_runs/ inside output/<case>/, and shows the piece of the report
# that proves what just happened. siesta_runs/ itself is never touched. Pauses between cases (`yes '' | bash
# example_3.5.sh` runs it straight through). Safe to re-run -- it always starts by wiping its own output/.
#
# siesta_runs/ holds three real, finished SIESTA 5.4 AIMD runs: 8 water molecules (Nose NVT, 350 K), 8-atom Si
# (Verlet NVE) and 8-atom Si (Nose-Parrinello-Rahman NPT, a variable cell). No SIESTA is needed. README section 6
# explains every number.

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

# a2t <case> <run> <args...>: a fresh copy of siesta_runs/<run> in output/<case>/, stb-ani2traj run there
a2t() {
    local case=$1 run=$2; shift 2
    rm -rf "$OUT/$case"; cp -a "$DIR/siesta_runs/$run" "$OUT/$case"
    echo "  \$ stb-ani2traj $*"
    (cd "$OUT/$case" && stb-ani2traj --no-intro "$@" < /dev/null 2>&1 | strip_colors > log.txt) || true
}
log() { echo "$OUT/$1/log.txt"; }
warnings() { section "$(log "$1")" "WARNINGS" | grep -E "^  [0-9]+\." | fold -s -w 116 | sed 's/^/  /'; }


banner "Welcome: a SIESTA AIMD trajectory, with its cell, time and thermodynamics"
cat <<'EOF'
SIESTA writes the frames of a molecular-dynamics run to <label>.ANI -- a plain multi-frame XYZ with no cell, no time
and no energies. stb-ani2traj writes the same frames as a trajectory that carries the cell (per frame for a
variable-cell run), and, in extended XYZ, each frame's time, temperature and energies; its report says what the run
was, which MD step each frame is, how the cell, temperature and energy behaved and how the atoms moved.

The tour:
   1  water, NVT: the whole report           5  equilibration and stride (pdb)
   2  which MD step is each frame?            6  what the report catches
   3  Si, NPT: a cell that changes            7  the interactive menu
   4  Si, NVE: energy conservation; --unwrap
EOF
pause


banner "Case 1 -- 8 water molecules, Nose NVT: the whole report"
a2t case1 water_nvt -l water
section "$(log case1)" "INPUTS" | sed 's/^/  /' | cut -c1-140
section "$(log case1)" "SYSTEM" | sed 's/^/  /'
section "$(log case1)" "THERMODYNAMICS" | sed 's/^/  /' | cut -c1-140
section "$(log case1)" "MOTION" | sed 's/^/  /'
section "$(log case1)" "OUTPUT" | sed 's/^/  /' | cut -c1-140
cat <<'EOF'
  -> The .ANI has no cell: it comes from each step's 'outcell' block in the log. The temperature averages the
     350 K target, but in 144 fs it still swings more than a 24-atom canonical ensemble would (the Nose thermostat
     has not settled). The xsf opens in OVITO, VMD or XCrySDen with its periodic box.
EOF
pause


banner "Case 2 -- Which MD step is each frame?"
section "$(log case1)" "FRAMES AND MD STEPS" | fold -s -w 116 | sed 's/^/  /'
cat <<'EOF'
  -> These runs have WriteCoorStep true, so the log prints the geometry of every step and the tool can match each
     frame to it: SIESTA writes a frame AFTER it moves the atoms, so frame k is the geometry of step k+1 -- and the
     last frame is the .XV. The energy and temperature written next to frame k are those of step k+1, not step k.
EOF
pause


banner "Case 3 -- Si, Nose-Parrinello-Rahman NPT: a cell that changes"
a2t case3 si_npt -l si -of xyz
section "$(log case3)" "LATTICE" | sed 's/^/  /'
section "$(log case3)" "THERMODYNAMICS" | grep -E "^ *(\||T |V |P |E_)|unit" | sed 's/^/  /' | cut -c1-140
echo "  the header line of the first frame:"; sed -n 2p "$OUT/case3/si_traj.xyz" | fold -w 112 | sed 's/^/    /'
cat <<'EOF'
  -> Every frame carries its own cell (Lattice="..."), so OVITO shows the box growing. The extended XYZ also
     carries each frame's MD step, time, T, E_KS, E_tot, volume and pressure (figures/thermo.png reads them back).
EOF
pause


banner "Case 4 -- Si, NVE: energy conservation; positions and --unwrap"
a2t case4 si_nve -l si -of xyz
section "$(log case4)" "THERMODYNAMICS" | grep -E "E_tot|E_KS|T  |drift" | sed 's/^/  /' | cut -c1-140
section "$(log case4)" "MOTION" | grep -E "Boundary|Outside|unwrap|Largest" | sed 's/^/  /'
cat <<'EOF'
  -> In NVE the kinetic and potential energies trade places while E_tot stays flat: its drift per atom and ps is
     the check of the time step and the SCF tolerance. The Si atoms that start on a cell face vibrate to negative
     fractional coordinates and stay there: SIESTA does not wrap .ANI positions, so there are no jumps to remove and
     --unwrap changes nothing; it matters for trajectories that were wrapped by another program.
EOF
pause


banner "Case 5 -- Equilibration and stride, pdb"
a2t case5 water_nvt -l water -of pdb --skip 100 --stride 5
section "$(log case5)" "FRAMES AND MD STEPS" | grep "Kept" | sed 's/^/  /'
section "$(log case5)" "OUTPUT" | sed 's/^/  /' | cut -c1-140
cat <<'EOF'
  -> --skip drops the first frames (equilibration), --stride keeps every Nth. PDB stores coordinates with 3
     decimals, which the read-back check states; it is VMD's own format.
EOF
pause


banner "Case 6 -- What the report catches"
rm -rf "$OUT/case6_nolog"; cp -a "$DIR/siesta_runs/si_npt" "$OUT/case6_nolog"; rm "$OUT/case6_nolog/calc.out"
echo "  \$ stb-ani2traj -l si --dry-run          (an NPT run whose log is missing)"
(cd "$OUT/case6_nolog" && stb-ani2traj --no-intro -l si --dry-run < /dev/null 2>&1 | strip_colors > log.txt) || true
warnings case6_nolog
rm -rf "$OUT/case6_app"; cp -a "$DIR/siesta_runs/si_nve" "$OUT/case6_app"
cat "$OUT/case6_app/si.ANI" "$OUT/case6_app/si.ANI" > "$OUT/case6_app/x" && mv "$OUT/case6_app/x" "$OUT/case6_app/si.ANI"
echo "  \$ stb-ani2traj -l si --dry-run          (a .ANI with a second run appended)"
(cd "$OUT/case6_app" && stb-ani2traj --no-intro -l si --dry-run < /dev/null 2>&1 | strip_colors > log.txt) || true
warnings case6_app
cat <<'EOF'
  -> Without the log a variable-cell run gets one static cell -- wrong, and said so. SIESTA appends to an existing
     .ANI when a run is restarted in the same folder: more frames than MD steps is flagged (and a last frame that is
     not the .XV would be too -- here the appended copy ends on the same geometry).
EOF
pause


banner "Case 7 -- The interactive menu (stb-suite, code 3.5)"
rm -rf "$OUT/case7"; cp -a "$DIR/siesta_runs/si_npt" "$OUT/case7"
# label (Enter: the only *.ANI), format 3 = xyz, skip 0, stride 1, no unwrap, default output, no dry run, no report
(cd "$OUT/case7" && printf '3.5\n\n3\n0\n1\nn\n\nn\nn\n\n0\n' | stb-suite 2>&1 | strip_colors > log.txt) || true
grep -E "^Input  |^Output  |^Read-back" "$OUT/case7/log.txt" | sed 's/^/  /'
cat <<'EOF'
  -> The menu asks for the label (Enter for the only *.ANI), the format, the frames to keep, --unwrap, the output
     name, dry run and report, and runs the same command.
EOF
pause


banner "Done"
cat <<EOF
Everything is in $OUT/. What to take from the tour:
  - .ANI has no cell; the converted trajectory has it per frame, and the extended XYZ also time and thermodynamics.
  - Frame k is the geometry AFTER MD step k; each frame gets the values of its own step.
  - Read THERMODYNAMICS (T vs target, E_tot drift) and MOTION (largest step, crossings) before using the run.
EOF
