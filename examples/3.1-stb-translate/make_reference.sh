#!/bin/bash
# =====================================================================================================
#  make_reference.sh -- regenerates reference_results/ of the 3.1 example (stb-translate)
#
#  Every case runs in reference_results/<case>/: the screen report (colours stripped) goes to report.txt, the
#  converted files stay next to it. Seconds on any machine; no model, no SIESTA. The inputs are in this folder
#  (python3 make_inputs.py rewrites the hand-made ones).
#
#      source ~/.anaconda3/etc/profile.d/conda.sh && conda activate stb    # any environment with stb_suite
#      bash make_reference.sh
# =====================================================================================================
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
OUT="$HERE/reference_results"
rm -rf "$OUT"
mkdir -p "$OUT"
export MPLBACKEND=Agg
nocolor() { sed 's/\x1b\[[0-9;]*m//g'; }

run() {   # run <case> <stb-translate args...>  (paths relative to this folder; outputs land in the case folder)
    local case="$1"; shift
    mkdir -p "$OUT/$case"
    echo "  [run]  $case: stb-translate $*"
    (cd "$OUT/$case" && stb-translate --no-intro "$@" < /dev/null 2>&1 | nocolor > report.txt) || true
}

echo "== 1. A relaxed SIESTA structure to CIF (formats from the names)"
run struct_out_to_cif -i ../../si111_adatom_relaxed.STRUCT_OUT -o si111_adatom.cif --view-image si111_adatom.png

run struct_out_to_cif_p1 -i ../../si111_adatom_relaxed.STRUCT_OUT -o si111_adatom_p1.cif --cif-symprec 0
run cif_symmetrize -i ../../si_noisy_0.002.fdf -o si_noisy_symmetrized.cif

echo "== 2. One crystal, every output format"
for f in si.cif si.vasp si.xyz si.xsf geometry.in si.gen si.data si_copy.fdf; do
    run "every_format" -i ../../si_bulk.fdf -o "$f"
    mv "$OUT/every_format/report.txt" "$OUT/every_format/report_${f}.txt"
done

echo "== 3. Fractional and Cartesian coordinates (-cf) on an oblique cell"
run coords_direct -i ../../si111_2x2.fdf -o slab_direct.fdf -cf direct
run coords_cartesian -i ../../si111_2x2.fdf -o slab_cartesian.fdf -cf cartesian

echo "== 4. Units: a lattice in Bohr"
run bohr -i ../../si_bohr.fdf -o si_primitive.vasp

echo "== 5. Symmetry and tolerance: the same crystal with noise"
for a in 0.002 0.02 0.2; do run "noisy_$a" -i "../../si_noisy_$a.fdf" --dry-run; done

echo "== 6. A CIF with its symmetry operations (asymmetric unit) to fdf"
run cif_asymmetric -i ../../nacl_asym.cif -o nacl.fdf

echo "== 7. A molecule: XYZ has no cell"
run molecule_xyz -i ../../water.xyz -o water.fdf --lattice-vectors 15 0 0 0 15 0 0 0 15

echo "== 8. What the report catches"
run mixup -i ../../si_mixup.vasp --dry-run
run wrongscale -i ../../si_wrongscale.vasp --dry-run

echo "== 9. Batch: every .fdf of this folder to CIF"
run batch -i '../../si*.fdf' -of cif --out-dir cif

echo "== 10. A slab to LAMMPS (triclinic cell) and back through the read-back check"
run lammps -i ../../si111_2x2.fdf -o slab.data

echo "done -- reference_results/<case>/report.txt"
