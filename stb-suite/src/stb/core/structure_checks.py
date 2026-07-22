"""Shared structure-malformation checklist, used by stb-inputfile, stb-fetch,
and stb-mlrelax -- the same 3 checks (atom proximity, lattice handedness,
atomic density) were independently duplicated in all three, each with a
comment acknowledging the thresholds are meant to stay identical. Extracted
here once mlrelax became the third near-identical copy, same "extract on
second/third use" policy as the rest of core/.
"""

import numpy as np

from stb.core import structure_io
from stb.core.cli import print_table

MIN_ATOM_DISTANCE_ANG = 0.5
MIN_DENSITY_ATOMS_PER_ANG3 = 0.01
MAX_DENSITY_ATOMS_PER_ANG3 = 0.15


def run_malformation_checks(pmg_structure, vacuum_axes, f_out=None,
                             min_atom_distance=MIN_ATOM_DISTANCE_ANG,
                             min_density=MIN_DENSITY_ATOMS_PER_ANG3,
                             max_density=MAX_DENSITY_ATOMS_PER_ANG3,
                             extra_checks=None):
    """Prints an explicit per-check table (Check | Result | Status), each row
    marked [OK]/[WARNING]/[SKIPPED] -- atom proximity (minimum-image pairwise
    distance), lattice handedness (cell determinant sign), and atomic density
    (bulk/non-vacuum-padded structures only; skipped otherwise, since a
    vacuum-padded cell's volume is dominated by the artificial vacuum by
    design).

    `extra_checks`, if given, is a list of (label, status, detail) tuples
    -- status one of "OK"/"WARNING"/"SKIPPED" -- appended as additional rows,
    e.g. stb-inputfile's own declared-atom-count header check, so every
    caller's checks render in one consistent table instead of mixing styles.
    """
    rows = []

    min_dist = structure_io.min_pairwise_distance(pmg_structure)
    if min_dist is None:
        rows.append((["Atom proximity", "only 1 atom -- not applicable", "[SKIPPED]"], None))
    elif min_dist < min_atom_distance:
        rows.append((["Atom proximity",
                       f"atoms unusually close: min. distance = {min_dist:.3f} Ang "
                       f"(< {min_atom_distance} Ang threshold)", "[WARNING]"], 'yellow'))
    else:
        rows.append((["Atom proximity",
                       f"min. distance = {min_dist:.3f} Ang (>= {min_atom_distance} Ang threshold)",
                       "[OK]"], None))

    det = np.linalg.det(pmg_structure.lattice.matrix)
    if det < 0:
        rows.append((["Lattice handedness",
                       "left-handed lattice (negative determinant) -- check vector order/signs",
                       "[WARNING]"], 'yellow'))
    else:
        rows.append((["Lattice handedness", "right-handed (positive determinant)", "[OK]"], None))

    if sum(vacuum_axes) > 0:
        rows.append((["Atomic density",
                       "skipped -- vacuum-padded structure (density metric not meaningful)",
                       "[SKIPPED]"], None))
    elif abs(det) > 0:
        density = len(pmg_structure) / abs(det)
        if min_density <= density <= max_density:
            rows.append((["Atomic density",
                           f"{density:.4f} atoms/Ang^3 (within [{min_density}, {max_density}])",
                           "[OK]"], None))
        else:
            rows.append((["Atomic density",
                           f"Unusual atomic density: {density:.4f} atoms/Ang^3 "
                           f"(outside [{min_density}, {max_density}])", "[WARNING]"], 'yellow'))

    for label, status, detail in (extra_checks or []):
        rows.append(([label, detail, f"[{status}]"], 'yellow' if status == "WARNING" else None))

    print_table(["Check", "Result", "Status"], rows, f_out)
