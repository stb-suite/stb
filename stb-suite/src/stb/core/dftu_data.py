"""Single source of truth for DFT+U/Hubbard-U reference data.

Consolidates SHELL_NAMES/DEFAULT_SHELL and REFERENCE_U, which used to be
duplicated identically in hubbardu.py and hubbardu_analysis.py, now also
used by dftu.py -- third consumer, past the "extract on second use" point
stated in CLAUDE.md.
"""

from __future__ import annotations

# Standard correlated shell per element (periodic-table block classification --
# uncontroversial chemistry, unlike the Hubbard U value itself).
SHELL_NAMES = {"3d": (3, 2), "4d": (4, 2), "5d": (5, 2), "4f": (4, 3), "5f": (5, 3)}

_3D = ["Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn"]
_4D = ["Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd"]
_5D = ["La", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg"]
_4F = ["Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]
_5F = ["Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"]

DEFAULT_SHELL = {
    **{el: "3d" for el in _3D}, **{el: "4d" for el in _4D}, **{el: "5d" for el in _5D},
    **{el: "4f" for el in _4F}, **{el: "5f" for el in _5F},
}

# Reference GGA+U values (eV), for comparison ONLY -- never substituted for a
# user-supplied or linear-response-computed U. Materials Project's oxide
# calibration (docs.materialsproject.org/methodology/materials-methodology/
# calculation-details/gga+u-calculations/hubbard-u-values), itself derived
# from Wang, Maxisch & Ceder, Phys. Rev. B 73, 195107 (2006). Calibrated for
# VASP GGA+U on oxides -- a starting-point sanity check, not a validation:
# the actual U depends on functional, pseudopotential, and basis set.
REFERENCE_U = {
    "V": 3.25, "Cr": 3.7, "Mn": 3.9, "Fe": 5.3,
    "Co": 3.32, "Ni": 6.2, "Mo": 4.38, "W": 6.2,
}


def ldau_proj_block(entries):
    """Formats a %block LDAU.proj snippet for one or more species.

    `entries` is a list of dicts, each with keys 'species', 'n', 'l', 'u',
    'j' (eV), and optionally 'rc'/'omega' (default 0.0/0.0 -- SIESTA computes
    its own defaults then). One 4-line stanza per entry, verified against
    SIESTA's own tutorial (siesta-project/tutorials, DFT+U/MnO/MnO.fdf) and
    documentation.
    """
    lines = ["%block LDAU.proj\n"]
    for e in entries:
        rc = e.get('rc', 0.0)
        omega = e.get('omega', 0.0)
        lines.append(f"{e['species']}   1\n")
        lines.append(f"n={e['n']}    {e['l']}\n")
        lines.append(f"{e['u']:.3f}    {e['j']:.3f}\n")
        lines.append(f"{rc:.3f}    {omega:.3f}\n")
    lines.append("%endblock LDAU.proj\n")
    return "".join(lines)
