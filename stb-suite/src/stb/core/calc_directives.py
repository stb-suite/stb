"""Regex substitute-or-append helpers for forcing a specific directive onto
a calc.fdf template's text (kept as a plain string throughout the suite,
since it's user-supplied and its exact formatting/comments are otherwise
preserved verbatim) -- shared by every WORKFLOW_TOOLS prep tool that needs
to write an INDEPENDENT (non-coupled) folder per candidate/image, where
letting SIESTA relax on its own would silently move the geometry away from
the specific point being sampled.
"""

import re

_MD_RUNTYPE_RE = re.compile(r'MD\.TypeOfRun\s+\S+', re.IGNORECASE)
_MD_NUMCGSTEPS_RE = re.compile(r'MD\.NumCGsteps\s+\d+', re.IGNORECASE)


def force_single_point(calc_text):
    """Substitutes/appends 'MD.TypeOfRun CG' and 'MD.NumCGsteps 0' onto
    calc_text (regex .subn, append if the tag isn't present at all rather
    than erroring, since an absent tag is a normal template, not a
    malformed one). Moved here from neb.py once stackingfault.py became a
    second consumer -- same reasoning in both callers: each folder (a NEB
    image, a stacking-fault grid point) is an independent, single-point
    evaluation with no inter-folder coupling once written to disk, so
    letting SIESTA relax one on its own would silently move it away from
    the specific point/geometry being sampled. Forced unconditionally, no
    opt-out flag, in both callers.
    """
    new_text, count = _MD_RUNTYPE_RE.subn('MD.TypeOfRun          CG', calc_text)
    if count == 0:
        new_text += "\nMD.TypeOfRun          CG\n"
    new_text, count = _MD_NUMCGSTEPS_RE.subn('MD.NumCGsteps         0', new_text)
    if count == 0:
        new_text += "MD.NumCGsteps         0\n"
    return new_text
