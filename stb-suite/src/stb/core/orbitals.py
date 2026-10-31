"""Shared SIESTA orbital-naming vocabulary (angular momentum l / magnetic
quantum number m -> orbital name, e.g. 'p', 'px', 'dxy').

Extracted from dos.py once stb-fatbands became a second consumer -- both
tools categorize per-orbital data (DOS weight in dos.py, wavefunction
projection weight in fatbands.py) into the same s/p/d/f (or detailed
px/py/pz/dxy/...) buckets, so there is exactly one naming scheme to keep
consistent rather than two independently maintained copies.
"""

# This defines the standard SIESTA order for real spherical harmonics
ORBITAL_MAP = {
    0: {0: 's'},
    1: {-1: 'py', 0: 'pz', 1: 'px'},
    2: {-2: 'dxy', -1: 'dyz', 0: 'dz2', 1: 'dxz', 2: 'dx2-y2'},
    3: {-3: 'f-3', -2: 'f-2', -1: 'f-1', 0: 'f0', 1: 'f1', 2: 'f2', 3: 'f3'} # Using simple f names
}

ORBITAL_ORDER = [
    's',
    'py', 'pz', 'px',
    'p', # For 'l' mode
    'dxy', 'dyz', 'dz2', 'dxz', 'dx2-y2',
    'd', # For 'l' mode
    'f-3', 'f-2', 'f-1', 'f0', 'f1', 'f2', 'f3',
    'f' # For 'l' mode
]


def get_orbital_name(l_val):
    """Maps angular momentum number 'l' to its name (s, p, d, f)."""
    l_map = {0: 's', 1: 'p', 2: 'd', 3: 'f'}
    return l_map.get(l_val, None) # Return None if not s,p,d,f


def get_detailed_orbital_name(l_val, m_val):
    """Maps (l, m) to orbital name (s, px, py, pz, dxy, ...)."""
    if l_val in ORBITAL_MAP:
        return ORBITAL_MAP[l_val].get(m_val, None) # Return None if m is invalid
    return None # Return None if l is invalid
