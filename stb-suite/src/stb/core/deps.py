"""Optional third-party dependency guards shared across CLI tools, plus a couple of
tiny sisl-based helpers shared by the tools that need sisl's own Geometry object
directly (not this suite's own FdfStructure -- see core/structure_io.py for that)."""

import os
import sys


def require_sisl():
    """Imports and returns sisl, or exits with a clear error if it's missing."""
    try:
        import sisl
    except ImportError:
        print("\n\033[91m[CRITICAL ERROR] sisl library not found.\033[0m")
        print("Please install it using: pip install sisl")
        sys.exit(1)
    return sisl


def require_icet():
    """Imports and returns icet, or exits with a clear error if it's missing."""
    try:
        import icet
    except ImportError:
        print("\n\033[91m[CRITICAL ERROR] icet library not found.\033[0m")
        print("Please install it using: pip install icet")
        sys.exit(1)
    return icet


def require_pyxtal():
    """Imports and returns pyxtal, or exits with a clear error if it's missing."""
    try:
        import pyxtal
    except ImportError:
        print("\n\033[91m[CRITICAL ERROR] pyxtal library not found.\033[0m")
        print("Please install it using: pip install pyxtal")
        sys.exit(1)
    return pyxtal


def require_pybader():
    """Imports and returns pybader.interface, or exits with a clear error if
    it's missing."""
    try:
        import pybader.interface
    except ImportError:
        print("\n\033[91m[CRITICAL ERROR] pybader library not found.\033[0m")
        print("Please install it using: pip install pybader")
        sys.exit(1)
    return pybader.interface


def read_sisl_geometry_xv_or_fdf(label):
    """Reads <label>.XV via sisl if present, else <label>.fdf, else returns
    (None, None). Shared by stb-bader and stb-workfunction, both of which
    need sisl's own Geometry object (for cube-writing / vacuum-axis
    detection) rather than this suite's lighter FdfStructure. Returns
    (geometry, source_path) so callers can name the file in messages;
    raises whatever sisl itself raises on a present-but-corrupt file --
    callers are expected to catch that themselves, since what to do about
    it (abort vs. degrade gracefully) differs by tool.
    """
    import sisl
    file_xv = f"{label}.XV"
    file_fdf = f"{label}.fdf"
    if os.path.exists(file_xv):
        return sisl.get_sile(file_xv).read_geometry(), file_xv
    if os.path.exists(file_fdf):
        return sisl.get_sile(file_fdf).read_geometry(), file_fdf
    return None, None


def require_mace():
    """Imports and returns the mace.calculators module, or exits with a clear
    error if it (or its PyTorch dependency) is missing. This is the suite's
    only heavy/optional dependency -- not part of the core install.
    """
    try:
        import mace.calculators
    except ImportError:
        print("\n\033[91m[CRITICAL ERROR] mace-torch (and/or PyTorch) not found.\033[0m")
        print("This tool needs the optional 'ml' extra. Install it with:")
        print("  pip install stb_suite[ml]")
        print("(or directly: pip install torch mace-torch)")
        sys.exit(1)
    return mace.calculators
