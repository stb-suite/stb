"""Optional third-party dependency guards shared across CLI tools."""

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
