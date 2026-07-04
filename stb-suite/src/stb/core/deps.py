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
