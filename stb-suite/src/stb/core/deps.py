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


def _pkg_resources_stand_in():
    """The two names pybader's __init__ takes from `pkg_resources`, built on the
    standard library. setuptools 81 dropped `pkg_resources`, and pybader 0.3.12
    still imports it, so on a fresh install `import pybader` fails although the
    package is there."""
    import types
    from importlib import metadata, resources

    def iter_entry_points(group, name=None):
        for entry_point in metadata.entry_points(group=group):
            if name is None or entry_point.name == name:
                yield entry_point

    def resource_listdir(package, resource_name):
        folder = resources.files(package).joinpath(resource_name)
        return [item.name for item in folder.iterdir()]

    module = types.ModuleType("pkg_resources")
    module.iter_entry_points = iter_entry_points
    module.resource_listdir = resource_listdir
    return module


def require_pybader():
    """Imports and returns pybader.interface, or exits with a clear error:
    'not installed' when the package is missing, the original exception when it
    is installed but fails to import."""
    stand_in = None
    try:
        import pkg_resources  # noqa: F401  (pybader needs it)
    except ImportError:
        # Only while pybader loads: another library must not find this
        # two-function module in place of the real one.
        stand_in = sys.modules["pkg_resources"] = _pkg_resources_stand_in()
    try:
        import pybader.interface
    except ImportError as error:
        if getattr(error, "name", None) == "pybader":
            print("\n\033[91m[CRITICAL ERROR] pybader library not found.\033[0m")
            print("Please install it using: pip install pybader")
        else:
            print("\n\033[91m[CRITICAL ERROR] pybader is installed but could not be imported:\033[0m")
            print(f"  {type(error).__name__}: {error}")
            print("Reinstalling pybader will not fix this; the missing module above is what it needs.")
        sys.exit(1)
    finally:
        if stand_in is not None and sys.modules.get("pkg_resources") is stand_in:
            del sys.modules["pkg_resources"]
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
