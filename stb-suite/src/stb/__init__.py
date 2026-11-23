"""STB-SUITE, the Siesta Toolbox Suite."""

from stb._input_header import install as _install_input_header
from stb._version import __version__

_install_input_header()

__all__ = ["__version__"]
