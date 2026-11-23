"""First line of every SIESTA input file the suite writes.

Tools write .fdf files in many ways: their own open() calls, shutil.copy of a template, phonopy's
SIESTA writer. An audit hook records every .fdf opened for writing by this process, and at exit each
one gets INPUT_HEADER as its first line. A comment changes nothing for SIESTA, and the header is never
repeated when a generated file is used as the template of a later stage.
"""

import atexit
import os
import sys
import tempfile

INPUT_HEADER = "# Created/modified with STB-SUITE (https://stb-suite.github.io/stb/)"
_OLD_HEADER = b"# automatic create using stb-translate"

_written = set()
_stamping = False
_installed = False


def mark_input(path):
    """Registers a file the suite wrote, for an input whose name does not end in .fdf."""
    _written.add(os.path.abspath(os.fsdecode(path)))


def stamped(data):
    """Returns the bytes of an input file with the header as its only copy, on the first line."""
    newline = b"\r\n" if b"\r\n" in data else b"\n"
    header = INPUT_HEADER.encode()
    lines = [ln for ln in data.split(newline) if ln.rstrip() != header and not ln.startswith(_OLD_HEADER)]
    return header + newline + newline.join(lines)


def _stamp_file(path):
    target = os.path.realpath(path)
    with open(target, "rb") as f:
        data = f.read()
    new = stamped(data)
    if new == data:
        return
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(target), prefix=".stb_header_")
    try:
        with os.fdopen(fd, "wb") as f:
            f.write(new)
        os.chmod(tmp, os.stat(target).st_mode & 0o7777)
        os.replace(tmp, target)
    except BaseException:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise


def stamp_written():
    """Stamps every recorded file; a file that is gone or cannot be rewritten is left alone."""
    global _stamping
    _stamping = True
    try:
        for path in sorted(_written):
            try:
                _stamp_file(path)
            except OSError:
                pass
        _written.clear()
    finally:
        _stamping = False


def _is_fdf(path):
    if isinstance(path, int) or path is None:
        return None
    name = os.fsdecode(path)
    return name if name.lower().endswith(".fdf") else None


def _audit(event, args):
    if _stamping or event not in ("open", "os.rename"):
        return
    try:
        if event == "open":
            path, mode, flags = args
            writes = (any(c in mode for c in "wax+") if isinstance(mode, str)
                      else bool(flags & (os.O_WRONLY | os.O_RDWR)))
            name = _is_fdf(path) if writes else None
        else:
            name = _is_fdf(args[1])
        if name:
            _written.add(os.path.abspath(name))
    except Exception:           # an audit hook that raises would make the open() itself fail
        pass


def install():
    global _installed
    if not _installed:
        _installed = True
        sys.addaudithook(_audit)
        atexit.register(stamp_written)
