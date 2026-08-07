"""Shared terminal UI helpers for the STB-SUITE command-line tools.

Consolidates COLORS/color_text()/show_intro(), which used to be duplicated
(with minor, purposeless drift) in 26 of the 27 modules under stb/, and
get_input()/get_float_input()/get_int_input(), which used to live only in
stb_suite.py.

Banner content (title/subtitle/footer lines) is NOT hardcoded here -- each
tool keeps its own text and passes it to show_intro(); only the printing
mechanics (ASCII logo, colors, borders, line-by-line reveal) are shared.
"""

from __future__ import annotations

import os
import io
import re
import sys
import time
import warnings
import threading
import contextlib
from time import sleep

COLORS = {
    'reset': '\033[0m',
    'cyan': '\033[96m',
    'blue': '\033[94m',
    'green': '\033[92m',
    'yellow': '\033[93m',
    'red': '\033[91m',
    'bold': '\033[1m',
    'underline': '\033[4m',
    'magenta': '\033[95m',
    'bg_red': '\033[41m',
    'white': '\033[97m',
}

_LOGO = r"""
.----------------.  .----------------.  .----------------.
| .--------------. || .--------------. || .--------------. |
| |    _______   | || |  _________   | || |   ______     | |
| |   /  ___  |  | || | |  _   _  |  | || |  |_   _ \    | |
| |  |  (__ \_|  | || | |_/ | | \_|  | || |    | |_) |   | |
| |   '.___`-.   | || |     | |      | || |    |  __'.   | |
| |  |`\____) |  | || |    _| |_     | || |   _| |__) |  | |
| |  |_______.'  | || |   |_____|    | || |  |_______/   | |
| |              | || |              | || |              | |
| '--------------' || '--------------' || '--------------' |
 '----------------'  '----------------'  '----------------'
"""


def color_text(text: str, color: str) -> str:
    """Returns text formatted with an ANSI color code, falling back to no color."""
    return f"{COLORS.get(color, COLORS['reset'])}{text}{COLORS['reset']}"


def show_intro(lines: list[str], delay: float = 0.2) -> None:
    """Clears the screen and prints the shared STB-SUITE ASCII banner.

    `lines` are centered inside '='*60 borders, one every `delay` seconds
    (pass delay=0 to disable the reveal animation).
    """
    os.system('cls' if os.name == 'nt' else 'clear')
    print(color_text(_LOGO, 'cyan'))
    print("\n" + "=" * 60)
    for line in lines:
        print(line.center(60))
        if delay:
            sleep(delay)
    print("=" * 60 + "\n")


def get_input(prompt: str, color: str = 'green') -> str:
    """Gets user input with a colored prompt (supports Tab-completion via readline)."""
    return input(color_text(prompt, color))


def get_float_input(prompt: str, default: float | None = None) -> float:
    """Gets a float from the user, retrying on invalid input."""
    while True:
        try:
            value_str = get_input(prompt)
            if value_str == "" and default is not None:
                return default
            return float(value_str)
        except ValueError:
            print(color_text("Please enter a valid number", 'red'))


def get_int_input(prompt: str, default: int | None = None) -> int:
    """Gets an integer from the user, retrying on invalid input."""
    while True:
        try:
            value_str = get_input(prompt)
            if value_str == "" and default is not None:
                return default
            return int(value_str)
        except ValueError:
            print(color_text("Please enter a valid integer", 'red'))


def print_dual(text: str, file_handle=None) -> None:
    """Prints to stdout with color, writes to file without color."""
    print(text)
    if file_handle:
        clean_text = re.sub(r'\x1b\[[0-9;]*m', '', text)
        file_handle.write(clean_text + "\n")


def print_section(label: str, f_out=None) -> None:
    """Prints one '[N] SECTION TITLE' header + separator, matching the
    suite's stage-report style (e.g. stb-ramanAnalysis's numbered sections).
    Moved here once stb-clean became a second consumer alongside
    stb-translate.
    """
    print_dual(f"\n{color_text(label, 'magenta')}", f_out)
    print_dual("-" * 60, f_out)


def print_table(headers, rows, f_out=None) -> None:
    """Prints a simple aligned ' | '-separated table via print_dual.

    `headers` is a list of column-name strings. `rows` is a list of
    (cells, color) tuples, where `cells` is a list of plain (uncolored)
    strings matching `headers` and `color` is an optional color name (e.g.
    'yellow') applied to the whole row, or None. Column widths are computed
    from the actual content, not fixed, so it reads well regardless of how
    long any one cell's text is (e.g. a full space-group label vs. a short
    number). Shared by core/structure_checks.py's per-check validation table
    and stb-mlrelax's before/after comparison table.
    """
    widths = [len(h) for h in headers]
    for cells, _ in rows:
        for i, cell in enumerate(cells):
            widths[i] = max(widths[i], len(cell))
    header_line = " | ".join(f"{h:<{w}}" for h, w in zip(headers, widths))
    print_dual(color_text(header_line, 'blue'), f_out)
    print_dual("-" * len(header_line), f_out)
    for cells, color in rows:
        line = " | ".join(f"{c:<{w}}" for c, w in zip(cells, widths))
        print_dual(color_text(line, color) if color else line, f_out)


@contextlib.contextmanager
def capture_library_noise(collector, label):
    """Captures stdout prints and warnings.warn() calls made by external
    libraries (MACE/torch/phonopy/spglib/pymatgen) during the wrapped block,
    instead of letting them interleave with a tool's own numbered report --
    appended to `collector` (a list of strings) and printed together, once,
    in a final LIBRARY WARNINGS section instead. A caller's own print_dual
    output never goes through here (only third-party calls are wrapped), so
    nothing from the report itself is ever captured/delayed. Moved here from
    phonons_create.py once xrdsearch.py became a second consumer of the
    exact same "collect noisy third-party output, report it once at the
    end" need.
    """
    buf = io.StringIO()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        with contextlib.redirect_stdout(buf):
            yield
    text = buf.getvalue().strip()
    if text:
        collector.append(f"[{label}]\n{text}")
    # De-duplicated by message text -- torch/mace commonly re-emit the exact
    # same DeprecationWarning once per call site internally (e.g. one per
    # torch.jit.load), which would otherwise print a dozen+ identical lines
    # here for a single underlying issue.
    seen = {}
    for w in caught:
        key = f"{w.category.__name__}: {w.message}"
        seen[key] = seen.get(key, 0) + 1
    for key, count in seen.items():
        suffix = f" (x{count})" if count > 1 else ""
        collector.append(f"[{label}] {key}{suffix}")


def run_with_spinner(func, *args, label: str = "Working", **kwargs):
    """Runs `func(*args, **kwargs)` in a background thread, showing an indeterminate
    CLI spinner with elapsed time while it blocks. For a slow call with no step-based
    progress callback to hook into -- e.g. icet's SQS search via pymatgen's
    SQSTransformation (stb-sqs), or pyxtal's XRD pattern-similarity calculation
    (stb-xrd --compare-to) -- this is elapsed-time feedback, not a real progress
    bar. Falls back to periodic plain-text lines when stderr isn't a terminal
    (e.g. output redirected to a file/log). Re-raises whatever `func` raised.
    """
    done = threading.Event()
    outcome = {}

    def worker():
        try:
            outcome["result"] = func(*args, **kwargs)
        except BaseException as exc:
            outcome["error"] = exc
        finally:
            done.set()

    thread = threading.Thread(target=worker, daemon=True)
    start = time.monotonic()
    thread.start()

    is_tty = sys.stderr.isatty()
    spinner_frames = "|/-\\"
    frame = 0
    last_line_at = 0.0
    try:
        while not done.wait(timeout=0.2):
            elapsed = time.monotonic() - start
            if is_tty:
                sys.stderr.write(f"\r  {spinner_frames[frame % len(spinner_frames)]} {label}... {elapsed:6.1f}s elapsed  ")
                sys.stderr.flush()
                frame += 1
            elif elapsed - last_line_at >= 5.0:
                sys.stderr.write(f"  ... {label} ({elapsed:.0f}s elapsed)\n")
                sys.stderr.flush()
                last_line_at = elapsed
    finally:
        if is_tty:
            sys.stderr.write("\r" + " " * 60 + "\r")
            sys.stderr.flush()

    if "error" in outcome:
        raise outcome["error"]
    return outcome["result"]
