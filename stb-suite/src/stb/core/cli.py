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
import re
import sys
import time
import threading
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
    """Prints to stdout with color, writes to file without color.

    Moved here from phonons_pos.py when stb-convergenceAnalysis needed the
    same helper; phonons_create.py, phonons_qha.py, phonons_ml.py and
    elastic_analysis.py still carry their own identical, undeduplicated
    copies of this function -- pre-existing before this move, left alone here
    since they belong to unrelated workflows.
    """
    print(text)
    if file_handle:
        clean_text = re.sub(r'\x1b\[[0-9;]*m', '', text)
        file_handle.write(clean_text + "\n")


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
