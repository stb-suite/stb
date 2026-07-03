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
