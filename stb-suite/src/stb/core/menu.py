"""Helpers the interactive menu shares with its extensions.

`run_tool` is the single way a menu entry runs a suite tool (it shells out to the
installed console command by name), and `prompt_pseudo_source` is the shared
pseudopotential-source prompt. They live here, not in `stb_suite.py`, so that an
extension package (see `stb_suite._load_plugins`) can use them without importing
the whole menu module.
"""

import os
import subprocess
from typing import List

from stb.core.cli import color_text, get_input
from stb.core.pseudopotentials import BANKS


def prompt_pseudo_source(optional: bool = True) -> str:
    """Shared pseudopotential-source prompt for every wrapper below that
    needs one (phonons, cohesive energy, input file, Hubbard U prep):
    a bundled bank (see core/pseudopotentials.py) or a custom path. Returns
    the raw string to pass straight through as the tool's -p/--pp-path/
    --pseudo-dir value (each tool resolves it itself); empty string only if
    `optional` and the user skips."""
    bank_list = list(BANKS.items())
    print(f"\n{color_text('Pseudopotential source:', 'yellow')}")
    for i, (name, info) in enumerate(bank_list, 1):
        print(f"  {color_text(str(i), 'cyan')} = Bundled: {info['description']} ({name})")
    print(f"  {color_text(str(len(bank_list) + 1), 'cyan')} = Custom path")
    prompt = f"Select option (1-{len(bank_list) + 1}"
    prompt += ", or Enter to skip): " if optional else "): "
    while True:
        choice = get_input(prompt).strip()
        if not choice and optional:
            return ""
        if choice.isdigit() and 1 <= int(choice) <= len(bank_list):
            return bank_list[int(choice) - 1][0]
        if choice == str(len(bank_list) + 1):
            path = os.path.expanduser(get_input("Custom pseudopotentials folder path: ").strip())
            if os.path.isdir(path):
                return path
            print(color_text(f"Path not found: '{path}'", 'red'))
            continue
        print(color_text("Invalid choice.", 'red'))


def prompt_model_card(model_path: str) -> str:
    """Asks for the model card (the mlff_model_card.json stb-mlffAnalysis writes next to
    each fine-tuned model) that must accompany a custom model file. Offers the one sitting
    next to the model as the default when it exists; asks again until a file is given."""
    from stb.core.model_card import CARD_FILE
    default = os.path.join(os.path.dirname(os.path.abspath(model_path)), CARD_FILE)
    default = default if os.path.isfile(default) else ""
    while True:
        suffix = f" [default: {default}]" if default else ""
        path = os.path.expanduser(get_input(f"  Model card (JSON) of this model{suffix}: ").strip()) or default
        if path and os.path.isfile(path):
            return path
        print(color_text("  Model card not found -- a custom model is always used with its card "
                         f"({CARD_FILE}, written by stb-mlffAnalysis).", 'red'))


def run_tool(tool_name: str, args: List[str], pause: bool = True) -> None:
    """Executes a suite tool as a subprocess.

    `pause=False` skips the "Press Enter to continue..." block -- for callers
    that invoke run_tool() several times in a row (e.g. run_strain_generator
    looping over multiple symmetry-equivalent directions in one go), so the
    user isn't interrupted after every individual subprocess call, only once
    at the very end.
    """
    try:
        cmd = [f"{tool_name}"] + args
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(color_text(f"\nError running {tool_name}: {e}", 'red'))
    except FileNotFoundError:
        print(color_text(f"\nTool {tool_name} not found!", 'red'))
        print(color_text(f"Make sure {tool_name} is in your system's PATH.", 'yellow'))
    if pause:
        input("\nPress Enter to continue...")
