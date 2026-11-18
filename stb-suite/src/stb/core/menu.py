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
from stb.core.pseudopotentials import BANKS, bank_info, list_bank_elements, short_bank_name


def _choose(title, options, default=1):
    """Numbered choice among `options` [(value, text)]; Enter takes `default` (1-based)."""
    print(f"\n  {color_text(title, 'yellow')}")
    for i, (_value, text) in enumerate(options, 1):
        print(f"    {color_text(str(i), 'cyan')} = {text}{' (default)' if i == default else ''}")
    while True:
        choice = get_input(f"  Select (1-{len(options)}) [default: {default}]: ").strip() or str(default)
        if choice.isdigit() and 1 <= int(choice) <= len(options):
            return options[int(choice) - 1][0]
        print(color_text("  Invalid choice.", 'red'))


def prompt_dojo_table() -> str:
    """Asks relativity / functional / accuracy and returns the PseudoDojo bank name in its
    shortest form (e.g. 'dojo-pbesol'; 'dojo' for the defaults); only combinations that exist
    are offered."""
    rel = _choose("Relativity:", [
        ("sr", "Scalar-relativistic"),
        ("fr", "Full-relativistic -- required for 'Spin spin-orbit'"),
        ("sr3plus", "Scalar, lanthanides with 4f in the core (Ce-Lu, Ln3+ compounds only)"),
    ])
    if rel == "sr3plus":
        print(color_text("  The 4f-in-core table is PBE standard and holds only Ce-Lu: give the other elements "
                         "a fallback (stb-pseudo --fallback-dir dojo), or a folder holding both.", 'yellow'))
        xc, acc = "pbe", "standard"
    else:
        xc_options = [("pbe", "PBE"), ("pbesol", "PBEsol")]
        if rel == "sr":
            xc_options.append(("pw", "LDA (PW92)"))
        xc = _choose("Functional (must match XC.Functional / XC.Authors of the calculation):", xc_options)
        acc = _choose("Accuracy:", [
            ("standard", "standard"),
            ("stringent", "stringent -- more semicore: groups 13-16 add the s,p shell below the d "
                          "(Ge/Sn z 14->22); Hf/Ta add 4f (z 12->26, 13->27). More expensive."),
        ])
    name = short_bank_name(f"dojo-{rel}-{xc}-{acc}")
    info = bank_info(name)
    print(color_text(f"  -> {name} = {info['name']}: {info['description']}, "
                     f"{len(list_bank_elements(name))} elements", 'green'))
    return name


def prompt_pseudo_source(optional: bool = True) -> str:
    """Shared pseudopotential-source prompt for every wrapper that needs one: a
    bundled bank (see core/pseudopotentials.py; the PseudoDojo bank asks for
    relativity, functional and accuracy) or a custom path. Returns the raw string
    to pass straight through as the tool's -p/--pp-path/--pseudo-dir value (each
    tool resolves it itself); empty string only if `optional` and the user skips."""
    print(f"\n{color_text('Pseudopotential source:', 'yellow')}")
    print(f"  {color_text('1', 'cyan')} = Bundled: PseudoDojo (PSML) -- choose relativity / functional / accuracy")
    print(f"  {color_text('2', 'cyan')} = Bundled: {BANKS['virtual_vault']['description']} (virtual_vault)")
    print(f"  {color_text('3', 'cyan')} = Custom path")
    prompt = "Select option (1-3"
    prompt += ", or Enter to skip): " if optional else "): "
    while True:
        choice = get_input(prompt).strip()
        if not choice and optional:
            return ""
        if choice == "1":
            return prompt_dojo_table()
        if choice == "2":
            return "virtual_vault"
        if choice == "3":
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
