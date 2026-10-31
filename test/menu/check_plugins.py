#!/usr/bin/env python3
"""Checks stb_suite._load_plugins, the hook through which an extension package
adds entries to the interactive menu (entry-point group "stb.plugins").

The real entry points are replaced by fake ones, so no package needs installing:
a plugin's entries are merged in and kept in numeric order (the menu lists them in
insertion order); a number already taken keeps the original entry and warns; a
plugin that fails to load is skipped with a warning; none of this raises. Also
checks that `run_tool` and `prompt_pseudo_source` live in stb.core.menu (the API an
extension imports) and that the menu module still re-exports them.

Usage: check_plugins.py   (exits 1 on any failure). Needs the stb-suite package.
"""
import contextlib
import importlib.metadata
import io
import sys
from types import SimpleNamespace

failures = []


def check(condition, message):
    print(f"[ OK ] {message}" if condition else f"[FAIL] {message}")
    if not condition:
        failures.append(message)


def entries(number, title):
    return {number: {"title": title, "description": "fixture",
                     "stages": {1: {"title": "Stage 1", "description": "fixture", "func": lambda: None}}}}


class FakeEntryPoints(list):
    def select(self, group):
        return FakeEntryPoints(self if group == "stb.plugins" else [])


def main():
    import stb.stb_suite as menu
    from stb.core import menu as core_menu

    check(menu.run_tool is core_menu.run_tool and menu.prompt_pseudo_source is core_menu.prompt_pseudo_source,
          "the menu re-exports run_tool and prompt_pseudo_source from stb.core.menu")

    taken = 1                                    # 5.1 exists in every install
    original_title = menu.WORKFLOW_TOOLS[taken]["title"]

    def boom():
        raise RuntimeError("boom")

    plugins = FakeEntryPoints([
        SimpleNamespace(name="demo", load=lambda: SimpleNamespace(WORKFLOW_TOOLS=entries(90, "Demo A"),
                                                          CATALYSIS_TOOLS={**entries(2, "Demo C"), **entries(1, "Demo B")})),
        SimpleNamespace(name="clash", load=lambda: SimpleNamespace(WORKFLOW_TOOLS=entries(taken, "Intruder"))),
        SimpleNamespace(name="broken", load=boom),
    ])
    real = importlib.metadata.entry_points
    importlib.metadata.entry_points = lambda *args, **kwargs: plugins
    output = io.StringIO()
    try:
        with contextlib.redirect_stdout(output):
            menu._load_plugins()
    except Exception as error:
        check(False, f"_load_plugins raised {error!r}")
    finally:
        importlib.metadata.entry_points = real
    text = output.getvalue()

    try:
        check(90 in menu.WORKFLOW_TOOLS and {1, 2} <= set(menu.CATALYSIS_TOOLS), "a plugin's entries are merged into the menu")
        keys = list(menu.CATALYSIS_TOOLS)
        check(keys == sorted(keys), "the menu is kept in numeric order (1 before 2, though the plugin gave 2 first)")
        check(menu.WORKFLOW_TOOLS[taken]["title"] == original_title, "a number already taken keeps the original entry")
        check("already taken" in text and "clash" in text, "a taken number produces a warning naming the plugin")
        check("could not be loaded" in text and "broken" in text, "a plugin that fails to load produces a warning")
        codes = menu._flatten_tool_codes()
        check("6.1.1" in codes and callable(codes["6.1.1"]), "a merged entry is reachable by its dotted code (6.1.1)")
        check("5.90.1" in codes, "an entry merged into an existing section is reachable too (5.90.1)")
    finally:
        menu.WORKFLOW_TOOLS.pop(90, None)        # leave the module as it was found
        for number in (1, 2):
            menu.CATALYSIS_TOOLS.pop(number, None)

    print(f"{len(failures)} failure(s)")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
