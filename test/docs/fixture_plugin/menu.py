"""Menu entries of the fictional test extension: workflow 5.99 with one stage (it has an example,
so a real guide) and tool 5.98 (no example: the site gives it a generated stub page)."""

from stb.core.menu import run_tool


def run_demo() -> None:
    run_tool("stb-demoext", [])


WORKFLOW_TOOLS = {
    98: {'title': "Demo Stub Tool (stb-demoext)",
         'description': "A fixture item with no example folder, so its guide page is a generated stub.",
         'func': run_demo},
    99: {'title': "Demo Extension Workflow",
         'description': "A fixture workflow that only exists to test the documentation build.",
         'stages': {
             1: {'title': "Stage 1 - Demo (stb-demoext)",
                 'description': "The only stage of the demo workflow.",
                 'func': run_demo},
         }},
}
