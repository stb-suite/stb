"""Menu entries of the fictional test extension: workflow 4.99 with one stage."""

from stb.core.menu import run_tool


def run_demo() -> None:
    run_tool("stb-demoext", [])


WORKFLOW_TOOLS = {
    99: {'title': "Demo Extension Workflow",
         'description': "A fixture workflow that only exists to test the documentation build.",
         'stages': {
             1: {'title': "Stage 1 - Demo (stb-demoext)",
                 'description': "The only stage of the demo workflow.",
                 'func': run_demo},
         }},
}
