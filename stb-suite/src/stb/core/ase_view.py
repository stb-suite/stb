"""Shared interactive 3D structure viewer (ASE), used by stb-translate and
stb-inputfile.
"""

from stb.core.cli import color_text


def view_structure_interactive(atoms):
    """Opens ASE's interactive 3D structure viewer. Needs a display (local X11,
    or `ssh -X`/`-Y` to a remote machine) and a working Tk installation."""
    try:
        from ase.visualize import view
        # block=True: waits for the GUI subprocess so a failure (e.g. no display)
        # raises here and is reported, instead of failing silently in the background.
        view(atoms, block=True)
    except Exception as e:
        print(color_text(f"[FAIL] Could not open the interactive 3D viewer: {e}", 'red'))
        print(color_text(
            "       This needs a display (local X11, or `ssh -X`/`-Y` to a remote "
            "machine) and a working Tk installation.", 'yellow'))
