"""Figures of the 3.3 example (stb-cube), every one drawn from the CUBE files stb-cube writes (read back with ASE),
so they also show that the cubes hold what the report says.

    python3 theory_figures.py        # figures/mesh_cutoff.png, si_density.png, h2o.png, o2_spin.png (~20 s)
"""
import os
import subprocess
import tempfile
import shutil
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import map_coordinates
from ase.io.cube import read_cube_data

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
BOHR = 0.529177210903


def cube(work, run, gtype, spin="total"):
    """Runs stb-cube on siesta_runs/<run> (copied into `work`) and returns (data, atoms) of the cube."""
    out = f"{run}_{gtype}_{spin}.cube"
    subprocess.run(["stb-cube", "--no-intro", "-l", f"{run}/{run}", "-t", gtype, "--spin", spin, "-o", out],
                   cwd=work, check=True, stdout=subprocess.DEVNULL)
    return read_cube_data(os.path.join(work, out))


def plane(data, cell, origin, u, v, nu=240, nv=240):
    """Values on the plane origin + s u + t v (s, t in [0, 1]), periodic trilinear interpolation of the cube."""
    s, t = np.meshgrid(np.linspace(0, 1, nu), np.linspace(0, 1, nv), indexing="ij")
    pts = origin + s[..., None] * u + t[..., None] * v
    frac = pts.reshape(-1, 3) @ np.linalg.inv(cell)
    idx = (frac * np.array(data.shape)).T
    return map_coordinates(data, idx, order=1, mode="grid-wrap").reshape(nu, nv)


# --------------------------------------------------------------------------------- SIESTA's mesh vs MeshCutoff --
def smooth(n):
    for p in (2, 3, 5):
        while n % p == 0:
            n //= p
    return n == 1


def siesta_mesh(cutoff, cell):
    """n_i: the smallest EVEN number with factors 2, 3, 5 only, at least 2 sqrt(Ecut)/|b_i| (b in 1/Bohr)."""
    b = np.linalg.norm(2 * np.pi * np.linalg.inv(cell / BOHR).T, axis=1)
    out = []
    for bi in b:
        n = int(np.ceil(2 * np.sqrt(cutoff) / bi - 1e-9))
        while n % 2 or not smooth(n):
            n += 1
        out.append(n)
    return np.array(out), b


def used_cutoff(n, b):
    return float(np.min((n * b / 2) ** 2))


runs = {"si": (np.array([[0, .5, .5], [.5, 0, .5], [.5, .5, 0]]) * 5.476, 200.0, 30, 248.851),
        "h2o": (np.eye(3) * 9.0, 100.0, 60, 122.834),
        "o2": (np.eye(3) * 8.0, 120.0, 54, None)}
for name, (cell, req, n_log, used_log) in runs.items():          # the rule must reproduce the three real runs
    n, b = siesta_mesh(req, cell)
    assert n[0] == n_log, (name, n, n_log)
    if used_log:
        assert abs(used_cutoff(n, b) - used_log) < 1e-2, (name, used_cutoff(n, b), used_log)

os.makedirs(FIG, exist_ok=True)
cell_si = runs["si"][0]
req = np.linspace(50, 400, 1401)
used, npts = [], []
for e in req:
    n, b = siesta_mesh(e, cell_si)
    used.append(used_cutoff(n, b))
    npts.append(n[0])
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
ax[0].plot(req, used, color="#1f6f8b", lw=1.6, label="used (effective)")
ax[0].plot(req, req, "k--", lw=0.8, label="required = used")
ax[0].plot([200], [248.851], "o", color="#c0392b", label="this example's Si run (log)")
ax[0].set_xlabel("Mesh.Cutoff required (Ry)")
ax[0].set_ylabel("cutoff the mesh really has (Ry)")
ax[0].legend(frameon=False, fontsize=8)
ax[0].set_title("(a) SIESTA rounds the mesh up")
ax[1].step(req, npts, where="post", color="#1f6f8b", lw=1.6)
ax[1].plot([200], [30], "o", color="#c0392b")
ax[1].set_xlabel("Mesh.Cutoff required (Ry)")
ax[1].set_ylabel("points per cell vector")
ax[1].set_title("(b) n: even, factors 2, 3, 5 only")
for a in ax:
    a.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "mesh_cutoff.png"), dpi=150)
plt.close(fig)

work = tempfile.mkdtemp()
try:
    for r in ("si", "h2o", "o2"):
        shutil.copytree(os.path.join(HERE, "siesta_runs", r), os.path.join(work, r))

    # ------------------------------------------------------------------- Si: density and deformation density --
    a = 5.476
    u = np.array([1.0, 1.0, 0.0]) * a          # [110], one conventional face diagonal
    v = np.array([0.0, 0.0, 1.0]) * a          # [001]
    # diamond atoms (conventional cell, translated) that lie on the (1-10) plane x = y, as (along u, along v)
    basis = [(0, 0, 0), (0, .5, .5), (.5, 0, .5), (.5, .5, 0), (.25, .25, .25), (.25, .75, .75), (.75, .25, .75),
             (.75, .75, .25)]
    IN_PLANE = sorted({(round((x + i) * np.sqrt(2) * a, 3), round((z + k) * a, 3))
                       for x, y, z in basis for i in (0, 1) for k in (0, 1)
                       if abs(x - y) < 1e-9 and x + i <= 1 and z + k <= 1})
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.3))
    for k, (gtype, title, cmap) in enumerate([("RHO", "(a) RHO: valence density", "viridis"),
                                              ("DRHO", "(b) DRHO: density - free atoms", "RdBu_r")]):
        data, atoms = cube(work, "si", gtype)
        p = plane(data, atoms.cell.array, np.zeros(3), u, v)
        lim = np.abs(p).max()
        im = ax[k].imshow(p.T, origin="lower", extent=[0, np.linalg.norm(u), 0, a], cmap=cmap,
                          vmin=-lim if gtype == "DRHO" else 0, vmax=lim)
        ax[k].contour(np.linspace(0, np.linalg.norm(u), p.shape[0]), np.linspace(0, a, p.shape[1]), p.T,
                      levels=8, colors="k", linewidths=0.4)
        for x, y in IN_PLANE:
            ax[k].plot(x, y, "o", ms=5, mfc="white", mec="k")
        ax[k].set_title(title)
        ax[k].set_xlabel("along [110] (Ang)")
        ax[k].set_ylabel("along [001] (Ang)")
        fig.colorbar(im, ax=ax[k], label="e/Bohr$^3$ (as in the cube)")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "si_density.png"), dpi=150)
    plt.close(fig)

    # ------------------------------------------------------------ water: DRHO in its plane, VH plane averages --
    data, atoms = cube(work, "h2o", "DRHO")
    cell = atoms.cell.array
    p = plane(data, cell, np.array([4.5, 2.0, 2.0]), np.array([0, 5.0, 0]), np.array([0, 0, 5.0]))
    vh, _ = cube(work, "h2o", "VH")
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    lim = 0.02
    im = ax[0].imshow(p.T, origin="lower", extent=[2, 7, 2, 7], cmap="RdBu_r", vmin=-lim, vmax=lim)
    for at in atoms:
        ax[0].plot(at.position[1], at.position[2], "o", ms=7 if at.symbol == "O" else 5, mfc="white", mec="k")
    ax[0].set_xlabel("y (Ang)")
    ax[0].set_ylabel("z (Ang)")
    ax[0].set_title("(a) DRHO in the molecular plane")
    fig.colorbar(im, ax=ax[0], label="e/Bohr$^3$ (clipped at $\\pm$0.02)")
    z = np.arange(vh.shape[2]) / vh.shape[2] * cell[2, 2]
    y = np.arange(vh.shape[1]) / vh.shape[1] * cell[1, 1]
    ax[1].plot(z, vh.mean(axis=(0, 1)), color="#c0392b", label="along z (the dipole axis)")
    ax[1].plot(y, vh.mean(axis=(0, 2)), color="#1f6f8b", label="along y")
    for at in atoms:
        ax[1].axvline(at.position[2], color="#c0392b", lw=0.5, ls=":")
    ax[1].set_ylim(-0.6, 0.3)
    ax[1].set_xlabel("position (Ang)")
    ax[1].set_ylabel("plane average of VH (eV)")
    ax[1].set_title("(b) VH: flat along y, sloped along z")
    ax[1].legend(frameon=False, fontsize=8)
    ax[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "h2o.png"), dpi=150)
    plt.close(fig)

    # --------------------------------------------------------------------------- O2: spin density, up and down --
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.9))
    for k, (spin, title) in enumerate([("up", "(a) up"), ("down", "(b) down"), ("diff", "(c) up - down")]):
        data, atoms = cube(work, "o2", "RHO", spin)
        p = plane(data, atoms.cell.array, np.array([1.5, 4.0, 1.5]), np.array([5.0, 0, 0]), np.array([0, 0, 5.0]))
        lim = 0.3 if spin != "diff" else 0.1
        im = ax[k].imshow(p.T, origin="lower", extent=[1.5, 6.5, 1.5, 6.5],
                          cmap="RdBu_r" if spin == "diff" else "viridis", vmin=-lim if spin == "diff" else 0,
                          vmax=lim)
        for at in atoms:
            ax[k].plot(at.position[0], at.position[2], "o", ms=6, mfc="white", mec="k")
        ax[k].set_title(title)
        ax[k].set_xlabel("x (Ang)")
        ax[k].set_ylabel("z (Ang)")
        fig.colorbar(im, ax=ax[k], label="e/Bohr$^3$ (clipped)")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "o2_spin.png"), dpi=150)
    plt.close(fig)
finally:
    shutil.rmtree(work)
print("figures/mesh_cutoff.png, si_density.png, h2o.png, o2_spin.png")
