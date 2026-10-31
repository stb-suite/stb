"""Figures of the 3.4 example (stb-siesta2wtb), drawn from the tb-*.ham files the tool writes:

  figures/hbn_bands.png   hBN: bands rebuilt from tb-NP.ham against SIESTA's own hbn.bands
  figures/fe_bands.png    Fe: up and down bands rebuilt from tb-SP.ham against fe.bands
  figures/precision.png   how much the file's 6 decimals change each eigenvalue, against its energy

    python3 theory_figures.py        # ~1 min (writes the two tb files in a scratch folder)
"""
import os
import shutil
import subprocess
import tempfile
import numpy as np
import scipy.linalg as sl
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
BOHR = 0.529177210903


def read_ham(path):
    with open(path) as f:
        head = [next(f) for _ in range(9)]
    cell = np.array([[float(x) for x in head[i].split()] for i in (3, 4, 5)])
    try:
        import pandas as pd
        data = pd.read_csv(path, skiprows=9, sep=r"\s+", header=None).to_numpy()
    except ImportError:
        data = np.loadtxt(path, skiprows=9)
    return float(head[2]), cell, int(head[6]), data


def hk_sk(data, dim, cell, kcart):
    phase = np.exp(1j * (data[:, :3] @ kcart))
    i, j = data[:, 3].astype(int) - 1, data[:, 4].astype(int) - 1
    hk, sk = np.zeros((dim, dim), complex), np.zeros((dim, dim), complex)
    np.add.at(hk, (i, j), (data[:, 5] + 1j * data[:, 6]) * phase)
    np.add.at(sk, (i, j), data[:, 7] * phase)
    return (hk + hk.conj().T) / 2, (sk + sk.conj().T) / 2


def read_bands(path):
    """SIESTA <label>.bands: E_F, path length x (1/Bohr), energies [nk, nspin, nbands] (eV)."""
    t = open(path).read().split()
    ef = float(t[0])
    nb, ns, nk = int(t[5]), int(t[6]), int(t[7])
    v = np.array(t[8:8 + nk * (1 + nb * ns)], float).reshape(nk, 1 + nb * ns)
    return ef, v[:, 0], v[:, 1:].reshape(nk, ns, nb)


def path_kpoints(cell, corners, x):
    """Cartesian k (1/Ang) at each path length x (1/Bohr) along the polyline through `corners` (reduced)."""
    b = 2 * np.pi * np.linalg.inv(cell).T                 # 1/Ang
    pts = np.array(corners) @ b
    seg = np.linalg.norm(np.diff(pts, axis=0), axis=1) * BOHR    # lengths in 1/Bohr, as SIESTA's x
    cum = np.concatenate([[0], np.cumsum(seg)])
    out = []
    for xi in x:
        s = min(np.searchsorted(cum, xi, side="right") - 1, len(seg) - 1)
        f = (xi - cum[s]) / seg[s] if seg[s] > 0 else 0.0
        out.append(pts[s] + f * (pts[s + 1] - pts[s]))
    return np.array(out), cum


def run_tool(work, run):
    subprocess.run(["stb-siesta2wtb", "-i", "calc.fdf", "-o", "calc.out", "--no-intro", "--no-check"],
                   cwd=os.path.join(work, run), check=True, stdout=subprocess.DEVNULL)


os.makedirs(FIG, exist_ok=True)
work = tempfile.mkdtemp()
try:
    for r in ("hbn", "fe"):
        shutil.copytree(os.path.join(HERE, "siesta_runs", r), os.path.join(work, r))
        run_tool(work, r)
    summary = []

    # ------------------------------------------------------------------------------------------------- hBN bands --
    ef_h, cell, dim, data = read_ham(os.path.join(work, "hbn", "tb-NP.ham"))
    ef, x, eb = read_bands(os.path.join(HERE, "siesta_runs", "hbn", "hbn.bands"))
    corners = [[0, 0, 0], [0.5, 0, 0], [1 / 3, 1 / 3, 0], [0, 0, 0]]
    kc, cum = path_kpoints(cell, corners, x)
    ef_file = np.array([sl.eigh(*hk_sk(data, dim, cell, k), eigvals_only=True) for k in kc])   # referred to E_F
    dev = np.abs(ef_file - (eb[:, 0, :] - ef))
    near = np.abs(eb[:, 0, :] - ef) < 5
    summary.append(f"hBN: max |dE| along the path within 5 eV of E_F {dev[near].max() * 1000:.3f} meV")
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    for n in range(eb.shape[2]):
        ax.plot(x, eb[:, 0, n] - ef, color="#1f6f8b", lw=1.4, label="SIESTA hbn.bands" if n == 0 else None)
    for n in range(dim):
        ax.plot(x[::4], ef_file[::4, n], "o", ms=2.6, mfc="none", mec="#c0392b",
                label="rebuilt from tb-NP.ham" if n == 0 else None)
    for c in cum:
        ax.axvline(c, color="0.7", lw=0.6)
    ax.set_xticks(cum)
    ax.set_xticklabels(["Γ", "M", "K", "Γ"])
    ax.axhline(0, color="0.4", ls="--", lw=0.6)
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(-22, 12)
    ax.set_ylabel("E − E$_F$ (eV)")
    ax.set_title("hBN monolayer: SIESTA vs the WanTiBEXOS file")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "hbn_bands.png"), dpi=150)
    plt.close(fig)

    # --------------------------------------------------------------------------------------------- Fe up / down --
    ef_h, cell, dim, data = read_ham(os.path.join(work, "fe", "tb-SP.ham"))
    no = dim // 2
    ef, x, eb = read_bands(os.path.join(HERE, "siesta_runs", "fe", "fe.bands"))
    corners = [[0, 0, 0], [0.5, 0.5, -0.5], [0, 0, 0.5], [0, 0, 0], [0.25, 0.25, 0.25]]
    kc, cum = path_kpoints(cell, corners, x)
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    worst = 0.0
    for spin, (color, name) in enumerate([("#c0392b", "up"), ("#1f6f8b", "down")]):
        blk = slice(spin * no, (spin + 1) * no)
        e_file = []
        for k in kc:
            hk, sk = hk_sk(data, dim, cell, k)
            e_file.append(sl.eigh(hk[blk, blk], sk[blk, blk], eigvals_only=True))
        e_file = np.array(e_file)
        ref = eb[:, spin, :] - ef
        near = np.abs(ref) < 5
        worst = max(worst, np.abs(e_file - ref)[near].max())
        for n in range(eb.shape[2]):
            ax.plot(x, ref[:, n], color=color, lw=1.3, label=f"SIESTA, spin {name}" if n == 0 else None)
        for n in range(no):
            ax.plot(x[::4], e_file[::4, n], "o", ms=2.4, mfc="none", mec="k",
                    label="rebuilt from tb-SP.ham" if (n == 0 and spin == 0) else None)
    summary.append(f"Fe: max |dE| along the path within 5 eV of E_F {worst * 1000:.3f} meV")
    for c in cum:
        ax.axvline(c, color="0.7", lw=0.6)
    ax.set_xticks(cum)
    ax.set_xticklabels(["Γ", "H", "N", "Γ", "P"])
    ax.axhline(0, color="0.4", ls="--", lw=0.6)
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(-9, 6)
    ax.set_ylabel("E − E$_F$ (eV)")
    ax.set_title("bcc Fe: up and down blocks of tb-SP.ham")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fe_bands.png"), dpi=150)
    plt.close(fig)

    # ----------------------------------------------------------------------------- precision of the 6 decimals --
    import sisl
    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    rng = np.random.default_rng(1)
    for run, code, color, label in (("hbn", "NP", "#1f6f8b", "hBN (NP)"), ("fe", "SP", "#c0392b", "Fe (SP)")):
        _, cell, dim, data = read_ham(os.path.join(work, run, f"tb-{code}.ham"))
        ham = sisl.get_sile(os.path.join(HERE, "siesta_runs", run, "calc.fdf")).read_hamiltonian()
        es, ds = [], []
        for _ in range(12):
            kf = rng.random(3)
            if run == "hbn":
                kf[2] = 0.0
            e = np.sort(sl.eigh(*hk_sk(data, dim, cell, kf @ (2 * np.pi * np.linalg.inv(cell).T)),
                                eigvals_only=True))
            ref = (np.sort(np.concatenate([ham.eigh(k=kf, spin=0), ham.eigh(k=kf, spin=1)]))
                   if ham.spin.is_polarized else np.sort(ham.eigh(k=kf)))
            es.append(ref)
            ds.append(np.abs(e - ref))
        es, ds = np.concatenate(es), np.concatenate(ds)
        ax.semilogy(es, np.maximum(ds * 1000, 1e-6), "o", ms=2.5, color=color, alpha=0.6, label=label)
    ax.axvspan(-5, 5, color="0.9", zorder=0, label="±5 eV of E$_F$")
    ax.set_xlabel("E − E$_F$ (eV), sisl on the SIESTA H, S")
    ax.set_ylabel("|E(file) − E(sisl)| (meV)")
    ax.set_title("What the file's 6 decimals change (12 random k each)")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "precision.png"), dpi=150)
    plt.close(fig)
finally:
    shutil.rmtree(work)
print("\n".join(summary))
print("figures/hbn_bands.png, fe_bands.png, precision.png")
