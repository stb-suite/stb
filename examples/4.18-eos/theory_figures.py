#!/usr/bin/env python3
"""Draws figures/eos_forms.png of README section 1.4: the energy, pressure and bulk modulus of the
equation-of-state forms stb-eosAnalysis can fit (Murnaghan, Birch-Murnaghan, Vinet), for the
parameters of bulk Si from the all-electron Delta-test reference (V0 = 20.453 A^3/atom,
B0 = 88.545 GPa, B0' = 4.31), compared with the plain harmonic approximation. The grey band is the
+/-8 % volume window of the example's scan: inside it the forms are indistinguishable, which is why
a narrow scan is enough for V0, E0 and B0 and not for B0' or for high pressures.

    python3 theory_figures.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ase.eos import birchmurnaghan, murnaghan, vinet

GPA = 160.21766208            # eV/A^3 -> GPa
V0, B0, BP = 20.4530, 88.545 / GPA, 4.31
V = np.linspace(0.62 * V0, 1.45 * V0, 800)


def harmonic(v):
    return 0.5 * B0 * (v - V0) ** 2 / V0


forms = {"harmonic (B0 only)": harmonic,
         "Murnaghan": lambda v: murnaghan(v, 0.0, B0, BP, V0),
         "Birch–Murnaghan (3rd order)": lambda v: birchmurnaghan(v, 0.0, B0, BP, V0),
         "Vinet": lambda v: vinet(v, 0.0, B0, BP, V0)}
style = {"harmonic (B0 only)": ("#888888", ":"), "Murnaghan": ("#CC79A7", "--"),
         "Birch–Murnaghan (3rd order)": ("#D55E00", "-"), "Vinet": ("#009E73", "-.")}

fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
for name, f in forms.items():
    e = f(V)
    p = -np.gradient(e, V) * GPA
    b = -V * np.gradient(p, V)
    color, ls = style[name]
    ax[0].plot(V / V0, 1000 * e, color=color, ls=ls, lw=2, label=name)
    ax[1].plot(V / V0, p, color=color, ls=ls, lw=2, label=name)
    if not name.startswith("harmonic"):
        ax[2].plot(p, b, color=color, ls=ls, lw=2, label=name)
for a in ax[:2]:
    a.axvspan(0.92, 1.08, color="#cccccc", alpha=0.5, label="±8 % scan window")
    a.set_xlabel("V / V0")
ax[0].set_ylim(-20, 700)
ax[0].set_ylabel("E − E0 (meV/atom)")
ax[0].set_title("Energy: the forms agree near V0")
ax[1].set_ylim(-12, 60)
ax[1].set_ylabel("P = −dE/dV (GPa)")
ax[1].axhline(0, color="k", lw=0.8)
ax[1].set_title("Pressure")
ax[2].set_xlim(-6, 60)
ax[2].set_ylim(60, 400)
ax[2].set_xlabel("P (GPa)")
ax[2].set_ylabel("B = −V dP/dV (GPa)")
ax[2].set_title("Bulk modulus along the curve (slope at P=0 is B0')")
ax[0].legend(fontsize=8, loc="upper right")
fig.suptitle("Equation-of-state forms for Si (B0 = 88.5 GPa, B0' = 4.31, V0 = 20.453 Å³/atom)", fontsize=13)
fig.tight_layout()
fig.savefig("figures/eos_forms.png", dpi=130)
print("figures/eos_forms.png written")
