"""Direct checks of core.hirshfeld's numerics, run by test.sh (section 6d).
Prints one 'OK <name>' line per passed check; any failure raises (exit 1).

- radial_density_profile on an analytic Gaussian (6 e-, sigma 0.6 Ang): the
  profile must hold the Gaussian's electrons with the default bins AND with
  bins far finer than the grid -- the old fixed-bin code left fine bins empty
  near the nucleus and read them as rho = 0 (5.13 of 6 e- on a real O atom);
- profile_radius: 99.9 % of a 3D Gaussian lies within ~4.03 sigma;
- estimate_remaining_error on a geometric history (known answer);
- unwrap_vacuum_axes on a molecule cut by the cell boundary.
"""
import numpy as np

from stb.core import hirshfeld as hf
from stb.hirshfeld_analysis import unwrap_vacuum_axes

box, n, sigma, n_e = 12.0, 96, 0.6, 6.0
lattice = np.eye(3) * box
x = np.arange(n) / n * box - box / 2
X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
rho = n_e / ((2 * np.pi) ** 1.5 * sigma ** 3) * np.exp(-(X**2 + Y**2 + Z**2) / (2 * sigma**2))
center = (0.5, 0.5, 0.5)

r, p = hf.radial_density_profile(rho, lattice, center)
n_default = hf.profile_electron_count(r, p)
print(f"default bins: N = {n_default:.4f} (exact {n_e})")
assert abs(n_default - n_e) < 0.01, n_default
assert r[0] == 0.0 and abs(p[0] - rho.max()) / rho.max() < 0.05, (r[0], p[0], rho.max())
print("OK profile electron count (default bins)")

r, p = hf.radial_density_profile(rho, lattice, center, n_bins=4000)
n_fine = hf.profile_electron_count(r, p)
print(f"4000 bins (~0.003 Ang, grid 0.125 Ang): N = {n_fine:.4f}")
assert abs(n_fine - n_e) < 0.02, n_fine
assert np.all(p > 0)
print("OK no empty bins (fine bins)")

r99 = hf.profile_radius(r, p)
print(f"r(99.9 %) = {r99:.3f} Ang = {r99 / sigma:.2f} sigma")
assert 3.6 * sigma < r99 < 4.4 * sigma, r99
assert hf.profile_radius(r, np.zeros_like(p)) is None
print("OK profile radius")

hist = [(0, None)] + [(k, 0.1 * 0.7 ** k) for k in range(1, 9)]
ratio, remaining = hf.estimate_remaining_error(hist)
assert abs(ratio - 0.7) < 1e-9, ratio
assert abs(remaining - 0.1 * 0.7 ** 8 * 0.7 / 0.3) < 1e-12, remaining
assert hf.estimate_remaining_error([(0, None), (1, 0.1)]) == (None, None)
ratio, remaining = hf.estimate_remaining_error([(0, None), (1, 0.1), (2, 0.2), (3, 0.4)])
assert ratio > 1 and remaining is None, (ratio, remaining)
print("OK remaining-error estimate")

frac = np.array([[0.5, 0.5, 0.98], [0.5, 0.5, 0.06]])   # a molecule cut by the c boundary
u = unwrap_vacuum_axes(frac, [False, False, True])
assert abs(abs(u[1, 2] - u[0, 2]) - 0.08) < 1e-9, u
assert np.allclose(u[:, :2], 0.5)
print("OK unwrap across the boundary")
