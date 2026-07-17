"""Generalized Quasi-Chemical Approximation (GQCA) mass-action solver for
substitutional alloy configurational thermodynamics.

Physics: A. Sher, M. van Schilfgaarde, A.-B. Chen, W. Chen, Phys. Rev. B 36,
4279 (1987); A.-B. Chen and A. Sher, "Semiconductor Alloys" (Plenum Press,
1995). An alloy A(1-x)Bx is treated as an ensemble of N-site clusters (here
N=2, the classic Guggenheim pair/"bond" approximation -- see stb-gqca's own
docstring for why larger clusters are out of scope), each in one of M=N+1
compositions j=0..N (j = number of B atoms in the cluster), statistically
independent "molecules" in a mass-action equilibrium.

Deliberately independent of icet -- the physics doesn't need it, keeping
this module light and unit-testable without the heavy optional import chain
sqs.py/gqca.py carry. Written generally in N (not hardcoded to N=2) since
the math is no simpler for a fixed N -- leaves the door open for larger
clusters later without rewriting the solver.

The mass-action law and mixing-energy reference (Delta_j interpolated in
the LOCAL cluster composition n_j/N, not the target alloy composition x)
are confirmed against Ferreira et al., Mater. Today Phys. 2024 ("Ab initio
modeling of superconducting alloys", Extended GQCA), which restates the
original Sher et al. 1987 formalism (their Eq. 3, 4, 13-15) before
extending it. The exact configurational entropy expression could NOT be
pinned to a primary-source closed form during this tool's development --
see mixing_entropy's own docstring, flagged [UNVERIFIED] there.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.optimize import brentq

_BOLTZMANN_EV_K = 8.617333262e-5  # eV/K


def degeneracies(N: int) -> list[int]:
    """g_j = C(N, j) for j=0..N -- Sher et al., PRB 36, 4279 (1987). For the
    pair cluster (N=2) this is trivially [1, 2, 1] (AA, AB, BB), since a
    single bond's own 2-fold symmetry makes both "arrangements" of one A
    and one B atom identical -- no ambiguity, unlike larger clusters (see
    module docstring).
    """
    return [math.comb(N, j) for j in range(N + 1)]


def mixing_energies(e: np.ndarray, N: int) -> np.ndarray:
    """Delta_j = e_j - (1 - n_j/N)*e_0 - (n_j/N)*e_N -- the excess/mixing
    energy of cluster composition j relative to a LINEAR interpolation, in
    the LOCAL cluster composition n_j/N (not the global target alloy
    composition x), between the two pure end-member cluster energies e_0
    (n_j=0) and e_N (n_j=N). CONFIRMED against Ferreira et al. 2024 Eq. 4
    ("Delta_j = H_j - (n-n_j)/n * H_A - n_j/n * H_B", n there == N here).

    `e` must have length N+1, ordered by increasing n_j (e[0] is the pure-A
    cluster, e[N] is the pure-B cluster).
    """
    e = np.asarray(e, dtype=float)
    if len(e) != N + 1:
        raise ValueError(f"Expected {N + 1} cluster energies (N={N}), got {len(e)}.")
    n = np.arange(N + 1, dtype=float)
    return e - (1.0 - n / N) * e[0] - (n / N) * e[N]


def cluster_fractions(eta: float, T: float, g: list[int], delta_e: np.ndarray, N: int) -> np.ndarray:
    """x_j(eta,T) = g_j * eta^n_j * exp(-Delta_j / kB*T) / Z, the mass-action
    equilibrium cluster-composition distribution (CONFIRMED against
    Ferreira et al. 2024 Eq. 13, same functional form) -- exposed
    separately from solve_eta so callers can inspect x_j(eta) at a trial
    eta, not just the self-consistent solution.

    eta=0 and eta=inf are handled directly (return the pure-A/pure-B
    limits [1,0,...,0]/[0,...,0,1]) rather than evaluating 0**0 or
    overflowing exp(), matching solve_eta's own x_target in {0,1}
    special-casing. T must be strictly positive (a physical temperature --
    the CLI never passes T<=0; the T->0 limit of the mass-action law is a
    genuinely different, unused degenerate case and isn't handled here).
    """
    if T <= 0:
        raise ValueError(f"T must be a positive temperature in Kelvin, got {T}.")
    g = np.asarray(g, dtype=float)
    delta_e = np.asarray(delta_e, dtype=float)
    n = np.arange(N + 1, dtype=float)
    if eta <= 0.0:
        result = np.zeros(N + 1)
        result[0] = 1.0
        return result
    if math.isinf(eta):
        result = np.zeros(N + 1)
        result[N] = 1.0
        return result
    kT = _BOLTZMANN_EV_K * T
    log_weights = np.log(g) + n * np.log(eta) - delta_e / kT
    weights = np.exp(log_weights - log_weights.max())
    return weights / weights.sum()


def _mean_composition(eta: float, T: float, g: list[int], delta_e: np.ndarray, N: int) -> float:
    """Sum_j (n_j/N)*x_j(eta,T) -- the mean alloy composition a given eta
    (at fixed T) implies. solve_eta finds the eta making this equal the
    caller's x_target.
    """
    x_j = cluster_fractions(eta, T, g, delta_e, N)
    n = np.arange(N + 1, dtype=float)
    return float(np.sum((n / N) * x_j))


def solve_eta(x_target: float, T: float, g: list[int], delta_e: np.ndarray, N: int) -> float:
    """Root-finds eta (scipy.optimize.brentq) such that
    Sum_j (n_j/N)*x_j(eta,T) == x_target -- the self-consistency condition
    Ferreira et al. 2024's Eq. 14-15 express as an order-N polynomial in
    eta. eta is a Lagrange-multiplier "activity" of species B; the mass-
    action law (cluster_fractions) is monotonic in eta (higher eta favors
    higher-n_j clusters), so the root is unique and brentq is safe here.

    x_target in {0, 1} is special-cased (returns 0.0 / math.inf directly)
    since the true root sits at eta's domain boundary, where brentq would
    either fail to bracket or blow up numerically evaluating eta**N there.
    """
    if x_target <= 0.0:
        return 0.0
    if x_target >= 1.0:
        return math.inf

    def f(log_eta):
        return _mean_composition(math.exp(log_eta), T, g, delta_e, N) - x_target

    # Bracket in log(eta) space (eta ranges over many orders of magnitude
    # for x_target near 0 or 1) -- expand geometrically until f changes
    # sign, same "don't assume a fixed bracket" caution as any other
    # root-find over an unknown-scale variable in this suite.
    lo, hi = -50.0, 50.0
    while f(lo) > 0:
        lo -= 50.0
    while f(hi) < 0:
        hi += 50.0
    log_eta = brentq(f, lo, hi, xtol=1e-12, rtol=1e-12)
    return math.exp(log_eta)


def mixing_enthalpy(x_j: np.ndarray, delta_e: np.ndarray) -> float:
    """Delta_H_mix = Sum_j x_j * Delta_j (Ferreira et al. 2024 Eq. 3, up to
    an overall normalization constant M their paper applies and we don't
    -- see NORMALIZATION note below). Units: eV per alloy sublattice SITE.

    NORMALIZATION: this is per sublattice site, not per N-site cluster and
    not per formula unit of the input structure. Every report line in
    stb-gqcaAnalysis states this explicitly -- getting the normalization
    wrong doesn't fail loudly, it just silently scales every result by a
    constant factor.
    """
    return float(np.sum(np.asarray(x_j, dtype=float) * np.asarray(delta_e, dtype=float)))


def mixing_entropy(x_j: np.ndarray, g: list[int]) -> float:
    """Delta_S_mix = -kB * Sum_j x_j * ln(x_j / g_j), eV/K per alloy
    sublattice site (same normalization as mixing_enthalpy).

    [UNVERIFIED] This is the standard "gas of non-interacting clusters"
    Boltzmann/Stirling entropy (same spirit as Rao et al., Acta Mater.
    2022, and the reference open-source implementation WMD-group/
    GQCA_alloys). Ferreira et al. 2024's own restatement of the entropy
    ("Delta_S = kB ln W", Eq. 5-6) could not be pinned to a fully closed
    form during this tool's development -- it may include an additional
    term beyond the cluster-mixing entropy used here (their EGQCA extends
    the base GQCA formalism, and their W may combine both a species-
    distribution term and this cluster-configuration term). Flagged here
    exactly as HER/OER flag their own unverified literature constants
    (see e.g. oer_analysis.py's own [UNVERIFIED] docstring notes) --
    documented, not hidden, pending a primary-source (Sher et al. 1987)
    confirmation.

    This specific functional form is also WHY Delta_G_mix(x) always
    comes out convex, structurally unable to predict a miscibility gap
    -- see common_tangent_gaps' own [STRUCTURAL LIMITATION] docstring
    note. Any future correction to this entropy formula that adds
    cluster-cluster correlation terms (rather than staying a pure
    independent-species Boltzmann/Stirling form) is the most likely path
    to eventually letting this module capture genuine phase separation.
    """
    x_j = np.asarray(x_j, dtype=float)
    g = np.asarray(g, dtype=float)
    total = 0.0
    for xj, gj in zip(x_j, g):
        if xj > 0:
            total += xj * math.log(xj / gj)
    return -_BOLTZMANN_EV_K * total


def free_energy(H: float, S: float, T: float) -> float:
    """Delta_G_mix = Delta_H_mix - T*Delta_S_mix."""
    return H - T * S


def random_reference_fractions(x: float, N: int) -> np.ndarray:
    """Binomial(N, x) distribution -- the SRO=0/fully-random reference
    x_j = C(N,j) * x^j * (1-x)^(N-j). Also the closed-form solution of the
    mass-action law when every cluster energy e_j is equal (Delta_j == 0
    for all j): the weights collapse to g_j*eta^n_j alone, and matching
    the composition constraint gives eta = x/(1-x) -- substituting back
    gives exactly this binomial distribution. Used both as the SRO
    reference in sro_parameter and as the analytic gate in this module's
    own test suite (see test/4-workflow/15-gqca/analysis/test.sh).
    """
    g = degeneracies(N)
    return np.array([g[j] * x**j * (1.0 - x) ** (N - j) for j in range(N + 1)], dtype=float)


def _cross(o, a, b) -> float:
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def _lower_hull(points: list[tuple[float, float]]) -> list[tuple[float, float]]:
    """Andrew's monotone-chain lower convex hull of `points` (already
    sorted by x ascending) -- standard O(n) pass once sorted, popping any
    vertex that would make a non-left ("cross <= 0") turn. Verified live
    during this function's development against a textbook regular-
    solution model (Delta_G_mix = Omega*x*(1-x) + RT*(x*ln(x) +
    (1-x)*ln(1-x)), Omega/RT=3, above the 2.0 critical ratio) --
    correctly recovers a symmetric gap at x~0.076/0.924, matching the
    known analytic behavior of that model.
    """
    hull: list[tuple[float, float]] = []
    for p in points:
        while len(hull) >= 2 and _cross(hull[-2], hull[-1], p) <= 0:
            hull.pop()
        hull.append(p)
    return hull


def common_tangent_gaps(x_grid: np.ndarray, g_grid: np.ndarray,
                         tol: float = 1e-9) -> list[tuple[float, float]]:
    """Miscibility gap(s) at one fixed T, found via the common-tangent
    (lower convex hull) construction on a DISCRETE Delta_G_mix(x) grid.

    Physics: at fixed T, the equilibrium (lowest free energy) state at an
    overall composition x inside a miscibility gap is NOT a single
    homogeneous alloy but a phase-separated mixture of two alloys at the
    gap's boundary compositions x_alpha/x_beta (the classic common-
    tangent / lever-rule construction) -- the boundary compositions are
    exactly the two points where a single straight line is tangent to
    Delta_G_mix(x) at both of them simultaneously, i.e. the endpoints of
    whichever lower-hull edge "skips over" one or more of the original
    grid points (those skipped points lie ABOVE the hull -- their
    homogeneous free energy is higher than the phase-separated mixture's).

    [STRUCTURAL LIMITATION -- verified, not just suspected] Within the
    independent-cluster (Sher-GQCA) formalism this whole module
    implements, Delta_G_mix(x) is PROVABLY CONVEX at every composition,
    for ANY choice of Delta_j and ANY N -- so in practice this function
    will essentially never find a gap on real GQCA output, no matter how
    unfavorable the cluster energies are. This isn't a numerical
    coincidence: Delta_G_mix(x) is obtained by minimizing a linear
    energy term minus T times a concave entropy over the cluster
    probabilities x_j, subject to a LINEAR constraint (mean composition
    = x) -- a textbook convex-analysis result (partial minimization of a
    jointly convex function over a linear constraint is itself convex in
    the remaining variable) guarantees the result is convex in x
    regardless of Delta_j/g_j/N. Verified numerically here across dozens
    of Delta_j/T/N combinations (including N=4, deliberately unfavorable
    AB-type energies up to 10 eV, and asymmetric Delta_0 != Delta_2) --
    not one produced a non-convex curve. Physically: this "gas of
    independent, non-interacting clusters" treatment captures SHORT-
    RANGE order (x_j deviating from the random/binomial distribution)
    but not genuine phase separation, which needs correlations BETWEEN
    clusters (the same reason a chain of decoupled 1D Ising spins never
    orders while a coupled 2D/3D lattice does) -- that requires a
    different, harder theory (a cluster variation method with self-
    consistent cluster-cluster coupling, or a full cluster-expansion +
    Monte Carlo treatment), out of scope here. A "no gap found" result
    is therefore an INTRINSIC property of this model, not a physical
    statement about the real alloy -- stb-gqcaAnalysis's own report
    states this explicitly every time this function runs, so a silent
    "no gap" is never misread as "this alloy is fully miscible."

    Returns a list of (x_alpha, x_beta) tuples, one per contiguous gap
    found (usually 0 or 1 for a simple binary system, but a more complex
    Delta_j set could in principle produce more than one). Always
    includes the (0.0, Delta_G_mix=0) and (1.0, Delta_G_mix=0) endpoints
    in the hull construction even if not present in `x_grid` -- Delta_G_
    mix is 0 by definition at the pure end-members regardless of what the
    grid's own min/max happen to be.

    [RESOLUTION CAVEAT] This is a discrete-grid estimate, not a
    continuous-curve one: x_alpha/x_beta can only land exactly ON a
    --x-points grid point, so the true (continuous) binodal compositions
    are only bracketed to within one grid spacing. Refine by increasing
    --x-points for a tighter estimate, not by fitting/extrapolating
    (Delta_G_mix's shape near the gap edges isn't guaranteed smooth
    enough for that, especially so close to the [UNVERIFIED] entropy
    term -- see mixing_entropy's own docstring).
    """
    x_grid = np.asarray(x_grid, dtype=float)
    g_grid = np.asarray(g_grid, dtype=float)
    points = sorted(zip(x_grid.tolist(), g_grid.tolist()))
    if points[0][0] > tol:
        points = [(0.0, 0.0)] + points
    if points[-1][0] < 1.0 - tol:
        points = points + [(1.0, 0.0)]

    hull = _lower_hull(points)
    xs_all = [p[0] for p in points]
    gaps = []
    for i in range(len(hull) - 1):
        xa, xb = hull[i][0], hull[i + 1][0]
        ia = xs_all.index(xa)
        ib = xs_all.index(xb)
        if ib - ia > 1:
            gaps.append((xa, xb))
    return gaps


def sro_parameter(x_j: np.ndarray, x: float, N: int) -> float:
    """[UNVERIFIED] Provisional deviation-from-random measure, pending a
    primary-source (Sher et al. 1987) confirmation of the exact Warren-
    Cowley-style normalization for a cluster-probability-based SRO
    parameter (Cowley's own original 1950 definition, Phys. Rev. 77, 669,
    is expressed in terms of pair-correlation shell probabilities, which
    the N=2 pair cluster's x_1 (mixed AB fraction) directly corresponds
    to -- but the precise normalization convention, e.g. whether it's
    bounded to [-1,1] the way Cowley's alpha_i parameters are, has not
    been independently confirmed here).

    Provisional definition: fractional deviation of the mixed-cluster
    (AB, j=1..N-1 lumped) probability from its random-mixing value,
    normalized to be dimensionless and (for the pair case) analogous in
    sign convention to Cowley's alpha (negative = ordering/AB-preferring,
    positive = clustering/like-neighbor-preferring).
    """
    x_j = np.asarray(x_j, dtype=float)
    x_random = random_reference_fractions(x, N)
    # "Mixed" clusters are every j strictly between the two pure end
    # members (0 < j < N) -- for N=2 this is just j=1 (the AB bond).
    p_mixed = float(np.sum(x_j[1:N]))
    p_mixed_random = float(np.sum(x_random[1:N]))
    if p_mixed_random == 0.0:
        return 0.0
    return (p_mixed_random - p_mixed) / p_mixed_random
