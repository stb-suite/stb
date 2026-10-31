"""Model-independent atomic-environment fingerprints and a set of them, used
to tell how new a configuration is relative to others: stb-mlffActiveLearning's
novelty/diversity filter and stb-mlff's near-duplicate check.
"""

import math

import numpy as np

# Atomic-environment fingerprint used for the novelty/diversity filter: the
# sorted distances from an atom to its FP_NEIGHBOURS nearest neighbours
# (periodic images included) within FP_CUTOFF Ang, padded with the cutoff.
# Model-independent (MACECalculator.get_descriptors does not work on the
# compiled model stb-mlffAnalysis also writes) and comparable across cell
# sizes, which a supercell candidate against an 8-atom training set needs.
FP_NEIGHBOURS = 16
FP_CUTOFF = 6.0


def environment_fingerprints(atoms, n_neighbours=FP_NEIGHBOURS, cutoff=FP_CUTOFF):
    """(n_atoms, n_neighbours) array: sorted distances (Ang) from each atom to
    its nearest neighbours within `cutoff`, padded with `cutoff`."""
    from ase.neighborlist import neighbor_list
    i_idx, dists = neighbor_list("id", atoms, cutoff)
    fp = np.full((len(atoms), n_neighbours), cutoff)
    order = np.lexsort((dists, i_idx))
    i_idx, dists = i_idx[order], dists[order]
    starts = np.searchsorted(i_idx, np.arange(len(atoms)))
    ends = np.searchsorted(i_idx, np.arange(len(atoms)), side="right")
    for a, (s, e) in enumerate(zip(starts, ends)):
        d = dists[s:min(e, s + n_neighbours)]
        fp[a, :len(d)] = d
    return fp


class EnvironmentSet:
    """Atomic-environment fingerprints of a set of configurations, kept per
    centre species (an environment is only compared to environments of the
    same element) and tagged with the configuration they came from."""

    def __init__(self):
        self.fps = {}      # symbol -> list of (n, d) arrays
        self.cfg = {}      # symbol -> list of config-id arrays
        self.n_configs = 0
        self._trees = {}

    def add(self, atoms):
        fp = environment_fingerprints(atoms)
        symbols = np.array(atoms.get_chemical_symbols())
        for sym in set(symbols):
            mask = symbols == sym
            self.fps.setdefault(sym, []).append(fp[mask])
            self.cfg.setdefault(sym, []).append(np.full(mask.sum(), self.n_configs))
            self._trees.pop(sym, None)
        self.n_configs += 1

    def _tree(self, sym):
        from scipy.spatial import cKDTree
        if sym not in self._trees:
            self._trees[sym] = (cKDTree(np.vstack(self.fps[sym])), np.concatenate(self.cfg[sym]))
        return self._trees[sym]

    def typical_spacing(self, k=64):
        """Median distance from an environment to the nearest one belonging
        to a DIFFERENT configuration -- the natural unit of 'how far apart
        two independent samples of this set are'. None with < 2 configs."""
        if self.n_configs < 2:
            return None
        nearest = []
        for sym in self.fps:
            tree, cfg = self._tree(sym)
            data = tree.data
            kk = min(k, len(data))
            dist, idx = tree.query(data, k=kk)
            dist, idx = np.atleast_2d(dist), np.atleast_2d(idx)
            for row in range(len(data)):
                other = cfg[idx[row]] != cfg[row]
                if other.any():
                    nearest.append(dist[row][other][0])
        positive = [d for d in nearest if d > 0]
        return float(np.median(positive)) if positive else None

    def novelty(self, atoms, scale):
        """Max over the atoms of `atoms` of the distance from its environment
        to the nearest one in this set (same element), in units of `scale`.
        inf when an element of `atoms` is absent from the set."""
        if self.n_configs == 0:
            return math.inf
        fp = environment_fingerprints(atoms)
        symbols = np.array(atoms.get_chemical_symbols())
        worst = 0.0
        for sym in set(symbols):
            if sym not in self.fps:
                return math.inf
            tree, _ = self._tree(sym)
            dist, _ = tree.query(fp[symbols == sym], k=1)
            worst = max(worst, float(np.max(dist)))
        return worst / scale
