"""Connectivity of every gene to the seed set, and per-gene edge statistics."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from search.config import ConnectivityConfig

FloatArray = NDArray[np.float64]
BoolArray = NDArray[np.bool_]


def compute_connectivity(
    gamma: FloatArray, seed_idx: NDArray[np.intp], cfg: ConnectivityConfig
) -> FloatArray:
    """C_raw = sum_j w_j gamma_j / sum_j w_j, with w = 0 for a seed's own column when
    ``exclude_self``. With binary membership and ``fraction`` this is the D of SEARCH 2014:
    seeds are divided by K - 1 and other genes by K."""
    n_genes, n_seeds = gamma.shape
    w = np.ones((n_genes, n_seeds))
    if cfg.exclude_self:
        w[seed_idx, np.arange(n_seeds)] = 0.0
    numerator = (w * gamma).sum(axis=1)
    if cfg.normalise == "count":
        return numerator
    denominator = w.sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        c = np.where(denominator > 0, numerator / denominator, np.nan)
    return np.asarray(c, dtype=np.float64)


def edge_statistics(
    delta: BoolArray, stat: FloatArray, gamma: FloatArray
) -> dict[str, NDArray[np.generic]]:
    """Edge counts, sign consistency and membership-weighted association per gene."""
    with np.errstate(invalid="ignore"):
        positive = delta & (stat > 0)
        negative = delta & (stat < 0)
    n_edges = delta.sum(axis=1)
    n_pos = positive.sum(axis=1)
    n_neg = negative.sum(axis=1)
    positive_consistent = (n_edges > 0) & (n_neg == 0)
    gsum = gamma.sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        weighted = np.where(gsum > 0, np.nansum(gamma * stat, axis=1) / gsum, np.nan)
    return {
        "n_edges": n_edges.astype(np.int64),
        "n_positive_edges": n_pos.astype(np.int64),
        "n_negative_edges": n_neg.astype(np.int64),
        "positive_consistent": positive_consistent,
        "weighted_seed_rho": np.asarray(weighted, dtype=np.float64),
    }


def nodes_degree(delta: BoolArray, seed_idx: NDArray[np.intp]) -> FloatArray:
    """For every seed, the fraction of seeds (itself included) it is connected to."""
    n_seeds = delta.shape[1]
    return delta[seed_idx, :].sum(axis=1) / float(n_seeds)
