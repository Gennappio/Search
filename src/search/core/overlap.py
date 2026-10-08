"""Seed-to-seed overlap: binary Jaccard O (SEARCH 2014), weighted S and sign concordance F
(Buffa 2010, equations 3 and 4)."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]
BoolArray = NDArray[np.bool_]


def jaccard_binary(delta: BoolArray) -> FloatArray:
    """O[a, b] = |cluster_a and cluster_b| / |cluster_a or cluster_b|, all genes counted."""
    d = delta.astype(np.float64)
    inter = d.T @ d
    sizes = d.sum(axis=0)
    union = sizes[:, None] + sizes[None, :] - inter
    with np.errstate(invalid="ignore", divide="ignore"):
        o = np.where(union > 0, inter / union, np.nan)
    return np.asarray(o, dtype=np.float64)


def omax(o: FloatArray) -> FloatArray:
    """Largest overlap of every seed with any other seed."""
    n = o.shape[0]
    out = np.full(n, np.nan)
    if n < 2:
        return out
    masked = o.copy()
    np.fill_diagonal(masked, np.nan)
    has_value = np.isfinite(masked).any(axis=1)
    out[has_value] = np.nanmax(masked[has_value], axis=1)
    return out


def shared_neighbourhood(
    gamma: FloatArray, stat: FloatArray, seed_idx: NDArray[np.intp]
) -> tuple[FloatArray, FloatArray]:
    """Weighted Jaccard S and sign concordance F between seeds, excluding the two seeds' rows.

    S[a, b] = sum_g min(gamma_ga, gamma_gb) / sum_g max(gamma_ga, gamma_gb)
    F[a, b] = sum_g min(gamma_ga, gamma_gb) sign(stat_ga stat_gb) / sum_g min(gamma_ga, gamma_gb)
    """
    n_seeds = gamma.shape[1]
    s = np.full((n_seeds, n_seeds), np.nan)
    f = np.full((n_seeds, n_seeds), np.nan)
    np.fill_diagonal(s, 1.0)
    np.fill_diagonal(f, 1.0)
    with np.errstate(invalid="ignore"):
        sign = np.sign(stat)
    sign = np.nan_to_num(sign, nan=0.0)
    for a in range(n_seeds - 1):
        ga = gamma[:, a]
        sa = sign[:, a]
        ra = seed_idx[a]
        for b in range(a + 1, n_seeds):
            gb = gamma[:, b]
            rb = seed_idx[b]
            lo = np.minimum(ga, gb)
            hi = np.maximum(ga, gb)
            prod_sign = sa * sign[:, b]
            sum_lo = lo.sum() - lo[ra] - lo[rb]
            sum_hi = hi.sum() - hi[ra] - hi[rb]
            weighted = lo * prod_sign
            sum_signed = weighted.sum() - weighted[ra] - weighted[rb]
            s_ab = sum_lo / sum_hi if sum_hi > 0 else np.nan
            f_ab = sum_signed / sum_lo if sum_lo > 0 else np.nan
            s[a, b] = s[b, a] = s_ab
            f[a, b] = f[b, a] = f_ab
    return s, f
