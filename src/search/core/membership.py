"""Membership of every gene to every seed cluster (the gamma matrix)."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from search.config import MembershipConfig

FloatArray = NDArray[np.float64]
BoolArray = NDArray[np.bool_]


def compute_membership(
    delta: BoolArray, stat: FloatArray, cfg: MembershipConfig
) -> tuple[FloatArray, NDArray[np.int64]]:
    """Return (gamma, cluster_size).

    - ``binary``: gamma = delta (SEARCH 2014).
    - ``cluster_size``: delta divided by the number of genes in the seed's cluster (GBM script).
    - ``gene_normalised``: delta divided by the number of seeds the gene is connected to
      (Buffa 2010, equation 2 as printed).
    - ``sigmoid``: delta times a logistic function of |stat| (soft membership).
    """
    d = delta.astype(np.float64)
    cluster_size = delta.sum(axis=0).astype(np.int64)
    if cfg.mode == "binary":
        gamma = d
    elif cfg.mode == "cluster_size":
        denom = np.where(cluster_size > 0, cluster_size, 1).astype(np.float64)
        gamma = d / denom[None, :]
    elif cfg.mode == "gene_normalised":
        degree = delta.sum(axis=1)
        denom = np.where(degree > 0, degree, 1).astype(np.float64)
        gamma = d / denom[:, None]
    elif cfg.mode == "sigmoid":
        with np.errstate(over="ignore", invalid="ignore"):
            soft = 1.0 / (1.0 + np.exp(-cfg.sigmoid_k * (np.abs(stat) - cfg.sigmoid_x0)))
        gamma = d * np.nan_to_num(soft, nan=0.0)
    else:  # pragma: no cover
        raise ValueError(f"unknown membership mode {cfg.mode}")
    return gamma, cluster_size
