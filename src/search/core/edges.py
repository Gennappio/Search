"""From association statistics to edges of the seed network (the delta matrix)."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from search.config import EdgeConfig
from search.core.association import AssociationResult
from search.errors import ConfigError
from search.stats import p_adjust

BoolArray = NDArray[np.bool_]
DECIMALS = 12


def compute_edges(
    assoc: AssociationResult, cfg: EdgeConfig, seed_idx: NDArray[np.intp]
) -> BoolArray:
    """genes x seeds boolean matrix: True where the gene is in the cluster of the seed.

    ``abs_threshold`` keeps |stat| > threshold. ``p_threshold`` keeps p < alpha, with alpha
    divided by the number of genes when ``bonferroni`` is ``n_genes``; ``p_adjust`` can instead
    adjust the p-values of each seed before thresholding. Missing statistics never make an edge.
    The seed's own entry is forced to ``self_edge``.
    """
    stat = assoc.stat
    n_genes, n_seeds = stat.shape
    # decisions are taken on values rounded to 12 decimals, so that floating-point noise that
    # depends on the order of samples or genes cannot flip a statistic sitting on the threshold
    if cfg.rule == "abs_threshold":
        with np.errstate(invalid="ignore"):
            delta = np.round(np.abs(stat), DECIMALS) > cfg.threshold
    elif cfg.rule == "p_threshold":
        if assoc.p is None:
            raise ConfigError(
                f"edge rule p_threshold needs p-values; method {assoc.method} gives none"
            )
        p = assoc.p
        if cfg.p_adjust != "none":
            p = np.column_stack([p_adjust(p[:, j], cfg.p_adjust) for j in range(n_seeds)])
        threshold = cfg.alpha / n_genes if cfg.bonferroni == "n_genes" else cfg.alpha
        with np.errstate(invalid="ignore"):
            delta = np.round(p, DECIMALS) < threshold
    else:  # pragma: no cover - pydantic restricts the literal
        raise ConfigError(f"unknown edge rule {cfg.rule}")
    delta = np.asarray(delta, dtype=bool)
    delta[np.isnan(stat)] = False
    delta[seed_idx, np.arange(n_seeds)] = cfg.self_edge
    return delta
