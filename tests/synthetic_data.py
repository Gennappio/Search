"""Synthetic expression data with planted structure, shared by several test layers."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from search.data import ExpressionSet


@dataclass
class Planted:
    expr: ExpressionSet
    seeds: list[str]
    module: list[str]
    background: list[str]


def make_planted(
    n_genes: int = 300,
    n_samples: int = 60,
    module_size: int = 50,
    n_seeds: int = 5,
    strength: float = 1.0,
    noise: float = 1.0,
    seed: int = 0,
    second_module: int = 0,
) -> Planted:
    """A module of ``module_size`` genes driven by one hidden factor, the rest pure noise.

    The first ``n_seeds`` module genes are the seeds. ``second_module`` adds an unrelated
    module driven by a second factor (useful for pruning tests).
    """
    rng = np.random.default_rng(seed)
    factor = rng.normal(size=n_samples)
    values = rng.normal(scale=noise, size=(n_genes, n_samples))
    values[:module_size] += strength * factor[None, :]
    if second_module:
        other = rng.normal(size=n_samples)
        start = module_size
        values[start : start + second_module] += strength * other[None, :]
    genes = [f"g{i:04d}" for i in range(n_genes)]
    samples = [f"s{j:03d}" for j in range(n_samples)]
    expr = ExpressionSet(values, pd.Index(genes), pd.Index(samples))
    module = genes[:module_size]
    return Planted(
        expr=expr,
        seeds=module[:n_seeds],
        module=module,
        background=genes[module_size + second_module :],
    )


def auc(scores_pos: np.ndarray, scores_neg: np.ndarray) -> float:
    """Probability that a positive scores above a negative (Mann-Whitney)."""
    from scipy.stats import mannwhitneyu

    u = mannwhitneyu(scores_pos, scores_neg, alternative="two-sided").statistic
    return float(u / (len(scores_pos) * len(scores_neg)))
