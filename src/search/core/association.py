"""Association between every seed and every gene: the first step of the core pass.

Correlation methods are vectorised: rows are ranked (Spearman) or not (Pearson), centred and
scaled to unit norm, and the genes x seeds matrix is one matrix product. Missing values are not
imputed: a gene with a missing value gets no association (NaN), and the QC stage warns about it.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from search.config import AssociationConfig
from search.stats import rank_rows, t_approx_pvalue

FloatArray = NDArray[np.float64]


@dataclass
class AssociationResult:
    stat: FloatArray  # genes x seeds: correlation or regression coefficient
    p: FloatArray | None  # genes x seeds: two-sided p-value, if available
    method: str
    n_samples: int
    coef_table: pd.DataFrame | None = None  # full coefficient tables (ols/glm)


def standardise_rows(x: FloatArray) -> FloatArray:
    """Centre every row and scale it to unit norm. Rows with NaN or zero variance become NaN."""
    x = np.asarray(x, dtype=np.float64)
    centred = x - x.mean(axis=1, keepdims=True)
    norm = np.sqrt((centred * centred).sum(axis=1, keepdims=True))
    with np.errstate(invalid="ignore", divide="ignore"):
        z = centred / norm
    bad = ~np.isfinite(z).all(axis=1)
    if bad.any():
        z[bad] = np.nan
    return z


def pearson_block(values: FloatArray, seed_idx: NDArray[np.intp]) -> FloatArray:
    """Pearson correlation of every row with the seed rows: genes x seeds."""
    z = standardise_rows(values)
    rho = z @ z[seed_idx].T
    return np.clip(rho, -1.0, 1.0)


def spearman_block(values: FloatArray, seed_idx: NDArray[np.intp]) -> FloatArray:
    """Spearman correlation of every row with the seed rows: genes x seeds."""
    return pearson_block(rank_rows(values), seed_idx)


def associate(
    values: FloatArray,
    seed_idx: NDArray[np.intp],
    cfg: AssociationConfig,
    sample_table: pd.DataFrame | None = None,
    seed_names: list[str] | None = None,
) -> AssociationResult:
    """Dispatch on the configured method."""
    n = int(values.shape[1])
    if cfg.method in ("spearman", "pearson"):
        rho = (
            spearman_block(values, seed_idx)
            if cfg.method == "spearman"
            else pearson_block(values, seed_idx)
        )
        return AssociationResult(
            stat=rho, p=t_approx_pvalue(rho, n), method=cfg.method, n_samples=n
        )
    if cfg.method in ("ols", "glm"):
        from search.core.linear import associate_linear

        return associate_linear(values, seed_idx, cfg, sample_table, seed_names)
    raise ValueError(f"unknown association method {cfg.method}")
