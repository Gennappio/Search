"""Per-sample signature score: the median of the absolute expression of the signature genes,
ranked between 0 and 1 across the samples of a cohort (Buffa 2010)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from search.config import ScoreConfig
from search.data import ExpressionSet
from search.errors import DataError
from search.stats import rank01


def median_abs_expression(
    expr: ExpressionSet, genes: list[str], min_genes: int = 3
) -> tuple[np.ndarray, list[str]]:
    present = [g for g in dict.fromkeys(genes) if g in expr.genes]
    if len(present) < min_genes:
        raise DataError(
            f"only {len(present)} of {len(set(genes))} signature genes present, "
            f"fewer than {min_genes}"
        )
    rows = expr.gene_positions(present)
    values = np.nanmedian(np.abs(expr.values[rows]), axis=0)
    return values, present


def median_abs_rank01(
    expr: ExpressionSet, genes: list[str], min_genes: int = 3, dataset: str = ""
) -> pd.DataFrame:
    e, present = median_abs_expression(expr, genes, min_genes)
    return pd.DataFrame(
        {
            "dataset": dataset,
            "sample": np.asarray(expr.samples),
            "E_median_abs_expression": e,
            "HS_rank01": rank01(e),
            "n_signature_genes_present": len(present),
            "n_signature_genes_total": len(set(genes)),
        }
    )


def score_samples(
    expr: ExpressionSet, genes: list[str], cfg: ScoreConfig, dataset: str = ""
) -> pd.DataFrame:
    if cfg.method != "median_abs_rank01":  # pragma: no cover - the config rejects other values
        raise DataError(f"unknown score method {cfg.method!r}")
    return median_abs_rank01(expr, genes, cfg.min_genes, dataset)
