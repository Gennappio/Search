"""Per-sample signature score. The score of Buffa 2010 is the median of the absolute expression
of the signature genes, ranked between 0 and 1 across samples. Other scorers come from
pysigscore when it is installed."""

from __future__ import annotations

import tempfile
from pathlib import Path

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


def pysigscore_scores(
    expr: ExpressionSet, genes: list[str], scorer: str = "Median", workdir: Path | None = None
) -> pd.Series:
    """One score per sample from pysigscore's direct scorers (Median, Mean, Z, ...)."""
    try:
        from pysigscore import Scorer
    except ImportError as exc:  # pragma: no cover - depends on the optional package
        raise ImportError(
            "pysigscore is not installed; it is not on PyPI, install it into the environment "
            "by hand: uv pip install --python .venv/bin/python <path or git url>"
        ) from exc
    present = [g for g in dict.fromkeys(genes) if g in expr.genes]
    frame = pd.DataFrame(
        np.abs(expr.values).T,
        index=[str(s) for s in expr.samples],
        columns=[str(g) for g in expr.genes],
    )
    base = (
        Path(workdir)
        if workdir is not None
        else Path(tempfile.mkdtemp(prefix="search_pysigscore_"))
    )
    engine = Scorer(name="search", path=str(base))
    engine.import_adata_as_dataframe(frame)
    engine.import_geneset_dict({"signature": present})
    results = engine.compute_scores(
        transform_list_data=["raw"],
        transform_list_enrich=["raw"],
        enrich_grid={scorer: {}},
        p_value_grid={scorer: 0},
        return_results=True,
    )
    key = next((k for k in results if k.startswith(f"{scorer}_score")), None)
    if key is None:
        raise DataError(f"pysigscore returned no '{scorer}_score' table; keys: {sorted(results)}")
    table = results[key]
    column = "signature" if "signature" in table.columns else table.columns[0]
    scores = table[column]
    if not scores.index.astype(str).isin(frame.index).all():
        scores = pd.Series(scores.to_numpy(), index=frame.index[: len(scores)])
    return scores.reindex(frame.index).astype(float)


def score_samples(
    expr: ExpressionSet, genes: list[str], cfg: ScoreConfig, dataset: str = ""
) -> pd.DataFrame:
    if cfg.method == "median_abs_rank01":
        return median_abs_rank01(expr, genes, cfg.min_genes, dataset)
    scorer = cfg.method.removeprefix("pysigscore_").capitalize()
    scores = pysigscore_scores(expr, genes, scorer)
    present = [g for g in dict.fromkeys(genes) if g in expr.genes]
    return pd.DataFrame(
        {
            "dataset": dataset,
            "sample": np.asarray(expr.samples),
            f"E_{scorer.lower()}": scores.to_numpy(),
            "HS_rank01": rank01(scores.to_numpy(dtype=float)),
            "n_signature_genes_present": len(present),
            "n_signature_genes_total": len(set(genes)),
        }
    )
