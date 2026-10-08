"""Under pure noise the empirical p-values of the null model must be uniform."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import stats as sps

from search.core.corepass import run_core
from search.data import ExpressionSet
from search.null import run_internal_null
from search.presets import BUFFA2010
from search.resampling import derive_rng, random_seed_set


@pytest.mark.slow
@pytest.mark.parametrize("pooling", ["pooled", "per_gene"])
def test_null_pvalues_are_uniform_on_noise(tmp_path: Path, pooling: str) -> None:
    rng = np.random.default_rng(0)
    n_genes, n_samples = 120, 30
    values = rng.normal(size=(n_genes, n_samples))
    genes = pd.Index([f"g{i}" for i in range(n_genes)])
    expr = ExpressionSet(values, genes, pd.Index([f"s{j}" for j in range(n_samples)]))
    # a lenient edge rule so that the statistic is not degenerate at zero
    cfg = BUFFA2010.with_overrides(
        {
            "edge.rule": "abs_threshold",
            "edge.threshold": 0.3,
            "null.n_rep": 200,
            "null.chunk_size": 50,
            "null.pooling": pooling,
        }
    )
    all_p = []
    for trial in range(1, 16):
        seeds = random_seed_set(np.asarray(genes), 4, derive_rng(500, trial))
        observed = run_core(expr, seeds, cfg).C_raw
        result = run_internal_null(
            expr,
            seeds,
            cfg,
            universe=None,
            observed=observed,
            ckpt_dir=tmp_path / f"{pooling}_{trial}",
            input_hash="x",
        )
        p = result.p
        mask = np.ones(n_genes, dtype=bool)
        mask[expr.gene_positions(seeds)] = False
        all_p.append(p[mask])
    pvals = np.concatenate(all_p)
    pvals = pvals[np.isfinite(pvals)]
    # the statistic is discrete and most genes have no edge, so p-values are conservative:
    # what must never happen is an excess of small p-values (anti-conservative behaviour)
    assert np.mean(pvals < 0.05) <= 0.07
    assert np.mean(pvals < 0.5) <= 0.55
    assert sps.kstest(pvals, "uniform", alternative="greater").pvalue > 0.01
