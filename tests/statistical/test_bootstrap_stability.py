"""Two halves of the bootstrap replicates must tell the same story about the genes."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import stats as sps

from search.bootstrap import run_bootstrap
from search.checkpoint import rep_column
from search.presets import BUFFA2010
from tests.synthetic_data import make_planted


@pytest.mark.slow
def test_split_half_bootstrap_means_agree(tmp_path: Path) -> None:
    data = make_planted(n_genes=200, n_samples=60, module_size=40, n_seeds=5, strength=1.5, seed=9)
    cfg = BUFFA2010.with_overrides(
        {"resampling.n_rep": 60, "resampling.chunk_size": 30, "null.n_rep": 0}
    )
    summary = run_bootstrap(
        data.expr, data.seeds, cfg, universe=None, ckpt_dir=tmp_path, input_hash="x"
    )
    chunks = [
        pd.read_parquet(tmp_path / f"bootstrap_chunk_{i:04d}.parquet").set_index("gene")
        for i in (1, 2)
    ]
    first = chunks[0][[rep_column(r) for r in range(1, 31)]].mean(axis=1)
    second = chunks[1][[rep_column(r) for r in range(31, 61)]].mean(axis=1)
    # the means agree on the whole (Pearson) and on who is at the top; background genes with
    # tied, near-constant ranks make an all-gene rank correlation meaningless
    pearson = sps.pearsonr(first, second).statistic
    top_a = set(first.sort_values(ascending=False).head(30).index)
    top_b = set(second.sort_values(ascending=False).head(30).index)
    overlap = len(top_a & top_b) / 30
    module = [g for g in data.module if g not in data.seeds]
    rho_module = sps.spearmanr(first[module], second[module]).statistic
    assert pearson > 0.95, f"split-half Pearson {pearson:.3f}"
    assert overlap >= 0.8, f"top-30 overlap {overlap:.2f}"
    assert rho_module > 0.5, f"module rank correlation {rho_module:.3f}"
    np.testing.assert_allclose(summary.mean, pd.concat(chunks, axis=1).mean(axis=1).to_numpy())
