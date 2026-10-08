import time

import numpy as np
import pandas as pd
import pytest

from search.core.corepass import run_core
from search.data import ExpressionSet
from search.presets import BUFFA2010


@pytest.mark.bench
def test_core_pass_under_one_second_on_cohort_sized_matrix() -> None:
    rng = np.random.default_rng(0)
    n_genes, n_samples = 16548, 95
    values = rng.normal(size=(n_genes, n_samples))
    expr = ExpressionSet(
        values, pd.Index([f"g{i}" for i in range(n_genes)]), pd.Index(range(n_samples))
    )
    seeds = [f"g{i}" for i in range(10)]
    run_core(expr, seeds, BUFFA2010)  # warm-up
    start = time.perf_counter()
    for _ in range(3):
        run_core(expr, seeds, BUFFA2010)
    elapsed = (time.perf_counter() - start) / 3
    print(f"\ncore pass on {n_genes}x{n_samples}, 10 seeds: {elapsed:.3f} s")
    assert elapsed < 1.0
