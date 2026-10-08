"""Properties that must hold for any input."""

import numpy as np
import pandas as pd
from hypothesis import given, settings
from hypothesis import strategies as st

from search.core.corepass import run_core
from search.data import ExpressionSet
from search.presets import BUFFA2010, SEARCH2014


def random_expr(seed: int, n_genes: int, n_samples: int) -> ExpressionSet:
    rng = np.random.default_rng(seed)
    values = rng.normal(size=(n_genes, n_samples))
    # a planted block so that some edges exist
    values[: n_genes // 2] += rng.normal(size=n_samples)[None, :]
    genes = [f"g{i}" for i in range(n_genes)]
    return ExpressionSet(values, pd.Index(genes), pd.Index([f"s{j}" for j in range(n_samples)]))


@settings(max_examples=20, deadline=None)
@given(
    seed=st.integers(0, 10**6),
    n_genes=st.integers(8, 25),
    n_samples=st.integers(6, 20),
    n_seeds=st.integers(2, 4),
    preset=st.sampled_from([BUFFA2010, SEARCH2014]),
)
def test_gene_and_sample_order_do_not_matter(seed, n_genes, n_samples, n_seeds, preset) -> None:
    expr = random_expr(seed, n_genes, n_samples)
    seeds = [f"g{i}" for i in range(n_seeds)]
    rng = np.random.default_rng(seed + 1)
    gene_perm = rng.permutation(n_genes)
    sample_perm = rng.permutation(n_samples)
    shuffled = ExpressionSet(
        expr.values[gene_perm][:, sample_perm], expr.genes[gene_perm], expr.samples[sample_perm]
    )
    a = run_core(expr, seeds, preset).gene_table().set_index("gene").sort_index()
    b = run_core(shuffled, seeds, preset).gene_table().set_index("gene").sort_index()
    for column in ["C_raw", "C_rank", "n_edges", "weighted_seed_rho"]:
        np.testing.assert_allclose(a[column].to_numpy(float), b[column].to_numpy(float), atol=1e-10)
    assert (a["positive_consistent"] == b["positive_consistent"]).all()


@settings(max_examples=20, deadline=None)
@given(
    seed=st.integers(0, 10**6),
    n_genes=st.integers(8, 25),
    n_samples=st.integers(6, 20),
    n_seeds=st.integers(2, 4),
    preset=st.sampled_from([BUFFA2010, SEARCH2014]),
)
def test_ranges_and_symmetries(seed, n_genes, n_samples, n_seeds, preset) -> None:
    expr = random_expr(seed, n_genes, n_samples)
    res = run_core(expr, [f"g{i}" for i in range(n_seeds)], preset)
    assert np.nanmin(res.C_raw) >= 0 and np.nanmax(res.C_raw) <= 1
    if res.C_rank is not None:
        assert np.nanmin(res.C_rank) >= 0 and np.nanmax(res.C_rank) <= 1
    assert np.allclose(res.O, res.O.T, equal_nan=True) and np.allclose(np.diag(res.O), 1)
    assert res.S is not None and res.F is not None
    assert np.allclose(res.S, res.S.T, equal_nan=True) and np.allclose(np.diag(res.S), 1)
    finite_f = res.F[np.isfinite(res.F)]
    assert (finite_f >= -1 - 1e-12).all() and (finite_f <= 1 + 1e-12).all()
    assert (res.delta[res.seed_idx, np.arange(n_seeds)]).all()  # self edges
    assert (res.gamma >= 0).all() and (res.gamma[~res.delta] == 0).all()
    assert (res.n_edges >= res.n_positive_edges + res.n_negative_edges).all()
