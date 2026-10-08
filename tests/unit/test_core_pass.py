import numpy as np
import pytest

from search.config import RunConfig
from search.core.corepass import run_core
from search.errors import SeedError
from search.presets import BUFFA2010, SEARCH2014
from tests.synthetic_data import make_planted


def test_core_result_shapes_and_tables() -> None:
    data = make_planted(n_genes=40, n_samples=30, module_size=10, n_seeds=3, seed=4)
    res = run_core(data.expr, data.seeds, BUFFA2010)
    assert res.n_genes == 40 and res.n_seeds == 3 and res.n_samples == 30
    assert res.stat.shape == (40, 3) and res.delta.shape == (40, 3)
    assert res.C_rank is not None and 0 <= np.nanmin(res.C_rank) <= np.nanmax(res.C_rank) <= 1
    assert res.is_seed.sum() == 3 and list(np.asarray(res.genes)[res.seed_idx]) == data.seeds
    table = res.gene_table()
    assert list(table.columns) == [
        "gene",
        "C_raw",
        "C_rank",
        "n_edges",
        "n_positive_edges",
        "n_negative_edges",
        "positive_consistent",
        "weighted_seed_rho",
        "is_seed",
    ]
    hoods = res.neighbourhoods()
    for seed in data.seeds:
        assert seed in hoods[seed].index  # self edge
    pairs = res.seed_pairs()
    assert len(pairs) == 9 and pairs["O"].between(0, 1).all()
    assert res.S is not None and res.F is not None and res.config_hash == BUFFA2010.core_hash()


def test_search2014_core_has_no_rank_and_binary_gamma() -> None:
    data = make_planted(n_genes=40, n_samples=30, module_size=10, n_seeds=3, seed=5)
    res = run_core(data.expr, data.seeds, SEARCH2014)
    assert res.C_rank is None
    assert set(np.unique(res.gamma)) <= {0.0, 1.0}
    assert res.C_raw.max() <= 1.0


def test_seed_errors_and_universe() -> None:
    data = make_planted(n_genes=30, n_samples=20, module_size=8, n_seeds=3, seed=6)
    with pytest.raises(SeedError):
        run_core(data.expr, data.seeds[:1], BUFFA2010)
    with pytest.raises(SeedError):
        run_core(data.expr, data.seeds + ["nope"], BUFFA2010)
    universe = list(data.expr.genes[:20])
    res = run_core(data.expr, data.seeds, BUFFA2010, universe=universe)
    assert res.n_genes == 20
    with pytest.raises(SeedError):
        run_core(data.expr, data.seeds, BUFFA2010, universe=list(data.expr.genes[10:]))


def test_explicit_sample_indices_and_column_order_invariance() -> None:
    data = make_planted(n_genes=30, n_samples=20, module_size=8, n_seeds=3, seed=7)
    perm = np.random.default_rng(0).permutation(20)
    full = run_core(data.expr, data.seeds, BUFFA2010)
    permuted = run_core(data.expr, data.seeds, BUFFA2010, sample_idx=perm)
    np.testing.assert_allclose(full.C_raw, permuted.C_raw)
    boot = np.random.default_rng(1).integers(0, 20, 20)
    direct = run_core(data.expr.take_samples(boot), data.seeds, BUFFA2010)
    via_idx = run_core(data.expr, data.seeds, BUFFA2010, sample_idx=boot)
    np.testing.assert_allclose(direct.C_raw, via_idx.C_raw)


def test_shared_neighbourhood_can_be_switched_off() -> None:
    data = make_planted(n_genes=30, n_samples=20, module_size=8, n_seeds=3, seed=8)
    cfg = RunConfig().with_overrides({"overlap.shared_neighbourhood": "never"})
    res = run_core(data.expr, data.seeds, cfg)
    assert res.S is None and res.F is None and res.O.shape == (3, 3)
