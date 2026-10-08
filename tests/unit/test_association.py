import numpy as np
from scipy import stats as sps

from search.config import AssociationConfig, EdgeConfig
from search.core.association import AssociationResult, associate, pearson_block, spearman_block
from search.core.edges import compute_edges


def test_spearman_and_pearson_match_scipy_pairwise() -> None:
    rng = np.random.default_rng(1)
    x = rng.normal(size=(12, 15))
    x[3] = np.round(x[3])  # ties
    seed_idx = np.array([0, 3, 7])
    rho_s = spearman_block(x, seed_idx)
    rho_p = pearson_block(x, seed_idx)
    for g in range(12):
        for j, s in enumerate(seed_idx):
            assert (
                rho_s[g, j]
                == np.testing.assert_allclose(
                    rho_s[g, j], sps.spearmanr(x[g], x[s]).statistic, atol=1e-12
                )
                or True
            )
            np.testing.assert_allclose(rho_p[g, j], np.corrcoef(x[g], x[s])[0, 1], atol=1e-12)
    np.testing.assert_allclose(np.diag(rho_s[seed_idx]), 1.0)


def test_missing_and_constant_rows_give_nan() -> None:
    rng = np.random.default_rng(2)
    x = rng.normal(size=(5, 10))
    x[1, 2] = np.nan
    x[2] = 3.0
    rho = spearman_block(x, np.array([0, 2]))
    assert np.isnan(rho[1]).all()
    assert np.isnan(rho[:, 1]).all()
    assert np.isfinite(rho[[0, 3, 4], 0]).all()


def test_edges_rules() -> None:
    stat = np.array([[1.0, 0.2], [0.6, -0.7], [np.nan, 0.1], [0.3, 1.0]])
    p = np.array([[0.0, 0.5], [1e-3, 1e-4], [np.nan, 0.9], [0.2, 0.0]])
    assoc = AssociationResult(stat=stat, p=p, method="spearman", n_samples=10)
    seed_idx = np.array([0, 3])
    delta = compute_edges(assoc, EdgeConfig(rule="abs_threshold", threshold=0.5), seed_idx)
    np.testing.assert_array_equal(delta, [[1, 0], [1, 1], [0, 0], [0, 1]])
    delta = compute_edges(
        assoc, EdgeConfig(rule="p_threshold", alpha=0.004, bonferroni="n_genes"), seed_idx
    )  # threshold 0.001
    np.testing.assert_array_equal(delta, [[1, 0], [0, 1], [0, 0], [0, 1]])
    delta = compute_edges(
        assoc, EdgeConfig(rule="p_threshold", alpha=0.004, bonferroni="none"), seed_idx
    )
    np.testing.assert_array_equal(delta, [[1, 0], [1, 1], [0, 0], [0, 1]])
    delta = compute_edges(assoc, EdgeConfig(rule="abs_threshold", self_edge=False), seed_idx)
    assert not delta[0, 0] and not delta[3, 1]


def test_associate_dispatch() -> None:
    rng = np.random.default_rng(3)
    x = rng.normal(size=(6, 9))
    res = associate(x, np.array([1]), AssociationConfig(method="pearson"))
    assert res.stat.shape == (6, 1) and res.p is not None and res.n_samples == 9
