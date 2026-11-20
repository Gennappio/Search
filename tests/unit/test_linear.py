import numpy as np
import pandas as pd
import pytest
from scipy import stats as sps

from search.config import AssociationConfig, RunConfig
from search.core.association import associate
from search.core.corepass import run_core
from search.core.linear import design_matrix, ols_fwl
from search.data import ExpressionSet
from search.errors import ConfigError, DataError
from tests.synthetic_data import auc


def manual_ols(y: np.ndarray, x: np.ndarray, design: np.ndarray) -> tuple[float, float]:
    full = np.column_stack([design, x])
    coef, *_ = np.linalg.lstsq(full, y, rcond=None)
    resid = y - full @ coef
    df = len(y) - full.shape[1]
    sigma2 = resid @ resid / df
    cov = sigma2 * np.linalg.inv(full.T @ full)
    beta = coef[-1]
    t = beta / np.sqrt(cov[-1, -1])
    return float(beta), float(2 * sps.t.sf(abs(t), df))


def test_fwl_matches_least_squares_pair_by_pair() -> None:
    rng = np.random.default_rng(0)
    n_genes, n = 15, 40
    values = rng.normal(size=(n_genes, n))
    table = pd.DataFrame({"age": rng.normal(size=n), "batch": rng.choice(["a", "b", "c"], n)})
    design, names = design_matrix(table, ["age", "batch"], n)
    assert names == ["intercept", "age", "batch_b", "batch_c"]
    seed_idx = np.array([0, 4, 9])
    res = ols_fwl(values, seed_idx, design)
    for g in range(n_genes):
        for j, s in enumerate(seed_idx):
            if g == s:
                assert res.p[g, j] == 0.0 and res.stat[g, j] == pytest.approx(1.0)
                continue
            beta, p = manual_ols(values[s], values[g], design)
            assert res.stat[g, j] == pytest.approx(beta, rel=1e-9)
            assert res.p[g, j] == pytest.approx(p, rel=1e-8)


def test_fwl_without_covariates_is_simple_regression() -> None:
    rng = np.random.default_rng(1)
    values = rng.normal(size=(6, 25))
    design, _ = design_matrix(None, [], 25)
    res = ols_fwl(values, np.array([2]), design)
    slope, intercept, r, p, se = sps.linregress(values[3], values[2])
    assert res.stat[3, 0] == pytest.approx(slope) and res.p[3, 0] == pytest.approx(p)


def test_design_matrix_errors() -> None:
    table = pd.DataFrame({"age": [1.0, np.nan, 3.0]})
    with pytest.raises(DataError, match="missing"):
        design_matrix(table, ["age"], 3)
    with pytest.raises(DataError, match="not in"):
        design_matrix(table, ["sex"], 3)
    with pytest.raises(DataError, match="need a sample table"):
        design_matrix(None, ["age"], 3)
    collinear = pd.DataFrame({"a": [1.0, 2.0, 3.0, 4.0], "b": [2.0, 4.0, 6.0, 8.0]})
    with pytest.raises(DataError, match="rank deficient"):
        design_matrix(collinear, ["a", "b"], 4)


def test_statsmodels_agrees_with_fwl() -> None:
    sm = pytest.importorskip("statsmodels.api")
    rng = np.random.default_rng(2)
    values = rng.normal(size=(5, 30))
    table = pd.DataFrame({"batch": rng.choice(["x", "y"], 30)})
    design, _ = design_matrix(table, ["batch"], 30)
    res = ols_fwl(values, np.array([1]), design)
    x = np.column_stack([design, values[4]])
    fit = sm.OLS(values[1], x).fit()
    assert res.stat[4, 0] == pytest.approx(fit.params[-1]) and res.p[4, 0] == pytest.approx(
        fit.pvalues[-1]
    )


def test_glm_gaussian_equals_ols_and_formula_requires_glm() -> None:
    pytest.importorskip("statsmodels")
    rng = np.random.default_rng(3)
    values = rng.normal(size=(4, 20))
    table = pd.DataFrame({"age": rng.normal(size=20)})
    ols = associate(
        values, np.array([0]), AssociationConfig(method="ols", covariates=["age"]), table
    )
    glm = associate(
        values, np.array([0]), AssociationConfig(method="glm", covariates=["age"]), table
    )
    np.testing.assert_allclose(ols.stat[1:, 0], glm.stat[1:, 0], rtol=1e-6)
    np.testing.assert_allclose(ols.p[1:, 0], glm.p[1:, 0], rtol=1e-4)
    with pytest.raises(ConfigError):
        associate(
            values, np.array([0]), AssociationConfig(method="ols", formula="age:x_node"), table
        )
    inter = associate(
        values,
        np.array([0]),
        AssociationConfig(method="glm", covariates=["age"], formula="x_node:age"),
        table,
    )
    assert np.isfinite(inter.stat[1:, 0]).all()


def test_glm_reads_covariates_whatever_their_names() -> None:
    pytest.importorskip("statsmodels")
    rng = np.random.default_rng(5)
    values = rng.normal(size=(5, 30))
    # names that are not valid Python names, as they come from R tables and spreadsheets
    table = pd.DataFrame(
        {
            "t.rfs": rng.normal(size=30),
            "tumour grade": rng.choice(["low", "high", "mid"], 30),
            "2nd-batch": rng.integers(0, 2, 30).astype(float),
        }
    )
    covariates = ["t.rfs", "tumour grade", "2nd-batch"]
    seeds = np.array([0, 1])
    ols = associate(values, seeds, AssociationConfig(method="ols", covariates=covariates), table)
    glm = associate(values, seeds, AssociationConfig(method="glm", covariates=covariates), table)
    assert np.isfinite(glm.stat).all()
    off_diagonal = np.ones_like(glm.stat, dtype=bool)
    off_diagonal[seeds, [0, 1]] = False  # a seed against itself: a perfect fit, p = 0
    np.testing.assert_allclose(ols.stat[off_diagonal], glm.stat[off_diagonal], rtol=1e-6)
    np.testing.assert_allclose(ols.p[off_diagonal], glm.p[off_diagonal], rtol=1e-4)


# the fits of this test fail on purpose, and statsmodels warns on its way to each failure
@pytest.mark.filterwarnings("ignore")
def test_glm_failures_are_errors_not_empty_networks(caplog: pytest.LogCaptureFixture) -> None:
    pytest.importorskip("statsmodels")
    rng = np.random.default_rng(6)
    values = rng.normal(size=(4, 25))
    table = pd.DataFrame({"age": rng.normal(size=25)})
    seeds = np.array([0])
    # a formula term that names no column is a configuration error
    with pytest.raises(ConfigError, match="no_such_column"):
        associate(values, seeds, AssociationConfig(method="glm", formula="no_such_column"), table)
    # a family that cannot fit these values at all (poisson on negative values): the run stops
    with pytest.raises(DataError, match="every one of the 4 GLM fits failed"):
        associate(values - 5.0, seeds, AssociationConfig(method="glm", glm_family="poisson"), table)
    # a family that fails on some pairs only: those get no association, and the log says so
    with caplog.at_level("WARNING", logger="search"):
        partial = associate(
            values, seeds, AssociationConfig(method="glm", glm_family="gamma"), table
        )
    assert np.isnan(partial.stat).sum() == 1 and np.isfinite(partial.stat).sum() == 3
    assert "1 of 4 fits failed" in caplog.text
    # the checks of the ols path hold for glm too: a covariate with a missing value
    holed = table.copy()
    holed.loc[3, "age"] = np.nan
    with pytest.raises(DataError, match="missing"):
        associate(values, seeds, AssociationConfig(method="glm", covariates=["age"]), holed)


def test_glm_gives_no_association_to_genes_the_engine_cannot_use() -> None:
    pytest.importorskip("statsmodels")
    rng = np.random.default_rng(7)
    values = rng.normal(size=(5, 25))
    values[2, 4] = np.nan  # a missing value
    values[3] = 1.5  # no variance
    for method in ("ols", "glm"):
        result = associate(values, np.array([0]), AssociationConfig(method=method), None)
        assert np.isnan(result.stat[[2, 3], 0]).all(), method
        assert np.isnan(result.p[[2, 3], 0]).all(), method
        assert np.isfinite(result.stat[[1, 4], 0]).all(), method


def test_planted_batch_effect_is_removed_by_the_covariate() -> None:
    rng = np.random.default_rng(4)
    n_genes, n = 120, 80
    batch = np.repeat([0, 1], n // 2)
    values = rng.normal(size=(n_genes, n))
    factor = rng.normal(size=n)
    values[:20] += 1.5 * factor[None, :]  # a real module, seeds inside
    values[20:60] += 3.0 * batch[None, :]  # genes shifted by batch only
    values[:5] += 3.0 * batch[None, :]  # the seeds carry the batch shift too
    table = pd.DataFrame({"batch": batch.astype(str)}, index=[f"s{j}" for j in range(n)])
    expr = ExpressionSet(
        values,
        pd.Index([f"g{i}" for i in range(n_genes)]),
        pd.Index(table.index),
        sample_table=table,
    )
    seeds = [f"g{i}" for i in range(5)]
    base = {
        "edge.rule": "p_threshold",
        "edge.alpha": 0.05,
        "edge.bonferroni": "n_genes",
        "membership.mode": "binary",
    }
    naive = run_core(expr, seeds, RunConfig().with_overrides({**base, "association.method": "ols"}))
    adjusted = run_core(
        expr,
        seeds,
        RunConfig().with_overrides(
            {**base, "association.method": "ols", "association.covariates": ["batch"]}
        ),
    )
    batch_only = [f"g{i}" for i in range(20, 60)]
    module = [f"g{i}" for i in range(5, 20)]
    naive_table = naive.gene_table().set_index("gene")
    adjusted_table = adjusted.gene_table().set_index("gene")
    assert naive_table.loc[batch_only, "n_edges"].mean() > 2  # spurious edges through the batch
    assert adjusted_table.loc[batch_only, "n_edges"].mean() < 0.5
    assert adjusted_table.loc[module, "n_edges"].mean() > 3
    background = [f"g{i}" for i in range(60, 120)]
    assert (
        auc(
            adjusted_table.loc[module, "C_raw"].to_numpy(),
            adjusted_table.loc[background, "C_raw"].to_numpy(),
        )
        > 0.95
    )
    assert adjusted.coef_table is not None and "batch_1" in adjusted.coef_table.columns


def test_bootstrap_keeps_covariates_aligned() -> None:
    rng = np.random.default_rng(5)
    n = 30
    table = pd.DataFrame({"age": rng.normal(size=n)}, index=[f"s{j}" for j in range(n)])
    expr = ExpressionSet(
        rng.normal(size=(8, n)),
        pd.Index([f"g{i}" for i in range(8)]),
        pd.Index(table.index),
        sample_table=table,
    )
    cfg = RunConfig().with_overrides(
        {"association.method": "ols", "association.covariates": ["age"]}
    )
    idx = rng.integers(0, n, n)
    via_idx = run_core(expr, ["g0", "g1"], cfg, sample_idx=idx)
    direct = run_core(expr.take_samples(idx), ["g0", "g1"], cfg)
    np.testing.assert_allclose(via_idx.stat, direct.stat, equal_nan=True)
