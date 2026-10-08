import numpy as np
import pytest
from scipy import stats as sps

from search.stats import p_adjust, rank01, row_sds, t_approx_pvalue


def test_rank01_basic_and_ties() -> None:
    np.testing.assert_allclose(rank01(np.array([3.0, 1.0, 2.0])), [1.0, 0.0, 0.5])
    np.testing.assert_allclose(rank01(np.array([1.0, 1.0, 2.0])), [0.25, 0.25, 1.0])
    np.testing.assert_allclose(
        rank01(np.array([3.0, 1.0, 2.0]), higher_is_better=False), [0.0, 1.0, 0.5]
    )


def test_rank01_missing_and_degenerate() -> None:
    out = rank01(np.array([np.nan, 2.0, 1.0]))
    assert np.isnan(out[0]) and out[1] == 1.0 and out[2] == 0.0
    assert np.isnan(rank01(np.array([5.0]))).all()
    assert np.isnan(rank01(np.array([np.nan, 7.0]))).all()


def test_t_approx_pvalue_matches_manual_formula() -> None:
    rho = np.array([0.0, 0.5, -0.9, np.nan])
    n = 20
    p = t_approx_pvalue(rho, n)
    t = rho[:3] * np.sqrt((n - 2) / (1 - rho[:3] ** 2))
    expected = 2 * sps.t.sf(np.abs(t), df=n - 2)
    np.testing.assert_allclose(p[:3], expected, rtol=1e-12)
    assert p[0] == pytest.approx(1.0)
    assert np.isnan(p[3])
    assert np.isnan(t_approx_pvalue(np.array([0.5]), 2)).all()


def test_p_adjust_known_values() -> None:
    p = np.array([0.01, 0.02, 0.03, 0.04])
    np.testing.assert_allclose(p_adjust(p, "bonferroni"), [0.04, 0.08, 0.12, 0.16])
    np.testing.assert_allclose(p_adjust(p, "bh"), [0.04, 0.04, 0.04, 0.04])
    np.testing.assert_allclose(p_adjust(np.array([0.5, 0.001, 0.02]), "bh"), [0.5, 0.003, 0.03])
    with_nan = p_adjust(np.array([0.01, np.nan, 0.5]), "bonferroni")
    np.testing.assert_allclose(with_nan[[0, 2]], [0.02, 1.0])
    assert np.isnan(with_nan[1])
    with pytest.raises(ValueError):
        p_adjust(p, "holm")


def test_row_sds() -> None:
    x = np.array([[1.0, 2.0, 3.0], [5.0, 5.0, 5.0], [np.nan, 1.0, 3.0]])
    np.testing.assert_allclose(row_sds(x), [1.0, 0.0, np.sqrt(2.0)])
