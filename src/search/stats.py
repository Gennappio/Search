"""Small statistical helpers shared by every stage."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from scipy import stats as sps

FloatArray = NDArray[np.float64]


def rank_rows(x: FloatArray) -> FloatArray:
    """Average ranks along each row; NaN entries stay NaN."""
    x = np.asarray(x, dtype=np.float64)
    ranks: FloatArray = sps.rankdata(x, axis=1, method="average", nan_policy="omit")
    return ranks


def rank01(x: FloatArray, higher_is_better: bool = True) -> FloatArray:
    """Fractional rank in [0, 1]: lowest value 0, highest 1, ties get the average rank.

    Non-finite values are left as NaN and do not take part in the ranking. With one or zero
    finite values every entry is NaN, as in the original script.
    """
    x = np.asarray(x, dtype=np.float64)
    out = np.full(x.shape, np.nan)
    ok = np.isfinite(x)
    n_ok = int(ok.sum())
    if n_ok <= 1:
        return out
    values = x[ok] if higher_is_better else -x[ok]
    r = sps.rankdata(values, method="average")
    out[ok] = (r - 1.0) / (n_ok - 1.0)
    return out


def t_approx_pvalue(rho: FloatArray, n: int) -> FloatArray:
    """Two-sided p-value of a correlation via the t approximation, n - 2 degrees of freedom."""
    rho = np.asarray(rho, dtype=np.float64)
    if n - 2 <= 0:
        return np.full(rho.shape, np.nan)
    r = np.clip(rho, -0.999999, 0.999999)
    with np.errstate(invalid="ignore", divide="ignore"):
        t = r * np.sqrt((n - 2) / np.maximum(1e-12, 1.0 - r * r))
    p = 2.0 * sps.t.sf(np.abs(t), df=n - 2)
    p = np.where(np.isfinite(rho), p, np.nan)
    return np.asarray(p, dtype=np.float64)


def p_adjust(p: FloatArray, method: str) -> FloatArray:
    """Multiple-testing adjustment of a vector of p-values, mirroring R's p.adjust.

    ``bonferroni`` multiplies by the number of non-missing tests; ``bh`` is Benjamini-Hochberg.
    Missing values are ignored and returned as NaN.
    """
    p = np.asarray(p, dtype=np.float64)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    n = int(ok.sum())
    if n == 0:
        return out
    values = p[ok]
    if method == "none":
        out[ok] = values
    elif method == "bonferroni":
        out[ok] = np.minimum(1.0, values * n)
    elif method == "bh":
        order = np.argsort(values)[::-1]  # from the largest p downwards
        ranked = values[order]
        factors = n / np.arange(n, 0, -1)
        adjusted = np.minimum.accumulate(np.minimum(1.0, factors * ranked))
        restored = np.empty(n)
        restored[order] = adjusted
        out[ok] = restored
    else:
        raise ValueError(f"unknown p-value adjustment: {method}")
    return out


def row_sds(x: FloatArray) -> FloatArray:
    """Sample standard deviation of each row, NaN-aware; non-finite results become 0."""
    x = np.asarray(x, dtype=np.float64)
    with np.errstate(invalid="ignore"):
        out = np.nanstd(x, axis=1, ddof=1)
    out[~np.isfinite(out)] = 0.0
    return np.asarray(out, dtype=np.float64)
