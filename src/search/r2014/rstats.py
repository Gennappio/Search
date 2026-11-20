"""R's ``sd``, ``cor``, ``round`` and ``p.adjust``, with the arithmetic in R's order and in R's
precision.

The aim is the same double, bit for bit, as R gives. Two programs that both compute "a
correlation" differ in the last digit as soon as they add the terms in a different order or
keep the running sum in a different precision, and the package then compares the result with a
hard threshold. So the sums here are accumulated one sample at a time, as the loops of
``src/library/stats/src/cov.c`` do, never through a matrix product or ``numpy.sum`` (which
adds in pairs), and in the accumulator of the platform being reproduced
(:mod:`search.r2014.platform`): a double where R's ``long double`` is one, the 80-bit type of
x86 (:mod:`search.r2014.x87`) where it is that.

Spearman and Kendall add up whole numbers or multiples of 1/4, which are exact in either
precision; for them only the last steps (a division, a square root) depend on the platform.
"""

from __future__ import annotations

import math
from typing import Literal

import numpy as np
from numpy.typing import NDArray
from scipy.stats import rankdata

from search.r2014 import platform, x87
from search.r2014.rsem import NA_REAL, FloatArray, RError, UnportedState

CorMethod = Literal["pearson", "spearman", "kendall"]
CorUse = Literal["everything", "all.obs", "complete.obs", "na.or.complete", "pairwise.complete.obs"]


def _mean_two_pass(a: FloatArray, ranks: bool = False) -> FloatArray:
    """The ``MEAN`` macro of cov.c for every row of ``a``: the mean, then a correction by the
    mean of the residuals. Rows are variables, columns are observations.

    ``ranks`` says that the rows are ranks: their mean, (n + 1) / 2, and every sum of them are
    exact in double precision, so the accumulator of the platform makes no difference.
    """
    n = a.shape[1]
    if platform.current().extended and not ranks:
        out = np.empty(a.shape[0])
        count = x87.from_int(n)
        for i, row in enumerate(a):
            if not np.isfinite(row).all():
                out[i] = float(np.sum(row))  # an infinite value: the mean is it, or NaN
                continue
            values = [x87.from_double(float(v)) for v in row]
            total_ext = x87.ZERO
            for v in values:
                total_ext = x87.add(total_ext, v)
            tmp_ext = x87.div(total_ext, count)
            resid_ext = x87.ZERO
            for v in values:
                resid_ext = x87.add(resid_ext, x87.sub(v, tmp_ext))
            out[i] = x87.to_double(x87.add(tmp_ext, x87.div(resid_ext, count)))
        return out
    total = np.zeros(a.shape[0])
    for k in range(n):
        total = total + a[:, k]
    with np.errstate(invalid="ignore", divide="ignore"):
        tmp = total / n
        resid = np.zeros(a.shape[0])
        for k in range(n):
            resid = resid + (a[:, k] - tmp)
        refined = tmp + resid / n
    return np.where(np.isfinite(tmp), refined, tmp)


def _sum_products(a: FloatArray, am: FloatArray, b: FloatArray, bm: FloatArray) -> FloatArray:
    """For every row i: sum over k of (a[i, k] - am[i]) * (b[i, k] - bm[i]), k in order, as a
    double. Used where the sum is exact or the platform accumulates in doubles."""
    total = np.zeros(a.shape[0])
    with np.errstate(invalid="ignore", over="ignore"):
        for k in range(a.shape[1]):
            total = total + (a[:, k] - am) * (b[:, k] - bm)
    return total


def _sum_products_extended(a: FloatArray, am: float, b: FloatArray, bm: float) -> x87.Ext:
    """The same sum for one pair of rows, in 80-bit precision, not yet rounded to a double."""
    am_ext = x87.from_double(am)
    bm_ext = x87.from_double(bm)
    total = x87.ZERO
    for p, q in zip(a, b, strict=True):
        dp = x87.sub(x87.from_double(float(p)), am_ext)
        dq = x87.sub(x87.from_double(float(q)), bm_ext)
        total = x87.add(total, x87.mul(dp, dq))
    return total


def _finish(total: FloatArray, n1: int, root: bool) -> FloatArray:
    """``(double)(sum / n1)`` or ``(double)sqrtl(sum / n1)`` for sums held exactly in doubles.

    In doubles this is one division (and a square root). In 80 bits the quotient is rounded
    to 64 bits first and to 53 afterwards, which is not always the same double.
    """
    with np.errstate(invalid="ignore"):
        if not platform.current().extended:
            return np.asarray(np.sqrt(total / n1) if root else total / n1, dtype=np.float64)
        out = np.empty(total.shape)
        cache: dict[float, float] = {}
        divisor = x87.from_int(n1)
        for index, value in np.ndenumerate(total):
            v = float(value)
            if v not in cache:
                if not math.isfinite(v) or (root and v < 0):
                    cache[v] = v if math.isinf(v) and v > 0 else float("nan")
                else:
                    quotient = x87.div(x87.from_double(v), divisor)
                    cache[v] = x87.to_double(x87.sqrt(quotient) if root else quotient)
            out[index] = cache[v]
        return out


def row_sd(values: FloatArray) -> FloatArray:
    """``apply(values, 1, sd)``. A row with a missing value, or a single observation, gives
    ``NA``."""
    g, n = values.shape
    out = np.full(g, NA_REAL)
    if n <= 1:
        return out
    ok = ~np.isnan(values).any(axis=1)
    if not ok.any():
        return out
    v = values[ok]
    m = _mean_two_pass(v)
    with np.errstate(invalid="ignore"):
        if platform.current().extended:
            n1 = x87.from_int(n - 1)
            variance = np.array(
                [
                    x87.to_double(
                        x87.div(_sum_products_extended(row, float(mean), row, float(mean)), n1)
                    )
                    if np.isfinite(row).all()
                    else float("nan")
                    for row, mean in zip(v, m, strict=True)
                ]
            )
            out[ok] = np.sqrt(variance)  # sd() takes the root of the double variance
        else:
            out[ok] = np.sqrt(_finish(_sum_products(v, m, v, m), n - 1, root=False))
    return out


def rank_average(values: FloatArray) -> FloatArray:
    """``rank(values, na.last = "keep")``: average ranks for ties, missing values stay
    missing."""
    out = np.full(len(values), NA_REAL)
    ok = ~np.isnan(values)
    if ok.any():
        out[ok] = rankdata(values[ok], method="average")
    return out


def _rank_rows(a: FloatArray) -> FloatArray:
    """``Rank`` inside ``cor``: every variable ranked over the observations; with a single
    observation R uses ``row(u)``, a matrix of ones."""
    if a.shape[1] <= 1:
        return np.ones_like(a)
    return np.vstack([rank_average(row) for row in a]) if a.shape[0] else a.copy()


def _clamp(v: FloatArray) -> FloatArray:
    """A correlation that rounding pushed past 1 is 1. R 3.1.0 did not do the same at -1."""
    with np.errstate(invalid="ignore"):
        if platform.current().clamps_at_minus_one:
            return np.where(v >= 1.0, 1.0, np.where(v <= -1.0, -1.0, v))
        return np.where(v > 1.0, 1.0, v)


def _kendall_sign_sum(x: FloatArray, ys: FloatArray) -> FloatArray:
    """Sum over all ordered pairs of observations of sign(dx) * sign(dy), for one ``x`` and
    every row of ``ys``. Whole numbers: exact."""
    sx = np.sign(x[:, None] - x[None, :])
    out = np.empty(ys.shape[0])
    for start in range(0, ys.shape[0], 512):
        block = ys[start : start + 512]
        sy = np.sign(block[:, :, None] - block[:, None, :])
        out[start : start + 512] = np.einsum("kl,gkl->g", sx, sy)
    return out


def _kendall_untied(a: FloatArray) -> FloatArray:
    """For every row: the number of ordered pairs of observations with different values."""
    out = np.empty(a.shape[0])
    for start in range(0, a.shape[0], 512):
        block = a[start : start + 512]
        out[start : start + 512] = (block[:, :, None] != block[:, None, :]).sum(axis=(1, 2))
    return out


def _cor_no_missing_rule(
    xs: FloatArray, ys: FloatArray, kendall: bool, ranks: bool = False
) -> FloatArray:
    """``cov_na_2`` of cov.c with ``cor = TRUE``: the result for a variable with a missing
    value is ``NA``; no observation is dropped. ``xs`` is K x n, ``ys`` is G x n. ``ranks``:
    the rows are ranks (Spearman), whose sums are exact in double precision."""
    k_x, n = xs.shape
    g = ys.shape[0]
    ans = np.full((k_x, g), NA_REAL)
    if n <= 1:
        return ans
    ok_x = ~np.isnan(xs).any(axis=1)
    ok_y = ~np.isnan(ys).any(axis=1)
    if not ok_x.any() or not ok_y.any():
        return ans
    yv = ys[ok_y]
    xv = xs[ok_x]
    n1 = n - 1
    extended_sums = platform.current().extended and not kendall and not ranks
    cov = np.empty((xv.shape[0], yv.shape[0]))
    with np.errstate(invalid="ignore", divide="ignore"):
        if kendall:
            # whole numbers: the covariance is the sum itself, the deviations its square root
            sd_x = _finish(_kendall_untied(xv), 1, root=True)
            sd_y = _finish(_kendall_untied(yv), 1, root=True)
            for i in range(xv.shape[0]):
                cov[i] = _kendall_sign_sum(xv[i], yv)
        elif extended_sums:
            xm = _mean_two_pass(xv)
            ym = _mean_two_pass(yv)
            divisor = x87.from_int(n1)

            def deviation(row: FloatArray, mean: float) -> float:
                if not np.isfinite(row).all():
                    return float("nan")
                total = _sum_products_extended(row, mean, row, mean)
                return x87.to_double(x87.sqrt(x87.div(total, divisor)))

            sd_x = np.array([deviation(row, float(mean)) for row, mean in zip(xv, xm, strict=True)])
            sd_y = np.array([deviation(row, float(mean)) for row, mean in zip(yv, ym, strict=True)])
            for i in range(xv.shape[0]):
                for j in range(yv.shape[0]):
                    if np.isfinite(xv[i]).all() and np.isfinite(yv[j]).all():
                        total = _sum_products_extended(xv[i], float(xm[i]), yv[j], float(ym[j]))
                        cov[i, j] = x87.to_double(x87.div(total, divisor))
                    else:
                        cov[i, j] = float("nan")
        else:
            xm = _mean_two_pass(xv, ranks)
            ym = _mean_two_pass(yv, ranks)
            sd_x = _finish(_sum_products(xv, xm, xv, xm), n1, root=True)
            sd_y = _finish(_sum_products(yv, ym, yv, ym), n1, root=True)
            for i in range(xv.shape[0]):
                xi = np.broadcast_to(xv[i], yv.shape)
                products = _sum_products(xi, np.full(yv.shape[0], xm[i]), yv, ym)
                cov[i] = _finish(products, n1, root=False)
        r = _clamp(cov / (sd_x[:, None] * sd_y[None, :]))
        # a standard deviation of zero: NA, with a warning in R
        r = np.where((sd_x[:, None] == 0) | (sd_y[None, :] == 0), NA_REAL, r)
    ans[np.ix_(ok_x, ok_y)] = r
    return ans


def _cor_pair(x: FloatArray, y: FloatArray, method: str) -> float:
    """One pair under ``pairwise.complete.obs``.

    Pearson goes through ``cov_pairwise2``: the observations present in both, and a
    single-pass mean. Spearman and Kendall are done by R code: the complete pairs are ranked
    and passed to the C routine with ``all.obs``.
    """
    ok = ~(np.isnan(x) | np.isnan(y))
    nobs = int(ok.sum())
    if method != "pearson":
        if nobs == 0:
            return NA_REAL
        xr = rank_average(x[ok])[None, :]
        yr = rank_average(y[ok])[None, :]
        rule = _cor_no_missing_rule(xr, yr, kendall=method == "kendall", ranks=method == "spearman")
        return float(rule[0, 0])
    if nobs < 2:
        return NA_REAL
    xv = x[ok]
    yv = y[ok]
    if platform.current().extended:
        return _cor_pair_extended(xv, yv)
    xsd = ysd = xysum = 0.0
    xmean = ymean = 0.0
    for a, b in zip(xv, yv, strict=True):
        xmean += float(a)
        ymean += float(b)
    xmean /= nobs
    ymean /= nobs
    for a, b in zip(xv, yv, strict=True):
        xm = float(a) - xmean
        ym = float(b) - ymean
        xysum += xm * ym
        xsd += xm * xm
        ysd += ym * ym
    if xsd == 0.0 or ysd == 0.0:
        return NA_REAL
    n1 = nobs - 1
    xsd /= n1
    ysd /= n1
    xysum /= n1
    xsd = math.sqrt(xsd)
    ysd = math.sqrt(ysd)
    value = xysum / (xsd * ysd)
    if value >= 1.0:
        return 1.0
    return -1.0 if value <= -1.0 and platform.current().clamps_at_minus_one else value


def _cor_pair_extended(xv: FloatArray, yv: FloatArray) -> float:
    """``COV_PAIRWISE_BODY`` for Pearson with every variable in 80 bits, as on x86: the means,
    the sums, the division by the two deviations, and one rounding to a double at the end."""
    if not (np.isfinite(xv).all() and np.isfinite(yv).all()):
        raise UnportedState("infinite values in a pairwise Pearson correlation")
    nobs = len(xv)
    xs_ext = [x87.from_double(float(v)) for v in xv]
    ys_ext = [x87.from_double(float(v)) for v in yv]
    xmean = ymean = x87.ZERO
    for a, b in zip(xs_ext, ys_ext, strict=True):
        xmean = x87.add(xmean, a)
        ymean = x87.add(ymean, b)
    count = x87.from_int(nobs)
    xmean = x87.div(xmean, count)
    ymean = x87.div(ymean, count)
    xsd = ysd = total = x87.ZERO
    for a, b in zip(xs_ext, ys_ext, strict=True):
        xm = x87.sub(a, xmean)
        ym = x87.sub(b, ymean)
        total = x87.add(total, x87.mul(xm, ym))
        xsd = x87.add(xsd, x87.mul(xm, xm))
        ysd = x87.add(ysd, x87.mul(ym, ym))
    if xsd[0] == 0 or ysd[0] == 0:
        return NA_REAL
    n1 = x87.from_int(nobs - 1)
    xsd = x87.div(xsd, n1)
    ysd = x87.div(ysd, n1)
    total = x87.div(total, n1)
    total = x87.div(total, x87.mul(x87.sqrt(xsd), x87.sqrt(ysd)))
    one = x87.from_int(1)
    if x87.less(one, total):
        total = one
    if platform.current().clamps_at_minus_one and x87.less(total, x87.neg(one)):
        total = x87.neg(one)
    return x87.to_double(total)


def r_cor(
    xs: FloatArray, ys: FloatArray, method: str = "spearman", use: str = "everything"
) -> FloatArray:
    """``cor(x, y, method, use)`` where the variables are the rows of ``xs`` (K x n) and of
    ``ys`` (G x n). Returns K x G.

    The method and the ``use`` are matched as R matches them: a unique prefix is enough.
    """
    method_full = _match_arg(method, ("pearson", "kendall", "spearman"), "method")
    use_full = _match_arg(
        use,
        ("all.obs", "complete.obs", "pairwise.complete.obs", "everything", "na.or.complete"),
        "use",
    )
    kendall = method_full == "kendall"
    if use_full == "pairwise.complete.obs":
        if method_full != "pearson" and (xs.size == 0 or ys.size == 0):
            raise RError("both 'x' and 'y' must be non-empty")
        out = np.empty((xs.shape[0], ys.shape[0]))
        for i in range(xs.shape[0]):
            for j in range(ys.shape[0]):
                out[i, j] = _cor_pair(xs[i], ys[j], method_full)
        return out
    if use_full in ("complete.obs", "na.or.complete"):
        complete = ~(np.isnan(xs).any(axis=0) | np.isnan(ys).any(axis=0))
        if not complete.any():
            if use_full == "complete.obs":
                raise RError("no complete element pairs")
            return np.full((xs.shape[0], ys.shape[0]), NA_REAL)
        xs = xs[:, complete]
        ys = ys[:, complete]
    elif use_full == "all.obs" and (np.isnan(xs).any() or np.isnan(ys).any()):
        raise RError("missing observations in cov/cor")
    if method_full == "spearman":
        xs = _rank_rows(xs)
        ys = _rank_rows(ys)
    return _cor_no_missing_rule(xs, ys, kendall, ranks=method_full == "spearman")


def cor_keeps_names(xs: FloatArray, ys: FloatArray, method: str, use: str) -> bool:
    """Whether the matrix ``cor(x, y)`` returns carries the names of the variables of ``y``.

    It does not in one case. For Spearman and Kendall under ``complete.obs`` and
    ``na.or.complete`` R drops the incomplete observations and ranks what is left with its
    ``Rank`` helper, which for a matrix of a single row returns ``row(u)``: a matrix of ones
    without names. With exactly one complete observation the correlations are all ``NA`` and
    the names are gone.
    """
    method_full = _match_arg(method, ("pearson", "kendall", "spearman"), "method")
    use_full = _match_arg(
        use,
        ("all.obs", "complete.obs", "pairwise.complete.obs", "everything", "na.or.complete"),
        "use",
    )
    if method_full == "pearson" or use_full not in ("complete.obs", "na.or.complete"):
        return True
    complete = ~(np.isnan(xs).any(axis=0) | np.isnan(ys).any(axis=0))
    return int(complete.sum()) != 1


def _match_arg(value: str, choices: tuple[str, ...], what: str) -> str:
    """``pmatch`` / ``match.arg``: an exact match, or a unique prefix."""
    if value in choices:
        return value
    hits = [c for c in choices if c.startswith(value)] if value else []
    if len(hits) != 1:
        if what == "use":
            raise RError("invalid 'use' argument")
        # dQuote(), in a UTF-8 locale; in the C locale R writes plain double quotes
        raise RError("'arg' should be one of " + ", ".join(f"“{c}”" for c in choices))
    return hits[0]


def r_round(values: FloatArray, digits: int) -> FloatArray:
    """``round(values, digits)`` of the platform's R.

    From R 4.0.0 (``private_fround`` in ``nmath/fround.c``): of the two neighbouring multiples
    of 10^-digits, the one nearer to the value, and the even one when they are equally near.
    Before: the whole part plus ``rint(fraction * 10^digits) / 10^digits``, computed in ``long
    double``. The two disagree on values that sit close to a half-way point, such as 0.003125,
    which is 0.00313 for the first and 0.00312 for the second.
    """
    modern = platform.current().round_to_representable
    if not modern and not platform.current().extended:
        raise UnportedState("the round() of R before 4.0.0 without extended precision")
    out = np.array(values, dtype=np.float64, copy=True)
    for i, x in enumerate(out):
        v = float(x)
        if math.isnan(v) or math.isinf(v) or v == 0.0:
            continue
        if modern:
            magnitude = _round_positive(abs(v), digits)
        else:
            magnitude = _round_positive_old(abs(v), digits)
        out[i] = -magnitude if v < 0 else magnitude
    return out


_DBL_DIG = 15


def _round_positive(x: float, dig: int) -> float:
    # more digits than a double holds: nothing to round
    e10 = int(math.floor(math.log10(x)))
    if dig > _DBL_DIG - 1 - e10 + 1:
        return x
    p10 = 10.0**dig if dig >= 0 else 1.0 / (10.0 ** (-dig))
    x10 = p10 * x
    i10 = math.floor(x10)
    xd = i10 / p10
    xu = math.ceil(x10) / p10
    du = xu - x
    dd = x - xd
    if dd < du or (dd == du and math.fmod(i10, 2.0) == 0.0):
        return xd
    return xu


def _round_positive_old(x: float, dig: int) -> float:
    """``fround`` of R 3.1.0 for ``dig > 0``, in 80-bit precision:
    ``(double)(intx + rint((double)((x - intx) * pow10)) / pow10)``."""
    if dig <= 0:
        raise UnportedState("round() to no decimals or to tens, as R before 4.0.0 did it")
    pow10 = x87.from_double(10.0**dig)
    whole = x87.from_double(float(math.floor(x)))
    fraction = x87.to_double(x87.mul(x87.sub(x87.from_double(x), whole), pow10))
    nearest = float(round(fraction))  # nearbyint: to nearest, ties to even
    return x87.to_double(x87.add(whole, x87.div(x87.from_double(nearest), pow10)))


def p_adjust(p: FloatArray, method: str) -> FloatArray:
    """``p.adjust(p, method)``. The number of comparisons is the number of p-values that are
    not missing: the default ``n = length(p)`` is evaluated after R has dropped them."""
    methods = ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
    method = _match_arg(method, methods, "method")
    if method == "fdr":
        method = "BH"
    p0 = np.array(p, dtype=np.float64, copy=True)
    nna = ~np.isnan(p0)
    pv = p0[nna]
    lp = n = len(pv)
    if n <= 1:
        return p0
    if n == 2 and method == "hommel":
        method = "hochberg"
    if method == "bonferroni":
        adj = np.minimum(1.0, n * pv)
    elif method == "holm":
        o = np.argsort(pv, kind="stable")
        i = np.arange(1, lp + 1)
        adj = _unorder(np.minimum(1.0, np.maximum.accumulate((n + 1 - i) * pv[o])), o)
    elif method in ("hochberg", "BH", "BY"):
        o = _order_decreasing(pv)
        i = np.arange(lp, 0, -1)
        if method == "hochberg":
            scaled = (n + 1 - i) * pv[o]
        elif method == "BH":
            scaled = n / i * pv[o]
        else:
            # q <- sum(1/(1L:n)): sum() accumulates in long double
            if platform.current().extended:
                total_ext = x87.ZERO
                for k in range(1, n + 1):
                    total_ext = x87.add(total_ext, x87.from_double(1.0 / k))
                harmonic = x87.to_double(total_ext)
            else:
                harmonic = 0.0
                for k in range(1, n + 1):
                    harmonic += 1.0 / k
            scaled = harmonic * n / i * pv[o]
        adj = _unorder(np.minimum(1.0, np.minimum.accumulate(scaled)), o)
    elif method == "hommel":
        o = np.argsort(pv, kind="stable")
        ps = pv[o]
        i = np.arange(1, n + 1)
        q = np.full(n, np.min(n * ps / i))
        pa = q.copy()
        for m in range(n - 1, 1, -1):
            i1 = np.arange(0, n - m + 1)
            i2 = np.arange(n - m + 1, n)
            q1 = np.min(m * ps[i2] / np.arange(2, m + 1))
            q[i1] = np.minimum(m * ps[i1], q1)
            q[i2] = q[n - m]
            pa = np.maximum(pa, q)
        adj = _unorder(np.maximum(pa, ps), o)
    else:
        adj = pv
    p0[nna] = adj
    return p0


def _order_decreasing(v: FloatArray) -> NDArray[np.intp]:
    """``order(v, decreasing = TRUE)``: ties keep their original order."""
    return np.argsort(-v, kind="stable")


def _unorder(sorted_values: FloatArray, order: NDArray[np.intp]) -> FloatArray:
    out = np.empty_like(sorted_values)
    out[order] = sorted_values
    return out
