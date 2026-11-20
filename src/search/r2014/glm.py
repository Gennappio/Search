"""``nodes.glm`` and ``nodes.fastLm`` of SEARCH 0.0.2.

Both fit, for every seed and every gene of the universe, the model ``seed ~ gene`` (plus the
terms of ``formula.string`` when covariates are given), and put a gene in the neighbourhood of
the seed when the p-value of its coefficient is below ``significance.threshold``.

What is reproduced exactly is the logic around the fits, where the defects are: which genes
enter a neighbourhood, a neighbourhood of one member losing its name (R5), a seed without
neighbours vanishing (R6), covariates used only through ``formula.string`` (D3). The fits
themselves follow R's algorithms (``glm.fit`` with the LINPACK QR of ``lm``, ``summary.glm``;
``fastLm`` of RcppArmadillo) but run on other linear-algebra and distribution routines, so
coefficients and p-values agree with R to rounding error, not bit for bit.

``formula.string`` is R code. The port reads sums of covariate names and their interactions
(``a + b``, ``a:b``, ``a*b``, also with ``x.node``); anything else is refused.
"""

from __future__ import annotations

import math
import re
import warnings
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy import special

from search.r2014.association import AssociationData, universe_rows
from search.r2014.rsem import (
    NA_REAL,
    FloatArray,
    Name,
    NamedVector,
    Names,
    NotReproducible,
    RError,
    RList,
    RMatrix,
    UnportedState,
)
from search.r2014.rstats import p_adjust

_EPS = float(np.finfo(np.float64).eps)


@dataclass
class Covariates:
    """The ``covariates`` data frame: one value per sample for every column. A column of
    numbers is numeric; a column of strings is what R turns into a factor. Missing values are
    ``NA_REAL`` and ``None``."""

    columns: dict[str, list[Any]]

    def rows(self, index: np.ndarray) -> Covariates:
        """``covariates[index, , drop = FALSE]`` with 0-based ``index``."""
        return Covariates(
            {name: [values[i] for i in index] for name, values in self.columns.items()}
        )

    def n_rows(self) -> int:
        return len(next(iter(self.columns.values()))) if self.columns else 0

    def is_numeric(self, name: str) -> bool:
        return all(
            isinstance(v, (int, float)) and not isinstance(v, bool) for v in self.columns[name]
        )


# ---- model formula ----------------------------------------------------------------------------
_VARIABLE = r"[A-Za-z.][A-Za-z0-9._]*|`[^`]+`"


def _parse_terms(
    formula_string: str | None, covariates: Covariates | None
) -> list[tuple[str, ...]]:
    """The terms of ``seed.exp ~ x.node + <formula.string>`` as tuples of variable names, in
    R's order: by number of variables, then as written."""
    text = (
        "x.node" if formula_string is None or covariates is None else f"x.node + {formula_string}"
    )
    variables: list[str] = []
    terms: list[tuple[str, ...]] = []

    def variable(token: str) -> str:
        name = token.strip().strip("`")
        if not re.fullmatch(_VARIABLE, token.strip()):
            raise UnportedState(
                f"formula.string term {token.strip()!r}: only names, ':' and '*' are read"
            )
        if name != "x.node" and (covariates is None or name not in covariates.columns):
            raise RError(f"object '{name}' not found", "eval(predvars, data, env)")
        if name not in variables:
            variables.append(name)
        return name

    def add(term: Sequence[str]) -> None:
        ordered = tuple(sorted(set(term), key=variables.index))
        if ordered not in terms:
            terms.append(ordered)

    for part in text.split("+"):
        if part.strip() == "1":
            continue
        if any(c in part for c in "()-^/|%~"):
            raise UnportedState(
                f"formula.string term {part.strip()!r}: only names, ':' and '*' are read"
            )
        factors = [[variable(v) for v in piece.split(":")] for piece in part.split("*")]
        # a*b*c: every non-empty subset of the factors, joined
        for mask in range(1, 2 ** len(factors)):
            add([v for k, group in enumerate(factors) if mask >> k & 1 for v in group])
    return sorted(terms, key=len)


def _model_matrix(
    y: FloatArray, x_node: FloatArray, covariates: Covariates | None, terms: list[tuple[str, ...]]
) -> tuple[FloatArray, FloatArray, list[str]]:
    """``model.frame`` and ``model.matrix``: the samples with a missing value in a variable of
    the model are dropped (``na.omit``); strings become factors with treatment contrasts."""
    used = [v for v in dict.fromkeys(v for term in terms for v in term) if v != "x.node"]
    keep = ~(np.isnan(y) | np.isnan(x_node))
    data: dict[str, Any] = {}
    for name in used:
        assert covariates is not None
        column = covariates.columns[name]
        if covariates.is_numeric(name):
            values = np.array(column, dtype=np.float64)
            keep &= ~np.isnan(values)
            data[name] = values
        else:
            keep &= np.array([v is not None for v in column], dtype=bool)
            data[name] = column
    yk, xk = y[keep], x_node[keep]
    n = int(keep.sum())
    columns: list[FloatArray] = [np.ones(n)]
    names: list[str] = ["(Intercept)"]
    term_set = set(terms)

    def coding(name: str, by_contrast: bool) -> list[tuple[str, FloatArray]]:
        if name == "x.node":
            return [("x.node", xk)]
        values = data[name]
        if isinstance(values, np.ndarray):
            return [(name, values[keep])]
        kept = [v for v, k in zip(values, keep, strict=True) if k]
        levels = sorted(set(kept))  # sort(unique(x)) in the C locale
        if by_contrast:
            if len(levels) < 2:
                raise RError("contrasts can be applied only to factors with 2 or more levels")
            levels = levels[1:]
        return [
            (f"{name}{level}", np.array([1.0 if v == level else 0.0 for v in kept]))
            for level in levels
        ]

    for term in terms:
        block: list[tuple[str, FloatArray]] = [("", np.ones(n))]
        for name in term:
            rest = tuple(v for v in term if v != name)
            pieces = coding(name, by_contrast=len(rest) == 0 or rest in term_set)
            # the first variable of a term varies fastest
            block = [
                (f"{label}:{piece_label}" if label else piece_label, column * piece)
                for piece_label, piece in pieces
                for label, column in block
            ]
        for label, term_column in block:
            names.append(label)
            columns.append(term_column)
    return yk, np.column_stack(columns), names


# ---- least squares as lm and glm do them ------------------------------------------------------
def _dqrls(
    x: FloatArray, y: FloatArray, tol: float
) -> tuple[FloatArray, FloatArray, int, list[int]]:
    """LINPACK ``dqrdc2`` and the solution step of ``dqrls``: Householder QR with the limited
    pivoting of R (a column whose norm falls below ``tol`` times its original norm is moved to
    the end). Returns the coefficients in pivoted order, the factor, the rank and the pivot."""
    a = np.array(x, dtype=np.float64, copy=True)
    n, p = a.shape
    qraux = np.array([np.linalg.norm(a[:, j]) for j in range(p)])
    work1 = qraux.copy()
    work2 = np.where(qraux == 0.0, 1.0, qraux)
    jpvt = list(range(p))
    lup = min(n, p)
    k = p + 1
    for l in range(1, lup + 1):  # noqa: E741  (the index names of the Fortran source)
        while l < k and qraux[l - 1] < work2[l - 1] * tol:
            a[:, l - 1 :] = np.roll(a[:, l - 1 :], -1, axis=1)
            for array in (qraux, work1, work2):
                array[l - 1 :] = np.roll(array[l - 1 :], -1)
            jpvt = [*jpvt[: l - 1], *jpvt[l:], jpvt[l - 1]]
            k -= 1
        if l == n:
            break
        nrmxl = float(np.linalg.norm(a[l - 1 :, l - 1]))
        if nrmxl != 0.0:
            if a[l - 1, l - 1] != 0.0:
                nrmxl = math.copysign(nrmxl, a[l - 1, l - 1])
            a[l - 1 :, l - 1] /= nrmxl
            a[l - 1, l - 1] += 1.0
            for j in range(l + 1, p + 1):
                t = -float(a[l - 1 :, l - 1] @ a[l - 1 :, j - 1]) / a[l - 1, l - 1]
                a[l - 1 :, j - 1] += t * a[l - 1 :, l - 1]
                if qraux[j - 1] != 0.0:
                    tt = max(1.0 - (abs(a[l - 1, j - 1]) / qraux[j - 1]) ** 2, 0.0)
                    if abs(tt) < 1e-6:
                        qraux[j - 1] = float(np.linalg.norm(a[l:, j - 1]))
                        work1[j - 1] = qraux[j - 1]
                    else:
                        qraux[j - 1] *= math.sqrt(tt)
            qraux[l - 1] = a[l - 1, l - 1]
            a[l - 1, l - 1] = -nrmxl
    rank = min(k - 1, n)
    qty = np.array(y, dtype=np.float64, copy=True)
    for j in range(1, min(rank, n - 1) + 1):
        if qraux[j - 1] != 0.0:
            saved = a[j - 1, j - 1]
            a[j - 1, j - 1] = qraux[j - 1]
            t = -float(a[j - 1 :, j - 1] @ qty[j - 1 :]) / a[j - 1, j - 1]
            qty[j - 1 :] += t * a[j - 1 :, j - 1]
            a[j - 1, j - 1] = saved
    coefficients = np.zeros(p)
    coefficients[:rank] = qty[:rank]
    for j in range(rank, 0, -1):
        if a[j - 1, j - 1] == 0.0:
            break
        coefficients[j - 1] /= a[j - 1, j - 1]
        if j != 1:
            coefficients[: j - 1] -= coefficients[j - 1] * a[: j - 1, j - 1]
    return coefficients, a, rank, jpvt


@dataclass
class _Family:
    name: str
    estimate_dispersion: bool
    linkfun: Any
    linkinv: Any
    mu_eta: Any
    variance: Any
    dev_resids: Any
    mustart: Any
    validmu: Any
    valideta: Any


def _poisson_dev(y: FloatArray, mu: FloatArray) -> FloatArray:
    r = mu.copy()
    positive = y > 0
    with np.errstate(divide="ignore", invalid="ignore"):
        r[positive] = (y * np.log(y / mu) - (y - mu))[positive]
    return np.asarray(2.0 * r, dtype=np.float64)


def _gamma_dev(y: FloatArray, mu: FloatArray) -> FloatArray:
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.asarray(
            -2.0 * (np.log(np.where(y == 0, 1.0, y / mu)) - (y - mu) / mu), dtype=np.float64
        )


def _family(name: str, y: FloatArray) -> _Family:
    """The families ``nodes.glm`` accepts, with their default links, as R defines them."""
    true = lambda v: True  # noqa: E731
    identity = lambda v: v  # noqa: E731
    ones = lambda v: np.ones_like(v)  # noqa: E731
    square = lambda y_, mu: (y_ - mu) ** 2  # noqa: E731
    log_inv = lambda eta: np.maximum(np.exp(eta), _EPS)  # noqa: E731
    positive = lambda mu: bool(np.all(np.isfinite(mu)) and np.all(mu > 0))  # noqa: E731
    if name in ("gaussian", "quasi"):
        return _Family(
            name, True, identity, identity, ones, ones, square, lambda: y.copy(), true, true
        )
    if name in ("poisson", "quasipoisson"):
        if np.any(y < 0):
            label = "Poisson" if name == "poisson" else "quasiPoisson"
            raise RError(f"negative values not allowed for the '{label}' family")
        return _Family(
            name, name == "quasipoisson", np.log, log_inv, log_inv, identity, _poisson_dev,
            lambda: y + 0.1, positive, true,
        )  # fmt: skip
    if name == "Gamma":
        if np.any(y <= 0):
            raise RError("non-positive values not allowed for the 'Gamma' family")
        return _Family(
            name, True, lambda mu: 1.0 / mu, lambda eta: 1.0 / eta, lambda eta: -1.0 / eta**2,
            lambda mu: mu**2, _gamma_dev, lambda: y.copy(), lambda mu: bool(np.all(mu > 0)),
            lambda eta: bool(np.all(np.isfinite(eta)) and np.all(eta != 0)),
        )  # fmt: skip
    if name == "inverse.gaussian":
        if np.any(y <= 0):
            raise RError("positive values only are allowed for the 'inverse.gaussian' family")
        return _Family(
            name, True, lambda mu: 1.0 / mu**2, lambda eta: 1.0 / np.sqrt(eta),
            lambda eta: -1.0 / (2.0 * eta**1.5), lambda mu: mu**3,
            lambda y_, mu: ((y_ - mu) ** 2) / (y_ * mu**2), lambda: y.copy(), true,
            lambda eta: bool(np.all(np.isfinite(eta)) and np.all(eta > 0)),
        )  # fmt: skip
    raise RError("glm.family not supported!")


def _glm_coef_table(
    y: FloatArray, x: FloatArray, family_name: str
) -> tuple[FloatArray, FloatArray]:
    """``glm.fit`` followed by ``summary.glm``: the estimate of every coefficient (NA when
    aliased) and its p-value (NA when aliased)."""
    nobs, nvars = x.shape
    if nobs == 0:
        raise UnportedState("a model without any complete sample")
    fam = _family(family_name, y)
    epsilon, maxit = 1e-8, 25
    mustart = fam.mustart()
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        eta = fam.linkfun(mustart)
        mu = fam.linkinv(eta)
        if not (fam.validmu(mu) and fam.valideta(eta)):
            raise RError("cannot find valid starting values: please specify some")
        devold = float(np.sum(fam.dev_resids(y, mu)))
        coefold: FloatArray | None = None
        coef: FloatArray | None = None
        start = np.zeros(nvars)
        w = np.ones(nobs)
        good = np.ones(nobs, dtype=bool)
        factor = np.zeros((nobs, nvars))
        rank = 0
        pivot = list(range(nvars))
        for _ in range(maxit):
            varmu = fam.variance(mu)
            if np.isnan(varmu).any():
                raise RError("NAs in V(mu)")
            if np.any(varmu == 0):
                raise RError("0s in V(mu)")
            mu_eta_val = fam.mu_eta(eta)
            if np.isnan(mu_eta_val).any():
                raise RError("NAs in d(mu)/d(eta)")
            good = mu_eta_val != 0
            if not good.any():
                break
            z = eta[good] + (y - mu)[good] / mu_eta_val[good]
            w = np.sqrt(mu_eta_val[good] ** 2 / varmu[good])
            fitted, factor, rank, pivot = _dqrls(
                x[good] * w[:, None], z * w, min(1e-7, epsilon / 1000)
            )
            if not np.all(np.isfinite(fitted)):
                break
            if nobs < rank:
                raise RError(f"X matrix has rank {rank}, but only {nobs} observations")
            start = np.zeros(nvars)
            start[pivot] = fitted
            eta = x @ start
            mu = fam.linkinv(eta)
            dev = float(np.sum(fam.dev_resids(y, mu)))
            if not math.isfinite(dev) or not (fam.valideta(eta) and fam.validmu(mu)):
                if coefold is None:
                    raise RError(
                        "no valid set of coefficients has been found: please supply starting values"
                    )
                halvings = 1
                while not math.isfinite(dev) or not (fam.valideta(eta) and fam.validmu(mu)):
                    if halvings > maxit:
                        raise RError(
                            "inner loop 1; cannot correct step size"
                            if not math.isfinite(dev)
                            else "inner loop 2; cannot correct step size"
                        )
                    halvings += 1
                    start = (start + coefold) / 2.0
                    eta = x @ start
                    mu = fam.linkinv(eta)
                    dev = float(np.sum(fam.dev_resids(y, mu)))
            if abs(dev - devold) / (abs(dev) + 0.1) < epsilon:
                coef = start
                break
            devold = dev
            coef = coefold = start
        if coef is None:
            raise UnportedState("a GLM fit that stops before its first coefficients")
        estimate = coef.copy()
        aliased = [pivot[k] for k in range(rank, nvars)]
        estimate[aliased] = NA_REAL
        residuals = (y - mu) / fam.mu_eta(eta)
        weights = np.zeros(nobs)
        weights[good] = w**2
        df_residual = nobs - rank
        if not fam.estimate_dispersion:
            dispersion = 1.0
        elif df_residual > 0:
            dispersion = float(np.sum((weights * residuals**2)[weights > 0])) / df_residual
        else:
            dispersion = float("nan")
        p_values = np.full(nvars, NA_REAL)
        if rank > 0:
            r_factor = np.triu(factor[:rank, :rank])
            r_inverse = np.linalg.solve(r_factor, np.eye(rank))
            unscaled = np.einsum("ij,ij->i", r_inverse, r_inverse)  # diag(chol2inv(R))
            s_err = np.sqrt(dispersion * unscaled)
            t_value = coef[pivot[:rank]] / s_err
            if not fam.estimate_dispersion:
                p = 2.0 * special.ndtr(-np.abs(t_value))
            elif df_residual > 0:
                p = 2.0 * special.stdtr(df_residual, -np.abs(t_value))
            else:
                p = np.full(rank, np.nan)
            p_values[pivot[:rank]] = p
    return estimate, p_values


def _fastlm_coef_table(y: FloatArray, x: FloatArray) -> FloatArray:
    """``coef(summary(fastLm(...)))``: Estimate, StdErr, t.value and p.value of every
    coefficient (``fastLm`` of RcppArmadillo: a least-squares solve, and the standard errors
    from the pseudo-inverse of X'X)."""
    n, k = x.shape
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        coef = np.linalg.lstsq(x, y, rcond=None)[0]
        res = y - x @ coef
        s2 = np.float64(res @ res) / np.float64(n - k)
        std_err = np.sqrt(s2 * np.diag(np.linalg.pinv(x.T @ x)))
        t_value = coef / std_err
        p_value = 2.0 * special.stdtr(n - k, -np.abs(t_value)) if n - k > 0 else np.full(k, np.nan)
    return np.column_stack([coef, std_err, t_value, p_value])


# ---- the two association functions ------------------------------------------------------------
def _regression_nodes(
    seeds: Sequence[Name],
    x: RMatrix,
    Blist: Sequence[Name],
    covariates: Covariates | None,
    significance_threshold: float,
    adjust_p: str | None,
    glm_family: str,
    formula_string: str | None,
    fast: bool,
) -> AssociationData:
    families = (
        ("gaussian",)
        if fast
        else ("gaussian", "Gamma", "inverse.gaussian", "poisson", "quasi", "quasipoisson")
    )
    if glm_family not in families:
        raise RError("glm.family not supported!")
    rows = universe_rows(x, Blist)
    assert x.rownames is not None
    if len(rows) < 2:
        # x[rownames(x) %in% Blist, ] is a plain vector (or has no row): x[seed, ] fails
        raise RError(
            "incorrect number of dimensions" if len(rows) == 1 else "subscript out of bounds",
            "x[seed, ]",
        )
    genes = [x.rownames[i] for i in rows]
    if len(set(genes)) != len(genes):
        raise UnportedState("duplicated row names in the data matrix")
    values = x.values[rows]
    position = {name: i for i, name in enumerate(genes)}
    if covariates is not None and formula_string is not None and covariates.n_rows() != x.ncol:
        raise RError(
            f"arguments imply differing number of rows: {x.ncol}, {covariates.n_rows()}",
            "data.frame(x.node = rep(NA, ncol(x)), covariates, stringsAsFactors = FALSE)",
        )
    # R defect D3: the covariates enter the model only when formula.string names them
    terms = _parse_terms(formula_string, covariates)
    value_columns = ["Estimate", "StdErr", "t.value", "p.value"] if fast else ["Estimate", "P"]

    items: list[NamedVector] = []
    item_names: Names = []
    blocks: list[FloatArray] = []
    block_rows: Names = []
    for seed in seeds:
        if seed is None or seed not in position:
            raise RError("subscript out of bounds", "x[seed, ]")
        seed_exp = values[position[seed]]
        present = seed_exp[~np.isnan(seed_exp)]
        if len(present) > 0 and np.ptp(present) == 0:
            # R defect R11: a constant is fitted exactly by the intercept, the residuals are
            # rounding noise, and so are the p-values that decide the neighbourhood
            warnings.warn(
                f"seed {seed!r} has the same value in every sample: its neighbourhood is decided "
                "by rounding error, in R and here, and the two do not agree",
                NotReproducible,
                stacklevel=3,
            )
        beta = np.full(len(genes), NA_REAL)
        p = np.full(len(genes), NA_REAL)
        tables: list[FloatArray] = []
        coefficient_names: list[str] | None = None
        for g in range(len(genes)):
            y, design, names = _model_matrix(seed_exp, values[g], covariates, terms)
            if fast:
                table = _fastlm_coef_table(y, design)
                p_column = 3
            else:
                table = np.column_stack(_glm_coef_table(y, design, glm_family))
                p_column = 1
            if coefficient_names is None:
                coefficient_names = names
            elif len(names) != len(coefficient_names):
                raise RError("number of rows of matrices must match (see arg 2)", "cbind")
            node = names.index("x.node")
            beta[g] = table[node, 0]
            p[g] = table[node, p_column]
            tables.append(table)
        if adjust_p is not None:
            p = p_adjust(p, adjust_p)
        with np.errstate(invalid="ignore"):
            selected = np.flatnonzero(p < significance_threshold)  # which(): a missing p is skipped
        if len(selected) == 1:
            # R defect R5: one row and one column of the table: the name is dropped
            items.append(NamedVector(beta[selected], None))
            item_names.append(seed)
        elif len(selected) > 1:
            items.append(NamedVector(beta[selected], [genes[i] for i in selected]))
            item_names.append(seed)
        # R defect R6: no gene selected: glm.res.beta[[seed]] <- NULL adds nothing to the list
        assert coefficient_names is not None
        blocks.append(np.hstack(tables))
        block_rows += [f"{seed}.{name}" for name in coefficient_names]
    other = RMatrix(
        np.vstack(blocks),
        block_rows,
        [f"{gene}.{column}" for gene in genes for column in value_columns],
    )
    return AssociationData(RList(items, item_names if items else None), other)


def nodes_glm(
    seeds: Sequence[Name],
    x: RMatrix,
    Blist: Sequence[Name],
    covariates: Covariates | None = None,
    significance_threshold: float = 0.01,
    adjust_p: str | None = None,
    glm_family: str = "gaussian",
    formula_string: str | None = None,
    mc_cores: int = 1,
) -> AssociationData:
    """``nodes.glm``. ``mc_cores`` only spreads the fits over processes in R and is ignored."""
    return _regression_nodes(
        seeds,
        x,
        Blist,
        covariates,
        significance_threshold,
        adjust_p,
        glm_family,
        formula_string,
        fast=False,
    )


def nodes_fastLm(
    seeds: Sequence[Name],
    x: RMatrix,
    Blist: Sequence[Name],
    covariates: Covariates | None = None,
    significance_threshold: float = 0.01,
    adjust_p: str | None = None,
    glm_family: str = "gaussian",
    formula_string: str | None = None,
    mc_cores: int = 1,
) -> AssociationData:
    """``nodes.fastLm``: the same with least squares, gaussian family only."""
    return _regression_nodes(
        seeds,
        x,
        Blist,
        covariates,
        significance_threshold,
        adjust_p,
        glm_family,
        formula_string,
        fast=True,
    )
