"""Linear-model association: the expression of a seed explained by the expression of a gene
plus covariates. The gaussian case is vectorised with the Frisch-Waugh-Lovell theorem: seed and
genes are residualised on the covariates once, then every seed-gene coefficient, its standard
error and p-value come from dot products of residuals. Other GLM families, and formulas with
interactions, go through statsmodels one pair at a time.

Both paths follow the rule of the engine: a gene with a missing value, or with the same value in
every sample, gets no association (NaN) and so never an edge.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from scipy import stats as sps

from search.config import AssociationConfig
from search.core.association import AssociationResult
from search.errors import ConfigError, DataError

FloatArray = NDArray[np.float64]
log = logging.getLogger("search")


def design_matrix(
    sample_table: pd.DataFrame | None, covariates: list[str], n_samples: int
) -> tuple[FloatArray, list[str]]:
    """Intercept plus covariate columns; categoricals one-hot encoded without the first level."""
    columns: list[FloatArray] = [np.ones(n_samples)]
    names = ["intercept"]
    if covariates:
        if sample_table is None or sample_table.shape[1] == 0:
            raise DataError(f"covariates {covariates} need a sample table")
        missing = [c for c in covariates if c not in sample_table.columns]
        if missing:
            raise DataError(f"covariates not in the sample table: {missing}")
        if len(sample_table) != n_samples:
            raise DataError("sample table and matrix have different numbers of samples")
        for name in covariates:
            column = sample_table[name]
            if column.isna().any():
                raise DataError(f"covariate {name!r} has missing values")
            if pd.api.types.is_bool_dtype(column):
                columns.append(column.to_numpy(dtype=np.float64))
                names.append(name)
            elif pd.api.types.is_numeric_dtype(column):
                columns.append(column.to_numpy(dtype=np.float64))
                names.append(name)
            else:
                dummies = pd.get_dummies(
                    column.astype(str), prefix=name, drop_first=True, dtype=float
                )
                for dummy in dummies.columns:
                    columns.append(dummies[dummy].to_numpy(dtype=np.float64))
                    names.append(str(dummy))
    x = np.column_stack(columns)
    rank = np.linalg.matrix_rank(x)
    if rank < x.shape[1]:
        raise DataError("the covariate design is rank deficient (collinear covariates)")
    return x, names


def usable_rows(values: FloatArray) -> NDArray[np.bool_]:
    """Rows a regression can use: no missing value and at least two different values.

    A constant row is told by its largest value being its smallest: its variance, computed in
    floating point, is often a tiny number and not zero.
    """
    with np.errstate(invalid="ignore"):
        varying = np.ptp(values, axis=1) > 0
    return np.asarray(np.isfinite(values).all(axis=1) & varying, dtype=bool)


def ols_fwl(
    values: FloatArray, seed_idx: NDArray[np.intp], design: FloatArray
) -> AssociationResult:
    """seed ~ gene + covariates for every seed and gene, through residuals on the covariates."""
    n = values.shape[1]
    q, _ = np.linalg.qr(design)
    residual = values - (values @ q) @ q.T  # genes x samples
    seeds = residual[seed_idx]  # K x N
    sxy = seeds @ residual.T  # K x G
    sxx = (residual * residual).sum(axis=1)  # G
    syy = (seeds * seeds).sum(axis=1)  # K
    df = n - design.shape[1] - 1
    with np.errstate(invalid="ignore", divide="ignore"):
        beta = sxy / sxx[None, :]
        rss = np.maximum(syy[:, None] - beta * sxy, 0.0)
        se = np.sqrt(rss / df / sxx[None, :]) if df > 0 else np.full_like(beta, np.nan)
        t = beta / se
    p = 2.0 * sps.t.sf(np.abs(t), df=df) if df > 0 else np.full_like(beta, np.nan)
    usable = usable_rows(values)
    # a constant gene leaves residuals of rounding noise, not of zeros: it is told beforehand
    no_gene = ~usable | (sxx == 0)
    beta[:, no_gene] = np.nan
    p[:, no_gene] = np.nan
    no_seed = ~usable[seed_idx]
    beta[no_seed, :] = np.nan
    p[no_seed, :] = np.nan
    for j, row in enumerate(seed_idx):
        if np.isfinite(beta[j, row]):
            p[j, row] = 0.0  # the seed explains itself exactly
    beta_gk = np.ascontiguousarray(beta.T)
    p_gk = np.ascontiguousarray(p.T)
    return AssociationResult(stat=beta_gk, p=p_gk, method="ols", n_samples=n)


def covariate_coefficients(
    values: FloatArray,
    seed_idx: NDArray[np.intp],
    design: FloatArray,
    names: list[str],
    seeds: list[str],
) -> pd.DataFrame:
    """Coefficients of the covariates alone for every seed (seed ~ covariates)."""
    coef, *_ = np.linalg.lstsq(design, values[seed_idx].T, rcond=None)
    table = pd.DataFrame(coef.T, columns=names)
    table.insert(0, "seed", seeds)
    return table


def _statsmodels_family(name: str) -> object:
    import statsmodels.api as sm

    families = {
        "gaussian": sm.families.Gaussian,
        "poisson": sm.families.Poisson,
        "gamma": sm.families.Gamma,
        "inverse_gaussian": sm.families.InverseGaussian,
    }
    return families[name]()


def quoted(name: str) -> str:
    """A column name as a term of a model formula.

    Written bare, a name is parsed as Python code, and one such as ``t.rfs`` or ``tumour grade``
    fails; inside ``Q()`` the column is read whatever characters its name holds.
    """
    return f"Q({name!r})"


def glm_loop(
    values: FloatArray,
    seed_idx: NDArray[np.intp],
    cfg: AssociationConfig,
    sample_table: pd.DataFrame | None,
) -> AssociationResult:
    """One statsmodels GLM per seed-gene pair: seed_exp ~ x_node + covariates (+ formula terms).

    A model that cannot be built is a configuration error. A fit that fails gives no
    association to that pair and is counted: the count is logged, and a run in which every fit
    fails stops, so that a network cannot come out empty without a word.
    """
    try:
        import statsmodels.formula.api as smf
    except ImportError as exc:  # pragma: no cover
        raise ImportError("the glm association needs statsmodels: pip install search[glm]") from exc
    n_genes, n = values.shape
    reserved = sorted({"x_node", "seed_exp"} & set(cfg.covariates))
    if reserved:
        raise ConfigError(f"covariates cannot be named {reserved}: the model uses these names")
    # the checks of the ols path: columns present, no missing value, no collinear covariates
    design_matrix(sample_table, cfg.covariates, n)
    frame = pd.DataFrame(index=range(n))
    for name in cfg.covariates:
        assert sample_table is not None
        frame[name] = sample_table[name].to_numpy()
    terms = ["x_node", *(quoted(name) for name in cfg.covariates)]
    if cfg.formula:
        terms.append(cfg.formula)
    formula = "seed_exp ~ " + " + ".join(terms)
    family = _statsmodels_family(cfg.glm_family)
    use_t = cfg.glm_family != "poisson"
    try:
        # built once on placeholder columns, before any fit
        placeholders = frame.assign(seed_exp=1.0, x_node=np.arange(n, dtype=np.float64))
        smf.glm(formula, data=placeholders, family=family)
    except Exception as exc:  # patsy raises its own errors for a term it cannot read
        raise ConfigError(f"the model {formula!r} cannot be built: {exc}") from exc
    beta = np.full((n_genes, len(seed_idx)), np.nan)
    p = np.full((n_genes, len(seed_idx)), np.nan)
    usable = usable_rows(values)
    genes = np.flatnonzero(usable)
    log.info(
        "glm association: %d x %d models, family %s", len(genes), len(seed_idx), cfg.glm_family
    )
    n_fits = 0
    failed = 0
    first_error = ""
    for j, row in enumerate(seed_idx):
        if not usable[row]:
            continue
        frame["seed_exp"] = values[row]
        for g in genes:
            frame["x_node"] = values[g]
            n_fits += 1
            try:
                # t statistics where the dispersion is estimated, as R's summary.glm does
                fit = smf.glm(formula, data=frame, family=family).fit(use_t=use_t)
            except Exception as exc:  # a fit this family cannot make on these values
                failed += 1
                first_error = first_error or f"{type(exc).__name__}: {exc}"
                continue
            if "x_node" in fit.params.index:
                beta[g, j] = float(fit.params["x_node"])
                p[g, j] = float(fit.pvalues["x_node"])
    if n_fits and failed == n_fits:
        raise DataError(
            f"every one of the {n_fits} GLM fits failed (family {cfg.glm_family}); "
            f"the first error was {first_error}"
        )
    if failed:
        log.warning(
            "glm association: %d of %d fits failed and give no association; the first error was %s",
            failed,
            n_fits,
            first_error,
        )
    return AssociationResult(stat=beta, p=p, method="glm", n_samples=n)


def associate_linear(
    values: FloatArray,
    seed_idx: NDArray[np.intp],
    cfg: AssociationConfig,
    sample_table: pd.DataFrame | None,
    seed_names: list[str] | None = None,
) -> AssociationResult:
    if cfg.method == "ols":
        if cfg.formula:
            raise ConfigError("formula terms need method glm; ols handles additive covariates only")
        design, names = design_matrix(sample_table, cfg.covariates, values.shape[1])
        result = ols_fwl(values, seed_idx, design)
        if cfg.covariates:
            result.coef_table = covariate_coefficients(
                values, seed_idx, design, names, seed_names or [str(i) for i in seed_idx]
            )
        return result
    if cfg.method == "glm":
        return glm_loop(values, seed_idx, cfg, sample_table)
    raise ConfigError(f"unknown linear association method {cfg.method}")
