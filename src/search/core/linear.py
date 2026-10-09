"""Linear-model association: the expression of a seed explained by the expression of a gene
plus covariates. The gaussian case is vectorised with the Frisch-Waugh-Lovell theorem: seed and
genes are residualised on the covariates once, then every seed-gene coefficient, its standard
error and p-value come from dot products of residuals. Other GLM families, and formulas with
interactions, go through statsmodels one pair at a time.
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
    beta[:, sxx == 0] = np.nan
    p[:, sxx == 0] = np.nan
    bad = ~np.isfinite(values).all(axis=1)
    beta[:, bad] = np.nan
    p[:, bad] = np.nan
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


def glm_loop(
    values: FloatArray,
    seed_idx: NDArray[np.intp],
    cfg: AssociationConfig,
    sample_table: pd.DataFrame | None,
) -> AssociationResult:
    """One statsmodels GLM per seed-gene pair: seed_exp ~ x_node + covariates (+ formula terms)."""
    try:
        import statsmodels.formula.api as smf
    except ImportError as exc:  # pragma: no cover
        raise ImportError("the glm association needs statsmodels: pip install search[glm]") from exc
    n_genes, n = values.shape
    if cfg.covariates and (sample_table is None or sample_table.shape[1] == 0):
        raise DataError(f"covariates {cfg.covariates} need a sample table")
    frame = pd.DataFrame(index=range(n))
    for name in cfg.covariates:
        if sample_table is None or name not in sample_table.columns:
            raise DataError(f"covariate {name!r} not in the sample table")
        frame[name] = sample_table[name].to_numpy()
    terms = ["x_node", *cfg.covariates]
    if cfg.formula:
        terms.append(cfg.formula)
    formula = "seed_exp ~ " + " + ".join(terms)
    family = _statsmodels_family(cfg.glm_family)
    use_t = cfg.glm_family != "poisson"
    beta = np.full((n_genes, len(seed_idx)), np.nan)
    p = np.full((n_genes, len(seed_idx)), np.nan)
    log.info("glm association: %d x %d models, family %s", n_genes, len(seed_idx), cfg.glm_family)
    for j, row in enumerate(seed_idx):
        frame["seed_exp"] = values[row]
        for g in range(n_genes):
            frame["x_node"] = values[g]
            try:
                # t statistics where the dispersion is estimated, as R's summary.glm does
                fit = smf.glm(formula, data=frame, family=family).fit(use_t=use_t)
            except Exception:  # singular fits, convergence failures
                continue
            if "x_node" in fit.params.index:
                beta[g, j] = float(fit.params["x_node"])
                p[g, j] = float(fit.pvalues["x_node"])
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
